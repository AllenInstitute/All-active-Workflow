import os,sys
import glob
from ateamopt.nwb_extractor import NWB_Extractor
from ateamopt.model_parameters import AllActive_Model_Parameters
from ateamopt.utils import utility
from ateamopt.optim_config_rules import filter_feat_proto_active,\
                                adjust_param_bounds
from ateamopt.jobscript.jobmodule import test_JobModule,\
            PBS_JobModule,Slurm_JobModule,ChainSubJob

import shutil
import logging
logger = logging.getLogger(__name__)

def main():

    parent_dir = os.path.abspath(os.path.join('.', os.pardir))
    path_to_cell_metadata = glob.glob(parent_dir+'/cell_metadata*.json')[0]
    cell_metadata=utility.load_json(path_to_cell_metadata)

    acceptable_stimtypes = ['Long Square','Ramp', 'Square - 2s Suprathreshold',
                        'Short Square - Triple','Noise 1','Noise 2']
    cell_id = cell_metadata['Cell_id']
    species = cell_metadata['Species']
    wasabi_bucket = 's3://aibs.snmo.01/'
    wasabi_bucket += '%s/%s'%(species.replace(' ',''),cell_id)

    # Extract data and get the features for the stage
    nwb_handler = NWB_Extractor(cell_id)
    ephys_data_path,stimmap_filename = nwb_handler.save_cell_data\
                                    (acceptable_stimtypes)
    feature_path = 'parameters/feature_set_stage2.json'
    filter_rule_func = filter_feat_proto_active
    select_dict = {'spike_proto': 2,
                   'nospike_proto' :0}
    add_fi_kink = True

    features_write_path,untrained_features_write_path,\
        protocols_write_path,all_protocols_write_path = \
        nwb_handler.get_ephys_features(feature_path,ephys_data_path,
                       stimmap_filename,filter_rule_func,select_dict,add_fi_kink)


    # Create the parameter bounds for the optimization
    model_params_handler = AllActive_Model_Parameters(cell_id)
    morph_path = model_params_handler.swc_path
    param_bounds_path = 'parameters/param_bounds_stage2.json'
    param_rule_func = adjust_param_bounds
    model_params,model_params_release= model_params_handler.get_opt_params\
                                (param_bounds_path,param_rule_func)
    param_write_path,release_param_write_path,release_params=\
                        model_params_handler.write_params_opt(model_params,model_params_release)
    mech_write_path,mech_release_write_path = model_params_handler.write_mechanisms_opt(model_params,\
                                    model_params_release,param_bounds_path)

    # Config file with all the necessary paths to feed into the optimization
    model_params_handler.write_opt_config_file(morph_path,param_write_path,
                                  mech_write_path,mech_release_write_path,
                                  features_write_path,untrained_features_write_path,
                                  protocols_write_path,all_protocols_write_path,
                                  release_params,release_param_write_path)

    # Copy the optimization files in the current directory

    optimizer_script=utility.locate_script_file('Optim_Main.py')
    stage_cwd = os.getcwd()

    for script_path in [optimizer_script]:
        shutil.copy(script_path,stage_cwd)

    # Get the conda environment
    conda_env = sys.argv[-1]

    # Create batch jobscript
    machine = cell_metadata['Machine']
    if 'hpc-login' in machine:
        jobtemplate_path = 'job_templates/Stage2_pbs.sh'
        batch_job = PBS_JobModule(jobtemplate_path,machine,conda_env=conda_env)
        batch_job.script_generator()
        analysis_jobtemplate_path = 'job_templates/Stage2_analyze_template_pbs.sh'

    elif any(substring in machine for substring in ['cori', 'bbp']):
        jobtemplate_path = 'job_templates/Stage2_slurm.sh'
        batch_job = Slurm_JobModule(jobtemplate_path,machine,conda_env=conda_env)
        batch_job.script_generator()
        analysis_jobtemplate_path = 'job_templates/Stage2_analyze_template.sh'

    else:
        cp_dir = 'checkpoints'
        testJob = test_JobModule(machine,'batch_job.sh','%s/seed1.pkl'%cp_dir,
                                 2,2)
        testJob.script_generator()
        analysis_cmd = 'python analysis_stage2.py -vv --cp_dir  checkpoints \n'
        analysis_cmd += 'aws s3 cp %s %s --recursive --profile wasabi\n'%(parent_dir,wasabi_bucket)
        testJob.adjust_template('sh chain_job.sh',analysis_cmd)


    # Create Analysis job for final stage
    if any(substring in machine for substring in ['cori', 'bbp','hpc-login']):
        analysis_job = ChainSubJob(analysis_jobtemplate_path,machine,\
                        script_name = 'analyze_results.sh',conda_env=conda_env)
        analysis_job.script_generator()
        s3_transfer_cmd = 'aws s3 cp %s %s --recursive --profile wasabi\n'%(parent_dir,wasabi_bucket)
        analysis_job.adjust_template('rm -rf $IPYTHONDIR',s3_transfer_cmd,add=True)

if __name__ == '__main__':
    main()
