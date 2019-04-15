import os,sys
import glob
from ateamopt.nwb_extractor import NWB_Extractor
from ateamopt.model_parameters import AllActive_Model_Parameters
from ateamopt.utils import utility
from ateamopt.optim_config_rules import filter_feat_proto_passive
from ateamopt.jobscript.jobmodule import test_JobModule,\
            PBS_JobModule,Slurm_JobModule,ChainSubJob
import shutil
import logging

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

def main():

    parent_dir = os.path.abspath(os.path.join('.', os.pardir))
    script_repo_dirname = 'Script_Repo'
    script_repo_dir = os.path.join(parent_dir, script_repo_dirname)
    path_to_cell_metadata = glob.glob(parent_dir+'/cell_metadata*.json')[0]
    cell_metadata=utility.load_json(path_to_cell_metadata)
    acceptable_stimtypes = ['Long Square']
    cell_id = cell_metadata['Cell_id']

    # Extract data and get the features for the stage
    nwb_handler = NWB_Extractor(cell_id)
    ephys_data_path,stimmap_filename = nwb_handler.save_cell_data(acceptable_stimtypes)
    feature_path = utility.locate_template_file(os.path.join('parameters',\
                        'feature_set_stage0.json'))
    filter_rule_func = filter_feat_proto_passive
    
    train_features,test_features,all_features,train_protocols,all_protocols = \
        nwb_handler.get_ephys_features(feature_path,ephys_data_path,
                                       stimmap_filename,filter_rule_func)
    features_write_path,untrained_features_write_path,all_features_write_path,\
        protocols_write_path,all_protocols_write_path = \
        nwb_handler.write_ephys_features(train_features,test_features,\
                             all_features,train_protocols,all_protocols)
    
    # Create the parameter bounds for the optimization
    model_params_handler = AllActive_Model_Parameters(cell_id)
    morph_path = model_params_handler.swc_path

    param_bounds_file = 'param_bounds_stage0.json'
    param_bounds_repo = os.path.abspath(os.path.join(script_repo_dir,param_bounds_file))
    param_bounds_repo = param_bounds_repo \
            if os.path.exists(param_bounds_repo) else None
    param_bounds_default_template = utility.locate_template_file(os.path.join('parameters',\
                                            param_bounds_file))
    param_bounds_path = param_bounds_repo or param_bounds_default_template
    model_params,model_params_release= model_params_handler.get_opt_params(param_bounds_path)
    param_write_path,release_param_write_path,release_params=\
                        model_params_handler.write_params_opt(model_params,model_params_release)
    
    model_mechs,model_mechs_release = model_params_handler.get_opt_mechanism(model_params,\
                        model_params_release,param_bounds_path)
    mech_write_path,mech_release_write_path = model_params_handler.write_mechanisms_opt(model_mechs,\
                        model_mechs_release)

    # Config file with all the necessary paths to feed into the optimization
    model_params_handler.write_opt_config_file(morph_path,param_write_path,
                                  mech_write_path,mech_release_write_path,
                                  features_write_path,untrained_features_write_path,
                                  all_features_write_path,
                                  protocols_write_path,all_protocols_write_path,
                                  release_params,release_param_write_path)

    # Copy the optimizer scripts in the current directory

    optimizer_script=utility.locate_script_file('Optim_Main.py')
    stage_cwd = os.getcwd()

    for script_path in [optimizer_script]:
        shutil.copy(script_path,stage_cwd)

    # Get the conda environment
    conda_env = sys.argv[-1]

    # Create batch jobscript
    machine = cell_metadata['Machine']
    if 'hpc-login' in machine:
        jobtemplate_path = 'job_templates/Stage0_pbs.sh'
        batch_job = PBS_JobModule(jobtemplate_path,machine,conda_env=conda_env)
        batch_job.script_generator()
        chain_jobtemplate_path = 'job_templates/Stage1_chainjob_template_pbs.sh'

    elif any(substring in machine for substring in ['cori', 'bbp']):
        jobtemplate_path = 'job_templates/Stage0_slurm.sh'
        batch_job = Slurm_JobModule(jobtemplate_path,machine,conda_env=conda_env)
        batch_job.script_generator()
        chain_jobtemplate_path = 'job_templates/Stage1_chainjob_template.sh'
    else:
        cp_dir = 'checkpoints'
        testJob = test_JobModule(machine,'batch_job.sh','%s/seed1.pkl'%cp_dir,
                                 2,2)
        testJob.script_generator()
        chain_jobtemplate_path = 'job_templates/Stage1_chainjob_template_pbs.sh'

    # Create Chain job for next stage

    chain_job = ChainSubJob(chain_jobtemplate_path,machine,conda_env=conda_env)
    chain_job.script_generator()


if __name__ == '__main__':
    main()












