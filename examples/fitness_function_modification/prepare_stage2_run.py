import os,sys
import glob
from ateamopt.nwb_extractor import NWB_Extractor
from ateamopt.model_parameters import AllActive_Model_Parameters
from ateamopt.utils import utility
from ateamopt.optim_config_rules import correct_voltage_feat_std,\
            adjust_param_bounds,filter_feat_proto_basic,entries_to_remove
from ateamopt.jobscript.jobmodule import test_JobModule,\
            PBS_JobModule,Slurm_JobModule,SGE_JobModule,ChainSubJob
import shutil
import logging
from collections import OrderedDict

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

def filter_feat_proto_active(features_dict,protocols_dict,all_protocols_dict,
                             select_dict, add_fi_kink = True,
                             add_DB_check = True):
    """
        Filter the features and protocols for the final
        stage of optimization
    """
    features_dict = correct_voltage_feat_std(features_dict)
    spiking_proto_dict = OrderedDict()
    non_spiking_proto_dict =  OrderedDict()
    training_stimtype_reject = ['LongDCSupra','Ramp','Short_Square_Triple','Noise']
    
    for feat_key,feat_val in features_dict.items():
        if any(reject_stim in feat_key for reject_stim in training_stimtype_reject):
            continue
        stim_amp = protocols_dict[feat_key]['stimuli'][0]['amp']
        if feat_val['soma']['Spikecount'][0] > 0:
            spiking_proto_dict[feat_key] = stim_amp
        else:
            non_spiking_proto_dict[feat_key] = stim_amp
            if 'depol_block' in feat_val['soma'].keys():
                del feat_val['soma']['depol_block']
    
    
    # Ignoring spiking protocol which are followed by non-spiking stim protocol
    max_nospiking_amp = max(non_spiking_proto_dict.values()) 
    f_key_list = []           
    for spike_stim,spike_amp in spiking_proto_dict.items():
        if spike_amp < max_nospiking_amp:
           f_key_list.append(spike_stim)
    spiking_proto_dict = entries_to_remove(f_key_list,spiking_proto_dict)       
           
           
    spiking_proto_sorted = sorted(spiking_proto_dict,
                           key=spiking_proto_dict.__getitem__)
    non_spiking_proto_sorted = sorted(non_spiking_proto_dict,
                           key=non_spiking_proto_dict.__getitem__)

    num_spike = select_dict['spike_proto'] # In descending order
    num_nospike = select_dict['nospike_proto'] # in descending order

    # Select spiking proto
    try:
        spiking_proto_select = spiking_proto_sorted[len(spiking_proto_sorted)-num_spike:\
                                                    len(spiking_proto_sorted)]
    except:
        logger.debug('Number of spiking protocols requested exceeds data')
        spiking_proto_select = spiking_proto_sorted

    # Select nospiking proto
    try:
        nonspiking_proto_select = non_spiking_proto_sorted[len(non_spiking_proto_sorted)-\
                                       num_nospike:len(non_spiking_proto_sorted)]
    except:
        logger.debug('Number of nonspiking protocols requested exceeds data')
        nonspiking_proto_select = non_spiking_proto_sorted

    if add_fi_kink:
        spiking_proto_select.append(spiking_proto_sorted[0]) # the fist spiking stim
        nonspiking_proto_select.append(non_spiking_proto_sorted[-1]) # the last non-spiking stim
        spiking_proto_select = list(set(spiking_proto_select))
        nonspiking_proto_select = list(set(nonspiking_proto_select))

    features_dict_filtered = {key:val for key,val in features_dict.items() \
                              if key in spiking_proto_select+
                                      nonspiking_proto_select}
    protocols_dict_filtered = {key:val for key,val in protocols_dict.items() \
                              if key in spiking_proto_select+
                                          nonspiking_proto_select}
    untrained_features_dict = {key:val for key,val in features_dict.items() \
                              if key not in spiking_proto_select+
                                  nonspiking_proto_select}
        
                
    if add_DB_check:
        max_proto_key = spiking_proto_sorted[-1]
        max_amp = max([proto['amp'] for proto in protocols_dict[max_proto_key]['stimuli']])
        DB_proto_delay = max([proto['delay'] for proto in protocols_dict[max_proto_key]['stimuli']])
        DB_proto_duration = min([proto['duration'] for proto in protocols_dict[max_proto_key]['stimuli']])
        DB_proto_stimend = min([proto['stim_end'] for proto in protocols_dict[max_proto_key]['stimuli']])
        DB_proto_totdur = min([proto['totduration'] for proto in protocols_dict[max_proto_key]['stimuli']])

        DB_holding_proto = [proto for proto in protocols_dict[max_proto_key]['stimuli'] \
                        if proto['delay'] == 0]
        DB_proto_dict = [{
                    'type' :'SquarePulse',
                    'amp' : max_amp + 0.01,
                    'delay' : DB_proto_delay,
                    'duration' : DB_proto_duration,
                    'stim_end' : DB_proto_stimend,
                    'totduration': DB_proto_totdur
                    }]
        if bool(DB_holding_proto):
            DB_proto_dict.append(DB_holding_proto[0])
        DB_feature_dict = {
                'depol_block' : [1.0,
                                 0.05]
                          }
        features_dict_filtered['DB_check_DC'] = {'soma': DB_feature_dict}
        protocols_dict_filtered['DB_check_DC'] = {'stimuli':DB_proto_dict}
        all_protocols_dict['DB_check_DC'] = {'stimuli':DB_proto_dict}

    return features_dict_filtered,untrained_features_dict,\
         protocols_dict_filtered,all_protocols_dict

def main():

    parent_dir = os.path.abspath(os.path.join('.', os.pardir))
    script_repo_dirname = 'Script_Repo'
    script_repo_dir = os.path.join(parent_dir, script_repo_dirname)
    path_to_cell_metadata = glob.glob(parent_dir+'/cell_metadata*.json')[0]
    cell_metadata=utility.load_json(path_to_cell_metadata)

    acceptable_stimtypes = ['Long Square','Ramp', 'Square - 2s Suprathreshold',
                        'Short Square - Triple','Noise 1','Noise 2']
    cell_id = cell_metadata['Cell_id']
    perisomatic_model_id = cell_metadata['Perisomatic_id']
    dend_type = cell_metadata['Dendrite_type']
    species = cell_metadata['Species']
    me_type = cell_metadata.pop('ME_type',None)

    # Remove ME type from metadata for backward compatibility
    utility.save_json(path_to_cell_metadata,cell_metadata)

    if me_type:
        cell_type = 'exc' if 'Exc' in me_type else 'inh'
    else:
        if dend_type == 'spiny':
            cell_type = 'exc'
        elif dend_type == 'aspiny':
            cell_type = 'inh'
        else:
            raise Exception('cell-type ambiguous')

    wasabi_bucket = 's3://aibs.test.ani/'
    wasabi_bucket += '%s/%s'%(species.replace(' ',''),cell_id)

    # Get the conda environment and nwb processing type
    if sys.argv[-1] == 'non_standard_nwb':
        non_standard_nwb = True
        conda_env = sys.argv[-2]
    else:
        conda_env = sys.argv[-1]
        non_standard_nwb = False

    # Extract data and get the features for the stage
    nwb_handler = NWB_Extractor(cell_id)
    ephys_data_path,stimmap_filename = nwb_handler.save_cell_data\
                (acceptable_stimtypes,non_standard_nwb=non_standard_nwb)
    feature_file = 'feature_set_stage2_exc.json' \
                if cell_type == 'exc' else \
                    'feature_set_stage2_inh.json'

    feature_set_repo = os.path.abspath(os.path.join(script_repo_dir,feature_file))
    feature_set_repo = feature_set_repo \
            if os.path.exists(feature_set_repo) else None
    feature_default_template = utility.locate_template_file(os.path.join('parameters',\
                        'feature_set_stage2.json')) # default includes check_AIS_initiation
    feature_path = feature_set_repo or feature_default_template

    filter_rule_func = filter_feat_proto_active
    select_dict = {'spike_proto': 2,
                   'nospike_proto' :0}
    add_fi_kink = True

    # Don't get features from these stim_types
    feature_reject_stim_type = ['Ramp','Short_Square_Triple']
    spiketimes_exp_path = 'Validation_Responses/spiketimes_exp_noise.pkl'
    train_features,test_features,all_features,train_protocols,all_protocols = \
        nwb_handler.get_ephys_features(feature_path,ephys_data_path,
                   stimmap_filename,filter_rule_func,select_dict,
                   add_fi_kink,feature_reject_stim_type= feature_reject_stim_type,
                   spiketimes_exp_path=spiketimes_exp_path)
    
    
    features_write_path,untrained_features_write_path,all_features_write_path,\
        protocols_write_path,all_protocols_write_path = \
        nwb_handler.write_ephys_features(train_features,test_features,\
                             all_features,train_protocols,all_protocols)
    
    
    # Get the basic ephys features for first round of training
    train_features_basic,train_protocols_basic = \
       filter_feat_proto_basic(train_features,train_protocols)
        
    features_basic_path = os.path.join('config/{}'.format(cell_id),
                                       'features_basic.json')
    protocols_basic_path = os.path.join('config/{}'.format(cell_id),
                                       'protocols_basic.json')
    
    utility.save_json(features_basic_path,train_features_basic)    
    utility.save_json(protocols_basic_path,train_protocols_basic)  
    
    
    # Create the parameter bounds for the optimization
    model_params_handler = AllActive_Model_Parameters(cell_id)
    morph_path = model_params_handler.swc_path

    if species == 'Mus musculus':
        param_bounds_file = 'param_bounds_stage2_mouse_exc.json' \
            if cell_type == 'exc' else \
                'param_bounds_stage2_mouse_inh.json'
    else:
        param_bounds_file = 'param_bounds_stage2_human_exc.json' \
            if cell_type == 'exc' else \
                'param_bounds_stage2_human_inh.json'

    param_bounds_repo = os.path.abspath(os.path.join(script_repo_dir,param_bounds_file))
    param_bounds_repo = param_bounds_repo \
            if os.path.exists(param_bounds_repo) else None
    param_bounds_default_template = utility.locate_template_file(os.path.join('parameters',\
                    'param_bounds_stage2.json')) #default parameters: Mouse spiny
    param_bounds_path = param_bounds_repo or param_bounds_default_template
    param_rule_func = adjust_param_bounds
    model_params,model_params_release= model_params_handler.get_opt_params\
                                (param_bounds_path,param_rule_func)
    param_write_path,release_param_write_path,release_params=\
                        model_params_handler.write_params_opt(model_params,model_params_release)
    model_mechs,model_mechs_release = model_params_handler.get_opt_mechanism(model_params,\
                        model_params_release,param_bounds_path)
    mech_write_path,mech_release_write_path = model_params_handler.write_mechanisms_opt(model_mechs,\
                        model_mechs_release)

    props = {}

    # Perisomatic model
    if perisomatic_model_id != '':
        perisomatic_dir = 'peri_model'
        peri_model_path = os.path.join(perisomatic_dir,'*fit_peri*.json')
        peri_param_path = glob.glob(peri_model_path)[0]
        peri_params_write_path, peri_mech_write_path = \
                model_params_handler.aibs_peri_to_bpopt(peri_param_path)
        props['peri_parameters'] = peri_params_write_path
        props['peri_mechanism'] = peri_mech_write_path

    # Config file with all the necessary paths to feed into the optimization
    model_params_handler.write_opt_config_file(morph_path,param_write_path,
                                  mech_write_path,mech_release_write_path,
                                  features_write_path,untrained_features_write_path,
                                  all_features_write_path,
                                  protocols_write_path,all_protocols_write_path,
                                  release_params,release_param_write_path,
                                  **props)

    # Basic config file
    model_params_handler.write_opt_config_file(morph_path,param_write_path,
                                  mech_write_path,mech_release_write_path,
                                  features_basic_path,untrained_features_write_path,
                                  all_features_write_path,
                                  protocols_basic_path,all_protocols_write_path,
                                  release_params,release_param_write_path,
                                  opt_config_filename = 'config_file_basic.json',
                                  **props)

    # Copy the optimization files in the current directory

    optimizer_script=utility.locate_script_file('Optim_Main.py')
    stage_cwd = os.getcwd()

    for script_path in [optimizer_script]:
        shutil.copy(script_path,stage_cwd)


    # Create batch jobscript
    machine = cell_metadata['Machine']
    if 'hpc-login' in machine:
        jobtemplate_path = 'job_templates/Stage2_pbs.sh'
        batch_job = PBS_JobModule(jobtemplate_path,machine,conda_env=conda_env)
        batch_job.script_generator()
        analysis_jobtemplate_path = 'job_templates/Stage2_analyze_template_pbs.sh'

    elif any(substring in machine for substring in ['cori', 'bbp']):
        jobtemplate_path = 'job_templates/Stage2_slurm_parallel_seed.sh'
        batch_job = Slurm_JobModule(jobtemplate_path,machine,conda_env=conda_env)
        batch_job.script_generator()
        analysis_jobtemplate_path = 'job_templates/Stage2_analyze_template.sh'

    elif 'aws' in machine:
        jobtemplate_path = 'job_templates/Stage2_sge.sh'
        batch_job = SGE_JobModule(jobtemplate_path,machine,conda_env=conda_env)
        batch_job.script_generator()
        analysis_jobtemplate_path = 'job_templates/Stage2_analyze_template_sge.sh'

    else:
        cp_dir = 'checkpoints'
        testJob = test_JobModule(machine,'batch_job.sh','%s/seed1.pkl'%cp_dir,
                                 2,2)
        testJob.script_generator()
        analysis_cmd = 'python analysis_stage2.py -vv --cp_dir  checkpoints \n'
        testJob.adjust_template('bash chain_job.sh',analysis_cmd)

    config_replace_cmd = 'sed -i -e "s/config_path config_file.json/' \
        'config_path config_file_basic.json/g" {}'.format(batch_job.script_name)
    os.system(config_replace_cmd)
    
    # fitness function modification
    ffmod_script = os.path.join(script_repo_dir,'ffmod_script.sh')
    with open(ffmod_script,'r') as templ:
            ffmod_cmd = templ.read()
    ffmod_cmd = ffmod_cmd.replace('qsub',batch_job.submit_verb)
    batch_job.adjust_template('analyze_results.sh','',partial_match = True)
    batch_job.adjust_template('mv ${CHECKPOINTS_DIR}/*',ffmod_cmd,
                              partial_match = True)
    
    # Create Analysis job for final stage
    if any(substring in machine for substring in ['cori','bbp',\
                                                  'hpc-login','aws']):
        analysis_job = ChainSubJob(analysis_jobtemplate_path,machine,\
                        script_name = 'analyze_results.sh',conda_env=conda_env)
        analysis_job.script_generator()

    

    # Transfer data to Wasabi only for AWS
    if 'aws' in machine:
        s3_transfer_cmd = '# aws s3 mv %s %s --recursive --profile wasabi\n'%(parent_dir,wasabi_bucket)
        analysis_job.adjust_template('python analysis_stage2.py',s3_transfer_cmd,
                                     add=True,partial_match=True)

if __name__ == '__main__':
    main()
