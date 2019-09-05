#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 12:39:10 2019

@author: anirbannandi
"""

from ateamopt.utils import utility


job_config_dict = {
        
    # Cell level config
        
    'cty_config' : { 
        'cell_id': '548494556',
        'e_type' : 'Exc',
        'me_type' : '',
        },
    
     # Job level config        
     
    'job_config' : {

        'highlevel_jobconfig' : 
                {  
                    'job_dir': 'Human_trial',
                    'conda_env' : 'ateam_opt',
                    'nwb_path' : '',
                    'swc_path' : '',
                    'nwb_type' : 'standard',
                    'axon_type' : 'stub_axon', # biological
                    'ephys_dir' : 'ephys_data',
                    'non_standard_nwb': False,
                    'acceptable_stimtypes' : ['Long Square', 'Noise 1', 'Noise 2'],
                    'feature_names_path' : 'feature_set_all.json',
                    'dryrun' : True,
#                    'script_repo_dir' : 'Script_Repo'
#                    'modfiles_dir' : 'modfiles'
                    'compiled_modfiles_dir' : 'x86_64'
                },
            
        'stage_jobconfig' : [        
                 {
                    'stage_name' : 'Stage0',
                    'stage_stimtypes': ['Long Square'],
                    'stage_features':'feature_set_stage0.json',
                    'stage_parameters' : 'param_bounds_stage0.json',
                    'cp_dir' : 'checkpoint',
                    'filter_rule' : 'filter_feat_proto_passive',
                    'nengines' : 256,
                    'nnodes' : 16,
                    'qos' : 'celltypes',
                    'nprocs': 16,
                    'error_stream' : 'job.err',
                    'output_stream' : 'job.out',
                    'jobmem' : '100g',
                    'jobtime' : '5:00:00',
                    'ipyp_db' : 'nodb',
                    'seed' : [1]
                 },
                 {
                    'stage_name' : 'Stage1',
                    'stage_stimtypes': ['Long Square'],
                    'stage_features':'feature_set_stage1.json',
                    'stage_parameters' : 'param_bounds_stage1.json',
                    'filter_rule' : 'filter_feat_proto_passive',
                    'nengines' : 256,
                    'nnodes' : 16,
                    'qos' : 'celltypes',
                    'nprocs': 16,
                    'error_stream' : 'job.err',
                    'output_stream' : 'job.out',
                    'jobmem' : '100g',
                    'jobtime' : '5:00:00',
                    'ipyp_db' : 'nodb',
                    'seed' : [1]
                 },
                 {
                    'stage_name' : 'Stage2',
                    'stage_stimtypes': ['Long Square'],
                    'stage_features':'feature_set_stage2.json',
                    'stage_parameters' : 'param_bounds_stage2.json',
                    'filter_rule' : 'filter_feat_proto_active',
                    'nengines' : 256,
                    'nnodes' : 16,
                    'qos' : 'celltypes',
                    'nprocs': 16,
                    'error_stream' : '/dev/null',
                    'output_stream' : '/dev/null',
                    'jobmem' : '100g',
                    'jobtime' : '5:00:00',
                    'nengines_analysis' : 10,
                    'nnodes_analysis' : 1,
                    'nprocs_analysis': 10,
                    'jobtime_analysis' : '30:00',
                    'ipyp_db' : 'nodb',
                    'seed' : [1,2],
                    'run_hof_analysis' : True,
                    'run_peri_comparison' :True,
                    'depol_block_check' : True,
                    'add_fi_kink' : True,
                    'calc_model_perf' : True,
                    'ipyp_analysis' : True,
                    'error_stream_analysis' : 'analysis.err',
                    'output_stream_analysis' : 'analysis.out',
                    'model_postprocess' : True,
                    'calc_time_statistics' : True
                 }
                ]
            }
    }
utility.save_json('job_config.json',job_config_dict)
