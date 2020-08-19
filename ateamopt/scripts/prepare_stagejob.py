import os
import sys
import glob
from ateamopt.nwb_extractor import NwbExtractor
from ateamopt.model_parameters import AllActive_Model_Parameters
from ateamopt.utils import utility
import ateamopt.optim_config_rules as filter_rules
from ateamopt.jobscript.jobmodule import test_JobModule,\
    PBS_JobModule, Slurm_JobModule, SGE_JobModule, ChainSubJob
import shutil
import logging
import argschema as ags
from ateamopt.optim_schema import Stage_Launch_Config
from collections import defaultdict


logger = logging.getLogger()


def main(args):
    # Job config
    job_config_path = sys.argv[-1]
    stage_jobconfig = args['stage_jobconfig']
    highlevel_job_props = args['highlevel_jobconfig']

    logging.basicConfig(level=highlevel_job_props['log_level'])

    job_dir = highlevel_job_props['job_dir']
    path_to_cell_metadata = glob.glob(
        os.path.join(job_dir, 'cell_metadata*.json'))[0]
    stage_tracker_path = os.path.join(job_dir, 'stage_tracker_config.json')

    cell_metadata = utility.load_json(path_to_cell_metadata)
    cell_id = cell_metadata['cell_id']
    peri_model_id = cell_metadata.get('peri_model_id')
    released_aa_model_path = cell_metadata.get('model_path_all_active')
    released_aa_model_id = cell_metadata.get('released_aa_model_id')

    nwb_path = highlevel_job_props['nwb_path']
    swc_path = highlevel_job_props['swc_path']
    all_features_path = highlevel_job_props['all_features_path']
    all_protocols_path = highlevel_job_props['all_protocols_path']

    stage_stimtypes = stage_jobconfig['stage_stimtypes']
    stage_feature_names_path = stage_jobconfig['stage_features']
    param_bounds_path = stage_jobconfig['stage_parameters']
    ap_init_flag = stage_jobconfig['AP_initiation_zone']
    ap_init_feature = 'check_AISInitiation'
    script_repo_dir = stage_jobconfig.get('script_repo_dir')
    depol_block_check = stage_jobconfig.get('depol_block_check')
    add_fi_kink = stage_jobconfig.get('add_fi_kink')
    analysis_parallel = (stage_jobconfig['analysis_config'].get('ipyparallel') and
                         stage_jobconfig['run_hof_analysis'])  # analysis batch job only for hof analysis
    param_bound_tolerance = stage_jobconfig.get('adjust_param_bounds_prev')
    prev_stage_path = stage_jobconfig.get('prev_stage_path')

    filter_rule_func = getattr(filter_rules, stage_jobconfig['filter_rule'])

    all_features = utility.load_json(all_features_path)
    all_protocols = utility.load_json(all_protocols_path)

    stage_feature_names = utility.load_json(
        stage_feature_names_path)['features']
    # AP init flag is prioritized over feature set file
    if ap_init_flag == 'soma':
        if ap_init_feature in stage_feature_names:
            stage_feature_names.remove(ap_init_feature)
    elif ap_init_flag == 'axon':
        if ap_init_feature not in stage_feature_names:
            stage_feature_names.append(ap_init_feature)

    select_stim_names = []
    for stim_name in all_features.keys():
        stim_type = stim_name.rsplit('_', 1)[0]
        stim_type_aibs = utility.aibs_stimname_map_inv[stim_type]
        if stim_type_aibs in stage_stimtypes:
            select_stim_names.append(stim_name)

    features_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    for stim_name, stim_dict in all_features.items():
        if stim_name in select_stim_names:
            for loc, loc_features in stim_dict.items():
                for feat, val in loc_features.items():
                    if feat in stage_feature_names:
                        features_dict[stim_name][loc][feat] = [val[0], val[1]]

    protocols_dict = {proto_key: proto_val for proto_key, proto_val in all_protocols.
                      items() if proto_key in select_stim_names}
    nwb_handler = NwbExtractor(cell_id, nwb_path=nwb_path)

    kwargs = {'depol_block_check': depol_block_check,
              'add_fi_kink': add_fi_kink}
    if depol_block_check:
        train_features, test_features, train_protocols, DB_proto_dict = filter_rule_func(
            features_dict, protocols_dict, **kwargs)
        # also append DB check info to all_protocols json and save
        all_protocols['DB_check_DC'] = {'stimuli': DB_proto_dict}
        utility.save_json(all_protocols_path, all_protocols)
    else:
        train_features, test_features, train_protocols = filter_rule_func(
            features_dict, protocols_dict, **kwargs)

    train_features_path, test_features_path, train_protocols_path = \
        nwb_handler.write_ephys_features(train_features, test_features,
                                         train_protocols)

    # Create the parameter bounds for the optimization
    if prev_stage_path:
        prev_stage_model_path = os.path.join(prev_stage_path, 'fitted_params', 'optim_param_%s_compact.json'
                                             % cell_id)
    else:
        prev_stage_model_path = None
    model_params_handler = AllActive_Model_Parameters(cell_id, swc_path=swc_path,
                                                      prev_stage_model_path=prev_stage_model_path,
                                                      released_aa_model_path=released_aa_model_path)

    model_params, model_params_release = model_params_handler.get_opt_params(param_bounds_path,
                                                                             prev_stage_tolerance=param_bound_tolerance)
    param_write_path, released_aa_param_write_path, released_aa_params =\
        model_params_handler.write_params_opt(model_params, model_params_release)

    model_mechs, model_mechs_release = model_params_handler.get_opt_mechanism(model_params,
                                                                              model_params_release, param_bounds_path)
    mech_write_path, mech_release_write_path = model_params_handler.write_mechanisms_opt(model_mechs,
                                                                                         model_mechs_release)

    props = {}
    if peri_model_id:
        peri_model_path = cell_metadata['model_path_perisomatic']
        peri_params_write_path, peri_mech_write_path = \
            model_params_handler.aibs_peri_to_bpopt(peri_model_path)
        props['released_peri_model'] = peri_params_write_path
        props['released_peri_mechanism'] = peri_mech_write_path

    # Config file with all the necessary paths to feed into the optimization
    # TODO: clarify how this fits into schema
    model_params_handler.write_opt_config_file(param_write_path,
                                               mech_write_path, mech_release_write_path,
                                               train_features_path, test_features_path,
                                               train_protocols_path,
                                               released_aa_params, released_aa_param_write_path,
                                               opt_config_filename=job_config_path,
                                               **props)

    # Copy the optimizer scripts in the current directory
    optimizer_script = stage_jobconfig['optim_config']['main_script']
    analysis_script = stage_jobconfig['analysis_config']['main_script']
    if script_repo_dir:
        optimizer_script_repo = os.path.abspath(os.path.join(script_repo_dir,
                                                             optimizer_script))
        optimizer_script_repo = optimizer_script_repo if os.path.exists(optimizer_script_repo)\
            else None
    else:
        optimizer_script_repo = None
    optimizer_script_default = utility.locate_script_file(optimizer_script)
    optimizer_script_path = optimizer_script_repo or optimizer_script_default
    stage_cwd = os.getcwd()
    shutil.copy(optimizer_script_path, stage_cwd)

    next_stage_job_props = utility.load_json(stage_tracker_path)

    machine = highlevel_job_props['machine']
    machine_match_patterns = ['hpc-login', 'aws', 'cori', 'bbp5']

    next_stage_jobconfig = {}
    try:
        next_stage_jobconfig['stage_jobconfig'] = next_stage_job_props.pop(0)
        next_stage_jobconfig['highlevel_jobconfig'] = highlevel_job_props
        next_stage_jobconfig['stage_jobconfig']['prev_stage_path'] = os.getcwd()

        chainjobtemplate_path = 'job_templates/chainjob_template.sh'
    except:
        pass

    utility.save_json(stage_tracker_path, next_stage_job_props)

    # Create batch jobscript
    if not any(substr in machine for substr in machine_match_patterns):
        testJob = test_JobModule(
            'batch_job.sh', job_config_path=job_config_path)

        testJob.script_generator(next_stage_job_config=next_stage_jobconfig)

    elif any(pattern in machine for pattern in ['hpc-login', 'aws']):
        jobtemplate_path = 'job_templates/pbs_jobtemplate.sh'
        batch_job = PBS_JobModule(jobtemplate_path, job_config_path)
        if analysis_parallel:
            batch_job.script_generator(analysis_jobname='analyze_job.sh')
            # A separate batch job needs to be created in this case
            analysis_job = PBS_JobModule(jobtemplate_path, job_config_path,
                                         script_name='analyze_job.sh')
            analysis_job.script_generator(analysis=True,
                                          next_stage_job_config=next_stage_jobconfig)
        else:
            batch_job.script_generator(next_stage_job_config=next_stage_jobconfig)

    elif any(pattern in machine for pattern in ['cori', 'bbp5']):

        if 'cori' in machine:
            jobtemplate_path = 'job_templates/nersc_slurm_jobtemplate.sh'
        else:
            jobtemplate_path = 'job_templates/bbp_slurm_jobtemplate.sh'

        batch_job = Slurm_JobModule(jobtemplate_path, job_config_path)
        if analysis_parallel:
            batch_job.script_generator(analysis_jobname='analyze_job.sh')
            # A separate batch job needs to be created in this case
            analysis_job = Slurm_JobModule(jobtemplate_path, job_config_path,
                                           script_name='analyze_job.sh')
            analysis_job.script_generator(analysis=True,
                                          next_stage_job_config=next_stage_jobconfig)
        else:
            batch_job.script_generator(next_stage_job_config=next_stage_jobconfig)

    if next_stage_jobconfig:

        stage_jobdir = os.path.join(highlevel_job_props['job_dir'],
                                    next_stage_jobconfig['stage_jobconfig']['stage_name'])

        next_stage_jobconfig_path = os.path.join(
            stage_jobdir, 'stage_job_config.json')
        utility.create_filepath(next_stage_jobconfig_path)
        utility.save_json(next_stage_jobconfig_path, next_stage_jobconfig)
        prepare_jobscript_default = utility.locate_script_file(
            'prepare_stagejob.py')
        analyze_jobscript_default = utility.locate_script_file(analysis_script)
        shutil.copy(prepare_jobscript_default, stage_jobdir)
        shutil.copy(analyze_jobscript_default, stage_jobdir)

        chain_job = ChainSubJob(chainjobtemplate_path,
                                next_stage_jobconfig_path)
        chain_job.script_generator()


if __name__ == '__main__':
    mod = ags.ArgSchemaParser(schema_type=Stage_Launch_Config)
    main(mod.args)
