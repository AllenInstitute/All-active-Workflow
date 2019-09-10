import os
import ateamopt
from ateamopt.nwb_extractor import NwbExtractor
import ateamopt.cell_data as cell_data
from ateamopt.utils import utility
import logging
import glob
import shutil
from ateamopt.jobscript.jobmodule import ChainSubJob
import argschema as ags
from ateamopt.optim_schema import Launch_Config
from ateamopt.optim_config_rules import correct_feat_statistics
from ateamopt.morph_handler import MorphHandler
import subprocess
import bluepyopt
import pandas as pd

logger = logging.getLogger()


def convert_paths(prop_dict):
    job_dict_abs = {}
    for j_key, j_val in prop_dict.items():
        try:
            if os.path.exists(j_val):
                job_dict_abs[j_key] = os.path.abspath(j_val)
            else:
                job_dict_abs[j_key] = j_val
        except:
            job_dict_abs[j_key] = j_val
    return job_dict_abs


def create_optim_job(args):
    level = logging.getLevelName(args['log_level'])
    logger.setLevel(level)

    cty_props = args['cty_config']
    cell_id = cty_props['cell_id']
    highlevel_job_props = args['job_config']['highlevel_jobconfig']
    stage_job_props = args['job_config']['stage_jobconfig']

    # Change any paths to absolute path
    for ii, stage_job_prop in enumerate(stage_job_props):
        stage_job_props[ii] = convert_paths(stage_job_prop)
    highlevel_job_props = convert_paths(highlevel_job_props)

    try:
        job_dir = os.path.join(os.getcwd(), highlevel_job_props['job_dir'])
    except:
        job_dir = os.path.join(os.getcwd(), str(cell_id))
    highlevel_job_props['job_dir'] = job_dir

    utility.create_dirpath(job_dir)
    os.chdir(job_dir)  # Change Working directory

    cty_config_path = 'user_config/cell_config.json'
    job_config_path = 'user_config/job_config.json'
    highlevel_jobconfig_path = 'high_level_job_config.json'
    stage_tracker_path = 'stage_tracker_config.json'

    utility.create_filepath(cty_config_path)
    utility.create_filepath(job_config_path)

    # Save a copy of the config files
    utility.save_json(cty_config_path, cty_props)
    utility.save_json(job_config_path, args['job_config'])

    try:
        ateamopt_dir = os.path.join(
            os.path.dirname(ateamopt.__file__), os.pardir)
        ateamopt_commitID = subprocess.check_output(
            ["git", "describe", "--tags"], cwd=ateamopt_dir).strip()
        ateamopt_commitID = ateamopt_commitID.decode() if isinstance(
            ateamopt_commitID, bytes) else ateamopt_commitID
        cty_props['ateamopt_tag'] = ateamopt_commitID
    except Exception as e:
        logger.debug(e)
    try:
        bluepyopt_dir = os.path.join(
            os.path.dirname(bluepyopt.__file__), os.pardir)
        bpopt_commitID = subprocess.check_output(
            ["git", "describe", "--tags"], cwd=bluepyopt_dir).strip()
        bpopt_commitID = bpopt_commitID.decode() if isinstance(
            bpopt_commitID, bytes) else bpopt_commitID
        cty_props['bluepyopt_tag'] = bpopt_commitID
    except:
        pass
    
    # pickling consistency depends on pandas version
    pd_version = pd.__version__
    cty_props['pandas_version'] = pd_version

    cell_metadata_path = glob.glob('cell_metadata*.json')

    if len(cell_metadata_path) == 0:
        cty_props.update(highlevel_job_props)
        cell_metadata, cell_metadata_path = cell_data.save_cell_metadata(
            **cty_props)
        morph_stats_filename = 'morph_stats_%s.json' % cell_id
        morph_handler = MorphHandler(
            cell_metadata['swc_path'], cell_id=cell_id)
        morph_handler.save_morph_data(morph_stats_filename)

    elif len(cell_metadata_path) == 1:
        cell_metadata_path = cell_metadata_path[0]
        cell_metadata = utility.load_json(cell_metadata_path)

    else:
        raise Exception('More than one metadata files found')

    # Extract ephys data
    ephys_dir = highlevel_job_props['ephys_dir']
    non_standard_nwb = highlevel_job_props['non_standard_nwb']
    acceptable_stimtypes = highlevel_job_props['acceptable_stimtypes']

    highlevel_job_props['nwb_path'] = cell_metadata['nwb_path']
    highlevel_job_props['swc_path'] = cell_metadata['swc_path']
    nwb_handler = NwbExtractor(
        cell_id, nwb_path=highlevel_job_props['nwb_path'])
    ephys_data_path, stimmap_filename = nwb_handler.save_cell_data(acceptable_stimtypes, non_standard_nwb=non_standard_nwb,
                                                                   ephys_dir=ephys_dir)
    feature_names_path = highlevel_job_props['feature_names_path']
    protocol_dict,feature_dict = nwb_handler.get_efeatures_all(feature_names_path,
                                          ephys_data_path,stimmap_filename)
    
    feature_dict = correct_feat_statistics(feature_dict,protocol_dict)
    all_protocols_filename = os.path.join(ephys_data_path,'all_protocols.json')
    all_features_filename = os.path.join(ephys_data_path,'all_features.json')
    utility.save_json(all_protocols_filename,protocol_dict)
    utility.save_json(all_features_filename,feature_dict)
    highlevel_job_props['stimmap_file'] = os.path.abspath(stimmap_filename)
    highlevel_job_props['machine'] = cell_metadata['machine']
    highlevel_job_props['log_level'] = args['log_level']
    highlevel_job_props['all_features_path'] = all_features_filename
    highlevel_job_props['all_protocols_path'] = all_protocols_filename

    highlevel_job_props = convert_paths(highlevel_job_props)
    utility.save_json(highlevel_jobconfig_path, highlevel_job_props)

    stage_level_jobconfig = {}
    stage_level_jobconfig['stage_jobconfig'] = stage_job_props.pop(0)
    stage_level_jobconfig['highlevel_jobconfig'] = highlevel_job_props

    utility.save_json(stage_tracker_path, stage_job_props)

    stage_jobdir = os.path.join(highlevel_job_props['job_dir'],
                                stage_level_jobconfig['stage_jobconfig']['stage_name'])
    stage_level_jobconfig_path = os.path.join(
        stage_jobdir, 'stage_job_config.json')
    utility.create_filepath(stage_level_jobconfig_path)

    utility.save_json(stage_level_jobconfig_path, stage_level_jobconfig)
    prepare_jobscript_default = utility.locate_script_file(
        'prepare_stagejob.py')
    analyze_jobscript_default = utility.locate_script_file(
        'analyze_stagejob.py')
    shutil.copy(prepare_jobscript_default, stage_jobdir)
    shutil.copy(analyze_jobscript_default, stage_jobdir)

    jobtemplate_path = 'job_templates/chainjob_template.sh'

    chain_job = ChainSubJob(jobtemplate_path, stage_level_jobconfig_path)
    chain_job.script_generator()

    chain_job.run_job()


def main():
    mod = ags.ArgSchemaParser(schema_type=Launch_Config)
    create_optim_job(mod.args)
