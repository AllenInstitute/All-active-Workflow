import os
import socket
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.biophysical_api import BiophysicalApi
import allensdk.api.queries.rma_api
import shutil
import logging
from ateamopt.utils import utility
import numpy as np
import glob

try:
    from ateam.data import lims
except:
    pass

logger = logging.getLogger(__name__)

def save_cell_metadata(**cell_metadata):
    cell_id = cell_metadata["cell_id"]
    metadata_filename = 'cell_metadata_%s.json' % cell_id

    # TODO: maybe move this to launch_optimjob?
    machine_name = socket.gethostname()
    machine_name = 'aws' if machine_name == 'master' else machine_name
    cell_metadata['machine'] = machine_name

    data_source = cell_metadata["data_source"]
    if data_source=="web":
        cell_metadata.update(get_data_web(cell_id))
    elif data_source=="lims":
        cell_metadata.update(get_data_lims(cell_id))

    cell_metadata.update(cell_props(cell_id))
    cell_metadata.update(model_props(cell_id))

    utility.save_json(metadata_filename, cell_metadata)
        
    return cell_metadata, metadata_filename

def get_data_lims(cell_id):
    lr = lims.LimsReader()
    nwb_path = lr.get_nwb_path_from_lims(cell_id)
    swc_path = lr.get_swc_path_from_lims(cell_id)
    return {"nwb_path": nwb_path, "swc_path": swc_path}


def get_data_web(cell_id):
    metadata = {}
    ctc.get_ephys_data(cell_id)
    nwb_path = os.path.abspath(
            utility.get_filepath_for_exten('.nwb')[0])
    ctc.get_reconstruction(cell_id)
    swc_path = os.path.abspath(
            utility.get_filepath_for_exten('.swc')[0])
    return {"nwb_path": nwb_path, "swc_path": swc_path}

# TODO: pull from lims, get_cells_df
def cell_props(cell_id):
    cell_metadata = {}
    ctc = CellTypesCache(manifest_file='cell_types/manifest.json')
    cells = ctc.get_cells()
    metadata_list = list(filter(lambda x: x['id'] == cell_id, cells))
    if metadata_list:
        cell_metadata.update(metadata_list[0])
    return cell_metadata


def model_props(cell_id):
    cell_metadata = {}
    # download existing model (if exists)
    bp = BiophysicalApi()
    try:
        model_list = bp.get_neuronal_models(cell_id)
        model_dict = {key: '' for key in template_model_dict.keys()}

        for model_type, template_id in template_model_dict.items():
            for model_meta in model_list:
                if model_meta['neuronal_model_template_id'] == template_id:
                    model_dict[model_type] = model_meta['id']

    except:
        logger.debug('No biophysical model available')

    api = allensdk.api.queries.rma_api.RmaApi()  # Get the model metadata

    for model_type, model_id in model_dict.items():
        if model_id != '':
            bp.cache_data(
                model_id, working_directory=model_dir[model_type])
            model_metadata = api.model_query("NeuronalModelRun",
                                                criteria="[neuronal_model_id$eq%s]" % model_id)
            model_metadata_select = [model_metadata[0]['rheobase_feature_average'],
                                        model_metadata[0]['explained_variance_ratio']]

            model_path = glob.glob('%s/*fit*.json' %
                                    model_dir[model_type])[0]
            for file_ in glob.glob('%s/*' % model_dir[model_type]):
                if file_ != model_path:
                    try:
                        os.remove(file_)
                    except:
                        shutil.rmtree(file_, ignore_errors=True)

            cell_metadata['model_path_{}'.format(
                model_type)] = os.path.abspath(model_path)
            cell_metadata['{}_id'.format(
                model_dir[model_type])] = str(model_id)
            for i, metadata_key in enumerate(opt_metric_dict[model_type]):
                cell_metadata[metadata_key] = model_metadata_select[i]
    return cell_metadata


# Formatting specs
template_model_dict = {'all_active': 491455321,
                        'perisomatic': 329230710}

opt_metric_dict = {'all_active': ['Feature_avg_Released_AllActive',
                                    'Explained_variance_Released_AllActive'],
                    'perisomatic': ['Feature_avg_Peri', 'Explained_variance_Peri']}

model_dir = {'all_active': 'released_aa_model',
                'perisomatic': 'peri_model'}
