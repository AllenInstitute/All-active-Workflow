
import os
import socket
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.biophysical_api import BiophysicalApi
import allensdk.api.queries.rma_api
import shutil
import logging
from ateamopt.utils import utility
import neurom as nm
from neurom import morphmath as mm
import numpy as np
from neurom.core.types import tree_type_checker, NEURITES


logger = logging.getLogger(__name__)

class Allactive_Optim(object):

    def __init__(self,nwb_path = None, swc_path = None):
        self.cell_id = os.path.basename(os.getcwd())
        self._nwb_path = nwb_path
        self._swc_path = swc_path

    @property
    def nwb_path(self):
        return self._nwb_path

    @nwb_path.setter
    def nwb_path(self,path):
        self._nwb_path = path

    @property
    def swc_path(self):
        return self._swc_path

    @swc_path.setter
    def swc_path(self,path):
        self._swc_path = path

    def save_morph_data(self):
        morph_stats_filename = 'morph_stats_%s.json'%self.cell_id

        if not os.path.exists(morph_stats_filename):
            swc_filename = utility.get_filepath_for_exten('.swc')
            swc_filename = [swc_file_ for swc_file_ in swc_filename if 'cell_types' in swc_file_][0]

            #  load a neuron from an SWC file
            nrn = nm.load_neuron(swc_filename)

            morph_stats = {}

            morph_stats['Cell_id'] = self.cell_id
            morph_stats['soma_suface'] = nm.get('soma_surface_areas', nrn)[0]
            morph_stats['soma_radius'] = np.mean(nm.get('soma_radii', nrn))

            # Morph stats
            for nrn_type_ in NEURITES:

                morph_stats['length' + '.' + str(nrn_type_).split('.')[1]] = np.sum(nm.get('segment_lengths',
                                                           nrn, neurite_type=nrn_type_))
                morph_stats['area' + '.' + str(nrn_type_).split('.')[1]] = sum(mm.segment_area(s)
                           for s in nm.iter_segments(nrn,neurite_filter=tree_type_checker(nrn_type_)))
                morph_stats['volume' + '.' + str(nrn_type_).split('.')[1]] = sum(mm.segment_volume(s)
                            for s in nm.iter_segments(nrn,neurite_filter=tree_type_checker(nrn_type_)))
                morph_stats['taper_rate' + '.' + str(nrn_type_).split('.')[1]] = \
                            np.mean([mm.segment_taper_rate(s) for s in nm.iter_segments(nrn, neurite_filter=tree_type_checker(nrn_type_))])

            utility.save_json(morph_stats_filename,morph_stats)


    def save_cell_metadata(self, get_data=True, **props):

        metadata_keys = ['Cell_id', 'Released_AllActive_id', 'Perisomatic_id',
                 'Dendrite_type','Species','Area','Cre_line','Layer',
                 'Feature_avg_Released_AllActive', 'Explained_variance_Released_AllActive',
                 'Feature_avg_Peri','Explained_variance_Peri','Machine']


        template_model_dict = {'all_active' :491455321,
             'perisomatic' : 329230710}

        opt_metric_dict = {'all_active' : ['Feature_avg_Released_AllActive',
                                      'Explained_variance_Released_AllActive'],
                    'perisomatic' : ['Feature_avg_Peri','Explained_variance_Peri']}
        model_dir = {'all_active' :'neuronal_model', 'perisomatic' : 'peri_model'}

        cell_id = self.cell_id

        metadata_filename='cell_metadata_%s.json'%cell_id

        if get_data:
            cell_metadata = {meta_field : '' for meta_field in metadata_keys}
            cell_metadata['Cell_id'] = cell_id
            try:
                cell_id = int(float(cell_metadata['Cell_id']))
            except:
                logger.debug('Cell not part of online product')

            ctc = CellTypesCache(manifest_file='cell_types/manifest.json')
            cells = ctc.get_cells()
            metadata_list = list(filter(lambda x: x['id'] == cell_id, cells))


            # download the ephys data and sweep metadata
            if not self.nwb_path and not 'nwb_path' in props:
                ctc.get_ephys_data(cell_id)
            else:
                src_nwb = props['nwb_path']
                self.nwb_path = src_nwb
                dest_nwb = 'cell_types/'
                shutil.copy(src_nwb, dest_nwb)
                set_cell_metadata = True

            # download morphology
            if not self.swc_path and not 'swc_path' in props:
                ctc.get_reconstruction(cell_id)
            else:
                src_swc = props['swc_path']
                self.swc_path = src_swc
                dest_swc = 'cell_types/'
                shutil.copy(src_swc, dest_swc)

            # download all-active model (if exists)
            bp = BiophysicalApi()
            try:
                model_list = bp.get_neuronal_models(cell_id)
                model_dict = {key : '' for key in template_model_dict.keys()}

                for model_type,template_id in template_model_dict.items():
                    for model_meta in model_list:
                        if model_meta['neuronal_model_template_id'] == template_id:
                            model_dict[model_type] = model_meta['id']

            except:
                logger.debug('No biophysical model available')


            if metadata_list:
                set_cell_metadata = False
                metadata_cell =  metadata_list[0]
                cell_metadata['Released_AllActive_id'] = str(model_dict['all_active'])
                cell_metadata['Perisomatic_id'] = str(model_dict['perisomatic'])
                cell_metadata['Dendrite_type'] = metadata_cell['dendrite_type']
                cell_metadata['Species'] = metadata_cell['species']
                cell_metadata['Area'] = metadata_cell['structure_area_abbrev']
                cell_metadata['Cre_line'] = metadata_cell['transgenic_line']
                cell_metadata['Layer'] = metadata_cell['structure_layer_name']


            api = allensdk.api.queries.rma_api.RmaApi() # Get the model metadata

            for model_type, model_id in model_dict.items():
                if model_id != '':
                    bp.cache_data(model_id,working_directory = model_dir[model_type])
                    model_metadata = api.model_query("NeuronalModelRun",
                                                     criteria="[neuronal_model_id$eq%s]"%model_id)
                    model_metadata_select = [model_metadata[0]['rheobase_feature_average'],
                                             model_metadata[0]['explained_variance_ratio']]

                    for i,metadata_key in enumerate(opt_metric_dict[model_type]):
                        cell_metadata[metadata_key] = model_metadata_select[i]

            machine_name = socket.gethostname()
            cell_metadata['Machine'] =  'aws' if machine_name == 'master' else machine_name
            if 'cori' in cell_metadata['Machine'] and 'premium_queue' in props:
                filename= 'nersc_queue.txt'
                content = 'premium'
                utility.save_file(filename,content)

            # Manually set metadata

            if set_cell_metadata:
                cell_metadata_set= set(cell_metadata)
                props_set = set(props)
                for metadata_key in cell_metadata_set.intersection(props_set):
                    cell_metadata[metadata_key] = props[metadata_key]


            utility.save_json(metadata_filename,cell_metadata)
        else:
            cell_metadata = utility.load_json(metadata_filename)
        return cell_metadata

