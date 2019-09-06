import os
import socket
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.biophysical_api import BiophysicalApi
import allensdk.api.queries.rma_api
import shutil
import logging
from ateamopt.utils import utility
from ateam.data import lims
import glob

logger = logging.getLogger(__name__)

class Allactive_Optim(object):

    def __init__(self,cell_id):
        self.cell_id = int(float(cell_id))


    def get_ephys_morphology(self,ctc,metadata,**kwargs):
        
        # download/retrieve ephys data and sweep metadata
        if not kwargs.get('nwb_path'):
            if 'hpc' in metadata['machine']: 
                lr = lims.LimsReader()
                logger.debug('Retrieving files from network')
                nwb_path = lr.get_nwb_path_from_lims(int(self.cell_id),get_sdk_version=True)

            else:
                ctc.get_ephys_data(self.cell_id)
                nwb_path = os.path.abspath(utility.get_filepath_for_exten('.nwb')[0])
            metadata['nwb_path'] = nwb_path

        # download/retrieve morphology
        if not kwargs.get('swc_path'):
            if 'hpc' in metadata['machine']:
                lr = lims.LimsReader()
                logger.debug('Retrieving files from network')
                swc_path = lr.get_swc_path_from_lims(int(self.cell_id))

            else:
                ctc.get_reconstruction(self.cell_id)
                swc_path =  os.path.abspath(utility.get_filepath_for_exten('.swc')[0])
            metadata['swc_path'] = swc_path
        return metadata

    def save_cell_metadata(self, **props):

#        metadata_keys = ['Cell_id', 'Released_AllActive_id', 'Perisomatic_id',
#                 'Dendrite_type','Species','Area','Cre_line','Layer',
#                 'Feature_avg_Released_AllActive', 'Explained_variance_Released_AllActive',
#                 'Feature_avg_Peri','Explained_variance_Peri','machine','Axon_type']

        cell_metadata = {prop_key:props.get(prop_key) for prop_key in ['cell_id','ateamopt_tag','bluepyopt_tag']}
        
        
        template_model_dict = {'all_active' :491455321,
             'perisomatic' : 329230710}

        opt_metric_dict = {'all_active' : ['Feature_avg_Released_AllActive',
                                      'Explained_variance_Released_AllActive'],
                    'perisomatic' : ['Feature_avg_Peri','Explained_variance_Peri']}
        model_dir = {'all_active' :'released_aa_model', 'perisomatic' : 'peri_model'}

        cell_id = self.cell_id

        metadata_filename='cell_metadata_%s.json'%cell_id
        machine_name = socket.gethostname()
        machine_name = 'aws' if machine_name == 'master' else machine_name
        cell_metadata['machine']= machine_name

        ctc = CellTypesCache(manifest_file='cell_types/manifest.json')
        cells = ctc.get_cells()
        metadata_list = list(filter(lambda x: x['id'] == cell_id, cells))
        if metadata_list:
            cell_metadata.update(metadata_list[0])
        cell_metadata = self.get_ephys_morphology(ctc,cell_metadata,**props)

        

        # download existing model (if exists)
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


        
#                set_cell_metadata = False
#                metadata_cell =  metadata_list[0]
#                cell_metadata['Released_AllActive_id'] = str(model_dict['all_active'])
#                cell_metadata['Perisomatic_id'] = str(model_dict['perisomatic'])
#                cell_metadata['Dendrite_type'] = metadata_cell['dendrite_type']
#                cell_metadata['Species'] = metadata_cell['species']
#                cell_metadata['Area'] = metadata_cell['structure_area_abbrev']
#                cell_metadata['Cre_line'] = metadata_cell['transgenic_line']
#                cell_metadata['Layer'] = metadata_cell['structure_layer_name']
                

        api = allensdk.api.queries.rma_api.RmaApi() # Get the model metadata

        for model_type, model_id in model_dict.items():
            if model_id != '':
                bp.cache_data(model_id,working_directory = model_dir[model_type])
                model_metadata = api.model_query("NeuronalModelRun",
                                                 criteria="[neuronal_model_id$eq%s]"%model_id)
                model_metadata_select = [model_metadata[0]['rheobase_feature_average'],
                                         model_metadata[0]['explained_variance_ratio']]
                
                model_path = glob.glob('%s/*fit*.json'%model_dir[model_type])[0]
                for file_ in glob.glob('%s/*'%model_dir[model_type]):
                    if file_ != model_path:
                        try:
                            os.remove(file_)
                        except:
                            shutil.rmtree(file_, ignore_errors=True)

                cell_metadata['model_path_{}'.format(model_type)] = os.path.abspath(model_path)
                cell_metadata['{}_id'.format(model_dir[model_type])] = str(model_id)
                for i,metadata_key in enumerate(opt_metric_dict[model_type]):
                    cell_metadata[metadata_key] = model_metadata_select[i]

        utility.save_json(metadata_filename,cell_metadata)
        
        return cell_metadata,metadata_filename

