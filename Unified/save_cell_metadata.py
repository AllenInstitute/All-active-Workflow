#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 12:57:55 2018

@author: anin
"""
import json
import os
import argparse
import socket
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.biophysical_api import BiophysicalApi
import allensdk.api.queries.rma_api
import shutil
import subprocess

template_model_dict = {'all_active' :491455321,
             'perisomatic' : 329230710}

opt_metric_dict = {'all_active' : ['Feature_avg_Released_AllActive', 
                                      'Explained_variance_Released_AllActive'],
                    'perisomatic' : ['Feature_avg_Peri','Explained_variance_Peri']}
model_dir = {'all_active' :'neuronal_model', 'perisomatic' : 'peri_model'}

def save_cell_info(cell_metadata,metadata_path):
    """ Create and save a dictionary with all the metadata """
        
    with open(metadata_path,'w') as handle:
        json.dump(cell_metadata,handle,indent=4)      
        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--swc_path', required=False, default=None,
                        help='Morphology path for unpublished cells')
    parser.add_argument('--nwb_path', required=False, default=None,
                        help='Ephys path for unpublished cells')
    parser.add_argument('--Species', required=False, default='Homo Sapiens',
                        choices=['Homo Sapiens','Mus musculus'],
                        help='Set Species of the cell in metadata')
    parser.add_argument('--Area', required=False, default='DG',
                        help='Set Area of the cell in metadata')
    parser.add_argument('--Dendrite_type', required=False, default='spiny',
                        help='Set dendrite type of the cell in metadata')
    parser.add_argument('--Layer', required=False, default=None,
                        help='Set Layer of the cell in metadata')
    parser.add_argument('--launch_job', action="store_true", default=False,
                        help='launch the Optimization')
    parser.add_argument('--premium_queue', action="store_true", default=False,
                        help='For premium queue on NERSC')
    
    args = parser.parse_args()
    
    
    
    metadata_keys = ['Cell_id', 'Released_AllActive_id', 'Perisomatic_id',
                     'Dendrite_type','Species','Area','Cre_line','Layer',
                     'Feature_avg_Released_AllActive', 'Explained_variance_Released_AllActive',
                     'Feature_avg_Peri','Explained_variance_Peri','Machine']
    
    #cell_id_str = raw_input('Enter the cell_id : ')
    cell_id_str = os.path.basename(os.getcwd())
    metadata_path = 'cell_metadata_' + cell_id_str +'.json'
    
    if not os.path.exists(metadata_path):
        cell_info_dict = {meta_field : '' for meta_field in metadata_keys}
        cell_info_dict['Cell_id'] = cell_id_str
        
        cell_id = int(float(cell_info_dict['Cell_id']))
        ctc = CellTypesCache(manifest_file='cell_types/manifest.json')
        cells = ctc.get_cells(require_reconstruction=True)
        metadata_list = list(filter(lambda x: x['id'] == cell_id, cells))
        set_cell_metadata = False
                
        # download the ephys data and sweep metadata
        if not args.nwb_path:
            ctc.get_ephys_data(cell_id)
        else:
            src_nwb = args.nwb_path
            dest_nwb = 'cell_types/'
            shutil.copy(src_nwb, dest_nwb)
            set_cell_metadata = True
        
        # download morphology
        if not args.swc_path:
            ctc.get_reconstruction(cell_id)
        else:
            src_swc = args.swc_path
            dest_swc = 'cell_types/'
            shutil.copy(src_swc, dest_swc)
        
        # download all-active model
        bp = BiophysicalApi()
        model_list = bp.get_neuronal_models(cell_id)
        model_dict = {key : '' for key in template_model_dict.keys()}
        
        for model_type,template_id in template_model_dict.items():
            for model_meta in model_list:
                if model_meta['neuronal_model_template_id'] == template_id:
                    model_dict[model_type] = model_meta['id']
        
        if metadata_list:
            metadata_cell =  metadata_list[0]
            cell_info_dict['Released_AllActive_id'] = str(model_dict['all_active'])
            cell_info_dict['Perisomatic_id'] = str(model_dict['perisomatic'])
            cell_info_dict['Dendrite_type'] = metadata_cell['dendrite_type']
            cell_info_dict['Species'] = metadata_cell['species']
            cell_info_dict['Area'] = metadata_cell['structure_area_abbrev']
            cell_info_dict['Cre_line'] = metadata_cell['transgenic_line']
            cell_info_dict['Layer'] = metadata_cell['structure_layer_name']

        
        api = allensdk.api.queries.rma_api.RmaApi() # Get the model metadata
        
        for model_type, model_id in model_dict.items():
            if model_id != '':
                bp.cache_data(model_id,working_directory = model_dir[model_type])
                model_metadata = api.model_query("NeuronalModelRun", 
                                                 criteria="[neuronal_model_id$eq%s]"%model_id) 
                model_metadata_select = [model_metadata[0]['rheobase_feature_average'],
                                         model_metadata[0]['explained_variance_ratio']]
                
                for i,metadata_key in enumerate(opt_metric_dict[model_type]):
                    cell_info_dict[metadata_key] = model_metadata_select[i]
        

        cell_info_dict['Machine'] =  socket.gethostname()
        
        if 'cori' in cell_info_dict['Machine'] and args.premium_queue:
            with open('nersc_queue.txt', 'a') as handle:
                handle.write('premium')
            
        
        # Set Metadata (Layer,Area,Species,Dendrite_type) from argparse
        if set_cell_metadata:
            for arg in vars(args):
                if arg in cell_info_dict.keys() and getattr(args, arg):
                    cell_info_dict[arg] = getattr(args, arg)
        

        save_cell_info(cell_info_dict,metadata_path)
    
    if args.launch_job:    
        try:
            subprocess.call(['./start_job.sh'])
        except:
            subprocess.call(['./start_job_bbp.sh'])

# Run this file in the parent directory (cell_id) to save the metadata
    
if __name__ == '__main__':
    main()
    