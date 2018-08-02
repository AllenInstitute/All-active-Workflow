#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 16:14:10 2018

@author: anin
"""

from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.biophysical_api import BiophysicalApi
import json

with open('cell_metadata.json','r') as info_dict:
    cell_info = json.load(info_dict)



def main():
    cell_id = int(float(cell_info['Cell_id']))
    ctc = CellTypesCache(manifest_file='cell_types/manifest.json')
    
    # download the ephys data and sweep metadata
    ctc.get_ephys_data(cell_id)
    
    # download morphology
    ctc.get_reconstruction(cell_id)
    
    # download all-active model
    if cell_info['Model_id'] != '':
        model_id = int(float(cell_info['Model_id']))
        bp = BiophysicalApi()
        bp.cache_data(model_id,working_directory='neuronal_model')
        
        

if __name__ == '__main__':
    main()