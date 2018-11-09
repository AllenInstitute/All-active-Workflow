#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 12:57:55 2018

@author: anin
"""

import json
import os

metadata_path = 'cell_metadata.json'

def save_cell_info(*args):
    """ Create and save a dictionary with all the metadata """
    
    cell_keys = ['Cell_id', 'Model_id', 'Perisomatic_id','Species',
                  'Cre_line','Area', 'Layer',
                  'Dendrite_type','Feature_avg','Explained_variance',
                     'Machine']
    cell_info_dict = {}
    for i,arg in enumerate(args):
        cell_info_dict[cell_keys[i]] = arg
        
    with open(metadata_path,'w') as cell_info:
        json.dump(cell_info_dict,cell_info,indent=4)      
        
def main():
    if not os.path.exists(metadata_path):
        cell_id = raw_input('Enter the cell_id : ')
        model_id = raw_input('Enter the model_id : ')
        perisomatic_id = raw_input('Enter the perisomatic model_id : ')
        species = raw_input('Enter the species : ')
        cre_line = raw_input('Enter the Cre-line : ')
        area = raw_input('Enter the Area : ')
        layer = raw_input('Enter the layer number : ')
        dendrite_type = raw_input('Enter dendrite type : ')
        feat_avg_released = raw_input('Enter the feature average of the released perisomatic model : ')
        exp_var_released = raw_input('Enter the explained variance of the released perisomatic model : ')
        machine =  raw_input('Enter the resource facility name : ')
        save_cell_info(cell_id, model_id,perisomatic_id,species, cre_line, area, 
                         layer, dendrite_type,
                         feat_avg_released,exp_var_released,
                         machine)
        
    
# Run this file in the parent directory (cell_id) to save the metadata
    
if __name__ == '__main__':
    main()
    