#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 12:57:55 2018

@author: anin
"""

import json



def save_cell_info(*args):
    """Create and save a dictionary with all the metadata"""
    
    cell_keys = ['Cell_id', 'Model_id', 'Species', 'Cre_line', 'Area', 'Layer']
    cell_info_dict = {}
    for i,arg in enumerate(args):
        cell_info_dict[cell_keys[i]] = arg
        
    with open('cell_metadata.json','w') as cell_info:
        json.dump(cell_info_dict,cell_info,indent=4)      
        
def main():
    cell_id = raw_input('Enter the cell_id : ')
    model_id = raw_input('Enter the model_id : ')
    species = raw_input('Enter the species : ')
    Cre_line = raw_input('Enter the Cre-line : ')
    Area = raw_input('Enter the Area : ')
    layer = raw_input('Enter the layer number: ')
    
    save_cell_info(cell_id, model_id, species, Cre_line, Area, layer)
    
    
    
if __name__ == '__main__':
    main()
    