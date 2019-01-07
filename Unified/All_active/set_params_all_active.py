#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 13:54:54 2018

@author: anin
"""

import json
import os

path_to_cell_metadata = os.path.abspath(os.path.join('.', os.pardir)) + '/cell_metadata.json'        
with open(path_to_cell_metadata,'r') as metadata:
    cell_metadata = json.load(metadata)

species = cell_metadata['Species']
dendrite_type = cell_metadata['Dendrite_type']

if species == 'Mouse' and dendrite_type == 'spiny':
    all_params = json.load(open('param_set1.json','r'))
elif species == 'Mouse' and dendrite_type == 'aspiny':
    all_params = json.load(open('param_set2.json','r'))
elif species == 'Human' and dendrite_type == 'spiny':
    all_params = json.load(open('param_set3.json','r'))
elif species == 'Human' and dendrite_type == 'aspiny':
    all_params = json.load(open('param_set4.json','r'))
    
    
def main():
    with open('all_param_bounds.json','w') as bound_file:
        json.dump(all_params,bound_file,indent=4)   
        
        
if __name__ == '__main__':
    main()