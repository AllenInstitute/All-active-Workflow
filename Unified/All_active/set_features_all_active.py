#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:13:54 2018

@author: anin
"""

import json
import os
import glob

parent_dir = os.path.abspath(os.path.join('.', os.pardir))
path_to_cell_metadata = glob.glob(parent_dir+'/*.json')[0]        
with open(path_to_cell_metadata,'r') as metadata:
        cell_metadata = json.load(metadata)

feature_set = {
        'somatic_features' : [
                              'voltage_base',
                              'steady_state_voltage',
                              'mean_frequency',
                              'time_to_first_spike',
                              'AP_amplitude_from_voltagebase',
                              'ISI_CV',
                              'AP_width',
                              'adaptation_index2',
                              'AHP_depth',
                              'depol_block',
                              'check_AISInitiation']
        
#       'dendritic_features' :['AP_amplitude_from_voltagebase',
#                              'AP_width']
            }

if cell_metadata['Dendrite_type'] == 'aspiny':
  feature_set['somatic_features'].remove('check_AISInitiation')


def main():
    with open('feature_set.json','w') as feature_file:
        json.dump(feature_set,feature_file,indent=4)


if __name__ == '__main__':
    main()
