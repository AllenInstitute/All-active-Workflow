#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:13:54 2018

@author: anin
"""

import json

feature_set = {
        'somatic_features' : [
                              'voltage_base',
                              'steady_state_voltage',
                              'mean_frequency',
                              'time_to_first_spike',
                              'AP_height',
                              'doublet_ISI',
                              'ISI_CV',
                              'AP_width',
                              'adaptation_index2',
                              'AHP_depth_abs',
                              'AHP_depth_abs_slow',
                              'check_AISInitiation'],
#       'dendritic_features' :['AP_amplitude_from_voltagebase',
#                              'AP_width']
            }


def main():
    with open('feature_set.json','w') as feature_file:
        json.dump(feature_set,feature_file,indent=4)


if __name__ == '__main__':
    main()
