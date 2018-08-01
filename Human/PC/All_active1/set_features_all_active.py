#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 14:13:54 2018

@author: anin
"""

import json

feature_set = {
        'somatic_features' : [
#                              'voltage_base',
                              'steady_state_voltage',
#                              'voltage_deflection_vb_ssse',
#                              'sag_amplitude',
                              'mean_frequency',
                              'time_to_first_spike',
#                              'time_to_last_spike',
                              'AP_amplitude_from_voltagebase',
#                              'min_voltage_between_spikes',
                              'ISI_log_slope',
                              'AP_width',
#                              'AP_duration_half_width',
                              'adaptation_index2',
                              'AHP_depth',
                              'check_AISInitiation'],
#       'dendritic_features' :['AP_amplitude_from_voltagebase',
#                              'AP_width']
            }


def main():
    with open('feature_set.json','w') as feature_file:
        json.dump(feature_set,feature_file,indent=4)


if __name__ == '__main__':
    main()
