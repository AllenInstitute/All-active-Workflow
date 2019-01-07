#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 10:26:26 2018

@author: anin
"""

import json

feature_set = {
        "subthreshold_features": [
                    "voltage_base", 
                    "steady_state_voltage", 
                    "voltage_deflection_vb_ssse", 
                    "decay_time_constant_after_stim",
            	    "sag_amplitude",
                    "Spikecount"
                    ]
            }


def main():
    with open('feature_set.json','w') as feature_file:
        json.dump(feature_set,feature_file,indent=4)   
        
        
if __name__ == '__main__':
    main()