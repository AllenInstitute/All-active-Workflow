#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 10:26:26 2018

@author: anin
"""

from ateamopt.utils import utility

feature_set = {
        'features' : ['voltage_base',
                       'steady_state_voltage',
                       'voltage_deflection_vb_ssse',
                       'decay_time_constant_after_stim',
                       'Spikecount']
        }


def main():
    path = 'feature_set_stage0.json'
    utility.save_json(path,feature_set)
    
        
if __name__ == '__main__':
    main()