#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 23:36:54 2018

@author: anirbannandi
"""

"""Script to get feature from abi traces"""

# pylint: disable=R0914, F0401, R0912

import os
import json
import numpy as np
import math
import collections
import efel



def all_features_path(cell_map, train_protocols_path):
    
    """Get feature values"""
    cell_name = cell_map.keys()[0]
    ephys_location = cell_map[cell_name]['ephys']
#    v_init_model = cell_map[cell_name]['v_init']
    stim_map = get_stim_map(os.path.join(ephys_location, 'StimMapReps.csv'))
    feature_set_map = get_feature_set_map(cell_map[cell_name]['feature_set_map'])
    stim_features = feature_set_map['somatic_features']
    features_meanstd = collections.defaultdict(
            lambda: collections.defaultdict(
                lambda: collections.defaultdict(dict)))
    all_features_json_filename = 'config/'+ cell_name +'/all_features.json'
    trained_features_json_filename = 'config/'+ cell_name +'/trained_features.json'
    untrained_features_json_filename = 'config/'+ cell_name +'/untrained_features.json'
    for stim_name,stim_params in stim_map.items():
                
                
#        print "\n### Getting features from %s of cell %s ###\n" \
#            % (stim_name, cell_name)
                
    
        sweeps = []
        for sweep_filename in stim_map[stim_name]['stimuli'][0]['sweep_filenames']:
            sweep_fullpath = os.path.join(
                ephys_location,
                sweep_filename)
    
            data = np.loadtxt(sweep_fullpath)
            time = data[:, 0]
            voltage = data[:, 1]
#            v_init_cell = voltage[0]
#            v_init_correction = v_init_cell - v_init_model 
            
            # Correct LJP
    #        voltage = voltage - specs['junctionpotential']
            time = time
    
            # Prepare sweep for eFEL
            sweep = {}
            sweep['T'] = time
            sweep['V'] = voltage
            sweep['stim_start'] = [stim_map[stim_name]['stimuli'][0]['delay']]
            sweep['stim_end'] = [stim_map[stim_name]['stimuli'][0]['stim_end']]
            sweep['T;location_AIS'] = time
            sweep['V;location_AIS'] = voltage
            sweep['stim_start;location_AIS'] = [stim_map[stim_name]['stimuli'][0]['delay']]
            sweep['stim_end;location_AIS'] = [stim_map[stim_name]['stimuli'][0]['stim_end']]

            sweeps.append(sweep)
    
        # Do the actual feature extraction
        feature_results = efel.getFeatureValues(sweeps, stim_features)
        
        temp_spike_count = 0
        for feature_temp_list in feature_results:
            if feature_temp_list['mean_frequency']:
                temp_spike_count += feature_temp_list['mean_frequency'][0]
        
        if temp_spike_count == 0: #only select spiking protocols
            continue
        for feature_name in stim_features:
            # For one feature, a list with values for every repetition
            feature_values = [np.mean(trace_dict[feature_name])
                              for trace_dict in feature_results
                              if trace_dict[feature_name] is not None]
            if len(feature_values) == 0:
               
                continue
            elif len(feature_values) == 1:
                mean = feature_values[0]
                std = 0.05 * abs(mean)
            elif len(feature_values) > 1:
                mean = np.mean(feature_values)
                std = np.std(feature_values)
            
            if std== 0 and len(feature_values) != 1:
               std = 0.05 * abs(mean)/math.sqrt(len(feature_values)) 
            
            if math.isnan(mean) or math.isnan(std):
                continue
            if mean == 0:
                std = 0.05
#            if feature_name in ['voltage_base', 'steady_state_voltage']:
#                    mean -= v_init_correction
            features_meanstd[stim_name]['soma'][
                feature_name] = [mean , std]


    save_json(features_meanstd, all_features_json_filename)
    train_protocols = load_json(train_protocols_path).keys()
    trained_features = dict()
    untrained_features = dict()    
    for key in features_meanstd:
        if key in train_protocols:
            trained_features[key] = features_meanstd[key]
        else:
            untrained_features[key] = features_meanstd[key]
    save_json(trained_features, trained_features_json_filename)
    save_json(untrained_features, untrained_features_json_filename)
    return all_features_json_filename, trained_features_json_filename,\
                    untrained_features_json_filename

def load_json(filename):
    """Load json file"""

    with open(filename) as file_h:
        return json.load(file_h)


def save_json(content, filename):
    """Load json file"""

    with open(filename, 'w') as file_h:
        return json.dump(
            content,
            file_h,
            sort_keys=True,
            indent=4,
            separators=(
                ',',
                ': '))




def get_stim_map(stim_map_filename):
    """Get stim map"""

    stim_map = collections.defaultdict(dict)

    with open(stim_map_filename, 'r') as stim_map_file:
        stim_map_content = stim_map_file.read()

    for line in stim_map_content.split('\n')[1:-1]:
        if line is not '':
            stim_name, stim_type, holding_current, amplitude_start, amplitude_end, \
                stim_start, stim_end, duration, sweeps = line.split(',')
            iter_dict= dict()
            iter_dict['type'] = stim_type.strip()
            iter_dict['hypamp'] = 1e9 * float(holding_current)
            iter_dict['amp'] = 1e9 * float(amplitude_start)
            iter_dict['amp_end'] = 1e9 * float(amplitude_end)
            iter_dict['delay'] = float(stim_start)
            iter_dict['duration'] = float(stim_end) - float(stim_start)
            iter_dict['stim_end'] = float(stim_end)
            iter_dict['totduration'] = float(duration)
            iter_dict['sweep_filenames'] = [
                x.strip() for x in sweeps.split('|')]
            
            iter_list = [iter_dict]
            stim_map[stim_name]['stimuli'] = iter_list
        
    return stim_map


def get_specs(specs_filename):
    """Get specs"""

    specs = {}
    with open(specs_filename, 'r') as specs_file:
        specs_content = specs_file.read()

    for line in specs_content.split('\n')[:-1]:
        var_name, var_value = line.split('=')
        specs[var_name] = float(var_value)

    return specs


def get_feature_set_map(feature_map_filename):
    """Get feature set map"""

    with open(feature_map_filename, 'r') as feature_map_file:
        feature_map = json.load(feature_map_file)

    return feature_map
