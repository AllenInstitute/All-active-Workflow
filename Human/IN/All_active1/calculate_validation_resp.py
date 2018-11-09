#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May  3 13:34:44 2018

@author: anirbannandi
"""

import allensdk
from allensdk.core.cell_types_cache import CellTypesCache
import os
import numpy as np
import bluepyopt.ephys as ephys
import pickle
import bluepyopt as bpopt
import json
import efel
import logging
import collections
import errno

import evaluator_helper


logging.basicConfig(level=logging.DEBUG) 
logger = logging.getLogger()

path_to_cell_metadata = os.path.abspath(os.path.join('.', os.pardir)) + '/cell_metadata.json'        
with open(path_to_cell_metadata,'r') as metadata:
        cell_metadata = json.load(metadata)    
cell_id = cell_metadata['Cell_id']


with open('config_file.json') as json_file:  
    path_data = json.load(json_file)

junction_potential = -14
nwb_path = path_data['ephys']
nwb_file = allensdk.core.nwb_data_set.NwbDataSet(nwb_path)
output_dir = os.getcwd() +'/preprocessed'

### Calculate the response only for Noise and Tri-blip stim

acceptable_stimtypes = [
    'Noise 1',
    'Noise 2',
    'Short Square - Triple',
    
    ]


distinct_id_map = {
    'Noise 1': 'Noise_1',
    'Noise 2': 'Noise_2',
    'Short Square - Triple' : 'Short_Square_Triple',
}


# Get the hall-of-famers  
 
optimized_params_path = 'analysis_params/plot_diversity_params.pkl'
hof = pickle.load(open(optimized_params_path,'r'))['hof']

### load the configuration files for Neuron simulation under Bluepyopt

release_params = path_data['release_params']                
morph_path = path_data['morphology']
protocol_path = path_data['protocols']
mech_path = path_data['mechanism']
feature_path = path_data['features']
param_path = path_data['parameters']
    
evaluator = evaluator_helper.create(protocol_path, feature_path, morph_path, 
                                    param_path, mech_path)

opt = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator)

validation_response_path = 'Validation_Responses/validation_responses.pkl'
    
if not os.path.exists(os.path.dirname(validation_response_path)):
    try:
        os.makedirs(os.path.dirname(validation_response_path))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
spiketimes_model_path = 'Validation_Responses/spiketimes_model.pkl'
spiketimes_exp_path = 'Validation_Responses/spiketimes_exp.pkl'

def Main(hof_index = 0):
    validation_responses = {}
    peaktimes_exp = {}
    peaktimes_model = {}
    
    best_individual = hof[hof_index]
    individual_dict = evaluator.param_dict(best_individual)

    ### Create a dictionary to keep track of which Noise type has already been simulated
    
    sim_run = collections.defaultdict(str)
    
    ## Loop through the stimulus sweeps to save the stim-response from nwb file and 
    ## also calculate the model response for the sweeps
    
    for sweep_number in nwb_file.get_sweep_numbers():
        sweep_data = nwb_file.get_sweep_metadata(sweep_number)
        stim_type = sweep_data['aibs_stimulus_name']
        
        if stim_type in acceptable_stimtypes:
            sweep = nwb_file.get_sweep(sweep_number)
    
            start_idx, stop_idx = sweep['index_range']
    
            stimulus_trace = sweep['stimulus'][start_idx:stop_idx]
            response_trace = sweep['response'][start_idx:stop_idx]
    
            sampling_rate = sweep['sampling_rate']
    
            time = np.arange(0, len(stimulus_trace)) / sampling_rate
            trace_name = '%s_%d' % (
                distinct_id_map[stim_type], sweep_number)
            
            if 'Noise' in stim_type:
                stim_type_sweep = distinct_id_map[stim_type]
            else:
                stim_type_sweep = trace_name
            
            # Check if the response for a particular simulus sweep is already calculated
            
            if os.path.exists(validation_response_path):
                validation_responses = pickle.load(open(validation_response_path, "r"))
                peaktimes_model = pickle.load(open(spiketimes_model_path, "r"))
                peaktimes_exp = pickle.load(open(spiketimes_exp_path, "r"))
                
                if trace_name in validation_responses.keys():
                    logger.debug('Response and Features are already calculated for %s'%trace_name)
                    sim_run[stim_type_sweep] = 'Done'
                
            
            time *= 1000.0
            response_trace *= 1000
            stimulus_trace *= 1e9 
            
            time_end = time[-1]
            response_end = response_trace[-1]
            stimulus_end = stimulus_trace[-1]
            
            # downsampling
            time = time[::5]
            response_trace = response_trace[::5] + junction_potential
            stimulus_trace =stimulus_trace[::5]
            
            if time_end != time[-1]:
                time = np.append(time,time_end)
                response_trace = np.append(response_trace,response_end)
                stimulus_trace = np.append(stimulus_trace,stimulus_end)
        
            nonzero_stim_indices = np.where(stimulus_trace != 0)[0]
            stim_start = time[nonzero_stim_indices[0]]
            stim_stop = time[nonzero_stim_indices[-1]]
            
            if trace_name not in peaktimes_exp.keys():
        
                response_trace_short_filename = '%s.%s' % (trace_name, 'txt')
        
                response_trace_filename = os.path.join(
                    output_dir, response_trace_short_filename)
                
                
                with open(response_trace_filename, 'w') as response_trace_file:
                    np.savetxt(response_trace_file,
                                  np.transpose([time, stimulus_trace, response_trace]))
                 
                
                trace = {}
                trace['T'] = time
                trace['V'] = response_trace
                trace['stim_start'] = [stim_start]
                trace['stim_end'] = [stim_stop]
#                if trace_name not in peaktimes_exp.keys():
                peaktimes_exp[trace_name] = efel.getFeatureValues(
                    [trace],
                    ['peak_time'])[0]['peak_time']
                
                pickle.dump(peaktimes_exp, open(spiketimes_exp_path, 'w'))
                logger.debug('Getting response features from %s'%trace_name)
            
            else:
                
                logger.debug('Response features are already calculated for %s'%trace_name)
            
            if sim_run[stim_type_sweep] != '':
                continue
            
            logger.debug('Calculating response for %s'%trace_name)
    
            soma_loc = ephys.locations.NrnSeclistCompLocation(
            name='soma',
            seclist_name='somatic',
            sec_index=0,
            comp_x=0.5)
            
            validation_recording = ephys.recordings.CompRecording(
            name='validation.soma.v',
            location=soma_loc,
            variable='v')
                        
            validation_time = time
            validation_current = stimulus_trace
            
            validation_stimulus = ephys.stimuli.NrnCurrentPlayStimulus(
                current_points=validation_current,
                time_points=validation_time,
                location=soma_loc)
            
            validation_protocol = ephys.protocols.SweepProtocol(
                'validation',
                [validation_stimulus],
                [validation_recording])
            
            validation_responses[trace_name] = opt.evaluator.run_protocols(
                    [validation_protocol],
                    param_values=individual_dict)
            
            resp_time = validation_responses[
                trace_name]['validation.soma.v']['time']
            resp_voltage = validation_responses[trace_name][
                'validation.soma.v']['voltage']
            
            logger.debug('Finished calculation for %s'%trace_name)
            
            trace = {}
            trace['T'] = resp_time
            trace['V'] = resp_voltage
            trace['stim_start'] = [stim_start]
            trace['stim_end'] = [stim_stop]
            peaktimes_model[trace_name] = efel.getFeatureValues(
                [trace],
                ['peak_time'])[0]['peak_time']
            
            # Save the responses
            
            pickle.dump(peaktimes_model, open(spiketimes_model_path, 'w'))
            pickle.dump(validation_responses, open(validation_response_path, 'w'))
            
            if 'Noise' in stim_type:
                sim_run[stim_type_sweep] = 'Done'
        

if __name__ == '__main__':
    Main()