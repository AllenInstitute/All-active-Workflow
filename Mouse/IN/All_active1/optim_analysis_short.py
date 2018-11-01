#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 14:37:29 2018

@author: anin
"""

#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 09:54:58 2018

@author: anin
"""

import os
import json
import pickle
import numpy as np

import errno
import logging
import bluepyopt.ephys as ephys



logger = logging.getLogger(__name__)

section_map_inv = {'somatic':'soma', 'axonal':'axon', 'apical':'apic',
               'basal':'dend', 'all':'all'}

with open('config_file.json') as json_file:  
    data = json.load(json_file)


release_params = data['release_params']
release_params_original = {} #to compare against the model released on the website
        
fit_json_path = data['fit_json']
morph_path = data['morphology']
protocol_path = data['protocols']
mech_path = data['mechanism']
feature_path = data['features']
param_path = data['parameters']
all_protocol_path = data['all_protocols']
original_release_param = data['original_parameters']


with open(protocol_path) as json_file:  
    train_protocols = json.load(json_file)
    
with open(param_path) as json_file:  
    parameters_optim = json.load(json_file)
    
path_to_cell_metadata = os.path.abspath(os.path.join('.', os.pardir)) + '/cell_metadata.json'        
with open(path_to_cell_metadata,'r') as metadata:
        cell_metadata = json.load(metadata)
    

#####################################################################
 

# Extracting the parameters
        
def plot_diversity(opt, checkpoint_file, param_names,hof_index):

    checkpoint = pickle.load(open(checkpoint_file, "r"))
    
    optimized_individual = [checkpoint['halloffame'][hof_index]]
    hof_list = checkpoint['halloffame'][1:]
    hof = checkpoint['halloffame']
    
    plot_diversity_params = {'optimized_individual': optimized_individual,
                             'hof_list' : hof_list,
                             'hof' : hof}
    checkpoint = None
    plot_diversity_params_path = 'analysis_params/plot_diversity_params.pkl'
    
    if not os.path.exists(os.path.dirname(plot_diversity_params_path)):
        try:
            os.makedirs(os.path.dirname(plot_diversity_params_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    
    with open(plot_diversity_params_path,'wb') as handle:
        pickle.dump(plot_diversity_params,handle)
    
    logger.debug('Saving the plot diversity parameters')    
    
#####################################################################

# GA evolution

def plot_GA_evolution(checkpoint_file):
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    
    log = checkpoint['logbook']
    gen_numbers = log.select('gen')
    mean = np.array(log.select('avg'))
    std = np.array(log.select('std'))
    minimum = np.array(log.select('min'))
    
    plot_GA_evolution_params = {'gen_numbers': gen_numbers,
                             'mean' : mean,
                             'std' : std,
                             'minimum' : minimum}
    checkpoint = None
    plot_GA_evolution_params_path = 'analysis_params/plot_GA_evolution_params.pkl'
    
    if not os.path.exists(os.path.dirname(plot_GA_evolution_params_path)):
        try:
            os.makedirs(os.path.dirname(plot_GA_evolution_params_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    
    with open(plot_GA_evolution_params_path,'wb') as handle:
        pickle.dump(plot_GA_evolution_params,handle)
    
    logger.debug('Saving the plot_GA_evolution parameters')    
    
    
#########################################################################


def get_responses(cell_evaluator, individuals, filename):
    responses = []
    response_status = 'Optimized Responses'
        
    if filename and os.path.exists(filename):
        logger.debug('Retrieving %s \n'%response_status)
        with open(filename) as fd:
            return pickle.load(fd)
        
    for i,individual in enumerate(individuals):
        if i == 0:
           logger.debug('Calculating %s \n'%response_status) 
        individual_dict = cell_evaluator.param_dict(individual)
        responses.append(
            cell_evaluator.run_protocols(
                cell_evaluator.fitness_protocols.values(),
                param_values=individual_dict)) 
        
    if filename:
        with open(filename, 'w') as fd:
            pickle.dump(responses, fd)

    return responses



#############################################################################

## Calculating responses


def plot_Response(opt,opt_release,checkpoint_file, responses_filename,
                  response_release_filename,hof_index):
    
    with open(checkpoint_file,'r') as handle:
        checkpoint = pickle.load(handle)
    hof = checkpoint['halloffame']
    checkpoint = None
    get_responses(opt.evaluator, [hof[hof_index]], responses_filename)
    
    if response_release_filename and os.path.exists(response_release_filename):
        responses_release = pickle.load(open(response_release_filename, "r"))
        logger.debug('Retrieving Released Responses')
    else:
        logger.debug('Calculating Released Responses')
        fitness_protocols = opt.evaluator.fitness_protocols
        responses_release = {}
        if original_release_param:
            nrn = ephys.simulators.NrnSimulator()
            
            for protocol in fitness_protocols.values():
                response_release = protocol.run(
                        cell_model=opt_release.evaluator.cell_model,
                        param_values=release_params_original,
                        sim=nrn)
                responses_release.update(response_release)
            
        if response_release_filename:    
            with open(response_release_filename, 'w') as fd:
                pickle.dump(responses_release, fd)
    
    
                    
