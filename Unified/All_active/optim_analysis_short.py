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
import glob
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
    
parent_dir = os.path.abspath(os.path.join('.', os.pardir))
path_to_cell_metadata = glob.glob(parent_dir+'/*.json')[0]
with open(path_to_cell_metadata,'r') as metadata:
        cell_metadata = json.load(metadata)
    

#####################################################################
 

# Extracting the parameters
        
def get_hof_params_responses(opt_gen,opt_train_feat,opt_untrain_feat,
                             checkpoint_file, cp_dir, hof_index, responses_filename):

    checkpoint = pickle.load(open(checkpoint_file, "r"))
    
    best_individual = [checkpoint['halloffame'][hof_index]]
        
    checkpoint = None
    
    logger.debug('Saving the hall of fame parameters')  
    
    
    ## Save hall-of-fame parameters for all seeds

    checkpoint_dir  = cp_dir + '/*.pkl'
    cp_list = glob.glob(checkpoint_dir)
    
    hof_params = []    
    seed_indices = []
    
    for i,cp_file in enumerate(cp_list):
        hof_i = pickle.load(open(cp_file, "r"))['halloffame']
        hof_params.extend(hof_i)
        seed = [cp_file.split('/')[1].split('.')[0]]*len(hof_i)
        seed_indices.extend(seed)

    
    
    plot_diversity_params_path = 'analysis_params/plot_diversity_params.pkl'
    hof_responses_filename = 'analysis_params/hof_response_all.pkl'    
    obj_list_gen_path = 'analysis_params/hof_obj.pkl'    
    obj_list_train_path = 'analysis_params/hof_obj_train.pkl'
    obj_list_untrain_path = 'analysis_params/hof_obj_untrain.pkl'
    seed_indices_path = 'analysis_params/seed_indices.pkl'
    
    if not os.path.exists(os.path.dirname(hof_responses_filename)):
        try:
            os.makedirs(os.path.dirname(hof_responses_filename))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
                
                
    # Calculate responses for the hall-of-fame parameters
    
    if not os.path.exists(hof_responses_filename):
        
        logger.debug('Calculating Hall of Fame Responses')
        
        response_list = opt_gen.toolbox.map(opt_gen.toolbox.save_sim_response,hof_params)
        pickle.dump(response_list,open(hof_responses_filename,'wb'))
        
    else:
        
        logger.debug('Retrieving Hall of Fame Responses')
        response_list = pickle.load(open(hof_responses_filename,'r'))
         
    logger.debug('Calculating scores for Optimized Responses')
    
    obj_list_gen = opt_gen.toolbox.map(opt_gen.toolbox.evaluate_response,response_list)
    obj_list_train = opt_train_feat.toolbox.map(opt_train_feat.toolbox.evaluate_response,response_list)
    obj_list_train_total = [np.sum(obj_dict_train.values()) for obj_dict_train in obj_list_train]
    obj_list_untrain = opt_untrain_feat.toolbox.map(opt_untrain_feat.toolbox.evaluate_response,response_list)
    
    
    arrange_lists = [hof_params, response_list, obj_list_gen, obj_list_train,
                                         obj_list_untrain,seed_indices]
    f_name_lists = [plot_diversity_params_path,hof_responses_filename,obj_list_gen_path,
                    obj_list_train_path,obj_list_untrain_path,seed_indices_path]
    
    
    for kk,arrange_list in enumerate(arrange_lists):
        arrange_list = [x for _,x in sorted(zip(obj_list_train_total,arrange_list))]
        arrange_lists[kk] = arrange_list
            
    arrange_lists[0] = {'optimized_individual': best_individual,
                     'hof_list' : arrange_lists[0]}
    
    for f_name,data in zip(f_name_lists,arrange_lists):
        with open(f_name,'wb') as handle:
            pickle.dump(data,handle)            
            
    with open(responses_filename, 'w') as fd:
        pickle.dump(response_list[0], fd)
    
    
        
    
#####################################################################

# GA evolution

def get_GA_evolution_params(checkpoint_file):
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
    
    


#############################################################################

## Calculating released responses


def get_release_responses(opt,opt_release,response_release_filename):
    

    
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
    
    
