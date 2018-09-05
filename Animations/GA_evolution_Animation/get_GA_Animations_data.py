#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 11:17:08 2018

@author: anin
"""

import numpy as np
import pickle
import bluepyopt.ephys as ephys
import bluepyopt as bpopt
import json
import math

import evaluator_helper

with open('config_file.json') as json_file:  
    data = json.load(json_file)
release_params = data['release_params']
morph_path = data['morphology']
protocol_path = data['protocols']
mech_path = data['mechanism']
feature_path = data['features']
param_path = data['parameters']
evaluator = evaluator_helper.create(protocol_path, feature_path, morph_path, 
                                    param_path, mech_path)

opt = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator)

checkpoint_file = 'checkpoints/seed1.pkl'
checkpoint = pickle.load(open(checkpoint_file, "r"))

off_spring_size = 1024
pop_size_accumulated = len(checkpoint['history'].genealogy_history.values())
gen_size = (pop_size_accumulated - off_spring_size)/(2*off_spring_size) +1

min_fitness_inds = list()
min_fitness_pop = list()

for gen in range(gen_size):
    if gen == 0:
        batch_size = off_spring_size
    else:
        batch_size = 2*off_spring_size
    
    if gen == 1:
        prev_batch_size = off_spring_size
    else:
        prev_batch_size = 2*off_spring_size
        
    batch_pop = checkpoint['history'].genealogy_history.values()[gen*prev_batch_size:\
                          (gen+1)*batch_size]
    batch_fitness = [batch_pop[i].fitness.sum for i in range(len(batch_pop))]
    min_fitness_batch = [x for _,x in sorted(zip(batch_fitness,batch_pop))]
    min_fitness = min(batch_fitness)
    
    min_fitness_inds.append(min_fitness_batch[0])
    min_fitness_pop.append(min_fitness)

checkpoint = None 


nrn = ephys.simulators.NrnSimulator()
fitness_protocols = opt.evaluator.fitness_protocols

frame_num = 60
frame_num1 = int(60*.67) # Select frame1 no. of frames from the initial generation
frame_num2 = frame_num - frame_num1 # Select frame2 no. of frames from the latter generation

mid_index = int(len(min_fitness_inds)/2)
init_generation_inds = min_fitness_inds[:mid_index]
min_fitness_pop1 =  min_fitness_pop[:mid_index]
ind1_space = int(math.ceil(len(init_generation_inds)/frame_num1))
selected_init_generation_inds = init_generation_inds[::ind1_space][:frame_num1]
gen1 = range(len(min_fitness_pop1))[::ind1_space][:frame_num1]
min_fitness1 = min_fitness_pop1[::ind1_space][:frame_num1]


latter_generation_inds = min_fitness_inds[mid_index:]
min_fitness_pop2 =  min_fitness_pop[mid_index:]
ind2_space = int(math.ceil(len(latter_generation_inds)/frame_num2))
selected_latter_generation_inds = latter_generation_inds[::ind2_space][:frame_num2]
gen2 = range(len(min_fitness_pop1),len(min_fitness_pop1)+len(min_fitness_pop2))[::ind2_space][:frame_num2]
min_fitness2 = min_fitness_pop2[::ind2_space][:frame_num2]

selected_individuals = selected_init_generation_inds + selected_latter_generation_inds
selected_min_fitness = min_fitness1 + min_fitness2
gen_numbers = gen1 + gen2
selected_fitness_gen = {'individual_fitness':selected_min_fitness,
                        'selected_gens' : gen_numbers}
 
selected_individuals_responses = []
   
key ='LongDC_46'

for i in range(len(selected_individuals)):
    ga_individual = dict()
    for index, param_name in enumerate(opt.evaluator.param_names):
        ga_individual[param_name] = selected_individuals[i][index]
    
    response_ind = fitness_protocols[key].run(
        cell_model=opt.evaluator.cell_model,
        param_values=ga_individual,
        sim=nrn)
    selected_individuals_responses.append(response_ind)

with open('opt_resp_ind.pkl', 'wb') as handle:
    pickle.dump(selected_individuals_responses, handle)
    
with open('selected_fitness.pkl', 'wb') as handle:
    pickle.dump(selected_fitness_gen, handle)    
    
