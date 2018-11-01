#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May  9 13:58:39 2018

@author: anirbannandi
"""

import pickle
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import glob
import numpy as np
import matplotlib.cm as cm
import json
import bluepyopt as bpopt
import pandas as pd
import os
import collections

import evaluator_helper


def Main():
    exp_dir = './preprocessed/*.txt'
    file_list = glob.glob(exp_dir)
    
    plt.style.use('fivethirtyeight')
    
    
    path_to_cell_metadata = os.path.abspath(os.path.join('.', os.pardir)) + '/cell_metadata.json'        
    with open(path_to_cell_metadata,'r') as metadata:
            cell_metadata = json.load(metadata) 
            
    cell_id = cell_metadata['Cell_id']
    layer = cell_metadata['Layer']
    area = cell_metadata['Area']
    species = cell_metadata['Species']
    cre_line = cell_metadata['Cre_line']
    dendrite_type = cell_metadata['Dendrite_type']
    feature_avg_released = cell_metadata['Feature_avg']
    explained_variance_released = cell_metadata['Explained_variance']
    
    validation_response_path = 'Validation_Responses/validation_responses.pkl'
    spiketimes_model_path = 'Validation_Responses/spiketimes_model.pkl'
    spiketimes_exp_path = 'Validation_Responses/spiketimes_exp.pkl'
    
    
    if os.path.exists(validation_response_path):
    
        model_response_pkl = pickle.load(open(validation_response_path,'r')) 
        
        noise_stimtypes = [
            'Noise_1',
            'Noise_2'
            ]        
        
        for key,model_resp in model_response_pkl.items():
            stim_type = key.rsplit('_',1)[0]
            if stim_type in noise_stimtypes:
                repetition_list = []
                for filename in file_list:
                    if stim_type in filename:
                        repetition_list.append(filename)
                len_plots = len(repetition_list)
                fig,ax = plt.subplots(len_plots+2,1,figsize = (8,10), dpi = 80, sharex = True)
                ax_index = 0
                color_index = np.linspace(0,.8,len_plots +1)
                for filename in repetition_list:
                    sweep = 'Sweep ' + filename.split('_')[2].replace('.txt', '')
                    exp_data = np.loadtxt(filename)
                    exp_data_time = exp_data[:,0]
                    exp_data_stim = exp_data[:,1]
                    exp_data_voltage = exp_data[:,2]
                    ax[ax_index].plot(exp_data_time, exp_data_voltage, 
                      color = cm.hot(color_index[ax_index]), lw =1, label = sweep)
                    ax[ax_index].set_ylabel('Voltage (mV)',fontsize = 14)
                    ax[-1].plot(exp_data_time, exp_data_stim, 
                      color = 'tomato', lw =1, label = 'Noise Stim')
                    
                    ax_index += 1
                ax[-1].set_ylabel('Stimulation (nA)',fontsize = 14)    
                ax[-1].set_xlabel('Time (ms)',fontsize = 14)
                model_response = model_resp['validation.soma.v']
                
                model_response_time = model_response['time'].as_matrix()
                model_response_voltage = model_response['voltage'].as_matrix()
                
            
                ax[ax_index].plot(model_response_time, model_response_voltage, 
                                  color = 'blue', lw = 1, label = 'Model')
                ax[ax_index].set_ylabel('Voltage (mV)',fontsize = 14)
                for ax_i in ax:
                    ax_i.legend(fontsize = 12)
                fig.suptitle('Validation for %s'%stim_type)
                fig.tight_layout(rect=[0, 0.03, 1, 0.95])
                fig.savefig('Validation_Responses/'+stim_type + '.png')
            
                plt.close(fig)
    
    ############### Explained variance calculation ###############
    

    if os.path.exists(spiketimes_exp_path):    
        spike_times_model_pkl = pickle.load(open(spiketimes_model_path,'r'))
        spike_times_exp_pkl = pickle.load(open(spiketimes_exp_path,'r'))    
        import exp_var_metric
        import math
        
        dt = 1/200.0 # ms
        sigma = [10] # ms
        
        exp_variance = {}
        for stim_type in noise_stimtypes:
            expt_trains = []
            for key,val in spike_times_exp_pkl.items():
                if stim_type in key:
                    expt_trains.append(val)
            for i,exp_train in enumerate(expt_trains):
                exp_train = [int(math.ceil(sp_time/dt)) for sp_time in exp_train]
                expt_trains[i] = np.asarray(exp_train)
            
            model_train = None
            for key,val in spike_times_model_pkl.items():
                if stim_type in key:
                    model_key = key
                    model_train = val
                    
            if model_train is not None:
                for i,sp_time in enumerate(model_train):
                    model_train[i] = int(math.ceil(sp_time/dt))
        
                model_train = model_train.astype(int)
                sweep_filename = 'preprocessed/'+model_key+'.txt'
                exp_data = np.loadtxt(sweep_filename)
                exp_data_time = exp_data[:,0]
                total_length = int(math.ceil(exp_data_time[-1]/dt))
                exp_variance[stim_type] = exp_var_metric.calculate_spike_time_metrics(expt_trains,
                            model_train, total_length, dt, sigma)
    else:
          exp_variance = {}      
    
    
    
    ################## Feature average calculation #####################
        
    
    
    with open('config_file.json') as json_file:  
        path_data = json.load(json_file)
    
    all_protocol_path = path_data['all_protocols']
    train_feature_path = path_data['features']
    untrained_feature_path = path_data['untrained_spiking_features']
    
    morph_path = path_data['morphology']
    mech_path = path_data['mechanism']
    param_path = path_data['parameters']
    
    evaluator_train = evaluator_helper.create(all_protocol_path, train_feature_path, morph_path, 
                                            param_path, mech_path)
    evaluator_validation = evaluator_helper.create(all_protocol_path, untrained_feature_path, morph_path, 
                                            param_path, mech_path)
    
    opt_train = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator_train)
    opt_validation = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator_validation)
    
    responses_filename = 'resp_opt.txt'
    with open(responses_filename, 'r') as fd:
        opt_responses = pickle.load(fd)
    opt_responses = opt_responses[0]
    
    objectives_train = opt_train.evaluator.fitness_calculator.calculate_scores(opt_responses)
    objectives_validation = opt_validation.evaluator.fitness_calculator.calculate_scores(opt_responses)
    
    objectives_train = collections.OrderedDict(sorted(objectives_train.iteritems()))
    objectives_validation = collections.OrderedDict(sorted(objectives_validation.iteritems()))
    
    feature_avg_train = np.mean(objectives_train.values())
    feature_avg_validation = np.mean(objectives_validation.values())    
    
    if bool(exp_variance):    
        if 'Noise_2' in exp_variance.keys():
            exp_variance_avg = exp_variance['Noise_2'][0]
        else:
            exp_variance_vals = exp_variance.values()
            exp_variance_avg = reduce(lambda x, y: x + y, exp_variance_vals) / len(exp_variance_vals)
            exp_variance_avg = exp_variance_avg[0]
    else:
        exp_variance_avg = None
    
    
    fitness_metrics = pd.DataFrame({'species' : [species],
                        'cell_id' : [cell_id],
                        'layer' : [layer],
                        'area' : [area],
                        'cre_line' : [cre_line],
                        'dendrite_type' : [dendrite_type],
                        'feature_avg_train' : [feature_avg_train],
                        'feature_avg_validation' : [feature_avg_validation],
                        'feature_avg_released' : [feature_avg_released],
                        'explained_variance' : [exp_variance_avg],
                        'explained_variance_released' : [explained_variance_released]
                        })
    
    fitness_filename = 'Validation_Responses/fitness_metrics_'+cell_id+'.csv' 
    fitness_metrics.to_csv(fitness_filename)    
    
if __name__ == '__main__':
    Main()
