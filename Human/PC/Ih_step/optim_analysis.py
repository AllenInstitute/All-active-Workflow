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
import pandas as pd
import matplotlib.pyplot as plt
import errno
import logging
from matplotlib.backends.backend_pdf import PdfPages
import math


logger = logging.getLogger(__name__)


with open('config_file.json') as json_file:  
    data = json.load(json_file)


release_params = data['release_params']        
fit_json_path = data['fit_json']
fit_protocol_path = data['protocols']
param_path = data['parameters']

section_map = {'somatic':'soma', 'axonal':'axon', 'apical':'apic',
               'basal':'dend', 'all':'all'}


with open(fit_protocol_path) as json_file:  
    train_protocols = json.load(json_file)

path_to_cell_metadata = os.path.abspath(os.path.join('.', os.pardir)) + '/cell_metadata.json'        
with open(path_to_cell_metadata,'r') as metadata:
        cell_metadata = json.load(metadata)
    
cell_id = cell_metadata['Cell_id']
layer = cell_metadata['Layer']
area = cell_metadata['Area']
species = cell_metadata['Species']

analysis_write_path = cell_id + '_analysis.pdf'
pdf_pages =  PdfPages(analysis_write_path)



def plot_diversity(opt, checkpoint_file, param_names):

    plt.style.use('ggplot')
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    import copy
    
    release_individual = copy.deepcopy(checkpoint['halloffame'][0])
    optimized_individual = [checkpoint['halloffame'][0]]


    param_values = opt.evaluator.params
    param_split_names = [name.split('.')[0] for name in param_names]
    unique = np.unique(np.asarray(param_split_names))
    ix = []
    for u in unique:
        for i,param in enumerate(param_split_names):
            if param == u:
                ix.append(i)
    param_names_arranged = [param_names[k] for k in ix]
    sect_names_arranged = [sect.split('.')[1] for sect in param_names_arranged]

    param_values_arranged = [param_values[k] for k in ix]
    optimized_individual_arranged = [optimized_individual[0][k] for k in ix]

    
    hof_list = checkpoint['halloffame'][1:]
    hof_list_arranged = []
    hof_df = pd.DataFrame([])
    for i in range(len(hof_list)):
        arranged_hof_item = [hof_list[i][j] for j in ix]
        hof_list_arranged.append(arranged_hof_item)
        temp_df = pd.DataFrame({'param_name' : param_names_arranged,
        'section':sect_names_arranged,
        'value' : arranged_hof_item,
        'label' : 'hall of fame',
        'fitness' : 2,
        'cell_id' : cell_id,
        'layer' : layer,
        'area' : area,
        'species' : species})
        hof_df=hof_df.append(temp_df) 
    

    if release_params:    
        for index, param_name in enumerate(param_names_arranged):
            release_individual[index] = release_params[param_name]

    fig, ax = plt.subplots(1,figsize=(6,6))
    
    param_count = len(param_names_arranged)
    x = np.arange(param_count)
    
    def add_line(ax, xpos, ypos):
        line = plt.Line2D([xpos, xpos], [ypos + .1, ypos],
                          transform=ax.transAxes, color='black')
        line.set_clip_on(False)
        ax.add_line(line)
        
        
    
            
    for k in range(len(hof_list_arranged)):
        if k == 0:
            ax.scatter(x, map(abs,hof_list_arranged[k]), marker = 'o', 
               alpha = 0.2, s=20, color= 'green', edgecolor='black',
               label='hall of fame')
        else:
            ax.scatter(x, map(abs,hof_list_arranged[k]), marker = 'o', 
               alpha = 0.2, s=20, color= 'green')
    
    abs_optimized_individual_arranged  = map(abs,optimized_individual_arranged)      
    ax.scatter(x, abs_optimized_individual_arranged, marker = 'x', 
               alpha = 1, s=200, color= 'blue', edgecolor='black',
               label='optimized')
    if release_params:
        abs_release_individual = map(abs,release_individual)
        ax.scatter(x, abs_release_individual, marker = 'x', 
                   alpha = 0.8, s=100, color= 'red', edgecolor='black'
                   , label = 'released')
    
    tick_labels = param_names_arranged
    
    def plot_tick(column, y):
                col_width = 0.25
                x = [column - col_width,
                     column + col_width]
                y = [y, y]
                ax.plot(x, y, color='black')
                
    min_list = list()            
    for i, parameter in enumerate(param_values_arranged):
        min_value = abs(parameter.lower_bound)
        max_value = abs(parameter.upper_bound)
        min_list.append(min_value)
        plot_tick(i, min_value)
        plot_tick(i, max_value)

    plt.xticks(x, tick_labels, rotation=90, ha='center')
    for xline in x:
        ax.axvline(xline, linewidth=1, color='white',linestyle=':')  
    ax.set_yscale('log')
    if min(min_list)<1e-10:
        ax.set_ylim((1e-10, 1e4))
    ax.set_xlim((-1, x[-1]+1))
    ax.set_ylabel('Parameter Values (Absolute)')
    ax.legend(prop={'size': 8}, frameon= True, shadow=True, fancybox=True)
    ax.set_title('Parameters')
    plt.tight_layout()
    pdf_pages.savefig(fig)
    plt.close(fig)


    # save optimized parameters in fit.json format
    
    fit_json_write_path = 'fitted_params/optim_param_'+cell_id+ '.json'
    if not os.path.exists(os.path.dirname(fit_json_write_path)):
        try:
            os.makedirs(os.path.dirname(fit_json_write_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    
    with open(param_path) as json_file:  
        params = json.load(json_file)
    
    optimized_param_dict = {key:optimized_individual_arranged[i] for i,key in \
                            enumerate(param_names_arranged)} 
    
    param_dict_final = {key.split('.')[0]+'.'+
                     section_map[key.split('.')[1]] : optimized_param_dict[key] 
                                        for key in optimized_param_dict.keys()} 
    with open(fit_json_path) as json_file:  
        model_data = json.load(json_file)
        
    param_dict_final_keys = param_dict_final.keys()
    for key in param_dict_final.keys():
        opt_name,opt_sect = key.split('.')
        data_key = 'genome'
        for j in range(len(model_data[data_key])):
            if model_data[data_key][j]['name'] == opt_name and \
                model_data[data_key][j]['section'] == opt_sect:
                model_data[data_key][j]['value'] = str(param_dict_final[key])
                param_dict_final_keys.remove(key)

    
    for key in param_dict_final_keys:
        param_name,sect = key.split('.') 
        model_data['genome'].append(
                {
                  'section' : sect,
                  'name'    : param_name,
                  'value'   : str(param_dict_final[key]),
                  'mechanism': 'Ih'      
                })

    model_data['passive'] = [{'ra' : param_dict_final['Ra.all']}]
    model_data['conditions'][0]['v_init'] = (item['value'] for item in params if \
                                item["param_name"] == "v_init").next()
    
    with open(fit_json_write_path, 'w') as outfile:
        json.dump(model_data, outfile,indent=4)
    
    print 'Saving the parameters in fit.json format'
   
    optim_param_write_path = 'optim_param_unformatted.json'
    with open(optim_param_write_path, 'w') as outfile:
        json.dump(optimized_param_dict, outfile,indent=4)
        
        
   # save parameters in a csv file to later plot in R
    if release_params:        
        released_df = pd.DataFrame({ 'param_name' : param_names_arranged,
            'section':sect_names_arranged,                        
            'value' : release_individual,
            'label' : 'Released',
            'fitness': 5,
            'cell_id' : cell_id,
            'layer' : layer,
            'area' : area,
            'species' : species})
    else:
          released_df = pd.DataFrame([])  
    
    optimized_df = pd.DataFrame({'param_name' : param_names_arranged,
        'section':sect_names_arranged,                         
        'value' : optimized_individual_arranged,
        'label' : 'Optimized',
        'fitness': 10,
        'cell_id' : cell_id,
        'layer' : layer,
        'area' : area,
        'species' : species})
    param_df = [optimized_df, released_df,hof_df] 
    param_df = pd.concat(param_df)   
    param_df.to_csv('params_'+cell_id+'.csv')
    
    print 'Saving the parameters in .csv for plotting in R'    
    
#####################################################################

# GA evolution

def plot_GA_evolution(checkpoint_file):
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    from matplotlib.ticker import MaxNLocator
    plt.style.use('ggplot')
    log = checkpoint['logbook']
    gen_numbers = log.select('gen')
    mean = np.array(log.select('avg'))
    std = np.array(log.select('std'))
    minimum = np.array(log.select('min'))
    stdminus = mean - std
    stdplus = mean + std
   
    fig, ax = plt.subplots(1, figsize=(8,8))

    ax.plot(gen_numbers,
            mean,
            color="white",
            linewidth=2,
            label='population average')
    
    
    
    ax.fill_between(gen_numbers,
        stdminus,
        stdplus,
        color='#3F5D7D',
        alpha = 0.5,
        label='population standard deviation')
    
    ax.plot(gen_numbers,
            minimum,
            color='red',
            linewidth=2,
            alpha = 0.8,
            label='population minimum')
    
    ax.legend(prop={'size': 12}, frameon= True, 
                    shadow=True, fancybox=True)
    
    left, bottom, width, height = [0.67, 0.6, 0.2, 0.15]
    ax2 = fig.add_axes([left, bottom, width, height])

    ax2.plot(gen_numbers,
            minimum,
            linewidth=2,
            color='red',
            alpha = 0.8)

    ax2.set_facecolor('#EAECEE')
    ax2.xaxis.set_major_locator(MaxNLocator(5))
    ax2.yaxis.set_major_locator(MaxNLocator(4))
    ax.set_xlim((1,gen_numbers[-1]))
    ax.set_xlabel('Generation #')
    ax.set_ylabel('Sum of objectives')
    ax.set_title('Evolution of the Objective')
    pdf_pages.savefig(fig)
    plt.close(fig)
    pdf_pages.close()
    
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


def feature_comp(opt, checkpoint_file,responses_filename):
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    hof = checkpoint['halloffame']

    
    responses = get_responses(opt.evaluator, hof, responses_filename)
    responses = responses[0] #select the optimized response     
    
    # objectives
    objectives = opt.evaluator.fitness_calculator.calculate_scores(responses)
    
    import collections
    objectives = collections.OrderedDict(sorted(objectives.iteritems()))
    
    feature_split_names = [name.split('.')[-1] for name in objectives.keys()]
    features = np.unique(np.asarray(feature_split_names))
    

    
    bar_width = 0.35
    opacity = 0.4

    plt.style.use('ggplot') 
    for i, feature in enumerate(features):
        fig, ax = plt.subplots(1, figsize=(8,8))    
        iter_dict = dict()
        iter_dict_release = dict()
        tick_label = []
        for key in objectives.keys():
            if key.split('.')[-1] == feature:
                iter_dict[key] = objectives[key]
                amp = train_protocols[key.split('.')[0]]['stimuli'][0]['amp']
                amp_reduced = round(amp,2)
                tick_label.append(amp_reduced)    
        index = np.arange(len(iter_dict.keys()))
       
        xtick_pos = index
        iter_dict = collections.OrderedDict(sorted(iter_dict.iteritems()))
        iter_dict_release = collections.OrderedDict(sorted(iter_dict_release.iteritems()))
        ax.bar(index,
              iter_dict.values(),
              bar_width,
              align='center',
              color='b',
              alpha=opacity,
              label='Optimized')
              
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(tick_label, fontsize= 8)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(8) 
        ax.set_ylabel('Objective value (# std)',fontsize= 12)
        ax.set_xlabel('Stimulus Amplitude (nA)',fontsize= 12)        
        ax.legend(prop={'size': 12},loc ='best')
        ax.set_title(feature, fontsize= 14)
        fig.tight_layout(rect=[0.05, 0.05, .95, 0.95])
        pdf_pages.savefig(fig)
        plt.close(fig)    
               
        logger.debug('Plotting comparisons for %s \n'%feature)
    
    stim_file = 'preprocessed/StimMapReps.csv'
    stim_df = pd.read_csv(stim_file, sep='\s*,\s*',
                           header=0, encoding='ascii', engine='python')
    
    csv_filename = 'error' + '_' + cell_id + '.csv'        

    feature_df = pd.DataFrame([])
    for key,val in objectives.iteritems():
        proto = key.split('.')[0]
        temp_df = pd.DataFrame({'error_fit' : val,
            'stim_amp': stim_df[stim_df['DistinctID'] == proto]['Amplitude_Start'],
            'cell_id' : cell_id,
            'layer' : layer,
            'area' : area,
            'species' : species,
            'feature' : key.split('.')[-1],
            'label' : 'Optimized'
            })

        
        feature_df=feature_df.append(temp_df) 
    feature_df.to_csv(csv_filename)

#############################################################################

## Plotting responses


def plot_Response(opt,checkpoint_file, responses_filename):
    stim_file = 'preprocessed/StimMapReps.csv'
    stim_df = pd.read_csv(stim_file, sep='\s*,\s*',
                           header=0, encoding='ascii', engine='python')
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    hof = checkpoint['halloffame']

    responses = get_responses(opt.evaluator, [hof[0]], responses_filename)
    response = responses[0]

    
    plt.style.use('ggplot') 
    training_plots = 0
    validation_plots = 0
    protocol_names = stim_df['DataPath'].sort_values()
    for i, trace_rep in enumerate(protocol_names):
        rep_id = trace_rep.split('|')[0]
        if rep_id.split('.')[0] in opt.evaluator.fitness_protocols.keys():
            if rep_id.split('.')[0] in train_protocols.keys():
                training_plots += len(trace_rep.split('|'))
            else:
                validation_plots += len(trace_rep.split('|'))
    n_col = 3   
    n_row_train =  int(math.ceil(training_plots/float(n_col)))
    n_row_val =  int(math.ceil(validation_plots/float(n_col)))
    fig_empty_index_train = range(training_plots,n_row_train*n_col)
    fig_empty_index_val = range(validation_plots,n_row_val*n_col)
    fig_train,ax_train = plt.subplots(n_row_train,n_col, figsize=(10,10))
    fig_val,ax_val = plt.subplots(n_row_val,n_col, figsize=(10,10))

    for ax_c,empty_index in zip([ax_val,ax_train],[fig_empty_index_val,fig_empty_index_train]):
        if len(empty_index) != 0:
            for ind in empty_index:
                ax_c[ind/n_col,ind%n_col].axis('off')
        
    
    index_train = 0
    index_val = 0
    for i, trace_rep in enumerate(protocol_names):
        rep_id = trace_rep.split('|')[0]
        name_loc = rep_id.split('.')[0] +'.soma.v'
        if rep_id.split('.')[0] in opt.evaluator.fitness_protocols.keys():
            if rep_id.split('.')[0] in train_protocols.keys():
                analysis_state = '(T)'
                ax = ax_train
                index = index_train
            else:
                analysis_state = '(V)'
                ax = ax_val
                index = index_val
            for name in trace_rep.split('|'):
                
                response_time = response[name_loc]['time']
                response_voltage = response[name_loc]['voltage']
                color = 'blue'
                l1, = ax[index/n_col,index%n_col].plot(response_time,
                        response_voltage,
                        color=color,
                        linewidth=1,
                        label= 'Model',
                        alpha = 0.8)

                FileName = 'preprocessed/' + name
                data = np.loadtxt(FileName) 
                l3, = ax[index/n_col,index%n_col].plot(data[:,0],
                            data[:,1],
                            color='black',
                            linewidth=1,
                            label = 'Experiment',
                            alpha = 0.5)  
                
                
                if (analysis_state == '(T)' and index/n_col == n_row_train-1) or \
                                (analysis_state == '(V)' and index/n_col == n_row_val-1): 
                    ax[index/n_col,index%n_col].set_xlabel('Time (ms)')
                
                if (analysis_state == '(T)' and index%n_col == 0) or \
                                (analysis_state == '(V)' and index%n_col == 0): 
                    ax[index/n_col,index%n_col].set_ylabel('Voltage (mV)')
                
                ax[index/n_col,index%n_col].set_title(name.split('.')[0], fontsize=8)
                
                if 'LongDCSupra' in name:
                    ax[index/n_col,index%n_col].set_xlim([0, 3000])
                elif 'LongDC' in name:
                    ax[index/n_col,index%n_col].set_xlim([0, 1600])
                    
                logger.debug('Plotting response comparisons for %s \n'%name.split('.')[0])
                index += 1
            if analysis_state == '(T)':
                index_train = index
                if index_train == training_plots:
                    fig_train.legend(handles = (l1,l3), loc = 'lower center', ncol=2)
            else:
                index_val = index
                if index_val == validation_plots:
                    fig_val.legend(handles = (l1,l3), loc = 'lower center', ncol=2)
            
    fig_train.suptitle('Training Set',fontsize=16) 
    fig_val.suptitle('Test Set',fontsize=16)       
    fig_val.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig_train.tight_layout(rect=[0, 0.03, 1, 0.95])
    pdf_pages.savefig(fig_train)
    plt.close(fig_train)
    pdf_pages.savefig(fig_val)
    plt.close(fig_val)
    
    