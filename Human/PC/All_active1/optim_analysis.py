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
import bluepyopt as bpopt
import copy
from collections import defaultdict
import efel
from scipy import interpolate


import evaluator_helper


logger = logging.getLogger(__name__)

with open('config_file.json') as json_file:  
    data = json.load(json_file)

section_map = {'somatic':'soma', 'axonal':'axon', 'apical':'apic',
               'basal':'dend', 'all':'all'}

release_params = data['release_params']
release_params_original = {}        
fit_json_path = data['fit_json']

morph_path = data['morphology']
protocol_path = data['protocols']
mech_path = data['mechanism']
feature_path = data['features']
param_path = data['parameters']
all_protocol_path = data['all_protocols']

with open(protocol_path) as json_file:  
    train_protocols = json.load(json_file)
    
with open(param_path) as json_file:  
    parameters_optim = json.load(json_file)


path_to_cell_metadata = os.path.abspath(os.path.join('.', os.pardir)) + '/cell_metadata.json'        
with open(path_to_cell_metadata,'r') as metadata:
        cell_metadata = json.load(metadata)
    
cell_id = cell_metadata['Cell_id']
layer = cell_metadata['Layer']
area = cell_metadata['Area']
species = cell_metadata['Species']
cre_line = cell_metadata['Cre_line']
dendrite_type = cell_metadata['Dendrite_type']

analysis_write_path = cell_id + '_analysis_Stage2.pdf'
pdf_pages =  PdfPages(analysis_write_path)

def plot_diversity(opt, checkpoint_file, param_names, hof_index):

    plt.style.use('ggplot')
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    
    release_individual = list()
    optimized_individual = [checkpoint['halloffame'][hof_index]]


    param_values = opt.evaluator.params
    param_names_arranged = sorted(param_names)
    sect_names_arranged = [sect.split('.')[1] for sect in param_names_arranged]
    ix = sorted(range(len(param_names)), key=lambda k: param_names[k])
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
        'area' : cell_metadata['Area'],
        'species' : cell_metadata['Species'],
        'cre_line' : cell_metadata['Cre_line'],
        'dendrite_type' : cell_metadata['Dendrite_type']})
        hof_df=hof_df.append(temp_df) 
    

    all_params = param_names_arranged + list(set(release_params.keys()) - 
                                             set(param_names_arranged)) 
    all_params_arranged = sorted(all_params)
    release_params_arranged = sorted(release_params.keys())
    sect_release_names_arranged = [sect.split('.')[1] for sect in release_params_arranged]
    
    for param_name in release_params_arranged:
        release_individual.append(release_params[param_name])

    fig, ax = plt.subplots(1,figsize=(6,6))
    
    param_count = len(all_params_arranged)
    x = np.arange(param_count)
    
    x_opt = []
    x_release = []
    for i in range(param_count):
        if all_params_arranged[i] in param_names_arranged:
            x_opt.append(i)
        if all_params_arranged[i] in release_params_arranged:
            x_release.append(i)
    
    def add_line(ax, xpos, ypos):
        line = plt.Line2D([xpos, xpos], [ypos + .1, ypos],
                          transform=ax.transAxes, color='black')
        line.set_clip_on(False)
        ax.add_line(line)
        
        
    
            
    for k in range(len(hof_list_arranged)):
        if k == 0:
            ax.scatter(x_opt, map(abs,hof_list_arranged[k]), marker = 'o', 
               alpha = 0.2, s=10, color= 'green', edgecolor='black',
               label='hall of fame')
        else:
            ax.scatter(x_opt, map(abs,hof_list_arranged[k]), marker = 'o', 
               alpha = 0.2, s=10, color= 'green')
    
    abs_optimized_individual_arranged  = map(abs,optimized_individual_arranged)      
    ax.scatter(x_opt, abs_optimized_individual_arranged, marker = 'x', 
               alpha = 1, s=100, color= 'blue', edgecolor='black',
               label='optimized')
    abs_release_individual = map(abs,release_individual)
    ax.scatter(x_release, abs_release_individual, marker = 'x', 
               alpha = 0.8, s=50, color= 'red', edgecolor='black',
               label = 'released')
    
    tick_labels = all_params_arranged
    
    def plot_tick(column, y):
                col_width = 0.3
                x = [column - col_width,
                     column + col_width]
                y = [y, y]
                ax.plot(x, y, color='black',linewidth = 1)
                
    min_list = list()            
    for i, parameter in zip(x_opt,param_values_arranged):
        min_value = abs(parameter.lower_bound)
        max_value = abs(parameter.upper_bound)
        min_list.append(min_value)
        plot_tick(i, min_value)
        plot_tick(i, max_value)

    plt.xticks(x, tick_labels, rotation=90, ha='center', fontsize = 6)
    plt.yticks(fontsize = 8)
    for xline in x:
        ax.axvline(xline, linewidth=.5, color='white',linestyle=':')  
    ax.set_yscale('log')
    if min(min_list)<1e-10:
        ax.set_ylim((1e-10, 1e4))
    ax.set_xlim((-1, x[-1]+1))
    ax.set_ylabel('Parameter Values (Absolute)', fontsize = 8)
    ax.legend(prop={'size': 8}, frameon= True, shadow=True, fancybox=True)
    ax.set_title('Parameters')
    plt.tight_layout()
    pdf_pages.savefig(fig)
    plt.close(fig)



    # save optimized parameters in fit.json format
    
    fit_json_write_path = 'fitted_params/optim_param_'+str(cell_id)+ '.json'
    if not os.path.exists(os.path.dirname(fit_json_write_path)):
        try:
            os.makedirs(os.path.dirname(fit_json_write_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    
    optimized_param_dict = {key:optimized_individual_arranged[i] for i,key in \
                            enumerate(param_names_arranged)} 
    
    param_dict_final = {key.split('.')[0]+'.'+
                     section_map[key.split('.')[1]] : optimized_param_dict[key] 
                                        for key in optimized_param_dict.keys()} 
    with open(fit_json_path) as json_file:  
        model_data = json.load(json_file)
    
    repeat_list = list()
    data_key = 'genome'
    model_data[data_key] = []
    
    for key in param_dict_final.keys():
        opt_name,opt_sect = key.split('.')
        mech_list = opt_name.split('_',1)
        
        if len(mech_list) > 1 and 'pas' not in mech_list:
            model_data[data_key].append({'name' : opt_name,
                                         'mechanism' : mech_list[1],
                                         'section' : opt_sect,
                                         'value' : str(param_dict_final[key])})
              
        else:
            repeat_list.append(key)
            
    for passive_param_key in repeat_list:
        opt_name,_ = passive_param_key.split('.')
        for sect in ['soma','axon','apic','dend']:
            model_data[data_key].append({'section' : sect,
                                         'name' : opt_name,
                                         'value' : str(param_dict_final[passive_param_key]),
                                         'mechanism' : ''})
    
    model_data['passive'] = [{'ra' : param_dict_final['Ra.all']}] 
    model_data['conditions'][0]['v_init'] = [parameter_optim['value'] for parameter_optim in \
               parameters_optim if parameter_optim['param_name'] == 'v_init'][0]
    
    with open(fit_json_write_path, 'w') as outfile:
        json.dump(model_data, outfile,indent=4)
    
    optim_param_write_path = 'optim_param_unformatted.json'    
    with open(optim_param_write_path, 'w') as outfile:
        json.dump(optimized_param_dict, outfile,indent=4)
        
    logger.debug('Saving the parameters in fit.json format')
   
    
   # save parameters in a csv file to later plot in R
            
    released_df = pd.DataFrame({ 'param_name' : release_params_arranged,
        'section':sect_release_names_arranged,                        
        'value' : release_individual,
        'label' : 'Released',
        'fitness': 5,
        'cell_id' : cell_id,
        'layer' : layer,
        'area' : cell_metadata['Area'],
        'species' : cell_metadata['Species'],
        'cre_line' : cell_metadata['Cre_line'],
        'dendrite_type' : cell_metadata['Dendrite_type']})
    
    optimized_df = pd.DataFrame({'param_name' : param_names_arranged,
        'section':sect_names_arranged,                         
        'value' : optimized_individual_arranged,
        'label' : 'Optimized',
        'fitness': 10,
        'cell_id' : cell_id,
        'layer' : layer,
        'area' : cell_metadata['Area'],
        'species' : cell_metadata['Species'],
        'cre_line' : cell_metadata['Cre_line'],
        'dendrite_type' : cell_metadata['Dendrite_type']})
        
    param_df = [optimized_df, released_df,hof_df] 
    param_df = pd.concat(param_df)  
    csv_filename = 'params_'+str(cell_id)+ '.csv'
    param_df.to_csv(csv_filename)
    
    logger.debug('Saving the parameters in .csv format')    
    
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
#    pdf_pages.close()
    
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

    # objectives
    responses = get_responses(opt.evaluator, hof, responses_filename)
    responses = responses[0] #select the optimized response 
    
    objectives = opt.evaluator.fitness_calculator.calculate_scores(responses)
    
    import collections
    objectives = collections.OrderedDict(sorted(objectives.iteritems()))
    
    feature_split_names = [name.split('.',1)[-1] for name in objectives.keys()]
    features = np.unique(np.asarray(feature_split_names))
    
    bar_width = 0.35
    opacity = 0.4    
    
    plt.style.use('ggplot') 
    for i, feature in enumerate(features):
        fig, ax = plt.subplots(1, figsize=(8,8))    
        iter_dict = collections.defaultdict(list)
        tick_label = []
        for key in objectives.keys():
            if key.split('.',1)[-1] == feature:
                
                amp = train_protocols[key.split('.')[0]]['stimuli'][0]['amp']
                amp_reduced = round(amp,3)
                iter_dict[str(amp_reduced)].append(objectives[key])
        index = np.arange(len(iter_dict.keys()))
        xtick_pos = index 
        iter_dict ={key:np.mean(val) for key,val in iter_dict.items()}
        iter_dict = collections.OrderedDict(sorted(iter_dict.iteritems()))
        tick_label = iter_dict.keys()
        ax.bar(index,
              iter_dict.values(),
              bar_width,
              align='center',
              color='b',
              alpha=opacity,
              label='Optimized')
  
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(tick_label, fontsize= 8)
        plt.xticks(rotation=90)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(8) 
        ax.set_ylabel('Objective value (# std)',fontsize= 12)
        ax.set_xlabel('Stimulus Amplitude (nA)',fontsize= 12)        
        ax.legend(prop={'size': 10},loc ='best')
        if 'soma' in feature:
            feature_title = feature.split('.')[1]
        else:
            feature_title = feature.split('.')[1] + ' at ' + feature.split('.')[0]
        ax.set_title(feature_title, fontsize= 12)
        fig.tight_layout(rect=[0.05, 0.05, .95, 0.95])
        pdf_pages.savefig(fig)
        plt.close(fig)    
        logger.debug('Plotting comparisons for %s \n'%feature)
    
    stim_file = 'preprocessed/StimMapReps.csv'
    stim_df = pd.read_csv(stim_file, sep='\s*,\s*',
                           header=0, encoding='ascii', engine='python')
    csv_filename = 'error' + '_' + str(cell_id) + '.csv'        

    feature_df = pd.DataFrame([])
    for key,val in objectives.iteritems():
        proto = key.split('.')[0]
        feature_name = key.split('.',1)[-1]
        
        # Remove dend1,2... specifications
        
        for i in range(10):
            if str(i) in feature_name:
                feature_name = feature_name.replace(str(i),'')
                
        temp_df = pd.DataFrame({'error_fit' : val,
            'stim_amp': stim_df[stim_df['DistinctID'] == proto]['Amplitude_Start'],
            'cell_id' : cell_id,
            'layer' : layer,
            'feature' : key.split('.')[-1],
            'label' : 'Optimized',
            'area' : cell_metadata['Area'],
            'species' : cell_metadata['Species'],
            'cre_line' : cell_metadata['Cre_line'],
            'dendrite_type' : cell_metadata['Dendrite_type']
            })
       
        
        feature_df=feature_df.append(temp_df) 
    feature_df.to_csv(csv_filename)

#############################################################################

## Plotting responses


def plot_Response(opt,checkpoint_file, responses_filename,hof_index):
    stim_file = 'preprocessed/StimMapReps.csv'
    stim_df = pd.read_csv(stim_file, sep='\s*,\s*',
                           header=0, encoding='ascii', engine='python')
    
    with open(checkpoint_file,'r') as handle:
        checkpoint = pickle.load(handle)
    hof = checkpoint['halloffame']

    responses = get_responses(opt.evaluator, [hof[hof_index]], responses_filename)
    
    
    plt.style.use('ggplot') 
    training_plots = 0
    
    protocol_names_original = stim_df['DataPath'].tolist()
    amp_start_original = stim_df['Stim_Start'].tolist()
    amp_end_original = stim_df['Stim_End'].tolist()

    protocol_names = sorted(protocol_names_original)
    idx = np.argsort(protocol_names_original)
    amp_start_list = [amp_start_original[i] for i in idx]
    amp_end_list = [amp_end_original[i] for i in idx]
    
    for i, trace_rep in enumerate(protocol_names):
        rep_id = trace_rep.split('|')[0]
        if rep_id.split('.')[0] in opt.evaluator.fitness_protocols.keys():
            training_plots += len(trace_rep.split('|'))
    n_col = 3  
    n_row = 5
    fig_per_page = n_col * n_row
    fig_pages =  int(math.ceil(training_plots/float(fig_per_page)))
    fig_mat = list()
    ax_mat = list()
    for page_i in range(fig_pages):
        
        fig,ax = plt.subplots(n_row,n_col, figsize=(10,10),squeeze = False)
        
        if page_i == fig_pages-1:
            remain_plots = training_plots - n_row*n_col*(fig_pages-1)
            fig_empty_index_train = range(remain_plots,n_row*n_col)
    

            if len(fig_empty_index_train) != 0:
                for ind in fig_empty_index_train:
                    ax[ind/n_col,ind%n_col].axis('off')
                    
        fig_mat.append(fig)
        ax_mat.append(ax)   
        
        
    index = 0
    index_plot = 0
    fig_index = 0
    for ix,trace_rep in enumerate(protocol_names):
        rep_id = trace_rep.split('|')[0]
        if rep_id.split('.')[0] in train_protocols.keys():
            state = ' (Train)'
        else:
            state = ' (Test)'
        name_loc = rep_id.split('.')[0] +'.soma.v'
        if rep_id.split('.')[0] in opt.evaluator.fitness_protocols.keys():
            
            for name in trace_rep.split('|'):
                ax_comp = ax_mat[fig_index]
                fig_comp = fig_mat[fig_index]
                response = responses[0]
                response_time = response[name_loc]['time']
                response_voltage = response[name_loc]['voltage']

                color = 'blue'
                l1, = ax_comp[index/n_col,index%n_col].plot(response_time,
                        response_voltage,
                        color=color,
                        linewidth=1,
                        label= 'Optimized',
                        alpha = 0.8)                    

                FileName = 'preprocessed/' + name
                data = np.loadtxt(FileName) 
                l3, = ax_comp[index/n_col,index%n_col].plot(data[:,0],
                            data[:,1],
                            color='black',
                            linewidth=1,
                            label = 'Cell Response',
                            alpha = 0.6)  

        
#                if index/n_col == n_row-1 and index%n_col == 0:
#                    fig_comp.legend(handles = (l1,l3),  loc = 'lower center', ncol=2)
                    
                if index/n_col == n_row-1: 
                    ax_comp[index/n_col,index%n_col].set_xlabel('Time (ms)')
                if index%n_col == 0: 
                    ax_comp[index/n_col,index%n_col].set_ylabel('Voltage (mV)')
                
                if 'LongDC' in name:
                    ax_comp[index/n_col,index%n_col].set_xlim([amp_start_list[ix]-200,\
                                                              amp_end_list[ix] +200])
                
                ax_comp[index/n_col,index%n_col].set_title(name.split('.')[0] + state, fontsize=8)
                logger.debug('Plotting response comparisons for %s \n'%name.split('.')[0])
                index += 1
                index_plot +=1
                if index%fig_per_page == 0 or index_plot == training_plots:
                    fig_comp.suptitle('Response Comparisons',fontsize=16)
                    fig_comp.legend(handles = (l1,l3),  loc = 'lower center', ncol=2)
                    fig_comp.tight_layout(rect=[0, 0.03, 1, 0.95])
                    pdf_pages.savefig(fig_comp)
                    plt.close(fig_comp)
                    index = 0
                    fig_index += 1
                    

def check_for_block(t_vec,v_vec,stim_params):

    stim_delay = stim_params['delay']
    stim_dur = stim_params['duration']
    v = v_vec.as_matrix()
    t = t_vec.as_matrix()
    stim_start_idx = np.flatnonzero(t >= stim_delay)[0]
    stim_end_idx = np.flatnonzero(t >= stim_delay + stim_dur)[0]
    depol_block_threshold = -50.0 # mV
    block_min_duration = 30.0 # ms
    long_hyperpol_threshold = -75.0 # mV

    bool_v = np.array(v > depol_block_threshold, dtype=int)
    up_indexes = np.flatnonzero(np.diff(bool_v) == 1)
    down_indexes = np.flatnonzero(np.diff(bool_v) == -1)
    if len(up_indexes) > len(down_indexes):
        down_indexes = np.append(down_indexes, [stim_end_idx])

    if len(up_indexes) == 0:
        # if it never gets high enough, that's not a good sign (meaning no spikes)
        return True
    else:
        max_depol_duration = np.max([t[down_indexes[k]] - t[up_idx] for k, up_idx in enumerate(up_indexes)])
        if max_depol_duration > block_min_duration:
            return True
    bool_v = np.array(v > long_hyperpol_threshold, dtype=int)
    up_indexes = np.flatnonzero(np.diff(bool_v) == 1)
    down_indexes = np.flatnonzero(np.diff(bool_v) == -1)
    down_indexes = down_indexes[(down_indexes > stim_start_idx) & (down_indexes < stim_end_idx)]
    if len(down_indexes) != 0:
        up_indexes = up_indexes[(up_indexes > stim_start_idx) & (up_indexes < stim_end_idx) & (up_indexes > down_indexes[0])]
        if len(up_indexes) < len(down_indexes):
            up_indexes = np.append(up_indexes, [stim_end_idx])
        max_hyperpol_duration = np.max([t[up_indexes[k]] - t[down_idx] for k, down_idx in enumerate(down_indexes)])
        if max_hyperpol_duration > block_min_duration:
            return True
    return False


def DB_check(checkpoint_file, DB_protocol_path, DB_response_path):
    
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    hof = checkpoint['halloffame'] 
    
    if os.path.exists(DB_response_path):
        responses = pickle.load(open(DB_response_path, 'r'))
       
        with open(DB_protocol_path) as json_file:  
            DB_protocols = json.load(json_file)
        select_protocols = [key for key in DB_protocols.keys() if 'DB' in key]
        
    else:
    
        # find the max amplitude used in training
        max_amp = 0
        for key, val in train_protocols.items():
            amp = val['stimuli'][0]['amp']
            if amp >= max_amp:
                max_key = key
                max_amp = amp
        

        # Increment of 10 and 20% of the maximum stimulation 
        # amplitude on which the model is trained        
        DB_protocols = copy.deepcopy(train_protocols)
        DB_protocol_increment = [.1*max_amp,.2*max_amp]
        select_protocols = []
        
        for i,inc in enumerate(DB_protocol_increment):
            DB_protocol_name = 'DB_'+str(i)
            temp_copy = copy.deepcopy(train_protocols)
            DB_protocols[DB_protocol_name] = temp_copy[max_key]
            DB_protocols[DB_protocol_name]['stimuli'][0]['amp'] += inc
            DB_protocols[DB_protocol_name]['stimuli'][0]['amp_end'] += inc
            select_protocols.append(DB_protocol_name)
        
        with open(DB_protocol_path, 'w') as json_file:
            json.dump(DB_protocols,json_file,indent=4)
            
        evaluator = evaluator_helper.create(DB_protocol_path, feature_path, morph_path, 
                                            param_path, mech_path)
            
        opt = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator)
    
        
        
        fitness_protocols = opt.evaluator.fitness_protocols
        select_run_protocols = []
        for protocol in fitness_protocols.values():
            if protocol.name in select_protocols:
                select_run_protocols.append(protocol)
        
        
        responses = defaultdict(list)        
        for i in range(len(hof)):
            individual = hof[i]
            individual_dict = evaluator.param_dict(individual)
            responses[i].append(evaluator.run_protocols(
                select_run_protocols,
                param_values=individual_dict)) 
    
    DB_dict = defaultdict(dict)
    
    # stim params same for both DB protocols
    stim_params = DB_protocols['DB_0']['stimuli'][0]
    locations = ['soma.v']
    
    hof_indices_DB_checked  = []
    
    for i in range(len(hof)):

        for proto in select_protocols:
            for loc in locations:
                name_loc = proto+'.'+loc
                t_vec = responses[i][0][name_loc]['time']
                v_vec = responses[i][0][name_loc]['voltage']
    
                depol_block = check_for_block(t_vec,v_vec,stim_params)
                DB_dict[i][name_loc] = depol_block
        
        depol_block_list = DB_dict[i].values()
        
        # if both the increased stimulus passes DB check add it to the hof list
        if not(depol_block_list[0]) and not(depol_block_list[0]):
            hof_indices_DB_checked.append(i)
            
    
    pickle.dump(responses, open(DB_response_path, 'w'))
    
    DB_check_path = 'DB_check.pkl'
    pickle.dump(hof_indices_DB_checked, open(DB_check_path, 'w'))
    

    logger.debug('%s out of %s passed the DB check'%(len(hof_indices_DB_checked), 
                          len(hof)))
    
    
    try:
       hof_index = hof_indices_DB_checked[0]            
       logger.debug('The population index which passed DB check with minimum error is ' 
                  'index = %s'%hof_index) 
       return hof_indices_DB_checked[0]
    except IndexError:
       return None
    

                    
def post_processing(checkpoint_file, responses_filename):
    
    
###################### Spikerate comparison between model and experiment ######################    
    
    # Reading the stimulus set from the data
    # Calculating Spikerate only for LongDC (1s step currents) 
    
    stim_map_filename = 'preprocessed/StimMapReps.csv'
    accept_stimtype_list = ['LongDC']
    stim_map = defaultdict(dict)
    with open(stim_map_filename, 'r') as stim_map_file:
        stim_map_content = stim_map_file.read()
        
    for line in stim_map_content.split('\n')[1:-1]:
        if line is not '':
            stim_name, stim_type, holding_current, amplitude_start, amplitude_end, \
                stim_start, stim_end, duration, sweeps = line.split(',')
            if any(stim_type in stim_name for stim_type in accept_stimtype_list):
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
    
    feature_mean_exp = defaultdict()
    stim_name_exp = defaultdict()
    
    somatic_features = ['Spikecount']
    spike_stim_keys_dict = {}
    
    for stim_name in stim_map.keys():
        stim_start = stim_map[stim_name]['stimuli'][0]['delay']
        stim_end = stim_map[stim_name]['stimuli'][0]['stim_end']
        sweeps = []
        for sweep_filename in stim_map[stim_name]['stimuli'][0]['sweep_filenames']:
            sweep_fullpath = os.path.join(
                'preprocessed',
                sweep_filename)
    
            data = np.loadtxt(sweep_fullpath)
            time = data[:, 0]
            voltage = data[:, 1]
    
            # Correct LJP
            voltage = voltage 
            time = time
    
            # Prepare sweep for eFEL
            sweep = {}
            sweep['T'] = time
            sweep['V'] = voltage
            sweep['stim_start'] = [stim_map[stim_name]['stimuli'][0]['delay']]
            sweep['stim_end'] = [stim_map[stim_name]['stimuli'][0]['stim_end']]
    
            sweeps.append(sweep)
        feature_results = efel.getFeatureValues(sweeps, somatic_features)
        for feature_name in somatic_features:
            feature_values = [np.mean(trace_dict[feature_name])
                              for trace_dict in feature_results
                              if trace_dict[feature_name] is not None]
            
        if feature_values:
            feature_mean = feature_values[0]
            
        else:
            feature_mean = 0
        if feature_mean != 0:
            spike_stim_keys_dict[stim_name] = stim_map[stim_name]['stimuli'][0]['amp']
        
        stim_amp =stim_map[stim_name]['stimuli'][0]['amp']
        stim_dur = (float(stim_end) - float(stim_start))/1e3  # in seconds
        feature_mean_exp[stim_amp] = feature_mean/stim_dur
        stim_name_exp[stim_amp] = stim_name
    
    stim_exp = sorted(stim_name_exp.keys())
    mean_freq_exp = [feature_mean_exp[amp] for amp in stim_exp]
        
    # Calculating the spikerate for the model
    
    with open(responses_filename, 'r') as fd:
        opt_responses = pickle.load(fd)
    opt_responses = opt_responses[0]
        
    feature_mean_model_dict = defaultdict()
    stim_model_dict = defaultdict()
    
    
    for key,val in opt_responses.items():
        if 'soma' in key:
            if any(stim_type in key for stim_type in accept_stimtype_list):
                stim_name = key.split('.')[0]
                resp_time = val['time'].as_matrix()
                resp_voltage = val['voltage'].as_matrix()
                stim_amp = stim_map[stim_name]['stimuli'][0]['amp']
                trace1 = {}
                trace1['T'] = resp_time
                trace1['V'] = resp_voltage
                trace1['stim_start'] = [stim_map[stim_name]['stimuli'][0]['delay']]
                trace1['stim_end'] = [stim_map[stim_name]['stimuli'][0]['stim_end']]
                model_traces = [trace1]
                feature_results_model = efel.getFeatureValues(model_traces, somatic_features)
                for feature_name in somatic_features:
                    feature_values_model = [np.mean(trace_dict[feature_name])
                                      for trace_dict in feature_results_model
                                      if trace_dict[feature_name] is not None]
                
                if feature_values_model:
                    feature_mean_model = feature_values_model[0]
                else:
                    feature_mean_model = 0
                
                stim_start = stim_map[stim_name]['stimuli'][0]['delay']
                stim_end = stim_map[stim_name]['stimuli'][0]['stim_end']
                stim_dur = (float(stim_end) - float(stim_start))/1e3
                feature_mean_model_dict[stim_amp] = feature_mean_model/stim_dur
                stim_model_dict[stim_amp] = stim_name
            
    select_stim_keys_list = sorted(spike_stim_keys_dict, key=spike_stim_keys_dict.__getitem__)
    if len(select_stim_keys_list) > 3:
        select_stim_keys =  [select_stim_keys_list[1], select_stim_keys_list[-2]]
    elif len(select_stim_keys_list) > 2:
        select_stim_keys =  [select_stim_keys_list[0], select_stim_keys_list[-2]]  
        
    stim_model = sorted(stim_model_dict.keys())
    mean_freq_model = [feature_mean_model_dict[amp] for amp in stim_model]
    
    max_freq = max([max(mean_freq_model), max(mean_freq_exp)])        
    
    plt.style.use('ggplot')
    fig,ax= plt.subplots(1,figsize=(5,5),dpi =80)   
    
    for key,val in stim_model_dict.items():
        if val in select_stim_keys:
            if val == select_stim_keys[0]:
                shift = (max_freq-feature_mean_model_dict[key])/2 
            else:
                shift = -feature_mean_model_dict[key]/2
            ax.annotate(val,xy=(key,feature_mean_model_dict[key]), 
                        xytext =(key,feature_mean_model_dict[key]+shift), rotation=90,
                        fontsize = 10,arrowprops=dict(facecolor='red', shrink=0.01,alpha =.5))
    ax.scatter(stim_exp, mean_freq_exp,color = 'black',s=50, alpha = .8, label='Experiment')
    ax.plot(stim_exp, mean_freq_exp,color = 'black',lw =.1, alpha = .1)

    ax.scatter(stim_model, mean_freq_model,color = 'blue',s=100,alpha = .9, marker = '*',label='Model')
    ax.plot(stim_model, mean_freq_model,color = 'blue',lw = .1, alpha = .1)

    ax.set_xlabel('Stimulation Amplitude (nA)',fontsize = 10)
    ax.set_ylabel('Spikes/sec',fontsize = 10)
    ax.legend()

    pdf_pages.savefig(fig)
    plt.close(fig)
    
    
 ################### Spike shape comparison between model and experiment ######################

    spike_features = ['peak_time']
    
    fig,ax= plt.subplots(1,2,figsize=(8,5),dpi=80)   
    for kk,stim_name in enumerate(select_stim_keys):
        sweeps = []
        for sweep_filename in stim_map[stim_name]['stimuli'][0]['sweep_filenames']:
            sweep_fullpath = os.path.join(
                'preprocessed',
                sweep_filename)
    
            data = np.loadtxt(sweep_fullpath)
            time = data[:, 0]
            voltage = data[:, 1]
    
            # Prepare sweep for eFEL
            sweep = {}
            sweep['T'] = time
            sweep['V'] = voltage
            sweep['stim_start'] = [stim_map[stim_name]['stimuli'][0]['delay']]
            sweep['stim_end'] = [stim_map[stim_name]['stimuli'][0]['stim_end']]               
            sweeps.append(sweep)
        
        
        model_sweeps = []    
        model_sweep = {}
        name_loc = stim_name+'.soma.v'
        resp_time = opt_responses[name_loc]['time'].as_matrix()
        resp_voltage = opt_responses[name_loc]['voltage'].as_matrix()
        model_sweep['T'] = resp_time
        model_sweep['V'] = resp_voltage
        model_sweep['stim_start'] = [stim_map[stim_name]['stimuli'][0]['delay']]
        model_sweep['stim_end'] = [stim_map[stim_name]['stimuli'][0]['stim_end']]
        model_sweeps.append(model_sweep)
        
        AP_shape_time = np.arange(-5,10,.005) #5ms and 10ms
        AP_shape_exp = np.zeros(AP_shape_time.size)
        AP_shape_model = np.zeros(AP_shape_time.size)
        feature_results = efel.getFeatureValues(sweeps, spike_features)
        feature_results_model = efel.getFeatureValues(model_sweeps, spike_features)
        
        num_spikes_exp = 0
        for k,sweep_filename in enumerate(stim_map[stim_name]['stimuli'][0]['sweep_filenames']):
            sweep_fullpath = os.path.join(
                'preprocessed',
                sweep_filename)
    
            data = np.loadtxt(sweep_fullpath)
            time = data[:, 0]
            voltage = data[:, 1]
            spike_times_exp = feature_results[k]['peak_time']
            for i,spike_time in enumerate(spike_times_exp):
                min_index = np.argmax(time >= spike_time - 10.5) #5ms
                max_index = np.argmax(time >= spike_time +15.5) #10ms
                AP_shape = voltage[min_index:max_index]
                t_shape  = time[min_index:max_index] - spike_time
                f_shape_exp = interpolate.interp1d(t_shape, AP_shape)
                AP_shape_exp += f_shape_exp(AP_shape_time) - f_shape_exp(AP_shape_time[0])
    
            num_spikes_exp += len(spike_times_exp)
        AP_shape_exp /= num_spikes_exp
        
        spike_times_model = feature_results_model[0]['peak_time']
    
        for i,spike_time_model in enumerate(spike_times_model):
            min_index = np.argmax(resp_time >= spike_time_model - 10.5) #5ms
            max_index = np.argmax(resp_time >= spike_time_model +15.5) #10ms
            AP_shape = resp_voltage[min_index:max_index]
            t_shape  = resp_time[min_index:max_index] - spike_time_model
            f_shape_model = interpolate.interp1d(t_shape, AP_shape)
            AP_shape_model += f_shape_model(AP_shape_time)- f_shape_model(AP_shape_time[0])
    
           
        num_spikes_model = len(spike_times_model)
        AP_shape_model /= num_spikes_model
            
            
    
        ax[kk].plot(AP_shape_time, AP_shape_exp,lw = 2, color = 'black',label = 'Experiment')
        ax[kk].plot(AP_shape_time, AP_shape_model,lw = 2, color = 'blue',label = 'Model')
        ax[kk].legend(prop={'size': 10})
        ax[kk].set_title(stim_name,fontsize = 12)
    
    pdf_pages.savefig(fig)
    plt.close(fig)
    pdf_pages.close()

    
     
 