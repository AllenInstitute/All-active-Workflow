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
import matplotlib

matplotlib.use('Agg')

import matplotlib.pyplot as plt
import errno
import logging
from matplotlib.backends.backend_pdf import PdfPages
import math
from shutil import copyfile
import glob
from collections import OrderedDict,defaultdict



logger = logging.getLogger(__name__)


with open('config_file.json') as json_file:  
    data = json.load(json_file)

section_map_inv = {'somatic':'soma', 'axonal':'axon', 'apical':'apic',
               'basal':'dend', 'all':'all'}


release_params = data['release_params']        
fit_json_path = data['fit_json']
fit_protocol_path = data['protocols']
param_path = data['parameters']

with open(fit_protocol_path) as json_file:  
    train_protocols = json.load(json_file)
    
if fit_json_path:
    with open(fit_json_path) as json_file:  
        model_data = json.load(json_file)

with open('passive_param_bounds.json','r') as bound_file:
        passive_params_dict = json.load(bound_file)
 
with open(param_path) as json_file:  
    params = json.load(json_file)
    
parent_dir = os.path.abspath(os.path.join('.', os.pardir))
path_to_cell_metadata = glob.glob(parent_dir+'/*.json')[0] 
with open(path_to_cell_metadata,'r') as metadata:
        cell_metadata = json.load(metadata)
    
cell_id = cell_metadata['Cell_id']
layer = cell_metadata['Layer']
area = cell_metadata['Area']
species = cell_metadata['Species']    
cre_line = cell_metadata['Cre_line']
dendrite_type = cell_metadata['Dendrite_type']


analysis_write_path = cell_id + '_analysis_Stage0.pdf'
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
        'cell_id' : cell_id,
        'layer' : layer,
        'area' : area,
        'species' : species,
        'cre_line' : cell_metadata['Cre_line'],
        'dendrite_type' : cell_metadata['Dendrite_type']})
    
        hof_df=hof_df.append(temp_df) 
    
    if release_params:    
        for index, param_name in enumerate(param_names_arranged):
            release_individual[index] = release_params[param_name]
        abs_release_individual = map(abs,release_individual)
            

    fig, ax = plt.subplots(1,figsize=(6,6))
    
    x = np.arange(len(optimized_individual_arranged))
    
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
    
    optimized_param_dict = {key:optimized_individual_arranged[i] for i,key in \
                            enumerate(param_names_arranged)} 
    
    
    fit_json_write_path = 'fitted_params/optim_param_'+cell_id+ '.json'
    fit_json_write_path_2 = './fit_opt.json'
    
    if not os.path.exists(os.path.dirname(fit_json_write_path)):
        try:
            os.makedirs(os.path.dirname(fit_json_write_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    optim_param_write_path = 'fitted_params/optim_param_unformatted_' + cell_id +'.json'
   
    copyfile(fit_json_write_path_2, fit_json_write_path)
        
    with open(optim_param_write_path, 'w') as outfile:
        json.dump(optimized_param_dict, outfile,indent=4)
        
    logger.debug('Saving the parameters in fit.json format')
    
   # save parameters in a csv file to later plot in R
            
    

    optimized_df = pd.DataFrame({'param_name' : param_names_arranged,
        'section':sect_names_arranged,                         
        'value' : optimized_individual_arranged,
        'label' : 'Optimized',
        'fitness': 10,
        'cell_id' : cell_id,
        'layer' : layer,
        'area' : area,
        'species' : species,
        'cre_line' : cell_metadata['Cre_line'],
        'dendrite_type' : cell_metadata['Dendrite_type']})
    
    if release_params:
        released_df = pd.DataFrame({ 'param_name' : param_names_arranged,
            'section':sect_names_arranged,                        
            'value' : release_individual,
            'label' : 'Released',
            'fitness': 5,
            'cell_id' : cell_id,
            'layer' : layer,
            'area' : area,
            'species' : species,
            'cre_line' : cell_metadata['Cre_line'],
            'dendrite_type' : cell_metadata['Dendrite_type']})
            
        param_df = [optimized_df, released_df,hof_df] 
    else:
        param_df = [optimized_df,hof_df] 
        
    param_df = pd.concat(param_df)   
    param_df.to_csv('params_'+cell_id+'.csv')
    
    logger.debug('Saving the parameters in .csv for plotting in R')    
    
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

    responses = get_responses(opt.evaluator, [hof[0]], responses_filename)
    responses = responses[0] #select the optimized response     
    
    # objectives       
    objectives = opt.evaluator.fitness_calculator.calculate_scores(responses)
    
    objectives = OrderedDict(sorted(objectives.iteritems()))
    
    feature_split_names = [name.split('.')[-1] for name in objectives.keys()]
    features = np.unique(np.asarray(feature_split_names))
    
    
    bar_width = 0.35
    opacity = 0.4

    plt.style.use('ggplot') 
    for i, feature in enumerate(features):
        fig, ax = plt.subplots(1, figsize=(8,8))    
        iter_dict = defaultdict(list)
#        iter_dict_release = defaultdict(list)
        for key in objectives.keys():
            if key.split('.')[-1] == feature:
                amp = train_protocols[key.split('.')[0]]['stimuli'][0]['amp']
                amp_reduced = round(amp,2)
                iter_dict[str(amp_reduced)].append(objectives[key])
        index = np.arange(len(iter_dict.keys()))
        xtick_pos = index
        iter_dict ={float(key):np.mean(val) for key,val in iter_dict.items()}
#        iter_dict_release ={key:np.mean(val) for key,val in iter_dict_release.items()}
        iter_dict = OrderedDict(sorted(iter_dict.iteritems()))
#        iter_dict_release = OrderedDict(sorted(iter_dict_release.iteritems()))
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
            'cre_line' : cell_metadata['Cre_line'],
            'dendrite_type' : cell_metadata['Dendrite_type'],
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
    all_plots = 0

    
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
            all_plots += len(trace_rep.split('|'))
            
            
    n_col = 3  
    n_row = 5
    fig_per_page = n_col * n_row
    fig_pages =  int(math.ceil(all_plots/float(fig_per_page)))
    fig_mat = list()
    ax_mat = list()
    for page_i in range(fig_pages):
        
        fig,ax = plt.subplots(n_row,n_col, figsize=(10,10),squeeze = False)
        
        if page_i == fig_pages-1:
            remain_plots = all_plots - n_row*n_col*(fig_pages-1)
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
                response_time = response[name_loc]['time']
                response_voltage = response[name_loc]['voltage']
                
                color = 'blue'
                l1, = ax_comp[index/n_col,index%n_col].plot(response_time,
                        response_voltage,
                        color=color,
                        linewidth=1,
                        label= 'Model',
                        alpha = 0.8)  

                FileName = 'preprocessed/' + name
                data = np.loadtxt(FileName)
                exp_time = data[:,0]
                exp_voltage = data[:,1]
                l3, = ax_comp[index/n_col,index%n_col].plot(exp_time,
                            exp_voltage,
                            color='black',
                            linewidth=1,
                            label = 'Experiment',
                            alpha = 0.8)
                
                if index/n_col == n_row-1: 
                    ax_comp[index/n_col,index%n_col].set_xlabel('Time (ms)')
                if index%n_col == 0: 
                    ax_comp[index/n_col,index%n_col].set_ylabel('Voltage (mV)')
                ax_comp[index/n_col,index%n_col].set_title(name.split('.')[0] + state, fontsize=8)
                
                if 'LongDC' in name:
                    ax_comp[index/n_col,index%n_col].set_xlim([amp_start_list[ix]-200,\
                                                              amp_end_list[ix]+200])
                        
                logger.debug('Plotting response comparisons for %s \n'%name.split('.')[0])
                index += 1
                index_plot +=1
                if index%fig_per_page == 0 or index_plot == all_plots:
                    fig_comp.suptitle('Response Comparisons',fontsize=16)
                    
                    handles = [l1, l3]
                    ncol = 2
                    labels = [h.get_label() for h in handles]
                    fig_comp.legend(handles = handles, labels=labels, loc = 'lower center', ncol=ncol)
                    fig_comp.tight_layout(rect=[0, 0.03, 1, 0.95])
                    pdf_pages.savefig(fig_comp)
                    plt.close(fig_comp)
                    index = 0
                    fig_index += 1
           
            

    
    
