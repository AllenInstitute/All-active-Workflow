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
import uncertainpy as un
import bluepyopt.ephys as ephys

logger = logging.getLogger(__name__)


with open('config_file.json') as json_file:  
    data = json.load(json_file)


release_params = data['release_params']        
fit_json_path = data['fit_json']
fit_protocol_path = data['protocols']

with open(fit_protocol_path) as json_file:  
    train_protocols = json.load(json_file)
    
cell_id = 468193142
layer = 'Layer 4'
analysis_write_path = str(cell_id) + '_analysis.pdf'
pdf_pages =  PdfPages(analysis_write_path)

def plot_diversity(opt, checkpoint_file, param_names):

    plt.style.use('ggplot')
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    import copy
    
    release_individual = copy.deepcopy(checkpoint['halloffame'][0])
    optimized_individual = [checkpoint['halloffame'][0]]
#    optimized_fitness = optimized_individual[0].fitness.values


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
        'layer' : layer})
        hof_df=hof_df.append(temp_df) 
        
    for index, param_name in enumerate(param_names_arranged):
        release_individual[index] = release_params[param_name]

#    param_history = checkpoint['history'].genealogy_history.values()
    fig, ax = plt.subplots(1,figsize=(6,6))
    
    param_count = len(release_individual)
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
    
    optimized_param_dict = {key:optimized_individual_arranged[i] for i,key in \
                            enumerate(param_names_arranged)} 
    section_map = {'somatic':'soma', 'axonal':'axon', 'apical':'apic', 'basal':'dend', 'all':'all'}
    
    param_dict_final = {key.split('.')[0]+'.'+
                     section_map[key.split('.')[1]] : optimized_param_dict[key] 
                                            for key in optimized_param_dict.keys()} 
    
    with open(fit_json_path) as json_file:  
        model_data = json.load(json_file)
    
    for key in param_dict_final.keys():
        opt_name,opt_sect = key.split('.')
        data_key = 'genome'
        repeat_list = list()
        remove_indices = list()
        for j in range(len(model_data[data_key])):

            if model_data[data_key][j]['name'] == opt_name:
               if opt_name not in repeat_list:                     
                   model_data[data_key][j]['value'] = str(param_dict_final[key])
                   model_data[data_key][j]['section'] = 'all'
                   repeat_list.append(opt_name)
               else:
                   remove_indices.append(j)
        model_data[data_key] = [i for j, i in enumerate(model_data[data_key]) if j not in remove_indices]
    fit_json_write_path = 'optim_param.json'
    optim_param_write_path = 'optim_param_unformatted.json'
    with open(fit_json_write_path, 'w') as outfile:
        json.dump(model_data, outfile,indent=4)
        
    with open(optim_param_write_path, 'w') as outfile:
        json.dump(optimized_param_dict, outfile,indent=4)
        
    print 'Saving the parameters in fit.json format'
    
   # save parameters in a csv file to later plot in R
            
    released_df = pd.DataFrame({ 'param_name' : param_names_arranged,
        'section':sect_names_arranged,                        
        'value' : release_individual,
        'label' : 'Released',
        'fitness': 5,
        'cell_id' : cell_id,
        'layer' : layer})
    
    optimized_df = pd.DataFrame({'param_name' : param_names_arranged,
        'section':sect_names_arranged,                         
        'value' : optimized_individual_arranged,
        'label' : 'Optimized',
        'fitness': 10,
        'cell_id' : cell_id,
        'layer' : layer})
    param_df = [optimized_df, released_df,hof_df] 
    param_df = pd.concat(param_df)   
    param_df.to_csv('params.csv')
    
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
    if filename and os.path.exists(filename):
        with open(filename) as fd:
            return pickle.load(fd)

    for individual in individuals:
        individual_dict = cell_evaluator.param_dict(individual)
        responses.append(
            cell_evaluator.run_protocols(
                cell_evaluator.fitness_protocols.values(),
                param_values=individual_dict)) 

    if filename:
        with open(filename, 'w') as fd:
            pickle.dump(responses, fd)

    return responses



def feature_comp(opt, checkpoint_file):
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    hof = checkpoint['halloffame']

    # objectives
    parameter_values = opt.evaluator.param_dict(hof[0])
    fitness_protocols = opt.evaluator.fitness_protocols
    responses = {}
    responses_release = {}
    
    nrn = ephys.simulators.NrnSimulator()
    
    for protocol in fitness_protocols.values():
        response = protocol.run(
            cell_model=opt.evaluator.cell_model,
            param_values=parameter_values,
            sim=nrn)
        response_release = protocol.run(
                cell_model=opt.evaluator.cell_model,
                param_values=release_params,
                sim=nrn)
        responses.update(response)
        responses_release.update(response_release)
    
    with open('release_response.txt', 'w') as fd:
        pickle.dump(responses_release, fd)
        
    objectives = opt.evaluator.fitness_calculator.calculate_scores(responses)
    objectives_release = opt.evaluator.fitness_calculator.calculate_scores(responses_release)
    
    import collections
    objectives = collections.OrderedDict(sorted(objectives.iteritems()))
    objectives_release = collections.OrderedDict(sorted(objectives_release.iteritems()))
    
    feature_split_names = [name.split('.')[-1] for name in objectives.keys()]
    features = np.unique(np.asarray(feature_split_names))
    
    
    bar_width = 0.35
    opacity = 0.4
    feature_draw_path = 'figures/Features/'
    try:
        os.makedirs(os.path.dirname(feature_draw_path))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
    plt.style.use('ggplot') 
#    fig, ax = plt.subplots(1, 3, figsize=(10,4))       
    for i, feature in enumerate(features):
        fig, ax = plt.subplots(1, figsize=(8,8))    
        iter_dict = dict()
        iter_dict_release = dict()
        tick_label = []
        for key in objectives.keys():
            if key.split('.')[-1] == feature:
                iter_dict[key] = objectives[key]
                iter_dict_release[key] = objectives_release[key]
                amp = train_protocols[key.split('.')[0]]['stimuli'][0]['amp']
                amp_reduced = round(amp,2)
                tick_label.append(amp_reduced)    
        index = np.arange(len(iter_dict.keys()))
#        ytick_pos = index + bar_width / 2
#        index = tick_label        
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
              
#        ax.barh(index + bar_width,
#              iter_dict_release.values(),
#              bar_width,
#              align='center',
#              color='r',
#              alpha=opacity,
#              label='Released')  
        ax.set_xticks(xtick_pos)
        ax.set_xticklabels(tick_label, fontsize= 8)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(8) 
        ax.set_ylabel('Objective value (# std)',fontsize= 12)
        ax.set_xlabel('Stimulus Amplitude (nA)',fontsize= 12)        
        ax.legend(prop={'size': 12},loc ='best')
        ax.set_title(feature, fontsize= 14)
#        fig.suptitle('Feature Comparison')
        fig.tight_layout(rect=[0.05, 0.05, .95, 0.95])
        pdf_pages.savefig(fig)
        plt.close(fig)    
       
#        fig.savefig(feature_draw_path+feature+'.png')
        
        logger.debug('Plotting comparisons for %s \n'%feature)
    
    stim_file = 'preprocessed/StimMapReps.csv'
    stim_df = pd.read_csv(stim_file, sep='\s*,\s*',
                           header=0, encoding='ascii', engine='python')
    voltage_deflection_df = pd.DataFrame([])
    for key,val in objectives.iteritems():
        if key.split('.')[-1] == 'voltage_deflection_vb_ssse':
            proto = key.split('.')[0]
            temp_df = pd.DataFrame({'deflection_fit' : val,
                                'stim_amp': stim_df[stim_df['DistinctID'] == proto]['Amplitude_Start'],
                                'cell_id' : cell_id
                                })
            voltage_deflection_df=voltage_deflection_df.append(temp_df) 

    voltage_deflection_df.to_csv('voltage_deflection_fit.csv')

#############################################################################

## Plotting responses


def plot_Response(opt,checkpoint_file, responses_filename):
    stim_file = 'preprocessed/StimMapReps.csv'
    stim_df = pd.read_csv(stim_file, sep='\s*,\s*',
                           header=0, encoding='ascii', engine='python')
    checkpoint = pickle.load(open(checkpoint_file, "r"))
    hof = checkpoint['halloffame']
#    responses_release = pickle.load(open('release_response.txt', "r"))

    response_draw_path = 'figures/Responses/'
    try:
        os.makedirs(os.path.dirname(response_draw_path))
    except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST:
            raise
    responses = get_responses(opt.evaluator, hof, responses_filename)
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
                        label= 'Optimized',
                        alpha = 0.8)
                    

                FileName = 'preprocessed/' + name
                data = np.loadtxt(FileName)
                exp_time = data[:,0]
                exp_voltage = data[:,1]
                l3, = ax[index/n_col,index%n_col].plot(exp_time,
                            exp_voltage,
                            color='black',
                            linewidth=1,
                            label = 'Cell Response',
                            alpha = 0.5)  
                data = None
                ax[index/n_col,index%n_col].legend(handles = [l1,l3],prop={'size': 7})
                if (analysis_state == '(T)' and index/n_col == n_row_train-1) or \
                                (analysis_state == '(V)' and index/n_col == n_row_val-1): 
                    ax[index/n_col,index%n_col].set_xlabel('Time (ms)')
                if (analysis_state == '(T)' and index%n_col == 0) or \
                                (analysis_state == '(V)' and index%n_col == 0): 
                    ax[index/n_col,index%n_col].set_ylabel('Voltage (mV)')
                ax[index/n_col,index%n_col].set_title(name.split('.')[0], fontsize=8)
                logger.debug('Plotting response comparisons for %s \n'%name.split('.')[0])
                index += 1
            if analysis_state == '(T)':
                index_train = index
            else:
                index_val = index
            
    fig_train.suptitle('Training Set',fontsize=16) 
    fig_val.suptitle('Test Set',fontsize=16)       
    fig_val.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig_train.tight_layout(rect=[0, 0.03, 1, 0.95])
    pdf_pages.savefig(fig_train)
    plt.close(fig_train)
    pdf_pages.savefig(fig_val)
    plt.close(fig_val)
    
    
def param_Sensitivity(opt):
    nrn = ephys.simulators.NrnSimulator()
    fitness_protocols = opt.evaluator.fitness_protocols
    key = 'LongDC_56'
    def nrnsim_bpopt(g_pas,e_pas,cm,Ra):
        sensitivity_params = {'g_pas.all' : g_pas,
                              'e_pas.all' : e_pas,
                              'cm.all' : cm,
                              'Ra.all' : Ra
                              }
        sensitivity_response = fitness_protocols[key].run(\
                    cell_model=opt.evaluator.cell_model,
                    param_values=sensitivity_params,
                    sim=nrn)
        name_loc = key + '.soma.v'
        time = np.asarray(sensitivity_response[name_loc]['time'])
        value = np.asarray(sensitivity_response[name_loc]['voltage'])
        info = {'stimulus_start':train_protocols[key]['stimuli'][0]['delay'], 
                'stimulus_end':train_protocols[key]['stimuli'][0]['stim_end']}
        return time, value, info
    
    optim_param_write_path = 'optim_param_unformatted.json'
    with open(optim_param_write_path) as read_file:
        optim_param = json.load(read_file)
    parameters ={ 'g_pas' : optim_param['g_pas.all'],
                  'e_pas' : optim_param['e_pas.all'],
                  'cm' : optim_param['cm.all'],
                  'Ra' : optim_param['Ra.all']                 
                 }
    features_to_run = ["voltage_deflection_vb_ssse",
                   "decay_time_constant_after_stim",
                   "steady_state_voltage",
                   "voltage_base"]
    features = un.EfelFeatures(features_to_run=features_to_run)

    model = un.Model(run=nrnsim_bpopt,adaptive=True,
                 labels=["Time (ms)", "Membrane potential (mV)"])
    parameters = un.Parameters(parameters)
    parameters.set_all_distributions(un.uniform(.1))
    
    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(model,
                                      parameters=parameters,
                                      features=features)
    UQ.quantify()
    
    import matplotlib as mpl
    mpl.rcParams.update(mpl.rcParamsDefault)
    l5pc =  un.Data("data/nrnsim_bpopt.h5")
    time = l5pc["nrnsim_bpopt"].time
    evaluations = l5pc["nrnsim_bpopt"].evaluations
#    sobol_first = l5pc["nrnsim_bpopt"].sobol_first
    variance = l5pc["nrnsim_bpopt"].variance
    mean = l5pc["nrnsim_bpopt"].mean
    percentile_95 = l5pc["nrnsim_bpopt"].percentile_95
    percentile_5 = l5pc["nrnsim_bpopt"].percentile_5
    
    ylabel = 'Membrane Potential (mV)'
    xlabel = 'Time (ms)'
    fig,ax1 = plt.subplots(1,figsize=(6,6))
        
    plt.style.use("ggplot")
    
    l2 = ax1.fill_between(time,
            percentile_5,
            percentile_95,
            color = 'steelblue',
            alpha = 0.8,
            label = '90% prediction interval')
    
    ax2 = ax1.twinx()
    l3, = ax2.plot(time, variance, label = 'Variance',
             alpha = 0.8)
    l1, = ax1.plot(time, mean,color='black',
                   label='Mean')
    ax1.set_xlabel(xlabel, color='black')
    ax1.tick_params('y', colors='black')
     
    ax1.set_ylabel(ylabel)
    
    ax2.set_ylabel('Variance $(\mathrm{mV}^2)$', color='r')
    ax2.tick_params('y', colors='r')
    ax2.grid(False)
    ax1.grid(True)
    
    # Shrink current axis's height by 10% on the bottom
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width*.9, box.height * 0.8])
    box2 = ax2.get_position()
    ax2.set_position([box2.x0, box2.y0 + box2.height * 0.2,
                 box2.width*.9, box2.height * 0.8])    
    ax1.legend(handles = [l1,l2,l3], bbox_to_anchor=(0.5, -0.16),
                   loc='upper center', ncol = 3)
    fig.suptitle('Sensitivity Analysis')
    pdf_pages.savefig(fig)
    plt.close(fig)
    
    import random

    random.seed(1)
    len_evaluations = np.shape(evaluations)[0]
    plot_num = 9
    rand_indices = random.sample(range(len_evaluations),plot_num)
    plot_col = 3
    plot_row = int(math.ceil(plot_num/float(plot_col)))
    fig, axes = plt.subplots(plot_row, plot_col, figsize=(8,6), sharex=True, sharey=True)
    for i in range(plot_num):
        ind_1 = i/plot_col
        ind_2 = i%plot_col
        axes[ind_1,ind_2].plot(time,evaluations[rand_indices[i]], lw =1)
        if ind_2 == 0:
            axes[ind_1,ind_2].set_ylabel('Voltage')
        if ind_1 == 2:
            axes[ind_1,ind_2].set_xlabel(xlabel)   
    fig.suptitle('Example membrane potentials',fontsize=12)    
    fig.tight_layout(rect=[0, 0.05, .95, 0.95])
    pdf_pages.savefig(fig)
    plt.close(fig)
    
    features = l5pc.keys()
    nr_plots = len(features)
    bar_width = 0.75
    opacity = 0.8
    x_labels =l5pc.uncertain_parameters
    
    x = np.arange(len(x_labels))
    plot_row = int(math.ceil(nr_plots/float(plot_col)))
    fig, axes = plt.subplots(plot_row, plot_col, figsize=(8,6))
    fig_empty_index = range(nr_plots,plot_row*plot_col)
    if len(fig_empty_index) != 0:
        for ind in fig_empty_index:
            axes[ind/plot_col,ind%plot_col].axis('off')
    
    axes[-1,-1].axis('off') 
    cmap = plt.get_cmap('RdBu')
    colors = cmap(np.linspace(0, 1, len(x)))
    
    for i in range(nr_plots):
        ind_1 = i/plot_col
        ind_2 = i%plot_col
        axes[ind_1,ind_2].bar(x, l5pc[features[i]].sobol_first_sum, 
            bar_width,align='center',alpha=opacity,color=colors)
    
        axes[ind_1,ind_2].set_xticklabels(x_labels, rotation=90)
        axes[ind_1,ind_2].set_xticks(x)
        feature = features[i]
        if feature == 'nrnsim_bpopt':
            feature = 'l5pc'
        axes[ind_1,ind_2].set_title(feature,fontsize=10)
    fig.suptitle('Normalized Sum of First Order Sobol Indices')        
    fig.tight_layout(rect=[0, 0.03, .95, 0.95])
    pdf_pages.savefig(fig)
    plt.close(fig)

#    pdf_pages.close()