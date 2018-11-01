#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 16:36:19 2018

@author: anin
"""

from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.biophysical_api import BiophysicalApi
from collections import defaultdict
import json
import glob
import os
import errno
import pickle
import bluepyopt as bpopt
import logging
import bluepyopt.ephys as ephys
import pandas as pd
import numpy as np
import matplotlib
import math
import efel
from scipy import interpolate



matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages




logging.basicConfig(level=logging.DEBUG) 
logger = logging.getLogger()
import evaluator_helper

section_map = {'soma' : 'somatic',
                         'apic':'apical',
                         'dend':'basal',
                         'axon':'axonal',
                         'all' : 'all'}

rev_potential = {'ena' : 53, 'ek' : -107}

with open('config_file.json') as json_file:  
    path_data = json.load(json_file)       

morph_path = path_data['morphology']
all_protocol_path = path_data['all_protocols']
protocol_path = path_data['protocols']
feature_path = path_data['features']

with open(protocol_path) as json_file:  
    train_protocols = json.load(json_file)

def plot_Response_Peri(opt,responses_filename,response_peri_filename,pdf_pages):
    
    stim_file = 'preprocessed/StimMapReps.csv'
    stim_df = pd.read_csv(stim_file, sep='\s*,\s*',
                           header=0, encoding='ascii', engine='python')

    responses = pickle.load(open(responses_filename, "r"))
    response = responses[0]  # get the response with minimum trainin error
    responses_release = pickle.load(open(response_peri_filename, "r"))
    logger.debug('Retrieving Optimized and Released Responses')
    
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
                        label= 'All-active',
                        alpha = 0.8)                    

                FileName = 'preprocessed/' + name
                data = np.loadtxt(FileName) 
                l3, = ax_comp[index/n_col,index%n_col].plot(data[:,0],
                            data[:,1],
                            color='black',
                            linewidth=1,
                            label = 'Experiment',
                            alpha = 0.8)
                
                responses_release_time = responses_release[name_loc]['time']
                responses_release_voltage = responses_release[name_loc]['voltage']
                l4,=ax_comp[index/n_col,index%n_col].plot(responses_release_time,
                        responses_release_voltage,
                        color='r',
                        linewidth=.1,
                        label = 'Perisomatic',
                        alpha = 0.4)  
            

                    
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
                    handles = [l1, l3, l4]
                    
                    labels = [h.get_label() for h in handles] 
                    fig_comp.legend(handles = handles, labels=labels, loc = 'lower center', ncol=3)
                    fig_comp.tight_layout(rect=[0, 0.03, 1, 0.95])
                    pdf_pages.savefig(fig_comp)
                    plt.close(fig_comp)
                    index = 0
                    fig_index += 1


def fI_curve_generator(responses_filename,response_peri_filename,perisomatic_model_id,pdf_pages):
    
    # Reading the stimulus set from the data
    # Calculating Spikerate only for LongDC (1s step currents) 
    
    stim_map_filename = 'preprocessed/StimMapReps.csv'
    reject_stimtype_list = ['LongDCSupra','Ramp', 'ShortDC']
    stim_map = defaultdict(dict)
    
    fI_curve_data = {}
    
    with open(stim_map_filename, 'r') as stim_map_file:
        stim_map_content = stim_map_file.read()
        
    for line in stim_map_content.split('\n')[1:-1]:
        if line is not '':
            stim_name, stim_type, holding_current, amplitude_start, amplitude_end, \
                stim_start, stim_end, duration, sweeps = line.split(',')
            if not any(stim_type in stim_name for stim_type in reject_stimtype_list):
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
    
    fI_curve_data['stim_exp'] = stim_exp
    fI_curve_data['freq_exp'] = mean_freq_exp
    
    # Calculating the spikerate for the models
    
    response_dict = {responses_filename:'All_active', 
                                 response_peri_filename : 'Perisomatic'}
    max_freq_model = 0
    
    for model_output,model_type in response_dict.items():
        
        if model_type == 'Perisomatic' and perisomatic_model_id == '':
            continue
        with open(model_output, 'r') as fd:
            opt_responses = pickle.load(fd)
        
        if model_type == 'All_active':
            opt_responses = opt_responses[0]
        
            
        feature_mean_model_dict = defaultdict()
        stim_model_dict = defaultdict()
        
        
        for key,val in opt_responses.items():
            if 'soma' in key:
                if not any(stim_type in key for stim_type in reject_stimtype_list):
                    stim_name = key.split('.')[0]
                    if 'DB' in stim_name:
                        continue
                    resp_time = val['time'].values
                    resp_voltage = val['voltage'].values
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
        if len(select_stim_keys_list) > 2:
            select_stim_keys =  [select_stim_keys_list[0], select_stim_keys_list[-2]]
            
        stim_model = sorted(stim_model_dict.keys())
        mean_freq_model = [feature_mean_model_dict[amp] for amp in stim_model]
        
        fI_curve_data['stim_' + model_type] = stim_model
        fI_curve_data['freq_' + model_type] = mean_freq_model
        
        if max(mean_freq_model) > max_freq_model:
            max_freq_model = max(mean_freq_model)

    
    max_freq = max([max_freq_model, max(mean_freq_exp)])        
    
    plt.style.use('ggplot')
    fig,ax= plt.subplots(1,figsize=(5,5),dpi =80)   
    
    for key,val in stim_model_dict.items():
        if val in select_stim_keys:
            if val == select_stim_keys[0]:
                shift = (max_freq+feature_mean_model_dict[key])/2 
            else:
                shift = -feature_mean_model_dict[key]/2
            ax.annotate(val,xy=(key,feature_mean_model_dict[key]), 
                        xytext =(key,feature_mean_model_dict[key]+shift), rotation=90,
                        fontsize = 10,arrowprops=dict(facecolor='red', shrink=0.01,alpha =.5))
    ax.scatter(stim_exp, mean_freq_exp,color = 'black',s=50, alpha = .8, label='Experiment')
    ax.plot(stim_exp, mean_freq_exp,color = 'black',lw =.1, alpha = .1)
    
    color_dict = {'All_active':'blue','Perisomatic':'magenta'}
    
    for model_type in response_dict.values():
        if model_type == 'Perisomatic' and perisomatic_model_id == '':
            continue
        ax.scatter(fI_curve_data['stim_' + model_type], fI_curve_data['freq_' + model_type],
                   color = color_dict[model_type],s=100,alpha = .9,
                   marker = '*',label= model_type)
        ax.plot(fI_curve_data['stim_' + model_type], fI_curve_data['freq_' + model_type],
                color = color_dict[model_type],lw = .1, alpha = .1)

    ax.set_xlabel('Stimulation Amplitude (nA)',fontsize = 10)
    ax.set_ylabel('Spikes/sec',fontsize = 10)
    ax.legend()

    pdf_pages.savefig(fig)
    plt.close(fig)
    return fI_curve_data, select_stim_keys,stim_map,response_dict


def get_spike_shape(select_stim_keys,stim_map,response_dict,perisomatic_model_id,pdf_pages):

    spike_features = ['peak_time']
    spike_shape_data = {}
    
    fig,ax= plt.subplots(1,2,figsize=(8,5),dpi=80) 
    color_dict = {'All_active':'blue','Perisomatic':'magenta'}
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
            
        
        feature_results = efel.getFeatureValues(sweeps, spike_features)
        AP_shape_time = np.arange(-2,5,.005) #2ms and 5ms
        AP_shape_exp = np.zeros(AP_shape_time.size)
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
                min_index = np.argmax(time >= spike_time - 10.5) #2ms
                max_index = np.argmax(time >= spike_time +15.5) #10ms
                AP_shape = voltage[min_index:max_index]
                t_shape  = time[min_index:max_index] - spike_time
                f_shape_exp = interpolate.interp1d(t_shape, AP_shape)
                AP_shape_exp += f_shape_exp(AP_shape_time)- \
                        f_shape_exp(AP_shape_time[0])
    
            num_spikes_exp += len(spike_times_exp)
        AP_shape_exp /= num_spikes_exp
        
        spike_shape_data[stim_name+'_exp'] = AP_shape_exp
        
        ax[kk].plot(AP_shape_time, AP_shape_exp,lw = 2, color = 'black',label = 'Experiment')

        
        for model_output,model_type in response_dict.items():
            if model_type == 'Perisomatic' and perisomatic_model_id == '':
                continue
            
            with open(model_output, 'r') as fd:
                opt_responses = pickle.load(fd)
            
            if model_type == 'All_active':
                opt_responses = opt_responses[0]
        
            model_sweeps = []    
            model_sweep = {}
            name_loc = stim_name+'.soma.v'
            resp_time = opt_responses[name_loc]['time'].values
            resp_voltage = opt_responses[name_loc]['voltage'].values
            model_sweep['T'] = resp_time
            model_sweep['V'] = resp_voltage
            model_sweep['stim_start'] = [stim_map[stim_name]['stimuli'][0]['delay']]
            model_sweep['stim_end'] = [stim_map[stim_name]['stimuli'][0]['stim_end']]
            model_sweeps.append(model_sweep)
        
           
            AP_shape_model = np.zeros(AP_shape_time.size)
            
            feature_results_model = efel.getFeatureValues(model_sweeps, spike_features)
        
        
        
            spike_times_model = feature_results_model[0]['peak_time']
        
            for i,spike_time_model in enumerate(spike_times_model):
                min_index = np.argmax(resp_time >= spike_time_model - 10.5) #2ms
                max_index = np.argmax(resp_time >= spike_time_model +15.5) #10ms
                AP_shape = resp_voltage[min_index:max_index]
                t_shape  = resp_time[min_index:max_index] - spike_time_model
                f_shape_model = interpolate.interp1d(t_shape, AP_shape)
                AP_shape_model += f_shape_model(AP_shape_time)- \
                                    f_shape_model(AP_shape_time[0])
    
           
            num_spikes_model = len(spike_times_model)
            AP_shape_model /= num_spikes_model
            spike_shape_data[stim_name+'_' +model_type] = AP_shape_model
    
            ax[kk].plot(AP_shape_time, AP_shape_model,lw = 2, 
              color = color_dict[model_type],label = model_type)
        ax[kk].legend(prop={'size': 10})
        ax[kk].set_title(stim_name,fontsize = 12)
        spike_shape_data[stim_name+'_time'] = AP_shape_time
    
    pdf_pages.savefig(fig)
    plt.close(fig)
    pdf_pages.close()
    return spike_shape_data


def Main():
    path_to_cell_metadata = os.path.abspath(os.path.join('.', os.pardir)) + '/cell_metadata.json'        
    with open(path_to_cell_metadata,'r') as metadata:
        cell_metadata = json.load(metadata)
    cell_id = cell_metadata['Cell_id']
    perisomatic_model_id = cell_metadata['Perisomatic_id']
    peri_response_filename = 'resp_peri.txt'
    all_active_response_filename = 'resp_opt.txt'
    analysis_write_path = cell_id + '_peri_comparison.pdf'
    pdf_pages =  PdfPages(analysis_write_path)    
    
    if perisomatic_model_id != '':
        perisomatic_model_id = int(float(cell_metadata['Perisomatic_id']))
        
            
        bp = BiophysicalApi()
        bp.cache_data(perisomatic_model_id,working_directory='peri_model')
        
        
        
        
        peri_param_path = glob.glob('./peri_model/*_fit.json')[0]
        
        
        with open(peri_param_path) as json_file:  
            peri_params = json.load(json_file)
                    
        peri_params_release = list()
        peri_mechs_release = defaultdict(list)
        peri_mechs_release['all'].append('pas')
        
        
        for key, values in peri_params.iteritems():            
            if key == 'genome':
                for j in range(len(values)):
                     iter_dict_release = {'param_name':peri_params[key][j]['name']}
                     iter_dict_release['sectionlist'] = section_map[peri_params[key][j]['section']]
                     iter_dict_release['type'] = 'section'
                     iter_dict_release['value'] = float(peri_params[key][j]['value'])
                     iter_dict_release['dist_type'] = 'uniform'
                     if peri_params[key][j]['mechanism'] != '':
                            iter_dict_release['mech'] = peri_params[key][j]['mechanism']
                            iter_dict_release['type'] = 'range'
                     peri_params_release.append(iter_dict_release)
                         
            elif key == 'passive':
                 for key_pas,val_pas in values[0].items():
                     if key_pas == 'cm':
                         for pas_param in val_pas:
                             iter_dict_release ={'param_name':'cm',
                                                 'sectionlist' : section_map[pas_param['section']],
                                                 'value' : pas_param['cm'],
                                                 'dist_type' : 'uniform',
                                                 'type' : 'section'}
                             peri_params_release.append(iter_dict_release)
                             
                     elif key_pas == 'ra':
                          iter_dict_release ={'param_name':'Ra',
                                                 'sectionlist' : 'all',
                                                 'value' : val_pas,
                                                 'dist_type' : 'uniform',
                                                 'type' : 'section'}
                          peri_params_release.append(iter_dict_release)
                     else:
                          iter_dict_release ={'param_name':key_pas,
                                                 'sectionlist' : 'all',
                                                 'value' : val_pas,
                                                 'dist_type' : 'uniform',
                                                 'type' : 'section'}
                          peri_params_release.append(iter_dict_release)
                             
                             
                         
                     
        
        for rev in rev_potential:
            iter_dict_release =  {'param_name':rev, 'sectionlist':'somatic', 
                                  'dist_type': 'uniform', 'type':'section'}
            if rev == 'ena':
                iter_dict_release['value'] = rev_potential[rev]
            elif rev == 'ek':
                iter_dict_release['value'] = rev_potential[rev]
            peri_params_release.append(iter_dict_release) 
        
        peri_params_release.append({"param_name": "celsius","type": "global","value": 34})     
        peri_params_release.append({"param_name": "v_init","type": "global",
                                     "value": peri_params['conditions'][0]["v_init"]})
                
        for param_dict in peri_params_release:
            if 'mech' in param_dict.keys():
                if param_dict['mech'] not in peri_mechs_release[param_dict['sectionlist']]:
                    peri_mechs_release[param_dict['sectionlist']].append(param_dict['mech'])         
                
                
                
        mechanism_write_path = 'perisomatic_config/' + 'mechanism.json'
        parameters_write_path = 'perisomatic_config/' + 'parameters.json'
        
        if not os.path.exists(os.path.dirname(mechanism_write_path)):
            try:
                os.makedirs(os.path.dirname(mechanism_write_path))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        with open(mechanism_write_path, 'w') as outfile:
            json.dump(peri_mechs_release, outfile,indent=4) 
        
        with open(parameters_write_path, 'w') as outfile:
            json.dump(peri_params_release, outfile,indent=4)   
        
        
        
        evaluator_peri =  evaluator_helper.create(all_protocol_path, feature_path, morph_path, 
                                                    parameters_write_path, mechanism_write_path)
        
        opt_peri = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator_peri)
        
        
        
        
        if os.path.exists(peri_response_filename):
            logger.debug('Retrieving Released Responses')
        else:
            logger.debug('Calculating Released Responses')
            fitness_protocols = opt_peri.evaluator.fitness_protocols
            responses_peri = {}
            nrn = ephys.simulators.NrnSimulator()
            release_params_original = {} #to compare against the model released on the website
            
            for protocol in fitness_protocols.values():
                response_release = protocol.run(
                        cell_model=opt_peri.evaluator.cell_model,
                        param_values=release_params_original,
                        sim=nrn)
                responses_peri.update(response_release)
            
            if peri_response_filename:    
                with open(peri_response_filename, 'w') as fd:
                    pickle.dump(responses_peri, fd)
        
        logger.debug('Plotting Response Comparisons')
        
        plot_Response_Peri(opt_peri,all_active_response_filename,peri_response_filename,
                                               pdf_pages)
        
    fI_curve_data, select_stim_keys,\
                stim_map,response_dict = fI_curve_generator(all_active_response_filename,
                                    peri_response_filename,perisomatic_model_id,pdf_pages)
    
    spike_shape_data = get_spike_shape(select_stim_keys,stim_map,
                                       response_dict,perisomatic_model_id,pdf_pages)  
    
    fI_curve_path = 'perisomatic_comparison/fI_curve_' +str(cell_id)+'.pkl'
    spike_shape_path = 'perisomatic_comparison/spike_shape_' +str(cell_id)+'.pkl'
    
    if not os.path.exists(os.path.dirname(fI_curve_path)):
        try:
            os.makedirs(os.path.dirname(fI_curve_path))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    with open(fI_curve_path, 'w') as outfile:
        pickle.dump(fI_curve_data, outfile) 
    
    with open(spike_shape_path, 'w') as outfile:
        pickle.dump(spike_shape_data, outfile)
        
            
if __name__ == '__main__':
    Main()