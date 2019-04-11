from collections import OrderedDict
from collections import defaultdict
import numpy as np

import logging
logger = logging.getLogger(__name__)


def filter_feat_proto_active(features_dict,protocols_dict,all_protocols_dict,
                             select_dict, add_fi_kink = True,
                             add_DB_check = True):
    """
        Filter the features and protocols for the final
        stage of optimization
    """
    features_dict = correct_voltage_feat_std(features_dict)
    spiking_proto_dict = OrderedDict()
    non_spiking_proto_dict =  OrderedDict()
    training_stimtype_reject = ['LongDCSupra','Ramp','Short_Square_Triple','Noise']
    
    for feat_key,feat_val in features_dict.items():
        if any(reject_stim in feat_key for reject_stim in training_stimtype_reject):
            continue
        stim_amp = protocols_dict[feat_key]['stimuli'][0]['amp']
        if feat_val['soma']['Spikecount'][0] > 0:
            spiking_proto_dict[feat_key] = stim_amp
            del feat_val['soma']['Spikecount']
        else:
            non_spiking_proto_dict[feat_key] = stim_amp
            if 'depol_block' in feat_val['soma'].keys():
                del feat_val['soma']['depol_block']
    
    
    # Ignoring spiking protocol which are followed by non-spiking stim protocol
    max_nospiking_amp = max(non_spiking_proto_dict.values())            
    for spike_stim,spike_amp in spiking_proto_dict.items():
        if spike_amp < max_nospiking_amp:
           del spiking_proto_dict[spike_stim] 
    spiking_proto_sorted = sorted(spiking_proto_dict,
                           key=spiking_proto_dict.__getitem__)
    non_spiking_proto_sorted = sorted(non_spiking_proto_dict,
                           key=non_spiking_proto_dict.__getitem__)

    num_spike = select_dict['spike_proto'] # In descending order
    num_nospike = select_dict['nospike_proto'] # in descending order

    # Select spiking proto
    try:
        spiking_proto_select = spiking_proto_sorted[len(spiking_proto_sorted)-num_spike:\
                                                    len(spiking_proto_sorted)]
    except:
        logger.debug('Number of spiking protocols requested exceeds data')
        spiking_proto_select = spiking_proto_sorted

    # Select nospiking proto
    try:
        nonspiking_proto_select = non_spiking_proto_sorted[len(non_spiking_proto_sorted)-\
                                       num_nospike:len(non_spiking_proto_sorted)]
    except:
        logger.debug('Number of nonspiking protocols requested exceeds data')
        nonspiking_proto_select = non_spiking_proto_sorted

    if add_fi_kink:
        spiking_proto_select.append(spiking_proto_sorted[0]) # the fist spiking stim
        nonspiking_proto_select.append(non_spiking_proto_sorted[-1]) # the last non-spiking stim
        spiking_proto_select = list(set(spiking_proto_select))
        nonspiking_proto_select = list(set(nonspiking_proto_select))

    features_dict_filtered = {key:val for key,val in features_dict.items() \
                              if key in spiking_proto_select+
                                      nonspiking_proto_select}
    protocols_dict_filtered = {key:val for key,val in protocols_dict.items() \
                              if key in spiking_proto_select+
                                          nonspiking_proto_select}
    untrained_features_dict = {key:val for key,val in features_dict.items() \
                              if key not in spiking_proto_select+
                                  nonspiking_proto_select}
    
    # For fI kink spiking proto only allow the following features
    for u_key,u_val in untrained_features_dict.items():
        u_val['soma'] = entries_to_remove(['depol_block',\
                       'check_AISInitiation','Spikecount'],u_val['soma'])
            
    for f_key,f_val in features_dict_filtered[spiking_proto_sorted[0]]['soma'].items():
        if f_key not in ['mean_frequency', 'depol_block',\
                       'check_AISInitiation']:
            features_dict_filtered[spiking_proto_sorted[0]]['soma'].pop(f_key)
    
    if add_DB_check:
        max_proto_key = spiking_proto_sorted[-1]
        max_amp = max([proto['amp'] for proto in protocols_dict[max_proto_key]['stimuli']])
        DB_proto_delay = max([proto['delay'] for proto in protocols_dict[max_proto_key]['stimuli']])
        DB_proto_duration = min([proto['duration'] for proto in protocols_dict[max_proto_key]['stimuli']])
        DB_proto_stimend = min([proto['stim_end'] for proto in protocols_dict[max_proto_key]['stimuli']])
        DB_proto_totdur = min([proto['totduration'] for proto in protocols_dict[max_proto_key]['stimuli']])

        DB_holding_proto = [proto for proto in protocols_dict[max_proto_key]['stimuli'] \
                        if proto['delay'] == 0]
        DB_proto_dict = [{
                    'type' :'SquarePulse',
                    'amp' : max_amp + 0.01,
                    'delay' : DB_proto_delay,
                    'duration' : DB_proto_duration,
                    'stim_end' : DB_proto_stimend,
                    'totduration': DB_proto_totdur
                    }]
        if bool(DB_holding_proto):
            DB_proto_dict.append(DB_holding_proto[0])
        DB_feature_dict = {
                'depol_block' : [1.0,
                                 0.05]
                          }
        features_dict_filtered['DB_check_DC'] = {'soma': DB_feature_dict}
        protocols_dict_filtered['DB_check_DC'] = {'stimuli':DB_proto_dict}
        all_protocols_dict['DB_check_DC'] = {'stimuli':DB_proto_dict}

    return features_dict_filtered,untrained_features_dict,\
         protocols_dict_filtered,all_protocols_dict


def filter_feat_proto_passive(features_dict,protocols_dict,all_protocols_dict,
                                  *args):

    features_dict = correct_voltage_feat_std(features_dict)
    spiking_proto = []
    for feat_key,feat_val in features_dict.items():
        if feat_val['soma']['Spikecount'][0] > 0:
            spiking_proto.append(feat_key)
        del feat_val['soma']['Spikecount']

    features_dict_filtered = {key:val for key,val in features_dict.items() \
                              if key not in spiking_proto}
    protocols_dict_filtered = {key:val for key,val in protocols_dict.items() \
                              if key not in spiking_proto}
    untrained_features_dict = {key:val for key,val in features_dict.items() \
                              if key not in features_dict_filtered.keys()}
    return features_dict_filtered,untrained_features_dict,\
                 protocols_dict_filtered,all_protocols_dict

def correct_voltage_feat_std(features_dict):
    feature_correct_list = ['voltage_base', 'steady_state_voltage']
    feature_stat = defaultdict(list)
    feature_key = []
    for feat_name in feature_correct_list:
        for key,val in features_dict.items():
            if val['soma'][feat_name][1] == 0:
                feature_stat[feat_name].append(val['soma'][feat_name][0])
                feature_key.append(key)

    feature_key = list(set(feature_key))

    for feat_key in feature_key:
        for feat_name in feature_correct_list:
            features_dict[feat_key]['soma'][feat_name][1] = \
                    np.std(feature_stat[feat_name])

    return features_dict

def adjust_param_bounds(model_param, model_param_prev,tolerance=0.5):
    lb_,ub_ = model_param['bounds']
    value = model_param_prev['value']
    lb = max(value - tolerance*abs(value),lb_)
    ub = min(value + tolerance*abs(value),ub_)
    adjusted_bound = [lb,ub]
    model_param['bounds'] = adjusted_bound
    return model_param

def entries_to_remove(entries, the_dict):
    for key in entries:
        if key in the_dict.keys():
            del the_dict[key]
    return the_dict