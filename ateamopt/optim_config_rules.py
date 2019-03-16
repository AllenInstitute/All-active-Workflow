
def filter_feat_proto_active(features_dict,protocols_dict):
    """ 
        Filter the features and protocols for the final
        stage of optimization
    """
    pass



def filter_feat_proto_passive(features_dict,protocols_dict):
    
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
    return features_dict_filtered, untrained_features_dict, protocols_dict_filtered


def adjust_param_bounds(model_param, model_param_prev,tolerance=0.5):
    lb_,ub_ = model_param['bounds']
    value = model_param_prev['value']
    lb = max(value - tolerance*abs(value),lb_)
    ub = min(value + tolerance*abs(value),ub_)
    adjusted_bound = [lb,ub]
    model_param['bounds'] = adjusted_bound
    return model_param