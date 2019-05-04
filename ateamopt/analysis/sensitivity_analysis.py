import uncertainpy as un
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ateamopt.utils import utility
import pandas as pd
import os
import math
import seaborn as sns


class SA_helper(object):
    
    def __init__(self, optim_param_path, sens_param_path,
                 param_mod_range,config_file):
        self.config = utility.load_json(config_file)
        self.sens_parameters = utility.load_json(sens_param_path)
        self.optim_param = utility.load_json(optim_param_path)
        self.param_range = param_mod_range
    
    def load_config(self,model_basepath=None):
        morph_path = os.path.normpath(self.config['morphology'])
        protocol_path = os.path.normpath(self.config['all_protocols'])
        mech_path = os.path.normpath(self.config['mechanism'])
        feature_path = os.path.normpath(self.config['features'])
        param_path = os.path.normpath(self.config['parameters'])
        
        if model_basepath:
           morph_path = os.path.join(model_basepath,morph_path)
           protocol_path = os.path.join(model_basepath,protocol_path)
           mech_path = os.path.join(model_basepath,mech_path)
           feature_path = os.path.join(model_basepath,feature_path)
           param_path = os.path.join(model_basepath,param_path)
           
        return morph_path,protocol_path,mech_path,feature_path,param_path

    
    def create_sens_param_dict(self):
        param_dict_uc = {}
        for key,val in self.sens_parameters.items():
            for sect in val:
                param_dict_uc['%s_%s'%(key,sect)] = '%s.%s'%(key,sect)
        return param_dict_uc
    
    def create_sa_bound(self,bpopt_param_bounds_path,sens_param_bounds_path,
                        max_bound = .5):
        
        # For parameter sensitivity create a new set of bounds because the permutations may 
        # fall outside the original bounds
        
        max_bound = max(max_bound,self.param_range+.1)
        
        bpopt_section_map = utility.bpopt_section_map
        param_bounds = utility.load_json(bpopt_param_bounds_path)
        optim_param = self.optim_param
        
        optim_param_bpopt_format = {}
        for key,val in optim_param.items():
            key_param,key_sect = key.split('.')
            try:
                key_sect = bpopt_section_map[key_sect]
            except:
                print('Already in bluepyopt format')
            optim_param_bpopt_format[key_param+'.'+key_sect]=val
            
        param_sens_list = list()
        for i,param_dict in enumerate(param_bounds):
            bound = param_dict.get('bounds')
            if bound:
                name_loc = param_dict['param_name'] + '.' + param_dict['sectionlist']
                lb = min(param_dict['bounds'][0], optim_param_bpopt_format[name_loc]-\
                         max_bound*abs(optim_param_bpopt_format[name_loc]))
                ub = max(param_dict['bounds'][1], optim_param_bpopt_format[name_loc]+\
                         max_bound*abs(optim_param_bpopt_format[name_loc]))
                param_dict['bounds'] = [lb,ub]
            param_sens_list.append(param_dict)
        
        utility.save_json(sens_param_bounds_path,param_sens_list)
        return optim_param_bpopt_format
    
    @staticmethod
    def plot_sobol_analysis(cell_data,features_to_run,
                analysis_path = 'figures/sensitivity_analysis_all_active.pdf'):
    
        nr_plots = len(features_to_run)
        bar_width = 0.75
        opacity = 1
        x_labels = cell_data.uncertain_parameters
        x = np.arange(len(x_labels))  
        
        n_cols = 3
        n_rows = int(math.ceil(nr_plots/float(n_cols)))
        
        utility.create_filepath(analysis_path)
        sns.set()
        fig, axes = plt.subplots(n_rows,n_cols,figsize=(12,8),dpi = 80, 
                             sharex= 'all',sharey = 'row',squeeze=False)
        
        # Turn off blank axis
        fig_empty_index = range(nr_plots,n_rows*n_cols)
        if len(fig_empty_index) != 0:
            for ind in fig_empty_index:
                axes[ind//n_cols,ind%n_cols].axis('off')
        
        current_palette = sns.color_palette()
        colors = sns.color_palette(current_palette,len(x))
        
        for i in range(nr_plots):
            feature = features_to_run[i]
            if feature == 'nrnsim_bpopt': # Voltage feature
                continue
           
            axes[i//n_cols,i%n_cols].bar(x, cell_data[feature].sobol_first_average, 
                bar_width,align='center',alpha=opacity,color=colors)
            
            axes[i//n_cols,i%n_cols].set_xticks(x)
            axes[i//n_cols,i%n_cols].set_xticklabels(x_labels, rotation=45,
                        fontsize=8)
            axes[i//n_cols,i%n_cols].set_title(feature,fontsize=9)
        fig.suptitle('Average First Order Sobol Indices',fontsize=12)        
        fig.tight_layout(rect=[0, 0.03, .95, 0.95])
        fig.savefig(analysis_path, bbox_inches = 'tight')
        plt.close(fig)
        
        
    @staticmethod
    def save_analysis_data(ucdata_path,**kwargs):
    
        uc_data =  un.Data(ucdata_path)
        model_name = uc_data.model_name
        features = list(uc_data.keys())
        features.remove(model_name)
        sens_datalist = []
        
        uc_params = uc_data.uncertain_parameters
        for i,param in enumerate(uc_params):
            sens_datadict = {}
            for feature in features:
                sens_datadict[feature] = uc_data[feature].sobol_first_average[i]
                sens_datadict['param_name'] = param
            for key,val in kwargs.items():
                sens_datadict.update({key:val})
                
            sens_datalist.append(sens_datadict)
        
        sens_datadf = pd.DataFrame(sens_datalist)
        
        return sens_datadf
    
    
    