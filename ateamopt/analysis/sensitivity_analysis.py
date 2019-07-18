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
from matplotlib.cm import ScalarMappable
from matplotlib.colors import ListedColormap


class SA_helper(object):
    
    def __init__(self, optim_param_path, sens_param_path,
                 param_mod_range,config_file):
        self.config = utility.load_json(config_file) if \
            sens_param_path else None
        self.sens_parameters = utility.load_json(sens_param_path) if \
            sens_param_path else None
        if optim_param_path:
            self.optim_param = utility.load_json(optim_param_path)
        else:
            self.optim_param = None
        self.param_range = param_mod_range if param_mod_range else 0
    
    def load_config(self,model_basepath=None,perisomatic=False):
        morph_path = os.path.normpath(self.config['morphology'])
        protocol_path = os.path.normpath(self.config['all_protocols'])
        feature_path = os.path.normpath(self.config['features'])
        mech_path = os.path.normpath(self.config['mechanism'])
        param_path = os.path.normpath(self.config['parameters'])
        
        if perisomatic:
            mech_path = os.path.normpath(self.config['peri_mechanism'])
            param_path = os.path.normpath(self.config['peri_parameters'])
            
        
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
        
        if 'genome' in optim_param.keys():
            print('The parameter file is in AIBS format')
            for aibs_param_dict in optim_param['genome']:
                param_name,param_sect = aibs_param_dict['name'],\
                                            aibs_param_dict['section']
                if param_name in ['e_pas','g_pas','Ra']:
                    param_sect = 'all'
                                                
                param_sect = bpopt_section_map[param_sect]
                optim_param_bpopt_format[param_name+'.'+param_sect]=\
                                float(aibs_param_dict['value'])
        else:
            print('The parameter file is in bluepyopt format')
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
    
    
    def create_sa_bound_peri(self,bpopt_param_bounds_path,sens_param_bounds_path,
                        max_bound = .5):
        
        # For parameter sensitivity create a new set of bounds because the permutations may 
        # fall outside the original bounds
        
        max_bound = max(max_bound,self.param_range+.1)        
        param_bounds = utility.load_json(bpopt_param_bounds_path)
        
        optim_param_bpopt_format = {}
        
        param_sens_list = list()
        for i,param_dict in enumerate(param_bounds):
            
            if 'sectionlist' in param_dict.keys():
                name_loc = param_dict['param_name'] + '.' + \
                                param_dict['sectionlist']
                
                if param_dict['param_name'] not in ['ena','ek']:
                    optim_param_bpopt_format[name_loc] = param_dict['value']
                    lb =  param_dict['value']-\
                         max_bound*abs(param_dict['value'])
                    ub = param_dict['value']+\
                         max_bound*abs(param_dict['value'])
                    param_dict['bounds'] = [lb,ub]
                    del param_dict['value']
            param_sens_list.append(param_dict)
        
        utility.save_json(sens_param_bounds_path,param_sens_list)
        return optim_param_bpopt_format
    
    @staticmethod
    def plot_sobol_analysis(cell_data,
                analysis_path = 'figures/sensitivity_analysis_all_active.pdf',
                palette='Set1'):
    
        features_to_run = list(cell_data.keys())
        features_to_run.remove(cell_data.model_name) # Voltage feature
        
        nr_plots = len(features_to_run)
        bar_width = 0.75
        opacity = 1
        x_labels = cell_data.uncertain_parameters
        x = np.arange(len(x_labels))  
        
        n_cols = 8
        n_rows = int(math.ceil(nr_plots/float(n_cols)))
        
        if n_rows == 1:
            n_cols = nr_plots
         
        utility.create_filepath(analysis_path)
        sns.set(style='whitegrid')
        fig, axes = plt.subplots(n_rows,n_cols,dpi = 80,figsize=(13,2.5), 
                             sharex= 'all',sharey = 'row',squeeze=False)
        # Turn off blank axis
        fig_empty_index = range(nr_plots,n_rows*n_cols)
        if len(fig_empty_index) != 0:
            for ind in fig_empty_index:
                axes[ind//n_cols,ind%n_cols].axis('off')
        
        current_palette = sns.color_palette(palette)
        colors = sns.color_palette(current_palette,len(x))
        my_cmap = ListedColormap(colors)
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0,len(x)-1))
        sm.set_array([])
        
        for i in range(nr_plots):
            feature = features_to_run[i]
            
            axes[i//n_cols,i%n_cols].bar(x, cell_data[feature].sobol_first_average, 
                bar_width,align='center',alpha=opacity,color=colors)
            
            axes[i//n_cols,i%n_cols].set_xticks(x)
            axes[i//n_cols,i%n_cols].set_xticklabels(x_labels, rotation=60,
                        ha='right',fontsize=8)
            axes[i//n_cols,i%n_cols].set_title(feature,fontsize=9)
            axes[i//n_cols,i%n_cols].grid(axis='x')
            if i%n_cols == 0:
                axes[i//n_cols,i%n_cols].set_ylabel('Sobol First Average',
                    fontsize=9)
                plt.setp(axes[i//n_cols,i%n_cols].get_yticklabels(),fontsize=8)

#        fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        
        cbar = plt.colorbar(sm,boundaries=np.arange(len(x)+1)-0.5)
        cbar.set_ticks(np.arange(len(x)))
        cbar.ax.set_yticklabels(x_labels, fontsize=8)
        fig.savefig(analysis_path, bbox_inches = 'tight')
        plt.close(fig)
        
    

    @staticmethod
    def plot_sobol_analysis_from_df(sa_data_df,
                analysis_path = 'figures/sens_analysis.pdf',palette='Set1'):
        
        param_names = sorted(list(sa_data_df.param_name.unique()))

        
        current_palette = sns.color_palette(palette)
        colors = sns.color_palette(current_palette,len(param_names))
        my_cmap = ListedColormap(colors)
        sm = ScalarMappable(cmap=my_cmap, norm=plt.Normalize(0,
                                     len(param_names)-1))
        sm.set_array([])
        
        g = sns.FacetGrid(sa_data_df,col='feature',
              sharex=True,sharey='row', height=6, 
              aspect=.5)
        g = g.map(sns.barplot,'param_name','sobol_index',
                  order=param_names,palette=palette,errwidth=1)        
        axes = g.axes.flatten()
        for ax_ in axes:
            title_ = ax_.get_title()
            title_ = title_.split('|')[-1].split('=')[-1]
            ax_.set_title(title_,fontsize=14)
            xticklabels = ax_.get_xticklabels()
            ax_.set_xticklabels(xticklabels,rotation=90,ha='center',fontsize=12)
            ax_.set_xlabel(None)
        axes[0].set_ylabel('sobol index',fontsize=14)
        g.fig.tight_layout(rect=[0, 0.03, 1, 0.95])
        cbar = plt.colorbar(sm,boundaries=np.arange(len(param_names)+1)-0.5)
        cbar.set_ticks(np.arange(len(param_names)))
        cbar.ax.set_yticklabels(param_names, fontsize=12)    
        g.fig.savefig(analysis_path,bbox_inches='tight')
        plt.close(g.fig)
        

    
    @staticmethod
    def save_analysis_data(ucdata_path,**kwargs):
    
        uc_data =  un.Data(ucdata_path)
        model_name = uc_data.model_name
        features = list(uc_data.keys())
        features.remove(model_name)
        sens_datalist = []
        
        filepath = kwargs.pop('filepath',None)
        
        uc_params = uc_data.uncertain_parameters
        for i,param in enumerate(uc_params):
            for feature in features:
                sens_datadict = {}
                sens_datadict['feature'] = feature
                sens_datadict['sobol_index'] = uc_data[feature].sobol_first_average[i]
                sens_datadict['param_name'] = param
                for key,val in kwargs.items():
                    sens_datadict.update({key:val})
                
                sens_datalist.append(sens_datadict.copy())
        
        sens_datadf = pd.DataFrame(sens_datalist)
        
        if filepath:
            utility.create_filepath(filepath)
            sens_datadf.to_csv(filepath, index=False)
        
        return sens_datadf
    
    
    