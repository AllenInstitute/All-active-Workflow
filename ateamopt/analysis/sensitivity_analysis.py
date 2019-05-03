import uncertainpy as un
import json
import bluepyopt as bpopt
import numpy as np
import bluepyopt.ephys as ephys
import matplotlib.pyplot as plt
from ateamopt.utils import utility
import copy
import errno
import os
import pandas as pd


class Sensititvity_Analysis(object):
    
    def __init__(self, optim_param_path, sens_param_path, config_file,stim_name):
        self.config = utility.load_json(config_file)
        self.sens_param_path = sens_param_path
        self.stim_name = stim_name
        self.optim_param_path = optim_param_path
    
    
    def load_config(self):
        morph_path = self.config['morphology']
        protocol_path = self.config['protocols']
        mech_path = self.config['mechanism']
        feature_path = self.config['features']
        param_path = self.config['parameters']
        return morph_path,protocol_path,mech_path,feature_path,param_path

    
    def create_sens_param_dict(self):
        param_dict_uc = {}
        for key,val in self.sens_parameters.items():
            for sect in val:
                param_dict_uc['%s_%s'%(key,sect)] = '%s.%s'%(key,sect)
        return param_dict_uc
    
    def create_sa_bound(self,bpopt_param_path, max_bound = .5):
        
        # For parameter sensitivity create a new set of bounds because the permutations may 
        # fall outside the original bounds

        param_bound = utility.load_json(bpopt_param_path)
        optim_param = utility.load_json(self.optim_param_path)
       
        param_sens_list = list()
        for i,param_dict in enumerate(param_bound):
            bound = param_dict.get('bounds')
            if bound:
                name_loc = param_dict['param_name'] + '.' + param_dict['sectionlist']
                lb = min(param_dict['bounds'][0], optim_param[name_loc]-\
                         max_bound*abs(optim_param[name_loc]))
                ub = max(param_dict['bounds'][1], optim_param[name_loc]+\
                         max_bound*abs(optim_param[name_loc]))
                param_dict['bounds'] = [lb,ub]
            param_sens_list.append(param_dict)
        
        utility.save_json(self.sens_param_path,param_sens_list)
            
        return optim_param
    
    @staticmethod
    def plot_sobol_analysis(cell_data,features_to_run,
                analysis_path = 'figures/sensitivity_analysis_all_active.pdf'):
    
        nr_plots = len(features_to_run)
        bar_width = 0.75
        opacity = 1
        x_labels = cell_data.uncertain_parameters
        x = np.arange(len(x_labels))  
        
        utility.create_filepath(analysis_path)
        plt.style.use('ggplot')
        fig, axes = plt.subplots(1,nr_plots,figsize=(12,4),dpi = 80, 
                                 sharey = True)
        
        cmap = plt.get_cmap('seismic')
        colors = cmap(np.linspace(0, 1, len(x)))
        
        for i in range(nr_plots):
            feature = features_to_run[i]
            if feature == 'nrnsim_bpopt': # Voltage feature
                continue
           
            axes[i].bar(x, cell_data[feature].sobol_first_average, 
                bar_width,align='center',alpha=opacity,color=colors)
           
            axes[i].set_xticklabels(x_labels, rotation=45,ha = 'right',
                    fontsize=8)
            axes[i].set_xticks(x)
            axes[i].grid(axis = 'x')
            
            
            axes[i].set_title(feature,fontsize=9)
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
    
    def nrnsim_bpopt(**kwargs):
        
        opt = kwargs.pop('opt')
        fitness_protocols = opt.evaluator.fitness_protocols
        nrn = ephys.simulators.NrnSimulator()
        sensitivity_params  = copy.deepcopy(kwargs.pop('optim_param'))
        for key,val in kwargs.items():
            sensitivity_params[param_dict_uc[key]] = val
            
        sensitivity_response = fitness_protocols[self.stim_name].run(\
                    cell_model=opt.evaluator.cell_model,
                    param_values=sensitivity_params,
                    sim=nrn)
        name_loc = self.stim_name + '.soma.v'
        time = sensitivity_response[name_loc]['time']
        value = sensitivity_response[name_loc]['voltage']
        info = {'stimulus_start':stim_protocols[self.stim_name]['stimuli'][0]['delay'], 
                'stimulus_end':stim_protocols[self.stim_name]['stimuli'][0]['stim_end']}
        
        return time, value, info
    