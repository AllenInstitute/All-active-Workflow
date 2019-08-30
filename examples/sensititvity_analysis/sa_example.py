import os
import bluepyopt.ephys as ephys
import bluepyopt as bpopt
import copy
from ateamopt.analysis.sensitivity_analysis import SA_helper
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
from ateamopt.utils import utility
import uncertainpy as un
import operator
from ateam.data import lims
import multiprocessing as mp
import sys,shutil


def nrnsim_bpopt(**kwargs):
    opt = kwargs.pop('opt')
    stim_name = kwargs.pop('stim_name')
    stim_protocols = kwargs.pop('stim_protocols')
    param_dict_uc = kwargs.pop('param_dict_uc')
    optim_param = kwargs.pop('optim_param')
    sensitivity_params  = copy.deepcopy(optim_param)
    
    fitness_protocols = opt.evaluator.fitness_protocols
    nrn = ephys.simulators.NrnSimulator()
    sensitivity_params  = copy.deepcopy(optim_param)
    
    for key,val in kwargs.items():
        sensitivity_params[param_dict_uc[key]] = val
        
    sensitivity_response = fitness_protocols[stim_name].run(\
                cell_model=opt.evaluator.cell_model,
                param_values=sensitivity_params,
                sim=nrn)
    name_loc = stim_name + '.soma.v'
    time = sensitivity_response[name_loc]['time']
    value = sensitivity_response[name_loc]['voltage']
    info = {'stimulus_start':stim_protocols[stim_name]['stimuli'][0]['delay'], 
            'stimulus_end':stim_protocols[stim_name]['stimuli'][0]['stim_end']}
    
    return time, value, info

def main():
    
    # Read sensitivity analysis config file
    sens_config_file = sys.argv[-1]
    sens_config_dict = utility.load_json(sens_config_file)
    cell_id = sens_config_dict['Cell_id']
    cpu_count = sens_config_dict['cpu_count'] if 'cpu_count'\
            in sens_config_dict.keys() else mp.cpu_count()
    perisomatic_sa = sens_config_dict.get('run_peri_analysis',False)
    
    # Parameters to vary (All-active) 
    select_aa_param_path = sens_config_dict['select_aa_param_path'] # knobs
    
    # Parameters to vary (Perisomatic) 
    if perisomatic_sa:
        select_peri_param_path = sens_config_dict['select_peri_param_path'] # knobs
    
    select_feature_path = sens_config_dict['select_feature_path'] # knobs
    param_mod_range = sens_config_dict.get('param_mod_range',.1) # knobs
    mechanism_path = sens_config_dict['mechanism']
    
    # config files with all the paths for Bluepyopt sim    
    lr = lims.LimsReader()
    morph_path = lr.get_swc_path_from_lims(int(cell_id))
    
    model_base_path='/allen/aibs/mat/ateam_shared/' \
                         'Mouse_Model_Fit_Metrics/{}'.format(cell_id)
                         
    opt_config_file = os.path.join(model_base_path,'config_file.json')
    if not os.path.exists(opt_config_file):
        opt_config = {
                "morphology": "",
                "parameters": "config/{}/parameters.json".format(cell_id),
                "mechanism": "config/{}/mechanism.json".format(cell_id),
                "protocols": "config/{}/protocols.json".format(cell_id),
                "all_protocols": "config/{}/all_protocols.json".format(cell_id),
                "features": "config/{}/features.json".format(cell_id),
                "peri_parameters": "config/{}/peri_parameters.json".format(cell_id),
                "peri_mechanism": "config/{}/peri_mechanism.json".format(cell_id)
                }
        opt_config_file = os.path.join(os.getcwd(),'config_file.json')
        utility.save_json(opt_config_file,opt_config)
    
    # optimized parameters around which select parameters are varied
    optim_param_path_aa = '/allen/aibs/mat/ateam_shared/Mouse_Model_Fit_Metrics/'\
    '{cell_id}/fitted_params/optim_param_unformatted_{cell_id}.json'.\
                    format(cell_id = cell_id)
    if not os.path.exists(optim_param_path_aa):
        optim_param_path_aa = '/allen/aibs/mat/ateam_shared/Mouse_Model_Fit_Metrics/'\
            '{cell_id}/fitted_params/optim_param_{cell_id}_bpopt.json'.\
                    format(cell_id = cell_id)
    
    SA_obj_aa = SA_helper(optim_param_path_aa,select_aa_param_path,param_mod_range,
                       opt_config_file)
    
    _,protocol_path,mech_path,feature_path,\
        param_bound_path = SA_obj_aa.load_config(model_base_path)
        
    # Make sure to get the parameter bounds big enough for BluePyOpt sim
    sens_param_bound_write_path_aa = "param_sensitivity_aa.json"
    optim_param_aa = SA_obj_aa.create_sa_bound(param_bound_path,
                                         sens_param_bound_write_path_aa)    
    param_dict_uc_aa = SA_obj_aa.create_sens_param_dict()
    parameters_aa ={key:optim_param_aa[val] for key,val in param_dict_uc_aa.items()}
    eval_handler_aa = Bpopt_Evaluator(protocol_path, feature_path,
                                   morph_path, sens_param_bound_write_path_aa,
                                   mech_path,
                                   ephys_dir=None,
                                   timed_evaluation = False)
    evaluator_aa = eval_handler_aa.create_evaluator()
    opt_aa = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator_aa)
    
    
    stim_protocols = utility.load_json(protocol_path)
    stim_protocols = {key:val for key,val in stim_protocols.items() \
                      if 'LongDC' in key}
    stim_dict = {key:val['stimuli'][0]['amp'] \
                     for key,val in stim_protocols.items()}
    sorted_stim_tuple= sorted(stim_dict.items(), key=operator.itemgetter(1))
    
    stim_name= sorted_stim_tuple[-1][0] # knobs (the max amp)
    
    # Copy compiled modfiles
    if not os.path.isdir('x86_64'):
        raise Exception('Compiled modfiles do not exist')
    
    efel_features = utility.load_json(select_feature_path)
    un_features = un.EfelFeatures(features_to_run=efel_features)
    
    un_parameters_aa = un.Parameters(parameters_aa)
    un_parameters_aa.set_all_distributions(un.uniform(param_mod_range))
    un_model_aa = un.Model(run=nrnsim_bpopt, interpolate=True,
                 labels=["Time (ms)", "Membrane potential (mV)"],
                 opt=opt_aa,stim_protocols =stim_protocols,
                 param_dict_uc = param_dict_uc_aa,
                 stim_name=stim_name,
                 optim_param=optim_param_aa)
    
    
    # Perform the uncertainty quantification
    UQ_aa = un.UncertaintyQuantification(un_model_aa,
                                      parameters=un_parameters_aa,
                                      features=un_features)    
    data_folder = 'sensitivity_data'
    sa_filename_aa = 'sa_allactive_%s.h5'%cell_id
    sa_filename_aa_csv = 'sa_allactive_%s.csv'%cell_id
    sa_data_path_aa = os.path.join(data_folder,sa_filename_aa)
    sa_aa_csv_path = os.path.join(data_folder,sa_filename_aa_csv)
    
    UQ_aa.quantify(seed=0,CPUs=cpu_count,data_folder=data_folder,
                   filename= sa_filename_aa)
    _ = SA_obj_aa.save_analysis_data(sa_data_path_aa,
                                filepath=sa_aa_csv_path)
        
    cell_data_aa =  un.Data(sa_data_path_aa)
    SA_obj_aa.plot_sobol_analysis(cell_data_aa,analysis_path = \
                          'figures/sa_analysis_aa_%s.pdf'%cell_id,
                          palette='Set1')
    
    # Perisomatic model
    
    if perisomatic_sa:
    
        try:
            optim_param_path_peri = None
            SA_obj_peri = SA_helper(optim_param_path_peri,select_peri_param_path,param_mod_range,
                                   opt_config_file)
            _,_,mech_path_peri,_,\
                    param_bound_path_peri = SA_obj_peri.load_config(model_base_path,
                                                                perisomatic=True)
            
            sens_param_bound_write_path_peri = "param_sensitivity_peri.json"
            optim_param_peri = SA_obj_peri.create_sa_bound_peri(param_bound_path_peri,
                                                     sens_param_bound_write_path_peri)
            
            param_dict_uc_peri = SA_obj_peri.create_sens_param_dict()
            parameters_peri ={key:optim_param_peri[val] for key,val in param_dict_uc_peri.items()}
            eval_handler_peri = Bpopt_Evaluator(protocol_path, feature_path,
                                               morph_path, sens_param_bound_write_path_peri,
                                               mech_path_peri,
                                               ephys_dir=None,
                                               timed_evaluation = False)
            evaluator_peri = eval_handler_peri.create_evaluator()
            opt_peri = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator_peri)
            un_parameters_peri= un.Parameters(parameters_peri)
            un_parameters_peri.set_all_distributions(un.uniform(param_mod_range))
            un_model_peri = un.Model(run=nrnsim_bpopt, interpolate=True,
                             labels=["Time (ms)", "Membrane potential (mV)"],
                             opt=opt_peri,stim_protocols =stim_protocols,
                             param_dict_uc = param_dict_uc_peri,
                             stim_name=stim_name,
                             optim_param=optim_param_peri)
            UQ_peri = un.UncertaintyQuantification(un_model_peri,
                                                  parameters=un_parameters_peri,
                                                  features=un_features)
            sa_filename_peri = 'sa_perisomatic_%s.h5'%cell_id
            sa_filename_peri_csv = 'sa_perisomatic_%s.csv'%cell_id
            sa_data_path_peri = os.path.join(data_folder,sa_filename_peri)
            sa_peri_csv_path = os.path.join(data_folder,sa_filename_peri_csv)
            
            UQ_peri.quantify(seed=0,CPUs=cpu_count,data_folder=data_folder,
                           filename= sa_filename_peri)
            _ = SA_obj_peri.save_analysis_data(sa_data_path_peri,
                                            filepath=sa_peri_csv_path)
            cell_data_peri =  un.Data(sa_data_path_peri)    
            SA_obj_peri.plot_sobol_analysis(cell_data_peri,analysis_path = \
                                      'figures/sa_analysis_peri_%s.pdf'%cell_id,
                                      palette='Set2')
        except Exception as e:
            print(e)
      
        
#    shutil.rmtree('x86_64')
        
if __name__ == '__main__':
    main()
