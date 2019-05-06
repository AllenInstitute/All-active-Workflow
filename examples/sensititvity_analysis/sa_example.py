import os
import bluepyopt.ephys as ephys
import bluepyopt as bpopt
import copy
from ateamopt.analysis.sensitivity_analysis import SA_helper
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
from ateamopt.utils import utility
import uncertainpy as un

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
    
    # config files with all the paths for Bluepyopt sim
    model_base_path='/allen/aibs/mat/anin/hpc_trials/' \
                         'BluePyOpt_speedup/483101699/Stage2/'
                         
    ephys_dir = os.path.join(model_base_path, 'preprocessed')
    opt_config_file = '/allen/aibs/mat/anin/hpc_trials/' \
                         'BluePyOpt_speedup/483101699/Stage2/config_file.json'
    
    # optimized parameters around which select parameters are varied
    optim_param_path = '/allen/aibs/mat/anin/hpc_trials/BluePyOpt_speedup/'\
    '483101699/Stage2/fitted_params/optim_param_483101699_bpopt.json'
    
    # Parameters to vary
    
    sens_params_path = 'select_params.json' # knobs
    if os.path.exists(sens_params_path):
        select_parameters = utility.load_json(sens_params_path)
    else:
        select_parameters = {'gbar_NaTs2_t':['somatic'],
                             'gbar_NaTa_t' : ['axonal'],
                             'gbar_Nap_Et2' : ['axonal','somatic'],
                             'gbar_Kv3_1' : ['axonal','somatic'],
                             'gbar_K_Tst' :  ['axonal','somatic']
                             }
        utility.save_json(sens_params_path,select_parameters)                         
    
    # Features 
    efel_features_path = 'features_sobol.json' # knobs
    
    if os.path.exists(efel_features_path):
        efel_features = utility.load_json(efel_features_path)
    else:    
        efel_features = ["AP_amplitude","AP_width",
                         "mean_frequency", "adaptation_index2", 
                         "AHP_depth", "time_to_first_spike","Spikecount"]
        utility.save_json(efel_features_path,efel_features)    
    
    param_mod_range = .1 # knobs
    SA_obj = SA_helper(optim_param_path,sens_params_path,param_mod_range,
                       opt_config_file)
    morph_path,protocol_path,mech_path,feature_path,\
        param_bound_path = SA_obj.load_config(model_base_path)
    
    # Make sure to get the parameter bounds big enough for BluePyOpt sim
    sens_param_bound_write_path = "parameters_sensitivity.json"

    optim_param = SA_obj.create_sa_bound(param_bound_path,sens_param_bound_write_path)
    param_dict_uc = SA_obj.create_sens_param_dict()
    parameters ={key:optim_param[val] for key,val in param_dict_uc.items()}
    
    eval_handler = Bpopt_Evaluator(protocol_path, feature_path,
                                   morph_path, param_bound_path,
                                   mech_path,
                                   ephys_dir=ephys_dir,
                                   timed_evaluation = False)
    evaluator = eval_handler.create_evaluator()
    opt = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator)
    stim_protocols = utility.load_json(protocol_path)
    
    stim_name= 'LongDC_55' # knobs
    
    # Check for compiled modfiles
    if not os.path.isdir('x86_64'):
        raise Exception('Compiled modfiles do not exist')
    
    un_features = un.EfelFeatures(features_to_run=efel_features)
    un_parameters = un.Parameters(parameters)
    un_parameters.set_all_distributions(un.uniform(param_mod_range))
    
    un_model = un.Model(run=nrnsim_bpopt, interpolate=True,
                 labels=["Time (ms)", "Membrane potential (mV)"],
                 opt=opt,stim_protocols =stim_protocols,
                 param_dict_uc = param_dict_uc,
                 stim_name=stim_name,
                 optim_param=optim_param)
    
    # Perform the uncertainty quantification
    UQ = un.UncertaintyQuantification(un_model,
                                      parameters=un_parameters,
                                      features=un_features)
    UQ.quantify(seed=0,CPUs=32)
    cell_data =  un.Data("data/nrnsim_bpopt.h5")
    SA_obj.plot_sobol_analysis(cell_data,efel_features)
    
if __name__ == '__main__':
    main()
