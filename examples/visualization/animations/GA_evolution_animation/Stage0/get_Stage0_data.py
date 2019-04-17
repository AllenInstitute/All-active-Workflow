import os
import bluepyopt.ephys as ephys
import bluepyopt as bpopt
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
from ateamopt.utils import utility
import numpy as np
from ateamopt.animation import animation_module

cell_id = '483101699'

opt_config_filename = 'config_file.json'
opt_config = utility.load_json(opt_config_filename)

all_protocols_write_path = opt_config['all_protocols']
features_write_path = opt_config['features']
morph_path = opt_config['morphology']
param_write_path = opt_config['parameters']
mech_write_path = opt_config['mechanism']


eval_handler = Bpopt_Evaluator(all_protocols_write_path, features_write_path,
                               morph_path, param_write_path,
                               mech_write_path)
evaluator = eval_handler.create_evaluator()
opt = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator)

checkpoint_file = 'checkpoints/seed1.pkl'
checkpoint = utility.load_pickle(checkpoint_file)

off_spring_size = 512
pop_size_accumulated = len(checkpoint['history'].genealogy_history.values())
gen_size = int((pop_size_accumulated - off_spring_size)/(2*off_spring_size) \
               +1)

min_fitness_individuals = list()
min_fitness_pop = list()


for gen in range(gen_size):
    batch_size = off_spring_size if gen == 0 else 2*off_spring_size
    prev_batch_size = off_spring_size if gen == 1 else  2*off_spring_size
    batch_pop = list(checkpoint['history'].genealogy_history.values())[gen*prev_batch_size:\
                          (gen+1)*batch_size]
    batch_fitness = [batch_pop[i].fitness.sum for i in range(len(batch_pop))]
    min_fitness_batch = [x for _,x in sorted(zip(batch_fitness,batch_pop))]
    min_fitness = min(batch_fitness)
    
    min_fitness_individuals.append(min_fitness_batch[0])
    min_fitness_pop.append(min_fitness)

min_fitness_individual_sorted =[x for _,x in \
            sorted(zip(min_fitness_pop,min_fitness_individuals),reverse=True)]
min_fitness_pop=sorted(min_fitness_pop,reverse=True)

nrn = ephys.simulators.NrnSimulator()
fitness_protocols = opt.evaluator.fitness_protocols
param_names= opt.evaluator.param_names
selected_individuals_responses = []  

select_stim_key ='LongDC_46'
protocol_dict = utility.load_json(all_protocols_write_path)\
                [select_stim_key]['stimuli'][0]
ephys_path = os.path.join(os.getcwd(),'ephys/%s.txt'%select_stim_key)
ephys_data_select_stim = np.loadtxt(ephys_path)

GA_responses_file = 'GA_responses_stage0.pkl'

if os.path.exists(GA_responses_file):
    selected_individuals_responses=utility.load_pickle(GA_responses_file)    
else:
    for i,individual in enumerate(min_fitness_individual_sorted):
        ga_individual = dict()
        for j, param_name in enumerate(param_names):
            ga_individual[param_name] = individual[j]
        
        response_ind = fitness_protocols[select_stim_key].run(
            cell_model=opt.evaluator.cell_model,
            param_values=ga_individual,
            sim=nrn)
        selected_individuals_responses.append(response_ind)
    utility.save_pickle(GA_responses_file,selected_individuals_responses)
    
abs_individual_list = []
for individual in min_fitness_individual_sorted:
    abs_individual_list.append(list(map(abs,individual)))
    
from matplotlib import gridspec
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
    
plt.style.use('ggplot')

anim_prefix = 'GA_anim'
files = []

for ii in range(gen_size):
    fname = '%s%03d.jpeg'%(anim_prefix,ii)
    fig,ax = plt.subplots(1,3,figsize=(12,4),gridspec_kw = \
                                  {'width_ratios':[9, .8, 8]})
#    plt.subplots_adjust(left=0.05, right=0.98)
    ax[0].scatter(np.arange(len(param_names)), abs_individual_list[ii], marker = 'o',
                       alpha = 0.8, s=100, color= 'red')
    ax[0].set_yscale('log')
    ax[0].set_ylim([1e-6, 1e3])
    ax[0].set_xticks(np.arange(len(param_names)))
    ax[0].set_xticklabels(param_names,rotation=45,fontsize =12,
              horizontalalignment='right')
    ax[0].set_ylabel('Absolute Parameter values',fontsize =12)
    
    ax[1].get_xaxis().set_visible(False)
    ax[1].yaxis.set_label_position("left")
    ax[1].set_ylabel('Error',fontsize =12)
    ax[1].grid('on')
    ax[1].set_facecolor('white')
    ax[1].yaxis.set_major_locator(MaxNLocator(2))
    ax[1].set_ylim([0, 200])
#    ax[1].set_yticks([])
    bar_width = 0.2
    opacity = 0.5
    ax[1].bar([0], [min_fitness_pop[ii]], bar_width, alpha=opacity,\
              color='r')
    
    resp_time = selected_individuals_responses[ii][select_stim_key+'.soma.v']['time']
    resp_voltage = selected_individuals_responses[ii][select_stim_key+'.soma.v']['voltage']
    ax[2].plot(resp_time,resp_voltage, color= 'b')
    ax[2].plot(ephys_data_select_stim[:,0],ephys_data_select_stim[:,1],color='k')
    ax[2].set_xlim([protocol_dict['delay']-200, protocol_dict['stim_end']+200])
    ax[2].set_ylabel('Voltage (mV)',fontsize =12)
    ax[2].set_xlabel('Time (ms)',fontsize =12)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    file_path = os.path.join('GA_figures',fname)
    utility.create_filepath(file_path)
    fig.savefig(file_path)
    plt.close(fig)
    files.append(file_path)

anim_handler = animation_module.Animation(movie_name = 'Stage0_GA.gif')
anim_handler.make_gif(files)