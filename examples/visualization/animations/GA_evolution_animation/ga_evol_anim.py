import numpy as np
from ateamopt.utils import utility
import seaborn as sns
import matplotlib.pyplot as plt
from ateamopt.animation import animation_module
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
import glob,os
import pandas as pd
from ipyparallel import Client

model_path = os.getcwd()

checkpoint_path = os.path.join(model_path,'checkpoints_final/seed4.pkl')
checkpoint = utility.load_pickle(checkpoint_path)
hof_inds = checkpoint['halloffame']
config_dir= os.path.join(model_path,'config/483101699')
protocol_path = os.path.join(config_dir,'protocols.json')
feature_path = os.path.join(config_dir,'features.json')
parameter_path = os.path.join(config_dir,'parameters.json')
mechanism_path = os.path.join(config_dir,'mechanism.json')
morph_path = os.path.join(model_path,'reconstruction.swc')
eval_handler = Bpopt_Evaluator(protocol_path,feature_path,morph_path,
                               parameter_path,mechanism_path,timeout=60)
evaluator = eval_handler.create_evaluator()
param_names = evaluator.param_names

offspring_size = 512                        
max_ngen = (len(checkpoint['history'].genealogy_history.values())-
               offspring_size)/(2*offspring_size)+1

max_ngen = int(max_ngen)
all_pop = list(checkpoint['history'].genealogy_history.values())

prefix= 'frames_ga_evol/evol_'
select_stim = 'LongDC_55'

gen_vector = []
error_vector=[]
max_fitness_vector = []

rc = Client(profile=os.getenv('IPYTHON_PROFILE'))
lview = rc.load_balanced_view()
nselect_inds = 128

def objective_evol(ax,gen_vector,error_vector,max_fitness_vector,gen,max_ngen):
    ax.scatter(gen_vector,error_vector,s=5,color='salmon',alpha=0.01,
               linewidth=0)
    ax.plot(range(1,gen+2),max_fitness_vector,lw=1,color='b')
    ax.grid(False)
    ax.set_xlim([1,max_ngen])
    ax.set_ylim([0,5000])
    ax.set_xlabel('Generations')
    ax.set_ylabel('Error (sum of z-scores)')

    h,l = ax.get_legend_handles_labels()
    dummy_h = ax.scatter([],[],color='salmon',alpha=1,s=20,label='All Models')
    h.append(dummy_h)
    l.append('All Models')
    ax.legend(h,l,loc='upper center', bbox_to_anchor=(0.5, -.15),frameon=False,
              ncol=2)
    sns.despine(ax=ax)    
    ax.set_title('Evolution')
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    return ax

def run_simulation(sim_tuple):
    print('Simulating Model')
    evaluator,param_dict = sim_tuple
    fitness_protocol = list(evaluator.fitness_protocols.values())[0]
    ind_response = evaluator.run_protocol(protocol=fitness_protocol, 
                                            param_values=param_dict)
    return ind_response


def plot_responses(ax,ind_responses,exp_data,gen_error,error_vector,
                   error_span):
    hof_num = 1
    for jj,ind_response in enumerate(ind_responses):
        if bool(ind_response):
            ind_time = ind_response['%s.soma.v'%select_stim]['time']
            ind_voltage = ind_response['%s.soma.v'%select_stim]['voltage']
            alpha_fitness = np.exp(-(gen_error[jj]- min(error_vector))/error_span)
            # the last indexed individual has highest fitness
            lw_fitness = 1 if jj >= nselect_inds-hof_num else .1   
            color_fitness = 'b' if jj >= nselect_inds-hof_num else 'lightsteelblue'
            if jj == nselect_inds-1:
                ax.plot(ind_time,ind_voltage,color=color_fitness,alpha=1,
                        lw=lw_fitness,label='Best Model')
            else:
                ax.plot(ind_time,ind_voltage,color=color_fitness,alpha=alpha_fitness,
                        lw=lw_fitness)
    
    ax.plot(exp_data[:,0],exp_data[:,1],color='k',lw=1,label='Experiment')
            
    ax.set_xlim([200,1350])
    ax.set_ylim([-105,50])
    ax.grid(False)
    sns.despine(ax=ax)
    ax.set_xlabel('time $(ms)$')
    ax.set_ylabel('voltage $(mV)$',labelpad=-5)
    ax.set_title('Fitness')
    h,l = ax.get_legend_handles_labels()
    if len(h) > 2:
        h,l = h[-2:],l[-2:]
    ax.legend(h,l,loc='upper center', bbox_to_anchor=(0.5, -.15),frameon=False,
              ncol=2)
    return ax

def plot_param_selection(ax,param_df,gen_error):
    hof_num = 1
    cond_params= list(param_df)
    param_df['fitness'] = [1/err_ for err_ in gen_error]
    param_df = param_df.melt(id_vars=['fitness'],value_vars=cond_params,
                 var_name='param_name',value_name='value')
    param_df = param_df.sort_values(by=['fitness','param_name'],axis=0)
    param_df = param_df.reset_index(drop=True)
    # the last indexed individual has highest fitness
    dummy_fitness = [1 if idx >= (offspring_size-hof_num)*len(cond_params) 
                else .8 for idx in param_df.index]
    param_df['dummy_fitness'] = dummy_fitness
    ax = sns.scatterplot(x='param_name',y='value',size='dummy_fitness',
                 data=param_df,hue='dummy_fitness',palette='Blues',
                 linewidth=.05,ax=ax,alpha=.7)
    ax.set(yscale='log')
    plt.setp(ax.get_xticklabels(),rotation=60,ha='right',fontsize=8)
    ax.xaxis.set_major_locator(plt.MaxNLocator(20))
    ax.set_xlabel('')
    ax.set_ylim([1e-6,1e3])
    ax.grid(False)
    sns.despine(ax=ax)   
    ax.get_legend().remove()
    ax.set_title('Selection')
    return ax

def plot_num_models(ax,gen_num):
    bar_height = [gen_num*offspring_size]
    ax.bar(range(len(bar_height)),bar_height,color='teal')
    ax.set_ylim([0,1e5])
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    sns.despine(ax=ax)
    ax.grid(False)
    ax.yaxis.set_major_locator(plt.MaxNLocator(1))
    ax.set_ylabel('# Models Evaluated')
    ax.set_xticklabels([])              
    return ax
    
hof_error = [ind.fitness.sum for ind in hof_inds]
hof_gen_pop = []
for gen in range(max_ngen):
    if gen == 0:
        gen_offspring_idx = 0
    else:
        gen_offspring_idx = offspring_size+(2*gen-1)*offspring_size
    gen_pop = all_pop[gen_offspring_idx:gen_offspring_idx+offspring_size]

    for hof_ind_ in hof_inds:
        if hof_ind_ in gen_pop:
            hof_gen_pop.append(hof_ind_)
            gen_pop.remove(hof_ind_)
    gen_pop += hof_gen_pop
    
    gen_error = [pop_.fitness.sum for pop_ in gen_pop]
    
    # Sort the individuals in descending order of error (ascending order of fitness)
    gen_pop = [ind for _,ind in sorted(zip(gen_error,gen_pop),reverse=True)][-offspring_size:]
    gen_error = sorted(gen_error,reverse=True)[-offspring_size:]
    error_vector.extend(gen_error)
    gen_vector.extend([gen+1]*offspring_size)
    max_fitness_vector.append(min(gen_error))
    error_span = max(error_vector)- min(error_vector)
    
    if (gen+1)%3 == 0:
        filename_ = '%s%03d.jpeg'%(prefix,gen)
        utility.create_filepath(filename_) 
        
        param_dict_list=[{param_name:param_val for param_name,param_val in 
                      zip(param_names,ind_params)} for ind_params in gen_pop]
        param_df = pd.DataFrame(param_dict_list)
        exp_data = np.loadtxt('preprocessed/%s.txt'%select_stim)
        
        sim_tuple_list = [(evaluator,param_dict) for param_dict in param_dict_list]
        ind_responses = lview.map_sync(run_simulation, sim_tuple_list[-nselect_inds:])
        
        sns.set(style='whitegrid')
        fig,ax = plt.subplots(1,4,figsize=(12,4),gridspec_kw={'width_ratios': [2.5,3,2.5,.3]})
        ax[0] = plot_param_selection(ax[0],param_df,gen_error)
        
        ax[1] = plot_responses(ax[1],ind_responses,exp_data,gen_error,error_vector,
                   error_span)
        ax[2] = objective_evol(ax[2],gen_vector,error_vector,max_fitness_vector,
                                      gen,max_ngen)
        
        ax[3] = plot_num_models(ax[3],gen+1)
        fig.subplots_adjust(wspace=.4)
        fig.savefig(filename_,dpi=400,bbox_inches='tight')
        plt.close(fig)
        
snaps = glob.glob('frames_ga_evol/*.jpeg')
snaps = sorted(list(snaps))
anim = animation_module.Animation(movie_name = 'movie.gif')
anim.make_gif(snaps,delay=50,repeat=False)