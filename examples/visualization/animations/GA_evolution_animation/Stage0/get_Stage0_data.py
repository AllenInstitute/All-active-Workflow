import numpy as np
import bluepyopt.ephys as ephys
import bluepyopt as bpopt
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
from ateamopt.utils import utility


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
opt_train = bpopt.optimisations.DEAPOptimisation(evaluator=evaluator)
checkpoint_file = 'checkpoints/seed1.pkl'
checkpoint = utility.load_pickle(checkpoint_file)

off_spring_size = 512
pop_size_accumulated = len(checkpoint['history'].genealogy_history.values())
gen_size = int((pop_size_accumulated - off_spring_size)/(2*off_spring_size) \
               +1)

min_fitness_inds = list()
min_fitness_pop = list()


for gen in range(gen_size):
    batch_size = off_spring_size if gen == 0 else 2*off_spring_size
    prev_batch_size = off_spring_size if gen == 1 else  2*off_spring_size
    batch_pop = list(checkpoint['history'].genealogy_history.values())[gen*prev_batch_size:\
                          (gen+1)*batch_size]
    batch_fitness = [batch_pop[i].fitness.sum for i in range(len(batch_pop))]
    min_fitness_batch = [x for _,x in sorted(zip(batch_fitness,batch_pop))]
    min_fitness = min(batch_fitness)
    
    min_fitness_inds.append(min_fitness_batch[0])
    min_fitness_pop.append(min_fitness)

frame_num = 60