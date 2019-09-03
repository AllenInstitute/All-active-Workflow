
import bluepyopt as bpopt
import logging
import os
import argparse
from datetime import datetime
from ateamopt.utils import utility
import shutil
from ateamopt.bpopt_evaluator import Bpopt_Evaluator


logger = logging.getLogger()


def create_optimizer(job_config,seed):
    '''returns configured bluepyopt.optimisations.DEAPOptimisation'''

    stage_jobconfig = job_config['stage_jobconfig']
    highlevel_job_props = job_config['highlevel_jobconfig']
    

    if stage_jobconfig['ipyp_optim']:
        from ipyparallel import Client
        rc = Client(profile=os.getenv('IPYTHON_PROFILE'))

        logger.debug('Using ipyparallel with %d engines', len(rc))
        lview = rc.load_balanced_view()

        def mapper(func, it):
            start_time = datetime.now()
            ret = lview.map_sync(func, it)
            logger.debug('Generation took %s', datetime.now() - start_time)

            # Save timing information for each generation
            f =  open('time_info.txt','a')
            f.write('%s\n'%(datetime.now() - start_time))
            f.close()

            return ret

        map_function = mapper
    else:
        map_function = None

    seed = os.getenv('BLUEPYOPT_SEED', seed)

    # load the configuration paths    
    
    ephys_dir = highlevel_job_props['ephys_dir']
    morph_path = highlevel_job_props['swc_path']
    protocol_path = job_config['train_protocols']
    mech_path = job_config['mechanism']
    feature_path = job_config['train_features']
    param_path = job_config['parameters']
    
    axon_type = highlevel_job_props.get('axon_type')
    
    props = {}
    for prop in ['timeout','learn_eval_trend']:
        if stage_jobconfig.get(prop):
            props[prop] = stage_jobconfig.get(prop)
    
    
    eval_handler = Bpopt_Evaluator(protocol_path, feature_path, morph_path,
                                    param_path, mech_path, axon_type=axon_type,
                                    ephys_dir=ephys_dir,**props)
    evaluator = eval_handler.create_evaluator()

    opt = bpopt.optimisations.DEAPOptimisation(
            evaluator=evaluator,
            map_function=map_function,
            seed=seed)
    return opt

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--job_config', required=True, 
                        help='Job config path')
    parser.add_argument('--seed', type=int, default=1,
                        help='Seed to use for optimization')
    return parser

def main():
    """Main"""
    args = get_parser().parse_args()
    
    seed = args.seed
    job_config_path = args.job_config
    job_config = utility.load_json(job_config_path)
    
    stage_jobconfig = job_config['stage_jobconfig']
    highlevel_job_props = job_config['highlevel_jobconfig']
    
    logging.basicConfig(level=highlevel_job_props['log_level'])
    
    opt = create_optimizer(job_config,seed)
    
    
    cp_file = os.path.join(stage_jobconfig['cp_dir'],'seed%s.pkl'%seed)
    utility.create_filepath(cp_file)
    if stage_jobconfig.get('cp_backup_dir'):
        cp_backup_file = os.path.join(stage_jobconfig['cp_backup_dir'],'seed%s.pkl'%seed)
        utility.create_filepath(cp_backup_file)
    else:
        cp_backup_file = None
    cp_backup_frequency = stage_jobconfig['cp_backup_frequency']
    max_ngen = stage_jobconfig['max_ngen']
    offspring_size = stage_jobconfig['offspring_size']
    
    continue_cp = os.path.exists(cp_file)
    logger.debug('Doing start or continue')
    
    if os.path.exists(cp_file):
        try:
            _ = utility.load_pickle(cp_file)
        except:
            logger.debug('Checkpoint file is corrupt! Looking for backup')
            if cp_backup_file and os.path.exists(cp_backup_file):
                shutil.copyfile(cp_backup_file,cp_file)
        
    opt.run(max_ngen=max_ngen,
            offspring_size=offspring_size,
            continue_cp=continue_cp,
            cp_filename=cp_file,
            cp_backup=cp_backup_file,
            cp_backup_frequency=cp_backup_frequency)


if __name__ == '__main__':
    main()
