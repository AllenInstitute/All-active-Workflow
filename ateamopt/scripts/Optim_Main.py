
import bluepyopt as bpopt
import logging
import os
import textwrap
from datetime import datetime
from ateamopt.utils import utility
import shutil
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
from ateamopt.optim_schema import Optim_Config
import argschema as ags



logger = logging.getLogger()


def create_optimizer(args):
    '''returns configured bluepyopt.optimisations.DEAPOptimisation'''

    stage_jobconfig = args['stage_jobconfig']
    highlevel_job_props = args['highlevel_jobconfig']
    seed = args['seed']

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
    
    morph_path = highlevel_job_props['swc_path']
    axon_type = highlevel_job_props['axon_type']
    ephys_dir = highlevel_job_props['ephys_dir']
    
    protocol_path = args['train_protocols']
    mech_path = args['mechanism']
    feature_path = args['train_features']
    param_path = args['parameters']
    
    timeout = stage_jobconfig['timeout']
    learn_eval_trend = stage_jobconfig['learn_eval_trend']
    
    
    eval_handler = Bpopt_Evaluator(protocol_path, feature_path, morph_path,
                                    param_path, mech_path, timeout = timeout,
                                    learn_eval_trend = learn_eval_trend,
                                    axon_type=axon_type,ephys_dir=ephys_dir)
    evaluator = eval_handler.create_evaluator()

    opt = bpopt.optimisations.DEAPOptimisation(
            evaluator=evaluator,
            map_function=map_function,
            seed=seed)
    return opt


def main():
    """Main"""
    mod = ags.ArgSchemaParser(schema_type=Optim_Config)
    args = mod.args
    seed = args['seed']
    stage_jobconfig = args['stage_jobconfig']
    highlevel_job_props = args['highlevel_jobconfig']
    
    logging.basicConfig(level=highlevel_job_props['log_level'])
    
    opt = create_optimizer(args)
    
    
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
