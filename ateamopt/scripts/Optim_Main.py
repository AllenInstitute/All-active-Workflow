
import bluepyopt as bpopt
import argparse
import logging
import os
import sys
import textwrap
from datetime import datetime
from ateamopt.utils import utility
import shutil
from ateamopt.bpopt_evaluator import Bpopt_Evaluator


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()



def create_optimizer(args):
    '''returns configured bluepyopt.optimisations.DEAPOptimisation'''

    if args.ipyparallel:
        from ipyparallel import Client
        rc = Client(profile=os.getenv('IPYTHON_PROFILE'))

        logger.debug('Using ipyparallel with %d engines', len(rc))
        lview = rc.load_balanced_view()

        def mapper(func, it):
            start_time = datetime.now()
            ret = lview.map_sync(func, it)
            if args.start or args.continu:
                logger.debug('Generation took %s', datetime.now() - start_time)

            # Save timing information for each generation
            f =  open('time_info.txt','a')
            f.write('%s\n'%(datetime.now() - start_time))
            f.close()

            return ret

        map_function = mapper
    else:
        map_function = None

    seed = os.getenv('BLUEPYOPT_SEED', args.seed)

    # load the configuration paths
    path_data=utility.load_json(args.config_path)

    morph_path = path_data['morphology']
    protocol_path = path_data['protocols']
    mech_path = path_data['mechanism']
    feature_path = path_data['features']
    param_path = path_data['parameters']
    eval_handler = Bpopt_Evaluator(protocol_path, feature_path, morph_path,
                                    param_path, mech_path, timeout = args.timeout,
                                    learn_eval_trend = args.learn_eval_trend)
    evaluator = eval_handler.create_evaluator()

    opt = bpopt.optimisations.DEAPOptimisation(
            evaluator=evaluator,
            map_function=map_function,
            seed=seed)
    return opt


def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='AIBS single neuron optimization',
        epilog=textwrap.dedent('''\
    The folling environment variables are considered:
    IPYTHON_PROFILE: if set, used as the path to the ipython profile
    BLUEPYOPT_SEED: The seed used for initial randomization
        '''))
    parser.add_argument('--start', action="store_true")
    parser.add_argument('--continu', action="store_true", default=False)
    parser.add_argument('--checkpoint', required=False, default=None,
                        help='Checkpoint pickle to avoid recalculation')
    parser.add_argument('--cp_backup', required=False, default=None,
                        help='Checkpoint backup to avoid corrupt pickling')
    parser.add_argument('--cp_backup_frequency', type=int, required=False, default=5,
                        help='Checkpoint backup frequency')
    parser.add_argument('--config_path', required=False, default='config_file.json',
                        help='For user defined configuration path for optimization')
    parser.add_argument('--offspring_size', type=int, required=False, default=2,
                        help='number of individuals in offspring')
    parser.add_argument('--max_ngen', type=int, required=False, default=2,
                        help='maximum number of generations')
    parser.add_argument('--seed', type=int, default=1,
                        help='Seed to use for optimization')
    parser.add_argument('--ipyparallel', action="store_true", default=False,
                        help='Use ipyparallel')
    parser.add_argument('-v', '--verbose', action='count', dest='verbose',
                        default=0, help='-v for INFO, -vv for DEBUG')
    parser.add_argument('--timeout', type=int, required=False, default=900,
                        help='Simulation cut-off time in seconds')
    parser.add_argument('--learn_eval_trend', action="store_true", default=False,
                        help='Modify the timeout based on evaluation times of previous generation')
    return parser



def main():
    """Main"""
    args = get_parser().parse_args()

    if args.verbose > 2:
        sys.exit('cannot be more verbose than -vv')
    logging.basicConfig(level=(logging.WARNING,
                               logging.INFO,
                               logging.DEBUG)[args.verbose],
                               stream=sys.stdout)
    opt = create_optimizer(args)

    if args.start or args.continu:
        logger.debug('Doing start or continue')
        
        if os.path.exists(args.checkpoint):
            try:
                _ = utility.load_pickle(args.checkpoint)
            except:
                logger.debug('Checkpoint file is corrupt! Looking for backup')
                shutil.copyfile(args.cp_backup,args.checkpoint)
            
        opt.run(max_ngen=args.max_ngen,
                offspring_size=args.offspring_size,
                continue_cp=args.continu,
                cp_filename=args.checkpoint,
                cp_backup=args.cp_backup,
                cp_backup_frequency=args.cp_backup_frequency)

       


if __name__ == '__main__':
    main()
