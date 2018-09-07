#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 15:40:38 2018

@author: anin
"""

import bluepyopt as bpopt
import argparse
import logging
import os
import sys
import textwrap
import json
from datetime import datetime
import shutil
from shutil import copyfile


import evaluator_helper
import checkpoint_decider

cp_backup = 'checkpoints_backup'
cp_source = 'checkpoints'

if os.path.exists('time_info_back_up.txt'):
    copyfile('time_info_back_up.txt', 'time_info.txt')
    

logging.basicConfig(level=logging.DEBUG) 
logger = logging.getLogger()


with open('config_file.json') as json_file:  
    path_data = json.load(json_file)
    
morph_path = path_data['morphology']
all_protocol_path = path_data['all_protocols']
protocol_path = path_data['protocols']
mech_path = path_data['mechanism']
feature_path = path_data['features']
param_path = path_data['parameters']


def create_optimizer(args):
    '''returns configured bluepyopt.optimisations.DEAPOptimisation'''
    
    if args.ipyparallel or os.getenv('CELLBENCHMARK_USEIPYP'):
        from ipyparallel import Client
        rc = Client(profile=os.getenv('IPYTHON_PROFILE'))
        
        logger.debug('Using ipyparallel with %d engines', len(rc))
        lview = rc.load_balanced_view()

        def mapper(func, it):
            start_time = datetime.now()
            ret = lview.map_sync(func, it)
            logger.debug('Generation took %s', datetime.now() - start_time)
            
            # Create a back-up checkpoint directory 
            # (to save optimization results in case checkpoint file is corrupted)
            
            if os.path.exists(cp_backup):
                shutil.rmtree(cp_backup)
            shutil.copytree(cp_source, cp_backup)
            
            # Save timing information for each generation
            f =  open('time_info.txt','a')
            f.write('%s\n'%(datetime.now() - start_time))
            f.close()
            
            # Create a back-up of the timing information
            copyfile('time_info.txt', 'time_info_back_up.txt')
            return ret

        map_function = mapper
    else:
        map_function = None
               
    seed = os.getenv('BLUEPYOPT_SEED', args.seed)    
    if args.analyse or args.short_analyse:
        
        evaluator = evaluator_helper.create(all_protocol_path, feature_path, morph_path, 
                                        param_path, mech_path)

    else:
            
        evaluator = evaluator_helper.create(protocol_path, feature_path, morph_path, 
                                      param_path, mech_path)
        
       
    
    opt = bpopt.optimisations.DEAPOptimisation(
        evaluator=evaluator,
        map_function=map_function,
        seed=seed)
    

    return opt
    


def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Cell Optimization Example',
        epilog=textwrap.dedent('''\
The folling environment variables are considered:
    CELLBENCHMARK_USEIPYP: if set, will use ipyparallel
    IPYTHON_PROFILE: if set, used as the path to the ipython profile
    BLUEPYOPT_SEED: The seed used for initial randomization
        '''))
    parser.add_argument('--start', action="store_true")
    parser.add_argument('--continu', action="store_true", default=False)
    parser.add_argument('--checkpoint', required=False, default=None,
                        help='Checkpoint pickle to avoid recalculation')
    parser.add_argument('--offspring_size', type=int, required=False, default=2,
                        help='number of individuals in offspring')
    parser.add_argument('--max_ngen', type=int, required=False, default=2,
                        help='maximum number of generations')
    parser.add_argument('--responses', required=False, default=None,
                        help='Response pickle file to avoid recalculation')
    parser.add_argument('--analyse', action="store_true")
    parser.add_argument('--short_analyse', action="store_true")
    parser.add_argument('--compile', action="store_true")
    parser.add_argument('--seed', type=int, default=1,
                        help='Seed to use for optimization')
    parser.add_argument('--ipyparallel', action="store_true", default=False,
                        help='Use ipyparallel')
    parser.add_argument('-v', '--verbose', action='count', dest='verbose',
                        default=0, help='-v for INFO, -vv for DEBUG')

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

    if args.compile:
        logger.debug('Doing compile')
        import commands
        commands.getstatusoutput('nrnivmodl modfiles/')

    if args.start or args.continu:
        logger.debug('Doing start or continue')
        opt.run(max_ngen=args.max_ngen,
                offspring_size=args.offspring_size,
                continue_cp=args.continu,
                cp_filename=args.checkpoint)
    
    if args.short_analyse:
        
        logger.debug('Doing analyse to save the results')
        from optim_analysis_short import plot_Response,plot_diversity,\
                                            plot_GA_evolution, get_responses
        if '/' in args.checkpoint:
            cp_dir = args.checkpoint.split('/')[0]
        else:
            cp_dir = '.'  # check in current directory
        args.checkpoint = checkpoint_decider.best_seed(cp_dir)
        
        
        if args.checkpoint is not None and os.path.isfile(args.checkpoint):
            hof_index = 0
            plot_Response(opt,args.checkpoint,args.responses,hof_index)
            
        else:
            print('No checkpoint file available run optimization '
                  'first with --start')

        logger.debug('Plotting Parameters - Optimized and Released')

        if not os.path.exists(args.checkpoint):
            raise Exception('Need a pickle file to plot the parameter diversity')

        plot_diversity(opt, args.checkpoint,
                                     opt.evaluator.param_names,hof_index)

        plot_GA_evolution(args.checkpoint)
        
    elif args.analyse:
        
        from optim_analysis import plot_Response,feature_comp,plot_diversity,\
                                                plot_GA_evolution,post_processing
        logger.debug('Plotting Response Comparisons')
        plot_Response(opt,args.checkpoint,args.responses)
        
        logger.debug('Plotting Feature Comparisons')
        feature_comp(opt,args.checkpoint,args.responses)
        
        logger.debug('Plotting Parameters - Optimized and Released')
        plot_diversity(opt, args.checkpoint,opt.evaluator.param_names)
        
        logger.debug('Plotting Evolution of the Objective')
        plot_GA_evolution(args.checkpoint)
        
        logger.debug('Plotting Spike shapes and mean frequency comparison')
        post_processing(args.checkpoint,args.responses)
        

if __name__ == '__main__':
    main()
