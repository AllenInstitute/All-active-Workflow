import os,glob
from ateamopt.utils import utility
from ateamopt.analysis.optim_analysis import Optim_Analyzer
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
import bluepyopt as bpopt
from matplotlib.backends.backend_pdf import PdfPages
import argparse


import logging
logger = logging.getLogger(__name__)


def analyzer_map(parallel=True):
    '''returns configured bluepyopt.optimisations.DEAPOptimisation'''

    if parallel:
        from ipyparallel import Client
        rc = Client(profile=os.getenv('IPYTHON_PROFILE'))
        
        logger.debug('Using ipyparallel with %d engines', len(rc))
        lview = rc.load_balanced_view()

        def mapper(func, it):
            ret = lview.map_sync(func, it)
            return ret

        map_function = mapper
    else:
        map_function = None
        
    return map_function

def get_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument('--cp_dir', required=False, default=None,
                        help='Directory for checkpoint pickle files')
    parser.add_argument('--config_path', required=False, default='config_file.json',
                        help='For user defined configuration path for optimization')

    parser.add_argument('--ipyparallel', action="store_true", default=False,
                        help='Use ipyparallel')
    parser.add_argument('-v', '--verbose', action='count', dest='verbose',
                        default=0, help='-v for INFO, -vv for DEBUG')

    return parser

def main():

    args = get_parser().parse_args()

    parent_dir = os.path.abspath(os.path.join('.', os.pardir))
    path_to_cell_metadata = glob.glob(parent_dir+'/cell_metadata*.json')[0] 
    cell_metadata=utility.load_json(path_to_cell_metadata)
    cell_id = cell_metadata['Cell_id']
    
    opt_config_filename = args.config_path
    opt_config = utility.load_json(opt_config_filename)
    
    all_protocols_write_path = opt_config['all_protocols']
    features_write_path = opt_config['features']
    morph_path = opt_config['morphology']
    param_write_path = opt_config['parameters']
    mech_write_path = opt_config['mechanism']
    release_param_write_path = opt_config['released_model']
    mech_release_write_path = opt_config['mechanism_release']
    
    eval_handler = Bpopt_Evaluator(all_protocols_write_path, features_write_path,
                                   morph_path, param_write_path,
                                   mech_write_path)
    evaluator = eval_handler.create_evaluator()
    map_function = analyzer_map(args.ipyparallel)
    opt_train = bpopt.optimisations.DEAPOptimisation(
            evaluator=evaluator,map_function=map_function)

    cp_dir = args.cp_dir   
    analysis_handler = Optim_Analyzer(opt_train,cp_dir)
    best_model = analysis_handler.get_best_model()
    
    aibs_params_modelname = 'fitted_params/optim_param_%s.json'%cell_id
    analysis_handler.save_params_aibs_format(aibs_params_modelname,
                                    best_model[0],expand_params = False)
    
    aibs_params_modelname = 'fitted_params/optim_param_%s_bpopt.json'%cell_id
    analysis_handler.save_params_bpopt_format(aibs_params_modelname,
                                    best_model[0])
    
    
    hof_model_params,_ = analysis_handler.get_all_models()
    hof_params_filename = 'analysis_params/hof_model_params.pkl'
    analysis_handler.save_hof_output_params(hof_model_params,hof_params_filename)
    
    
    GA_evol_path = 'analysis_params/GA_evolution_params.pkl'
    analysis_handler.save_GA_evolultion_info(GA_evol_path)
    
    response_list = analysis_handler.get_model_responses(best_model)
    resp_filename = './resp_opt.txt'
    analysis_handler.save_best_response(response_list[0], resp_filename)
    
    if release_param_write_path:
        eval_handler_release = Bpopt_Evaluator(all_protocols_write_path, 
                                   features_write_path,
                                   morph_path, release_param_write_path,
                                   mech_release_write_path,
                                   do_replace_axon = False,
                                   do_replace_axon_swc = True)
        evaluator_release = eval_handler_release.create_evaluator()
        opt_release = bpopt.optimisations.DEAPOptimisation(
                            evaluator=evaluator_release)
    else:
        opt_release = None
    
    resp_release_filename = './resp_release.txt'
    analysis_handler.get_release_responses(opt_release,resp_release_filename)    
    
    
    stim_mapfile = 'preprocessed/StimMapReps.csv'
    analysis_write_path = cell_id + '_Stage2.pdf'
    pdf_pages =  PdfPages(analysis_write_path)
        
    pdf_pages= analysis_handler.plot_grid_Response(resp_filename,
                                        resp_release_filename,
                                        stim_mapfile,
                                        pdf_pages)    
        
    pdf_pages= analysis_handler.plot_feature_comp(resp_filename,
                         resp_release_filename, pdf_pages)
    
    pdf_pages = analysis_handler.plot_GA_evol(GA_evol_path,pdf_pages)
    pdf_pages = analysis_handler.plot_param_diversity(hof_params_filename,
                                 pdf_pages)
    pdf_pages.close()
    
    
if __name__ == '__main__':
    main()    
    