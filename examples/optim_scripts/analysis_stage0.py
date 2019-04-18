import os,glob
from ateamopt.utils import utility
from ateamopt.analysis.optim_analysis import Optim_Analyzer
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
import bluepyopt as bpopt
from matplotlib.backends.backend_pdf import PdfPages

import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

def main():

    parent_dir = os.path.abspath(os.path.join('.', os.pardir))
    path_to_cell_metadata = glob.glob(parent_dir+'/*.json')[0]
    cell_metadata=utility.load_json(path_to_cell_metadata)
    cell_id = cell_metadata['Cell_id']

    opt_config_filename = 'config_file.json'
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
                                   mech_write_path,
                                   'timed_evaluation' = False)
    evaluator = eval_handler.create_evaluator()

    opt_train = bpopt.optimisations.DEAPOptimisation(
                evaluator=evaluator)

    cp_dir = 'checkpoints/'
    analysis_handler = Optim_Analyzer(opt_train,cp_dir)
    best_model = analysis_handler.get_best_model()

    aibs_params_modelname = 'fitted_params/optim_param_%s.json'%cell_id
    analysis_handler.save_params_aibs_format(aibs_params_modelname,
                                    best_model[0])

    aibs_params_modelname = 'fitted_params/optim_param_%s_bpopt.json'%cell_id
    analysis_handler.save_params_bpopt_format(aibs_params_modelname,
                                    best_model[0])


    hof_model_params,_ = analysis_handler.get_all_models()
    hof_params_filename = 'analysis_params/hof_model_params.pkl'
    analysis_handler.save_hof_output_params(hof_model_params,hof_params_filename)


    GA_evol_path = 'analysis_params/GA_evolution_params.pkl'
    analysis_handler.save_GA_evolultion_info(GA_evol_path)

    hof_responses_filename = 'analysis_params/hof_response_all.pkl'
    response_list = analysis_handler.get_model_responses(best_model,hof_responses_filename)
    resp_filename = './resp_opt.txt'
    analysis_handler.save_best_response(response_list[0], resp_filename)

    if release_param_write_path:
        eval_handler_release = Bpopt_Evaluator(all_protocols_write_path,
                                   features_write_path,
                                   morph_path, release_param_write_path,
                                   mech_release_write_path,
                                   do_replace_axon = False,
                                   do_replace_axon_swc = True,
                                   'timed_evaluation' = False)
        evaluator_release = eval_handler_release.create_evaluator()
        opt_release = bpopt.optimisations.DEAPOptimisation(
                            evaluator=evaluator_release)
    else:
        opt_release = None

    resp_release_filename = './resp_release.txt'
    analysis_handler.get_release_responses(opt_release,resp_release_filename)


    stim_mapfile = 'preprocessed/StimMapReps.csv'
    analysis_write_path = cell_id + '_Stage0.pdf'
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
