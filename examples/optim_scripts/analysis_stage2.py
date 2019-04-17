import os,glob
from ateamopt.utils import utility
from ateamopt.analysis.optim_analysis import Optim_Analyzer
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
import bluepyopt as bpopt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import argparse
from ateamopt.analysis import analysis_module
import logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger()

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

def get_opt_obj(protocols_path,features_path,morph_path, param_path,
                mech_path,map_function):
    eval_handler = Bpopt_Evaluator(protocols_path, features_path,
                                   morph_path, param_path,
                                   mech_path,timed_evaluation = False)
    evaluator = eval_handler.create_evaluator()
    
    opt = bpopt.optimisations.DEAPOptimisation(
            evaluator=evaluator,map_function=map_function)
    return opt


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
    perisomatic_model_id = cell_metadata['Perisomatic_id']

    opt_config_filename = args.config_path
    opt_config = utility.load_json(opt_config_filename)

    all_protocols_write_path = opt_config['all_protocols']
    features_write_path = opt_config['features']
    all_feature_path = opt_config['all_features']
    untrained_feature_path = opt_config['untrained_features']
    morph_path = opt_config['morphology']
    param_write_path = opt_config['parameters']
    mech_write_path = opt_config['mechanism']
    release_param_write_path = opt_config['released_model']
    mech_release_write_path = opt_config['mechanism_release']

    map_function = analyzer_map(args.ipyparallel)
    opt_train = get_opt_obj(all_protocols_write_path, features_write_path,
                                   morph_path, param_write_path,
                                   mech_write_path,map_function)
    opt_gen = get_opt_obj(all_protocols_write_path, all_feature_path,
                                   morph_path, param_write_path,
                                   mech_write_path,map_function)
    opt_untrain = get_opt_obj(all_protocols_write_path,untrained_feature_path,
                                   morph_path, param_write_path,
                                   mech_write_path,map_function)

    analysis_handler = Optim_Analyzer(opt_train,args.cp_dir)
    best_model = analysis_handler.get_best_model() # Model with least training error

    aibs_params_modelname = 'fitted_params/optim_param_%s.json'%cell_id
    analysis_handler.save_params_aibs_format(aibs_params_modelname,
                                    best_model[0],expand_params = True)

    bpopt_params_modelname = 'fitted_params/optim_param_%s_bpopt.json'%cell_id
    analysis_handler.save_params_bpopt_format(bpopt_params_modelname,
                                    best_model[0])

    hof_model_params,seed_indices = analysis_handler.get_all_models()
    hof_params_filename = 'analysis_params/hof_model_params.pkl'
    hof_responses_filename = 'analysis_params/hof_response_all.pkl'
    obj_list_train_filename = 'analysis_params/hof_obj_train.pkl'
    obj_list_all_filename = 'analysis_params/hof_obj_all.pkl'
    obj_list_untrain_filename = 'analysis_params/hof_obj_untrain.pkl'
    seed_indices_filename = 'analysis_params/seed_indices.pkl'

    # Response for the entire hall of fame not arranged
    hof_response_list = analysis_handler.get_model_responses(hof_model_params,
                                                     hof_responses_filename)

    # Sort everything with respect to training error
    obj_list_train = analysis_handler.get_response_scores(hof_response_list)
    
    analysis_handler._opt = opt_gen
    obj_list_gen = analysis_handler.get_response_scores(hof_response_list)
    
    analysis_handler._opt = opt_untrain
    obj_list_untrain = analysis_handler.get_response_scores(hof_response_list)
    
    analysis_handler._opt = opt_train
    score_list_train = [np.sum(list(obj_dict_train.values())) \
                        for obj_dict_train in obj_list_train]

    hof_response_sorted = analysis_handler.organize_models(hof_response_list,
                                                                score_list_train)
    seed_indices_sorted = analysis_handler.organize_models(seed_indices,
                                                                score_list_train)
    obj_list_train_sorted = analysis_handler.organize_models(obj_list_train,
                                                                score_list_train)
    obj_list_gen_sorted = analysis_handler.organize_models(obj_list_gen,
                                                                score_list_train)
    obj_list_untrain_sorted = analysis_handler.organize_models(obj_list_untrain,
                                                                score_list_train)

    # Save the sorted hall of fame output in .pkl
    hof_model_params_sorted = analysis_handler.save_hof_output_params(hof_model_params,\
                                                              hof_params_filename)
    utility.save_pickle(hof_responses_filename, hof_response_sorted)
    utility.save_pickle(obj_list_train_filename, obj_list_train_sorted)
    utility.save_pickle(obj_list_all_filename, obj_list_gen_sorted)
    utility.save_pickle(obj_list_untrain_filename, obj_list_untrain_sorted)
    utility.save_pickle(seed_indices_filename, seed_indices_sorted)

    # Save the entire hall of fame parameters
    for i,hof_param in enumerate(hof_model_params_sorted):
        utility.save_json('fitted_params/hof_param_%s_%s.json'%(cell_id,i),hof_param)

    GA_evol_path = 'analysis_params/GA_evolution_params.pkl'
    analysis_handler.save_GA_evolultion_info(GA_evol_path)

    resp_filename = './resp_opt.txt'
    analysis_handler.save_best_response(hof_response_sorted[0], resp_filename)

    if release_param_write_path:
        eval_handler_release = Bpopt_Evaluator(all_protocols_write_path,
                                   features_write_path,
                                   morph_path, release_param_write_path,
                                   mech_release_write_path,
                                   stub_axon = False,
                                   do_replace_axon = True,
                                   timed_evaluation = False)
        evaluator_release = eval_handler_release.create_evaluator()
        opt_release = bpopt.optimisations.DEAPOptimisation(
                            evaluator=evaluator_release)
    else:
        opt_release = None

    resp_release_filename = './resp_release.txt'
    analysis_handler.get_release_responses(opt_release,resp_release_filename)

    stim_mapfile = 'preprocessed/StimMapReps.csv'
    analysis_write_path = '%s_Stage2.pdf'%cell_id
    pdf_pages =  PdfPages(analysis_write_path)
    model_type = 'All-active'
    pdf_pages= analysis_handler.plot_grid_Response(resp_filename,
                                        resp_release_filename,
                                        stim_mapfile,
                                        pdf_pages)
    
    pdf_pages= analysis_handler.plot_feature_comp(resp_filename,
                         resp_release_filename, pdf_pages)
    pdf_pages = analysis_handler.plot_GA_evol(GA_evol_path,pdf_pages)
    pdf_pages = analysis_handler.plot_param_diversity(hof_params_filename,
                                 pdf_pages)
    exp_fi_path = 'Validation_Responses/fI_exp_%s.pkl'%cell_id 
    model_fi_path = 'Validation_Responses/fI_aa_%s.pkl'%cell_id
    exp_AP_shape_path = 'Validation_Responses/AP_shape_exp_%s.pkl'%cell_id
    model_AP_shape_path = 'Validation_Responses/AP_shape_aa_%s.pkl'%cell_id
                    
    pdf_pages = analysis_handler.postprocess(stim_mapfile,resp_filename,pdf_pages,\
                         exp_fi_path, model_fi_path,exp_AP_shape_path,model_AP_shape_path,
                         model_type)
    
    # Perisomatic model
    
    if perisomatic_model_id != '':
        resp_peri_filename = './resp_peri.txt'
        peri_param_write_path = opt_config['peri_parameters']
        peri_mech_write_path = opt_config['peri_mechanism']
        eval_handler_peri = Bpopt_Evaluator(all_protocols_write_path,
                                   features_write_path,
                                   morph_path, peri_param_write_path,
                                   peri_mech_write_path,
                                   timed_evaluation = False)
        evaluator_peri = eval_handler_peri.create_evaluator()
        opt_peri = bpopt.optimisations.DEAPOptimisation(
                            evaluator=evaluator_peri)
        model_type = 'Perisomatic'
        analysis_handler.get_release_responses(opt_peri,resp_peri_filename)
        pdf_pages= analysis_handler.plot_grid_Response(resp_filename,
                                        resp_peri_filename,
                                        stim_mapfile,pdf_pages,resp_comparison = model_type)    
        model_fi_path = 'Validation_Responses/fI_peri_%s.pkl'%cell_id
        model_AP_shape_path = 'Validation_Responses/AP_shape_peri_%s.pkl'%cell_id
        pdf_pages = analysis_handler.postprocess(stim_mapfile,resp_peri_filename,pdf_pages,
                             exp_fi_path, model_fi_path,exp_AP_shape_path,model_AP_shape_path,
                             model_type)
        
    spiketimes_exp_path ='Validation_Responses/spiketimes_exp_noise.pkl'
    spiketimes_hof_path = 'Validation_Responses/spiketimes_model_noise.pkl'
    exp_variance_hof_path = 'Validation_Responses/exp_variance_hof.pkl'
    model_perf_filename = 'Validation_Responses/fitness_metrics_'+cell_id+'.csv'
    pdf_pages = analysis_handler.hof_statistics(stim_mapfile, pdf_pages,
                       obj_list_all_filename, hof_responses_filename,
                       obj_list_train_filename,obj_list_untrain_filename,
                       seed_indices_filename,spiketimes_exp_path,spiketimes_hof_path,
                       exp_variance_hof_path,cell_metadata,model_perf_filename)
    pdf_pages.close()
    time_by_gen_filename = 'time_info.txt'
    if os.path.exists(time_by_gen_filename):
        time_metrics_filename = 'time_metrics_%s.csv'%cell_id
        analysis_module.save_optimization_time(time_by_gen_filename,
                        time_metrics_filename,cell_metadata)
    
if __name__ == '__main__':
    main()
