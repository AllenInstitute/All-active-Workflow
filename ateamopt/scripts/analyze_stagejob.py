import os
import glob
from ateamopt.utils import utility
from ateamopt.analysis.optim_analysis import Optim_Analyzer
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
import bluepyopt as bpopt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from ateamopt.analysis import analysis_module
from ateamopt.optim_schema import Optim_Config
import argschema as ags
import logging

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


def get_opt_obj(protocols_path, features_path, morph_path, param_path,
                mech_path, map_function, **kwargs):
    eval_handler = Bpopt_Evaluator(protocols_path, features_path,
                                   morph_path, param_path,
                                   mech_path, **kwargs)
    evaluator = eval_handler.create_evaluator()

    opt = bpopt.optimisations.DEAPOptimisation(
        evaluator=evaluator, map_function=map_function)
    return opt


def main(args):
    stage_jobconfig = args['stage_jobconfig']
    highlevel_job_props = args['highlevel_jobconfig']

    logging.basicConfig(level=highlevel_job_props['log_level'])

    parent_dir = os.path.abspath(os.path.join('.', os.pardir))
    path_to_cell_metadata = glob.glob(parent_dir+'/cell_metadata*.json')[0]
    cell_metadata = utility.load_json(path_to_cell_metadata)
    cell_id = cell_metadata['cell_id']
    peri_model_id = cell_metadata.get('peri_model_id')
    released_aa_model_id = cell_metadata.get('released_aa_model_id')

    all_protocols_path = highlevel_job_props['all_protocols_path']
    all_features_path = highlevel_job_props['all_features_path']
    morph_path = highlevel_job_props['swc_path']
    axon_type = highlevel_job_props['axon_type']
    ephys_dir = highlevel_job_props['ephys_dir']

    train_features_path = args['train_features']
    test_features_path = args['test_features']
    param_write_path = args['parameters']
    mech_write_path = args['mechanism']
    release_param_write_path = args['released_aa_model']
    mech_release_write_path = args['released_aa_mechanism']

    analysis_parallel = (stage_jobconfig['analysis_config'].get('ipyparallel') and 
            stage_jobconfig['run_hof_analysis'])

    props = dict(axon_type=axon_type, ephys_dir=ephys_dir)

    map_function = analyzer_map(analysis_parallel)
    opt_train = get_opt_obj(all_protocols_path, train_features_path,
                            morph_path, param_write_path,
                            mech_write_path, map_function, **props)
    opt_all = get_opt_obj(all_protocols_path, all_features_path,
                          morph_path, param_write_path,
                          mech_write_path, map_function, **props)
    opt_test = get_opt_obj(all_protocols_path, test_features_path,
                           morph_path, param_write_path,
                           mech_write_path, map_function, **props)

    analysis_handler = Optim_Analyzer(args, opt_all)
    best_model = analysis_handler.get_best_model()  # Model with least training error

    aibs_params_modelname = 'fitted_params/optim_param_%s.json' % cell_id
    analysis_handler.save_params_aibs_format(aibs_params_modelname,
                                             best_model[0],expand_params = True)
    aibs_params_compact_modelname = 'fitted_params/optim_param_%s_compact.json'%cell_id
    analysis_handler.save_params_aibs_format(aibs_params_compact_modelname,best_model[0])
    bpopt_params_modelname = 'fitted_params/optim_param_%s_bpopt.json'%cell_id
    analysis_handler.save_params_bpopt_format(bpopt_params_modelname,
                                              best_model[0])
    
    # Export hoc model
    if stage_jobconfig['hoc_export']:
        hoc_export_path = 'fitted_params/model_template_%s.hoc'%cell_id
        utility.create_filepath(hoc_export_path)
        best_param_dict = {key:best_model[0][i] for i,key in \
                            enumerate(opt_train.evaluator.param_names)}
        model_string = opt_train.evaluator.cell_model.create_hoc(best_param_dict)
        with open(hoc_export_path, "w") as hoc_template:
            hoc_template.write(model_string)
    
    hof_model_params, seed_indices = analysis_handler.get_all_models()

    if not stage_jobconfig.get('run_hof_analysis'):
        hof_model_params, seed_indices = best_model, [seed_indices[0]]

    hof_params_filename = 'analysis_params/hof_model_params.pkl'
    hof_responses_filename = 'analysis_params/hof_response_all.pkl'
    obj_list_train_filename = 'analysis_params/hof_obj_train.pkl'
    obj_list_all_filename = 'analysis_params/hof_obj_all.pkl'
    feat_list_all_filename = 'analysis_params/hof_features_all.pkl'
    obj_list_test_filename = 'analysis_params/hof_obj_test.pkl'
    seed_indices_filename = 'analysis_params/seed_indices.pkl'

    score_list_train_filename = 'analysis_params/score_list_train.pkl'

    # Response for the entire hall of fame not arranged
    hof_response_list = analysis_handler.get_model_responses(hof_model_params,
                                                             hof_responses_filename)

    analysis_handler._opt = opt_train
    obj_list_train = analysis_handler.get_response_scores(hof_response_list)

    # Sort everything with respect to training error

    if not os.path.exists(score_list_train_filename):
        score_list_train = [np.sum(list(obj_dict_train.values()))
                            for obj_dict_train in obj_list_train]
        utility.create_filepath(score_list_train_filename)
        utility.save_pickle(score_list_train_filename, score_list_train)
    else:
        score_list_train = utility.load_pickle(score_list_train_filename)

    if not os.path.exists(seed_indices_filename):
        seed_indices_sorted = analysis_handler.organize_models(seed_indices,
                                                               score_list_train)
        utility.save_pickle(seed_indices_filename, seed_indices_sorted)

    if not os.path.exists(obj_list_train_filename):
        obj_list_train_sorted = analysis_handler.organize_models(obj_list_train,
                                                                 score_list_train)
        utility.save_pickle(obj_list_train_filename, obj_list_train_sorted)

    if not os.path.exists(obj_list_all_filename):
        analysis_handler._opt = opt_all
        obj_list_gen = analysis_handler.get_response_scores(hof_response_list)
        obj_list_gen_sorted = analysis_handler.organize_models(obj_list_gen,
                                                               score_list_train)
        utility.save_pickle(obj_list_all_filename, obj_list_gen_sorted)

    if not os.path.exists(feat_list_all_filename):
        analysis_handler._opt = opt_all
        feat_list_gen = analysis_handler.get_response_features(
            hof_response_list)
        feat_list_gen_sorted = analysis_handler.organize_models(feat_list_gen,
                                                                score_list_train)
        utility.save_pickle(feat_list_all_filename, feat_list_gen_sorted)

    if not os.path.exists(obj_list_test_filename):
        analysis_handler._opt = opt_test
        obj_list_test = analysis_handler.get_response_scores(hof_response_list)

        obj_list_test_sorted = analysis_handler.organize_models(obj_list_test,
                                                                score_list_train)
        utility.save_pickle(obj_list_test_filename, obj_list_test_sorted)

    analysis_handler._opt = opt_train

    # Save the sorted hof responses at the end
    hof_response_sorted = analysis_handler.organize_models(hof_response_list,
                                                           score_list_train)
    utility.save_pickle(hof_responses_filename, hof_response_sorted)

    # Save the sorted hall of fame output in .pkl
    hof_model_params_sorted = analysis_handler.save_hof_output_params(hof_model_params,
                                                  hof_params_filename, score_list_train)
    # Save the entire hall of fame parameters
    for i, hof_param in enumerate(hof_model_params_sorted):
        aibs_params_modelname = os.path.join('fitted_params','hof_param_%s_%s.json' % (
            cell_id,i))
        analysis_handler.save_params_aibs_format(aibs_params_modelname,
                                                 hof_param, expand_params=True)

    # Now save the sorted score
    utility.save_pickle(score_list_train_filename, sorted(score_list_train))

    GA_evol_path = os.path.join('analysis_params','GA_evolution_params.pkl')
    analysis_handler.save_GA_evolultion_info(GA_evol_path)

    resp_filename = os.path.join(os.getcwd(),'resp_opt.txt')
    analysis_handler.save_best_response(hof_response_sorted[0], resp_filename)

    if release_param_write_path:
        eval_handler_release = Bpopt_Evaluator(all_protocols_path,
                                               all_features_path,
                                               morph_path, release_param_write_path,
                                               mech_release_write_path,
                                               stub_axon=False,
                                               do_replace_axon=True,
                                               ephys_dir=ephys_dir)
        evaluator_release = eval_handler_release.create_evaluator()
        opt_release = bpopt.optimisations.DEAPOptimisation(
            evaluator=evaluator_release)
        resp_release_filename = os.path.join(os.getcwd(),'resp_release.txt')
        analysis_handler.get_release_responses(
            opt_release, resp_release_filename)
        resp_release_aa = utility.load_pickle(resp_release_filename)[0]
        features_release_aa = opt_release.evaluator.fitness_calculator.\
            calculate_features(resp_release_aa)
        features_aa_filename = os.path.join('Validation_Responses','Features_released_aa_%s.pkl'% cell_id)
        utility.create_filepath(features_aa_filename)
        utility.save_pickle(features_aa_filename, features_release_aa)
    else:
        resp_release_filename = None

    stim_mapfile = highlevel_job_props['stimmap_file']
    analysis_write_path = '%s_%s.pdf' % (
        cell_id, stage_jobconfig['stage_name'])
    pdf_pages = PdfPages(analysis_write_path)
    model_type = 'All-active'
    pdf_pages = analysis_handler.plot_grid_Response(resp_filename,
                                                    resp_release_filename,
                                                    stim_mapfile,
                                                    pdf_pages)

    pdf_pages = analysis_handler.plot_feature_comp(resp_filename,
                                                   resp_release_filename, pdf_pages)
    pdf_pages = analysis_handler.plot_GA_evol(GA_evol_path, pdf_pages)
    pdf_pages = analysis_handler.plot_param_diversity(hof_params_filename,
                                                      pdf_pages)

    if stage_jobconfig['model_postprocess']:
        exp_fi_path = os.path.join('Validation_Responses','fI_exp_%s.pkl'%cell_id)
        model_fi_path = os.path.join('Validation_Responses','fI_aa_%s.pkl' % cell_id)
        exp_AP_shape_path = os.path.join('Validation_Responses','AP_shape_exp_%s.pkl' % cell_id)
        model_AP_shape_path = os.path.join('Validation_Responses','AP_shape_aa_%s.pkl' % cell_id)

        pdf_pages = analysis_handler.postprocess(stim_mapfile, resp_filename, pdf_pages,
                                                 exp_fi_path, model_fi_path, exp_AP_shape_path, model_AP_shape_path,
                                                 model_type)

    # Perisomatic model

    if peri_model_id and stage_jobconfig['run_peri_comparison']:
        resp_peri_filename = os.path.join(os.getcwd(),'resp_peri.txt')
        peri_param_path = args['released_peri_model']
        peri_mech_path = args['released_peri_mechanism']

        props_peri = props.copy()
        props_peri['axon_type'] = 'stub_axon'
        eval_handler_peri = Bpopt_Evaluator(all_protocols_path,
                                            all_features_path,
                                            morph_path, peri_param_path,
                                            peri_mech_path, **props_peri)
        evaluator_peri = eval_handler_peri.create_evaluator()
        opt_peri = bpopt.optimisations.DEAPOptimisation(
            evaluator=evaluator_peri)
        model_type = 'Perisomatic'
        analysis_handler.get_release_responses(opt_peri, resp_peri_filename)
        resp_peri = utility.load_pickle(resp_peri_filename)[0]
        features_peri = opt_peri.evaluator.fitness_calculator.calculate_features(
            resp_peri)
        features_peri_filename = os.path.join('Validation_Responses','Features_peri_%s.pkl' % cell_id)
        utility.create_filepath(features_peri_filename)
        utility.save_pickle(features_peri_filename, features_peri)

        pdf_pages = analysis_handler.plot_grid_Response(resp_filename,
                        resp_peri_filename,stim_mapfile, pdf_pages, 
                        resp_comparison=model_type)
        if stage_jobconfig['model_postprocess']:
            model_fi_path = os.path.join('Validation_Responses','fI_peri_%s.pkl' % cell_id)
            model_AP_shape_path = os.path.join('Validation_Responses','AP_shape_peri_%s.pkl'%cell_id)
            pdf_pages = analysis_handler.postprocess(stim_mapfile, resp_peri_filename, 
                        pdf_pages, exp_fi_path, model_fi_path, exp_AP_shape_path, 
                        model_AP_shape_path, model_type)

    if stage_jobconfig.get('calc_model_perf'):
        spiketimes_exp_path = os.path.join('Validation_Responses','spiketimes_exp_noise.pkl')
        all_features = utility.load_json(all_features_path)
        spiketimes_noise_exp = {}
        for stim_, feat in all_features.items():
            if 'Noise' in stim_:
                if 'peak_time' in feat['soma'].keys():
                    spiketimes_noise_exp[stim_] = feat['soma']['peak_time'][2]
        if bool(spiketimes_noise_exp):
            utility.create_filepath(spiketimes_exp_path)
            utility.save_pickle(spiketimes_exp_path, spiketimes_noise_exp)

        spiketimes_hof_path = os.path.join('Validation_Responses','spiketimes_model_noise.pkl')
        exp_variance_hof_path = os.path.join('Validation_Responses','exp_variance_hof.pkl')
        model_perf_filename = os.path.join('Validation_Responses','fitness_metrics_%s.csv'%cell_id)
        pdf_pages = analysis_handler.hof_statistics(stim_mapfile, pdf_pages,
                                                    obj_list_all_filename, hof_responses_filename,
                                                    obj_list_train_filename, obj_list_test_filename,
                                                    seed_indices_filename, spiketimes_exp_path, spiketimes_hof_path,
                                                    exp_variance_hof_path, cell_metadata, model_perf_filename)
    pdf_pages.close()

    if stage_jobconfig.get('calc_time_statistics'):
#        time_by_gen_filename = 'time_info.txt'
#        if os.path.exists(time_by_gen_filename):
#            time_metrics_filename = 'time_metrics_%s.csv' % cell_id
#            analysis_module.save_optimization_time(time_by_gen_filename,
#                                                   time_metrics_filename, cell_metadata)
        compute_statistics_filename = 'compute_metrics_%s.csv' % cell_id
        opt_logbook = 'logbook_info.txt'
        analysis_module.save_compute_statistics(opt_logbook,compute_statistics_filename)

if __name__ == '__main__':
    mod = ags.ArgSchemaParser(schema_type=Optim_Config)
    main(mod.args)
