import os,sys
import glob
from ateamopt.nwb_extractor import NWB_Extractor
from ateamopt.model_parameters import AllActive_Model_Parameters
from ateamopt.utils import utility
from ateamopt.optim_config_rules import filter_feat_proto_passive
from ateamopt.analysis.optim_analysis import Optim_Analyzer
from ateamopt.bpopt_evaluator import Bpopt_Evaluator
import bluepyopt as bpopt
from ateamopt.jobscript.jobmodule import test_JobModule,\
            PBS_JobModule,Slurm_JobModule,ChainSubJob
from matplotlib.backends.backend_pdf import PdfPages
import shutil
import logging

logger = logging.getLogger(__name__)


def main():
    
    parent_dir = os.path.abspath(os.path.join('.', os.pardir))
    path_to_cell_metadata = glob.glob(parent_dir+'/*.json')[0] 
    cell_metadata=utility.load_json(path_to_cell_metadata)
    
    acceptable_stimtypes = ['Long Square']
    
    cell_id = cell_metadata['Cell_id']
    
    # Extract data and get the features for the stage
    nwb_handler = NWB_Extractor(cell_id)
    ephys_data_path,stimmap_filename = nwb_handler.save_cell_data(acceptable_stimtypes,
                                                    non_standard_nwb = True)
    feature_path = 'parameters/feature_set_stage0.json'
    filter_rule_func = filter_feat_proto_passive 
    features_write_path,untrained_features_write_path,\
        protocols_write_path,all_protocols_write_path = \
        nwb_handler.get_ephys_features(feature_path,ephys_data_path,
                                       stimmap_filename,filter_rule_func)

    # Create the parameter bounds for the optimization
    model_params_handler = AllActive_Model_Parameters(cell_id)
    morph_path = model_params_handler.swc_path
    param_bounds_path = 'parameters/param_bounds_stage0.json'
    model_params,model_params_release= model_params_handler.get_opt_params(param_bounds_path)
    param_write_path,release_param_write_path,release_params=\
                        model_params_handler.write_params_opt(model_params,model_params_release)
    mech_write_path,mech_release_write_path = model_params_handler.write_mechanisms_opt(model_params,\
                                        model_params_release,param_bounds_path)
    
    # Config file with all the necessary paths to feed into the optimization
    model_params_handler.write_opt_config_file(morph_path,param_write_path,
                                  mech_write_path,mech_release_write_path,
                                  features_write_path,untrained_features_write_path,
                                  protocols_write_path,all_protocols_write_path,
                                  release_params,release_param_write_path)

    # Copy the optimer scripts in the current directory
    
    optimizer_script=utility.locate_script_file('Optim_Main.py')
    stage_cwd = os.getcwd()
    
    for script_path in [optimizer_script]:
        shutil.copy(script_path,stage_cwd)

    # Create batch jobscript
    machine = cell_metadata['Machine']
    if 'hpc-login' in machine:
        jobtemplate_path = 'job_templates/Stage0_pbs.sh'
        batch_job = PBS_JobModule(jobtemplate_path,machine)
        batch_job.script_generator()
    elif any(substring in machine for substring in ['cori', 'bbp']):
        jobtemplate_path = 'job_templates/Stage0_slurm.sh'
        batch_job = Slurm_JobModule(jobtemplate_path,machine)
        batch_job.script_generator()
    
    
    # Create Chain job for next stage
    chain_jobtemplate_path = 'job_templates/Stage1_chainjob_template.sh'
    chain_job = ChainSubJob(chain_jobtemplate_path,machine)
    chain_job.script_generator()
    
    if 'cori' in machine:
        chain_job.adjust_for_NERSC('#SBATCH -p prod', '#SBATCH -q regular')
        chain_job.adjust_for_NERSC('#SBATCH -C cpu|nvme', '#SBATCH -C haswell')                      
        chain_job.adjust_for_NERSC('#SBATCH -A proj36','#SBATCH -L SCRATCH')
        chain_job.adjust_for_NERSC('#SBATCH -n 256', '#SBATCH -N 8') 
        Path_append ='export PATH="/global/common/software/m2043/AIBS_Opt/software/x86_64/bin:$PATH"'
        chain_job.adjust_for_NERSC('source activate %s'%chain_job.conda_env, Path_append, 
                                  add = True)    
    
    
    # Trial run
    
    if sys.argv[-1] == 'dryrun':
        
        ### Create Jobscript
        
        cp_dir = './' 
        machine = cell_metadata['Machine'] 
        testJob = test_JobModule(machine,'test_job.sh','%s/seed1.pkl'%cp_dir,
                                 2,2)
        testJob.script_generator()
        testJob.run_job()
    
        ### Analysis
            
        eval_handler = Bpopt_Evaluator(all_protocols_write_path, features_write_path,
                                       morph_path, param_write_path,
                                       mech_write_path)
        evaluator = eval_handler.create_evaluator()
            
        opt_train = bpopt.optimisations.DEAPOptimisation(
                    evaluator=evaluator)
        
           
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
        
        response_list = analysis_handler.get_model_responses(best_model)
        resp_filename = './resp_opt.txt'
        analysis_handler.save_best_response(response_list[0], resp_filename)  
        
        if release_param_write_path:
            eval_handler_release = Bpopt_Evaluator(all_protocols_write_path, 
                                       features_write_path,
                                       morph_path, release_param_write_path,
                                       mech_write_path,
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
    
    
    
    
    
    
    
    
    
    
    
    
    