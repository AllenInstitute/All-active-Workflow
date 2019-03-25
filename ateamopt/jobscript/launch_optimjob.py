import os
import argparse
from ateamopt.allactive_optim import Allactive_Optim
from ateamopt.utils import utility
from ateamopt.jobscript.jobmodule import ChainSubJob
import shutil
import logging

logger = logging.getLogger(__name__)


def create_optim_job(args):
    
    optim_dir = os.path.join(os.getcwd(),str(args.cell_id))
    utility.create_dirpath(optim_dir)
    os.chdir(optim_dir) # Change Working directory
    
    optim_model = Allactive_Optim()
    me_props = {}
    
    if args.swc_path:  
        me_props['swc_path'] = args.swc_path
        
    if args.nwb_path:
        me_props['nwb_path'] = args.nwb_path
    
    script_dir = os.path.join(os.getcwd(),'Script_Repo')
    if os.path.exists(script_dir):
        shutil.rmtree(script_dir)
    try:    
        shutil.copytree(args.ext_scripts, script_dir)

    except Exception as e:
        print(e)
    
    cell_metadata = optim_model.save_cell_metadata(**me_props)
    machine = cell_metadata['Machine']
    
    if args.qos and 'cori' in machine:
       with open('qos.txt', 'a') as handle:
           handle.write(args.qos) 
    jobtemplate_path = 'job_templates/Stage0_chainjob_template.sh'
    chain_job = ChainSubJob(jobtemplate_path,machine)
    chain_job.script_generator()
    chain_job.run_job()
    


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cell_id', required=True, default=None,
                        help='Feed the Cell id')
    parser.add_argument('--ext_scripts', required=True, default=None,
                        help='Directory for runtime scripts')
    parser.add_argument('--swc_path', required=False, default=None,
                        help='Morphology path for unpublished cells')
    parser.add_argument('--nwb_path', required=False, default=None,
                        help='Ephys path for unpublished cells')
    parser.add_argument('--qos', required=False, default=None,
                        help='Specify queue for NERSC')
    args = parser.parse_args()
    create_optim_job(args)
