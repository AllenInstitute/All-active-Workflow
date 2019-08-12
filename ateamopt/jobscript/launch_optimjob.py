import os
import argparse
from ateamopt.cell_data import CellDataTool
from ateamopt.utils import utility
from ateamopt.jobscript.jobmodule import ChainSubJob
import shutil
import logging
# import pkg_resources

logger = logging.getLogger(__name__)


def create_optim_job(args):

    optim_dir = os.path.join(os.getcwd(), str(args.cell_id))
    utility.create_dirpath(optim_dir)
    os.chdir(optim_dir)  # Change Working directory

    cell_data = CellDataTool()
    me_props = {}

    if args.swc_path:
        me_props['swc_path'] = args.swc_path

    if args.nwb_path:
        me_props['nwb_path'] = args.nwb_path

    if args.me_type:
        me_props['me_type'] = args.me_type

    script_dir = os.path.join(os.getcwd(), 'Script_Repo')
    if os.path.exists(script_dir):
        shutil.rmtree(script_dir)
    try:
        # if args.ext_scripts is None:
            # TODO: use package resources for scripts
            # args.ext_scripts = pkg_resources.resource_filename() 
        shutil.copytree(args.ext_scripts, script_dir)

    except Exception as e:
        print(e)

    # TODO: select subset of args or structure better; for now, pass all through
    cell_metadata = cell_data.save_cell_metadata(**vars(args))
    cell_data.save_morph_data()

    machine = cell_metadata['Machine']

    if args.qos and 'cori' in machine:
        with open('qos.txt', 'w') as handle:
            handle.write(args.qos)

    jobtemplate_path = 'job_templates/Stage0_chainjob_template.sh'
    chain_job = ChainSubJob(jobtemplate_path, machine, conda_env=args.conda_env,
                            non_standard_nwb=args.non_std_nwb)
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
    parser.add_argument('--me_type', required=False, default=None,
                        help='ME type from the clustering')
    parser.add_argument('--qos', required=False, default=None,
                        help='Specify queue for NERSC')
    parser.add_argument('--conda_env', required=False, default='ateam_opt',
                        help='Specify the conda environment')
    parser.add_argument('--non_std_nwb', action="store_true", default=False,
                        help='Non standard processing of .nwb files')
    parser.add_argument('--from_lims', action="store_true", default=False,
                        help='Pull cell data from LIMS')
    args = parser.parse_args()
    create_optim_job(args)
