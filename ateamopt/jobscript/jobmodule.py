import os
import logging
from subprocess import Popen, PIPE
from ateamopt.utils import utility
import re

logger = logging.getLogger(__name__)


def script_decorator(func):
    def func_wrapper(*args, **kwargs):
         job_string= func(*args, **kwargs)
         return job_string
    return func_wrapper


dryrun_config = dict(offspring_size = 2,max_ngen= 2,cp_backup_dir= None,
            nengines= 2,nnodes=1,nnodes_analysis=1,nprocs=2,
            nprocs_analysis=2,nengines_analysis= 2,
            jobtime='10:00',jobtime_analysis='10:00')

class JobModule(object):

    def __init__(self, script_name = 'Jobscript.sh'):
        
        self.script_name = script_name

    def adjust_template(self,match_line, replace_line, add = False,
                        partial_match = False, add_in_place = False):
        with open(self.script_name, "r") as in_file:
            buf = in_file.readlines()

            with open(self.script_name, "w") as out_file:
                for line in buf:
                    match_eval = line == "%s\n"%match_line \
                        if not partial_match else match_line in line

                    if match_eval:
                        if add:
                            line +=  "%s\n"%replace_line
                        elif add_in_place:
                            line = line.rstrip()
                            line +=  " %s\n"%replace_line
                        else:
                            line = "%s\n"%replace_line

                    out_file.write(line)
                    

class ChainSubJob(JobModule):

    def __init__(self, script_template, job_config_path,
                     script_name = 'chain_job.sh'):
        
        super(ChainSubJob,self).__init__(script_name)
        self.script_template = utility.locate_template_file(script_template)
        self.job_config_path = job_config_path
        
  
    def script_generator(self):
        all_config = utility.load_json(self.job_config_path)
        stage_jobconfig = all_config['stage_jobconfig']
        highlevel_job_props = all_config['highlevel_jobconfig']
        
        
        machine = highlevel_job_props['machine']
        if any(substring in machine for substring in ['aws', 'hpc-login']):
            self.submit_cmd = 'qsub'            
        elif any(substring in machine for substring in ['cori', 'bbp']):
            self.submit_cmd = 'sbatch'
        else:
            self.submit_cmd = 'bash'
        conda_env = highlevel_job_props.get('conda_env','ateam_opt')
        with open(self.script_template,'r') as job_template:
            subjob_string = job_template.read()

        subjob_string = subjob_string.replace('conda_env',conda_env)
    
        if stage_jobconfig.get('stage_name'):
            stage_path = os.path.join(highlevel_job_props['job_dir'],
                                      stage_jobconfig['stage_name'])
            subjob_string = subjob_string.replace('stage_jobdir',stage_path)
            subjob_string = subjob_string.replace(' Stage ',' %s '%stage_jobconfig['stage_name'])
            subjob_string = subjob_string.replace('submit_cmd',self.submit_cmd)
        
        if highlevel_job_props.get('modfiles_dir') and not highlevel_job_props.get\
                                ('compiled_modfiles_dir'):
            subjob_string = subjob_string.replace('modfiles_dir_abs',
                          os.path.basename(highlevel_job_props['modfiles_dir']))
            subjob_string = subjob_string.replace('modfiles_dir',
                                      highlevel_job_props['modfiles_dir'])
            
            subjob_string = re.sub('# Copy compiled[\S\s]*fi','',
                                   subjob_string) 
        elif highlevel_job_props.get('compiled_modfiles_dir'):
            subjob_string = re.sub('# Copy modfiles[\S\s]*modfiles_dir_abs; fi','',
                                   subjob_string)  
            subjob_string = subjob_string.replace('compiled_modfiles_dir',
                          highlevel_job_props['compiled_modfiles_dir'])
        else:
            subjob_string = re.sub('# Copy modfiles[\S\s]*fi','',
                                   subjob_string) 
            
 
        with open(self.script_name, "w") as chainsubjob_script:
            chainsubjob_script.write(subjob_string)
            
            
        if bool(stage_jobconfig) and self.submit_cmd == 'bash':
            self.adjust_template('RES=$(%s batch_job.sh)'%self.submit_cmd, 
                                 '%s batch_job.sh'%self.submit_cmd)
            self.adjust_template('echo ${RES##* }', '',partial_match = True)
        

#    def script_generator(self):
#        # Adjusting the job based on machine
#        job_config = utility.load_json(self.job_config_path)
#        machine = job_config['machine']
#        if 'cori' in machine:
#            self.adjust_template('#SBATCH -p prod', '#SBATCH -q regular')
#            self.adjust_template('#SBATCH -C cpu|nvme', '#SBATCH -C haswell')
#            self.adjust_template('#SBATCH -A proj36','#SBATCH -L SCRATCH')
#            Path_append ='export PATH="/global/common/software/m2043/AIBS_Opt/software/x86_64/bin:$PATH"'
#            self.adjust_template('source activate %s'%self.conda_env, Path_append,
#                                      add = True)
#            # HDF5 file locking
#            HDF5_cmd = "export HDF5_USE_FILE_LOCKING=FALSE"
#            self.adjust_template('source activate %s'%self.conda_env,
#                                  HDF5_cmd,add = True)
#
#        elif 'hpc-login' in machine:
#            self.adjust_template('cp -r $SCRIPT_REPO/modfiles $STAGE_DIR/',
#                             'cp -r $SCRIPT_REPO/x86_64 $STAGE_DIR/',add = True)
#            self.adjust_template('nrnivmodl modfiles/',
#                     '\techo "Loading compiled modfiles"',partial_match = True)
#            
#
#        elif self.submit_cmd == 'bash':
#            self.adjust_template('RES=$(sh batch_job.sh)', 'sh batch_job.sh')
#            self.adjust_template('echo ${RES##* }', '',partial_match = True)


    def run_job(self):

        os.system('chmod +x %s'%self.script_name)
        os.system('bash %s'%self.script_name)


class test_JobModule(JobModule):

    def __init__(self,script_name,job_config_path):
        
        super(test_JobModule,self).__init__(script_name)
        
        self.job_config_path = job_config_path

    def script_generator(self,chain_job='chain_job.sh',**kwargs):
        job_config = utility.load_json(self.job_config_path)
        stage_jobconfig = job_config['stage_jobconfig']
        highlevel_job_props = job_config['highlevel_jobconfig']
        for option,option_val in dryrun_config.items():
            if option in stage_jobconfig:
                stage_jobconfig[option] = option_val
        stage_jobconfig['ipyp_optim'] = False
        if isinstance(stage_jobconfig.get('seed'),list):
            stage_jobconfig['seed'] = 1
        utility.save_json(self.job_config_path,job_config)
        
        
        testjob_string = '#!/bin/bash\n'
        testjob_string +='set -ex\n'
        testjob_string +='source activate %s\n'%highlevel_job_props['conda_env']
        testjob_string += 'python %s --job_config %s\n'%\
                (stage_jobconfig['main_script'],self.job_config_path)
        testjob_string += 'python %s --input_json %s\n'\
                        %(stage_jobconfig['analysis_script'],self.job_config_path)
        
        if 'next_stage_job_config' in kwargs.keys():
            if bool(kwargs['next_stage_job_config']):
                testjob_string += 'bash %s\n'%chain_job
        with open(self.script_name, "w") as shell_script:
            shell_script.write(testjob_string)

    def run_job(self):
        utility.create_filepath(self.cp_file)
        os.system('chmod +x %s'%self.script_name)
        process = Popen(['sh', '%s'%self.script_name], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        logger.debug(stderr)




class Slurm_JobModule(JobModule):
    def __init__(self, script_template,job_config_path,script_name = 'batch_job.sh'):

        super(Slurm_JobModule,self).__init__(script_name)
        self.job_config_path = job_config_path
        self.script_template = utility.locate_template_file(script_template)
        self.submit_cmd = 'sbatch'

    def script_generator(self):
        job_config = utility.load_json(self.job_config_path)
        with open(self.script_template,'r') as job_template:
            batchjob_string = job_template.read()

        batchjob_string = batchjob_string.replace('conda_env',self.conda_env)
        with open(self.script_name, "w") as batchjob_script:
            batchjob_script.write(batchjob_string)

        if 'cori' in self.machine:
            self.adjust_template('#SBATCH -p prod', '#SBATCH -q regular')
            self.adjust_template('#SBATCH -C cpu|nvme', '#SBATCH -C haswell')
            self.adjust_template('#SBATCH -A proj36','#SBATCH -L SCRATCH')
            self.adjust_template('#SBATCH -n 256', '#SBATCH -N 8')
            Path_append ='export PATH="/global/common/software/m2043/AIBS_Opt/software/x86_64/bin:$PATH"'
            self.adjust_template('source activate %s'%self.conda_env, Path_append,
                                  add = True)
            self.adjust_template('#SBATCH --mail-user','',partial_match = True)

    def submit_job(self):
        os.system('chmod +x %s'%self.script_name)
        process = Popen(['%s', '%s'%(self.submit_cmd,self.script_name)],
                         stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        logger.debug(stderr)


class PBS_JobModule(JobModule):
    def __init__(self, script_template,job_config_path,script_name = 'batch_job.sh'):
        self.job_config_path = job_config_path
        super(PBS_JobModule,self).__init__(script_name)
        self.script_template = utility.locate_template_file(script_template)
        self.submit_cmd = 'qsub'


    def script_generator(self,chain_job='chain_job.sh',**kwargs):
        job_config = utility.load_json(self.job_config_path)
        stage_jobconfig = job_config['stage_jobconfig']
        
        highlevel_job_props = job_config['highlevel_jobconfig']
        
        if highlevel_job_props['dryrun'] or (kwargs.get('analysis') and \
                             not stage_jobconfig.get('run_hof_analysis')):
            for option,option_val in dryrun_config.items():
                if stage_jobconfig.get(option):
                    stage_jobconfig[option] = option_val
            
        utility.save_json(self.job_config_path,job_config)
        
        
        with open(self.script_template,'r') as job_template:
            batchjob_string = job_template.read()

        jobname = '%s.%s'%(os.path.basename(highlevel_job_props['job_dir']),
                       stage_jobconfig['stage_name'])
        seed_string = ''.join(['%s '%seed_ for seed_ in stage_jobconfig['seed']])
        
        # High level job config
        batchjob_string = batchjob_string.replace('conda_env',
                          highlevel_job_props['conda_env'])
        batchjob_string = batchjob_string.replace('jobname',jobname)
        
        # Stage Job config
        batchjob_string = batchjob_string.replace('jobscript_name',self.script_name)
        batchjob_string = batchjob_string.replace('jobmem',
                          stage_jobconfig['jobmem'])
        batchjob_string = batchjob_string.replace('ipyp_db',
                          stage_jobconfig['ipyp_db'])
        batchjob_string = batchjob_string.replace('qos',
                          job_config['stage_jobconfig']['qos'])
        batchjob_string = batchjob_string.replace('main_script',
                          stage_jobconfig['main_script'])
        batchjob_string = batchjob_string.replace('job_config_path',
                          self.job_config_path)
        batchjob_string = batchjob_string.replace('seed_list',seed_string)
        batchjob_string = batchjob_string.replace('analysis_script',
                          stage_jobconfig['analysis_script'])
        
        # Job config analysis vs optimization 
        batchjob_string = (batchjob_string.replace('jobtime',
          stage_jobconfig['jobtime_analysis'])
          if kwargs.get('analysis') else batchjob_string.replace('jobtime',
          stage_jobconfig['jobtime']))
        batchjob_string = (batchjob_string.replace('error_stream',
             stage_jobconfig['error_stream_analysis']) if kwargs.get('analysis')
            else batchjob_string.replace('error_stream',stage_jobconfig['error_stream']))
        batchjob_string = (batchjob_string.replace('output_stream',
             stage_jobconfig['output_stream_analysis']) if kwargs.get('analysis')
            else batchjob_string.replace('output_stream',stage_jobconfig['output_stream']))
        batchjob_string = (batchjob_string.replace('nnodes',
          str(stage_jobconfig['nnodes_analysis']))
          if kwargs.get('analysis') else batchjob_string.replace('nnodes',
          str(stage_jobconfig['nnodes'])))
        batchjob_string = (batchjob_string.replace('nprocs',
          str(stage_jobconfig['nprocs_analysis']))
          if kwargs.get('analysis') else batchjob_string.replace('nprocs',
          str(stage_jobconfig['nprocs'])))
        batchjob_string = (batchjob_string.replace('nengines',
          str(stage_jobconfig['nprocs_analysis']))
          if kwargs.get('analysis') else batchjob_string.replace('nengines',
          str(stage_jobconfig['nengines'])))
        
         
        if kwargs.get('analysis'):
            batchjob_string = re.sub('# Run[\S\s]*pids','',
                                   batchjob_string) 
            if not stage_jobconfig['run_hof_analysis']:
                batchjob_string = re.sub('# Configure[\S\s]*pids','',
                                   batchjob_string) 
         
        if 'next_stage_job_config' in kwargs.keys():
            if bool(kwargs['next_stage_job_config']):
                if (kwargs.get('analysis') and stage_jobconfig.get('ipyp_analysis'))\
                    or not kwargs.get('analysis'):
                        batchjob_string += 'bash %s\n'%chain_job
        
        if not kwargs.get('analysis') and stage_jobconfig.get('ipyp_analysis'):
            batchjob_string = re.sub('# Analyze[\S\s]*.json','qsub analyze_job.sh',
                                   batchjob_string) 
            
            
        with open(self.script_name, "w") as batchjob_script:
            batchjob_script.write(batchjob_string)

    def submit_job(self):

        os.system('chmod +x %s'%self.script_name)
        process = Popen(['%s', '%s'%(self.submit_cmd,self.script_name)],
                         stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        logger.debug(stderr)

class SGE_JobModule(JobModule):
    def __init__(self, script_template, machine,script_name = 'batch_job.sh',
                     conda_env='ateam_opt'):

        super(SGE_JobModule,self).__init__(machine,script_name)
        self.conda_env = conda_env
        self.script_template = utility.locate_template_file(script_template)
        self.submit_cmd = 'qsub'


    def script_generator(self):
        with open(self.script_template,'r') as job_template:
            batchjob_string = job_template.read()

        batchjob_string = batchjob_string.replace('conda_env',self.conda_env)
        with open(self.script_name, "w") as batchjob_script:
            batchjob_script.write(batchjob_string)

    def submit_job(self):

        os.system('chmod +x %s'%self.script_name)
        process = Popen(['%s', '%s'%(self.submit_cmd,self.script_name)],
                         stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        logger.debug(stderr)
