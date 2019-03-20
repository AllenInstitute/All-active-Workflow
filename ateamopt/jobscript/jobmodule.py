import os
import logging
from subprocess import Popen, PIPE
from ateamopt.utils import utility


logger = logging.getLogger(__name__)


class JobModule(object):
    
    def __init__(self, machine,script_name = 'Jobscript.sh'):
        
        self.machine = machine
        self.script_name = script_name
        

class ChainSubJob(JobModule):
    
    def __init__(self, script_template, machine,script_name = 'chain_job.sh',
                     conda_env='ateam_opt'):
        
        super().__init__(machine,script_name)
        self.conda_env = conda_env
        self.script_template = utility.locate_template_file(script_template)
        
        
    def script_generator(self):
        with open(self.script_template,'r') as job_template:
            subjob_string = job_template.read()
        
        if self.conda_env != 'ateam_opt':
            subjob_string = subjob_string.replace('ateam_opt',self.conda_env)
            
        if 'hpc-login' in self.machine:
            submit_cmd = 'qsub'
            subjob_string = subjob_string.replace('submit_cmd',submit_cmd)
            
        chainsubjob_script = open(self.script_name, "w")
        chainsubjob_script.write(subjob_string)
     
    def run_job(self):
        
        os.system('chmod +x %s'%self.script_name)
        process = Popen(['sh', '%s'%self.script_name], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        print(stderr)        
            

class test_JobModule(JobModule):
    
    def __init__(self,machine,script_name,cp_file,
                 offspring, max_ngen,
                 optim_script= 'Optim_Main.py',
                 parallel = False,
                 job_status = 'start'):
        super().__init__(machine,script_name)
        self.cp_file = cp_file
        self.offspring = offspring
        self.max_ngen = max_ngen
        self.optim_script = optim_script
        self.parallel = parallel
        self.job_status = job_status
        
        
    def script_generator(self):
        testjob_string = 'python %s -vv --checkpoint %s'%(self.optim_script, 
                                                          self.cp_file)
            
        if self.parallel:
           testjob_string += ' --ipyparallel' 
            
        testjob_string += ' --offspring_size=%s --max_ngen=%s --%s'%(self.offspring,
                                                    self.max_ngen,self.job_status)
        
        shell_script = open(self.script_name, "w")
        shell_script.write(testjob_string)
     
    def run_job(self):
        utility.create_filepath(self.cp_file)
        os.system('chmod +x %s'%self.script_name)
        process = Popen(['sh', '%s'%self.script_name], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        print(stderr)
        
        
        

class Slurm_JobModule(JobModule):
    pass

class PBS_JobModule(JobModule):
    def __init__(self, script_template, machine,script_name = 'batch_job.sh',
                     conda_env='ateam_opt'):
    
        super().__init__(machine,script_name)
        self.conda_env = conda_env
        self.script_template = utility.locate_template_file(script_template)
        self.submit_verb = 'qsub'
        
        
    def script_generator(self):
        with open(self.script_template,'r') as job_template:
            batchjob_string = job_template.read()
        
        if self.conda_env != 'ateam_opt':
            batchjob_string = batchjob_string.replace('ateam_opt',self.conda_env)
            
         
        batchjob_script = open(self.script_name, "w")
        batchjob_script.write(batchjob_string)
        
    def submit_job(self):
        
        os.system('chmod +x %s'%self.script_name)
        process = Popen(['%s', '%s'%(self.submit_verb,self.script_name)],
                         stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        print(stderr)        