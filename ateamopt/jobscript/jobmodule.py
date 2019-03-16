import os
import logging
from subprocess import Popen, PIPE
from ateamopt.utils import utility


logger = logging.getLogger(__name__)


class JobModule(object):
    
    def __init__(self, machine,script_name = 'Jobscript.sh'):
        
        self.machine = machine
        self.script_name = script_name
        

class test_JobModule(JobModule):
    
    def __init__(self,machine,script_name,cp_file,
                 offspring, max_ngen,
                 optim_script= 'Optim_Main.py',
                 compile_stat = False,
                 parallel = False,
                 job_status = 'start'):
        super().__init__(machine,script_name)
        self.cp_file = cp_file
        self.offspring = offspring
        self.max_ngen = max_ngen
        self.optim_script = optim_script
        self.parallel = parallel
        self.compile_stat = compile_stat
        self.job_status = job_status
        
        
    def script_generator(self):
        testjob_string = 'python %s -vv --checkpoint %s'%(self.optim_script, 
                                                          self.cp_file)
        if self.compile_stat:
            testjob_string += ' --compile' 
            
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
    pass