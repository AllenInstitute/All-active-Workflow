import os
import logging
from subprocess import Popen, PIPE
from ateamopt.utils import utility


logger = logging.getLogger(__name__)


class JobModule(object):

    def __init__(self, machine,script_name = 'Jobscript.sh'):

        self.machine = machine
        self.script_name = script_name


    def adjust_template(self,match_line, replace_line, add = False):
        with open(self.script_name, "r") as in_file:
            buf = in_file.readlines()

        with open(self.script_name, "w") as out_file:
            for line in buf:
                if line == "%s\n"%match_line:
                    if add:
                        line +=  "%s\n"%replace_line
                    else:
                        line = "%s\n"%replace_line
                out_file.write(line)

class ChainSubJob(JobModule):

    def __init__(self, script_template, machine,script_name = 'chain_job.sh',
                     conda_env='ateam_opt'):

        super(ChainSubJob,self).__init__(machine,script_name)

        self.conda_env = conda_env
        self.script_template = utility.locate_template_file(script_template)


    def script_generator(self):
        with open(self.script_template,'r') as job_template:
            subjob_string = job_template.read()

        subjob_string = subjob_string.replace('conda_env',self.conda_env)

        if 'hpc-login' in self.machine:
            submit_cmd = 'qsub'
            subjob_string = subjob_string.replace('submit_cmd',submit_cmd)
        elif any(substring in self.machine for substring in ['cori', 'bbp']):
            submit_cmd = 'sbatch'
            subjob_string = subjob_string.replace('submit_cmd',submit_cmd)
        else:
            submit_cmd = 'sh'
            subjob_string = subjob_string.replace('submit_cmd',submit_cmd)

        with open(self.script_name, "w") as chainsubjob_script:
            chainsubjob_script.write(subjob_string)

        # Adjusting the job based on machine

        if 'cori' in self.machine:
            self.adjust_template('#SBATCH -p prod', '#SBATCH -q regular')
            self.adjust_template('#SBATCH -C cpu|nvme', '#SBATCH -C haswell')
            self.adjust_template('#SBATCH -A proj36','#SBATCH -L SCRATCH')
            self.adjust_template('#SBATCH -n 256', '#SBATCH -N 8')
            Path_append ='export PATH="/global/common/software/m2043/AIBS_Opt/software/x86_64/bin:$PATH"'
            self.adjust_template('source activate %s'%self.conda_env, Path_append,
                                      add = True)
            # HDF5 file locking
            HDF5_cmd = "export HDF5_USE_FILE_LOCKING=FALSE"
            self.adjust_template('source activate %s'%self.conda_env,
                                  HDF5_cmd,add = True)

        elif 'hpc-login' in self.machine:
            self.adjust_template('nrnivmodl modfiles/',
                                 'echo "Loading compiled modfiles"')



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
        super(test_JobModule,self).__init__(machine,script_name)
        self.cp_file = cp_file
        self.offspring = offspring
        self.max_ngen = max_ngen
        self.optim_script = optim_script
        self.parallel = parallel
        self.job_status = job_status


    def script_generator(self):
        utility.create_filepath(self.cp_file)
        testjob_string = 'python %s -vv --checkpoint %s'%(self.optim_script,
                                                          self.cp_file)

        if self.parallel:
           testjob_string += ' --ipyparallel'

        testjob_string += ' --offspring_size=%s --max_ngen=%s --%s\n'%(self.offspring,
                                                    self.max_ngen,self.job_status)
        testjob_string += 'sh chain_job.sh\n'
        with open(self.script_name, "w") as shell_script:
            shell_script.write(testjob_string)

    def run_job(self):
        utility.create_filepath(self.cp_file)
        os.system('chmod +x %s'%self.script_name)
        process = Popen(['sh', '%s'%self.script_name], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        print(stderr)




class Slurm_JobModule(JobModule):
    def __init__(self, script_template, machine,script_name = 'batch_job.sh'):

        super(Slurm_JobModule,self).__init__(machine,script_name)

        if 'cori' in self.machine:
            self.conda_env = 'ateam'
        elif 'bbp' in self.machine:
            self.conda_env = 'CompNeuro'

        self.script_template = utility.locate_template_file(script_template)
        self.submit_verb = 'sbatch'

    def script_generator(self):
        with open(self.script_template,'r') as job_template:
            batchjob_string = job_template.read()

        batchjob_string = batchjob_string.replace('conda_env',self.conda_env)
        with open(self.script_name, "w") as batchjob_script:
            batchjob_script.write(batchjob_string)

        self.adjust_template('#SBATCH -p prod', '#SBATCH -q regular')
        self.adjust_template('#SBATCH -C cpu|nvme', '#SBATCH -C haswell')
        self.adjust_template('#SBATCH -A proj36','#SBATCH -L SCRATCH')
        self.adjust_template('#SBATCH -n 256', '#SBATCH -N 8')
        Path_append ='export PATH="/global/common/software/m2043/AIBS_Opt/software/x86_64/bin:$PATH"'
        self.adjust_template('source activate %s'%self.conda_env, Path_append,
                              add = True)

    def submit_job(self):

        os.system('chmod +x %s'%self.script_name)
        process = Popen(['%s', '%s'%(self.submit_verb,self.script_name)],
                         stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        print(stderr)


class PBS_JobModule(JobModule):
    def __init__(self, script_template, machine,script_name = 'batch_job.sh',
                     conda_env='ateam_opt'):

        super(PBS_JobModule,self).__init__(machine,script_name)
        self.conda_env = conda_env
        self.script_template = utility.locate_template_file(script_template)
        self.submit_verb = 'qsub'


    def script_generator(self):
        with open(self.script_template,'r') as job_template:
            batchjob_string = job_template.read()

        batchjob_string = batchjob_string.replace('ateam_opt',self.conda_env)
        with open(self.script_name, "w") as batchjob_script:
            batchjob_script.write(batchjob_string)

    def submit_job(self):

        os.system('chmod +x %s'%self.script_name)
        process = Popen(['%s', '%s'%(self.submit_verb,self.script_name)],
                         stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        print(stderr)
