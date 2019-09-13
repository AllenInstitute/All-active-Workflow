#!/bin/bash
set -e

# Activate conda environment
source activate conda_env

# Change to job directory
cd stage_jobdir

# Run scripts to prepare for the batch-job
python prepare_stagejob.py --input_json stage_job_config.json

# Copy modfiles
if [ -d modfiles_dir ]; then cp -r modfiles_dir stage_jobdir; fi

# Compile mechanisms
if [ -d modfiles_dir ]; then nrnivmodl modfiles_dir_abs; fi

# Copy compiled modfiles
if [ -d compiled_modfiles_dir ]; then cp -r compiled_modfiles_dir stage_jobdir; fi

# Launch batch job
echo "Launching Stage Opimization"
RES=$(submit_cmd batch_job.sh)
echo ${RES##* } > Job.txt

