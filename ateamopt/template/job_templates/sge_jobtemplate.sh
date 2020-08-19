#!/bin/bash

#$ -pe orte nengines
#$ -N jobname
#$ -cwd
#$ -V
#$ -j y
#$ -m bes
#$ -M email
#$ -S /bin/bash

set -ex

source activate conda_env

# Configure ipython profile
PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

export IPYTHONDIR=$PWD/.ipython
export IPYTHON_PROFILE=sge.$JOB_ID

# Start ipcontroller and engines
ipcontroller --init --ip='*' --ipyparallel_db --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
mpiexec -np nengines --output-filename $LOGS/engine ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

# Run Optimization
pids=""
for seed in seed_list; do
    python optimizer_script             \
        --seed=${seed}                     \
        --job_config job_config_path &
    pids+="$! "
done

wait $pids

# Analyze results
python analysis_script --input_json stage_job_config.json