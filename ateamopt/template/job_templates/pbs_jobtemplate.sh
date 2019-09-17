#!/bin/bash

#PBS -q qos
#PBS -l walltime=jobtime
#PBS -l nodes=nnodes:ppn=nprocs
#PBS -l mem=jobmem
#PBS -N jobname
#PBS -e error_stream
#PBS -o output_stream
#PBS -r n
#PBS -m bea

cd $PBS_O_WORKDIR

set -e

source activate conda_env

# Relaunch batch job if not finished
qsub -W depend=afternotok:$PBS_JOBID jobscript_name


# Configure ipython profile
PWD=$(pwd)
export IPYTHONDIR=$PWD/.ipython
ipython profile create
file $IPYTHONDIR
export IPYTHON_PROFILE=pbs.$PBS_JOBID

# Start ipcontroller and engines
ipcontroller --init --ip='*' --ipyp_db --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 30
file $IPYTHONDIR/profile_$IPYTHON_PROFILE
mpiexec -n nengines ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 30


# Run Optimization
pids=""
for seed in seed_list; do
    python main_script             \
        --seed ${seed}             \
        --input_json job_config_path &
    pids+="$! "
done

wait $pids


# Analyze results
python analysis_script --input_json stage_job_config.json
