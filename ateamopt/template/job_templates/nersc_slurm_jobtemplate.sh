#!/bin/bash

#SBATCH -q qos
#SBATCH -t jobtime
#SBATCH -N nnodes
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH --mail-type=ALL
#SBATCH -J jobname
#SBATCH --signal=B:USR1@60

run_dependent_script() {
func="$1" ; shift
for sig ; do
trap "$func $sig" "$sig"
done
}

# trap function to relaunch the optimization
func_trap() {
sbatch jobscript_name
}

# submit launch script upon signal USR1
run_dependent_script func_trap USR1

set -ex

source activate conda_env
export PATH="/global/common/software/m2043/AIBS_Opt/software/x86_64/bin:$PATH"


PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

# Configure ipython profile
export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=slurm.${SLURM_JOBID}

# Start ipcontroller and engines
ipcontroller --init --ip='*' --ipyp_db --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun -n nengines --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10


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
