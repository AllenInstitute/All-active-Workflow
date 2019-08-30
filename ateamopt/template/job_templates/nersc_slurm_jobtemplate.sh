#!/bin/bash

#SBATCH -q qos
#SBATCH -t jobtime
#SBATCH -n nengines
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
sbatch batch_job
}

# submit launch script upon signal USR1
run_dependent_script func_trap USR1

set -ex

source activate conda_env
export PATH="/global/common/software/m2043/AIBS_Opt/software/x86_64/bin:$PATH"


PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS
seed=1

# Configure ipython profile
export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=slurm.${SLURM_JOBID}

# Start ipcontroller and engines
ipcontroller --init --ip='*' --ipyp_db --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun -n nengines --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10


python optimizer_script             \
    --seed=${seed}                     \
    --job_config job_config_path &


# Check if the job for final seed is finished
if [[ $seed = max_seed ]]; then
    sbatch analyze_results.sh
else
    seed_new=$(($seed+1))
    sed -i -e "s/seed=$seed/seed=$seed_new/g" batch_job
    sbatch batch_job
fi
