#!/bin/sh

#SBATCH -p prod
#SBATCH -t 12:00:00
#SBATCH -n 256
#SBATCH -C cpu|nvme
#SBATCH -A proj36
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anin@alleninstitute.org
#SBATCH -J Stage2
#SBATCH --signal=B:USR1@120

run_dependent_script() {
func="$1" ; shift
for sig ; do
trap "$func $sig" "$sig"
done
}

# trap function to relaunch the optimization
func_trap() {
sbatch batch_job.sh
}

# submit launch script upon signal USR1
run_dependent_script func_trap USR1

set -ex

source activate conda_env

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

OFFSPRING_SIZE=512
MAX_NGEN=200
seed=1

export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=benchmark.${SLURM_JOBID}

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun -n 256 --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

CHECKPOINTS_DIR="checkpoints"
mkdir -p ${CHECKPOINTS_DIR}


# Check the job status : Start or continue
if [ "$(ls -A $CHECKPOINTS_DIR)" ]; then
    JOB_STATUS=continu
else
    JOB_STATUS=start
fi


python Optim_Main.py             \
    -vv                                \
    --offspring_size=${OFFSPRING_SIZE} \
    --max_ngen=${MAX_NGEN}             \
    --seed=${seed}                     \
    --ipyparallel                      \
    --$JOB_STATUS                         \
    --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" &

pid=$!
wait $pid


# If job finishes in time analyze result
mv ${CHECKPOINTS_DIR}/seed${seed}.pkl checkpoints_final/


# check if the job with 4th seed is finished

if [[ $seed = 4 ]]; then
    sbatch analyze_results.sh
else
    seed_new=$(($seed+1))
    sed -i -e "s/seed=$seed/seed=$seed_new/g" batch_job.sh
    sbatch batch_job.sh
fi