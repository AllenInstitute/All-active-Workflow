#!/bin/sh

#SBATCH -p prod
#SBATCH -t 12:00:00
#SBATCH -n 256
#SBATCH -C cpu|nvme
#SBATCH -A proj36
#SBATCH --mail-type=ALL
#SBATCH --mail-user=anin@alleninstitute.org
#SBATCH -J Stage2
#SBATCH --signal=B:USR1@60

run_dependent_script() {
func="$1" ; shift
for sig ; do
trap "$func $sig" "$sig"
done
}

# trap function to relaunch the optimization
func_trap() {
scancel ${SLURM_JOBID}
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
# timeout=900
seed=1

export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=slurm.${SLURM_JOBID}

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun -n 256 --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

CHECKPOINTS_DIR="checkpoints"
CHECKPOINTS_BACKUP="checkpoints_backup"
mkdir -p $CHECKPOINTS_DIR
mkdir -p $CHECKPOINTS_BACKUP
mkdir -p checkpoints_final


# Check the job status : Start or continue
if [ "$(ls -A $CHECKPOINTS_DIR)" ]; then
    JOB_STATUS=continu
else
    JOB_STATUS=start
fi

pids=""
for seed in {1..4}; do
    python Optim_Main.py             \
        -vv                                \
        --offspring_size=${OFFSPRING_SIZE} \
        --max_ngen=${MAX_NGEN}             \
        --seed=${seed}                     \
        --ipyparallel                      \
        --$JOB_STATUS                         \
        # --timeout=$timeout              \
        --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" \
        --cp_backup "${CHECKPOINTS_BACKUP}/seed${seed}.pkl" &
    pids+="$! "
done

wait $pids


# If job finishes in time analyze result
mv ${CHECKPOINTS_DIR}/* checkpoints_final/


# check if the job with 4th seed is finished

# if [[ $seed = 4 ]]; then
sbatch analyze_results.sh
# else
#     seed_new=$(($seed+1))
#     sed -i -e "s/seed=$seed/seed=$seed_new/g" batch_job.sh
#     sbatch batch_job.sh
# fi
