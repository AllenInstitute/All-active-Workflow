#!/bin/bash

#$ -pe orte 128
#$ -N Stage2
#$ -cwd
#$ -V
#$ -j y
#$ -m bes
#$ -M anin@alleninstitute.org
#$ -S /bin/bash

set -ex

source activate conda_env

# Relaunch batch job if not finished
#qsub -W depend=afternotok:$JOB_ID batch_job.sh

OFFSPRING_SIZE=512
MAX_NGEN=200
timeout=300


PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

export IPYTHONDIR=$PWD/.ipython
export IPYTHON_PROFILE=sge.$JOB_ID

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
mpiexec -np 128 --output-filename $LOGS/engine ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
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
        --$JOB_STATUS                      \
        --timeout=$timeout                 \
        --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" \
        --cp_backup "${CHECKPOINTS_BACKUP}/seed${seed}.pkl" &
    pids+="$! "
done

wait $pids

# If job finishes in time analyze result
mv ${CHECKPOINTS_DIR}/* checkpoints_final/

qsub analyze_results.sh

