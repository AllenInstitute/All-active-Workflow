#!/bin/bash

#$ -pe orte 128
#$ -N Stage0
#$ -cwd
#$ -V
#$ -j y
#$ -m bes
#$ -M anin@alleninstitute.org
#$ -S /bin/bash

set -ex

source activate conda_env

OFFSPRING_SIZE=512
MAX_NGEN=50
seed=1
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
    --$JOB_STATUS                      \
    --timeout=$timeout                 \
    --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" &

pid=$!
wait $pid

# Launch the Stage 1 optimization
bash chain_job.sh
