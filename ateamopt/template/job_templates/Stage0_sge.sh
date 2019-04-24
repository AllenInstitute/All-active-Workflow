#!/bin/sh

#$ -pe orte 256
#$ -N Stage0
#$ -cwd
#$ -V
#$ -j y
#$ -m bes
#$ -M anin@alleninstitute.org
#$ -S /bin/sh

set -ex

#source activate conda_env

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

python Optim_Main.py             \
    -vv                                \
    --offspring_size=${OFFSPRING_SIZE} \
    --max_ngen=${MAX_NGEN}             \
    --seed=${seed}                     \
    --ipyparallel                      \
    --start                            \
    --timeout=$timeout                 \
    --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" &

pid=$!
wait $pid

# Launch the Stage 1 optimization
sh chain_job.sh
