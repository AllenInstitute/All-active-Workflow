#!/bin/sh

#PBS -q celltypes
#PBS -l walltime=4:00:00
#PBS -l nodes=16:ppn=16
#PBS -l mem=100g
#PBS -N Stage1
#PBS -r n
#PBS -m n

cd $PBS_O_WORKDIR

set -ex

source activate conda_env

OFFSPRING_SIZE=512
MAX_NGEN=50
seed=1


PWD=$(pwd)
export IPYTHONDIR=$PWD/.ipython
file $IPYTHONDIR
export IPYTHON_PROFILE=pbs.$PBS_JOBID

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 30
file $IPYTHONDIR/$IPYTHON_PROFILE
mpiexec -n 256 ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 30

CHECKPOINTS_DIR="checkpoints"
mkdir -p ${CHECKPOINTS_DIR}

python Optim_Main.py             \
    -vv                                \
    --offspring_size=${OFFSPRING_SIZE} \
    --max_ngen=${MAX_NGEN}             \
    --seed=${seed}                     \
    --ipyparallel                      \
    --start                        \
    --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" &

pid=$!
wait $pid

# Launch the Stage 2 optimization
sh chain_job.sh
