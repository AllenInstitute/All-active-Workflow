#!/bin/sh

#PBS -q celltypes
#PBS -l walltime=4:00:00
#PBS -l nodes=16:ppn=16
#PBS -l mem=10g
#PBS -N Stage0
#PBS -r n
#PBS -m n

cd $PBS_O_WORKDIR

set -ex

source activate ateam_opt

OFFSPRING_SIZE=512
MAX_NGEN=50

PWD=$(pwd)
export IPYTHONDIR=$PWD/.ipython
export IPYTHON_PROFILE=pbs.$PBS_JOBID

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
mpiexec -n 256 ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

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

# Launch the Stage 1 optimization 
#qsub launch_stage1_hpc.pbs
