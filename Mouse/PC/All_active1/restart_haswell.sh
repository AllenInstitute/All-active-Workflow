#!/bin/bash

#SBATCH -q premium
#SBATCH -N 10
#SBATCH -t 12:00:00
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH --mail-user=anirban.nandi@wustl.edu
#SBATCH --mail-type=ALL
#SBATCH -J Stage2



set -e
set -x

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

OFFSPRING_SIZE=640
MAX_NGEN=100

export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=benchmark.${SLURM_JOBID}

ipcontroller --init --ip='*' --sqlitedb --profile=${IPYTHON_PROFILE} &
sleep 10
srun -n 320 -c 2 --cpu_bind=cores --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

CHECKPOINTS_DIR="checkpoints"
mkdir -p ${CHECKPOINTS_DIR}

pids=""
for seed in 1; do
    python Optim_Main.py             \
        -vv                                \
        --compile                          \
        --offspring_size=${OFFSPRING_SIZE} \
        --max_ngen=${MAX_NGEN}             \
        --seed=${seed}                     \
        --ipyparallel                      \
        --continu                         \
        --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" &
    pids+="$! "
done

wait $pids
