#!/bin/bash

#SBATCH -q premium
#SBATCH -N 8
#SBATCH -t 4:00:00
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH --mail-user=anirban.nandi@wustl.edu
#SBATCH --mail-type=ALL
#SBATCH -J Stage0



set -e
set -x

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

OFFSPRING_SIZE=512
MAX_NGEN=50

export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=benchmark.${SLURM_JOBID}

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun -n 256 --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
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
        --start                         \
        --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" &
    pids+="$! "
done

wait $pids

# Launch the passive+Ih optimization (Stage 1)
sh launch_stage1.sh