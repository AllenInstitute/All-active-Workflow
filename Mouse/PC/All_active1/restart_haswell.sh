#!/bin/bash

#SBATCH -q premium
#SBATCH -N 16
#SBATCH -t 12:00:00
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH --mail-user=anirban.nandi@wustl.edu
#SBATCH --mail-type=ALL
#SBATCH -J Stage2
#SBATCH --signal=B:USR1@60

run_dependent_script() { 
func="$1" ; shift 
for sig ; do 
trap "$func $sig" "$sig" 
done 
} 

# trap function to launch the passive+Ih optimization (Stage 2)
func_trap() {
mv checkpoints checkpoints_old
mv checkpoints_backup checkpoints 
sbatch restart_haswell.sh
#echo Launching Stage2 through xfer job because of timeout
} 

#submit launch script upon signal USR1 
run_dependent_script func_trap USR1 

set -e
set -x

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

OFFSPRING_SIZE=1024
MAX_NGEN=500

export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=benchmark.${SLURM_JOBID}

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun -n 512 --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
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

# If job finishes in time analyze result
mv checkpoints checkpoints_old
mv checkpoints_backup checkpoints 
sh analyse.sh
