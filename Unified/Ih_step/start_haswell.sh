#!/bin/bash

#SBATCH -q regular
#SBATCH -N 8
#SBATCH -t 4:00:00
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH --mail-type=ALL
#SBATCH -J Stage1
#SBATCH --signal=B:USR1@60

run_dependent_script() { 
func="$1" ; shift 
for sig ; do 
trap "$func $sig" "$sig" 
done 
} 

# trap function to launch the passive+Ih optimization (Stage 2)
func_trap() { 
sbatch launch_stage2.slurm
} 

#submit launch script upon signal USR1 
run_dependent_script func_trap USR1 


set -e
set -x

source activate ateam
export LD_LIBRARY_PATH="/global/common/software/m2043/AIBS_Opt/software/x86_64/lib:$LD_LIBRARY_PATH"
export PATH="/global/common/software/m2043/AIBS_Opt/software/x86_64/bin:$PATH"

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


# Launch the passive+Ih optimization (Stage 1) once the job finishes successfully
sbatch launch_stage2.slurm
