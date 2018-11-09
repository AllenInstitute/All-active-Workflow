#!/bin/bash

#SBATCH -q regular
#SBATCH -N 8
#SBATCH -t 36:00:00
#SBATCH -C haswell
#SBATCH -L SCRATCH
#SBATCH --mail-user=anirban.nandi@wustl.edu
#SBATCH --mail-type=ALL
#SBATCH -J Stage2
#SBATCH --signal=B:USR1@120

run_dependent_script() { 
func="$1" ; shift 
for sig ; do 
trap "$func $sig" "$sig" 
done 
} 

# trap function to launch the passive+Ih optimization (Stage 2)
func_trap() {
new_cp=checkpoints.${SLURM_JOBID}
mv checkpoints $new_cp
mv checkpoints_backup checkpoints 
sbatch restart_batchjob_stage2.slurm
} 

#submit launch script upon signal USR1 
run_dependent_script func_trap USR1 

set -e
set -x

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

OFFSPRING_SIZE=512
MAX_NGEN=200

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


# If job finishes in time analyze result

mv ${CHECKPOINTS_DIR}/seed${seed}.pkl checkpoints_final/

#new_cp=checkpoints.${SLURM_JOBID}
#mv checkpoints $new_cp
#mv checkpoints_backup checkpoints

# check if the job with 4th seed is finished

if [[ $seed = 4 ]]; then
    sbatch analyse_stage2.slurm
else
    seed_new=$(($seed+1))
    sed -i -e "s/seed in $seed/seed in $seed_new/g" start_haswell.sh 
    sed -i -e "s/seed in $seed/seed in $seed_new/g" restart_haswell.sh
    sbatch start_batchjob_stage2.slurm
fi
