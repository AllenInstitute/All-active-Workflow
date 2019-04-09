#!/bin/sh

#SBATCH -p prod
#SBATCH -t 8:00:00
#SBATCH -n 40
#SBATCH -C cpu|nvme
#SBATCH -A proj36
#SBATCH --mail-type=ALL
#SBATCH -J analyze_Stage2

set -ex

source activate conda_env

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=benchmark.${SLURM_JOBID}

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

python analysis_stage2.py -vv --checkpoint  checkpoints_final  --ipyparallel
# python save_time.py
# python perisomatic_sim.py

rm -rf $IPYTHONDIR
