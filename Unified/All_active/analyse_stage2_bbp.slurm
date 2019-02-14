#!/bin/bash

#SBATCH -p prod 
#SBATCH -t 5:00:00 
#SBATCH -n 40
#SBATCH -C cpu|nvme
#SBATCH -A proj36
#SBATCH --mail-user=anin@alleninstitute.org
#SBATCH --mail-type=ALL
#SBATCH -J analyze_Stage2 


set -e
set -x

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=benchmark.${SLURM_JOBID}

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

source activate CompNeuro
python Optim_Main.py -vv --checkpoint  checkpoints_final/seed1.pkl   --responses resp_opt.txt   --response_release resp_release.txt  --short_analyse  --ipyparallel
python Optim_Main.py -vv --responses resp_opt.txt   --response_release resp_release.txt --analyse
python save_time.py
python perisomatic_sim.py 
