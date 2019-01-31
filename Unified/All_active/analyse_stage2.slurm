#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -t 4:00:00 
#SBATCH -L SCRATCH
#SBATCH --mail-user=anirban.nandi@wustl.edu
#SBATCH --mail-type=ALL
#SBATCH -J analyze_Stage2 
#SBATCH -C haswell 
#SBATCH -n 1 


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

python Optim_Main.py -vv --checkpoint  checkpoints_final/seed1.pkl   --responses resp_opt.txt   --response_release resp_release.txt  --short_analyse  --ipyparallel
python Optim_Main.py -vv --responses resp_opt.txt   --response_release resp_release.txt --analyse
python save_time.py
python perisomatic_sim.py 
#python get_validation_metrics.py
