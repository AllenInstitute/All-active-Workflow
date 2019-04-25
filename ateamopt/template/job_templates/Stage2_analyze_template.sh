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

PARENT_DIR=$(<pwd.txt)
CELL_ID=$(<cell_id.txt)

export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=slurm.${SLURM_JOBID}

ipcontroller --init --ip='*' --nodb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

python analysis_stage2.py -vv --cp_dir  checkpoints_final  --ipyparallel

# Cleaning up

rm -rf $IPYTHONDIR $LOGS
rm -rf checkpoints_backup



