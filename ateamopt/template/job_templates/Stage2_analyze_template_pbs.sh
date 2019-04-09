#!/bin/sh

#PBS -q celltypes
#PBS -l walltime=8:00:00
#PBS -l nodes=4:ppn=10
#PBS -l mem=100g
#PBS -N analyze_Stage2
#PBS -r n
#PBS -m n

cd $PBS_O_WORKDIR

set -ex

source activate conda_env

PWD=$(pwd)
PARENT_DIR=$(<pwd.txt)
CELL_ID=$(<cell_id.txt)

LOGS=$PWD/logs
mkdir -p $LOGS

export IPYTHONDIR=${PWD}/.ipython
export IPYTHON_PROFILE=benchmark.$PBS_JOBID

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
srun --output="${LOGS}/engine_%j_%2t.out" ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

python analysis_stage2.py -vv --checkpoint  checkpoints_final  --ipyparallel
# python save_time.py
# python perisomatic_sim.py

rm -rf $IPYTHONDIR
