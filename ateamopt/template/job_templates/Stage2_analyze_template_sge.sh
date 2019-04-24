#!/bin/bash

#$ -pe orte 40
#$ -N analyze_Stage2
#$ -cwd
#$ -V
#$ -j y
#$ -m bes
#$ -M anin@alleninstitute.org
#$ -S /bin/bash

set -ex

#source activate conda_env

PWD=$(pwd)
PARENT_DIR=$(<pwd.txt)
CELL_ID=$(<cell_id.txt)

export IPYTHONDIR=$PWD/.ipython
export IPYTHON_PROFILE=sge.$JOB_ID

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 10
mpiexec -np 40 ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 10

python analysis_stage2.py -vv --cp_dir  checkpoints_final  --ipyparallel

# Cleaning up

rm -rf $IPYTHONDIR $LOGS
rm -rf checkpoints_backup



