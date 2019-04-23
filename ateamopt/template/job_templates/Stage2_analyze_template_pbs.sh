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

export IPYTHONDIR=$PWD/.ipython
file $IPYTHONDIR
export IPYTHON_PROFILE=pbs.$PBS_JOBID

ipcontroller --init --ip='*' --sqlitedb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 30
file $IPYTHONDIR/$IPYTHON_PROFILE
mpiexec -n 40 ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 30

python analysis_stage2.py -vv --cp_dir  checkpoints_final  --ipyparallel

# Cleaning up

rm -rf $IPYTHONDIR $LOGS
rm -rf checkpoints_backup



