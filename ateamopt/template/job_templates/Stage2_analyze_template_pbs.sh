#!/bin/sh

#PBS -q celltypes
#PBS -l walltime=8:00:00
#PBS -l nodes=4:ppn=10
#PBS -l mem=100g
#PBS -N analyze_Stage2
#PBS -e /dev/null
#PBS -o /dev/null
#PBS -r n
#PBS -m bea

cd $PBS_O_WORKDIR

set -ex

source activate conda_env

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

PARENT_DIR=$(<pwd.txt)
CELL_ID=$(<cell_id.txt)

export IPYTHONDIR=$PWD/.ipython
ipython profile create
file $IPYTHONDIR
export IPYTHON_PROFILE=pbs.$PBS_JOBID

ipcontroller --init --ip='*' --nodb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 30
file $IPYTHONDIR/$IPYTHON_PROFILE
mpiexec -n 40 ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 30

python analysis_stage2.py -vv --cp_dir  checkpoints_final  --ipyparallel
# pid="$! "
# wait $pid

# Cleaning up

# rm -rf $LOGS
# rm -rf checkpoints_backup



