#!/bin/sh

#PBS -q celltypes
#PBS -l walltime=30:00
#PBS -l nodes=10:ppn=10
#PBS -l mem=100g
#PBS -N FI_Data_Parallel
#PBS -j oe
#PBS -r n
#PBS -m n

cd $PBS_O_WORKDIR
set -ex
source activate ateam_opt

PWD=$(pwd)
export IPYTHONDIR=$PWD/.ipython
file $IPYTHONDIR
export IPYTHON_PROFILE=pbs.$PBS_JOBID

ipcontroller --init --ip='*' --nodb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 30
file $IPYTHONDIR/$IPYTHON_PROFILE
mpiexec -n 100 ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 30

python avg_fi_parallel.py
