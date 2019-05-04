#!/bin/sh

#PBS -q celltypes
#PBS -l walltime=10:00
#PBS -l nodes=1:ppn=32
#PBS -l mem=100g
#PBS -N sens_analysis
#PBS -j oe
#PBS -r n
#PBS -m n

cd $PBS_O_WORKDIR
set -ex
source activate ateam_opt

python sa_example.py
