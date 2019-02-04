#!/bin/bash -l 
#SBATCH -q debug 
#SBATCH -t 0:10:00 
#SBATCH -L SCRATCH
#SBATCH -J relaunch_Stage2 
#SBATCH -C haswell 
#SBATCH -N 1 


sbatch restart_haswell.sh