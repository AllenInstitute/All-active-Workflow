#!/bin/bash -l 
#SBATCH -q shared 
#SBATCH -t 0:30:00 
#SBATCH -J relaunch_Stage2 
#SBATCH -C haswell 
#SBATCH -n 1 


sbatch start_haswell.sh