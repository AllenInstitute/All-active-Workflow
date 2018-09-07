#!/bin/bash -l 
#SBATCH -q shared 
#SBATCH -t 5:00:00 
#SBATCH -J analyze_Stage2 
#SBATCH -C haswell 
#SBATCH -n 1 

sh short_analyse.sh
python calculate_validation_resp.py
