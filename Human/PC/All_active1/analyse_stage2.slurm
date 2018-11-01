#!/bin/bash -l 
#SBATCH -q regular 
#SBATCH -t 5:00:00 
#SBATCH -L SCRATCH
#SBATCH --mail-user=anirban.nandi@wustl.edu
#SBATCH --mail-type=ALL
#SBATCH -J analyze_Stage2 
#SBATCH -C haswell 
#SBATCH -n 1 

sh short_analyse.sh
python calculate_validation_resp.py
python save_time.py
sh analyse.sh
python get_validation_metrics.py

