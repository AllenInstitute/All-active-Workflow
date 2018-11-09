#!/bin/bash

#SBATCH -p prod 
#SBATCH -t 5:00:00 
#SBATCH -n 1
#SBATCH -C cpu|nvme
#SBATCH -A proj36
#SBATCH --mail-user=anirban.nandi@wustl.edu
#SBATCH --mail-type=ALL
#SBATCH -J analyze_Stage2 

source activate CompNeuro
sh short_analyse.sh
python calculate_validation_resp.py
python save_time.py
sh analyse.sh
python perisomatic_sim.py
python get_validation_metrics.py
