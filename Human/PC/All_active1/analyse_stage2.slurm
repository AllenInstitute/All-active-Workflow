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
<<<<<<< HEAD
python perisomatic_sim.py
=======
<<<<<<< HEAD
python perisomatic_sim.py
=======
>>>>>>> b673469fd583adba569777e00d7e8f6f81d54529
>>>>>>> ba68a140446b6740bc3d3e024fb7a7e94724bfe1
python get_validation_metrics.py

