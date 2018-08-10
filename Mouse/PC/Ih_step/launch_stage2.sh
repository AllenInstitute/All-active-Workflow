#!/bin/bash

export ALL_ACTIV_REPO="/project/projectdirs/m2043/AIBS/ani/Mouse/PC/All_Active1_Repo"

JOBID_1=$(<Job_1.txt)
PARENT_DIR=$(<pwd.txt)

export ALL_ACTIV_DIR=$PARENT_DIR/All_Active1
mkdir $ALL_ACTIV_DIR


STATUS_1=$(sacct -j ${JOBID_1} -o State| sed -n '3 p'| xargs) # get the status of the job 
if [[ $STATUS_1 = "COMPLETED" ]]; then
    echo "Stage 1 finished successfully" > Stage1_status.txt
else
    echo "Stage 1 did NOT finish successfully" > Stage1_status.txt
fi
python Optim_Main.py --checkpoint checkpoints/seed1.pkl --short_analyse
echo "Saving the Optimized parameters for the next stage"
#rm -rf preprocessed/ 


# Launch the All-active optimization (Stage 2)

cp -r cell_types/ $ALL_ACTIV_DIR/
#rm -rf cell_types
mv fit_opt.json $ALL_ACTIV_DIR/cell_types/
cp -r $ALL_ACTIV_REPO/* $ALL_ACTIV_DIR/
cd $ALL_ACTIV_DIR
python starter_optim.py
nrnivmodl modfiles/
echo "Launching Stage 2 Opimization"
RES_2=$(sbatch start_haswell.sh) && RES_3=$(sbatch --dependency=afternotok:${RES_2##* } restart_haswell.sh)  # sbatch command goes here
