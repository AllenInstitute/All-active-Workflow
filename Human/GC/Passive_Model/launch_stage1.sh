#!/bin/bash

export PASS_IH_REPO="/project/projectdirs/m2043/AIBS/ani/Human/GC/Passive_Ih_Repo"

JOBID_0=$(<Job_0.txt)
PARENT_DIR=$(<pwd.txt)


export PASS_IH_DIR=$PARENT_DIR/Ih_Step
mkdir $PASS_IH_DIR
cp pwd.txt $PASS_IH_DIR/


STATUS_0=$(sacct -j ${JOBID_0} -o State| sed -n '3 p'| xargs) # get the status of the job 
if [[ $STATUS_0 = "COMPLETED" ]]; then
    echo "Stage 0 finished successfully" > Stage0_status.txt
else
    echo "Stage 0 did NOT finish successfully" > Stage0_status.txt
fi
python Optim_Main.py --checkpoint checkpoints/seed1.pkl --short_analyse
echo "Saving the Optimized parameters for the next stage"
rm -rf preprocessed/ 


# Launch the passive+Ih optimization (Stage 1)

cp -r cell_types/ $PASS_IH_DIR/
cp cell_id.txt $PASS_IH_DIR/
rm -rf cell_types/
mv fit_opt.json $PASS_IH_DIR/cell_types/
cp -r $PASS_IH_REPO/* $PASS_IH_DIR/
cd $PASS_IH_DIR
python starter_optim.py
nrnivmodl modfiles/
echo "Launching Stage 1 Opimization"
RES_1=$(sbatch start_haswell.sh)  # sbatch command goes here
echo ${RES_1##* } > Job_1.txt


