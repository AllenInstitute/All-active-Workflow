#!/bin/bash -l 
#SBATCH -q debug
#SBATCH -t 0:30:00 
#SBATCH -L SCRATCH
#SBATCH -J launch_Stage_1 
#SBATCH -C haswell 
#SBATCH -N 1 

export PASS_IH_REPO="/project/projectdirs/m2043/AIBS/ani/Human/PC/Passive_Ih_Repo"

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
python Optim_Main.py --checkpoint checkpoints_backup/seed1.pkl --short_analyse
echo "Saving the Optimized parameters for the next stage"
#rm -rf preprocessed/ 


# Launch the passive+Ih optimization (Stage 1)

cp -r cell_types/ $PASS_IH_DIR/
cp cell_id.txt $PASS_IH_DIR/
#rm -rf cell_types/
mv fit_opt.json $PASS_IH_DIR/cell_types/
cp -r $PASS_IH_REPO/* $PASS_IH_DIR/
cd $PASS_IH_DIR
python set_features_passive_and_Ih.py
python set_params_passive_and_Ih.py
python starter_optim.py
nrnivmodl modfiles/
STAGE="_STAGE1"
STAGE_NEXT="_STAGE2"
CELL_ID=$(<cell_id.txt)
JOBNAME=$CELL_ID$STAGE
LAUNCH_JOBNAME=$CELL_ID$STAGE_NEXT
sed -i -e "s/Stage1/$JOBNAME/g" start_haswell.sh
sed -i -e "s/Stage_2/$LAUNCH_JOBNAME/g" launch_stage2.slurm
echo "Launching Stage 1 Opimization"
RES_1=$(sbatch start_haswell.sh)  # sbatch command goes here
echo ${RES_1##* } > Job_1.txt
