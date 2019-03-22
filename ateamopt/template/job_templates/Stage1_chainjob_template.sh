#!/bin/sh

#SBATCH -q debug
#SBATCH -t 0:30:00
#SBATCH -L SCRATCH
#SBATCH -J launch_Stage_1
#SBATCH -C haswell
#SBATCH -N 1


set -ex


source activate ateam

JOBID_0=$(<Job_0.txt)
PARENT_DIR=$(<pwd.txt)
export STAGE_DIR=$PARENT_DIR/Stage1
export SCRIPT_REPO=$PARENT_DIR/Script_Repo

mkdir $STAGE_DIR
cp pwd.txt $STAGE_DIR/

STATUS_0=$(sacct -j ${JOBID_0} -o State| sed -n '3 p'| xargs) # get the status of the job
if [[ $STATUS_0 = "COMPLETED" ]]; then
    echo "Stage 0 finished successfully" > Stage0_status.txt
else
    echo "Stage 0 did NOT finish successfully" > Stage0_status.txt
fi
python Optim_Main.py --checkpoint checkpoints_backup/seed1.pkl --short_analyse
python Optim_Main.py -vv --checkpoint checkpoints_backup/seed1.pkl  --responses resp_opt.txt   --analyse

echo "Saving the Optimized parameters for the next stage"
#rm -rf preprocessed/


# Launch the passive+Ih optimization (Stage 1)

cp -r cell_types/ $STAGE_DIR/
cp cell_id.txt $STAGE_DIR/
mv fit_opt.json $STAGE_DIR/cell_types/
if [ -d "peri_model" ]; then mv peri_model/ $STAGE_DIR/; fi
cp -r $PASS_IH_REPO/* $STAGE_DIR/
if [ -f nersc_queue.txt ]; then cp nersc_queue.txt $STAGE_DIR/ ; fi
cd $STAGE_DIR
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
if [ -f nersc_queue.txt ]; then sed -i -e "s/regular/premium/g" start_haswell.sh ; fi
sed -i -e "s/Stage_2/$LAUNCH_JOBNAME/g" launch_stage2.slurm
echo "Launching Stage 1 Opimization"
RES_1=$(sbatch start_haswell.sh)  # sbatch command goes here
echo ${RES_1##* } > Job_1.txt


