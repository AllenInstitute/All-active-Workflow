#!/bin/bash

PARENT_DIR=$(pwd)
CELL_ID=${PARENT_DIR##*/}
export PASS_DIR="Passive_Model"
export PASSIVE_REPO="/gpfs/bbp.cscs.ch/home/anirban/Unified/Passive_Repo"

if [ ! -d "$PASS_DIR" ]; then
	
	mkdir $PASS_DIR && mv cell_types $PASS_DIR/
	model_filename=$(find . -maxdepth 2 -name "*fit.json")
	if [ ! -z "$model_filename" ]; then
		mv $model_filename $PASS_DIR/cell_types/fit_parameters.json && rm -rf neuronal_model
	fi
	
else
	echo "All data already downloaded" 
fi

# Launch the passive optimization (Stage 0)

cp -r $PASSIVE_REPO/* $PASS_DIR/
cd $PASS_DIR
python set_features_passive.py
python set_params_passive.py
python starter_optim.py
STAGE="_STAGE0"
STAGE_NEXT="_STAGE1"
JOBNAME=$CELL_ID$STAGE
LAUNCH_JOBNAME=$CELL_ID$STAGE_NEXT
sed -i -e "s/Stage0/$JOBNAME/g" start_bbp.sh
sed -i -e "s/Stage_1/$LAUNCH_JOBNAME/g" launch_stage1_bbp.slurm
echo "Launching Stage 0 Opimization"
RES_0=$(sbatch start_bbp.sh)  # sbatch command goes here
echo ${RES_0##* } > Job_0.txt
echo $PARENT_DIR > pwd.txt
echo $CELL_ID > cell_id.txt

