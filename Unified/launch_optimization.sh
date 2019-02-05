#!/bin/bash

PARENT_DIR=$(pwd)
CELL_ID=${PARENT_DIR##*/}
export PASS_DIR="Passive_Model"
export PASSIVE_REPO="/project/projectdirs/m2043/AIBS/ani/Unified/Passive_Repo"

source activate ateam
export HDF5_USE_FILE_LOCKING=FALSE

if [ ! -d "$PASS_DIR" ]; then
	
	mkdir $PASS_DIR && mv cell_types $PASS_DIR/
	all_active_model_file=$(find neuronal_model -maxdepth 1 -name "*fit.json")
	peri_model_file=$(find peri_model -maxdepth 1 -name "*fit.json")
	if [ ! -z "$all_active_model_file" ]; then
		mv $all_active_model_file $PASS_DIR/cell_types/fit_parameters.json && rm -rf neuronal_model
	fi
	if [ ! -z "$peri_model_file" ]; then
		mkdir -p $PASS_DIR/peri_model && mv $peri_model_file $PASS_DIR/peri_model/ && rm -rf peri_model
	fi
	
else
	echo "All data already moved for Stage 0" 
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
sed -i -e "s/Stage0/$JOBNAME/g" start_haswell.sh
sed -i -e "s/Stage_1/$LAUNCH_JOBNAME/g" launch_stage1.slurm
echo "Launching Stage 0 Opimization"
RES_0=$(sbatch start_haswell.sh)  # sbatch command goes here
echo ${RES_0##* } > Job_0.txt
echo $PARENT_DIR > pwd.txt
echo $CELL_ID > cell_id.txt

