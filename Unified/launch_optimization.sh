#!/bin/bash

set -ex

PARENT_DIR=$(pwd)
CELL_ID=${PARENT_DIR##*/}
export PASS_DIR="Passive_Model"
export PASSIVE_REPO="/project/projectdirs/m2043/AIBS/ani/Unified/Passive_Repo"

source activate ateam
export HDF5_USE_FILE_LOCKING=FALSE

if [ ! -d "$PASS_DIR" ]; then
	
	mkdir $PASS_DIR && mv cell_types $PASS_DIR/
	all_active_model_dir=neuronal_model
	peri_model_dir=peri_model
	
	if [ -d "$all_active_model_dir" ]; then
		all_active_model_file=$(find $all_active_model_dir -maxdepth 1 -name "*fit.json")
		mv $all_active_model_file $PASS_DIR/cell_types/fit_parameters.json && rm -rf $all_active_model_dir
	fi
	if [ -d "$peri_model_dir" ]; then
		peri_model_file=$(find $peri_model_dir -maxdepth 1 -name "*fit.json")
		mkdir -p $PASS_DIR/$peri_model_dir && mv $peri_model_file $PASS_DIR/$peri_model_dir/ && rm -rf $peri_model_dir
	fi
	
else
	echo "All data already moved for Stage 0" 
fi

# Launch the passive optimization (Stage 0)

cp -r $PASSIVE_REPO/* $PASS_DIR/
if [ -f nersc_queue.txt ]; then cp nersc_queue.txt $PASS_DIR/ ; fi
cd $PASS_DIR
python set_features_passive.py
python set_params_passive.py
python starter_optim.py
STAGE="_STAGE0"
STAGE_NEXT="_STAGE1"
JOBNAME=$CELL_ID$STAGE
LAUNCH_JOBNAME=$CELL_ID$STAGE_NEXT
sed -i -e "s/Stage0/$JOBNAME/g" start_haswell.sh
if [ -f nersc_queue.txt ]; then sed -i -e "s/regular/premium/g" start_haswell.sh ; fi
sed -i -e "s/Stage_1/$LAUNCH_JOBNAME/g" launch_stage1.slurm
echo "Launching Stage 0 Opimization"
RES_0=$(sbatch start_haswell.sh)  # sbatch command goes here
echo ${RES_0##* } > Job_0.txt
echo $PARENT_DIR > pwd.txt
echo $CELL_ID > cell_id.txt

