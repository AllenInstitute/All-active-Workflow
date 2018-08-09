#!/bin/bash


PARENT_DIR=$(pwd)
export PASS_DIR="Passive_Model"

export PASSIVE_REPO="/project/projectdirs/m2043/AIBS/ani/Mouse/PC/Passive_Repo"

if [ ! -d "$PASS_DIR" ]; then
	
	python get_ephys_morphology_model.py && mkdir $PASS_DIR && mv cell_types $PASS_DIR/
	echo "Downloading ephys data, morphology and (possibly) model"
	model_filename=$(find . -maxdepth 2 -name "*fit.json")
	if [ ! -z "$model_filename" ]; then
		mv $model_filename $PASS_DIR/cell_types/fit_parameters.json && mv neuronal_model/modfiles $PASS_DIR && rm -rf neuronal_model
	fi
	
else
	echo "All data already downloaded" 
fi

# Launch the passive optimization (Stage 0)

cp -r $PASSIVE_REPO/* $PASS_DIR/
cd $PASS_DIR
python starter_optim.py
echo "Launching Stage 0 Opimization"
RES_0=$(sbatch start_haswell.sh)  # sbatch command goes here
echo ${RES_0##* } > Job_0.txt
echo $PARENT_DIR > pwd.txt

