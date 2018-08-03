#!/bin/bash


PWD=$(pwd)
export PASS_DIR="Passive_Model"
export PASS_IH_DIR="Ih_Step"
export ALL_ACTIV_DIR="All_Active1"

export PASSIVE_REPO="/project/projectdirs/m2043/AIBS/ani/Human/PC/Passive_Repo"
export PASS_IH_REPO="/project/projectdirs/m2043/AIBS/ani/Human/PC/Passive_Ih_Repo"
export ALL_ACTIV_REPO="/project/projectdirs/m2043/AIBS/ani/Human/PC/All_Active1_Repo"

if [ ! -d "$PASS_DIR" ]; then
	
	python get_ephys_morphology_model.py && mkdir $PASS_DIR && mv cell_types $PASS_DIR/
	echo "Downloading ephys data, morphology and (possibly) model"
	model_filename=$(find . -maxdepth 2 -name "*fit.json")
	if [ ! -z "$model_filename" ]; then
		mv $model_filename $PASS_DIR/cell_types/fit_parameters.json && mv neuronal_model/modfiles $PASS_DIR && rm -rf neuronal_model
	fi
	
else
	echo "All data already downloded" 
fi

# Launch the passive optimization (Stage 0)

cp -r $PASSIVE_REPO/* $PASS_DIR/
cd $PASS_DIR
python starter_optim.py
echo "Launching Stage 0 Opimization"
RES0=$(sbatch start_haswell.sh)  # sbatch command goes here
STATUS_0=$(sacct -j ${RES##* } -o State| sed -n '3 p'| xargs) # get the status of the job 
if [[ $STATUS_0 = "COMPLETED" ]]; then
    echo "Stage 0 finished successfully" > Stage0_status.txt
else
    echo "Stage 0 did NOT finish successfully" > Stage0_status.txt
fi
python Optim_Main.py --checkpoint seed1.pkl --short_analyse
echo "Saving the Optimized parameters for the next stage"
rm -rf preprocessed/ 
cd ..


# Launch the passive+Ih optimization (Stage 1)

mkdir $PASS_IH_DIR
cp -r $PASS_DIR/cell_types/ $PASS_IH_DIR
rm -rf $PASS_DIR/cell_types
mv $PASS_DIR/fit_opt.json $PASS_IH_DIR/cell_types
cp -r $PASS_IH_REPO/* $PASS_IH_DIR/
cd $PASS_IH_DIR
python starter_optim.py
nrnivmodl modfiles/
echo "Launching Stage 1 Opimization"
RES1=$(sbatch start_haswell.sh)  # sbatch command goes here
STATUS_1=$(sacct -j ${RES1##* } -o State| sed -n '3 p'| xargs) # get the status of the job 
if [[ $STATUS_1 = "COMPLETED" ]]; then
    echo "Stage 1 finished successfully" > Stage1_status.txt
else
    echo "Stage 1 did NOT finish successfully" > Stage1_status.txt
fi
python Optim_Main.py --checkpoint seed1.pkl --short_analyse
echo "Saving the Optimized parameters for the next stage"
cd ..

# Launch the All-active optimization (Stage 2)

mkdir $ALL_ACTIV_DIR
cp -r $PASS_IH_DIR/cell_types/ $ALL_ACTIV_DIR
rm -rf $PASS_IH_DIR/cell_types
mv $PASS_IH_DIR/fit_opt.json $ALL_ACTIV_DIR/cell_types
cp -r $ALL_ACTIV_REPO/* $ALL_ACTIV_DIR/
cd $ALL_ACTIV_DIR
python starter_optim.py
nrnivmodl modfiles/
echo "Launching Stage 2 Opimization"
RES2=$(sbatch start_haswell.sh) && RES3=$(sbatch --dependency=afternotok:${RES2##* } restart_haswell.sh)  # sbatch command goes here
STATUS_3=$(sacct -j ${RES3##* } -o State| sed -n '3 p'| xargs) # get the status of the job 
if [[ $STATUS_3 = "COMPLETED" ]]; then
    echo "Stage 2 finished successfully" > Stage1_status.txt
    echo "This is the end"

else
    echo "Stage 2 did NOT finish successfully. Run restart_haswell. " > Stage1_status.txt
fi
cd ..