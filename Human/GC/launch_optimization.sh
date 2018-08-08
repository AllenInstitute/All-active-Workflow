#!/bin/bash


PARENT_DIR=$(pwd)
CELL_ID=${PARENT_DIR##*/}
export PASS_DIR="Passive_Model"
export PASSIVE_REPO="/project/projectdirs/m2043/AIBS/ani/Human/GC/Passive_Repo"

echo "All data already downloaded" 
mkdir $PASS_DIR

# Launch the passive optimization (Stage 0)

cp -r cell_types/ $PASS_DIR/
cp -r $PASSIVE_REPO/* $PASS_DIR/
cd $PASS_DIR
python starter_optim.py
echo $CELL_ID > cell_id.txt
echo "Launching Stage 0 Opimization"
RES_0=$(sbatch start_haswell.sh)  # sbatch command goes here
echo ${RES_0##* } > Job_0.txt
echo $PARENT_DIR > pwd.txt
