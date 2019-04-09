#!/bin/sh

set -ex

PARENT_DIR=$(pwd)
CELL_ID=${PARENT_DIR##*/}
export STAGE_DIR=$PARENT_DIR/Stage0
export SCRIPT_REPO=$PARENT_DIR/Script_Repo


# Activate conda environment

source activate conda_env

# Move the data

if [ ! -d "$STAGE_DIR" ]; then

	mkdir $STAGE_DIR && mv cell_types $STAGE_DIR/
	all_active_model_dir=neuronal_model
	peri_model_dir=peri_model

	if [ -d "$all_active_model_dir" ]; then
		all_active_model_file=$(find $all_active_model_dir -maxdepth 1 -name "*fit.json")
		mv $all_active_model_file $STAGE_DIR/cell_types/fit_parameters.json && rm -rf $all_active_model_dir
	fi
	if [ -d "$peri_model_dir" ]; then
		peri_model_file=$(find $peri_model_dir -maxdepth 1 -name "*fit.json")
		mkdir -p $STAGE_DIR/$peri_model_dir && mv $peri_model_file $STAGE_DIR/$peri_model_dir/ && rm -rf $peri_model_dir
	fi

else
	echo "All data already moved"
fi


# Run scripts to prepare for the batch-job

cp $SCRIPT_REPO/{prepare_stage0_run.py,analysis_stage0.py} $STAGE_DIR/
cp -r $SCRIPT_REPO/modfiles $STAGE_DIR/
if [ -f qos.txt ]; then cp qos.txt $STAGE_DIR/ ; fi # Specific to Cori
cd $STAGE_DIR

python prepare_stage0_run.py conda_env
STAGE="_STAGE0"
STAGE_NEXT="_STAGE1"
JOBNAME=$CELL_ID$STAGE
LAUNCH_JOBNAME=$CELL_ID$STAGE_NEXT
sed -i -e "s/Stage0/$JOBNAME/g" batch_job.sh
if [ -f qos.txt ]; then
	queue=$(<qos.txt)
	sed -i -e "s/regular/$queue/g" batch_job.sh # Specific to Cori
fi
sed -i -e "s/Stage_1/$LAUNCH_JOBNAME/g" chain_job.sh
if [ -d modfiles ]; then # Compile mechanisms
	nrnivmodl modfiles/
fi
echo $PARENT_DIR > pwd.txt
echo $CELL_ID > cell_id.txt

# Launch batch job

echo "Launching Stage 0 Opimization"
RES=$(submit_cmd batch_job.sh)
echo ${RES##* } > Job_0.txt

