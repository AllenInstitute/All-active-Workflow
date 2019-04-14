#!/bin/sh

#SBATCH -p prod
#SBATCH -t 30:00
#SBATCH -n 1
#SBATCH -C cpu|nvme
#SBATCH -A proj36
#SBATCH -J launch_Stage_1

set -ex

# Activate conda environment

source activate conda_env

JOBID_0=$(<Job_0.txt)
PARENT_DIR=$(<pwd.txt)
CELL_ID=$(<cell_id.txt)
IPYTHONDIR=.ipython

export STAGE_DIR=$PARENT_DIR/Stage1
export SCRIPT_REPO=$PARENT_DIR/Script_Repo

mkdir -p $STAGE_DIR

# Check if the batch job was completed

STATUS_0=$(sacct -j ${JOBID_0} -o State| sed -n '3 p'| xargs) # get the status of the job
if [[ $STATUS_0 = "COMPLETED" ]]; then
    echo "Stage 0 finished successfully" > Stage0_status.txt
else
    echo "Stage 0 did NOT finish successfully" > Stage0_status.txt
fi

# Run analysis

python analysis_stage0.py
echo "Saving the Optimized parameters for the next stage"

# Cleaning up large files and Moving data

# rm -rf $IPYTHONDIR
cp -r cell_types $STAGE_DIR/
mv fitted_params/fit_opt.json $STAGE_DIR/cell_types/
if [ -d "peri_model" ]; then mv peri_model/ $STAGE_DIR/; fi
if [ -f qos.txt ]; then cp qos.txt $STAGE_DIR/ ; fi
cp $SCRIPT_REPO/{prepare_stage1_run.py,analysis_stage1.py} $STAGE_DIR/
cp -r $SCRIPT_REPO/modfiles $STAGE_DIR/

# Run scripts to prepare for the batch-job

cd $STAGE_DIR
python prepare_stage1_run.py conda_env
if [ -d modfiles ]; then nrnivmodl modfiles/ ; fi # Compile mechanisms
STAGE="_STAGE1"
STAGE_NEXT="_STAGE2"
JOBNAME=$CELL_ID$STAGE
LAUNCH_JOBNAME=$CELL_ID$STAGE_NEXT
sed -i -e "s/Stage1/$JOBNAME/g" batch_job.sh
if [ -f qos.txt ]; then
    queue=$(<qos.txt)
    sed -i -e "s/regular/$queue/g" batch_job.sh # Specific to Cori
fi
sed -i -e "s/Stage_2/$LAUNCH_JOBNAME/g" chain_job.sh
echo $PARENT_DIR > pwd.txt
echo $CELL_ID > cell_id.txt

# Launch the batch job (Stage 1)
echo "Launching Stage 1 Opimization"
RES=$(submit_cmd batch_job.sh)
echo ${RES##* } > Job_1.txt


