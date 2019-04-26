#!/bin/sh

#PBS -q celltypes
#PBS -l walltime=24:00:00
#PBS -l nodes=16:ppn=16
#PBS -l mem=100g
#PBS -N Stage2
#PBS -e /dev/null
#PBS -r n
#PBS -m n

cd $PBS_O_WORKDIR

set -ex

source activate conda_env

# Relaunch batch job if not finished
qsub -W depend=afternotok:$PBS_JOBID batch_job.sh

OFFSPRING_SIZE=512
MAX_NGEN=200
timeout=300


PWD=$(pwd)
export IPYTHONDIR=$PWD/.ipython
ipython profile create
file $IPYTHONDIR
export IPYTHON_PROFILE=pbs.$PBS_JOBID

ipcontroller --init --ip='*' --nodb --ping=30000 --profile=${IPYTHON_PROFILE} &
sleep 30
file $IPYTHONDIR/$IPYTHON_PROFILE
mpiexec -n 256 ipengine --timeout=3000 --profile=${IPYTHON_PROFILE} &
sleep 30

CHECKPOINTS_DIR="checkpoints"
CHECKPOINTS_BACKUP="checkpoints_backup"

mkdir -p $CHECKPOINTS_DIR
mkdir -p $CHECKPOINTS_BACKUP
mkdir -p checkpoints_final


# Check the job status : Start or continue
if [ "$(ls -A $CHECKPOINTS_DIR)" ]; then
    JOB_STATUS=continu
else
    JOB_STATUS=start
fi

pids=""
for seed in {1..4}; do
    python Optim_Main.py             \
        -vv                                \
        --offspring_size=${OFFSPRING_SIZE} \
        --max_ngen=${MAX_NGEN}             \
        --seed=${seed}                     \
        --ipyparallel                      \
        --$JOB_STATUS                      \
        --timeout=$timeout                 \
        --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" \
        --cp_backup "${CHECKPOINTS_BACKUP}/seed${seed}.pkl" &
    pids+="$! "
done

wait $pids

# If job finishes in time analyze result
mv ${CHECKPOINTS_DIR}/* checkpoints_final/

# check if the job with 4th seed is finished
# if [[ $seed = 4 ]]; then
qsub analyze_results.sh
# else
#     seed_new=$(($seed+1))
#     sed -i -e "s/seed=$seed/seed=$seed_new/g" batch_job.sh
#     qsub batch_job.sh
# fi
