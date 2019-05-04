#!/bin/sh

#PBS -q celltypes
#PBS -l walltime=5:00:00
#PBS -l nodes=16:ppn=16
#PBS -l mem=100g
#PBS -N Stage1
#PBS -e /dev/null
#PBS -o /dev/null
#PBS -r n
#PBS -m bea

cd $PBS_O_WORKDIR

set -ex

source activate conda_env

PWD=$(pwd)
LOGS=$PWD/logs
mkdir -p $LOGS

OFFSPRING_SIZE=512
MAX_NGEN=50
seed=1
timeout=300

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
mkdir -p ${CHECKPOINTS_DIR}

# Check the job status : Start or continue
if [ "$(ls -A $CHECKPOINTS_DIR)" ]; then
    JOB_STATUS=continu
else
    JOB_STATUS=start
fi

python Optim_Main.py             \
    -vv                                \
    --offspring_size=${OFFSPRING_SIZE} \
    --max_ngen=${MAX_NGEN}             \
    --seed=${seed}                     \
    --ipyparallel                      \
    --$JOB_STATUS                      \
    --timeout=$timeout                 \
    --checkpoint "${CHECKPOINTS_DIR}/seed${seed}.pkl" &

pid=$!
wait $pid

# Launch the Stage 2 optimization
sh chain_job.sh
