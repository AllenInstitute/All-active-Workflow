#!/bin/bash

csv_file=$1
script_repo=$2
conda=$3
max_jobs=$4
index=0
pids=""

while IFS=, read -r cell_id
do
    echo "Launching optimization for $cell_id"
    nohup launch_optimjob --cell_id $cell_id --ext_scripts $script_repo --conda_env $conda > out$index.log 2> err$index.log &
    pids+="$! "
    index=$(($index+1))
    # echo $index
    # if [[ $index = 10 ]]; then
    #     break
    # fi
done < $csv_file

wait $pids

rm -rf err*.log out*.log
echo "# of jobs submitted = $index"

