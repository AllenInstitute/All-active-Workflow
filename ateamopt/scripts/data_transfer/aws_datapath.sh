#!/bin/bash

set -ex

aws_mouse=/aibs_data

# To be used with starcluster
# Mouse cells
cd $aws_mouse
for dir_mouse in */
  do
      if [ -f $dir_mouse/Stage2/time_metrics_*.csv ]; then
          (cd $dir_mouse && echo $(pwd))
      fi
  done
