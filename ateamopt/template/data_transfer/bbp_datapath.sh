#!/bin/bash

set -ex

bbp_human=/gpfs/bbp.cscs.ch/home/anirban/Optimizations/Human_Cells/running_cells
bbp_mouse=/gpfs/bbp.cscs.ch/home/anirban/Optimizations/Mouse_Cells/running_cells

# Human cells
cd $bbp_human
for dir_human in */
    do
        if [ -f $dir_human/All_Active1/time_metrics_*.csv ]; then
            (cd $dir_human && echo $(pwd))
        fi
    done


# Mouse cells
cd $bbp_mouse
for dir_mouse in */
  do
      if [ -f $dir_mouse/Stage2/time_metrics_*.csv ]; then
          (cd $dir_mouse && echo $(pwd))
      fi
  done
