#!/bin/bash

set -ex

nersc_human=/global/project/projectdirs/m2043/AIBS/ani/Optimizations/Human_Cells/running_cells
nersc_mouse=/global/project/projectdirs/m2043/AIBS/ani/Optimizations/Mouse_Cells/running_cells

#echo CELL_ID,SPECIES,STATUS,MACHINE
# Human cells
if [ -d $nersc_human ]; then
    cd $nersc_human
    for dir_human in */
           do
                   if [ -f $dir_human/cell_metadata*.json ]; then
                           (cd $dir_human && echo $(pwd))
                   fi
           done
fi

# Mouse cells
if [ -d $nersc_mouse ]; then
    cd $nersc_mouse
    for dir_mouse in */
           do
                   if [ -f $dir_mouse/Stage2/time_metrics_*.csv  ]; then
                           (cd $dir_mouse && echo $(pwd))
                   fi
           done
fi
