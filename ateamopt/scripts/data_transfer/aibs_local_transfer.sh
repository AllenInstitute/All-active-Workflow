#!/bin/sh

set -ex

aibs_datapath='/allen/programs/celltypes/workgroups/humancolumn_ephysmodeling/anin/Optimizations_HPC/Mouse'
aibs_metricpath='/allen/aibs/mat/ateam_shared/Mouse_Model_Fit_Metrics'

cd $aibs_datapath
for dir_mouse in */
    do
        if [ -f ${dir_mouse}Stage2/time_metrics_*.csv ]; then
            cd $dir_mouse
            path=$(pwd)
            CELL_ID=${path##*/}
            echo $CELL_ID

            target_path=$aibs_metricpath/$CELL_ID
            mkdir -p $target_path
            cp -r Stage2/*.pdf $target_path/
            cp -r Stage2/fitted_params $target_path/
            cp -r Stage2/config $target_path/
            cp -r Stage2/Validation_Responses/exp* $target_path/
            cp -r Stage2/Validation_Responses/fitness* $target_path/
            cp -r cell_metadata* $target_path/
            cp -r morph_stats* $target_path/
            cp -r Stage2/time_metrics* $target_path/
            cp -r Stage2/Validation_Responses/fI* $target_path/
            cp -r Stage2/Validation_Responses/AP* $target_path/

            cd $aibs_datapath
        fi
    done


