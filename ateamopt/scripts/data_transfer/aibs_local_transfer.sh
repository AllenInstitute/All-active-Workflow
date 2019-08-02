#!/bin/sh

set -ex

# Can be run from any location on CentOS
aibs_datapath='/allen/programs/celltypes/workgroups/humancolumn_ephysmodeling/anin/Optimizations_HPC/Mouse'
aibs_runpath='/allen/aibs/mat/anin/hpc_trials'
aibs_metricpath='/allen/aibs/mat/ateam_shared/Mouse_Model_Fit_Metrics'
cell_complete_path=""

cd $aibs_runpath
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
            cp -r Stage2/config_file.json $target_path/
            cp -r Stage2/analysis_params/hof_features_all.pkl $target_path/
            cp -r Stage2/analysis_params/hof_obj*.pkl $target_path/
            cp -r Stage2/analysis_params/score_list_train.pkl $target_path/
            cp -r Stage2/analysis_params/seed_indices.pkl $target_path/
            cp -r Stage2/Validation_Responses/exp*.csv $target_path/
            cp -r Stage2/Validation_Responses/fitness* $target_path/
            cp -r Stage2/Validation_Responses/*.pkl $target_path/
            cp -r cell_metadata* $target_path/
            cp -r morph_stats* $target_path/
            cp -r Stage2/time_metrics* $target_path/

            cell_complete_path+="$path "

            cd $aibs_runpath
        fi
    done


for complete_cell in $cell_complete_path
    do
        CELL_ID=${complete_cell##*/}
        if [ ! -d $aibs_datapath/$CELL_ID ]; then
            mv $complete_cell $aibs_datapath
            sleep 5
        else
            echo "$CELL_ID is optimized twice"
        fi
    done
