#!/bin/sh

set -ex

bbp_host=anirban@bbpv1.epfl.ch

scp bbp_datapath.sh $bbp_host:~/
ssh $bbp_host "sh bbp_datapath.sh" > bbp_log
ssh $bbp_host "rm bbp_datapath.sh"

for line in $(<bbp_log);
    do
        bbp_path=$bbp_host:$line
        CELL_ID=${bbp_path##*/}
        mkdir -p $CELL_ID
        scp -r $bbp_path/Stage2/*.pdf $CELL_ID/
        scp -r $bbp_path/Stage2/fitted_params $CELL_ID/
        scp -r $bbp_path/Stage2/config $CELL_ID/
        scp -r $bbp_path/Stage2/analysis_params/hof_features_all.pkl $CELL_ID/
        scp -r $bbp_path/Stage2/analysis_params/hof_obj*.pkl $CELL_ID/
        scp -r $bbp_path/Stage2/analysis_params/score_list_train.pkl $CELL_ID/
        scp -r $bbp_path/Stage2/analysis_params/seed_indices.pkl $CELL_ID/
        scp -r $bbp_path/Stage2/Validation_Responses/exp*.csv $CELL_ID/
        scp -r $bbp_path/Stage2/Validation_Responses/fitness* $CELL_ID/
        scp -r $bbp_path/Stage2/Validation_Responses/*.pkl $CELL_ID/
        scp -r $bbp_path/cell_metadata* $CELL_ID/
        scp -r $bbp_path/morph_stats* $CELL_ID/
        scp -r $bbp_path/Stage2/time_metrics* $CELL_ID/
    done
