#!/bin/sh

set -ex

nersc_host=ani@cori.nersc.gov

scp -i /local1/nersc_ssh/mynersc nersc_datapath.sh $nersc_host:~/
ssh -i /local1/nersc_ssh/mynersc $nersc_host "sh nersc_datapath.sh" > nersc_log
ssh -i /local1/nersc_ssh/mynersc $nersc_host "rm nersc_datapath.sh"


for line in $(<nersc_log);
    do
        nersc_path=$nersc_host:$line
        CELL_ID=${nersc_path##*/}
        mkdir -p $CELL_ID
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/*.pdf $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/fitted_params $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/config $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/analysis_params/hof_features_all.pkl $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/analysis_params/hof_obj*.pkl $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/analysis_params/score_list_train.pkl $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/analysis_params/seed_indices.pkl $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/Validation_Responses/exp* $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/Validation_Responses/fitness* $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/Validation_Responses/Features* $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/cell_metadata* $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/morph_stats* $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/time_metrics* $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/Validation_Responses/fI* $CELL_ID/
        scp -i /local1/nersc_ssh/mynersc -r $nersc_path/Stage2/Validation_Responses/AP* $CELL_ID/
    done
