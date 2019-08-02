#!/bin/sh

set -ex

# Run this script from where the metadata is saved.

# Mouse cells
s3_mouse_bucket=s3://aibs.test.ani/Musmusculus/
aws s3 ls $s3_mouse_bucket --profile wasabi > aws_log

for line in $(<aws_log);
    do
        if [ $line != "PRE" ]; then
            CELL_ID=$line
            echo $CELL_ID
            if [ ! -d $CELL_ID ] ; then
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/ $CELL_ID --exclude "*" --include "*.pdf" --recursive --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/fitted_params ${CELL_ID}fitted_params --recursive --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/config $CELL_ID/config --recursive --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/config_file.json $CELL_ID/ --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/analysis_params/hof_features_all.pkl ${CELL_ID} --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/analysis_params/ ${CELL_ID} --exclude "*" --include "hof_obj*" --recursive --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/analysis_params/score_list_train.pkl ${CELL_ID} --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/analysis_params/seed_indices.pkl ${CELL_ID} --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/Validation_Responses/ ${CELL_ID} --exclude "*" --include "*.pkl" --recursive --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/Validation_Responses/ ${CELL_ID} --exclude "*" --include "exp*.csv" --recursive --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/Validation_Responses/ ${CELL_ID} --exclude "*" --include "fitness_metrics*" --recursive --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID} $CELL_ID --exclude "*" --include "cell_metadata*.json" --recursive --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID} $CELL_ID --exclude "*" --include "morph_stats*.json" --recursive --profile wasabi
                aws s3 cp ${s3_mouse_bucket}${CELL_ID}Stage2/ $CELL_ID --exclude "*" --include "time_metrics*" --recursive --profile wasabi
            fi
        fi
    done
