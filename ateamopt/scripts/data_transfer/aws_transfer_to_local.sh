#!/bin/sh

set -ex

starcluster put ani_cluster aws_datapath.sh /root
starcluster sshmaster ani_cluster "sh aws_datapath.sh" > aws_log
starcluster sshmaster ani_cluster "rm aws_datapath.sh"


for line in $(<aws_log);
    do
        aws_path=$line
        CELL_ID=${aws_path##*/}
        mkdir -p $CELL_ID
        starcluster get ani_cluster $aws_path/Stage2/*.pdf $CELL_ID/
        starcluster get ani_cluster $aws_path/Stage2/fitted_params $CELL_ID/
        starcluster get ani_cluster $aws_path/Stage2/config $CELL_ID/
        starcluster get ani_cluster $aws_path/Stage2/Validation_Responses/exp* $CELL_ID/
        starcluster get ani_cluster $aws_path/Stage2/Validation_Responses/fitness* $CELL_ID/
        starcluster get ani_cluster $aws_path/cell_metadata* $CELL_ID/
        starcluster get ani_cluster $aws_path/morph_stats* $CELL_ID/
        starcluster get ani_cluster $aws_path/Stage2/time_metrics* $CELL_ID/
        starcluster get ani_cluster $aws_path/Stage2/Validation_Responses/fI* $CELL_ID/
        starcluster get ani_cluster $aws_path/Stage2/Validation_Responses/AP* $CELL_ID/
    done
