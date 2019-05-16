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
                echo "$CELL_ID does not exist"
            fi
        fi
    done
