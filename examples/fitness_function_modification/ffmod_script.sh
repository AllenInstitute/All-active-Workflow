if [ $MAX_NGEN = 100 ]; then
    sed -i -e "s/MAX_NGEN=100/MAX_NGEN=200/g" batch_job.sh
    sed -i -e "s/config_file_basic.json/config_file.json/g" batch_job.sh
    qsub batch_job.sh
else
    mv ${CHECKPOINTS_DIR}/* checkpoints_final/
    qsub analyze_results.sh
fi
