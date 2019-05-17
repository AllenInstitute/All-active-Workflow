if [ $MAX_NGEN = 200 ]; then
    sed -i -e "s/MAX_NGEN=200/MAX_NGEN=400/g" batch_job.sh
    sed -i -e "s/config_file_basic.json/config_file.json/g" batch_job.sh
    sed -i -e "s/Optim_Main.py/Optim_Main_cont.py/g" batch_job.sh
    sleep 30
    qsub batch_job.sh
else
    mv ${CHECKPOINTS_DIR}/* checkpoints_final/
    qsub analyze_results.sh
fi
