## Instructions for Launching Chained Jobs in NERSC machines

> Make sure you have the following directory structure:

     .      
     ├── ...
     ├── Parent Directory (Named by the cell id)
                     ├── ...
                     ├── launch_optimization.sh
                     ├── start_job.sh
                     ├── get_ephys_morphology_model.py
                     ├── save_cell_metadata.py
                     ├── cell_metadata.json
                     │
                     ├── Passive_Model (Stage 0) # sh launch_stage1.sh is added at the end of the batch script 
                     ├── Ih_Step (Stage 1)       # sh launch_stage2.sh is added at the end of the batch script
                     ├── All-Active1 (Stage 2)
    

### Step 1
> Create a directory with the cell-id and make sure launch_optimization.sh, get_ephys_morphology_model.py and save_cell_metadata.py are in this directory.

### Step 2
> Run save_cell_metadata.py to input cell-id, model-id (if available), dendrite type, species etc. This will create the cell_metadata.json. Just hit enter for fields you have no information on.

### Step 3
> Create a repo within the NERSC filesystem from which the necessary scripts will be copied into each cell directory. These files are common to each cell type (e.g., starter_optim.py, get_features.py). Look at the README.md at Human/PC for the list of necessary files.

### Step 4
> Now either execute start_job.sh or copy and paste the command at the terminal and hit enter and that should be it. This will daemonize(run in background) the process. Here first the cell data is downloaded using the allensdk api. It will also create the directory Passive Model and copy the necessary scripts from a repo location within the NERSC filesystem. After the necessary files are produced, this will launch the Stage 0 job and so on. 
