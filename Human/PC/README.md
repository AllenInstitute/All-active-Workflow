### Directory structure for the scripts to work:

* Put the .nwb and .swc file in a directory named 'cell_types'
* Then run starter_optim.py

> Make sure you have the following directory structure:

     .      
     ├── ...
     ├── cell_metadata.json
     ├── Passive_Model/Ih_step/All_Active1 (Stage 0, Stage 1, Stage 2)
     │               ├── ...
     │               ├── cell_types                             
     │               │   ├── .nwb                               # stimulus response file
     │               │   ├── .swc                               # morphology file
     │               │   └── fit_parameters.json                # parameters of a previously fitted model (if one exists)
     │               │── config ── cell_id                                 
     │               │              ├── parameters.json         # parameters to be fitted
     │               │              ├── mechanism.json          # mechanisms to be added in each section 
     │               │              ├── features.json           # features on which the model is fitted
     │               │              └── protocols.json          # protocols used (stimulus and recording)
     │               │
     │               │── starter_optim.py
     │               │── get_features.py
     │               │── Optim_Main.py
     │               │── model_helper.py
     │               │── evaluator_helper.py
     │               │── feature_set.json
     │               │── *bounds.json
     │               │── modfiles/ 
     │               ├── ...
     │
     ├── ...
