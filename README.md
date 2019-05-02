[![Build Status](https://travis-ci.com/anirban6908/All-active-Workflow.svg?token=93Twb9jDYFzVNoM9gSjr&branch=master)](https://travis-ci.com/anirban6908/All-active-Workflow)

# All-active-Workflow
Creating the code base for All-active Model generation and analysis written on top of Bluepyopt

#### Launching optimization jobs
* AIBS hpc
```sh
$ source activate ateam_opt
$ launch_optimjob --help # Look at the options
$ launch_optimjob --cell_id xyz --ext_scripts /allen/aibs/mat/anin/software/All-active-Workflow/examples/optim_scripts 
$  script_repo='/allen/aibs/mat/anin/software/All-active-Workflow/examples/optim_scripts'
$ submit_opt_jobs -f cell_data.csv -r $script_repo -c ateam_opt -m 2 # Launching multiple jobs from a csv file
```
* NERSC Cori - using shared conda environment
```sh
$ source activate ateam
$ script_repo='/global/homes/a/ani/shared_software/All-active-Workflow/examples/optim_scripts'
$ launch_optimjob --cell_id xyz --conda_env ateam --ext_scripts $script_repo 
$ submit_opt_jobs -f cell_data.csv -r $script_repo -c ateam -m 2 # Launching multiple jobs from a csv file
```
* NERSC Cori - using docker image
```sh
$ 
``` 



