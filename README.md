[![Build Status](https://travis-ci.com/anirban6908/All-active-Workflow.svg?token=93Twb9jDYFzVNoM9gSjr&branch=master)](https://travis-ci.com/anirban6908/All-active-Workflow)

# All-active-Workflow
Creating the code base for All-active Model generation and analysis written on top of Bluepyopt

#### Launching optimization jobs
* AIBS hpc
```sh
$ launch_optimjob --help # Look at the options
$ launch_optimjob --cell_id xyz --ext_scripts /allen/aibs/mat/anin/software/All-active-Workflow/examples/optim_scripts 
$ submit_opt_jobs -f cell_data.csv -r /allen/aibs/mat/anin/software/All-active-Workflow/examples/optim_scripts/ -c ateam_opt -m 2 # Launching multiple jobs from a csv file
```
* NERSC Cori
```
```



