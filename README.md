[![Build Status](https://travis-ci.com/anirban6908/All-active-Workflow.svg?token=93Twb9jDYFzVNoM9gSjr&branch=master)](https://travis-ci.com/anirban6908/All-active-Workflow)

# All-active-Workflow
Creating the code base for All-active Model generation and analysis written on top of Bluepyopt

#### Launching optimization jobs
* AIBS hpc
```sh
$ source activate ateam_opt
$ launch_optimjob --help # Look at the options
$ launch_optimjob --cell_id xyz --ext_scripts /allen/aibs/mat/anin/software/All-active-Workflow/examples/optim_scripts 
$ launch_optimjob --cell_id xyz --ext_scripts /allen/aibs/mat/anin/software/All-active-Workflow/examples/optim_scripts --me_type ME_Exc_1 # launch jobs by passing me type
$ script_repo='/allen/aibs/mat/anin/software/All-active-Workflow/examples/optim_scripts'
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
# On your machine 
$ docker build -t=ateam_opt . # From inside the directory containing dockerfile
$ docker run -it ateam_opt # Run the container interactively
$ docker container ls # get the container id
$ docker container stop $container_id
$ docker login # login to dockerhub
$ docker tag ateam_opt anirban6908/aibs_software_ani:ateam_opt # tag the image for upload
$ docker push anirban6908/aibs_software_ani:ateam_opt # push the image

# On NERSC
$ shifterimg -v pull docker:anirban6908/aibs_software_ani:ateam_opt # pull the image (only needs to be done once)
$ salloc -N 1 -C haswell -q debug --image anirban6908/aibs_software_ani:ateam_opt --volume="/global/homes/a/ani/docker/bmtk_docker_sims:/app" # Run the image interactively
$ shifter /bin/bash # interactive bash shell in your shifter image
``` 
* BBP5
```sh
$ source activate CompNeuro
$ script_repo=~/All-active-Workflow/examples/optim_scripts
$ launch_optimjob --cell_id xyz --conda_env CompNeuro --ext_scripts $script_repo 
$ submit_opt_jobs -f cell_data.csv -r $script_repo -c CompNeuro -m 2 # Launching multiple jobs from a csv file
```


