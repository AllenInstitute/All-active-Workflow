[![Build Status](https://travis-ci.com/anirban6908/All-active-Workflow.svg?token=93Twb9jDYFzVNoM9gSjr&branch=master)](https://travis-ci.com/anirban6908/All-active-Workflow)
[![Generic badge](https://img.shields.io/badge/License-Allen_Institute-yellow.svg)](https://alleninstitute.org/legal/terms-use/)


# All-active model generation and evaluation workflow
Code base for all-active model generation, evaluation and analysis (main engine: bluepyopt). Method developed and present in Nandi et al. (2020): https://www.biorxiv.org/content/10.1101/2020.04.09.030239v1

Main developers: Ani Nandi, Werner van Geit, Tom Chartrand and Costas A. Anastassiou

* Genetic algorithm in action: **selection** + **evaluation** + **evolution**

![alt text](examples/visualization/animations/GA_evolution_animation/movie.gif "all-active model optimization") 

## Dependencies

* Python 3.5+
* Neuron 7.5 (compiled with Python support)
* BluePyOpt
* eFEL
* AllenSDK
* NeuroM
* Pandas
* Seaborn


## Installation
```bash
pip install git+https://github.com/AllenInstitute/All-active-Workflow # install directly from the repository
```
or
```bash
git clone https://github.com/AllenInstitute/All-active-Workflow # clone repository from github
cd All-active-Workflow
pip install -e . # install in editable mode
```


## Quickstart: Launching optimization jobs
In a Unix environment :
```bash
$ source activate conda-env # conda environment with all dependencies
$ launch_optimjob --help # Look at the options
$ launch_optimjob --input_json job_config.json
```

* Configurable staged optimization 
    <details> <summary>job_config.json</summary>

    ```json
    {
        "cty_config": {
            "cell_id": "483101699"
        },
        "job_config": {
            "highlevel_jobconfig": {
                "conda_env": "ateam_opt",
                "axon_type": "stub_axon",
                "data_source": "web",
                "ephys_dir": "ephys_data",
                "non_standard_nwb": false,
                "feature_stimtypes": [
                    "Long Square"
                ],
                "feature_names_path": "feature_set_all.json",
                "compiled_modfiles_dir": "x86_64",
                "job_dir": "483101699_benchmark_timeout"
            },
            "stage_jobconfig": [
                {
                    "stage_name": "Stage0",
                    "stage_stimtypes": [
                        "Long Square"
                    ],
                    "stage_features": "feature_set_stage0.json",
                    "stage_parameters": "param_bounds_stage0.json",
                    "filter_rule": "filter_feat_proto_passive",
                    "offspring_size": 512,
                    "max_ngen": 50,
                    "optim_config":{
                        "nengines": 256,
                        "nnodes": 16,
                        "qos": "celltypes",
                        "nprocs": 16,
                        "error_stream": "job.err",
                        "output_stream": "job.out",
                        "jobmem": "100g",
                        "jobtime": "5:00:00",
                        "ipyparallel": true,
                        "ipyparallel_db": "sqlitedb",
                        "main_script": "Optim_Main.py"
                    },
                    "analysis_config":{
                        "main_script": "analyze_stagejob.py"
                    },
                    "seed": [
                        1
                    ]
                },
                {
                    "stage_name": "Stage1",
                    "stage_stimtypes": [
                        "Long Square"
                    ],
                    "stage_features": "feature_set_stage1.json",
                    "stage_parameters": "param_bounds_stage1.json",
                    "filter_rule": "filter_feat_proto_passive",
                    "offspring_size": 512,
                    "max_ngen": 50,
                    "optim_config":{
                        "nengines": 256,
                        "nnodes": 16,
                        "qos": "celltypes",
                        "nprocs": 16,
                        "error_stream": "job.err",
                        "output_stream": "job.out",
                        "jobmem": "100g",
                        "jobtime": "5:00:00",
                        "ipyparallel": true,
                        "ipyparallel_db": "sqlitedb",
                        "main_script": "Optim_Main.py"
                    },
                    "analysis_config":{
                        "main_script": "analyze_stagejob.py"
                    },
                    "seed": [
                        1
                    ]
                },
                {
                    "stage_name": "Stage2",
                    "stage_stimtypes": [
                        "Long Square"
                    ],
                    "stage_features": "feature_set_stage2.json",
                    "stage_parameters": "param_bounds_stage2_mouse_spiny.json",
                    "filter_rule": "filter_feat_proto_active",
                    "AP_initiation_zone": "axon",
                    "offspring_size": 512,
                    "cp_backup_dir": "checkpoints_backup",
                    "max_ngen": 200,
                    "optim_config":{
                        "nengines": 256,
                        "nnodes": 16,
                        "qos": "celltypes",
                        "nprocs": 16,
                        "error_stream": "job.err",
                        "output_stream": "job.out",
                        "jobmem": "150g",
                        "jobtime": "12:00:00",
                        "ipyparallel": true,
                        "ipyparallel_db": "sqlitedb",
                        "main_script": "Optim_Main.py"
                    },
                    "analysis_config":{
                        "main_script": "analyze_stagejob.py",
                        "ipyparallel": true,
                        "ipyparallel_db": "nodb",
                        "error_stream": "analysis.err",
                        "output_stream": "analysis.out",
                        "nengines": 40,
                        "nnodes": 4,
                        "nprocs": 10,
                        "jobtime": "10:00:00",
                        "jobmem": "100g",
                        "qos": "celltypes"
                    },
                    "seed": [
                        1,
                        2,
                        3,
                        4
                    ],
                    "run_hof_analysis": true,
                    "run_peri_comparison": false,
                    "depol_block_check": true,
                    "add_fi_kink": true,
                    "calc_model_perf": true,
                    "model_postprocess": true,
                    "calc_time_statistics": true,
                    "timeout": 300,
                    "hoc_export": true
                }
            ]
        }
    }

    ```
    </details>

## HPC usage support
Instructions are included to run these compute intensive jobs on 4 different architectures (Allen Institute HPC cluster, NERSC, BBP5 and AWS).


## Level of support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.

