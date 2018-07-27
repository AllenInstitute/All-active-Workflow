### Directory structure for the scripts to work:

* Put the .nwb and .swc file in a directory named 'cell_types'
* Then run starter_optim.py

> Make sure you have the following directory structure:


    .
    ├── ...
    ├── cell_types                             # Test files (alternatively `spec` or `tests`)
    │   ├── .nwb                               # Load and stress tests
    │   ├── .swc                               # End-to-end, integration tests (alternatively `e2e`)
    │   └── fit_parameters.json                # Unit tests
    │── config ── cell_id                                 # Test files (alternatively `spec` or `tests`)
    │              ├── .nwb                               # Load and stress tests
    │              ├── .swc                               # End-to-end, integration tests (alternatively `e2e`)
    │              └── fit_parameters.json
