### Running the optimization jobs using shared conda environment on NERSC

The location of the main code base on NERSC filesystem is

`/project/projectdirs/m2043/AIBS/ani/Unified/`

#### Shared python environment

Add the following line to your .bashrc.ext on the home directory 

`export PATH="/global/common/software/m2043/AIBS_Opt/software:$PATH"`

Log out and log back in or use

`source ~/.bashrc.ext`

for changes to take effect. The shared conda environment ateam with all necessary packages should be available at this time.

`source activate ateam`

#### Create a directory with the same name as cell id 

Note that for cells part of the online product the name has to match the cell id exactly.

`mkdir <cell_id>`

`cd <cell_id>`

Copy the file `move_opt_files.sh` (available @ `/project/projectdirs/m2043/AIBS/ani/Unified/move_opt_files.sh`) in this directory and run this file

`sh move_opt_files.sh #this will copy the necessary starter files to the current directory`

#### Run the save_cell_metadata.py with suitable options

For cells part of the online product run

`python save_cell_metadata.py --launch_job`

For DG cells provide the .nwb and .swc path using --nwb_path and --swc_path options. For manually setting other fields in the metadata (this should be populated automatically for online cells) run 

`python save_cell_metadata.py --help`

for help.

