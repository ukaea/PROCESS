## A short summary of what can you find in this directory:

#### manual_start - contains examples of stellarator runs. Some of them should even work. You can run them by run_me.py script. To select specific case you have to change the 'prefix' variable in the script to match the prefic of the input and stella_config file (output will be given the same prefix)

#### autostart - contains scripts to generate a scan over bt  values. To generate a new case, you can use start.py script. It will execute the remaining scripts in the right order, keeping the names consistant. You can use other scripts as well to performe specific actions.

#### templates and config_files - contain the files use by autostart scripts to generate run subdirectories.