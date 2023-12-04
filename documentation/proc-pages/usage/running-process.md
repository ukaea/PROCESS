# Running PROCESS

Process can be run in two main mode:

1. SingleRun (where Process runs once)
2. VaryRun (where PROCESS runs multiple times, varying iteration parameters until a solution is found).

For a SingleRun:

```bash
# Use an IN.DAT file in the current directory
process

# Use an IN.DAT file outside of the current directory
process -i path/to/IN.DAT 
```

For a VaryRun:

```bash
# Use a configuration file called run_process.conf in the current directory
process -v

# Use a conf file outside of the current directory
process -v -c path/to/my_conf_file.conf
```

The available arguments are:

```bash
process [--input, -i input_file_path] [--varyiterparams, -v] [--varyiterparamsconfig, -c config_file_path] [--help, -h]
```

Help is available with:

```bash
process --help
```

## Config file for VaryRun

### Configuration File

The configuration file `.conf` has the following style:

```
* CONFIG FILE FOR RUNNING PROCESS WITH MODIFIED IN.DAT
* (Does not allow for inline comments!)
* 

* Path to working directory in which PROCESS is run.
WDIR = .

* original IN.DAT name (should not be called IN.DAT!)
ORIGINAL_IN_DAT = ref_IN.DAT

* Not used any more, Python package used instead
* PATH to PROCESS binary
* PROCESS= ../../../process.exe

* ONE line comment to be put into README.txt
COMMENT = 

* Max no. iterations
NITER = 1000

* integer seed for random number generator; use None for random seed
SEED = None

* factor within which the iteration variables are changed
FACTOR = 1.9

* Number of allowed unfeasible points that do not trigger rerunning.
NO_ALLOWED_UNFEASIBLE = 0

* include a summary file with the iteration variables at each stage.
INCLUDE_ITERVAR_DIFF = True

* add iteration variables - comma separated
ADD_IXC = 

* remove iteration variables - comma separated
DEL_IXC = 

* add constraint equations  - comma separated 
ADD_ICC = 

* remove constraint equations - comma separated
DEL_ICC = 
```
