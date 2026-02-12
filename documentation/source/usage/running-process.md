# Running PROCESS

There are a number of ways to run PROCESS.  The first two are determined by the value of the switch `ioptimz` in the input file:

`ioptimz = -2` for evaluation. The physics and engineering models are evaluated and the equality (e.g. model consistency) constraints will be solved using `scipy`'s `fsolve`. The input file and default values together define the input values of all variables. The values of the parameters used to solve the equalities (solution parameters) will be different in the solution, however. This is used when evaluating a set of input parameters (e.g. a "point") whilst ensuring that the models are self-consistent.

`ioptimz = 1` for optimisation.  Those variables specified as iteration variables (specified by equations such as `ixc = 1` in the input file) are automatically varied during the iteration process, between the bounds given by the arrays `boundl` (lower bounds) and `boundu` (upper bounds).  The input file contains the *initial* values of the iteration variables, but the final values will not be same.  The iteration process continues until convergence, or until the maximum number of iterations (`maxcal`) is reached.  If the code converges, the constraints will be satisfied and the Figure of Merit will be maximised or minimised. 

If the optimisation fails to converge, a third option is available by using the command line option `-v` (VaryRun), together with a configuration file.  PROCESS is run repeatedly in optimisation mode (`ioptimz` = 1 must be set), but a new input file is written each time, with different, randomly selected *initial* values of the iteration variables.  The factor within which the initial values of the iteration variables are changed is `FACTOR` (in the configuration file).  For example, `FACTOR = 1.1` will vary the initial values randomly by up to 10%.  This is repeated until PROCESS converges, or until the maximum number of PROCESS runs (`NITER`) is reached.  Sometimes this procedure will generate a converged solution when a single optimisation run does not.

A SCAN is available in any of these modes.  One input variable can be scanned (`scan_dim = 1`) or two input variables (`scan_dim = 2`).  A scan variable must not be an iteration variable.  For details, see [scan_module](source/reference/process/scan.md).

--------------

## To run PROCESS

The default PROCESS input file name is IN.DAT. If no input file name is given as an argument in the command line, it assumes an IN.DAT file is present in the current directory:
```bash
# Use an IN.DAT file in the current directory
process
```
In this case the output files will be `OUT.DAT` and `MFILE.DAT`.

The user can provide a named input file, provided the last 6 characters of the input file name are IN.DAT.
```bash
process -i path/to/my_file_name_IN.DAT
```
will produce the following output files in the same directory as the input file:
```
    my_file_name_OUT.DAT
    my_file_name_MFILE.DAT
```

It may be convenient to automatically generate the summary plots after the `PROCESS` run has finished.
Setting the `--full-output` flag:

```bash
process -i path/to/my_file_name_IN.DAT --full-output
```

will produce the following output files in the same directory as the input file:

```
    my_file_name_OUT.DAT
    my_file_name_MFILE.DAT
    SankeyPowerFlow.pdf
    my_file_name.MFILE.DATSUMMARY.pdf
    my_file_name.MFILE_radial_build.pdf
```

---------------

### VaryRun 

The default VaryRun configuration filename is `run_process.conf`. If no configuration filename is given as an argument in the command line, `run_process.conf` is assumed to be present in the current directory:
```bash
# Use a configuration file called run_process.conf in the current directory
process -v
```
The user can provide a named configuration file:
```
process -v -c path/to/my_conf_file.conf
```
When using VaryRun, the input filename is defined inside the configuration file.  The `-i` argument should not be used.

-----------------

### Command line arguments

The full set of command line arguments is available with:

```bash
process --help
```

------------

## Configuration file for VaryRun

The configuration file has the following format:

```
* CONFIG FILE FOR RUNNING PROCESS WITH MODIFIED IN.DAT
* (Does not allow for inline comments!)
* 

* Path to working directory in which PROCESS is run.
WDIR = .

* original IN.DAT name (should not be called IN.DAT!)
ORIGINAL_IN_DAT = large_tokamak_IN.DAT

* ONE line comment to be put into README.txt
COMMENT = 

* Maximum number of runs
NITER = 1000

* integer seed for random number generator; use None for random seed
SEED = 5

* factor within which the initial values of the iteration variables are changed
FACTOR = 1.5

* Number of allowed unconverged scan points in a single PROCESS scan.
* Leave this at 0 to stop the vary run only when a converged solution is found
NO_ALLOWED_UNFEASIBLE = 0

* include a summary file with the iteration variables at each stage.
INCLUDE_ITERVAR_DIFF = True

* The following parameters are only required if you wish to add or remove iteration variables
* or constraints from the input file supplied.
* Otherwise leave these lines unchanged.

* List of additional iteration variables to be added to IN.DAT - comma separated
ADD_IXC = 

* List of iteration variables to be deleted from IN.DAT - comma separated
DEL_IXC = 

* List of additional constraint equations to be added to IN.DAT  - comma separated 
ADD_ICC = 

* List of constrained equations to be deleted from IN.DAT - comma separated
DEL_ICC = 
```
