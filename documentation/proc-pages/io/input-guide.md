# Input File

The input file `IN.DAT` is used to change the values of the physics, engineering 
and other code parameters from their default values, and to set up the constraint equations, 
iteration variables etc. 

If the code encounters a problem reading the input file, it will stop immediately 
with an error message. The last line of the output file `OUT.DAT` may give an 
indication of where in the input file the problem lies.

## File Naming Convention

The default PROCESS input file name is `IN.DAT`. The user can provide a named 
input file, that will produce named output files, provided the last 6 characters 
of the input file name are `IN.DAT`.

```bash
process -i [path-to-file]/my_file_name_IN.DAT
```

Will produce output files named:

- `my_file_name_OUT.DAT`
- `my_file_name_MFILE.DAT`

If no input file name is given as the first argument to the code, it assumes an 
`IN.DAT` file is present in the current directory.

## Constraints

PROCESS permits a large number of constraint equations, all of which are formulated 
in the source file `constraint equations.f90`. 

In the input file constraint equations are specified as in the following example:

```
icc = 2 
```

where `icc` is the constraints array in PROCESS and the user is requesting constraint 
equation 2. A comment on the same line is recommended:

```
icc = 2 * Global power balance (consistency equation)
```

Some constraints have `f-value` variables. These must be set as iteration variables, 
which are discussed below.


!!! Info "Constraints"  
    A full list of constraints is given on the variable description page in the row labelled 
    `lablcc` [here](../vardes/#numerics).  
    See [solver](../solver/solver-guide) page for more info

## Iteration Variables

Variables that are adjusted by PROCESS in order to satisfy the constraints and 
optimise the figure of merit are referred to as iteration variables. Successive calls 
are made to the physics and engineering routines, with slightly different values for 
the iteration variables on each call, and the equation solver determines the effect on the 
output due to these small changes to the input.

Iteration variables must never be initialised to zero. The code will not be able to adjust 
the variable’s value if this is done, and it will stop with an error message.

An iteration variable can be specified in the input file as in the following example:

```
ixc = 3
```

where `ixc` is the iteration variable array in PROCESS. An in-line comment is recommended:

```
ixc = 3 * Plasma major radius [m]
```

For example, the major radius is available as an iteration variable, and appears in the variable
description file as `rmajor /8.14/ : plasma major radius (m) (iteration variable 3)`. If it
is selected as an iteration variable, it will be adjusted by the code. The value input by the user (or
the default, if no value is specified), will be used as the starting value.

!!! Note "Iteration Variables"  
    A full list of iteration variables is given on the variable description page in the row labelled 
    `lablxc` [here](../vardes/#numerics).  
    (See [solver](../solver/solver-guide) page for more info)

## Bounds

The upper and lower bounds of an iteration variable can be set by the user in 
the input file as in the following example, 
where `boundl` is the lower bound and `boundu` is the upper bound:

```
boundl(3) = 8 
boundu(3) = 12
```

where `3` is the iteration variable number (in this case the major radius). It is good practice to
place the iteration variable and its bounds in a block:

```
ixc = 3 * Plasma major radius [m]
boundl(3) = 8
boundu(3) = 12
```
If bounds are not specified default values are used.

## Numerics

The user can select which solver to use, but only one solver is available at present (VMCON).

```
ioptimz  = 1 * for optimisation VMCON only
```

The user can select the figure of merit to be used:

```
minmax   = 1 * Switch for figure-of-merit (see lablmm for descriptions)
```

In this case the user is choosing option `1`, which is major radius. For `minmax`

* a **positive** value means **minimise** the figure of merit
* a **negative** value means **maximise** the figure of merit

The user can also input the allowed error tolerance on the solver solution:

```
epsvmc   = 1.0e-8 * Error tolerance for vmcon
```

!!! Info "Figure of Merit"  
    A full list of figures of merit is given on the variable description page in the row labelled 
    `lablmm` [here](../vardes/#numerics).  

## Input Variables

One can enter an input into the `IN.DAT` by:

```
rmajor = 8.90 * Plasma major radius [m]
```

The `*` is for adding comments to the input file. To comment out an entire line 
one can add a `*` to the beginning of the line, as below:

```
*rmajor = 8.90 * Plasma major radius [m]
```

!!! Info "Variable Descriptions"  
    A full list of inputs variables is given in the PROCESS `html` documentation 
    file `vardes.html` and on the variable description page [here](../vardes).

## Scan

PROCESS can scan either one or two variables within a single run.  
This provides a method of determining the sensitivity of the
results to different input assumptions. 

1-D scan.  The user specifies which variable is to be scanned as in the following example:

```
nsweep = 1 
isweep = 4
sweep = 2.8, 2.9, 3.0, 3.1
```

where `nsweep` is the scan variable chosen (see [variable descriptions](../vardes)),
`isweep` is the number of scan points and `sweep` is the array of scan values. In this example, 
PROCESS runs for each of the four values given of the plasma aspect ratio variable `aspect`. 

2-D scan.  Two variables are scanned over a rectangular grid of values, 
using the switch `scan_dim = 2` as below

```
scan_dim = 2

nweep = 1
isweep = 4
sweep = 2.8, 2.9, 3.0, 3.1

nweep_2 = 4
isweep_2 = 3
sweep_2 = 1.0, 1.1, 1.2
```

where the scan parameters have duplicate names with `_2` for the second scan 
dimension.  In this example, PROCESS runs for each of the four values given 
of the plasma aspect ratio variable `aspect` (specified by `nsweep`=1) and the three values given for
the energy confinement time enhancement factor `hfact` (specified by `nsweep_2`=4).

The results from the previous scan point are used as the input to the next
scan point. The output files contain all of the scan points for a given run.

Scan variables should not be confused with iteration variables.

!!! Info "Scanning"    
    For obvious reasons, the active scanning variable or variables must not also be active
    iteration variables.   
    A full list of scan variables is [here](../vardes/#scan_module).  

## Examples

Example `IN.DAT` files are available in the repository in the folder `/tests/regression/scenarios/`.
