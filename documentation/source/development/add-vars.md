# Guide for adding Variables & Constraints 

Specific instructions must be followed to add an input, iteration variable,
optimisation figure of merit and constraints to the `PROCESS` code.

  **At all times the [`PROCESS` style guide](../development/standards.md) must be used.**

!!! note

    As the code is quickly converging towards a wholly Python codebase the respective files may change in type from `.f90` to `.py`.

-----------------

## Add an input

To add a `PROCESS` input, please follow below:

1. Choose the most relevant module `XX` and add the variable in the `XX_variables` defined in `XX_variables.f90`.
 
2. Add a description of the input variable below the declaration, using the FORD      formatting described in the standards section specifying the units.
  
3. Specify a sensible default value in the `init_XX_variables()` function within the corresponding model `.py` main file
  
4. Add the parameter to the `INPUT_VARIABLES` dictionary in `input.py`.  

Here is an example of the code to add:
  

Variable definition example in `tfcoil_variables.f90`:
```fortran
  real(dp) :: rho_tf_joints
  !! TF joints surfacic resistivity [ohm.m]
  !! Feldmetal joints assumed.
```

Variable initialization example in `tf_coil.py`:
```python
    def init_tfcoil_variables():
    ...
      tfv.rho_tf_joints = 2.5e-10
```

Code example in the `input.py` file:

```python
  INPUT_VARIABLES = {
  ...
    "rho_tf_joints": InputVariable(data_structure.tfcoil_variables, float, range=(0.0, 0.01)),
```

-----------------

## Add an iteration variable

To add a `PROCESS` iteration variable please follow the steps below, in addition to the instructions for adding an input variable:


1. The parameter `ipnvars` in module `numerics` of `numerics.f90` will normally be greater than the actual number of iteration variables, and does not need to be changed.
2. Append a new iteration number key to the end of the `ITERATION_VARIABLES` dictionary  in `iteration_variables.py`. The associated variable is the corresponding key value.
3. Set the variable origin file and then the associated lower and upper bounds
4. Update the `lablxc` derscription in `numerics.f90`.
  
It should be noted that iteration variables must not be reset elsewhere in the
code. That is, they may only be assigned new values when originally
initialised (in the relevant module, or in the input file if required).
Otherwise, the numerical procedure cannot adjust the value as it requires, and
the program will fail.

Here is a code snippet showing how `rmajor` is defined in `iteration_variables.py`

```python
ITERATION_VARIABLES = {
  ...
  3: IterationVariable("rmajor", fortran.physics_variables, 0.1, 50.00),
```

-----------------

## Add a figure of merit

New figures of merit are added to `PROCESS` in the following way:

1. Increment the parameter `ipnfoms` in module `numerics` in source file `numerics.f90` to accommodate the new figure of merit.
  
2. Assign a description of the new figure of merit to the relevant element of array `lablmm` in module `numerics` in the source file `numerics.f90`.
  
3. Add the new figure of merit equation to `objective_function()` in `objectives.py`, following the method used in the existing examples. The value of figure of merit case should be of order unity, so select a reasonable scaling factor if necessary. 
  
4. Add the new figure of merit description to the `OBJECTIVES_NAMES` dictionary in `objectives.py`
  
An example can be found below:


```python
objective_function():
  ...
  match figure_of_merit:
  ...  
  case 1:
        objective_metric = 0.2 * physics_variables.rmajor
```

-----------

## Add a scan variable

After following the instruction to add an input variable, you can make the variable a scan variable by following these steps:

1. Increment the parameter `IPNSCNV` defined in `scan_variables.py` in the data_structure directory, to accommodate the new scanning variable. The incremented value will identify your scan variable.
  
2. Add a short description of the new scanning variable in the `nsweep` comment in `scan_variables.py`, alongside its identification number.
  
3. Update the `ScanVariables` enum in the `scan.py` file by adding a new case statement connecting the variable to the scan integer switch, the variable name and a short description.
  
4. Add a comment in the corresponding variable file in the data_structure directory, eg, `data_structure/[XX]_variables.py`, to add the variable description indicating the scan switch number.
  

`nsweep` comment example:
```fortran

  integer :: nsweep = 1
  !! nsweep /1/ : switch denoting quantity to scan:<UL>
  !!         <LI> 1  aspect
  !!         <LI> 2  pflux_div_heat_load_max_mw
  ...
  !!         <LI> 54 GL_nbti upper critical field at 0 Kelvin
  !!         <LI> 55 `dr_shld_inboard` : Inboard neutron shield thickness </UL>
```

`SCAN_VARIABLES` case example:

```python
  class ScanVariables(Enum):
    aspect: ScanVariable("aspect", "Aspect_ratio", 1),
    pflux_div_heat_load_max_mw: ScanVariable("pflux_div_heat_load_max_mw", "Div_heat_limit_(MW/m2)", 2),
    ...
    Bc2_0K: ScanVariable("Bc2(0K)", "GL_NbTi Bc2(0K)", 54),
    dr_shld_inboard : ScanVariable("dr_shld_inboard", "Inboard neutronic shield", 55),
```

---------------

## Add a constraint equation

Constraint equations are added to *PROCESS* in the `process/constraints.py` file. They are registered with the `ConstraintManager` whenever the application is run. Each equation has a unique name that is currently an integer, however upgrades to the input file format in the future will allow arbitrary hashable constraint names. 

A constraint is simply added by registering the constraint to the manager using a decorator.

```python
@ConstraintManager.register_constraint(1234, "m", "=")
def my_constraint_function(): ...
```
The arguments to the `register_constraint` function are:

- Name (again, currently an integer)
- Unit (for output reporting purposes)
- Symbol (e.g. =, >=, <=. Again, for output reporting purposes)


`my_constraint_function` should be named appropriately and return a `ConstraintResult` which contains the:

- Normalised residual error
- Constraint value
- Constraint error

```python
@ConstraintManager.register_constraint(1234, "m", "=")
def my_constraint_function():
  normalised_residual = ...
  value = ...
  error = ...
  return ConstraintResult(normalised_residual, value, error)
```
