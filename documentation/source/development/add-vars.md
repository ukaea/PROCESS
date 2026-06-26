# Guide for adding Variables & Constraints 

Specific instructions must be followed to add an input, iteration variable,
optimisation figure of merit and constraints to the `PROCESS` code.

  **At all times the [`PROCESS` style guide](../development/standards.md) must be used.**


-----------------

## Add an input

To add a `PROCESS` input, please follow below:

1. Choose the most relevant module `XX` and add the variable in the `XX_variables` defined in `XX_variables.py`.
 
2. Add a description of the input variable below the declaration, using the formatting described in the standards section specifying the units.

3. Assign a sensible default initial value, and a type.
  
4. Add the parameter to the `INPUT_VARIABLES` dictionary in `input.py`.  

Here is an example of the code to add:
  

Variable definition and initial value setting example in `tfcoil_variables.py`:
```python
  e_tf_coil_magnetic_stored: float = 0.0
  """Stored magnetic energy in a single TF coil (J)"""
```

Code example in the `input.py` file:

```python
  INPUT_VARIABLES = {
  ...
    "rho_tf_joints": InputVariable("tfcoil", float, range=(0.0, 0.01)),
```

-----------------

## Add an iteration variable

To add a `PROCESS` iteration variable please follow the steps below, in addition to the instructions for adding an input variable:


1. The parameter `IPNVARS` in module `numerics` of `numerics.py` will normally be greater than the actual number of iteration variables, and does not need to be changed.
2. Append a new iteration number key to the end of the `ITERATION_VARIABLES` dictionary  in `iteration_variables.py`. The associated variable is the corresponding key value.
3. Set the variable origin file and then the associated lower and upper bounds
4. Update the `lablxc` description in `numerics.py`.
  
It should be noted that iteration variables must not be reset elsewhere in the
code. That is, they may only be assigned new values when originally
initialised (in the relevant module, or in the input file if required).
Otherwise, the numerical procedure cannot adjust the value as it requires, and
the program will fail.

Here is a code snippet showing how `rmajor` is defined in `iteration_variables.py`

```python
ITERATION_VARIABLES = {
  ...
  3: IterationVariable("rmajor", "physics", 0.1, 50.00),
```

-----------------

## Add a figure of merit

New figures of merit are added to `PROCESS` in the following way:

1. Increment the parameter `IPNFOMS` in module `numerics` in source file `numerics.py` to accommodate the new figure of merit.
  
2. Assign the new integer value and description string of the new figure of merit to the `FiguresOfMerit` enumerator in `numerics.py`.
  
3. Add the new figure of merit equation to `objective_function()` in `objectives.py`, following the method used in the existing examples. The value of figure of merit case should be of order unity, so select a reasonable scaling factor if necessary. 
  
An example can be found below:


```python
objective_function():
  ...
  try:
      figure_of_merit = FiguresOfMerit(abs(minmax))
  ...  
  if figure_of_merit == FiguresOfMerit.MAJOR_RADIUS:
        objective_metric = 0.2 * data.physics.rmajor
```

-----------

## Add a scan variable

After following the instruction to add an input variable, you can make the variable a scan variable by following these steps:

1. Increment the parameter `IPNSCNV` defined in `scan_variables.py` in the data_structure directory, to accommodate the new scanning variable. The incremented value will identify your scan variable.
  
2. Add a short description of the new scanning variable in the `nsweep` comment in `scan_variables.py`, alongside its identification number.
  
3. Update the `ScanVariables` enum in the `scan.py` file by adding a new case statement connecting the variable to the scan integer switch, the variable name and a short description.
  
4. Add a comment in the corresponding variable file in the data_structure directory, eg, `data_structure/[XX]_variables.py`, to add the variable description indicating the scan switch number.
  

`nsweep` comment example:
```python
  nsweep: int = None
  """Switch denoting quantity to scan:<UL>
  <LI> 1  aspect
  <LI> 2  pflux_div_heat_load_max_mw
  <LI> 3  p_plant_electric_net_required_mw
  <LI> 4  hfact
```

`nsweep_dict` example in `scans.py`
```python
nsweep_dict = {
        1: "aspect",
        2: "pflux_div_heat_load_max_mw",
        3: "p_plant_electric_net_mw",
        4: "hfact",
        ...
}
```


`ScanVariables` case example:

```python
  class ScanVariables(Enum):
    
    ...

    aspect = ScanVariable("aspect", "Aspect_ratio", 1)
    pflux_div_heat_load_max_mw = ScanVariable(
        "pflux_div_heat_load_max_mw", "Div_heat_limit_(MW/m2)", 2
    )
    p_plant_electric_net_required_mw = ScanVariable(
        "p_plant_electric_net_required_mw", "Net_electric_power_(MW)", 3
    )
```



---------------

## Add a constraint equation

Constraint equations are added to *PROCESS* in the `process/core/solver/constraints.py` file. They are registered with the `ConstraintManager` whenever the application is run. Each equation has a unique name that is currently an integer, however upgrades to the input file format in the future will allow arbitrary hashable constraint names. 

A constraint is simply added by registering the constraint to the manager using a decorator.

```python
@ConstraintManager.register_constraint(1234, "m", "=")
def my_constraint_function(constraint_registration): ...
```
The arguments to the `register_constraint` function are:

- Name (again, currently an integer)
- Unit (for output reporting purposes)
- Symbol (e.g. =, >=, <=. Again, for output reporting purposes)


`my_constraint_function` should be named appropriately and return a `ConstraintResult` which contains the:

- Normalised residual error
- Constraint value
- Constraint bound
- Constraint residual

The recommended way to do this is using one of the functions `geq`, `leq`, or `eq` depending on whether the constraint is desired to be $v\geq b$, $v\leq b$, or $v=b$, respectively.

```python
@ConstraintManager.register_constraint(1234, "m", "=")
def my_constraint_function(constraint_registration):
  return geq(value, bound, constraint_registration)
```
