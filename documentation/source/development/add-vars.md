# Guide for adding variables and constraints

A guide for how to add new inputs, constraints, figures of merit, and scan variables. 

**At all times ensure new variables and functions adhere to the [`PROCESS` style guide](../development/standards.md).**


All of these features rely on 'variables' which belong to a 'data structure'. All of the data structures can be found in `process/data_structure`.

In general, a variable within a data structure could act as an:

- Input variable: is specified by the user in an `IN.DAT` and its value is not changed once PROCESS is running.
- Iteration variable: is modified by the solver to try and optimise for some figure of merit.
- Scan variable: is sequentially modified by the `Scan` class to some `IN.DAT`-defined values.
- Intermediate variable: is calculated within a model and then used within other models.
- Output variable: is calculated within a model and then written out to the `MFILE.DAT` and `OUT.DAT`.

It is advised that a variable is either used to define a particular PROCESS run (input variable, iteration variable, or scan variable) or mutated within a PROCESS run (output variable or intermediate variable). Mixing the two classes of variable (e.g. having a variable that can be input but is also mutated within a model) will lead to confusing, dangerous, and incorrect results.


-----------------

## Add a new variable
You may need to add a variable to PROCESS when changing or creating models. In most cases, you will want to add your variable to an existing data structure. Creating an entirely new data structure is beyond the scope of this guide, so please seek support from the PROCESS maintainers.

For example, if you are adding a new variable that relates to the blanket model, you would add the variable to `process/data_structure/blanket_variables.py` as part of the `BlanketData` dataclass. 

```python
@dataclass(slots=True)
class BlanketData:
  <... existing variables ...>

  my_new_blanket_variable: float = 0.0
  """my variable description [m]"""
```

Here, `[m]` is the units and should be replaced with the appropriate units for the variable being added.

This variable could then be used within a model

```python
self.data.blanket.my_new_blanket_variable = 1.0
...
another_variable = self.data.blanket.my_new_blanket_variable / 2.0
```

-----------------

## Add a new input
Adding an input in PROCESS means that some variable in a data structure can be set from the `IN.DAT`. Inputs are defined in the `process/core/input.py` file in the `INPUT_VARIABLES` dictionary. Adding a new entry to this dictionary will create a new input.

Continuing with the example from the previous section:
```python
INPUT_VARIABLES = {
  ...
  "my_new_blanket_variable": InputVariable("blanket", float),
}
```

`InputVariable` has several additional fields to support validation and the parsing of arrays; please consult the dataclass for these additional arguments.

You would replace `"blanket"` with the name of the data structure your specific variable belongs to (found by looking at `DataStructure` in `process/core/model.py`).

Now, in the `IN.DAT`, you could set an initial value for `my_new_blanket_variable` by writing:

```
my_new_blanket_variable = 1.0
```

Equally, you could choose not to specify the variable in the `IN.DAT` and the value would be intialised to its default value (`0.0` in this example, as we specified [earlier](#add-a-new-variable)).

-----------------

## Add an iteration variable

Adding an iteration variable allows the PROCESS solver to change the variable as part of the optimisation/solving loop. Iteration variables are defined in `process/core/solver/iteration_variables.py` in the `ITERATION_VARIABLES` dictionary. You would add a new entry to this dictionary to create a new iteration variable:


```python
ITERATION_VARIABLES = {
    123: IterationVariable("my_new_blanket_variable", "blanket", 0.1, 1.0),
}
```

In this example:

- `123` is the identifier of the iteration variable, and must be unique.
- `"blanket"` is the data structure the variable will be set on.
- `0.1` is the default lower bound of the variable.
- `1.0` is the default upper bound of the variable.

You will often want to [add a variable to the input file](#add-a-new-input) if it is an iteration variable. That way, you can specify the initial value of the iteration variable in the `IN.DAT`.

!!! note
    Iteration variables are the exception to the best-practice of not allowing inputs to be mutated. The solver will obviously need to modify iteration variables away from their default value, however it remains true that _models_ should not mutate inputs. Furthermore, models should not change the values of iteration variables set by the solver.

The iteration variable can be enabled in the `IN.DAT` by:
```
ixc = 123

my_new_blanket_variable = 0.5
```

Again, note that omitting the line specifying the initial value for an iteration variable means the initial value would be set to that variable's default value that we specified [earlier](#add-a-new-variable).

-----------------

## Add a figure of merit

A figure of merit is the scalar that the optimiser (e.g. VMCON) will try and minimise or maximise. The figures of merit are specified in `process/core/solver/objectives.py` in the `objective_function()` function.

To add a new figure of merit, first create a new entry in the `FiguresOfMerit` enum in `process/data_structure/numerics.py`:

```python
class FiguresOfMerit(IntEnum):
    ...
    BLANKET_FIGURE_OF_MERIT = (20, "my FOM description")
```
Here `20` will be the identifier of the figure of merit, and **must** be unique. 

Finally, add the equation to `process/core/solver/objectives.py`:
```python
elif figure_of_merit == FiguresOfMerit.BLANKET_FIGURE_OF_MERIT:
  objective_metric = data.blanket.my_new_blanket_variable
```

Note that you will want to scale the `objective_metric` such that it is of the order unity if the variable is not already. This is recommended (but not necessary) because most optimisers are tuned to work well on data of this scale.

The figure of merit can be selected in the `IN.DAT`:
```
minmax = 20
```
Remember, setting `minmax = -20` would minimise instead of maximise our new variable.

-----------------

## Add a scan variable

After following the instructions to add an input variable, you can then make a scan variable.

First, add the variable to the `ScanVariables` enum in `process/core/scan.py`. 
```python
class ScanVariables(Enum):
    ...
    blanket_scan_variable = ScanVariable(
        "my_new_blanket_variable", "A blanket variable", 82
    )
```
Here, `82` is the identifier of the scan variable and must be unique. 

Next, increment the parameter `IPNSCNV` in `process/data_structure/scan_variables.py` and be sure to add a description of the scan variable in the docstring of the `nsweep` variable.

Finally, in `process/core/scan.py`, add the scan variable to the `Scan.scan_select()` method.
```python
match nwp:
  ...
  case 82:
    self.data.tfcoil.my_new_blanket_variable = swp[iscn - 1]
```

Please see the [scan documentation](../io/input-guide.md#input-guide) for how to set up a scan `IN.DAT`

-----------------

## Add a constraint equation

Constraint equations are added to PROCESS in the `process/core/solver/constraints.py` file. They are registered with the `ConstraintManager` whenever the application is run. Each equation has a unique name that is currently an integer.

A constraint is simply added by registering the constraint to the manager using a decorator.

```python
@ConstraintManager.register_constraint(1234, "m", "=")
def my_constraint_function(constraint_registration): ...
```

The arguments to the `register_constraint()` function are:

- Name (again, currently an integer): in the example above it is `1234`
- Unit (for output reporting purposes): in the example above it is `m`
- Symbol (e.g. `=`, `>=`, `<=`. Again, for output reporting purposes): in the example above it is `=`


`my_constraint_function()` should be named appropriately and return a `ConstraintResult` which contains the:

- Constraint residual
- Normalised residual
- Constraint value
- Constraint bound

The recommended way to do this is using one of the functions `geq`, `leq`, or `eq` depending on whether the constraint is desired to be $v\geq b$, $v\leq b$, or $v=b$, respectively, where $v$ is the value (something that PROCESS calculates) and $b$ is the bound (the limit set by the user in the `IN.DAT`).

```python
@ConstraintManager.register_constraint(1234, "m", "=")
def my_constraint_function(constraint_registration):
    return eq(value, bound, constraint_registration)
```
