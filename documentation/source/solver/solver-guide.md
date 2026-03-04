
# Constraint Equations

Any computer program naturally contains many equations. The built-in equation 
solvers within PROCESS act on a special class, known as constraint equations, 
all of which are formulated in the source file `process/core/optimisation/constraints.py`. These 
can be split into two types:
 
**Equality constraints (consistency equations)** that enforce consistency between the physics and 
engineering parameters in the various models.

**Inequality constraints (limit equations)** that enforce a value be greater/less than or equal to some bound (limit).

A PROCESS input file will, for example, define which constraint equations are being used as follows:

```
...
neqns = 3

* Equalities
icc = 1 * Beta
icc = 2 * Global power balance
icc = 11 * Radial build

* Inequalities
icc = 9 * Fusion power upper limit
icc = 5 * Density upper limit
icc = 24 * Beta upper limit
icc = 15 * LH power threshold limit
...
```

Here each `icc=n` statement tells PROCESS to activate a constraint with the name `n`. A list of the constraints and 
their corresponding names can be found [here](../../source/reference/process/data_structure/numerics/#process.data_structure.numerics.lablcc).

The `neqns = 3` statement is telling PROCESS to treat the first `3` equations as equality constraints, and the rest as inequality constraints. Therefore, it is imperative that all equality constraints are stated before any inequality constraints.


In both types of equations, an optimiser/solver uses the normalised residuals $c_i$ of the constraints (and sometimes its gradient, depending on the solver/optimiser) to guide the solution towards one that satisfies all of the constraints.

## Consistency Equations

Consistency equations are equalities that ensure that the machine produced by PROCESS is 
self-consistent. This means, therefore, that many of these constraint equations should 
always be used, namely equations 1, 2, 10 and 11. Equation 7 should also be activated 
if neutral beam injection is used. The other consistency equations can be activated if 
required. A typical consistency equation ensures that two functions $g$ and $h$ are equal:

$$
g(x, y, z, \cdots) = h(x, y, z, \cdots)
$$

$$
c_i = 1 - \frac{g}{h}
$$

The optimiser/solver will attempt to find a solution that produces $c_i = 0$ for all equality constraints.

## Limit Equations

The limit equations are inequalities that ensure that various physics or engineering 
limits are not exceeded. 

As with the consistency equations, the general form of the limit equations is

$$
c_i = 1 - \frac{h_{max}}{h}
$$

where $h_{max}$ is the maximum allowed value of the quantity $h$, or

$$
c_i = 1 - f.\frac{h}{h_{min}}
$$

where $h_{min}$ is the minimum allowed value of the quantity $h$.

For example, to set the net electric power to a certain value, the following 
should be carried out:

1. Activate `constraint 16` (net electric power lower limit) by including it in the `icc` array
2. Set `p_plant_electric_net_required_mw` to the required net electric power.

The optimiser/solver will attempt to find a solution that produces $c_i \geq 0$ for all inequality constraints.

## Iteration Variables

...

## Figure of Merit

In optimisation mode, PROCESS finds the self-consistent set of iteration 
variable values that maximises or minimises a certain function of them, 
known as the figure of merit. 

Several possible figures of merit are available, all of which are in the 
source file `evaluators.f90`. 

Switch `minmax` is used to control which figure of merit is to be used. If the 
figure of merit is to be minimised, `minmax` should be **positive**, and if a 
maximised figure of merit is desired, `minmax` should be **negative**.

## Convergence

...

## Optimisation mode

Switch `ioptimz` should be set to 1 for optimisation mode.

If `ioptimz = 0`, a non-optimisation pass is performed first. Occasionally this provides a feasible set of initial conditions that aids convergence of the optimiser, but it is recommended to use `ioptimz = 1`.

Enable all the relevant consistency equations, and it is advisable to enable the corresponding iterations variables. A number of limit equations (inequality constraints) can also be activated. In optimisation mode, the number of iteration variables is unlimited.

It may still be difficult, if not impossible, to reconcile the fusion power and the net electric power with the required values. This may well be due to the power conversion efficiency values being used.

If scan of a given variable are to be made over a large range of values, it is often a good idea to start the scan in the middle of the desired range, and to split the in two - one going downwards from the initial value,  and the other upwards. This ensures that the whole range of the scan produces well-converged machines (assuming a "good" initial point), without sharp changes in gradient in the parameter values.

It should be remembered that the value of the scan variable is set in the array `sweep`, and this overrules any value set for the variable elsewhere in the input file.

The output from an optimisation run contains an indication as to which iteration variables lie at their limit values.

## Evaluation mode

Evaluation mode is used to evaluate models for a given set of input parameters (a "point") whilst ensuring that the models are self-consistent. It can also be used to perform benchmark comparison, whereby the machine size, output power etc. are known and one only wishes to find the calculated stresses, beta values and fusion powers, for example.

Running `PROCESS` in evaluation mode requires few changes to be made to the input file from the optimisation case. The main differences between optimisation mode and evaluation mode are:

1. Evaluation mode does not apply lower or upper bounds to the iteration variables. 

1. Inequality constraints are not enforced.

2. The number of iteration variables must be equal to the number of equality constraints being solved.

3. An objective function is not available, as no optimisation is taking place.

As before, the user must decide which constraint equations and iteration variables to activate. For example, an extract from an input file might look like:
```
* Evaluation problem: evaluate models consistently by solving equality constraints only
ioptimz  = -2 * evaluation mode

*---------------Constraint Equations---------------*
* Define number of equality constraints
neqns = 2

* Equalities
icc = 1 * Beta
icc = 2 * Global power balance

* Inequalities
* Not enforced, but values reported
icc = 5 * Density upper limit
icc = 8 * Neutron wall load upper limit
icc = 9 * Fusion power upper limit
...

*---------------Iteration Variables----------------*
* Used to solve equality constraints
ixc = 4 * temp_plasma_electron_vol_avg_kev
ixc = 6 * nd_plasma_electrons_vol_avg
...
```
The beta and plasma power balance equality constraints are recommended to ensure model consistency, and the `temp_plasma_electron_vol_avg_kev` and `nd_plasma_electrons_vol_avg` iteration variables used to solve them.
