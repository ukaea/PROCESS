
# Constraint Equations

Any computer program naturally contains many equations. The built-in equation 
solvers within PROCESS act on a special class, known as constraint equations, 
all of which are formulated in the source file `constraint equations.f90`. These 
can be split into two types:
 
**Equality constraints (consistency equations)** -- that enforce consistency between the physics and 
engineering parameters in the various models

**Inequality constraints (limit equations)** -- that enforce various parameters to lie within their allowed 
limits. The `neqns` constraint equations that the user chooses for a given run are 
activated by including the equation numbers in the first `neqns` elements of 
array `icc`.

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

The equation solvers VMCON and HYBRD need the constraint equations $c_i$ to be given in 
the form shown, since they adjust the iteration variables so as to obtain $c_i = 0$, 
thereby ensuring that $g = h$.

## Limit Equations

The limit equations are inequalities that ensure that various physics or engineering 
limits are not exceeded. 

The f-values are used as follows. In general, limit equations have the form

$$
\mathrm{calculated\ quantity} = f \times \mathrm{maximum\ allowable\ value}
$$

where $f$ is the `f-value`. In optimisation mode, all iteration variables have prescribed 
lower and upper bounds. If $f$ is chosen to be an iteration variable and is given a 
lower bound of zero and an upper bound of one, then the limit equation does indeed 
constrain the calculated quantity to lie between zero and its maximum allowable value, 
as required. 

As with the consistency equations, the general form of the limit equations is

$$
c_i = 1 - f.\frac{h_{max}}{h}
$$

where $h_{max}$ is the maximum allowed value of the quantity $h$. Sometimes, the limit 
equation and `f-value` are used to ensure that quantity $h$ is larger than its minimum
value $h_{min}$. In this case, $0 ≤ f ≤ 1$ (as before), but the equation takes the form

$$
c_i = 1 - f.\frac{h}{h_{min}}
$$

By fixing the `f-value` (i.e. not including it in the `ixc` array), the limit equations 
can be used as equality constraints. 

For example, to set the net electric power to a certain value, the following 
should be carried out:

1. Activate `constraint 16` (net electric power lower limit) by including it in the `icc` array
2. Set `p_plant_electric_net_required_mw` to the required net electric power.

Limit equations are not restricted to optimisation mode. In non-optimisation mode, the iteration
variables are not bounded, but the `f-values` can still be used to provide information about 
how calculated values compare with limiting values, without having to change the characteristics 
of the device being benchmarked to find a solution.

It is for this reason that all the constraint equations used in PROCESS are formulated as equalities,
despite the fact that equation solver VMCON can solve for inequalities as well. The use of `f-values`
precludes this need, and allows the non-optimising equation solver HYBRD to use the same constraint
equations.

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

Enable all the relevant consistency equations, and it is advisable to enable the corresponding iterations variables shown first as corresponding itvs in the vardes.html file. A number of limit equations (inequality constraints) can also be activated. For limit equations, the corresponding f-value must be selected as an iteration variable. In optimisation more, the number of iteration variables is unlimited.

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