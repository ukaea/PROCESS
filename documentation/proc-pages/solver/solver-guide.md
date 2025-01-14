
# Constraint Equations

Any computer program naturally contains many equations. The built-in equation 
solvers within PROCESS act on a special class, known as constraint equations, 
all of which are formulated in the source file `constraint equations.f90`. These 
can be split into two types:
 
**Consistency equations** -- that enforce consistency between the physics and 
engineering parameters

**limit equations** -- that enforce various parameters to lie within their allowed 
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
limits are not exceeded. Each of these equations has an associated `f-value`, which allow 
them to be coded as equalities.

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
2. Set the corresponding `f-value` `fpnetel = 1.0D0`
3. Ensure that `fpnetel` (iteration variable no. 25) **IS NOT** selected as an iteration variable.
4. Set `pnetelin` to the required net electric power.

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

## Non-optimisation mode

Non-optimisation mode is sometimes used to perform benchmark comparison, whereby the machine size, output power etc. are known and one only wishes to find the calculated stresses, beta values and fusion powers, for example.

Running `PROCESS` in non-optimisation mode requires few changes to be made to the input file from the optimisation case. The main differences between optimisation mode and non-optimisation mode are:

1. Non-optimisation mode does NOT apply lower or upper bounds to the iteration variables. It follows that limit equations are no enforced.

2. In non-optimisation mode the number of active iteration variables must be equal to the number of constraints.

3. A figure of merit is not available in non-optimisation mode.

4. Scans cannot be performed in non-optimisation mode.

As before, the user must decide which constraint equations and iteration variables to activate.