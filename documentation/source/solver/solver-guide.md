# Solver PROCESS docs

## Constraint Equations

PROCESS models a fusion power plant at a **self-consistent operating point**. It contains a large set of coupled physics and engineering models that describe different aspects of the plant. For a design to constitute a PROCESS solution, quantities calculated by different parts of the code must satisfy the required relationships between them.

PROCESS represents these relationships using **constraint equations**, which are defined in `process/core/solver/constraints.py`.

Constraint equations fall into two classes:

- **Equality constraints (consistency equations)** define relationships that must be satisfied for the calculated plant to be self-consistent.
- **Inequality constraints (limit equations)** define physics or engineering limits that determine whether a solution is feasible.

PROCESS uses user-defined **iteration variables** as the degrees of freedom available to the solver when solving the selected constraint problem.

## PROCESS Solutions

At the heart of a PROCESS solution is a plasma operating in a self-consistent equilibrium state. Its state cannot be chosen arbitrarily: plasma temperature, density, heating, fusion reactions, radiation, transport losses and fuelling are all coupled.

The principal plasma state is described by the electron temperature and density,

$$
T_e,\qquad n_e.
$$

To determine these two quantities, two independent governing equations are required.

### Governing plasma equations

In PROCESS these are:

1. **Global plasma power balance - constraint 2**
2. **Fuel-ion equilibrium - constraint 93**

Together these equations close the plasma system and define the equilibrium plasma state solved by PROCESS.

For a given set of design parameters, a plasma solution is therefore a pair

$$
(T_e,n_e)
$$

for which both governing equations are satisfied simultaneously.

#### Global plasma power balance

For the plasma to remain in a steady operating state, the power supplied to the plasma must balance the power lost from it.

Schematically,

$$
P_{\mathrm{heating}} = P_{\mathrm{loss}}.
$$

The heating and loss terms depend on the plasma state and on the selected PROCESS physics models. The power balance therefore provides one relationship between plasma temperature and density.

PROCESS represents this condition using **constraint 2**.

#### Fuel-ion equilibrium

The plasma must also satisfy particle balance for the fuel ions.

Fusion reactions consume fuel ions. To maintain the plasma state, this consumption must be consistent with the calculated fuelling rate and fuel burn-up.

Schematically,

$$
\text{fuel supplied} = \text{fuel required to sustain the fusion reaction rate}.
$$

PROCESS represents this condition using **constraint 93**.

Together, global plasma power balance and fuel-ion equilibrium form a coupled system in $T_e$ and $n_e$, whose simultaneous solution determines the self-consistent plasma state.

### Other consistency equations

The governing plasma equations are not necessarily the only equality constraints required for a complete PROCESS calculation.

Other physics and engineering models may introduce additional consistency relationships, for example for the machine radial build.

A **PROCESS solution** is therefore a point at which the governing plasma equations and all other required equality constraints are satisfied.

### Feasible solutions

A self-consistent PROCESS solution is not necessarily feasible.

PROCESS therefore also defines **inequality constraints**, or **limit equations**, which represent physics and engineering limits on the design.

Examples include limits on plasma density, beta, fusion power, neutron wall load, divertor power loading, magnet stresses and burn time.

This gives an important distinction:

$$
\text{PROCESS solution}
=
\text{required equality constraints satisfied}
$$

whereas

$$
\text{feasible PROCESS solution}
=
\text{solution that also satisfies the required inequality constraints}
$$

A point may therefore be a valid solution of the PROCESS governing equations while still being infeasible because one or more design limits are violated.

## Specifying constraint equations

Constraint equations are selected using the `icc` array. For example:

```text
n_equality_constraints = 3

* Equalities
icc = 2   * Global plasma power balance
icc = 11  * Radial build
icc = 93  * Fuel-ion equilibrium

* Inequalities
icc = 5   * Density upper limit
icc = 9   * Fusion power upper limit
icc = 15  * L-H power threshold limit
icc = 24  * Beta upper limit
```

Each `icc = n` statement activates constraint `n`.

`n_equality_constraints` specifies how many entries at the beginning of the `icc` array are treated as equality constraints. All equality constraints must therefore be listed before any inequality constraints.

A list of the available constraints and their corresponding names can be found [here](../../source/reference/process/data_structure/numerics/#process.data_structure.numerics.lablcc).

For a typical equality constraint,

$$
g(\mathbf{x}) = h(\mathbf{x}),
$$

PROCESS forms the normalised residual

$$
c_i = 1 - \frac{g}{h},
$$

which is satisfied when

$$
c_i = 0.
$$

For inequality constraints, the residual is formulated so that the permitted region satisfies

$$
c_i \geq 0.
$$

## Iteration Variables

Iteration variables are the quantities that the solver is allowed to vary in order to satisfy the active constraints.

They are selected using the `ixc` array, for example:

```text
ixc = 4  * temp_plasma_electron_vol_avg_kev
ixc = 6  * nd_plasma_electrons_vol_avg
```

The equations are coupled, so iteration variables should not generally be thought of as corresponding one-to-one with individual constraints. Changing one iteration variable may affect several constraint residuals.

## Figure of Merit

In optimisation mode there may be many feasible solutions. PROCESS selects between them using a **figure of merit**, or objective function.

The switch `i_figure_merit` selects the objective. A positive value indicates minimisation, while a negative value indicates maximisation.

The optimisation problem can therefore be viewed as finding the feasible PROCESS solution that gives the best value of the selected figure of merit.

## Convergence

PROCESS solves a nonlinear, coupled constrained problem iteratively. Numerical convergence indicates that the solver has satisfied its convergence criteria. The constraint residuals should then be inspected to determine whether the required consistency equations and limits have been satisfied.

For equality constraints,

$$
c_i \approx 0.
$$

Numerical convergence should not be confused with feasibility: a converged calculation may still violate one or more inequality constraints.

## Optimisation mode

Optimisation mode is selected using

```text
i_process_run_mode = 1
```

In this mode, PROCESS varies the active iteration variables in order to:

1. satisfy the equality constraints;
2. satisfy the inequality constraints; and
3. optimise the selected figure of merit.

The number of iteration variables may exceed the number of equality constraints because the additional degrees of freedom allow the optimiser to explore the feasible design space.

## Evaluation mode

Evaluation mode is selected using

```text
i_process_run_mode = -2
```

It is intended for evaluating a specified design while maintaining model consistency, rather than searching for an optimum.

In evaluation mode:

1. equality constraints are solved;
2. inequality constraints may be reported but are not enforced;
3. the number of iteration variables must equal the number of equality constraints;
4. no figure of merit is optimised; and
5. iteration-variable optimisation bounds are not applied.

For example, for a specified design point the plasma state can be obtained by solving the two governing plasma equations for $T_e$ and $n_e$:

```text
i_process_run_mode = -2

n_equality_constraints = 2

* Equalities
icc = 2   * Global plasma power balance
icc = 93  * Fuel-ion equilibrium

* Inequalities
* Reported but not enforced
icc = 5   * Density upper limit
icc = 8   * Neutron wall load upper limit
icc = 9   * Fusion power upper limit

* Iteration variables
ixc = 4   * temp_plasma_electron_vol_avg_kev
ixc = 6   * nd_plasma_electrons_vol_avg
```

This is useful for benchmark calculations or for evaluating a prescribed machine design without performing an optimisation.