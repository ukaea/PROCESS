# Culham Lower Hybrid | `cullhy()`

- `i_hcd_primary/i_hcd_secondary` = 6

This routine calculates the current drive parameters for a
lower hybrid system, based on the AEA FUS 172 model.
AEA FUS 251: A User's Guide to the PROCESS Systems Code
AEA FUS 172: Physics Assessment for the European Reactor Study[^1]

1. Call the [`lhrad()`](#lower-hybrid-wave-absorption-radius--lhrad) method to calculate the lower hybrid wave absorption radius, `rratio`.
2. Calculate the penetration radius, `rpenet`, by multiplying `rratio` with the minor radius of the plasma.
3. Calculate the local density, `dlocal`, using the `neprofile()` function from the `profiles_module` module. This function takes into account various plasma parameters such as the density profile, electron density at the edge, pedestal density, separatrix density, and the value of the parameter `alphan`.
4. Similarly, calculate the local temperature, `tlocal`, using the `teprofile()` function from the `profiles_module` module. This function considers parameters such as the temperature profile, electron temperature at the edge, pedestal temperature, separatrix temperature, `alphat`, and `tbeta`.
5. Calculate the local toroidal magnetic field, `blocal`, using the formula `b_plasma_toroidal_on_axis * rmajor / (rmajor - rpenet)`. Here, `b_plasma_toroidal_on_axis` is the toroidal magnetic field at the magnetic axis, and `rmajor` is the major radius of the plasma.
6. Calculate the parallel refractive index, `nplacc`, which is needed for plasma access. It uses the local density `dlocal` and the local magnetic field `blocal` to calculate a fraction `frac`. `nplacc` is then obtained by adding `frac` to the square root of `1.0 + frac * frac`.
7. Calculate the local inverse aspect ratio, `epslh`, by dividing `rpenet` by `rmajor`.
8. Calculate the LH normalised efficiency, `x`, using the formula `24.0 / (nplacc * sqrt(tlocal))`.
9. Calculate several intermediate terms, `term01`, `term02`, `term03`, and `term04`, using different formulas involving `nplacc`, `physics_variables.zeff`, `tlocal`, `epslh`, and `x`.
10. Calculate the current drive efficiency, `gamlh`, using the formula `term01 * term02 * (1.0e0 - term03 / term04)`.
11. Return the current drive efficiency normalised by the product of `0.1e0 * dlocal` and `physics_variables.rmajor`.

[^1]: T. C. Hender, M. K. Bevir, M. Cox, R. J. Hastie, P. J. Knight, C. N. Lashmore-Davies, B. Lloyd, G. P. Maddison, A. W. Morris, M. R. O'Brien, M.F. Turner abd H. R. Wilson, *"Physics Assessment for the European Reactor Study"*, AEA Fusion Report AEA FUS 172 (1992)

## Lower Hybrid wave absorption radius | `lhrad`()

This routine determines numerically the minor radius at which the damping of Lower Hybrid waves occurs, using a Newton-Raphson method to establish the correct minor radius ratio. The required minor radius ratio has been found when the difference between the results of the two formulae for the energy E given in AEA FUS 172, p.58 is sufficiently close to zero.

Correction to refractive index (kept within valid bounds)
  $\mathtt{drfind} = \min\left(0.7, \max\left(0.1, \frac{12.5}{\text{te0}}\right)\right)$

Use Newton-Raphson method to establish the correct minor radius ratio. The required minor radius ratio has been found when the difference between the results of the two formulae for the energy E given in AEA FUS 172, p.58 is sufficiently close to zero.

Iterate over the following steps to find the minor radius ratio:

1. Set an initial guess for the minor radius ratio, $\mathtt{rat0}$, to 0.8.
2. Repeat the following steps for a maximum of 100 iterations:
    - Calculate the minor radius ratios, $r1$ and $r2$, by subtracting and adding 0.1% of $\mathtt{rat0}$, respectively.
    - Evaluate the function $g$ at $\mathtt{rat0}$, $r1$, and $r2$ using the method `lheval(drfind, rat)`.
    - Calculate the gradient of $g$ with respect to the minor radius ratio, $\frac{{dg}}{{dr}}$, using the formula $\frac{{g2 - g1}}{{r2 - r1}}$.
    - Calculate a new approximation for the minor radius ratio, $\mathtt{rat1}$, using the formula $\mathtt{rat0} - \frac{{g0}}{{\frac{{dg}}{{dr}}}}$.
    - Ensure that $\mathtt{rat1}$ is within the bounds of 0.0001 and 0.9999.
    - If the absolute value of $g0$ is less than or equal to 0.01, exit the loop.
    - Update $\mathtt{rat0}$ with the new approximation, $\mathtt{rat1}$.
3. If the loop completes all 100 iterations without finding a satisfactory solution, report an error and set $\mathtt{rat0}$ to 0.8.
4. Return the final value of $\mathtt{rat0}$.
