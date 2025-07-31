# Blanket

## General Blanket methods | `BlanketLibrary`

### Coolant pressure drop | `pressure_drop()`

The pressure drop in the coolant is given by the [Darcy-Weisbach Equation](https://en.wikipedia.org/wiki/Darcy%E2%80%93Weisbach_equation)

For a cylindrical pipe of uniform diameter the pressure loss due to viscous effects can be characterized by:

$$
\Delta P = L\left[f_{\text{D}}\frac{\rho}{2}\frac{\langle v \rangle^2}{D_{\text{H}}}\right]
$$

where $L$ is the pipe length, $f_{\text{D}}$ is the [Darcy friction factor](https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae), $\rho$ is the coolant density, $\langle v \rangle$ is the mean flow coolant velocity and $D_{\text{H}}$ is the hydraulic diameter or the pipe diameter in this case.

To find the Darcy friction factor we need to know the Reynolds number given by:

$$
\text{Re} = \frac{\rho v L}{\mu}
$$

here $L$ is the characteristic length which we set to be the pipe diameter and $\mu$ is the coolant dynamic viscosity.

Using the Reynolds number we calculate the Darcy friction factor using the Haaland approximation calculated by `darcy_friction_haaland()`.

For the radius of the pipe bend we assume it to be 3 times the radius of the coolant channel.

The elbow coefficients for the 90 and 180 degree bends $\left(f_{\text{90,elbow}}, f_{\text{180,elbow}}\right)$ are clalculated via [`elbow_coeff()`](#pipe-bend-elbow-coefficient--elbow_coeff).

The pressure drop for the straights along the entire pipe length is the same as above:

$$
\Delta P = L\left[f_{\text{D}}\frac{\rho}{2}\frac{\langle v \rangle^2}{D_{\text{H}}}\right]
$$

where we define $\frac{f_{\text{D}}L}{D_{\text{H}}}$ as our straight section coefficient.

The pressure drop for the 90 and 180 degree bends are:

$$
\Delta P = N_{\text{90}} \left[f_{\text{90,elbow}} \frac{\rho \langle v \rangle^2}{2}\right]
$$

$$
\Delta P = N_{\text{180}} \left[f_{\text{180,elbow}} \frac{\rho \langle v \rangle^2}{2}\right]
$$

where $N_{\text{90}}$ and $N_{\text{180}}$ are the number of 90 and 180 degree bends in the system.

The total returned pressure drop is simply:

$$
\Delta P = L\left[f_{\text{D}}\frac{\rho}{2}\frac{\langle v \rangle^2}{D_{\text{H}}}\right] + N_{\text{90}} \left[f_{\text{90,elbow}} \frac{\rho \langle v \rangle^2}{2}\right] + N_{\text{180}} \left[f_{\text{180,elbow}} \frac{\rho \langle v \rangle^2}{2}\right]
$$


### Pipe bend elbow coefficient | `elbow_coeff()`