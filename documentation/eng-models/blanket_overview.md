# Blanket

## General blanket methods | `BlanketLibrary`


--------------------

### Coolant mechanical pumping power | `coolant_pumping_power()`

To calculate the coolant pumping power we use the change in enthalpies of the coolant as it goes through the pump. 
**We assume the pump is isentropic so the entropy change of the coolant is 0**. 

The mechanical pumping power is defined as:

$$
P = \frac{\frac{\dot{m} \times \left(H_{\text{out}}-H_{\text{in}}\right)}{\eta}}{\left(1-fp\right)}
$$

where $\dot{m}$ is the coolant mass flow rate, $H$ is the coolant enthalpy, $\eta$ is the isentropic efficiency of the pump and $\gamma$ is the adiabatic index of the coolant.

$$
fp = \frac{T_{\text{pump,out}}\left(\frac{P_{\text{pump,out}}}{P_{\text{pump,in}}}\right)^{-\frac{\gamma -1}{\gamma}}}{\eta \left(T_{\text{pump,in}}-T_{\text{pump,out}}\right)}
$$

------------------

### Coolant pressure drop | `coolant_friction_pressure_drop()`

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

-------------------

### Pipe bend elbow coefficient | `elbow_coeff()`

This function calculates the elbow bend coefficients for pressure drop calculations.

$$
a = 1.0 \quad \text{if} \ \theta = 90^{\circ} \\
a = 0.9 \times \sin{\left(\frac{\theta \pi}{180^{\circ}}\right)} \quad \text{if} \ \theta < 70^{\circ} \\
a = 0.7 + 0.35 \times \sin{\left(\frac{\theta}{90^{\circ}} \times \frac{\pi}{180^{\circ}}\right)} \quad \text{if} \ \theta > 90^{\circ} \\
$$

where $\theta$ is the angle of the pipe bend.

$$
b = \frac{0.21}{\sqrt{\frac{R_{\text{elbow}}}{D_{\text{pipe}}}}}\quad \text{if} \ \frac{R_{\text{elbow}}}{D_{\text{pipe}}} \ge 1 \\
b = \frac{0.21}{\left(\frac{R_{\text{elbow}}}{D_{\text{pipe}}}\right)^{2.5}}\quad \text{if} \ \frac{R_{\text{elbow}}}{D_{\text{pipe}}} \le 1 \\
\text{else} \quad b =0.21
$$

The elbow coefficient is given by:

$$
ab + \left( f_{\text{D}} \times \frac{R_{\text{elbow}}}{D_{\text{pipe}}}\right) \times \theta \times \left(\frac{\pi}{180^{\circ}}\right)
$$

