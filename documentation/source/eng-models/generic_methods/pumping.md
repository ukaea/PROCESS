# Pumping Methods

## Pumping coolant friction | `darcy_friction_haaland()`

 The pressure drop is based on the Darcy friction factor, using the [Haaland equation](https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation), an approximation to the implicit Colebrook–White equation. 

$$
\frac{1}{\sqrt{f}} = -1.8 \log_{10}{\left[ \left(\frac{\epsilon / D}{3.7}\right)^{1.11} + \frac{6.9}{\text{Re}} \right]}
$$

-------------------

## Reynolds number | `calculate_reynolds_number()`

$$
\mathrm{Re} = \frac{\rho v \left(2r_{\text{channel}}\right)}{\mu}
$$

where $\rho$ is the coolant density and $\mu$ is the coolant viscosity.

-------------------

## Gnielinski heat transfer | `gnielinski_heat_transfer_coefficient()`

1. **Calculate the Reynolds number:**

    $$
    \mathrm{Re} = \frac{\rho v \left(2r_{\text{channel}}\right)}{\mu}
    $$

    where $\rho$ is the coolant density and $\mu$ is the coolant viscosity.

2. **Calculate the Prandtl number:**

    $$
    \mathrm{Pr} = \frac{c_{\text{p}}\mu}{k}
    $$

    were $c_{\text{p}}$ is the coolant heat capacity and $k$ is the coolant thermal conductivity.

3. **Calculate the Darcy friction factor using the [`darcy_friction_haaland()`](#pumping-coolant-friction--darcy_friction_haaland) method:**

    $$
    f = \texttt{darcy_friction_haaland()}
    $$

4. **Calculate the Nusselt number using the [Gnielinski correlation](https://en.wikipedia.org/wiki/Nusselt_number#Gnielinski_correlation):**

    $$
    \mathrm{Nu_D}  = \frac{\left(f/8\right)\left(\mathrm{Re}-1000\right)\mathrm{Pr}}{1+12.7\left(f/8\right)^{0.5}\left(\mathrm{Pr}^{2/3}-1\right)}
    $$

    The relation is valid for:

    $$
    0.5 \le \mathrm{Pr} \le 2000 \\
    3000 \le \mathrm{Re} \le 5 \times 10^6
    $$

5. **Calculate the heat transfer coefficient with the Nusselt number:**

    $$
    h = \frac{\mathrm{Nu_D}k}{2r_{\text{channel}}}
    $$

-------------------------

## Pipe bend elbow coefficient | `elbow_coeff()`

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