# Plasma Beta

## Overview

The efficiency of confinement of plasma pressure by the magnetic field is represented by the ratio:

$$
\beta = \frac{2\mu_0p}{B^2}
$$

There are several different measures of this type, arising from different choices of definition and from the need to quantify different equilibrium properties.

In its expanded form of magnetic field components:

$$
\beta = \frac{2\mu_0p}{B^2_{\text{T}}+B^2_{\text{P}}}
$$

This is often broken down into separate related quantities knows as the toroidal $\beta_{\text{t}}$ and the poloidal beta $\beta_{\text{p}}$ whose definitions are:

$$
\beta_{\text{t}} = \frac{2\mu_0p}{B^2_{\text{T}}} \\
\beta_{\text{p}} = \frac{2\mu_0p}{B^2_{\text{P}}}
$$

The relationship between these quantities can be written as:

$$
\frac{1}{\beta} = \frac{1}{\beta_{\text{t}}} + \frac{1}{\beta_{\text{p}}}
$$

implying that the total $\beta$ is dominated by the smaller of the two contributions.Observe that by definition $\beta \le 1$. However,either $\beta_{\text{t}}$ or $\beta_{\text{p}}$, but not both, can be greater than unity.

The above in its simplest form is known as the total thermal beta $\beta_{\text{th}}$ as it only takes into account the pressure produced by the ions and electrons such that:

$$
\beta_{\text{th}} = \frac{2\mu_0 \langle n_{\text{e}}T_{\text{e}}+n_{\text{i}}T_{\text{i}}\rangle}{B^2}
$$

In future reactors there will be a notable pressure contribution due to the productions of fast alpha particles and that from plasma-beam interactions if used. So the true total beta can be defined as:

$$
\beta_{\text{tot}} = \frac{2\mu_0 \langle n_{\text{e}}T_{\text{e}}+n_{\text{i}}T_{\text{i}} + p_{\alpha}+ p_{\text{beam}}\rangle}{B^2} \\
\beta_{\text{tot}} = \beta_{\text{th}}+\beta_{\alpha}+ \beta_{\text{beam}}
$$

Models for the fast alpha particle pressure contribution can be found [here](plasma_alpha_beta_contribution.md).

----------------

## Derivation of plasma beta parameter

From the lowest order momentum balance equation:

$$
\nabla p = \mathbf{j} \times \mathbf{B}
$$

as it pertains to plasma equilibrium. The current and the field must also satisfy Maxwell’s equations

$$
\mu_0\mathbf{j} = \nabla \times \mathbf{B},  \quad \nabla \cdot \mathbf{B} = 0
$$

A consequence of the above is that the current and  magnetic field lie on isobaric surfaces. This follows directly from the observation that $\mathbf{j} \cdot \nabla p = \mathbf{B} \cdot \nabla p = 0$ and $\nabla p$ is everywhere normal to the surface $p = \text{const}$. The force is everywhere normal to the isobaric surface and just balances the pressure gradient force, $-\nabla p$. Although the current and the magnetic field lie in a common flux surface, they can only be parallel in regions where the pressure gradient vanishes. Since currents that are parallel to the field do not contribute to the $\mathbf{j} \times \mathbf{B}$ force, they are known as “force-free currents.”

Taking the divergence of $\mu_0\mathbf{j} = \nabla \times \mathbf{B}$  yields $\nabla \cdot \mathbf{j} = 0$ which is consistent with the quasineutrality assumption. The vector product of $\mathbf{B}$ with the momentum balance equation leads to an expression for the current perpendicular to $\mathbf{B}$,

$$
\mathbf{j}_{\perp} = \frac{\mathbf{B} \times \nabla p}{B^2}
$$

We see that there can be no perpendicular current in the absence of a pressure gradient. This leads to a momentum balance expressed in terms of the pressure and magnetic field.

$$
\nabla \left(p+\frac{B^2}{2\mu_0}\right) = \frac{1}{\mu_0}\left(\mathbf{B} \cdot \nabla \right)\mathbf{B}
$$

The right hand side vanishes when the field lines are straight and parallel, in which case above reduces to a simple statement that the total (kinetic plus magnetic) pressure is constant everywhere within a confined plasma, 

$$
p + \frac{B^2}{2\mu_0} \equiv \frac{B^2}{2\mu_0}\left(1+\beta \right) = \text{const}
$$

where we have introduced the quantity

$$
\beta \equiv \frac{2\mu_0p}{B^2}
$$

------------------------

## Beta Limit

The plasma beta limit[^1] is given by 

$$\begin{aligned}
\beta < 0.01\, g \, \frac{I(\mbox{MA})}{a(\mbox{m}) \, B_0(\mbox{T})}
\end{aligned}$$

where $B_0$ is the axial vacuum toroidal field. The beta
coefficient $g$ is set using input parameter `dnbeta`. To apply the beta limit, 
constraint equation 24 should be turned on with iteration variable 36
(`fbetatry`). 

By default, $\beta$ is defined with respect to the total equilibrium B-field [^2]. 

| `i_beta_component` | Description |
| :-: | - |
| 0 (default) | Apply the $\beta$ limit to the total plasma beta (including the contribution from fast ions) |
| 1 | Apply the $\beta$ limit to only the thermal component of beta |
| 2 | Apply the $\beta$ limit to only the thermal plus neutral beam contributions to beta |
| 3 | Apply the $\beta$ limit to the total beta (including the contribution from fast ions), calculated using only the toroidal field |

### Setting the Beta $g$ Coefficient

Switch `iprofile` determines how the beta $g$ coefficient `dnbeta` should 
be calculated.

| `iprofile` | Description |
| :-: | - |
| 0 | `alphaj`, `rli` and `dnbeta` are inputs. |
| 1 (default) | `alphaj`, `rli` and `dnbeta` are calulcated consistently. `dnbeta` calculated using $g=4l_i$ [^3].  This is only recommended for high aspect ratio tokamaks.|
| 2 | `alphaj` and `rli` are inputs. `dnbeta` calculated using $g=2.7(1+5\epsilon^{3.5})$ (which gives g = 3.0 for aspect ratio = 3) |
| 3 | `alphaj` and `rli` are inputs. `dnbeta` calculated using $g=3.12+3.5\epsilon^{1.7}$ [^4]|
| 4 | `alphaj` and `dnbeta` are inputs. `rli` calculated from elongation [^4]. This is only recommended for spherical tokamaks.|
| 5 | `alphaj` is an input.  `rli` calculated from elongation and `dnbeta` calculated using $g=3.12+3.5\epsilon^{1.7}$ [^4]. This is only recommended for spherical tokamaks.|
| 6 | `alphaj` and `c_beta` are inputs.  `rli` calculated from elongation and `dnbeta` calculated using $C_{\beta}=(g-3.7)F_p / 12.5-3.5F_p$, where $F_p$ is the pressure peaking and $C_{\beta}$ is the destabilisation papermeter (default 0.5)[^5]. See Section 2.4 of Tholerus et al. (2024) for a more detailed description.  <u> This is only recommended for spherical tokamaks <u>.|

Further details on the calculation of `alphaj` and `rli` is given in [Plasma Current](./plasma_current.md).

----------------------

## Key Constraints

### Beta consistency

This constraint can be activated by stating `icc = 1` in the input file.

Relationship between beta, temperature (keV) and density

**It is highly recommended to always have this constraint on as it is a global consistency checker**

----------------

### Poloidal beta and inverse aspect upper limit

This constraint can be activated by stating `icc = 6` in the input file.

To apply a limit to the value of $\epsilon\beta_p$, where $\epsilon = a/R$ is
the inverse aspect ratio and $\beta_p$ is the poloidal $\beta$, constraint equation no. 6 should be 
turned on with iteration variable no. 8 (`fbeta`). The limiting value of $\epsilon\beta_p$ 
is be set using input parameter `epbetmax`.

--------------------

### Beta upper limit

This constraint can be activated by stating `icc = 24` in the input file.

--------------------

### Poloidal upper limit

This constraint can be activated by stating `icc = 48` in the input file.

-------------------

### Beta lower limit

This constraint can be activated by stating `icc = 84` in the input file.

[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',

[^2]: D.J. Ward, 'PROCESS Fast Alpha Pressure', Work File Note F/PL/PJK/PROCESS/CODE/050

[^3]: Tokamaks 4th Edition, Wesson, page 116

[^4]: Menard et al. (2016), Nuclear Fusion, 56, 106023

[^5]: Tholerus et al. (2024), arXiv:2403.09460
