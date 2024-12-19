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

The calculation of the beta component given by neutral beams is calculated in the neutral beam fusion calculations in [`beam_fusion()`](../fusion_reactions/beam_reactions.md#beam-slowing-down-properties--beam_fusion). A description can be found [here](../fusion_reactions/beam_reactions.md#derivation-of-beam-slowing-down-rate-and-critical-energy).

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

## Troyon Beta Limit

The Troyon plasma beta limit is given by[^0][^1]:

$$\begin{aligned}
\beta < 0.01\, g \, \frac{I \ [\mbox{MA}]}{a \ [\mbox{m}] \, B_0 \ [\mbox{T}]}
\end{aligned}$$

where $B_0$ is the axial vacuum toroidal field. The beta
coefficient $g$ is set using input parameter `beta_norm_max`. To apply the Troyon limit please see [Beta upper limit](#beta-upper-limit).

By default, $\beta$ is defined with respect to the total equilibrium B-field. This can be changed depending on the setting of `i_beta_component` seen below:

| `i_beta_component` | Description |
| :-: | - |
| 0 (default) | Apply the $\beta$ limit to the total plasma beta (including the contribution from fast alphas and  neutral beams) |
| 1 | Apply the $\beta$ limit to only the thermal component of beta |
| 2 | Apply the $\beta$ limit to only the thermal plus neutral beam contributions to beta |
| 3 | Apply the $\beta$ limit to the total toroidal beta (including the contribution from fast alphas and neutral beams) |

------------

### Setting the Beta $g$ Coefficient

Switch `iprofile` determines how the beta $g$ coefficient `beta_norm_max` should
be calculated. The following switch options are available below:

#### User Input

This can be activated by stating `iprofile = 0` in the input file.

`alphaj`, `rli` and `beta_norm_max` are inputs.

---------

#### Wesson Relation

This can be activated by stating `iprofile = 1` in the input file.

`alphaj`, `rli` and `beta_norm_max` are calculated consistently. 

`beta_norm_max` is calculated using:  

$$
g = 4l_i
$$

This relation is based off of data taken from DIII-D shots[^7].

This is only recommended for high aspect ratio tokamaks[^3].

---------

#### Original Scaling Law

This can be activated by stating `iprofile = 2` in the input file.

`alphaj` and `rli` are inputs. `beta_norm_max` calculated using $g=2.7(1+5\epsilon^{3.5})$ (which gives g = 3.0 for aspect ratio = 3)

---------

#### Menard Beta Relation

This can be activated by stating `iprofile = 3` in the input file.

 `alphaj` and `rli` are inputs. `beta_norm_max` calculated using $g=3.12+3.5\epsilon^{1.7}$ [^4]

---------

#### Menard Inductance Relation

This can be activated by stating `iprofile = 4` in the input file.

`alphaj` and `beta_norm_max` are inputs. `rli` calculated from elongation [^4]. This is only recommended for spherical tokamaks.

---------

#### Menard Beta & Inductance Relation

This can be activated by stating `iprofile = 5` in the input file.

`alphaj` is an input.  `rli` calculated from elongation and `beta_norm_max` calculated using $g=3.12+3.5\epsilon^{1.7}$ [^4]. This is only recommended for spherical tokamaks.

---------

#### Tholerus Relation

This can be activated by stating `iprofile = 6` in the input file.

`alphaj` and `c_beta` are inputs.  `rli` calculated from elongation and `beta_norm_max` calculated using 

$$
C_{\beta}\approx\frac{(g-3.7)F_p}{12.5-3.5 F_p}
$$

where $F_p$ is the pressure peaking, $F_p = p_{\text{ax}} / \langle p \rangle$ and $C_{\beta}$ is the destabilization parameter (default 0.5)[^5].  

<u> This is only recommended for spherical tokamaks </u>

---------

Further details on the calculation of `alphaj` and `rli` is given in [Plasma Current](./plasma_current.md).

----------------------

## Key Constraints

### Beta consistency

This constraint can be activated by stating `icc = 1` in the input file.

Ensures the relationship between $\beta$, density, temperature and total magnetic field is withheld by checking the fixed input or iteration variable $\mathtt{beta}$ is consistent in value with the rest of the physics parameters

$$
\mathtt{beta} \equiv \frac{2\mu_0 \langle n_{\text{e}}T_{\text{e}}+n_{\text{i}}T_{\text{i}}\rangle}{B^2} + \beta_{\alpha} + \beta_{\text{beam}}
$$

**It is highly recommended to always have this constraint on as it is a global consistency checker**

----------------

### Poloidal beta and inverse aspect upper limit

This constraint can be activated by stating `icc = 6` in the input file [^6].

The limiting value of $\epsilon\beta_p$ is be set using input parameter `beta_poloidal_eps_max`.

The scaling value `fbeta_poloidal_eps` can be varied also.

!!! note "Origin of the $\epsilon\beta_p$ limit"

    High poloidal beta shots in TFTR were performed[^6] and it was found that as $\beta_p$,
    exceeds approximately 1.2 times the aspect ratio, a separatrix with
    an inside poloidal field null is observed to limit the outer boundary
    of the plasma. Since the curvature of TFTR’s applied vertical field
    is constant, the appearance of the poloidal field null corresponds to
    the equilibrium poloidal beta limit.

--------------------

### Beta upper limit

This constraint can be activated by stating `icc = 24` in the input file.

It is the general setting of the $\beta$ limit depending on the $\beta_{\text{N}}$ value calculated in the [beta limit](#beta-limit) calculations.

The upper limit value of beta is calculated by `calculate_beta_limit()`. The beta
coefficient $g$ can be set using `beta_norm_max`, depending on the setting of [`iprofile`](#setting-the-beta--coefficient). It can be set directly or follow some relation.

The scaling value `fbeta_max` can be varied also.

**It is recommended to have this constraint on as it is a plasma stability model**

--------------------

### Poloidal upper limit

This constraint can be activated by stating `icc = 48` in the input file.

The value of `beta_poloidal_max` can be set to the desired maximum poloidal beta. The scaling value `fbeta_poloidal` can be varied also.

-------------------

### Beta lower limit

This constraint can be activated by stating `icc = 84` in the input file.

The value of `beta_max` can be set to the desired minimum total beta. The scaling value `fbeta_min` can be varied also.

[^0]: F. Troyon et.al,  “Beta limit in tokamaks. Experimental and computational status,” Plasma Physics and Controlled Fusion, vol. 30, no. 11, pp. 1597–1609, Oct. 1988, doi: https://doi.org/10.1088/0741-3335/30/11/019.

[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',

[^2]: D.J. Ward, 'PROCESS Fast Alpha Pressure', Work File Note F/PL/PJK/PROCESS/CODE/050

[^3]: Tokamaks 4th Edition, Wesson, page 116

[^4]: J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,” Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016, doi: https://doi.org/10.1088/0029-5515/56/10/106023.

[^5]: E. Tholerus et al., “Flat-top plasma operational space of the STEP power plant,” Nuclear Fusion, Aug. 2024, doi: https://doi.org/10.1088/1741-4326/ad6ea2.
‌
[^6]: M. E. Mauel et al., “Operation at the tokamak equilibrium poloidal beta-limit in TFTR,” Nuclear Fusion, vol. 32, no. 8, pp. 1468–1473, Aug. 1992. doi:https://dx.doi.org/10.1088/0029-5515/32/8/I14

[^7]: T. T. S et al., “Profile Optimization and High Beta Discharges and Stability of High Elongation Plasmas in the DIII-D Tokamak,” Osti.gov, Oct. 1990. https://www.osti.gov/biblio/6194284 (accessed Dec. 19, 2024).
‌
