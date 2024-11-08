# Beta Limit

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

## Setting the Beta $g$ Coefficient

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

## Limiting $\epsilon\beta_p$

To apply a limit to the value of $\epsilon\beta_p$, where $\epsilon = a/R$ is
the inverse aspect ratio and $\beta_p$ is the poloidal $\beta$, constraint equation no. 6 should be 
turned on with iteration variable no. 8 (`fbeta`). The limiting value of $\epsilon\beta_p$ 
is be set using input parameter `epbetmax`.

## Key Constraints



[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',

[^2]: D.J. Ward, 'PROCESS Fast Alpha Pressure', Work File Note F/PL/PJK/PROCESS/CODE/050

[^3]: Tokamaks 4th Edition, Wesson, page 116

[^4]: Menard et al. (2016), Nuclear Fusion, 56, 106023

[^5]: Tholerus et al. (2024), arXiv:2403.09460
