The plasma beta limit[^1] is given by 

$$\begin{aligned}
\beta < g \, \frac{I(\mbox{MA})}{a(\mbox{m}) \, B_0(\mbox{T})}
\end{aligned}$$

where $B_0$ is the axial vacuum toroidal field. The beta
coefficient $g$ is set using input parameter `dnbeta`. To apply the beta limit, 
constraint equation 24 should be turned on with iteration variable 36
(`fbetatry`). 

By default, $\beta$ is defined with respect to the total equilibrium B-field [^2]. 

| `iculbl` | Description |
| :-: | - |
| 0 (default) | Apply the $\beta$ limit to the total plasma beta (including the contribution from fast ions) |
| 1 | Apply the $\beta$ limit to only the thermal component of beta |
| 2 | Apply the $\beta$ limit to only the thermal plus neutral beam contributions to beta |
| 3 | Apply the $\beta$ limit to the total beta (including the contribution from fast ions), calculated using only the toroidal field |

### Scaling of beta $g$ coefficient

Switch `gtscale` determines how the beta $g$ coefficient `dnbeta` should 
be calculated, using the inverse aspect ratio $\epsilon = a/R$.

| `gtscale` | Description |
| :-: | - |
| 0 | `dnbeta` is an input. |
| 1 | $g=2.7(1+5\epsilon^{3.5})$ (which gives g = 3.0 for aspect ratio = 3) |
| 2 | $g=3.12+3.5\epsilon^{1.7}$ (based on Menard et al. "Fusion Nuclear Science Facilities and Pilot Plants Based on the Spherical Tokamak", Nucl. Fusion, 2016, 44)  |

!!! Note 
    `gtscale` is over-ridden if `iprofile` = 1.

### Limiting $\epsilon\beta_p$

To apply a limit to the value of $\epsilon\beta_p$, where $\epsilon = a/R$ is
the inverse aspect ratio and $\beta_p$ is the poloidal $\beta$, constraint equation no. 6 should be 
turned on with iteration variable no. 8 (`fbeta`). The limiting value of $\epsilon\beta_p$ 
is be set using input parameter `epbetmax`.

[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
[^2]: D.J. Ward, 'PROCESS Fast Alpha Pressure', Work File Note F/PL/PJK/PROCESS/CODE/050