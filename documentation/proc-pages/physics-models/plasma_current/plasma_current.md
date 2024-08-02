# Plasma Current 

## Overview







## Plasma Current Calculation | `calculate_plasma_current()`

A number of plasma current scaling laws are available in PROCESS [^1]. These are calculated in 
routine `calculate_plasma_current()`, in `physics.py`. The safety factor $q_{95}$ required to prevent 
disruptive MHD instabilities dictates the plasma current $I_{\text{p}}$:

$$\begin{aligned}
I_{\text{p}} = f_q \frac{2\pi}{\mu_0}  \frac{a^2 B_{\text{T}}}{R \ q_{95}}
\end{aligned}$$

The factor $f_q$ makes allowance for toroidal effects and plasma shaping (elongation and 
triangularity). Several formulae for this factor are available depending on the value of 
the switch `i_plasma_current`, as follows:

---------------

### Peng analytic fit

Switch value: `i_plasma_current = 1`

The formula for calculating `fq` is:

$$f_q = \left(\frac{{1.22 - 0.68  \epsilon}}{{(1.0 - \epsilon^2)^2}}\right)  \mathtt{{sf}}^2$$

Where $\epsilon$ is the inverse [aspect ratio](../plasma_geometry.md) ($\mathtt{eps}$) and $\mathtt{sf}$ is the shaping factor calculated in the [poloidal perimeter](../plasma_geometry.md#poloidal-perimeter) function in `plasma_geometry.py`

-----------

### 2.Peng double null divertor scaling (ST)
[^4]

### 3. Simple ITER scaling

### 4. Revised ITER scaling 
[^5]

### 5. Todd empirical scaling, I

### 6. Todd empirical scaling, II

### 7. Connor-Hastie model

### 8. Sauter model, allows negative $\delta$

### 9. Scaling for spherical tokamaks, based on a fit to a family of equilibria derived by Fiesta: 

## Plasma Current Profile Consistency

A limited degree of self-consistency between the plasma current profile and other parameters [^6] can be 
enforced by setting switch `iprofile = 1`. This sets the current 
profile peaking factor $\alpha_J$ (`alphaj`),  the normalised internal inductance $l_i$ (`rli`) and beta limit $g$-factor (`dnbeta`) using the 
safety factor on axis `q0` and the cylindrical safety factor $q*$ (`qstar`):   

$$\begin{aligned}
\alpha_J = \frac{q*}{q_0} - 1
\end{aligned}$$

$$\begin{aligned}
l_i = \rm{ln}(1.65+0.89\alpha_J)
\end{aligned}$$

$$\begin{aligned}
g = 4 l_i
\end{aligned}$$

It is recommended that current scaling law `i_plasma_current = 4` is used if `iprofile = 1`. 
This relation is only applicable to large aspect ratio tokamaks.

For spherical tokamaks, the normalised internal inductance can be set from the elongation using `iprofile = 4` or `iprofile = 5`:

$$\begin{aligned}
l_i = 3.4 - \kappa_x
\end{aligned}$$

Further desciption of `iprofile` is given in [Beta Limit](./plasma_beta.md).

[^1]: D.J. Ward, 'PROCESS Fast Alpha Pressure', Work File Note F/PL/PJK/PROCESS/CODE/050
[^2]:  Albajar, Nuclear Fusion **41** (2001) 665
[^3]: M. Kovari, R. Kemp, H. Lux, P. Knight, J. Morris, D.J. Ward, '“PROCESS”: A systems code for fusion power plants—Part 1: Physics' Fusion Engineering and Design 89 (2014) 3054–3069
[^4]: J.D. Galambos, 'STAR Code : Spherical Tokamak Analysis and Reactor Code',
Unpublished internal Oak Ridge document.
[^5]: W.M. Nevins, 'Summary Report: ITER Specialists' Meeting on Heating and
Current Drive', ITER-TN-PH-8-4, 13--17 June 1988, Garching, FRG
[^6]: Y. Sakamoto, 'Recent progress in vertical stability analysis in JA',
Task meeting EU-JA #16, Fusion for Energy, Garching, 24--25 June 2014