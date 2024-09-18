# Plasma Current 

## Overview







## Plasma Current Calculation | `calculate_plasma_current()`

This function calculates the plasma current shaping factor ($f_q$), then plasma current ($I_{\text{p}}$) then qstar ($q^*$) then normalized beta ($\beta_{\text{N}}$) then poloidal field and the profile settings for $\mathtt{alphaj}$ ($\alpha_J$) and $\mathtt{rli}$ ($l_{\mathtt{i}}$)

A number of plasma current scaling laws are available in PROCESS. These are calculated in 
routine `calculate_plasma_current()`, in `physics.py`. The safety factor $q_{95}$ required to prevent disruptive MHD instabilities dictates the plasma current $I_{\text{p}}$:

$$\begin{aligned}
I_{\text{p}} = f_q \frac{2\pi}{\mu_0}  \frac{a^2 B_{\text{T}}}{R \ q_{95}}
\end{aligned}$$

$$
\mu_0 I = \mu_0 \int j_{\phi}\cdot \text{d}s = \int \nabla \times B \cdot \text{d}s = \oint B \cdot \text{d}l_{\text{p}} = 2\pi a l B_{\text{p}}(a)
$$

Where $l$ is the ratio of the poloidal plasma circumference to the circumference of the inscribed circle of radius $a$. The function $q$ becomes:

$$
q(r) = \frac{RB_{\phi}}{2\pi}\oint \frac{\text{d}l_{\text{p}}}{R^2B_{\text{p}}} = \frac{rlB_{\phi}}{RB_{\text{p}(r)}}
$$

The factor $f_q$ makes allowance for toroidal effects and plasma shaping (elongation and 
triangularity). Several formulae for this factor are available depending on the value of 
the switch `i_plasma_current`, as follows:

---------------

### 1. Calculate plasma current shaping function $f_q$

------------

#### Peng analytic fit | `calculate_current_coefficient_peng()`

Switch value: `i_plasma_current = 1`

The formula for calculating `fq` is:

$$f_q = \left(\frac{{1.22 - 0.68  \epsilon}}{{(1.0 - \epsilon^2)^2}}\right)  \mathtt{{sf}}^2$$

Where $\epsilon$ is the inverse [aspect ratio](../plasma_geometry.md) ($\mathtt{eps}$) and $\mathtt{sf}$ is the shaping factor calculated in the [poloidal perimeter](../plasma_geometry.md#poloidal-perimeter) function in `plasma_geometry.py`

-----------

#### STAR, Peng double null divertor scaling (ST)

Switch value: `i_plasma_current = 2` [^3] [^4]

This is currently the only scaling in which the calculated plasma current does not follow the form of that show above. The bounds of applicability is for aspect ratios $\le 3$
It uses the `calculate_plasma_current_peng()` function. 

##### `calculate_plasma_current_peng()`

The plasma current is given by:

$$
I_{\text{p}} = \frac{5a B_{\text{T}}\kappa}{2\pi^2 \bar{q}}(F_1+F_2)\left(\frac{\arcsin{E_1}}{E_1}+\frac{\arcsin{E_2}}{E_2}\right)
$$

The values of $F_1$, $F_2$, $d_1$ & $d_2$ are first calculated from the `_plascar_bpol` function.

The values of $E_1$ & $E_2$ are then calculated such as

$$
E_1 = \frac{2\kappa}{d_1(1+\delta)}
$$

$$
E_2 = \frac{2\kappa}{d_2(1.0 - \delta)}
$$    

$I_{\text{p}}$ from above is then calculated from these values

----------------

#### Simple ITER scaling

Switch value: `i_plasma_current = 3` [^5]

The simple cyclindrical case is assumed so:

$$
f_q = 1
$$

-----------------
#### ITER IPDG89 scaling | `calculate_current_coefficient_ipdg89()`

Switch value: `i_plasma_current = 4`[^5] [^7]

The formula for calculating `fq` is:

$$
f_q = \left(\frac{{0.5 \cdot (1.17 - 0.65 \cdot \epsilon)}}{{(1.0 - \epsilon^2)^2}} \cdot \left(1.0 + \kappa_{95}^2 \cdot \left(1.0 + 2.0 \cdot \delta_{95}^2 - 1.2 \cdot \delta_{95}^3\right)\right)\right)
$$





--------------

#### Todd empirical scaling, I | `calculate_current_coefficient_todd()`

Switch value: `i_plasma_current = 5`[^6] [^7]

The formula for calculating `fq` is:

$$
F_{T1}  |  f_q = \left(
                (1.0 + 2.0\epsilon^2)
                \cdot 0.5
                \cdot (1.0 + \kappa_{95}^2)
                \cdot (
                    1.24
                    - 0.54 \cdot \kappa_{95}
                    + 0.3 \cdot (\kappa_{95}^2 + \delta_{95}^2)
                    + 0.125 \cdot \delta_{95}
                )
            \right)
$$


-------------------

#### Todd empirical scaling, II 

Switch value: `i_plasma_current = 6` [^6] [^7]

This function is similar to the previous [Todd scaling](#todd-empirical-scaling-i) except it is mltiplied by a new elongation dependant term


$$
f_q = F_{T1} \times \left(1+[\kappa-1.2]^3\right) 
$$

------------------


#### Connor-Hastie model | `calculate_current_coefficient_hastie()`

Switch value: `i_plasma_current = 7` [^7] [^8]

Asymptotically correct in the range of: 

- $\epsilon \ll 1$
- $\delta \ll 1$ 
- $(\kappa - 1) \ll 1$

Assumes a parabolic profile as seen [here](../profiles/plasma_profiles.md#parabolic-profile-l-mode)

$$
f_q = \left(\frac{(\kappa+1)^2}{2}\right)\left(1+\left(\frac{\kappa+1}{2}\right)^2\epsilon^2+\frac{1}{2}\Delta^{\prime 2}+2\frac{\Delta}{R_0} \\
+ \frac{1}{2}\left(E^{\prime 2}+\frac{E^2}{r^2}\right)+\frac{1}{2}\left(T^{\prime 2}+\frac{4T^2}{r^2}\right)\right)
$$

where:

$$
\frac{T}{r} = \frac{\kappa \delta}{(\kappa+1)^2}  \ \ \ ; \ \ \ \frac{E}{r} = \frac{\kappa-1}{\kappa+1}
$$

$$
T^{\prime} = 2 \left(\frac{1+\lambda}{1+\frac{\lambda}{2}}\right)\frac{T}{r} \ \ \ ; \ \ \ E^{\prime} = 2 \left(\frac{1+\lambda}{1+\frac{\lambda}{3}}\right)\frac{E}{r}
$$

$$
\Delta^{\prime} = \frac{\kappa+1}{2}\epsilon\frac{l_i}{2}+\frac{\beta_0 (1+\lambda)^2}{(\frac{\kappa+1}{2}\epsilon)(1+\nu)}
$$

with 

$$
l_i = \frac{1+\lambda}{\lambda}\left(\frac{1+\lambda}{\lambda}\text{ln}(1-\lambda)-1\right) 
$$

$$
\frac{\Delta}{R_0} = \frac{\beta_0}{6}\left[1+\frac{5}{6}\lambda+\frac{1}{4}\lambda^2\right]+\left(\frac{\kappa+1}{2}\epsilon\right)^2 \frac{1}{8}\left(1-\frac{\lambda^2}{3}\right)
$$

and the parameters $\lambda$ and $\nu$ characterise the current profile $J_{\phi} = \frac{j_0}{(1+\lambda r^2)^2}$ and pressure profile $p = p_0(1-r^2)^{\nu}$



---------------

#### Sauter model, allows negative triangularity

Switch value: `i_plasma_current = 8`[^9]

Assumes zero squareness with the $w_{07}$ parameter ($w_{07}$ = 1). 

The values for $\kappa$ and $\delta$ is taken at the separatrix and is not the 95% values.

!!! note "$w_{07}$ setting"

    There is currently no parameter to set in the input file in order to change $w_{07}$ directly. This can be done directly through the code in `physics.py`

$$
f_q = \frac{4.1 \times 10^6}{\frac{2\pi}{\mu_0}}[1+1.2(\kappa-1)+0.56(\kappa-1)^2] \\
\times (1+0.09\delta +0.16\delta^2)\frac{1+0.45\delta \epsilon}{1-0.74\epsilon}[1+0.55(w_{07}-1)]
$$

----------------------

#### FIESTA for ST's: 

Switch value: `i_plasma_current = 9`[^10]



Assumptions:

 - D-shaped plasmas with $A < 3$:
 - X-pont values of $\kappa$ & $\delta$ used

A set of eqlibria from FIESTA were created and compared to the calculatons from the [Peng double null divertor scaling](plasma_current.md#peng-double-null-divertor-scaling-st) For the low elongation equilibria, the calculated values for
the plasma current were close to those from FIESTA, however moving to
higher elongations causes an underestimate by up to 20%.

Given the parameter dependencies a new plasma current relation based on fits to FIESTA
equilibria was created. It showed no bias with any parameter fitted and that the fit is accurate to 10%.
The linear relation between these and the 95% values expressed in does not hold at high values of elongation and triangularity as per the [`ishape = 0`](../plasma_geometry.md#plasma-geometry-parameters-geomty) relation.


$$
f_q = 0.538 (1.0 + 2.440 \epsilon^{2.736}) \kappa^{2.154}\delta^{0.060}
$$

-----------------------------

### 2. Calculate the cylidrical safety factor

$$
q^* = \frac{5 \times 10^6a^2B_T}{RI_{\text{p}}}\frac{(1+\kappa^2(1+2\delta^2-1.2\delta^3))}{2}
$$

--------------

### 3. Caclulate the normalized beta

The total normlaized beta is calculated as per:

$$
\beta_N = \beta\frac{1\times10^8  a B_{\text{T}}}{I_{\text{P}}}
$$

### 4. Plasma Current Poloidal Field

For calculating the poloidal magnetic field created due to the presence of the plasma current, [Ampere's law](https://en.wikipedia.org/wiki/Amp%C3%A8re%27s_circuital_law) can simply be used. In this case the poloidal field is simply returned as:

$$
B_{\text{p}} = \frac{\mu_0 I_{\text{p}}}{\mathtt{pperim}}
$$

Where `pperim` is the plasma poloidal perimieter calculated [here](../plasma_geometry.md#poloidal-perimeter).

In the case of using the Peng double null scaling ([`i_plasma_current = 2`](plasma_current.md#star-peng-double-null-divertor-scaling-st)), the values $F_1$ and $F_2$ are calculated from [_plasc_bpol](plasma_current.md#_plasc_bpol) and used to calculated the poloidal field from the toroidal field as per:

$$
B_{\text{p}} = B_{\text{T}}\frac{F_1 + F_2}{2\pi \overline{q}}
$$

------------

### 5. Plasma Current Profile Consistency

A limited degree of self-consistency between the plasma current profile and other parameters can be 
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

For spherical tokamaks, the normalised internal inductance can be set from the elongation using `iprofile = 4` or `iprofile = 5` or `iprofile = 6`[^11]:

$$\begin{aligned}
l_i = 3.4 - \kappa_x
\end{aligned}$$

Further desciption of `iprofile` is given in [Beta Limit](../plasma_beta.md).

## _plasc_bpol

This intenral function is used to calculate the plasma shape and poloidal coefficients for calculating the plasma current in the [Peng double null scaling from the STAR code](plasma_current.md#star-peng-double-null-divertor-scaling-st). If this scaling is selected the coefficents are also used to calculate the [poloidal field from the plasma current](plasma_current.md#plasma-current-poloidal-field).


where if $A < \frac{\kappa^2}{(1+\delta)}+\delta$:

$$
F_1 = f_1\left(g-h_1\ln \left(\frac{1+y_1}{1-y_1}\right)\right)
$$

$$
f_1 = \frac{d_1(1+\delta)\epsilon}{(1+\epsilon)(c_1\epsilon -1)}
$$

$$
h_1 = \frac{1+(1-c_1)\frac{\epsilon}{2}}{(1+\epsilon)(c_1 \epsilon -1)}
$$

$$
y_1 = \frac{\sqrt{c_1\epsilon-1}}{1+\epsilon}\frac{1+\delta}{\kappa}
$$

and if $A > \frac{\kappa^2}{(1+\delta)}+\delta$:

$$
F_1 = f_1\left(-g +2h_1 \arctan(y_1)\right)
$$

$$
f_1 = -\frac{d_1(1+\delta)\epsilon}{(1+\epsilon)(c_1\epsilon -1)}
$$

$$
h_1 = \frac{1+(1-c_1)\frac{\epsilon}{2}}{(1+\epsilon)(1-c_1 \epsilon)}
$$

$$
y_1 = \frac{\sqrt{1-c_1\epsilon}}{1+\epsilon}\frac{1+\delta}{\kappa}
$$

where both conditions share:

$$
g = \frac{\epsilon \kappa}{(1-\epsilon \delta)}
$$

$$
F_2 = f_2\left(g +2h_2 \arctan(y_2)\right)
$$

$$
f_2 = -\frac{d_2(1-\delta)\epsilon}{(1-\epsilon)(c_2\epsilon -1)}
$$

$$
h_2 = \frac{1+(c_2-1)\frac{\epsilon}{2}}{(1-\epsilon)(c_2 \epsilon+1)}
$$

$$
E_1 = \frac{2\kappa}{d_1(1+\delta)} ; \ \  E_2 = \frac{2\kappa}{d_2(1-\delta)}
$$

$$
c_1 = \frac{\kappa^2}{(1+\delta)}+\delta \ ; \ \ c_2 = \frac{\kappa^2}{(1-\delta)}+\delta
$$

$$
d_1 = \left(\frac{\kappa}{1+\delta}\right)^2+1 \ ; \ \ d_2 = \left(\frac{\kappa}{1-\delta}\right)^2+1
$$

$$
y_2 = \frac{\sqrt{c_2\epsilon+1}}{1-\epsilon}\frac{1-\delta}{\kappa}
$$

[^3]: Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992). 'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A), 1729–1738. https://doi.org/10.13182/FST92-A29971
[^4]: J.D. Galambos, 'STAR Code : Spherical Tokamak Analysis and Reactor Code',
Unpublished internal Oak Ridge document.
[^5]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
[^6]: D.C.Robinson and T.N.Todd, Plasma and Contr Fusion 28 (1986) 1181
[^7]: T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
[^8]: J.W.Connor and R.J.Hastie, Culham Lab Report CLM-M106 (1985). https://scientific-publications.ukaea.uk/wp-content/uploads/CLM-M106-1.pdf
[^9]: O. Sauter, Geometric formulas for system codes including the effect of negative triangularity, Fusion Engineering and Design, Volume 112, 2016, Pages 633-645, ISSN 0920-3796,
https://doi.org/10.1016/j.fusengdes.2016.04.033.
[^10]: Stuart I. Muldrew, Hanni Lux, Geof Cunningham, Tim C. Hender, Sebastien Kahn, Peter J. Knight, Bhavin Patel, Garry M. Voss, Howard R. Wilson, “PROCESS”: Systems studies of spherical tokamaks, Fusion Engineering and Design, Volume 154, 2020, 111530, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2020.111530.
[^11]: J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,” Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016, doi: https://doi.org/10.1088/0029-5515/56/10/106023.

