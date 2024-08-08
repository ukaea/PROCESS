# Plasma Current Scaling Laws

A number of plasma current scaling laws are available in PROCESS [^1]. These are calculated in 
routine `culcur`, which is called by `physics`. The safety factor $q_{95}$ required to prevent 
disruptive MHD instabilities dictates the plasma current Ip:

$$\begin{aligned}
I_p = \frac{2\pi}{\mu_0} B_t \frac{a^2 f_q}{Rq_{95}}
\end{aligned}$$

The factor $f_q$ makes allowance for toroidal effects and plasma shaping (elongation and 
triangularity). Several formulae for this factor are available [2,3] depending on the value of 
the switch `icurr`, as follows:

| `icurr` | Description |
| :-: | - |
| 1 | Peng analytic fit | 
| 2 | Peng double null divertor scaling (ST)[^4] | 
| 3 | Simple ITER scaling | 
| 4 | Revised ITER scaling[^5]  $f_q = \frac{1.17-0.65\epsilon}{2(1-\epsilon^2)^2} (1 + \kappa_{95}^2 (1+2\delta_{95}^2 - 1.2\delta_{95}^3) )$| 
| 5 | Todd empirical scaling, I | 
| 6 | Todd empirical scaling, II | 
| 7 | Connor-Hastie model | 
| 8 | Sauter model, allows negative $\delta$ | 
| 9 | Scaling for spherical tokamaks, based on a fit to a family of equilibria derived by Fiesta: $f_q = 0.538 (1 + 2.44\epsilon^{2.736}) \kappa^{2.154} \delta^{0.06}$|

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

It is recommended that current scaling law `icurr = 4` is used if `iprofile = 1`. 
This relation is only applicable to large aspect ratio tokamaks.

For spherical tokamaks, the normalised internal inductance can be set from the elongation using `iprofile = 4` or `iprofile = 5`:

$$\begin{aligned}
l_i = 3.4 - \kappa_x
\end{aligned}$$

Further desciption of `iprofile` is given in [Beta Limit](./plasma_beta.md).

## Bootstrap, Diamagnetic and Pfirsch-Schlüter Current Scalings

The fraction of the plasma current provided by the bootstrap effect
can be either input into the code directly, or calculated using one of four
methods, as summarised here. Note that methods `ibss = 1-3` do not take into account the 
existence of pedestals, whereas the Sauter et al. scaling 
(`ibss = 4`) allows general profiles to be used. 

| `ibss` | Description |
| :-: | - |
| 1 | ITER scaling -- To use the ITER scaling method for the bootstrap current fraction.  Set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$). This method is valid at high aspect ratio only.
| 2 | General scaling -- To use a more general scaling method, set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$).
| 3 | Numerically fitted scaling [^8] -- To use a numerically fitted scaling method, valid for all aspect ratios, set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$).
| 4 | Sauter, Angioni and Lin-Liu scaling [^9] [^10] -- Set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$).
| 5 | Sakai, Fujita and Okamoto scaling [^11] -- Set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$). The model includes the toroidal diamagnetic current in the calculation due to the dataset, so `idia = 0` can only be used with it

!!! Note "Fixed Bootstrap Current"
    Direct input -- To input the bootstrap current fraction directly, set `bscfmax` 
    to $(-1)$ times the required value (e.g. -0.73 sets the bootstrap faction to 0.73).

The diamagnetic current fraction $f_{dia}$ is strongly related to $\beta$ and is typically small,
hence it is usually neglected.  For high $\beta$ plasmas, such as those at tight
aspect ratio, it should be included and two scalings are offered.  If the diamagnetic
current is expected to be above one per cent of the plasma current, a warning
is issued to calculate it.

`idia = 0` Diamagnetic current fraction is zero.

`idia = 1` Diamagnetic current fraction is calculated using a fit to spherical tokamak calculations by Tim Hender:

$$f_{dia} = \frac{\beta}{2.8}$$

`idia = 2` Diamagnetic current fraction is calculated using a SCENE fit for all aspect ratios:

$$f_{dia} = 0.414 \space \beta \space (\frac{0.1 q_{95}}{q_0} + 0.44)$$

A similar scaling is available for the Pfirsch-Schlüter current fraction $f_{PS}$.  This is
typically smaller than the diamagnetic current, but is negative.

`ips = 0` Pfirsch-Schlüter current fraction is set to zero.

`ips = 1` Pfirsch-Schlüter current fraction is calculated using a SCENE fit for all aspect ratios:

$$ f_{PS} = -0.09 \beta $$

There is no ability to input the diamagnetic and Pfirsch-Schlüter current
directly.  In this case, it is recommended to turn off these two scalings 
and to use the method of fixing the bootstrap current fraction.

[^1]: D.J. Ward, 'PROCESS Fast Alpha Pressure', Work File Note F/PL/PJK/PROCESS/CODE/050
[^2]:  Albajar, Nuclear Fusion **41** (2001) 665
[^3]: M. Kovari, R. Kemp, H. Lux, P. Knight, J. Morris, D.J. Ward, '“PROCESS”: A systems code for fusion power plants—Part 1: Physics' Fusion Engineering and Design 89 (2014) 3054–3069
[^4]: J.D. Galambos, 'STAR Code : Spherical Tokamak Analysis and Reactor Code',
Unpublished internal Oak Ridge document.
[^5]: W.M. Nevins, 'Summary Report: ITER Specialists' Meeting on Heating and
Current Drive', ITER-TN-PH-8-4, 13--17 June 1988, Garching, FRG
[^6]: Y. Sakamoto, 'Recent progress in vertical stability analysis in JA',
Task meeting EU-JA #16, Fusion for Energy, Garching, 24--25 June 2014
[^7]: Menard et al. (2016), Nuclear Fusion, 56, 106023
[^8]: H.R. Wilson, Nuclear Fusion **32** (1992) 257
[^9]: O. Sauter, C. Angioni and Y.R. Lin-Liu, Physics of Plasmas **6** (1999) 2834 
[^10]: O. Sauter, C. Angioni and Y.R. Lin-Liu, Physics of Plasmas **9** (2002) 5140
[^11]: Ryosuke Sakai, Takaaki Fujita, Atsushi Okamoto, Derivation of bootstrap current fraction scaling formula for 0-D system code analysis, Fusion Engineering and Design, Volume 149, 2019, 111322, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2019.111322.
