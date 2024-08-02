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

[^1]: D.J. Ward, 'PROCESS Fast Alpha Pressure', Work File Note F/PL/PJK/PROCESS/CODE/050
[^2]:  Albajar, Nuclear Fusion **41** (2001) 665
[^3]: M. Kovari, R. Kemp, H. Lux, P. Knight, J. Morris, D.J. Ward, '“PROCESS”: A systems code for fusion power plants—Part 1: Physics' Fusion Engineering and Design 89 (2014) 3054–3069
[^4]: J.D. Galambos, 'STAR Code : Spherical Tokamak Analysis and Reactor Code',
Unpublished internal Oak Ridge document.
[^5]: W.M. Nevins, 'Summary Report: ITER Specialists' Meeting on Heating and
Current Drive', ITER-TN-PH-8-4, 13--17 June 1988, Garching, FRG
[^6]: Y. Sakamoto, 'Recent progress in vertical stability analysis in JA',
Task meeting EU-JA #16, Fusion for Energy, Garching, 24--25 June 2014