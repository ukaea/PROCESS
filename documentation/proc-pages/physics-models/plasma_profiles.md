# Plasma Profiles

???+ Note

    Profile sizes are set to 501 point by default. this can be varied in the `__init__` of `PlasmaProfile` 
All is calculated from `plasma_profiles.py` and `profiles.py`

if physics_variables.tratio > 0.0e0:
            physics_variables.ti = physics_variables.tratio * physics_variables.te    
## Initialization
Plasma profile class is `PlasmaProfile`. Initialization sets profile class size and `neprofile` and `teprofile` to `NProfile` & `TProfile` from `profiles`

| Profile parameter     | Density   |                | Temperature |                |
|-----------------------|-----------|----------------|-------------|----------------|
| Pedestal radius (r/a) | `rhopedn` | $\rho_{ped,n}$ | `rhopedt`   | $\rho_{ped,T}$ |
| Plasma centre value   | `ne0`     | $n_0$          | `te0`       | $T_0$          |
| Pedestal value        | `neped`   | $n_{ped}$      | `teped`     | $T_{ped}$      |
| Separatrix value      | `nesep`   | $n_{sep}$      | `tesep`     | $T_{sep}$      |
| Profile index         | `alphan`  | $\alpha_n$     | `alphat`    | $\alpha_T$     |
| Profile index $\beta$ |           |                | `tbeta`     | $\beta_T$      |


### `NProfile`

def run(self):
        """_summary_
        Subroutine which calls functions and stores nprofile data.
        """
        self.normalise_profile_x()
        self.calculate_profile_dx()
        self.`set_physics_variables`()
        self.calculate_profile_y(
            self.profile_x,
            physics_variables.rhopedn,
            physics_variables.ne0,
            physics_variables.neped,
            physics_variables.nesep,
            physics_variables.alphan,
        )
        self.integrate_profile_y()


Firstly the profile x-dimension is normalised in `normalise_profile_x` by simply dividing the profile size by its max value

The steps between the normalized points is then done by `calculate_profile_dx` which divided the max x-dimension by the number of points.

`set_physics_variables` is then ran which performs `ncore`

 physics_variables.ne0 = self.ncore(
            physics_variables.rhopedn,
            physics_variables.neped,
            physics_variables.nesep,
            physics_variables.dene,
            physics_variables.alphan,
        )
        physics_variables.ni0 = (
            physics_variables.dnitot / physics_variables.dene * physics_variables.ne0
        )

$$
ncore=\frac{1}{3{rhopedn}^2} \left(3nav \left(1+alphan\right)+nsep \left(1+alphan\right) \left(-2+rhopedn+{rhopedn}^2\right)-nped \left(\left(1+alphan\right) \left(1+rhopedn\right)+\left(alphan-2\right) {rhopedn}^2\right)\right)
$$

$$
physics_variables.ni0 = (
            physics_variables.dnitot / physics_variables.dene * physics_variables.ne0
        )
$$

The y profile is then calculated using `calculate_profile_y`

calculate_profile_y(self, rho, rhopedn, n0, nped, nsep, alphan):
        """This routine calculates the density at each normalised minor radius position
        rho for a ELIOS-type density pedestal profile (nprofile).
        Authors:
            R Kemp, CCFE, Culham Science Centre
            H Lux, CCFE, Culham Science Centre
            P J Knight, CCFE, Culham Science Centre
        References:
            J.Johner, Fusion Science and Technology 59 (2011), pp 308-349

        :param rho: normalised minor radius vector
        :type rho: float
        :param rhopedn: normalised minor radius pedestal position
        :type rhopedn: float
        :param n0: central density (/m3)
        :type n0: float
        :param nped: oedestal desnity (/m3)
        :type nped: float
        :param nsep: separatrix density (/m3)
        :type nsep: float
        :param alphan: density peaking parameter
        :type alphan: float
        """

        if physics_variables.ipedestal == 0:
            self.profile_y = n0 * (1 - rho**2) ** alphan

        #  Error trap; shouldn't happen unless volume-averaged density has
        #  been allowed to drop below nped. This may happen during a HYBRD case,
        #  but should have been prevented for optimisation runs.

        #  Input checks

        if n0 < nped:
            logger.info(
                f"NPROFILE: density pedestal is higher than core density. {nped = }, {n0 = }"
            )
        rho_index = rho <= rhopedn
        self.profile_y[rho_index] = (
            nped + (n0 - nped) * (1 - (rho[rho_index] / rhopedn) ** 2) ** alphan
        )
        # Invert the rho_index
        self.profile_y[~rho_index] = nsep + (nped - nsep) * (1 - rho[~rho_index]) / (
            1 - rhopedn
        )

Profile is then integrated with `integrate_profile_y` using simpsons integration

"""
        Integrate profile_y values using scipy.integrate.simpson() function.
        """
        self.profile_integ = integrate.simpson(
            self.profile_y, x=self.profile_x, dx=self.profile_dx
        )

## No pedestal, `(ipedestal = 0)`



if physics_variables.ipedestal == 0:
            self.parabolic_paramterisation()
            self.calculate_profile_factors()
            self.calculate_parabolic_profile_factors()
## Pedestal, `(ipedestal = 1)`
If switch `ipedestal = 0`, no pedestal is present.  The plasma profiles are assumed to be of the form

$$\begin{aligned}
\mbox{Density : } n(\rho) & = n_0 \left( 1 - \rho^2 \right)^{\alpha_n} \\
\mbox{Temperature : } T(\rho) & = T_0 \left( 1 - \rho^2 \right)^{\alpha_T} \\
\mbox{Current : } J(r) & = J_0 \left( 1 - \rho^2 \right)^{\alpha_J}
\end{aligned}$$

where $\rho = r/a$, and $a$ is the plasma minor radius. This gives
volume-averaged values $\langle n \rangle = n_0 / (1+\alpha_n)$, and
line-averaged values $\bar{n} \sim n_0 / \sqrt{(1+\alpha_n)}$, etc.  These
volume- and line-averages are used throughout the code along with the profile
indices $\alpha$, in the various physics models, many of which are fits to
theory-based or empirical scalings. Thus, the plasma model in PROCESS may
be described as 1/2-D.  The relevant profile index variables are
`alphan`, `alphat` and `alphaj`, respectively.

If `ipedestal` = 1 the density and temperature profiles include a pedestal.  
If `ipedestal` = 1 the density and temperature profiles use the forms given below [^1].  

$$\begin{aligned}
\mbox{density:} \qquad n(\rho) = \left\{ 
\begin{aligned}
    & n_{ped} + (n_0 - n_{ped}) \left( 1 -
    \frac{\rho^2}{\rho_{ped,n}^2}\right)^{\alpha_n}
   & \qquad 0 \leq \rho \leq \rho_{ped,n} \\
   & n_{sep} + (n_{ped} - n_{sep})\left( \frac{1- \rho}{1-\rho_{ped,n}}\right)
   & \qquad \rho_{ped,n} < \rho \leq 1
\end{aligned}
\right.
\end{aligned}$$

$$\begin{aligned}
\mbox{temperature:} \qquad T(\rho) = \left\{ 
\begin{aligned}
   & T_{ped} + (T_0 - T_{ped}) \left( 1 - \frac{\rho^{\beta_T}}
    {\rho_{ped,T}^{\beta_T}}\right)^{\alpha_T}  & \qquad 0 \leq \rho \leq \rho_{ped,T} \\
   & T_{sep} + (T_{ped} - T_{sep})\left( \frac{1- \rho}{1-\rho_{ped,T}}\right)
   & \qquad \rho_{ped,T} < \rho \leq 1
\end{aligned}
\right.
\end{aligned}$$

Subscripts $0$, $ped$ and $sep$, denote values at the centre ($\rho = 0$), the
pedestal ($\rho = \rho_{ped}$) and the separatrix ($\rho=1$),
respectively. The density and temperature peaking parameters $\alpha_n$ and a
$\alpha_T$ as well as the second exponent $\beta_T$ (input parameter
`tbeta`, not to be confused with the plasma beta) in the temperature
profile can be chosen by the user, as can the pedestal heights and the values
at the separatrix (`neped, nesep` for the electron density, and
`teped, tesep` for the electron temperature); the ion equivalents are
scaled from the electron values by the ratio of the volume-averaged values).

The density at the centre is given by:

$$\begin{aligned}
  \nonumber
  n_0 & = & \frac{1}{3\rho_{ped,n}^2} \left[3\langle n\rangle (1+\alpha_n)
    + n_{sep} (1+\alpha_n) (-2 + \rho_{ped,n} + \rho_{ped,n}^2) \right.\\
   & & \left. - n_{ped}\left( (1 + \alpha_n)(1+ \rho_{ped,n}) + (\alpha_n -2)
    \rho_{ped,n}^2 \right) \right]
\end{aligned}$$

where $\langle n \rangle$ is the volume-averaged density. The temperature at
the centre is given by

$$\begin{aligned}
T_0 = T_{ped} + \gamma \left[ T_{ped}\, \rho_{ped,T}^2 - \langle T \rangle +
  \frac{1}{3}(1 - \rho_{ped,T}) \left[ \, (1 + 2\rho_{ped,T}) \, T_{ped} + ( 2 +
    \rho_{ped,T}) \, T_{sep} \, \right] \right]
\end{aligned}$$

with 

$$\begin{aligned}
\gamma = \left\{
\begin{aligned}
  & \frac{ -\Gamma(1+\alpha_T+2/\beta_T)}
  {\rho_{ped,T}^2 \, \Gamma(1+\alpha_T) \, \Gamma((2+\beta_T)/\beta_T)}
  \qquad \text{for integer} \, \alpha_T \\
  &\frac{\Gamma(-\alpha_T)\sin(\pi\alpha)\, \Gamma(1+\alpha_T+2/\beta_T)}
  {\pi\rho_{ped,T}^2 \, \Gamma((2+\beta_T)/\beta_T)}
  \qquad \text{for non-integer} \, \alpha_T
\end{aligned}
\right.
\end{aligned}$$

where $\Gamma$ is the gamma function.

Note that density and temperature can have different pedestal positions
$\rho_{ped,n}$ (`rhopedn`) and $\rho_{ped,T}$ (`rhopedt`) in agreement with 
simulations.

If `ipedestal` = 1 or 2 then the pedestal density `neped` is set as a fraction `fgwped` of the 
Greenwald density (providing `fgwped` >= 0).  The default value of `fgwped` is 0.8[^2].

[^1]: M. Bernert et al. Plasma Phys. Control. Fus. **57** (2015) 014038
[^2]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',