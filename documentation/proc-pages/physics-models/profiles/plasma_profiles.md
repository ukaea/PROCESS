# Plasma Profiles | `PlasmaProfile`

!!! warning " Un-realistic profiles"

    If `ipedestal >= 1` it is highly recommended to use constraint equation 81 (icc=81). This enforces solutions in which $n_0$ has to be greater than $n_{ped}$. 
    Negative $n_0$ values can also arise during iteration, so it is important to be weary on how low the lower bound for $n_e (\mathtt{dene})$ is set.

In `PROCESS` the density, temperature and current profiles of the plasma for electrons and ions can take two forms depending on the switch value for `ipedestal`. Either without a [pedestal](http://fusionwiki.ciemat.es/wiki/Pedestal), `ipedestal = 0` or with a pedestal `ipedestal = 1`. ADD MORE DESCIPTION IN HERE

The files responsible for calculting and storing the profiles are `plasma_profiles.py` and `profiles.py`. A central plasma profile object is created from the `PlasmaProfile` class that contains attributes for the plasma density and temperature. The density and temperature profiles are in themselves objects of the [`Profile`](./plasma_profiles_abstract_class.md) abstract base class. [`Profile`](./plasma_profiles_abstract_class.md), `NProfile` and `TProfile` are all defined in `profiles.py`. `PlasmaProfile` is exclusively in `plasma_profiles.py` 

<figure markdown>
![UML of profiles](../../images/profiles_uml.png){ width="100%"}
<figcaption>Figure 1: UML class breakdown of the plasma profiles</figcaption>
</figure>





------------

If `ipedestal = 0` then no pedestal is present and the function describing the profiles is given by:

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

| Profile parameter                | Density   | Temperature | Current  |
|----------------------------------|-----------|-------------|----------------|
| Plasma centre value              | `ne0`, $n_0$         | `te0`, $T_0$        |  `plascur`, $J_0$        |
| Profile index/ peaking parameter | `alphan`, $\alpha_n$ | `alphat`, $\alpha_T$    |  `alphaj`, $\alpha_J$    |



If `ipedestal=1` there is now a pedestal present in the profile and they follow the form shown below:


$$\begin{aligned}
\mbox{Density:} \qquad n(\rho) = \left\{ 
\begin{aligned}
    & \text{n}_{\text{ped}} + (n_0 - \text{n}_{\text{ped}}) \left( 1 -
    \frac{\rho^2}{\rho_{\text{ped},n}^2}\right)^{\alpha_n}
   & \qquad 0 \leq \rho \leq \rho_{\text{ped},n} \\
   & \text{n}_{\text{sep}} + (\text{n}_{\text{ped}} - \text{n}_{\text{sep}})\left( \frac{1- \rho}{1-\rho_{\text{ped},n}}\right)
   & \qquad \rho_{\text{ped},n} < \rho \leq 1
\end{aligned}
\right.
\end{aligned}$$


$$\begin{aligned}
\mbox{Temperature:} \qquad T(\rho) = \left\{ 
\begin{aligned}
   & \text{T}_{\text{ped}} + (T_0 - \text{T}_{\text{ped}}) \left( 1 - \frac{\rho^{\beta_T}}
    {\rho_{\text{ped},T}^{\beta_T}}\right)^{\alpha_T}  & \qquad 0 \leq \rho \leq \rho_{\text{ped},T} \\
   & \text{T}_{\text{sep}} + (\text{T}_{\text{ped}} - \text{T}_{\text{sep}})\left( \frac{1- \rho}{1-\rho_{\text{ped},T}}\right)
   & \qquad \rho_{\text{ped},T} < \rho \leq 1
\end{aligned}
\right.
\end{aligned}$$ 

| Profile parameter                | Density   |        Temperature |                
|----------------------------------|-----------|-------------|
| Pedestal radius (r/a)            | `rhopedn`,$\rho_{ped,n}$ |   `rhopedt`, $\rho_{ped,T}$   |  
| Plasma centre value              | `ne0`, $n_0$      |           `te0`, $T_0$       |           
| Pedestal value                   | `neped`, $n_{ped}$    |       `teped`, $T_{ped}$     |       
| Separatrix value                 | `nesep`, $n_{sep}$   |        `tesep`, $T_{sep}$     |       
| Profile index/ peaking parameter | `alphan`, $\alpha_n$  |       `alphat`, $\alpha_T$    |      
| Profile index $\beta$            |           |                 `tbeta`, $\beta_T$     |       






## Initialization
The parent plasma profile class is `PlasmaProfile`. Initialization sets the profile class size and `neprofile` and `teprofile` to `NProfile` & `TProfile` from `profiles`

???+ Note

    Profile sizes are set to 501 point by default. this can be varied in the `__init__` of `PlasmaProfile`. Changing this will affect the values when doing Simpsons rule integration on the profiles. 



| Profile parameter                | Density   |                | Temperature |                |
|----------------------------------|-----------|----------------|-------------|----------------|
| Pedestal radius (r/a)            | `rhopedn` | $\rho_{ped,n}$ | `rhopedt`   | $\rho_{ped,T}$ |
| Plasma centre value              | `ne0`     | $n_0$          | `te0`       | $T_0$          |
| Pedestal value                   | `neped`   | $n_{ped}$      | `teped`     | $T_{ped}$      |
| Separatrix value                 | `nesep`   | $n_{sep}$      | `tesep`     | $T_{sep}$      |
| Profile index/ peaking parameter | `alphan`  | $\alpha_n$     | `alphat`    | $\alpha_T$     |
| Profile index $\beta$            |           |                | `tbeta`     | $\beta_T$      |




### Temperature `TProfle()`

1. Firstly the profile x-dimension is normalised in `normalise_profile_x()` by simply dividing the profile size by its max value

2. The steps between the normalized points is then done by `calculate_profile_dx()` which divided the max x-dimension by the number of points.


3. `set_physics_variables()` is then ran which performs `tcore()` which calculates the central electron density. The ion central density is then calculated by the ratio from this scaling.

-------------------------------------

#### `tcore`

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


-------------------------------------------------


$$
T_{i,0} = \left(
           \frac{T_i}{T_e}T_{e,0} 
        \right)
$$

4. The y profile is then calculated using `calculate_profile_y()`. This routine calculates the temperature at each normalised minor radius position $\rho$ for a HELIOS-type density pedestal profile (tprofile)[^3]

If `ipedestal == 0` then the original parabolic profile form is used.

$$
T(\rho) = T_0 \left( 1 - \rho^2 \right)^{\alpha_T}  
$$

The central temperature ($T_0$) is then checked to make sure it is not less than the pedestal temperature, $n_{ped}$    


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

5. Profile is then integrated with `integrate_profile_y()` using simpsons integration




## Plasma Parameterization
Ion temperature is set with `tratio` which just takes $T_i = \mathtt{tratio}\times T_e$

if physics_variables.tratio > 0.0e0:
            physics_variables.ti = physics_variables.tratio * physics_variables.te 

### ipedestal = 0

#### `parabolic_paramterisation()`
`parabolic_paramterisation()` is ran
 Pedestal values resent to agree with original parabolinc profiles.

 `Nprofile()` and `TProfile()` is re-ran


profile factor is calculated:

$$
\mathtt{pcoef} = \frac{(1+\alpha_n)(1+\alpha_T)}{(1+\alpha_T + \alpha_n)}
$$

line avergaed electron density is calculated

$$
\mathtt{dnla} = n_e \times (1+\alpha_n) \times 0.886227 \times \mathtt{gamfun}(\alpha_n+1) / \mathtt{gamfun}(\alpha_n+1.5)
$$

where $\mathtt{gamfun}$ is a gamma function calculator from `maths_library.f90` and is found below.

```fortran
recursive function gamfun(x) result(gamma)
    !! Calculates the gamma function for arbitrary real x
    !! author: P J Knight, CCFE, Culham Science Centre
    !! x : input real : gamma function argument
    !! This routine evaluates the gamma function, using an
    !! asymptotic expansion based on Stirling's approximation.
    !! http://en.wikipedia.org/wiki/Gamma_function
    !! T&amp;M/PKNIGHT/LOGBOOK24, p.5
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    !  Arguments
    real(dp), intent(in) :: x
    real(dp) :: gamma
    !  Local variables
    real(dp), parameter :: sqtwopi = 2.5066282746310005D0
    real(dp), parameter :: c1 = 8.3333333333333333D-2  !  1/12
    real(dp), parameter :: c2 = 3.4722222222222222D-3  !  1/288
    real(dp), parameter :: c3 = 2.6813271604938272D-3  !  139/51840
    real(dp), parameter :: c4 = 2.2947209362139918D-4  !  571/2488320
    real(dp) :: summ, denom
    integer :: i,n
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (x > 1.0D0) then
       summ = 1.0D0 + c1/x + c2/x**2 - c3/x**3 - c4/x**4
       gamma = exp(-x) * x**(x-0.5D0) * sqtwopi * summ
    else
       !  Use recurrence formula to shift the argument to >1
       !  gamma(x) = gamma(x+n) / (x*(x+1)*(x+2)*...*(x+n-1))
       !  where n is chosen to make x+n > 1
       n = int(-x) + 2
       denom = x
       do i = 1,n-1
          denom = denom*(x+i)
       end do
       gamma = gamfun(x+n)/denom
    end if
  end function gamfun
```
Set the density weighted temperatures

$$\begin{aligned}
\mathtt{ten} = \mathtt{pcoef}T_e \\
\mathtt{tin} = \mathtt{pcoef}T_i
\end{aligned}$$

Calculate central values for temperature and density

$$\begin{aligned}
\mathtt{te0} = T_e \times (1+\alpha_T) \\
\mathtt{ti0} = T_i \times (1+\alpha_T)
\end{aligned}$$

$$\begin{aligned}
\mathtt{ne0} = n_e \times (1+\alpha_n) \\
\mathtt{ni0} = \mathtt{dnitot} \times (1+\alpha_n)
\end{aligned}$$


#### `calculate_profile_factors()`

The central plasma pressure is calculated from the ideal gas law.

$$
p_0 = (\mathtt{ne0} \times \mathtt{te0}+\mathtt{ni0}\times \mathtt{ti0})\times 1000 \times e
$$

Pressure profile index (N.B. no pedestal effects included here)
N.B. p0 is NOT equal to <p> * (1 + alphap), but p(rho) = n(rho)*T(rho)
and <p> = <n>.T_n where <...> denotes volume-averages and T_n is the
density-weighted temperature

$$
\alpha_p = \alpha_n + \alpha_T
$$

#### `calculate_parabolic_profile_factor()`

$$dtdrho_{max}=\left(-2^{alphat} \left(-1+alphat\right)^{-1+alphat} alphat\times-1+2alphat\right)^{0.5\times{10}^0-alphat} te0$$

```python
def calculate_parabolic_profile_factors():
        """The gradient information for ipedestal = 0:
        All formulas can be obtained from the analytical parametric form of the ipedestal profiles
        rho_max is obtained by equalling the second derivative to zero e.g.
        """
        if physics_variables.ipedestal == 0:
            if physics_variables.alphat > 1.0:
                # Rho (normalized radius), where temperature derivative is largest
                rho_te_max = 1.0 / np.sqrt(-1.0 + 2.0 * physics_variables.alphat)
                dtdrho_max = (
                    -(2.0**physics_variables.alphat)
                    * (-1.0 + physics_variables.alphat)
                    ** (-1.0 + physics_variables.alphat)
                    * physics_variables.alphat
                    * (-1.0 + 2.0 * physics_variables.alphat)
                    ** (0.5e0 - physics_variables.alphat)
                    * physics_variables.te0
                )
                te_max = (
                    physics_variables.te0
                    * (1 - rho_te_max**2) ** physics_variables.alphat
                )

            elif physics_variables.alphat <= 1.0 and physics_variables.alphat > 0.0:
                # This makes the profiles very 'boxy'
                # The gradient diverges here at the edge so define some 'wrong' value of 0.9
                # to approximate the gradient
                rho_te_max = 0.9
                dtdrho_max = (
                    -2.0
                    * physics_variables.alphat
                    * rho_te_max
                    * (1 - rho_te_max**2) ** (-1.0 + physics_variables.alphat)
                    * physics_variables.te0
                )
                te_max = (
                    physics_variables.te0
                    * (1 - rho_te_max**2) ** physics_variables.alphat
                )
            else:
                raise ValueError(f"alphat is negative: { physics_variables.alphat}")

            # Same for density
            if physics_variables.alphan > 1.0:
                rho_ne_max = 1.0 / np.sqrt(-1.0 + 2.0 * physics_variables.alphan)
                dndrho_max = (
                    -(2.0**physics_variables.alphan)
                    * (-1.0 + physics_variables.alphan)
                    ** (-1.0 + physics_variables.alphan)
                    * physics_variables.alphan
                    * (-1.0 + 2.0 * physics_variables.alphan)
                    ** (0.5 - physics_variables.alphan)
                    * physics_variables.ne0
                )
                ne_max = (
                    physics_variables.ne0
                    * (1e0 - rho_ne_max**2) ** physics_variables.alphan
                )
            elif physics_variables.alphan <= 1.0 and physics_variables.alphan > 0.0:
                # This makes the profiles very 'boxy'
                # The gradient diverges here at the edge so define some 'wrong' value of 0.9
                # to approximate the gradient
                rho_ne_max = 0.9
                dndrho_max = (
                    -2.0
                    * physics_variables.alphan
                    * rho_ne_max
                    * (1 - rho_ne_max**2) ** (-1.0 + physics_variables.alphan)
                    * physics_variables.ne0
                )
                ne_max = (
                    physics_variables.ne0
                    * (1 - rho_ne_max**2) ** physics_variables.alphan
                )
            else:
                raise ValueError(f"alphan is negative: { physics_variables.alphan}")

            # set normalized gradient length
            # te at rho_te_max
            physics_variables.gradient_length_te = (
                -dtdrho_max * physics_variables.rminor * rho_te_max / te_max
            )
            # same for density:
            physics_variables.gradient_length_ne = (
                -dndrho_max * physics_variables.rminor * rho_ne_max / ne_max
            )
```

### ipedestal = 1

#### `pedestal_parameterisation()`

 `Nprofile()` and `TProfile()` is re-ran

 Perform integrations to calculate ratio of density-weighted to volume-averaged temperature, etc. Density-weighted temperature = $\frac{\int{nT \ dV}}{\int{n \ dV}}$,  which is approximately equal to the ratio $\frac{\int{\rho \ n(\rho) T(\rho) \ d\rho}}{\int{\rho \ n(\rho) \ d\rho}}$

 Density weighted temperatures are thus set as such

$$
\mathtt{ten} = \frac{\int{\rho \ n(\rho) T(\rho) \  d\rho}}{\int{\rho \ n(\rho) \ d\rho}} \\
\mathtt{tin} = \mathtt{ten}\frac{T_i}{T_e}
$$

Set profile factor:

$$
\mathtt{pcoef} = \frac{\mathtt{ten}}{T_e}
$$

Caclulate the line avaerged electron density:

$$
\mathtt{dnla} = \int{n(\rho) \ d\rho}
$$

#### `calculate_profile_factors()`

The central plasma pressure is calculated from the ideal gas law.

$$
p_0 = (\mathtt{ne0} \times \mathtt{te0}+\mathtt{ni0}\times \mathtt{ti0})\times 1000 \times e
$$

Pressure profile index (N.B. no pedestal effects included here)
N.B. p0 is NOT equal to <p> * (1 + alphap), but p(rho) = n(rho)*T(rho)
and <p> = <n>.T_n where <...> denotes volume-averages and T_n is the
density-weighted temperature

$$
\alpha_p = \alpha_n + \alpha_T
$$

If `ipedestal` = 1 or 2 then the pedestal density `neped` is set as a fraction `fgwped` of the 
Greenwald density (providing `fgwped` >= 0).  The default value of `fgwped` is 0.8[^2].

[^1]: M. Bernert et al. Plasma Phys. Control. Fus. **57** (2015) 014038
[^2]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
[^3]: Johner Jean (2011) HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies, Fusion Science and Technology, 59:2, 308-349, DOI: 10.13182/FST11-A11650


# Plasma Profiles

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

If `ipedestal` = 1, 2 or 3 the density and temperature profiles include a pedestal.  
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

!!! warning " Un-realistic profiles"

    If `ipedestal >= 1` it is highly recommended to use constraint equation 81 (icc=81). This enforces solutions in which $n_0$ has to be greater than $n_{ped}$. 
    Negative $n_0$ values can also arise during iteration, so it is important to be weary on how low the lower bound for $n_e (\mathtt{dene})$ is set.

[^1]: M. Bernert et al. Plasma Phys. Control. Fus. **57** (2015) 014038    

[^2]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',