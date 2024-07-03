# Temperature Profile | `TProfile(Profile)`

The temperature profile class is organised around a central runner function that is called each time the plasma is parameterised by the parent `PlasmaProfile()` class. It is called by `pedestal_parameterisation()` and `parabolic parameterisation()`. The sequence of the runner function can be seen below along with explanation of the following calculations.

## Runner function | `run()`

1. Firstly the profile x-dimension is normalised in [`normalise_profile_x()`](./plasma_profiles_abstract_class.md/#normalise-the-profile-in-x--normalise_profile_x) by simply dividing the profile size by its max value.

2. The steps between the normalized points is then done by [`calculate_profile_dx()`](./plasma_profiles_abstract_class.md#calculate-the-profile-steps-in-x--calculate_profile_dx) which divides the max x-dimension by the number of points.

3. The core electron and ion temperatures are claculated via [`set_physics_variables()`]() 


    ### Calculate core values | `set_physics_variables()`

    The core electron temperature is calculated using the `tcore` method.
    
    #### Electron core density of a pedestalised profile | `tcore()`

    This function calculates the core electron density for a pedestalsied profile in $\text{keV}$. The inclusion of a new $\beta_T$ exponent term allows a more accurate description of temperature profiles with a triangular shape or a strong gradient near the pedestal (characteristic of regimes with an internal transport barrier).

    | Profile parameter / Input               | Temperature   |
    |----------------------------------|-----------|
    | Pedestal radius (r/a)            | `rhopedn`, $\rho_{\text{ped,T}}$ |
    | Pedestal value                   | `neped`, $T_{\text{ped}}$ |
    | Separatrix value                 | `nesep`, $T_{\text{sep}}$ |
    | Average temperature             | `dene`, $\langle T_e \rangle$ |
    | Profile index/ peaking parameter | `alphan`, $\alpha_T$ |
    | Profile index/ peaking parameter | `tbeta`, $\beta_T$ |


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

    If `ncore` is returned as being less that 0, it is forced into a state of `ncore = 1E-6` in order to help convergance. This will also give a warning to the user to raise the lower bound of the average electron density `dene`.

    The core ion density is then set from the output out `dnitot` which is the total ion density such as:

    $$
    n_{i,0} = \left(\frac{n_i}{n_e}\right)n_{e,0}
    $$


4. The y profile is then calculated using `calculate_profile_y()`. This routine calculates the density at each normalised minor radius position $\rho$ for a HELIOS-type density pedestal profile (nprofile)[^1]

    ### Calculate desnity at each radius position | `calculate_profile_y()`

    Normalized 

    | Profile parameter / Input               | Density   |
    |----------------------------------|-----------|
    | Normalized plasma radii            | `profile_x` |
    | Pedestal radius (r/a)            | `rhopedn`, $\rho_{\text{ped,n}}$ |
    | Core density                | `ne0`, $n_{\text{e0}}$ |
    | Pedestal value                   | `neped`, $n_{\text{ped}}$ |
    | Separatrix value                 | `nesep`, $n_{\text{sep}}$ |
    | Profile index/ peaking parameter | `alphan`, $\alpha_n$ |

    If `ipedestal == 0` then the original parabolic profile form is used

    $$
    n(\rho) = n_0(1 - \rho^2)^{\alpha_n} 
    $$

    The central density ($n_0$) is then checked to make sure it is not less than the pedestal density, $n_{\text{ped}}$.
    If it is less then a logger warning is pushed to the terminal at runtime.

    Values of the profile density are then assigned based on the desnity function below across bounds from 0 to `rhopedn` and `rhopedn` to 1.  



    $$\begin{aligned}
    \mbox{Density:} \ n(\rho) = \left\{ 
    \begin{aligned}
        & n_{\text{ped}} + (n_0 - n_{\text{ped}}) \left( 1 -
        \frac{\rho^2}{\rho_{\text{ped,n}}^2}\right)^{\alpha_n}
    & \ 0 \leq \rho \leq \rho_{\text{ped,n}} \\
    & n_{\text{sep}} + (n_{\text{ped}} - n_{\text{sep}})\left( \frac{1- \rho}{1-\rho_{\text{ped,n}}}\right)
    & \ \rho_{\text{ped,n}} < \rho \leq 1
    \end{aligned}
    \right.
    \end{aligned}$$
        

5. Profile is then integrated with `integrate_profile_y()` using Simpsons integration from the profile abstract base class

[^1]: Jean, J. (2011). *HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies*. Fusion Science and Technology, 59(2), 308–349. https://doi.org/10.13182/FST11-A11650