# Temperature Profile | `TProfile(Profile)`

The temperature profile class is organised around a central runner function that is called each time the plasma is parameterised by the parent `PlasmaProfile()` class. It is called by `pedestal_parameterisation()` and `parabolic parameterisation()`. The sequence of the runner function can be seen below along with explanation of the following calculations.

## Runner function | `run()`

1. Firstly the profile x-dimension is normalised in [`normalise_profile_x()`](./plasma_profiles_abstract_class.md/#normalise-the-profile-in-x--normalise_profile_x) by simply dividing the profile size by its max value.

2. The steps between the normalized points is then done by [`calculate_profile_dx()`](./plasma_profiles_abstract_class.md#calculate-the-profile-steps-in-x--calculate_profile_dx) which divides the max x-dimension by the number of points.

3. The core electron and ion temperatures are claculated via [`set_physics_variables()`]() 


    ### Calculate core values | `set_physics_variables()`

    The core electron temperature is calculated using the `tcore` method.
    
    #### Electron core density of a pedestalised profile | `tcore()`

    This function calculates the core electron density for a pedestalsied profile in $\text{keV}$. The inclusion of a new $\beta_T$ exponent term allows a more accurate description of temperature profiles with a triangular shape or a strong gradient near the pedestal (characteristic of regimes with an [internal transport barrier](https://wiki.fusion.ciemat.es/wiki/Internal_Transport_Barrier)).

    | Profile parameter / Input               | Temperature   |
    |----------------------------------|-----------|
    | Pedestal radius (r/a)            | `rhopedn`, $\rho_{\text{ped,T}}$ |
    | Pedestal value                   | `neped`, $T_{\text{ped}}$ |
    | Separatrix value                 | `nesep`, $T_{\text{sep}}$ |
    | Average temperature             | `dene`, $\langle T_\text{e} \rangle$ |
    | Profile index/ peaking parameter | `alphan`, $\alpha_T$ |
    | Profile index/ peaking parameter | `tbeta`, $\beta_T$ |


    $$\begin{aligned}
    T_0 = T_{\text{ped}} + \gamma \left[ T_{\text{ped}}\, \rho_{\text{ped},T}^2 - \langle T \rangle +
    \frac{1}{3}(1 - \rho_{\text{ped},T}) \left[ \, (1 + 2\rho_{\text{ped},T}) \, T_{\text{ped}} + ( 2 +
        \rho_{\text{ped},T}) \, T_{\text{sep}} \, \right] \right]
    \end{aligned}$$

    with

    $$\begin{aligned}
    \gamma = \left\{
    \begin{aligned}
    & \frac{ -\Gamma(1+\alpha_T+2/\beta_T)}
    {\rho_{\text{ped},T}^2 \, \Gamma(1+\alpha_T) \, \Gamma((2+\beta_T)/\beta_T)}
    \qquad \text{for integer} \, \alpha_T \\
    &\frac{\Gamma(-\alpha_T)\sin(\pi\alpha)\, \Gamma(1+\alpha_T+2/\beta_T)}
    {\pi\rho_{\text{ped},T}^2 \, \Gamma((2+\beta_T)/\beta_T)}
    \qquad \text{for non-integer} \, \alpha_T
    \end{aligned}
    \right.
    \end{aligned}$$


    where $\Gamma$ is the gamma function.

    ##### Derivation

    Taking the standard form of the profile both in the core region and the pedestal region and multipling each function by $2\pi\rho$ and integrating within their bounds of appliability to get their volume of integration. 
    
    $$\begin{aligned}
	\qquad T(\rho) =  
    & T_{ped} + (T_0 - T_{ped}) \left( 1 - \frac{\rho^{\beta_T}}
    {\rho_{ped,T}^{\beta_T}}\right)^{\alpha_T}  & \qquad 0 \leq \rho \leq \rho_{ped,T} \\
    & T_{sep} + (T_{ped} - T_{sep})\left( \frac{1- \rho}{1-\rho_{ped,T}}\right)
    & \qquad \rho_{ped,T} < \rho \leq 1
	\end{aligned}$$

    $$
    \begin{split}
        & \int_{0}^{\rho_{ped,T}} 2\pi \rho \left( T_{ped} + (T_0 - T_{ped}) \left( 1 - \frac{\rho^{\beta_T}}{\rho_{ped,T}^{\beta_T}}\right)^{\alpha_T} \right) \, d\rho  \\ 
        & = \pi  \left(\rho_{ped,T}^2 T_{ped} + \frac{\pi  (-1)^{2/\beta } \csc (\pi  \alpha ) \left(-\rho_{ped,T}^{-\beta }\right)^{-2/\beta } \Gamma \left(\frac{\beta +2}{\beta }\right) (T_{ped}-T_0)}{\Gamma (-\alpha ) \Gamma \left(\alpha +\frac{2}{\beta }+1\right)}+\right)
    \end{split}
    $$

    $$
    \begin{equation}\label{eq:3}
	\begin{split}
		& \int_{\rho_{ped,T}}^{1} 2\pi \rho \left( T_{sep} + (T_{ped} - T_{sep})\left( \frac{1- \rho}{1-\rho_{ped,T}}\right) \right) \, d\rho  \\ 
		& = -\frac{1}{3} \pi  (\rho_{ped,T}-1) (2 \rho_{ped,T} T_{ped}+T_{ped}+(\rho_{ped,T}+2) T_{sep})
	\end{split}
    \end{equation}
    $$

    Summing these two together equates to the volume averaged temperature, $\langle T \rangle$. This can then be substituted back in 

    (2) + (3) $= \langle T \rangle$
    
    
   
    $$
    \begin{aligned}
    \nonumber
    T_0 = T_{\text{ped}} + \frac{1}{3 \pi  \Gamma \left(\frac{\beta +2}{\beta}\right)}(-1)^{-2/\beta}\sin (\pi  \alpha) \left(-\rho_{ped,T}^{-\beta}\right)\left(\langle T \rangle-\rho_{ped,T}^2 T_{ped}\right)^{2/\beta}
    & \quad (\rho_{ped,T}-1)(2 \rho_{ped,T} T_{ped}+T_{ped}+(\rho_{ped,T}+2) T_{sep}) -3 \Gamma (-\alpha ) \Gamma \left(\alpha +\frac{2}{\beta }+1\right)
    \end{aligned} 
    $$
    

The core ion temperature is then set such as:

$$
T_{i,0} = \left(\frac{T_i}{T_e}\right)T_{e,0}
$$


4. The y profile is then calculated using `calculate_profile_y()`. This routine calculates the temperature at each normalised minor radius position $\rho$ for a HELIOS-type temperature pedestal profile (nprofile)[^1]

    ### Calculate temperature at each radius position | `calculate_profile_y()`

    Normalized 

    | Profile parameter / Input               | Density   |
    |----------------------------------|-----------|
    | Normalized plasma radii            | `profile_x` |
    | Pedestal radius (r/a)            | `rhopedt`, $\rho_{\text{ped,T}}$ |
    | Core density                | `te0`, $T_{\text{e0}}$ |
    | Pedestal value                   | `teped`, $T_{\text{ped}}$ |
    | Separatrix value                 | `tesep`, $T_{\text{sep}}$ |
    | Profile index/ peaking parameter | `alphat`, $\alpha_T$ |
    | 2nd profile index/ peaking parameter | `tbeta`, $\beta_T$ |

    If `ipedestal == 0` then the original parabolic profile form is used

    $$
    T(\rho) = T_0(1 - \rho^2)^{\alpha_T} 
    $$

    The central tmeprature ($T_0$) is then checked to make sure it is not less than the pedestal temperature, $T_{\text{ped}}$.
    If it is less then a logger warning is pushed to the terminal at runtime.

    Values of the profile temperature are then assigned based on the desnity function below across bounds from 0 to `rhopedn` and `rhopedn` to 1.  



    $$\begin{aligned}
    \mbox{Temperature:} \ \ T(\rho) = \left\{ 
    \begin{aligned}
    & \text{T}_{\text{ped}} + (T_0 - \text{T}_{\text{ped}}) \left( 1 - \frac{\rho^{\beta_T}}
    {\rho_{\text{ped},T}^{\beta_T}}\right)^{\alpha_T}  &  0 \leq \rho \leq \rho_{\text{ped},T} \\
    & \text{T}_{\text{sep}} + (\text{T}_{\text{ped}} - \text{T}_{\text{sep}})\left( \frac{1- \rho}{1-\rho_{\text{ped},T}}\right)
    &  \rho_{\text{ped},T} < \rho \leq 1
    \end{aligned}
    \right.
    \end{aligned}$$ 
        

5. Profile is then integrated with `integrate_profile_y()` using Simpsons integration from the profile abstract base class

[^1]: Jean, J. (2011). *HELIOS: A Zero-Dimensional Tool for Next Step and Reactor Studies*. Fusion Science and Technology, 59(2), 308–349. https://doi.org/10.13182/FST11-A11650