# Density Profile | `NProfile(Profile)`

The desnity profile class is organised around a central runner function that is called each time the plasma is parameterised by the parent [`PlasmaProfile()`](./plasma_profiles.md) class. It is called by `pedestal_parameterisation()` and `parabolic parameterisation()`. The sequence of the runner function can be seen below along with explanation of the following calculations.

## Runner function | `run()`

1. Firstly the profile x-dimension is normalised in [`normalise_profile_x()`](./plasma_profiles_abstract_class.md/#normalise-the-profile-in-x--normalise_profile_x) by simply dividing the profile size by its max value.

2. The steps between the normalized points is then done by [`calculate_profile_dx()`](./plasma_profiles_abstract_class.md#calculate-the-profile-steps-in-x--calculate_profile_dx) which divides the max x-dimension by the number of points.

3. The core electron and ion temperatures are claculated via [`set_physics_variables()`]() 


    ### Calculate core values | `set_physics_variables()`

    The core electron density is calculated using the `ncore` method.
    The core ion density is then set from the output out `dnitot` which is the total ion density such as:

    $$
    n_{i,0} = \left(\frac{n_i}{n_e}\right)n_{e,0}
    $$
    
    #### Electron core density of a pedestalised profile | `ncore()`

    This function calculates the core electron density for a pedestalsied profile (`ipedestal == 1`). It takes in values of 

    | Profile parameter / Input               | Density   |
    |----------------------------------|-----------|
    | Pedestal radius (r/a)            | `rhopedn`, $\rho_{\text{ped,n}}$ |
    | Pedestal value                   | `neped`, $n_{\text{ped}}$ |
    | Separatrix value                 | `nesep`, $n_{\text{sep}}$ |
    | Average density             | `dene`, $\langle n \rangle$ |
    | Profile index/ peaking parameter | `alphan`, $\alpha_n$ |


    $$\begin{aligned}
        \nonumber
        n_0 & = & \frac{1}{3\rho_{\text{ped,n}}^2} \left[3\langle n\rangle (1+\alpha_n)
            + n_{\text{sep}} (1+\alpha_n) (-2 + \rho_{\text{ped,n}} + \rho_{\text{ped,n}}^2) \right.\\
        & & \left. - n_{\text{ped}}\left( (1 + \alpha_n)(1+ \rho_{\text{ped,n}}) + (\alpha_n -2)
            \rho_{\text{ped,n}}^2 \right) \right]
    \end{aligned}$$

    If `ncore` is returned as being less that 0, it is forced into a state of `ncore = 1E-6` in order to help convergance. This will also give a warning to the user to raise the lower bound of the average electron density `dene`.

    ##### Derivation

    Taking the standard model for the density profile and calculating the sums of the [volumes of integration](https://en.wikipedia.org/wiki/Solid_of_revolution) for the main profile and pedestal regions will give us a density equal to the average plasma density, $\langle n \rangle$ (`dene`). Calculating this and then re-arranging to find $n_0$ allows us to find the central density
    
    
    $$
    \langle n \rangle = \int_{0}^{\rho_{\text{ped},n}}2\pi\rho \left[n_{\text{ped}} + (n_0 - n_{\text{ped}}) \left( 1 -
    \frac{\rho^2}{\rho_{\text{ped},n}^2}\right)^{\alpha_n}\right] d\rho \\
    + \int_{\rho_{\text{ped},n}}^{1} 2\pi\rho \left[n_{\text{sep}} + (n_{\text{ped}} - n_{\text{sep}})\left( \frac{1- \rho}{1-\rho_{\text{ped},n}}\right)\right] d\rho
    $$

    $$
    \langle n \rangle = \frac{(n_0+n_{\text{ped}}\alpha_\text{n})\rho_{\text{ped}}^2}{1+\alpha_n}+\frac{1}{3}(1-\rho_{\text{ped}})(n_{\text{ped}}+2n_{\text{ped}}\rho_{\text{ped}}+n_{\text{sep}}(2+\rho_{\text{ped}}))
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