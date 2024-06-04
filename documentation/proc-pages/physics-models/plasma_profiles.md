# Plasma Profiles

!!! warning " Un-realistic profiles"

    If `ipedestal >= 1` it is highly recommended to use constraint equation 81 (icc=81). This enforces solutions in which $n_0$ has to be greater than $n_{ped}$. 
    Negative $n_0$ values can also arise during iteration, so it is important to be weary on how low the lower bound for $n_e (\mathtt{dene})$ is set.

The density and temperature profile of the plasma for electrons and ions can take two forms depending on the switch value for `ipedestal`. Either without a [pedestal](http://fusionwiki.ciemat.es/wiki/Pedestal) ,`ipedestal = 0` or with a pedestal `ipedestal = 1`.

The files responsible for calculting and storing the profiles are `plasma_profiles.py` and `profiles.py`.

## Initialization
The plasma profile class is `PlasmaProfile`. Initialization sets the profile class size and `neprofile` and `teprofile` to `NProfile` & `TProfile` from `profiles`

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


### Density `NProfile()`

1. Firstly the profile x-dimension is normalised in `normalise_profile_x()` by simply dividing the profile size by its max value

2. The steps between the normalized points is then done by `calculate_profile_dx()` which divided the max x-dimension by the number of points.

3. `set_physics_variables()` is then ran which performs `ncore()` which calculates the central electron density. The ion central density is then calculated by the ratio from this scaling.

-------------------------------------

#### `ncore`

$$\begin{aligned}
  \nonumber
  n_0 & = & \frac{1}{3\rho_{ped,n}^2} \left[3\langle n\rangle (1+\alpha_n)
    + n_{sep} (1+\alpha_n) (-2 + \rho_{ped,n} + \rho_{ped,n}^2) \right.\\
   & & \left. - n_{ped}\left( (1 + \alpha_n)(1+ \rho_{ped,n}) + (\alpha_n -2)
    \rho_{ped,n}^2 \right) \right]
\end{aligned}$$

where $\rho = r/a$, and $a$ is the plasma minor radius. This gives
volume-averaged values $\langle n \rangle = n_0 / (1+\alpha_n)$, and
line-averaged values $\bar{n} \sim n_0 / \sqrt{(1+\alpha_n)}$, etc.  These
volume- and line-averages are used throughout the code along with the profile
indices $\alpha$, in the various physics models, many of which are fits to
theory-based or empirical scalings. Thus, the plasma model in PROCESS may
be described as 1/2-D.  The relevant profile index variables are
`alphan`, `alphat` and `alphaj`, respectively.

Subscripts $0$, $ped$ and $sep$, denote values at the centre ($\rho = 0$), the
pedestal ($\rho = \rho_{ped}$) and the separatrix ($\rho=1$),
respectively. The density and temperature peaking parameters $\alpha_n$ and a
$\alpha_T$ as well as the second exponent $\beta_T$ (input parameter
`tbeta`, not to be confused with the plasma beta) in the temperature
profile can be chosen by the user, as can the pedestal heights and the values
at the separatrix (`neped, nesep` for the electron density, and
`teped, tesep` for the electron temperature); the ion equivalents are
scaled from the electron values by the ratio of the volume-averaged values).

-------------------------------------------------


$$
n_{i,0} = \left(
            \frac{n_{i,tot}}{n_e} n_{e,0}
        \right)
$$

4. The y profile is then calculated using `calculate_profile_y()`. This routine calculates the density at each normalised minor radius position $\rho$ for a HELIOS-type density pedestal profile (nprofile)[^3]

If `ipedestal == 0` then the original parabolic profile form is used





<!DOCTYPE html>
<html lang="en">

  <head>

<meta charset="utf-8">
<title>Bokeh Plot</title>







<script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-2.4.0.min.js"></script>
<script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-widgets-2.4.0.min.js"></script>
<script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-mathjax-2.4.0.min.js"></script>
<script type="text/javascript">
    Bokeh.set_log_level("info");
</script>
        
      
      

  </head>


  <body>

      
        
          
          
            
        <div class="bk-root" id="8fce6150-bf8e-4367-9dd2-7b71a9f9c6f2" data-root-id="1040"></div>





<script type="application/json" id="1151">
    {"1ddcbc33-1892-4eba-ae9b-4e27540773e4":{"defs":[],"roots":{"references":[{"attributes":{"coordinates":null,"group":null,"text":"Parabolic plasma profile"},"id":"1004","type":"Title"},{"attributes":{},"id":"1043","type":"BasicTickFormatter"},{"attributes":{},"id":"1024","type":"WheelZoomTool"},{"attributes":{"below":[{"id":"1014"}],"center":[{"id":"1017"},{"id":"1021"}],"height":400,"left":[{"id":"1018"}],"renderers":[{"id":"1034"}],"title":{"id":"1004"},"toolbar":{"id":"1026"},"toolbar_location":"below","width":400,"x_range":{"id":"1006"},"x_scale":{"id":"1010"},"y_range":{"id":"1008"},"y_scale":{"id":"1012"}},"id":"1003","subtype":"Figure","type":"Plot"},{"attributes":{},"id":"1044","type":"AllLabels"},{"attributes":{"source":{"id":"1002"}},"id":"1035","type":"CDSView"},{"attributes":{"line_alpha":0.2,"line_color":"#1f77b4","line_width":3,"x":{"field":"x"},"y":{"field":"y"}},"id":"1033","type":"Line"},{"attributes":{},"id":"1046","type":"BasicTickFormatter"},{"attributes":{"tools":[{"id":"1022"},{"id":"1023"},{"id":"1024"}]},"id":"1026","type":"Toolbar"},{"attributes":{"bottom_units":"screen","coordinates":null,"fill_alpha":0.5,"fill_color":"lightgrey","group":null,"left_units":"screen","level":"overlay","line_alpha":1.0,"line_color":"black","line_dash":[4,4],"line_width":2,"right_units":"screen","syncable":false,"top_units":"screen"},"id":"1025","type":"BoxAnnotation"},{"attributes":{},"id":"1047","type":"AllLabels"},{"attributes":{},"id":"1019","type":"BasicTicker"},{"attributes":{"axis":{"id":"1018"},"coordinates":null,"dimension":1,"group":null,"ticker":null},"id":"1021","type":"Grid"},{"attributes":{"axis_label":"$$n_e$$","coordinates":null,"formatter":{"id":"1043"},"group":null,"major_label_policy":{"id":"1044"},"ticker":{"id":"1019"}},"id":"1018","type":"LinearAxis"},{"attributes":{"axis":{"id":"1014"},"coordinates":null,"group":null,"ticker":null},"id":"1017","type":"Grid"},{"attributes":{},"id":"1012","type":"LinearScale"},{"attributes":{},"id":"1010","type":"LinearScale"},{"attributes":{},"id":"1015","type":"BasicTicker"},{"attributes":{"children":[{"id":"1036"},{"id":"1037"}]},"id":"1039","type":"Column"},{"attributes":{"line_alpha":0.1,"line_color":"#1f77b4","line_width":3,"x":{"field":"x"},"y":{"field":"y"}},"id":"1032","type":"Line"},{"attributes":{},"id":"1048","type":"UnionRenderers"},{"attributes":{},"id":"1049","type":"Selection"},{"attributes":{"coordinates":null,"data_source":{"id":"1002"},"glyph":{"id":"1031"},"group":null,"hover_glyph":null,"muted_glyph":{"id":"1033"},"nonselection_glyph":{"id":"1032"},"view":{"id":"1035"}},"id":"1034","type":"GlyphRenderer"},{"attributes":{"data":{"x":{"__ndarray__":"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","dtype":"float64","order":"little","shape":[500]},"y":{"__ndarray__":"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","dtype":"float64","order":"little","shape":[500]}},"selected":{"id":"1049"},"selection_policy":{"id":"1048"}},"id":"1002","type":"ColumnDataSource"},{"attributes":{"end":10.0,"js_property_callbacks":{"change:value":[{"id":"1038"}]},"start":0.05,"step":0.025,"title":"Decay index","value":0.5},"id":"1037","type":"Slider"},{"attributes":{"args":{"alphan":{"id":"1037"},"n0":{"id":"1036"},"source":{"id":"1002"}},"code":"\n    const A = n0.value\n    const B = alphan.value\n    const x = source.data.x\n    const y = Array.from(x, (x) =&gt;  A*(1-x**2)**B)\n    source.data = { x, y }\n"},"id":"1038","type":"CustomJS"},{"attributes":{"overlay":{"id":"1025"}},"id":"1023","type":"BoxZoomTool"},{"attributes":{"axis_label":"$$\\rho / a$$","coordinates":null,"formatter":{"id":"1046"},"group":null,"major_label_policy":{"id":"1047"},"ticker":{"id":"1015"}},"id":"1014","type":"LinearAxis"},{"attributes":{"end":7},"id":"1008","type":"Range1d"},{"attributes":{"line_alpha":0.6,"line_color":"#1f77b4","line_width":3,"x":{"field":"x"},"y":{"field":"y"}},"id":"1031","type":"Line"},{"attributes":{"children":[{"id":"1003"},{"id":"1039"}]},"id":"1040","type":"Row"},{"attributes":{},"id":"1022","type":"PanTool"},{"attributes":{"end":6.0,"js_property_callbacks":{"change:value":[{"id":"1038"}]},"start":0.1,"step":0.05,"title":"$$\\delta \\text{ (damping factor, 1/s)}$$","value":1.0},"id":"1036","type":"Slider"},{"attributes":{"end":1.1},"id":"1006","type":"Range1d"}],"root_ids":["1040"]},"title":"Bokeh Application","version":"2.4.0"}}
</script>
<script type="text/javascript">
    (function() {
    const fn = function() {
        Bokeh.safely(function() {
        (function(root) {
            function embed_document(root) {

            const docs_json = document.getElementById('1151').textContent;
            const render_items = [{"docid":"1ddcbc33-1892-4eba-ae9b-4e27540773e4","root_ids":["1040"],"roots":{"1040":"8fce6150-bf8e-4367-9dd2-7b71a9f9c6f2"}}];
            root.Bokeh.embed.embed_items(docs_json, render_items);

            }
            if (root.Bokeh !== undefined) {
            embed_document(root);
            } else {
            let attempts = 0;
            const timer = setInterval(function(root) {
                if (root.Bokeh !== undefined) {
                clearInterval(timer);
                embed_document(root);
                } else {
                attempts++;
                if (attempts > 100) {
                    clearInterval(timer);
                    console.log("Bokeh: ERROR: Unable to run BokehJS code because BokehJS library is missing");
                }
                }
            }, 10, root)
            }
        })(window);
        });
    };
    if (document.readyState != "loading") fn();
    else document.addEventListener("DOMContentLoaded", fn);
    })();
</script>

  </body>

</html>



$$
n(\rho) = n_0(1 - \rho^2)^{\alpha_n} 
$$

The central density ($n_0$) is then checked to make sure it is not less than the pedestal density, $n_{ped}$    

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
        

5. Profile is then integrated with `integrate_profile_y()` using simpsons integration


Integrate profile_y values using scipy.integrate.simpson() function.
"""
self.profile_integ = integrate.simpson(
    self.profile_y, x=self.profile_x, dx=self.profile_dx
)

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