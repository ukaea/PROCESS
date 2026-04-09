# Plasma Radiation

## Overview

In `PROCESS` we have two distinct plasma radiation contributions that are calculated separately. These are the [synchrotron radiation](#synchrotron-radiation--psync_albajar_fidone) and then the [line and Bremsstrahlung radiation](#line--bremsstrahlung-radiation) combined.

By changing the input parameter `radius_plasma_core_norm`, the user may set the normalised radius defining the 'core' region of the plasma. Anything past this normalised radius is known as the plasma 'edge'. **The synchrotron radiation loss is always deemed to come from the 'core' region.**


----------------


## Line & Bremsstrahlung Radiation

The radiation per unit volume is determined using loss functions computed by the `ADAS405` code [^1].
The effective collisional–radiative coefficients necessary to determine the ionization state and radiative losses of each ionic species, 
assuming equilibrium ionization balance in an optically thin plasma, were sourced from `ADF11`-derived data files [^2]. 
These coefficients utilize the generalized collisional-radiative approach [^3] for elements such as $\text{H, He, Be, C, N, O, Ne,}$ and $\text{Si}$. 

For $\text{Ni}$, the data is based on [^4], for $\text{Fe}$ it is derived from [^5]; and for $\text{W,}$ the data is obtained from [^6].
The $\text{Ni}$ and $\text{Fe}$ rates incorporate a density dependence as described in [^7]. For $\text{Kr}$ and $\text{Xe}$, data is taken from the ADAS baseline.

The computed loss functions exhibit a weak dependence on density but are evaluated at a fixed electron density of $10^{19} \text{m}^{-3}$.
This differs from strict coronal equilibrium, which assumes density independence. 
In practice, non-local effects arising from density and temperature gradients are significant but are not considered here. 
The loss functions account for Bremsstrahlung, line radiation, and recombination radiation, represented by:

$$
P_i = n_i n_e L_Z (Z_i, T)
$$

where $P_i$ is the radiation per unit volume (excluding synchrotron radiation),
$L_Z (Z_i, T)$ is the loss function for ion species $i$ at temperature $T$, and $n_i$ is the density of ion species $i$.

The source data comprises of 200 separate data points with a fitted temperature range of $1 \ \text{eV}$ tp $40 \ \text{keV}$. The fitted radiation profiles for each species against temperature can be found in the figure below.




<figure markdown>
![ADAS Radiation](../images/adas_radiation.png){ width = "100"}
<figcaption>Figure 1: Radiation loss functions as a function of temperature, at 10^19 electrons m^−3.
The lowest line is H, with the other lines in the order listed. The dashed lines show
the bremsstrahlung calculated using a separate method.</figcaption>
</figure>

For the regime above $40 \ \text{keV}$ whichc an be seen by the dashed lines in Figure 1 above, the radiation is assumed to be Bhremsstrahlung dominated so the values are extrapolated via a power law.


!!! note "Location of impurities"

    All species/impurities are currently assumed to be distributed at a fixed relative value to the electron density profile throughout the plasma. So the electron density and temperature profiles are integrated using the values of `f_nd_impurity_electron()` to get the line and Bhremmsstrahlung radiation of a species across the plasma profile. 
    

!!! note "Particle confinement"
    For the loss function $L_Z (Z_i, T)$, the particles are assumed to be in equailibrium and thus have infinite confinement time. (This is the data found at the bottom of the ADAS derived files in `process/data/lz_non_corona_14_elements`)

    This assumes:

    - There are no particle losses to the walls of the confinement vessel or through magnetic field lines.

    - The system exists in a perfectly steady state where the ionization and recombination rates precisely balance, allowing the coronal equilibrium model to apply perfectly. 

    - The radiation is "optically thin," meaning that the emitted photons can escape the plasma without being re-absorbed by other particles. 

The radiation emission is numerically integrated over the plasma profile, 
using the corresponding temperature and density distributions. Emission is only considered from within the separatrix, 
as the PROCESS model does not account for the scrape-off layer. 
The plasma inside the separatrix is divided into two regions: the “core” and the “edge,” separated by a normalised minor radius defined by the user. 
Radiation is calculated separately for the core and edge regions, except for synchrotron radiation, which is assumed to originate solely from the core.

--------------

### Reduction of core radiation

Work by H.Lux et al.[^8] has showed that reducing the total core contribution of the radiation from that predicted by the raw ADAS data allows better fitting to the $\text{IPB98(y,2)}$ scaling for actual plasma stored energies calculated in the `ASTRA/TGLF` codes.

Setting of the value `f_p_plasma_core_rad_reduction` will cause each radiated power density value in the core profile to be multiplied by this value. This in essence allows a scaling of the core radiated power to be higher or lower than the raw ADAS data.

------------

## Synchrotron radiation | `psync_albajar_fidone()`

The formula below is the current synchrotron radiation loss power implemented in `PROCESS`. It is a combination of the general form given by Albajar et al. [^9] and a wall reflectivity correction given by Fidone et al. [^10]

$$
\begin{aligned}
P_{\text{syn}}(\mathrm{MW})= & 3.84 \times 10^{-8}\times \frac{(1-f_{\text{reflect}})^{0.62}}{\left[1+0.12 \frac{T_{\text{e0}}}{p_{a_0}^{0.41}}\left(1.0 - f_{\text{reflect}}\right)^{0.41}\right]^{1.51}}  \\
& \times R a^{1.38} \kappa^{0.79}B_{\text{T}}^{2.62} n_{\text{e0,20}}^{0.38} T_{\text{e0}}\left(16+T_{\text{e0}}\right)^{2.61} \\
& \times K\left(\alpha_n, \alpha_T, \beta_T\right) \\
& \times G(A)
\end{aligned}
$$

$$
K(\alpha_n, \alpha_T, \beta_T) = \left(\alpha_n +3.87\alpha_T +1.46\right)^{-0.79} \\
\times \left(1.98+\alpha_T\right)^{1.36}\beta_T^{2.14} \\
\times \left(\beta_T^{1.53}+1.87\alpha_T-0.16\right)^{-1.33}
$$

$$
p_{a_0} = 6.04 \times 10^3 \frac{a n_{\text{e0(20)}}}{B_{\text{T}}}
$$

where $T_{\text{e0}}$ is the central electron temperature in keV, $R$ is the plasma major radius, $a$ is the plasma minor radius, $\kappa$ is the plasma separatrix elongation, $B_{\text{T}}$ is the on axis toroidal magnetic field, $n_{\text{e0,20}}$ is the central electron density in units of $10^{20} \text{m}^{-3}$, $\alpha_n$ is the density profile peaking parameter, $\alpha_T$ is the temperature profile peaking parameter, $\beta_T$ is the [secondary temperature profile peaking parameter](./profiles/plasma_profiles.md#pedestal-profile--h-mode) and $A$ is the plasma aspect ratio.

The original form of the synchrotron radiation formula presented by Albajar et al. [^9] is based off using multiple non-linear regression from a database consisting of 3000 complete computations of the analytical expression for synchrotron power derived by Albajar et al.[^9] to include aspect ratio and temperature profile dependence.

The fitting variable range is:

- $10 < T_{\text{e0}} < 100 \ \text{keV}$
- $10^2 < p_{\text{a0}} < 1 \times 10^4$
- $1 < \kappa < 2.5$
- $0 < \alpha_n < 2$
- $0 < \alpha_T < 8$
- $1 < \beta_T < 8$


The root mean square error of the fit is found to be 5.8%

The original form also uses a fair estimation of the effect of a wall with
a reflection coefficient $f_{\text{reflect}}$ obtained from the Trubnikov approach [^11].

$$
P_{\text{syn,r}} \propto \left(1 − f_{\text{reflect}}\right)^{1/2} P_{\text{syn}}
$$

It should also be noted that the wall reflection coefficient for the synchrotron radiation is poorly known.

- The wall reflection factor ($f_{\text{reflect}}$) may be set by the user by inputting `f_sync_reflect = <value>`.

The wall reflectivity correction presented later by Fidone et al. [^10] is of the form:

$$
P_{\text{syn,r}} \propto \frac{\left(1 − f_{\text{reflect}}\right)^{0.62}}{\left[1+\left(1 − f_{\text{reflect}}\right)^{0.41}\right]^{1.51}} P_{\text{syn}}
$$

A comparison of the different reflection functions can be seen in the graph below:

<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>Normalised synchrotron radiation loss vs wall reflectivity</title>
    <style>
      html, body {
        box-sizing: border-box;
        display: flow-root;
        height: 100%;
        margin: 0;
        padding: 0;
      }
    </style>
    <script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-3.6.0.min.js"></script>
    <script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-mathjax-3.6.0.min.js"></script>
    <script type="text/javascript">
        Bokeh.set_log_level("info");
    </script>
  </head>
  <body>
    <div id="eb800ba2-3e84-4f46-932d-71bda79eab80" data-root-id="p1007" style="display: contents;"></div>
  
    <script type="application/json" id="cbf88062-06c9-46fe-8dfe-451e6ae8e999">
      {"48247fbd-6fdc-41d7-b5a8-5e857f114784":{"version":"3.6.0","title":"Bokeh Application","roots":[{"type":"object","name":"Figure","id":"p1007","attributes":{"width":400,"height":400,"x_range":{"type":"object","name":"Range1d","id":"p1017"},"y_range":{"type":"object","name":"Range1d","id":"p1018"},"x_scale":{"type":"object","name":"LinearScale","id":"p1019"},"y_scale":{"type":"object","name":"LinearScale","id":"p1020"},"title":{"type":"object","name":"Title","id":"p1010","attributes":{"text":"Normalised synchrotron loss vs wall reflectivity"}},"renderers":[{"type":"object","name":"GlyphRenderer","id":"p1050","attributes":{"data_source":{"type":"object","name":"ColumnDataSource","id":"p1001","attributes":{"selected":{"type":"object","name":"Selection","id":"p1002","attributes":{"indices":[],"line_indices":[]}},"selection_policy":{"type":"object","name":"UnionRenderers","id":"p1003"},"data":{"type":"map","entries":[["x",{"type":"ndarray","array":{"type":"bytes","data":"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"},"shape":[500],"dtype":"float64","order":"little"}],["y",{"type":"ndarray","array":{"type":"bytes","data":"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"},"shape":[500],"dtype":"float64","order":"little"}]]}}},"view":{"type":"object","name":"CDSView","id":"p1051","attributes":{"filter":{"type":"object","name":"AllIndices","id":"p1052"}}},"glyph":{"type":"object","name":"Line","id":"p1047","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"blue","line_alpha":0.6,"line_width":3}},"nonselection_glyph":{"type":"object","name":"Line","id":"p1048","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"blue","line_alpha":0.1,"line_width":3}},"muted_glyph":{"type":"object","name":"Line","id":"p1049","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"blue","line_alpha":0.2,"line_width":3}}}},{"type":"object","name":"GlyphRenderer","id":"p1061","attributes":{"data_source":{"type":"object","name":"ColumnDataSource","id":"p1004","attributes":{"selected":{"type":"object","name":"Selection","id":"p1005","attributes":{"indices":[],"line_indices":[]}},"selection_policy":{"type":"object","name":"UnionRenderers","id":"p1006"},"data":{"type":"map","entries":[["x",{"type":"ndarray","array":{"type":"bytes","data":"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"},"shape":[500],"dtype":"float64","order":"little"}],["y2",{"type":"ndarray","array":{"type":"bytes","data":"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"},"shape":[500],"dtype":"float64","order":"little"}]]}}},"view":{"type":"object","name":"CDSView","id":"p1062","attributes":{"filter":{"type":"object","name":"AllIndices","id":"p1063"}}},"glyph":{"type":"object","name":"Line","id":"p1058","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y2"},"line_color":"red","line_alpha":0.6,"line_width":3}},"nonselection_glyph":{"type":"object","name":"Line","id":"p1059","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y2"},"line_color":"red","line_alpha":0.1,"line_width":3}},"muted_glyph":{"type":"object","name":"Line","id":"p1060","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y2"},"line_color":"red","line_alpha":0.2,"line_width":3}}}}],"toolbar":{"type":"object","name":"Toolbar","id":"p1016","attributes":{"tools":[{"type":"object","name":"PanTool","id":"p1031"},{"type":"object","name":"WheelZoomTool","id":"p1032","attributes":{"renderers":"auto"}},{"type":"object","name":"BoxZoomTool","id":"p1033","attributes":{"overlay":{"type":"object","name":"BoxAnnotation","id":"p1034","attributes":{"syncable":false,"line_color":"black","line_alpha":1.0,"line_width":2,"line_dash":[4,4],"fill_color":"lightgrey","fill_alpha":0.5,"level":"overlay","visible":false,"left":{"type":"number","value":"nan"},"right":{"type":"number","value":"nan"},"top":{"type":"number","value":"nan"},"bottom":{"type":"number","value":"nan"},"left_units":"canvas","right_units":"canvas","top_units":"canvas","bottom_units":"canvas","handles":{"type":"object","name":"BoxInteractionHandles","id":"p1040","attributes":{"all":{"type":"object","name":"AreaVisuals","id":"p1039","attributes":{"fill_color":"white","hover_fill_color":"lightgray"}}}}}}}},{"type":"object","name":"SaveTool","id":"p1041"},{"type":"object","name":"ResetTool","id":"p1042"},{"type":"object","name":"HelpTool","id":"p1043"}]}},"left":[{"type":"object","name":"LinearAxis","id":"p1026","attributes":{"ticker":{"type":"object","name":"BasicTicker","id":"p1027","attributes":{"mantissas":[1,2,5]}},"formatter":{"type":"object","name":"BasicTickFormatter","id":"p1028"},"axis_label":"Normalized synchrotron radiation loss","major_label_policy":{"type":"object","name":"AllLabels","id":"p1029"}}}],"below":[{"type":"object","name":"LinearAxis","id":"p1021","attributes":{"ticker":{"type":"object","name":"BasicTicker","id":"p1022","attributes":{"mantissas":[1,2,5]}},"formatter":{"type":"object","name":"BasicTickFormatter","id":"p1023"},"axis_label":"Wall reflectivity, \\  $$[f_{\\text{reflect}}]$$","major_label_policy":{"type":"object","name":"AllLabels","id":"p1024"}}}],"center":[{"type":"object","name":"Grid","id":"p1025","attributes":{"axis":{"id":"p1021"}}},{"type":"object","name":"Grid","id":"p1030","attributes":{"dimension":1,"axis":{"id":"p1026"}}},{"type":"object","name":"Legend","id":"p1053","attributes":{"location":"bottom_left","title":"Legend","items":[{"type":"object","name":"LegendItem","id":"p1054","attributes":{"label":{"type":"value","value":"Fidone - Correction"},"renderers":[{"id":"p1050"}]}},{"type":"object","name":"LegendItem","id":"p1064","attributes":{"label":{"type":"value","value":"Albajar - Original"},"renderers":[{"id":"p1061"}]}}]}}]}}]}}
</script>
<script type="text/javascript">
    (function() {
    const fn = function() {
        Bokeh.safely(function() {
        (function(root) {
            function embed_document(root) {
            const docs_json = document.getElementById('cbf88062-06c9-46fe-8dfe-451e6ae8e999').textContent;
            const render_items = [{"docid":"48247fbd-6fdc-41d7-b5a8-5e857f114784","roots":{"p1007":"eb800ba2-3e84-4f46-932d-71bda79eab80"},"root_ids":["p1007"]}];
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

The power loss due to synchrotron radiation grows as the aspect ratio decreases. At high aspect ratios $(A > 6)$, although $P_{\text{syn}}$ increases with $R$, the normalised synchrotron loss saturates. This is due to the fact that the magnetic field inhomogeneity vanishes for large $A$. For the above reason a $G$ correction factor is implemented.This gives an root mean square error of 6.2% with respect to a secondary
dataset consisting of 640 complete computations for the same range of variables variables above and for an aspect ratio interval of $1.5 <A< 15$.

$$
G\left(A\right) = 0.93\left[1+0.85 e^{-0.82A}\right]
$$



-------------

## Impurity Radiation Class | `ImpurityRadiation`

### Initialization | `__init__()`

Initialize the FusionReactionRate class with the given plasma profile.

#### Parameters:
- `plasma_profile (PlasmaProfile)`: The parameterized temperature and density profiles of the plasma. Taken from the plasma_profile object.

#### Attributes:
- `plasma_profile` (`PlasmaProfile`): Plasma profile instance.
- `rho` (`numpy.ndarray`): Density profile along the x-axis.
- `rhodx` (`numpy.ndarray`): Density profile step size along the x-axis.
- `imp` (`numpy.ndarray`): Indices of impurities with a fraction greater than 1.0e-30.
- `pimp_profile` (`numpy.ndarray`): Impurity profile array initialized to zeros.
- `pden_impurity_rad_profile` (`numpy.ndarray`): Total radiation profile array initialized to zeros.
- `pden_impurity_core_rad_profile` (`numpy.ndarray`): Core radiation profile array initialized to zeros.
- `pden_impurity_rad_total` (`float`): Total radiation initialized to 0.0.
- `pden_impurity_core_rad_total` (`float`): Core radiation initialized to 0.0.



--------------------



## Key Constraints

### Radiation power density upper limit

This constraint can be activated by stating `icc = 17` in the input file.

$$
\frac{P_{\text{rad}}}{V_{\text{plasma}}} < \frac{P_{\text{inj}} + P_{\alpha}f_{\alpha, \text{plasma}} + P_{\text{charged}} + P_{\text{ohmic}} }{V_{\text{plasma}}}
$$

Ensures that the calculated total radiation power density does not exceed the total heating power to the plasma (i.e. the sum of the fusion
alpha power coupled to the plasma, other charged particle fusion power, auxiliary injected power and
the ohmic heating power). If the radiation power is higher than the heating power then the change in the stored plasma thermal energy is negative, $-\frac{\mathrm{d}W}{\mathrm{d}T}$. This means the plasma is cooling and is thus not a viable plasma solution.

$f_{\alpha, \text{plasma}}$ is the fraction of alpha power that is coupled to the plasma (`f_p_alpha_plasma_deposited`).

**It is recommended to have this constraint on as it is a plasma stability model**

----------------

### Radiation wall load upper limit

This constraint can be activated by stating `icc = 67` in the input file.

The limiting value of $q_{\text{fw,rad}}$ in $\mathrm {MWm^{-2}}$ is be set using input parameter `pflux_fw_rad_max`.

[^1]: “ADAS: Documentation,” Adas.ac.uk, 2024. https://www.adas.ac.uk/manual.php
[^2]: “OPEN-ADAS,” Adas.ac.uk, 2025. https://open.adas.ac.uk/adf11 (accessed Jan. 15, 2025).
[^3]: H. P. Summers et al., “Ionization state, excited populations and emission of impurities in dynamic finite density plasmas: I. The generalized collisional–radiative model for light elements,” Plasma Physics and Controlled Fusion, vol. 48, no. 2, pp. 263–293, Jan. 2006, doi: https://doi.org/10.1088/0741-3335/48/2/007.
[^4]: M. Arnaud and R. Rothenflug, “An updated evaluation of recombination and ionization rates,” Astronomy & Astrophysics Supplement Series, vol. 60, no. 3, pp. 425–457, Jun. 1985.
[^5]: M. Arnaud, R. Rothenflug, An updated evaluation of recombination and ionization rates, Astron. & Astrophys. Supp. Ser. 60 (1985) 425–457.
[^6]: T. Pütterich, R. Neu, R. Dux, A. D. Whiteford, and M. G. O’Mullane, “Modelling of measured tungsten spectra from ASDEX Upgrade and predictions for ITER,” Plasma Physics and Controlled Fusion, vol. 50, no. 8, p. 085016, Jun. 2008, doi: https://doi.org/10.1088/0741-3335/50/8/085016.
[^7]: Summers, H. P. (1974). Tables and Graphs of Collisional Dielectronic Recombination and Ionisation Coefficients and Ionisation Equilibria of H-like to A-like Ions of Elements. Appleton Laboratory.
[^8]: H. Lux, R. Kemp, E. Fable, and R. Wenninger, “Radiation and confinement in 0D fusion systems codes,” Plasma Physics and Controlled Fusion, vol. 58, no. 7, pp. 075001–075001, May 2016, doi: https://doi.org/10.1088/0741-3335/58/7/075001.
[^9]: F. Albajar, J. Johner, and G. Granata, “Improved calculation of synchrotron radiation losses in realistic tokamak plasmas,"Nuclear Fusion, vol. 41, no. 6, pp. 665–678, Jun. 2001, doi: https://doi.org/10.1088/0029-5515/41/6/301.
[^10]: I. Fidone, G Giruzzi, and G. Granata, “Synchrotron radiation loss in tokamaks of arbitrary geometry,” Nuclear Fusion, vol. 41, no. 12, pp. 1755–1758, Dec. 2001, doi: https://doi.org/10.1088/0029-5515/41/12/102.
[^11]: Trubnikov, B.A., in Reviews of Plasma Physics, Vol. 7 (Leontovich, M.A., Ed.), Consultants Bureau, New York (1979) 345.
