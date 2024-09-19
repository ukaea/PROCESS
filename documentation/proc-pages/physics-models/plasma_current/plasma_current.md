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

The parameters $\lambda$ and $\nu$ characterise the current profile $J_{\phi} = \frac{j_0}{(1+\lambda r^2)^2}$ and pressure profile $p = p_0(1-r^2)^{\nu}$. The pressure profile used by Connor-Hastie is still of the normal parabolic type.


!!! warning Assumed current profiles

    Since $\lambda$ in this case is treated as the current profile index `alphaj` within PROCESS will use its assumed standard [parabolic profile](../profiles/plasma_profiles.md`). Even though the profile given above by Connor-Hastie is different. The differenc between these two assumed current profiles can be experimented with below.

<!DOCTYPE html>
<html lang="en">
  
  <head>
    
<meta charset="utf-8">
<title>My Bokeh Plot</title>



    


    
<script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-2.4.0.min.js"></script>
<script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-widgets-2.4.0.min.js"></script>
<script type="text/javascript">
    Bokeh.set_log_level("info");
</script>
        
      
      
    
  </head>
  
  
  <body>
    
      
        
    
    
    
    <div class="bk-root" id="f1598c7d-4e49-4714-a81f-6697714facb5" data-root-id="1075"></div>
    
    



<script type="application/json" id="1196">
    {"028b7c7b-7de9-4b2e-80bd-57279ba813bc":{"defs":[],"roots":{"references":[{"attributes":{"line_alpha":0.1,"line_color":"red","line_width":3,"x":{"field":"x"},"y":{"field":"y2"}},"id":"1056","type":"Line"},{"attributes":{"overlay":{"id":"1028"}},"id":"1024","type":"BoxZoomTool"},{"attributes":{},"id":"1025","type":"SaveTool"},{"attributes":{},"id":"1045","type":"AllLabels"},{"attributes":{"coordinates":null,"data_source":{"id":"1002"},"glyph":{"id":"1055"},"group":null,"hover_glyph":null,"muted_glyph":{"id":"1057"},"nonselection_glyph":{"id":"1056"},"view":{"id":"1059"}},"id":"1058","type":"GlyphRenderer"},{"attributes":{"axis":{"id":"1014"},"coordinates":null,"group":null,"ticker":null},"id":"1017","type":"Grid"},{"attributes":{"args":{"alpha":{"id":"1072"},"n0":{"id":"1071"},"source":{"id":"1002"}},"code":"\n    const A = n0.value\n    const B = alpha.value\n\n    const x = source.data.x\n    const y1 = Array.from(x, (x) =&gt; A * (1 - x**2)**B)\n    const y2 = Array.from(x, (x) =&gt; A / (1+B*x**2)**2) // Example transformation for the second line\n    source.data = { x, y1, y2 }\n"},"id":"1073","type":"CustomJS"},{"attributes":{"children":[{"id":"1071"},{"id":"1072"}]},"id":"1074","type":"Column"},{"attributes":{},"id":"1027","type":"HelpTool"},{"attributes":{},"id":"1050","type":"Selection"},{"attributes":{"children":[{"id":"1003"},{"id":"1074"}]},"id":"1075","type":"Row"},{"attributes":{},"id":"1006","type":"Range1d"},{"attributes":{},"id":"1047","type":"BasicTickFormatter"},{"attributes":{"below":[{"id":"1014"}],"center":[{"id":"1017"},{"id":"1021"},{"id":"1052"}],"height":400,"left":[{"id":"1018"}],"renderers":[{"id":"1040"},{"id":"1058"}],"title":{"id":"1004"},"toolbar":{"id":"1029"},"width":400,"x_range":{"id":"1006"},"x_scale":{"id":"1010"},"y_range":{"id":"1008"},"y_scale":{"id":"1012"}},"id":"1003","subtype":"Figure","type":"Plot"},{"attributes":{"axis_label":"Current Density, J","coordinates":null,"formatter":{"id":"1044"},"group":null,"major_label_policy":{"id":"1045"},"ticker":{"id":"1019"}},"id":"1018","type":"LinearAxis"},{"attributes":{"line_alpha":0.6,"line_color":"blue","line_width":3,"x":{"field":"x"},"y":{"field":"y1"}},"id":"1037","type":"Line"},{"attributes":{"line_alpha":0.2,"line_color":"blue","line_width":3,"x":{"field":"x"},"y":{"field":"y1"}},"id":"1039","type":"Line"},{"attributes":{"end":10,"js_property_callbacks":{"change:value":[{"id":"1073"}]},"start":0.1,"step":0.1,"title":"Plasma centre value | j0","value":5},"id":"1071","type":"Slider"},{"attributes":{},"id":"1012","type":"LinearScale"},{"attributes":{"source":{"id":"1002"}},"id":"1041","type":"CDSView"},{"attributes":{},"id":"1044","type":"BasicTickFormatter"},{"attributes":{},"id":"1048","type":"AllLabels"},{"attributes":{"bottom_units":"screen","coordinates":null,"fill_alpha":0.5,"fill_color":"lightgrey","group":null,"left_units":"screen","level":"overlay","line_alpha":1.0,"line_color":"black","line_dash":[4,4],"line_width":2,"right_units":"screen","syncable":false,"top_units":"screen"},"id":"1028","type":"BoxAnnotation"},{"attributes":{"coordinates":null,"data_source":{"id":"1002"},"glyph":{"id":"1037"},"group":null,"hover_glyph":null,"muted_glyph":{"id":"1039"},"nonselection_glyph":{"id":"1038"},"view":{"id":"1041"}},"id":"1040","type":"GlyphRenderer"},{"attributes":{"line_alpha":0.6,"line_color":"red","line_width":3,"x":{"field":"x"},"y":{"field":"y2"}},"id":"1055","type":"Line"},{"attributes":{"end":10,"js_property_callbacks":{"change:value":[{"id":"1073"}]},"start":0.01,"step":0.01,"title":"Profile Index  | alphaj","value":2},"id":"1072","type":"Slider"},{"attributes":{},"id":"1022","type":"PanTool"},{"attributes":{"coordinates":null,"group":null,"items":[{"id":"1053"},{"id":"1070"}]},"id":"1052","type":"Legend"},{"attributes":{"source":{"id":"1002"}},"id":"1059","type":"CDSView"},{"attributes":{"tools":[{"id":"1022"},{"id":"1023"},{"id":"1024"},{"id":"1025"},{"id":"1026"},{"id":"1027"}]},"id":"1029","type":"Toolbar"},{"attributes":{},"id":"1026","type":"ResetTool"},{"attributes":{"label":{"value":"PROCESS Parabolic Current Profile"},"renderers":[{"id":"1040"}]},"id":"1053","type":"LegendItem"},{"attributes":{"axis":{"id":"1018"},"coordinates":null,"dimension":1,"group":null,"ticker":null},"id":"1021","type":"Grid"},{"attributes":{},"id":"1019","type":"BasicTicker"},{"attributes":{"data":{"x":{"__ndarray__":"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","dtype":"float64","order":"little","shape":[500]},"y1":{"__ndarray__":"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","dtype":"float64","order":"little","shape":[500]},"y2":{"__ndarray__":"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","dtype":"float64","order":"little","shape":[500]}},"selected":{"id":"1050"},"selection_policy":{"id":"1049"}},"id":"1002","type":"ColumnDataSource"},{"attributes":{"line_alpha":0.1,"line_color":"blue","line_width":3,"x":{"field":"x"},"y":{"field":"y1"}},"id":"1038","type":"Line"},{"attributes":{"axis_label":"Normalized Plasma Radius","coordinates":null,"formatter":{"id":"1047"},"group":null,"major_label_policy":{"id":"1048"},"ticker":{"id":"1015"}},"id":"1014","type":"LinearAxis"},{"attributes":{"end":10},"id":"1008","type":"Range1d"},{"attributes":{"label":{"value":"Hastie Current Profile"},"renderers":[{"id":"1058"}]},"id":"1070","type":"LegendItem"},{"attributes":{},"id":"1049","type":"UnionRenderers"},{"attributes":{},"id":"1023","type":"WheelZoomTool"},{"attributes":{"line_alpha":0.2,"line_color":"red","line_width":3,"x":{"field":"x"},"y":{"field":"y2"}},"id":"1057","type":"Line"},{"attributes":{},"id":"1010","type":"LinearScale"},{"attributes":{},"id":"1015","type":"BasicTicker"},{"attributes":{"coordinates":null,"group":null,"text":"PROCESS vs Hastie Current Profile"},"id":"1004","type":"Title"}],"root_ids":["1075"]},"title":"Bokeh Application","version":"2.4.0"}}
</script>
<script type="text/javascript">
    (function() {
    const fn = function() {
        Bokeh.safely(function() {
        (function(root) {
            function embed_document(root) {
            
            const docs_json = document.getElementById('1196').textContent;
            const render_items = [{"docid":"028b7c7b-7de9-4b2e-80bd-57279ba813bc","root_ids":["1075"],"roots":{"1075":"f1598c7d-4e49-4714-a81f-6697714facb5"}}];
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

---------------

#### Sauter model, allows negative triangularity | `calculate_current_coefficient_sauter()`

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

#### FIESTA for ST's | `calculate_current_coefficient_fiesta()` 

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

-----------------------

## Key Constraints

-----------------------

### Plasma current ramp-up time lower limit

This constraint can be activated by stating `icc = 41` in the input file.

The value of `tohsm` can be set to the required minimum plasma current ramp up time at the start of a pulse. The scaling value `ftohs` can be varied also

The calculated plasma current ramp up time `tohs` is dictated by the [pulsed plant operation configuration](../pulsed-plant.md).

This constraint will ensure that the value of `tohs` is always greater than or equal to `tohsm`

--------------------

### Plasma and Rod current fraction upper limit

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

