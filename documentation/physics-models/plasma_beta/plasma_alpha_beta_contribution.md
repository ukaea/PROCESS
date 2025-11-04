# Fast Alpha Pressure Contribution | `fast_alpha_beta()`

The pressure contribution from the fast alpha particles can be controlled using switch `i_beta_fast_alpha`.
This sets the value of the physics variable, `beta_fast_alpha`.

A contribution from fast alphas to the total plasma pressure can be relatively high (~ 10-30%) for fusion temperatures of interest Because the maximum volume-averaged beta achievable in a tokamak is limited by  magnetohydrrdynamic (MHD) instabilities (e.g.,ballooning and kink modes), the presence of fast alphas, if $\langle \beta_{\text{tot}} \rangle$ is constant, reduces the background thermal plasma pressure. Furthermore, the energetic alpha population can influence (favorably or unfavorably) the bulk plasma ballooning mode stability boundaries[^2].

A comaprison between the two avaialbe scalings can be experimented with below:

<!DOCTYPE html>
<html lang="en">

<head>

<meta charset="utf-8">
<title>My Bokeh Plot</title>

<script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-2.4.0.min.js"></script>
<script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-widgets-2.4.0.min.js"></script>
<script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-mathjax-2.4.0.min.js"></script>
<script type="text/javascript">
    Bokeh.set_log_level("info");
</script>




</head>


<body>
    
      
        
          
          
            
        <div class="bk-root" id="51d9fa83-acf5-40ed-a2e9-10dad5bd7ad5" data-root-id="1074"></div>
    
    



<script type="application/json" id="1185">
    {"68cf82d8-8d08-40cf-8418-5dcbd5f6848c":{"defs":[],"roots":{"references":[{"attributes":{},"id":"1022","type":"PanTool"},{"attributes":{"overlay":{"id":"1028"}},"id":"1024","type":"BoxZoomTool"},{"attributes":{"source":{"id":"1002"}},"id":"1041","type":"CDSView"},{"attributes":{},"id":"1025","type":"SaveTool"},{"attributes":{},"id":"1023","type":"WheelZoomTool"},{"attributes":{},"id":"1050","type":"Selection"},{"attributes":{"line_alpha":0.2,"line_color":"blue","line_width":3,"x":{"field":"x"},"y":{"field":"y1"}},"id":"1039","type":"Line"},{"attributes":{"axis_label":"\\[\\left(\\frac{n_{\\text{DT}}}{n_{\\text{e}}}\\right)^2\\]","coordinates":null,"formatter":{"id":"1047"},"group":null,"major_label_policy":{"id":"1048"},"ticker":{"id":"1015"}},"id":"1014","type":"LinearAxis"},{"attributes":{},"id":"1019","type":"BasicTicker"},{"attributes":{},"id":"1049","type":"UnionRenderers"},{"attributes":{"axis":{"id":"1018"},"coordinates":null,"dimension":1,"group":null,"ticker":null},"id":"1021","type":"Grid"},{"attributes":{"children":[{"id":"1071"}]},"id":"1073","type":"Column"},{"attributes":{},"id":"1012","type":"LinearScale"},{"attributes":{"data":{"x":{"__ndarray__":"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","dtype":"float64","order":"little","shape":[500]},"y1":{"__ndarray__":"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","dtype":"float64","order":"little","shape":[500]},"y2":{"__ndarray__":"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","dtype":"float64","order":"little","shape":[500]}},"selected":{"id":"1050"},"selection_policy":{"id":"1049"}},"id":"1002","type":"ColumnDataSource"},{"attributes":{"tools":[{"id":"1022"},{"id":"1023"},{"id":"1024"},{"id":"1025"},{"id":"1026"},{"id":"1027"}]},"id":"1029","type":"Toolbar"},{"attributes":{},"id":"1048","type":"AllLabels"},{"attributes":{},"id":"1015","type":"BasicTicker"},{"attributes":{"line_alpha":0.6,"line_color":"red","line_width":3,"x":{"field":"x"},"y":{"field":"y2"}},"id":"1055","type":"Line"},{"attributes":{"axis_label":"\\[\\ \\beta_{\\alpha} / \\beta_{\\text{th}}\\]","coordinates":null,"formatter":{"id":"1044"},"group":null,"major_label_policy":{"id":"1045"},"ticker":{"id":"1019"}},"id":"1018","type":"LinearAxis"},{"attributes":{"line_alpha":0.1,"line_color":"blue","line_width":3,"x":{"field":"x"},"y":{"field":"y1"}},"id":"1038","type":"Line"},{"attributes":{"line_alpha":0.1,"line_color":"red","line_width":3,"x":{"field":"x"},"y":{"field":"y2"}},"id":"1056","type":"Line"},{"attributes":{"line_alpha":0.2,"line_color":"red","line_width":3,"x":{"field":"x"},"y":{"field":"y2"}},"id":"1057","type":"Line"},{"attributes":{},"id":"1047","type":"BasicTickFormatter"},{"attributes":{"children":[{"id":"1003"},{"id":"1073"}]},"id":"1074","type":"Row"},{"attributes":{"axis":{"id":"1014"},"coordinates":null,"group":null,"ticker":null},"id":"1017","type":"Grid"},{"attributes":{"coordinates":null,"group":null,"items":[{"id":"1053"},{"id":"1070"}]},"id":"1052","type":"Legend"},{"attributes":{"args":{"source":{"id":"1002"},"te":{"id":"1071"}},"code":"\n    const A = te.value\n\n    const x = source.data.x\n    const y1 = Array.from(x, (x) =&gt; 0.29 * (A/20 - 0.37) * (x)**2) // Example transformation for the first line\n    const y2 = Array.from(x, (x) =&gt; 0.26 * (A/20 - 0.65)**0.5 * (x)**2) // Example transformation for the second line\n    source.data = { x, y1, y2 }\n"},"id":"1072","type":"CustomJS"},{"attributes":{"end":40,"js_property_callbacks":{"change:value":[{"id":"1072"}]},"start":0.1,"step":0.1,"title":"Density weighted temperature [keV] | T","value":15},"id":"1071","type":"Slider"},{"attributes":{"source":{"id":"1002"}},"id":"1059","type":"CDSView"},{"attributes":{"label":{"value":"Ward Model"},"renderers":[{"id":"1058"}]},"id":"1070","type":"LegendItem"},{"attributes":{},"id":"1026","type":"ResetTool"},{"attributes":{"label":{"value":"IPDG89 Model"},"renderers":[{"id":"1040"}]},"id":"1053","type":"LegendItem"},{"attributes":{"coordinates":null,"data_source":{"id":"1002"},"glyph":{"id":"1055"},"group":null,"hover_glyph":null,"muted_glyph":{"id":"1057"},"nonselection_glyph":{"id":"1056"},"view":{"id":"1059"}},"id":"1058","type":"GlyphRenderer"},{"attributes":{},"id":"1045","type":"AllLabels"},{"attributes":{},"id":"1044","type":"BasicTickFormatter"},{"attributes":{"end":5},"id":"1006","type":"Range1d"},{"attributes":{},"id":"1008","type":"Range1d"},{"attributes":{"below":[{"id":"1014"}],"center":[{"id":"1017"},{"id":"1021"},{"id":"1052"}],"height":400,"left":[{"id":"1018"}],"renderers":[{"id":"1040"},{"id":"1058"}],"title":{"id":"1004"},"toolbar":{"id":"1029"},"width":400,"x_range":{"id":"1006"},"x_scale":{"id":"1010"},"y_range":{"id":"1008"},"y_scale":{"id":"1012"}},"id":"1003","subtype":"Figure","type":"Plot"},{"attributes":{"bottom_units":"screen","coordinates":null,"fill_alpha":0.5,"fill_color":"lightgrey","group":null,"left_units":"screen","level":"overlay","line_alpha":1.0,"line_color":"black","line_dash":[4,4],"line_width":2,"right_units":"screen","syncable":false,"top_units":"screen"},"id":"1028","type":"BoxAnnotation"},{"attributes":{"coordinates":null,"data_source":{"id":"1002"},"glyph":{"id":"1037"},"group":null,"hover_glyph":null,"muted_glyph":{"id":"1039"},"nonselection_glyph":{"id":"1038"},"view":{"id":"1041"}},"id":"1040","type":"GlyphRenderer"},{"attributes":{},"id":"1027","type":"HelpTool"},{"attributes":{},"id":"1010","type":"LinearScale"},{"attributes":{"coordinates":null,"group":null,"text":"Fast Alpha Beta Fractions"},"id":"1004","type":"Title"},{"attributes":{"line_alpha":0.6,"line_color":"blue","line_width":3,"x":{"field":"x"},"y":{"field":"y1"}},"id":"1037","type":"Line"}],"root_ids":["1074"]},"title":"Bokeh Application","version":"2.4.0"}}
</script>
<script type="text/javascript">
    (function() {
    const fn = function() {
        Bokeh.safely(function() {
        (function(root) {
            function embed_document(root) {
            
            const docs_json = document.getElementById('1185').textContent;
            const render_items = [{"docid":"68cf82d8-8d08-40cf-8418-5dcbd5f6848c","root_ids":["1074"],"roots":{"1074":"51d9fa83-acf5-40ed-a2e9-10dad5bd7ad5"}}];
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

Both available models are truncated so that they return a value for $\frac{\beta_{\alpha}}{\beta_{\text{th}}}$ no less than 0.3. This value is then multiplied by the internally calcualted value fo the thermal beta $\beta_{\text{th}}$, using the ion and electron desnity weighted temperatures. The value is then scaled by the fraction of total alpha power vs alpha power produced internally by the plasma to account for additional alpha pressure produced by beam-plasma interactions.

----------------------

## IPDG89 model

This model can be used by setting: `i_beta_fast_alpha` = 0[^1]

### Derivation

Below is the derivation given by Uckan, N. A. et.al. [^2].

$$
\beta_{\alpha} = \frac{2\mu_0 \left(\frac{2n_{\alpha}\langle E_{\alpha}\rangle }{3}\right)}{B^2}
$$

Normalizing to the plasma thermal beta $\beta_{\text{th}} = \beta_{\text{i}}+ \beta_{\text{e}}$

$$
\frac{\beta_{\alpha}}{\beta_{\text{th}}} = \frac{\left(\frac{2E_{\alpha,0}}{3T_{\text{e}}}\right) \left(\frac{n_{\alpha}}{n_{\text{e}}}\right) \left(\frac{\langle E_{\alpha} \rangle}{E_{\alpha,0}}\right)}{1+f_{\text{nT}}}
$$

where $f_{\text{nT}} = f_{\text{n}}f_{\text{T}} = \left(\frac{n_{\text{i}}}{n_{\text{e}}}\right)\left(\frac{T_{\text{i}}}{T_{\text{e}}}\right)$. For $T_{\text{i}} \approx T_{\text{e}}$ and $Z_{\text{eff}} \approx 1.5$ (with $Z = 6$, carbon), $f_{\text{nT}} \approx 0.9$ ad typical local values of fractional fast alpha density , beta and energy are given in the table below.

| $T \ [\mathrm{keV}]$ | $\langle\sigma v\rangle_{{D T}} \ \left[\mathrm{m}^3 / \mathrm{s}\right]$ | $n_{\mathrm{\alpha}} / n_{\text{e}} \ [\%]$ | $\beta_{\mathrm{\alpha}} / \beta_{\mathrm{th}} \ [\%]$ | $\langle E_{\alpha}\rangle/ E_{\alpha,0}$ |
| :---: | :---: | :---: | :---: | :---: |
| $5$ | $1.35 \times 10^{-23}$ | $0.01$ | $0.73$ | $0.3$ |
| $10$ | $1.13 \times 10^{-22}$ | $0.9$ | $4.2$ | $0.34$ |
| $20$ | $4.31 \times 10^{-22}$ | $0.8$ | $19$ | $0.39$ |
| $30$ | $6.65 \times 10^{-22}$ | $1.8$ | $31$ | $0.41$ |
| $40$ | $7.93 \times 10^{-22}$ | $2.7$ | $34$ | $0.41$ |
| $50$ | $8.54 \times 10^{-22}$ | $3.45$ | $34$ | $0.4$ |

Assuming a parabolic profile for temperature and density where $\alpha_{\text{n}} \approx 0.0 - 0.5$ (relatively flat density profile) and $\alpha_{\text{T}} \approx 1$ the above beta ratio becomes:

$$
\gamma_{\alpha} = \frac{\beta_{\alpha}}{\beta_{\text{th}}} \approx\frac{ 0.32 f_{\text{DT}}^2\left(\frac{T_{\text{i}}}{T_{\text{e}}}\right) \langle T_{\text{e,10}} \rangle^{5/2} \langle U_{\alpha \text{e}} \rangle}{1+ f_{\text{nT}}} \\
= 0.32 f_{\text{DT}}^2\left(\frac{T_{\text{i}}}{T_{\text{e}}}\right) \langle T_{\text{e,10}} \rangle^{5/2} \left[\frac{2^{5/2}}{(1+f_{\text{nT}})^{7/2}}\right]
$$

where $\langle U_{\alpha \text{e}} \rangle$ is the fraction of alpha energy given to the electrons, $\langle T \rangle = \langle n_{\text{e}}T_{\text{e}}+n_{\text{i}}T_{\text{i}} / \langle 2n_{\text{e}} \rangle = \langle T_{\text{e}} \rangle \left(1+f_{\text{nT}}\right)/2$ is the density-weighted average  temperature. For analytical simplicity, $\langle \sigma v \rangle_{\text{DT}}$ (the fusion reaction-rate parameter) is approximated as $\langle \sigma v \rangle_{\text{DT}} \approx 1.1 x 10^{-22} \left(T_{\text{i,10}}\right)^2$, which is accurate enough for $T \approx 7-20 \  \mathrm{keV}$. 

For the chosen profiles and $Z_{\text{eff}} \approx 1.5$, the average pressure contribution from fast alphas is $\gamma_{\alpha} \approx 5-20 \%$ for $\langle T \rangle \approx 6-15 \  \mathrm{keV}.$ Direct comparison between the predictions and a large number of 1-1/2-D $\mathtt{WHIST}$ transport code calculations (having similar profile shapes and $Z_{\text{eff}}$ values) shows good agreement (within ±15%) over the temperature range $\left(\langle T \rangle \approx 5-20 \  \mathrm{keV}\right)$ considered. A benchmark between above and $\mathtt{WHIST}$ has resulted in a simple functional fit that is more convenient to use in global analyses:

$$
\gamma_{\alpha} \approx 0.2\left(T_{\text{10}}-0.37\right), \ \text{for} \  Z_{\text{eff}} \approx 1.5, T_{\text{i}}/T_{\text{e}} \approx 1, \langle T \rangle \approx 5-20 \text{keV}
$$

!!! quote "Fit validity"

    "*To zeroth order, the assumption of different profiles $\left(\alpha_{\text{n}} \approx 0-1.0, \alpha_{\text{T}} \approx 0.5-2.0\right)$ did 
    not appear to have any significant effect on this simple fit. As expected, significant 
    deviations from the $\gamma_{\alpha}$ were seen in simulations for anomalous fast alpha diffusion 
    and energy relaxation. In such cases, however, global analysis is not adequate to describe 
    the fast alpha behavior*[^2]"

-------------------

The above derivation form is not that explicitly given in IPDG89[^1] and `PROCESS`. Terms for the electron and DT fuel ion species have been added back in. The $\left(\langle T_{\text{10}}\rangle-0.37\right)$ term still remains.

$$
\beta_{\alpha} =\beta_{\text{th}}\left[0.29 \, \left( \frac{\langle T_{10} \rangle}{20} -
0.37 \right) \, \left( \frac{n_{\text{DT}}}{n_{\text{e}}} \right)^2\right]
$$

For $Z_{\text{eff}} \approx 1.5, T_{\text{i}}/T_{\text{e}} \approx 1, \langle T \rangle \approx 5-20 \text{keV}$


-----------------------

## Ward model

This model can be used by setting: `i_beta_fast_alpha` = 1 (default)[^3]

IPDG89[^1] estimates the pressure in the fast alpha particles. Based on transport code modelling however, looking at the data for local values in IPDG89[^1] this expression overestimates at high temperatures. A better fit to the result at higher temperature is:

$$
\beta_{\alpha} =\beta_{\text{th}} \left[0.26 \, \left( \frac{\langle T_{10} \rangle}{20} -
0.65 \right)^{0.5} \, \left( \frac{n_{\text{DT}}}{n_{\text{e}}} \right)^2\right]
$$

The new fit is much better at higher temperature because it more reliably captures the effect of saturation of fusion cross section with increasing temperature. The new fit cannot be used below peak temperatures of around 15 keV however.


[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989', https://inis.iaea.org/collection/NCLCollectionStore/_Public/21/068/21068960.pdf

[^2]: Uckan, N. A., Tolliver, J. S., Houlberg, W. A., and Attenberger, S. E. Influence of fast alpha diffusion and thermal alpha buildup on tokamak reactor performance. United States: N. p., 1987. Web.https://www.osti.gov/servlets/purl/5611706

[^3]: D.J. Ward, 'PROCESS Fast Alpha Pressure', Work File Note F/PL/PJK/PROCESS/CODE/050