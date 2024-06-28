# Plasma Profiles | `PlasmaProfile`

!!! warning " Un-realistic profiles"

    If `ipedestal >= 1` it is highly recommended to use constraint equation 81 (icc=81). This enforces solutions in which $n_0$ has to be greater than $n_{ped}$. 
    Negative $n_0$ values can also arise during iteration, so it is important to be weary on how low the lower bound for $n_e (\mathtt{dene})$ is set.

In `PROCESS` the density, temperature and current profiles of the plasma for electrons and ions can take two forms depending on the switch value for `ipedestal`. Either without a [pedestal](http://fusionwiki.ciemat.es/wiki/Pedestal) ,`ipedestal = 0` or with a pedestal `ipedestal = 1`.

The files responsible for calculting and storing the profiles are `plasma_profiles.py` and `profiles.py`.

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



    
    
    
        <div class="bk-root" id="2743d8bf-3974-44a5-accc-fa6f529a38a2" data-root-id="1046"></div>
    
    



<script type="application/json" id="1224">
    {"f8ad2f4a-b793-4e84-8627-1dfb70865666":{"defs":[],"roots":{"references":[{"attributes":{},"id":"1054","type":"BasicTickFormatter"},{"attributes":{},"id":"1025","type":"WheelZoomTool"},{"attributes":{},"id":"1024","type":"PanTool"},{"attributes":{},"id":"1052","type":"AllLabels"},{"attributes":{"overlay":{"id":"1030"}},"id":"1026","type":"BoxZoomTool"},{"attributes":{"end":10,"js_property_callbacks":{"change:value":[{"id":"1044"}]},"start":0.01,"step":0.01,"title":"\u03b1\u2099: Density profile exponent","value":1},"id":"1005","type":"Slider"},{"attributes":{},"id":"1027","type":"SaveTool"},{"attributes":{},"id":"1028","type":"ResetTool"},{"attributes":{"coordinates":null,"group":null,"text":"Parabolic Plasma Density Profile"},"id":"1047","type":"Title"},{"attributes":{},"id":"1055","type":"AllLabels"},{"attributes":{"line_alpha":0.6,"line_color":"#1f77b4","line_width":3,"x":{"field":"x"},"y":{"field":"y"}},"id":"1039","type":"Line"},{"attributes":{"data":{"x":{"__ndarray__":"AAAAAAAAAAD7EnmctWpgP/sSeZy1anA/eJy1ahCgeD/7EnmctWqAP7pXlwNjhYQ/eJy1ahCgiD834dPRvbqMP/sSeZy1apA/WjUIUAx4kj+6V5cDY4WUPxl6Jre5kpY/eJy1ahCgmD/YvkQeZ62aPzfh09G9upw/lwNjhRTInj/7EnmctWqgPyukQPZgcaE/WjUIUAx4oj+Kxs+pt36jP7pXlwNjhaQ/6eheXQ6MpT8Zeia3uZKmP0kL7hBlmac/eJy1ahCgqD+oLX3Eu6apP9i+RB5nrao/CFAMeBK0qz834dPRvbqsP2dymytpwa0/lwNjhRTIrj/GlCrfv86vP/sSeZy1arA/k9tcSQvusD8rpED2YHGxP8NsJKO29LE/WjUIUAx4sj/y/ev8YfuyP4rGz6m3frM/Io+zVg0CtD+6V5cDY4W0P1Ige7C4CLU/6eheXQ6MtT+BsUIKZA+2Pxl6Jre5krY/sUIKZA8Wtz9JC+4QZZm3P+HT0b26HLg/eJy1ahCguD8QZZkXZiO5P6gtfcS7prk/QPZgcREquj/YvkQeZ626P3CHKMu8MLs/CFAMeBK0uz+fGPAkaDe8Pzfh09G9urw/z6m3fhM+vT9ncpsracG9P/86f9i+RL4/lwNjhRTIvj8uzEYyaku/P8aUKt+/zr8/ry4HxgopwD/7EnmctWrAP0f36nJgrMA/k9tcSQvuwD/fv84fti/BPyukQPZgccE/d4iyzAuzwT/DbCSjtvTBPw5RlnlhNsI/WjUIUAx4wj+mGXomt7nCP/L96/xh+8I/PuJd0ww9wz+Kxs+pt37DP9aqQYBiwMM/Io+zVg0CxD9ucyUtuEPEP7pXlwNjhcQ/BjwJ2g3HxD9SIHuwuAjFP54E7YZjSsU/6eheXQ6MxT81zdAzuc3FP4GxQgpkD8Y/zZW04A5Rxj8Zeia3uZLGP2VemI1k1MY/sUIKZA8Wxz/9Jnw6ulfHP0kL7hBlmcc/le9f5w/bxz/h09G9uhzIPy24Q5RlXsg/eJy1ahCgyD/EgCdBu+HIPxBlmRdmI8k/XEkL7hBlyT+oLX3Eu6bJP/QR75pm6Mk/QPZgcREqyj+M2tJHvGvKP9i+RB5nrco/JKO29BHvyj9whyjLvDDLP7xrmqFncss/CFAMeBK0yz9TNH5OvfXLP58Y8CRoN8w/6/xh+xJ5zD834dPRvbrMP4PFRaho/Mw/z6m3fhM+zT8bjilVvn/NP2dymytpwc0/s1YNAhQDzj//On/YvkTOP0sf8a5phs4/lwNjhRTIzj/j59RbvwnPPy7MRjJqS88/erC4CBWNzz/GlCrfv87PP4k8zlo1CNA/ry4Hxgop0D/VIEAx4EnQP/sSeZy1atA/IQWyB4uL0D9H9+pyYKzQP23pI941zdA/k9tcSQvu0D+5zZW04A7RP9+/zh+2L9E/BbIHi4tQ0T8rpED2YHHRP1GWeWE2ktE/d4iyzAuz0T+deus34dPRP8NsJKO29NE/6F5dDowV0j8OUZZ5YTbSPzRDz+Q2V9I/WjUIUAx40j+AJ0G74ZjSP6YZeia3udI/zAuzkYza0j/y/ev8YfvSPxjwJGg3HNM/PuJd0ww90z9k1JY+4l3TP4rGz6m3ftM/sLgIFY2f0z/WqkGAYsDTP/yceus34dM/Io+zVg0C1D9IgezB4iLUP25zJS24Q9Q/lGVemI1k1D+6V5cDY4XUP+BJ0G44ptQ/BjwJ2g3H1D8sLkJF4+fUP1Ige7C4CNU/eBK0G44p1T+eBO2GY0rVP8P2JfI4a9U/6eheXQ6M1T8P25fI46zVPzXN0DO5zdU/W78Jn47u1T+BsUIKZA/WP6eje3U5MNY/zZW04A5R1j/zh+1L5HHWPxl6Jre5ktY/P2xfIo+z1j9lXpiNZNTWP4tQ0fg59dY/sUIKZA8W1z/XNEPP5DbXP/0mfDq6V9c/Ixm1pY941z9JC+4QZZnXP2/9Jnw6utc/le9f5w/b1z+74ZhS5fvXP+HT0b26HNg/B8YKKZA92D8tuEOUZV7YP1OqfP86f9g/eJy1ahCg2D+eju7V5cDYP8SAJ0G74dg/6nJgrJAC2T8QZZkXZiPZPzZX0oI7RNk/XEkL7hBl2T+CO0RZ5oXZP6gtfcS7ptk/zh+2L5HH2T/0Ee+aZujZPxoEKAY8Cdo/QPZgcREq2j9m6Jnc5kraP4za0ke8a9o/sswLs5GM2j/YvkQeZ63aP/6wfYk8zto/JKO29BHv2j9Kle9f5w/bP3CHKMu8MNs/lnlhNpJR2z+8a5qhZ3LbP+Jd0ww9k9s/CFAMeBK02z8uQkXj59TbP1M0fk699ds/eSa3uZIW3D+fGPAkaDfcP8UKKZA9WNw/6/xh+xJ53D8R75pm6JncPzfh09G9utw/XdMMPZPb3D+DxUWoaPzcP6m3fhM+Hd0/z6m3fhM+3T/1m/Dp6F7dPxuOKVW+f90/QYBiwJOg3T9ncpsracHdP41k1JY+4t0/s1YNAhQD3j/ZSEZt6SPeP/86f9i+RN4/JS24Q5Rl3j9LH/GuaYbeP3ERKho/p94/lwNjhRTI3j+99Zvw6ejeP+Pn1Fu/Cd8/CdoNx5Qq3z8uzEYyakvfP1S+f50/bN8/erC4CBWN3z+govFz6q3fP8aUKt+/zt8/7IZjSpXv3z+JPM5aNQjgP5y1ahCgGOA/ry4Hxgop4D/Cp6N7dTngP9UgQDHgSeA/6Jnc5kpa4D/7EnmctWrgPw6MFVIge+A/IQWyB4uL4D80fk699ZvgP0f36nJgrOA/WnCHKMu84D9t6SPeNc3gP4BiwJOg3eA/k9tcSQvu4D+mVPn+df7gP7nNlbTgDuE/zEYyaksf4T/fv84fti/hP/I4a9UgQOE/BbIHi4tQ4T8YK6RA9mDhPyukQPZgceE/Ph3dq8uB4T9RlnlhNpLhP2QPFhehouE/d4iyzAuz4T+KAU+CdsPhP5166zfh0+E/sPOH7Uvk4T/DbCSjtvThP9blwFghBeI/6F5dDowV4j/71/nD9iXiPw5RlnlhNuI/IcoyL8xG4j80Q8/kNlfiP0e8a5qhZ+I/WjUIUAx44j9trqQFd4jiP4AnQbvhmOI/k6DdcEyp4j+mGXomt7niP7mSFtwhyuI/zAuzkYza4j/fhE9H9+riP/L96/xh++I/BXeIsswL4z8Y8CRoNxzjPytpwR2iLOM/PuJd0ww94z9RW/qId03jP2TUlj7iXeM/d00z9Exu4z+Kxs+pt37jP50/bF8ij+M/sLgIFY2f4z/DMaXK96/jP9aqQYBiwOM/6SPeNc3Q4z/8nHrrN+HjPw8WF6Gi8eM/Io+zVg0C5D81CFAMeBLkP0iB7MHiIuQ/W/qId00z5D9ucyUtuEPkP4HsweIiVOQ/lGVemI1k5D+n3vpN+HTkP7pXlwNjheQ/zdAzuc2V5D/gSdBuOKbkP/PCbCSjtuQ/BjwJ2g3H5D8ZtaWPeNfkPywuQkXj5+Q/P6fe+k345D9SIHuwuAjlP2WZF2YjGeU/eBK0G44p5T+Li1DR+DnlP54E7YZjSuU/sH2JPM5a5T/D9iXyOGvlP9Zvwqeje+U/6eheXQ6M5T/8YfsSeZzlPw/bl8jjrOU/IlQ0fk695T81zdAzuc3lP0hGbekj3uU/W78Jn47u5T9uOKZU+f7lP4GxQgpkD+Y/lCrfv84f5j+no3t1OTDmP7ocGCukQOY/zZW04A5R5j/gDlGWeWHmP/OH7UvkceY/BgGKAU+C5j8Zeia3uZLmPyzzwmwko+Y/P2xfIo+z5j9S5fvX+cPmP2VemI1k1OY/eNc0Q8/k5j+LUNH4OfXmP57Jba6kBec/sUIKZA8W5z/Eu6YZeibnP9c0Q8/kNuc/6q3fhE9H5z/9Jnw6ulfnPxCgGPAkaOc/Ixm1pY945z82klFb+ojnP0kL7hBlmec/XISKxs+p5z9v/SZ8OrrnP4J2wzGlyuc/le9f5w/b5z+oaPyceuvnP7vhmFLl++c/zlo1CFAM6D/h09G9uhzoP/RMbnMlLeg/B8YKKZA96D8aP6fe+k3oPy24Q5RlXug/QDHgSdBu6D9Tqnz/On/oP2YjGbWlj+g/eJy1ahCg6D+LFVIge7DoP56O7tXlwOg/sQeLi1DR6D/EgCdBu+HoP9f5w/Yl8ug/6nJgrJAC6T/96/xh+xLpPxBlmRdmI+k/I941zdAz6T82V9KCO0TpP0nQbjimVOk/XEkL7hBl6T9vwqeje3XpP4I7RFnmhek/lbTgDlGW6T+oLX3Eu6bpP7umGXomt+k/zh+2L5HH6T/hmFLl+9fpP/QR75pm6Ok/B4uLUNH46T8aBCgGPAnqPy19xLumGeo/QPZgcREq6j9Tb/0mfDrqP2bomdzmSuo/eWE2klFb6j+M2tJHvGvqP59Tb/0mfOo/sswLs5GM6j/FRaho/JzqP9i+RB5nreo/6zfh09G96j/+sH2JPM7qPxEqGj+n3uo/JKO29BHv6j83HFOqfP/qP0qV71/nD+s/XQ6MFVIg6z9whyjLvDDrP4MAxYAnQes/lnlhNpJR6z+p8v3r/GHrP7xrmqFncus/z+Q2V9KC6z/iXdMMPZPrP/XWb8Kno+s/CFAMeBK06z8byagtfcTrPy5CRePn1Os/QbvhmFLl6z9TNH5OvfXrP2atGgQoBuw/eSa3uZIW7D+Mn1Nv/SbsP58Y8CRoN+w/spGM2tJH7D/FCimQPVjsP9iDxUWoaOw/6/xh+xJ57D/+df6wfYnsPxHvmmbomew/JGg3HFOq7D834dPRvbrsP0pacIcoy+w/XdMMPZPb7D9wTKny/evsP4PFRaho/Ow/lj7iXdMM7T+pt34TPh3tP7wwG8moLe0/z6m3fhM+7T/iIlQ0fk7tP/Wb8OnoXu0/CBWNn1Nv7T8bjilVvn/tPy4HxgopkO0/QYBiwJOg7T9U+f51/rDtP2dymytpwe0/eus34dPR7T+NZNSWPuLtP6DdcEyp8u0/s1YNAhQD7j/Gz6m3fhPuP9lIRm3pI+4/7MHiIlQ07j//On/YvkTuPxK0G44pVe4/JS24Q5Rl7j84plT5/nXuP0sf8a5phu4/XpiNZNSW7j9xESoaP6fuP4SKxs+pt+4/lwNjhRTI7j+qfP86f9juP731m/Dp6O4/0G44plT57j/j59RbvwnvP/ZgcREqGu8/CdoNx5Qq7z8bU6p8/zrvPy7MRjJqS+8/QUXj59Rb7z9Uvn+dP2zvP2c3HFOqfO8/erC4CBWN7z+NKVW+f53vP6Ci8XPqre8/sxuOKVW+7z/GlCrfv87vP9kNx5Qq3+8/7IZjSpXv7z8AAAAAAADwPw==","dtype":"float64","order":"little","shape":[500]},"y":{"__ndarray__":"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","dtype":"float64","order":"little","shape":[500]}},"selected":{"id":"1057"},"selection_policy":{"id":"1056"}},"id":"1006","type":"ColumnDataSource"},{"attributes":{"line_alpha":0.1,"line_color":"#1f77b4","line_width":3,"x":{"field":"x"},"y":{"field":"y"}},"id":"1040","type":"Line"},{"attributes":{"bottom_units":"screen","coordinates":null,"fill_alpha":0.5,"fill_color":"lightgrey","group":null,"left_units":"screen","level":"overlay","line_alpha":1.0,"line_color":"black","line_dash":[4,4],"line_width":2,"right_units":"screen","syncable":false,"top_units":"screen"},"id":"1030","type":"BoxAnnotation"},{"attributes":{"line_alpha":0.2,"line_color":"#1f77b4","line_width":3,"x":{"field":"x"},"y":{"field":"y"}},"id":"1041","type":"Line"},{"attributes":{"coordinates":null,"data_source":{"id":"1006"},"glyph":{"id":"1039"},"group":null,"hover_glyph":null,"muted_glyph":{"id":"1041"},"nonselection_glyph":{"id":"1040"},"view":{"id":"1043"}},"id":"1042","type":"GlyphRenderer"},{"attributes":{"end":10,"js_property_callbacks":{"change:value":[{"id":"1044"}]},"start":0.01,"step":0.01,"title":"n\u2080: On-axis density","value":5},"id":"1004","type":"Slider"},{"attributes":{"args":{"alphan":{"id":"1005"},"n0":{"id":"1004"},"source":{"id":"1006"}},"code":"\n    const data = source.data;\n    const x = data['x']\n    const y = data['y']\n    for (var i = 0; i &lt; x.length; i++) {\n        y[i] = n0.value * (1 - x[i]**2) ** alphan.value;\n    }\n    source.change.emit();\n"},"id":"1044","type":"CustomJS"},{"attributes":{},"id":"1051","type":"BasicTickFormatter"},{"attributes":{},"id":"1056","type":"UnionRenderers"},{"attributes":{"below":[{"id":"1016"}],"center":[{"id":"1019"},{"id":"1023"}],"height":400,"left":[{"id":"1020"}],"renderers":[{"id":"1042"}],"sizing_mode":"stretch_width","title":{"id":"1047"},"toolbar":{"id":"1031"},"width":400,"x_range":{"id":"1008"},"x_scale":{"id":"1012"},"y_range":{"id":"1010"},"y_scale":{"id":"1014"}},"id":"1007","subtype":"Figure","type":"Plot"},{"attributes":{},"id":"1057","type":"Selection"},{"attributes":{"end":10},"id":"1010","type":"Range1d"},{"attributes":{"source":{"id":"1006"}},"id":"1043","type":"CDSView"},{"attributes":{"tools":[{"id":"1024"},{"id":"1025"},{"id":"1026"},{"id":"1027"},{"id":"1028"},{"id":"1029"}]},"id":"1031","type":"Toolbar"},{"attributes":{},"id":"1017","type":"BasicTicker"},{"attributes":{},"id":"1008","type":"Range1d"},{"attributes":{},"id":"1012","type":"LinearScale"},{"attributes":{},"id":"1014","type":"LinearScale"},{"attributes":{"axis_label":"$$\\rho/a$$","coordinates":null,"formatter":{"id":"1054"},"group":null,"major_label_policy":{"id":"1055"},"ticker":{"id":"1017"}},"id":"1016","type":"LinearAxis"},{"attributes":{"axis_label":"$$n_{e}$$","coordinates":null,"formatter":{"id":"1051"},"group":null,"major_label_policy":{"id":"1052"},"ticker":{"id":"1021"}},"id":"1020","type":"LinearAxis"},{"attributes":{},"id":"1029","type":"HelpTool"},{"attributes":{"axis":{"id":"1016"},"coordinates":null,"group":null,"ticker":null},"id":"1019","type":"Grid"},{"attributes":{"axis":{"id":"1020"},"coordinates":null,"dimension":1,"group":null,"ticker":null},"id":"1023","type":"Grid"},{"attributes":{},"id":"1021","type":"BasicTicker"},{"attributes":{"children":[{"id":"1007"},{"id":"1045"}],"sizing_mode":"stretch_width"},"id":"1046","type":"Column"},{"attributes":{"children":[{"id":"1004"},{"id":"1005"}],"sizing_mode":"stretch_width"},"id":"1045","type":"Column"}],"root_ids":["1046"]},"title":"Bokeh Application","version":"2.4.0"}}
</script>
<script type="text/javascript">
    (function() {
    const fn = function() {
        Bokeh.safely(function() {
        (function(root) {
            function embed_document(root) {
            
            const docs_json = document.getElementById('1224').textContent;
            const render_items = [{"docid":"f8ad2f4a-b793-4e84-8627-1dfb70865666","root_ids":["1046"],"roots":{"1046":"2743d8bf-3974-44a5-accc-fa6f529a38a2"}}];
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






<div class="bk-root" id="7c8d29d7-49da-474e-b90a-27f627d3c601" data-root-id="1053"></div>





<script type="application/json" id="1164">
{"8410e422-b8ca-495a-9647-24a4cbe0a3bb":{"defs":[],"roots":{"references":[{"attributes":{"end":10,"format":"0.00","js_property_callbacks":{"change:value":[{"id":"1050"}]},"start":0.01,"step":0.01,"title":"\ud835\udc47\u209b\u2091\u209a: Seperatrix temperature","value":0.5},"id":"1047","type":"Slider"},{"attributes":{},"id":"1027","type":"HelpTool"},{"attributes":{"end":1,"start":0},"id":"1006","type":"DataRange1d"},{"attributes":{"end":25,"format":"0.00","js_property_callbacks":{"change:value":[{"id":"1050"}]},"start":0.01,"step":0.01,"title":"\ud835\udc47\u2080: Central plasma temperature","value":15},"id":"1046","type":"Slider"},{"attributes":{},"id":"1019","type":"BasicTicker"},{"attributes":{},"id":"1015","type":"BasicTicker"},{"attributes":{"axis":{"id":"1018"},"coordinates":null,"dimension":1,"group":null,"ticker":null},"id":"1021","type":"Grid"},{"attributes":{"axis_label":"$$T_{e}$$","coordinates":null,"formatter":{"id":"1056"},"group":null,"major_label_policy":{"id":"1057"},"ticker":{"id":"1019"}},"id":"1018","type":"LinearAxis"},{"attributes":{"below":[{"id":"1014"}],"center":[{"id":"1017"},{"id":"1021"}],"left":[{"id":"1018"}],"renderers":[{"id":"1040"}],"sizing_mode":"stretch_width","title":{"id":"1004"},"toolbar":{"id":"1029"},"toolbar_location":"below","x_range":{"id":"1006"},"x_scale":{"id":"1010"},"y_range":{"id":"1008"},"y_scale":{"id":"1012"}},"id":"1003","subtype":"Figure","type":"Plot"},{"attributes":{"start":0},"id":"1008","type":"DataRange1d"},{"attributes":{"end":10,"format":"0.00","js_property_callbacks":{"change:value":[{"id":"1050"}]},"start":0.01,"step":0.01,"title":"\ud835\udefc\u209c","value":2.0},"id":"1044","type":"Slider"},{"attributes":{},"id":"1023","type":"WheelZoomTool"},{"attributes":{"coordinates":null,"group":null,"text":"Pedestal Plasma Temeprature Profile"},"id":"1004","type":"Title"},{"attributes":{},"id":"1059","type":"BasicTickFormatter"},{"attributes":{},"id":"1012","type":"LinearScale"},{"attributes":{"data":{"x":{"__ndarray__":"AAAAAAAAAABbv1Kg1q+EP1u/UqDWr5Q/CB988MEHnz9bv1Kg1q+kPzJvZ0jM26k/CB988MEHrz9wZ0jM2xmyP1u/UqDWr7Q/RhdddNFFtz8yb2dIzNu5Px3HcRzHcbw/CB988MEHvz96O0Ni3s7AP3BnSMzbGcI/ZZNNNtlkwz9bv1Kg1q/EP1HrVwrU+sU/RhdddNFFxz88Q2LezpDIPzJvZ0jM28k/J5tssskmyz8dx3Ecx3HMPxPzdobEvM0/CB988MEHzz9/pUCtXynQP3o7Q2LeztA/ddFFF1100T9wZ0jM2xnSP2r9SoFav9I/ZZNNNtlk0z9gKVDrVwrUP1u/UqDWr9Q/VlVVVVVV1T9R61cK1PrVP0yBWr9SoNY/RhdddNFF1z9BrV8pUOvXPzxDYt7OkNg/N9lkk0022T8yb2dIzNvZPy0Fav1Kgdo/J5tssskm2z8iMW9nSMzbPx3HcRzHcdw/GF100UUX3T8T83aGxLzdPw6JeTtDYt4/CB988MEH3z8DtX6lQK3fP3+lQK1fKeA/ffDBBx984D96O0Ni3s7gP3eGxLydIeE/ddFFF1104T9yHMdxHMfhP3BnSMzbGeI/bbLJJpts4j9q/UqBWr/iP2hIzNsZEuM/ZZNNNtlk4z9j3s6QmLfjP2ApUOtXCuQ/XnTRRRdd5D9bv1Kg1q/kP1gK1PqVAuU/VlVVVVVV5T9ToNavFKjlP1HrVwrU+uU/TjbZZJNN5j9MgVq/UqDmP0nM2xkS8+Y/RhdddNFF5z9EYt7OkJjnP0GtXylQ6+c/P/jggw8+6D88Q2LezpDoPzmO4ziO4+g/N9lkk0026T80JObtDInpPzJvZ0jM2+k/L7rooosu6j8tBWr9SoHqPypQ61cK1Oo/J5tssskm6z8l5u0MiXnrPyIxb2dIzOs/IHzwwQcf7D8dx3Ecx3HsPxsS83aGxOw/GF100UUX7T8VqPUrBWrtPxPzdobEvO0/ED744IMP7j8OiXk7Q2LuPwvU+pUCte4/CB988MEH7z8Gav1KgVrvPwO1fqVAre8/AAAAAAAA8D8=","dtype":"float64","order":"little","shape":[100]},"y":{"__ndarray__":"AAAAAAAALkCzuzXHc/4tQHhfgWnP+S1ATD3izBPyLUB9iKxwQuctQK1ViW1d2S1A0pp2dWfILUA1L8fTY7QtQHHLIm1WnS1AdQmGv0ODLUCAZELiMGYtQCk5/oUjRi1AVsW09CEjLUBAKLYRM/0sQHdip1le1CxA2lWC4quoLECcxZVbJHosQEJWhQ3RSCxApo1J2rsULED20i89790rQK9u2kp2pCtApIpAsVxoK0D6Ma63rikrQChRxD556CpA+7V4wMmkKkCODxZQrl4qQFbuO5o1FipAFMTe5G7LKUDf40cPan4pQCKCFZI3LylAmbQ6f+jdKEBUcv+BjoooQLiTAN87NShAedIvdAPeJ0CfydO4+IQnQIf1h70vKidA4LM8LL3NJkCqQzdItm8mQDzFEe4wECZAPDq7k0OvJUCmhXdIBU0lQMZr37SN6SRAQJLgGvWEJEADgL1VVB8kQFmdDdrEuCNA2jO9tWBRI0Bybg2QQukiQGRZlKmFgCJAP+I83EUXIkDp10abn60hQJzqRvOvQyFA4asmipTZIECZjiSfa28gQPbm0wpUBSBA79Q5fto2H0D3X3mOrmMeQEpfio9lkR1A8IU6/UHAHECESQGGh/AbQEzi/wp7IhtANksBoGJWGkDSQXqLhYwZQF5GiUYsxRhAsZv2fKAAGEBWRzQNLT8XQHQRXggegRZA2YQ5ssDGFUAB7zWBYxAVQANgbB5WXhRApKqfZemwE0BKZDxlbwgTQAjlWF47ZRJAjEe1xKHHEUAyabs++C8RQPzpfqWVnhBAjCy9BNITEEBmrLo1DSAPQLyd4LEbJw5AToJhxYY9DUCYNAYdCGQMQFwU68pbmwtAtgaARkDkCkD+dYhsdj8KQORRG3/BrQlAXg+jJecvCUCrqN1sr8YIQF+d3MbkcghAUfIEC1Q1CECqMQ92zA4IQNtqB6ofAAhAK7rooosuBkB+pUCtXykEQMqQmLczJAJAHXzwwQcfAEDSzpCYtzP8P3ilQK1fKfg/HXzwwQcf9D+2UqDWrxTwP7ZSoNavFOg/AAAAAAAA4D8=","dtype":"float64","order":"little","shape":[100]}},"selected":{"id":"1062"},"selection_policy":{"id":"1061"}},"id":"1002","type":"ColumnDataSource"},{"attributes":{},"id":"1061","type":"UnionRenderers"},{"attributes":{"children":[{"id":"1052"}],"sizing_mode":"stretch_width"},"id":"1053","type":"Column"},{"attributes":{},"id":"1062","type":"Selection"},{"attributes":{"end":1,"format":"0.00","js_property_callbacks":{"change:value":[{"id":"1050"}]},"start":0.01,"step":0.01,"title":"\ud835\udf0c\u209a\u2091d,T: Pedestal position","value":0.9},"id":"1042","type":"Slider"},{"attributes":{},"id":"1060","type":"AllLabels"},{"attributes":{"source":{"id":"1002"}},"id":"1041","type":"CDSView"},{"attributes":{"axis_label":"$$\\rho/a$$","coordinates":null,"formatter":{"id":"1059"},"group":null,"major_label_policy":{"id":"1060"},"ticker":{"id":"1015"}},"id":"1014","type":"LinearAxis"},{"attributes":{"tools":[{"id":"1022"},{"id":"1023"},{"id":"1024"},{"id":"1025"},{"id":"1026"},{"id":"1027"}]},"id":"1029","type":"Toolbar"},{"attributes":{"line_alpha":0.2,"line_color":"#1f77b4","line_width":2,"x":{"field":"x"},"y":{"field":"y"}},"id":"1039","type":"Line"},{"attributes":{"line_color":"#1f77b4","line_width":2,"x":{"field":"x"},"y":{"field":"y"}},"id":"1037","type":"Line"},{"attributes":{},"id":"1056","type":"BasicTickFormatter"},{"attributes":{"children":[{"id":"1042"},{"id":"1043"},{"id":"1044"},{"id":"1045"},{"id":"1046"},{"id":"1047"}]},"id":"1051","type":"Column"},{"attributes":{"coordinates":null,"data_source":{"id":"1002"},"glyph":{"id":"1037"},"group":null,"hover_glyph":null,"muted_glyph":{"id":"1039"},"nonselection_glyph":{"id":"1038"},"view":{"id":"1041"}},"id":"1040","type":"GlyphRenderer"},{"attributes":{},"id":"1057","type":"AllLabels"},{"attributes":{"end":10,"format":"0.00","js_property_callbacks":{"change:value":[{"id":"1050"}]},"start":0.01,"step":0.01,"title":"\ud835\udefd\u209c","value":2},"id":"1043","type":"Slider"},{"attributes":{"bottom_units":"screen","coordinates":null,"fill_alpha":0.5,"fill_color":"lightgrey","group":null,"left_units":"screen","level":"overlay","line_alpha":1.0,"line_color":"black","line_dash":[4,4],"line_width":2,"right_units":"screen","syncable":false,"top_units":"screen"},"id":"1028","type":"BoxAnnotation"},{"attributes":{"children":[{"id":"1003"},{"id":"1051"}],"sizing_mode":"stretch_width"},"id":"1052","type":"Row"},{"attributes":{"end":10,"format":"0.00","js_property_callbacks":{"change:value":[{"id":"1050"}]},"start":0.01,"step":0.01,"title":"\ud835\udc47\u209a\u2091d: Pedestal temperature","value":3},"id":"1045","type":"Slider"},{"attributes":{"line_alpha":0.1,"line_color":"#1f77b4","line_width":2,"x":{"field":"x"},"y":{"field":"y"}},"id":"1038","type":"Line"},{"attributes":{},"id":"1026","type":"ResetTool"},{"attributes":{"args":{"alphat":{"id":"1044"},"rhopedt":{"id":"1042"},"source":{"id":"1002"},"t0":{"id":"1046"},"tbeta":{"id":"1043"},"teped":{"id":"1045"},"tesep":{"id":"1047"}},"code":"\n    const data = source.data;\n    const x = data['x'];\n    const y = data['y'];\n    const rho_index = x.map((value) =&gt; value &lt;= rhopedt.value);\n    \n    for (let i = 0; i &lt; x.length; i++) {\n        if (rho_index[i]) {\n            y[i] = teped.value + (t0.value - teped.value) * Math.pow(1 - (x[i] / rhopedt.value) ** tbeta.value, alphat.value);\n        } else {\n            y[i] = tesep.value + (teped.value - tesep.value) * (1 - x[i]) / (1 - rhopedt.value);\n        }\n    }\n    \n    source.change.emit();\n"},"id":"1050","type":"CustomJS"},{"attributes":{},"id":"1025","type":"SaveTool"},{"attributes":{"axis":{"id":"1014"},"coordinates":null,"group":null,"ticker":null},"id":"1017","type":"Grid"},{"attributes":{"overlay":{"id":"1028"}},"id":"1024","type":"BoxZoomTool"},{"attributes":{},"id":"1010","type":"LinearScale"},{"attributes":{},"id":"1022","type":"PanTool"}],"root_ids":["1053"]},"title":"Bokeh Application","version":"2.4.0"}}
</script>
<script type="text/javascript">
(function() {
const fn = function() {
Bokeh.safely(function() {
(function(root) {
function embed_document(root) {

const docs_json = document.getElementById('1164').textContent;
const render_items = [{"docid":"8410e422-b8ca-495a-9647-24a4cbe0a3bb","root_ids":["1053"],"roots":{"1053":"7c8d29d7-49da-474e-b90a-27f627d3c601"}}];
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
















<figure markdown>
![UML of profiles](../images/profiles_uml.png){ width="100%"}
<figcaption>Figure 1: UML class breakdown of the plasma profiles</figcaption>
</figure>


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