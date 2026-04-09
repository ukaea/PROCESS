# Divertor

The divertor provides a means of removing plasma reaching the scrape-off layer. 
The principal outputs from the code are the divertor heat load, used to 
determine its lifetime, and its peak temperature. The divertor is cooled either 
by gaseous helium or by pressurised water.

Switch `i_single_null` controls the overall plasma configuration. Setting `i_single_null= 0` 
corresponds to an up-down symmetric, double null configuration, while 
`i_single_null= 1` (the default) assumes a single null plasma with the divertor at the 
bottom of the machine. The vertical build and PF coil current scaling 
algorithms take the value of this switch into account, although not the plasma 
geometry at present.

The Harrison-Kukushkin-Hotston divertor model[^1] developed for ITER is available, but is unlikely to be relevant for a reactor.

The divertor heat flux `pflux_div_heat_load_mw` can be calculated or it can be input by the user. Options are selected using the switch `i_div_heat_load`:

| `i_div_heat_load` | Description |
| :-: | - |
| 0 | the user inputs the value for `pflux_div_heat_load_mw` |
| 1 | the Peng chamber model (`divtart`) is called to calculate `pflux_div_heat_load_mw` |
| 2 | the Wade heat flux model (`divwade`) is called to calculate `pflux_div_heat_load_mw` |

---------------

## Peng Chamber model | `divtart()`

!!! Note ""
    `i_div_heat_load == 1`

The tight aspect ratio tokamak divertor model (`divtart()`) [^5] calculates the divertor heat flux by 
assuming that the power is evenly spread around the divertor chamber by the action of a gaseous 
target. Each divertor is assumed to be approximately triangular in the R,Z plane.

The inner radius of the divertor region is given by:

$$
r_{\text{inner}} = R_0 - \left(a \times \delta \right) - \left(3 \times \Delta r_{\text{FW, plasma}}\right) + \mathtt{drtop}
$$

This is treated as the same as the centrepost and first wall thickness at the divertor height.

The outer radius of the divertor region is simply

$$
r_{\text{outer}} = R_0 + a
$$

The height of the divertor box is taken to be the same as the gap between the plasma X-point and the divertor structure.

Therefore the vertical plate surface area in the divertor box is given by:

$$
A_{\text{div,vertical}} = 2\pi r_{\text{inner}} \times \Delta z_{\text{div}}
$$

The top horizontal plate area is given by:

$$
A_{\text{div,horizontal}} = \pi \left(r_{\text{outer}}^2 - r_{\text{inner}}^2 \right)
$$

We define an angle $\theta$ such that:

$$
\theta = \arctan\left({\frac{\Delta z_{\text{plasma,div}}}{r_{\text{outer}} - r_{\text{inner}}}}\right)
$$

The diagonal plate area is then given by:

$$
A_{\text{div,diagonal}} = \frac{A_{\text{div,horizontal}}}{\cos\left({\theta}\right)^2}
$$

Therefore the total divertor surface area is given by:

$$
A_{\text{div}} = \left(A_{\text{div,vertical}} + A_{\text{div,horizontal}} + A_{\text{div,diagonal}}\right)
$$

If a double null machine is set up with `i_single_null = 0` then the total divertor area is simply:

$$
A_{\text{div}} = 2 \times \left(A_{\text{div,vertical}} + A_{\text{div,horizontal}} + A_{\text{div,diagonal}}\right)
$$

The divertor heat load is then found as:

$$
\mathtt{pflux_div_heat_load_mw} = \frac{\mathtt{p_plasma_separatrix_mw}}{A_{\text{div}}}
$$

!!! warning "Radiated power area"

    The main assumption of the Peng gaseous divertor model is that the power radiated to the divertor is equally radiated in the divertor box across all three surfaces. This may not truly be the case in reality.

The interactive graph below can be used to investigate how changing the key parameters changes the divertor configuration. The grey box represents the first wall, the far right red line represents the right hand edge of the divertor region ($r_{\text{outer}}$), the far left red line represents the left hand edge of the divertor region ($r_{\text{inner}}$) and the blue line represents the bottom of the divertor region ($\Delta z_{\text{plasma,div}}$).

-------------

<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>Bokeh Plot</title>
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
    <script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-widgets-3.6.0.min.js"></script>
    <script type="text/javascript">
        Bokeh.set_log_level("info");
    </script>
  </head>
  <body>
    <div id="d9948402-939a-4dd0-8cec-f597292e2690" data-root-id="p1080" style="display: contents;"></div>
  
    <script type="application/json" id="dd3f126d-9395-4057-9e70-86889cdc9528">
      {"de14171d-5750-472f-9777-8fa7c5b4f122":{"version":"3.6.0","title":"Bokeh Application","roots":[{"type":"object","name":"Column","id":"p1080","attributes":{"children":[{"type":"object","name":"Figure","id":"p1011","attributes":{"width":650,"height":800,"x_range":{"type":"object","name":"Range1d","id":"p1021","attributes":{"end":10}},"y_range":{"type":"object","name":"Range1d","id":"p1022","attributes":{"end":10}},"x_scale":{"type":"object","name":"LinearScale","id":"p1023"},"y_scale":{"type":"object","name":"LinearScale","id":"p1024"},"title":{"type":"object","name":"Title","id":"p1014","attributes":{"text":"divtart divertor box"}},"renderers":[{"type":"object","name":"GlyphRenderer","id":"p1054","attributes":{"data_source":{"type":"object","name":"ColumnDataSource","id":"p1008","attributes":{"selected":{"type":"object","name":"Selection","id":"p1009","attributes":{"indices":[],"line_indices":[]}},"selection_policy":{"type":"object","name":"UnionRenderers","id":"p1010"},"data":{"type":"map","entries":[["x",{"type":"ndarray","array":{"type":"bytes","data":"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"},"shape":[256],"dtype":"float64","order":"little"}],["y",{"type":"ndarray","array":{"type":"bytes","data":"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"},"shape":[256],"dtype":"float64","order":"little"}],["linspace",{"type":"ndarray","array":{"type":"bytes","data":"GC1EVPshCcAgjy7nhO8IwCjxGHoOvQjAMFMDDZiKCMA3te2fIVgIwD8X2DKrJQjAR3nCxTTzB8BP26xYvsAHwFc9l+tHjgfAX5+BftFbB8BmAWwRWykHwG5jVqTk9gbAdsVAN27EBsB+JyvK95EGwIaJFV2BXwbAjuv/7wotBsCWTeqClPoFwJ2v1BUeyAXApRG/qKeVBcCtc6k7MWMFwLXVk866MAXAvTd+YUT+BMDEmWj0zcsEwMz7UodXmQTA1F09GuFmBMDcvyetajQEwOQhEkD0AQTA7IP80n3PA8D05eZlB50DwPxH0fiQagPAA6q7ixo4A8ALDKYepAUDwBNukLEt0wLAG9B6RLegAsAjMmXXQG4CwCqUT2rKOwLAMvY5/VMJAsA6WCSQ3dYBwEK6DiNnpAHAShz5tfBxAcBSfuNIej8BwFrgzdsDDQHAYkK4bo3aAMBqpKIBF6gAwHEGjZSgdQDAeWh3JypDAMCBymG6sxAAwBFZmJp6vP+/IR1twI1X/78x4UHmoPL+v0ClFgy0jf6/UGnrMcco/r9gLcBX2sP9v2/xlH3tXv2/f7VpowD6/L+PeT7JE5X8v549E+8mMPy/rgHoFDrL+7++xbw6TWb7v86JkWBgAfu/3U1mhnOc+r/tETushjf6v/3VD9KZ0vm/DJrk96xt+b8cXrkdwAj5vywijkPTo/i/O+ZiaeY++L9LqjeP+dn3v1tuDLUMdfe/ajLh2h8Q97969rUAM6v2v4q6iiZGRva/mn5fTFnh9b+pQjRybHz1v7kGCZh/F/W/ycrdvZKy9L/YjrLjpU30v+hShwm56PO/+BZcL8yD878H2zBV3x7zvxefBXvyufK/J2PaoAVV8r82J6/GGPDxv0brg+wri/G/Vq9YEj8m8b9mcy04UsHwv3Y3Al5lXPC/CPetB/Hu778of1dTFyXvv0gHAZ89W+6/aI+q6mOR7b+IF1Q2isfsv6if/YGw/eu/xCenzdYz67/kr1AZ/WnqvwQ4+mQjoOm/JMCjsEnW6L9ESE38bwzov2TQ9keWQue/hFigk7x45r+g4Enf4q7lv8Bo8yoJ5eS/4PCcdi8b5L8AeUbCVVHjvyAB8A18h+K/QImZWaK94b9cEUOlyPPgv3yZ7PDuKeC/OEMseSrA3r94U38Qdyzdv7hj0qfDmNu/+HMlPxAF2r8whHjWXHHYv3CUy22p3da/sKQeBfZJ1b/wtHGcQrbTvzDFxDOPItK/cNUXy9uO0L9gy9XEUPbNv9Dre/Ppzsq/UAwiIoOnx7/QLMhQHIDEv1BNbn+1WMG/oNsoXJ1ivL+gHHW5zxO2vwC7gi0Eiq+/AD0b6Gjsor8A/M6KNjuJvwD8zoo2O4k/AD0b6Gjsoj8Au4ItBIqvP4AcdbnPE7Y/oNsoXJ1ivD9QTW5/tVjBP9AsyFAcgMQ/UAwiIoOnxz/Q63vz6c7KP1DL1cRQ9s0/cNUXy9uO0D8wxcQzjyLSP/C0cZxCttM/sKQeBfZJ1T9wlMttqd3WPzCEeNZccdg/8HMlPxAF2j+4Y9Knw5jbP3hTfxB3LN0/OEMseSrA3j98mezw7ingP1wRQ6XI8+A/PImZWaK94T8gAfANfIfiPwB5RsJVUeM/4PCcdi8b5D/AaPMqCeXkP6DgSd/iruU/gFigk7x45j9k0PZHlkLnP0RITfxvDOg/JMCjsEnW6D8EOPpkI6DpP+SvUBn9aeo/xCenzdYz6z+on/2BsP3rP4gXVDaKx+w/aI+q6mOR7T9IBwGfPVvuPyh/V1MXJe8/CPetB/Hu7z90NwJeZVzwP2RzLThSwfA/VK9YEj8m8T9E64PsK4vxPzgnr8YY8PE/KGPaoAVV8j8YnwV78rnyPwjbMFXfHvM/+BZcL8yD8z/oUocJuejzP9iOsuOlTfQ/yMrdvZKy9D+4BgmYfxf1P6hCNHJsfPU/mH5fTFnh9T+IuoomRkb2P3j2tQAzq/Y/bDLh2h8Q9z9cbgy1DHX3P0yqN4/52fc/POZiaeY++D8sIo5D06P4PxxeuR3ACPk/DJrk96xt+T/81Q/SmdL5P+wRO6yGN/o/3E1mhnOc+j/MiZFgYAH7P7zFvDpNZvs/rAHoFDrL+z+gPRPvJjD8P5B5PskTlfw/gLVpowD6/D9w8ZR97V79P2AtwFfaw/0/UGnrMcco/j9ApRYMtI3+PzDhQeag8v4/IB1twI1X/z8QWZiaerz/P4DKYbqzEABAeGh3JypDAEBwBo2UoHUAQGqkogEXqABAYkK4bo3aAEBa4M3bAw0BQFJ+40h6PwFAShz5tfBxAUBCug4jZ6QBQDpYJJDd1gFAMvY5/VMJAkAqlE9qyjsCQCIyZddAbgJAGtB6RLegAkASbpCxLdMCQAwMph6kBQNABKq7ixo4A0D8R9H4kGoDQPTl5mUHnQNA7IP80n3PA0DkIRJA9AEEQNy/J61qNARA1F09GuFmBEDM+1KHV5kEQMSZaPTNywRAvDd+YUT+BEC01ZPOujAFQKxzqTsxYwVAphG/qKeVBUCer9QVHsgFQJZN6oKU+gVAjuv/7wotBkCGiRVdgV8GQH4nK8r3kQZAdsVAN27EBkBuY1ak5PYGQGYBbBFbKQdAXp+BftFbB0BWPZfrR44HQE7brFi+wAdARnnCxTTzB0BAF9gyqyUIQDi17Z8hWAhAMFMDDZiKCEAo8Rh6Dr0IQCCPLueE7whAGC1EVPshCUA="},"shape":[256],"dtype":"float64","order":"little"}]]}}},"view":{"type":"object","name":"CDSView","id":"p1055","attributes":{"filter":{"type":"object","name":"AllIndices","id":"p1056"}}},"glyph":{"type":"object","name":"Patch","id":"p1051","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_width":0,"fill_color":"pink"}},"nonselection_glyph":{"type":"object","name":"Patch","id":"p1052","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.1,"line_width":0,"fill_color":"pink","fill_alpha":0.1,"hatch_alpha":0.1}},"muted_glyph":{"type":"object","name":"Patch","id":"p1053","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.2,"line_width":0,"fill_color":"pink","fill_alpha":0.2,"hatch_alpha":0.2}}}},{"type":"object","name":"GlyphRenderer","id":"p1066","attributes":{"data_source":{"type":"object","name":"ColumnDataSource","id":"p1063","attributes":{"selected":{"type":"object","name":"Selection","id":"p1064","attributes":{"indices":[],"line_indices":[]}},"selection_policy":{"type":"object","name":"UnionRenderers","id":"p1065"},"data":{"type":"map"}}},"view":{"type":"object","name":"CDSView","id":"p1067","attributes":{"filter":{"type":"object","name":"AllIndices","id":"p1068"}}},"glyph":{"type":"object","name":"Quad","id":"p1062","attributes":{"left":{"type":"value","value":0.0},"right":{"type":"value","value":2.7},"bottom":{"type":"value","value":-4.0},"top":{"type":"value","value":10.0},"fill_color":{"type":"value","value":"grey"}}}}},{"type":"object","name":"GlyphRenderer","id":"p1075","attributes":{"data_source":{"type":"object","name":"ColumnDataSource","id":"p1069","attributes":{"selected":{"type":"object","name":"Selection","id":"p1070","attributes":{"indices":[],"line_indices":[]}},"selection_policy":{"type":"object","name":"UnionRenderers","id":"p1071"},"data":{"type":"map","entries":[["x",[3.1,3.1,7.0]],["y",[5.0,6.0,6.0]]]}}},"view":{"type":"object","name":"CDSView","id":"p1076","attributes":{"filter":{"type":"object","name":"AllIndices","id":"p1077"}}},"glyph":{"type":"object","name":"Patch","id":"p1072","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_width":2,"fill_color":"green","fill_alpha":0.6}},"nonselection_glyph":{"type":"object","name":"Patch","id":"p1073","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.1,"line_width":2,"fill_color":"green","fill_alpha":0.1,"hatch_alpha":0.1}},"muted_glyph":{"type":"object","name":"Patch","id":"p1074","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.2,"line_width":2,"fill_color":"green","fill_alpha":0.2,"hatch_alpha":0.2}}}}],"toolbar":{"type":"object","name":"Toolbar","id":"p1020","attributes":{"tools":[{"type":"object","name":"PanTool","id":"p1035"},{"type":"object","name":"WheelZoomTool","id":"p1036","attributes":{"renderers":"auto"}},{"type":"object","name":"BoxZoomTool","id":"p1037","attributes":{"overlay":{"type":"object","name":"BoxAnnotation","id":"p1038","attributes":{"syncable":false,"line_color":"black","line_alpha":1.0,"line_width":2,"line_dash":[4,4],"fill_color":"lightgrey","fill_alpha":0.5,"level":"overlay","visible":false,"left":{"type":"number","value":"nan"},"right":{"type":"number","value":"nan"},"top":{"type":"number","value":"nan"},"bottom":{"type":"number","value":"nan"},"left_units":"canvas","right_units":"canvas","top_units":"canvas","bottom_units":"canvas","handles":{"type":"object","name":"BoxInteractionHandles","id":"p1044","attributes":{"all":{"type":"object","name":"AreaVisuals","id":"p1043","attributes":{"fill_color":"white","hover_fill_color":"lightgray"}}}}}}}},{"type":"object","name":"SaveTool","id":"p1045"},{"type":"object","name":"ResetTool","id":"p1046"},{"type":"object","name":"HelpTool","id":"p1047"}]}},"left":[{"type":"object","name":"LinearAxis","id":"p1030","attributes":{"ticker":{"type":"object","name":"BasicTicker","id":"p1031","attributes":{"mantissas":[1,2,5]}},"formatter":{"type":"object","name":"BasicTickFormatter","id":"p1032"},"axis_label":"Height","major_label_policy":{"type":"object","name":"AllLabels","id":"p1033"}}}],"below":[{"type":"object","name":"LinearAxis","id":"p1025","attributes":{"ticker":{"type":"object","name":"BasicTicker","id":"p1026","attributes":{"mantissas":[1,2,5]}},"formatter":{"type":"object","name":"BasicTickFormatter","id":"p1027"},"axis_label":"Radius","major_label_policy":{"type":"object","name":"AllLabels","id":"p1028"}}}],"center":[{"type":"object","name":"Grid","id":"p1029","attributes":{"axis":{"id":"p1025"}}},{"type":"object","name":"Grid","id":"p1034","attributes":{"dimension":1,"axis":{"id":"p1030"}}},{"type":"object","name":"Legend","id":"p1057","attributes":{"location":"top_left","items":[{"type":"object","name":"LegendItem","id":"p1058","attributes":{"label":{"type":"value","value":"Plasma"},"renderers":[{"id":"p1054"}]}},{"type":"object","name":"LegendItem","id":"p1078","attributes":{"label":{"type":"value","value":"Divertor box"},"renderers":[{"id":"p1075"}]}}]}},{"type":"object","name":"Span","id":"p1059","attributes":{"location":3.1,"dimension":"height","line_color":"red","line_width":2}},{"type":"object","name":"Span","id":"p1060","attributes":{"location":7.0,"dimension":"height","line_color":"red","line_width":2}},{"type":"object","name":"Span","id":"p1061","attributes":{"location":5.0,"line_color":"blue","line_width":2}}]}},{"type":"object","name":"Slider","id":"p1001","attributes":{"js_property_callbacks":{"type":"map","entries":[["change:value",[{"type":"object","name":"CustomJS","id":"p1079","attributes":{"args":{"type":"map","entries":[["source",{"id":"p1008"}],["r0",{"id":"p1001"}],["a",{"type":"object","name":"Slider","id":"p1002","attributes":{"js_property_callbacks":{"type":"map","entries":[["change:value",[{"id":"p1079"}]]]},"title":"Minor radius","start":0.01,"end":2.5,"value":2.0,"step":0.01}}],["delta",{"type":"object","name":"Slider","id":"p1003","attributes":{"js_property_callbacks":{"type":"map","entries":[["change:value",[{"id":"p1079"}]]]},"title":"Triangularity | delta","start":0.0,"end":1,"value":0.5,"step":0.01}}],["kappa",{"type":"object","name":"Slider","id":"p1004","attributes":{"js_property_callbacks":{"type":"map","entries":[["change:value",[{"id":"p1079"}]]]},"title":"Elongation | kappa","start":0.0,"end":2.5,"value":2.0,"step":0.01}}],["square_value",0.0],["vline",{"id":"p1060"}],["vline_inner",{"id":"p1059"}],["hline",{"id":"p1061"}],["first_wall_gap",{"type":"object","name":"Slider","id":"p1005","attributes":{"js_property_callbacks":{"type":"map","entries":[["change:value",[{"id":"p1079"}]]]},"title":"First Wall Gap","start":0.01,"end":0.5,"value":0.3,"step":0.01}}],["plasma_divertor_gap",{"type":"object","name":"Slider","id":"p1006","attributes":{"js_property_callbacks":{"type":"map","entries":[["change:value",[{"id":"p1079"}]]]},"title":"Plasma to Divertor Gap","start":0.01,"end":1.0,"value":1.0,"step":0.01}}],["divertor_height",{"type":"object","name":"Slider","id":"p1007","attributes":{"js_property_callbacks":{"type":"map","entries":[["change:value",[{"id":"p1079"}]]]},"title":"Divertor vertical height","start":0.01,"end":1.5,"value":1.0,"step":0.01}}],["left_rect",{"id":"p1062"}],["divertor_triangle",{"id":"p1075"}]]},"code":"\n    const A = r0.value;\n    const B = a.value;\n    const D = delta.value;\n    const K = kappa.value;\n    const S = square_value;\n    const F = first_wall_gap.value;\n    const P = plasma_divertor_gap.value;\n    const H = divertor_height.value;\n    const x = source.data['linspace'];\n    const R = source.data['x'];\n    const Z = source.data['y'];\n    for (let i = 0; i &lt; x.length; i++) {\n        R[i] = A + B * Math.cos(x[i] + D * Math.sin(x[i]) - S * Math.sin(2 * x[i]));\n        Z[i] = K * B * Math.sin(x[i] + S * Math.sin(2 * x[i]));\n    }\n    vline.location = A + B;\n    vline_inner.location = A - B*D - 3.0*F;\n    hline.location = K*B + P;\n    left_rect.left = 0.0;\n    left_rect.right = A - B -F;\n    divertor_triangle.data_source.data['x'] = [\n        A - B*D - 3.0*F,\n        A - B*D - 3.0*F,\n        A + B\n    ];\n    divertor_triangle.data_source.data['y'] = [\n        K*B + P,\n        K*B + P + H,\n        K*B + P + H\n    ];\n    source.change.emit();\n    divertor_triangle.data_source.change.emit();\n"}}]]]},"title":"Major radius","start":0.1,"end":10,"value":5.0,"step":0.1}},{"id":"p1002"},{"id":"p1003"},{"id":"p1004"},{"id":"p1005"},{"id":"p1006"},{"id":"p1007"}]}}]}}
  </script>
  <script type="text/javascript">
    (function() {
      const fn = function() {
        Bokeh.safely(function() {
          (function(root) {
            function embed_document(root) {
            const docs_json = document.getElementById('dd3f126d-9395-4057-9e70-86889cdc9528').textContent;
            const render_items = [{"docid":"de14171d-5750-472f-9777-8fa7c5b4f122","roots":{"p1080":"d9948402-939a-4dd0-8cec-f597292e2690"},"root_ids":["p1080"]}];
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

--------------------

## Wade Heat Flux Model

!!! Note ""
    `i_div_heat_load == 2`

A divertor heat flux model is provided in Appendix A.II. of [^2].  This uses the Eich scaling 
[^3] and S-factor [^4] to calculate the SOL width at the outboard divertor, mapped to the midplane:

$$
\lambda_{int} = \lambda_{q,Eich} + 1.64S
$$

where

$$
\lambda_{q,Eich} = 1.35 \, P_{\mathrm{SOL}}^{-0.02} \, R_{o}^{0.04} \, B_{p}^{-0.92} \, \epsilon^{0.42}
$$

$$
S = 0.12(n_{e,mid}/10^{19})^{-0.02} \, P_{\mathrm{SOL}}^{-0.21} \, R_{o}^{0.71} \, B_{p}^{-0.82}.
$$

This is then used to calculate the wetted area in the divertor

$$
A_{wetted} = 2\pi N_{div} R \lambda_{int} F_{exp} \sin(\theta_{div})
$$

where $N_{div}$ is the number of divertors (1 or 2), $F_{exp}$ is the relevant flux expansion, and 
$\theta_{div}$ is the tilt of the separatrix relative to the target in the poloidal plane, and has the form

$$
\theta_{div} = \sin^{-1} [(1+1/\alpha_{div}^{2})\sin\beta_{div}],
$$

where

$$
\alpha_{div} = F_{exp}\alpha_{mid}
$$

$$
\alpha_{mid} = \tan^{-1}\frac{B_{p,mid}}{B_{T,mid}}
$$

where $B_{p,mid}$ and $B_{T,mid}$ are the poloidal and toroidal fields on the outer midplane. The 
parameter $\beta_{div}$ is the angle of incidence between the field line and the target.

The divertor heat flux in $\mathrm{MW}/\mathrm{m^{2}}$ is then 

$$
q_{div} = P_{\mathrm{SOL}}(1-f_{rad,div})/A_{wetted}
$$

where $f_{rad,div}$ is the SOL radiative fraction.

For the purposes of this model, the following are inputs:

- Flux expansion $F_{exp}$  (`f_div_flux_expansion`, default = 2)  
- Field line angle with respect to divertor target plate (degrees) $\beta_{div}$ (`deg_div_field_plate`), also 
  available as an iteration variable (170)  
- SOL radiative fraction, $f_{rad,div}$ (`rad_fraction_sol`).

[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)

[^2]: M.R. Wade & J.A. Leuer, 'Cost Drivers for a Tokamak-Based Compact Pilot Plant, Fusion Science and Technology, 77:2, 119-143 (2021)

[^3]: T. Eich et al, 'Scaling of the tokamak near the scrape-off layer H-mode power width and implications for ITER', Nucl. Fusion 53 093031 (2013)

[^4]: A. Scarabosio et al, 'Scaling of the divertor power spreading (S-factor) in open and closed divertor operation in JET and ASDEX Upgrade, Journal of Nuclear Materials, Vol. 463, 49-54 (2015)

[^5]: Y.-K. M. Peng, J. B. Hicksand AEA Fusion, Culham (UK), "Engineering feasibility of tight aspect ratio Tokamak (spherical torus) reactors". 1990. https://inis.iaea.org/records/ey2rf-dah04
