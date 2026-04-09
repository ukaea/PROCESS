import numpy as np
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import figure, output_file, save

T0 = Slider(start=0.1, end=10, value=10.0, step=0.1, title="Plasma centre value | T0")
alpha = Slider(start=0.01, end=10, value=2.0, step=0.01, title="Profile Index  | alphan")
Trho = Slider(start=0.1, end=1, value=0.9, step=0.01, title="Pedestal position | Trho")
Tped = Slider(start=0.01, end=10, value=2.0, step=0.01, title="Pedestal value  | Tped")
Tsep = Slider(start=0.1, end=10, value=0.5, step=0.1, title="Separatrix value | Tsep")
Tbeta = Slider(start=0.01, end=10, value=1.0, step=0.01, title="Beta value  | Tbeta")


output_file("profile_pedestal_plot.html")
x = np.linspace(0, Trho.value, 50000)
y = (
    Tped.value
    + (T0.value - Tped.value)
    * (1 - ((x**Tbeta.value) / (Trho.value**Tbeta.value))) ** alpha.value
)
x2 = np.linspace(Trho.value, 1, 50000)
y2 = Tsep.value + (Tped.value - Tsep.value) * ((1 - x2) / (1 - Trho.value))
source = ColumnDataSource(data={"x": x, "y": y, "x2": x2, "y2": y2})

plot = figure(
    x_range=(0, 1),
    y_range=(0, 10),
    width=400,
    height=400,
    title="Pedestal Profile | H-mode",
)
plot.xaxis.axis_label = r"Normalised Radius, $$ \rho $$"
plot.yaxis.axis_label = r"Density, $$n_e$$"

plot.line("x", "y", source=source, line_width=3, line_alpha=0.6)
plot.line("x2", "y2", source=source, line_width=3, line_alpha=0.6)

callback = CustomJS(
    args={
        "source": source,
        "T0": T0,
        "alpha": alpha,
        "Trho": Trho,
        "Tped": Tped,
        "Tsep": Tsep,
        "Tbeta": Tbeta,
    },
    code="""
   const A = T0.value
    const B = alpha.value
    const C = Trho.value
    const D = Tped.value
    const E = Tsep.value
    const F = Tbeta.value

    const x = Array.from({length: 50000}, (_, i) => i * C/50000 )
    const y = Array.from(x, (x) => D + (A - D)*(1-((x**F)/(C**F)))**B)
    const x2 = Array.from({length: 50000}, (_, i) => C + i * (1 - C)/50000)
    const y2 = Array.from(x2, (x) => E + (D - E)*((1-(x))/(1-C)))

    source.data = { x, y, x2, y2 }
""",
)


T0.js_on_change("value", callback)
alpha.js_on_change("value", callback)
Trho.js_on_change("value", callback)
Tped.js_on_change("value", callback)
Tsep.js_on_change("value", callback)
Tbeta.js_on_change("value", callback)

# Save the plot as HTML
output_file("profile_pedestal_plot.html", title="Pedestal Profile")
save(
    row(plot, column(T0, alpha, Trho, Tped, Tsep, Tbeta)),
    filename="profile_pedestal_plot.html",
)
