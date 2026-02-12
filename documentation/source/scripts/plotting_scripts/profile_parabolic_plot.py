import numpy as np
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import figure, output_file, save

x = np.linspace(0, 1, 500)
y = 5.0 * (1 - x**2) ** 2.0

source = ColumnDataSource(data={"x": x, "y": y})

plot = figure(
    x_range=(0, 1),
    y_range=(0, 10),
    width=400,
    height=400,
    title="Parabolic Profile | L-mode",
)
plot.xaxis.axis_label = r"Normalized Radius, $$ \rho $$"
plot.yaxis.axis_label = r"Density, $$n_e$$"

plot.line("x", "y", source=source, line_width=3, line_alpha=0.6)

n0 = Slider(start=0.1, end=10, value=5, step=0.1, title="Plasma centre value | n0")
alpha = Slider(start=0.01, end=10, value=2, step=0.01, title="Profile Index  | alphan")


callback = CustomJS(
    args={"source": source, "n0": n0, "alpha": alpha},
    code="""
    const A = n0.value
    const B = alpha.value

    const x = source.data.x
    const y = Array.from(x, (x) => A*(1-x**2)**B)
    source.data = { x, y }
""",
)

n0.js_on_change("value", callback)
alpha.js_on_change("value", callback)

# Save the plot as HTML
output_file("profile_parabolic_plot.html", title="Parabolic Profile")
save(
    row(plot, column(n0, alpha)),
    filename="profile_parabolic_plot.html",
)
