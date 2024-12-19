import numpy as np
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import figure, output_file, save

x = np.linspace(0, 1, 500)
y1 = 5.0 * (1 - x**2) ** 2.0
y2 = 5.0 / (1 + 2.0 * x**2) ** 2  # Initial data for the second line

source = ColumnDataSource(data=dict(x=x, y1=y1, y2=y2))

plot = figure(
    x_range=(0, 1),
    y_range=(0, 10),
    width=400,
    height=400,
    title="PROCESS vs Hastie Current Profile",
)
plot.xaxis.axis_label = "Normalized Plasma Radius"
plot.yaxis.axis_label = "Current Density, J"

plot.line(
    "x",
    "y1",
    source=source,
    line_width=3,
    line_alpha=0.6,
    color="blue",
    legend_label="PROCESS Parabolic Current Profile",
)
plot.line(
    "x",
    "y2",
    source=source,
    line_width=3,
    line_alpha=0.6,
    color="red",
    legend_label="Hastie Current Profile",
)

n0 = Slider(start=0.1, end=10, value=5, step=0.1, title="Plasma centre value | j0")
alpha = Slider(start=0.01, end=10, value=2, step=0.01, title="Profile Index  | alphaj")

callback = CustomJS(
    args=dict(source=source, n0=n0, alpha=alpha),
    code="""
    const A = n0.value
    const B = alpha.value

    const x = source.data.x
    const y1 = Array.from(x, (x) => A * (1 - x**2)**B)
    const y2 = Array.from(x, (x) => A / (1+B*x**2)**2) // Example transformation for the second line
    source.data = { x, y1, y2 }
""",
)

n0.js_on_change("value", callback)
alpha.js_on_change("value", callback)

# Save the plot as HTML
output_file("profile_hastie.html", title="Hastie Current Profile")
save(row(plot, column(n0, alpha)), filename="profile_hastie.html")
