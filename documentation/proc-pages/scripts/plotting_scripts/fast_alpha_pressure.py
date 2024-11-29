import numpy as np
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import figure, output_file, save

x = np.linspace(0, 10, 500)
y1 = 0.29 * (x) ** 2(0.75 - 0.37)
y2 = 0.26(0.75 - 0.65) ** 0.5 * (x) ** 2  # Initial data for the second line

source = ColumnDataSource(data=dict(x=x, y1=y1, y2=y2))

plot = figure(
    x_range=(0, 10),
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

te = Slider(start=0.1, end=40, value=15, step=0.1, title="Plasma centre value | j0")

callback = CustomJS(
    args=dict(source=source, te=te),
    code="""
    const A = te.value

    const x = source.data.x
    const y1 = Array.from(x, (x) => 0.29 * (A/20 - 0.37) * (x)**2) // Example transformation for the first line
    const y2 = Array.from(x, (x) => A / (1+B*x**2)**2) // Example transformation for the second line
    source.data = { x, y1, y2 }
""",
)

te.js_on_change("value", callback)

# Save the plot as HTML
output_file("my_bokeh_plot.html", title="My Bokeh Plot")
save(row(plot, column(te)), filename="PROCESS/profile_alpha_beta.html")
