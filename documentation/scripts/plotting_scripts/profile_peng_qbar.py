import numpy as np
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import figure, output_file, save

x = np.linspace(1.0, 3.0, 500)
y1 = 5.0 * 1.3 * (1.0 - (1.0 / x) ** 0.6)
y2 = 5.0 * (1.0 + 2.6 * (1.0 / x) ** 2.8)  # Initial data for the second line

source = ColumnDataSource(data={"x": x, "y1": y1, "y2": y2})

plot = figure(
    x_range=(1.0, 3.0),
    y_range=(0, 10),
    width=400,
    height=400,
    title="PROCESS vs STAR code, average safety factor",
)
plot.xaxis.axis_label = "Aspect ratio"
plot.yaxis.axis_label = "Average safety factor "

plot.line(
    "x",
    "y1",
    source=source,
    line_width=3,
    line_alpha=0.6,
    color="blue",
    legend_label="PROCESS",
)
plot.line(
    "x",
    "y2",
    source=source,
    line_width=3,
    line_alpha=0.6,
    color="red",
    legend_label="STAR",
)

q95 = Slider(start=1.0, end=10, value=5, step=0.1, title=r"95% safety factor, q95")
qbar_0 = Slider(start=1.0, end=10, value=5, step=0.1, title="STAR value")

callback = CustomJS(
    args={"source": source, "q95": q95, "qbar": qbar_0},
    code="""
    const A = q95.value;
    const B = qbar.value;

    const x = source.data['x'];
    const y1 = x.map(xi => A * 1.3 * (1.0 - (1.0 / xi) ** 0.6));
    const y2 = x.map(xi => B * (1.0 + 2.6 * (1.0 / xi) ** 2.8)); // Example transformation for the second line

    source.data['y1'] = y1;
    source.data['y2'] = y2;
    source.change.emit();
""",
)

q95.js_on_change("value", callback)
qbar_0.js_on_change("value", callback)

# Save the plot as HTML
output_file("profile_star_qbar.html", title="Peng q Profile")
save(
    row(plot, column(q95, qbar_0)),
    filename="profile_star_qbar.html",
)
