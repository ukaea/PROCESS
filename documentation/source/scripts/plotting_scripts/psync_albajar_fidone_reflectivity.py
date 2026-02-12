import numpy as np
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, output_file, save

x = np.linspace(0.0, 1, 500)
y = ((1 - x) ** 0.62 / (1 + (1 - x) ** 0.41) ** 1.51) / 0.3511112189
y2 = (1 - x) ** 0.5
source = ColumnDataSource(data={"x": x, "y": y})
source2 = ColumnDataSource(data={"x": x, "y2": y2})

plot = figure(
    x_range=(0, 1),
    y_range=(0, 1),
    width=400,
    height=400,
    title=r"Normalised synchrotron loss vs wall reflectivity",
)
plot.xaxis.axis_label = r"Wall reflectivity, \  $$[f_{\text{reflect}}]$$"
plot.yaxis.axis_label = "Normalized synchrotron radiation loss"

plot.line(
    "x",
    "y",
    source=source,
    line_width=3,
    line_alpha=0.6,
    line_color="blue",
    legend_label="Fidone - Correction",
)
plot.line(
    "x",
    "y2",
    source=source2,
    line_width=3,
    line_alpha=0.6,
    line_color="red",
    legend_label="Albajar - Original",
)

plot.legend.title = "Legend"
plot.legend.location = "bottom_left"

# Save the plot as HTML
output_file(
    "psync_albajar_fidone_reflectivity.html",
    title=r"Normalised synchrotron radiation loss vs wall reflectivity",
)
save(plot)
