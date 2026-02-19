import numpy as np
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, output_file, save

x = np.linspace(1.0, 5, 500)
y = 3.12 + 3.5 * (1 / x) ** 1.7
source = ColumnDataSource(data={"x": x, "y": y})

plot = figure(
    x_range=(1, 5),
    y_range=(2, 8),
    width=400,
    height=400,
    title="Menard Normalised Beta Limit",
)
plot.xaxis.axis_label = r"Aspect ratio, \  $$[A]$$"
plot.yaxis.axis_label = r"Normalised beta limit, \  $$[\beta_N]$$"

plot.line("x", "y", source=source, line_width=3, line_alpha=0.6)

# Save the plot as HTML
output_file("menard_beta_norm.html", title="Menard Normalised Beta Limit")
save(plot)
