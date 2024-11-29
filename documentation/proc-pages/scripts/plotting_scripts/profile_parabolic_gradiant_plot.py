import numpy as np
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import figure
from bokeh.io import save


x = np.linspace(0, 2, 500)
y = -2 * x * (1 / np.sqrt(2 * x - 1)) * (1 - (1 / np.sqrt(2 * x - 1)) ^ 2) ^ (x - 1)

source = ColumnDataSource(data=dict(x=x, y=y))

plot = figure(x_range=(0, 2), y_range=(0, -5), width=400, height=400, title="")
plot.xaxis.axis_label = r"$$ \alpha $$"
plot.yaxis.axis_label = r"Gradiant Value"

plot.line("x", "y", source=source, line_width=3, line_alpha=0.6)

n0 = Slider(start=0.1, end=10, value=5, step=0.1, title="Plasma centre value | n0")
alpha = Slider(start=0.01, end=10, value=2, step=0.01, title="Profile Index  | alphan")


callback = CustomJS(
    args=dict(source=source, n0=n0, alpha=alpha),
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
save(
    row(plot, column(n0, alpha)),
    filename="/PROCESS/profile_parabolic_gradiant_plot.html",
)
