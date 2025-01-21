import math

import numpy as np
from bokeh.io import output_file, save
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Slider
from bokeh.plotting import figure

square = Slider(start=-1.0, end=1.0, value=0.0, step=0.01, title="Squareness")
r0 = Slider(start=0.1, end=10, value=5.0, step=0.1, title="Major radius")
a = Slider(start=0.01, end=10, value=3.0, step=0.01, title="Minor radius")
delta = Slider(start=0.0, end=1, value=0.5, step=0.01, title="Triangularity | delta")
kappa = Slider(start=0.0, end=10, value=2.0, step=0.01, title="Elongation | kappa")

output_file("plass_shape_plot.html")
x = np.linspace(-np.pi, np.pi, 256)

# Sauter
R = r0.value + a.value * np.cos(
    x + delta.value * np.sin(x) - square.value * np.sin(2 * x)
)
Z = kappa.value * a.value * np.sin(x + square.value * np.sin(2 * x))

# Mirror the Sauter plot
R_mirror = -R

# PROCESS
x1 = (
    2.0 * r0.value * (1.0 + delta.value)
    - a.value * (delta.value**2 + kappa.value**2 - 1.0)
) / (2.0 * (1.0 + delta.value))
x2 = (
    2.0 * r0.value * (delta.value - 1.0)
    - a.value * (delta.value**2 + kappa.value**2 - 1.0)
) / (2.0 * (delta.value - 1.0))
r1 = 0.5 * math.sqrt(
    (a.value**2 * ((delta.value + 1.0) ** 2 + kappa.value**2) ** 2)
    / ((delta.value + 1.0) ** 2)
)
r2 = 0.5 * math.sqrt(
    (a.value**2 * ((delta.value - 1.0) ** 2 + kappa.value**2) ** 2)
    / ((delta.value - 1.0) ** 2)
)
theta1 = np.arcsin((kappa.value * a.value) / r1)
theta2 = np.arcsin((kappa.value * a.value) / r2)

inang = 1.0 / r1
outang = 1.5 / r2

angs1 = np.linspace(-theta1 + np.pi, (inang + theta1) + np.pi, 256, endpoint=True)
angs2 = np.linspace(-(outang + theta2), theta2, 256, endpoint=True)

xs1 = -(r1 * np.cos(angs1) - x1)
ys1 = r1 * np.sin(angs1)
xs2 = -(r2 * np.cos(angs2) - x2)
ys2 = r2 * np.sin(angs2)

# Mirror the PROCESS plot
xs1_mirror = -xs1
xs2_mirror = -xs2

source = ColumnDataSource(
    data={
        "x": R,
        "y": Z,
        "x_mirror": R_mirror,
        "y_mirror": Z,
        "x1": xs1,
        "y1": ys1,
        "x2": xs2,
        "y2": ys2,
        "x1_mirror": xs1_mirror,
        "y1_mirror": ys1,
        "x2_mirror": xs2_mirror,
        "y2_mirror": ys2,
        "linspace": x,
    }
)

plot = figure(
    x_range=(-10, 10), y_range=(-10, 10), width=800, height=800, title="Plasma Shape"
)
plot.xaxis.axis_label = r"Radius"
plot.yaxis.axis_label = r"Height"

plot.line("x", "y", source=source, line_width=3, line_alpha=0.6, legend_label="Sauter")
plot.line("x_mirror", "y_mirror", source=source, line_width=3, line_alpha=0.6)
plot.line(
    "x1",
    "y1",
    source=source,
    line_width=3,
    line_alpha=0.6,
    color="red",
    legend_label="PROCESS",
)
plot.line(
    "x1_mirror", "y1_mirror", source=source, line_width=3, line_alpha=0.6, color="red"
)
plot.line("x2", "y2", source=source, line_width=3, line_alpha=0.6, color="red")
plot.line(
    "x2_mirror", "y2_mirror", source=source, line_width=3, line_alpha=0.6, color="red"
)

plot.legend.location = "top_left"

callback = CustomJS(
    args={
        "source": source,
        "r0": r0,
        "a": a,
        "delta": delta,
        "kappa": kappa,
        "square": square,
    },
    code="""
    const A = r0.value;
    const B = a.value;
    const D = delta.value;
    const K = kappa.value;
    const S = square.value;
    const x = source.data['linspace'];
    const R = source.data['x'];
    const Z = source.data['y'];
    const R_mirror = source.data['x_mirror'];
    const Z_mirror = source.data['y_mirror'];
    const x1 = (2.0 * A * (1.0 + D) - B * (D**2 + K**2 - 1.0)) / (2.0 * (1.0 + D));
    const x2 = (2.0 * A * (D - 1.0) - B * (D**2 + K**2 - 1.0)) / (2.0 * (D - 1.0));
    const r1 = 0.5 * Math.sqrt((B**2 * ((D + 1.0) ** 2 + K**2) ** 2) / ((D + 1.0) ** 2));
    const r2 = 0.5 * Math.sqrt((B**2 * ((D - 1.0) ** 2 + K**2) ** 2) / ((D - 1.0) ** 2));
    const theta1 = Math.asin((K * B) / r1);
    const theta2 = Math.asin((K * B) / r2);
    const inang = 1.0 / r1;
    const outang = 1.5 / r2;
    const angs1 = Array.from({length: 256}, (_, i) => -theta1 + Math.PI + i * ((inang + 2*theta1) / 255));
    const angs2 = Array.from({length: 256}, (_, i) => -(outang + theta2) + i * ((theta2 + (outang + theta2)) / 255));
    const xs1 = source.data['x1'];
    const ys1 = source.data['y1'];
    const xs2 = source.data['x2'];
    const ys2 = source.data['y2'];
    const xs1_mirror = source.data['x1_mirror'];
    const ys1_mirror = source.data['y1_mirror'];
    const xs2_mirror = source.data['x2_mirror'];
    const ys2_mirror = source.data['y2_mirror'];
    for (let i = 0; i < x.length; i++) {
        R[i] = A + B * Math.cos(x[i] + D * Math.sin(x[i]) - S * Math.sin(2 * x[i]));
        Z[i] = K * B * Math.sin(x[i] + S * Math.sin(2 * x[i]));
        R_mirror[i] = -R[i];
        Z_mirror[i] = Z[i];
        xs1[i] = -(r1 * Math.cos(angs1[i]) - x1);
        ys1[i] = r1 * Math.sin(angs1[i]);
        xs2[i] = -(r2 * Math.cos(angs2[i]) - x2);
        ys2[i] = r2 * Math.sin(angs2[i]);
        xs1_mirror[i] = -xs1[i];
        ys1_mirror[i] = ys1[i];
        xs2_mirror[i] = -xs2[i];
        ys2_mirror[i] = ys2[i];
    }
    source.change.emit();
""",
)

r0.js_on_change("value", callback)
a.js_on_change("value", callback)
delta.js_on_change("value", callback)
kappa.js_on_change("value", callback)
square.js_on_change("value", callback)

# Save the plot as HTML
save(row(plot, column(r0, a, delta, kappa, square)), filename="./plasma_plot.html")
