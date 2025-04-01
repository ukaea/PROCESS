import numpy as np
from bokeh.io import output_file, save
from bokeh.layouts import column
from bokeh.models import ColumnDataSource, CustomJS, Quad, Slider, Span
from bokeh.plotting import figure

# Fixed value for square
square_value = 0.0

r0 = Slider(start=0.1, end=10, value=5.0, step=0.1, title="Major radius")
a = Slider(start=0.01, end=2.5, value=2.0, step=0.01, title="Minor radius")
delta = Slider(start=0.0, end=1, value=0.5, step=0.01, title="Triangularity | delta")
kappa = Slider(start=0.0, end=2.5, value=2.0, step=0.01, title="Elongation | kappa")
first_wall_gap_inboard = Slider(
    start=0.01, end=0.5, value=0.3, step=0.01, title="First Wall Gap"
)
plasma_divertor_gap = Slider(
    start=0.01, end=1.0, value=1.0, step=0.01, title="Plasma to Divertor Gap"
)
divertor_height = Slider(
    start=0.01, end=1.5, value=1.0, step=0.01, title="Divertor vertical height"
)

output_file("divtart_divertor.html")
x = np.linspace(-np.pi, np.pi, 256)

# Sauter
R = r0.value + a.value * np.cos(
    x + delta.value * np.sin(x) - square_value * np.sin(2 * x)
)
Z = kappa.value * a.value * np.sin(x + square_value * np.sin(2 * x))

source = ColumnDataSource(
    data={
        "x": R,
        "y": Z,
        "linspace": x,
    }
)

plot = figure(
    x_range=(0, 10),
    y_range=(0, 10),
    width=650,
    height=800,
    title="divtart divertor box",
)
plot.xaxis.axis_label = r"Radius"
plot.yaxis.axis_label = r"Height"

plot.patch(
    "x",
    "y",
    source=source,
    fill_color="pink",
    line_width=0,
    fill_alpha=1.0,
    legend_label="Plasma",
)

# Add vertical line for inner radius of divertor region
vline_inner = Span(
    location=r0.value - (a.value * delta.value) - (3.0 * first_wall_gap_inboard.value),
    dimension="height",
    line_color="red",
    line_width=2,
)
plot.add_layout(vline_inner)

# Add vertical line for outer radius of divertor region
vline = Span(
    location=r0.value + a.value, dimension="height", line_color="red", line_width=2
)
plot.add_layout(vline)

# Add horizontal line for divertor to plasma gap
hline = Span(
    location=kappa.value * a.value + plasma_divertor_gap.value,
    dimension="width",
    line_color="blue",
    line_width=2,
)
plot.add_layout(hline)

# Add rectangle to the left of the plasma minus the wall gap
left_rect = Quad(
    left=0.0,
    right=r0.value - a.value - first_wall_gap_inboard.value,
    bottom=-kappa.value * a.value,
    top=10.0,
    fill_color="grey",
    fill_alpha=1.0,
)
plot.add_glyph(left_rect)

# Add triangle for the divertor
divertor_triangle = plot.patch(
    [
        r0.value - (a.value * delta.value) - (3.0 * first_wall_gap_inboard.value),
        r0.value - (a.value * delta.value) - (3.0 * first_wall_gap_inboard.value),
        r0.value + a.value,
    ],
    [
        kappa.value * a.value + plasma_divertor_gap.value,
        kappa.value * a.value + plasma_divertor_gap.value + divertor_height.value,
        kappa.value * a.value + plasma_divertor_gap.value + divertor_height.value,
    ],
    fill_color="green",
    fill_alpha=0.6,
    line_width=2,
    legend_label="Divertor box",
)

plot.legend.location = "top_left"

callback = CustomJS(
    args={
        "source": source,
        "r0": r0,
        "a": a,
        "delta": delta,
        "kappa": kappa,
        "square_value": square_value,
        "vline": vline,
        "vline_inner": vline_inner,
        "hline": hline,
        "first_wall_gap": first_wall_gap_inboard,
        "plasma_divertor_gap": plasma_divertor_gap,
        "divertor_height": divertor_height,
        "left_rect": left_rect,
        "divertor_triangle": divertor_triangle,
    },
    code="""
    const A = r0.value;
    const B = a.value;
    const D = delta.value;
    const K = kappa.value;
    const S = square_value;
    const F = first_wall_gap.value;
    const P = plasma_divertor_gap.value;
    const H = divertor_height.value;
    const x = source.data['linspace'];
    const R = source.data['x'];
    const Z = source.data['y'];
    for (let i = 0; i < x.length; i++) {
        R[i] = A + B * Math.cos(x[i] + D * Math.sin(x[i]) - S * Math.sin(2 * x[i]));
        Z[i] = K * B * Math.sin(x[i] + S * Math.sin(2 * x[i]));
    }
    vline.location = A + B;
    vline_inner.location = A - B*D - 3.0*F;
    hline.location = K*B + P;
    left_rect.left = 0.0;
    left_rect.right = A - B -F;
    divertor_triangle.data_source.data['x'] = [
        A - B*D - 3.0*F,
        A - B*D - 3.0*F,
        A + B
    ];
    divertor_triangle.data_source.data['y'] = [
        K*B + P,
        K*B + P + H,
        K*B + P + H
    ];
    source.change.emit();
    divertor_triangle.data_source.change.emit();
""",
)

r0.js_on_change("value", callback)
a.js_on_change("value", callback)
delta.js_on_change("value", callback)
kappa.js_on_change("value", callback)
first_wall_gap_inboard.js_on_change("value", callback)
plasma_divertor_gap.js_on_change("value", callback)
divertor_height.js_on_change("value", callback)
# Save the plot as HTML
save(
    column(
        plot,
        r0,
        a,
        delta,
        kappa,
        first_wall_gap_inboard,
        plasma_divertor_gap,
        divertor_height,
    ),
    filename="./divtart_divertor.html",
)
