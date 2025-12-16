import numpy as np
from bokeh.io import output_file, save
from bokeh.layouts import column, row
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
dr_first_wall = Slider(
    start=0.05, end=0.3, value=0.1, step=0.01, title="First Wall Thickness"
)
plasma_divertor_gap = Slider(
    start=0.01, end=1.0, value=1.0, step=0.01, title="Plasma to Divertor Gap"
)
divertor_height = Slider(
    start=0.01, end=1.5, value=1.0, step=0.01, title="Divertor vertical height"
)
dz_blkt_upper = Slider(
    start=0.01, end=1.5, value=0.5, step=0.01, title="Blanket upper height"
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
    y_range=(-8, 8),
    width=650,
    height=600,
    title="D-shaped first wall with divertor",
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

z_fw_inboard_half = (
    (kappa.value * a.value)
    + plasma_divertor_gap.value
    + divertor_height.value
    - dr_first_wall.value
    - dz_blkt_upper.value
)

# Add rectangle to the left of the plasma minus the wall gap
left_rect = Quad(
    left=r0.value - a.value - first_wall_gap_inboard.value - dr_first_wall.value,
    right=r0.value - a.value - first_wall_gap_inboard.value,
    bottom=-z_fw_inboard_half,
    top=z_fw_inboard_half,
    fill_color="grey",
    fill_alpha=1.0,
)
plot.add_glyph(left_rect)


# Outboard first wall as a half-ellipse with thickness 0.1
theta = np.linspace(-np.pi / 2, np.pi / 2, 128)

# Outer ellipse parameters
a_outer = 2 * a.value + 2 * first_wall_gap_inboard.value + dr_first_wall.value
b_outer = (
    (kappa.value * a.value)
    + plasma_divertor_gap.value
    + divertor_height.value
    - dz_blkt_upper.value
)

# Inner ellipse parameters
a_inner = 2 * a.value + 2 * first_wall_gap_inboard.value
b_inner = (
    (kappa.value * a.value)
    + plasma_divertor_gap.value
    + divertor_height.value
    - dr_first_wall.value
    - dz_blkt_upper.value
)

# Parametric equations for half-ellipses
x_outer = (r0.value - a.value - first_wall_gap_inboard.value) + a_outer * np.cos(theta)
y_outer = b_outer * np.sin(theta)
x_inner = (r0.value - a.value - first_wall_gap_inboard.value) + a_inner * np.cos(theta)
y_inner = b_inner * np.sin(theta)

# Concatenate to form a closed path (half-elliptical annulus)
x_fw = np.concatenate([x_outer, x_inner[::-1]])
y_fw = np.concatenate([y_outer, y_inner[::-1]])

fw_source = ColumnDataSource(data={"x": x_fw, "y": y_fw})

fw_patch = plot.patch(
    "x",
    "y",
    source=fw_source,
    fill_color="grey",
    fill_alpha=1.0,
    line_width=0,
    legend_label="FW",
)
# Add a transparent patch between the two blue lines to indicate the divertor region
# Define ColumnDataSource for divertor patch
divertor_source = ColumnDataSource(
    data={
        "x": [
            plot.x_range.start,
            plot.x_range.end,
            plot.x_range.end,
            plot.x_range.start,
        ],
        "y": [
            -(kappa.value * a.value + plasma_divertor_gap.value),
            -(kappa.value * a.value + plasma_divertor_gap.value),
            -(kappa.value * a.value + plasma_divertor_gap.value + divertor_height.value),
            -(kappa.value * a.value + plasma_divertor_gap.value + divertor_height.value),
        ],
    }
)

plot.patch(
    "x",
    "y",
    source=divertor_source,
    fill_color="blue",
    fill_alpha=0.15,
    line_alpha=0,
    legend_label="Divertor region",
)

# Add horizontal line for divertor to plasma gap
hline = Span(
    location=-(kappa.value * a.value + plasma_divertor_gap.value),
    dimension="width",
    line_color="blue",
    line_width=2,
)
plot.add_layout(hline)

# Add horizontal line for divertor to plasma gap
hline_1 = Span(
    location=-(
        kappa.value * a.value + plasma_divertor_gap.value + divertor_height.value
    ),
    dimension="width",
    line_color="blue",
    line_width=2,
)
plot.add_layout(hline_1)

plot.legend.location = "top_left"

callback = CustomJS(
    args={
        "source": source,
        "fw_source": fw_source,
        "divertor_source": divertor_source,
        "r0": r0,
        "a": a,
        "delta": delta,
        "kappa": kappa,
        "square_value": square_value,
        "hline": hline,
        "hline_1": hline_1,
        "first_wall_gap": first_wall_gap_inboard,
        "dr_first_wall": dr_first_wall,
        "plasma_divertor_gap": plasma_divertor_gap,
        "divertor_height": divertor_height,
        "dz_blkt_upper": dz_blkt_upper,
        "left_rect": left_rect,
        "plot": plot,
    },
    code="""
    const A = r0.value;
    const B = a.value;
    const D = delta.value;
    const K = kappa.value;
    const S = square_value;
    const F = first_wall_gap.value;
    const DR = dr_first_wall.value;
    const P = plasma_divertor_gap.value;
    const H = divertor_height.value;
    const DZ = dz_blkt_upper.value;
    const x = source.data['linspace'];
    const R = source.data['x'];
    const Z = source.data['y'];
    for (let i = 0; i < x.length; i++) {
        R[i] = A + B * Math.cos(x[i] + D * Math.sin(x[i]) - S * Math.sin(2 * x[i]));
        Z[i] = K * B * Math.sin(x[i] + S * Math.sin(2 * x[i]));
    }
    const z_fw_inboard_half = (K * B) + P + H - DR - DZ;
    hline.location = -(K*B + P);
    hline_1.location = -(K*B + P + H);
    left_rect.left = A - B - F - DR;
    left_rect.right = A - B - F;
    left_rect.bottom = -z_fw_inboard_half;
    left_rect.top = z_fw_inboard_half;

    // Update outboard FW ellipse
    const n = 128;
    const theta = [];
    for (let i = 0; i < n; i++) {
        theta.push(-Math.PI/2 + i * Math.PI/(n-1));
    }
    const a_outer = 2*B + 2*F + DR;
    const b_outer = (K*B) + P + H - DZ;
    const a_inner = 2*B + 2*F;
    const b_inner = (K*B) + P + H - DR - DZ;
    const x0 = A - B - F;
    const x_outer = [];
    const y_outer = [];
    const x_inner = [];
    const y_inner = [];
    for (let i = 0; i < n; i++) {
        x_outer.push(x0 + a_outer * Math.cos(theta[i]));
        y_outer.push(b_outer * Math.sin(theta[i]));
        x_inner.push(x0 + a_inner * Math.cos(theta[i]));
        y_inner.push(b_inner * Math.sin(theta[i]));
    }
    const x_fw = x_outer.concat(x_inner.slice().reverse());
    const y_fw = y_outer.concat(y_inner.slice().reverse());
    fw_source.data['x'] = x_fw;
    fw_source.data['y'] = y_fw;
    fw_source.change.emit();

    // Update divertor patch
    const x_start = plot.x_range.start;
    const x_end = plot.x_range.end;
    const y_top = -(K*B + P);
    const y_bottom = -(K*B + P + H);
    divertor_source.data['x'] = [x_start, x_end, x_end, x_start];
    divertor_source.data['y'] = [y_top, y_top, y_bottom, y_bottom];
    divertor_source.change.emit();

    source.change.emit();
""",
)

r0.js_on_change("value", callback)
a.js_on_change("value", callback)
delta.js_on_change("value", callback)
kappa.js_on_change("value", callback)
first_wall_gap_inboard.js_on_change("value", callback)
dr_first_wall.js_on_change("value", callback)
plasma_divertor_gap.js_on_change("value", callback)
divertor_height.js_on_change("value", callback)
dz_blkt_upper.js_on_change("value", callback)


# Arrange sliders in two columns with a small gap between them
slider_col1 = column(r0, a, delta, kappa, width=300)
slider_col2 = column(
    first_wall_gap_inboard,
    dr_first_wall,
    plasma_divertor_gap,
    divertor_height,
    dz_blkt_upper,
    width=300,
)
sliders_row = row(slider_col1, slider_col2, spacing=20)

# Save the plot as HTML
save(
    column(
        plot,
        sliders_row,
    ),
    filename="./d_shaped_fw.html",
)
