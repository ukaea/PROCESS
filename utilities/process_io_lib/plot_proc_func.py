"""

  Functions for creating one-page summary

  James Morris
  26/03/2014

  CCFE

  Compatible with PROCESS version ???

"""

import math
import scipy as sp
import numpy as np

try:
    import process_io_lib.process_dicts as proc_dict
except ImportError:
    print(
        "The Python dictionaries have not yet been created. Please run \
'make dicts'!"
    )
    exit()

RADIAL_BUILD = [
    "bore",
    "ohcth",
    "gapoh",
    "tfcth",
    "gapds",
    "d_vv_in",
    "shldith",
    "blnkith",
    "fwith",
    "scrapli",
    "rminori",
    "rminoro",
    "scraplo",
    "fwoth",
    "blnkoth",
    "shldoth",
    "d_vv_out",
    "gapsto",
    "tfthko",
]

VERTICAL_BUILD = [
    "tfcth",
    "vgap2",
    "d_vv_top",
    "shldtth",
    "blnktth",
    "top first wall vertical thickness (m)",
    "top scrape-off vertical thickness (m)",
    "rminor*kappa",
    "midplane",
    "rminor*kappa",
    "vgap",
    "divfix",
    "shldtth",
    "d_vv_top",
    "vgap2",
    "tfcth",
]


FILLCOLS = ("0.8", "0", "#00ff00", "#00ffff", "#ff0000", "#ff00ff")

ANIMATION_INFO = [
    ("rmajor", "Major radius", "m"),
    ("rminor", "Minor radius", "m"),
    ("aspect", "Aspect ratio", ""),
]
# ("", "", ""),
# ("", "", "")]


def plotdh(axis, r0, a, delta, kap):
    """Plots half a thin D-section.

    Arguments:
        axis --> axis object to plot to
        r0 --> major radius
        a --> minor radius
        delta --> triangularity
        kap --> elongation

    Returns:
        rs --> radial coordinates of D-section
        zs --> vertical coordinates of D-section

    """
    angs = np.linspace(0, np.pi, 50, endpoint=True)
    rs = r0 + a * np.cos(angs + delta * np.sin(1.0 * angs))
    zs = kap * a * np.sin(angs)
    axis.plot(rs, zs, color="black")
    return rs, zs


def plotdhgap(axis, inpt, outpt, inthk, outthk, toppt, topthk, delta, col):
    """Plots half a thick D-section with a gap.

    Arguments:
        axis --> axis object to plot to
        inpt --> inner points
        outpt --> outer points
        inthk --> inner thickness
        outthk --> outer thickness
        toppt --> top points
        topthk --> top thickness
        delta --> triangularity
        col --> color for fill

    """
    arc = np.pi / 4.0
    r01 = (inpt + outpt) / 2.0
    r02 = (inpt + inthk + outpt - outthk) / 2.0
    a1 = r01 - inpt
    a2 = r02 - inpt - inthk
    kap1 = toppt / a1
    kap2 = (toppt - topthk) / a2
    # angs = ((np.pi/2.) - arc/2.) * findgen(50)/49.
    angs = np.linspace(0.0, (np.pi / 2.0) - arc / 2.0, 50, endpoint=True)
    rs1 = r01 + a1 * np.cos(angs + delta * np.sin(angs))
    zs1 = kap1 * a1 * np.sin(angs)
    rs2 = r02 + a2 * np.cos(angs + delta * np.sin(angs))
    zs2 = kap2 * a2 * np.sin(angs)
    # angs = !pi + ((!pi/2.) - arc) * findgen(50)/49.
    angs = np.linspace(np.pi, np.pi + ((np.pi / 2.0) - arc), 50, endpoint=True)
    rs3 = r01 + a1 * np.cos(angs + delta * np.sin(angs))
    zs3 = kap1 * a1 * np.sin(angs)
    rs4 = r02 + a2 * np.cos(angs + delta * np.sin(angs))
    zs4 = kap2 * a2 * np.sin(angs)

    axis.plot(
        np.concatenate([rs1, rs2[::-1]]),
        np.concatenate([zs1, zs2[::-1]]),
        color="black",
    )
    axis.plot(
        np.concatenate([rs3, rs4[::-1]]),
        -np.concatenate([zs3, zs4[::-1]]),
        color="black",
    )
    axis.fill(
        np.concatenate([rs1, rs2[::-1]]), np.concatenate([zs1, zs2[::-1]]), color=col
    )
    axis.fill(
        np.concatenate([rs3, rs4[::-1]]), -np.concatenate([zs3, zs4[::-1]]), color=col
    )


def plot_plasma(axis, mfile_data, scan):
    """Plots the plasma boundary arcs.

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use

    """

    r0 = mfile_data.data["rmajor"].get_scan(scan)
    a = mfile_data.data["rminor"].get_scan(scan)
    delta = 1.5 * mfile_data.data["triang95"].get_scan(scan)
    kappa = (1.1 * mfile_data.data["kappa95"].get_scan(scan)) + 0.04
    i_single_null = mfile_data.data["i_single_null"].get_scan(scan)

    x1 = (2.0 * r0 * (1.0 + delta) - a * (delta**2 + kappa**2 - 1.0)) / (
        2.0 * (1.0 + delta)
    )
    x2 = (2.0 * r0 * (delta - 1.0) - a * (delta**2 + kappa**2 - 1.0)) / (
        2.0 * (delta - 1.0)
    )
    r1 = 0.5 * math.sqrt(
        (a**2 * ((delta + 1.0) ** 2 + kappa**2) ** 2) / ((delta + 1.0) ** 2)
    )
    r2 = 0.5 * math.sqrt(
        (a**2 * ((delta - 1.0) ** 2 + kappa**2) ** 2) / ((delta - 1.0) ** 2)
    )
    theta1 = sp.arcsin((kappa * a) / r1)
    theta2 = sp.arcsin((kappa * a) / r2)
    inang = 1.0 / r1
    outang = 1.5 / r2
    if i_single_null == 0:
        angs1 = np.linspace(
            -(inang + theta1) + np.pi, (inang + theta1) + np.pi, 256, endpoint=True
        )
        angs2 = np.linspace(-(outang + theta2), (outang + theta2), 256, endpoint=True)
    elif i_single_null < 0:
        angs1 = np.linspace(
            -(inang + theta1) + np.pi, theta1 + np.pi, 256, endpoint=True
        )
        angs2 = np.linspace(-theta2, (outang + theta2), 256, endpoint=True)
    else:
        angs1 = np.linspace(
            -theta1 + np.pi, (inang + theta1) + np.pi, 256, endpoint=True
        )
        angs2 = np.linspace(-(outang + theta2), theta2, 256, endpoint=True)
    xs1 = -(r1 * np.cos(angs1) - x1)
    ys1 = r1 * np.sin(angs1)
    xs2 = -(r2 * np.cos(angs2) - x2)
    ys2 = r2 * np.sin(angs2)
    axis.plot(xs1, ys1, color="black")
    axis.plot(xs2, ys2, color="black")


def plot_centre_cross(axis, mfile_data, scan):
    """Function to plot centre cross on plot

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use

    """
    rmajor = mfile_data.data["rmajor"].get_scan(scan)
    axis.plot(
        [rmajor - 0.25, rmajor + 0.25, rmajor, rmajor, rmajor],
        [0, 0, 0, 0.25, -0.25],
        color="black",
    )


def cumulative_radial_build(section, mfile_data, scan):
    """Function for calculating the cumulative radial build up to and
    including the given section.

    Arguments:
        section --> section of the radial build to go up to
        mfile_data --> MFILE data object
        scan --> scan number to use

    Returns:
        cumulative_build --> cumulative radial build up to section given

    """

    cumulative_build = 0
    for item in RADIAL_BUILD:
        if item == "rminori" or item == "rminoro":
            cumulative_build += mfile_data.data["rminor"].get_scan(scan)
        elif "d_vv_" in item:
            cumulative_build += mfile_data.data["d_vv_in"].get_scan(scan)
        else:
            cumulative_build += mfile_data.data[item].get_scan(scan)
        if item == section:
            break
    return cumulative_build


def cumulative_vertical_build(mfile_data, scan):
    """Function for calculating the cumulative vertical build

    Arguments:
        mfile_data --> MFILE data object
        scan --> scan number to use

    Returns:
        cumulative_build --> vertical build list

    """

    cumulative_build = list()
    cumulative_build.append(0.0)

    # Top
    top_build = 0
    start = False
    for i in range(len(VERTICAL_BUILD) - 1, -1, -1):
        if not start:
            if VERTICAL_BUILD[i] == "midplane":
                start = True
        else:
            top_build += mfile_data.data[VERTICAL_BUILD[i]].get_scan(scan)
            cumulative_build.insert(0, top_build)

    # bottom
    start = False
    bottom_build = 0
    for j in range(0, len(VERTICAL_BUILD), 1):
        if not start:
            if VERTICAL_BUILD[j] == "midplane":
                start = True
        else:
            bottom_build -= mfile_data.data[VERTICAL_BUILD[j]].get_scan(scan)
            cumulative_build.append(bottom_build)

    return cumulative_build


def plot_machine_pic(axis, mfile_data, scan=-1):
    """Function to plot machine picture

    Arguments:
      axis --> axis object to add plot to
      mfile_data --> MFILE data object
      scan --> scan number to use

    """

    xmin = 0
    xmax = 20
    ymin = -15
    ymax = 15
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)
    axis.set_xlabel("R / m")
    axis.set_ylabel("Z / m")
    axis.set_title("Radial build")

    if mfile_data.data["i_single_null"].get_scan(scan):
        plot_i_single_null_cryo(axis, mfile_data, scan)
        plot_shield_i_single_null(axis, mfile_data, scan)

    plot_plasma(axis, mfile_data, scan)
    plot_centre_cross(axis, mfile_data, scan)


def plot_i_single_null_cryo(axis, mfile_data, scan):
    """Function to plot top of cryostat

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE data object
        scan --> scan number to use

    """
    triang = mfile_data.data["triang95"].get_scan(scan)

    # cryostat
    temp_array_1 = ()
    temp_array_2 = ()
    # Outer side
    radx = (
        cumulative_radial_build("ddwo", mfile_data, scan)
        + cumulative_radial_build("gapds", mfile_data, scan)
    ) / 2.0
    rminx = (
        cumulative_radial_build("ddwo", mfile_data, scan)
        - cumulative_radial_build("gapds", mfile_data, scan)
    ) / 2.0

    a = cumulative_vertical_build(mfile_data, scan)
    kapx = a[2] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_1 = temp_array_1 + ((rs, zs))

    kapx = a[13] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_2 = temp_array_2 + ((rs, zs))

    # Inner side
    radx = (
        cumulative_radial_build("shldoth", mfile_data, scan)
        + cumulative_radial_build("d_vv_out", mfile_data, scan)
    ) / 2.0
    rminx = (
        cumulative_radial_build("shldoth", mfile_data, scan)
        - cumulative_radial_build("d_vv_out", mfile_data, scan)
    ) / 2.0

    a = cumulative_vertical_build(mfile_data, scan)
    kapx = a[3] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_1 = temp_array_1 + ((rs, zs))
    kapx = a[12] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_2 = temp_array_2 + ((rs, zs))

    rs = np.concatenate([temp_array_1[0], temp_array_1[2][::-1]])
    zs = np.concatenate([temp_array_1[1], temp_array_1[3][::-1]])
    axis.fill(rs, zs, color=FILLCOLS[1])

    rs = np.concatenate([temp_array_2[0], temp_array_2[2][::-1]])
    zs = np.concatenate([temp_array_2[1], temp_array_2[3][::-1]])
    axis.fill(rs, zs, color=FILLCOLS[1])


def plot_shield_i_single_null(axis, mfile_data, scan):
    """Function to plot single null case of shield

    Arguments:
      axis --> axis object to plot to
      mfile_data --> MFILE.DAT object
      scan --> scan number to use

    """
    triang = mfile_data.data["triang95"].get_scan(scan)
    blnkith = mfile_data.data["blnkith"].get_scan(scan)
    blnkoth = mfile_data.data["blnkoth"].get_scan(scan)
    fwith = mfile_data.data["fwith"].get_scan(scan)
    fwoth = mfile_data.data["fwoth"].get_scan(scan)
    blnktth = mfile_data.data["blnktth"].get_scan(scan)
    tfwvt = mfile_data.data["top first wall vertical thickness (m)"].get_scan(scan)
    c_shldith = cumulative_radial_build("shldith", mfile_data, scan)
    c_blnkoth = cumulative_radial_build("blnkoth", mfile_data, scan)
    c_blnkith = cumulative_radial_build("blnkith", mfile_data, scan)
    c_fwoth = cumulative_radial_build("fwoth", mfile_data, scan)
    a = cumulative_vertical_build(mfile_data, scan)

    temp_array_1 = ()
    temp_array_2 = ()

    radx = (
        cumulative_radial_build("shldoth", mfile_data, scan)
        + cumulative_radial_build("d_vv_out", mfile_data, scan)
    ) / 2.0
    rminx = (
        cumulative_radial_build("shldoth", mfile_data, scan)
        - cumulative_radial_build("d_vv_out", mfile_data, scan)
    ) / 2.0
    kapx = a[3] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_1 = temp_array_1 + ((rs, zs))
    kapx = a[12] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_2 = temp_array_2 + ((rs, zs))

    radx = (
        cumulative_radial_build("blnkoth", mfile_data, scan)
        + cumulative_radial_build("shldith", mfile_data, scan)
    ) / 2.0
    rminx = (
        cumulative_radial_build("blnkoth", mfile_data, scan)
        - cumulative_radial_build("shldith", mfile_data, scan)
    ) / 2.0
    kapx = a[4] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_1 = temp_array_1 + ((rs, zs))
    kapx = a[11] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_2 = temp_array_2 + ((rs, zs))

    radx = (
        cumulative_radial_build("fwoth", mfile_data, scan)
        + cumulative_radial_build("blnkith", mfile_data, scan)
    ) / 2.0
    rminx = (
        cumulative_radial_build("fwoth", mfile_data, scan)
        - cumulative_radial_build("blnkith", mfile_data, scan)
    ) / 2.0
    kapx = a[5] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_1 = temp_array_1 + ((rs, zs))

    radx = (
        cumulative_radial_build("scraplo", mfile_data, scan)
        + cumulative_radial_build("fwith", mfile_data, scan)
    ) / 2.0
    rminx = (
        cumulative_radial_build("scraplo", mfile_data, scan)
        - cumulative_radial_build("fwith", mfile_data, scan)
    ) / 2.0
    kapx = a[6] / rminx
    (rs, zs) = plotdh(axis, radx, rminx, triang, kapx)
    temp_array_1 = temp_array_1 + ((rs, zs))

    # Top
    rs = np.concatenate([temp_array_1[0], temp_array_1[2][::-1]])
    zs = np.concatenate([temp_array_1[1], temp_array_1[3][::-1]])
    axis.fill(rs, zs, color=FILLCOLS[2])

    # for i in range(len(rs)):
    # print(rs[i], zs[i])

    rs = np.concatenate([temp_array_1[0 + 2 * 1], temp_array_1[2 + 2 * 1][::-1]])
    zs = np.concatenate([temp_array_1[1 + 2 * 1], temp_array_1[3 + 2 * 1][::-1]])
    axis.fill(rs, zs, color=FILLCOLS[1 + 2])

    rs = np.concatenate([temp_array_1[0 + 2 * 2], temp_array_1[2 + 2 * 2][::-1]])
    zs = np.concatenate([temp_array_1[1 + 2 * 2], temp_array_1[3 + 2 * 2][::-1]])
    axis.fill(rs, zs, color=FILLCOLS[2 + 2])

    # Bottom
    rs = np.concatenate([temp_array_2[0], temp_array_2[2][::-1]])
    zs = np.concatenate([temp_array_2[1], temp_array_2[3][::-1]])
    axis.fill(rs, zs, color=FILLCOLS[2])

    plotdhgap(
        axis,
        c_shldith,
        c_blnkoth,
        blnkith,
        blnkoth,
        a[11],
        -blnktth,
        triang,
        FILLCOLS[3],
    )
    plotdhgap(
        axis,
        c_blnkith,
        c_fwoth,
        fwith,
        fwoth,
        a[11] + blnktth,
        -tfwvt,
        triang,
        FILLCOLS[4],
    )


def angle_check(angle1, angle2):
    """Function to perform TF coil angle check"""
    if angle1 > 1:
        angle1 = 1
    if angle1 < -1:
        angle1 = -1
    if angle2 > 1:
        angle2 = 1
    if angle2 < -1:
        angle2 = -1
    return angle1, angle2


def plot_tf_coils(axis, mfile_data, scan):
    """Function to plot TF coils

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """

    vert_build = cumulative_vertical_build(mfile_data, scan)

    # Arc points
    x_1 = mfile_data.data["xarc(1)"].get_scan(scan)
    y_1 = mfile_data.data["yarc(1)"].get_scan(scan)
    x_2 = mfile_data.data["xarc(2)"].get_scan(scan)
    y_2 = mfile_data.data["yarc(2)"].get_scan(scan)
    x_3 = mfile_data.data["xarc(3)"].get_scan(scan)
    y_3 = mfile_data.data["yarc(3)"].get_scan(scan)
    x_4 = mfile_data.data["xarc(4)"].get_scan(scan)
    y_4 = mfile_data.data["yarc(4)"].get_scan(scan)
    x_5 = mfile_data.data["xarc(5)"].get_scan(scan)
    # y_5 = mfile_data.data["yarc(5)"].get_scan(scan)

    # Arc centres
    x_c_1 = mfile_data.data["xctfc(1)"].get_scan(scan)
    y_c_1 = mfile_data.data["yctfc(1)"].get_scan(scan)
    x_c_2 = mfile_data.data["xctfc(2)"].get_scan(scan)
    y_c_2 = mfile_data.data["yctfc(2)"].get_scan(scan)
    x_c_3 = mfile_data.data["xctfc(3)"].get_scan(scan)
    y_c_3 = mfile_data.data["yctfc(3)"].get_scan(scan)
    x_c_4 = mfile_data.data["xctfc(4)"].get_scan(scan)
    y_c_4 = mfile_data.data["yctfc(4)"].get_scan(scan)

    # Arc radii
    rad_1 = np.sqrt((x_1 - x_c_1) ** 2 + (y_1 - y_c_1) ** 2)
    rad_2 = np.sqrt((x_2 - x_c_2) ** 2 + (y_2 - y_c_2) ** 2)
    rad_3 = np.sqrt((x_3 - x_c_3) ** 2 + (y_3 - y_c_3) ** 2)
    rad_4 = np.sqrt((x_4 - x_c_4) ** 2 + (y_4 - y_c_4) ** 2)

    # Arc angles
    a_1_1 = (x_1 - x_c_1) / rad_1
    a_1_2 = (x_2 - x_c_1) / rad_1
    a_1_1, a_1_2 = angle_check(a_1_1, a_1_2)
    angle_1_1 = sp.arcsin(a_1_1)
    angle_1_2 = sp.arcsin(a_1_2)

    a_2_1 = (x_2 - x_c_2) / rad_2
    a_2_2 = (x_3 - x_c_2) / rad_2
    a_2_1, a_2_2 = angle_check(a_2_1, a_2_2)
    angle_2_1 = sp.arcsin(a_2_1)
    angle_2_2 = sp.arcsin(a_2_2)

    a_3_1 = (x_3 - x_c_3) / rad_3
    a_3_2 = (x_4 - x_c_3) / rad_3
    a_3_1, a_3_2 = angle_check(a_3_1, a_3_2)
    angle_3_1 = sp.arcsin(a_3_1)
    angle_3_2 = sp.arcsin(a_3_2)

    a_4_1 = (x_4 - x_c_4) / rad_4
    a_4_2 = (x_5 - x_c_4) / rad_4
    a_4_1, a_4_2 = angle_check(a_4_1, a_4_2)
    angle_4_1 = sp.arcsin(a_4_1)
    angle_4_2 = sp.arcsin(a_4_2)

    in_x_1 = (
        rad_1 * np.sin(np.linspace(angle_1_1, angle_1_2, 30, endpoint=True)) + x_c_1
    )
    in_x_2 = (
        rad_2 * np.sin(np.linspace(angle_2_1, angle_2_2, 30, endpoint=True)) + x_c_2
    )
    in_x_3 = (
        rad_3 * np.sin(np.linspace(angle_3_1, angle_3_2, 30, endpoint=True)) + x_c_3
    )
    in_x_4 = (
        rad_4 * np.sin(np.linspace(angle_4_1, angle_4_2, 30, endpoint=True)) + x_c_4
    )
    in_x = np.concatenate((in_x_1, in_x_2, in_x_3, in_x_4))

    in_y_1 = (
        rad_1 * np.cos(np.linspace(angle_1_1, angle_1_2, 30, endpoint=True)) + y_c_1
    )
    in_y_2 = (
        rad_2 * np.cos(np.linspace(angle_2_1, angle_2_2, 30, endpoint=True)) + y_c_2
    )
    in_y_3 = (
        rad_3 * np.cos(np.linspace(angle_3_1, angle_3_2, 30, endpoint=True)) + y_c_3
    )
    in_y_4 = (
        rad_4 * np.cos(np.linspace(angle_4_1, angle_4_2, 30, endpoint=True)) + y_c_4
    )
    in_y = np.concatenate((in_y_1, in_y_2, in_y_3, in_y_4))

    in_x = np.concatenate([in_x, in_x[::-1]])
    in_y = np.concatenate([in_y, -in_y[::-1]])

    # in_width = cumulative_radial_build("gapsto", mfile_data, scan) - \
    #    cumulative_radial_build("tfcth", mfile_data, scan)
    # out_width = in_width + mfile_data.data["tfcth"].get_scan(scan) + \
    #    mfile_data.data["tfthko"].get_scan(scan)
    # out_x = ((in_x - cumulative_radial_build("tfcth", mfile_data, scan)) *
    #         (out_width/in_width))

    centre_in_x = (
        cumulative_radial_build("tfcth", mfile_data, scan)
        + cumulative_radial_build("gapsto", mfile_data, scan)
    ) / 2.0
    centre_out_x = (
        cumulative_radial_build("gapoh", mfile_data, scan)
        + cumulative_radial_build("tfthko", mfile_data, scan)
    ) / 2.0
    in_width = cumulative_radial_build(
        "gapsto", mfile_data, scan
    ) - cumulative_radial_build("tfcth", mfile_data, scan)
    out_width = cumulative_radial_build(
        "tfthko", mfile_data, scan
    ) - cumulative_radial_build("gapoh", mfile_data, scan)

    out_x = ((in_x - centre_in_x) * (out_width / in_width)) + centre_out_x
    extern = vert_build[7] / vert_build[6]
    if vert_build[-1]:
        extern = (vert_build[0] - vert_build[15]) / (vert_build[1] - vert_build[14])
    out_y = in_y * extern

    shift = (vert_build[0] + vert_build[15]) / 2.0
    out_y += shift
    in_y += shift

    ins_x = np.concatenate((in_x, out_x[::-1]))
    ins_y = np.concatenate((in_y, out_y[::-1]))

    axis.fill(ins_x, ins_y, color=FILLCOLS[0])


def plot_pf_coils(axis, mfile_data, scan):
    """Function to plot PF coils

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """

    coils_r = []
    coils_z = []
    coils_dr = []
    coils_dz = []
    coil_text = []

    # Number of coils (1 is OH coil)
    number_of_coils = 0
    for item in mfile_data.data.keys():
        if "rpf(" in item:
            number_of_coils += 1

    # Get OH coil
    coils_r.append(mfile_data.data["rpf(nohc)"].get_scan(scan))
    coils_z.append(mfile_data.data["zpf(nohc)"].get_scan(scan))
    coils_dr.append(mfile_data.data["ohdr"].get_scan(scan))
    coils_dz.append(mfile_data.data["ohdz"].get_scan(scan))
    coil_text.append("OH")

    # Rest of the coils
    for coil in range(1, number_of_coils):
        coils_r.append(mfile_data.data["rpf({:02})".format(coil)].get_scan(scan))
        coils_z.append(mfile_data.data["zpf({:02})".format(coil)].get_scan(scan))
        coils_dr.append(mfile_data.data["pfdr{:02}".format(coil)].get_scan(scan))
        coils_dz.append(mfile_data.data["pfdz{:02}".format(coil)].get_scan(scan))
        coil_text.append(str(coil))

    for i in range(len(coils_r)):
        r_1 = coils_r[i] - 0.5 * coils_dr[i]
        z_1 = coils_z[i] - 0.5 * coils_dz[i]
        r_2 = coils_r[i] - 0.5 * coils_dr[i]
        z_2 = coils_z[i] + 0.5 * coils_dz[i]
        r_3 = coils_r[i] + 0.5 * coils_dr[i]
        z_3 = coils_z[i] + 0.5 * coils_dz[i]
        r_4 = coils_r[i] + 0.5 * coils_dr[i]
        z_4 = coils_z[i] - 0.5 * coils_dz[i]
        r_5 = coils_r[i] - 0.5 * coils_dr[i]
        z_5 = coils_z[i] - 0.5 * coils_dz[i]

        r_points = [r_1, r_2, r_3, r_4, r_5]
        z_points = [z_1, z_2, z_3, z_4, z_5]
        axis.plot(r_points, z_points, color="black")
        axis.text(
            coils_r[i],
            coils_z[i],
            coil_text[i],
            ha="center",
            va="center",
            fontsize="smaller",
        )


def plot_info(axis, data, mfile_data, scan):
    """Function to plot data in written form on a matplotlib plot.

    Arguments:
        axis --> axis object to plot to
        data --> plot information
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """
    eqpos = 0.7
    for i in range(len(data)):
        axis.text(0, -i, data[i][1], ha="left", va="center")
        if isinstance(data[i][0], str):
            if data[i][0] == "":
                axis.text(eqpos, -i, "\n", ha="left", va="center")
            elif data[i][0][0] == "#":
                axis.text(
                    -0.05, -i, "{}\n".format(data[i][0][1:]), ha="left", va="center"
                )
            elif data[i][0][0] == "!":
                value = data[i][0][1:]
                axis.text(
                    0.4,
                    -i,
                    "-->  " + str(value.replace('"', "")) + " " + data[i][2],
                    ha="left",
                    va="center",
                )
            else:
                if mfile_data.data[data[i][0]].exists:
                    dat = mfile_data.data[data[i][0]].get_scan(scan)
                    if isinstance(dat, str):
                        value = dat
                    else:
                        value = "{:.4g}".format(
                            mfile_data.data[data[i][0]].get_scan(scan)
                        )
                    if "alpha" in data[i][0]:
                        value = str(float(value) + 1.0)
                    axis.text(
                        eqpos,
                        -i,
                        "= " + value + " " + data[i][2],
                        ha="left",
                        va="center",
                    )
                else:
                    mfile_data.data[data[i][0]].get_scan(-1)
                    axis.text(
                        eqpos, -i, "=" + "ERROR! Var missing", ha="left", va="center"
                    )
        else:
            dat = data[i][0]
            if isinstance(dat, str):
                value = dat
            else:
                value = "{:.4g}".format(data[i][0])
            axis.text(
                eqpos, -i, "= " + value + " " + data[i][2], ha="left", va="center"
            )


def plot_geometry_info(axis, mfile_data, scan):
    """Function to plot geometry info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """
    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.text(-0.05, 1, "Geometry:", ha="left", va="center")
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    in_blanket_thk = mfile_data.data["shldith"].get_scan(scan) + mfile_data.data[
        "blnkith"
    ].get_scan(scan)
    out_blanket_thk = mfile_data.data["shldoth"].get_scan(scan) + mfile_data.data[
        "blnkoth"
    ].get_scan(scan)

    data = [
        ("rmajor", "$R_0$", "m"),
        ("rminor", "a", "m"),
        ("aspect", "A", ""),
        ("kappa95", r"$\kappa_{95}$", ""),
        ("triang95", r"$\delta_{95}$", ""),
        ("sarea", "Surface area", "m$^2$"),
        ("vol", "Plasma volume", "m$^3$"),
        ("n_tf", "No. of TF coils", ""),
        (in_blanket_thk, "i/b blkt/shld", "m"),
        (out_blanket_thk, "o/b blkt/shld", "m"),
        ("powfmw", "Fusion power", "MW"),
        ("", "", ""),
        ("#User Info", "", ""),
        ("!" + str(mfile_data.data["procver"].get_scan(scan)), "PROCESS Version", ""),
        ("!" + mfile_data.data["date"].get_scan(scan), "Date:", ""),
        ("!" + mfile_data.data["time"].get_scan(scan), "Time:", ""),
        ("!" + mfile_data.data["username"].get_scan(scan), "User:", ""),
        (
            "!"
            + proc_dict.DICT_OPTIMISATION_VARS[
                abs(int(mfile_data.data["minmax"].get_scan(scan)))
            ],
            "Optimising:",
            "",
        ),
    ]

    plot_info(axis, data, mfile_data, scan)


def plot_physics_info(axis, mfile_data, scan):
    """Function to plot geometry info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """
    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.text(-0.05, 1, "Physics:", ha="left", va="center")
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    nong = mfile_data.data["dnla"].get_scan(scan) / mfile_data.data[
        "dlimit(7)"
    ].get_scan(scan)

    dnz = mfile_data.data["dnz"].get_scan(scan) / mfile_data.data["dene"].get_scan(scan)

    data = [
        ("plascur/1d6", "$I_p$", "MA"),
        ("bt", "Vacuum $B_T$ as $R_0$", "T"),
        ("q", "$q_{edge}$", ""),
        ("normalised thermal beta", r"$\beta_N$, thermal", "% m T MA$^{-1}$"),
        ("normalised total beta", r"$\beta_N$, total", "% m T MA$^{-1}$"),
        ("thermal poloidal beta", r"$\beta_P$, thermal", ""),
        ("betap", r"$\beta_P$, total", ""),
        ("te", r"$< t_e >$", "keV"),
        ("dene", r"$< n_e >$", "m$^{-3}$"),
        (nong, r"$< n_{\mathrm{e,line}} >/n_G$", ""),
        ("alphat", r"$T_{e0}/ < T_e >$", ""),
        ("alphan", r"$n_{e0}/ < n_{\mathrm{e, vol}} >$", ""),
        ("zeff", r"$Z_{\mathrm{eff}}$", ""),
        ("zeffso", r"$Z_{\mathrm{eff, SoL}}$", ""),
        (dnz, r"$n_Z/ < n_{\mathrm{e, vol}} >$", ""),
        ("taueff", r"$\tau_e$", "s"),
        ("hfact", "H-factor", ""),
        ("tauelaw", "Scaling law", ""),
    ]

    plot_info(axis, data, mfile_data, scan)


def plot_magnetics_info(axis, mfile_data, scan):
    """Function to plot magnet info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """

    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.text(-0.05, 1, "Coil currents etc:", ha="left", va="center")
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    # Number of coils (1 is OH coil)
    number_of_coils = 0
    for item in mfile_data.data.keys():
        if "rpf(" in item:
            number_of_coils += 1

    pf_info = []
    for i in range(1, number_of_coils):
        if i % 2 != 0:
            pf_info.append(
                (
                    mfile_data.data["ric({:02})".format(i)].get_scan(scan),
                    "PF {}".format(i),
                )
            )

    tburn = mfile_data.data["tburn"].get_scan(scan) / 3600.0
    tftype = proc_dict.DICT_TF_TYPE[mfile_data.data["i_tf_sc_mat"].get_scan(scan)]
    vssoft = mfile_data.data["vsres"].get_scan(scan) + mfile_data.data[
        "vsind"
    ].get_scan(scan)

    data = [
        (pf_info[0][0], pf_info[0][1], "MA"),
        (pf_info[1][0], pf_info[1][1], "MA"),
        (pf_info[2][0], pf_info[2][1], "MA"),
        (vssoft, "Startup flux swing", "Wb"),
        ("vstot", "Available flux swing", "Wb"),
        (tburn, "Burn time", "hrs"),
        ("", "", ""),
        ("#TF coil type is {}".format(tftype), "", ""),
        ("bmaxtf", "Peak field at conductor", "T"),
        ("iooic", r"I/I$_{\mathrm{crit}}$", ""),
        ("tmarg", "Temperature margin", "K"),
        ("sig_tf_case", "TF case maximum shear stress (Tresca criterion)", "Pa"),
        ("sig_tf_wp", "TF conduit maximum shear stress (Tresca criterion)", "Pa"),
        (
            "sig_tf_case_max",
            "Allowable maximum shear stress in the TF case (Tresca criterion)",
            "Pa",
        ),
        (
            "sig_tf_wp_max",
            "Allowable maximum shear stress in the TF conduit (Tresca criterion)",
            "Pa",
        ),
        ("", "", ""),
        ("#Costs", "", ""),
        ("coe", "Cost of electricity", r"\$MWh"),
        ("concost", "Constructed cost", r"M\$"),
        ("capcost", "Total capex", r"M\$"),
    ]

    plot_info(axis, data, mfile_data, scan)


def plot_power_info(axis, mfile_data, scan):
    """Function to plot power info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """

    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.text(-0.05, 1, "Power flows/economics:", ha="left", va="center")
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    dnla = mfile_data.data["dnla"].get_scan(scan) / 1.0e20
    bt = mfile_data.data["bt"].get_scan(scan)
    surf = mfile_data.data["sarea"].get_scan(scan)
    pthresh = 0.0488 * dnla**0.717 * bt**0.803 * surf**0.941 * 0.8
    err = (
        0.057**2
        + (0.035 * sp.log(dnla)) ** 2
        + (0.032 * sp.log(bt)) ** 2
        + (0.019 * sp.log(surf)) ** 2
    )
    err = np.sqrt(err) * pthresh

    gross_eff = 100.0 * (
        mfile_data.data["pgrossmw"].get_scan(scan)
        / mfile_data.data["pthermmw"].get_scan(scan)
    )

    net_eff = 100.0 * (
        (
            mfile_data.data["pgrossmw"].get_scan(scan)
            - mfile_data.data["htpmw"].get_scan(scan)
        )
        / (
            mfile_data.data["pthermmw"].get_scan(scan)
            - mfile_data.data["htpmw"].get_scan(scan)
        )
    )

    plant_eff = 100.0 * (
        mfile_data.data["pnetelmw"].get_scan(scan)
        / mfile_data.data["powfmw"].get_scan(scan)
    )

    data = [
        ("wallmw", "Av. neutron wall load", "MW m$^{-2}$"),
        ("pinnerzoneradmw", "inner zone radiation", "MW"),
        ("psyncpv*vol", "Synchrotron radiation", "MW"),
        ("pouterzoneradmw", "outer zone radiation", "MW"),
        ("pnucblkt", "Nuclear heating in blanket", "MW"),
        ("pnucshld", "Nuclear heating in shield", "MW"),
        ("pdivt", "Psep / Pdiv", "MW"),
        (pthresh, "H-mode threshold (M=2.5)", r"$\pm${:.3f} MW".format(err)),
        ("fwbllife", "FW/Blanket life", "years"),
        ("divlife", "Divertor life", "years"),
        ("pthermmw", "Thermal Power", "MW"),
        (gross_eff, "Gross cycle efficiency", "%"),
        (net_eff, "Net cycle efficiency", "%"),
        ("pgrossmw", "Gross electric power", "MW"),
        ("pnetelmw", "Net electric power", "MW"),
        (
            plant_eff,
            "Plant efficiency " + r"$\frac{P_{\mathrm{e,net}}}{P_{\mathrm{fus}}}$",
            "%",
        ),
    ]

    plot_info(axis, data, mfile_data, scan)


def plot_current_drive_info(axis, mfile_data, scan):
    """Function to plot current drive info

    Arguments:
        axis --> axis object to plot to
        mfile_data --> MFILE.DAT object
        scan --> scan number to use

    """

    xmin = 0
    xmax = 1
    ymin = -16
    ymax = 1

    axis.text(-0.05, 1, "Neutral Beam Current Drive:", ha="left", va="center")
    axis.set_ylim([ymin, ymax])
    axis.set_xlim([xmin, xmax])
    axis.set_axis_off()
    axis.set_autoscaley_on(False)
    axis.set_autoscalex_on(False)

    pinjie = mfile_data.data["pinjmw"].get_scan(scan)

    pdivt = mfile_data.data["pdivt"].get_scan(scan)
    pdivr = pdivt / mfile_data.data["rmajor"].get_scan(scan)

    pdivnr = (
        10.0e20
        * mfile_data.data["pdivt"].get_scan(scan)
        / (
            mfile_data.data["rmajor"].get_scan(scan)
            * mfile_data.data["dene"].get_scan(scan)
        )
    )

    dnla = mfile_data.data["dnla"].get_scan(scan) / 1.0e20
    bt = mfile_data.data["bt"].get_scan(scan)
    surf = mfile_data.data["sarea"].get_scan(scan)
    pthresh = 0.0488 * dnla**0.717 * bt**0.803 * surf**0.941 * 0.8
    flh = pdivt / pthresh

    powerht = mfile_data.data["powerht"].get_scan(scan)
    psync = mfile_data.data["psyncpv*vol"].get_scan(scan)
    pbrem = mfile_data.data["pinnerzoneradmw"].get_scan(scan)
    hfact = mfile_data.data["hfact"].get_scan(scan)
    hstar = hfact * (powerht / (powerht + psync + pbrem)) ** 0.31

    data = [
        (pinjie, "SS auxiliary power", "MW"),
        ("pheat", "Power for heating only", "MW"),
        ("bootipf", "Bootstrap fraction", ""),
        ("faccd", "Auxiliary fraction", ""),
        ("facoh", "Ohmic fraction", ""),
        ("gamnb", "NB gamma", "$10^{20}$ A W$^{-1}$ m$^{-2}$"),
        ("enbeam", "NB energy", "keV"),
        ("powerht", "Assumed heating power", "MW"),
        (pdivr, r"$\frac{P_{\mathrm{div}}}{R_{0}}$", "MW m$^{-1}$"),
        (
            pdivnr,
            r"$\frac{P_{\mathrm{div}}}{<n> R_{0}}$",
            r"$\times 10^{-20}$ MW m$^{2}$",
        ),
        (flh, r"$\frac{P_{\mathrm{div}}}{P_{\mathrm{LH}}}$", ""),
        (hstar, "H* (non-rad. corr.)", ""),
    ]

    plot_info(axis, data, mfile_data, scan)


# def plot_video_info(axis, data, mfile_data, scan):
#     """Function to plot data to accompany video in written form on a
#     matplotlib plot.
#     """
#     eqpos = 0.7
#     for i in range(len(data)):
#         axis.text(0, -i, data[i][1], ha='left', va='center')
#         if isinstance(data[i][0], str):
#             if data[i][0] == "":
#                 axis.text(eqpos, -i, "\n",
#                           ha='left', va='center')
#             elif data[i][0][0] == "#":
#                 axis.text(-0.05, -i, "%s\n" % data[i][0][1:],
#                           ha='left', va='center')
#             elif data[i][0][0] == "!":
#                 value = data[i][0][1:]
#                 axis.text(0.4, -i, "-->  " + str(value.replace('"', '')) +
#                           " " + data[i][2], ha='left', va='center')
#             else:
#                 dat = mfile_data.data[data[i][0]].get_scan(scan)
#                 if isinstance(dat, str):
#                     value = dat
#                 else:
#                     value = "%.4g" % mfile_data.data[data[i][0]].get_scan(scan)
#                 if "alpha" in data[i][0]:
#                     value = str(float(value) + 1.0)
#                 axis.text(eqpos, -i, '= ' + value + ' ' + data[i][2],
#                           ha='left', va='center')
#         else:
#             dat = data[i][0]
#             if isinstance(dat, str):
#                 value = dat
#             else:
#                 value = "%.4g" % data[i][0]
#             axis.text(eqpos, -i, '= ' + value + ' ' + data[i][2],
#                       ha='left', va='center')
#
#
# def plot_animation_info(axis, mfile_data, scan, animation_info=ANIMATION_INFO):
#     """Function to plot animation information. Where animation information
#     is defined by the user. There is a default animation_info that can be
#     passed as an argument.
#     """
#
#     xmin = 0
#     xmax = 1
#     ymin = -16
#     ymax = 1
#
#     axis.text(-0.05, 1, 'Plant Information:', ha='left', va='center')
#     axis.set_ylim([ymin, ymax])
#     axis.set_xlim([xmin, xmax])
#     axis.set_axis_off()
#     axis.set_autoscaley_on(False)
#     axis.set_autoscalex_on(False)
#
#     plot_video_info(axis, animation_info, mfile_data, scan)
#
#
# def plot_animation_graph(axis, mfile_data, scan, x, y, xlabel, ylabel):
#     """Function to plot animation information. Where animation information
#     is defined by the user. There is a default animation_info that can be
#     passed as an argument.
#     """
#
#     #axis.set_ylim([ymin, ymax])
#     #axis.set_xlim([xmin, xmax])
#     #axis.set_axis_off()
#     axis.set_autoscaley_on(False)
#     axis.set_autoscalex_on(False)
#     axis.set_xlabel(xlabel)
#     axis.set_ylabel(ylabel)
#     axis.plot(x, y)
#     a = random.uniform(0, 1)
#     b = random.uniform(0, 1)
#     axis.plot([a, a],
#               [b, b],
#               color="red", marker="o", alpha=0.5, markersize=20)
