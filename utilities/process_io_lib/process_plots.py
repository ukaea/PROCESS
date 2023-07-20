"""

  PROCESS plotting library

  James Morris
  27/07/17
  CCFE

"""
# flake8: noqa

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

color_sequence = [
    "#1f77b4",
    "#aec7e8",
    "#ff7f0e",
    "#ffbb78",
    "#2ca02c",
    "#98df8a",
    "#d62728",
    "#ff9896",
    "#9467bd",
    "#c5b0d5",
    "#8c564b",
    "#c49c94",
    "#e377c2",
    "#f7b6d2",
    "#7f7f7f",
    "#c7c7c7",
    "#bcbd22",
    "#dbdb8d",
    "#17becf",
    "#9edae5",
]

plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "Ubuntu"
plt.rcParams["font.monospace"] = "Ubuntu Mono"
plt.rcParams["font.size"] = 10
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["xtick.labelsize"] = 9
plt.rcParams["ytick.labelsize"] = 9
plt.rcParams["legend.fontsize"] = 10

# PROCESS libraries import
try:
    import process.io.mfile as mfile
except ImportError:
    print("The PROCESS python libraries are not in your path. Please check PYTHONPATH")
    exit()


def align_yaxis(ax1, v1, ax2, v2):
    """adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1"""
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    adjust_yaxis(ax2, (y1 - y2) / 2, v2)
    adjust_yaxis(ax1, (y2 - y1) / 2, v1)


def adjust_yaxis(ax, ydif, v):
    """shift axis ax by ydiff, maintaining point v at the same location"""
    inv = ax.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, ydif))
    miny, maxy = ax.get_ylim()
    miny, maxy = miny - v, maxy - v
    if -miny > maxy or (-miny == maxy and dy > 0):
        nminy = miny
        nmaxy = miny * (maxy + dy) / (miny + dy)
    else:
        nmaxy = maxy
        nminy = maxy * (miny + dy) / (maxy + dy)
    ax.set_ylim(nminy + v, nmaxy + v)


def plot_pulse_timings(
    mfile_obj, save=False, show=False, return_fig=False, save_path=""
):
    """
    Plot the PROCESS pulse timings plot

    Will plot Ip, CS and PF currents over time

    mfile_obj  : MFILE.DAT object
    save       : save to "process_times.png"
    show       : show to screen on running
    return_fig : return fig object

    """

    # setup figure
    fig = plt.figure(dpi=100)
    fig.subplots_adjust(hspace=0.1)

    # # Setup plot axes
    gs = gridspec.GridSpec(4, 1, height_ratios=[1, 4, 4, 4])
    axis_0 = fig.add_subplot(gs[0])
    axis_1 = fig.add_subplot(gs[1])
    axis_2 = fig.add_subplot(gs[2])
    axis_3 = fig.add_subplot(gs[3])

    # Number, n, of PF circuits including CS coil (n-1) and plasma n
    number_of_pf_circuits = int(mfile_obj.data["ncirt"].get_scan(-1))

    # Number of time points
    n_time = 6

    # Time array
    time_list = ["tramp", "tohs", "theat", "tburn", "tqnch", "tdwell"]
    times_array = np.array([mfile_obj.data[val].get_scan(-1) for val in time_list])
    time = np.array([sum(times_array[:i]) for i in range(len(times_array))])

    # Current arrays
    current_arrays = np.zeros(shape=(number_of_pf_circuits, n_time))

    # PF coils
    for i in range(number_of_pf_circuits - 2):
        for j in range(n_time):
            variable_name = "pfc{0:02}t{1:02}".format(i + 1, j + 1)
            current_arrays[i][j] = mfile_obj.data[variable_name].get_scan(-1) / 1e6

    # Central solenoid
    for m in range(n_time):
        variable_name = "cst{0:02}".format(m + 1)
        current_arrays[number_of_pf_circuits - 2][m] = (
            mfile_obj.data[variable_name].get_scan(-1) / 1e6
        )

    # Plasma
    for p in range(n_time):
        variable_name = "plasmat{0:02}".format(p + 1)
        current_arrays[number_of_pf_circuits - 1][p] = (
            mfile_obj.data[variable_name].get_scan(-1) / 1e6
        )

    # Plotting

    # Axis 1
    axis_1.axhline(y=0, color="k")
    axis_1.set_ylabel("Current [MA]")
    axis_1.legend(loc="upper left", bbox_to_anchor=(1, 1.05))

    # limits
    ymin_1 = 1.2 * min(current_arrays[-2])
    ymax_1 = 1.2 * max(current_arrays[-2])
    axis_1.set_ylim([ymin_1, ymax_1])

    for ti in range(len(time)):
        axis_1.axvline(x=time[ti], color="k", linestyle="--", linewidth=0.5, alpha=0.5)

    # Axis 2
    axis_2.axhline(y=0, color="k")
    axis_2.set_ylabel("Current [MA]")
    axis_2.set_ylim([0, 25])
    axis_2.legend(loc="upper left", bbox_to_anchor=(1, 1.05))

    # limits
    ymin_2 = 1.2 * min(current_arrays[-1])
    ymax_2 = 1.2 * max(current_arrays[-1])
    axis_2.set_ylim([ymin_2, ymax_2])

    for ti in range(len(time)):
        axis_2.axvline(x=time[ti], color="k", linestyle="--", linewidth=0.5, alpha=0.5)

    # Axis 3
    ymin_3 = 0
    ymax_3 = 0
    for s in range(number_of_pf_circuits - 2):
        axis_3.plot(
            time,
            current_arrays[s],
            label="PF circuit {0}".format(s + 1),
            color=color_sequence[s],
        )

        # limits
        if 1.2 * min(current_arrays[s]) < ymin_3:
            ymin_3 = 1.2 * min(current_arrays[s])
        if 1.2 * max(current_arrays[s]) > ymax_3:
            ymax_3 = 1.2 * max(current_arrays[s])

    axis_3.axhline(y=0, color="k")
    axis_3.legend(loc="upper left", bbox_to_anchor=(1, 1.05))
    axis_3.set_xlabel("Time [s]")
    axis_3.set_ylabel("Current [MA]")

    # limits
    axis_3.set_ylim([ymin_3, ymax_3])

    for ti in range(len(time)):
        axis_3.axvline(x=time[ti], color="k", linestyle="--", linewidth=0.5, alpha=0.5)

    # Time bar plot
    y_pos = 0
    left = 0
    patch_handles = list()

    for i in range(len(time_list)):
        time_length = mfile_obj.data[time_list[i]].get_scan(-1)
        if time_length != 0:
            patch_handles.append(
                axis_0.barh(
                    y_pos,
                    time_length,
                    color=color_sequence[i],
                    align="center",
                    left=left,
                    height=5,
                    label=time_list[i],
                )
            )

            # accumulate the left-hand offsets
            left += time_length

    axis_0.legend(
        bbox_to_anchor=(0, 1.1, 1, 0.2),
        loc="lower left",
        mode="expand",
        borderaxespad=0,
        ncol=15,
        labelspacing=2,
    )
    axis_0.axes.get_yaxis().set_visible(False)
    axis_0.tick_params(axis="both", which="both", length=0)

    xticklabels = (
        axis_0.get_xticklabels() + axis_1.get_xticklabels() + axis_2.get_xticklabels()
    )
    plt.setp(xticklabels, visible=False)

    # Options
    if save:
        plt.savefig(
            "{0}process_times.png".format(save_path), dpi=300, bbox_inches="tight"
        )

    if show:

        plt.show()

    if return_fig:
        return fig

    return


def radial_bar_plot(
    build_params, mfile_obj, save=False, show=False, return_fig=False, save_path=""
):
    """
    Plot the PROCESS radial build bar plot

    Will plot radial build thicknesses

    build_params : list of radial build entry names
    mfile_obj    : MFILE.DAT object
    save         : save to "process_times.png"
    show         : show to screen on running
    return_fig   : return fig object
    """

    fig = plt.figure(figsize=(18, 2))
    ax = fig.add_subplot(111)

    y_pos = 0
    left = 0
    patch_handles = list()

    for i in range(len(build_params)):
        item_thickness = mfile_obj.data[build_params[i]].get_scan(-1)
        if item_thickness != 0:
            patch_handles.append(
                ax.barh(
                    y_pos,
                    item_thickness,
                    color=color_sequence[i],
                    align="center",
                    left=left,
                    height=5,
                    label=build_params[i],
                )
            )

            # accumulate the left-hand offsets
            left += item_thickness

    ax.legend(
        bbox_to_anchor=(0, 1.02, 1, 0.2),
        loc="lower left",
        mode="expand",
        borderaxespad=0,
        ncol=15,
        labelspacing=2,
    )
    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlabel("Radial Build [m]")
    last_tick = left
    minor_ticks = np.arange(0, last_tick, 0.1)
    ax.set_xticks(minor_ticks, minor=True)

    # Options
    if save:
        plt.savefig(
            "{0}process_radial_plot.png".format(save_path), dpi=300, bbox_inches="tight"
        )

    if show:
        plt.show()

    if return_fig:
        return fig

    return


def plasma_profiles_plot(
    mfile_obj, save=False, show=False, return_fig=False, save_path=""
):
    """
    Plot the PROCESS plasma profiles

    mfile_obj    : MFILE.DAT object
    save         : save to "process_times.png"
    show         : show to screen on running
    return_fig   : return fig object
    """

    fig = plt.figure(dpi=100)
    fig.subplots_adjust(wspace=0.5)
    axis_1 = fig.add_subplot(131, aspect=0.05)
    axis_2 = fig.add_subplot(132, aspect=1 / 35)
    axis_3 = fig.add_subplot(133, aspect=1 / 10)

    # Parameters
    ipedestal = mfile_obj.data["ipedestal"].get_scan(-1)
    rhopedn = mfile_obj.data["rhopedn"].get_scan(-1)
    rhopedt = mfile_obj.data["rhopedt"].get_scan(-1)
    ne0 = mfile_obj.data["ne0"].get_scan(-1)
    neped = mfile_obj.data["neped"].get_scan(-1)
    teped = mfile_obj.data["teped"].get_scan(-1)
    alphan = mfile_obj.data["alphan"].get_scan(-1)
    alphat = mfile_obj.data["alphat"].get_scan(-1)
    tesep = mfile_obj.data["tesep"].get_scan(-1)
    te0 = mfile_obj.data["te0"].get_scan(-1)
    tbeta = mfile_obj.data["tbeta"].get_scan(-1)
    q0 = mfile_obj.data["q0"].get_scan(-1)
    q95 = mfile_obj.data["q95"].get_scan(-1)

    # density profile plot

    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 20
    axis_1.set_ylim([ymin, ymax])
    axis_1.set_xlim([xmin, xmax])
    axis_1.set_autoscaley_on(False)
    axis_1.set_xlabel("r/a")
    axis_1.set_ylabel("ne / 1e19 m-3")
    axis_1.set_title("Density profile")

    if mfile_obj.data["ipedestal"] == 1:
        rhocore1 = np.linspace(0, 0.95 * rhopedn)
        rhocore2 = np.linspace(0.95 * rhopedn, rhopedn)
        rhocore = np.append(rhocore1, rhocore2)
        ncore = neped + (ne0 - neped) * (1 - rhocore**2 / rhopedn**2) ** alphan

        rhosep = np.linspace(rhopedn, 1)
        nsep = nesep + (neped - nesep) * (1 - rhosep) / (1 - min(0.9999, rhopedn))

        rho = np.append(rhocore, rhosep)
        ne = np.append(ncore, nsep)
    else:
        rho1 = np.linspace(0, 0.95)
        rho2 = np.linspace(0.95, 1)
        rho = np.append(rho1, rho2)
        ne = ne0 * (1 - rho**2) ** alphan
    ne = ne / 1e19
    axis_1.plot(rho, ne)

    # temperature profile plot

    xmin = 0
    xmax = 1
    ymin = 0
    ymax = 35
    axis_2.set_ylim([ymin, ymax])
    axis_2.set_xlim([xmin, xmax])
    axis_2.set_autoscaley_on(False)
    axis_2.set_xlabel("r/a")
    axis_2.set_ylabel("Te / keV")
    axis_2.set_title("Temperature profile")

    if ipedestal == 1:
        rhocore1 = np.linspace(0, 0.9 * rhopedn)
        rhocore2 = np.linspace(0.9 * rhopedn, rhopedn)
        rhocore = np.append(rhocore1, rhocore2)
        tcore = teped + (te0 - teped) * (1 - (rhocore / rhopedn) ** tbeta) ** alphat

        rhosep = np.linspace(rhopedn, 1)
        tsep = tesep + (teped - tesep) * (1 - rhosep) / (1 - min(0.9999, rhopedt))

        rho = np.append(rhocore, rhosep)
        te = np.append(tcore, tsep)
    else:
        rho = np.linspace(0, 1)
        te = te0 * (1 - rho**2) ** alphat
    axis_2.plot(rho, te)

    # q profile plot

    XMIN = 0
    XMAX = 1
    YMIN = 0
    YMAX = 10
    axis_3.set_ylim([YMIN, YMAX])
    axis_3.set_xlim([XMIN, XMAX])
    axis_3.set_autoscaley_on(False)
    axis_3.set_xlabel("r/a")
    axis_3.set_ylabel("q(r)")
    axis_3.set_title("q profile")

    rho = np.linspace(0, 1)
    q_r_nevin = q0 + (q95 - q0) * (rho + rho * rho + rho**3) / (3.0)
    q_r_sauter = q0 + (q95 - q0) * (rho * rho)

    axis_3.plot(rho, q_r_nevin, label="Nevin")
    axis_3.plot(rho, q_r_sauter, label="Sauter")
    axis_3.legend()

    # Options
    if save:
        plt.savefig(
            "{0}process_plasma_profiles.png".format(save_path),
            dpi=300,
            bbox_inches="tight",
        )

    if show:
        plt.show()

    if return_fig:
        return fig

    return


if __name__ == "__main__":
    print("Running in test mode")
    mf = mfile.MFile(filename="MFILE.DAT")
    plot_pulse_timings(mf, save=True, show=False, save_path="figures/")
    radial_build = [
        "bore",
        "ohcth",
        "precomp",
        "gapoh",
        "tfcth",
        "deltf",
        "thshield_ib",
        "gapds",
        "d_vv_in",
        "shldith",
        "vvblgap",
        "blnkith",
        "fwith",
        "scrapli",
        "rminor",
    ]
    radial_bar_plot(radial_build, mf, show=False, save=True)
    plasma_profiles_plot(mf, show=False, save=True)
