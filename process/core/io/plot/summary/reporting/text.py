"""Text-based report panel helpers."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import matplotlib.pyplot as plt

    from process.core.io.mfile import MFile


def plot_info(axis: plt.Axes, data, mfile: MFile, scan: int):
    """Function to plot data in written form on a matplotlib plot.

    Parameters
    ----------
    axis :
        axis object to plot to
    data :
        plot information
    mfile :
        MFILE
    scan :
        scan number to use
    """
    eqpos = 0.75
    for i in range(len(data)):
        colorflag = "black"
        if mfile.data[data[i][0]].exists:
            if mfile.data[data[i][0]].var_flag == "ITV":
                colorflag = "red"
            elif mfile.data[data[i][0]].var_flag == "OP":
                colorflag = "blue"
        axis.text(0, -i, data[i][1], color=colorflag, ha="left", va="center")
        if isinstance(data[i][0], str):
            if not data[i][0]:
                axis.text(eqpos, -i, "\n", ha="left", va="center")
            elif data[i][0][0] == "#":
                axis.text(
                    -0.05,
                    -i,
                    f"{data[i][0][1:]}\n",
                    ha="left",
                    va="center",
                )
            elif data[i][0][0] == "!":
                value = data[i][0][1:].replace('"', "")
                axis.text(
                    0.4,
                    -i,
                    f"-->  {value} {data[i][2]}",
                    ha="left",
                    va="center",
                )
            elif mfile.data[data[i][0]].exists:
                dat = mfile.get(data[i][0], scan=scan)
                if isinstance(dat, str):
                    value = dat
                else:
                    value = f"{mfile.get(data[i][0], scan=scan):.4g}"
                if "alpha" in data[i][0]:
                    value = str(float(value) + 1.0)
                    axis.text(
                        eqpos,
                        -i,
                        f"= {value} {data[i][2]}",
                        color=colorflag,
                        ha="left",
                        va="center",
                    )
            else:
                mfile.get(data[i][0], scan=-1)
                axis.text(
                    eqpos,
                    -i,
                    "= ERROR! Var missing",
                    color=colorflag,
                    ha="left",
                    va="center",
                )
        else:
            dat = data[i][0]
            value = dat if isinstance(dat, str) else f"{data[i][0]:.4g}"
            axis.text(
                eqpos,
                -i,
                f"= {value} {data[i][2]}",
                color=colorflag,
                ha="left",
                va="center",
            )


__all__ = ["plot_info"]
