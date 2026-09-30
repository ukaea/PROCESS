"""Radial build lookup utilities used by geometry and magnet plots."""

from __future__ import annotations

from typing import TYPE_CHECKING

from process.core.io.plot.summary.constants import RADIAL_BUILD

if TYPE_CHECKING:
    from process.core.io.mfile import MFile


def cumulative_radial_build(section, mfile: MFile, scan: int):
    """Function for calculating the cumulative radial build up to and
    including the given section.

    Parameters
    ----------
    section :
        section of the radial build to go up to
    mfile :
        MFILE data object
    scan :
        scan number to use

    Returns
    -------
    :
        cumulative_build:cumulative radial build up to section given
    """
    complete = False
    cumulative_build = 0
    for item in RADIAL_BUILD:
        if item in {"rminori", "rminoro"}:
            cumulative_build += mfile.get("rminor", scan=scan)
        elif item in {"vvblgapi", "vvblgapo"}:
            cumulative_build += mfile.get("dr_shld_blkt_gap", scan=scan)
        elif "dr_vv_inboard" in item:
            cumulative_build += mfile.get("dr_vv_inboard", scan=scan)
        elif "dr_vv_outboard" in item:
            cumulative_build += mfile.get("dr_vv_outboard", scan=scan)
        else:
            cumulative_build += mfile.get(item, scan=scan)
        if item == section:
            complete = True
            break

    if complete is False:
        print("radial build parameter ", section, " not found")
    return cumulative_build


def cumulative_radial_build2(section, mfile: MFile, scan: int):
    """Function for calculating the cumulative radial build up to and
    including the given section.

    Parameters
    ----------
    section :
        section of the radial build to go up to
    mfile :
        MFILE data object
    scan :
        scan number to use

    Returns
    -------
    :
        cumulative_build --> cumulative radial build up to and including
        section given
        previous         --> cumulative radial build up to section given
    """
    cumulative_build = 0
    build = 0
    for item in RADIAL_BUILD:
        if item in {"rminori", "rminoro"}:
            build = mfile.get("rminor", scan=scan)
        elif item in {"vvblgapi", "vvblgapo"}:
            build = mfile.get("dr_shld_blkt_gap", scan=scan)
        elif "dr_vv_inboard" in item:
            build = mfile.get("dr_vv_inboard", scan=scan)
        elif "dr_vv_outboard" in item:
            build = mfile.get("dr_vv_outboard", scan=scan)
        else:
            build = mfile.get(item, scan=scan)
        cumulative_build += build
        if item == section:
            break
    previous = cumulative_build - build
    return (cumulative_build, previous)


__all__ = ["cumulative_radial_build", "cumulative_radial_build2"]
