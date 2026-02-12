"""
Module to hold plotting functions, used in plot_proc.py, which are common to multiple reactor components
"""

import numpy as np


def dh_vertices(
    r0: float, a: float, triang: float, kap: float
) -> tuple[np.ndarray, np.ndarray]:
    """Returns the radial and vertical coordinates which, when plotted, plots half a thin D-section, centred on z = 0

    Parameters
    ----------
    r0:
        major radius of centre
    a:
        horizontal radius
    triang:
        plasma triangularity
    kap:
        plasma elongation

    Returns
    -------
    :
        tuple containing radial and vertical coordinates which, when plotted, plots a half thin D-section with a gap
    """
    angs = np.linspace(0, np.pi, 50, endpoint=True)
    rs = r0 + a * np.cos(angs + triang * np.sin(1.0 * angs))
    zs = kap * a * np.sin(angs)
    return rs, zs


def dhgap_vertices(
    inpt: float,
    outpt: float,
    inthk: float,
    outthk: float,
    toppt: float,
    topthk: float,
    triang: float,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Returns the radial and vertical coordinates which, when plotted, plots a half thick D-section with a gap

    Parameters
    ----------
    inpt:
        inner point
    outpt:
        outer point
    inthk:
        inner thickness
    outthk:
        outer thickness
    toppt:
        top point
    topthk:
        top thickness
    triang:
        plasma triangularity

    Returns
    -------
    :
        tuple containing radial and vertical coordinates which, when plotted, plots a half thick D-section with a gap
    """
    arc = np.pi / 4.0
    r01 = (inpt + outpt) / 2.0
    r02 = (inpt + inthk + outpt - outthk) / 2.0
    a1 = r01 - inpt
    a2 = r02 - inpt - inthk
    kap1 = toppt / a1
    kap2 = (toppt - topthk) / a2
    angs = np.linspace(0.0, (np.pi / 2.0) - arc / 2.0, 50, endpoint=True)
    rs1 = r01 + a1 * np.cos(angs + triang * np.sin(angs))
    zs1 = kap1 * a1 * np.sin(angs)
    rs2 = r02 + a2 * np.cos(angs + triang * np.sin(angs))
    zs2 = kap2 * a2 * np.sin(angs)
    angs = np.linspace(np.pi, np.pi + ((np.pi / 2.0) - arc), 50, endpoint=True)
    rs3 = r01 + a1 * np.cos(angs + triang * np.sin(angs))
    zs3 = kap1 * a1 * np.sin(angs)
    rs4 = r02 + a2 * np.cos(angs + triang * np.sin(angs))
    zs4 = kap2 * a2 * np.sin(angs)

    return rs1, rs2, rs3, rs4, zs1, zs2, zs3, zs4


def ellips_fill_vertices(
    a1: float = 0,
    a2: float = 0,
    b1: float = 0,
    b2: float = 0,
    x0: float = 0,
    y0: float = 0,
    ang1: float = 0,
    ang2: float = np.pi / 2,
) -> list[tuple[float, float]]:
    """Returns the vertices of a shape which, when filled, fills the space between two concentric ellipse sectors

    Parameters
    ----------
    a1:
        horizontal radius to be filled, defaults to 0
    a2:
        horizontal radius to be filled, defaults to 0
    b1:
        vertical radius to be filled, defaults to 0
    b2:
        vertical radius to be filled, defaults to 0
    x0:
        x coordinate of centre of ellipses, defaults to 0
    y0:
        y coordinate of centre of ellipses, defaults to 0
    ang1:
        polar angle at start, defaults to 0
    ang2:
        polar angle at end, defaults to np.pi/2

    Returns
    -------
    :
        list containing (R,Z) coordinates which, when plotted, fill space between ellipses
    """
    angs = np.linspace(ang1, ang2, endpoint=True)
    r1 = ((np.cos(angs) / a1) ** 2 + (np.sin(angs) / b1) ** 2) ** (-0.5)
    xs1 = r1 * np.cos(angs) + x0
    ys1 = r1 * np.sin(angs) + y0
    angs = np.linspace(ang2, ang1, endpoint=True)
    r2 = ((np.cos(angs) / a2) ** 2 + (np.sin(angs) / b2) ** 2) ** (-0.5)
    xs2 = r2 * np.cos(angs) + x0
    ys2 = r2 * np.sin(angs) + y0
    verts = list(zip(xs1, ys1, strict=False))
    verts.extend(list(zip(xs2, ys2, strict=False)))
    endpoint = verts[-1:]
    verts.extend(endpoint)

    return verts
