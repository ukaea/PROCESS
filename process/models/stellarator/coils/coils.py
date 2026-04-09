import logging

import numpy as np

import process.models.superconductors as superconductors
from process.core.exceptions import ProcessValueError
from process.data_structure import (
    stellarator_configuration,
)

logger = logging.getLogger(__name__)


def jcrit_from_material(
    b_max,
    t_helium,
    i_tf_sc_mat,
    b_crit_upper_nbti,
    b_crit_sc,
    f_a_tf_turn_cable_copper,
    f_hts,
    t_crit_nbti,
    t_crit_sc,
    f_a_tf_turn_cable_space_extra_void,
    j_wp,
):
    strain = -0.005  # for now a small value
    f_he = f_a_tf_turn_cable_space_extra_void  # this is helium fraction in the superconductor (set it to the fixed global variable here)

    f_tf_conductor_copper = (
        f_a_tf_turn_cable_copper  # fcutfsu is a global variable. Is the copper fraction
    )
    # of a cable conductor.

    if i_tf_sc_mat == 1:  # ITER Nb3Sn critical surface parameterization
        bc20m = 32.97  # these are values taken from sctfcoil.f90
        tc0m = 16.06

        #  j_crit_sc returned by itersc is the critical current density in the
        #  superconductor - not the whole strand, which contains copper
        if b_max > bc20m:
            j_crit_sc = 1.0e-9  # Set to a small nonzero value
        else:
            (
                j_crit_sc,
                _bcrit,
                _tcrit,
            ) = superconductors.itersc(t_helium, b_max, strain, bc20m, tc0m)

        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1.0 - f_tf_conductor_copper) * (1.0e0 - f_he)

        # This is needed right now. Can we change it later?
        j_crit_sc = max(1.0e-9, j_crit_sc)
        j_crit_cable = max(1.0e-9, j_crit_cable)

    elif i_tf_sc_mat == 2:
        # Bi-2212 high temperature superconductor parameterization
        #  Current density in a strand of Bi-2212 conductor
        #  N.B. jcrit returned by bi2212 is the critical current density
        #  in the strand, not just the superconducting portion.
        #  The parameterization for j_crit_cable assumes a particular strand
        #  composition that does not require a user-defined copper fraction,
        #  so this is irrelevant in this model

        jstrand = j_wp / (1 - f_he)
        #  jstrand = 0  # as far as I can tell this will always be 0
        #  because jwp was never set in fortran (so 0)

        j_crit_cable, tmarg = superconductors.bi2212(
            b_max, jstrand, t_helium, f_hts
        )  # bi2212 outputs j_crit_cable
        j_crit_sc = j_crit_cable / (1 - f_tf_conductor_copper)
        _tcrit = t_helium + tmarg
    elif i_tf_sc_mat == 3:  # NbTi data
        bc20m = 15.0
        tc0m = 9.3
        c0 = 1.0

        if b_max > bc20m:
            j_crit_sc = 1.0e-9  # Set to a small nonzero value
        else:
            j_crit_sc, _tcrit = superconductors.jcrit_nbti(
                t_helium,
                b_max,
                c0,
                bc20m,
                tc0m,
            )
            # I dont need tcrit here so dont use it.

        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)

        # This is needed right now. Can we change it later?
        j_crit_sc = max(1.0e-9, j_crit_sc)
        j_crit_cable = max(1.0e-9, j_crit_cable)
    elif i_tf_sc_mat == 4:  # As (1), but user-defined parameters
        bc20m = b_crit_sc
        tc0m = t_crit_sc
        j_crit_sc, _bcrit, _tcrit = superconductors.itersc(
            t_helium, b_max, strain, bc20m, tc0m
        )
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)
    elif i_tf_sc_mat == 5:  # WST Nb3Sn parameterisation
        bc20m = 32.97
        tc0m = 16.06

        #  j_crit_sc returned by itersc is the critical current density in the
        #  superconductor - not the whole strand, which contains copper

        j_crit_sc, _bcrit, _tcrit = superconductors.western_superconducting_nb3sn(
            t_helium,
            b_max,
            strain,
            bc20m,
            tc0m,
        )
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)
    elif i_tf_sc_mat == 6:  # ! "REBCO" 2nd generation HTS superconductor in CrCo strand
        j_crit_sc, _validity = superconductors.jcrit_rebco(t_helium, b_max, 0)
        j_crit_sc = max(1.0e-9, j_crit_sc)
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)

    elif i_tf_sc_mat == 7:  # Durham Ginzburg-Landau Nb-Ti parameterisation
        bc20m = b_crit_upper_nbti
        tc0m = t_crit_nbti
        j_crit_sc, _bcrit, _tcrit = superconductors.gl_nbti(
            t_helium, b_max, strain, bc20m, tc0m
        )
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)
    elif i_tf_sc_mat == 8:
        bc20m = 429
        tc0m = 185
        j_crit_sc, _bcrit, _tcrit = superconductors.gl_rebco(
            t_helium, b_max, strain, bc20m, tc0m
        )
        # A0 calculated for tape cross section already
        # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
        j_crit_cable = j_crit_sc * (1 - f_tf_conductor_copper) * (1 - f_he)
    else:
        raise ProcessValueError(
            "Illegal value for i_pf_superconductor", i_tf_sc_mat=i_tf_sc_mat
        )

    return j_crit_sc * 1e-6


def intersect(x1, y1, x2, y2, xin):
    """Routine to find the x (abscissa) intersection point of two curves
    each defined by tabulated (x,y) values

    This routine estimates the x point (abscissa) at which two curves
    defined by tabulated (x,y) values intersect, using simple
    linear interpolation and the Newton-Raphson method.
    The routine will stop with an error message if no crossing point
    is found within the x ranges of the two curves.

    Parameters
    ----------
    x1 :
        x values for first curve
    y1 :
        y values for first curve
    x2 :
        x values for first curve
    y2 :
        y values for first curve
    xin :
        initial guess for intersection point

    Returns
    -------
    :
        x value at point of intersection on exit

    """
    x = xin
    n1 = len(x1)
    n2 = len(x2)

    xmin = max(np.amin(x1), np.amin(x2))
    xmax = min(np.max(x1), np.amax(x2))

    if xmin >= xmax:
        logger.error(
            f"X ranges not overlapping. {np.amin(x1)=} {np.amin(x2)=} "
            f"{np.amax(x1)=} {np.amax(x2)=}"
        )

    #  Ensure input guess for x is within this range
    if x < xmin:
        x = xmin
    elif x > xmax:
        x = xmax

    #  Find overall y range, and set tolerance
    #  in final difference in y values

    ymin = min(np.amin(y1), np.amin(y2))
    ymax = max(np.max(y1), np.max(y2))

    epsy = 1.0e-6 * (ymax - ymin)

    #  Finite difference dx

    dx = 0.01e0 / max(n1, n2) * (xmax - xmin)

    for _i in range(100):
        #  Find difference in y values at x

        y01 = np.interp(x, x1, y1)
        y02 = np.interp(x, x2, y2)
        y = y01 - y02

        if abs(y) < epsy:
            break

        #  Find difference in y values at x+dx

        y01 = np.interp(x + dx, x1, y1)
        y02 = np.interp(x + dx, x2, y2)
        yright = y01 - y02

        #  Find difference in y values at x-dx

        y01 = np.interp(x - dx, x1, y1)
        y02 = np.interp(x - dx, x2, y2)
        yleft = y01 - y02

        #  Adjust x using Newton-Raphson method

        x = x - 2.0e0 * dx * y / (yright - yleft)

        if x < xmin:
            logger.error(
                f"X has dropped below Xmin; X={x} has been set equal to Xmin={xmin}"
            )
            x = xmin
            break

        if x > xmax:
            logger.error(
                f"X has risen above Xmax; X={x} has been set equal to Xmax={xmin}"
            )
            x = xmax
            break
    else:
        logger.error("Convergence too slow; X may be wrong...")

    return x


def bmax_from_awp(wp_width_radial, current, n_tf_coils, r_coil_major, r_coil_minor):
    """Returns a fitted function for bmax for stellarators

    Returns a fitted function for bmax in dependence
    of the winding pack. The stellarator type config
    is taken from the parent scope.

    Parameters
    ----------
    wp_width_radial :

    current :

    n_tf_coils :

    r_coil_major :

    r_coil_minor :


    """

    return (
        2e-1  # this is mu x 1e6, to use current in MA
        * current
        * n_tf_coils
        / (r_coil_major - r_coil_minor)
        * (
            stellarator_configuration.stella_config_a1
            + stellarator_configuration.stella_config_a2 * r_coil_major / wp_width_radial
        )
    )
