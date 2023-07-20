"""Contains some methods from the now-removed kallenbach_module.f90 file.

Running this file will run the kallenbach_paper routine, which should probably
become a test at some point.
"""

import logging
from pathlib import Path
import numpy

from process.fortran import constants
from process.fortran import physics_variables
from process.fortran import div_kal_vars
from process.fortran import impurity_radiation_module
from process.fortran import divertor_ode

from process.fortran import read_and_get_atomic_data
from process.fortran import plot_radiation

from process.fortran import global_variables
from process.utilities import f2py_string_patch
from process.main import SingleRun

logger = logging.getLogger(__name__)
# Logging handler for console output
s_handler = logging.StreamHandler()
s_handler.setLevel(logging.INFO)
logger.addHandler(s_handler)

CURRENT_PATH = Path(__file__).resolve().parent


def kallenbach_paper():
    # TODO: this should probably be turned into a pytest at some point

    SingleRun.init_module_vars()
    SingleRun.initialise()

    rmajor = 8.0e0
    rminor = 2.75e0
    bt = 4.00972e0 * (rmajor + rminor) / rmajor
    plascur = 1.33542e0 * (2.0e0 * numpy.pi * rminor) / constants.rmu0
    q = 3.0e0
    t_target = 2.3e0
    q_target_total = 4.175e6
    target_angle = 30.0e0
    b_pol = 0.956e0
    physics_variables.tesep = 0.298
    div_kal_vars.target_spread = 7.0e-3
    div_kal_vars.netau_sol = 0.5e0
    div_kal_vars.lambda_q_omp = 0.002e0

    # MDK Issue #494.
    # Invert the equation for lcon to ensure lcon=100 for this test.
    # lcon = lcon_factor * 0.395d0*pi*q*rmajor/lambda_omp**0.196
    div_kal_vars.lcon_factor = 100.0e0 / (
        0.395e0 * numpy.pi * 3.0e0 * 8.0e0 / 0.002e0**0.196
    )

    global_variables.fileprefix = f2py_string_patch.string_to_f2py_compatible(
        global_variables.fileprefix, str(CURRENT_PATH)
    )

    logger.info(
        f"""Divertor: Kallenbach 1D Model

    Major radius [m]: {rmajor}
    Minor radius [m]: {rminor}
    Toroidal field [T]: {bt}
    Plasma current [A]: {plascur}
    q95 [A]: {q}
    """
    )

    for i in range(1, impurity_radiation_module.nimp):
        impurity_radiation_module.impurity_arr_frac[i] = 0

    # Set the impurity array fraction of Nitrogen
    # gives 0.04 in SOL, as in Kallenbach paper
    impurity_radiation_module.impurity_arr_frac[4] = 8.0e-3

    divertor_ode.divertor_kallenbach(
        rmajor=rmajor,
        rminor=rminor,
        bt=bt,
        plascur=plascur,
        q=q,
        verboseset=False,
        ttarget=t_target,
        qtargettotal=q_target_total,
        targetangle=target_angle,
        unit_test=False,
        bp=b_pol,
        outfile=constants.nout,
        iprint=1,
    )

    read_and_get_atomic_data.plot_rates()
    logger.info(
        """Rate coefficients for deuterium - saved in "rate_coefficients.txt"
    Compare to Figure 2 in Kallenbach 2016"""
    )

    plot_radiation.plot_lz()
    logger.info(
        """Radiative loss functions - saved in "radiative_loss_functions.txt"
    Compare to Figure 3 in Kallenbach 2016."""
    )

    plot_radiation.plot_z()
    logger.info(
        """Reads mean Z and mean Z^2 - saved in "mean_Z.tx"
    Compare to plots such as He_z.ps etc in /home/mkovari/sol/kallenbach/divertor_ode/LZ_NON_CORONA."""
    )


if __name__ == "__main__":
    kallenbach_paper()
