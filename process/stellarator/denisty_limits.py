import numpy as np
import logging
logger = logging.getLogger(__name__)

from process import process_output as po
from process.exceptions import ProcessValueError


from process.fortran import (
    physics_variables,
    stellarator_variables,
)

def st_denisty_limits(stellarator, output):
    """Routine to reiterate the physics loop
    author: J Lion, IPP Greifswald
    None
    This routine reiterates some physics modules.
    """

    physics_variables.dnelimt = st_sudo_density_limit(
        physics_variables.bt,
        physics_variables.p_plasma_loss_mw,
        physics_variables.rmajor,
        physics_variables.rminor,
    )

    # Calculates the ECRH parameters

    ne0_max_ECRH, bt_ecrh = st_d_limit_ecrh(
        stellarator_variables.max_gyrotron_frequency, physics_variables.bt
    )

    ne0_max_ECRH = min(physics_variables.ne0, ne0_max_ECRH)
    bt_ecrh = min(physics_variables.bt, bt_ecrh)

    if output:
        print_output(
            stellarator,
            bt_ecrh,
            ne0_max_ECRH,
        )


def st_sudo_density_limit(bt, powht, rmajor, rminor):
        """Routine to calculate the Sudo density limit in a stellarator
        author: P J Knight, CCFE, Culham Science Centre
        bt     : input real : Toroidal field on axis (T)
        powht  : input real : Absorbed heating power (MW)
        rmajor : input real : Plasma major radius (m)
        rminor : input real : Plasma minor radius (m)
        dlimit : output real : Maximum volume-averaged plasma density (/m3)
        This routine calculates the density limit for a stellarator.
        S.Sudo, Y.Takeiri, H.Zushi et al., Scalings of Energy Confinement
        and Density Limit in Stellarator/Heliotron Devices, Nuclear Fusion
        vol.30, 11 (1990).
        """
        arg = powht * bt / (rmajor * rminor * rminor)

        if arg <= 0.0e0:
            raise ProcessValueError(
                "Negative square root imminent",
                arg=arg,
                powht=powht,
                bt=bt,
                rmajor=rmajor,
                rminor=rminor,
            )

        #  Maximum line-averaged electron density

        dnlamx = 0.25e20 * np.sqrt(arg)

        #  Scale the result so that it applies to the volume-averaged
        #  electron density

        dlimit = dnlamx * physics_variables.dene / physics_variables.nd_electron_line

        #  Set the required value for icc=5

        physics_variables.dnelimt = dlimit

        return dlimit


def st_d_limit_ecrh(gyro_frequency_max, bt_input):
        """Routine to calculate the density limit due to an ECRH heating scheme on axis
        depending on an assumed maximal available gyrotron frequency.
        author: J Lion, IPP Greifswald
        gyro_frequency_max     : input real : Maximal available Gyrotron frequency (1/s) NOT (rad/s)
        bt  : input real : Maximal magnetic field on axis (T)
        dlimit_ecrh : output real : Maximum peak plasma density by ECRH constraints (/m3)
        bt_max : output real : Maximum allowable b field for ecrh heating (T)
        This routine calculates the density limit due to an ECRH heating scheme on axis
        """
        gyro_frequency = min(1.76e11 * bt_input, gyro_frequency_max * 2.0e0 * np.pi)

        # Restrict b field to the maximal available gyrotron frequency
        bt_max = (gyro_frequency_max * 2.0e0 * np.pi) / 1.76e11

        #                      me*e0/e^2       * w^2
        ne0_max = max(0.0e0, 3.142077e-4 * gyro_frequency**2)

        # Check if parabolic profiles are used:
        if physics_variables.ipedestal == 0:
            # Parabolic profiles used, use analytical formula:
            dlimit_ecrh = ne0_max
        else:
            logger.warning(
                "It was used physics_variables.ipedestal = 1 in a stellarator routine. PROCESS will pretend it got parabolic profiles (physics_variables.ipedestal = 0)."
            )
            dlimit_ecrh = ne0_max

        return dlimit_ecrh, bt_max


def print_output(stellarator, bt_ecrh, ne0_max_ECRH):
    po.oheadr(stellarator.outfile, "ECRH Ignition at lower values. Information:")

    po.ovarre(
        stellarator.outfile,
        "Maximal available gyrotron freq (input)",
        "(max_gyro_frequency)",
        stellarator_variables.max_gyrotron_frequency,
    )

    po.ovarre(
        stellarator.outfile,
        "Operating point: bfield",
        "(bt)",
        physics_variables.bt
    )

    po.ovarre(
        stellarator.outfile,
        "Operating point: Peak density",
        "(ne0)",
        physics_variables.ne0,
    )
    po.ovarre(
        stellarator.outfile,
        "Operating point: Peak temperature",
        "(te0)",
        physics_variables.te0,
    )

    po.ovarre(stellarator.outfile, "Ignition point: bfield (T)", "(bt_ecrh)", bt_ecrh)
    po.ovarre(
        stellarator.outfile,
        "Ignition point: density (/m3)",
        "(ne0_max_ECRH)",
        ne0_max_ECRH,
    )
    po.ovarre(
        stellarator.outfile,
        "Maximum reachable ECRH temperature (pseudo) (KEV)",
        "(te0_ecrh_achievable)",
        stellarator_variables.te0_ecrh_achievable,
    )

    powerht_local, pscalingmw_local = stellarator.power_at_ignition_point(
        stellarator_variables.max_gyrotron_frequency, stellarator_variables.te0_ecrh_achievable
    )
    po.ovarre(
        stellarator.outfile,
        "Ignition point: Heating Power (MW)",
        "(powerht_ecrh)",
        powerht_local,
    )
    po.ovarre(
        stellarator.outfile,
        "Ignition point: Loss Power (MW)",
        "(pscalingmw_ecrh)",
        pscalingmw_local,
    )

    if powerht_local >= pscalingmw_local:
        po.ovarin(stellarator.outfile, "Operation point ECRH ignitable?", "(ecrh_bool)", 1)
    else:
        po.ovarin(stellarator.outfile, "Operation point ECRH ignitable?", "(ecrh_bool)", 0)

