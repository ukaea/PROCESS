from copy import copy
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

    #  Set the required value for icc=5
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
        if physics_variables.ipedestal != 0:
            logger.warning(
                "It was used physics_variables.ipedestal != 0 in a stellarator routine. PROCESS will pretend it got parabolic profiles (physics_variables.ipedestal = 0)."
            )
        # Assume parabolic profiles anyway, use analytical formula:  
        dlimit_ecrh = ne0_max

        return dlimit_ecrh, bt_max


def power_at_ignition_point(stellarator, gyro_frequency_max, te0_available):
        """Routine to calculate if the plasma is ignitable with the current values for the B field. Assumes
        current ECRH achievable peak temperature (which is inaccurate as the cordey pass should be calculated)
        author: J Lion, IPP Greifswald
        gyro_frequency_max : input real : Maximal available Gyrotron frequency (1/s) NOT (rad/s)
        te0_available : input real : Reachable peak electron temperature, reached by ECRH (KEV)
        powerht_out : output real: Heating Power at ignition point (MW)
        pscalingmw_out : output real: Heating Power loss at ignition point (MW)
        This routine calculates the density limit due to an ECRH heating scheme on axis
        Assumes current peak temperature (which is inaccurate as the cordey pass should be calculated)
        Maybe use this: https://doi.org/10.1088/0029-5515/49/8/085026
        """
        te_old = copy(physics_variables.te)
        # Volume averaged physics_variables.te from te0_achievable
        physics_variables.te = te0_available / (1.0e0 + physics_variables.alphat)
        ne0_max, bt_ecrh_max = st_d_limit_ecrh(
            gyro_frequency_max, physics_variables.bt
        )
        # Now go to point where ECRH is still available
        # In density..
        dene_old = copy(physics_variables.dene)
        physics_variables.dene = min(
            dene_old, ne0_max / (1.0e0 + physics_variables.alphan)
        )

        # And B-field..
        bt_old = copy(physics_variables.bt)
        physics_variables.bt = min(bt_ecrh_max, physics_variables.bt)

        stellarator.st_phys(False)
        stellarator.st_phys(
            False
        )  # The second call seems to be necessary for all values to "converge" (and is sufficient)

        powerht_out = max(
            copy(physics_variables.p_plasma_loss_mw), 0.00001e0
        )  # the radiation module sometimes returns negative heating power
        pscalingmw_out = copy(physics_variables.pscalingmw)

        # Reverse it and do it again because anything more efficiently isn't suitable with the current implementation
        # This is bad practice but seems to be necessary as of now:
        physics_variables.te = te_old
        physics_variables.dene = dene_old
        physics_variables.bt = bt_old

        # The second call seems to be necessary for all values to "converge" (and is sufficient)
        stellarator.st_phys(False)
        stellarator.st_phys(False)

        return powerht_out, pscalingmw_out


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

    powerht_local, pscalingmw_local = power_at_ignition_point(stellarator, 
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

