"""author: J. Morris, M. Kovari (UKAEA)
Module containing global variables relating to the PF coil power conversion system
"""

import numpy as np

acptmax: float = None
"""average of currents in PF circuits (kA)"""

ensxpfm: float = None
"""maximum stored energy in the PF circuits (MJ)"""

i_pf_energy_storage_source: int = None
"""Switch for PF coil energy storage option:
 - =1 all power from MGF (motor-generator flywheel) units
 - =2 all pulsed power from line
 - =3 PF power from MGF, heating from line
   (In fact, options 1 and 3 are not treated differently)
"""

pfckts: float = None
"""number of PF coil circuits"""

spfbusl: float = None
"""total PF coil circuit bus length (m)"""

spsmva: float = None
"""sum of PF power supply ratings (MVA)"""

srcktpm: float = None
"""sum of resistive PF coil power (kW)"""

vpfskv: float = None
"""PF coil voltage (kV)"""

peakpoloidalpower: float = None
"""Peak absolute rate of change of stored energy in poloidal field (MW)"""

maxpoloidalpower: float = None
"""Maximum permitted absolute rate of change of stored energy in poloidal field (MW)"""

poloidalpower: list[float] = None
"""Poloidal power usage at time t (MW)"""

f_p_pf_energy_store_loss: float = None
"""Fraction of PF magnetic energy moved into/out of storage that is lost each time"""

f_p_pf_psu_loss: float = None
"""Fraction of inductive power flow lost in the PF power supplies/converters."""


def init_pf_power_variables():
    """Initialise PF coil power variables"""
    global \
        acptmax, \
        ensxpfm, \
        i_pf_energy_storage_source, \
        pfckts, \
        spfbusl, \
        spsmva, \
        srcktpm, \
        vpfskv, \
        peakpoloidalpower, \
        maxpoloidalpower, \
        poloidalpower, \
        f_p_pf_energy_store_loss, \
        f_p_pf_psu_loss

    acptmax = 0.0
    ensxpfm = 0.0
    i_pf_energy_storage_source = 2
    pfckts = 0.0
    spfbusl = 0.0
    spsmva = 0.0
    srcktpm = 0.0
    vpfskv = 0.0
    peakpoloidalpower = 0.0
    maxpoloidalpower = 1000.0
    poloidalpower = np.zeros(5)
    f_p_pf_energy_store_loss = 0.1
    f_p_pf_psu_loss = 0.1
