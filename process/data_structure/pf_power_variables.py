"""Module containing variables for the PF coil power models"""

from dataclasses import dataclass, field

import numpy as np

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@dataclass(slots=True)
class PFPowerData(PROCESSModelData):
    """Dataclass holding PF power variables"""

    acptmax: Parameter[float] = 0.0
    """average of currents in PF circuits (kA)"""

    ensxpfm: Parameter[float] = 0.0
    """maximum stored energy in the PF circuits (MJ)"""

    i_pf_energy_storage_source: int = 2
    """Switch for PF coil energy storage option:
    - =1 all power from MGF (motor-generator flywheel) units
    - =2 all pulsed power from line
    - =3 PF power from MGF, heating from line
    (In fact, options 1 and 3 are not treated differently)
    """

    pfckts: Parameter[float] = 0.0
    """number of PF coil circuits"""

    spfbusl: Parameter[float] = 0.0
    """total PF coil circuit bus length (m)"""

    spsmva: Parameter[float] = 0.0
    """sum of PF power supply ratings (MVA)"""

    srcktpm: Parameter[float] = 0.0
    """sum of resistive PF coil power (kW)"""

    vpfskv: Parameter[float] = 0.0
    """PF coil voltage (kV)"""

    peakpoloidalpower: Parameter[float] = 0.0
    """Peak absolute rate of change of stored energy in poloidal field (MW)"""

    maxpoloidalpower: Parameter[float] = 1000.0
    """Maximum permitted absolute rate of change of stored energy in poloidal field (MW)"""

    poloidalpower: Parameter[list[float]] = field(default_factory=lambda: np.zeros(5))
    """Poloidal power usage at time t (MW)"""

    f_p_pf_energy_store_loss: Parameter[float] = 0.1
    """Fraction of PF magnetic energy moved into/out of storage that is lost each time"""

    f_p_pf_psu_loss: Parameter[float] = 0.1
    """Fraction of inductive power flow lost in the PF power supplies/converters."""


CREATE_DICTS_FROM_DATACLASS = PFPowerData
