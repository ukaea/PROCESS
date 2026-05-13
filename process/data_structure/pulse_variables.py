from dataclasses import dataclass


@dataclass
class PulseData:
    bctmp: float = 320.0
    """first wall bulk coolant temperature (C)"""

    dtstor: float = 300.0
    """maximum allowable temperature change in stainless steel thermal storage block (K) (`istore=3`)"""

    istore: int = 1
    """Switch for thermal storage method:
    - =1 option 1 of Electrowatt report, AEA FUS 205
    - =2 option 2 of Electrowatt report, AEA FUS 205
    - =3 stainless steel block
    """

    itcycl: int = 1
    """Switch for first wall axial stress model:
    - =1 total axial constraint, no bending
    - =2 no axial constraint, no bending
    - =3 no axial constraint, bending
    """

    i_pulsed_plant: int = 0
    """Switch for reactor model:
    - =0 continuous operation
    - =1 pulsed operation
    """


CREATE_DICTS_FROM_DATACLASS = PulseData
