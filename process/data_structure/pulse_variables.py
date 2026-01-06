bctmp: float = None
"""first wall bulk coolant temperature (C)"""

dtstor: float = None
"""maximum allowable temperature change in stainless steel thermal storage block (K) (`istore=3`)"""

istore: int = None
"""Switch for thermal storage method:
- =1 option 1 of Electrowatt report, AEA FUS 205
- =2 option 2 of Electrowatt report, AEA FUS 205
- =3 stainless steel block
"""

itcycl: int = None
"""Switch for first wall axial stress model:
- =1 total axial constraint, no bending
- =2 no axial constraint, no bending
- =3 no axial constraint, bending
"""

i_pulsed_plant: int = None
"""Switch for reactor model:
- =0 continuous operation
- =1 pulsed operation
"""


def init_pulse_variables():
    """Initialise the pulse variables"""
    global bctmp, dtstor, istore, itcycl, i_pulsed_plant

    bctmp = 320.0
    dtstor = 300.0
    istore = 1
    itcycl = 1
    i_pulsed_plant = 0
