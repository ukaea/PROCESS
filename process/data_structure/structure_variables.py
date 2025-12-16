aintmass: float = None
"""intercoil structure mass (kg)"""

clgsmass: float = None
"""gravity support structure for TF coil, PF coil and intercoil support systems (kg)"""

coldmass: float = None
"""total mass of components at cryogenic temperatures (kg)"""

fncmass: float = None
"""PF coil outer support fence mass (kg)"""

gsmass: float = None
"""reactor core gravity support mass (kg)"""


def init_structure_variables():
    """Initialise structure variables"""
    global aintmass, clgsmass, coldmass, fncmass, gsmass

    aintmass = 0.0
    clgsmass = 0.0
    coldmass = 0.0
    fncmass = 0.0
    gsmass = 0.0
