"""Module containing tokamak plasma physics routines
author: P J Knight, CCFE, Culham Science Centre
N/A
This module contains all the primary plasma physics routines
for a tokamak device.
"""

iscz: int = None


err242: int = None


err243: int = None


rad_fraction_lcfs: float = None


e_plasma_beta: float = None
"""[J]"""


total_loss_power: float = None
"""[W]"""


t_energy_confinement_beta: float = None
"""[s]"""


ptarmw: float = None


lambdaio: float = None


drsep: float = None


fio: float = None


fli: float = None


flo: float = None


fui: float = None


fuo: float = None


plimw: float = None


plomw: float = None


puimw: float = None


puomw: float = None


rho_star: float = None

nu_star: float = None

beta_mcdonald: float = None

itart_r: float = None

# Var in subroutine plasma_composition which requires re-initialisation on
# each new run:
first_call: int = None


def init_physics_module():
    """Initialise the physics module"""
    global first_call
    global iscz
    global err242
    global err243
    global rad_fraction_lcfs
    global e_plasma_beta
    global total_loss_power
    global t_energy_confinement_beta
    global ptarmw
    global lambdaio
    global drsep
    global fio
    global fli
    global flo
    global fui
    global fuo
    global plimw
    global plomw
    global puimw
    global puomw
    global rho_star
    global nu_star
    global beta_mcdonald
    global itart_r

    first_call = 1
    iscz = 0
    err242 = 0
    err243 = 0
    rad_fraction_lcfs = 0.0
    e_plasma_beta = 0.0
    total_loss_power = 0.0
    t_energy_confinement_beta = 0.0
    ptarmw = 0.0
    lambdaio = 0.0
    drsep = 0.0
    fio = 0.0
    fli = 0.0
    flo = 0.0
    fui = 0.0
    fuo = 0.0
    plimw = 0.0
    plomw = 0.0
    puimw = 0.0
    puomw = 0.0
    rho_star = 0.0
    nu_star = 0.0
    beta_mcdonald = 0.0
    itart_r = 0.0
