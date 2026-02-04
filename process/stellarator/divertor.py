import numpy as np

from process import constants
from process import process_output as po
from process.data_structure import (
    build_variables,
    divertor_variables,
    fwbs_variables,
    physics_variables,
    stellarator_variables,
)


def st_div(stellarator, output: bool):
    """Routine to call the stellarator divertor model
    author: P J Knight, CCFE, Culham Science Centre
    author: F Warmer, IPP Greifswald
    outfile : input integer : output file unit
    iprint : input integer : switch for writing to output file (1=yes)
    This routine calls the divertor model for a stellarator,
    developed by Felix Warmer.
    Stellarator Divertor Model for the Systems
    Code PROCESS, F. Warmer, 21/06/2013
    """
    Theta = stellarator_variables.flpitch  # ~bmn [rad] field line pitch
    r = physics_variables.rmajor
    p_div = physics_variables.p_plasma_separatrix_mw
    alpha = divertor_variables.anginc
    xi_p = divertor_variables.xpertin
    T_scrape = divertor_variables.tdiv

    #  Scrape-off temperature in Joules

    e = T_scrape * constants.ELECTRON_CHARGE

    #  Sound speed of particles (m/s)

    c_s = np.sqrt(e / (physics_variables.m_fuel_amu * constants.UMASS))

    #  Island size (m)

    w_r = 4.0e0 * np.sqrt(
        stellarator_variables.bmn
        * r
        / (stellarator_variables.shear * stellarator_variables.n_res)
    )

    #  Perpendicular (to plate) distance from X-point to divertor plate (m)

    Delta = stellarator_variables.f_w * w_r

    #  Length 'along' plasma (m)

    l_p = 2 * np.pi * r * (stellarator_variables.m_res) / stellarator_variables.n_res

    #  Connection length from X-point to divertor plate (m)

    l_x_t = Delta / Theta

    #  Power decay length (m)

    l_q = np.sqrt(xi_p * (l_x_t / c_s))

    #  Channel broadening length (m)

    l_b = np.sqrt(xi_p * l_p / (c_s))

    #  Channel broadening factor

    f_x = 1.0e0 + (l_b / (l_p * Theta))

    #  Length of a single divertor plate (m)

    l_d = f_x * l_p * (Theta / alpha)

    #  Total length of divertor plates (m)

    l_t = 2.0e0 * stellarator_variables.n_res * l_d

    #  Wetted area (m2)

    a_eff = l_t * l_q

    #  Divertor plate width (m): assume total area is wetted area/stellarator_variables.fdivwet

    darea = a_eff / stellarator_variables.fdivwet
    l_w = darea / l_t

    #  Divertor heat load (MW/m2)

    q_div = stellarator_variables.f_asym * (p_div / a_eff)

    #  Transfer to global variables

    divertor_variables.pflux_div_heat_load_mw = q_div
    divertor_variables.a_div_surface_total = darea

    fwbs_variables.f_ster_div_single = darea / build_variables.a_fw_total

    if output:
        output(stellarator, a_eff, l_d, l_w, f_x, l_q, w_r, Delta)


def output(stellarator, a_eff, l_d, l_w, f_x, l_q, w_r, Delta):
    """
    Outputs a summary of divertor-related parameters and results to the stellartor object.

    Parameters:
        stellarator: An object containing stellarator configuration and output handle.
        a_eff (float): Effective divertor wetted area (m²).
        l_d (float): Divertor plate length (m).
        l_w (float): Divertor plate width (m).
        f_x (float): Flux channel broadening factor.
        l_q (float): Power decay width (m).
        w_r (float): Island width (m).
        Delta (float): Perpendicular distance from X-point to plate (m).

    The function writes various physical and geometric parameters related to the divertor,
    including power, angles, heat transport coefficients, resonance numbers, field perturbations,
    and other relevant quantities, to the output file associated with the stellarator object.
    """

    po.oheadr(stellarator.outfile, "Divertor")

    po.ovarre(
        stellarator.outfile,
        "Power to divertor (MW)",
        "(p_plasma_separatrix_mw.)",
        physics_variables.p_plasma_separatrix_mw,
    )
    po.ovarre(
        stellarator.outfile,
        "Angle of incidence (deg)",
        "(anginc)",
        divertor_variables.anginc * 180.0e0 / np.pi,
    )
    po.ovarre(
        stellarator.outfile,
        "Perp. heat transport coefficient (m2/s)",
        "(xpertin)",
        divertor_variables.xpertin,
    )
    po.ovarre(
        stellarator.outfile,
        "Divertor plasma temperature (eV)",
        "(tdiv)",
        divertor_variables.tdiv,
    )
    po.ovarre(
        stellarator.outfile,
        "Radiated power fraction in SOL",
        "(f_rad)",
        stellarator_variables.f_rad,
    )
    po.ovarre(
        stellarator.outfile,
        "Heat load peaking factor",
        "(f_asym)",
        stellarator_variables.f_asym,
    )
    po.ovarin(
        stellarator.outfile,
        "Poloidal resonance number",
        "(m_res)",
        stellarator_variables.m_res,
    )
    po.ovarin(
        stellarator.outfile,
        "Toroidal resonance number",
        "(n_res)",
        stellarator_variables.n_res,
    )
    po.ovarre(
        stellarator.outfile,
        "Relative radial field perturbation",
        "(bmn)",
        stellarator_variables.bmn,
    )
    po.ovarre(
        stellarator.outfile,
        "Field line pitch (rad)",
        "(flpitch)",
        stellarator_variables.flpitch,
    )
    po.ovarre(
        stellarator.outfile,
        "Island size fraction factor",
        "(f_w)",
        stellarator_variables.f_w,
    )
    po.ovarre(
        stellarator.outfile,
        "Magnetic stellarator_variables.shear (/m)",
        "(shear)",
        stellarator_variables.shear,
    )
    po.ovarre(stellarator.outfile, "Divertor wetted area (m2)", "(A_eff)", a_eff)
    po.ovarre(
        stellarator.outfile,
        "Wetted area fraction of total plate area",
        "(fdivwet)",
        stellarator_variables.fdivwet,
    )
    po.ovarre(stellarator.outfile, "Divertor plate length (m)", "(L_d)", l_d)
    po.ovarre(stellarator.outfile, "Divertor plate width (m)", "(L_w)", l_w)
    po.ovarre(stellarator.outfile, "Flux channel broadening factor", "(F_x)", f_x)
    po.ovarre(stellarator.outfile, "Power decay width (cm)", "(100*l_q)", 100.0e0 * l_q)
    po.ovarre(stellarator.outfile, "Island width (m)", "(w_r)", w_r)
    po.ovarre(
        stellarator.outfile,
        "Perp. distance from X-point to plate (m)",
        "(Delta)",
        Delta,
    )
    po.ovarre(
        stellarator.outfile,
        "Peak heat load (MW/m2)",
        "(pflux_div_heat_load_mw)",
        divertor_variables.pflux_div_heat_load_mw,
    )
