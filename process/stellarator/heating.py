import logging

from process import process_output as po
from process.data_structure import (
    current_drive_variables,
    heat_transport_variables,
    physics_variables,
    stellarator_variables,
)
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


def st_heat(stellarator, f_output: bool):
    """Routine to calculate the auxiliary heating power
    in a stellarator
    author: P J Knight, CCFE, Culham Science Centre
    outfile : input integer : output file unit
    iprint : input integer : switch for writing to output file (1=yes)
    This routine calculates the auxiliary heating power for
    a stellarator device.
    AEA FUS 172: Physics Assessment for the European Reactor Study
    """
    f_p_beam_injected_ions = None
    if stellarator_variables.isthtr == 1:
        current_drive_variables.p_hcd_ecrh_injected_total_mw = (
            current_drive_variables.p_hcd_primary_extra_heat_mw
        )
        current_drive_variables.p_hcd_injected_ions_mw = 0
        current_drive_variables.p_hcd_injected_electrons_mw = (
            current_drive_variables.p_hcd_ecrh_injected_total_mw
        )
        current_drive_variables.eta_hcd_primary_injector_wall_plug = (
            current_drive_variables.eta_ecrh_injector_wall_plug
        )
        current_drive_variables.p_hcd_electric_total_mw = (
            current_drive_variables.p_hcd_injected_ions_mw
            + current_drive_variables.p_hcd_injected_electrons_mw
        ) / current_drive_variables.eta_hcd_primary_injector_wall_plug

    elif stellarator_variables.isthtr == 2:
        current_drive_variables.p_hcd_lowhyb_injected_total_mw = (
            current_drive_variables.p_hcd_primary_extra_heat_mw
        )
        current_drive_variables.p_hcd_injected_ions_mw = 0
        current_drive_variables.p_hcd_injected_electrons_mw = (
            current_drive_variables.p_hcd_lowhyb_injected_total_mw
        )
        current_drive_variables.eta_hcd_primary_injector_wall_plug = (
            current_drive_variables.eta_lowhyb_injector_wall_plug
        )
        heat_transport_variables.p_hcd_electric_total_mw = (
            current_drive_variables.p_hcd_injected_ions_mw
            + current_drive_variables.p_hcd_injected_electrons_mw
        ) / current_drive_variables.eta_hcd_primary_injector_wall_plug

    elif stellarator_variables.isthtr == 3:
        (
            _effnbss,
            f_p_beam_injected_ions,
            current_drive_variables.f_p_beam_shine_through,
        ) = stellarator.current_drive.culnbi()
        current_drive_variables.p_hcd_beam_injected_total_mw = (
            current_drive_variables.p_hcd_primary_extra_heat_mw
            * (1 - current_drive_variables.f_p_beam_orbit_loss)
        )
        current_drive_variables.p_beam_orbit_loss_mw = (
            current_drive_variables.p_hcd_primary_extra_heat_mw
            * current_drive_variables.f_p_beam_orbit_loss
        )
        current_drive_variables.p_hcd_injected_ions_mw = (
            current_drive_variables.p_hcd_beam_injected_total_mw * f_p_beam_injected_ions
        )
        current_drive_variables.p_hcd_injected_electrons_mw = (
            current_drive_variables.p_hcd_beam_injected_total_mw
            * (1 - f_p_beam_injected_ions)
        )
        current_drive_variables.eta_hcd_primary_injector_wall_plug = (
            current_drive_variables.eta_beam_injector_wall_plug
        )
        current_drive_variables.p_hcd_electric_total_mw = (
            current_drive_variables.p_hcd_injected_ions_mw
            + current_drive_variables.p_hcd_injected_electrons_mw
        ) / current_drive_variables.eta_hcd_primary_injector_wall_plug
    else:
        logger.error(f"isthtr {stellarator_variables.isthtr}")
        logger.error(f"isthtr type {type(stellarator_variables.isthtr)}")
        raise ProcessValueError(
            "Illegal value for isthtr", isthtr=stellarator_variables.isthtr
        )

    #  Total injected power

    current_drive_variables.p_hcd_injected_total_mw = (
        current_drive_variables.p_hcd_injected_electrons_mw
        + current_drive_variables.p_hcd_injected_ions_mw
    )

    #  Calculate neutral beam current

    if abs(current_drive_variables.p_hcd_beam_injected_total_mw) > 1e-8:
        current_drive_variables.c_beam_total = (
            1e-3
            * (current_drive_variables.p_hcd_beam_injected_total_mw * 1e6)
            / current_drive_variables.e_beam_kev
        )
    else:
        current_drive_variables.c_beam_total = 0

    #  Ratio of fusion to input (injection+ohmic) power

    if (
        abs(
            current_drive_variables.p_hcd_injected_total_mw
            + current_drive_variables.p_beam_orbit_loss_mw
            + physics_variables.p_plasma_ohmic_mw
        )
        < 1e-6
    ):
        current_drive_variables.big_q_plasma = 1e18
    else:
        current_drive_variables.big_q_plasma = physics_variables.p_fusion_total_mw / (
            current_drive_variables.p_hcd_injected_total_mw
            + current_drive_variables.p_beam_orbit_loss_mw
            + physics_variables.p_plasma_ohmic_mw
        )

    if f_output:
        output(stellarator, f_p_beam_injected_ions)


def output(stellarator, f_p_beam_injected_ions=None):
    po.oheadr(stellarator.outfile, "Auxiliary Heating System")

    if stellarator_variables.isthtr == 1:
        po.ocmmnt(stellarator.outfile, "Electron Cyclotron Resonance Heating")
    elif stellarator_variables.isthtr == 2:
        po.ocmmnt(stellarator.outfile, "Lower Hybrid Heating")
    elif stellarator_variables.isthtr == 3:
        po.ocmmnt(stellarator.outfile, "Neutral Beam Injection Heating")

    if physics_variables.i_plasma_ignited == 1:
        po.ocmmnt(
            stellarator.outfile,
            "Ignited plasma; injected power only used for start-up phase",
        )

    po.oblnkl(stellarator.outfile)

    po.ovarre(
        stellarator.outfile,
        "Auxiliary power supplied to plasma (MW)",
        "(p_hcd_primary_extra_heat_mw)",
        current_drive_variables.p_hcd_primary_extra_heat_mw,
    )
    po.ovarre(
        stellarator.outfile,
        "Fusion gain factor Q",
        "(big_q_plasma)",
        current_drive_variables.big_q_plasma,
    )

    if abs(current_drive_variables.p_hcd_beam_injected_total_mw) > 1e-8:
        po.ovarre(
            stellarator.outfile,
            "Neutral beam energy (KEV)",
            "(enbeam)",
            current_drive_variables.enbeam,
        )
        po.ovarre(
            stellarator.outfile,
            "Neutral beam current (A)",
            "(c_beam_total)",
            current_drive_variables.c_beam_total,
        )
        po.ovarre(
            stellarator.outfile,
            "Fraction of beam energy to ions",
            "(f_p_beam_injected_ions)",
            f_p_beam_injected_ions,
        )
        po.ovarre(
            stellarator.outfile,
            "Neutral beam shine-through fraction",
            "(f_p_beam_shine_through)",
            current_drive_variables.f_p_beam_shine_through,
        )
        po.ovarre(
            stellarator.outfile,
            "Neutral beam orbit loss power (MW)",
            "(p_beam_orbit_loss_mw)",
            current_drive_variables.p_beam_orbit_loss_mw,
        )
        po.ovarre(
            stellarator.outfile,
            "Beam duct shielding thickness (m)",
            "(dx_beam_shield)",
            current_drive_variables.dx_beam_shield,
        )
        po.ovarre(
            stellarator.outfile,
            "R injection tangent / R-major",
            "(f_radius_beam_tangency_rmajor)",
            current_drive_variables.f_radius_beam_tangency_rmajor,
        )
        po.ovarre(
            stellarator.outfile,
            "Beam centreline tangency radius (m)",
            "(radius_beam_tangency)",
            current_drive_variables.radius_beam_tangency,
        )
        po.ovarre(
            stellarator.outfile,
            "Maximum possible tangency radius (m)",
            "(radius_beam_tangency_max)",
            current_drive_variables.radius_beam_tangency_max,
        )
        po.ovarre(
            stellarator.outfile,
            "Beam decay lengths to centre",
            "(n_beam_decay_lengths_core)",
            current_drive_variables.n_beam_decay_lengths_core,
        )
