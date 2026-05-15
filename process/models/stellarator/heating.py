import logging

from process.core import process_output as po
from process.core.exceptions import ProcessValueError
from process.core.model import DataStructure
from process.data_structure import physics_variables

logger = logging.getLogger(__name__)


def st_heat(stellarator, f_output: bool, data: DataStructure):
    """Routine to calculate the auxiliary heating power
    in a stellarator

    This routine calculates the auxiliary heating power for
    a stellarator device.

    Parameters
    ----------
    stellarator :

    f_output:

    data: DataStructure
        data structure object

    """
    f_p_beam_injected_ions = None
    if data.stellarator.isthtr == 1:
        data.current_drive.p_hcd_ecrh_injected_total_mw = (
            data.current_drive.p_hcd_primary_extra_heat_mw
        )
        data.current_drive.p_hcd_injected_ions_mw = 0
        data.current_drive.p_hcd_injected_electrons_mw = (
            data.current_drive.p_hcd_ecrh_injected_total_mw
        )
        data.current_drive.eta_hcd_primary_injector_wall_plug = (
            data.current_drive.eta_ecrh_injector_wall_plug
        )
        data.current_drive.p_hcd_electric_total_mw = (
            data.current_drive.p_hcd_injected_ions_mw
            + data.current_drive.p_hcd_injected_electrons_mw
        ) / data.current_drive.eta_hcd_primary_injector_wall_plug

    elif data.stellarator.isthtr == 2:
        data.current_drive.p_hcd_lowhyb_injected_total_mw = (
            data.current_drive.p_hcd_primary_extra_heat_mw
        )
        data.current_drive.p_hcd_injected_ions_mw = 0
        data.current_drive.p_hcd_injected_electrons_mw = (
            data.current_drive.p_hcd_lowhyb_injected_total_mw
        )
        data.current_drive.eta_hcd_primary_injector_wall_plug = (
            data.current_drive.eta_lowhyb_injector_wall_plug
        )
        data.heat_transport.p_hcd_electric_total_mw = (
            data.current_drive.p_hcd_injected_ions_mw
            + data.current_drive.p_hcd_injected_electrons_mw
        ) / data.current_drive.eta_hcd_primary_injector_wall_plug

    elif data.stellarator.isthtr == 3:
        (
            _effnbss,
            f_p_beam_injected_ions,
            data.current_drive.f_p_beam_shine_through,
        ) = stellarator.current_drive.culnbi()
        data.current_drive.p_hcd_beam_injected_total_mw = (
            data.current_drive.p_hcd_primary_extra_heat_mw
            * (1 - data.current_drive.f_p_beam_orbit_loss)
        )
        data.current_drive.p_beam_orbit_loss_mw = (
            data.current_drive.p_hcd_primary_extra_heat_mw
            * data.current_drive.f_p_beam_orbit_loss
        )
        data.current_drive.p_hcd_injected_ions_mw = (
            data.current_drive.p_hcd_beam_injected_total_mw * f_p_beam_injected_ions
        )
        data.current_drive.p_hcd_injected_electrons_mw = (
            data.current_drive.p_hcd_beam_injected_total_mw
            * (1 - f_p_beam_injected_ions)
        )
        data.current_drive.eta_hcd_primary_injector_wall_plug = (
            data.current_drive.eta_beam_injector_wall_plug
        )
        data.current_drive.p_hcd_electric_total_mw = (
            data.current_drive.p_hcd_injected_ions_mw
            + data.current_drive.p_hcd_injected_electrons_mw
        ) / data.current_drive.eta_hcd_primary_injector_wall_plug
    else:
        logger.error(f"isthtr {data.stellarator.isthtr}")
        logger.error(f"isthtr type {type(data.stellarator.isthtr)}")
        raise ProcessValueError(
            "Illegal value for isthtr", isthtr=data.stellarator.isthtr
        )

    #  Total injected power

    data.current_drive.p_hcd_injected_total_mw = (
        data.current_drive.p_hcd_injected_electrons_mw
        + data.current_drive.p_hcd_injected_ions_mw
    )

    #  Calculate neutral beam current

    if abs(data.current_drive.p_hcd_beam_injected_total_mw) > 1e-8:
        data.current_drive.c_beam_total = (
            1e-3
            * (data.current_drive.p_hcd_beam_injected_total_mw * 1e6)
            / data.current_drive.e_beam_kev
        )
    else:
        data.current_drive.c_beam_total = 0

    #  Ratio of fusion to input (injection+ohmic) power

    if (
        abs(
            data.current_drive.p_hcd_injected_total_mw
            + data.current_drive.p_beam_orbit_loss_mw
            + physics_variables.p_plasma_ohmic_mw
        )
        < 1e-6
    ):
        data.current_drive.big_q_plasma = 1e18
    else:
        data.current_drive.big_q_plasma = physics_variables.p_fusion_total_mw / (
            data.current_drive.p_hcd_injected_total_mw
            + data.current_drive.p_beam_orbit_loss_mw
            + physics_variables.p_plasma_ohmic_mw
        )

    if f_output:
        output(stellarator, data, f_p_beam_injected_ions)


def output(stellarator, data: DataStructure, f_p_beam_injected_ions=None):
    po.oheadr(stellarator.outfile, "Auxiliary Heating System")

    if data.stellarator.isthtr == 1:
        po.ocmmnt(stellarator.outfile, "Electron Cyclotron Resonance Heating")
    elif data.stellarator.isthtr == 2:
        po.ocmmnt(stellarator.outfile, "Lower Hybrid Heating")
    elif data.stellarator.isthtr == 3:
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
        data.current_drive.p_hcd_primary_extra_heat_mw,
    )
    po.ovarre(
        stellarator.outfile,
        "Fusion gain factor Q",
        "(big_q_plasma)",
        data.current_drive.big_q_plasma,
    )

    if abs(data.current_drive.p_hcd_beam_injected_total_mw) > 1e-8:
        po.ovarre(
            stellarator.outfile,
            "Neutral beam energy (KEV)",
            "(enbeam)",
            data.current_drive.enbeam,
        )
        po.ovarre(
            stellarator.outfile,
            "Neutral beam current (A)",
            "(c_beam_total)",
            data.current_drive.c_beam_total,
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
            data.current_drive.f_p_beam_shine_through,
        )
        po.ovarre(
            stellarator.outfile,
            "Neutral beam orbit loss power (MW)",
            "(p_beam_orbit_loss_mw)",
            data.current_drive.p_beam_orbit_loss_mw,
        )
        po.ovarre(
            stellarator.outfile,
            "Beam duct shielding thickness (m)",
            "(dx_beam_shield)",
            data.current_drive.dx_beam_shield,
        )
        po.ovarre(
            stellarator.outfile,
            "R injection tangent / R-major",
            "(f_radius_beam_tangency_rmajor)",
            data.current_drive.f_radius_beam_tangency_rmajor,
        )
        po.ovarre(
            stellarator.outfile,
            "Beam centreline tangency radius (m)",
            "(radius_beam_tangency)",
            data.current_drive.radius_beam_tangency,
        )
        po.ovarre(
            stellarator.outfile,
            "Maximum possible tangency radius (m)",
            "(radius_beam_tangency_max)",
            data.current_drive.radius_beam_tangency_max,
        )
        po.ovarre(
            stellarator.outfile,
            "Beam decay lengths to centre",
            "(n_beam_decay_lengths_core)",
            data.current_drive.n_beam_decay_lengths_core,
        )
