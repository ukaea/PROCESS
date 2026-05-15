import numpy as np

from process.core.exceptions import ProcessValueError
from process.core.model import DataStructure
from process.data_structure import (
    divertor_variables,
    pf_power_variables,
    physics_variables,
    tfcoil_variables,
)

OBJECTIVE_NAMES = {
    1: "Plasma major radius",
    3: "neutron wall load",
    4: "total TF + PF coil power",
    5: "ratio fusion power:injection power",
    6: "cost of electricity",
    7: "constructed cost",
    8: "aspect ratio",
    9: "divertor heat load",
    10: "toroidal field on axis",
    11: "injection power",
    14: "pulse length",
    15: "plant availability factor",
    16: "Major radius/burn time",
    17: "net electrical output",
    18: "NULL",
    19: "Major radius/burn time",
}


def objective_function(minmax: int, data: DataStructure) -> float:
    """Calculate the specified objective function

    Parameters
    ----------
    minimax : int
        the ID and sign of the figure of merit to evaluate.
        A negative value indicates maximisation.
        A positive value indicates minimisation.
        * 1: Major radius
        * 3: Neutron wall load
        * 4: TF coil + PF coil power
        * 5: Fusion gain
        * 6: Cost of electricity
        * 7: Direct/constructed/capital cost
        * 8: Aspect ratio
        * 9: Divertor heat load
        * 10: Toroidal field on axis
        * 11: Injected power
        * 14: Pulse length
        * 15: Plant availability
        * 16: Major radius/burn time
        * 17: Net electrical output
        * 18: NULL, f(x) = 1
        * 19: Major radius/burn time
    data: DataStructure
        data structure object for providing data to the
        objective function
    """
    figure_of_merit = abs(minmax)

    # -1 = maximise
    # +1 = minimise
    objective_sign = np.sign(minmax)

    match figure_of_merit:
        case 1:
            objective_metric = 0.2 * physics_variables.rmajor
        case 3:
            objective_metric = physics_variables.pflux_fw_neutron_mw
        case 4:
            objective_metric = (
                tfcoil_variables.tfcmw + 1e-3 * pf_power_variables.srcktpm
            ) / 10.0
        case 5:
            objective_metric = physics_variables.p_fusion_total_mw / (
                data.current_drive.p_hcd_injected_total_mw
                + data.current_drive.p_beam_orbit_loss_mw
                + physics_variables.p_plasma_ohmic_mw
            )
        case 6:
            objective_metric = data.costs.coe / 100.0
        case 7:
            objective_metric = (
                data.costs.cdirt / 1.0e3
                if data.costs.ireactor == 0
                else data.costs.concost / 1.0e4
            )
        case 8:
            objective_metric = physics_variables.aspect
        case 9:
            objective_metric = divertor_variables.pflux_div_heat_load_mw
        case 10:
            objective_metric = physics_variables.b_plasma_toroidal_on_axis
        case 11:
            objective_metric = data.current_drive.p_hcd_injected_total_mw
        case 14:
            objective_metric = data.times.t_plant_pulse_burn / 2.0e4
        case 15:
            if data.costs.i_plant_availability != 1:
                raise ProcessValueError("minmax=15 requires i_plant_availability=1")
            objective_metric = data.costs.f_t_plant_available
        case 16:
            objective_metric = 0.95 * (physics_variables.rmajor / 9.0) - 0.05 * (
                data.times.t_plant_pulse_burn / 7200.0
            )
        case 17:
            objective_metric = data.heat_transport.p_plant_electric_net_mw / 500.0
        case 18:
            objective_metric = 1.0
        case 19:
            objective_metric = -0.5 * (data.current_drive.big_q_plasma / 20.0) - 0.5 * (
                data.times.t_plant_pulse_burn / 7200.0
            )

    return objective_sign * objective_metric
