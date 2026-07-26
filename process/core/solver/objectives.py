import numpy as np

from process.core.exceptions import ProcessValueError
from process.core.model import DataStructure
from process.data_structure.numerics import FiguresOfMerit


def objective_function(minmax: int, data: DataStructure) -> float:
    """Calculate the specified objective function

    Parameters
    ----------
    minmax : int
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
    try:
        figure_of_merit = FiguresOfMerit(abs(minmax))
    except ValueError as err:
        raise ProcessValueError(f"Invalid minmax value: {minmax}") from err

    # -1 = maximise
    # +1 = minimise
    objective_sign = np.sign(minmax)

    if figure_of_merit == FiguresOfMerit.MAJOR_RADIUS:
        objective_metric = 0.2 * data.physics.rmajor
    elif figure_of_merit == FiguresOfMerit.NEUTRON_WALL_LOAD:
        objective_metric = data.physics.pflux_fw_neutron_mw
    elif figure_of_merit == FiguresOfMerit.P_TF_PLUS_P_PF:
        objective_metric = (data.tfcoil.tfcmw + 1e-3 * data.pf_power.srcktpm) / 10.0
    elif figure_of_merit == FiguresOfMerit.FUSION_GAIN_Q:
        objective_metric = data.current_drive.big_q_plasma
    elif figure_of_merit == FiguresOfMerit.COST_OF_ELECTRICITY:
        objective_metric = data.costs.coe / 100.0
    elif figure_of_merit == FiguresOfMerit.CAPITAL_COST:
        objective_metric = (
            data.costs.cdirt / 1.0e3
            if data.costs.ireactor == 0
            else data.costs.concost / 1.0e4
        )
    elif figure_of_merit == FiguresOfMerit.ASPECT_RATIO:
        objective_metric = data.physics.aspect
    elif figure_of_merit == FiguresOfMerit.DIVERTOR_HEAT_LOAD:
        objective_metric = data.divertor.pflux_div_heat_load_mw
    elif figure_of_merit == FiguresOfMerit.TOROIDAL_FIELD:
        objective_metric = data.physics.b_plasma_toroidal_on_axis
    elif figure_of_merit == FiguresOfMerit.TOTAL_INJECTED_POWER:
        objective_metric = data.current_drive.p_hcd_injected_total_mw
    elif figure_of_merit == FiguresOfMerit.PULSE_LENGTH:
        objective_metric = data.times.t_plant_pulse_burn / 2.0e4
    elif figure_of_merit == FiguresOfMerit.PLANT_AVAILABILITY_FACTOR:
        if data.costs.i_plant_availability != 1:
            raise ProcessValueError("minmax=15 requires i_plant_availability=1")
        objective_metric = data.costs.f_t_plant_available
    elif figure_of_merit == FiguresOfMerit.MIN_R0_MAX_TAU_BURN:
        objective_metric = 0.95 * (data.physics.rmajor / 9.0) - 0.05 * (
            data.times.t_plant_pulse_burn / 7200.0
        )
    elif figure_of_merit == FiguresOfMerit.NET_ELECTRICAL_OUTPUT:
        objective_metric = data.heat_transport.p_plant_electric_net_mw / 500.0
    elif figure_of_merit == FiguresOfMerit.NULL_FIGURE_OF_MERIT:
        objective_metric = 1.0
    elif figure_of_merit == FiguresOfMerit.MAX_Q_MAX_T_PLANT_PULSE_BURN:
        objective_metric = -0.5 * (data.current_drive.big_q_plasma / 20.0) - 0.5 * (
            data.times.t_plant_pulse_burn / 7200.0
        )

    return objective_sign * objective_metric
