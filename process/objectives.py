import numpy as np

from process.fortran import (
    physics_variables,
    tfcoil_variables,
    pf_power_variables,
    current_drive_variables,
    cost_variables,
    divertor_variables,
    times_variables,
    heat_transport_variables,
)


def objective_function(minmax: int) -> float:
    """Calculate the specified objective function

    :param minimax: the ID and sign of the figure of merit to evaluate.
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
    :type minimax: int
    """
    figure_of_merit = abs(minmax)

    # -1 = maximise
    # +1 = minimise
    objective_sign = np.sign(minmax)

    match figure_of_merit:
        case 1:
            objective_metric = 0.2 * physics_variables.rmajor
        case 3:
            objective_metric = physics_variables.wallmw
        case 4:
            objective_metric = (
                tfcoil_variables.tfcmw + 1e-3 * pf_power_variables.srcktpm
            ) / 10.0
        case 5:
            objective_metric = physics_variables.fusion_power / (
                current_drive_variables.pinjmw
                + current_drive_variables.porbitlossmw
                + physics_variables.pohmmw
            )
        case 6:
            objective_metric = cost_variables.coe / 100.0
        case 7:
            objective_metric = (
                cost_variables.cdirt / 1.0e3
                if cost_variables.ireactor == 0
                else cost_variables.concost / 1.0e4
            )
        case 8:
            objective_metric = physics_variables.aspect
        case 9:
            objective_metric = divertor_variables.hldiv
        case 10:
            objective_metric = physics_variables.bt
        case 11:
            objective_metric = current_drive_variables.pinjmw
        case 14:
            objective_metric = times_variables.t_burn / 2.0e4
        case 15:
            if cost_variables.iavail != 1:
                raise ValueError("minmax=15 requires iavail=1")
            objective_metric = cost_variables.cfactr
        case 16:
            objective_metric = 0.95 * (physics_variables.rmajor / 9.0) - 0.05 * (
                times_variables.t_burn / 7200.0
            )
        case 17:
            objective_metric = heat_transport_variables.pnetelmw / 500.0
        case 18:
            objective_metric = 1.0
        case 19:
            objective_metric = -0.5 * (current_drive_variables.bigq / 20.0) - 0.5 * (
                times_variables.t_burn / 7200.0
            )

    return objective_sign * objective_metric
