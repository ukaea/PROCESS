from dataclasses import astuple, dataclass

import numpy as np
from tabulate import tabulate

import process.process_output as process_output
from process.caller import write_output_files
from process.exceptions import ProcessValueError
from process.fortran import (
    build_variables,
    constants,
    constraint_variables,
    constraints,
    cost_variables,
    cs_fatigue_variables,
    current_drive_variables,
    divertor_variables,
    error_handling,
    fwbs_variables,
    global_variables,
    heat_transport_variables,
    impurity_radiation_module,
    numerics,
    pfcoil_variables,
    physics_variables,
    rebco_variables,
    scan_module,
    tfcoil_variables,
)
from process.optimiser import Optimiser
from process.utilities.f2py_string_patch import (
    f2py_compatible_to_string,
    string_to_f2py_compatible,
)


@dataclass
class ScanVariable:
    variable_name: str
    variable_description: str

    def __iter__(self):
        return iter(astuple(self))


SCAN_VARIABLES = {
    1: ScanVariable("aspect", "Aspect_ratio"),
    2: ScanVariable("pflux_div_heat_load_max_mw", "Div_heat_limit_(MW/m2)"),
    3: ScanVariable("pnetelin", "Net_electric_power_(MW)"),
    4: ScanVariable("hfact", "Confinement_H_factor"),
    5: ScanVariable("oacdcp", "TF_inboard_leg_J_(MA/m2)"),
    6: ScanVariable("walalw", "Allow._wall_load_(MW/m2)"),
    7: ScanVariable("beamfus0", "Beam_bkgrd_multiplier"),
    8: ScanVariable("fqval", "Big_Q_f-value"),
    9: ScanVariable("te", "Electron_temperature_keV"),
    10: ScanVariable("boundu(15)", "Volt-second_upper_bound"),
    11: ScanVariable("beta_norm_max", "Beta_coefficient"),
    12: ScanVariable("f_c_plasma_bootstrap_max", "Bootstrap_fraction"),
    13: ScanVariable("boundu(10)", "H_factor_upper_bound"),
    14: ScanVariable("fiooic", "TFC_Iop_/_Icrit_f-value"),
    15: ScanVariable("fjprot", "TFC_Jprot_limit_f-value"),
    16: ScanVariable("rmajor", "Plasma_major_radius_(m)"),
    17: ScanVariable("bmxlim", "Max_toroidal_field_(T)"),
    18: ScanVariable("gammax", "Maximum_CD_gamma"),
    19: ScanVariable("boundl(16)", "CS_thickness_lower_bound"),
    20: ScanVariable("t_burn_min", "Minimum_burn_time_(s)"),
    22: ScanVariable("cfactr", "Plant_availability_factor"),
    23: ScanVariable("boundu(72)", "Ip/Irod_upper_bound"),
    24: ScanVariable("p_fusion_total_max_mw", "Fusion_power_limit_(MW)"),
    25: ScanVariable("kappa", "Plasma_elongation"),
    26: ScanVariable("triang", "Plasma_triangularity"),
    27: ScanVariable("tbrmin", "Min_tritium_breed._ratio"),
    28: ScanVariable("bt", "Tor._field_on_axis_(T)"),
    29: ScanVariable("coreradius", "Core_radius"),
    31: ScanVariable(
        "f_alpha_energy_confinement_min", "t_alpha_confinement/taueff_lower_limit"
    ),
    32: ScanVariable("epsvmc", "VMCON error tolerance"),
    38: ScanVariable("boundu(129)", " Neon upper limit"),
    39: ScanVariable("boundu(131)", " Argon upper limit"),
    40: ScanVariable("boundu(135)", " Xenon upper limit"),
    41: ScanVariable("dr_blkt_outboard", "Outboard blanket thick."),
    42: ScanVariable("fimp(9)", "Argon fraction"),
    44: ScanVariable("sig_tf_case_max", "Allowable_stress_in_tf_coil_case_Tresca_(pa)"),
    45: ScanVariable("tmargmin_tf", "Minimum_allowable_temperature_margin"),
    46: ScanVariable("boundu(152)", "Max allowable fgwsep"),
    48: ScanVariable("n_pancake", "TF Coil - n_pancake"),
    49: ScanVariable("n_layer", "TF Coil - n_layer"),
    50: ScanVariable("fimp(13)", "Xenon fraction"),
    51: ScanVariable("f_p_div_lower", "lower_divertor_power_fraction"),
    52: ScanVariable("rad_fraction_sol", "SoL radiation fraction"),
    53: ScanVariable("boundu(157)", "Max allowable fvssu"),
    54: ScanVariable("Bc2(0K)", "GL_NbTi Bc2(0K)"),
    55: ScanVariable("dr_shld_inboard", "Inboard neutronic shield"),
    56: ScanVariable("crypmw_max", "max allowable crypmw"),
    57: ScanVariable("boundl(2)", "bt minimum"),
    58: ScanVariable("dr_fw_plasma_gap_inboard", "Inboard FW-plasma sep gap"),
    59: ScanVariable("dr_fw_plasma_gap_outboard", "Outboard FW-plasma sep gap"),
    60: ScanVariable(
        "sig_tf_wp_max", "Allowable_stress_in_tf_coil_conduit_Tresca_(pa)"
    ),
    61: ScanVariable("copperaoh_m2_max", "Max CS coil current / copper area"),
    62: ScanVariable("coheof", "CS coil current density at EOF (A/m2)"),
    63: ScanVariable("dr_cs", "CS coil thickness (m)"),
    64: ScanVariable("ohhghf", "CS height (m)"),
    65: ScanVariable("n_cycle_min", "CS stress cycles min"),
    66: ScanVariable("oh_steel_frac", "CS steel fraction"),
    67: ScanVariable("t_crack_vertical", "Initial crack vertical size (m)"),
    68: ScanVariable(
        "inlet_temp_liq", "Inlet Temperature Liquid Metal Breeder/Coolant (K)"
    ),
    69: ScanVariable(
        "outlet_temp_liq", "Outlet Temperature Liquid Metal Breeder/Coolant (K)"
    ),
    70: ScanVariable(
        "blpressure_liq", "Blanket liquid metal breeder/coolant pressure (Pa)"
    ),
    71: ScanVariable(
        "n_liq_recirc", "Selected number of liquid metal breeder recirculations per day"
    ),
    72: ScanVariable(
        "bz_channel_conduct_liq",
        "Conductance of liquid metal breeder duct walls (A V-1 m-1)",
    ),
    73: ScanVariable(
        "pnuc_fw_ratio_dcll", "Ratio of FW nuclear power as fraction of total (FW+BB)"
    ),
    74: ScanVariable(
        "f_nuc_pow_bz_struct",
        "Fraction of BZ power cooled by primary coolant for dual-coolant balnket",
    ),
    75: ScanVariable("dx_fw_module", "dx_fw_module of first wall cooling channels (m)"),
    76: ScanVariable("etath", "Thermal conversion eff."),
    77: ScanVariable("startupratio", "Gyrotron redundancy"),
    78: ScanVariable("fkind", "Multiplier for Nth of a kind costs"),
    79: ScanVariable(
        "eta_ecrh_injector_wall_plug", "ECH wall plug to injector efficiency"
    ),
    80: ScanVariable("fcoolcp", "Coolant fraction of TF"),
    81: ScanVariable("n_tf_turn", "Number of turns in TF"),
}


class Scan:
    """Perform a parameter scan using the Fortran scan module."""

    def __init__(self, models, solver):
        """Immediately run the run_scan() method.

        :param models: physics and engineering model objects
        :type models: process.main.Models
        :param solver: which solver to use, as specified in solver.py
        :type solver: str
        """
        self.models = models
        self.optimiser = Optimiser(models, solver)
        self.run_scan()

    def run_scan(self):
        """Call VMCON over a range of values of one of the variables.

        This method calls the optimisation routine VMCON a number of times, by
        performing a sweep over a range of values of a particular variable. A
        number of output variable values are written to the MFILE.DAT file at
        each scan point, for plotting or other post-processing purposes.
        """
        # Turn off error reporting (until next output)
        error_handling.errors_on = False

        if scan_module.isweep == 0:
            ifail = self.doopt()
            write_output_files(models=self.models, ifail=ifail)
            error_handling.show_errors()
            return

        if scan_module.isweep > scan_module.ipnscns:
            raise ProcessValueError(
                "Illegal value of isweep",
                isweep=scan_module.isweep,
                ipnscns=scan_module.ipnscns,
            )

        if scan_module.scan_dim == 2:
            self.scan_2d()
        else:
            self.scan_1d()

    def doopt(self):
        """Run the optimiser."""
        # If no optimisation is required, leave the method
        if numerics.ioptimz < 0:
            return None

        ifail = self.optimiser.run()
        self.post_optimise(ifail)

        return ifail

    def post_optimise(self, ifail: int):
        """Called after calling the optimising equation solver from Python.
        author: P J Knight, CCFE, Culham Science Centre
        ifail   : input integer : error flag
        """
        numerics.sqsumsq = (numerics.rcm[: numerics.neqns] ** 2).sum() ** 0.5

        error_handling.errors_on = True

        process_output.oheadr(constants.nout, "Numerics")
        process_output.ocmmnt(
            constants.nout, "PROCESS has performed a VMCON (optimisation) run."
        )
        if ifail != 1:
            process_output.ovarin(constants.nout, "VMCON error flag", "(ifail)", ifail)
            process_output.oheadr(
                constants.iotty, "PROCESS COULD NOT FIND A FEASIBLE SOLUTION"
            )
            process_output.oblnkl(constants.iotty)

            error_handling.idiags[0] = ifail
            error_handling.report_error(132)

            self.verror(ifail)
            process_output.oblnkl(constants.nout)
            process_output.oblnkl(constants.iotty)
        else:
            process_output.ocmmnt(
                constants.nout, "and found a feasible set of parameters."
            )
            process_output.oblnkl(constants.nout)
            process_output.ovarin(constants.nout, "VMCON error flag", "(ifail)", ifail)
            process_output.oheadr(constants.iotty, "PROCESS found a feasible solution")

            if numerics.sqsumsq >= 1.0e-2:
                process_output.oblnkl(constants.nout)
                process_output.ocmmnt(
                    constants.nout,
                    "WARNING: Constraint residues are HIGH; consider re-running",
                )
                process_output.ocmmnt(
                    constants.nout,
                    "   with lower values of EPSVMC to confirm convergence...",
                )
                process_output.ocmmnt(
                    constants.nout,
                    "   (should be able to get down to about 1.0E-8 okay)",
                )
                process_output.oblnkl(constants.nout)
                process_output.ocmmnt(
                    constants.iotty,
                    "WARNING: Constraint residues are HIGH; consider re-running",
                )
                process_output.ocmmnt(
                    constants.iotty,
                    "   with lower values of EPSVMC to confirm convergence...",
                )
                process_output.ocmmnt(
                    constants.iotty,
                    "   (should be able to get down to about 1.0E-8 okay)",
                )
                process_output.oblnkl(constants.iotty)

                error_handling.fdiags[0] = numerics.sqsumsq
                error_handling.report_error(134)

        process_output.ovarin(
            constants.nout, "Number of iteration variables", "(nvar)", numerics.nvar
        )
        process_output.ovarin(
            constants.nout,
            "Number of constraints (total)",
            "(neqns+nineqns)",
            numerics.neqns + numerics.nineqns,
        )
        process_output.ovarin(
            constants.nout, "Optimisation switch", "(ioptimz)", numerics.ioptimz
        )
        process_output.ovarin(
            constants.nout, "Figure of merit switch", "(minmax)", numerics.minmax
        )

        objf_name = string_to_f2py_compatible(
            numerics.objf_name,
            f'"{f2py_compatible_to_string(numerics.lablmm[abs(numerics.minmax) - 1])}"',
        )

        numerics.objf_name = objf_name

        process_output.ovarst(
            constants.nout, "Objective function name", "(objf_name)", numerics.objf_name
        )
        process_output.ovarre(
            constants.nout,
            "Normalised objective function",
            "(norm_objf)",
            numerics.norm_objf,
            "OP ",
        )
        process_output.ovarre(
            constants.nout,
            "Square root of the sum of squares of the constraint residuals",
            "(sqsumsq)",
            numerics.sqsumsq,
            "OP ",
        )
        process_output.ovarre(
            constants.nout,
            "VMCON convergence parameter",
            "(convergence_parameter)",
            global_variables.convergence_parameter,
            "OP ",
        )
        process_output.ovarin(
            constants.nout,
            "Number of VMCON iterations",
            "(nviter)",
            numerics.nviter,
            "OP ",
        )
        process_output.oblnkl(constants.nout)

        if ifail == 1:
            string1 = "PROCESS has successfully optimised"
        else:
            string1 = "PROCESS has tried to optimise"

        string2 = "minimise" if numerics.minmax > 0 else "maximise"

        process_output.write(
            constants.nout,
            f"{string1} the iteration variables to {string2} the figure of merit: {objf_name}\n",
        )

        written_warning = False

        solution_vector_table = []
        for i in range(numerics.nvar):
            numerics.xcs[i] = numerics.xcm[i] * numerics.scafc[i]

            name = f2py_compatible_to_string(numerics.lablxc[numerics.ixc[i] - 1])
            solution_vector_table.append([name, numerics.xcs[i], numerics.xcm[i]])

            xminn = 1.01 * numerics.itv_scaled_lower_bounds[i]
            xmaxx = 0.99 * numerics.itv_scaled_upper_bounds[i]

            if numerics.xcm[i] < xminn or numerics.xcm[i] > xmaxx:
                if not written_warning:
                    written_warning = True
                    process_output.ocmmnt(
                        constants.nout,
                        (
                            "Certain operating limits have been reached,"
                            "\n as shown by the following iteration variables that are"
                            "\n at or near to the edge of their prescribed range :\n"
                        ),
                    )

                xcval = numerics.xcm[i] * numerics.scafc[i]

                if numerics.xcm[i] < xminn:
                    location, bound = "below", "lower"
                else:
                    location, bound = "above", "upper"
                process_output.write(
                    constants.nout,
                    f"   {name:<30}= {xcval} is at or {location} its {bound} bound:"
                    f" {numerics.itv_scaled_upper_bounds[i] * numerics.scafc[i]}",
                )

            process_output.ovarre(
                constants.mfile,
                f2py_compatible_to_string(numerics.lablxc[numerics.ixc[i] - 1]),
                f"(itvar{i + 1:03d})",
                numerics.xcs[i],
            )

            if numerics.boundu[i] == numerics.boundl[i]:
                xnorm = 1.0
            else:
                xnorm = min(
                    max(
                        (numerics.xcm[i] - numerics.itv_scaled_lower_bounds[i])
                        / (
                            numerics.itv_scaled_upper_bounds[i]
                            - numerics.itv_scaled_lower_bounds[i]
                        ),
                        0.0,
                    ),
                    1.0,
                )

            process_output.ovarre(
                constants.mfile,
                f"{name} (final value/initial value)",
                f"(xcm{i + 1:03d})",
                numerics.xcm[i],
            )
            process_output.ovarre(
                constants.mfile,
                f"{name} (range normalised)",
                f"(nitvar{i + 1:03d})",
                xnorm,
            )

        process_output.osubhd(
            constants.nout, "The solution vector is comprised as follows :"
        )
        process_output.write(
            constants.nout,
            tabulate(
                solution_vector_table,
                headers=["", "Final value", "Final / initial"],
                numalign="left",
            ),
        )

        process_output.osubhd(
            constants.nout,
            "The following equality constraint residues should be close to zero :",
        )

        con1, con2, err, sym, lab = constraints.constraint_eqns(
            numerics.neqns + numerics.nineqns, -1
        )

        equality_constraint_table = []
        for i in range(numerics.neqns):
            name = f2py_compatible_to_string(numerics.lablcc[numerics.icc[i] - 1])

            equality_constraint_table.append([
                name,
                sym[i],
                f"{con2[i]} {f2py_compatible_to_string(lab[i])}",
                f"{err[i]} {f2py_compatible_to_string(lab[i])}",
                con1[i],
            ])
            process_output.ovarre(
                constants.mfile,
                f"{name:<33} normalised residue",
                f"(eq_con{numerics.icc[i]:03d})",
                con1[i],
            )

        process_output.write(
            constants.nout,
            tabulate(
                equality_constraint_table,
                headers=[
                    "",
                    "",
                    "Physical constraint",
                    "Constraint residue",
                    "Normalised residue",
                ],
                numalign="left",
            ),
        )

        if numerics.nineqns > 0:
            inequality_constraint_table = []
            process_output.osubhd(
                constants.nout,
                "The following inequality constraint residues should be "
                "greater than or approximately equal to zero :",
            )

            for i in range(numerics.neqns, numerics.neqns + numerics.nineqns):
                name = f2py_compatible_to_string(numerics.lablcc[numerics.icc[i] - 1])
                inequality_constraint_table.append([
                    name,
                    sym[i],
                    f"{con2[i]} {f2py_compatible_to_string(lab[i])}",
                    f"{err[i]} {f2py_compatible_to_string(lab[i])}",
                ])
                process_output.ovarre(
                    constants.mfile,
                    f"{name} normalised residue",
                    f"(ineq_con{numerics.icc[i]:03d})",
                    numerics.rcm[i],
                )

            process_output.write(
                constants.nout,
                tabulate(
                    inequality_constraint_table,
                    headers=[
                        "",
                        "",
                        "Physical constraint",
                        "Constraint residue",
                    ],
                    numalign="left",
                ),
            )

    def verror(self, ifail: int):
        """Routine to print out relevant messages in the case of an
        unfeasible result from a VMCON (optimisation) run
        author: P J Knight, CCFE, Culham Science Centre
        ifail  : input integer : error flag
        This routine prints out relevant messages in the case of
        an unfeasible result from a VMCON (optimisation) run.
        """
        if ifail == -1:
            process_output.ocmmnt(constants.nout, "User-terminated execution of VMCON.")
            process_output.ocmmnt(
                constants.iotty, "User-terminated execution of VMCON."
            )
        elif ifail == 0:
            process_output.ocmmnt(
                constants.nout, "Improper input parameters to the VMCON routine."
            )
            process_output.ocmmnt(constants.nout, "PROCESS coding must be checked.")

            process_output.ocmmnt(
                constants.iotty, "Improper input parameters to the VMCON routine."
            )
            process_output.ocmmnt(constants.iotty, "PROCESS coding must be checked.")
        elif ifail == 2:
            process_output.ocmmnt(
                constants.nout,
                "The maximum number of calls has been reached without solution.",
            )
            process_output.ocmmnt(
                constants.nout,
                "The code may be stuck in a minimum in the residual space that is significantly above zero.",
            )
            process_output.oblnkl(constants.nout)
            process_output.ocmmnt(
                constants.nout, "There is either no solution possible, or the code"
            )
            process_output.ocmmnt(
                constants.nout, "is failing to escape from a deep local minimum."
            )
            process_output.ocmmnt(
                constants.nout,
                "Try changing the variables in IXC, or modify their initial values.",
            )

            process_output.ocmmnt(
                constants.iotty,
                "The maximum number of calls has been reached without solution.",
            )
            process_output.ocmmnt(
                constants.iotty,
                "The code may be stuck in a minimum in the residual space that is significantly above zero.",
            )
            process_output.oblnkl(constants.nout)
            process_output.oblnkl(constants.iotty)
            process_output.ocmmnt(
                constants.iotty, "There is either no solution possible, or the code"
            )
            process_output.ocmmnt(
                constants.iotty, "is failing to escape from a deep local minimum."
            )
            process_output.ocmmnt(
                constants.iotty,
                "Try changing the variables in IXC, or modify their initial values.",
            )
        elif ifail == 3:
            process_output.ocmmnt(
                constants.nout, "The line search required the maximum of 10 calls."
            )
            process_output.ocmmnt(
                constants.nout, "A feasible solution may be difficult to achieve."
            )
            process_output.ocmmnt(
                constants.nout, "Try changing or adding variables to IXC."
            )

            process_output.ocmmnt(
                constants.iotty, "The line search required the maximum of 10 calls."
            )
            process_output.ocmmnt(
                constants.iotty, "A feasible solution may be difficult to achieve."
            )
            process_output.ocmmnt(
                constants.iotty, "Try changing or adding variables to IXC."
            )
        elif ifail == 4:
            process_output.ocmmnt(
                constants.nout, "An uphill search direction was found."
            )
            process_output.ocmmnt(
                constants.nout, "Try changing the equations in ICC, or"
            )
            process_output.ocmmnt(constants.nout, "adding new variables to IXC.")

            process_output.ocmmnt(
                constants.iotty, "An uphill search direction was found."
            )
            process_output.ocmmnt(
                constants.iotty, "Try changing the equations in ICC, or"
            )
            process_output.ocmmnt(constants.iotty, "adding new variables to IXC.")
        elif ifail == 5:
            process_output.ocmmnt(
                constants.nout, "The quadratic programming technique was unable to"
            )
            process_output.ocmmnt(constants.nout, "find a feasible point.")
            process_output.oblnkl(constants.nout)
            process_output.ocmmnt(
                constants.nout, "Try changing or adding variables to IXC, or modify"
            )
            process_output.ocmmnt(
                constants.nout,
                "their initial values (especially if only 1 optimisation",
            )
            process_output.ocmmnt(constants.nout, "iteration was performed).")

            process_output.ocmmnt(
                constants.iotty, "The quadratic programming technique was unable to"
            )
            process_output.ocmmnt(constants.iotty, "find a feasible point.")
            process_output.oblnkl(constants.iotty)
            process_output.ocmmnt(
                constants.iotty, "Try changing or adding variables to IXC, or modify"
            )
            process_output.ocmmnt(
                constants.iotty,
                "their initial values (especially if only 1 optimisation",
            )
            process_output.ocmmnt(constants.iotty, "iteration was performed).")
        elif ifail == 6:
            process_output.ocmmnt(
                constants.nout, "The quadratic programming technique was restricted"
            )
            process_output.ocmmnt(
                constants.nout, "by an artificial bound, or failed due to a singular"
            )
            process_output.ocmmnt(constants.nout, "matrix.")
            process_output.ocmmnt(
                constants.nout, "Try changing the equations in ICC, or"
            )
            process_output.ocmmnt(constants.nout, "adding new variables to IXC.")

            process_output.ocmmnt(
                constants.iotty, "The quadratic programming technique was restricted"
            )
            process_output.ocmmnt(
                constants.iotty, "by an artificial bound, or failed due to a singular"
            )
            process_output.ocmmnt(constants.iotty, "matrix.")
            process_output.ocmmnt(
                constants.iotty, "Try changing the equations in ICC, or"
            )
            process_output.ocmmnt(constants.iotty, "adding new variables to IXC.")

    def scan_1d(self):
        """Run a 1-D scan."""
        # initialise dict which will contain ifail values for each scan point
        scan_1d_ifail_dict = {}

        for iscan in range(1, scan_module.isweep + 1):
            self.scan_1d_write_point_header(iscan)
            ifail = self.doopt()
            scan_1d_ifail_dict[iscan] = ifail
            write_output_files(models=self.models, ifail=ifail)

            error_handling.show_errors()
            error_handling.init_error_handling()

        # outvar now contains results
        self.scan_1d_write_plot()
        print(
            " ****************************************** Scan Convergence Summary ****************************************** \n"
        )
        sweep_values = scan_module.sweep[: scan_module.isweep]
        nsweep_var_name, _ = self.scan_select(
            scan_module.nsweep, scan_module.sweep, scan_module.isweep
        )
        converged_count = 0
        # offsets for aligning the converged/unconverged column
        max_sweep_value_length = len(str(np.max(sweep_values)).replace(".", ""))
        offsets = [
            max_sweep_value_length - len(str(sweep_val).replace(".", ""))
            for sweep_val in sweep_values
        ]
        for iscan in range(1, scan_module.isweep + 1):
            if scan_1d_ifail_dict[iscan] == 1:
                converged_count += 1
                print(
                    f"Scan {iscan:02d}: {nsweep_var_name} = {sweep_values[iscan - 1]} "
                    + " " * offsets[iscan - 1]
                    + "\u001b[32mCONVERGED \u001b[0m"
                )
            else:
                print(
                    f"Scan {iscan:02d}: {nsweep_var_name} = {sweep_values[iscan - 1]} "
                    + " " * offsets[iscan - 1]
                    + "\u001b[31mUNCONVERGED \u001b[0m"
                )
        converged_percentage = converged_count / scan_module.isweep * 100
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")

    def scan_2d(self):
        """Run a 2-D scan."""
        # Initialise intent(out) arrays
        self.scan_2d_init()
        iscan = 1

        # initialise array which will contain ifail values for each scan point
        scan_2d_ifail_list = np.zeros(
            (scan_module.noutvars, scan_module.ipnscns), dtype=np.float64, order="F"
        )
        for iscan_1 in range(1, scan_module.isweep + 1):
            for iscan_2 in range(1, scan_module.isweep_2 + 1):
                self.scan_2d_write_point_header(iscan, iscan_1, iscan_2)
                ifail = self.doopt()

                write_output_files(models=self.models, ifail=ifail)

                error_handling.show_errors()
                error_handling.init_error_handling()
                scan_2d_ifail_list[iscan_1][iscan_2] = ifail
                iscan = iscan + 1

        print(
            " ****************************************** Scan Convergence Summary ****************************************** \n"
        )
        sweep_1_values = scan_module.sweep[: scan_module.isweep]
        sweep_2_values = scan_module.sweep_2[: scan_module.isweep_2]
        nsweep_var_name, _ = self.scan_select(
            scan_module.nsweep, scan_module.sweep, scan_module.isweep
        )
        nsweep_2_var_name, _ = self.scan_select(
            scan_module.nsweep_2, scan_module.sweep_2, scan_module.isweep_2
        )
        converged_count = 0
        scan_point = 1
        # offsets for aligning the converged/unconverged column
        max_sweep1_value_length = len(str(np.max(sweep_1_values)).replace(".", ""))
        max_sweep2_value_length = len(str(np.max(sweep_2_values)).replace(".", ""))
        offsets = np.zeros(
            (scan_module.isweep, scan_module.isweep_2), dtype=int, order="F"
        )
        for count1, sweep1 in enumerate(sweep_1_values):
            for count2, sweep2 in enumerate(sweep_2_values):
                offsets[count1][count2] = (
                    max_sweep1_value_length
                    - len(str(sweep1).replace(".", ""))
                    + max_sweep2_value_length
                    - len(str(sweep2).replace(".", ""))
                )

        for iscan_1 in range(1, scan_module.isweep + 1):
            for iscan_2 in range(1, scan_module.isweep_2 + 1):
                if scan_2d_ifail_list[iscan_1][iscan_2] == 1:
                    converged_count += 1
                    print(
                        f"Scan {scan_point:02d}: ({nsweep_var_name} = {sweep_1_values[iscan_1 - 1]}, {nsweep_2_var_name} = {sweep_2_values[iscan_2 - 1]}) "
                        + " " * offsets[iscan_1 - 1][iscan_2 - 1]
                        + "\u001b[32mCONVERGED \u001b[0m"
                    )
                    scan_point += 1
                else:
                    print(
                        f"Scan {scan_point:02d}: ({nsweep_var_name} = {sweep_1_values[iscan_1 - 1]}, {nsweep_2_var_name} = {sweep_2_values[iscan_2 - 1]}) "
                        + " " * offsets[iscan_1 - 1][iscan_2 - 1]
                        + "\u001b[31mUNCONVERGED \u001b[0m"
                    )
                    scan_point += 1
        converged_percentage = (
            converged_count / (scan_module.isweep * scan_module.isweep_2) * 100
        )
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")

    def scan_2d_init(self):
        process_output.ovarin(
            constants.mfile,
            "Number of first variable scan points",
            "(isweep)",
            scan_module.isweep,
        )
        process_output.ovarin(
            constants.mfile,
            "Number of second variable scan points",
            "(isweep_2)",
            scan_module.isweep_2,
        )
        process_output.ovarin(
            constants.mfile,
            "Scanning first variable number",
            "(nsweep)",
            scan_module.nsweep,
        )
        process_output.ovarin(
            constants.mfile,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_module.nsweep_2,
        )
        process_output.ovarin(
            constants.mfile,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_module.nsweep_2,
        )
        process_output.ovarin(
            constants.mfile,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_module.nsweep_2,
        )

    def scan_1d_write_point_header(self, iscan: int):
        global_variables.iscan_global = iscan
        global_variables.vlabel, global_variables.xlabel = self.scan_select(
            scan_module.nsweep, scan_module.sweep, iscan
        )

        process_output.oblnkl(constants.nout)
        process_output.ostars(constants.nout, 110)

        process_output.write(
            constants.nout,
            f"***** Scan point {iscan} of {scan_module.isweep} : {f2py_compatible_to_string(global_variables.xlabel)}"
            f", {f2py_compatible_to_string(global_variables.vlabel)} = {scan_module.sweep[iscan - 1]} "
            "*****",
        )
        process_output.ostars(constants.nout, 110)
        process_output.oblnkl(constants.mfile)
        process_output.ovarin(constants.mfile, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan} of {scan_module.isweep} : "
            f"{f2py_compatible_to_string(global_variables.xlabel)} , {f2py_compatible_to_string(global_variables.vlabel)}"
            f" = {scan_module.sweep[iscan - 1]}"
        )

    def scan_2d_write_point_header(self, iscan, iscan_1, iscan_2):
        iscan_r = scan_module.isweep_2 - iscan_2 + 1 if iscan_1 % 2 == 0 else iscan_2

        # Makes iscan available globally (read-only)
        global_variables.iscan_global = iscan

        global_variables.vlabel, global_variables.xlabel = self.scan_select(
            scan_module.nsweep, scan_module.sweep, iscan_1
        )
        global_variables.vlabel_2, global_variables.xlabel_2 = self.scan_select(
            scan_module.nsweep_2, scan_module.sweep_2, iscan_r
        )

        process_output.oblnkl(constants.nout)
        process_output.ostars(constants.nout, 110)

        process_output.write(
            constants.nout,
            f"***** 2D Scan point {iscan} of {scan_module.isweep * scan_module.isweep_2} : "
            f"{f2py_compatible_to_string(global_variables.vlabel)} = {scan_module.sweep[iscan - 1]} and"
            f" {f2py_compatible_to_string(global_variables.vlabel_2)} = {scan_module.sweep_2[iscan_r]} "
            "*****",
        )
        process_output.ostars(constants.nout, 110)
        process_output.oblnkl(constants.mfile)
        process_output.ovarin(constants.mfile, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan}:  {f2py_compatible_to_string(global_variables.xlabel)}, "
            f"{f2py_compatible_to_string(global_variables.vlabel)} = {scan_module.sweep[iscan - 1]}"
            f" and {f2py_compatible_to_string(global_variables.xlabel_2)}, "
            f"{f2py_compatible_to_string(global_variables.vlabel_2)} = {scan_module.sweep_2[iscan_r]} "
        )

        return iscan_r

    def scan_1d_write_plot(self):
        if scan_module.first_call_1d:
            process_output.ovarin(
                constants.mfile, "Number of scan points", "(isweep)", scan_module.isweep
            )
            process_output.ovarin(
                constants.mfile,
                "Scanning variable number",
                "(nsweep)",
                scan_module.nsweep,
            )

            scan_module.first_call_1d = False

    def scan_select(self, nwp, swp, iscn):
        match nwp:
            case 1:
                physics_variables.aspect = swp[iscn - 1]
            case 2:
                divertor_variables.pflux_div_heat_load_max_mw = swp[iscn - 1]
            case 3:
                constraint_variables.pnetelin = swp[iscn - 1]
            case 4:
                physics_variables.hfact = swp[iscn - 1]
            case 5:
                tfcoil_variables.oacdcp = swp[iscn - 1]
            case 6:
                constraint_variables.walalw = swp[iscn - 1]
            case 7:
                physics_variables.beamfus0 = swp[iscn - 1]
            case 8:
                constraint_variables.fqval = swp[iscn - 1]
            case 9:
                physics_variables.te = swp[iscn - 1]
            case 10:
                numerics.boundu[14] = swp[iscn - 1]
            case 11:
                physics_variables.beta_norm_max = swp[iscn - 1]
            case 12:
                current_drive_variables.f_c_plasma_bootstrap_max = swp[iscn - 1]
            case 13:
                numerics.boundu[9] = swp[iscn - 1]
            case 14:
                constraint_variables.fiooic = swp[iscn - 1]
            case 15:
                constraint_variables.fjprot = swp[iscn - 1]
            case 16:
                physics_variables.rmajor = swp[iscn - 1]
            case 17:
                constraint_variables.bmxlim = swp[iscn - 1]
            case 18:
                constraint_variables.gammax = swp[iscn - 1]
            case 19:
                numerics.boundl[15] = swp[iscn - 1]
            case 20:
                constraint_variables.t_burn_min = swp[iscn - 1]
            case 22:
                if cost_variables.iavail == 1:
                    raise ProcessValueError("Do not scan cfactr if iavail=1")
                cost_variables.cfactr = swp[iscn - 1]
            case 23:
                numerics.boundu[71] = swp[iscn - 1]
            case 24:
                constraint_variables.p_fusion_total_max_mw = swp[iscn - 1]
            case 25:
                physics_variables.kappa = swp[iscn - 1]
            case 26:
                physics_variables.triang = swp[iscn - 1]
            case 27:
                constraint_variables.tbrmin = swp[iscn - 1]
            case 28:
                physics_variables.bt = swp[iscn - 1]
            case 29:
                impurity_radiation_module.coreradius = swp[iscn - 1]
            case 31:
                constraint_variables.f_alpha_energy_confinement_min = swp[iscn - 1]
            case 32:
                numerics.epsvmc = swp[iscn - 1]
            case 38:
                numerics.boundu[128] = swp[iscn - 1]
            case 39:
                numerics.boundu[130] = swp[iscn - 1]
            case 40:
                numerics.boundu[134] = swp[iscn - 1]
            case 41:
                build_variables.dr_blkt_outboard = swp[iscn - 1]
            case 42:
                impurity_radiation_module.fimp[8] = swp[iscn - 1]
                impurity_radiation_module.impurity_arr_frac[8] = (
                    impurity_radiation_module.fimp[8]
                )
            case 44:
                tfcoil_variables.sig_tf_case_max = swp[iscn - 1]
            case 45:
                tfcoil_variables.tmargmin_tf = swp[iscn - 1]
            case 46:
                numerics.boundu[151] = swp[iscn - 1]
            case 48:
                tfcoil_variables.n_pancake = int(swp[iscn - 1])
            case 49:
                tfcoil_variables.n_layer = int(swp[iscn - 1])
            case 50:
                impurity_radiation_module.fimp[12] = swp[iscn - 1]
                impurity_radiation_module.impurity_arr_frac[12] = (
                    impurity_radiation_module.fimp[12]
                )
            case 51:
                physics_variables.f_p_div_lower = swp[iscn - 1]
            case 52:
                physics_variables.rad_fraction_sol = swp[iscn - 1]
            case 53:
                numerics.boundu[156] = swp[iscn - 1]
            case 54:
                tfcoil_variables.b_crit_upper_nbti = swp[iscn - 1]
            case 55:
                build_variables.dr_shld_inboard = swp[iscn - 1]
            case 56:
                heat_transport_variables.crypmw_max = swp[iscn - 1]
            case 57:
                numerics.boundl[1] = swp[iscn - 1]
            case 58:
                build_variables.dr_fw_plasma_gap_inboard = swp[iscn - 1]
            case 59:
                build_variables.dr_fw_plasma_gap_outboard = swp[iscn - 1]
            case 60:
                tfcoil_variables.sig_tf_wp_max = swp[iscn - 1]
            case 61:
                rebco_variables.copperaoh_m2_max = swp[iscn - 1]
            case 62:
                pfcoil_variables.coheof = swp[iscn - 1]
            case 63:
                build_variables.dr_cs = swp[iscn - 1]
            case 64:
                pfcoil_variables.ohhghf = swp[iscn - 1]
            case 65:
                cs_fatigue_variables.n_cycle_min = swp[iscn - 1]
            case 66:
                pfcoil_variables.oh_steel_frac = swp[iscn - 1]
            case 67:
                cs_fatigue_variables.t_crack_vertical = swp[iscn - 1]
            case 68:
                fwbs_variables.inlet_temp_liq = swp[iscn - 1]
            case 69:
                fwbs_variables.outlet_temp_liq = swp[iscn - 1]
            case 70:
                fwbs_variables.blpressure_liq = swp[iscn - 1]
            case 71:
                fwbs_variables.n_liq_recirc = swp[iscn - 1]
            case 72:
                fwbs_variables.bz_channel_conduct_liq = swp[iscn - 1]
            case 73:
                fwbs_variables.pnuc_fw_ratio_dcll = swp[iscn - 1]
            case 74:
                fwbs_variables.f_nuc_pow_bz_struct = swp[iscn - 1]
            case 75:
                fwbs_variables.dx_fw_module = swp[iscn - 1]
            case 76:
                heat_transport_variables.etath = swp[iscn - 1]
            case 77:
                cost_variables.startupratio = swp[iscn - 1]
            case 78:
                cost_variables.fkind = swp[iscn - 1]
            case 79:
                current_drive_variables.eta_ecrh_injector_wall_plug = swp[iscn - 1]
            case 80:
                tfcoil_variables.fcoolcp = swp[iscn - 1]
            case 81:
                tfcoil_variables.n_tf_turn = swp[iscn - 1]
            case _:
                raise ProcessValueError("Illegal scan variable number", nwp=nwp)

        return SCAN_VARIABLES[int(nwp)]


def init_scan_module():
    """Initialise the scan module"""
    scan_module.scan_dim = 1
    scan_module.isweep = 0
    scan_module.isweep_2 = 0
    scan_module.nsweep = 1
    scan_module.nsweep_2 = 3
    scan_module.sweep[:] = 0.0
    scan_module.sweep_2[:] = 0.0
    scan_module.first_call_1d = True
    scan_module.first_call_2d = True
