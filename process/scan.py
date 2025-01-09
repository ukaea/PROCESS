import numpy as np

from process.caller import write_output_files
from process.fortran import (
    constants,
    error_handling,
    global_variables,
    numerics,
    process_output,
    scan_module,
)
from process.optimiser import Optimiser
from process.utilities.f2py_string_patch import f2py_compatible_to_string


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
            error_handling.idiags[1] = scan_module.isweep
            error_handling.idiags[2] = scan_module.ipnscns
            error_handling.report_error(94)

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
        scan_module.post_optimise(ifail)

        return ifail

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
        nsweep_var_name, _ = scan_module.scan_select(
            scan_module.nsweep, scan_module.sweep, scan_module.isweep
        )
        converged_count = 0
        nsweep_var_name = nsweep_var_name.decode("utf-8")
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
        outvar = np.zeros(
            (scan_module.noutvars, scan_module.ipnscns), dtype=np.float64, order="F"
        )
        sweep_1_vals = np.ndarray(scan_module.ipnscns, dtype=np.float64, order="F")
        sweep_2_vals = np.ndarray(scan_module.ipnscns, dtype=np.float64, order="F")

        scan_module.scan_2d_init()
        iscan = 1

        # initialise array which will contain ifail values for each scan point
        scan_2d_ifail_list = np.zeros(
            (scan_module.noutvars, scan_module.ipnscns), dtype=np.float64, order="F"
        )
        for iscan_1 in range(1, scan_module.isweep + 1):
            for iscan_2 in range(1, scan_module.isweep_2 + 1):
                iscan_r = self.scan_2d_write_point_header(iscan, iscan_1, iscan_2)
                ifail = self.doopt()

                write_output_files(models=self.models, ifail=ifail)

                outvar, sweep_1_vals, sweep_2_vals = scan_module.scan_2d_store_output(
                    ifail,
                    iscan_1,
                    iscan_r,
                    iscan,
                    scan_module.noutvars,
                    scan_module.ipnscns,
                )
                error_handling.show_errors()
                error_handling.init_error_handling()
                scan_2d_ifail_list[iscan_1][iscan_2] = ifail
                iscan = iscan + 1

        scan_module.scan_2d_write_plot(iscan, outvar, sweep_1_vals, sweep_2_vals)
        print(
            " ****************************************** Scan Convergence Summary ****************************************** \n"
        )
        sweep_1_values = scan_module.sweep[: scan_module.isweep]
        sweep_2_values = scan_module.sweep_2[: scan_module.isweep_2]
        nsweep_var_name, _ = scan_module.scan_select(
            scan_module.nsweep, scan_module.sweep, scan_module.isweep
        )
        nsweep_2_var_name, _ = scan_module.scan_select(
            scan_module.nsweep_2, scan_module.sweep_2, scan_module.isweep_2
        )
        converged_count = 0
        scan_point = 1
        nsweep_var_name = nsweep_var_name.decode("utf-8")
        nsweep_2_var_name = nsweep_2_var_name.decode("utf-8")
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

    def scan_1d_write_point_header(self, iscan: int):
        global_variables.iscan_global = iscan
        global_variables.vlabel, global_variables.xlabel = scan_module.scan_select(
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

        global_variables.vlabel, global_variables.xlabel = scan_module.scan_select(
            scan_module.nsweep, scan_module.sweep, iscan_1
        )
        global_variables.vlabel_2, global_variables.xlabel_2 = scan_module.scan_select(
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
            self.plabel = [
                "Ifail____________________",
                "Sqsumsq__________________",
                "Electric_cost_(mil/kwh)__",
                "Capital_cost_(mil/kwh)___",
                "Fuel_cost_(mil/kwh)______",
                "Operations_cost_(mil/kwh)",
                "Capital_cost_(millions)__",
                "Core_costs_(millions)____",
                "Direct_cost_(billions)___",
                "Major_Radius_(m)_________",
                "Aspect_Ratio_____________",
                "Plasma_Current_(MA)______",
                "B_Toroidal_Axis_(T)______",
                "B_total_on_axis_(T)______",
                "Safety_Factor____________",
                "q95_min_(zero_if_i_plasma_geometry=0)",
                "Beta_____________________",
                "Beta_Limit_______________",
                "Epsilon_Beta_Poloidal____",
                "Dens.weight_Te_(10keV)___",
                "Average_Dens_(10^20/m^3)_",
                "H-fact_Iter_Power________",
                "H-fact_Iter_Offset_______",
                "Fusion_Power_(MW)________",
                "nb_Fusion_Power_(MW)_____",
                "Wall_Load_(MW/m^2)_______",
                "Injection_Power_(MW)_____",
                "Inject_Pwr_Wall_Plug_(MW)",
                "Heating_Power_(MW)_______",
                "Current_Drive_(MW)_______",
                "Big_Q____________________",
                "Bootstrap_Fraction_______",
                "Neutral_Beam_Energy_(MeV)",
                "Divertor_Heat_(MW/m^2)___",
                "TF_coil_Power_(MW)_______",
                "TF_coil_weight_(kg)______",
                "vM_stress_in_TF_case_(Pa)",
                "J_TF_inboard_leg_(MA/m^2)",
                "Centrepost_max_T_(TART)__",
                "Res_TF_inbrd_leg_Pwr_(MW)",
                "Coolant_Fraction_Ctr.____",
                "C/P_coolant_radius_(m)___",
                "C/P_coolant_velocity(m/s)",
                "C/P_pump_power_(MW)______",
                "PF_coil_Power_(MW)_______",
                "PF_coil_weight_(kg)______",
                "Gross_Elect_Pwr_(MW)_____",
                "Net_electric_Pwr_(MW)____",
                "Recirculating_Fraction___",
                "Psep/R___________________",
                "Tot._radiation_power_(MW)",
                "First_wall_peak_temp_(K)_",
                "Cu_frac_TFC_conductor____",
                "Winding_pack_area_TFC(m2)",
                "Conductor_area_TFC_(m2)__",
                "Area_TF_inboard_leg_(m2)_",
                "Taup/taueff_lower_limit__",
                "Plasma_temp_at_sep_[keV]_",
                "SOL_density_at_OMP_______",
                "Power_through__separatrix",
                "neomp/nesep______________",
                "qtargettotal_____________",
                "Total_pressure_at_target_",
                "Temperature_at_target____",
                "Helium_fraction__________",
                "Momentum_loss_factor_____",
                "totalpowerlost_[W]_______",
                "H__concentration_________",
                "He_concentration_________",
                "Be_concentration_________",
                "C__concentration_________",
                "N__concentration_________",
                "O__concentration_________",
                "Ne_concentration_________",
                "Si_concentration_________",
                "Ar_concentration_________",
                "Fe_concentration_________",
                "Ni_concentration_________",
                "Kr_concentration_________",
                "Xe_concentration_________",
                "W__concentration_________",
                "teped____________________",
                "Max_field_on_TF_coil_____",
            ]
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
