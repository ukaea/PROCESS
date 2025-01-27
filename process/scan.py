import numpy as np

from process.caller import write_output_files
from process.fortran import error_handling, numerics, scan_module
from process.optimiser import Optimiser


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
        # f90wrap requires that arrays output from Fortran subroutines are
        # defined in Python and passed in as an argument, in similar style to
        # an intent(inout) argument. They are modified, but not returned.
        # Initialise intent(out) array outvar
        outvar = np.zeros(
            (scan_module.noutvars, scan_module.ipnscns), dtype=np.float64, order="F"
        )

        # initialise dict which will contain ifail values for each scan point
        scan_1d_ifail_dict = {}

        for iscan in range(1, scan_module.isweep + 1):
            scan_module.scan_1d_write_point_header(iscan)
            ifail = self.doopt()
            scan_1d_ifail_dict[iscan] = ifail
            write_output_files(models=self.models, ifail=ifail)

            # outvar is an intent(out) of scan_1d_store_output()
            outvar = scan_module.scan_1d_store_output(
                iscan,
                ifail,
                scan_module.noutvars,
                scan_module.ipnscns,
            )
            error_handling.show_errors()
            error_handling.init_error_handling()

        # outvar now contains results
        scan_module.scan_1d_write_plot(iscan, outvar)
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
                iscan_R = scan_module.scan_2d_write_point_header(
                    iscan, iscan_1, iscan_2
                )
                ifail = self.doopt()

                write_output_files(models=self.models, ifail=ifail)

                outvar, sweep_1_vals, sweep_2_vals = scan_module.scan_2d_store_output(
                    ifail,
                    iscan_1,
                    iscan_R,
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
