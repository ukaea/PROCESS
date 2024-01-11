from process.fortran import error_handling
from process.fortran import scan_module
from process.fortran import numerics
from process.optimiser import Optimiser
from process import final
import numpy as np


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
            final.finalise(self.models, ifail)
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
            return

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

        for iscan in range(1, scan_module.isweep + 1):
            scan_module.scan_1d_write_point_header(iscan)
            ifail = self.doopt()

            final.finalise(self.models, ifail)

            # outvar is an intent(out) of scan_1d_store_output()
            outvar = scan_module.scan_1d_store_output(
                iscan,
                ifail,
                scan_module.noutvars,
                scan_module.ipnscns,
            )

        # outvar now contains results
        scan_module.scan_1d_write_plot(iscan, outvar)

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

        for iscan_1 in range(1, scan_module.isweep + 1):
            for iscan_2 in range(1, scan_module.isweep_2 + 1):
                iscan_R = scan_module.scan_2d_write_point_header(
                    iscan, iscan_1, iscan_2
                )
                ifail = self.doopt()

                final.finalise(self.models, ifail)

                outvar, sweep_1_vals, sweep_2_vals = scan_module.scan_2d_store_output(
                    ifail,
                    iscan_1,
                    iscan_R,
                    iscan,
                    scan_module.noutvars,
                    scan_module.ipnscns,
                )

                iscan = iscan + 1

        scan_module.scan_2d_write_plot(iscan, outvar, sweep_1_vals, sweep_2_vals)
