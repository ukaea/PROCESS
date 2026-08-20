"""Scanning mechanics"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from enum import Enum
from types import DynamicClassAttribute
from typing import TYPE_CHECKING

import numpy as np

from process.core import constants, process_output
from process.core.caller import write_output_files
from process.core.exceptions import ProcessValueError
from process.core.io.data_structure_dicts import get_dicts
from process.core.log import logging_model_handler, show_errors
from process.core.solver import constraints
from process.core.solver.solver_handler import SolverHandler
from process.data_structure.scan_variables import IPNSCNS, NOUTVARS, ScanData
from process.models.availability import AvailabilityModel

if TYPE_CHECKING:
    from process.core.model import DataStructure, Model


logger = logging.getLogger(__name__)


@dataclass
class ScanVariable:
    """Scan variable container"""

    number: int
    area: Area = field(repr=False)
    _out_name_: str | None = None


class Area(Enum):
    """Scan variable data structure area shorthand

    this mirrors data_structure.area
    """

    P = "physics"
    D = "divertor"
    C = "constraints"
    T = "tfcoil"
    TR = "rebco"
    CD = "current_drive"
    NUM = "numerics"
    CST = "costs"
    IR = "impurity_radiation"
    B = "build"
    HT = "heat_transport"
    PF = "pf_coil"
    CS = "cs_fatigue"
    FWBS = "fwbs"


class ScanVariables(ScanVariable, Enum):
    """Scan variable options"""

    @classmethod
    def _missing_(cls, value):
        if isinstance(value, int):
            for sv in cls:
                if sv.number == value:
                    return sv
            raise ProcessValueError("Illegal scan variable number", nwp=value)

        if isinstance(value, str):
            return cls[value.replace("(", "__").replace(")", "")]

        return super()._missing_(value)

    @DynamicClassAttribute
    def fname(self):
        """Full name"""
        if "__" in self.name:
            return self.name.replace("__", "(") + ")"
        return self.name

    @DynamicClassAttribute
    def out_name(self):
        """Output name"""
        if self._out_name_ is None:
            return self.fname
        return self._out_name_

    def set(self, data: DataStructure, sweep_val: float):
        """Set value of scan variable

        Raises
        ------
        ProcessValueError
            Scanning f_t_plant_available if i_plant_availability=1"

        """
        var_area = getattr(data, self.area.value)

        if (
            self.number == 22
            and AvailabilityModel(var_area.i_plant_availability)
            != AvailabilityModel.USER_INPUT
        ):
            raise ProcessValueError(
                "Do not scan f_t_plant_available if i_plant_availability=1"
            )

        if "__" in self.name:
            name, index = self.name.split("__")
            getattr(var_area, name)[int(index) - 1] = sweep_val
            if name == "f_nd_impurity_electrons":
                var_area.f_nd_impurity_electron_array[int(index - 1)] = sweep_val
        else:
            setattr(var_area, self.name, sweep_val)
            name = self.name

        self._data_ = getattr(var_area, name)
        self._description_ = get_dicts()["DICT_DESCRIPTIONS"][name]

    @DynamicClassAttribute
    def data(self):
        """Variable value

        Raises
        ------
        ValueError
            data not set
        """
        if hasattr(self, "_data_"):
            return self._data_
        raise ValueError("Data not available")

    @DynamicClassAttribute
    def description(self):
        """Variable description

        Raises
        ------
        ValueError
            description not set
        """
        if hasattr(self, "_description_"):
            return self._description_
        raise ValueError("Description not available")

    def get_val(self, mfile, scan):
        """Get value from mfile"""
        # TODO this will fail for boundu/l we should write the scan variable to
        # the mfile and use that directly (also replacign output names)
        return mfile.get(self.out_name, scan=scan)

    aspect = (1, Area.P)
    pflux_div_heat_load_max_mw = (2, Area.D)
    p_plant_electric_net_required_mw = (3, Area.C)
    hfact = (4, Area.P)
    j_tf_coil_full_area = (5, Area.T)
    pflux_fw_neutron_max_mw = (6, Area.C)
    beamfus0 = (7, Area.P)
    temp_plasma_electron_vol_avg_kev = (9, Area.P)
    boundu__15 = (10, Area.NUM)
    beta_norm_max = (11, Area.P)
    f_c_plasma_bootstrap_max = (12, Area.CD)
    boundu__10 = (13, Area.NUM)
    f_j_tf_wp_critical_max = (14, Area.C)  # TODO is this needed
    rmajor = (16, Area.P)
    b_tf_inboard_max = (17, Area.C, "b_tf_inboard_peak_symmetric")
    eta_cd_norm_hcd_primary_max = (18, Area.C)
    boundl__16 = (19, Area.NUM)
    t_burn_min = (20, Area.C)
    f_t_plant_available = (22, Area.CST)
    p_fusion_total_max_mw = (24, Area.C)
    kappa = (25, Area.P)
    triang = (26, Area.P)
    tbrmin = (27, Area.C)
    b_plasma_toroidal_on_axis = (28, Area.P)
    coreradius = (29, Area.IR)
    f_alpha_energy_confinement_min = (31, Area.C)
    epsvmc = (32, Area.NUM)
    boundu__129 = (38, Area.NUM)
    boundu__131 = (39, Area.NUM)
    boundu__135 = (40, Area.NUM)
    dr_blkt_outboard = (41, Area.B)
    f_nd_impurity_electrons__9 = (42, Area.IR)
    sig_tf_case_max = (44, Area.T)
    temp_tf_superconductor_margin_min = (45, Area.T)
    boundu__152 = (46, Area.NUM)
    n_tf_wp_pancakes = (48, Area.T)
    n_tf_wp_layers = (49, Area.T)
    f_nd_impurity_electrons__13 = (50, Area.IR)
    f_p_div_lower = (51, Area.P)
    rad_fraction_sol = (52, Area.P)
    boundu__157 = (53, Area.NUM)
    b_crit_upper_nbti = (54, Area.T)
    dr_shld_inboard = (55, Area.B)
    p_cryo_plant_electric_max_mw = (56, Area.HT)
    boundl__2 = (57, Area.NUM)
    dr_fw_plasma_gap_inboard = (58, Area.B)
    dr_fw_plasma_gap_outboard = (59, Area.B)
    sig_tf_wp_max = (60, Area.T)
    copperaoh_m2_max = (61, Area.TR)
    coheof = (62, Area.PF)
    dr_cs = (63, Area.B)
    ohhghf = (64, Area.PF)
    n_cycle_min = (65, Area.CS)
    oh_steel_frac = (66, Area.PF)
    dz_cs_turn_crack_initial = (67, Area.CS)
    inlet_temp_liq = (68, Area.FWBS)
    outlet_temp_liq = (69, Area.FWBS)
    blpressure_liq = (70, Area.FWBS)
    n_liq_recirc = (71, Area.FWBS)
    bz_channel_conduct_liq = (72, Area.FWBS)
    pnuc_fw_ratio_dcll = (73, Area.FWBS)
    f_nuc_pow_bz_struct = (74, Area.FWBS)
    dx_fw_module = (75, Area.FWBS)
    eta_turbine = (76, Area.HT)
    startupratio = (77, Area.CST)
    fkind = (78, Area.CST)
    eta_ecrh_injector_wall_plug = (79, Area.CD)
    fcoolcp = (80, Area.T)
    n_tf_coil_turns = (81, Area.T)


class Scan:
    """Perform a parameter scan using the Fortran scan module."""

    def __init__(self, models: Model, solver: str, data: DataStructure):
        """Immediately run the run_scan() method.

        Parameters
        ----------
        models :
            Physics and engineering model objects
        solver :
            Which solver to use, as specified in solver.py
        data :
            Data structure object
        """
        self.models = models
        self.solver = solver
        self.data = data
        self.solver_handler = SolverHandler(models, solver, data)
        self.run_scan()

    def run_scan(self):
        """Call a solver over a range of values of one of the variables.

        This method calls the optimisation routine VMCON a number of times, by
        performing a sweep over a range of values of a particular variable. A
        number of output variable values are written to the MFILE.DAT file at
        each scan point, for plotting or other post-processing purposes.

        Raises
        ------
        ProcessValueError
            isweep value greater than IPNSCNS
        """
        if self.data.scan.isweep == 0:
            # Solve single problem, rather than an array of problems (scan)
            # doopt() can also run just an evaluation
            start_time = time.time()
            ifail = self.doopt()
            write_output_files(
                models=self.models,
                data=self.data,
                ifail=ifail,
                runtime=time.time() - start_time,
            )
            show_errors(constants.NOUT)
            return

        if self.data.scan.isweep > IPNSCNS:
            raise ProcessValueError(
                "Illegal value of isweep",
                isweep=self.data.scan.isweep,
                IPNSCNS=IPNSCNS,
            )

        if self.data.scan.scan_dim == 2:
            self.scan_2d()
        else:
            self.scan_1d()

    def doopt(self):
        """Run the optimiser or solver."""
        ifail = self.solver_handler.run()
        constraints.constraints_output(self.data, self.solver)

        return ifail

    def scan_1d(self):
        """Run a 1-D scan."""
        # initialise dict which will contain ifail values for each scan point
        scan_1d_ifail_dict = {}

        for iscan in range(1, self.data.scan.isweep + 1):
            self.scan_1d_write_point_header(iscan)
            start_time = time.time()
            ifail = self.doopt()
            scan_1d_ifail_dict[iscan] = ifail
            write_output_files(
                models=self.models,
                data=self.data,
                ifail=ifail,
                runtime=time.time() - start_time,
            )

            show_errors(constants.NOUT)
            logging_model_handler.clear_logs()

        # outvar now contains results
        self.scan_1d_write_plot(self.data.scan)
        print("Scan Convergence Summary \n")
        sweep_values = self.data.scan.sweep[: self.data.scan.isweep]
        nsweep_var = self.scan_select(
            self.data.scan.nsweep, self.data.scan.sweep, self.data.scan.isweep
        )
        converged_count = 0
        # offsets for aligning the converged/unconverged column
        max_sweep_value_length = len(str(np.max(sweep_values)).replace(".", ""))
        offsets = [
            max_sweep_value_length - len(str(sweep_val).replace(".", ""))
            for sweep_val in sweep_values
        ]
        for iscan in range(1, self.data.scan.isweep + 1):
            if scan_1d_ifail_dict[iscan] == 1:
                converged_count += 1
                print(
                    f"Scan {iscan:02d}: {nsweep_var.fname} = {sweep_values[iscan - 1]} "
                    + " " * offsets[iscan - 1]
                    + "\u001b[32mCONVERGED \u001b[0m"
                )
            else:
                print(
                    f"Scan {iscan:02d}: {nsweep_var.fname} = {sweep_values[iscan - 1]} "
                    + " " * offsets[iscan - 1]
                    + "\u001b[31mUNCONVERGED \u001b[0m"
                )
        converged_percentage = converged_count / self.data.scan.isweep * 100
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")

    def scan_2d(self):
        """Run a 2-D scan."""
        # Initialise intent(out) arrays
        self.scan_2d_init(self.data.scan)
        iscan = 1

        # initialise array which will contain ifail values for each scan point
        scan_2d_ifail_list = np.zeros(
            (NOUTVARS, IPNSCNS),
            dtype=np.float64,
            order="F",
        )
        for iscan_1 in range(1, self.data.scan.isweep + 1):
            for iscan_2 in range(1, self.data.scan.isweep_2 + 1):
                self.scan_2d_write_point_header(iscan, iscan_1, iscan_2)
                start_time = time.time()
                ifail = self.doopt()
                write_output_files(
                    models=self.models,
                    data=self.data,
                    ifail=ifail,
                    runtime=time.time() - start_time,
                )

                show_errors(constants.NOUT)
                logging_model_handler.clear_logs()
                scan_2d_ifail_list[iscan_1][iscan_2] = ifail
                iscan += 1

        print("Scan Convergence Summary\n")
        sweep_1_values = self.data.scan.sweep[: self.data.scan.isweep]
        sweep_2_values = self.data.scan.sweep_2[: self.data.scan.isweep_2]
        nsweep_var = self.scan_select(
            self.data.scan.nsweep, self.data.scan.sweep, self.data.scan.isweep
        )
        nsweep_2_var = self.scan_select(
            self.data.scan.nsweep_2, self.data.scan.sweep_2, self.data.scan.isweep_2
        )
        converged_count = 0
        scan_point = 1
        # offsets for aligning the converged/unconverged column
        max_sweep1_value_length = len(str(np.max(sweep_1_values)).replace(".", ""))
        max_sweep2_value_length = len(str(np.max(sweep_2_values)).replace(".", ""))
        offsets = np.zeros(
            (self.data.scan.isweep, self.data.scan.isweep_2), dtype=int, order="F"
        )
        for count1, sweep1 in enumerate(sweep_1_values):
            for count2, sweep2 in enumerate(sweep_2_values):
                offsets[count1][count2] = (
                    max_sweep1_value_length
                    - len(str(sweep1).replace(".", ""))
                    + max_sweep2_value_length
                    - len(str(sweep2).replace(".", ""))
                )

        for iscan_1 in range(1, self.data.scan.isweep + 1):
            for iscan_2 in range(1, self.data.scan.isweep_2 + 1):
                if scan_2d_ifail_list[iscan_1][iscan_2] == 1:
                    converged_count += 1
                    print(
                        (
                            f"Scan {scan_point:02d}: ({nsweep_var.fname} = "
                            f"{sweep_1_values[iscan_1 - 1]}, {nsweep_2_var.fname} "
                            f"= {sweep_2_values[iscan_2 - 1]}) "
                        )
                        + " " * offsets[iscan_1 - 1][iscan_2 - 1]
                        + "\u001b[32mCONVERGED \u001b[0m"
                    )
                    scan_point += 1
                else:
                    print(
                        (
                            f"Scan {scan_point:02d}: ({nsweep_var.fname} = "
                            f"{sweep_1_values[iscan_1 - 1]}, {nsweep_2_var.fname} = "
                            f"{sweep_2_values[iscan_2 - 1]}) "
                        )
                        + " " * offsets[iscan_1 - 1][iscan_2 - 1]
                        + "\u001b[31mUNCONVERGED \u001b[0m"
                    )
                    scan_point += 1
        converged_percentage = (
            converged_count / (self.data.scan.isweep * self.data.scan.isweep_2) * 100
        )
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")

    @staticmethod
    def scan_2d_init(scan_data: ScanData):
        """Scan 2d initialisation"""
        process_output.ovarre(
            constants.MFILE,
            "Number of first variable scan points",
            "(isweep)",
            scan_data.isweep,
        )
        process_output.ovarre(
            constants.MFILE,
            "Number of second variable scan points",
            "(isweep_2)",
            scan_data.isweep_2,
        )
        process_output.ovarre(
            constants.MFILE,
            "Scanning first variable number",
            "(nsweep)",
            scan_data.nsweep,
        )
        process_output.ovarre(
            constants.MFILE,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_data.nsweep_2,
        )
        process_output.ovarre(
            constants.MFILE,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_data.nsweep_2,
        )
        process_output.ovarre(
            constants.MFILE,
            "Scanning second variable number",
            "(nsweep_2)",
            scan_data.nsweep_2,
        )

    def scan_1d_write_point_header(self, iscan: int):
        """Scan 1d header"""
        self.data.globals.iscan_global = iscan
        sv = self.scan_select(self.data.scan.nsweep, self.data.scan.sweep, iscan)

        self.data.globals.vlabel = sv.fname
        self.data.globals.xlabel = sv.description

        process_output.oblnkl(constants.NOUT)
        process_output.ostars(constants.NOUT, 110)

        process_output.write(
            constants.NOUT,
            f"***** Scan point {iscan} of {self.data.scan.isweep} : "
            f"{self.data.globals.xlabel}"
            f", {self.data.globals.vlabel} = {self.data.scan.sweep[iscan - 1]} "
            "*****",
        )
        process_output.ostars(constants.NOUT, 110)
        process_output.oblnkl(constants.MFILE)
        process_output.ovarre(constants.MFILE, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan} of {self.data.scan.isweep} : "
            f"{self.data.globals.xlabel} , {self.data.globals.vlabel}"
            f" = {self.data.scan.sweep[iscan - 1]}"
        )

    def scan_2d_write_point_header(self, iscan, iscan_1, iscan_2):
        """Scan 2d header"""
        iscan_r = self.data.scan.isweep_2 - iscan_2 + 1 if iscan_1 % 2 == 0 else iscan_2

        # Makes iscan available globally (read-only)
        self.data.globals.iscan_global = iscan
        sv_1 = self.scan_select(self.data.scan.nsweep, self.data.scan.sweep, iscan_1)

        self.data.globals.vlabel = sv_1.fname
        self.data.globals.xlabel = sv_1.data.description

        sv_2 = self.scan_select(self.data.scan.nsweep_2, self.data.scan.sweep_2, iscan_r)

        self.data.globals.vlabel_2 = sv_2.fname
        self.data.globals.xlabel_2 = sv_2.data.description

        process_output.oblnkl(constants.NOUT)
        process_output.ostars(constants.NOUT, 110)

        process_output.write(
            constants.NOUT,
            f"***** 2D Scan point {iscan} of "
            f"{self.data.scan.isweep * self.data.scan.isweep_2} : "
            f"{self.data.globals.vlabel} = {self.data.scan.sweep[iscan_1 - 1]} and"
            f" {self.data.globals.vlabel_2} = {self.data.scan.sweep_2[iscan_r - 1]} "
            "*****",
        )
        process_output.ostars(constants.NOUT, 110)
        process_output.oblnkl(constants.MFILE)
        process_output.ovarre(constants.MFILE, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan}:  {self.data.globals.xlabel}, "
            f"{self.data.globals.vlabel} = {self.data.scan.sweep[iscan_1 - 1]}"
            f" and {self.data.globals.xlabel_2}, "
            f"{self.data.globals.vlabel_2} = {self.data.scan.sweep_2[iscan_r - 1]} "
        )

        return iscan_r

    @staticmethod
    def scan_1d_write_plot(scan_data: ScanData):
        """Scan 1d plotter"""
        if scan_data.first_call_1d:
            process_output.ovarre(
                constants.MFILE,
                "Number of scan points",
                "(isweep)",
                scan_data.isweep,
            )
            process_output.ovarre(
                constants.MFILE,
                "Scanning variable number",
                "(nsweep)",
                scan_data.nsweep,
            )

            scan_data.first_call_1d = False

    def scan_select(self, nsweep, sweep, iscan) -> ScanVariables:
        """Select a scan"""
        sv = ScanVariables(nsweep)
        sv.set(self.data, sweep[iscan - 1])
        return sv
