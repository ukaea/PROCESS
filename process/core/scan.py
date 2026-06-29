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
from process.models.availability import AvailabilityModel

if TYPE_CHECKING:
    from process.core.model import DataStructure
    from process.main import Models


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
    t_crack_vertical = (67, Area.CS)
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


@dataclass
class ScanRes:
    iscan: int
    ifail: int
    solver: SolverHandler


class Scan:
    """Perform a parameter scan

    Parameters
    ----------
    models :
        Physics and engineering model objects
    solver :
        Which solver to use, as specified in solver.py
    data :
        Data structure object
    """

    def __init__(self, models: Models, solver: str, data: DataStructure):
        self.models = models
        self.solver = solver
        self.data = data

    def _run(self, iscan, nsweep, sweep, data):
        """
        This method calls the optimisation routine VMCON a number of times, by
        performing a sweep over a range of values of a particular variable. A
        number of output variable values are written to the MFILE.DAT file at
        each scan point, for plotting or other post-processing purposes.

        Raises
        ------
        ProcessValueError
            isweep value greater than IPNSCNS
        """
        sh = SolverHandler(self.models, self.solver, data)
        # TODO queue the output to avoid race condition (?)
        if data.scan.nsweep is not None:
            self.write_point_header(iscan)
        start_time = time.time()
        ifail = sh.run()
        end_time = time.time() - start_time
        write_output_files(models=self.models, data=data, ifail=ifail, runtime=end_time)
        nums = data.numerics
        nums.sqsumsq = sum(r**2 for r in nums.rcm[: nums.neqns]) ** 0.5

        if self.data.scan.scan_dim == 2:
            self.scan_2d()
        else:
            self.scan_1d()

    def write_outputs(self):
        write_output_files(
            models=self.models,
            data=self.data,
            ifail=self._ifail,
            runtime=self._finish_time - self._start_time,
        )
        show_errors(constants.NOUT)

        logging_model_handler.clear_logs()
        optimisation_output(data)
        constraints.constraints_output(data, self.solver)

        return ScanRes(iscan, ifail, sh)

    def _set_v_x_label(self, iscan: list[int]):
        sv = [
            self.scan_select(self.data.scan.nsweep, self.data.scan.sweep, isc)
            for isc in iscan
        ]
        self.data.globals.vlabel = [s.fname for s in sv]
        self.data.globals.xlabel = [s.data.description for s in sv]

    def write_point_header(self, iscan):
        self._set_v_x_label(iscan)

        process_output.oblnkl(constants.NOUT)
        process_output.oblnkl(constants.MFILE)

        process_output.write(
            constants.NOUT,
            f"Scan point {iscan} of {np.prod(self.data.scan.isweep)} : \n".join(
                f"{v} = {self.data.scan.sweep[iscan[no] - 1]}"
                for no, v in enumerate(self.data.globals.vlabel)
            ),
        )
        process_output.ovarin(constants.MFILE, "Scan point number", "(iscan)", iscan)

        print(
            f"Starting scan point {iscan}:  {self.data.globals.xlabel}, \n".join(
                f"{v} = {self.data.scan.sweep[iscan[no] - 1]}"
                for no, v in enumerate(self.data.globals.vlabel)
            )
        )

    def scan_select(self, nsweep, sweep, iscan):
        sv = ScanVariables(nsweep)
        sv.set(self.data, sweep[iscan - 1])
        return sv

    def run(self):
        """Call a solver over a range of values of one of the variables.

        This method calls the optimisation routine VMCON a number of times, by
        performing a sweep over a range of values of a particular variable. A
        number of output variable values are written to the MFILE.DAT file at
        each scan point, for plotting or other post-processing purposes.
        """
        # vectorise running of self._run
        if self.data.scan.nsweep is not None:
            for d, n, v in (
                ("Number of scan points", "(isweep)", self.data.scan.isweep),
                ("Scanning variable number", "(nsweep)", self.data.scan.nsweep),
            ):
                process_output.ovarin(constants.MFILE, d, n, v)

        # TODO copy of self.data for each vectorised run (?)
        scan_res = np.vectorise(self._run)(
            self.data.scan.isweep, self.data.scan.nsweep, self.data.scan.sweep, self.data
        )

        if self.data.scan.nsweep is not None:
            self.summary(scan_res)

    def summary(self, scan_res):
        print("Scan Convergence Summary\n")
        sweep_values = self.data.scan.sweep
        nsweep_var = [ScanVariables(nsw) for nsw in self.data.scan.nsweep]

        conv_list = []
        converged_count = 0
        conv_str = "\u001b[3{}CONVERGED \u001b[0m"
        for no, sr in enumerate(scan_res):
            if sr.ifail == 1:
                converged_count += 1
                conv = conv_str.format("2m")
            else:
                conv = conv_str.format("1mUN")
            conv_list.append([
                "{sr.iscan:02d}",
                nsweep_var[no].fname,
                sweep_values[sr.iscan],
                conv,
            ])

        print(
            tabulate(conv_list, headers=["Iscan", "Sweep Var", "Sweep Val", "Converged"])
        )

        converged_percentage = converged_count / np.prod(self.data.scan.isweep) * 100
        print(f"\nConvergence Percentage: {converged_percentage:.2f}%")


def optimisation_output(data: DataStructure):
    nums = data.numerics

    written_warning = False

    # Output optimisation parameters
    solution_vector_table = []
    for i in range(nums.nvar):
        nums.xcs[i] = nums.xcm[i] * nums.scafc[i]

        name = nums.lablxc[nums.ixc[i] - 1]
        solution_vector_table.append([name, nums.xcs[i], nums.xcm[i]])

        xminn = 1.01 * nums.itv_scaled_lower_bounds[i]
        xmaxx = 0.99 * nums.itv_scaled_upper_bounds[i]

        # Write to output file if close to optimisation parameter bounds
        if nums.xcm[i] < xminn or nums.xcm[i] > xmaxx:
            if not written_warning:
                written_warning = True
                process_output.ocmmnt(
                    constants.NOUT,
                    (
                        "Certain operating limits have been reached,"
                        "\n as shown by the following optimisation parameters that are"
                        "\n at or near to the edge of their prescribed range :\n"
                    ),
                )

            xcval = nums.xcm[i] * nums.scafc[i]

            if nums.xcm[i] < xminn:
                location, bound = "below", "lower"
                bounds = nums.itv_scaled_lower_bounds
            else:
                location, bound = "above", "upper"
                bounds = nums.itv_scaled_upper_bounds
            process_output.write(
                constants.NOUT,
                f"   {name:<30}= {xcval} is at or {location} its {bound} bound:"
                f" {bounds[i] * nums.scafc[i]}",
            )

        xnorm = (
            1.0
            if nums.boundu[i] == nums.boundl[i]
            else min(
                max(
                    (nums.xcm[i] - nums.itv_scaled_lower_bounds[i])
                    / (
                        nums.itv_scaled_upper_bounds[i] - nums.itv_scaled_lower_bounds[i]
                    ),
                    0.0,
                ),
                1.0,
            )
        )

        # Write optimisation parameters to mfile
        for d, var, v in (
            (nums.lablxc[nums.ixc[i] - 1], f"(itvar{i + 1:03d})", nums.xcs[i]),
            (f"{name} (final value/initial value)", f"(xcm{i + 1:03d})", nums.xcm[i]),
            (f"{name} (range normalised)", f"(nitvar{i + 1:03d})", xnorm),
            (
                f"{name} (upper bound)",
                f"(boundu{i + 1:03d})",
                nums.itv_scaled_upper_bounds[i] * nums.scafc[i],
            ),
            (
                f"{name} (lower bound)",
                f"(boundl{i + 1:03d})",
                nums.itv_scaled_lower_bounds[i] * nums.scafc[i],
            ),
        ):
            process_output.ovarre(constants.MFILE, d, var, v)

    # Write optimisation parameter headings to output file
    process_output.osubhd(
        constants.NOUT, "The solution vector is comprised as follows :"
    )
    process_output.write(
        constants.NOUT,
        tabulate(
            solution_vector_table,
            headers=["", "Final value", "Final / initial"],
            numalign="left",
        ),
    )
