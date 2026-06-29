from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from enum import Enum
from types import DynamicClassAttribute
from typing import TYPE_CHECKING

import numpy as np
from tabulate import tabulate

from process.core import constants, process_output
from process.core.caller import write_output_files
from process.core.exceptions import ProcessValueError
from process.core.log import logging_model_handler, show_errors
from process.core.solver import constraints
from process.core.solver.solver_handler import SolverHandler

if TYPE_CHECKING:
    from process.core.model import DataStructure
    from process.main import Models

logger = logging.getLogger(__name__)


@dataclass
class ScanVariable:
    number: int
    area: SVE = field(repr=False)


class SVE(Enum):
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
    @classmethod
    def _missing_(cls, value):
        if isinstance(value, int):
            for sv in cls:
                if sv.number == value:
                    return sv
        raise ProcessValueError("Illegal scan variable number", nwp=value)

    def fname(self):
        if "__" in self.name:
            return self.name.replace("__", "(") + ")"
        return self.name

    def set(self, data: DataStructure, sweep_val: float):
        var_area = getattr(data, self.area.value)

        if self.number == 22 and var_area.i_plant_availability == 1:
            raise ProcessValueError(
                "Do not scan f_t_plant_available if i_plant_availability=1"
            )

        if "__" in self.name:
            name, index = self.name.split("__")
            getattr(var_area, name)[int(index)] = sweep_val
            if name == "f_nd_impurity_electrons":
                var_area.f_nd_impurity_electron_array[int(index - 1)] = sweep_val
        else:
            setattr(var_area, self.name, sweep_val)
            name = self.name

        self._data_ = getattr(var_area, name)

    @DynamicClassAttribute
    def data(self):
        if hasattr(self, "_data_"):
            return self._data_
        raise ValueError("Data not available")

    def get_val(self, mfile, scan):
        return mfile.get(self.name, scan=scan)

    aspect = (1, SVE.P)
    pflux_div_heat_load_max_mw = (2, SVE.D)
    p_plant_electric_net_required_mw = (3, SVE.C)
    hfact = (4, SVE.P)
    j_tf_coil_full_area = (5, SVE.T)
    pflux_fw_neutron_max_mw = (6, SVE.C)
    beamfus0 = (7, SVE.P)
    temp_plasma_electron_vol_avg_kev = (9, SVE.P)
    boundu__14 = (10, SVE.NUM)
    beta_norm_max = (11, SVE.P)
    f_c_plasma_bootstrap_max = (12, SVE.CD)
    boundu__10 = (13, SVE.NUM)
    rmajor = (16, SVE.P)
    b_tf_inboard_max = (17, SVE.C)
    eta_cd_norm_hcd_primary_max = (18, SVE.C)
    boundl__16 = (19, SVE.NUM)
    t_burn_min = (20, SVE.C)
    f_t_plant_available = (22, SVE.CST)
    p_fusion_total_max_mw = (24, SVE.C)
    kappa = (25, SVE.P)
    triang = (26, SVE.P)
    tbrmin = (27, SVE.C)
    b_plasma_toroidal_on_axis = (28, SVE.P)
    coreradius = (29, SVE.IR)
    f_alpha_energy_confinement_min = (31, SVE.C)
    epsvmc = (32, SVE.NUM)
    boundu__129 = (38, SVE.NUM)
    boundu__131 = (39, SVE.NUM)
    boundu__135 = (40, SVE.NUM)
    dr_blkt_outboard = (41, SVE.B)
    f_nd_impurity_electrons__9 = (42, SVE.IR)
    sig_tf_case_max = (44, SVE.T)
    temp_tf_superconductor_margin_min = (45, SVE.T)
    boundu__152 = (46, SVE.NUM)
    n_tf_wp_pancakes = (48, SVE.T)
    n_tf_wp_layers = (49, SVE.T)
    f_nd_impurity_electrons__13 = (50, SVE.IR)
    f_p_div_lower = (51, SVE.P)
    rad_fraction_sol = (52, SVE.P)
    boundu__157 = (53, SVE.NUM)
    b_crit_upper_nbti = (54, SVE.T)
    dr_shld_inboard = (55, SVE.B)
    p_cryo_plant_electric_max_mw = (56, SVE.HT)
    boundl__2 = (57, SVE.NUM)
    dr_fw_plasma_gap_inboard = (58, SVE.B)
    dr_fw_plasma_gap_outboard = (59, SVE.B)
    sig_tf_wp_max = (60, SVE.T)
    copperaoh_m2_max = (61, SVE.TR)
    coheof = (62, SVE.PF)
    dr_cs = (63, SVE.B)
    ohhghf = (64, SVE.PF)
    n_cycle_min = (65, SVE.CS)
    oh_steel_frac = (66, SVE.PF)
    t_crack_vertical = (67, SVE.CS)
    inlet_temp_liq = (68, SVE.FWBS)
    outlet_temp_liq = (69, SVE.FWBS)
    blpressure_liq = (70, SVE.FWBS)
    n_liq_recirc = (71, SVE.FWBS)
    bz_channel_conduct_liq = (72, SVE.FWBS)
    pnuc_fw_ratio_dcll = (73, SVE.FWBS)
    f_nuc_pow_bz_struct = (74, SVE.FWBS)
    dx_fw_module = (75, SVE.FWBS)
    eta_turbine = (76, SVE.HT)
    startupratio = (77, SVE.CST)
    fkind = (78, SVE.CST)
    eta_ecrh_injector_wall_plug = (79, SVE.CD)
    fcoolcp = (80, SVE.T)
    n_tf_coil_turns = (81, SVE.T)


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
