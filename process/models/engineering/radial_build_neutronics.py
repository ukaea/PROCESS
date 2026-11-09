# imports
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import openmc
from matproplib.conditions import OperationalConditions
from matproplib.library.fluids import Air, Water
from matproplib.library.steel import SS316_L
from matproplib.material import Material, material, mixture
from matproplib.nucleides import Elements, atomic_fraction_to_mass_fraction
from matproplib.properties.dependent import (
    Density,
)
from matproplib.properties.group import props
from tabulate import tabulate

from process.core import constants

# from process.data_structure import fwbs_variables as fwbs
from process.core import process_output as po
from process.core.repository import get_process_root

data_dir = get_process_root() / "data" / "process_openmc_data"


def check_isotope_files(folder, elements_list):
    """
    Checks whether there are nuclear data files corresponding with each element
    in the provided list. Ensure files in folder location use naming convention
    C012 rather than C12.

    If any desired elements aren't present in the folder it will add these to a
    json file that can be used to download just the missing elements. This will
    result in an error being produced until these files are downloaded to the given
    folder.
    """
    folderpath = Path(folder)
    tendl_path = folderpath / "tendl-2021-ace" / "tendl21c"
    return_list = []
    for el in elements_list:
        el_file = el.element.symbol + str(el.mass_number).zfill(3)
        filepath = tendl_path / el_file
        if not filepath.is_file():
            return_list += [el_file]
    if len(return_list) > 0:
        # remove duplicates
        return_list = list(set(return_list))
        element_dict = {"tendl": {"neutron": return_list}}
        file = folderpath / "missing_elemenets.json"
        with open(file, "w", encoding="utf-8") as fh:
            json.dump(element_dict, fh)
        raise ValueError(
            f"Required neutron data for material in radial build not found in {tendl_path}. {len(return_list)} instances"
            f" identified and listed in file: '{file}'. Run 'scripts/openmc_data_download/neutronics_data_downloader'"
            f" with '--isotope-file {file}'."
        )


def removal_xs_from_datafile(
    keyname, folder: str | None = None, energy: str = "fusion", factor=1.0
):
    """
    Obtains the total cross section of the given element and then calculates the removal
    cross section for it. Values can be produced at either fusion (14MeV) or fission (2MeV)
    energies.
    """
    if energy == "fusion":
        ev = 14e6
    elif energy == "fission":
        ev = 2e6
    else:
        raise ValueError(
            f"Invalid input for energy ({energy}). Must be either 'fusion' or 'fission'."
        )

    folderpath = Path(folder)
    tendl_path = folderpath / "tendl-2021-ace" / "tendl21c"
    data_file = keyname
    filepath = tendl_path / data_file
    if filepath.is_file():
        neutron_data = openmc.data.IncidentNeutron.from_ace(filepath)
        total = neutron_data.reactions[1].xs["294K"]
        if isinstance(total, openmc.data.function.ResonancesWithBackground):
            for res in total.resonances:
                if res.energy_min <= ev <= res.energy_max:
                    print(
                        f"Warning: {ev} is within the resonance of {res} found in {neutron_data}"
                    )
            total = total.background
        total_xs = total(x=ev)
        return soliman_cross_section_fit(total_xs, energy, factor)
    raise ValueError(f"Invalid filepath: {filepath}, no data found at this location.")


def fw_fluence_fpy(flux, energy):
    """
    Function to calculate the neutron fluence at the FW for 1 FPY.

    inputs
    ------
    flux:
        neutron flux on FW measured in MW/m2
    energy:
        neutron energy from fusion reaction in MeV

    Returns
    -------
    fluence:
        neutron fluence measured in n/cm2 at FW
    """
    neutron_flux_w_m2 = flux * 1e6
    neutron_energy_joules = energy * 1e6 * 1.6e-19
    fpy_seconds = 3600 * 24 * 365
    neutron_fluence_m2 = neutron_flux_w_m2 * fpy_seconds / neutron_energy_joules
    return neutron_fluence_m2 / 1e4  # convert to cm2 from m2


# needed functions
def top_hat(x, a, b):
    """
    Top hat function used for neutron fluence propagation.
    """
    alpha = np.min([a, b])
    beta = np.max([a, b])
    return np.heaviside(x - alpha, 1) - np.heaviside(x - beta, 0)


def component_mapping(label):
    """
    Mapping function to get component from the radial build label.
    """
    mapping_dict = {
        "VV": ["dr_vv_inboard", "dr_vv_outboard"],
        "TF": ["dr_tf_inboard", "dr_tf_outboard"],
        "TS": ["dr_shld_thermal_inboard", "dr_shld_thermal_outboard"],
        "RS": ["dr_shld_inboard", "dr_shld_outboard"],
        "BB": ["dr_blkt_inboard", "dr_blkt_outboard"],
        "FW": ["dr_fw_inboard", "dr_fw_outboard"],
        "CS": ["dr_cs", "dr_cs_precomp"],
        "Gap": [
            "dr_bore",
            "dr_cs_tf_gap",
            "dr_tf_shld_gap",
            "dr_shld_vv_gap_inboard",
            "dr_shld_blkt_gap",
            "dr_shld_vv_gap_outboard",
            "dr_tf_shld_gap",
        ],
    }
    for key, rb_list in mapping_dict.items():
        if label in rb_list:
            return key
    raise ValueError(f"Invalid input label: {label}. Not found in mapping dict.")


def _convert_micro_macro(sigma, a):
    return 0.6022 * sigma / a


def _convert_macro_micro(sigma_rho, a):
    return sigma_rho * a / 0.6022


def _material_elements_sigma(
    elements: Elements,
    path: str | None = None,
    input_factor: float = 1.0,
    model_factor: float = 1.0,
):
    elements_breakdown = {}
    for e_n in elements:
        a = elements[e_n].element.mass_number
        keyname = e_n + str(a).zfill(3)
        sigma = removal_xs_from_datafile(keyname, folder=path, factor=model_factor)
        sigma *= input_factor
        sigma = _convert_micro_macro(sigma, a)
        elements_breakdown[e_n] = {
            "fraction": elements[e_n].fraction,
            "sigma": sigma,
        }
    return elements_breakdown


def _compound_sigma(elements_dict):
    sigma_c = 0
    for el in elements_dict:
        sigma_c += elements_dict[el]["fraction"] * elements_dict[el]["sigma"]
    return sigma_c


def soliman_cross_section_fit(xs, energy: str = "fusion", factor=1.0):
    """
    Fit to estimate removal cross section from total cross section at fission or 14 MeV

    From A. Soliman 2025 paper: Energy-dependent neutron removal cross-section.
    """
    fission_params = [0.954, 2.435, 0.395, 2.137, 4.421, -0.385, -0.333]
    fission_params = [factor * i for i in fission_params]
    fusion_params = [0.542, 11.388, 0.077, 1.659, 7.979, -0.206, -1.43]
    fusion_params = [factor * i for i in fusion_params]
    if energy == "fusion":
        a = fusion_params[2] / (1 + np.exp((xs - fusion_params[3]) / fusion_params[5]))
        b = (1 - fusion_params[2]) / (
            1 + np.exp((xs - fusion_params[4]) / fusion_params[6])
        )
        return fusion_params[0] + fusion_params[1] * (a + b)
    if energy == "fission":
        a = fission_params[2] / (
            1 + np.exp((xs - fission_params[3]) / fission_params[5])
        )
        b = (1 - fission_params[2]) / (
            1 + np.exp((xs - fission_params[4]) / fission_params[6])
        )
        return fission_params[0] + fission_params[1] * (a + b)
    raise ValueError(
        f"Invalid energy input: {energy}. Should be either 'fusion' or 'fission'."
    )


def calculate_shielding_sigma_wrapper(
    compound: Material,
    path: str | None = None,
    input_factor: float = 1.0,
    model_factor: float = 1.0,
):
    """
    Function to calculate the mass removal cross section for the input material.

    Makes use of the materials elements properties (mass and atomic number) to
    approximate its microscopic cross section. This is converted to macroscopic,
    which is then used to calculate the compound cross section making use of the
    mass fraction weighting of each element.

    Works for matproplib materials and mixtures.
    """
    elements = atomic_fraction_to_mass_fraction(compound.elements.root)
    elements_breakdown = _material_elements_sigma(
        elements, path=path, input_factor=input_factor, model_factor=model_factor
    )
    return _compound_sigma(elements_breakdown)


def fluence_decay(r, r0, f0, mu):
    """
    Exponential decay that is cyclidrically corrected.

    Inputs are in m, converted to cm.
    """
    r *= 100
    r0 *= 100
    factor1 = (r0 / r) * f0
    factor2 = np.exp(-mu * np.abs(r - r0))
    return factor1 * factor2


def _default_material(component_name):
    steel_defaults = ["VV", "BB", "TF", "RS", "Steel"]
    air_defaults = ["Gap", "Air"]
    water_defaults = ["Water"]
    if component_name in steel_defaults:
        material = SS316_L()
    elif component_name in air_defaults:
        material = Air()
    elif component_name in water_defaults:
        material = Water()
    else:
        raise ValueError(f"Component name {component_name} has no recognised defaults")
    return material


def _build_mixture(component_name, mixture_dict):
    mixture_list = []
    for mat_item in mixture_dict:
        mat = _build_material(mat_item, mixture_dict[mat_item])
        mixture_list += [(mat, mixture_dict[mat_item]["fraction"])]
    return mixture(
        component_name,
        mixture_list,
        mix_condition=OperationalConditions(temperature=300, pressure=(1, "atm")),
        fraction_type="mass",
    )


def _build_material(component_name, material_dict):
    if material_dict["default"]:
        mat = _default_material(component_name=component_name)
    else:
        mat = material(
            component_name,
            elements=material_dict["elements"],
            properties=props(
                density=Density(
                    value=material_dict["density"]["value"],
                    unit=material_dict["density"]["unit"],
                )
            ),
        )()
    return mat


# class
class RadialBuildNeutronics:
    """
    class to handle merging of radial build inputs and component compositions to produce
    neutron fluence across the radial build
    """

    def __init__(
        self,
        ib: bool,
        component_dict: dict,
        incident_fluence: float | None = None,
        datapath: str = data_dir,
        input_factor: float = 0,
        model_factor: float = 0,
    ):
        self.datapath = datapath
        self.xs_input_factor = 1 + input_factor
        self.xs_model_factor = 1 + model_factor
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.ib = ib
        self.component_inputs = component_dict
        self.incident_fluence = incident_fluence
        self.op_cond = OperationalConditions(temperature=300, pressure=(1, "atm"))
        self.components, self.list_of_elements = self._construct_materials()
        check_isotope_files(data_dir, self.list_of_elements)
        self.component_sigmas = self._cross_section_wrapper()
        self.rb_fluence_dict = None

    def output(self):
        if self.rb_fluence_dict is None:
            self.run()
        # print(self.rb_fluence_dict)
        if self.ib:
            po.oheadr(self.outfile, "Inboard Radial Build Fluence")
            po.oblnkl(self.outfile)
            po.ovarin(
                self.outfile,
                "Fluence at inboard TF outer winding pack edge n/cm^2",
                "(tf_wp_fluence)",
                self.data.fwbs.tf_wp_fluence,
            )
            label = "ib"
            po.ovarre(
                self.mfile,
                "Number of inboard radial components",
                "(no_ib_comp)",
                f"{len(self.rb_fluence_dict)}",
            )
        else:
            po.oheadr(self.outfile, "Outboard Radial Build Fluence")
            label = "ob"
            po.ovarre(
                self.mfile,
                "Number of outboard radial components",
                "(no_ob_comp)",
                f"{len(self.rb_fluence_dict)}",
            )
        po.oblnkl(self.outfile)
        rows = []
        for i, key in enumerate(self.rb_fluence_dict):
            entry = self.rb_fluence_dict[key]
            if i == 0:
                r = entry["region"][0]
                f = entry["fluence"][0]
                sigma = entry["sigma"]
                rows += [
                    [
                        f"{self.initial_label}",
                        f"{sigma}",
                        f"{r}",
                        f"{np.format_float_scientific(f, precision=4)}",
                    ]
                ]
            r = entry["region"][1]
            f = entry["fluence"][1]
            sigma = entry["sigma"]
            rows += [
                [
                    f"{key}",
                    f"{sigma}",
                    f"{r}",
                    f"{np.format_float_scientific(f, precision=4)}",
                ]
            ]

        for line in tabulate(
            rows,
            headers=["Label", "Cross Section (/cm)", "Radius(m)", "Fluence(n/cm^2)"],
            tablefmt="plain",
        ).splitlines():
            po.write(self.outfile, line)

        for i, key in enumerate(self.rb_fluence_dict):
            entry = self.rb_fluence_dict[key]
            sigma = entry["sigma"]
            if i == 0:
                r = entry["region"][0]
                f = entry["fluence"][0]
                po.ovarre(
                    self.mfile,
                    f"Component {i} radius (m)",
                    f"({label}_r_comp{i})",
                    f"{r}",
                )
                po.ovarre(
                    self.mfile,
                    f"Component {i} cross section (/cm)",
                    f"({label}_xs_comp{i})",
                    f"{r}",
                )
                po.ovarre(
                    self.mfile,
                    f"Component {i} fluence (n/cm2)",
                    f"({label}_fluence_comp{i})",
                    f"{f:.6e}",
                )

            r = entry["region"][1]
            f = entry["fluence"][1]
            po.ovarre(
                self.mfile,
                f"Component {i} radius (m)",
                f"({label}_r_comp{i})",
                f"{r}",
            )
            po.ovarre(
                self.mfile,
                f"Component {i} cross section (/cm)",
                f"({label}_xs_comp{i})",
                f"{r}",
            )
            po.ovarre(
                self.mfile,
                f"Component {i} fluence (n/cm2)",
                f"({label}_fluence_comp{i})",
                f"{f:.6e}",
            )

    def run(self, radial_build=None, flux=None, test=False):
        if radial_build is None:
            radial_build = self.data.build.radial_build
        if self.incident_fluence is None:
            if flux is None:
                flux = self.data.physics.pflux_fw_neutron_mw
            self.incident_fluence = fw_fluence_fpy(flux, 14.1)
        if test:
            rb_labels, rb_values = radial_build
        else:
            rb_labels, rb_values = self.get_rb(radial_build=radial_build)
        self.initial_label = rb_labels[0]
        self.bounds = [rb_values[0], rb_values[-1]]
        self.rb_fluence_dict = self._radial_build_conversion(
            rb_labels[1:], rb_values[1:]
        )
        self._merge_component_inputs()
        self._calculate_fluence()
        if self.ib and not test:
            self._tf_wp_fluence()

    def _radial_build_conversion(self, labels, values):
        # print("conv", self.bounds[0], values)
        output_dict = {}
        cv = self.bounds[0]
        for label, value in zip(labels, values, strict=True):
            output_dict[label] = {
                "component": component_mapping(label),
                "region": [cv, value],
            }
            # print(label, cv, value)
            cv = value
        return output_dict

    def _construct_materials(self):
        els = []
        comp_dict = {}
        for comp in self.component_inputs:
            input_vals = self.component_inputs[comp]
            if input_vals["type"] == "cross_section":
                continue
            if input_vals["type"] == "material":
                mat = _build_material(comp, input_vals)
                comp_dict[comp] = mat
                elements = mat.elements.root
                for el in elements:
                    els += [elements[el].element]
            elif input_vals["type"] == "mixture":
                mix = _build_mixture(comp, input_vals["materials"])
                comp_dict[comp] = mix
                elements = mix.elements.root
                for el in elements:
                    els += [elements[el].element]
            elif input_vals["type"] != "cross_section":
                raise ValueError(f"Input type {input_vals['type']} is not recognised")
        els = list(els)
        return comp_dict, els

    def _cross_section_wrapper(self):
        sigma_dict = {}
        for comp in self.components:
            input_vals = self.component_inputs[comp]
            if input_vals["type"] == "cross_section":
                sigma_dict[comp] = input_vals["sigma"] * self.xs_factor
            elif input_vals["type"] in {"material", "mixture"}:
                sigma = calculate_shielding_sigma_wrapper(
                    self.components[comp],
                    input_factor=self.xs_input_factor,
                    path=self.datapath,
                    model_factor=self.xs_model_factor,
                ) * self.components[comp].density.value_as(self.op_cond, "g/cm^3")
                sigma_dict[comp] = sigma
            else:
                raise ValueError(f"Input type {input_vals['type']} is not recognised")
        return sigma_dict

    def _merge_component_inputs(self):
        for key in self.rb_fluence_dict:
            entry = self.rb_fluence_dict[key]
            component = entry["component"]
            self.rb_fluence_dict[key]["sigma"] = self.component_sigmas[component]

    def _calculate_fluence(self):
        f0 = self.incident_fluence
        for key in self.rb_fluence_dict:
            entry = self.rb_fluence_dict[key]
            f = fluence_decay(entry["region"][1], entry["region"][0], f0, entry["sigma"])
            self.rb_fluence_dict[key]["fluence"] = [f0, f]
            f0 = f

    def _tf_wp_fluence(self):
        radius = self.data.superconducting_tfcoil.r_tf_wp_inboard_outer
        self.data.fwbs.tf_wp_fluence = self.fluence_at_radius(radius)

    def fluence_at_radius(self, radius):
        if np.min(self.bounds) <= radius <= np.max(self.bounds):
            fluence = 0
            for key in self.rb_fluence_dict:
                entry = self.rb_fluence_dict[key]
                tval = top_hat(radius, entry["region"][0], entry["region"][1])
                fval = fluence_decay(
                    radius, entry["region"][0], entry["fluence"][0], entry["sigma"]
                )
                fluence += tval * fval
            return fluence
        raise ValueError(f"radius {radius} not in valid region {self.bounds}")

    def plot_fluence(self, n, ax=None, show=True, line="r-"):
        radii = np.linspace(self.bounds[0], self.bounds[1], n)
        f = []
        for r in radii:
            f += [self.fluence_at_radius(r)]
        if ax is None:
            _, ax = plt.subplots()
        ax.plot(radii, f, line)
        ax.plot(self.bounds[0], self.incident_fluence, "go")
        ax.vlines(self.bounds[0], 0, self.incident_fluence, colors="blue")
        for key in self.rb_fluence_dict:
            entry = self.rb_fluence_dict[key]
            ax.plot(entry["region"][1], entry["fluence"][1], "go")
            ax.vlines(entry["region"][1], 0, self.incident_fluence, colors="blue")
            ax.text(
                0.5 * (entry["region"][0] + entry["region"][1]),
                0.9 * self.incident_fluence,
                entry["component"],
            )
        ax.set_xlabel("Radius [m]")
        ax.set_ylabel("Fluence [n/cm^2]")
        ax.set_yscale("log")
        if show:
            plt.show()

    def get_rb(self, radial_build):
        """
        Function to extract key radial build labels and values in
        valid range for ib and ob for purpose of shielding.

        IB: dr_bore -> dr_fw_plasma_gap_inboard
        OB: dr_fw_plasma_gap_outboard -> dr_tf_outboard

        ib is bool where true returns IB labels and values whilst false returns OB.
        """
        if isinstance(radial_build, list):
            radial_build = np.array(radial_build)

        if self.ib:
            start_idx = 0  # centreline
            end_idx = (
                np.where(radial_build[:, 1] == "dr_fw_plasma_gap_inboard")[0][0] + 1
            )

        else:
            start_idx = np.where(radial_build[:, 1] == "dr_fw_plasma_gap_outboard")[0][0]
            end_idx = len(radial_build[:, 1])  # TF

        labels = radial_build[start_idx:end_idx, 1]
        values = radial_build[start_idx:end_idx, 3]
        if self.ib:
            if np.abs(values[0]) <= 1e-6:
                values[0] = 1e-6
            labels = labels[-1:0:-1]  # invert and remove final label
            values = values[-2::-1]  # invert and remove initial value
        # else:
        # print("ob", labels, values)
        return labels, values
