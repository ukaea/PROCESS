import json
from pathlib import Path

from process.exceptions import ProcessValueError
from process.fortran import stellarator_configuration

HELIAS5B = {
    "name": "Helias 5b",
    "rmajor_ref": 22.2,
    "rminor_ref": 1.80,
    "aspect_ref": 12.33,
    "coil_rmajor": 22.44,
    "coil_rminor": 4.76,
    "bt_ref": 5.6,
    "WP_area": 0.8 * 0.6,
    "WP_bmax": 11.44,
    "symmetry": 5,
    "coilspermodule": 10,
    "a1": 0.688,
    "a2": 0.025,
    "vol_plasma": 1422.63,  # This value is for Helias 5
    "dmin": 0.84,
    "max_portsize_width": 2.12,
    "plasma_surface": 1960.0,  # Plasma Surface
    "maximal_coil_height": 12.7,  # [m] Full height max point to min point
    "coilsurface": 4817.7,  # Coil surface, dimensionfull. At reference point
    "coillength": 1680.0,  # Central filament length of machine with outer radius 1m.
    "I0": 13.06,  # Coil Current needed to produce 1T on axis in [MA] at outer radius 1m
    "inductance": 1655.76e-6,  # inductance in muH
    "WP_ratio": 1.2,  # The fit values in stellarator config class should be calculated using this value.
    "max_force_density": 120.0,  # [MN/m^3]
    "max_force_density_MNm": 98.0,  # [MN/m]
    "max_lateral_force_density": 92.4,  # [MN/m^3]
    "max_radial_force_density": 113.5,  # [MN/m^3]
    "centering_force_max_MN": 189.5,
    "centering_force_min_MN": -55.7,
    "centering_force_avg_MN": 93.0,
    "min_plasma_coil_distance": 1.9,
    "derivative_min_LCFS_coils_dist": -1.0,  # this is approximated for now
    "min_bend_radius": 1.0,  # [m]
    "neutron_peakfactor": 1.6,
    "epseff": 0.015,
}


HELIAS4 = {
    "name": "Helias 4",
    # Reference point where all the other variables are determined from
    # Plasma outer radius
    "rmajor_ref": 17.6,
    "rminor_ref": 2.0,
    "aspect_ref": 8.8,
    # Coil radii
    "coil_rmajor": 18.39,
    "coil_rminor": 4.94,
    "bt_ref": 5.6,
    "WP_area": 0.8 * 0.6,
    "WP_bmax": 11.51,
    "symmetry": 4,
    "coilspermodule": 10,
    "a1": 0.676,
    "a2": 0.029,
    "vol_plasma": 1380.0,
    "dmin": 1.08,
    "max_portsize_width": 3.24,
    "plasma_surface": 1900.0,
    "maximal_coil_height": 13.34,  # [m] Full height max point to min point
    "coilsurface": 4100.0,  # Coil surface, dimensionfull. At reference point
    "coillength": 1435.07,  # Central filament length of machine with outer radius 1m.
    "I0": 13.146,  # Coil Current needed to produce b0 on axis in [MA] at reference point
    "inductance": 1290.4e-6,  # inductance/R*A^2 in muH
    "WP_ratio": 1.3,
    "max_force_density": 120.0,  # [MN/m^3]
    "max_force_density_MNm": 98.0,  # [MN/m]
    "max_lateral_force_density": 87.9,  # [MN/m^3]
    "max_radial_force_density": 109.9,  # [MN/m^3]
    "centering_force_max_MN": 226.0,
    "centering_force_min_MN": -35.3,
    "centering_force_avg_MN": 125.8,
    "min_plasma_coil_distance": 1.7,
    "derivative_min_LCFS_coils_dist": -1.0,  # this is approximated for now
    "min_bend_radius": 0.86,  # [m]
    "neutron_peakfactor": 1.6,
    "epseff": 0.015,
}

HELIAS3 = {
    "name": "Helias 3",
    # Reference point where all the other variables are determined from
    # Plasma outer radius
    "rmajor_ref": 13.86,
    "rminor_ref": 2.18,
    "aspect_ref": 6.36,
    # Coil radii
    "coil_rmajor": 14.53,
    "coil_rminor": 6.12,
    "bt_ref": 5.6,
    "WP_bmax": 12.346,
    "WP_area": 0.8 * 0.6,
    "symmetry": 3,
    "coilspermodule": 10,
    # Bmax fit parameters
    "a1": 0.56,
    "a2": 0.030,
    "vol_plasma": 1300.8,
    "dmin": 1.145,
    "max_portsize_width": 3.24,  # ??? guess. not ready yet
    "plasma_surface": 1600.00,
    "maximal_coil_height": 17.74,  # [m] Full height max point to min point
    "coilsurface": 4240.0,  # Coil surface, dimensionfull. At reference point
    "coillength": 1287.3,  # Central filament length of machine with outer radius 1m.
    "I0": 14.23,  # Coil Current needed to produce 1T on axis in [MA] at outer radius 1m
    "inductance": 1250.7e-6,  # inductance in muH
    "WP_ratio": 1.3,
    "max_force_density": 120.0,  # [MN/m]
    "max_force_density_MNm": 98.0,
    "max_lateral_force_density": 96.6,  # [MN/m^3]
    "max_radial_force_density": 130.5,  # [MN/m^3]
    "centering_force_max_MN": 428.1,
    "centering_force_min_MN": -70.3,
    "centering_force_avg_MN": 240.9,
    "min_plasma_coil_distance": 1.78,
    "derivative_min_LCFS_coils_dist": -1.0,  # this is approximated for now
    "min_bend_radius": 1.145,  # [m]
    "neutron_peakfactor": 1.6,
    "epseff": 0.015,
}

W7X30 = {
    "name": "W7X-30",
    # Reference point where all the other variables are determined from
    # Plasma outer radius
    "rmajor_ref": 5.50,
    "rminor_ref": 0.49,
    "aspect_ref": 11.2,
    # Coil radii
    "coil_rmajor": 5.62,
    "coil_rminor": 1.36,
    "bt_ref": 3.0,
    "WP_area": 0.18 * 0.15,
    "WP_bmax": 10.6,
    "symmetry": 5,
    "coilspermodule": 6,
    "a1": 0.98,
    "a2": 0.041,
    "vol_plasma": 26.4,
    "dmin": 0.21,
    "max_portsize_width": 0.5,
    "plasma_surface": 128.3,
    "maximal_coil_height": 3.6,  # [m] Full height max point to min point
    "coilsurface": 370.0,  # Coil surface, dimensionfull. At reference point
    "coillength": 303.4,  # Central filament length of machine with outer radius 1m.
    "I0": 2.9,  # Coil Current needed to produce b0 on axis in [MA] at reference point
    "inductance": 252.7e-6,  # inductance/R*A^2 in muH
    "WP_ratio": 1.2,
    "max_force_density": 350.0,  # [MN/m^3]
    "max_force_density_MNm": 98.0,  # [MN/m]
    "max_lateral_force_density": 271.1,  # [MN/m^3]
    "max_radial_force_density": 305.2,  # [MN/m^3]
    "centering_force_max_MN": 7.95,
    "centering_force_min_MN": -2.15,
    "centering_force_avg_MN": 3.46,
    "min_plasma_coil_distance": 0.45,
    "derivative_min_LCFS_coils_dist": -1.0,  # this is approximated for now
    "min_bend_radius": 0.186,  # [m]
    "neutron_peakfactor": 1.6,
    "epseff": 0.015,
}

W7X50 = {
    "name": "W7X-50",
    # Reference point where all the other variables are determined from
    # Plasma outer radius
    "rmajor_ref": 5.5,
    "rminor_ref": 0.49,
    "aspect_ref": 11.2,
    # Coil radii
    "coil_rmajor": 5.62,
    "coil_rminor": 1.18,
    "bt_ref": 3.0,
    "WP_area": 0.18 * 0.15,
    "WP_bmax": 6.3,
    "symmetry": 5,
    "coilspermodule": 10,
    "a1": 0.66,
    "a2": 0.025,
    "vol_plasma": 26.4,
    "dmin": 0.28,
    "max_portsize_width": 0.3,
    "plasma_surface": 128.3,
    "maximal_coil_height": 3.1,  # [m] Full height max point to min point
    "coilsurface": 299.85,  # Coil surface, dimensionfull. At reference point
    "coillength": 420.67,  # Central filament length of machine with outer radius 1m.
    "I0": 1.745,  # Coil Current needed to produce b0 on axis in [MA] at reference point
    "inductance": 412.4e-6,  # inductance/R*A^2 in muH
    "WP_ratio": 1.2,
    "max_force_density": 250.0,  # [MN/m^3]
    "max_force_density_MNm": 98.0,  # [MN/m]
    "max_lateral_force_density": 116.4,  # [MN/m^3]
    "max_radial_force_density": 148.0,  # [MN/m^3]
    "centering_force_max_MN": 2.99,
    "centering_force_min_MN": -1.29,
    "centering_force_avg_MN": 1.61,
    "min_plasma_coil_distance": 0.39,
    "derivative_min_LCFS_coils_dist": -1.0,  # this is approximated for now
    "min_bend_radius": 0.39,  # [m]
    "neutron_peakfactor": 1.6,
    "epseff": 0.015,
}


def load_stellarator_config(istell: int, config_file: Path | None):
    """Load the appropriate Stellarator machine configuration
    given the `istell` switch:

    istell = 1: Helias5 machine
    istell = 2: Helias4 machine
    istell = 3: Helias3 machine
    istell = 4: w7x30 machine
    istell = 5: w7x50 machine
    istell = 6: Init from json
    """
    match istell:
        case 1:
            machine_config = HELIAS5B
        case 2:
            machine_config = HELIAS4
        case 3:
            machine_config = HELIAS3
        case 4:
            machine_config = W7X30
        case 5:
            machine_config = W7X50
        case 6:
            if config_file is None:
                raise ProcessValueError("Stellarator config file is None but istell=6")

            with open(config_file) as f:
                machine_config = json.load(f)
        case _:
            raise ProcessValueError(f"{istell=} is not an integer in the range [1, 6]")

    for variable_name, variable_value in machine_config.items():
        setattr(
            stellarator_configuration,
            f"stella_config_{variable_name.lower()}",
            variable_value,
        )
