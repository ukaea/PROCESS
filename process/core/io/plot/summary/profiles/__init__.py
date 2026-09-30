"""Public API for this summary plotting concern."""

from __future__ import annotations

import process.core.io.plot.summary.profiles.atomic as _atomic
import process.core.io.plot.summary.profiles.misc as _misc
import process.core.io.plot.summary.profiles.plasma as _plasma
import process.core.io.plot.summary.profiles.radiation as _radiation
import process.core.io.plot.summary.profiles.stress as _stress

_MODULES = (_atomic, _misc, _plasma, _radiation, _stress)
_REGISTRY = {}
for _module in _MODULES:
    _REGISTRY.update({name: getattr(_module, name) for name in _module.__all__})
for _module in _MODULES:
    _module.__dict__.update(_REGISTRY)
interp1d_profile = _REGISTRY["interp1d_profile"]
plot_beta_profiles = _REGISTRY["plot_beta_profiles"]
plot_collision_frequency_profile = _REGISTRY["plot_collision_frequency_profile"]
plot_collision_time_profile = _REGISTRY["plot_collision_time_profile"]
plot_cs_hoop_stress_contour_profile = _REGISTRY["plot_cs_hoop_stress_contour_profile"]
plot_cs_hoop_stress_profile = _REGISTRY["plot_cs_hoop_stress_profile"]
plot_cs_radial_stress_contour_profile = _REGISTRY[
    "plot_cs_radial_stress_contour_profile"
]
plot_cs_radial_stress_profile = _REGISTRY["plot_cs_radial_stress_profile"]
plot_cs_stress_time_profile = _REGISTRY["plot_cs_stress_time_profile"]
plot_cs_tresca_2d_contour = _REGISTRY["plot_cs_tresca_2d_contour"]
plot_cs_vertical_stress_profile = _REGISTRY["plot_cs_vertical_stress_profile"]
plot_cs_von_mises_2d_contour = _REGISTRY["plot_cs_von_mises_2d_contour"]
plot_cumulative_plasma_thermal_energy_profiles = _REGISTRY[
    "plot_cumulative_plasma_thermal_energy_profiles"
]
plot_debye_length_profile = _REGISTRY["plot_debye_length_profile"]
plot_electron_frequency_profile = _REGISTRY["plot_electron_frequency_profile"]
plot_fusion_rate_contours = _REGISTRY["plot_fusion_rate_contours"]
plot_fusion_rate_profiles = _REGISTRY["plot_fusion_rate_profiles"]
plot_ion_charge_profile = _REGISTRY["plot_ion_charge_profile"]
plot_ion_frequency_profile = _REGISTRY["plot_ion_frequency_profile"]
plot_ion_slowing_down_time_profile = _REGISTRY["plot_ion_slowing_down_time_profile"]
plot_jprofile = _REGISTRY["plot_jprofile"]
plot_larmor_radius_profile = _REGISTRY["plot_larmor_radius_profile"]
plot_line_brem_loss_function_profile = _REGISTRY["plot_line_brem_loss_function_profile"]
plot_line_brem_power_density_profile = _REGISTRY["plot_line_brem_power_density_profile"]
plot_mean_free_path_profile = _REGISTRY["plot_mean_free_path_profile"]
plot_n_profiles = _REGISTRY["plot_n_profiles"]
plot_plasma_effective_charge_profile = _REGISTRY["plot_plasma_effective_charge_profile"]
plot_plasma_poloidal_pressure_contours = _REGISTRY[
    "plot_plasma_poloidal_pressure_contours"
]
plot_plasma_pressure_gradient_profiles = _REGISTRY[
    "plot_plasma_pressure_gradient_profiles"
]
plot_plasma_pressure_profiles = _REGISTRY["plot_plasma_pressure_profiles"]
plot_plasma_thermal_energy_profiles = _REGISTRY["plot_plasma_thermal_energy_profiles"]
plot_qprofile = _REGISTRY["plot_qprofile"]
plot_rad_contour = _REGISTRY["plot_rad_contour"]
plot_resistivity_profile = _REGISTRY["plot_resistivity_profile"]
plot_t_profiles = _REGISTRY["plot_t_profiles"]
plot_velocity_profile = _REGISTRY["plot_velocity_profile"]
plot_vertical_stress_contour_profile = _REGISTRY["plot_vertical_stress_contour_profile"]
profiles_with_pedestal = _REGISTRY["profiles_with_pedestal"]
read_imprad_data = _REGISTRY["read_imprad_data"]
__all__ = [
    "interp1d_profile",
    "plot_beta_profiles",
    "plot_collision_frequency_profile",
    "plot_collision_time_profile",
    "plot_cs_hoop_stress_contour_profile",
    "plot_cs_hoop_stress_profile",
    "plot_cs_radial_stress_contour_profile",
    "plot_cs_radial_stress_profile",
    "plot_cs_stress_time_profile",
    "plot_cs_tresca_2d_contour",
    "plot_cs_vertical_stress_profile",
    "plot_cs_von_mises_2d_contour",
    "plot_cumulative_plasma_thermal_energy_profiles",
    "plot_debye_length_profile",
    "plot_electron_frequency_profile",
    "plot_fusion_rate_contours",
    "plot_fusion_rate_profiles",
    "plot_ion_charge_profile",
    "plot_ion_frequency_profile",
    "plot_ion_slowing_down_time_profile",
    "plot_jprofile",
    "plot_larmor_radius_profile",
    "plot_line_brem_loss_function_profile",
    "plot_line_brem_power_density_profile",
    "plot_mean_free_path_profile",
    "plot_n_profiles",
    "plot_plasma_effective_charge_profile",
    "plot_plasma_poloidal_pressure_contours",
    "plot_plasma_pressure_gradient_profiles",
    "plot_plasma_pressure_profiles",
    "plot_plasma_thermal_energy_profiles",
    "plot_qprofile",
    "plot_rad_contour",
    "plot_resistivity_profile",
    "plot_t_profiles",
    "plot_velocity_profile",
    "plot_vertical_stress_contour_profile",
    "profiles_with_pedestal",
    "read_imprad_data",
]
