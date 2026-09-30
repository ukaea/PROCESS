"""Public API for this summary plotting concern."""

from __future__ import annotations

import process.core.io.plot.summary.plasma.confinement as _confinement
import process.core.io.plot.summary.plasma.current_drive as _current_drive
import process.core.io.plot.summary.plasma.overview as _overview
import process.core.io.plot.summary.plasma.physics as _physics

_MODULES = (_confinement, _current_drive, _overview, _physics)
_REGISTRY = {}
for _module in _MODULES:
    _REGISTRY.update({name: getattr(_module, name) for name in _module.__all__})
for _module in _MODULES:
    _module.__dict__.update(_REGISTRY)
plot_bootstrap_comparison = _REGISTRY["plot_bootstrap_comparison"]
plot_brunner_divertor_power_split_comparison_stackplot = _REGISTRY[
    "plot_brunner_divertor_power_split_comparison_stackplot"
]
plot_confinement_time_comparison = _REGISTRY["plot_confinement_time_comparison"]
plot_current_drive_info = _REGISTRY["plot_current_drive_info"]
plot_detailed_plasma_parameters = _REGISTRY["plot_detailed_plasma_parameters"]
plot_magnetic_fields_in_plasma = _REGISTRY["plot_magnetic_fields_in_plasma"]
plot_main_plasma_information = _REGISTRY["plot_main_plasma_information"]
plot_max_normalised_beta_comparison = _REGISTRY["plot_max_normalised_beta_comparison"]
plot_plasma = _REGISTRY["plot_plasma"]
plot_plasma_coloumb_logarithms = _REGISTRY["plot_plasma_coloumb_logarithms"]
plot_plasma_current_comparison = _REGISTRY["plot_plasma_current_comparison"]
plot_plasma_outboard_toroidal_ripple_map = _REGISTRY[
    "plot_plasma_outboard_toroidal_ripple_map"
]
plot_sol_power_decay_length_comparison = _REGISTRY[
    "plot_sol_power_decay_length_comparison"
]
reaction_plot_grid = _REGISTRY["reaction_plot_grid"]
__all__ = [
    "plot_bootstrap_comparison",
    "plot_brunner_divertor_power_split_comparison_stackplot",
    "plot_confinement_time_comparison",
    "plot_current_drive_info",
    "plot_detailed_plasma_parameters",
    "plot_magnetic_fields_in_plasma",
    "plot_main_plasma_information",
    "plot_max_normalised_beta_comparison",
    "plot_plasma",
    "plot_plasma_coloumb_logarithms",
    "plot_plasma_current_comparison",
    "plot_plasma_outboard_toroidal_ripple_map",
    "plot_sol_power_decay_length_comparison",
    "reaction_plot_grid",
]
