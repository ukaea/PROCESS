"""Public API for this summary plotting concern."""

from __future__ import annotations

import process.core.io.plot.summary.reporting.constraints as _constraints
import process.core.io.plot.summary.reporting.layouts as _layouts
import process.core.io.plot.summary.reporting.misc as _misc
import process.core.io.plot.summary.reporting.panels as _panels

_MODULES = (_constraints, _layouts, _misc, _panels)
_REGISTRY = {}
for _module in _MODULES:
    _REGISTRY.update({name: getattr(_module, name) for name in _module.__all__})
for _module in _MODULES:
    _module.__dict__.update(_REGISTRY)
RadialBuild = _REGISTRY["RadialBuild"]
draw_bend = _REGISTRY["draw_bend"]
plot_centre_cross = _REGISTRY["plot_centre_cross"]
plot_cover_page = _REGISTRY["plot_cover_page"]
plot_density_limit_comparison = _REGISTRY["plot_density_limit_comparison"]
plot_ebw_ecrh_coupling_graph = _REGISTRY["plot_ebw_ecrh_coupling_graph"]
plot_equality_constraint_equations = _REGISTRY["plot_equality_constraint_equations"]
plot_fw_90_deg_pipe_bend = _REGISTRY["plot_fw_90_deg_pipe_bend"]
plot_h_threshold_comparison = _REGISTRY["plot_h_threshold_comparison"]
plot_header = _REGISTRY["plot_header"]
plot_inequality_constraint_equations = _REGISTRY["plot_inequality_constraint_equations"]
plot_info = _REGISTRY["plot_info"]
plot_iteration_variables = _REGISTRY["plot_iteration_variables"]
plot_lower_vertical_build = _REGISTRY["plot_lower_vertical_build"]
plot_separatrix_power_split = _REGISTRY["plot_separatrix_power_split"]
plot_upper_vertical_build = _REGISTRY["plot_upper_vertical_build"]
__all__ = [
    "RadialBuild",
    "draw_bend",
    "plot_centre_cross",
    "plot_cover_page",
    "plot_density_limit_comparison",
    "plot_ebw_ecrh_coupling_graph",
    "plot_equality_constraint_equations",
    "plot_fw_90_deg_pipe_bend",
    "plot_h_threshold_comparison",
    "plot_header",
    "plot_inequality_constraint_equations",
    "plot_info",
    "plot_iteration_variables",
    "plot_lower_vertical_build",
    "plot_separatrix_power_split",
    "plot_upper_vertical_build",
]
