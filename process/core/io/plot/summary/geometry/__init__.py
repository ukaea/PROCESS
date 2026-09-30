"""Public API for this summary plotting concern."""

from __future__ import annotations

import process.core.io.plot.summary.geometry.build as _build
import process.core.io.plot.summary.geometry.misc as _misc
import process.core.io.plot.summary.geometry.poloidal as _poloidal
import process.core.io.plot.summary.geometry.toroidal as _toroidal

_MODULES = (_build, _misc, _poloidal, _toroidal)
_REGISTRY = {}
for _module in _MODULES:
    _REGISTRY.update({name: getattr(_module, name) for name in _module.__all__})
for _module in _MODULES:
    _module.__dict__.update(_REGISTRY)
arc = _REGISTRY["arc"]
arc_fill = _REGISTRY["arc_fill"]
cumulative_radial_build = _REGISTRY["cumulative_radial_build"]
cumulative_radial_build2 = _REGISTRY["cumulative_radial_build2"]
plot_blanket = _REGISTRY["plot_blanket"]
plot_blkt_pipe_bends = _REGISTRY["plot_blkt_pipe_bends"]
plot_blkt_structure = _REGISTRY["plot_blkt_structure"]
plot_cryostat = _REGISTRY["plot_cryostat"]
plot_first_wall_poloidal_cross_section = _REGISTRY[
    "plot_first_wall_poloidal_cross_section"
]
plot_first_wall_top_down_cross_section = _REGISTRY[
    "plot_first_wall_top_down_cross_section"
]
plot_firstwall = _REGISTRY["plot_firstwall"]
plot_full_machine_poloidal_cross_section = _REGISTRY[
    "plot_full_machine_poloidal_cross_section"
]
plot_geometry_info = _REGISTRY["plot_geometry_info"]
plot_radial_build = _REGISTRY["plot_radial_build"]
plot_shield = _REGISTRY["plot_shield"]
plot_vacuum_vessel_and_divertor = _REGISTRY["plot_vacuum_vessel_and_divertor"]
poloidal_cross_section = _REGISTRY["poloidal_cross_section"]
toroidal_cross_section = _REGISTRY["toroidal_cross_section"]
__all__ = [
    "arc",
    "arc_fill",
    "cumulative_radial_build",
    "cumulative_radial_build2",
    "plot_blanket",
    "plot_blkt_pipe_bends",
    "plot_blkt_structure",
    "plot_cryostat",
    "plot_first_wall_poloidal_cross_section",
    "plot_first_wall_top_down_cross_section",
    "plot_firstwall",
    "plot_full_machine_poloidal_cross_section",
    "plot_geometry_info",
    "plot_radial_build",
    "plot_shield",
    "plot_vacuum_vessel_and_divertor",
    "poloidal_cross_section",
    "toroidal_cross_section",
]
