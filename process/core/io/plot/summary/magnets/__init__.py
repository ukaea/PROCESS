"""Public API for this summary plotting concern."""

from __future__ import annotations

import process.core.io.plot.summary.magnets.cables as _cables
import process.core.io.plot.summary.magnets.cs as _cs
import process.core.io.plot.summary.magnets.pf as _pf
import process.core.io.plot.summary.magnets.tf as _tf

_MODULES = (_cables, _cs, _pf, _tf)
_REGISTRY = {}
for _module in _MODULES:
    _REGISTRY.update({name: getattr(_module, name) for name in _module.__all__})
for _module in _MODULES:
    _module.__dict__.update(_REGISTRY)
TF_outboard = _REGISTRY["TF_outboard"]
plot_cable_in_conduit_cable = _REGISTRY["plot_cable_in_conduit_cable"]
plot_corc_cable_geometry = _REGISTRY["plot_corc_cable_geometry"]
plot_cs_coil_structure = _REGISTRY["plot_cs_coil_structure"]
plot_cs_turn_structure = _REGISTRY["plot_cs_turn_structure"]
plot_hts_tape_geometry = _REGISTRY["plot_hts_tape_geometry"]
plot_magnetics_info = _REGISTRY["plot_magnetics_info"]
plot_pf_coils = _REGISTRY["plot_pf_coils"]
plot_pf_cs_plasma_mutual_inductance = _REGISTRY["plot_pf_cs_plasma_mutual_inductance"]
plot_pf_dimensions = _REGISTRY["plot_pf_dimensions"]
plot_physics_info = _REGISTRY["plot_physics_info"]
plot_quench_time_evolution = _REGISTRY["plot_quench_time_evolution"]
plot_resistive_tf_info = _REGISTRY["plot_resistive_tf_info"]
plot_resistive_tf_wp = _REGISTRY["plot_resistive_tf_wp"]
plot_superconducting_tf_wp = _REGISTRY["plot_superconducting_tf_wp"]
plot_tf_cable_in_conduit_turn = _REGISTRY["plot_tf_cable_in_conduit_turn"]
plot_tf_coil_structure = _REGISTRY["plot_tf_coil_structure"]
plot_tf_coils = _REGISTRY["plot_tf_coils"]
plot_tf_corc_cable_summary_box = _REGISTRY["plot_tf_corc_cable_summary_box"]
plot_tf_croco_turn = _REGISTRY["plot_tf_croco_turn"]
plot_tf_stress = _REGISTRY["plot_tf_stress"]
secs_to_hms = _REGISTRY["secs_to_hms"]
__all__ = [
    "TF_outboard",
    "plot_cable_in_conduit_cable",
    "plot_corc_cable_geometry",
    "plot_cs_coil_structure",
    "plot_cs_turn_structure",
    "plot_hts_tape_geometry",
    "plot_magnetics_info",
    "plot_pf_coils",
    "plot_pf_cs_plasma_mutual_inductance",
    "plot_pf_dimensions",
    "plot_physics_info",
    "plot_quench_time_evolution",
    "plot_resistive_tf_info",
    "plot_resistive_tf_wp",
    "plot_superconducting_tf_wp",
    "plot_tf_cable_in_conduit_turn",
    "plot_tf_coil_structure",
    "plot_tf_coils",
    "plot_tf_corc_cable_summary_box",
    "plot_tf_croco_turn",
    "plot_tf_stress",
    "secs_to_hms",
]
