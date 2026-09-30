"""PROCESS summary plotting package."""

from __future__ import annotations

import process.core.io.plot.summary.api as _api
import process.core.io.plot.summary.common as _common
import process.core.io.plot.summary.geometry as _geometry
import process.core.io.plot.summary.magnets as _magnets
import process.core.io.plot.summary.plasma as _plasma
import process.core.io.plot.summary.power_flow as _power_flow
import process.core.io.plot.summary.profiles as _profiles
import process.core.io.plot.summary.reporting as _reporting
import process.core.io.plot.summary.time_profiles as _time_profiles

_MODULES = (
    _api,
    _common,
    _geometry,
    _magnets,
    _plasma,
    _power_flow,
    _profiles,
    _reporting,
    _time_profiles,
)
_REGISTRY = {}
for _module in _MODULES:
    _REGISTRY.update({name: getattr(_module, name) for name in _module.__all__})
for _module in _MODULES:
    _module.__dict__.update(_REGISTRY)


plot_summary = _api.plot_summary
main_plot = _api.main_plot
create_thickness_builds = _api.create_thickness_builds

__all__ = ["create_thickness_builds", "main_plot", "plot_summary"]
