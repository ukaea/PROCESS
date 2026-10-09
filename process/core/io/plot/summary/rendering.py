"""Reusable rendering primitives for summary plots."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from collections.abc import Iterable

    from matplotlib.axes import Axes
    from matplotlib.text import Annotation, Text
    from matplotlib.transforms import Transform


@dataclass(frozen=True)
class TextPanel:
    """Declarative text panel rendered through a Matplotlib axes."""

    x: float
    y: float
    text: str
    options: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class ArrowSpec:
    """Declarative annotation arrow."""

    start: tuple[float, float]
    end: tuple[float, float]
    transform: Transform
    options: dict[str, Any] = field(default_factory=dict)


def draw_text_panels(axis: Axes, panels: Iterable[TextPanel]) -> list[Text]:
    """Render a sequence of declarative text panels."""
    return [axis.text(panel.x, panel.y, panel.text, **panel.options) for panel in panels]


def draw_arrows(axis: Axes, arrows: Iterable[ArrowSpec]) -> list[Annotation]:
    """Render a sequence of declarative arrows."""
    return [
        axis.annotate(
            "",
            xy=arrow.end,
            xytext=arrow.start,
            xycoords=arrow.transform,
            arrowprops=arrow.options,
        )
        for arrow in arrows
    ]


__all__ = [
    "ArrowSpec",
    "TextPanel",
    "draw_arrows",
    "draw_text_panels",
]
