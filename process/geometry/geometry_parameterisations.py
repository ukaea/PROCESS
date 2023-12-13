"""
Module to hold geometry parameterisation dataclasses common to multiple reactor components
"""
from dataclasses import dataclass


@dataclass
class RectangleGeometry:
    """Holds data for rectangular geometries"""

    center_x: float
    """rectangle x coordinate anchor point"""
    center_z: float
    """rectangle z coordinate anchor point"""
    width: float
    """rectangle width"""
    height: float
    """rectangle height"""
