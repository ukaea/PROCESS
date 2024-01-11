"""
Module to hold geometry parameterisation dataclasses common to multiple reactor components
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class RectangleGeometry:
    """Holds data for rectangular geometries"""

    anchor_x: float
    """rectangle x coordinate anchor point"""
    anchor_z: float
    """rectangle z coordinate anchor point"""
    width: float
    """rectangle width"""
    height: float
    """rectangle height"""


@dataclass
class ArbitraryGeometry:
    """Holds radial and vertical coordinates for arbitrary reactor component shapes

    Example: a triangular shaped component with vertices [(0,0), (1,0), (0,1)] would be represented by ArbitraryGeometry(rs=[0,1,0], zs=[0,0,1])"""

    rs: np.ndarray
    """outboard and inboard radial coordinates"""
    zs: np.ndarray
    """outboard and inboard vertical coordinates"""
