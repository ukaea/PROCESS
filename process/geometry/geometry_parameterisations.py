from dataclasses import dataclass


@dataclass
class RectangleGeometry:
    center_x: float
    center_z: float
    width: float
    height: float
