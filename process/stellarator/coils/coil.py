"""
This is a conceptual draft of the coil class, to calculated space constraints for each coil separately
"""

from dataclasses import dataclass
from typing import List

@dataclass
class Coil:

    B_max_ref: float
    B_max: float = None

    def calculate_winding_pack():
        pass

    def calculate_B_max():
        pass

    def calculate_number_of_turns():
        pass

    def calculate_dump_voltage(self):
        pass

    def check_coil_plasma_distance(self):
        pass



@dataclass
class Coil_set():
    coils: List[Coil]

    def check_coil_coil_distance(self):
        pass

    def check_plasma_coil_distances(self):
        pass