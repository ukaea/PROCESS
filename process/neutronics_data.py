import numpy as np

class ZeroContinuousFunc:
    """A dummy class that returns 0 whenever called."""
    def __call__(self, x):
        """Return 0 for al cases."""
        if np.isscalar(x):
            return 0.0
        return np.zeros_like(x)

def get_avg_atomic_mass(composition: dict[str, float]) -> float:
    """Calculate the average atomic mass number.
    Parameters
    ----------
    composition:
        a dictionary showing the fraction that each species makes up.
    """
    total_fraction = sum(composition.values())
    return sum(extract_atomic_mass(species) * fraction / total_fraction for species, fraction in composition.items())

def extract_atomic_mass(str) -> float:
    return

material_density_data_bank = ...
material_composition_data_bank = ...
xs_data_bank = ...
breeding_xs_data_bank = ...