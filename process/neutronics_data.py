from dataclasses import dataclass
from itertools import islice
from pathlib import Path

import numpy as np
from numpy import typing as npt
from scipy.constants import Avogadro

from process.exceptions import ProcessValidationError

BARNS_CM2 = 1e-24
N_A = Avogadro
N2N_Q_VALUE = ...
_ATOMIC_MASS = {}


def extract_atomic_mass(isotope: str) -> float:
    """Copied from openmc 0.15.2: openmc.data.data"""
    if not _ATOMIC_MASS:
        mass_file = Path(__file__) / "data" / "mass_1.mas20.txt"
        with open(mass_file) as ame:
            # Read lines in file starting at line 37
            for line in islice(ame, 36, None):
                name = f"{line[20:22].strip()}{int(line[16:19])}"
                mass = float(line[106:109]) + 1e-6 * float(
                    line[110:116] + "." + line[117:123]
                )
                _ATOMIC_MASS[name.lower()] = mass
    return _ATOMIC_MASS[isotope.lower()]


def get_avg_atomic_mass(composition: dict[str, float]) -> float:
    """Calculate the average atomic mass number.
    Parameters
    ----------
    composition:
        a dictionary showing the fraction that each species makes up.
    """
    total_fraction = sum(composition.values())
    return sum(
        extract_atomic_mass(species) * (fraction / total_fraction)
        for species, fraction in composition.items()
    )


material_density_data_bank = ...
material_composition_data_bank = ...
xs_data_bank = ...
breeding_xs_data_bank = ...


def calculate_average_macro_xs(
    composition: dict[str, float],
    micro_xs: dict[str, float | npt.NDArray[np.float64]],
    density: float,
) -> float | npt.NDArray[np.float64]:
    r"""
    Calculate the macroscopic cross-section for a specific energy group and a specific
    reaction (a scalar value), when given the microscopic cross-section values of its
    components and its density.

    Parameters
    ----------
    composition:
        Fraction of each species of atoms of the medium. Does not have to be normalised.
        Given in the format {'species': float(fraction)}.
    micro_xs:
        A dictionary of each species, and their microscopic cross-sections, given in
        [barns].
        Given in the format {'species': float(microscopic cross-section)}.
    density:
        Density of the medium, given in [g/cm^3]

    Notes
    -----
    .. math::
        \Sigma &= N_d\sigma
        &= \frac{N_A}{A} \rho \sigma
    """

    total_fraction = sum(composition.values())
    weighted_micro_xs, weighted_atomic_mass = [], []
    if composition.keys() != micro_xs.keys():
        raise KeyError(
            "The two dictionaries 'composition' and 'micro_xs' must have matching keys."
        )
    for species, fraction in composition.items():
        normalized_fraction = fraction / total_fraction
        weighted_atomic_mass.append(
            normalized_fraction * extract_atomic_mass(species)
        )
        weighted_micro_xs.append(normalized_fraction * micro_xs[species])
    avg_sigma = np.sum(weighted_micro_xs, axis=0)
    avg_mass_amu = sum(weighted_atomic_mass)
    return (
        (BARNS_CM2 * N_A) / avg_mass_amu * avg_sigma * density
    )  # N_A/A * rho * sigma


def discretize_xs(
    continuous_xs, group_structure: list[float]
) -> npt.NDArray[np.float64]:
    """
    Discretise a continuous cross-section function into a group structure of n discrete
    flaots.

    Parameters
    ----------
    continuous_xs:
        continuous cross-section function to be discretized.
    group_structure:
        group structure of neutron energies, the n+1 energy bin boundaries for the n
        neutron groups, in descending energies.

    Returns
    -------
    :
        microscopic cross-section values discretized according to the group structure.
    """
    return group_structure, continuous_xs


def scattering_weight_matrix(
    group_structure: list[float], atomic_mass: float
) -> npt.NDArray:
    """
    Parameters
    ----------
    group_structure:
        the n+1 energy bin boundaries for the n neutron groups, in descending energies.
    atomic_mass:
        atomic mass of the medium

    Returns
    -------
    :
        A lower triangular matrix, where the i-th row contains the normalized weights for
        scattering down from the i-th bin to the j-th bin. The main-diagonal contains
        the self-scattering cross-section.
        e.g. [2,1] would be the fraction of elastic scattering reactions that causes
        neutrons from group 3 to scatter into group 2.
        The upper triangle must be all zeros.
        np.sum(axis=1) == np.ones(len(group_structure)-1).
    """
    return atomic_mass, group_structure


def expand_macro_neutron_multiplication_xs_into_matrix(
    discrete_n2n_xs: npt.NDArray[np.float64],
    group_structure: list[float],
    q_value: float,
) -> npt.NDArray:
    """Instead of only having the macroscopic cross-section values for the (n,2n)
    reaction for each group, calculate the macroscopic cross-section of neutron in the
    [i]-th bin producing a neutron in the [j]-bin, recorded as the [i,j] element in the
    matrix.

    Parameters
    ----------
    discrete_n2n_xs:
        The group-wise macroscopic cross-section for the n2n reaction.
        A 1D numpy array, with len==number of neutron groups.
    group_structure:
        The n+1 energy bin boundaries for the n neutron groups, in descending energies.
    q_value:
        The difference in q-value.

    Returns
    -------
    :
        A macroscopic neutron multiplication cross-section matrix, such that each row
        should sum to = 2 * discrete_n2n_xs (since two neutrons should be produced per
        neutron consumed in the (n,2n) reaction).

    Notes
    -----
    This is a three body problem. TODO: further investigation is needed to figure out
    how to distribute the two outputted neutron's energies!
    """


class ZeroContinuousFunc:
    """A dummy class that returns 0 whenever called."""

    def __call__(self, x):
        """Return 0 for al cases."""
        if np.isscalar(x):
            return 0.0
        return np.zeros_like(x)


@dataclass
class MaterialMacroInfo:
    """Material information"""

    sigma_t: npt.NDArray[np.float64]  # Sigma_total, 1D array of len = n
    sigma_s: npt.NDArray  # Sigma_scattering from group i to j, 2D array of n*n
    group_structure: npt.NDArray  # energy bin edges, 1D array of len = n+1
    avg_atomic_mass: float  # average atomic mass (weighted by fraction)

    def __post_init__(self):
        """Validation to confirm the shape is correct."""
        # force into float or numpy arrays.
        self.sigma_t = np.array(self.sigma_t, dtype=float)
        self.sigma_s = np.array(self.sigma_s, dtype=float)
        self.group_structure = np.array(self.group_structure, dtype=float)
        self.avg_atomic_mass = float(self.avg_atomic_mass)

        if np.diff(self.group_structure)>=0:
            raise ValueError(
                "The group structure must be defined beginning from the highest energy "
                "bin (i.e. lowest lethargy bin) edge, descending to the lowest energy. "
                "Similarly the cross-section must be arranged with the highest energy "
                "group first, and the lowest energy group last."
            )

        if np.shape(self.sigma_t) != (self.n_groups,):
            raise ProcessValidationError(
                f"total group-wise cross-sections should have {self.n_groups} "
                "groups as specified by the group_structure."
            )
        if np.shape(self.sigma_s) != (self.n_groups, self.n_groups):
            raise ProcessValidationError(
                "Group-wise scattering cross-sections be a square matrix of "
                f"shape n*n, where n= number of groups = {self.n_groups}."
            )
        if (self.sigma_s.sum(axis=1) > self.sigma_t).any():
            raise ProcessValidationError(
                "Each group's scattering cross-section should be smaller than "
                "or equal to its total cross-section."
            )

    @property
    def n_groups(self):
        """Store this attribute upon first retrieval."""
        if not hasattr(self, "_n_groups"):
            self._n_groups = len(self.group_structure) - 1
        return self._n_groups


def get_material_nuclear_data(
    material: str, group_structure: list[float]
) -> MaterialMacroInfo:
    """
    The constants that is directly used.

    Parameters
    ----------
    materials:
        Which material that we want to get the data out of.
    group_structure:
        the n+1 energy bin boundaries for the n neutron groups, in descending energies.

    Returns
    -------
    discrete_macro_total_xs:
        All group-wise total cross-sections, given in ascending order.
    macroscopic scattering cross-section matrix:
        A lower-triangular matrix of total scattering cross-sections,
        e.g. [2,1] would be the scattering cross-section from neutron group 3 to group 2.
        The upper triangle MUST be zeros.
    avg_atomic_mass:
        average atomic mass, used for further calcuations

    Notes
    -----
    When calculating the group-wise diffusion coefficients, the i-th group's scattering
    cross-section :math:`\\Sigma_{s}` is the i-th element on the main diagonal of the
    macroscopic scattering cross-section matrix.
    """
    density = material_density_data_bank[material]
    composition = material_composition_data_bank[material]
    avg_atomic_mass = get_avg_atomic_mass(composition)

    # dicts of {"isotope": npt.NDArray[np.float64] 1D arrays}
    micro_total_xs = {}
    micro_scattering_xs = {}
    micro_n2n_xs = {}
    for species in composition:
        total_xs_continuous, elastic_scattering_xs_continuous = xs_data_bank[
            species
        ]
        n2n_xs_continuous = breeding_xs_data_bank.get(
            species, ZeroContinuousFunc
        )
        micro_total_xs[species] = discretize_xs(
            total_xs_continuous, group_structure
        )
        micro_scattering_xs[species] = discretize_xs(
            elastic_scattering_xs_continuous, group_structure
        )
        micro_n2n_xs[species] = discretize_xs(
            n2n_xs_continuous, group_structure
        )

    discrete_macro_total_xs = calculate_average_macro_xs(
        composition, micro_total_xs, density
    )
    discrete_macro_scattering_xs = calculate_average_macro_xs(
        composition, micro_scattering_xs, density
    )
    source_matrix = (
        scattering_weight_matrix(group_structure, avg_atomic_mass).T
        * discrete_macro_scattering_xs
    ).T
    n2n_susceptible_species = ...
    for (
        modified_composition_file,
        modified_density,
        q_value,
    ) in n2n_susceptible_species:
        source_matrix += expand_macro_neutron_multiplication_xs_into_matrix(
            calculate_average_macro_xs(composition, micro_n2n_xs, density),
            group_structure,
            q_value,
        )

    return MaterialMacroInfo(
        discrete_macro_total_xs,
        source_matrix,
        group_structure,
        avg_atomic_mass,
    )
