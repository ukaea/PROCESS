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


material_density_data_bank = {"stainless steel": {"Fe": 6.0}}
material_composition_data_bank = {"stainless steel": {"Fe60": 1.0}}
xs_data_bank = ...
breeding_xs_data_bank = ...
fission_xs_data_bank = ...


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
        Possible to have some missing values here.
    density:
        Density of the medium, given in [g/cm^3]

    Returns
    -------
    :
        macroscopic cross-section of the material.

    Notes
    -----
    .. math::
        \Sigma &= N_d\sigma
        &= \frac{N_A}{A} \rho \sigma
    """

    total_fraction = sum(composition.values())
    weighted_micro_xs, weighted_atomic_mass = [], []
    if extra_species:=set(micro_xs.keys()).difference(set(composition.keys())):
        raise KeyError(
            f"micro_xs contains species not specified by the composition: {extra_species}"
        )
    for species, fraction in composition.items():
        if species in micro_xs:
            frac = fraction / total_fraction

            weighted_atomic_mass.append(frac * extract_atomic_mass(species))
            weighted_micro_xs.append(frac * micro_xs[species])

    avg_sigma = np.sum(weighted_micro_xs, axis=0)
    avg_mass_amu = sum(weighted_atomic_mass)

    # N_A/A * rho * sigma
    return N_A / avg_mass_amu * density * (avg_sigma * BARNS_CM2)


def discretize_xs(
    continuous_xs, group_structure: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """
    Discretise a continuous cross-section function into a group structure of n discrete
    floats. This is done by taking the average xs per bin. (assuming constant-lethargy within bin.)

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
        1D array.
    """
    lethargy = np.log10(group_structure)
    lower_leth, upper_leth = lethargy[:-1], lethargy[1:]
    leth_spacing = np.diff(lethargy)
    continuous_xs_lethargy_input = lambda e: continuous_xs(10**-e)

    return integrate(continuous_xs_lethargy_input, lower_leth, upper_leth)/leth_spacing

def discretize_scattering_xs(
    continuous_xs, group_structure: npt.NDArray[np.float64], atomic_mass: float
) -> npt.NDArray:
    """
    Discretize a scattering macroscopic cross-section from a continuous
    function into a numpy array of len==len(group_structure)-1.

    Parameters
    ----------
    continuous_xs:
        continuous microscopic scattering cross-section function to be discretized.
    group_structure:
        group structure of neutron energies, the n+1 energy bin boundaries for the n
        neutron groups, in descending energies.

    Returns
    -------
    :
        microscopic scattering cross-section values discretized according to the
        group structure. 2D array, where element [i,j] is proportional to the
        probability of scattering neutrons from bin i to bin j.
    """
    return (
        discretize_xs(continuous_xs, group_structure) *
        scattering_weight_matrix(group_structure, atomic_mass).T
    ).T

def discretize_n2n_xs(
    continuous_xs, group_structure: npt.NDArray[np.float64], q_value: float
) -> npt.NDArray:
    """
    Discretize a n2n reaction macroscopic cross-section from a continuous
    function into a numpy array of len==len(group_structure)-1.

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
        microscopic n2n reaction cross-section values discretized according to the
        group structure. 2D array, where element [i,j] is proportional to the
        probability of n2n neutrons being born into bin j due to a reaction
        caused by a neutron in bin i.
    """
    return (
        discretize_xs(continuous_xs, group_structure) * 2 *
        n2n_weight_matrix(group_structure, atomic_mass).T
    ).T

def discretize_fission_xs(
    continuous_xs, group_structure: npt.NDArray[np.float64],
    fission_spectrum_continuous, num_neutrons: float,
) -> npt.NDArray:
    """
    Discretize a fission macroscopic cross-section from a continuous
    function into a numpy array of len==len(group_structure)-1.

    Parameters
    ----------
    continuous_xs:
        continuous cross-section function to be discretized.
    group_structure:
        group structure of neutron energies, the n+1 energy bin boundaries for the n
        neutron groups, in descending energies.
    fission_spectrum_continuous:
        The neutron spectrum formed by the fission neutrons (e.g. the Watt spectrum)
        as a continuous, integrable function.

    Returns
    -------
    :
        microscopic fission cross-section values discretized according to the
        group structure. 2D array, where element [i,j] is proportional to the
        probability of a fission neutron being born into bin j due to a reaction
        caused by a neutron in bin i.
    """
    fiss_spec = integrate(fission_spectrum_continuous, group_structure)

    scaled_fission_spectrum = fiss_spec/fiss_spec.sum() * num_neutrons
    return np.outer(
        discretize_xs(continuous_xs, group_structure),
        scaled_fission_spectrum,
    )

def scattering_weight_matrix(
    group_structure: npt.NDArray[np.float64], atomic_mass: float
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
        An upper triangular matrix, where the i-th row contains the normalized weights for
        scattering down from the i-th bin to the j-th bin. The main-diagonal contains
        the self-scattering cross-section.
        e.g. [2,1] would be the fraction of elastic scattering reactions that causes
        neutrons from group 3 to scatter into group 2.
        The lower triangle must be all zeros.
        np.sum(axis=1) == np.ones(len(group_structure)-1).
    """
    return atomic_mass, group_structure


def n2n_weight_matrix(
    group_structure: npt.NDArray[np.float64], q_value: float
) -> npt.NDArray:
    """
    Parameters
    ----------
    group_structure:
        the n+1 energy bin boundaries for the n neutron groups, in descending energies,
        in eV.
    q_value:
        the q-value of the reaction in eV.

    Returns
    -------
    :
        A macroscopic cross-section matrix, where the j-th column of the i-th row
        expresses the probability of one of the n2n neutrons trigged by a i-th
        bin neutron ends up in the j-th bin.
        np.sum(axis=1) == np.ones(len(group_structure)-1).
    """
    return


@dataclass
class MaterialMacroInfo:
    """
    Material information.

    Parameters
    ----------
    sigma_t:
        total macroscopic cross-section, 1D array of len = n.
    sigma_s:
        Source matrix, 2D array of shape (n, n). Includes a sum of the
        scattering matrix and n2n matrix. Scattering matrix is the
        macroscopic scattering cross-section from group i to j, and the
        n2n matrix is the macroscopic cross-section for production of
        group j neutrons due to group i neutrons. It should be
        mostly-upper-triangular, i.e. the lower triangle must have small values
        compared to the average macroscopic cross-section value of the matrix.

        e.g. [0,3] would be the cross-section for proudction of group 4 neutrons
        due to reactions (scattering and n2n) caused by group 1 neutrons.
    group_structure:
        energy bin edges, 1D array of len = n+1, in eV.
    avg_atomic_mass:
        average atomic mass (weighted by fraction)
    """

    sigma_t: npt.NDArray[np.float64]
    sigma_s: npt.NDArray
    group_structure: npt.NDArray
    avg_atomic_mass: float

    def __post_init__(self):
        """Validation to confirm the shape is correct."""
        # force into float or numpy arrays of floats.
        self.sigma_t = np.array(self.sigma_t, dtype=float)
        self.sigma_s = np.array(self.sigma_s, dtype=float)
        self.group_structure = np.clip(self.group_structure, 1E-9, np.inf)
        self.avg_atomic_mass = float(self.avg_atomic_mass)

        if np.diff(self.group_structure) >= 0:
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

    @property
    def n_groups(self):
        """Store this attribute upon first retrieval."""
        if not hasattr(self, "_n_groups"):
            self._n_groups = len(self.group_structure) - 1
        return self._n_groups

    @property
    def downscatter_only(self):
        """
        Checks if the source matrix suggests that neutrons from each group can
        only cause neutrons to be born in higher-lethargy (lower energy) groups.

        If True, the transport equation can be solved acyclically, from low
        to high lethargy groups in a single-pass.

        If False, and the transport equation will need to be solved iteratively,
        as neutron fluxes in higher-lethargy groups in turn affects the neutron
        flux in lower-lethargy groups.
        """
        return ~self.contains_upscatter

    @property
    def contains_upscatter(self):
        return np.tril(self.sigma_s, k=-1).any()

def get_material_nuclear_data(
    material: str, group_structure: npt.NDArray[np.float64]
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
    :
        MaterialMacroInfo encapsulating the total macroscopic cross-section,
        source (i.e. scattering + n2n) macroscopic cross-section, group structure
        (same as the input group structure) and average atomic mass number.

    Notes
    -----
    When calculating the group-wise diffusion coefficients, the i-th group's scattering
    cross-section :math:`\\Sigma_{s}` is the i-th element on the main diagonal of the
    macroscopic scattering cross-section matrix.
    """
    density = material_density_data_bank[material]
    composition = material_composition_data_bank[material]
    avg_atomic_mass = get_avg_atomic_mass(composition)

    # dicts of {"isotope": npt.NDArray[np.float64] 1D/2D arrays}
    micro_total_xs, micro_scattering_xs, micro_n2n_xs = {}, {}, {}

    for species in composition:
        total_xs_continuous, elastic_scattering_xs_continuous = xs_data_bank[species]

        micro_total_xs[species] = discretize_xs(
            total_xs_continuous, group_structure
        )
        micro_scattering_xs[species] = discretize_scattering_xs(
            elastic_scattering_xs_continuous,
            group_structure,
            extract_atomic_mass(species),

        )
        if species in breeding_xs_data_bank:
            n2n_xs_continuous, q_value = breeding_xs_data_bank[species]
            micro_n2n_xs[species] = discretize_n2n_xs(
                n2n_xs_continuous, group_structure, q_value,
            )
        if species in fission_xs_data_bank:
            fission_xs_continuous, fission_spectrum_continuous, num_neutrons_per_fission = fission_xs_data_bank[species]
            micro_fiss_xs[species] = discretize_fission_xs(
                fission_xs_continuous, group_structure,
                fission_spectrum_continuous, num_neutrons_per_fission,
            )

    discrete_macro_total_xs = calculate_average_macro_xs(
        composition, micro_total_xs, density
    )
    discrete_macro_scattering_xs = calculate_average_macro_xs(
        composition, micro_scattering_xs, density
    )
    discrete_macro_n2n_xs = calculate_average_macro_xs(
        composition, micro_n2n_xs, density
    )
    discrete_macro_fission_xs = calculate_average_macro_xs(
        composition, micro_fiss_xs, density
    )
    source_matrix = discrete_macro_scattering_xs + discrete_macro_n2n_xs + discrete_macro_fission_xs

    return MaterialMacroInfo(
        discrete_macro_total_xs,
        source_matrix,
        group_structure,
        avg_atomic_mass,
    )
