import warnings
from dataclasses import dataclass
from itertools import islice, pairwise
from pathlib import Path

import numpy as np
from numpy import typing as npt
from scipy.constants import Avogadro

from process.exceptions import ProcessValidationError

BARNS_TO_M2 = 1e-28
N_A = Avogadro
N2N_Q_VALUE = ...
_ATOMIC_MASS = {}
EV_TO_J = 1.602e-19
DT_NEUTRON_E = 14.06e6 * EV_TO_J


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
        Density of the medium, given in [kg/m^3]

    Returns
    -------
    :
        macroscopic cross-section of the material. [m]

    Notes
    -----
    .. math::
        \Sigma &= N_d\sigma
        &= \frac{N_A}{A} \rho \sigma
    """

    total_fraction = sum(composition.values())
    weighted_micro_xs, weighted_atomic_mass = [], []
    if extra_species := set(micro_xs.keys()).difference(
        set(composition.keys())
    ):
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

    return N_A / (avg_mass_amu / 1000) * density * (avg_sigma * BARNS_TO_M2)


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
    lethargy = -np.log(group_structure)
    lower_leth, upper_leth = lethargy[:-1], lethargy[1:]
    leth_spacing = np.diff(lethargy)

    def continuous_xs_lethargy_input(lethargy):
        return continuous_xs(np.exp(-lethargy))

    return (
        integrate(continuous_xs_lethargy_input, lower_leth, upper_leth)
        / leth_spacing
    )


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
        discretize_xs(continuous_xs, group_structure)
        * scattering_weight_matrix(group_structure, atomic_mass).T
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
        discretize_xs(continuous_xs, group_structure)
        * n2n_weight_matrix(group_structure, q_value).T
    ).T


def discretize_fission_xs(
    continuous_xs,
    group_structure: npt.NDArray[np.float64],
    fission_spectrum_continuous,
    num_neutrons: float,
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
    e_max, e_min = group_structure[:-1], group_structure[1:]
    fiss_spec = integrate(fission_spectrum_continuous, e_max, e_min)

    scaled_fission_spectrum = fiss_spec / fiss_spec.sum() * num_neutrons
    return np.outer(
        discretize_xs(continuous_xs, group_structure),
        scaled_fission_spectrum,
    )


def _get_alpha(atomic_mass: float):
    return ((atomic_mass - 1) / (atomic_mass + 1)) ** 2


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

    alpha = _get_alpha(atomic_mass)
    n_groups = len(group_structure) - 1
    matrix = np.zeros([n_groups, n_groups], dtype=float)
    for i, in_group in enumerate(pairwise(group_structure)):
        for g, out_group in enumerate(pairwise(group_structure)):
            if i == g:
                for energy_limits, case_descr in _split_into_energy_limits(
                    alpha, in_group
                ):
                    matrix[i, g] += _convolved_scattering_fraction(
                        energy_limits, alpha, in_group, out_group, case_descr
                    )
            elif i < g:
                for energy_limits, case_descr in _split_into_energy_limits(
                    alpha, in_group, out_group
                ):
                    matrix[i, g] += _convolved_scattering_fraction(
                        energy_limits, alpha, in_group, out_group, case_descr
                    )
    return matrix


def _split_into_energy_limits(
    alpha,
    in_group: tuple[float, float],
    out_group: tuple[float, float] | None = None,
):
    """
    Spit out a list of integration limits to be used by _convolved_scattering_fraction.

    Parameters
    ----------
    alpha:
        Parameters
    in_group:
        Descending energy bounds of the input group. Scattering neutrons are
        generated from this energy group. (e_i1, e_i)
    out_group:
        Descending energy bounds of the output group. Scattering neutrons are
        born into this energy group. (e_g1, e_g)

    Returns
    -------
    :
        list of 2-tuples, each two tuple contains the following items:
    energy_limits:
        A tuple of two floats in descending order, describing the integration
        limits.
    case_descr:
        A string to describe which case it belongs to.
    """
    e_i1, e_i = in_group
    max_min_scattered_energy = alpha * e_i1
    min_min_scattered_energy = alpha * e_i

    if not out_group:
        # min_min_scattered_energy < e_i < max_min_scattered_energy < e_i1
        if (alpha > 0) and (e_i < max_min_scattered_energy):
            mid_e = min(e_i / alpha, e_i1)
            return [
                ((e_i1, mid_e), "self, complete"),
                ((mid_e, e_i), "self, upper-half"),
            ]
        # min_min_scattered_energy < max_min_scattered_energy < e_i < e_i1
        return [((e_i1, e_i), "self, upper-half")]
    e_g1, e_g = out_group
    if alpha == 0:
        # min_min_scattered_energy < max_min_scattered_energy < e_g < e_g1
        return [((e_i1, e_i), "down, middle")]

    top_e = min(e_g1 / alpha, e_i1)
    mid_e = min(e_g / alpha, e_i1)

    # e_g < e_g1 < min_min_scattered_energy < max_min_scattered_energy
    if e_g1 <= min_min_scattered_energy:
        return []

    # e_g < min_min_scattered_energy < e_g1 < max_min_scattered_energy
    # e_g < min_min_scattered_energy < max_min_scattered_energy < e_g1
    if e_g <= min_min_scattered_energy:
        return [((top_e, e_i), "down, upper-half")]

    # min_min_scattered_energy < e_g < e_g1 < max_min_scattered_energy
    # min_min_scattered_energy < e_g < max_min_scattered_energy < e_g1
    if e_g < max_min_scattered_energy:
        return [
            ((top_e, mid_e), "down, upper-half"),
            ((mid_e, e_i), "down, middle"),
        ]

    # min_min_scattered_energy < max_min_scattered_energy < e_g < e_g1
    return [((top_e, e_i), "down, middle")]


def _convolved_scattering_fraction(
    energy_limits: tuple[float, float],
    alpha: float,
    in_group: tuple[float, float],
    out_group: tuple[float, float],
    case_description: str,
):
    r"""
    The fraction of neutron flux from bin i that would get scattered into bin g
    is calculated as $M_{ig} = \int_{E_{min}}^{E_{max}} flux_i(E) dist(E) dE$,
    where $flux_i=$ normalized flux (we assume the neutron flux in bin i has
    constant in per-unit-lethargy space), i.e. follows a 1/E distribution.
    Hence after intergration, frac is always accompanied by a factor of
    $1/(ln(E_{i-1}) - ln(E_i))$.
    """
    e_max, e_min = energy_limits
    e_i1, e_i = in_group
    e_g1, e_g = out_group

    _const = 1 / (np.log(e_i1) - np.log(e_i))
    _am1i = 1 / (1 - alpha)
    _diff_log_e = np.log(e_max) - np.log(e_min)
    _diff_inv_e = 1 / e_min - 1 / e_max
    match case_description:
        case "self, complete":
            return _const * _diff_log_e
        case "self, upper-half":
            return _const * _am1i * (_diff_log_e - e_g * _diff_inv_e)
        case "down, upper-half":
            return _const * _am1i * -(alpha * _diff_log_e - e_g1 * _diff_inv_e)
        case "down, middle":
            return _const * _am1i * (e_g1 - e_g) * _diff_inv_e


def n2n_weight_matrix(
    group_structure: npt.NDArray[np.float64], q_value: float
) -> npt.NDArray:
    """
    Parameters
    ----------
    group_structure:
        the n+1 energy bin boundaries for the n neutron groups, in descending energies,
        in J.
    q_value:
        the q-value of the reaction in J.

    Returns
    -------
    :
        A macroscopic cross-section matrix, where the j-th column of the i-th row
        expresses the probability of one of the n2n neutrons trigged by a i-th
        bin neutron ends up in the j-th bin.
        np.sum(axis=1) <= np.ones(len(group_structure)-1) * 2, i.e. the
        probability distribution in each row (i.e. i-th bin) is normalized to
        2, i.e. the number of neutrons.
    """
    # Assume that the two neutrons would share the resulting energy evenly, i.e.
    # each take half of the neutron
    # To make things even simpler, we'll assume the neutron flux is
    shift_e = -q_value / 2

    e_i1, e_i = group_structure[:-1], group_structure[1:]
    weight = 1 / (np.log(e_i1) - np.log(e_i))

    n_groups = len(group_structure) - 1
    e_g1 = np.broadcast_to(e_i1, [n_groups, n_groups])
    e_g = np.broadcast_to(e_i, [n_groups, n_groups])

    e_min = np.clip((e_g + shift_e).T, e_i, e_i1).T
    e_max = np.clip((e_g1 + shift_e).T, e_min.T, e_i1).T

    matrix = np.log(e_max) - np.log(e_min)
    return (weight * matrix.T).T


@dataclass
class MaterialMacroInfo:
    """
    Material information.

    Parameters
    ----------
    sigma_t:
        total macroscopic cross-section, 1D array of len = n.
    sigma_s:
        Macroscopic scattering cross-section from group i to j, forming a 2D
        array of shape (n, n). It should be mostly-upper-triangular, i.e. the
        lower triangle (excluding the main diagonal) must have small values
        compared to the average macroscopic cross-section value of the matrix.
        Neutrons fluxes are assumed to be isotropic before and after scattering.

        e.g. [0,3] would be the cross-section for group 4 neutrons scattered-in
        from group 1 neutrons.
    sigma_in:
        In-source matrix: for now, it includes a sum of the matrix of (n,2n)
        reactions and fission reactions. Same logic as the scattering matrix,
        i.e. probability of group j neutrons produced (presumed to be
        isotropic) per unit flux of group i neutrons.

        e.g. [0,3] would be the cross-section for the proudction of group 4
        neutrons due to n,2n and fission reactions caused by group 1 neutrons.
    group_structure:
        energy bin edges, 1D array of len = n+1, in [J].
    avg_atomic_mass:
        average atomic mass (weighted by fraction)
    """

    group_structure: npt.NDArray
    avg_atomic_mass: float
    sigma_t: npt.NDArray[np.float64]
    sigma_s: npt.NDArray
    sigma_in: npt.NDArray | None = None
    name: str = ""

    def __post_init__(self):
        """
        Validation of group_structure, sigma_s and sigma_t to confirm their
        shapes are correct.
        """
        # force into float or numpy arrays of floats.
        self.group_structure = np.array(self.group_structure, dtype=float)
        self.avg_atomic_mass = float(self.avg_atomic_mass)
        self.sigma_t = np.array(self.sigma_t, dtype=float)
        self.sigma_s = np.array(self.sigma_s, dtype=float)
        if self.sigma_in:
            self.sigma_in = np.array(self.sigma_in, dtype=float)
        else:
            self.sigma_in = np.zeros_like(self.sigma_s)

        if (self.group_structure <= 0).any():
            warnings.warn("Zero energy (inf. lethargy) not allowed.")
            self.group_structure = np.clip(
                self.group_structure, 1e-9 * EV_TO_J, np.inf
            )
        if (np.diff(self.group_structure) >= 0).any():
            raise ValueError(
                "The group structure must be defined descendingly, from the "
                "highest energy bin (i.e. lowest lethargy bin) edge to the "
                "lowest energy bin edge, which can't be zero (infinite "
                "lethargy). Similarly the cross-section must be arranged "
                "according to these bin edges."
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
                "Total cross-section should include the scattering cross-section."
            )
        if np.tril(self.sigma_s, k=-1).any():
            warnings.warn(
                "Elastic up-scattering seems unlikely in this model! "
                "Check if the group structure is chosen correctly?",
            )

    @property
    def n_groups(self):
        """
        Number of groups in the group structure.
        Store this attribute upon first retrieval.
        """
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
        return ~(
            np.tril(self.sigma_s, k=-1).any()
            or np.tril(self.sigma_in, k=-1).any()
        )


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
    density = material_density_data_bank[material]  # [kg/m^3]
    composition = material_composition_data_bank[material]
    avg_atomic_mass = get_avg_atomic_mass(composition)
    n_groups = len(group_structure) - 1

    # dicts of {"isotope": npt.NDArray[np.float64] 1D/2D arrays}
    micro_total_xs = {}
    micro_scattering_xs = {}
    micro_n2n_xs = {}
    micro_fiss_xs = {}

    for species in composition:
        total_xs_continuous, elastic_scattering_xs_continuous = xs_data_bank[
            species
        ]

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
                n2n_xs_continuous,
                group_structure,
                q_value,
            )
        if species in fission_xs_data_bank:
            (
                fission_xs_continuous,
                fission_spectrum_continuous,
                num_neutrons_per_fission,
            ) = fission_xs_data_bank[species]
            micro_fiss_xs[species] = discretize_fission_xs(
                fission_xs_continuous,
                group_structure,
                fission_spectrum_continuous,
                num_neutrons_per_fission,
            )

    discrete_macro_total_xs = calculate_average_macro_xs(
        composition, micro_total_xs, density
    )
    discrete_macro_scattering_xs = calculate_average_macro_xs(
        composition, micro_scattering_xs, density
    )
    source_matrix = np.zeros([n_groups, n_groups], dtype=float)
    if micro_n2n_xs:
        source_matrix += calculate_average_macro_xs(
            composition, micro_n2n_xs, density
        )
    if micro_fiss_xs:
        source_matrix += calculate_average_macro_xs(
            composition, micro_fiss_xs, density
        )

    return MaterialMacroInfo(
        group_structure,
        avg_atomic_mass,
        discrete_macro_total_xs,
        discrete_macro_scattering_xs,
        source_matrix,
    )
