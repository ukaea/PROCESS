"""
Most of derivation is basd on the book Reactor Analysis, Duderstadt and Hamilton, 1976,
ISBN:9780471223634.
"""

import functools

import numpy as np
from numpy import testing as npt
from scipy.constants import Avogadro
from neutronics_data import (ZeroContinuousFunc,
    extract_atomic_mass,
    get_avg_atomic_mass,
    material_density_data_bank,
    material_composition_data_bank,
    xs_data_bank,
    breeding_xs_data_bank
)

BARNS_CM2 = 1e-24
N_A = Avogadro
N2N_Q_VALUE = ...


def groupwise(func):
    """Rename the current func as groupwise_func, and over-write the current method
    as one that sums up all group-wise func.
    """
    method_name = func.__name__
    groupwise_name = f"groupwise_{method_name}"

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        groupwise_func = getattr(self, groupwise_name)
        return np.sum(
            [
                groupwise_func(n, *args, **kwargs)
                for n in range(1, self.n_groups + 1)
            ],
            axis=0,
        )

    def wrapper_setattr(cls):
        """Save the decorated (groupwise) function under a different name after it has
        been created, then overwrite the current method with the np.sum implementation.
        """
        setattr(cls, groupwise_name, func)
        setattr(cls, method_name, wrapper)
        return cls

    # Instead of returning a function, we return a descriptor that registers itself later
    return RegisterLater(wrapper_setattr)


class RegisterLater:
    """Descriptor class"""

    def __init__(self, installer):
        """Modifies the class AFTER it has been created."""
        self.installer = installer

    def __set_name__(self, owner, name):
        """Re-write method name"""
        self.installer(owner)


def get_diffusion_coefficient_and_length(
    total_xs: float, scattering_xs: float, avg_atomic_mass: float
) -> tuple[float, np.complex128]:
    r"""
    Calculate the diffusion coefficient for a given scattering and total macro-scopic
    cross-section in a given medium.

    Parameters
    ----------
    total_xs:
        macroscopic total cross-section `\sigma_{total}`, i.e. any reaction between
        nuclei and neutrons, that either changes the neutron's path or remove it from
        that energy group.
        Unit: cm^-1
    scattering_xs:
        macroscopic total cross-section `\sigma_{scatter}`, i.e. number of reactions per
        unit distance travelled by the neutron that leads to it being scattered (without
        getting absorbed).
        Unit: cm^-1
    avg_atomic_mass:
        Average atomic mass in [amu]. This can be approximated by the atomic number 'A'
        of the medium that the neutron passes through. The effect of the more-anisotropic
        scattering due to smaller atomic number can be accounted for by increasing the
        diffusion coefficient (which decreases the transport macroscopic cross-section,
        :math:`\Sigma_{tr}=\frac{1}{3D}`).

    Returns
    -------
    diffusion_coef:
        The diffusion coefficient as given by Reactor Analysis, Duderstadt and Hamilton.
        unit: [cm]
    diffusion_len:
        The characteristic diffusion length as given by Reactor Analysis, Duderstadt and
        Hamilton.
        unit: [cm]
    """

    transport_xs = total_xs - 2 / (3 * avg_atomic_mass) * scattering_xs
    diffusion_coef = 1 / 3 / transport_xs
    diffusion_len = np.sqrt(
        complex(diffusion_coef / (total_xs - scattering_xs), 0.0)
    )
    return diffusion_coef, diffusion_len


def extrapolation_length(diffusion_coefficient: float) -> float:
    """Get the extrapolation length of the final medium :math:`\\delta`.

    Notes
    -----
    Diffusion theory breaks down at the vacuum boundary, where once the neutron exits,
    it will travel indefinitely into free space, never to return. To counteract this
    problem, we can approximate the neutron profile quite closely by assuming that the
    flux goes to 0 at an extended boundary, rather than at the vacuum boundary.
    THis yields a very close approximation. All of this equation is provided by
    Duderstadt and Hamilton.
    """
    return 0.7104 * 3 * diffusion_coefficient


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
    return


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
    return


def expand_macro_neutron_multiplication_xs_into_matrix(
    discrete_n2n_xs: npt.NDArray[np.float64],
    group_structure: list[float],
    q_value: float,
) -> npt.NDArray:
    """Instead of only having the macroscopic cross-section values for the (n,2n)
    reaction for each group, calculate the macroscopic cross-section of neutron in the
    [i]-th bin producing a neutron in the [j]-bin, recorded as the [i,j] element in the
    matrix.

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


def get_material_nuclear_data(
    material: str, group_structure: list[float]
) -> tuple[npt.NDArray[np.float64], npt.NDArray]:
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
    discrete_macro_scattering_xs_matrix = (
        scattering_weight_matrix(group_structure, avg_atomic_mass).T
        * discrete_macro_scattering_xs
    ).T
    discrete_n2n_xs = calculate_average_macro_xs(
        composition, micro_n2n_xs, density
    )
    discrete_macro_mult_xs_matrix = (
        expand_macro_neutron_multiplication_xs_into_matrix(
            discrete_n2n_xs, group_structure, N2N_Q_VALUE
            # TODO: need to restructure this as we may have more than one types of (n,2n)
            # reactions, requiring different values of N2N_Q_VALUE.
        )
    )
    source_matrix = (
        discrete_macro_scattering_xs_matrix + discrete_macro_mult_xs_matrix
    )

    return discrete_macro_total_xs, source_matrix, avg_atomic_mass


class NeutronFluxProfile:
    """Neutron flux in the first wall or the blanket."""

    def __init__(
        self,
        flux: float,
        fw_mat: str,
        x_fw: float,
        bz_mat: str,
        x_bz: float,
        n_groups: int = 1,
    ):
        """Initialize a particular FW-BZ geometry and neutron flux.

        Parameters
        ----------
        flux:
            Nuetron flux directly emitted by the plasma, incident on the first wall.
            unit: arbitrary, but whichever unit used here will be the same unit used at
            the output neutron_flux.
        fw_mat:
            first wall material
        x_fw:
            thickness of the first wall [m]. It will be converted and stored in [cm].
        bz_mat:
            blanket material
        x_bz:
            thickness of the blanket [m]. It will be converted and stored in [cm].
        n_groups:
            number of groups used to approximate this.
        """
        self.flux = flux
        if not (0 < x_fw < x_bz):
            raise ValueError(
                f"Cannot construct a first-wall+blanket module where{x_fw=}, {x_bz=}."
            )
        self.fw_mat = fw_mat
        self.x_fw = x_fw * 100
        self.bz_mat = bz_mat
        self.x_bz = x_bz * 100
        self.n_groups = n_groups

        self.fw_sigma_t, self.fw_sigma_s, self.fw_A = (
            get_material_nuclear_data(self.fw_mat, n_groups)
        )
        self.bz_sigma_t, self.bz_sigma_s, self.bz_A = (
            get_material_nuclear_data(self.bz_mat, n_groups)
        )
        # Dictionaries indexed by integers, so that we don't have to worry about ordering
        self.integration_constants = {}
        self.extended_boundary = {}
        self.l_fw, self.l_bz = {}, {}

    def solve_one_group(self) -> None:
        i = 0
        if i in self.integration_constants:
            return  # skip if it has already been solved.
        c1 = self.flux * ...
        c2 = self.flux * ...
        c3 = self.flux * ...
        c4 = self.flux * ...
        self.l_fw[i], d_fw = get_diffusion_coefficient_and_length(
            self.fw_sigma_t[i],
            self.fw_sigma_s[i, i],
            self.fw_A,
        )
        self.l_bz[i], d_bz = get_diffusion_coefficient_and_length(
            self.bz_sigma_t[i],
            self.bz_sigma_s[i, i],
            self.bz_A,
        )
        self.extended_boundary[i] = self.x_bz + extrapolation_length(d_bz)
        self.integration_constants[i] = [c1, c2, c3, c4]

    def solve_group_n(self, n: int) -> None:
        if n not in range(1, self.n_groups):
            raise ValueError(
                f"n must be a positive integer between 1 and {self.n_groups}!"
            )
        if n == 1:
            self.solve_one_group()
        for k in range(n - 1):
            if k not in self.integration_constants:
                self.solve_group_n(k)
        i = n - 1
        if i in self.integration_constants:
            return  # skip if it has already been solved.
        c1 = self.flux * ...
        c2 = self.flux * ...
        c3 = self.flux * ...
        c4 = self.flux * ...
        self.l_fw[i], d_fw = get_diffusion_coefficient_and_length(
            self.fw_sigma_t[i],
            self.fw_sigma_s[i, i],
            self.fw_A,
        )
        self.l_bz[i], d_bz = get_diffusion_coefficient_and_length(
            self.bz_sigma_t[i],
            self.bz_sigma_s[i, i],
            self.bz_A,
        )
        self.extended_boundary[i] = self.x_bz + extrapolation_length(d_bz)
        self.integration_constants[i] = [c1, c2, c3, c4]

    @groupwise
    def neutron_flux_fw(self, n: int, x: float | npt.NDArray) -> npt.NDArray:
        """Neutron flux at the first wall."""
        i = n - 1
        c1, c2 = self.integration_constants[i][:2]
        return np.real(
            c1 * np.sinh(x / self.l_fw[i]) + c2 * np.cosh(x / self.l_fw[i])
        )

    @groupwise
    def neutron_flux_bz(self, n: int, x: float | npt.NDArray) -> npt.NDArray:
        """Neutron flux at the blanket."""
        i = n - 1
        c3, c4 = self.integration_constants[i][2:]
        return np.real(
            c3 * np.sinh(x / self.l_bz[i]) + c4 * np.cosh(x / self.l_bz[i])
        )

    @groupwise
    def neutron_flux_at(self, n: int, x: float | npt.NDArray) -> npt.NDArray:
        """
        Neutron flux
        Parameters
        ----------
        n:
            Neutron group number
        x:
            The depth where we want the neutron flux (m).
        """
        if np.isscalar(x):
            return self.groupwise_neutron_flux_at(n, [x])[0]
        x = np.asarray(x)
        in_fw = abs(x) <= self.x_fw
        in_bz = np.logical_and(self.x_fw < abs(x), abs(x) <= self.x_bz)
        if (~np.logical_or(in_fw, in_bz)).any():
            raise ValueError(
                f"for neutron group {n}, neutron flux can only be calculated up to "
                f"{self.extended_boundary[n]}, which {x} violates!"
            )

        out_flux = np.zeros_like(x)
        out_flux[in_fw] = self.groupwise_neutron_flux_fw(x[in_fw])
        out_flux[in_bz] = self.groupwise_neutron_flux_bz(x[in_bz])
        return out_flux
