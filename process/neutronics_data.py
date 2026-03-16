from abc import ABC, abstractmethod, abstractproperty
import warnings
import json
from dataclasses import dataclass
from itertools import islice, pairwise
from pathlib import Path

import numpy as np
from numpy import typing as npt
from scipy.constants import Avogadro

from process.exceptions import ProcessValidationError

BARNS_TO_M2 = 1e-28
N_A = Avogadro
EV_TO_J = 1.602e-19
DT_NEUTRON_E = 14.06e6 * EV_TO_J

with open(Path(__file__).parent / "data" / "atomic_data.json") as j:
    _data = json.load(j)
    _ATOMIC_MASS_SOURCE = _data["ATOMIC_MASS_SOURCE"]
    ATOMIC_MASS = _data["ATOMIC_MASS"]
    _NATURAL_ABUNDANCE_SOURCE = _data["NATURAL_ABUNDANCE_SOURCE"]
    NATURAL_ABUNDANCE = _data["NATURAL_ABUNDANCE"]

def get_isotopic_composition(element: str):
    """Get the isotopic composition of a naturally occurring element."""
    return {iso:frac for iso,frac in NATURAL_ABUNDANCE.items() if iso.startswith(element)}


def elem_to_isotopic_comp(composition: dict[str, float]) -> dict[str, float]:
    """
    Convert an element composition dictionary into an isotope composition
    dictionary using NATURAL_ABUNDANCE.
    Parameters
    ----------
    composition:
        A dictionary showing the atomic fraction that each relevant element
        makes up.

    Returns
    -------
    new_comp_dict:
        A dictionary showing the atomic fraction that each relevant isotope
        makes up.
    """
    new_comp_dict = {}
    for elem, overall_fraction in composition.items():
        for isotope, fraction in get_isotopic_composition(elem).items():
            new_comp_dict[isotope] = overall_fraction * fraction
    return new_comp_dict


def get_avg_atomic_mass(composition: dict[str, float]) -> float:
    """Calculate the average atomic mass number.
    Parameters
    ----------
    composition:
        A dictionary showing the atomic fraction that each species makes up.
    """
    total_fraction = sum(composition.values())
    return sum(
        ATOMIC_MASS[species] * (fraction / total_fraction)
        for species, fraction in composition.items()
    )

class ENDFRecord(ABC):
    @abstractmethod
    def get_xs(
            self, mt: int, group_structure: npt.NDArray[np.float64]
        ) -> npt.NDArray[np.float64]:
        """
        A method that returns the formatted (integrated, weighted by lethargy)
        cross-section.
        """
        pass

@dataclass
class MTTreeNode:
    reaction_name: str
    main_number: int
    constituents: tuple[int]
    auxillary: dict[int, str]

    def resolve_xs(self, endf_record: ENDFRecord, group_structure: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """
        Check if the endf record is sesnsible by following the correct summation rules,
        and then extract the correct cross-section.

        Parameters
        ----------
        endf_record:
            an instance of a class with the method with a signature of (self, group_structure)
        group_structure:
            A 1D array of numbers representing the group boundaries, from high
            energy to low energy (low lethargy to high lethargy).

        Returns
        -------
        xs:
            formatted to the correct group structure.
        """
        main_xs = endf_record.get_xs(main_number, group_structure)
        constituent_xs = np.sum([endf_record.get_xs(mt, group_structure) for mt in constituents], axis=0)
        if np.sum(main_xs) < np.sum(constituent_xs):
            raise ValueError(f"MT={main_number} does not match the constituents")
        if len(main_xs)!=(len(group_structure)-1):
            raise ValueError(f"{self.mt_tuple} resolves to zero!")
        return main_xs

class MTResolutionRule():
    def __init__(self, mt_tuple: tuple[MTTreeNode], redundant_reaction: str, redundant_mt: int):
        self.redundant_reaction = str(redundant_reaction)
        self.redundant_mt = int(redundant_mt)
        self.mt_tuple = mt_tuple

    def resolve_xs(self, endf_record: ENDFRecord, group_structure: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """
        Check if the endf record is sesnsible by following the correct summation rules,
        and then extract the correct cross-section.

        Parameters
        ----------
        endf_record:
            an instance of a class with the method with a signature of (self, group_structure)
        group_structure:
            A 1D array of numbers representing the group boundaries, from high
            energy to low energy (low lethargy to high lethargy).

        Returns
        -------
        xs:
            formatted to the correct group structure.
        """
        redundant_xs = endf_record.get_xs(redundant_mt, group_structure)
        main_xs = np.sum([mt.resolve_xs(endf_record, group_structure) for mt in self.mt_tuple], axis=0)
        if redundant_xs:
            if np.sum(redundant_xs) < np.sum(main_xs):
                raise ValueError(f"mt={redundant_mt} failed to include everything.")
            return redundant_xs
        return main_xs

MT_N2N = MTTreeNode("n,2n", 16, range(875, 892), {11: "n,2nd", 21: "n,2nf", 24: "n,2na", 30: "n,2n2a", 41: "n,2np"}, )
MT_N3N = MTTreeNode("n,3n", 17, (), {25: "n,3na", 38: "n,3nf", 42: "n,3np"})
MT_N4N = MTTreeNode("n,4n", 37, (), {}, )
MT_TRITON = MTTreeNode("n,t", 105, range(700, 750), {33: "n,nt", 113:"n,t2a", 116:"n,pt", })
MT_SCAT_INELASTIC = MTTreeNode("n,n'", 4, range(50, 92),
    {20: "n,nf", 22: "n,na", 23: "n,n3a", 28: "n,np", 29: "n,n2a", 32: "n,nd",
    33: "n,nt", 34: "n,n3He", 35: "n,nd2a", 36: "n,nt2a", 44: "n,n2p",
    45: "n,npa"}
)


class ExtractedNuclearData:
    """A class for storing the useful nuclear data extracted out of the ENDFRecord."""
    rules = {
        "total" : MTResolutionRule((), "n,tot", 1),
        "elastic_scattering" : MTResolutionRule((), "n,n", 2),
        "inelastic_scattering" : MTResolutionRule((MT_SCAT_INELASTIC,), "n,n'", 4),
        "neutron_producing" : MTResolutionRule((MT_N2N, MT_N3N, MT_N4N), "n,Xn", 201),
        "triton_producing" : MTResolutionRule((MT_TRITON,), "n,Xt", 205),
    }
    def __init__(self,
        group_structure: npt.NDArray[np.float64],
        atomic_mass: float,
        endf_record: ENDFRecord,
    ):
        self.group_structure = group_structure
        self.atomic_mass = atomic_mass
        self.endf_record = endf_record

        self.sigma_total = self.rules["total"].resolve(self.endf_record, self.group_structure)
        self.sigma_triton = self.rules["triton_producing"].resolve(self.endf_record, self.group_structure)

        self.sigma_scatter = (
            self.rules["elastic_scattering"].resolve(
                self.endf_record, self.group_structure
            )
            * scattering_weight_matrix(
                self.group_structure, self.atomic_mass
            ).T
        ).T
        self.sigma_scatter +=(
            self.rules["inelastic_scattering"].resolve(
                self.endf_record, self.group_structure
            ) * scattering_weight_matrix(
                self.group_structure, self.atomic_mass,
                self.endf_record.q_values.get("inelastic_scattering", 0.0)
            ).T  # TODO: calculate with an inelastic scattering weight matrix instead.
        ).T

        neutron_producing = self.rules["neutron_producing"].resolve(self.endf_record, self.group_structure)
        if np.sum(neutron_producing) == 0:
            return  # skip the whole thing below to save time if 0.
        n2n = self.endf_record.get_xs(16, self.group_structure)
        n3n = self.endf_record.get_xs(17, self.group_structure)
        n4n = self.endf_record.get_xs(37, self.group_structure)
        all_n2n = MT_N2N.resolve(self.endf_record, self.group_structure)
        all_n3n = MT_N3N.resolve(self.endf_record, self.group_structure)
        all_n4n = MT_N4N.resolve(self.endf_record, self.group_structure)
        mt201 = self.endf_record.get_xs(201, self.group_structure)
        if np.isclose(n2n + n3n + n4n, neutron_producing, atol=0.0):
            self.sigma_in = (
                n2n * nXn(self.group_structure, self.endf_record.q_values.get("n,2n", 0.0), 2).T
                +n3n * nXn(self.group_structure, self.endf_record.q_values.get("n,3n", 0.0), 3).T
                +n4n * nXn(self.group_structure, self.endf_record.q_values.get("n,4n", 0.0), 4).T
            ).T
        elif np.isclose(n2n + n3n + n4n, neutron_producing, atol=0.0):
            self.sigma_in = (
                all_n2n * nXn(self.group_structure, self.endf_record.q_values.get("n,2n", 0.0), 2).T
                +all_n3n * nXn(self.group_structure, self.endf_record.q_values.get("n,3n", 0.0), 3).T
                +all_n4n * nXn(self.group_structure, self.endf_record.q_values.get("n,4n", 0.0), 4).T
            ).T
        else:
            self.sigma_in = (
                neutron_producing * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get("n,2n", 0.0), 2).T
            ).T
        self.sigma_in = (
            _neutron_producing.resolve(self.endf_record, self.group_structure)
            * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get("n,2n", 0.0), 2).T
        ).T

material_density_data_bank = {"stainless steel": {"Fe": 6.0}}
material_composition_data_bank = {"stainless steel": {"Fe60": 1.0}}
xs_data_bank = ...
breeding_xs_data_bank = {}
fission_xs_data_bank = {}

def read_mass(iso_name: str):
    length = len(iso_name)
    name, mass = [], []
    for i in range(length):
        name.append(iso_name[i])
        if iso_name[i+1].isnumeric():
            break
    for j in range(i+1, length):
        mass.append(iso_name[j])
    return "".join(name), int("".join(mass))

def is_natural(iso_name: str):
    iso_symbol, iso_num = read_mass(iso_name)
    return iso_num==0



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


def discretize_nXn_xs(
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
        * nXn_weight_matrix(group_structure, q_value, X).T
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
    group_structure: npt.NDArray[np.float64], atomic_mass: float, energy_lost=0.0
) -> npt.NDArray:
    """
    Parameters
    ----------
    group_structure:
        the n+1 energy bin boundaries for the n neutron groups, in descending energies.
    atomic_mass:
        atomic mass of the medium
    energy_lost:
        The energy lost in J. Currently a placeholder variable that isn't used.
        TODO: make the inelastic scattering weight matrix account for this loss.

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


def nXn_weight_matrix(
    group_structure: npt.NDArray[np.float64], q_value: float, X: int,
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
        1, i.e. the number of neutrons.
    """
    # Assume that the two neutrons would share the resulting energy evenly, i.e.
    # each take half of the neutron
    # To make things even simpler, we'll assume the neutron flux is
    shift_e = -q_value / X

    e_i1, e_i = group_structure[:-1], group_structure[1:]
    weight = 1 / (np.log(e_i1) - np.log(e_i))

    n_groups = len(group_structure) - 1
    e_g1 = np.broadcast_to(e_i1, [n_groups, n_groups])
    e_g = np.broadcast_to(e_i, [n_groups, n_groups])

    e_min = np.clip((e_g + shift_e).T, e_i, e_i1).T
    e_max = np.clip((e_g1 + shift_e).T, e_min.T, e_i1).T

    matrix = np.log(e_max) - np.log(e_min)
    return (weight * matrix.T).T


class MaterialMacroInfo():
    def __init__(self,
        group_structure: npt.NDArray[np.float64],
        density: float,
        elements: dict[str, float],
        name: str="",
        source: str="",
        comment: str="",
    ):
        """
        Parameters
        ----------
        group_structure:
            energy bin edges, 1D array of len = n+1, in [J].
        density:
            density of the material in kg/m3
        name:
            name of the material (optional)
        source:
            data source described by a string.
        comment:
            comment on the data (string)
        """
        self.group_structure = np.asarray(group_structure)
        if (np.diff(self.group_structure) >= 0).any():
            raise ValueError(
                "The group structure must be defined descendingly, from the "
                "highest energy bin (i.e. lowest lethargy bin) edge to the "
                "lowest energy bin edge, which can't be zero (infinite "
                "lethargy). Similarly the cross-section must be arranged "
                "according to these bin edges."
            )
        if (self.group_structure <= 0).any():
            warnings.warn("Zero energy (inf. lethargy) not allowed.")
            self.group_structure = np.clip(
                self.group_structure, 1e-9 * EV_TO_J, np.inf
            )
        self.density = float(density)
        self.name = name
        self._populated = False
        self.elements = elements
        self.avg_atomic_mass = get_avg_atomic_mass(elem_to_isotopic_comp(self.elements))
        self.number_density = N_A / self.avg_atomic_mass * 1000

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

    def __repr__(self):
        return super().__repr__().replace(" at ", f" '{self.name}' at ")

    def _set_sigma(self, sigma_t, sigma_s, sigma_in=None, sigma_triton=None) -> None:
        """Populate the values directly. Mainly to make unit-testing easier."""
        self._sigma_total = np.asarray(sigma_t)
        self._sigma_scatter = np.asarray(sigma_s)
        if sigma_in is not None:
            self._sigma_in = np.asarray(sigma_in)
        else:
            self._sigma_in = np.zeros([self.n_groups, self.n_groups])
        if sigma_triton is not None:
            self._sigma_triton = np.asarray(sigma_triton)
        else:
            self._sigma_triton = np.zeros(self.n_groups)
        self._confirm_sigma()

    def _confirm_sigma(self) -> None:
        if np.shape(self._sigma_total) != (self.n_groups,):
            raise ProcessValidationError(
                f"total group-wise cross-sections should have {self.n_groups} "
                "groups as specified by the group_structure."
            )
        if np.shape(self._sigma_scatter) != (self.n_groups, self.n_groups):
            raise ProcessValidationError(
                "Group-wise scattering cross-sections be a square matrix of "
                f"shape n*n, where n= number of groups = {self.n_groups}."
            )
        if (self._sigma_scatter.sum(axis=1) > self._sigma_total).any():
            raise ProcessValidationError(
                "Total cross-section should include the scattering cross-section."
            )
        if np.tril(self._sigma_scatter, k=-1).any():
            warnings.warn(
                "Elastic up-scattering seems unlikely in this model! "
                "Check if the group structure is chosen correctly?",
            )
        self._populated = True
        return

    def populate_from_data_library(self, xs_dict: [str, ExtractedNuclearData]) -> None:
        """Populate the cross-section values according to the nuclear data library."""
        self._sigma_total = np.zeros(self.n_groups)
        self._sigma_scatter = np.zeros([self.n_groups, self.n_groups])
        self._sigma_in = np.zeros([self.n_groups, self.n_groups])
        self._sigma_triton = np.zeros(self.n_groups)
        for element, elem_frac in self.elements.items():
            natural = element + "0"
            if natural in xs_dict:
                self._add_data_from_single_record(
                    xs_dict[natural], elem_frac * self.number_density
                )
            else:
                for isotope, natural_abundance in get_isotopic_composition():
                    self._add_data_from_single_record(
                        xs_dict[isotope], natural_abundance * elem_frac * self.number_density
                    )
        self._confirm_sigma()
        return

    def _add_data_from_single_record(
            self, xs_data: ExtractedNuclearData, partial_number_density: float
        ) -> None:
        if not np.isclose(
            self.group_structure, xs_data.group_structure, atol=0.0
        ):
            raise ValueError(f"Mismatched group structure with {xs_data}.")
        self._sigma_total += xs_data.sigma_total * partial_number_density * BARNS_TO_M2
        self._sigma_scatter += xs_data.sigma_scatter * partial_number_density * BARNS_TO_M2
        self._sigma_in += xs_data.sigma_in * partial_number_density * BARNS_TO_M2
        self._sigma_triton += xs_data.sigma_triton * partial_number_density * BARNS_TO_M2

    @property
    def sigma_t(self):
        """total macroscopic cross-section, 1D array of len = n."""
        if not self._populated:
            raise ValueError("Empty cross-section data!")
        return self._sigma_total

    @property
    def sigma_s(self):
        """
        Macroscopic scattering cross-section from group i to j, forming a 2D
        array of shape (n, n). It should be mostly-upper-triangular, i.e. the
        lower triangle (excluding the main diagonal) must have small values
        compared to the average macroscopic cross-section value of the matrix.
        Neutrons fluxes are assumed to be isotropic before and after scattering.

        e.g. [0,3] would be the cross-section for group 4 neutrons scattered-in
        from group 1 neutrons.
        """
        if not self._populated:
            raise ValueError("Empty cross-section data!")
        return self._sigma_scatter

    @property
    def sigma_in(self):
        """
        In-source matrix: for now, it includes a sum of the matrix of (n,2n)
        reactions and fission reactions. Same logic as the scattering matrix,
        i.e. probability of group j neutrons produced (presumed to be
        isotropic) per unit flux of group i neutrons.

        e.g. [0,3] would be the cross-section for the proudction of group 4
        neutrons due to n,2n and fission reactions caused by group 1 neutrons.
        """
        if not self._populated:
            raise ValueError("Empty cross-section data!")
        return self._sigma_in

    @property
    def sigma_triton(self):
        if not self._populated:
            raise ValueError("Empty cross-section data!")
        return self._sigma_triton

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

    @property
    def element_set(self):
        if not hasattr(self, "_element_set"):
            self._element_set = set(self.elements.keys())
        return self._element_set

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
            ATOMIC_MASS[species],
        )
        if species in breeding_xs_data_bank:
            n2n_xs_continuous, q_value = breeding_xs_data_bank[species]
            micro_n2n_xs[species] = discretize_nXn_xs(
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
