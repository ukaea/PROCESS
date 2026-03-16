from abc import ABC, abstractmethod, abstractproperty
import warnings
import json
from dataclasses import dataclass
from itertools import islice, pairwise
from typing import Iterable, Callable
from pathlib import Path

import numpy as np
from numpy import typing as npt
from scipy.constants import Avogadro

try:
    from process.exceptions import ProcessValidationError
except ImportError as e:
    ProcessValidationError = ValueError

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
    return {iso:frac for iso,frac in NATURAL_ABUNDANCE.items() if read_elem_and_mass(iso)[0]==element}


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
        main_xs = endf_record.get_xs(self.main_number, group_structure)
        constituent_xs = np.sum([endf_record.get_xs(mt, group_structure) for mt in self.constituents], axis=0)
        if np.sum(main_xs) < np.sum(constituent_xs) and not np.isclose(np.sum(main_xs), np.sum(constituent_xs), atol=0).all():
            raise ValueError(f"MT={self.main_number} does not match the constituents")
        if len(main_xs)!=(len(group_structure)-1):
            raise ValueError(f"{self.mt_tuple} resolves to zero!")
        auxillary = np.sum([endf_record.get_xs(mt, group_structure) for mt in self.auxillary], axis=0)
        return main_xs + auxillary


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
        redundant_xs = endf_record.get_xs(self.redundant_mt, group_structure)
        main_xs = np.sum([mt.resolve_xs(endf_record, group_structure) for mt in self.mt_tuple], axis=0)
        if np.sum(redundant_xs):
            if np.sum(redundant_xs) < np.sum(main_xs) and not np.isclose(np.sum(redundant_xs), np.sum(main_xs), atol=1E-10, rtol=0.01):
                raise ValueError(f"mt={self.redundant_mt} failed to include everything.")
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
        self.group_structure = np.asarray(group_structure)
        self.atomic_mass = atomic_mass
        self.endf_record = endf_record

        self.sigma_total = self.rules["total"].resolve_xs(self.endf_record, self.group_structure)
        self.sigma_triton = self.rules["triton_producing"].resolve_xs(self.endf_record, self.group_structure)

        self.sigma_scatter = (
            self.rules["elastic_scattering"].resolve_xs(
                self.endf_record, self.group_structure
            )
            * scattering_weight_matrix(
                self.group_structure, self.atomic_mass
            ).T
        ).T
        self.sigma_scatter +=(
            self.rules["inelastic_scattering"].resolve_xs(
                self.endf_record, self.group_structure
            ) * scattering_weight_matrix(
                self.group_structure, self.atomic_mass,
                self.endf_record.q_values.get(4, 0.0)
            ).T  # TODO: calculate with an inelastic scattering weight matrix instead.
        ).T

        neutron_producing = self.rules["neutron_producing"].resolve_xs(self.endf_record, self.group_structure)
        if np.sum(neutron_producing) == 0:
            self.sigma_in = np.zeros(len(self.group_structure)-1)
            return
        n2n = self.endf_record.get_xs(16, self.group_structure)
        n3n = self.endf_record.get_xs(17, self.group_structure)
        n4n = self.endf_record.get_xs(37, self.group_structure)
        all_n2n = MT_N2N.resolve_xs(self.endf_record, self.group_structure)
        all_n3n = MT_N3N.resolve_xs(self.endf_record, self.group_structure)
        all_n4n = MT_N4N.resolve_xs(self.endf_record, self.group_structure)
        mt201 = self.endf_record.get_xs(201, self.group_structure)
        if np.isclose(n2n + n3n + n4n, neutron_producing, atol=0.0).all():
            self.sigma_in = (
                n2n * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get(16, 0.0), 2).T
                +n3n * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get(17, 0.0), 3).T
                +n4n * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get(37, 0.0), 4).T
            ).T
        elif np.isclose(n2n + n3n + n4n, neutron_producing, atol=0.0).all():
            self.sigma_in = (
                all_n2n * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get(16, 0.0), 2).T
                +all_n3n * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get(17, 0.0), 3).T
                +all_n4n * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get(37, 0.0), 4).T
            ).T
        else:
            self.sigma_in = (
                neutron_producing * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get(16, 0.0), 2).T
            ).T
        self.sigma_in = (
            neutron_producing
            * nXn_weight_matrix(self.group_structure, self.endf_record.q_values.get(16, 0.0), 2).T
        ).T
        return


class DummyExtractedNuclearData(ExtractedNuclearData):
    """A dummy class that returns zeros for all of the relevant cross-sections."""
    def __init__(self, group_structure, atomic_mass):
        self.group_structure = np.asarray(group_structure)
        self.atomic_mass = atomic_mass
        n_groups = len(self.group_structure) - 1
        self.sigma_total = np.zeros(n_groups)
        self.sigma_scatter = np.zeros([n_groups, n_groups])
        self.sigma_in = np.zeros([n_groups, n_groups])
        self.sigma_triton = np.zeros(n_groups)


def read_elem_and_mass(iso_name: str) -> tuple[str, str]:
    """
    Read the element and mass number of a string representation of an isotope name.

    Parameters
    ----------
    iso_name:
        an isotope name e.g. "Fe56"

    Returns
    -------
    element_name:
        A string
    mass_number:
        An integer
    """
    length = len(iso_name)
    name, mass = [], []
    for i in range(length):
        name.append(iso_name[i])
        if iso_name[i+1].isnumeric():
            break
    for j in range(i+1, length):
        mass.append(iso_name[j])
    return "".join(name), "".join(mass)


def is_natural(iso_name: str):
    element, mass = read_elem_and_mass(iso_name)
    return mass=="0"


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

    def _add_data_from_single_record(
        self, xs_data: ExtractedNuclearData, partial_number_density: float
    ) -> None:
        if not np.isclose(
            self.group_structure, xs_data.group_structure, atol=0.0
        ).all():
            raise ValueError(f"Mismatched group structure with {xs_data}.")
        self._sigma_total += xs_data.sigma_total * partial_number_density * BARNS_TO_M2
        self._sigma_scatter += xs_data.sigma_scatter * partial_number_density * BARNS_TO_M2
        self._sigma_in += xs_data.sigma_in * partial_number_density * BARNS_TO_M2
        self._sigma_triton += xs_data.sigma_triton * partial_number_density * BARNS_TO_M2

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
                for isotope, natural_abundance in get_isotopic_composition(element).items():
                    self._add_data_from_single_record(
                        xs_dict[isotope], natural_abundance * elem_frac * self.number_density
                    )
        self._confirm_sigma()
        return

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


def accumulate_data_requirements(mat_list: Iterable[MaterialMacroInfo]) -> tuple[set[str], set[str]]:
    """
    Get the set of isotopes and elements that needs to be extracted from the
    nuclear data library.

    Parameters
    ----------
    mat_list:
        an Iterable of MaterialMacroInfo
    
    Returns
    -------
    required_element_set:
        a set containing all of the elements that has appeared on the mat_list
    required_isotope_set:
        a set containing all of the isotopes that has appeared on the mat_list
    """
    required_element_set, required_isotope_set = set(), set()
    for mat in mat_list:
        required_element_set = required_element_set.union(mat.element_set)
    for elem in required_element_set:
        for isotope in get_isotopic_composition(elem):
            required_isotope_set.add(isotope)
    return required_element_set, required_isotope_set

def populate_from_nuclear_data_library(
    endf_record_path_generator: Iterable[Path],
    mat_list: Iterable[MaterialMacroInfo],
    quick_isotope_checker: Callable[[Path], tuple[str, str, bool]],
    nuclear_data_extractor: Callable[[Path, npt.NDArray[np.float64]], ExtractedNuclearData],
    group_structure: npt.NDArray[np.float64],
    *,
    silence_missing_error: bool=False,
) -> None:
    """
    Read only the relevant files from endf_record_path_generator, and then
    save these data as ExtractedNuclearData. These data is then used to mutate
    the MaterialMacroInfo in mat_list, so that they're now populated with the
    required data.

    Parameters
    ----------
    endf_record_path_generator:
        a function that generates the list of paths that will be checked for
        matching nuclear data
    mat_list:
        an iterable of MaterialMacroInfo from which the required isotope set
        will be observed.
    quick_isotope_checker:
        A function to quickly get the isotope name using only the endf file 
        name. This should preferably not involve parsing the entire endf file.
    nuclear_data_extractor:
        A function to make the ExtractedNuclearData when given an endf_record
        Path.

    Returns
    -------
    :
    """
    # scrape for all isotopes required.
    required_element_set, required_isotope_set = accumulate_data_requirements(
        mat_list
    )

    # Gather the records
    data_extracted = {}
    for endf_record in endf_record_path_generator:
        isotope = quick_isotope_checker(endf_record)
        element, mass = read_elem_and_mass(isotope)
        if is_natural(isotope) and element in required_element_set:
            required_element_set.remove(read_elem_and_mass(isotope)[0])
            data_extracted[isotope] = nuclear_data_extractor(endf_record, group_structure)
        elif isotope in required_isotope_set:
            required_isotope_set.remove(isotope)
            data_extracted[isotope] = nuclear_data_extractor(endf_record, group_structure)

    # Check no data is missing
    for isotope in list(required_isotope_set):
        if read_elem_and_mass(isotope)[0] not in required_element_set:
            required_isotope_set.remove(isotope)
    if required_isotope_set:
        if not silence_missing_error:
            raise RuntimeError(
                f"Missing nuclear data records for {required_isotope_set}"
            )
        # fill with dummy
        for isotope in required_isotope_set:
            data_extracted[isotope] = DummyExtractedNuclearData(
                group_structure, read_elem_and_mass(isotope)[1]
            )

    for mat in mat_list:
        mat.populate_from_data_library(data_extracted)
    return