"""
Most of derivation is basd on the book Reactor Analysis, Duderstadt and Hamilton, 1976,
ISBN:9780471223634.
"""

import functools

import numpy as np
from neutronics_data import MaterialMacroInfo
from numpy import typing as npt
from scipy.special import expm1

from process.exceptions import ProcessValidationError


def groupwise(func):
    """Rename the current func as groupwise_func, and over-write the current method
    as one that sums up all group-wise func.
    """
    method_name = func.__name__
    groupwise_name = f"groupwise_{method_name}"

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        groupwise_func = getattr(self, groupwise_name)
        return np.array(
            [
                groupwise_func(n, *args, **kwargs)
                for n in range(1, self.n_groups + 1)
            ],
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
) -> tuple[float, float]:
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
    diffusion_len_2:
        The square of the characteristic diffusion length as given by Reactor Analysis,
        Duderstadt and Hamilton.
        unit: [cm]
    """

    transport_xs = total_xs - 2 / (3 * avg_atomic_mass) * scattering_xs
    diffusion_coef = 1 / 3 / transport_xs
    diffusion_len_2 = diffusion_coef / (total_xs - scattering_xs)
    return diffusion_coef, diffusion_len_2


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


class NeutronFluxProfile:
    """Neutron flux in the first wall or the blanket."""

    def __init__(
        self,
        flux: float,
        x_fw: float,
        x_bz: float,
        fw_mat: MaterialMacroInfo,
        bz_mat: MaterialMacroInfo,
    ):
        """Initialize a particular FW-BZ geometry and neutron flux.

        Parameters
        ----------
        flux:
            Nuetron flux directly emitted by the plasma, incident on the first wall.
            unit: cm^-2 s^-1
        x_fw:
            thickness of the first wall [m]. It will be converted and stored in [cm].
        x_bz:
            thickness of the blanket [m]. It will be converted and stored in [cm].
        fw_mat:
            first wall material information
        bz_mat:
            blanket material information

        Attributes
        ----------
        x_fw:
            thickness of the first wall, converted from [m] into [cm].
        x_bz:
            thickness of the blanket, converted from [m] into [cm].
        n_groups:
            number of groups in the group structure
        group_structure:
            energy bin edges, 1D array of len = n_groups+1
        """
        self.flux = flux  # flux incident on the first wall.
        if not (0 < x_fw < x_bz):
            raise ValueError(
                f"Cannot construct a first-wall+blanket module where{x_fw=}, {x_bz=}."
            )
        self.x_fw = x_fw * 100
        self.x_bz = x_bz * 100
        self.fw_mat = fw_mat
        self.bz_mat = bz_mat
        self.n_groups = self.fw_mat.n_groups
        self.group_structure = self.fw_mat.group_structure
        if not np.isclose(
            self.fw_mat.group_structure, self.bz_mat.group_structure
        ):
            raise ProcessValidationError(
                "The first-wall material info and breeding zone material info"
                "must have the same group structure!"
            )

        # macroscopic cross-sections Sigma saved here.
        # Dictionaries indexed by integers, so that we can create these values
        # out of sequence.
        self.integration_constants = {}
        self.l_fw_2, self.l_bz_2 = {}, {}  # diffusion lengths squared
        self.extended_boundary = {}

    def solve_lowest_group(self) -> None:
        """
        Solve the highest-energy (lowest-lethargy)-group's neutron diffusion equation.
        Store the solved constants in self.extended_boundary[0], self.l_fw_2[0],
        self.l_bz_2[0], and self.integration_constants[0].
        """
        i = 0
        if i in self.integration_constants:
            return  # skip if it has already been solved.
        self.l_fw_2[i], d_fw = get_diffusion_coefficient_and_length(
            self.fw_mat.sigma_t[i],
            self.fw_mat.sigma_s[i, i],
            self.fw_mat.avg_atomic_mass,
        )
        self.l_bz_2[i], d_bz = get_diffusion_coefficient_and_length(
            self.bz_mat.sigma_t[i],
            self.bz_mat.sigma_s[i, i],
            self.bz_mat.avg_atomic_mass,
        )
        l_fw = np.sqrt(abs(self.l_fw_2[i]))
        l_bz = np.sqrt(abs(self.l_bz_2[i]))
        x_fw, x_bz = self.x_fw, self.x_bz
        self.extended_boundary[i] = x_bz + extrapolation_length(d_bz)
        if self.l_fw_2[i] > 0:
            sinh_fw = np.sinh(x_fw / l_fw)
            cosh_fw = np.cosh(x_fw / l_fw)
            tanh_fw = np.tanh(x_fw / l_fw)
        else:
            sinh_fw = np.sin(x_fw / l_fw)
            cosh_fw = np.cos(x_fw / l_fw)
            tanh_fw = np.tan(x_fw / l_fw)
        if self.l_bz_2[i] > 0:
            sinh_bz = np.sinh((self.extended_boundary[i] - x_fw) / l_bz)
            cosh_bz = np.cosh((self.extended_boundary[i] - x_fw) / l_bz)
            tanh_bz = np.tanh((self.extended_boundary[i] - x_fw) / l_bz)
        else:
            sinh_bz = np.sin((self.extended_boundary[i] - x_fw) / l_bz)
            cosh_bz = np.cos((self.extended_boundary[i] - x_fw) / l_bz)
            tanh_bz = np.tan((self.extended_boundary[i] - x_fw) / l_bz)

        c2 = (
            self.flux
            * np.exp(x_fw / l_fw)
            / 2
            * ((l_fw / d_fw) + (l_bz / d_bz) * tanh_bz)
            / (cosh_fw + sinh_fw * tanh_bz * (d_fw / l_fw) * (l_bz / d_bz))
        )
        c1 = c2 - l_fw / d_fw * self.flux

        c3_c4_common_factor = (
            self.flux
            * np.exp(x_fw / l_fw)
            / 2
            * (1 - tanh_fw)
            / ((d_bz / l_bz) * cosh_bz + (d_fw / l_fw) * tanh_fw * sinh_bz)
        )
        c3 = c3_c4_common_factor * np.exp(self.extended_boundary[i] / l_bz)
        c4 = -c3_c4_common_factor * np.exp(-self.extended_boundary[i] / l_bz)
        self.integration_constants[i] = [c1, c2, c3, c4]

    def solve_group_n(self, n: int) -> None:
        """
        Solve the n-th group of neutron's diffusion equation.
        Store the solved constants in self.extended_boundary[n-1], self.l_fw_2[n-1],
        self.l_bz_2[n-1], and self.integration_constants[n-1].

        Parameters
        ----------
        n:
            The index of the neutron group whose constants are being solved.
            allowed range: [1, self.n_groups]
            This gets translated to the index i=n-1 used for the dictionaries of
            constants attached to self.
        """
        if n not in range(1, self.n_groups):
            raise ValueError(
                f"n must be a positive integer between 1 and {self.n_groups}!"
            )
        if n == 1:
            return self.solve_lowest_group()
        for k in range(n - 1):
            if k not in self.integration_constants:
                self.solve_group_n(k)
        i = n - 1
        if i in self.integration_constants:
            return None  # skip if it has already been solved.
        self.l_fw_2[i], d_fw = get_diffusion_coefficient_and_length(
            self.fw_mat.sigma_t[i],
            self.fw_mat.sigma_s[i, i],
            self.fw_mat.avg_atomic_mass,
        )
        self.l_bz_2[i], d_bz = get_diffusion_coefficient_and_length(
            self.bz_mat.sigma_t[i],
            self.bz_mat.sigma_s[i, i],
            self.bz_mat.avg_atomic_mass,
        )
        l_fw = np.sqrt(abs(self.l_fw_2[i]))
        l_bz = np.sqrt(abs(self.l_bz_2[i]))
        c1 = ...
        c2 = ...
        c3 = ...
        c4 = ...
        self.extended_boundary[i] = self.x_bz + extrapolation_length(d_bz)
        self.integration_constants[i] = [c1, c2, c3, c4]
        return None

    @groupwise
    def neutron_flux_fw(self, n: int, x: float | npt.NDArray) -> npt.NDArray:
        """Neutron flux of the n-th group at the first wall, at location x [m].

        Parameters
        ----------
        n:
            The index of the neutron group whose flux is being evaluated.
        x:
            The position where the neutron flux has to be evaluated.
            Note that this function does not enforce a check for x=inside the first-wall,
            thus if x is outside the first-wall, an extrapolated first-wall flux value up
            to that point will be given, this flux is still guaranteed to be non-singular
            i.e. finite, but not guaranteed to be positive.

        Returns
        -------
        flux:
            Neutron flux at x meter from the first wall.
        """
        i = n - 1
        c1, c2 = self.integration_constants[i][:2]
        l_fw = np.sqrt(abs(self.l_fw_2[i]))
        x_l_fw = abs(x) / l_fw
        return c1 * np.exp(x_l_fw) + c2 * np.exp(-x_l_fw)

    @groupwise
    def neutron_flux_bz(self, n: int, x: float | npt.NDArray) -> npt.NDArray:
        """Neutron flux of the n-th groupat the blanket, at location x [m].

        Parameters
        ----------
        n:
            The index of the neutron group whose flux is being evaluated.
        x:
            The position where the neutron flux has to be evaluated. [m]
            Note that this function does not enforce a check for x=inside the blanket,
            thus if x is outside the blanket, an extrapolated blanket flux value up to
            that point will be given, this flux is still guaranteed to be non-singular
            i.e. finite, but not guaranteed to be positive.

        Returns
        -------
        flux:
            Neutron flux at x meter from the first wall.
        """
        i = n - 1
        c3, c4 = self.integration_constants[i][2:]
        l_bz = np.sqrt(abs(self.l_bz_2[i]))
        x_l_bz = abs(x) / l_bz
        return c3 * np.exp(x_l_bz) + c4 * np.exp(-x_l_bz)

    @groupwise
    def neutron_flux_at(self, n: int, x: float | npt.NDArray) -> npt.NDArray:
        """
        Neutron flux anywhere within the valid range of x,
        i.e. from -self.x_bz to self.x_bz.

        Parameters
        ----------
        n:
            Neutron group number
        x:
            The depth where we want the neutron flux (m).
            Valid only between x= -extended boundary to +-extended boundary of that group

        Raises
        ------
        ValueError
            The inputted x
        """
        if np.isscalar(x):
            return self.groupwise_neutron_flux_at(n, [x])[0]
        x = np.asarray(x)
        in_fw = abs(x) <= self.x_fw
        in_bz = np.logical_and(self.x_fw < abs(x), abs(x) <= self.x_bz)
        if (~np.logical_or(in_fw, in_bz)).any():
            raise ValueError(
                f"for neutron group {n}, neutron flux can only be calculated up to "
                f"{self.extended_boundary[n]} cm, which {x * 100} cm violates!"
            )

        out_flux = np.zeros_like(x)
        out_flux[in_fw] = self.groupwise_neutron_flux_fw(n, x[in_fw])
        out_flux[in_bz] = self.groupwise_neutron_flux_bz(n, x[in_bz])
        return out_flux

    # scalar values (one or two such floats per neutron group.)
    @groupwise
    def reaction_rate_fw(self, n: int, reaction_type: str) -> float:
        """Calculate the reaction rate in the first wall.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved.
        reaction_type:
            Two options: "total" or "non-scattering".

        """
        self.solve_group_n(n)
        i = n - 1
        l_fw = np.sqrt(abs(self.l_fw_2[i]))
        c1, c2, c3, c4 = self.integration_constants[i]
        if reaction_type == "non-scattering":
            sigma = self.fw_mat.sigma_t[i] - self.fw_mat.sigma_s[i, i]
        elif reaction_type == "total":
            sigma = self.fw_mat.sigma_t[i]
        else:
            raise NotImplementedError(
                f"Not yet implemented the reaction type {reaction_type}"
            )
        return sigma * (
            c1 * expm1(self.x_fw / l_fw) - c2 * expm1(-self.x_fw / l_fw)
        )

    @groupwise
    def reaction_rate_bz(self, n: int, reaction_type: str) -> float:
        """Calculate the reaction rate in the blanket.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved.
        reaction_type:
            Two options: "total" or "non-scattering".

        """
        self.solve_group_n(n)
        i = n - 1
        l_bz = np.sqrt(abs(self.l_bz_2[i]))
        c1, c2, c3, c4 = self.integration_constants[i]
        if reaction_type == "non-scattering":
            sigma = self.bz_mat.sigma_t[i] - self.bz_mat.sigma_s[i, i]
        elif reaction_type == "total":
            sigma = self.bz_mat.sigma_t[i]
        else:
            raise NotImplementedError(
                f"Not yet implemented the reaction type {reaction_type}"
            )
        # thicknesses in terms of bz path lengths
        bz_thick = (self.x_bz - self.x_fw) / l_bz
        fw_thick = self.x_fw / l_bz
        return sigma * (
            c3 * np.exp(fw_thick) * expm1(bz_thick)
            - c4 * np.exp(-fw_thick) * expm1(-bz_thick)
        )

    @groupwise
    def flux_fw2bz(self, n: int) -> float:
        self.solve_group_n(n)
        i = n - 1
        c1, c2, c3, c4 = self.integration_constants[i]
        l_bz_2, d_bz = get_diffusion_coefficient_and_length(
            self.bz_mat.sigma_t[i],
            self.bz_mat.sigma_s[i, i],
            self.bz_mat.avg_atomic_mass,
        )
        l_bz = np.sqrt(abs(l_bz_2))
        return (
            -d_bz
            / l_bz
            * (c3 * np.exp(self.x_fw / l_bz) - c4 * np.exp(-self.x_fw / l_bz))
        )
        # equivalent definition below: (should yield the same answer)
        # self.flux_fw2bz[i] = - d_fw / l_fw * (c1 * np.exp(x_fw/l_fw) - c2 * np.exp(-x_fw/l_fw))

    @groupwise
    def flux_escaped(self, n: int) -> float:
        self.solve_group_n(n)
        i = n - 1
        c1, c2, c3, c4 = self.integration_constants[i]
        l_bz_2, d_bz = get_diffusion_coefficient_and_length(
            self.bz_mat.sigma_t[i],
            self.bz_mat.sigma_s[i, i],
            self.bz_mat.avg_atomic_mass,
        )
        l_bz = np.sqrt(abs(l_bz_2))
        self.flux_escaped[i] = (
            -d_bz
            / l_bz
            * (c3 * np.exp(self.x_bz / l_bz) - c4 * np.exp(-self.x_bz / l_bz))
        )
