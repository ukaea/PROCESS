"""
Most of derivation is basd on the book Reactor Analysis, Duderstadt and Hamilton, 1976,
ISBN:9780471223634.
"""

import functools
import inspect

import numpy as np
from numpy import typing as npt
from scipy.special import expm1

from process.exceptions import ProcessValidationError
from process.neutronics_data import MaterialMacroInfo


def summarize_values(func):
    """
    Keep groupwise_func unchanged, but create a new method under a similar name
    (but with the prefix "groupwise_" removed) which outputs the sum of every
    groupwise value.
    """
    summary_method_name = func.__name__[10:]
    # confirm this is a groupwise method
    func_params = inspect.signature(func).parameters
    if not (func.__name__.startswith("groupwise_") and "n" in func_params):
        raise ValueError(
            "The decorated method is designed to turn groupwise methods into "
            "summed flux/reaction/current methods."
        )

    @functools.wraps(func)
    def wrapper(self, *args, **kwargs):
        groupwise_func = getattr(self, func.__name__)
        return np.sum(
            [groupwise_func(n, *args, **kwargs) for n in range(self.n_groups)],
            axis=0,
        )

    def wrapper_setattr(cls):
        """
        Attach the method that outputs the summed values of all of the
        groupwise values to the same parent class.
        """
        setattr(cls, func.__name__, func)
        setattr(cls, summary_method_name, wrapper)
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
            thickness of the blanket + first-wall [m]. It will be converted and stored in [cm].
        fw_mat:
            first wall material information
        bz_mat:
            blanket material information

        Attributes
        ----------
        x_fw_cm:
            thickness of the first wall, converted from [m] into [cm].
        x_bz_cm:
            thickness of the blanket + first-wall, converted from [m] into [cm].
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
        self.x_fw_cm = x_fw * 100
        self.x_bz_cm = x_bz * 100
        self.fw_mat = fw_mat
        self.bz_mat = bz_mat
        self.n_groups = self.fw_mat.n_groups
        self.group_structure = self.fw_mat.group_structure
        if not np.allclose(
            self.fw_mat.group_structure,
            self.bz_mat.group_structure,
            atol=0,
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
        n = 0
        if n in self.integration_constants:
            return  # skip if it has already been solved.
        self.l_fw_2[n], d_fw = get_diffusion_coefficient_and_length(
            self.fw_mat.sigma_t[n],
            self.fw_mat.sigma_s[n, n],
            self.fw_mat.avg_atomic_mass,
        )
        self.l_bz_2[n], d_bz = get_diffusion_coefficient_and_length(
            self.bz_mat.sigma_t[n],
            self.bz_mat.sigma_s[n, n],
            self.bz_mat.avg_atomic_mass,
        )
        l_fw = np.sqrt(abs(self.l_fw_2[n]))
        l_bz = np.sqrt(abs(self.l_bz_2[n]))
        x_fw, x_bz = self.x_fw_cm, self.x_bz_cm
        self.extended_boundary[n] = x_bz + extrapolation_length(d_bz)
        if self.l_fw_2[n] > 0:
            sinh_fw = np.sinh(x_fw / l_fw)
            cosh_fw = np.cosh(x_fw / l_fw)
            tanh_fw = np.tanh(x_fw / l_fw)
        else:
            sinh_fw = np.sin(x_fw / l_fw)
            cosh_fw = np.cos(x_fw / l_fw)
            tanh_fw = np.tan(x_fw / l_fw)
        if self.l_bz_2[n] > 0:
            sinh_bz = np.sinh((self.extended_boundary[n] - x_fw) / l_bz)
            cosh_bz = np.cosh((self.extended_boundary[n] - x_fw) / l_bz)
            tanh_bz = np.tanh((self.extended_boundary[n] - x_fw) / l_bz)
        else:
            sinh_bz = np.sin((self.extended_boundary[n] - x_fw) / l_bz)
            cosh_bz = np.cos((self.extended_boundary[n] - x_fw) / l_bz)
            tanh_bz = np.tan((self.extended_boundary[n] - x_fw) / l_bz)

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
        c3 = c3_c4_common_factor * np.exp(self.extended_boundary[n] / l_bz)
        c4 = -c3_c4_common_factor * np.exp(-self.extended_boundary[n] / l_bz)
        self.integration_constants[n] = [c1, c2, c3, c4]

    def solve_group_n(self, n: int) -> None:
        """
        Solve the n-th group of neutron's diffusion equation, where n<=n_groups-1.
        Store the solved constants in self.extended_boundary[n-1], self.l_fw_2[n-1],
        self.l_bz_2[n-1], and self.integration_constants[n-1].

        Parameters
        ----------
        n:
            The index of the neutron group whose constants are being solved. n <= n_groups - 1.
            Therefore n=0 shows the reaction rate for group 1, n=1 for group 2, etc.
        """
        if n not in range(self.n_groups):
            raise ValueError(
                f"n must be a positive integer between 1 and {self.n_groups}!"
            )
        if n == 0:
            return self.solve_lowest_group()
        for k in range(n):
            if k not in self.integration_constants:
                self.solve_group_n(k)
        if n in self.integration_constants:
            return None  # skip if it has already been solved.
        self.l_fw_2[n], d_fw = get_diffusion_coefficient_and_length(
            self.fw_mat.sigma_t[n],
            self.fw_mat.sigma_s[n, n],
            self.fw_mat.avg_atomic_mass,
        )
        self.l_bz_2[n], d_bz = get_diffusion_coefficient_and_length(
            self.bz_mat.sigma_t[n],
            self.bz_mat.sigma_s[n, n],
            self.bz_mat.avg_atomic_mass,
        )
        l_fw = np.sqrt(abs(self.l_fw_2[n]))
        l_bz = np.sqrt(abs(self.l_bz_2[n]))
        c1 = ...
        c2 = ...
        c3 = ...
        c4 = ...
        self.extended_boundary[n] = self.x_bz_cm + extrapolation_length(d_bz)
        self.integration_constants[n] = [c1, c2, c3, c4]
        return None

    @summarize_values
    def groupwise_neutron_flux_fw(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """Neutron flux of the n-th group at the first wall, at location x [m].

        Parameters
        ----------
        n:
            The index of the neutron group whose flux is being evaluated. n <= n_groups - 1.
            Therefore n=0 shows the reaction rate for group 1, n=1 for group 2, etc.
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
        c1, c2 = self.integration_constants[n][:2]
        l_fw = np.sqrt(abs(self.l_fw_2[n]))
        x_l_fw = abs(x) / l_fw
        return c1 * np.exp(x_l_fw) + c2 * np.exp(-x_l_fw)

    @summarize_values
    def groupwise_neutron_flux_bz(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """Neutron flux of the n-th groupat the blanket, at location x [m].

        Parameters
        ----------
        n:
            The index of the neutron group whose flux is being evaluated. n <= n_groups - 1.
            Therefore n=0 shows the reaction rate for group 1, n=1 for group 2, etc.
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
        c3, c4 = self.integration_constants[n][2:]
        l_bz = np.sqrt(abs(self.l_bz_2[n]))
        x_l_bz = abs(x) / l_bz
        return c3 * np.exp(x_l_bz) + c4 * np.exp(-x_l_bz)

    @summarize_values
    def groupwise_neutron_flux_at(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """
        Neutron flux anywhere within the valid range of x,
        i.e. from -self.x_bz_cm to self.x_bz_cm.

        Parameters
        ----------
        n:
            Neutron group index. n <= n_groups - 1.
            Therefore n=0 shows the reaction rate for group 1, n=1 for group 2, etc.
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
        in_fw = abs(x) <= self.x_fw_cm
        in_bz = np.logical_and(self.x_fw_cm < abs(x), abs(x) <= self.x_bz_cm)
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
    @summarize_values
    def groupwise_reaction_rate_fw(self, n: int, reaction_type: str) -> float:
        """Calculate the reaction rate in the first wall.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved. n <= n_groups - 1.
            Therefore n=0 shows the reaction rate for group 1, n=1 for group 2, etc.
        reaction_type:
            Two options: "total" or "non-scattering".

        """
        self.solve_group_n(n)
        l_fw = np.sqrt(abs(self.l_fw_2[n]))
        c1, c2, c3, c4 = self.integration_constants[n]
        if reaction_type == "non-scattering":
            sigma = self.fw_mat.sigma_t[n] - self.fw_mat.sigma_s[n, n]
        elif reaction_type == "total":
            sigma = self.fw_mat.sigma_t[n]
        else:
            raise NotImplementedError(
                f"Not yet implemented the reaction type {reaction_type}"
            )
        return sigma * (
            c1 * expm1(self.x_fw_cm / l_fw) - c2 * expm1(-self.x_fw_cm / l_fw)
        )

    @summarize_values
    def groupwise_reaction_rate_bz(self, n: int, reaction_type: str) -> float:
        """Calculate the reaction rate in the blanket.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved. n <= n_groups - 1.
            Therefore n=0 shows the reaction rate for group 1, n=1 for group 2, etc.
        reaction_type:
            Two options: "total" or "non-scattering".

        """
        self.solve_group_n(n)
        l_bz = np.sqrt(abs(self.l_bz_2[n]))
        c1, c2, c3, c4 = self.integration_constants[n]
        if reaction_type == "non-scattering":
            sigma = self.bz_mat.sigma_t[n] - self.bz_mat.sigma_s[n, n]
        elif reaction_type == "total":
            sigma = self.bz_mat.sigma_t[n]
        else:
            raise NotImplementedError(
                f"Not yet implemented the reaction type {reaction_type}"
            )
        # thicknesses in terms of bz path lengths
        bz_thick = (self.x_bz_cm - self.x_fw_cm) / l_bz
        fw_thick = self.x_fw_cm / l_bz
        return sigma * (
            c3 * np.exp(fw_thick) * expm1(bz_thick)
            - c4 * np.exp(-fw_thick) * expm1(-bz_thick)
        )

    @summarize_values
    def groupwise_current_fw2bz(self, n: int) -> float:
        """
        Net current from the first wall to breeding zone.
        Parameters
        ----------
        n:
            The index of the neutron group that we want the current for. n <= n_groups - 1.
            Therefore n=0 shows the reaction rate for group 1, n=1 for group 2, etc.

        Returns
        -------
        :
            current in cm^-2
        """
        self.solve_group_n(n)
        c1, c2, c3, c4 = self.integration_constants[n]
        l_bz_2, d_bz = get_diffusion_coefficient_and_length(
            self.bz_mat.sigma_t[n],
            self.bz_mat.sigma_s[n, n],
            self.bz_mat.avg_atomic_mass,
        )
        l_bz = np.sqrt(abs(l_bz_2))
        return (
            -d_bz
            / l_bz
            * (
                c3 * np.exp(self.x_fw_cm / l_bz)
                - c4 * np.exp(-self.x_fw_cm / l_bz)
            )
        )
        # equivalent definition below: (should yield the same answer)
        # return - d_fw / l_fw * (c1 * np.exp(x_fw/l_fw) - c2 * np.exp(-x_fw/l_fw))

    @summarize_values
    def groupwise_current_escaped(self, n: int) -> float:
        """
        Neutron current escaped from the breeding zone to outside the reactor.
        Parameters
        ----------
        n:
            The index of the neutron group that we want the current for. n <= n_groups - 1.
            Therefore n=0 shows the reaction rate for group 1, n=1 for group 2, etc.

        Returns
        -------
        :
            current in cm^-2
        """
        self.solve_group_n(n)
        c1, c2, c3, c4 = self.integration_constants[n]
        l_bz_2, d_bz = get_diffusion_coefficient_and_length(
            self.bz_mat.sigma_t[n],
            self.bz_mat.sigma_s[n, n],
            self.bz_mat.avg_atomic_mass,
        )
        l_bz = np.sqrt(abs(l_bz_2))
        return (
            -d_bz
            / l_bz
            * (
                c3 * np.exp(self.x_bz_cm / l_bz)
                - c4 * np.exp(-self.x_bz_cm / l_bz)
            )
        )
