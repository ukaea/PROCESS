"""
Most of derivation is basd on the book Reactor Analysis, Duderstadt and Hamilton, 1976,
ISBN:9780471223634.
"""

import functools
import inspect
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass

import numpy as np
from matplotlib import pyplot as plt
from numpy import typing as npt
from scipy import optimize

from process.exceptions import ProcessValidationError, ProcessValueError
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
    avg_atomic_mass: float,
    total_xs_cm: float,
    scattering_xs_cm: float,
    in_source_xs_cm: float,
) -> tuple[float, float]:
    r"""
    Calculate the diffusion coefficient for a given scattering and total macro-scopic
    cross-section in a given medium.

    Parameters
    ----------
    total_xs_cm:
        macroscopic total cross-section `\sigma_{total}`, i.e. any reaction between
        nuclei and neutrons, that either changes the neutron's path or remove it from
        that energy group.
        Unit: cm^-1
    scattering_xs_cm:
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

    transport_xs = total_xs_cm - 2 / (3 * avg_atomic_mass) * scattering_xs_cm
    diffusion_coef = 1 / 3 / transport_xs
    diffusion_len_2 = diffusion_coef / (
        total_xs_cm - scattering_xs_cm - in_source_xs_cm
    )
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


@dataclass
class IntegrationConstants:
    """
    Inside each material, there are two new hyperbolic trig funcs per group, i.e.
    group n=0 has cosh(x/L[0]) and sinh(x/L[0]),
    group n=1 has cosh(x/L[0]), sinh(x/L[0]), cosh(x/L[1]) and sinh(x/L[1]),
    etc.
    To get the neutron flux, each trig func has to be scaled by an integration
    constant. E.g.
    group n=0: fw: fw_c[0][0] * cosh(x/L[0]) + fw_s[0][0] * sinh(x/L[0])
    group n=1: fw: fw_c[1][0] * cosh(x/L[0]) + fw_c[1][1] * cosh(x/L[1])
                + fw_s[1][0] * sinh(x/L[0]) + fw_s[1][1] * sinh(x/L[1])

    """

    fw_c: Iterable[float]
    fw_s: Iterable[float]
    bz_c: Iterable[float]
    bz_s: Iterable[float]

    def validate_length(self, expected_length: int):
        """Validate that all fields has the correct length."""
        for const_name, const_value in asdict(self).items():
            if len(const_value) != expected_length:
                raise ProcessValueError(
                    f"Expected {const_name} to have len=={expected_length}, "
                    f"got {len(const_value)} instead."
                )


class AutoPopulatingDict:
    """
    Class that behaves like a dictionary, but if the required key does not
    exist in the dictionary, it will call the populating_method to populate
    that specific key.
    """

    def __init__(self, populating_method: Callable[[int], None], name: str):
        """
        Attributes
        ----------
        _dict:
            A dictionary indexed by integers, so that we can populate its
            values out of sequence.
        populating_method:
            The method to be called if the requested index n is currently
            unpopulated in the dictionary. This method should populate the
            dictionary.
        """
        self._dict = {}
        self._name = name
        self._attempting_to_access = set()
        self.populating_method = populating_method

    def __getitem__(self, i: int):
        """Check if index i is in the dictionary or not. If not, populate it."""
        if i in self._attempting_to_access:
            raise RecursionError(
                f"retrieving the value of {self._name}[{i}] requires the "
                f"value of {self._name}[{i}]."
            )
        if i not in self._dict:
            self._attempting_to_access.add(i)
            self.populating_method(i)
            if i not in self._dict:
                raise RuntimeError(
                    f"{self.populating_method}({i}) failed to populate key {i} "
                    "in the dictionary!"
                )
            self._attempting_to_access.discard(i)
        return self._dict[i]

    def __contains__(self, i: int):
        """Check if key 'i' is in the dictionary or not."""
        return i in self._dict

    def __setitem__(self, i: int, value: float):
        """Check if dict i is in the index or not."""
        if hasattr(value, "validate_length"):
            value.validate_length(i + 1)
        self._dict[i] = value
        self._attempting_to_access.discard(i)

    def values(self):
        return self._dict.values()

    def __repr__(self):
        return f"AutoPopulatingDict{self._dict}"


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
        fw_mat:
            first wall material information
        bz_mat:
            blanket material information
        n_groups:
            number of groups in the group structure
        group_structure:
            energy bin edges, 1D array of len = n_groups+1

        integration_constants:
            Integration constants that determine the flux shape (and therefore
            reaction rates and neutron current) of each group. A set of four
            constants, two for fw and two for bz; each with unit: [cm^-2 s^-1]
        l_fw_2:
            square of the characteristic diffusion length as given by Reactor Analysis,
            Duderstadt and Hamilton, applied to the fw. unit: [cm]
        l_bz_2:
            square of the characteristic diffusion length as given by Reactor Analysis,
            Duderstadt and Hamilton, applied to the bz. unit: [cm]
        d_fw_cm:
            diffusion constants in the fw. unit: [cm]
        d_bz_cm:
            diffusion constants in the bz. unit: [cm]
        extended_boundary_cm:
            extended boundary (outside the bz) for each group, stored in cm.
            This value should be larger than x_bz for each of them.
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

        self.integration_constants = AutoPopulatingDict(
            self.solve_group_n, "integration_constants"
        )
        # diffusion lengths squared
        self.l_fw_2 = AutoPopulatingDict(self.solve_group_n, "l_fw_2")
        self.l_bz_2 = AutoPopulatingDict(self.solve_group_n, "l_bz_2")
        self.d_fw_cm = AutoPopulatingDict(self.solve_group_n, "d_fw_cm")
        self.d_bz_cm = AutoPopulatingDict(self.solve_group_n, "d_bz_cm")
        self.extended_boundary_cm = AutoPopulatingDict(
            self.solve_group_n, "extended_boundary_cm"
        )

    def solve_lowest_group(self) -> None:
        """
        Solve the highest-energy (lowest-lethargy)-group's neutron diffusion equation.
        Store the solved constants in self.extended_boundary_cm[0], self.l_fw_2[0],
        self.l_bz_2[0], and self.integration_constants[0].
        integration_constants have units of [cm^-2 s^-1].
        """
        n = 0
        if n in self.integration_constants:
            return  # skip if it has already been solved.
        self.d_fw_cm[n], self.l_fw_2[n] = get_diffusion_coefficient_and_length(
            self.fw_mat.avg_atomic_mass,
            self.fw_mat.sigma_t_cm[n],
            self.fw_mat.sigma_s_cm[n, n],
            self.fw_mat.sigma_in_cm[n, n],
        )
        self.d_bz_cm[n], self.l_bz_2[n] = get_diffusion_coefficient_and_length(
            self.bz_mat.avg_atomic_mass,
            self.bz_mat.sigma_t_cm[n],
            self.bz_mat.sigma_s_cm[n, n],
            self.bz_mat.sigma_in_cm[n, n],
        )
        l_fw = np.sqrt(abs(self.l_fw_2[n]))
        l_bz = np.sqrt(abs(self.l_bz_2[n]))
        x_fw, x_bz = self.x_fw_cm, self.x_bz_cm
        d_fw = self.d_fw_cm[n]
        d_bz = self.d_bz_cm[n]
        self.extended_boundary_cm[n] = x_bz + extrapolation_length(d_bz)
        if self.l_fw_2[n] > 0:
            s_fw = np.sinh(x_fw / l_fw)
            c_fw = np.cosh(x_fw / l_fw)
            t_fw = np.tanh(x_fw / l_fw)
        else:
            s_fw = np.sin(x_fw / l_fw)
            c_fw = np.cos(x_fw / l_fw)
            t_fw = np.tan(x_fw / l_fw)
        if self.l_bz_2[n] > 0:
            c_bz = np.cosh(self.extended_boundary_cm[n] / l_bz)
            s_bz = np.sinh(self.extended_boundary_cm[n] / l_bz)
            c_bz_mod = np.cosh((self.extended_boundary_cm[n] - x_fw) / l_bz)
            s_bz_mod = np.sinh((self.extended_boundary_cm[n] - x_fw) / l_bz)
            t_bz_mod = np.tanh((self.extended_boundary_cm[n] - x_fw) / l_bz)
        else:
            c_bz = np.cos(self.extended_boundary_cm[n] / l_bz)
            s_bz = np.sin(self.extended_boundary_cm[n] / l_bz)
            c_bz_mod = np.cos((self.extended_boundary_cm[n] - x_fw) / l_bz)
            s_bz_mod = np.sin((self.extended_boundary_cm[n] - x_fw) / l_bz)
            t_bz_mod = np.tan((self.extended_boundary_cm[n] - x_fw) / l_bz)

        fw_c_factor = - self.flux * (
            l_fw / d_fw 
            - np.exp(x_fw / l_fw) * (
                (l_fw / d_fw) + (l_bz / d_bz) * t_bz_mod
            )/(
                c_fw + s_fw * t_bz_mod * (d_fw / l_fw) * (l_bz / d_bz)
            )
        )
        fw_s_factor = - self.flux * l_fw / d_fw

        bz_common_factor = (
            self.flux
            * np.exp(x_fw / l_fw)
            * (1 - t_fw)
            / ((d_bz / l_bz) * c_bz_mod + (d_fw / l_fw) * t_fw * s_bz_mod)
        )
        bz_c_factor = bz_common_factor * s_bz
        bz_s_factor = -bz_common_factor * c_bz

        self.integration_constants[n] = IntegrationConstants(
            [fw_c_factor], [fw_s_factor], [bz_c_factor], [bz_s_factor]
        )
        return

    def solve_group_n(self, n: int) -> None:
        """
        Solve the n-th group of neutron's diffusion equation, where n<=n_groups-1.
        Store the solved constants in self.extended_boundary_cm[n-1], self.l_fw_2[n-1],
        self.l_bz_2[n-1], and self.integration_constants[n-1].

        Parameters
        ----------
        n:
            The index of the neutron group whose constants are being solved.
            The allowed range of values = [0, self.n_groups-1]. Therefore,
            n=0 shows the reaction rate for group 1, n=1 for group 2, etc.
        """
        if n not in range(self.n_groups):
            raise ValueError(
                f"n must be a positive integer between 0 and {self.n_groups}!"
            )
        if n == 0:
            return self.solve_lowest_group()
        # ensure all lower groups are solved.
        for k in range(n):
            if k not in self.integration_constants:
                self.solve_group_n(k)
        if not (self.fw_mat.downscatter_only and self.bz_mat.downscatter_only):
            raise NotImplementedError(
                "Will implement solve_group_n in a loop later..."
            )
        if n in self.integration_constants:
            return None  # skip if it has already been solved.
        else:
            # Parameter to be changed later, to allow solving non-down-scatter
            # only systems by iterating.
            first_iteration = True 
        self.d_fw_cm[n], self.l_fw_2[n] = get_diffusion_coefficient_and_length(
            self.fw_mat.avg_atomic_mass,
            self.fw_mat.sigma_t_cm[n],
            self.fw_mat.sigma_s_cm[n, n],
            self.fw_mat.sigma_in_cm[n, n],
        )
        self.d_bz_cm[n], self.l_bz_2[n] = get_diffusion_coefficient_and_length(
            self.bz_mat.avg_atomic_mass,
            self.bz_mat.sigma_t_cm[n],
            self.bz_mat.sigma_s_cm[n, n],
            self.bz_mat.sigma_in_cm[n, n],
        )
        self.extended_boundary_cm[n] = self.x_bz_cm + extrapolation_length(
            self.d_bz_cm[n]
        )

        # Setting up aliases for shorter code
        l_fw_2, l_bz_2 = self.l_fw_2[n], self.l_bz_2[n]
        d_fw, d_bz = self.d_fw_cm[n], self.d_bz_cm[n]
        src_fw = self.fw_mat.sigma_s_cm + self.fw_mat.sigma_in_cm
        src_bz = self.bz_mat.sigma_s_cm + self.bz_mat.sigma_in_cm

        self.integration_constants[n] = IntegrationConstants([], [], [], [])
        for g in range(n):
            # if the characteristic length of group [g] coincides with the
            # characteristic length of group [n], then that particular
            # cosh/sinh would be indistinguishable from group [n]'s
            # cosh/sinh anyways, therefore we can set the coefficient to 0.
            if np.isclose(diff_fw := self.l_fw_2[g] - l_fw_2, 0):
                scale_factor_fw = 0.0
            else:
                scale_factor_fw = (l_fw_2 * self.l_fw_2[g]) / d_fw / diff_fw
            if np.isclose(diff_bz := self.l_bz_2[g] - l_bz_2, 0):
                scale_factor_bz = 0.0
            else:
                scale_factor_bz = (l_bz_2 * self.l_bz_2[g]) / d_bz / diff_bz
            self.integration_constants[n].fw_c.append(
                sum(src_fw[i, n] * self.integration_constants[i].fw_c[g] for i in range(g, n))
                * scale_factor_fw
            )
            self.integration_constants[n].fw_s.append(
                sum(src_fw[i, n] * self.integration_constants[i].fw_s[g] for i in range(g, n))
                * scale_factor_fw
            )
            self.integration_constants[n].bz_c.append(
                sum(src_bz[i, n] * self.integration_constants[i].bz_c[g] for i in range(g, n))
                * scale_factor_bz
            )
            self.integration_constants[n].bz_s.append(
                sum(src_bz[i, n] * self.integration_constants[i].bz_s[g] for i in range(g, n))
                * scale_factor_bz
            )

        fw_c_factor_guess = 0.0
        fw_s_factor_guess = 0.0
        bz_c_factor_guess = 0.0
        bz_s_factor_guess = 0.0
        self.integration_constants[n].fw_c.append(fw_c_factor_guess)
        self.integration_constants[n].fw_s.append(fw_s_factor_guess)
        self.integration_constants[n].bz_c.append(bz_c_factor_guess)
        self.integration_constants[n].bz_s.append(bz_s_factor_guess)

        def _set_constants(input_vector: Iterable[float]):
            _ic_n.fw_c[n] = input_vector[0]
            _ic_n.fw_s[n] = input_vector[1]
            _ic_n.bz_c[n] = input_vector[2]
            _ic_n.bz_s[n] = input_vector[3]

        def _evaluate_fit():
            flux_continuity = self.groupwise_neutron_flux_fw(n, self.x_fw_cm/100) - self.groupwise_neutron_flux_bz(n, self.x_fw_cm/100)
            flux_at_boundary = self.groupwise_neutron_flux_bz(n, self.extended_boundary_cm[n])
            current_continuity = self.groupwise_neutron_current_fw(n, self.x_fw_cm/100) - self.groupwise_neutron_current_bz(n, self.x_fw_cm/100)
            influx_fw, influx_bz = 0.0, 0.0
            for g in range(self.n_groups):
                if g>n and not first_iteration:
                    if not self.fw_mat.downscatter_only:
                        influx_fw += (self.fw_mat.sigma_s_cm[g, n] + self.fw_mat.sigma_in_cm[g, n]) * self.groupwise_integrated_flux_fw(g)
                    if not self.bz_mat.downscatter_only:
                        influx_bz += (self.bz_mat.sigma_s_cm[g, n] + self.bz_mat.sigma_in_cm[g, n]) * self.groupwise_integrated_flux_bz(g)
                elif g<=n:
                    influx_fw += (self.fw_mat.sigma_s_cm[g, n] + self.fw_mat.sigma_in_cm[g, n]) * self.groupwise_integrated_flux_fw(g)
                    influx_bz += (self.bz_mat.sigma_s_cm[g, n] + self.bz_mat.sigma_in_cm[g, n]) * self.groupwise_integrated_flux_bz(g)
            removal_fw = self.fw_mat.sigma_t_cm[n] * self.groupwise_integrated_flux_fw(n)
            removal_bz = self.bz_mat.szsigma_t_cm[n] * self.groupwise_integrated_flux_bz(n)
            # conservation_fw = influx_fw - removal_fw - current_fw2bz
            # conservation_bz = current_fw2bz + influx_bz - removal_bz - escaped_bz
            total_neutron_conservation = (
                influx_fw
                + influx_bz
                - removal_fw
                - removal_bz
                - self.groupwise_neutron_current_escaped(_ic_n, self, n)
            )
            return np.array([
                flux_continuity,
                flux_at_boundary,
                current_continuity,
                total_neutron_conservation,
            ])

        def objective(four_integration_constants_vector):
            _set_constants(four_integration_constants_vector)
            return _evaluate_fit()

        results = optimize.root(objective, x0=[
            fw_c_factor_guess,
            fw_s_factor_guess,
            bz_c_factor_guess,
            bz_s_factor_guess,
        ])
        _set_constants(results.res)
        return None

    @summarize_values
    def groupwise_neutron_flux_fw(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """
        Neutron flux[cm^-2 s^-1] of the n-th group at the first wall, at location x [m].

        Parameters
        ----------
        n:
            The index of the neutron group whose flux is being evaluated. n <= n_groups - 1.
            Therefore n=0 shows the flux for group 1, n=1 for group 2, etc.
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
        x_cm = x * 100

        trig_funcs = []
        for g in range(n + 1):
            if self.l_fw_2[g]>0:
                c, s = np.cosh, np.sinh
                l_fw = np.sqrt(self.l_fw_2[g])
            else:
                c, s = np.cos, np.sin
                l_fw = np.sqrt(-self.l_fw_2[g])
            trig_funcs.append(self.integration_constants[n].fw_c[g] * c(abs(x_cm) / l_fw))
            trig_funcs.append(self.integration_constants[n].fw_s[g] * s(abs(x_cm) / l_fw))
        return np.sum(trig_funcs, axis=0)

    @summarize_values
    def groupwise_neutron_flux_bz(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """Neutron flux[cm^-2 s^-1] of the n-th groupat the blanket, at location x [m].

        Parameters
        ----------
        n:
            The index of the neutron group whose flux is being evaluated. n <= n_groups - 1.
            Therefore n=0 shows the flux for group 1, n=1 for group 2, etc.
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
        x_cm = x * 100

        trig_funcs = []
        for g in range(n + 1):
            if self.l_bz_2[g]>0:
                c, s = np.cosh, np.sinh
                l_bz = np.sqrt(self.l_bz_2[g])
            else:
                c, s = np.cos, np.sin
                l_bz = np.sqrt(-self.l_bz_2[g])
            trig_funcs.append(self.integration_constants[n].bz_c[g] * c(abs(x_cm) / l_bz))
            trig_funcs.append(self.integration_constants[n].bz_s[g] * s(abs(x_cm) / l_bz))
        return np.sum(trig_funcs, axis=0)

    @summarize_values
    def groupwise_neutron_flux_at(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """
        Neutron flux [cm^-2 s^-1] anywhere within the valid range of x,
        i.e. [-self.extended_boundary_cm[n], self.extended_boundary_cm[n]].

        Parameters
        ----------
        n:
            Neutron group index. n <= n_groups - 1.
            Therefore n=0 shows the flux for group 1, n=1 for group 2, etc.
        x:
            The depth where we want the neutron flux [m].
            Valid only between x= -extended boundary to +-extended boundary of that group

        Raises
        ------
        ValueError
            The inputted x
        """
        if np.isscalar(x):
            return self.groupwise_neutron_flux_at(n, [x])[0]
        x = np.asarray(x)
        in_fw = abs(x * 100) <= self.x_fw_cm
        in_bz = np.logical_and(
            self.x_fw_cm < abs(x * 100),
            abs(x * 100) <= self.extended_boundary_cm[n],
        )
        if (~np.logical_or(in_fw, in_bz)).any():
            raise ValueError(
                f"for neutron group {n}, neutron flux can only be calculated "
                f"up to {self.extended_boundary_cm[n]} cm, which {x * 100} cm violates!"
            )

        out_flux = np.zeros_like(x)
        out_flux[in_fw] = self.groupwise_neutron_flux_fw(n, x[in_fw])
        out_flux[in_bz] = self.groupwise_neutron_flux_bz(n, x[in_bz])
        return out_flux

    # scalar values (one such float per neutron group.)
    @summarize_values
    def groupwise_integrated_flux_fw(self, n: int) -> float:
        """
        Calculate the integrated flux[cm s^-1], which can be mulitplied to any
        macroscopic cross-section [cm^-1] to get the reaction rate [s^-1] in
        the first wall.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved. n <= n_groups - 1.
            Therefore n=0 shows the integrated flux for group 1, n=1 for group 2, etc.
        """
        integrals = []
        for g in range(n + 1):
            if self.l_fw_2[g]>0:
                l_fw = np.sqrt(self.l_fw_2[g])
                integrals.append(
                    l_fw * self.integration_constants[n].fw_c[g] * np.sinh(self.x_fw_cm / l_fw)
                )
                integrals.append(
                    l_fw * self.integration_constants[n].fw_s[g] * (np.cosh(self.x_fw_cm / l_fw) - 1)
                )
            else:
                l_fw = np.sqrt(-self.l_fw_2[g])
                integrals.append(
                    l_fw * self.integration_constants[n].fw_c[g] * np.sin(self.x_fw_cm / l_fw)
                )
                integrals.append(
                    l_fw * self.integration_constants[n].fw_s[g] * (1 - np.cos(self.x_fw_cm / l_fw))
                )

        return np.sum(integrals, axis=0)

    @summarize_values
    def groupwise_integrated_flux_bz(self, n: int) -> float:
        """
        Calculate the integrated flux[cm s^-1], which can be mulitplied to any
        macroscopic cross-section [cm^-1] to get the reaction rate [s^-1] in
        the blanket.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved. n <= n_groups - 1.
            Therefore n=0 shows the integrated flux for group 1, n=1 for group 2, etc.
        """
        integrals = []
        for g in range(n + 1):
            if self.l_bz_2[g]>0:
                l_bz = np.sqrt(self.l_bz_2[g])
                integrals.append(
                    l_bz * self.integration_constants[n].bz_c[g] *
                    (np.sinh(self.x_bz_cm / l_bz) - np.sinh(self.x_fw_cm / l_bz))
                )
                integrals.append(
                    l_bz * self.integration_constants[n].bz_s[g] *
                    (np.cosh(self.x_bz_cm / l_bz) - np.cosh(self.x_fw_cm / l_bz))
                )
            else:
                l_bz = np.sqrt(-self.l_bz_2[g])
                integrals.append(
                    l_bz * self.integration_constants[n].bz_c[g] *
                    (np.sin(self.x_bz_cm / l_bz) - np.sin(self.x_fw_cm / l_bz))
                )
                integrals.append(
                    -l_bz * self.integration_constants[n].bz_s[g] *
                    (np.cos(self.x_bz_cm / l_bz) - np.cos(self.x_fw_cm / l_bz))
                )
        return np.sum(integrals, axis=0)

    @summarize_values
    def groupwise_neutron_current_fw(self, n: int, x: float | npt.NDArray) -> float:
        """Get the neutron current (in the outward direction) in the fw"""
        x_cm = x * 100
        differentials = []
        for g in range(n + 1):
            if self.l_fw_2[g]>0:
                l_fw = np.sqrt(self.l_fw_2[g])
                differentials.append(
                    self.integration_constants[n].fw_c[g] / l_fw * np.sinh(x_cm / l_fw)
                )
                differentials.append(
                    self.integration_constants[n].fw_s[g] / l_fw * np.cosh(x_cm / l_fw)
                )
            else:
                l_fw = np.sqrt(-self.l_fw_2[g])
                differentials.append(
                    - self.integration_constants[n].fw_c[g] / l_fw * np.sin(x_cm / l_fw)
                )
                differentials.append(
                    self.integration_constants[n].fw_s[g] / l_fw * np.cos(x_cm / l_fw)
                )
        return -self.d_fw_cm[n] * np.sum(differentials, axis=0)

    @summarize_values
    def groupwise_neutron_current_bz(self, n: int, x: float | npt.NDArray) -> float:
        """Get the neutron current (in the outward direction) in the bz."""
        x_cm = x * 100
        differentials = []
        for g in range(n + 1):
            if self.l_bz_2[g]>0:
                l_bz = np.sqrt(self.l_bz_2[g])
                differentials.append(
                    self.integration_constants[n].bz_c[g] / l_bz * np.sinh(x_cm / l_bz)
                )
                differentials.append(
                    self.integration_constants[n].bz_s[g] / l_bz * np.cosh(x_cm / l_bz)
                )
            else:
                l_bz = np.sqrt(-self.l_bz_2[g])
                differentials.append(
                    - self.integration_constants[n].bz_c[g] / l_bz * np.sin(x_cm / l_bz)
                )
                differentials.append(
                    self.integration_constants[n].bz_s[g] / l_bz * np.cos(x_cm / l_bz)
                )
        return -self.d_bz_cm[n] * np.sum(differentials, axis=0)

    @summarize_values
    def groupwise_neutron_current_at(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """
        Neutron current [cm^-2 s^-1] anywhere within the valid range of x,
        i.e. from -self.x_bz_cm to self.x_bz_cm.

        Parameters
        ----------
        n:
            Neutron group index. n <= n_groups - 1.
            Therefore n=0 shows the neutron current for group 1, n=1 for group 2, etc.
        x:
            The depth where we want the neutron flux [m].
            Valid only between x= -extended boundary to +-extended boundary of that group

        Raises
        ------
        ValueError
            The inputted x
        """
        if np.isscalar(x):
            return self.groupwise_neutron_flux_at(n, [x])[0]
        x = np.asarray(x)
        in_fw = abs(x * 100) <= self.x_fw_cm
        in_bz = np.logical_and(
            self.x_fw_cm < abs(x * 100),
            abs(x * 100) <= self.x_bz_cm,
        )
        if (~np.logical_or(in_fw, in_bz)).any():
            raise ValueError(
                f"for neutron group {n}, neutron flux can only be calculated "
                f"up to {self.x_bz_cm} cm, which {x * 100} cm violates!"
            )

        out_current = np.zeros_like(x)
        out_current[in_fw] = self.groupwise_neutron_current_fw(n, x[in_fw])
        out_current[in_bz] = self.groupwise_neutron_current_bz(n, x[in_bz])
        return out_current

    @summarize_values
    def groupwise_neutron_current_fw2bz(self, n: int) -> float:
        """
        Net current from the first wall to breeding zone.
        Parameters
        ----------
        n:
            The index of the neutron group that we want the current for. n <= n_groups - 1.
            Therefore n=0 shows the neutron current for group 1, n=1 for group 2, etc.

        Returns
        -------
        :
            current in cm^-2
        """
        return self.groupwise_neutron_current_bz(n, self.x_fw_cm/100)

    @summarize_values
    def groupwise_neutron_current_escaped(self, n: int) -> float:
        """
        Neutron current escaped from the breeding zone to outside the reactor.
        Parameters
        ----------
        n:
            The index of the neutron group that we want the current for. n <= n_groups - 1.
            Therefore n=0 shows the neutron current for group 1, n=1 for group 2, etc.

        Returns
        -------
        :
            current in cm^-2
        """
        return self.groupwise_neutron_current_bz(n, self.x_bz_cm/100)

    def plot(
        self,
        ax: plt.Axes | None = None,
        n_points: int = 100,
        symmetric: bool = True,
        plot_groups: bool = True,
    ):
        """
        Make a rough plot of the neutron flux.

        Parameters
        ----------
        ax:
            A plt.Axes object to plot on.
        n_points:
            number of points to be used for plotting.
        symmetric:
            Whether to plot from -x to x (symmetric), or from 0 to x (right side only.)
        plot_groups:
            Whether to plot each individual group's neutron flux.
            If True, a legend will be added to help label the groups.
        """
        self.solve_group_n(self.n_groups - 1)
        ax = ax or plt.axes()

        x_bz_left, x_fw, x_bz_right = _generate_x_range(
            min(self.extended_boundary_cm.values()),
            n_points,
            symmetric,
            self.x_fw_cm,
        )
        ax.plot(
            x_bz_left,
            self.neutron_flux_bz(x_bz_left),
            color="black",
            ls=(0, (3, 1, 1, 1)),
        )
        ax.plot(
            x_fw,
            self.neutron_flux_fw(x_fw),
            label="total (FW)",
            color="black",
            ls="solid",
        )
        ax.plot(
            x_bz_right,
            self.neutron_flux_bz(x_bz_right),
            color="black",
            label="total (BZ)",
            ls=(0, (3, 1, 1, 1)),
        )

        if plot_groups:
            for n in range(self.n_groups):
                x_bz_left, x_fw, x_bz_right = _generate_x_range(
                    self.extended_boundary_cm[n],
                    n_points,
                    symmetric,
                    self.x_fw_cm,
                )
                ax.plot(
                    x_bz_left,
                    self.neutron_flux_bz(x_bz_left),
                    color=f"C{n}",
                    ls=(0, (3, 1, 1, 1)),
                )
                ax.plot(
                    x_fw,
                    self.neutron_flux_fw(x_fw),
                    label=f"group {n + 1} (FW)",
                    color=f"C{n}",
                    ls="solid",
                )
                ax.plot(
                    x_bz_right,
                    self.neutron_flux_bz(x_bz_right),
                    label=f"group {n + 1} (BZ)",
                    color=f"C{n}",
                    ls=(0, (3, 1, 1, 1)),
                )
        ax.legend()
        ax.set_title("Neutron flux profile")
        ax.set_xlabel("Distance from the plasma-fw interface [m]")
        ax.set_ylabel("Neutron flux [cm^-2 s^-1]")
        return ax


def _generate_x_range(
    x_max_cm: float,
    approx_n_points: int,
    symmetric: bool,
    fw_bz_split_point_cm: float | None = None,
):
    """Helper function for finding the range of x-values to be plotted.

    Parameters
    ----------
    x_max_cm:
        absolute value of the maximum x that we want to plot.
        This is typically obtained by min(extended_boundary_cm), or
        extended_boundary_cm[n] [cm]
    symmetric:
        Whether we want to plot the negative side of the x-axis as well, forming
        a symmetric plot.
    approx_n_points:
        number of points to be plotted.
    fw_bz_split_point_cm:
        FW and BZ region splits at this distance. If provided, we generate separate
        x ranges for the FW and BZ. [cm]

    Returns
    -------
    if (fw_bz_split_point_cm, symmetric) = (float value, True):
        x_range_bz_left:
            The array of x-values to be used for plotting the left (negative)
            side of bz. [m]
        x_range_fw:
            The array of x-values to be used for plotting both the left and
            right sides of fw. [m]
        x_range_bz_right:
            The array of x-values to be used for plotting the right (positive)
            side of bz. [m]
    elif (fw_bz_split_point_cm, symmetric) = (float value, False):
        x_range_fw:
            The array of x-values to be used for plotting the right (positive)
            side of fw [m]
        x_range_bz:
            The array of x-values to be used for plotting the right (positive)
            side of bz [m]
    elif (fw_bz_split_point_cm, symmetric) = (None, True):
        x_range:
            The array of x-values to be used for plotting both sides of the
            neutron flux. [m]
    else (fw_bz_split_point_cm, symmetric) = (None, False):
        x_range:
            The array of x-values to be used for plotting the right (positive)
            side neutron flux. [m]
    """

    x_max = x_max_cm / 100
    if fw_bz_split_point_cm is not None:
        x_range = np.linspace(0, x_max_cm, approx_n_points)
        n_points_fw = (x_range < fw_bz_split_point_cm).sum() + 1
        n_points_bz = (x_range >= fw_bz_split_point_cm).sum() + 1

        split_point = fw_bz_split_point_cm / 100
        if symmetric:
            return (
                np.linspace(-x_max, -split_point, n_points_bz),
                np.linspace(-split_point, split_point, n_points_fw * 2 - 1),
                np.linspace(split_point, x_max, n_points_bz),
            )
        return (
            np.linspace(0, split_point, n_points_fw),
            np.linspace(split_point, x_max, n_points_bz),
        )

    if symmetric:
        return np.linspace(-x_max, x_max, approx_n_points * 2 - 1)
    return np.linspace(0, x_max, approx_n_points)
