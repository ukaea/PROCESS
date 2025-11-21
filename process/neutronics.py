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
    total_xs: float,
    scattering_xs: float,
    in_source_xs: float,
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
        Unit: m^-1
    scattering_xs:
        macroscopic total cross-section `\sigma_{scatter}`, i.e. number of reactions per
        unit distance travelled by the neutron that leads to it being scattered (without
        getting absorbed).
        Unit: m^-1
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
        unit: [m]
    diffusion_len_2:
        The square of the characteristic diffusion length as given by Reactor Analysis,
        Duderstadt and Hamilton.
        unit: [m^2]
    """

    transport_xs = total_xs - 2 / (3 * avg_atomic_mass) * scattering_xs
    diffusion_coef = 1 / 3 / transport_xs
    diffusion_len_2 = diffusion_coef / (
        total_xs - scattering_xs - in_source_xs
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
class Coefficients:
    """
    Inside each material, there are two new hyperbolic trig funcs per group, i.e.
    group n=0 has cosh(x/L[0]) and sinh(x/L[0]),
    group n=1 has cosh(x/L[0]), sinh(x/L[0]), cosh(x/L[1]) and sinh(x/L[1]),
    etc.
    To get the neutron flux, each trig func has to be scaled by a coeffient.
    E.g.
    group n=0: fw: fw_c[0][0] * cosh(x/L[0]) + fw_s[0][0] * sinh(x/L[0])
    group n=1: fw: fw_c[1][0] * cosh(x/L[0]) + fw_c[1][1] * cosh(x/L[1])
                + fw_s[1][0] * sinh(x/L[0]) + fw_s[1][1] * sinh(x/L[1])

    """

    fw_c: Iterable[float]
    fw_s: Iterable[float]
    bz_c: Iterable[float]
    bz_s: Iterable[float]

    def validate_length(self, exp_len: int, parent_name: str):
        """Validate that all fields has the correct length."""
        for const_name, const_value in asdict(self).items():
            if len(const_value) != exp_len:
                raise ProcessValueError(
                    f"{parent_name}'s [{exp_len-1}]-th item is expected to "
                    "have .{const_name} of length={exp_len}, but instead got "
                    f"{const_value}."
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
        self.name = name
        self._attempting_to_access = set()
        self.populating_method = populating_method

    def __getitem__(self, i: int):
        """Check if index i is in the dictionary or not. If not, populate it."""
        if i in self._attempting_to_access:
            raise RecursionError(
                f"retrieving the value of {self.name}[{i}] requires the "
                f"value of {self.name}[{i}]."
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
            value.validate_length(i + 1, parent_name=self.name)
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
            unit: m^-2 s^-1
        x_fw:
            thickness of the first wall [m].
        x_bz:
            thickness of the blanket + first-wall [m].
        fw_mat:
            first wall material information
        bz_mat:
            blanket material information

        Attributes
        ----------
        fw_mat:
            first wall material information
        bz_mat:
            blanket material information
        n_groups:
            number of groups in the group structure
        group_structure:
            energy bin edges, 1D array of len = n_groups+1

        coefficients:
            Coefficients that determine the flux shape (and therefore
            reaction rates and neutron current) of each group. A set of four
            constants, two for fw and two for bz; each with unit: [m^-2 s^-1]
        l_fw_2:
            square of the characteristic diffusion length as given by Reactor Analysis,
            Duderstadt and Hamilton, applied to the fw. unit: [m^2]
        l_bz_2:
            square of the characteristic diffusion length as given by Reactor Analysis,
            Duderstadt and Hamilton, applied to the bz. unit: [m^2]
        d_fw:
            diffusion constants in the fw. unit: [m]
        d_bz:
            diffusion constants in the bz. unit: [m]
        extended_boundary:
            extended boundary (outside the bz) for each group.
            This value should be larger than x_bz for each of them.
        """
        self.flux = flux  # flux incident on the first wall.
        if not (0 < x_fw < x_bz):
            raise ValueError(
                f"Cannot construct a first-wall+blanket module where{x_fw=}, {x_bz=}."
            )
        self.x_fw, self.x_bz = x_fw, x_bz
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

        self.coefficients = AutoPopulatingDict(
            self.solve_group_n, "coefficients"
        )
        # diffusion lengths squared
        self.l_fw_2 = AutoPopulatingDict(self.solve_group_n, "l_fw_2")
        self.l_bz_2 = AutoPopulatingDict(self.solve_group_n, "l_bz_2")
        self.d_fw = AutoPopulatingDict(self.solve_group_n, "d_fw")
        self.d_bz = AutoPopulatingDict(self.solve_group_n, "d_bz")
        self.extended_boundary = AutoPopulatingDict(
            self.solve_group_n, "extended_boundary"
        )

    def solve_lowest_group(self) -> None:
        """
        Solve the highest-energy (lowest-lethargy)-group's neutron diffusion equation.
        Store the solved constants in self.extended_boundary[0], self.l_fw_2[0],
        self.l_bz_2[0], and self.coefficients[0].
        coefficients have units of [m^-2 s^-1].
        """
        n = 0
        if n in self.coefficients:
            return  # skip if it has already been solved.
        self.d_fw[n], self.l_fw_2[n] = get_diffusion_coefficient_and_length(
            self.fw_mat.avg_atomic_mass,
            self.fw_mat.sigma_t[n],
            self.fw_mat.sigma_s[n, n],
            self.fw_mat.sigma_in[n, n],
        )
        self.d_bz[n], self.l_bz_2[n] = get_diffusion_coefficient_and_length(
            self.bz_mat.avg_atomic_mass,
            self.bz_mat.sigma_t[n],
            self.bz_mat.sigma_s[n, n],
            self.bz_mat.sigma_in[n, n],
        )
        l_fw = np.sqrt(abs(self.l_fw_2[n]))
        l_bz = np.sqrt(abs(self.l_bz_2[n]))
        x_fw, x_bz = self.x_fw, self.x_bz
        d_fw = self.d_fw[n]
        d_bz = self.d_bz[n]
        self.extended_boundary[n] = x_bz + extrapolation_length(d_bz)
        if self.l_fw_2[n] > 0:
            s_fw = np.sinh(x_fw / l_fw)
            c_fw = np.cosh(x_fw / l_fw)
            t_fw = np.tanh(x_fw / l_fw)
        else:
            s_fw = np.sin(x_fw / l_fw)
            c_fw = np.cos(x_fw / l_fw)
            t_fw = np.tan(x_fw / l_fw)
        if self.l_bz_2[n] > 0:
            c_bz = np.cosh(self.extended_boundary[n] / l_bz)
            s_bz = np.sinh(self.extended_boundary[n] / l_bz)
            c_bz_mod = np.cosh((self.extended_boundary[n] - x_fw) / l_bz)
            s_bz_mod = np.sinh((self.extended_boundary[n] - x_fw) / l_bz)
            t_bz_mod = np.tanh((self.extended_boundary[n] - x_fw) / l_bz)
        else:
            c_bz = np.cos(self.extended_boundary[n] / l_bz)
            s_bz = np.sin(self.extended_boundary[n] / l_bz)
            c_bz_mod = np.cos((self.extended_boundary[n] - x_fw) / l_bz)
            s_bz_mod = np.sin((self.extended_boundary[n] - x_fw) / l_bz)
            t_bz_mod = np.tan((self.extended_boundary[n] - x_fw) / l_bz)

        fw_c_factor = -self.flux * (
            l_fw / d_fw
            - np.exp(x_fw / l_fw)
            * ((l_fw / d_fw) + (l_bz / d_bz) * t_bz_mod)
            / (c_fw + s_fw * t_bz_mod * (d_fw / l_fw) * (l_bz / d_bz))
        )
        fw_s_factor = -self.flux * l_fw / d_fw

        bz_common_factor = (
            self.flux
            * np.exp(x_fw / l_fw)
            * (1 - t_fw)
            / ((d_bz / l_bz) * c_bz_mod + (d_fw / l_fw) * t_fw * s_bz_mod)
        )
        bz_c_factor = bz_common_factor * s_bz
        bz_s_factor = -bz_common_factor * c_bz

        self.coefficients[n] = Coefficients(
            [fw_c_factor], [fw_s_factor], [bz_c_factor], [bz_s_factor]
        )
        return

    def solve_group_n(self, n: int) -> None:
        """
        Solve the n-th group of neutron's diffusion equation, where n<=n_groups-1.
        Store the solved constants in self.extended_boundary[n-1], self.l_fw_2[n-1],
        self.l_bz_2[n-1], and self.coefficients[n-1].

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
            if k not in self.coefficients:
                self.solve_group_n(k)
        if not (self.fw_mat.downscatter_only and self.bz_mat.downscatter_only):
            raise NotImplementedError(
                "Will implement solve_group_n in a loop later..."
            )
        if n in self.coefficients:
            return None  # skip if it has already been solved.
        # Parameter to be changed later, to allow solving non-down-scatter
        # only systems by iterating.
        first_iteration = True
        self.d_fw[n], self.l_fw_2[n] = get_diffusion_coefficient_and_length(
            self.fw_mat.avg_atomic_mass,
            self.fw_mat.sigma_t[n],
            self.fw_mat.sigma_s[n, n],
            self.fw_mat.sigma_in[n, n],
        )
        self.d_bz[n], self.l_bz_2[n] = get_diffusion_coefficient_and_length(
            self.bz_mat.avg_atomic_mass,
            self.bz_mat.sigma_t[n],
            self.bz_mat.sigma_s[n, n],
            self.bz_mat.sigma_in[n, n],
        )
        self.extended_boundary[n] = self.x_bz + extrapolation_length(
            self.d_bz[n]
        )

        # Setting up aliases for shorter code
        l_fw_2, l_bz_2 = self.l_fw_2[n], self.l_bz_2[n]
        d_fw, d_bz = self.d_fw[n], self.d_bz[n]
        src_fw = self.fw_mat.sigma_s + self.fw_mat.sigma_in
        src_bz = self.bz_mat.sigma_s + self.bz_mat.sigma_in

        self.coefficients[n] = Coefficients([], [], [], [])
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
            self.coefficients[n].fw_c.append(
                sum(
                    src_fw[i, n] * self.coefficients[i].fw_c[g]
                    for i in range(g, n)
                )
                * scale_factor_fw
            )
            self.coefficients[n].fw_s.append(
                sum(
                    src_fw[i, n] * self.coefficients[i].fw_s[g]
                    for i in range(g, n)
                )
                * scale_factor_fw
            )
            self.coefficients[n].bz_c.append(
                sum(
                    src_bz[i, n] * self.coefficients[i].bz_c[g]
                    for i in range(g, n)
                )
                * scale_factor_bz
            )
            self.coefficients[n].bz_s.append(
                sum(
                    src_bz[i, n] * self.coefficients[i].bz_s[g]
                    for i in range(g, n)
                )
                * scale_factor_bz
            )

        fw_c_factor_guess = 0.0
        fw_s_factor_guess = 0.0
        bz_c_factor_guess = 0.0
        bz_s_factor_guess = 0.0
        self.coefficients[n].fw_c.append(fw_c_factor_guess)
        self.coefficients[n].fw_s.append(fw_s_factor_guess)
        self.coefficients[n].bz_c.append(bz_c_factor_guess)
        self.coefficients[n].bz_s.append(bz_s_factor_guess)

        def _set_constants(input_vector: Iterable[float]):
            self.coefficients[n].fw_c[n] = input_vector[0]
            self.coefficients[n].fw_s[n] = input_vector[1]
            self.coefficients[n].bz_c[n] = input_vector[2]
            self.coefficients[n].bz_s[n] = input_vector[3]

        def _evaluate_fit():
            flux_continuity = self.groupwise_neutron_flux_fw(
                n, self.x_fw
            ) - self.groupwise_neutron_flux_bz(n, self.x_fw)
            flux_at_boundary = self.groupwise_neutron_flux_bz(
                n, self.extended_boundary[n]
            )
            current_continuity = self.groupwise_neutron_current_fw(
                n, self.x_fw
            ) - self.groupwise_neutron_current_bz(n, self.x_fw)
            influx_fw, influx_bz = 0.0, 0.0
            for g in range(self.n_groups):
                if g > n and not first_iteration:
                    if not self.fw_mat.downscatter_only:
                        influx_fw += (
                            self.fw_mat.sigma_s[g, n]
                            + self.fw_mat.sigma_in[g, n]
                        ) * self.groupwise_integrated_flux_fw(g)
                    if not self.bz_mat.downscatter_only:
                        influx_bz += (
                            self.bz_mat.sigma_s[g, n]
                            + self.bz_mat.sigma_in[g, n]
                        ) * self.groupwise_integrated_flux_bz(g)
                elif g <= n:
                    influx_fw += (
                        self.fw_mat.sigma_s[g, n] + self.fw_mat.sigma_in[g, n]
                    ) * self.groupwise_integrated_flux_fw(g)
                    influx_bz += (
                        self.bz_mat.sigma_s[g, n] + self.bz_mat.sigma_in[g, n]
                    ) * self.groupwise_integrated_flux_bz(g)
            removal_fw = self.fw_mat.sigma_t[
                n
            ] * self.groupwise_integrated_flux_fw(n)
            removal_bz = self.bz_mat.sigma_t[
                n
            ] * self.groupwise_integrated_flux_bz(n)
            # conservation_fw = influx_fw - removal_fw - current_fw2bz
            # conservation_bz = current_fw2bz + influx_bz - removal_bz - escaped_bz
            total_neutron_conservation = (
                influx_fw
                + influx_bz
                - removal_fw
                - removal_bz
                - self.groupwise_neutron_current_escaped(n)
            )
            return np.array([
                flux_continuity,
                flux_at_boundary,
                current_continuity,
                total_neutron_conservation,
            ])

        def objective(four_coefficients_vector):
            _set_constants(four_coefficients_vector)
            return _evaluate_fit()

        results = optimize.root(
            objective,
            x0=[
                fw_c_factor_guess,
                fw_s_factor_guess,
                bz_c_factor_guess,
                bz_s_factor_guess,
            ],
        )
        _set_constants(results.res)
        return None

    @summarize_values
    def groupwise_neutron_flux_fw(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """
        Neutron flux[m^-2 s^-1] of the n-th group at the first wall, at location x [m].

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
        trig_funcs = []
        for g in range(n + 1):
            if self.l_fw_2[g] > 0:
                c, s = np.cosh, np.sinh
                l_fw = np.sqrt(self.l_fw_2[g])
            else:
                c, s = np.cos, np.sin
                l_fw = np.sqrt(-self.l_fw_2[g])
            trig_funcs.append(
                self.coefficients[n].fw_c[g] * c(abs(x) / l_fw)
            )
            trig_funcs.append(
                self.coefficients[n].fw_s[g] * s(abs(x) / l_fw)
            )
        return np.sum(trig_funcs, axis=0)

    @summarize_values
    def groupwise_neutron_flux_bz(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """Neutron flux[m^-2 s^-1] of the n-th groupat the blanket, at location x [m].

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
        trig_funcs = []
        for g in range(n + 1):
            if self.l_bz_2[g] > 0:
                c, s = np.cosh, np.sinh
                l_bz = np.sqrt(self.l_bz_2[g])
            else:
                c, s = np.cos, np.sin
                l_bz = np.sqrt(-self.l_bz_2[g])
            trig_funcs.append(
                self.coefficients[n].bz_c[g] * c(abs(x) / l_bz)
            )
            trig_funcs.append(
                self.coefficients[n].bz_s[g] * s(abs(x) / l_bz)
            )
        return np.sum(trig_funcs, axis=0)

    @summarize_values
    def groupwise_neutron_flux_at(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """
        Neutron flux [m^-2 s^-1] anywhere within the valid range of x,
        i.e. [-self.extended_boundary[n], self.extended_boundary[n]].

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
        in_fw = abs(x) <= self.x_fw
        in_bz = np.logical_and(
            self.x_fw < abs(x),
            abs(x) <= self.extended_boundary[n],
        )
        if (~np.logical_or(in_fw, in_bz)).any():
            raise ValueError(
                f"for neutron group {n}, neutron flux can only be calculated "
                f"up to {self.extended_boundary[n]} m, which {x} m violates!"
            )

        out_flux = np.zeros_like(x)
        out_flux[in_fw] = self.groupwise_neutron_flux_fw(n, x[in_fw])
        out_flux[in_bz] = self.groupwise_neutron_flux_bz(n, x[in_bz])
        return out_flux

    # scalar values (one such float per neutron group.)
    @summarize_values
    def groupwise_integrated_flux_fw(self, n: int) -> float:
        """
        Calculate the integrated flux[m^-1 s^-1], which can be mulitplied to any
        macroscopic cross-section [m^-1] to get the reaction rate [s^-1] in
        the first wall.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved. n <= n_groups - 1.
            Therefore n=0 shows the integrated flux for group 1, n=1 for group 2, etc.
        """
        integrals = []
        for g in range(n + 1):
            if self.l_fw_2[g] > 0:
                l_fw = np.sqrt(self.l_fw_2[g])
                integrals.append(
                    l_fw
                    * self.coefficients[n].fw_c[g]
                    * np.sinh(self.x_fw / l_fw)
                )
                integrals.append(
                    l_fw
                    * self.coefficients[n].fw_s[g]
                    * (np.cosh(self.x_fw / l_fw) - 1)
                )
            else:
                l_fw = np.sqrt(-self.l_fw_2[g])
                integrals.append(
                    l_fw
                    * self.coefficients[n].fw_c[g]
                    * np.sin(self.x_fw / l_fw)
                )
                integrals.append(
                    l_fw
                    * self.coefficients[n].fw_s[g]
                    * (1 - np.cos(self.x_fw / l_fw))
                )

        return np.sum(integrals, axis=0)

    @summarize_values
    def groupwise_integrated_flux_bz(self, n: int) -> float:
        """
        Calculate the integrated flux[m^-1 s^-1], which can be mulitplied to any
        macroscopic cross-section [m^-1] to get the reaction rate [s^-1] in
        the blanket.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved. n <= n_groups - 1.
            Therefore n=0 shows the integrated flux for group 1, n=1 for group 2, etc.
        """
        integrals = []
        for g in range(n + 1):
            if self.l_bz_2[g] > 0:
                l_bz = np.sqrt(self.l_bz_2[g])
                integrals.append(
                    l_bz
                    * self.coefficients[n].bz_c[g]
                    * (np.sinh(self.x_bz / l_bz) - np.sinh(self.x_fw / l_bz))
                )
                integrals.append(
                    l_bz
                    * self.coefficients[n].bz_s[g]
                    * (np.cosh(self.x_bz / l_bz) - np.cosh(self.x_fw / l_bz))
                )
            else:
                l_bz = np.sqrt(-self.l_bz_2[g])
                integrals.append(
                    l_bz
                    * self.coefficients[n].bz_c[g]
                    * (np.sin(self.x_bz / l_bz) - np.sin(self.x_fw / l_bz))
                )
                integrals.append(
                    -l_bz
                    * self.coefficients[n].bz_s[g]
                    * (np.cos(self.x_bz / l_bz) - np.cos(self.x_fw / l_bz))
                )
        return np.sum(integrals, axis=0)

    @summarize_values
    def groupwise_neutron_current_fw(
        self, n: int, x: float | npt.NDArray
    ) -> float:
        """Get the neutron current (in the outward direction) in the fw"""
        differentials = []
        for g in range(n + 1):
            if self.l_fw_2[g] > 0:
                l_fw = np.sqrt(self.l_fw_2[g])
                differentials.append(
                    self.coefficients[n].fw_c[g]
                    / l_fw
                    * np.sinh(x / l_fw)
                )
                differentials.append(
                    self.coefficients[n].fw_s[g]
                    / l_fw
                    * np.cosh(x / l_fw)
                )
            else:
                l_fw = np.sqrt(-self.l_fw_2[g])
                differentials.append(
                    -self.coefficients[n].fw_c[g]
                    / l_fw
                    * np.sin(x / l_fw)
                )
                differentials.append(
                    self.coefficients[n].fw_s[g]
                    / l_fw
                    * np.cos(x / l_fw)
                )
        return -self.d_fw[n] * np.sum(differentials, axis=0)

    @summarize_values
    def groupwise_neutron_current_bz(
        self, n: int, x: float | npt.NDArray
    ) -> float:
        """Get the neutron current (in the outward direction) in the bz."""
        differentials = []
        for g in range(n + 1):
            if self.l_bz_2[g] > 0:
                l_bz = np.sqrt(self.l_bz_2[g])
                differentials.append(
                    self.coefficients[n].bz_c[g]
                    / l_bz
                    * np.sinh(x / l_bz)
                )
                differentials.append(
                    self.coefficients[n].bz_s[g]
                    / l_bz
                    * np.cosh(x / l_bz)
                )
            else:
                l_bz = np.sqrt(-self.l_bz_2[g])
                differentials.append(
                    -self.coefficients[n].bz_c[g]
                    / l_bz
                    * np.sin(x / l_bz)
                )
                differentials.append(
                    self.coefficients[n].bz_s[g]
                    / l_bz
                    * np.cos(x / l_bz)
                )
        return -self.d_bz[n] * np.sum(differentials, axis=0)

    @summarize_values
    def groupwise_neutron_current_at(
        self, n: int, x: float | npt.NDArray
    ) -> npt.NDArray:
        """
        Neutron current [m^-2 s^-1] anywhere within the valid range of x,
        i.e. from -self.x_bz to self.x_bz.

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
        in_fw = abs(x) <= self.x_fw
        in_bz = np.logical_and(
            self.x_fw < abs(x),
            abs(x) <= self.x_bz,
        )
        if (~np.logical_or(in_fw, in_bz)).any():
            raise ValueError(
                f"for neutron group {n}, neutron flux can only be calculated "
                f"up to {self.x_bz} m, which {x} m violates!"
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
            current in m^-2
        """
        return self.groupwise_neutron_current_bz(n, self.x_fw)

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
            current in m^-2
        """
        return self.groupwise_neutron_current_bz(n, self.x_bz)

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
            min(self.extended_boundary.values()),
            n_points,
            symmetric,
            self.x_fw,
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
                    self.extended_boundary[n],
                    n_points,
                    symmetric,
                    self.x_fw,
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
        ax.set_ylabel("Neutron flux [m^-2 s^-1]")
        return ax


def _generate_x_range(
    x_max: float,
    approx_n_points: int,
    symmetric: bool,
    fw_bz_split_point: float | None = None,
):
    """Helper function for finding the range of x-values to be plotted.

    Parameters
    ----------
    x_max:
        absolute value of the maximum x that we want to plot.
        This is typically obtained by min(extended_boundary), or
        extended_boundary[n] [m]
    symmetric:
        Whether we want to plot the negative side of the x-axis as well, forming
        a symmetric plot.
    approx_n_points:
        number of points to be plotted.
    fw_bz_split_point:
        FW and BZ region splits at this distance. If provided, we generate separate
        x ranges for the FW and BZ. [m]

    Returns
    -------
    if (fw_bz_split_point, symmetric) = (float value, True):
        x_range_bz_left:
            The array of x-values to be used for plotting the left (negative)
            side of bz. [m]
        x_range_fw:
            The array of x-values to be used for plotting both the left and
            right sides of fw. [m]
        x_range_bz_right:
            The array of x-values to be used for plotting the right (positive)
            side of bz. [m]
    elif (fw_bz_split_point, symmetric) = (float value, False):
        x_range_fw:
            The array of x-values to be used for plotting the right (positive)
            side of fw [m]
        x_range_bz:
            The array of x-values to be used for plotting the right (positive)
            side of bz [m]
    elif (fw_bz_split_point, symmetric) = (None, True):
        x_range:
            The array of x-values to be used for plotting both sides of the
            neutron flux. [m]
    else (fw_bz_split_point, symmetric) = (None, False):
        x_range:
            The array of x-values to be used for plotting the right (positive)
            side neutron flux. [m]
    """

    if fw_bz_split_point is not None:
        x_range = np.linspace(0, x_max, approx_n_points)
        n_points_fw = (x_range < fw_bz_split_point).sum() + 1
        n_points_bz = (x_range >= fw_bz_split_point).sum() + 1

        if symmetric:
            return (
                np.linspace(-x_max, -fw_bz_split_point, n_points_bz),
                np.linspace(
                    -fw_bz_split_point, fw_bz_split_point, n_points_fw * 2 - 1
                ),
                np.linspace(fw_bz_split_point, x_max, n_points_bz),
            )
        return (
            np.linspace(0, fw_bz_split_point, n_points_fw),
            np.linspace(fw_bz_split_point, x_max, n_points_bz),
        )

    if symmetric:
        return np.linspace(-x_max, x_max, approx_n_points * 2 - 1)
    return np.linspace(0, x_max, approx_n_points)
