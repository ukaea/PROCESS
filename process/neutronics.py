"""
Most of derivation is basd on the book Reactor Analysis, Duderstadt and Hamilton, 1976,
ISBN:9780471223634.
"""

import functools
import inspect
from itertools import pairwise
from collections.abc import Callable, Iterable
from dataclasses import asdict, dataclass

import numpy as np
from matplotlib import pyplot as plt
from numpy import typing as npt
from scipy import optimize

from process.exceptions import ProcessValidationError, ProcessValueError
from process.neutronics_data import MaterialMacroInfo, DT_NEUTRON_E


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
            "flux/integrated flux/current methods."
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
    diffusion_const:
        The diffusion coefficient as given by Reactor Analysis, Duderstadt and Hamilton.
        unit: [m]
    diffusion_len_2:
        The square of the characteristic diffusion length as given by Reactor Analysis,
        Duderstadt and Hamilton.
        unit: [m^2]
    """

    transport_xs = total_xs - 2 / (3 * avg_atomic_mass) * scattering_xs
    diffusion_const = 1 / 3 / transport_xs
    diffusion_len_2 = diffusion_const / (
        total_xs - scattering_xs - in_source_xs
    )
    return diffusion_const, diffusion_len_2


def extrapolation_length(diffusion_constficient: float) -> float:
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
    return 0.7104 * 3 * diffusion_constficient


@dataclass
class Coefficients:
    """
    Inside each material, there are two new hyperbolic trig funcs per group, i.e.
    group n=0 has cosh(x/L[0]) and sinh(x/L[0]),
    group n=1 has cosh(x/L[0]), sinh(x/L[0]), cosh(x/L[1]) and sinh(x/L[1]),
    etc.
    To get the neutron flux, each trig func has to be scaled by a coefficient.
    E.g. Let's say for the first wall, which is num_layer=0,
    group n=0: NeutronFluxProfile.coefficients[0, 0] = Coefficients(...)
        fw_grp0 = NeutronFluxProfile.coefficients[0, 0]
        flux = fw_grp0.c[0] * cosh(x/L[0]) + fw_grp0.s[0] * sinh(x/L[0])
    
    group n=1: NeutronFluxProfile.coefficients[0, 1] = Coefficients(...)
        fw_grp1 = NeutronFluxProfile.coefficients[0, 1]
        flux = fw_grp1.c[0] * cosh(x/L[0]) + fw_grp1.c[1] * cosh(x/L[1])
                + fw_grp1.s[0] * sinh(x/L[0]) + fw_grp1.s[1] * sinh(x/L[1])

    """

    c: Iterable[float]
    s: Iterable[float]

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
        return f"AutoPopulatingDict({self.name}):{self._dict}"

class LayerSpecificGroupwiseConstants:
    """An object containing multiple AutoPopulatingDict"""
    def __init__(self, populating_method: Callable[[int], None], layer_names: list[str], quantity_description: str):
        """
        Create an object that contains as many AutoPopulatingDict as there are
        items in layer_names.

        Parameters
        ----------
        populating_method:
            The method to be called if the requested index n in any one of the
            AutoPopulatingDict is currently unpopulated. This method should
            populate that dictionary.
        layer_names:
            A list of strings, each of which is the descriptive name for that
            layer. While the actual content in each string could be empty,
            the length of this list MUST be equal to the total number
            of dictionaries required.
        quantity_description:
            A name to be given to this specific instance of the class, to help
            label what quantity is being stored.
        """
        self._name = quantity_description
        layer_dicts = []
        for num_layer, _layer_name in enumerate(layer_names):
            name = f"{self._name} for layer {num_layer}"
            if _layer_name:
                name += f":{_layer_name}"
            layer_dicts.append(AutoPopulatingDict(populating_method, name))
        self._dicts = tuple(layer_dicts)
        self.n_layers = len(self._dicts)

    def __iter__(self):
        return self._dicts.__iter__()

    def __len__(self) -> int:
        return len(self._dicts)
        
    def __setitem__(self, index: int | tuple[int, int], value):
        """
        Act as if this is a 2D array, where the first-axis is the layer and the
        second axis is the group. support slice of the thing, 
        """
        if isinstance(index, tuple) and len(index)>=2:
            if len(index)>2:
                raise IndexError("2D array indexed with more than 2 indices!")
            layer_index, group_index = index
            self._dicts[layer_index][group_index] = value
        else:
            super().__setitem__(index, value)

    def __getitem__(self, index: int | tuple[int, int]):
        """
        Act as if this is a 2D array, where the first-axis is the layer and the
        second axis is the group. Handle slices as well.
        """
        if isinstance(index, tuple) and len(index)>=2:
            if len(index)>2:
                raise IndexError("2D array indexed with more than 2 indices!")
            layer_index, group_index = index
            if isinstance(layer_index, slice):
                return tuple(
                    [_dict[group_index] for _dict in self._dicts[layer_index]]
                )
            return self._dicts[layer_index][group_index]
        return self._dicts[index]

    def has_populated(self, n: int) -> bool:
        """
        Check if group n's constants are populated for every layer's dict.

        Parameter
        ---------
        n:
            group number to check the population status of every layer's dict.


        Returns
        -------
        :
            Whether group n's constants are poplated across all dictionaries
            stored by the current instance.
        """
        return all(n in layer_dict for layer_dict in self._dicts)

    def __repr__(self):
        return f"A tuple of {self.n_layers} layers of {self._name}"

UNIT_LOOKUP = {
    "integrated_flux": "m^-1 s^-1",
    "reaction_rate": "m^-2 s^-1",
    "flux": "m^-2 s^-1",
    "current": "m^-2 s^-1",
    "heating": "W m^-3",
}

class NeutronFluxProfile:
    """
    Calculate the neutron flux, neutron current, and neutron heating in the
    mirrored infinite-slab model, where each layer extend infinitely in y- and
    z-directions, but has finite width in the x-direction. Each layer's
    thickness is defined along the positive x-axis starting at 0, and then
    reflected along x=0 to fill out the negative x-axis.
    """

    def __init__(
        self,
        flux: float,
        layer_x: npt.NDArray[np.float64],
        materials: Iterable[MaterialMacroInfo],
        init_neutron_energy: float = DT_NEUTRON_E,
    ):
        """Initialize a particular FW-BZ geometry and neutron flux.

        Parameters
        ----------
        flux:
            Neutron flux directly emitted by the plasma, incident on the first
            wall. unit: m^-2 s^-1
        init_neutron_energy:
            Neutron's initial energy when it first exit the plasma, before any
            downscattering or reactions. unit: J.
        layer_x:
            The x-coordinates of the right side of every layers. By definition,
            the plasma is situated at x=0, so all values in layer_x must be >0.
            E.g. layer_x[0] is the thickness of the first wall,
            layer_x[1] is the thickness of the first wall + breeding zone,
            etc.
        materials:
            Every layer's material information.

        Attributes
        ----------
        interface_x:
            The x-coordinates of every plasma-layer interfaces/ layer-layer
            interface/ layer-void interface. For n_layers,
            there will be the interface between the first layer and the plasma,
            plus (n_layers - 1) interfaces between layers, plus the interface
            between the final layer and the void into which neutrons are lost.
            E.g. interface_x[0] = 0.0 = the plasma-fw interface.
        n_layers:
            Number of layers
        n_groups:
            Number of groups in the group structure
        group_structure:
            Energy bin edges, 1D array of len = n_groups+1
        energy:
            The mean neutron energy of each group.

        coefficients:
            Coefficients that determine the flux shape (and therefore reaction
            rates, neutron current, etc.) of each group. Each coefficient has
            unit: [m^-2 s^-1]
        l2:
            Square of the characteristic diffusion length of each layer as
            given by Reactor Analysis, Duderstadt and Hamilton. unit: [m^2]
        diffusion_const:
            Diffusion coefficient of each layer. unit: [m]
        extended_boundary:
            Extended boundary for each group. These values should be larger
            than layer_x[-1].
        """
        # flux incident on the first wall at the highest energy.
        self.flux = flux
        self.init_neutron_energy = init_neutron_energy

        # layers
        self.layer_x = np.array(layer_x).ravel()
        if not (np.diff(self.layer_x)>0).all():
            raise ValueError(
                f"Model cannot have non-positive layer thicknesses."
            )

        self.layer_x.flags.writeable = False
        self.interface_x = np.array([0.0, *self.layer_x])
        self.interface_x.flags.writeable = False

        self.materials = tuple(materials)
        if len(self.layer_x) != len(self.materials):
            raise ProcessValidationError(
                "The number of layers specified by self.materials must match "
                "the number of x-positions specified by layer_x."
            )
        self.n_layers = len(self.materials)

        # groups
        fw_mat = self.materials[0]
        for mat in self.materials[1:]:
            if not np.allclose(
                fw_mat.group_structure, mat.group_structure, atol=0,
            ):
                raise ProcessValidationError(
                    "All material info must have the same group structure!"
                )
        self.n_groups = fw_mat.n_groups
        self.group_structure = fw_mat.group_structure
        self.group_energy = self._calculate_mean_energy(
            self.group_structure, self.init_neutron_energy
        )

        mat_name_list = [mat.name for mat in self.materials]
        self.coefficients = LayerSpecificGroupwiseConstants(
            self.solve_group_n, mat_name_list, "Coefficients"
        )
        self.l2 = LayerSpecificGroupwiseConstants(
            self.solve_group_n, mat_name_list,
            "Characteristic diffusion length squared"
        )
        self.diffusion_const = LayerSpecificGroupwiseConstants(
            self.solve_group_n, mat_name_list, "Diffusion coefficient D"
        )
        self.extended_boundary = AutoPopulatingDict(
            self.solve_group_n, "extended_boundary"
        )

    @staticmethod
    def _calculate_mean_energy(
        group_structure: npt.NDArray[np.float64], init_neutron_e: float,
    ) -> npt.NDArray[np.float64]:
        """
        Calculate the mean energy of each neutron group in Joule.
        When implementing this method, we can choose either a weighted mean or
        an unweighted mean. The weighted mean (where neutron flux is assumed
        constant w.r.t. neutorn lethargy within the bin) is more accurate, but
        may end up being incorrect if the bins in the group structure are too
        wide, as it heavily biases towards lower energy. In contrast, the
        simple unweighted mean does not have such problem, but is inconsistent
        with the rest of the program's assumption (const. flux w.r.t. lethargy),
        therefore the former is chosen.

        Parameters
        ----------
        group_structure:
            The neutron energy bin's edges, in descending order, with len=
            n_groups + 1.
        init_neutron_e:
            The neutrons entering from the plasma to the FW is assumed to be
            monoenergetic, with this energy.

        Returns
        -------
        mean_neutron_e:
            Mean energy of neutrons. The bin containing init_neutron_e
            (likely the highest energy bin, i.e. bin[0]) is assumed to be
            dominated by the unscattered neutrons entering from the plasma,
            therefore it is 
        """
        high, low = group_structure[:-1], group_structure[1:]
        weighted_mean = (high - low)/(np.log(high)- np.log(low))
        # unweighted_mean = np.mean([high, low], axis=0)
        first_bin = np.logical_and(high>init_neutron_e, init_neutron_e>=low)
        if first_bin.sum()<1:
            raise ValueError(
                "The energy of neutrons incident from the plasma is not "
                "captured by the group structure!"
            )
        elif first_bin.sum()>1:
            raise ValueError(
                "The energy of neutrons incident from the plasma is NOT captured by the group "
                "structure!"
            )
        weighted_mean[first_bin] = init_neutron_e
        return weighted_mean

    def solve_lowest_group(self) -> None:
        """
        Solve the highest-energy (lowest-lethargy)-group's neutron diffusion equation.
        Store the solved constants in self.extended_boundary[0], self.l2[:, 0],
        and self.coefficients[:, 0].
        self.coefficients each have units of [m^-2 s^-1].
        """
        n = 0
        if self.coefficients.has_populated(n):
            return  # skip if it has already been solved.
        for num_layer, mat in enumerate(self.materials):
            self.diffusion_const[num_layer, n], self.l2[num_layer, n] = get_diffusion_coefficient_and_length(
                mat.avg_atomic_mass,
                mat.sigma_t[n],
                mat.sigma_s[n, n],
                mat.sigma_in[n, n],
            )
        self.extended_boundary[n] = self.layer_x[-1] + extrapolation_length(
            self.diffusion_const[-1, n]
        )
        if self.n_layers == 2:
            l_fw, l_bz = np.sqrt(abs(self.l2[0, n])), np.sqrt(abs(self.l2[1, n]))
            x_fw, = self.layer_x[:-1]
            d_fw, d_bz = self.diffusion_const[0, n], self.diffusion_const[1, n]
            if self.l2[0, n] > 0:
                s_fw = np.sinh(x_fw / l_fw)
                c_fw = np.cosh(x_fw / l_fw)
                t_fw = np.tanh(x_fw / l_fw)
            else:
                s_fw = np.sin(x_fw / l_fw)
                c_fw = np.cos(x_fw / l_fw)
                t_fw = np.tan(x_fw / l_fw)
            if self.l2[1, n] > 0:
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

            self.coefficients[0, n] = Coefficients([fw_c_factor], [fw_s_factor])
            self.coefficients[1, n] = Coefficients([bz_c_factor], [bz_s_factor])
        else:
            raise NotImplementedError("Only implemented 2 groups so far.")
        return

    def solve_group_n(self, n: int) -> None:
        """
        Solve the n-th group of neutron's diffusion equation, where n <=
        n_groups-1. Store the solved constants in self.extended_boundary[n-1],
        self.l2[:, n-1], and self.coefficients[:, n-1].

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
            if not self.coefficients.has_populated(k):
                self.solve_group_n(k)
        if not all(mat.downscatter_only for mat in self.materials):
            raise NotImplementedError(
                "Will implement solve_group_n in a loop later..."
            )
        if self.coefficients.has_populated(n):
            return None  # skip if it has already been solved.
        # Parameter to be changed later, to allow solving non-down-scatter
        # only systems by iterating.
        first_iteration = True
        for num_layer, mat in enumerate(self.materials):
            self.diffusion_const[num_layer, n], self.l2[num_layer, n] = get_diffusion_coefficient_and_length(
                mat.avg_atomic_mass,
                mat.sigma_t[n],
                mat.sigma_s[n, n],
                mat.sigma_in[n, n],
            )
        self.extended_boundary[n] = self.layer_x[-1] + extrapolation_length(
            self.diffusion_const[-1, n]
        )


        for num_layer in range(self.n_layers):
            # Setting up aliases for shorter code
            _coefs = Coefficients([], [])
            mat = self.materials[num_layer]
            src_matrix = mat.sigma_s + mat.sigma_in
            diffusion_const_n = self.diffusion_const[num_layer, n]
            for g in range(n):
                # if the characteristic length of group [g] coincides with the
                # characteristic length of group [n], then that particular
                # cosh/sinh would be indistinguishable from group [n]'s
                # cosh/sinh anyways, therefore we can set the coefficient to 0.
                scale_factor = 0.0
                l2n, l2g = self.l2[num_layer, n], self.l2[num_layer, g]
                if not np.isclose(l2_diff:=(l2g - l2n), 0):
                    scale_factor = (l2n * l2g) / diffusion_const_n / l2_diff
                _coefs.c.append(
                    sum(
                        src_matrix[i, n] * self.coefficients[num_layer, i].c[g]
                        for i in range(g, n)
                    )
                    * scale_factor
                )
                _coefs.s.append(
                    sum(
                        src_matrix[i, n] * self.coefficients[num_layer, i].s[g]
                        for i in range(g, n)
                    )
                    * scale_factor
                )

            c_factor_guess = 0.0
            s_factor_guess = 0.0
            _coefs.c.append(c_factor_guess)
            _coefs.s.append(s_factor_guess)
            self.coefficients[num_layer, n] = _coefs


        def _set_coefficients(input_vector: Iterable[float]):
            for num_layer in range(self.n_layers):
                i = num_layer * 2
                self.coefficients[num_layer, n].c[n] = input_vector[i]
                self.coefficients[num_layer, n].s[n] = input_vector[i+1]

        def _evaluate_fit():
            # zero flux condition 
            flux_at_boundary = self.groupwise_neutron_flux_in_layer(
                n, self.n_layers - 1, self.extended_boundary[n]
            )
            # conservation condition
            influx, removal = 0.0, 0.0
            for g in range(self.n_groups):
                if g > n and not first_iteration:
                    for num_layer, mat in self.materials:
                        if not mat.downscatter_only:
                            influx += (mat.sigma_s[g, n] + mat.sigma_in[g, n]) * self.groupwise_integrated_flux_in_layer(g, num_layer)
                elif g <= n:
                    for num_layer, mat in self.materials:
                        influx += (mat.sigma_s[g, n] + mat.sigma_in[g, n]) * self.groupwise_integrated_flux_in_layer(g, num_layer)
            for num_layer, mat in self.materials:
                removal += mat.sigma_t[n] * self.groupwise_integrated_flux_in_layer(n, num_layer)
            # conservation_fw = influx_fw - removal_fw - current_fw2bz
            # conservation_bz = current_fw2bz + influx_bz - removal_bz - escaped_bz
            total_neutron_conservation = influx - removal - self.groupwise_neutron_current_escaped(n)

            conditions = [flux_at_boundary, total_neutron_conservation]

            # enforce continuity conditions
            for num_layer in range(self.n_layers - 1):
                x = self.layer_x[num_layer]
                flux_continuity = self.groupwise_neutron_flux_in_layer(n, num_layer, x) - self.groupwise_neutron_flux_in_layer(n, num_layer + 1, x)
                current_continuity = self.groupwise_neutron_current_in_layer(n, num_layer, x) - self.groupwise_neutron_current_in_layer(n, num_layer + 1, x)
            conditions.append(flux_continuity)
            conditions.append(current_continuity)
            # TODO: there may be one more condition that I should impose: the slope at x=0 should be 0.
            return np.array(conditions)

        def objective(coefficients_vector):
            _set_coefficients(coefficients_vector)
            return _evaluate_fit()

        x0 = np.flatten([[layer_coefs[n].c[n], layer_coefs[n].s[n]] for layer_coefs in self.coefficients])
        results = optimize.root(objective, x0=x0)
        _set_coefficients(results.res)
        return None

    @summarize_values
    def groupwise_neutron_flux_in_layer(
        self, n: int, num_layer: int, x: float | npt.NDArray
    ) -> float | npt.NDArray:
        """
        Neutron flux[m^-2 s^-1] of the n-th group int h specified layer,
        at location x [m].

        Parameters
        ----------
        n:
            The index of the neutron group whose flux is being evaluated.
            n <= n_groups - 1.
            Therefore n=0 shows the flux for group 1, n=1 for group 2, etc.
        num_layer:
            The index of the layer that we want to get the neutron flux for.
        x:
            The position where the neutron flux has to be evaluated.

        Returns
        -------
        flux:
            Neutron flux at x meter from the first wall.
        """
        trig_funcs = []
        for g in range(n + 1):
            if self.l2[num_layer, g] > 0:
                c, s = np.cosh, np.sinh
                l = np.sqrt(self.l2[num_layer, g])
            else:
                c, s = np.cos, np.sin
                l = np.sqrt(-self.l2[num_layer, g])
            trig_funcs.append(self.coefficients[num_layer, n].c[g] * c(abs(x) / l))
            trig_funcs.append(self.coefficients[num_layer, n].s[g] * s(abs(x) / l))
        return np.sum(trig_funcs, axis=0)


    @summarize_values
    def groupwise_neutron_heating_in_layer(
        self, n: int, num_layer: int, x: float | npt.NDArray
    ) -> float | npt.NDArray:
        return (
            self.groupwise_linear_heating_density_in_layer(n, num_layer)
            * self.groupwise_neutron_flux_in_layer(n, num_layer, x)
        )

    @summarize_values
    def groupwise_neutron_flux_at(
        self, n: int, x: float | npt.NDArray
    ) -> float | npt.NDArray:
        """
        Neutron flux [m^-2 s^-1] anywhere. Neutron flux is assumed to be
        unperturbed once it leaves the final layer.

        Parameters
        ----------
        n:
            Neutron group index. n <= n_groups - 1.
            Therefore n=0 shows the flux for group 1, n=1 for group 2, etc.
        x:
            The depth where we want the neutron flux [m]. Out-of-bounds x would
            be clipped back to bound.
        """
        if np.isscalar(x):
            return self.groupwise_neutron_flux_at(n, [x])[0]
        x = np.clip(np.asarray(x), -self.layer_x[-1], self.layer_x[-1])

        out_flux = np.zeros_like(x)
        abs_x = abs(x)
        for num_layer, (xmin, xmax) in enumerate(pairwise(self.interface_x)):
            in_layer = xmin <= abs_x < xmax
            if num_layer==(self.n_layers-1):
                in_layer = xmin <= abs_x <= xmax
            if in_layer.any():
                out_flux[in_layer] = self.groupwise_neutron_flux_in_layer(
                    n, num_layer, x[in_layer]
                )
        return out_flux

    @summarize_values
    def groupwise_neutron_current_in_layer(
        self, n: int, num_layer: int, x: float | npt.NDArray
    ) -> float | npt.NDArray:
        """
        Get the neutron current (right=positive, left=negative) in any layer.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved. n <= n_groups - 1.
            Therefore n=0 shows the integrated flux for group 1, n=1 for group 2, etc.
        num_layer:
            The index of the layer that we want to get the neutron flux for.
        x:
            The depth where we want the neutron flux [m].
        """
        differentials = []
        for g in range(n + 1):
            l2g = self.l2[num_layer, g]
            if l2g > 0:
                l = np.sqrt(l2g)
                differentials.append(
                    self.coefficients[num_layer, n].c[g]
                    / l
                    * np.sinh(abs(x) / l)
                )
                differentials.append(
                    self.coefficients[num_layer, n].s[g]
                    / l
                    * np.cosh(abs(x) / l)
                )
            else:
                l = np.sqrt(-l2g)
                differentials.append(
                    -self.coefficients[num_layer, n].c[g]
                    / l
                    * np.sin(abs(x) / l)
                )
                differentials.append(
                    self.coefficients[num_layer, n].s[g]
                    / l
                    * np.cos(abs(x) / l)
                )
    
        return -self.diffusion_const[num_layer, n] * np.sum(differentials, axis=0) * _get_sign_of(x)

    @summarize_values
    def groupwise_neutron_current_at(
        self, n: int, x: float | npt.NDArray
    ) -> float | npt.NDArray:
        """
        Neutron current [m^-2 s^-1]. Neutron current is assumed to be
        unperturbed once it leaves the final layer.

        Parameters
        ----------
        n:
            Neutron group index. n <= n_groups - 1.
            Therefore n=0 shows the neutron current for group 1, n=1 for group 2, etc.
        x:
            The depth where we want the neutron current [m]. Out-of-bounds x
            would be clipped back to bound.
        """
        if np.isscalar(x):
            return self.groupwise_neutron_flux_at(n, [x])[0]
        x = np.clip(np.asarray(x), -self.layer_x[-1], self.layer_x[-1])

        current = np.zeros_like(x)
        abs_x = abs(x)
        for num_layer, (xmin, xmax) in enumerate(pairwise(self.interface_x)):
            in_layer = xmin <= abs_x < xmax
            if num_layer==(self.n_layers-1):
                in_layer = xmin <= abs_x <= xmax
            if in_layer.any():
                current[in_layer] = self.groupwise_neutron_current_in_layer(
                    n, num_layer, x[in_layer]
                )
        return current

    @summarize_values
    def groupwise_neutron_current_through_interface(
        self, n: int, n_interface: int, *, default_to_inner_layer : bool=True
    ) -> float:
        """
        Net current from left to right on the positive side of the model, at
        the specified interface number.

        Parameters
        ----------
        n:
            The index of the neutron group that we want the current for. n <= n_groups - 1.
            Therefore n=0 shows the neutron current for group 1, n=1 for group 2, etc.
        n_interface:
            The index of the interface that we want the net neutron current
            through for.
            E.g. for n_interface=0, that would be getting the current at the
            plasma-fw interface. For n_interface=n_layers, that would be
            getting the neutron current leaking from the final layer into the void.

        Returns
        -------
        :
            current in m^-2
        """
        x = self.interface_x[n_interface]

        if default_to_inner_layer:
            if n_interface==0:
                return self.groupwise_neutron_current_in_layer(n, 0, x)
            return self.groupwise_neutron_current_in_layer(n, n_interface-1, x)
        else:
            if n_interface==self.n_layers:
                return self.groupwise_neutron_current_in_layer(
                    n, self.n_layers-1, x
                )
            return self.groupwise_neutron_current_in_layer(n, n_interface, x)

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
        return self.groupwise_neutron_current_through_interface(
            n, self.n_layers
        )

    # scalar values (one such float per neutron group, and per layer.)
    @summarize_values
    def groupwise_integrated_flux_in_layer(
        self, n: int, num_layer: int
    ) -> float:
        """
        Calculate the integrated flux[m^-1 s^-1], which can be mulitplied to any
        macroscopic cross-section [m^-1] to get the reaction rate [s^-1] in
        any layer specified.

        Parameters
        ----------
        n:
            The index of the neutron group that needs to be solved. n <= n_groups - 1.
            Therefore n=0 shows the integrated flux for group 1, n=1 for group 2, etc.
        num_layer:
            The index of the layer that we want to get the neutron flux for.
        """
        integrals = []
        # set integration limits
        x_start = self.layer_x[num_layer-1]
        if num_layer==0:
            x_start = 0.0
        x_end = self.layer_x[num_layer]

        for g in range(n + 1):
            l2g = self.l2[num_layer, g]
            if l2g > 0:
                l = np.sqrt(l2g)
                integrals.append(
                    l * self.coefficients[num_layer, n].c[g]
                    * (np.sinh(x_end / l) - np.sinh(x_start / l))
                )
                integrals.append(
                    l * self.coefficients[num_layer, n].s[g]
                    * (np.cosh(x_end / l) - np.cosh(x_start / l))
                )
            else:
                l = np.sqrt(-l2g)
                integrals.append(
                    l * self.coefficients[num_layer, n].c[g]
                    * (np.sin(x_end / l) - np.sin(x_start / l))
                )
                integrals.append(
                    -l * self.coefficients[num_layer, n].s[g]
                    * (np.cos(x_end / l) - np.cos(x_start / l))
                )
        return np.sum(integrals, axis=0)

    @summarize_values
    def groupwise_integrated_heating_in_layer(
        self, n: int, num_layer: int,
    ) -> float:
        """
        How much neutron heating is 
        """
        return (
            self.groupwise_linear_heating_density_in_layer(n, num_layer)
            * self.groupwise_integrated_flux_in_layer(n, num_layer)
        )

    # Do NOT add a summarize_values decorator, as you can't add cross-sections
    # from different groups together without first multiplying by flux to get reaction rate.
    def groupwise_linear_heating_density_in_layer(
        self, n: int, num_layer: int,
    ) -> float:
        """
        All reactions that does not lead to scattering are assumed to have
        the full energy of the neutron deposited into the material.
        Obviously this contradicts the assumption of neutrons retaining some of
        its energy in the n,2n reaction, but we hope this is a small enough
        error that we can overlook it.
        """
        mat = self.materials[num_layer]
        non_scatter_xs = mat.sigma_t[n] - mat.sigma_s[n,:].sum()
        lost_energy = (
            (self.group_energy[n] - self.group_energy) * mat.sigma_s[n,:]
        ).sum()
        return self.group_energy[n] * non_scatter_xs + lost_energy

    @classmethod
    def get_output_unit(cls, method) -> str | None:
        """
        Check a method's outputted quantity's unit
        Parameters
        ----------
        method:
            A method whose name we shall be inspecting and comparing against
            UNIT_LOOKUP.
        Returns
        -------
        :
            If a match is found, return the unit as a string. Otherwise, return
            None.
        """
        for quantity, unit in UNIT_LOOKUP.items():
            if quantity in method.__name__:
                return unit

    def plot(
        self,
        quantity: str = "flux",
        ax: plt.Axes | None = None,
        *,
        plot_groups: bool = True,
        symmetric: bool = True,
        n_points: int = 100,
    ):
        """
        Make a rough plot of the neutron flux.

        Parameters
        ----------
        quantity:
            Options of plotting which quantity: {"flux", "current", "heating"}.
        ax:
            A plt.Axes object to plot on.
        n_points:
            Number of points to be used for plotting.
        symmetric:
            Whether to plot from -x to x (symmetric), or from 0 to x
            (right side only.)
        plot_groups:
            Whether to plot each individual group's neutron flux.
            If True, a legend will be added to help label the groups.
        """
        self.solve_group_n(self.n_groups - 1)
        ax = ax or plt.axes()
        method_name = f"neutron_{quantity}_in_layer"
        total_function = getattr(self, method_name)
        unit = self.get_output_unit(total_function)
        ylabel = f"{quantity}({unit})"

        if plot_groups:
            groupwise_function = getattr(self, f"groupwise_{method_name}")
        x_ranges = _generate_x_range(
            self.interface_x.copy(),
            max(self.extended_boundary.values()),
            min_total_num_points=n_points,
            symmetric=symmetric
        )
        for num_layer in range(self.n_layers):
            if symmetric:
                neg_x = next(x_ranges)
                ax.plot(neg_x, total_function(num_layer, neg_x), color="black")
            pos_x = next(x_ranges)
            plot_dict = {"label": "total"} if num_layer==0 else {}
            ax.plot(
                pos_x, total_function(num_layer, pos_x),
                color="black", **plot_dict,
            )
            if plot_groups:
                for n in range(self.n_groups):
                    if symmetric:
                        ax.plot(
                            neg_x, groupwise_function(n, num_layer, neg_x),
                            color=f"C{n}"
                        )
                    plot_dict = {"label": f"group {n}"} if num_layer==0 else {}
                    ax.plot(
                        pos_x, groupwise_function(n, num_layer, pos_x),
                        color=f"C{n}", **plot_dict,
                    )
        ax.legend()
        ax.set_title(f"Neutron {quantity} profile")
        ax.set_xlabel("Distance from the plasma-fw interface [m]")
        ax.set_ylabel(ylabel)
        
        # plotting the interfaces for ease of comprehension.
        ylims = ax.get_ylim()
        for (xmin, xmax), mat in zip(pairwise(self.interface_x), self.materials):
            _plot_vertical_dotted_line(ax, xmin, ylims, symmetric=symmetric)
            ax.text(np.mean([xmin, xmax]), 0, mat.name, ha="center", va="center")
            if symmetric:
                ax.text(-np.mean([xmin, xmax]), 0, mat.name, ha="center", va="center")
        _plot_vertical_dotted_line(ax, xmax, ylims, symmetric=symmetric)
        return ax

def _get_sign_of(x_values):
    """
    Get sign of any real number, but also forces 0.0 to be +ve and -0.0 to be -ve.
    The neutron current for the first group (in a non-breeding/weakly breeding
    scenario) is strongest at x=0, but have different signs when limit x-> 0^+
    and limit x-> 0^-. This function allows the input x to behave like 0^+ when
    it's =0.0 and like 0^- when it's =-0.0, giving the correct neutron current
    at those locations, rather than setting the neutron current to zer0.
    """
    negatives = np.signbit(x_values)
    return np.array(negatives, dtype=float) * -2 + 1

def _plot_vertical_dotted_line(ax, x, ylims, *, symmetric: bool=True):
    if symmetric:
        ax.plot([-x, -x], ylims, color="black", ls="--")
    ax.plot([x, x], ylims, color="black", ls="--")
    return None

def _generate_x_range(
    interface_x: npt.NDArray[np.float64],
    extended_boundary: float | None,
    *,
    min_total_num_points: int,
    symmetric: bool,
):
    """Helper generator for finding the range of x-values to be plotted.

    Parameters
    ----------
    interface_x:
        The x-coordinates of each interface. It should be all-positive, and
        ascending only.
    extended_boundary:
        extended boundary
    min_total_num_points:
        Approximate number of points to be plotted. Increase this number to increase the resolution
    symmetric:
        Whether to return two copies (one negative, one positive) of x-ranges
        per layer.

    Yields
    -------
    :
        A generator of x-coordinates used for plotting, each of which falls
        within the limit of the xmin and xmax of that layer, forming a total
        of a minimum of min_total_num_points.
        If symmetric=True, then the number of numpy arrays in the list
        = 2*n_layers, where out_x_range[0] and out_x_range[1] are for the
        negative and positive sides of the first layer respectively;
        out_x_range[2] and out_x_range[3] are for the negative and positive
        sides of the second layer respectively.
        Otherwise, the number of numpy arrays in the list = n_layers.
    """
    full_x_range = np.linspace(
        interface_x.min(), interface_x.max(), min_total_num_points
    )
    for num_layer, (xmin, xmax) in enumerate(pairwise(interface_x)):
        num_points = np.logical_and(xmin<full_x_range, full_x_range<xmax).sum() + 2
        layer_x_range = np.linspace(xmin, xmax, num_points)
        if num_layer==(len(interface_x)-1):
            layer_x_range = np.append(layer_x_range, abs(extended_boundary))
        if symmetric:
            yield -layer_x_range[::-1]
        yield layer_x_range.copy()
