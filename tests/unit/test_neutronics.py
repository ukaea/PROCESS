import numpy as np
import pytest

from process.exceptions import ProcessValidationError
from process.neutronics import (
    AutoPopulatingDict,
    LayerSpecificGroupwiseConstants,
    NeutronFluxProfile,
    _get_sign_of,
)
from process.neutronics_data import DT_NEUTRON_E, EV_TO_J, MaterialMacroInfo

MAX_E = DT_NEUTRON_E * 1.01
MIN_E = 1e-9 * EV_TO_J


def test_group_structure_0_energy():
    with pytest.warns():
        MaterialMacroInfo(
            [1.0, 0.0],
            0.1,
            [1.0],
            [[1.0]],
        )


def test_group_structure_too_short():
    with pytest.raises(ProcessValidationError):
        MaterialMacroInfo(
            [1.0],
            0.1,
            [1.0],
            [[1.0]],
        )


def test_sigma_s_incorrect_shape():
    with pytest.raises(ProcessValidationError):
        MaterialMacroInfo(
            [1000, 10, 1.0],
            0.1,
            [1.0, 2.0],
            [1.0, 1.0],
        )


def test_sigma_s_too_large():
    with pytest.raises(ProcessValidationError):
        MaterialMacroInfo(
            [1000, 10, 1.0],
            0.1,
            [1.0, 2.0],
            [[1.0, 1.0], [1.0, 1.0]],
        )


def test_warn_up_elastic_scatter():
    with pytest.warns():
        MaterialMacroInfo(
            [1000, 10, 1.0],
            0.1,
            [1.0, 2.0],
            [[0.5, 0.5], [1.0, 1.0]],
        )


def test_throw_index_error():
    layer_specific_const = LayerSpecificGroupwiseConstants(
        lambda x: x, ["", ""], ["Dummy constants"]
    )
    with pytest.raises(IndexError):
        layer_specific_const[0, 0, 0]
    with pytest.raises(IndexError):
        layer_specific_const[0, 0, 0] = 1


def test_iter_and_len():
    layer_specific_const = LayerSpecificGroupwiseConstants(
        lambda x: x, ["", ""], ["Dummy constants"]
    )
    assert len(layer_specific_const) == 2
    as_list = list(layer_specific_const)
    assert len(as_list) == 2
    assert isinstance(as_list[0], AutoPopulatingDict)


def test_has_local_fluxes():
    """Test that the groupwise decorator has worked on the local fluxes methods."""
    assert hasattr(NeutronFluxProfile, "neutron_flux_at")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_flux_at")
    assert hasattr(NeutronFluxProfile, "neutron_flux_in_layer")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_flux_in_layer")


def test_has_boundary_current():
    """Test that the groupwise decorator has worked on the boundary fluxes methods."""
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_current_in_layer")
    assert hasattr(NeutronFluxProfile, "neutron_current_in_layer")
    assert hasattr(
        NeutronFluxProfile, "groupwise_neutron_current_through_interface"
    )
    assert hasattr(NeutronFluxProfile, "neutron_current_through_interface")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_current_escaped")
    assert hasattr(NeutronFluxProfile, "neutron_current_escaped")


def test_has_reaction_rates():
    """Test that the groupwise decorator has worked on the reactions methods."""
    assert hasattr(NeutronFluxProfile, "groupwise_integrated_flux_in_layer")
    assert hasattr(NeutronFluxProfile, "integrated_flux_in_layer")
    assert hasattr(NeutronFluxProfile, "groupwise_integrated_heating_in_layer")
    assert hasattr(NeutronFluxProfile, "integrated_heating_in_layer")


def test_has_volumetric_heating():
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_heating_in_layer")
    assert hasattr(NeutronFluxProfile, "neutron_heating_in_layer")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_heating_at")
    assert hasattr(NeutronFluxProfile, "neutron_heating_at")


def test_has_plot():
    assert hasattr(NeutronFluxProfile, "plot")


def test_units():
    nfp = NeutronFluxProfile
    assert (
        nfp.get_output_unit(nfp.groupwise_integrated_flux_in_layer)
        == "m^-1 s^-1"
    )
    assert (
        nfp.get_output_unit(nfp.groupwise_integrated_heating_in_layer)
        == "W m^-2"
    )
    assert (
        nfp.get_output_unit(nfp.groupwise_linear_heating_density_in_layer)
        == "J m^-1"
    )
    assert nfp.get_output_unit(nfp.groupwise_neutron_current_at) == "m^-2 s^-1"
    assert (
        nfp.get_output_unit(nfp.groupwise_neutron_current_escaped)
        == "m^-2 s^-1"
    )
    assert (
        nfp.get_output_unit(nfp.groupwise_neutron_current_in_layer)
        == "m^-2 s^-1"
    )
    assert (
        nfp.get_output_unit(nfp.groupwise_neutron_current_through_interface)
        == "m^-2 s^-1"
    )
    assert nfp.get_output_unit(nfp.groupwise_neutron_flux_at) == "m^-2 s^-1"
    assert (
        nfp.get_output_unit(nfp.groupwise_neutron_flux_in_layer) == "m^-2 s^-1"
    )
    assert (
        nfp.get_output_unit(nfp.groupwise_neutron_heating_in_layer) == "W m^-3"
    )
    assert nfp.get_output_unit(nfp.groupwise_neutron_heating_at) == "W m^-3"

    assert nfp.get_output_unit(nfp.integrated_flux_in_layer) == "m^-1 s^-1"
    assert nfp.get_output_unit(nfp.integrated_heating_in_layer) == "W m^-2"
    assert nfp.get_output_unit(nfp.neutron_current_at) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_current_escaped) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_current_in_layer) == "m^-2 s^-1"
    assert (
        nfp.get_output_unit(nfp.neutron_current_through_interface)
        == "m^-2 s^-1"
    )
    assert nfp.get_output_unit(nfp.neutron_flux_at) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_flux_in_layer) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_heating_in_layer) == "W m^-3"
    assert nfp.get_output_unit(nfp.neutron_heating_at) == "W m^-3"


def test_get_sign_func():
    signs = _get_sign_of(
        np.array([-np.inf, -2, -1, -0.0, 0.0, 1.0, 2.0, np.inf])
    )
    np.testing.assert_equal(signs, [-1, -1, -1, -1, 1, 1, 1, 1])


def test_three_group():
    """Can use degenerate cross-sections values between groups?"""


def test_two_group():
    """Ensure continuity at interface for both groups."""

def test_two_identical_materials():
    """
    A 2-layer model (both layers being made of material A) should have the same
    neutron spectrum and flux profiles as a one-layer model.
    """