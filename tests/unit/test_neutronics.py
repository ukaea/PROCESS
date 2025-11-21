import pytest

from process.exceptions import ProcessValidationError
from process.neutronics import NeutronFluxProfile
from process.neutronics_data import MaterialMacroInfo


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


def test_has_local_fluxes():
    """Test that the groupwise decorator has worked on the local fluxes methods."""
    assert hasattr(NeutronFluxProfile, "neutron_flux_at")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_flux_at")
    assert hasattr(NeutronFluxProfile, "neutron_flux_fw")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_flux_fw")
    assert hasattr(NeutronFluxProfile, "neutron_flux_bz")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_flux_bz")


def test_has_boundary_current():
    """Test that the groupwise decorator has worked on the boundary fluxes methods."""
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_current_fw2bz")
    assert hasattr(NeutronFluxProfile, "neutron_current_fw2bz")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_current_escaped")
    assert hasattr(NeutronFluxProfile, "neutron_current_escaped")


def test_has_reactions():
    """Test that the groupwise decorator has worked on the reactions methods."""
    assert hasattr(NeutronFluxProfile, "groupwise_integrated_flux_fw")
    assert hasattr(NeutronFluxProfile, "integrated_flux_fw")
    assert hasattr(NeutronFluxProfile, "groupwise_integrated_flux_bz")
    assert hasattr(NeutronFluxProfile, "integrated_flux_bz")


def test_three_group():
    # 1. No negative flux
    dummy = [10000, 1000, 100, 1]
    # fw_mat = MaterialMacroInfo()
    bz_mat = MaterialMacroInfo
    # 2. same L_1 and L_3 shoudl yield coefficients[2].fw_c[0] and
    # coefficients[2].fw_s[0] = 0.0
    # self.l_fw_2 == self.l_fw_2[n] case.
    #


def test_two_group():
    """Ensure continuity at interface for both groups."""
