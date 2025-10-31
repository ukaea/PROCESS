import pytest

from process.exceptions import ProcessValidationError
from process.neutronics import NeutronFluxProfile
from process.neutronics_data import MaterialMacroInfo


def test_group_structure_too_short():
    with pytest.raises(ProcessValidationError):
        MaterialMacroInfo(
            [1.0],
            [[1.0]],
            [0.0],
            0.1,
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
    assert hasattr(NeutronFluxProfile, "groupwise_reaction_rate_fw")
    assert hasattr(NeutronFluxProfile, "reaction_rate_fw")
    assert hasattr(NeutronFluxProfile, "groupwise_reaction_rate_bz")
    assert hasattr(NeutronFluxProfile, "reaction_rate_bz")


# 2 groups
# fw_material = MaterialMacroInfo(
    # [1.0, 0.5], [[0.5, 0.1], [0, 0.1]], [1000, 100, 0], 1.0
# )
# bz_material = MaterialMacroInfo(
    # [1.0, 0.5], [[0.5, 0.1], [0, 0.1]], [1000, 100, 0], 2.0
# )
# 1 group
fw_material = MaterialMacroInfo(
    [1.0], [[0.5]], [100, 0], 1.0
)
bz_material = MaterialMacroInfo(
    [1.0], [[0.5]], [100, 0], 2.0
)
profile = NeutronFluxProfile(
    1.0,  # flux
    0.01,  # 1cm
    0.11,  # 10cm thick blanket
    fw_material,
    bz_material,
)
