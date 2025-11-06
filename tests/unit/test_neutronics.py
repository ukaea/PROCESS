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


def test_against_desmos_number():
    dummy = [100, 1]  # dummy group structure
    # translate from mean-free-path lengths (mfp) to macroscopic cross-sections
    # [cm]
    mfp_fw_s = 118 * 0.01
    mfp_fw_t = 16.65 * 0.01
    sigma_fw_t = 1 / mfp_fw_t  # [1/cm]
    sigma_fw_s = 1 / mfp_fw_s  # [1/cm]
    a_fw = 52
    fw_material = MaterialMacroInfo([sigma_fw_t], [[sigma_fw_s]], dummy, a_fw)

    mfp_bz_s = 97 * 0.01
    mfp_bz_t = 35.8 * 0.01
    sigma_bz_s = 1 / mfp_bz_s  # [1/cm]
    sigma_bz_t = 1 / mfp_bz_t  # [1/cm]
    a_bz = 71
    bz_material = MaterialMacroInfo([sigma_bz_t], [[sigma_bz_s]], dummy, a_bz)

    x_fw, x_bz = 5.72 * 0.01, 85 * 0.01
    incoming_flux = 41
    neutron_profile = NeutronFluxProfile(
        incoming_flux, x_fw, x_bz, fw_material, bz_material
    )
    neutron_profile.plot()
