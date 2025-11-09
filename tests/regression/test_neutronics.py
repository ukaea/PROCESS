import pytest
import numpy as np

from process.neutronics import NeutronFluxProfile
from process.neutronics_data import MaterialMacroInfo

def test_against_desmos_number():
    """
    Regression test against Desmos snapshot:
    https://www.desmos.com/calculator/18xojespuo
    """
    dummy = [100, 1]  # dummy group structure
    # translate from mean-free-path lengths (mfp) to macroscopic cross-sections
    # [m]
    mfp_fw_s = 118 * 0.01
    mfp_fw_t = 16.65 * 0.01
    sigma_fw_t = 1 / mfp_fw_t  # [1/m]
    sigma_fw_s = 1 / mfp_fw_s  # [1/m]
    a_fw = 52
    fw_material = MaterialMacroInfo([sigma_fw_t], [[sigma_fw_s]], dummy, a_fw)

    mfp_bz_s = 97 * 0.01
    mfp_bz_t = 35.8 * 0.01
    sigma_bz_s = 1 / mfp_bz_s  # [1/m]
    sigma_bz_t = 1 / mfp_bz_t  # [1/m]
    a_bz = 71
    bz_material = MaterialMacroInfo([sigma_bz_t], [[sigma_bz_s]], dummy, a_bz)

    x_fw, x_bz = 5.72 * 0.01, 85 * 0.01
    incoming_flux = 41
    neutron_profile = NeutronFluxProfile(
        incoming_flux, x_fw, x_bz, fw_material, bz_material
    )
    
    const = neutron_profile.integration_constants[0]  # alias to fit line width
    assert np.isclose(const.fw_pos, 1.98923249017), "c1"
    assert np.isclose(const.fw_neg, 78.5454445887), "c2"
    assert np.isclose(const.bz_pos, -0.0126020377605), "c3"
    assert np.isclose(const.bz_neg, 60.6997676395), "c4"

    assert np.isclose(neutron_profile.neutron_flux_fw(x_fw), 48.72444)
    assert np.isclose(neutron_profile.neutron_flux_bz(x_fw), 48.72444)
    assert np.isclose(neutron_profile.neutron_flux_at(x_fw), 48.72444)

    assert np.isclose(neutron_profile.neutron_current_fw2bz(), 22.3980214162)
    assert np.isclose(neutron_profile.neutron_current_escaped(), 1.22047369356)

    assert np.isclose(neutron_profile.flux,
        neutron_profile.neutron_current_escaped()
        + neutron_profile.reaction_rate_fw("removal")
        + neutron_profile.reaction_rate_bz("removal")
    )
def test_one_group_with_fission():
    """Expecting a cosine-shape (dome shape!) of neutron flux profile"""
