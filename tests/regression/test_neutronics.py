import numpy as np

from process.neutronics import NeutronFluxProfile
from process.neutronics_data import MaterialMacroInfo


def test_one_group():
    """
    Regression test against Desmos snapshot:
    https://www.desmos.com/calculator/18xojespuo
    """
    dummy = [100, 1]  # dummy group structure
    # translate from mean-free-path lengths (mfp) to macroscopic cross-sections
    mfp_fw_s = 118 * 0.01  # [m]
    mfp_fw_t = 16.65 * 0.01  # [m]
    sigma_fw_t = 1 / mfp_fw_t  # [1/m]
    sigma_fw_s = 1 / mfp_fw_s  # [1/m]
    a_fw = 52
    fw_material = MaterialMacroInfo(
        dummy, a_fw, [sigma_fw_t], [[sigma_fw_s]], name="fw"
    )

    mfp_bz_s = 97 * 0.01  # [m]
    mfp_bz_t = 35.8 * 0.01  # [m]
    sigma_bz_s = 1 / mfp_bz_s  # [1/m]
    sigma_bz_t = 1 / mfp_bz_t  # [1/m]
    a_bz = 71
    bz_material = MaterialMacroInfo(
        dummy, a_bz, [sigma_bz_t], [[sigma_bz_s]], name="bz"
    )

    x_fw, x_bz = 5.72 * 0.01, 85 * 0.01
    incoming_flux = 41
    neutron_profile = NeutronFluxProfile(
        incoming_flux, x_fw, x_bz, fw_material, bz_material
    )

    const = neutron_profile.coefficients[0]  # alias to fit line width
    assert np.isclose(const.fw_c[0], 80.5346770788), "c5"
    assert np.isclose(const.fw_s[0], -76.5562120985), "c6"
    assert np.isclose(const.bz_c[0], 60.6871656017), "c7"
    assert np.isclose(const.bz_s[0], -60.7123696772), "c8"

    assert np.isclose(neutron_profile.neutron_flux_fw(x_fw), 48.72444)
    assert np.isclose(neutron_profile.neutron_flux_bz(x_fw), 48.72444)
    assert np.isclose(neutron_profile.neutron_flux_at(x_fw), 48.72444)

    assert np.isclose(neutron_profile.neutron_current_fw2bz(), 22.3980214162)
    assert np.isclose(neutron_profile.neutron_current_escaped(), 1.22047369356)

    fw_removal = sigma_fw_t - sigma_fw_s - fw_material.sigma_in[0, 0]
    bz_removal = sigma_bz_t - sigma_bz_s - bz_material.sigma_in[0, 0]
    assert np.isclose(
        neutron_profile.flux,
        neutron_profile.neutron_current_escaped()
        + fw_removal * neutron_profile.integrated_flux_fw()
        + bz_removal * neutron_profile.integrated_flux_bz(),
    ), "Conservation of neutrons"


def test_one_group_with_fission():
    """
    Regression test against Desmos snapshot with fission involved:
    https://www.desmos.com/calculator/cd830add9c
    Expecting a cosine-shape (dome shape!) of neutron flux profile.
    """
    dummy = [100, 1]
    mfp_fw_s = 118 * 0.01  # [m]
    mfp_fw_t = 16.65 * 0.01  # [m]
    sigma_fw_t = 1 / mfp_fw_t  # [1/m]
    sigma_fw_s = 1 / mfp_fw_s  # [1/m]
    a_fw = 52
    fw_material = MaterialMacroInfo(
        dummy, a_fw, [sigma_fw_t], [[sigma_fw_s]], name="fw"
    )

    mfp_bz_s = 97 * 0.01  # [m]
    mfp_bz_t = 35.8 * 0.01  # [m]
    sigma_bz_s = 1 / mfp_bz_s  # [1/m]
    sigma_bz_t = 1 / mfp_bz_t  # [1/m]
    a_bz = 71

    g = 1.2
    nu_sigma_bz_f = g * (sigma_bz_t - sigma_bz_s)
    bz_material = MaterialMacroInfo(
        dummy, a_bz, [sigma_bz_t], [[sigma_bz_s]], [[nu_sigma_bz_f]], name="bz"
    )

    x_fw, x_bz = 5.72 * 0.01, 85 * 0.01
    incoming_flux = 41
    neutron_profile = NeutronFluxProfile(
        incoming_flux, x_fw, x_bz, fw_material, bz_material
    )
    assert np.isclose(neutron_profile.l_bz_2[0], -((58.2869567709 / 100) ** 2))
    assert np.isclose(
        neutron_profile.neutron_flux_at(-4.79675 / 100), 159.9434
    ), "Minimum flux in FW"
    assert np.isclose(
        neutron_profile.neutron_flux_at(4.79675 / 100), 159.9434
    ), "Minimum flux in FW"
    assert np.isclose(
        neutron_profile.neutron_flux_at(18.96382 / 100), 164.81245
    ), "Maximum flux in BZ"
    assert np.isclose(
        neutron_profile.neutron_flux_at(-18.96382 / 100), 164.81245
    ), "Maximum flux in BZ"
    assert np.isclose(
        neutron_profile.neutron_flux_fw(x_fw),
        neutron_profile.neutron_flux_bz(x_fw),
    ), "Flux continuity assurance"
    assert np.isclose(
        neutron_profile.neutron_current_fw2bz(), -7.6275782637960745
    ), "Negative current because BZ (breeding) is backflowing into the FW"
    assert np.isclose(
        neutron_profile.neutron_current_escaped(), 30.665951670177186
    ), "positive escaped current."
    fw_removal = sigma_fw_t - sigma_fw_s - fw_material.sigma_in[0, 0]
    bz_removal = sigma_bz_t - sigma_bz_s - bz_material.sigma_in[0, 0]
    assert np.isclose(
        neutron_profile.flux,
        neutron_profile.neutron_current_escaped()
        + fw_removal * neutron_profile.integrated_flux_fw()
        + bz_removal * neutron_profile.integrated_flux_bz(),
    ), "Conservation of neutrons"


def test_two_groups():
    """
    Same group n=0 values as test_one_group. Second group can have a massive removal cross-section
    """
