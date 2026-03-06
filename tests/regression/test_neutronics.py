import numpy as np
from scipy.integrate import trapezoid

from process.neutronics import NeutronFluxProfile
from process.neutronics_data import DT_NEUTRON_E, EV_TO_J, MaterialMacroInfo, scattering_weight_matrix

MAX_E = DT_NEUTRON_E * 1.01
MIN_E = 1/40 * EV_TO_J

def _diffusion_equation_in_layer(test_profile, n, num_layer, x):
    """
    Get the three terms in the diffusion equation (equation 5 in the paper.
    """
    diffusion_out = test_profile.diffusion_const[num_layer, n] * test_profile._groupwise_flux_curvature_in_layer(n, num_layer, x)
    total_removal = test_profile.materials[num_layer].sigma_t[n] * test_profile.groupwise_neutron_flux_in_layer(n, num_layer, x)

    source_in_terms = []
    in_matrix = test_profile.materials[num_layer].sigma_s + test_profile.materials[num_layer].sigma_in
    for g, all_sources_entering_from_g in enumerate(in_matrix[:, n]):
        source_in_terms.append(all_sources_entering_from_g * test_profile.groupwise_neutron_flux_in_layer(g, num_layer, x))

    return diffusion_out, total_removal, np.sum(source_in_terms)


def test_one_group():
    """
    Regression test against Desmos snapshot:
    https://www.desmos.com/calculator/18xojespuo
    """
    dummy = [MAX_E, MIN_E]  # dummy group structure
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
        incoming_flux, [x_fw, x_bz], [fw_material, bz_material]
    )

    layer_group_coefs = neutron_profile.coefficients
    assert np.isclose(layer_group_coefs[0, 0].c[0], 80.5346770788), "c5"
    assert np.isclose(layer_group_coefs[0, 0].s[0], -76.5562120985), "c6"
    assert np.isclose(layer_group_coefs[1, 0].c[0], 60.6871656017), "c7"
    assert np.isclose(layer_group_coefs[1, 0].s[0], -60.7123696772), "c8"

    assert np.isclose(neutron_profile.neutron_flux_in_layer(0, x_fw), 48.72444)
    assert np.isclose(neutron_profile.neutron_flux_in_layer(1, x_fw), 48.72444)
    assert np.isclose(neutron_profile.neutron_flux_at(x_fw), 48.72444)

    assert np.isclose(
        neutron_profile.neutron_current_through_interface(1), 22.3980214162
    )
    assert np.isclose(neutron_profile.neutron_current_escaped(), 1.22047369356)

    fw_removal = sigma_fw_t - sigma_fw_s - fw_material.sigma_in[0, 0]
    bz_removal = sigma_bz_t - sigma_bz_s - bz_material.sigma_in[0, 0]
    assert np.isclose(
        neutron_profile.fluxes[0],
        neutron_profile.neutron_current_escaped()
        + fw_removal * neutron_profile.integrated_flux_in_layer(0)
        + bz_removal * neutron_profile.integrated_flux_in_layer(1),
    ), "Conservation of neutrons"

    x_fw = np.linspace(*neutron_profile.interface_x[0:2], 100000)
    manually_integrated_heating_fw = trapezoid(
        neutron_profile.neutron_heating_in_layer(0, x_fw),
        x_fw,
    )
    x_bz = np.linspace(*neutron_profile.interface_x[1:3], 100000)
    manually_integrated_heating_bz = trapezoid(
        neutron_profile.neutron_heating_in_layer(1, x_bz),
        x_bz,
    )
    assert np.isclose(
        neutron_profile.integrated_heating_in_layer(0),
        manually_integrated_heating_fw,
        atol=0,
        rtol=1e-8,
    ), "Correctly integrated heating in FW"
    assert np.isclose(
        neutron_profile.integrated_heating_in_layer(1),
        manually_integrated_heating_bz,
        atol=0,
        rtol=1e-8,
    ), "Correctly integrated heating in BZ"


def test_one_group_with_fission():
    """
    Regression test against Desmos snapshot with fission involved:
    https://www.desmos.com/calculator/cd830add9c
    Expecting a cosine-shape (dome shape!) of neutron flux profile.
    """
    dummy = [MAX_E, MIN_E]
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
        incoming_flux,
        [x_fw, x_bz],
        [fw_material, bz_material],
    )
    assert np.isclose(neutron_profile.l2[1, 0], -((58.2869567709 / 100) ** 2))
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
        neutron_profile.neutron_flux_in_layer(0, x_fw),
        neutron_profile.neutron_flux_in_layer(1, x_fw),
    ), "Flux continuity assurance"
    assert np.isclose(
        neutron_profile.neutron_current_through_interface(1),
        -7.6275782637960745,
    ), "Negative current because BZ (breeding) is backflowing into the FW"
    assert np.isclose(
        neutron_profile.neutron_current_escaped(), 30.665951670177186
    ), "positive escaped current."
    fw_removal = sigma_fw_t - sigma_fw_s - fw_material.sigma_in[0, 0]
    bz_removal = sigma_bz_t - sigma_bz_s - bz_material.sigma_in[0, 0]

    assert np.isclose(
        neutron_profile.fluxes[0],
        neutron_profile.neutron_current_escaped()
        + fw_removal * neutron_profile.integrated_flux_in_layer(0)
        + bz_removal * neutron_profile.integrated_flux_in_layer(1),
    ), "Conservation of neutrons"
    x_fw = np.linspace(*neutron_profile.interface_x[0:2], 100000)
    manually_integrated_heating_fw = trapezoid(
        neutron_profile.neutron_heating_in_layer(0, x_fw),
        x_fw,
    )
    x_bz = np.linspace(*neutron_profile.interface_x[1:3], 100000)
    manually_integrated_heating_bz = trapezoid(
        neutron_profile.neutron_heating_in_layer(1, x_bz),
        x_bz,
    )
    assert np.isclose(
        neutron_profile.integrated_heating_in_layer(0),
        manually_integrated_heating_fw,
        atol=0,
        rtol=1e-8,
    ), "Correctly integrated heating in FW"
    assert np.isclose(
        neutron_profile.integrated_heating_in_layer(1),
        manually_integrated_heating_bz,
        atol=0,
        rtol=1e-8,
    ), "Correctly integrated heating in BZ"

def test_two_identical_materials():
    """
    A 2-layer model (both layers being made of material A) should have the same
    neutron spectrum and flux profiles as a one-layer model.
    """
    
def test_5_5():
    """Create an arbitrary 5-layer 5-group model. Check for continuity and conformity to the equation."""
    dummy_group_structure = np.geomspace(MAX_E, MIN_E, 5+1)
    mat_list = []
    at_masses = np.geomspace(1, 100, 5)[[3,0,2,4,1]]
    sigma_t_lists = [  # arbitrarily chosen and rearranged numbers
        1/(80 + np.linspace(-40, 40, 5)[[0,2,4,1,3]]),
        1/(200 + np.linspace(-30, 30, 5)[[0,3,1,4,2]]),
        1/(300 + np.linspace(-20, 20, 5)[[1,3,0,2,4]]),
        1/(100 + np.linspace(-40, 40, 5)[[1,4,2,0,3]]),
        1/(50 + np.linspace(-10, 10, 5)[[2,0,3,1,4]]),
    ]
    sigma_s_list = [
        sigma_t_lists[0] * [0.8, 0.7, 0.6, 0.6, 0.5],
        sigma_t_lists[1] * [0.4, 0.3, 0.2, 0.1, 0.0],
        sigma_t_lists[2] * [0.8, 0.7, 0.6, 0.6, 0.5],
        sigma_t_lists[3] * [0.9, 0.9, 0.6, 0.4, 0.1],
        sigma_t_lists[4] * [0.6, 0.5, 0.3, 0.3, 0.2],
    ]
    sigma_in_list = [
        [0.001, 0.002, 0, 0, 0.005],
        np.where([0,0,1,1,0], sigma_t_lists[1]*1.0, [0,0,0,0,0]),
        np.where([0,1,0,1,0], sigma_t_lists[2]*0.6, [0,0,0,0,0]),
        np.where([1,0,1,1,1], sigma_t_lists[3]*2.2, [0,0,0,0,0]),
        [0,0,0,0,0],
    ]
    for i in range(5):
        mat_list.append(MaterialMacroInfo(
            dummy_group_structure, at_masses[i],
            sigma_t=sigma_t_lists[i],
            sigma_s=(sigma_s_list[i] * scattering_weight_matrix(
                dummy_group_structure, at_masses[i]
            ).T).T,
            sigma_in=(
                sigma_in_list[i] * scattering_weight_matrix(
                dummy_group_structure, at_masses[i]
            ).T).T,
            name=f"mat{i}"
        ))
    neutron_profile = NeutronFluxProfile(
        1.0, [5, 10, 15, 20, 25], mat_list
    )
    neutron_profile.solve_group_n(4)
