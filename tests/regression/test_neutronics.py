import numpy as np
import pytest
from scipy.integrate import trapezoid

from process.models.neutronics.base import NeutronFluxProfile
from process.models.neutronics.data import (
    DT_NEUTRON_E,
    EV_TO_J,
    MaterialMacroInfo,
    calculate_mean_energy_and_incident_bin,
    scattering_weight_matrix,
)

MAX_E = DT_NEUTRON_E * 1.01
MIN_E = 1 / 40 * EV_TO_J


def _diffusion_equation_in_layer(test_profile, n, num_layer, x):
    """
    Get the three terms in the diffusion equation (equation 5 in the paper.
    """
    diffusion_out = test_profile.materials[num_layer].diffusion_const[
        n
    ] * test_profile._groupwise_flux_curvature_in_layer(n, num_layer, x)  # noqa: SLF001
    total_removal = test_profile.materials[num_layer].sigma_t[
        n
    ] * test_profile.groupwise_neutron_flux_in_layer(n, num_layer, x)

    source_in_terms = []
    in_matrix = (
        test_profile.materials[num_layer].sigma_s
        + test_profile.materials[num_layer].sigma_in
    )
    for g, all_sources_entering_from_g in enumerate(in_matrix[:, n]):
        source_in_terms.append(
            all_sources_entering_from_g
            * test_profile.groupwise_neutron_flux_in_layer(g, num_layer, x)
        )

    return diffusion_out, total_removal, np.sum(source_in_terms)


def test_1_group_1_layer():
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
    fw_material = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="fw")
    fw_material.avg_atomic_mass = a_fw
    fw_material._set_sigma([sigma_fw_t], [[sigma_fw_s]])  # noqa: SLF001

    mfp_bz_s = 97 * 0.01  # [m]
    mfp_bz_t = 35.8 * 0.01  # [m]
    sigma_bz_s = 1 / mfp_bz_s  # [1/m]
    sigma_bz_t = 1 / mfp_bz_t  # [1/m]
    a_bz = 71
    bz_material = MaterialMacroInfo(dummy, 1.0, {"Lu": 1.0}, name="bz")
    bz_material.avg_atomic_mass = a_bz
    bz_material._set_sigma([sigma_bz_t], [[sigma_bz_s]])  # noqa: SLF001

    x_fw, x_bz = 5.72 * 0.01, 85 * 0.01
    incoming_flux = 41
    neutron_profile = NeutronFluxProfile(
        incoming_flux, [x_fw, x_bz], [fw_material, bz_material]
    )

    layer_group_coefs = neutron_profile.coefficients
    assert np.isclose(layer_group_coefs[0, 0].c[0], 78.5454445887), "c1"
    assert np.isclose(layer_group_coefs[0, 0].s[0], 1.98923249017), "c2"
    assert np.isclose(layer_group_coefs[1, 0].c[0], 60.6997676395), "c3"
    assert np.isclose(layer_group_coefs[1, 0].s[0], -0.0126020377605), "c4"

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
    assert np.isclose(neutron_profile.neutron_current_at(0), incoming_flux)


def test_1_group_1_layer_with_fission():
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
    fw_material = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="fw")
    fw_material.avg_atomic_mass = a_fw
    fw_material._set_sigma([sigma_fw_t], [[sigma_fw_s]])  # noqa: SLF001

    mfp_bz_s = 97 * 0.01  # [m]
    mfp_bz_t = 35.8 * 0.01  # [m]
    sigma_bz_s = 1 / mfp_bz_s  # [1/m]
    sigma_bz_t = 1 / mfp_bz_t  # [1/m]
    a_bz = 71

    g = 1.2
    nu_sigma_bz_f = g * (sigma_bz_t - sigma_bz_s)
    bz_material = MaterialMacroInfo(dummy, 1.0, {"Lu": 1.0}, name="bz")
    bz_material.avg_atomic_mass = a_bz
    bz_material._set_sigma([sigma_bz_t], [[sigma_bz_s]], [[nu_sigma_bz_f]])  # noqa: SLF001

    x_fw, x_bz = 5.72 * 0.01, 85 * 0.01
    incoming_flux = 41
    neutron_profile = NeutronFluxProfile(
        incoming_flux,
        [x_fw, x_bz],
        [fw_material, bz_material],
    )
    assert np.isclose(neutron_profile.materials[1].l2[0], -((58.2869567709 / 100) ** 2))
    assert np.isclose(neutron_profile.neutron_flux_at(-4.79675 / 100), 159.9434), (
        "Minimum flux in FW"
    )
    assert np.isclose(neutron_profile.neutron_flux_at(4.79675 / 100), 159.9434), (
        "Minimum flux in FW"
    )
    assert np.isclose(neutron_profile.neutron_flux_at(18.96382 / 100), 164.81245), (
        "Maximum flux in BZ"
    )
    assert np.isclose(neutron_profile.neutron_flux_at(-18.96382 / 100), 164.81245), (
        "Maximum flux in BZ"
    )
    assert np.isclose(
        neutron_profile.neutron_flux_in_layer(0, x_fw),
        neutron_profile.neutron_flux_in_layer(1, x_fw),
    ), "Flux continuity assurance"
    assert np.isclose(
        neutron_profile.neutron_current_through_interface(1),
        -7.6275782637960745,
    ), "Negative current because BZ (breeding) is backflowing into the FW"
    assert np.isclose(neutron_profile.neutron_current_escaped(), 30.665951670177186), (
        "positive escaped current."
    )
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


@pytest.mark.filterwarnings("ignore:FigureCanvasAgg is non-interactive, and thus cannot be shown")
def test_fission_plot():
    """Same regression test as test_one_group_with_fission
    https://www.desmos.com/calculator/cd830add9c
    But we also plot.
    """
    dummy = [MAX_E, MIN_E]
    mfp_fw_s = 118 * 0.01  # [m]
    mfp_fw_t = 16.65 * 0.01  # [m]
    sigma_fw_t = 1 / mfp_fw_t  # [1/m]
    sigma_fw_s = 1 / mfp_fw_s  # [1/m]
    a_fw = 52
    fw_material = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="fw")
    fw_material.avg_atomic_mass = a_fw
    fw_material._set_sigma([sigma_fw_t], [[sigma_fw_s]])  # noqa: SLF001

    mfp_bz_s = 97 * 0.01  # [m]
    mfp_bz_t = 35.8 * 0.01  # [m]
    sigma_bz_s = 1 / mfp_bz_s  # [1/m]
    sigma_bz_t = 1 / mfp_bz_t  # [1/m]
    a_bz = 71

    g = 1.2
    nu_sigma_bz_f = g * (sigma_bz_t - sigma_bz_s)
    bz_material = MaterialMacroInfo(dummy, 1.0, {"Lu": 1.0}, name="bz")
    bz_material.avg_atomic_mass = a_bz
    bz_material._set_sigma([sigma_bz_t], [[sigma_bz_s]], [[nu_sigma_bz_f]])  # noqa: SLF001

    x_fw, x_bz = 5.72 * 0.01, 85 * 0.01
    incoming_flux = 41
    neutron_profile = NeutronFluxProfile(
        incoming_flux,
        [x_fw, x_bz],
        [fw_material, bz_material],
    )
    neutron_profile.plot("flux")
    neutron_profile.plot("heating")
    neutron_profile.plot("tritium_production")
    neutron_profile.plot("current")


def test_2_groups_2_layers():
    """Create a 2-layer 2-group model."""
    dummy = np.geomspace(MAX_E, MIN_E, 3)
    mat1 = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="mat1")
    mat1._set_sigma([100.0, 200], [[90, 1.0], [0.0, 80.0]])
    mat2 = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="mat2")
    mat2._set_sigma([10.0, 20], [[9, 1.0], [0.0, 8.0]])
    neutron_profile = NeutronFluxProfile(
        1.0,
        [0.2, 0.3],
        [mat1, mat2],
    )
    # Continuity
    for n in range(neutron_profile.n_groups):
        for num_layer in range(neutron_profile.n_layers-1):
            x = neutron_profile.layer_x[num_layer]
            np.testing.assert_almost_equal(
                neutron_profile.groupwise_neutron_flux_in_layer(n, num_layer, x),
                neutron_profile.groupwise_neutron_flux_in_layer(n, num_layer+1, x),
            )
            np.testing.assert_almost_equal(
                neutron_profile.groupwise_neutron_current_in_layer(n, num_layer, x),
                neutron_profile.groupwise_neutron_current_in_layer(n, num_layer+1, x),
            )


def test_2_groups_1_layer():
    dummy = np.geomspace(MAX_E, MIN_E, 3)  # dummy group structure
    # translate from mean-free-path lengths (mfp) to macroscopic cross-sections
    mfp_fw_s = 118 * 0.01  # [m]
    mfp_fw_t = 16.65 * 0.01  # [m]
    sigma_fw_s = 1 / mfp_fw_s  # [1/m]
    x_fw = 5.72 * 0.1
    fw_material = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="fw")
    fw_material._set_sigma(  # noqa: SLF001
        [1 / mfp_fw_t, 1 / (mfp_fw_t + 0.5)],
        [[sigma_fw_s, sigma_fw_s], [0.0, sigma_fw_s]],
    )
    incoming_flux = 100.0
    neutron_profile = NeutronFluxProfile(incoming_flux, [x_fw], [fw_material])
    neutron_profile.solve()
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            0, 0, neutron_profile.extended_boundary[0]
        ),
        0,
    ), "Extended boundary condition check for group 0"
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            1, 0, neutron_profile.extended_boundary[1]
        ),
        0,
    ), "Extended boundary condition check for group 1"
    num_layer = 0
    mid_point = np.mean(neutron_profile.interface_x[num_layer : num_layer + 2])
    for n in range(neutron_profile.n_groups):
        diffusion_out, total_removal, source_in = _diffusion_equation_in_layer(
            neutron_profile, n, 0, mid_point
        )
        assert np.isclose(diffusion_out, total_removal - source_in), (
            "Check that the diffusion equation holds up at an arbitrary point."
        )
    removal_xs = [
        mat.sigma_t - mat.sigma_s.sum(axis=1) for mat in neutron_profile.materials
    ]
    assert np.isclose(
        sum(neutron_profile.fluxes),
        neutron_profile.neutron_current_escaped()
        + sum(
            sum(
                removal_xs[num_layer][n]
                * neutron_profile.groupwise_integrated_flux_in_layer(n, num_layer)
                for n in range(neutron_profile.n_groups)
            )
            for num_layer in range(neutron_profile.n_layers)
        ),
    ), "Conservation of neutrons"
    assert np.isclose(neutron_profile.neutron_current_at(0), incoming_flux)

def test_2_groups_2_layers_two_identical_materials():
    """
    A 2-layer model (both layers being made of material A) should have the same
    neutron spectrum and flux profiles as a one-layer model.
    """
    dummy = np.geomspace(MAX_E, MIN_E, 3)
    mat1 = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="mat1")
    mat1._set_sigma([10.0, 20], [[9, 1.0], [0.0, 8.0]])
    mat2 = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="mat2")
    mat2._set_sigma([10.0, 20], [[9, 1.0], [0.0, 8.0]])
    neutron_profile1 = NeutronFluxProfile(
        1.0,
        [0.2, 0.4],
        [mat1, mat2],
    )
    neutron_profile2 = NeutronFluxProfile(
        1.0,
        [0.4],
        [mat1],
    )
    x = np.linspace(0, neutron_profile2.extended_boundary[0])
    np.testing.assert_array_almost_equal(
        neutron_profile1.groupwise_neutron_flux_at(0, x),
        neutron_profile2.groupwise_neutron_flux_at(0, x)
    )
    x = np.linspace(0, neutron_profile2.extended_boundary[1])
    np.testing.assert_array_almost_equal(
        neutron_profile1.groupwise_neutron_flux_at(1, x),
        neutron_profile2.groupwise_neutron_flux_at(1, x)
    )

def test_3_groups_1_layer():
    dummy = np.geomspace(MAX_E, MIN_E, 4)  # dummy group structure
    # translate from mean-free-path lengths (mfp) to macroscopic cross-sections
    mfp_fw_s = 118 * 0.01  # [m]
    mfp_fw_t = 16.65 * 0.01  # [m]
    sigma_fw_s = 1 / mfp_fw_s  # [1/m]
    x_fw = 5.72 * 0.1
    fw_material = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="fw")
    fw_material._set_sigma(  # noqa: SLF001
        [1 / mfp_fw_t, 1 / (mfp_fw_t + 0.25), 1 / (mfp_fw_t + 0.5)],
        [
            [sigma_fw_s, sigma_fw_s, sigma_fw_s],
            [0, sigma_fw_s, sigma_fw_s],
            [0, 0, sigma_fw_s],
        ],
    )
    incoming_flux = 100.0
    neutron_profile = NeutronFluxProfile(incoming_flux, [x_fw], [fw_material])
    neutron_profile.solve()
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            0, 0, neutron_profile.extended_boundary[0]
        ),
        0,
    ), "Extended boundary condition check for group 0"
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            1, 0, neutron_profile.extended_boundary[1]
        ),
        0,
    ), "Extended boundary condition check for group 1"
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            2, 0, neutron_profile.extended_boundary[2]
        ),
        0,
    ), "Extended boundary condition check for group 2"
    num_layer = 0
    mid_point = np.mean(neutron_profile.interface_x[num_layer : num_layer + 2])
    for n in range(neutron_profile.n_groups):
        diffusion_out, total_removal, source_in = _diffusion_equation_in_layer(
            neutron_profile, n, 0, mid_point
        )
        assert np.isclose(diffusion_out, total_removal - source_in), (
            "Check that the diffusion equation holds up at an arbitrary point."
        )
    removal_xs = [
        mat.sigma_t - mat.sigma_s.sum(axis=1) for mat in neutron_profile.materials
    ]
    assert np.isclose(
        sum(neutron_profile.fluxes),
        neutron_profile.neutron_current_escaped()
        + sum(
            sum(
                removal_xs[num_layer][n]
                * neutron_profile.groupwise_integrated_flux_in_layer(n, num_layer)
                for n in range(neutron_profile.n_groups)
            )
            for num_layer in range(neutron_profile.n_layers)
        ),
    ), "Conservation of neutrons"
    assert np.isclose(neutron_profile.neutron_current_at(0), incoming_flux)


def test_4_groups_1_layer():
    dummy = np.geomspace(MAX_E, MIN_E, 5)  # dummy group structure
    # translate from mean-free-path lengths (mfp) to macroscopic cross-sections
    mfp_fw_s = 118 * 0.01  # [m]
    mfp_fw_t = 16.65 * 0.01  # [m]

    sigma_fw_s = 1 / mfp_fw_s  # [1/m]
    x_fw = 5.72 * 0.1
    fw_material = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="fw")
    fw_material._set_sigma(  # noqa: SLF001
        [
            1 / mfp_fw_t,
            1 / (mfp_fw_t + 0.25),
            1 / (mfp_fw_t + 0.5),
            1 / (mfp_fw_t + 0.75),
        ],
        [
            [sigma_fw_s / 4, sigma_fw_s / 4, sigma_fw_s, sigma_fw_s],
            [0, sigma_fw_s / 3, sigma_fw_s / 3, sigma_fw_s],
            [0, 0, sigma_fw_s / 3, sigma_fw_s],
            [0, 0, 0, 0.3],
        ],
        # [
        #     [0,0.001, 0.001, 0.001],
        #     [0,0.01,0, 0],
        #     [0,0,0,0],
        #     [0,0,0,0],
        # ],
    )
    incoming_flux = 100.0
    neutron_profile = NeutronFluxProfile(incoming_flux, [x_fw], [fw_material])
    neutron_profile.solve()
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            0, 0, neutron_profile.extended_boundary[0]
        ),
        0,
    ), "Extended boundary condition check for group 0"
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            1, 0, neutron_profile.extended_boundary[1]
        ),
        0,
    ), "Extended boundary condition check for group 1"
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            2, 0, neutron_profile.extended_boundary[2]
        ),
        0,
    ), "Extended boundary condition check for group 2"
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            3, 0, neutron_profile.extended_boundary[3]
        ),
        0,
    ), "Extended boundary condition check for group 3"

    num_layer = 0
    mid_point = np.mean(neutron_profile.interface_x[num_layer : num_layer + 2])
    for n in range(neutron_profile.n_groups):
        diffusion_out, total_removal, source_in = _diffusion_equation_in_layer(
            neutron_profile, n, 0, mid_point
        )
        assert np.isclose(diffusion_out, total_removal - source_in), (
            "Check that the diffusion equation holds up at an arbitrary point."
        )
    assert np.isclose(neutron_profile.neutron_current_at(0), incoming_flux)
    removal_xs = [
        mat.sigma_t - mat.sigma_s.sum(axis=1) - mat.sigma_in.sum(axis=1)
        for mat in neutron_profile.materials
    ]
    assert np.isclose(
        sum(neutron_profile.fluxes),
        neutron_profile.neutron_current_escaped()
        + sum(
            sum(
                removal_xs[num_layer][n]
                * neutron_profile.groupwise_integrated_flux_in_layer(n, num_layer)
                for n in range(neutron_profile.n_groups)
            )
            for num_layer in range(neutron_profile.n_layers)
        ),
    ), "Conservation of neutrons"

def test_4_groups_4_layers():
    dummy = np.geomspace(MAX_E, MIN_E, 5)  # dummy group structure
    # translate from mean-free-path lengths (mfp) to macroscopic cross-sections
    tungsten = MaterialMacroInfo(dummy, 19300.0, {"Te": 1.0}, name="tungsten")
    tungsten._set_sigma(  # noqa: SLF001
        [36.70793446, 82.75605723, 557.10852765, 68.99371192],
        [
            [3.56412195e00, 3.10469266e-03, 1.86354209e-07, 9.59501248e-11],
            [0.00000000e00, 7.17103718e00, 2.01902003e-02, 3.70852436e-06],
            [0.00000000e00, 0.00000000e00, 7.95579530e00, 2.28258235e-02],
            [0.00000000e00, 0.00000000e00, 0.00000000e00, 7.98194881e00],
        ],
    )
    lithium = MaterialMacroInfo(dummy, 534.0, {"Te": 1.0}, name="lithium")
    lithium._set_sigma(  # noqa: SLF001
        [8.53058076, 5.00130054, 8.6628814, 54.01467005],
        [
            [6.2826794, 0.01282172, 0.0, 0.0],
            [0.0, 4.41033136, 0.02001536, 0.0],
            [0.0, 0.0, 4.42975314, 0.02010321],
            [0.0, 0.0, 0.0, 4.22604695],
        ],
    )
    ss316 = MaterialMacroInfo(dummy, 7930.0, {"Te": 1.0}, name="ss316")
    ss316._set_sigma(  # noqa: SLF001
        [58.75691059, 274.00378023, 3562.92253327, 86.75244242],
        [
            [3.48853683e01, 1.01273247e-02, 0.00000000e00, 0.00000000e00],
            [0.00000000e00, 2.70925014e02, 1.66745243e-01, 0.00000000e00],
            [0.00000000e00, 0.00000000e00, 3.46806397e03, 2.03609085e00],
            [0.00000000e00, 0.00000000e00, 0.00000000e00, 4.68049763e01],
        ],
    )
    concrete = MaterialMacroInfo(dummy, 3600.0, {"Te": 1.0}, name="concrete")
    concrete._set_sigma(  # noqa: SLF001
        [4.12742084, 7.42553891, 8.24183468, 8.29305646],
        [
            [3.56412195e00, 3.10469266e-03, 1.86354209e-07, 9.59501248e-11],
            [0.00000000e00, 7.17103718e00, 2.01902003e-02, 3.70852436e-06],
            [0.00000000e00, 0.00000000e00, 7.95579530e00, 2.28258235e-02],
            [0.00000000e00, 0.00000000e00, 0.00000000e00, 7.98194881e00],
        ],
    )
    incoming_flux = 100.0
    neutron_profile = NeutronFluxProfile(
        incoming_flux,
        np.cumsum([0.05, 0.3, 0.2, 0.4]),
        [tungsten, lithium, ss316, concrete],
    )
    neutron_profile.solve()
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            0, 0, neutron_profile.extended_boundary[0]
        ),
        0,
    ), "Extended boundary condition check for group 0"
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            1, 0, neutron_profile.extended_boundary[1]
        ),
        0,
    ), "Extended boundary condition check for group 1"
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            2, 0, neutron_profile.extended_boundary[2]
        ),
        0,
    ), "Extended boundary condition check for group 2"
    assert np.isclose(
        neutron_profile.groupwise_neutron_flux_in_layer(
            3, 0, neutron_profile.extended_boundary[3]
        ),
        0,
    ), "Extended boundary condition check for group 3"

    num_layer = 0
    mid_point = np.mean(neutron_profile.interface_x[num_layer : num_layer + 2])
    for n in range(neutron_profile.n_groups):
        diffusion_out, total_removal, source_in = _diffusion_equation_in_layer(
            neutron_profile, n, 0, mid_point
        )
        assert np.isclose(diffusion_out, total_removal - source_in), (
            "Check that the diffusion equation holds up at an arbitrary point."
        )
    assert np.isclose(neutron_profile.neutron_current_at(0), incoming_flux)
    removal_xs = [
        mat.sigma_t - mat.sigma_s.sum(axis=1) - mat.sigma_in.sum(axis=1)
        for mat in neutron_profile.materials
    ]
    assert np.isclose(
        sum(neutron_profile.fluxes),
        neutron_profile.neutron_current_escaped()
        + sum(
            sum(
                removal_xs[num_layer][n]
                * neutron_profile.groupwise_integrated_flux_in_layer(n, num_layer)
                for n in range(neutron_profile.n_groups)
            )
            for num_layer in range(neutron_profile.n_layers)
        ),
    ), "Conservation of neutrons"


@pytest.mark.filterwarnings("ignore:Calculation of flux")
def test_5_groups_5_layers():
    """Create an arbitrary 5-layer 5-group model. Check for continuity and conformity to the equation."""
    dummy_group_structure = np.geomspace(MAX_E, MIN_E, 5 + 1)
    mat_list = []
    at_masses = np.geomspace(1, 100, 5)[[3, 0, 2, 4, 1]]
    sigma_t_lists = [  # arbitrarily chosen and rearranged numbers
        1 / (80 + np.linspace(-40, 40, 5)[[0, 2, 4, 1, 3]]),
        1 / (200 + np.linspace(-30, 30, 5)[[0, 3, 1, 4, 2]]),
        1 / (300 + np.linspace(-20, 20, 5)[[1, 3, 0, 2, 4]]),
        1 / (100 + np.linspace(-40, 40, 5)[[1, 4, 2, 0, 3]]),
        1 / (50 + np.linspace(-10, 10, 5)[[2, 0, 3, 1, 4]]),
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
        np.where([0, 0, 1, 1, 0], sigma_t_lists[1] * 1.0, [0, 0, 0, 0, 0]),
        np.where([0, 1, 0, 1, 0], sigma_t_lists[2] * 0.6, [0, 0, 0, 0, 0]),
        np.where([1, 0, 1, 1, 1], sigma_t_lists[3] * 2.2, [0, 0, 0, 0, 0]),
        [0, 0, 0, 0, 0],
    ]
    incident_energy = np.mean(sorted(dummy_group_structure)[-2:])
    dummy_group_energy, _ = calculate_mean_energy_and_incident_bin(
        dummy_group_structure, incident_energy
    )
    for i in range(5):
        mat = MaterialMacroInfo(dummy_group_structure, 1.0, {"C": 1.0}, name=f"mat{i}")
        mat.avg_atomic_mass = at_masses[i]
        mat._set_sigma(  # noqa: SLF001
            sigma_t=sigma_t_lists[i],
            sigma_s=(
                sigma_s_list[i]
                * scattering_weight_matrix(
                    dummy_group_structure, dummy_group_energy, at_masses[i]
                ).T
            ).T,
            sigma_in=(
                sigma_in_list[i]
                * scattering_weight_matrix(
                    dummy_group_structure, dummy_group_energy, at_masses[i]
                ).T
            ).T,
        )
        mat_list.append(mat)
    incoming_flux = 100.0
    neutron_profile = NeutronFluxProfile(incoming_flux, [5, 10, 15, 20, 25], mat_list)
    neutron_profile.solve()
    for num_layer in range(neutron_profile.n_layers):
        mid_point = np.mean(neutron_profile.interface_x[num_layer : num_layer + 2])
        layer_x = neutron_profile.layer_x[num_layer]
        for n in range(neutron_profile.n_groups):
            # Check for conformity with the diffusion equation
            diffusion_out, total_removal, source_in = _diffusion_equation_in_layer(
                neutron_profile, n, num_layer, mid_point
            )
            assert np.isclose(diffusion_out, total_removal - source_in), (
                "Check that the diffusion equation holds up at an arbitrary point."
            )
            if num_layer == neutron_profile.n_layers - 1:
                continue
            # Check for continuity of flux and current
            assert np.isclose(
                neutron_profile.groupwise_neutron_flux_in_layer(n, num_layer, layer_x),
                neutron_profile.groupwise_neutron_flux_in_layer(
                    n, num_layer + 1, layer_x
                ),
            )
            assert np.isclose(
                neutron_profile.groupwise_neutron_current_in_layer(
                    n, num_layer, layer_x
                ),
                neutron_profile.groupwise_neutron_current_in_layer(
                    n, num_layer + 1, layer_x
                ),
            )
        assert np.isclose(
            neutron_profile.neutron_flux_in_layer(num_layer, layer_x),
            neutron_profile.neutron_flux_in_layer(num_layer + 1, layer_x),
        )
        assert np.isclose(
            neutron_profile.neutron_current_in_layer(num_layer, layer_x),
            neutron_profile.neutron_current_in_layer(num_layer + 1, layer_x),
        )
    # Check for extended boundary flux = 0
    num_layer = neutron_profile.n_layers - 1
    for n in range(neutron_profile.n_groups):
        assert np.isclose(
            neutron_profile.groupwise_neutron_flux_in_layer(
                n, num_layer, neutron_profile.extended_boundary[n]
            ),
            0,
        ), f"flux at Extended boundary of group {n} should = 0"

    no_incident_flux_err_msg = (
        "Expected no incident neutron flux from the plasma except in energy group 0."
    )
    for n in range(neutron_profile.n_groups):
        assert np.isclose(
            neutron_profile.groupwise_neutron_current_at(n, 0),
            incoming_flux * int(n == 0),
        ), no_incident_flux_err_msg

    sigma_t = np.array([mat.sigma_t for mat in neutron_profile.materials])
    sigma_s = np.array([mat.sigma_s for mat in neutron_profile.materials])
    sigma_in = np.array([mat.sigma_in for mat in neutron_profile.materials])
    shape = np.array([neutron_profile.n_layers, neutron_profile.n_groups])
    in_flow = np.zeros(shape)
    in_scatter = np.zeros(shape)
    removal = np.zeros(shape)
    for n in range(neutron_profile.n_groups):
        for num_layer in range(neutron_profile.n_layers):
            in_flow[num_layer, n] = neutron_profile.groupwise_neutron_current_in_layer(
                n, num_layer, neutron_profile.interface_x[num_layer]
            ) - neutron_profile.groupwise_neutron_current_in_layer(
                n, num_layer, neutron_profile.interface_x[num_layer + 1]
            )
            in_scatter[num_layer, n] = sum(
                (sigma_s[num_layer, in_group, n] + sigma_in[num_layer, in_group, n])
                * neutron_profile.groupwise_integrated_flux_in_layer(in_group, num_layer)
                for in_group in range(neutron_profile.n_groups)
                if in_group < n
            )
            removal[num_layer, n] = (
                sigma_t[num_layer, n]
                - sigma_s[num_layer, n, n]
                - sigma_in[num_layer, n, n]
            ) * neutron_profile.groupwise_integrated_flux_in_layer(n, num_layer)
        assert np.isclose(
            neutron_profile.groupwise_neutron_current_through_interface(
                n, num_layer + 1
            ),
            neutron_profile.groupwise_neutron_current_escaped(n),
        )
        assert np.isclose(
            neutron_profile.groupwise_neutron_current_at(
                n, neutron_profile.layer_x[num_layer]
            ),
            neutron_profile.groupwise_neutron_current_through_interface(
                n, num_layer + 1
            ),
        )
    assert np.isclose(in_flow + in_scatter, removal, atol=0, rtol=1e-9).all(), (
        f"Mismatch between {num_layer} group {n} influx and outflux"
    )

    removal_xs = np.zeros(shape)
    int_flux = np.zeros(shape)
    for num_layer in range(neutron_profile.n_layers):
        removal_xs[num_layer] = (
            neutron_profile.materials[num_layer].sigma_t
            # Count all neutrons consumed, but not if they popped up anywhere else.
            - neutron_profile.materials[num_layer].sigma_s.sum(axis=1)
            - neutron_profile.materials[num_layer].sigma_in.sum(axis=1)
        )
        int_flux[num_layer] = [
            neutron_profile.groupwise_integrated_flux_in_layer(n, num_layer)
            for n in range(neutron_profile.n_groups)
        ]

    assert np.isclose(
        sum(neutron_profile.fluxes),
        neutron_profile.neutron_current_escaped() + (removal_xs * int_flux).sum(),
    ), "Conservation of neutrons"

