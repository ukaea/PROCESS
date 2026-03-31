import numpy as np
import pytest

from process.core.exceptions import ProcessValidationError
from process.models.neutronics.base import (
    AutoPopulatingDict,
    LayerSpecificGroupwiseConstants,
    NeutronFluxProfile,
    _get_sign_of,
)
from process.models.neutronics.data import DT_NEUTRON_E, EV_TO_J, MaterialMacroInfo

MAX_E = DT_NEUTRON_E * 1.01
MIN_E = 1 / 40 * EV_TO_J


def test_group_structure_0_energy():
    with pytest.warns():
        mat = MaterialMacroInfo([1.0, 0.0], 1.0, {"C": 1.0})
        mat._set_sigma([1.0], [[1.0]])  # noqa: SLF001


def test_group_structure_too_short():
    with pytest.raises(ProcessValidationError):
        mat = MaterialMacroInfo([1.0], 1.0, {"C": 1.0})
        mat._set_sigma([1.0], [[1.0]])  # noqa: SLF001


def test_sigma_s_incorrect_shape():
    with pytest.raises(ProcessValidationError):
        mat = MaterialMacroInfo([1000, 10, 1.0], 1.0, {"C": 1.0})
        mat._set_sigma([1.0, 2.0], [1.0, 1.0])  # noqa: SLF001


def test_sigma_s_too_large():
    with pytest.raises(ProcessValidationError):
        mat = MaterialMacroInfo([1000, 10, 1.0], 1.0, {"C": 1.0})
        mat._set_sigma([1.0, 2.0], [[1.0, 1.0], [1.0, 1.0]])  # noqa: SLF001


def test_warn_up_elastic_scatter():
    with pytest.warns():
        mat = MaterialMacroInfo([1000, 10, 1.0], 1.0, {"C": 1.0})
        mat._set_sigma([1.0, 2.0], [[0.5, 0.5], [1.0, 1.0]])  # noqa: SLF001


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
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_current_through_interface")
    assert hasattr(NeutronFluxProfile, "neutron_current_through_interface")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_current_escaped")
    assert hasattr(NeutronFluxProfile, "neutron_current_escaped")


def test_has_reaction_rates():
    """Test that the groupwise decorator has worked on the reactions methods."""
    assert hasattr(NeutronFluxProfile, "groupwise_integrated_flux_in_layer")
    assert hasattr(NeutronFluxProfile, "integrated_flux_in_layer")
    assert hasattr(NeutronFluxProfile, "groupwise_integrated_heating_in_layer")
    assert hasattr(NeutronFluxProfile, "integrated_heating_in_layer")
    assert hasattr(
        NeutronFluxProfile, "groupwise_integrated_tritium_production_in_layer"
    )
    assert hasattr(NeutronFluxProfile, "integrated_tritium_production_in_layer")


def test_has_volumetric_heating():
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_heating_in_layer")
    assert hasattr(NeutronFluxProfile, "neutron_heating_in_layer")
    assert hasattr(NeutronFluxProfile, "groupwise_neutron_heating_at")
    assert hasattr(NeutronFluxProfile, "neutron_heating_at")


def test_has_tritium_production():
    assert hasattr(NeutronFluxProfile, "groupwise_tritium_production_in_layer")
    assert hasattr(NeutronFluxProfile, "tritium_production_in_layer")
    assert hasattr(NeutronFluxProfile, "groupwise_tritium_production_at")
    assert hasattr(NeutronFluxProfile, "tritium_production_at")


def test_has_plot():
    assert hasattr(NeutronFluxProfile, "plot")


def test_units():
    nfp = NeutronFluxProfile
    assert nfp.get_output_unit(nfp.groupwise_integrated_flux_in_layer) == "m^-1 s^-1"
    assert nfp.get_output_unit(nfp.groupwise_integrated_heating_in_layer) == "W m^-2"
    assert (
        nfp.get_output_unit(nfp.groupwise_integrated_tritium_production_in_layer)
        == "mole m^-2"
    )
    assert nfp.get_output_unit(nfp.groupwise_linear_heating_density_in_layer) == "J m^-1"
    assert nfp.get_output_unit(nfp.groupwise_neutron_current_at) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.groupwise_neutron_current_escaped) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.groupwise_neutron_current_in_layer) == "m^-2 s^-1"
    assert (
        nfp.get_output_unit(nfp.groupwise_neutron_current_through_interface)
        == "m^-2 s^-1"
    )
    assert nfp.get_output_unit(nfp.groupwise_neutron_flux_at) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.groupwise_neutron_flux_in_layer) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.groupwise_neutron_heating_at) == "W m^-3"
    assert nfp.get_output_unit(nfp.groupwise_neutron_heating_in_layer) == "W m^-3"
    assert nfp.get_output_unit(nfp.groupwise_tritium_production_at) == "mole m^-3"
    assert nfp.get_output_unit(nfp.groupwise_tritium_production_in_layer) == "mole m^-3"

    assert nfp.get_output_unit(nfp.integrated_flux_in_layer) == "m^-1 s^-1"
    assert nfp.get_output_unit(nfp.integrated_heating_in_layer) == "W m^-2"
    assert nfp.get_output_unit(nfp.integrated_tritium_production_in_layer) == "mole m^-2"
    assert nfp.get_output_unit(nfp.neutron_current_at) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_current_escaped) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_current_in_layer) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_current_through_interface) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_flux_at) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_flux_in_layer) == "m^-2 s^-1"
    assert nfp.get_output_unit(nfp.neutron_heating_in_layer) == "W m^-3"
    assert nfp.get_output_unit(nfp.neutron_heating_at) == "W m^-3"
    assert nfp.get_output_unit(nfp.tritium_production_in_layer) == "mole m^-3"


def test_get_sign_func():
    signs = _get_sign_of(np.array([-np.inf, -2, -1, -0.0, 0.0, 1.0, 2.0, np.inf]))
    np.testing.assert_equal(signs, [-1, -1, -1, -1, 1, 1, 1, 1])


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


def test_same_l_in_2_groups_warns():
    dummy = np.geomspace(MAX_E, MIN_E, 3)  # dummy group structure
    # translate from mean-free-path lengths (mfp) to macroscopic cross-sections
    mfp_fw_s = 118 * 0.01  # [m]
    mfp_fw_t = 16.65 * 0.01  # [m]
    sigma_fw_t = 1 / mfp_fw_t  # [1/m]
    sigma_fw_s = 1 / mfp_fw_s  # [1/m]
    x_fw = 5.72 * 0.01
    fw_material = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="fw")
    fw_material._set_sigma(  # noqa: SLF001
        [sigma_fw_t, sigma_fw_t], [[sigma_fw_s, sigma_fw_s], [0.0, sigma_fw_s]]
    )
    incoming_flux = 100.0
    neutron_profile = NeutronFluxProfile(incoming_flux, [x_fw], [fw_material])
    with pytest.warns(UserWarning):
        neutron_profile.solve()


@pytest.mark.filterwarnings(
    "ignore:Group 0 and group 1 has the same neutron diffusion lengths"
)
@pytest.mark.filterwarnings("error")
@pytest.mark.xfail
def test_same_l_in_2_groups_calculate_flux():
    dummy = np.geomspace(MAX_E, MIN_E, 3)  # dummy group structure
    # translate from mean-free-path lengths (mfp) to macroscopic cross-sections
    mfp_fw_s = 118 * 0.01  # [m]
    mfp_fw_t = 16.65 * 0.01  # [m]
    sigma_fw_t = 1 / mfp_fw_t  # [1/m]
    sigma_fw_s = 1 / mfp_fw_s  # [1/m]
    x_fw = 5.72 * 0.01
    fw_material = MaterialMacroInfo(dummy, 1.0, {"Te": 1.0}, name="fw")
    fw_material._set_sigma(  # noqa: SLF001
        [sigma_fw_t, sigma_fw_t], [[sigma_fw_s, sigma_fw_s], [0.0, sigma_fw_s]]
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


def test_two_group():
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


@pytest.mark.filterwarnings("ignore:Calculation of flux")
def test_three_group():
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


def test_four_group():
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
