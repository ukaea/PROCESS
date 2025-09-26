import pytest

from process.current_drive import (
    CurrentDrive,
    ElectronBernstein,
    ElectronCyclotron,
    IonCyclotron,
    LowerHybrid,
    NeutralBeam,
)
from process.data_structure import (
    current_drive_variables,
    heat_transport_variables,
    physics_variables,
)
from process.plasma_profiles import PlasmaProfile


@pytest.fixture
def current_drive():
    """Provides CurrentDrive object for testing.

    :returns current_drive: initialised CurrentDrive object
    :rtype: process.current_drive.CurrentDrive
    """
    return CurrentDrive(
        PlasmaProfile(),
        electron_cyclotron=ElectronCyclotron(plasma_profile=PlasmaProfile()),
        electron_bernstein=ElectronBernstein(plasma_profile=PlasmaProfile()),
        neutral_beam=NeutralBeam(plasma_profile=PlasmaProfile()),
        lower_hybrid=LowerHybrid(plasma_profile=PlasmaProfile()),
        ion_cyclotron=IonCyclotron(plasma_profile=PlasmaProfile()),
    )


def test_cudriv_primary_lower_hybrid(current_drive):
    current_drive_variables.i_hcd_primary = 1  # Lower Hybrid
    current_drive_variables.i_hcd_secondary = 0
    current_drive_variables.i_hcd_calculations = 1
    current_drive_variables.p_hcd_primary_extra_heat_mw = 0.0
    current_drive_variables.eta_cd_hcd_secondary = 0.0
    physics_variables.nd_plasma_electrons_vol_avg = 1e20
    physics_variables.temp_plasma_electron_vol_avg_keV = 10
    physics_variables.rmajor = 6.2
    physics_variables.plasma_current = 15e6
    physics_variables.f_c_plasma_auxiliary = 0.2
    current_drive_variables.eta_lowhyb_injector_wall_plug = 0.4
    current_drive.cudriv()

    assert current_drive_variables.eta_cd_hcd_primary == pytest.approx(
        0.07812311, rel=1e-6
    )
    assert current_drive_variables.eta_cd_norm_hcd_primary == pytest.approx(
        0.48436326, rel=1e-6
    )
    assert current_drive_variables.f_c_plasma_hcd_secondary == pytest.approx(
        0.0, rel=1e-6
    )
    assert current_drive_variables.p_hcd_primary_injected_mw == pytest.approx(
        38.4009309, rel=1e-6
    )
    assert current_drive_variables.c_hcd_primary_driven == pytest.approx(
        3000000.0, rel=1e-6
    )
    assert current_drive_variables.f_c_plasma_hcd_primary == pytest.approx(
        0.2, rel=1e-6
    )
    assert heat_transport_variables.p_hcd_primary_electric_mw == pytest.approx(
        96.00232726, rel=1e-6
    )
    assert heat_transport_variables.p_hcd_secondary_electric_mw == pytest.approx(
        0.0, rel=1e-6
    )
    assert current_drive_variables.p_hcd_injected_total_mw == pytest.approx(
        38.4009309, rel=1e-6
    )
    assert current_drive_variables.p_hcd_injected_electrons_mw == pytest.approx(
        38.4009309, rel=1e-6
    )


def test_cudriv_primary_lower_hybrid_with_heat(current_drive):
    current_drive_variables.i_hcd_primary = 1  # Lower Hybrid
    current_drive_variables.i_hcd_secondary = 0
    current_drive_variables.i_hcd_calculations = 1
    current_drive_variables.p_hcd_primary_extra_heat_mw = 5.0  # Adding primary heat
    current_drive_variables.eta_cd_hcd_secondary = 0.0
    physics_variables.nd_plasma_electrons_vol_avg = 1e20
    physics_variables.temp_plasma_electron_vol_avg_keV = 10
    physics_variables.rmajor = 6.2
    physics_variables.plasma_current = 15e6
    physics_variables.f_c_plasma_auxiliary = 0.2
    current_drive_variables.eta_lowhyb_injector_wall_plug = 0.4
    current_drive.cudriv()

    assert current_drive_variables.eta_cd_hcd_primary == pytest.approx(
        0.07812311, rel=1e-6
    )
    assert current_drive_variables.eta_cd_norm_hcd_primary == pytest.approx(
        0.48436326, rel=1e-6
    )
    assert current_drive_variables.f_c_plasma_hcd_secondary == pytest.approx(
        0.0, rel=1e-6
    )
    assert current_drive_variables.p_hcd_primary_injected_mw == pytest.approx(
        38.4009309, rel=1e-6
    )  # Adjusted for extra heat
    assert current_drive_variables.c_hcd_primary_driven == pytest.approx(
        3000000.0, rel=1e-6
    )
    assert current_drive_variables.f_c_plasma_hcd_primary == pytest.approx(
        0.2, rel=1e-6
    )
    assert heat_transport_variables.p_hcd_primary_electric_mw == pytest.approx(
        108.50232725752629, rel=1e-6
    )  # Adjusted for extra heat
    assert heat_transport_variables.p_hcd_secondary_electric_mw == pytest.approx(
        0.0, rel=1e-6
    )
    assert current_drive_variables.p_hcd_injected_total_mw == pytest.approx(
        43.4009309, rel=1e-6
    )  # Adjusted for extra heat
    assert current_drive_variables.p_hcd_injected_electrons_mw == pytest.approx(
        43.4009309, rel=1e-6
    )
    assert current_drive_variables.p_hcd_injected_ions_mw == pytest.approx(
        0.0, rel=1e-6
    )


def test_sigbeam(current_drive):
    assert current_drive.neutral_beam.sigbeam(
        1e3, 13.07, 8.0e-1, 0.1, 1e-4, 1e-4, 1e-4
    ) == pytest.approx(2.013589662302492e-11)
