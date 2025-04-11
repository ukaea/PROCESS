from typing import Any, NamedTuple

import pytest

from process.current_drive import (
    CurrentDrive,
    ElectronBernstein,
    ElectronCyclotron,
    IonCyclotron,
    LowerHybrid,
    NeutralBeam,
)
from process.fortran import (
    cost_variables,
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


class CudrivParam(NamedTuple):
    p_hcd_secondary_electric_mw: Any = None

    p_hcd_electric_total_mw: Any = None

    p_hcd_ecrh_injected_total_mw: Any = None

    p_hcd_beam_injected_total_mw: Any = None

    p_hcd_lowhyb_injected_total_mw: Any = None

    c_beam_total: Any = None

    p_beam_orbit_loss_mw: Any = None

    i_hcd_primary: Any = None

    i_hcd_secondary: Any = None

    p_hcd_primary_extra_heat_mw: Any = None

    p_hcd_secondary_extra_heat_mw: Any = None

    p_hcd_secondary_injected_mw: Any = None

    i_hcd_calculations: Any = None

    feffcd: Any = None

    f_p_beam_injected_ions: Any = None

    f_p_beam_shine_through: Any = None

    eta_cd_norm_hcd_primary: Any = None

    eta_cd_norm_ecrh: Any = None

    eta_lowhyb_injector_wall_plug: Any = None

    eta_hcd_primary_injector_wall_plug: Any = None

    eta_hcd_secondary_injector_wall_plug: Any = None

    eta_ecrh_injector_wall_plug: Any = None

    f_p_beam_orbit_loss: Any = None

    p_hcd_injected_total_mw: Any = None

    pwpnb: Any = None

    eta_beam_injector_wall_plug: Any = None

    e_beam_kev: Any = None

    eta_cd_hcd_primary: Any = None

    p_hcd_lowhyb_electric_mw: Any = None

    p_hcd_ecrh_electric_mw: Any = None

    p_beam_injected_mw: Any = None

    p_beam_shine_through_mw: Any = None

    p_hcd_injected_electrons_mw: Any = None

    p_hcd_injected_ions_mw: Any = None

    bigq: Any = None

    f_c_plasma_bootstrap: Any = None

    f_c_plasma_bootstrap_max: Any = None

    n_beam_decay_lengths_core: Any = None

    p_hcd_injected_max: Any = None

    dx_beam_shield: Any = None

    frbeam: Any = None

    rtanbeam: Any = None

    rtanmax: Any = None

    f_c_plasma_diamagnetic: Any = None

    f_c_plasma_pfirsch_schluter: Any = None

    f_c_plasma_internal: Any = None

    n_ecrh_harmonic: Any = None

    xi_ebw: Any = None

    dene: Any = None

    te: Any = None

    rmajor: Any = None

    ten: Any = None

    zeff: Any = None

    dlamee: Any = None

    beta: Any = None

    rhopedt: Any = None

    rhopedn: Any = None

    te0: Any = None

    teped: Any = None

    tesep: Any = None

    alphat: Any = None

    alphan: Any = None

    ne0: Any = None

    nesep: Any = None

    neped: Any = None

    bt: Any = None

    rminor: Any = None

    tbeta: Any = None

    plasma_current: Any = None

    ipedestal: Any = None

    f_c_plasma_auxiliary: Any = None

    ignite: Any = None

    p_plasma_ohmic_mw: Any = None

    fusion_power: Any = None

    f_c_plasma_inductive: Any = None

    f_c_plasma_non_inductive: Any = None

    startupratio: Any = None

    iprint: Any = None

    outfile: Any = None

    expected_p_hcd_electric_total_mw: Any = None

    expected_p_hcd_ecrh_injected_total_mw: Any = None

    expected_gamcd: Any = None

    expected_etacd: Any = None

    expected_p_hcd_injected_total_mw: Any = None

    expected_effcd: Any = None

    expected_p_hcd_ecrh_electric_mw: Any = None

    expected_p_hcd_injected_electrons_mw: Any = None

    expected_bigq: Any = None


@pytest.mark.parametrize(
    "cudrivparam",
    (
        CudrivParam(
            p_hcd_secondary_electric_mw=0,
            p_hcd_electric_total_mw=0,
            p_hcd_ecrh_injected_total_mw=0,
            p_hcd_beam_injected_total_mw=0,
            p_hcd_lowhyb_injected_total_mw=0,
            c_beam_total=0,
            p_beam_orbit_loss_mw=0,
            i_hcd_primary=10,
            i_hcd_secondary=0,
            p_hcd_primary_extra_heat_mw=75,
            p_hcd_secondary_extra_heat_mw=0,
            p_hcd_secondary_injected_mw=0,
            i_hcd_calculations=1,
            feffcd=1,
            f_p_beam_injected_ions=0.5,
            f_p_beam_shine_through=0,
            eta_cd_norm_hcd_primary=0,
            eta_cd_norm_ecrh=0.30000000000000004,
            eta_lowhyb_injector_wall_plug=0.29999999999999999,
            eta_hcd_primary_injector_wall_plug=0,
            eta_hcd_secondary_injector_wall_plug=0,
            eta_ecrh_injector_wall_plug=0.5,
            f_p_beam_orbit_loss=0,
            p_hcd_injected_total_mw=0,
            pwpnb=0,
            eta_beam_injector_wall_plug=0.29999999999999999,
            e_beam_kev=1000,
            eta_cd_hcd_primary=0,
            p_hcd_lowhyb_electric_mw=0,
            p_hcd_ecrh_electric_mw=0,
            p_beam_injected_mw=0,
            p_beam_shine_through_mw=0,
            p_hcd_injected_electrons_mw=0,
            p_hcd_injected_ions_mw=0,
            bigq=0,
            f_c_plasma_bootstrap=0.27635918746616817,
            f_c_plasma_bootstrap_max=0.95000000000000007,
            n_beam_decay_lengths_core=0,
            p_hcd_injected_max=200,
            dx_beam_shield=0.5,
            frbeam=1.05,
            rtanbeam=0,
            rtanmax=0,
            f_c_plasma_diamagnetic=0,
            f_c_plasma_pfirsch_schluter=0,
            f_c_plasma_internal=0.27635918746616817,
            n_ecrh_harmonic=1,
            xi_ebw=0.80000000000000004,
            dene=7.5e19,
            te=12,
            rmajor=8,
            ten=12.626131115905864,
            zeff=2.0909945616489103,
            dlamee=17.510652035055571,
            beta=0.030000000000000006,
            rhopedt=0.94000000000000006,
            rhopedn=0.94000000000000006,
            te0=24.402321098330372,
            teped=5.5,
            tesep=0.10000000000000001,
            alphat=1.45,
            alphan=1,
            ne0=8.515060981068918e19,
            nesep=4.1177885154594193e19,
            neped=7.000240476281013e19,
            bt=5.7000000000000002,
            rminor=2.6666666666666665,
            tbeta=2,
            plasma_current=18398455.678867526,
            ipedestal=1,
            f_c_plasma_auxiliary=0.12364081253383186,
            ignite=0,
            p_plasma_ohmic_mw=0,
            fusion_power=0,
            f_c_plasma_inductive=0.59999999999999998,
            f_c_plasma_non_inductive=0.40000000000000002,
            startupratio=1,
            iprint=0,
            outfile=11,
            expected_p_hcd_electric_total_mw=240.99200038011492,
            expected_p_hcd_ecrh_injected_total_mw=120.49600019005746,
            expected_gamcd=0.30000000000000004,
            expected_etacd=0.5,
            expected_p_hcd_injected_total_mw=120.49600019005746,
            expected_effcd=0.05000000000000001,
            expected_p_hcd_ecrh_electric_mw=240.99200038011492,
            expected_p_hcd_injected_electrons_mw=120.49600019005746,
            expected_bigq=0,
        ),
        CudrivParam(
            p_hcd_secondary_electric_mw=0,
            p_hcd_electric_total_mw=240.99200038011492,
            p_hcd_ecrh_injected_total_mw=120.49600019005746,
            p_hcd_beam_injected_total_mw=0,
            p_hcd_lowhyb_injected_total_mw=0,
            c_beam_total=0,
            p_beam_orbit_loss_mw=0,
            i_hcd_primary=10,
            i_hcd_secondary=0,
            p_hcd_primary_extra_heat_mw=75,
            p_hcd_secondary_extra_heat_mw=0,
            p_hcd_secondary_injected_mw=0,
            i_hcd_calculations=1,
            feffcd=1,
            f_p_beam_injected_ions=0.5,
            f_p_beam_shine_through=0,
            eta_cd_norm_hcd_primary=0.30000000000000004,
            eta_cd_norm_ecrh=0.30000000000000004,
            eta_lowhyb_injector_wall_plug=0.29999999999999999,
            eta_hcd_primary_injector_wall_plug=0.5,
            eta_hcd_secondary_injector_wall_plug=0,
            eta_ecrh_injector_wall_plug=0.5,
            f_p_beam_orbit_loss=0,
            p_hcd_injected_total_mw=120.49600019005746,
            pwpnb=0,
            eta_beam_injector_wall_plug=0.29999999999999999,
            e_beam_kev=1000,
            eta_cd_hcd_primary=0.05000000000000001,
            p_hcd_lowhyb_electric_mw=0,
            p_hcd_ecrh_electric_mw=240.99200038011492,
            p_beam_injected_mw=0,
            p_beam_shine_through_mw=0,
            p_hcd_injected_electrons_mw=120.49600019005746,
            p_hcd_injected_ions_mw=0,
            bigq=0,
            f_c_plasma_bootstrap=0.27635918746616817,
            f_c_plasma_bootstrap_max=0.95000000000000007,
            n_beam_decay_lengths_core=0,
            p_hcd_injected_max=200,
            dx_beam_shield=0.5,
            frbeam=1.05,
            rtanbeam=8.4000000000000004,
            rtanmax=13.179564451855533,
            f_c_plasma_diamagnetic=0,
            f_c_plasma_pfirsch_schluter=0,
            f_c_plasma_internal=0.27635918746616817,
            n_ecrh_harmonic=1,
            xi_ebw=0.80000000000000004,
            dene=7.5e19,
            te=12,
            rmajor=8,
            ten=12.626131115905864,
            zeff=2.0909945616489103,
            dlamee=17.510652035055571,
            beta=0.030000000000000006,
            rhopedt=0.94000000000000006,
            rhopedn=0.94000000000000006,
            te0=24.402321098330372,
            teped=5.5,
            tesep=0.10000000000000001,
            alphat=1.45,
            alphan=1,
            ne0=8.515060981068918e19,
            nesep=4.1177885154594193e19,
            neped=7.000240476281013e19,
            bt=5.7000000000000002,
            rminor=2.6666666666666665,
            tbeta=2,
            plasma_current=18398455.678867526,
            ipedestal=1,
            f_c_plasma_auxiliary=0.12364081253383186,
            ignite=0,
            p_plasma_ohmic_mw=0.76707314489379119,
            fusion_power=1051.6562748933977,
            f_c_plasma_inductive=0.59999999999999998,
            f_c_plasma_non_inductive=0.40000000000000002,
            startupratio=1,
            iprint=0,
            outfile=11,
            expected_p_hcd_electric_total_mw=240.99200038011492,
            expected_p_hcd_ecrh_injected_total_mw=120.49600019005746,
            expected_gamcd=0.30000000000000004,
            expected_etacd=0.5,
            expected_p_hcd_injected_total_mw=120.49600019005746,
            expected_effcd=0.05000000000000001,
            expected_p_hcd_ecrh_electric_mw=240.99200038011492,
            expected_p_hcd_injected_electrons_mw=120.49600019005746,
            expected_bigq=8.6725187311435423,
        ),
    ),
)
def test_cudriv(cudrivparam, monkeypatch, current_drive):
    """
    Automatically generated Regression Unit Test for cudriv.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param cudrivparam: the data used to mock and assert in this test.
    :type cudrivparam: cudrivparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_secondary_electric_mw",
        cudrivparam.p_hcd_secondary_electric_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_total_mw",
        cudrivparam.p_hcd_electric_total_mw,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_ecrh_injected_total_mw",
        cudrivparam.p_hcd_ecrh_injected_total_mw,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_beam_injected_total_mw",
        cudrivparam.p_hcd_beam_injected_total_mw,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_lowhyb_injected_total_mw",
        cudrivparam.p_hcd_lowhyb_injected_total_mw,
    )

    monkeypatch.setattr(
        current_drive_variables, "c_beam_total", cudrivparam.c_beam_total
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_beam_orbit_loss_mw",
        cudrivparam.p_beam_orbit_loss_mw,
    )

    monkeypatch.setattr(
        current_drive_variables, "i_hcd_primary", cudrivparam.i_hcd_primary
    )

    monkeypatch.setattr(
        current_drive_variables, "i_hcd_secondary", cudrivparam.i_hcd_secondary
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_primary_extra_heat_mw",
        cudrivparam.p_hcd_primary_extra_heat_mw,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_secondary_extra_heat_mw",
        cudrivparam.p_hcd_secondary_extra_heat_mw,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_secondary_injected_mw",
        cudrivparam.p_hcd_secondary_injected_mw,
    )

    monkeypatch.setattr(
        current_drive_variables, "i_hcd_calculations", cudrivparam.i_hcd_calculations
    )

    monkeypatch.setattr(current_drive_variables, "feffcd", cudrivparam.feffcd)

    monkeypatch.setattr(
        current_drive_variables,
        "f_p_beam_injected_ions",
        cudrivparam.f_p_beam_injected_ions,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "f_p_beam_shine_through",
        cudrivparam.f_p_beam_shine_through,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "eta_cd_norm_hcd_primary",
        cudrivparam.eta_cd_norm_hcd_primary,
    )

    monkeypatch.setattr(
        current_drive_variables, "eta_cd_norm_ecrh", cudrivparam.eta_cd_norm_ecrh
    )

    monkeypatch.setattr(
        current_drive_variables,
        "eta_lowhyb_injector_wall_plug",
        cudrivparam.eta_lowhyb_injector_wall_plug,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "eta_hcd_primary_injector_wall_plug",
        cudrivparam.eta_hcd_primary_injector_wall_plug,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "eta_hcd_secondary_injector_wall_plug",
        cudrivparam.eta_hcd_secondary_injector_wall_plug,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "eta_ecrh_injector_wall_plug",
        cudrivparam.eta_ecrh_injector_wall_plug,
    )

    monkeypatch.setattr(
        current_drive_variables, "f_p_beam_orbit_loss", cudrivparam.f_p_beam_orbit_loss
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_injected_total_mw",
        cudrivparam.p_hcd_injected_total_mw,
    )

    monkeypatch.setattr(current_drive_variables, "pwpnb", cudrivparam.pwpnb)

    monkeypatch.setattr(
        current_drive_variables,
        "eta_beam_injector_wall_plug",
        cudrivparam.eta_beam_injector_wall_plug,
    )

    monkeypatch.setattr(current_drive_variables, "e_beam_kev", cudrivparam.e_beam_kev)

    monkeypatch.setattr(
        current_drive_variables, "eta_cd_hcd_primary", cudrivparam.eta_cd_hcd_primary
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_lowhyb_electric_mw",
        cudrivparam.p_hcd_lowhyb_electric_mw,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_ecrh_electric_mw",
        cudrivparam.p_hcd_ecrh_electric_mw,
    )

    monkeypatch.setattr(
        current_drive_variables, "p_beam_injected_mw", cudrivparam.p_beam_injected_mw
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_beam_shine_through_mw",
        cudrivparam.p_beam_shine_through_mw,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_injected_electrons_mw",
        cudrivparam.p_hcd_injected_electrons_mw,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_injected_ions_mw",
        cudrivparam.p_hcd_injected_ions_mw,
    )

    monkeypatch.setattr(current_drive_variables, "bigq", cudrivparam.bigq)

    monkeypatch.setattr(
        current_drive_variables,
        "f_c_plasma_bootstrap",
        cudrivparam.f_c_plasma_bootstrap,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "f_c_plasma_bootstrap_max",
        cudrivparam.f_c_plasma_bootstrap_max,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "n_beam_decay_lengths_core",
        cudrivparam.n_beam_decay_lengths_core,
    )

    monkeypatch.setattr(
        current_drive_variables, "p_hcd_injected_max", cudrivparam.p_hcd_injected_max
    )

    monkeypatch.setattr(
        current_drive_variables, "dx_beam_shield", cudrivparam.dx_beam_shield
    )

    monkeypatch.setattr(current_drive_variables, "frbeam", cudrivparam.frbeam)

    monkeypatch.setattr(current_drive_variables, "rtanbeam", cudrivparam.rtanbeam)

    monkeypatch.setattr(current_drive_variables, "rtanmax", cudrivparam.rtanmax)

    monkeypatch.setattr(
        current_drive_variables,
        "f_c_plasma_diamagnetic",
        cudrivparam.f_c_plasma_diamagnetic,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "f_c_plasma_pfirsch_schluter",
        cudrivparam.f_c_plasma_pfirsch_schluter,
    )

    monkeypatch.setattr(
        current_drive_variables,
        "f_c_plasma_internal",
        cudrivparam.f_c_plasma_internal,
    )

    monkeypatch.setattr(
        current_drive_variables, "n_ecrh_harmonic", cudrivparam.n_ecrh_harmonic
    )

    monkeypatch.setattr(current_drive_variables, "xi_ebw", cudrivparam.xi_ebw)

    monkeypatch.setattr(physics_variables, "dene", cudrivparam.dene)

    monkeypatch.setattr(physics_variables, "te", cudrivparam.te)

    monkeypatch.setattr(physics_variables, "rmajor", cudrivparam.rmajor)

    monkeypatch.setattr(physics_variables, "ten", cudrivparam.ten)

    monkeypatch.setattr(physics_variables, "zeff", cudrivparam.zeff)

    monkeypatch.setattr(physics_variables, "dlamee", cudrivparam.dlamee)

    monkeypatch.setattr(physics_variables, "beta", cudrivparam.beta)

    monkeypatch.setattr(physics_variables, "rhopedt", cudrivparam.rhopedt)

    monkeypatch.setattr(physics_variables, "rhopedn", cudrivparam.rhopedn)

    monkeypatch.setattr(physics_variables, "te0", cudrivparam.te0)

    monkeypatch.setattr(physics_variables, "teped", cudrivparam.teped)

    monkeypatch.setattr(physics_variables, "tesep", cudrivparam.tesep)

    monkeypatch.setattr(physics_variables, "alphat", cudrivparam.alphat)

    monkeypatch.setattr(physics_variables, "alphan", cudrivparam.alphan)

    monkeypatch.setattr(physics_variables, "ne0", cudrivparam.ne0)

    monkeypatch.setattr(physics_variables, "nesep", cudrivparam.nesep)

    monkeypatch.setattr(physics_variables, "neped", cudrivparam.neped)

    monkeypatch.setattr(physics_variables, "bt", cudrivparam.bt)

    monkeypatch.setattr(physics_variables, "rminor", cudrivparam.rminor)

    monkeypatch.setattr(physics_variables, "tbeta", cudrivparam.tbeta)

    monkeypatch.setattr(physics_variables, "plasma_current", cudrivparam.plasma_current)

    monkeypatch.setattr(physics_variables, "ipedestal", cudrivparam.ipedestal)

    monkeypatch.setattr(
        physics_variables, "f_c_plasma_auxiliary", cudrivparam.f_c_plasma_auxiliary
    )

    monkeypatch.setattr(physics_variables, "ignite", cudrivparam.ignite)

    monkeypatch.setattr(
        physics_variables, "p_plasma_ohmic_mw", cudrivparam.p_plasma_ohmic_mw
    )

    monkeypatch.setattr(physics_variables, "fusion_power", cudrivparam.fusion_power)

    monkeypatch.setattr(
        physics_variables,
        "f_c_plasma_inductive",
        cudrivparam.f_c_plasma_inductive,
    )

    monkeypatch.setattr(
        physics_variables,
        "f_c_plasma_non_inductive",
        cudrivparam.f_c_plasma_non_inductive,
    )

    monkeypatch.setattr(cost_variables, "startupratio", cudrivparam.startupratio)

    current_drive.cudriv()

    assert heat_transport_variables.p_hcd_electric_total_mw == pytest.approx(
        cudrivparam.expected_p_hcd_electric_total_mw
    )

    assert current_drive_variables.p_hcd_ecrh_injected_total_mw == pytest.approx(
        cudrivparam.expected_p_hcd_ecrh_injected_total_mw
    )

    assert current_drive_variables.eta_cd_norm_hcd_primary == pytest.approx(
        cudrivparam.expected_gamcd
    )

    assert current_drive_variables.eta_hcd_primary_injector_wall_plug == pytest.approx(
        cudrivparam.expected_etacd
    )

    assert current_drive_variables.p_hcd_injected_total_mw == pytest.approx(
        cudrivparam.expected_p_hcd_injected_total_mw
    )

    assert current_drive_variables.eta_cd_hcd_primary == pytest.approx(
        cudrivparam.expected_effcd
    )

    assert current_drive_variables.p_hcd_ecrh_electric_mw == pytest.approx(
        cudrivparam.expected_p_hcd_ecrh_electric_mw
    )

    assert current_drive_variables.p_hcd_injected_electrons_mw == pytest.approx(
        cudrivparam.expected_p_hcd_injected_electrons_mw
    )

    assert current_drive_variables.bigq == pytest.approx(cudrivparam.expected_bigq)


def test_sigbeam(current_drive):
    assert current_drive.neutral_beam.sigbeam(
        1e3, 13.07, 8.0e-1, 0.1, 1e-4, 1e-4, 1e-4
    ) == pytest.approx(2.013589662302492e-11)
