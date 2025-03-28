from typing import Any, NamedTuple

import pytest

from process.build import Build
from process.fortran import (
    build_variables,
    current_drive_variables,
    divertor_variables,
    physics_variables,
    tfcoil_variables,
)


@pytest.fixture
def build():
    """Provides Build object for testing.

    :returns build: initialised Build object
    :rtype: process.build.Build
    """
    return Build()


class DivgeomParam(NamedTuple):
    rspo: Any = None

    plleno: Any = None

    tfoffset: Any = None

    plsepi: Any = None

    plleni: Any = None

    plsepo: Any = None

    betao: Any = None

    betai: Any = None

    itart: Any = None

    rmajor: Any = None

    rminor: Any = None

    idivrt: Any = None

    kappa: Any = None

    triang: Any = None

    iprint: Any = None

    outfile: Any = None

    expected_divht: Any = None


class RippleAmplitudeParam(NamedTuple):
    rminor: Any = None

    rmajor: Any = None

    tinstf: Any = None

    n_tf_coils: Any = None

    dx_tf_inboard_out_toroidal: Any = None

    casths: Any = None

    dr_tf_wp: Any = None

    thkcas: Any = None

    casths_fraction: Any = None

    i_tf_sup: Any = None

    i_tf_wp_geom: Any = None

    tfinsgap: Any = None

    tfc_sidewall_is_fraction: Any = None

    r_tf_inboard_in: Any = None

    ripmax: Any = None

    r_tf_outboard_mid: Any = None

    expected_ripple: Any = None

    expected_r_tf_outboard_midmin: Any = None

    expected_flag: Any = None


@pytest.mark.parametrize(
    "divgeomparam",
    (
        DivgeomParam(
            rspo=8.2125352340518898,
            plleno=1,
            tfoffset=0,
            plsepi=1,
            plleni=1,
            plsepo=1.5,
            betao=1,
            betai=1,
            itart=0,
            rmajor=8.8901000000000003,
            rminor=2.8677741935483869,
            idivrt=1,
            kappa=1.8480000000000001,
            triang=0.5,
            iprint=0,
            outfile=11,
            expected_divht=2.002443311884611,
        ),
        DivgeomParam(
            rspo=8.2125352340518898,
            plleno=1,
            tfoffset=0,
            plsepi=1,
            plleni=1,
            plsepo=1.5,
            betao=1,
            betai=1,
            itart=0,
            rmajor=8.8901000000000003,
            rminor=2.8677741935483869,
            idivrt=1,
            kappa=1.8480000000000001,
            triang=0.5,
            iprint=0,
            outfile=11,
            expected_divht=2.002443311884611,
        ),
    ),
)
def test_divgeom(divgeomparam, monkeypatch, build):
    """
    Automatically generated Regression Unit Test for divgeom.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param divgeomparam: the data used to mock and assert in this test.
    :type divgeomparam: divgeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param build: fixture containing an initialised `Build` object
    :type build: tests.unit.test_build.build (functional fixture)
    """

    monkeypatch.setattr(build_variables, "rspo", divgeomparam.rspo)

    monkeypatch.setattr(build_variables, "plleno", divgeomparam.plleno)

    monkeypatch.setattr(build_variables, "tfoffset", divgeomparam.tfoffset)

    monkeypatch.setattr(build_variables, "plsepi", divgeomparam.plsepi)

    monkeypatch.setattr(build_variables, "plleni", divgeomparam.plleni)

    monkeypatch.setattr(build_variables, "plsepo", divgeomparam.plsepo)

    monkeypatch.setattr(divertor_variables, "betao", divgeomparam.betao)

    monkeypatch.setattr(divertor_variables, "betai", divgeomparam.betai)

    monkeypatch.setattr(physics_variables, "itart", divgeomparam.itart)

    monkeypatch.setattr(physics_variables, "rmajor", divgeomparam.rmajor)

    monkeypatch.setattr(physics_variables, "rminor", divgeomparam.rminor)

    monkeypatch.setattr(physics_variables, "idivrt", divgeomparam.idivrt)

    monkeypatch.setattr(physics_variables, "kappa", divgeomparam.kappa)

    monkeypatch.setattr(physics_variables, "triang", divgeomparam.triang)

    divht = build.divgeom(output=False)

    assert divht == pytest.approx(divgeomparam.expected_divht)


@pytest.mark.parametrize(
    "rippleamplitudeparam",
    (
        RippleAmplitudeParam(
            rminor=2.8677741935483869,
            rmajor=8.8901000000000003,
            tinstf=0.0080000000000000019,
            n_tf_coils=16,
            dx_tf_inboard_out_toroidal=1,
            casths=0.05000000000000001,
            dr_tf_wp=0.54261087836601019,
            thkcas=0.52465000000000006,
            casths_fraction=0.059999999999999998,
            i_tf_sup=1,
            i_tf_wp_geom=0,
            tfinsgap=0.01,
            tfc_sidewall_is_fraction=False,
            r_tf_inboard_in=2.9939411851091102,
            ripmax=0.60000000000000009,
            r_tf_outboard_mid=14.988874193548387,
            expected_ripple=2.3850014198003961,
            expected_r_tf_outboard_midmin=16.519405859443332,
            expected_flag=0,
        ),
        RippleAmplitudeParam(
            rminor=2.8677741935483869,
            rmajor=8.8901000000000003,
            tinstf=0.0080000000000000019,
            n_tf_coils=16,
            dx_tf_inboard_out_toroidal=1,
            casths=0.05000000000000001,
            dr_tf_wp=0.54261087836601019,
            thkcas=0.52465000000000006,
            casths_fraction=0.059999999999999998,
            i_tf_sup=1,
            i_tf_wp_geom=0,
            tfinsgap=0.01,
            tfc_sidewall_is_fraction=False,
            r_tf_inboard_in=2.9939411851091102,
            ripmax=0.60000000000000009,
            r_tf_outboard_mid=16.519405859443332,
            expected_ripple=0.59999999999999987,
            expected_r_tf_outboard_midmin=16.519405859443332,
            expected_flag=0,
        ),
    ),
)
def test_ripple_amplitude(rippleamplitudeparam, monkeypatch, build):
    """
    Automatically generated Regression Unit Test for ripple_amplitude.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param rippleamplitudeparam: the data used to mock and assert in this test.
    :type rippleamplitudeparam: rippleamplitudeparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param build: fixture containing an initialised `Build` object
    :type build: tests.unit.test_build.build (functional fixture)
    """

    monkeypatch.setattr(physics_variables, "rminor", rippleamplitudeparam.rminor)

    monkeypatch.setattr(physics_variables, "rmajor", rippleamplitudeparam.rmajor)

    monkeypatch.setattr(tfcoil_variables, "tinstf", rippleamplitudeparam.tinstf)

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", rippleamplitudeparam.n_tf_coils)

    monkeypatch.setattr(
        tfcoil_variables,
        "dx_tf_inboard_out_toroidal",
        rippleamplitudeparam.dx_tf_inboard_out_toroidal,
    )

    monkeypatch.setattr(tfcoil_variables, "casths", rippleamplitudeparam.casths)

    monkeypatch.setattr(tfcoil_variables, "dr_tf_wp", rippleamplitudeparam.dr_tf_wp)

    monkeypatch.setattr(tfcoil_variables, "thkcas", rippleamplitudeparam.thkcas)

    monkeypatch.setattr(
        tfcoil_variables, "casths_fraction", rippleamplitudeparam.casths_fraction
    )

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", rippleamplitudeparam.i_tf_sup)

    monkeypatch.setattr(
        tfcoil_variables, "i_tf_wp_geom", rippleamplitudeparam.i_tf_wp_geom
    )

    monkeypatch.setattr(tfcoil_variables, "tfinsgap", rippleamplitudeparam.tfinsgap)

    monkeypatch.setattr(
        tfcoil_variables,
        "tfc_sidewall_is_fraction",
        rippleamplitudeparam.tfc_sidewall_is_fraction,
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", rippleamplitudeparam.r_tf_inboard_in
    )

    ripple, r_tf_outboard_midmin, flag = build.ripple_amplitude(
        ripmax=rippleamplitudeparam.ripmax,
        r_tf_outboard_mid=rippleamplitudeparam.r_tf_outboard_mid,
    )

    assert ripple == pytest.approx(rippleamplitudeparam.expected_ripple)

    assert r_tf_outboard_midmin == pytest.approx(
        rippleamplitudeparam.expected_r_tf_outboard_midmin
    )

    assert flag == pytest.approx(rippleamplitudeparam.expected_flag)


class PortszParam(NamedTuple):
    r_tf_outboard_mid: Any = None

    dr_tf_outboard: Any = None

    rtanbeam: Any = None

    rtanmax: Any = None

    nbshield: Any = None

    beamwd: Any = None

    frbeam: Any = None

    rmajor: Any = None

    dx_tf_inboard_out_toroidal: Any = None

    n_tf_coils: Any = None

    expected_rtanbeam: Any = None

    expected_rtanmax: Any = None


@pytest.mark.parametrize(
    "portszparam",
    (
        PortszParam(
            r_tf_outboard_mid=16.519405859443332,
            dr_tf_outboard=1.208,
            rtanbeam=0,
            rtanmax=0,
            nbshield=0.5,
            beamwd=0.57999999999999996,
            frbeam=1.05,
            rmajor=8.8901000000000003,
            dx_tf_inboard_out_toroidal=1.6395161177915356,
            n_tf_coils=16,
            expected_rtanbeam=9.3346050000000016,
            expected_rtanmax=14.735821603386416,
        ),
        PortszParam(
            r_tf_outboard_mid=16.519405859443332,
            dr_tf_outboard=1.208,
            rtanbeam=9.3346050000000016,
            rtanmax=14.735821603386416,
            nbshield=0.5,
            beamwd=0.57999999999999996,
            frbeam=1.05,
            rmajor=8.8901000000000003,
            dx_tf_inboard_out_toroidal=1.6395161177915356,
            n_tf_coils=16,
            expected_rtanbeam=9.3346050000000016,
            expected_rtanmax=14.735821603386416,
        ),
    ),
)
def test_portsz(portszparam, monkeypatch, build):
    """
    Automatically generated Regression Unit Test for portsz.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param portszparam: the data used to mock and assert in this test.
    :type portszparam: portszparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param build: fixture containing an initialised `Build` object
    :type build: tests.unit.test_build.build (functional fixture)
    """

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", portszparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(build_variables, "dr_tf_outboard", portszparam.dr_tf_outboard)

    monkeypatch.setattr(current_drive_variables, "rtanbeam", portszparam.rtanbeam)

    monkeypatch.setattr(current_drive_variables, "rtanmax", portszparam.rtanmax)

    monkeypatch.setattr(current_drive_variables, "nbshield", portszparam.nbshield)

    monkeypatch.setattr(current_drive_variables, "beamwd", portszparam.beamwd)

    monkeypatch.setattr(current_drive_variables, "frbeam", portszparam.frbeam)

    monkeypatch.setattr(physics_variables, "rmajor", portszparam.rmajor)

    monkeypatch.setattr(
        tfcoil_variables,
        "dx_tf_inboard_out_toroidal",
        portszparam.dx_tf_inboard_out_toroidal,
    )

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", portszparam.n_tf_coils)

    build.portsz()

    assert current_drive_variables.rtanbeam == pytest.approx(
        portszparam.expected_rtanbeam
    )

    assert current_drive_variables.rtanmax == pytest.approx(
        portszparam.expected_rtanmax
    )
