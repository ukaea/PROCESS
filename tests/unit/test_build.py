from typing import Any, NamedTuple

import pytest

from process.build import Build
from process.data_structure import (
    build_variables,
    divertor_variables,
    physics_variables,
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

    dz_tf_plasma_centre_offset: Any = None

    plsepi: Any = None

    plleni: Any = None

    plsepo: Any = None

    betao: Any = None

    betai: Any = None

    itart: Any = None

    rmajor: Any = None

    rminor: Any = None

    n_divertors: Any = None

    kappa: Any = None

    triang: Any = None

    iprint: Any = None

    outfile: Any = None

    expected_divht: Any = None


class RippleAmplitudeParam(NamedTuple):
    rminor: Any = None

    rmajor: Any = None

    dx_tf_wp_insulation: Any = None

    n_tf_coils: Any = None

    dx_tf_inboard_out_toroidal: Any = None

    dr_tf_wp_with_insulation: Any = None

    r_tf_wp_inboard_inner: Any = None

    r_tf_wp_inboard_centre: Any = None

    r_tf_wp_inboard_outer: Any = None

    i_tf_sup: Any = None

    i_tf_wp_geom: Any = None

    dx_tf_wp_insertion_gap: Any = None

    tfc_sidewall_is_fraction: Any = None

    ripple_b_tf_plasma_edge_max: Any = None

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
            dz_tf_plasma_centre_offset=0,
            plsepi=1,
            plleni=1,
            plsepo=1.5,
            betao=1,
            betai=1,
            itart=0,
            rmajor=8.8901000000000003,
            rminor=2.8677741935483869,
            n_divertors=1,
            kappa=1.8480000000000001,
            triang=0.5,
            iprint=0,
            outfile=11,
            expected_divht=2.002443311884611,
        ),
        DivgeomParam(
            rspo=8.2125352340518898,
            plleno=1,
            dz_tf_plasma_centre_offset=0,
            plsepi=1,
            plleni=1,
            plsepo=1.5,
            betao=1,
            betai=1,
            itart=0,
            rmajor=8.8901000000000003,
            rminor=2.8677741935483869,
            n_divertors=1,
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

    monkeypatch.setattr(
        build_variables,
        "dz_tf_plasma_centre_offset",
        divgeomparam.dz_tf_plasma_centre_offset,
    )

    monkeypatch.setattr(build_variables, "plsepi", divgeomparam.plsepi)

    monkeypatch.setattr(build_variables, "plleni", divgeomparam.plleni)

    monkeypatch.setattr(build_variables, "plsepo", divgeomparam.plsepo)

    monkeypatch.setattr(divertor_variables, "betao", divgeomparam.betao)

    monkeypatch.setattr(divertor_variables, "betai", divgeomparam.betai)

    monkeypatch.setattr(physics_variables, "itart", divgeomparam.itart)

    monkeypatch.setattr(physics_variables, "rmajor", divgeomparam.rmajor)

    monkeypatch.setattr(physics_variables, "rminor", divgeomparam.rminor)

    monkeypatch.setattr(physics_variables, "n_divertors", divgeomparam.n_divertors)

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
            dx_tf_wp_insulation=0.0080000000000000019,
            n_tf_coils=16,
            dx_tf_inboard_out_toroidal=1,
            dr_tf_wp_with_insulation=0.54261087836601019,
            r_tf_wp_inboard_inner=2.9939411851091102,
            r_tf_wp_inboard_centre=3.0939411851091102,
            r_tf_wp_inboard_outer=3.1939411851091102,
            i_tf_sup=1,
            i_tf_wp_geom=0,
            dx_tf_wp_insertion_gap=0.01,
            ripple_b_tf_plasma_edge_max=0.60000000000000009,
            r_tf_outboard_mid=14.988874193548387,
            expected_ripple=2.6693801986341232,
            expected_r_tf_outboard_midmin=16.62750502597848,
            expected_flag=0,
        ),
        RippleAmplitudeParam(
            rminor=2.8677741935483869,
            rmajor=8.8901000000000003,
            dx_tf_wp_insulation=0.0080000000000000019,
            n_tf_coils=16,
            dx_tf_inboard_out_toroidal=1,
            dr_tf_wp_with_insulation=0.54261087836601019,
            r_tf_wp_inboard_inner=0.6,
            r_tf_wp_inboard_centre=0.7,
            r_tf_wp_inboard_outer=0.8,
            i_tf_sup=1,
            i_tf_wp_geom=0,
            dx_tf_wp_insertion_gap=0.01,
            ripple_b_tf_plasma_edge_max=0.60000000000000009,
            r_tf_outboard_mid=16.519405859443332,
            expected_ripple=0.6590302350151329,
            expected_r_tf_outboard_midmin=16.62750502597848,
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

    ripple, r_tf_outboard_midmin, flag = build.plasma_outboard_edge_toroidal_ripple(
        ripple_b_tf_plasma_edge_max=rippleamplitudeparam.ripple_b_tf_plasma_edge_max,
        r_tf_outboard_mid=rippleamplitudeparam.r_tf_outboard_mid,
        n_tf_coils=rippleamplitudeparam.n_tf_coils,
        rmajor=rippleamplitudeparam.rmajor,
        rminor=rippleamplitudeparam.rminor,
        r_tf_wp_inboard_inner=0.6,
        r_tf_wp_inboard_centre=0.7,
        r_tf_wp_inboard_outer=0.8,
        dx_tf_wp_primary_toroidal=rippleamplitudeparam.dx_tf_wp_insulation,
        i_tf_sup=rippleamplitudeparam.i_tf_sup,
        i_tf_shape=0,
        dx_tf_wp_insulation=rippleamplitudeparam.dx_tf_wp_insulation,
        dx_tf_wp_insertion_gap=rippleamplitudeparam.dx_tf_wp_insertion_gap,
    )

    assert ripple == pytest.approx(rippleamplitudeparam.expected_ripple)

    assert r_tf_outboard_midmin == pytest.approx(
        rippleamplitudeparam.expected_r_tf_outboard_midmin
    )

    assert flag == pytest.approx(rippleamplitudeparam.expected_flag)


class PortszParam(NamedTuple):
    r_tf_outboard_mid: Any = None

    dr_tf_outboard: Any = None

    radius_beam_tangency: Any = None

    radius_beam_tangency_max: Any = None

    dx_beam_shield: Any = None

    dx_beam_duct: Any = None

    f_radius_beam_tangency_rmajor: Any = None

    rmajor: Any = None

    dx_tf_inboard_out_toroidal: Any = None

    n_tf_coils: Any = None

    expected_radius_beam_tangency: Any = None

    expected_radius_beam_tangency_max: Any = None


@pytest.mark.parametrize(
    "portszparam",
    (
        PortszParam(
            r_tf_outboard_mid=16.519405859443332,
            dr_tf_outboard=1.208,
            radius_beam_tangency=0,
            radius_beam_tangency_max=0,
            dx_beam_shield=0.5,
            dx_beam_duct=0.57999999999999996,
            f_radius_beam_tangency_rmajor=1.05,
            rmajor=8.8901000000000003,
            dx_tf_inboard_out_toroidal=1.6395161177915356,
            n_tf_coils=16,
            expected_radius_beam_tangency=9.3346050000000016,
            expected_radius_beam_tangency_max=14.735821603386416,
        ),
        PortszParam(
            r_tf_outboard_mid=16.519405859443332,
            dr_tf_outboard=1.208,
            radius_beam_tangency=9.3346050000000016,
            radius_beam_tangency_max=14.735821603386416,
            dx_beam_shield=0.5,
            dx_beam_duct=0.57999999999999996,
            f_radius_beam_tangency_rmajor=1.05,
            rmajor=8.8901000000000003,
            dx_tf_inboard_out_toroidal=1.6395161177915356,
            n_tf_coils=16,
            expected_radius_beam_tangency=9.3346050000000016,
            expected_radius_beam_tangency_max=14.735821603386416,
        ),
    ),
)
def test_calculate_beam_port_size(portszparam, build):
    """
    Regression Unit Test for calculate_beam_port_size with explicit inputs.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param portszparam: the data used to mock and assert in this test.
    :type portszparam: portszparam

    :param build: fixture containing an initialised `Build` object
    :type build: tests.unit.test_build.build (functional fixture)
    """

    radius_beam_tangency, radius_beam_tangency_max = build.calculate_beam_port_size(
        r_tf_outboard_mid=portszparam.r_tf_outboard_mid,
        dr_tf_outboard=portszparam.dr_tf_outboard,
        dx_beam_shield=portszparam.dx_beam_shield,
        dx_beam_duct=portszparam.dx_beam_duct,
        f_radius_beam_tangency_rmajor=portszparam.f_radius_beam_tangency_rmajor,
        rmajor=portszparam.rmajor,
        dx_tf_inboard_out_toroidal=portszparam.dx_tf_inboard_out_toroidal,
        n_tf_coils=portszparam.n_tf_coils,
    )

    assert radius_beam_tangency == pytest.approx(
        portszparam.expected_radius_beam_tangency
    )

    assert radius_beam_tangency_max == pytest.approx(
        portszparam.expected_radius_beam_tangency_max
    )
