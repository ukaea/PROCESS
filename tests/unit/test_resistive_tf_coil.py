from typing import Any, NamedTuple

import pytest

from process.fortran import (
    build_variables,
    physics_variables,
    tfcoil_variables,
)
from process.resistive_tf_coil import ResistiveTFCoil


@pytest.fixture
def resistive_tf_coil():
    """Provides SuperconductingTFCoil object for testing.

    :returns: initialised SuperconductingTFCoil object
    :rtype: process.sctfcoil.SuperconductingTFCoil
    """
    return ResistiveTFCoil()


class ResTfInternalGeomParam(NamedTuple):
    n_tf_turn: Any = None

    thicndut: Any = None

    thkcas: Any = None

    dr_tf_wp: Any = None

    tftort: Any = None

    tfareain: Any = None

    c_tf_total: Any = None

    fcoolcp: Any = None

    cpttf: Any = None

    cdtfleg: Any = None

    casthi: Any = None

    aiwp: Any = None

    acasetf: Any = None

    tinstf: Any = None

    n_tf_coils: Any = None

    dr_tf_outboard: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    r_cp_top: Any = None

    itart: Any = None

    expected_n_tf_turn: Any = None

    expected_cpttf: Any = None

    expected_cdtfleg: Any = None

    expected_aiwp: Any = None

    expected_acasetf: Any = None


@pytest.mark.parametrize(
    "restfinternalgeomparam",
    (
        ResTfInternalGeomParam(
            n_tf_turn=0,
            thicndut=0.00080000000000000004,
            thkcas=0,
            dr_tf_wp=0.15483000000000002,
            tftort=0.45367650933034859,
            tfareain=0.0753112923616783,
            c_tf_total=25500000,
            fcoolcp=0.12725,
            cpttf=70000,
            cdtfleg=0,
            casthi=0.0077415000000000019,
            aiwp=0,
            acasetf=0,
            tinstf=0,
            n_tf_coils=12,
            dr_tf_outboard=0.15483000000000002,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.87643571428571443,
            itart=1,
            expected_n_tf_turn=1,
            expected_cpttf=2125000,
            expected_cdtfleg=421788350.27812088,
            expected_aiwp=0.00030678028680367151,
            expected_acasetf=0.00061190425043863676,
        ),
        ResTfInternalGeomParam(
            n_tf_turn=1,
            thicndut=0.00080000000000000004,
            thkcas=0,
            dr_tf_wp=0.14708850000000001,
            tftort=0.44435902370665786,
            tfareain=0.0753112923616783,
            c_tf_total=25500000,
            fcoolcp=0.12725,
            cpttf=2125000,
            cdtfleg=421788350.27812088,
            casthi=0.0077415000000000019,
            aiwp=0.00030678028680367151,
            acasetf=0.00061190425043863676,
            tinstf=0,
            n_tf_coils=12,
            dr_tf_outboard=0.15483000000000002,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.85843571428571441,
            itart=1,
            expected_n_tf_turn=1,
            expected_cpttf=2125000,
            expected_cdtfleg=430664525.98439038,
            expected_aiwp=0.00029439388680367086,
            expected_acasetf=0.00061190425043863676,
        ),
    ),
)
def test_res_tf_internal_geom(restfinternalgeomparam, monkeypatch, resistive_tf_coil):
    """
    Automatically generated Regression Unit Test for res_tf_internal_geom.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param restfinternalgeomparam: the data used to mock and assert in this test.
    :type restfinternalgeomparam: restfinternalgeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param sctfcoil: initialised SuperconductingTFCoil object
    :type sctfcoil: process.sctfcoil.SuperconductingTFCoil
    """

    monkeypatch.setattr(tfcoil_variables, "n_tf_turn", restfinternalgeomparam.n_tf_turn)

    monkeypatch.setattr(tfcoil_variables, "thicndut", restfinternalgeomparam.thicndut)

    monkeypatch.setattr(tfcoil_variables, "thkcas", restfinternalgeomparam.thkcas)

    monkeypatch.setattr(tfcoil_variables, "dr_tf_wp", restfinternalgeomparam.dr_tf_wp)

    monkeypatch.setattr(tfcoil_variables, "tftort", restfinternalgeomparam.tftort)

    monkeypatch.setattr(tfcoil_variables, "tfareain", restfinternalgeomparam.tfareain)

    monkeypatch.setattr(
        tfcoil_variables, "c_tf_total", restfinternalgeomparam.c_tf_total
    )

    monkeypatch.setattr(tfcoil_variables, "fcoolcp", restfinternalgeomparam.fcoolcp)

    monkeypatch.setattr(tfcoil_variables, "cpttf", restfinternalgeomparam.cpttf)

    monkeypatch.setattr(tfcoil_variables, "cdtfleg", restfinternalgeomparam.cdtfleg)

    monkeypatch.setattr(tfcoil_variables, "casthi", restfinternalgeomparam.casthi)

    monkeypatch.setattr(tfcoil_variables, "aiwp", restfinternalgeomparam.aiwp)

    monkeypatch.setattr(tfcoil_variables, "acasetf", restfinternalgeomparam.acasetf)

    monkeypatch.setattr(tfcoil_variables, "tinstf", restfinternalgeomparam.tinstf)

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_coils", restfinternalgeomparam.n_tf_coils
    )

    monkeypatch.setattr(
        build_variables, "dr_tf_outboard", restfinternalgeomparam.dr_tf_outboard
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", restfinternalgeomparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", restfinternalgeomparam.r_tf_inboard_out
    )

    monkeypatch.setattr(build_variables, "r_cp_top", restfinternalgeomparam.r_cp_top)

    monkeypatch.setattr(physics_variables, "itart", restfinternalgeomparam.itart)

    resistive_tf_coil.res_tf_internal_geom()

    assert tfcoil_variables.n_tf_turn == pytest.approx(
        restfinternalgeomparam.expected_n_tf_turn
    )

    assert tfcoil_variables.cpttf == pytest.approx(
        restfinternalgeomparam.expected_cpttf
    )

    assert tfcoil_variables.cdtfleg == pytest.approx(
        restfinternalgeomparam.expected_cdtfleg
    )

    assert tfcoil_variables.aiwp == pytest.approx(restfinternalgeomparam.expected_aiwp)

    assert tfcoil_variables.acasetf == pytest.approx(
        restfinternalgeomparam.expected_acasetf
    )

    assert tfcoil_variables.n_tf_turn == pytest.approx(
        restfinternalgeomparam.expected_n_tf_turn
    )

    assert tfcoil_variables.cpttf == pytest.approx(
        restfinternalgeomparam.expected_cpttf
    )

    assert tfcoil_variables.cdtfleg == pytest.approx(
        restfinternalgeomparam.expected_cdtfleg
    )

    assert tfcoil_variables.aiwp == pytest.approx(restfinternalgeomparam.expected_aiwp)

    assert tfcoil_variables.acasetf == pytest.approx(
        restfinternalgeomparam.expected_acasetf
    )
