"""Unit and Integration tests for tfcoil.f90."""

from typing import Any, NamedTuple

import numpy as np
import pytest

from process.build import Build
from process.fortran import (
    build_variables,
    fwbs_variables,
    physics_variables,
    sctfcoil_module,
    tfcoil_variables,
)
from process.fortran import build_variables as bv
from process.fortran import fwbs_variables as fwbsv
from process.fortran import tfcoil_variables as tfv
from process.sctfcoil import SuperconductingTFCoil
from process.tfcoil import TFCoil


@pytest.fixture
def tfcoil():
    """Provides TFCoil object for testing.

    :param monkeypatch: pytest mocking fixture
    :type monkeypatch: MonkeyPatch

    :return tfcoil: initialised TFCoil object
    :type tfcoil: process.tfcoil.TFCoil
    """
    return TFCoil(build=Build(), sctfcoil=SuperconductingTFCoil())


@pytest.mark.parametrize(
    "temperature, expected_density",
    [(24.6, 130.02434313053487), (30.2, 113.09539723009078), (43.6, 85.26924709595201)],
)
def test_he_density(temperature, expected_density, tfcoil):
    """Tests `he_density` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param expected_density: expected result of the routine.
    :type expected_density: float

    :param tfcoil: fixture containing an initialised `TFCoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    density = tfcoil.he_density(temperature)

    assert pytest.approx(density) == expected_density


@pytest.mark.parametrize(
    "temperature, expected_cp",
    [(24.6, 5674.909063980127), (30.2, 5798.42049712345), (43.6, 5673.218322000001)],
)
def test_he_cp(temperature, expected_cp, tfcoil):
    """Tests `he_cp` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param expected_cp: expected result of the routine.
    :type expected_cp: float

    :param tfcoil: fixture containing an initialised `TFCoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """

    cp = tfcoil.he_cp(temperature)

    assert pytest.approx(cp) == expected_cp


@pytest.mark.parametrize(
    "temperature, expected_visco",
    [
        (20.6, 6.889108080243641e-06),
        (26.2, 6.859929884441028e-06),
        (43.6, 7.717393982e-06),
    ],
)
def test_he_visco(temperature, expected_visco, tfcoil):
    """Tests `he_visco` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param expected_visco: expected result of the routine.
    :type expected_visco: float

    :param tfcoil: fixture containing an initialised `TFCoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    visco = tfcoil.he_visco(temperature)

    assert pytest.approx(visco) == expected_visco


# error module needs to be initialised here because the temperature ranges are (in some cases) greater than the error conditions
@pytest.mark.parametrize(
    "temperature, expected_th_cond",
    [
        (20.6, 0.0585183573711527),
        (24.2, 0.05720100686027678),
        (43.6, 0.061189437089717184),
        (50.6, 0.06409264503),
        (54.4, 0.065706872),
    ],
)
def test_he_th_cond(temperature, expected_th_cond, reinitialise_error_module, tfcoil):
    """Tests `he_th_cond` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param expected_th_cond: expected result of the routine.
    :type expected_th_cond: float

    :param reinitialise_error_module: teardown any error side-effects
    :type reinitialise_error_module: None

    :param tfcoil: fixture containing an initialised `TFCoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    th_cond = tfcoil.he_th_cond(temperature)

    assert pytest.approx(th_cond) == expected_th_cond


@pytest.mark.parametrize(
    "temperature, expected_th_cond",
    [
        (54.4, 844.9049012800042),
        (66.9, 571.151543384937),
        (109.5, 233.66333020125035),
        (151, 250.4911087866094),
    ],
)
def test_al_th_cond(temperature, expected_th_cond, tfcoil):
    """Tests `he_th_cond` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param al_th_cond: expected result of the routine.
    :type al_th_cond: float

    :param tfcoil: fixture containing an initialised `TFCoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    th_cond = tfcoil.al_th_cond(temperature)

    assert pytest.approx(th_cond) == expected_th_cond


class CntrpstTestAsset(NamedTuple):
    """Test asset for a test case of cntrpst

    :i_tf_sup: value for tfcoil_variables.i_tf_sup to be mocked with (0=Copper, 2=Cryogenic aluminium)
    :type i_tf_sup: integer
    :tcoolin: value for tfcoil_variables.tcoolin to be mocked with (centrepost coolant inlet temperature)
    :type tcoolin: float

    :expected_dtiocool: expected value of tfcoil_variables.dtiocool after tfcoil.cntrpst routine has run
    :type expected_dtiocool: float
    :expected_tcpav2: expected value of tfcoil_variables.tcpav2 after tfcoil.cntrpst routine has run
    :type expected_tcpav2: float
    :expected_tcpmax: expected value of tfcoil_variables.tcpmax after tfcoil.cntrpst routine has run
    :type expected_tcpmax: float
    :expected_ppump: expected value of tfcoil_variables.ppump after tfcoil.cntrpst routine has run
    :type expected_ppump: float
    """

    i_tf_sup: int
    tcoolin: float
    expected_dtiocool: float
    expected_tcpav2: float
    expected_tcpmax: float
    expected_ppump: float


@pytest.mark.parametrize(
    "cntrpst_asset",
    [
        CntrpstTestAsset(
            0, 100.0, 0.00075899, 100.00109611, 100.00147829, 7.05905966e08
        ),
        # CntrpstTestAsset(1, 45.0),
        CntrpstTestAsset(
            2, 43.6, 0.00645998, 43.60678774, 43.61001841, 80926408.5501315
        ),
    ],
)
def test_cntrpst(cntrpst_asset, monkeypatch, reinitialise_error_module, tfcoil):
    """Integration test for cntrpst

    Testing tfcoil module variables being set:
        - dtiocool
        - tcpav2
        - tcpmax
        - ppump

    :param cntrpst_asset: test asset containing values to mock and expected results for the represented test case
    :type cntrpst_asset: CntrpstTestAsset

    :param monkeypatch: Mock fixture
    :type monkeypatch: object

    :param reinitialise_error_module: teardown any error side-effects
    :type reinitialise_error_module: None

    :param tfcoil: fixture containing an initialised `TFCoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    monkeypatch.setattr(tfv, "a_cp_cool", 1)
    monkeypatch.setattr(tfv, "n_tf_coils", 16)
    monkeypatch.setattr(tfv, "rcool", 0.005)
    monkeypatch.setattr(tfv, "vcool", 20.0)
    monkeypatch.setattr(tfv, "vol_cond_cp", 2)
    monkeypatch.setattr(tfv, "p_cp_resistive", 1)
    monkeypatch.setattr(tfv, "i_tf_sup", cntrpst_asset.i_tf_sup)
    monkeypatch.setattr(tfv, "tcoolin", cntrpst_asset.tcoolin)
    monkeypatch.setattr(fwbsv, "pnuc_cp_tf", 1)
    monkeypatch.setattr(bv, "hmax", 1)
    monkeypatch.setattr(bv, "dr_tf_outboard", 0.5)

    tfcoil.cntrpst()

    # appears to be the same for all cases?
    assert pytest.approx(tfv.ncool) == 203718.3271576

    assert pytest.approx(tfv.dtiocool, abs=1e-8) == cntrpst_asset.expected_dtiocool
    assert pytest.approx(tfv.tcpav2) == cntrpst_asset.expected_tcpav2
    assert pytest.approx(tfv.tcpmax) == cntrpst_asset.expected_tcpmax
    assert pytest.approx(tfv.ppump) == cntrpst_asset.expected_ppump


class TfGlobalGeometryParam(NamedTuple):
    r_tf_outboard_mid: Any = None

    r_cp_top: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    dr_tf_outboard: Any = None

    tfareain: Any = None

    c_tf_total: Any = None

    tftort: Any = None

    n_tf_coils: Any = None

    a_tf_leg_outboard: Any = None

    i_tf_sup: Any = None

    dztop: Any = None

    i_tf_case_geom: Any = None

    itart: Any = None

    kappa: Any = None

    rminor: Any = None

    h_cp_top: Any = None

    r_tf_outboard_in: Any = None

    r_tf_outboard_out: Any = None

    theta_coil: Any = None

    tan_theta_coil: Any = None

    expected_tfareain: Any = None

    expected_tftort: Any = None

    expected_a_tf_leg_outboard: Any = None

    expected_r_tf_outboard_in: Any = None

    expected_r_tf_outboard_out: Any = None

    expected_theta_coil: Any = None

    expected_tan_theta_coil: Any = None


@pytest.mark.parametrize(
    "tfglobalgeometryparam",
    (
        TfGlobalGeometryParam(
            r_tf_outboard_mid=16.519405859443332,
            r_cp_top=4.20194118510911,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            dr_tf_outboard=1.208,
            tfareain=0,
            c_tf_total=0,
            tftort=1,
            n_tf_coils=16,
            a_tf_leg_outboard=0,
            i_tf_sup=1,
            dztop=0,
            i_tf_case_geom=0,
            itart=0,
            kappa=1.8480000000000001,
            rminor=2.8677741935483869,
            h_cp_top=0,
            r_tf_outboard_in=0,
            r_tf_outboard_out=0,
            theta_coil=0,
            tan_theta_coil=0,
            expected_tfareain=27.308689677971632,
            expected_tftort=1.6395161177915356,
            expected_a_tf_leg_outboard=1.9805354702921749,
            expected_r_tf_outboard_in=15.915405859443332,
            expected_r_tf_outboard_out=17.123405859443331,
            expected_theta_coil=0.19634954084936207,
            expected_tan_theta_coil=0.19891236737965801,
        ),
        TfGlobalGeometryParam(
            r_tf_outboard_mid=17.063351291812893,
            r_cp_top=4.4822055399518357,
            r_tf_inboard_in=2.9538679176819831,
            r_tf_inboard_out=4.4822055399518357,
            dr_tf_outboard=1.5283376222698528,
            tfareain=35.703669036223495,
            c_tf_total=241812532.66279837,
            tftort=1.7488698442633552,
            n_tf_coils=16,
            a_tf_leg_outboard=2.6728635794409041,
            i_tf_sup=1,
            dztop=0,
            i_tf_case_geom=0,
            itart=0,
            kappa=1.8480000000000001,
            rminor=2.9620024998595755,
            h_cp_top=0,
            r_tf_outboard_in=16.299182480677967,
            r_tf_outboard_out=17.827520102947819,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            expected_tfareain=35.703669036223495,
            expected_tftort=1.7488698442633552,
            expected_a_tf_leg_outboard=2.6728635794409041,
            expected_r_tf_outboard_in=16.299182480677967,
            expected_r_tf_outboard_out=17.827520102947819,
            expected_theta_coil=0.19634954084936207,
            expected_tan_theta_coil=0.19891236737965801,
        ),
    ),
)
def test_tf_global_geometry(tfglobalgeometryparam, monkeypatch, tfcoil):
    """
    Automatically generated Regression Unit Test for tf_global_geometry.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfglobalgeometryparam: the data used to mock and assert in this test.
    :type tfglobalgeometryparam: tfglobalgeometryparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param sctfcoil: initialised SuperconductingTFCoil object
    :type sctfcoil: process.sctfcoil.SuperconductingTFCoil
    """

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", tfglobalgeometryparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(build_variables, "r_cp_top", tfglobalgeometryparam.r_cp_top)

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfglobalgeometryparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfglobalgeometryparam.r_tf_inboard_out
    )

    monkeypatch.setattr(
        build_variables, "dr_tf_outboard", tfglobalgeometryparam.dr_tf_outboard
    )

    monkeypatch.setattr(tfcoil_variables, "tfareain", tfglobalgeometryparam.tfareain)

    monkeypatch.setattr(
        tfcoil_variables, "c_tf_total", tfglobalgeometryparam.c_tf_total
    )

    monkeypatch.setattr(tfcoil_variables, "tftort", tfglobalgeometryparam.tftort)

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_coils", tfglobalgeometryparam.n_tf_coils
    )

    monkeypatch.setattr(
        tfcoil_variables, "a_tf_leg_outboard", tfglobalgeometryparam.a_tf_leg_outboard
    )

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tfglobalgeometryparam.i_tf_sup)

    monkeypatch.setattr(tfcoil_variables, "dztop", tfglobalgeometryparam.dztop)

    monkeypatch.setattr(
        tfcoil_variables, "i_tf_case_geom", tfglobalgeometryparam.i_tf_case_geom
    )

    monkeypatch.setattr(physics_variables, "itart", tfglobalgeometryparam.itart)

    monkeypatch.setattr(physics_variables, "kappa", tfglobalgeometryparam.kappa)

    monkeypatch.setattr(physics_variables, "rminor", tfglobalgeometryparam.rminor)

    monkeypatch.setattr(sctfcoil_module, "h_cp_top", tfglobalgeometryparam.h_cp_top)

    monkeypatch.setattr(
        sctfcoil_module, "r_tf_outboard_in", tfglobalgeometryparam.r_tf_outboard_in
    )

    monkeypatch.setattr(
        sctfcoil_module, "r_tf_outboard_out", tfglobalgeometryparam.r_tf_outboard_out
    )

    monkeypatch.setattr(sctfcoil_module, "theta_coil", tfglobalgeometryparam.theta_coil)

    monkeypatch.setattr(
        sctfcoil_module, "tan_theta_coil", tfglobalgeometryparam.tan_theta_coil
    )

    tfcoil.tf_global_geometry()

    assert tfcoil_variables.tfareain == pytest.approx(
        tfglobalgeometryparam.expected_tfareain
    )

    assert tfcoil_variables.tftort == pytest.approx(
        tfglobalgeometryparam.expected_tftort
    )

    assert tfcoil_variables.a_tf_leg_outboard == pytest.approx(
        tfglobalgeometryparam.expected_a_tf_leg_outboard
    )

    assert sctfcoil_module.r_tf_outboard_in == pytest.approx(
        tfglobalgeometryparam.expected_r_tf_outboard_in
    )

    assert sctfcoil_module.r_tf_outboard_out == pytest.approx(
        tfglobalgeometryparam.expected_r_tf_outboard_out
    )

    assert sctfcoil_module.theta_coil == pytest.approx(
        tfglobalgeometryparam.expected_theta_coil
    )

    assert sctfcoil_module.tan_theta_coil == pytest.approx(
        tfglobalgeometryparam.expected_tan_theta_coil
    )


class TfCurrentParam(NamedTuple):
    casthi: Any = None

    c_tf_total: Any = None

    rbmax: Any = None

    i_tf_sup: Any = None

    casths_fraction: Any = None

    tinstf: Any = None

    tftort: Any = None

    bmaxtf: Any = None

    tfinsgap: Any = None

    tfc_sidewall_is_fraction: Any = None

    casths: Any = None

    casthi_is_fraction: Any = None

    casthi_fraction: Any = None

    n_tf_coils: Any = None

    thicndut: Any = None

    thkcas: Any = None

    oacdcp: Any = None

    tfareain: Any = None

    r_tf_inboard_out: Any = None

    r_tf_inboard_in: Any = None

    dr_tf_inboard: Any = None

    bt: Any = None

    rmajor: Any = None

    tfc_current: Any = None

    theta_coil: Any = None

    expected_c_tf_total: Any = None

    expected_rbmax: Any = None

    expected_bmaxtf: Any = None

    expected_oacdcp: Any = None

    expected_tfc_current: Any = None


@pytest.mark.parametrize(
    "tfcurrentparam",
    (
        TfCurrentParam(
            casthi=0.060000000000000012,
            c_tf_total=0,
            rbmax=0,
            i_tf_sup=1,
            casths_fraction=0.059999999999999998,
            tinstf=0.0080000000000000019,
            tftort=1.6395161177915356,
            bmaxtf=0,
            tfinsgap=0.01,
            tfc_sidewall_is_fraction=False,
            casths=0.05000000000000001,
            casthi_is_fraction=False,
            casthi_fraction=0.050000000000000003,
            n_tf_coils=16,
            thicndut=0.002,
            thkcas=0.52465000000000006,
            oacdcp=8673900,
            tfareain=27.308689677971632,
            r_tf_inboard_out=4.20194118510911,
            r_tf_inboard_in=2.9939411851091102,
            dr_tf_inboard=1.208,
            bt=5.3292000000000002,
            rmajor=8.8901000000000003,
            tfc_current=0,
            theta_coil=0.19634954084936207,
            expected_c_tf_total=236885604.60000002,
            expected_rbmax=4.0432020634751211,
            expected_bmaxtf=11.717722779177526,
            expected_oacdcp=8674367.2945641987,
            expected_tfc_current=14805350.287500001,
        ),
    ),
)
def test_tf_current(tfcurrentparam, monkeypatch, tfcoil):
    """
    Automatically generated Regression Unit Test for tf_current.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcurrentparam: the data used to mock and assert in this test.
    :type tfcurrentparam: tfcurrentparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param sctfcoil: initialised SuperconductingTFCoil object
    :type sctfcoil: process.sctfcoil.SuperconductingTFCoil
    """

    monkeypatch.setattr(tfcoil_variables, "casthi", tfcurrentparam.casthi)

    monkeypatch.setattr(tfcoil_variables, "c_tf_total", tfcurrentparam.c_tf_total)

    monkeypatch.setattr(tfcoil_variables, "rbmax", tfcurrentparam.rbmax)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tfcurrentparam.i_tf_sup)

    monkeypatch.setattr(
        tfcoil_variables, "casths_fraction", tfcurrentparam.casths_fraction
    )

    monkeypatch.setattr(tfcoil_variables, "tinstf", tfcurrentparam.tinstf)

    monkeypatch.setattr(tfcoil_variables, "tftort", tfcurrentparam.tftort)

    monkeypatch.setattr(tfcoil_variables, "bmaxtf", tfcurrentparam.bmaxtf)

    monkeypatch.setattr(tfcoil_variables, "tfinsgap", tfcurrentparam.tfinsgap)

    monkeypatch.setattr(
        tfcoil_variables,
        "tfc_sidewall_is_fraction",
        tfcurrentparam.tfc_sidewall_is_fraction,
    )

    monkeypatch.setattr(tfcoil_variables, "casths", tfcurrentparam.casths)

    monkeypatch.setattr(
        tfcoil_variables, "casthi_is_fraction", tfcurrentparam.casthi_is_fraction
    )

    monkeypatch.setattr(
        tfcoil_variables, "casthi_fraction", tfcurrentparam.casthi_fraction
    )

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", tfcurrentparam.n_tf_coils)

    monkeypatch.setattr(tfcoil_variables, "thicndut", tfcurrentparam.thicndut)

    monkeypatch.setattr(tfcoil_variables, "thkcas", tfcurrentparam.thkcas)

    monkeypatch.setattr(tfcoil_variables, "oacdcp", tfcurrentparam.oacdcp)

    monkeypatch.setattr(tfcoil_variables, "tfareain", tfcurrentparam.tfareain)

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfcurrentparam.r_tf_inboard_out
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfcurrentparam.r_tf_inboard_in
    )

    monkeypatch.setattr(build_variables, "dr_tf_inboard", tfcurrentparam.dr_tf_inboard)

    monkeypatch.setattr(physics_variables, "bt", tfcurrentparam.bt)

    monkeypatch.setattr(physics_variables, "rmajor", tfcurrentparam.rmajor)

    monkeypatch.setattr(sctfcoil_module, "tfc_current", tfcurrentparam.tfc_current)

    monkeypatch.setattr(sctfcoil_module, "theta_coil", tfcurrentparam.theta_coil)

    tfcoil.tf_current()

    assert tfcoil_variables.c_tf_total == pytest.approx(
        tfcurrentparam.expected_c_tf_total
    )

    assert tfcoil_variables.rbmax == pytest.approx(tfcurrentparam.expected_rbmax)

    assert tfcoil_variables.bmaxtf == pytest.approx(tfcurrentparam.expected_bmaxtf)

    assert tfcoil_variables.oacdcp == pytest.approx(tfcurrentparam.expected_oacdcp)

    assert sctfcoil_module.tfc_current == pytest.approx(
        tfcurrentparam.expected_tfc_current
    )


@pytest.mark.parametrize(
    "a, b, expected_circumference",
    (
        (2.667950e9, 6.782819e8, 11464316399.111176),
        (4.7186039761812131, 3.6192586838709673, 26.308134540723429),
    ),
    ids=["johndcook", "baseline_2018"],
)
def test_circumference(a, b, expected_circumference, tfcoil):
    """Unit test for the sctfcoil circumference routine.

    This unit test uses values from an external blog referenced in the
    routine header (https://www.johndcook.com/blog/2013/05/05/ramanujan-circumference-ellipse/)
    as well as results obtained from baseline 2018.

    :param a: the value of a (x/a)^2...  in the formula of an ellipse
    :type a: float

    :param b: the value of b ...(y/b)^2...  in the formula of an ellipse
    :type b: float

    :param expected_circumference: the expected result of the routine given inputs a and b
    :type expected_circumference: float

    :param sctfcoil: initialised SuperconductingTFCoil object
    :type sctfcoil: process.sctfcoil.SuperconductingTFCoil
    """
    assert tfcoil.circumference(a, b) == pytest.approx(expected_circumference)


class TfFieldAndForceParam(NamedTuple):
    rminor: Any = None

    rmajor: Any = None

    bt: Any = None

    itart: Any = None

    r_tf_outboard_mid: Any = None

    r_vv_inboard_out: Any = None

    r_tf_inboard_mid: Any = None

    r_cp_top: Any = None

    vforce: Any = None

    n_tf_coils: Any = None

    taucq: Any = None

    sigvvall: Any = None

    cforce: Any = None

    c_tf_total: Any = None

    bmaxtf: Any = None

    i_tf_sup: Any = None

    f_vforce_inboard: Any = None

    vforce_outboard: Any = None

    tinstf: Any = None

    thicndut: Any = None

    dr_tf_wp: Any = None

    tfinsgap: Any = None

    i_cp_joints: Any = None

    casthi: Any = None

    r_tf_outboard_in: Any = None

    r_wp_inner: Any = None

    r_wp_outer: Any = None

    vforce_inboard_tot: Any = None

    expected_vforce: Any = None

    expected_cforce: Any = None

    expected_f_vforce_inboard: Any = None

    expected_vforce_outboard: Any = None

    expected_vforce_inboard_tot: Any = None


@pytest.mark.parametrize(
    "tffieldandforceparam",
    (
        TfFieldAndForceParam(
            rminor=0.97142857142857153,
            rmajor=1.7000000000000002,
            bt=3,
            itart=1,
            r_tf_outboard_mid=4.1688435714285719,
            r_vv_inboard_out=0.20483000000000001,
            r_tf_inboard_mid=0.077415000000000012,
            r_cp_top=0.87643571428571443,
            vforce=0,
            n_tf_coils=12,
            taucq=30,
            sigvvall=93000000,
            cforce=0,
            c_tf_total=25500000,
            bmaxtf=34.862617362267024,
            i_tf_sup=0,
            f_vforce_inboard=0.5,
            vforce_outboard=0,
            tinstf=0,
            thicndut=0.00080000000000000004,
            dr_tf_wp=0.15483000000000002,
            tfinsgap=0.01,
            i_cp_joints=1,
            casthi=0.0077415000000000019,
            r_tf_outboard_in=4.0914285714285716,
            r_wp_inner=0,
            r_wp_outer=0.14708850000000001,
            vforce_inboard_tot=0,
            expected_vforce=12380916.66459452,
            expected_cforce=37041530.947408713,
            expected_f_vforce_inboard=0.59539634897566385,
            expected_vforce_outboard=8413494.7991220243,
            expected_vforce_inboard_tot=148570999.97513425,
        ),
        TfFieldAndForceParam(
            rminor=0.97142857142857153,
            rmajor=1.7000000000000002,
            bt=3,
            itart=1,
            r_tf_outboard_mid=4.1868435714285717,
            r_vv_inboard_out=0.20483000000000001,
            r_tf_inboard_mid=0.077415000000000012,
            r_cp_top=0.85843571428571441,
            vforce=12380916.66459452,
            n_tf_coils=12,
            taucq=30,
            sigvvall=93000000,
            cforce=37041530.947408713,
            c_tf_total=25500000,
            bmaxtf=34.862617362267024,
            i_tf_sup=0,
            f_vforce_inboard=0.59539634897566385,
            vforce_outboard=8413494.7991220243,
            tinstf=0,
            thicndut=0.00080000000000000004,
            dr_tf_wp=0.14708850000000001,
            tfinsgap=0.01,
            i_cp_joints=1,
            casthi=0.0077415000000000019,
            r_tf_outboard_in=4.1094285714285714,
            r_wp_inner=0,
            r_wp_outer=0.14708850000000001,
            vforce_inboard_tot=148570999.97513425,
            expected_vforce=12268469.138442248,
            expected_cforce=37041530.947408713,
            expected_f_vforce_inboard=0.58932254522566518,
            expected_vforce_outboard=8549450.0771621168,
            expected_vforce_inboard_tot=147221629.66130698,
        ),
    ),
)
def test_tf_field_and_force(tffieldandforceparam, monkeypatch, tfcoil):
    """
    Automatically generated Regression Unit Test for tf_field_and_force.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param tffieldandforceparam: the data used to mock and assert in this test.
    :type tffieldandforceparam: tffieldandforceparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "rminor", tffieldandforceparam.rminor)

    monkeypatch.setattr(physics_variables, "rmajor", tffieldandforceparam.rmajor)

    monkeypatch.setattr(physics_variables, "bt", tffieldandforceparam.bt)

    monkeypatch.setattr(physics_variables, "itart", tffieldandforceparam.itart)

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", tffieldandforceparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(
        build_variables, "r_vv_inboard_out", tffieldandforceparam.r_vv_inboard_out
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_mid", tffieldandforceparam.r_tf_inboard_mid
    )

    monkeypatch.setattr(build_variables, "r_cp_top", tffieldandforceparam.r_cp_top)

    monkeypatch.setattr(tfcoil_variables, "vforce", tffieldandforceparam.vforce)

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", tffieldandforceparam.n_tf_coils)

    monkeypatch.setattr(tfcoil_variables, "cforce", tffieldandforceparam.cforce)

    monkeypatch.setattr(tfcoil_variables, "c_tf_total", tffieldandforceparam.c_tf_total)

    monkeypatch.setattr(tfcoil_variables, "bmaxtf", tffieldandforceparam.bmaxtf)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tffieldandforceparam.i_tf_sup)

    monkeypatch.setattr(
        tfcoil_variables, "f_vforce_inboard", tffieldandforceparam.f_vforce_inboard
    )

    monkeypatch.setattr(
        tfcoil_variables, "vforce_outboard", tffieldandforceparam.vforce_outboard
    )

    monkeypatch.setattr(tfcoil_variables, "tinstf", tffieldandforceparam.tinstf)

    monkeypatch.setattr(tfcoil_variables, "thicndut", tffieldandforceparam.thicndut)

    monkeypatch.setattr(tfcoil_variables, "dr_tf_wp", tffieldandforceparam.dr_tf_wp)

    monkeypatch.setattr(tfcoil_variables, "tfinsgap", tffieldandforceparam.tfinsgap)

    monkeypatch.setattr(
        tfcoil_variables, "i_cp_joints", tffieldandforceparam.i_cp_joints
    )

    monkeypatch.setattr(tfcoil_variables, "casthi", tffieldandforceparam.casthi)

    monkeypatch.setattr(
        sctfcoil_module, "r_tf_outboard_in", tffieldandforceparam.r_tf_outboard_in
    )

    monkeypatch.setattr(sctfcoil_module, "r_wp_inner", tffieldandforceparam.r_wp_inner)

    monkeypatch.setattr(sctfcoil_module, "r_wp_outer", tffieldandforceparam.r_wp_outer)

    monkeypatch.setattr(
        sctfcoil_module, "vforce_inboard_tot", tffieldandforceparam.vforce_inboard_tot
    )

    tfcoil.tf_field_and_force()

    assert tfcoil_variables.vforce == pytest.approx(
        tffieldandforceparam.expected_vforce
    )

    assert tfcoil_variables.cforce == pytest.approx(
        tffieldandforceparam.expected_cforce
    )

    assert tfcoil_variables.f_vforce_inboard == pytest.approx(
        tffieldandforceparam.expected_f_vforce_inboard
    )

    assert tfcoil_variables.vforce_outboard == pytest.approx(
        tffieldandforceparam.expected_vforce_outboard
    )

    assert sctfcoil_module.vforce_inboard_tot == pytest.approx(
        tffieldandforceparam.expected_vforce_inboard_tot
    )


class TfcindParam(NamedTuple):
    yarc: Any = None

    xarc: Any = None

    tfind: Any = None

    tfthk: Any = None

    expected_yarc: Any = None

    expected_tfind: Any = None


@pytest.mark.parametrize(
    "tfcindparam",
    (
        TfcindParam(
            yarc=np.array(
                (
                    4.5228880258064512,
                    7.5381467096774184,
                    0,
                    -9.0730900215620327,
                    -5.4438540129372193,
                ),
                order="F",
            ),
            xarc=np.array(
                (
                    4.20194118510911,
                    8.316545161290323,
                    15.915405859443332,
                    8.316545161290323,
                    4.20194118510911,
                ),
                order="F",
            ),
            tfind=0,
            tfthk=1.208,
            expected_tfind=5.4453892599192845e-06,
        ),
        TfcindParam(
            yarc=np.array(
                (
                    4.5336880258064509,
                    7.5561467096774191,
                    0,
                    -9.0730900215620327,
                    -5.4438540129372193,
                ),
                order="F",
            ),
            xarc=np.array(
                (
                    4.20194118510911,
                    8.316545161290323,
                    15.915405859443332,
                    8.316545161290323,
                    4.20194118510911,
                ),
                order="F",
            ),
            tfind=5.4524893280368181e-06,
            tfthk=1.208,
            expected_tfind=5.4524893280368181e-06,
        ),
    ),
)
def test_tfcind(tfcindparam, monkeypatch, tfcoil):
    """
    Automatically generated Regression Unit Test for tfcind.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcindparam: the data used to mock and assert in this test.
    :type tfcindparam: tfcindparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(tfcoil_variables, "tfind", tfcindparam.tfind)

    tfind = tfcoil.tfcind(
        tfthk=tfcindparam.tfthk, xarc=tfcindparam.xarc, yarc=tfcindparam.yarc
    )

    assert tfind == pytest.approx(tfcindparam.expected_tfind)


class TfCoilAreaAndMassesParam(NamedTuple):
    hr1: Any = None

    r_tf_outboard_mid: Any = None

    dr_tf_inboard: Any = None

    r_tf_inboard_mid: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    hmax: Any = None

    denstl: Any = None

    whtconsh: Any = None

    whttf: Any = None

    whtcas: Any = None

    tficrn: Any = None

    tfcryoarea: Any = None

    whtgw: Any = None

    tfocrn: Any = None

    whtconsc: Any = None

    whtconcu: Any = None

    whtcon: Any = None

    whtconin: Any = None

    vftf: Any = None

    dcond: Any = None

    dcondins: Any = None

    len_tf_coil: Any = None

    dcase: Any = None

    acndttf: Any = None

    n_tf_turn: Any = None

    n_tf_coils: Any = None

    aiwp: Any = None

    acasetfo: Any = None

    acasetf: Any = None

    fcutfsu: Any = None

    awphec: Any = None

    acstf: Any = None

    whttflgs: Any = None

    whtcp: Any = None

    whtconal: Any = None

    vol_cond_cp: Any = None

    i_tf_sup: Any = None

    i_tf_sc_mat: Any = None

    a_tf_leg_outboard: Any = None

    thkcas: Any = None

    voltfleg: Any = None

    cplen: Any = None

    itart: Any = None

    awpc: Any = None

    awptf: Any = None

    vol_ins_cp: Any = None

    vol_gr_ins_cp: Any = None

    vol_case_cp: Any = None

    a_leg_ins: Any = None

    a_leg_gr_ins: Any = None

    a_leg_cond: Any = None

    theta_coil: Any = None

    tan_theta_coil: Any = None

    expected_whtconsh: Any = None

    expected_whtcas: Any = None

    expected_tficrn: Any = None

    expected_tfcryoarea: Any = None

    expected_whtgw: Any = None

    expected_tfocrn: Any = None

    expected_whtconsc: Any = None

    expected_whtconcu: Any = None

    expected_whtcon: Any = None

    expected_whtconin: Any = None

    expected_cplen: Any = None


@pytest.mark.parametrize(
    "tfcoilareaandmassesparam",
    (
        TfCoilAreaAndMassesParam(
            hr1=0,
            r_tf_outboard_mid=16.519405859443332,
            dr_tf_inboard=1.208,
            r_tf_inboard_mid=3.5979411851091103,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            hmax=9.0730900215620327,
            denstl=7800,
            whtconsh=0,
            whttf=0,
            whtcas=0,
            tficrn=0,
            tfcryoarea=0,
            whtgw=0,
            tfocrn=0,
            whtconsc=0,
            whtconcu=0,
            whtcon=0,
            whtconin=0,
            vftf=0.30000000000000004,
            dcond=np.array(
                np.array(
                    (6080, 6080, 6070, 6080, 6080, 8500, 6070, 8500, 8500), order="F"
                ),
                order="F",
            ).transpose(),
            dcondins=1800,
            len_tf_coil=50.483843027201402,
            dcase=8000,
            acndttf=0.0014685061538103825,
            n_tf_turn=200,
            n_tf_coils=16,
            aiwp=0.087880174466980876,
            acasetfo=1.2752592893394648,
            acasetf=1.0015169239205168,
            fcutfsu=0.80884,
            awphec=0.015707963267948974,
            acstf=0.001293323051622732,
            whttflgs=0,
            whtcp=0,
            whtconal=0,
            vol_cond_cp=0,
            i_tf_sup=1,
            i_tf_sc_mat=5,
            a_tf_leg_outboard=1.9805354702921749,
            thkcas=0.52465000000000006,
            voltfleg=0,
            cplen=0,
            itart=0,
            awpc=0.70527618095271016,
            awptf=0.64024601555360383,
            vol_ins_cp=0,
            vol_gr_ins_cp=0,
            vol_case_cp=0,
            a_leg_ins=0,
            a_leg_gr_ins=0,
            a_leg_cond=0,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            expected_whtconsh=115651.90127937049,
            expected_whtcas=1034021.9996272125,
            expected_tficrn=0.8197580588957678,
            expected_tfcryoarea=6381.2092203414386,
            expected_whtgw=5909.3507916745702,
            expected_tfocrn=0.59553192892551199,
            expected_whtconsc=5802.5700395134345,
            expected_whtconcu=58744.465423173802,
            expected_whtcon=188184.68882144717,
            expected_whtconin=7985.7520793894437,
            expected_cplen=20.562180043124066,
        ),
        TfCoilAreaAndMassesParam(
            hr1=0,
            r_tf_outboard_mid=16.519405859443332,
            dr_tf_inboard=1.208,
            r_tf_inboard_mid=3.5979411851091103,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            hmax=9.0730900215620327,
            denstl=7800,
            whtconsh=115651.90127937049,
            whttf=19649856.627845347,
            whtcas=1034021.9996272125,
            tficrn=0.8197580588957678,
            tfcryoarea=6381.2092203414386,
            whtgw=5909.3507916745702,
            tfocrn=0.59553192892551199,
            whtconsc=5802.5700395134345,
            whtconcu=58744.465423173802,
            whtcon=0,
            whtconin=0,
            vftf=0.30000000000000004,
            dcond=np.array(
                np.array(
                    (6080, 6080, 6070, 6080, 6080, 8500, 6070, 8500, 8500), order="F"
                ),
                order="F",
            ).transpose(),
            dcondins=1800,
            len_tf_coil=50.514015976170839,
            dcase=8000,
            acndttf=0.0014685061538103825,
            n_tf_turn=200,
            n_tf_coils=16,
            aiwp=0.087880174466980876,
            acasetfo=1.2752592893394648,
            acasetf=1.0015169239205168,
            fcutfsu=0.80884,
            awphec=0.015707963267948974,
            acstf=0.001293323051622732,
            whttflgs=0,
            whtcp=0,
            whtconal=0,
            vol_cond_cp=0,
            i_tf_sup=1,
            i_tf_sc_mat=5,
            a_tf_leg_outboard=1.9805354702921749,
            thkcas=0.52465000000000006,
            voltfleg=0,
            cplen=20.562180043124066,
            itart=0,
            awpc=0.70527618095271016,
            awptf=0.64024601555360383,
            vol_ins_cp=0,
            vol_gr_ins_cp=0,
            vol_case_cp=0,
            a_leg_ins=0,
            a_leg_gr_ins=0,
            a_leg_cond=0,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            expected_whtconsh=115721.02357090525,
            expected_whtcas=1034699.2182961091,
            expected_tficrn=0.8197580588957678,
            expected_tfcryoarea=6385.0231118485681,
            expected_whtgw=5912.8826650262808,
            expected_tfocrn=0.59553192892551199,
            expected_whtconsc=5806.038092640837,
            expected_whtconcu=58779.575542593491,
            expected_whtcon=188297.16217276,
            expected_whtconin=7990.5249666247555,
            expected_cplen=20.562180043124066,
        ),
    ),
)
def test_tf_coil_area_and_masses(tfcoilareaandmassesparam, monkeypatch, tfcoil):
    """
    Automatically generated Regression Unit Test for tf_coil_area_and_masses.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcoilareaandmassesparam: the data used to mock and assert in this test.
    :type tfcoilareaandmassesparam: tfcoilareaandmassesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "hr1", tfcoilareaandmassesparam.hr1)

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", tfcoilareaandmassesparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(
        build_variables, "dr_tf_inboard", tfcoilareaandmassesparam.dr_tf_inboard
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_mid", tfcoilareaandmassesparam.r_tf_inboard_mid
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfcoilareaandmassesparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfcoilareaandmassesparam.r_tf_inboard_out
    )

    monkeypatch.setattr(build_variables, "hmax", tfcoilareaandmassesparam.hmax)

    monkeypatch.setattr(fwbs_variables, "denstl", tfcoilareaandmassesparam.denstl)

    monkeypatch.setattr(tfcoil_variables, "whtconsh", tfcoilareaandmassesparam.whtconsh)

    monkeypatch.setattr(tfcoil_variables, "whttf", tfcoilareaandmassesparam.whttf)

    monkeypatch.setattr(tfcoil_variables, "whtcas", tfcoilareaandmassesparam.whtcas)

    monkeypatch.setattr(tfcoil_variables, "tficrn", tfcoilareaandmassesparam.tficrn)

    monkeypatch.setattr(
        tfcoil_variables, "tfcryoarea", tfcoilareaandmassesparam.tfcryoarea
    )

    monkeypatch.setattr(tfcoil_variables, "whtgw", tfcoilareaandmassesparam.whtgw)

    monkeypatch.setattr(tfcoil_variables, "tfocrn", tfcoilareaandmassesparam.tfocrn)

    monkeypatch.setattr(tfcoil_variables, "whtconsc", tfcoilareaandmassesparam.whtconsc)

    monkeypatch.setattr(tfcoil_variables, "whtconcu", tfcoilareaandmassesparam.whtconcu)

    monkeypatch.setattr(tfcoil_variables, "whtcon", tfcoilareaandmassesparam.whtcon)

    monkeypatch.setattr(tfcoil_variables, "whtconin", tfcoilareaandmassesparam.whtconin)

    monkeypatch.setattr(
        tfcoil_variables, "tfcryoarea", tfcoilareaandmassesparam.tfcryoarea
    )

    monkeypatch.setattr(tfcoil_variables, "vftf", tfcoilareaandmassesparam.vftf)

    monkeypatch.setattr(tfcoil_variables, "dcond", tfcoilareaandmassesparam.dcond)

    monkeypatch.setattr(tfcoil_variables, "dcondins", tfcoilareaandmassesparam.dcondins)

    monkeypatch.setattr(
        tfcoil_variables, "len_tf_coil", tfcoilareaandmassesparam.len_tf_coil
    )

    monkeypatch.setattr(tfcoil_variables, "dcase", tfcoilareaandmassesparam.dcase)

    monkeypatch.setattr(tfcoil_variables, "acndttf", tfcoilareaandmassesparam.acndttf)

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_turn", tfcoilareaandmassesparam.n_tf_turn
    )

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_coils", tfcoilareaandmassesparam.n_tf_coils
    )

    monkeypatch.setattr(tfcoil_variables, "aiwp", tfcoilareaandmassesparam.aiwp)

    monkeypatch.setattr(tfcoil_variables, "acasetfo", tfcoilareaandmassesparam.acasetfo)

    monkeypatch.setattr(tfcoil_variables, "acasetf", tfcoilareaandmassesparam.acasetf)

    monkeypatch.setattr(tfcoil_variables, "fcutfsu", tfcoilareaandmassesparam.fcutfsu)

    monkeypatch.setattr(tfcoil_variables, "awphec", tfcoilareaandmassesparam.awphec)

    monkeypatch.setattr(tfcoil_variables, "acstf", tfcoilareaandmassesparam.acstf)

    monkeypatch.setattr(tfcoil_variables, "whttflgs", tfcoilareaandmassesparam.whttflgs)

    monkeypatch.setattr(tfcoil_variables, "whtcp", tfcoilareaandmassesparam.whtcp)

    monkeypatch.setattr(tfcoil_variables, "whtconal", tfcoilareaandmassesparam.whtconal)

    monkeypatch.setattr(
        tfcoil_variables, "vol_cond_cp", tfcoilareaandmassesparam.vol_cond_cp
    )

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tfcoilareaandmassesparam.i_tf_sup)

    monkeypatch.setattr(
        tfcoil_variables, "i_tf_sc_mat", tfcoilareaandmassesparam.i_tf_sc_mat
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "a_tf_leg_outboard",
        tfcoilareaandmassesparam.a_tf_leg_outboard,
    )

    monkeypatch.setattr(tfcoil_variables, "thkcas", tfcoilareaandmassesparam.thkcas)

    monkeypatch.setattr(tfcoil_variables, "voltfleg", tfcoilareaandmassesparam.voltfleg)

    monkeypatch.setattr(tfcoil_variables, "cplen", tfcoilareaandmassesparam.cplen)

    monkeypatch.setattr(physics_variables, "itart", tfcoilareaandmassesparam.itart)

    monkeypatch.setattr(sctfcoil_module, "awpc", tfcoilareaandmassesparam.awpc)

    monkeypatch.setattr(sctfcoil_module, "awptf", tfcoilareaandmassesparam.awptf)

    monkeypatch.setattr(
        sctfcoil_module, "vol_ins_cp", tfcoilareaandmassesparam.vol_ins_cp
    )

    monkeypatch.setattr(
        sctfcoil_module, "vol_gr_ins_cp", tfcoilareaandmassesparam.vol_gr_ins_cp
    )

    monkeypatch.setattr(
        sctfcoil_module, "vol_case_cp", tfcoilareaandmassesparam.vol_case_cp
    )

    monkeypatch.setattr(
        sctfcoil_module, "a_leg_ins", tfcoilareaandmassesparam.a_leg_ins
    )

    monkeypatch.setattr(
        sctfcoil_module, "a_leg_gr_ins", tfcoilareaandmassesparam.a_leg_gr_ins
    )

    monkeypatch.setattr(
        sctfcoil_module, "a_leg_cond", tfcoilareaandmassesparam.a_leg_cond
    )

    monkeypatch.setattr(
        sctfcoil_module, "theta_coil", tfcoilareaandmassesparam.theta_coil
    )

    monkeypatch.setattr(
        sctfcoil_module, "tan_theta_coil", tfcoilareaandmassesparam.tan_theta_coil
    )

    tfcoil.tf_coil_area_and_masses()

    assert tfcoil_variables.whtconsh == pytest.approx(
        tfcoilareaandmassesparam.expected_whtconsh
    )

    assert tfcoil_variables.whtcas == pytest.approx(
        tfcoilareaandmassesparam.expected_whtcas
    )

    assert tfcoil_variables.tficrn == pytest.approx(
        tfcoilareaandmassesparam.expected_tficrn
    )

    assert tfcoil_variables.tfcryoarea == pytest.approx(
        tfcoilareaandmassesparam.expected_tfcryoarea
    )

    assert tfcoil_variables.whtgw == pytest.approx(
        tfcoilareaandmassesparam.expected_whtgw
    )

    assert tfcoil_variables.tfocrn == pytest.approx(
        tfcoilareaandmassesparam.expected_tfocrn
    )

    assert tfcoil_variables.whtconsc == pytest.approx(
        tfcoilareaandmassesparam.expected_whtconsc
    )

    assert tfcoil_variables.whtconcu == pytest.approx(
        tfcoilareaandmassesparam.expected_whtconcu
    )

    assert tfcoil_variables.whtcon == pytest.approx(
        tfcoilareaandmassesparam.expected_whtcon
    )

    assert tfcoil_variables.whtconin == pytest.approx(
        tfcoilareaandmassesparam.expected_whtconin
    )

    assert tfcoil_variables.cplen == pytest.approx(
        tfcoilareaandmassesparam.expected_cplen
    )
