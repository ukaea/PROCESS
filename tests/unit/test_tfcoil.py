"""Unit and Integration tests for tfcoil.f90."""

import json
from collections.abc import Sequence
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Any, NamedTuple

import numpy as np
import pytest

import process.models.tfcoil.base as tfcoil_module
from process.data_structure import (
    build_variables,
    fwbs_variables,
    superconducting_tf_coil_variables,
    tfcoil_variables,
)
from process.models.build import Build
from process.models.tfcoil.base import TFCoil


@pytest.fixture
def tfcoil():
    """Provides TFCoil object for testing.

    :param monkeypatch: pytest mocking fixture
    :type monkeypatch: MonkeyPatch

    :return tfcoil: initialised TFCoil object
    :type tfcoil: process.tfcoil.TFCoil
    """
    return TFCoil(build=Build())


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
    :temp_cp_coolant_inlet: value for tfcoil_variables.temp_cp_coolant_inlet to be mocked with (centrepost coolant inlet temperature)
    :type temp_cp_coolant_inlet: float

    :expected_dtiocool: expected value of tfcoil_variables.dtemp_cp_coolant after tfcoil.cntrpst routine has run
    :type expected_dtiocool: float
    :expected_tcpav2: expected value of tfcoil_variables.tcpav2 after tfcoil.cntrpst routine has run
    :type expected_tcpav2: float
    :expected_temp_cp_peak: expected value of tfcoil_variables.temp_cp_peak after tfcoil.cntrpst routine has run
    :type expected_temp_cp_peak: float
    :expected_ppump: expected value of tfcoil_variables.p_cp_coolant_pump_elec after tfcoil.cntrpst routine has run
    :type expected_ppump: float
    """

    i_tf_sup: int
    temp_cp_coolant_inlet: float
    expected_dtiocool: float
    expected_tcpav2: float
    expected_temp_cp_peak: float
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
        - dtemp_cp_coolant
        - tcpav2
        - temp_cp_peak
        - p_cp_coolant_pump_elec

    :param cntrpst_asset: test asset containing values to mock and expected results for the represented test case
    :type cntrpst_asset: CntrpstTestAsset

    :param monkeypatch: Mock fixture
    :type monkeypatch: object

    :param reinitialise_error_module: teardown any error side-effects
    :type reinitialise_error_module: None

    :param tfcoil: fixture containing an initialised `TFCoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    monkeypatch.setattr(tfcoil_variables, "a_cp_cool", 1)
    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", 16)
    monkeypatch.setattr(tfcoil_variables, "radius_cp_coolant_channel", 0.005)
    monkeypatch.setattr(tfcoil_variables, "vel_cp_coolant_midplane", 20.0)
    monkeypatch.setattr(tfcoil_variables, "vol_cond_cp", 2)
    monkeypatch.setattr(tfcoil_variables, "p_cp_resistive", 1)
    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", cntrpst_asset.i_tf_sup)
    monkeypatch.setattr(
        tfcoil_variables, "temp_cp_coolant_inlet", cntrpst_asset.temp_cp_coolant_inlet
    )
    monkeypatch.setattr(fwbs_variables, "pnuc_cp_tf", 1)
    monkeypatch.setattr(build_variables, "z_tf_inside_half", 1)
    monkeypatch.setattr(build_variables, "dr_tf_outboard", 0.5)

    tfcoil.cntrpst()

    # appears to be the same for all cases?
    assert pytest.approx(tfcoil_variables.n_cp_coolant_channels_total) == 203718.3271576

    assert (
        pytest.approx(tfcoil_variables.dtemp_cp_coolant, abs=1e-8)
        == cntrpst_asset.expected_dtiocool
    )
    assert pytest.approx(tfcoil_variables.tcpav2) == cntrpst_asset.expected_tcpav2
    assert (
        pytest.approx(tfcoil_variables.temp_cp_peak)
        == cntrpst_asset.expected_temp_cp_peak
    )
    assert (
        pytest.approx(tfcoil_variables.p_cp_coolant_pump_elec)
        == cntrpst_asset.expected_ppump
    )


@pytest.mark.parametrize(
    "i_tf_case_geom, i_f_dr_tf_plasma_case, f_dr_tf_plasma_case, tfc_sidewall_is_fraction, casths_fraction, n_tf_coils, dr_tf_inboard, dr_tf_nose_case, r_tf_inboard_out, r_tf_inboard_in, r_tf_outboard_mid, dr_tf_outboard, expected",
    [
        (
            0,  # Circular plasma-facing front case
            False,
            0.0,
            False,
            0.0,
            16,
            0.5,
            0.1,
            2.0,
            1.5,
            5.0,
            0.3,
            (
                pytest.approx(0.19634954084936207),  # rad_tf_coil_inboard_toroidal_half
                pytest.approx(0.198912367379658),  # tan_theta_coil
                pytest.approx(5.497787143782138),  # a_tf_inboard_total
                pytest.approx(4.85),  # r_tf_outboard_in
                pytest.approx(5.15),  # r_tf_outboard_out
                pytest.approx(0.780361288064513),  # dx_tf_inboard_out_toroidal
                pytest.approx(0.2341083864193539),  # a_tf_leg_outboard
                pytest.approx(0.03842943919353914),  # dr_tf_plasma_case
                pytest.approx(0.0),  # dx_tf_side_case_min
            ),
        ),
        (
            1,  # Straight plasma-facing front case
            True,
            0.1,
            True,
            0.05,
            12,
            0.4,
            0.05,
            1.8,
            1.4,
            4.5,
            0.25,
            (
                pytest.approx(0.2617993877991494),  # rad_tf_coil_inboard_toroidal_half
                pytest.approx(0.2679491924311227),  # tan_theta_coil
                pytest.approx(3.562478398964007),  # a_tf_inboard_total
                pytest.approx(4.375),  # r_tf_outboard_in
                pytest.approx(4.625),  # r_tf_outboard_out
                pytest.approx(0.9317485623690747),  # dx_tf_inboard_out_toroidal
                pytest.approx(0.23293714059226867),  # a_tf_leg_outboard
                pytest.approx(0.04),  # dr_tf_plasma_case
                pytest.approx(0.019426316451256392),  # dx_tf_side_case_min
            ),
        ),
    ],
)
def test_tf_global_geometry(
    tfcoil,
    i_tf_case_geom,
    i_f_dr_tf_plasma_case,
    f_dr_tf_plasma_case,
    tfc_sidewall_is_fraction,
    casths_fraction,
    n_tf_coils,
    dr_tf_inboard,
    dr_tf_nose_case,
    r_tf_inboard_out,
    r_tf_inboard_in,
    r_tf_outboard_mid,
    dr_tf_outboard,
    expected,
):
    """Test the tf_global_geometry method."""
    result = tfcoil.tf_global_geometry(
        i_tf_case_geom,
        i_f_dr_tf_plasma_case,
        f_dr_tf_plasma_case,
        tfc_sidewall_is_fraction,
        casths_fraction,
        n_tf_coils,
        dr_tf_inboard,
        dr_tf_nose_case,
        r_tf_inboard_out,
        r_tf_inboard_in,
        r_tf_outboard_mid,
        dr_tf_outboard,
    )
    assert result == expected


@pytest.mark.parametrize(
    "n_tf_coils, b_plasma_toroidal_on_axis, rmajor, r_b_tf_inboard_peak, a_tf_inboard_total, expected",
    [
        (
            16,  # Number of TF coils
            5.0,  # Toroidal magnetic field at plasma center [T]
            6.2,  # Major radius of the plasma [m]
            2.5,  # Radius at which the peak inboard B field occurs [m]
            0.8,  # Cross-sectional area of the inboard leg of the TF coil [m²]
            (
                pytest.approx(12.4),  # b_tf_inboard_peak_symmetric
                pytest.approx(154999999.93042317),  # c_tf_total
                pytest.approx(9687499.995651448),  # c_tf_coil
                pytest.approx(193749999.91302896),  # oacdcp
            ),
        ),
        (
            12,  # Number of TF coils
            3.0,  # Toroidal magnetic field at plasma center [T]
            5.0,  # Major radius of the plasma [m]
            1.8,  # Radius at which the peak inboard B field occurs [m]
            0.5,  # Cross-sectional area of the inboard leg of the TF coil [m²]
            (
                pytest.approx(8.333333333),  # b_tf_inboard_peak_symmetric
                pytest.approx(74999999.9663338),  # c_tf_total
                pytest.approx(6249999.997194484),  # c_tf_coil
                pytest.approx(149999999.9326676),  # oacdcp
            ),
        ),
    ],
)
def test_tf_current(
    tfcoil,
    n_tf_coils,
    b_plasma_toroidal_on_axis,
    rmajor,
    r_b_tf_inboard_peak,
    a_tf_inboard_total,
    expected,
):
    """Test the tf_current method."""
    result = tfcoil.tf_current(
        n_tf_coils,
        b_plasma_toroidal_on_axis,
        rmajor,
        r_b_tf_inboard_peak,
        a_tf_inboard_total,
    )
    assert result == expected


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

    b_plasma_toroidal_on_axis: Any = None

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

    b_tf_inboard_peak_symmetric: Any = None

    i_tf_sup: Any = None

    f_vforce_inboard: Any = None

    vforce_outboard: Any = None

    dx_tf_wp_insulation: Any = None

    dx_tf_turn_insulation: Any = None

    dr_tf_wp_with_insulation: Any = None

    dx_tf_wp_insertion_gap: Any = None

    i_cp_joints: Any = None

    dr_tf_plasma_case: Any = None

    r_tf_outboard_in: Any = None

    r_tf_wp_inboard_inner: Any = None

    r_tf_wp_inboard_outer: Any = None

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
            b_plasma_toroidal_on_axis=3,
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
            b_tf_inboard_peak_symmetric=34.862617362267024,
            i_tf_sup=0,
            f_vforce_inboard=0.5,
            vforce_outboard=0,
            dx_tf_wp_insulation=0,
            dx_tf_turn_insulation=0.00080000000000000004,
            dr_tf_wp_with_insulation=0.15483000000000002,
            dx_tf_wp_insertion_gap=0.01,
            i_cp_joints=1,
            dr_tf_plasma_case=0.0077415000000000019,
            r_tf_outboard_in=4.0914285714285716,
            r_tf_wp_inboard_inner=0,
            r_tf_wp_inboard_outer=0.14708850000000001,
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
            b_plasma_toroidal_on_axis=3,
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
            b_tf_inboard_peak_symmetric=34.862617362267024,
            i_tf_sup=0,
            f_vforce_inboard=0.59539634897566385,
            vforce_outboard=8413494.7991220243,
            dx_tf_wp_insulation=0,
            dx_tf_turn_insulation=0.00080000000000000004,
            dr_tf_wp_with_insulation=0.14708850000000001,
            dx_tf_wp_insertion_gap=0.01,
            i_cp_joints=1,
            dr_tf_plasma_case=0.0077415000000000019,
            r_tf_outboard_in=4.1094285714285714,
            r_tf_wp_inboard_inner=0,
            r_tf_wp_inboard_outer=0.14708850000000001,
            vforce_inboard_tot=148570999.97513425,
            expected_vforce=12268469.138442248,
            expected_cforce=37041530.947408713,
            expected_f_vforce_inboard=0.58932254522566518,
            expected_vforce_outboard=8549450.0771621168,
            expected_vforce_inboard_tot=147221629.66130698,
        ),
    ),
)
def test_tf_field_and_force(tffieldandforceparam, tfcoil):
    """
    Automatically generated Regression Unit Test for tf_field_and_force.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param tffieldandforceparam: the data used to mock and assert in this test.
    :type tffieldandforceparam: tffieldandforceparam
    """

    cforce, vforce, vforce_outboard, vforce_inboard_tot, f_vforce_inboard = (
        tfcoil.tf_field_and_force(
            i_tf_sup=tffieldandforceparam.i_tf_sup,
            r_tf_wp_inboard_outer=tffieldandforceparam.r_tf_wp_inboard_outer,
            r_tf_wp_inboard_inner=tffieldandforceparam.r_tf_wp_inboard_inner,
            r_tf_outboard_in=tffieldandforceparam.r_tf_outboard_in,
            dx_tf_wp_insulation=tffieldandforceparam.dx_tf_wp_insulation,
            dx_tf_wp_insertion_gap=tffieldandforceparam.dx_tf_wp_insertion_gap,
            b_tf_inboard_peak_symmetric=tffieldandforceparam.b_tf_inboard_peak_symmetric,
            c_tf_total=tffieldandforceparam.c_tf_total,
            n_tf_coils=tffieldandforceparam.n_tf_coils,
            dr_tf_plasma_case=tffieldandforceparam.dr_tf_plasma_case,
            rmajor=tffieldandforceparam.rmajor,
            b_plasma_toroidal_on_axis=tffieldandforceparam.b_plasma_toroidal_on_axis,
            r_cp_top=tffieldandforceparam.r_cp_top,
            itart=tffieldandforceparam.itart,
            i_cp_joints=tffieldandforceparam.i_cp_joints,
            f_vforce_inboard=tffieldandforceparam.f_vforce_inboard,
        )
    )

    assert vforce == pytest.approx(tffieldandforceparam.expected_vforce)

    assert cforce == pytest.approx(tffieldandforceparam.expected_cforce)

    assert f_vforce_inboard == pytest.approx(
        tffieldandforceparam.expected_f_vforce_inboard
    )

    assert vforce_outboard == pytest.approx(
        tffieldandforceparam.expected_vforce_outboard
    )

    assert vforce_inboard_tot == pytest.approx(
        tffieldandforceparam.expected_vforce_inboard_tot
    )


class TfcindParam(NamedTuple):
    z_tf_arc: Any
    r_tf_arc: Any
    ind_tf_coil: Any
    dr_tf_inboard: Any
    itart: Any
    i_tf_shape: Any
    z_tf_inside_half: Any
    dr_tf_outboard: Any
    r_tf_outboard_mid: Any
    r_tf_inboard_mid: Any
    expected_ind_tf_coil: Any


@pytest.mark.parametrize(
    "tfcindparam",
    (
        TfcindParam(
            z_tf_arc=np.array((
                4.5228880258064512,
                7.5381467096774184,
                0,
                -9.0730900215620327,
                -5.4438540129372193,
            )),
            r_tf_arc=np.array((
                4.20194118510911,
                8.316545161290323,
                15.915405859443332,
                8.316545161290323,
                4.20194118510911,
            )),
            ind_tf_coil=0,
            dr_tf_inboard=1.208,
            itart=0,
            i_tf_shape=1,
            # The following 4 params are not used by tf_coil_self_inductance because
            # this tests the D-shaped coil branch. However they are provided because
            # the function is Numba compiled and cannot be called with None arguments.
            z_tf_inside_half=0.0,
            dr_tf_outboard=0.0,
            r_tf_outboard_mid=0.0,
            r_tf_inboard_mid=0.0,
            expected_ind_tf_coil=5.4453892599192845e-06,
        ),
        TfcindParam(
            z_tf_arc=np.array((
                4.5336880258064509,
                7.5561467096774191,
                0,
                -9.0730900215620327,
                -5.4438540129372193,
            )),
            r_tf_arc=np.array((
                4.20194118510911,
                8.316545161290323,
                15.915405859443332,
                8.316545161290323,
                4.20194118510911,
            )),
            ind_tf_coil=5.4524893280368181e-06,
            dr_tf_inboard=1.208,
            itart=0,
            i_tf_shape=1,
            # following 4 params are not used by needed for numba to be happy
            z_tf_inside_half=0.0,
            dr_tf_outboard=0.0,
            r_tf_outboard_mid=0.0,
            r_tf_inboard_mid=0.0,
            expected_ind_tf_coil=5.4524893280368181e-06,
        ),
        TfcindParam(
            dr_tf_inboard=1.208,
            itart=0,
            i_tf_shape=0,
            z_tf_inside_half=9.0730900215620327,
            dr_tf_outboard=1.208,
            r_tf_outboard_mid=16.519405859443332,
            r_tf_inboard_mid=3.5979411851091103,
            # following 3 params are not used by needed for numba to be happy
            z_tf_arc=np.zeros(3),
            r_tf_arc=np.zeros(3),
            ind_tf_coil=0.0,
            expected_ind_tf_coil=6.26806810007207e-06,
        ),
    ),
)
def test_tf_coil_self_inductance(tfcindparam, monkeypatch, tfcoil):
    """
    Automatically generated Regression Unit Test for tf_coil_self_inductance().

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcindparam: the data used to mock and assert in this test.
    :type tfcindparam: tfcindparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(tfcoil_variables, "ind_tf_coil", tfcindparam.ind_tf_coil)

    ind_tf_coil = tfcoil.tf_coil_self_inductance(
        dr_tf_inboard=tfcindparam.dr_tf_inboard,
        r_tf_arc=tfcindparam.r_tf_arc,
        z_tf_arc=tfcindparam.z_tf_arc,
        itart=tfcindparam.itart,
        i_tf_shape=tfcindparam.i_tf_shape,
        z_tf_inside_half=tfcindparam.z_tf_inside_half,
        dr_tf_outboard=tfcindparam.dr_tf_outboard,
        r_tf_outboard_mid=tfcindparam.r_tf_outboard_mid,
        r_tf_inboard_mid=tfcindparam.r_tf_inboard_mid,
    )

    assert ind_tf_coil == pytest.approx(tfcindparam.expected_ind_tf_coil)


class TfCoilAreaAndMassesParam(NamedTuple):
    r_tf_outboard_mid: Any = None

    r_tf_inboard_mid: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    tficrn: Any = None

    tfcryoarea: Any = None

    len_tf_coil: Any = None

    tfocrn: Any = None

    rad_tf_coil_inboard_toroidal_half: Any = None

    tan_theta_coil: Any = None

    expected_tficrn: Any = None

    expected_tfcryoarea: Any = None

    expected_tfocrn: Any = None


@pytest.mark.parametrize(
    "tfcoilareaandmassesparam",
    (
        TfCoilAreaAndMassesParam(
            r_tf_outboard_mid=16.519405859443332,
            r_tf_inboard_mid=3.5979411851091103,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            len_tf_coil=50.483843027201402,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            expected_tficrn=0.8197580588957678,
            expected_tfcryoarea=6381.2092203414386,
            expected_tfocrn=0.59553192892551199,
        ),
    ),
)
def test_generic_tf_coil_area_and_masses(tfcoilareaandmassesparam, monkeypatch, tfcoil):
    """
    Automatically generated Regression Unit Test for tf_coil_area_and_masses.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcoilareaandmassesparam: the data used to mock and assert in this test.
    :type tfcoilareaandmassesparam: tfcoilareaandmassesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", tfcoilareaandmassesparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(
        tfcoil_variables, "len_tf_coil", tfcoilareaandmassesparam.len_tf_coil
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

    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "rad_tf_coil_inboard_toroidal_half",
        tfcoilareaandmassesparam.rad_tf_coil_inboard_toroidal_half,
    )
    monkeypatch.setattr(
        superconducting_tf_coil_variables,
        "tan_theta_coil",
        tfcoilareaandmassesparam.tan_theta_coil,
    )

    tfcoil.generic_tf_coil_area_and_masses()

    assert tfcoil_variables.tficrn == pytest.approx(
        tfcoilareaandmassesparam.expected_tficrn
    )

    assert tfcoil_variables.tfcryoarea == pytest.approx(
        tfcoilareaandmassesparam.expected_tfcryoarea
    )

    assert tfcoil_variables.tfocrn == pytest.approx(
        tfcoilareaandmassesparam.expected_tfocrn
    )


class StressclParam(NamedTuple):
    dr_tf_inboard: Any = None

    r_tf_inboard_mid: Any = None

    dr_bore: Any = None

    dr_cs: Any = None

    i_tf_inside_cs: Any = None

    dr_cs_tf_gap: Any = None

    z_tf_inside_half: Any = None

    r_tf_inboard_in: Any = None

    casestr: Any = None

    n_tf_coil_turns: Any = None

    dr_tf_wp_with_insulation: Any = None

    i_tf_tresca: Any = None

    a_tf_turn_cable_space_no_void: Any = None

    vforce: Any = None

    c_tf_total: Any = None

    j_tf_wp: Any = None

    sig_tf_cs_bucked: Any = None

    sig_tf_case: Any = None

    sig_tf_wp: Any = None

    dx_tf_turn_steel: Any = None

    insstrain: Any = None

    dx_tf_wp_insulation: Any = None

    dx_tf_turn_insulation: Any = None

    a_tf_turn_steel: Any = None

    dx_tf_wp_insertion_gap: Any = None

    a_tf_coil_inboard_case: Any = None

    sig_tf_case_max: Any = None

    poisson_steel: Any = None

    poisson_copper: Any = None

    poisson_al: Any = None

    n_tf_graded_layers: Any = None

    i_tf_sup: Any = None

    i_tf_bucking: Any = None

    fcoolcp: Any = None

    eyoung_cond_axial: Any = None

    eyoung_steel: Any = None

    eyoung_res_tf_buck: Any = None

    eyoung_ins: Any = None

    eyoung_al: Any = None

    eyoung_copper: Any = None

    a_tf_coil_wp_turn_insulation: Any = None

    a_tf_wp_steel: Any = None

    c_tf_turn: Any = None

    n_tf_coils: Any = None

    i_tf_stress_model: Any = None

    sig_tf_wp_max: Any = None

    i_tf_turns_integer: Any = None

    dr_tf_plasma_case: Any = None

    a_tf_wp_conductor: Any = None

    a_tf_wp_extra_void: Any = None

    a_tf_wp_coolant_channels: Any = None

    poisson_ins: Any = None

    eyoung_cond_trans: Any = None

    poisson_cond_axial: Any = None

    poisson_cond_trans: Any = None

    dia_tf_turn_coolant_channel: Any = None

    f_a_tf_turn_cable_copper: Any = None

    str_wp: Any = None

    n_tf_wp_stress_layers: Any = None

    i_pf_conductor: Any = None

    f_a_cs_turn_steel: Any = None

    f_z_cs_tf_internal: Any = None

    j_cs_flat_top_end: Any = None

    j_cs_pulse_start: Any = None

    n_pf_coils_in_group: Any = None

    c_pf_coil_turn_peak_input: Any = None

    a_tf_wp_with_insulation: Any = None

    a_tf_coil_inboard_steel: Any = None

    a_tf_coil_inboard_insulation: Any = None

    r_tf_wp_inboard_inner: Any = None

    r_tf_wp_inboard_outer: Any = None

    dx_tf_wp_toroidal_min: Any = None

    dx_tf_wp_toroidal_average: Any = None

    dx_tf_side_case_average: Any = None

    a_tf_plasma_case: Any = None

    a_tf_coil_nose_case: Any = None

    rad_tf_coil_inboard_toroidal_half: Any = None

    tan_theta_coil: Any = None

    dr_tf_turn_cable_space: Any = None

    dx_tf_turn_cable_space_average: Any = None

    vforce_inboard_tot: Any = None

    iprint: Any = None

    outfile: Any = None

    n_radial_array: Any = None

    n_tf_layer: Any = None

    expected_casestr: Any = None

    expected_sig_tf_case: Any = None

    expected_sig_tf_wp: Any = None

    expected_insstrain: Any = None

    expected_str_wp: Any = None


@pytest.mark.parametrize(
    "stressclparam",
    (
        StressclParam(
            dr_tf_inboard=1.208,
            i_tf_inside_cs=0,
            dr_cs_tf_gap=0.01,
            r_tf_inboard_mid=3.5979411851091103,
            dr_bore=2.3322000000000003,
            dr_cs=0.55242000000000002,
            z_tf_inside_half=9.0730900215620327,
            r_tf_inboard_in=2.9939411851091102,
            casestr=0,
            n_tf_coil_turns=200,
            dr_tf_wp_with_insulation=0.54261087836601019,
            i_tf_tresca=0,
            a_tf_turn_cable_space_no_void=0.001293323051622732,
            vforce=250545611.13801825,
            c_tf_total=236885604.60000002,
            j_tf_wp=23124470.793774806,
            sig_tf_cs_bucked=0,
            sig_tf_case=0,
            sig_tf_wp=0,
            dx_tf_turn_steel=0.0080000000000000002,
            insstrain=0,
            dx_tf_wp_insulation=0.0080000000000000019,
            dx_tf_turn_insulation=0.002,
            a_tf_turn_steel=0.0014685061538103825,
            dx_tf_wp_insertion_gap=0.01,
            a_tf_coil_inboard_case=1.0015169239205168,
            sig_tf_case_max=580000000,
            poisson_steel=0.29999999999999999,
            poisson_copper=0.34999999999999998,
            poisson_al=0.34999999999999998,
            n_tf_graded_layers=1,
            i_tf_sup=1,
            i_tf_bucking=1,
            fcoolcp=0.29999999999999999,
            eyoung_cond_axial=0,
            eyoung_steel=205000000000,
            eyoung_res_tf_buck=150000000000,
            eyoung_ins=20000000000,
            eyoung_al=69000000000.0,
            eyoung_copper=117000000000.0,
            a_tf_coil_wp_turn_insulation=0.087880174466980876,
            a_tf_wp_steel=0.29370123076207649,
            c_tf_turn=74026.751437500003,
            n_tf_coils=16,
            i_tf_stress_model=1,
            sig_tf_wp_max=580000000,
            i_tf_turns_integer=1,
            dr_tf_plasma_case=0.060000000000000012,
            a_tf_wp_conductor=0.1653572639592335,
            a_tf_wp_extra_void=0.07759938309736393,
            a_tf_wp_coolant_channels=0.015707963267948974,
            poisson_ins=0.34000000000000002,
            eyoung_cond_trans=0,
            poisson_cond_axial=0.30000001192092896,
            poisson_cond_trans=0.30000001192092896,
            dia_tf_turn_coolant_channel=0.010000000000000002,
            f_a_tf_turn_cable_copper=0.80884,
            str_wp=0,
            n_tf_wp_stress_layers=5,
            i_pf_conductor=0,
            f_a_cs_turn_steel=0.57874999999999999,
            f_z_cs_tf_internal=0.90000000000000002,
            j_cs_flat_top_end=20726000,
            j_cs_pulse_start=0,
            n_pf_coils_in_group=np.array((1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0)),
            c_pf_coil_turn_peak_input=np.array(
                (
                    42200,
                    42200,
                    42200,
                    42200,
                    43000,
                    43000,
                    43000,
                    43000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                    40000,
                ),
            ),
            a_tf_wp_with_insulation=0.70527618095271016,
            a_tf_coil_inboard_steel=1.2952181546825934,
            a_tf_coil_inboard_insulation=0.11646247019991701,
            r_tf_wp_inboard_inner=3.5185911851091101,
            r_tf_wp_inboard_outer=4.06120206347512,
            dx_tf_wp_toroidal_min=1.299782604942499,
            dx_tf_wp_toroidal_average=1.299782604942499,
            dx_tf_side_case_average=0.10396600719086938,
            a_tf_plasma_case=0.18607458590131154,
            a_tf_coil_nose_case=0.70261616505511615,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            dr_tf_turn_cable_space=0.030661087836601014,
            dx_tf_turn_cable_space_average=0.036389912284773368,
            vforce_inboard_tot=4008729778.208292,
            iprint=0,
            outfile=11,
            n_radial_array=100,
            n_tf_layer=3,
            expected_casestr=0.00094360452596334093,
            expected_sig_tf_case=543381805.25001633,
            expected_sig_tf_wp=397005702.35272157,
            expected_insstrain=-0.006687152422925652,
            expected_str_wp=0.0015619754370069119,
        ),
        StressclParam(
            dr_tf_inboard=1.208,
            i_tf_inside_cs=0,
            dr_cs_tf_gap=0.01,
            r_tf_inboard_mid=3.5979411851091103,
            dr_bore=2.3322000000000003,
            dr_cs=0.55242000000000002,
            z_tf_inside_half=9.0730900215620327,
            r_tf_inboard_in=2.9939411851091102,
            casestr=0.00094360452596334093,
            n_tf_coil_turns=200,
            dr_tf_wp_with_insulation=0.54261087836601019,
            i_tf_tresca=0,
            a_tf_turn_cable_space_no_void=0.001293323051622732,
            vforce=250545611.13801825,
            c_tf_total=236885604.60000002,
            j_tf_wp=23124470.793774806,
            sig_tf_cs_bucked=0,
            sig_tf_case=543381805.25001633,
            sig_tf_wp=397005702.35272157,
            dx_tf_turn_steel=0.0080000000000000002,
            insstrain=0,
            dx_tf_wp_insulation=0.0080000000000000019,
            dx_tf_turn_insulation=0.002,
            a_tf_turn_steel=0.0014685061538103825,
            dx_tf_wp_insertion_gap=0.01,
            a_tf_coil_inboard_case=1.0015169239205168,
            sig_tf_case_max=580000000,
            poisson_steel=0.29999999999999999,
            poisson_copper=0.34999999999999998,
            poisson_al=0.34999999999999998,
            n_tf_graded_layers=1,
            i_tf_sup=1,
            i_tf_bucking=1,
            fcoolcp=0.29999999999999999,
            eyoung_cond_axial=0,
            eyoung_steel=205000000000,
            eyoung_res_tf_buck=150000000000,
            eyoung_ins=20000000000,
            eyoung_al=69000000000.0,
            eyoung_copper=117000000000.0,
            a_tf_coil_wp_turn_insulation=0.087880174466980876,
            a_tf_wp_steel=0.29370123076207649,
            c_tf_turn=74026.751437500003,
            n_tf_coils=16,
            i_tf_stress_model=1,
            sig_tf_wp_max=580000000,
            i_tf_turns_integer=1,
            dr_tf_plasma_case=0.060000000000000012,
            a_tf_wp_conductor=0.1653572639592335,
            a_tf_wp_extra_void=0.07759938309736393,
            a_tf_wp_coolant_channels=0.015707963267948974,
            poisson_ins=0.34000000000000002,
            eyoung_cond_trans=0,
            poisson_cond_axial=0.30000001192092896,
            poisson_cond_trans=0.30000001192092896,
            dia_tf_turn_coolant_channel=0.010000000000000002,
            f_a_tf_turn_cable_copper=0.80884,
            str_wp=0.0015619754370069119,
            n_tf_wp_stress_layers=5,
            i_pf_conductor=0,
            f_a_cs_turn_steel=0.57874999999999999,
            f_z_cs_tf_internal=0.90000000000000002,
            j_cs_flat_top_end=20726000,
            j_cs_pulse_start=19311657.760000002,
            n_pf_coils_in_group=np.array((1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0)),
            c_pf_coil_turn_peak_input=np.array((
                42200,
                42200,
                42200,
                42200,
                43000,
                43000,
                43000,
                43000,
                40000,
                40000,
                40000,
                40000,
                40000,
                40000,
                40000,
                40000,
                40000,
                40000,
                40000,
                40000,
                40000,
                40000,
            )),
            a_tf_wp_with_insulation=0.70527618095271016,
            a_tf_coil_inboard_steel=1.2952181546825934,
            a_tf_coil_inboard_insulation=0.11646247019991701,
            r_tf_wp_inboard_inner=3.5185911851091101,
            r_tf_wp_inboard_outer=4.06120206347512,
            dx_tf_wp_toroidal_min=1.299782604942499,
            dx_tf_wp_toroidal_average=1.299782604942499,
            dx_tf_side_case_average=0.10396600719086938,
            a_tf_plasma_case=0.18607458590131154,
            a_tf_coil_nose_case=0.70261616505511615,
            rad_tf_coil_inboard_toroidal_half=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            dr_tf_turn_cable_space=0.030661087836601014,
            dx_tf_turn_cable_space_average=0.036389912284773368,
            vforce_inboard_tot=4008729778.208292,
            iprint=0,
            outfile=11,
            n_radial_array=100,
            n_tf_layer=3,
            expected_casestr=0.00094360452596334093,
            expected_sig_tf_case=543381805.25001633,
            expected_sig_tf_wp=397005702.35272157,
            expected_insstrain=-0.006687152422925652,
            expected_str_wp=0.0015619754370069119,
        ),
    ),
)
def test_stresscl(stressclparam, monkeypatch, tfcoil):
    """
    Automatically generated Regression Unit Test for stresscl.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param stressclparam: the data used to mock and assert in this test.
    :type stressclparam: stressclparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    (
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        sig_tf_wp,
        sig_tf_case,
        _,
        str_wp,
        casestr,
        insstrain,
        _,
    ) = tfcoil.stresscl(
        stressclparam.n_tf_layer,
        stressclparam.n_radial_array,
        stressclparam.n_tf_wp_stress_layers,
        stressclparam.i_tf_bucking,
        stressclparam.r_tf_inboard_in,
        stressclparam.dr_bore,
        stressclparam.z_tf_inside_half,
        stressclparam.f_z_cs_tf_internal,
        stressclparam.dr_cs,
        stressclparam.i_tf_inside_cs,
        stressclparam.dr_tf_inboard,
        stressclparam.dr_cs_tf_gap,
        stressclparam.i_pf_conductor,
        stressclparam.j_cs_flat_top_end,
        stressclparam.j_cs_pulse_start,
        stressclparam.c_pf_coil_turn_peak_input,
        stressclparam.n_pf_coils_in_group,
        70 / 22,
        3e-3,
        stressclparam.f_a_cs_turn_steel,
        stressclparam.eyoung_steel,
        stressclparam.poisson_steel,
        stressclparam.eyoung_cond_axial,
        stressclparam.poisson_cond_axial,
        stressclparam.eyoung_cond_trans,
        stressclparam.poisson_cond_trans,
        stressclparam.eyoung_ins,
        stressclparam.poisson_ins,
        stressclparam.dx_tf_turn_insulation,
        stressclparam.eyoung_copper,
        stressclparam.poisson_copper,
        stressclparam.i_tf_sup,
        stressclparam.eyoung_res_tf_buck,
        stressclparam.r_tf_wp_inboard_inner,
        stressclparam.tan_theta_coil,
        stressclparam.rad_tf_coil_inboard_toroidal_half,
        stressclparam.r_tf_wp_inboard_outer,
        stressclparam.a_tf_coil_inboard_steel,
        stressclparam.a_tf_plasma_case,
        stressclparam.a_tf_coil_nose_case,
        stressclparam.dx_tf_wp_insertion_gap,
        stressclparam.dx_tf_wp_insulation,
        stressclparam.n_tf_coil_turns,
        stressclparam.i_tf_turns_integer,
        stressclparam.dx_tf_turn_cable_space_average,
        stressclparam.dr_tf_turn_cable_space,
        stressclparam.dia_tf_turn_coolant_channel,
        stressclparam.f_a_tf_turn_cable_copper,
        stressclparam.dx_tf_turn_steel,
        stressclparam.dx_tf_side_case_average,
        stressclparam.dx_tf_wp_toroidal_average,
        stressclparam.a_tf_coil_inboard_insulation,
        stressclparam.a_tf_wp_steel,
        stressclparam.a_tf_wp_conductor,
        stressclparam.a_tf_wp_with_insulation,
        stressclparam.eyoung_al,
        stressclparam.poisson_al,
        stressclparam.fcoolcp,
        stressclparam.n_tf_graded_layers,
        stressclparam.c_tf_total,
        stressclparam.dr_tf_plasma_case,
        stressclparam.i_tf_stress_model,
        stressclparam.vforce_inboard_tot,
        stressclparam.i_tf_tresca,
        stressclparam.a_tf_coil_inboard_case,
        stressclparam.vforce,
        stressclparam.a_tf_turn_steel,
    )

    assert casestr == pytest.approx(stressclparam.expected_casestr, rel=0.01)

    assert sig_tf_case == pytest.approx(stressclparam.expected_sig_tf_case, rel=0.01)

    assert sig_tf_wp == pytest.approx(stressclparam.expected_sig_tf_wp, rel=0.01)

    assert insstrain == pytest.approx(stressclparam.expected_insstrain, rel=0.01)

    assert str_wp == pytest.approx(stressclparam.expected_str_wp, rel=0.01)


class PlaneStressParam(NamedTuple):
    linesolv: Any = None

    n_radial_array: Any = None

    nlayers: Any = None

    nu: Any = None

    rad: Any = None

    ey: Any = None

    j: Any = None

    expected_data_key: str = ""


@dataclass
class PlaneStressExpectedData:
    expected_sigr: Sequence[float]

    expected_sigt: Sequence[float]

    expected_r_deflect: Sequence[float]

    expected_rradius: Sequence[float]

    def __post_init__(self):
        for f in fields(self):
            setattr(self, f.name, np.asarray(getattr(self, f.name), dtype=float))


@pytest.mark.parametrize(
    "planestressparam",
    (
        PlaneStressParam(
            n_radial_array=100,
            nlayers=3,
            nu=np.array(
                (0.29999999999999999, 0.30904421667064924, 0.29999999999999999),
            ),
            rad=np.array((
                2.9939411851091102,
                3.5414797139565706,
                4.0876202904571599,
                4.1476202904571595,
            )),
            ey=np.array((205000000000, 43126670035.025253, 205000000000)),
            j=np.array((0, 18097185.781970859, 0)),
            expected_data_key="test1",
        ),
        PlaneStressParam(
            n_radial_array=100,
            nlayers=3,
            nu=np.array(
                np.array(
                    (0.29999999999999999, 0.30904421667064924, 0.29999999999999999),
                )
            ),
            rad=np.array((
                2.9939411851091102,
                3.5414797139565706,
                4.0876202904571599,
                4.1476202904571595,
            )),
            ey=np.array((205000000000, 43126670035.025253, 205000000000)),
            j=np.array((0, 18097185.781970859, 0)),
            expected_data_key="test2",
        ),
        PlaneStressParam(
            linesolv=None,
            n_radial_array=100,
            nlayers=3,
            nu=np.array([0.3, 0.34006912702297704, 0.3]),
            rad=np.array([
                3.6732023601326333,
                3.7688101124061717,
                3.7649909451102674,
                3.8249909451102675,
            ]),
            ey=np.array([2.05000000e11, 21085960915.80571, 2.05000000e11]),
            j=np.array([0.00000000e00, -2245759961.294637, 0.00000000e00]),
            expected_data_key="test3",
        ),
    ),
)
def test_plane_stress(planestressparam, skip_if_incompatible_system, request):
    """
    Automatically generated Regression Unit Test for plane_stress.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param planestressparam: the data used to mock and assert in this test.
    :type planestressparam: planestressparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    sigr, sigt, r_deflect, rradius = tfcoil_module.plane_stress(
        n_radial_array=planestressparam.n_radial_array,
        nlayers=planestressparam.nlayers,
        nu=planestressparam.nu,
        rad=planestressparam.rad,
        ey=planestressparam.ey,
        j=planestressparam.j,
    )
    with open(Path(__file__).parent / "tf_plane_stress_expected_data.json") as data_file:
        result = PlaneStressExpectedData(
            **json.load(data_file)[planestressparam.expected_data_key]
        )

    assert sigr == pytest.approx(result.expected_sigr)
    assert sigt == pytest.approx(result.expected_sigt)
    assert r_deflect == pytest.approx(result.expected_r_deflect)
    assert rradius == pytest.approx(result.expected_rradius)


class ExtendedPlaneStrainParam(NamedTuple):
    n_radial_array: Any = None

    nlayers: Any = None

    i_tf_bucking: Any = None

    nu_t: Any = None

    nu_zt: Any = None

    ey_t: Any = None

    ey_z: Any = None

    d_curr: Any = None

    rad: Any = None

    v_force: Any = None

    expected_data_key: str = ""

    nan_init: list[str] | None = None


@dataclass
class ExtendedPlaneStressExpectedData:
    expected_sigr: Sequence[float]

    expected_sigt: Sequence[float]

    expected_sigz: Sequence[float]

    expected_str_r: Sequence[float]

    expected_str_t: Sequence[float]

    expected_str_z: Sequence[float]

    expected_r_deflect: Sequence[float]

    expected_rradius: Sequence[float]

    def __post_init__(self):
        for f in fields(self):
            setattr(self, f.name, np.asarray(getattr(self, f.name), dtype=float))

    def insert_0th_element(self, keys):
        if keys is None:
            return
        for f in fields(self):
            for k in keys:
                if f.name.endswith(k):
                    setattr(
                        self,
                        f.name,
                        np.asarray([None, *getattr(self, f.name)], dtype=float),
                    )


@pytest.mark.parametrize(
    "extendedplanestrainparam",
    (
        ExtendedPlaneStrainParam(
            n_radial_array=100,
            nlayers=2,
            i_tf_bucking=0,
            nu_t=np.array((0.34999999999999998, 0.29999999999999999)),
            nu_zt=np.array((0.34948024015688117, 0.29999999999999999)),
            ey_t=np.array((117000000000, 205000000000)),
            ey_z=np.array((97843910970.178864, 205000000000.00003)),
            d_curr=np.array((375174117.44492483, 0)),
            rad=np.array((0, 0.14708850000000001, 0.15483000000000002)),
            v_force=147221629.66130698,
            expected_data_key="test1",
            nan_init=["sigr", "sigt", "sigz", "str_r", "str_t", "r_deflect"],
        ),
        ExtendedPlaneStrainParam(
            n_radial_array=100,
            nlayers=5,
            i_tf_bucking=3,
            nu_t=np.array((
                0.30000000502169133,
                0.34000000000000002,
                0.29999999999999999,
                0.30901178507421895,
                0.29999999999999999,
            )),
            nu_zt=np.array((
                0.31163570564277626,
                0.34000000000000002,
                0.29999999999999999,
                0.31377709779186291,
                0.29999999999999999,
            )),
            ey_t=np.array((
                118643750000,
                2500000000,
                205000000000,
                43163597776.087654,
                205000000000,
            )),
            ey_z=np.array(
                (
                    48005309351.198608,
                    2500000000,
                    205000000000,
                    124208626934.75433,
                    390554854819.81116,
                ),
            ),
            d_curr=np.array((0, 0, 0, 18343613.563061949, 0)),
            rad=np.array((
                2.3322000000000003,
                2.8846200000000004,
                2.9346200000000002,
                3.4817726429672304,
                4.0290604740948242,
                4.0890604740948238,
            )),
            v_force=4051971733.3410816,
            expected_data_key="test2",
        ),
    ),
)
def test_extended_plane_strain(extendedplanestrainparam):
    """
    Automatically generated Regression Unit Test for extended_plane_strain.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param extendedplanestrainparam: the data used to mock and assert in this test.
    :type extendedplanestrainparam: extendedplanestrainparam
    """
    (
        rradius,
        sigr,
        sigt,
        sigz,
        str_r,
        str_t,
        str_z,
        r_deflect,
    ) = tfcoil_module.extended_plane_strain(
        n_radial_array=extendedplanestrainparam.n_radial_array,
        nlayers=extendedplanestrainparam.nlayers,
        i_tf_bucking=extendedplanestrainparam.i_tf_bucking,
        nu_t=extendedplanestrainparam.nu_t,
        nu_zt=extendedplanestrainparam.nu_zt,
        ey_t=extendedplanestrainparam.ey_t,
        ey_z=extendedplanestrainparam.ey_z,
        d_curr=extendedplanestrainparam.d_curr,
        rad=extendedplanestrainparam.rad,
        v_force=extendedplanestrainparam.v_force,
    )

    with open(
        Path(__file__).parent / "tf_extended_plane_stress_expected_data.json"
    ) as data_file:
        data = json.load(data_file)

    result = ExtendedPlaneStressExpectedData(
        **data[extendedplanestrainparam.expected_data_key]
    )

    result.insert_0th_element(extendedplanestrainparam.nan_init)

    # assert sigr == pytest.approx(result.expected_sigr, rel=0.01)

    np.testing.assert_array_almost_equal(sigr, np.array(result.expected_sigr), decimal=3)

    np.testing.assert_array_almost_equal(sigt, np.array(result.expected_sigt), decimal=3)

    np.testing.assert_array_almost_equal(sigz, np.array(result.expected_sigz), decimal=3)

    np.testing.assert_array_almost_equal(str_r, np.array(result.expected_str_r))

    np.testing.assert_array_almost_equal(str_t, np.array(result.expected_str_t))

    np.testing.assert_array_almost_equal(str_z, np.array(result.expected_str_z))

    np.testing.assert_array_almost_equal(r_deflect, np.array(result.expected_r_deflect))

    np.testing.assert_array_almost_equal(rradius, np.array(result.expected_rradius))


class EyoungParallelParam(NamedTuple):
    eyoung_j_1: Any = None

    eyoung_j_2: Any = None

    a_1: Any = None

    a_2: Any = None

    poisson_j_perp_1: Any = None

    poisson_j_perp_2: Any = None

    expected_eyoung_j_3: Any = None

    expected_a_3: Any = None

    expected_poisson_j_perp_3: Any = None


@pytest.mark.parametrize(
    "eyoungparallelparam",
    (
        EyoungParallelParam(
            eyoung_j_1=0,
            eyoung_j_2=0,
            a_1=0.010000000000000002,
            a_2=0,
            poisson_j_perp_1=0.30000001192092896,
            poisson_j_perp_2=0,
            expected_eyoung_j_3=0,
            expected_a_3=0.010000000000000002,
            expected_poisson_j_perp_3=0.30000001192092896,
        ),
        EyoungParallelParam(
            eyoung_j_1=0,
            eyoung_j_2=0,
            a_1=0.020661087836601012,
            a_2=0.010000000000000002,
            poisson_j_perp_1=0.30000001192092896,
            poisson_j_perp_2=0.30000001192092896,
            expected_eyoung_j_3=0,
            expected_a_3=0.030661087836601014,
            expected_poisson_j_perp_3=0.30000001192092896,
        ),
    ),
)
def test_eyoung_parallel(eyoungparallelparam, monkeypatch):
    """
    Automatically generated Regression Unit Test for eyoung_parallel.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param eyoungparallelparam: the data used to mock and assert in this test.
    :type eyoungparallelparam: eyoungparallelparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    eyoung_j_3, a_3, poisson_j_perp_3 = tfcoil_module.eyoung_parallel(
        eyoung_j_1=eyoungparallelparam.eyoung_j_1,
        eyoung_j_2=eyoungparallelparam.eyoung_j_2,
        a_1=eyoungparallelparam.a_1,
        a_2=eyoungparallelparam.a_2,
        poisson_j_perp_1=eyoungparallelparam.poisson_j_perp_1,
        poisson_j_perp_2=eyoungparallelparam.poisson_j_perp_2,
    )

    assert eyoung_j_3 == pytest.approx(eyoungparallelparam.expected_eyoung_j_3)

    assert a_3 == pytest.approx(eyoungparallelparam.expected_a_3)

    assert poisson_j_perp_3 == pytest.approx(
        eyoungparallelparam.expected_poisson_j_perp_3
    )


class EyoungTNestedSquaresParam(NamedTuple):
    n: Any = None

    eyoung_j_in: Any = None

    l_in: Any = None

    poisson_j_perp_in: Any = None

    expected_eyoung_j_out: Any = None

    expected_l_out: Any = None

    expected_poisson_j_perp_out: Any = None

    expected_eyoung_stiffest: Any = None


class EyoungSeriesParam(NamedTuple):
    eyoung_j_1: Any = None

    eyoung_j_2: Any = None

    l_1: Any = None

    l_2: Any = None

    poisson_j_perp_1: Any = None

    poisson_j_perp_2: Any = None

    expected_eyoung_j_3: Any = None

    expected_l_3: Any = None

    expected_poisson_j_perp_3: Any = None


class EyoungParallelArrayParam(NamedTuple):
    n: Any = None

    eyoung_j_in: Any = None

    a_in: Any = None

    poisson_j_perp_in: Any = None

    expected_eyoung_j_out: Any = None

    expected_a_out: Any = None

    expected_poisson_j_perp_out: Any = None


@pytest.mark.parametrize(
    "eyoungtnestedsquaresparam",
    (
        EyoungTNestedSquaresParam(
            n=4,
            eyoung_j_in=np.array((0, 0, 205000000000, 20000000000, 0)),
            l_in=np.array((
                0.010000000000000002,
                0.020661087836601012,
                0.016,
                0.0041799999999999997,
                0,
            )),
            poisson_j_perp_in=np.array((
                0.29999999999999999,
                0.30000001192092896,
                0.29999999999999999,
                0.34000000000000002,
                0,
            )),
            expected_eyoung_j_out=38289891367.115105,
            expected_l_out=0.050841087836601018,
            expected_poisson_j_perp_out=0.30931445806415137,
            expected_eyoung_stiffest=116443733140.5881,
        ),
        EyoungTNestedSquaresParam(
            n=4,
            eyoung_j_in=np.array((0, 0, 205000000000, 20000000000, 0)),
            l_in=np.array((
                0.010000000000000002,
                0.020661087836601012,
                0.016,
                0.0041799999999999997,
                0,
            )),
            poisson_j_perp_in=np.array((
                0.29999999999999999,
                0.30000001192092896,
                0.29999999999999999,
                0.34000000000000002,
                0,
            )),
            expected_eyoung_j_out=38289891367.115105,
            expected_l_out=0.050841087836601018,
            expected_poisson_j_perp_out=0.30931445806415137,
            expected_eyoung_stiffest=116443733140.5881,
        ),
    ),
)
def test_eyoung_t_nested_squares(eyoungtnestedsquaresparam, monkeypatch):
    """
    Automatically generated Regression Unit Test for eyoung_t_nested_squares.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param eyoungtnestedsquaresparam: the data used to mock and assert in this test.
    :type eyoungtnestedsquaresparam: eyoungtnestedsquaresparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    (
        _eyoung_j_out,
        _l_out,
        _poisson_j_perp_out,
        eyoung_stiffest,
    ) = tfcoil_module.eyoung_t_nested_squares(
        n=eyoungtnestedsquaresparam.n,
        eyoung_j_in=eyoungtnestedsquaresparam.eyoung_j_in,
        l_in=eyoungtnestedsquaresparam.l_in,
        poisson_j_perp_in=eyoungtnestedsquaresparam.poisson_j_perp_in,
    )

    # assert eyoung_j_out == pytest.approx(
    #     eyoungtnestedsquaresparam.expected_eyoung_j_out
    # )

    # assert l_out == pytest.approx(eyoungtnestedsquaresparam.expected_l_out)

    # assert poisson_j_perp_out == pytest.approx(
    #     eyoungtnestedsquaresparam.expected_poisson_j_perp_out
    # )

    assert eyoung_stiffest == pytest.approx(
        eyoungtnestedsquaresparam.expected_eyoung_stiffest
    )


@pytest.mark.parametrize(
    "eyoungseriesparam",
    (
        EyoungSeriesParam(
            eyoung_j_1=0,
            eyoung_j_2=117000000000,
            l_1=0.003949573550844649,
            l_2=0.016711514285756363,
            poisson_j_perp_1=0.30000001192092896,
            poisson_j_perp_2=0.34999999999999998,
            expected_eyoung_j_3=0,
            expected_l_3=0.020661087836601012,
            expected_poisson_j_perp_3=0.30000001192092896,
        ),
        EyoungSeriesParam(
            eyoung_j_1=0,
            eyoung_j_2=0,
            l_1=0.020661087836601012,
            l_2=0.010000000000000002,
            poisson_j_perp_1=0.30000001192092896,
            poisson_j_perp_2=0.29999999999999999,
            expected_eyoung_j_3=0,
            expected_l_3=0.030661087836601014,
            expected_poisson_j_perp_3=0.30000001192092896,
        ),
    ),
)
def test_eyoung_series(eyoungseriesparam, monkeypatch):
    """
    Automatically generated Regression Unit Test for eyoung_series.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param eyoungseriesparam: the data used to mock and assert in this test.
    :type eyoungseriesparam: eyoungseriesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    eyoung_j_3, l_3, poisson_j_perp_3 = tfcoil_module.eyoung_series(
        eyoung_j_1=eyoungseriesparam.eyoung_j_1,
        eyoung_j_2=eyoungseriesparam.eyoung_j_2,
        l_1=eyoungseriesparam.l_1,
        l_2=eyoungseriesparam.l_2,
        poisson_j_perp_1=eyoungseriesparam.poisson_j_perp_1,
        poisson_j_perp_2=eyoungseriesparam.poisson_j_perp_2,
    )

    assert eyoung_j_3 == pytest.approx(eyoungseriesparam.expected_eyoung_j_3)

    assert l_3 == pytest.approx(eyoungseriesparam.expected_l_3)

    assert poisson_j_perp_3 == pytest.approx(eyoungseriesparam.expected_poisson_j_perp_3)


@pytest.mark.parametrize(
    "eyoungparallelarrayparam",
    (
        EyoungParallelArrayParam(
            n=5,
            eyoung_j_in=np.array((205000000000, 20000000000, 117000000000, 0, 0)),
            a_in=np.array((
                0.29370123076207649,
                0.11646247019991701,
                0.13374756938078641,
                0.031609694578447076,
                0.1297552160314831,
            )),
            poisson_j_perp_in=np.array((
                0.29999999999999999,
                0.34000000000000002,
                0.34999999999999998,
                0.30000001192092896,
                0.29999999999999999,
            )),
            expected_eyoung_j_out=110859361820.72557,
            expected_a_out=0.70527618095271016,
            expected_poisson_j_perp_out=0.31608714140682664,
        ),
        EyoungParallelArrayParam(
            n=5,
            eyoung_j_in=np.array((205000000000, 20000000000, 117000000000, 0, 0)),
            a_in=np.array((
                0.29370123076207649,
                0.11646247019991701,
                0.13374756938078641,
                0.031609694578447076,
                0.1297552160314831,
            )),
            poisson_j_perp_in=np.array((
                0.29999999999999999,
                0.34000000000000002,
                0.34999999999999998,
                0.30000001192092896,
                0.29999999999999999,
            )),
            expected_eyoung_j_out=110859361820.72557,
            expected_a_out=0.70527618095271016,
            expected_poisson_j_perp_out=0.31608714140682664,
        ),
    ),
)
def test_eyoung_parallel_array(eyoungparallelarrayparam, monkeypatch):
    """
    Automatically generated Regression Unit Test for eyoung_parallel_array.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param eyoungparallelarrayparam: the data used to mock and assert in this test.
    :type eyoungparallelarrayparam: eyoungparallelarrayparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    eyoung_j_out, a_out, poisson_j_perp_out = tfcoil_module.eyoung_parallel_array(
        n=eyoungparallelarrayparam.n,
        eyoung_j_in=eyoungparallelarrayparam.eyoung_j_in,
        a_in=eyoungparallelarrayparam.a_in,
        poisson_j_perp_in=eyoungparallelarrayparam.poisson_j_perp_in,
    )

    assert eyoung_j_out == pytest.approx(eyoungparallelarrayparam.expected_eyoung_j_out)

    assert a_out == pytest.approx(eyoungparallelarrayparam.expected_a_out)

    assert poisson_j_perp_out == pytest.approx(
        eyoungparallelarrayparam.expected_poisson_j_perp_out
    )


@pytest.mark.parametrize(
    "sx, sy, sz, expected",
    (
        (0, -3.2e8, 2.4e8, 486621002.42385757),
        (-2.8e8, 0, 2.4e8, 450777106.7833858),
    ),
)
def test_sigvm(sx, sy, sz, expected):
    # could not find an example of a use in PROCESS where
    # tx, ty, or tz were anything other than 0
    ret = tfcoil_module.sigvm(sx, sy, sz, 0, 0, 0)

    assert ret == pytest.approx(expected)


@pytest.mark.parametrize(
    "ind_tf_coil, c_tf_total, n_tf_coils, expected_total, expected_total_gj, expected_single",
    [
        (1.0, 2.0, 3, 2.0, 2.0e-9, 2 / 3),
        (0.5, 4.0, 2, 4.0, 4.0e-9, 2.0),
        (2.0, 5.0, 4, 25.0, 25.0e-9, 6.25),
        (0.0, 5.0, 1, 0.0, 0.0, 0.0),
        (1.0, 0.0, 10, 0.0, 0.0, 0.0),
    ],
)
def test_tf_stored_magnetic_energy(
    ind_tf_coil,
    c_tf_total,
    n_tf_coils,
    expected_total,
    expected_total_gj,
    expected_single,
):
    tfc = TFCoil(build=None)
    result = tfc.tf_stored_magnetic_energy(
        ind_tf_coil=ind_tf_coil, c_tf_total=c_tf_total, n_tf_coils=n_tf_coils
    )
    assert pytest.approx(result[0]) == expected_total
    assert pytest.approx(result[1]) == expected_total_gj
    assert pytest.approx(result[2]) == expected_single
