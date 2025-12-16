from typing import Any, NamedTuple

import numpy as np
import pytest

from process.data_structure import (
    build_variables,
    buildings_variables,
    constraint_variables,
    cost_variables,
    current_drive_variables,
    fwbs_variables,
    heat_transport_variables,
    numerics,
    pf_power_variables,
    pfcoil_variables,
    physics_variables,
    power_variables,
    tfcoil_variables,
    times_variables,
)
from process.data_structure import primary_pumping_variables as ppv
from process.power import Power


@pytest.fixture
def power():
    """Provides power object for testing.

    :returns: initialised power object
    :rtype: process.power.Power
    """
    return Power()


class CryoParam(NamedTuple):
    qnuc: Any = None

    inuclear: Any = None

    qss: Any = None

    qac: Any = None

    qcl: Any = None

    qmisc: Any = None

    i_tf_sup: Any = None

    coldmass: Any = None

    c_tf_turn: Any = None

    ensxpfm: Any = None

    p_tf_nuclear_heat_mw: Any = None

    n_tf_coils: Any = None

    tfcryoarea: Any = None

    t_plant_pulse_plasma_present: Any = None

    expected_qss: Any = None

    expected_qac: Any = None

    expected_qcl: Any = None

    expected_qmisc: Any = None

    expected_helpow: Any = None


@pytest.mark.parametrize(
    "cryoparam",
    (
        CryoParam(
            qnuc=12920,
            inuclear=1,
            qss=0,
            qac=0,
            qcl=0,
            qmisc=0,
            i_tf_sup=1,
            coldmass=47352637.039762333,
            c_tf_turn=74026.751437500003,
            ensxpfm=37429.525515086898,
            p_tf_nuclear_heat_mw=0.044178296011112193,
            n_tf_coils=16,
            tfcryoarea=0,
            t_plant_pulse_plasma_present=10364.426139387357,
            expected_qss=20361.633927097802,
            expected_qac=3611.3456752656607,
            expected_qcl=16108.2211128,
            expected_qmisc=23850.540321823562,
            expected_helpow=76851.741036987034,
        ),
        CryoParam(
            qnuc=12920,
            inuclear=1,
            qss=20361.633927097802,
            qac=3611.3456752656607,
            qcl=16108.2211128,
            qmisc=23850.540321823562,
            i_tf_sup=1,
            coldmass=47308985.527808741,
            c_tf_turn=74026.751437500003,
            ensxpfm=37427.228965055205,
            p_tf_nuclear_heat_mw=0.045535131445547841,
            n_tf_coils=16,
            tfcryoarea=0,
            t_plant_pulse_plasma_present=364.42613938735633,
            expected_qss=20342.863776957758,
            expected_qac=102701.82327748176,
            expected_qcl=16108.2211128,
            expected_qmisc=68432.80867525778,
            expected_helpow=220505.71684249729,
        ),
    ),
)
def test_cryo(cryoparam, monkeypatch, power):
    """
    Automatically generated Regression Unit Test for cryo.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param cryoparam: the data used to mock and assert in this test.
    :type cryoparam: cryoparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(fwbs_variables, "qnuc", cryoparam.qnuc)

    monkeypatch.setattr(fwbs_variables, "inuclear", cryoparam.inuclear)

    monkeypatch.setattr(power_variables, "qss", cryoparam.qss)

    monkeypatch.setattr(power_variables, "qac", cryoparam.qac)

    monkeypatch.setattr(power_variables, "qcl", cryoparam.qcl)

    monkeypatch.setattr(power_variables, "qmisc", cryoparam.qmisc)

    helpow = power.cryo(
        i_tf_sup=cryoparam.i_tf_sup,
        coldmass=cryoparam.coldmass,
        c_tf_turn=cryoparam.c_tf_turn,
        ensxpfm=cryoparam.ensxpfm,
        p_tf_nuclear_heat_mw=cryoparam.p_tf_nuclear_heat_mw,
        n_tf_coils=cryoparam.n_tf_coils,
        tfcryoarea=cryoparam.tfcryoarea,
        t_plant_pulse_plasma_present=cryoparam.t_plant_pulse_plasma_present,
    )

    assert power_variables.qss == pytest.approx(cryoparam.expected_qss)

    assert power_variables.qac == pytest.approx(cryoparam.expected_qac)

    assert power_variables.qcl == pytest.approx(cryoparam.expected_qcl)

    assert power_variables.qmisc == pytest.approx(cryoparam.expected_qmisc)

    assert helpow == pytest.approx(cryoparam.expected_helpow)


class PfpwrParam(NamedTuple):
    iohcl: Any = None

    peakmva: Any = None

    pfckts: Any = None

    maxpoloidalpower: Any = None

    peakpoloidalpower: Any = None

    spfbusl: Any = None

    poloidalpower: Any = None

    spsmva: Any = None

    vpfskv: Any = None

    ensxpfm: Any = None

    acptmax: Any = None

    srcktpm: Any = None

    n_pf_coil_groups: Any = None

    c_pf_coil_turn: Any = None

    p_pf_electric_supplies_mw: Any = None

    rho_pf_coil: Any = None

    n_pf_cs_plasma_circuits: Any = None

    n_pf_coils_in_group: Any = None

    c_pf_cs_coils_peak_ma: Any = None

    etapsu: Any = None

    c_pf_coil_turn_peak_input: Any = None

    c_pf_cs_coil_pulse_end_ma: Any = None

    ind_pf_cs_plasma_mutual: Any = None

    n_pf_coil_turns: Any = None

    f_a_pf_coil_void: Any = None

    j_pf_coil_wp_peak: Any = None

    r_pf_coil_middle: Any = None

    p_plasma_ohmic_mw: Any = None

    rmajor: Any = None

    active_constraints: Any = None

    ioptimz: Any = None

    t_pulse_cumulative: Any = None

    intervallabel: Any = None

    timelabel: Any = None

    t_plant_pulse_plasma_current_ramp_up: Any = None

    outfile: Any = None

    iprint: Any = None

    expected_peakmva: Any = None

    expected_pfckts: Any = None

    expected_peakpoloidalpower: Any = None

    expected_spfbusl: Any = None

    expected_poloidalpower: Any = None

    expected_spsmva: Any = None

    expected_vpfskv: Any = None

    expected_ensxpfm: Any = None

    expected_acptmax: Any = None

    expected_srcktpm: Any = None


@pytest.mark.parametrize(
    "pfpwrparam",
    (
        PfpwrParam(
            iohcl=1,
            peakmva=0,
            pfckts=0,
            maxpoloidalpower=1000,
            peakpoloidalpower=0,
            spfbusl=0,
            poloidalpower=np.array(
                np.array((0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            spsmva=0,
            vpfskv=0,
            ensxpfm=0,
            acptmax=0,
            srcktpm=0,
            n_pf_coil_groups=4,
            c_pf_coil_turn=np.array(
                (
                    (
                        0,
                        0,
                        -0,
                        -0,
                        -0,
                        -0,
                        -0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        42200,
                        42200,
                        3020.1587941721036,
                        3020.1587941721036,
                        3300.7614790391262,
                        3300.7614790391262,
                        40065.680000000008,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        192.99998911734409,
                        -6144.2544496188857,
                        -42200,
                        -42200,
                        -43000,
                        -43000,
                        -43000,
                        17721306.969367817,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        192.99998911734409,
                        -6144.2544496188857,
                        -42200,
                        -42200,
                        -43000,
                        -43000,
                        -43000,
                        17721306.969367817,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        192.99998911734409,
                        -6144.2544496188857,
                        -42200,
                        -42200,
                        -43000,
                        -43000,
                        -43000,
                        17721306.969367817,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0,
                        0,
                        -0,
                        -0,
                        -0,
                        -0,
                        -0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                ),
                order="F",
            ).transpose(),
            p_pf_electric_supplies_mw=0,
            rho_pf_coil=0,
            n_pf_cs_plasma_circuits=8,
            n_pf_coils_in_group=np.array(
                np.array((1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            c_pf_cs_coils_peak_ma=np.array(
                np.array(
                    (
                        14.742063826112622,
                        20.032681634901664,
                        -8.1098913365453491,
                        -8.1098913365453491,
                        -5.5984385047179153,
                        -5.5984385047179153,
                        -186.98751599968145,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            etapsu=0.90000000000000002,
            c_pf_coil_turn_peak_input=np.array(
                np.array(
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
                    order="F",
                ),
                order="F",
            ).transpose(),
            c_pf_cs_coil_pulse_end_ma=np.array(
                np.array(
                    (
                        0.067422231232391661,
                        -2.9167273287450968,
                        -8.1098913365453491,
                        -8.1098913365453491,
                        -5.5984385047179153,
                        -5.5984385047179153,
                        -186.98751599968148,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            ind_pf_cs_plasma_mutual=np.array(
                (
                    (
                        2.4933245328128875,
                        0.044628616646610005,
                        0.23809409972275392,
                        0.15765363220324294,
                        0.21869592803714374,
                        0.066200200497513878,
                        0.88106839153571348,
                        0.0008151322258474719,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.044628616646610005,
                        4.3316688790171352,
                        0.18920709024491231,
                        0.2933332275969987,
                        0.078421246973137196,
                        0.283752898388758,
                        0.85440319548278287,
                        0.00086087843592316565,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.23809409972275389,
                        0.18920709024491231,
                        3.128721334334037,
                        1.1084361059087036,
                        0.72476925375751233,
                        0.39082336057406458,
                        0.54626354354859585,
                        0.0017044090640384037,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.15765363220324294,
                        0.2933332275969987,
                        1.1084361059087036,
                        3.128721334334037,
                        0.39082336057406458,
                        0.72476925375751233,
                        0.54626354354859585,
                        0.0017044090640384037,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.21869592803714374,
                        0.078421246973137196,
                        0.72476925375751244,
                        0.39082336057406464,
                        1.3966126540799821,
                        0.15016488330980787,
                        0.32769603485124171,
                        0.00088156051922038358,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.066200200497513864,
                        0.28375289838875795,
                        0.39082336057406464,
                        0.72476925375751244,
                        0.15016488330980787,
                        1.3966126540799821,
                        0.32769603485124171,
                        0.00088156051922038358,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.88106839153571348,
                        0.85440319548278287,
                        0.54626354354859585,
                        0.54626354354859585,
                        0.32769603485124171,
                        0.32769603485124171,
                        25.013930780082362,
                        0.0049030712741391239,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.0008151322258474719,
                        0.00086087843592316565,
                        0.0017044090640384037,
                        0.0017044090640384037,
                        0.00088156051922038358,
                        0.00088156051922038358,
                        0.0049030712741391239,
                        1.6039223939491056e-05,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                ),
                order="F",
            ).transpose(),
            n_pf_coil_turns=np.array(
                np.array(
                    (
                        349.33800535811901,
                        474.70809561378354,
                        192.17751982334951,
                        192.17751982334951,
                        130.19624429576547,
                        130.19624429576547,
                        4348.5468837135222,
                        1,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            f_a_pf_coil_void=np.array(
                np.array(
                    (
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            j_pf_coil_wp_peak=np.array(
                np.array(
                    (
                        11000000,
                        11000000,
                        6000000,
                        6000000,
                        8000000,
                        8000000,
                        8000000,
                        8000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            r_pf_coil_middle=np.array(
                np.array(
                    (
                        6.2732560483870969,
                        6.2732560483870969,
                        18.401280308184159,
                        18.401280308184159,
                        16.803394770584916,
                        16.803394770584916,
                        2.6084100000000001,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            p_plasma_ohmic_mw=0.61391840981850698,
            rmajor=8.8901000000000003,
            active_constraints=(
                True,
                True,
                False,
                False,
                True,
                False,
                False,
                True,
                False,
                False,
                True,
                False,
                True,
                False,
                True,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                True,
                True,
                True,
                True,
                False,
                False,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                True,
                False,
                True,
                False,
                False,
                True,
                False,
                False,
                True,
                False,
                False,
                False,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ),
            ioptimz=1,
            t_pulse_cumulative=np.array(
                np.array(
                    (
                        0,
                        500,
                        677.21306969367811,
                        687.21306969367811,
                        10687.213069693678,
                        10864.426139387357,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            intervallabel=(
                "t_plant_pulse_coil_precharge      ",
                "t_plant_pulse_plasma_current_ramp_up       ",
                "t_plant_pulse_fusion_ramp      ",
                "t_plant_pulse_burn      ",
                "t_plant_pulse_plasma_current_ramp_down      ",
            ),
            timelabel=(
                "Start      ",
                "BOP        ",
                "EOR        ",
                "BOF        ",
                "EOF        ",
                "EOP        ",
            ),
            t_plant_pulse_plasma_current_ramp_up=177.21306969367816,
            outfile=11,
            iprint=0,
            expected_peakmva=736.39062584245937,
            expected_pfckts=12,
            expected_peakpoloidalpower=211.21199231967319,
            expected_spfbusl=2533.4495999999999,
            expected_poloidalpower=np.array(
                np.array(
                    (59332953.082890816, 43806300.444207191, 0, 0, -211211992.31967318),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_spsmva=845.66824574150155,
            expected_vpfskv=20,
            expected_ensxpfm=37429.525515086898,
            expected_acptmax=24.816666666666666,
            expected_srcktpm=1071.1112934857531,
        ),
        PfpwrParam(
            iohcl=1,
            peakmva=736.39062584245937,
            pfckts=12,
            maxpoloidalpower=1000,
            peakpoloidalpower=211.21199231967319,
            spfbusl=2533.4495999999999,
            poloidalpower=np.array(
                np.array(
                    (59332953.082890816, 43806300.444207191, 0, 0, -211211992.31967318),
                    order="F",
                ),
                order="F",
            ).transpose(),
            spsmva=845.66824574150155,
            vpfskv=20,
            ensxpfm=37429.525515086898,
            acptmax=24.816666666666666,
            srcktpm=1071.1112934857531,
            n_pf_coil_groups=4,
            c_pf_coil_turn=np.array(
                (
                    (
                        0,
                        0,
                        -0,
                        -0,
                        -0,
                        -0,
                        -0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        33663.946773824558,
                        38185.429487079651,
                        3066.1011211106556,
                        3066.1011211106556,
                        3142.8828598960072,
                        3142.8828598960072,
                        40065.680000000008,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        42200,
                        42200,
                        -38360.428378196812,
                        -38360.428378196812,
                        -39064.277281521267,
                        -39064.277281521267,
                        7172.8553168274502,
                        17721306.969367817,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        42200,
                        42200,
                        -38360.428378196812,
                        -38360.428378196812,
                        -39064.277281521267,
                        -39064.277281521267,
                        7172.8553168274502,
                        17721306.969367817,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        43.81218847453755,
                        -5618.2831008025678,
                        -42200,
                        -42200,
                        -43000,
                        -43000,
                        -43000,
                        17721306.969367817,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0,
                        0,
                        -0,
                        -0,
                        -0,
                        -0,
                        -0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                ),
                order="F",
            ).transpose(),
            p_pf_electric_supplies_mw=0.89998039031509891,
            rho_pf_coil=0,
            n_pf_cs_plasma_circuits=8,
            n_pf_coils_in_group=np.array(
                np.array((1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            c_pf_cs_coils_peak_ma=np.array(
                np.array(
                    (
                        18.579095475129442,
                        22.175439215004367,
                        -8.1210132461605742,
                        -8.1210132461605742,
                        -5.575080047168135,
                        -5.575080047168135,
                        -186.98751599968145,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            etapsu=0.90000000000000002,
            c_pf_coil_turn_peak_input=np.array(
                np.array(
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
                    order="F",
                ),
                order="F",
            ).transpose(),
            c_pf_cs_coil_pulse_end_ma=np.array(
                np.array(
                    (
                        0.019288882290113718,
                        -2.9523197960789949,
                        -8.1210132461605742,
                        -8.1210132461605742,
                        -5.575080047168135,
                        -5.575080047168135,
                        -186.98751599968148,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            ind_pf_cs_plasma_mutual=np.array(
                (
                    (
                        3.7834082671748859,
                        0.062121647727783093,
                        0.30015331189839162,
                        0.19867383883991577,
                        0.27436487704364948,
                        0.082948031292997063,
                        1.1061712527993555,
                        0.0010241850221498481,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.0621216477277831,
                        5.1972808039781917,
                        0.20973249911052369,
                        0.32515436295986727,
                        0.086447229668541736,
                        0.3127934446710578,
                        0.94579280357174933,
                        0.00095296065575475799,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.30015331189839167,
                        0.20973249911052369,
                        3.136721879042204,
                        1.1114784104069368,
                        0.7227350889143983,
                        0.38972646092520974,
                        0.54701268968450789,
                        0.0017067464916032077,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.19867383883991577,
                        0.32515436295986722,
                        1.1114784104069368,
                        3.136721879042204,
                        0.38972646092520974,
                        0.7227350889143983,
                        0.54701268968450789,
                        0.0017067464916032077,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.27436487704364942,
                        0.08644722966854175,
                        0.72273508891439842,
                        0.38972646092520968,
                        1.385724786008854,
                        0.14891442656236412,
                        0.32632878326620535,
                        0.00087788236968843478,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.082948031292997063,
                        0.31279344467105774,
                        0.38972646092520968,
                        0.72273508891439842,
                        0.14891442656236412,
                        1.385724786008854,
                        0.32632878326620535,
                        0.00087788236968843478,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        1.1061712527993555,
                        0.94579280357174933,
                        0.54701268968450789,
                        0.54701268968450789,
                        0.32632878326620535,
                        0.32632878326620535,
                        25.013930780082362,
                        0.0049030712741391239,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (
                        0.0010241850221498481,
                        0.00095296065575475799,
                        0.0017067464916032077,
                        0.0017067464916032077,
                        0.00087788236968843478,
                        0.00087788236968843478,
                        0.0049030712741391239,
                        1.6039223939491056e-05,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
                ),
                order="F",
            ).transpose(),
            n_pf_coil_turns=np.array(
                np.array(
                    (
                        440.26292595093463,
                        525.48434158778116,
                        192.44107218389988,
                        192.44107218389988,
                        129.65302435274731,
                        129.65302435274731,
                        4348.5468837135222,
                        1,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                        100,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            f_a_pf_coil_void=np.array(
                np.array(
                    (
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                        0.29999999999999999,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            j_pf_coil_wp_peak=np.array(
                np.array(
                    (
                        11000000,
                        11000000,
                        6000000,
                        6000000,
                        8000000,
                        8000000,
                        8000000,
                        8000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                        30000000,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            r_pf_coil_middle=np.array(
                np.array(
                    (
                        6.2732560483870969,
                        6.2732560483870969,
                        18.401280308184159,
                        18.401280308184159,
                        16.803394770584916,
                        16.803394770584916,
                        2.6084100000000001,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            p_plasma_ohmic_mw=0.61391840981850698,
            rmajor=8.8901000000000003,
            active_constraints=(
                True,
                True,
                False,
                False,
                True,
                False,
                False,
                True,
                False,
                False,
                True,
                False,
                True,
                False,
                True,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                True,
                True,
                True,
                True,
                False,
                False,
                True,
                True,
                True,
                True,
                True,
                True,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                True,
                False,
                True,
                False,
                False,
                True,
                False,
                False,
                True,
                False,
                False,
                False,
                True,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
                False,
            ),
            ioptimz=1,
            t_pulse_cumulative=np.array(
                np.array(
                    (
                        0,
                        500,
                        677.21306969367811,
                        687.21306969367811,
                        687.21306969367811,
                        864.42613938735622,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            intervallabel=(
                "t_plant_pulse_coil_precharge      ",
                "t_plant_pulse_plasma_current_ramp_up       ",
                "t_plant_pulse_fusion_ramp      ",
                "t_plant_pulse_burn      ",
                "t_plant_pulse_plasma_current_ramp_down      ",
            ),
            timelabel=(
                "Start      ",
                "BOP        ",
                "EOR        ",
                "BOF        ",
                "EOF        ",
                "EOP        ",
            ),
            t_plant_pulse_plasma_current_ramp_up=177.21306969367816,
            outfile=11,
            iprint=0,
            expected_peakmva=90.673341440806112,
            expected_pfckts=12,
            expected_peakpoloidalpower=9900,
            expected_spfbusl=2533.4495999999999,
            expected_poloidalpower=np.array(
                np.array(
                    (
                        59043243.553314812,
                        -69656470.894853994,
                        0,
                        9900000000,
                        -211199033.0608803,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_spsmva=354.86210489492782,
            expected_vpfskv=20,
            expected_ensxpfm=37427.228965055205,
            expected_acptmax=24.816666666666666,
            expected_srcktpm=1069.8879533693198,
        ),
    ),
)
def test_pfpwr(pfpwrparam, monkeypatch, power):
    """
    Automatically generated Regression Unit Test for pfpwr.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param pfpwrparam: the data used to mock and assert in this test.
    :type pfpwrparam: pfpwrparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "iohcl", pfpwrparam.iohcl)

    monkeypatch.setattr(heat_transport_variables, "peakmva", pfpwrparam.peakmva)

    monkeypatch.setattr(pf_power_variables, "pfckts", pfpwrparam.pfckts)

    monkeypatch.setattr(
        pf_power_variables, "maxpoloidalpower", pfpwrparam.maxpoloidalpower
    )

    monkeypatch.setattr(
        pf_power_variables, "peakpoloidalpower", pfpwrparam.peakpoloidalpower
    )

    monkeypatch.setattr(pf_power_variables, "spfbusl", pfpwrparam.spfbusl)

    monkeypatch.setattr(pf_power_variables, "poloidalpower", pfpwrparam.poloidalpower)

    monkeypatch.setattr(pf_power_variables, "spsmva", pfpwrparam.spsmva)

    monkeypatch.setattr(pf_power_variables, "vpfskv", pfpwrparam.vpfskv)

    monkeypatch.setattr(pf_power_variables, "ensxpfm", pfpwrparam.ensxpfm)

    monkeypatch.setattr(pf_power_variables, "acptmax", pfpwrparam.acptmax)

    monkeypatch.setattr(pf_power_variables, "srcktpm", pfpwrparam.srcktpm)

    monkeypatch.setattr(
        pfcoil_variables, "n_pf_coil_groups", pfpwrparam.n_pf_coil_groups
    )

    monkeypatch.setattr(pfcoil_variables, "c_pf_coil_turn", pfpwrparam.c_pf_coil_turn)

    monkeypatch.setattr(
        pfcoil_variables,
        "p_pf_electric_supplies_mw",
        pfpwrparam.p_pf_electric_supplies_mw,
    )

    monkeypatch.setattr(pfcoil_variables, "rho_pf_coil", pfpwrparam.rho_pf_coil)

    monkeypatch.setattr(
        pfcoil_variables, "n_pf_cs_plasma_circuits", pfpwrparam.n_pf_cs_plasma_circuits
    )

    monkeypatch.setattr(
        pfcoil_variables, "n_pf_coils_in_group", pfpwrparam.n_pf_coils_in_group
    )

    monkeypatch.setattr(
        pfcoil_variables, "c_pf_cs_coils_peak_ma", pfpwrparam.c_pf_cs_coils_peak_ma
    )

    monkeypatch.setattr(pfcoil_variables, "etapsu", pfpwrparam.etapsu)

    monkeypatch.setattr(
        pfcoil_variables,
        "c_pf_coil_turn_peak_input",
        pfpwrparam.c_pf_coil_turn_peak_input,
    )

    monkeypatch.setattr(
        pfcoil_variables,
        "c_pf_cs_coil_pulse_end_ma",
        pfpwrparam.c_pf_cs_coil_pulse_end_ma,
    )

    monkeypatch.setattr(
        pfcoil_variables, "ind_pf_cs_plasma_mutual", pfpwrparam.ind_pf_cs_plasma_mutual
    )

    monkeypatch.setattr(pfcoil_variables, "n_pf_coil_turns", pfpwrparam.n_pf_coil_turns)

    monkeypatch.setattr(
        pfcoil_variables, "f_a_pf_coil_void", pfpwrparam.f_a_pf_coil_void
    )

    monkeypatch.setattr(
        pfcoil_variables, "j_pf_coil_wp_peak", pfpwrparam.j_pf_coil_wp_peak
    )

    monkeypatch.setattr(
        pfcoil_variables, "r_pf_coil_middle", pfpwrparam.r_pf_coil_middle
    )

    monkeypatch.setattr(
        physics_variables, "p_plasma_ohmic_mw", pfpwrparam.p_plasma_ohmic_mw
    )

    monkeypatch.setattr(physics_variables, "rmajor", pfpwrparam.rmajor)

    monkeypatch.setattr(numerics, "active_constraints", pfpwrparam.active_constraints)

    monkeypatch.setattr(numerics, "ioptimz", pfpwrparam.ioptimz)

    monkeypatch.setattr(
        times_variables, "t_pulse_cumulative", pfpwrparam.t_pulse_cumulative
    )

    monkeypatch.setattr(
        times_variables,
        "t_plant_pulse_plasma_current_ramp_up",
        pfpwrparam.t_plant_pulse_plasma_current_ramp_up,
    )

    power.pfpwr(output=False)

    assert heat_transport_variables.peakmva == pytest.approx(pfpwrparam.expected_peakmva)

    assert pf_power_variables.pfckts == pytest.approx(pfpwrparam.expected_pfckts)

    assert pf_power_variables.peakpoloidalpower == pytest.approx(
        pfpwrparam.expected_peakpoloidalpower
    )

    assert pf_power_variables.spfbusl == pytest.approx(pfpwrparam.expected_spfbusl)

    assert pf_power_variables.poloidalpower == pytest.approx(
        pfpwrparam.expected_poloidalpower
    )

    assert pf_power_variables.spsmva == pytest.approx(pfpwrparam.expected_spsmva)

    assert pf_power_variables.vpfskv == pytest.approx(pfpwrparam.expected_vpfskv)

    assert pf_power_variables.ensxpfm == pytest.approx(pfpwrparam.expected_ensxpfm)

    assert pf_power_variables.acptmax == pytest.approx(pfpwrparam.expected_acptmax)

    assert pf_power_variables.srcktpm == pytest.approx(pfpwrparam.expected_srcktpm)


class AcpowParam(NamedTuple):
    a_plant_floor_effective: Any = None

    p_plant_electric_base: Any = None

    p_cryo_plant_electric_mw: Any = None

    vachtmw: Any = None

    p_tf_electric_supplies_mw: Any = None

    p_tritium_plant_electric_mw: Any = None

    p_hcd_electric_total_mw: Any = None

    tlvpmw: Any = None

    peakmva: Any = None

    p_plant_electric_base_total_mw: Any = None

    fmgdmw: Any = None

    pflux_plant_floor_electric: Any = None

    p_coolant_pump_elec_total_mw: Any = None

    pacpmw: Any = None

    i_pf_energy_storage_source: Any = None

    srcktpm: Any = None

    iprint: Any = None

    outfile: Any = None

    expected_pacpmw: Any = None


@pytest.mark.parametrize(
    "acpowparam",
    (
        AcpowParam(
            a_plant_floor_effective=379218.8908858358,
            p_plant_electric_base=5000000,
            p_cryo_plant_electric_mw=37.900388528497025,
            vachtmw=0.5,
            p_tf_electric_supplies_mw=9.1507079104675704,
            p_tritium_plant_electric_mw=15,
            p_hcd_electric_total_mw=129.94611930107126,
            tlvpmw=0,
            peakmva=736.39062584245937,
            p_plant_electric_base_total_mw=0,
            fmgdmw=0,
            pflux_plant_floor_electric=150,
            p_coolant_pump_elec_total_mw=234.28554165620102,
            pacpmw=0,
            i_pf_energy_storage_source=2,
            srcktpm=1071.1112934857531,
            iprint=0,
            outfile=11,
            expected_pacpmw=1164.244494532182,
        ),
        AcpowParam(
            a_plant_floor_effective=381580.9594357388,
            p_plant_electric_base=5000000,
            p_cryo_plant_electric_mw=108.74512702403499,
            vachtmw=0.5,
            p_tf_electric_supplies_mw=9.1507079104675704,
            p_tritium_plant_electric_mw=15,
            p_hcd_electric_total_mw=129.94611930107126,
            tlvpmw=699.34943812129745,
            peakmva=90.673341440806112,
            p_plant_electric_base_total_mw=62.23714391536082,
            fmgdmw=0,
            pflux_plant_floor_electric=150,
            p_coolant_pump_elec_total_mw=234.2162627659944,
            pacpmw=1226.1273281650574,
            i_pf_energy_storage_source=2,
            srcktpm=1069.8879533693198,
            iprint=0,
            outfile=11,
            expected_pacpmw=589.3014463957436,
        ),
    ),
)
def test_acpow(acpowparam, monkeypatch, power):
    """
    Automatically generated Regression Unit Test for acpow.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acpowparam: the data used to mock and assert in this test.
    :type acpowparam: acpowparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        buildings_variables,
        "a_plant_floor_effective",
        acpowparam.a_plant_floor_effective,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_base",
        acpowparam.p_plant_electric_base,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_cryo_plant_electric_mw",
        acpowparam.p_cryo_plant_electric_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "vachtmw", acpowparam.vachtmw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_tf_electric_supplies_mw",
        acpowparam.p_tf_electric_supplies_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_tritium_plant_electric_mw",
        acpowparam.p_tritium_plant_electric_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_total_mw",
        acpowparam.p_hcd_electric_total_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "tlvpmw", acpowparam.tlvpmw)

    monkeypatch.setattr(heat_transport_variables, "peakmva", acpowparam.peakmva)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_base_total_mw",
        acpowparam.p_plant_electric_base_total_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "fmgdmw", acpowparam.fmgdmw)

    monkeypatch.setattr(
        heat_transport_variables,
        "pflux_plant_floor_electric",
        acpowparam.pflux_plant_floor_electric,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_coolant_pump_elec_total_mw",
        acpowparam.p_coolant_pump_elec_total_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "pacpmw", acpowparam.pacpmw)

    monkeypatch.setattr(
        pf_power_variables,
        "i_pf_energy_storage_source",
        acpowparam.i_pf_energy_storage_source,
    )

    monkeypatch.setattr(pf_power_variables, "srcktpm", acpowparam.srcktpm)

    power.acpow(output=False)

    assert heat_transport_variables.pacpmw == pytest.approx(acpowparam.expected_pacpmw)


class Power2Param(NamedTuple):
    p_plant_electric_net_required_mw: Any = None

    ipnet: Any = None

    ireactor: Any = None

    p_hcd_injected_total_mw: Any = None

    p_blkt_multiplication_mw: Any = None

    inuclear: Any = None

    p_blkt_nuclear_heat_total_mw: Any = None

    p_fw_rad_total_mw: Any = None

    qnuc: Any = None

    eta_coolant_pump_electric: Any = None

    f_p_blkt_multiplication: Any = None

    p_div_rad_total_mw: Any = None

    f_ster_div_single: Any = None

    f_a_fw_outboard_hcd: Any = None

    i_thermal_electric_conversion: Any = None

    pnuc_cp: Any = None

    p_div_nuclear_heat_total_mw: Any = None

    i_p_coolant_pumping: Any = None

    p_tf_nuclear_heat_mw: Any = None

    p_fw_hcd_nuclear_heat_mw: Any = None

    p_shld_nuclear_heat_mw: Any = None

    p_fw_hcd_rad_total_mw: Any = None

    p_fw_nuclear_heat_total_mw: Any = None

    p_shld_coolant_pump_mw: Any = None

    p_blkt_coolant_pump_mw: Any = None

    p_shld_secondary_heat_mw: Any = None

    f_p_shld_coolant_pump_total_heat: Any = None

    temp_turbine_coolant_in: Any = None

    p_plant_electric_net_mw: Any = None

    f_p_div_coolant_pump_total_heat: Any = None

    f_p_blkt_coolant_pump_total_heat: Any = None

    vachtmw: Any = None

    p_div_coolant_pump_mw: Any = None

    n_primary_heat_exchangers: Any = None

    helpow: Any = None

    p_fw_coolant_pump_mw: Any = None

    p_plant_electric_recirc_mw: Any = None

    p_plant_primary_heat_mw: Any = None

    f_p_fw_coolant_pump_total_heat: Any = None

    p_plant_electric_base_total_mw: Any = None

    i_shld_primary_heat: Any = None

    p_hcd_electric_total_mw: Any = None

    fachtmw: Any = None

    p_plant_electric_gross_mw: Any = None

    p_plant_secondary_heat_mw: Any = None

    p_tritium_plant_electric_mw: Any = None

    p_hcd_secondary_heat_mw: Any = None

    p_tf_electric_supplies_mw: Any = None

    p_coolant_pump_elec_total_mw: Any = None

    eta_turbine: Any = None

    p_cryo_plant_electric_mw: Any = None

    p_div_secondary_heat_mw: Any = None

    p_hcd_electric_loss_mw: Any = None

    p_coolant_pump_loss_total_mw: Any = None

    helpow_cryal: Any = None

    p_pf_electric_supplies_mw: Any = None

    p_alpha_total_mw: Any = None

    i_plasma_ignited: Any = None

    p_plasma_inner_rad_mw: Any = None

    p_plasma_rad_mw: Any = None

    itart: Any = None

    p_plasma_separatrix_mw: Any = None

    p_fw_alpha_mw: Any = None

    n_divertors: Any = None

    p_plasma_ohmic_mw: Any = None

    i_rad_loss: Any = None

    p_fusion_total_mw: Any = None

    p_non_alpha_charged_mw: Any = None

    pscalingmw: Any = None

    f_p_alpha_plasma_deposited: Any = None

    p_cp_coolant_pump_elec: Any = None

    i_tf_sup: Any = None

    tfcmw: Any = None

    temp_tf_cryo: Any = None

    tcoolin: Any = None

    eff_tf_cryo: Any = None

    p_fw_blkt_coolant_pump_mw: Any = None

    p_shld_coolant_pump_elec_mw: Any = None

    p_div_coolant_pump_elec_mw: Any = None

    p_coolant_pump_total_mw: Any = None

    p_fw_blkt_heat_deposited_mw: Any = None

    p_fw_blkt_coolant_pump_elec_mw: Any = None

    p_div_heat_deposited_mw: Any = None

    p_fw_heat_deposited_mw: Any = None

    p_shld_heat_deposited_mw: Any = None

    p_cp_coolant_pump_elec_mw: Any = None

    p_plant_core_systems_elec_mw: Any = None

    f_p_div_primary_heat: Any = None

    qss: Any = None

    qac: Any = None

    qcl: Any = None

    qmisc: Any = None

    outfile: Any = None

    iprint: Any = None

    expected_p_plant_electric_net_mw: Any = None

    expected_precircmw: Any = None

    expected_fachtmw: Any = None

    expected_p_plant_electric_gross_mw: Any = None

    expected_p_plant_secondary_heat_mw: Any = None

    expected_pcoresystems: Any = None


@pytest.mark.parametrize(
    "power2param",
    (
        Power2Param(
            p_plant_electric_net_required_mw=500,
            ipnet=0,
            ireactor=1,
            p_hcd_injected_total_mw=51.978447720428512,
            p_blkt_multiplication_mw=377.93233088402548,
            inuclear=1,
            p_blkt_nuclear_heat_total_mw=1504.711566619962,
            p_fw_rad_total_mw=254.87601794907812,
            qnuc=12920,
            eta_coolant_pump_electric=0.87000000000000011,
            f_p_blkt_multiplication=1.2690000534057617,
            p_div_rad_total_mw=33.119482558354782,
            f_ster_div_single=0.115,
            f_a_fw_outboard_hcd=0,
            i_thermal_electric_conversion=2,
            pnuc_cp=0,
            p_div_nuclear_heat_total_mw=182.69222981118057,
            i_p_coolant_pumping=3,
            p_tf_nuclear_heat_mw=0.044178296011112193,
            p_fw_hcd_nuclear_heat_mw=0,
            p_shld_nuclear_heat_mw=1.3609360176065353,
            p_fw_hcd_rad_total_mw=0,
            p_fw_nuclear_heat_total_mw=276.76827393356979,
            p_shld_coolant_pump_mw=0.0068046800880326762,
            p_blkt_coolant_pump_mw=0,
            p_shld_secondary_heat_mw=0,
            f_p_shld_coolant_pump_total_heat=0.0050000000000000001,
            temp_turbine_coolant_in=0,
            p_plant_electric_net_mw=0,
            f_p_div_coolant_pump_total_heat=0.0050000000000000001,
            f_p_blkt_coolant_pump_total_heat=0.0050000000000000001,
            vachtmw=0.5,
            p_div_coolant_pump_mw=1.7942175899286208,
            n_primary_heat_exchangers=3,
            helpow=76851.741036987034,
            p_fw_coolant_pump_mw=0,
            p_plant_electric_recirc_mw=0,
            p_plant_primary_heat_mw=2620.2218111502593,
            f_p_fw_coolant_pump_total_heat=0.0050000000000000001,
            p_plant_electric_base_total_mw=62.23714391536082,
            i_shld_primary_heat=1,
            p_hcd_electric_total_mw=129.94611930107126,
            fachtmw=0,
            p_plant_electric_gross_mw=0,
            p_plant_secondary_heat_mw=0,
            p_tritium_plant_electric_mw=15,
            p_hcd_secondary_heat_mw=0,
            p_tf_electric_supplies_mw=9.1507079104675704,
            p_coolant_pump_elec_total_mw=234.28554165620102,
            eta_turbine=0.37500000000000006,
            p_cryo_plant_electric_mw=37.900388528497025,
            p_div_secondary_heat_mw=0,
            p_hcd_electric_loss_mw=77.967671580642758,
            p_coolant_pump_loss_total_mw=30.457120415306122,
            helpow_cryal=0,
            p_pf_electric_supplies_mw=0.89998039031509891,
            p_alpha_total_mw=396.66154806848488,
            i_plasma_ignited=0,
            p_plasma_inner_rad_mw=113.53817859231452,
            p_plasma_rad_mw=287.99550050743289,
            itart=0,
            p_plasma_separatrix_mw=143.03180561618876,
            p_fw_alpha_mw=19.833077403424262,
            n_divertors=1,
            p_plasma_ohmic_mw=0.61391840981850698,
            i_rad_loss=1,
            p_fusion_total_mw=1985.785106643267,
            p_non_alpha_charged_mw=1.6064693283140403,
            pscalingmw=325.08626176539281,
            f_p_alpha_plasma_deposited=0.94999999999999996,
            p_cp_coolant_pump_elec=0,
            i_tf_sup=1,
            tfcmw=0,
            temp_tf_cryo=4.5,
            tcoolin=313.14999999999998,
            eff_tf_cryo=0.13,
            p_fw_blkt_coolant_pump_mw=202.02739897087824,
            p_shld_coolant_pump_elec_mw=0.0078214713655548,
            p_div_coolant_pump_elec_mw=2.0623190688834718,
            p_coolant_pump_total_mw=203.8284212408949,
            p_fw_blkt_heat_deposited_mw=2258.2163348769122,
            p_fw_blkt_coolant_pump_elec_mw=232.21540111595198,
            p_div_heat_deposited_mw=360.63773557565275,
            p_fw_heat_deposited_mw=0,
            p_shld_heat_deposited_mw=1.3677406976945679,
            p_cp_coolant_pump_elec_mw=0,
            p_plant_core_systems_elec_mw=0,
            f_p_div_primary_heat=0.13763633828287813,
            qss=20361.633927097802,
            qac=3611.3456752656607,
            qcl=16108.2211128,
            qmisc=23850.540321823562,
            outfile=11,
            iprint=0,
            expected_p_plant_electric_net_mw=549.90044139,
            expected_precircmw=432.68273779,
            expected_fachtmw=5.0,
            expected_p_plant_electric_gross_mw=982.58317918134742,
            expected_p_plant_secondary_heat_mw=176.92004712,
            expected_pcoresystems=68.45107682927969,
        ),
        Power2Param(
            p_plant_electric_net_required_mw=500,
            ipnet=0,
            ireactor=1,
            p_hcd_injected_total_mw=51.978447720428512,
            p_blkt_multiplication_mw=377.8143718115644,
            inuclear=1,
            p_blkt_nuclear_heat_total_mw=1549.9285082739402,
            p_fw_rad_total_mw=254.87601794907812,
            qnuc=12920,
            eta_coolant_pump_electric=0.87000000000000011,
            f_p_blkt_multiplication=1.2690000534057617,
            p_div_rad_total_mw=33.119482558354782,
            f_ster_div_single=0.115,
            f_a_fw_outboard_hcd=0,
            i_thermal_electric_conversion=2,
            pnuc_cp=0,
            p_div_nuclear_heat_total_mw=182.6352084763719,
            i_p_coolant_pumping=3,
            p_tf_nuclear_heat_mw=0.045535131445547841,
            p_fw_hcd_nuclear_heat_mw=0,
            p_shld_nuclear_heat_mw=1.4036212304705389,
            p_fw_hcd_rad_total_mw=0,
            p_fw_nuclear_heat_total_mw=230.95082168283884,
            p_shld_coolant_pump_mw=0.0070181061523526943,
            p_blkt_coolant_pump_mw=0,
            p_shld_secondary_heat_mw=0,
            f_p_shld_coolant_pump_total_heat=0.0050000000000000001,
            temp_turbine_coolant_in=0,
            p_plant_electric_net_mw=493.01760776192009,
            f_p_div_coolant_pump_total_heat=0.0050000000000000001,
            f_p_blkt_coolant_pump_total_heat=0.0050000000000000001,
            vachtmw=0.5,
            p_div_coolant_pump_mw=1.7933419035282543,
            n_primary_heat_exchangers=3,
            helpow=220505.71684249729,
            p_fw_coolant_pump_mw=0,
            p_plant_electric_recirc_mw=489.9198817019128,
            p_plant_primary_heat_mw=2619.4223856129224,
            f_p_fw_coolant_pump_total_heat=0.0050000000000000001,
            p_plant_electric_base_total_mw=62.237143915360818,
            i_shld_primary_heat=1,
            p_hcd_electric_total_mw=129.94611930107126,
            fachtmw=61.882833632875375,
            p_plant_electric_gross_mw=982.58317918134742,
            p_plant_secondary_heat_mw=234.15719103660052,
            p_tritium_plant_electric_mw=15,
            p_hcd_secondary_heat_mw=0,
            p_tf_electric_supplies_mw=9.1507079104675704,
            p_coolant_pump_elec_total_mw=234.2162627659944,
            eta_turbine=0.37500000000000006,
            p_cryo_plant_electric_mw=108.74512702403499,
            p_div_secondary_heat_mw=0,
            p_hcd_electric_loss_mw=77.967671580642758,
            p_coolant_pump_loss_total_mw=30.448114159579291,
            helpow_cryal=0,
            p_pf_electric_supplies_mw=0.068213156646500808,
            p_alpha_total_mw=396.53774329057228,
            i_plasma_ignited=0,
            p_plasma_inner_rad_mw=113.53817859231452,
            p_plasma_rad_mw=287.99550050743289,
            itart=0,
            p_plasma_separatrix_mw=142.91368967092416,
            p_fw_alpha_mw=19.826887164528632,
            n_divertors=1,
            p_plasma_ohmic_mw=0.61391840981850698,
            i_rad_loss=1,
            p_fusion_total_mw=1985.1653095257811,
            p_non_alpha_charged_mw=1.6059679220663614,
            pscalingmw=325.00280675287695,
            f_p_alpha_plasma_deposited=0.94999999999999996,
            p_cp_coolant_pump_elec=0,
            i_tf_sup=1,
            tfcmw=0,
            temp_tf_cryo=4.5,
            tcoolin=313.14999999999998,
            eff_tf_cryo=0.13,
            p_fw_blkt_coolant_pump_mw=201.96778859673452,
            p_shld_coolant_pump_elec_mw=0.0080667886808651647,
            p_div_coolant_pump_elec_mw=2.0613125327910966,
            p_coolant_pump_total_mw=203.76814860641511,
            p_fw_blkt_heat_deposited_mw=2257.5500236671205,
            p_fw_blkt_coolant_pump_elec_mw=232.14688344452242,
            p_div_heat_deposited_mw=360.46172260917911,
            p_fw_heat_deposited_mw=0,
            p_shld_heat_deposited_mw=1.4106393366228915,
            p_cp_coolant_pump_elec_mw=0,
            p_plant_core_systems_elec_mw=125.68822074464052,
            f_p_div_primary_heat=0.13761114839248584,
            qss=20342.863776957758,
            qac=102701.82327748176,
            qcl=16108.2211128,
            qmisc=68432.80867525778,
            outfile=11,
            iprint=0,
            expected_p_plant_electric_net_mw=479.65696445,
            expected_precircmw=502.62643016,
            expected_fachtmw=5.0,
            expected_p_plant_electric_gross_mw=982.28339460484608,
            expected_p_plant_secondary_heat_mw=246.92536896,
            expected_pcoresystems=138.46404809114904,
        ),
    ),
)
def test_power2(power2param, monkeypatch, power):
    """
    Automatically generated Regression Unit Test for power2.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param power2param: the data used to mock and assert in this test.
    :type power2param: power2param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        constraint_variables,
        "p_plant_electric_net_required_mw",
        power2param.p_plant_electric_net_required_mw,
    )

    monkeypatch.setattr(cost_variables, "ipnet", power2param.ipnet)

    monkeypatch.setattr(cost_variables, "ireactor", power2param.ireactor)

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_injected_total_mw",
        power2param.p_hcd_injected_total_mw,
    )

    monkeypatch.setattr(
        fwbs_variables, "p_blkt_multiplication_mw", power2param.p_blkt_multiplication_mw
    )

    monkeypatch.setattr(fwbs_variables, "inuclear", power2param.inuclear)

    monkeypatch.setattr(
        fwbs_variables,
        "p_blkt_nuclear_heat_total_mw",
        power2param.p_blkt_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(
        fwbs_variables, "p_fw_rad_total_mw", power2param.p_fw_rad_total_mw
    )

    monkeypatch.setattr(fwbs_variables, "qnuc", power2param.qnuc)

    monkeypatch.setattr(
        fwbs_variables,
        "eta_coolant_pump_electric",
        power2param.eta_coolant_pump_electric,
    )

    monkeypatch.setattr(
        fwbs_variables, "f_p_blkt_multiplication", power2param.f_p_blkt_multiplication
    )

    monkeypatch.setattr(
        fwbs_variables, "p_div_rad_total_mw", power2param.p_div_rad_total_mw
    )

    monkeypatch.setattr(
        fwbs_variables, "f_ster_div_single", power2param.f_ster_div_single
    )

    monkeypatch.setattr(
        fwbs_variables, "f_a_fw_outboard_hcd", power2param.f_a_fw_outboard_hcd
    )

    monkeypatch.setattr(
        fwbs_variables,
        "i_thermal_electric_conversion",
        power2param.i_thermal_electric_conversion,
    )

    monkeypatch.setattr(fwbs_variables, "pnuc_cp", power2param.pnuc_cp)

    monkeypatch.setattr(
        fwbs_variables,
        "p_div_nuclear_heat_total_mw",
        power2param.p_div_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(
        fwbs_variables, "i_p_coolant_pumping", power2param.i_p_coolant_pumping
    )

    monkeypatch.setattr(
        fwbs_variables, "p_tf_nuclear_heat_mw", power2param.p_tf_nuclear_heat_mw
    )

    monkeypatch.setattr(
        fwbs_variables, "p_fw_hcd_nuclear_heat_mw", power2param.p_fw_hcd_nuclear_heat_mw
    )

    monkeypatch.setattr(
        fwbs_variables, "p_shld_nuclear_heat_mw", power2param.p_shld_nuclear_heat_mw
    )

    monkeypatch.setattr(
        fwbs_variables, "p_fw_hcd_rad_total_mw", power2param.p_fw_hcd_rad_total_mw
    )

    monkeypatch.setattr(
        fwbs_variables,
        "p_fw_nuclear_heat_total_mw",
        power2param.p_fw_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_shld_coolant_pump_mw",
        power2param.p_shld_coolant_pump_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_blkt_coolant_pump_mw",
        power2param.p_blkt_coolant_pump_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_shld_secondary_heat_mw",
        power2param.p_shld_secondary_heat_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "f_p_shld_coolant_pump_total_heat",
        power2param.f_p_shld_coolant_pump_total_heat,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "temp_turbine_coolant_in",
        power2param.temp_turbine_coolant_in,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_net_mw",
        power2param.p_plant_electric_net_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "f_p_div_coolant_pump_total_heat",
        power2param.f_p_div_coolant_pump_total_heat,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "f_p_blkt_coolant_pump_total_heat",
        power2param.f_p_blkt_coolant_pump_total_heat,
    )

    monkeypatch.setattr(heat_transport_variables, "vachtmw", power2param.vachtmw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_div_coolant_pump_mw",
        power2param.p_div_coolant_pump_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "n_primary_heat_exchangers",
        power2param.n_primary_heat_exchangers,
    )

    monkeypatch.setattr(heat_transport_variables, "helpow", power2param.helpow)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_fw_coolant_pump_mw",
        power2param.p_fw_coolant_pump_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_recirc_mw",
        power2param.p_plant_electric_recirc_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_primary_heat_mw",
        power2param.p_plant_primary_heat_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "f_p_fw_coolant_pump_total_heat",
        power2param.f_p_fw_coolant_pump_total_heat,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_base_total_mw",
        power2param.p_plant_electric_base_total_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables, "i_shld_primary_heat", power2param.i_shld_primary_heat
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_total_mw",
        power2param.p_hcd_electric_total_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "fachtmw", power2param.fachtmw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_gross_mw",
        power2param.p_plant_electric_gross_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_secondary_heat_mw",
        power2param.p_plant_secondary_heat_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_tritium_plant_electric_mw",
        power2param.p_tritium_plant_electric_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_secondary_heat_mw",
        power2param.p_hcd_secondary_heat_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_tf_electric_supplies_mw",
        power2param.p_tf_electric_supplies_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_coolant_pump_elec_total_mw",
        power2param.p_coolant_pump_elec_total_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "eta_turbine", power2param.eta_turbine)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_cryo_plant_electric_mw",
        power2param.p_cryo_plant_electric_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_div_secondary_heat_mw",
        power2param.p_div_secondary_heat_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_loss_mw",
        power2param.p_hcd_electric_loss_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_coolant_pump_loss_total_mw",
        power2param.p_coolant_pump_loss_total_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables, "helpow_cryal", power2param.helpow_cryal
    )

    monkeypatch.setattr(
        pfcoil_variables,
        "p_pf_electric_supplies_mw",
        power2param.p_pf_electric_supplies_mw,
    )

    monkeypatch.setattr(
        physics_variables, "p_alpha_total_mw", power2param.p_alpha_total_mw
    )

    monkeypatch.setattr(
        physics_variables, "i_plasma_ignited", power2param.i_plasma_ignited
    )

    monkeypatch.setattr(
        physics_variables, "p_plasma_inner_rad_mw", power2param.p_plasma_inner_rad_mw
    )

    monkeypatch.setattr(
        physics_variables, "p_plasma_rad_mw", power2param.p_plasma_rad_mw
    )

    monkeypatch.setattr(physics_variables, "itart", power2param.itart)

    monkeypatch.setattr(
        physics_variables, "p_plasma_separatrix_mw", power2param.p_plasma_separatrix_mw
    )

    monkeypatch.setattr(physics_variables, "p_fw_alpha_mw", power2param.p_fw_alpha_mw)

    monkeypatch.setattr(physics_variables, "n_divertors", power2param.n_divertors)

    monkeypatch.setattr(
        physics_variables, "p_plasma_ohmic_mw", power2param.p_plasma_ohmic_mw
    )

    monkeypatch.setattr(physics_variables, "i_rad_loss", power2param.i_rad_loss)

    monkeypatch.setattr(
        physics_variables, "p_fusion_total_mw", power2param.p_fusion_total_mw
    )

    monkeypatch.setattr(
        physics_variables,
        "p_non_alpha_charged_mw",
        power2param.p_non_alpha_charged_mw,
    )

    monkeypatch.setattr(physics_variables, "pscalingmw", power2param.pscalingmw)

    monkeypatch.setattr(
        physics_variables,
        "f_p_alpha_plasma_deposited",
        power2param.f_p_alpha_plasma_deposited,
    )

    monkeypatch.setattr(
        tfcoil_variables, "p_cp_coolant_pump_elec", power2param.p_cp_coolant_pump_elec
    )

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", power2param.i_tf_sup)

    monkeypatch.setattr(tfcoil_variables, "tfcmw", power2param.tfcmw)

    monkeypatch.setattr(tfcoil_variables, "temp_tf_cryo", power2param.temp_tf_cryo)

    monkeypatch.setattr(tfcoil_variables, "tcoolin", power2param.tcoolin)

    monkeypatch.setattr(tfcoil_variables, "eff_tf_cryo", power2param.eff_tf_cryo)

    monkeypatch.setattr(
        ppv, "p_fw_blkt_coolant_pump_mw", power2param.p_fw_blkt_coolant_pump_mw
    )

    monkeypatch.setattr(
        power_variables,
        "p_shld_coolant_pump_elec_mw",
        power2param.p_shld_coolant_pump_elec_mw,
    )

    monkeypatch.setattr(
        power_variables,
        "p_div_coolant_pump_elec_mw",
        power2param.p_div_coolant_pump_elec_mw,
    )

    monkeypatch.setattr(
        power_variables, "p_coolant_pump_total_mw", power2param.p_coolant_pump_total_mw
    )

    monkeypatch.setattr(
        power_variables,
        "p_fw_blkt_heat_deposited_mw",
        power2param.p_fw_blkt_heat_deposited_mw,
    )

    monkeypatch.setattr(
        power_variables,
        "p_fw_blkt_coolant_pump_elec_mw",
        power2param.p_fw_blkt_coolant_pump_elec_mw,
    )

    monkeypatch.setattr(
        power_variables, "p_div_heat_deposited_mw", power2param.p_div_heat_deposited_mw
    )

    monkeypatch.setattr(
        power_variables, "p_fw_heat_deposited_mw", power2param.p_fw_heat_deposited_mw
    )

    monkeypatch.setattr(
        power_variables,
        "p_shld_heat_deposited_mw",
        power2param.p_shld_heat_deposited_mw,
    )

    monkeypatch.setattr(
        power_variables,
        "p_cp_coolant_pump_elec_mw",
        power2param.p_cp_coolant_pump_elec_mw,
    )

    monkeypatch.setattr(
        power_variables,
        "p_plant_core_systems_elec_mw",
        power2param.p_plant_core_systems_elec_mw,
    )

    monkeypatch.setattr(
        power_variables, "f_p_div_primary_heat", power2param.f_p_div_primary_heat
    )

    monkeypatch.setattr(power_variables, "qss", power2param.qss)

    monkeypatch.setattr(power_variables, "qac", power2param.qac)

    monkeypatch.setattr(power_variables, "qcl", power2param.qcl)

    monkeypatch.setattr(power_variables, "qmisc", power2param.qmisc)

    power.plant_electric_production()

    assert heat_transport_variables.p_plant_electric_net_mw == pytest.approx(
        power2param.expected_p_plant_electric_net_mw
    )

    assert heat_transport_variables.p_plant_electric_recirc_mw == pytest.approx(
        power2param.expected_precircmw
    )

    assert heat_transport_variables.fachtmw == pytest.approx(
        power2param.expected_fachtmw
    )

    assert heat_transport_variables.p_plant_electric_gross_mw == pytest.approx(
        power2param.expected_p_plant_electric_gross_mw
    )

    assert heat_transport_variables.p_plant_secondary_heat_mw == pytest.approx(
        power2param.expected_p_plant_secondary_heat_mw
    )

    assert power_variables.p_plant_core_systems_elec_mw == pytest.approx(
        power2param.expected_pcoresystems
    )
