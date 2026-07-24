from typing import Any, NamedTuple

import numpy as np
import pytest


@pytest.fixture
def power(process_models):
    """Fixture to get the Power instance from process_models.

    :returns: initialised power object
    :rtype: process.power.Power
    """
    return process_models.power


class CryoParam(NamedTuple):
    qnuc: Any = None

    inuclear: Any = None

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
    [
        CryoParam(
            qnuc=12920,
            inuclear=1,
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
    ],
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

    monkeypatch.setattr(power.data.fwbs, "qnuc", cryoparam.qnuc)

    monkeypatch.setattr(power.data.fwbs, "inuclear", cryoparam.inuclear)

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

    assert power.data.power.qss == pytest.approx(cryoparam.expected_qss)

    assert power.data.power.qac == pytest.approx(cryoparam.expected_qac)

    assert power.data.power.qcl == pytest.approx(cryoparam.expected_qcl)

    assert power.data.power.qmisc == pytest.approx(cryoparam.expected_qmisc)

    assert helpow == pytest.approx(cryoparam.expected_helpow)


class PfpwrParam(NamedTuple):
    iohcl: Any = None

    poloidalpower: Any = None

    n_pf_coil_groups: Any = None

    c_pf_coil_turn: Any = None

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
    [
        PfpwrParam(
            iohcl=1,
            poloidalpower=np.zeros(5),
            n_pf_coil_groups=4,
            c_pf_coil_turn=np.array(
                (
                    np.zeros(22),
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
                    np.zeros(22),
                ),
                order="F",
            ).transpose(),
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
                    *np.zeros((14, 22)),
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
            f_a_pf_coil_void=np.full(22, 0.29999999999999999),
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
            poloidalpower=np.array(
                np.array(
                    (59332953.082890816, 43806300.444207191, 0, 0, -211211992.31967318),
                    order="F",
                ),
                order="F",
            ).transpose(),
            n_pf_coil_groups=4,
            c_pf_coil_turn=np.array(
                (
                    np.zeros(22),
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
                    np.zeros(22),
                ),
                order="F",
            ).transpose(),
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
                    *np.zeros((14, 22)),
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
            f_a_pf_coil_void=np.full(22, 0.29999999999999999),
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
    ],
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

    monkeypatch.setattr(power.data.build, "iohcl", pfpwrparam.iohcl)

    monkeypatch.setattr(power.data.pf_power, "poloidalpower", pfpwrparam.poloidalpower)

    for field in [
        "n_pf_coil_groups",
        "c_pf_coil_turn",
        "rho_pf_coil",
        "n_pf_cs_plasma_circuits",
        "n_pf_coils_in_group",
        "c_pf_cs_coils_peak_ma",
        "etapsu",
        "c_pf_coil_turn_peak_input",
        "c_pf_cs_coil_pulse_end_ma",
        "ind_pf_cs_plasma_mutual",
        "n_pf_coil_turns",
        "f_a_pf_coil_void",
        "j_pf_coil_wp_peak",
        "r_pf_coil_middle",
    ]:
        monkeypatch.setattr(power.data.pf_coil, field, getattr(pfpwrparam, field))

    monkeypatch.setattr(
        power.data.physics, "p_plasma_ohmic_mw", pfpwrparam.p_plasma_ohmic_mw
    )

    monkeypatch.setattr(power.data.physics, "rmajor", pfpwrparam.rmajor)

    monkeypatch.setattr(
        power.data.numerics, "active_constraints", pfpwrparam.active_constraints
    )

    monkeypatch.setattr(power.data.numerics, "ioptimz", pfpwrparam.ioptimz)

    monkeypatch.setattr(
        power.data.times, "t_pulse_cumulative", pfpwrparam.t_pulse_cumulative
    )

    monkeypatch.setattr(
        power.data.times,
        "t_plant_pulse_plasma_current_ramp_up",
        pfpwrparam.t_plant_pulse_plasma_current_ramp_up,
    )

    power.pfpwr(output=False)

    assert power.data.heat_transport.peakmva == pytest.approx(
        pfpwrparam.expected_peakmva
    )

    assert power.data.pf_power.pfckts == pytest.approx(pfpwrparam.expected_pfckts)

    assert power.data.pf_power.peakpoloidalpower == pytest.approx(
        pfpwrparam.expected_peakpoloidalpower
    )

    assert power.data.pf_power.spfbusl == pytest.approx(pfpwrparam.expected_spfbusl)

    assert power.data.pf_power.poloidalpower == pytest.approx(
        pfpwrparam.expected_poloidalpower
    )

    assert power.data.pf_power.spsmva == pytest.approx(pfpwrparam.expected_spsmva)

    assert power.data.pf_power.vpfskv == pytest.approx(pfpwrparam.expected_vpfskv)

    assert power.data.pf_power.ensxpfm == pytest.approx(pfpwrparam.expected_ensxpfm)

    assert power.data.pf_power.acptmax == pytest.approx(pfpwrparam.expected_acptmax)

    assert power.data.pf_power.srcktpm == pytest.approx(pfpwrparam.expected_srcktpm)


class AcpowParam(NamedTuple):
    p_cryo_plant_electric_mw: Any = None

    vachtmw: Any = None

    p_tf_electric_supplies_mw: Any = None

    p_tritium_plant_electric_mw: Any = None

    p_hcd_electric_total_mw: Any = None

    peakmva: Any = None

    p_plant_electric_base_total_mw: Any = None

    fmgdmw: Any = None

    p_coolant_pump_elec_total_mw: Any = None

    i_pf_energy_storage_source: Any = None

    srcktpm: Any = None

    expected_pacpmw: Any = None


@pytest.mark.parametrize(
    "acpowparam",
    [
        AcpowParam(
            p_cryo_plant_electric_mw=37.900388528497025,
            vachtmw=0.5,
            p_tf_electric_supplies_mw=9.1507079104675704,
            p_tritium_plant_electric_mw=15,
            p_hcd_electric_total_mw=129.94611930107126,
            peakmva=736.39062584245937,
            p_plant_electric_base_total_mw=0,
            fmgdmw=0,
            p_coolant_pump_elec_total_mw=234.28554165620102,
            i_pf_energy_storage_source=2,
            srcktpm=1071.1112934857531,
            expected_pacpmw=1164.244494532182,
        ),
        AcpowParam(
            p_cryo_plant_electric_mw=108.74512702403499,
            vachtmw=0.5,
            p_tf_electric_supplies_mw=9.1507079104675704,
            p_tritium_plant_electric_mw=15,
            p_hcd_electric_total_mw=129.94611930107126,
            peakmva=90.673341440806112,
            p_plant_electric_base_total_mw=62.23714391536082,
            fmgdmw=0,
            p_coolant_pump_elec_total_mw=234.2162627659944,
            i_pf_energy_storage_source=2,
            srcktpm=1069.8879533693198,
            expected_pacpmw=589.3014463957436,
        ),
    ],
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

    for field in [
        "p_cryo_plant_electric_mw",
        "vachtmw",
        "p_tf_electric_supplies_mw",
        "p_tritium_plant_electric_mw",
        "p_hcd_electric_total_mw",
        "peakmva",
        "p_plant_electric_base_total_mw",
        "fmgdmw",
        "p_coolant_pump_elec_total_mw",
    ]:
        monkeypatch.setattr(
            power.data.heat_transport,
            field,
            getattr(acpowparam, field),
        )

    for field in [
        "i_pf_energy_storage_source",
        "srcktpm",
    ]:
        monkeypatch.setattr(
            power.data.pf_power,
            field,
            getattr(acpowparam, field),
        )

    power.acpow(output=False)

    assert power.data.heat_transport.pacpmw == pytest.approx(acpowparam.expected_pacpmw)


class PlantElectricProductionParam(NamedTuple):
    ireactor: Any = None

    i_p_coolant_pumping: Any = None

    p_tf_nuclear_heat_mw: Any = None

    p_shld_secondary_heat_mw: Any = None

    vachtmw: Any = None

    p_plant_primary_heat_mw: Any = None

    p_hcd_electric_total_mw: Any = None

    p_tritium_plant_electric_mw: Any = None

    p_hcd_secondary_heat_mw: Any = None

    p_tf_electric_supplies_mw: Any = None

    p_coolant_pump_elec_total_mw: Any = None

    eta_turbine: Any = None

    p_cryo_plant_electric_mw: Any = None

    p_div_secondary_heat_mw: Any = None

    p_hcd_electric_loss_mw: Any = None

    p_coolant_pump_loss_total_mw: Any = None

    p_pf_electric_supplies_mw: Any = None

    itart: Any = None

    p_fusion_total_mw: Any = None

    p_cp_coolant_pump_elec: Any = None

    i_tf_sup: Any = None

    expected_p_plant_electric_net_mw: Any = None

    expected_precircmw: Any = None

    expected_fachtmw: Any = None

    expected_p_plant_electric_gross_mw: Any = None

    expected_p_plant_secondary_heat_mw: Any = None

    expected_pcoresystems: Any = None


@pytest.mark.parametrize(
    "plantelecprodparam",
    [
        PlantElectricProductionParam(
            ireactor=1,
            i_p_coolant_pumping=3,
            p_tf_nuclear_heat_mw=0.044178296011112193,
            p_shld_secondary_heat_mw=0,
            vachtmw=0.5,
            p_plant_primary_heat_mw=2620.2218111502593,
            p_hcd_electric_total_mw=129.94611930107126,
            p_tritium_plant_electric_mw=15,
            p_hcd_secondary_heat_mw=0,
            p_tf_electric_supplies_mw=9.1507079104675704,
            p_coolant_pump_elec_total_mw=234.28554165620102,
            eta_turbine=0.37500000000000006,
            p_cryo_plant_electric_mw=37.900388528497025,
            p_div_secondary_heat_mw=0,
            p_hcd_electric_loss_mw=77.967671580642758,
            p_coolant_pump_loss_total_mw=30.457120415306122,
            p_pf_electric_supplies_mw=0.89998039031509891,
            itart=0,
            p_fusion_total_mw=1985.785106643267,
            p_cp_coolant_pump_elec=0,
            i_tf_sup=1,
            expected_p_plant_electric_net_mw=549.90044139,
            expected_precircmw=432.68273779,
            expected_fachtmw=5.0,
            expected_p_plant_electric_gross_mw=982.58317918134742,
            expected_p_plant_secondary_heat_mw=176.92004712,
            expected_pcoresystems=68.45107682927969,
        ),
        PlantElectricProductionParam(
            ireactor=1,
            i_p_coolant_pumping=3,
            p_tf_nuclear_heat_mw=0.045535131445547841,
            p_shld_secondary_heat_mw=0,
            vachtmw=0.5,
            p_plant_primary_heat_mw=2619.4223856129224,
            p_hcd_electric_total_mw=129.94611930107126,
            p_tritium_plant_electric_mw=15,
            p_hcd_secondary_heat_mw=0,
            p_tf_electric_supplies_mw=9.1507079104675704,
            p_coolant_pump_elec_total_mw=234.2162627659944,
            eta_turbine=0.37500000000000006,
            p_cryo_plant_electric_mw=108.74512702403499,
            p_div_secondary_heat_mw=0,
            p_hcd_electric_loss_mw=77.967671580642758,
            p_coolant_pump_loss_total_mw=30.448114159579291,
            p_pf_electric_supplies_mw=0.068213156646500808,
            itart=0,
            p_fusion_total_mw=1985.1653095257811,
            p_cp_coolant_pump_elec=0,
            i_tf_sup=1,
            expected_p_plant_electric_net_mw=479.65696445,
            expected_precircmw=502.62643016,
            expected_fachtmw=5.0,
            expected_p_plant_electric_gross_mw=982.28339460484608,
            expected_p_plant_secondary_heat_mw=246.92536896,
            expected_pcoresystems=138.46404809114904,
        ),
    ],
)
def test_plant_electric_production(plantelecprodparam, monkeypatch, power):
    """
    Automatically generated Regression Unit Test for plant_electric_production.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param power2param: the data used to mock and assert in this test.
    :type power2param: power2param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(power.data.costs, "ireactor", plantelecprodparam.ireactor)

    for field in [
        "i_p_coolant_pumping",
        "p_tf_nuclear_heat_mw",
    ]:
        monkeypatch.setattr(power.data.fwbs, field, getattr(plantelecprodparam, field))

    for field in [
        "p_shld_secondary_heat_mw",
        "vachtmw",
        "p_plant_primary_heat_mw",
        "p_hcd_electric_total_mw",
        "p_tritium_plant_electric_mw",
        "p_hcd_secondary_heat_mw",
        "p_tf_electric_supplies_mw",
        "p_coolant_pump_elec_total_mw",
        "eta_turbine",
        "p_cryo_plant_electric_mw",
        "p_div_secondary_heat_mw",
        "p_hcd_electric_loss_mw",
        "p_coolant_pump_loss_total_mw",
    ]:
        monkeypatch.setattr(
            power.data.heat_transport,
            field,
            getattr(plantelecprodparam, field),
        )

    monkeypatch.setattr(
        power.data.pf_coil,
        "p_pf_electric_supplies_mw",
        plantelecprodparam.p_pf_electric_supplies_mw,
    )

    monkeypatch.setattr(power.data.physics, "itart", plantelecprodparam.itart)

    monkeypatch.setattr(
        power.data.physics, "p_fusion_total_mw", plantelecprodparam.p_fusion_total_mw
    )

    for field in [
        "p_cp_coolant_pump_elec",
        "i_tf_sup",
    ]:
        monkeypatch.setattr(
            power.data.tfcoil,
            field,
            getattr(plantelecprodparam, field),
        )

    power.plant_electric_production()

    assert power.data.heat_transport.p_plant_electric_net_mw == pytest.approx(
        plantelecprodparam.expected_p_plant_electric_net_mw
    )

    assert power.data.heat_transport.p_plant_electric_recirc_mw == pytest.approx(
        plantelecprodparam.expected_precircmw
    )

    assert power.data.heat_transport.fachtmw == pytest.approx(
        plantelecprodparam.expected_fachtmw
    )

    assert power.data.heat_transport.p_plant_electric_gross_mw == pytest.approx(
        plantelecprodparam.expected_p_plant_electric_gross_mw
    )

    assert power.data.heat_transport.p_plant_secondary_heat_mw == pytest.approx(
        plantelecprodparam.expected_p_plant_secondary_heat_mw
    )

    assert power.data.power.p_plant_core_systems_elec_mw == pytest.approx(
        plantelecprodparam.expected_pcoresystems
    )
