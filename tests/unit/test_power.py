from typing import Any, NamedTuple

import numpy
import pytest

from process.fortran import (
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
    tfcoil_variables,
    times_variables,
)
from process.fortran import primary_pumping_variables as ppv
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

    m_components_cryo_cooled: Any = None

    cpttf: Any = None

    ensxpfm: Any = None

    p_tf_nuclear_heat_mw: Any = None

    n_tf: Any = None

    a_tf_cryo: Any = None

    t_pulse_repetition: Any = None

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
            m_components_cryo_cooled=47352637.039762333,
            cpttf=74026.751437500003,
            ensxpfm=37429.525515086898,
            p_tf_nuclear_heat_mw=0.044178296011112193,
            n_tf=16,
            a_tf_cryo=0,
            t_pulse_repetition=10364.426139387357,
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
            m_components_cryo_cooled=47308985.527808741,
            cpttf=74026.751437500003,
            ensxpfm=37427.228965055205,
            p_tf_nuclear_heat_mw=0.045535131445547841,
            n_tf=16,
            a_tf_cryo=0,
            t_pulse_repetition=364.42613938735633,
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

    monkeypatch.setattr(power, "qss", cryoparam.qss)

    monkeypatch.setattr(power, "qac", cryoparam.qac)

    monkeypatch.setattr(power, "qcl", cryoparam.qcl)

    monkeypatch.setattr(power, "qmisc", cryoparam.qmisc)

    helpow = power.cryo(
        i_tf_sup=cryoparam.i_tf_sup,
        m_components_cryo_cooled=cryoparam.m_components_cryo_cooled,
        cpttf=cryoparam.cpttf,
        ensxpfm=cryoparam.ensxpfm,
        p_tf_nuclear_heat_mw=cryoparam.p_tf_nuclear_heat_mw,
        n_tf=cryoparam.n_tf,
        a_tf_cryo=cryoparam.a_tf_cryo,
        t_pulse_repetition=cryoparam.t_pulse_repetition,
    )

    assert power.qss == pytest.approx(cryoparam.expected_qss)

    assert power.qac == pytest.approx(cryoparam.expected_qac)

    assert power.qcl == pytest.approx(cryoparam.expected_qcl)

    assert power.qmisc == pytest.approx(cryoparam.expected_qmisc)

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

    p_pf_resisitve_total_kw: Any = None

    ngrp: Any = None

    cpt: Any = None

    pfwpmw: Any = None

    pfclres: Any = None

    ncirt: Any = None

    ncls: Any = None

    ric: Any = None

    etapsu: Any = None

    cptdin: Any = None

    curpfs: Any = None

    sxlg: Any = None

    turns: Any = None

    vf: Any = None

    rjconpf: Any = None

    rpf: Any = None

    p_plasma_ohmic_mw: Any = None

    rmajor: Any = None

    active_constraints: Any = None

    ioptimz: Any = None

    tim: Any = None

    intervallabel: Any = None

    timelabel: Any = None

    t_current_ramp_up: Any = None

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

    expected_p_pf_resisitve_total_kw: Any = None


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
            poloidalpower=numpy.array(
                numpy.array((0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            spsmva=0,
            vpfskv=0,
            ensxpfm=0,
            acptmax=0,
            p_pf_resisitve_total_kw=0,
            ngrp=4,
            cpt=numpy.array(
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
            pfwpmw=0,
            pfclres=0,
            ncirt=8,
            ncls=numpy.array(
                numpy.array((1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            ric=numpy.array(
                numpy.array(
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
            cptdin=numpy.array(
                numpy.array(
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
            curpfs=numpy.array(
                numpy.array(
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
            sxlg=numpy.array(
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
            turns=numpy.array(
                numpy.array(
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
            vf=numpy.array(
                numpy.array(
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
            rjconpf=numpy.array(
                numpy.array(
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
            rpf=numpy.array(
                numpy.array(
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
            ),
            ioptimz=1,
            tim=numpy.array(
                numpy.array(
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
                "t_precharge      ",
                "t_current_ramp_up       ",
                "t_fusion_ramp      ",
                "t_burn      ",
                "t_ramp_down      ",
            ),
            timelabel=(
                "Start      ",
                "BOP        ",
                "EOR        ",
                "BOF        ",
                "EOF        ",
                "EOP        ",
            ),
            t_current_ramp_up=177.21306969367816,
            outfile=11,
            iprint=0,
            expected_peakmva=736.39062584245937,
            expected_pfckts=12,
            expected_peakpoloidalpower=211.21199231967319,
            expected_spfbusl=2533.4495999999999,
            expected_poloidalpower=numpy.array(
                numpy.array(
                    (59332953.082890816, 43806300.444207191, 0, 0, -211211992.31967318),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_spsmva=845.66824574150155,
            expected_vpfskv=20,
            expected_ensxpfm=37429.525515086898,
            expected_acptmax=24.816666666666666,
            expected_p_pf_resisitve_total_kw=1071.1112934857531,
        ),
        PfpwrParam(
            iohcl=1,
            peakmva=736.39062584245937,
            pfckts=12,
            maxpoloidalpower=1000,
            peakpoloidalpower=211.21199231967319,
            spfbusl=2533.4495999999999,
            poloidalpower=numpy.array(
                numpy.array(
                    (59332953.082890816, 43806300.444207191, 0, 0, -211211992.31967318),
                    order="F",
                ),
                order="F",
            ).transpose(),
            spsmva=845.66824574150155,
            vpfskv=20,
            ensxpfm=37429.525515086898,
            acptmax=24.816666666666666,
            p_pf_resisitve_total_kw=1071.1112934857531,
            ngrp=4,
            cpt=numpy.array(
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
            pfwpmw=0.89998039031509891,
            pfclres=0,
            ncirt=8,
            ncls=numpy.array(
                numpy.array((1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            ric=numpy.array(
                numpy.array(
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
            cptdin=numpy.array(
                numpy.array(
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
            curpfs=numpy.array(
                numpy.array(
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
            sxlg=numpy.array(
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
            turns=numpy.array(
                numpy.array(
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
            vf=numpy.array(
                numpy.array(
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
            rjconpf=numpy.array(
                numpy.array(
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
            rpf=numpy.array(
                numpy.array(
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
            ),
            ioptimz=1,
            tim=numpy.array(
                numpy.array(
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
                "t_precharge      ",
                "t_current_ramp_up       ",
                "t_fusion_ramp      ",
                "t_burn      ",
                "t_ramp_down      ",
            ),
            timelabel=(
                "Start      ",
                "BOP        ",
                "EOR        ",
                "BOF        ",
                "EOF        ",
                "EOP        ",
            ),
            t_current_ramp_up=177.21306969367816,
            outfile=11,
            iprint=0,
            expected_peakmva=90.673341440806112,
            expected_pfckts=12,
            expected_peakpoloidalpower=9900,
            expected_spfbusl=2533.4495999999999,
            expected_poloidalpower=numpy.array(
                numpy.array(
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
            expected_p_pf_resisitve_total_kw=1069.8879533693198,
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

    monkeypatch.setattr(
        pf_power_variables,
        "p_pf_resisitve_total_kw",
        pfpwrparam.p_pf_resisitve_total_kw,
    )

    monkeypatch.setattr(pfcoil_variables, "ngrp", pfpwrparam.ngrp)

    monkeypatch.setattr(pfcoil_variables, "cpt", pfpwrparam.cpt)

    monkeypatch.setattr(pfcoil_variables, "pfwpmw", pfpwrparam.pfwpmw)

    monkeypatch.setattr(pfcoil_variables, "pfclres", pfpwrparam.pfclres)

    monkeypatch.setattr(pfcoil_variables, "ncirt", pfpwrparam.ncirt)

    monkeypatch.setattr(pfcoil_variables, "ncls", pfpwrparam.ncls)

    monkeypatch.setattr(pfcoil_variables, "ric", pfpwrparam.ric)

    monkeypatch.setattr(pfcoil_variables, "etapsu", pfpwrparam.etapsu)

    monkeypatch.setattr(pfcoil_variables, "cptdin", pfpwrparam.cptdin)

    monkeypatch.setattr(pfcoil_variables, "curpfs", pfpwrparam.curpfs)

    monkeypatch.setattr(pfcoil_variables, "sxlg", pfpwrparam.sxlg)

    monkeypatch.setattr(pfcoil_variables, "turns", pfpwrparam.turns)

    monkeypatch.setattr(pfcoil_variables, "vf", pfpwrparam.vf)

    monkeypatch.setattr(pfcoil_variables, "rjconpf", pfpwrparam.rjconpf)

    monkeypatch.setattr(pfcoil_variables, "rpf", pfpwrparam.rpf)

    monkeypatch.setattr(
        physics_variables, "p_plasma_ohmic_mw", pfpwrparam.p_plasma_ohmic_mw
    )

    monkeypatch.setattr(physics_variables, "rmajor", pfpwrparam.rmajor)

    monkeypatch.setattr(numerics, "active_constraints", pfpwrparam.active_constraints)

    monkeypatch.setattr(numerics, "ioptimz", pfpwrparam.ioptimz)

    monkeypatch.setattr(times_variables, "tim", pfpwrparam.tim)

    monkeypatch.setattr(
        times_variables, "t_current_ramp_up", pfpwrparam.t_current_ramp_up
    )

    power.pfpwr(output=False)

    assert heat_transport_variables.peakmva == pytest.approx(
        pfpwrparam.expected_peakmva
    )

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

    assert pf_power_variables.p_pf_resisitve_total_kw == pytest.approx(
        pfpwrparam.expected_p_pf_resisitve_total_kw
    )


class AcpowParam(NamedTuple):
    a_floor_total: Any = None

    p_baseload_electrical: Any = None

    p_cryo_plant_mw: Any = None

    p_vacuum_pumps_mw: Any = None

    p_tf_electrical_mw: Any = None

    trithtmw: Any = None

    p_hcd_electrical_mw: Any = None

    tlvpmw: Any = None

    peakmva: Any = None

    p_baseload_electrical_total_mw: Any = None

    fmgdmw: Any = None

    pwpm2: Any = None

    p_pump_cool_elec_total_mw: Any = None

    p_pulsed_power_total_mw: Any = None

    i_pf_power_source: Any = None

    p_pf_resisitve_total_kw: Any = None

    iprint: Any = None

    outfile: Any = None

    expected_p_pulsed_power_total_mw: Any = None


@pytest.mark.parametrize(
    "acpowparam",
    (
        AcpowParam(
            a_floor_total=379218.8908858358,
            p_baseload_electrical=5000000,
            p_cryo_plant_mw=37.900388528497025,
            p_vacuum_pumps_mw=0.5,
            p_tf_electrical_mw=9.1507079104675704,
            trithtmw=15,
            p_hcd_electrical_mw=129.94611930107126,
            tlvpmw=0,
            peakmva=736.39062584245937,
            p_baseload_electrical_total_mw=0,
            fmgdmw=0,
            pwpm2=150,
            p_pump_cool_elec_total_mw=234.28554165620102,
            p_pulsed_power_total_mw=0,
            i_pf_power_source=2,
            p_pf_resisitve_total_kw=1071.1112934857531,
            iprint=0,
            outfile=11,
            expected_p_pulsed_power_total_mw=1164.244494532182,
        ),
        AcpowParam(
            a_floor_total=381580.9594357388,
            p_baseload_electrical=5000000,
            p_cryo_plant_mw=108.74512702403499,
            p_vacuum_pumps_mw=0.5,
            p_tf_electrical_mw=9.1507079104675704,
            trithtmw=15,
            p_hcd_electrical_mw=129.94611930107126,
            tlvpmw=699.34943812129745,
            peakmva=90.673341440806112,
            p_baseload_electrical_total_mw=61.882833632875375,
            fmgdmw=0,
            pwpm2=150,
            p_pump_cool_elec_total_mw=234.2162627659944,
            p_pulsed_power_total_mw=1226.1273281650574,
            i_pf_power_source=2,
            p_pf_resisitve_total_kw=1069.8879533693198,
            iprint=0,
            outfile=11,
            expected_p_pulsed_power_total_mw=589.3014463957436,
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

    monkeypatch.setattr(buildings_variables, "a_floor_total", acpowparam.a_floor_total)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_baseload_electrical",
        acpowparam.p_baseload_electrical,
    )

    monkeypatch.setattr(
        heat_transport_variables, "p_cryo_plant_mw", acpowparam.p_cryo_plant_mw
    )

    monkeypatch.setattr(
        heat_transport_variables, "p_vacuum_pumps_mw", acpowparam.p_vacuum_pumps_mw
    )

    monkeypatch.setattr(
        heat_transport_variables, "p_tf_electrical_mw", acpowparam.p_tf_electrical_mw
    )

    monkeypatch.setattr(heat_transport_variables, "trithtmw", acpowparam.trithtmw)

    monkeypatch.setattr(
        heat_transport_variables, "p_hcd_electrical_mw", acpowparam.p_hcd_electrical_mw
    )

    monkeypatch.setattr(heat_transport_variables, "tlvpmw", acpowparam.tlvpmw)

    monkeypatch.setattr(heat_transport_variables, "peakmva", acpowparam.peakmva)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_baseload_electrical_total_mw",
        acpowparam.p_baseload_electrical_total_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "fmgdmw", acpowparam.fmgdmw)

    monkeypatch.setattr(heat_transport_variables, "pwpm2", acpowparam.pwpm2)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_pump_cool_elec_total_mw",
        acpowparam.p_pump_cool_elec_total_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_pulsed_power_total_mw",
        acpowparam.p_pulsed_power_total_mw,
    )

    monkeypatch.setattr(
        pf_power_variables, "i_pf_power_source", acpowparam.i_pf_power_source
    )

    monkeypatch.setattr(
        pf_power_variables,
        "p_pf_resisitve_total_kw",
        acpowparam.p_pf_resisitve_total_kw,
    )

    power.acpow(output=False)

    assert heat_transport_variables.p_pulsed_power_total_mw == pytest.approx(
        acpowparam.expected_p_pulsed_power_total_mw
    )


class Power2Param(NamedTuple):
    pnetelin: Any = None

    ipnet: Any = None

    ireactor: Any = None

    pinjmw: Any = None

    emultmw: Any = None

    inuclear: Any = None

    p_blanket_nuclear_heat_mw: Any = None

    p_fw_radiation_mw: Any = None

    qnuc: Any = None

    eta_pump_coolant_electrical: Any = None

    emult: Any = None

    p_div_radiation_mw: Any = None

    fdiv: Any = None

    fhcd: Any = None

    secondary_cycle: Any = None

    pnuc_cp: Any = None

    p_div_nuclear_heat_mw: Any = None

    primary_pumping: Any = None

    p_tf_nuclear_heat_mw: Any = None

    p_hcd_nuclear_heat_mw: Any = None

    p_shield_nuclear_heat_mw: Any = None

    p_hcd_radiation_mw: Any = None

    p_fw_nuclear_heat_mw: Any = None

    p_shield_pump_cool_mw: Any = None

    p_blanket_pumping_mw: Any = None

    p_shield_thermal_secondary_mw: Any = None

    fpumpshld: Any = None

    tturb: Any = None

    p_net_electrical_mw: Any = None

    fpumpdiv: Any = None

    fpumpblkt: Any = None

    p_vacuum_pumps_mw: Any = None

    p_div_pump_cool_mw: Any = None

    nphx: Any = None

    helpow: Any = None

    p_fw_pumping_mw: Any = None

    p_recirc_electrical_mw: Any = None

    p_thermal_primary_mw: Any = None

    fpumpfw: Any = None

    p_baseload_electrical_total_mw: Any = None

    i_shield_power_generation: Any = None

    p_hcd_electrical_mw: Any = None

    fachtmw: Any = None

    p_gross_electrical: Any = None

    p_thermal_secondary_mw: Any = None

    trithtmw: Any = None

    p_hcd_thermal_secondary_mw: Any = None

    p_tf_electrical_mw: Any = None

    p_pump_cool_elec_total_mw: Any = None

    eta_thermal_electric: Any = None

    p_cryo_plant_mw: Any = None

    p_div_thermal_secondary_mw: Any = None

    p_hcd_electrical_loss_mw: Any = None

    p_pump_cool_loss_total_mw: Any = None

    helpow_cryal: Any = None

    pfwpmw: Any = None

    alpha_power_total: Any = None

    i_ignited: Any = None

    pinnerzoneradmw: Any = None

    pradmw: Any = None

    itart: Any = None

    pdivt: Any = None

    p_fw_alpha_mw: Any = None

    idivrt: Any = None

    p_plasma_ohmic_mw: Any = None

    iradloss: Any = None

    fusion_power: Any = None

    non_alpha_charged_power: Any = None

    pscalingmw: Any = None

    f_alpha_plasma: Any = None

    p_cp_pump_cool: Any = None

    i_tf_sup: Any = None

    tfcmw: Any = None

    temp_tf_coil_cryo: Any = None

    tcoolin: Any = None

    eff_tf_cryo: Any = None

    p_fw_blkt_pump_cool_mw: Any = None

    p_shield_pump_cool_elec_mw: Any = None

    p_div_pump_cool_elec_mw: Any = None

    p_pump_coolant_total_mw: Any = None

    p_fw_blkt_coolant_thermal_mw: Any = None

    p_fw_blkt_pump_cool_elec_mw: Any = None

    p_div_coolant_thermal_mw: Any = None

    p_fw_coolant_thermal_mw: Any = None

    p_shield_coolant_thermal_mw: Any = None

    ppumpmw: Any = None

    p_core_electrical_mw: Any = None

    f_div_thermal_primary: Any = None

    qss: Any = None

    qac: Any = None

    qcl: Any = None

    qmisc: Any = None

    outfile: Any = None

    iprint: Any = None

    expected_p_net_electrical_mw: Any = None

    expected_p_recirc_electrical_mw: Any = None

    expected_fachtmw: Any = None

    expected_p_gross_electrical: Any = None

    expected_p_thermal_secondary_mw: Any = None

    expected_p_core_electrical_mw: Any = None


@pytest.mark.parametrize(
    "power2param",
    (
        Power2Param(
            pnetelin=500,
            ipnet=0,
            ireactor=1,
            pinjmw=51.978447720428512,
            emultmw=377.93233088402548,
            inuclear=1,
            p_blanket_nuclear_heat_mw=1504.711566619962,
            p_fw_radiation_mw=254.87601794907812,
            qnuc=12920,
            eta_pump_coolant_electrical=0.87000000000000011,
            emult=1.2690000534057617,
            p_div_radiation_mw=33.119482558354782,
            fdiv=0.115,
            fhcd=0,
            secondary_cycle=2,
            pnuc_cp=0,
            p_div_nuclear_heat_mw=182.69222981118057,
            primary_pumping=3,
            p_tf_nuclear_heat_mw=0.044178296011112193,
            p_hcd_nuclear_heat_mw=0,
            p_shield_nuclear_heat_mw=1.3609360176065353,
            p_hcd_radiation_mw=0,
            p_fw_nuclear_heat_mw=276.76827393356979,
            p_shield_pump_cool_mw=0.0068046800880326762,
            p_blanket_pumping_mw=0,
            p_shield_thermal_secondary_mw=0,
            fpumpshld=0.0050000000000000001,
            tturb=0,
            p_net_electrical_mw=0,
            fpumpdiv=0.0050000000000000001,
            fpumpblkt=0.0050000000000000001,
            p_vacuum_pumps_mw=0.5,
            p_div_pump_cool_mw=1.7942175899286208,
            nphx=3,
            helpow=76851.741036987034,
            p_fw_pumping_mw=0,
            p_recirc_electrical_mw=0,
            p_thermal_primary_mw=2620.2218111502593,
            fpumpfw=0.0050000000000000001,
            p_baseload_electrical_total_mw=61.882833632875375,
            i_shield_power_generation=1,
            p_hcd_electrical_mw=129.94611930107126,
            fachtmw=0,
            p_gross_electrical=0,
            p_thermal_secondary_mw=0,
            trithtmw=15,
            p_hcd_thermal_secondary_mw=0,
            p_tf_electrical_mw=9.1507079104675704,
            p_pump_cool_elec_total_mw=234.28554165620102,
            eta_thermal_electric=0.37500000000000006,
            p_cryo_plant_mw=37.900388528497025,
            p_div_thermal_secondary_mw=0,
            p_hcd_electrical_loss_mw=77.967671580642758,
            p_pump_cool_loss_total_mw=30.457120415306122,
            helpow_cryal=0,
            pfwpmw=0.89998039031509891,
            alpha_power_total=396.66154806848488,
            i_ignited=0,
            pinnerzoneradmw=113.53817859231452,
            pradmw=287.99550050743289,
            itart=0,
            pdivt=143.03180561618876,
            p_fw_alpha_mw=19.833077403424262,
            idivrt=1,
            p_plasma_ohmic_mw=0.61391840981850698,
            iradloss=1,
            fusion_power=1985.785106643267,
            non_alpha_charged_power=1.6064693283140403,
            pscalingmw=325.08626176539281,
            f_alpha_plasma=0.94999999999999996,
            p_cp_pump_cool=0,
            i_tf_sup=1,
            tfcmw=0,
            temp_tf_coil_cryo=4.5,
            tcoolin=313.14999999999998,
            eff_tf_cryo=0.13,
            p_fw_blkt_pump_cool_mw=202.02739897087824,
            p_shield_pump_cool_elec_mw=0.0078214713655548,
            p_div_pump_cool_elec_mw=2.0623190688834718,
            p_pump_coolant_total_mw=203.8284212408949,
            p_fw_blkt_coolant_thermal_mw=2258.2163348769122,
            p_fw_blkt_pump_cool_elec_mw=232.21540111595198,
            p_div_coolant_thermal_mw=360.63773557565275,
            p_fw_coolant_thermal_mw=0,
            p_shield_coolant_thermal_mw=1.3677406976945679,
            ppumpmw=0,
            p_core_electrical_mw=0,
            f_div_thermal_primary=0.13763633828287813,
            qss=20361.633927097802,
            qac=3611.3456752656607,
            qcl=16108.2211128,
            qmisc=23850.540321823562,
            outfile=11,
            iprint=0,
            expected_p_net_electrical_mw=493.01760776192009,
            expected_p_recirc_electrical_mw=489.56557141942733,
            expected_fachtmw=61.882833632875375,
            expected_p_gross_electrical=982.58317918134742,
            expected_p_thermal_secondary_mw=233.80288075411508,
            expected_p_core_electrical_mw=125.33391046215507,
        ),
        Power2Param(
            pnetelin=500,
            ipnet=0,
            ireactor=1,
            pinjmw=51.978447720428512,
            emultmw=377.8143718115644,
            inuclear=1,
            p_blanket_nuclear_heat_mw=1549.9285082739402,
            p_fw_radiation_mw=254.87601794907812,
            qnuc=12920,
            eta_pump_coolant_electrical=0.87000000000000011,
            emult=1.2690000534057617,
            p_div_radiation_mw=33.119482558354782,
            fdiv=0.115,
            fhcd=0,
            secondary_cycle=2,
            pnuc_cp=0,
            p_div_nuclear_heat_mw=182.6352084763719,
            primary_pumping=3,
            p_tf_nuclear_heat_mw=0.045535131445547841,
            p_hcd_nuclear_heat_mw=0,
            p_shield_nuclear_heat_mw=1.4036212304705389,
            p_hcd_radiation_mw=0,
            p_fw_nuclear_heat_mw=230.95082168283884,
            p_shield_pump_cool_mw=0.0070181061523526943,
            p_blanket_pumping_mw=0,
            p_shield_thermal_secondary_mw=0,
            fpumpshld=0.0050000000000000001,
            tturb=0,
            p_net_electrical_mw=493.01760776192009,
            fpumpdiv=0.0050000000000000001,
            fpumpblkt=0.0050000000000000001,
            p_vacuum_pumps_mw=0.5,
            p_div_pump_cool_mw=1.7933419035282543,
            nphx=3,
            helpow=220505.71684249729,
            p_fw_pumping_mw=0,
            p_recirc_electrical_mw=489.56557141942733,
            p_thermal_primary_mw=2619.4223856129224,
            fpumpfw=0.0050000000000000001,
            p_baseload_electrical_total_mw=62.237143915360818,
            i_shield_power_generation=1,
            p_hcd_electrical_mw=129.94611930107126,
            fachtmw=61.882833632875375,
            p_gross_electrical=982.58317918134742,
            p_thermal_secondary_mw=233.80288075411508,
            trithtmw=15,
            p_hcd_thermal_secondary_mw=0,
            p_tf_electrical_mw=9.1507079104675704,
            p_pump_cool_elec_total_mw=234.2162627659944,
            eta_thermal_electric=0.37500000000000006,
            p_cryo_plant_mw=108.74512702403499,
            p_div_thermal_secondary_mw=0,
            p_hcd_electrical_loss_mw=77.967671580642758,
            p_pump_cool_loss_total_mw=30.448114159579291,
            helpow_cryal=0,
            pfwpmw=0.068213156646500808,
            alpha_power_total=396.53774329057228,
            i_ignited=0,
            pinnerzoneradmw=113.53817859231452,
            pradmw=287.99550050743289,
            itart=0,
            pdivt=142.91368967092416,
            p_fw_alpha_mw=19.826887164528632,
            idivrt=1,
            p_plasma_ohmic_mw=0.61391840981850698,
            iradloss=1,
            fusion_power=1985.1653095257811,
            non_alpha_charged_power=1.6059679220663614,
            pscalingmw=325.00280675287695,
            f_alpha_plasma=0.94999999999999996,
            p_cp_pump_cool=0,
            i_tf_sup=1,
            tfcmw=0,
            temp_tf_coil_cryo=4.5,
            tcoolin=313.14999999999998,
            eff_tf_cryo=0.13,
            p_fw_blkt_pump_cool_mw=201.96778859673452,
            p_shield_pump_cool_elec_mw=0.0080667886808651647,
            p_div_pump_cool_elec_mw=2.0613125327910966,
            p_pump_coolant_total_mw=203.76814860641511,
            p_fw_blkt_coolant_thermal_mw=2257.5500236671205,
            p_fw_blkt_pump_cool_elec_mw=232.14688344452242,
            p_div_coolant_thermal_mw=360.46172260917911,
            p_fw_coolant_thermal_mw=0,
            p_shield_coolant_thermal_mw=1.4106393366228915,
            ppumpmw=0,
            p_core_electrical_mw=125.33391046215507,
            f_div_thermal_primary=0.13761114839248584,
            qss=20342.863776957758,
            qac=102701.82327748176,
            qcl=16108.2211128,
            qmisc=68432.80867525778,
            outfile=11,
            iprint=0,
            expected_p_net_electrical_mw=422.4198205312706,
            expected_p_recirc_electrical_mw=559.86357407357548,
            expected_fachtmw=62.237143915360818,
            expected_p_gross_electrical=982.28339460484608,
            expected_p_thermal_secondary_mw=304.16251287817744,
            expected_p_core_electrical_mw=195.70119200650984,
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

    monkeypatch.setattr(constraint_variables, "pnetelin", power2param.pnetelin)

    monkeypatch.setattr(cost_variables, "ipnet", power2param.ipnet)

    monkeypatch.setattr(cost_variables, "ireactor", power2param.ireactor)

    monkeypatch.setattr(current_drive_variables, "pinjmw", power2param.pinjmw)

    monkeypatch.setattr(fwbs_variables, "emultmw", power2param.emultmw)

    monkeypatch.setattr(fwbs_variables, "inuclear", power2param.inuclear)

    monkeypatch.setattr(
        fwbs_variables,
        "p_blanket_nuclear_heat_mw",
        power2param.p_blanket_nuclear_heat_mw,
    )

    monkeypatch.setattr(
        fwbs_variables, "p_fw_radiation_mw", power2param.p_fw_radiation_mw
    )

    monkeypatch.setattr(fwbs_variables, "qnuc", power2param.qnuc)

    monkeypatch.setattr(
        fwbs_variables,
        "eta_pump_coolant_electrical",
        power2param.eta_pump_coolant_electrical,
    )

    monkeypatch.setattr(fwbs_variables, "emult", power2param.emult)

    monkeypatch.setattr(
        fwbs_variables, "p_div_radiation_mw", power2param.p_div_radiation_mw
    )

    monkeypatch.setattr(fwbs_variables, "fdiv", power2param.fdiv)

    monkeypatch.setattr(fwbs_variables, "fhcd", power2param.fhcd)

    monkeypatch.setattr(fwbs_variables, "secondary_cycle", power2param.secondary_cycle)

    monkeypatch.setattr(fwbs_variables, "pnuc_cp", power2param.pnuc_cp)

    monkeypatch.setattr(
        fwbs_variables, "p_div_nuclear_heat_mw", power2param.p_div_nuclear_heat_mw
    )

    monkeypatch.setattr(fwbs_variables, "primary_pumping", power2param.primary_pumping)

    monkeypatch.setattr(
        fwbs_variables, "p_tf_nuclear_heat_mw", power2param.p_tf_nuclear_heat_mw
    )

    monkeypatch.setattr(
        fwbs_variables, "p_hcd_nuclear_heat_mw", power2param.p_hcd_nuclear_heat_mw
    )

    monkeypatch.setattr(
        fwbs_variables, "p_shield_nuclear_heat_mw", power2param.p_shield_nuclear_heat_mw
    )

    monkeypatch.setattr(
        fwbs_variables, "p_hcd_radiation_mw", power2param.p_hcd_radiation_mw
    )

    monkeypatch.setattr(
        fwbs_variables, "p_fw_nuclear_heat_mw", power2param.p_fw_nuclear_heat_mw
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_shield_pump_cool_mw",
        power2param.p_shield_pump_cool_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_blanket_pumping_mw",
        power2param.p_blanket_pumping_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_shield_thermal_secondary_mw",
        power2param.p_shield_thermal_secondary_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "fpumpshld", power2param.fpumpshld)

    monkeypatch.setattr(heat_transport_variables, "tturb", power2param.tturb)

    monkeypatch.setattr(
        heat_transport_variables, "p_net_electrical_mw", power2param.p_net_electrical_mw
    )

    monkeypatch.setattr(heat_transport_variables, "fpumpdiv", power2param.fpumpdiv)

    monkeypatch.setattr(heat_transport_variables, "fpumpblkt", power2param.fpumpblkt)

    monkeypatch.setattr(
        heat_transport_variables, "p_vacuum_pumps_mw", power2param.p_vacuum_pumps_mw
    )

    monkeypatch.setattr(
        heat_transport_variables, "p_div_pump_cool_mw", power2param.p_div_pump_cool_mw
    )

    monkeypatch.setattr(heat_transport_variables, "nphx", power2param.nphx)

    monkeypatch.setattr(heat_transport_variables, "helpow", power2param.helpow)

    monkeypatch.setattr(
        heat_transport_variables, "p_fw_pumping_mw", power2param.p_fw_pumping_mw
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_recirc_electrical_mw",
        power2param.p_recirc_electrical_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_thermal_primary_mw",
        power2param.p_thermal_primary_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "fpumpfw", power2param.fpumpfw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_baseload_electrical_total_mw",
        power2param.p_baseload_electrical_total_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "i_shield_power_generation",
        power2param.i_shield_power_generation,
    )

    monkeypatch.setattr(
        heat_transport_variables, "p_hcd_electrical_mw", power2param.p_hcd_electrical_mw
    )

    monkeypatch.setattr(heat_transport_variables, "fachtmw", power2param.fachtmw)

    monkeypatch.setattr(
        heat_transport_variables, "p_gross_electrical", power2param.p_gross_electrical
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_thermal_secondary_mw",
        power2param.p_thermal_secondary_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "trithtmw", power2param.trithtmw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_thermal_secondary_mw",
        power2param.p_hcd_thermal_secondary_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables, "p_tf_electrical_mw", power2param.p_tf_electrical_mw
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_pump_cool_elec_total_mw",
        power2param.p_pump_cool_elec_total_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "eta_thermal_electric",
        power2param.eta_thermal_electric,
    )

    monkeypatch.setattr(
        heat_transport_variables, "p_cryo_plant_mw", power2param.p_cryo_plant_mw
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_div_thermal_secondary_mw",
        power2param.p_div_thermal_secondary_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electrical_loss_mw",
        power2param.p_hcd_electrical_loss_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_pump_cool_loss_total_mw",
        power2param.p_pump_cool_loss_total_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables, "helpow_cryal", power2param.helpow_cryal
    )

    monkeypatch.setattr(pfcoil_variables, "pfwpmw", power2param.pfwpmw)

    monkeypatch.setattr(
        physics_variables, "alpha_power_total", power2param.alpha_power_total
    )

    monkeypatch.setattr(physics_variables, "i_ignited", power2param.i_ignited)

    monkeypatch.setattr(
        physics_variables, "pinnerzoneradmw", power2param.pinnerzoneradmw
    )

    monkeypatch.setattr(physics_variables, "pradmw", power2param.pradmw)

    monkeypatch.setattr(physics_variables, "itart", power2param.itart)

    monkeypatch.setattr(physics_variables, "pdivt", power2param.pdivt)

    monkeypatch.setattr(physics_variables, "p_fw_alpha_mw", power2param.p_fw_alpha_mw)

    monkeypatch.setattr(physics_variables, "idivrt", power2param.idivrt)

    monkeypatch.setattr(
        physics_variables, "p_plasma_ohmic_mw", power2param.p_plasma_ohmic_mw
    )

    monkeypatch.setattr(physics_variables, "iradloss", power2param.iradloss)

    monkeypatch.setattr(physics_variables, "fusion_power", power2param.fusion_power)

    monkeypatch.setattr(
        physics_variables,
        "non_alpha_charged_power",
        power2param.non_alpha_charged_power,
    )

    monkeypatch.setattr(physics_variables, "pscalingmw", power2param.pscalingmw)

    monkeypatch.setattr(physics_variables, "f_alpha_plasma", power2param.f_alpha_plasma)

    monkeypatch.setattr(tfcoil_variables, "p_cp_pump_cool", power2param.p_cp_pump_cool)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", power2param.i_tf_sup)

    monkeypatch.setattr(tfcoil_variables, "tfcmw", power2param.tfcmw)

    monkeypatch.setattr(
        tfcoil_variables, "temp_tf_coil_cryo", power2param.temp_tf_coil_cryo
    )

    monkeypatch.setattr(tfcoil_variables, "tcoolin", power2param.tcoolin)

    monkeypatch.setattr(tfcoil_variables, "eff_tf_cryo", power2param.eff_tf_cryo)

    monkeypatch.setattr(
        ppv, "p_fw_blkt_pump_cool_mw", power2param.p_fw_blkt_pump_cool_mw
    )

    monkeypatch.setattr(
        power, "p_shield_pump_cool_elec_mw", power2param.p_shield_pump_cool_elec_mw
    )

    monkeypatch.setattr(
        power, "p_div_pump_cool_elec_mw", power2param.p_div_pump_cool_elec_mw
    )

    monkeypatch.setattr(
        power, "p_pump_coolant_total_mw", power2param.p_pump_coolant_total_mw
    )

    monkeypatch.setattr(
        power, "p_fw_blkt_coolant_thermal_mw", power2param.p_fw_blkt_coolant_thermal_mw
    )

    monkeypatch.setattr(
        power, "p_fw_blkt_pump_cool_elec_mw", power2param.p_fw_blkt_pump_cool_elec_mw
    )

    monkeypatch.setattr(
        power, "p_div_coolant_thermal_mw", power2param.p_div_coolant_thermal_mw
    )

    monkeypatch.setattr(
        power, "p_fw_coolant_thermal_mw", power2param.p_fw_coolant_thermal_mw
    )

    monkeypatch.setattr(
        power, "p_shield_coolant_thermal_mw", power2param.p_shield_coolant_thermal_mw
    )

    monkeypatch.setattr(power, "ppumpmw", power2param.ppumpmw)

    monkeypatch.setattr(power, "p_core_electrical_mw", power2param.p_core_electrical_mw)

    monkeypatch.setattr(
        power, "f_div_thermal_primary", power2param.f_div_thermal_primary
    )

    monkeypatch.setattr(power, "qss", power2param.qss)

    monkeypatch.setattr(power, "qac", power2param.qac)

    monkeypatch.setattr(power, "qcl", power2param.qcl)

    monkeypatch.setattr(power, "qmisc", power2param.qmisc)

    power.power2(output=False)

    assert heat_transport_variables.p_net_electrical_mw == pytest.approx(
        power2param.expected_p_net_electrical_mw
    )

    assert heat_transport_variables.p_recirc_electrical_mw == pytest.approx(
        power2param.expected_p_recirc_electrical_mw
    )

    assert heat_transport_variables.fachtmw == pytest.approx(
        power2param.expected_fachtmw
    )

    assert heat_transport_variables.p_gross_electrical == pytest.approx(
        power2param.expected_p_gross_electrical
    )

    assert heat_transport_variables.p_thermal_secondary_mw == pytest.approx(
        power2param.expected_p_thermal_secondary_mw
    )

    assert power.p_core_electrical_mw == pytest.approx(
        power2param.expected_p_core_electrical_mw
    )


class Power3Param(NamedTuple):
    etacd: Any = None

    p_pump_cool_elec_total_mw: Any = None

    pinjmax: Any = None

    p_cryo_plant_mw: Any = None

    p_vacuum_pumps_mw: Any = None

    p_tf_electrical_mw: Any = None

    trithtmw: Any = None

    p_hcd_electrical_mw: Any = None

    fachtmw: Any = None

    p_gross_electrical: Any = None

    poloidalpower: Any = None

    t_precharge: Any = None

    t_burn: Any = None

    t_fusion_ramp: Any = None

    t_between_pulse: Any = None

    t_ramp_down: Any = None

    t_current_ramp_up: Any = None

    outfile: Any = None

    iprint: Any = None


@pytest.mark.parametrize(
    "power3param",
    (
        Power3Param(
            etacd=0.40000000000000002,
            p_pump_cool_elec_total_mw=234.28554165620102,
            pinjmax=120,
            p_cryo_plant_mw=37.900388528497025,
            p_vacuum_pumps_mw=0.5,
            p_tf_electrical_mw=9.1507079104675704,
            trithtmw=15,
            p_hcd_electrical_mw=129.94611930107126,
            fachtmw=61.882833632875375,
            p_gross_electrical=982.58317918134742,
            poloidalpower=numpy.array(
                numpy.array(
                    (59332953.082890816, 43806300.444207191, 0, 0, -211211992.31967318),
                    order="F",
                ),
                order="F",
            ).transpose(),
            t_precharge=500,
            t_burn=0,
            t_fusion_ramp=10,
            t_between_pulse=0,
            t_ramp_down=177.21306969367816,
            t_current_ramp_up=177.21306969367816,
            outfile=11,
            iprint=0,
        ),
        Power3Param(
            etacd=0.40000000000000002,
            p_pump_cool_elec_total_mw=234.2162627659944,
            pinjmax=120,
            p_cryo_plant_mw=108.74512702403499,
            p_vacuum_pumps_mw=0.5,
            p_tf_electrical_mw=9.1507079104675704,
            trithtmw=15,
            p_hcd_electrical_mw=129.94611930107126,
            fachtmw=62.237143915360818,
            p_gross_electrical=982.28339460484608,
            poloidalpower=numpy.array(
                numpy.array(
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
            t_precharge=500,
            t_burn=10230.533336387543,
            t_fusion_ramp=10,
            t_between_pulse=0,
            t_ramp_down=177.21306969367816,
            t_current_ramp_up=177.21306969367816,
            outfile=11,
            iprint=0,
        ),
    ),
)
def test_power3(power3param, monkeypatch, power):
    """
    Automatically generated Regression Unit Test for power3.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param power3param: the data used to mock and assert in this test.
    :type power3param: power3param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(current_drive_variables, "etacd", power3param.etacd)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_pump_cool_elec_total_mw",
        power3param.p_pump_cool_elec_total_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "pinjmax", power3param.pinjmax)

    monkeypatch.setattr(
        heat_transport_variables, "p_cryo_plant_mw", power3param.p_cryo_plant_mw
    )

    monkeypatch.setattr(
        heat_transport_variables, "p_vacuum_pumps_mw", power3param.p_vacuum_pumps_mw
    )

    monkeypatch.setattr(
        heat_transport_variables, "p_tf_electrical_mw", power3param.p_tf_electrical_mw
    )

    monkeypatch.setattr(heat_transport_variables, "trithtmw", power3param.trithtmw)

    monkeypatch.setattr(
        heat_transport_variables, "p_hcd_electrical_mw", power3param.p_hcd_electrical_mw
    )

    monkeypatch.setattr(heat_transport_variables, "fachtmw", power3param.fachtmw)

    monkeypatch.setattr(
        heat_transport_variables, "p_gross_electrical", power3param.p_gross_electrical
    )

    monkeypatch.setattr(pf_power_variables, "poloidalpower", power3param.poloidalpower)

    monkeypatch.setattr(times_variables, "t_precharge", power3param.t_precharge)

    monkeypatch.setattr(times_variables, "t_burn", power3param.t_burn)

    monkeypatch.setattr(times_variables, "t_fusion_ramp", power3param.t_fusion_ramp)

    monkeypatch.setattr(times_variables, "t_between_pulse", power3param.t_between_pulse)

    monkeypatch.setattr(times_variables, "t_ramp_down", power3param.t_ramp_down)

    monkeypatch.setattr(
        times_variables, "t_current_ramp_up", power3param.t_current_ramp_up
    )

    power.power3(output=False)
