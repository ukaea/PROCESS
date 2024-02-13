from typing import NamedTuple, Any
import pytest
import numpy


from process.fortran import fwbs_variables
from process.fortran import heat_transport_variables
from process.fortran import pfcoil_variables
from process.fortran import numerics
from process.fortran import physics_variables
from process.fortran import build_variables
from process.fortran import pf_power_variables
from process.fortran import times_variables
from process.fortran import buildings_variables
from process.fortran import primary_pumping_variables as ppv
from process.fortran import constraint_variables
from process.fortran import cost_variables
from process.fortran import current_drive_variables
from process.fortran import tfcoil_variables
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

    cpttf: Any = None

    ensxpfm: Any = None

    ptfnuc: Any = None

    n_tf: Any = None

    tfsai: Any = None

    tpulse: Any = None

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
            cpttf=74026.751437500003,
            ensxpfm=37429.525515086898,
            ptfnuc=0.044178296011112193,
            n_tf=16,
            tfsai=0,
            tpulse=10364.426139387357,
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
            cpttf=74026.751437500003,
            ensxpfm=37427.228965055205,
            ptfnuc=0.045535131445547841,
            n_tf=16,
            tfsai=0,
            tpulse=364.42613938735633,
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
        coldmass=cryoparam.coldmass,
        cpttf=cryoparam.cpttf,
        ensxpfm=cryoparam.ensxpfm,
        ptfnuc=cryoparam.ptfnuc,
        n_tf=cryoparam.n_tf,
        tfsai=cryoparam.tfsai,
        tpulse=cryoparam.tpulse,
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

    srcktpm: Any = None

    ngrp: Any = None

    cpt: Any = None

    pfwpmw: Any = None

    pfclres: Any = None

    ncirt: Any = None

    ncls: Any = None

    ric: Any = None

    etapsu: Any = None

    cptdin: Any = None

    curpfb: Any = None

    sxlg: Any = None

    turns: Any = None

    vf: Any = None

    rjconpf: Any = None

    rpf: Any = None

    pohmmw: Any = None

    rmajor: Any = None

    active_constraints: Any = None

    ioptimz: Any = None

    tim: Any = None

    intervallabel: Any = None

    timelabel: Any = None

    tohs: Any = None

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
            poloidalpower=numpy.array(
                numpy.array((0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            spsmva=0,
            vpfskv=0,
            ensxpfm=0,
            acptmax=0,
            srcktpm=0,
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
            curpfb=numpy.array(
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
            pohmmw=0.61391840981850698,
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
                "tramp      ",
                "tohs       ",
                "t_fusion_ramp      ",
                "tburn      ",
                "tqnch      ",
            ),
            timelabel=(
                "Start      ",
                "BOP        ",
                "EOR        ",
                "BOF        ",
                "EOF        ",
                "EOP        ",
            ),
            tohs=177.21306969367816,
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
            expected_srcktpm=1071.1112934857531,
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
            srcktpm=1071.1112934857531,
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
            curpfb=numpy.array(
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
            pohmmw=0.61391840981850698,
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
                "tramp      ",
                "tohs       ",
                "t_fusion_ramp      ",
                "tburn      ",
                "tqnch      ",
            ),
            timelabel=(
                "Start      ",
                "BOP        ",
                "EOR        ",
                "BOF        ",
                "EOF        ",
                "EOP        ",
            ),
            tohs=177.21306969367816,
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

    monkeypatch.setattr(pfcoil_variables, "ngrp", pfpwrparam.ngrp)

    monkeypatch.setattr(pfcoil_variables, "cpt", pfpwrparam.cpt)

    monkeypatch.setattr(pfcoil_variables, "pfwpmw", pfpwrparam.pfwpmw)

    monkeypatch.setattr(pfcoil_variables, "pfclres", pfpwrparam.pfclres)

    monkeypatch.setattr(pfcoil_variables, "ncirt", pfpwrparam.ncirt)

    monkeypatch.setattr(pfcoil_variables, "ncls", pfpwrparam.ncls)

    monkeypatch.setattr(pfcoil_variables, "ric", pfpwrparam.ric)

    monkeypatch.setattr(pfcoil_variables, "etapsu", pfpwrparam.etapsu)

    monkeypatch.setattr(pfcoil_variables, "cptdin", pfpwrparam.cptdin)

    monkeypatch.setattr(pfcoil_variables, "curpfb", pfpwrparam.curpfb)

    monkeypatch.setattr(pfcoil_variables, "sxlg", pfpwrparam.sxlg)

    monkeypatch.setattr(pfcoil_variables, "turns", pfpwrparam.turns)

    monkeypatch.setattr(pfcoil_variables, "vf", pfpwrparam.vf)

    monkeypatch.setattr(pfcoil_variables, "rjconpf", pfpwrparam.rjconpf)

    monkeypatch.setattr(pfcoil_variables, "rpf", pfpwrparam.rpf)

    monkeypatch.setattr(physics_variables, "pohmmw", pfpwrparam.pohmmw)

    monkeypatch.setattr(physics_variables, "rmajor", pfpwrparam.rmajor)

    monkeypatch.setattr(numerics, "active_constraints", pfpwrparam.active_constraints)

    monkeypatch.setattr(numerics, "ioptimz", pfpwrparam.ioptimz)

    monkeypatch.setattr(times_variables, "tim", pfpwrparam.tim)

    monkeypatch.setattr(times_variables, "tohs", pfpwrparam.tohs)

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

    assert pf_power_variables.srcktpm == pytest.approx(pfpwrparam.expected_srcktpm)


class AcpowParam(NamedTuple):

    efloor: Any = None

    baseel: Any = None

    crypmw: Any = None

    vachtmw: Any = None

    tfacpd: Any = None

    trithtmw: Any = None

    pinjwp: Any = None

    tlvpmw: Any = None

    peakmva: Any = None

    fcsht: Any = None

    fmgdmw: Any = None

    pwpm2: Any = None

    htpmw: Any = None

    pacpmw: Any = None

    iscenr: Any = None

    srcktpm: Any = None

    iprint: Any = None

    outfile: Any = None

    expected_pacpmw: Any = None


@pytest.mark.parametrize(
    "acpowparam",
    (
        AcpowParam(
            efloor=379218.8908858358,
            baseel=5000000,
            crypmw=37.900388528497025,
            vachtmw=0.5,
            tfacpd=9.1507079104675704,
            trithtmw=15,
            pinjwp=129.94611930107126,
            tlvpmw=0,
            peakmva=736.39062584245937,
            fcsht=0,
            fmgdmw=0,
            pwpm2=150,
            htpmw=234.28554165620102,
            pacpmw=0,
            iscenr=2,
            srcktpm=1071.1112934857531,
            iprint=0,
            outfile=11,
            expected_pacpmw=1226.1273281650574,
        ),
        AcpowParam(
            efloor=381580.9594357388,
            baseel=5000000,
            crypmw=108.74512702403499,
            vachtmw=0.5,
            tfacpd=9.1507079104675704,
            trithtmw=15,
            pinjwp=129.94611930107126,
            tlvpmw=699.34943812129745,
            peakmva=90.673341440806112,
            fcsht=61.882833632875375,
            fmgdmw=0,
            pwpm2=150,
            htpmw=234.2162627659944,
            pacpmw=1226.1273281650574,
            iscenr=2,
            srcktpm=1069.8879533693198,
            iprint=0,
            outfile=11,
            expected_pacpmw=651.53859031110437,
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

    monkeypatch.setattr(buildings_variables, "efloor", acpowparam.efloor)

    monkeypatch.setattr(heat_transport_variables, "baseel", acpowparam.baseel)

    monkeypatch.setattr(heat_transport_variables, "crypmw", acpowparam.crypmw)

    monkeypatch.setattr(heat_transport_variables, "vachtmw", acpowparam.vachtmw)

    monkeypatch.setattr(heat_transport_variables, "tfacpd", acpowparam.tfacpd)

    monkeypatch.setattr(heat_transport_variables, "trithtmw", acpowparam.trithtmw)

    monkeypatch.setattr(heat_transport_variables, "pinjwp", acpowparam.pinjwp)

    monkeypatch.setattr(heat_transport_variables, "tlvpmw", acpowparam.tlvpmw)

    monkeypatch.setattr(heat_transport_variables, "peakmva", acpowparam.peakmva)

    monkeypatch.setattr(heat_transport_variables, "fcsht", acpowparam.fcsht)

    monkeypatch.setattr(heat_transport_variables, "fmgdmw", acpowparam.fmgdmw)

    monkeypatch.setattr(heat_transport_variables, "pwpm2", acpowparam.pwpm2)

    monkeypatch.setattr(heat_transport_variables, "htpmw", acpowparam.htpmw)

    monkeypatch.setattr(heat_transport_variables, "pacpmw", acpowparam.pacpmw)

    monkeypatch.setattr(pf_power_variables, "iscenr", acpowparam.iscenr)

    monkeypatch.setattr(pf_power_variables, "srcktpm", acpowparam.srcktpm)

    power.acpow(output=False)

    assert heat_transport_variables.pacpmw == pytest.approx(acpowparam.expected_pacpmw)


class Power2Param(NamedTuple):

    pnetelin: Any = None

    ipnet: Any = None

    ireactor: Any = None

    pinjmw: Any = None

    emultmw: Any = None

    inuclear: Any = None

    pnucblkt: Any = None

    pradfw: Any = None

    qnuc: Any = None

    etahtp: Any = None

    emult: Any = None

    praddiv: Any = None

    fdiv: Any = None

    fhcd: Any = None

    secondary_cycle: Any = None

    pnuc_cp: Any = None

    pnucdiv: Any = None

    primary_pumping: Any = None

    ptfnuc: Any = None

    pnuchcd: Any = None

    pnucshld: Any = None

    pradhcd: Any = None

    pnucfw: Any = None

    htpmw_shld: Any = None

    htpmw_blkt: Any = None

    psecshld: Any = None

    fpumpshld: Any = None

    tturb: Any = None

    pnetelmw: Any = None

    fpumpdiv: Any = None

    fpumpblkt: Any = None

    vachtmw: Any = None

    htpmw_div: Any = None

    nphx: Any = None

    helpow: Any = None

    htpmw_fw: Any = None

    precircmw: Any = None

    pthermmw: Any = None

    fpumpfw: Any = None

    fcsht: Any = None

    iprimshld: Any = None

    pinjwp: Any = None

    fachtmw: Any = None

    pgrossmw: Any = None

    psechtmw: Any = None

    trithtmw: Any = None

    psechcd: Any = None

    tfacpd: Any = None

    htpmw: Any = None

    etath: Any = None

    crypmw: Any = None

    psecdiv: Any = None

    pinjht: Any = None

    htpsecmw: Any = None

    helpow_cryal: Any = None

    pfwpmw: Any = None

    palpmw: Any = None

    ignite: Any = None

    pinnerzoneradmw: Any = None

    pradmw: Any = None

    itart: Any = None

    pdivt: Any = None

    palpfwmw: Any = None

    idivrt: Any = None

    pohmmw: Any = None

    iradloss: Any = None

    powfmw: Any = None

    pchargemw: Any = None

    pscalingmw: Any = None

    falpha: Any = None

    ppump: Any = None

    i_tf_sup: Any = None

    tfcmw: Any = None

    tmpcry: Any = None

    tcoolin: Any = None

    eff_tf_cryo: Any = None

    htpmw_fw_blkt: Any = None

    htpmwe_shld: Any = None

    htpmwe_div: Any = None

    htpmw_mech: Any = None

    pthermfw_blkt: Any = None

    htpmwe_fw_blkt: Any = None

    pthermdiv: Any = None

    pthermfw: Any = None

    pthermshld: Any = None

    ppumpmw: Any = None

    pcoresystems: Any = None

    pdivfraction: Any = None

    qss: Any = None

    qac: Any = None

    qcl: Any = None

    qmisc: Any = None

    outfile: Any = None

    iprint: Any = None

    expected_pnetelmw: Any = None

    expected_precircmw: Any = None

    expected_fachtmw: Any = None

    expected_pgrossmw: Any = None

    expected_psechtmw: Any = None

    expected_pcoresystems: Any = None


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
            pnucblkt=1504.711566619962,
            pradfw=254.87601794907812,
            qnuc=12920,
            etahtp=0.87000000000000011,
            emult=1.2690000534057617,
            praddiv=33.119482558354782,
            fdiv=0.115,
            fhcd=0,
            secondary_cycle=2,
            pnuc_cp=0,
            pnucdiv=182.69222981118057,
            primary_pumping=3,
            ptfnuc=0.044178296011112193,
            pnuchcd=0,
            pnucshld=1.3609360176065353,
            pradhcd=0,
            pnucfw=276.76827393356979,
            htpmw_shld=0.0068046800880326762,
            htpmw_blkt=0,
            psecshld=0,
            fpumpshld=0.0050000000000000001,
            tturb=0,
            pnetelmw=0,
            fpumpdiv=0.0050000000000000001,
            fpumpblkt=0.0050000000000000001,
            vachtmw=0.5,
            htpmw_div=1.7942175899286208,
            nphx=3,
            helpow=76851.741036987034,
            htpmw_fw=0,
            precircmw=0,
            pthermmw=2620.2218111502593,
            fpumpfw=0.0050000000000000001,
            fcsht=61.882833632875375,
            iprimshld=1,
            pinjwp=129.94611930107126,
            fachtmw=0,
            pgrossmw=0,
            psechtmw=0,
            trithtmw=15,
            psechcd=0,
            tfacpd=9.1507079104675704,
            htpmw=234.28554165620102,
            etath=0.37500000000000006,
            crypmw=37.900388528497025,
            psecdiv=0,
            pinjht=77.967671580642758,
            htpsecmw=30.457120415306122,
            helpow_cryal=0,
            pfwpmw=0.89998039031509891,
            palpmw=396.66154806848488,
            ignite=0,
            pinnerzoneradmw=113.53817859231452,
            pradmw=287.99550050743289,
            itart=0,
            pdivt=143.03180561618876,
            palpfwmw=19.833077403424262,
            idivrt=1,
            pohmmw=0.61391840981850698,
            iradloss=1,
            powfmw=1985.785106643267,
            pchargemw=1.6064693283140403,
            pscalingmw=325.08626176539281,
            falpha=0.94999999999999996,
            ppump=0,
            i_tf_sup=1,
            tfcmw=0,
            tmpcry=4.5,
            tcoolin=313.14999999999998,
            eff_tf_cryo=0.13,
            htpmw_fw_blkt=202.02739897087824,
            htpmwe_shld=0.0078214713655548,
            htpmwe_div=2.0623190688834718,
            htpmw_mech=203.8284212408949,
            pthermfw_blkt=2258.2163348769122,
            htpmwe_fw_blkt=232.21540111595198,
            pthermdiv=360.63773557565275,
            pthermfw=0,
            pthermshld=1.3677406976945679,
            ppumpmw=0,
            pcoresystems=0,
            pdivfraction=0.13763633828287813,
            qss=20361.633927097802,
            qac=3611.3456752656607,
            qcl=16108.2211128,
            qmisc=23850.540321823562,
            outfile=11,
            iprint=0,
            expected_pnetelmw=493.01760776192009,
            expected_precircmw=489.56557141942733,
            expected_fachtmw=61.882833632875375,
            expected_pgrossmw=982.58317918134742,
            expected_psechtmw=233.80288075411508,
            expected_pcoresystems=125.33391046215507,
        ),
        Power2Param(
            pnetelin=500,
            ipnet=0,
            ireactor=1,
            pinjmw=51.978447720428512,
            emultmw=377.8143718115644,
            inuclear=1,
            pnucblkt=1549.9285082739402,
            pradfw=254.87601794907812,
            qnuc=12920,
            etahtp=0.87000000000000011,
            emult=1.2690000534057617,
            praddiv=33.119482558354782,
            fdiv=0.115,
            fhcd=0,
            secondary_cycle=2,
            pnuc_cp=0,
            pnucdiv=182.6352084763719,
            primary_pumping=3,
            ptfnuc=0.045535131445547841,
            pnuchcd=0,
            pnucshld=1.4036212304705389,
            pradhcd=0,
            pnucfw=230.95082168283884,
            htpmw_shld=0.0070181061523526943,
            htpmw_blkt=0,
            psecshld=0,
            fpumpshld=0.0050000000000000001,
            tturb=0,
            pnetelmw=493.01760776192009,
            fpumpdiv=0.0050000000000000001,
            fpumpblkt=0.0050000000000000001,
            vachtmw=0.5,
            htpmw_div=1.7933419035282543,
            nphx=3,
            helpow=220505.71684249729,
            htpmw_fw=0,
            precircmw=489.56557141942733,
            pthermmw=2619.4223856129224,
            fpumpfw=0.0050000000000000001,
            fcsht=62.237143915360818,
            iprimshld=1,
            pinjwp=129.94611930107126,
            fachtmw=61.882833632875375,
            pgrossmw=982.58317918134742,
            psechtmw=233.80288075411508,
            trithtmw=15,
            psechcd=0,
            tfacpd=9.1507079104675704,
            htpmw=234.2162627659944,
            etath=0.37500000000000006,
            crypmw=108.74512702403499,
            psecdiv=0,
            pinjht=77.967671580642758,
            htpsecmw=30.448114159579291,
            helpow_cryal=0,
            pfwpmw=0.068213156646500808,
            palpmw=396.53774329057228,
            ignite=0,
            pinnerzoneradmw=113.53817859231452,
            pradmw=287.99550050743289,
            itart=0,
            pdivt=142.91368967092416,
            palpfwmw=19.826887164528632,
            idivrt=1,
            pohmmw=0.61391840981850698,
            iradloss=1,
            powfmw=1985.1653095257811,
            pchargemw=1.6059679220663614,
            pscalingmw=325.00280675287695,
            falpha=0.94999999999999996,
            ppump=0,
            i_tf_sup=1,
            tfcmw=0,
            tmpcry=4.5,
            tcoolin=313.14999999999998,
            eff_tf_cryo=0.13,
            htpmw_fw_blkt=201.96778859673452,
            htpmwe_shld=0.0080667886808651647,
            htpmwe_div=2.0613125327910966,
            htpmw_mech=203.76814860641511,
            pthermfw_blkt=2257.5500236671205,
            htpmwe_fw_blkt=232.14688344452242,
            pthermdiv=360.46172260917911,
            pthermfw=0,
            pthermshld=1.4106393366228915,
            ppumpmw=0,
            pcoresystems=125.33391046215507,
            pdivfraction=0.13761114839248584,
            qss=20342.863776957758,
            qac=102701.82327748176,
            qcl=16108.2211128,
            qmisc=68432.80867525778,
            outfile=11,
            iprint=0,
            expected_pnetelmw=422.4198205312706,
            expected_precircmw=559.86357407357548,
            expected_fachtmw=62.237143915360818,
            expected_pgrossmw=982.28339460484608,
            expected_psechtmw=304.16251287817744,
            expected_pcoresystems=195.70119200650984,
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

    monkeypatch.setattr(fwbs_variables, "pnucblkt", power2param.pnucblkt)

    monkeypatch.setattr(fwbs_variables, "pradfw", power2param.pradfw)

    monkeypatch.setattr(fwbs_variables, "qnuc", power2param.qnuc)

    monkeypatch.setattr(fwbs_variables, "etahtp", power2param.etahtp)

    monkeypatch.setattr(fwbs_variables, "emult", power2param.emult)

    monkeypatch.setattr(fwbs_variables, "praddiv", power2param.praddiv)

    monkeypatch.setattr(fwbs_variables, "fdiv", power2param.fdiv)

    monkeypatch.setattr(fwbs_variables, "fhcd", power2param.fhcd)

    monkeypatch.setattr(fwbs_variables, "secondary_cycle", power2param.secondary_cycle)

    monkeypatch.setattr(fwbs_variables, "pnuc_cp", power2param.pnuc_cp)

    monkeypatch.setattr(fwbs_variables, "pnucdiv", power2param.pnucdiv)

    monkeypatch.setattr(fwbs_variables, "primary_pumping", power2param.primary_pumping)

    monkeypatch.setattr(fwbs_variables, "ptfnuc", power2param.ptfnuc)

    monkeypatch.setattr(fwbs_variables, "pnuchcd", power2param.pnuchcd)

    monkeypatch.setattr(fwbs_variables, "pnucshld", power2param.pnucshld)

    monkeypatch.setattr(fwbs_variables, "pradhcd", power2param.pradhcd)

    monkeypatch.setattr(fwbs_variables, "pnucfw", power2param.pnucfw)

    monkeypatch.setattr(heat_transport_variables, "htpmw_shld", power2param.htpmw_shld)

    monkeypatch.setattr(heat_transport_variables, "htpmw_blkt", power2param.htpmw_blkt)

    monkeypatch.setattr(heat_transport_variables, "psecshld", power2param.psecshld)

    monkeypatch.setattr(heat_transport_variables, "fpumpshld", power2param.fpumpshld)

    monkeypatch.setattr(heat_transport_variables, "tturb", power2param.tturb)

    monkeypatch.setattr(heat_transport_variables, "pnetelmw", power2param.pnetelmw)

    monkeypatch.setattr(heat_transport_variables, "fpumpdiv", power2param.fpumpdiv)

    monkeypatch.setattr(heat_transport_variables, "fpumpblkt", power2param.fpumpblkt)

    monkeypatch.setattr(heat_transport_variables, "vachtmw", power2param.vachtmw)

    monkeypatch.setattr(heat_transport_variables, "htpmw_div", power2param.htpmw_div)

    monkeypatch.setattr(heat_transport_variables, "nphx", power2param.nphx)

    monkeypatch.setattr(heat_transport_variables, "helpow", power2param.helpow)

    monkeypatch.setattr(heat_transport_variables, "htpmw_fw", power2param.htpmw_fw)

    monkeypatch.setattr(heat_transport_variables, "precircmw", power2param.precircmw)

    monkeypatch.setattr(heat_transport_variables, "pthermmw", power2param.pthermmw)

    monkeypatch.setattr(heat_transport_variables, "fpumpfw", power2param.fpumpfw)

    monkeypatch.setattr(heat_transport_variables, "fcsht", power2param.fcsht)

    monkeypatch.setattr(heat_transport_variables, "iprimshld", power2param.iprimshld)

    monkeypatch.setattr(heat_transport_variables, "pinjwp", power2param.pinjwp)

    monkeypatch.setattr(heat_transport_variables, "fachtmw", power2param.fachtmw)

    monkeypatch.setattr(heat_transport_variables, "pgrossmw", power2param.pgrossmw)

    monkeypatch.setattr(heat_transport_variables, "psechtmw", power2param.psechtmw)

    monkeypatch.setattr(heat_transport_variables, "trithtmw", power2param.trithtmw)

    monkeypatch.setattr(heat_transport_variables, "psechcd", power2param.psechcd)

    monkeypatch.setattr(heat_transport_variables, "tfacpd", power2param.tfacpd)

    monkeypatch.setattr(heat_transport_variables, "htpmw", power2param.htpmw)

    monkeypatch.setattr(heat_transport_variables, "etath", power2param.etath)

    monkeypatch.setattr(heat_transport_variables, "crypmw", power2param.crypmw)

    monkeypatch.setattr(heat_transport_variables, "psecdiv", power2param.psecdiv)

    monkeypatch.setattr(heat_transport_variables, "pinjht", power2param.pinjht)

    monkeypatch.setattr(heat_transport_variables, "htpsecmw", power2param.htpsecmw)

    monkeypatch.setattr(
        heat_transport_variables, "helpow_cryal", power2param.helpow_cryal
    )

    monkeypatch.setattr(pfcoil_variables, "pfwpmw", power2param.pfwpmw)

    monkeypatch.setattr(physics_variables, "palpmw", power2param.palpmw)

    monkeypatch.setattr(physics_variables, "ignite", power2param.ignite)

    monkeypatch.setattr(
        physics_variables, "pinnerzoneradmw", power2param.pinnerzoneradmw
    )

    monkeypatch.setattr(physics_variables, "pradmw", power2param.pradmw)

    monkeypatch.setattr(physics_variables, "itart", power2param.itart)

    monkeypatch.setattr(physics_variables, "pdivt", power2param.pdivt)

    monkeypatch.setattr(physics_variables, "palpfwmw", power2param.palpfwmw)

    monkeypatch.setattr(physics_variables, "idivrt", power2param.idivrt)

    monkeypatch.setattr(physics_variables, "pohmmw", power2param.pohmmw)

    monkeypatch.setattr(physics_variables, "iradloss", power2param.iradloss)

    monkeypatch.setattr(physics_variables, "powfmw", power2param.powfmw)

    monkeypatch.setattr(physics_variables, "pchargemw", power2param.pchargemw)

    monkeypatch.setattr(physics_variables, "pscalingmw", power2param.pscalingmw)

    monkeypatch.setattr(physics_variables, "falpha", power2param.falpha)

    monkeypatch.setattr(tfcoil_variables, "ppump", power2param.ppump)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", power2param.i_tf_sup)

    monkeypatch.setattr(tfcoil_variables, "tfcmw", power2param.tfcmw)

    monkeypatch.setattr(tfcoil_variables, "tmpcry", power2param.tmpcry)

    monkeypatch.setattr(tfcoil_variables, "tcoolin", power2param.tcoolin)

    monkeypatch.setattr(tfcoil_variables, "eff_tf_cryo", power2param.eff_tf_cryo)

    monkeypatch.setattr(ppv, "htpmw_fw_blkt", power2param.htpmw_fw_blkt)

    monkeypatch.setattr(power, "htpmwe_shld", power2param.htpmwe_shld)

    monkeypatch.setattr(power, "htpmwe_div", power2param.htpmwe_div)

    monkeypatch.setattr(power, "htpmw_mech", power2param.htpmw_mech)

    monkeypatch.setattr(power, "pthermfw_blkt", power2param.pthermfw_blkt)

    monkeypatch.setattr(power, "htpmwe_fw_blkt", power2param.htpmwe_fw_blkt)

    monkeypatch.setattr(power, "pthermdiv", power2param.pthermdiv)

    monkeypatch.setattr(power, "pthermfw", power2param.pthermfw)

    monkeypatch.setattr(power, "pthermshld", power2param.pthermshld)

    monkeypatch.setattr(power, "ppumpmw", power2param.ppumpmw)

    monkeypatch.setattr(power, "pcoresystems", power2param.pcoresystems)

    monkeypatch.setattr(power, "pdivfraction", power2param.pdivfraction)

    monkeypatch.setattr(power, "qss", power2param.qss)

    monkeypatch.setattr(power, "qac", power2param.qac)

    monkeypatch.setattr(power, "qcl", power2param.qcl)

    monkeypatch.setattr(power, "qmisc", power2param.qmisc)

    power.power2(output=False)

    assert heat_transport_variables.pnetelmw == pytest.approx(
        power2param.expected_pnetelmw
    )

    assert heat_transport_variables.precircmw == pytest.approx(
        power2param.expected_precircmw
    )

    assert heat_transport_variables.fachtmw == pytest.approx(
        power2param.expected_fachtmw
    )

    assert heat_transport_variables.pgrossmw == pytest.approx(
        power2param.expected_pgrossmw
    )

    assert heat_transport_variables.psechtmw == pytest.approx(
        power2param.expected_psechtmw
    )

    assert power.pcoresystems == pytest.approx(power2param.expected_pcoresystems)


class Power3Param(NamedTuple):

    etacd: Any = None

    htpmw: Any = None

    pinjmax: Any = None

    crypmw: Any = None

    vachtmw: Any = None

    tfacpd: Any = None

    trithtmw: Any = None

    pinjwp: Any = None

    fachtmw: Any = None

    pgrossmw: Any = None

    poloidalpower: Any = None

    tramp: Any = None

    tburn: Any = None

    t_fusion_ramp: Any = None

    tdwell: Any = None

    tqnch: Any = None

    tohs: Any = None

    outfile: Any = None

    iprint: Any = None


@pytest.mark.parametrize(
    "power3param",
    (
        Power3Param(
            etacd=0.40000000000000002,
            htpmw=234.28554165620102,
            pinjmax=120,
            crypmw=37.900388528497025,
            vachtmw=0.5,
            tfacpd=9.1507079104675704,
            trithtmw=15,
            pinjwp=129.94611930107126,
            fachtmw=61.882833632875375,
            pgrossmw=982.58317918134742,
            poloidalpower=numpy.array(
                numpy.array(
                    (59332953.082890816, 43806300.444207191, 0, 0, -211211992.31967318),
                    order="F",
                ),
                order="F",
            ).transpose(),
            tramp=500,
            tburn=0,
            t_fusion_ramp=10,
            tdwell=0,
            tqnch=177.21306969367816,
            tohs=177.21306969367816,
            outfile=11,
            iprint=0,
        ),
        Power3Param(
            etacd=0.40000000000000002,
            htpmw=234.2162627659944,
            pinjmax=120,
            crypmw=108.74512702403499,
            vachtmw=0.5,
            tfacpd=9.1507079104675704,
            trithtmw=15,
            pinjwp=129.94611930107126,
            fachtmw=62.237143915360818,
            pgrossmw=982.28339460484608,
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
            tramp=500,
            tburn=10230.533336387543,
            t_fusion_ramp=10,
            tdwell=0,
            tqnch=177.21306969367816,
            tohs=177.21306969367816,
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

    monkeypatch.setattr(heat_transport_variables, "htpmw", power3param.htpmw)

    monkeypatch.setattr(heat_transport_variables, "pinjmax", power3param.pinjmax)

    monkeypatch.setattr(heat_transport_variables, "crypmw", power3param.crypmw)

    monkeypatch.setattr(heat_transport_variables, "vachtmw", power3param.vachtmw)

    monkeypatch.setattr(heat_transport_variables, "tfacpd", power3param.tfacpd)

    monkeypatch.setattr(heat_transport_variables, "trithtmw", power3param.trithtmw)

    monkeypatch.setattr(heat_transport_variables, "pinjwp", power3param.pinjwp)

    monkeypatch.setattr(heat_transport_variables, "fachtmw", power3param.fachtmw)

    monkeypatch.setattr(heat_transport_variables, "pgrossmw", power3param.pgrossmw)

    monkeypatch.setattr(pf_power_variables, "poloidalpower", power3param.poloidalpower)

    monkeypatch.setattr(times_variables, "tramp", power3param.tramp)

    monkeypatch.setattr(times_variables, "tburn", power3param.tburn)

    monkeypatch.setattr(times_variables, "t_fusion_ramp", power3param.t_fusion_ramp)

    monkeypatch.setattr(times_variables, "tdwell", power3param.tdwell)

    monkeypatch.setattr(times_variables, "tqnch", power3param.tqnch)

    monkeypatch.setattr(times_variables, "tohs", power3param.tohs)

    power.power3(output=False)
