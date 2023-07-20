"""Unit tests for physics_functions.f90."""
from typing import Any, NamedTuple
from process.fortran import physics_functions_module as pf
from process.fortran import physics_variables as pv
import pytest
from pytest import approx
import numpy


def test_plasma_elongation_ipb(monkeypatch):
    """Test plasma_elongation_IPB().
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(pv, "rmajor", 9.137)
    monkeypatch.setattr(pv, "rminor", 2.947)
    monkeypatch.setattr(pv, "vol", 2634.0)
    kappaa_ipb = pf.plasma_elongation_ipb()
    assert kappaa_ipb == approx(1.682, abs=0.001)


def test_total_mag_field(monkeypatch):
    """Test total_mag_field().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(pv, "bt", 5.278)
    monkeypatch.setattr(pv, "bp", 0.852)
    btot = pf.total_mag_field()
    assert btot == approx(5.347, abs=0.001)


def test_beta_poloidal(monkeypatch):
    """Test beta_poloidal().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(pv, "btot", 5.347)
    monkeypatch.setattr(pv, "beta", 0.0307)
    monkeypatch.setattr(pv, "bp", 0.852)
    betap = pf.beta_poloidal()
    assert betap == approx(1.209, abs=0.001)


def test_res_diff_time(monkeypatch):
    """Test res_diff_time().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(pv, "rmajor", 9.137)
    monkeypatch.setattr(pv, "rplas", 2.909e-9)
    monkeypatch.setattr(pv, "kappa95", 1.650)
    res_time = pf.res_diff_time()
    assert res_time == approx(4784.3, abs=0.1)


class Palph2Param(NamedTuple):
    falpha: Any = None

    fdeut: Any = None

    ifalphap: Any = None

    bp: Any = None

    bt: Any = None

    dene: Any = None

    deni: Any = None

    dnitot: Any = None

    falpe: Any = None

    falpi: Any = None

    palpnb: Any = None

    pchargepv: Any = None

    ten: Any = None

    tin: Any = None

    vol: Any = None

    palppv: Any = None

    pneutpv: Any = None

    expected_palppv: Any = None

    expected_pneutpv: Any = None

    expected_palpmw: Any = None

    expected_pneutmw: Any = None

    expected_pchargemw: Any = None

    expected_betaft: Any = None

    expected_palpepv: Any = None

    expected_palpipv: Any = None

    expected_pfuscmw: Any = None

    expected_powfmw: Any = None


@pytest.mark.parametrize(
    "palph2param",
    (
        Palph2Param(
            falpha=0.95,
            fdeut=0.5,
            ifalphap=1,
            bp=0.86,
            bt=5.3292,
            dene=7.432e19,
            deni=6.226e19,
            dnitot=6.743e19,
            falpe=0.68,
            falpi=0.32,
            palpnb=0,
            pchargepv=0.00066,
            ten=13.84,
            tin=13.84,
            vol=2426.25,
            palppv=0.163,
            pneutpv=0.654,
            expected_palppv=0.163,
            expected_pneutpv=0.654,
            expected_palpmw=395.47875,
            expected_pneutmw=1586.7675,
            expected_pchargemw=1.601325,
            expected_betaft=0.00423788,
            expected_palpipv=0.049552,
            expected_palpepv=0.105298,
            expected_pfuscmw=397.080075,
            expected_powfmw=1983.847575,
        ),
        Palph2Param(
            falpha=0.95,
            fdeut=0.5,
            ifalphap=1,
            bp=0.86,
            bt=5.3292,
            dene=7.432e19,
            deni=6.226e19,
            dnitot=6.743e19,
            falpe=0.68,
            falpi=0.32,
            palpnb=100.5,
            pchargepv=0.00066,
            ten=13.84,
            tin=13.84,
            vol=2426.25,
            palppv=0.163,
            pneutpv=0.654,
            expected_palppv=0.20442195,
            expected_pneutpv=0.81968779,
            expected_palpmw=495.97875,
            expected_pneutmw=1988.7675,
            expected_pchargemw=1.601325,
            expected_betaft=0.00531482,
            expected_palpipv=0.062144272,
            expected_palpepv=0.132056578,
            expected_pfuscmw=497.580075,
            expected_powfmw=2486.347575,
        ),
        Palph2Param(
            falpha=0.95,
            fdeut=0.5,
            ifalphap=0,
            bp=0.86,
            bt=5.3292,
            dene=7.432e19,
            deni=6.226e19,
            dnitot=6.743e19,
            falpe=0.68,
            falpi=0.32,
            palpnb=100.5,
            pchargepv=0.00066,
            ten=13.84,
            tin=13.84,
            vol=2426.25,
            palppv=0.163,
            pneutpv=0.654,
            expected_palppv=0.20442195,
            expected_pneutpv=0.81968779,
            expected_palpmw=495.97875,
            expected_pneutmw=1988.7675,
            expected_pchargemw=1.601325,
            expected_betaft=0.00701622,
            expected_palpipv=0.062144272,
            expected_palpepv=0.132056578,
            expected_pfuscmw=497.580075,
            expected_powfmw=2486.347575,
        ),
        Palph2Param(
            falpha=0.95,
            fdeut=2.5,
            ifalphap=0,
            bp=0.86,
            bt=5.3292,
            dene=7.432e19,
            deni=6.226e19,
            dnitot=6.743e19,
            falpe=0.68,
            falpi=0.32,
            palpnb=100.5,
            pchargepv=0.00066,
            ten=13.84,
            tin=13.84,
            vol=2426.25,
            palppv=0.163,
            pneutpv=0.654,
            expected_palppv=0.20442195,
            expected_pneutpv=0.81968779,
            expected_palpmw=495.97875,
            expected_pneutmw=1988.7675,
            expected_pchargemw=1.601325,
            expected_betaft=0.0,
            expected_palpipv=0.062144272,
            expected_palpepv=0.132056578,
            expected_pfuscmw=497.580075,
            expected_powfmw=2486.347575,
        ),
    ),
)
def test_palph2(palph2param, monkeypatch):
    """
    Automatically generated Regression Unit Test for palph2.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param palph2param: the data used to mock and assert in this test.
    :type palph2param: palph2param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(pv, "falpha", palph2param.falpha)
    monkeypatch.setattr(pv, "fdeut", palph2param.fdeut)

    # allow inout params to be mutated
    palppv = numpy.array(palph2param.palppv)
    pneutpv = numpy.array(palph2param.pneutpv)

    (palpmw, pneutmw, pchargemw, betaft, palpipv, palpepv, pfuscmw, powfmw) = pf.palph2(
        ifalphap=palph2param.ifalphap,
        bp=palph2param.bp,
        bt=palph2param.bt,
        dene=palph2param.dene,
        deni=palph2param.deni,
        dnitot=palph2param.dnitot,
        falpe=palph2param.falpe,
        falpi=palph2param.falpi,
        palpnb=palph2param.palpnb,
        pchargepv=palph2param.pchargepv,
        ten=palph2param.ten,
        tin=palph2param.tin,
        vol=palph2param.vol,
        palppv=palppv,
        pneutpv=pneutpv,
    )

    assert palppv == pytest.approx(palph2param.expected_palppv)
    assert pneutpv == pytest.approx(palph2param.expected_pneutpv)
    assert palpmw == pytest.approx(palph2param.expected_palpmw)
    assert pneutmw == pytest.approx(palph2param.expected_pneutmw)
    assert pchargemw == pytest.approx(palph2param.expected_pchargemw)
    assert betaft == pytest.approx(palph2param.expected_betaft)
    assert palpepv == pytest.approx(palph2param.expected_palpepv)
    assert palpipv == pytest.approx(palph2param.expected_palpipv)
    assert pfuscmw == pytest.approx(palph2param.expected_pfuscmw)
    assert powfmw == pytest.approx(palph2param.expected_powfmw)


@pytest.mark.parametrize(
    "t, reaction, expected_bosch_hale",
    (
        (0.0, -1, 0.0),
        (55.73, 1, 8.832857074192583e-22),
        (55.73, 2, 7.067916724597656e-23),
        (55.73, 3, 1.3127277533210717e-23),
        (55.73, 4, 1.1329338540436287e-23),
    ),
    ids=["t_0", "DT", "DHE3", "DD1", "DD2"],
)
def test_bosch_hale(t, reaction, expected_bosch_hale):
    """
    Unit test for the bosch_hale function.

    :param t: input Maxwellian density-weighted ion temperature
    :type t: float
    :param reaction: input flag for fusion reaction to use
    :type reaction: int
    :param expected_bosch_hale: expected return value from the bosch_hale function
    :type expected_bosch_hale: float
    """
    bosch_hale = pf.bosch_hale(t, reaction)

    assert bosch_hale == approx(expected_bosch_hale, abs=1e-23)


class PhalphParams(NamedTuple):
    alphan: float
    alphat: float
    deni: float
    fdeut: float
    fhe3: float
    ftrit: float
    ti: float

    ipedestal: int
    te: float
    rhopedt: float
    te0: float
    teped: float
    tesep: float
    tbeta: float
    dene: float
    rhopedn: float
    ne0: float
    neped: float
    nesep: float

    expected_palppv: float
    expected_pchargepv: float
    expected_pneutpv: float
    expected_sigvdt: float
    expected_fusionrate: float
    expected_alpharate: float
    expected_protonrate: float
    expected_pdtpv: float
    expected_pdhe3pv: float
    expected_pddpv: float


@pytest.mark.parametrize(
    "phalphparams",
    (
        PhalphParams(
            1.0,
            1.45,
            6.2262793637240177e19,
            0.5,
            0.0,
            0.5,
            13.84,
            1,
            12.33,
            0.94,
            28.089723663920328,
            3.7775374842470044,
            0.1,
            2.0,
            7.4321e19,
            0.94,
            9.7756974320342041e19,
            5.8300851381352219e19,
            3.4294618459618943e19,
            0.19030547335201128,
            0.000813064815368787,
            0.7616652065443165,
            3.48375533764882e-22,
            3.3979157349166214e17,
            3.3763297977777126e17,
            1030380967354505.6,
            0.9515273667600563,
            0.0,
            0.0012563779516400874,
        ),
    ),
)
def test_phalph(phalphparams, monkeypatch):
    """
    Automatically generated Integration for palph.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param palphparam: the data used to mock and assert in this test.
    :type palphparam: palphparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(pv, "ipedestal", phalphparams.ipedestal)
    monkeypatch.setattr(pv, "te", phalphparams.te)
    monkeypatch.setattr(pv, "rhopedt", phalphparams.rhopedt)
    monkeypatch.setattr(pv, "te0", phalphparams.te0)
    monkeypatch.setattr(pv, "teped", phalphparams.teped)
    monkeypatch.setattr(pv, "tesep", phalphparams.tesep)
    monkeypatch.setattr(pv, "tbeta", phalphparams.tbeta)
    monkeypatch.setattr(pv, "dene", phalphparams.dene)
    monkeypatch.setattr(pv, "rhopedn", phalphparams.rhopedn)
    monkeypatch.setattr(pv, "ne0", phalphparams.ne0)
    monkeypatch.setattr(pv, "neped", phalphparams.neped)
    monkeypatch.setattr(pv, "nesep", phalphparams.nesep)

    (
        palppv,
        pchargepv,
        pneutpv,
        sigvdt,
        fusionrate,
        alpharate,
        protonrate,
        pdtpv,
        pdhe3pv,
        pddpv,
    ) = pf.palph(
        phalphparams.alphan,
        phalphparams.alphat,
        phalphparams.deni,
        phalphparams.fdeut,
        phalphparams.fhe3,
        phalphparams.ftrit,
        phalphparams.ti,
    )

    assert palppv == pytest.approx(phalphparams.expected_palppv)
    assert pchargepv == pytest.approx(phalphparams.expected_pchargepv)
    assert pneutpv == pytest.approx(phalphparams.expected_pneutpv)
    assert sigvdt == pytest.approx(phalphparams.expected_sigvdt)
    assert fusionrate == pytest.approx(phalphparams.expected_fusionrate)
    assert alpharate == pytest.approx(phalphparams.expected_alpharate)
    assert protonrate == pytest.approx(phalphparams.expected_protonrate)
    assert pdtpv == pytest.approx(phalphparams.expected_pdtpv)
    assert pdhe3pv == pytest.approx(phalphparams.expected_pdhe3pv)
    assert pddpv == pytest.approx(phalphparams.expected_pddpv)
