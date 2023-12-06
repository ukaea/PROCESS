"""Unit tests for physics.f90."""
from typing import Any, NamedTuple
from process.fortran import physics_variables
import pytest
from process.physics import Physics
from process.plasma_profiles import PlasmaProfile


@pytest.fixture
def physics():
    """Provides Physics object for testing.

    :returns: initialised Physics object
    :rtype: process.physics.Physics
    """
    return Physics(PlasmaProfile())


def test_diamagnetic_fraction_hender(physics):
    """Test diamagnetic_fraction_hender()."""
    beta = 0.14
    diacf = physics.diamagnetic_fraction_hender(beta)
    assert diacf == pytest.approx(0.05, abs=0.0001)


def test_diamagnetic_fraction_scene(physics):
    """Test diamagnetic_fraction_scene."""
    beta = 0.15
    q95 = 3.0
    q0 = 1.0
    diacf = physics.diamagnetic_fraction_scene(beta, q95, q0)
    assert diacf == pytest.approx(0.0460, abs=0.0001)


def test_ps_fraction_scene(physics):
    """Test ps_fraction_scene."""
    beta = 0.15
    pscf = physics.ps_fraction_scene(beta)
    assert pscf == pytest.approx(-0.0135, abs=0.0001)


class BootstrapFractionIter89Param(NamedTuple):
    aspect: Any = None

    beta: Any = None

    bt: Any = None

    cboot: Any = None

    plascur: Any = None

    q95: Any = None

    q0: Any = None

    rmajor: Any = None

    vol: Any = None

    expected_bootipf: Any = None


@pytest.mark.parametrize(
    "bootstrapfractioniter89param",
    (
        BootstrapFractionIter89Param(
            aspect=3,
            beta=0.030000000000000006,
            bt=5.7802910787445487,
            cboot=1,
            plascur=18398455.678867526,
            q95=3.5,
            q0=1,
            rmajor=8,
            vol=1888.1711539956691,
            expected_bootipf=0.30255906256775245,
        ),
    ),
)
def test_bootstrap_fraction_iter89(bootstrapfractioniter89param, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_iter89.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param bootstrapfractioniter89param: the data used to mock and assert in this test.
    :type bootstrapfractioniter89param: bootstrapfractioniter89param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    bootipf = physics.bootstrap_fraction_iter89(
        aspect=bootstrapfractioniter89param.aspect,
        beta=bootstrapfractioniter89param.beta,
        bt=bootstrapfractioniter89param.bt,
        cboot=bootstrapfractioniter89param.cboot,
        plascur=bootstrapfractioniter89param.plascur,
        q95=bootstrapfractioniter89param.q95,
        q0=bootstrapfractioniter89param.q0,
        rmajor=bootstrapfractioniter89param.rmajor,
        vol=bootstrapfractioniter89param.vol,
    )

    assert bootipf == pytest.approx(bootstrapfractioniter89param.expected_bootipf)


class BootstrapFractionNevinsParam(NamedTuple):
    te0: Any = None

    ne0: Any = None

    alphan: Any = None

    betat: Any = None

    bt: Any = None

    dene: Any = None

    plascur: Any = None

    q0: Any = None

    q95: Any = None

    alphat: Any = None

    rmajor: Any = None

    rminor: Any = None

    ten: Any = None

    zeff: Any = None

    expected_fibs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionnevinsparam",
    (
        BootstrapFractionNevinsParam(
            te0=24.402321098330372,
            ne0=8.515060981068918e19,
            alphan=1.0,
            betat=0.03,
            bt=5.7,
            dene=18398455.678867526,
            plascur=18398455.678867526,
            q0=1,
            q95=3.5,
            alphat=1.45,
            rmajor=8,
            rminor=2.6666666666666665,
            ten=12.626131115905864,
            zeff=2.0909945616489103,
            expected_fibs=889258771342.7881,
        ),
    ),
)
def test_bootstrap_fraction_nevins(bootstrapfractionnevinsparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_nevins.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param bootstrapfractionnevinsparam: the data used to mock and assert in this test.
    :type bootstrapfractionnevinsparam: bootstrapfractionnevinsparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "te0", bootstrapfractionnevinsparam.te0)

    monkeypatch.setattr(physics_variables, "ne0", bootstrapfractionnevinsparam.ne0)

    fibs = physics.bootstrap_fraction_nevins(
        alphan=bootstrapfractionnevinsparam.alphan,
        alphat=bootstrapfractionnevinsparam.alphat,
        betat=bootstrapfractionnevinsparam.betat,
        bt=bootstrapfractionnevinsparam.bt,
        dene=bootstrapfractionnevinsparam.dene,
        plascur=bootstrapfractionnevinsparam.plascur,
        q0=bootstrapfractionnevinsparam.q0,
        q95=bootstrapfractionnevinsparam.q95,
        rmajor=bootstrapfractionnevinsparam.rmajor,
        rminor=bootstrapfractionnevinsparam.rminor,
        ten=bootstrapfractionnevinsparam.ten,
        zeff=bootstrapfractionnevinsparam.zeff,
    )

    assert fibs == pytest.approx(bootstrapfractionnevinsparam.expected_fibs)


class BootstrapFractionWilsonParam(NamedTuple):
    alphaj: Any = None

    alphap: Any = None

    alphat: Any = None

    betpth: Any = None

    q0: Any = None

    qpsi: Any = None

    rmajor: Any = None

    rminor: Any = None

    expected_bfw: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionwilsonparam",
    (
        BootstrapFractionWilsonParam(
            alphaj=1.9008029008029004,
            alphap=2.4500000000000002,
            alphat=1.45,
            betpth=1.0874279209664601,
            q0=1,
            qpsi=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
            expected_bfw=0.42321339288758714,
        ),
        BootstrapFractionWilsonParam(
            alphaj=1.9008029008029004,
            alphap=2.4500000000000002,
            alphat=1.45,
            betpth=0.99075943086768326,
            q0=1,
            qpsi=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
            expected_bfw=0.38559122143951252,
        ),
    ),
)
def test_bootstrap_fraction_wilson(bootstrapfractionwilsonparam, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_wilson.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param bootstrapfractionwilsonparam: the data used to mock and assert in this test.
    :type bootstrapfractionwilsonparam: bootstrapfractionwilsonparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    bfw = physics.bootstrap_fraction_wilson(
        alphaj=bootstrapfractionwilsonparam.alphaj,
        alphap=bootstrapfractionwilsonparam.alphap,
        alphat=bootstrapfractionwilsonparam.alphat,
        betpth=bootstrapfractionwilsonparam.betpth,
        q0=bootstrapfractionwilsonparam.q0,
        qpsi=bootstrapfractionwilsonparam.qpsi,
        rmajor=bootstrapfractionwilsonparam.rmajor,
        rminor=bootstrapfractionwilsonparam.rminor,
    )

    assert bfw == pytest.approx(bootstrapfractionwilsonparam.expected_bfw)


class BootstrapFractionSauterParam(NamedTuple):
    dnitot: Any = None

    rminor: Any = None

    tesep: Any = None

    ti: Any = None

    triang: Any = None

    q0: Any = None

    afuel: Any = None

    zeff: Any = None

    rhopedn: Any = None

    bt: Any = None

    plascur: Any = None

    xarea: Any = None

    fhe3: Any = None

    teped: Any = None

    dene: Any = None

    te: Any = None

    rmajor: Any = None

    q: Any = None

    nesep: Any = None

    te0: Any = None

    neped: Any = None

    tbeta: Any = None

    ne0: Any = None

    alphan: Any = None

    rhopedt: Any = None

    alphat: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionsauterparam",
    (
        BootstrapFractionSauterParam(
            dnitot=6.6125550702454276e19,
            rminor=2.6666666666666665,
            tesep=0.10000000000000001,
            ti=12,
            triang=0.5,
            q0=1,
            afuel=2.5,
            zeff=2.0909945616489103,
            rhopedn=0.94000000000000006,
            bt=5.7000000000000002,
            plascur=18398455.678867526,
            xarea=38.39822223637151,
            fhe3=0,
            teped=5.5,
            dene=7.5e19,
            te=12,
            rmajor=8,
            q=3.5,
            nesep=4.1177885154594193e19,
            te0=24.402321098330372,
            neped=7.000240476281013e19,
            tbeta=2,
            ne0=8.515060981068918e19,
            alphan=1,
            rhopedt=0.94000000000000006,
            alphat=1.45,
            expected_bfs=0.27635918746616817,
        ),
    ),
)
def test_bootstrap_fraction_sauter(bootstrapfractionsauterparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_sauter.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        physics_variables, "dnitot", bootstrapfractionsauterparam.dnitot
    )

    monkeypatch.setattr(
        physics_variables, "rminor", bootstrapfractionsauterparam.rminor
    )

    monkeypatch.setattr(physics_variables, "tesep", bootstrapfractionsauterparam.tesep)

    monkeypatch.setattr(physics_variables, "ti", bootstrapfractionsauterparam.ti)

    monkeypatch.setattr(
        physics_variables, "triang", bootstrapfractionsauterparam.triang
    )

    monkeypatch.setattr(physics_variables, "q0", bootstrapfractionsauterparam.q0)

    monkeypatch.setattr(physics_variables, "afuel", bootstrapfractionsauterparam.afuel)

    monkeypatch.setattr(physics_variables, "zeff", bootstrapfractionsauterparam.zeff)

    monkeypatch.setattr(
        physics_variables, "rhopedn", bootstrapfractionsauterparam.rhopedn
    )

    monkeypatch.setattr(physics_variables, "bt", bootstrapfractionsauterparam.bt)

    monkeypatch.setattr(
        physics_variables, "plascur", bootstrapfractionsauterparam.plascur
    )

    monkeypatch.setattr(physics_variables, "xarea", bootstrapfractionsauterparam.xarea)

    monkeypatch.setattr(physics_variables, "fhe3", bootstrapfractionsauterparam.fhe3)

    monkeypatch.setattr(physics_variables, "teped", bootstrapfractionsauterparam.teped)

    monkeypatch.setattr(physics_variables, "dene", bootstrapfractionsauterparam.dene)

    monkeypatch.setattr(physics_variables, "te", bootstrapfractionsauterparam.te)

    monkeypatch.setattr(
        physics_variables, "rmajor", bootstrapfractionsauterparam.rmajor
    )

    monkeypatch.setattr(physics_variables, "q", bootstrapfractionsauterparam.q)

    monkeypatch.setattr(physics_variables, "nesep", bootstrapfractionsauterparam.nesep)

    monkeypatch.setattr(physics_variables, "te0", bootstrapfractionsauterparam.te0)

    monkeypatch.setattr(physics_variables, "neped", bootstrapfractionsauterparam.neped)

    monkeypatch.setattr(physics_variables, "tbeta", bootstrapfractionsauterparam.tbeta)

    monkeypatch.setattr(physics_variables, "ne0", bootstrapfractionsauterparam.ne0)

    monkeypatch.setattr(
        physics_variables, "alphan", bootstrapfractionsauterparam.alphan
    )

    monkeypatch.setattr(
        physics_variables, "rhopedt", bootstrapfractionsauterparam.rhopedt
    )

    monkeypatch.setattr(
        physics_variables, "alphat", bootstrapfractionsauterparam.alphat
    )

    bfs = physics.bootstrap_fraction_sauter()

    assert bfs == pytest.approx(bootstrapfractionsauterparam.expected_bfs)


class CulcurParam(NamedTuple):
    normalised_total_beta: Any = None

    beta: Any = None

    icurr: Any = None

    iprofile: Any = None

    alphaj: Any = None

    rli: Any = None

    alphap: Any = None

    bt: Any = None

    eps: Any = None

    kappa: Any = None

    kappa95: Any = None

    p0: Any = None

    pperim: Any = None

    q0: Any = None

    qpsi: Any = None

    rmajor: Any = None

    rminor: Any = None

    sf: Any = None

    triang: Any = None

    triang95: Any = None

    expected_normalised_total_beta: Any = None

    expected_alphaj: Any = None

    expected_rli: Any = None

    expected_bp: Any = None

    expected_qstar: Any = None

    expected_plascur: Any = None


@pytest.mark.parametrize(
    "culcurparam",
    (
        CulcurParam(
            normalised_total_beta=0,
            beta=0.030000000000000006,
            icurr=4,
            iprofile=1,
            alphaj=1,
            rli=0.90000000000000002,
            alphap=0,
            bt=5.7000000000000002,
            eps=0.33333333333333331,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p0=0,
            pperim=24.081367139525412,
            q0=1,
            qpsi=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
            sf=1.4372507312498271,
            triang=0.5,
            triang95=0.33333333333333331,
            expected_normalised_total_beta=2.4784688886891844,
            expected_alphaj=1.9008029008029004,
            expected_rli=1.2064840230894305,
            expected_bp=0.96008591022564971,
            expected_qstar=2.9008029008029004,
            expected_plascur=18398455.678867526,
        ),
        CulcurParam(
            normalised_total_beta=2.4784688886891844,
            beta=0.030000000000000006,
            icurr=4,
            iprofile=1,
            alphaj=1.9008029008029004,
            rli=1.2064840230894305,
            alphap=2.4500000000000002,
            bt=5.7000000000000002,
            eps=0.33333333333333331,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p0=626431.90482713911,
            pperim=24.081367139525412,
            q0=1,
            qpsi=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
            sf=1.4372507312498271,
            triang=0.5,
            triang95=0.33333333333333331,
            expected_normalised_total_beta=2.4784688886891844,
            expected_alphaj=1.9008029008029004,
            expected_rli=1.2064840230894305,
            expected_bp=0.96008591022564971,
            expected_qstar=2.9008029008029004,
            expected_plascur=18398455.678867526,
        ),
    ),
)
def test_culcur(culcurparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for culcur.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param culcurparam: the data used to mock and assert in this test.
    :type culcurparam: culcurparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        physics_variables, "normalised_total_beta", culcurparam.normalised_total_beta
    )

    monkeypatch.setattr(physics_variables, "beta", culcurparam.beta)

    _, _, bp, qstar, plascur = physics.culcur(
        icurr=culcurparam.icurr,
        iprofile=culcurparam.iprofile,
        alphaj=culcurparam.alphaj,
        rli=culcurparam.rli,
        alphap=culcurparam.alphap,
        bt=culcurparam.bt,
        eps=culcurparam.eps,
        kappa=culcurparam.kappa,
        kappa95=culcurparam.kappa95,
        p0=culcurparam.p0,
        pperim=culcurparam.pperim,
        q0=culcurparam.q0,
        qpsi=culcurparam.qpsi,
        rmajor=culcurparam.rmajor,
        rminor=culcurparam.rminor,
        sf=culcurparam.sf,
        triang=culcurparam.triang,
        triang95=culcurparam.triang95,
    )

    assert physics_variables.normalised_total_beta == pytest.approx(
        culcurparam.expected_normalised_total_beta
    )

    assert bp == pytest.approx(culcurparam.expected_bp)

    assert qstar == pytest.approx(culcurparam.expected_qstar)

    assert plascur == pytest.approx(culcurparam.expected_plascur)
