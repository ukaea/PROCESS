"""Unit tests for physics.f90."""

from typing import Any, NamedTuple
from process.fortran import (
    constants,
    physics_variables,
    physics_module,
    current_drive_variables,
    impurity_radiation_module,
)
import numpy
import pytest
from process.physics import (
    Physics,
    calculate_poloidal_field,
    diamagnetic_fraction_scene,
    diamagnetic_fraction_hender,
    ps_fraction_scene,
    calculate_plasma_current_peng,
    culblm,
    calculate_current_coefficient_hastie,
    vscalc,
    rether,
    beta_poloidal,
    res_diff_time,
)
from process.plasma_profiles import PlasmaProfile
from process.current_drive import CurrentDrive
from process.impurity_radiation import initialise_imprad


@pytest.fixture
def physics():
    """Provides Physics object for testing.

    :returns: initialised Physics object
    :rtype: process.physics.Physics
    """
    return Physics(PlasmaProfile(), CurrentDrive(PlasmaProfile()))


def test_beta_poloidal():
    """Test beta_poloidal()"""
    betap = beta_poloidal(5.347, 0.852, 0.0307)
    assert betap == pytest.approx(1.209, abs=0.001)


def test_res_diff_time():
    """Test res_diff_time()"""
    res_time = res_diff_time(9.137, 2.909e-9, 1.65)
    assert res_time == pytest.approx(4784.3, abs=0.1)


def test_diamagnetic_fraction_hender():
    """Test diamagnetic_fraction_hender()."""
    beta = 0.14
    diacf = diamagnetic_fraction_hender(beta)
    assert diacf == pytest.approx(0.05, abs=0.0001)


def test_diamagnetic_fraction_scene():
    """Test diamagnetic_fraction_scene."""
    beta = 0.15
    q95 = 3.0
    q0 = 1.0
    diacf = diamagnetic_fraction_scene(beta, q95, q0)
    assert diacf == pytest.approx(0.0460, abs=0.0001)


def test_ps_fraction_scene():
    """Test ps_fraction_scene."""
    beta = 0.15
    pscf = ps_fraction_scene(beta)
    assert pscf == pytest.approx(-0.0135, abs=0.0001)


class BootstrapFractionIter89Param(NamedTuple):
    aspect: Any = None

    beta: Any = None

    bt: Any = None

    plasma_current: Any = None

    q95: Any = None

    q0: Any = None

    rmajor: Any = None

    plasma_volume: Any = None

    expected_bootipf: Any = None


@pytest.mark.parametrize(
    "bootstrapfractioniter89param",
    (
        BootstrapFractionIter89Param(
            aspect=3,
            beta=0.030000000000000006,
            bt=5.7802910787445487,
            plasma_current=18398455.678867526,
            q95=3.5,
            q0=1,
            rmajor=8,
            plasma_volume=1888.1711539956691,
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

    bootstrap_current_fraction = physics.bootstrap_fraction_iter89(
        aspect=bootstrapfractioniter89param.aspect,
        beta=bootstrapfractioniter89param.beta,
        bt=bootstrapfractioniter89param.bt,
        plasma_current=bootstrapfractioniter89param.plasma_current,
        q95=bootstrapfractioniter89param.q95,
        q0=bootstrapfractioniter89param.q0,
        rmajor=bootstrapfractioniter89param.rmajor,
        plasma_volume=bootstrapfractioniter89param.plasma_volume,
    )

    assert bootstrap_current_fraction == pytest.approx(
        bootstrapfractioniter89param.expected_bootipf
    )


class BootstrapFractionNevinsParam(NamedTuple):
    te0: Any = None

    ne0: Any = None

    alphan: Any = None

    betat: Any = None

    bt: Any = None

    dene: Any = None

    plasma_current: Any = None

    q0: Any = None

    q95: Any = None

    alphat: Any = None

    rmajor: Any = None

    rminor: Any = None

    te: Any = None

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
            plasma_current=18398455.678867526,
            q0=1,
            q95=3.5,
            alphat=1.45,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.626131115905864,
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
        plasma_current=bootstrapfractionnevinsparam.plasma_current,
        q0=bootstrapfractionnevinsparam.q0,
        q95=bootstrapfractionnevinsparam.q95,
        rmajor=bootstrapfractionnevinsparam.rmajor,
        rminor=bootstrapfractionnevinsparam.rminor,
        te=bootstrapfractionnevinsparam.te,
        zeff=bootstrapfractionnevinsparam.zeff,
    )

    assert fibs == pytest.approx(bootstrapfractionnevinsparam.expected_fibs)


class BootstrapFractionWilsonParam(NamedTuple):
    alphaj: Any = None

    alphap: Any = None

    alphat: Any = None

    betpth: Any = None

    q0: Any = None

    q95: Any = None

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
            q95=3.5,
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
            q95=3.5,
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
        q95=bootstrapfractionwilsonparam.q95,
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

    plasma_current: Any = None

    xarea: Any = None

    f_helium3: Any = None

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
            dnitot=7.1297522422781575e19,
            rminor=2.6666666666666665,
            tesep=0.10000000000000001,
            ti=12.570861186498382,
            triang=0.5,
            q0=1,
            afuel=2.5,
            zeff=2.5211399464385624,
            rhopedn=0.9400000000000001,
            bt=5.326133750416047,
            plasma_current=16528278.760008096,
            xarea=38.39822223637151,
            f_helium3=0,
            teped=5.5,
            dene=8.016748468651018e19,
            te=12.570861186498382,
            rmajor=8,
            q=3.5,
            nesep=3.6992211545476006e19,
            te0=25.986118047669795,
            neped=6.2886759627309195e19,
            tbeta=2,
            ne0=1.054474759840606e20,
            alphan=1,
            rhopedt=0.9400000000000001,
            alphat=1.45,
            expected_bfs=0.4110838247346975,
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
        physics_variables, "plasma_current", bootstrapfractionsauterparam.plasma_current
    )

    monkeypatch.setattr(physics_variables, "xarea", bootstrapfractionsauterparam.xarea)

    monkeypatch.setattr(
        physics_variables, "f_helium3", bootstrapfractionsauterparam.f_helium3
    )

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
    physics.plasma_profile.run()
    bfs = physics.bootstrap_fraction_sauter(physics.plasma_profile)

    assert bfs == pytest.approx(bootstrapfractionsauterparam.expected_bfs)


class BootstrapFractionSakaiParam(NamedTuple):
    betap: Any = None

    q95: Any = None

    q0: Any = None

    alphan: Any = None

    alphat: Any = None

    eps: Any = None

    rli: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionsakaiparam",
    (
        BootstrapFractionSakaiParam(
            betap=1.3184383457774960,
            q95=3.5151046634673557,
            q0=1.0,
            alphan=1.0,
            alphat=1.45,
            eps=1 / 3,
            rli=1.2098126022585098,
            expected_bfs=0.3501274900057279,
        ),
        BootstrapFractionSakaiParam(
            betap=1.1701245502231756,
            q95=5.1746754543339177,
            q0=2.0,
            alphan=0.9,
            alphat=1.3999999999999999,
            eps=1 / 1.8,
            rli=0.3,
            expected_bfs=0.81877746774625,
        ),
    ),
)
def test_bootstrap_fraction_sakai(bootstrapfractionsakaiparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_sakai.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.
    and tests/regression/input_files/st_regression.IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "betap", bootstrapfractionsakaiparam.betap)

    monkeypatch.setattr(physics_variables, "q95", bootstrapfractionsakaiparam.q95)

    monkeypatch.setattr(physics_variables, "q0", bootstrapfractionsakaiparam.q0)

    monkeypatch.setattr(physics_variables, "alphan", bootstrapfractionsakaiparam.alphan)

    monkeypatch.setattr(physics_variables, "alphat", bootstrapfractionsakaiparam.alphat)

    monkeypatch.setattr(physics_variables, "eps", bootstrapfractionsakaiparam.eps)

    monkeypatch.setattr(physics_variables, "rli", bootstrapfractionsakaiparam.rli)

    bfs = physics.bootstrap_fraction_sakai(
        betap=bootstrapfractionsakaiparam.betap,
        q95=bootstrapfractionsakaiparam.q95,
        q0=bootstrapfractionsakaiparam.q0,
        alphan=bootstrapfractionsakaiparam.alphan,
        alphat=bootstrapfractionsakaiparam.alphat,
        eps=bootstrapfractionsakaiparam.eps,
        rli=bootstrapfractionsakaiparam.rli,
    )

    assert bfs == pytest.approx(bootstrapfractionsakaiparam.expected_bfs)


class BootstrapFractionAriesParam(NamedTuple):
    betap: Any = None

    rli: Any = None

    core_density: Any = None

    average_density: Any = None

    inverse_aspect: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionariesparam",
    (
        BootstrapFractionAriesParam(
            betap=1.2708883332338736,
            rli=1.4279108047138775,
            core_density=1.0695994460047332e20,
            average_density=8.1317358967210131e19,
            inverse_aspect=1 / 3,
            expected_bfs=4.3237405809568441e-01,
        ),
    ),
)
def test_bootstrap_fraction_aries(bootstrapfractionariesparam, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_aries.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam
    """

    bfs = physics.bootstrap_fraction_aries(
        betap=bootstrapfractionariesparam.betap,
        rli=bootstrapfractionariesparam.rli,
        core_density=bootstrapfractionariesparam.core_density,
        average_density=bootstrapfractionariesparam.average_density,
        inverse_aspect=bootstrapfractionariesparam.inverse_aspect,
    )

    assert bfs == pytest.approx(bootstrapfractionariesparam.expected_bfs)


class BootstrapFractionAndradeParam(NamedTuple):
    betap: Any = None

    core_pressure: Any = None

    average_pressure: Any = None

    inverse_aspect: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionandradeparam",
    (
        BootstrapFractionAndradeParam(
            betap=1.2708883332338736,
            core_pressure=8.3049163275475602e05,
            average_pressure=2.4072221239268288e05,
            inverse_aspect=1 / 3,
            expected_bfs=4.6240007834873120e-01,
        ),
    ),
)
def test_bootstrap_fraction_andrade(bootstrapfractionandradeparam, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_andrade.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam
    """

    bfs = physics.bootstrap_fraction_andrade(
        betap=bootstrapfractionandradeparam.betap,
        core_pressure=bootstrapfractionandradeparam.core_pressure,
        average_pressure=bootstrapfractionandradeparam.average_pressure,
        inverse_aspect=bootstrapfractionandradeparam.inverse_aspect,
    )

    assert bfs == pytest.approx(bootstrapfractionandradeparam.expected_bfs)


class BootstrapFractionHoangParam(NamedTuple):
    betap: Any = None

    pressure_index: Any = None

    current_index: Any = None

    inverse_aspect: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionhoangparam",
    (
        BootstrapFractionHoangParam(
            betap=1.2708883332338736,
            pressure_index=2.4500000000000002e00,
            current_index=2.8314361644755763e00,
            inverse_aspect=1 / 3,
            expected_bfs=0.27190974213794156,
        ),
    ),
)
def test_bootstrap_fraction_hoang(bootstrapfractionhoangparam, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_hoang.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam
    """

    bfs = physics.bootstrap_fraction_hoang(
        betap=bootstrapfractionhoangparam.betap,
        pressure_index=bootstrapfractionhoangparam.pressure_index,
        current_index=bootstrapfractionhoangparam.current_index,
        inverse_aspect=bootstrapfractionhoangparam.inverse_aspect,
    )

    assert bfs == pytest.approx(bootstrapfractionhoangparam.expected_bfs)


class BootstrapFractionWongParam(NamedTuple):
    betap: Any = None

    density_index: Any = None

    temperature_index: Any = None

    inverse_aspect: Any = None

    elongation: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionwongparam",
    (
        BootstrapFractionWongParam(
            betap=1.2708883332338736,
            density_index=1.0000000000000000e00,
            temperature_index=1.4500000000000000e00,
            inverse_aspect=1 / 3,
            elongation=1.8500000000000001e00,
            expected_bfs=7.0706527916080808e-01,
        ),
    ),
)
def test_bootstrap_fraction_wong(bootstrapfractionwongparam, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_wong.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam
    """

    bfs = physics.bootstrap_fraction_wong(
        betap=bootstrapfractionwongparam.betap,
        density_index=bootstrapfractionwongparam.density_index,
        temperature_index=bootstrapfractionwongparam.temperature_index,
        inverse_aspect=bootstrapfractionwongparam.inverse_aspect,
        elongation=bootstrapfractionwongparam.elongation,
    )

    assert bfs == pytest.approx(bootstrapfractionwongparam.expected_bfs)


class BootstrapFractionGiIParam(NamedTuple):
    betap: Any = None

    pressure_index: Any = None

    temperature_index: Any = None

    inverse_aspect: Any = None

    effective_charge: Any = None

    q95: Any = None

    q0: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractiongiiparam",
    (
        BootstrapFractionGiIParam(
            betap=1.2708883332338736,
            pressure_index=2.4500000000000002e00,
            temperature_index=1.4500000000000000e00,
            inverse_aspect=1 / 3,
            effective_charge=2.5368733516769737e00,
            q95=3.4656394133756647e00,
            q0=1.0,
            expected_bfs=7.9639753138719782e-01,
        ),
    ),
)
def test_bootstrap_fraction_gi_I(bootstrapfractiongiiparam, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_gi.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam
    """

    bfs = physics.bootstrap_fraction_gi_I(
        betap=bootstrapfractiongiiparam.betap,
        pressure_index=bootstrapfractiongiiparam.pressure_index,
        temperature_index=bootstrapfractiongiiparam.temperature_index,
        inverse_aspect=bootstrapfractiongiiparam.inverse_aspect,
        effective_charge=bootstrapfractiongiiparam.effective_charge,
        q95=bootstrapfractiongiiparam.q95,
        q0=bootstrapfractiongiiparam.q0,
    )

    assert bfs == pytest.approx(bootstrapfractiongiiparam.expected_bfs)


class BootstrapFractionGiIIParam(NamedTuple):
    betap: Any = None

    pressure_index: Any = None

    temperature_index: Any = None

    inverse_aspect: Any = None

    effective_charge: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractiongiiiparam",
    (
        BootstrapFractionGiIIParam(
            betap=1.2708883332338736,
            pressure_index=2.4500000000000002e00,
            temperature_index=1.4500000000000000e00,
            inverse_aspect=1 / 3,
            effective_charge=2.5368733516769737e00,
            expected_bfs=8.8502865710180589e-01,
        ),
    ),
)
def test_bootstrap_fraction_gi_II(bootstrapfractiongiiiparam, physics):
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_gi.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam
    """

    bfs = physics.bootstrap_fraction_gi_II(
        betap=bootstrapfractiongiiiparam.betap,
        pressure_index=bootstrapfractiongiiiparam.pressure_index,
        temperature_index=bootstrapfractiongiiiparam.temperature_index,
        inverse_aspect=bootstrapfractiongiiiparam.inverse_aspect,
        effective_charge=bootstrapfractiongiiiparam.effective_charge,
    )

    assert bfs == pytest.approx(bootstrapfractiongiiiparam.expected_bfs)


class PlasmaCurrentParam(NamedTuple):
    normalised_total_beta: Any = None

    beta: Any = None

    i_plasma_current: Any = None

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

    q95: Any = None

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

    expected_plasma_current: Any = None


@pytest.mark.parametrize(
    "plasmacurrentparam",
    (
        PlasmaCurrentParam(
            normalised_total_beta=0,
            beta=0.030000000000000006,
            i_plasma_current=4,
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
            q95=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
            sf=1.4372507312498271,
            triang=0.5,
            triang95=0.33333333333333331,
            expected_normalised_total_beta=2.4784688886891844,
            expected_alphaj=1.9008029008029004,
            expected_rli=1.2064840230894305,
            expected_bp=0.96008591022564971,
            expected_qstar=3.869423496255382,
            expected_plasma_current=18398455.678867526,
        ),
        PlasmaCurrentParam(
            normalised_total_beta=2.4784688886891844,
            beta=0.030000000000000006,
            i_plasma_current=4,
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
            q95=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
            sf=1.4372507312498271,
            triang=0.5,
            triang95=0.33333333333333331,
            expected_normalised_total_beta=2.4784688886891844,
            expected_alphaj=1.9008029008029004,
            expected_rli=1.2064840230894305,
            expected_bp=0.96008591022564971,
            expected_qstar=3.869423496255382,
            expected_plasma_current=18398455.678867526,
        ),
    ),
)
def test_calculate_plasma_current(plasmacurrentparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for calculate_plasma_current().

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param plasmacurrentparam: the data used to mock and assert in this test.
    :type plasmacurrentparam: plasmacurrentparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        physics_variables,
        "normalised_total_beta",
        plasmacurrentparam.normalised_total_beta,
    )

    monkeypatch.setattr(physics_variables, "beta", plasmacurrentparam.beta)

    _, _, bp, qstar, plasma_current = physics.calculate_plasma_current(
        i_plasma_current=plasmacurrentparam.i_plasma_current,
        iprofile=plasmacurrentparam.iprofile,
        alphaj=plasmacurrentparam.alphaj,
        rli=plasmacurrentparam.rli,
        alphap=plasmacurrentparam.alphap,
        bt=plasmacurrentparam.bt,
        eps=plasmacurrentparam.eps,
        kappa=plasmacurrentparam.kappa,
        kappa95=plasmacurrentparam.kappa95,
        p0=plasmacurrentparam.p0,
        pperim=plasmacurrentparam.pperim,
        q0=plasmacurrentparam.q0,
        q95=plasmacurrentparam.q95,
        rmajor=plasmacurrentparam.rmajor,
        rminor=plasmacurrentparam.rminor,
        sf=plasmacurrentparam.sf,
        triang=plasmacurrentparam.triang,
        triang95=plasmacurrentparam.triang95,
    )

    assert physics_variables.normalised_total_beta == pytest.approx(
        plasmacurrentparam.expected_normalised_total_beta
    )

    assert bp == pytest.approx(plasmacurrentparam.expected_bp)

    assert qstar == pytest.approx(plasmacurrentparam.expected_qstar)

    assert plasma_current == pytest.approx(plasmacurrentparam.expected_plasma_current)


@pytest.mark.parametrize(
    ("arguments", "expected"),
    (
        (
            {
                "q95": 2.5,
                "aspect": 2.7,
                "eps": 0.37037037,
                "rminor": 1.5,
                "bt": 12,
                "kappa": 1.85,
                "delta": 0.5,
            },
            46.84050744522757,
        ),
        (
            {
                "q95": 2.5,
                "aspect": 3.0,
                "eps": 0.33333333,
                "rminor": 1.5,
                "bt": 12,
                "kappa": 1.85,
                "delta": 0.5,
            },
            39.90862317467305,
        ),
    ),
)
def test_calculate_plasma_current_peng(arguments, expected):
    assert calculate_plasma_current_peng(**arguments) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("arguments", "expected"),
    (
        (
            {
                "i_plasma_current": 2,
                "ip": 1.6e7,
                "q95": 2.5,
                "aspect": 2.7,
                "eps": 0.37037037,
                "bt": 12,
                "kappa": 1.85,
                "delta": 0.5,
                "perim": 24,
                "rmu0": constants.rmu0,
            },
            4.3453802853633166,
        ),
        (
            {
                "i_plasma_current": 2,
                "ip": 1.6e7,
                "q95": 2.5,
                "aspect": 3.0,
                "eps": 0.33333333,
                "bt": 12,
                "kappa": 1.85,
                "delta": 0.5,
                "perim": 24,
                "rmu0": constants.rmu0,
            },
            3.702311392804667,
        ),
        (
            {
                "i_plasma_current": 3,
                "ip": 1.6e7,
                "q95": 2.5,
                "aspect": 3.0,
                "eps": 0.33333333,
                "bt": 12,
                "kappa": 1.85,
                "delta": 0.5,
                "perim": 24,
                "rmu0": constants.rmu0,
            },
            0.8377580413333333,
        ),
    ),
)
def test_calculate_poloidal_field(arguments, expected):
    assert calculate_poloidal_field(**arguments) == pytest.approx(expected)


def test_culblm():
    assert culblm(12, 4.879, 18300000, 2.5) == pytest.approx(0.0297619)


def test_conhas():
    assert calculate_current_coefficient_hastie(
        5, 5, 12, 0.5, 0.33, 1.85, 2e3, constants.rmu0
    ) == pytest.approx(2.518876726889116)


class PlasmaCompositionParam(NamedTuple):
    f_tritium_beam: Any = None

    impurity_arr_frac: Any = None

    impurity_arr_z: Any = None

    impurity_arr_amass: Any = None

    alphat: Any = None

    ignite: Any = None

    f_alpha_electron: Any = None

    afuel: Any = None

    f_tritium: Any = None

    deni: Any = None

    aion: Any = None

    dnitot: Any = None

    protium: Any = None

    zeffai: Any = None

    rncne: Any = None

    rnone: Any = None

    f_alpha_ion: Any = None

    ralpne: Any = None

    dlamee: Any = None

    rnbeam: Any = None

    zeff: Any = None

    dnz: Any = None

    pcoef: Any = None

    alpha_rate_density_total: Any = None

    rnfene: Any = None

    abeam: Any = None

    dlamie: Any = None

    te: Any = None

    proton_rate_density: Any = None

    f_deuterium: Any = None

    alphan: Any = None

    dnbeam: Any = None

    f_helium3: Any = None

    dnalp: Any = None

    dene: Any = None

    dnprot: Any = None

    iscz: Any = None

    err242: Any = None

    err243: Any = None

    ptarmw: Any = None

    lambdaio: Any = None

    drsep: Any = None

    fio: Any = None

    rho_star: Any = None

    nu_star: Any = None

    beta_mcdonald: Any = None

    itart_r: Any = None

    first_call: Any = None

    expected_impurity_arr_frac: Any = None

    expected_f_alpha_electron: Any = None

    expected_afuel: Any = None

    expected_deni: Any = None

    expected_aion: Any = None

    expected_dnitot: Any = None

    expected_zeffai: Any = None

    expected_f_alpha_ion: Any = None

    expected_dlamee: Any = None

    expected_zeff: Any = None

    expected_dnz: Any = None

    expected_abeam: Any = None

    expected_dlamie: Any = None

    expected_dnalp: Any = None

    expected_dnprot: Any = None

    expected_first_call: Any = None


@pytest.mark.parametrize(
    "plasmacompositionparam",
    (
        PlasmaCompositionParam(
            f_tritium_beam=9.9999999999999995e-07,
            impurity_arr_frac=[
                0.90000000000000002,
                0.10000000000000001,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.00038000000000000008,
                5.0000000000000021e-06,
            ],
            impurity_arr_z=[1, 2, 4, 6, 7, 8, 10, 14, 18, 26, 28, 36, 54, 74],
            impurity_arr_amass=[
                1.01,
                4.0030000000000001,
                9.0099999999999998,
                12.01,
                14.01,
                15.999000000000001,
                20.18,
                28.09,
                39.950000000000003,
                55.850000000000001,
                58.700000000000003,
                83.799999999999997,
                131.30000000000001,
                183.84999999999999,
            ],
            alphat=1.45,
            ignite=0,
            f_alpha_electron=0,
            afuel=0,
            f_tritium=0.5,
            deni=0,
            aion=0,
            dnitot=0,
            protium=0,
            zeffai=0,
            rncne=0,
            rnone=0,
            f_alpha_ion=0,
            ralpne=0.10000000000000001,
            dlamee=0,
            rnbeam=0,
            zeff=0,
            dnz=0,
            pcoef=0,
            alpha_rate_density_total=0,
            rnfene=0,
            abeam=0,
            dlamie=0,
            te=12,
            proton_rate_density=0,
            f_deuterium=0.5,
            alphan=1,
            dnbeam=0,
            f_helium3=0,
            dnalp=0,
            dene=7.5e19,
            dnprot=0,
            iscz=0,
            err242=0,
            err243=0,
            ptarmw=0,
            lambdaio=0,
            drsep=0,
            fio=0,
            rho_star=0,
            nu_star=0,
            beta_mcdonald=0,
            itart_r=0,
            first_call=1,
            expected_impurity_arr_frac=[
                0.78128900936605694,
                0.10000000000000001,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.00038000000000000008,
                5.0000000000000021e-06,
            ],
            expected_f_alpha_electron=0.6845930883190634,
            expected_afuel=2.5,
            expected_deni=5.8589175702454272e19,
            expected_aion=2.7265017998473029,
            expected_dnitot=6.6125550702454276e19,
            expected_zeffai=0.43248858851447464,
            expected_f_alpha_ion=0.3154069116809366,
            expected_dlamee=17.510652035055571,
            expected_zeff=2.0909945616489103,
            expected_dnz=28875000000000004,
            expected_abeam=2.0000010000000001,
            expected_dlamie=17.810652035055568,
            expected_dnalp=7.5e18,
            expected_dnprot=7500000000000000,
            expected_first_call=0,
        ),
        PlasmaCompositionParam(
            f_tritium_beam=9.9999999999999995e-07,
            impurity_arr_frac=(
                0.78128900936605694,
                0.10000000000000001,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.00038000000000000008,
                5.0000000000000021e-06,
            ),
            impurity_arr_z=numpy.array(
                numpy.array(
                    (1, 2, 4, 6, 7, 8, 10, 14, 18, 26, 28, 36, 54, 74), order="F"
                ),
                order="F",
            ).transpose(),
            impurity_arr_amass=numpy.array(
                numpy.array(
                    (
                        1.01,
                        4.0030000000000001,
                        9.0099999999999998,
                        12.01,
                        14.01,
                        15.999000000000001,
                        20.18,
                        28.09,
                        39.950000000000003,
                        55.850000000000001,
                        58.700000000000003,
                        83.799999999999997,
                        131.30000000000001,
                        183.84999999999999,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            alphat=1.45,
            ignite=0,
            f_alpha_electron=0.6845930883190634,
            afuel=2.5,
            f_tritium=0.5,
            deni=5.8589175702454272e19,
            aion=2.7265017998473029,
            dnitot=6.6125550702454276e19,
            protium=0,
            zeffai=0.43248858851447464,
            rncne=0,
            rnone=0,
            f_alpha_ion=0.3154069116809366,
            ralpne=0.10000000000000001,
            dlamee=17.510652035055571,
            rnbeam=0,
            zeff=2.0909945616489103,
            dnz=28875000000000004,
            pcoef=1.0521775929921553,
            alpha_rate_density_total=1.973996644759543e17,
            rnfene=0,
            abeam=2.0000010000000001,
            dlamie=17.810652035055568,
            te=12,
            proton_rate_density=540072280299564.38,
            f_deuterium=0.5,
            alphan=1,
            dnbeam=0,
            f_helium3=0,
            dnalp=7.5e18,
            dene=7.5e19,
            dnprot=7500000000000000,
            iscz=0,
            err242=0,
            err243=0,
            ptarmw=33.990985729118783,
            lambdaio=0.00157,
            drsep=-0.014999999999999999,
            fio=0.40999999999999998,
            rho_star=0,
            nu_star=0,
            beta_mcdonald=0,
            itart_r=0,
            first_call=0,
            expected_impurity_arr_frac=(
                0.78128900936605694,
                0.10000000000000001,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0.00038000000000000008,
                5.0000000000000021e-06,
            ),
            expected_f_alpha_electron=0.73096121787894142,
            expected_afuel=2.5,
            expected_deni=5.8576156204039725e19,
            expected_aion=2.7262064639685937,
            expected_dnitot=6.6125550702454276e19,
            expected_zeffai=0.43258985127992111,
            expected_f_alpha_ion=0.26903878212105858,
            expected_dlamee=17.510652035055571,
            expected_zeff=2.0909945616489103,
            expected_dnz=28875000000000004,
            expected_abeam=2.0000010000000001,
            expected_dlamie=17.810652035055568,
            expected_dnalp=7.5e18,
            expected_dnprot=20519498414548412,
            expected_first_call=0,
        ),
    ),
)
def test_plasma_composition(plasmacompositionparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for plasma_composition.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param plasmacompositionparam: the data used to mock and assert in this test.
    :type plasmacompositionparam: plasmacompositionparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    impurity_radiation_module.init_impurity_radiation_module()
    initialise_imprad()

    monkeypatch.setattr(
        current_drive_variables, "f_tritium_beam", plasmacompositionparam.f_tritium_beam
    )

    monkeypatch.setattr(
        impurity_radiation_module,
        "impurity_arr_frac",
        plasmacompositionparam.impurity_arr_frac,
    )

    monkeypatch.setattr(
        impurity_radiation_module,
        "impurity_arr_z",
        plasmacompositionparam.impurity_arr_z,
    )

    monkeypatch.setattr(
        impurity_radiation_module,
        "impurity_arr_amass",
        plasmacompositionparam.impurity_arr_amass,
    )

    monkeypatch.setattr(physics_variables, "alphat", plasmacompositionparam.alphat)

    monkeypatch.setattr(physics_variables, "ignite", plasmacompositionparam.ignite)

    monkeypatch.setattr(
        physics_variables, "f_alpha_electron", plasmacompositionparam.f_alpha_electron
    )

    monkeypatch.setattr(physics_variables, "afuel", plasmacompositionparam.afuel)

    monkeypatch.setattr(
        physics_variables, "f_tritium", plasmacompositionparam.f_tritium
    )

    monkeypatch.setattr(physics_variables, "deni", plasmacompositionparam.deni)

    monkeypatch.setattr(physics_variables, "aion", plasmacompositionparam.aion)

    monkeypatch.setattr(physics_variables, "dnitot", plasmacompositionparam.dnitot)

    monkeypatch.setattr(physics_variables, "protium", plasmacompositionparam.protium)

    monkeypatch.setattr(physics_variables, "zeffai", plasmacompositionparam.zeffai)

    monkeypatch.setattr(physics_variables, "rncne", plasmacompositionparam.rncne)

    monkeypatch.setattr(physics_variables, "rnone", plasmacompositionparam.rnone)

    monkeypatch.setattr(
        physics_variables, "f_alpha_ion", plasmacompositionparam.f_alpha_ion
    )

    monkeypatch.setattr(physics_variables, "ralpne", plasmacompositionparam.ralpne)

    monkeypatch.setattr(physics_variables, "dlamee", plasmacompositionparam.dlamee)

    monkeypatch.setattr(physics_variables, "rnbeam", plasmacompositionparam.rnbeam)

    monkeypatch.setattr(physics_variables, "zeff", plasmacompositionparam.zeff)

    monkeypatch.setattr(physics_variables, "dnz", plasmacompositionparam.dnz)

    monkeypatch.setattr(physics_variables, "pcoef", plasmacompositionparam.pcoef)

    monkeypatch.setattr(
        physics_variables,
        "alpha_rate_density_total",
        plasmacompositionparam.alpha_rate_density_total,
    )

    monkeypatch.setattr(physics_variables, "rnfene", plasmacompositionparam.rnfene)

    monkeypatch.setattr(physics_variables, "abeam", plasmacompositionparam.abeam)

    monkeypatch.setattr(physics_variables, "dlamie", plasmacompositionparam.dlamie)

    monkeypatch.setattr(physics_variables, "te", plasmacompositionparam.te)

    monkeypatch.setattr(
        physics_variables,
        "proton_rate_density",
        plasmacompositionparam.proton_rate_density,
    )

    monkeypatch.setattr(
        physics_variables, "f_deuterium", plasmacompositionparam.f_deuterium
    )

    monkeypatch.setattr(physics_variables, "alphan", plasmacompositionparam.alphan)

    monkeypatch.setattr(physics_variables, "dnbeam", plasmacompositionparam.dnbeam)

    monkeypatch.setattr(
        physics_variables, "f_helium3", plasmacompositionparam.f_helium3
    )

    monkeypatch.setattr(physics_variables, "dnalp", plasmacompositionparam.dnalp)

    monkeypatch.setattr(physics_variables, "dene", plasmacompositionparam.dene)

    monkeypatch.setattr(physics_variables, "dnprot", plasmacompositionparam.dnprot)

    monkeypatch.setattr(physics_module, "iscz", plasmacompositionparam.iscz)

    monkeypatch.setattr(physics_module, "err242", plasmacompositionparam.err242)

    monkeypatch.setattr(physics_module, "err243", plasmacompositionparam.err243)

    monkeypatch.setattr(physics_module, "ptarmw", plasmacompositionparam.ptarmw)

    monkeypatch.setattr(physics_module, "lambdaio", plasmacompositionparam.lambdaio)

    monkeypatch.setattr(physics_module, "drsep", plasmacompositionparam.drsep)

    monkeypatch.setattr(physics_module, "fio", plasmacompositionparam.fio)

    monkeypatch.setattr(physics_module, "rho_star", plasmacompositionparam.rho_star)

    monkeypatch.setattr(physics_module, "nu_star", plasmacompositionparam.nu_star)

    monkeypatch.setattr(
        physics_module, "beta_mcdonald", plasmacompositionparam.beta_mcdonald
    )

    monkeypatch.setattr(physics_module, "itart_r", plasmacompositionparam.itart_r)

    monkeypatch.setattr(physics_module, "first_call", plasmacompositionparam.first_call)

    physics.plasma_composition()

    assert impurity_radiation_module.impurity_arr_frac == pytest.approx(
        plasmacompositionparam.expected_impurity_arr_frac
    )

    assert physics_variables.f_alpha_electron == pytest.approx(
        plasmacompositionparam.expected_f_alpha_electron
    )

    assert physics_variables.afuel == pytest.approx(
        plasmacompositionparam.expected_afuel
    )

    assert physics_variables.deni == pytest.approx(plasmacompositionparam.expected_deni)

    assert physics_variables.aion == pytest.approx(plasmacompositionparam.expected_aion)

    assert physics_variables.dnitot == pytest.approx(
        plasmacompositionparam.expected_dnitot
    )

    assert physics_variables.zeffai == pytest.approx(
        plasmacompositionparam.expected_zeffai
    )

    assert physics_variables.f_alpha_ion == pytest.approx(
        plasmacompositionparam.expected_f_alpha_ion
    )

    assert physics_variables.dlamee == pytest.approx(
        plasmacompositionparam.expected_dlamee
    )

    assert physics_variables.zeff == pytest.approx(plasmacompositionparam.expected_zeff)

    assert physics_variables.dnz == pytest.approx(plasmacompositionparam.expected_dnz)

    assert physics_variables.abeam == pytest.approx(
        plasmacompositionparam.expected_abeam
    )

    assert physics_variables.dlamie == pytest.approx(
        plasmacompositionparam.expected_dlamie
    )

    assert physics_variables.dnalp == pytest.approx(
        plasmacompositionparam.expected_dnalp
    )

    assert physics_variables.dnprot == pytest.approx(
        plasmacompositionparam.expected_dnprot
    )

    assert physics_module.first_call == pytest.approx(
        plasmacompositionparam.expected_first_call
    )


class VscalcParam(NamedTuple):
    csawth: Any = None

    eps: Any = None

    inductive_current_fraction: Any = None

    gamma: Any = None

    kappa: Any = None

    plasma_current: Any = None

    rli: Any = None

    rmajor: Any = None

    rplas: Any = None

    t_burn: Any = None

    t_fusion_ramp: Any = None

    expected_phiint: Any = None

    expected_rlp: Any = None

    expected_vsbrn: Any = None

    expected_vsind: Any = None

    expected_vsres: Any = None

    expected_vsstt: Any = None


@pytest.mark.parametrize(
    "vscalcparam",
    (
        VscalcParam(
            csawth=1,
            eps=0.33333333333333331,
            inductive_current_fraction=0.59999999999999998,
            gamma=0.30000000000000004,
            kappa=1.8500000000000001,
            plasma_current=18398455.678867526,
            rli=1.2064840230894305,
            rmajor=8,
            rplas=3.7767895536275952e-09,
            t_burn=1000,
            t_fusion_ramp=10,
            expected_phiint=111.57651734747576,
            expected_rlp=1.4075705307248088e-05,
            expected_vsbrn=42.109179697761263,
            expected_vsind=258.97124024420435,
            expected_vsres=55.488435095110333,
            expected_vsstt=356.56885503707593,
        ),
        VscalcParam(
            csawth=1,
            eps=0.33333333333333331,
            inductive_current_fraction=0.59999999999999998,
            gamma=0.30000000000000004,
            kappa=1.8500000000000001,
            plasma_current=18398455.678867526,
            rli=1.2064840230894305,
            rmajor=8,
            rplas=3.7767895536275952e-09,
            t_burn=0,
            t_fusion_ramp=10,
            expected_phiint=111.57651734747576,
            expected_rlp=1.4075705307248088e-05,
            expected_vsbrn=0.41692257126496302,
            expected_vsind=258.97124024420435,
            expected_vsres=55.488435095110333,
            expected_vsstt=314.87659791057968,
        ),
    ),
)
def test_vscalc(vscalcparam):
    """
    Automatically generated Regression Unit Test for vscalc.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param vscalcparam: the data used to mock and assert in this test.
    :type vscalcparam: vscalcparam
    """

    phiint, rlp, vsbrn, vsind, vsres, vsstt = vscalc(
        csawth=vscalcparam.csawth,
        eps=vscalcparam.eps,
        inductive_current_fraction=vscalcparam.inductive_current_fraction,
        gamma=vscalcparam.gamma,
        kappa=vscalcparam.kappa,
        plasma_current=vscalcparam.plasma_current,
        rli=vscalcparam.rli,
        rmajor=vscalcparam.rmajor,
        rplas=vscalcparam.rplas,
        t_burn=vscalcparam.t_burn,
        t_fusion_ramp=vscalcparam.t_fusion_ramp,
        rmu0=constants.rmu0,
    )

    assert phiint == pytest.approx(vscalcparam.expected_phiint)

    assert rlp == pytest.approx(vscalcparam.expected_rlp)

    assert vsbrn == pytest.approx(vscalcparam.expected_vsbrn)

    assert vsind == pytest.approx(vscalcparam.expected_vsind)

    assert vsres == pytest.approx(vscalcparam.expected_vsres)

    assert vsstt == pytest.approx(vscalcparam.expected_vsstt)


class PhyauxParam(NamedTuple):
    tauratio: Any = None

    burnup_in: Any = None

    aspect: Any = None

    dene: Any = None

    deni: Any = None

    dnalp: Any = None

    fusion_rate_density_total: Any = None

    alpha_rate_density_total: Any = None

    plasma_current: Any = None

    sbar: Any = None

    taueff: Any = None

    plasma_volume: Any = None

    expected_burnup: Any = None

    expected_dntau: Any = None

    expected_figmer: Any = None

    expected_fusrat: Any = None

    expected_qfuel: Any = None

    expected_rndfuel: Any = None

    expected_taup: Any = None


@pytest.mark.parametrize(
    "phyauxparam",
    (
        PhyauxParam(
            tauratio=1,
            burnup_in=0,
            aspect=3,
            dene=7.5e19,
            deni=5.8589175702454272e19,
            dnalp=7.5e18,
            fusion_rate_density_total=1.9852091609123786e17,
            alpha_rate_density_total=1.973996644759543e17,
            plasma_current=18398455.678867526,
            sbar=1,
            taueff=3.401323521525641,
            plasma_volume=1888.1711539956691,
            expected_burnup=0.20383432558954095,
            expected_dntau=2.5509926411442307e20,
            expected_figmer=55.195367036602576,
            expected_fusrat=3.7484146722826997e20,
            expected_qfuel=1.8389516394951276e21,
            expected_rndfuel=3.7484146722826997e20,
            expected_taup=37.993985551650177,
        ),
        PhyauxParam(
            tauratio=1,
            burnup_in=0,
            aspect=3,
            dene=7.5e19,
            deni=5.8576156204039725e19,
            dnalp=7.5e18,
            fusion_rate_density_total=1.9843269653375773e17,
            alpha_rate_density_total=1.9731194318497056e17,
            plasma_current=18398455.678867526,
            sbar=1,
            taueff=3.402116961408892,
            plasma_volume=1888.1711539956691,
            expected_burnup=0.20387039462081086,
            expected_dntau=2.5515877210566689e20,
            expected_figmer=55.195367036602576,
            expected_fusrat=3.7467489360461772e20,
            expected_qfuel=1.8378092331723546e21,
            expected_rndfuel=3.7467489360461772e20,
            expected_taup=38.010876984618747,
        ),
    ),
)
def test_phyaux(phyauxparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for phyaux.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param phyauxparam: the data used to mock and assert in this test.
    :type phyauxparam: phyauxparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "tauratio", phyauxparam.tauratio)

    monkeypatch.setattr(physics_variables, "burnup_in", phyauxparam.burnup_in)

    burnup, dntau, figmer, fusrat, qfuel, rndfuel, taup = physics.phyaux(
        aspect=phyauxparam.aspect,
        dene=phyauxparam.dene,
        deni=phyauxparam.deni,
        dnalp=phyauxparam.dnalp,
        fusion_rate_density_total=phyauxparam.fusion_rate_density_total,
        alpha_rate_density_total=phyauxparam.alpha_rate_density_total,
        plasma_current=phyauxparam.plasma_current,
        sbar=phyauxparam.sbar,
        taueff=phyauxparam.taueff,
        plasma_volume=phyauxparam.plasma_volume,
    )

    assert burnup == pytest.approx(phyauxparam.expected_burnup)

    assert dntau == pytest.approx(phyauxparam.expected_dntau)

    assert figmer == pytest.approx(phyauxparam.expected_figmer)

    assert fusrat == pytest.approx(phyauxparam.expected_fusrat)

    assert qfuel == pytest.approx(phyauxparam.expected_qfuel)

    assert rndfuel == pytest.approx(phyauxparam.expected_rndfuel)

    assert taup == pytest.approx(phyauxparam.expected_taup)


def test_rether():
    assert rether(
        1.0, 1.45, 7.5e19, 17.81065204, 12.0, 13.0, 0.43258985
    ) == pytest.approx(0.028360489673597476)


class PohmParam(NamedTuple):
    aspect: Any = None

    plasma_res_factor: Any = None

    inductive_current_fraction: Any = None

    kappa95: Any = None

    plasma_current: Any = None

    rmajor: Any = None

    rminor: Any = None

    ten: Any = None

    plasma_volume: Any = None

    zeff: Any = None

    expected_pohmpv: Any = None

    expected_pohmmw: Any = None

    expected_rpfac: Any = None

    expected_rplas: Any = None


@pytest.mark.parametrize(
    "pohmparam",
    (
        PohmParam(
            aspect=3,
            plasma_res_factor=0.70000000000000007,
            inductive_current_fraction=0.59999999999999998,
            kappa95=1.6517857142857142,
            plasma_current=18398455.678867526,
            rmajor=8,
            rminor=2.6666666666666665,
            ten=12.626131115905864,
            plasma_volume=1888.1711539956691,
            zeff=2.0909945616489103,
            expected_pohmpv=0.0004062519138005805,
            expected_pohmmw=0.7670731448937912,
            expected_rpfac=2.5,
            expected_rplas=3.7767895536275952e-09,
        ),
    ),
)
def test_pohm(pohmparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for pohm.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param pohmparam: the data used to mock and assert in this test.
    :type pohmparam: pohmparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "aspect", pohmparam.aspect)

    monkeypatch.setattr(
        physics_variables, "plasma_res_factor", pohmparam.plasma_res_factor
    )

    pohmpv, pohmmw, rpfac, rplas = physics.pohm(
        inductive_current_fraction=pohmparam.inductive_current_fraction,
        kappa95=pohmparam.kappa95,
        plasma_current=pohmparam.plasma_current,
        rmajor=pohmparam.rmajor,
        rminor=pohmparam.rminor,
        ten=pohmparam.ten,
        plasma_volume=pohmparam.plasma_volume,
        zeff=pohmparam.zeff,
    )

    assert pohmpv == pytest.approx(pohmparam.expected_pohmpv)

    assert pohmmw == pytest.approx(pohmparam.expected_pohmmw)

    assert rpfac == pytest.approx(pohmparam.expected_rpfac)

    assert rplas == pytest.approx(pohmparam.expected_rplas)


class CalculateDensityLimitParam(NamedTuple):
    i_density_limit: Any = None

    bt: Any = None

    pdivt: Any = None

    pinjmw: Any = None

    plasma_current: Any = None

    prn1: Any = None

    q95: Any = None

    qcyl: Any = None

    rmajor: Any = None

    rminor: Any = None

    sarea: Any = None

    zeff: Any = None

    expected_dnelimt: Any = None

    expected_dlimit: Any = None


@pytest.mark.parametrize(
    "calculatedensitylimitparam",
    (
        CalculateDensityLimitParam(
            i_density_limit=7,
            bt=5.1847188735686647,
            pdivt=162.32943903093374,
            pinjmw=79.928763793309031,
            plasma_current=16702766.338258133,
            prn1=0.4614366315228275,
            q95=3.5068029786872268,
            qcyl=3.8769445264202052,
            rmajor=8,
            rminor=2.6666666666666665,
            sarea=1173.8427771245592,
            zeff=2.5668755115791471,
            expected_dnelimt=7.4765470107450917e19,
            expected_dlimit=(
                5.2955542598288974e19,
                1.0934080161360552e20,
                4.3286395478282289e19,
                1.9109162908046821e21,
                4.2410183109151497e20,
                5.0149533075302982e19,
                7.4765470107450917e19,
                8.7406037163890049e20,
            ),
        ),
    ),
)
def test_calculate_density_limit(calculatedensitylimitparam, physics):
    """
    Automatically generated Regression Unit Test for calculate_density_limit().

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param calculatedensitylimitparam: the data used to mock and assert in this test.
    :type calculatedensitylimitparam: calculatedensitylimitparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    dlimit, dnelimt = physics.calculate_density_limit(
        i_density_limit=calculatedensitylimitparam.i_density_limit,
        bt=calculatedensitylimitparam.bt,
        pdivt=calculatedensitylimitparam.pdivt,
        pinjmw=calculatedensitylimitparam.pinjmw,
        plasma_current=calculatedensitylimitparam.plasma_current,
        prn1=calculatedensitylimitparam.prn1,
        q95=calculatedensitylimitparam.q95,
        qcyl=calculatedensitylimitparam.qcyl,
        rmajor=calculatedensitylimitparam.rmajor,
        rminor=calculatedensitylimitparam.rminor,
        sarea=calculatedensitylimitparam.sarea,
        zeff=calculatedensitylimitparam.zeff,
    )

    assert dnelimt == pytest.approx(calculatedensitylimitparam.expected_dnelimt)

    assert dlimit == pytest.approx(calculatedensitylimitparam.expected_dlimit)


class PcondParam(NamedTuple):
    iradloss: Any = None

    tauee_in: Any = None

    pradpv: Any = None

    kappaa_ipb: Any = None

    pohmmw: Any = None

    f_alpha_plasma: Any = None

    iinvqd: Any = None

    isc: Any = None

    ignite: Any = None

    afuel: Any = None

    alpha_power_total: Any = None

    aspect: Any = None

    bt: Any = None

    dene: Any = None

    dnitot: Any = None

    dnla: Any = None

    eps: Any = None

    hfact: Any = None

    kappa: Any = None

    kappa95: Any = None

    non_alpha_charged_power: Any = None

    pinjmw: Any = None

    plasma_current: Any = None

    pcoreradpv: Any = None

    q: Any = None

    qstar: Any = None

    rmajor: Any = None

    rminor: Any = None

    te: Any = None

    ten: Any = None

    tin: Any = None

    plasma_volume: Any = None

    xarea: Any = None

    zeff: Any = None

    expected_kappaa_ipb: Any = None

    expected_kappaa: Any = None

    expected_powerht: Any = None

    expected_ptrepv: Any = None

    expected_ptripv: Any = None

    expected_tauee: Any = None

    expected_taueff: Any = None

    expected_tauei: Any = None


@pytest.mark.parametrize(
    "pcondparam",
    (
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=32,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=6.1946504123280199,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.012570664670798823,
            expected_ptripv=0.011156836223364355,
            expected_tauee=21.17616899712392,
            expected_tauei=21.17616899712392,
            expected_taueff=21.17616899712392,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=33,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.96163409847948844,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.081458458765440875,
            expected_ptripv=0.072296788376262536,
            expected_tauee=3.2679051814806361,
            expected_tauei=3.2679051814806361,
            expected_taueff=3.2679051814806366,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=34,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.1960655150953154,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.081327750398391824,
            expected_ptripv=0.072180780839479111,
            expected_tauee=3.2731572946627923,
            expected_tauei=3.2731572946627923,
            expected_taueff=3.2731572946627923,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=35,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.78453691772934719,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.12077931279593926,
            expected_ptripv=0.10719520783694242,
            expected_tauee=2.2040075681235445,
            expected_tauei=2.2040075681235445,
            expected_taueff=2.2040075681235445,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=36,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.1619079679077555,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.081309182573405442,
            expected_ptripv=0.072164301346324039,
            expected_tauee=3.2739047552801135,
            expected_tauei=3.2739047552801135,
            expected_taueff=3.2739047552801135,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=37,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.7340642483550435,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.08142612910014517,
            expected_ptripv=0.072268094843318462,
            expected_tauee=3.269202679985145,
            expected_tauei=3.269202679985145,
            expected_taueff=3.2692026799851455,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=38,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.1117392853827024,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.072709146314778414,
            expected_ptripv=0.064531515128155068,
            expected_tauee=3.6611421391548524,
            expected_tauei=3.6611421391548524,
            expected_taueff=3.6611421391548529,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=39,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.84477227311274361,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.078529089520063822,
            expected_ptripv=0.069696886639614319,
            expected_tauee=3.3898077909969717,
            expected_tauei=3.3898077909969717,
            expected_taueff=3.3898077909969717,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=40,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.6096667103064701,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.081359821043062386,
            expected_ptripv=0.072209244483967094,
            expected_tauee=3.2718670722507683,
            expected_tauei=3.2718670722507683,
            expected_taueff=3.2718670722507692,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=41,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.67053301699102119,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.081513009259203975,
            expected_ptripv=0.072345203550858064,
            expected_tauee=3.2657182196344126,
            expected_tauei=3.2657182196344126,
            expected_taueff=3.2657182196344126,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=42,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=2.1212580310552207,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.073097992274882811,
            expected_ptripv=0.064876627403966491,
            expected_tauee=3.6416666339340682,
            expected_tauei=3.6416666339340682,
            expected_taueff=3.6416666339340686,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=43,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=50.095480115636271,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.081423406537449478,
            expected_ptripv=0.072265678488503571,
            expected_tauee=3.2693119926464509,
            expected_tauei=3.2693119926464509,
            expected_taueff=3.2693119926464513,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=44,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=87.869318916638761,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.081419881443596701,
            expected_ptripv=0.072262549863580605,
            expected_tauee=3.2694535383156871,
            expected_tauei=3.2694535383156871,
            expected_taueff=3.2694535383156871,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=45,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=28.562137619592285,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.081421142032658531,
            expected_ptripv=0.072263668673609824,
            expected_tauee=3.2694029195542003,
            expected_tauei=3.2694029195542003,
            expected_taueff=3.2694029195542003,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=46,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.50082256309019457,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.079599509500323962,
            expected_ptripv=0.070646915991500636,
            expected_tauee=3.3442231132583498,
            expected_tauei=3.3442231132583498,
            expected_taueff=3.3442231132583502,
            expected_powerht=290.18368660937881,
        ),
        PcondParam(
            iradloss=1,
            tauee_in=0,
            pradpv=0.11824275660100725,
            kappaa_ipb=1.68145080681586,
            pohmmw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            iinvqd=1,
            isc=47,
            ignite=0,
            afuel=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            dnitot=7.1529510234203251e19,
            dnla=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.77961193402355955,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            non_alpha_charged_power=1.2453296074483358,
            pinjmw=75.397788712812741,
            plasma_current=16616203.759182997,
            pcoreradpv=0.047757569353246924,
            q=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            ten=13.745148298980761,
            tin=13.745148298980761,
            plasma_volume=1888.1711539956691,
            xarea=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappaa_ipb=1.68145080681586,
            expected_kappaa=1.7187938085542791,
            expected_ptrepv=0.070108167457759038,
            expected_ptripv=0.062223069561580698,
            expected_tauee=3.7969687288631331,
            expected_tauei=3.7969687288631331,
            expected_taueff=3.7969687288631335,
            expected_powerht=290.18368660937881,
        ),
    ),
)
def test_pcond(pcondparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for pcond.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param pcondparam: the data used to mock and assert in this test.
    :type pcondparam: pcondparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "iradloss", pcondparam.iradloss)

    monkeypatch.setattr(physics_variables, "tauee_in", pcondparam.tauee_in)

    monkeypatch.setattr(physics_variables, "pradpv", pcondparam.pradpv)

    monkeypatch.setattr(physics_variables, "kappaa_ipb", pcondparam.kappaa_ipb)

    monkeypatch.setattr(physics_variables, "pohmmw", pcondparam.pohmmw)

    monkeypatch.setattr(physics_variables, "f_alpha_plasma", pcondparam.f_alpha_plasma)

    kappaa, ptrepv, ptripv, tauee, tauei, taueff, powerht = physics.pcond(
        iinvqd=pcondparam.iinvqd,
        isc=pcondparam.isc,
        ignite=pcondparam.ignite,
        afuel=pcondparam.afuel,
        alpha_power_total=pcondparam.alpha_power_total,
        aspect=pcondparam.aspect,
        bt=pcondparam.bt,
        dene=pcondparam.dene,
        dnitot=pcondparam.dnitot,
        dnla=pcondparam.dnla,
        eps=pcondparam.eps,
        hfact=pcondparam.hfact,
        kappa=pcondparam.kappa,
        kappa95=pcondparam.kappa95,
        non_alpha_charged_power=pcondparam.non_alpha_charged_power,
        pinjmw=pcondparam.pinjmw,
        plasma_current=pcondparam.plasma_current,
        pcoreradpv=pcondparam.pcoreradpv,
        q=pcondparam.q,
        qstar=pcondparam.qstar,
        rmajor=pcondparam.rmajor,
        rminor=pcondparam.rminor,
        te=pcondparam.te,
        ten=pcondparam.ten,
        tin=pcondparam.tin,
        plasma_volume=pcondparam.plasma_volume,
        xarea=pcondparam.xarea,
        zeff=pcondparam.zeff,
    )

    assert physics_variables.kappaa_ipb == pytest.approx(pcondparam.expected_kappaa_ipb)

    assert kappaa == pytest.approx(pcondparam.expected_kappaa)

    assert powerht == pytest.approx(pcondparam.expected_powerht)

    assert ptrepv == pytest.approx(pcondparam.expected_ptrepv)

    assert ptripv == pytest.approx(pcondparam.expected_ptripv)

    assert tauee == pytest.approx(pcondparam.expected_tauee)

    assert taueff == pytest.approx(pcondparam.expected_taueff)

    assert tauei == pytest.approx(pcondparam.expected_tauei)
