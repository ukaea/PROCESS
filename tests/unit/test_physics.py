"""Unit tests for physics.f90."""

from typing import Any, NamedTuple

import numpy as np
import pytest

from process.current_drive import CurrentDrive
from process.fortran import (
    constants,
    current_drive_variables,
    impurity_radiation_module,
    physics_module,
    physics_variables,
)
from process.impurity_radiation import initialise_imprad
from process.physics import (
    Physics,
    calculate_beta_limit,
    calculate_current_coefficient_hastie,
    calculate_plasma_current_peng,
    calculate_poloidal_beta,
    calculate_poloidal_field,
    calculate_volt_second_requirements,
    diamagnetic_fraction_hender,
    diamagnetic_fraction_scene,
    ps_fraction_scene,
    res_diff_time,
    rether,
)
from process.plasma_profiles import PlasmaProfile


@pytest.fixture
def physics():
    """Provides Physics object for testing.

    :returns: initialised Physics object
    :rtype: process.physics.Physics
    """
    return Physics(PlasmaProfile(), CurrentDrive(PlasmaProfile()))


def test_calculate_poloidal_beta():
    """Test calculate_poloidal_beta()"""
    beta_poloidal = calculate_poloidal_beta(5.347, 0.852, 0.0307)
    assert beta_poloidal == pytest.approx(1.209, abs=0.001)


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

    vol_plasma: Any = None

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
            vol_plasma=1888.1711539956691,
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
        vol_plasma=bootstrapfractioniter89param.vol_plasma,
    )

    assert bootstrap_current_fraction == pytest.approx(
        bootstrapfractioniter89param.expected_bootipf
    )


class BootstrapFractionNevinsParam(NamedTuple):
    te0: Any = None

    ne0: Any = None

    alphan: Any = None

    beta_toroidal: Any = None

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
            beta_toroidal=0.03,
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
        beta_toroidal=bootstrapfractionnevinsparam.beta_toroidal,
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
    nd_ions_total: Any = None

    rminor: Any = None

    tesep: Any = None

    ti: Any = None

    triang: Any = None

    q0: Any = None

    m_fuel_amu: Any = None

    zeff: Any = None

    rhopedn: Any = None

    bt: Any = None

    plasma_current: Any = None

    a_plasma_poloidal: Any = None

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
            nd_ions_total=7.1297522422781575e19,
            rminor=2.6666666666666665,
            tesep=0.10000000000000001,
            ti=12.570861186498382,
            triang=0.5,
            q0=1,
            m_fuel_amu=2.5,
            zeff=2.5211399464385624,
            rhopedn=0.9400000000000001,
            bt=5.326133750416047,
            plasma_current=16528278.760008096,
            a_plasma_poloidal=38.39822223637151,
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
        physics_variables, "nd_ions_total", bootstrapfractionsauterparam.nd_ions_total
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

    monkeypatch.setattr(
        physics_variables, "m_fuel_amu", bootstrapfractionsauterparam.m_fuel_amu
    )

    monkeypatch.setattr(physics_variables, "zeff", bootstrapfractionsauterparam.zeff)

    monkeypatch.setattr(
        physics_variables, "rhopedn", bootstrapfractionsauterparam.rhopedn
    )

    monkeypatch.setattr(physics_variables, "bt", bootstrapfractionsauterparam.bt)

    monkeypatch.setattr(
        physics_variables, "plasma_current", bootstrapfractionsauterparam.plasma_current
    )

    monkeypatch.setattr(
        physics_variables,
        "a_plasma_poloidal",
        bootstrapfractionsauterparam.a_plasma_poloidal,
    )

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
    beta_poloidal: Any = None

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
            beta_poloidal=1.3184383457774960,
            q95=3.5151046634673557,
            q0=1.0,
            alphan=1.0,
            alphat=1.45,
            eps=1 / 3,
            rli=1.2098126022585098,
            expected_bfs=0.3501274900057279,
        ),
        BootstrapFractionSakaiParam(
            beta_poloidal=1.1701245502231756,
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

    monkeypatch.setattr(
        physics_variables, "beta_poloidal", bootstrapfractionsakaiparam.beta_poloidal
    )

    monkeypatch.setattr(physics_variables, "q95", bootstrapfractionsakaiparam.q95)

    monkeypatch.setattr(physics_variables, "q0", bootstrapfractionsakaiparam.q0)

    monkeypatch.setattr(physics_variables, "alphan", bootstrapfractionsakaiparam.alphan)

    monkeypatch.setattr(physics_variables, "alphat", bootstrapfractionsakaiparam.alphat)

    monkeypatch.setattr(physics_variables, "eps", bootstrapfractionsakaiparam.eps)

    monkeypatch.setattr(physics_variables, "rli", bootstrapfractionsakaiparam.rli)

    bfs = physics.bootstrap_fraction_sakai(
        beta_poloidal=bootstrapfractionsakaiparam.beta_poloidal,
        q95=bootstrapfractionsakaiparam.q95,
        q0=bootstrapfractionsakaiparam.q0,
        alphan=bootstrapfractionsakaiparam.alphan,
        alphat=bootstrapfractionsakaiparam.alphat,
        eps=bootstrapfractionsakaiparam.eps,
        rli=bootstrapfractionsakaiparam.rli,
    )

    assert bfs == pytest.approx(bootstrapfractionsakaiparam.expected_bfs)


class BootstrapFractionAriesParam(NamedTuple):
    beta_poloidal: Any = None

    rli: Any = None

    core_density: Any = None

    average_density: Any = None

    inverse_aspect: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionariesparam",
    (
        BootstrapFractionAriesParam(
            beta_poloidal=1.2708883332338736,
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
        beta_poloidal=bootstrapfractionariesparam.beta_poloidal,
        rli=bootstrapfractionariesparam.rli,
        core_density=bootstrapfractionariesparam.core_density,
        average_density=bootstrapfractionariesparam.average_density,
        inverse_aspect=bootstrapfractionariesparam.inverse_aspect,
    )

    assert bfs == pytest.approx(bootstrapfractionariesparam.expected_bfs)


class BootstrapFractionAndradeParam(NamedTuple):
    beta_poloidal: Any = None

    core_pressure: Any = None

    average_pressure: Any = None

    inverse_aspect: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionandradeparam",
    (
        BootstrapFractionAndradeParam(
            beta_poloidal=1.2708883332338736,
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
        beta_poloidal=bootstrapfractionandradeparam.beta_poloidal,
        core_pressure=bootstrapfractionandradeparam.core_pressure,
        average_pressure=bootstrapfractionandradeparam.average_pressure,
        inverse_aspect=bootstrapfractionandradeparam.inverse_aspect,
    )

    assert bfs == pytest.approx(bootstrapfractionandradeparam.expected_bfs)


class BootstrapFractionHoangParam(NamedTuple):
    beta_poloidal: Any = None

    pressure_index: Any = None

    current_index: Any = None

    inverse_aspect: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionhoangparam",
    (
        BootstrapFractionHoangParam(
            beta_poloidal=1.2708883332338736,
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
        beta_poloidal=bootstrapfractionhoangparam.beta_poloidal,
        pressure_index=bootstrapfractionhoangparam.pressure_index,
        current_index=bootstrapfractionhoangparam.current_index,
        inverse_aspect=bootstrapfractionhoangparam.inverse_aspect,
    )

    assert bfs == pytest.approx(bootstrapfractionhoangparam.expected_bfs)


class BootstrapFractionWongParam(NamedTuple):
    beta_poloidal: Any = None

    density_index: Any = None

    temperature_index: Any = None

    inverse_aspect: Any = None

    elongation: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionwongparam",
    (
        BootstrapFractionWongParam(
            beta_poloidal=1.2708883332338736,
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
        beta_poloidal=bootstrapfractionwongparam.beta_poloidal,
        density_index=bootstrapfractionwongparam.density_index,
        temperature_index=bootstrapfractionwongparam.temperature_index,
        inverse_aspect=bootstrapfractionwongparam.inverse_aspect,
        elongation=bootstrapfractionwongparam.elongation,
    )

    assert bfs == pytest.approx(bootstrapfractionwongparam.expected_bfs)


class BootstrapFractionGiIParam(NamedTuple):
    beta_poloidal: Any = None

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
            beta_poloidal=1.2708883332338736,
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
def test_bootstrap_fraction_gi_I(bootstrapfractiongiiparam, physics):  # noqa: N802
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_gi.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam
    """

    bfs = physics.bootstrap_fraction_gi_I(
        beta_poloidal=bootstrapfractiongiiparam.beta_poloidal,
        pressure_index=bootstrapfractiongiiparam.pressure_index,
        temperature_index=bootstrapfractiongiiparam.temperature_index,
        inverse_aspect=bootstrapfractiongiiparam.inverse_aspect,
        effective_charge=bootstrapfractiongiiparam.effective_charge,
        q95=bootstrapfractiongiiparam.q95,
        q0=bootstrapfractiongiiparam.q0,
    )

    assert bfs == pytest.approx(bootstrapfractiongiiparam.expected_bfs)


class BootstrapFractionGiIIParam(NamedTuple):
    beta_poloidal: Any = None

    pressure_index: Any = None

    temperature_index: Any = None

    inverse_aspect: Any = None

    effective_charge: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractiongiiiparam",
    (
        BootstrapFractionGiIIParam(
            beta_poloidal=1.2708883332338736,
            pressure_index=2.4500000000000002e00,
            temperature_index=1.4500000000000000e00,
            inverse_aspect=1 / 3,
            effective_charge=2.5368733516769737e00,
            expected_bfs=8.8502865710180589e-01,
        ),
    ),
)
def test_bootstrap_fraction_gi_II(bootstrapfractiongiiiparam, physics):  # noqa: N802
    """
    Automatically generated Regression Unit Test for bootstrap_fraction_gi.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.

    :param bootstrapfractionsauterparam: the data used to mock and assert in this test.
    :type bootstrapfractionsauterparam: bootstrapfractionsauterparam
    """

    bfs = physics.bootstrap_fraction_gi_II(
        beta_poloidal=bootstrapfractiongiiiparam.beta_poloidal,
        pressure_index=bootstrapfractiongiiiparam.pressure_index,
        temperature_index=bootstrapfractiongiiiparam.temperature_index,
        inverse_aspect=bootstrapfractiongiiiparam.inverse_aspect,
        effective_charge=bootstrapfractiongiiiparam.effective_charge,
    )

    assert bfs == pytest.approx(bootstrapfractiongiiiparam.expected_bfs)


class PlasmaCurrentParam(NamedTuple):
    beta_norm_total: Any = None

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

    len_plasma_poloidal: Any = None

    q0: Any = None

    q95: Any = None

    rmajor: Any = None

    rminor: Any = None

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
            beta_norm_total=0,
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
            len_plasma_poloidal=24.081367139525412,
            q0=1,
            q95=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
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
            beta_norm_total=2.4784688886891844,
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
            len_plasma_poloidal=24.081367139525412,
            q0=1,
            q95=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
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
        "beta_norm_total",
        plasmacurrentparam.beta_norm_total,
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
        len_plasma_poloidal=plasmacurrentparam.len_plasma_poloidal,
        q0=plasmacurrentparam.q0,
        q95=plasmacurrentparam.q95,
        rmajor=plasmacurrentparam.rmajor,
        rminor=plasmacurrentparam.rminor,
        triang=plasmacurrentparam.triang,
        triang95=plasmacurrentparam.triang95,
    )

    assert physics_variables.beta_norm_total == pytest.approx(
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


def test_calculate_beta_limit():
    assert calculate_beta_limit(12, 4.879, 18300000, 2.5) == pytest.approx(0.0297619)


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

    m_fuel_amu: Any = None

    f_tritium: Any = None

    nd_fuel_ions: Any = None

    m_ions_total_amu: Any = None

    nd_ions_total: Any = None

    f_nd_protium_electrons: Any = None

    zeffai: Any = None

    rncne: Any = None

    rnone: Any = None

    f_alpha_ion: Any = None

    f_nd_alpha_electron: Any = None

    dlamee: Any = None

    f_nd_beam_electron: Any = None

    zeff: Any = None

    nd_impurities: Any = None

    pcoef: Any = None

    alpha_rate_density_total: Any = None

    rnfene: Any = None

    m_beam_amu: Any = None

    te: Any = None

    proton_rate_density: Any = None

    f_deuterium: Any = None

    alphan: Any = None

    nd_beam_ions: Any = None

    f_helium3: Any = None

    nd_alphas: Any = None

    dene: Any = None

    nd_protons: Any = None

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

    expected_m_fuel_amu: Any = None

    expected_nd_fuel_ions: Any = None

    expected_m_ions_total_amu: Any = None

    expected_nd_ions_total: Any = None

    expected_zeffai: Any = None

    expected_f_alpha_ion: Any = None

    expected_dlamee: Any = None

    expected_zeff: Any = None

    expected_nd_impurities: Any = None

    expected_m_beam_amu: Any = None

    expected_nd_alphas: Any = None

    expected_nd_protons: Any = None

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
            m_fuel_amu=0,
            f_tritium=0.5,
            nd_fuel_ions=0,
            m_ions_total_amu=0,
            nd_ions_total=0,
            f_nd_protium_electrons=0,
            zeffai=0,
            rncne=0,
            rnone=0,
            f_alpha_ion=0,
            f_nd_alpha_electron=0.10000000000000001,
            dlamee=0,
            f_nd_beam_electron=0,
            zeff=0,
            nd_impurities=0,
            pcoef=0,
            alpha_rate_density_total=0,
            rnfene=0,
            m_beam_amu=0,
            te=12,
            proton_rate_density=0,
            f_deuterium=0.5,
            alphan=1,
            nd_beam_ions=0,
            f_helium3=0,
            nd_alphas=0,
            dene=7.5e19,
            nd_protons=0,
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
            expected_m_fuel_amu=2.5145269632339478,
            expected_nd_fuel_ions=5.8589175702454272e19,
            expected_m_ions_total_amu=2.7395439636787726,
            expected_nd_ions_total=6.6125550702454276e19,
            expected_zeffai=0.43046641789338563,
            expected_f_alpha_ion=0.3154069116809366,
            expected_dlamee=17.510652035055571,
            expected_zeff=2.0909945616489103,
            expected_nd_impurities=28875000000000004,
            expected_m_beam_amu=2.01355414,
            expected_nd_alphas=7.5e18,
            expected_nd_protons=7500000000000000,
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
            impurity_arr_z=np.array(
                np.array((1, 2, 4, 6, 7, 8, 10, 14, 18, 26, 28, 36, 54, 74), order="F"),
                order="F",
            ).transpose(),
            impurity_arr_amass=np.array(
                np.array(
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
            m_fuel_amu=2.5,
            f_tritium=0.5,
            nd_fuel_ions=5.8589175702454272e19,
            m_ions_total_amu=2.7395439636787726,
            nd_ions_total=6.6125550702454276e19,
            f_nd_protium_electrons=0,
            zeffai=0.43046641789338563,
            rncne=0,
            rnone=0,
            f_alpha_ion=0.3154069116809366,
            f_nd_alpha_electron=0.10000000000000001,
            dlamee=17.510652035055571,
            f_nd_beam_electron=0,
            zeff=2.0909945616489103,
            nd_impurities=28875000000000004,
            pcoef=1.0521775929921553,
            alpha_rate_density_total=1.973996644759543e17,
            rnfene=0,
            m_beam_amu=2.01355414,
            te=12,
            proton_rate_density=540072280299564.38,
            f_deuterium=0.5,
            alphan=1,
            nd_beam_ions=0,
            f_helium3=0,
            nd_alphas=7.5e18,
            dene=7.5e19,
            nd_protons=7500000000000000,
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
            expected_m_fuel_amu=2.5145269632339478,
            expected_nd_fuel_ions=5.8576156204039725e19,
            expected_m_ions_total_amu=2.739245767577763,
            expected_nd_ions_total=6.6125550702454276e19,
            expected_zeffai=0.43056686748101997,
            expected_f_alpha_ion=0.26903878212105858,
            expected_dlamee=17.510652035055571,
            expected_zeff=2.0909945616489103,
            expected_nd_impurities=28875000000000004,
            expected_m_beam_amu=2.01355414,
            expected_nd_alphas=7.5e18,
            expected_nd_protons=20519498414548412,
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

    monkeypatch.setattr(
        physics_variables, "m_fuel_amu", plasmacompositionparam.m_fuel_amu
    )

    monkeypatch.setattr(
        physics_variables, "f_tritium", plasmacompositionparam.f_tritium
    )

    monkeypatch.setattr(
        physics_variables, "nd_fuel_ions", plasmacompositionparam.nd_fuel_ions
    )

    monkeypatch.setattr(
        physics_variables, "m_ions_total_amu", plasmacompositionparam.m_ions_total_amu
    )

    monkeypatch.setattr(
        physics_variables, "nd_ions_total", plasmacompositionparam.nd_ions_total
    )

    monkeypatch.setattr(
        physics_variables,
        "f_nd_protium_electrons",
        plasmacompositionparam.f_nd_protium_electrons,
    )

    monkeypatch.setattr(physics_variables, "zeffai", plasmacompositionparam.zeffai)

    monkeypatch.setattr(physics_variables, "rncne", plasmacompositionparam.rncne)

    monkeypatch.setattr(physics_variables, "rnone", plasmacompositionparam.rnone)

    monkeypatch.setattr(
        physics_variables, "f_alpha_ion", plasmacompositionparam.f_alpha_ion
    )

    monkeypatch.setattr(
        physics_variables,
        "f_nd_alpha_electron",
        plasmacompositionparam.f_nd_alpha_electron,
    )

    monkeypatch.setattr(physics_variables, "dlamee", plasmacompositionparam.dlamee)

    monkeypatch.setattr(
        physics_variables,
        "f_nd_beam_electron",
        plasmacompositionparam.f_nd_beam_electron,
    )

    monkeypatch.setattr(physics_variables, "zeff", plasmacompositionparam.zeff)

    monkeypatch.setattr(
        physics_variables, "nd_impurities", plasmacompositionparam.nd_impurities
    )

    monkeypatch.setattr(physics_variables, "pcoef", plasmacompositionparam.pcoef)

    monkeypatch.setattr(
        physics_variables,
        "alpha_rate_density_total",
        plasmacompositionparam.alpha_rate_density_total,
    )

    monkeypatch.setattr(physics_variables, "rnfene", plasmacompositionparam.rnfene)

    monkeypatch.setattr(
        physics_variables, "m_beam_amu", plasmacompositionparam.m_beam_amu
    )

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

    monkeypatch.setattr(
        physics_variables, "nd_beam_ions", plasmacompositionparam.nd_beam_ions
    )

    monkeypatch.setattr(
        physics_variables, "f_helium3", plasmacompositionparam.f_helium3
    )

    monkeypatch.setattr(
        physics_variables, "nd_alphas", plasmacompositionparam.nd_alphas
    )

    monkeypatch.setattr(physics_variables, "dene", plasmacompositionparam.dene)

    monkeypatch.setattr(
        physics_variables, "nd_protons", plasmacompositionparam.nd_protons
    )

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

    assert physics_variables.m_fuel_amu == pytest.approx(
        plasmacompositionparam.expected_m_fuel_amu
    )

    assert physics_variables.nd_fuel_ions == pytest.approx(
        plasmacompositionparam.expected_nd_fuel_ions
    )

    assert physics_variables.m_ions_total_amu == pytest.approx(
        plasmacompositionparam.expected_m_ions_total_amu
    )

    assert physics_variables.nd_ions_total == pytest.approx(
        plasmacompositionparam.expected_nd_ions_total
    )

    assert physics_variables.zeffai == pytest.approx(
        plasmacompositionparam.expected_zeffai
    )

    assert physics_variables.f_alpha_ion == pytest.approx(
        plasmacompositionparam.expected_f_alpha_ion
    )

    assert physics_variables.zeff == pytest.approx(plasmacompositionparam.expected_zeff)

    assert physics_variables.nd_impurities == pytest.approx(
        plasmacompositionparam.expected_nd_impurities
    )

    assert physics_variables.m_beam_amu == pytest.approx(
        plasmacompositionparam.expected_m_beam_amu
    )

    assert physics_variables.nd_alphas == pytest.approx(
        plasmacompositionparam.expected_nd_alphas
    )

    assert physics_variables.nd_protons == pytest.approx(
        plasmacompositionparam.expected_nd_protons
    )

    assert physics_module.first_call == pytest.approx(
        plasmacompositionparam.expected_first_call
    )


class VoltSecondReqParam(NamedTuple):
    csawth: Any = None

    eps: Any = None

    inductive_current_fraction: Any = None

    ejima_coeff: Any = None

    kappa: Any = None

    plasma_current: Any = None

    rli: Any = None

    rmajor: Any = None

    res_plasma: Any = None

    t_burn: Any = None

    t_fusion_ramp: Any = None

    expected_vs_plasma_internal: Any = None

    expected_rlp: Any = None

    expected_vsbrn: Any = None

    expected_vsind: Any = None

    expected_vsres: Any = None

    expected_vsstt: Any = None


@pytest.mark.parametrize(
    "voltsecondreqparam",
    (
        VoltSecondReqParam(
            csawth=1,
            eps=0.33333333333333331,
            inductive_current_fraction=0.59999999999999998,
            ejima_coeff=0.30000000000000004,
            kappa=1.8500000000000001,
            plasma_current=18398455.678867526,
            rli=1.2064840230894305,
            rmajor=8,
            res_plasma=3.7767895536275952e-09,
            t_burn=1000,
            t_fusion_ramp=10,
            expected_vs_plasma_internal=111.57651734747576,
            expected_rlp=1.4075705307248088e-05,
            expected_vsbrn=42.109179697761263,
            expected_vsind=258.97124024420435,
            expected_vsres=55.488435095110333,
            expected_vsstt=356.56885503707593,
        ),
        VoltSecondReqParam(
            csawth=1,
            eps=0.33333333333333331,
            inductive_current_fraction=0.59999999999999998,
            ejima_coeff=0.30000000000000004,
            kappa=1.8500000000000001,
            plasma_current=18398455.678867526,
            rli=1.2064840230894305,
            rmajor=8,
            res_plasma=3.7767895536275952e-09,
            t_burn=0,
            t_fusion_ramp=10,
            expected_vs_plasma_internal=111.57651734747576,
            expected_rlp=1.4075705307248088e-05,
            expected_vsbrn=0.41692257126496302,
            expected_vsind=258.97124024420435,
            expected_vsres=55.488435095110333,
            expected_vsstt=314.87659791057968,
        ),
    ),
)
def test_vscalc(voltsecondreqparam):
    """
    Automatically generated Regression Unit Test for calculate_volt_second_requirements().

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param voltsecondreqparam: the data used to mock and assert in this test.
    :type voltsecondreqparam: voltsecondreqparam
    """

    vs_plasma_internal, rlp, vsbrn, vsind, vsres, vsstt = (
        calculate_volt_second_requirements(
            csawth=voltsecondreqparam.csawth,
            eps=voltsecondreqparam.eps,
            inductive_current_fraction=voltsecondreqparam.inductive_current_fraction,
            ejima_coeff=voltsecondreqparam.ejima_coeff,
            kappa=voltsecondreqparam.kappa,
            plasma_current=voltsecondreqparam.plasma_current,
            rli=voltsecondreqparam.rli,
            rmajor=voltsecondreqparam.rmajor,
            res_plasma=voltsecondreqparam.res_plasma,
            t_burn=voltsecondreqparam.t_burn,
            t_fusion_ramp=voltsecondreqparam.t_fusion_ramp,
        )
    )

    assert vs_plasma_internal == pytest.approx(
        voltsecondreqparam.expected_vs_plasma_internal
    )

    assert rlp == pytest.approx(voltsecondreqparam.expected_rlp)

    assert vsbrn == pytest.approx(voltsecondreqparam.expected_vsbrn)

    assert vsind == pytest.approx(voltsecondreqparam.expected_vsind)

    assert vsres == pytest.approx(voltsecondreqparam.expected_vsres)

    assert vsstt == pytest.approx(voltsecondreqparam.expected_vsstt)


class PhyauxParam(NamedTuple):
    tauratio: Any = None

    burnup_in: Any = None

    aspect: Any = None

    dene: Any = None

    te: Any = None

    nd_fuel_ions: Any = None

    nd_alphas: Any = None

    fusion_rate_density_total: Any = None

    alpha_rate_density_total: Any = None

    plasma_current: Any = None

    sbar: Any = None

    t_energy_confinement: Any = None

    vol_plasma: Any = None

    expected_burnup: Any = None

    expected_ntau: Any = None

    expected_nTtau: Any = None

    expected_figmer: Any = None

    expected_fusrat: Any = None

    expected_qfuel: Any = None

    expected_rndfuel: Any = None

    expected_t_alpha_confinement: Any = None


@pytest.mark.parametrize(
    "phyauxparam",
    (
        PhyauxParam(
            tauratio=1,
            burnup_in=0,
            aspect=3,
            dene=7.5e19,
            te=12.569,
            nd_fuel_ions=5.8589175702454272e19,
            nd_alphas=7.5e18,
            fusion_rate_density_total=1.9852091609123786e17,
            alpha_rate_density_total=1.973996644759543e17,
            plasma_current=18398455.678867526,
            sbar=1,
            t_energy_confinement=3.401323521525641,
            vol_plasma=1888.1711539956691,
            expected_burnup=0.20383432558954095,
            expected_ntau=2.5509926411442307e20,
            expected_nTtau=3.253e21,
            expected_figmer=55.195367036602576,
            expected_fusrat=3.7484146722826997e20,
            expected_qfuel=1.8389516394951276e21,
            expected_rndfuel=3.7484146722826997e20,
            expected_t_alpha_confinement=37.993985551650177,
        ),
        PhyauxParam(
            tauratio=1,
            burnup_in=0,
            aspect=3,
            dene=7.5e19,
            te=12.569,
            nd_fuel_ions=5.8576156204039725e19,
            nd_alphas=7.5e18,
            fusion_rate_density_total=1.9843269653375773e17,
            alpha_rate_density_total=1.9731194318497056e17,
            plasma_current=18398455.678867526,
            sbar=1,
            t_energy_confinement=3.402116961408892,
            vol_plasma=1888.1711539956691,
            expected_burnup=0.20387039462081086,
            expected_ntau=2.5515877210566689e20,
            expected_nTtau=3.253e21,
            expected_figmer=55.195367036602576,
            expected_fusrat=3.7467489360461772e20,
            expected_qfuel=1.8378092331723546e21,
            expected_rndfuel=3.7467489360461772e20,
            expected_t_alpha_confinement=38.010876984618747,
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

    burnup, ntau, nTtau, figmer, fusrat, qfuel, rndfuel, t_alpha_confinement, _ = (
        physics.phyaux(
            aspect=phyauxparam.aspect,
            dene=phyauxparam.dene,
            te=phyauxparam.te,
            nd_fuel_ions=phyauxparam.nd_fuel_ions,
            nd_alphas=phyauxparam.nd_alphas,
            fusion_rate_density_total=phyauxparam.fusion_rate_density_total,
            alpha_rate_density_total=phyauxparam.alpha_rate_density_total,
            plasma_current=phyauxparam.plasma_current,
            sbar=phyauxparam.sbar,
            t_energy_confinement=phyauxparam.t_energy_confinement,
            vol_plasma=phyauxparam.vol_plasma,
        )
    )

    assert burnup == pytest.approx(phyauxparam.expected_burnup)

    assert ntau == pytest.approx(phyauxparam.expected_ntau)

    assert figmer == pytest.approx(phyauxparam.expected_figmer)

    assert fusrat == pytest.approx(phyauxparam.expected_fusrat)

    assert qfuel == pytest.approx(phyauxparam.expected_qfuel)

    assert rndfuel == pytest.approx(phyauxparam.expected_rndfuel)

    assert t_alpha_confinement == pytest.approx(
        phyauxparam.expected_t_alpha_confinement
    )


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

    vol_plasma: Any = None

    zeff: Any = None

    expected_pden_plasma_ohmic_mw: Any = None

    expected_p_plasma_ohmic_mw: Any = None

    expected_rpfac: Any = None

    expected_res_plasma: Any = None


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
            vol_plasma=1888.1711539956691,
            zeff=2.0909945616489103,
            expected_pden_plasma_ohmic_mw=0.0004062519138005805,
            expected_p_plasma_ohmic_mw=0.7670731448937912,
            expected_rpfac=2.5,
            expected_res_plasma=3.7767895536275952e-09,
        ),
    ),
)
def test_pohm(pohmparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for plasma_ohmic_heating.

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

    (
        pden_plasma_ohmic_mw,
        p_plasma_ohmic_mw,
        rpfac,
        res_plasma,
    ) = physics.plasma_ohmic_heating(
        inductive_current_fraction=pohmparam.inductive_current_fraction,
        kappa95=pohmparam.kappa95,
        plasma_current=pohmparam.plasma_current,
        rmajor=pohmparam.rmajor,
        rminor=pohmparam.rminor,
        ten=pohmparam.ten,
        vol_plasma=pohmparam.vol_plasma,
        zeff=pohmparam.zeff,
    )

    assert pden_plasma_ohmic_mw == pytest.approx(
        pohmparam.expected_pden_plasma_ohmic_mw
    )

    assert p_plasma_ohmic_mw == pytest.approx(pohmparam.expected_p_plasma_ohmic_mw)

    assert rpfac == pytest.approx(pohmparam.expected_rpfac)

    assert res_plasma == pytest.approx(pohmparam.expected_res_plasma)


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

    a_plasma_surface: Any = None

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
            a_plasma_surface=1173.8427771245592,
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
        a_plasma_surface=calculatedensitylimitparam.a_plasma_surface,
        zeff=calculatedensitylimitparam.zeff,
    )

    assert dnelimt == pytest.approx(calculatedensitylimitparam.expected_dnelimt)

    assert dlimit == pytest.approx(calculatedensitylimitparam.expected_dlimit)


class ConfinementTimeParam(NamedTuple):
    i_rad_loss: Any = None

    tauee_in: Any = None

    pden_plasma_rad_mw: Any = None

    kappa_ipb: Any = None

    p_plasma_ohmic_mw: Any = None

    f_alpha_plasma: Any = None

    i_confinement_time: Any = None

    ignite: Any = None

    m_fuel_amu: Any = None

    alpha_power_total: Any = None

    aspect: Any = None

    bt: Any = None

    dene: Any = None

    nd_ions_total: Any = None

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

    vol_plasma: Any = None

    a_plasma_poloidal: Any = None

    zeff: Any = None

    expected_kappa_ipb: Any = None

    expected_kappaa: Any = None

    expected_p_plasma_loss_mw: Any = None

    expected_pden_electron_transport_loss_mw: Any = None

    expected_pden_ion_transport_loss_mw: Any = None

    expected_tauee: Any = None

    expected_t_energy_confinement: Any = None

    expected_t_ion_energy_confinement: Any = None


@pytest.mark.parametrize(
    "confinementtimeparam",
    (
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=32,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.012572050692511346,
            expected_pden_ion_transport_loss_mw=0.011158066358576262,
            expected_tauee=21.17616899712392,
            expected_t_ion_energy_confinement=21.17616899712392,
            expected_t_energy_confinement=21.17616899712392,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=33,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.08146744024696746,
            expected_pden_ion_transport_loss_mw=0.07230475970642361,
            expected_tauee=3.2679051814806361,
            expected_t_ion_energy_confinement=3.2679051814806361,
            expected_t_energy_confinement=3.2679051814806366,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=34,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.0813367174682195,
            expected_pden_ion_transport_loss_mw=0.07218873937883169,
            expected_tauee=3.2731572946627923,
            expected_t_ion_energy_confinement=3.2731572946627923,
            expected_t_energy_confinement=3.2731572946627923,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=35,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.12079262973297819,
            expected_pden_ion_transport_loss_mw=0.10720702701193681,
            expected_tauee=2.2040075681235445,
            expected_t_ion_energy_confinement=2.2040075681235445,
            expected_t_energy_confinement=2.2040075681235445,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=36,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.08131814759597392,
            expected_pden_ion_transport_loss_mw=0.07217225806867361,
            expected_tauee=3.2739047552801135,
            expected_t_ion_energy_confinement=3.2739047552801135,
            expected_t_energy_confinement=3.2739047552801135,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=37,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.08143510701705373,
            expected_pden_ion_transport_loss_mw=0.07227606300977574,
            expected_tauee=3.269202679985145,
            expected_t_ion_energy_confinement=3.269202679985145,
            expected_t_energy_confinement=3.2692026799851455,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=38,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.07271716311087716,
            expected_pden_ion_transport_loss_mw=0.06453863027150285,
            expected_tauee=3.6611421391548524,
            expected_t_ion_energy_confinement=3.6611421391548524,
            expected_t_energy_confinement=3.6611421391548529,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=39,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.07853774801394538,
            expected_pden_ion_transport_loss_mw=0.06970457130869961,
            expected_tauee=3.3898077909969717,
            expected_t_ion_energy_confinement=3.3898077909969717,
            expected_t_energy_confinement=3.3898077909969717,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=40,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.0840021318362596,
            expected_pden_ion_transport_loss_mw=0.07455437336481374,
            expected_tauee=3.169298972363837,
            expected_t_ion_energy_confinement=3.169298972363837,
            expected_t_energy_confinement=3.169298972363837,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=41,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.08311313602000579,
            expected_pden_ion_transport_loss_mw=0.07376536331761714,
            expected_tauee=3.203198469625145,
            expected_t_ion_energy_confinement=3.203198469625145,
            expected_t_energy_confinement=3.203198469625145,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=42,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.07310605194452542,
            expected_pden_ion_transport_loss_mw=0.0648837805988509,
            expected_tauee=3.6416666339340682,
            expected_t_ion_energy_confinement=3.6416666339340682,
            expected_t_energy_confinement=3.6416666339340686,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=43,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.08143238415417252,
            expected_pden_ion_transport_loss_mw=0.07227364638853734,
            expected_tauee=3.2693119926464509,
            expected_t_ion_energy_confinement=3.2693119926464509,
            expected_t_energy_confinement=3.2693119926464513,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=44,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.08142885867164847,
            expected_pden_ion_transport_loss_mw=0.07227051741865713,
            expected_tauee=3.2694535383156871,
            expected_t_ion_energy_confinement=3.2694535383156871,
            expected_t_energy_confinement=3.2694535383156871,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=45,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.08143011939694184,
            expected_pden_ion_transport_loss_mw=0.07227163634959588,
            expected_tauee=3.2694029195542003,
            expected_t_ion_energy_confinement=3.2694029195542003,
            expected_t_energy_confinement=3.2694029195542003,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=46,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.07960828601702878,
            expected_pden_ion_transport_loss_mw=0.07065470540932789,
            expected_tauee=3.3442231132583498,
            expected_t_ion_energy_confinement=3.3442231132583498,
            expected_t_energy_confinement=3.3442231132583502,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
        ConfinementTimeParam(
            i_rad_loss=1,
            tauee_in=0,
            pden_plasma_rad_mw=0.11824275660100725,
            kappa_ipb=1.68145080681586,
            p_plasma_ohmic_mw=0.63634001890069991,
            f_alpha_plasma=0.94999999999999996,
            i_confinement_time=47,
            ignite=0,
            m_fuel_amu=2.5,
            alpha_power_total=319.03020327154269,
            aspect=3,
            bt=5.2375830857646202,
            dene=8.0593948787884524e19,
            nd_ions_total=7.1529510234203251e19,
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
            vol_plasma=1888.1711539956691,
            a_plasma_poloidal=38.39822223637151,
            zeff=2.4987360098030775,
            expected_kappa_ipb=1.68145080681586,
            expected_pden_electron_transport_loss_mw=0.07148441348179191,
            expected_pden_ion_transport_loss_mw=0.06344452856118785,
            expected_tauee=3.7242785823911264,
            expected_t_ion_energy_confinement=3.7242785823911264,
            expected_t_energy_confinement=3.7242785823911264,
            expected_p_plasma_loss_mw=290.18368660937881,
        ),
    ),
)
def test_calculate_confinement_time(confinementtimeparam, monkeypatch, physics):
    """
    Automatically generated Regression Unit Test for calculate_confinement_time().

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param confinementtimeparam: the data used to mock and assert in this test.
    :type confinementtimeparam: confinementtimeparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        physics_variables, "i_rad_loss", confinementtimeparam.i_rad_loss
    )

    monkeypatch.setattr(physics_variables, "tauee_in", confinementtimeparam.tauee_in)

    monkeypatch.setattr(
        physics_variables, "pden_plasma_rad_mw", confinementtimeparam.pden_plasma_rad_mw
    )

    monkeypatch.setattr(physics_variables, "kappa_ipb", confinementtimeparam.kappa_ipb)

    monkeypatch.setattr(
        physics_variables, "p_plasma_ohmic_mw", confinementtimeparam.p_plasma_ohmic_mw
    )

    monkeypatch.setattr(
        physics_variables, "f_alpha_plasma", confinementtimeparam.f_alpha_plasma
    )

    (
        pden_electron_transport_loss_mw,
        pden_ion_transport_loss_mw,
        t_electron_energy_confinement,
        t_ion_energy_confinement,
        t_energy_confinement,
        p_plasma_loss_mw,
    ) = physics.calculate_confinement_time(
        i_confinement_time=confinementtimeparam.i_confinement_time,
        ignite=confinementtimeparam.ignite,
        m_fuel_amu=confinementtimeparam.m_fuel_amu,
        alpha_power_total=confinementtimeparam.alpha_power_total,
        aspect=confinementtimeparam.aspect,
        bt=confinementtimeparam.bt,
        dene=confinementtimeparam.dene,
        nd_ions_total=confinementtimeparam.nd_ions_total,
        dnla=confinementtimeparam.dnla,
        eps=confinementtimeparam.eps,
        hfact=confinementtimeparam.hfact,
        kappa=confinementtimeparam.kappa,
        kappa95=confinementtimeparam.kappa95,
        non_alpha_charged_power=confinementtimeparam.non_alpha_charged_power,
        pinjmw=confinementtimeparam.pinjmw,
        plasma_current=confinementtimeparam.plasma_current,
        pcoreradpv=confinementtimeparam.pcoreradpv,
        q=confinementtimeparam.q,
        qstar=confinementtimeparam.qstar,
        rmajor=confinementtimeparam.rmajor,
        rminor=confinementtimeparam.rminor,
        ten=confinementtimeparam.ten,
        tin=confinementtimeparam.tin,
        vol_plasma=confinementtimeparam.vol_plasma,
        zeff=confinementtimeparam.zeff,
    )

    assert physics_variables.kappa_ipb == pytest.approx(
        confinementtimeparam.expected_kappa_ipb
    )

    assert p_plasma_loss_mw == pytest.approx(
        confinementtimeparam.expected_p_plasma_loss_mw
    )

    assert pden_electron_transport_loss_mw == pytest.approx(
        confinementtimeparam.expected_pden_electron_transport_loss_mw
    )

    assert pden_ion_transport_loss_mw == pytest.approx(
        confinementtimeparam.expected_pden_ion_transport_loss_mw
    )

    assert t_electron_energy_confinement == pytest.approx(
        confinementtimeparam.expected_tauee
    )

    assert t_energy_confinement == pytest.approx(
        confinementtimeparam.expected_t_energy_confinement
    )

    assert t_ion_energy_confinement == pytest.approx(
        confinementtimeparam.expected_t_ion_energy_confinement
    )
