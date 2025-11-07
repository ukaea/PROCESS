"""Unit tests for physics.f90."""

from typing import Any, NamedTuple

import numpy as np
import pytest

from process import constants
from process.current_drive import (
    CurrentDrive,
    ElectronBernstein,
    ElectronCyclotron,
    IonCyclotron,
    LowerHybrid,
    NeutralBeam,
)
from process.data_structure import (
    current_drive_variables,
    impurity_radiation_module,
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
    return Physics(
        PlasmaProfile(),
        CurrentDrive(
            PlasmaProfile(),
            electron_cyclotron=ElectronCyclotron(plasma_profile=PlasmaProfile()),
            ion_cyclotron=IonCyclotron(plasma_profile=PlasmaProfile()),
            neutral_beam=NeutralBeam(plasma_profile=PlasmaProfile()),
            electron_bernstein=ElectronBernstein(plasma_profile=PlasmaProfile()),
            lower_hybrid=LowerHybrid(plasma_profile=PlasmaProfile()),
        ),
    )


def test_calculate_poloidal_beta():
    """Test calculate_poloidal_beta()"""
    beta_poloidal = calculate_poloidal_beta(5.347, 0.852, 0.0307)
    assert beta_poloidal == pytest.approx(1.209, abs=0.001)


def test_res_diff_time():
    """Test res_diff_time()"""
    t_plasma_res_diffusion = res_diff_time(9.137, 2.909e-9, 1.65)
    assert t_plasma_res_diffusion == pytest.approx(4784.3, abs=0.1)


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

    b_plasma_toroidal_on_axis: Any = None

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
            b_plasma_toroidal_on_axis=5.7802910787445487,
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

    f_c_plasma_bootstrap = physics.bootstrap_fraction_iter89(
        aspect=bootstrapfractioniter89param.aspect,
        beta=bootstrapfractioniter89param.beta,
        b_plasma_toroidal_on_axis=bootstrapfractioniter89param.b_plasma_toroidal_on_axis,
        plasma_current=bootstrapfractioniter89param.plasma_current,
        q95=bootstrapfractioniter89param.q95,
        q0=bootstrapfractioniter89param.q0,
        rmajor=bootstrapfractioniter89param.rmajor,
        vol_plasma=bootstrapfractioniter89param.vol_plasma,
    )

    assert f_c_plasma_bootstrap == pytest.approx(
        bootstrapfractioniter89param.expected_bootipf
    )


class BootstrapFractionNevinsParam(NamedTuple):
    temp_plasma_electron_on_axis_kev: Any = None

    nd_plasma_electron_on_axis: Any = None

    alphan: Any = None

    beta_toroidal: Any = None

    b_plasma_toroidal_on_axis: Any = None

    nd_plasma_electrons_vol_avg: Any = None

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
            temp_plasma_electron_on_axis_kev=24.402321098330372,
            nd_plasma_electron_on_axis=8.515060981068918e19,
            alphan=1.0,
            beta_toroidal=0.03,
            b_plasma_toroidal_on_axis=5.7,
            nd_plasma_electrons_vol_avg=18398455.678867526,
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

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_electron_on_axis_kev",
        bootstrapfractionnevinsparam.temp_plasma_electron_on_axis_kev,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_electron_on_axis",
        bootstrapfractionnevinsparam.nd_plasma_electron_on_axis,
    )

    fibs = physics.bootstrap_fraction_nevins(
        alphan=bootstrapfractionnevinsparam.alphan,
        alphat=bootstrapfractionnevinsparam.alphat,
        beta_toroidal=bootstrapfractionnevinsparam.beta_toroidal,
        b_plasma_toroidal_on_axis=bootstrapfractionnevinsparam.b_plasma_toroidal_on_axis,
        nd_plasma_electrons_vol_avg=bootstrapfractionnevinsparam.nd_plasma_electrons_vol_avg,
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
    nd_plasma_ions_total_vol_avg: Any = None

    rminor: Any = None

    temp_plasma_separatrix_kev: Any = None

    temp_plasma_ion_vol_avg_kev: Any = None

    triang: Any = None

    q0: Any = None

    m_fuel_amu: Any = None

    zeff: Any = None

    radius_plasma_pedestal_density_norm: Any = None

    b_plasma_toroidal_on_axis: Any = None

    plasma_current: Any = None

    a_plasma_poloidal: Any = None

    f_plasma_fuel_helium3: Any = None

    temp_plasma_pedestal_kev: Any = None

    nd_plasma_electrons_vol_avg: Any = None

    te: Any = None

    rmajor: Any = None

    q95: Any = None

    nd_plasma_separatrix_electron: Any = None

    temp_plasma_electron_on_axis_kev: Any = None

    nd_plasma_pedestal_electron: Any = None

    tbeta: Any = None

    nd_plasma_electron_on_axis: Any = None

    alphan: Any = None

    radius_plasma_pedestal_temp_norm: Any = None

    alphat: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionsauterparam",
    (
        BootstrapFractionSauterParam(
            nd_plasma_ions_total_vol_avg=7.1297522422781575e19,
            rminor=2.6666666666666665,
            temp_plasma_separatrix_kev=0.10000000000000001,
            temp_plasma_ion_vol_avg_kev=12.570861186498382,
            triang=0.5,
            q0=1,
            m_fuel_amu=2.5,
            zeff=2.5211399464385624,
            radius_plasma_pedestal_density_norm=0.9400000000000001,
            b_plasma_toroidal_on_axis=5.326133750416047,
            plasma_current=16528278.760008096,
            a_plasma_poloidal=38.39822223637151,
            f_plasma_fuel_helium3=0,
            temp_plasma_pedestal_kev=5.5,
            nd_plasma_electrons_vol_avg=8.016748468651018e19,
            te=12.570861186498382,
            rmajor=8,
            q95=3.5,
            nd_plasma_separatrix_electron=3.6992211545476006e19,
            temp_plasma_electron_on_axis_kev=25.986118047669795,
            nd_plasma_pedestal_electron=6.2886759627309195e19,
            tbeta=2,
            nd_plasma_electron_on_axis=1.054474759840606e20,
            alphan=1,
            radius_plasma_pedestal_temp_norm=0.9400000000000001,
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
        physics_variables,
        "nd_plasma_ions_total_vol_avg",
        bootstrapfractionsauterparam.nd_plasma_ions_total_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables, "rminor", bootstrapfractionsauterparam.rminor
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_separatrix_kev",
        bootstrapfractionsauterparam.temp_plasma_separatrix_kev,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_ion_vol_avg_kev",
        bootstrapfractionsauterparam.temp_plasma_ion_vol_avg_kev,
    )

    monkeypatch.setattr(
        physics_variables, "triang", bootstrapfractionsauterparam.triang
    )

    monkeypatch.setattr(physics_variables, "q0", bootstrapfractionsauterparam.q0)

    monkeypatch.setattr(
        physics_variables, "m_fuel_amu", bootstrapfractionsauterparam.m_fuel_amu
    )

    monkeypatch.setattr(
        physics_variables,
        "n_charge_plasma_effective_vol_avg",
        bootstrapfractionsauterparam.zeff,
    )

    monkeypatch.setattr(
        physics_variables,
        "radius_plasma_pedestal_density_norm",
        bootstrapfractionsauterparam.radius_plasma_pedestal_density_norm,
    )

    monkeypatch.setattr(
        physics_variables,
        "b_plasma_toroidal_on_axis",
        bootstrapfractionsauterparam.b_plasma_toroidal_on_axis,
    )

    monkeypatch.setattr(
        physics_variables, "plasma_current", bootstrapfractionsauterparam.plasma_current
    )

    monkeypatch.setattr(
        physics_variables,
        "a_plasma_poloidal",
        bootstrapfractionsauterparam.a_plasma_poloidal,
    )

    monkeypatch.setattr(
        physics_variables,
        "f_plasma_fuel_helium3",
        bootstrapfractionsauterparam.f_plasma_fuel_helium3,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_pedestal_kev",
        bootstrapfractionsauterparam.temp_plasma_pedestal_kev,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_electrons_vol_avg",
        bootstrapfractionsauterparam.nd_plasma_electrons_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_electron_vol_avg_kev",
        bootstrapfractionsauterparam.te,
    )

    monkeypatch.setattr(
        physics_variables, "rmajor", bootstrapfractionsauterparam.rmajor
    )

    monkeypatch.setattr(physics_variables, "q95", bootstrapfractionsauterparam.q95)

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_separatrix_electron",
        bootstrapfractionsauterparam.nd_plasma_separatrix_electron,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_electron_on_axis_kev",
        bootstrapfractionsauterparam.temp_plasma_electron_on_axis_kev,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_pedestal_electron",
        bootstrapfractionsauterparam.nd_plasma_pedestal_electron,
    )

    monkeypatch.setattr(physics_variables, "tbeta", bootstrapfractionsauterparam.tbeta)

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_electron_on_axis",
        bootstrapfractionsauterparam.nd_plasma_electron_on_axis,
    )

    monkeypatch.setattr(
        physics_variables, "alphan", bootstrapfractionsauterparam.alphan
    )

    monkeypatch.setattr(
        physics_variables,
        "radius_plasma_pedestal_temp_norm",
        bootstrapfractionsauterparam.radius_plasma_pedestal_temp_norm,
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

    ind_plasma_internal_norm: Any = None

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
            ind_plasma_internal_norm=1.2098126022585098,
            expected_bfs=0.3501274900057279,
        ),
        BootstrapFractionSakaiParam(
            beta_poloidal=1.1701245502231756,
            q95=5.1746754543339177,
            q0=2.0,
            alphan=0.9,
            alphat=1.3999999999999999,
            eps=1 / 1.8,
            ind_plasma_internal_norm=0.3,
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
        physics_variables,
        "beta_poloidal_vol_avg",
        bootstrapfractionsakaiparam.beta_poloidal,
    )

    monkeypatch.setattr(physics_variables, "q95", bootstrapfractionsakaiparam.q95)

    monkeypatch.setattr(physics_variables, "q0", bootstrapfractionsakaiparam.q0)

    monkeypatch.setattr(physics_variables, "alphan", bootstrapfractionsakaiparam.alphan)

    monkeypatch.setattr(physics_variables, "alphat", bootstrapfractionsakaiparam.alphat)

    monkeypatch.setattr(physics_variables, "eps", bootstrapfractionsakaiparam.eps)

    monkeypatch.setattr(
        physics_variables,
        "ind_plasma_internal_norm",
        bootstrapfractionsakaiparam.ind_plasma_internal_norm,
    )

    bfs = physics.bootstrap_fraction_sakai(
        beta_poloidal=bootstrapfractionsakaiparam.beta_poloidal,
        q95=bootstrapfractionsakaiparam.q95,
        q0=bootstrapfractionsakaiparam.q0,
        alphan=bootstrapfractionsakaiparam.alphan,
        alphat=bootstrapfractionsakaiparam.alphat,
        eps=bootstrapfractionsakaiparam.eps,
        ind_plasma_internal_norm=bootstrapfractionsakaiparam.ind_plasma_internal_norm,
    )

    assert bfs == pytest.approx(bootstrapfractionsakaiparam.expected_bfs)


class BootstrapFractionAriesParam(NamedTuple):
    beta_poloidal: Any = None

    ind_plasma_internal_norm: Any = None

    core_density: Any = None

    average_density: Any = None

    inverse_aspect: Any = None

    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionariesparam",
    (
        BootstrapFractionAriesParam(
            beta_poloidal=1.2708883332338736,
            ind_plasma_internal_norm=1.4279108047138775,
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
        ind_plasma_internal_norm=bootstrapfractionariesparam.ind_plasma_internal_norm,
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


class BootstrapFractionSugiyamaLModeParam(NamedTuple):
    eps: Any = None
    beta_poloidal: Any = None
    alphan: Any = None
    alphat: Any = None
    zeff: Any = None
    q95: Any = None
    q0: Any = None
    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionsugiyamalparam",
    (
        BootstrapFractionSugiyamaLModeParam(
            eps=0.33333333,
            beta_poloidal=1.2708883332338736,
            alphan=1.0,
            alphat=1.45,
            zeff=2.5,
            q95=3.5,
            q0=1.0,
            expected_bfs=0.5700433347072669,
        ),
        BootstrapFractionSugiyamaLModeParam(
            eps=0.25,
            beta_poloidal=1.1,
            alphan=0.9,
            alphat=1.3,
            zeff=2.0,
            q95=4.0,
            q0=1.2,
            expected_bfs=0.42806341892858024,
        ),
    ),
)
def test_bootstrap_fraction_sugiyama_l_mode(bootstrapfractionsugiyamalparam, physics):
    """
    Test bootstrap_fraction_sugiyama_l_mode function.

    This test validates the Sugiyama L-mode bootstrap fraction calculation
    using predefined parameters and expected results.

    :param bootstrapfractionsugiyamalparam: Parameters for the test case.
    :type bootstrapfractionsugiyamalparam: BootstrapFractionSugiyamaLModeParam
    """
    bfs = physics.bootstrap_fraction_sugiyama_l_mode(
        eps=bootstrapfractionsugiyamalparam.eps,
        beta_poloidal=bootstrapfractionsugiyamalparam.beta_poloidal,
        alphan=bootstrapfractionsugiyamalparam.alphan,
        alphat=bootstrapfractionsugiyamalparam.alphat,
        zeff=bootstrapfractionsugiyamalparam.zeff,
        q95=bootstrapfractionsugiyamalparam.q95,
        q0=bootstrapfractionsugiyamalparam.q0,
    )

    assert bfs == pytest.approx(bootstrapfractionsugiyamalparam.expected_bfs)


class BootstrapFractionSugiyamaHModeParam(NamedTuple):
    eps: Any = None
    beta_poloidal: Any = None
    alphan: Any = None
    alphat: Any = None
    tbeta: Any = None
    zeff: Any = None
    q95: Any = None
    q0: Any = None
    radius_plasma_pedestal_density_norm: Any = None
    nd_plasma_pedestal_electron: Any = None
    n_greenwald: Any = None
    temp_plasma_pedestal_kev: Any = None
    expected_bfs: Any = None


@pytest.mark.parametrize(
    "bootstrapfractionsugiyamahparam",
    (
        BootstrapFractionSugiyamaHModeParam(
            eps=0.33333333,
            beta_poloidal=1.2708883332338736,
            alphan=1.0,
            alphat=1.45,
            tbeta=2.0,
            zeff=2.5,
            q95=3.5,
            q0=1.0,
            radius_plasma_pedestal_density_norm=0.9,
            nd_plasma_pedestal_electron=6.0e19,
            n_greenwald=8.0e19,
            temp_plasma_pedestal_kev=5.0,
            expected_bfs=0.5875359328840783,
        ),
        BootstrapFractionSugiyamaHModeParam(
            eps=0.25,
            beta_poloidal=1.1,
            alphan=0.9,
            alphat=1.3,
            tbeta=1.8,
            zeff=2.0,
            q95=4.0,
            q0=1.2,
            radius_plasma_pedestal_density_norm=0.85,
            nd_plasma_pedestal_electron=5.5e19,
            n_greenwald=7.5e19,
            temp_plasma_pedestal_kev=4.5,
            expected_bfs=0.40154857221044604,
        ),
        # JA-DEMO steady state case from the paper
        BootstrapFractionSugiyamaHModeParam(
            eps=0.2847058823529412,
            beta_poloidal=1.52,
            alphan=1.7,
            alphat=4.5,
            tbeta=2.65,
            zeff=1.84,
            q95=4.09,
            q0=1.0,
            radius_plasma_pedestal_density_norm=0.91,
            nd_plasma_pedestal_electron=0.98e20,
            n_greenwald=1e20,
            temp_plasma_pedestal_kev=6.0,
            expected_bfs=0.5634482876932788,
        ),
        # ITER 15MA case from the paper
        BootstrapFractionSugiyamaHModeParam(
            eps=0.3225806451612903,
            beta_poloidal=0.691,
            alphan=0.3,
            alphat=4.0,
            tbeta=2.0,
            zeff=1.8,
            q95=6.26,
            q0=1.0,
            radius_plasma_pedestal_density_norm=0.93,
            nd_plasma_pedestal_electron=0.75e20,
            n_greenwald=1e20,
            temp_plasma_pedestal_kev=6.0,
            expected_bfs=0.2770187998673241,
        ),
    ),
)
def test_bootstrap_fraction_sugiyama_h_mode(bootstrapfractionsugiyamahparam, physics):
    """
    Test bootstrap_fraction_sugiyama_h_mode function.

    This test validates the Sugiyama H-mode bootstrap fraction calculation
    using predefined parameters and expected results.

    :param bootstrapfractionsugiyamahparam: Parameters for the test case.
    :type bootstrapfractionsugiyamahparam: BootstrapFractionSugiyamaHModeParam
    """
    bfs = physics.bootstrap_fraction_sugiyama_h_mode(
        eps=bootstrapfractionsugiyamahparam.eps,
        beta_poloidal=bootstrapfractionsugiyamahparam.beta_poloidal,
        alphan=bootstrapfractionsugiyamahparam.alphan,
        alphat=bootstrapfractionsugiyamahparam.alphat,
        tbeta=bootstrapfractionsugiyamahparam.tbeta,
        zeff=bootstrapfractionsugiyamahparam.zeff,
        q95=bootstrapfractionsugiyamahparam.q95,
        q0=bootstrapfractionsugiyamahparam.q0,
        radius_plasma_pedestal_density_norm=bootstrapfractionsugiyamahparam.radius_plasma_pedestal_density_norm,
        nd_plasma_pedestal_electron=bootstrapfractionsugiyamahparam.nd_plasma_pedestal_electron,
        n_greenwald=bootstrapfractionsugiyamahparam.n_greenwald,
        temp_plasma_pedestal_kev=bootstrapfractionsugiyamahparam.temp_plasma_pedestal_kev,
    )

    assert bfs == pytest.approx(bootstrapfractionsugiyamahparam.expected_bfs)


class PlasmaCurrentParam(NamedTuple):
    beta_norm_total: Any = None

    beta: Any = None

    i_plasma_current: Any = None

    alphaj: Any = None

    alphap: Any = None

    b_plasma_toroidal_on_axis: Any = None

    eps: Any = None

    kappa: Any = None

    kappa95: Any = None

    pres_plasma_on_axis: Any = None

    len_plasma_poloidal: Any = None

    q95: Any = None

    rmajor: Any = None

    rminor: Any = None

    triang: Any = None

    triang95: Any = None

    expected_normalised_total_beta: Any = None

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
            alphaj=1,
            alphap=0,
            b_plasma_toroidal_on_axis=5.7000000000000002,
            eps=0.33333333333333331,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            pres_plasma_on_axis=0,
            len_plasma_poloidal=24.081367139525412,
            q95=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            triang95=0.33333333333333331,
            expected_normalised_total_beta=2.4784688886891844,
            expected_bp=0.96008591022564971,
            expected_qstar=2.900802902105021,
            expected_plasma_current=18398455.678867526,
        ),
        PlasmaCurrentParam(
            beta_norm_total=2.4784688886891844,
            beta=0.030000000000000006,
            i_plasma_current=4,
            alphaj=1.9008029008029004,
            alphap=2.4500000000000002,
            b_plasma_toroidal_on_axis=5.7000000000000002,
            eps=0.33333333333333331,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            pres_plasma_on_axis=626431.90482713911,
            len_plasma_poloidal=24.081367139525412,
            q95=3.5,
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            triang95=0.33333333333333331,
            expected_normalised_total_beta=2.4784688886891844,
            expected_bp=0.96008591022564971,
            expected_qstar=2.900802902105021,
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

    monkeypatch.setattr(
        physics_variables, "beta_total_vol_avg", plasmacurrentparam.beta
    )

    b_plasma_poloidal_average, qstar, plasma_current = physics.calculate_plasma_current(
        i_plasma_current=plasmacurrentparam.i_plasma_current,
        alphaj=plasmacurrentparam.alphaj,
        alphap=plasmacurrentparam.alphap,
        b_plasma_toroidal_on_axis=plasmacurrentparam.b_plasma_toroidal_on_axis,
        eps=plasmacurrentparam.eps,
        kappa=plasmacurrentparam.kappa,
        kappa95=plasmacurrentparam.kappa95,
        pres_plasma_on_axis=plasmacurrentparam.pres_plasma_on_axis,
        len_plasma_poloidal=plasmacurrentparam.len_plasma_poloidal,
        q95=plasmacurrentparam.q95,
        rmajor=plasmacurrentparam.rmajor,
        rminor=plasmacurrentparam.rminor,
        triang=plasmacurrentparam.triang,
        triang95=plasmacurrentparam.triang95,
    )

    assert physics_variables.beta_norm_total == pytest.approx(
        plasmacurrentparam.expected_normalised_total_beta
    )

    assert b_plasma_poloidal_average == pytest.approx(plasmacurrentparam.expected_bp)

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
                "b_plasma_toroidal_on_axis": 12,
                "kappa": 1.85,
                "delta": 0.5,
            },
            37.08247253247967,
        ),
        (
            {
                "q95": 2.5,
                "aspect": 3.0,
                "eps": 0.33333333,
                "rminor": 1.5,
                "b_plasma_toroidal_on_axis": 12,
                "kappa": 1.85,
                "delta": 0.5,
            },
            31.594671010223617,
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
                "b_plasma_toroidal_on_axis": 12,
                "kappa": 1.85,
                "delta": 0.5,
                "perim": 24,
                "rmu0": constants.RMU0,
            },
            3.4401302177092803,
        ),
        (
            {
                "i_plasma_current": 2,
                "ip": 1.6e7,
                "q95": 2.5,
                "aspect": 3.0,
                "eps": 0.33333333,
                "b_plasma_toroidal_on_axis": 12,
                "kappa": 1.85,
                "delta": 0.5,
                "perim": 24,
                "rmu0": constants.RMU0,
            },
            2.9310284627233205,
        ),
        (
            {
                "i_plasma_current": 3,
                "ip": 1.6e7,
                "q95": 2.5,
                "aspect": 3.0,
                "eps": 0.33333333,
                "b_plasma_toroidal_on_axis": 12,
                "kappa": 1.85,
                "delta": 0.5,
                "perim": 24,
                "rmu0": constants.RMU0,
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
        5, 5, 12, 0.5, 0.33, 1.85, 2e3, constants.RMU0
    ) == pytest.approx(2.518876726889116)


class PlasmaCompositionParam(NamedTuple):
    f_beam_tritium: Any = None

    f_nd_impurity_electron_array: Any = None

    impurity_arr_z: Any = None

    m_impurity_amu_array: Any = None

    alphat: Any = None

    i_plasma_ignited: Any = None

    f_alpha_electron: Any = None

    m_fuel_amu: Any = None

    f_plasma_fuel_tritium: Any = None

    nd_plasma_fuel_ions_vol_avg: Any = None

    m_ions_total_amu: Any = None

    nd_plasma_ions_total_vol_avg: Any = None

    f_nd_protium_electrons: Any = None

    zeffai: Any = None

    f_nd_plasma_carbon_electron: Any = None

    f_nd_plasma_oxygen_electron: Any = None

    f_alpha_ion: Any = None

    f_nd_alpha_electron: Any = None

    dlamee: Any = None

    f_nd_beam_electron: Any = None

    zeff: Any = None

    nd_plasma_impurities_vol_avg: Any = None

    f_temp_plasma_electron_density_vol_avg: Any = None

    fusden_alpha_total: Any = None

    f_nd_plasma_iron_argon_electron: Any = None

    m_beam_amu: Any = None

    te: Any = None

    proton_rate_density: Any = None

    f_plasma_fuel_deuterium: Any = None

    alphan: Any = None

    nd_beam_ions: Any = None

    f_plasma_fuel_helium3: Any = None

    nd_plasma_alphas_vol_avg: Any = None

    nd_plasma_electrons_vol_avg: Any = None

    nd_plasma_protons_vol_avg: Any = None

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
            f_beam_tritium=9.9999999999999995e-07,
            f_nd_impurity_electron_array=np.array([
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
            ]),
            impurity_arr_z=np.array([
                1.0,
                2.0,
                4.0,
                6.0,
                7.0,
                8.0,
                10.0,
                14.0,
                18.0,
                26.0,
                28.0,
                36.0,
                54.0,
                74.0,
            ]),
            m_impurity_amu_array=np.array([
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
            ]),
            alphat=1.45,
            i_plasma_ignited=0,
            f_alpha_electron=0,
            m_fuel_amu=0,
            f_plasma_fuel_tritium=0.5,
            nd_plasma_fuel_ions_vol_avg=0,
            m_ions_total_amu=0,
            nd_plasma_ions_total_vol_avg=0,
            f_nd_protium_electrons=0,
            zeffai=0,
            f_nd_plasma_carbon_electron=0,
            f_nd_plasma_oxygen_electron=0,
            f_alpha_ion=0,
            f_nd_alpha_electron=0.10000000000000001,
            dlamee=0,
            f_nd_beam_electron=0,
            zeff=0,
            nd_plasma_impurities_vol_avg=0,
            f_temp_plasma_electron_density_vol_avg=0,
            fusden_alpha_total=0,
            f_nd_plasma_iron_argon_electron=0,
            m_beam_amu=0,
            te=12,
            proton_rate_density=0,
            f_plasma_fuel_deuterium=0.5,
            alphan=1,
            nd_beam_ions=0,
            f_plasma_fuel_helium3=0,
            nd_plasma_alphas_vol_avg=0,
            nd_plasma_electrons_vol_avg=7.5e19,
            nd_plasma_protons_vol_avg=0,
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
            expected_impurity_arr_frac=np.array([
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
            ]),
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
            f_beam_tritium=9.9999999999999995e-07,
            f_nd_impurity_electron_array=np.array([
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
            ]),
            impurity_arr_z=np.array(
                np.array((1, 2, 4, 6, 7, 8, 10, 14, 18, 26, 28, 36, 54, 74), order="F"),
                order="F",
            ).transpose(),
            m_impurity_amu_array=np.array(
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
            i_plasma_ignited=0,
            f_alpha_electron=0.6845930883190634,
            m_fuel_amu=2.5,
            f_plasma_fuel_tritium=0.5,
            nd_plasma_fuel_ions_vol_avg=5.8589175702454272e19,
            m_ions_total_amu=2.7395439636787726,
            nd_plasma_ions_total_vol_avg=6.6125550702454276e19,
            f_nd_protium_electrons=0,
            zeffai=0.43046641789338563,
            f_nd_plasma_carbon_electron=0,
            f_nd_plasma_oxygen_electron=0,
            f_alpha_ion=0.3154069116809366,
            f_nd_alpha_electron=0.10000000000000001,
            dlamee=17.510652035055571,
            f_nd_beam_electron=0,
            zeff=2.0909945616489103,
            nd_plasma_impurities_vol_avg=28875000000000004,
            f_temp_plasma_electron_density_vol_avg=1.0521775929921553,
            fusden_alpha_total=1.973996644759543e17,
            f_nd_plasma_iron_argon_electron=0,
            m_beam_amu=2.01355414,
            te=12,
            proton_rate_density=540072280299564.38,
            f_plasma_fuel_deuterium=0.5,
            alphan=1,
            nd_beam_ions=0,
            f_plasma_fuel_helium3=0,
            nd_plasma_alphas_vol_avg=7.5e18,
            nd_plasma_electrons_vol_avg=7.5e19,
            nd_plasma_protons_vol_avg=7500000000000000,
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
            expected_impurity_arr_frac=np.array([
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
            ]),
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
    initialise_imprad()

    monkeypatch.setattr(
        current_drive_variables, "f_beam_tritium", plasmacompositionparam.f_beam_tritium
    )

    monkeypatch.setattr(
        impurity_radiation_module,
        "f_nd_impurity_electron_array",
        plasmacompositionparam.f_nd_impurity_electron_array,
    )

    monkeypatch.setattr(
        impurity_radiation_module,
        "impurity_arr_z",
        plasmacompositionparam.impurity_arr_z,
    )

    monkeypatch.setattr(
        impurity_radiation_module,
        "m_impurity_amu_array",
        plasmacompositionparam.m_impurity_amu_array,
    )

    monkeypatch.setattr(physics_variables, "alphat", plasmacompositionparam.alphat)

    monkeypatch.setattr(
        physics_variables, "i_plasma_ignited", plasmacompositionparam.i_plasma_ignited
    )

    monkeypatch.setattr(
        physics_variables, "f_alpha_electron", plasmacompositionparam.f_alpha_electron
    )

    monkeypatch.setattr(
        physics_variables, "m_fuel_amu", plasmacompositionparam.m_fuel_amu
    )

    monkeypatch.setattr(
        physics_variables,
        "f_plasma_fuel_tritium",
        plasmacompositionparam.f_plasma_fuel_tritium,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_fuel_ions_vol_avg",
        plasmacompositionparam.nd_plasma_fuel_ions_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables, "m_ions_total_amu", plasmacompositionparam.m_ions_total_amu
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_ions_total_vol_avg",
        plasmacompositionparam.nd_plasma_ions_total_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables,
        "f_nd_protium_electrons",
        plasmacompositionparam.f_nd_protium_electrons,
    )

    monkeypatch.setattr(physics_variables, "zeffai", plasmacompositionparam.zeffai)

    monkeypatch.setattr(
        physics_variables,
        "f_nd_plasma_carbon_electron",
        plasmacompositionparam.f_nd_plasma_carbon_electron,
    )

    monkeypatch.setattr(
        physics_variables,
        "f_nd_plasma_oxygen_electron",
        plasmacompositionparam.f_nd_plasma_oxygen_electron,
    )

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

    monkeypatch.setattr(
        physics_variables,
        "n_charge_plasma_effective_vol_avg",
        plasmacompositionparam.zeff,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_impurities_vol_avg",
        plasmacompositionparam.nd_plasma_impurities_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables,
        "f_temp_plasma_electron_density_vol_avg",
        plasmacompositionparam.f_temp_plasma_electron_density_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables,
        "fusden_alpha_total",
        plasmacompositionparam.fusden_alpha_total,
    )

    monkeypatch.setattr(
        physics_variables,
        "f_nd_plasma_iron_argon_electron",
        plasmacompositionparam.f_nd_plasma_iron_argon_electron,
    )

    monkeypatch.setattr(
        physics_variables, "m_beam_amu", plasmacompositionparam.m_beam_amu
    )

    monkeypatch.setattr(
        physics_variables, "temp_plasma_electron_vol_avg_kev", plasmacompositionparam.te
    )

    monkeypatch.setattr(
        physics_variables,
        "proton_rate_density",
        plasmacompositionparam.proton_rate_density,
    )

    monkeypatch.setattr(
        physics_variables,
        "f_plasma_fuel_deuterium",
        plasmacompositionparam.f_plasma_fuel_deuterium,
    )

    monkeypatch.setattr(physics_variables, "alphan", plasmacompositionparam.alphan)

    monkeypatch.setattr(
        physics_variables, "nd_beam_ions", plasmacompositionparam.nd_beam_ions
    )

    monkeypatch.setattr(
        physics_variables,
        "f_plasma_fuel_helium3",
        plasmacompositionparam.f_plasma_fuel_helium3,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_alphas_vol_avg",
        plasmacompositionparam.nd_plasma_alphas_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_electrons_vol_avg",
        plasmacompositionparam.nd_plasma_electrons_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_protons_vol_avg",
        plasmacompositionparam.nd_plasma_protons_vol_avg,
    )

    monkeypatch.setattr(physics_variables, "iscz", plasmacompositionparam.iscz)

    monkeypatch.setattr(physics_variables, "err242", plasmacompositionparam.err242)

    monkeypatch.setattr(physics_variables, "err243", plasmacompositionparam.err243)

    monkeypatch.setattr(physics_variables, "ptarmw", plasmacompositionparam.ptarmw)

    monkeypatch.setattr(physics_variables, "lambdaio", plasmacompositionparam.lambdaio)

    monkeypatch.setattr(physics_variables, "drsep", plasmacompositionparam.drsep)

    monkeypatch.setattr(physics_variables, "fio", plasmacompositionparam.fio)

    monkeypatch.setattr(physics_variables, "rho_star", plasmacompositionparam.rho_star)

    monkeypatch.setattr(physics_variables, "nu_star", plasmacompositionparam.nu_star)

    monkeypatch.setattr(
        physics_variables, "beta_mcdonald", plasmacompositionparam.beta_mcdonald
    )

    monkeypatch.setattr(physics_variables, "itart_r", plasmacompositionparam.itart_r)

    monkeypatch.setattr(
        physics_variables, "first_call", plasmacompositionparam.first_call
    )

    physics.plasma_composition()

    assert impurity_radiation_module.f_nd_impurity_electron_array == pytest.approx(
        plasmacompositionparam.expected_impurity_arr_frac
    )

    assert physics_variables.f_alpha_electron == pytest.approx(
        plasmacompositionparam.expected_f_alpha_electron
    )

    assert physics_variables.m_fuel_amu == pytest.approx(
        plasmacompositionparam.expected_m_fuel_amu
    )

    assert physics_variables.nd_plasma_fuel_ions_vol_avg == pytest.approx(
        plasmacompositionparam.expected_nd_fuel_ions
    )

    assert physics_variables.m_ions_total_amu == pytest.approx(
        plasmacompositionparam.expected_m_ions_total_amu
    )

    assert physics_variables.nd_plasma_ions_total_vol_avg == pytest.approx(
        plasmacompositionparam.expected_nd_ions_total
    )

    assert physics_variables.zeffai == pytest.approx(
        plasmacompositionparam.expected_zeffai
    )

    assert physics_variables.f_alpha_ion == pytest.approx(
        plasmacompositionparam.expected_f_alpha_ion
    )

    assert physics_variables.n_charge_plasma_effective_vol_avg == pytest.approx(
        plasmacompositionparam.expected_zeff
    )

    assert physics_variables.nd_plasma_impurities_vol_avg == pytest.approx(
        plasmacompositionparam.expected_nd_impurities
    )

    assert physics_variables.m_beam_amu == pytest.approx(
        plasmacompositionparam.expected_m_beam_amu
    )

    assert physics_variables.nd_plasma_alphas_vol_avg == pytest.approx(
        plasmacompositionparam.expected_nd_alphas
    )

    assert physics_variables.nd_plasma_protons_vol_avg == pytest.approx(
        plasmacompositionparam.expected_nd_protons
    )

    assert physics_variables.first_call == pytest.approx(
        plasmacompositionparam.expected_first_call
    )


class VoltSecondReqParam(NamedTuple):
    csawth: Any = None

    eps: Any = None

    f_c_plasma_inductive: Any = None

    ejima_coeff: Any = None

    kappa: Any = None

    plasma_current: Any = None

    ind_plasma_internal_norm: Any = None

    rmajor: Any = None

    res_plasma: Any = None

    t_plant_pulse_burn: Any = None

    t_plant_pulse_fusion_ramp: Any = None

    expected_vs_plasma_internal: Any = None

    expected_ind_plasma: Any = None

    expected_vs_plasma_burn_required: Any = None

    expected_vs_plasma_ramp_required: Any = None

    expected_vs_plasma_ind_ramp: Any = None

    expected_vs_plasma_res_ramp: Any = None

    expected_vs_plasma_total_required: Any = None

    expected_v_plasma_loop_burn: Any = None


@pytest.mark.parametrize(
    "voltsecondreqparam",
    (
        VoltSecondReqParam(
            csawth=1,
            eps=0.33333333333333331,
            f_c_plasma_inductive=0.59999999999999998,
            ejima_coeff=0.30000000000000004,
            kappa=1.8500000000000001,
            plasma_current=18398455.678867526,
            ind_plasma_internal_norm=1.2064840230894305,
            rmajor=8,
            res_plasma=3.7767895536275952e-09,
            t_plant_pulse_burn=1000,
            t_plant_pulse_fusion_ramp=10,
            expected_vs_plasma_internal=111.57651734747576,
            expected_ind_plasma=1.4075705307248088e-05,
            expected_vs_plasma_burn_required=42.109179697761263,
            expected_vs_plasma_ramp_required=314.4596753393147,
            expected_vs_plasma_ind_ramp=258.97124024420435,
            expected_vs_plasma_res_ramp=55.488435095110333,
            expected_vs_plasma_total_required=356.56885503707593,
            expected_v_plasma_loop_burn=0.0416922571264963,
        ),
        VoltSecondReqParam(
            csawth=1,
            eps=0.33333333333333331,
            f_c_plasma_inductive=0.59999999999999998,
            ejima_coeff=0.30000000000000004,
            kappa=1.8500000000000001,
            plasma_current=18398455.678867526,
            ind_plasma_internal_norm=1.2064840230894305,
            rmajor=8,
            res_plasma=3.7767895536275952e-09,
            t_plant_pulse_burn=0,
            t_plant_pulse_fusion_ramp=10,
            expected_vs_plasma_internal=111.57651734747576,
            expected_ind_plasma=1.4075705307248088e-05,
            expected_vs_plasma_burn_required=0.41692257126496302,
            expected_vs_plasma_ramp_required=314.4596753393147,
            expected_vs_plasma_ind_ramp=258.97124024420435,
            expected_vs_plasma_res_ramp=55.488435095110333,
            expected_vs_plasma_total_required=314.87659791057968,
            expected_v_plasma_loop_burn=0.0416922571264963,
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

    (
        vs_plasma_internal,
        ind_plasma,
        vs_plasma_burn_required,
        vs_plasma_ramp_required,
        vs_plasma_ind_ramp,
        vs_plasma_res_ramp,
        vs_plasma_total_required,
        v_plasma_loop_burn,
    ) = calculate_volt_second_requirements(
        csawth=voltsecondreqparam.csawth,
        eps=voltsecondreqparam.eps,
        f_c_plasma_inductive=voltsecondreqparam.f_c_plasma_inductive,
        ejima_coeff=voltsecondreqparam.ejima_coeff,
        kappa=voltsecondreqparam.kappa,
        plasma_current=voltsecondreqparam.plasma_current,
        ind_plasma_internal_norm=voltsecondreqparam.ind_plasma_internal_norm,
        rmajor=voltsecondreqparam.rmajor,
        res_plasma=voltsecondreqparam.res_plasma,
        t_plant_pulse_burn=voltsecondreqparam.t_plant_pulse_burn,
        t_plant_pulse_fusion_ramp=voltsecondreqparam.t_plant_pulse_fusion_ramp,
    )

    assert vs_plasma_internal == pytest.approx(
        voltsecondreqparam.expected_vs_plasma_internal
    )

    assert ind_plasma == pytest.approx(voltsecondreqparam.expected_ind_plasma)

    assert vs_plasma_burn_required == pytest.approx(
        voltsecondreqparam.expected_vs_plasma_burn_required
    )

    assert vs_plasma_ramp_required == pytest.approx(
        voltsecondreqparam.expected_vs_plasma_ramp_required
    )

    assert vs_plasma_ind_ramp == pytest.approx(
        voltsecondreqparam.expected_vs_plasma_ind_ramp
    )

    assert vs_plasma_res_ramp == pytest.approx(
        voltsecondreqparam.expected_vs_plasma_res_ramp
    )

    assert vs_plasma_total_required == pytest.approx(
        voltsecondreqparam.expected_vs_plasma_total_required
    )

    assert v_plasma_loop_burn == pytest.approx(
        voltsecondreqparam.expected_v_plasma_loop_burn
    )


class PhyauxParam(NamedTuple):
    tauratio: Any = None

    burnup_in: Any = None

    aspect: Any = None

    nd_plasma_electrons_vol_avg: Any = None

    te: Any = None

    nd_plasma_fuel_ions_vol_avg: Any = None

    nd_plasma_alphas_vol_avg: Any = None

    fusden_total: Any = None

    fusden_alpha_total: Any = None

    plasma_current: Any = None

    sbar: Any = None

    t_energy_confinement: Any = None

    vol_plasma: Any = None

    expected_burnup: Any = None

    expected_ntau: Any = None

    expected_nTtau: Any = None

    expected_figmer: Any = None

    expected_fusrat: Any = None

    expected_molflow_plasma_fuelling_required: Any = None

    expected_rndfuel: Any = None

    expected_t_alpha_confinement: Any = None


@pytest.mark.parametrize(
    "phyauxparam",
    (
        PhyauxParam(
            tauratio=1,
            burnup_in=0,
            aspect=3,
            nd_plasma_electrons_vol_avg=7.5e19,
            te=12.569,
            nd_plasma_fuel_ions_vol_avg=5.8589175702454272e19,
            nd_plasma_alphas_vol_avg=7.5e18,
            fusden_total=1.9852091609123786e17,
            fusden_alpha_total=1.973996644759543e17,
            plasma_current=18398455.678867526,
            sbar=1,
            t_energy_confinement=3.401323521525641,
            vol_plasma=1888.1711539956691,
            expected_burnup=0.20383432558954095,
            expected_ntau=2.5509926411442307e20,
            expected_nTtau=3.253e21,
            expected_figmer=55.195367036602576,
            expected_fusrat=3.7484146722826997e20,
            expected_molflow_plasma_fuelling_required=1.8389516394951276e21,
            expected_rndfuel=3.7484146722826997e20,
            expected_t_alpha_confinement=37.993985551650177,
        ),
        PhyauxParam(
            tauratio=1,
            burnup_in=0,
            aspect=3,
            nd_plasma_electrons_vol_avg=7.5e19,
            te=12.569,
            nd_plasma_fuel_ions_vol_avg=5.8576156204039725e19,
            nd_plasma_alphas_vol_avg=7.5e18,
            fusden_total=1.9843269653375773e17,
            fusden_alpha_total=1.9731194318497056e17,
            plasma_current=18398455.678867526,
            sbar=1,
            t_energy_confinement=3.402116961408892,
            vol_plasma=1888.1711539956691,
            expected_burnup=0.20387039462081086,
            expected_ntau=2.5515877210566689e20,
            expected_nTtau=3.253e21,
            expected_figmer=55.195367036602576,
            expected_fusrat=3.7467489360461772e20,
            expected_molflow_plasma_fuelling_required=1.8378092331723546e21,
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

    (
        burnup,
        ntau,
        nTtau,
        figmer,
        fusrat,
        molflow_plasma_fuelling_required,
        rndfuel,
        t_alpha_confinement,
        _,
    ) = physics.phyaux(
        aspect=phyauxparam.aspect,
        nd_plasma_electrons_vol_avg=phyauxparam.nd_plasma_electrons_vol_avg,
        te=phyauxparam.te,
        nd_plasma_fuel_ions_vol_avg=phyauxparam.nd_plasma_fuel_ions_vol_avg,
        nd_plasma_alphas_vol_avg=phyauxparam.nd_plasma_alphas_vol_avg,
        fusden_total=phyauxparam.fusden_total,
        fusden_alpha_total=phyauxparam.fusden_alpha_total,
        plasma_current=phyauxparam.plasma_current,
        sbar=phyauxparam.sbar,
        t_energy_confinement=phyauxparam.t_energy_confinement,
        vol_plasma=phyauxparam.vol_plasma,
    )

    assert burnup == pytest.approx(phyauxparam.expected_burnup)

    assert ntau == pytest.approx(phyauxparam.expected_ntau)

    assert figmer == pytest.approx(phyauxparam.expected_figmer)

    assert fusrat == pytest.approx(phyauxparam.expected_fusrat)

    assert molflow_plasma_fuelling_required == pytest.approx(
        phyauxparam.expected_molflow_plasma_fuelling_required
    )

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

    f_c_plasma_inductive: Any = None

    kappa95: Any = None

    plasma_current: Any = None

    rmajor: Any = None

    rminor: Any = None

    temp_plasma_electron_density_weighted_kev: Any = None

    vol_plasma: Any = None

    zeff: Any = None

    expected_pden_plasma_ohmic_mw: Any = None

    expected_p_plasma_ohmic_mw: Any = None

    expected_f_res_plasma_neo: Any = None

    expected_res_plasma: Any = None


@pytest.mark.parametrize(
    "pohmparam",
    (
        PohmParam(
            aspect=3,
            plasma_res_factor=0.70000000000000007,
            f_c_plasma_inductive=0.59999999999999998,
            kappa95=1.6517857142857142,
            plasma_current=18398455.678867526,
            rmajor=8,
            rminor=2.6666666666666665,
            temp_plasma_electron_density_weighted_kev=12.626131115905864,
            vol_plasma=1888.1711539956691,
            zeff=2.0909945616489103,
            expected_pden_plasma_ohmic_mw=0.0004062519138005805,
            expected_p_plasma_ohmic_mw=0.7670731448937912,
            expected_f_res_plasma_neo=2.5,
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
        f_res_plasma_neo,
        res_plasma,
    ) = physics.plasma_ohmic_heating(
        f_c_plasma_inductive=pohmparam.f_c_plasma_inductive,
        kappa95=pohmparam.kappa95,
        plasma_current=pohmparam.plasma_current,
        rmajor=pohmparam.rmajor,
        rminor=pohmparam.rminor,
        temp_plasma_electron_density_weighted_kev=pohmparam.temp_plasma_electron_density_weighted_kev,
        vol_plasma=pohmparam.vol_plasma,
        zeff=pohmparam.zeff,
    )

    assert pden_plasma_ohmic_mw == pytest.approx(
        pohmparam.expected_pden_plasma_ohmic_mw
    )

    assert p_plasma_ohmic_mw == pytest.approx(pohmparam.expected_p_plasma_ohmic_mw)

    assert f_res_plasma_neo == pytest.approx(pohmparam.expected_f_res_plasma_neo)

    assert res_plasma == pytest.approx(pohmparam.expected_res_plasma)


class CalculateDensityLimitParam(NamedTuple):
    i_density_limit: Any = None

    b_plasma_toroidal_on_axis: Any = None

    p_plasma_separatrix_mw: Any = None

    p_hcd_injected_total_mw: Any = None

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
            b_plasma_toroidal_on_axis=5.1847188735686647,
            p_plasma_separatrix_mw=162.32943903093374,
            p_hcd_injected_total_mw=79.928763793309031,
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

    nd_plasma_electron_max_array, nd_plasma_electrons_max = (
        physics.calculate_density_limit(
            i_density_limit=calculatedensitylimitparam.i_density_limit,
            b_plasma_toroidal_on_axis=calculatedensitylimitparam.b_plasma_toroidal_on_axis,
            p_plasma_separatrix_mw=calculatedensitylimitparam.p_plasma_separatrix_mw,
            p_hcd_injected_total_mw=calculatedensitylimitparam.p_hcd_injected_total_mw,
            plasma_current=calculatedensitylimitparam.plasma_current,
            prn1=calculatedensitylimitparam.prn1,
            q95=calculatedensitylimitparam.q95,
            qcyl=calculatedensitylimitparam.qcyl,
            rmajor=calculatedensitylimitparam.rmajor,
            rminor=calculatedensitylimitparam.rminor,
            a_plasma_surface=calculatedensitylimitparam.a_plasma_surface,
            zeff=calculatedensitylimitparam.zeff,
        )
    )

    assert nd_plasma_electrons_max == pytest.approx(
        calculatedensitylimitparam.expected_dnelimt
    )

    assert nd_plasma_electron_max_array == pytest.approx(
        calculatedensitylimitparam.expected_dlimit
    )


class ConfinementTimeParam(NamedTuple):
    i_rad_loss: Any = None

    tauee_in: Any = None

    pden_plasma_rad_mw: Any = None

    kappa_ipb: Any = None

    p_plasma_ohmic_mw: Any = None

    f_p_alpha_plasma_deposited: Any = None

    i_confinement_time: Any = None

    i_plasma_ignited: Any = None

    m_fuel_amu: Any = None

    p_alpha_total_mw: Any = None

    aspect: Any = None

    b_plasma_toroidal_on_axis: Any = None

    nd_plasma_electrons_vol_avg: Any = None

    nd_plasma_ions_total_vol_avg: Any = None

    nd_plasma_electron_line: Any = None

    eps: Any = None

    hfact: Any = None

    kappa: Any = None

    kappa95: Any = None

    p_non_alpha_charged_mw: Any = None

    p_hcd_injected_total_mw: Any = None

    plasma_current: Any = None

    pden_plasma_core_rad_mw: Any = None

    q95: Any = None

    qstar: Any = None

    rmajor: Any = None

    rminor: Any = None

    te: Any = None

    temp_plasma_electron_density_weighted_kev: Any = None

    temp_plasma_ion_density_weighted_kev: Any = None

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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=32,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=6.1946504123280199,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=33,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.96163409847948844,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=34,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.1960655150953154,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=35,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.78453691772934719,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=36,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.1619079679077555,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=37,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.7340642483550435,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=38,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.1117392853827024,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=39,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.84477227311274361,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=40,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=1.6096667103064701,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=41,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.67053301699102119,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=42,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=2.1212580310552207,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=43,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=50.095480115636271,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=44,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=87.869318916638761,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=45,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=28.562137619592285,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=46,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.50082256309019457,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
            f_p_alpha_plasma_deposited=0.94999999999999996,
            i_confinement_time=47,
            i_plasma_ignited=0,
            m_fuel_amu=2.5,
            p_alpha_total_mw=319.03020327154269,
            aspect=3,
            b_plasma_toroidal_on_axis=5.2375830857646202,
            nd_plasma_electrons_vol_avg=8.0593948787884524e19,
            nd_plasma_ions_total_vol_avg=7.1529510234203251e19,
            nd_plasma_electron_line=8.925359201116491e19,
            eps=0.33333333333333331,
            hfact=0.77961193402355955,
            kappa=1.8500000000000001,
            kappa95=1.6517857142857142,
            p_non_alpha_charged_mw=1.2453296074483358,
            p_hcd_injected_total_mw=75.397788712812741,
            plasma_current=16616203.759182997,
            pden_plasma_core_rad_mw=0.047757569353246924,
            q95=3.5610139569387185,
            qstar=2.9513713188821282,
            rmajor=8,
            rminor=2.6666666666666665,
            te=12.437097421317889,
            temp_plasma_electron_density_weighted_kev=13.745148298980761,
            temp_plasma_ion_density_weighted_kev=13.745148298980761,
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
        physics_variables,
        "f_p_alpha_plasma_deposited",
        confinementtimeparam.f_p_alpha_plasma_deposited,
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
        i_plasma_ignited=confinementtimeparam.i_plasma_ignited,
        m_fuel_amu=confinementtimeparam.m_fuel_amu,
        p_alpha_total_mw=confinementtimeparam.p_alpha_total_mw,
        aspect=confinementtimeparam.aspect,
        b_plasma_toroidal_on_axis=confinementtimeparam.b_plasma_toroidal_on_axis,
        nd_plasma_electrons_vol_avg=confinementtimeparam.nd_plasma_electrons_vol_avg,
        nd_plasma_ions_total_vol_avg=confinementtimeparam.nd_plasma_ions_total_vol_avg,
        nd_plasma_electron_line=confinementtimeparam.nd_plasma_electron_line,
        eps=confinementtimeparam.eps,
        hfact=confinementtimeparam.hfact,
        kappa=confinementtimeparam.kappa,
        kappa95=confinementtimeparam.kappa95,
        p_non_alpha_charged_mw=confinementtimeparam.p_non_alpha_charged_mw,
        p_hcd_injected_total_mw=confinementtimeparam.p_hcd_injected_total_mw,
        plasma_current=confinementtimeparam.plasma_current,
        pden_plasma_core_rad_mw=confinementtimeparam.pden_plasma_core_rad_mw,
        q95=confinementtimeparam.q95,
        qstar=confinementtimeparam.qstar,
        rmajor=confinementtimeparam.rmajor,
        rminor=confinementtimeparam.rminor,
        temp_plasma_electron_density_weighted_kev=confinementtimeparam.temp_plasma_electron_density_weighted_kev,
        temp_plasma_ion_density_weighted_kev=confinementtimeparam.temp_plasma_ion_density_weighted_kev,
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


def test_calculate_plasma_masses():
    """Test calculate_plasma_masses()"""
    m_fuel_amu = 2.5
    m_ions_total_amu = 3.0
    nd_plasma_ions_total_vol_avg = 1.0e20
    nd_plasma_fuel_ions_vol_avg = 0.8e20
    nd_plasma_alphas_vol_avg = 0.1e20
    vol_plasma = 100.0
    nd_plasma_electrons_vol_avg = 1.0e20

    (
        m_plasma_fuel_ions,
        m_plasma_ions_total,
        m_plasma_alpha,
        m_plasma_electron,
        m_plasma,
    ) = Physics.calculate_plasma_masses(
        m_fuel_amu=m_fuel_amu,
        m_ions_total_amu=m_ions_total_amu,
        nd_plasma_ions_total_vol_avg=nd_plasma_ions_total_vol_avg,
        nd_plasma_fuel_ions_vol_avg=nd_plasma_fuel_ions_vol_avg,
        nd_plasma_alphas_vol_avg=nd_plasma_alphas_vol_avg,
        vol_plasma=vol_plasma,
        nd_plasma_electrons_vol_avg=nd_plasma_electrons_vol_avg,
    )

    assert m_plasma_fuel_ions == pytest.approx(3.32107813784e-05, abs=1e-30)
    assert m_plasma_ions_total == pytest.approx(4.9816172067599995e-05, abs=1e-30)
    assert m_plasma_alpha == pytest.approx(6.644657345e-06, abs=1e-30)
    assert m_plasma_electron == pytest.approx(9.1093837139e-09, abs=1e-34)
    assert m_plasma == pytest.approx(4.982528145131389e-05, abs=1e-30)


def test_calculate_current_profile_index_wesson():
    """Test calculate_current_profile_index_wesson()."""
    qstar = 3.5
    q0 = 1.5
    result = Physics.calculate_current_profile_index_wesson(qstar, q0)
    assert result == pytest.approx(1.33333, abs=0.0001)


def test_calculate_internal_inductance_wesson():
    """Test calculate_internal_inductance_wesson()."""
    alphaj = 0.8
    result = Physics.calculate_internal_inductance_wesson(alphaj)
    assert result == pytest.approx(0.8595087177751706, abs=0.0001)


def test_calculate_internal_inductance_menard():
    """Test calculate_internal_inductance_menard()."""
    kappa = 2.8
    result = Physics.calculate_internal_inductance_menard(kappa)
    assert result == pytest.approx(0.6, abs=0.001)


def test_calculate_beta_norm_max_wesson():
    """Test calculate_beta_norm_max_wesson()."""
    ind_plasma_internal_norm = 1.5
    result = Physics.calculate_beta_norm_max_wesson(ind_plasma_internal_norm)
    assert result == pytest.approx(6.0, abs=0.001)


def test_calculate_beta_norm_max_original():
    """Test calculate_beta_norm_max_original()"""
    eps = 0.5
    result = Physics.calculate_beta_norm_max_original(eps)
    assert result == pytest.approx(3.8932426932522994, abs=0.00001)


def test_calculate_beta_norm_max_menard():
    """Test calculate_beta_norm_max_menard()."""
    eps = 0.5
    result = Physics.calculate_beta_norm_max_menard(eps)
    assert result == pytest.approx(4.197251361676802, abs=0.000001)


def test_calculate_beta_norm_max_thloreus():
    """Test calculate_beta_norm_max_thloreus()"""
    c_beta = 0.5
    pres_plasma_on_axis = 2.0
    pres_plasma_vol_avg = 1.0
    result = Physics.calculate_beta_norm_max_thloreus(
        c_beta, pres_plasma_on_axis, pres_plasma_vol_avg
    )
    assert result == pytest.approx(5.075, abs=0.00001)


def test_calculate_beta_norm_max_stambaugh():
    """Test calculate_beta_norm_max_thloreus()"""
    f_c_plasma_bootstrap = 0.7
    kappa = 2.0
    aspect = 2.5
    result = Physics.calculate_beta_norm_max_stambaugh(
        f_c_plasma_bootstrap, kappa, aspect
    )
    assert result == pytest.approx(3.840954484207041, abs=0.00001)


def test_calculate_internal_inductance_iter_3():
    """Test calculate_normalised_internal_inductance_iter_3."""
    result = Physics.calculate_normalised_internal_inductance_iter_3(
        b_plasma_poloidal_vol_avg=1.0, c_plasma=1.5e7, vol_plasma=1000.0, rmajor=6.2
    )
    assert result == pytest.approx(0.9078959099585583, abs=0.00001)
