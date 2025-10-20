"""Unit tests for costs.f90."""

from typing import Any, NamedTuple

import numpy as np
import pytest
from pytest import approx

from process import data_structure
from process.costs import Costs
from process.data_structure import (
    build_variables,
    buildings_variables,
    cost_variables,
    current_drive_variables,
    divertor_variables,
    fwbs_variables,
    heat_transport_variables,
    ife_variables,
    pf_power_variables,
    pfcoil_variables,
    physics_variables,
    pulse_variables,
    structure_variables,
    tfcoil_variables,
    times_variables,
    vacuum_variables,
)


@pytest.fixture
def costs():
    """Provides Pulse object for testing.

    :returns: initialised Pulse object
    :rtype: process.pulse.Pulse
    """
    return Costs()


def acc2261_param(**kwargs):
    """Make parameters for a single acc2261() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {"i_blkt_coolant_type": 1, "expected": approx(49.68, abs=0.01)}

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def acc2261_params():
    """Create a list of parameter dicts for the acc2261 fixture.

    Case 1: Reactor cooling system (He)
    Case 2: Reactor cooling system (H2O)

    :return: List of parameter dicts
    :rtype: list
    """
    return [
        acc2261_param(),
        acc2261_param(i_blkt_coolant_type=2, expected=approx(53.85, abs=0.01)),
    ]


@pytest.fixture(params=acc2261_params(), ids=["he", "h2o"])
def acc2261_fix(costs, request, monkeypatch):
    """Fixture for the acc2261() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of acc2261() for the parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock variables used by acc2261()
    monkeypatch.setattr(cost_variables, "fkind", 1)
    monkeypatch.setattr(cost_variables, "lsa", 1)
    monkeypatch.setattr(heat_transport_variables, "p_fw_div_heat_deposited_mw", 0.0)
    monkeypatch.setattr(fwbs_variables, "p_blkt_nuclear_heat_total_mw", 1558.0)
    monkeypatch.setattr(fwbs_variables, "p_shld_nuclear_heat_mw", 1.478)
    monkeypatch.setattr(heat_transport_variables, "p_plant_primary_heat_mw", 2647.0)
    monkeypatch.setattr(heat_transport_variables, "n_primary_heat_exchangers", 3)
    monkeypatch.setattr(cost_variables, "c2261", 0)

    # Parameterised mocks
    monkeypatch.setattr(
        fwbs_variables, "i_blkt_coolant_type", param["i_blkt_coolant_type"]
    )

    # Return the expected result for the given parameter list
    return param["expected"]


def test_acc2261(acc2261_fix, costs):
    """Test acc2261.

    :param acc2261_fix: Expected value of acc2261()
    :type acc2261_fix: ApproxScalar
    """
    # Run acc2261() with the current fixture,
    # then assert the result (c2261) is the expected one
    costs.acc2261()
    assert cost_variables.c2261 == acc2261_fix


def test_acc2262(monkeypatch, costs):
    """Test acc2262()."""
    # Mock module variables
    monkeypatch.setattr(cost_variables, "fkind", 1)
    monkeypatch.setattr(cost_variables, "lsa", 4)
    monkeypatch.setattr(heat_transport_variables, "p_hcd_electric_loss_mw", 76.5)
    monkeypatch.setattr(heat_transport_variables, "p_cryo_plant_electric_mw", 39.936)
    monkeypatch.setattr(heat_transport_variables, "vachtmw", 0.5)
    monkeypatch.setattr(heat_transport_variables, "p_tritium_plant_electric_mw", 15.0)
    monkeypatch.setattr(heat_transport_variables, "fachtmw", 64.835)
    monkeypatch.setattr(cost_variables, "c2262", 0)

    costs.acc2262()
    assert cost_variables.c2262 == approx(29.408, abs=0.01)


def test_acc2263(monkeypatch, costs):
    """Test acc2263().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(cost_variables, "fkind", 1)
    monkeypatch.setattr(cost_variables, "lsa", 4)
    monkeypatch.setattr(cost_variables, "uccry", 9.3e4)
    monkeypatch.setattr(data_structure.tfcoil_variables, "tftmp", 4.5)
    monkeypatch.setattr(heat_transport_variables, "helpow", 80.980e3)
    monkeypatch.setattr(cost_variables, "c2263", 0)

    costs.acc2263()
    assert cost_variables.c2263 == approx(180.76, abs=0.01)


def test_acc2271(monkeypatch, costs):
    """Test acc2271().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(cost_variables, "ucf1", 2.23e7)
    monkeypatch.setattr(cost_variables, "fkind", 1)
    monkeypatch.setattr(cost_variables, "c2271", 0)

    costs.acc2271()
    assert cost_variables.c2271 == approx(22.3, abs=0.01)


def test_acc2272(monkeypatch, costs):
    """Test acc2272().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(physics_variables, "rndfuel", 7.158e20)
    monkeypatch.setattr(physics_variables, "m_fuel_amu", 2.5)
    monkeypatch.setattr(cost_variables, "fkind", 1)
    monkeypatch.setattr(cost_variables, "c2271", 0)

    costs.acc2272()
    assert cost_variables.c2272 == approx(114.707, abs=0.01)


def acc2273_param(**kwargs):
    """Make parameters for a single acc2273() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {
        "f_plasma_fuel_tritium": 0.0001,
        "volrci": data_structure.buildings_variables.volrci,
        "wsvol": data_structure.buildings_variables.wsvol,
        "expected": approx(0.0, abs=0.00001),
    }

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def acc2273_params():
    """Create a list of parameter dicts for the acc2273 fixture.

    Case 1: ttrit low
    Case 2: ttrit high

    :return: List of parameter dicts
    :rtype: list
    """
    return [
        acc2273_param(),
        acc2273_param(
            f_plasma_fuel_tritium=0.5,
            volrci=1299783.4,
            wsvol=132304.1,
            expected=approx(74.12, abs=0.01),
        ),
    ]


@pytest.fixture(params=acc2273_params(), ids=["ttrit_low", "ttrit_high"])
def acc2273_fix(request, monkeypatch, costs):
    """Fixture for the acc2273() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of acc2273() for the parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock variables used by acc2273()
    # Some may be parameterised
    monkeypatch.setattr(data_structure.buildings_variables, "wsvol", param["wsvol"])
    monkeypatch.setattr(data_structure.buildings_variables, "volrci", param["volrci"])
    monkeypatch.setattr(
        physics_variables, "f_plasma_fuel_tritium", param["f_plasma_fuel_tritium"]
    )

    # Mock result var as negative, as an expected result is 0
    # Otherwise could get false positive result
    monkeypatch.setattr(cost_variables, "c2273", -1)

    # Return the expected result for the given parameter list
    return param["expected"]


def test_acc2273(acc2273_fix, costs):
    """Test acc2273.

    :param acc2273_fix: Expected value of acc2273()
    :type acc2273_fix: ApproxScalar
    """
    # Run acc2273() with the current fixture,
    # then assert the result is the expected one
    costs.acc2273()
    assert cost_variables.c2273 == acc2273_fix


def test_acc2274(monkeypatch, costs):
    """Test acc2274().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(data_structure.buildings_variables, "wsvol", 132304.1)
    monkeypatch.setattr(data_structure.buildings_variables, "volrci", 1299783.4)
    monkeypatch.setattr(cost_variables, "fkind", 1)
    monkeypatch.setattr(cost_variables, "c2274", 0)

    costs.acc2274()
    assert cost_variables.c2274 == approx(84.10, abs=0.01)


def acc228_param(**kwargs):
    """Make parameters for a single acc228() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {"fkind": 1.0, "expected": approx(150.0, abs=0.01)}

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def acc228_params():
    """Create a list of parameter dicts for the acc228 fixture.

    Case 1: fkind 1
    Case 2: fkind 0.5

    :return: List of parameter dicts
    :rtype: list
    """
    return [acc228_param(), acc228_param(fkind=0.5, expected=approx(75.0, abs=0.01))]


@pytest.fixture(params=acc228_params(), ids=["fkind_1", "fkind_0p5"])
def acc228_fix(request, monkeypatch, costs):
    """Fixture for the acc228() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of acc228() for the parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock variables used by acc228()
    # Some may be parameterised
    monkeypatch.setattr(cost_variables, "fkind", param["fkind"])
    monkeypatch.setattr(cost_variables, "c228", 0)

    # Return the expected result for the given parameter list
    return param["expected"]


def test_acc228(acc228_fix, costs):
    """Test acc228.

    :param acc228_fix: Expected value of acc228()
    :type acc228_fix: ApproxScalar
    """
    # Run acc228() with the current fixture,
    # then assert the result is the expected one
    costs.acc228()
    assert cost_variables.c228 == acc228_fix


def acc229_param(**kwargs):
    """Make parameters for a single acc229() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {"fkind": 1.0, "expected": approx(125.0, abs=0.01)}

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def acc229_params():
    """Create a list of parameter dicts for the acc229 fixture.

    Case 1: fkind 1
    Case 2: fkind 0.5

    :return: List of parameter dicts
    :rtype: list
    """
    return [acc229_param(), acc229_param(fkind=0.5, expected=approx(62.5, abs=0.01))]


@pytest.fixture(params=acc229_params(), ids=["fkind_1", "fkind_0p5"])
def acc229_fix(request, monkeypatch, costs):
    """Fixture for the acc229() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of acc229() for the parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock variables used by acc229()
    # Some may be parameterised
    monkeypatch.setattr(cost_variables, "fkind", param["fkind"])
    monkeypatch.setattr(cost_variables, "c229", 0)

    # Return the expected result for the given parameter list
    return param["expected"]


def test_acc229(acc229_fix, costs):
    """Test acc229.

    :param acc229_fix: Expected value of acc229()
    :type acc229_fix: ApproxScalar
    """
    # Run acc229() with the current fixture,
    # then assert the result is the expected one
    costs.acc229()
    assert cost_variables.c229 == acc229_fix


def acc23_param(**kwargs):
    """Make parameters for a single acc23() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {"i_blkt_coolant_type": 1, "expected": approx(230, abs=0.01)}

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def acc23_params():
    """Create a list of parameter dicts for the acc23 fixture.

    Case 1: he coolant
    Case 2: h2o coolant

    :return: List of parameter dicts
    :rtype: list
    """
    return [
        acc23_param(),
        acc23_param(i_blkt_coolant_type=2, expected=approx(245, abs=0.01)),
    ]


@pytest.fixture(params=acc23_params(), ids=["he", "h2o"])
def acc23_fix(request, monkeypatch, costs):
    """Fixture for the acc23() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of acc23() for the parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock variables used by acc23()
    # Some may be parameterised
    monkeypatch.setattr(
        fwbs_variables, "i_blkt_coolant_type", param["i_blkt_coolant_type"]
    )
    monkeypatch.setattr(heat_transport_variables, "p_plant_electric_gross_mw", 1200.0)
    monkeypatch.setattr(cost_variables, "c23", 0)

    # Return the expected result for the given parameter list
    return param["expected"]


def test_acc23(acc23_fix, costs):
    """Test acc23.

    :param acc23_fix: Expected value of acc23()
    :type acc23_fix: ApproxScalar
    """
    # Run acc23() with the current fixture,
    # then assert the result is the expected one
    costs.acc23()
    assert cost_variables.c23 == acc23_fix


def test_acc241(monkeypatch, costs):
    """Test acc241().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(cost_variables, "lsa", 4)
    monkeypatch.setattr(cost_variables, "c241", 0)

    costs.acc241()
    assert cost_variables.c241 == approx(18.4, abs=0.01)


def test_acc242(monkeypatch, costs):
    """Test acc242().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(cost_variables, "lsa", 4)
    monkeypatch.setattr(heat_transport_variables, "pacpmw", 630.0)
    monkeypatch.setattr(
        heat_transport_variables, "p_plant_electric_base_total_mw", 65.0
    )
    monkeypatch.setattr(cost_variables, "c242", 0)

    costs.acc242()
    assert cost_variables.c242 == approx(9.06, abs=0.01)


def test_acc243(monkeypatch, costs):
    """Test acc243().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(cost_variables, "lsa", 4)
    monkeypatch.setattr(heat_transport_variables, "tlvpmw", 403.8)
    monkeypatch.setattr(cost_variables, "c243", 0)

    costs.acc243()
    assert cost_variables.c243 == approx(8.08, abs=0.01)


def test_acc244(monkeypatch, costs):
    """Test acc244().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(cost_variables, "lsa", 4)
    monkeypatch.setattr(cost_variables, "c244", 0)

    costs.acc244()
    assert cost_variables.c244 == approx(6.80, abs=0.01)


def test_acc245(monkeypatch, costs):
    """Test acc245().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(cost_variables, "lsa", 4)
    monkeypatch.setattr(cost_variables, "c245", 0)

    costs.acc245()
    assert cost_variables.c245 == approx(1.5, abs=0.01)


def acc25_param(**kwargs):
    """Make parameters for a single acc25() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {"lsa": 4, "expected": approx(25, abs=0.01)}

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def acc25_params():
    """Create a list of parameter dicts for the acc25 fixture.

    Case 1: lsa_4
    Case 2: lsa_1

    :return: List of parameter dicts
    :rtype: list
    """
    return [acc25_param(), acc25_param(lsa=1, expected=approx(19.25, abs=0.01))]


@pytest.fixture(params=acc25_params(), ids=["lsa_4", "lsa_1"])
def acc25_fix(request, monkeypatch, costs):
    """Fixture for the acc25() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of acc25() for the parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock variables used by acc25()
    # Some may be parameterised
    monkeypatch.setattr(cost_variables, "ucmisc", 2.5e7)
    monkeypatch.setattr(cost_variables, "lsa", param["lsa"])
    monkeypatch.setattr(cost_variables, "c25", 0)

    # Return the expected result for the given parameter list
    return param["expected"]


def test_acc25(acc25_fix, costs):
    """Test acc25.

    :param acc25_fix: Expected value of acc25()
    :type acc25_fix: ApproxScalar
    """
    # Run acc25() with the current fixture,
    # then assert the result is the expected one
    costs.acc25()
    assert cost_variables.c25 == acc25_fix


def acc26_param(**kwargs):
    """Make parameters for a single acc26() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {
        "ireactor": 0,
        "p_fusion_total_mw": 2000.0,
        "p_hcd_electric_total_mw": 250.0,
        "tfcmw": 50.0,
        "p_plant_primary_heat_mw": heat_transport_variables.p_plant_primary_heat_mw,
        "p_plant_electric_gross_mw": heat_transport_variables.p_plant_electric_gross_mw,
        "expected": approx(87.9, abs=0.01),
    }

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def acc26_params():
    """Create a list of parameter dicts for the acc26 fixture.

    Case 1: ireactor = 0
    Case 2: ireactor = 1

    :return: List of parameter dicts
    :rtype: list
    """
    return [
        acc26_param(),
        acc26_param(
            ireactor=1,
            p_fusion_total_mw=physics_variables.p_fusion_total_mw,
            p_hcd_electric_total_mw=heat_transport_variables.p_hcd_electric_total_mw,
            tfcmw=data_structure.tfcoil_variables.tfcmw,
            p_plant_primary_heat_mw=3000.0,
            p_plant_electric_gross_mw=700.0,
        ),
    ]


@pytest.fixture(params=acc26_params(), ids=["ireactor_0", "ireactor_1"])
def acc26_fix(request, monkeypatch, costs):
    """Fixture for the acc26() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of acc26() for the parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock variables used by acc26()
    # Some may be parameterised
    monkeypatch.setattr(cost_variables, "lsa", 4)
    monkeypatch.setattr(cost_variables, "ireactor", param["ireactor"])
    monkeypatch.setattr(
        physics_variables,
        "p_fusion_total_mw",
        param["p_fusion_total_mw"],
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_total_mw",
        param["p_hcd_electric_total_mw"],
    )
    monkeypatch.setattr(data_structure.tfcoil_variables, "tfcmw", param["tfcmw"])
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_primary_heat_mw",
        param["p_plant_primary_heat_mw"],
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_gross_mw",
        param["p_plant_electric_gross_mw"],
    )
    monkeypatch.setattr(cost_variables, "c26", 0)

    # Return the expected result for the given parameter list
    return param["expected"]


def test_acc26(acc26_fix, costs):
    """Test acc26.

    :param acc26_fix: Expected value of acc26()
    :type acc26_fix: ApproxScalar
    """
    # Run acc26() with the current fixture,
    # then assert the result is the expected one
    costs.acc26()
    assert cost_variables.c26 == acc26_fix


def test_acc9(monkeypatch, costs):
    """Test acc9().

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    monkeypatch.setattr(cost_variables, "lsa", 4)
    monkeypatch.setattr(cost_variables, "cdirt", 30e3)
    monkeypatch.setattr(cost_variables, "cowner", 0.15)
    monkeypatch.setattr(cost_variables, "fcontng", 0.195)
    monkeypatch.setattr(cost_variables, "cindrt", 0)

    costs.acc9()
    assert cost_variables.cindrt == approx(10005.0, abs=0.1)
    assert cost_variables.ccont == approx(7800.98, abs=0.1)


class Acc21Param(NamedTuple):
    shovol: Any = None

    triv: Any = None

    elevol: Any = None

    rbvol: Any = None

    cryvol: Any = None

    rmbvol: Any = None

    admvol: Any = None

    convol: Any = None

    wsvol: Any = None

    ucrb: Any = None

    ireactor: Any = None

    cturbb: Any = None

    lsa: Any = None

    csi: Any = None

    cland: Any = None

    c21: Any = None

    c211: Any = None

    c212: Any = None

    c213: Any = None

    c214: Any = None

    c2141: Any = None

    c2142: Any = None

    c215: Any = None

    c216: Any = None

    c217: Any = None

    c2171: Any = None

    c2172: Any = None

    c2173: Any = None

    c2174: Any = None

    expected_c21: Any = None

    expected_c211: Any = None

    expected_c212: Any = None

    expected_c213: Any = None

    expected_c214: Any = None

    expected_c2141: Any = None

    expected_c2142: Any = None

    expected_c215: Any = None

    expected_c216: Any = None

    expected_c217: Any = None

    expected_c2171: Any = None

    expected_c2172: Any = None

    expected_c2173: Any = None

    expected_c2174: Any = None


@pytest.mark.parametrize(
    "acc21param",
    (
        Acc21Param(
            shovol=100000,
            triv=40000,
            elevol=51601.097615432001,
            rbvol=1356973.2891062023,
            cryvol=15247.180612719381,
            rmbvol=421473.52130148414,
            admvol=100000,
            convol=60000,
            wsvol=130018.25667917728,
            ucrb=400,
            ireactor=1,
            cturbb=38,
            lsa=2,
            csi=16,
            cland=19.199999999999999,
            c21=0,
            c211=0,
            c212=0,
            c213=0,
            c214=0,
            c2141=0,
            c2142=0,
            c215=0,
            c216=0,
            c217=0,
            c2171=0,
            c2172=0,
            c2173=0,
            c2174=0,
            expected_c21=740.00647752036286,
            expected_c211=32.640000000000001,
            expected_c212=455.94302513968393,
            expected_c213=31.919999999999998,
            expected_c214=142.28887143307821,
            expected_c2141=92.049817052244123,
            expected_c2142=50.239054380834098,
            expected_c215=12.431999999999999,
            expected_c216=16.471070358845893,
            expected_c217=48.311510588754764,
            expected_c2171=15.119999999999999,
            expected_c2172=17.640000000000001,
            expected_c2173=9.6599999999999984,
            expected_c2174=5.8915105887547687,
        ),
        Acc21Param(
            shovol=100000,
            triv=40000,
            elevol=51609.268177478581,
            rbvol=1358540.6868905292,
            cryvol=25826.919937316459,
            rmbvol=423252.94369581528,
            admvol=100000,
            convol=60000,
            wsvol=130255.93791329287,
            ucrb=400,
            ireactor=1,
            cturbb=38,
            lsa=2,
            csi=16,
            cland=19.199999999999999,
            c21=740.00647752036286,
            c211=32.640000000000001,
            c212=455.94302513968393,
            c213=31.919999999999998,
            c214=142.28887143307821,
            c2141=92.049817052244123,
            c2142=50.239054380834098,
            c215=12.431999999999999,
            c216=16.471070358845893,
            c217=48.311510588754764,
            c2171=15.119999999999999,
            c2172=17.640000000000001,
            c2173=9.6599999999999984,
            c2174=5.8915105887547687,
            expected_c21=745.10420837411039,
            expected_c211=32.640000000000001,
            expected_c212=456.46967079521778,
            expected_c213=31.919999999999998,
            expected_c214=142.76933731286238,
            expected_c2141=92.438442903166035,
            expected_c2142=50.330894409696356,
            expected_c215=12.431999999999999,
            expected_c216=16.47367840225116,
            expected_c217=52.399521863779071,
            expected_c2171=15.119999999999999,
            expected_c2172=17.640000000000001,
            expected_c2173=9.6599999999999984,
            expected_c2174=9.9795218637790786,
        ),
    ),
)
def test_acc21(acc21param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc21.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc21param: the data used to mock and assert in this test.
    :type acc21param: acc21param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(buildings_variables, "shovol", acc21param.shovol)

    monkeypatch.setattr(buildings_variables, "triv", acc21param.triv)

    monkeypatch.setattr(buildings_variables, "elevol", acc21param.elevol)

    monkeypatch.setattr(buildings_variables, "rbvol", acc21param.rbvol)

    monkeypatch.setattr(buildings_variables, "cryvol", acc21param.cryvol)

    monkeypatch.setattr(buildings_variables, "rmbvol", acc21param.rmbvol)

    monkeypatch.setattr(buildings_variables, "admvol", acc21param.admvol)

    monkeypatch.setattr(buildings_variables, "convol", acc21param.convol)

    monkeypatch.setattr(buildings_variables, "wsvol", acc21param.wsvol)

    monkeypatch.setattr(cost_variables, "ucrb", acc21param.ucrb)

    monkeypatch.setattr(cost_variables, "ireactor", acc21param.ireactor)

    monkeypatch.setattr(cost_variables, "cturbb", acc21param.cturbb)

    monkeypatch.setattr(cost_variables, "lsa", acc21param.lsa)

    monkeypatch.setattr(cost_variables, "csi", acc21param.csi)

    monkeypatch.setattr(cost_variables, "cland", acc21param.cland)

    monkeypatch.setattr(cost_variables, "c21", acc21param.c21)

    monkeypatch.setattr(cost_variables, "c211", acc21param.c211)

    monkeypatch.setattr(cost_variables, "c212", acc21param.c212)

    monkeypatch.setattr(cost_variables, "c213", acc21param.c213)

    monkeypatch.setattr(cost_variables, "c214", acc21param.c214)

    monkeypatch.setattr(cost_variables, "c2141", acc21param.c2141)

    monkeypatch.setattr(cost_variables, "c2142", acc21param.c2142)

    monkeypatch.setattr(cost_variables, "c215", acc21param.c215)

    monkeypatch.setattr(cost_variables, "c216", acc21param.c216)

    monkeypatch.setattr(cost_variables, "c217", acc21param.c217)

    monkeypatch.setattr(cost_variables, "c2171", acc21param.c2171)

    monkeypatch.setattr(cost_variables, "c2172", acc21param.c2172)

    monkeypatch.setattr(cost_variables, "c2173", acc21param.c2173)

    monkeypatch.setattr(cost_variables, "c2174", acc21param.c2174)

    costs.acc21()

    assert cost_variables.c21 == pytest.approx(acc21param.expected_c21)

    assert cost_variables.c211 == pytest.approx(acc21param.expected_c211)

    assert cost_variables.c212 == pytest.approx(acc21param.expected_c212)

    assert cost_variables.c213 == pytest.approx(acc21param.expected_c213)

    assert cost_variables.c214 == pytest.approx(acc21param.expected_c214)

    assert cost_variables.c2141 == pytest.approx(acc21param.expected_c2141)

    assert cost_variables.c2142 == pytest.approx(acc21param.expected_c2142)

    assert cost_variables.c215 == pytest.approx(acc21param.expected_c215)

    assert cost_variables.c216 == pytest.approx(acc21param.expected_c216)

    assert cost_variables.c217 == pytest.approx(acc21param.expected_c217)

    assert cost_variables.c2171 == pytest.approx(acc21param.expected_c2171)

    assert cost_variables.c2172 == pytest.approx(acc21param.expected_c2172)

    assert cost_variables.c2173 == pytest.approx(acc21param.expected_c2173)

    assert cost_variables.c2174 == pytest.approx(acc21param.expected_c2174)


class Acc2211Param(NamedTuple):
    a_fw_total: Any = None

    ucblss: Any = None

    fkind: Any = None

    fwallcst: Any = None

    ucblli2o: Any = None

    ifueltyp: Any = None

    lsa: Any = None

    fwmatm: Any = None

    uccarb: Any = None

    ife: Any = None

    ucconc: Any = None

    c22: Any = None

    c2211: Any = None

    expected_fwallcst: Any = None


@pytest.mark.parametrize(
    "acc2211param",
    (
        Acc2211Param(
            a_fw_total=1601.1595634509963,
            ucblss=90,
            fkind=1,
            fwallcst=0,
            ucblli2o=600,
            ifueltyp=1,
            lsa=2,
            fwmatm=np.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            uccarb=50,
            ife=0,
            ucconc=0.10000000000000001,
            c22=0,
            c2211=0,
            expected_fwallcst=143.19827300247195,
        ),
        Acc2211Param(
            a_fw_total=1891.2865102700493,
            ucblss=90,
            fkind=1,
            fwallcst=143.19827300247195,
            ucblli2o=600,
            ifueltyp=1,
            lsa=2,
            fwmatm=np.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            uccarb=50,
            ife=0,
            ucconc=0.10000000000000001,
            c22=3474.7391916096453,
            c2211=0,
            expected_fwallcst=167.7865317453867,
        ),
    ),
)
def test_acc2211(acc2211param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2211.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2211param: the data used to mock and assert in this test.
    :type acc2211param: acc2211param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "a_fw_total", acc2211param.a_fw_total)

    monkeypatch.setattr(cost_variables, "ucblss", acc2211param.ucblss)

    monkeypatch.setattr(cost_variables, "fkind", acc2211param.fkind)

    monkeypatch.setattr(cost_variables, "fwallcst", acc2211param.fwallcst)

    monkeypatch.setattr(cost_variables, "ucblli2o", acc2211param.ucblli2o)

    monkeypatch.setattr(cost_variables, "ifueltyp", acc2211param.ifueltyp)

    monkeypatch.setattr(cost_variables, "lsa", acc2211param.lsa)

    monkeypatch.setattr(ife_variables, "fwmatm", acc2211param.fwmatm)

    monkeypatch.setattr(ife_variables, "uccarb", acc2211param.uccarb)

    monkeypatch.setattr(ife_variables, "ife", acc2211param.ife)

    monkeypatch.setattr(ife_variables, "ucconc", acc2211param.ucconc)

    monkeypatch.setattr(cost_variables, "c22", acc2211param.c22)

    monkeypatch.setattr(cost_variables, "c2211", acc2211param.c2211)

    costs.acc2211()

    assert cost_variables.fwallcst == pytest.approx(acc2211param.expected_fwallcst)


class Acc2212Param(NamedTuple):
    ucblss: Any = None

    ucblbreed: Any = None

    ucblbe: Any = None

    ucblli: Any = None

    ucblvd: Any = None

    ucblli2o: Any = None

    blkcst: Any = None

    ucbllipb: Any = None

    ifueltyp: Any = None

    lsa: Any = None

    fkind: Any = None

    i_blanket_type: Any = None

    m_blkt_lithium: Any = None

    m_blkt_li2o: Any = None

    whtblbreed: Any = None

    m_blkt_vanadium: Any = None

    m_blkt_beryllium: Any = None

    m_blkt_steel_total: Any = None

    wtbllipb: Any = None

    ucflib: Any = None

    blmatm: Any = None

    ife: Any = None

    ucconc: Any = None

    mflibe: Any = None

    uccarb: Any = None

    c22: Any = None

    c2212: Any = None

    c22121: Any = None

    c22122: Any = None

    c22123: Any = None

    c22124: Any = None

    c22125: Any = None

    c22126: Any = None

    c22127: Any = None

    c22128: Any = None

    expected_blkcst: Any = None

    expected_c22121: Any = None

    expected_c22122: Any = None

    expected_c22123: Any = None


@pytest.mark.parametrize(
    "acc2212param",
    (
        Acc2212Param(
            ucblss=90,
            ucblbreed=875,
            ucblbe=260,
            ucblli=875,
            ucblvd=280,
            ucblli2o=600,
            blkcst=0,
            ucbllipb=10.300000000000001,
            ifueltyp=1,
            lsa=2,
            fkind=1,
            i_blanket_type=1,
            m_blkt_lithium=0,
            m_blkt_li2o=1258110.2710352642,
            whtblbreed=0,
            m_blkt_vanadium=0,
            m_blkt_beryllium=1184720.5052248738,
            m_blkt_steel_total=1058196.5489677608,
            wtbllipb=0,
            ucflib=84,
            blmatm=np.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            ife=0,
            ucconc=0.10000000000000001,
            mflibe=0,
            uccarb=50,
            c22=0,
            c2212=0,
            c22121=0,
            c22122=0,
            c22123=0,
            c22124=0,
            c22125=0,
            c22126=0,
            c22127=0,
            c22128=0,
            expected_blkcst=868.59838754004318,
            expected_c22121=231.02049851885039,
            expected_c22122=566.14962196586885,
            expected_c22123=71.428267055323843,
        ),
        Acc2212Param(
            ucblss=90,
            ucblbreed=875,
            ucblbe=260,
            ucblli=875,
            ucblvd=280,
            ucblli2o=600,
            blkcst=868.59838754004318,
            ucbllipb=10.300000000000001,
            ifueltyp=1,
            lsa=2,
            fkind=1,
            i_blanket_type=1,
            m_blkt_lithium=0,
            m_blkt_li2o=1260437.468838267,
            whtblbreed=0,
            m_blkt_vanadium=0,
            m_blkt_beryllium=1186911.9498227015,
            m_blkt_steel_total=1060153.955039866,
            wtbllipb=0,
            ucflib=84,
            blmatm=np.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            ife=0,
            ucconc=0.10000000000000001,
            mflibe=0,
            uccarb=50,
            c22=3474.7391916096453,
            c2212=0,
            c22121=231.02049851885039,
            c22122=566.14962196586885,
            c22123=71.428267055323843,
            c22124=0,
            c22125=0,
            c22126=0,
            c22127=0,
            c22128=0,
            expected_blkcst=870.20508315783786,
            expected_c22121=231.44783021542679,
            expected_c22122=567.19686097722013,
            expected_c22123=71.560391965190959,
        ),
    ),
)
def test_acc2212(acc2212param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2212.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2212param: the data used to mock and assert in this test.
    :type acc2212param: acc2212param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucblss", acc2212param.ucblss)

    monkeypatch.setattr(cost_variables, "ucblbreed", acc2212param.ucblbreed)

    monkeypatch.setattr(cost_variables, "ucblbe", acc2212param.ucblbe)

    monkeypatch.setattr(cost_variables, "ucblli", acc2212param.ucblli)

    monkeypatch.setattr(cost_variables, "ucblvd", acc2212param.ucblvd)

    monkeypatch.setattr(cost_variables, "ucblli2o", acc2212param.ucblli2o)

    monkeypatch.setattr(cost_variables, "blkcst", acc2212param.blkcst)

    monkeypatch.setattr(cost_variables, "ucbllipb", acc2212param.ucbllipb)

    monkeypatch.setattr(cost_variables, "ifueltyp", acc2212param.ifueltyp)

    monkeypatch.setattr(cost_variables, "lsa", acc2212param.lsa)

    monkeypatch.setattr(cost_variables, "fkind", acc2212param.fkind)

    monkeypatch.setattr(fwbs_variables, "i_blanket_type", acc2212param.i_blanket_type)

    monkeypatch.setattr(fwbs_variables, "m_blkt_lithium", acc2212param.m_blkt_lithium)

    monkeypatch.setattr(fwbs_variables, "m_blkt_li2o", acc2212param.m_blkt_li2o)

    monkeypatch.setattr(fwbs_variables, "whtblbreed", acc2212param.whtblbreed)

    monkeypatch.setattr(fwbs_variables, "m_blkt_vanadium", acc2212param.m_blkt_vanadium)

    monkeypatch.setattr(
        fwbs_variables, "m_blkt_beryllium", acc2212param.m_blkt_beryllium
    )

    monkeypatch.setattr(
        fwbs_variables, "m_blkt_steel_total", acc2212param.m_blkt_steel_total
    )

    monkeypatch.setattr(fwbs_variables, "wtbllipb", acc2212param.wtbllipb)

    monkeypatch.setattr(ife_variables, "ucflib", acc2212param.ucflib)

    monkeypatch.setattr(ife_variables, "blmatm", acc2212param.blmatm)

    monkeypatch.setattr(ife_variables, "ife", acc2212param.ife)

    monkeypatch.setattr(ife_variables, "ucconc", acc2212param.ucconc)

    monkeypatch.setattr(ife_variables, "mflibe", acc2212param.mflibe)

    monkeypatch.setattr(ife_variables, "uccarb", acc2212param.uccarb)

    monkeypatch.setattr(cost_variables, "c22", acc2212param.c22)

    monkeypatch.setattr(cost_variables, "c2212", acc2212param.c2212)

    monkeypatch.setattr(cost_variables, "c22121", acc2212param.c22121)

    monkeypatch.setattr(cost_variables, "c22122", acc2212param.c22122)

    monkeypatch.setattr(cost_variables, "c22123", acc2212param.c22123)

    monkeypatch.setattr(cost_variables, "c22124", acc2212param.c22124)

    monkeypatch.setattr(cost_variables, "c22125", acc2212param.c22125)

    monkeypatch.setattr(cost_variables, "c22126", acc2212param.c22126)

    monkeypatch.setattr(cost_variables, "c22127", acc2212param.c22127)

    monkeypatch.setattr(cost_variables, "c22128", acc2212param.c22128)

    costs.acc2212()

    assert cost_variables.blkcst == pytest.approx(acc2212param.expected_blkcst)

    assert cost_variables.c22121 == pytest.approx(acc2212param.expected_c22121)

    assert cost_variables.c22122 == pytest.approx(acc2212param.expected_c22122)

    assert cost_variables.c22123 == pytest.approx(acc2212param.expected_c22123)


class Acc2213Param(NamedTuple):
    ucpens: Any = None

    ucshld: Any = None

    fkind: Any = None

    ucblli2o: Any = None

    lsa: Any = None

    wpenshld: Any = None

    whtshld: Any = None

    shmatm: Any = None

    uccarb: Any = None

    ife: Any = None

    ucconc: Any = None

    c22: Any = None

    c2213: Any = None

    c22131: Any = None

    c22132: Any = None

    expected_c2213: Any = None

    expected_c22131: Any = None

    expected_c22132: Any = None


@pytest.mark.parametrize(
    "acc2213param",
    (
        Acc2213Param(
            ucpens=32,
            ucshld=32,
            fkind=1,
            ucblli2o=600,
            lsa=2,
            wpenshld=2294873.8131476026,
            whtshld=2294873.8131476026,
            shmatm=np.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            uccarb=50,
            ife=0,
            ucconc=0.10000000000000001,
            c22=0,
            c2213=0,
            c22131=0,
            c22132=0,
            expected_c2213=110.15394303108492,
            expected_c22131=55.076971515542461,
            expected_c22132=55.076971515542461,
        ),
        Acc2213Param(
            ucpens=32,
            ucshld=32,
            fkind=1,
            ucblli2o=600,
            lsa=2,
            wpenshld=2297808.3935174868,
            whtshld=2297808.3935174868,
            shmatm=np.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            uccarb=50,
            ife=0,
            ucconc=0.10000000000000001,
            c22=3474.7391916096453,
            c2213=110.15394303108492,
            c22131=55.076971515542461,
            c22132=55.076971515542461,
            expected_c2213=110.29480288883934,
            expected_c22131=55.147401444419671,
            expected_c22132=55.147401444419671,
        ),
    ),
)
def test_acc2213(acc2213param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2213.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2213param: the data used to mock and assert in this test.
    :type acc2213param: acc2213param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucpens", acc2213param.ucpens)

    monkeypatch.setattr(cost_variables, "ucshld", acc2213param.ucshld)

    monkeypatch.setattr(cost_variables, "fkind", acc2213param.fkind)

    monkeypatch.setattr(cost_variables, "ucblli2o", acc2213param.ucblli2o)

    monkeypatch.setattr(cost_variables, "lsa", acc2213param.lsa)

    monkeypatch.setattr(fwbs_variables, "wpenshld", acc2213param.wpenshld)

    monkeypatch.setattr(fwbs_variables, "whtshld", acc2213param.whtshld)

    monkeypatch.setattr(ife_variables, "shmatm", acc2213param.shmatm)

    monkeypatch.setattr(ife_variables, "uccarb", acc2213param.uccarb)

    monkeypatch.setattr(ife_variables, "ife", acc2213param.ife)

    monkeypatch.setattr(ife_variables, "ucconc", acc2213param.ucconc)

    monkeypatch.setattr(cost_variables, "c22", acc2213param.c22)

    monkeypatch.setattr(cost_variables, "c2213", acc2213param.c2213)

    monkeypatch.setattr(cost_variables, "c22131", acc2213param.c22131)

    monkeypatch.setattr(cost_variables, "c22132", acc2213param.c22132)

    costs.acc2213()

    assert cost_variables.c2213 == pytest.approx(acc2213param.expected_c2213)

    assert cost_variables.c22131 == pytest.approx(acc2213param.expected_c22131)

    assert cost_variables.c22132 == pytest.approx(acc2213param.expected_c22132)


class Acc2214Param(NamedTuple):
    fkind: Any = None

    lsa: Any = None

    gsmass: Any = None

    c22: Any = None

    c2214: Any = None

    expected_c2214: Any = None


@pytest.mark.parametrize(
    "acc2214param",
    (
        Acc2214Param(
            fkind=1,
            lsa=2,
            gsmass=1631228.030796848,
            c22=0,
            c2214=0,
            expected_c2214=47.672639200037878,
        ),
        Acc2214Param(
            fkind=1,
            lsa=2,
            gsmass=1626877.8363395864,
            c22=3474.7391916096453,
            c2214=47.672639200037878,
            expected_c2214=47.545504767024411,
        ),
    ),
)
def test_acc2214(acc2214param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2214.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2214param: the data used to mock and assert in this test.
    :type acc2214param: acc2214param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "fkind", acc2214param.fkind)

    monkeypatch.setattr(cost_variables, "lsa", acc2214param.lsa)

    monkeypatch.setattr(structure_variables, "gsmass", acc2214param.gsmass)

    monkeypatch.setattr(cost_variables, "c22", acc2214param.c22)

    monkeypatch.setattr(cost_variables, "c2214", acc2214param.c2214)

    costs.acc2214()

    assert cost_variables.c2214 == pytest.approx(acc2214param.expected_c2214)


class Acc2215Param(NamedTuple):
    ifueltyp: Any = None

    divcst: Any = None

    fkind: Any = None

    ucdiv: Any = None

    a_div_surface_total: Any = None

    ife: Any = None

    c22: Any = None

    c2215: Any = None

    expected_divcst: Any = None


@pytest.mark.parametrize(
    "acc2215param",
    (
        Acc2215Param(
            ifueltyp=1,
            divcst=0,
            fkind=1,
            ucdiv=500000,
            a_div_surface_total=177.80928909705162,
            ife=0,
            c22=0,
            c2215=0,
            expected_divcst=88.904644548525795,
        ),
        Acc2215Param(
            ifueltyp=1,
            divcst=88.904644548525795,
            fkind=1,
            ucdiv=500000,
            a_div_surface_total=177.80928909705162,
            ife=0,
            c22=3474.7391916096453,
            c2215=0,
            expected_divcst=88.904644548525795,
        ),
    ),
)
def test_acc2215(acc2215param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2215.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2215param: the data used to mock and assert in this test.
    :type acc2215param: acc2215param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ifueltyp", acc2215param.ifueltyp)

    monkeypatch.setattr(cost_variables, "divcst", acc2215param.divcst)

    monkeypatch.setattr(cost_variables, "fkind", acc2215param.fkind)

    monkeypatch.setattr(cost_variables, "ucdiv", acc2215param.ucdiv)

    monkeypatch.setattr(
        divertor_variables, "a_div_surface_total", acc2215param.a_div_surface_total
    )

    monkeypatch.setattr(ife_variables, "ife", acc2215param.ife)

    monkeypatch.setattr(cost_variables, "c22", acc2215param.c22)

    monkeypatch.setattr(cost_variables, "c2215", acc2215param.c2215)

    costs.acc2215()

    assert cost_variables.divcst == pytest.approx(acc2215param.expected_divcst)


class Acc2221Param(NamedTuple):
    uccpclb: Any = None

    uccase: Any = None

    uccu: Any = None

    fkind: Any = None

    cconshtf: Any = None

    ucsc: Any = None

    ifueltyp: Any = None

    uccpcl1: Any = None

    ucwindtf: Any = None

    cpstcst: Any = None

    lsa: Any = None

    cconfix: Any = None

    itart: Any = None

    clgsmass: Any = None

    aintmass: Any = None

    m_tf_coil_copper: Any = None

    m_tf_coil_superconductor: Any = None

    m_tf_coil_case: Any = None

    n_tf_coils: Any = None

    whttflgs: Any = None

    whtcp: Any = None

    i_tf_sup: Any = None

    supercond_cost_model: Any = None

    j_crit_str_tf: Any = None

    n_tf_coil_turns: Any = None

    len_tf_coil: Any = None

    i_tf_sc_mat: Any = None

    c22: Any = None

    c2221: Any = None

    c22211: Any = None

    c22212: Any = None

    c22213: Any = None

    c22214: Any = None

    c22215: Any = None

    expected_c22211: Any = None

    expected_c22212: Any = None


@pytest.mark.parametrize(
    "acc2221param",
    (
        Acc2221Param(
            uccpclb=150,
            uccase=50,
            uccu=75,
            fkind=1,
            cconshtf=75,
            ucsc=np.array(
                np.array((600, 600, 300, 600, 600, 600, 300, 1200, 1200), order="F"),
                order="F",
            ).transpose(),
            ifueltyp=1,
            uccpcl1=250,
            ucwindtf=480,
            cpstcst=0,
            lsa=2,
            cconfix=80,
            itart=0,
            clgsmass=1953582.3684708222,
            aintmass=5829865.436088616,
            m_tf_coil_copper=58744.465423173802,
            m_tf_coil_superconductor=5802.5700395134345,
            m_tf_coil_case=1034021.9996272125,
            n_tf_coils=16,
            whttflgs=0,
            whtcp=0,
            i_tf_sup=1,
            supercond_cost_model=0,
            j_crit_str_tf=300.0,
            n_tf_coil_turns=200,
            len_tf_coil=50.483843027201402,
            i_tf_sc_mat=5,
            c22=0,
            c2221=0,
            c22211=0,
            c22212=0,
            c22213=0,
            c22214=0,
            c22215=0,
            expected_c22211=127.79612438919186,
            expected_c22212=65.523989541865234,
        ),
        Acc2221Param(
            uccpclb=150,
            uccase=50,
            uccu=75,
            fkind=1,
            cconshtf=75,
            ucsc=np.array(
                np.array((600, 600, 300, 600, 600, 600, 300, 1200, 1200), order="F"),
                order="F",
            ).transpose(),
            ifueltyp=1,
            uccpcl1=250,
            ucwindtf=480,
            cpstcst=0,
            lsa=2,
            cconfix=80,
            itart=0,
            clgsmass=1951781.4798732549,
            aintmass=5829865.436088616,
            m_tf_coil_copper=58779.575542593491,
            m_tf_coil_superconductor=5806.038092640837,
            m_tf_coil_case=1034699.2182961091,
            n_tf_coils=16,
            whttflgs=0,
            whtcp=0,
            i_tf_sup=1,
            supercond_cost_model=0,
            j_crit_str_tf=300.0,
            n_tf_coil_turns=200,
            len_tf_coil=50.514015976170839,
            i_tf_sc_mat=5,
            c22=3474.7391916096453,
            c2221=1122.5144544988982,
            c22211=127.79612438919186,
            c22212=65.523989541865234,
            c22213=698.99887174799562,
            c22214=172.4182702723208,
            c22215=57.77719854752457,
            expected_c22211=127.87250498362496,
            expected_c22212=65.563151615791654,
        ),
        Acc2221Param(
            uccpclb=150,
            uccase=50,
            uccu=75,
            fkind=1,
            cconshtf=75,
            ucsc=np.array(
                np.array((600, 600, 300, 600, 600, 600, 300, 1200, 1200), order="F"),
                order="F",
            ).transpose(),
            ifueltyp=1,
            uccpcl1=250,
            ucwindtf=480,
            cpstcst=0,
            lsa=2,
            cconfix=80,
            itart=0,
            clgsmass=1951781.4798732549,
            aintmass=5829865.436088616,
            m_tf_coil_copper=58779.575542593491,
            m_tf_coil_superconductor=5806.038092640837,
            m_tf_coil_case=1034699.2182961091,
            n_tf_coils=16,
            whttflgs=0,
            whtcp=0,
            i_tf_sup=1,
            supercond_cost_model=1,
            j_crit_str_tf=300.0,
            n_tf_coil_turns=200,
            len_tf_coil=50.514015976170839,
            i_tf_sc_mat=5,
            c22=3474.7391916096453,
            c2221=1122.5144544988982,
            c22211=127.79612438919186,
            c22212=65.523989541865234,
            c22213=698.99887174799562,
            c22214=172.4182702723208,
            c22215=57.77719854752457,
            expected_c22211=1462760.833721748,
            expected_c22212=65.563151615791654,
        ),
    ),
)
def test_acc2221(acc2221param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2221.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2221param: the data used to mock and assert in this test.
    :type acc2221param: acc2221param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "uccpclb", acc2221param.uccpclb)

    monkeypatch.setattr(cost_variables, "uccase", acc2221param.uccase)

    monkeypatch.setattr(cost_variables, "uccu", acc2221param.uccu)

    monkeypatch.setattr(cost_variables, "fkind", acc2221param.fkind)

    monkeypatch.setattr(cost_variables, "cconshtf", acc2221param.cconshtf)

    monkeypatch.setattr(cost_variables, "ucsc", acc2221param.ucsc)

    monkeypatch.setattr(cost_variables, "ifueltyp", acc2221param.ifueltyp)

    monkeypatch.setattr(cost_variables, "uccpcl1", acc2221param.uccpcl1)

    monkeypatch.setattr(cost_variables, "ucwindtf", acc2221param.ucwindtf)

    monkeypatch.setattr(cost_variables, "cpstcst", acc2221param.cpstcst)

    monkeypatch.setattr(cost_variables, "lsa", acc2221param.lsa)

    monkeypatch.setattr(cost_variables, "cconfix", acc2221param.cconfix)

    monkeypatch.setattr(physics_variables, "itart", acc2221param.itart)

    monkeypatch.setattr(structure_variables, "clgsmass", acc2221param.clgsmass)

    monkeypatch.setattr(structure_variables, "aintmass", acc2221param.aintmass)

    monkeypatch.setattr(
        tfcoil_variables, "m_tf_coil_copper", acc2221param.m_tf_coil_copper
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "m_tf_coil_superconductor",
        acc2221param.m_tf_coil_superconductor,
    )

    monkeypatch.setattr(tfcoil_variables, "m_tf_coil_case", acc2221param.m_tf_coil_case)

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", acc2221param.n_tf_coils)

    monkeypatch.setattr(tfcoil_variables, "whttflgs", acc2221param.whttflgs)

    monkeypatch.setattr(tfcoil_variables, "whtcp", acc2221param.whtcp)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", acc2221param.i_tf_sup)

    monkeypatch.setattr(
        cost_variables, "supercond_cost_model", acc2221param.supercond_cost_model
    )

    monkeypatch.setattr(tfcoil_variables, "j_crit_str_tf", acc2221param.j_crit_str_tf)

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_coil_turns", acc2221param.n_tf_coil_turns
    )

    monkeypatch.setattr(tfcoil_variables, "len_tf_coil", acc2221param.len_tf_coil)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sc_mat", acc2221param.i_tf_sc_mat)

    monkeypatch.setattr(cost_variables, "c22", acc2221param.c22)

    monkeypatch.setattr(cost_variables, "c2221", acc2221param.c2221)

    monkeypatch.setattr(cost_variables, "c22211", acc2221param.c22211)

    monkeypatch.setattr(cost_variables, "c22212", acc2221param.c22212)

    monkeypatch.setattr(cost_variables, "c22213", acc2221param.c22213)

    monkeypatch.setattr(cost_variables, "c22214", acc2221param.c22214)

    monkeypatch.setattr(cost_variables, "c22215", acc2221param.c22215)

    costs.acc2221()

    assert cost_variables.c22211 == pytest.approx(acc2221param.expected_c22211)

    assert cost_variables.c22212 == pytest.approx(acc2221param.expected_c22212)


class Acc2222Param(NamedTuple):
    iohcl: Any = None

    uccase: Any = None

    uccu: Any = None

    cconshpf: Any = None

    ucfnc: Any = None

    cconfix: Any = None

    ucsc: Any = None

    ucwindpf: Any = None

    lsa: Any = None

    fkind: Any = None

    j_pf_coil_wp_peak: Any = None

    supercond_cost_model: Any = None

    j_crit_str_cs: Any = None

    j_crit_str_pf: Any = None

    i_pf_conductor: Any = None

    f_a_cs_void: Any = None

    n_cs_pf_coils: Any = None

    n_pf_coil_turns: Any = None

    i_pf_superconductor: Any = None

    m_pf_coil_structure_total: Any = None

    c_pf_cs_coils_peak_ma: Any = None

    r_pf_coil_middle: Any = None

    i_cs_superconductor: Any = None

    fcupfsu: Any = None

    fcuohsu: Any = None

    f_a_pf_coil_void: Any = None

    awpoh: Any = None

    fncmass: Any = None

    dcond: Any = None

    c22: Any = None

    c2222: Any = None

    c22221: Any = None

    c22222: Any = None

    c22223: Any = None

    c22224: Any = None

    expected_c2222: Any = None

    expected_c22221: Any = None

    expected_c22222: Any = None

    expected_c22223: Any = None

    expected_c22224: Any = None


@pytest.mark.parametrize(
    "acc2222param",
    (
        Acc2222Param(
            iohcl=1,
            uccase=50,
            uccu=75,
            cconshpf=70,
            ucfnc=35,
            cconfix=80,
            ucsc=np.array(
                np.array((600, 600, 300, 600, 600, 600, 300, 1200, 1200), order="F"),
                order="F",
            ).transpose(),
            ucwindpf=465,
            lsa=2,
            fkind=1,
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
            supercond_cost_model=0,
            j_crit_str_cs=100.0,
            j_crit_str_pf=200.0,
            i_pf_conductor=0,
            f_a_cs_void=0.29999999999999999,
            n_cs_pf_coils=7,
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
            i_pf_superconductor=3,
            m_pf_coil_structure_total=2695737.563343476,
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
            i_cs_superconductor=5,
            fcupfsu=0.68999999999999995,
            fcuohsu=0.70000000000000007,
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
            awpoh=3.8004675824985918,
            fncmass=310716.52923547616,
            dcond=np.array(
                np.array(
                    (6080, 6080, 6070, 6080, 6080, 8500, 6070, 8500, 8500),
                    order="F",
                ),
                order="F",
            ).transpose(),
            c22=0,
            c2222=0,
            c22221=0,
            c22222=0,
            c22223=0,
            c22224=0,
            expected_c2222=626.57984594974835,
            expected_c22221=434.46640986938519,
            expected_c22222=69.02908267696219,
            expected_c22223=113.89491205126185,
            expected_c22224=9.1894413521392071,
        ),
        Acc2222Param(
            iohcl=1,
            uccase=50,
            uccu=75,
            cconshpf=70,
            ucfnc=35,
            cconfix=80,
            ucsc=np.array(
                np.array((600, 600, 300, 600, 600, 600, 300, 1200, 1200), order="F"),
                order="F",
            ).transpose(),
            ucwindpf=465,
            lsa=2,
            fkind=1,
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
            supercond_cost_model=0,
            j_crit_str_cs=100.0,
            j_crit_str_pf=200.0,
            i_pf_conductor=0,
            f_a_cs_void=0.29999999999999999,
            n_cs_pf_coils=7,
            n_pf_coil_turns=np.array(
                np.array(
                    (
                        440.26292595093469,
                        525.4843415877815,
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
            i_pf_superconductor=3,
            m_pf_coil_structure_total=2510424.9065680322,
            c_pf_cs_coils_peak_ma=np.array(
                np.array(
                    (
                        18.579095475129446,
                        22.175439215004378,
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
            i_cs_superconductor=5,
            fcupfsu=0.68999999999999995,
            fcuohsu=0.70000000000000007,
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
            awpoh=3.8004675824985918,
            fncmass=310716.52923547616,
            dcond=np.array(
                np.array(
                    (6080, 6080, 6070, 6080, 6080, 8500, 6070, 8500, 8500),
                    order="F",
                ),
                order="F",
            ).transpose(),
            c22=3474.7391916096453,
            c2222=626.57984594974835,
            c22221=434.46640986938519,
            c22222=69.02908267696219,
            c22223=113.89491205126185,
            c22224=9.1894413521392071,
            expected_c2222=634.503192513881,
            expected_c22221=448.04573758127646,
            expected_c22222=71.202561277966055,
            expected_c22223=106.06545230249935,
            expected_c22224=9.1894413521392071,
        ),
        Acc2222Param(
            iohcl=1,
            uccase=50,
            uccu=75,
            cconshpf=70,
            ucfnc=35,
            cconfix=80,
            ucsc=np.array(
                np.array((600, 600, 300, 600, 600, 600, 300, 1200, 1200), order="F"),
                order="F",
            ).transpose(),
            ucwindpf=465,
            lsa=2,
            fkind=1,
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
            supercond_cost_model=1,
            j_crit_str_cs=100.0,
            j_crit_str_pf=200.0,
            i_pf_conductor=0,
            f_a_cs_void=0.29999999999999999,
            n_cs_pf_coils=7,
            n_pf_coil_turns=np.array(
                np.array(
                    (
                        440.26292595093469,
                        525.4843415877815,
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
            i_pf_superconductor=3,
            m_pf_coil_structure_total=2510424.9065680322,
            c_pf_cs_coils_peak_ma=np.array(
                np.array(
                    (
                        18.579095475129446,
                        22.175439215004378,
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
            i_cs_superconductor=5,
            fcupfsu=0.68999999999999995,
            fcuohsu=0.70000000000000007,
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
            awpoh=3.8004675824985918,
            fncmass=310716.52923547616,
            dcond=np.array(
                np.array(
                    (6080, 6080, 6070, 6080, 6080, 8500, 6070, 8500, 8500),
                    order="F",
                ),
                order="F",
            ).transpose(),
            c22=3474.7391916096453,
            c2222=626.57984594974835,
            c22221=434.46640986938519,
            c22222=69.02908267696219,
            c22223=113.89491205126185,
            c22224=9.1894413521392071,
            expected_c2222=2271626.1414324627,
            expected_c22221=2271439.6839775303,
            expected_c22222=71.202561277966055,
            expected_c22223=106.06545230249935,
            expected_c22224=9.1894413521392071,
        ),
    ),
)
def test_acc2222(acc2222param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2222.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2222param: the data used to mock and assert in this test.
    :type acc2222param: acc2222param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "iohcl", acc2222param.iohcl)

    monkeypatch.setattr(cost_variables, "uccase", acc2222param.uccase)

    monkeypatch.setattr(cost_variables, "uccu", acc2222param.uccu)

    monkeypatch.setattr(cost_variables, "cconshpf", acc2222param.cconshpf)

    monkeypatch.setattr(cost_variables, "ucfnc", acc2222param.ucfnc)

    monkeypatch.setattr(cost_variables, "cconfix", acc2222param.cconfix)

    monkeypatch.setattr(cost_variables, "ucsc", acc2222param.ucsc)

    monkeypatch.setattr(cost_variables, "ucwindpf", acc2222param.ucwindpf)

    monkeypatch.setattr(cost_variables, "lsa", acc2222param.lsa)

    monkeypatch.setattr(cost_variables, "fkind", acc2222param.fkind)

    monkeypatch.setattr(
        pfcoil_variables, "j_pf_coil_wp_peak", acc2222param.j_pf_coil_wp_peak
    )

    monkeypatch.setattr(
        cost_variables, "supercond_cost_model", acc2222param.supercond_cost_model
    )

    monkeypatch.setattr(pfcoil_variables, "j_crit_str_cs", acc2222param.j_crit_str_cs)

    monkeypatch.setattr(pfcoil_variables, "j_crit_str_pf", acc2222param.j_crit_str_pf)

    monkeypatch.setattr(pfcoil_variables, "i_pf_conductor", acc2222param.i_pf_conductor)

    monkeypatch.setattr(pfcoil_variables, "f_a_cs_void", acc2222param.f_a_cs_void)

    monkeypatch.setattr(pfcoil_variables, "n_cs_pf_coils", acc2222param.n_cs_pf_coils)

    monkeypatch.setattr(
        pfcoil_variables, "n_pf_coil_turns", acc2222param.n_pf_coil_turns
    )

    monkeypatch.setattr(
        pfcoil_variables, "i_pf_superconductor", acc2222param.i_pf_superconductor
    )

    monkeypatch.setattr(
        pfcoil_variables,
        "m_pf_coil_structure_total",
        acc2222param.m_pf_coil_structure_total,
    )

    monkeypatch.setattr(
        pfcoil_variables, "c_pf_cs_coils_peak_ma", acc2222param.c_pf_cs_coils_peak_ma
    )

    monkeypatch.setattr(
        pfcoil_variables, "r_pf_coil_middle", acc2222param.r_pf_coil_middle
    )

    monkeypatch.setattr(
        pfcoil_variables, "i_cs_superconductor", acc2222param.i_cs_superconductor
    )

    monkeypatch.setattr(pfcoil_variables, "fcupfsu", acc2222param.fcupfsu)

    monkeypatch.setattr(pfcoil_variables, "fcuohsu", acc2222param.fcuohsu)

    monkeypatch.setattr(
        pfcoil_variables, "f_a_pf_coil_void", acc2222param.f_a_pf_coil_void
    )

    monkeypatch.setattr(pfcoil_variables, "awpoh", acc2222param.awpoh)

    monkeypatch.setattr(structure_variables, "fncmass", acc2222param.fncmass)

    monkeypatch.setattr(tfcoil_variables, "dcond", acc2222param.dcond)

    monkeypatch.setattr(cost_variables, "c22", acc2222param.c22)

    monkeypatch.setattr(cost_variables, "c2222", acc2222param.c2222)

    monkeypatch.setattr(cost_variables, "c22221", acc2222param.c22221)

    monkeypatch.setattr(cost_variables, "c22222", acc2222param.c22222)

    monkeypatch.setattr(cost_variables, "c22223", acc2222param.c22223)

    monkeypatch.setattr(cost_variables, "c22224", acc2222param.c22224)

    costs.acc2222()

    assert cost_variables.c2222 == pytest.approx(acc2222param.expected_c2222)

    assert cost_variables.c22221 == pytest.approx(acc2222param.expected_c22221)

    assert cost_variables.c22222 == pytest.approx(acc2222param.expected_c22222)

    assert cost_variables.c22223 == pytest.approx(acc2222param.expected_c22223)

    assert cost_variables.c22224 == pytest.approx(acc2222param.expected_c22224)


class Acc2223Param(NamedTuple):
    uccryo: Any = None

    lsa: Any = None

    fkind: Any = None

    m_vv: Any = None

    c22: Any = None

    c2223: Any = None

    expected_c2223: Any = None


@pytest.mark.parametrize(
    "acc2223param",
    (
        Acc2223Param(
            uccryo=32,
            lsa=2,
            fkind=1,
            m_vv=9043937.8018644415,
            c22=0,
            c2223=0,
            expected_c2223=244.54807816241447,
        ),
        Acc2223Param(
            uccryo=32,
            lsa=2,
            fkind=1,
            m_vv=9056931.558219457,
            c22=3474.7391916096453,
            c2223=244.54807816241447,
            expected_c2223=244.89942933425411,
        ),
    ),
)
def test_acc2223(acc2223param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2223.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2223param: the data used to mock and assert in this test.
    :type acc2223param: acc2223param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "uccryo", acc2223param.uccryo)

    monkeypatch.setattr(cost_variables, "lsa", acc2223param.lsa)

    monkeypatch.setattr(cost_variables, "fkind", acc2223param.fkind)

    monkeypatch.setattr(fwbs_variables, "m_vv", acc2223param.m_vv)

    monkeypatch.setattr(cost_variables, "c22", acc2223param.c22)

    monkeypatch.setattr(cost_variables, "c2223", acc2223param.c2223)

    costs.acc2223()

    assert cost_variables.c2223 == pytest.approx(acc2223param.expected_c2223)


class Acc223Param(NamedTuple):
    ucich: Any = None

    fkind: Any = None

    ucnbi: Any = None

    ucech: Any = None

    uclh: Any = None

    ifueltyp: Any = None

    cdcost: Any = None

    fcdfuel: Any = None

    p_hcd_lowhyb_injected_total_mw: Any = None

    i_hcd_primary: Any = None

    p_hcd_ecrh_injected_total_mw: Any = None

    p_beam_injected_mw: Any = None

    dcdrv2: Any = None

    mcdriv: Any = None

    cdriv2: Any = None

    dcdrv0: Any = None

    edrive: Any = None

    etadrv: Any = None

    ifedrv: Any = None

    ife: Any = None

    dcdrv1: Any = None

    cdriv1: Any = None

    cdriv3: Any = None

    cdriv0: Any = None

    c22: Any = None

    c223: Any = None

    c2231: Any = None

    c2232: Any = None

    c2233: Any = None

    c2234: Any = None

    expected_cdcost: Any = None

    expected_c223: Any = None

    expected_c2231: Any = None


@pytest.mark.parametrize(
    "acc223param",
    (
        Acc223Param(
            ucich=3,
            fkind=1,
            ucnbi=3.2999999999999998,
            ucech=3,
            uclh=3.2999999999999998,
            ifueltyp=1,
            cdcost=0,
            fcdfuel=0.10000000000000001,
            p_hcd_lowhyb_injected_total_mw=0,
            i_hcd_primary=10,
            p_hcd_ecrh_injected_total_mw=51.978447720428512,
            p_beam_injected_mw=0,
            dcdrv2=59.899999999999999,
            mcdriv=1,
            cdriv2=244.90000000000001,
            dcdrv0=111.40000000000001,
            edrive=5000000,
            etadrv=0,
            ifedrv=2,
            ife=0,
            dcdrv1=78,
            cdriv1=163.19999999999999,
            cdriv3=1.4630000000000001,
            cdriv0=154.30000000000001,
            c22=0,
            c223=0,
            c2231=0,
            c2232=0,
            c2233=0,
            c2234=0,
            expected_cdcost=140.341808845157,
            expected_c223=140.341808845157,
            expected_c2231=140.341808845157,
        ),
        Acc223Param(
            ucich=3,
            fkind=1,
            ucnbi=3.2999999999999998,
            ucech=3,
            uclh=3.2999999999999998,
            ifueltyp=1,
            cdcost=140.341808845157,
            fcdfuel=0.10000000000000001,
            p_hcd_lowhyb_injected_total_mw=0,
            i_hcd_primary=10,
            p_hcd_ecrh_injected_total_mw=51.978447720428512,
            p_beam_injected_mw=0,
            dcdrv2=59.899999999999999,
            mcdriv=1,
            cdriv2=244.90000000000001,
            dcdrv0=111.40000000000001,
            edrive=5000000,
            etadrv=0,
            ifedrv=2,
            ife=0,
            dcdrv1=78,
            cdriv1=163.19999999999999,
            cdriv3=1.4630000000000001,
            cdriv0=154.30000000000001,
            c22=3474.7391916096453,
            c223=140.341808845157,
            c2231=140.341808845157,
            c2232=0,
            c2233=0,
            c2234=0,
            expected_cdcost=140.341808845157,
            expected_c223=140.341808845157,
            expected_c2231=140.341808845157,
        ),
    ),
)
def test_acc223(acc223param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc223.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc223param: the data used to mock and assert in this test.
    :type acc223param: acc223param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucich", acc223param.ucich)

    monkeypatch.setattr(cost_variables, "fkind", acc223param.fkind)

    monkeypatch.setattr(cost_variables, "ucnbi", acc223param.ucnbi)

    monkeypatch.setattr(cost_variables, "ucech", acc223param.ucech)

    monkeypatch.setattr(cost_variables, "uclh", acc223param.uclh)

    monkeypatch.setattr(cost_variables, "ifueltyp", acc223param.ifueltyp)

    monkeypatch.setattr(cost_variables, "cdcost", acc223param.cdcost)

    monkeypatch.setattr(cost_variables, "fcdfuel", acc223param.fcdfuel)

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_lowhyb_injected_total_mw",
        acc223param.p_hcd_lowhyb_injected_total_mw,
    )

    monkeypatch.setattr(
        current_drive_variables, "i_hcd_primary", acc223param.i_hcd_primary
    )

    monkeypatch.setattr(
        current_drive_variables,
        "p_hcd_ecrh_injected_total_mw",
        acc223param.p_hcd_ecrh_injected_total_mw,
    )

    monkeypatch.setattr(
        current_drive_variables, "p_beam_injected_mw", acc223param.p_beam_injected_mw
    )

    monkeypatch.setattr(ife_variables, "dcdrv2", acc223param.dcdrv2)

    monkeypatch.setattr(ife_variables, "mcdriv", acc223param.mcdriv)

    monkeypatch.setattr(ife_variables, "cdriv2", acc223param.cdriv2)

    monkeypatch.setattr(ife_variables, "dcdrv0", acc223param.dcdrv0)

    monkeypatch.setattr(ife_variables, "edrive", acc223param.edrive)

    monkeypatch.setattr(ife_variables, "etadrv", acc223param.etadrv)

    monkeypatch.setattr(ife_variables, "ifedrv", acc223param.ifedrv)

    monkeypatch.setattr(ife_variables, "ife", acc223param.ife)

    monkeypatch.setattr(ife_variables, "dcdrv1", acc223param.dcdrv1)

    monkeypatch.setattr(ife_variables, "cdriv1", acc223param.cdriv1)

    monkeypatch.setattr(ife_variables, "cdriv3", acc223param.cdriv3)

    monkeypatch.setattr(ife_variables, "cdriv0", acc223param.cdriv0)

    monkeypatch.setattr(cost_variables, "c22", acc223param.c22)

    monkeypatch.setattr(cost_variables, "c223", acc223param.c223)

    monkeypatch.setattr(cost_variables, "c2231", acc223param.c2231)

    monkeypatch.setattr(cost_variables, "c2232", acc223param.c2232)

    monkeypatch.setattr(cost_variables, "c2233", acc223param.c2233)

    monkeypatch.setattr(cost_variables, "c2234", acc223param.c2234)

    costs.acc223()

    assert cost_variables.cdcost == pytest.approx(acc223param.expected_cdcost)

    assert cost_variables.c223 == pytest.approx(acc223param.expected_c223)

    assert cost_variables.c2231 == pytest.approx(acc223param.expected_c2231)


class Acc224Param(NamedTuple):
    fkind: Any = None

    dlscal: Any = None

    m_vv_vacuum_duct_shield: Any = None

    n_vac_pumps_high: Any = None

    dia_vv_vacuum_ducts: Any = None

    i_vacuum_pump_type: Any = None

    n_vv_vacuum_ducts: Any = None

    c22: Any = None

    c224: Any = None

    c2241: Any = None

    c2242: Any = None

    c2243: Any = None

    c2244: Any = None

    c2245: Any = None

    c2246: Any = None

    expected_c224: Any = None

    expected_c2241: Any = None

    expected_c2242: Any = None

    expected_c2243: Any = None

    expected_c2244: Any = None

    expected_c2246: Any = None


@pytest.mark.parametrize(
    "acc224param",
    (
        Acc224Param(
            fkind=1,
            dlscal=4.9196133171476717,
            m_vv_vacuum_duct_shield=0,
            n_vac_pumps_high=46,
            dia_vv_vacuum_ducts=0.57081858183821432,
            i_vacuum_pump_type=1,
            n_vv_vacuum_ducts=16,
            c22=0,
            c224=0,
            c2241=0,
            c2242=0,
            c2243=0,
            c2244=0,
            c2245=0,
            c2246=0,
            expected_c224=34.593599813216727,
            expected_c2241=17.940000000000001,
            expected_c2242=4.6799999999999997,
            expected_c2243=3.3256586023918255,
            expected_c2244=7.3479412108249003,
            expected_c2246=1.3,
        ),
        Acc224Param(
            fkind=1,
            dlscal=4.9184638394909044,
            m_vv_vacuum_duct_shield=0,
            n_vac_pumps_high=46,
            dia_vv_vacuum_ducts=0.57072331228476758,
            i_vacuum_pump_type=1,
            n_vv_vacuum_ducts=16,
            c22=3474.7391916096453,
            c224=34.593599813216727,
            c2241=17.940000000000001,
            c2242=4.6799999999999997,
            c2243=3.3256586023918255,
            c2244=7.3479412108249003,
            c2245=0,
            c2246=1.3,
            expected_c224=34.591105904913036,
            expected_c2241=17.940000000000001,
            expected_c2242=4.6799999999999997,
            expected_c2243=3.3248815554958511,
            expected_c2244=7.346224349417187,
            expected_c2246=1.3,
        ),
    ),
)
def test_acc224(acc224param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc224.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc224param: the data used to mock and assert in this test.
    :type acc224param: acc224param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "fkind", acc224param.fkind)

    monkeypatch.setattr(vacuum_variables, "dlscal", acc224param.dlscal)

    monkeypatch.setattr(
        vacuum_variables, "m_vv_vacuum_duct_shield", acc224param.m_vv_vacuum_duct_shield
    )

    monkeypatch.setattr(
        vacuum_variables, "n_vac_pumps_high", acc224param.n_vac_pumps_high
    )

    monkeypatch.setattr(
        vacuum_variables, "dia_vv_vacuum_ducts", acc224param.dia_vv_vacuum_ducts
    )

    monkeypatch.setattr(
        vacuum_variables, "i_vacuum_pump_type", acc224param.i_vacuum_pump_type
    )

    monkeypatch.setattr(
        vacuum_variables, "n_vv_vacuum_ducts", acc224param.n_vv_vacuum_ducts
    )

    monkeypatch.setattr(cost_variables, "c22", acc224param.c22)

    monkeypatch.setattr(cost_variables, "c224", acc224param.c224)

    monkeypatch.setattr(cost_variables, "c2241", acc224param.c2241)

    monkeypatch.setattr(cost_variables, "c2242", acc224param.c2242)

    monkeypatch.setattr(cost_variables, "c2243", acc224param.c2243)

    monkeypatch.setattr(cost_variables, "c2244", acc224param.c2244)

    monkeypatch.setattr(cost_variables, "c2245", acc224param.c2245)

    monkeypatch.setattr(cost_variables, "c2246", acc224param.c2246)

    costs.acc224()

    assert cost_variables.c224 == pytest.approx(acc224param.expected_c224)

    assert cost_variables.c2241 == pytest.approx(acc224param.expected_c2241)

    assert cost_variables.c2242 == pytest.approx(acc224param.expected_c2242)

    assert cost_variables.c2243 == pytest.approx(acc224param.expected_c2243)

    assert cost_variables.c2244 == pytest.approx(acc224param.expected_c2244)

    assert cost_variables.c2246 == pytest.approx(acc224param.expected_c2246)


class Acc2251Param(NamedTuple):
    uctfsw: Any = None

    fkind: Any = None

    ucbus: Any = None

    uctfbr: Any = None

    uctfps: Any = None

    uctfbus: Any = None

    v_tf_coil_dump_quench_kv: Any = None

    tfcmw: Any = None

    len_tf_bus: Any = None

    e_tf_magnetic_stored_total_gj: Any = None

    i_tf_sup: Any = None

    m_tf_bus: Any = None

    tfckw: Any = None

    n_tf_coils: Any = None

    c_tf_turn: Any = None

    c22: Any = None

    c225: Any = None

    c2251: Any = None

    c22511: Any = None

    c22512: Any = None

    c22513: Any = None

    c22514: Any = None

    c22515: Any = None

    expected_c2251: Any = None

    expected_c22511: Any = None

    expected_c22512: Any = None

    expected_c22513: Any = None

    expected_c22514: Any = None

    expected_c22515: Any = None


@pytest.mark.parametrize(
    "acc2251param",
    (
        Acc2251Param(
            uctfsw=1,
            fkind=1,
            ucbus=0.123,
            uctfbr=1.22,
            uctfps=24,
            uctfbus=100,
            v_tf_coil_dump_quench_kv=9.9882637896807953,
            tfcmw=0,
            len_tf_bus=3397.0129827974288,
            e_tf_magnetic_stored_total_gj=152.78343648685947,
            i_tf_sup=1,
            m_tf_bus=0,
            tfckw=32474.753636211804,
            n_tf_coils=16,
            c_tf_turn=74026.751437500003,
            c22=0,
            c225=0,
            c2251=0,
            c22511=0,
            c22512=0,
            c22513=0,
            c22514=0,
            c22515=0,
            expected_c2251=98.457845594540643,
            expected_c22511=4.3480381629432125,
            expected_c22512=31.601916254373826,
            expected_c22513=26.777101385200407,
            expected_c22514=4.7999999999999998,
            expected_c22515=30.930789792023205,
        ),
        Acc2251Param(
            uctfsw=1,
            fkind=1,
            ucbus=0.123,
            uctfbr=1.22,
            uctfps=24,
            uctfbus=100,
            v_tf_coil_dump_quench_kv=10.001287165953382,
            tfcmw=0,
            len_tf_bus=3397.0129827974288,
            e_tf_magnetic_stored_total_gj=152.98264590137683,
            i_tf_sup=1,
            m_tf_bus=0,
            tfckw=32505.257577809778,
            n_tf_coils=16,
            c_tf_turn=74026.751437500003,
            c22=3474.7391916096453,
            c225=185.05656643685359,
            c2251=98.457845594540643,
            c22511=4.3480381629432125,
            c22512=31.601916254373826,
            c22513=26.777101385200407,
            c22514=4.7999999999999998,
            c22515=30.930789792023205,
            expected_c2251=98.524335872804144,
            expected_c22511=4.3508966768725132,
            expected_c22512=31.630686371167478,
            expected_c22513=26.811963032740941,
            expected_c22514=4.7999999999999998,
            expected_c22515=30.930789792023205,
        ),
    ),
)
def test_acc2251(acc2251param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2251.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2251param: the data used to mock and assert in this test.
    :type acc2251param: acc2251param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "uctfsw", acc2251param.uctfsw)

    monkeypatch.setattr(cost_variables, "fkind", acc2251param.fkind)

    monkeypatch.setattr(cost_variables, "ucbus", acc2251param.ucbus)

    monkeypatch.setattr(cost_variables, "uctfbr", acc2251param.uctfbr)

    monkeypatch.setattr(cost_variables, "uctfps", acc2251param.uctfps)

    monkeypatch.setattr(cost_variables, "uctfbus", acc2251param.uctfbus)

    monkeypatch.setattr(
        tfcoil_variables,
        "v_tf_coil_dump_quench_kv",
        acc2251param.v_tf_coil_dump_quench_kv,
    )

    monkeypatch.setattr(tfcoil_variables, "tfcmw", acc2251param.tfcmw)

    monkeypatch.setattr(tfcoil_variables, "len_tf_bus", acc2251param.len_tf_bus)

    monkeypatch.setattr(
        tfcoil_variables,
        "e_tf_magnetic_stored_total_gj",
        acc2251param.e_tf_magnetic_stored_total_gj,
    )

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", acc2251param.i_tf_sup)

    monkeypatch.setattr(tfcoil_variables, "m_tf_bus", acc2251param.m_tf_bus)

    monkeypatch.setattr(tfcoil_variables, "tfckw", acc2251param.tfckw)

    monkeypatch.setattr(tfcoil_variables, "n_tf_coils", acc2251param.n_tf_coils)

    monkeypatch.setattr(tfcoil_variables, "c_tf_turn", acc2251param.c_tf_turn)

    monkeypatch.setattr(cost_variables, "c22", acc2251param.c22)

    monkeypatch.setattr(cost_variables, "c225", acc2251param.c225)

    monkeypatch.setattr(cost_variables, "c2251", acc2251param.c2251)

    monkeypatch.setattr(cost_variables, "c22511", acc2251param.c22511)

    monkeypatch.setattr(cost_variables, "c22512", acc2251param.c22512)

    monkeypatch.setattr(cost_variables, "c22513", acc2251param.c22513)

    monkeypatch.setattr(cost_variables, "c22514", acc2251param.c22514)

    monkeypatch.setattr(cost_variables, "c22515", acc2251param.c22515)

    costs.acc2251()

    assert cost_variables.c2251 == pytest.approx(acc2251param.expected_c2251)

    assert cost_variables.c22511 == pytest.approx(acc2251param.expected_c22511)

    assert cost_variables.c22512 == pytest.approx(acc2251param.expected_c22512)

    assert cost_variables.c22513 == pytest.approx(acc2251param.expected_c22513)

    assert cost_variables.c22514 == pytest.approx(acc2251param.expected_c22514)

    assert cost_variables.c22515 == pytest.approx(acc2251param.expected_c22515)


class Acc2252Param(NamedTuple):
    ucpfcb: Any = None

    ucpfbk: Any = None

    fkind: Any = None

    ucpfb: Any = None

    ucpfdr1: Any = None

    ucpfic: Any = None

    ucpfbs: Any = None

    ucpfps: Any = None

    peakmva: Any = None

    ensxpfm: Any = None

    spfbusl: Any = None

    pfckts: Any = None

    srcktpm: Any = None

    vpfskv: Any = None

    acptmax: Any = None

    c22: Any = None

    c225: Any = None

    c2252: Any = None

    c22521: Any = None

    c22522: Any = None

    c22523: Any = None

    c22524: Any = None

    c22525: Any = None

    c22526: Any = None

    c22527: Any = None

    expected_c22521: Any = None

    expected_c22522: Any = None

    expected_c22523: Any = None

    expected_c22524: Any = None


@pytest.mark.parametrize(
    "acc2252param",
    (
        Acc2252Param(
            ucpfcb=75000,
            ucpfbk=16600,
            fkind=1,
            ucpfb=210,
            ucpfdr1=150,
            ucpfic=10000,
            ucpfbs=4900,
            ucpfps=35000,
            peakmva=736.39062584245937,
            ensxpfm=37429.525515086898,
            spfbusl=2533.4495999999999,
            pfckts=12,
            srcktpm=1071.1112934857531,
            vpfskv=20,
            acptmax=24.816666666666666,
            c22=0,
            c225=0,
            c2252=0,
            c22521=0,
            c22522=0,
            c22523=0,
            c22524=0,
            c22525=0,
            c22526=0,
            c22527=0,
            expected_c22521=25.773671904486076,
            expected_c22522=3.5999999999999996,
            expected_c22523=13.203072590399998,
            expected_c22524=1.36406376579542,
        ),
        Acc2252Param(
            ucpfcb=75000,
            ucpfbk=16600,
            fkind=1,
            ucpfb=210,
            ucpfdr1=150,
            ucpfic=10000,
            ucpfbs=4900,
            ucpfps=35000,
            peakmva=90.673341440806084,
            ensxpfm=37427.228965055205,
            spfbusl=2533.4495999999999,
            pfckts=12,
            srcktpm=1069.8879533693198,
            vpfskv=20,
            acptmax=24.816666666666666,
            c22=3474.7391916096453,
            c225=185.05656643685359,
            c2252=65.813098499070378,
            c22521=25.773671904486076,
            c22522=3.5999999999999996,
            c22523=13.203072590399998,
            c22524=1.36406376579542,
            c22525=15.357861411125848,
            c22526=5.6144288272630343,
            c22527=0.89999999999999991,
            expected_c22521=3.1735669504282127,
            expected_c22522=3.5999999999999996,
            expected_c22523=13.203072590399998,
            expected_c22524=1.3629730294999658,
        ),
    ),
)
def test_acc2252(acc2252param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2252.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2252param: the data used to mock and assert in this test.
    :type acc2252param: acc2252param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucpfcb", acc2252param.ucpfcb)

    monkeypatch.setattr(cost_variables, "ucpfbk", acc2252param.ucpfbk)

    monkeypatch.setattr(cost_variables, "fkind", acc2252param.fkind)

    monkeypatch.setattr(cost_variables, "ucpfb", acc2252param.ucpfb)

    monkeypatch.setattr(cost_variables, "ucpfdr1", acc2252param.ucpfdr1)

    monkeypatch.setattr(cost_variables, "ucpfic", acc2252param.ucpfic)

    monkeypatch.setattr(cost_variables, "ucpfbs", acc2252param.ucpfbs)

    monkeypatch.setattr(cost_variables, "ucpfps", acc2252param.ucpfps)

    monkeypatch.setattr(heat_transport_variables, "peakmva", acc2252param.peakmva)

    monkeypatch.setattr(pf_power_variables, "ensxpfm", acc2252param.ensxpfm)

    monkeypatch.setattr(pf_power_variables, "spfbusl", acc2252param.spfbusl)

    monkeypatch.setattr(pf_power_variables, "pfckts", acc2252param.pfckts)

    monkeypatch.setattr(pf_power_variables, "srcktpm", acc2252param.srcktpm)

    monkeypatch.setattr(pf_power_variables, "vpfskv", acc2252param.vpfskv)

    monkeypatch.setattr(pf_power_variables, "acptmax", acc2252param.acptmax)

    monkeypatch.setattr(cost_variables, "c22", acc2252param.c22)

    monkeypatch.setattr(cost_variables, "c225", acc2252param.c225)

    monkeypatch.setattr(cost_variables, "c2252", acc2252param.c2252)

    monkeypatch.setattr(cost_variables, "c22521", acc2252param.c22521)

    monkeypatch.setattr(cost_variables, "c22522", acc2252param.c22522)

    monkeypatch.setattr(cost_variables, "c22523", acc2252param.c22523)

    monkeypatch.setattr(cost_variables, "c22524", acc2252param.c22524)

    monkeypatch.setattr(cost_variables, "c22525", acc2252param.c22525)

    monkeypatch.setattr(cost_variables, "c22526", acc2252param.c22526)

    monkeypatch.setattr(cost_variables, "c22527", acc2252param.c22527)

    costs.acc2252()

    assert cost_variables.c22521 == pytest.approx(acc2252param.expected_c22521)

    assert cost_variables.c22522 == pytest.approx(acc2252param.expected_c22522)

    assert cost_variables.c22523 == pytest.approx(acc2252param.expected_c22523)

    assert cost_variables.c22524 == pytest.approx(acc2252param.expected_c22524)


class Acc2253Param(NamedTuple):
    ucblss: Any = None

    fkind: Any = None

    p_plant_primary_heat_mw: Any = None

    p_plant_electric_net_mw: Any = None

    i_pulsed_plant: Any = None

    dtstor: Any = None

    istore: Any = None

    t_plant_pulse_no_burn: Any = None

    c22: Any = None

    c225: Any = None

    c2253: Any = None

    expected_c2253: Any = None


@pytest.mark.parametrize(
    "acc2253param",
    (
        Acc2253Param(
            ucblss=90,
            fkind=1,
            p_plant_primary_heat_mw=2620.2218111502593,
            p_plant_electric_net_mw=493.01760776192009,
            i_pulsed_plant=1,
            dtstor=300,
            istore=1,
            t_plant_pulse_no_burn=854.42613938735622,
            c22=0,
            c225=0,
            c2253=0,
            expected_c2253=20.785622343242554,
        ),
        Acc2253Param(
            ucblss=90,
            fkind=1,
            p_plant_primary_heat_mw=2619.4223856129224,
            p_plant_electric_net_mw=422.4198205312706,
            i_pulsed_plant=1,
            dtstor=300,
            istore=1,
            t_plant_pulse_no_burn=854.42613938735622,
            c22=3474.7391916096453,
            c225=185.05656643685359,
            c2253=20.785622343242554,
            expected_c2253=17.809219633598371,
        ),
    ),
)
def test_acc2253(acc2253param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2253.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2253param: the data used to mock and assert in this test.
    :type acc2253param: acc2253param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucblss", acc2253param.ucblss)

    monkeypatch.setattr(cost_variables, "fkind", acc2253param.fkind)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_primary_heat_mw",
        acc2253param.p_plant_primary_heat_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_net_mw",
        acc2253param.p_plant_electric_net_mw,
    )

    monkeypatch.setattr(pulse_variables, "i_pulsed_plant", acc2253param.i_pulsed_plant)

    monkeypatch.setattr(pulse_variables, "dtstor", acc2253param.dtstor)

    monkeypatch.setattr(pulse_variables, "istore", acc2253param.istore)

    monkeypatch.setattr(
        times_variables, "t_plant_pulse_no_burn", acc2253param.t_plant_pulse_no_burn
    )

    monkeypatch.setattr(cost_variables, "c22", acc2253param.c22)

    monkeypatch.setattr(cost_variables, "c225", acc2253param.c225)

    monkeypatch.setattr(cost_variables, "c2253", acc2253param.c2253)

    costs.acc2253()

    assert cost_variables.c2253 == pytest.approx(acc2253param.expected_c2253)


class Acc226Param(NamedTuple):
    c226: Any = None

    c2261: Any = None

    c2262: Any = None

    c2263: Any = None

    c22: Any = None

    expected_c226: Any = None


@pytest.mark.parametrize(
    "acc226param",
    (
        Acc226Param(
            c226=0,
            c2261=85.82488824875719,
            c2262=20.313088941037051,
            c2263=122.17123799205466,
            c22=0,
            expected_c226=228.30921518184891,
        ),
        Acc226Param(
            c226=228.30921518184891,
            c2261=86.412964519098367,
            c2262=25.118525150548585,
            c2263=247.55533515524576,
            c22=3474.7391916096453,
            expected_c226=359.08682482489269,
        ),
    ),
)
def test_acc226(acc226param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc226.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc226param: the data used to mock and assert in this test.
    :type acc226param: acc226param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "c226", acc226param.c226)

    monkeypatch.setattr(cost_variables, "c2261", acc226param.c2261)

    monkeypatch.setattr(cost_variables, "c2262", acc226param.c2262)

    monkeypatch.setattr(cost_variables, "c2263", acc226param.c2263)

    monkeypatch.setattr(cost_variables, "c22", acc226param.c22)

    costs.acc226()

    assert cost_variables.c226 == pytest.approx(acc226param.expected_c226)


class Acc2261Param(NamedTuple):
    uchts: Any = None

    lsa: Any = None

    fkind: Any = None

    i_blkt_coolant_type: Any = None

    p_shld_nuclear_heat_mw: Any = None

    p_blkt_nuclear_heat_total_mw: Any = None

    p_plant_primary_heat_mw: Any = None

    p_fw_div_heat_deposited_mw: Any = None

    n_primary_heat_exchangers: Any = None

    c226: Any = None

    c2261: Any = None

    c22: Any = None

    chx: Any = None

    cpp: Any = None

    expected_c2261: Any = None

    expected_chx: Any = None

    expected_cpp: Any = None


@pytest.mark.parametrize(
    "acc2261param",
    (
        Acc2261Param(
            uchts=np.array(
                np.array((15.300000000000001, 19.100000000000001), order="F"),
                order="F",
            ).transpose(),
            lsa=2,
            fkind=1,
            i_blkt_coolant_type=1,
            p_shld_nuclear_heat_mw=1.3609360176065353,
            p_blkt_nuclear_heat_total_mw=1504.711566619962,
            p_plant_primary_heat_mw=2620.2218111502593,
            p_fw_div_heat_deposited_mw=0,
            n_primary_heat_exchangers=3,
            c226=0,
            c2261=0,
            c22=0,
            chx=0,
            cpp=0,
            expected_c2261=85.82488824875719,
            expected_chx=57.169226428813381,
            expected_cpp=28.655661819943806,
        ),
        Acc2261Param(
            uchts=np.array(
                np.array((15.300000000000001, 19.100000000000001), order="F"),
                order="F",
            ).transpose(),
            lsa=2,
            fkind=1,
            i_blkt_coolant_type=1,
            p_shld_nuclear_heat_mw=1.4036212304705389,
            p_blkt_nuclear_heat_total_mw=1549.9285082739402,
            p_plant_primary_heat_mw=2619.4223856129224,
            p_fw_div_heat_deposited_mw=0,
            n_primary_heat_exchangers=3,
            c226=228.30921518184891,
            c2261=85.82488824875719,
            c22=3474.7391916096453,
            chx=57.169226428813381,
            cpp=28.655661819943806,
            expected_c2261=86.412964519098367,
            expected_chx=57.157016301470911,
            expected_cpp=29.255948217627452,
        ),
    ),
)
def test_acc2261_rut(acc2261param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2261.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2261param: the data used to mock and assert in this test.
    :type acc2261param: acc2261param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "uchts", acc2261param.uchts)

    monkeypatch.setattr(cost_variables, "lsa", acc2261param.lsa)

    monkeypatch.setattr(cost_variables, "fkind", acc2261param.fkind)

    monkeypatch.setattr(
        fwbs_variables, "i_blkt_coolant_type", acc2261param.i_blkt_coolant_type
    )

    monkeypatch.setattr(
        fwbs_variables, "p_shld_nuclear_heat_mw", acc2261param.p_shld_nuclear_heat_mw
    )

    monkeypatch.setattr(
        fwbs_variables,
        "p_blkt_nuclear_heat_total_mw",
        acc2261param.p_blkt_nuclear_heat_total_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_primary_heat_mw",
        acc2261param.p_plant_primary_heat_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_fw_div_heat_deposited_mw",
        acc2261param.p_fw_div_heat_deposited_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "n_primary_heat_exchangers",
        acc2261param.n_primary_heat_exchangers,
    )

    monkeypatch.setattr(cost_variables, "c226", acc2261param.c226)

    monkeypatch.setattr(cost_variables, "c2261", acc2261param.c2261)

    monkeypatch.setattr(cost_variables, "c22", acc2261param.c22)

    monkeypatch.setattr(cost_variables, "chx", acc2261param.chx)

    monkeypatch.setattr(cost_variables, "cpp", acc2261param.cpp)

    costs.acc2261()

    assert cost_variables.c2261 == pytest.approx(acc2261param.expected_c2261)

    assert cost_variables.chx == pytest.approx(acc2261param.expected_chx)

    assert cost_variables.cpp == pytest.approx(acc2261param.expected_cpp)


class Acc2262Param(NamedTuple):
    lsa: Any = None

    fkind: Any = None

    tfacmw: Any = None

    ife: Any = None

    tdspmw: Any = None

    p_hcd_electric_loss_mw: Any = None

    vachtmw: Any = None

    p_tritium_plant_electric_mw: Any = None

    fachtmw: Any = None

    p_cryo_plant_electric_mw: Any = None

    c226: Any = None

    c2262: Any = None

    c22: Any = None

    cpp: Any = None

    cppa: Any = None

    expected_c2262: Any = None

    expected_cppa: Any = None


@pytest.mark.parametrize(
    "acc2262param",
    (
        Acc2262Param(
            lsa=2,
            fkind=1,
            tfacmw=0,
            ife=0,
            tdspmw=0.01,
            p_hcd_electric_loss_mw=77.967671580642758,
            vachtmw=0.5,
            p_tritium_plant_electric_mw=15,
            fachtmw=61.882833632875375,
            p_cryo_plant_electric_mw=37.900388528497025,
            c226=0,
            c2262=0,
            c22=0,
            cpp=28.655661819943806,
            cppa=0,
            expected_c2262=20.313088941037051,
            expected_cppa=20.313088941037051,
        ),
        Acc2262Param(
            lsa=2,
            fkind=1,
            tfacmw=0,
            ife=0,
            tdspmw=0.01,
            p_hcd_electric_loss_mw=77.967671580642758,
            vachtmw=0.5,
            p_tritium_plant_electric_mw=15,
            fachtmw=62.237143915360818,
            p_cryo_plant_electric_mw=108.74512702403499,
            c226=228.30921518184891,
            c2262=20.313088941037051,
            c22=3474.7391916096453,
            cpp=29.255948217627452,
            cppa=20.313088941037051,
            expected_c2262=25.118525150548585,
            expected_cppa=25.118525150548585,
        ),
    ),
)
def test_acc2262_rut(acc2262param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2262.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2262param: the data used to mock and assert in this test.
    :type acc2262param: acc2262param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "lsa", acc2262param.lsa)

    monkeypatch.setattr(cost_variables, "fkind", acc2262param.fkind)

    monkeypatch.setattr(ife_variables, "tfacmw", acc2262param.tfacmw)

    monkeypatch.setattr(ife_variables, "ife", acc2262param.ife)

    monkeypatch.setattr(ife_variables, "tdspmw", acc2262param.tdspmw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_loss_mw",
        acc2262param.p_hcd_electric_loss_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "vachtmw", acc2262param.vachtmw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_tritium_plant_electric_mw",
        acc2262param.p_tritium_plant_electric_mw,
    )

    monkeypatch.setattr(heat_transport_variables, "fachtmw", acc2262param.fachtmw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_cryo_plant_electric_mw",
        acc2262param.p_cryo_plant_electric_mw,
    )

    monkeypatch.setattr(cost_variables, "c226", acc2262param.c226)

    monkeypatch.setattr(cost_variables, "c2262", acc2262param.c2262)

    monkeypatch.setattr(cost_variables, "c22", acc2262param.c22)

    monkeypatch.setattr(cost_variables, "cpp", acc2262param.cpp)

    monkeypatch.setattr(cost_variables, "cppa", acc2262param.cppa)

    costs.acc2262()

    assert cost_variables.c2262 == pytest.approx(acc2262param.expected_c2262)

    assert cost_variables.cppa == pytest.approx(acc2262param.expected_cppa)


class Acc2263Param(NamedTuple):
    uccry: Any = None

    lsa: Any = None

    fkind: Any = None

    helpow: Any = None

    temp_tf_cryo: Any = None

    c226: Any = None

    c2263: Any = None

    c22: Any = None

    expected_c2263: Any = None


@pytest.mark.parametrize(
    "acc2263param",
    (
        Acc2263Param(
            uccry=93000,
            lsa=2,
            fkind=1,
            helpow=76851.741036987034,
            temp_tf_cryo=4.5,
            c226=0,
            c2263=0,
            c22=0,
            expected_c2263=122.17123799205466,
        ),
        Acc2263Param(
            uccry=93000,
            lsa=2,
            fkind=1,
            helpow=220505.71684249729,
            temp_tf_cryo=4.5,
            c226=228.30921518184891,
            c2263=122.17123799205466,
            c22=3474.7391916096453,
            expected_c2263=247.55533515524576,
        ),
    ),
)
def test_acc2263_rut(acc2263param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2263.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2263param: the data used to mock and assert in this test.
    :type acc2263param: acc2263param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "uccry", acc2263param.uccry)

    monkeypatch.setattr(cost_variables, "lsa", acc2263param.lsa)

    monkeypatch.setattr(cost_variables, "fkind", acc2263param.fkind)

    monkeypatch.setattr(heat_transport_variables, "helpow", acc2263param.helpow)

    monkeypatch.setattr(tfcoil_variables, "temp_tf_cryo", acc2263param.temp_tf_cryo)

    monkeypatch.setattr(cost_variables, "c226", acc2263param.c226)

    monkeypatch.setattr(cost_variables, "c2263", acc2263param.c2263)

    monkeypatch.setattr(cost_variables, "c22", acc2263param.c22)

    costs.acc2263()

    assert cost_variables.c2263 == pytest.approx(acc2263param.expected_c2263)


class Acc227Param(NamedTuple):
    c227: Any = None

    c2271: Any = None

    c2272: Any = None

    c2273: Any = None

    c2274: Any = None

    c22: Any = None

    expected_c227: Any = None


@pytest.mark.parametrize(
    "acc227param",
    (
        Acc227Param(
            c227=0,
            c2271=22.300000000000001,
            c2272=114.02873340990777,
            c2273=69.115208498727412,
            c2274=79.525098581749191,
            c22=0,
            expected_c227=284.96904049038437,
        ),
        Acc227Param(
            c227=284.96904049038437,
            c2271=22.300000000000001,
            c2272=114.00948752346841,
            c2273=69.202425860597359,
            c2274=79.60537144364551,
            c22=3474.7391916096453,
            expected_c227=285.11728482771127,
        ),
    ),
)
def test_acc227(acc227param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc227.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc227param: the data used to mock and assert in this test.
    :type acc227param: acc227param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "c227", acc227param.c227)

    monkeypatch.setattr(cost_variables, "c2271", acc227param.c2271)

    monkeypatch.setattr(cost_variables, "c2272", acc227param.c2272)

    monkeypatch.setattr(cost_variables, "c2273", acc227param.c2273)

    monkeypatch.setattr(cost_variables, "c2274", acc227param.c2274)

    monkeypatch.setattr(cost_variables, "c22", acc227param.c22)

    costs.acc227()

    assert cost_variables.c227 == pytest.approx(acc227param.expected_c227)


class Acc2271Param(NamedTuple):
    ucf1: Any = None

    fkind: Any = None

    c227: Any = None

    c2271: Any = None

    c22: Any = None

    expected_c2271: Any = None


@pytest.mark.parametrize(
    "acc2271param",
    (
        Acc2271Param(
            ucf1=22300000,
            fkind=1,
            c227=0,
            c2271=0,
            c22=0,
            expected_c2271=22.300000000000001,
        ),
        Acc2271Param(
            ucf1=22300000,
            fkind=1,
            c227=284.96904049038437,
            c2271=22.300000000000001,
            c22=3474.7391916096453,
            expected_c2271=22.300000000000001,
        ),
    ),
)
def test_acc2271_rut(acc2271param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2271.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2271param: the data used to mock and assert in this test.
    :type acc2271param: acc2271param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucf1", acc2271param.ucf1)

    monkeypatch.setattr(cost_variables, "fkind", acc2271param.fkind)

    monkeypatch.setattr(cost_variables, "c227", acc2271param.c227)

    monkeypatch.setattr(cost_variables, "c2271", acc2271param.c2271)

    monkeypatch.setattr(cost_variables, "c22", acc2271param.c22)

    costs.acc2271()

    assert cost_variables.c2271 == pytest.approx(acc2271param.expected_c2271)


class Acc2272Param(NamedTuple):
    fkind: Any = None

    fburn: Any = None

    reprat: Any = None

    ife: Any = None

    gain: Any = None

    edrive: Any = None

    wtgpd: Any = None

    rndfuel: Any = None

    m_fuel_amu: Any = None

    c227: Any = None

    c2272: Any = None

    c22: Any = None

    expected_wtgpd: Any = None

    expected_c2272: Any = None


@pytest.mark.parametrize(
    "acc2272param",
    (
        Acc2272Param(
            fkind=1,
            fburn=0.33329999999999999,
            reprat=0,
            ife=0,
            gain=0,
            edrive=5000000,
            wtgpd=0,
            rndfuel=7.0799717510383796e20,
            m_fuel_amu=2.5,
            c227=0,
            c2272=0,
            c22=0,
            expected_wtgpd=507.88376577416528,
            expected_c2272=114.02873340990777,
        ),
        Acc2272Param(
            fkind=1,
            fburn=0.33329999999999999,
            reprat=0,
            ife=0,
            gain=0,
            edrive=5000000,
            wtgpd=507.88376577416528,
            rndfuel=7.0777619721108953e20,
            m_fuel_amu=2.5,
            c227=284.96904049038437,
            c2272=114.02873340990777,
            c22=3474.7391916096453,
            expected_wtgpd=507.72524666099866,
            expected_c2272=114.00948752346841,
        ),
    ),
)
def test_acc2272_rut(acc2272param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2272.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2272param: the data used to mock and assert in this test.
    :type acc2272param: acc2272param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "fkind", acc2272param.fkind)

    monkeypatch.setattr(ife_variables, "fburn", acc2272param.fburn)

    monkeypatch.setattr(ife_variables, "reprat", acc2272param.reprat)

    monkeypatch.setattr(ife_variables, "ife", acc2272param.ife)

    monkeypatch.setattr(ife_variables, "gain", acc2272param.gain)

    monkeypatch.setattr(ife_variables, "edrive", acc2272param.edrive)

    monkeypatch.setattr(physics_variables, "wtgpd", acc2272param.wtgpd)

    monkeypatch.setattr(physics_variables, "rndfuel", acc2272param.rndfuel)

    monkeypatch.setattr(physics_variables, "m_fuel_amu", acc2272param.m_fuel_amu)

    monkeypatch.setattr(cost_variables, "c227", acc2272param.c227)

    monkeypatch.setattr(cost_variables, "c2272", acc2272param.c2272)

    monkeypatch.setattr(cost_variables, "c22", acc2272param.c22)

    costs.acc2272()

    assert physics_variables.wtgpd == pytest.approx(acc2272param.expected_wtgpd)

    assert cost_variables.c2272 == pytest.approx(acc2272param.expected_c2272)


class Acc2273Param(NamedTuple):
    wsvol: Any = None

    volrci: Any = None

    fkind: Any = None

    f_plasma_fuel_tritium: Any = None

    c227: Any = None

    c2273: Any = None

    c22: Any = None

    expected_c2273: Any = None


@pytest.mark.parametrize(
    "acc2273param",
    (
        Acc2273Param(
            wsvol=130018.25667917728,
            volrci=1205439.8543893537,
            fkind=1,
            f_plasma_fuel_tritium=0.5,
            c227=0,
            c2273=0,
            c22=0,
            expected_c2273=69.115208498727412,
        ),
        Acc2273Param(
            wsvol=130255.93791329287,
            volrci=1206887.4047542624,
            fkind=1,
            f_plasma_fuel_tritium=0.5,
            c227=284.96904049038437,
            c2273=69.115208498727412,
            c22=3474.7391916096453,
            expected_c2273=69.202425860597359,
        ),
    ),
)
def test_acc2273_rut(acc2273param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2273.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2273param: the data used to mock and assert in this test.
    :type acc2273param: acc2273param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(buildings_variables, "wsvol", acc2273param.wsvol)

    monkeypatch.setattr(buildings_variables, "volrci", acc2273param.volrci)

    monkeypatch.setattr(cost_variables, "fkind", acc2273param.fkind)

    monkeypatch.setattr(
        physics_variables, "f_plasma_fuel_tritium", acc2273param.f_plasma_fuel_tritium
    )

    monkeypatch.setattr(cost_variables, "c227", acc2273param.c227)

    monkeypatch.setattr(cost_variables, "c2273", acc2273param.c2273)

    monkeypatch.setattr(cost_variables, "c22", acc2273param.c22)

    costs.acc2273()

    assert cost_variables.c2273 == pytest.approx(acc2273param.expected_c2273)


class Acc2274Param(NamedTuple):
    wsvol: Any = None

    volrci: Any = None

    fkind: Any = None

    c227: Any = None

    c2274: Any = None

    c22: Any = None

    expected_c2274: Any = None


@pytest.mark.parametrize(
    "acc2274param",
    (
        Acc2274Param(
            wsvol=130018.25667917728,
            volrci=1205439.8543893537,
            fkind=1,
            c227=0,
            c2274=0,
            c22=0,
            expected_c2274=79.525098581749191,
        ),
        Acc2274Param(
            wsvol=130255.93791329287,
            volrci=1206887.4047542624,
            fkind=1,
            c227=284.96904049038437,
            c2274=79.525098581749191,
            c22=3474.7391916096453,
            expected_c2274=79.60537144364551,
        ),
    ),
)
def test_acc2274_rut(acc2274param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2274.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2274param: the data used to mock and assert in this test.
    :type acc2274param: acc2274param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(buildings_variables, "wsvol", acc2274param.wsvol)

    monkeypatch.setattr(buildings_variables, "volrci", acc2274param.volrci)

    monkeypatch.setattr(cost_variables, "fkind", acc2274param.fkind)

    monkeypatch.setattr(cost_variables, "c227", acc2274param.c227)

    monkeypatch.setattr(cost_variables, "c2274", acc2274param.c2274)

    monkeypatch.setattr(cost_variables, "c22", acc2274param.c22)

    costs.acc2274()

    assert cost_variables.c2274 == pytest.approx(acc2274param.expected_c2274)


class Acc228Param(NamedTuple):
    uciac: Any = None

    fkind: Any = None

    c228: Any = None

    c22: Any = None

    expected_c228: Any = None


@pytest.mark.parametrize(
    "acc228param",
    (
        Acc228Param(
            uciac=150000000,
            fkind=1,
            c228=0,
            c22=0,
            expected_c228=150,
        ),
        Acc228Param(
            uciac=150000000,
            fkind=1,
            c228=150,
            c22=3474.7391916096453,
            expected_c228=150,
        ),
    ),
)
def test_acc228_rut(acc228param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc228.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc228param: the data used to mock and assert in this test.
    :type acc228param: acc228param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "uciac", acc228param.uciac)

    monkeypatch.setattr(cost_variables, "fkind", acc228param.fkind)

    monkeypatch.setattr(cost_variables, "c228", acc228param.c228)

    monkeypatch.setattr(cost_variables, "c22", acc228param.c22)

    costs.acc228()

    assert cost_variables.c228 == pytest.approx(acc228param.expected_c228)


class Acc229Param(NamedTuple):
    ucme: Any = None

    fkind: Any = None

    c229: Any = None

    c22: Any = None

    expected_c229: Any = None


@pytest.mark.parametrize(
    "acc229param",
    (
        Acc229Param(
            ucme=300000000,
            fkind=1,
            c229=0,
            c22=0,
            expected_c229=300,
        ),
        Acc229Param(
            ucme=300000000,
            fkind=1,
            c229=300,
            c22=3474.7391916096453,
            expected_c229=300,
        ),
    ),
)
def test_acc229_rut(acc229param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc229.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc229param: the data used to mock and assert in this test.
    :type acc229param: acc229param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucme", acc229param.ucme)

    monkeypatch.setattr(cost_variables, "fkind", acc229param.fkind)

    monkeypatch.setattr(cost_variables, "c229", acc229param.c229)

    monkeypatch.setattr(cost_variables, "c22", acc229param.c22)

    costs.acc229()

    assert cost_variables.c229 == pytest.approx(acc229param.expected_c229)


class Acc23Param(NamedTuple):
    ucturb: Any = None

    ireactor: Any = None

    i_blkt_coolant_type: Any = None

    p_plant_electric_gross_mw: Any = None

    c23: Any = None

    expected_c23: Any = None


@pytest.mark.parametrize(
    "acc23param",
    (
        Acc23Param(
            ucturb=np.array(
                np.array((230000000, 245000000), order="F"), order="F"
            ).transpose(),
            ireactor=1,
            i_blkt_coolant_type=1,
            p_plant_electric_gross_mw=982.58317918134742,
            c23=0,
            expected_c23=194.83812507173698,
        ),
        Acc23Param(
            ucturb=np.array(
                np.array((230000000, 245000000), order="F"), order="F"
            ).transpose(),
            ireactor=1,
            i_blkt_coolant_type=1,
            p_plant_electric_gross_mw=982.28339460484608,
            c23=194.83812507173698,
            expected_c23=194.78878460447092,
        ),
    ),
)
def test_acc23_rut(acc23param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc23.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc23param: the data used to mock and assert in this test.
    :type acc23param: acc23param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucturb", acc23param.ucturb)

    monkeypatch.setattr(cost_variables, "ireactor", acc23param.ireactor)

    monkeypatch.setattr(
        fwbs_variables, "i_blkt_coolant_type", acc23param.i_blkt_coolant_type
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_gross_mw",
        acc23param.p_plant_electric_gross_mw,
    )

    monkeypatch.setattr(cost_variables, "c23", acc23param.c23)

    costs.acc23()

    assert cost_variables.c23 == pytest.approx(acc23param.expected_c23)


class Acc24Param(NamedTuple):
    c24: Any = None

    c241: Any = None

    c242: Any = None

    c243: Any = None

    c244: Any = None

    c245: Any = None

    expected_c24: Any = None


@pytest.mark.parametrize(
    "acc24param",
    (
        Acc24Param(
            c24=0,
            c241=14.443999999999999,
            c242=12.196675853540341,
            c243=10.979786178504369,
            c244=5.3380000000000001,
            c245=1.1775,
            expected_c24=44.135962032044716,
        ),
        Acc24Param(
            c24=44.135962032044716,
            c241=14.443999999999999,
            c242=7.2671621358073066,
            c243=6.4715020827802281,
            c244=5.3380000000000001,
            c245=1.1775,
            expected_c24=34.698164218587536,
        ),
    ),
)
def test_acc24(acc24param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc24.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc24param: the data used to mock and assert in this test.
    :type acc24param: acc24param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "c24", acc24param.c24)

    monkeypatch.setattr(cost_variables, "c241", acc24param.c241)

    monkeypatch.setattr(cost_variables, "c242", acc24param.c242)

    monkeypatch.setattr(cost_variables, "c243", acc24param.c243)

    monkeypatch.setattr(cost_variables, "c244", acc24param.c244)

    monkeypatch.setattr(cost_variables, "c245", acc24param.c245)

    costs.acc24()

    assert cost_variables.c24 == pytest.approx(acc24param.expected_c24)


class Acc241Param(NamedTuple):
    lsa: Any = None

    c24: Any = None

    c241: Any = None

    expected_c241: Any = None


@pytest.mark.parametrize(
    "acc241param",
    (
        Acc241Param(
            lsa=2,
            c24=0,
            c241=0,
            expected_c241=14.443999999999999,
        ),
        Acc241Param(
            lsa=2,
            c24=44.135962032044716,
            c241=14.443999999999999,
            expected_c241=14.443999999999999,
        ),
    ),
)
def test_acc241_rut(acc241param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc241.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc241param: the data used to mock and assert in this test.
    :type acc241param: acc241param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "lsa", acc241param.lsa)

    monkeypatch.setattr(cost_variables, "c24", acc241param.c24)

    monkeypatch.setattr(cost_variables, "c241", acc241param.c241)

    costs.acc241()

    assert cost_variables.c241 == pytest.approx(acc241param.expected_c241)


class Acc242Param(NamedTuple):
    lsa: Any = None

    pacpmw: Any = None

    p_plant_electric_base_total_mw: Any = None

    c24: Any = None

    c242: Any = None

    cpp: Any = None

    expected_c242: Any = None


@pytest.mark.parametrize(
    "acc242param",
    (
        Acc242Param(
            lsa=2,
            pacpmw=1226.1273281650574,
            p_plant_electric_base_total_mw=61.882833632875375,
            c24=0,
            c242=0,
            cpp=28.655661819943806,
            expected_c242=12.196675853540341,
        ),
        Acc242Param(
            lsa=2,
            pacpmw=651.53859031110449,
            p_plant_electric_base_total_mw=62.237143915360818,
            c24=44.135962032044716,
            c242=12.196675853540341,
            cpp=29.255948217627452,
            expected_c242=7.2671621358073075,
        ),
    ),
)
def test_acc242_rut(acc242param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc242.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc242param: the data used to mock and assert in this test.
    :type acc242param: acc242param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "lsa", acc242param.lsa)

    monkeypatch.setattr(heat_transport_variables, "pacpmw", acc242param.pacpmw)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_base_total_mw",
        acc242param.p_plant_electric_base_total_mw,
    )

    monkeypatch.setattr(cost_variables, "c24", acc242param.c24)

    monkeypatch.setattr(cost_variables, "c242", acc242param.c242)

    monkeypatch.setattr(cost_variables, "cpp", acc242param.cpp)

    costs.acc242()

    assert cost_variables.c242 == pytest.approx(acc242param.expected_c242)


class Acc243Param(NamedTuple):
    lsa: Any = None

    tlvpmw: Any = None

    c24: Any = None

    c243: Any = None

    expected_c243: Any = None


@pytest.mark.parametrize(
    "acc243param",
    (
        Acc243Param(
            lsa=2,
            tlvpmw=699.34943812129745,
            c24=0,
            c243=0,
            expected_c243=10.979786178504369,
        ),
        Acc243Param(
            lsa=2,
            tlvpmw=412.19758489046046,
            c24=44.135962032044716,
            c243=10.979786178504369,
            expected_c243=6.471502082780229,
        ),
    ),
)
def test_acc243_rut(acc243param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc243.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc243param: the data used to mock and assert in this test.
    :type acc243param: acc243param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "lsa", acc243param.lsa)

    monkeypatch.setattr(heat_transport_variables, "tlvpmw", acc243param.tlvpmw)

    monkeypatch.setattr(cost_variables, "c24", acc243param.c24)

    monkeypatch.setattr(cost_variables, "c243", acc243param.c243)

    costs.acc243()

    assert cost_variables.c243 == pytest.approx(acc243param.expected_c243)


class Acc244Param(NamedTuple):
    lsa: Any = None

    c24: Any = None

    c244: Any = None

    expected_c244: Any = None


@pytest.mark.parametrize(
    "acc244param",
    (
        Acc244Param(
            lsa=2,
            c24=0,
            c244=0,
            expected_c244=5.3380000000000001,
        ),
        Acc244Param(
            lsa=2,
            c24=44.135962032044716,
            c244=5.3380000000000001,
            expected_c244=5.3380000000000001,
        ),
    ),
)
def test_acc244_rut(acc244param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc244.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc244param: the data used to mock and assert in this test.
    :type acc244param: acc244param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "lsa", acc244param.lsa)

    monkeypatch.setattr(cost_variables, "c24", acc244param.c24)

    monkeypatch.setattr(cost_variables, "c244", acc244param.c244)

    costs.acc244()

    assert cost_variables.c244 == pytest.approx(acc244param.expected_c244)


class Acc245Param(NamedTuple):
    lsa: Any = None

    c24: Any = None

    c245: Any = None

    expected_c245: Any = None


@pytest.mark.parametrize(
    "acc245param",
    (
        Acc245Param(
            lsa=2,
            c24=0,
            c245=0,
            expected_c245=1.1775,
        ),
        Acc245Param(
            lsa=2,
            c24=44.135962032044716,
            c245=1.1775,
            expected_c245=1.1775,
        ),
    ),
)
def test_acc245_rut(acc245param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc245.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc245param: the data used to mock and assert in this test.
    :type acc245param: acc245param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "lsa", acc245param.lsa)

    monkeypatch.setattr(cost_variables, "c24", acc245param.c24)

    monkeypatch.setattr(cost_variables, "c245", acc245param.c245)

    costs.acc245()

    assert cost_variables.c245 == pytest.approx(acc245param.expected_c245)


class Acc25Param(NamedTuple):
    ucmisc: Any = None

    lsa: Any = None

    c25: Any = None

    expected_c25: Any = None


@pytest.mark.parametrize(
    "acc25param",
    (
        Acc25Param(
            ucmisc=25000000,
            lsa=2,
            c25=0,
            expected_c25=22.125,
        ),
        Acc25Param(
            ucmisc=25000000,
            lsa=2,
            c25=22.125,
            expected_c25=22.125,
        ),
    ),
)
def test_acc25_rut(acc25param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc25.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc25param: the data used to mock and assert in this test.
    :type acc25param: acc25param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucmisc", acc25param.ucmisc)

    monkeypatch.setattr(cost_variables, "lsa", acc25param.lsa)

    monkeypatch.setattr(cost_variables, "c25", acc25param.c25)

    costs.acc25()

    assert cost_variables.c25 == pytest.approx(acc25param.expected_c25)


class Acc26Param(NamedTuple):
    ireactor: Any = None

    uchrs: Any = None

    lsa: Any = None

    p_plant_primary_heat_mw: Any = None

    p_hcd_electric_total_mw: Any = None

    p_plant_electric_gross_mw: Any = None

    p_fusion_total_mw: Any = None

    tfcmw: Any = None

    c26: Any = None

    expected_c26: Any = None


@pytest.mark.parametrize(
    "acc26param",
    (
        Acc26Param(
            ireactor=1,
            uchrs=87900000,
            lsa=2,
            p_plant_primary_heat_mw=2620.2218111502593,
            p_hcd_electric_total_mw=129.94611930107126,
            p_plant_electric_gross_mw=982.58317918134742,
            p_fusion_total_mw=1985.785106643267,
            tfcmw=0,
            c26=0,
            expected_c26=56.327648771765475,
        ),
        Acc26Param(
            ireactor=1,
            uchrs=87900000,
            lsa=2,
            p_plant_primary_heat_mw=2619.4223856129224,
            p_hcd_electric_total_mw=129.94611930107126,
            p_plant_electric_gross_mw=982.28339460484608,
            p_fusion_total_mw=1985.1653095257811,
            tfcmw=0,
            c26=56.327648771765475,
            expected_c26=56.310463295064743,
        ),
    ),
)
def test_acc26_rut(acc26param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc26.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc26param: the data used to mock and assert in this test.
    :type acc26param: acc26param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ireactor", acc26param.ireactor)

    monkeypatch.setattr(cost_variables, "uchrs", acc26param.uchrs)

    monkeypatch.setattr(cost_variables, "lsa", acc26param.lsa)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_primary_heat_mw",
        acc26param.p_plant_primary_heat_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_total_mw",
        acc26param.p_hcd_electric_total_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_gross_mw",
        acc26param.p_plant_electric_gross_mw,
    )

    monkeypatch.setattr(
        physics_variables,
        "p_fusion_total_mw",
        acc26param.p_fusion_total_mw,
    )

    monkeypatch.setattr(tfcoil_variables, "tfcmw", acc26param.tfcmw)

    monkeypatch.setattr(cost_variables, "c26", acc26param.c26)

    costs.acc26()

    assert cost_variables.c26 == pytest.approx(acc26param.expected_c26)


class Acc9Param(NamedTuple):
    fcontng: Any = None

    lsa: Any = None

    cowner: Any = None

    cdirt: Any = None

    cfind: Any = None

    cindrt: Any = None

    ccont: Any = None

    expected_cindrt: Any = None

    expected_ccont: Any = None


@pytest.mark.parametrize(
    "acc9param",
    (
        Acc9Param(
            fcontng=0.15000000000000002,
            lsa=2,
            cowner=0.14999999999999999,
            cdirt=4532.1724050055554,
            cfind=np.array(
                np.array(
                    (
                        0.24399999999999999,
                        0.24399999999999999,
                        0.24399999999999999,
                        0.28999999999999998,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            cindrt=0,
            ccont=0,
            expected_cindrt=1271.7275768445588,
            expected_ccont=870.58499727751723,
        ),
        Acc9Param(
            fcontng=0.15000000000000002,
            lsa=2,
            cowner=0.14999999999999999,
            cdirt=4641.9862239386794,
            cfind=np.array(
                np.array(
                    (
                        0.24399999999999999,
                        0.24399999999999999,
                        0.24399999999999999,
                        0.28999999999999998,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            cindrt=1271.7275768445588,
            ccont=870.58499727751723,
            expected_cindrt=1302.5413344371934,
            expected_ccont=891.67913375638113,
        ),
    ),
)
def test_acc9_rut(acc9param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc9.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc9param: the data used to mock and assert in this test.
    :type acc9param: acc9param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "fcontng", acc9param.fcontng)

    monkeypatch.setattr(cost_variables, "lsa", acc9param.lsa)

    monkeypatch.setattr(cost_variables, "cowner", acc9param.cowner)

    monkeypatch.setattr(cost_variables, "cdirt", acc9param.cdirt)

    monkeypatch.setattr(cost_variables, "cfind", acc9param.cfind)

    monkeypatch.setattr(cost_variables, "cindrt", acc9param.cindrt)

    monkeypatch.setattr(cost_variables, "ccont", acc9param.ccont)

    costs.acc9()

    assert cost_variables.cindrt == pytest.approx(acc9param.expected_cindrt)

    assert cost_variables.ccont == pytest.approx(acc9param.expected_ccont)


class Acc2253Param(NamedTuple):
    ucblss: Any = None

    fkind: Any = None

    p_plant_primary_heat_mw: Any = None

    p_plant_electric_net_mw: Any = None

    i_pulsed_plant: Any = None

    dtstor: Any = None

    istore: Any = None

    t_plant_pulse_no_burn: Any = None

    c22: Any = None

    c225: Any = None

    c2253: Any = None

    expected_c2253: Any = None


@pytest.mark.parametrize(
    "acc2253param",
    (
        Acc2253Param(
            ucblss=90,
            fkind=1,
            p_plant_primary_heat_mw=2620.2218111502593,
            p_plant_electric_net_mw=493.01760776192009,
            i_pulsed_plant=1,
            dtstor=300,
            istore=1,
            t_plant_pulse_no_burn=854.42613938735622,
            c22=0,
            c225=0,
            c2253=0,
            expected_c2253=20.785622343242554,
        ),
        Acc2253Param(
            ucblss=90,
            fkind=1,
            p_plant_primary_heat_mw=2619.4223856129224,
            p_plant_electric_net_mw=422.4198205312706,
            i_pulsed_plant=1,
            dtstor=300,
            istore=1,
            t_plant_pulse_no_burn=854.42613938735622,
            c22=3474.7391916096453,
            c225=185.05656643685359,
            c2253=20.785622343242554,
            expected_c2253=17.809219633598371,
        ),
    ),
)
def test_acc2253_urt(acc2253param, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for acc2253.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param acc2253param: the data used to mock and assert in this test.
    :type acc2253param: acc2253param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "ucblss", acc2253param.ucblss)

    monkeypatch.setattr(cost_variables, "fkind", acc2253param.fkind)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_primary_heat_mw",
        acc2253param.p_plant_primary_heat_mw,
    )

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_net_mw",
        acc2253param.p_plant_electric_net_mw,
    )

    monkeypatch.setattr(pulse_variables, "i_pulsed_plant", acc2253param.i_pulsed_plant)

    monkeypatch.setattr(pulse_variables, "dtstor", acc2253param.dtstor)

    monkeypatch.setattr(pulse_variables, "istore", acc2253param.istore)

    monkeypatch.setattr(
        times_variables, "t_plant_pulse_no_burn", acc2253param.t_plant_pulse_no_burn
    )

    monkeypatch.setattr(cost_variables, "c22", acc2253param.c22)

    monkeypatch.setattr(cost_variables, "c225", acc2253param.c225)

    monkeypatch.setattr(cost_variables, "c2253", acc2253param.c2253)

    costs.acc2253()

    assert cost_variables.c2253 == pytest.approx(acc2253param.expected_c2253)


class CoelcParam(NamedTuple):
    fcdfuel: Any = None

    uche3: Any = None

    life_plant: Any = None

    ifueltyp: Any = None

    cpstcst: Any = None

    coeoam: Any = None

    coecap: Any = None

    output_costs: Any = None

    coe: Any = None

    lsa: Any = None

    f_t_plant_available: Any = None

    divcst: Any = None

    ucfuel: Any = None

    divlife: Any = None

    divlife_cal: Any = None

    coefuelt: Any = None

    moneyint: Any = None

    life_hcd_fpy: Any = None

    cdrlife_cal: Any = None

    capcost: Any = None

    cplife: Any = None

    cplife_cal: Any = None

    fwallcst: Any = None

    fcr0: Any = None

    discount_rate: Any = None

    decomf: Any = None

    cdcost: Any = None

    fcap0: Any = None

    fcap0cp: Any = None

    ucwst: Any = None

    ucoam: Any = None

    dtlife: Any = None

    blkcst: Any = None

    dintrt: Any = None

    concost: Any = None

    cfind: Any = None

    life_blkt_fpy: Any = None

    life_blkt: Any = None

    uctarg: Any = None

    ife: Any = None

    reprat: Any = None

    p_plant_electric_net_mw: Any = None

    itart: Any = None

    wtgpd: Any = None

    f_plasma_fuel_helium3: Any = None

    t_plant_pulse_total: Any = None

    t_plant_pulse_burn: Any = None

    outfile: Any = None

    expected_coeoam: Any = None

    expected_coecap: Any = None

    expected_coe: Any = None

    expected_coefuelt: Any = None

    expected_moneyint: Any = None

    expected_capcost: Any = None


@pytest.mark.parametrize(
    "coelcparam",
    (
        CoelcParam(
            fcdfuel=0.10000000000000001,
            uche3=1000000,
            life_plant=40,
            ifueltyp=1,
            cpstcst=0,
            coeoam=0,
            coecap=0,
            output_costs=0,
            coe=0,
            lsa=2,
            f_t_plant_available=0.75000000000000011,
            divcst=88.904644548525795,
            ucfuel=3.4500000000000002,
            divlife=6.1337250397740126,
            divlife_cal=6.1337250397740126,
            coefuelt=0,
            moneyint=0,
            life_hcd_fpy=19.216116010620578,
            cdrlife_cal=19.216116010620578,
            capcost=0,
            cplife=0,
            cplife_cal=0,
            fwallcst=143.19827300247195,
            fcr0=0.065000000000000016,
            discount_rate=0.060000000000000012,
            decomf=0.10000000000000001,
            cdcost=140.341808845157,
            fcap0=1.1499999999999999,
            fcap0cp=1.0600000000000001,
            ucwst=np.array(
                np.array(
                    (0, 3.9399999999999999, 5.9100000000000001, 7.8799999999999999),
                    order="F",
                ),
                order="F",
            ).transpose(),
            ucoam=np.array(
                np.array(
                    (
                        68.799999999999997,
                        68.799999999999997,
                        68.799999999999997,
                        74.400000000000006,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            dtlife=0,
            blkcst=868.59838754004318,
            dintrt=0,
            concost=6674.484979127632,
            cfind=np.array(
                np.array(
                    (
                        0.24399999999999999,
                        0.24399999999999999,
                        0.24399999999999999,
                        0.28999999999999998,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            life_blkt_fpy=19.216116010620578,
            life_blkt=19.216116010620578,
            uctarg=0.29999999999999999,
            ife=0,
            reprat=0,
            p_plant_electric_net_mw=493.01760776192009,
            itart=0,
            wtgpd=507.88376577416528,
            f_plasma_fuel_helium3=0,
            t_plant_pulse_total=10864.426139387357,
            t_plant_pulse_burn=0,
            outfile=11,
            expected_coeoam=4.4099029328740929e20,
            expected_coecap=4.9891775218979061e21,
            expected_coe=6.95253391e21,
            expected_coefuelt=1.48018708e21,
            expected_moneyint=1001.1727468691442,
            expected_capcost=7675.6577259967762,
        ),
        CoelcParam(
            fcdfuel=0.10000000000000001,
            uche3=1000000,
            life_plant=40,
            ifueltyp=1,
            cpstcst=0,
            coeoam=4.4099029328740929e20,
            coecap=4.9891775218979061e21,
            output_costs=0,
            coe=6.9525339143363677e21,
            lsa=2,
            f_t_plant_available=0.75000000000000011,
            divcst=88.904644548525795,
            ucfuel=3.4500000000000002,
            divlife=6.145510750914414,
            divlife_cal=6.145510750914414,
            coefuelt=1.4801870771036603e21,
            moneyint=1001.1727468691442,
            life_hcd_fpy=19.222115557991025,
            cdrlife_cal=19.222115557991025,
            capcost=7675.6577259967762,
            cplife=0,
            cplife_cal=0,
            fwallcst=167.7865317453867,
            fcr0=0.065000000000000016,
            discount_rate=0.060000000000000012,
            decomf=0.10000000000000001,
            cdcost=140.341808845157,
            fcap0=1.1499999999999999,
            fcap0cp=1.0600000000000001,
            ucwst=np.array(
                np.array(
                    (0, 3.9399999999999999, 5.9100000000000001, 7.8799999999999999),
                    order="F",
                ),
                order="F",
            ).transpose(),
            ucoam=np.array(
                np.array(
                    (
                        68.799999999999997,
                        68.799999999999997,
                        68.799999999999997,
                        74.400000000000006,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            dtlife=0,
            blkcst=870.20508315783786,
            dintrt=0,
            concost=6836.2066921322539,
            cfind=np.array(
                np.array(
                    (
                        0.24399999999999999,
                        0.24399999999999999,
                        0.24399999999999999,
                        0.28999999999999998,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            life_blkt_fpy=19.222115557991025,
            life_blkt=19.222115557991025,
            uctarg=0.29999999999999999,
            ife=0,
            reprat=0,
            p_plant_electric_net_mw=422.4198205312706,
            itart=0,
            wtgpd=507.72524666099866,
            f_plasma_fuel_helium3=0,
            t_plant_pulse_total=864.42613938735622,
            t_plant_pulse_burn=10230.533336387549,
            outfile=11,
            expected_coeoam=1.2419424614419636,
            expected_coecap=15.547404530833255,
            expected_coe=21.50420973,
            expected_coefuelt=4.58342338,
            expected_moneyint=1025.4310038198375,
            expected_capcost=7861.6376959520912,
        ),
    ),
)
def test_coelc(coelcparam, monkeypatch, costs):
    """
    Automatically generated Regression Unit Test for coelc.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param coelcparam: the data used to mock and assert in this test.
    :type coelcparam: coelcparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(cost_variables, "fcdfuel", coelcparam.fcdfuel)

    monkeypatch.setattr(cost_variables, "uche3", coelcparam.uche3)

    monkeypatch.setattr(cost_variables, "life_plant", coelcparam.life_plant)

    monkeypatch.setattr(cost_variables, "ifueltyp", coelcparam.ifueltyp)

    monkeypatch.setattr(cost_variables, "cpstcst", coelcparam.cpstcst)

    monkeypatch.setattr(cost_variables, "coeoam", coelcparam.coeoam)

    monkeypatch.setattr(cost_variables, "coecap", coelcparam.coecap)

    monkeypatch.setattr(cost_variables, "output_costs", coelcparam.output_costs)

    monkeypatch.setattr(cost_variables, "coe", coelcparam.coe)

    monkeypatch.setattr(cost_variables, "lsa", coelcparam.lsa)

    monkeypatch.setattr(
        cost_variables, "f_t_plant_available", coelcparam.f_t_plant_available
    )

    monkeypatch.setattr(cost_variables, "divcst", coelcparam.divcst)

    monkeypatch.setattr(cost_variables, "ucfuel", coelcparam.ucfuel)

    monkeypatch.setattr(cost_variables, "divlife", coelcparam.divlife)

    monkeypatch.setattr(cost_variables, "divlife_cal", coelcparam.divlife_cal)

    monkeypatch.setattr(cost_variables, "coefuelt", coelcparam.coefuelt)

    monkeypatch.setattr(cost_variables, "moneyint", coelcparam.moneyint)

    monkeypatch.setattr(cost_variables, "life_hcd_fpy", coelcparam.life_hcd_fpy)

    monkeypatch.setattr(cost_variables, "cdrlife_cal", coelcparam.cdrlife_cal)

    monkeypatch.setattr(cost_variables, "capcost", coelcparam.capcost)

    monkeypatch.setattr(cost_variables, "cplife", coelcparam.cplife)

    monkeypatch.setattr(cost_variables, "cplife_cal", coelcparam.cplife_cal)

    monkeypatch.setattr(cost_variables, "fwallcst", coelcparam.fwallcst)

    monkeypatch.setattr(cost_variables, "fcr0", coelcparam.fcr0)

    monkeypatch.setattr(cost_variables, "discount_rate", coelcparam.discount_rate)

    monkeypatch.setattr(cost_variables, "decomf", coelcparam.decomf)

    monkeypatch.setattr(cost_variables, "cdcost", coelcparam.cdcost)

    monkeypatch.setattr(cost_variables, "fcap0", coelcparam.fcap0)

    monkeypatch.setattr(cost_variables, "fcap0cp", coelcparam.fcap0cp)

    monkeypatch.setattr(cost_variables, "ucwst", coelcparam.ucwst)

    monkeypatch.setattr(cost_variables, "ucoam", coelcparam.ucoam)

    monkeypatch.setattr(cost_variables, "dtlife", coelcparam.dtlife)

    monkeypatch.setattr(cost_variables, "blkcst", coelcparam.blkcst)

    monkeypatch.setattr(cost_variables, "dintrt", coelcparam.dintrt)

    monkeypatch.setattr(cost_variables, "concost", coelcparam.concost)

    monkeypatch.setattr(cost_variables, "cfind", coelcparam.cfind)

    monkeypatch.setattr(fwbs_variables, "life_blkt_fpy", coelcparam.life_blkt_fpy)

    monkeypatch.setattr(fwbs_variables, "life_blkt", coelcparam.life_blkt)

    monkeypatch.setattr(ife_variables, "uctarg", coelcparam.uctarg)

    monkeypatch.setattr(ife_variables, "ife", coelcparam.ife)

    monkeypatch.setattr(ife_variables, "reprat", coelcparam.reprat)

    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_net_mw",
        coelcparam.p_plant_electric_net_mw,
    )

    monkeypatch.setattr(physics_variables, "itart", coelcparam.itart)

    monkeypatch.setattr(physics_variables, "wtgpd", coelcparam.wtgpd)

    monkeypatch.setattr(
        physics_variables, "f_plasma_fuel_helium3", coelcparam.f_plasma_fuel_helium3
    )

    monkeypatch.setattr(
        times_variables, "t_plant_pulse_total", coelcparam.t_plant_pulse_total
    )

    monkeypatch.setattr(
        times_variables, "t_plant_pulse_burn", coelcparam.t_plant_pulse_burn
    )

    costs.coelc()

    assert cost_variables.coeoam == pytest.approx(coelcparam.expected_coeoam)

    assert cost_variables.coecap == pytest.approx(coelcparam.expected_coecap)

    assert cost_variables.coe == pytest.approx(coelcparam.expected_coe)

    assert cost_variables.coefuelt == pytest.approx(coelcparam.expected_coefuelt)

    assert cost_variables.moneyint == pytest.approx(coelcparam.expected_moneyint)

    assert cost_variables.capcost == pytest.approx(coelcparam.expected_capcost)
