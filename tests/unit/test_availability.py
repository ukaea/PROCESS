"""Unit tests for availability.f90."""

import pytest
from pytest import approx

from process import data_structure
from process.availability import Availability
from process.data_structure import constraint_variables as ctv
from process.data_structure import cost_variables as cv
from process.data_structure import divertor_variables as dv
from process.data_structure import fwbs_variables as fwbsv
from process.data_structure import ife_variables as ifev
from process.data_structure import physics_variables as pv
from process.data_structure import tfcoil_variables as tfv
from process.data_structure import times_variables as tv
from process.init import init_all_module_vars


@pytest.fixture
def availability():
    """Provides Availability object for testing.

    :return availability: initialised Availability object
    :type availability: process.availability.Availability
    """
    return Availability()


@pytest.mark.parametrize(
    "life_fw_fpy, ibkt_life, bktlife_exp_param",
    ((0.0000001, 0, 0.5), (0.0000001, 1, 2.5), (1.0, 0, 0.5), (1.0, 1, 1.25)),
)
def test_avail_0(monkeypatch, availability, life_fw_fpy, ibkt_life, bktlife_exp_param):
    """Test avail for i_plant_availability = 0

    :param monkeypatch: mocking fixture
    :type monkeypatch: MonkeyPatch
    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """

    # Mock module vars
    monkeypatch.setattr(ifev, "ife", 0)
    monkeypatch.setattr(pv, "p_fusion_total_mw", 4.0e3)
    monkeypatch.setattr(fwbsv, "life_fw_fpy", life_fw_fpy)
    monkeypatch.setattr(cv, "ibkt_life", ibkt_life)
    monkeypatch.setattr(cv, "abktflnc", 4.0)
    monkeypatch.setattr(pv, "pflux_fw_neutron_mw", 10.0)
    monkeypatch.setattr(cv, "tlife", 30.0)
    monkeypatch.setattr(cv, "life_dpa", 40.0)
    monkeypatch.setattr(cv, "adivflnc", 8.0)
    monkeypatch.setattr(dv, "pflux_div_heat_load_mw", 10.0)
    monkeypatch.setattr(tv, "t_plant_pulse_total", 5.0)
    monkeypatch.setattr(cv, "i_plant_availability", 0)
    monkeypatch.setattr(cv, "f_t_plant_available", 0.8)
    monkeypatch.setattr(tv, "t_plant_pulse_burn", 500.0)
    monkeypatch.setattr(pv, "itart", 1)

    availability.avail(output=False)
    cpfact_obs = cv.cpfact
    cpfact_exp = 80.0
    assert pytest.approx(cpfact_obs) == cpfact_exp

    bktlife_obs = fwbsv.life_blkt_fpy
    bktlife_exp = bktlife_exp_param
    assert pytest.approx(bktlife_obs) == bktlife_exp

    divlife_obs = cv.divlife
    divlife_exp = 1.0
    assert pytest.approx(divlife_obs) == divlife_exp

    cplife_obs = cv.cplife
    cplife_exp = 30.0
    assert pytest.approx(cplife_obs) == cplife_exp


def test_avail_1(monkeypatch, availability):
    """Test avail for i_plant_availability = 1

    :param monkeypatch: mocking fixture
    :type monkeypatch: MonkeyPatch
    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """
    # Initialise fortran variables to keep test isolated from others
    init_all_module_vars()

    # Mock module vars
    monkeypatch.setattr(cv, "i_plant_availability", 1)
    monkeypatch.setattr(cv, "divlife", 1.0)
    monkeypatch.setattr(fwbsv, "life_blkt_fpy", 7.0)
    monkeypatch.setattr(cv, "t_div_replace_yrs", 0.1)
    monkeypatch.setattr(cv, "t_blkt_replace_yrs", 0.2)
    monkeypatch.setattr(cv, "tcomrepl", 0.3)
    monkeypatch.setattr(cv, "uubop", 0.4)
    monkeypatch.setattr(cv, "uucd", 0.5)
    monkeypatch.setattr(cv, "uudiv", 0.6)
    monkeypatch.setattr(cv, "uufuel", 0.7)
    monkeypatch.setattr(cv, "uufw", 0.8)
    monkeypatch.setattr(cv, "uumag", 0.9)
    monkeypatch.setattr(cv, "uuves", 0.11)

    availability.avail(output=False)
    cfactr_obs = cv.f_t_plant_available
    cfactr_exp = 0.0006344554455445239
    assert pytest.approx(cfactr_exp) == cfactr_obs

    # Initialise fortran variables again to reset for other tests
    init_all_module_vars()


def test_calc_u_unplanned_hcd(availability):
    """Test calc_u_unplanned_hcd.

    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """
    expected = 0.02
    result = availability.calc_u_unplanned_hcd()
    assert result == expected


def test_calc_u_unplanned_bop(monkeypatch, availability):
    """Test calc_u_unplanned_bop.

    :param monkeypatch: Mock fixture
    :type monkeypatch: object

    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """

    # Mock module variable t_plant_operational_total_yrs
    monkeypatch.setattr(cv, "t_plant_operational_total_yrs", 25.0)

    # Call subroutine and check result is within an absolute tolerance
    result = availability.calc_u_unplanned_bop(output=False)
    assert result == approx(0.009, abs=0.0005)


def calc_u_planned_param(**kwargs):
    """Create a dict of parameters for a single calc_u_planned() test.

    :return: Dict of parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {
        "abktflnc": 5.0,
        "adivflnc": 10.0,
        "cpstflnc": 0.0,
        "tlife": 30.0,
        "num_rh_systems": 5,
        "pflux_fw_neutron_mw": 1.0,
        "pflux_div_heat_load_mw": 10.0,
        "itart": 0,
        "expected": approx(0.3, abs=0.05),
    }

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def calc_u_planned_params():
    """Create a list of parameter dicts for the calc_u_planned fixture.

    This consists of the nominal and nominal ST cases.

    :return: List of parameter dicts
    :rtype: list
    """
    return [
        calc_u_planned_param(),  # Nominal
        calc_u_planned_param(
            abktflnc=20.0,
            adivflnc=25.0,
            cpstflnc=20.0,
            tlife=30.0,
            num_rh_systems=4,
            pflux_fw_neutron_mw=1.0,
            pflux_div_heat_load_mw=1.0,
            itart=1,
            expected=approx(0.03, abs=0.005),
        ),  # Nominal ST
    ]


@pytest.fixture(params=calc_u_planned_params(), ids=["nominal", "ST"])
def calc_u_planned_fix(request, monkeypatch):
    """Fixture for the calc_u_planned() variables.

    This fixture is parameterised using calc_u_planned_params(). This runs the
    fixture for two different parameter sets, testing the nominal and ST cases.
    These two cases are identified using the "ids" list, which determines
    the test ID for the fixture instance with each parameter value.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of calc_u_planned() for the parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock all module variables used by calc_u_planned()
    # Some are parameterised
    monkeypatch.setattr(
        dv,
        "pflux_div_heat_load_mw",
        param["pflux_div_heat_load_mw"],
    )
    monkeypatch.setattr(data_structure.fwbs_variables, "life_blkt_fpy", 0.0)
    monkeypatch.setattr(pv, "pflux_fw_neutron_mw", param["pflux_fw_neutron_mw"])
    monkeypatch.setattr(pv, "itart", param["itart"])
    monkeypatch.setattr(cv, "tlife", param["tlife"])
    monkeypatch.setattr(cv, "divlife", 0.0)
    monkeypatch.setattr(cv, "adivflnc", param["adivflnc"])
    monkeypatch.setattr(cv, "abktflnc", param["abktflnc"])
    monkeypatch.setattr(cv, "cdrlife", 0.0)
    monkeypatch.setattr(cv, "cplife", 0.0)
    monkeypatch.setattr(cv, "cpstflnc", param["cpstflnc"])
    monkeypatch.setattr(cv, "num_rh_systems", param["num_rh_systems"])

    # Return the expected result for the given parameter list
    return param["expected"]


def test_calc_u_planned(calc_u_planned_fix, availability):
    """Test calc_u_planned.

    :param calc_u_planned_fix: Expected value of calc_u_planned()
    :type calc_u_planned_fix: ApproxScalar

    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """
    # Arguments

    # Run calc_u_planned() with the current fixture, then assert the result
    # is the expected one
    result = availability.calc_u_planned(output=False)
    assert result == calc_u_planned_fix


def calc_u_unplanned_magnets_param(**kwargs):
    """Create a parameter dict for the calc_u_unplanned_magnets fixture.

    :return: Parameter dict for the fixture
    :rtype: dict
    """
    defaults = {
        "temp_margin": 1.5,
        "temp_tf_superconductor_margin_min": 1.5,
        "temp_cs_superconductor_margin_min": 1.5,
        "t_plant_operational_total_yrs": 30,
        "conf_mag": 1.0,
        "expected": approx(0.02, abs=0.005),
    }

    # Merge default values with optional keyword args to override them
    return {**defaults, **kwargs}


def calc_u_unplanned_magnets_params():
    """Generate list of parameter dicts to parameterise the fixture.

    Parameter dicts for cases:
    no_degradation
    no_degradation_conf
    degradation_conf

    :return: List of parameter dicts
    :rtype: list
    """
    return [
        calc_u_unplanned_magnets_param(),
        calc_u_unplanned_magnets_param(
            temp_margin=2.0,
            temp_tf_superconductor_margin_min=1.6,
            temp_cs_superconductor_margin_min=1.6,
            conf_mag=0.8,
        ),
        calc_u_unplanned_magnets_param(
            temp_margin=1.8,
            temp_tf_superconductor_margin_min=1.6,
            temp_cs_superconductor_margin_min=1.6,
            conf_mag=0.8,
            expected=approx(0.03, abs=0.005),
        ),
    ]


@pytest.fixture(
    params=calc_u_unplanned_magnets_params(),
    ids=["no_degredation", "no_degradation_conf", "degradation_conf"],
)
def calc_u_unplanned_magnets_fix(request, monkeypatch):
    """Fixture for the calc_u_unplanned_magnets() test.

    :param request: Allows access to parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    :return: Expected return value of calc_u_unplanned_magnets()
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock module variables
    monkeypatch.setattr(
        cv, "t_plant_operational_total_yrs", param["t_plant_operational_total_yrs"]
    )
    monkeypatch.setattr(cv, "conf_mag", param["conf_mag"])
    monkeypatch.setattr(
        tfv,
        "temp_cs_superconductor_margin_min",
        param["temp_cs_superconductor_margin_min"],
    )
    monkeypatch.setattr(
        tfv,
        "temp_tf_superconductor_margin_min",
        param["temp_tf_superconductor_margin_min"],
    )
    monkeypatch.setattr(tfv, "temp_margin", param["temp_margin"])

    return param["expected"]


def test_calc_u_unplanned_magnets(calc_u_unplanned_magnets_fix, availability):
    """Test function for calc_u_unplanned_magnets().

    :param calc_u_unplanned_magnets_fix: Expected return value
    :type calc_u_unplanned_magnets_fix: ApproxScalar
    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """

    result = availability.calc_u_unplanned_magnets(output=False)
    assert result == calc_u_unplanned_magnets_fix


def calc_u_unplanned_divertor_param(**kwargs):
    """Make parameters for a single calc_u_unplanned_divertor() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {
        "divlife": 1.99,
        "t_plant_pulse_total": 9000,
        "expected": approx(0.02, abs=0.005),
    }

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def calc_u_unplanned_divertor_params():
    """Create a list of parameter dicts for the calc_u_planned fixture.

    Case 1: below nref
    Case 2: above nu
    Case 3: between

    :return: List of parameter dicts
    :rtype: list
    """
    return [
        calc_u_unplanned_divertor_param(),
        calc_u_unplanned_divertor_param(divlife=4, expected=approx(1, abs=0)),
        calc_u_unplanned_divertor_param(divlife=3, expected=approx(0.1, abs=0.05)),
    ]


@pytest.fixture(
    params=calc_u_unplanned_divertor_params(), ids=["below_nref", "above_nu", "between"]
)
def calc_u_unplanned_divertor_fix(request, monkeypatch):
    """Fixture for the calc_u_unplanned_divertor() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of calc_u_unplanned_divertor() for the
    parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock variables used by calc_u_unplanned_divertor()
    # Some may be parameterised
    monkeypatch.setattr(tv, "t_plant_pulse_total", param["t_plant_pulse_total"])
    monkeypatch.setattr(cv, "divlife", param["divlife"])

    # Return the expected result for the given parameter list
    return param["expected"]


def test_calc_u_unplanned_divertor(calc_u_unplanned_divertor_fix, availability):
    """Test calc_u_unplanned_divertor.

    :param calc_u_unplanned_divertor_fix: Expected value of
    calc_u_unplanned_divertor()
    :type calc_u_unplanned_divertor_fix: ApproxScalar

    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """

    # Run calc_u_unplanned_divertor() with the current fixture,
    # then assert the result is the expected one
    result = availability.calc_u_unplanned_divertor(output=False)
    assert result == calc_u_unplanned_divertor_fix


def calc_u_unplanned_fwbs_param(**kwargs):
    """Make parameters for a single calc_u_unplanned_fwbs() test.

    :return: Parameters, including expected result
    :rtype: dict
    """
    # Default parameters
    defaults = {
        "life_blkt_fpy": 5,
        "t_plant_pulse_total": 9000,
        "expected": approx(0.02, abs=0.005),
    }

    # Merge default dict with any optional keyword arguments to override values
    return {**defaults, **kwargs}


def calc_u_unplanned_fwbs_params():
    """Create a list of parameter dicts for the calc_u_planned fixture.

    Case 1: below nref
    Case 2: above nu
    Case 3: between

    :return: List of parameter dicts
    :rtype: list
    """
    return [
        calc_u_unplanned_fwbs_param(),
        calc_u_unplanned_fwbs_param(life_blkt_fpy=15, expected=approx(1, abs=0)),
        calc_u_unplanned_fwbs_param(life_blkt_fpy=8.5, expected=approx(0.1, abs=0.005)),
    ]


@pytest.fixture(
    params=calc_u_unplanned_fwbs_params(), ids=["below_nref", "above_nu", "between"]
)
def calc_u_unplanned_fwbs_fix(request, monkeypatch):
    """Fixture for the calc_u_unplanned_fwbs() variables.

    :param request: Request object for accessing parameters
    :type request: object
    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    :return: Expected return value of calc_u_unplanned_fwbs() for the parameter list
    :rtype: ApproxScalar
    """
    param = request.param

    # Mock variables used by calc_u_unplanned_fwbs()
    # Some may be parameterised
    monkeypatch.setattr(tv, "t_plant_pulse_total", param["t_plant_pulse_total"])
    monkeypatch.setattr(
        data_structure.fwbs_variables, "life_blkt_fpy", param["life_blkt_fpy"]
    )

    # Return the expected result for the given parameter list
    return param["expected"]


def test_calc_u_unplanned_fwbs(calc_u_unplanned_fwbs_fix, availability):
    """Test calc_u_unplanned_fwbs.

    :param calc_u_unplanned_fwbs_fix: Expected value of calc_u_unplanned_fwbs()
    :type calc_u_unplanned_fwbs_fix: ApproxScalar

    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """

    # Run calc_u_unplanned_fwbs() with the current fixture,
    # then assert the result is the expected one
    result = availability.calc_u_unplanned_fwbs(output=False)
    assert result == calc_u_unplanned_fwbs_fix


def test_avail_2(monkeypatch, availability):
    """Test avail_2 routine

    :param monkeypatch: Mock fixture
    :type monkeypatch: object

    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """

    # Mock return values for for functions called in avail_2
    def mock_calc_u_planned(*args, **kwargs):
        return 0.01

    def mock_calc_u_unplanned_magnets(*args, **kwargs):
        return 0.02

    def mock_calc_u_unplanned_divertor(*args, **kwargs):
        return 0.03

    def mock_calc_u_unplanned_fwbs(*args, **kwargs):
        return 0.04

    def mock_calc_u_unplanned_bop(*args, **kwargs):
        return 0.05

    def mock_calc_u_unplanned_hcd(*args, **kwargs):
        return 0.06

    def mock_calc_u_unplanned_vacuum(*args, **kwargs):
        return 0.07

    # Mock module functions
    monkeypatch.setattr(availability, "calc_u_planned", mock_calc_u_planned)
    monkeypatch.setattr(
        availability, "calc_u_unplanned_magnets", mock_calc_u_unplanned_magnets
    )
    monkeypatch.setattr(
        availability, "calc_u_unplanned_divertor", mock_calc_u_unplanned_divertor
    )
    monkeypatch.setattr(
        availability, "calc_u_unplanned_fwbs", mock_calc_u_unplanned_fwbs
    )
    monkeypatch.setattr(availability, "calc_u_unplanned_bop", mock_calc_u_unplanned_bop)
    monkeypatch.setattr(availability, "calc_u_unplanned_hcd", mock_calc_u_unplanned_hcd)
    monkeypatch.setattr(
        availability, "calc_u_unplanned_vacuum", mock_calc_u_unplanned_vacuum
    )

    # Mock module variables
    monkeypatch.setattr(tv, "t_plant_pulse_burn", 5.0)
    monkeypatch.setattr(tv, "t_plant_pulse_total", 50.0)
    monkeypatch.setattr(ifev, "ife", 0)
    monkeypatch.setattr(pv, "itart", 1)
    monkeypatch.setattr(fwbsv, "life_blkt_fpy", 5.0)
    monkeypatch.setattr(cv, "divlife", 10.0)
    monkeypatch.setattr(cv, "cplife", 15.0)

    availability.avail_2(False)

    cfactr_obs = cv.f_t_plant_available
    cfactr_exp = 0.7173
    assert pytest.approx(cfactr_obs) == cfactr_exp

    cpfact_obs = cv.cpfact
    cpfact_exp = 0.07173
    assert pytest.approx(cpfact_obs) == cpfact_exp

    bktlife_obs = fwbsv.life_blkt_fpy
    bktlife_exp = 6.97058413
    assert pytest.approx(bktlife_obs) == bktlife_exp

    divlife_obs = cv.divlife
    divlife_exp = 13.94116827
    assert pytest.approx(divlife_obs) == divlife_exp

    cplife_obs = cv.cplife
    cplife_exp = 20.9117524
    assert pytest.approx(cplife_obs) == cplife_exp


def test_avail_st(monkeypatch, availability):
    """Test avail_st routine

    :param monkeypatch: Mock fixture
    :type monkeypatch: object

    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """
    # Initialise fortran variables to keep test isolated from others
    init_all_module_vars()
    monkeypatch.setattr(cv, "tmain", 1.0)
    monkeypatch.setattr(cv, "tlife", 30.0)
    monkeypatch.setattr(cv, "u_unplanned_cp", 0.05)
    monkeypatch.setattr(tv, "t_plant_pulse_burn", 5.0)
    monkeypatch.setattr(tv, "t_plant_pulse_total", 9000.0)
    monkeypatch.setattr(cv, "adivflnc", 10.0)
    monkeypatch.setattr(dv, "pflux_div_heat_load_mw", 10.0)
    monkeypatch.setattr(cv, "ibkt_life", 0)
    monkeypatch.setattr(cv, "abktflnc", 10.0)
    monkeypatch.setattr(pv, "pflux_fw_neutron_mw", 10.0)
    monkeypatch.setattr(cv, "cplife", 5.0)
    monkeypatch.setattr(cv, "cdrlife", 15.0)

    availability.avail_st(output=False)

    assert pytest.approx(cv.t_plant_operational_total_yrs) == 15.0
    assert pytest.approx(cv.f_t_plant_available) == 0.27008858
    assert pytest.approx(cv.cpfact, abs=1.0e-8) == 0.00015005


@pytest.mark.parametrize("i_tf_sup, exp", ((1, 6.337618), (0, 4)))
def test_cp_lifetime(monkeypatch, availability, i_tf_sup, exp):
    """Test cp_lifetime routine

    :param monkeypatch: Mock fixture
    :type monkeypatch: object

    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """

    monkeypatch.setattr(tfv, "i_tf_sup", i_tf_sup)
    monkeypatch.setattr(ctv, "nflutfmax", 1.0e23)
    monkeypatch.setattr(fwbsv, "neut_flux_cp", 5.0e14)
    monkeypatch.setattr(cv, "cpstflnc", 20.0)
    monkeypatch.setattr(pv, "pflux_fw_neutron_mw", 5.0)
    monkeypatch.setattr(cv, "tlife", 30.0)

    cplife = availability.cp_lifetime()

    assert pytest.approx(cplife) == exp


def test_divertor_lifetime(monkeypatch, availability):
    """Test divertor_lifetime routine

    :param monkeypatch: Mock fixture
    :type monkeypatch: object

    :param availability: fixture containing an initialised `Availability` object
    :type availability: tests.unit.test_availability.availability (functional fixture)
    """

    monkeypatch.setattr(cv, "adivflnc", 100.0)
    monkeypatch.setattr(dv, "pflux_div_heat_load_mw", 10.0)
    monkeypatch.setattr(cv, "tlife", 30.0)

    divlife_obs = availability.divertor_lifetime()
    divlife_exp = 10.0

    assert pytest.approx(divlife_obs) == divlife_exp
