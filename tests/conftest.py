"""pytest config file for all tests.

Defines fixtures that will be shared across all test modules.
"""

import os
import warnings

import pytest
from _pytest.fixtures import SubRequest
from system_check import system_compatible

from process.log import logging_model_handler


def pytest_addoption(parser):
    """Add custom CLI options to pytest.

    Add regression tolerance and overwrite options to pytest.
    :param parser: pytest's CLI arg parser
    :type parser: _pytest.config.argparsing.Parser
    """
    parser.addoption(
        "--reg-tolerance",
        default=5.0,
        type=float,
        help="Percentage tolerance for regression tests",
    )
    parser.addoption(
        "--solver",
        default="vmcon",
        type=str,
        help="Solver to use in regression tests, e.g. vmcon",
    )
    parser.addoption(
        "--opt-params-only",
        action="store_true",
        default=False,
        help="Only regression test optimisation parameters: useful for solver comparisons",
    )


@pytest.fixture
def reg_tolerance(request):
    """Return the value of the --reg-tolerance CLI option.

    Gets the percentage tolerance to use in regression tests.
    e.g. "pytest --reg-tolerance=5" returns "5" here.
    :param request: request fixture to access CLI args
    :type request: SubRequest
    :return: value of --reg-tolerance option
    :rtype: float
    """
    return request.config.getoption("--reg-tolerance")


@pytest.fixture
def overwrite_refs_opt(request):
    """Return the value of the --overwrite CLI option.

    e.g. "pytest --overwrite" returns True here.
    :param request: request fixture to access CLI args
    :type request: SubRequest
    :return: True if --overwrite option present
    :rtype: bool
    """

    return request.config.getoption("--overwrite")


@pytest.fixture
def solver_name(request):
    """Return the value of the --solver CLI option.

    :param request: request fixture to access CLI args
    :type request: SubRequest
    :return: name of solver to use
    :rtype: str
    """
    return request.config.getoption("--solver")


@pytest.fixture
def opt_params_only(request: SubRequest) -> bool:
    """Return the value of the --opt-params-only CLI option.

    :param request: request fixture to access CLI args
    :type request: _pytest.fixtures.SubRequest
    :return: True if --opt-params-only option present
    :rtype: bool
    """
    return request.config.getoption("--opt-params-only")


@pytest.fixture(scope="session", autouse=True)
def precondition(request):
    """Check a user for an outdated system
    Warn a user if their system is outdated and could have
    regression test floating point issues.

    Exits the test suite if trying to overwrite tests
    to stop inaccurate test assets being written.

    e.g. "pytest --overwrite" returns True here.
    :param request: request fixture to access CLI args
    :type request: SubRequest
    """
    compatible = system_compatible()
    if compatible:
        return
    warnings.warn(
        """
        \u001b[33m\033[1mYou are running the PROCESS test suite on an outdated system.\033[0m
        This can cause floating point rounding errors in regression tests.

        Please see documentation for information on running PROCESS (and tests)
        using a Docker/Singularity container.
        """,
        UserWarning,
        stacklevel=2,
    )


@pytest.fixture(scope="session", autouse=True)
def initialise_error_module():
    """pytest fixture to initialise the error module before tests run.

    Initialise the error module initially otherwise segmentation faults can
    occur when tested subroutines raise errors.
    """
    logging_model_handler.clear_logs()


@pytest.fixture
def reinitialise_error_module():
    """Re-initialise the error module.

    If a subroutine raises an error and writes to error variables, this should
    be cleaned up when the test finishes to prevent any side-effects.

    """
    # TODO Perhaps this should be autoused by all tests? Specify use explicitly
    # for now for known error-raisers
    yield
    logging_model_handler.clear_logs()


@pytest.fixture(autouse=True)
def return_to_root():
    """Various parts of PROCESS change directories and do not always change back.
    This fixture ensures that, at the end of each test, the cwd is reset to what it
    was at the beginning of the test.
    """
    cwd = os.getcwd()
    yield
    os.chdir(cwd)
