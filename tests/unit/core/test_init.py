"""Unit tests for input sanity checks in process.core.init.check_process."""

import pytest

from process.core.data_structure.base import DataStructure
from process.core.exceptions import ProcessValidationError
from process.core.init import check_process
from process.models.tfcoil.base import TFConductorModel


def _validation_error_message(data):
    """Run check_process and return any validation error message.

    Later, unrelated checks may still fire on an otherwise-default
    DataStructure, so callers assert on the message content rather than
    on whether an error was raised.
    """
    try:
        check_process(data, None)
    except ProcessValidationError as error:
        return str(error)
    return ""


def test_zero_thickness_superconducting_tf_is_rejected():
    """SC TF with dr_tf_inboard left at 0 and neither ixc=13 nor ixc=140
    active must fail validation instead of silently building a machine
    with no inboard TF coil.
    """
    data = DataStructure()
    data.tfcoil.i_tf_sup = TFConductorModel.SUPERCONDUCTING
    data.build.dr_tf_inboard = 0.0

    with pytest.raises(ProcessValidationError, match="dr_tf_inboard"):
        check_process(data, None)


def test_explicit_tf_thickness_is_accepted():
    data = DataStructure()
    data.tfcoil.i_tf_sup = TFConductorModel.SUPERCONDUCTING
    data.build.dr_tf_inboard = 0.5

    assert "dr_tf_inboard" not in _validation_error_message(data)


def test_wp_thickness_iteration_variable_is_accepted():
    """With ixc = 140 active, dr_tf_inboard is derived in the build model,
    so a zero input value is legitimate.
    """
    data = DataStructure()
    data.tfcoil.i_tf_sup = TFConductorModel.SUPERCONDUCTING
    data.build.dr_tf_inboard = 0.0
    data.numerics.n_iteration_variables = 1
    data.numerics.ixc[0] = 140

    assert "dr_tf_inboard" not in _validation_error_message(data)


def test_zero_thickness_resistive_tf_is_rejected():
    """Resistive TF coils use dr_tf_inboard through the same radial-build
    path as superconducting ones, so a zero thickness is equally invalid.
    """
    data = DataStructure()
    data.tfcoil.i_tf_sup = TFConductorModel.WATER_COOLED_COPPER
    data.build.dr_tf_inboard = 0.0

    with pytest.raises(ProcessValidationError, match="dr_tf_inboard"):
        check_process(data, None)


def test_explicit_thickness_resistive_tf_is_accepted():
    data = DataStructure()
    data.tfcoil.i_tf_sup = TFConductorModel.WATER_COOLED_COPPER
    data.build.dr_tf_inboard = 0.5

    assert "dr_tf_inboard" not in _validation_error_message(data)


@pytest.mark.parametrize("bad_value", [0.0, -0.5])
def test_thickness_iteration_variable_does_not_exempt(bad_value):
    """ixc = 13 does not exempt a non-positive input value: it seeds the
    first model evaluation, and an exactly-zero value is only rejected
    later by the generic iteration-variable check.  A negative value
    cannot come from an input file (the parser bounds dr_tf_inboard to
    [0, 10]), but check_process guards the data structure however it was
    populated.
    """
    data = DataStructure()
    data.tfcoil.i_tf_sup = TFConductorModel.SUPERCONDUCTING
    data.build.dr_tf_inboard = bad_value
    data.numerics.n_iteration_variables = 1
    data.numerics.ixc[0] = 13

    with pytest.raises(ProcessValidationError, match="dr_tf_inboard"):
        check_process(data, None)


def test_stellarator_is_not_checked():
    """Stellarators calculate dr_tf_inboard during the model run."""
    data = DataStructure()
    data.stellarator.istell = 1
    data.tfcoil.i_tf_sup = TFConductorModel.SUPERCONDUCTING
    data.build.dr_tf_inboard = 0.0

    assert "dr_tf_inboard" not in _validation_error_message(data)
