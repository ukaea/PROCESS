import copy
from dataclasses import dataclass
from typing import Annotated

import numpy as np
import pytest

import process.core.data_structure.parameter
from process.core.data_structure.parameter import (
    Parameter,
    ParameterMetadata,
    PROCESSModelData,
)


@pytest.fixture
def turn_off_access_records(monkeypatch):
    monkeypatch.setattr(
        process.core.data_structure.parameter, "KEEP_EDIT_USE_RECORDS", False
    )
    monkeypatch.setattr(
        process.core.data_structure.parameter,
        "FILTER_EDIT_USE_RECORDS_PATH",
        lambda frame: "/tests/" in frame.filename,
    )


@pytest.fixture
def turn_on_access_records(monkeypatch):
    monkeypatch.setattr(
        process.core.data_structure.parameter, "KEEP_EDIT_USE_RECORDS", True
    )
    monkeypatch.setattr(
        process.core.data_structure.parameter,
        "FILTER_EDIT_USE_RECORDS_PATH",
        lambda frame: "/tests/" in frame.filename,
    )


@dataclass
class ShouldFailMutatableDefault(PROCESSModelData):
    mutable_default: Parameter[float] = Parameter("mutable_default", 123.0)  # noqa: RUF009


@dataclass
class ShouldFailMutatableDefaultAnnotated(PROCESSModelData):
    mutable_default: Annotated[Parameter[float], ParameterMetadata()] = Parameter(  # noqa: RUF009
        "mutable_default", 123.0
    )


@pytest.mark.parametrize(
    "dataclass_class",
    [ShouldFailMutatableDefault, ShouldFailMutatableDefaultAnnotated],
)
def test_mutable_default_fails(dataclass_class):
    with pytest.raises(TypeError, match="Initialise the field as a bare constant"):
        dataclass_class()


@dataclass
class ShouldFailBareParameter(PROCESSModelData):
    my_param: Parameter = 123.0


@dataclass
class ShouldFailBareParameterAnnotated(PROCESSModelData):
    my_param: Annotated[Parameter, ParameterMetadata()] = 123.0


@pytest.mark.parametrize(
    "dataclass_class",
    [ShouldFailBareParameter, ShouldFailBareParameterAnnotated],
)
def test_bare_parameter_fails(dataclass_class):
    with pytest.raises(TypeError, match="is typed as a bare"):
        dataclass_class()


@dataclass
class ShouldFailWrongAnnotation(PROCESSModelData):
    my_param: Annotated[Parameter, "not_a_ParameterMetadata"] = 123.0


def test_wrong_annotation_fails():
    with pytest.raises(TypeError, match="is annotated with the wrong type of data"):
        ShouldFailWrongAnnotation()


@dataclass
class ExampleModelDataclass(PROCESSModelData):
    my_param: Annotated[
        Parameter[float],
        ParameterMetadata(description="My parameter", long_name="my_parameter"),
    ] = 42.0


@pytest.fixture
def example_model_dataclass():
    return ExampleModelDataclass()


def test_process_model_data_parameter(example_model_dataclass):
    assert example_model_dataclass.my_param == 42.0  # noqa: RUF069
    assert example_model_dataclass.my_param.value == 42.0  # noqa: RUF069
    assert example_model_dataclass.my_param.description == "My parameter"
    assert example_model_dataclass.my_param.long_name == "my_parameter"


def test_process_model_data_parameter_mutation(example_model_dataclass):
    example_model_dataclass.my_param *= 2
    assert isinstance(example_model_dataclass.my_param, Parameter)
    assert example_model_dataclass.my_param == 84

    example_model_dataclass.my_param = example_model_dataclass.my_param * 2  # noqa: PLR6104
    assert isinstance(example_model_dataclass.my_param, Parameter)
    assert example_model_dataclass.my_param == 168


def test_error_when_records_disabled(turn_off_access_records, example_model_dataclass):
    with pytest.raises(RuntimeError, match="Usage records are disabled"):
        _use = example_model_dataclass.my_param.usage_records

    with pytest.raises(RuntimeError, match="Edit records are disabled"):
        _edit = example_model_dataclass.my_param.edit_records


def test_no_records_when_disabled(turn_off_access_records, example_model_dataclass):
    example_model_dataclass.my_param = 100.0
    _some_param = example_model_dataclass.my_param * 2

    assert example_model_dataclass.my_param._used == []
    assert example_model_dataclass.my_param._edited == []


def test_parameter_use_record(turn_on_access_records, example_model_dataclass):
    assert len(example_model_dataclass.my_param.usage_records) == 0

    _some_param = example_model_dataclass.my_param * 2

    assert len(example_model_dataclass.my_param.usage_records) == 1
    assert example_model_dataclass.my_param.usage_records[0].value == 42.0  # noqa: RUF069

    example_model_dataclass.my_param = 7.0
    _some_other_param = example_model_dataclass.my_param + 4.0

    assert len(example_model_dataclass.my_param.usage_records) == 2
    assert example_model_dataclass.my_param.usage_records[1].value == 7.0  # noqa: RUF069


def test_parameter_use_record_not_created_on_edit(
    turn_on_access_records, example_model_dataclass
):
    assert len(example_model_dataclass.my_param.usage_records) == 0

    example_model_dataclass.my_param = 2

    # Used again to get the usage records
    assert len(example_model_dataclass.my_param.usage_records) == 0


def test_parameter_edit_record(turn_on_access_records, example_model_dataclass):
    assert len(example_model_dataclass.my_param.edit_records) == 0

    example_model_dataclass.my_param = 2.0

    assert example_model_dataclass.my_param == 2.0  # noqa: RUF069
    assert len(example_model_dataclass.my_param.edit_records) == 1
    assert example_model_dataclass.my_param.edit_records[0].value == 42.0  # noqa: RUF069
    assert example_model_dataclass.my_param.edit_records[0].new_value == 2.0  # noqa: RUF069

    example_model_dataclass.my_param = 4.0

    assert len(example_model_dataclass.my_param.edit_records) == 2
    assert example_model_dataclass.my_param.edit_records[1].value == 2.0  # noqa: RUF069
    assert example_model_dataclass.my_param.edit_records[1].new_value == 4.0  # noqa: RUF069


def test_parameter_edit_inplace_record(turn_on_access_records, example_model_dataclass):
    example_model_dataclass.my_param *= 2.0

    assert example_model_dataclass.my_param == 84.0  # noqa: RUF069
    assert len(example_model_dataclass.my_param.edit_records) == 1
    assert example_model_dataclass.my_param.edit_records[0].value == 42.0  # noqa: RUF069
    assert example_model_dataclass.my_param.edit_records[0].new_value == 84.0  # noqa: RUF069


def test_parameter_edit_record_another_parameter(
    turn_on_access_records, example_model_dataclass
):
    example_model_dataclass.my_param = Parameter("another_param", 7.0)

    assert example_model_dataclass.my_param == 7.0  # noqa: RUF069
    assert len(example_model_dataclass.my_param.edit_records) == 1
    assert example_model_dataclass.my_param.edit_records[0].value == 42.0  # noqa: RUF069
    assert example_model_dataclass.my_param.edit_records[0].new_value == 7.0  # noqa: RUF069


def test_parameter_edit_record_no_frame_filter(
    turn_on_access_records, example_model_dataclass, monkeypatch
):
    monkeypatch.setattr(
        process.core.data_structure.parameter,
        "FILTER_EDIT_USE_RECORDS_PATH",
        lambda _: False,
    )

    example_model_dataclass.my_param = 72.0

    assert example_model_dataclass.my_param.edit_records[0].frame_code is None
    assert example_model_dataclass.my_param.edit_records[0].frame_file is None
    assert example_model_dataclass.my_param.edit_records[0].frame_function is None
    assert example_model_dataclass.my_param.edit_records[0].frame_lineno is None


@pytest.mark.parametrize("value", [1.0, Parameter("a_param", 7.0), [1.0, 2.0]])
def test_set_field(example_model_dataclass, value):
    example_model_dataclass.set_field("my_param", value)
    assert example_model_dataclass.my_param is value


def test_arrayify_parameter():
    param = Parameter("param", 42.0)
    array = np.array(param)

    assert isinstance(array, np.ndarray)
    np.testing.assert_array_equal(array, np.array(42.0))


@pytest.mark.parametrize(
    "val1",
    [
        np.float64(1.0),
        np.array(1.0),
        np.array([1.0]),
        np.array([[1.0]]),
        np.int_(1),
        np.array([[1]]),
    ],
)
@pytest.mark.parametrize(
    "val2",
    [
        np.float64(2.0),
        np.array(2.0),
        np.array([2.0]),
        np.array([[2.0]]),
        np.int_(2),
        np.array([[2]]),
    ],
)
def test_ops(val1, val2):
    result_normal = val1 + val2

    result_parameter = Parameter("val1", val1) + Parameter("val2", val2)
    assert result_normal == result_parameter
    # non-in-place ops should return the correct type e.g. float, not a Parameter
    # python types WILL be coerced into a Numpy type
    assert type(result_normal) is type(result_parameter)


@pytest.mark.parametrize(
    "val1",
    [
        np.float64(1.0),
        np.array(1.0),
        np.array([1.0]),
        np.array([[1.0]]),
        np.int_(1),
        np.array([[1]]),
    ],
)
@pytest.mark.parametrize(
    "val2",
    [
        np.float64(2.0),
        np.array(2.0),
        np.array([2.0]),
        np.array([[2.0]]),
        np.int_(2),
        np.array([[2]]),
    ],
)
def test_ops_in_place(val1, val2):

    val1_copy = copy.deepcopy(val1)
    val2_copy = copy.deepcopy(val2)
    try:
        val1_copy += val2_copy
    except Exception:  # noqa: BLE001
        pytest.skip("Not a valid in-place addition type combination.")

    param1 = Parameter("param1", val1)
    param2 = Parameter("param2", val2)

    param1 += param2

    assert val1_copy == param1
    # non-in-place ops should return the correct type e.g. float, not a Parameter
    # python types WILL be coerced into a Numpy type
    assert type(val1_copy) is type(param1.value)


def test_one_inplace_param():
    a = np.array([1.0, 2.0, 3.0])
    b = np.array([1.0, 2.0, 3.0])

    result = np.add(a, b)

    x = Parameter("x", np.zeros_like(a))
    id_x = id(x._value)

    result_param = np.add(a, b, out=x)

    assert id(x._value) == id_x
    assert (x == result).all()
    assert (result_param == result).all()
    assert result_param is x


def test_in_place_many_parameters():
    x = np.array([1.5, 2.7, -3.2])
    fractional = np.zeros_like(x)
    integral = np.zeros_like(x)

    result = np.modf(x, out=(fractional, integral))

    fractional_param = Parameter("fractional_param", np.zeros_like(x))
    integral_param = Parameter("integral_param", np.zeros_like(x))

    id_fractional_param = id(fractional_param._value)
    id_integral_param = id(integral_param._value)

    result_param = np.modf(x, out=(fractional_param, integral_param))

    assert len(result_param) == 2
    assert (result[0] == result_param[0]).all()
    assert (result[1] == result_param[1]).all()
    assert (fractional == fractional_param).all()
    assert (integral == integral_param).all()
    assert id(fractional_param._value) == id_fractional_param
    assert id(integral_param._value) == id_integral_param


def test_in_place_parameter_and_mutable():
    x = np.array([1.5, 2.7, -3.2])
    fractional = np.zeros_like(x)
    integral = Parameter("integral", np.float64(0.0))

    with pytest.raises(
        TypeError, match="Numpy ufunc out keyword has a mixture of Parameter"
    ):
        np.modf(x, out=(fractional, integral))
