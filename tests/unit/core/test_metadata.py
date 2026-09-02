from dataclasses import dataclass
from typing import Annotated

import numpy as np
import pytest

import process.core.metadata
from process.core.metadata import Parameter, ParameterMetadata, PROCESSModelData


@pytest.fixture
def turn_off_access_records(monkeypatch):
    monkeypatch.setattr(process.core.metadata, "KEEP_EDIT_USE_RECORDS", False)
    monkeypatch.setattr(process.core.metadata, "FILTER_EDIT_USE_RECORDS_PATH", "/tests/")


@pytest.fixture
def turn_on_access_records(monkeypatch):
    monkeypatch.setattr(process.core.metadata, "KEEP_EDIT_USE_RECORDS", True)
    monkeypatch.setattr(process.core.metadata, "FILTER_EDIT_USE_RECORDS_PATH", "/tests/")


@dataclass
class ShouldFailMutatableDefault(PROCESSModelData):
    mutable_default: Parameter[float] = Parameter("mutable_default", 123.0)  # ruff:ignore[RUF009]


@dataclass
class ShouldFailMutatableDefaultAnnotated(PROCESSModelData):
    mutable_default: Annotated[Parameter[float], ParameterMetadata()] = Parameter(  # ruff:ignore[RUF009]
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
    assert example_model_dataclass.my_param == 42.0  # ruff:ignore[RUF069]
    assert example_model_dataclass.my_param.value == 42.0  # ruff:ignore[RUF069]
    assert example_model_dataclass.my_param.description == "My parameter"
    assert example_model_dataclass.my_param.long_name == "my_parameter"


def test_process_model_data_parameter_mutation(example_model_dataclass):
    example_model_dataclass.my_param *= 2
    assert isinstance(example_model_dataclass.my_param, Parameter)
    assert example_model_dataclass.my_param == 84

    example_model_dataclass.my_param = example_model_dataclass.my_param * 2  # ruff:ignore[PLR6104]
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
    # Used when the usage records are accessed hence its not 0
    assert len(example_model_dataclass.my_param.usage_records) == 1
    assert example_model_dataclass.my_param.usage_records[0].value == 42.0  # ruff:ignore[RUF069]

    _some_param = example_model_dataclass.my_param * 2

    # Used once in the calculation and thrice to access the usage records
    assert len(example_model_dataclass.my_param.usage_records) == 4
    for record in example_model_dataclass.my_param.usage_records:
        assert record.value == 42.0  # ruff:ignore[RUF069]


def test_parameter_use_record_not_created_on_edit(
    turn_on_access_records, example_model_dataclass
):
    # Only used when accessing the usage records, hence 1 not 0
    assert len(example_model_dataclass.my_param.usage_records) == 1

    example_model_dataclass.my_param = 2

    # Used again to get the usage records
    assert len(example_model_dataclass.my_param.usage_records) == 2


def test_parameter_edit_record(turn_on_access_records, example_model_dataclass):
    assert len(example_model_dataclass.my_param.edit_records) == 0

    example_model_dataclass.my_param = 2.0

    assert example_model_dataclass.my_param == 2.0  # ruff:ignore[RUF069]
    assert len(example_model_dataclass.my_param.edit_records) == 1
    assert example_model_dataclass.my_param.edit_records[0].value == 42.0  # ruff:ignore[RUF069]
    assert example_model_dataclass.my_param.edit_records[0].new_value == 2.0  # ruff:ignore[RUF069]


def test_parameter_edit_inplace_record(turn_on_access_records, example_model_dataclass):
    example_model_dataclass.my_param *= 2.0

    assert example_model_dataclass.my_param == 84.0  # ruff:ignore[RUF069]
    assert len(example_model_dataclass.my_param.edit_records) == 1
    assert example_model_dataclass.my_param.edit_records[0].value == 42.0  # ruff:ignore[RUF069]
    assert example_model_dataclass.my_param.edit_records[0].new_value == 84.0  # ruff:ignore[RUF069]


def test_parameter_edit_record_another_parameter(
    turn_on_access_records, example_model_dataclass
):
    example_model_dataclass.my_param = Parameter("another_param", 7.0)

    assert example_model_dataclass.my_param == 7.0  # ruff:ignore[RUF069]
    assert len(example_model_dataclass.my_param.edit_records) == 1
    assert example_model_dataclass.my_param.edit_records[0].value == 42.0  # ruff:ignore[RUF069]
    assert example_model_dataclass.my_param.edit_records[0].new_value == 7.0  # ruff:ignore[RUF069]


@pytest.mark.parametrize("value", [1.0, Parameter("a_param", 7.0), [1.0, 2.0]])
def test_set_field(example_model_dataclass, value):
    example_model_dataclass.set_field("my_param", value)
    assert example_model_dataclass.my_param is value


def test_arrayify_parameter():
    param = Parameter("param", 42.0)
    array = np.array(param)

    assert isinstance(array, np.ndarray)
    np.testing.assert_array_equal(array, np.array(42.0))
