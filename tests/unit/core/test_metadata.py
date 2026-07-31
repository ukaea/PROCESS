from dataclasses import dataclass
from typing import Annotated

import pytest

from process.core.metadata import Parameter, ParameterMetadata, PROCESSModelData


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
    mutable_default: Parameter = 123.0


@dataclass
class ShouldFailBareParameterAnnotated(PROCESSModelData):
    mutable_default: Annotated[Parameter, ParameterMetadata()] = 123.0


@pytest.mark.parametrize(
    "dataclass_class",
    [ShouldFailBareParameter, ShouldFailBareParameterAnnotated],
)
def test_bare_parameter_fails(dataclass_class):
    with pytest.raises(TypeError, match="is typed as a bare"):
        dataclass_class()


@dataclass
class ShouldFailWrongAnnotation(PROCESSModelData):
    mutable_default: Annotated[Parameter, "not_a_ParameterMetadata"] = 123.0


def test_wrong_annotation_fails():
    with pytest.raises(TypeError, match="is annotated with the wrong type of data"):
        ShouldFailWrongAnnotation()
