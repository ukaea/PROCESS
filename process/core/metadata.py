import inspect
from collections.abc import Generator
from copy import deepcopy
from dataclasses import asdict, dataclass, fields
from typing import Annotated, Generic, get_args, get_origin

import numpy as np
from parameter_frame import Parameter as DefaultParameter
from parameter_frame import ParameterValueType

KEEP_EDIT_USE_RECORDS = True


@dataclass(slots=True, kw_only=True, frozen=True)
class UseRecord:
    value: ParameterValueType
    frame_file: str
    frame_lineno: int
    frame_function: str
    frame_code: list[str] | None


@dataclass(slots=True, kw_only=True, frozen=True)
class EditRecord(UseRecord):
    new_value: ParameterValueType


class Parameter(DefaultParameter, Generic[ParameterValueType]):
    def __init__(
        self,
        name: str,
        value: ParameterValueType,
        unit: str = "",
        source: str = "",
        description: str = "",
        long_name: str = "",
        symbol: str = "",
        latex_symbol: str = "",
        _value_types: tuple[type, ...] | None = None,
    ):
        self._latext_symbol = latex_symbol
        self._symbol = symbol

        self._edited = []
        self._used = []
        super().__init__(name, value, unit, source, description, long_name, _value_types)

    @property
    def symbol(self) -> str:
        return self._symbol

    @property
    def latex_symbol(self) -> str:
        return self._latex_symbol

    def __eq__(self, o, /):
        """Check if this parameter is equal to something.

        Parameters are equal if their names and values (with matching
        units) are equal.

        In PROCESS, a parameter is equal to a non-Parameter 'o' if the
        parameter value equals 'o'.

        Returns
        -------
        :
            True if the parameters are equal, False otherwise.
        """
        if not isinstance(o, DefaultParameter):
            return self.value == o
        return super().__eq__(o)

    def __hash__(self):
        return super().__hash__()

    def reset_edit_use_records(self):
        self._edited = []
        self._used = []

    @property
    def edit_records(self) -> list[EditRecord]:
        return deepcopy(self._edited)

    @property
    def usage_records(self) -> list[UseRecord]:
        return deepcopy(self._used)


@dataclass(slots=True, kw_only=True)
class ParameterMetadata:
    unit: str = ""
    source: str = ""
    description: str = ""
    long_name: str = ""


class PROCESSModelData:
    __slots__ = []

    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, "__dataclass_fields__"):
            raise TypeError(f"{cls.__name__} must be a dataclass!")

        return super().__new__(cls, *args, **kwargs)

    def __post_init__(self):
        for f in fields(self):
            current_value = getattr(self, f.name)
            # Check that the Parameter has not been instantiated yet (this will cause
            # issues and is bad with dataclasses)
            if isinstance(current_value, Parameter):
                error_msg = (
                    f"Field {f.name} is initialised as a {type(current_value).__name__}."
                    " This is dangerous as it is mutable!"
                    f" Initialise the field as a bare constant e.g. {f.name}: {f.type!r}"
                    " = 0.0"
                )

                raise TypeError(error_msg)

            field_type = f.type
            # Extract the metadata
            origin_type = get_origin(field_type)
            metadata = {}
            if isinstance(origin_type, type) and issubclass(origin_type, Annotated):
                field_type, metad = get_args(field_type)

                if not isinstance(metad, ParameterMetadata):
                    error_msg = (
                        f"{f.name} is annotated with the wrong type of data"
                        f" ({type(metad).__name__}), expected ParameterMetadata."
                    )
                    raise TypeError(error_msg)

                metadata = asdict(metad)

            # Check for non-generic types that are Parameters
            if isinstance(field_type, type) and issubclass(field_type, Parameter):
                error_msg = (
                    f"{f.name} is typed as a bare {field_type.__name__}"
                    f" on dataclass {self.__class__}."
                    f" You must specify a generic e.g."
                    f" {field_type.__name__}[{type(f.default).__name__}]."
                )
                raise TypeError(error_msg)

            # Get it again in case the field_type has changed (when Annotated)
            origin_type = get_origin(field_type)

            # Make the field a Parameter if it:
            # 1. Is generic (origin_type is None for concrete types)
            # 2. Is a Parameter
            if (origin_type is not None) and (issubclass(origin_type, Parameter)):
                parameter = Parameter(f.name, getattr(self, f.name), **metadata)
                setattr(self, f.name, parameter)

    def __setattr__(self, name, value):
        # we are setting this attribute for the first time (e.g. creating the dataclass)
        if not hasattr(self, name):
            super().__setattr__(name, value)

        current_value = getattr(self, name)
        if KEEP_EDIT_USE_RECORDS and isinstance(current_value, Parameter):
            try:
                called_from = next(
                    filter(lambda frame: "/models/" in frame.filename, inspect.stack())
                )
            except StopIteration:
                pass
            else:
                current_value._edited.append(
                    EditRecord(
                        value=np.copy(current_value.value),
                        new_value=np.copy(value.value)
                        if isinstance(value, Parameter)
                        else np.copy(value),
                        frame_file=called_from.filename,
                        frame_lineno=called_from.lineno,
                        frame_function=called_from.function,
                        frame_code=called_from.code_context,
                    )
                )

        # Not everything is a Parameter in PROCESS
        if isinstance(current_value, Parameter) and not isinstance(value, Parameter):
            current_value.set_value(value)
            return
        super().__setattr__(name, value)

    def __getattribute__(self, name):
        if KEEP_EDIT_USE_RECORDS and (
            isinstance(current_value := super().__getattribute__(name), Parameter)
        ):
            try:
                called_from = next(
                    filter(lambda frame: "/models/" in frame.filename, inspect.stack())
                )
            except StopIteration:
                return current_value

            current_value._used.append(
                UseRecord(
                    value=np.copy(current_value.value),
                    frame_file=called_from.filename,
                    frame_lineno=called_from.lineno,
                    frame_function=called_from.function,
                    frame_code=called_from.code_context,
                )
            )
            return current_value

        return super().__getattribute__(name)

    def parameters(self) -> Generator[tuple[str, Parameter], None, None]:
        return (
            (field.name, param)
            for field in fields(self)
            if isinstance(param := getattr(self, field.name), Parameter)
        )

    def reset_edit_use_records(self):
        for field in fields(self):
            value = getattr(self, field.name)

            if isinstance(value, Parameter):
                value.reset_edit_use_records()
