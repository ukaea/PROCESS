"""Classes related to handling and collecting metadata on the PROCESS data structure."""

import inspect
import logging
from collections.abc import Callable, Generator
from copy import deepcopy
from dataclasses import asdict, dataclass, fields
from typing import Annotated, Any, Generic, get_args, get_origin

import numpy as np
from parameter_frame import Parameter as DefaultParameter
from parameter_frame import ParameterValueType

logger = logging.getLogger(__name__)

KEEP_EDIT_USE_RECORDS = False
FILTER_EDIT_USE_RECORDS_PATH: Callable[[inspect.FrameInfo], bool] = lambda frame: (  # noqa: E731
    "/models/" in frame.filename
)
"""An optional filter to find the frame to display in the edit/use record.

If FILTER_EDIT_USE_RECORDS_PATH(frame) = True the first frame is returned
(probably not a useful frame because it will be one of the functions in this file)

If FILTER_EDIT_USE_RECORDS_PATH(frame) = False no frame information is recorded.

By default, the functions selects the first frame to come from a PROCESS model.
I.e. process/models/**/*.py
"""


@dataclass(slots=True, kw_only=True, frozen=True)
class UseRecord:
    """A dataclass which records the location where a Parameter is used.

    Notes
    -----
    Use records are only created when the variable is accessed in a file
    in the `model` subdirectory.
    """

    value: ParameterValueType
    """The current value of the Parameter when it is accessed."""
    frame_file: str
    """The file path of the file where the parameter was accessed."""
    frame_lineno: int
    """The line number of `frame_file` where the parameter was accessed."""
    frame_function: str
    """The name of the function in `frame_file` where the parameter was accessed."""
    frame_code: list[str] | None
    """A copy of the Python code (usually a single line) which accesses the parameter."""


@dataclass(slots=True, kw_only=True, frozen=True)
class EditRecord(UseRecord):
    """A dataclass which records the location where a Parameter is edited."""

    new_value: Any
    """The value that the Parameter is being updated to.

    Notes
    -----
    If the new value is also a Parameter only its `.value` is copied over.
    """


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
    def value(self):
        if KEEP_EDIT_USE_RECORDS:
            try:
                called_from = next(filter(FILTER_EDIT_USE_RECORDS_PATH, inspect.stack()))
            except StopIteration:
                called_from = inspect.FrameInfo(None, None, None, None, None, None)

            self._used.append(
                UseRecord(
                    value=np.copy(self._value),
                    frame_file=called_from.filename,
                    frame_lineno=called_from.lineno,
                    frame_function=called_from.function,
                    frame_code=called_from.code_context,
                )
            )

        return super().value

    def set_value(self, new_value, source=""):
        if KEEP_EDIT_USE_RECORDS:
            try:
                called_from = next(
                    filter(
                        FILTER_EDIT_USE_RECORDS_PATH,
                        inspect.stack(),
                    )
                )
            except StopIteration:
                called_from = inspect.FrameInfo(None, None, None, None, None, None)

            self._edited.append(
                EditRecord(
                    value=np.copy(self.history()[0].value),
                    new_value=np.copy(new_value._value)
                    if isinstance(new_value, Parameter)
                    else np.copy(new_value),
                    frame_file=called_from.filename,
                    frame_lineno=called_from.lineno,
                    frame_function=called_from.function,
                    frame_code=called_from.code_context,
                )
            )
        return super().set_value(new_value, source)

    @property
    def edit_records(self) -> list[EditRecord]:
        if not KEEP_EDIT_USE_RECORDS:
            raise RuntimeError(
                f"Edit records are disabled because {KEEP_EDIT_USE_RECORDS = }"
            )
        return deepcopy(self._edited)

    @property
    def usage_records(self) -> list[UseRecord]:
        """The usage records for this Parameter.

        Raises
        ------
        RuntimeError
            KEEP_EDIT_USE_RECORDS is false meaning no uses of this Parameter
            would have been recorded.

        Notes
        -----
        If `dataclass.my_param.usage_records` is called in a file that matches the
        FILTER_EDIT_USE_RECORDS_PATH then this action will create a new use record
        which will be included in the return from this method.
        """
        if not KEEP_EDIT_USE_RECORDS:
            raise RuntimeError(
                f"Usage records are disabled because {KEEP_EDIT_USE_RECORDS = }"
            )
        return deepcopy(self._used)

    def __deepcopy__(self, memo):
        """Create a copy of this Parameter.

        Only copies across the value, not the history or use/edit records etc.
        """
        return self.__class__(name=self._name, value=deepcopy(self._value))


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

        # Do not want a use record to be created here because we editing it
        current_value = getattr(self, name)

        # Not everything is a Parameter in PROCESS
        if isinstance(current_value, Parameter):
            logger.debug(
                f"Doing self.{name} = {value!r} only copies {value} into "
                f"{name}.value. Use set_field to exactly set the dataclass field."
            )
            if isinstance(value, Parameter):
                current_value.set_value(value.value, source=value.name)
                return
            current_value.set_value(value)
            return

        super().__setattr__(name, value)

    def set_field(self, name, value):
        """Forcibly set self.name to value.

        This bypasses the Parameter logic when doing self.name = value which maintains
        the Parameterness of self.name.
        """
        super().__setattr__(name, value)

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


def unwrap_parameter(func):
    def wrapper(*args, **kwargs):
        return func(
            *[arg.value if isinstance(arg, Parameter) else arg for arg in args],
            **{
                k: (v.value if isinstance(v, Parameter) else v)
                for k, v in kwargs.items()
            },
        )

    return wrapper
