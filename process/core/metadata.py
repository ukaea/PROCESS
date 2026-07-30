from collections.abc import Generator
from dataclasses import asdict, dataclass, fields
from typing import Annotated, get_args, get_origin

from parameter_frame import Parameter


@dataclass(slots=True, kw_only=True)
class ParameterMetadata:
    unit: str = ""
    source: str = ""
    description: str = ""
    long_name: str = ""


class PROCESSModelData:
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
        if isinstance(current_value, Parameter) and not isinstance(
            value, type(current_value)
        ):
            current_value.set_value(value)
            return
        super().__setattr__(name, value)

    def parameters(self) -> Generator[tuple[str, Parameter], None, None]:
        return (
            (field.name, param)
            for field in fields(self)
            if isinstance(param := getattr(self, field.name), Parameter)
        )
