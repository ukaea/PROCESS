from dataclasses import Field, field, fields
from typing import Any, get_origin

from parameter_frame import Parameter


def parameter_field_factory(
    *parameter_args,
    parameter_class=Parameter,
    field_kwargs: dict[str, Any] | None = None,
    **parameter_kwargs,
) -> Field:
    return field(
        default_factory=lambda: parameter_class(*parameter_args, **parameter_kwargs),
        **(field_kwargs or {}),
    )


class PROCESSModelData:
    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, "__dataclass_fields__"):
            raise TypeError(f"{cls.__name__} must be a dataclass!")

        return super().__new__(cls, *args, **kwargs)

    def __post_init__(self):
        for f in fields(self):
            current_value = getattr(self, f.name)
            if isinstance(current_value, Parameter) and f.name != current_value.name:
                error_msg = (
                    f"Field '{f.name}' of {self.__class__} is a {type(current_value).__name__} "
                    f"with a different name ('{current_value.name}'). "
                    f'\nChange the instantiation to read: {type(current_value).__name__}("{f.name}", ...)'
                )
                raise ValueError(error_msg)

            # Check for non-generic types that are Parameters
            if isinstance(f.type, type) and issubclass(f.type, Parameter):
                error_msg = (
                    f"{f.name} is typed as a bare {f.type.__name__} on dataclass {self.__class__}. "
                    f"You must specify a generic e.g. {f.type.__name__}[float]."
                )
                raise TypeError(error_msg)

            origin_type = get_origin(f.type)

            # Make the field a Parameter if it:
            # 1. Is generic (origin_type is None for concrete types)
            # 2. Is a Parameter
            if (origin_type is not None) and (issubclass(origin_type, Parameter)):
                setattr(self, f.name, Parameter(f.name, getattr(self, f.name)))

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
