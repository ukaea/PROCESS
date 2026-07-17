from dataclasses import fields
from typing import get_origin

from parameter_frame import Parameter


class PROCESSModelData:
    def __new__(cls, *args, **kwargs):
        if not hasattr(cls, "__dataclass_fields__"):
            raise TypeError(f"{cls.__name__} must be a dataclass!")

        return super().__new__(cls, *args, **kwargs)

    def __post_init__(self):
        for f in fields(self):
            # Check for non-generic types that are Parameters
            if isinstance(f.type, type) and issubclass(f.type, Parameter):
                error_msg = (
                    f"{f.name} is typed as a bare {f.type.__name__} on dataclass {self.__class__.__name__}. "
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
