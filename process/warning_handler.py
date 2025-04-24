import inspect
from dataclasses import dataclass
from typing import ClassVar

import process.process_output as process_output


class ProcessUserWarning(UserWarning):
    pass


@dataclass
class ProcessWarning:
    msg: str
    location: str


class WarningManager:
    _warnings: ClassVar[list[ProcessWarning]] = []

    def __init__(self):
        raise NotImplementedError(f"{self.__class__.__name__} cannot be instantiated.")

    @classmethod
    def create_warning(cls, msg: str, **diagnostics):
        caller = inspect.stack()[1]
        warning = ProcessWarning(
            f"{msg}"
            + ("\n\t" if diagnostics else "")
            + "\n\t".join([f"\t{d}: {v!r}" for d, v in [*diagnostics.items()]]),
            f"{caller.filename} @ line {caller.lineno}",
        )

        # will not add the warning again if it has the same message and line number
        if warning not in cls._warnings:
            cls._warnings.append(warning)

    @classmethod
    def warnings(cls):
        return cls._warnings

    @classmethod
    def reinitialise(cls):
        cls._warnings = []

    @classmethod
    def show_errors(cls, file_unit: int):
        warning_string = (
            "******************************************** Errors and Warnings *********************************************"
            f"\n{cls.warning_string()}"
        )
        print(warning_string)
        process_output.write(file_unit, warning_string)

    @classmethod
    def warning_string(cls):
        return "\n\n".join([f"({w.location}) {w.msg}" for w in cls._warnings])
