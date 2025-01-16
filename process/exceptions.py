class ProcessException(Exception):
    """A base Exception to derive other PROCESS exceptions from"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args)
        self._diagnostics = kwargs

    def __str__(self):
        exception_message = super().__str__()
        diagnostics_message = "\n".join([
            f"\t{d}: {repr(v)}" for d, v in self._diagnostics.items()
        ])

        if diagnostics_message:
            return f"{exception_message}\n{diagnostics_message}"

        return exception_message


class ProcessValidationError(ProcessException):
    """Exception raised when validating PROCESS input.

    E.g. initial values, constraint/variable combinations, switch combinations"""


class ProcessValueError(ProcessException, ValueError):
    pass
