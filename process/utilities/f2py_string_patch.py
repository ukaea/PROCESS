import numpy as np
import re
import warnings


def string_to_f2py_compatible(
    target: np.ndarray, string: str = None, except_length: bool = False
) -> str:
    """Return a string that is compatible with f2py interface.
    Fortran character(len=x) will throw an error (by f2py) if the passed string is
    not exactly of length x. This means truncating and lengthening of the strings
    must occur before the string is passed to the Fortran.

    :param target: the Fortran character variable/parameter to be targeted.
    :type target: np.ndarray
    :param string: the string to be passed to the Fortran interface.
    :type string: str
    :default string: None (which becomes '')
    :param except_length: if the string is too long an exception can be raised or the string truncated with a warning.
                          True: raise runtime exception on string too long
                          False: truncate string and raise a runtime warning
    :type except_length: bool
    :default except_length: False

    :returns new_string: string with artificially added whitespace to make string correct length
    :rtype new_string: str
    """
    target_info = re.findall(r"^\|([A-Z]+)([0-9]+)$", str(target.dtype))

    target_type, target_size = target_info[0]

    target_size = int(target_size)

    if string is None:
        string = ""

    if target_type != "S":
        raise TypeError(
            f'{target} is not of type string ("|S<>") is instead {str(target.dtype)}'
        )

    if string:
        if len(string) > target_size and except_length:
            raise RuntimeError(
                f"String string of length {len(string)} is trying to initiate as {target} with length {target_size}"
            )
        elif len(string) > target_size:
            warnings.warn(
                f"String string of length {len(string)} is trying to initiate as {target} with length \
                {target_size}. String string will be truncated!",
            )

        string = string[0:target_size]

        difference = target_size - len(string)
        new_string = string + " " * difference
    else:
        new_string = " " * target_size

    return new_string


def f2py_compatible_to_string(result: np.ndarray) -> str:
    """Strings from the Fortran interface are numpy.ndarray's
    with the string held as a memory buffer. This buffer
    needs to be converted to Python bytes which can then be decoded.
    Finally, the artificial whitespace can be stripped.

    :param result: a character(len=x) variable from the Fortran interface
    :type result: np.ndarray

    :returns string: string with artificially added whitespace removed
    :rtype new_string: str
    """
    return result.tobytes().decode().strip()
