from contextlib import suppress
from pathlib import Path

import numpy as np

from process.core import constants
from process.data_structure import global_variables, numerics


class OutputFileManager:
    @classmethod
    def open_files(cls, *, mode="w"):
        cls._outfile = open(  # noqa: SIM115
            Path(global_variables.output_prefix + "OUT.DAT"), mode
        )
        cls._mfile = open(  # noqa: SIM115
            Path(global_variables.output_prefix + "MFILE.DAT"), mode
        )

    @classmethod
    def open_idempotence_files(cls):
        cls._outfile.close()
        cls._mfile.close()

        cls._outfile = open(  # noqa: SIM115
            Path(global_variables.output_prefix + "IDEM_OUT.DAT"), "w"
        )
        cls._mfile = open(  # noqa: SIM115
            Path(global_variables.output_prefix + "IDEM_MFILE.DAT"), "w"
        )

    @classmethod
    def close_idempotence_files(cls):
        Path(cls._outfile.name).unlink()
        Path(cls._mfile.name).unlink()
        cls._outfile.close()
        cls._mfile.close()
        cls.open_files(mode="a")

    @classmethod
    def finish(cls):
        cls._outfile.close()
        cls._mfile.close()


def write(file, string: str):
    if file == constants.MFILE:
        OutputFileManager._mfile.write(f"{string}\n")  # noqa: SLF001
    elif file == constants.NOUT:
        OutputFileManager._outfile.write(f"{string}\n")  # noqa: SLF001
    elif file == constants.IOTTY:
        print(string)


def ocentr(file, string: str, width: int, *, character="*"):
    """Write a centred header within a line of characters to a file

    Parameters
    ----------
    file :
        the integer unit of the file
    string :
        the heading text
    width :
        the desired with of the header
    character :
        the character to pad the heading with (*) (Default value = "*")

    """
    write(file, f"{f' {string} ':{character}^{width}}")
    write(constants.MFILE, f"# {string} #")


def ostars(file, width: int, *, character="*"):
    """Write a line of characters to a file

    Parameters
    ----------
    file :
        the integer unit of the file
    width :
        the desired with of the line
    character :
        the character to fill the line with (*) (Default value = "*")

    """
    write(file, character * width)


def oheadr(file, string: str, *, width: int = 110, character="*"):
    """Write a centred header within a line of characters between two blank lines

    Parameters
    ----------
    file :
        the integer unit of the file
    string :
        the heading text
    width :
        the desired with of the header
    character :
        the character to pad the heading with (*) (Default value = "*")
    """
    oblnkl(file)
    ocentr(file, string, width, character=character)
    oblnkl(file)


def oshead(file, string: str, *, width: int = 80, character="*"):
    """Write a short centred header within a line of characters between two blank lines

    Parameters
    ----------
    file :
        the integer unit of the file
    string :
        the heading text
    width :
        the desired with of the header
    character :
        the character to pad the heading with (*) (Default value = "*")
    """
    oheadr(file, string, width=width, character=character)


def oblnkl(file):
    """Write a blank line to a file

    Parameters
    ----------
    file :
        the integer unit of the file
    """
    write(file, " ")


def osubhd(file, string):
    """Write a subheading between two blank lines

    Parameters
    ----------
    file :
        the integer unit of the file
    string :
        the heading text
    """
    oblnkl(file)
    write(file, string)
    oblnkl(file)


def ocmmnt(file, string: str):
    """Write a comment to a file

    Parameters
    ----------
    file :
        the integer unit of the file
    string :
        the comment text
    """
    write(file, string)


def ovarre(file, descr: str, varnam: str, value, output_flag: str = ""):
    replacement_character = "_"
    if file != constants.MFILE:
        replacement_character = " "

    description = f"{descr:<72}".replace(" ", replacement_character)
    varname = f"{varnam:<30}".replace(" ", replacement_character)

    if isinstance(value, np.ndarray):
        value = value.item()
    if isinstance(value, bytes):
        # TODO: remove when Fortran is gone
        value = value.decode().strip()
    if isinstance(value, str):
        # try and convert the value to a float
        # if it fails, leave as a string
        with suppress(ValueError):
            value = float(value)

    format_value = f"{value:.17e}" if isinstance(value, float) else f"{value: >12}"

    if varnam.strip("()") in numerics.name_xc:
        # MDK add ITV label if it is an iteration variable
        # The ITV flag overwrites the output_flag
        output_flag = "ITV"

    line = f"{description}{replacement_character} {varname}{replacement_character} {format_value} {output_flag}"
    write(file, line)
    if file != constants.MFILE:
        ovarre(constants.MFILE, descr, varnam, value, output_flag)


def ocosts(file, varnam: str, descr: str, value):
    ovarre(file, descr, varnam, value)


def ovarrf(file, descr: str, varnam: str, value, output_flag: str = ""):
    ovarre(file, descr, varnam, value, output_flag)


def ovarin(file, descr: str, varnam: str, value, output_flag: str = ""):
    ovarre(file, descr, varnam, value, output_flag)


def ovarst(file, descr: str, varnam: str, value, output_flag: str = ""):
    ovarre(file, descr, varnam, value, output_flag)


def obuild(file, descr: str, thick: float, total: float, variable_name: str = ""):
    write(file, f"{descr:<50}{thick:.3e}{' ':<10}{total:.3e}  {variable_name}")
