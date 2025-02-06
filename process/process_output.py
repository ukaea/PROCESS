import numpy as np

from process.fortran import constants, process_output_fortran

# necessary to avoid using process_output in the code through
# two different interfaces
ocentr = process_output_fortran.ocentr
ostars = process_output_fortran.ostars
oheadr = process_output_fortran.oheadr
oshead = process_output_fortran.oshead
oblnkl = process_output_fortran.oblnkl
osubhd = process_output_fortran.osubhd
ocmmnt = process_output_fortran.ocmmnt
write = process_output_fortran.write
dblcol = process_output_fortran.dblcol
ovarin = process_output_fortran.ovarin
ovarst = process_output_fortran.ovarst
obuild = process_output_fortran.obuild


def ovarre(file, descr: str, varnam: str, value, output_flag: str = ""):
    replacement_character = "_"
    if file != constants.mfile:
        replacement_character = " "

    description = f"{descr:<72}".replace(" ", replacement_character)
    varname = f"{varnam:<30}".replace(" ", replacement_character)

    if isinstance(value, np.ndarray):
        value = value.item()
    elif isinstance(value, str):
        value = float(value)

    line = f"{description}{replacement_character} {varname}{replacement_character} {value:.17e} {output_flag}"
    write(file, line)
    if file != constants.mfile:
        ovarre(constants.mfile, descr, varnam, value, output_flag)


def ocosts(file, varnam: str, descr: str, value):
    ovarre(file, descr, varnam, value)


def ovarrf(file, descr: str, varnam: str, value, output_flag: str = ""):
    ovarre(file, descr, varnam, value, output_flag)
