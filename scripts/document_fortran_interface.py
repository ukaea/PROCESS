"""
Documentation tools cannot document the fortan shared object file
due. However, f2py auto-documents the entire interface.

This file processes that interface into a dummy-Python format,
_fortran.py with all the same names and docstrings, however
this is a format that can be documented.

fortran (python module)
├─ fortran module (python class)
│  ├─ subroutine (python class method)
│  ├─ function (python class method)
│  ├─ module variable (python class attribute)
"""

import inspect
from pathlib import Path
from typing import Any, NamedTuple

from process import fortran


class FortranModuleMember(NamedTuple):
    """A container for data about a specific subroutine/function/module contained within a wrapped module.

    :name: the name of the function/subroutine/module variable this object holds data for
    :type name: str

    :docstring: the docstring (autogenerated by f2py) for this function/subroutine/module variable
    :type docstring: str
    """

    name: str
    docstring: str


class FortranModuleVariable(NamedTuple):
    """A container for data about a specific variable contained within a wrapped module.

    :name: the name of the function/subroutine/module variable this object holds data for
    :type name: str

    :docstring: the docstring (autogenerated by f2py) for this function/subroutine/module variable
    :type docstring: str

    :var_type: the type of this variable
    :type var_type: Any

    :value: the default value, if any, of this member
    :type value: Any
    """

    name: str
    docstring: str
    value: Any = None
    var_type: Any = None
    module: Any = None


class FortranModule(NamedTuple):
    """A container for data about a specific wrapped fortran module

    :name: the name of the module this object holds data for
    :type name: str

    :docstring: the docstring (autogenerated by f2py) for this module
    :type docstring: str

    :members: a list of members, ie subroutines/functions/module variables
    :type members: List[FortranModuleMember]
    """

    name: str
    docstring: str
    members: list[FortranModuleMember]


def get_modules(ftrn) -> list[FortranModule]:
    """Returns a list of modules in the wrapped Python-Fortran interface

    :param ftrn: top-level wrapped fortran interface (the name of the module output from f2py) e.g. f2py wrapping PROCESS uses `-m fortran` so this top-level interface is process.fortran
    :type ftrn: f2py-generated Python-Fortran interface

    :return classes: the list of modules
    :type classes: List[FortranModule]
    """
    classes = []

    for name, module in inspect.getmembers(ftrn):
        if type(module) == type(fortran.physics_variables):  # noqa: E721
            classes.append(
                FortranModule(
                    name=name, docstring=module.__doc__, members=get_members(module)
                )
            )

    return classes


def get_members(
    fortran_module,
) -> list[FortranModuleMember | FortranModuleVariable]:
    """Returns a list of members (subroutine, function, module variables) of the module

    :param fortran_module: the Fortran module to get the members of
    :type fortran_module: a wrapped Fortran module

    :return members: a list of subroutines, functions, and variables of the fortran_module
    :type members: List[Union[FortranModuleMember, FortranModuleVariable]]
    """
    members = []

    for name, member in inspect.getmembers(fortran_module):
        if name[0:2] == "__":
            continue

        if type(member) == type(fortran.constants.init_constants):  # noqa: E721
            docstring = member.__doc__
            if is_variable(member):
                members.append(
                    FortranModuleVariable(
                        name=name,
                        docstring=docstring,
                        value=member,
                        var_type=type(member).__qualname__,
                        module=type(member).__module__,
                    )
                )
            else:
                members.append(FortranModuleMember(name=name, docstring=docstring))

    return members


def is_variable(member) -> bool:
    """Checks if the member is a variable or not."""
    # should likely be kept as a function since f2py's interface keeps changing
    return member.__doc__ is None


def create_module_signature(mod: FortranModule) -> str:
    """Creates the signature for a wrapped fortran module, corresponding to one class.
    Manages the generation of import statements, variable signatures, and function signatures.
    """

    docstring = f"Abstract representation of the F2Py-generated wrapper around the {mod.name} module"
    # f2py gives modules unhelpful docstrings, so a more accurate and concise docstring is
    # created for use in documentation

    functions: list[str] = []  # subroutines/functions
    variables: list[str] = []  # module variables
    imports: set[str] = set()  # import non-builtin types
    # set to avoid import duplication (at a class level)
    # duplication could still occur between classes, although
    # neither matter

    for i in mod.members:
        if isinstance(i, FortranModuleVariable):
            variables.append(create_variable_signature(i))
            imports.add(create_import(i))
        else:
            functions.append(create_function_signature(i))

    header = "\n".join([i for i in imports if i]) + "\n\n"
    body = "\n".join(variables) + "\n\n" + "\n\n".join(functions)

    return f'{header}class {mod.name}:\n\t"""{docstring}"""\n{body}'


def create_import(var: FortranModuleVariable) -> str:
    """Creates the import statement for a var's type if not a builtin."""
    # ignore builtins
    if var.module == "builtins":
        return None

    return f"from {var.module} import {var.var_type}"


def create_variable_signature(var: FortranModuleVariable) -> str:
    """Creates the signature (abstract code) for module variables.
    In its abstract representation, it is a class variable.
    """
    base_string = [var.name]
    # if a type is declared, type hint as such
    if var.var_type is not None:
        base_string.append(f": {var.var_type}")
    # if a default is given, show as such
    if var.value is not None:
        base_string.append(f" = {var.value}")

    return "\t" + "".join(base_string)


def create_function_signature(func: FortranModuleMember) -> str:
    """Creates the signature (abstract code) for module functions/subroutines.
    In its abstract representation, it is a class method.
    """
    # sometimes f2py puts odd characters in the docstring
    # this assert statement will show that this has happened
    # and explain a crash of the script
    assert all(ord(c) != 0 for c in func.docstring)

    docstring = func.docstring.replace("\n", "\n\t\t\t").strip("\n\t")
    return f'\t@classmethod\n\tdef {func.name}(cls, *args, **kwargs):\n\t\t"""{docstring}"""\n\t\tpass'


if __name__ == "__main__":
    fortran_module_definitions = [
        create_module_signature(i) for i in get_modules(fortran)
    ]

    string = (
        '"""Abstract definitions of all wrapped modules, including automatically generated docstrings for subroutines/functions created by f2py"""\n\n'
        + "".join(fortran_module_definitions)
    )

    # write _fortran.py to the process package
    # so it will be autodocumented as if it
    # were a regular source file
    current_dir = Path(__file__).resolve().parent
    target_dir = current_dir / "../process"
    with open(target_dir / "_fortran.py", "w") as file:
        file.write(string)
