"""
As the Python conversion continues, we want to maintain backwards compatibility with our utility tools.
These tools rely on the so-called 'dictionaries' to provide some of the meta data used in their function.

This extension to the existing `create_dicts.py` provides a means to allow data structures in Python to
still be included in the 'dictionaries'.
"""

import importlib.util
import inspect
import pathlib
import sys
from typing import Any, List, NamedTuple

# the directory which this script resides (scripts/)
CURRENT_DIR = pathlib.Path(__file__).resolve().parent

# the following code allows us to import the Process directory as a package
# ignoring current installs of Process and meaning we get access to the latest
# code changes
# This means the installation of process does not need to be a build dependency
# of this script and a potential race condition can be removed

# find import information for the process directory sitting above this scripts directory
process_spec = importlib.util.spec_from_file_location(
    "process", CURRENT_DIR.parent / "process/__init__.py"
)
# get the module from the above import information
process_module = importlib.util.module_from_spec(process_spec)
# add the process_module (from the directory NOT the possibly pre-installed version)
# to the modules list overwritting any other Process install for this "session"
# of the interpreter
sys.modules[process_spec.name] = process_module
# execute the module so we can use it when we come to import it later
process_spec.loader.exec_module(process_module)

# now that "process" is temporarily available to import
# as a package (despite being a directory) we can install what we need

from process.main import Models  # noqa: E402
from process.variables import AnnotatedVariable  # noqa: E402


class AnnotatedVariableData(NamedTuple):
    """Holds data about a variable being processed before insertion into the dictionaries.

    :param parent: the name of the class this variable was declared within
    :type parent: str

    :param name: the name of the variable
    :type name: str

    :param obj: a pointer to the actual variable, therefore holding the initial value of this annotated variable.
    :type obj: _Variable (hidden within AnnotatedVariable)

    :param docstring: provide a docstring, description, to the variable.
    :type docstring: str

    :param units: the units of this variable
    :type units: str
    """

    parent: str
    name: str
    obj: Any
    docstring: str
    units: str


# _Variable is a hidden class so cannot be used
# to check type. This hackery is needed
# to provde the type checking.
_Variable_module_name = AnnotatedVariable(object).__class__.__module__
_Variable_name = AnnotatedVariable(object).__class__.__name__


def get_non_dunder_class_members(_object: object) -> dict:
    """Given a class, this function returns all members that are not dunder methods/variables.
    It should be noted that member is class variables and class methods.

    :param _object: the object (instantiated class) to be interogated
    :type _object: object

    :return non_dunder_members: a dictionary mapping a members name to the underlying member
    :type non_dunder_members: dict
    """
    non_dunder_members = {}
    for name, value in inspect.getmembers(_object):
        # exclude methods or variables starting double-underscore
        if name[0:2] != "__":
            non_dunder_members[name] = value

    return non_dunder_members


def get_annotated_variables(parent_name: str, _object) -> List[AnnotatedVariableData]:
    """Given a physics and engineering module, this function extracts and returns all of the
    annotated variables as a named tuple with important information contained.

    :param parent_name: the name of the parent (physics and engineering) module
    :type parent_name: str

    :param _object: an instance physics and engineering module to gather annotated variables from
    :type _object: PROCESS physics and engineering module

    :return annotated_variables: array of AnnotatedVariableData
    :type annotated_variables: List[AnnotatedVariableData]
    """
    annotated_variables = []
    for name, member in inspect.getmembers(_object):
        # hackery beget hackery
        # since the underlying _Variable class is hidden,
        # this is the only way to check the underlying class
        if (
            member.__class__.__module__ == _Variable_module_name
            and member.__class__.__name__ == _Variable_name
        ):
            annotated_variables.append(
                AnnotatedVariableData(
                    parent=parent_name,
                    name=name,
                    obj=member,
                    docstring=member.__doc__,
                    units=member.__units__,
                )
            )

    return annotated_variables


def get_python_variables() -> List[AnnotatedVariableData]:
    """Drives the discovery of AnnotatedVariables, returning them back to be
    parsed into the dicts by `create_dicts.py`.
    """
    models_obj = Models()
    # non-dunder methods of Models are all the models
    # that we want to get the variables from
    models_dict = get_non_dunder_class_members(models_obj)

    variables = []
    for name, model in models_dict.items():
        for annotated_variable in get_annotated_variables(name, model):
            variables.append(annotated_variable)

    return variables
