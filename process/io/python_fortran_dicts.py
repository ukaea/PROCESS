"""Access Python-Fortran dictionaries for Fortran source information.

This ultimately provides Process Python with the ability to access variable
information in the Process Fortran source code.

Process Python can call process.io.python_fortran_dicts.get_dicts() to load
the dicts from the JSON file created and saved at build-time and use them.
"""
from pkg_resources import resource_filename
import json
import logging

logger = logging.getLogger(__name__)


def get_dicts():
    """Return dicts loaded from the JSON file for use in Python.

    :return: Python-Fortran dicts
    :rtype: dict
    """
    # Get the dicts filename from the package data
    dicts_filename = resource_filename("process", "io/python_fortran_dicts.json")

    # Return loaded dicts
    try:
        with open(dicts_filename, "r") as dicts_file:
            return json.load(dicts_file)
    except FileNotFoundError as error:
        logger.exception("Can't find the dicts JSON file")
        raise error
