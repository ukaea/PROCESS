"""Contains method related to the PROCESS repository/installation"""

from pathlib import Path

import process


def get_process_root() -> Path:
    """Returns the root directory of PROCESS.

    E.g. '/home/user/process'
    """

    return Path(process.__file__).resolve().parent
