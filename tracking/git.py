"""Simple submodule to provide access to some git attributes about the repository"""

import subprocess  # noqa: S404
from pathlib import Path


def git_commit_message(directory=str | bytes | None) -> str:
    """Get the commit message for `directory` or the current directory if not provided

    Parameters
    ----------
    :param directory: the directory to get the commit message at

    Returns
    -------
    :
        the commit message for the git repository

    Raises
    ------
    subprocess.CalledProcessError
        if the directory checked is not a git repository.

    """
    if directory is None:
        directory = Path.cwd()

    commit_message = subprocess.run(  # noqa: S602
        "git log -1 --pretty=%B",  # noqa: S607
        shell=True,
        capture_output=True,
        cwd=directory,
        check=True,
    )

    return commit_message.stdout.decode()


def git_commit_hash(directory=None) -> str:
    """Get the commit hash for `directory` or the current directory if not provided

    Parameters
    ----------
    :param directory: the directory to get the commit hash at

    Returns
    -------
    :
        the commit hash for the git repository

    Raises
    ------
    subprocess.CalledProcessError
        if the directory checked is not a git repository.
    """
    if directory is None:
        directory = Path.cwd()

    commit_hash = subprocess.run(  # noqa: S602
        'git log --format="%H" -n 1',  # noqa: S607
        shell=True,
        capture_output=True,
        cwd=directory,
        check=True,
    )

    return commit_hash.stdout.decode()
