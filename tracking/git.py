"""Simple submodule to provide access to some git attributes about the current repository"""

from pathlib import Path
import subprocess


def git_commit_message(directory=None) -> str:
    """Get the commit message for `directory` or, if `directory` is not provided, the current working directory.

    Will raise a subprocess.CalledProcessError if the directory checked is not a git repository.

    :param directory: the directory to get the commit message at
    :type directory: str | bytes | path-like

    :return commit_message: the commit message for the git repository
    :type commit_message: str
    """
    if directory is None:
        directory = Path.cwd()

    commit_message = subprocess.run(
        "git log -1 --pretty=%B",
        shell=True,
        capture_output=True,
        cwd=directory,
        check=True,
    )

    return commit_message.stdout.decode()


def git_commit_hash(directory=None) -> str:
    """Get the commit hash for `directory` or, if `directory` is not provided, the current working directory.

    Will raise a subprocess.CalledProcessError if the directory checked is not a git repository.

    :param directory: the directory to get the commit hash for
    :type directory: str | bytes | path-like

    :return commit_hash: the commit hash for the git repository
    :type commit_hash: str
    """
    if directory is None:
        directory = Path.cwd()

    commit_hash = subprocess.run(
        'git log --format="%H" -n 1',
        shell=True,
        capture_output=True,
        cwd=directory,
        check=True,
    )

    return commit_hash.stdout.decode()
