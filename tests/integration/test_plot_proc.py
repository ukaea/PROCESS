"""Integration tests for plot_proc.py."""

from shutil import copy

from process.core.io.plot.cli import plot_proc


def test_input_file(temp_data, mfile_name, cli_runner):
    """Run plot_proc on an input MFILE and check for an output.

    :param temp_data: temporary data dir
    :type temp_data: Path
    :param mfile_name: name of the mfile in the data dir
    :type mfile_name: str
    """
    mfile = temp_data / mfile_name
    mfile_str = str(mfile)
    mfile_str = str(mfile)
    cli_runner(plot_proc, args=["-f", mfile_str])

    # Assert a pdf has been created
    assert len(list(temp_data.glob("*.pdf")))


def test_input_file_cwd(temp_data_cwd, mfile_name, cli_runner):
    """Run plot_proc on an MFILE in the cwd.

    :param temp_data_cwd: temporary data dir, which is also the cwd
    :type temp_data_cwd: Path
    :param mfile_name: name of the mfile in the data dir
    :type mfile_name: str
    """
    # Copy the mfile to its default name
    copy(temp_data_cwd / mfile_name, temp_data_cwd / "MFILE.DAT")

    # Run plot_proc with no args, which will look for the default-named mfile
    cli_runner(plot_proc, args=["-f", (temp_data_cwd / "MFILE.DAT").as_posix()])

    # Assert a pdf has been created
    assert len(list(temp_data_cwd.glob("*.pdf")))
