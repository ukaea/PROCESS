"""Integration tests for plot_sankey.py."""

from process.core.io.plot.cli import sankey


def test_plot_sankey(temp_data_cwd, mfile_name, cli_runner):
    """Assert plot_sankey can make a pdf in the cwd from an mfile.

    :param temp_data_cwd: temp path to data dir, which is also the cwd
    :type temp_data_cwd: Path
    :param mfile_name: name of the mfile to be used
    :type mfile_name: str
    """
    mfile_path = temp_data_cwd / mfile_name
    mfile_path_str = str(mfile_path)
    cli_runner(sankey, args=["-f", mfile_path_str])

    # Assert a pdf has been created
    assert len(list(temp_data_cwd.glob("*.pdf"))) > 0
