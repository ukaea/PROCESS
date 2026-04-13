"""Integration tests for costs plotting"""

import pytest

from process.core.io.plot.costs.cli import bar_plot, pie_plot


@pytest.mark.parametrize("cmd", [bar_plot, pie_plot])
def test_input_file(cmd, temp_data_cwd, mfile_name, cli_runner):
    """Run costs plots on an input mfile and check a pdf is produced"""
    mfile = temp_data_cwd / mfile_name
    mfile_str = str(mfile)
    cli_runner(cmd, args=[mfile_str, "-s"])

    # Assert a pdf has been created
    assert list(temp_data_cwd.glob("*.pdf"))
