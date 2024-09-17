"""Integration tests for costs_bar.py."""

from process.io import costs_bar


def test_input_file(temp_data_cwd, mfile_name):
    """Run costs_bar on an input MFILE and check for a pdf output.

    :param temp_data: temporary data dir, which is also the cwd
    :type temp_data: Path
    :param mfile_name: name of the mfile in the data dir
    :type mfile_name: str
    """
    mfile = temp_data_cwd / mfile_name
    mfile_str = str(mfile)
    costs_bar.main(args=["-f", mfile_str, "-s"])

    # Assert a pdf has been created
    assert len(list(temp_data_cwd.glob("*.pdf")))
