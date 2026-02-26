"""Integration tests for mfile_to_csv.py."""

from process.core.io import mfile_to_csv


def test_mfile_to_csv(temp_data, mfile_name):
    """Run mfile_to_csv via CLI on an MFILE and check for a CSV output.

    varlist.txt defines which variables to extract to CSV.
    :param temp_data: temporary data dir
    :type temp_data: Path
    :param mfile_name: name of the mfile in the data dir
    :type mfile_name: str
    """
    mfile = temp_data / mfile_name
    varlist = temp_data / "mfile_to_csv_varlist.json"
    mfile_to_csv.main(args=["-f", str(mfile), "-v", str(varlist)])

    # Assert a .csv has been produced
    assert len(list(temp_data.glob("*.csv")))
