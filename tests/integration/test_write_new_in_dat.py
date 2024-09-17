"""Integration tests for write_new_in_dat.py."""

from process.io import write_new_in_dat
from process.io.mfile import MFile
from process.io.in_dat import InDat
from pytest import approx


def test_write_new_in_dat(temp_data, mfile_name):
    """Ensure solution vector from MFILE.DAT is copied to new IN.DAT.

    :param temp_data: temporary data dir
    :type temp_data: Path
    :param mfile_name: name of the mfile in the data dir
    :type mfile_name: str
    """
    mfile_path = temp_data / mfile_name
    in_dat_path = temp_data / "ref_IN.DAT"
    new_in_dat_path = temp_data / "new_IN.DAT"
    # Get final value of te optimisation parameter
    mfile = MFile(mfile_path)
    te_exp = mfile.data["te"].get_scan(-1)

    # Write new IN.DAT then inspect value in new input file
    write_new_in_dat.main(
        args=["-f", str(mfile_path), "-i", str(in_dat_path), "-o", str(new_in_dat_path)]
    )
    in_dat = InDat(str(new_in_dat_path))
    te_obs = in_dat.data["te"].get_value

    # Assert mfile value is now the same as IN.DAT value
    assert te_obs == approx(te_exp)
