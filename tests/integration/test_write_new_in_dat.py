"""Integration tests for write_new_in_dat.py."""

from click.testing import CliRunner
from pytest import approx

from process.core.io.in_dat import InDat
from process.core.io.in_dat.cli import new_indat
from process.core.io.mfile import MFile


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
    # Get final value of te and f_nd_impurity_electrons(13) optimisation parameters
    mfile = MFile(mfile_path)
    te_exp = mfile.data["temp_plasma_electron_vol_avg_kev"].get_scan(-1)
    fimp13_exp = mfile.data["f_nd_impurity_electrons(13)"].get_scan(-1)

    # Write new IN.DAT then inspect value in new input file
    runner = CliRunner()
    runner.invoke(
        new_indat,
        args=["-f", str(mfile_path), "-i", str(in_dat_path), "-o", str(new_in_dat_path)],
    )
    in_dat = InDat(new_in_dat_path)
    te_obs = in_dat.data["temp_plasma_electron_vol_avg_kev"].get_value
    fimp13_obs = in_dat.data["f_nd_impurity_electrons"].get_value[12]

    # Assert mfile values are now the same as IN.DAT value
    assert te_obs == approx(te_exp)
    assert fimp13_obs == approx(fimp13_exp)
