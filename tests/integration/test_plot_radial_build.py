"""Integration tests for plot_radial_build.py."""
from process.io import plot_radial_build


def test_plot_radial_build(temp_data, mfile_name):
    """Run plot_scans script on a scan MFILE.DAT and check for a PDF output.

    :param temp_data: temporary data dir
    :type temp_data: Path
    :param scan_mfile_name: name of the mfile in the data dir
    :type scan_mfile_name: str
    """
    mfile = temp_data / mfile_name

    plot_radial_build.main(args=["-f", str(mfile), "--outputdir", str(temp_data)])

    assert len(list(temp_data.glob("*.pdf")))
