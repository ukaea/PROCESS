"""Integration tests for plot_scans.py."""
from process.io import plot_scans


def test_plot_scans(temp_data, scan_mfile_name):
    """Run plot_scans script on a scan MFILE.DAT and check for a PDF output.

    :param temp_data: temporary data dir
    :type temp_data: Path
    :param scan_mfile_name: name of the scan mfile in the data dir
    :type scan_mfile_name: str
    """
    mfile = temp_data / scan_mfile_name

    plot_scans.main(
        args=["-f", str(mfile), "-yv", "pnetelmw", "--outputdir", str(temp_data)]
    )

    assert len(list(temp_data.glob("*.pdf")))
