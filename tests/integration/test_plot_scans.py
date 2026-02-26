"""Integration tests for plot_scans.py."""

from process.core.io import plot_scans


def test_plot_scans(temp_data, scan_mfile_name):
    """Run plot_scans script on a scan MFILE.DAT and check for a PDF output.

    :param temp_data: temporary data dir
    :type temp_data: Path
    :param scan_mfile_name: name of the scan mfile in the data dir
    :type scan_mfile_name: str
    """
    mfile = temp_data / scan_mfile_name

    plot_scans.main(
        args=[
            "-f",
            str(mfile),
            "-yv",
            "p_plant_electric_net_mw",
            "--outputdir",
            str(temp_data),
        ]
    )

    assert len(list(temp_data.glob("*.pdf")))


def test_plot_scans_stack(temp_data, scan_mfile_name):
    """Run plot_scans script with stacked plots switch on a scan MFILE.DAT and check for a PDF output.

    :param temp_data: temporary data dir
    :type temp_data: Path
    :param scan_mfile_name: name of the scan mfile in the data dir
    :type scan_mfile_name: str
    """
    mfile = temp_data / scan_mfile_name

    plot_scans.main(
        args=[
            "-f",
            str(mfile),
            "-yv",
            "p_plant_electric_net_mw b_plasma_toroidal_on_axis rmajor",
            "-stc",
            "--outputdir",
            str(temp_data),
        ]
    )

    assert len(list(temp_data.glob("*.pdf")))


def test_plot_scans_2d_contour(temp_data, scan_2d_mfile_name):
    """Run plot_scans script with 2D contour plot switch on a scan MFILE.DAT and check for a PDF output.

    :param temp_data: temporary data dir
    :type temp_data: Path
    :param scan_mfile_name: name of the scan mfile in the data dir
    :type scan_mfile_name: str
    """
    mfile = temp_data / scan_2d_mfile_name

    plot_scans.main(
        args=[
            "-f",
            str(mfile),
            "-yv",
            "beta_total_vol_avg",
            "-2DC",
            "--outputdir",
            str(temp_data),
        ]
    )

    assert len(list(temp_data.glob("*.pdf")))
