"""Integration tests for plot_summary.py."""

from shutil import copy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from process.core.io.mfile import MFile
from process.core.io.plot.cli import plot_summary_cli
from process.core.io.plot.summary.profiles.bootstrap import plot_sauter_bootstrap_summary
from process.models.geometry.plasma import plasma_geometry


def test_input_file(temp_data, mfile_name, cli_runner):
    """Run plot_summary on an input MFILE and check for an output.

    :param temp_data: temporary data dir
    :type temp_data: Path
    :param mfile_name: name of the mfile in the data dir
    :type mfile_name: str
    """
    mfile = temp_data / mfile_name
    mfile_str = str(mfile)
    mfile_str = str(mfile)
    cli_runner(plot_summary_cli, args=["-f", mfile_str])

    # Assert the pdf has been created, named after the MFILE
    assert (temp_data / f"{mfile_name}.SUMMARY.pdf").exists()


def test_input_file_cwd(temp_data_cwd, mfile_name, cli_runner):
    """Run plot_summary on an MFILE in the cwd.

    :param temp_data_cwd: temporary data dir, which is also the cwd
    :type temp_data_cwd: Path
    :param mfile_name: name of the mfile in the data dir
    :type mfile_name: str
    """
    # Copy the mfile to its default name
    copy(temp_data_cwd / mfile_name, temp_data_cwd / "MFILE.DAT")

    # Run plot_summary with no args, which will look for the default-named mfile
    cli_runner(plot_summary_cli, args=["-f", (temp_data_cwd / "MFILE.DAT").as_posix()])

    # Assert the pdf has been created, named after the MFILE
    assert (temp_data_cwd / "MFILE.DAT.SUMMARY.pdf").exists()


@pytest.mark.parametrize("shape", [0, 1])
def test_sauter_summary_elliptical_surfaces(temp_data, mfile_name, monkeypatch, shape):
    mfile = MFile(temp_data / mfile_name)
    get_value = mfile.get
    monkeypatch.setattr(
        mfile,
        "get",
        lambda variable, *, scan=-1: (
            shape if variable == "i_plasma_shape" else get_value(variable, scan=scan)
        ),
    )
    fig = plt.figure(figsize=(12, 9), dpi=80)
    try:
        plot_sauter_bootstrap_summary(fig, mfile, scan=-1)
        fig.canvas.draw()
        geometry = next(
            axis
            for axis in fig.axes
            if axis.get_title().startswith("Assumed flux surfaces in plasma")
        )
        rmajor = mfile.get("rmajor")
        rminor = mfile.get("rminor")
        area = mfile.get("a_plasma_poloidal")
        kappa = area / (np.pi * rminor**2)
        assert geometry.get_title() == (
            "Assumed flux surfaces in plasma\n"
            rf"Effective elongation $\kappa_{{\rm eff}} = {kappa:.3f}$"
        )
        boundary = plasma_geometry(
            rmajor=rmajor,
            rminor=rminor,
            triang=mfile.get("triang"),
            kappa=mfile.get("kappa"),
            i_single_null=mfile.get("i_single_null"),
            i_plasma_shape=shape,
            square=mfile.get("plasma_square"),
        )
        expected_r = boundary.rs if shape == 0 else [boundary.rs]
        expected_z = boundary.zs if shape == 0 else [boundary.zs]
        boundary_lines = geometry.lines[: len(expected_r)]
        for line, rs, zs in zip(boundary_lines, expected_r, expected_z, strict=True):
            np.testing.assert_allclose(line.get_xdata(), rs)
            np.testing.assert_allclose(line.get_ydata(), zs)
            assert line.get_zorder() < 2
        assert boundary_lines[0].get_label() == "PROCESS plasma outline"
        assert geometry.collections or geometry.patches
        assert all(
            fill.get_zorder() < boundary_lines[0].get_zorder()
            for fill in (*geometry.collections, *geometry.patches)
        )
        theta = np.linspace(0.0, 2.0 * np.pi, 361)
        ellipse_lines = geometry.lines[len(expected_r) : len(expected_r) + 4]
        for line, fraction in zip(ellipse_lines, (0.25, 0.5, 0.75, 1.0), strict=True):
            assert line.get_label() == f"Assumed $r/a={fraction:g}$"
            np.testing.assert_allclose(
                line.get_xdata(), rmajor + fraction * rminor * np.cos(theta)
            )
            np.testing.assert_allclose(
                line.get_ydata(), kappa * fraction * rminor * np.sin(theta)
            )
        outer_r = ellipse_lines[3].get_xdata()
        outer_z = ellipse_lines[3].get_ydata()
        assert (np.max(outer_r) - np.min(outer_r)) / 2 == pytest.approx(rminor)
        assert np.pi * rminor * np.max(outer_z) == pytest.approx(area)
        assert geometry.get_aspect() == pytest.approx(1.0)
        renderer = fig.canvas.get_renderer()
        geometry_bounds = geometry.get_tightbbox(renderer)
        assert all(
            not axis.get_tightbbox(renderer).overlaps(geometry_bounds)
            for axis in fig.axes
            if axis is not geometry
        )
    finally:
        plt.close(fig)
