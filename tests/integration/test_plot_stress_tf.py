"""Integration tests for plot_stress_tf.py."""

from process.core.io.plot.cli import plot_tf_stress


def test_input_file(temp_data_cwd, cli_runner):
    """Run plot_stress_tf on an input MFILE and check for a pdf output.

    :param temp_data: temporary data dir, which is also the cwd
    :type temp_data: Path
    """
    mfile = temp_data_cwd / "SIG_TF.json"
    mfile_str = str(mfile)
    cli_runner(plot_tf_stress, args=["-f", mfile_str])

    # Assert a pdf has been created
    assert len(list(temp_data_cwd.glob("*.pdf")))
