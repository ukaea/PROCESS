"""Test the plot_solutions tool."""

import pytest
from process.io.plot_solutions import RunMetadata, plot_mfile_solutions
from typing import List, Sequence


@pytest.fixture
def run_metadata(temp_data) -> List[RunMetadata]:
    """Return runs metadata.

    :param temp_data: temporary dir containing data files
    :type temp_data: Path
    :return: list of solutions and their tags
    :rtype: List[RunMetadata]
    """
    run_metadata = [
        RunMetadata(temp_data / "large_tokamak_1_MFILE.DAT", "large tokamak 1"),
        RunMetadata(temp_data / "large_tokamak_2_MFILE.DAT", "large tokamak 2"),
        RunMetadata(temp_data / "large_tokamak_3_MFILE.DAT", "large tokamak 3"),
    ]

    return run_metadata


def test_plot_mfile_solutions(run_metadata: Sequence[RunMetadata]):
    """Plot three mfile solutions and their differences using RunMetadata objects.

    :param run_metadata: run_metadata for plotting
    :type run_metadata: Sequence[RunMetadata]
    """

    _, results_df = plot_mfile_solutions(
        runs_metadata=run_metadata, plot_title="3 large tokamak solutions"
    )

    # Check shape of solution df as evidence of plotting without error
    assert results_df.shape == (3, 93)
