from pathlib import Path

from process.io.mfile_utils import get_mfile_initial_ixc_values


def test_get_mfile_initial_ixc_values(input_file):
    iteration_variable_names, iteration_variable_values = get_mfile_initial_ixc_values(
        Path(input_file)
    )

    assert iteration_variable_names[0] == "b_plasma_toroidal_on_axis"
    assert iteration_variable_values[0] == 5.7

    assert iteration_variable_names[1] == "rmajor"
    assert iteration_variable_values[1] == 8.0

    assert iteration_variable_names[-1] == "dr_tf_wp_with_insulation"
    assert iteration_variable_values[-1] == 0.5

    # A default not provided in the MFile
    assert iteration_variable_names[-4] == "f_nd_alpha_electron"
    assert iteration_variable_values[-4] == 0.1
