import pytest

from process.models.physics.exhaust import PlasmaExhaust


def test_calculate_brunner_divertor_power_splits_at_zero_separatrix_separation():
    result = PlasmaExhaust().calculate_brunner_divertor_power_splits(
        dr_outboard_midplane_sep=0.0,
        len_plasma_sol_outboard_power_decay=0.01,
        len_plasma_sol_inboard_power_decay=0.01,
    )

    assert result.f_p_div_inboard_separatrix == pytest.approx(0.16)
    assert result.f_p_div_outboard_separatrix == pytest.approx(0.84)
    assert result.f_p_div_inboard_lower_separatrix == pytest.approx(0.08)
    assert result.f_p_div_inboard_upper_separatrix == pytest.approx(0.08)
    assert result.f_p_div_outboard_lower_separatrix == pytest.approx(0.42)
    assert result.f_p_div_outboard_upper_separatrix == pytest.approx(0.42)


def test_calculate_brunner_divertor_power_splits_asymmetric_case():
    result = PlasmaExhaust().calculate_brunner_divertor_power_splits(
        dr_outboard_midplane_sep=0.01,
        len_plasma_sol_outboard_power_decay=0.005,
        len_plasma_sol_inboard_power_decay=0.0075,
    )

    assert result.f_p_div_inboard_separatrix == pytest.approx(0.4010068950189542)
    assert result.f_p_div_outboard_separatrix == pytest.approx(0.5989931049810458)
    assert result.f_p_div_inboard_lower_separatrix == pytest.approx(0.08365345781749392)
    assert result.f_p_div_inboard_upper_separatrix == pytest.approx(0.31735343720146025)
    assert result.f_p_div_outboard_lower_separatrix == pytest.approx(0.07140172838484167)
    assert result.f_p_div_outboard_upper_separatrix == pytest.approx(0.527591376596204)

    assert (
        result.f_p_div_inboard_lower_separatrix + result.f_p_div_inboard_upper_separatrix
    ) == pytest.approx(result.f_p_div_inboard_separatrix)
    assert (
        result.f_p_div_outboard_lower_separatrix
        + result.f_p_div_outboard_upper_separatrix
    ) == pytest.approx(result.f_p_div_outboard_separatrix)
    assert (
        result.f_p_div_inboard_separatrix + result.f_p_div_outboard_separatrix
    ) == pytest.approx(1.0)
