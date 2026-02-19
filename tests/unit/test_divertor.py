"""Unit tests for divertor.f90 subroutines/functions"""

import pytest

from process.data_structure import divertor_variables as dv
from process.data_structure import tfcoil_variables as tfv
from process.models.divertor import Divertor


@pytest.fixture
def divertor():
    """Provides Divertor object for testing.

    :returns: initialised Divertor object
    :rtype: process.divertor.Divertor
    """
    return Divertor()


class TestDivertor:
    def test_divtart(self, monkeypatch, divertor):
        """Test the divtart subroutine.

        Uses test data from the second call of this subroutine by FNSF regression test.

        :param monkeypatch: pytest mocking fixture
        :type monkeypatch: object

        :param divertor: fixture containing an initialised `Divertor` object
        :type divertor: tests.unit.test_divertor.divertor (functional fixture)
        """

        monkeypatch.setattr(tfv, "drtop", 0)

        rmajor = 1.7
        rminor = 0.97142857142857153
        triang = 0.5
        dr_fw_plasma_gap_inboard = 0.09595
        dz_xpoint_divertor = 0.5
        p_plasma_separatrix_mw = 7.7197999809272062
        i_single_null = 0
        dz_divertor = 0.5
        monkeypatch.setattr(dv, "i_div_heat_load", 1)

        expected_pflux_div_heat_load_mw = 0.087770426974167357

        pflux_div_heat_load_mw = divertor.divtart(
            rmajor,
            rminor,
            triang,
            dr_fw_plasma_gap_inboard,
            dz_xpoint_divertor,
            p_plasma_separatrix_mw,
            False,
            i_single_null,
            dz_divertor,
        )

        assert pflux_div_heat_load_mw == pytest.approx(expected_pflux_div_heat_load_mw)

    def test_divwade(self, monkeypatch, divertor):
        """Test the divwade subroutine.

        Uses test data from the second call of this subroutine by FNSF regression test.

        :param monkeypatch: pytest mocking fixture
        :type monkeypatch: object

        :param divertor: fixture containing an initialised `Divertor` object
        :type divertor: tests.unit.test_divertor.divertor (functional fixture)
        """

        monkeypatch.setattr(tfv, "drtop", 0)

        rmajor = 2.0
        rminor = 1.0
        aspect = 2.0
        b_plasma_toroidal_on_axis = 0.5
        b_plasma_poloidal_average = 0.09595
        p_plasma_separatrix_mw = 1.0e2
        f_div_flux_expansion = 2
        nd_plasma_separatrix_electron = 1.0e19
        deg_div_field_plate = 5.0
        rad_fraction_sol = 8.0e-1
        f_p_div_lower = 1.0

        expected_pflux_div_heat_load_mw = 0.58898578

        pflux_div_heat_load_mw = divertor.divwade(
            rmajor,
            rminor,
            aspect,
            b_plasma_toroidal_on_axis,
            b_plasma_poloidal_average,
            p_plasma_separatrix_mw,
            f_div_flux_expansion,
            nd_plasma_separatrix_electron,
            deg_div_field_plate,
            rad_fraction_sol,
            f_p_div_lower,
            False,
        )

        assert pflux_div_heat_load_mw == pytest.approx(expected_pflux_div_heat_load_mw)


@pytest.mark.parametrize(
    "p_plasma_rad_mw, f_ster_div_single, n_divertors, expected",
    [
        (10.0, 0.5, 2, 10.0),
        (0.0, 1.0, 1, 0.0),
        (5.5, 0.2, 3, 3.3),
        (100.0, 0.0, 5, 0.0),
        (7.0, 0.25, 4, 7.0),
    ],
)
def test_set_incident_radiation_power(
    divertor, p_plasma_rad_mw, f_ster_div_single, n_divertors, expected
):
    result = divertor.incident_radiation_power(
        p_plasma_rad_mw, f_ster_div_single, n_divertors
    )
    assert result == pytest.approx(expected)


@pytest.mark.parametrize(
    "p_plasma_neutron_mw, f_ster_div_single, n_divertors, expected",
    [
        (20.0, 0.5, 2, 20.0),
        (0.0, 1.0, 1, 0.0),
        (8.0, 0.25, 4, 8.0),
        (50.0, 0.0, 3, 0.0),
        (12.5, 0.4, 2, 10.0),
    ],
)
def test_set_incident_neutron_power(
    divertor, p_plasma_neutron_mw, f_ster_div_single, n_divertors, expected
):
    result = divertor.incident_neutron_power(
        p_plasma_neutron_mw, f_ster_div_single, n_divertors
    )
    assert result == pytest.approx(expected)
