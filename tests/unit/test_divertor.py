"""Unit tests for divertor.f90 subroutines/functions"""

import pytest

from process.divertor import Divertor
from process.fortran import divertor_variables as dv
from process.fortran import tfcoil_variables as tfv


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
        bt = 0.5
        bp = 0.09595
        p_plasma_separatrix_mw = 1.0e2
        f_div_flux_expansion = 2
        nesep = 1.0e19
        deg_div_field_plate = 5.0
        rad_fraction_sol = 8.0e-1
        f_p_div_lower = 1.0

        expected_pflux_div_heat_load_mw = 0.58898578

        pflux_div_heat_load_mw = divertor.divwade(
            rmajor,
            rminor,
            aspect,
            bt,
            bp,
            p_plasma_separatrix_mw,
            f_div_flux_expansion,
            nesep,
            deg_div_field_plate,
            rad_fraction_sol,
            f_p_div_lower,
            False,
        )

        assert pflux_div_heat_load_mw == pytest.approx(expected_pflux_div_heat_load_mw)
