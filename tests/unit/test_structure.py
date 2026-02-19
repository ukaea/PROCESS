"""Unit tests for structure.py"""

import pytest

from process.models.structure import Structure


@pytest.fixture
def structure():
    """Provides Structure object for testing.

    :returns: initialised Structure object
    :rtype: process.structure.Structure
    """
    return Structure()


class TestStructure:
    def test_structure(self, structure):
        """Tests the structure subroutine"""
        ai: float = 17721306.969367817
        r0: float = 8.8901
        a: float = 2.8677741935483869
        akappa: float = 1.848
        b0: float = 5.3292
        tf_h_width: float = 15.337464674334223
        tfhmax: float = 9.0730900215620327
        shldmass: float = 2294873.8131476026
        dvrtmass: float = 43563.275828777645
        pfmass: float = 5446188.2481440185
        tfmass: float = 21234909.756419446
        m_fw_total: float = 224802.80270851994
        blmass: float = 3501027.3252278985
        m_fw_blkt_div_coolant_total: float = 1199.6389920083477
        dewmass: float = 16426726.727684354
        i_tf_sup: int = 1
        i_pf_conductor: int = 0

        expected_fncmass: float = 310716.52923547616
        expected_aintmass: float = 5829865.436088617
        expected_clgsmass: float = 2018975.3864451263
        expected_coldmass: float = 48937690.168336436
        expected_gsm: float = 1685092.6111564008

        fncmass, aintmass, clgsmass, coldmass, gsm = structure.structure(
            ai,
            r0,
            a,
            akappa,
            b0,
            i_tf_sup,
            i_pf_conductor,
            tf_h_width,
            tfhmax,
            shldmass,
            dvrtmass,
            pfmass,
            tfmass,
            m_fw_total,
            blmass,
            m_fw_blkt_div_coolant_total,
            dewmass,
            output=False,
        )

        assert fncmass == pytest.approx(expected_fncmass)
        assert aintmass == pytest.approx(expected_aintmass)
        assert clgsmass == pytest.approx(expected_clgsmass)
        assert coldmass == pytest.approx(expected_coldmass)
        assert gsm == pytest.approx(expected_gsm)
