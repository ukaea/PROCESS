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

        expected_fncmass = 310716.52923547616
        expected_aintmass = 5829865.436088617
        expected_clgsmass = 2018975.3864451263
        expected_coldmass = 48937690.168336436
        expected_gsm = 1685092.6111564008

        fncmass, aintmass, clgsmass, coldmass, gsm = structure.structure(
            ai=17721306.969367817,
            r0=8.8901,
            a=2.8677741935483869,
            akappa=1.848,
            b0=5.3292,
            i_tf_sup=1,
            i_pf_conductor=0,
            tf_h_width=15.337464674334223,
            tfhmax=9.0730900215620327,
            shldmass=2294873.8131476026,
            dvrtmass=43563.275828777645,
            pfmass=5446188.2481440185,
            tfmass=21234909.756419446,
            m_fw_total=224802.80270851994,
            blmass=3501027.3252278985,
            m_fw_blkt_div_coolant_total=1199.6389920083477,
            dewmass=16426726.727684354,
            output=False,
        )

        assert fncmass == pytest.approx(expected_fncmass)
        assert aintmass == pytest.approx(expected_aintmass)
        assert clgsmass == pytest.approx(expected_clgsmass)
        assert coldmass == pytest.approx(expected_coldmass)
        assert gsm == pytest.approx(expected_gsm)
