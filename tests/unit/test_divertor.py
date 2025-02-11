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
    def test_divert(self, divertor):
        """Test the divert subroutine.

        Uses test data from the first call of this subroutine by baseline 2018.

        :param divertor: fixture containing an initialised `Divertor` object
        :type divertor: tests.unit.test_divertor.divertor (functional fixture)
        """
        adas = 0.052617908173833536
        m_ions_total_amu = 2.5
        anginc = 0.262
        c1div = 0.45
        c2div = -7
        c3div = 0.54
        c4div = -3.6
        c5div = 0.7
        delld = 1
        delne = 0.34294618459618942
        fdfs = 10
        fififi = 0.004
        frgd = 0.33822011777067068
        frrp = 0.4
        minstang = 0.5
        omegan = 1
        qdiv = 0.21095478076086721
        pdiv = 100.12226393133213
        rbpbt = 0.067124908423896915
        rconl = 0.48318885267795408
        rmaj = 8.8901
        rsrd = 1.2262681695075397
        tconl = 78.118065070110475
        xpara = 600
        xperp = 2.66

        expected_delta = 1.9073843835398801e-12
        expected_delw = 0.021077486208183741
        expected_dendiv = 3.7696864807168677
        expected_densin = 19.136765450021365
        expected_gamdt = 6663.1765788767079
        expected_lamp = 0.51210900896428913
        expected_omlarg = 1.2535646368097362
        expected_ppdiv = 4.2920738515450783
        expected_ppdivr = 4.3475484299911091
        expected_ptpdiv = 5.4674204940995681
        expected_tdiv = 27.70537887577682
        expected_tsep = 305.08793598453673

        (
            delta,
            delw,
            dendiv,
            densin,
            gamdt,
            lamp,
            omlarg,
            ppdiv,
            ppdivr,
            ptpdiv,
            tdiv,
            tsep,
        ) = divertor.divert(
            adas,
            m_ions_total_amu,
            anginc,
            delne,
            c1div,
            c2div,
            c3div,
            c4div,
            c5div,
            delld,
            fdfs,
            fififi,
            frgd,
            frrp,
            minstang,
            omegan,
            qdiv,
            pdiv,
            rbpbt,
            rconl,
            rmaj,
            rsrd,
            tconl,
            xpara,
            xperp,
        )

        assert delta == pytest.approx(expected_delta)
        assert delw == pytest.approx(expected_delw)
        assert dendiv == pytest.approx(expected_dendiv)
        assert densin == pytest.approx(expected_densin)
        assert gamdt == pytest.approx(expected_gamdt)
        assert lamp == pytest.approx(expected_lamp)
        assert omlarg == pytest.approx(expected_omlarg)
        assert ppdiv == pytest.approx(expected_ppdiv)
        assert ppdivr == pytest.approx(expected_ppdivr)
        assert ptpdiv == pytest.approx(expected_ptpdiv)
        assert tdiv == pytest.approx(expected_tdiv)
        assert tsep == pytest.approx(expected_tsep)

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
        pdivt = 7.7197999809272062
        monkeypatch.setattr(dv, "i_hldiv", 1)

        expected_hldiv = 0.087770426974167357

        hldiv = divertor.divtart(
            rmajor,
            rminor,
            triang,
            dr_fw_plasma_gap_inboard,
            dz_xpoint_divertor,
            pdivt,
            False,
        )

        assert hldiv == pytest.approx(expected_hldiv)

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
        pdivt = 1.0e2
        flux_exp = 2
        nesep = 1.0e19
        beta_div = 5.0
        rad_fraction_sol = 8.0e-1
        ftar = 1.0

        expected_hldiv = 0.58898578

        hldiv = divertor.divwade(
            rmajor,
            rminor,
            aspect,
            bt,
            bp,
            pdivt,
            flux_exp,
            nesep,
            beta_div,
            rad_fraction_sol,
            ftar,
            False,
        )

        assert hldiv == pytest.approx(expected_hldiv)

    def test_erprcy(self, divertor):
        """Test the erprcy subroutine.

        Uses test data from the first call of this subroutine by baseline 2018.

        :param divertor: fixture containing an initialised `Divertor` object
        :type divertor: tests.unit.test_divertor.divertor (functional fixture)
        """
        tdiv = 150
        ndiv = 0.38105131621798821

        expected_erprcy = 24.949836803003997

        erprcy = divertor.erprcy(tdiv, ndiv)

        assert erprcy == pytest.approx(expected_erprcy)

    def test_ftdiv(self, divertor):
        """Test the ftdiv subroutine.

        Uses test data from the first call of this subroutine by baseline 2018.

        :param divertor: fixture containing an initialised `Divertor` object
        :type divertor: tests.unit.test_divertor.divertor (functional fixture)
        """
        m_ions_total_amu = 2.5
        coefl = 1.3854518853592164
        delne = 0.34294618459618942
        fififi = 0.0040000000000000001
        omegan = 1
        omlarg = 1.2535646368097362
        qdiv = 0.21095478076086721
        tconl = 78.118065070110475
        xpara = 600
        xperp = 2.6600000000000001
        xx = 0.90000000000000002
        yy = 150

        expected_ftdiv = 140.7905577701393

        ftdiv = divertor.ftdiv(
            m_ions_total_amu,
            coefl,
            delne,
            fififi,
            omegan,
            omlarg,
            qdiv,
            tconl,
            xpara,
            xperp,
            xx,
            yy,
        )

        assert ftdiv == pytest.approx(expected_ftdiv)

    def test_ftpts(self, divertor):
        """Test the ftpts subroutine.

        Uses test data from the first call of this subroutine by baseline 2018.

        :param divertor: fixture containing an initialised `Divertor` object
        :type divertor: tests.unit.test_divertor.divertor (functional fixture)
        """
        m_ions_total_amu = 2.5
        coefl = 1.3854518853592164
        delne = 0.34294618459618942
        fififi = 0.0040000000000000001
        omegan = 1
        omlarg = 1.2535646368097362
        qdiv = 0.21095478076086721
        tconl = 78.118065070110475
        xpara = 600
        xperp = 2.6600000000000001
        xx = 0.90000000000000002
        yy = 150

        expected_ftpts = -7.962899353691732

        ftpts = divertor.ftpts(
            m_ions_total_amu,
            coefl,
            delne,
            fififi,
            omegan,
            omlarg,
            qdiv,
            tconl,
            xpara,
            xperp,
            xx,
            yy,
        )

        assert ftpts == pytest.approx(expected_ftpts)

    def test_gammash(self, divertor):
        """Test the gammash subroutine.

        Uses test data from the first call of this subroutine by baseline 2018.

        :param divertor: fixture containing an initialised `Divertor` object
        :type divertor: tests.unit.test_divertor.divertor (functional fixture)
        """
        gcoef = 0.0040000000000000001
        tdiv = 150

        expected_gammash = 8.155694305511572

        gammash = divertor.gammash(gcoef, tdiv)

        assert gammash == pytest.approx(expected_gammash)
