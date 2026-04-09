from typing import Any, NamedTuple

import pytest

import process.models.superconductors as superconductors


class IterscParam(NamedTuple):
    thelium: Any = None

    bmax: Any = None

    strain: Any = None

    bc20max: Any = None

    tc0max: Any = None

    expected_jcrit: Any = None

    expected_bcrit: Any = None

    expected_tcrit: Any = None


@pytest.mark.parametrize(
    "iterscparam",
    (
        IterscParam(
            thelium=4.75,
            bmax=13.008974843466492,
            strain=0.001601605753441172,
            bc20max=32.969999999999999,
            tc0max=16.059999999999999,
            expected_jcrit=692348194.774593,
            expected_bcrit=27.092853296363597,
            expected_tcrit=11.338458919718571,
        ),
        IterscParam(
            thelium=6.2510000000000003,
            bmax=13.008974843466492,
            strain=0.001601605753441172,
            bc20max=32.969999999999999,
            tc0max=16.059999999999999,
            expected_jcrit=495889332.08959526,
            expected_bcrit=24.442648486388464,
            expected_tcrit=11.338458919718571,
        ),
    ),
)
def test_itersc(iterscparam):
    """
    Automatically generated Regression Unit Test for itersc.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param iterscparam: the data used to mock and assert in this test.
    :type iterscparam: iterscparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    jcrit, bcrit, tcrit = superconductors.itersc(
        temp_conductor=iterscparam.thelium,
        b_conductor=iterscparam.bmax,
        strain=iterscparam.strain,
        b_c20max=iterscparam.bc20max,
        temp_c0max=iterscparam.tc0max,
    )

    assert jcrit == pytest.approx(iterscparam.expected_jcrit)

    assert bcrit == pytest.approx(iterscparam.expected_bcrit)

    assert tcrit == pytest.approx(iterscparam.expected_tcrit)


class JcritNbtiParam(NamedTuple):
    temperature: Any = None

    bmax: Any = None

    c0: Any = None

    bc20max: Any = None

    tc0max: Any = None

    expected_jcrit: Any = None

    expected_tcrit: Any = None


@pytest.mark.parametrize(
    "jcritnbtiparam",
    (
        JcritNbtiParam(
            temperature=4.75,
            bmax=8.0517923638507547,
            c0=10000000000,
            bc20max=15,
            tc0max=9.3000000000000007,
            expected_jcrit=906668274.04561484,
            expected_tcrit=5.9060082696285683,
        ),
        JcritNbtiParam(
            temperature=6,
            bmax=8.0517923638507547,
            c0=10000000000,
            bc20max=15,
            tc0max=9.3000000000000007,
            expected_jcrit=-73718607.547511846,
            expected_tcrit=5.9060082696285683,
        ),
    ),
)
def test_jcrit_nbti(jcritnbtiparam):
    """
    Automatically generated Regression Unit Test for jcrit_nbti.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param jcritnbtiparam: the data used to mock and assert in this test.
    :type jcritnbtiparam: jcritnbtiparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    jcrit, tcrit = superconductors.jcrit_nbti(
        temp_conductor=jcritnbtiparam.temperature,
        b_conductor=jcritnbtiparam.bmax,
        c0=jcritnbtiparam.c0,
        b_c20max=jcritnbtiparam.bc20max,
        temp_c0max=jcritnbtiparam.tc0max,
    )

    assert jcrit == pytest.approx(jcritnbtiparam.expected_jcrit)

    assert tcrit == pytest.approx(jcritnbtiparam.expected_tcrit)


def test_jcrit_rebco():
    jcrit_rebco, validity = superconductors.jcrit_rebco(4.75, 7.0)

    assert jcrit_rebco == pytest.approx(55870234414.171684)
    assert validity


def test_current_sharing_rebco():
    assert superconductors.current_sharing_rebco(11.0, 2e7) == pytest.approx(
        71.28702697627514
    )


def test_bi2212():
    jcrit, tmarg = superconductors.bi2212(7.0, 2e7, 4.75, 0.2)

    assert jcrit == pytest.approx(174017403.16041547)
    assert tmarg == pytest.approx(13.750991122745397)


def test_gl_nbti():
    jcrit, tcrit, bcrit = superconductors.gl_nbti(4.75, 7.0, 2, 9.5, 13.75)

    assert jcrit == pytest.approx(2551683055.6511745)
    assert tcrit == pytest.approx(7.277374792835339)
    assert bcrit == pytest.approx(13.662161361675887)


def test_wstsc():
    jcrit, bcrit, tcrit = superconductors.western_superconducting_nb3sn(
        4.75, 27.0, 0.001, 30.0, 25.0
    )

    assert jcrit == pytest.approx(195513.0673058944)
    assert bcrit == pytest.approx(27.329369840368482)
    assert tcrit == pytest.approx(5.170678992915718)


def test_gl_rebco():
    jcrit, bcrit, tcrit = superconductors.gl_rebco(4.75, 7.0, 2, 30.0, 25.0)

    assert jcrit == pytest.approx(14527.765708690296)
    assert bcrit == pytest.approx(9.439350824747793)
    assert tcrit == pytest.approx(24.66989137698065)


def test_hijc_rebco():
    jcrit, bcrit, tcrit = superconductors.hijc_rebco(
        temp_conductor=4.75,
        b_conductor=7.0,
        b_c20max=30.0,
        t_c0=25.0,
        dr_hts_tape=4.0e-3,
        dx_hts_tape_rebco=1.0e-6,
        dx_hts_tape_total=6.5e-5,
    )

    assert jcrit == pytest.approx(111046017.5)
    assert bcrit == pytest.approx(22.335736687814954)
    assert tcrit == pytest.approx(24.999125)


@pytest.mark.parametrize(
    "dia_croco_strand, dx_croco_strand_copper, dx_hts_tape_rebco, dx_hts_tape_copper, dx_hts_tape_hastelloy, expected",
    [
        (
            0.010,  # 10 mm
            0.001,  # 1 mm
            1e-6,  # 1 um
            2e-6,  # 2 um
            3e-6,  # 3 um
            (
                0.008,  # dia_croco_strand_tape_region
                pytest.approx(959.3950997769347, rel=1e-3),  # n_croco_strand_hts_tapes
                pytest.approx(
                    3.8934279435385194e-05, rel=1e-3
                ),  # a_croco_strand_copper_total
                pytest.approx(
                    1.5989918329615573e-05, rel=1e-3
                ),  # a_croco_strand_hastelloy
                pytest.approx(1.8285645798205533e-05, rel=1e-3),  # a_croco_strand_solder
                pytest.approx(5.329972776538525e-06, rel=1e-3),  # a_croco_strand_rebco
                pytest.approx(7.85398e-5, rel=1e-3),  # croco_strand_area
                pytest.approx(5.5556e-3, rel=1e-3),  # dr_hts_tape
            ),
        ),
        (
            0.0054,  # baseline diameter
            0.0005,  # baseline copper thickness
            1e-6,
            2e-6,
            3e-6,
            (
                0.0044,
                pytest.approx(527.6673048773141, rel=1e-6),
                pytest.approx(1.0921535531100803e-05, rel=1e-6),
                pytest.approx(4.836950294708712e-06, rel=1e-6),
                pytest.approx(5.531407853957174e-06, rel=1e-6),
                pytest.approx(1.612316764902904e-06, rel=1e-6),
                pytest.approx(2.2902210444669593e-05, rel=1e-6),
                pytest.approx(0.0030555555555555553, rel=1e-6),
            ),
        ),
    ],
)
def test_calculate_croco_cable_geometry(
    dia_croco_strand,
    dx_croco_strand_copper,
    dx_hts_tape_rebco,
    dx_hts_tape_copper,
    dx_hts_tape_hastelloy,
    expected,
):
    result = superconductors.calculate_croco_cable_geometry(
        dia_croco_strand,
        dx_croco_strand_copper,
        dx_hts_tape_rebco,
        dx_hts_tape_copper,
        dx_hts_tape_hastelloy,
    )
    for r, e in zip(result, expected, strict=True):
        assert r == e
