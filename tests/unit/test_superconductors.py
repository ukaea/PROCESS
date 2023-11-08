import pytest
from typing import NamedTuple, Any

import process.superconductors as superconductors


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
        thelium=iterscparam.thelium,
        bmax=iterscparam.bmax,
        strain=iterscparam.strain,
        bc20max=iterscparam.bc20max,
        tc0max=iterscparam.tc0max,
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
        temperature=jcritnbtiparam.temperature,
        bmax=jcritnbtiparam.bmax,
        c0=jcritnbtiparam.c0,
        bc20max=jcritnbtiparam.bc20max,
        tc0max=jcritnbtiparam.tc0max,
    )

    assert jcrit == pytest.approx(jcritnbtiparam.expected_jcrit)

    assert tcrit == pytest.approx(jcritnbtiparam.expected_tcrit)


def test_jcrit_rebco():
    jcrit_rebco, validity = superconductors.jcrit_rebco(4.75, 7.0)

    assert jcrit_rebco == pytest.approx(55870234414.171684)
    assert validity


def test_current_sharing_rebco():
    assert superconductors.current_sharing_rebco(7.0, 2e7) == pytest.approx(
        75.76286550648135
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
    jcrit, bcrit, tcrit = superconductors.wstsc(4.75, 27.0, 0.001, 30.0, 25.0)

    assert jcrit == pytest.approx(195513.0673058944)
    assert bcrit == pytest.approx(27.329369840368482)
    assert tcrit == pytest.approx(5.170678992915718)


def test_gl_rebco():
    jcrit, bcrit, tcrit = superconductors.gl_rebco(4.75, 7.0, 2, 30.0, 25.0)

    assert jcrit == pytest.approx(14527.765708690296)
    assert bcrit == pytest.approx(9.439350824747793)
    assert tcrit == pytest.approx(24.66989137698065)


def test_hijc_rebco():
    jcrit, bcrit, tcrit = superconductors.hijc_rebco(4.75, 7.0, 2, 30.0, 25.0)

    assert jcrit == pytest.approx(111046017.5)
    assert bcrit == pytest.approx(22.335736687814954)
    assert tcrit == pytest.approx(24.999125)
