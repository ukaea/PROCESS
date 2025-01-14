"""Unit tests for plasma_geometry.f90."""

from typing import Any, NamedTuple

import pytest

from process.plasma_geometry import PlasmaGeom


@pytest.fixture
def plasma():
    """Fixture to create a PFCoil object.

    :return: an instance of PFCoil
    :rtype: process.pfcoil.PFCoil
    """
    plasma = PlasmaGeom()

    return plasma


class XparamParam(NamedTuple):
    a: Any = None

    kap: Any = None

    tri: Any = None

    expected_xi: Any = None

    expected_thetai: Any = None

    expected_xo: Any = None

    expected_thetao: Any = None


@pytest.mark.parametrize(
    "xparamparam",
    (
        XparamParam(
            a=2.8677741935483869,
            kap=1.8480000000000001,
            tri=0.5,
            expected_xi=10.510690667870968,
            expected_thetai=0.52847258461252744,
            expected_xo=5.4154130183225808,
            expected_thetao=1.3636548755403939,
        ),
        XparamParam(
            a=2.8677741935483869,
            kap=1.8480000000000001,
            tri=0.5,
            expected_xi=10.510690667870968,
            expected_thetai=0.52847258461252744,
            expected_xo=5.4154130183225808,
            expected_thetao=1.3636548755403939,
        ),
    ),
)
def test_xparam(xparamparam, monkeypatch, plasma):
    """
    Automatically generated Regression Unit Test for xparam.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param xparamparam: the data used to mock and assert in this test.
    :type xparamparam: xparamparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    xi, thetai, xo, thetao = plasma.xparam(
        a=xparamparam.a, kap=xparamparam.kap, tri=xparamparam.tri
    )

    assert xi == pytest.approx(xparamparam.expected_xi)

    assert thetai == pytest.approx(xparamparam.expected_thetai)

    assert xo == pytest.approx(xparamparam.expected_xo)

    assert thetao == pytest.approx(xparamparam.expected_thetao)


@pytest.mark.parametrize(
    "a, kap, tri, expected_perim",
    [
        (
            2.8677741935483869,
            1.8480000000000001,
            0.5,
            25.8787324576261,
        )
    ],
)
def test_perim(a, kap, tri, expected_perim, plasma):
    """Tests `perim` function.

    :param a: test asset passed to the routine representing the plasma minor radius, in meters.
    :type a: float

    :param kap: test asset passed to the routine representing the plasma separatrix elongation.
    :type kap: float

    :param tri: test asset passed to the routine representing the plasma separatrix triangularity.
    :type tri: float

    :param expected_perim: expected result of the function.
    :type expected_perim: float
    """
    perim = plasma.perim(a, kap, tri)

    assert pytest.approx(perim) == expected_perim


@pytest.mark.parametrize(
    "rmajor, rminor, xi, thetai, xo, thetao, expected_xvol",
    [
        (
            9.2995201822511735,
            2.9998452200810237,
            10.261919050584332,
            0.54748563700358688,
            5.4205364969154601,
            1.4001019213417263,
            2694.44683912,
        )
    ],
)
def test_xvol(rmajor, rminor, xi, thetai, xo, thetao, expected_xvol, plasma):
    """Tests `xvol` function.
    :param rmajor: test asset passed to the routine representing the plasma major radius (m).
    :type rmajor: float

    :param rminor: test asset passed to the routine representing the plasma minor radius (m).
    :type rminor: float

    :param xi: test asset passed to the routine representing the radius of arc describing inboard surface (m).
    :type xi: float

    :param thetai: test asset passed to the routine representing the half-angle of arc describing inboard surface.
    :type thetai: float

    :param xo: test asset passed to the routine representing the radius of arc describing outboard surface (m).
    :type xo: float

    :param thetao: test asset passed to the routine representing the half-angle of arc describing outboard surface.
    :type thetao: float

    :param expected_xvol: expected result of the function.
    :type expected_xvol: float
    """
    xvol = plasma.xvol(rmajor, rminor, xi, thetai, xo, thetao)

    assert pytest.approx(xvol) == expected_xvol


@pytest.mark.parametrize(
    "xi, thetai, xo ,thetao, expected_xsecta",
    [
        (
            10.261919050584332,
            0.54748563700358688,
            5.4205364969154601,
            1.4001019213417263,
            47.069149087374726,
        )
    ],
)
def test_xsecta(xi, thetai, xo, thetao, expected_xsecta, plasma):
    """Tests `xsecta` function.
    :param xi: test asset passed to the routine representing the radius of arc describing inboard surface (m).
    :type xi: float

    :param thetai: test asset passed to the routine representing the half-angle of arc describing inboard surface.
    :type thetai: float

    :param xo: test asset passed to the routine representing the radius of arc describing outboard surface (m).
    :type xo: float

    :param thetao: test asset passed to the routine representing the half-angle of arc describing outboard surface.
    :type thetao: float

    :param expected_xsecta: expected result of the function.
    :type expected_xsecta: float
    """
    xsecta = plasma.xsecta(xi, thetai, xo, thetao)

    assert pytest.approx(xsecta) == expected_xsecta


@pytest.mark.parametrize(
    "r, a, kap, tri, expected_fvol",
    [
        (
            9.2995201822511735,
            2.8677741935483869,
            1.8480000000000001,
            0.5,
            2540.3858087477,
        )
    ],
)
def test_fvol(r, a, kap, tri, expected_fvol, plasma):
    """Tests `fvol` function.
    :param r: test asset passed to the routine representing the plasma major radius, in meters.
    :type r: float

    :param a: test asset passed to the routine representing the plasma minor radius, in meters.
    :type a: float

    :param kap: test asset passed to the routine representing the plasma separatrix elongation.
    :type kap: float

    :param tri: test asset passed to the routine representing the plasma separatrix triangularity.
    :type tri: float

    :param expected_fvol: expected result of the function.
    :type expected_fvol: float

    """
    fvol = plasma.fvol(r, a, kap, tri)

    assert pytest.approx(fvol) == expected_fvol


@pytest.mark.parametrize(
    "a, kap, tri, expected_xsect0",
    [
        (
            2.8677741935483869,
            1.8480000000000001,
            0.5,
            44.367959142766466,
        )
    ],
)
def test_xsect0(a, kap, tri, expected_xsect0, plasma):
    """Tests `xsect0` function.

    :param a: test asset passed to the routine representing the plasma minor radius, in meters.
    :type a: float

    :param kap: test asset passed to the routine representing the plasma separatrix elongation.
    :type kap: float

    :param tri: test asset passed to the routine representing the plasma separatrix triangularity.
    :type tri: float

    :param expected_xsect0: expected result of the function.
    :type expected_xsect0: float

    """
    xsect0 = plasma.xsect0(a, kap, tri)

    assert pytest.approx(xsect0) == expected_xsect0


class XsurfParam(NamedTuple):
    rmajor: Any = None

    rminor: Any = None

    xi: Any = None

    thetai: Any = None

    xo: Any = None

    thetao: Any = None

    expected_xsi: Any = None

    expected_xso: Any = None


@pytest.mark.parametrize(
    "xsurfparam",
    (
        XsurfParam(
            rmajor=8.8901000000000003,
            rminor=2.8677741935483869,
            xi=10.510690667870968,
            thetai=0.52847258461252744,
            xo=5.4154130183225808,
            thetao=1.3636548755403939,
            expected_xsi=454.0423505329922,
            expected_xso=949.22962703393853,
        ),
        XsurfParam(
            rmajor=8.8901000000000003,
            rminor=2.8677741935483869,
            xi=10.510690667870968,
            thetai=0.52847258461252744,
            xo=5.4154130183225808,
            thetao=1.3636548755403939,
            expected_xsi=454.0423505329922,
            expected_xso=949.22962703393853,
        ),
    ),
)
def test_xsurf(xsurfparam, monkeypatch, plasma):
    """
    Automatically generated Regression Unit Test for xsurf.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param xsurfparam: the data used to mock and assert in this test.
    :type xsurfparam: xsurfparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    xsi, xso = plasma.xsurf(
        rmajor=xsurfparam.rmajor,
        rminor=xsurfparam.rminor,
        xi=xsurfparam.xi,
        thetai=xsurfparam.thetai,
        xo=xsurfparam.xo,
        thetao=xsurfparam.thetao,
    )

    assert xsi == pytest.approx(xsurfparam.expected_xsi)

    assert xso == pytest.approx(xsurfparam.expected_xso)


class SurfaParam(NamedTuple):
    a: Any = None

    r: Any = None

    k: Any = None

    d: Any = None

    expected_sa: Any = None

    expected_so: Any = None


@pytest.mark.parametrize(
    "surfaparam",
    (
        SurfaParam(
            a=0.97142857142857153,
            r=1.7000000000000002,
            k=2.75,
            d=0.5,
            expected_sa=117.06185446474473,
            expected_so=86.489845577107758,
        ),
        SurfaParam(
            a=0.97142857142857153,
            r=1.7000000000000002,
            k=2.75,
            d=0.5,
            expected_sa=117.06185446474473,
            expected_so=86.489845577107758,
        ),
    ),
)
def test_surfa(surfaparam, monkeypatch, plasma):
    """
    Automatically generated Regression Unit Test for surfa.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param surfaparam: the data used to mock and assert in this test.
    :type surfaparam: surfaparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    sa, so = plasma.surfa(
        a=surfaparam.a, r=surfaparam.r, k=surfaparam.k, d=surfaparam.d
    )

    assert sa == pytest.approx(surfaparam.expected_sa)

    assert so == pytest.approx(surfaparam.expected_so)
