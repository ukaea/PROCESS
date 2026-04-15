"""Unit tests for plasma_geometry.f90."""

from typing import Any, NamedTuple

import pytest

import process.models.physics.plasma_geometry as pg
from process.models.physics.plasma_geometry import PlasmaGeom


@pytest.fixture
def plasma():
    """Fixture to create a PFCoil object.

    :return: an instance of PFCoil
    :rtype: process.pfcoil.PFCoil
    """

    return PlasmaGeom()


class PlasmaAnglesArcsParam(NamedTuple):
    a: Any = None

    kap: Any = None

    tri: Any = None

    expected_xi: Any = None

    expected_thetai: Any = None

    expected_xo: Any = None

    expected_thetao: Any = None


@pytest.mark.parametrize(
    "plasmaanglesarcsparam",
    [
        PlasmaAnglesArcsParam(
            a=2.8677741935483869,
            kap=1.8480000000000001,
            tri=0.5,
            expected_xi=10.510690667870968,
            expected_thetai=0.52847258461252744,
            expected_xo=5.4154130183225808,
            expected_thetao=1.3636548755403939,
        ),
    ],
)
def test_plasma_angles_arcs(
    plasmaanglesarcsparam: PlasmaAnglesArcsParam, monkeypatch, plasma
):
    """
    Automatically generated Regression Unit Test for plasma_angles_arcs().

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    Parameters
    ----------
    plasmaanglesarcsparam:
        the data used to mock and assert in this test.
    monkeypatch:
        pytest fixture used to mock module/class variables
    """

    xi, thetai, xo, thetao = plasma.plasma_angles_arcs(
        a=plasmaanglesarcsparam.a,
        kappa=plasmaanglesarcsparam.kap,
        triang=plasmaanglesarcsparam.tri,
    )

    assert xi == pytest.approx(plasmaanglesarcsparam.expected_xi)

    assert thetai == pytest.approx(plasmaanglesarcsparam.expected_thetai)

    assert xo == pytest.approx(plasmaanglesarcsparam.expected_xo)

    assert thetao == pytest.approx(plasmaanglesarcsparam.expected_thetao)


@pytest.mark.parametrize(
    ("a", "kap", "tri", "expected_perim"),
    [
        (
            2.8677741935483869,
            1.8480000000000001,
            0.5,
            25.8787324576261,
        )
    ],
)
def test_perim(a: float, kap: float, tri: float, expected_perim: float):
    """Tests `perim` function.

    Parameters
    ----------
    a:
        the plasma minor radius, in meters.
    kap:
        the plasma separatrix elongation.
    tri:
        the plasma separatrix triangularity.
    expected_perim:
        expected result of the function.
    """
    perim = pg.perim(a, kap, tri)

    assert pytest.approx(perim) == expected_perim


@pytest.mark.parametrize(
    ("rmajor", "rminor", "xi", "thetai", "xo", "thetao", "expected_plasma_volume"),
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
def test_plasma_volume(
    rmajor: float,
    rminor: float,
    xi: float,
    thetai: float,
    xo: float,
    thetao: float,
    expected_plasma_volume: float,
    plasma: float,
):
    """Tests `plasma_volume()` function.

    Parameters
    ----------
    rmajor:
        the plasma major radius (m).
    rminor:
        the plasma minor radius (m).
    xi:
        the radius of arc describing inboard surface (m).
    thetai:
        the half-angle of arc describing inboard surface.
    xo:
        the radius of arc describing outboard surface (m).
    thetao:
        the half-angle of arc describing outboard surface.
    expected_plasma_volume:
        expected result of the function.
    """
    plasma_volume = plasma.plasma_volume(rmajor, rminor, xi, thetai, xo, thetao)

    assert pytest.approx(plasma_volume) == expected_plasma_volume


@pytest.mark.parametrize(
    ("xi", "thetai", "xo", "thetao", "expected_plasma_cross_section"),
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
def test_plasma_cross_section(
    xi, thetai, xo, thetao, expected_plasma_cross_section, plasma
):
    """Tests `plasma_cross_section()` function.
    xi: the radius of arc describing inboard surface (m).
    :type xi: float

    thetai: the half-angle of arc describing inboard surface.
    :type thetai: float

    xo: the radius of arc describing outboard surface (m).
    :type xo: float

    thetao: the half-angle of arc describing outboard surface.
    :type thetao: float

    expected_plasma_cross_section: expected result of the function.
    :type expected_plasma_cross_section: float
    """
    plasma_cross_section = plasma.plasma_cross_section(xi, thetai, xo, thetao)

    assert pytest.approx(plasma_cross_section) == expected_plasma_cross_section


@pytest.mark.parametrize(
    ("r", "a", "kap", "tri", "expected_fvol"),
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
def test_fvol(r: float, a: float, kap: float, tri: float, expected_fvol: float):
    """Tests `fvol` function.

    Parameters
    ----------
    r:
        the plasma major radius, in meters.
    a:
        the plasma minor radius, in meters.
    kap:
        the plasma separatrix elongation.
    tri:
        test asset passed to the routine
        representing the plasma separatrix triangularity.
    expected_fvol: expected result of the function.
    """
    fvol = pg.fvol(r, a, kap, tri)

    assert pytest.approx(fvol) == expected_fvol


@pytest.mark.parametrize(
    ("a", "kap", "tri", "expected_xsect0"),
    [
        (
            2.8677741935483869,
            1.8480000000000001,
            0.5,
            44.367959142766466,
        )
    ],
)
def test_xsect0(a: float, kap: float, tri: float, expected_xsect0: float):
    """Tests `xsect0` function.

    Parameters
    ----------
    a:
        the plasma minor radius, in meters.
    kap:
        the plasma separatrix elongation.
    tri:
        test asset passed to the routine
        representing the plasma separatrix triangularity.
    expected_xsect0:
        expected result of the function.
    """
    xsect0 = pg.xsect0(a, kap, tri)

    assert pytest.approx(xsect0) == expected_xsect0


class SurfaceAreaParam(NamedTuple):
    rmajor: Any = None

    rminor: Any = None

    xi: Any = None

    thetai: Any = None

    xo: Any = None

    thetao: Any = None

    expected_xsi: Any = None

    expected_xso: Any = None


@pytest.mark.parametrize(
    "surfaceareaparam",
    [
        SurfaceAreaParam(
            rmajor=8.8901000000000003,
            rminor=2.8677741935483869,
            xi=10.510690667870968,
            thetai=0.52847258461252744,
            xo=5.4154130183225808,
            thetao=1.3636548755403939,
            expected_xsi=454.0423505329922,
            expected_xso=949.22962703393853,
        ),
    ],
)
def test_plasma_surface_area(surfaceareaparam: SurfaceAreaParam, monkeypatch, plasma):
    """
    Automatically generated Regression Unit Test for plasma_surface_area().

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    Parameters
    ----------
    surfaceareaparam:
        the data used to mock and assert in this test.
    monkeypatch:
        pytest fixture used to mock module/class variables
    """

    xsi, xso = plasma.plasma_surface_area(
        rmajor=surfaceareaparam.rmajor,
        rminor=surfaceareaparam.rminor,
        xi=surfaceareaparam.xi,
        thetai=surfaceareaparam.thetai,
        xo=surfaceareaparam.xo,
        thetao=surfaceareaparam.thetao,
    )

    assert xsi == pytest.approx(surfaceareaparam.expected_xsi)

    assert xso == pytest.approx(surfaceareaparam.expected_xso)


class SurfaParam(NamedTuple):
    a: Any = None

    r: Any = None

    k: Any = None

    d: Any = None

    expected_sa: Any = None

    expected_so: Any = None


@pytest.mark.parametrize(
    "surfaparam",
    [
        SurfaParam(
            a=0.97142857142857153,
            r=1.7000000000000002,
            k=2.75,
            d=0.5,
            expected_sa=117.06185446474473,
            expected_so=86.489845577107758,
        ),
    ],
)
def test_surfa(surfaparam: SurfaParam, monkeypatch):
    """
    Automatically generated Regression Unit Test for surfa.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    Parameters
    ----------
    surfaparam:
        the data used to mock and assert in this test.
    monkeypatch:
        pytest fixture used to mock module/class variables
    """

    sa, so = pg.surfa(a=surfaparam.a, r=surfaparam.r, k=surfaparam.k, d=surfaparam.d)

    assert sa == pytest.approx(surfaparam.expected_sa)

    assert so == pytest.approx(surfaparam.expected_so)
