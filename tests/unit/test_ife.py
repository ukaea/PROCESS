"""Unit tests for ife."""
from typing import NamedTuple, Any

import pytest
import numpy

from process.ife import IFE
from process.availability import Availability
from process.costs import Costs
from process.fortran import ife_module, build_variables, ife_variables


@pytest.fixture
def ife():
    """Provides IFE object for testing.

    :returns: initialised IFE object
    :rtype: process.ife.IFE
    """
    return IFE(Availability(), Costs())


def test_ifetgt(monkeypatch):
    """Test ifetgt.

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    # Mock module variables
    # Repetition Rate (Hz)
    monkeypatch.setattr(ife_variables, "reprat", 4.0)
    # IFE target factory power at 6 Hz repetition rate
    monkeypatch.setattr(ife_variables, "ptargf", 2.0)
    monkeypatch.setattr(ife_variables, "tfacmw", 0)

    ife_module.ifetgt()
    assert ife_variables.tfacmw == pytest.approx(1.506, abs=0.001)


class SombldParam(NamedTuple):
    fwarea: Any = None
    ifetyp: Any = None
    chrad: Any = None
    r1: Any = None
    fwdr: Any = None
    r2: Any = None
    v1dr: Any = None
    r3: Any = None
    bldr: Any = None
    r4: Any = None
    v2dr: Any = None
    r5: Any = None
    shdr: Any = None
    r6: Any = None
    v3dr: Any = None
    r7: Any = None
    zl7: Any = None
    v3dzl: Any = None
    zl6: Any = None
    shdzl: Any = None
    zl5: Any = None
    v2dzl: Any = None
    zl4: Any = None
    bldzl: Any = None
    zl3: Any = None
    v1dzl: Any = None
    zl2: Any = None
    fwdzl: Any = None
    zl1: Any = None
    chdzl: Any = None
    chdzu: Any = None
    zu1: Any = None
    fwdzu: Any = None
    zu2: Any = None
    v1dzu: Any = None
    zu3: Any = None
    bldzu: Any = None
    zu4: Any = None
    v2dzu: Any = None
    zu5: Any = None
    shdzu: Any = None
    zu6: Any = None
    v3dzu: Any = None
    zu7: Any = None
    fwmatv: Any = None
    v1matv: Any = None
    blmatv: Any = None
    v2matv: Any = None
    shmatv: Any = None
    v3matv: Any = None
    chmatv: Any = None
    chvol: Any = None
    fwvol: Any = None
    v1vol: Any = None
    blvol: Any = None
    v2vol: Any = None
    somtdr: Any = None
    sombdr: Any = None
    shvol: Any = None
    v3vol: Any = None
    chmatf: Any = None
    fwmatf: Any = None
    blmatf: Any = None
    v1matf: Any = None
    v2matf: Any = None
    shmatf: Any = None
    v3matf: Any = None
    expected_fwarea: Any = None
    expected_r1: Any = None
    expected_r2: Any = None
    expected_r3: Any = None
    expected_r4: Any = None
    expected_r5: Any = None
    expected_r6: Any = None
    expected_r7: Any = None
    expected_zl7: Any = None
    expected_zl6: Any = None
    expected_zl5: Any = None
    expected_zl4: Any = None
    expected_zl3: Any = None
    expected_zl2: Any = None
    expected_zl1: Any = None
    expected_zu1: Any = None
    expected_zu2: Any = None
    expected_zu3: Any = None
    expected_zu4: Any = None
    expected_zu5: Any = None
    expected_zu6: Any = None
    expected_zu7: Any = None
    expected_fwmatv: Any = None
    expected_v1matv: Any = None
    expected_blmatv: Any = None
    expected_v2matv: Any = None
    expected_shmatv: Any = None
    expected_v3matv: Any = None
    expected_chmatv: Any = None
    expected_chvol: Any = None
    expected_fwvol: Any = None
    expected_v1vol: Any = None
    expected_blvol: Any = None
    expected_v2vol: Any = None
    expected_shvol: Any = None
    expected_v3vol: Any = None


@pytest.mark.parametrize(
    "sombldparam",
    (
        SombldParam(
            fwarea=0,
            ifetyp=2,
            chrad=3.5,
            r1=0,
            fwdr=0.055,
            r2=0,
            v1dr=0.0050000000000000001,
            r3=0,
            bldr=0.55000000000000004,
            r4=0,
            v2dr=2.1899999999999999,
            r5=0,
            shdr=0.20000000000000001,
            r6=0,
            v3dr=3.5,
            r7=0,
            zl7=0,
            v3dzl=0,
            zl6=0,
            shdzl=0.35000000000000003,
            zl5=0,
            v2dzl=0,
            zl4=0,
            bldzl=0.65000000000000013,
            zl3=0,
            v1dzl=1.75,
            zl2=0,
            fwdzl=0,
            zl1=0,
            chdzl=3.1000000000000001,
            chdzu=3.7000000000000002,
            zu1=0,
            fwdzu=0.055,
            zu2=0,
            v1dzu=0.0050000000000000001,
            zu3=0,
            bldzu=0.55000000000000004,
            zu4=0,
            v2dzu=0.85000000000000009,
            zu5=0,
            shdzu=0.20000000000000001,
            zu6=0,
            v3dzu=13.640000000000001,
            zu7=0,
            fwmatv=numpy.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            v1matv=numpy.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            blmatv=numpy.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            v2matv=numpy.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            shmatv=numpy.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            v3matv=numpy.array(
                (
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            chmatv=numpy.array(
                numpy.array((0, 0, 0, 0, 0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            chvol=0,
            fwvol=numpy.array(numpy.array((0, 0, 0), order="F"), order="F").transpose(),
            v1vol=numpy.array(numpy.array((0, 0, 0), order="F"), order="F").transpose(),
            blvol=numpy.array(numpy.array((0, 0, 0), order="F"), order="F").transpose(),
            v2vol=numpy.array(numpy.array((0, 0, 0), order="F"), order="F").transpose(),
            somtdr=2.7000000000000002,
            sombdr=2.7000000000000002,
            shvol=numpy.array(numpy.array((0, 0, 0), order="F"), order="F").transpose(),
            v3vol=numpy.array(numpy.array((0, 0, 0), order="F"), order="F").transpose(),
            chmatf=numpy.array(
                numpy.array((1, 0, 0, 0, 0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            fwmatf=numpy.array(
                (
                    (0.050000000000000003, 0, 1),
                    (0, 0, 0),
                    (0.086400000000000005, 0.090899999999999995, 0),
                    (0.86360000000000015, 0.90910000000000002, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            blmatf=numpy.array(
                (
                    (0.050000000000000003, 0, 0),
                    (0, 0, 0),
                    (0.0086, 0.0091000000000000004, 0),
                    (0.94140000000000001, 0.9909, 1),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            v1matf=numpy.array(
                (
                    (0.050000000000000003, 0, 1),
                    (0, 0, 0),
                    (0.95000000000000007, 1, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            v2matf=numpy.array(
                (
                    (0.95000000000000007, 0.95000000000000007, 0.95000000000000007),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0.050000000000000003, 0.050000000000000003, 0.5),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            shmatf=numpy.array(
                (
                    (0, 0.30000000000000004, 0.30000000000000004),
                    (1, 0.70000000000000007, 0.70000000000000007),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            v3matf=numpy.array(
                (
                    (1, 1, 1),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            expected_fwarea=258.39044229575416,
            expected_r1=3.5,
            expected_r2=3.5550000000000002,
            expected_r3=3.5600000000000001,
            expected_r4=4.1100000000000003,
            expected_r5=6.3000000000000007,
            expected_r6=6.5000000000000009,
            expected_r7=10,
            expected_zl7=5.8499999999999996,
            expected_zl6=5.8499999999999996,
            expected_zl5=5.5,
            expected_zl4=5.5,
            expected_zl3=4.8499999999999996,
            expected_zl2=3.1000000000000001,
            expected_zl1=3.1000000000000001,
            expected_zu1=3.7000000000000002,
            expected_zu2=3.7550000000000003,
            expected_zu3=3.7600000000000002,
            expected_zu4=4.3100000000000005,
            expected_zu5=5.1600000000000001,
            expected_zu6=5.3600000000000003,
            expected_zu7=19,
            expected_fwmatv=numpy.array(
                (
                    (0, 0, 1.4221859043107128),
                    (-0, 0, 0),
                    (0, 0.1954426256072351, 0),
                    (0, 1.954641264461358, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                ),
                order="F",
            ).transpose(),
            expected_v1matv=numpy.array(
                (
                    (0, 0, 23.355974233572464),
                    (-0, 0, 0),
                    (0, 0.19879697242619909, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                ),
                order="F",
            ).transpose(),
            expected_blmatv=numpy.array(
                (
                    (0, 0, 0),
                    (-0, 0, 0),
                    (0, 0.60678614525284247, 0),
                    (0, 66.073010036378193, 93.868477240185996),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                ),
                order="F",
            ).transpose(),
            expected_v2matv=numpy.array(
                (
                    (0, 479.3074876785567, 543.51842785324482),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (0, 25.226709877818774, 286.06233044907617),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                    (-0, 0, 0),
                ),
                order="F",
            ).transpose(),
            expected_shmatv=numpy.array(
                (
                    (0, 7.4813887452587426, 13.092430304202773),
                    (90.156169335658547, 17.456573738937063, 30.549004043139806),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            expected_v3matv=numpy.array(
                (
                    (4508.4603472585413, 1810.4684303372624, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                    (0, 0, 0),
                ),
                order="F",
            ).transpose(),
            expected_chmatv=numpy.array(
                numpy.array((82.10028801381327, 0, 0, 0, 0, 0, 0, 0, 0), order="F"),
                order="F",
            ).transpose(),
            expected_chvol=82.10028801381327,
            expected_fwvol=numpy.array(
                numpy.array(
                    (-0.24380329788183547, 2.150083890068593, 1.4221859043107128),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_v1vol=numpy.array(
                numpy.array(
                    (-0.022352431730291104, 0.19879697242619909, 23.355974233572464),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_blvol=numpy.array(
                numpy.array(
                    (-2.6505617218337005, 66.679796181631033, 93.868477240185996),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_v2vol=numpy.array(
                numpy.array(
                    (-14.324343031454905, 504.53419755637543, 572.12466089815234),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_shvol=numpy.array(
                numpy.array(
                    (90.156169335658547, 24.937962484195804, 43.641434347342575),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_v3vol=numpy.array(
                numpy.array((4508.4603472585413, 1810.4684303372624, 0), order="F"),
                order="F",
            ).transpose(),
        ),
    ),
)
def test_sombld(sombldparam, monkeypatch, ife):
    """
    Automatically generated Regression Unit Test for sombld.

    This test was generated using data from IFE/new_IFE2.IN.DAT.

    :param sombldparam: the data used to mock and assert in this test.
    :type sombldparam: sombldparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(build_variables, "fwarea", sombldparam.fwarea)
    monkeypatch.setattr(ife_variables, "ifetyp", sombldparam.ifetyp)
    monkeypatch.setattr(ife_variables, "chrad", sombldparam.chrad)
    monkeypatch.setattr(ife_variables, "r1", sombldparam.r1)
    monkeypatch.setattr(ife_variables, "fwdr", sombldparam.fwdr)
    monkeypatch.setattr(ife_variables, "r2", sombldparam.r2)
    monkeypatch.setattr(ife_variables, "v1dr", sombldparam.v1dr)
    monkeypatch.setattr(ife_variables, "r3", sombldparam.r3)
    monkeypatch.setattr(ife_variables, "bldr", sombldparam.bldr)
    monkeypatch.setattr(ife_variables, "r4", sombldparam.r4)
    monkeypatch.setattr(ife_variables, "v2dr", sombldparam.v2dr)
    monkeypatch.setattr(ife_variables, "r5", sombldparam.r5)
    monkeypatch.setattr(ife_variables, "shdr", sombldparam.shdr)
    monkeypatch.setattr(ife_variables, "r6", sombldparam.r6)
    monkeypatch.setattr(ife_variables, "v3dr", sombldparam.v3dr)
    monkeypatch.setattr(ife_variables, "r7", sombldparam.r7)
    monkeypatch.setattr(ife_variables, "zl7", sombldparam.zl7)
    monkeypatch.setattr(ife_variables, "v3dzl", sombldparam.v3dzl)
    monkeypatch.setattr(ife_variables, "zl6", sombldparam.zl6)
    monkeypatch.setattr(ife_variables, "shdzl", sombldparam.shdzl)
    monkeypatch.setattr(ife_variables, "zl5", sombldparam.zl5)
    monkeypatch.setattr(ife_variables, "v2dzl", sombldparam.v2dzl)
    monkeypatch.setattr(ife_variables, "zl4", sombldparam.zl4)
    monkeypatch.setattr(ife_variables, "bldzl", sombldparam.bldzl)
    monkeypatch.setattr(ife_variables, "zl3", sombldparam.zl3)
    monkeypatch.setattr(ife_variables, "v1dzl", sombldparam.v1dzl)
    monkeypatch.setattr(ife_variables, "zl2", sombldparam.zl2)
    monkeypatch.setattr(ife_variables, "fwdzl", sombldparam.fwdzl)
    monkeypatch.setattr(ife_variables, "zl1", sombldparam.zl1)
    monkeypatch.setattr(ife_variables, "chdzl", sombldparam.chdzl)
    monkeypatch.setattr(ife_variables, "chdzu", sombldparam.chdzu)
    monkeypatch.setattr(ife_variables, "zu1", sombldparam.zu1)
    monkeypatch.setattr(ife_variables, "fwdzu", sombldparam.fwdzu)
    monkeypatch.setattr(ife_variables, "zu2", sombldparam.zu2)
    monkeypatch.setattr(ife_variables, "v1dzu", sombldparam.v1dzu)
    monkeypatch.setattr(ife_variables, "zu3", sombldparam.zu3)
    monkeypatch.setattr(ife_variables, "bldzu", sombldparam.bldzu)
    monkeypatch.setattr(ife_variables, "zu4", sombldparam.zu4)
    monkeypatch.setattr(ife_variables, "v2dzu", sombldparam.v2dzu)
    monkeypatch.setattr(ife_variables, "zu5", sombldparam.zu5)
    monkeypatch.setattr(ife_variables, "shdzu", sombldparam.shdzu)
    monkeypatch.setattr(ife_variables, "zu6", sombldparam.zu6)
    monkeypatch.setattr(ife_variables, "v3dzu", sombldparam.v3dzu)
    monkeypatch.setattr(ife_variables, "zu7", sombldparam.zu7)
    monkeypatch.setattr(ife_variables, "fwmatv", sombldparam.fwmatv)
    monkeypatch.setattr(ife_variables, "v1matv", sombldparam.v1matv)
    monkeypatch.setattr(ife_variables, "blmatv", sombldparam.blmatv)
    monkeypatch.setattr(ife_variables, "v2matv", sombldparam.v2matv)
    monkeypatch.setattr(ife_variables, "shmatv", sombldparam.shmatv)
    monkeypatch.setattr(ife_variables, "v3matv", sombldparam.v3matv)
    monkeypatch.setattr(ife_variables, "chmatv", sombldparam.chmatv)
    monkeypatch.setattr(ife_variables, "chvol", sombldparam.chvol)
    monkeypatch.setattr(ife_variables, "fwvol", sombldparam.fwvol)
    monkeypatch.setattr(ife_variables, "v1vol", sombldparam.v1vol)
    monkeypatch.setattr(ife_variables, "blvol", sombldparam.blvol)
    monkeypatch.setattr(ife_variables, "v2vol", sombldparam.v2vol)
    monkeypatch.setattr(ife_variables, "somtdr", sombldparam.somtdr)
    monkeypatch.setattr(ife_variables, "sombdr", sombldparam.sombdr)
    monkeypatch.setattr(ife_variables, "shvol", sombldparam.shvol)
    monkeypatch.setattr(ife_variables, "v3vol", sombldparam.v3vol)
    monkeypatch.setattr(ife_variables, "chmatf", sombldparam.chmatf)
    monkeypatch.setattr(ife_variables, "fwmatf", sombldparam.fwmatf)
    monkeypatch.setattr(ife_variables, "blmatf", sombldparam.blmatf)
    monkeypatch.setattr(ife_variables, "v1matf", sombldparam.v1matf)
    monkeypatch.setattr(ife_variables, "v2matf", sombldparam.v2matf)
    monkeypatch.setattr(ife_variables, "shmatf", sombldparam.shmatf)
    monkeypatch.setattr(ife_variables, "v3matf", sombldparam.v3matf)

    ife.sombld()

    assert build_variables.fwarea == pytest.approx(sombldparam.expected_fwarea)
    assert ife_variables.r1 == pytest.approx(sombldparam.expected_r1)
    assert ife_variables.r2 == pytest.approx(sombldparam.expected_r2)
    assert ife_variables.r3 == pytest.approx(sombldparam.expected_r3)
    assert ife_variables.r4 == pytest.approx(sombldparam.expected_r4)
    assert ife_variables.r5 == pytest.approx(sombldparam.expected_r5)
    assert ife_variables.r6 == pytest.approx(sombldparam.expected_r6)
    assert ife_variables.r7 == pytest.approx(sombldparam.expected_r7)
    assert ife_variables.zl7 == pytest.approx(sombldparam.expected_zl7)
    assert ife_variables.zl6 == pytest.approx(sombldparam.expected_zl6)
    assert ife_variables.zl5 == pytest.approx(sombldparam.expected_zl5)
    assert ife_variables.zl4 == pytest.approx(sombldparam.expected_zl4)
    assert ife_variables.zl3 == pytest.approx(sombldparam.expected_zl3)
    assert ife_variables.zl2 == pytest.approx(sombldparam.expected_zl2)
    assert ife_variables.zl1 == pytest.approx(sombldparam.expected_zl1)
    assert ife_variables.zu1 == pytest.approx(sombldparam.expected_zu1)
    assert ife_variables.zu2 == pytest.approx(sombldparam.expected_zu2)
    assert ife_variables.zu3 == pytest.approx(sombldparam.expected_zu3)
    assert ife_variables.zu4 == pytest.approx(sombldparam.expected_zu4)
    assert ife_variables.zu5 == pytest.approx(sombldparam.expected_zu5)
    assert ife_variables.zu6 == pytest.approx(sombldparam.expected_zu6)
    assert ife_variables.zu7 == pytest.approx(sombldparam.expected_zu7)
    assert ife_variables.fwmatv == pytest.approx(sombldparam.expected_fwmatv)
    assert ife_variables.v1matv == pytest.approx(sombldparam.expected_v1matv)
    assert ife_variables.blmatv == pytest.approx(sombldparam.expected_blmatv)
    assert ife_variables.v2matv == pytest.approx(sombldparam.expected_v2matv)
    assert ife_variables.shmatv == pytest.approx(sombldparam.expected_shmatv)
    assert ife_variables.v3matv == pytest.approx(sombldparam.expected_v3matv)
    assert ife_variables.chmatv == pytest.approx(sombldparam.expected_chmatv)
    assert ife_variables.chvol == pytest.approx(sombldparam.expected_chvol)
    assert ife_variables.fwvol == pytest.approx(sombldparam.expected_fwvol)
    assert ife_variables.v1vol == pytest.approx(sombldparam.expected_v1vol)
    assert ife_variables.blvol == pytest.approx(sombldparam.expected_blvol)
    assert ife_variables.v2vol == pytest.approx(sombldparam.expected_v2vol)
    assert ife_variables.shvol == pytest.approx(sombldparam.expected_shvol)
    assert ife_variables.v3vol == pytest.approx(sombldparam.expected_v3vol)
