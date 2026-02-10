"""Unit tests for ife."""

from typing import Any, NamedTuple

import numpy as np
import pytest

from process.availability import Availability
from process.costs import Costs
from process.data_structure import (
    buildings_variables,
    cost_variables,
    first_wall_variables,
    fwbs_variables,
    heat_transport_variables,
    ife_variables,
    physics_variables,
)
from process.ife import IFE


@pytest.fixture
def ife():
    """Provides IFE object for testing.

    :returns: initialised IFE object
    :rtype: process.ife.IFE
    """
    return IFE(Availability(), Costs())


def test_ifetgt(monkeypatch, ife):
    """Test ifetgt.

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    # Mock module variables
    # Repetition Rate (Hz)
    monkeypatch.setattr(ife_variables, "reprat", 4.0)
    # IFE target factory power at 6 Hz repetition rate
    monkeypatch.setattr(ife_variables, "ptargf", 2.0)
    monkeypatch.setattr(ife_variables, "tfacmw", 0.0)

    ife.ifetgt()
    assert ife_variables.tfacmw == pytest.approx(1.506, abs=0.001)


class SombldParam(NamedTuple):
    a_fw_total: Any = None
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
    expected_a_fw_total: Any = None
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
            a_fw_total=0.0,
            ifetyp=2,
            chrad=3.5,
            r1=0.0,
            fwdr=0.055,
            r2=0.0,
            v1dr=0.0050000000000000001,
            r3=0.0,
            bldr=0.55000000000000004,
            r4=0.0,
            v2dr=2.1899999999999999,
            r5=0.0,
            shdr=0.20000000000000001,
            r6=0.0,
            v3dr=3.5,
            r7=0.0,
            zl7=0.0,
            v3dzl=0.0,
            zl6=0.0,
            shdzl=0.35000000000000003,
            zl5=0.0,
            v2dzl=0.0,
            zl4=0.0,
            bldzl=0.65000000000000013,
            zl3=0.0,
            v1dzl=1.75,
            zl2=0.0,
            fwdzl=0.0,
            zl1=0.0,
            chdzl=3.1000000000000001,
            chdzu=3.7000000000000002,
            zu1=0.0,
            fwdzu=0.055,
            zu2=0.0,
            v1dzu=0.0050000000000000001,
            zu3=0.0,
            bldzu=0.55000000000000004,
            zu4=0.0,
            v2dzu=0.85000000000000009,
            zu5=0.0,
            shdzu=0.20000000000000001,
            zu6=0.0,
            v3dzu=13.640000000000001,
            zu7=0.0,
            fwmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v1matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            blmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v2matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            shmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v3matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            chmatv=np.array(
                np.array(
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                ),
            ).T,
            chvol=0.0,
            fwvol=np.array(
                np.array(
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v1vol=np.array(
                np.array(
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            blvol=np.array(
                np.array(
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v2vol=np.array(
                np.array(
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            somtdr=2.7000000000000002,
            sombdr=2.7000000000000002,
            shvol=np.array(
                np.array(
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v3vol=np.array(
                np.array(
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            chmatf=np.array(
                np.array(
                    (1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                ),
            ).T,
            fwmatf=np.array(
                (
                    (0.050000000000000003, 0.0, 1.0),
                    (0.0, 0.0, 0.0),
                    (0.086400000000000005, 0.090899999999999995, 0.0),
                    (0.86360000000000015, 0.90910000000000002, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            blmatf=np.array(
                (
                    (0.050000000000000003, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0086, 0.0091000000000000004, 0.0),
                    (0.94140000000000001, 0.9909, 1.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v1matf=np.array(
                (
                    (0.050000000000000003, 0.0, 1.0),
                    (0.0, 0.0, 0.0),
                    (0.95000000000000007, 1.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v2matf=np.array(
                (
                    (0.95000000000000007, 0.95000000000000007, 0.95000000000000007),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.050000000000000003, 0.050000000000000003, 0.5),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            shmatf=np.array(
                (
                    (0.0, 0.30000000000000004, 0.30000000000000004),
                    (1.0, 0.70000000000000007, 0.70000000000000007),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v3matf=np.array(
                (
                    (1.0, 1.0, 1.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_a_fw_total=258.39044229575416,
            expected_r1=3.5,
            expected_r2=3.5550000000000002,
            expected_r3=3.5600000000000001,
            expected_r4=4.1100000000000003,
            expected_r5=6.3000000000000007,
            expected_r6=6.5000000000000009,
            expected_r7=10.0,
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
            expected_zu7=19.0,
            expected_fwmatv=np.array(
                (
                    (0.0, 0.0, 1.4221859043107128),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.1954426256072351, 0.0),
                    (0.0, 1.954641264461358, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_v1matv=np.array(
                (
                    (0.0, 0.0, 23.355974233572464),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.19879697242619909, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_blmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.60678614525284247, 0.0),
                    (0.0, 66.073010036378193, 93.868477240185996),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_v2matv=np.array(
                (
                    (0.0, 479.3074876785567, 543.51842785324482),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 25.226709877818774, 286.06233044907617),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_shmatv=np.array(
                (
                    (0.0, 7.4813887452587426, 13.092430304202773),
                    (90.156169335658547, 17.456573738937063, 30.549004043139806),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_v3matv=np.array(
                (
                    (4508.4603472585413, 1810.4684303372624, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_chmatv=np.array(
                np.array(
                    (82.10028801381327, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                ),
            ).T,
            expected_chvol=82.10028801381327,
            expected_fwvol=np.array(
                np.array(
                    (-0.24380329788183547, 2.150083890068593, 1.4221859043107128),
                ),
            ).T,
            expected_v1vol=np.array(
                np.array(
                    (-0.022352431730291104, 0.19879697242619909, 23.355974233572464),
                ),
            ).T,
            expected_blvol=np.array(
                np.array(
                    (-2.6505617218337005, 66.679796181631033, 93.868477240185996),
                ),
            ).T,
            expected_v2vol=np.array(
                np.array(
                    (-14.324343031454905, 504.53419755637543, 572.12466089815234),
                ),
            ).T,
            expected_shvol=np.array(
                np.array(
                    (90.156169335658547, 24.937962484195804, 43.641434347342575),
                ),
            ).T,
            expected_v3vol=np.array(
                np.array(
                    (4508.4603472585413, 1810.4684303372624, 0.0),
                ),
            ).T,
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
    monkeypatch.setattr(first_wall_variables, "a_fw_total", sombldparam.a_fw_total)
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

    assert first_wall_variables.a_fw_total == pytest.approx(
        sombldparam.expected_a_fw_total
    )
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


class DriverParam(NamedTuple):
    edrive: Any = None
    etave: Any = None
    gainve: Any = None
    expected_etadrv: Any = None
    expected_gain: Any = None


@pytest.mark.parametrize(
    "driverparam",
    (
        DriverParam(
            edrive=5000000.0,
            etave=np.array(
                (
                    0.082000000000000003,
                    0.079000000000000001,
                    0.075999999999999998,
                    0.072999999999999995,
                    0.069000000000000006,
                    0.066000000000000003,
                    0.062,
                    0.058999999999999997,
                    0.055,
                    0.050999999999999997,
                ),
            ).T,
            gainve=np.array(
                (60.0, 95.0, 115.0, 125.0, 133.0, 141.0, 152.0, 160.0, 165.0, 170.0),
            ).T,
            expected_etadrv=0.069000000000000006,
            expected_gain=133.0,
        ),
    ),
)
def test_driver(driverparam, ife):
    """
    Automatically generated Regression Unit Test for driver.

    This test was generated using data from ife/normal_driver.IN.DAT.

    :param driverparam: the data used to mock and assert in this test.
    :type driverparam: driverparam
    """

    gain, etadrv = ife.driver(
        edrive=driverparam.edrive, etave=driverparam.etave, gainve=driverparam.gainve
    )

    assert etadrv == pytest.approx(driverparam.expected_etadrv)
    assert gain == pytest.approx(driverparam.expected_gain)


class LasdrvParam(NamedTuple):
    edrive: Any = None
    expected_etadrv: Any = None
    expected_gain: Any = None


@pytest.mark.parametrize(
    "lasdrvparam",
    (
        LasdrvParam(
            edrive=5000000.0,
            expected_etadrv=0.069000000000000006,
            expected_gain=136.0,
        ),
    ),
)
def test_lasdrv(lasdrvparam, ife):
    """
    Automatically generated Regression Unit Test for lasdrv.

    This test was generated using data from ife/laser_drive.IN.DAT.

    :param lasdrvparam: the data used to mock and assert in this test.
    :type lasdrvparam: lasdrvparam
    """

    gain, etadrv = ife.lasdrv(edrive=lasdrvparam.edrive)

    assert etadrv == pytest.approx(lasdrvparam.expected_etadrv)
    assert gain == pytest.approx(lasdrvparam.expected_gain)


class HylbldParam(NamedTuple):
    a_fw_total: Any = None
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
    expected_a_fw_total: Any = None
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


@pytest.mark.parametrize(
    "hylbldparam",
    (
        HylbldParam(
            a_fw_total=0.0,
            ifetyp=3,
            chrad=3.5,
            r1=0.0,
            fwdr=0.055,
            r2=0.0,
            v1dr=0.0050000000000000001,
            r3=0.0,
            bldr=0.55000000000000004,
            r4=0.0,
            v2dr=2.1899999999999999,
            r5=0.0,
            shdr=0.20000000000000001,
            r6=0.0,
            v3dr=3.5,
            r7=0.0,
            zl7=0.0,
            v3dzl=0.0,
            zl6=0.0,
            shdzl=0.35000000000000003,
            zl5=0.0,
            v2dzl=0.0,
            zl4=0.0,
            bldzl=0.65000000000000013,
            zl3=0.0,
            v1dzl=1.75,
            zl2=0.0,
            fwdzl=0.0,
            zl1=0.0,
            chdzl=3.1000000000000001,
            chdzu=3.7000000000000002,
            zu1=0.0,
            fwdzu=0.055,
            zu2=0.0,
            v1dzu=0.0050000000000000001,
            zu3=0.0,
            bldzu=0.55000000000000004,
            zu4=0.0,
            v2dzu=0.85000000000000009,
            zu5=0.0,
            shdzu=0.20000000000000001,
            zu6=0.0,
            v3dzu=13.640000000000001,
            zu7=0.0,
            fwmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v1matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            blmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v2matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            shmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v3matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            chmatv=np.array(
                np.array(
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                ),
            ).T,
            expected_a_fw_total=281.91872215483801,
            expected_r1=3.5,
            expected_r2=3.5550000000000002,
            expected_r3=3.5600000000000001,
            expected_r4=4.1100000000000003,
            expected_r5=6.3000000000000007,
            expected_r6=6.5000000000000009,
            expected_r7=10.0,
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
            expected_zu7=19.0,
            expected_fwmatv=np.array([
                [
                    0.564099880474099,
                    0.564099880474099,
                    0.564099880474099,
                    0.0,
                    0.0,
                    0.0,
                    10.717897729007879,
                    10.717897729007879,
                    10.717897729007879,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_v1matv=np.array([
                [
                    0.9617133751957783,
                    0.9617133751957783,
                    0.9617133751957783,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_blmatv=np.array([
                [
                    5.70202090409477,
                    5.70202090409477,
                    5.70202090409477,
                    0.0,
                    0.0,
                    0.0,
                    51.318188136852925,
                    51.318188136852925,
                    51.318188136852925,
                ],
                [
                    0.0,
                    0.0,
                    0.0,
                    5.627242752258617,
                    5.627242752258617,
                    5.627242752258617,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    10.348278966422315,
                    10.348278966422315,
                    10.348278966422315,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ]),
            expected_v2matv=np.array([
                [
                    702.6090256928657,
                    702.6090256928657,
                    702.6090256928657,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_shmatv=np.array([
                [
                    4.286640343970205,
                    4.286640343970205,
                    4.286640343970205,
                    16.289233307086782,
                    16.289233307086782,
                    16.289233307086782,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    17.653394518684465,
                    17.653394518684465,
                    17.653394518684465,
                ],
                [
                    4.413348629671108,
                    4.413348629671108,
                    4.413348629671108,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ]),
            expected_v3matv=np.array([
                [
                    2033.7963980993259,
                    2033.7963980993259,
                    2033.7963980993259,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_chmatv=np.array([
                323.26988405438965,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ]),
        ),
    ),
)
def test_hylbld(hylbldparam, monkeypatch, ife):
    """
    Automatically generated Regression Unit Test for hylbld.

    This test was generated using data from ife/hylbld.IN.DAT.

    :param hylbldparam: the data used to mock and assert in this test.
    :type hylbldparam: hylbldparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(first_wall_variables, "a_fw_total", hylbldparam.a_fw_total)
    monkeypatch.setattr(ife_variables, "ifetyp", hylbldparam.ifetyp)
    monkeypatch.setattr(ife_variables, "chrad", hylbldparam.chrad)
    monkeypatch.setattr(ife_variables, "r1", hylbldparam.r1)
    monkeypatch.setattr(ife_variables, "fwdr", hylbldparam.fwdr)
    monkeypatch.setattr(ife_variables, "r2", hylbldparam.r2)
    monkeypatch.setattr(ife_variables, "v1dr", hylbldparam.v1dr)
    monkeypatch.setattr(ife_variables, "r3", hylbldparam.r3)
    monkeypatch.setattr(ife_variables, "bldr", hylbldparam.bldr)
    monkeypatch.setattr(ife_variables, "r4", hylbldparam.r4)
    monkeypatch.setattr(ife_variables, "v2dr", hylbldparam.v2dr)
    monkeypatch.setattr(ife_variables, "r5", hylbldparam.r5)
    monkeypatch.setattr(ife_variables, "shdr", hylbldparam.shdr)
    monkeypatch.setattr(ife_variables, "r6", hylbldparam.r6)
    monkeypatch.setattr(ife_variables, "v3dr", hylbldparam.v3dr)
    monkeypatch.setattr(ife_variables, "r7", hylbldparam.r7)
    monkeypatch.setattr(ife_variables, "zl7", hylbldparam.zl7)
    monkeypatch.setattr(ife_variables, "v3dzl", hylbldparam.v3dzl)
    monkeypatch.setattr(ife_variables, "zl6", hylbldparam.zl6)
    monkeypatch.setattr(ife_variables, "shdzl", hylbldparam.shdzl)
    monkeypatch.setattr(ife_variables, "zl5", hylbldparam.zl5)
    monkeypatch.setattr(ife_variables, "v2dzl", hylbldparam.v2dzl)
    monkeypatch.setattr(ife_variables, "zl4", hylbldparam.zl4)
    monkeypatch.setattr(ife_variables, "bldzl", hylbldparam.bldzl)
    monkeypatch.setattr(ife_variables, "zl3", hylbldparam.zl3)
    monkeypatch.setattr(ife_variables, "v1dzl", hylbldparam.v1dzl)
    monkeypatch.setattr(ife_variables, "zl2", hylbldparam.zl2)
    monkeypatch.setattr(ife_variables, "fwdzl", hylbldparam.fwdzl)
    monkeypatch.setattr(ife_variables, "zl1", hylbldparam.zl1)
    monkeypatch.setattr(ife_variables, "chdzl", hylbldparam.chdzl)
    monkeypatch.setattr(ife_variables, "chdzu", hylbldparam.chdzu)
    monkeypatch.setattr(ife_variables, "zu1", hylbldparam.zu1)
    monkeypatch.setattr(ife_variables, "fwdzu", hylbldparam.fwdzu)
    monkeypatch.setattr(ife_variables, "zu2", hylbldparam.zu2)
    monkeypatch.setattr(ife_variables, "v1dzu", hylbldparam.v1dzu)
    monkeypatch.setattr(ife_variables, "zu3", hylbldparam.zu3)
    monkeypatch.setattr(ife_variables, "bldzu", hylbldparam.bldzu)
    monkeypatch.setattr(ife_variables, "zu4", hylbldparam.zu4)
    monkeypatch.setattr(ife_variables, "v2dzu", hylbldparam.v2dzu)
    monkeypatch.setattr(ife_variables, "zu5", hylbldparam.zu5)
    monkeypatch.setattr(ife_variables, "shdzu", hylbldparam.shdzu)
    monkeypatch.setattr(ife_variables, "zu6", hylbldparam.zu6)
    monkeypatch.setattr(ife_variables, "v3dzu", hylbldparam.v3dzu)
    monkeypatch.setattr(ife_variables, "zu7", hylbldparam.zu7)
    monkeypatch.setattr(ife_variables, "fwmatv", hylbldparam.fwmatv)
    monkeypatch.setattr(ife_variables, "v1matv", hylbldparam.v1matv)
    monkeypatch.setattr(ife_variables, "blmatv", hylbldparam.blmatv)
    monkeypatch.setattr(ife_variables, "v2matv", hylbldparam.v2matv)
    monkeypatch.setattr(ife_variables, "shmatv", hylbldparam.shmatv)
    monkeypatch.setattr(ife_variables, "v3matv", hylbldparam.v3matv)
    monkeypatch.setattr(ife_variables, "chmatv", hylbldparam.chmatv)

    ife.hylbld()

    assert first_wall_variables.a_fw_total == pytest.approx(
        hylbldparam.expected_a_fw_total
    )
    assert ife_variables.r1 == pytest.approx(hylbldparam.expected_r1)
    assert ife_variables.r2 == pytest.approx(hylbldparam.expected_r2)
    assert ife_variables.r3 == pytest.approx(hylbldparam.expected_r3)
    assert ife_variables.r4 == pytest.approx(hylbldparam.expected_r4)
    assert ife_variables.r5 == pytest.approx(hylbldparam.expected_r5)
    assert ife_variables.r6 == pytest.approx(hylbldparam.expected_r6)
    assert ife_variables.r7 == pytest.approx(hylbldparam.expected_r7)
    assert ife_variables.zl7 == pytest.approx(hylbldparam.expected_zl7)
    assert ife_variables.zl6 == pytest.approx(hylbldparam.expected_zl6)
    assert ife_variables.zl5 == pytest.approx(hylbldparam.expected_zl5)
    assert ife_variables.zl4 == pytest.approx(hylbldparam.expected_zl4)
    assert ife_variables.zl3 == pytest.approx(hylbldparam.expected_zl3)
    assert ife_variables.zl2 == pytest.approx(hylbldparam.expected_zl2)
    assert ife_variables.zl1 == pytest.approx(hylbldparam.expected_zl1)
    assert ife_variables.zu1 == pytest.approx(hylbldparam.expected_zu1)
    assert ife_variables.zu2 == pytest.approx(hylbldparam.expected_zu2)
    assert ife_variables.zu3 == pytest.approx(hylbldparam.expected_zu3)
    assert ife_variables.zu4 == pytest.approx(hylbldparam.expected_zu4)
    assert ife_variables.zu5 == pytest.approx(hylbldparam.expected_zu5)
    assert ife_variables.zu6 == pytest.approx(hylbldparam.expected_zu6)
    assert ife_variables.zu7 == pytest.approx(hylbldparam.expected_zu7)
    assert ife_variables.fwmatv == pytest.approx(hylbldparam.expected_fwmatv)
    assert ife_variables.v1matv == pytest.approx(hylbldparam.expected_v1matv)
    assert ife_variables.blmatv == pytest.approx(hylbldparam.expected_blmatv)
    assert ife_variables.v2matv == pytest.approx(hylbldparam.expected_v2matv)
    assert ife_variables.shmatv == pytest.approx(hylbldparam.expected_shmatv)
    assert ife_variables.v3matv == pytest.approx(hylbldparam.expected_v3matv)
    assert ife_variables.chmatv == pytest.approx(hylbldparam.expected_chmatv)


class IondrvParam(NamedTuple):
    edrive: Any = None
    expected_gain: Any = None
    expected_etadrv: Any = None


@pytest.mark.parametrize(
    "iondrvparam",
    (
        IondrvParam(
            edrive=5000000.0,
            expected_gain=87.0,
            expected_etadrv=0.28199999999999997,
        ),
    ),
)
def test_iondrv(iondrvparam, ife):
    """
    Automatically generated Regression Unit Test for iondrv.

    This test was generated using data from ife/IFE.IN.DAT.

    :param iondrvparam: the data used to mock and assert in this test.
    :type iondrvparam: iondrvparam
    """

    gain, etadrv = ife.iondrv(edrive=iondrvparam.edrive)

    assert gain == pytest.approx(iondrvparam.expected_gain)
    assert etadrv == pytest.approx(iondrvparam.expected_etadrv)


class IfefbsParam(NamedTuple):
    a_fw_total: Any = None
    life_plant: Any = None
    abktflnc: Any = None
    f_t_plant_available: Any = None
    den_steel: Any = None
    m_fw_total: Any = None
    m_blkt_total: Any = None
    whtshld: Any = None
    m_blkt_beryllium: Any = None
    m_blkt_vanadium: Any = None
    m_blkt_steel_total: Any = None
    m_blkt_li2o: Any = None
    m_blkt_lithium: Any = None
    life_blkt_fpy: Any = None
    life_fw_fpy: Any = None
    chmatm: Any = None
    chmatv: Any = None
    fwmatm: Any = None
    fwmatv: Any = None
    v1matm: Any = None
    v1matv: Any = None
    blmatm: Any = None
    blmatv: Any = None
    v2matm: Any = None
    v2matv: Any = None
    shmatm: Any = None
    shmatv: Any = None
    v3matm: Any = None
    v3matv: Any = None
    mflibe: Any = None
    fbreed: Any = None
    ifetyp: Any = None
    pflux_fw_neutron_mw: Any = None
    expected_m_fw_total: Any = None
    expected_m_blkt_total: Any = None
    expected_whtshld: Any = None
    expected_life_blkt_fpy: Any = None
    expected_life_fw_fpy: Any = None
    expected_fwmatm: Any = None
    expected_v1matm: Any = None
    expected_blmatm: Any = None
    expected_v2matm: Any = None
    expected_shmatm: Any = None
    expected_mflibe: Any = None


@pytest.mark.parametrize(
    "ifefbsparam",
    (
        IfefbsParam(
            a_fw_total=188.02432031734912,
            life_plant=30.0,
            abktflnc=20.0,
            f_t_plant_available=0.75000000000000011,
            den_steel=7800.0,
            m_fw_total=0.0,
            m_blkt_total=0.0,
            whtshld=0.0,
            m_blkt_beryllium=0.0,
            m_blkt_vanadium=0.0,
            m_blkt_steel_total=0.0,
            m_blkt_li2o=0.0,
            m_blkt_lithium=0.0,
            life_blkt_fpy=0.0,
            life_fw_fpy=0.0,
            chmatm=np.array(
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            ).T,
            chmatv=np.array(
                (261.69466804402975, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            ).T,
            fwmatm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            fwmatv=np.array(
                (
                    (0.41446560639912183, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.71619656785768249, 0.19849778071617336, 0.0),
                    (7.1586499537256332, 1.9851961765574611, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v1matm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v1matv=np.array(
                (
                    (0.038306479877786521, 0.0, 69.676755145437298),
                    (0.0, 0.0, 0.0),
                    (0.72782311767794394, 0.19907644327267379, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            blmatm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            blmatv=np.array(
                (
                    (5.7053341062470615, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.98131746627449445, 0.26560582680483946, 0.0),
                    (107.42003055241966, 28.921847668232459, 34.494263221407721),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v2matm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v2matv=np.array(
                (
                    (667.47857440822247, 100.68702352994043, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (35.130451284643286, 5.2993170278916013, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            shmatm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            shmatv=np.array(
                (
                    (0.0, 7.9639373768501356, 13.936890409487711),
                    (85.732806879404109, 18.582520545983648, 32.519410955471322),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v3matm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v3matv=np.array(
                (
                    (2033.7963980993259, 4285.1323794964783, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            mflibe=0.0,
            fbreed=0.52600000000000002,
            ifetyp=1,
            pflux_fw_neutron_mw=8.8876851857005388,
            expected_m_fw_total=20574.366184891722,
            expected_m_blkt_total=347956.92928704334,
            expected_whtshld=1067310.9593707009,
            expected_life_blkt_fpy=3.000406304846492,
            expected_life_fw_fpy=3.000406304846492,
            expected_fwmatm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (1647.2521060726697, 456.5448956471987, 0.0),
                    (14460.47290652578, 4010.0962766460711, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_v1matm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (1673.9931706592711, 457.8758195271497, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_blmatm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (2257.0301724313372, 610.89340165113072, 0.0),
                    (216988.46171588771, 58422.132289829569, 69678.411707243591),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_v2matm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (70963.511594979442, 10704.620396341035, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_shmatm=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (668715.89365935209, 144943.66025867246, 253651.40545267632),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            expected_mflibe=939298.9596781712,
        ),
    ),
)
def test_ifefbs(ifefbsparam, monkeypatch, ife):
    """
    Automatically generated Regression Unit Test for ifefbs.

    This test was generated using data from ife/IFE.IN.DAT.

    :param ifefbsparam: the data used to mock and assert in this test.
    :type ifefbsparam: ifefbsparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(first_wall_variables, "a_fw_total", ifefbsparam.a_fw_total)
    monkeypatch.setattr(cost_variables, "life_plant", ifefbsparam.life_plant)
    monkeypatch.setattr(cost_variables, "abktflnc", ifefbsparam.abktflnc)
    monkeypatch.setattr(
        cost_variables, "f_t_plant_available", ifefbsparam.f_t_plant_available
    )
    monkeypatch.setattr(fwbs_variables, "den_steel", ifefbsparam.den_steel)
    monkeypatch.setattr(fwbs_variables, "m_fw_total", ifefbsparam.m_fw_total)
    monkeypatch.setattr(fwbs_variables, "m_blkt_total", ifefbsparam.m_blkt_total)
    monkeypatch.setattr(fwbs_variables, "whtshld", ifefbsparam.whtshld)
    monkeypatch.setattr(fwbs_variables, "m_blkt_beryllium", ifefbsparam.m_blkt_beryllium)
    monkeypatch.setattr(fwbs_variables, "m_blkt_vanadium", ifefbsparam.m_blkt_vanadium)
    monkeypatch.setattr(
        fwbs_variables, "m_blkt_steel_total", ifefbsparam.m_blkt_steel_total
    )
    monkeypatch.setattr(fwbs_variables, "m_blkt_li2o", ifefbsparam.m_blkt_li2o)
    monkeypatch.setattr(fwbs_variables, "m_blkt_lithium", ifefbsparam.m_blkt_lithium)
    monkeypatch.setattr(fwbs_variables, "life_blkt_fpy", ifefbsparam.life_blkt_fpy)
    monkeypatch.setattr(fwbs_variables, "life_fw_fpy", ifefbsparam.life_fw_fpy)
    monkeypatch.setattr(ife_variables, "chmatm", ifefbsparam.chmatm)
    monkeypatch.setattr(ife_variables, "chmatv", ifefbsparam.chmatv)
    monkeypatch.setattr(ife_variables, "fwmatm", ifefbsparam.fwmatm)
    monkeypatch.setattr(ife_variables, "fwmatv", ifefbsparam.fwmatv)
    monkeypatch.setattr(ife_variables, "v1matm", ifefbsparam.v1matm)
    monkeypatch.setattr(ife_variables, "v1matv", ifefbsparam.v1matv)
    monkeypatch.setattr(ife_variables, "blmatm", ifefbsparam.blmatm)
    monkeypatch.setattr(ife_variables, "blmatv", ifefbsparam.blmatv)
    monkeypatch.setattr(ife_variables, "v2matm", ifefbsparam.v2matm)
    monkeypatch.setattr(ife_variables, "v2matv", ifefbsparam.v2matv)
    monkeypatch.setattr(ife_variables, "shmatm", ifefbsparam.shmatm)
    monkeypatch.setattr(ife_variables, "shmatv", ifefbsparam.shmatv)
    monkeypatch.setattr(ife_variables, "v3matm", ifefbsparam.v3matm)
    monkeypatch.setattr(ife_variables, "v3matv", ifefbsparam.v3matv)
    monkeypatch.setattr(ife_variables, "mflibe", ifefbsparam.mflibe)
    monkeypatch.setattr(ife_variables, "fbreed", ifefbsparam.fbreed)
    monkeypatch.setattr(ife_variables, "ifetyp", ifefbsparam.ifetyp)
    monkeypatch.setattr(
        physics_variables, "pflux_fw_neutron_mw", ifefbsparam.pflux_fw_neutron_mw
    )

    ife.ifefbs(output=False)

    assert fwbs_variables.m_fw_total == pytest.approx(ifefbsparam.expected_m_fw_total)
    assert fwbs_variables.m_blkt_total == pytest.approx(
        ifefbsparam.expected_m_blkt_total
    )
    assert fwbs_variables.whtshld == pytest.approx(ifefbsparam.expected_whtshld)
    assert fwbs_variables.life_blkt_fpy == pytest.approx(
        ifefbsparam.expected_life_blkt_fpy
    )
    assert fwbs_variables.life_fw_fpy == pytest.approx(ifefbsparam.expected_life_fw_fpy)
    assert ife_variables.fwmatm == pytest.approx(ifefbsparam.expected_fwmatm)
    assert ife_variables.v1matm == pytest.approx(ifefbsparam.expected_v1matm)
    assert ife_variables.blmatm == pytest.approx(ifefbsparam.expected_blmatm)
    assert ife_variables.v2matm == pytest.approx(ifefbsparam.expected_v2matm)
    assert ife_variables.shmatm == pytest.approx(ifefbsparam.expected_shmatm)
    assert ife_variables.mflibe == pytest.approx(ifefbsparam.expected_mflibe)


class GenbldParam(NamedTuple):
    a_fw_total: Any = None
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
    expected_a_fw_total: Any = None
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


@pytest.mark.parametrize(
    "genbldparam",
    (
        GenbldParam(
            a_fw_total=0.0,
            ifetyp=0.0,
            chrad=3.5,
            r1=0.0,
            fwdr=0.055,
            r2=0.0,
            v1dr=0.0050000000000000001,
            r3=0.0,
            bldr=0.55000000000000004,
            r4=0.0,
            v2dr=2.1899999999999999,
            r5=0.0,
            shdr=0.20000000000000001,
            r6=0.0,
            v3dr=3.5,
            r7=0.0,
            zl7=0.0,
            v3dzl=0.0,
            zl6=0.0,
            shdzl=0.35000000000000003,
            zl5=0.0,
            v2dzl=0.0,
            zl4=0.0,
            bldzl=0.65000000000000013,
            zl3=0.0,
            v1dzl=1.75,
            zl2=0.0,
            fwdzl=0.0,
            zl1=0.0,
            chdzl=3.1000000000000001,
            chdzu=3.7000000000000002,
            zu1=0.0,
            fwdzu=0.055,
            zu2=0.0,
            v1dzu=0.0050000000000000001,
            zu3=0.0,
            bldzu=0.55000000000000004,
            zu4=0.0,
            v2dzu=0.85000000000000009,
            zu5=0.0,
            shdzu=0.20000000000000001,
            zu6=0.0,
            v3dzu=13.640000000000001,
            zu7=0.0,
            fwmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v1matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            blmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v2matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            shmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v3matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            chmatv=np.array(
                np.array(
                    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
                ),
            ).T,
            expected_a_fw_total=226.5088303238241,
            expected_r1=3.5,
            expected_r2=3.5550000000000002,
            expected_r3=3.5600000000000001,
            expected_r4=4.1100000000000003,
            expected_r5=6.3000000000000007,
            expected_r6=6.5000000000000009,
            expected_r7=10.0,
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
            expected_zu7=19.0,
            expected_fwmatv=np.array([
                [
                    0.41446560639912183,
                    0.41446560639912183,
                    0.41446560639912183,
                    0.0,
                    0.0,
                    0.0,
                    7.874846521583314,
                    7.874846521583314,
                    7.874846521583314,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_v1matv=np.array([
                [
                    0.7661295975557304,
                    0.7661295975557304,
                    0.7661295975557304,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_blmatv=np.array([
                [
                    5.7053341062470615,
                    5.7053341062470615,
                    5.7053341062470615,
                    0.0,
                    0.0,
                    0.0,
                    51.34800695622355,
                    51.34800695622355,
                    51.34800695622355,
                ],
                [
                    0.0,
                    0.0,
                    0.0,
                    5.837490699007461,
                    5.837490699007461,
                    5.837490699007461,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    10.348278966422315,
                    10.348278966422315,
                    10.348278966422315,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ]),
            expected_v2matv=np.array([
                [
                    702.6090256928657,
                    702.6090256928657,
                    702.6090256928657,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_shmatv=np.array([
                [
                    4.286640343970205,
                    4.286640343970205,
                    4.286640343970205,
                    16.289233307086782,
                    16.289233307086782,
                    16.289233307086782,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    17.653394518684465,
                    17.653394518684465,
                    17.653394518684465,
                ],
                [
                    4.413348629671108,
                    4.413348629671108,
                    4.413348629671108,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ]),
            expected_v3matv=np.array([
                [
                    2033.7963980993259,
                    2033.7963980993259,
                    2033.7963980993259,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_chmatv=np.array([
                261.69466804402975,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ]),
        ),
    ),
)
def test_genbld(genbldparam, monkeypatch, ife):
    """
    Automatically generated Regression Unit Test for genbld.

    This test was generated using data from ife/generic.IN.DAT.

    :param genbldparam: the data used to mock and assert in this test.
    :type genbldparam: genbldparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(first_wall_variables, "a_fw_total", genbldparam.a_fw_total)
    monkeypatch.setattr(ife_variables, "ifetyp", genbldparam.ifetyp)
    monkeypatch.setattr(ife_variables, "chrad", genbldparam.chrad)
    monkeypatch.setattr(ife_variables, "r1", genbldparam.r1)
    monkeypatch.setattr(ife_variables, "fwdr", genbldparam.fwdr)
    monkeypatch.setattr(ife_variables, "r2", genbldparam.r2)
    monkeypatch.setattr(ife_variables, "v1dr", genbldparam.v1dr)
    monkeypatch.setattr(ife_variables, "r3", genbldparam.r3)
    monkeypatch.setattr(ife_variables, "bldr", genbldparam.bldr)
    monkeypatch.setattr(ife_variables, "r4", genbldparam.r4)
    monkeypatch.setattr(ife_variables, "v2dr", genbldparam.v2dr)
    monkeypatch.setattr(ife_variables, "r5", genbldparam.r5)
    monkeypatch.setattr(ife_variables, "shdr", genbldparam.shdr)
    monkeypatch.setattr(ife_variables, "r6", genbldparam.r6)
    monkeypatch.setattr(ife_variables, "v3dr", genbldparam.v3dr)
    monkeypatch.setattr(ife_variables, "r7", genbldparam.r7)
    monkeypatch.setattr(ife_variables, "zl7", genbldparam.zl7)
    monkeypatch.setattr(ife_variables, "v3dzl", genbldparam.v3dzl)
    monkeypatch.setattr(ife_variables, "zl6", genbldparam.zl6)
    monkeypatch.setattr(ife_variables, "shdzl", genbldparam.shdzl)
    monkeypatch.setattr(ife_variables, "zl5", genbldparam.zl5)
    monkeypatch.setattr(ife_variables, "v2dzl", genbldparam.v2dzl)
    monkeypatch.setattr(ife_variables, "zl4", genbldparam.zl4)
    monkeypatch.setattr(ife_variables, "bldzl", genbldparam.bldzl)
    monkeypatch.setattr(ife_variables, "zl3", genbldparam.zl3)
    monkeypatch.setattr(ife_variables, "v1dzl", genbldparam.v1dzl)
    monkeypatch.setattr(ife_variables, "zl2", genbldparam.zl2)
    monkeypatch.setattr(ife_variables, "fwdzl", genbldparam.fwdzl)
    monkeypatch.setattr(ife_variables, "zl1", genbldparam.zl1)
    monkeypatch.setattr(ife_variables, "chdzl", genbldparam.chdzl)
    monkeypatch.setattr(ife_variables, "chdzu", genbldparam.chdzu)
    monkeypatch.setattr(ife_variables, "zu1", genbldparam.zu1)
    monkeypatch.setattr(ife_variables, "fwdzu", genbldparam.fwdzu)
    monkeypatch.setattr(ife_variables, "zu2", genbldparam.zu2)
    monkeypatch.setattr(ife_variables, "v1dzu", genbldparam.v1dzu)
    monkeypatch.setattr(ife_variables, "zu3", genbldparam.zu3)
    monkeypatch.setattr(ife_variables, "bldzu", genbldparam.bldzu)
    monkeypatch.setattr(ife_variables, "zu4", genbldparam.zu4)
    monkeypatch.setattr(ife_variables, "v2dzu", genbldparam.v2dzu)
    monkeypatch.setattr(ife_variables, "zu5", genbldparam.zu5)
    monkeypatch.setattr(ife_variables, "shdzu", genbldparam.shdzu)
    monkeypatch.setattr(ife_variables, "zu6", genbldparam.zu6)
    monkeypatch.setattr(ife_variables, "v3dzu", genbldparam.v3dzu)
    monkeypatch.setattr(ife_variables, "zu7", genbldparam.zu7)
    monkeypatch.setattr(ife_variables, "fwmatv", genbldparam.fwmatv)
    monkeypatch.setattr(ife_variables, "v1matv", genbldparam.v1matv)
    monkeypatch.setattr(ife_variables, "blmatv", genbldparam.blmatv)
    monkeypatch.setattr(ife_variables, "v2matv", genbldparam.v2matv)
    monkeypatch.setattr(ife_variables, "shmatv", genbldparam.shmatv)
    monkeypatch.setattr(ife_variables, "v3matv", genbldparam.v3matv)
    monkeypatch.setattr(ife_variables, "chmatv", genbldparam.chmatv)

    ife.genbld()

    assert first_wall_variables.a_fw_total == pytest.approx(
        genbldparam.expected_a_fw_total
    )
    assert ife_variables.r1 == pytest.approx(genbldparam.expected_r1)
    assert ife_variables.r2 == pytest.approx(genbldparam.expected_r2)
    assert ife_variables.r3 == pytest.approx(genbldparam.expected_r3)
    assert ife_variables.r4 == pytest.approx(genbldparam.expected_r4)
    assert ife_variables.r5 == pytest.approx(genbldparam.expected_r5)
    assert ife_variables.r6 == pytest.approx(genbldparam.expected_r6)
    assert ife_variables.r7 == pytest.approx(genbldparam.expected_r7)
    assert ife_variables.zl7 == pytest.approx(genbldparam.expected_zl7)
    assert ife_variables.zl6 == pytest.approx(genbldparam.expected_zl6)
    assert ife_variables.zl5 == pytest.approx(genbldparam.expected_zl5)
    assert ife_variables.zl4 == pytest.approx(genbldparam.expected_zl4)
    assert ife_variables.zl3 == pytest.approx(genbldparam.expected_zl3)
    assert ife_variables.zl2 == pytest.approx(genbldparam.expected_zl2)
    assert ife_variables.zl1 == pytest.approx(genbldparam.expected_zl1)
    assert ife_variables.zu1 == pytest.approx(genbldparam.expected_zu1)
    assert ife_variables.zu2 == pytest.approx(genbldparam.expected_zu2)
    assert ife_variables.zu3 == pytest.approx(genbldparam.expected_zu3)
    assert ife_variables.zu4 == pytest.approx(genbldparam.expected_zu4)
    assert ife_variables.zu5 == pytest.approx(genbldparam.expected_zu5)
    assert ife_variables.zu6 == pytest.approx(genbldparam.expected_zu6)
    assert ife_variables.zu7 == pytest.approx(genbldparam.expected_zu7)
    assert ife_variables.fwmatv == pytest.approx(genbldparam.expected_fwmatv)
    assert ife_variables.v1matv == pytest.approx(genbldparam.expected_v1matv)
    assert ife_variables.blmatv == pytest.approx(genbldparam.expected_blmatv)
    assert ife_variables.v2matv == pytest.approx(genbldparam.expected_v2matv)
    assert ife_variables.shmatv == pytest.approx(genbldparam.expected_shmatv)
    assert ife_variables.v3matv == pytest.approx(genbldparam.expected_v3matv)
    assert ife_variables.chmatv == pytest.approx(genbldparam.expected_chmatv)


class Ifepw1Param(NamedTuple):
    f_p_blkt_multiplication: Any = None
    fhole: Any = None
    p_blkt_nuclear_heat_total_mw: Any = None
    p_shld_nuclear_heat_mw: Any = None
    pnucloss: Any = None
    priheat: Any = None
    p_plant_primary_heat_mw: Any = None
    p_fw_div_heat_deposited_mw: Any = None
    n_primary_heat_exchangers: Any = None
    p_hcd_electric_total_mw: Any = None
    p_hcd_electric_loss_mw: Any = None
    p_cryo_plant_electric_mw: Any = None
    helpow: Any = None
    pdrive: Any = None
    ifetyp: Any = None
    etadrv: Any = None
    pifecr: Any = None
    p_fusion_total_mw: Any = None
    expected_p_blkt_nuclear_heat_total_mw: Any = None
    expected_priheat: Any = None
    expected_p_plant_primary_heat_mw: Any = None
    expected_p_fw_div_heat_deposited_mw: Any = None
    expected_nphx: Any = None
    expected_p_hcd_electric_total_mw: Any = None
    expected_p_hcd_electric_loss_mw: Any = None
    expected_p_cryo_plant_electric_mw: Any = None
    expected_helpow: Any = None


@pytest.mark.parametrize(
    "ifepw1param",
    (
        Ifepw1Param(
            f_p_blkt_multiplication=1.26,
            fhole=0.0,
            p_blkt_nuclear_heat_total_mw=0.0,
            p_shld_nuclear_heat_mw=0.0,
            pnucloss=0.0,
            priheat=0.0,
            p_plant_primary_heat_mw=0.0,
            p_fw_div_heat_deposited_mw=0.0,
            n_primary_heat_exchangers=0.0,
            p_hcd_electric_total_mw=0.0,
            p_hcd_electric_loss_mw=0.0,
            p_cryo_plant_electric_mw=0.0,
            helpow=0.0,
            pdrive=23100000,
            ifetyp=1,
            etadrv=0.28199999999999997,
            pifecr=10.0,
            p_fusion_total_mw=2009.6999999999998,
            expected_p_blkt_nuclear_heat_total_mw=1924.4887199999998,
            expected_priheat=2532.2219999999998,
            expected_p_plant_primary_heat_mw=2532.2219999999998,
            expected_p_fw_div_heat_deposited_mw=607.73327999999992,
            expected_nphx=3.0,
            expected_p_hcd_electric_total_mw=81.914893617021278,
            expected_p_hcd_electric_loss_mw=58.814893617021283,
            expected_p_cryo_plant_electric_mw=10.0,
            expected_helpow=20266.75905075,
        ),
    ),
)
def test_ifepw1(ifepw1param, monkeypatch, ife):
    """
    Automatically generated Regression Unit Test for ifepw1.

    This test was generated using data from ife/IFE.IN.DAT.

    :param ifepw1param: the data used to mock and assert in this test.
    :type ifepw1param: ifepw1param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        fwbs_variables, "f_p_blkt_multiplication", ifepw1param.f_p_blkt_multiplication
    )
    monkeypatch.setattr(fwbs_variables, "fhole", ifepw1param.fhole)
    monkeypatch.setattr(
        fwbs_variables,
        "p_blkt_nuclear_heat_total_mw",
        ifepw1param.p_blkt_nuclear_heat_total_mw,
    )
    monkeypatch.setattr(
        fwbs_variables, "p_shld_nuclear_heat_mw", ifepw1param.p_shld_nuclear_heat_mw
    )
    monkeypatch.setattr(fwbs_variables, "pnucloss", ifepw1param.pnucloss)
    monkeypatch.setattr(heat_transport_variables, "priheat", ifepw1param.priheat)
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_primary_heat_mw",
        ifepw1param.p_plant_primary_heat_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_fw_div_heat_deposited_mw",
        ifepw1param.p_fw_div_heat_deposited_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "n_primary_heat_exchangers",
        ifepw1param.n_primary_heat_exchangers,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_total_mw",
        ifepw1param.p_hcd_electric_total_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_loss_mw",
        ifepw1param.p_hcd_electric_loss_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_cryo_plant_electric_mw",
        ifepw1param.p_cryo_plant_electric_mw,
    )
    monkeypatch.setattr(heat_transport_variables, "helpow", ifepw1param.helpow)
    monkeypatch.setattr(ife_variables, "pdrive", ifepw1param.pdrive)
    monkeypatch.setattr(ife_variables, "ifetyp", ifepw1param.ifetyp)
    monkeypatch.setattr(ife_variables, "etadrv", ifepw1param.etadrv)
    monkeypatch.setattr(ife_variables, "pifecr", ifepw1param.pifecr)
    monkeypatch.setattr(
        physics_variables, "p_fusion_total_mw", ifepw1param.p_fusion_total_mw
    )

    ife.ifepw1()

    assert fwbs_variables.p_blkt_nuclear_heat_total_mw == pytest.approx(
        ifepw1param.expected_p_blkt_nuclear_heat_total_mw
    )
    assert heat_transport_variables.priheat == pytest.approx(
        ifepw1param.expected_priheat
    )
    assert heat_transport_variables.p_plant_primary_heat_mw == pytest.approx(
        ifepw1param.expected_p_plant_primary_heat_mw
    )
    assert heat_transport_variables.p_fw_div_heat_deposited_mw == pytest.approx(
        ifepw1param.expected_p_fw_div_heat_deposited_mw
    )
    assert heat_transport_variables.n_primary_heat_exchangers == pytest.approx(
        ifepw1param.expected_nphx
    )
    assert heat_transport_variables.p_hcd_electric_total_mw == pytest.approx(
        ifepw1param.expected_p_hcd_electric_total_mw
    )
    assert heat_transport_variables.p_hcd_electric_loss_mw == pytest.approx(
        ifepw1param.expected_p_hcd_electric_loss_mw
    )
    assert heat_transport_variables.p_cryo_plant_electric_mw == pytest.approx(
        ifepw1param.expected_p_cryo_plant_electric_mw
    )
    assert heat_transport_variables.helpow == pytest.approx(ifepw1param.expected_helpow)


class Bld2019Param(NamedTuple):
    a_fw_total: Any = None
    trcl: Any = None
    stcl: Any = None
    tbr: Any = None
    f_p_blkt_multiplication: Any = None
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
    expected_a_fw_total: Any = None
    expected_tbr: Any = None
    expected_emult: Any = None
    expected_r1: Any = None
    expected_r2: Any = None
    expected_r3: Any = None
    expected_bldr: Any = None
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
    expected_v3dzu: Any = None
    expected_zu7: Any = None
    expected_blmatv: Any = None
    expected_v2matv: Any = None
    expected_shmatv: Any = None
    expected_v3matv: Any = None
    expected_chmatv: Any = None


@pytest.mark.parametrize(
    "bld2019param",
    (
        Bld2019Param(
            a_fw_total=0.0,
            trcl=1.0,
            stcl=3.0,
            tbr=0.0,
            f_p_blkt_multiplication=1.26,
            ifetyp=4,
            chrad=3.5,
            r1=0.0,
            fwdr=0.0,
            r2=0.0,
            v1dr=0.0,
            r3=0.0,
            bldr=0.55000000000000004,
            r4=0.0,
            v2dr=2.1899999999999999,
            r5=0.0,
            shdr=0.20000000000000001,
            r6=0.0,
            v3dr=3.5,
            r7=0.0,
            zl7=0.0,
            v3dzl=0.0,
            zl6=0.0,
            shdzl=0.35000000000000003,
            zl5=0.0,
            v2dzl=0.0,
            zl4=0.0,
            bldzl=0.65000000000000013,
            zl3=0.0,
            v1dzl=0.0,
            zl2=0.0,
            fwdzl=0.0,
            zl1=0.0,
            chdzl=3.1000000000000001,
            chdzu=3.7000000000000002,
            zu1=0.0,
            fwdzu=0.0,
            zu2=0.0,
            v1dzu=0.0,
            zu3=0.0,
            bldzu=0.55000000000000004,
            zu4=0.0,
            v2dzu=0.0,
            zu5=0.0,
            shdzu=0.20000000000000001,
            zu6=0.0,
            v3dzu=13.640000000000001,
            zu7=0.0,
            fwmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v1matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            blmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v2matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            shmatv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            v3matv=np.array(
                (
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                    (0.0, 0.0, 0.0),
                ),
            ).T,
            chmatv=np.array(
                (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            ).T,
            expected_a_fw_total=36.573165036030936,
            expected_tbr=1.3991938274222466,
            expected_emult=0.87599606521170326,
            expected_r1=3.5,
            expected_r2=3.5,
            expected_r3=3.5,
            expected_bldr=3.0547221599995269,
            expected_r4=6.5547221599995265,
            expected_r5=8.744722159999526,
            expected_r6=8.9447221599995252,
            expected_r7=12.444722159999525,
            expected_zl7=4.0999999999999996,
            expected_zl6=4.0999999999999996,
            expected_zl5=3.75,
            expected_zl4=3.75,
            expected_zl3=3.1000000000000001,
            expected_zl2=3.1000000000000001,
            expected_zl1=3.1000000000000001,
            expected_zu1=3.7000000000000002,
            expected_zu2=3.7000000000000002,
            expected_zu3=3.7000000000000002,
            expected_zu4=4.25,
            expected_zu5=4.25,
            expected_zu6=4.4500000000000002,
            expected_v3dzu=18.590999999999998,
            expected_zu7=23.040999999999997,
            expected_blmatv=np.array([
                [
                    441.3488374179947,
                    32.80730799732719,
                    32.80730799732719,
                    0.0,
                    0.0,
                    0.0,
                    295.2657719759447,
                    295.2657719759447,
                    214.79732252854913,
                ],
                [
                    0.0,
                    0.0,
                    0.0,
                    22.192896602036015,
                    22.192896602036015,
                    22.192896602036015,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    46.84643251850827,
                    46.84643251850827,
                    46.84643251850827,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ]),
            expected_v2matv=np.array([
                [
                    715.7783572702281,
                    715.7783572702281,
                    715.7783572702281,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_shmatv=np.array([
                [
                    4.445842309995582,
                    4.445842309995582,
                    4.445842309995582,
                    16.89420077798321,
                    16.89420077798321,
                    16.89420077798321,
                    0.0,
                    0.0,
                    0.0,
                ],
                [
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    34.73839406353055,
                    34.73839406353055,
                    34.73839406353055,
                ],
                [
                    8.357477981712734,
                    8.357477981712734,
                    8.357477981712734,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
            ]),
            expected_v3matv=np.array([
                [
                    2010.8678816698186,
                    2010.8678816698186,
                    2010.8678816698186,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ]),
            expected_chmatv=np.array([
                261.69466804402975,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ]),
        ),
    ),
)
def test_bld2019(bld2019param, monkeypatch, ife):
    """
    Automatically generated Regression Unit Test for bld2019.

    This test was generated using data from ife/2019.IN.DAT.

    :param bld2019param: the data used to mock and assert in this test.
    :type bld2019param: bld2019param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(first_wall_variables, "a_fw_total", bld2019param.a_fw_total)
    monkeypatch.setattr(buildings_variables, "trcl", bld2019param.trcl)
    monkeypatch.setattr(buildings_variables, "stcl", bld2019param.stcl)
    monkeypatch.setattr(fwbs_variables, "tbr", bld2019param.tbr)
    monkeypatch.setattr(
        fwbs_variables, "f_p_blkt_multiplication", bld2019param.f_p_blkt_multiplication
    )
    monkeypatch.setattr(ife_variables, "ifetyp", bld2019param.ifetyp)
    monkeypatch.setattr(ife_variables, "chrad", bld2019param.chrad)
    monkeypatch.setattr(ife_variables, "r1", bld2019param.r1)
    monkeypatch.setattr(ife_variables, "fwdr", bld2019param.fwdr)
    monkeypatch.setattr(ife_variables, "r2", bld2019param.r2)
    monkeypatch.setattr(ife_variables, "v1dr", bld2019param.v1dr)
    monkeypatch.setattr(ife_variables, "r3", bld2019param.r3)
    monkeypatch.setattr(ife_variables, "bldr", bld2019param.bldr)
    monkeypatch.setattr(ife_variables, "r4", bld2019param.r4)
    monkeypatch.setattr(ife_variables, "v2dr", bld2019param.v2dr)
    monkeypatch.setattr(ife_variables, "r5", bld2019param.r5)
    monkeypatch.setattr(ife_variables, "shdr", bld2019param.shdr)
    monkeypatch.setattr(ife_variables, "r6", bld2019param.r6)
    monkeypatch.setattr(ife_variables, "v3dr", bld2019param.v3dr)
    monkeypatch.setattr(ife_variables, "r7", bld2019param.r7)
    monkeypatch.setattr(ife_variables, "zl7", bld2019param.zl7)
    monkeypatch.setattr(ife_variables, "v3dzl", bld2019param.v3dzl)
    monkeypatch.setattr(ife_variables, "zl6", bld2019param.zl6)
    monkeypatch.setattr(ife_variables, "shdzl", bld2019param.shdzl)
    monkeypatch.setattr(ife_variables, "zl5", bld2019param.zl5)
    monkeypatch.setattr(ife_variables, "v2dzl", bld2019param.v2dzl)
    monkeypatch.setattr(ife_variables, "zl4", bld2019param.zl4)
    monkeypatch.setattr(ife_variables, "bldzl", bld2019param.bldzl)
    monkeypatch.setattr(ife_variables, "zl3", bld2019param.zl3)
    monkeypatch.setattr(ife_variables, "v1dzl", bld2019param.v1dzl)
    monkeypatch.setattr(ife_variables, "zl2", bld2019param.zl2)
    monkeypatch.setattr(ife_variables, "fwdzl", bld2019param.fwdzl)
    monkeypatch.setattr(ife_variables, "zl1", bld2019param.zl1)
    monkeypatch.setattr(ife_variables, "chdzl", bld2019param.chdzl)
    monkeypatch.setattr(ife_variables, "chdzu", bld2019param.chdzu)
    monkeypatch.setattr(ife_variables, "zu1", bld2019param.zu1)
    monkeypatch.setattr(ife_variables, "fwdzu", bld2019param.fwdzu)
    monkeypatch.setattr(ife_variables, "zu2", bld2019param.zu2)
    monkeypatch.setattr(ife_variables, "v1dzu", bld2019param.v1dzu)
    monkeypatch.setattr(ife_variables, "zu3", bld2019param.zu3)
    monkeypatch.setattr(ife_variables, "bldzu", bld2019param.bldzu)
    monkeypatch.setattr(ife_variables, "zu4", bld2019param.zu4)
    monkeypatch.setattr(ife_variables, "v2dzu", bld2019param.v2dzu)
    monkeypatch.setattr(ife_variables, "zu5", bld2019param.zu5)
    monkeypatch.setattr(ife_variables, "shdzu", bld2019param.shdzu)
    monkeypatch.setattr(ife_variables, "zu6", bld2019param.zu6)
    monkeypatch.setattr(ife_variables, "v3dzu", bld2019param.v3dzu)
    monkeypatch.setattr(ife_variables, "zu7", bld2019param.zu7)
    monkeypatch.setattr(ife_variables, "fwmatv", bld2019param.fwmatv)
    monkeypatch.setattr(ife_variables, "v1matv", bld2019param.v1matv)
    monkeypatch.setattr(ife_variables, "blmatv", bld2019param.blmatv)
    monkeypatch.setattr(ife_variables, "v2matv", bld2019param.v2matv)
    monkeypatch.setattr(ife_variables, "shmatv", bld2019param.shmatv)
    monkeypatch.setattr(ife_variables, "v3matv", bld2019param.v3matv)
    monkeypatch.setattr(ife_variables, "chmatv", bld2019param.chmatv)

    ife.bld2019()

    assert first_wall_variables.a_fw_total == pytest.approx(
        bld2019param.expected_a_fw_total
    )
    assert fwbs_variables.tbr == pytest.approx(bld2019param.expected_tbr)
    assert fwbs_variables.f_p_blkt_multiplication == pytest.approx(
        bld2019param.expected_emult
    )
    assert ife_variables.r1 == pytest.approx(bld2019param.expected_r1)
    assert ife_variables.r2 == pytest.approx(bld2019param.expected_r2)
    assert ife_variables.r3 == pytest.approx(bld2019param.expected_r3)
    assert ife_variables.bldr == pytest.approx(bld2019param.expected_bldr)
    assert ife_variables.r4 == pytest.approx(bld2019param.expected_r4)
    assert ife_variables.r5 == pytest.approx(bld2019param.expected_r5)
    assert ife_variables.r6 == pytest.approx(bld2019param.expected_r6)
    assert ife_variables.r7 == pytest.approx(bld2019param.expected_r7)
    assert ife_variables.zl7 == pytest.approx(bld2019param.expected_zl7)
    assert ife_variables.zl6 == pytest.approx(bld2019param.expected_zl6)
    assert ife_variables.zl5 == pytest.approx(bld2019param.expected_zl5)
    assert ife_variables.zl4 == pytest.approx(bld2019param.expected_zl4)
    assert ife_variables.zl3 == pytest.approx(bld2019param.expected_zl3)
    assert ife_variables.zl2 == pytest.approx(bld2019param.expected_zl2)
    assert ife_variables.zl1 == pytest.approx(bld2019param.expected_zl1)
    assert ife_variables.zu1 == pytest.approx(bld2019param.expected_zu1)
    assert ife_variables.zu2 == pytest.approx(bld2019param.expected_zu2)
    assert ife_variables.zu3 == pytest.approx(bld2019param.expected_zu3)
    assert ife_variables.zu4 == pytest.approx(bld2019param.expected_zu4)
    assert ife_variables.zu5 == pytest.approx(bld2019param.expected_zu5)
    assert ife_variables.zu6 == pytest.approx(bld2019param.expected_zu6)
    assert ife_variables.v3dzu == pytest.approx(bld2019param.expected_v3dzu)
    assert ife_variables.zu7 == pytest.approx(bld2019param.expected_zu7)
    assert ife_variables.blmatv == pytest.approx(bld2019param.expected_blmatv)
    assert ife_variables.v2matv == pytest.approx(bld2019param.expected_v2matv)
    assert ife_variables.shmatv == pytest.approx(bld2019param.expected_shmatv)
    assert ife_variables.v3matv == pytest.approx(bld2019param.expected_v3matv)
    assert ife_variables.chmatv == pytest.approx(bld2019param.expected_chmatv)


class IfeacpParam(NamedTuple):
    a_plant_floor_effective: Any = None
    p_plant_electric_base: Any = None
    pflux_plant_floor_electric: Any = None
    pacpmw: Any = None
    p_cryo_plant_electric_mw: Any = None
    vachtmw: Any = None
    p_tritium_plant_electric_mw: Any = None
    p_hcd_electric_total_mw: Any = None
    p_plant_electric_base_total_mw: Any = None
    tlvpmw: Any = None
    tdspmw: Any = None
    tfacmw: Any = None
    htpmw_ife: Any = None
    reprat: Any = None
    lipmw: Any = None
    ifetyp: Any = None
    expected_pacpmw: Any = None
    expected_fcsht: Any = None
    expected_tlvpmw: Any = None


@pytest.mark.parametrize(
    "ifeacpparam",
    (
        IfeacpParam(
            a_plant_floor_effective=128814.70697706047,
            p_plant_electric_base=5000000.0,
            pflux_plant_floor_electric=150.0,
            pacpmw=0.0,
            p_cryo_plant_electric_mw=10.0,
            vachtmw=0.5,
            p_tritium_plant_electric_mw=15.0,
            p_hcd_electric_total_mw=81.914893617021278,
            p_plant_electric_base_total_mw=0.0,
            tlvpmw=0.0,
            tdspmw=0.01,
            tfacmw=1.6656107044913124,
            htpmw_ife=10.0,
            reprat=4.6200000000000001,
            lipmw=0.0,
            ifetyp=1,
            expected_pacpmw=141.11271036807165,
            expected_fcsht=24.322206046559071,
            expected_tlvpmw=54.187816751050391,
        ),
    ),
)
def test_ifeacp(ifeacpparam, monkeypatch, ife):
    """
    Automatically generated Regression Unit Test for ifeacp.

    This test was generated using data from tests/regression/input_files/IFE.IN.DAT.

    :param ifeacpparam: the data used to mock and assert in this test.
    :type ifeacpparam: ifeacpparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        buildings_variables,
        "a_plant_floor_effective",
        ifeacpparam.a_plant_floor_effective,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_base",
        ifeacpparam.p_plant_electric_base,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "pflux_plant_floor_electric",
        ifeacpparam.pflux_plant_floor_electric,
    )
    monkeypatch.setattr(heat_transport_variables, "pacpmw", ifeacpparam.pacpmw)
    monkeypatch.setattr(
        heat_transport_variables,
        "p_cryo_plant_electric_mw",
        ifeacpparam.p_cryo_plant_electric_mw,
    )
    monkeypatch.setattr(heat_transport_variables, "vachtmw", ifeacpparam.vachtmw)
    monkeypatch.setattr(
        heat_transport_variables,
        "p_tritium_plant_electric_mw",
        ifeacpparam.p_tritium_plant_electric_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_total_mw",
        ifeacpparam.p_hcd_electric_total_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_base_total_mw",
        ifeacpparam.p_plant_electric_base_total_mw,
    )
    monkeypatch.setattr(heat_transport_variables, "tlvpmw", ifeacpparam.tlvpmw)
    monkeypatch.setattr(ife_variables, "tdspmw", ifeacpparam.tdspmw)
    monkeypatch.setattr(ife_variables, "tfacmw", ifeacpparam.tfacmw)
    monkeypatch.setattr(ife_variables, "htpmw_ife", ifeacpparam.htpmw_ife)
    monkeypatch.setattr(ife_variables, "reprat", ifeacpparam.reprat)
    monkeypatch.setattr(ife_variables, "lipmw", ifeacpparam.lipmw)
    monkeypatch.setattr(ife_variables, "ifetyp", ifeacpparam.ifetyp)

    ife.ifeacp()

    assert heat_transport_variables.pacpmw == pytest.approx(ifeacpparam.expected_pacpmw)
    assert heat_transport_variables.p_plant_electric_base_total_mw == pytest.approx(
        ifeacpparam.expected_fcsht
    )
    assert heat_transport_variables.tlvpmw == pytest.approx(ifeacpparam.expected_tlvpmw)


class IfebdgParam(NamedTuple):
    wrbi: Any = None
    rbwt: Any = None
    rbrt: Any = None
    fndt: Any = None
    trcl: Any = None
    hcwt: Any = None
    hccl: Any = None
    wgt2: Any = None
    stcl: Any = None
    pibv: Any = None
    a_plant_floor_effective: Any = None
    triv: Any = None
    conv: Any = None
    admv: Any = None
    shov: Any = None
    admvol: Any = None
    convol: Any = None
    elevol: Any = None
    rbvol: Any = None
    rmbvol: Any = None
    shovol: Any = None
    volrci: Any = None
    wsvol: Any = None
    volnucb: Any = None
    whtshld: Any = None
    helpow: Any = None
    zl7: Any = None
    zu7: Any = None
    r7: Any = None
    zl6: Any = None
    zu6: Any = None
    r6: Any = None
    expected_wrbi: Any = None
    expected_a_plant_floor_effective: Any = None
    expected_admvol: Any = None
    expected_convol: Any = None
    expected_elevol: Any = None
    expected_rbvol: Any = None
    expected_rmbvol: Any = None
    expected_shovol: Any = None
    expected_volrci: Any = None
    expected_wsvol: Any = None
    expected_volnucb: Any = None


@pytest.mark.parametrize(
    "ifebdgparam",
    (
        IfebdgParam(
            wrbi=0.0,
            rbwt=3.2000000000000002,
            rbrt=3.2000000000000002,
            fndt=2.0,
            trcl=1.0,
            hcwt=1.5,
            hccl=5.0,
            wgt2=100000.0,
            stcl=3.0,
            pibv=40000.0,
            a_plant_floor_effective=0.0,
            triv=40000.0,
            conv=60000.0,
            admv=100000.0,
            shov=100000.0,
            admvol=0.0,
            convol=0.0,
            elevol=0.0,
            rbvol=0.0,
            rmbvol=0.0,
            shovol=0.0,
            volrci=0.0,
            wsvol=0.0,
            volnucb=0.0,
            whtshld=1067310.9593707009,
            helpow=20277.29636048527,
            zl7=5.8499999999999996,
            zu7=19.0,
            r7=10.0,
            zl6=5.8499999999999996,
            zu6=5.3600000000000003,
            r6=6.5000000000000009,
            expected_wrbi=10.0,
            expected_a_plant_floor_effective=128814.70697706047,
            expected_admvol=100000.0,
            expected_convol=60000.0,
            expected_elevol=40000.0,
            expected_rbvol=20943.647999999997,
            expected_rmbvol=294010.79749999999,
            expected_shovol=100000.0,
            expected_volrci=9940.0,
            expected_wsvol=110101.88589999999,
            expected_volnucb=461884.59386236279,
        ),
    ),
)
def test_ifebdg(ifebdgparam, monkeypatch, ife):
    """
    Automatically generated Regression Unit Test for ifebdg.

    This test was generated using data from tests/regression/input_files/IFE.IN.DAT.

    :param ifebdgparam: the data used to mock and assert in this test.
    :type ifebdgparam: ifebdgparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(buildings_variables, "wrbi", ifebdgparam.wrbi)
    monkeypatch.setattr(buildings_variables, "rbwt", ifebdgparam.rbwt)
    monkeypatch.setattr(buildings_variables, "rbrt", ifebdgparam.rbrt)
    monkeypatch.setattr(buildings_variables, "fndt", ifebdgparam.fndt)
    monkeypatch.setattr(buildings_variables, "trcl", ifebdgparam.trcl)
    monkeypatch.setattr(buildings_variables, "hcwt", ifebdgparam.hcwt)
    monkeypatch.setattr(buildings_variables, "hccl", ifebdgparam.hccl)
    monkeypatch.setattr(buildings_variables, "wgt2", ifebdgparam.wgt2)
    monkeypatch.setattr(buildings_variables, "stcl", ifebdgparam.stcl)
    monkeypatch.setattr(buildings_variables, "pibv", ifebdgparam.pibv)
    monkeypatch.setattr(
        buildings_variables,
        "a_plant_floor_effective",
        ifebdgparam.a_plant_floor_effective,
    )
    monkeypatch.setattr(buildings_variables, "triv", ifebdgparam.triv)
    monkeypatch.setattr(buildings_variables, "conv", ifebdgparam.conv)
    monkeypatch.setattr(buildings_variables, "admv", ifebdgparam.admv)
    monkeypatch.setattr(buildings_variables, "shov", ifebdgparam.shov)
    monkeypatch.setattr(buildings_variables, "admvol", ifebdgparam.admvol)
    monkeypatch.setattr(buildings_variables, "convol", ifebdgparam.convol)
    monkeypatch.setattr(buildings_variables, "elevol", ifebdgparam.elevol)
    monkeypatch.setattr(buildings_variables, "rbvol", ifebdgparam.rbvol)
    monkeypatch.setattr(buildings_variables, "rmbvol", ifebdgparam.rmbvol)
    monkeypatch.setattr(buildings_variables, "shovol", ifebdgparam.shovol)
    monkeypatch.setattr(buildings_variables, "volrci", ifebdgparam.volrci)
    monkeypatch.setattr(buildings_variables, "wsvol", ifebdgparam.wsvol)
    monkeypatch.setattr(buildings_variables, "volnucb", ifebdgparam.volnucb)
    monkeypatch.setattr(fwbs_variables, "whtshld", ifebdgparam.whtshld)
    monkeypatch.setattr(heat_transport_variables, "helpow", ifebdgparam.helpow)
    monkeypatch.setattr(ife_variables, "zl7", ifebdgparam.zl7)
    monkeypatch.setattr(ife_variables, "zu7", ifebdgparam.zu7)
    monkeypatch.setattr(ife_variables, "r7", ifebdgparam.r7)
    monkeypatch.setattr(ife_variables, "zl6", ifebdgparam.zl6)
    monkeypatch.setattr(ife_variables, "zu6", ifebdgparam.zu6)
    monkeypatch.setattr(ife_variables, "r6", ifebdgparam.r6)

    ife.ifebdg()

    assert buildings_variables.wrbi == pytest.approx(ifebdgparam.expected_wrbi)
    assert buildings_variables.a_plant_floor_effective == pytest.approx(
        ifebdgparam.expected_a_plant_floor_effective
    )
    assert buildings_variables.admvol == pytest.approx(ifebdgparam.expected_admvol)
    assert buildings_variables.convol == pytest.approx(ifebdgparam.expected_convol)
    assert buildings_variables.elevol == pytest.approx(ifebdgparam.expected_elevol)
    assert buildings_variables.rbvol == pytest.approx(ifebdgparam.expected_rbvol)
    assert buildings_variables.rmbvol == pytest.approx(ifebdgparam.expected_rmbvol)
    assert buildings_variables.shovol == pytest.approx(ifebdgparam.expected_shovol)
    assert buildings_variables.volrci == pytest.approx(ifebdgparam.expected_volrci)
    assert buildings_variables.wsvol == pytest.approx(ifebdgparam.expected_wsvol)
    assert buildings_variables.volnucb == pytest.approx(ifebdgparam.expected_volnucb)


class Ifepw2Param(NamedTuple):
    ireactor: Any = None
    pnucloss: Any = None
    f_p_blkt_multiplication: Any = None
    tbr: Any = None
    p_blkt_nuclear_heat_total_mw: Any = None
    fachtmw: Any = None
    p_plant_electric_base_total_mw: Any = None
    p_plant_secondary_heat_mw: Any = None
    p_hcd_electric_loss_mw: Any = None
    vachtmw: Any = None
    p_tritium_plant_electric_mw: Any = None
    p_cryo_plant_electric_mw: Any = None
    p_plant_electric_gross_mw: Any = None
    p_plant_primary_heat_mw: Any = None
    eta_turbine: Any = None
    fgrosbop: Any = None
    p_plant_electric_recirc_mw: Any = None
    pacpmw: Any = None
    p_plant_electric_net_mw: Any = None
    p_hcd_electric_total_mw: Any = None
    p_fw_div_heat_deposited_mw: Any = None
    n_primary_heat_exchangers: Any = None
    tdspmw: Any = None
    tfacmw: Any = None
    htpmw_ife: Any = None
    fauxbop: Any = None
    ifetyp: Any = None
    taufall: Any = None
    expected_fachtmw: Any = None
    expected_p_plant_secondary_heat_mw: Any = None
    expected_p_plant_electric_gross_mw: Any = None
    expected_precircmw: Any = None
    expected_p_plant_electric_net_mw: Any = None


@pytest.mark.parametrize(
    "ifepw2param",
    (
        Ifepw2Param(
            ireactor=1,
            pnucloss=0.0,
            f_p_blkt_multiplication=1.26,
            tbr=0.0,
            p_blkt_nuclear_heat_total_mw=1924.4887199999998,
            fachtmw=0.0,
            p_plant_electric_base_total_mw=24.322206046559071,
            p_plant_secondary_heat_mw=0.0,
            p_hcd_electric_loss_mw=58.814893617021283,
            vachtmw=0.5,
            p_tritium_plant_electric_mw=15.0,
            p_cryo_plant_electric_mw=10.0,
            p_plant_electric_gross_mw=0.0,
            p_plant_primary_heat_mw=2532.2219999999998,
            eta_turbine=0.45000000000000001,
            fgrosbop=0.0,
            p_plant_electric_recirc_mw=0.0,
            pacpmw=141.11271036807165,
            p_plant_electric_net_mw=0.0,
            p_hcd_electric_total_mw=81.914893617021278,
            p_fw_div_heat_deposited_mw=607.73327999999992,
            n_primary_heat_exchangers=3,
            tdspmw=0.01,
            tfacmw=1.6656107044913124,
            htpmw_ife=10.0,
            fauxbop=0.0,
            ifetyp=1,
            taufall=0.0,
            expected_fachtmw=24.322206046559071,
            expected_p_plant_secondary_heat_mw=120.31271036807168,
            expected_p_plant_electric_gross_mw=1139.4999,
            expected_precircmw=141.11271036807165,
            expected_p_plant_electric_net_mw=998.38718963192832,
        ),
    ),
)
def test_ifepw2(ifepw2param, monkeypatch, ife):
    """
    Automatically generated Regression Unit Test for ifepw2.

    This test was generated using data from tests/regression/input_files/IFE.IN.DAT.

    :param ifepw2param: the data used to mock and assert in this test.
    :type ifepw2param: ifepw2param

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(cost_variables, "ireactor", ifepw2param.ireactor)
    monkeypatch.setattr(fwbs_variables, "pnucloss", ifepw2param.pnucloss)
    monkeypatch.setattr(
        fwbs_variables, "f_p_blkt_multiplication", ifepw2param.f_p_blkt_multiplication
    )
    monkeypatch.setattr(fwbs_variables, "tbr", ifepw2param.tbr)
    monkeypatch.setattr(
        fwbs_variables,
        "p_blkt_nuclear_heat_total_mw",
        ifepw2param.p_blkt_nuclear_heat_total_mw,
    )
    monkeypatch.setattr(heat_transport_variables, "fachtmw", ifepw2param.fachtmw)
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_base_total_mw",
        ifepw2param.p_plant_electric_base_total_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_secondary_heat_mw",
        ifepw2param.p_plant_secondary_heat_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_loss_mw",
        ifepw2param.p_hcd_electric_loss_mw,
    )
    monkeypatch.setattr(heat_transport_variables, "vachtmw", ifepw2param.vachtmw)
    monkeypatch.setattr(
        heat_transport_variables,
        "p_tritium_plant_electric_mw",
        ifepw2param.p_tritium_plant_electric_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_cryo_plant_electric_mw",
        ifepw2param.p_cryo_plant_electric_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_gross_mw",
        ifepw2param.p_plant_electric_gross_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_primary_heat_mw",
        ifepw2param.p_plant_primary_heat_mw,
    )
    monkeypatch.setattr(heat_transport_variables, "eta_turbine", ifepw2param.eta_turbine)
    monkeypatch.setattr(heat_transport_variables, "fgrosbop", ifepw2param.fgrosbop)
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_recirc_mw",
        ifepw2param.p_plant_electric_recirc_mw,
    )
    monkeypatch.setattr(heat_transport_variables, "pacpmw", ifepw2param.pacpmw)
    monkeypatch.setattr(
        heat_transport_variables,
        "p_plant_electric_net_mw",
        ifepw2param.p_plant_electric_net_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_hcd_electric_total_mw",
        ifepw2param.p_hcd_electric_total_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "p_fw_div_heat_deposited_mw",
        ifepw2param.p_fw_div_heat_deposited_mw,
    )
    monkeypatch.setattr(
        heat_transport_variables,
        "n_primary_heat_exchangers",
        ifepw2param.n_primary_heat_exchangers,
    )
    monkeypatch.setattr(ife_variables, "tdspmw", ifepw2param.tdspmw)
    monkeypatch.setattr(ife_variables, "tfacmw", ifepw2param.tfacmw)
    monkeypatch.setattr(ife_variables, "htpmw_ife", ifepw2param.htpmw_ife)
    monkeypatch.setattr(ife_variables, "fauxbop", ifepw2param.fauxbop)
    monkeypatch.setattr(ife_variables, "ifetyp", ifepw2param.ifetyp)
    monkeypatch.setattr(ife_variables, "taufall", ifepw2param.taufall)

    ife.ifepw2()

    assert heat_transport_variables.fachtmw == pytest.approx(
        ifepw2param.expected_fachtmw
    )
    assert heat_transport_variables.p_plant_secondary_heat_mw == pytest.approx(
        ifepw2param.expected_p_plant_secondary_heat_mw
    )
    assert heat_transport_variables.p_plant_electric_gross_mw == pytest.approx(
        ifepw2param.expected_p_plant_electric_gross_mw
    )
    assert heat_transport_variables.p_plant_electric_recirc_mw == pytest.approx(
        ifepw2param.expected_precircmw
    )
    assert heat_transport_variables.p_plant_electric_net_mw == pytest.approx(
        ifepw2param.expected_p_plant_electric_net_mw
    )
