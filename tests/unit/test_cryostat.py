from typing import Any, NamedTuple

import numpy as np
import pytest

from process.cryostat import Cryostat
from process.data_structure import (
    blanket_library,
    build_variables,
    buildings_variables,
    fwbs_variables,
    pfcoil_variables,
)


@pytest.fixture
def cryostat_fixture():
    """Provides BlanketLibrary object for testing.

    :returns: initialised BlanketLibrary object
    :rtype: process.blanket_library.BlanketLibrary
    """
    return Cryostat()


class ExternalCryoGeometryParam(NamedTuple):
    f_z_cryostat: Any = None
    z_tf_inside_half: Any = None
    dr_tf_inboard: Any = None
    dr_cryostat: Any = None
    r_cryostat_inboard: Any = None
    dr_pf_cryostat: Any = None
    z_cryostat_half_inside: Any = None
    vol_cryostat: Any = None
    m_vv: Any = None
    vol_vv: Any = None
    den_steel: Any = None
    dewmkg: Any = None
    r_pf_coil_outer: Any = None
    z_pf_coil_upper: Any = None
    dz_tf_cryostat: Any = None
    dz_pf_cryostat: Any = None
    expected_r_cryostat_inboard: Any = None
    expected_z_cryostat_half_inside: Any = None
    expected_vol_cryostat: Any = None
    expected_vvmass: Any = None
    expected_dewmkg: Any = None
    expected_dz_tf_cryostat: Any = None
    expected_dz_pf_cryostat: Any = None


@pytest.mark.parametrize(
    "externalcryogeometryparam",
    (
        ExternalCryoGeometryParam(
            f_z_cryostat=4.2679999999999998,
            z_tf_inside_half=8.8182171641274945,
            dr_tf_inboard=0.92672586247397692,
            dr_cryostat=0.15000000000000002,
            r_cryostat_inboard=0,
            dr_pf_cryostat=0.5,
            z_cryostat_half_inside=0,
            vol_cryostat=0,
            m_vv=0,
            vol_vv=1016.2876250857248,
            den_steel=7800,
            dewmkg=0,
            r_pf_coil_outer=np.array(
                np.array(
                    (
                        6.1290994712971543,
                        6.2110624909068086,
                        17.305470903073743,
                        17.305470903073743,
                        15.620546715016166,
                        15.620546715016166,
                        2.5506597842255361,
                        10.666666666666666,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            z_pf_coil_upper=np.array(
                np.array(
                    (
                        9.9154920004377978,
                        -11.249338850841614,
                        3.2350365669570316,
                        -3.2350365669570316,
                        7.8723998771612473,
                        -7.8723998771612473,
                        7.9363954477147454,
                        4.9333333333333336,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            dz_tf_cryostat=2.5,
            dz_pf_cryostat=0,
            expected_r_cryostat_inboard=17.805470903073743,
            expected_z_cryostat_half_inside=15.259637557000296,
            expected_vol_cryostat=818.1630389343372,
            expected_vvmass=7927043.4756686538,
            expected_dewmkg=14308715.179356484,
            expected_dz_tf_cryostat=5.514694530398824,
            expected_dz_pf_cryostat=5.3441455565624985,
        ),
    ),
)
def test_external_cryo_geometry(
    externalcryogeometryparam, monkeypatch, cryostat_fixture
):
    """
    Automatically generated Regression Unit Test for external_cryo_geometry.

    This test was generated using data from tests/regression/input_files/large_tokamak_eval.IN.DAT.

    :param externalcryogeometryparam: the data used to mock and assert in this test.
    :type externalcryogeometryparam: externalcryogeometryparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        build_variables, "f_z_cryostat", externalcryogeometryparam.f_z_cryostat
    )
    monkeypatch.setattr(
        build_variables, "z_tf_inside_half", externalcryogeometryparam.z_tf_inside_half
    )
    monkeypatch.setattr(
        build_variables, "dr_tf_inboard", externalcryogeometryparam.dr_tf_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_cryostat", externalcryogeometryparam.dr_cryostat
    )
    monkeypatch.setattr(
        fwbs_variables,
        "r_cryostat_inboard",
        externalcryogeometryparam.r_cryostat_inboard,
    )
    monkeypatch.setattr(
        fwbs_variables, "dr_pf_cryostat", externalcryogeometryparam.dr_pf_cryostat
    )
    monkeypatch.setattr(
        fwbs_variables,
        "z_cryostat_half_inside",
        externalcryogeometryparam.z_cryostat_half_inside,
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_cryostat", externalcryogeometryparam.vol_cryostat
    )
    monkeypatch.setattr(fwbs_variables, "m_vv", externalcryogeometryparam.m_vv)
    monkeypatch.setattr(fwbs_variables, "vol_vv", externalcryogeometryparam.vol_vv)
    monkeypatch.setattr(
        fwbs_variables, "den_steel", externalcryogeometryparam.den_steel
    )
    monkeypatch.setattr(fwbs_variables, "dewmkg", externalcryogeometryparam.dewmkg)
    monkeypatch.setattr(
        pfcoil_variables, "r_pf_coil_outer", externalcryogeometryparam.r_pf_coil_outer
    )
    monkeypatch.setattr(
        pfcoil_variables, "z_pf_coil_upper", externalcryogeometryparam.z_pf_coil_upper
    )
    monkeypatch.setattr(
        buildings_variables, "dz_tf_cryostat", externalcryogeometryparam.dz_tf_cryostat
    )
    monkeypatch.setattr(
        blanket_library, "dz_pf_cryostat", externalcryogeometryparam.dz_pf_cryostat
    )

    cryostat_fixture.external_cryo_geometry()

    assert fwbs_variables.r_cryostat_inboard == pytest.approx(
        externalcryogeometryparam.expected_r_cryostat_inboard
    )
    assert fwbs_variables.z_cryostat_half_inside == pytest.approx(
        externalcryogeometryparam.expected_z_cryostat_half_inside
    )
    assert fwbs_variables.vol_cryostat == pytest.approx(
        externalcryogeometryparam.expected_vol_cryostat
    )
    assert fwbs_variables.m_vv == pytest.approx(
        externalcryogeometryparam.expected_vvmass
    )
    assert fwbs_variables.dewmkg == pytest.approx(
        externalcryogeometryparam.expected_dewmkg
    )
    assert buildings_variables.dz_tf_cryostat == pytest.approx(
        externalcryogeometryparam.expected_dz_tf_cryostat
    )
    assert blanket_library.dz_pf_cryostat == pytest.approx(
        externalcryogeometryparam.expected_dz_pf_cryostat
    )
