from typing import Any, NamedTuple

import numpy as np
import pytest

from process.blanket_library import BlanketLibrary
from process.fortran import (
    blanket_library,
    build_variables,
    divertor_variables,
    fwbs_variables,
    physics_variables,
)
from process.fw import Fw


@pytest.fixture
def blanket_library_fixture():
    """Provides BlanketLibrary object for testing.

    :returns: initialised BlanketLibrary object
    :rtype: process.blanket_library.BlanketLibrary
    """
    return BlanketLibrary(Fw())


class PrimaryCoolantPropertiesParam(NamedTuple):
    i_fw_coolant_type: Any = None

    temp_fw_coolant_in: Any = None

    temp_fw_coolant_out: Any = None

    pres_fw_coolant: Any = None

    den_fw_coolant: Any = None

    cp_fw: Any = None

    cv_fw: Any = None

    i_blkt_coolant_type: Any = None

    temp_blkt_coolant_in: Any = None

    temp_blkt_coolant_out: Any = None

    pres_blkt_coolant: Any = None

    den_blkt_coolant: Any = None

    i_blkt_dual_coolant: Any = None

    visc_blkt_coolant: Any = None

    cp_bl: Any = None

    cv_bl: Any = None

    visc_fw_coolant: Any = None

    i_fw_blkt_shared_coolant: Any = None

    expected_den_fw_coolant: Any = None

    expected_cp_fw: Any = None

    expected_cv_fw: Any = None

    expected_den_blkt_coolant: Any = None

    expected_visc_blkt_coolant: Any = None

    expected_cp_bl: Any = None

    expected_cv_bl: Any = None

    expected_visc_fw_coolant: Any = None


@pytest.mark.parametrize(
    "primarycoolantpropertiesparam",
    (
        PrimaryCoolantPropertiesParam(
            i_fw_coolant_type="helium",
            temp_fw_coolant_in=573,
            temp_fw_coolant_out=773,
            pres_fw_coolant=8000000,
            den_fw_coolant=0,
            cp_fw=0,
            cv_fw=0,
            i_blkt_coolant_type=1,
            temp_blkt_coolant_in=573,
            temp_blkt_coolant_out=773,
            pres_blkt_coolant=8000000,
            den_blkt_coolant=0,
            i_blkt_dual_coolant=2,
            visc_blkt_coolant=0,
            cp_bl=0,
            cv_bl=0,
            visc_fw_coolant=0,
            i_fw_blkt_shared_coolant=0,
            expected_den_fw_coolant=5.6389735407435868,
            expected_cp_fw=5188.5588430173211,
            expected_cv_fw=3123.5687263525392,
            expected_den_blkt_coolant=5.6389735407435868,
            expected_visc_blkt_coolant=3.5036293160410249e-05,
            expected_cp_bl=5188.5588430173211,
            expected_cv_bl=3123.5687263525392,
            expected_visc_fw_coolant=3.5036293160410249e-05,
        ),
        PrimaryCoolantPropertiesParam(
            i_fw_coolant_type="helium",
            temp_fw_coolant_in=573,
            temp_fw_coolant_out=773,
            pres_fw_coolant=8000000,
            den_fw_coolant=5.6389735407435868,
            cp_fw=5188.5588430173211,
            cv_fw=3123.5687263525392,
            i_blkt_coolant_type=1,
            temp_blkt_coolant_in=573,
            temp_blkt_coolant_out=773,
            pres_blkt_coolant=8000000,
            den_blkt_coolant=5.6389735407435868,
            i_blkt_dual_coolant=2,
            visc_blkt_coolant=3.5036293160410249e-05,
            cp_bl=5188.5588430173211,
            cv_bl=3123.5687263525392,
            visc_fw_coolant=3.5036293160410249e-05,
            i_fw_blkt_shared_coolant=0,
            expected_den_fw_coolant=5.6389735407435868,
            expected_cp_fw=5188.5588430173211,
            expected_cv_fw=3123.5687263525392,
            expected_den_blkt_coolant=5.6389735407435868,
            expected_visc_blkt_coolant=3.5036293160410249e-05,
            expected_cp_bl=5188.5588430173211,
            expected_cv_bl=3123.5687263525392,
            expected_visc_fw_coolant=3.5036293160410249e-05,
        ),
    ),
)
def test_primary_coolant_properties(
    primarycoolantpropertiesparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for primary_coolant_properties.

    This test was generated using data from dcll/dcll_mms_lt_IN.DAT.

    :param primarycoolantpropertiesparam: the data used to mock and assert in this test.
    :type primarycoolantpropertiesparam: primarycoolantpropertiesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    # monkeypatch doesnt work for strings
    # but helium is the default
    # monkeypatch.setattr(
    #     fwbs_variables, "i_fw_coolant_type", primarycoolantpropertiesparam.i_fw_coolant_type
    # )

    monkeypatch.setattr(
        fwbs_variables,
        "temp_fw_coolant_in",
        primarycoolantpropertiesparam.temp_fw_coolant_in,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "temp_fw_coolant_out",
        primarycoolantpropertiesparam.temp_fw_coolant_out,
    )

    monkeypatch.setattr(
        fwbs_variables, "pres_fw_coolant", primarycoolantpropertiesparam.pres_fw_coolant
    )

    monkeypatch.setattr(
        fwbs_variables, "den_fw_coolant", primarycoolantpropertiesparam.den_fw_coolant
    )

    monkeypatch.setattr(fwbs_variables, "cp_fw", primarycoolantpropertiesparam.cp_fw)

    monkeypatch.setattr(fwbs_variables, "cv_fw", primarycoolantpropertiesparam.cv_fw)

    monkeypatch.setattr(
        fwbs_variables,
        "i_blkt_coolant_type",
        primarycoolantpropertiesparam.i_blkt_coolant_type,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "temp_blkt_coolant_in",
        primarycoolantpropertiesparam.temp_blkt_coolant_in,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "temp_blkt_coolant_out",
        primarycoolantpropertiesparam.temp_blkt_coolant_out,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "pres_blkt_coolant",
        primarycoolantpropertiesparam.pres_blkt_coolant,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "den_blkt_coolant",
        primarycoolantpropertiesparam.den_blkt_coolant,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "i_blkt_dual_coolant",
        primarycoolantpropertiesparam.i_blkt_dual_coolant,
    )

    monkeypatch.setattr(
        fwbs_variables,
        "visc_blkt_coolant",
        primarycoolantpropertiesparam.visc_blkt_coolant,
    )

    monkeypatch.setattr(fwbs_variables, "cp_bl", primarycoolantpropertiesparam.cp_bl)

    monkeypatch.setattr(fwbs_variables, "cv_bl", primarycoolantpropertiesparam.cv_bl)

    monkeypatch.setattr(
        fwbs_variables, "visc_fw_coolant", primarycoolantpropertiesparam.visc_fw_coolant
    )

    monkeypatch.setattr(
        fwbs_variables,
        "i_fw_blkt_shared_coolant",
        primarycoolantpropertiesparam.i_fw_blkt_shared_coolant,
    )

    blanket_library_fixture.primary_coolant_properties(output=False)

    assert fwbs_variables.den_fw_coolant == pytest.approx(
        primarycoolantpropertiesparam.expected_den_fw_coolant, rel=1e-4
    )

    assert fwbs_variables.cp_fw == pytest.approx(
        primarycoolantpropertiesparam.expected_cp_fw, rel=1e-4
    )

    assert fwbs_variables.cv_fw == pytest.approx(
        primarycoolantpropertiesparam.expected_cv_fw, rel=1e-4
    )

    assert fwbs_variables.den_blkt_coolant == pytest.approx(
        primarycoolantpropertiesparam.expected_den_blkt_coolant, rel=1e-4
    )

    assert fwbs_variables.visc_blkt_coolant == pytest.approx(
        primarycoolantpropertiesparam.expected_visc_blkt_coolant, rel=1e-4
    )

    assert fwbs_variables.cp_bl == pytest.approx(
        primarycoolantpropertiesparam.expected_cp_bl, rel=1e-4
    )

    assert fwbs_variables.cv_bl == pytest.approx(
        primarycoolantpropertiesparam.expected_cv_bl, rel=1e-4
    )

    assert fwbs_variables.visc_fw_coolant == pytest.approx(
        primarycoolantpropertiesparam.expected_visc_fw_coolant, rel=1e-4
    )


def test_deltap_tot_inboard_first_wall(monkeypatch, blanket_library_fixture):
    monkeypatch.setattr(fwbs_variables, "radius_fw_channel", 0.006)
    monkeypatch.setattr(fwbs_variables, "a_bz_liq", 0.22481)

    data = {
        "icoolpump": 1,
        "flow_velocity": 15.9,
        "flleng": 4,
        "no90": 2,
        "no180": 0,
        "coolant_density": 5.6,
        "coolant_dynamic_viscosity": 3.5e-5,
        "coolant_electrical_conductivity": 0.0,
        "pol_channel_length": 1.89,
        "nopolchan": 0,
        "label": "Inboard first wall",
    }

    assert (
        pytest.approx(blanket_library_fixture.deltap_tot(False, **data))
        == 5885.192672142268
    )


def test_deltap_tot_outboard_blanket_breeder_liquid(
    monkeypatch, blanket_library_fixture
):
    monkeypatch.setattr(fwbs_variables, "radius_fw_channel", 0.006)
    monkeypatch.setattr(fwbs_variables, "a_bz_liq", 0.22481)
    monkeypatch.setattr(fwbs_variables, "i_blkt_liquid_breeder_channel_type", 1)
    monkeypatch.setattr(fwbs_variables, "b_bz_liq", 0.11625)
    monkeypatch.setattr(fwbs_variables, "b_mag_blkt", [8.393, 3.868])
    monkeypatch.setattr(fwbs_variables, "bz_channel_conduct_liq", 833000)
    monkeypatch.setattr(fwbs_variables, "th_wall_secondary", 0.0125)

    data = {
        "icoolpump": 2,
        "flow_velocity": 0.06,
        "flleng": 4.7,
        "no90": 2,
        "no180": 1,
        "coolant_density": 9753.2,
        "coolant_dynamic_viscosity": 0.0017,
        "coolant_electrical_conductivity": 861800.8,
        "pol_channel_length": 1.89,
        "nopolchan": 0,
        "label": "Outboard blanket breeder liquid",
    }

    assert (
        pytest.approx(blanket_library_fixture.deltap_tot(False, **data))
        == 56.962742615936264
    )


def test_pumppower_primary_helium(monkeypatch, blanket_library_fixture):
    monkeypatch.setattr(fwbs_variables, "etaiso", 0.9)
    monkeypatch.setattr(fwbs_variables, "etaiso_liq", 0.85)

    data = {
        "icoolpump": 2,
        "temp_in": 570,
        "temp_out": 720,
        "pressure": 1700000,
        "pdrop": 303517.3,
        "mf": 35677.7,
        "primary_coolant_switch": 1,
        "coolant_density": 9753.25,
        "label": "Liquid Metal Breeder/Coolant",
    }

    assert (
        pytest.approx(blanket_library_fixture.pumppower(False, **data))
        == 1.8251284651310427
    )


def test_pumppower_secondary_pb_li(monkeypatch, blanket_library_fixture):
    monkeypatch.setattr(fwbs_variables, "etaiso", 0.9)
    monkeypatch.setattr(fwbs_variables, "etaiso_liq", 0.85)

    data = {
        "icoolpump": 1,
        "temp_in": 573,
        "temp_out": 773,
        "pressure": 8000000,
        "pdrop": 20088.23,
        "mf": 956.3,
        "primary_coolant_switch": "Helium",
        "coolant_density": 5.64,
        "label": "First Wall and Blanket",
    }

    assert (
        pytest.approx(blanket_library_fixture.pumppower(False, **data), rel=1e-4)
        == 3.2374845432302464
    )


class ComponentHalfHeightParam(NamedTuple):
    hmax: Any = None
    dz_xpoint_divertor: Any = None
    dz_shld_vv_gap: Any = None
    dz_blkt_upper: Any = None
    dz_shld_upper: Any = None
    dr_fw_plasma_gap_inboard: Any = None
    dr_fw_plasma_gap_outboard: Any = None
    dr_fw_inboard: Any = None
    dr_fw_outboard: Any = None
    dz_vv_lower: Any = None
    dz_vv_upper: Any = None
    z_plasma_xpoint_lower: Any = None
    z_plasma_xpoint_upper: Any = None
    n_divertors: Any = None
    dz_divertor: Any = None
    icomponent: Any = None
    expected_icomponent: Any = None
    expected_half_height: Any = None


@pytest.mark.parametrize(
    "componenthalfheightparam",
    (
        ComponentHalfHeightParam(
            hmax=8.8182171641274945,
            dz_xpoint_divertor=2.0018838307941582,
            dz_shld_vv_gap=0.16300000000000001,
            dz_blkt_upper=0.85000000000000009,
            dz_shld_upper=0.59999999999999998,
            dr_fw_plasma_gap_inboard=0.25,
            dr_fw_plasma_gap_outboard=0.25,
            dr_fw_inboard=0.018000000000000002,
            dr_fw_outboard=0.018000000000000002,
            dz_vv_lower=0.30000000000000004,
            dz_vv_upper=0.30000000000000004,
            z_plasma_xpoint_lower=4.93333333333333333,
            z_plasma_xpoint_upper=4.93333333333333333,
            n_divertors=1,
            dz_divertor=0.62000000000000011,
            icomponent=0,
            expected_icomponent=0,
            expected_half_height=5.9532752487304119,
        ),
    ),
)
def test_component_half_height(
    componenthalfheightparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for component_half_height.

    This test was generated using data from tests/regression/input_files/large_tokamak.IN.DAT.

    :param componenthalfheightparam: the data used to mock and assert in this test.
    :type componenthalfheightparam: componenthalfheightparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(build_variables, "hmax", componenthalfheightparam.hmax)
    monkeypatch.setattr(
        build_variables,
        "dz_xpoint_divertor",
        componenthalfheightparam.dz_xpoint_divertor,
    )
    monkeypatch.setattr(
        build_variables,
        "dz_shld_vv_gap",
        componenthalfheightparam.dz_shld_vv_gap,
    )
    monkeypatch.setattr(
        build_variables, "dz_blkt_upper", componenthalfheightparam.dz_blkt_upper
    )
    monkeypatch.setattr(
        build_variables, "dz_shld_upper", componenthalfheightparam.dz_shld_upper
    )
    monkeypatch.setattr(
        build_variables,
        "dr_fw_plasma_gap_inboard",
        componenthalfheightparam.dr_fw_plasma_gap_inboard,
    )
    monkeypatch.setattr(
        build_variables,
        "dr_fw_plasma_gap_outboard",
        componenthalfheightparam.dr_fw_plasma_gap_outboard,
    )
    monkeypatch.setattr(
        build_variables, "dr_fw_inboard", componenthalfheightparam.dr_fw_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_fw_outboard", componenthalfheightparam.dr_fw_outboard
    )
    monkeypatch.setattr(
        build_variables, "dz_vv_lower", componenthalfheightparam.dz_vv_lower
    )
    monkeypatch.setattr(
        build_variables, "dz_vv_upper", componenthalfheightparam.dz_vv_upper
    )
    monkeypatch.setattr(
        build_variables,
        "z_plasma_xpoint_lower",
        componenthalfheightparam.z_plasma_xpoint_lower,
    )
    monkeypatch.setattr(
        build_variables,
        "z_plasma_xpoint_upper",
        componenthalfheightparam.z_plasma_xpoint_upper,
    )
    monkeypatch.setattr(
        physics_variables, "n_divertors", componenthalfheightparam.n_divertors
    )
    monkeypatch.setattr(
        divertor_variables, "dz_divertor", componenthalfheightparam.dz_divertor
    )

    half_height = blanket_library_fixture.component_half_height(
        componenthalfheightparam.icomponent
    )

    assert half_height == pytest.approx(componenthalfheightparam.expected_half_height)


class DshapedComponentParam(NamedTuple):
    rsldi: Any = None
    dr_shld_inboard: Any = None
    dr_blkt_inboard: Any = None
    dr_fw_inboard: Any = None
    dr_fw_plasma_gap_inboard: Any = None
    dr_fw_plasma_gap_outboard: Any = None
    dr_fw_outboard: Any = None
    blareaib: Any = None
    blareaob: Any = None
    blarea: Any = None
    dr_blkt_outboard: Any = None
    dz_blkt_upper: Any = None
    shareaib: Any = None
    shareaob: Any = None
    sharea: Any = None
    dr_shld_outboard: Any = None
    dz_shld_upper: Any = None
    rsldo: Any = None
    dr_vv_inboard: Any = None
    dr_vv_outboard: Any = None
    dz_vv_upper: Any = None
    dz_vv_lower: Any = None
    vol_blkt_inboard: Any = None
    vol_blkt_outboard: Any = None
    vol_blkt_total: Any = None
    volshld: Any = None
    vol_vv: Any = None
    rminor: Any = None
    volshldi: Any = None
    volshldo: Any = None
    vol_vv_inboard: Any = None
    vol_vv_outboard: Any = None
    hblnkt: Any = None
    hshld: Any = None
    hvv: Any = None
    icomponent: Any = None
    expected_blareaib: Any = None
    expected_blareaob: Any = None
    expected_blarea: Any = None
    expected_shareaib: Any = None
    expected_shareaob: Any = None
    expected_sharea: Any = None
    expected_vol_blkt_outboard: Any = None
    expected_volblkt: Any = None
    expected_volshld: Any = None
    expected_vol_vv: Any = None
    expected_volshldi: Any = None
    expected_volshldo: Any = None
    expected_vol_vv_inboard: Any = None
    expected_vol_vv_outboard: Any = None
    expected_icomponent: Any = None


@pytest.mark.parametrize(
    "dshapedcomponentparam",
    (
        DshapedComponentParam(
            rsldi=1.5,
            dr_shld_inboard=0.40000000000000002,
            dr_blkt_inboard=0,
            dr_fw_inboard=0.018000000000000002,
            dr_fw_plasma_gap_inboard=0.10000000000000001,
            dr_fw_plasma_gap_outboard=0.10000000000000001,
            dr_fw_outboard=0.018000000000000002,
            blareaib=0,
            blareaob=0,
            blarea=0,
            dr_blkt_outboard=1,
            dz_blkt_upper=0.5,
            shareaib=0,
            shareaob=0,
            sharea=0,
            dr_shld_outboard=0.30000000000000004,
            dz_shld_upper=0.60000000000000009,
            rsldo=8.4000000000000004,
            dr_vv_inboard=0.20000000000000001,
            dr_vv_outboard=0.30000000000000004,
            dz_vv_upper=0.30000000000000004,
            dz_vv_lower=0.30000000000000004,
            vol_blkt_inboard=0,
            vol_blkt_outboard=0,
            vol_blkt_total=0,
            volshld=0,
            vol_vv=0,
            rminor=2.5,
            volshldi=0,
            volshldo=0,
            vol_vv_inboard=0,
            vol_vv_outboard=0,
            hblnkt=8.25,
            hshld=8.75,
            hvv=9.4349999999999987,
            icomponent=0,
            expected_blareaib=196.97785938008002,
            expected_blareaob=852.24160940262459,
            expected_blarea=1049.2194687827046,
            expected_shareaib=0,
            expected_shareaob=0,
            expected_sharea=0,
            expected_vol_blkt_outboard=691.06561956756764,
            expected_volblkt=691.06561956756764,
            expected_volshld=0,
            expected_vol_vv=0,
            expected_volshldi=0,
            expected_volshldo=0,
            expected_vol_vv_inboard=0,
            expected_vol_vv_outboard=0,
            expected_icomponent=0,
        ),
        DshapedComponentParam(
            rsldi=1.5,
            dr_shld_inboard=0.40000000000000002,
            dr_blkt_inboard=0,
            dr_fw_inboard=0.018000000000000002,
            dr_fw_plasma_gap_inboard=0.10000000000000001,
            dr_fw_plasma_gap_outboard=0.10000000000000001,
            dr_fw_outboard=0.018000000000000002,
            blareaib=196.97785938008002,
            blareaob=852.24160940262459,
            blarea=1049.2194687827046,
            dr_blkt_outboard=1,
            dz_blkt_upper=0.5,
            shareaib=0,
            shareaob=0,
            sharea=0,
            dr_shld_outboard=0.30000000000000004,
            dz_shld_upper=0.60000000000000009,
            rsldo=8.4000000000000004,
            dr_vv_inboard=0.20000000000000001,
            dr_vv_outboard=0.30000000000000004,
            dz_vv_upper=0.30000000000000004,
            dz_vv_lower=0.30000000000000004,
            vol_blkt_inboard=0,
            vol_blkt_outboard=691.06561956756764,
            vol_blkt_total=691.06561956756764,
            volshld=0,
            vol_vv=0,
            rminor=2.5,
            volshldi=0,
            volshldo=0,
            vol_vv_inboard=0,
            vol_vv_outboard=0,
            hblnkt=8.25,
            hshld=8.75,
            hvv=9.4349999999999987,
            icomponent=0,
            expected_blareaib=196.97785938008002,
            expected_blareaob=852.24160940262459,
            expected_blarea=1049.2194687827046,
            expected_shareaib=208.91591146372122,
            expected_shareaob=1013.8483589087293,
            expected_sharea=1222.7642703724505,
            expected_vol_blkt_outboard=691.06561956756764,
            expected_volblkt=691.06561956756764,
            expected_volshld=450.46122947809488,
            expected_vol_vv=0,
            expected_volshldi=79.896984366095609,
            expected_volshldo=370.5642451119993,
            expected_vol_vv_inboard=0,
            expected_vol_vv_outboard=0,
            expected_icomponent=1,
        ),
        DshapedComponentParam(
            rsldi=1.5,
            dr_shld_inboard=0.40000000000000002,
            dr_blkt_inboard=0,
            dr_fw_inboard=0.018000000000000002,
            dr_fw_plasma_gap_inboard=0.10000000000000001,
            dr_fw_plasma_gap_outboard=0.10000000000000001,
            dr_fw_outboard=0.018000000000000002,
            blareaib=196.97785938008002,
            blareaob=852.24160940262459,
            blarea=1049.2194687827046,
            dr_blkt_outboard=1,
            dz_blkt_upper=0.5,
            shareaib=208.91591146372122,
            shareaob=1013.8483589087293,
            sharea=1222.7642703724505,
            dr_shld_outboard=0.30000000000000004,
            dz_shld_upper=0.60000000000000009,
            rsldo=8.4000000000000004,
            dr_vv_inboard=0.20000000000000001,
            dr_vv_outboard=0.30000000000000004,
            dz_vv_upper=0.30000000000000004,
            dz_vv_lower=0.30000000000000004,
            vol_blkt_inboard=0,
            vol_blkt_outboard=691.06561956756764,
            vol_blkt_total=691.06561956756764,
            volshld=450.46122947809488,
            vol_vv=0,
            rminor=2.5,
            volshldi=79.896984366095609,
            volshldo=370.5642451119993,
            vol_vv_inboard=0,
            vol_vv_outboard=0,
            hblnkt=8.25,
            hshld=8.75,
            hvv=9.4349999999999987,
            icomponent=1,
            expected_blareaib=196.97785938008002,
            expected_blareaob=852.24160940262459,
            expected_blarea=1049.2194687827046,
            expected_shareaib=208.91591146372122,
            expected_shareaob=1013.8483589087293,
            expected_sharea=1222.7642703724505,
            expected_vol_blkt_outboard=691.06561956756764,
            expected_volblkt=691.06561956756764,
            expected_volshld=450.46122947809488,
            expected_vol_vv=340.45369594344834,
            expected_volshldi=79.896984366095609,
            expected_volshldo=370.5642451119993,
            expected_vol_vv_inboard=34.253413020620215,
            expected_vol_vv_outboard=306.20028292282814,
            expected_icomponent=2,
        ),
    ),
)
def test_dshaped_component(dshapedcomponentparam, monkeypatch, blanket_library_fixture):
    """
    Automatically generated Regression Unit Test for dshaped_component.

    This test was generated using data from tests/regression/input_files/st_regression.IN.DAT.

    :param dshapedcomponentparam: the data used to mock and assert in this test.
    :type dshapedcomponentparam: dshapedcomponentparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(build_variables, "rsldi", dshapedcomponentparam.rsldi)
    monkeypatch.setattr(
        build_variables, "dr_shld_inboard", dshapedcomponentparam.dr_shld_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_blkt_inboard", dshapedcomponentparam.dr_blkt_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_fw_inboard", dshapedcomponentparam.dr_fw_inboard
    )
    monkeypatch.setattr(
        build_variables,
        "dr_fw_plasma_gap_inboard",
        dshapedcomponentparam.dr_fw_plasma_gap_inboard,
    )
    monkeypatch.setattr(
        build_variables,
        "dr_fw_plasma_gap_outboard",
        dshapedcomponentparam.dr_fw_plasma_gap_outboard,
    )
    monkeypatch.setattr(
        build_variables, "dr_fw_outboard", dshapedcomponentparam.dr_fw_outboard
    )
    monkeypatch.setattr(build_variables, "blareaib", dshapedcomponentparam.blareaib)
    monkeypatch.setattr(build_variables, "blareaob", dshapedcomponentparam.blareaob)
    monkeypatch.setattr(build_variables, "blarea", dshapedcomponentparam.blarea)
    monkeypatch.setattr(
        build_variables, "dr_blkt_outboard", dshapedcomponentparam.dr_blkt_outboard
    )
    monkeypatch.setattr(
        build_variables, "dz_blkt_upper", dshapedcomponentparam.dz_blkt_upper
    )
    monkeypatch.setattr(build_variables, "shareaib", dshapedcomponentparam.shareaib)
    monkeypatch.setattr(build_variables, "shareaob", dshapedcomponentparam.shareaob)
    monkeypatch.setattr(build_variables, "sharea", dshapedcomponentparam.sharea)
    monkeypatch.setattr(
        build_variables, "dr_shld_outboard", dshapedcomponentparam.dr_shld_outboard
    )
    monkeypatch.setattr(
        build_variables, "dz_shld_upper", dshapedcomponentparam.dz_shld_upper
    )
    monkeypatch.setattr(build_variables, "rsldo", dshapedcomponentparam.rsldo)
    monkeypatch.setattr(
        build_variables, "dr_vv_inboard", dshapedcomponentparam.dr_vv_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_vv_outboard", dshapedcomponentparam.dr_vv_outboard
    )
    monkeypatch.setattr(
        build_variables, "dz_vv_upper", dshapedcomponentparam.dz_vv_upper
    )
    monkeypatch.setattr(
        build_variables, "dz_vv_lower", dshapedcomponentparam.dz_vv_lower
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_inboard", dshapedcomponentparam.vol_blkt_inboard
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_outboard", dshapedcomponentparam.vol_blkt_outboard
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_total", dshapedcomponentparam.vol_blkt_total
    )
    monkeypatch.setattr(fwbs_variables, "volshld", dshapedcomponentparam.volshld)
    monkeypatch.setattr(fwbs_variables, "vol_vv", dshapedcomponentparam.vol_vv)
    monkeypatch.setattr(physics_variables, "rminor", dshapedcomponentparam.rminor)
    monkeypatch.setattr(blanket_library, "volshldi", dshapedcomponentparam.volshldi)
    monkeypatch.setattr(blanket_library, "volshldo", dshapedcomponentparam.volshldo)
    monkeypatch.setattr(
        blanket_library, "vol_vv_inboard", dshapedcomponentparam.vol_vv_inboard
    )
    monkeypatch.setattr(
        blanket_library, "vol_vv_outboard", dshapedcomponentparam.vol_vv_outboard
    )
    monkeypatch.setattr(blanket_library, "hblnkt", dshapedcomponentparam.hblnkt)
    monkeypatch.setattr(blanket_library, "hshld", dshapedcomponentparam.hshld)
    monkeypatch.setattr(blanket_library, "hvv", dshapedcomponentparam.hvv)

    blanket_library_fixture.dshaped_component(dshapedcomponentparam.icomponent)

    assert build_variables.blareaib == pytest.approx(
        dshapedcomponentparam.expected_blareaib
    )
    assert build_variables.blareaob == pytest.approx(
        dshapedcomponentparam.expected_blareaob
    )
    assert build_variables.blarea == pytest.approx(
        dshapedcomponentparam.expected_blarea
    )
    assert fwbs_variables.vol_blkt_outboard == pytest.approx(
        dshapedcomponentparam.expected_vol_blkt_outboard
    )
    assert fwbs_variables.vol_blkt_total == pytest.approx(
        dshapedcomponentparam.expected_volblkt
    )


class EllipticalComponentParam(NamedTuple):
    rsldi: Any = None
    dr_shld_inboard: Any = None
    dr_blkt_inboard: Any = None
    rsldo: Any = None
    dr_shld_outboard: Any = None
    dr_blkt_outboard: Any = None
    blareaib: Any = None
    blareaob: Any = None
    blarea: Any = None
    dz_blkt_upper: Any = None
    shareaib: Any = None
    shareaob: Any = None
    sharea: Any = None
    dz_shld_upper: Any = None
    dr_vv_inboard: Any = None
    dr_vv_outboard: Any = None
    dz_vv_upper: Any = None
    dz_vv_lower: Any = None
    vol_blkt_inboard: Any = None
    vol_blkt_outboard: Any = None
    vol_blkt_total: Any = None
    volshld: Any = None
    vol_vv: Any = None
    rmajor: Any = None
    rminor: Any = None
    triang: Any = None
    volshldi: Any = None
    volshldo: Any = None
    vol_vv_inboard: Any = None
    vol_vv_outboard: Any = None
    hblnkt: Any = None
    hshld: Any = None
    hvv: Any = None
    icomponent: Any = None
    expected_blareaib: Any = None
    expected_blareaob: Any = None
    expected_blarea: Any = None
    expected_shareaib: Any = None
    expected_shareaob: Any = None
    expected_sharea: Any = None
    expected_vol_blkt_inboard: Any = None
    expected_vol_blkt_outboard: Any = None
    expected_volblkt: Any = None
    expected_volshld: Any = None
    expected_vol_vv: Any = None
    expected_volshldi: Any = None
    expected_volshldo: Any = None
    expected_vol_vv_inboard: Any = None
    expected_vol_vv_outboard: Any = None
    expected_icomponent: Any = None


@pytest.mark.parametrize(
    "ellipticalcomponentparam",
    (
        EllipticalComponentParam(
            rsldi=4.0833333333333339,
            dr_shld_inboard=0.30000000000000004,
            dr_blkt_inboard=0.70000000000000007,
            rsldo=12.716666666666667,
            dr_shld_outboard=0.80000000000000004,
            dr_blkt_outboard=1,
            blareaib=0,
            blareaob=0,
            blarea=0,
            dz_blkt_upper=0.85000000000000009,
            shareaib=0,
            shareaob=0,
            sharea=0,
            dz_shld_upper=0.59999999999999998,
            dr_vv_inboard=0.30000000000000004,
            dr_vv_outboard=0.30000000000000004,
            dz_vv_upper=0.30000000000000004,
            dz_vv_lower=0.30000000000000004,
            vol_blkt_inboard=0,
            vol_blkt_outboard=0,
            vol_blkt_total=0,
            volshld=0,
            vol_vv=0,
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            volshldi=0,
            volshldo=0,
            vol_vv_inboard=0,
            vol_vv_outboard=0,
            hblnkt=5.9532752487304119,
            hshld=6.8032752487304133,
            hvv=7.5032752487304135,
            icomponent=0,
            expected_blareaib=664.9687712975541,
            expected_blareaob=1101.3666396424403,
            expected_blarea=1766.3354109399943,
            expected_shareaib=0,
            expected_shareaob=0,
            expected_sharea=0,
            expected_vol_blkt_inboard=315.83946385183026,
            expected_vol_blkt_outboard=1020.3677420460117,
            expected_volblkt=1336.207205897842,
            expected_volshld=0,
            expected_vol_vv=0,
            expected_volshldi=0,
            expected_volshldo=0,
            expected_vol_vv_inboard=0,
            expected_vol_vv_outboard=0,
            expected_icomponent=0,
        ),
        EllipticalComponentParam(
            rsldi=4.0833333333333339,
            dr_shld_inboard=0.30000000000000004,
            dr_blkt_inboard=0.70000000000000007,
            rsldo=12.716666666666667,
            dr_shld_outboard=0.80000000000000004,
            dr_blkt_outboard=1,
            blareaib=664.9687712975541,
            blareaob=1101.3666396424403,
            blarea=1766.3354109399943,
            dz_blkt_upper=0.85000000000000009,
            shareaib=0,
            shareaob=0,
            sharea=0,
            dz_shld_upper=0.59999999999999998,
            dr_vv_inboard=0.30000000000000004,
            dr_vv_outboard=0.30000000000000004,
            dz_vv_upper=0.30000000000000004,
            dz_vv_lower=0.30000000000000004,
            vol_blkt_inboard=315.83946385183026,
            vol_blkt_outboard=1020.3677420460117,
            vol_blkt_total=1336.207205897842,
            volshld=0,
            vol_vv=0,
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            volshldi=0,
            volshldo=0,
            vol_vv_inboard=0,
            vol_vv_outboard=0,
            hblnkt=5.9532752487304119,
            hshld=6.8032752487304133,
            hvv=7.5032752487304135,
            icomponent=1,
            expected_blareaib=664.9687712975541,
            expected_blareaob=1101.3666396424403,
            expected_blarea=1766.3354109399943,
            expected_shareaib=700.06731267447844,
            expected_shareaob=1344.1106481995357,
            expected_sharea=2044.1779608740142,
            expected_vol_blkt_inboard=315.83946385183026,
            expected_vol_blkt_outboard=1020.3677420460117,
            expected_volblkt=1336.207205897842,
            expected_volshld=1124.4621612595051,
            expected_vol_vv=0,
            expected_volshldi=177.89822933168091,
            expected_volshldo=946.56393192782434,
            expected_vol_vv_inboard=0,
            expected_vol_vv_outboard=0,
            expected_icomponent=1,
        ),
        EllipticalComponentParam(
            rsldi=4.0833333333333339,
            dr_shld_inboard=0.30000000000000004,
            dr_blkt_inboard=0.70000000000000007,
            rsldo=12.716666666666667,
            dr_shld_outboard=0.80000000000000004,
            dr_blkt_outboard=1,
            blareaib=664.9687712975541,
            blareaob=1101.3666396424403,
            blarea=1766.3354109399943,
            dz_blkt_upper=0.85000000000000009,
            shareaib=700.06731267447844,
            shareaob=1344.1106481995357,
            sharea=2044.1779608740142,
            dz_shld_upper=0.59999999999999998,
            dr_vv_inboard=0.30000000000000004,
            dr_vv_outboard=0.30000000000000004,
            dz_vv_upper=0.30000000000000004,
            dz_vv_lower=0.30000000000000004,
            vol_blkt_inboard=315.83946385183026,
            vol_blkt_outboard=1020.3677420460117,
            vol_blkt_total=1336.207205897842,
            volshld=1124.4621612595051,
            vol_vv=0,
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            volshldi=177.89822933168091,
            volshldo=946.56393192782434,
            vol_vv_inboard=0,
            vol_vv_outboard=0,
            hblnkt=5.9532752487304119,
            hshld=6.8032752487304133,
            hvv=7.5032752487304135,
            icomponent=2,
            expected_blareaib=664.9687712975541,
            expected_blareaob=1101.3666396424403,
            expected_blarea=1766.3354109399943,
            expected_shareaib=700.06731267447844,
            expected_shareaob=1344.1106481995357,
            expected_sharea=2044.1779608740142,
            expected_vol_blkt_inboard=315.83946385183026,
            expected_vol_blkt_outboard=1020.3677420460117,
            expected_volblkt=1336.207205897842,
            expected_volshld=1124.4621612595051,
            expected_vol_vv=584.07334775041659,
            expected_volshldi=177.89822933168091,
            expected_volshldo=946.56393192782434,
            expected_vol_vv_inboard=143.03162449152501,
            expected_vol_vv_outboard=441.04172325889158,
            expected_icomponent=2,
        ),
    ),
)
def test_elliptical_component(
    ellipticalcomponentparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for elliptical_component.

    This test was generated using data from tests/regression/input_files/large_tokamak_once_through.IN.DAT.

    :param ellipticalcomponentparam: the data used to mock and assert in this test.
    :type ellipticalcomponentparam: ellipticalcomponentparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(build_variables, "rsldi", ellipticalcomponentparam.rsldi)
    monkeypatch.setattr(
        build_variables, "dr_shld_inboard", ellipticalcomponentparam.dr_shld_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_blkt_inboard", ellipticalcomponentparam.dr_blkt_inboard
    )
    monkeypatch.setattr(build_variables, "rsldo", ellipticalcomponentparam.rsldo)
    monkeypatch.setattr(
        build_variables, "dr_shld_outboard", ellipticalcomponentparam.dr_shld_outboard
    )
    monkeypatch.setattr(
        build_variables, "dr_blkt_outboard", ellipticalcomponentparam.dr_blkt_outboard
    )
    monkeypatch.setattr(build_variables, "blareaib", ellipticalcomponentparam.blareaib)
    monkeypatch.setattr(build_variables, "blareaob", ellipticalcomponentparam.blareaob)
    monkeypatch.setattr(build_variables, "blarea", ellipticalcomponentparam.blarea)
    monkeypatch.setattr(
        build_variables, "dz_blkt_upper", ellipticalcomponentparam.dz_blkt_upper
    )
    monkeypatch.setattr(build_variables, "shareaib", ellipticalcomponentparam.shareaib)
    monkeypatch.setattr(build_variables, "shareaob", ellipticalcomponentparam.shareaob)
    monkeypatch.setattr(build_variables, "sharea", ellipticalcomponentparam.sharea)
    monkeypatch.setattr(
        build_variables, "dz_shld_upper", ellipticalcomponentparam.dz_shld_upper
    )
    monkeypatch.setattr(
        build_variables, "dr_vv_inboard", ellipticalcomponentparam.dr_vv_inboard
    )
    monkeypatch.setattr(
        build_variables, "dr_vv_outboard", ellipticalcomponentparam.dr_vv_outboard
    )
    monkeypatch.setattr(
        build_variables, "dz_vv_upper", ellipticalcomponentparam.dz_vv_upper
    )
    monkeypatch.setattr(
        build_variables, "dz_vv_lower", ellipticalcomponentparam.dz_vv_lower
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_inboard", ellipticalcomponentparam.vol_blkt_inboard
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_outboard", ellipticalcomponentparam.vol_blkt_outboard
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_total", ellipticalcomponentparam.vol_blkt_total
    )
    monkeypatch.setattr(fwbs_variables, "volshld", ellipticalcomponentparam.volshld)
    monkeypatch.setattr(fwbs_variables, "vol_vv", ellipticalcomponentparam.vol_vv)
    monkeypatch.setattr(physics_variables, "rmajor", ellipticalcomponentparam.rmajor)
    monkeypatch.setattr(physics_variables, "rminor", ellipticalcomponentparam.rminor)
    monkeypatch.setattr(physics_variables, "triang", ellipticalcomponentparam.triang)
    monkeypatch.setattr(blanket_library, "volshldi", ellipticalcomponentparam.volshldi)
    monkeypatch.setattr(blanket_library, "volshldo", ellipticalcomponentparam.volshldo)
    monkeypatch.setattr(
        blanket_library, "vol_vv_inboard", ellipticalcomponentparam.vol_vv_inboard
    )
    monkeypatch.setattr(
        blanket_library, "vol_vv_outboard", ellipticalcomponentparam.vol_vv_outboard
    )
    monkeypatch.setattr(blanket_library, "hblnkt", ellipticalcomponentparam.hblnkt)
    monkeypatch.setattr(blanket_library, "hshld", ellipticalcomponentparam.hshld)
    monkeypatch.setattr(blanket_library, "hvv", ellipticalcomponentparam.hvv)

    blanket_library_fixture.elliptical_component(ellipticalcomponentparam.icomponent)

    assert build_variables.blareaib == pytest.approx(
        ellipticalcomponentparam.expected_blareaib
    )
    assert build_variables.blareaob == pytest.approx(
        ellipticalcomponentparam.expected_blareaob
    )
    assert build_variables.blarea == pytest.approx(
        ellipticalcomponentparam.expected_blarea
    )
    assert build_variables.shareaib == pytest.approx(
        ellipticalcomponentparam.expected_shareaib
    )
    assert build_variables.shareaob == pytest.approx(
        ellipticalcomponentparam.expected_shareaob
    )
    assert build_variables.sharea == pytest.approx(
        ellipticalcomponentparam.expected_sharea
    )
    assert fwbs_variables.vol_blkt_inboard == pytest.approx(
        ellipticalcomponentparam.expected_vol_blkt_inboard
    )
    assert fwbs_variables.vol_blkt_outboard == pytest.approx(
        ellipticalcomponentparam.expected_vol_blkt_outboard
    )
    assert fwbs_variables.vol_blkt_total == pytest.approx(
        ellipticalcomponentparam.expected_volblkt
    )
    assert fwbs_variables.volshld == pytest.approx(
        ellipticalcomponentparam.expected_volshld
    )
    assert fwbs_variables.vol_vv == pytest.approx(
        ellipticalcomponentparam.expected_vol_vv
    )
    assert blanket_library.volshldi == pytest.approx(
        ellipticalcomponentparam.expected_volshldi
    )
    assert blanket_library.volshldo == pytest.approx(
        ellipticalcomponentparam.expected_volshldo
    )
    assert blanket_library.vol_vv_inboard == pytest.approx(
        ellipticalcomponentparam.expected_vol_vv_inboard
    )
    assert blanket_library.vol_vv_outboard == pytest.approx(
        ellipticalcomponentparam.expected_vol_vv_outboard
    )


class ApplyCoverageFactorsParam(NamedTuple):
    blareaob: Any = None
    blarea: Any = None
    blareaib: Any = None
    shareaib: Any = None
    shareaob: Any = None
    sharea: Any = None
    f_ster_div_single: Any = None
    f_a_fw_hcd: Any = None
    vol_blkt_outboard: Any = None
    vol_blkt_inboard: Any = None
    vol_blkt_total: Any = None
    fvolsi: Any = None
    fvolso: Any = None
    volshld: Any = None
    vol_vv: Any = None
    fvoldw: Any = None
    n_divertors: Any = None
    volshldi: Any = None
    volshldo: Any = None
    expected_blareaob: Any = None
    expected_blarea: Any = None
    expected_shareaob: Any = None
    expected_sharea: Any = None
    expected_vol_blkt_outboard: Any = None
    expected_volblkt: Any = None
    expected_volshld: Any = None
    expected_vol_vv: Any = None
    expected_volshldo: Any = None


@pytest.mark.parametrize(
    "applycoveragefactorsparam",
    (
        ApplyCoverageFactorsParam(
            blareaob=1101.3666396424403,
            blarea=1766.3354109399943,
            blareaib=664.9687712975541,
            shareaib=700.06731267447844,
            shareaob=1344.1106481995357,
            sharea=2044.1779608740142,
            f_ster_div_single=0.115,
            f_a_fw_hcd=0,
            vol_blkt_outboard=1020.3677420460117,
            vol_blkt_inboard=315.83946385183026,
            vol_blkt_total=1336.207205897842,
            fvolsi=1,
            fvolso=0.64000000000000001,
            volshld=1124.4621612595051,
            vol_vv=584.07334775041659,
            fvoldw=1.74,
            n_divertors=1,
            volshldi=177.89822933168091,
            volshldo=946.56393192782434,
            expected_blareaob=898.23806738434075,
            expected_blarea=1563.2068386818949,
            expected_shareaob=860.23081484770285,
            expected_sharea=1560.2981275221814,
            expected_vol_blkt_outboard=866.70391336775992,
            expected_volblkt=1182.5433772195902,
            expected_volshld=783.69914576548854,
            expected_vol_vv=1016.2876250857248,
            expected_volshldo=605.80091643380763,
        ),
    ),
)
def test_apply_coverage_factors(
    applycoveragefactorsparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for apply_coverage_factors.

    This test was generated using data from tests/regression/input_files/large_tokamak_once_through.IN.DAT.

    :param applycoveragefactorsparam: the data used to mock and assert in this test.
    :type applycoveragefactorsparam: applycoveragefactorsparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(build_variables, "blareaob", applycoveragefactorsparam.blareaob)
    monkeypatch.setattr(build_variables, "blarea", applycoveragefactorsparam.blarea)
    monkeypatch.setattr(build_variables, "blareaib", applycoveragefactorsparam.blareaib)
    monkeypatch.setattr(build_variables, "shareaib", applycoveragefactorsparam.shareaib)
    monkeypatch.setattr(build_variables, "shareaob", applycoveragefactorsparam.shareaob)
    monkeypatch.setattr(build_variables, "sharea", applycoveragefactorsparam.sharea)
    monkeypatch.setattr(
        fwbs_variables, "f_ster_div_single", applycoveragefactorsparam.f_ster_div_single
    )
    monkeypatch.setattr(
        fwbs_variables, "f_a_fw_hcd", applycoveragefactorsparam.f_a_fw_hcd
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_outboard", applycoveragefactorsparam.vol_blkt_outboard
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_inboard", applycoveragefactorsparam.vol_blkt_inboard
    )
    monkeypatch.setattr(
        fwbs_variables, "vol_blkt_total", applycoveragefactorsparam.vol_blkt_total
    )
    monkeypatch.setattr(fwbs_variables, "fvolsi", applycoveragefactorsparam.fvolsi)
    monkeypatch.setattr(fwbs_variables, "fvolso", applycoveragefactorsparam.fvolso)
    monkeypatch.setattr(fwbs_variables, "volshld", applycoveragefactorsparam.volshld)
    monkeypatch.setattr(fwbs_variables, "vol_vv", applycoveragefactorsparam.vol_vv)
    monkeypatch.setattr(fwbs_variables, "fvoldw", applycoveragefactorsparam.fvoldw)
    monkeypatch.setattr(
        physics_variables, "n_divertors", applycoveragefactorsparam.n_divertors
    )
    monkeypatch.setattr(blanket_library, "volshldi", applycoveragefactorsparam.volshldi)
    monkeypatch.setattr(blanket_library, "volshldo", applycoveragefactorsparam.volshldo)

    blanket_library_fixture.apply_coverage_factors()

    assert build_variables.blareaob == pytest.approx(
        applycoveragefactorsparam.expected_blareaob
    )
    assert build_variables.blarea == pytest.approx(
        applycoveragefactorsparam.expected_blarea
    )
    assert build_variables.shareaob == pytest.approx(
        applycoveragefactorsparam.expected_shareaob
    )
    assert build_variables.sharea == pytest.approx(
        applycoveragefactorsparam.expected_sharea
    )
    assert fwbs_variables.vol_blkt_outboard == pytest.approx(
        applycoveragefactorsparam.expected_vol_blkt_outboard
    )
    assert fwbs_variables.vol_blkt_total == pytest.approx(
        applycoveragefactorsparam.expected_volblkt
    )
    assert fwbs_variables.volshld == pytest.approx(
        applycoveragefactorsparam.expected_volshld
    )
    assert fwbs_variables.vol_vv == pytest.approx(
        applycoveragefactorsparam.expected_vol_vv
    )
    assert blanket_library.volshldo == pytest.approx(
        applycoveragefactorsparam.expected_volshldo
    )


class BlanketModPolHeightParam(NamedTuple):
    dr_fw_plasma_gap_inboard: Any = None
    dr_fw_plasma_gap_outboard: Any = None
    i_fw_blkt_vv_shape: Any = None
    n_blkt_inboard_modules_poloidal: Any = None
    f_ster_div_single: Any = None
    n_blkt_outboard_modules_poloidal: Any = None
    itart: Any = None
    rminor: Any = None
    n_divertors: Any = None
    rmajor: Any = None
    triang: Any = None
    bllengi: Any = None
    bllengo: Any = None
    hblnkt: Any = None
    expected_bllengi: Any = None
    expected_bllengo: Any = None


@pytest.mark.parametrize(
    "blanketmodpolheightparam",
    (
        BlanketModPolHeightParam(
            dr_fw_plasma_gap_inboard=0.25,
            dr_fw_plasma_gap_outboard=0.25,
            i_fw_blkt_vv_shape=2,
            n_blkt_inboard_modules_poloidal=7,
            f_ster_div_single=0.115,
            n_blkt_outboard_modules_poloidal=8,
            itart=0,
            rminor=2.6666666666666665,
            n_divertors=1,
            rmajor=8,
            triang=0.5,
            bllengi=0,
            bllengo=0,
            hblnkt=5.9532752487304119,
            expected_bllengi=1.6252823720672551,
            expected_bllengo=1.7853902013340495,
        ),
        BlanketModPolHeightParam(
            dr_fw_plasma_gap_inboard=0.10000000000000001,
            dr_fw_plasma_gap_outboard=0.10000000000000001,
            i_fw_blkt_vv_shape=1,
            n_blkt_inboard_modules_poloidal=7,
            f_ster_div_single=0.115,
            n_blkt_outboard_modules_poloidal=8,
            itart=1,
            rminor=2.5,
            n_divertors=2,
            rmajor=4.5,
            triang=0.5,
            bllengi=0,
            bllengo=0,
            hblnkt=8.25,
            expected_bllengi=2.3571428571428572,
            expected_bllengo=2.0597205347177807,
        ),
    ),
)
def test_blanket_mod_pol_height(
    blanketmodpolheightparam,
    monkeypatch,
    blanket_library_fixture,
):
    """
    Automatically generated Regression Unit Test for blanket_mod_pol_height.

    This test was generated using data from blanket_files/large_tokamak_primary_pumping2.IN.DAT.

    :param blanketmodpolheightparam: the data used to mock and assert in this test.
    :type blanketmodpolheightparam: blanketmodpolheightparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        build_variables,
        "dr_fw_plasma_gap_inboard",
        blanketmodpolheightparam.dr_fw_plasma_gap_inboard,
    )
    monkeypatch.setattr(
        build_variables,
        "dr_fw_plasma_gap_outboard",
        blanketmodpolheightparam.dr_fw_plasma_gap_outboard,
    )
    monkeypatch.setattr(
        fwbs_variables,
        "i_fw_blkt_vv_shape",
        blanketmodpolheightparam.i_fw_blkt_vv_shape,
    )
    monkeypatch.setattr(
        fwbs_variables,
        "n_blkt_inboard_modules_poloidal",
        blanketmodpolheightparam.n_blkt_inboard_modules_poloidal,
    )
    monkeypatch.setattr(
        fwbs_variables, "f_ster_div_single", blanketmodpolheightparam.f_ster_div_single
    )
    monkeypatch.setattr(
        fwbs_variables,
        "n_blkt_outboard_modules_poloidal",
        blanketmodpolheightparam.n_blkt_outboard_modules_poloidal,
    )
    monkeypatch.setattr(physics_variables, "itart", blanketmodpolheightparam.itart)
    monkeypatch.setattr(physics_variables, "rminor", blanketmodpolheightparam.rminor)
    monkeypatch.setattr(
        physics_variables, "n_divertors", blanketmodpolheightparam.n_divertors
    )
    monkeypatch.setattr(physics_variables, "rmajor", blanketmodpolheightparam.rmajor)
    monkeypatch.setattr(physics_variables, "triang", blanketmodpolheightparam.triang)
    monkeypatch.setattr(blanket_library, "bllengi", blanketmodpolheightparam.bllengi)
    monkeypatch.setattr(blanket_library, "bllengo", blanketmodpolheightparam.bllengo)
    monkeypatch.setattr(blanket_library, "hblnkt", blanketmodpolheightparam.hblnkt)

    blanket_library_fixture.blanket_mod_pol_height()

    assert blanket_library.bllengi == pytest.approx(
        blanketmodpolheightparam.expected_bllengi
    )
    assert blanket_library.bllengo == pytest.approx(
        blanketmodpolheightparam.expected_bllengo
    )


class LiquidBreederPropertiesParam(NamedTuple):
    inlet_temp_liq: Any = None
    outlet_temp_liq: Any = None
    a_bz_liq: Any = None
    b_bz_liq: Any = None
    den_liq: Any = None
    specific_heat_liq: Any = None
    thermal_conductivity_liq: Any = None
    dynamic_viscosity_liq: Any = None
    electrical_conductivity_liq: Any = None
    i_blkt_liquid_breeder_type: Any = None
    hartmann_liq: Any = None
    b_mag_blkt: Any = None
    i_blkt_inboard: Any = None
    i_blkt_dual_coolant: Any = None
    bt: Any = None
    aspect: Any = None
    rmajor: Any = None
    dr_blkt_inboard: Any = None
    dr_blkt_outboard: Any = None
    ip: Any = None
    ofile: Any = None
    expected_den_liq: Any = None
    expected_specific_heat_liq: Any = None
    expected_thermal_conductivity_liq: Any = None
    expected_dynamic_viscosity_liq: Any = None
    expected_electrical_conductivity_liq: Any = None
    expected_hartmann_liq: Any = None
    expected_b_mag_blkt: Any = None


@pytest.mark.parametrize(
    "liquidbreederpropertiesparam",
    (
        LiquidBreederPropertiesParam(
            inlet_temp_liq=570,
            outlet_temp_liq=720,
            a_bz_liq=0.20000000000000001,
            b_bz_liq=0.20000000000000001,
            den_liq=9500,
            specific_heat_liq=190,
            thermal_conductivity_liq=30,
            dynamic_viscosity_liq=0,
            electrical_conductivity_liq=0,
            i_blkt_liquid_breeder_type=0,
            hartmann_liq=np.array(np.array((0, 0), order="F"), order="F").transpose(),
            b_mag_blkt=np.array(np.array((5, 5), order="F"), order="F").transpose(),
            i_blkt_inboard=1,
            i_blkt_dual_coolant=0,
            bt=5.7000000000000002,
            aspect=3,
            rmajor=8,
            dr_blkt_inboard=0.70000000000000007,
            dr_blkt_outboard=1,
            ip=0,
            ofile=11,
            expected_den_liq=9753.2497999999996,
            expected_specific_heat_liq=189.12018,
            expected_thermal_conductivity_liq=9.238260167312621,
            expected_dynamic_viscosity_liq=0.0017477589255667965,
            expected_electrical_conductivity_liq=861800.80431007256,
            expected_hartmann_liq=np.array(
                np.array((20319.245984309102, 9067.8426109080938), order="F"),
                order="F",
            ).transpose(),
            expected_b_mag_blkt=np.array(
                np.array((9.1505016722408019, 4.0835820895522392), order="F"),
                order="F",
            ).transpose(),
        ),
        LiquidBreederPropertiesParam(
            inlet_temp_liq=570,
            outlet_temp_liq=720,
            a_bz_liq=0.20000000000000001,
            b_bz_liq=0.20000000000000001,
            den_liq=9500,
            specific_heat_liq=190,
            thermal_conductivity_liq=30,
            expected_thermal_conductivity_liq=30,  # doesn't change when i_blkt_liquid_breeder_type=1
            dynamic_viscosity_liq=0,
            electrical_conductivity_liq=0,
            i_blkt_liquid_breeder_type=1,
            hartmann_liq=np.array(np.array((0, 0), order="F"), order="F").transpose(),
            b_mag_blkt=np.array(np.array((5, 5), order="F"), order="F").transpose(),
            i_blkt_inboard=1,
            i_blkt_dual_coolant=0,
            bt=5.7000000000000002,
            aspect=3,
            rmajor=8,
            dr_blkt_inboard=0.70000000000000007,
            dr_blkt_outboard=1,
            ip=0,
            ofile=11,
            expected_den_liq=305.30702851374997,
            expected_specific_heat_liq=34.640761200690406,
            expected_dynamic_viscosity_liq=0.00037298826343426359,
            expected_electrical_conductivity_liq=596562356750.5,
            expected_hartmann_liq=np.array(
                np.array((36595294.326740541, 16331332.841336451), order="F"),
                order="F",
            ).transpose(),
            expected_b_mag_blkt=np.array(
                np.array((9.1505016722408019, 4.0835820895522392), order="F"),
                order="F",
            ).transpose(),
        ),
    ),
)
def test_liquid_breeder_properties(
    liquidbreederpropertiesparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for liquid_breeder_properties.

    This test was generated using data from blanket_files/large_tokamak_dcll.IN.DAT.

    :param liquidbreederpropertiesparam: the data used to mock and assert in this test.
    :type liquidbreederpropertiesparam: liquidbreederpropertiesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        fwbs_variables, "inlet_temp_liq", liquidbreederpropertiesparam.inlet_temp_liq
    )
    monkeypatch.setattr(
        fwbs_variables, "outlet_temp_liq", liquidbreederpropertiesparam.outlet_temp_liq
    )
    monkeypatch.setattr(
        fwbs_variables, "a_bz_liq", liquidbreederpropertiesparam.a_bz_liq
    )
    monkeypatch.setattr(
        fwbs_variables, "b_bz_liq", liquidbreederpropertiesparam.b_bz_liq
    )
    monkeypatch.setattr(fwbs_variables, "den_liq", liquidbreederpropertiesparam.den_liq)
    monkeypatch.setattr(
        fwbs_variables,
        "specific_heat_liq",
        liquidbreederpropertiesparam.specific_heat_liq,
    )
    monkeypatch.setattr(
        fwbs_variables,
        "thermal_conductivity_liq",
        liquidbreederpropertiesparam.thermal_conductivity_liq,
    )
    monkeypatch.setattr(
        fwbs_variables,
        "dynamic_viscosity_liq",
        liquidbreederpropertiesparam.dynamic_viscosity_liq,
    )
    monkeypatch.setattr(
        fwbs_variables,
        "electrical_conductivity_liq",
        liquidbreederpropertiesparam.electrical_conductivity_liq,
    )
    monkeypatch.setattr(
        fwbs_variables,
        "i_blkt_liquid_breeder_type",
        liquidbreederpropertiesparam.i_blkt_liquid_breeder_type,
    )
    monkeypatch.setattr(
        fwbs_variables, "hartmann_liq", liquidbreederpropertiesparam.hartmann_liq
    )
    monkeypatch.setattr(
        fwbs_variables, "b_mag_blkt", liquidbreederpropertiesparam.b_mag_blkt
    )
    monkeypatch.setattr(
        fwbs_variables, "i_blkt_inboard", liquidbreederpropertiesparam.i_blkt_inboard
    )
    monkeypatch.setattr(
        fwbs_variables,
        "i_blkt_dual_coolant",
        liquidbreederpropertiesparam.i_blkt_dual_coolant,
    )
    monkeypatch.setattr(physics_variables, "bt", liquidbreederpropertiesparam.bt)
    monkeypatch.setattr(
        physics_variables, "aspect", liquidbreederpropertiesparam.aspect
    )
    monkeypatch.setattr(
        physics_variables, "rmajor", liquidbreederpropertiesparam.rmajor
    )
    monkeypatch.setattr(
        build_variables, "dr_blkt_inboard", liquidbreederpropertiesparam.dr_blkt_inboard
    )
    monkeypatch.setattr(
        build_variables,
        "dr_blkt_outboard",
        liquidbreederpropertiesparam.dr_blkt_outboard,
    )

    blanket_library_fixture.liquid_breeder_properties()

    assert fwbs_variables.den_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_den_liq
    )
    assert fwbs_variables.specific_heat_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_specific_heat_liq
    )
    assert fwbs_variables.thermal_conductivity_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_thermal_conductivity_liq
    )
    assert fwbs_variables.dynamic_viscosity_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_dynamic_viscosity_liq
    )
    assert fwbs_variables.electrical_conductivity_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_electrical_conductivity_liq
    )
    assert fwbs_variables.hartmann_liq == pytest.approx(
        liquidbreederpropertiesparam.expected_hartmann_liq
    )
    assert fwbs_variables.b_mag_blkt == pytest.approx(
        liquidbreederpropertiesparam.expected_b_mag_blkt
    )


class PressureDropParam(NamedTuple):
    radius_fw_channel: Any = None
    a_bz_liq: Any = None
    b_bz_liq: Any = None
    roughness_fw_channel: Any = None
    ip: Any = None
    ofile: Any = None
    i_ps: Any = None
    num_90: Any = None
    num_180: Any = None
    l_pipe: Any = None
    den: Any = None
    vsc: Any = None
    vv: Any = None
    label: Any = None
    expected_pressure_drop_out: Any = None


@pytest.mark.parametrize(
    "pressuredropparam",
    (
        PressureDropParam(
            radius_fw_channel=0.0060000000000000001,
            a_bz_liq=0.20000000000000001,
            b_bz_liq=0.20000000000000001,
            roughness_fw_channel=9.9999999999999995e-07,
            ip=0,
            ofile=11,
            i_ps=1,
            num_90=2,
            num_180=0,
            l_pipe=4,
            den=10.405276820718059,
            vsc=3.604452999475736e-05,
            vv=32.753134225223164,
            label="Inboard first wall",
            expected_pressure_drop_out=36214.869527556766,
        ),
    ),
)
def test_pressure_drop(pressuredropparam, monkeypatch, blanket_library_fixture):
    """
    Automatically generated Regression Unit Test for pressure_drop.

    This test was generated using data from blanket_files/large_tokamak_primary_pumping2.IN.DAT.

    :param pressuredropparam: the data used to mock and assert in this test.
    :type pressuredropparam: pressuredropparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        fwbs_variables, "radius_fw_channel", pressuredropparam.radius_fw_channel
    )
    monkeypatch.setattr(fwbs_variables, "a_bz_liq", pressuredropparam.a_bz_liq)
    monkeypatch.setattr(fwbs_variables, "b_bz_liq", pressuredropparam.b_bz_liq)
    monkeypatch.setattr(
        fwbs_variables, "roughness_fw_channel", pressuredropparam.roughness_fw_channel
    )

    pressure_drop_out = blanket_library_fixture.pressure_drop(
        i_ps=pressuredropparam.i_ps,
        num_90=pressuredropparam.num_90,
        num_180=pressuredropparam.num_180,
        l_pipe=pressuredropparam.l_pipe,
        den=pressuredropparam.den,
        vsc=pressuredropparam.vsc,
        vv=pressuredropparam.vv,
        label=pressuredropparam.label,
    )

    assert pressure_drop_out == pytest.approx(
        pressuredropparam.expected_pressure_drop_out
    )


class LiquidBreederPressureDropMhdParam(NamedTuple):
    i_blkt_liquid_breeder_channel_type: Any = None
    a_bz_liq: Any = None
    b_bz_liq: Any = None
    b_mag_blkt: Any = None
    bz_channel_conduct_liq: Any = None
    th_wall_secondary: Any = None
    ip: Any = None
    ofile: Any = None
    vel: Any = None
    vsc: Any = None
    conduct_liq: Any = None
    l_channel: Any = None
    num_pol: Any = None
    label: Any = None
    expected_liquid_breeder_pressure_drop_mhd_out: Any = None


@pytest.mark.parametrize(
    "liquidbreederpressuredropmhdparam",
    (
        LiquidBreederPressureDropMhdParam(
            i_blkt_liquid_breeder_channel_type=0,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=np.array(
                np.array((8.3930173480023953, 3.8678755427951117), order="F"),
                order="F",
            ).transpose(),
            bz_channel_conduct_liq=833000,
            th_wall_secondary=0.012500000000000001,
            ip=0,
            ofile=11,
            vel=0.061438753831945352,
            vsc=0.0017477589503313536,
            conduct_liq=861800.80431007256,
            l_channel=1.8937231989768815,
            num_pol=3072,
            label="Outboard blanket breeder liquid",
            expected_liquid_breeder_pressure_drop_mhd_out=282697824.60502106,
        ),
        LiquidBreederPressureDropMhdParam(
            i_blkt_liquid_breeder_channel_type=1,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=np.array(
                np.array((8.3930173480023953, 3.8678755427951117), order="F"),
                order="F",
            ).transpose(),
            bz_channel_conduct_liq=833000,
            th_wall_secondary=0.012500000000000001,
            ip=0,
            ofile=11,
            vel=0.061438753831945352,
            vsc=0.0017477589503313536,
            conduct_liq=861800.80431007256,
            l_channel=1.8937231989768815,
            num_pol=3072,
            label="Outboard blanket breeder liquid",
            expected_liquid_breeder_pressure_drop_mhd_out=160029.28473931071,
        ),
        LiquidBreederPressureDropMhdParam(
            i_blkt_liquid_breeder_channel_type=2,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=np.array(
                np.array((8.3930173480023953, 3.8678755427951117), order="F"),
                order="F",
            ).transpose(),
            bz_channel_conduct_liq=833000,
            th_wall_secondary=0.012500000000000001,
            ip=0,
            ofile=11,
            vel=0.061438753831945352,
            vsc=0.0017477589503313536,
            conduct_liq=861800.80431007256,
            l_channel=1.8937231989768815,
            num_pol=3072,
            label="Outboard blanket breeder liquid",
            expected_liquid_breeder_pressure_drop_mhd_out=282697824.60502106,
        ),
        LiquidBreederPressureDropMhdParam(
            i_blkt_liquid_breeder_channel_type=2,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=np.array(
                np.array((8.3930173480023953, 3.8678755427951117), order="F"),
                order="F",
            ).transpose(),
            bz_channel_conduct_liq=833000,
            th_wall_secondary=0.012500000000000001,
            ip=0,
            ofile=11,
            vel=0.042500391943592931,
            vsc=0.0017477589503313536,
            conduct_liq=861800.80431007256,
            l_channel=1.7176027768600395,
            num_pol=1792,
            label="Inboard blanket breeder liquid",
            expected_liquid_breeder_pressure_drop_mhd_out=487177576.30010164,
        ),
    ),
)
def test_liquid_breeder_pressure_drop_mhd(
    liquidbreederpressuredropmhdparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for liquid_breeder_pressure_drop_mhd.

    This test was generated using data from blanket_files/dcll_mms_lt_IN.DAT.

    :param liquidbreederpressuredropmhdparam: the data used to mock and assert in this test.
    :type liquidbreederpressuredropmhdparam: liquidbreederpressuredropmhdparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        fwbs_variables,
        "i_blkt_liquid_breeder_channel_type",
        liquidbreederpressuredropmhdparam.i_blkt_liquid_breeder_channel_type,
    )
    monkeypatch.setattr(
        fwbs_variables, "a_bz_liq", liquidbreederpressuredropmhdparam.a_bz_liq
    )
    monkeypatch.setattr(
        fwbs_variables, "b_bz_liq", liquidbreederpressuredropmhdparam.b_bz_liq
    )
    monkeypatch.setattr(
        fwbs_variables, "b_mag_blkt", liquidbreederpressuredropmhdparam.b_mag_blkt
    )
    monkeypatch.setattr(
        fwbs_variables,
        "bz_channel_conduct_liq",
        liquidbreederpressuredropmhdparam.bz_channel_conduct_liq,
    )
    monkeypatch.setattr(
        fwbs_variables,
        "th_wall_secondary",
        liquidbreederpressuredropmhdparam.th_wall_secondary,
    )

    liquid_breeder_pressure_drop_mhd_out = (
        blanket_library_fixture.liquid_breeder_pressure_drop_mhd(
            vel=liquidbreederpressuredropmhdparam.vel,
            vsc=liquidbreederpressuredropmhdparam.vsc,
            conduct_liq=liquidbreederpressuredropmhdparam.conduct_liq,
            l_channel=liquidbreederpressuredropmhdparam.l_channel,
            num_pol=liquidbreederpressuredropmhdparam.num_pol,
            label=liquidbreederpressuredropmhdparam.label,
        )
    )

    assert liquid_breeder_pressure_drop_mhd_out == pytest.approx(
        liquidbreederpressuredropmhdparam.expected_liquid_breeder_pressure_drop_mhd_out
    )
