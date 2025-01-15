from typing import Any, NamedTuple

import numpy
import pytest

from process.blanket_library import BlanketLibrary
from process.fortran import (
    blanket_library,
    build_variables,
    buildings_variables,
    divertor_variables,
    fwbs_variables,
    pfcoil_variables,
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
    fwcoolant: Any = None

    temp_fw_in: Any = None

    temp_fw_out: Any = None

    pres_fw: Any = None

    rhof_fw: Any = None

    cp_fw: Any = None

    cv_fw: Any = None

    coolwh: Any = None

    temp_blkt_in: Any = None

    outlet_temp: Any = None

    blpressure: Any = None

    rhof_bl: Any = None

    i_blkt_dual_coolant: Any = None

    visc_bl: Any = None

    cp_bl: Any = None

    cv_bl: Any = None

    visc_fw: Any = None

    ipump: Any = None

    expected_rhof_fw: Any = None

    expected_cp_fw: Any = None

    expected_cv_fw: Any = None

    expected_rhof_bl: Any = None

    expected_visc_bl: Any = None

    expected_cp_bl: Any = None

    expected_cv_bl: Any = None

    expected_visc_fw: Any = None


@pytest.mark.parametrize(
    "primarycoolantpropertiesparam",
    (
        PrimaryCoolantPropertiesParam(
            fwcoolant="helium",
            temp_fw_in=573,
            temp_fw_out=773,
            pres_fw=8000000,
            rhof_fw=0,
            cp_fw=0,
            cv_fw=0,
            coolwh=1,
            temp_blkt_in=573,
            outlet_temp=773,
            blpressure=8000000,
            rhof_bl=0,
            i_blkt_dual_coolant=2,
            visc_bl=0,
            cp_bl=0,
            cv_bl=0,
            visc_fw=0,
            ipump=0,
            expected_rhof_fw=5.6389735407435868,
            expected_cp_fw=5188.5588430173211,
            expected_cv_fw=3123.5687263525392,
            expected_rhof_bl=5.6389735407435868,
            expected_visc_bl=3.5036293160410249e-05,
            expected_cp_bl=5188.5588430173211,
            expected_cv_bl=3123.5687263525392,
            expected_visc_fw=3.5036293160410249e-05,
        ),
        PrimaryCoolantPropertiesParam(
            fwcoolant="helium",
            temp_fw_in=573,
            temp_fw_out=773,
            pres_fw=8000000,
            rhof_fw=5.6389735407435868,
            cp_fw=5188.5588430173211,
            cv_fw=3123.5687263525392,
            coolwh=1,
            temp_blkt_in=573,
            outlet_temp=773,
            blpressure=8000000,
            rhof_bl=5.6389735407435868,
            i_blkt_dual_coolant=2,
            visc_bl=3.5036293160410249e-05,
            cp_bl=5188.5588430173211,
            cv_bl=3123.5687263525392,
            visc_fw=3.5036293160410249e-05,
            ipump=0,
            expected_rhof_fw=5.6389735407435868,
            expected_cp_fw=5188.5588430173211,
            expected_cv_fw=3123.5687263525392,
            expected_rhof_bl=5.6389735407435868,
            expected_visc_bl=3.5036293160410249e-05,
            expected_cp_bl=5188.5588430173211,
            expected_cv_bl=3123.5687263525392,
            expected_visc_fw=3.5036293160410249e-05,
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
    #     fwbs_variables, "fwcoolant", primarycoolantpropertiesparam.fwcoolant
    # )

    monkeypatch.setattr(
        fwbs_variables, "temp_fw_in", primarycoolantpropertiesparam.temp_fw_in
    )

    monkeypatch.setattr(
        fwbs_variables, "temp_fw_out", primarycoolantpropertiesparam.temp_fw_out
    )

    monkeypatch.setattr(
        fwbs_variables, "pres_fw", primarycoolantpropertiesparam.pres_fw
    )

    monkeypatch.setattr(
        fwbs_variables, "rhof_fw", primarycoolantpropertiesparam.rhof_fw
    )

    monkeypatch.setattr(fwbs_variables, "cp_fw", primarycoolantpropertiesparam.cp_fw)

    monkeypatch.setattr(fwbs_variables, "cv_fw", primarycoolantpropertiesparam.cv_fw)

    monkeypatch.setattr(fwbs_variables, "coolwh", primarycoolantpropertiesparam.coolwh)

    monkeypatch.setattr(
        fwbs_variables, "temp_blkt_in", primarycoolantpropertiesparam.temp_blkt_in
    )

    monkeypatch.setattr(
        fwbs_variables, "outlet_temp", primarycoolantpropertiesparam.outlet_temp
    )

    monkeypatch.setattr(
        fwbs_variables, "blpressure", primarycoolantpropertiesparam.blpressure
    )

    monkeypatch.setattr(
        fwbs_variables, "rhof_bl", primarycoolantpropertiesparam.rhof_bl
    )

    monkeypatch.setattr(
        fwbs_variables,
        "i_blkt_dual_coolant",
        primarycoolantpropertiesparam.i_blkt_dual_coolant,
    )

    monkeypatch.setattr(
        fwbs_variables, "visc_bl", primarycoolantpropertiesparam.visc_bl
    )

    monkeypatch.setattr(fwbs_variables, "cp_bl", primarycoolantpropertiesparam.cp_bl)

    monkeypatch.setattr(fwbs_variables, "cv_bl", primarycoolantpropertiesparam.cv_bl)

    monkeypatch.setattr(
        fwbs_variables, "visc_fw", primarycoolantpropertiesparam.visc_fw
    )

    monkeypatch.setattr(fwbs_variables, "ipump", primarycoolantpropertiesparam.ipump)

    blanket_library_fixture.primary_coolant_properties(output=False)

    assert fwbs_variables.rhof_fw == pytest.approx(
        primarycoolantpropertiesparam.expected_rhof_fw, rel=1e-4
    )

    assert fwbs_variables.cp_fw == pytest.approx(
        primarycoolantpropertiesparam.expected_cp_fw, rel=1e-4
    )

    assert fwbs_variables.cv_fw == pytest.approx(
        primarycoolantpropertiesparam.expected_cv_fw, rel=1e-4
    )

    assert fwbs_variables.rhof_bl == pytest.approx(
        primarycoolantpropertiesparam.expected_rhof_bl, rel=1e-4
    )

    assert fwbs_variables.visc_bl == pytest.approx(
        primarycoolantpropertiesparam.expected_visc_bl, rel=1e-4
    )

    assert fwbs_variables.cp_bl == pytest.approx(
        primarycoolantpropertiesparam.expected_cp_bl, rel=1e-4
    )

    assert fwbs_variables.cv_bl == pytest.approx(
        primarycoolantpropertiesparam.expected_cv_bl, rel=1e-4
    )

    assert fwbs_variables.visc_fw == pytest.approx(
        primarycoolantpropertiesparam.expected_visc_fw, rel=1e-4
    )


def test_deltap_tot_inboard_first_wall(monkeypatch, blanket_library_fixture):
    monkeypatch.setattr(fwbs_variables, "afw", 0.006)
    monkeypatch.setattr(fwbs_variables, "a_bz_liq", 0.22481)

    data = dict(
        icoolpump=1,
        flow_velocity=15.9,
        flleng=4,
        no90=2,
        no180=0,
        coolant_density=5.6,
        coolant_dynamic_viscosity=3.5e-5,
        coolant_electrical_conductivity=0.0,
        pol_channel_length=1.89,
        nopolchan=0,
        label="Inboard first wall",
    )

    assert (
        pytest.approx(blanket_library_fixture.deltap_tot(False, **data))
        == 5885.192672142268
    )


def test_deltap_tot_outboard_blanket_breeder_liquid(
    monkeypatch, blanket_library_fixture
):
    monkeypatch.setattr(fwbs_variables, "afw", 0.006)
    monkeypatch.setattr(fwbs_variables, "a_bz_liq", 0.22481)
    monkeypatch.setattr(fwbs_variables, "ifci", 1)
    monkeypatch.setattr(fwbs_variables, "b_bz_liq", 0.11625)
    monkeypatch.setattr(fwbs_variables, "b_mag_blkt", [8.393, 3.868])
    monkeypatch.setattr(fwbs_variables, "bz_channel_conduct_liq", 833000)
    monkeypatch.setattr(fwbs_variables, "th_wall_secondary", 0.0125)

    data = dict(
        icoolpump=2,
        flow_velocity=0.06,
        flleng=4.7,
        no90=2,
        no180=1,
        coolant_density=9753.2,
        coolant_dynamic_viscosity=0.0017,
        coolant_electrical_conductivity=861800.8,
        pol_channel_length=1.89,
        nopolchan=0,
        label="Outboard blanket breeder liquid",
    )

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
    vgap_xpoint_divertor: Any = None
    vgap_vv_thermalshield: Any = None
    blnktth: Any = None
    shldtth: Any = None
    scrapli: Any = None
    scraplo: Any = None
    fwith: Any = None
    fwoth: Any = None
    d_vv_bot: Any = None
    d_vv_top: Any = None
    rminor: Any = None
    kappa: Any = None
    idivrt: Any = None
    divfix: Any = None
    icomponent: Any = None
    expected_icomponent: Any = None
    expected_half_height: Any = None


@pytest.mark.parametrize(
    "componenthalfheightparam",
    (
        ComponentHalfHeightParam(
            hmax=8.8182171641274945,
            vgap_xpoint_divertor=2.0018838307941582,
            vgap_vv_thermalshield=0.16300000000000001,
            blnktth=0.85000000000000009,
            shldtth=0.59999999999999998,
            scrapli=0.25,
            scraplo=0.25,
            fwith=0.018000000000000002,
            fwoth=0.018000000000000002,
            d_vv_bot=0.30000000000000004,
            d_vv_top=0.30000000000000004,
            rminor=2.6666666666666665,
            kappa=1.8500000000000001,
            idivrt=1,
            divfix=0.62000000000000011,
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
        "vgap_xpoint_divertor",
        componenthalfheightparam.vgap_xpoint_divertor,
    )
    monkeypatch.setattr(
        build_variables,
        "vgap_vv_thermalshield",
        componenthalfheightparam.vgap_vv_thermalshield,
    )
    monkeypatch.setattr(build_variables, "blnktth", componenthalfheightparam.blnktth)
    monkeypatch.setattr(build_variables, "shldtth", componenthalfheightparam.shldtth)
    monkeypatch.setattr(build_variables, "scrapli", componenthalfheightparam.scrapli)
    monkeypatch.setattr(build_variables, "scraplo", componenthalfheightparam.scraplo)
    monkeypatch.setattr(build_variables, "fwith", componenthalfheightparam.fwith)
    monkeypatch.setattr(build_variables, "fwoth", componenthalfheightparam.fwoth)
    monkeypatch.setattr(build_variables, "d_vv_bot", componenthalfheightparam.d_vv_bot)
    monkeypatch.setattr(build_variables, "d_vv_top", componenthalfheightparam.d_vv_top)
    monkeypatch.setattr(physics_variables, "rminor", componenthalfheightparam.rminor)
    monkeypatch.setattr(physics_variables, "kappa", componenthalfheightparam.kappa)
    monkeypatch.setattr(physics_variables, "idivrt", componenthalfheightparam.idivrt)
    monkeypatch.setattr(divertor_variables, "divfix", componenthalfheightparam.divfix)

    half_height = blanket_library_fixture.component_half_height(
        componenthalfheightparam.icomponent
    )

    assert half_height == pytest.approx(componenthalfheightparam.expected_half_height)


class DshapedComponentParam(NamedTuple):
    rsldi: Any = None
    shldith: Any = None
    blnkith: Any = None
    fwith: Any = None
    scrapli: Any = None
    scraplo: Any = None
    fwoth: Any = None
    blareaib: Any = None
    blareaob: Any = None
    blarea: Any = None
    blnkoth: Any = None
    blnktth: Any = None
    shareaib: Any = None
    shareaob: Any = None
    sharea: Any = None
    shldoth: Any = None
    shldtth: Any = None
    rsldo: Any = None
    d_vv_in: Any = None
    d_vv_out: Any = None
    d_vv_top: Any = None
    d_vv_bot: Any = None
    volblkti: Any = None
    volblkto: Any = None
    volblkt: Any = None
    volshld: Any = None
    vdewin: Any = None
    rminor: Any = None
    volshldi: Any = None
    volshldo: Any = None
    volvvi: Any = None
    volvvo: Any = None
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
    expected_volblkto: Any = None
    expected_volblkt: Any = None
    expected_volshld: Any = None
    expected_vdewin: Any = None
    expected_volshldi: Any = None
    expected_volshldo: Any = None
    expected_volvvi: Any = None
    expected_volvvo: Any = None
    expected_icomponent: Any = None


@pytest.mark.parametrize(
    "dshapedcomponentparam",
    (
        DshapedComponentParam(
            rsldi=1.5,
            shldith=0.40000000000000002,
            blnkith=0,
            fwith=0.018000000000000002,
            scrapli=0.10000000000000001,
            scraplo=0.10000000000000001,
            fwoth=0.018000000000000002,
            blareaib=0,
            blareaob=0,
            blarea=0,
            blnkoth=1,
            blnktth=0.5,
            shareaib=0,
            shareaob=0,
            sharea=0,
            shldoth=0.30000000000000004,
            shldtth=0.60000000000000009,
            rsldo=8.4000000000000004,
            d_vv_in=0.20000000000000001,
            d_vv_out=0.30000000000000004,
            d_vv_top=0.30000000000000004,
            d_vv_bot=0.30000000000000004,
            volblkti=0,
            volblkto=0,
            volblkt=0,
            volshld=0,
            vdewin=0,
            rminor=2.5,
            volshldi=0,
            volshldo=0,
            volvvi=0,
            volvvo=0,
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
            expected_volblkto=691.06561956756764,
            expected_volblkt=691.06561956756764,
            expected_volshld=0,
            expected_vdewin=0,
            expected_volshldi=0,
            expected_volshldo=0,
            expected_volvvi=0,
            expected_volvvo=0,
            expected_icomponent=0,
        ),
        DshapedComponentParam(
            rsldi=1.5,
            shldith=0.40000000000000002,
            blnkith=0,
            fwith=0.018000000000000002,
            scrapli=0.10000000000000001,
            scraplo=0.10000000000000001,
            fwoth=0.018000000000000002,
            blareaib=196.97785938008002,
            blareaob=852.24160940262459,
            blarea=1049.2194687827046,
            blnkoth=1,
            blnktth=0.5,
            shareaib=0,
            shareaob=0,
            sharea=0,
            shldoth=0.30000000000000004,
            shldtth=0.60000000000000009,
            rsldo=8.4000000000000004,
            d_vv_in=0.20000000000000001,
            d_vv_out=0.30000000000000004,
            d_vv_top=0.30000000000000004,
            d_vv_bot=0.30000000000000004,
            volblkti=0,
            volblkto=691.06561956756764,
            volblkt=691.06561956756764,
            volshld=0,
            vdewin=0,
            rminor=2.5,
            volshldi=0,
            volshldo=0,
            volvvi=0,
            volvvo=0,
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
            expected_volblkto=691.06561956756764,
            expected_volblkt=691.06561956756764,
            expected_volshld=450.46122947809488,
            expected_vdewin=0,
            expected_volshldi=79.896984366095609,
            expected_volshldo=370.5642451119993,
            expected_volvvi=0,
            expected_volvvo=0,
            expected_icomponent=1,
        ),
        DshapedComponentParam(
            rsldi=1.5,
            shldith=0.40000000000000002,
            blnkith=0,
            fwith=0.018000000000000002,
            scrapli=0.10000000000000001,
            scraplo=0.10000000000000001,
            fwoth=0.018000000000000002,
            blareaib=196.97785938008002,
            blareaob=852.24160940262459,
            blarea=1049.2194687827046,
            blnkoth=1,
            blnktth=0.5,
            shareaib=208.91591146372122,
            shareaob=1013.8483589087293,
            sharea=1222.7642703724505,
            shldoth=0.30000000000000004,
            shldtth=0.60000000000000009,
            rsldo=8.4000000000000004,
            d_vv_in=0.20000000000000001,
            d_vv_out=0.30000000000000004,
            d_vv_top=0.30000000000000004,
            d_vv_bot=0.30000000000000004,
            volblkti=0,
            volblkto=691.06561956756764,
            volblkt=691.06561956756764,
            volshld=450.46122947809488,
            vdewin=0,
            rminor=2.5,
            volshldi=79.896984366095609,
            volshldo=370.5642451119993,
            volvvi=0,
            volvvo=0,
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
            expected_volblkto=691.06561956756764,
            expected_volblkt=691.06561956756764,
            expected_volshld=450.46122947809488,
            expected_vdewin=340.45369594344834,
            expected_volshldi=79.896984366095609,
            expected_volshldo=370.5642451119993,
            expected_volvvi=34.253413020620215,
            expected_volvvo=306.20028292282814,
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
    monkeypatch.setattr(build_variables, "shldith", dshapedcomponentparam.shldith)
    monkeypatch.setattr(build_variables, "blnkith", dshapedcomponentparam.blnkith)
    monkeypatch.setattr(build_variables, "fwith", dshapedcomponentparam.fwith)
    monkeypatch.setattr(build_variables, "scrapli", dshapedcomponentparam.scrapli)
    monkeypatch.setattr(build_variables, "scraplo", dshapedcomponentparam.scraplo)
    monkeypatch.setattr(build_variables, "fwoth", dshapedcomponentparam.fwoth)
    monkeypatch.setattr(build_variables, "blareaib", dshapedcomponentparam.blareaib)
    monkeypatch.setattr(build_variables, "blareaob", dshapedcomponentparam.blareaob)
    monkeypatch.setattr(build_variables, "blarea", dshapedcomponentparam.blarea)
    monkeypatch.setattr(build_variables, "blnkoth", dshapedcomponentparam.blnkoth)
    monkeypatch.setattr(build_variables, "blnktth", dshapedcomponentparam.blnktth)
    monkeypatch.setattr(build_variables, "shareaib", dshapedcomponentparam.shareaib)
    monkeypatch.setattr(build_variables, "shareaob", dshapedcomponentparam.shareaob)
    monkeypatch.setattr(build_variables, "sharea", dshapedcomponentparam.sharea)
    monkeypatch.setattr(build_variables, "shldoth", dshapedcomponentparam.shldoth)
    monkeypatch.setattr(build_variables, "shldtth", dshapedcomponentparam.shldtth)
    monkeypatch.setattr(build_variables, "rsldo", dshapedcomponentparam.rsldo)
    monkeypatch.setattr(build_variables, "d_vv_in", dshapedcomponentparam.d_vv_in)
    monkeypatch.setattr(build_variables, "d_vv_out", dshapedcomponentparam.d_vv_out)
    monkeypatch.setattr(build_variables, "d_vv_top", dshapedcomponentparam.d_vv_top)
    monkeypatch.setattr(build_variables, "d_vv_bot", dshapedcomponentparam.d_vv_bot)
    monkeypatch.setattr(fwbs_variables, "volblkti", dshapedcomponentparam.volblkti)
    monkeypatch.setattr(fwbs_variables, "volblkto", dshapedcomponentparam.volblkto)
    monkeypatch.setattr(fwbs_variables, "volblkt", dshapedcomponentparam.volblkt)
    monkeypatch.setattr(fwbs_variables, "volshld", dshapedcomponentparam.volshld)
    monkeypatch.setattr(fwbs_variables, "vdewin", dshapedcomponentparam.vdewin)
    monkeypatch.setattr(physics_variables, "rminor", dshapedcomponentparam.rminor)
    monkeypatch.setattr(blanket_library, "volshldi", dshapedcomponentparam.volshldi)
    monkeypatch.setattr(blanket_library, "volshldo", dshapedcomponentparam.volshldo)
    monkeypatch.setattr(blanket_library, "volvvi", dshapedcomponentparam.volvvi)
    monkeypatch.setattr(blanket_library, "volvvo", dshapedcomponentparam.volvvo)
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
    assert fwbs_variables.volblkto == pytest.approx(
        dshapedcomponentparam.expected_volblkto
    )
    assert fwbs_variables.volblkt == pytest.approx(
        dshapedcomponentparam.expected_volblkt
    )


class EllipticalComponentParam(NamedTuple):
    rsldi: Any = None
    shldith: Any = None
    blnkith: Any = None
    rsldo: Any = None
    shldoth: Any = None
    blnkoth: Any = None
    blareaib: Any = None
    blareaob: Any = None
    blarea: Any = None
    blnktth: Any = None
    shareaib: Any = None
    shareaob: Any = None
    sharea: Any = None
    shldtth: Any = None
    d_vv_in: Any = None
    d_vv_out: Any = None
    d_vv_top: Any = None
    d_vv_bot: Any = None
    volblkti: Any = None
    volblkto: Any = None
    volblkt: Any = None
    volshld: Any = None
    vdewin: Any = None
    rmajor: Any = None
    rminor: Any = None
    triang: Any = None
    volshldi: Any = None
    volshldo: Any = None
    volvvi: Any = None
    volvvo: Any = None
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
    expected_volblkti: Any = None
    expected_volblkto: Any = None
    expected_volblkt: Any = None
    expected_volshld: Any = None
    expected_vdewin: Any = None
    expected_volshldi: Any = None
    expected_volshldo: Any = None
    expected_volvvi: Any = None
    expected_volvvo: Any = None
    expected_icomponent: Any = None


@pytest.mark.parametrize(
    "ellipticalcomponentparam",
    (
        EllipticalComponentParam(
            rsldi=4.0833333333333339,
            shldith=0.30000000000000004,
            blnkith=0.70000000000000007,
            rsldo=12.716666666666667,
            shldoth=0.80000000000000004,
            blnkoth=1,
            blareaib=0,
            blareaob=0,
            blarea=0,
            blnktth=0.85000000000000009,
            shareaib=0,
            shareaob=0,
            sharea=0,
            shldtth=0.59999999999999998,
            d_vv_in=0.30000000000000004,
            d_vv_out=0.30000000000000004,
            d_vv_top=0.30000000000000004,
            d_vv_bot=0.30000000000000004,
            volblkti=0,
            volblkto=0,
            volblkt=0,
            volshld=0,
            vdewin=0,
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            volshldi=0,
            volshldo=0,
            volvvi=0,
            volvvo=0,
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
            expected_volblkti=315.83946385183026,
            expected_volblkto=1020.3677420460117,
            expected_volblkt=1336.207205897842,
            expected_volshld=0,
            expected_vdewin=0,
            expected_volshldi=0,
            expected_volshldo=0,
            expected_volvvi=0,
            expected_volvvo=0,
            expected_icomponent=0,
        ),
        EllipticalComponentParam(
            rsldi=4.0833333333333339,
            shldith=0.30000000000000004,
            blnkith=0.70000000000000007,
            rsldo=12.716666666666667,
            shldoth=0.80000000000000004,
            blnkoth=1,
            blareaib=664.9687712975541,
            blareaob=1101.3666396424403,
            blarea=1766.3354109399943,
            blnktth=0.85000000000000009,
            shareaib=0,
            shareaob=0,
            sharea=0,
            shldtth=0.59999999999999998,
            d_vv_in=0.30000000000000004,
            d_vv_out=0.30000000000000004,
            d_vv_top=0.30000000000000004,
            d_vv_bot=0.30000000000000004,
            volblkti=315.83946385183026,
            volblkto=1020.3677420460117,
            volblkt=1336.207205897842,
            volshld=0,
            vdewin=0,
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            volshldi=0,
            volshldo=0,
            volvvi=0,
            volvvo=0,
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
            expected_volblkti=315.83946385183026,
            expected_volblkto=1020.3677420460117,
            expected_volblkt=1336.207205897842,
            expected_volshld=1124.4621612595051,
            expected_vdewin=0,
            expected_volshldi=177.89822933168091,
            expected_volshldo=946.56393192782434,
            expected_volvvi=0,
            expected_volvvo=0,
            expected_icomponent=1,
        ),
        EllipticalComponentParam(
            rsldi=4.0833333333333339,
            shldith=0.30000000000000004,
            blnkith=0.70000000000000007,
            rsldo=12.716666666666667,
            shldoth=0.80000000000000004,
            blnkoth=1,
            blareaib=664.9687712975541,
            blareaob=1101.3666396424403,
            blarea=1766.3354109399943,
            blnktth=0.85000000000000009,
            shareaib=700.06731267447844,
            shareaob=1344.1106481995357,
            sharea=2044.1779608740142,
            shldtth=0.59999999999999998,
            d_vv_in=0.30000000000000004,
            d_vv_out=0.30000000000000004,
            d_vv_top=0.30000000000000004,
            d_vv_bot=0.30000000000000004,
            volblkti=315.83946385183026,
            volblkto=1020.3677420460117,
            volblkt=1336.207205897842,
            volshld=1124.4621612595051,
            vdewin=0,
            rmajor=8,
            rminor=2.6666666666666665,
            triang=0.5,
            volshldi=177.89822933168091,
            volshldo=946.56393192782434,
            volvvi=0,
            volvvo=0,
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
            expected_volblkti=315.83946385183026,
            expected_volblkto=1020.3677420460117,
            expected_volblkt=1336.207205897842,
            expected_volshld=1124.4621612595051,
            expected_vdewin=584.07334775041659,
            expected_volshldi=177.89822933168091,
            expected_volshldo=946.56393192782434,
            expected_volvvi=143.03162449152501,
            expected_volvvo=441.04172325889158,
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
    monkeypatch.setattr(build_variables, "shldith", ellipticalcomponentparam.shldith)
    monkeypatch.setattr(build_variables, "blnkith", ellipticalcomponentparam.blnkith)
    monkeypatch.setattr(build_variables, "rsldo", ellipticalcomponentparam.rsldo)
    monkeypatch.setattr(build_variables, "shldoth", ellipticalcomponentparam.shldoth)
    monkeypatch.setattr(build_variables, "blnkoth", ellipticalcomponentparam.blnkoth)
    monkeypatch.setattr(build_variables, "blareaib", ellipticalcomponentparam.blareaib)
    monkeypatch.setattr(build_variables, "blareaob", ellipticalcomponentparam.blareaob)
    monkeypatch.setattr(build_variables, "blarea", ellipticalcomponentparam.blarea)
    monkeypatch.setattr(build_variables, "blnktth", ellipticalcomponentparam.blnktth)
    monkeypatch.setattr(build_variables, "shareaib", ellipticalcomponentparam.shareaib)
    monkeypatch.setattr(build_variables, "shareaob", ellipticalcomponentparam.shareaob)
    monkeypatch.setattr(build_variables, "sharea", ellipticalcomponentparam.sharea)
    monkeypatch.setattr(build_variables, "shldtth", ellipticalcomponentparam.shldtth)
    monkeypatch.setattr(build_variables, "d_vv_in", ellipticalcomponentparam.d_vv_in)
    monkeypatch.setattr(build_variables, "d_vv_out", ellipticalcomponentparam.d_vv_out)
    monkeypatch.setattr(build_variables, "d_vv_top", ellipticalcomponentparam.d_vv_top)
    monkeypatch.setattr(build_variables, "d_vv_bot", ellipticalcomponentparam.d_vv_bot)
    monkeypatch.setattr(fwbs_variables, "volblkti", ellipticalcomponentparam.volblkti)
    monkeypatch.setattr(fwbs_variables, "volblkto", ellipticalcomponentparam.volblkto)
    monkeypatch.setattr(fwbs_variables, "volblkt", ellipticalcomponentparam.volblkt)
    monkeypatch.setattr(fwbs_variables, "volshld", ellipticalcomponentparam.volshld)
    monkeypatch.setattr(fwbs_variables, "vdewin", ellipticalcomponentparam.vdewin)
    monkeypatch.setattr(physics_variables, "rmajor", ellipticalcomponentparam.rmajor)
    monkeypatch.setattr(physics_variables, "rminor", ellipticalcomponentparam.rminor)
    monkeypatch.setattr(physics_variables, "triang", ellipticalcomponentparam.triang)
    monkeypatch.setattr(blanket_library, "volshldi", ellipticalcomponentparam.volshldi)
    monkeypatch.setattr(blanket_library, "volshldo", ellipticalcomponentparam.volshldo)
    monkeypatch.setattr(blanket_library, "volvvi", ellipticalcomponentparam.volvvi)
    monkeypatch.setattr(blanket_library, "volvvo", ellipticalcomponentparam.volvvo)
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
    assert fwbs_variables.volblkti == pytest.approx(
        ellipticalcomponentparam.expected_volblkti
    )
    assert fwbs_variables.volblkto == pytest.approx(
        ellipticalcomponentparam.expected_volblkto
    )
    assert fwbs_variables.volblkt == pytest.approx(
        ellipticalcomponentparam.expected_volblkt
    )
    assert fwbs_variables.volshld == pytest.approx(
        ellipticalcomponentparam.expected_volshld
    )
    assert fwbs_variables.vdewin == pytest.approx(
        ellipticalcomponentparam.expected_vdewin
    )
    assert blanket_library.volshldi == pytest.approx(
        ellipticalcomponentparam.expected_volshldi
    )
    assert blanket_library.volshldo == pytest.approx(
        ellipticalcomponentparam.expected_volshldo
    )
    assert blanket_library.volvvi == pytest.approx(
        ellipticalcomponentparam.expected_volvvi
    )
    assert blanket_library.volvvo == pytest.approx(
        ellipticalcomponentparam.expected_volvvo
    )


class ApplyCoverageFactorsParam(NamedTuple):
    blareaob: Any = None
    blarea: Any = None
    blareaib: Any = None
    shareaib: Any = None
    shareaob: Any = None
    sharea: Any = None
    fdiv: Any = None
    fhcd: Any = None
    volblkto: Any = None
    volblkti: Any = None
    volblkt: Any = None
    fvolsi: Any = None
    fvolso: Any = None
    volshld: Any = None
    vdewin: Any = None
    fvoldw: Any = None
    idivrt: Any = None
    volshldi: Any = None
    volshldo: Any = None
    expected_blareaob: Any = None
    expected_blarea: Any = None
    expected_shareaob: Any = None
    expected_sharea: Any = None
    expected_volblkto: Any = None
    expected_volblkt: Any = None
    expected_volshld: Any = None
    expected_vdewin: Any = None
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
            fdiv=0.115,
            fhcd=0,
            volblkto=1020.3677420460117,
            volblkti=315.83946385183026,
            volblkt=1336.207205897842,
            fvolsi=1,
            fvolso=0.64000000000000001,
            volshld=1124.4621612595051,
            vdewin=584.07334775041659,
            fvoldw=1.74,
            idivrt=1,
            volshldi=177.89822933168091,
            volshldo=946.56393192782434,
            expected_blareaob=898.23806738434075,
            expected_blarea=1563.2068386818949,
            expected_shareaob=860.23081484770285,
            expected_sharea=1560.2981275221814,
            expected_volblkto=866.70391336775992,
            expected_volblkt=1182.5433772195902,
            expected_volshld=783.69914576548854,
            expected_vdewin=1016.2876250857248,
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
    monkeypatch.setattr(fwbs_variables, "fdiv", applycoveragefactorsparam.fdiv)
    monkeypatch.setattr(fwbs_variables, "fhcd", applycoveragefactorsparam.fhcd)
    monkeypatch.setattr(fwbs_variables, "volblkto", applycoveragefactorsparam.volblkto)
    monkeypatch.setattr(fwbs_variables, "volblkti", applycoveragefactorsparam.volblkti)
    monkeypatch.setattr(fwbs_variables, "volblkt", applycoveragefactorsparam.volblkt)
    monkeypatch.setattr(fwbs_variables, "fvolsi", applycoveragefactorsparam.fvolsi)
    monkeypatch.setattr(fwbs_variables, "fvolso", applycoveragefactorsparam.fvolso)
    monkeypatch.setattr(fwbs_variables, "volshld", applycoveragefactorsparam.volshld)
    monkeypatch.setattr(fwbs_variables, "vdewin", applycoveragefactorsparam.vdewin)
    monkeypatch.setattr(fwbs_variables, "fvoldw", applycoveragefactorsparam.fvoldw)
    monkeypatch.setattr(physics_variables, "idivrt", applycoveragefactorsparam.idivrt)
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
    assert fwbs_variables.volblkto == pytest.approx(
        applycoveragefactorsparam.expected_volblkto
    )
    assert fwbs_variables.volblkt == pytest.approx(
        applycoveragefactorsparam.expected_volblkt
    )
    assert fwbs_variables.volshld == pytest.approx(
        applycoveragefactorsparam.expected_volshld
    )
    assert fwbs_variables.vdewin == pytest.approx(
        applycoveragefactorsparam.expected_vdewin
    )
    assert blanket_library.volshldo == pytest.approx(
        applycoveragefactorsparam.expected_volshldo
    )


class ExternalCryoGeometryParam(NamedTuple):
    clhsf: Any = None
    hmax: Any = None
    tfcth: Any = None
    ddwex: Any = None
    rdewex: Any = None
    rpf2dewar: Any = None
    zdewex: Any = None
    vdewex: Any = None
    vvmass: Any = None
    vdewin: Any = None
    denstl: Any = None
    dewmkg: Any = None
    rb: Any = None
    zh: Any = None
    clh1: Any = None
    hcryopf: Any = None
    expected_rdewex: Any = None
    expected_zdewex: Any = None
    expected_vdewex: Any = None
    expected_vvmass: Any = None
    expected_dewmkg: Any = None
    expected_clh1: Any = None
    expected_hcryopf: Any = None


@pytest.mark.parametrize(
    "externalcryogeometryparam",
    (
        ExternalCryoGeometryParam(
            clhsf=4.2679999999999998,
            hmax=8.8182171641274945,
            tfcth=0.92672586247397692,
            ddwex=0.15000000000000002,
            rdewex=0,
            rpf2dewar=0.5,
            zdewex=0,
            vdewex=0,
            vvmass=0,
            vdewin=1016.2876250857248,
            denstl=7800,
            dewmkg=0,
            rb=numpy.array(
                numpy.array(
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
            zh=numpy.array(
                numpy.array(
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
            clh1=2.5,
            hcryopf=0,
            expected_rdewex=17.805470903073743,
            expected_zdewex=15.259637557000296,
            expected_vdewex=810.95017445543203,
            expected_vvmass=7927043.4756686538,
            expected_dewmkg=14252454.836421024,
            expected_clh1=5.514694530398824,
            expected_hcryopf=5.3441455565624985,
        ),
    ),
)
def test_external_cryo_geometry(
    externalcryogeometryparam, monkeypatch, blanket_library_fixture
):
    """
    Automatically generated Regression Unit Test for external_cryo_geometry.

    This test was generated using data from tests/regression/input_files/large_tokamak_once_through.IN.DAT.

    :param externalcryogeometryparam: the data used to mock and assert in this test.
    :type externalcryogeometryparam: externalcryogeometryparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(build_variables, "clhsf", externalcryogeometryparam.clhsf)
    monkeypatch.setattr(build_variables, "hmax", externalcryogeometryparam.hmax)
    monkeypatch.setattr(build_variables, "tfcth", externalcryogeometryparam.tfcth)
    monkeypatch.setattr(build_variables, "ddwex", externalcryogeometryparam.ddwex)
    monkeypatch.setattr(fwbs_variables, "rdewex", externalcryogeometryparam.rdewex)
    monkeypatch.setattr(
        fwbs_variables, "rpf2dewar", externalcryogeometryparam.rpf2dewar
    )
    monkeypatch.setattr(fwbs_variables, "zdewex", externalcryogeometryparam.zdewex)
    monkeypatch.setattr(fwbs_variables, "vdewex", externalcryogeometryparam.vdewex)
    monkeypatch.setattr(fwbs_variables, "vvmass", externalcryogeometryparam.vvmass)
    monkeypatch.setattr(fwbs_variables, "vdewin", externalcryogeometryparam.vdewin)
    monkeypatch.setattr(fwbs_variables, "denstl", externalcryogeometryparam.denstl)
    monkeypatch.setattr(fwbs_variables, "dewmkg", externalcryogeometryparam.dewmkg)
    monkeypatch.setattr(pfcoil_variables, "rb", externalcryogeometryparam.rb)
    monkeypatch.setattr(pfcoil_variables, "zh", externalcryogeometryparam.zh)
    monkeypatch.setattr(buildings_variables, "clh1", externalcryogeometryparam.clh1)
    monkeypatch.setattr(blanket_library, "hcryopf", externalcryogeometryparam.hcryopf)

    blanket_library_fixture.external_cryo_geometry()

    assert fwbs_variables.rdewex == pytest.approx(
        externalcryogeometryparam.expected_rdewex
    )
    assert fwbs_variables.zdewex == pytest.approx(
        externalcryogeometryparam.expected_zdewex
    )
    assert fwbs_variables.vdewex == pytest.approx(
        externalcryogeometryparam.expected_vdewex
    )
    assert fwbs_variables.vvmass == pytest.approx(
        externalcryogeometryparam.expected_vvmass
    )
    assert fwbs_variables.dewmkg == pytest.approx(
        externalcryogeometryparam.expected_dewmkg
    )
    assert buildings_variables.clh1 == pytest.approx(
        externalcryogeometryparam.expected_clh1
    )
    assert blanket_library.hcryopf == pytest.approx(
        externalcryogeometryparam.expected_hcryopf
    )


class BlanketModPolHeightParam(NamedTuple):
    scrapli: Any = None
    scraplo: Any = None
    fwbsshape: Any = None
    nblktmodpi: Any = None
    fdiv: Any = None
    nblktmodpo: Any = None
    itart: Any = None
    rminor: Any = None
    idivrt: Any = None
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
            scrapli=0.25,
            scraplo=0.25,
            fwbsshape=2,
            nblktmodpi=7,
            fdiv=0.115,
            nblktmodpo=8,
            itart=0,
            rminor=2.6666666666666665,
            idivrt=1,
            rmajor=8,
            triang=0.5,
            bllengi=0,
            bllengo=0,
            hblnkt=5.9532752487304119,
            expected_bllengi=1.6252823720672551,
            expected_bllengo=1.7853902013340495,
        ),
        BlanketModPolHeightParam(
            scrapli=0.10000000000000001,
            scraplo=0.10000000000000001,
            fwbsshape=1,
            nblktmodpi=7,
            fdiv=0.115,
            nblktmodpo=8,
            itart=1,
            rminor=2.5,
            idivrt=2,
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
    monkeypatch.setattr(build_variables, "scrapli", blanketmodpolheightparam.scrapli)
    monkeypatch.setattr(build_variables, "scraplo", blanketmodpolheightparam.scraplo)
    monkeypatch.setattr(fwbs_variables, "fwbsshape", blanketmodpolheightparam.fwbsshape)
    monkeypatch.setattr(
        fwbs_variables, "nblktmodpi", blanketmodpolheightparam.nblktmodpi
    )
    monkeypatch.setattr(fwbs_variables, "fdiv", blanketmodpolheightparam.fdiv)
    monkeypatch.setattr(
        fwbs_variables, "nblktmodpo", blanketmodpolheightparam.nblktmodpo
    )
    monkeypatch.setattr(physics_variables, "itart", blanketmodpolheightparam.itart)
    monkeypatch.setattr(physics_variables, "rminor", blanketmodpolheightparam.rminor)
    monkeypatch.setattr(physics_variables, "idivrt", blanketmodpolheightparam.idivrt)
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
    i_bb_liq: Any = None
    hartmann_liq: Any = None
    b_mag_blkt: Any = None
    iblnkith: Any = None
    i_blkt_dual_coolant: Any = None
    bt: Any = None
    aspect: Any = None
    rmajor: Any = None
    blnkith: Any = None
    blnkoth: Any = None
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
            i_bb_liq=0,
            hartmann_liq=numpy.array(
                numpy.array((0, 0), order="F"), order="F"
            ).transpose(),
            b_mag_blkt=numpy.array(
                numpy.array((5, 5), order="F"), order="F"
            ).transpose(),
            iblnkith=1,
            i_blkt_dual_coolant=0,
            bt=5.7000000000000002,
            aspect=3,
            rmajor=8,
            blnkith=0.70000000000000007,
            blnkoth=1,
            ip=0,
            ofile=11,
            expected_den_liq=9753.2497999999996,
            expected_specific_heat_liq=189.12018,
            expected_thermal_conductivity_liq=9.238260167312621,
            expected_dynamic_viscosity_liq=0.0017477589255667965,
            expected_electrical_conductivity_liq=861800.80431007256,
            expected_hartmann_liq=numpy.array(
                numpy.array((20319.245984309102, 9067.8426109080938), order="F"),
                order="F",
            ).transpose(),
            expected_b_mag_blkt=numpy.array(
                numpy.array((9.1505016722408019, 4.0835820895522392), order="F"),
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
            expected_thermal_conductivity_liq=30,  # doesn't change when i_bb_liq=1
            dynamic_viscosity_liq=0,
            electrical_conductivity_liq=0,
            i_bb_liq=1,
            hartmann_liq=numpy.array(
                numpy.array((0, 0), order="F"), order="F"
            ).transpose(),
            b_mag_blkt=numpy.array(
                numpy.array((5, 5), order="F"), order="F"
            ).transpose(),
            iblnkith=1,
            i_blkt_dual_coolant=0,
            bt=5.7000000000000002,
            aspect=3,
            rmajor=8,
            blnkith=0.70000000000000007,
            blnkoth=1,
            ip=0,
            ofile=11,
            expected_den_liq=305.30702851374997,
            expected_specific_heat_liq=34.640761200690406,
            expected_dynamic_viscosity_liq=0.00037298826343426359,
            expected_electrical_conductivity_liq=596562356750.5,
            expected_hartmann_liq=numpy.array(
                numpy.array((36595294.326740541, 16331332.841336451), order="F"),
                order="F",
            ).transpose(),
            expected_b_mag_blkt=numpy.array(
                numpy.array((9.1505016722408019, 4.0835820895522392), order="F"),
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
        fwbs_variables, "i_bb_liq", liquidbreederpropertiesparam.i_bb_liq
    )
    monkeypatch.setattr(
        fwbs_variables, "hartmann_liq", liquidbreederpropertiesparam.hartmann_liq
    )
    monkeypatch.setattr(
        fwbs_variables, "b_mag_blkt", liquidbreederpropertiesparam.b_mag_blkt
    )
    monkeypatch.setattr(
        fwbs_variables, "iblnkith", liquidbreederpropertiesparam.iblnkith
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
        build_variables, "blnkith", liquidbreederpropertiesparam.blnkith
    )
    monkeypatch.setattr(
        build_variables, "blnkoth", liquidbreederpropertiesparam.blnkoth
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
    afw: Any = None
    a_bz_liq: Any = None
    b_bz_liq: Any = None
    roughness: Any = None
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
            afw=0.0060000000000000001,
            a_bz_liq=0.20000000000000001,
            b_bz_liq=0.20000000000000001,
            roughness=9.9999999999999995e-07,
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
    monkeypatch.setattr(fwbs_variables, "afw", pressuredropparam.afw)
    monkeypatch.setattr(fwbs_variables, "a_bz_liq", pressuredropparam.a_bz_liq)
    monkeypatch.setattr(fwbs_variables, "b_bz_liq", pressuredropparam.b_bz_liq)
    monkeypatch.setattr(fwbs_variables, "roughness", pressuredropparam.roughness)

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
    ifci: Any = None
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
            ifci=0,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=numpy.array(
                numpy.array((8.3930173480023953, 3.8678755427951117), order="F"),
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
            ifci=1,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=numpy.array(
                numpy.array((8.3930173480023953, 3.8678755427951117), order="F"),
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
            ifci=2,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=numpy.array(
                numpy.array((8.3930173480023953, 3.8678755427951117), order="F"),
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
            ifci=2,
            a_bz_liq=0.22481469639955909,
            b_bz_liq=0.11625000000000001,
            b_mag_blkt=numpy.array(
                numpy.array((8.3930173480023953, 3.8678755427951117), order="F"),
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
    monkeypatch.setattr(fwbs_variables, "ifci", liquidbreederpressuredropmhdparam.ifci)
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
