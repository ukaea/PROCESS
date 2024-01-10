import pytest
from typing import NamedTuple, Any

from process.blanket_library import BlanketLibrary
from process.fw import Fw
from process.fortran import fwbs_variables


@pytest.fixture
def blanket_library():
    """Provides BlanketLibrary object for testing.

    :returns: initialised BlanketLibrary object
    :rtype: process.blanket_library.BlanketLibrary
    """
    return BlanketLibrary(Fw())


class PrimaryCoolantPropertiesParam(NamedTuple):
    fwcoolant: Any = None

    fwinlet: Any = None

    fwoutlet: Any = None

    fwpressure: Any = None

    rhof_fw: Any = None

    cp_fw: Any = None

    cv_fw: Any = None

    coolwh: Any = None

    inlet_temp: Any = None

    outlet_temp: Any = None

    blpressure: Any = None

    rhof_bl: Any = None

    icooldual: Any = None

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
            fwinlet=573,
            fwoutlet=773,
            fwpressure=8000000,
            rhof_fw=0,
            cp_fw=0,
            cv_fw=0,
            coolwh=1,
            inlet_temp=573,
            outlet_temp=773,
            blpressure=8000000,
            rhof_bl=0,
            icooldual=2,
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
            fwinlet=573,
            fwoutlet=773,
            fwpressure=8000000,
            rhof_fw=5.6389735407435868,
            cp_fw=5188.5588430173211,
            cv_fw=3123.5687263525392,
            coolwh=1,
            inlet_temp=573,
            outlet_temp=773,
            blpressure=8000000,
            rhof_bl=5.6389735407435868,
            icooldual=2,
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
    primarycoolantpropertiesparam, monkeypatch, blanket_library
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
        fwbs_variables, "fwinlet", primarycoolantpropertiesparam.fwinlet
    )

    monkeypatch.setattr(
        fwbs_variables, "fwoutlet", primarycoolantpropertiesparam.fwoutlet
    )

    monkeypatch.setattr(
        fwbs_variables, "fwpressure", primarycoolantpropertiesparam.fwpressure
    )

    monkeypatch.setattr(
        fwbs_variables, "rhof_fw", primarycoolantpropertiesparam.rhof_fw
    )

    monkeypatch.setattr(fwbs_variables, "cp_fw", primarycoolantpropertiesparam.cp_fw)

    monkeypatch.setattr(fwbs_variables, "cv_fw", primarycoolantpropertiesparam.cv_fw)

    monkeypatch.setattr(fwbs_variables, "coolwh", primarycoolantpropertiesparam.coolwh)

    monkeypatch.setattr(
        fwbs_variables, "inlet_temp", primarycoolantpropertiesparam.inlet_temp
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
        fwbs_variables, "icooldual", primarycoolantpropertiesparam.icooldual
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

    blanket_library.primary_coolant_properties(output=False)

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


def test_deltap_tot_inboard_first_wall(monkeypatch, blanket_library):
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

    assert pytest.approx(blanket_library.deltap_tot(False, **data)) == 5885.192672142268


def test_deltap_tot_outboard_blanket_breeder_liquid(monkeypatch, blanket_library):
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
        pytest.approx(blanket_library.deltap_tot(False, **data)) == 56.962742615936264
    )


def test_pumppower_primary_helium(monkeypatch, blanket_library):
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

    assert pytest.approx(blanket_library.pumppower(False, **data)) == 1.8251284651310427


def test_pumppower_secondary_pb_li(monkeypatch, blanket_library):
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
        pytest.approx(blanket_library.pumppower(False, **data), rel=1e-4)
        == 3.2374845432302464
    )
