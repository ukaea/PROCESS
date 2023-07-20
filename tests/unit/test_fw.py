import pytest
from typing import NamedTuple, Any

from process.fw import Fw
from process.fortran import fwbs_variables


@pytest.fixture
def fw():
    """Provides Fw object for testing.

    :returns: initialised Fw object
    :rtype: process.fw.Fw
    """
    return Fw()


class FwTempParam(NamedTuple):
    fw_th_conductivity: Any = None

    fwcoolant: Any = None

    fwinlet: Any = None

    fwpressure: Any = None

    fwoutlet: Any = None

    pitch: Any = None

    fw_channel_length: Any = None

    tpeak: Any = None

    peaking_factor: Any = None

    fw_wall: Any = None

    pnuc_deposited: Any = None

    afw: Any = None

    thickness: Any = None

    area: Any = None

    prad_incident: Any = None

    label: str = None

    expected_tpeakfw: Any = None

    expected_cfmean: Any = None

    expected_rhofmean: Any = None

    expected_massrate: Any = None


@pytest.mark.parametrize(
    "fwtempparam",
    (
        FwTempParam(
            fw_th_conductivity=28.34,
            fwcoolant="helium",
            fwinlet=573,
            fwpressure=8000000,
            fwoutlet=773,
            pitch=0.005000000000000001,
            fw_channel_length=4,
            tpeak=873,
            peaking_factor=1,
            fw_wall=0.0030000000000000001,
            pnuc_deposited=75.219932653459054,
            afw=0.0060000000000000001,
            thickness=0.018000000000000002,
            area=612.23369764444396,
            prad_incident=97.271629070225231,
            label="Inboard first wall",
            expected_tpeakfw=827.22264979106728,
            expected_cfmean=5188.6731857247905,
            expected_rhofmean=5.7616377393642955,
            expected_massrate=0.0054299309578952218,
        ),
        FwTempParam(
            fw_th_conductivity=28.34,
            fwcoolant="helium",
            fwinlet=573,
            fwpressure=8000000,
            fwoutlet=773,
            pitch=0.005000000000000001,
            fw_channel_length=4,
            tpeak=873,
            peaking_factor=1,
            fw_wall=0.0030000000000000001,
            pnuc_deposited=121.50088652655793,
            afw=0.0060000000000000001,
            thickness=0.018000000000000002,
            area=988.92586580655245,
            prad_incident=176.95628839065773,
            label="Outboard first wall",
            expected_tpeakfw=829.53877268685892,
            expected_cfmean=5188.6731857247905,
            expected_rhofmean=5.7616377393642955,
            expected_massrate=0.005816503189155049,
        ),
    ),
)
def test_fw_temp(fwtempparam, monkeypatch, fw):
    """
    Automatically generated Regression Unit Test for fw_temp.

    This test was generated using data from dcll/dcll_mms_lt_IN.DAT.

    :param fwtempparam: the data used to mock and assert in this test.
    :type fwtempparam: fwtempparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        fwbs_variables, "fw_th_conductivity", fwtempparam.fw_th_conductivity
    )

    # monkeypatch doesnt work for strings
    # but helium is the default
    # monkeypatch.setattr(fwbs_variables, "fwcoolant", fwtempparam.fwcoolant)

    monkeypatch.setattr(fwbs_variables, "fwinlet", fwtempparam.fwinlet)

    monkeypatch.setattr(fwbs_variables, "fwpressure", fwtempparam.fwpressure)

    monkeypatch.setattr(fwbs_variables, "fwoutlet", fwtempparam.fwoutlet)

    monkeypatch.setattr(fwbs_variables, "pitch", fwtempparam.pitch)

    monkeypatch.setattr(
        fwbs_variables, "fw_channel_length", fwtempparam.fw_channel_length
    )

    monkeypatch.setattr(fwbs_variables, "tpeak", fwtempparam.tpeak)

    monkeypatch.setattr(fwbs_variables, "peaking_factor", fwtempparam.peaking_factor)

    monkeypatch.setattr(fwbs_variables, "fw_wall", fwtempparam.fw_wall)

    tpeakfw, cfmean, rhofmean, massrate = fw.fw_temp(
        False,
        pnuc_deposited=fwtempparam.pnuc_deposited,
        afw=fwtempparam.afw,
        thickness=fwtempparam.thickness,
        area=fwtempparam.area,
        prad_incident=fwtempparam.prad_incident,
        label=fwtempparam.label,
    )

    assert tpeakfw == pytest.approx(fwtempparam.expected_tpeakfw)

    assert cfmean == pytest.approx(fwtempparam.expected_cfmean, rel=1e-4)

    assert rhofmean == pytest.approx(fwtempparam.expected_rhofmean, rel=1e-4)

    assert massrate == pytest.approx(fwtempparam.expected_massrate, rel=1e-4)
