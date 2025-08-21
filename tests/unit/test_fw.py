from typing import Any, NamedTuple

import pytest

from process.data_structure import fwbs_variables
from process.fw import Fw


@pytest.fixture
def fw():
    """Provides Fw object for testing.

    :returns: initialised Fw object
    :rtype: process.fw.Fw
    """
    return Fw()


class FwTempParam(NamedTuple):
    fw_th_conductivity: Any = None

    i_fw_coolant_type: Any = None

    temp_fw_coolant_in: Any = None

    pres_fw_coolant: Any = None

    temp_fw_coolant_out: Any = None

    dx_fw_module: Any = None

    len_fw_channel: Any = None

    temp_fw_peak: Any = None

    f_fw_peak: Any = None

    dr_fw_wall: Any = None

    pnuc_deposited: Any = None

    radius_fw_channel: Any = None

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
            i_fw_coolant_type="helium",
            temp_fw_coolant_in=573,
            pres_fw_coolant=8000000,
            temp_fw_coolant_out=773,
            dx_fw_module=0.005000000000000001,
            len_fw_channel=4,
            temp_fw_peak=873,
            f_fw_peak=1,
            dr_fw_wall=0.0030000000000000001,
            pnuc_deposited=75.219932653459054,
            radius_fw_channel=0.0060000000000000001,
            thickness=0.018000000000000002,
            area=612.23369764444396,
            prad_incident=97.271629070225231,
            label="Inboard first wall",
            expected_tpeakfw=846.4965128389427,
            expected_cfmean=5188.6731857247905,
            expected_rhofmean=5.7616377393642955,
            expected_massrate=0.0054299309578952218,
        ),
        FwTempParam(
            fw_th_conductivity=28.34,
            i_fw_coolant_type="helium",
            temp_fw_coolant_in=573,
            pres_fw_coolant=8000000,
            temp_fw_coolant_out=773,
            dx_fw_module=0.005000000000000001,
            len_fw_channel=4,
            temp_fw_peak=873,
            f_fw_peak=1,
            dr_fw_wall=0.0030000000000000001,
            pnuc_deposited=121.50088652655793,
            radius_fw_channel=0.0060000000000000001,
            thickness=0.018000000000000002,
            area=988.92586580655245,
            prad_incident=176.95628839065773,
            label="Outboard first wall",
            expected_tpeakfw=850.851622254886,
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
    # monkeypatch.setattr(fwbs_variables, "i_fw_coolant_type", fwtempparam.i_fw_coolant_type)

    monkeypatch.setattr(
        fwbs_variables, "temp_fw_coolant_in", fwtempparam.temp_fw_coolant_in
    )

    monkeypatch.setattr(fwbs_variables, "pres_fw_coolant", fwtempparam.pres_fw_coolant)

    monkeypatch.setattr(
        fwbs_variables, "temp_fw_coolant_out", fwtempparam.temp_fw_coolant_out
    )

    monkeypatch.setattr(fwbs_variables, "dx_fw_module", fwtempparam.dx_fw_module)

    monkeypatch.setattr(fwbs_variables, "len_fw_channel", fwtempparam.len_fw_channel)

    monkeypatch.setattr(fwbs_variables, "temp_fw_peak", fwtempparam.temp_fw_peak)

    monkeypatch.setattr(fwbs_variables, "f_fw_peak", fwtempparam.f_fw_peak)

    monkeypatch.setattr(fwbs_variables, "dr_fw_wall", fwtempparam.dr_fw_wall)

    tpeakfw, cfmean, rhofmean, massrate = fw.fw_temp(
        False,
        pnuc_deposited=fwtempparam.pnuc_deposited,
        radius_fw_channel=fwtempparam.radius_fw_channel,
        dr_fw=fwtempparam.thickness,
        a_fw=fwtempparam.area,
        prad_incident=fwtempparam.prad_incident,
        label=fwtempparam.label,
    )

    assert tpeakfw == pytest.approx(fwtempparam.expected_tpeakfw)

    assert cfmean == pytest.approx(fwtempparam.expected_cfmean, rel=1e-4)

    assert rhofmean == pytest.approx(fwtempparam.expected_rhofmean, rel=1e-4)

    assert massrate == pytest.approx(fwtempparam.expected_massrate, rel=1e-4)


def test_darcy_friction_haaland(fw):
    assert fw.darcy_friction_haaland(5500, 1e-6, 0.1) == pytest.approx(
        0.0366668931278784
    )


def test_heat_transfer(fw):
    assert fw.heat_transfer(
        mflux_coolant=112.19853108876258,
        den_coolant=8.8673250601290707,
        radius_channel=0.0060000000000000001,
        heatcap_coolant=5184.9330299967578,
        visc_coolant=4.0416219836935569e-05,
        thermcond_coolant=0.3211653052986152,
        roughness_fw_channel=6e-8,
    ) == pytest.approx(1929.2042015869506)


def test_fw_thermal_conductivity(monkeypatch, fw):
    monkeypatch.setattr(fwbs_variables, "fw_th_conductivity", 28.9)

    assert fw.fw_thermal_conductivity(1900.0) == pytest.approx(326.70406785462256)
