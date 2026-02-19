from typing import Any, NamedTuple

import pytest

from process.data_structure import water_usage_variables
from process.models.water_use import WaterUse


@pytest.fixture
def water_use():
    """Provides WaterUse object for testing.

    :return water_use: initialised WaterUse object
    :type water_use: process.water_use.WaterUse
    """
    return WaterUse()


class CoolingTowersParam(NamedTuple):
    airtemp: Any
    waterdens: Any
    latentheat: Any
    volheat: Any
    evapratio: Any
    evapvol: Any
    energypervol: Any
    volperenergy: Any
    waterusetower: Any
    outfile: Any
    iprint: Any
    wastetherm: Any
    expected_volheat: Any
    expected_evapratio: Any
    expected_evapvol: Any
    expected_energypervol: Any
    expected_volperenergy: Any
    expected_waterusetower: Any


@pytest.mark.parametrize(
    "coolingtowersparam",
    (
        CoolingTowersParam(
            airtemp=15,
            waterdens=998.01999999999998,
            latentheat=2257000,
            volheat=0,
            evapratio=0,
            evapvol=0,
            energypervol=0,
            volperenergy=0,
            waterusetower=0,
            outfile=11,
            iprint=0,
            wastetherm=141491977.80211401,
            expected_volheat=2252531140,
            expected_evapratio=0.79171374999999999,
            expected_evapvol=49731.230059988622,
            expected_energypervol=2845133282.0732241,
            expected_volperenergy=0.000351477382905348,
            expected_waterusetower=69623.722083984059,
        ),
        CoolingTowersParam(
            airtemp=15,
            waterdens=998.01999999999998,
            latentheat=2257000,
            volheat=2252531140,
            evapratio=0.37103027579005016,
            evapvol=23306.140640523186,
            energypervol=6071017075.9073286,
            volperenergy=0.00016471704617146833,
            waterusetower=69623.722083984059,
            outfile=11,
            iprint=0,
            wastetherm=141448808.82309783,
            expected_volheat=2252531140,
            expected_evapratio=0.79171374999999999,
            expected_evapvol=49716.057140221317,
            expected_energypervol=2845133282.0732241,
            expected_volperenergy=0.000351477382905348,
            expected_waterusetower=69602.479996309834,
        ),
    ),
)
def test_cooling_towers(coolingtowersparam, monkeypatch, water_use):
    monkeypatch.setattr(water_usage_variables, "airtemp", coolingtowersparam.airtemp)
    monkeypatch.setattr(water_usage_variables, "waterdens", coolingtowersparam.waterdens)
    monkeypatch.setattr(
        water_usage_variables, "latentheat", coolingtowersparam.latentheat
    )
    monkeypatch.setattr(water_usage_variables, "volheat", coolingtowersparam.volheat)
    monkeypatch.setattr(water_usage_variables, "evapratio", coolingtowersparam.evapratio)
    monkeypatch.setattr(water_usage_variables, "evapvol", coolingtowersparam.evapvol)
    monkeypatch.setattr(
        water_usage_variables, "energypervol", coolingtowersparam.energypervol
    )
    monkeypatch.setattr(
        water_usage_variables, "volperenergy", coolingtowersparam.volperenergy
    )
    monkeypatch.setattr(
        water_usage_variables, "waterusetower", coolingtowersparam.waterusetower
    )

    water_use.cooling_towers(wastetherm=coolingtowersparam.wastetherm, output=False)

    assert water_usage_variables.volheat == pytest.approx(
        coolingtowersparam.expected_volheat
    )
    assert water_usage_variables.evapratio == pytest.approx(
        coolingtowersparam.expected_evapratio
    )
    assert water_usage_variables.evapvol == pytest.approx(
        coolingtowersparam.expected_evapvol
    )
    assert water_usage_variables.energypervol == pytest.approx(
        coolingtowersparam.expected_energypervol
    )
    assert water_usage_variables.volperenergy == pytest.approx(
        coolingtowersparam.expected_volperenergy
    )
    assert water_usage_variables.waterusetower == pytest.approx(
        coolingtowersparam.expected_waterusetower
    )


class CoolingWaterBodyParam(NamedTuple):
    watertemp: Any
    windspeed: Any
    waterdens: Any
    latentheat: Any
    volheat: Any
    evapratio: Any
    evapvol: Any
    energypervol: Any
    volperenergy: Any
    wateruserecirc: Any
    wateruseonethru: Any
    outfile: Any
    iprint: Any
    wastetherm: Any
    expected_evapratio: Any
    expected_evapvol: Any
    expected_energypervol: Any
    expected_volperenergy: Any
    expected_wateruserecirc: Any
    expected_wateruseonethru: Any


@pytest.mark.parametrize(
    "coolingwaterbodyparam",
    (
        CoolingWaterBodyParam(
            watertemp=5,
            windspeed=4,
            waterdens=998.01999999999998,
            latentheat=2257000,
            volheat=2252531140,
            evapratio=0.79171374999999999,
            evapvol=49731.230059988622,
            energypervol=2845133282.0732241,
            volperenergy=0.000351477382905348,
            wateruserecirc=0,
            wateruseonethru=0,
            outfile=11,
            iprint=0,
            wastetherm=141491977.80211401,
            expected_evapratio=0.37103027579005016,
            expected_evapvol=23306.140640523186,
            expected_energypervol=6071017075.9073286,
            expected_volperenergy=0.00016471704617146833,
            expected_wateruserecirc=23360.25002895,
            expected_wateruseonethru=2289304.50283738,
        ),
        CoolingWaterBodyParam(
            watertemp=5,
            windspeed=4,
            waterdens=998.01999999999998,
            latentheat=2257000,
            volheat=2252531140,
            evapratio=0.79171374999999999,
            evapvol=49716.057140221317,
            energypervol=2845133282.0732241,
            volperenergy=0.000351477382905348,
            wateruserecirc=23306.140640523186,
            wateruseonethru=2284001.7827712721,
            outfile=11,
            iprint=0,
            wastetherm=141448808.82309783,
            expected_evapratio=0.37103027579005016,
            expected_evapvol=23299.029973813402,
            expected_energypervol=6071017075.9073286,
            expected_volperenergy=0.00016471704617146833,
            expected_wateruserecirc=23353.12285355,
            expected_wateruseonethru=2288606.0396483,
        ),
    ),
)
def test_cooling_water_body(coolingwaterbodyparam, monkeypatch, water_use):
    monkeypatch.setattr(
        water_usage_variables, "watertemp", coolingwaterbodyparam.watertemp
    )
    monkeypatch.setattr(
        water_usage_variables, "windspeed", coolingwaterbodyparam.windspeed
    )
    monkeypatch.setattr(
        water_usage_variables, "waterdens", coolingwaterbodyparam.waterdens
    )
    monkeypatch.setattr(
        water_usage_variables, "latentheat", coolingwaterbodyparam.latentheat
    )
    monkeypatch.setattr(water_usage_variables, "volheat", coolingwaterbodyparam.volheat)
    monkeypatch.setattr(
        water_usage_variables, "evapratio", coolingwaterbodyparam.evapratio
    )
    monkeypatch.setattr(water_usage_variables, "evapvol", coolingwaterbodyparam.evapvol)
    monkeypatch.setattr(
        water_usage_variables, "energypervol", coolingwaterbodyparam.energypervol
    )
    monkeypatch.setattr(
        water_usage_variables, "volperenergy", coolingwaterbodyparam.volperenergy
    )
    monkeypatch.setattr(
        water_usage_variables, "wateruserecirc", coolingwaterbodyparam.wateruserecirc
    )
    monkeypatch.setattr(
        water_usage_variables, "wateruseonethru", coolingwaterbodyparam.wateruseonethru
    )

    water_use.cooling_water_body(
        wastetherm=coolingwaterbodyparam.wastetherm, output=False
    )

    assert water_usage_variables.evapratio == pytest.approx(
        coolingwaterbodyparam.expected_evapratio
    )
    assert water_usage_variables.evapvol == pytest.approx(
        coolingwaterbodyparam.expected_evapvol
    )
    assert water_usage_variables.energypervol == pytest.approx(
        coolingwaterbodyparam.expected_energypervol
    )
    assert water_usage_variables.volperenergy == pytest.approx(
        coolingwaterbodyparam.expected_volperenergy
    )
    assert water_usage_variables.wateruserecirc == pytest.approx(
        coolingwaterbodyparam.expected_wateruserecirc
    )
    assert water_usage_variables.wateruseonethru == pytest.approx(
        coolingwaterbodyparam.expected_wateruseonethru
    )
