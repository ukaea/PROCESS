from dataclasses import dataclass

import pytest


@pytest.fixture
def water_use(process_models):
    """Provides WaterUse object for testing.

    :return water_use: initialised WaterUse object
    :type water_use: process.water_use.WaterUse
    """
    return process_models.water_use


@dataclass
class CoolingTowersParam:
    airtemp: float
    waterdens: float
    latentheat: float
    volheat: float
    evapratio: float
    evapvol: float
    energypervol: float
    volperenergy: float
    waterusetower: float
    wastetherm: float
    expected_volheat: float
    expected_evapratio: float
    expected_evapvol: float
    expected_energypervol: float
    expected_volperenergy: float
    expected_waterusetower: float


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
    for field in [
        "airtemp",
        "waterdens",
        "latentheat",
        "volheat",
        "evapratio",
        "evapvol",
        "energypervol",
        "volperenergy",
        "waterusetower",
    ]:
        monkeypatch.setattr(
            water_use.data.water_use, field, getattr(coolingtowersparam, field)
        )

    water_use.cooling_towers(wastetherm=coolingtowersparam.wastetherm, output=False)

    assert water_use.data.water_use.volheat == pytest.approx(
        coolingtowersparam.expected_volheat
    )
    assert water_use.data.water_use.evapratio == pytest.approx(
        coolingtowersparam.expected_evapratio
    )
    assert water_use.data.water_use.evapvol == pytest.approx(
        coolingtowersparam.expected_evapvol
    )
    assert water_use.data.water_use.energypervol == pytest.approx(
        coolingtowersparam.expected_energypervol
    )
    assert water_use.data.water_use.volperenergy == pytest.approx(
        coolingtowersparam.expected_volperenergy
    )
    assert water_use.data.water_use.waterusetower == pytest.approx(
        coolingtowersparam.expected_waterusetower
    )


@dataclass
class CoolingWaterBodyParam:
    watertemp: float
    windspeed: float
    waterdens: float
    latentheat: float
    volheat: float
    evapratio: float
    evapvol: float
    energypervol: float
    volperenergy: float
    wateruserecirc: float
    wateruseonethru: float
    wastetherm: float
    expected_evapratio: float
    expected_evapvol: float
    expected_energypervol: float
    expected_volperenergy: float
    expected_wateruserecirc: float
    expected_wateruseonethru: float


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
    for field in [
        "watertemp",
        "windspeed",
        "waterdens",
        "latentheat",
        "volheat",
        "evapratio",
        "energypervol",
        "energypervol",
        "volperenergy",
        "wateruserecirc",
        "wateruseonethru",
    ]:
        monkeypatch.setattr(
            water_use.data.water_use,
            field,
            getattr(coolingwaterbodyparam, field),
        )

    water_use.cooling_water_body(
        wastetherm=coolingwaterbodyparam.wastetherm, output=False
    )

    assert water_use.data.water_use.evapratio == pytest.approx(
        coolingwaterbodyparam.expected_evapratio
    )
    assert water_use.data.water_use.evapvol == pytest.approx(
        coolingwaterbodyparam.expected_evapvol
    )
    assert water_use.data.water_use.energypervol == pytest.approx(
        coolingwaterbodyparam.expected_energypervol
    )
    assert water_use.data.water_use.volperenergy == pytest.approx(
        coolingwaterbodyparam.expected_volperenergy
    )
    assert water_use.data.water_use.wateruserecirc == pytest.approx(
        coolingwaterbodyparam.expected_wateruserecirc
    )
    assert water_use.data.water_use.wateruseonethru == pytest.approx(
        coolingwaterbodyparam.expected_wateruseonethru
    )
