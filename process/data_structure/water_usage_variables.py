"""Module containing global variables relating to the water usage

References
----------
https://www.usgs.gov/special-topic/water-science-school/science/water-density
https://www.thermal-engineering.org/what-is-latent-heat-of-vaporization-definition/
"""

airtemp: float = None
"""ambient air temperature (degrees Celsius)"""

watertemp: float = None
"""water temperature (degrees Celsius)"""

windspeed: float = None
"""wind speed (m/s)"""

waterdens: float = None
"""density of water (kg/m3)
for simplicity, set to static value applicable to water at 21 degC
"""

latentheat: float = None
"""latent heat of vaporization (J/kg)
for simplicity, set to static value applicable at 1 atm (100 kPa) air pressure
"""

volheat: float = None
"""volumetric heat of vaporization (J/m3)"""

evapratio: float = None
"""evaporation ratio: ratio of the heat used to evaporate water
to the total heat discharged through the tower
"""

evapvol: float = None
"""evaporated volume of water (m3)"""

energypervol: float = None
"""input waste (heat) energy cooled per evaporated volume (J/m3)"""

volperenergy: float = None
"""volume evaporated by units of heat energy (m3/MJ)"""

waterusetower: float = None
"""total volume of water used in cooling tower (m3)"""

wateruserecirc: float = None
"""total volume of water used in recirculating system (m3)"""

wateruseonethru: float = None
"""total volume of water used in once-through system (m3)"""


def init_watuse_variables():
    """Initialise water variables"""
    global \
        airtemp, \
        watertemp, \
        windspeed, \
        waterdens, \
        latentheat, \
        volheat, \
        evapratio, \
        evapvol, \
        energypervol, \
        volperenergy, \
        waterusetower, \
        wateruserecirc, \
        wateruseonethru

    airtemp = 15.0
    watertemp = 5.0
    windspeed = 4.0
    waterdens = 998.02
    latentheat = 2257000.0
    volheat = 0.0
    evapratio = 0.0
    evapvol = 0.0
    energypervol = 0.0
    volperenergy = 0.0
    waterusetower = 0.0
    wateruserecirc = 0.0
    wateruseonethru = 0.0
