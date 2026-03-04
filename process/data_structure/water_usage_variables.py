"""Module containing global variables relating to the water usage

References
----------
https://www.usgs.gov/special-topic/water-science-school/science/water-density
https://www.thermal-engineering.org/what-is-latent-heat-of-vaporization-definition/
"""

from dataclasses import dataclass


@dataclass
class WaterUseData:
    airtemp: float = 15.0
    """ambient air temperature (degrees Celsius)"""
    watertemp: float = 5.0
    """water temperature (degrees Celsius)"""
    windspeed: float = 4.0
    """wind speed (m/s)"""
    waterdens: float = 998.02
    """density of water (kg/m3)
    for simplicity, set to static value applicable to water at 21 degC
    """
    latentheat: float = 2257000.0
    """latent heat of vaporization (J/kg)
    for simplicity, set to static value applicable at 1 atm (100 kPa) air pressure
    """
    volheat: float = 0.0
    """volumetric heat of vaporization (J/m3)"""
    evapratio: float = 0.0
    """evaporation ratio: ratio of the heat used to evaporate water
    to the total heat discharged through the tower
    """

    evapvol: float = 0.0
    """evaporated volume of water (m3)"""

    energypervol: float = 0.0
    """input waste (heat) energy cooled per evaporated volume (J/m3)"""

    volperenergy: float = 0.0
    """volume evaporated by units of heat energy (m3/MJ)"""

    waterusetower: float = 0.0
    """total volume of water used in cooling tower (m3)"""

    wateruserecirc: float = 0.0
    """total volume of water used in recirculating system (m3)"""

    wateruseonethru: float = 0.0
    """total volume of water used in once-through system (m3)"""


# Another disgusting we may need to do in the transition period to support the dicts.
# Once all variables in the new data structure we can make the dicts from the DataStructure...
# and then in the long term put metadata on these classes and entierly remove the dicts
# In the meantime... the dicts will check each module for a '_CREATE_DICTS_FROM_DATACLASS' attribute
# and, if present, use this to create the dict ...
CREATE_DICTS_FROM_DATACLASS = WaterUseData
