"""Module containing variables for the water usage models

References
----------
https://www.usgs.gov/special-topic/water-science-school/science/water-density
https://www.thermal-engineering.org/what-is-latent-heat-of-vaporization-definition/
"""

from dataclasses import dataclass

from parameter_frame import Parameter

from process.core.metadata import PROCESSModelData


@dataclass(slots=True)
class WaterUseData(PROCESSModelData):
    """Dataclass holding water use variables"""

    airtemp: Parameter[float] = 15.0
    """ambient air temperature (degrees Celsius)"""
    watertemp: Parameter[float] = 5.0
    """water temperature (degrees Celsius)"""
    windspeed: Parameter[float] = 4.0
    """wind speed (m/s)"""
    waterdens: Parameter[float] = 998.02
    """density of water (kg/m3)
    for simplicity, set to static value applicable to water at 21 degC
    """
    latentheat: Parameter[float] = 2257000.0
    """latent heat of vaporization (J/kg)
    for simplicity, set to static value applicable at 1 atm (100 kPa) air pressure
    """
    volheat: Parameter[float] = 0.0
    """volumetric heat of vaporization (J/m3)"""
    evapratio: Parameter[float] = 0.0
    """evaporation ratio: ratio of the heat used to evaporate water
    to the total heat discharged through the tower
    """

    evapvol: Parameter[float] = 0.0
    """evaporated volume of water (m3/day)"""

    energypervol: Parameter[float] = 0.0
    """input waste (heat) energy cooled per evaporated volume (J/m3)"""

    volperenergy: Parameter[float] = 0.0
    """volume evaporated by units of heat energy (m3/MJ)"""

    waterusetower: Parameter[float] = 0.0
    """total volume of water used in cooling tower (m3/day)"""

    wateruserecirc: Parameter[float] = 0.0
    """total volume of water used in recirculating system (m3/day)"""

    wateruseonethru: Parameter[float] = 0.0
    """total volume of water used in once-through system (m3/day)"""


# Another disgusting we may need to do in the transition period to support the dicts.
# Once all variables in the new data structure we can make the dicts from the DataStructure...
# and then in the long term put metadata on these classes and entirely remove the dicts
# In the meantime... the dicts will check each module for a '_CREATE_DICTS_FROM_DATACLASS' attribute
# and, if present, use this to create the dict ...
CREATE_DICTS_FROM_DATACLASS = WaterUseData
