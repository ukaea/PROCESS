import dataclasses

from CoolProp.CoolProp import PropsSI


@dataclasses.dataclass
class FluidProperties:
    temperature: float
    """fluid temperature [K]"""
    pressure: float
    """fluid pressure [Pa]"""
    density: float
    """fluid density [kg/m3]"""
    enthalpy: float
    """fluid specific enthalpy [J/kg]"""
    entropy: float
    """fluid entropy [J/kg/K]"""
    specific_heat_const_p: float
    """fluid specific heat capacity at constant pressure [J/kg/K]"""
    specific_heat_const_v: float
    """fluid specific heat capacity at constant volume [J/kg/K]"""
    viscosity: float
    """fluid viscosity [Pa.s]"""
    thermal_conductivity: float
    """fluid thermal conductivity [W/m/K]"""

    @classmethod
    def of(
        cls,
        fluid_name: str,
        *,
        temperature: float | None = None,
        pressure: float | None = None,
        entropy: float | None = None,
        vapor_quality: float | None = None,
    ):
        """Calculates the fluid properties of a fluid given its temperature and pressure.

        :param fluid_name: the name of the fluid to calculate properties for, e.g. 'Helium' or 'Water'.
        :param temperature: the current temperature [K] of the fluid to calculate the properties with respect to.
        :param pressure: the current pressure [Pa] of the fluid to calculate the properties with respect to.
        :param entropy: the current entropy [J/kg/K] of the fluid to calculate the properties with respect to.
        :param vapor_quality: the molar vapor quality [mol/mol] of the fluid to calculate the properties with respect to.
        `[0, 1]`, where `0` is a saturated liquid and `1` is a saturated vapor.
        """
        coolprop_inputs = []

        if temperature is not None:
            coolprop_inputs += ["T", temperature]

        if pressure is not None:
            coolprop_inputs += ["P", pressure]

        if entropy is not None:
            coolprop_inputs += ["S", entropy]

        if vapor_quality is not None:
            coolprop_inputs += ["Q", vapor_quality]

        coolprop_inputs.append(fluid_name.title())

        return cls(
            temperature=PropsSI("T", *coolprop_inputs),
            pressure=PropsSI("P", *coolprop_inputs),
            density=PropsSI("D", *coolprop_inputs),
            enthalpy=PropsSI("H", *coolprop_inputs),
            entropy=PropsSI("S", *coolprop_inputs),
            specific_heat_const_p=PropsSI("C", *coolprop_inputs),
            specific_heat_const_v=PropsSI("CVMASS", *coolprop_inputs),
            viscosity=PropsSI("V", *coolprop_inputs),
            thermal_conductivity=PropsSI("CONDUCTIVITY", *coolprop_inputs),
        )
