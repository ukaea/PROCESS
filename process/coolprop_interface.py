from functools import cache

from CoolProp.CoolProp import PropsSI


class FluidProperties:
    def __init__(self, coolprop_inputs: list[str | float]):
        self._coolprop_inputs = tuple(coolprop_inputs)

    @property
    def temperature(self):
        """fluid temperature [K]"""
        return _temperature(self._coolprop_inputs)

    @property
    def pressure(self):
        """fluid pressure [Pa]"""
        return _pressure(self._coolprop_inputs)

    @property
    def density(self):
        """fluid density [kg/m3]"""
        return _density(self._coolprop_inputs)

    @property
    def enthalpy(self):
        """fluid specific enthalpy [J/kg]"""
        return _enthalpy(self._coolprop_inputs)

    @property
    def entropy(self):
        """fluid entropy [J/kg/K]"""
        return _entropy(self._coolprop_inputs)

    @property
    def specific_heat_const_p(self):
        """fluid specific heat capacity at constant pressure [J/kg/K]"""
        return _specific_heat_const_p(self._coolprop_inputs)

    @property
    def specific_heat_const_v(self):
        """fluid specific heat capacity at constant volume [J/kg/K]"""
        return _specific_heat_const_v(self._coolprop_inputs)

    @property
    def viscosity(self):
        """fluid viscosity [Pa.s]"""
        return _viscosity(self._coolprop_inputs)

    @property
    def thermal_conductivity(self):
        """fluid thermal conductivity [W/m/K]"""
        return _thermal_conductivity(self._coolprop_inputs)

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

        Parameters
        ----------
        fluid_name :
            the name of the fluid to calculate properties for, e.g. 'Helium' or 'Water'.
        temperature :
            the current temperature [K] of the fluid to calculate the properties with respect to.
        pressure :
            the current pressure [Pa] of the fluid to calculate the properties with respect to.
        entropy :
            the current entropy [J/kg/K] of the fluid to calculate the properties with respect to.
        vapor_quality :
            the molar vapor quality [mol/mol] of the fluid to calculate the properties with respect to.
            `[0, 1]`, where `0` is a saturated liquid and `1` is a saturated vapor.
        fluid_name: str :



        temperature: float | None :
             (Default value = None)
        pressure: float | None :
             (Default value = None)
        entropy: float | None :
             (Default value = None)
        vapor_quality: float | None :
             (Default value = None)
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

        return cls(coolprop_inputs)


@cache
def _temperature(coolprop_inputs: tuple[str | float]):
    return PropsSI("T", *coolprop_inputs)


@cache
def _pressure(coolprop_inputs: tuple[str | float]):
    return PropsSI("P", *coolprop_inputs)


@cache
def _density(coolprop_inputs: tuple[str | float]):
    return PropsSI("D", *coolprop_inputs)


@cache
def _enthalpy(coolprop_inputs: tuple[str | float]):
    return PropsSI("H", *coolprop_inputs)


@cache
def _entropy(coolprop_inputs: tuple[str | float]):
    return PropsSI("S", *coolprop_inputs)


@cache
def _specific_heat_const_p(coolprop_inputs: tuple[str | float]):
    return PropsSI("C", *coolprop_inputs)


@cache
def _specific_heat_const_v(coolprop_inputs: tuple[str | float]):
    return PropsSI("CVMASS", *coolprop_inputs)


@cache
def _viscosity(coolprop_inputs: tuple[str | float]):
    return PropsSI("V", *coolprop_inputs)


@cache
def _thermal_conductivity(coolprop_inputs: tuple[str | float]):
    return PropsSI("CONDUCTIVITY", *coolprop_inputs)
