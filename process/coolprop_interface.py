from CoolProp.CoolProp import PropsSI


class FluidProperties:
    def __init__(self, coolprop_inputs: list[str | float]):
        self._coolprop_inputs = coolprop_inputs

    @property
    def temperature(self):
        """fluid temperature [K]"""
        return PropsSI("T", *self._coolprop_inputs)

    @property
    def pressure(self):
        """fluid pressure [Pa]"""
        return PropsSI("P", *self._coolprop_inputs)

    @property
    def density(self):
        """fluid density [kg/m3]"""
        return PropsSI("D", *self._coolprop_inputs)

    @property
    def enthalpy(self):
        """fluid specific enthalpy [J/kg]"""
        return PropsSI("H", *self._coolprop_inputs)

    @property
    def entropy(self):
        """fluid entropy [J/kg/K]"""
        return PropsSI("S", *self._coolprop_inputs)

    @property
    def specific_heat_const_p(self):
        """fluid specific heat capacity at constant pressure [J/kg/K]"""
        return PropsSI("C", *self._coolprop_inputs)

    @property
    def specific_heat_const_v(self):
        """fluid specific heat capacity at constant volume [J/kg/K]"""
        return PropsSI("CVMASS", *self._coolprop_inputs)

    @property
    def viscosity(self):
        """fluid viscosity [Pa.s]"""
        return PropsSI("V", *self._coolprop_inputs)

    @property
    def thermal_conductivity(self):
        """fluid thermal conductivity [W/m/K]"""
        return PropsSI("CONDUCTIVITY", *self._coolprop_inputs)

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

        return cls(coolprop_inputs)
