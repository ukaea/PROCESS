"""Engineering models for pumping system analysis."""

import logging
from dataclasses import dataclass
from enum import IntEnum
from types import DynamicClassAttribute

import numpy as np

from process.core.coolprop_interface import FluidProperties
from process.core.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


class CoolantType(IntEnum):
    """Enum for coolant types."""

    HELIUM = (1, "Helium")
    WATER = (2, "Water")

    def __new__(cls, value: int, full_name: str):
        """Create a new CoolantType enum member with value and full_name.

        Parameters
        ----------
        value : int
            The numeric value of the enum member.
        full_name : str
            The full name description of the coolant type.
            This should match the CoolProp name for the coolant type.

        Returns
        -------
        CoolantType
            A new enum member with the specified value and full_name.
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._full_name_ = full_name
        return obj

    @DynamicClassAttribute
    def full_name(self):
        """The full name of the coolant type."""
        return self._full_name_


@dataclass
class CoolantFrictionLossParameters:
    """Parameters for calculating coolant friction losses."""

    dpres_total: float
    """Total pressure drop across the coolant channel (Pa)"""
    dpres_straight: float
    """Pressure drop due to straight length of the coolant channel (Pa)"""
    dpres_90: float
    """Pressure drop due to 90 degree bends in the coolant channel (Pa)"""
    dpres_90_total: float
    """Total pressure drop due to 90 degree bends in the coolant channel (Pa)"""
    dpres_180: float
    """Pressure drop due to 180 degree bends in the coolant channel (Pa)"""
    dpres_180_total: float
    """Total pressure drop due to 180 degree bends in the coolant channel (Pa)"""
    dpres_bends_total: float
    """Total pressure drop due to bends in the coolant channel (Pa)"""
    reynolds_number: float
    """Reynolds number of the coolant flow in the channel"""
    darcy_friction_factor: float
    """Darcy friction factor for the coolant flow in the channel"""
    f_straight: float
    """Friction factor for straight length of the coolant channel"""
    len_straight: float
    """Length of straight sections of the coolant channel (m)"""
    f_elbow_90: float
    """Friction factor for 90 degree bends in the coolant channel"""
    f_elbow_180: float
    """Friction factor for 180 degree bends in the coolant channel"""


def coolant_friction_pressure_drop(
    i_ps: int,
    radius_pipe_90_deg_bend: float,
    radius_pipe_180_deg_bend: float,
    n_pipe_90_deg_bends: float,
    n_pipe_180_deg_bends: float,
    len_pipe: float,
    den_coolant: float,
    visc_coolant: float,
    vel_coolant: float,
    roughness_channel: float,
    radius_channel: float,
    a_bz_liq: float,
    b_bz_liq: float,
) -> CoolantFrictionLossParameters:
    """Pressure drops are calculated for a pipe with a number of 90
    and 180 degree bends. The pressure drop due to frictional forces along
    the total straight length of the pipe is calculated, then the pressure
    drop due to the bends is calculated. The total pressure drop is the sum
    of all contributions.

    Parameters
    ----------
    i_ps :
        switch for primary or secondary coolant
    radius_pipe_90_deg_bend :
        radius of 90 degree bend in pipe [m]
    radius_pipe_180_deg_bend :
        radius of 180 degree bend in pipe [m]
    n_pipe_90_deg_bends :
        number of 90 degree bends in the pipe
    n_pipe_180_deg_bends :
        number of 180 degree bends in the pipe
    len_pipe :
        total flow length along pipe [m]
    den_coolant :
        coolant density [kg/m³]
    visc_coolant :
        coolant viscosity [Pa s]
    vel_coolant :
        coolant flow velocity [m/s]
    roughness_channel :
        roughness of the channel wall (ε) [m]
    radius_channel :
        radius of the channel [m]
    a_bz_liq :
        width of the breeding blanket coolant channel [m]
    b_bz_liq :
        height of the breeding blanket coolant channel [m]


    Returns
    -------
    :
        `CoolantFrictionLossParameters` dataclass containing:
        - Total pressure drop due to friction (Pa)
        - Pressure drop due to straight sections (Pa)
        - Pressure drop due to 90 degree bends (Pa)
        - Pressure drop due to 180 degree bends (Pa)
        - Reynolds number
        - Darcy friction factor
        - Pressure drop coefficient for straight sections
        - Pressure drop coefficient for 90 degree bends
        - Pressure drop coefficient for 180 degree bends

    Notes
    -----
        Darcy-Weisbach Equation (straight pipe):

        ΔP = λ * L/D * (p 〈v〉²) / 2

        λ - Darcy friction factor, L - pipe length, D - hydraulic diameter,
        p - fluid density, 〈v〉 - fluid flow average velocity

        This function also calculates pressure drop equations for elbow bends,
        with modified coefficients.

        N.B. Darcy friction factor is estimated from the Haaland approximation.
    """
    # Calculate hydraulic dimater for round or retancular pipe (m)
    dia_pipe = pipe_hydraulic_diameter(
        i_channel_shape=i_ps,
        radius_fw_channel=radius_channel,
        a_bz_liq=a_bz_liq,
        b_bz_liq=b_bz_liq,
    )

    # Reynolds number
    reynolds_number = calculate_reynolds_number(
        den_coolant=den_coolant,
        vel_coolant=vel_coolant,
        radius_channel=dia_pipe / 2,
        visc_coolant=visc_coolant,
    )

    # Calculate Darcy friction factor
    # N.B. friction function Uses Haaland approx. which assumes a filled
    # circular pipe.
    # Use dh which allows us to do fluid calculations for non-circular tubes
    # (dh is estimate appropriate for fully developed flow).

    darcy_friction_factor = darcy_friction_haaland(
        reynolds=reynolds_number,
        roughness_channel=roughness_channel,
        radius_channel=radius_channel,
    )

    # Pressure drop coefficient

    # Straight section
    f_straight = darcy_friction_factor * len_pipe / dia_pipe

    # 90 degree elbow pressure drop coefficient
    f_elbow_90 = elbow_coeff(
        radius_pipe_elbow=radius_pipe_90_deg_bend,
        deg_pipe_elbow=90.0,
        darcy_friction=darcy_friction_factor,
        dia_pipe=dia_pipe,
    )

    # 180 degree elbow pressure drop coefficient
    f_elbow_180 = elbow_coeff(
        radius_pipe_elbow=radius_pipe_180_deg_bend,
        deg_pipe_elbow=180.0,
        darcy_friction=darcy_friction_factor,
        dia_pipe=dia_pipe,
    )

    # Pressure drop due to friction in straight sections
    dpres_straight = f_straight * 0.5 * den_coolant * vel_coolant**2

    # Pressure drop due to 90 and 180 degree bends
    dpres_90 = f_elbow_90 * 0.5 * den_coolant * vel_coolant**2
    dpres_90_total = n_pipe_90_deg_bends * dpres_90
    dpres_180 = f_elbow_180 * 0.5 * den_coolant * vel_coolant**2
    dpres_180_total = n_pipe_180_deg_bends * dpres_180

    dpres_bends_total = dpres_90_total + dpres_180_total

    # Total pressure drop (Pa)
    dpres_total = dpres_straight + dpres_bends_total

    return CoolantFrictionLossParameters(
        dpres_total=dpres_total,
        dpres_straight=dpres_straight,
        dpres_90=dpres_90,
        dpres_90_total=dpres_90_total,
        dpres_180=dpres_180,
        dpres_180_total=dpres_180_total,
        dpres_bends_total=dpres_bends_total,
        reynolds_number=reynolds_number,
        darcy_friction_factor=darcy_friction_factor,
        f_straight=f_straight,
        len_straight=len_pipe,
        f_elbow_90=f_elbow_90,
        f_elbow_180=f_elbow_180,
    )


def pipe_hydraulic_diameter(
    i_channel_shape, radius_fw_channel: float, a_bz_liq: float, b_bz_liq: float
) -> float:
    """Caculate the hydraulic diameter (m) for a given coolant pipe size/shape.


    Parameters
    ----------
    i_channel_shape :
        switch for circular or rectangular channel crossection.
        Shape depends on whether primary or secondary coolant
    """
    # If primary coolant then circular channels assumed
    if i_channel_shape == 1:
        return 2.0 * radius_fw_channel

    # If secondary coolant then rectangular channels assumed
    if i_channel_shape == 2:
        return 2 * a_bz_liq * b_bz_liq / (a_bz_liq + b_bz_liq)

    raise ProcessValueError(f"i_channel_shape ={i_channel_shape} is an invalid option.")


@staticmethod
def darcy_friction_haaland(
    reynolds: float, roughness_channel: float, radius_channel: float
) -> float:
    """Calculate Darcy friction factor using the Haaland equation.

    Parameters
    ----------
    reynolds:
        Reynolds number.
    roughness_channel:
        Roughness of the first wall coolant channel (m).
    radius_channel:
        Radius of the first wall coolant channel (m).

    Returns
    -------
    :
        Darcy friction factor.

    Notes
    -----
    The Haaland equation is an approximation to the implicit Colebrook-White equation.
    It is used to calculate the Darcy friction factor for turbulent flow in pipes.

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation
    """
    # Bracketed term in Haaland equation
    bracket = (roughness_channel / radius_channel / 3.7) ** 1.11 + 6.9 / reynolds

    # Calculate Darcy friction factor
    return (1.8 * np.log10(bracket)) ** (-2)


def coolant_pumping_power(
    i_liquid_breeder: int,
    temp_coolant_pump_outlet: float,
    temp_coolant_pump_inlet: float,
    pres_coolant_pump_inlet: float,
    dpres_coolant: float,
    mflow_coolant_total: float,
    i_coolant_type: int,
    den_coolant: float,
    etaiso: float,
    etaiso_liq: float,
) -> float:
    """Calculate the coolant pumping power in MW for the first wall (FW) or breeding blanket (BZ) coolant.

    Parameters
    ----------
    i_liquid_breeder : int
        Switch for primary coolant or secondary coolant/breeder (1=primary He/H2O, 2=secondary PbLi/Li).
    temp_coolant_pump_outlet : float
        Pump outlet temperature (K).
    temp_coolant_pump_inlet : float
        Pump inlet temperature (K).
    pres_coolant_pump_inlet : float
        Pump inlet pressure (Pa).
    dpres_coolant : float
        Coolant pressure drop (Pa).
    mflow_coolant_total : float
        Total coolant mass flow rate in (kg/s).
    i_coolant_type : int
        Type of FW/blanket coolant (e.g., 1=Helium, 2=Water)
    den_coolant : float
        Density of coolant or liquid breeder (kg/m³).
    etaiso : float
        Isentropic efficiency of the pump for primary coolant.
    etaiso_liq : float
        Isentropic efficiency of the pump for secondary coolant/breeder.

    Returns
    -------
    float
        Pumping power in MW.

    References
    ----------
        - Idel'Cik, I. E. (1969), Memento des pertes de charge
        - S.P. Sukhatme (2005), A Textbook on Heat Transfer
    """
    # Pump outlet pressure (Pa)
    # The pump adds the pressure lost going through the coolant channels back
    pres_coolant_pump_outlet = pres_coolant_pump_inlet + dpres_coolant

    # Adiabatic index for helium or water
    gamma = (5 / 3) if i_coolant_type == CoolantType.HELIUM else (4 / 3)

    # If calculating for primary coolant
    if i_liquid_breeder == 1:
        # The pumping power is be calculated in the most general way,
        # using enthalpies before and after the pump.

        pump_outlet_fluid_properties = FluidProperties.of(
            fluid_name=CoolantType(i_coolant_type).full_name,
            temperature=temp_coolant_pump_outlet,
            pressure=pres_coolant_pump_outlet,
        )

        # Assume isentropic pump so that s1 = s2
        s1 = pump_outlet_fluid_properties.entropy

        # Get specific enthalpy at the outlet (J/kg) before pump using pressure and entropy s1
        pump_inlet_fluid_properties = FluidProperties.of(
            fluid_name=CoolantType(i_coolant_type).full_name,
            pressure=pres_coolant_pump_inlet,
            entropy=s1,
        )

        # Pumping power (MW) is given by enthalpy change, with a correction for
        # the isentropic efficiency of the pump.
        fp = (
            temp_coolant_pump_outlet
            * (
                1
                - (pres_coolant_pump_outlet / pres_coolant_pump_inlet)
                ** -((gamma - 1) / gamma)
            )
            / (etaiso * (temp_coolant_pump_inlet - temp_coolant_pump_outlet))
        )
        pumppower = (
            1e-6
            * mflow_coolant_total
            * (
                pump_outlet_fluid_properties.enthalpy
                - pump_inlet_fluid_properties.enthalpy
            )
            / etaiso
        ) / (1 - fp)

    # If calculating for secondary coolant/breeder...
    else:
        # Calculate specific volume
        spec_vol = 1 / den_coolant

        # Pumping power (MW) is given by pressure change, with a correction for
        # the isentropic efficiency of the pump.
        fp = (
            temp_coolant_pump_outlet
            * (
                1
                - (pres_coolant_pump_outlet / pres_coolant_pump_inlet)
                ** -((gamma - 1) / gamma)
            )
            / (etaiso_liq * (temp_coolant_pump_inlet - temp_coolant_pump_outlet))
        )
        pumppower = (
            1e-6 * mflow_coolant_total * spec_vol * dpres_coolant / etaiso_liq
        ) / (1 - fp)

    # Error for dpres_coolant too large
    if fp >= 1:
        raise ProcessValueError("Pressure drops in coolant are too large to be feasible")

    return pumppower


def gnielinski_heat_transfer_coefficient(
    mflux_coolant: float,
    den_coolant: float,
    radius_channel: float,
    heatcap_coolant: float,
    visc_coolant: float,
    thermcond_coolant: float,
    roughness_channel: float,
) -> float:
    """Calculate heat transfer coefficient using Gnielinski correlation.

    Parameters
    ----------
    mflux_coolant:
        Coolant mass flux in a single channel (kg/m²/s).
    den_coolant:
        Coolant density (average of inlet and outlet) (kg/m³).
    radius_channel:
        Coolant pipe radius (m).
    heatcap_coolant:
        Coolant specific heat capacity (average of inlet and outlet) (J/kg/K).
    visc_coolant:
        Coolant viscosity (average of inlet and outlet) (Pa.s).
    thermcond_coolant:
        Thermal conductivity of coolant (average of inlet and outlet) (W/m.K).
    roughness_channel:
        Roughness of the coolant channel (m).

    Returns
    -------
    :
        Heat transfer coefficient (W/m²K).

    Notes
    -----
    Gnielinski correlation. Ignore the distinction between wall and
    bulk temperatures. Valid for: 3000 < Re < 5e6, 0.5 < Pr < 2000

    References
    ----------
    [1] https://en.wikipedia.org/wiki/Nusselt_number#Gnielinski_correlation

    """
    # Calculate pipe diameter (m)
    diameter = 2 * radius_channel

    # Calculate flow velocity (m/s)
    vel_coolant = mflux_coolant / den_coolant

    # Calculate Reynolds number
    reynolds = calculate_reynolds_number(
        den_coolant=den_coolant,
        vel_coolant=vel_coolant,
        radius_channel=diameter / 2,
        visc_coolant=visc_coolant,
    )

    # Calculate Prandtl number
    pr = heatcap_coolant * visc_coolant / thermcond_coolant

    # Calculate Darcy friction factor, using Haaland equation
    f = darcy_friction_haaland(
        reynolds=reynolds,
        roughness_channel=roughness_channel,
        radius_channel=radius_channel,
    )

    # Calculate the Nusselt number
    nusselt = (
        (f / 8.0)
        * (reynolds - 1000.0)
        * pr
        / (1 + 12.7 * np.sqrt(f / 8.0) * (pr ** (2 / 3) - 1.0))
    )

    # Calculate the heat transfer coefficient (W/m^2K)
    heat_transfer_coefficient = nusselt * thermcond_coolant / (2.0 * radius_channel)

    # Check that Reynolds number is in valid range for the Gnielinski correlation
    if (reynolds <= 3000.0) or (reynolds > 5.0e6):
        logger.error("Reynolds number out of range : [3e3-5000e3]. %s", reynolds)

    # Check that Prandtl number is in valid range for the Gnielinski correlation
    if (pr < 0.5) or (pr > 2000.0):
        logger.error("Prandtl number out of range : [0.5-2000]. %s", pr)

    # Check that the Darcy friction factor is in valid range for the Gnielinski
    # correlation
    if f <= 0.0:
        logger.error("Negative Darcy friction factor (f). %s", f)

    return heat_transfer_coefficient


@staticmethod
def calculate_reynolds_number(
    den_coolant: float,
    vel_coolant: float,
    radius_channel: float,
    visc_coolant: float,
) -> float:
    """Calculate Reynolds number for flow in a pipe.

    Parameters
    ----------
    den_coolant:
        Coolant density (average of inlet and outlet) (kg/m³).
    vel_coolant:
        Coolant velocity in a single channel (m/s).
    radius_channel:
        Coolant pipe radius (m).
    visc_coolant:
        Coolant viscosity (average of inlet and outlet) (Pa.s).

    Returns
    -------
    :
        Reynolds number.

    """
    # Calculate pipe diameter (m)
    diameter = 2 * radius_channel

    # Calculate Reynolds number
    return den_coolant * vel_coolant * diameter / visc_coolant


@staticmethod
def elbow_coeff(
    radius_pipe_elbow: float,
    deg_pipe_elbow: float,
    darcy_friction: float,
    dia_pipe: float,
) -> float:
    """Calculates elbow bend coefficients for pressure drop calculations.

    Parameters
    ----------
    radius_pipe_elbow : float
        Pipe elbow radius (m)
    deg_pipe_elbow : float
        Pipe elbow angle (degrees)
    darcy_friction : float
        Darcy friction factor
    dia_pipe : float
        Pipe diameter (m)

    Returns
    -------
    float
        Elbow coefficient for pressure drop calculation

    References
    ----------
    [1] Idel'Cik, I. E. (1969), Memento des pertes de charge,
    Collection de la Direction des Etudes et Recherches d'Electricité de France.
    """
    if deg_pipe_elbow == 90:
        a = 1.0
    elif deg_pipe_elbow < 70:
        a = 0.9 * np.sin(deg_pipe_elbow * np.pi / 180.0)
    elif deg_pipe_elbow > 100:
        a = 0.7 + (0.35 * np.sin((deg_pipe_elbow / 90.0) * (np.pi / 180.0)))
    else:
        raise ProcessValueError(
            "No formula for 70 <= elbow angle(deg) <= 100, only 90 deg option available in this range."
        )

    r_ratio = radius_pipe_elbow / dia_pipe

    if r_ratio > 1:
        b = 0.21 / r_ratio**0.5
    elif r_ratio < 1:
        b = 0.21 / r_ratio**2.5
    else:
        b = 0.21

    # Singularity
    ximt = a * b

    # Friction
    xift = (
        (np.pi / 180.0)
        * darcy_friction
        * (radius_pipe_elbow / dia_pipe)
        * deg_pipe_elbow
    )

    # Elbow Coefficient
    return ximt + xift


@staticmethod
def calculate_required_mass_flow_rate(
    p_heat_total: float,
    heatcap_coolant: float,
    temp_in_coolant: float,
    temp_out_coolant: float,
) -> float:
    """Calculate the required mass flow rate of coolant to remove the specified heat
    load, due to the fundamental energy balance formula.

    Parameters
    ----------
    p_heat_total:
        Total heat load to be removed (W).
    heatcap_coolant:
        Specific heat capacity of the coolant (J/kg/K).
    temp_in_coolant:
        Inlet temperature of the coolant (K).
    temp_out_coolant:
        Outlet temperature of the coolant (K).

    Returns
    -------
    float
        Required mass flow rate of the coolant (kg/s).

    Notes
    -----
    The heat capacity is assumed to be constant over the temperature range of the
    coolant.

    """
    return p_heat_total / (heatcap_coolant * (temp_out_coolant - temp_in_coolant))
