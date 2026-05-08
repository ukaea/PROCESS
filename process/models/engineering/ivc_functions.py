"""Module containing functions common to multiple in vessel components"""

import numpy as np


def calculate_pipe_bend_radius(i_ps: int, radius_fw_channel: float, b_bz_liq: float):
    """Set the pipe bend radius based on the coolant type.

    Parameters
    ----------
    i_ps : int
        switch for primary or secondary coolant
    radius_fw_channel: float
        radius of first wall cooling channels [m]
    b_bz_liq: float
        radial width of the rectangular cooling channel [m]
        for long poloidal sections of blanket breeding zone

    """
    # If primary coolant or secondary coolant (See DCLL)
    radius_pipe_90_deg_bend = (3 * radius_fw_channel) if (i_ps == 1) else b_bz_liq
    radius_pipe_180_deg_bend = radius_pipe_90_deg_bend / 2

    return radius_pipe_90_deg_bend, radius_pipe_180_deg_bend


def pumping_powers_as_fractions(
    f_p_fw_coolant_pump_total_heat: float,
    f_p_blkt_coolant_pump_total_heat: float,
    f_p_shld_coolant_pump_total_heat: float,
    f_p_div_coolant_pump_total_heat: float,
    p_fw_nuclear_heat_total_mw: float,
    psurffwi: float,
    psurffwo: float,
    p_blkt_nuclear_heat_total_mw: float,
    p_shld_nuclear_heat_mw: float,
    p_cp_shield_nuclear_heat_mw: float,
    p_plasma_separatrix_mw: float,
    p_div_nuclear_heat_total_mw: float,
    p_div_rad_total_mw: float,
) -> tuple[float, float, float, float]:
    """Calculate mechanical pumping powers as fractions of thermal power in each component.

    Parameters
    ----------
    f_p_fw_coolant_pump_total_heat : float
        Fraction for FW coolant pump.
    f_p_blkt_coolant_pump_total_heat : float
        Fraction for blanket coolant pump.
    f_p_shld_coolant_pump_total_heat : float
        Fraction for shield coolant pump.
    f_p_div_coolant_pump_total_heat : float
        Fraction for divertor coolant pump.
    p_fw_nuclear_heat_total_mw : float
        Total FW nuclear heating (MW).
    psurffwi : float
        Inboard FW surface heating (MW).
    psurffwo : float
        Outboard FW surface heating (MW).
    p_blkt_nuclear_heat_total_mw : float
        Total blanket nuclear heating (MW).
    p_shld_nuclear_heat_mw : float
        Shield nuclear heating (MW).
    p_cp_shield_nuclear_heat_mw : float
        CP shield nuclear heating (MW).
    p_plasma_separatrix_mw : float
        Plasma separatrix power (MW).
    p_div_nuclear_heat_total_mw : float
        Divertor nuclear heating (MW).
    p_div_rad_total_mw : float
        Divertor radiative power (MW).

    Returns
    -------
    tuple[float, float, float, float]
        Tuple of pumping powers (MW) for FW, blanket, shield, and divertor.
    """
    p_fw_coolant_pump_mw = f_p_fw_coolant_pump_total_heat * (
        p_fw_nuclear_heat_total_mw + psurffwi + psurffwo
    )
    p_blkt_coolant_pump_mw = (
        f_p_blkt_coolant_pump_total_heat * p_blkt_nuclear_heat_total_mw
    )
    p_shld_coolant_pump_mw = f_p_shld_coolant_pump_total_heat * (
        p_shld_nuclear_heat_mw + p_cp_shield_nuclear_heat_mw
    )
    p_div_coolant_pump_mw = f_p_div_coolant_pump_total_heat * (
        p_plasma_separatrix_mw + p_div_nuclear_heat_total_mw + p_div_rad_total_mw
    )
    return (
        p_fw_coolant_pump_mw,
        p_blkt_coolant_pump_mw,
        p_shld_coolant_pump_mw,
        p_div_coolant_pump_mw,
    )


def eshellarea(rshell, rmini, rmino, zminor):
    """Routine to calculate the inboard, outboard and total surface areas
    of a toroidal shell comprising two elliptical sections

    Parameters
    ----------
    rshell : float
        major radius of centre of both ellipses [m]
    rmini : float
        horizontal distance from rshell to inboard elliptical shell [m]
    rmino : float
        horizontal distance from rshell to outboard elliptical shell [m]
    zminor : float
        vertical internal half-height of shell [m]

    Returns
    -------
    tuple[float, float, float]
        Tuple containing:
        - ain: Surface area of inboard straight section (m²)
        - aout: Surface area of outboard curved section (m²)
        - atot: Total surface area of shell (m²)
    """
    # Inboard section
    elong = zminor / rmini
    ain = 2.0 * np.pi * elong * (np.pi * rshell * rmini - 2.0 * rmini * rmini)

    # Outboard section
    elong = zminor / rmino
    aout = 2.0 * np.pi * elong * (np.pi * rshell * rmino + 2.0 * rmino * rmino)

    return ain, aout, ain + aout


def dshellarea(
    rmajor: float, rminor: float, zminor: float
) -> tuple[float, float, float]:
    """Calculate the inboard, outboard, and total surface areas of a D-shaped toroidal shell.

    The inboard section is assumed to be a cylinder.
    The outboard section is defined by a semi-ellipse, centred on the major radius of the inboard section.

    Parameters
    ----------
    rmajor : float
        Major radius of inboard straight section (m)
    rminor : float
        Horizontal width of shell (m)
    zminor : float
        Vertical half-height of shell (m)

    Returns
    -------
    tuple[float, float, float]
        Tuple containing:
        - ain: Surface area of inboard straight section (m²)
        - aout: Surface area of outboard curved section (m²)
        - atot: Total surface area of shell (m²)
    """
    # Area of inboard cylindrical shell
    ain = 4.0 * zminor * np.pi * rmajor

    # Area of elliptical outboard section
    elong = zminor / rminor
    aout = 2.0 * np.pi * elong * (np.pi * rmajor * rminor + 2.0 * rminor * rminor)

    return ain, aout, ain + aout


def eshellvol(rshell, rmini, rmino, zminor, drin, drout, dz):
    """Routine to calculate the inboard, outboard and total volumes
    of a toroidal shell comprising two elliptical sections

    This routine calculates the volume of the inboard and outboard sections
    of a toroidal shell defined by two co-centred semi-ellipses.
    Each section's internal and external surfaces are in turn defined
    by two semi-ellipses. The volumes of each section are calculated as
    the difference in those of the volumes of revolution enclosed by their
    inner and outer surfaces.

    Parameters
    ----------
    rshell :
        major radius of centre of both ellipses (m)
    rmini :
        horizontal distance from rshell to outer edge of inboard elliptical shell (m)
    rmino :
        horizontal distance from rshell to inner edge of outboard elliptical shell (m)
    zminor :
        vertical internal half-height of shell (m)
    drin :
        horiz. thickness of inboard shell at midplane (m)
    drout :
        horiz. thickness of outboard shell at midplane (m)
    dz :
        vertical thickness of shell at top/bottom (m)

    Returns
    -------
    vin :
        volume of inboard section (m3)
    vout:
        volume of outboard section (m3)
    vtot:
        total volume of shell (m3)
    -------
    """
    # Inboard section

    # Volume enclosed by outer (higher R) surface of elliptical section
    # and the vertical straight line joining its ends
    a = rmini
    b = zminor
    elong = b / a
    v1 = 2.0 * np.pi * elong * (0.5 * np.pi * rshell * a**2 - (2.0 / 3.0) * a**3)

    # Volume enclosed by inner (lower R) surface of elliptical section
    # and the vertical straight line joining its ends
    a = rmini + drin
    b = zminor + dz
    elong = b / a
    v2 = 2.0 * np.pi * elong * (0.5 * np.pi * rshell * a**2 - (2.0 / 3.0) * a**3)

    # Volume of inboard section of shell
    vin = v2 - v1

    # Outboard section

    # Volume enclosed by inner (lower R) surface of elliptical section
    # and the vertical straight line joining its ends
    a = rmino
    b = zminor
    elong = b / a
    v1 = 2.0 * np.pi * elong * (0.5 * np.pi * rshell * a**2 + (2.0 / 3.0) * a**3)

    # Volume enclosed by outer (higher R) surface of elliptical section
    # and the vertical straight line joining its ends
    a = rmino + drout
    b = zminor + dz
    elong = b / a
    v2 = 2.0 * np.pi * elong * (0.5 * np.pi * rshell * a**2 + (2.0 / 3.0) * a**3)

    # Volume of outboard section of shell
    vout = v2 - v1

    return vin, vout, vin + vout


def dshellvol(rmajor, rminor, zminor, drin, drout, dz):
    """Routine to calculate the inboard, outboard and total volumes
    of a D-shaped toroidal shell

    This routine calculates the volume of the inboard and outboard sections
    of a D-shaped toroidal shell defined by the above input parameters.
    The inboard section is assumed to be a cylinder of uniform thickness.
    The outboard section's internal and external surfaces are defined
    by two semi-ellipses, centred on the outer edge of the inboard section;
    its volume is calculated as the difference in those of the volumes of
    revolution enclosed by the two surfaces.

    Parameters
    ----------
    rmajor :
        major radius to outer point of inboardstraight section of shell (m)
    rminor :
        horizontal internal width of shell (m)
    zminor :
        vertical internal half-height of shell (m)
    drin :
        horiz. thickness of inboard shell at midplane (m)
    drout :
        horiz. thickness of outboard shell at midplane (m)
    dz :
        vertical thickness of shell at top/bottom (m)

    Returns
    -------
    vin :
        volume of inboard straight section (m3)
    vout:
        volume of outboard curved section (m3)
    vtot:
        total volume of shell (m3)

    """
    # Volume of inboard cylindrical shell
    vin = 2.0 * (zminor + dz) * np.pi * (rmajor**2 - (rmajor - drin) ** 2)

    # Volume enclosed by inner surface of elliptical outboard section
    # and the vertical straight line joining its ends
    a = rminor
    b = zminor
    elong = b / a
    v1 = 2.0 * np.pi * elong * (0.5 * np.pi * rmajor * a**2 + (2.0 / 3.0) * a**3)

    # Volume enclosed by outer surface of elliptical outboard section
    # and the vertical straight line joining its ends
    a = rminor + drout
    b = zminor + dz
    elong = b / a
    v2 = 2.0 * np.pi * elong * (0.5 * np.pi * rmajor * a**2 + (2.0 / 3.0) * a**3)

    # Volume of elliptical outboard shell
    vout = v2 - v1

    return vin, vout, vin + vout
