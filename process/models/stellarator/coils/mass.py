"""Module for coil mass calculations in stellarators."""

from process.core import constants
from process.data_structure import fwbs_variables, tfcoil_variables


def calculate_coils_mass(a_tf_wp_with_insulation: float, a_tf_wp_no_insulation: float):
    """Calculates the mass of stellarator coils by aggregating the masses of various coil components.

    This function computes the masses of conductor constituents (casing, ground insulation, superconductor, copper),
    conduit masses (steel and insulation), and then calculates the total conductor and coil masses.

    Parameters
    ----------
    a_tf_wp_with_insulation :
        Area of the toroidal field coil winding pack with insulation.
    a_tf_wp_no_insulation:
        Area of the toroidal field coil winding pack without insulation.

    Returns
    -------
    :
        The function performs calculations and updates external state.

    """

    #  Masses of conductor constituents
    casing()
    ground_insulation(a_tf_wp_with_insulation, a_tf_wp_no_insulation)
    superconductor()
    copper()

    # conduit masses
    conduit_steel()
    conduit_insulation()

    # Total masses
    total_conductor()
    total_coil()


def casing():
    """[kg] Mass of case
    (no need for correction factors as is the case for tokamaks)
    This is only correct if the winding pack is 'thin' (len_tf_coil>>sqrt(tfcoil_variables.a_tf_coil_inboard_case)).

    """
    tfcoil_variables.m_tf_coil_case = (
        tfcoil_variables.len_tf_coil
        * tfcoil_variables.a_tf_coil_inboard_case
        * tfcoil_variables.den_tf_coil_case
    )


def ground_insulation(a_tf_wp_with_insulation, a_tf_wp_no_insulation):
    """Mass of ground-wall insulation [kg]
    (assumed to be same density/material as conduit insulation)

    Parameters
    ----------
    a_tf_wp_with_insulation :

    a_tf_wp_no_insulation :

    """
    tfcoil_variables.m_tf_coil_wp_insulation = (
        tfcoil_variables.len_tf_coil
        * (a_tf_wp_with_insulation - a_tf_wp_no_insulation)
        * tfcoil_variables.den_tf_wp_turn_insulation
    )


def superconductor():
    """[kg] mass of Superconductor
    a_tf_wp_coolant_channels is 0 for a stellarator. but keep this term for now.

    """
    tfcoil_variables.m_tf_coil_superconductor = (
        tfcoil_variables.len_tf_coil
        * tfcoil_variables.n_tf_coil_turns
        * tfcoil_variables.a_tf_turn_cable_space_no_void
        * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
        * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_copper)
        - tfcoil_variables.len_tf_coil * tfcoil_variables.a_tf_wp_coolant_channels
    ) * tfcoil_variables.dcond[tfcoil_variables.i_tf_sc_mat - 1]


def copper():
    """[kg] mass of Copper in conductor"""
    tfcoil_variables.m_tf_coil_copper = (
        tfcoil_variables.len_tf_coil
        * tfcoil_variables.n_tf_coil_turns
        * tfcoil_variables.a_tf_turn_cable_space_no_void
        * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
        * tfcoil_variables.f_a_tf_turn_cable_copper
        - tfcoil_variables.len_tf_coil * tfcoil_variables.a_tf_wp_coolant_channels
    ) * constants.den_copper


def conduit_steel():
    """[kg] mass of Steel conduit (sheath)"""
    tfcoil_variables.m_tf_wp_steel_conduit = (
        tfcoil_variables.len_tf_coil
        * tfcoil_variables.n_tf_coil_turns
        * tfcoil_variables.a_tf_turn_steel
        * fwbs_variables.den_steel
    )
    # if (i_tf_sc_mat==6)   tfcoil_variables.m_tf_wp_steel_conduit = fcondsteel * a_tf_wp_no_insulation *tfcoil_variables.len_tf_coil* fwbs_variables.denstl


def conduit_insulation():
    """Conduit insulation mass [kg]
    (tfcoil_variables.a_tf_coil_wp_turn_insulation already contains tfcoil_variables.n_tf_coil_turns)

    """
    tfcoil_variables.m_tf_coil_wp_turn_insulation = (
        tfcoil_variables.len_tf_coil
        * tfcoil_variables.a_tf_coil_wp_turn_insulation
        * tfcoil_variables.den_tf_wp_turn_insulation
    )


def total_conductor():
    """[kg] Total conductor mass"""
    tfcoil_variables.m_tf_coil_conductor = (
        tfcoil_variables.m_tf_coil_superconductor
        + tfcoil_variables.m_tf_coil_copper
        + tfcoil_variables.m_tf_wp_steel_conduit
        + tfcoil_variables.m_tf_coil_wp_turn_insulation
    )


def total_coil():
    """[kg] Total coil mass"""
    tfcoil_variables.m_tf_coils_total = (
        tfcoil_variables.m_tf_coil_case
        + tfcoil_variables.m_tf_coil_conductor
        + tfcoil_variables.m_tf_coil_wp_insulation
    ) * tfcoil_variables.n_tf_coils
