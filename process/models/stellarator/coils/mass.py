"""Module for coil mass calculations in stellarators."""

from process.core import constants
from process.core.model import DataStructure


def calculate_coils_mass(
    a_tf_wp_with_insulation: float, a_tf_wp_no_insulation: float, data: DataStructure
):
    """Calculates the mass of stellarator coils by aggregating the masses of
    various coil components.

    This function computes the masses of conductor constituents
    (casing, ground insulation, superconductor, copper),
    conduit masses (steel and insulation),
    and then calculates the total conductor and coil masses.

    Parameters
    ----------
    a_tf_wp_with_insulation :
        Area of the toroidal field coil winding pack with insulation.
    a_tf_wp_no_insulation:
        Area of the toroidal field coil winding pack without insulation.
    data: DataStructure
        data structure object to provide model data
    """
    #  Masses of conductor constituents
    casing(data)
    ground_insulation(a_tf_wp_with_insulation, a_tf_wp_no_insulation, data)
    superconductor(data)
    copper(data)

    # conduit masses
    conduit_steel(data)
    conduit_insulation(data)

    # Total masses
    total_conductor(data)
    total_coil(data)


def casing(data: DataStructure):
    """[kg] Mass of case
    (no need for correction factors as is the case for tokamaks)
    This is only correct if the winding pack is 'thin'
    (len_tf_coil>>sqrt(data.tfcoil.a_tf_coil_inboard_case)).

    """
    data.tfcoil.m_tf_coil_case = (
        data.tfcoil.len_tf_coil
        * data.tfcoil.a_tf_coil_inboard_case
        * data.tfcoil.den_tf_coil_case
    )


def ground_insulation(
    a_tf_wp_with_insulation, a_tf_wp_no_insulation, data: DataStructure
):
    """Mass of ground-wall insulation [kg]
    (assumed to be same density/material as conduit insulation)

    Parameters
    ----------
    a_tf_wp_with_insulation :

    a_tf_wp_no_insulation :

    """
    data.tfcoil.m_tf_coil_wp_insulation = (
        data.tfcoil.len_tf_coil
        * (a_tf_wp_with_insulation - a_tf_wp_no_insulation)
        * data.tfcoil.den_tf_wp_turn_insulation
    )


def superconductor(data: DataStructure):
    """[kg] mass of Superconductor
    a_tf_wp_coolant_channels is 0 for a stellarator. but keep this term for now.

    """
    data.tfcoil.m_tf_coil_superconductor = (
        data.tfcoil.len_tf_coil
        * data.tfcoil.n_tf_coil_turns
        * data.tfcoil.a_tf_turn_cable_space_no_void
        * (1.0e0 - data.tfcoil.f_a_tf_turn_cable_space_extra_void)
        * (1.0e0 - data.tfcoil.f_a_tf_turn_cable_copper)
        - data.tfcoil.len_tf_coil * data.tfcoil.a_tf_wp_coolant_channels
    ) * data.tfcoil.dcond[data.tfcoil.i_tf_sc_mat - 1]


def copper(data: DataStructure):
    """[kg] mass of Copper in conductor"""
    data.tfcoil.m_tf_coil_copper = (
        data.tfcoil.len_tf_coil
        * data.tfcoil.n_tf_coil_turns
        * data.tfcoil.a_tf_turn_cable_space_no_void
        * (1.0e0 - data.tfcoil.f_a_tf_turn_cable_space_extra_void)
        * data.tfcoil.f_a_tf_turn_cable_copper
        - data.tfcoil.len_tf_coil * data.tfcoil.a_tf_wp_coolant_channels
    ) * constants.DEN_COPPER


def conduit_steel(data):
    """[kg] mass of Steel conduit (sheath)"""
    data.tfcoil.m_tf_wp_steel_conduit = (
        data.tfcoil.len_tf_coil
        * data.tfcoil.n_tf_coil_turns
        * data.tfcoil.a_tf_turn_steel
        * data.fwbs.den_steel
    )
    # if (i_tf_sc_mat==6)
    # data.tfcoil.m_tf_wp_steel_conduit = (
    #  fcondsteel * a_tf_wp_no_insulation *data.tfcoil.len_tf_coil* fwbs_variables.denstl
    # )


def conduit_insulation(data: DataStructure):
    """Conduit insulation mass [kg]
    (data.tfcoil.a_tf_coil_wp_turn_insulation
    already contains data.tfcoil.n_tf_coil_turns)

    """
    data.tfcoil.m_tf_coil_wp_turn_insulation = (
        data.tfcoil.len_tf_coil
        * data.tfcoil.a_tf_coil_wp_turn_insulation
        * data.tfcoil.den_tf_wp_turn_insulation
    )


def total_conductor(data: DataStructure):
    """[kg] Total conductor mass"""
    data.tfcoil.m_tf_coil_conductor = (
        data.tfcoil.m_tf_coil_superconductor
        + data.tfcoil.m_tf_coil_copper
        + data.tfcoil.m_tf_wp_steel_conduit
        + data.tfcoil.m_tf_coil_wp_turn_insulation
    )


def total_coil(data: DataStructure):
    """[kg] Total coil mass"""
    data.tfcoil.m_tf_coils_total = (
        data.tfcoil.m_tf_coil_case
        + data.tfcoil.m_tf_coil_conductor
        + data.tfcoil.m_tf_coil_wp_insulation
    ) * data.tfcoil.n_tf_coils
