"""Module for coil force calculations in stellarators."""

from process.data_structure import (
    stellarator_configuration,
    tfcoil_variables,
    stellarator_variables
)

def calculate_max_force_density(a_tf_wp_no_insulation):
    """Calculate the maximum force density in the TF coil winding pack from scaling. [MN/m3]"""

    tfcoil_variables.max_force_density = (
        stellarator_configuration.stella_config_max_force_density
        * stellarator_variables.f_st_i_total
        / stellarator_variables.f_st_n_coils
        * tfcoil_variables.b_tf_inboard_peak_symmetric
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_wp_area
        / a_tf_wp_no_insulation
    )


def calculate_max_force_density_mnm():
    """Calculate the maximum force per meter in the TF coil winding pack from scaling. [MN/m]"""
    return (
        stellarator_configuration.stella_config_max_force_density_mnm
        * stellarator_variables.f_st_i_total
        / stellarator_variables.f_st_n_coils
        * tfcoil_variables.b_tf_inboard_peak_symmetric
        / stellarator_configuration.stella_config_wp_bmax
    )


def calculate_maximum_stress():
    """Approximate, very simple maxiumum stress (needed for limitation of icc 32), in Pa"""
    tfcoil_variables.sig_tf_wp = (
        tfcoil_variables.max_force_density
        * tfcoil_variables.dr_tf_wp_with_insulation
        * 1.0e6
    )


def calculate_max_lateral_force_density(a_tf_wp_no_insulation):
    """Calculate the maximum lateral force density in the TF coil winding pack from scaling. [MN/m3]"""
    return (
        stellarator_configuration.stella_config_max_lateral_force_density
        * stellarator_variables.f_st_i_total
        / stellarator_variables.f_st_n_coils
        * tfcoil_variables.b_tf_inboard_peak_symmetric
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_wp_area
        / a_tf_wp_no_insulation
    )


def calculate_max_radial_force_density(a_tf_wp_no_insulation):
    """Calculate the maximum radial force density in the TF coil winding pack from scaling. [MN/m3]"""
    return (
        stellarator_configuration.stella_config_max_radial_force_density
        * stellarator_variables.f_st_i_total
        / stellarator_variables.f_st_n_coils
        * tfcoil_variables.b_tf_inboard_peak_symmetric
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_wp_area
        / a_tf_wp_no_insulation
    )


def calculate_centering_force_max_mn():
    """Calculate the maximum centering force in the TF coils from scaling. [MN]"""
    return (
    stellarator_configuration.stella_config_centering_force_max_mn
    * stellarator_variables.f_st_i_total
    / stellarator_variables.f_st_n_coils
    * tfcoil_variables.b_tf_inboard_peak_symmetric
    / stellarator_configuration.stella_config_wp_bmax
    * stellarator_configuration.stella_config_coillength
    / tfcoil_variables.n_tf_coils
    / tfcoil_variables.len_tf_coil
    )


def calculate_centering_force_min_mn():
    """Calculate the minimum centering force in the TF coils from scaling. [MN]"""
    return (
    stellarator_configuration.stella_config_centering_force_min_mn
    * stellarator_variables.f_st_i_total
    / stellarator_variables.f_st_n_coils
    * tfcoil_variables.b_tf_inboard_peak_symmetric
    / stellarator_configuration.stella_config_wp_bmax
    * stellarator_configuration.stella_config_coillength
    / tfcoil_variables.n_tf_coils
    / tfcoil_variables.len_tf_coil
    )


def calculate_centering_force_avg_mn():
    """Calculate the average centering force in the TF coils from scaling. [MN]"""
    return (
    stellarator_configuration.stella_config_centering_force_avg_mn
    * stellarator_variables.f_st_i_total
    / stellarator_variables.f_st_n_coils
    * tfcoil_variables.b_tf_inboard_peak_symmetric
    / stellarator_configuration.stella_config_wp_bmax
    * stellarator_configuration.stella_config_coillength
    / tfcoil_variables.n_tf_coils
    / tfcoil_variables.len_tf_coil
    )