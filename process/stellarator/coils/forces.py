"""Module for coil force calculations in stellarators."""

from process.fortran import (
    stellarator_configuration,
    tfcoil_variables,
)
from process.fortran import (
    stellarator_module as st,
)

def calculate_max_force_density(a_tf_wp_no_insulation):

    tfcoil_variables.max_force_density = (
        stellarator_configuration.stella_config_max_force_density
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_wp_area
        / a_tf_wp_no_insulation
    )


def calculate_max_force_density_mnm():
    return (
        stellarator_configuration.stella_config_max_force_density_mnm
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
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
    return (
        stellarator_configuration.stella_config_max_lateral_force_density
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_wp_area
        / a_tf_wp_no_insulation
    )


def calculate_max_radial_force_density(a_tf_wp_no_insulation):
    return (
        stellarator_configuration.stella_config_max_radial_force_density
        * st.f_i
        / st.f_n
        * tfcoil_variables.b_tf_inboard_peak
        / stellarator_configuration.stella_config_wp_bmax
        * stellarator_configuration.stella_config_wp_area
        / a_tf_wp_no_insulation
    )


def calculate_centering_force_max_mn():
    return (
    stellarator_configuration.stella_config_centering_force_max_mn
    * st.f_i
    / st.f_n
    * tfcoil_variables.b_tf_inboard_peak
    / stellarator_configuration.stella_config_wp_bmax
    * stellarator_configuration.stella_config_coillength
    / tfcoil_variables.n_tf_coils
    / tfcoil_variables.len_tf_coil
    )


def calculate_centering_force_min_mn():
    return (
    stellarator_configuration.stella_config_centering_force_min_mn
    * st.f_i
    / st.f_n
    * tfcoil_variables.b_tf_inboard_peak
    / stellarator_configuration.stella_config_wp_bmax
    * stellarator_configuration.stella_config_coillength
    / tfcoil_variables.n_tf_coils
    / tfcoil_variables.len_tf_coil
    )


def calculate_centering_force_avg_mn():
    return (
    stellarator_configuration.stella_config_centering_force_avg_mn
    * st.f_i
    / st.f_n
    * tfcoil_variables.b_tf_inboard_peak
    / stellarator_configuration.stella_config_wp_bmax
    * stellarator_configuration.stella_config_coillength
    / tfcoil_variables.n_tf_coils
    / tfcoil_variables.len_tf_coil
    )