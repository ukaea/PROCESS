"""Module for coil force calculations in stellarators."""

from process.core.model import DataStructure


def calculate_max_force_density(a_tf_wp_no_insulation, data: DataStructure):
    """Calculate the maximum force density in the TF coil winding pack from scaling. [MN/m3]

    Parameters
    ----------
    a_tf_wp_no_insulation :

    data: DataStructure
        data structure object
    """
    data.tfcoil.max_force_density = (
        data.stellarator_config.stella_config_max_force_density
        * data.stellarator.f_st_i_total
        / data.stellarator.f_st_n_coils
        * data.tfcoil.b_tf_inboard_peak_symmetric
        / data.stellarator_config.stella_config_wp_bmax
        * data.stellarator_config.stella_config_wp_area
        / a_tf_wp_no_insulation
    )


def calculate_max_force_density_mnm(data: DataStructure):
    """Calculate the maximum force per meter in the TF coil winding pack from scaling. [MN/m]"""
    return (
        data.stellarator_config.stella_config_max_force_density_mnm
        * data.stellarator.f_st_i_total
        / data.stellarator.f_st_n_coils
        * data.tfcoil.b_tf_inboard_peak_symmetric
        / data.stellarator_config.stella_config_wp_bmax
    )


def calculate_maximum_stress(data):
    """Approximate, very simple maxiumum stress (needed for limitation of icc 32), in Pa"""
    data.tfcoil.sig_tf_wp = (
        data.tfcoil.max_force_density * data.tfcoil.dr_tf_wp_with_insulation * 1.0e6
    )


def calculate_max_lateral_force_density(a_tf_wp_no_insulation, data: DataStructure):
    """Calculate the maximum lateral force density in the TF coil winding pack from scaling. [MN/m3]

    Parameters
    ----------
    a_tf_wp_no_insulation :

    data: DataStructure
        data structure object
    """
    return (
        data.stellarator_config.stella_config_max_lateral_force_density
        * data.stellarator.f_st_i_total
        / data.stellarator.f_st_n_coils
        * data.tfcoil.b_tf_inboard_peak_symmetric
        / data.stellarator_config.stella_config_wp_bmax
        * data.stellarator_config.stella_config_wp_area
        / a_tf_wp_no_insulation
    )


def calculate_max_radial_force_density(a_tf_wp_no_insulation, data):
    """Calculate the maximum radial force density in the TF coil winding pack from scaling. [MN/m3]

    Parameters
    ----------
    a_tf_wp_no_insulation :

    data: DataStructure
        data structure object
    """
    return (
        data.stellarator_config.stella_config_max_radial_force_density
        * data.stellarator.f_st_i_total
        / data.stellarator.f_st_n_coils
        * data.tfcoil.b_tf_inboard_peak_symmetric
        / data.stellarator_config.stella_config_wp_bmax
        * data.stellarator_config.stella_config_wp_area
        / a_tf_wp_no_insulation
    )


def calculate_centering_force_max_mn(data: DataStructure):
    """Calculate the maximum centering force in the TF coils from scaling. [MN]"""
    return (
        data.stellarator_config.stella_config_centering_force_max_mn
        * data.stellarator.f_st_i_total
        / data.stellarator.f_st_n_coils
        * data.tfcoil.b_tf_inboard_peak_symmetric
        / data.stellarator_config.stella_config_wp_bmax
        * data.stellarator_config.stella_config_coillength
        / data.tfcoil.n_tf_coils
        / data.tfcoil.len_tf_coil
    )


def calculate_centering_force_min_mn(data: DataStructure):
    """Calculate the minimum centering force in the TF coils from scaling. [MN]"""
    return (
        data.stellarator_config.stella_config_centering_force_min_mn
        * data.stellarator.f_st_i_total
        / data.stellarator.f_st_n_coils
        * data.tfcoil.b_tf_inboard_peak_symmetric
        / data.stellarator_config.stella_config_wp_bmax
        * data.stellarator_config.stella_config_coillength
        / data.tfcoil.n_tf_coils
        / data.tfcoil.len_tf_coil
    )


def calculate_centering_force_avg_mn(data: DataStructure):
    """Calculate the average centering force in the TF coils from scaling. [MN]"""
    return (
        data.stellarator_config.stella_config_centering_force_avg_mn
        * data.stellarator.f_st_i_total
        / data.stellarator.f_st_n_coils
        * data.tfcoil.b_tf_inboard_peak_symmetric
        / data.stellarator_config.stella_config_wp_bmax
        * data.stellarator_config.stella_config_coillength
        / data.tfcoil.n_tf_coils
        / data.tfcoil.len_tf_coil
    )
