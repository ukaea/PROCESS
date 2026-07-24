from dataclasses import dataclass


@dataclass(slots=True)
class StellaratorConfigData:
    stella_config_name: str = ""
    """Name of the configuration"""

    stella_config_symmetry: int = 0
    """Number of coils [1]"""

    stella_config_coilspermodule: int = 0
    """Coils per module [1]"""

    stella_config_rmajor_ref: float = 0.0
    """Reference Point for major radius where all the other variables are determined [m]"""

    stella_config_rminor_ref: float = 0.0
    """Reference Point for minor radius where all the other variables are determined [m]"""

    stella_config_coil_rmajor: float = 0.0
    """Reference Point for coil major radius [m]"""

    stella_config_coil_rminor: float = 0.0
    """Reference Point for coil minor radius [m]"""

    stella_config_aspect_ref: float = 0.0
    """Reference Point for aspect ratio where all the other variables are determined [1]"""

    stella_config_bt_ref: float = 0.0
    """Reference Point for toroidal b where all the other variables are determined [T]"""

    stella_config_wp_area: float = 0.0
    """Winding pack area at the reference point [m^2]"""

    stella_config_wp_bmax: float = 0.0
    """The maximal magnetic field in the winding pack at the reference size of the winding pack [T]"""

    stella_config_i0: float = 0.0
    """Coil current needed for b0 at the reference point [MA]"""

    stella_config_a1: float = 0.0
    """Magnetic field fit parameter a1 (for the maximal field on the coils) [1]"""

    stella_config_a2: float = 0.0
    """Magnetic field fit parameter a2 [1]"""

    stella_config_dmin: float = 0.0
    """Minimal intercoil distance at the reference point [m]"""

    stella_config_inductance: float = 0.0
    """inductance at the reference point [H]"""

    stella_config_coilsurface: float = 0.0
    """Coil surface at the reference point [m2]"""

    stella_config_coillength: float = 0.0
    """Total coil length at the reference point [m]"""

    stella_config_max_portsize_width: float = 0.0
    """Port size in toroidal direction at the reference point [m]"""

    stella_config_maximal_coil_height: float = 0.0
    """The maximal coil height at reference point. [m]"""

    stella_config_min_plasma_coil_distance: float = 0.0
    """The minimal distance between coil and plasma at the reference point [m]"""

    stella_config_derivative_min_lcfs_coils_dist: float = 0.0
    """The derivative of min_plasma_coil_distance wrt to the minor plasma radius at the reference point [1]"""

    stella_config_vol_plasma: float = 0.0
    """The plasma volume at the reference point. Scales as a*R^2. [m^3]"""

    stella_config_plasma_surface: float = 0.0
    """The plasma surface a the reference point. [m^2]"""

    stella_config_wp_ratio: float = 0.0
    """Ratio radial to toroidal length of the winding pack. (a1 and a2 should be calculated using this value) [1]"""

    stella_config_max_force_density: float = 0.0
    """Maximal toroidal and radially averaged force density at reference point in a WP cross section [MN/m^3]"""

    stella_config_max_force_density_mnm: float = 0.0
    """Maximal integrated force density at reference point in a WP cross section [MN/m]"""

    stella_config_min_bend_radius: float = 0.0
    """Minimal bending radius at reference point [m]"""

    stella_config_epseff: float = 0.0
    """Maximal epsilon effective in the core region [1]"""

    stella_config_max_lateral_force_density: float = 0.0
    """Maximal lateral force density of the coil set [MN/m]"""

    stella_config_max_radial_force_density: float = 0.0
    """Maximal radial force density of the coil set [MN/m]"""

    stella_config_centering_force_max_mn: float = 0.0
    """Maximal centering force of a coil in the coil set [MN]"""

    stella_config_centering_force_min_mn: float = 0.0
    """Minimal centering force of a coil in the coil set (negative means pointing outwards) [MN]"""

    stella_config_centering_force_avg_mn: float = 0.0
    """Average centering force the coils in the coil set [MN/coil]"""

    stella_config_neutron_peakfactor: float = 0.0
    """The neutron peaking factor determined through inhomogeneities on the stellarator wall (qmax/qavg) [1]"""


CREATE_DICTS_FROM_DATACLASS = StellaratorConfigData
