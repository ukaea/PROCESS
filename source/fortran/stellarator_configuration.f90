module stellarator_configuration
    !! author: J Lion, IPP Greifswald
    !! Module containing defining parameters for a stellarator
    !!
    !! This module contains a set of constants that defines a
    !! stellarator configuration. These parameters are based on external
    !! calculations and are hardcoded right now into this module. There will be
    !! the possibiltiy to set them via an input file in the future.
    !! The list below will be modified in further commits.
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
   use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    character (len = 20) :: stella_config_name
   ! Name of the configuration

   integer stella_config_symmetry
   !  Number of coils [1]

   integer stella_config_coilspermodule
   !  Coils per module [1]

   real(dp)  stella_config_rmajor_ref
   !  Reference Point for major radius where all the other variables are determined [m]

   real(dp)  stella_config_rminor_ref
   !  Reference Point for minor radius where all the other variables are determined [m]

   real(dp)  stella_config_coil_rmajor
   !  Reference Point for coil major radius [m]

   real(dp)  stella_config_coil_rminor
   !  Reference Point for coil minor radius [m]

   real(dp)  stella_config_aspect_ref
   !  Reference Point for aspect ratio where all the other variables are determined [1]

   real(dp)  stella_config_bt_ref
   !  Reference Point for toroidal b where all the other variables are determined [T]

   real(dp) stella_config_WP_area
   !  Winding pack area at the reference point [m^2]

   real(dp)  stella_config_WP_bmax
   !  The maximal magnetic field in the winding pack at the reference size of the winding pack [T]

   real(dp) stella_config_i0
   !  Coil current needed for b0 at the reference point [MA]

   real(dp) stella_config_a1
   !  Magnetic field fit parameter a1 (for the maximal field on the coils) [1]

   real(dp) stella_config_a2
   !  Magnetic field fit parameter a2 [1]

   real(dp) stella_config_dmin
   !  Minimal intercoil distance at the reference point [m]

   real(dp) stella_config_inductance
   !  inductance at the reference point [H]

   real(dp) stella_config_coilsurface
   !  Coil surface at the reference point [m2]

   real(dp) stella_config_coillength
   !  Total coil length at the reference point [m]

   real(dp) stella_config_max_portsize_width
   !  Port size in toroidal direction at the reference point [m]

   real(dp) stella_config_maximal_coil_height
   !  The maximal coil height at reference point. [m]

   real(dp) stella_config_min_plasma_coil_distance
   !  The minimal distance between coil and plasma at the reference point [m]

   real(dp) stella_config_derivative_min_LCFS_coils_dist
   !  The derivative of min_plasma_coil_distance wrt to the minor plasma radius at the reference point [1]

   real(dp) stella_config_vol_plasma
   !  The plasma volume at the reference point. Scales as a*R^2. [m^3]

   real(dp) stella_config_plasma_surface
   !  The plasma surface a the reference point. [m^2]

   real(dp) stella_config_WP_ratio
   !  Ratio radial to toroidal length of the winding pack. (a1 and a2 should be calculated using this value) [1]

   real(dp) stella_config_max_force_density
   !  Maximal toroidal and radially averaged force density at reference point in a WP cross section [MN/m^3]

   real(dp) stella_config_max_force_density_MNm
   !  Maximal integrated force density at reference point in a WP cross section [MN/m]

   real(dp) stella_config_min_bend_radius
   !  Minimal bending radius at reference point [m]

   real(dp) stella_config_epseff
   !  Maximal epsilon effective in the core region [1]

   real(dp) stella_config_max_lateral_force_density
   !  Maximal lateral force density of the coil set [MN/m]

   real(dp) stella_config_max_radial_force_density
   !  Maximal radial force density of the coil set [MN/m]

   real(dp) stella_config_centering_force_max_MN
   !  Maximal centering force of a coil in the coil set [MN]

   real(dp) stella_config_centering_force_min_MN
   !  Minimal centering force of a coil in the coil set (negative means pointing outwards) [MN]

   real(dp) stella_config_centering_force_avg_MN
   !  Average centering force the coils in the coil set [MN/coil]

   real(dp) :: stella_config_neutron_peakfactor
   !  The neutron peaking factor determined through inhomogeneities on the stellarator wall (qmax/qavg) [1]
 end module stellarator_configuration
