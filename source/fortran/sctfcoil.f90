module sctfcoil_module

!! Module containing superconducting TF coil routines
!! author: P J Knight, CCFE, Culham Science Centre
!! author: J Morris, CCFE, Culham Science Centre
!! author: S Kahn, CCFE, Culham Science Centre
!! N/A
!! This module contains routines for calculating the
!! parameters of a superconducting TF coil system for a
!! fusion power plant.
!! PROCESS Superconducting TF Coil Model, J. Morris, CCFE, 1st May 2014
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
   use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

   implicit none

! Module variables
!-----------------

   real(dp) :: tf_fit_t
!! Dimensionless winding pack width

   real(dp) :: tf_fit_z
!! Dimensionless winding pack radial thickness

   real(dp) :: tf_fit_y
!! Ratio of peak field with ripple to nominal axisymmetric peak field

   real(dp) :: c_tf_coil
!! Current in each TF coil

   real(dp) :: awpc
!! Total cross-sectional area of winding pack including
!! GW insulation and insertion gap [m2]

   real(dp) :: awptf
!! Total cross-sectional area of winding pack without
!! ground insulation and insertion gap [m2]

   real(dp) :: a_tf_steel
!! Inboard coil steel coil cross-sectional area [m2]

   real(dp) :: a_tf_ins
!! Inboard coil insulation cross-section per coil [m2]

   real(dp) :: f_tf_steel
!! Inboard coil steel fraction [-]

   real(dp) :: f_tf_ins
!! Inboard coil insulation fraction [-]

   real(dp) :: z_cp_top
!! Vertical distance from the midplane to the top of the tapered section [m]

   real(dp) :: r_tf_outboard_in
!! Radial position of plasma-facing edge of TF coil outboard leg [m]

   real(dp) :: r_tf_outboard_out
!! Radial position of outer edge of TF coil inboard leg [m]

   real(dp) :: r_wp_inner
!! Radial position of inner edge and centre of winding pack [m]

   real(dp) :: r_wp_outer
!! Radial position of outer edge and centre of winding pack [m]

   real(dp) :: r_wp_centre
!! Radial position of centre and centre of winding pack [m]

   real(dp) :: dr_tf_wp_top
!! Conductor layer radial thickness at centercollumn top [m]
!! Ground insulation layer included, only defined for itart = 1

   real(dp) :: vol_ins_cp
!! CP turn insulation volume [m3]

   real(dp) :: vol_gr_ins_cp
!! CP ground insulation volume [m3]

   real(dp) :: vol_case_cp
!! Volume of the CP outer casing cylinder

   real(dp) :: t_wp_toroidal
!! Minimal toroidal thickness of of winding pack [m]

   real(dp) :: t_wp_toroidal_av
!! Averaged toroidal thickness of of winding pack [m]

   real(dp) :: t_lat_case_av
!! Average lateral casing thickness [m]

   real(dp) :: a_case_front
!! Front casing area [m2]

   real(dp) :: a_case_nose
!! Nose casing area [m2]

   real(dp) :: a_ground_ins
!! Inboard mid-plane cross-section area of the WP ground insulation [m2]

   real(dp) :: a_leg_ins
!! TF ouboard leg turn insulation area per coil [m2]

   real(dp) :: a_leg_gr_ins
!! TF outboard leg ground insulation area per coil [m2]

   real(dp) :: a_leg_cond
!! Exact TF ouboard leg conductor area [m2]

   real(dp) :: rad_tf_coil_toroidal
!! Half toroidal angular extent of a single TF coil inboard leg

   real(dp) :: tan_theta_coil
!! Tan half toroidal angular extent of a single TF coil inboard leg

   real(dp) :: t_conductor_radial, t_conductor_toroidal
!! Conductor area radial and toroidal dimension (integer turn only) [m]

   real(dp) :: dr_tf_turn_cable_space, t_cable_toroidal
!! Cable area radial and toroidal dimension (integer turn only) [m]

   real(dp) :: dr_tf_turn, t_turn_toroidal
!! Turn radial and toroidal dimension (integer turn only) [m]

   real(dp) :: t_cable
!! Cable area averaged dimension (square shape) [m]

   real(dp) :: vforce_inboard_tot
!! Total inboard vertical tension (all coils) [N]

! Vacuum Vessel stress on TF coil quench

   real(dp) :: vv_stress_quench
   !! The Tresca stress experienced by the Vacuum Vessel when the SCTF coil quenches [Pa]

! croco_strand
   real(dp) :: croco_strand_area
   real(dp) :: croco_strand_critical_current

! conductor
   real(dp) :: conductor_copper_area,  conductor_copper_fraction
   real(dp) :: conductor_copper_bar_area
   real(dp) :: conductor_hastelloy_area, conductor_hastelloy_fraction
   real(dp) :: conductor_helium_area, conductor_helium_fraction
   real(dp) :: conductor_solder_area, conductor_solder_fraction
   real(dp) :: conductor_jacket_area, conductor_jacket_fraction
   real(dp) :: conductor_rebco_area,  conductor_rebco_fraction
   real(dp) :: conductor_critical_current
   real(dp) :: conductor_acs
   real(dp) :: conductor_area
!! Area of cable space inside jacket

   real(dp):: T1, time2, tau2, e_tf_magnetic_stored_total
! (OBSOLETE, but leave for moment)
! real (kind(1.0D0)) ::croco_quench_factor
! real(dp):: jwdgpro_1, jwdgpro_2,  etamax

! Var in tf_res_heating requiring re-initialisation on each new run
! Not sure what is really doing --> to be checked
   integer :: is_leg_cp_temp_same
end module sctfcoil_module
