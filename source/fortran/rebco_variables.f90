module rebco_variables
  !! author: M. Kovari
  !!
  !! Module for the REBCO HTS superconductor variables
  !!
  !! Variables relating to the REBCO HTS tape, strand and conductor
  !! Conduit information is in the modules relating to each coil.
  !!
  !!### References
  !!
  !! - Updated 13/11/18 using data from Lewandowska et al 2018.

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  real(dp) :: rebco_thickness
  !! thickness of REBCO layer in tape (m) (`iteration variable 138`)

  real(dp) :: copper_thick
  !! thickness of copper layer in tape (m) (`iteration variable 139`)

  real(dp) :: hastelloy_thickness
  !! thickness of Hastelloy layer in tape (m)

  real(dp) :: tape_width
  !! Mean width of tape (m)

  real(dp) :: tape_thickness
  !! thickness of tape, inc. all layers (hts, copper, substrate, etc.) (m)

  real(dp) :: croco_od
  !! Outer diameter of CroCo strand (m)

  real(dp) :: croco_id
  !! Inner diameter of CroCo copper tube (m)

  real(dp) :: croco_thick
  !! Thickness of CroCo copper tube (m) (`iteration variable 158`)

  real(dp) :: copper_rrr
  !! residual resistivity ratio copper in TF superconducting cable

  real(dp) :: copperA_m2
  !! TF coil current / copper area (A/m2)

  real(dp) :: coppera_m2_max
  !! Maximum TF coil current / copper area (A/m2)

  real(dp) :: f_coppera_m2
  !! f-value for constraint 75: TF coil current / copper area < copperA_m2_max

  real(dp) :: copperaoh_m2
  !! CS coil current / copper area (A/m2) (`sweep variable 61`)

  real(dp) :: copperaoh_m2_max
  !! Maximum CS coil current / copper area (A/m2)

  real(dp) :: f_copperaoh_m2
  !! f-value for constraint 88: CS coil current / copper area < copperA_m2_max

  !#TODO: variables need descriptions and units
  real(dp) :: stack_thickness
  real(dp) :: tapes
  real(dp) :: rebco_area
  real(dp) :: copper_area
  real(dp) :: hastelloy_area
  real(dp) :: solder_area
  real(dp) :: croco_area

  contains

  subroutine init_rebco_variables
    !! Initialise module variables
    implicit none

    rebco_thickness = 1.0D-6
    copper_thick = 100.0D-6
    hastelloy_thickness = 50.0D-6
    tape_width = 4.0D-3
    croco_od = 0.0D0
    croco_id = 0.0D0
    croco_thick = 2.5D-3
    copper_rrr = 100d0
    coppera_m2_max = 1D8
    f_coppera_m2 = 1d0
    tape_thickness = 6.5D-5
    stack_thickness = 0.0D0
    tapes = 0.0D0
    rebco_area = 0.0D0
    copper_area = 0.0D0
    hastelloy_area = 0.0D0
    solder_area = 0.0D0
    croco_area = 0.0D0
    copperA_m2 = 0.0D0
    copperaoh_m2_max = 1D8
    f_copperaoh_m2 = 1d0
    copperaoh_m2 = 0.0D0
  end subroutine init_rebco_variables
end module rebco_variables
