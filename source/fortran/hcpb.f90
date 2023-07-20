
module ccfe_hcpb_module
  !! author: J Morris (UKAEA)
  !!
  !! This module contains the PROCESS CCFE HCPB blanket model
  !! based on CCFE HCPB model from the PROCESS engineering paper
  !! PROCESS Engineering paper (M. Kovari et al.)
  !!
  !!### References
  !!
  !! - Kovari et al., Fusion Engineering and Design 104 (2016) 9-20

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  ! Variables for output to file
  integer :: ip, ofile

  ! Smeared densities of build sections
  real(dp) :: armour_density
  !! FW armour density [kg/m3]

  real(dp) :: fw_density
  !! FW density [kg/m3]

  real(dp) :: blanket_density
  !! Blanket density [kg/m3]

  real(dp) :: shield_density
  !! Shield density [kg/m3]

  real(dp) :: vv_density
  !! Vacuum vessel density [kg/m3]

  real(dp) :: x_blanket
  !! Blanket exponent (tonne/m2)

  real(dp) :: x_shield
  !! Shield exponent (tonne/m2)

  real(dp) :: tfc_nuc_heating
  !! Unit nuclear heating in TF coil (W per W of fusion power)

  real(dp) :: fw_armour_u_nuc_heating
  !! Unit heating of FW and armour in FW armour (W/kg per W of fusion power)

  real(dp) :: shld_u_nuc_heating
  !! Unit nuclear heating in shield (W per W of fusion power)

  real(dp) :: pnuc_tot_blk_sector
  !! Total nuclear power deposited in blanket covered sector (FW, BLKT, SHLD, TF) (MW)

  real(dp) :: exp_blanket, exp_shield1, exp_shield2
  !! Exponential factors in nuclear heating calcs

contains

  subroutine init_ccfe_hcpb_module
    !! Initialise module variables
    implicit none

    ip = 0
    ofile = 0
    armour_density = 0.0D0
    fw_density = 0.0D0
    blanket_density = 0.0D0
    shield_density = 0.0D0
    vv_density = 0.0D0
    x_blanket = 0.0D0
    x_shield = 0.0D0
    tfc_nuc_heating = 0.0D0
    fw_armour_u_nuc_heating = 0.0D0
    shld_u_nuc_heating = 0.0D0
    exp_blanket = 0.0D0
    exp_shield1 = 0.0D0
    exp_shield2 = 0.0D0
  end subroutine init_ccfe_hcpb_module

end module ccfe_hcpb_module
