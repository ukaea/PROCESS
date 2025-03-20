
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
end module ccfe_hcpb_module
