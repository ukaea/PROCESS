module superconductors
  !! Module containing superconducter critical surfaces and conductor data
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none
contains

!--------------------------------------------------------------------

subroutine croco(jcritsc, croco_strand_area, croco_strand_critical_current, &
    conductor_copper_area, conductor_copper_fraction, conductor_copper_bar_area, &
    conductor_hastelloy_area, conductor_hastelloy_fraction, conductor_helium_area, &
    conductor_helium_fraction, conductor_solder_area, conductor_solder_fraction, &
    conductor_rebco_area, conductor_rebco_fraction, conductor_critical_current, &
    conductor_area, croco_od,croco_thick)

    !! "CroCo" (cross-conductor) strand and cable design for
    !! "REBCO" 2nd generation HTS superconductor
    ! Updated 13/11/18 using data from Lewandowska et al 2018.

    use rebco_variables, only: copper_area, copper_thick, croco_id, &
      hastelloy_area, hastelloy_thickness, rebco_area, solder_area, &
      stack_thickness, tape_thickness, tape_width, tapes, rebco_thickness
    use resistive_materials, only: volume_fractions, supercon_strand
    use constants, only: pi
    implicit none
    real(dp), intent(in) ::jcritsc
    real(dp) :: d, scaling, croco_od, croco_thick

    ! conductor
    real(dp), intent(inout) :: conductor_copper_area,  conductor_copper_fraction
    real(dp), intent(inout) :: conductor_copper_bar_area
    real(dp), intent(inout) :: conductor_hastelloy_area, conductor_hastelloy_fraction
    real(dp), intent(inout) :: conductor_helium_area, conductor_helium_fraction
    real(dp), intent(inout) :: conductor_solder_area, conductor_solder_fraction
    real(dp), intent(inout) :: conductor_rebco_area,  conductor_rebco_fraction
    real(dp), intent(inout) :: conductor_critical_current
    real(dp), intent(in) :: conductor_area

    ! croco_strand
    real(dp), intent(inout) :: croco_strand_area
    real(dp), intent(inout) :: croco_strand_critical_current


    ! Define local alias
    d = croco_od
    !d = conductor_width / 3.0d0 - thwcndut * ( 2.0d0 / 3.0d0 )

    croco_id = d - 2.0d0 * croco_thick !scaling * 5.4d-3
    if (croco_id <= 0.0d0) then
        write(*,*) 'Warning: negitive inner croco diameter!'
        write(*,*)'croco_id =', croco_id, ',croco_thick = ', croco_thick, ', croco_od =', croco_od
    end if
    ! Define the scaling factor for the input REBCO variable
    ! Ratio of new croco inner diameter and fixed base line value
    scaling = croco_id / 5.4d-3
    tape_width = scaling * 3.75d-3
    ! Properties of a single strand
    tape_thickness = rebco_thickness + copper_thick + hastelloy_thickness
    stack_thickness = sqrt(croco_id**2 - tape_width**2)
    tapes = stack_thickness / tape_thickness

    copper_area = pi * croco_thick * d - pi * croco_thick**2 &  ! copper tube
                  + copper_thick*tape_width*tapes          ! copper in tape
    hastelloy_area = hastelloy_thickness * tape_width * tapes
    solder_area = pi / 4.0d0 * croco_id**2 - stack_thickness * tape_width

    rebco_area = rebco_thickness * tape_width * tapes
    croco_strand_area =  pi / 4.0d0 * d**2
    croco_strand_critical_current = jcritsc * rebco_area

    ! Conductor properties
    !conductor%number_croco = conductor%acs*(1d0-cable_helium_fraction-copper_bar)/croco_strand_area
    conductor_critical_current = croco_strand_critical_current * 6.0d0
    ! Area of core = area of strand
    conductor_copper_bar_area = croco_strand_area
    conductor_copper_area = copper_area * 6.0d0 + conductor_copper_bar_area
    conductor_copper_fraction = conductor_copper_area / conductor_area

    ! Helium area is set by the user.
    !conductor_helium_area = cable_helium_fraction * conductor_acs
    conductor_helium_area = pi / 2.0d0 * d**2
    conductor_helium_fraction = conductor_helium_area / conductor_area

    conductor_hastelloy_area = hastelloy_area * 6.0d0
    conductor_hastelloy_fraction = conductor_hastelloy_area / conductor_area

    conductor_solder_area = solder_area * 6.0d0
    conductor_solder_fraction = conductor_solder_area / conductor_area

    conductor_rebco_area = rebco_area * 6.0d0
    conductor_rebco_fraction = conductor_rebco_area / conductor_area

end subroutine croco

end module superconductors
