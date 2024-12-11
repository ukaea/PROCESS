module maths_library

  !! Library of mathematical and numerical routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains a large number of routines to enable
  !! PROCESS to perform a variety of numerical procedures, including
  !! linear algebra, zero finding, integration and minimisation.
  !! The module is an amalgamation of the contents of several
  !! different pre-existing PROCESS source files, which themselves
  !! were derived from a number of different numerical libraries
  !! including BLAS and MINPAC.
  !! http://en.wikipedia.org/wiki/Gamma_function
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  !use process_output

  implicit none

contains
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eshellvol(rshell,rmini,rmino,zminor,drin,drout,dz,vin,vout,vtot)

    !! Routine to calculate the inboard, outboard and total volumes
    !! of a toroidal shell comprising two elliptical sections
    !! author: P J Knight, CCFE, Culham Science Centre
    !! rshell : input real : major radius of centre of both ellipses (m)
    !! rmini  : input real : horizontal distance from rshell to outer edge
    !! of inboard elliptical shell (m)
    !! rmino  : input real : horizontal distance from rshell to inner edge
    !! of outboard elliptical shell (m)
    !! zminor : input real : vertical internal half-height of shell (m)
    !! drin   : input real : horiz. thickness of inboard shell at midplane (m)
    !! drout  : input real : horiz. thickness of outboard shell at midplane (m)
    !! dz     : input real : vertical thickness of shell at top/bottom (m)
    !! vin    : output real : volume of inboard section (m3)
    !! vout   : output real : volume of outboard section (m3)
    !! vtot   : output real : total volume of shell (m3)
    !! This routine calculates the volume of the inboard and outboard sections
    !! of a toroidal shell defined by two co-centred semi-ellipses.
    !! Each section's internal and external surfaces are in turn defined
    !! by two semi-ellipses. The volumes of each section are calculated as
    !! the difference in those of the volumes of revolution enclosed by their
    !! inner and outer surfaces.
    !! <P>See also <A HREF="eshellarea.html"><CODE>eshellarea</CODE></A>
    !! Internal CCFE note T&amp;M/PKNIGHT/PROCESS/009, P J Knight:
    !! Surface Area and Volume Calculations for Toroidal Shells
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use constants, only: pi, twopi
    implicit none

    !  Arguments
    real(kind=dp), intent(in) :: rshell, rmini, rmino, zminor, drin, drout, dz
    real(kind=dp), intent(out) :: vin, vout, vtot

    !  Local variables
    real(kind=dp) :: a, b, elong, v1, v2

    !  Global shared variables
    !  Input: pi,twopi
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! #TODO - Review both equations containing dz and attempt to separate
    !         top and bottom of vacuum vessel thickness
    !         See issue #433 for explanation

    !  Inboard section

    !  Volume enclosed by outer (higher R) surface of elliptical section
    !  and the vertical straight line joining its ends
    a = rmini ; b = zminor ; elong = b/a
    v1 = twopi * elong * (0.5D0*pi*rshell*a*a - 2.0D0/3.0D0*a*a*a)

    !  Volume enclosed by inner (lower R) surface of elliptical section
    !  and the vertical straight line joining its ends
    a = rmini+drin ; b = zminor+dz ; elong = b/a
    v2 = twopi * elong * (0.5D0*pi*rshell*a*a - 2.0D0/3.0D0*a*a*a)

    !  Volume of inboard section of shell
    vin = v2 - v1

    !  Outboard section

    !  Volume enclosed by inner (lower R) surface of elliptical section
    !  and the vertical straight line joining its ends
    a = rmino ; b = zminor ; elong = b/a
    v1 = twopi * elong * (0.5D0*pi*rshell*a*a + 2.0D0/3.0D0*a*a*a)

    !  Volume enclosed by outer (higher R) surface of elliptical section
    !  and the vertical straight line joining its ends
    a = rmino+drout ; b = zminor+dz ; elong = b/a
    v2 = twopi * elong * (0.5D0*pi*rshell*a*a + 2.0D0/3.0D0*a*a*a)

    !  Volume of outboard section of shell
    vout = v2 - v1

    !  Total shell volume
    vtot = vin + vout

  end subroutine eshellvol

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dshellvol(rmajor,rminor,zminor,drin,drout,dz,vin,vout,vtot)

    !! Routine to calculate the inboard, outboard and total volumes
    !! of a D-shaped toroidal shell
    !! author: P J Knight, CCFE, Culham Science Centre
    !! rmajor : input real : major radius to outer point of inboard
    !! straight section of shell (m)
    !! rminor : input real : horizontal internal width of shell (m)
    !! zminor : input real : vertical internal half-height of shell (m)
    !! drin   : input real : horiz. thickness of inboard shell at midplane (m)
    !! drout  : input real : horiz. thickness of outboard shell at midplane (m)
    !! dz     : input real : vertical thickness of shell at top/bottom (m)
    !! vin    : output real : volume of inboard straight section (m3)
    !! vout   : output real : volume of outboard curved section (m3)
    !! vtot   : output real : total volume of shell (m3)
    !! This routine calculates the volume of the inboard and outboard sections
    !! of a D-shaped toroidal shell defined by the above input parameters.
    !! The inboard section is assumed to be a cylinder of uniform thickness.
    !! The outboard section's internal and external surfaces are defined
    !! by two semi-ellipses, centred on the outer edge of the inboard section;
    !! its volume is calculated as the difference in those of the volumes of
    !! revolution enclosed by the two surfaces.
    !! <P>See also <A HREF="dshellarea.html"><CODE>dshellarea</CODE></A>
    !! Internal CCFE note T&amp;M/PKNIGHT/PROCESS/009, P J Knight:
    !! Surface Area and Volume Calculations for Toroidal Shells
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use constants, only: pi, twopi
    implicit none

    !  Arguments
    real(kind=dp), intent(in) :: rmajor, rminor, zminor, drin, drout, dz
    real(kind=dp), intent(out) :: vin, vout, vtot

    !  Local variables
    real(kind=dp) :: a, b, elong, v1, v2

    !  Global shared variables
    !  Input: pi,twopi
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! #TODO - Review both equations containing dz and attempt to separate
    !         top and bottom of vacuum vessel thickness
    !         See issue #433 for explanation

    !  Volume of inboard cylindrical shell
    vin = 2.0D0*(zminor+dz) * pi*(rmajor**2 - (rmajor-drin)**2)

    !  Volume enclosed by inner surface of elliptical outboard section
    !  and the vertical straight line joining its ends
    a = rminor ; b = zminor ; elong = b/a
    v1 = twopi * elong * (0.5D0*pi*rmajor*a*a + 2.0D0/3.0D0*a*a*a)

    !  Volume enclosed by outer surface of elliptical outboard section
    !  and the vertical straight line joining its ends
    a = rminor+drout ; b = zminor+dz ; elong = b/a
    v2 = twopi * elong * (0.5D0*pi*rmajor*a*a + 2.0D0/3.0D0*a*a*a)

    !  Volume of elliptical outboard shell
    vout = v2 - v1

    !  Total shell volume
    vtot = vin + vout

  end subroutine dshellvol

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! ------------------------------------------------------------------------
  pure function variable_error(variable)
      real(dp), intent(in) ::variable
      logical::variable_error

      if((variable/=variable).or.(variable<-9.99D99).or.(variable>9.99D99))then
          variable_error = .TRUE.
      else
          variable_error = .FALSE.
      end if

  end function variable_error

end module maths_library
