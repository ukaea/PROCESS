module fwbs_module

  !! Module containing first wall, blanket and shield routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains routines for calculating the
  !! parameters of the first wall, blanket and shield components
  !! of a fusion power plant.

  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  private
  public :: sctfcoil_nuclear_heating_iter90 !, blanket_neutronics

contains

!   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sctfcoil_nuclear_heating_iter90(coilhtmx,dpacop,htheci,nflutf, &
       pheci,pheco,ptfiwp,ptfowp,raddose,ptfnuc)

    !! Superconducting TF coil nuclear heating estimate
    !! author: P J Knight, CCFE, Culham Science Centre
    !! coilhtmx : output real : peak magnet heating (MW/m3)
    !! dpacop : output real : copper stabiliser displacements/atom
    !! htheci : output real : peak TF coil case heating (MW/m3)
    !! nflutf : output real : maximum neutron fluence (n/m2)
    !! pheci : output real : inboard coil case heating (MW)
    !! pheco : output real : outboard coil case heating (MW)
    !! ptfiwp : output real : inboard TF coil winding pack heating (MW)
    !! ptfowp : output real : outboard TF coil winding pack heating (MW)
    !! raddose : output real : insulator dose (rad)
    !! ptfnuc : output real : TF coil nuclear heating (MW)
    !! This subroutine calculates the nuclear heating in the
    !! superconducting TF coils, assuming an exponential neutron
    !! attenuation through the blanket and shield materials.
    !! The estimates are based on 1990 ITER data.
    !! <P>The arrays <CODE>coef(i,j)</CODE> and <CODE>decay(i,j)</CODE>
    !! are used for exponential decay approximations of the
    !! (superconducting) TF coil nuclear parameters.
    !! <UL><P><LI><CODE>j = 1</CODE> : stainless steel shield (assumed)
    !! <P><LI><CODE>j = 2</CODE> : tungsten shield (not used)</UL>
    !! Note: Costing and mass calculations elsewhere assume
    !! stainless steel only.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use build_variables, only: blnkith, blnkoth, fwith, fwoth, shldith, shldoth
    use cost_variables, only: cfactr, tlife
    use physics_variables, only: wallmw
    use tfcoil_variables, only: casthi, i_tf_sup, tfsai, tfsao, dr_tf_wp, &
      tinstf

		use maths_library, only: tril
    implicit none

    !  Arguments

    real(dp), intent(out) :: coilhtmx,dpacop,htheci,nflutf, &
         pheci,pheco,ptfiwp,ptfowp,raddose,ptfnuc

    !  Local variables

    integer, parameter :: ishmat = 1  !  stainless steel coil casing is assumed

    real(dp), dimension(5) :: fact
    real(dp), dimension(5,2) :: coef
    real(dp), dimension(7,2) :: decay
    real(dp) :: dshieq,dshoeq,fpsdt,fpydt,ptfi,ptfo,wpthk

    !  Global shared variables

    !  Input: blnkith,blnkoth,casthi,cfactr,fwith,fwoth,i_tf_sup,shldith
    !  Input: shldoth,tfsai,tfsao,dr_tf_wp,tinstf,tlife,wallmw

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (i_tf_sup /= 1) then  !  Resistive coils
       coilhtmx = 0.0D0
       ptfiwp = 0.0D0
       ptfowp = 0.0D0
       htheci = 0.0D0
       pheci = 0.0D0
       pheco = 0.0D0
       raddose = 0.0D0
       nflutf = 0.0D0
       dpacop = 0.0D0

       ptfnuc = 0.0D0

    else

       !  TF coil nuclear heating coefficients in region i (first element),
       !  assuming shield material j (second element where present)

       fact(1) = 8.0D0
       fact(2) = 8.0D0
       fact(3) = 6.0D0
       fact(4) = 4.0D0
       fact(5) = 4.0D0

       coef(1,1) = 10.3D0
       coef(2,1) = 11.6D0
       coef(3,1) = 7.08D5
       coef(4,1) = 2.19D18
       coef(5,1) = 3.33D-7
       coef(1,2) = 8.32D0
       coef(2,2) = 10.6D0
       coef(3,2) = 7.16D5
       coef(4,2) = 2.39D18
       coef(5,2) = 3.84D-7

       decay(1,1) = 10.05D0
       decay(2,1) = 17.61D0
       decay(3,1) = 13.82D0
       decay(4,1) = 13.24D0
       decay(5,1) = 14.31D0
       decay(6,1) = 13.26D0
       decay(7,1) = 13.25D0
       decay(1,2) = 10.02D0
       decay(2,2) = 3.33D0
       decay(3,2) = 15.45D0
       decay(4,2) = 14.47D0
       decay(5,2) = 15.87D0
       decay(6,2) = 15.25D0
       decay(7,2) = 17.25D0

       !  N.B. The vacuum vessel appears to be ignored

       dshieq = shldith + fwith + blnkith
       dshoeq = shldoth + fwoth + blnkoth

       !  Winding pack radial thickness, including groundwall insulation

       wpthk = dr_tf_wp + 2.0D0*tinstf

       !  Nuclear heating rate in inboard TF coil (MW/m**3)

       coilhtmx = fact(1) * wallmw * coef(1,ishmat) * &
            exp(-decay(6,ishmat) * (dshieq + casthi))

       !  Total nuclear heating (MW)

       ptfiwp = coilhtmx * tfsai * &
            (1.0D0-exp(-decay(1,ishmat)*wpthk)) / decay(1,ishmat)
       ptfowp = fact(1) * wallmw * coef(1,ishmat) * &
            exp(-decay(6,ishmat) * (dshoeq + casthi)) * tfsao * &
            (1.0D0 - exp(-decay(1,ishmat)*wpthk)) / decay(1,ishmat)

       !  Nuclear heating in plasma-side TF coil case (MW)

       htheci = fact(2) * wallmw * coef(2,ishmat) * &
            exp(-decay(7,ishmat) * dshieq)
       pheci = htheci * tfsai * (1.0D0-exp(-decay(2,ishmat)*casthi))/ &
            decay(2,ishmat)
       pheco = fact(2) * wallmw * coef(2,ishmat) * &
            exp(-decay(7,ishmat) * dshoeq) * tfsao * &
            (1.0D0-exp(-decay(2,ishmat)*casthi))/decay(2,ishmat)
       ptfi = ptfiwp + pheci
       ptfo = ptfowp + pheco

       ptfnuc = ptfi + ptfo

       !  Full power DT operation years for replacement of TF Coil
       !  (or plant life)

       fpydt = cfactr * tlife
       fpsdt = fpydt * 3.154D7  !  seconds

       !  Insulator dose (rad)

       raddose = coef(3,ishmat) * fpsdt * fact(3) * wallmw * &
            exp(-decay(3,ishmat) * (dshieq+casthi))

       !  Maximum neutron fluence in superconductor (n/m**2)

       nflutf = fpsdt * fact(4) * wallmw * coef(4,ishmat) * &
            exp(-decay(4,ishmat) * (dshieq+casthi))

       !  Atomic displacement in copper stabilizer

       dpacop = fpsdt * fact(5) * wallmw * coef(5,ishmat) * &
            exp(-decay(5,ishmat) * (dshieq + casthi) )

    end if

  end subroutine sctfcoil_nuclear_heating_iter90


end module fwbs_module
