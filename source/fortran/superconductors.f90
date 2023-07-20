module superconductors
  !! Module containing superconducter critical surfaces and conductor data
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none
contains

! -------------------------------------------------------------------

subroutine jcrit_rebco(temperature, b, jcrit, validity, iprint)

    !! Critical current density for "REBCO" 2nd generation HTS superconductor
    !! temperature : input real : superconductor temperature (K)
    !! b : input real : Magnetic field at superconductor (T)
    !! jcrit : output real : Critical current density in superconductor (A/m2)
    ! Will return a negative number if the temperature is greater than Tc0, the
    ! zero-field critical temperature.
    use constants, only: pi
    implicit none

    !  Arguments
    real(dp), intent(in) :: temperature, b
    real(dp), intent(out) :: jcrit
    logical, intent(out) :: validity
    integer, intent(in) :: iprint

    !  Local variables
    real(dp) :: birr, factor, tcb

    !  Parameters
    real(dp), parameter :: tc0 = 90.0d0        !  (K)
    real(dp), parameter :: birr0 = 132.5d0     !  (T)
    real(dp), parameter :: a = 1.82962d8       ! scaling constant
    real(dp), parameter :: p = 0.5875d0        ! exponent
    real(dp), parameter :: q = 1.7d0           ! exponent
    real(dp), parameter :: alpha =1.54121d0    ! exponent
    real(dp), parameter :: beta = 1.96679d0    ! exponent
    real(dp), parameter :: oneoveralpha = 1d0 / alpha  ! inverse

    validity = .true.
    if((temperature<4.2d0).or.(temperature>72.0d0) )validity = .false.
    if(temperature<65)then
        if((b<0.0d0).or.(b>15.0d0))validity = .false.
    else
        if((b<0.0d0).or.(b>11.5d0))validity = .false.
    endif

    if((iprint==1).and.(.not.validity))then
        write(*,*)'ERROR in subroutine jcrit_rebco: input out of range'
        write(*,*)'temperature = ', temperature, '   Field = ', b
    endif

    if(temperature<tc0)then
        ! Normal case
        birr = birr0 * (1 - temperature/tc0)**alpha
    else
        ! If temp is greater than critical temp, ensure result is real but negative.
        birr = birr0 * (1 - temperature/tc0)
    end if

    if(b<birr)then
        ! Normal case
        factor = (b/birr)**p * (1-b/birr)**q
        jcrit = (a/b) * (birr**beta) * factor
    else
        ! Field is too high
        ! Ensure result is real but negative, and varies with temperature.
        ! tcb = critical temperature at field b
        tcb = tc0 * (1-(b/birr0)**oneoveralpha)
        jcrit = -(temperature - tcb)
    end if

end subroutine jcrit_rebco
! -------------------------------------------------------------------

subroutine current_sharing_rebco(current_sharing_t, bfield, j)

    !! Current sharing temperature for "REBCO" 2nd generation HTS superconductor
    !! temperature : input real : superconductor temperature (K)
    !! b : input real : Magnetic field at superconductor (T)
    !! j : input real : Current density in superconductor (A/m2)
    !! current_sharing_t : output real : Current sharing temperature (K)
    ! Uses 'function jcrit_rebco' and the finds the temperature for given field and current density
    ! Will return a negative number if the current density is too high

    use maths_library, only: secant_solve
    implicit none

    !  Arguments
    real(dp), intent(in) :: j, bfield
    real(dp), intent(out) :: current_sharing_t
    real(dp)::x1,x2         ! Initial guesses for temperature
    logical::error                   ! True if the solver does not converge
    real(dp)::residual      ! Residual current density error
    real(dp), parameter ::opt_tol = 1d7 ! Tolerance in current density

    x1 = 4d0
    x2 = 20d0
    ! Solve for deltaj_rebco = 0
    call secant_solve(deltaj_rebco,x1,x2,current_sharing_t,error,residual,opt_tol)

contains

    function deltaj_rebco(temperature)
        ! Required because secant_solve requires a function not a subroutine
        ! All we do is throw away the second output, 'validity'!
        ! This needs to be 'contained' because 'bfield' must be available but cannot
        ! be passed, because secant_solve requires a function of one variable!
        ! Also j must be available.
        real(dp), intent(in) :: temperature
        real(dp)::deltaj_rebco, jcritical
        logical :: validity
        integer :: iprint = 0

        call jcrit_rebco(temperature, bfield, jcritical, validity, iprint)
        deltaj_rebco = jcritical - j
    end function deltaj_rebco

end subroutine current_sharing_rebco
!--------------------------------------------------------------------

function function_jcrit_rebco(temperature,b)
    ! Required because secant_solve requires a function not a subroutine
    ! All we do is throw away the second output, 'validity'!
    real(dp), intent(in) :: temperature, b
    real(dp)::function_jcrit_rebco
    logical :: validity
    integer :: iprint = 0

    call jcrit_rebco(temperature, b, function_jcrit_rebco, validity, iprint)
end function function_jcrit_rebco

!--------------------------------------------------------------------

subroutine test_quench()
    use resistive_materials, only: resistive_material
    implicit none

    real(dp) :: jcrit
    real(dp) :: B, delta_b,jcrit42,jcrit14,jcrit22, jcrit33
    real(dp) :: jcrit50, jcrit65, jcrit72, jcrit90
    real(dp) :: T, delta_t, copper0, copper10, hastelloyB, RRR
    real(dp) :: copper0B, copper10B

    integer::i
    logical :: validity
    type(resistive_material)::copper, jacket, solder, hastelloy
    type(resistive_material):: helium
    copper%rrr = 100.0d0
10  format(1x,a,1pe10.2,a,1pe10.2)
20  format(f10.1, 9(1pe10.2))
30  format(9(a10))

    write(*,*)'TEST of subroutine jcrit_rebco:'
    call jcrit_rebco(14.0d0, 15.0d0, jcrit, validity, 0)
    write(*,*)'Temperature = 14 K,   Field = 15 T, validity = ', validity
    write(*,10)'jcrit                                         = ', jcrit
    write(*,*) 'Critical current from Superpower plot (approx)=  2.5e10 A/m2'
    write(*,*)


    open(unit=1,file='REBCO_JC.DAT')
    B = 0.5d0
    delta_b = 0.5d0
    write(1,30)'Field (T)', '4.2', '14', '22', '33', '50', '65', '72', '90'
    do i=1,30
        call jcrit_rebco(4.2d0, B, jcrit42, validity, 0)
        call jcrit_rebco(14.0d0, B, jcrit14, validity, 0)
        call jcrit_rebco(22.0d0, B, jcrit22, validity, 0)
        call jcrit_rebco(33.0d0, B, jcrit33, validity, 0)
        call jcrit_rebco(50.0d0, B, jcrit50, validity, 0)
        call jcrit_rebco(65.0d0, B, jcrit65, validity, 0)
        call jcrit_rebco(72.0d0, B, jcrit72, validity, 0)
        ! Beyond the critical temperature at any non-zero field
        call jcrit_rebco(90.0d0, B, jcrit90, validity, 0)

        write(1,20) B, jcrit42,jcrit14,jcrit22, jcrit33, jcrit50, jcrit65, jcrit72, jcrit90
        B = B + delta_b
    end do
    close(unit=1)
    write(*,*) 'Critical current data written to REBCO_JC.DAT'

    open(unit=1,file='quench_data.DAT')

    T = 4d0
    delta_t = 1d0
    write(1,*)'Solids: Cp as a function of Temperature (J/kg/K), and Helium: Cp x density J/K/m3'
    write(1,30)'', 'solder', 'solder(B)','jacket', 'hastell', 'hastell(B)','copper', 'copper(B)', 'He Cp*rho'
    do i=1,297
        call solder_properties(T, solder%cp)
        call jacket_properties(T, jacket%cp)
        call hastelloy_properties(T, hastelloy%cp)
        call copper_properties(T,0.0d0, copper%resistivity, copper%rrr, copper%cp)
        call helium_properties(T, helium%cp_density)
        hastelloyB = cHastelloyC276(T)

        write(1,20) T, solder%cp, cSn40Pb(sngl(T)), jacket%cp, hastelloy%cp, &
                       hastelloyB,copper%cp, cCu(sngl(T)), helium%cp_density
        T = T + delta_t
    enddo
    write(1,*)

    write(1,*)'Resistivity of copper at different fields.  RRR = ', copper%rrr
    T = 4d0
    delta_t = 1d0
    RRR = copper%rrr
    write(1,30)'Temp (K)', '0 T', '10 T', '0 T (B)', '10 T (B)'
    do i=1,297
        call copper_properties(T,0.0d0, copper%resistivity, copper%rrr, copper%cp)
        copper0 = copper%resistivity
        call copper_properties(T,10.0d0, copper%resistivity, copper%rrr, copper%cp)
        copper10 = copper%resistivity
        copper0B = rCu(sngl(T),0.0,sngl(RRR))
        copper10B = rCu(sngl(T),10.0,sngl(RRR))
        write(1,20) T, copper0, copper10,copper0B, copper10B
        T = T + delta_t
    end do
    close(unit=1)
    write(*,*) 'Data written to quench_data.DAT'

end subroutine test_quench
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine itersc(thelium,bmax,strain,bc20max,tc0max,jcrit,bcrit,tcrit)

  !! Implementation of ITER Nb3Sn critical surface implementation
  !! author: R Kemp, CCFE, Culham Science Centre
  !! author: P J Knight, CCFE, Culham Science Centre
  !! thelium : input real : Coolant/SC temperature (K)
  !! bmax : input real : Magnetic field at conductor (T)
  !! strain : input real : Strain in superconductor
  !! bc20max : input real : Upper critical field (T) for superconductor
  !! at zero temperature and strain
  !! tc0max : input real : Critical temperature (K) at zero field and strain
  !! jcrit : output real : Critical current density in superconductor (A/m2)
  !! bcrit : output real : Critical field (T)
  !! tcrit : output real : Critical temperature (K)
  !! This routine calculates the critical current density and
  !! temperature in the superconducting TF coils using the
  !! ITER Nb3Sn critical surface model.
  !! $J_C(B,T,\epsilon)$ Parameterization for ITER Nb3Sn production,
  !! L. Bottura, CERN-ITER Collaboration Report, Version 2, April 2nd 2008
  !! (distributed by Arnaud Devred, ITER, 10th April 2008)
  !! ITER Nb3Sn critical surface parameterization (2MMF7J) (2008),
  !! https://user.iter.org/?uid=2MMF7J&action=get_document
  !! ITER DDD 11-7: Magnets - conductors (2NBKXY) (2009),
  !! https://user.iter.org/?uid=2NBKXY&action=get_document
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use error_handling, only: fdiags, report_error
  use constants, only: pi
  implicit none

  !  Arguments

  real(dp), intent(in) :: thelium, bmax, strain, bc20max, tc0max
  real(dp), intent(out) :: jcrit, bcrit, tcrit

  !  Local variables

  !  Parameters named in Bottura

  !real(dp), parameter :: csc = 16500.0D6 !  scaling constant C
  real(dp), parameter :: csc = 19922.0D6 !  scaling constant C
  real(dp), parameter :: p = 0.63D0      !  low field exponent p
  real(dp), parameter :: q = 2.1D0       !  high field exponent q
  !real(dp), parameter :: ca1 = 44.0D0    !  strain fitting constant C_{a1}
  real(dp), parameter :: ca1 = 44.48D0    !  strain fitting constant C_{a1}
  !real(dp), parameter :: ca2 = 4.0D0     !  strain fitting constant C_{a2}
  real(dp), parameter :: ca2 = 0.0D0     !  strain fitting constant C_{a2}
  real(dp), parameter :: eps0a = 0.00256D0  !  epsilon_{0,a}
  !real(dp), parameter :: epsmax = -0.003253075D0  !  epsilon_{max} (not used)
  !real(dp), parameter :: epsmax = -1.1D-3  !  epsilon_{max} (not used)

  real(dp), parameter :: diter = 0.82D0  !  ITER strand diameter (mm)
  real(dp), parameter :: cuiter = 0.5D0  !  ITER strand copper fraction

  real(dp) :: bred, epssh, t, bc20eps, &
       tc0eps, bzero, strfun, jc1, jc2, jc3, scalefac

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  $\epsilon_{sh}$
  epssh = (ca2*eps0a)/(sqrt(ca1**2 - ca2**2))

  !  Strain function $s(\epsilon)$
  !  0.83 < s < 1.0, for -0.005 < strain < 0.005
  strfun = sqrt(epssh**2 + eps0a**2) - sqrt((strain-epssh)**2 + eps0a**2)
  strfun = strfun*ca1 - ca2*strain
  strfun = 1.0D0 + (1.0D0/(1.0D0 - ca1*eps0a))*strfun

  !  $B^*_{C2} (0,\epsilon)$
  bc20eps = bc20max*strfun

  !  $T^*_C (0,\epsilon)$
  tc0eps = tc0max * strfun**(1.0D0/3.0D0)

  !  Reduced temperature, restricted to be < 1
  !  Should remain < 1 for thelium < 0.94*tc0max (i.e. 15 kelvin for i_tf_sc_mat=1)

  if (thelium/tc0eps >= 1.0D0) then
     fdiags(1) = thelium ; fdiags(2) = tc0eps
     call report_error(159)
  end if
  t = min(thelium/tc0eps, 0.9999D0)

  !  Reduced magnetic field at zero temperature
  !  Should remain < 1 for bmax < 0.83*bc20max (i.e. 27 tesla for i_tf_sc_mat=1)

  if (bmax/bc20eps >= 1.0D0) then
     fdiags(1) = bmax ; fdiags(2) = bc20eps
     call report_error(160)
  end if
  bzero = min(bmax/bc20eps, 0.9999D0)

  !  Critical temperature (K)
  tcrit = tc0eps * (1.0D0 - bzero)**(1.0D0/1.52D0)  !  bzero must be < 1 to avoid NaNs

  !  Critical field (T)
  bcrit = bc20eps * (1.0D0 - t**1.52D0)

  !  Reduced magnetic field, restricted to be < 1
  if (bmax/bcrit >= 1.0D0) then
     fdiags(1) = bmax ; fdiags(2) = bcrit
     call report_error(161)
  end if
  bred = min(bmax/bcrit, 0.9999D0)

  !  Critical current density in superconductor (A/m2)
  !  ITER parameterization is for the current in a single strand,
  !  not per unit area, so scalefac converts to current density

  scalefac = pi * (0.5D0*diter)**2 * (1.0D0-cuiter)

  jc1 = (csc/bmax)*strfun
  jc2 = (1.0D0-t**1.52D0) * (1.0D0-t**2)  !  t must be < 1 to avoid NaNs
  jc3 = bred**p * (1.0D0-bred)**q  !  bred must be < 1 to avoid NaNs

  jcrit = jc1 * jc2 * jc3 / scalefac

end subroutine itersc

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bi2212(bmax,jstrand,tsc,fhts,jcrit,tmarg)

    !! Fitted parameterization to Bi-2212 superconductor properties
    !! author: P J Knight, CCFE, Culham Science Centre
    !! author: M Kovari, CCFE, Culham Science Centre
    !! bmax    : input real : Magnetic field at conductor (T)
    !! jstrand : input real : Current density in strand (A/m2)
    !! tsc     : input real : Superconductor temperature (K)
    !! fhts    : input real : Adjustment factor (<= 1) to account for strain,
    !! radiation damage, fatigue or AC losses
    !! jcrit : output real : Critical current density in strand (A/m2)
    !! tmarg : output real : Temperature margin (K)
    !! This routine calculates the critical current density and
    !! the temperature margin for Bi-2212 superconductor in the TF coils
    !! using a fit by M. Kovari to measurements described in the reference,
    !! specifically from the points shown in Figure 6.
    !! <P>Bi-2212 (Bi<SUB>2</SUB>Sr<SUB>2</SUB>CaCu<SUB>2</SUB>O<SUB>8-x</SUB>)
    !! is a first-generation high temperature superconductor; it still needs
    !! to be operated below about 10K, but remains superconducting at much
    !! higher fields at that temperature than Nb3Sn etc.
    !! The model's range of validity is T &lt; 20K, adjusted field
    !! b &lt; 104 T, B &gt; 6 T.
    !! A transformative superconducting magnet technology for fields well
    !! above 30 T using isotropic round wire multifilament
    !! Bi2Sr2CaCu2O8-x conductor, D. C. Larbalestier et al., preprint,
    !! 9th April 2013
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use error_handling, only: fdiags, report_error
    implicit none

    !  Arguments

    real(dp), intent(in) :: bmax, jstrand, tsc, fhts
    real(dp), intent(out) :: jcrit, tmarg

    !  Local variables

    real(dp) :: b

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Adjusted field (T)

    b = bmax / exp(-0.168D0*(tsc-4.2D0))

    !  Engineering (i.e. strand) critical current density (A/m2)

    jcrit = fhts * (1.175D9*exp(-0.02115D0*b) - 1.288D8)

    !  Temperature margin (K)
    !  Simple inversion of above calculation, using actual current density
    !  in strand instead of jcrit

    tmarg = 1.0D0/0.168D0 * &
    log( log(1.175D9/(jstrand/fhts + 1.288D8)) / (0.02115*bmax) ) &
    + 4.2D0 - tsc

    !  Check if ranges of validity have been violated

    if ((tsc > 20.0D0).or.(bmax < 6.0D0).or.(b > 104.0D0)) then
        fdiags(1) = tsc ; fdiags(2) = bmax ; fdiags(3) = b
        call report_error(106)
        write(*,*) 'Warning in routine BI2212:'
        write(*,*) 'Range of validity of the HTS Bi-2212 model has been violated:'
        write(*,*) '   S/C temperature (K) = ',tsc, ' (should be < 20 K)'
        write(*,*) 'Field at conductor (T) = ',bmax, ' (should be > 6 T)'
        write(*,*) '    Adjusted field (T) = ',b, ' (should be < 104 T)'
        write(*,*) ' '
    end if

end subroutine bi2212
!------------------------------------------------------------------
subroutine jcrit_nbti(temperature,bmax,c0,bc20max,tc0max,jcrit,tcrit)

    !! Critical current density in a NbTi superconductor strand
    !! author: P J Knight, CCFE, Culham Science Centre
    !! temperature : input real : SC temperature (K)
    !! bmax : input real : Magnetic field at conductor (T)
    !! c0   : input real : Scaling constant (A/m2)
    !! bc20max : input real : Upper critical field (T) for superconductor
    !! at zero temperature and strain
    !! tc0max : input real : Critical temperature (K) at zero field and strain
    !! jcrit : output real : Critical current density in superconductor (A/m2)
    !! tcrit : output real : Critical temperature (K)
    !! This routine calculates the critical current density and
    !! temperature in superconducting TF coils using NbTi
    !! as the superconductor.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: pi
    implicit none

    !  Arguments
    real(dp), intent(in) :: temperature, bmax, c0, bc20max, tc0max
    real(dp), intent(out) :: jcrit, tcrit

    !  Local variables
    real(dp) :: bratio, tbar
    !-----------------------------------
    bratio = bmax/bc20max

    if (bmax < bc20max) then
        !  Critical temperature (K)
        tcrit = tc0max * (1.0D0 - bratio)**0.59D0
    else
        ! Allow bmax > bc20max but set error flag
        ! Fudge to give real (negative) value if bratio < 1
        tcrit = tc0max * (1.0D0 - bratio)
    end if

    ! Allow tbar to be negative but set error flag
    tbar = 1.0D0 - temperature/tcrit

    !  Critical current density (A/m2)
    jcrit = c0 * (1.0D0 - bratio) * tbar

    ! if ((temperature > tcrit).or.(bmax > bc20max))then
    !     write(*,*)'jcrit_nbti: out of range: ', 'bmax =', bmax, ' bc20max =', bc20max, &
    !               ' temperature =',temperature, ' tcrit =',tcrit
    ! end if

end subroutine jcrit_nbti
!--------------------------------------------------------------------
subroutine GL_nbti(thelium,bmax,strain,bc20max,t_c0,jcrit,bcrit,tcrit)

    !!  Author: S B L Chislett-McDonald Durham University
    !!  Category: subroutine
    !!
    !!  Critical current density of the superconductor in an ITER
    !!  Nb-Ti strand based on the Ginzburg-Landau theory of superconductivity
    !!
    !! \begin{equation}
    !!  J_{c,TS}(B,T,\epsilon_{I}) = A(\epsilon_{I}) \left[T_{c}(\epsilon_{I})*(1-t^2)\right]^2\left
    !!  [B_{c2}(\epsilon_I)*(1-t^\nu)\right]^{n-3}b^{p-1}(1-b)^q~.
    !! \end{equation}
    !!
    !!  - \( \thelium \) -- Coolant/SC temperature [K]
    !!  - \( \bmax \) -- Magnetic field at conductor [T]
    !!  - \( \\epsilon_{I} \) -- Intrinsic strain in superconductor [\%]
    !!  - \( \B_{c2}(\epsilon_I) \) -- Strain dependent upper critical field [T]
    !!  - \( \b \) -- Reduced field = bmax / \B_{c2}(\epsilon_I)*(1-t^\nu) [unitless]
    !!  - \( \T_{c}(\epsilon_{I}) \) -- Strain dependent critical temperature (K)
    !!  - \( \t \) -- Reduced temperature = thelium / \T_{c}(\epsilon_{I}) [unitless]
    !!  - \( \A(epsilon_{I}) \) -- Strain dependent Prefactor [A / ( m\(^2\) K\(^-2) T\(^n-3))]
    !!  - \( \J_{c,TS} \) --  Critical current density in superconductor [A / m\(^-2\)]

    implicit none

    !  Arguments
    real(dp), intent(in) :: thelium, bmax, strain, bc20max, t_c0
    real(dp), intent(out) :: jcrit, tcrit, bcrit

    !  Local variables
    real(dp) :: strain_func, T_e, A_e, b_reduced, t_reduced, epsilon_I
    !critical current density prefactor (strand non-copper J_c)
    real(dp), parameter :: A_0 = 1102D6
    !flux pinning field scaling parameters
    real(dp), parameter :: p = 0.49D0
    real(dp), parameter :: q = 0.56D0
    real(dp), parameter :: n = 1.83D0
    !temperatute scaling parameter
    real(dp), parameter :: v = 1.42D0
    !strain scaling parameters
    real(dp), parameter :: c2 = -0.0025D0
    real(dp), parameter :: c3 = -0.0003D0
    real(dp), parameter :: c4 = -0.0001D0
    real(dp), parameter :: em = -0.002D-2
    !strain conversion parameters
    real(dp), parameter :: u = 0.0D0
    real(dp), parameter :: w = 2.2D0

    epsilon_I = strain - em

    strain_func = 1 + c2*(epsilon_I)**2 + c3*(epsilon_I)**3 + c4*(epsilon_I)**4

    T_e = t_c0 * strain_func**(1 / w)

    t_reduced = thelium/T_e

    A_e = A_0 * strain_func**(u / w)

    !  Critical Field
    bcrit = bc20max * (1 - t_reduced**v) * strain_func

    b_reduced = bmax/bcrit

    !  Critical temperature (K)
    tcrit = T_e

    !  Critical current density (A/m2)
    if (b_reduced <= 1.0D0) then
        jcrit = A_e * (T_e*(1-t_reduced**2))**2 * bcrit**(n-3) * b_reduced**(p-1) * (1 - b_reduced)**q
    else !Fudge to yield negative single valued function of Jc for fields above Bc2
        jcrit = A_e * (T_e*(1-t_reduced**2))**2 * bcrit**(n-3) * b_reduced**(p-1) * (1 - b_reduced)**1.0
    end if

end subroutine GL_nbti

subroutine wstsc(temperature,bmax,strain,bc20max,tc0max,jcrit,bcrit,tcrit)

    !! Implementation of WST Nb3Sn critical surface implementation
    !! author: J Morris, CCFE, Culham Science Centre
    !! temperature : input real : SC temperature (K)
    !! bmax : input real : Magnetic field at conductor (T)
    !! strain : input real : Strain in superconductor
    !! bc20max : input real : Upper critical field (T) for superconductor
    !! at zero temperature and strain
    !! tc0max : input real : Critical temperature (K) at zero field and strain
    !! jcrit : output real : Critical current density in superconductor (A/m2)
    !! bcrit : output real : Critical field (T)
    !! tcrit : output real : Critical temperature (K)
    !! This routine calculates the critical current density and
    !! temperature in the superconducting TF coils using the
    !! WST Nb3Sn critical surface model.
    !! V. Corato et al, "Common operating values for DEMO magnets design for 2016",
    !! https://scipub.euro-fusion.org/wp-content/uploads/eurofusion/WPMAGREP16_16565_submitted.pdf
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use error_handling, only: fdiags, report_error
    use maths_library, only: variable_error
    implicit none

    ! Arguments
    real(dp), intent(in) :: temperature, bmax, strain, bc20max, tc0max
    real(dp), intent(out) :: jcrit, bcrit, tcrit

    ! Local variables

    ! Scaling constant C [AT/mm2]
    real(dp), parameter :: csc = 83075.0D0
    ! Low field exponent p
    real(dp), parameter :: p = 0.593D0
    ! High field exponent q
    real(dp), parameter :: q = 2.156D0
    ! Strain fitting constant C_{a1}
    real(dp), parameter :: ca1 = 50.06D0
    ! Strain fitting constant C_{a2}
    real(dp), parameter :: ca2 = 0.0D0
    ! epsilon_{0,a}
    real(dp), parameter :: eps0a = 0.00312D0

    !real(dp), parameter :: epsmax = -1.1D-3  !  epsilon_{max} (not used)

    real(dp) :: bred, epssh, t, bc20eps, &
    tc0eps, bzero, strfun, jc1, jc2, jc3, scalefac

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  $\epsilon_{sh}$
    epssh = (ca2*eps0a)/(sqrt(ca1**2 - ca2**2))

    !  Strain function $s(\epsilon)$
    !  0.83 < s < 1.0, for -0.005 < strain < 0.005
    strfun = sqrt(epssh**2 + eps0a**2) - sqrt((strain-epssh)**2 + eps0a**2)
    strfun = strfun*ca1 - ca2*strain
    strfun = 1.0D0 + (1.0D0/(1.0D0 - ca1*eps0a))*strfun
    if(strfun<0.d0)write(*,*)'subroutine wstsc: strfun<0.d0. strfun =', strfun, ', strain = ',strain

    !  $B^*_{C2} (0,\epsilon)$
    bc20eps = bc20max*strfun

    !  $T^*_C (0,\epsilon)$
    tc0eps = tc0max * strfun**(1.0D0/3.0D0)

    !  Reduced temperature
    !  Should remain < 1 for temperature < 0.94*tc0max (i.e. 15 kelvin for i_tf_sc_mat=1)

    if (temperature/tc0eps >= 1.0D0) then
        fdiags(1) = temperature ; fdiags(2) = tc0eps
        call report_error(159)
    end if
    ! t = min(thelium/tc0eps, 0.9999D0)
    t = temperature/tc0eps

    !  Reduced magnetic field at zero temperature
    !  Should remain < 1 for bmax < 0.83*bc20max (i.e. 27 tesla for i_tf_sc_mat=1)

    if (bmax/bc20eps >= 1.0D0) then
        fdiags(1) = bmax ; fdiags(2) = bc20eps
        call report_error(160)
    end if

    ! bzero = min(bmax/bc20eps, 0.9999D0)
    bzero = bmax/bc20eps

    if (bzero < 1.0d0) then
        !  Critical temperature (K)
        tcrit = tc0eps * (1.0D0 - bzero)**(1.0D0/1.52D0)
    else
        ! Allow bzero > 1, fudge to give real (negative) value of tcrit
        ! This generates a real (negative) and continuous (but not differentiable)
        ! function of bzero.
        tcrit = tc0eps
    end if



    !  Critical field (T). Negative if normalised temperature t>1
    if(t>0.0d0)then
        bcrit = bc20eps * (1.0D0 - t**1.52D0)
    else
        ! Allow t<0, fudge to give real value of bcrit
        bcrit = bc20eps * (1.0D0 - t)
    end if

    !  Reduced magnetic field, restricted to be < 1
    if (bmax/bcrit >= 1.0D0) then
        fdiags(1) = bmax ; fdiags(2) = bcrit
        call report_error(161)
    end if
    ! bred = min(bmax/bcrit, 0.9999D0)
    bred = bmax/bcrit

    if ((bred>0.0d0).and.(bred < 1.0d0)) then
        jc3 = bred**p * (1.0D0-bred)**q  !  bred must be < 1 to avoid NaNs
    else
        ! Allow bred > 1 or <0, fudge to give real (negative) value of jc3
        ! This generates a real (negative) and continuous (but not differentiable)
        ! function of bred.
        jc3 = bred * (1.0D0-bred)
        if(variable_error(jc3))then
            write(*,'(a24, 8(a12,es12.3))')'jc3 jcrit is NaN.',' bred=',bred, ' bmax=',bmax, ' bcrit=',bcrit, ' t=',t
            stop 1
        end if
    end if

    !  Critical current density in superconductor (A/m2)
    jc1 = (csc/bmax)*strfun

    if(t>0.0d0)then
        jc2 = (1.0D0-t**1.52D0) * (1.0D0-t**2)
    else
        ! Allow t<0, fudge to give real value of jc2
        ! This generates a real and continuous (but not differentiable) function of t.
        jc2 = (1.0D0-t) * (1.0D0-t**2)
    end if

    ! jc3 = bred**p * (1.0D0-bred)**q  !  bred must be < 1 to avoid NaNs

    ! scale from mm2 to m2
    scalefac = 1.0D6

    jcrit = jc1 * jc2 * jc3*scalefac
    if(variable_error(jcrit))then
        write(*,'(a24, 8(a12,es12.3))')'WST jcrit is NaN.',' jc1=',jc1, ' jc2=',jc2, ' jc3=',jc3, ' t=',t
        write(*,'(a24, 8(a12,es12.3))')'T=',T,' bmax=',bmax,' strain=',strain,' bc20max=',bc20max, &
                                       ' tc0max=',tc0max,'jcrit=',jcrit,' bcrit=',bcrit,' tcrit=', tcrit
        stop 1
    end if

end subroutine wstsc
!--------------------------------------------------------------------------

subroutine GL_REBCO(thelium,bmax,strain,bc20max,t_c0,jcrit,bcrit,tcrit) !SCM added 13/06/18

  !!  Author: S B L Chislett-McDonald Durham University
  !!  Category: subroutine
  !!
  !!  Critical current density of a SuperPower REBCO tape based on measurements by P. Branch
  !!  at Durham University
  !!  https://git.ccfe.ac.uk/process/process/uploads/e98c6ea13da782cdc6fe16daea92078a/20200707_Branch-Osamura-Hampshire_-_accepted_SuST.pdf
  !!  and fit to state-of-the-art measurements at 4.2 K published in SuST
  !!  http://dx.doi.org/10.1088/0953-2048/24/3/035001
  !!
  !! \begin{equation}
  !!  J_{c,TS}(B,T,\epsilon_{I}) = A(\epsilon_{I}) \left[T_{c}(\epsilon_{I})*(1-t^2)\right]^2\left
  !!  [B_{c2}(\epsilon_I)*(1-t)^s\right]^{n-3}b^{p-1}(1-b)^q~.
  !! \end{equation}
  !!
  !!  - \( \thelium \) -- Coolant/SC temperature [K]
  !!  - \( \bmax \) -- Magnetic field at conductor [T]
  !!  - \( \\epsilon_{I} \) -- Intrinsic strain in superconductor [\%]
  !!  - \( \B_{c2}(\epsilon_I) \) -- Strain dependent upper critical field [T]
  !!  - \( \b \) -- Reduced field = bmax / \B_{c2}(\epsilon_I)*(1-t^\nu) [unitless]
  !!  - \( \T_{c}(\epsilon_{I}) \) -- Strain dependent critical temperature (K)
  !!  - \( \t \) -- Reduced temperature = thelium / \T_{c}(\epsilon_{I}) [unitless]
  !!  - \( \A(epsilon_{I}) \) -- Strain dependent Prefactor [A / ( m\(^2\) K\(^-2) T\(^n-3))]
  !!  - \( \J_{c,TS} \) --  Critical current density in superconductor [A / m\(^-2\)]
  !!  - \( \\epsilon_{m} \) -- Strain at which peak in J_c occurs [\%]


  implicit none

  !  Arguments
  real(dp), intent(in) :: thelium, bmax, strain, bc20max, t_c0
  real(dp), intent(out) :: jcrit, tcrit, bcrit

  !  Local variables
  real(dp) :: strain_func, T_e, A_e, b_reduced, t_reduced, epsilon_I
  !critical current density prefactor
  real(dp), parameter :: A_0 = 2.95D2
  !flux pinning field scaling parameters
  real(dp), parameter :: p = 0.32D0
  real(dp), parameter :: q = 2.50D0
  real(dp), parameter :: n = 3.33D0
  !temperatute scaling parameter
  real(dp), parameter :: s = 5.27D0
  !strain scaling parameters
  real(dp), parameter :: c2 = -0.0191D0
  real(dp), parameter :: c3 = 0.0039D0
  real(dp), parameter :: c4 = 0.00103D0
  real(dp), parameter :: em = 0.058D0
  !strain conversion parameters
  real(dp), parameter :: u = 0.0D0
  real(dp), parameter :: w = 2.2D0

  epsilon_I = strain - em

  strain_func = 1 + c2*(epsilon_I)**2 + c3*(epsilon_I)**3 + c4*(epsilon_I)**4

  T_e = t_c0 * strain_func**(1 / w)

  t_reduced = thelium/T_e

  A_e = A_0 * strain_func**(u / w)

  !  Critical Field
  bcrit = bc20max * (1 - t_reduced)**s * strain_func

  b_reduced = bmax/bcrit

  !  Critical temperature (K)
  tcrit = T_e

  !  Critical current density (A/m2)
  jcrit = A_e * (T_e*(1-t_reduced**2))**2 * bcrit**(n-3) * b_reduced**(p-1) * (1 - b_reduced)**q

end subroutine GL_REBCO

!--------------------------------------------------------------------------

subroutine HIJC_REBCO(thelium,bmax,strain,bc20max,t_c0,jcrit,bcrit,tcrit)

    !! Implementation of High Current Density REBCO tape
    !! author: R Chapman, UKAEA
    !! thelium : input real : SC temperature (K)
    !! bmax : input real : Magnetic field at conductor (T)
    !! strain : input real : Strain in superconductor
    !! bc20max : input real : Upper critical field (T) for superconductor
    !! at zero temperature and strain
    !! t_c0 : input real : Critical temperature (K) at zero field and strain
    !! jcrit : output real : Critical current density in superconductor (A/m2)
    !! bcrit : output real : Critical field (T)
    !! tcrit : output real : Critical temperature (K)
    !!
    !! Returns the critical current of a REBCO tape based on a critical surface
    !! (field, temperature) parameterization. Based in part on the parameterization
    !! described in: M. J. Wolf, N. Bagrets, W. H. Fietz, C. Lange and K. Weiss,
    !! "Critical Current Densities of 482 A/mm2 in HTS CrossConductors at 4.2 K and 12 T,"
    !! in IEEE Transactions on Applied Superconductivity, vol. 28, no. 4, pp. 1-4,
    !! June 2018, Art no. 4802404, doi: 10.1109/TASC.2018.2815767. And on the experimental
    !! data presented here: "2G HTS Wire Development at SuperPower", Drew W. Hazelton,
    !! February 16, 2017 https://indico.cern.ch/event/588810/contributions/2473740/
    !! The high Ic parameterization is a result of modifications based on Ic values
    !! observed in: "Conceptual design of HTS magnets for fusion nuclear science facility",
    !! Yuhu Zhai, Danko van der Laan, Patrick Connolly, Charles Kessel, 2021,
    !! https://doi.org/10.1016/j.fusengdes.2021.112611
    !! The parameter A is transformed into a function A(T) based on a Newton polynomial fit
    !! considering A(4.2 K) = 2.2e8, A(20 K) = 2.3e8 and A(65 K) = 3.5e8. These values were
    !! selected manually. A good fit to the pubished data can be seen in the 4-10 T range
    !! but the fit deviates at very low or very high field.
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use rebco_variables, only: tape_thickness, tape_width, rebco_thickness

    implicit none

    !  Arguments
    real(dp), intent(in) :: thelium, bmax, strain, bc20max, t_c0
    real(dp), intent(out) :: jcrit, tcrit, bcrit

    !  Local variables
    real(dp) :: A_t
    real(dp), parameter :: a = 1.4D0
    real(dp), parameter :: b = 2.005D0
    !critical current density prefactor
    real(dp), parameter :: A_0 = 2.2D8
    !flux pinning field scaling parameters
    real(dp), parameter :: p = 0.39D0
    real(dp), parameter :: q = 0.9D0
    !strain conversion parameters
    real(dp), parameter :: u = 33450.0D0
    real(dp), parameter :: v = -176577.0D0

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Critical Field (T)
    !  B_crit(T) calculated using temperature and critical temperature
    bcrit = bc20max * (1.0D0 - thelium/t_c0)**a

    ! Critical temperature (K)
    !  scaled to match behaviour in GL_REBCO routine,
    !  ONLY TO BE USED until a better suggestion is received
    tcrit = 0.999965D0 * t_c0

    ! finding A(T); constants based on a Newton polynomial fit to pubished data
    A_t = A_0 + ( u * thelium**2 ) + ( v * thelium )

    ! Critical current density (A/m2)

    jcrit = ( A_t / bmax ) * bcrit**b * ( bmax / bcrit )**p * ( 1 - bmax/bcrit )**q

    ! Jc times HTS area: default area is width 4mm times HTS layer thickness 1 um,
    ! divided by the tape area to provide engineering Jc per tape, then multiplied by fraction 0.4
    ! to reach the level of current density expected in the space where the tapes are wound in A/m^2!
    jcrit = jcrit * (tape_width * rebco_thickness) / (tape_width * tape_thickness) * 0.4D0

end subroutine HIJC_REBCO

!----------------------------------------------------------------

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

subroutine copper_properties(T,B, copper_resistivity, copper_rrr, copper_cp)
    ! Review of ROXIE's Material Properties Database for Quench Simulation,
    ! Author: Giulio Manfreda, December 2011
    ! https://espace.cern.ch/roxie/Documentation/Materials.pdf
    ! Different models use different definitions for residual resistivity ratio RRR.
    ! CUDI: resistivity at 290 K / 4 K.
    ! The range of validity of this t is between 4 K and 300 K.
    implicit none

    real(dp), intent(in) :: T, B   ! temperature, field
    real(dp), intent(out) :: copper_resistivity
    real(dp), intent(in) :: copper_rrr
    real(dp), intent(out) :: copper_cp
    real(dp):: bracket, logt, sum
    ! Fitting constants: resistivity
    real(dp), parameter:: t5 = 2.32547d9
    real(dp), parameter:: t3 = 9.57137d5
    real(dp), parameter:: t1 = 1.62735d2
    real(dp), parameter:: mr = 5.0d-11     ! ohm.m/T
    real(dp), parameter:: a = 1.7          ! ohm.m
    ! Fitting constants: specific heat (p.13)
    ! Checked against
    ! http://cryogenics.nist.gov/MPropsMAY/OFHC%20Copper/OFHC_Copper_rev1.htm
    real(dp), parameter::a0 = -1.91844d0
    real(dp), parameter::a1 = -0.15973d0
    real(dp), parameter::a2 =  8.61013d0
    real(dp), parameter::a3 = -18.996d0
    real(dp), parameter::a4 =  21.9661d0
    real(dp), parameter::a5 = -12.7328d0
    real(dp), parameter::a6 =  3.54322d0
    real(dp), parameter::a7 = -0.3797d0

    ! page 5: Copper resistivity is computed in CUDI with the t function similar
    ! to the one of McAshan [McA88]
    ! Note this formula is much quicker to evaluate than the NIST formula.

    bracket = 1d0 / (t5/T**5 + t3/T**3 + t1/T)

    copper_resistivity = 1.d-8 * (a/copper_rrr + bracket) + mr*B

    ! Specific heat
    ! NIST typical polynomial interpolation, equation 4 page 3
    logt = log10(T)
    sum = a0 + a1*logt + a2*logt**2 + a3*logt**3 + a4*logt**4 + a5*logt**5 + a6*logt**6 +a7*logt**7
    copper_cp = 10**sum

end subroutine copper_properties
! -------------------------------------------------------------------------
subroutine hastelloy_properties(temperature, hastelloy_cp)
    implicit none

    real(dp), intent(in) :: temperature   ! temperature
    real(dp), intent(out) :: hastelloy_cp
    real(dp) :: T

    T = temperature
    if(temperature>300d0) T=300.0d0

    ! Reinhard Heller: obtained by fitting data published in :
    ! J. Lu, E. S. Choi, and H. D. Zhou,
    ! "Physical properties of Hastelloy C-276 at cryogenic temperatures,"
    ! J. Appl. Phys. 103(6) 2008 064908

    if(T<42.2571d0)then
        hastelloy_cp = 0.60796d0 + 0.15309d0*T - 0.00237d0*T**2 + 6.76732d-4*T**3
    else if(T.ge.42.2571d0)then
        hastelloy_cp = -147.06251d0 + 5.43432d0*T - 0.01937d0*T**2 + 2.71669d-5*T**3   &
                       -8.12438d-9*T**4
    endif

    ! if(T<48.7135d0)then
    !     hastelloy%resistivity = 1.23386d-6 - 1.40462d-9*T + 1.09943d-10*T**2 -  &
    !                             3.81875d-12*T**3 + 6.3866d-14*T**4 - 4.10322d-16*T**5
    ! else if(T.ge.48.7135d0)then
    !     hastelloy%resistivity = 1.22641d-6 + 1.19188d-10*T
    ! endif
end subroutine hastelloy_properties
! --------------------------------------------------------------------------
subroutine solder_properties(T, solder_cp)
    use constants, only: pi
    implicit none

    real(dp), intent(in) :: T   ! temperature
    real(dp), intent(out) :: solder_cp

    ! Reinhard Heller: obtained by fitting data
    ! Material Database from Cryodata Software Package, CRYOCOMP, version 3.0, Florence, SC

    if(T<20.0d0)then
        solder_cp = 4.88028d0 - 2.92865d0*T + 0.52736d0*T**2 - 0.01861d0*T**3 + 2.36019d-4*T**4
    else if((T.ge.20.0d0).and.(T<300.0d0))then
        solder_cp = -10.15269d0 + 3.70087d0*T - 0.02947d0*T**2 + 1.02933d-4*T**3 - 1.29518d-7*T**4
    else
        solder_cp = 181.29d0
    endif
end subroutine solder_properties
! -------------------------------------------------------------------------
subroutine jacket_properties(T, jacket_cp)
    implicit none

    real(dp), intent(in) :: T   ! temperature
    real(dp), intent(out) :: jacket_cp
    real(dp):: logt, sum
    ! Fitting constants: specific heat (p.13)
    ! http://cryogenics.nist.gov/MPropsMAY/304Stainless/304Stainless_rev.htm
    real(dp), parameter::a0 = 22.0061
    real(dp), parameter::a1 = -127.5528
    real(dp), parameter::a2 = 303.647
    real(dp), parameter::a3 = -381.0098
    real(dp), parameter::a4 = 274.0328
    real(dp), parameter::a5 = -112.9212
    real(dp), parameter::a6 = 24.7593
    real(dp), parameter::a7 = -2.239153

    logt = log10(T)
    if(T>300) logt = log10(300d0)
    sum = a0 + a1*logt + a2*logt**2 + a3*logt**3 + a4*logt**4 + a5*logt**5 + &
                         a6*logt**6 + a7*logt**7
    jacket_cp = 10**sum

end subroutine jacket_properties
! ----------------------------------------------------------------------------
subroutine helium_properties(T, helium_cp_density)
    implicit none

    ! Isobaric Data for P = 0.60000 MPa
    ! http://webbook.nist.gov/chemistry/fluid/
    ! See very approximate fits in quench_data.xlsx (Issue #522)

    real(dp), intent(in) :: T   ! temperature
    real(dp), intent(out) :: helium_cp_density
    ! Cp x density J/K/m3
    if((T>=4d0).and.(T<10d0))then
        helium_cp_density = 12.285d3*T**3 - 309.92d3*T**2 + 2394.6d3*T - 5044.8d3
    else if(T>=10d0)then
        ! This probably works OK for arbitrarily high temperature
        helium_cp_density = 1745.1d3*T**(-1.031d0)
    else
        write(*,*)'temperature is below the range of helium data fit (>4 K): ', T
    end if

end subroutine helium_properties

! ---------------------------------------------------------------------------
!#####################################################################
!
!                  Sn-40Pb PROPERTIES PACKAGE
!                  ------------------------
! Contains functions for the calculation of the thermo-physical
! properties of Sn-40Pb
!
!#####################################################################

!#####################################################################
real function dSn40Pb(T)
    !#####################################################################
    !
    ! Density of Sn-40Pb
    !
    ! Range: 0 <= T <= inf K
    !
    !                        References
    !                        ----------
    ! Cryocomp v3.0
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   dSn40Pb         x    density                           Kg/m**3
    !
    !
    ! Author : L.Bottura at Cryosoft
    ! Version: 1.0   June 2015
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T
    ! fit variables

    ! local variables
    dSn40Pb = 8400.0

    return
! SJP Issue #835
! For Intel compliance add "end function"
end function dSn40Pb
!#####################################################################
real function cSn40Pb(T)
    !#####################################################################
    !
    ! Specific heat of Sn-40Pb
    !
    ! Range: 1 <= T <= 300 K
    !
    !                        References
    !                        ----------
    ! MISSING
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   cSn40Pb         x    specific heat                     J/Kg K
    !
    !
    ! Author : L.Bottura at Cryosoft
    ! Version: 1.0   June 2015
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T
    ! fit variables
    real     Tmin,Tmax,T0
    real     AA,BB,CC,DD,a,b,c,d,na,nb,nc,nd

    data     Tmin / 1.0 /, Tmax / 300.0/
    data     T0 / 3.96249314 /

    save
    ! local variables
    real     TT

    !
    TT=T
    TT=min(TT,Tmax)
    TT=max(TT,Tmin)

    if(TT.ge.1.and.TT.le.T0) then

        AA = -0.07848998
        BB = 0.2150118
        CC = -0.42863307
        DD = 0.06224209
        a = 1.89412709
        b = 556.143928
        c = 0.20145037
        d = 0.00086726
        na = -0.10237772
        nb = -0.06338306
        nc = 0.65668303
        nd = 0.31793398

    else

        AA = -0.0002679
        BB = -0.39223444
        CC = 1.11558347
        DD = -1.06506974
        a = 18.7914061
        b = 43.8271896
        c = 6.10665608
        d = 3.14878533
        na = 4.50255658
        nb = 1.22694566
        nc = 2.24090761
        nd = 3.38478332

    endif

    cSn40Pb = AA*TT**1/(1+TT/a)**na + BB*TT**2/(1+TT/b)**nb + &
    CC*TT**3/(1+TT/c)**nc + DD*TT**4/(1+TT/d)**nd

    return
! SJP Issue #835
! For Intel compliance add "end function"
end function cSn40Pb
!#####################################################################
real function kSn40Pb(T)
    !#####################################################################
    !
    ! Thermal conductivity of Sn-40Pb
    !
    ! Range: 2 <= T <= 300 K
    !
    !                        References
    !                        ----------
    ! MISSING
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   kSn40Pb          x    thermal conductivity              W/m K
    !
    !
    ! Author : L.Bottura at Cryosoft
    ! Version: 1.0   June 2015
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T
    ! fit variables
    real     Tmin,Tmax,T0
    real     rho0,rho273,RRR
    real     alpha,alpha2,beta
    real     L0,K0
    real     m,n,arg1,arg2
    real     W0,Wi,Wi0
    real     l

    data     Tmin / 2.0/, Tmax / 300.0/
    data     T0 / 22.702953 /
    data     rho273 / 1.50E-07 / , RRR / 26.4214669 /
    data     alpha2 / 2.1156885E-08 /
    data     L0 / 2.443E-08 / , K0 / 5.74E+01 /
    data     m / 3.092523 /, n / 3.0255307 /
    data     l / 0.20328615 /
    save
    ! local variables
    real     TT

    !
    TT=T
    TT=min(TT,Tmax)
    TT=max(TT,Tmin)

    rho0 = rho273/(RRR-1)
    beta = rho0/L0
    arg1 = (m-n)*(m+l)
    arg2 = (TT-T0)/T0
    alpha = alpha2*(beta/n/alpha2)**arg1
    W0 = beta/TT
    Wi = alpha*TT**n
    Wi0 = Wi+W0

    if(TT.le.T0) then

        kSn40Pb = 1/Wi0

    else
        kSn40Pb = 1/Wi0 + K0*(1-exp(-arg2))
    endif

    return
! SJP Issue #835
! For Intel compliance add "end function"
end function kSn40Pb
!#####################################################################
real function rSn40Pb(T)
    !#####################################################################
    !
    ! Electrical resistivity of Sn-40Pb
    !
    ! Range: 0 <= T <= 300 K
    !
    !                        References
    !                        ----------
    ! MISSING
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   rSn40Pb         x    resistivity                       Ohm m
    !
    ! Author : L.Bottura at CryoSoft
    ! Version: 1.0   June 2015
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T
    ! fit variables
    real     RRR,rho273,rho0,rhoi,rhoi0
    real     p1,p2,p3,p4,p5,p6,p7
    real     arg
    real     Tmin,Tmax

    data     RRR / 45 /, rho273 / 1.50E-07 /
    data     p1 / 3.98307546E-14 /, p2 / 3.47905728 /, &
    p3 / 3.6576564E+9 /, p4 / -1.10 /, &
    p5 / 0.0 /, p6 / 1.0 / , &
    p7 / 0.333 /
    data     Tmin   / 0.0/, Tmax   / 300.0/
    save
    ! local variables
    real     TT

    !
    TT=T
    TT=min(TT,Tmax)
    TT=max(TT,Tmin)

    rho0 = rho273/(RRR-1)
    arg = (p5/TT)**p6
    rhoi = p1*TT**p2/(1.+p1*p3*TT**(p2+p4)*exp(-arg))
    rhoi0 = (p7*rhoi*rho0)/(rhoi+rho0)

    rSn40Pb = rho0+rhoi+rhoi0

    !
    return
! SJP Issue #835
! For Intel compliance add "end function"
end function rSn40Pb

!#####################################################################
!
!                  HASTELLOY C276 PROPERTIES PACKAGE
!                  ---------------------------------
! Contains functions for the calculation of the thermo-physical
! properties of Hastelloy C276 (Haynes superalloy)
!
!#####################################################################

!#####################################################################
real function dHastelloyC276(T)
    !#####################################################################
    !
    ! Density of Hastelloy C276
    !
    ! Range: 0 <= T <= inf K
    !
    !                        References
    !                        ----------
    ! J.Lu, E.S. Choi, H.D. Zhou, Physical Properties of Hastelloy(R)
    ! C-276 at cryogenic temperatures, J. Appl. Phys, 103, 064908, 2008
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   dHastelloyC276    x    density                           Kg/m**3
    !
    !
    ! Author : L.Bottura at Cryosoft
    ! Version: 1.0   December 2012
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T
    ! fit variables

    ! local variables

    dHastelloyC276 = 8890.0

    return
! SJP Issue #835
! For Intel compliance add "end function"
end function dHastelloyC276
!#####################################################################
real(dp) function cHastelloyC276(T)
    !#####################################################################
    !
    ! Specific heat (not density) of Hastelloy C276
    !
    ! Range: 1 <= T <= 300 K
    !
    !                        References
    !                        ----------
    ! J.Lu, E.S. Choi, H.D. Zhou, Physical Properties of Hastelloy(R)
    ! C-276 at cryogenic temperatures, J. Appl. Phys, 103, 064908, 2008
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   cHastelloyC276    x    specific heat                     J/Kg K
    !
    !
    ! Author : L.Bottura at Cryosoft
    ! Version: 1.0   December 2012
    !
    !#####################################################################
    implicit none
    ! external variables
    real(dp)     T
    ! fit variables
    real(dp)     AA,BB,CC,DD,a,b,c,d,na,nb,nc,nd
    real(dp)     Tmin,Tmax
    data     AA  / 18.46314493  / ,  BB  / 1298.042986 / ,&
    CC  /-1105.076534  / ,  DD  / 2.226310361 /
    data     a   /  2.634486768 / ,  b   /  13.0796954 / ,&
    c   / 23.29856632  / ,  d   / 0.932886135 /
    data     na  /   14.010908  / ,  nb  / 12.68689825 / ,&
    nc  /  14.64960608 / ,  nd  / 2.626448071 /
    data     Tmin / 1.0/, Tmax / 300.0/
    save
    ! local variables
    real(dp)     TT

    TT=T
    TT=min(TT,Tmax)
    TT=max(TT,Tmin)

    cHastelloyC276 = AA*TT**na/(a+TT)**na + BB*TT**nb/(b+TT)**nb + &
    CC*TT**nc/(c+TT)**nc + DD*TT**nd/(d+TT)**nd

    return
! SJP Issue #835
! For Intel compliance add "end function"
end function cHastelloyC276
!#####################################################################
real function kHastelloyC276(T)
    !#####################################################################
    !
    ! Thermal conductivity of Hastelloy C276
    !
    ! Range: 2 <= T <= 300 K
    !
    !                        References
    !                        ----------
    ! J.Lu, E.S. Choi, H.D. Zhou, Physical Properties of Hastelloy(R)
    ! C-276 at cryogenic temperatures, J. Appl. Phys, 103, 064908, 2008
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   kHastelloyC276    x    thermal conductivity              W/m K
    !
    !
    ! Author : L.Bottura at Cryosoft
    ! Version: 1.0   December 2012
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T
    ! fit variables
    real     p1,p2,p3,p4
    real     Tmin,Tmax
    data     p1  /   0.2093 / ,  p2  / 125.7155 / ,&
    p3  /  60.0344 / ,  p4  / 2.1669   /
    data     Tmin / 2.0/, Tmax / 400.0/
    save
    ! local variables
    real     TT

    TT=T
    TT=min(TT,Tmax)
    TT=max(TT,Tmin)
    kHastelloyC276 = p1*TT*(1+(TT/p2)**p4) / (1+(TT/p3)**p4)

    return
! SJP Issue #835
! For Intel compliance add "end function"
end function kHastelloyC276
!#####################################################################
real function rHastelloyC276(T)
    !#####################################################################
    !
    ! Electrical resistivity of Hastelloy C276
    !
    ! Range: 0 <= T <= 400 K
    !
    !                        References
    !                        ----------
    ! J.Lu, E.S. Choi, H.D. Zhou, Physical Properties of Hastelloy(R)
    ! C-276 at cryogenic temperatures, J. Appl. Phys, 103, 064908, 2008
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   rHastelloyC276    x    resistivity                       Ohm m
    !
    ! Author : L.Bottura at CryoSoft
    ! Version: 1.0   October 2012
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T
    ! fit variables
    real     A,B,C
    real     Tmin,Tmax
    data     A  /  1.22634e-6 / ,  B  / 1.26121e-10 / ,&
    C  / -2.2417e-14 /
    data     Tmin   / 0.0/, Tmax   / 400.0/
    save
    ! local variables
    real     TT

    TT=T
    TT=min(TT,Tmax)
    TT=max(TT,Tmin)
    rHastelloyC276 = A + B*TT + C*TT*TT
    return
! SJP Issue #835
! For Intel compliance add "end function"
end function rHastelloyC276
!#####################################################################
!
!                   COPPER PROPERTIES PACKAGE
!                   -------------------------
!
! Contains functions for the calculation of the thermo-physical
! properties of pure Copper (Cu)
!
!#####################################################################

!#####################################################################
real function dCu(T)
    !#####################################################################
    !
    ! Density of Copper
    !
    ! Range: 0 <= T <= inf K
    !
    !                        References
    !                        ----------
    ! http://en.wikipedia.org/wiki/Copper
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   dCu               x    density                           Kg/m**3
    !
    !
    ! Author : L.Bottura at Cryosoft
    ! Version: 1.0   October 2012
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T
    ! fit variables

    ! local variables

    !
    dCu = 8960.0
    !

    return
! SJP Issue #835
! For Intel compliance add "end function"
end function dCu
!#####################################################################
real function cCu(T)
    !#####################################################################
    !
    ! Specific heat of Copper
    !
    ! Range: 1 <= T <= 1000 K
    !
    !                        References
    !                        ----------
    ! CryoComp version 3.0
    ! Handbook on Materials for S.C. Machinery, NBS Boulder (yellow
    ! book), 1977
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   cCu               x    specific heat                     J/Kg K
    !
    ! Author : L.Bottura & C. Rosso at Cryosoft
    ! Version: 3.5   October 2012
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T
    ! fit variables
    real     T0,T1
    real     Tmin,Tmax
    real     a1,a2,a3
    real     b0,b1,b2,b3,b4
    real     AA,CC,DD,a,c,d,na,nc,nd
    data     T0  /   10.4529369 / ,  T1  /  48.26583891 /
    data     a1  /   0.01188007 / ,  a2  /  -0.00050323 / , &
    &         a3  /   0.00075762 /
    data     b0  /  -5.03283229 / ,  b1  /   1.27426834 / , &
    &         b2  /  -0.11610718 / ,  b3  /   0.00522402 / ,&
    &         b4  /  -5.2996E-05 /
    data     AA  / -65.07570094 / ,  &
    &         CC  /  624.7552517 / ,  DD  /  0.529512119 /
    data     a   /  1.833505318 / , &
    &         c   /  16.55124429 / ,  d   / -0.000101401 /
    data     na  /  0.518553624 / ,&
    &         nc  /  2.855560719 / ,  nd  /  2.983928329 /
    data     Tmin / 1.0/, Tmax / 1000.0/
    save

    ! local variables
    real     TT

    !
    TT=T
    TT=max(TT,Tmin)
    TT=min(TT,Tmax)

    if (TT.le.t0) then
        cCu = a1*TT + a2*TT**2 + a3*TT**3
    elseif (TT.gt.t0 .and. TT.le.t1) then
        cCu = b0 + b1*TT+ b2*TT**2 + b3*TT**3 + b4*TT**4
    elseif(TT.gt.t1) then
        cCu = AA*TT   /(a+TT)**na +&
        &         CC*TT**3/(c+TT)**nc + DD*TT**4/(d+TT)**nd
    endif

    !
    return
! SJP Issue #835
! For Intel compliance add "end function"
end function cCu
!#####################################################################
real function kCu(T,B,RRR)
    !#####################################################################
    !
    ! Thermal conductivity of Copper
    !
    ! Range: 0.1 <= T <= 1000 K, 0 <= B <= 30 T,  1.5 <= RRR <= 3000
    !
    !                        References
    !                        ----------
    ! CryoComp version 3.0
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   B               x      magnetic field                    T
    !   RRR             x      residual resistivity ratio        -
    !   kCu               x    thermal conductivity              W/m K
    !
    ! Author : L.Bottura & C. Rosso at CryoSoft
    ! Version: 3.5   October 2012
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T,B,RRR
    ! fit variables
    real     Tmin,Tmax
    parameter(Tmin   = 0.1, Tmax   = 1000.0)
    real     Bmin,Bmax
    parameter(Bmin   = 0.0, Bmax   = 30.0)
    real     RRRmin,RRRmax
    parameter(RRRmin = 1.5, RRRmax = 3000.0)
    real     rho273
    parameter(rho273=1.54e-8)
    real     L0
    parameter(L0=2.443e-8)
    real     p1,p2,p3,p4,p5,p6
    parameter(p1=1.754e-8, p2=2.763, p3=1102.,&
    &          p4=-0.165  , p5=70.0 , p6=1.765)
    ! local variables
    real     TT,BB,R
    real     rhozero,beta,p7,arga,argb,argc,argd
    real     wt,wi0,wi,w0,wc,magr
    !real     magrCu

    !
    TT=T
    TT=max(TT,Tmin)
    TT=min(TT,Tmax)

    BB=B
    BB=max(BB,Bmin)
    BB=min(BB,Bmax)

    R =RRR
    R =max(R,RRRmin)
    R =min(R,RRRmax)

    !
    rhozero = rho273/(R-1.)

    beta    = rhozero/L0

    p7      = 0.838/(beta/0.0003)**0.1661
    arga    = min((log(TT/470.0)/0.70)**2,30.0)
    argb    = min((log(TT/ 87.0)/0.45)**2,30.0)
    argc    = min((log(TT/ 21.0)/0.50)**2,30.0)
    argd    = min((p5/TT)**p6,30.0)

    wc      = -0.00012*log(TT/420.0)*exp(-arga) -&
    &           0.00016*log(TT/ 73.0)*exp(-argb) -&
    &           0.00002*log(TT/ 18.0)*exp(-argc)
    w0      = beta/TT
    wi      = p1*TT**p2/(1.+p1*p3*TT**(p2+p4)*exp(-argd))+wc
    wi0     = p7*wi*w0/(wi+w0)
    wt      = w0+wi+wi0

    magr    = magrCu(TT,BB,R)

    kCu     = 1.0/(wt*magr)
    return
! SJP Issue #835
! For Intel compliance add "end function"
end function kCu
!#####################################################################
real function rCu(T,B,RRR)
    !#####################################################################
    !
    ! Electrical resistivity of Copper
    !
    ! Range: 0.1 <= T <= 1000 K, 0 <= B <= 30 T,  1.5 <= RRR <= 3000
    !
    !                        References
    !                        ----------
    ! CryoComp version 3.0
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   B               x      magnetic field                    T
    !   RRR             x      residual resistivity ratio        -
    !   rCu               x    resistivity                       Ohm m
    !
    ! Author : L.Bottura & C. Rosso at CryoSoft
    ! Version: 3.5   October 2012
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T,B,RRR
    ! fit variables
    real     Tmin,Tmax
    parameter(Tmin   = 0.1, Tmax   = 1000.0)
    real     Bmin,Bmax
    parameter(Bmin   = 0.0, Bmax   = 30.0)
    real     RRRmin,RRRmax
    parameter(RRRmin = 1.5, RRRmax = 3000.0)
    real     rho273
    parameter(rho273=1.54e-8)
    real     p1,p2,p3,p4,p5,p6,p7
    parameter(p1=0.1171e-16, p2=4.49, p3=3.841e10,&
    &          p4=-1.14, p5=50.0, p6=6.428, p7=0.4531)
    ! local variables
    real     TT,BB,R,rhozero,arg,rhoi,rhoi0,rho0,magr
    !real     magrCu

    TT=T
    TT=max(TT,Tmin)
    TT=min(TT,Tmax)

    BB=B
    BB=max(BB,Bmin)
    BB=min(BB,Bmax)

    R =RRR
    R =max(R,RRRmin)
    R =min(R,RRRmax)

    ! resistivity at absolute zero
    rhozero = rho273/(R-1.0)

    ! resistivity at the temperature t
    arg     = min((p5/TT)**p6,30.0)
    rhoi    = p1*TT**p2/(1.+p1*p3*TT**(p2+p4)*exp(-arg))
    rhoi0   = p7*rhoi*rhozero/(rhoi+rhozero)
    rho0    = rhozero+rhoi+rhoi0

    ! transverse magneto-resistance factor
    magr    = magrCu(TT,BB,R)

    ! resistivity
    rCu     = magr * rho0

    return
! SJP Issue #835
! For Intel compliance add "end function"
end function rCu
!#####################################################################
!
! Auxiliary functions and calculations
!
!#####################################################################

!#####################################################################
real function magrCu(T,B,RRR)
    !#####################################################################
    !
    ! Magneto-resistivity factor of Copper, given as the ratio between
    ! resistivity in transverse magnetic field B to resistivity at zero
    ! field
    !
    ! Range: 0.1 <= T <= 1000 K, 0 <= B <= 30 T,  1.5 <= RRR <= 3000
    !
    !                        References
    !                        ----------
    ! CryoComp version 3.0
    !
    ! variable          I/O    meaning                           units
    ! --------------------------------------------------------------------
    !   T               x      absolute temperature              K
    !   B               x      magnetic field                    T
    !   RRR             x      residual resistivity ratio        -
    !   magrCu            x    magneto-resistivity factor        -
    !
    ! Author : L.Bottura and C. Rosso at Cryosoft
    ! Version: 1.1   October 2012
    !
    !#####################################################################
    implicit none
    ! external variables
    real     T,B,RRR
    ! fit variables
    real     Tmin,Tmax
    parameter(Tmin   = 0.1, Tmax   = 1000.0)
    real     RRRmin,RRRmax
    parameter(RRRmin = 1.5, RRRmax = 3000.0)
    real     Bmin,Bmax
    parameter(Bmin   = 0.0, Bmax   = 30.0)
    real     brrmin,brrmax
    parameter(brrmin = 0.0, brrmax = 40.0e3)
    real     rho273,rhorrr
    parameter(rho273=1.54e-8,rhorrr=2.37e-8)
    real     p1,p2,p3,p4,p5,p6,p7
    parameter(p1=0.1171e-16, p2=4.49, p3=3.841e10,&
    &          p4=-1.14, p5=50.0, p6=6.428, p7=0.4531)
    real     a1,a2,a3,a4
    parameter(a1=0.382806e-03, a2=0.132407e+01, a3=0.167634e-02,&
    &          a4=0.789953e+00)
    ! local variables
    real     TT,BB,R,rhozero,rhoice,arg,rhoi,rhoi0,rho0,brr,magr

    TT=T
    TT=max(TT,Tmin)
    TT=min(TT,Tmax)

    BB=B
    BB=max(BB,Bmin)
    BB=min(BB,Bmax)

    R =RRR
    R =max(R,RRRmin)
    R =min(R,RRRmax)

    ! resistivity at absolute zero
    rhozero = rho273/(R-1.0)

    ! fit for the resistivity at the ice temperature
    rhoice  = rho273 + rhorrr/R

    ! resistivity at the temperature t
    arg     = min((p5/TT)**p6,30.0)
    rhoi    = p1*TT**p2/(1.+p1*p3*TT**(p2+p4)*exp(-arg))
    rhoi0   = p7*rhoi*rhozero/(rhoi+rhozero)
    rho0    = rhozero+rhoi+rhoi0

    ! product of field and residual resistivity ratio
    brr     = BB*rhoice/rho0
    brr     = max(brr,brrmin)
    brr     = min(brr,brrmax)

    ! fit for the transverse magneto-resistance increase
    if(brr .gt. 1.) then
        magr = a1*brr**a2/(1.0+a3*brr**a4)
    else
        magr = 0.0
    endif

    ! transverse magneto-resistance factor
    magrCu  = magr+1.0

    return
! SJP Issue #835
! For Intel compliance add "end function"
end function magrCu

!-----------------------------------------------------------------


end module superconductors
