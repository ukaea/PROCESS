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

end module superconductors
