 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physics_module

  !! Module containing tokamak plasma physics routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains all the primary plasma physics routines
  !! for a tokamak device.
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  !  Module-level variables

  integer :: iscz
  integer :: err242, err243
  real(dp) :: rad_fraction_LCFS
  real(dp) :: total_plasma_internal_energy  ! [J]
  real(dp) :: total_loss_power        ! [W]
  real(dp) :: total_energy_conf_time  ! [s]
  real(dp) :: ptarmw, lambdaio, drsep
  real(dp) :: fio, fLI, fLO, fUI, fUO, pLImw, pLOmw, pUImw, pUOmw
  real(dp) :: rho_star
  real(dp) :: nu_star
  real(dp) :: beta_mcdonald
  real(dp) :: itart_r

  ! Var in subroutine plasma_composition which requires re-initialisation on
  ! each new run
  integer :: first_call

  contains

  subroutine init_physics_module
    !! Initialise module variables
    implicit none

    first_call = 1
    iscz = 0
    err242 = 0
    err243 = 0
    rad_fraction_LCFS = 0.0D0
    total_plasma_internal_energy = 0.0D0
    total_loss_power = 0.0D0
    total_energy_conf_time = 0.0D0
    ptarmw = 0.0D0
    lambdaio = 0.0D0
    drsep = 0.0D0
    fio = 0.0D0
    fLI = 0.0D0
    fLO = 0.0D0
    fUI = 0.0D0
    fUO = 0.0D0
    pLImw = 0.0D0
    pLOmw = 0.0D0
    pUImw = 0.0D0
    pUOmw = 0.0D0
    rho_star   = 0.0D0
    nu_star   = 0.0D0
    beta_mcdonald = 0.0D0
    itart_r = 0.0D0
  end subroutine init_physics_module

  subroutine subr(a, b)
     implicit none
     real, intent(in) :: a
     real, intent(out) :: b
     b = a
  end subroutine subr

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function eped_warning()
      ! Issue #413.  MDK 26/2/18: improved output
    use physics_variables, only: rminor, tesep, triang, rmajor, kappa, &
      normalised_total_beta, plascur
    implicit none
      character(len=100) :: eped_warning, info_string
      eped_warning=''
      info_string = ''

      if((triang<0.399d0).or.(triang>0.601d0)) then
          write(info_string , '(1pe13.4)') triang
          eped_warning='triang = '//trim(info_string)
      endif
      if((kappa<1.499d0).or.(kappa>2.001d0)) then
          write(info_string , '(1pe13.4)') kappa
          eped_warning=trim(eped_warning)//'  kappa = '//trim(info_string)
      endif
      if((plascur<9.99d6).or.(plascur>20.01d6)) then
          write(info_string , '(1pe13.4)') plascur
          eped_warning=trim(eped_warning)//'  plascur = '//trim(info_string)
      endif
      if((rmajor<6.99d0).or.(rmajor>11.01d0)) then
          write(info_string , '(1pe13.4)') rmajor
          eped_warning=trim(eped_warning)//'  rmajor = '//trim(info_string)
      endif
      if((rminor<1.99d0).or.(rminor>3.501d0))then
          write(info_string , '(1pe13.4)') rminor
          eped_warning=trim(eped_warning)//'  rminor = '//trim(info_string)
      endif
      if((normalised_total_beta<1.99d0).or.(normalised_total_beta>3.01d0))then
          write(info_string , '(1pe13.4)') normalised_total_beta
          eped_warning=trim(eped_warning)//'  normalised_total_beta = '//trim(info_string)
      endif
      if(tesep>0.5)then
          write(info_string , '(1pe13.4)') tesep
          eped_warning=trim(eped_warning)//'  tesep = '//trim(info_string)
      endif
  end function eped_warning

 ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function bootstrap_fraction_iter89(aspect,beta,bt,cboot,plascur,q95,q0,rmajor,vol)

    !! Original ITER calculation of bootstrap-driven fraction
    !! of the plasma current.
    !! author: P J Knight, CCFE, Culham Science Centre
    !! aspect  : input real : plasma aspect ratio
    !! beta    : input real : plasma total beta
    !! bt      : input real : toroidal field on axis (T)
    !! cboot   : input real : bootstrap current fraction multiplier
    !! plascur : input real : plasma current (A)
    !! q95     : input real : safety factor at 95% surface
    !! q0      : input real : central safety factor
    !! rmajor  : input real : plasma major radius (m)
    !! vol     : input real : plasma volume (m3)
    !! This routine performs the original ITER calculation of the
    !! plasma current bootstrap fraction.
    !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
    !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    use constants, only: pi, rmu0
    implicit none

    real(dp) :: bootstrap_fraction_iter89

    !  Arguments

    real(dp), intent(in) :: aspect, beta, bt, cboot, &
         plascur, q95, q0, rmajor, vol

    !  Local variables

    real(dp) :: betapbs, bpbs, cbs, xbs, bootipf

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xbs = min(10.0D0, q95/q0)
    cbs = cboot * (1.32D0 - 0.235D0*xbs + 0.0185D0*xbs**2)
    bpbs = rmu0*plascur/(2.0D0*pi*sqrt(vol/(2.0D0* pi**2 *rmajor)) )
    betapbs = beta*bt**2 / bpbs**2

    if (betapbs <= 0.0D0) then  !  only possible if beta <= 0.0
       bootipf = 0.0D0
    else
       bootipf = cbs * ( betapbs/sqrt(aspect) )**1.3D0
    end if

    bootstrap_fraction_iter89 = bootipf

  end function bootstrap_fraction_iter89

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function bootstrap_fraction_nevins(alphan,alphat,betat,bt,dene,plascur, &
       q95,q0,rmajor,rminor,ten,zeff)

    !! Bootstrap current fraction from Nevins et al scaling
    !! author: P J Knight, CCFE, Culham Science Centre
    !! alphan : input real :  density profile index
    !! alphat : input real :  temperature profile index
    !! betat  : input real :  total plasma beta (with respect to the toroidal
    !! field)
    !! bt     : input real :  toroidal field on axis (T)
    !! dene   : input real :  electron density (/m3)
    !! plascur: input real :  plasma current (A)
    !! q0     : input real :  central safety factor
    !! q95    : input real :  safety factor at 95% surface
    !! rmajor : input real :  plasma major radius (m)
    !! rminor : input real :  plasma minor radius (m)
    !! ten    : input real :  density weighted average plasma temperature (keV)
    !! zeff   : input real :  plasma effective charge
    !! This function calculates the bootstrap current fraction,
    !! using the Nevins et al method, 4/11/90.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use physics_variables, only: te0, ne0
      use constants, only: rmu0, echarge
      use maths_library, only: quanc8
    implicit none

    real(dp) :: bootstrap_fraction_nevins

    !  Arguments

    real(dp), intent(in) :: alphan,alphat,betat,bt,dene,plascur, &
         q0,q95,rmajor,rminor,ten,zeff

    !  Local variables

    integer :: nofun
    real(dp) :: aibs,ainteg,betae0,dum1,fibs,flag

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Calculate peak electron beta

    betae0 = ne0 * te0 * 1.0D3*echarge / ( bt**2 /(2.0D0*rmu0) )

    !  Call integration routine

    call quanc8(bsinteg,0.0D0,0.999D0,0.001D0,0.001D0,ainteg,dum1, &
         nofun,flag)

    !  Calculate bootstrap current and fraction

    aibs = 2.5D0 * betae0 * rmajor * bt * q95 * ainteg
    fibs = 1.0D6 * aibs / plascur

    bootstrap_fraction_nevins = fibs

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function bsinteg(y)

      !! Integrand function for Nevins et al bootstrap current scaling
      !! author: P J Knight, CCFE, Culham Science Centre
      !! y : input real : abscissa of integration, = normalised minor radius
      !! This function calculates the integrand function for the
      !! Nevins et al bootstrap current scaling, 4/11/90.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: bsinteg

      !  Arguments

      real(dp), intent(in) :: y

      !  Local variables

      real(dp) :: alphai,al1,al2,a1,a2,betae,c1,c2,c3, &
           d,del,pratio,q,x,z

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Constants for fit to q-profile

      c1 = 1.0D0
      c2 = 1.0D0
      c3 = 1.0D0

      !  Compute average electron beta

      betae = dene*ten*1.0D3*echarge/(bt**2/(2.0D0*rmu0))

      del = rminor*sqrt(y)/rmajor
      x = (1.46D0*sqrt(del) + 2.4D0*del)/(1.0D0 - del)**1.5D0
      z = zeff
      d = 1.414D0*z + z*z + x*(0.754D0 + 2.657D0*z + 2.0D0*z*z) &
           + x*x*(0.348D0 + 1.243D0*z + z*z)
      al2 = -x*(0.884D0 + 2.074D0*z)/d
      a2 = alphat*(1.0D0-y)**(alphan+alphat-1.0D0)
      alphai = -1.172D0/(1.0D0 + 0.462D0*x)
      a1 = (alphan+alphat)*(1.0D0-y)**(alphan+alphat-1.0D0)
      al1 = x*(0.754D0+2.21D0*z+z*z+x*(0.348D0+1.243D0*z+z*z))/d

      !  q-profile

      q = q0 + (q95-q0)*(c1*y + c2*y*y + c3*y**3)/(c1+c2+c3)

      pratio = (betat - betae) / betae

      bsinteg = (q/q95)*(al1*(a1 + pratio*(a1+alphai*a2) ) + al2*a2 )

    end function bsinteg

  end function bootstrap_fraction_nevins

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function bootstrap_fraction_wilson(alphaj,alphap,alphat,betpth, &
       q0,qpsi,rmajor,rminor)

    !! Bootstrap current fraction from Wilson et al scaling
    !! author: P J Knight, CCFE, Culham Science Centre
    !! alphaj  : input real :  current profile index
    !! alphap  : input real :  pressure profile index
    !! alphat  : input real :  temperature profile index
    !! beta    : input real :  total beta
    !! betpth  : input real :  thermal component of poloidal beta
    !! q0      : input real :  safety factor on axis
    !! qpsi    : input real :  edge safety factor
    !! rmajor  : input real :  major radius (m)
    !! rminor  : input real :  minor radius (m)
    !! This function calculates the bootstrap current fraction
    !! using the numerically fitted algorithm written by Howard Wilson.
    !! AEA FUS 172: Physics Assessment for the European Reactor Study
    !! H. R. Wilson, Nuclear Fusion <B>32</B> (1992) 257
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use error_handling, only: fdiags, report_error
		use maths_library, only: linesolv
    implicit none

    real(dp) :: bootstrap_fraction_wilson

    !  Arguments

    real(dp), intent(in) :: alphaj,alphap,alphat,betpth, &
         q0,qpsi,rmajor,rminor

    !  Local variables

    integer :: i
    real(dp), dimension(12) :: a, b
    real(dp) :: aj,alfpnw,alftnw,eps1,r1,r2, &
         saj,seps1,sss,termj,termp,termt,term1,term2,z

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  alphap, alphat and alphaj are indices relevant to profiles of
    !  the form
    !             p = p0.(1-(r/a)**2)**alphap, etc.
    !
    !  Convert these indices to those relevant to profiles of the form
    !             p = p0.psi**alfpnw, etc.

    term1 = log(0.5D0)
    term2 = log(q0/qpsi)

    termp = 1.0D0-0.5D0**(1.0D0/alphap)
    termt = 1.0D0-0.5D0**(1.0D0/alphat)
    termj = 1.0D0-0.5D0**(1.0D0/alphaj)

    alfpnw = term1/log( log( (q0+(qpsi-q0)*termp)/qpsi )/term2)
    alftnw = term1/log( log( (q0+(qpsi-q0)*termt)/qpsi )/term2)
    aj     = term1/log( log( (q0+(qpsi-q0)*termj)/qpsi )/term2)

    !  Crude check for NaN errors or other illegal values...

    if ((aj /= aj).or.(alfpnw /= alfpnw).or.(alftnw /= alftnw).or.(aj <= 0.0D0)) then
       fdiags(1) = aj ; fdiags(2) = alfpnw ; fdiags(3) = alftnw ; fdiags(4) = aj
       call report_error(76)
    end if

    !  Ratio of ionic charge to electron charge

    z = 1.0D0

    !  Inverse aspect ratio: r2 = maximum plasma radius, r1 = minimum

    r2 = rmajor+rminor
    r1 = rmajor-rminor
    eps1 = (r2-r1)/(r2+r1)

    !  Coefficients fitted using least squares techniques

    saj = sqrt(aj)

    a(1)  =    1.41D0*(1.0D0-0.28D0*saj)*(1.0D0+0.12D0/z)
    a(2)  =    0.36D0*(1.0D0-0.59D0*saj)*(1.0D0+0.8D0/z)
    a(3)  =   -0.27D0*(1.0D0-0.47D0*saj)*(1.0D0+3.0D0/z)
    a(4)  =  0.0053D0*(1.0D0+5.0D0/z)
    a(5)  =   -0.93D0*(1.0D0-0.34D0*saj)*(1.0D0+0.15D0/z)
    a(6)  =   -0.26D0*(1.0D0-0.57D0*saj)*(1.0D0-0.27D0*z)
    a(7)  =   0.064D0*(1.0D0-0.6D0*aj+0.15D0*aj*aj)*(1.0D0+7.6D0/z)
    a(8)  = -0.0011D0*(1.0D0+9.0D0/z)
    a(9)  =   -0.33D0*(1.0D0-aj+0.33D0*aj*aj)
    a(10) =   -0.26D0*(1.0D0-0.87D0/saj-0.16D0*aj)
    a(11) =   -0.14D0*(1.0D0-1.14D0/saj-0.45D0*saj)
    a(12) = -0.0069D0

    seps1 = sqrt(eps1)

    b(1)  = 1.0D0
    b(2)  = alfpnw
    b(3)  = alftnw
    b(4)  = alfpnw*alftnw
    b(5)  = seps1
    b(6)  = alfpnw*seps1
    b(7)  = alftnw*seps1
    b(8)  = alfpnw*alftnw*seps1
    b(9)  = eps1
    b(10) = alfpnw*eps1
    b(11) = alftnw*eps1
    b(12) = alfpnw*alftnw*eps1

    sss = 0.0D0
    do i = 1,12
       sss = sss + a(i)*b(i)
    end do

    !  Empirical bootstrap current fraction

    bootstrap_fraction_wilson = seps1 * betpth * sss

  end function bootstrap_fraction_wilson

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function bootstrap_fraction_sauter()

    !! Bootstrap current fraction from Sauter et al scaling
    !! author: P J Knight, CCFE, Culham Science Centre
    !! None
    !! This function calculates the bootstrap current fraction
    !! using the Sauter, Angioni and Lin-Liu scaling.
    !! <P>The code was supplied by Emiliano Fable, IPP Garching
    !! (private communication).
    !! O. Sauter, C. Angioni and Y. R. Lin-Liu,
    !! Physics of Plasmas <B>6</B> (1999) 2834
    !! O. Sauter, C. Angioni and Y. R. Lin-Liu, (ERRATA)
    !! Physics of Plasmas <B>9</B> (2002) 5140
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use physics_variables, only: dnitot, rminor, tesep, ti, triang, q0, afuel, &
      zeff, rhopedn, bt, plascur, xarea, fhe3, teped, dene, te, rmajor, q, &
      nesep, te0, neped, tbeta, ne0, alphan, rhopedt, alphat
		use profiles_module, only: tprofile, nprofile
		use constants, only: pi
    implicit none

    real(dp) :: bootstrap_fraction_sauter

    !  Arguments

    !  Local variables

    integer, parameter :: nr = 200
    integer :: ir
    real(dp) :: da,drho,iboot,jboot,roa
    real(dp) :: dlogne_drho,dlogte_drho,dlogti_drho
    real(dp), dimension(nr) :: amain,mu,ne,ni,rho,sqeps,tempe,tempi,zef,zmain

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Populate profile arrays

    do ir = 1,nr
       roa = dble(ir)/nr
       rho(ir) = sqrt(xarea/pi) * roa !  local circularised minor radius (m)
       sqeps(ir) = sqrt(roa * rminor/rmajor)

       ne(ir) = 1.0D-19 * nprofile(roa,rhopedn,ne0,neped,nesep,alphan)
       ni(ir) = dnitot/dene * ne(ir)
       tempe(ir) = tprofile(roa,rhopedt,te0,teped,tesep,alphat,tbeta)
       tempi(ir) = ti/te * tempe(ir)

       zef(ir) = zeff  !  Flat Zeff profile assumed

       !  mu = 1/safety factor
       !  Parabolic q profile assumed

       mu(ir) = 1.0D0 / (q0 + (q-q0)*roa**2)
       amain(ir) = afuel  !  fuel ion mass
       zmain(ir) = 1.0D0 + fhe3  !  sum(Zi.ni)/sum(ni) over fuel ions i
    end do

    !  Ensure that density and temperature values are not zero at edge

    if (ne(nr) == 0.0D0) then
       ne(nr) = 1.0D-4*ne(nr-1)
       ni(nr) = 1.0D-4*ni(nr-1)
    end if

    if (tempe(nr) == 0.0D0) then
       tempe(nr) = 1.0D-4*tempe(nr-1)
       tempi(nr) = 1.0D-4*tempi(nr-1)
    end if

    !  Calculate total bootstrap current (MA) by summing along profiles

    iboot = 0.0D0
    do ir = 1,nr

       if (ir == nr) then
          jboot = 0.0D0
          da = 0.0D0
       else
          drho = rho(ir+1) - rho(ir)
          da = 2.0D0*pi*rho(ir)*drho  !  area of annulus

          dlogte_drho = (log(tempe(ir+1)) - log(tempe(ir))) / drho
          dlogti_drho = (log(tempi(ir+1)) - log(tempi(ir))) / drho
          dlogne_drho = (log(ne(ir+1)) - log(ne(ir))) / drho

          !  The factor of 0.5 below arises because in ASTRA the logarithms
          !  are coded as (e.g.):  (Te(j+1)-Te(j))/(Te(j+1)+Te(j)), which
          !  actually corresponds to grad(log(Te))/2. So the factors dcsa etc.
          !  are a factor two larger than one might otherwise expect.

          jboot = 0.5D0 * ( dcsa(ir,nr) * dlogne_drho &
               + hcsa(ir,nr) * dlogte_drho &
               + xcsa(ir,nr) * dlogti_drho )
          jboot = -bt/(0.2D0*pi*rmajor) * rho(ir)*mu(ir) * jboot  !  MA/m2
       end if

       iboot = iboot + da*jboot  !  MA

    end do

    bootstrap_fraction_sauter = 1.0D6 * iboot/plascur

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function beta_poloidal_local(j,nr)

      !! Local beta poloidal calculation
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! nr : input integer : maximum value of j
      !! This function calculates the local beta poloidal.
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! <P>beta poloidal = 4*pi*ne*Te/Bpo**2
      !! Pereverzev, 25th April 1989 (?)
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: beta_poloidal_local

      !  Arguments

      integer, intent(in) :: j, nr

      !  Local variables

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (j /= nr)  then
         beta_poloidal_local = 1.6D-4*pi * (ne(j+1)+ne(j)) * (tempe(j+1)+tempe(j))
      else
         beta_poloidal_local = 6.4D-4*pi * ne(j)*tempe(j)
      end if

      beta_poloidal_local = beta_poloidal_local * &
           ( rmajor/(bt*rho(j)*abs(mu(j)+1.0D-4)) )**2

    end function beta_poloidal_local

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function beta_poloidal_local_total(j,nr)

      !! Local beta poloidal calculation, including ion pressure
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! nr : input integer : maximum value of j
      !! This function calculates the local total beta poloidal.
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! <P>beta poloidal = 4*pi*(ne*Te+ni*Ti)/Bpo**2
      !! where ni is the sum of all ion densities (thermal)
      !! Pereverzev, 25th April 1989 (?)
      !! E Fable, private communication, 15th May 2014
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: beta_poloidal_local_total

      !  Arguments

      integer, intent(in) :: j, nr

      !  Local variables

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (j /= nr)  then
         beta_poloidal_local_total = 1.6D-4*pi * ( &
              ( (ne(j+1)+ne(j)) * (tempe(j+1)+tempe(j)) ) + &
              ( (ni(j+1)+ni(j)) * (tempi(j+1)+tempi(j)) ) )
      else
         beta_poloidal_local_total = 6.4D-4*pi * (ne(j)*tempe(j) + ni(j)*tempi(j))
      end if

      beta_poloidal_local_total = beta_poloidal_local_total * &
           ( rmajor/(bt*rho(j)*abs(mu(j)+1.0D-4)) )**2

    end function beta_poloidal_local_total

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function nues(j)

      !! Relative frequency of electron collisions
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! This function calculates the relative frequency of electron
      !! collisions: <I>NU* = Nuei*q*Rt/eps**1.5/Vte</I>
      !! The electron-ion collision frequency NUEI=NUEE*1.4*ZEF is
      !! used.
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! Yushmanov, 30th April 1987 (?)
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: nues

      !  Arguments

      integer, intent(in) :: j

      !  Local variables

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nues = nuee(j) * 1.4D0*zef(j)*rmajor / &
           abs(mu(j)*(sqeps(j)**3)*sqrt(tempe(j))*1.875D7)

    end function nues

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function nuee(j)

      !! Frequency of electron-electron collisions
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! This function calculates the frequency of electron-electron
      !! collisions (Hz): <I>NUEE = 4*SQRT(pi)/3*Ne*e**4*lambd/
      !! SQRT(Me)/Te**1.5</I>
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! Yushmanov, 25th April 1987 (?),
      !! updated by Pereverzev, 9th November 1994 (?)
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: nuee

      !  Arguments

      integer, intent(in) :: j

      !  Local variables

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nuee = 670.0D0 * coulg(j) * ne(j) / ( tempe(j)*sqrt(tempe(j)) )

    end function nuee

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function coulg(j)

      !! Coulomb logarithm
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! This function calculates the Coulomb logarithm, valid
      !! for e-e collisions (T_e > 0.01 keV), and for
      !! e-i collisions (T_e > 0.01*Zeff^2) (Alexander, 9/5/1994).
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! C. A. Ordonez and M. I. Molina, Phys. Plasmas <B>1</B> (1994) 2515
      !! Rev. Mod. Phys., V.48, Part 1 (1976) 275
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: coulg

      !  Arguments

      integer, intent(in) :: j

      !  Local variables

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      coulg = 15.9D0 - 0.5D0*log(ne(j)) + log(tempe(j))

    end function coulg

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function nuis(j)

      !! Relative frequency of ion collisions
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! This function calculates the relative frequency of ion
      !! collisions: <I>NU* = Nui*q*Rt/eps**1.5/Vti</I>
      !! The full ion collision frequency NUI is used.
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! Yushmanov, 30th April 1987 (?)
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: nuis

      !  Arguments

      integer, intent(in) :: j

      !  Local variables

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      nuis = 3.2D-6 * nui(j)*rmajor / ( abs(mu(j)+1.0D-4) * &
           sqeps(j)**3 * sqrt(tempi(j)/amain(j)) )

    end function nuis

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function nui(j)

      !! Full frequency of ion collisions
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! This function calculates the full frequency of ion
      !! collisions (Hz).
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! None
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real(dp) :: nui

      !  Arguments

      integer, intent(in) :: j

      !  Local variables

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !	Coulomb logarithm = 15 is used

      nui = zmain(j)**4 * ni(j) * 322.0D0 / ( tempi(j)*sqrt(tempi(j)*amain(j)) )

    end function nui

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function dcsa(j,nr)

      !! Grad(ln(ne)) coefficient in the Sauter bootstrap scaling
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! nr : input integer : maximum value of j
      !! This function calculates the coefficient scaling grad(ln(ne))
      !! in the Sauter bootstrap current scaling.
      !! Code by Angioni, 29th May 2002.
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! O. Sauter, C. Angioni and Y. R. Lin-Liu,
      !! Physics of Plasmas <B>6</B> (1999) 2834
      !! O. Sauter, C. Angioni and Y. R. Lin-Liu, (ERRATA)
      !! Physics of Plasmas <B>9</B> (2002) 5140
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  DCSA $\equiv \mathcal{L}_{31}$, Eq.14a, Sauter et al, 1999

      implicit none

      real(dp) :: dcsa

      !  Arguments

      integer, intent(in) :: j,nr

      !  Local variables

      real(dp) :: zz,zft,zdf

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (j == 1) then
         dcsa = 0.0D0
      else
         zz = zef(j)
         zft = tpf(j)
         zdf = 1.0D0 + (1.0D0 - 0.1D0*zft)*sqrt(nues(j))
         zdf = zdf + 0.5D0*(1.0D0-zft)*nues(j)/zz
         zft = zft/zdf  !  $f^{31}_{teff}(\nu_{e*})$, Eq.14b
         dcsa = (1.0D0 + 1.4D0/(zz+1.0D0))*zft - 1.9D0/(zz+1.0D0)*zft*zft
         dcsa = dcsa + (0.3D0*zft*zft + 0.2D0*zft*zft*zft)*zft / (zz+1.0D0)

         !  Corrections suggested by Fable, 15/05/2015
         !dcsa = dcsa*beta_poloidal_local(j,nr) * (1.0D0+tempi(j)/(zz*tempe(j)))
         dcsa = dcsa*beta_poloidal_local_total(j,nr)
      end if

    end function dcsa

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function hcsa(j,nr)

      !! Grad(ln(Te)) coefficient in the Sauter bootstrap scaling
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! nr : input integer : maximum value of j
      !! This function calculates the coefficient scaling grad(ln(Te))
      !! in the Sauter bootstrap current scaling.
      !! Code by Angioni, 29th May 2002.
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! O. Sauter, C. Angioni and Y. R. Lin-Liu,
      !! Physics of Plasmas <B>6</B> (1999) 2834
      !! O. Sauter, C. Angioni and Y. R. Lin-Liu, (ERRATA)
      !! Physics of Plasmas <B>9</B> (2002) 5140
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  HCSA $\equiv ?$, Sauter et al, 1999

      implicit none

      real(dp) :: hcsa

      !  Arguments

      integer, intent(in) :: j,nr

      !  Local variables

      real(dp) :: zz,zft,zdf,zfte,zfte2,zfte3,zfte4
      real(dp) :: zfti,zfti2,zfti3,zfti4,hcee,hcei

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (j == 1) then
         hcsa = 0.0D0
      else
         zz = zef(j)
         zft = tpf(j)
         zdf = 1.0D0 + 0.26D0*(1.0D0-zft)*sqrt(nues(j))
         zdf = zdf + 0.18D0*(1.0D0-0.37D0*zft)*nues(j)/sqrt(zz)
         zfte = zft/zdf  !  $f^{32\_ee}_{teff}(\nu_{e*})$, Eq.15d
         zfte2 = zfte*zfte
         zfte3 = zfte*zfte2
         zfte4 = zfte2*zfte2

         zdf = 1.0D0 + (1.0D0 + 0.6D0*zft)*sqrt(nues(j))
         zdf = zdf + 0.85D0*(1.0D0 - 0.37D0*zft)*nues(j)*(1.0D0+zz)
         zfti = zft/zdf  !  $f^{32\_ei}_{teff}(\nu_{e*})$, Eq.15e
         zfti2 = zfti*zfti
         zfti3 = zfti*zfti2
         zfti4 = zfti2*zfti2

         hcee = (0.05D0 + 0.62D0*zz) / zz / (1.0D0 + 0.44D0*zz) * (zfte-zfte4)
         hcee = hcee + (zfte2 - zfte4 - 1.2D0*(zfte3-zfte4)) / (1.0D0 + 0.22D0*zz)
         hcee = hcee + 1.2D0/(1.0D0 + 0.5D0*zz)*zfte4  !  $F_{32\_ee}(X)$, Eq.15b

         hcei = -(0.56D0 + 1.93D0*zz) / zz / (1.0D0 + 0.44*zz) * (zfti-zfti4)
         hcei = hcei + 4.95D0/(1.0D0 + 2.48D0*zz) * &
              (zfti2 - zfti4 - 0.55D0*(zfti3-zfti4))
         hcei = hcei - 1.2D0/(1.0D0 + 0.5D0*zz)*zfti4  !  $F_{32\_ei}(Y)$, Eq.15c

         !  Corrections suggested by Fable, 15/05/2015
         !hcsa = beta_poloidal_local(j,nr)*(hcee + hcei) + dcsa(j,nr) &
         !     / (1.0D0 + tempi(j)/(zz*tempe(j)))
         hcsa = beta_poloidal_local(j,nr)*(hcee + hcei) + dcsa(j,nr) &
              * beta_poloidal_local(j,nr)/beta_poloidal_local_total(j,nr)
      end if

    end function hcsa

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function xcsa(j,nr)

      !! Grad(ln(Ti)) coefficient in the Sauter bootstrap scaling
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! nr : input integer : maximum value of j
      !! This function calculates the coefficient scaling grad(ln(Ti))
      !! in the Sauter bootstrap current scaling.
      !! Code by Angioni, 29th May 2002.
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! O. Sauter, C. Angioni and Y. R. Lin-Liu,
      !! Physics of Plasmas <B>6</B> (1999) 2834
      !! O. Sauter, C. Angioni and Y. R. Lin-Liu, (ERRATA)
      !! Physics of Plasmas <B>9</B> (2002) 5140
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: xcsa

      !  Arguments

      integer, intent(in) :: j,nr

      !  Local variables

      real(dp) :: zz,zft,zdf,a0,alp,a1,zfte

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (j == 1) then
         xcsa = 0.0D0
      else
         zz = zef(j)
         zft = tpf(j)
         zdf = 1.0D0 + (1.0D0 - 0.1D0*zft)*sqrt(nues(j))
         zdf = zdf + 0.5D0*(1.0D0 - 0.5D0*zft)*nues(j)/zz
         zfte = zft/zdf  !  $f^{34}_{teff}(\nu_{e*})$, Eq.16b

         xcsa = (1.0D0 + 1.4D0/(zz+1.0D0))*zfte - 1.9D0/(zz+1.0D0)*zfte*zfte
         xcsa = xcsa + (0.3D0*zfte*zfte + 0.2D0*zfte*zfte*zfte)*zfte &
              / (zz+1.0D0)  !  Eq.16a

         a0 = -1.17D0*(1.0D0-zft)
         a0 = a0 / (1.0D0 - 0.22D0*zft - 0.19D0*zft*zft)  !  $\alpha_0$, Eq.17a

         alp = (a0 + 0.25D0*(1.0D0 - zft*zft)*sqrt(nuis(j))) / &
              (1.0D0 + 0.5*sqrt(nuis(j)))
         a1 = nuis(j)*nuis(j) * zft**6
         alp = (alp + 0.315D0*a1) / (1.0D0 + 0.15D0*a1)  !  $\alpha(\nu_{i*})$, Eq.17b

         !  Corrections suggested by Fable, 15/05/2015
         !xcsa = beta_poloidal_local(j,nr) * (xcsa*alp)*tempi(j)/zz/tempe(j)
         !xcsa = xcsa + dcsa(j,nr) / (1.0D0 + zz*tempe(j)/tempi(j))

         xcsa = (beta_poloidal_local_total(j,nr)-beta_poloidal_local(j,nr)) &
              * (xcsa*alp)
         xcsa = xcsa + dcsa(j,nr) * &
              (1.0D0 - beta_poloidal_local(j,nr)/beta_poloidal_local_total(j,nr))
      end if

    end function xcsa

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function tpf(j)

      !! Trapped particle fraction
      !! author: P J Knight, CCFE, Culham Science Centre
      !! j  : input integer : radial element index in range 1 to nr
      !! This function calculates the trapped particle fraction at
      !! a given radius.
      !! <P>A number of different fits are provided, but the one
      !! to be used is hardwired prior to run-time.
      !! <P>The code was supplied by Emiliano Fable, IPP Garching
      !! (private communication).
      !! O. Sauter et al, Plasma Phys. Contr. Fusion <B>44</B> (2002) 1999
      !! O. Sauter, 2013:
      !! http://infoscience.epfl.ch/record/187521/files/lrp_012013.pdf
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: tpf

      !  Arguments

      integer, intent(in) :: j

      !  Local variables

      integer, parameter :: ASTRA=1, SAUTER2002=2, SAUTER2013=3

      real(dp) :: eps,epseff,g,s,zz

      integer, parameter :: fit = ASTRA

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      s = sqeps(j)
      eps = s*s

      select case (fit)

      case (ASTRA)

         !  ASTRA method, from Emiliano Fable, private communication
         !  (Excluding h term which dominates for inverse aspect ratios < 0.5,
         !  and tends to take the trapped particle fraction to 1)

         zz = 1.0D0 - eps

         g = 1.0D0 - zz*sqrt(zz) / (1.0D0 + 1.46D0*s)

         !  Advised by Emiliano to ignore ASTRA's h below
         !
         !h = 0.209D0 * (sqrt(tempi(j)*amain(j))/zmain(j)*mu(j)*rmajor*bt)**0.3333D0
         !tpf = min(1.0D0, max(g, h))

         tpf = g

      case (SAUTER2002)

         !  Equation 4 of Sauter 2002
         !  Similar to, but not quite identical to g above

         tpf = 1.0D0 - (1.0D0-eps)**2 / (1.0D0 + 1.46D0*s) / sqrt(1.0D0 - eps*eps)

      case (SAUTER2013)

         !  Includes correction for triangularity

         epseff = 0.67D0*(1.0D0 - 1.4D0*triang*abs(triang)) * eps

         tpf = 1.0D0 - sqrt( (1.0D0-eps)/(1.0D0+eps) ) * &
              (1.0D0 - epseff) / (1.0D0 + 2.0D0*sqrt(epseff))

      end select

    end function tpf

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !subroutine fast_alpha_bs()

      !  BSALP (local per index J) is in MA/m^2

      !  Required... before we can use this routine:
      !  fast alpha pressure profile
      !  poloidal flux profile vs local minor radius  (and grad(psi))
      !  Shafranov shift vs local minor radius

      !  all lengths in meters,
      !  temperatures in keV where j is the radial index,
      !  IPOL is R*Bphi / (R0*Bphi0)  (i.e. the normalized poloidal current integral)
      !  PFAST is the alpha pressure
      !  TE is the electron temperature
      !  SHIF is the Shafranov shift (defined with respect to the geom. major radius)
      !  AMETR is the minor radius
      !  RTOR = R0
      !  ZEF = Z effective
      !  FP = PSI (magnetic flux, poloidal) defined such that
      !    B_pol = grad(PSI) / (2*PI*R)

      ! ZBIRTH = 1.

      ! !MeV already included in PFAST ,convert PFAST from keV*1e19 to J
      ! ZDPDPSI = 3./2.*1.60218*1.e3* &
      !      (PFAST(J)-PFAST(J-1))/((FP(J)-FP(J-1))/GP2)
      ! ZSB=(0.027*ZEF(J-1)*(TE(J-1)/20.)**(3./2.))**(1./3.)
      ! ZSB1=(0.027*ZEF(J)*(TE(J)/20.)**(3./2.))**(1./3.)
      ! ZSC=(5./3.)**(1./3.)*ZSB
      ! ZSC1=(5./3.)**(1./3.)*ZSB1
      ! ZDSC3DPSI = 3./2.*1.60218*1.e3*PFAST(J)* &
      !      (ZSC1**3.-ZSC**3.)/((FP(J)-FP(J-1))/GP2)
      ! ZEPS=AMETR(J)/RTOR
      ! ZFP=1.-1.46*(1.+0.67/ZEF(J))*ZEPS**0.5+ 0.46*(1.+2.1/ZEF(J))*ZEPS

      ! ZDR0DR=(SHIF(J)-SHIF(J-1))/(AMETR(J)-AMETR(J-1))

      ! ZY=(1.-ZDR0DR/ZEPS*(1.-(1.-ZEPS**2.)**0.5)) &
      !      /(1.+ZEPS*ZDR0DR/2.)/(1.-ZEPS**2.)**.5

      ! ZA11=-ZSB**(3./2.)*(0.12+2.24*ZSB**(3./2.)-0.9*ZSB**3.) &
      !      /(0.7+18.4*ZSB**(3./2.)+ &
      !      23.5*ZSB**3.+101.*ZSB**(9./2.))*ZY

      ! ZA12=-2./3.*(0.5+0.8*ZSC**(3./2.)+0.4*ZSC**3.) &
      !      /(1.+2.3*ZSC**(3./2.)+4.*ZSC**3.)

      ! ZA21=(7.e-4+0.02*ZSB**(3./2.)+0.4*ZSB**3.)/ &
      !      (0.01-0.61*ZSB**(3./2.)+ &
      !      24.8*ZSB**3.-53.4*ZSB**(9./2.)+118.*ZSB**6.)*ZY

      ! ZA22=2./3.*(0.1+3.25*ZSC**(3./2.)-1.1*ZSC**3.)/ &
      !      (1.e-3+0.6*ZSC**(3./2.)+ &
      !      8.6*ZSC**3.+3.1*ZSC**(9./2.)+15.1*ZSC**6.)

      ! ZB1=(0.155+3.9*ZSB**(3./2.)-3.1*ZSB**3.+0.3*ZSB**6.) &
      !      /(0.1+3.*ZSB**(3./2.)-2.1*ZSB**3.)*ZY

      ! ZB2=(1.3-0.5*ZSB**(3./2.)+5.9*ZSB**3.)/ &
      !      (1.-0.34*ZSB**(3./2.)+4.9*ZSB**3.)*ZY

      ! ZA1 = -ZA11+(2.*ZA11+(2.*ZB1-3.)*ZA12)*ZEPS**.5 &
      !      -(ZA11+2.*(ZB1-1.)*ZA12)*ZEPS

      ! ZA2 = -ZA21+(2.*ZA21+(2.*ZB2-3.)*ZA22)*ZEPS**.5 &
      !      -(ZA21+2.*(ZB2-1.)*ZA22)*ZEPS

      ! !bootstrap current by alphas
      ! BSALP=-ZEPS**.5*(1.-2./ZEF(J)*ZFP)*IPOL(J)* &
      !      RTOR*ZBIRTH*(ZA1*ZDPDPSI+ZA2*ZDSC3DPSI)
      ! BSALP=BSALP/1.e6

    !end subroutine fast_alpha_bs

  end function bootstrap_fraction_sauter

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine diamagnetic_fraction_hender(beta,diacf)

    !! author: S.I. Muldrew, CCFE, Culham Science Centre
    !! Diamagnetic contribution at tight aspect ratio.
    !! Tim Hender fit
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   !  Arguments

   real(dp), intent(in) ::  beta
   real(dp), intent(out) :: diacf

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   diacf = beta / 2.8D0


  end subroutine diamagnetic_fraction_hender

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine diamagnetic_fraction_scene(beta,q95,q0,diacf)

    !! author: S.I. Muldrew, CCFE, Culham Science Centre
    !! Diamagnetic fraction based on SCENE fit by Tim Hender
    !! See Issue #992
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    real(dp), intent(in) ::  beta, q95, q0
    real(dp), intent(out) :: diacf

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    diacf = beta * (0.1D0*q95/q0+0.44D0) * 4.14D-1

  end subroutine diamagnetic_fraction_scene

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ps_fraction_scene(beta,pscf)

    !! author: S.I. Muldrew, CCFE, Culham Science Centre
    !! Pfirsch-Schlüter fraction based on SCENE fit by Tim Hender
    !! See Issue #992

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    real(dp), intent(in) ::  beta
    real(dp), intent(out) :: pscf

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    pscf = -9.0D-2 * beta

  end subroutine ps_fraction_scene

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine culcur(alphaj,alphap,bt,eps,icurr,iprofile,kappa,kappa95, &
       p0,pperim,q0,qpsi,rli,rmajor,rminor,sf,triang,triang95,bp,qstar,plascur)

    !! Routine to calculate the plasma current
    !! author: P J Knight, CCFE, Culham Science Centre
    !! alphaj   : input/output real : current profile index
    !! alphap   : input real :  pressure profile index
    !! bt       : input real :  toroidal field on axis (T)
    !! eps      : input real :  inverse aspect ratio
    !! icurr    : input integer : current scaling model to use
    !! 1 = Peng analytic fit
    !! 2 = Peng divertor scaling (TART)
    !! 3 = simple ITER scaling
    !! 4 = revised ITER scaling
    !! 5 = Todd empirical scaling I
    !! 6 = Todd empirical scaling II
    !! 7 = Connor-Hastie model
    !! iprofile : input integer : switch for current profile consistency
    !! 0 = use input alphaj, rli
    !! 1 = make these consistent with q, q0
    !! kappa    : input real :  plasma elongation
    !! kappa95  : input real :  plasma elongation at 95% surface
    !! p0       : input real :  central plasma pressure (Pa)
    !! pperim   : input real :  plasma perimeter length (m)
    !! q0       : input real :  plasma safety factor on axis
    !! qpsi     : input real :  plasma edge safety factor (= q-bar for icurr=2)
    !! rli      : input/output real : plasma normalised internal inductance
    !! rmajor   : input real :  major radius (m)
    !! rminor   : input real :  minor radius (m)
    !! sf       : input real :  shape factor for icurr=1 (=A/pi in documentation)
    !! triang   : input real :  plasma triangularity
    !! triang95 : input real :  plasma triangularity at 95% surface
    !! bp       : output real : poloidal field (T)
    !! qstar    : output real : equivalent cylindrical safety factor (shaped)
    !! plascur  : output real : plasma current (A)
    !! This routine performs the calculation of the
    !! plasma current, with a choice of formula for the edge
    !! safety factor. It will also make the current profile parameters
    !! consistent with the q-profile if required.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !! J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
    !! unpublished internal Oak Ridge document
    !! Y.-K. M. Peng, J. Galambos and P.C. Shipe, 1992,
    !! Fusion Technology, 21, 1729
    !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
    !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
    !! M. Kovari et al, 2014, "PROCESS": A systems code for fusion power plants -
    !! Part 1: Physics https://www.sciencedirect.com/science/article/pii/S0920379614005961
    !! H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO
    !! https://iopscience.iop.org/article/10.1088/0029-5515/53/7/073019
    !! T. Hartmann, 2013, Development of a modular systems code to analyse the
    !! implications of physics assumptions on the design of a demonstration fusion power plant
    !! https://inis.iaea.org/search/search.aspx?orig_q=RN:45031642
    !! Sauter, Geometric formulas for systems codes..., FED 2016
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use error_handling, only: idiags, report_error
		use physics_variables, only: normalised_total_beta, beta
		use global_variables, only: run_tests
		use constants, only: pi, rmu0
    implicit none

    !  Arguments

    integer, intent(in) :: icurr, iprofile
    real(dp), intent(inout) :: alphaj, rli
    real(dp), intent(in) :: alphap, bt, eps, kappa, &
         kappa95, p0, pperim, q0, qpsi, rmajor, rminor, sf, triang, triang95
    real(dp), intent(out) :: bp, qstar, plascur

    !  Local variables

    real(dp) :: asp, curhat, fq, w07

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Aspect ratio

    asp = 1.0D0/eps

    !  Calculate the function Fq that scales the edge q from the
    !  circular cross-section cylindrical case

    !  First check for negative triangularity using unsuitable current scaling

    if ((icurr.ne.8).and.(triang.lt.0.0)) then
     write(*,*) 'Triangularity is negative without icurr = 8.'
     write(*,*) 'Please check and try again.'
     write(*,*) 'PROCESS stopping'
     stop 1
    end if

    select case (icurr)

    case (1)  !  Peng analytical fit
       fq = (1.22D0-0.68D0*eps)/((1.0D0-eps*eps)**2) * sf**2

    case (2)  !  Peng scaling for double null divertor; TARTs [STAR Code]
       curhat = 1.0D6 * plasc(qpsi,asp,rminor,bt,kappa,triang)/bt

    case (3)  !  Simple ITER scaling (simply the cylindrical case)
       fq = 1.0D0

    case (4)  !  ITER formula (IPDG89)
       fq = 0.5D0 * (1.17D0-0.65D0*eps)/((1.0D0-eps*eps)**2) * &
            (1.0D0 + kappa95**2 * &
            (1.0D0 + 2.0D0*triang95**2 - 1.2D0*triang95**3) )

    case (5, 6) !  Todd empirical scalings

       fq = (1.0D0+2.0D0*eps*eps) * 0.5D0*(1.0D0+kappa95**2) * &
            (1.24D0-0.54D0*kappa95+0.3D0*(kappa95**2 + triang95**2) + &
            0.125D0*triang95)

       if (icurr == 6) fq = fq * (1.0D0 + ( abs(kappa95-1.2D0) )**3)

    case (7)  !  Connor-Hastie asymptotically-correct expression

       !  N.B. If iprofile=1, alphaj will be wrong during the first call (only)

       call conhas(alphaj,alphap,bt,triang95,eps,kappa95,p0,fq)

    case (8)  !  Sauter scaling allowing negative triangularity [FED May 2016]

        ! Assumes zero squareness, note takes kappa, delta at separatrix not _95

        w07 = 1.0d0    ! zero squareness - can be modified later if required

        fq = (1.0d0 + 1.2d0*(kappa - 1.0d0) + 0.56d0*(kappa-1.0d0)**2) * &
             (1.0d0 + 0.09d0 * triang + 0.16d0 * triang**2) * &
       (1.0d0 + 0.45d0 * triang * eps)/(1.0d0 - 0.74d0 * eps) * &
       (1.0d0 + 0.55d0 * (w07 - 1.0d0))

    case (9) ! FIESTA ST fit

       fq = 0.538D0 * (1.0D0 + 2.440D0*eps**2.736D0) * kappa**2.154D0 * triang**0.060D0

    case default
       idiags(1) = icurr ; call report_error(77)

    end select

    !  Calculate the ratio of plasma current to toroidal field

    if (icurr /= 2) then
       curhat = 5.0D6 * rminor**2 / (rmajor*qpsi) * fq
    end if
    if (icurr == 8) then
       curhat = 4.1d6 * rminor**2 / (rmajor*qpsi) * fq
    end if

    !  Calculate the equivalent edge safety factor (= qcyl)

    qstar = 5.0D6 * rminor**2 / (rmajor*curhat) * 0.5D0 * &
         (1.0D0 + kappa95**2 * &
         (1.0D0 + 2.0D0*triang95**2 - 1.2D0*triang95**3) )

    !  Calculate plasma current

    plascur = curhat * bt

    normalised_total_beta = 1.0D8*beta*rminor*bt/plascur

    !  Calculate the poloidal field

    bp = bpol(icurr,plascur,qpsi,asp,bt,kappa,triang,pperim)

    !  Ensure current profile consistency, if required
    !  This is as described in Hartmann and Zohm only if icurr = 4 as well...

    if (iprofile == 1) then
       alphaj = qstar/q0 - 1.0D0
       rli = log(1.65D0 + 0.89D0*alphaj)  !  Tokamaks 4th Edition, Wesson, page 116
    end if

  contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function plasc(qbar,aspect,rminor,bt,kappa,delta)

      !! Function to calculate plasma current (Peng scaling)
      !! author: J Galambos, FEDC/ORNL
      !! author: P J Knight, CCFE, Culham Science Centre
      !! aspect : input real :  plasma aspect ratio
      !! bt     : input real :  toroidal field on axis (T)
      !! delta  : input real :  plasma triangularity
      !! kappa  : input real :  plasma elongation
      !! qbar   : input real :  edge q-bar
      !! rminor : input real :  plasma minor radius (m)
      !! This function calculates the plasma current in MA,
      !! using a scaling from Peng, Galambos and Shipe (1992).
      !! It is primarily used for Tight Aspect Ratio Tokamaks and is
      !! selected via <CODE>icurr=2</CODE>.
      !! J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
      !! unpublished internal Oak Ridge document
      !! Y.-K. M. Peng, J. Galambos and P.C. Shipe, 1992,
      !! Fusion Technology, 21, 1729
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: plasc

      !  Arguments

      real(dp), intent(in) :: aspect,bt,delta,kappa,qbar,rminor

      !  Local variables

      real(dp) :: c1,c2,d1,d2,eps,e1,e2,f1,f2,ff1,ff2,g,h1,h2,y1,y2

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      eps = 1.0D0/aspect

      c1 = kappa**2/(1.0D0+delta) + delta
      c2 = kappa**2/(1.0D0-delta) - delta

      d1 = (kappa/(1.0D0+delta))**2 + 1.0D0
      d2 = (kappa/(1.0D0-delta))**2 + 1.0D0

      if (aspect < c1) then
         y1 = sqrt( (c1*eps - 1.0D0)/(1.0D0+eps) ) * (1.0D0 + delta)/kappa
      else
         y1 = sqrt( (1.0D0 - c1*eps)/(1.0D0+eps) ) * (1.0D0 + delta)/kappa
      end if
      y2 = sqrt( (c2*eps+1.0D0)/(1.0D0-eps) ) * (1.0D0-delta)/kappa

      e1 = 2.0D0*kappa/(d1*(1.0D0+delta))
      e2 = 2.0D0*kappa/(d2*(1.0D0-delta))

      h2 = (1.0D0 + (c2-1.0D0)*eps/2.0D0) / &
           sqrt( (1.0D0-eps)*(c2*eps+1.0D0) )
      f2 = (d2*(1.0D0-delta)*eps) / ( (1.0D0-eps)*(c2*eps+1.0D0) )
      g = eps*kappa / (1.0D0 - eps*delta)
      ff2 = f2 * (g + 2.0D0*h2*atan(y2) )

      if (aspect < c1) then
         h1 = (1.0D0 + (1.0D0-c1)*eps/2.0D0) / &
              sqrt( (1.0D0+eps)*(c1*eps-1.0D0) )
         f1 = (d1*(1.0D0+delta)*eps) / ( (1.0D0+eps)*(c1*eps-1.0D0) )
         ff1 = f1*(g - h1*log( (1.0D0+y1)/(1.0D0-y1) ) )
      else
         h1 = (1.0D0 + (1.0D0-c1)*eps/2.0D0) / &
              sqrt( (1.0D0+eps)*(1.0D0-c1*eps) )
         f1 = -(d1*(1.0D0+delta)*eps) / ( (1.0D0+eps)*(c1*eps-1.0D0) )
         ff1 = f1*( -g + 2.0D0*h1*atan(y1) )
      end if

      plasc = rminor*bt/qbar * 5.0D0*kappa/(2.0D0*pi**2) * &
           ( asin(e1)/e1 + asin(e2)/e2 ) * (ff1 + ff2)

    end function plasc

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine conhas(alphaj,alphap,bt,delta95,eps,kappa95,p0,fq)

      !! Routine to calculate the F coefficient used for scaling the
      !! plasma current
      !! author: P J Knight, CCFE, Culham Science Centre
      !! alphaj   : input real :  current profile index
      !! alphap   : input real :  pressure profile index
      !! bt       : input real :  toroidal field on axis (T)
      !! delta95  : input real :  plasma triangularity 95%
      !! eps      : input real :  inverse aspect ratio
      !! kappa95  : input real :  plasma elongation 95%
      !! p0       : input real :  central plasma pressure (Pa)
      !! fq       : output real : scaling for edge q from circular
      !! cross-section cylindrical case
      !! This routine calculates the F coefficient used for scaling the
      !! plasma current, using the Connor-Hastie scaling given in
      !! AEA FUS 172.
      !! AEA FUS 172: Physics Assessment for the European Reactor Study
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !  Arguments

      real(dp), intent(in) :: alphaj,alphap,bt,delta95,eps,kappa95,p0
      real(dp), intent(out) :: fq

      !  Local variables

      real(dp) :: beta0, deltap, deltar, eprime, er, kap1, &
           lambda, lamp1, li, nu, tprime, tr

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Exponent in Connor-Hastie current profile - matching total
      !  current gives the following trivial relation

      lambda = alphaj

      !  Exponent in Connor-Hastie pressure profile

      nu = alphap

      !  Central plasma beta

      beta0 = 2.0D0 * rmu0 * p0 / (bt**2)

      !  Plasma internal inductance

      lamp1 = 1.0D0 + lambda
      li = lamp1/lambda * (lamp1/lambda * log(lamp1) - 1.0D0)

      !  T/r in AEA FUS 172

      kap1 = kappa95 + 1.0D0
      tr = kappa95 * delta95 / kap1**2

      !  E/r in AEA FUS 172

      er = (kappa95-1.0D0)/kap1

      !  T primed in AEA FUS 172

      tprime = 2.0D0 * tr * lamp1/(1.0D0 + 0.5D0*lambda)

      !  E primed in AEA FUS 172

      eprime = er * lamp1/(1.0D0 + lambda/3.0D0)

      !  Delta primed in AEA FUS 172

      deltap = 0.5D0*kap1 * eps * 0.5D0*li + &
           beta0/(0.5D0*kap1*eps) * lamp1**2 / (1.0D0+nu)

      !  Delta/R0 in AEA FUS 172

      deltar = beta0/6.0D0 * (1.0D0 + 5.0D0*lambda/6.0D0 + 0.25D0*lambda**2) &
           + (0.5D0*kap1*eps)**2 * 0.125D0*(1.0D0-(lambda**2)/3.0D0)

      !  F coefficient

      fq = (0.5D0*kap1)**2 * &
           ( 1.0D0 + eps**2 * (0.5D0*kap1)**2 + 0.5D0*deltap**2 + &
           2.0D0*deltar + 0.5D0*(eprime**2 + er**2) + &
           0.5D0*(tprime**2 + 4.0D0*tr**2) )

    end subroutine conhas

  end subroutine culcur

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function bpol(icurr,ip,qbar,aspect,bt,kappa,delta,perim)

    !! Function to calculate poloidal field
    !! author: J Galambos, FEDC/ORNL
    !! author: P J Knight, CCFE, Culham Science Centre
    !! icurr  : input integer : current scaling model to use
    !! ip     : input real :  plasma current (A)
    !! qbar   : input real :  edge q-bar
    !! aspect : input real :  plasma aspect ratio
    !! bt     : input real :  toroidal field on axis (T)
    !! kappa  : input real :  plasma elongation
    !! delta  : input real :  plasma triangularity
    !! perim  : input real :  plasma perimeter (m)
    !! This function calculates the poloidal field in Tesla,
    !! using a simple calculation using Stoke's Law for conventional
    !! tokamaks, or for TARTs, a scaling from Peng, Galambos and
    !! Shipe (1992).
    !! J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
    !! unpublished internal Oak Ridge document
    !! Y.-K. M. Peng, J. Galambos and P.C. Shipe, 1992,
    !! Fusion Technology, 21, 1729
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use constants, only: pi, rmu0
    implicit none

    real(dp) :: bpol

    !  Arguments

    integer, intent(in) :: icurr
    real(dp), intent(in) :: aspect,bt,delta,ip,kappa,perim,qbar

    !  Local variables

    real(dp) :: c1,c2,d1,d2,eps,f1,f2,ff1,ff2,g,h1,h2,y1,y2

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (icurr /= 2) then

       !  Stoke's Law

       bpol = rmu0 * ip / perim

    else  !  Original coding, only suitable for TARTs [STAR Code]

       eps = 1.0D0/aspect

       c1 = kappa**2/(1.0D0+delta) + delta
       c2 = kappa**2/(1.0D0-delta) - delta

       d1 = (kappa/(1.0D0+delta))**2 + 1.0D0
       d2 = (kappa/(1.0D0-delta))**2 + 1.0D0

       if (aspect < c1) then
          y1 = sqrt( (c1*eps - 1.0D0)/(1.0D0+eps) ) * (1.0D0 + delta)/kappa
       else
          y1 = sqrt( (1.0D0 - c1*eps)/(1.0D0+eps) ) * (1.0D0 + delta)/kappa
       end if
       y2 = sqrt( (c2*eps+1.0D0)/(1.0D0-eps) ) * (1.0D0-delta)/kappa

       h2 = (1.0D0 + (c2-1.0D0)*eps/2.0D0) / &
            sqrt( (1.0D0-eps)*(c2*eps+1.0D0) )
       f2 = (d2*(1.0D0-delta)*eps) / ( (1.0D0-eps)*(c2*eps+1.0D0) )
       g = eps*kappa / (1.0D0 - eps*delta)
       ff2 = f2 * (g + 2.0D0*h2*atan(y2) )

       if (aspect < c1) then
          h1 = (1.0D0 + (1.0D0-c1)*eps/2.0D0) / &
               sqrt( (1.0D0+eps)*(c1*eps-1.0D0) )
          f1 = (d1*(1.0D0+delta)*eps) / ( (1.0D0+eps)*(c1*eps-1.0D0) )
          ff1 = f1*(g - h1*log( (1.0D0+y1)/(1.0D0-y1) ) )
       else
          h1 = (1.0D0 + (1.0D0-c1)*eps/2.0D0) / &
               sqrt( (1.0D0+eps)*(1.0D0-c1*eps) )
          f1 = -(d1*(1.0D0+delta)*eps) / ( (1.0D0+eps)*(c1*eps-1.0D0) )
          ff1 = f1*( -g + 2.0D0*h1*atan(y1) )
       end if

       bpol = bt * (ff1 + ff2) / (2.0D0 * pi * qbar)

    end if

  end function bpol

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine culblm(bt,dnbeta,plascur,rminor,betalim)

    !! Beta scaling limit
    !! author: P J Knight, CCFE, Culham Science Centre
    !! bt      : input real :  toroidal B-field on plasma axis (T)
    !! dnbeta  : input real :  Troyon-like g coefficient
    !! plascur : input real :  plasma current (A)
    !! rminor  : input real :  plasma minor axis (m)
    !! betalim : output real : beta limit as defined below
    !! This subroutine calculates the beta limit, using
    !! the algorithm documented in AEA FUS 172.
    !! <P>The limit applies to beta defined with respect to the total B-field.
    !! Switch ICULBL determines which components of beta to include (see
    !! routine <A HREF="constraints.html">constraints</A> for coding):
    !! <UL>
    !! <P><LI>If ICULBL = 0, then the limit is applied to the total beta
    !! <P><LI>If ICULBL = 1, then the limit is applied to the thermal beta only
    !! <P><LI>If ICULBL = 2, then the limit is applied to the thermal +
    !! neutral beam beta components
    !! </UL>
    !! The default value for the g coefficient is DNBETA = 3.5
    !! AEA FUS 172: Physics Assessment for the European Reactor Study
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    real(dp), intent(in) :: bt, dnbeta, plascur, rminor
    real(dp), intent(out) :: betalim

    !  Local variables

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    betalim = 0.01D0 * dnbeta * (plascur/1.0D6) / (rminor*bt)

  end subroutine culblm

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine plasma_composition

    !! Calculates various plasma component fractional makeups
    !! author: P J Knight, CCFE, Culham Science Centre
    !! None
    !! This subroutine determines the various plasma component
    !! fractional makeups. It is the replacement for the original
    !! It is the replacement for the original routine <CODE>betcom</CODE>,
    !! and is used in conjunction with the new impurity radiation model
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use current_drive_variables, only: ftritbm
		use error_handling, only: fdiags, report_error
    use impurity_radiation_module, only: nimp, impurity_arr_frac, element2index, &
      zav_of_te, impurity_arr_Z, impurity_arr_amass
    use physics_variables, only: alphat, ignite, falpe, afuel, ftrit, deni, &
      aion, dnitot, protium, zeffai, rncne, rnone, falpi, ralpne, dlamee, &
      rnbeam, zeff, dnz, pcoef, alpharate, rnfene, abeam, dlamie, te, &
      protonrate, fdeut, alphan, dnbeam, fhe3, dnalp, dene, dnprot
		use maths_library, only: secant_solve
    implicit none

    !  Arguments

    !  Local variables

    real(dp) :: znimp, pc, znfuel
    integer :: imp

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Ion density components
    !  ======================

    !  Alpha ash portion

    dnalp = dene * ralpne

    !  Protons
    !  This calculation will be wrong on the first call as the particle
    !  production rates are evaluated later in the calling sequence
    !  Issue #557 Allow protium impurity to be specified: 'protium'
    !  This will override the calculated value which is a minimum.

    if (alpharate < 1.0D-6) then  !  not calculated yet...
       dnprot = max(protium*dene, dnalp * (fhe3 + 1.0D-3)) !  rough estimate
    else
       dnprot = max(protium*dene, dnalp * protonrate/alpharate)
    end if

    !  Beam hot ion component
    !  If ignited, prevent beam fusion effects

    if (ignite == 0) then
       dnbeam = dene * rnbeam
    else
       dnbeam = 0.0D0
    end if

    !  Sum of Zi.ni for all impurity ions (those with charge > helium)

    znimp = 0.0D0
    do imp = 1,nimp
       if (impurity_arr_Z(imp) > 2) then
         ! znimp = znimp + impurity_arr(imp)%Z*(impurity_arr_frac(imp) * dene)
          znimp = znimp + Zav_of_te(imp,te)*(impurity_arr_frac(imp) * dene)
       end if
    end do

    !  Fuel portion - conserve charge neutrality
    !  znfuel is the sum of Zi.ni for the three fuel ions

    znfuel = dene - 2.0D0*dnalp - dnprot - dnbeam - znimp

    !  Fuel ion density, deni
    !  For D-T-He3 mix, deni = nD + nT + nHe3, while znfuel = nD + nT + 2*nHe3
    !  So deni = znfuel - nHe3 = znfuel - fhe3*deni

    deni = znfuel/(1.0D0+fhe3)

    !  Ensure that deni is never negative or zero

    if (deni < 0.0D0) then
       fdiags(1) = deni ; call report_error(78)
       deni = max(deni,1.0D0)
    end if

    !  Set hydrogen and helium impurity fractions for
    !  radiation calculations

    impurity_arr_frac(element2index('H_')) = &
         (dnprot + (fdeut+ftrit)*deni + dnbeam)/dene

    impurity_arr_frac(element2index('He')) = fhe3*deni/dene + ralpne

    !  Total impurity density

    dnz = 0.0D0
    do imp = 1,nimp
       if (impurity_arr_Z(imp) > 2) then
          dnz = dnz + impurity_arr_frac(imp)*dene
       end if
    end do

    !  Total ion density

    dnitot = deni + dnalp + dnprot + dnbeam + dnz

    !  Set some (obsolescent) impurity fraction variables
    !  for the benefit of other routines

    rncne = impurity_arr_frac(element2index('C_'))
    rnone = impurity_arr_frac(element2index('O_'))
    ! Issue #261 Remove zfear.  Use the sum of Fe and Ar concentrations
    ! if (zfear == 0) then
    !    rnfene = impurity_arr(element2index('Fe'))%frac
    ! else
    !    rnfene = impurity_arr(element2index('Ar'))%frac
    ! end if
    rnfene = impurity_arr_frac(element2index('Fe')) + impurity_arr_frac(element2index('Ar'))

    !  Effective charge
    !  Calculation should be sum(ni.Zi^2) / sum(ni.Zi),
    !  but ne = sum(ni.Zi) through quasineutrality

    zeff = 0.0D0
    do imp = 1,nimp
       !zeff = zeff + impurity_arr_frac(imp) * (impurity_arr(imp)%Z)**2
       zeff = zeff + impurity_arr_frac(imp) * Zav_of_te(imp,te)**2
    end do

    !  Define coulomb logarithm
    !  (collisions: ion-electron, electron-electron)

    dlamee = 31.0D0 - (log(dene)/2.0D0) + log(te*1000.0D0)
    dlamie = 31.3D0 - (log(dene)/2.0D0) + log(te*1000.0D0)

    !  Fraction of alpha energy to ions and electrons
    !  From Max Fenstermacher
    !  (used with electron and ion power balance equations only)
    !  No consideration of pchargepv here...

    !  pcoef now calculated in plasma_profiles, after the very first
    !  call of plasma_composition; use old parabolic profile estimate
    !  in this case

    if (first_call == 1) then
       pc = (1.0D0 + alphan)*(1.0D0 + alphat)/(1.0D0+alphan+alphat)
       first_call = 0
    else
       pc = pcoef
    end if

    falpe = 0.88155D0 * exp(-te*pc/67.4036D0)
    falpi = 1.0D0 - falpe

    !  Average atomic masses

    afuel = 2.0D0*fdeut + 3.0D0*ftrit + 3.0D0*fhe3
    abeam = 2.0D0*(1.0D0-ftritbm) + 3.0D0*ftritbm

    !  Density weighted mass

    aion = afuel*deni + 4.0D0*dnalp + dnprot + abeam*dnbeam
    do imp = 1,nimp
       if (impurity_arr_Z(imp) > 2) then
          aion = aion + dene*impurity_arr_frac(imp)*impurity_arr_amass(imp)
       end if
    end do
    aion = aion/dnitot

    !  Mass weighted plasma effective charge

    zeffai = ( fdeut*deni/2.0D0 + ftrit*deni/3.0D0 + 4.0D0*fhe3*deni/3.0D0 + &
         dnalp + dnprot + (1.0D0-ftritbm)*dnbeam/2.0D0 + ftritbm*dnbeam/3.0D0 &
         ) / dene
    do imp = 1,nimp
       if (impurity_arr_Z(imp) > 2) then
          zeffai = zeffai + impurity_arr_frac(imp) &
          !     * (impurity_arr(imp)%Z)**2 / impurity_arr_amass(imp)
               * Zav_of_te(imp,te)**2 / impurity_arr_amass(imp)
       end if
    end do

  end subroutine plasma_composition

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine culdlm(bt,idensl,pdivt,plascur,prn1,qcyl,q95, &
       rmajor,rminor,sarea,zeff,dlimit,dnelimt)

    !! Density limit calculation
    !! author: P J Knight, CCFE, Culham Science Centre
    !! bt       : input real :  toroidal field on axis (T)
    !! idensl   : input/output integer : switch denoting which formula to enforce
    !! pdivt    : input real :  power flowing to the edge plasma via
    !! charged particles (MW)
    !! plascur  : input real :  plasma current (A)
    !! prn1     : input real :  edge density / average plasma density
    !! qcyl     : input real :  equivalent cylindrical safety factor (qstar)
    !! q95      : input real :  safety factor at 95% surface
    !! rmajor   : input real :  plasma major radius (m)
    !! rminor   : input real :  plasma minor radius (m)
    !! sarea    : input real :  plasma surface area (m**2)
    !! zeff     : input real :  plasma effective charge
    !! dlimit(7): output real array : average plasma density limit using
    !! seven different models (m**-3)
    !! dnelimt  : output real : enforced average plasma density limit (m**-3)
    !! This routine calculates several different formulae for the
    !! density limit, and enforces the one chosen by the user.
    !! AEA FUS 172: Physics Assessment for the European Reactor Study
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use error_handling, only: fdiags, idiags, report_error
		use constants, only: pi
    implicit none

    !  Arguments

    integer, intent(inout) :: idensl
    real(dp), intent(in) :: bt, pdivt, plascur, prn1, q95, &
         qcyl, rmajor, rminor, sarea, zeff
    real(dp), intent(out) :: dnelimt
    real(dp), dimension(7), intent(out) :: dlimit

    !  Local variables

    real(dp) :: denom, dlim, qperp

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Check for illegal values of IDENSL

    if ((idensl < 1).or.(idensl > 7)) then
       idiags(1) = idensl ; call report_error(79)
    end if

    dlimit(:) = 0.0D0

    !  Power per unit area crossing the plasma edge
    !  (excludes radiation and neutrons)

      qperp = pdivt/sarea

    !  Old ASDEX density limit formula
    !  This applies to the density at the plasma edge, so must be scaled
    !  to give the density limit applying to the average plasma density.

    dlim = 1.54D20 * qperp**0.43D0 * bt**0.31D0 /(q95*rmajor)**0.45D0
    dlimit(1) = dlim/prn1

    !  Borrass density limit model for ITER (I)
    !  This applies to the density at the plasma edge, so must be scaled
    !  to give the density limit applying to the average plasma density.
    !  Borrass et al, ITER-TN-PH-9-6 (1989)

    dlim = 1.8D20 * qperp**0.53D0 * bt**0.31D0 /(q95*rmajor)**0.22D0
    dlimit(2) = dlim/prn1

    !  Borrass density limit model for ITER (II)
    !  This applies to the density at the plasma edge, so must be scaled
    !  to give the density limit applying to the average plasma density.
    !  This formula is (almost) identical to that in the original routine
    !  denlim (now deleted).

    dlim = 0.5D20 * qperp**0.57D0 * bt**0.31D0 /(q95*rmajor)**0.09D0
    dlimit(3) = dlim/prn1

    !  JET edge radiation density limit model
    !  This applies to the density at the plasma edge, so must be scaled
    !  to give the density limit applying to the average plasma density.
    !  qcyl=qstar here, but literature is not clear.

    denom = (zeff-1.0D0)*( 1.0D0-4.0D0/(3.0D0*qcyl) )
    if (denom <= 0.0D0) then
       if (idensl == 4) then
          fdiags(1) = denom ; fdiags(2) = qcyl
          call report_error(80)
          idensl = 5
       end if
       dlimit(4) = 0.0D0
    else
       dlim = 1.0D20 * sqrt(pdivt/denom)
       dlimit(4) = dlim/prn1
    end if

    !  JET simplified density limit model
    !  This applies to the density at the plasma edge, so must be scaled
    !  to give the density limit applying to the average plasma density.

    dlim = 0.237D20 * bt * sqrt(pdivt)/rmajor
    dlimit(5) = dlim/prn1

    !  Hugill-Murakami M.q limit
    !  qcyl=qstar here, which is okay according to the literature

    dlimit(6) = 3.0D20 * bt / (rmajor*qcyl)

    !  Greenwald limit

    dlimit(7) = 1.0D14 * plascur/(pi*rminor*rminor)

    !  Enforce the chosen density limit

    dnelimt = dlimit(idensl)

  end subroutine culdlm

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pcond(afuel,palpmw,aspect,bt,dnitot,dene,dnla,eps,hfact, &
       iinvqd,isc,ignite,kappa,kappa95,kappaa,pchargemw,pinjmw,&
       plascur,pcoreradpv,rmajor,rminor,te,ten,tin,q,qstar,vol, &
       xarea,zeff,ptrepv,ptripv,tauee,tauei,taueff,powerht)

    !! Routine to calculate the confinement times and
    !! the transport power loss terms.
    !! author: P J Knight, CCFE, Culham Science Centre
    !! afuel     : input real :  average mass of fuel (amu)
    !! palpmw    : input real :  alpha particle power (MW)
    !! aspect    : input real :  aspect ratio
    !! bt        : input real :  toroidal field on axis (T)
    !! dene      : input real :  volume averaged electron density (/m3)
    !! dnitot    : input real :  total ion density (/m3)
    !! dnla      : input real :  line-averaged electron density (/m3)
    !! eps       : input real :  inverse aspect ratio
    !! hfact     : input real :  H factor on energy confinement scalings
    !! iinvqd    : input integer :  switch for inverse quadrature
    !! isc       : input integer :  switch for energy confinement scaling to use
    !! ignite    : input integer :  switch for ignited calculation
    !! kappa     : input real :  plasma elongation
    !! kappa95   : input real :  plasma elongation at 95% surface
    !! kappaa    : output real : plasma elongation calculated using area ratio
    !! pchargemw : input real :  non-alpha charged particle fusion power (MW)
    !! pinjmw    : input real :  auxiliary power to ions and electrons (MW)
    !! plascur   : input real :  plasma current (A)
    !! pcoreradpv: input real :  total core radiation power (MW/m3)
    !! q         : input real :  edge safety factor (tokamaks), or
    !! rotational transform iotabar (stellarators)
    !! qstar     : input real :  equivalent cylindrical edge safety factor
    !! rmajor    : input real :  plasma major radius (m)
    !! rminor    : input real :  plasma minor radius (m)
    !! te        : input real :  average electron temperature (keV)
    !! ten       : input real :  density weighted average electron temp. (keV)
    !! tin       : input real :  density weighted average ion temperature (keV)
    !! vol       : input real :  plasma volume (m3)
    !! xarea     : input real :  plasma cross-sectional area (m2)
    !! zeff      : input real :  plasma effective charge
    !! ptrepv    : output real : electron transport power (MW/m3)
    !! ptripv    : output real : ion transport power (MW/m3)
    !! tauee     : output real : electron energy confinement time (s)
    !! taueff    : output real : global energy confinement time (s)
    !! tauei     : output real : ion energy confinement time (s)
    !! powerht   : output real : heating power (MW) assumed in calculation
    !! This subroutine calculates the energy confinement time
    !! using one of a large number of scaling laws, and the
    !! transport power loss terms.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !! N. A. Uckan and ITER Physics Group,
    !! "ITER Physics Design Guidelines: 1989",
    !! ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)
    !! ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)
    !! A. Murari et al 2015 Nucl. Fusion, 55, 073009
    !! C.C. Petty 2008 Phys. Plasmas, 15, 080501
    !! P.T. Lang et al. 2012 IAEA conference proceeding EX/P4-01
    !! ITER physics basis Chapter 2, 1999 Nuclear Fusion 39 2175
    !! Nuclear Fusion corrections, 2008 Nuclear Fusion 48 099801
    !! Menard 2019, Phil. Trans. R. Soc. A 377:20170440
    !! Kaye et al. 2006, Nucl. Fusion 46 848
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use error_handling, only: idiags, report_error
    use physics_variables, only: iradloss, tauee_in, pradpv, kappaa_ipb, &
      pohmmw, falpha
		use startup_variables, only: ptaue, gtaue, ftaue, rtaue, qtaue
		use constants, only: pi
    implicit none

    !  Arguments
    integer, intent(in) :: iinvqd, isc, ignite
    real(dp), intent(in) :: afuel, palpmw, aspect, bt, dene, &
         dnitot, dnla, eps, hfact, kappa, kappa95, pchargemw, pinjmw, &
         plascur, pcoreradpv, q, qstar, rmajor, rminor, te, &
         ten, tin, vol, xarea, zeff
    real(dp), intent(out) :: kappaa, powerht, ptrepv, ptripv, &
         tauee, taueff, tauei

    !  Local variables
    real(dp) :: chii,ck2,denfac,dnla19,dnla20,eps2,gjaeri,iotabar, &
         n20,pcur,qhat,ratio,rll,str2,str5,taueena,tauit1,tauit2, &
         term1,term2, h, qratio, nratio, nGW, taunstx,taupetty

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Neoclassical ion transport loss
    !  Calculate ion energy confinement time
    !
    !  N.B. This calculation is superseded later in the routine

    eps2 = eps/2.0D0
    str5 = (2.0D0/(1.0D0+(kappa**2)))
    ck2 = (0.66D0+(1.88D0*(sqrt(eps2)))-(1.54D0*eps2))* &
         (1.0D0+(1.5D0*(eps2**2)))
    chii = (6.5D-22)*ck2*zeff*(aspect**1.5D0)*dene*(q**2)*str5/ &
         ((sqrt(tin))*(bt**2))
    str2 = (2.0D0*(kappa**2)/(1.0D0+(kappa**2)))
    tauei = 0.375D0*rminor**2/chii*str2

    !  Calculate heating power (MW)
    powerht = falpha*palpmw + pchargemw + pohmmw

    !  If the device is not ignited, add the injected auxiliary power
    if (ignite == 0) powerht = powerht + pinjmw

    !  Include the radiation as a loss term if requested
    if (iradloss == 0) then
       powerht = powerht - pradpv*vol
    else if (iradloss == 1) then
       powerht = powerht - pcoreradpv*vol ! shouldn't this be vol_core instead of vol?
    else
       continue  !  do not adjust powerht for radiation
    end if

    !  Ensure heating power is positive (shouldn't be necessary)
    powerht = max(powerht,1.0D-3)

    !  Line averaged electron density in scaled units
    dnla20 = dnla * 1.0D-20
    dnla19 = dnla * 1.0D-19

    !  Volume averaged electron density in units of 10**20 m**-3
    n20 = dene / 1.0D20

    !  Plasma current in MA
    pcur = plascur / 1.0D6

    ! Separatrix kappa defined with X-section for general use
    kappaa = xarea/(pi*rminor*rminor)

    ! Separatrix kappa defined with plasma volume for IPB scalings
    kappaa_IPB = vol / ( 2.0D0 * pi*pi * rminor*rminor * rmajor )

    !  Calculate Neo-Alcator confinement time (used in several scalings)
    taueena = 0.07D0 * n20 * rminor * rmajor*rmajor * qstar

    !  For reference (see startup.f90):
    !  gtaue = offset term in tauee scaling
    !  ptaue = exponent for density term in tauee scaling
    !  qtaue = exponent for temperature term in tauee scaling
    !  rtaue = exponent for power term in tauee scaling

    !  Electron energy confinement times

    select case (isc)

    case (1)  !  Neo-Alcator scaling (ohmic)
       !tauee = taueena
       tauee = hfact * taueena
       gtaue = 0.0D0
       ptaue = 1.0D0
       qtaue = 0.0D0
       rtaue = 0.0D0

    case (2)  !  Mirnov scaling (H-mode)
       tauee = hfact * 0.2D0 * rminor * sqrt(kappa95) * pcur
       gtaue = 0.0D0
       ptaue = 0.0D0
       qtaue = 0.0D0
       rtaue = 0.0D0

    case (3)  !  Merezhkin-Muhkovatov scaling (L-mode)
       tauee = hfact * 3.5D-3 * rmajor**2.75D0 * rminor**0.25D0 * &
            kappa95**0.125D0 * qstar * dnla20 * sqrt(afuel) / &
            sqrt(ten/10.0D0)
       gtaue = 0.0D0
       ptaue = 1.0D0
       qtaue = -0.5D0
       rtaue = 0.0D0

    case (4)  !  Shimomura scaling (H-mode)
       tauee = hfact * 0.045D0 * rmajor * rminor * bt * sqrt(kappa95) &
            * sqrt(afuel)
       gtaue = 0.0D0
       ptaue = 0.0D0
       qtaue = 0.0D0
       rtaue = 0.0D0

    case (5)  !  Kaye-Goldston scaling (L-mode)
       tauee = hfact * 0.055D0 * kappa95**0.28D0 * pcur**1.24D0 * &
            n20**0.26D0 * rmajor**1.65D0 * sqrt(afuel/1.5D0) / &
            ( bt**0.09D0 * rminor**0.49D0 * powerht**0.58D0 )
       gtaue = 0.0D0
       ptaue = 0.26D0
       qtaue = 0.0D0
       rtaue = -0.58D0
       if (iinvqd /= 0) tauee = 1.0D0 / &
            sqrt(1.0D0/taueena**2 + 1.0D0/tauee**2)

    case (6)  !  ITER Power scaling - ITER 89-P (L-mode)
       tauee = hfact * 0.048D0 * pcur**0.85D0 * rmajor**1.2D0 * &
            rminor**0.3D0 * sqrt(kappa) * dnla20**0.1D0 * bt**0.2D0 * &
            sqrt(afuel) / sqrt(powerht)
       gtaue = 0.0D0
       ptaue = 0.1D0
       qtaue = 0.0D0
       rtaue = -0.5D0

    case (7)  !  ITER Offset linear scaling - ITER 89-O (L-mode)

       term1 = 0.04D0 * pcur**0.5D0 * rmajor**0.3D0 * &
            rminor**0.8D0 * kappa**0.6D0 * afuel**0.5D0
       term2 = 0.064D0 * pcur**0.8D0 * rmajor**1.6D0 * &
            rminor**0.6D0 * kappa**0.5D0 * dnla20**0.6D0 * &
            bt**0.35D0 * afuel**0.2D0 / powerht
       tauee = hfact * (term1 + term2)
       gtaue = hfact*term1
       ptaue = 0.6D0
       qtaue = 0.0D0
       rtaue = -1.0D0

    case (8)  !  Rebut-Lallia offset linear scaling (L-mode)
       rll = (rminor**2 * rmajor * kappa95)**0.333D0
       tauee = hfact * 1.65D0 * sqrt(afuel/2.0D0) * &
            ( 1.2D-2 * pcur * rll**1.5D0 / sqrt(zeff) + &
            0.146D0 * dnla20**0.75D0 * sqrt(pcur) * sqrt(bt) * &
            rll**2.75D0 * zeff**0.25D0 /powerht )
       gtaue = hfact * 1.65D0 * sqrt(afuel/2.0D0) * &
            (1.2D-2 * pcur * rll**1.5D0 / sqrt(zeff))
       ptaue = 0.75D0
       qtaue = 0.0D0
       rtaue = -1.0D0

    case (9)  !  Goldston scaling (L-mode)
       tauee = hfact * 0.037D0 * pcur * rmajor**1.75D0 * &
            rminor**(-0.37D0) * sqrt(kappa95) * sqrt(afuel/1.5D0) / &
            sqrt(powerht)
       gtaue = 0.0D0
       ptaue = 0.0D0
       qtaue = 0.0D0
       rtaue = -0.5D0
       if (iinvqd /= 0) tauee = 1.0D0 / &
            sqrt(1.0D0/taueena**2 + 1.0D0/tauee**2)

    case (10)  !  T10 scaling
       denfac = dnla20 * rmajor * qstar / (1.3D0*bt)
       denfac = min(1.0D0,denfac)
       tauee = hfact * 0.095D0 * rmajor * rminor * bt * &
            sqrt(kappa95) * denfac / powerht**0.4D0 * &
            ( zeff**2 * pcur**4 / &
            (rmajor * rminor * qstar**3 * kappa95**1.5D0) )**0.08D0
       gtaue = 0.0D0
       ptaue = 1.0D0
       qtaue = 0.0D0
       rtaue = -0.4D0

    case (11)  !  JAERI scaling
       gjaeri = zeff**0.4D0 * ((15.0D0-zeff)/20.0D0)**0.6D0 * &
            (3.0D0 * qstar * (qstar+5.0D0) / ((qstar+2.0D0) * &
            (qstar+7.0D0)))**0.6D0
       tauee = hfact * (0.085D0 * kappa95 * rminor**2 * sqrt(afuel) + &
            0.069D0 * n20**0.6D0 * pcur * bt**0.2D0 * rminor**0.4D0 * &
            rmajor**1.6D0 * sqrt(afuel) * gjaeri * kappa95**0.2D0 / &
            powerht)
       gtaue = hfact * 0.085D0 * kappa95 * rminor**2 * sqrt(afuel)
       ptaue = 0.6D0
       qtaue = 0.0D0
       rtaue = -1.0D0

    case (12)  !  Kaye-Big scaling
       tauee = hfact * 0.105D0 * sqrt(rmajor) * rminor**0.8D0 * &
            bt**0.3D0 * kappa95**0.25D0 * pcur**0.85D0 * &
            n20**0.1D0 * afuel**0.5D0 / powerht**0.5D0
       gtaue = 0.0D0
       ptaue = 0.1D0
       qtaue = 0.0D0
       rtaue = -0.5D0

    case (13)  !  ITER H-mode scaling - ITER H90-P
       tauee = hfact * 0.064D0 * pcur**0.87D0 * rmajor**1.82D0 * &
            rminor**(-0.12D0) * kappa**0.35D0 * dnla20**0.09D0 * &
            bt**0.15D0 * sqrt(afuel) / sqrt(powerht)
       gtaue = 0.0D0
       ptaue = 0.09D0
       qtaue = 0.0D0
       rtaue = -0.5D0

    case (14)  !  Minimum of ITER 89-P (isc=6) and ITER 89-O (isc=7)
       tauit1 = hfact * 0.048D0 * pcur**0.85D0 * rmajor**1.2D0 * &
            rminor**0.3D0 * sqrt(kappa) * dnla20**0.1D0 * bt**0.2D0 * &
            sqrt(afuel) / sqrt(powerht)
       term1 = 0.04D0 * pcur**0.5D0 * rmajor**0.3D0 * &
            rminor**0.8D0 * kappa**0.6D0 * afuel**0.5D0
       term2 = 0.064D0 * pcur**0.8D0 * rmajor**1.6D0 * &
            rminor**0.6D0 * kappa**0.5D0 * dnla20**0.6D0 * &
            bt**0.35D0 * afuel**0.2D0 / powerht
       tauit2 = hfact * (term1 + term2)
       tauee = min(tauit1,tauit2)

       if (tauit1 < tauit2) then
          gtaue = 0.0D0
          ptaue = 0.1D0
          qtaue = 0.0D0
          rtaue = -0.5D0
       else
          gtaue = hfact*term1
          ptaue = 0.6D0
          qtaue = 0.0D0
          rtaue = -1.0D0
       end if

    case (15)  !  Riedel scaling (L-mode)
       tauee = hfact * 0.044D0 * pcur**0.93D0 * rmajor**1.37D0 * &
            rminor**(-0.049D0) * kappa95**0.588D0 * dnla20**0.078D0 * &
            bt**0.152D0 / powerht**0.537D0
       gtaue = 0.0D0
       ptaue = 0.078D0
       qtaue = 0.0D0
       rtaue = -0.537D0

    case (16)  !  Christiansen et al scaling (L-mode)
       tauee = hfact * 0.24D0 * pcur**0.79D0 * rmajor**0.56D0 * &
            rminor**1.46D0 * kappa95**0.73D0 * dnla20**0.41D0 * &
            bt**0.29D0 / (powerht**0.79D0 * afuel**0.02D0)
       gtaue = 0.0D0
       ptaue = 0.41D0
       qtaue = 0.0D0
       rtaue = -0.79D0

    case (17)  !  Lackner-Gottardi scaling (L-mode)
       qhat = (1.0D0+kappa95**2) * rminor**2 * bt /(0.4D0 * pcur * rmajor)
       tauee = hfact * 0.12D0 * pcur**0.8D0 * rmajor**1.8D0 * &
            rminor**0.4D0 * kappa95 * (1.0D0+kappa95)**(-0.8D0) * &
            dnla20**0.6D0 * qhat**0.4D0 / powerht**0.6D0
       gtaue = 0.0D0
       ptaue = 0.6D0
       qtaue = 0.0D0
       rtaue = -0.6D0

    case (18)  !  Neo-Kaye scaling (L-mode)
       tauee = hfact * 0.063D0 * pcur**1.12D0 * rmajor**1.3D0 * &
            rminor**(-0.04D0) * kappa95**0.28D0 * dnla20**0.14D0 * &
            bt**0.04D0 * sqrt(afuel) / powerht**0.59D0
       gtaue = 0.0D0
       ptaue = 0.14D0
       qtaue = 0.0D0
       rtaue = -0.59D0

    case (19)  !  Riedel scaling (H-mode)
       tauee = hfact * 0.1D0 * sqrt(afuel) * pcur**0.884D0 * &
            rmajor**1.24D0 * rminor**(-0.23D0) * kappa95**0.317D0 * &
            bt**0.207D0 * dnla20**0.105D0 / powerht**0.486D0
       gtaue = 0.0D0
       ptaue = 0.105D0
       qtaue = 0.0D0
       rtaue = -0.486D0

    case (20)  !  Amended version of ITER H90-P law
       !  Nuclear Fusion 32 (1992) 318
       tauee = hfact * 0.082D0 * pcur**1.02D0 * &
            bt**0.15D0 * sqrt(afuel) * rmajor**1.60D0 / &
            (powerht**0.47D0 * kappa**0.19D0)
       gtaue = 0.0D0
       ptaue = 0.0D0
       qtaue = 0.0D0
       rtaue = -0.47D0

    case (21)  !  Large Helical Device scaling (stellarators)
       !  S.Sudo, Y.Takeiri, H.Zushi et al., Nuclear Fusion 30 (1990) 11
       tauee = hfact * 0.17D0 * rmajor**0.75D0 * rminor**2 * &
            dnla20**0.69D0 * bt**0.84D0 * powerht**(-0.58D0)
       gtaue = 0.0D0
       ptaue = 0.69D0
       qtaue = 0.0D0
       rtaue = 0.58D0

    case (22)  !  Gyro-reduced Bohm scaling
       !  R.J.Goldston, H.Biglari, G.W.Hammett et al., Bull.Am.Phys.Society,
       !  volume 34, 1964 (1989)
       tauee = hfact * 0.25D0 * bt**0.8D0 * dnla20**0.6D0 * &
            powerht**(-0.6D0) * rminor**2.4D0 * rmajor**0.6D0
       gtaue = 0.0D0
       ptaue = 0.6D0
       qtaue = 0.0D0
       rtaue = -0.6D0

    case (23)  !  Lackner-Gottardi stellarator scaling
       !  K.Lackner and N.A.O.Gottardi, Nuclear Fusion, 30, p.767 (1990)
       iotabar = q  !  dummy argument q is actual argument iotabar for stellarators
       tauee = hfact * 0.17D0 * rmajor * rminor**2 * dnla20**0.6D0 * &
            bt**0.8D0 * powerht**(-0.6D0) * iotabar**0.4D0
       gtaue = 0.0D0
       ptaue = 0.6D0
       qtaue = 0.0D0
       rtaue = -0.6D0

    case (24)  !  ITER-93H scaling (ELM-free; multiply by 0.85 for ELMy version)
       !  S.Kaye and the ITER Joint Central Team and Home Teams, in Plasma
       !  Physics and Controlled Nuclear Fusion Research (Proc. 15th
       !  Int. Conf., Seville, 1994) IAEA-CN-60/E-P-3
       tauee = hfact * 0.053D0 * pcur**1.06D0 * bt**0.32D0 * &
            powerht**(-0.67D0) * afuel**0.41D0 * rmajor**1.79D0 * &
            dnla20**0.17D0 * aspect**0.11D0 * kappa**0.66D0
       gtaue = 0.0D0
       ptaue = 0.17D0
       qtaue = 0.0D0
       rtaue = -0.67D0

   case (25)  !  Issue #508 Remove RFP option.

       !  Next two are ITER-97 H-mode scalings
       !  J. G. Cordey et al., EPS Berchtesgaden, 1997

    case (26)  !  ELM-free: ITERH-97P
       tauee = hfact * 0.031D0 * pcur**0.95D0 * bt**0.25D0 * &
            powerht**(-0.67D0) * dnla19**0.35D0 * &
            rmajor**1.92D0 * aspect**(-0.08D0) * kappa**0.63D0 * &
            afuel**0.42D0
       gtaue = 0.0D0
       ptaue = 0.35D0
       qtaue = 0.0D0
       rtaue = -0.67D0

    case (27)  !  ELMy: ITERH-97P(y)
       tauee = hfact * 0.029D0 * pcur**0.90D0 * bt**0.20D0 * &
            powerht**(-0.66D0) * dnla19**0.40D0 * &
            rmajor**2.03D0 * aspect**(-0.19D0) * kappa**0.92D0 * &
            afuel**0.2D0
       gtaue = 0.0D0
       ptaue = 0.4D0
       qtaue = 0.0D0
       rtaue = -0.66D0

    case (28)  !  ITER-96P (= ITER-97L) L-mode scaling
       !  S.M.Kaye and the ITER Confinement Database Working Group,
       !  Nuclear Fusion 37 (1997) 1303
       !  N.B. tau_th formula used
       tauee = hfact * 0.023D0 * pcur**0.96D0 * bt**0.03D0 * &
            kappa95**0.64D0 * rmajor**1.83D0 * aspect**0.06D0 * &
            dnla19**0.40D0 * afuel**0.20D0 * powerht**(-0.73D0)
       gtaue = 0.0D0
       ptaue = 0.4D0
       qtaue = 0.0D0
       rtaue = -0.73D0

    case (29)  !  Valovic modified ELMy-H mode scaling
       tauee = hfact * 0.067D0 * pcur**0.9D0 * bt**0.17D0 * &
            dnla19**0.45D0 * afuel**0.05D0 * rmajor**1.316D0 * &
            rminor**0.79D0 * kappa**0.56D0 * powerht**(-0.68D0)
       gtaue = 0.0D0
       ptaue = 0.45D0
       qtaue = 0.0D0
       rtaue = -0.68D0

    case (30)  !  Kaye PPPL Workshop April 1998 L-mode scaling
       tauee = hfact * 0.021D0 * pcur**0.81D0 * bt**0.14D0 * &
            kappa**0.7D0 * rmajor**2.01D0 * aspect**(-0.18D0) * &
            dnla19**0.47D0 * afuel**0.25D0 * powerht**(-0.73D0)
       gtaue = 0.0D0
       ptaue = 0.47D0
       qtaue = 0.0D0
       rtaue = -0.73D0

    case (31)  !  ITERH-PB98P(y), ELMy H-mode scaling
       tauee = hfact * 0.0615D0 * pcur**0.9D0 * bt**0.1D0 * &
            dnla19**0.4D0 * powerht**(-0.66D0) * rmajor**2 * &
            kappaa**0.75D0 * aspect**(-0.66D0) * afuel**0.2D0
       gtaue = 0.0D0
       ptaue = 0.4D0
       qtaue = 0.0D0
       rtaue = -0.66D0

    case (32)  !  IPB98(y), ELMy H-mode scaling
       !  Data selection : full ITERH.DB3
       !  Nuclear Fusion 39 (1999) 2175, Table 5
       tauee = hfact * 0.0365D0 * pcur**0.97D0 * bt**0.08D0 * &
            dnla19**0.41D0 * powerht**(-0.63D0) * rmajor**1.93D0 * &
            kappa**0.67D0 * aspect**(-0.23D0) * afuel**0.2D0
       gtaue = 0.0D0
       ptaue = 0.41D0
       qtaue = 0.0D0
       rtaue = -0.63D0

    case (33)  !  IPB98(y,1), ELMy H-mode scaling
       !  Data selection : full ITERH.DB3
       !  Nuclear Fusion 39 (1999) 2175, Table 5
       tauee = hfact * 0.0503D0 * pcur**0.91D0 * bt**0.15D0 * &
            dnla19**0.44D0 * powerht**(-0.65D0) * rmajor**2.05D0 * &
            kappaa_IPB**0.72D0 * aspect**(-0.57D0) * afuel**0.13D0
       gtaue = 0.0D0
       ptaue = 0.44D0
       qtaue = 0.0D0
       rtaue = -0.65D0

    case (34)  !  IPB98(y,2), ELMy H-mode scaling
       !  Data selection : ITERH.DB3, NBI only
       !  Nuclear Fusion 39 (1999) 2175, Table 5
       tauee = hfact * 0.0562D0 * pcur**0.93D0 * bt**0.15D0 * &
            dnla19**0.41D0 * powerht**(-0.69D0) * rmajor**1.97D0 * &
            kappaa_IPB**0.78D0 * aspect**(-0.58D0) * afuel**0.19D0
       gtaue = 0.0D0
       ptaue = 0.41D0
       qtaue = 0.0D0
       rtaue = -0.69D0

    case (35)  !  IPB98(y,3), ELMy H-mode scaling
       !  Data selection : ITERH.DB3, NBI only, no C-Mod
       !  Nuclear Fusion 39 (1999) 2175, Table 5
       tauee = hfact * 0.0564D0 * pcur**0.88D0 * bt**0.07D0 * &
            dnla19**0.40D0 * powerht**(-0.69D0) * rmajor**2.15D0 * &
            kappaa_IPB**0.78D0 * aspect**(-0.64D0) * afuel**0.20D0
       gtaue = 0.0D0
       ptaue = 0.4D0
       qtaue = 0.0D0
       rtaue = -0.69D0

    case (36)  !  IPB98(y,4), ELMy H-mode scaling
       !  Data selection : ITERH.DB3, NBI only, ITER like devices
       !  Nuclear Fusion 39 (1999) 2175, Table 5
       tauee = hfact * 0.0587D0 * pcur**0.85D0 * bt**0.29D0 * &
            dnla19**0.39D0 * powerht**(-0.70D0) * rmajor**2.08D0 * &
            kappaa_IPB**0.76D0 * aspect**(-0.69D0) * afuel**0.17D0
       gtaue = 0.0D0
       ptaue = 0.39D0
       qtaue = 0.0D0
       rtaue = -0.70D0

    case (37)  !  ISS95 stellarator scaling
       !  U. Stroth et al., Nuclear Fusion, 36, p.1063 (1996)
       !  Assumes kappa = 1.0, triang = 0.0
       iotabar = q  !  dummy argument q is actual argument iotabar for stellarators
       tauee = hfact * 0.079D0 * rminor**2.21D0 * rmajor**0.65D0 * dnla19**0.51D0 * &
            bt**0.83D0 * powerht**(-0.59D0) * iotabar**0.4D0
       gtaue = 0.0D0
       ptaue = 0.51D0
       qtaue = 0.0D0
       rtaue = -0.59D0

    case (38)  !  ISS04 stellarator scaling
       !  H. Yamada et al., Nuclear Fusion, 45, p.1684 (2005)
       !  Assumes kappa = 1.0, triang = 0.0
       iotabar = q  !  dummy argument q is actual argument iotabar for stellarators
       tauee = hfact * 0.134D0 * rminor**2.28D0 * rmajor**0.64D0 * dnla19**0.54D0 * &
            bt**0.84D0 * powerht**(-0.61D0) * iotabar**0.41D0
       gtaue = 0.0D0
       ptaue = 0.54D0
       qtaue = 0.0D0
       rtaue = -0.61D0

    case (39)  !  DS03 beta-independent H-mode scaling
       !  T. C. Luce, C. C. Petty and J. G. Cordey,
       !  Plasma Phys. Control. Fusion 50 (2008) 043001, eqn.4.13, p.67
       tauee = hfact * 0.028D0 * pcur**0.83D0 * bt**0.07D0 * &
            dnla19**0.49D0 * powerht**(-0.55D0) * rmajor**2.11D0 * &
            kappa95**0.75D0 * aspect**(-0.3D0) * afuel**0.14D0
       gtaue = 0.0D0
       ptaue = 0.49D0
       qtaue = 0.0D0
       rtaue = -0.55D0

    case (40)  !  "Non-power law" (NPL) Murari energy confinement scaling
       !   Based on the ITPA database of H-mode discharges
       !   A new approach to the formulation and validation of scaling expressions for plasma confinement in tokamaks
       !   A. Murari et al 2015 Nucl. Fusion 55 073009, doi:10.1088/0029-5515/55/7/073009
       !   Table 4.  (Issue #311)
       !  Note that aspect ratio and M (afuel) do not appear, and B (bt) only
       !  appears in the "saturation factor" h.
       h = dnla19**0.448D0 / (1.0D0 + exp(-9.403D0*(bt/dnla19)**1.365D0))
       tauee = hfact * 0.0367D0 * pcur**1.006D0 * rmajor**1.731D0 * kappaa**1.450D0 * &
               powerht**(-0.735D0) * h

       gtaue = 0.0D0
       ptaue = 0.448D0
       qtaue = 0.0D0
       rtaue = -0.735D0

    case (41) ! Beta independent dimensionless confinement scaling
       ! C.C. Petty 2008 Phys. Plasmas 15, 080501, equation 36
       ! Note that there is no dependence on the average fuel mass 'afuel'
       tauee = hfact * 0.052D0 * pcur**0.75D0 * bt**0.3D0 * &
            dnla19**0.32D0 * powerht**(-0.47D0) * rmajor**2.09D0 * &
            kappaa**0.88D0 * aspect**(-0.84D0)

       gtaue = 0.0D0
       ptaue = 0.32D0
       qtaue = 0.0D0
       rtaue = -0.47D0

    case (42) ! High density relevant confinement scaling
       ! P.T. Lang et al. 2012, IAEA conference proceeding EX/P4-01
       ! q should be q95: incorrect if icurr = 2 (ST current scaling)
       qratio = q/qstar
       ! Greenwald density in m^-3
       nGW = 1.0D14 * plascur/(pi*rminor*rminor)
       nratio = dnla/nGW
       tauee = hfact * 6.94D-7 * plascur**1.3678D0 * bt**0.12D0 * &
            dnla**0.032236D0 * (powerht*1.0D6)**(-0.74D0) * rmajor**1.2345D0 * &
            kappaa_IPB**0.37D0 * aspect**2.48205D0 * afuel**0.2D0 * &
            qratio**0.77D0 * aspect**(-0.9D0*log(aspect)) * &
            nratio**(-0.22D0*log(nratio))

       gtaue = 0.0D0
       ptaue = 0.032236D0 -0.22D0*log(nratio)
       qtaue = 0.0D0
       rtaue = -0.74D0

    case (43)  !  Hubbard et al. 2017 I-mode confinement time scaling - nominal
      tauee = hfact * 0.014D0 * (plascur/1.0D6)**0.68D0 * bt**0.77D0 * dnla20**0.02D0 &
              * powerht**(-0.29D0)
      gtaue = 0.0D0
      ptaue = 0.02D0
      qtaue = 0.0D0
      rtaue = -0.29D0

    case (44)  !  Hubbard et al. 2017 I-mode confinement time scaling - lower
      tauee = hfact * 0.014D0 * (plascur/1.0D6)**0.60D0 * bt**0.70D0 * dnla20**(-0.03D0) &
              * powerht**(-0.33D0)
      gtaue = 0.0D0
      ptaue = -0.03D0
      qtaue = 0.0D0
      rtaue = -0.33D0

    case (45)  !  Hubbard et al. 2017 I-mode confinement time scaling - upper
      tauee = hfact * 0.014D0 * (plascur/1.0D6)**0.76D0 * bt**0.84D0  * dnla20**0.07 &
              * powerht**(-0.25D0)
      gtaue = 0.0D0
      ptaue = 0.07D0
      qtaue = 0.0D0
      rtaue = -0.25D0

    case (46)  !  NSTX, ELMy H-mode scaling
      !  NSTX scaling with IPB98(y,2) for other variables
      !  Menard 2019, Phil. Trans. R. Soc. A 377:20170440
      !  Kaye et al. 2006, Nucl. Fusion 46 848
      tauee = hfact * 0.095D0 * pcur**0.57D0 * bt**1.08D0 * &
           dnla19**0.44D0 * powerht**(-0.73D0) * rmajor**1.97D0 * &
           kappaa_IPB**0.78D0 * aspect**(-0.58D0) * afuel**0.19D0
      gtaue = 0.0D0
      ptaue = 0.44D0
      qtaue = 0.0D0
      rtaue = -0.73D0

    case (47) ! NSTX-Petty08 Hybrid
      ! Linear interpolation between NSTX and Petty08 in eps
      ! Menard 2019, Phil. Trans. R. Soc. A 377:20170440
      if ((1.0D0/aspect).le.0.4D0) then
      ! Petty08, i.e. case (41)
        tauee = hfact * 0.052D0 * pcur**0.75D0 * bt**0.3D0 * &
              dnla19**0.32D0 * powerht**(-0.47D0) * rmajor**2.09D0 * &
              kappaa**0.88D0 * aspect**(-0.84D0)

        gtaue = 0.0D0
        ptaue = 0.32D0
        qtaue = 0.0D0
        rtaue = -0.47D0

      else if ((1.0D0/aspect).ge.0.6D0) then
        ! NSTX, i.e.case (46)
        tauee = hfact * 0.095D0 * pcur**0.57D0 * bt**1.08D0 * &
              dnla19**0.44D0 * powerht**(-0.73D0) * rmajor**1.97D0 * &
              kappaa_IPB**0.78D0 * aspect**(-0.58D0) * afuel**0.19D0

        gtaue = 0.0D0
        ptaue = 0.44D0
        qtaue = 0.0D0
        rtaue = -0.73D0

      else
        taupetty = 0.052D0 * pcur**0.75D0 * bt**0.3D0 * &
                dnla19**0.32D0 * powerht**(-0.47D0) * rmajor**2.09D0 * &
                kappaa**0.88D0 * aspect**(-0.84D0)
        taunstx= 0.095D0 * pcur**0.57D0 * bt**1.08D0 * &
                dnla19**0.44D0 * powerht**(-0.73D0) * rmajor**1.97D0 * &
                kappaa_IPB**0.78D0 * aspect**(-0.58D0) * afuel**0.19D0

        tauee = hfact*((((1.0D0/aspect)-0.4D0)/(0.6D0-0.4D0))*taunstx + &
                 ((0.6D0-(1.0D0/aspect))/(0.6D0-0.4D0))*taupetty)

        gtaue = 0.0D0
        ptaue = ((((1.0D0/aspect)-0.4D0)/(0.6D0-0.4D0))*0.32D0 + &
                ((0.6D0-(1.0D0/aspect))/(0.6D0-0.4D0))*0.44D0)
        qtaue = 0.0D0
        rtaue = ((((1.0D0/aspect)-0.4D0)/(0.6D0-0.4D0))*(-0.47D0) + &
                ((0.6D0-(1.0D0/aspect))/(0.6D0-0.4D0))*(-0.73D0))
      end if

    case (48) ! NSTX gyro-Bohm (Buxton)
      ! P F Buxton et al. 2019 Plasma Phys. Control. Fusion 61 035006
      tauee = hfact * 0.21D0 * pcur**0.54D0 * bt**0.91D0 * &
         powerht**(-0.38D0) * rmajor**2.14D0 * dnla20**(-0.05D0)

      gtaue = 0.0D0
      ptaue = -0.05D0
      qtaue = 0.0D0
      rtaue = -0.38D0

    case (49) ! tauee is an input
      tauee = hfact * tauee_in

      gtaue = 0.0D0
      ptaue = 0.0D0
      qtaue = 0.0D0
      rtaue = 0.0D0

    case default
       idiags(1) = isc ; call report_error(81)

    end select

    !  Ion energy confinement time
    !  N.B. Overwrites earlier calculation above

    tauei = tauee

    !  Calculation of the transport power loss terms
    !  Transport losses in Watts/m3 are 3/2 * n.e.T / tau , with T in eV
    !  (here, tin and ten are in keV, and ptrepv and ptripv are in MW/m3)

    ptripv = 2.403D-22 * dnitot*tin/tauei
    ptrepv = 2.403D-22 * dene*ten/tauee

    ratio = dnitot/dene * tin/ten

    !  Global energy confinement time

    taueff = ((ratio + 1.0D0)/(ratio/tauei + 1.0D0/tauee))

    ! This is used only in subroutine startup, which is currently (r400)
    ! not used.
    ftaue = (tauee-gtaue) / &
         (n20**ptaue * (te/10.0D0)**qtaue * powerht**rtaue)

  end subroutine pcond

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine vscalc(csawth,eps,facoh,gamma,kappa,rmajor,rplas, &
       plascur,theat,tburn,phiint,rli,rlp,vsbrn,vsind,vsres,vsstt)

    !! Volt-second requirements
    !! author: P J Knight, CCFE, Culham Science Centre
    !! csawth : input real :  coefficient for sawteeth effects
    !! eps    : input real :  inverse aspect ratio
    !! facoh  : input real :  fraction of plasma current produced inductively
    !! gamma  : input real :  Ejima coeff for resistive start-up V-s component
    !! kappa  : input real :  plasma elongation
    !! plascur: input real :  plasma current (A)
    !! rli    : input real :  plasma normalised inductivity
    !! rmajor : input real :  plasma major radius (m)
    !! rplas  : input real :  plasma resistance (ohm)
    !! theat  : input real :  heating time (s)
    !! tburn  : input real :  burn time (s)
    !! phiint : output real : internal plasma volt-seconds (Wb)
    !! rlp    : output real : plasma inductance (H)
    !! vsbrn  : output real : volt-seconds needed during flat-top (heat+burn) (Wb)
    !! vsind  : output real : internal and external plasma inductance V-s (Wb)
    !! vsres  : output real : resistive losses in start-up volt-seconds (Wb)
    !! vsstt  : output real : total volt-seconds needed (Wb)
    !! This subroutine calculates the volt-second requirements and some
    !! other related items.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use constants, only: rmu0
    implicit none

    !  Arguments

    real(dp), intent(in) :: csawth, eps, facoh, gamma, kappa, &
         plascur, rli, rmajor, rplas, tburn, theat
    real(dp), intent(out) :: phiint, rlp, vsbrn, vsind, vsres, vsstt

    !  Local variables

    real(dp) :: aeps,beps,rlpext,rlpint,vburn

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Internal inductance

    rlpint = rmu0 * rmajor * rli/2.0D0
    phiint = rlpint*plascur

    !  Start-up resistive component
    !  Uses ITER formula without the 10 V-s add-on

    vsres = gamma * rmu0*plascur*rmajor

    !  Hirshman, Neilson: Physics of Fluids, 29 (1986) p790
    !  fit for external inductance

    aeps = (1.0D0 + 1.81D0*sqrt(eps)+2.05D0*eps)*log(8.0D0/eps) &
         - (2.0D0 + 9.25D0*sqrt(eps)-1.21D0*eps)
    beps = 0.73D0 * sqrt(eps) *(1.0D0 + 2.0D0*eps**4-6.0D0*eps**5 &
         + 3.7D0*eps**6)
    rlpext = rmajor*rmu0 * aeps*(1.0D0-eps)/(1.0D0-eps+beps*kappa)

    rlp = rlpext + rlpint

    !  Inductive V-s component

    vsind = rlp * plascur
    vsstt = vsres + vsind

    !  Loop voltage during flat-top
    !  Include enhancement factor in flattop V-s requirement
    !  to account for MHD sawtooth effects.

    vburn = plascur * rplas * facoh * csawth

    !  N.B. tburn on first iteration will not be correct
    !  if the pulsed reactor option is used, but the value
    !  will be correct on subsequent calls.

    vsbrn = vburn*(theat + tburn)
    vsstt = vsstt + vsbrn

  end subroutine vscalc

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine phyaux(aspect,dene,deni,fusionrate,alpharate,plascur,sbar,dnalp, &
       taueff,vol,burnup,dntau,figmer,fusrat,qfuel,rndfuel,taup)

    !! Auxiliary physics quantities
    !! author: P J Knight, CCFE, Culham Science Centre
    !! aspect : input real :  plasma aspect ratio
    !! dene   : input real :  electron density (/m3)
    !! deni   : input real :  fuel ion density (/m3)
    !! dnalp  : input real :  alpha ash density (/m3)
    !! fusionrate : input real :  fusion reaction rate (/m3/s)
    !! alpharate  : input real :  alpha particle production rate (/m3/s)
    !! plascur: input real :  plasma current (A)
    !! sbar   : input real :  exponent for aspect ratio (normally 1)
    !! taueff : input real :  global energy confinement time (s)
    !! vol    : input real :  plasma volume (m3)
    !! burnup : output real : fractional plasma burnup
    !! dntau  : output real : plasma average n-tau (s/m3)
    !! figmer : output real : physics figure of merit
    !! fusrat : output real : number of fusion reactions per second
    !! qfuel  : output real : fuelling rate for D-T (nucleus-pairs/sec)
    !! rndfuel: output real : fuel burnup rate (reactions/s)
    !! taup   : output real : (alpha) particle confinement time (s)
    !! This subroutine calculates extra physics related items
    !! needed by other parts of the code
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use physics_variables, only: tauratio,burnup_in
    implicit none

    !  Arguments

    real(dp), intent(in) :: aspect, dene, deni, dnalp, &
         fusionrate, alpharate, plascur, sbar, taueff, vol
    real(dp), intent(out) :: burnup, dntau, figmer, fusrat, &
         qfuel, rndfuel, taup

    !  Local variables

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    figmer = 1.0D-6 * plascur * aspect**sbar

    dntau = taueff*dene

    !  Fusion reactions per second

    fusrat = fusionrate*vol

    !  Alpha particle confinement time (s)
    !  Number of alphas / alpha production rate

    if (alpharate /= 0.0D0) then
      taup = dnalp / alpharate
    else  !  only likely if DD is only active fusion reaction
      taup = 0.0D0
    end if

    !  Fractional burnup

    !  (Consider detailed model in: G. L. Jackson, V. S. Chan, R. D. Stambaugh,
    !  Fusion Science and Technology, vol.64, no.1, July 2013, pp.8-12)

    !  The ratio of ash to fuel particle confinement times is given by
    !  tauratio
    !  Possible logic...
    !  burnup = fuel ion-pairs burned/m3 / initial fuel ion-pairs/m3;
    !  fuel ion-pairs burned/m3 = alpha particles/m3 (for both D-T and D-He3 reactions)
    !  initial fuel ion-pairs/m3 = burnt fuel ion-pairs/m3 + unburnt fuel-ion pairs/m3
    !  Remember that unburnt fuel-ion pairs/m3 = 0.5 * unburnt fuel-ions/m3
    if (burnup_in <= 1.0D-9) then
    burnup = dnalp / (dnalp + 0.5D0*deni) / tauratio
    else
       burnup = burnup_in
    end if
    !  Fuel burnup rate (reactions/second) (previously Amps)

    rndfuel = fusrat

    !  Required fuelling rate (fuel ion pairs/second) (previously Amps)

    qfuel = rndfuel/burnup

  end subroutine phyaux

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine rether(alphan,alphat,dene,dlamie,te,ti,zeffai,piepv)

    !! Routine to find the equilibration power between the
    !! ions and electrons
    !! author: P J Knight, CCFE, Culham Science Centre
    !! alphan : input real :  density profile index
    !! alphat : input real :  temperature profile index
    !! dene   : input real :  electron density (/m3)
    !! dlamie : input real :  ion-electron coulomb logarithm
    !! te     : input real :  electron temperature (keV)
    !! ti     : input real :  ion temperature (keV)
    !! zeffai : input real :  mass weighted plasma effective charge
    !! piepv  : output real : ion/electron equilibration power (MW/m3)
    !! This routine calculates the equilibration power between the
    !! ions and electrons.
    !! Unknown origin
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    real(dp), intent(in) :: alphan, alphat, dene, dlamie, &
         te, ti, zeffai
    real(dp), intent(out) :: piepv

    !  Local variables

    real(dp) :: conie, profie

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    profie = (1.0D0+alphan)**2 / &
         ( (2.0D0*alphan - 0.5D0*alphat + 1.0D0) * sqrt(1.0D0+alphat) )

    conie = 2.42165D-41 * dlamie * dene**2 * zeffai * profie

    piepv = conie*(ti-te)/(te**1.5D0)

  end subroutine rether

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pohm(facoh,kappa95,plascur,rmajor,rminor,ten,vol, &
       zeff,pohmpv,pohmmw,rpfac,rplas)

    !! Ohmic power calculation
    !! author: P J Knight, CCFE, Culham Science Centre
    !! facoh  : input real :  fraction of plasma current produced inductively
    !! kappa95: input real :  plasma elongation at 95% flux
    !! plascur: input real :  plasma current (A)
    !! rmajor : input real :  plasma major radius (m)
    !! rminor : input real :  plasma minor radius (m)
    !! ten    : input real :  density weighted average electron temperature (keV)
    !! vol    : input real :  plasma volume (m3)
    !! zeff   : input real :  plasma effective charge
    !! pohmpv : output real : ohmic heating power per unit volume (MW/m3)
    !! pohmmw : output real : ohmic heating power (MW)
    !! rpfac  : output real : neoclassical resistivity enhancement factor
    !! rplas  : output real : plasma resistance (ohm)
    !! This routine finds the ohmic heating power per unit volume.
    !! The expression is a good fit for alphan = 0.5, alphat = 1.0,
    !! alphaj = 1.5, aspect = 2.5 -- 4.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
    !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use error_handling, only: fdiags, report_error
		use physics_variables, only: aspect, plasma_res_factor
    implicit none

    !  Arguments

    real(dp), intent(in) :: facoh, kappa95, plascur, rmajor, &
         rminor, ten, vol, zeff
    real(dp), intent(out) :: pohmpv, pohmmw, rpfac, rplas

    !  Local variables

    real(dp) :: t10

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Density weighted electron temperature in 10 keV units

    t10 = ten/10.0D0

    !  Plasma resistance, from loop voltage calculation in IPDG89

    rplas = plasma_res_factor * 2.15D-9 * zeff*rmajor / (kappa95*rminor**2 * t10**1.5D0)

    !  Neo-classical resistivity enhancement factor
    !  Taken from  N. A. Uckan et al, Fusion Technology 13 (1988) p.411.
    !  The expression is valid for aspect ratios in the range 2.5--4.

    rpfac = 4.3D0 - 0.6D0*rmajor/rminor
    rplas = rplas * rpfac

    !  Check to see if plasma resistance is negative
    !  (possible if aspect ratio is too high)

    if (rplas <= 0.0D0) then
       fdiags(1) = rplas ; fdiags(2) = aspect
       call report_error(83)
    end if

    !  Ohmic heating power per unit volume
    !  Corrected from: pohmpv = (facoh*plascur)**2 * ...

    pohmpv = facoh * plascur**2 * rplas * 1.0D-6/vol

    !  Total ohmic heating power

    pohmmw = pohmpv*vol

  end subroutine pohm

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine igmarcal(outfile)

    !! Routine to calculate ignition margin
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile   : input integer : Fortran output unit identifier
    !! This routine calculates the ignition margin at the final point
    !! with different scalings.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use current_drive_variables, only: pinjmw
    use physics_variables, only: vol, palpmw, zeff, pchargemw, hfac, xarea, &
      tin, tauscl, eps, kappaa, dnla, kappa95, ten, te, kappa, dnitot, dene, &
      iinvqd, rminor, bt, rmajor, ignite, aspect, qstar, q, afuel, plascur, &
      pcoreradpv
		use process_output, only: oheadr, oblnkl
    implicit none

    !  Arguments

    integer, intent(in) :: outfile

    !  Local variables

    integer :: iisc
    real(dp), parameter :: d1 = 1.0D0
    real(dp) :: powerhtz, ptrez, ptriz, &
         taueez, taueffz, taueiz

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call oheadr(outfile,'Energy confinement times, and required H-factors :')

    write(outfile,10)
10  format(t5,'scaling law', t30,'confinement time (s)', &
         t55,'H-factor for')

    write(outfile,20)
20  format(t34,'for H = 1',t54,'power balance')

    call oblnkl(outfile)

    !  Calculate power balances for all scaling laws assuming H = 1

    do iisc = 32,47
       call pcond(afuel,palpmw,aspect,bt,dnitot,dene,dnla,eps,d1, &
            iinvqd,iisc,ignite,kappa,kappa95,kappaa,pchargemw,pinjmw, &
            plascur,pcoreradpv,rmajor,rminor,te,ten,tin,q,qstar,vol, &
            xarea,zeff,ptrez,ptriz,taueez,taueiz,taueffz,powerhtz)
       hfac(iisc) = fhfac(iisc)

       write(outfile,30) tauscl(iisc),taueez,hfac(iisc)
    end do
30  format(t2,a24,t34,f7.3,t58,f7.3)

  end subroutine igmarcal

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function fhfac(is)

    !! Function to find H-factor for power balance
    !! author: P J Knight, CCFE, Culham Science Centre
    !! is : input integer : confinement time scaling law of interest
    !! This function calculates the H-factor required for power balance,
    !! using the given energy confinement scaling law.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use maths_library, only: zeroin
    implicit none

    real(dp) :: fhfac

    !  Arguments

    integer, intent(in) :: is

    !  Local variables

    real(dp), parameter :: abserr = 0.003D0  !  numerical tolerance
    real(dp), parameter :: xlow = 0.01D0     !  minimum bound on H-factor
    real(dp), parameter :: xhigh = 100.0D0   !  maximum bound on H-factor

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    iscz = is

    !  Find value of H-factor for which function FHZ is zero
    !  (this occurs at power balance)

    fhfac = zeroin(xlow,xhigh,fhz,abserr)

  end function fhfac

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function fhz(hhh)

    !! Function used to find power balance
    !! author: P J Knight, CCFE, Culham Science Centre
    !! hhh : input real : test value for confinement time H-factor
    !! This function is used to find power balance.
    !! <CODE>FHZ</CODE> is zero at power balance, which is achieved
    !! using routine <A HREF="zeroin.html">ZEROIN</A> to adjust the
    !! value of <CODE>hhh</CODE>, the confinement time H-factor.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use current_drive_variables, only: pinjmw
    use physics_variables, only: iradloss, vol, palpmw, pradpv, pchargemw, &
      zeff, pohmpv, pchargepv, xarea, tin, eps, kappaa, dnla, palppv, kappa95, &
      ten, te, kappa, falpha, dnitot, dene, iinvqd, rminor, bt, rmajor, &
      ignite, aspect, qstar, q, afuel, plascur, pcoreradpv
    implicit none

    real(dp) :: fhz

    !  Arguments

    real(dp), intent(in) :: hhh

    !  Local variables

    real(dp) :: powerhtz,ptrez,ptriz,taueezz,taueiz,taueffz

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call pcond(afuel,palpmw,aspect,bt,dnitot,dene,dnla,eps,hhh, &
         iinvqd,iscz,ignite,kappa,kappa95,kappaa,pchargemw,pinjmw, &
         plascur,pcoreradpv,rmajor,rminor,te,ten,tin,q,qstar,vol, &
         xarea,zeff,ptrez,ptriz,taueezz,taueiz,taueffz,powerhtz)

    ! MDK All the scaling laws now contain hfact, so this code no longer required.
    !if (iscz < 3) then  !  only laws 1 and 2 are affected???
    !   ptrez = ptrez/hhh
    !   ptriz = ptriz/hhh
    !end if

    !  At power balance, fhz is zero.

    fhz = ptrez + ptriz - falpha*palppv - pchargepv - pohmpv

    !  Take into account whether injected power is included in tau_e
    !  calculation (i.e. whether device is ignited)

    if (ignite == 0) fhz = fhz - pinjmw/vol

    !  Include the radiation power if requested

    if (iradloss == 0) then
       fhz = fhz + pradpv
    else if (iradloss == 1) then
       fhz = fhz + pcoreradpv
    else
       continue
    end if

  end function fhz

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine outplas(outfile)

    !! Subroutine to output the plasma physics information
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer : Fortran output unit identifier
    !! This routine writes the plasma physics information
    !! to a file, in a tidy format.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constraint_variables, only: maxradwallload, peakradwallload, fbetatry, &
      taulimit, peakfactrad
    use current_drive_variables, only: bscf_nevins, bscfmax, cboot, &
      bscf_wilson, bscf_sauter, pinjmw, bscf_iter89, bootipf, pinjimw, pinjemw, &
      psipf, pscf_scene, diacf_hender, diacf_scene, diaipf
		use error_handling, only: fdiags, idiags, report_error
    use impurity_radiation_module, only: nimp, coreradiationfraction, &
      coreradius, fimp, impurity_arr_frac, impurity_arr_Label
    use physics_variables, only: ieped, ftar, dnelimt, fgwped, kappaa, deni, &
      betap, iculbl, rad_fraction_total, palpnb, ten, falpi, iradloss, pthrmw, &
      ralpne, taueff, dntau, dene, rad_fraction_sol, iprofile, rhopedn, &
      xarea, itart, epbetmax, neped, te0, ptrimw, dnbeta, powerht, psyncpv, &
      res_time, ignite, vol, bvert, tbeta, photon_wall, burnup, kappaa_ipb, &
      hfact, ilhthresh, alphan, fkzohm, alpha_crit, pohmmw, pouterzoneradmw, qlim, &
      qfuel, triang95, rplas, zeff, pdhe3, plascur, pdt, pdd, pbrempv, &
      ipedestal, dlamie, vsres, falpe, rli, ptremw, alphat, rminor, isc, &
      teped, fdeut, gamma, dnprot, ftrit, aion, btot, vsbrn, betanb, protium, &
      pchargemw, wallmw, vsstt, aspect, ti, q0, pinnerzoneradmw, &
      normalised_total_beta, pdivmax, dnbeam, kappa95, nesep_crit, fhe3, &
      triang, pneutmw, tauee, betalim, rlp, te, dlimit, ne0, qstar, dnalp, &
      taup, sarea, ti0, plhthresh, bp, dnitot, pradmw, pradsolmw, csawth, rndfuel, q95, &
      rhopedt, tauratio, pperim, tesep, vsind, ibss, alphaj, dnz, q, ssync, &
      psolradmw, tauei, ishape, plinepv, palpmw, palpfwmw, icurr, pdivt, &
      gammaft, powfmw
    use physics_variables, only: betaft, tauscl, fgwsep, rmajor, falpha, &
      nesep, facoh, kappa, dlimit, beta, dlimit, eps, pthrmw, dnla, bt, &
      pthrmw, pthrmw, pthrmw, idivrt, ips, idia
    use process_output, only: int_to_string2, ovarre, ovarrf, oheadr, &
      oblnkl, ovarin, ocmmnt, osubhd, ovarst
    use numerics, only: active_constraints, boundu, icc, &
        boundl, ioptimz
    use reinke_variables, only: fzactual, impvardiv, fzmin
	 use constants, only: rmu0, mproton, mfile, echarge, pi, epsilon0
    use stellarator_variables, only: iotabar, istell
    implicit none

    !  Arguments

    integer, intent(in) :: outfile

    !  Local variables

    real(dp) :: betath, tot_power_plasma, normalised_toroidal_beta
    ! pinj
    integer :: imp
    character(len=30) :: tauelaw
    character(len=30) :: str1,str2
    real(dp) :: fgwped_out ! neped/dlimit(7)
    real(dp) :: fgwsep_out ! nesep/dlimit(7)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Dimensionless plasma parameters. See reference below.
    nu_star = 1/rmu0  * (15.d0*echarge**4 * dlamie) / (4.d0*pi**1.5d0 * epsilon0**2) * &
              vol**2 * rmajor**2 * bt * sqrt(eps) * dnla**3 * kappa           / &
              (total_plasma_internal_energy**2 * plascur)

   rho_star = sqrt(2.d0* mproton * aion * total_plasma_internal_energy / (3.d0 * vol * dnla) ) / &
              (echarge * bt * eps * rmajor)

   beta_mcdonald = 4.d0/3.d0 *rmu0 * total_plasma_internal_energy / (vol * bt**2)

   call oheadr(outfile,'Plasma')

   if (istell == 0) then
      select case (idivrt)
      case (0)
         call ocmmnt(outfile,'Plasma configuration = limiter')
      case (1)
         call ocmmnt(outfile,'Plasma configuration = single null divertor')
      case (2)
         call ocmmnt(outfile,'Plasma configuration = double null divertor')
      case default
         idiags(1) = idivrt ; call report_error(85)
      end select
   else
      call ocmmnt(outfile,'Plasma configuration = stellarator')
   end if

   if (istell == 0) then
      if (itart == 0) then
         itart_r = itart
         call ovarrf(outfile,'Tokamak aspect ratio = Conventional, itart = 0','(itart)',itart_r)
      else if (itart == 1) then
         itart_r = itart
         call ovarrf(outfile,'Tokamak aspect ratio = Spherical, itart = 1','(itart)',itart_r)
      end if
   end if

   call osubhd(outfile,'Plasma Geometry :')
   call ovarrf(outfile,'Major radius (m)','(rmajor)',rmajor)
   call ovarrf(outfile,'Minor radius (m)','(rminor)',rminor, 'OP ')
   call ovarrf(outfile,'Aspect ratio','(aspect)',aspect)

   if (istell == 0) then

       select case (ishape)
       case (0,6,8)
          call ovarrf(outfile,'Elongation, X-point (input value used)', '(kappa)',kappa, 'IP ')
       case (1)
          call ovarrf(outfile,'Elongation, X-point (TART scaling)', '(kappa)',kappa, 'OP ')
       case (2,3)
          call ovarrf(outfile,'Elongation, X-point (Zohm scaling)', '(kappa)',kappa, 'OP ')
          call ovarrf(outfile,'Zohm scaling adjustment factor', '(fkzohm)',fkzohm)
       case (4,5,7)
          call ovarrf(outfile,'Elongation, X-point (calculated from kappa95)', '(kappa)',kappa, 'OP ')
       case (9)
          call ovarrf(outfile,'Elongation, X-point (calculated from aspect ratio and li(3))', &
               '(kappa)',kappa, 'OP ')
       case (10)
         call ovarrf(outfile,'Elongation, X-point (calculated from aspect ratio and stability margin)', &
         '(kappa)',kappa, 'OP ')
       case (11)
         call ovarrf(outfile,'Elongation, X-point (calculated from aspect ratio via Menard 2016)', &
         '(kappa)',kappa, 'OP ')
       case default
          idiags(1) = ishape ; call report_error(86)
       end select

       select case (ishape)
       case (4,5,7)
          call ovarrf(outfile,'Elongation, 95% surface (input value used)', &
               '(kappa95)',kappa95, 'IP ')
       case default
          call ovarrf(outfile,'Elongation, 95% surface (calculated from kappa)', &
               '(kappa95)',kappa95, 'OP ')
       end select

       call ovarrf(outfile,'Elongation, area ratio calc.','(kappaa)',kappaa, 'OP ')

       select case (ishape)
       case (0,2,6,8,9,10,11)
          call ovarrf(outfile,'Triangularity, X-point (input value used)', &
               '(triang)',triang, 'IP ')
       case (1)
          call ovarrf(outfile,'Triangularity, X-point (TART scaling)', &
               '(triang)',triang, 'OP ')
       case (3,4,5,7)
          call ovarrf(outfile,'Triangularity, X-point (calculated from triang95)', &
               '(triang)',triang, 'OP ')
       end select

       select case (ishape)
       case (3,4,5,7)
          call ovarrf(outfile,'Triangularity, 95% surface (input value used)', &
               '(triang95)',triang95, 'IP ')
       case default
          call ovarrf(outfile,'Triangularity, 95% surface (calculated from triang)', &
               '(triang95)',triang95, 'OP ')
       end select

       call ovarrf(outfile,'Plasma poloidal perimeter (m)','(pperim)',pperim, 'OP ')
    end if

    call ovarrf(outfile,'Plasma cross-sectional area (m2)','(xarea)',xarea, 'OP ')
    call ovarre(outfile,'Plasma surface area (m2)','(sarea)',sarea, 'OP ')
    call ovarre(outfile,'Plasma volume (m3)','(vol)',vol, 'OP ')

    call osubhd(outfile,'Current and Field :')


      if (istell == 0) then
         if (iprofile == 0) then
            call ocmmnt(outfile, &
               'Consistency between q0,q,alphaj,rli,dnbeta is not enforced')
         else
            call ocmmnt(outfile, &
               'Consistency between q0,q,alphaj,rli,dnbeta is enforced')
         end if
         call oblnkl(outfile)
         call ovarin(outfile,'Plasma current scaling law used','(icurr)',icurr)


         call ovarrf(outfile,'Plasma current (MA)','(plascur/1D6)',plascur/1.0D6, 'OP ')
         !call ovarrf(outfile,'Plasma current (A)','(plascur)',plascur, 'OP ')
         if (iprofile == 1) then
            call ovarrf(outfile,'Current density profile factor','(alphaj)',alphaj, 'OP ')
         else
            call ovarrf(outfile,'Current density profile factor','(alphaj)',alphaj)
         end if

         call ovarrf(outfile,'Plasma internal inductance, li','(rli)',rli, 'OP ')
         call ovarrf(outfile,'Vertical field at plasma (T)','(bvert)',bvert, 'OP ')
      end if

    call ovarrf(outfile,'Vacuum toroidal field at R (T)','(bt)',bt)
    call ovarrf(outfile,'Average poloidal field (T)','(bp)',bp, 'OP ')

    call ovarrf(outfile,'Total field (sqrt(bp^2 + bt^2)) (T)','(btot)',btot, 'OP ')


   if (istell == 0) then
       call ovarrf(outfile,'Safety factor on axis','(q0)',q0)

       if (icurr == 2) then
          call ovarrf(outfile,'Mean edge safety factor','(q)',q)
       end if

       call ovarrf(outfile,'Safety factor at 95% flux surface','(q95)',q95)

       call ovarrf(outfile,'Cylindrical safety factor (qcyl)','(qstar)',qstar, 'OP ')

       if (ishape == 1) then
          call ovarrf(outfile,'Lower limit for edge safety factor q', '(qlim)',qlim, 'OP ')
       end if
    else
       call ovarrf(outfile,'Rotational transform','(iotabar)',iotabar)
    end if

    call osubhd(outfile,'Beta Information :')

    betath = beta-betaft-betanb
    gammaft = (betaft + betanb)/betath

    call ovarre(outfile,'Total plasma beta','(beta)',beta)
    call ovarre(outfile,'Total poloidal beta','(betap)',betap, 'OP ')
    call ovarre(outfile,'Total toroidal beta',' ',beta*(btot/bt)**2, 'OP ')
    call ovarre(outfile,'Fast alpha beta','(betaft)',betaft, 'OP ')
    call ovarre(outfile,'Beam ion beta','(betanb)',betanb, 'OP ')
    call ovarre(outfile,'(Fast alpha + beam beta)/(thermal beta)','(gammaft)',gammaft, 'OP ')

    call ovarre(outfile,'Thermal beta',' ',betath, 'OP ')
    call ovarre(outfile,'Thermal poloidal beta',' ',betath*(btot/bp)**2, 'OP ')
    call ovarre(outfile,'Thermal toroidal beta (= beta-exp)',' ', betath*(btot/bt)**2, 'OP ')

    call ovarrf(outfile,'2nd stability beta : beta_p / (R/a)', '(eps*betap)',eps*betap, 'OP ')
    call ovarrf(outfile,'2nd stability beta upper limit','(epbetmax)', epbetmax)


    if (istell == 0) then
       if (iprofile == 1) then
            call ovarrf(outfile,'Beta g coefficient','(dnbeta)',dnbeta, 'OP ')
       else
            call ovarrf(outfile,'Beta g coefficient','(dnbeta)',dnbeta)
       end if

       call ovarrf(outfile,'Normalised thermal beta',' ',1.0D8*betath*rminor*bt/plascur, 'OP ')
       !call ovarrf(outfile,'Normalised total beta',' ',1.0D8*beta*rminor*bt/plascur, 'OP ')
       call ovarrf(outfile,'Normalised total beta',' ',normalised_total_beta, 'OP ')
       !call ovarrf(outfile,'Normalised toroidal beta',' ',normalised_total_beta*(btot/bt)**2, 'OP ')
       normalised_toroidal_beta=normalised_total_beta*(btot/bt)**2
       call ovarrf(outfile,'Normalised toroidal beta','(normalised_toroidal_beta)',normalised_toroidal_beta, 'OP ')
    end if


    if (iculbl == 0) then
       call ovarrf(outfile,'Limit on total beta','(betalim)',betalim, 'OP ')
    else if (iculbl == 1) then
       call ovarrf(outfile,'Limit on thermal beta','(betalim)',betalim, 'OP ')
    else
       call ovarrf(outfile,'Limit on thermal + NB beta','(betalim)', betalim, 'OP ')
    end if

    call ovarre(outfile,'Plasma thermal energy (J)',' ', 1.5D0*betath*btot*btot/(2.0D0*rmu0)*vol, 'OP ')

	call ovarre(outfile,'Total plasma internal energy (J)','(total_plasma_internal_energy)', total_plasma_internal_energy, 'OP ')

    call osubhd(outfile,'Temperature and Density (volume averaged) :')
    call ovarrf(outfile,'Electron temperature (keV)','(te)',te)
    call ovarrf(outfile,'Electron temperature on axis (keV)','(te0)',te0, 'OP ')
    call ovarrf(outfile,'Ion temperature (keV)','(ti)',ti)
    call ovarrf(outfile,'Ion temperature on axis (keV)','(ti0)',ti0, 'OP ')
    call ovarrf(outfile,'Electron temp., density weighted (keV)','(ten)',ten, 'OP ')
    call ovarre(outfile,'Electron density (/m3)','(dene)',dene)
    call ovarre(outfile,'Electron density on axis (/m3)','(ne0)',ne0, 'OP ')
    call ovarre(outfile,'Line-averaged electron density (/m3)','(dnla)',dnla, 'OP ')

    if (istell == 0) then
     call ovarre(outfile,'Line-averaged electron density / Greenwald density', &
         '(dnla_gw)',dnla/dlimit(7), 'OP ')
    end if


    call ovarre(outfile,'Ion density (/m3)','(dnitot)',dnitot, 'OP ')
    call ovarre(outfile,'Fuel density (/m3)','(deni)',deni, 'OP ')
    call ovarre(outfile,'Total impurity density with Z > 2 (no He) (/m3)','(dnz)',dnz, 'OP ')
    call ovarre(outfile,'Helium ion density (thermalised ions only) (/m3)','(dnalp)',dnalp, 'OP ')
    call ovarre(outfile,'Proton density (/m3)','(dnprot)',dnprot, 'OP ')
    if(protium > 1.0d-10)then
        call ovarre(outfile,'Seeded protium density / electron density','(protium)',protium)
    end if

    call ovarre(outfile,'Hot beam density (/m3)','(dnbeam)',dnbeam, 'OP ')
    call ovarre(outfile,'Density limit from scaling (/m3)','(dnelimt)',dnelimt, 'OP ')
    if ((ioptimz > 0).and.(active_constraints(5))) then
        call ovarre(outfile,'Density limit (enforced) (/m3)','(boundu(9)*dnelimt)',boundu(9)*dnelimt, 'OP ')
    end if
    call ovarre(outfile,'Helium ion density (thermalised ions only) / electron density','(ralpne)',ralpne)
    call oblnkl(outfile)


   call ocmmnt(outfile,'Impurities')
   call oblnkl(outfile)
   call ocmmnt(outfile,'Plasma ion densities / electron density:')
   do imp = 1,nimp
      ! MDK Update fimp, as this will make the ITV output work correctly.
      fimp(imp) = impurity_arr_frac(imp)
      str1 = impurity_arr_Label(imp) // ' concentration'
      str2 = '(fimp('//int_to_string2(imp)//'))'
      ! MDK Add output flag for H which is calculated.
      if (imp==1) then
        !call ovarre(outfile,str1,str2,impurity_arr_frac(imp), 'OP ')
        call ovarre(outfile,str1,str2,fimp(imp), 'OP ')
      else
        call ovarre(outfile,str1,str2,fimp(imp))
      end if
   end do

    call ovarre(outfile,'Average mass of all ions (amu)','(aion)',aion, 'OP ')
    ! MDK Say which impurity is varied, if iteration variable fimpvar (102) is turned on
    !if (any(ixc == 102)) then
    !    call ovarst(outfile,'Impurity used as an iteration variable' , '', '"' // impurity_arr(impvar)%label // '"')
    !    call ovarre(outfile,'Fractional density of variable impurity (ion / electron density)','(fimpvar)',fimpvar)
    !end if
    call oblnkl(outfile)
    call ovarrf(outfile,'Effective charge','(zeff)',zeff, 'OP ')

    ! Issue #487.  No idea what zeffai is.
    ! I haven't removed it as it is used in subroutine rether,
    !   (routine to find the equilibration power between the ions and electrons)
    ! call ovarrf(outfile,'Mass weighted effective charge','(zeffai)',zeffai, 'OP ')

    call ovarrf(outfile,'Density profile factor','(alphan)',alphan)
    call ovarin(outfile,'Plasma profile model','(ipedestal)',ipedestal)

    if(ipedestal.ge.1)then
        if (ne0<neped) then
            call report_error(213)
        end if
        call ocmmnt(outfile,'Pedestal profiles are used.')
        call ovarrf(outfile,'Density pedestal r/a location','(rhopedn)',rhopedn)
        if(fgwped >= 0d0)then
            call ovarre(outfile,'Electron density pedestal height (/m3)','(neped)',neped, 'OP ')
        else
            call ovarre(outfile,'Electron density pedestal height (/m3)','(neped)',neped)
        end if

        ! This code is ODD! Don't change it! No explanation why fgwped and fgwsep
        ! must be assigned to their exisiting values!
        fgwped_out = neped/dlimit(7)
        fgwsep_out = nesep/dlimit(7)
        if(fgwped >= 0d0) fgwped = neped/dlimit(7)
        if(fgwsep >= 0d0) fgwsep = nesep/dlimit(7)

        call ovarre(outfile,'Electron density at pedestal / nGW','(fgwped_out)',fgwped_out)
        call ovarrf(outfile,'Temperature pedestal r/a location','(rhopedt)',rhopedt)
        ! Issue #413 Pedestal scaling
        call ovarin(outfile,'Pedestal scaling switch','(ieped)',ieped)
        if(ieped==1)then
            call ocmmnt(outfile,'Saarelma 6-parameter pedestal temperature scaling is ON')

            if(eped_warning() /= '')then
                call ocmmnt(outfile,'WARNING: Pedestal parameters are outside the range of applicability of the scaling:')
                call ocmmnt(outfile,'triang: 0.4 - 0.6; kappa: 1.5 - 2.0;   plascur: 10 - 20 MA, rmajor: 7 - 11 m;')
                call ocmmnt(outfile,'rminor: 2 - 3.5 m; tesep: 0 - 0.5 keV; normalised_total_beta: 2 - 3; ')
                write(*,*)'WARNING: Pedestal parameters are outside the range of applicability of the scaling:'
                write(*,*)'triang: 0.4 - 0.6; kappa: 1.5 - 2.0;   plascur: 10 - 20 MA, rmajor: 7 - 11 m;'
                write(*,*)'rminor: 2 - 3.5 m; tesep: 0 - 0.5 keV; normalised_total_beta: 2 - 3'
                write(*,*)trim(eped_warning())
            endif
        endif
        call ovarrf(outfile,'Electron temp. pedestal height (keV)','(teped)',teped)
        if (any(icc == 78)) then
           call ovarrf(outfile,'Electron temp. at separatrix (keV)','(tesep)',tesep, 'OP ')
        else
           call ovarrf(outfile,'Electron temp. at separatrix (keV)','(tesep)',tesep)
        endif
        call ovarre(outfile,'Electron density at separatrix (/m3)','(nesep)',nesep)
        call ovarre(outfile,'Electron density at separatrix / nGW','(fgwsep_out)',fgwsep_out)

    endif

    ! Issue 558 - addition of constraint 76 to limit the value of nesep, in proportion with the ballooning parameter and Greenwald density
    if(any(icc==76))then
       call ovarre(outfile,'Critical ballooning parameter value','(alpha_crit)',alpha_crit)
       call ovarre(outfile,'Critical electron density at separatrix (/m3)','(nesep_crit)',nesep_crit)
    endif

    call ovarrf(outfile,'Temperature profile index','(alphat)',alphat)
    call ovarrf(outfile,'Temperature profile index beta','(tbeta)',tbeta)

   if (istell == 0) then
      call osubhd(outfile,'Density Limit using different models :')
      call ovarre(outfile,'Old ASDEX model','(dlimit(1))',dlimit(1), 'OP ')
      call ovarre(outfile,'Borrass ITER model I','(dlimit(2))',dlimit(2), 'OP ')
      call ovarre(outfile,'Borrass ITER model II','(dlimit(3))',dlimit(3), 'OP ')
      call ovarre(outfile,'JET edge radiation model','(dlimit(4))',dlimit(4), 'OP ')
      call ovarre(outfile,'JET simplified model','(dlimit(5))',dlimit(5), 'OP ')
      call ovarre(outfile,'Hugill-Murakami Mq model','(dlimit(6))',dlimit(6), 'OP ')
      call ovarre(outfile,'Greenwald model','(dlimit(7))',dlimit(7), 'OP ')
   end if

    call osubhd(outfile,'Fuel Constituents :')
    call ovarrf(outfile,'Deuterium fuel fraction','(fdeut)',fdeut)
    call ovarrf(outfile,'Tritium fuel fraction','(ftrit)',ftrit)
    if (fhe3 > 1.0D-3) call ovarrf(outfile,'3-Helium fuel fraction','(fhe3)',fhe3)

    call osubhd(outfile,'Fusion Power :')
    call ovarre(outfile,'Total fusion power (MW)','(powfmw)',powfmw, 'OP ')
    call ovarre(outfile,' =    D-T fusion power (MW)','(pdt)',pdt, 'OP ')
    call ovarre(outfile,'  +   D-D fusion power (MW)','(pdd)',pdd, 'OP ')
    call ovarre(outfile,'  + D-He3 fusion power (MW)','(pdhe3)',pdhe3, 'OP ')
    call ovarre(outfile,'Alpha power: total (MW)','(palpmw)',palpmw, 'OP ')
    call ovarre(outfile,'Alpha power: beam-plasma (MW)','(palpnb)',palpnb, 'OP ')
    call ovarre(outfile,'Neutron power (MW)','(pneutmw)',pneutmw, 'OP ')
    call ovarre(outfile,'Charged particle power (excluding alphas) (MW)', '(pchargemw)',pchargemw, 'OP ')
    tot_power_plasma=falpha*palpmw+pchargemw+pohmmw+pinjmw
    call ovarre(outfile,'Total power deposited in plasma (MW)','(tot_power_plasma)',tot_power_plasma, 'OP ')
    !call ovarre(outfile,'Total power deposited in plasma (MW)','()',falpha*palpmw+pchargemw+pohmmw+pinjmw, 'OP ')

    call osubhd(outfile,'Radiation Power (excluding SOL):')
    call ovarre(outfile,'Bremsstrahlung radiation power (MW)','(pbrempv*vol)', pbrempv*vol, 'OP ')
    call ovarre(outfile,'Line radiation power (MW)','(plinepv*vol)', plinepv*vol, 'OP ')
    call ovarre(outfile,'Synchrotron radiation power (MW)','(psyncpv*vol)', psyncpv*vol, 'OP ')
    call ovarrf(outfile,'Synchrotron wall reflectivity factor','(ssync)',ssync)
    call ovarre(outfile,"Normalised minor radius defining 'core'", '(coreradius)',coreradius)
    call ovarre(outfile,"Fraction of core radiation subtracted from P_L", &
         '(coreradiationfraction)',coreradiationfraction)
    call ovarre(outfile,'Radiation power from inner zone (MW)', '(pinnerzoneradmw)',pinnerzoneradmw, 'OP ')
    call ovarre(outfile,'Radiation power from outer zone (MW)','(pouterzoneradmw)', pouterzoneradmw, 'OP ')

    if (istell/=0) then
      call ovarre(outfile,'SOL radiation power as imposed by f_rad (MW)','(psolradmw)', psolradmw, 'OP ')
    end if

    call ovarre(outfile,'Total radiation power from inside LCFS (MW)','(pradmw)',pradmw, 'OP ')
    call ovarre(outfile,'LCFS radiation fraction = total radiation in LCFS / total power deposited in plasma', &
        '(rad_fraction_LCFS)', rad_fraction_LCFS, 'OP ')
    call ovarre(outfile,'Nominal mean radiation load on inside surface of reactor (MW/m2)', &
        '(photon_wall)', photon_wall, 'OP ')
    call ovarre(outfile,'Peaking factor for radiation wall load', &
        '(peakfactrad)', peakfactrad, 'IP ')
    call ovarre(outfile,'Maximum permitted radiation wall load (MW/m^2)', &
        '(maxradwallload)', maxradwallload, 'IP ')
    call ovarre(outfile,'Peak radiation wall load (MW/m^2)', &
        '(peakradwallload)', peakradwallload, 'OP ')
    call ovarre(outfile,'Fast alpha particle power incident on the first wall (MW)', &
        '(palpfwmw)', palpfwmw, 'OP ')
    call ovarre(outfile,'Nominal mean neutron load on inside surface of reactor (MW/m2)', &
        '(wallmw)', wallmw, 'OP ')

    if (istell == 0) then
      call oblnkl(outfile)
      call ovarre(outfile,'Power incident on the divertor targets (MW)', &
         '(ptarmw)',ptarmw, 'OP ')
      call ovarre(outfile, 'Fraction of power to the lower divertor', &
         '(ftar)', ftar, 'IP ')
      call ovarre(outfile,'Outboard side heat flux decay length (m)', &
         '(lambdaio)',lambdaio, 'OP ')
      if (idivrt == 2) then
         call ovarre(outfile,'Midplane seperation of the two magnetic closed flux surfaces (m)', &
            '(drsep)',drsep, 'OP ')
      end if
      call ovarre(outfile,'Fraction of power on the inner targets', &
         '(fio)',fio, 'OP ')
      call ovarre(outfile,'Fraction of power incident on the lower inner target', &
         '(fLI)',fLI, 'OP ')
      call ovarre(outfile,'Fraction of power incident on the lower outer target', &
         '(fLO)',fLO, 'OP ')
      if (idivrt == 2 ) then
         call ovarre(outfile,'Fraction of power incident on the upper inner target', &
         '(fUI)',fUI, 'OP ')
         call ovarre(outfile,'Fraction of power incident on the upper outer target', &
         '(fUO)',fUO, 'OP ')
      end if
      call ovarre(outfile,'Power incident on the lower inner target (MW)', &
         '(pLImw)',pLImw, 'OP ')
      call ovarre(outfile,'Power incident on the lower outer target (MW)', &
         '(pLOmw)',pLOmw, 'OP ')
      if (idivrt == 2) then
         call ovarre(outfile,'Power incident on the upper innner target (MW)', &
            '(pUImw)',pUImw, 'OP ')
         call ovarre(outfile,'Power incident on the upper outer target (MW)', &
            '(pUOmw)',pUOmw, 'OP ')
      end if
    end if

    call oblnkl(outfile)
    call ovarre(outfile,'Ohmic heating power (MW)','(pohmmw)',pohmmw, 'OP ')
    call ovarrf(outfile,'Fraction of alpha power deposited in plasma','(falpha)',falpha, 'OP ')
    call ovarrf(outfile,'Fraction of alpha power to electrons','(falpe)',falpe, 'OP ')
    call ovarrf(outfile,'Fraction of alpha power to ions','(falpi)',falpi, 'OP ')
    call ovarre(outfile,'Ion transport (MW)','(ptrimw)',ptrimw, 'OP ')
    call ovarre(outfile,'Electron transport (MW)','(ptremw)',ptremw, 'OP ')
    call ovarre(outfile,'Injection power to ions (MW)','(pinjimw)',pinjimw, 'OP ')
    call ovarre(outfile,'Injection power to electrons (MW)','(pinjemw)',pinjemw, 'OP ')
    if (ignite == 1) then
       call ocmmnt(outfile,'  (Injected power only used for start-up phase)')
    end if
    call ovarin(outfile,'Ignited plasma switch (0=not ignited, 1=ignited)', '(ignite)',ignite)

    call oblnkl(outfile)
    call ovarre(outfile,'Power into divertor zone via charged particles (MW)','(pdivt)',pdivt, 'OP ')

    if (pdivt <= 0.001D0) then
       fdiags(1) = pdivt ; call report_error(87)
       call oblnkl(outfile)
       call ocmmnt(outfile,'  BEWARE: possible problem with high radiation power')
       call ocmmnt(outfile,'          Power into divertor zone is unrealistic;')
       call ocmmnt(outfile,'          divertor calculations will be nonsense!')
       call ocmmnt(outfile,'  Set constraint 17 (Radiation fraction upper limit).')
       call oblnkl(outfile)
    end if

    if (idivrt == 2) then
      ! Double null divertor configuration
      call ovarre(outfile,'Pdivt / R ratio (MW/m) (On peak divertor)','(pdivmax/rmajor)',pdivmax/rmajor, 'OP ')
      call ovarre(outfile,'Pdivt Bt / qAR ratio (MWT/m) (On peak divertor)','(pdivmaxbt/qar)', ((pdivmax*bt)/(q95*aspect*rmajor)), 'OP ')
    else
      ! Single null divertor configuration
      call ovarre(outfile,'Psep / R ratio (MW/m)','(pdivt/rmajor)',pdivt/rmajor, 'OP ')
      call ovarre(outfile,'Psep Bt / qAR ratio (MWT/m)','(pdivtbt/qar)', ((pdivt*bt)/(q95*aspect*rmajor)), 'OP ')
    end if

      if (istell == 0) then
         call osubhd(outfile,'H-mode Power Threshold Scalings :')

         call ovarre(outfile,'ITER 1996 scaling: nominal (MW)','(pthrmw(1))', pthrmw(1), 'OP ')
         call ovarre(outfile,'ITER 1996 scaling: upper bound (MW)','(pthrmw(2))', pthrmw(2), 'OP ')
         call ovarre(outfile,'ITER 1996 scaling: lower bound (MW)','(pthrmw(3))', pthrmw(3), 'OP ')
         call ovarre(outfile,'ITER 1997 scaling (1) (MW)','(pthrmw(4))',pthrmw(4), 'OP ')
         call ovarre(outfile,'ITER 1997 scaling (2) (MW)','(pthrmw(5))',pthrmw(5), 'OP ')
         call ovarre(outfile,'Martin 2008 scaling: nominal (MW)', '(pthrmw(6))',pthrmw(6), 'OP ')
         call ovarre(outfile,'Martin 2008 scaling: 95% upper bound (MW)', '(pthrmw(7))',pthrmw(7), 'OP ')
         call ovarre(outfile,'Martin 2008 scaling: 95% lower bound (MW)', '(pthrmw(8))',pthrmw(8), 'OP ')
         call ovarre(outfile,'Snipes 2000 scaling: nominal (MW)', '(pthrmw(9))',pthrmw(9), 'OP ')
         call ovarre(outfile,'Snipes 2000 scaling: upper bound (MW)', '(pthrmw(10))',pthrmw(10), 'OP ')
         call ovarre(outfile,'Snipes 2000 scaling: lower bound (MW)', '(pthrmw(11))',pthrmw(11), 'OP ')
         call ovarre(outfile,'Snipes 2000 scaling (closed divertor): nominal (MW)', '(pthrmw(12))',pthrmw(12), 'OP ')
         call ovarre(outfile,'Snipes 2000 scaling (closed divertor): upper bound (MW)', '(pthrmw(13))',pthrmw(13), 'OP ')
         call ovarre(outfile,'Snipes 2000 scaling (closed divertor): lower bound (MW)', '(pthrmw(14))',pthrmw(14), 'OP ')
         call ovarre(outfile,'Hubbard 2012 L-I threshold - nominal (MW)', '(pthrmw(15))',pthrmw(15), 'OP ')
         call ovarre(outfile,'Hubbard 2012 L-I threshold - lower bound (MW)', '(pthrmw(16))',pthrmw(16), 'OP ')
         call ovarre(outfile,'Hubbard 2012 L-I threshold - upper bound (MW)', '(pthrmw(17))',pthrmw(17), 'OP ')
         call ovarre(outfile,'Hubbard 2017 L-I threshold', '(pthrmw(18))',pthrmw(18), 'OP ')
         call ovarre(outfile,'Martin 2008 aspect ratio corrected scaling: nominal (MW)', '(pthrmw(19))',pthrmw(19), 'OP ')
         call ovarre(outfile,'Martin 2008 aspect ratio corrected scaling: 95% upper bound (MW)', '(pthrmw(20))',pthrmw(20), 'OP ')
         call ovarre(outfile,'Martin 2008 aspect ratio corrected scaling: 95% lower bound (MW)', '(pthrmw(21))',pthrmw(21), 'OP ')
         call oblnkl(outfile)
         if ((ilhthresh.eq.9).or.(ilhthresh.eq.10).or.(ilhthresh.eq.11)) then
            if ((bt < 0.78D0).or.(bt > 7.94D0)) then
               call ocmmnt(outfile,'(bt outside Snipes 2000 fitted range)')
               call report_error(201)
            end if
            if ((rminor < 0.15D0).or.(rminor > 1.15D0)) then
               call ocmmnt(outfile,'(rminor outside Snipes 2000 fitted range)')
               call report_error(202)
            end if
            if ((rmajor < 0.55D0).or.(rmajor > 3.37D0)) then
               call ocmmnt(outfile,'(rmajor outside Snipes 2000 fitted range)')
               call report_error(203)
            end if
            if ((dnla < 0.09D20).or.(dnla > 3.16D20)) then
               call ocmmnt(outfile,'(dnla outside Snipes 2000 fitted range)')
               call report_error(204)
            end if
            if ((kappa < 1.0D0).or.(kappa > 2.04D0)) then
               call ocmmnt(outfile,'(kappa outside Snipes 2000 fitted range)')
               call report_error(205)
            end if
            if ((triang < 0.07D0).or.(triang > 0.74D0)) then
               call ocmmnt(outfile,'(triang outside Snipes 2000 fitted range)')
               call report_error(206)
            end if
         call oblnkl(outfile)
         end if
         if ((ilhthresh.eq.12).or.(ilhthresh.eq.13).or.(ilhthresh.eq.14)) then
            call ocmmnt(outfile,'(L-H threshold for closed divertor only. Limited data used in Snipes fit)')
            call oblnkl(outfile)
            call report_error(207)
         end if
         if ((ioptimz > 0).and.(active_constraints(15))) then
            call ovarre(outfile,'L-H threshold power (enforced) (MW)', '(boundl(103)*plhthresh)',boundl(103)*plhthresh, 'OP ')
            call ovarre(outfile,'L-H threshold power (MW)', '(plhthresh)',plhthresh, 'OP ')
         else
            call ovarre(outfile,'L-H threshold power (NOT enforced) (MW)', '(plhthresh)',plhthresh, 'OP ')
         end if
      end if

    call osubhd(outfile,'Confinement :')

    if (ignite == 1) then
       call ocmmnt(outfile, &
            'Device is assumed to be ignited for the calculation of confinement time')
       call oblnkl(outfile)
    end if

    write(outfile,200) tauscl(isc)
200 format(' Confinement scaling law',T45,A24)

    if (index(tauscl(isc),'(') /= 0) then
       tauelaw = '"'//trim(tauscl(isc)(1:index(tauscl(isc),'(',.true.)-1))//'"'
    else
       tauelaw = '"'//trim(tauscl(isc))//'"'
    end if
    call ovarst(mfile,'Confinement scaling law','(tauelaw)',trim(tauelaw))

    call ovarrf(outfile,'Confinement H factor','(hfact)',hfact)
    call ovarrf(outfile,'Global thermal energy confinement time (s)','(taueff)',taueff, 'OP ')
    call ovarrf(outfile,'Ion energy confinement time (s)','(tauei)',tauei, 'OP ')
    call ovarrf(outfile,'Electron energy confinement time (s)','(tauee)',tauee, 'OP ')
    call ovarre(outfile,'n.tau = Volume-average electron density x Energy confinement time (s/m3)', &
        '(dntau)', dntau, 'OP ')
    call ocmmnt(outfile,'Triple product = Vol-average electron density x Vol-average&
        & electron temperature x Energy confinement time:')
    call ovarre(outfile,'Triple product  (keV s/m3)','(dntau*te)',dntau*te, 'OP ')
    call ovarre(outfile,'Transport loss power assumed in scaling law (MW)', '(powerht)',powerht, 'OP ')
    call ovarin(outfile,'Switch for radiation loss term usage in power balance', '(iradloss)',iradloss)
    if (iradloss == 0) then
       call ovarre(outfile,'Radiation power subtracted from plasma power balance (MW)', '',pradmw, 'OP ')
       call ocmmnt(outfile,'  (Radiation correction is total radiation power)')
    else if (iradloss == 1) then
       call ovarre(outfile,'Radiation power subtracted from plasma power balance (MW)', '',pinnerzoneradmw, 'OP ')
       call ocmmnt(outfile,'  (Radiation correction is core radiation power)')
    else
       call ovarre(outfile,'Radiation power subtracted from plasma power balance (MW)', '',0.0D0)
       call ocmmnt(outfile,'  (No radiation correction applied)')
    end if
    call ovarrf(outfile,'Alpha particle confinement time (s)','(taup)',taup, 'OP ')
    ! Note alpha confinement time is no longer equal to fuel particle confinement time.
    call ovarrf(outfile,'Alpha particle/energy confinement time ratio','(taup/taueff)',taup/taueff, 'OP ')
    call ovarrf(outfile,'Lower limit on taup/taueff','(taulimit)',taulimit)
    call ovarrf(outfile,'Total energy confinement time including radiation loss (s)', &
         '(total_energy_conf_time)', total_energy_conf_time, 'OP ')
    call ocmmnt(outfile,'  (= stored energy including fast particles / loss power including radiation')

   if (istell == 0) then
      ! Issues 363 Output dimensionless plasma parameters MDK
      call osubhd(outfile,'Dimensionless plasma parameters')
      call ocmmnt(outfile,'For definitions see')
      call ocmmnt(outfile,'Recent progress on the development and analysis of the ITPA global H-mode confinement database')
      call ocmmnt(outfile,'D.C. McDonald et al, 2007 Nuclear Fusion v47, 147. (nu_star missing 1/mu0)')
      call ovarre(outfile,'Normalized plasma pressure beta as defined by McDonald et al', '(beta_mcdonald)',beta_mcdonald,'OP ')
      call ovarre(outfile,'Normalized ion Larmor radius', '(rho_star)', rho_star,'OP ')
      call ovarre(outfile,'Normalized collisionality', '(nu_star)',nu_star,'OP ')
      call ovarre(outfile,'Volume measure of elongation','(kappaa_IPB)',kappaa_IPB,'OP ')


      call osubhd(outfile,'Plasma Volt-second Requirements :')
      call ovarre(outfile,'Total volt-second requirement (Wb)','(vsstt)',vsstt, 'OP ')
      call ovarre(outfile,'Inductive volt-seconds (Wb)','(vsind)',vsind, 'OP ')
      call ovarrf(outfile,'Ejima coefficient','(gamma)',gamma)
      call ovarre(outfile,'Start-up resistive (Wb)','(vsres)',vsres, 'OP ')
      call ovarre(outfile,'Flat-top resistive (Wb)','(vsbrn)',vsbrn, 'OP ')

      call ovarrf(outfile,'bootstrap current fraction multiplier', '(cboot)',cboot)
      call ovarrf(outfile,'Bootstrap fraction (ITER 1989)', '(bscf_iter89)',bscf_iter89, 'OP ')


      call ovarrf(outfile,'Bootstrap fraction (Sauter et al)', '(bscf_sauter)',bscf_sauter, 'OP ')


      call ovarrf(outfile,'Bootstrap fraction (Nevins et al)', '(bscf_nevins)',bscf_nevins, 'OP ')
      call ovarrf(outfile,'Bootstrap fraction (Wilson)', '(bscf_wilson)',bscf_wilson, 'OP ')
      call ovarrf(outfile,'Diamagnetic fraction (Hender)', '(diacf_hender)',diacf_hender, 'OP ')
      call ovarrf(outfile,'Diamagnetic fraction (SCENE)', '(diacf_scene)',diacf_scene, 'OP ')
      call ovarrf(outfile,'Pfirsch-Schlueter fraction (SCENE)', '(pscf_scene)',pscf_scene, 'OP ')
      ! Error to catch if bootstap fraction limit has been enforced
      if (err242==1)then
         call report_error(242)
      end if
      ! Error to catch if self-driven current fraction limit has been enforced
      if (err243==1)then
         call report_error(243)
      end if

      if (bscfmax < 0.0D0) then
         call ocmmnt(outfile,'  (User-specified bootstrap current fraction used)')
      else if (ibss == 1) then
         call ocmmnt(outfile,'  (ITER 1989 bootstrap current fraction model used)')
      else if (ibss == 2) then
         call ocmmnt(outfile,'  (Nevins et al bootstrap current fraction model used)')
      else if (ibss == 3) then
         call ocmmnt(outfile,'  (Wilson bootstrap current fraction model used)')
      else if (ibss == 4) then
         call ocmmnt(outfile,'  (Sauter et al bootstrap current fraction model used)')
      end if

      if (idia == 0) then
         call ocmmnt(outfile,'  (Diamagnetic current fraction not calculated)')
         ! Error to show if diamagnetic current is above 1% but not used
         if (diacf_scene.gt.0.01D0) then
         call report_error(244)
         end if
      else if (idia == 1) then
         call ocmmnt(outfile,'  (Hender diamagnetic current fraction scaling used)')
      else if (idia == 2) then
         call ocmmnt(outfile,'  (SCENE diamagnetic current fraction scaling used)')

      if (ips == 0) then
            call ocmmnt(outfile,'  (Pfirsch-Schlüter current fraction not calculated)')
      else if (ips == 1) then
            call ocmmnt(outfile,'  (SCENE Pfirsch-Schlüter current fraction scaling used)')
      end if

      endif

      call ovarrf(outfile,'Bootstrap fraction (enforced)','(bootipf.)',bootipf, 'OP ')
      call ovarrf(outfile,'Diamagnetic fraction (enforced)','(diaipf.)',diaipf, 'OP ')
      call ovarrf(outfile,'Pfirsch-Schlueter fraction (enforced)','(psipf.)',psipf, 'OP ')

      call ovarre(outfile,'Loop voltage during burn (V)','(vburn)', plascur*rplas*facoh, 'OP ')
      call ovarre(outfile,'Plasma resistance (ohm)','(rplas)',rplas, 'OP ')

      call ovarre(outfile,'Resistive diffusion time (s)','(res_time)',res_time, 'OP ')
      call ovarre(outfile,'Plasma inductance (H)','(rlp)',rlp, 'OP ')
      call ovarrf(outfile,'Coefficient for sawtooth effects on burn V-s requirement','(csawth)',csawth)
   end if

    call osubhd(outfile,'Fuelling :')
    call ovarre(outfile,'Ratio of He and pellet particle confinement times','(tauratio)',tauratio)
    call ovarre(outfile,'Fuelling rate (nucleus-pairs/s)','(qfuel)',qfuel, 'OP ')
    call ovarre(outfile,'Fuel burn-up rate (reactions/s)','(rndfuel)',rndfuel, 'OP ')
    call ovarrf(outfile,'Burn-up fraction','(burnup)',burnup, 'OP ')


    if (any(icc == 78)) then
       call osubhd(outfile,'Reinke Criterion :')
       call ovarin(outfile,'index of impurity to be iterated for divertor detachment', '(impvardiv)',impvardiv)
       call ovarre(outfile,'Minimum Impurity fraction from Reinke','(fzmin)',fzmin, 'OP ')
       call ovarre(outfile,'Actual Impurity fraction','(fzactual)',fzactual)
    endif
  end subroutine outplas

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine outtim(outfile)

    !! Routine to print out the times of the various stages
    !! during a single plant cycle
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer : Fortran output unit identifier
    !! This routine writes out the times of the various stages
    !! during a single plant cycle.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use process_output, only: ovarrf, ovarre, oheadr, oblnkl
		use times_variables, only: tramp, theat, tcycle, tohs, tdwell, tqnch, tburn
    implicit none

    !  Arguments

    integer, intent(in) :: outfile

    !  Local variables

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call oheadr(outfile,'Times')

    call ovarrf(outfile,'Initial charge time for CS from zero current (s)','(tramp)', tramp)
    call ovarrf(outfile,'Plasma current ramp-up time (s)','(tohs)',tohs)
    call ovarrf(outfile,'Heating time (s)','(theat)',theat)
    call ovarre(outfile,'Burn time (s)','(tburn)',tburn, 'OP ')
    call ovarrf(outfile,'Reset time to zero current for CS (s)','(tqnch)',tqnch)
    call ovarrf(outfile,'Time between pulses (s)','(tdwell)',tdwell)
    call oblnkl(outfile)
    !call ovarre(outfile,'Pulse time (s)','(tpulse)',tpulse, 'OP ')
    !call ovarrf(outfile,'Down time (s)','(tdown)',tdown, 'OP ')
    call ovarre(outfile,'Total plant cycle time (s)','(tcycle)',tcycle, 'OP ')

  end subroutine outtim

end module physics_module
