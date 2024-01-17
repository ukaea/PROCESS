module physics_functions_module

  !! Module containing physics subfunctions
  !! author: K Ellis, CCFE, Culham Science Centre
  !! N/A
  !! This module contains physics routines which can be called by physics or
  !! other modules
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  public :: beamfus, palph, palph2

  !  Module-level variables
  real(dp) :: vcritx

contains

  subroutine init_physics_functions
    !! Initialise module variables
    implicit none

    vcritx = 0.0D0
  end subroutine init_physics_functions

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine pthresh(dene,dnla,bt,rmajor,kappa,sarea,aion,aspect,pthrmw)
    !! author: P J Knight, CCFE, Culham Science Centre
    !! L-mode to H-mode power threshold calculation
    !! dene   : input real :  volume-averaged electron density (/m3)
    !! dnla   : input real :  line-averaged electron density (/m3)
    !! bt     : input real :  toroidal field on axis (T)
    !! rmajor : input real :  plasma major radius (m)
    !! kappa  : input real :  plasma elongation
    !! sarea  : input real :  plasma surface area (m**2)
    !! aion   : input real :  average mass of all ions (amu)
    !! aspect : input real :  aspect ratio
    !! pthrmw(17) : output real array : power threshold (different scalings)
    !! This routine calculates the power threshold for the L-mode to
    !! H-mode transition.
    !! ITER Physics Design Description Document, p.2-2
    !! ITER-FDR Plasma Performance Assessments, p.III-9
    !! Snipes, 24th EPS Conference, Berchtesgaden 1997, p.961
    !! Martin et al, 11th IAEA Tech. Meeting on H-mode Physics and
    !! Transport Barriers, Journal of Physics: Conference Series
    !! 123 (2008) 012033
    !! J A Snipes and the International H-mode Threshold Database
    !! Working Group, 2000, Plasma Phys. Control. Fusion, 42, A299
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use physics_variables, only: rminor, plascur
    implicit none

    !  Arguments

    real(dp), intent(in) :: dene,dnla,bt,rmajor,kappa,sarea,aion,aspect
    real(dp), dimension(21), intent(out) :: pthrmw

    !  Local variables

    real(dp) :: dene20,dnla20,marterr

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    dene20 = 1.0D-20*dene
    dnla20 = 1.0D-20*dnla

    !  ITER-DDD, D.Boucher
    !  Fit to 1996 H-mode power threshold database: nominal

    pthrmw(1) = 0.45D0 * dene20**0.75D0 * bt * rmajor**2

    !  Fit to 1996 H-mode power threshold database: upper bound

    pthrmw(2) = 0.37D0 * dene20 * bt * rmajor**2.5D0

    !  Fit to 1996 H-mode power threshold database: lower bound

    pthrmw(3) = 0.54D0 * dene20**0.5D0 * bt * rmajor**1.5D0

    !  J. A. Snipes, ITER H-mode Threshold Database Working Group,
    !  Controlled Fusion and Plasma Physics, 24th EPS Conference,
    !  Berchtesgaden, June 1997, vol.21A, part III, p.961

    pthrmw(4) = 0.65D0 * dnla20**0.93D0 * bt**0.86D0 * rmajor**2.15D0

    pthrmw(5) = 0.42D0 * dnla20**0.80D0 * bt**0.90D0 * rmajor**1.99D0 &
         * kappa**0.76D0

    !  Martin et al (2008) for recent ITER scaling, with mass correction
    !  and 95% confidence limits

    pthrmw(6) = 0.0488D0 * dnla20**0.717D0 * bt**0.803D0 &
         * sarea**0.941D0 * (2.0D0/aion)

    marterr = 0.057D0**2 + (0.035D0 * log(dnla20))**2 &
         + (0.032D0 * log(bt))**2 + (0.019D0 * log(sarea))**2
    marterr = sqrt(marterr) * pthrmw(6)

    pthrmw(7) = pthrmw(6) + 2.0D0*marterr
    pthrmw(8) = pthrmw(6) - 2.0D0*marterr

    ! Snipes et al (2000) scaling with mass correction
    ! Nominal, upper and lower

    pthrmw(9) = 1.42D0 * dnla20**0.58D0 * bt**0.82D0 * rmajor &
               * rminor**0.81D0 * (2.0D0/aion)

    pthrmw(10) = 1.547D0 * dnla20**0.615D0 * bt**0.851D0 &
              * rmajor**1.089D0 * rminor**0.876D0 * (2.0D0/aion)

    pthrmw(11) = 1.293D0 * dnla20**0.545D0 * bt**0.789D0 &
              * rmajor**0.911D0 * rminor**0.744D0 * (2.0D0/aion)

    ! Snipes et al (2000) scaling (closed divertor) with mass correction
    ! Nominal, upper and lower

    pthrmw(12) = 0.8D0 * dnla20**0.5D0 * bt**0.53D0 * rmajor**1.51D0 &
               * (2.0D0/aion)

    pthrmw(13) = 0.867D0 * dnla20**0.561D0 * bt**0.588D0 * rmajor**1.587D0 &
               * (2.0D0/aion)

    pthrmw(14) = 0.733D0 * dnla20**0.439D0 * bt**0.472D0 * rmajor**1.433D0 &
               * (2.0D0/aion)

    ! Hubbard et al. 2012 L-I threshold scaling

    ! Nominal
    pthrmw(15) = 2.11 * (plascur/1.0D6)**0.94 * dnla20**0.65

    ! Lower bound
    pthrmw(16) = 2.11 * (plascur/1.0D6)**0.70 * dnla20**0.47

    ! Upper bound
    pthrmw(17) = 2.11 * (plascur/1.0D6)**1.18 * dnla20**0.83

    ! Hubbard et al. 2017 L-I threshold scaling
    pthrmw(18) = 0.162 * dnla20 * sarea * (bt)**0.26

    !  Aspect ratio corrected Martin et al (2008)
    !  Correction: Takizuka 2004, Plasma Phys. Control Fusion 46 A227
    if (aspect.le.2.7D0) then
        pthrmw(19) = pthrmw(6) * (0.098D0 * aspect / (1.0D0 - (2.0D0/(1.0D0 + aspect))**0.5D0))
        pthrmw(20) = pthrmw(7) * (0.098D0 * aspect / (1.0D0 - (2.0D0/(1.0D0 + aspect))**0.5D0))
        pthrmw(21) = pthrmw(8) * (0.098D0 * aspect / (1.0D0 - (2.0D0/(1.0D0 + aspect))**0.5D0))
    else
        pthrmw(19) = pthrmw(6)
        pthrmw(20) = pthrmw(7)
        pthrmw(21) = pthrmw(8)
    end if

  end subroutine pthresh

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine palph(alphan,alphat,deni,fdeut,fhe3,ftrit,ti, &
       palppv,pchargepv,pneutpv,sigvdt,fusionrate,alpharate,protonrate, &
       pdtpv,pdhe3pv,pddpv)

    !! (Initial part of) fusion power and fast alpha pressure calculations
    !! author: P J Knight, CCFE, Culham Science Centre
    !! alphan     : input real :  density profile index
    !! alphat     : input real :  temperature profile index
    !! deni       : input real :  fuel ion density (/m3)
    !! fdeut      : input real :  deuterium fuel fraction
    !! fhe3       : input real :  helium-3 fuel fraction
    !! ftrit      : input real :  tritium fuel fraction
    !! ti         : input real :  ion temperature (keV)
    !! palppv     : output real : alpha particle fusion power per volume (MW/m3)
    !! pchargepv  : output real : other charged particle fusion power/volume (MW/m3)
    !! pneutpv    : output real : neutron fusion power per volume (MW/m3)
    !! sigvdt     : output real : profile averaged <sigma v DT> (m3/s)
    !! fusionrate : output real : fusion reaction rate (reactions/m3/s)
    !! alpharate  : output real : alpha particle production rate (/m3/s)
    !! protonrate : output real : proton production rate (/m3/s)
    !! pdtpv      : output real : D-T fusion power (MW/m3)
    !! pdhe3pv    : output real : D-He3 fusion power (MW/m3)
    !! pddpv      : output real : D-D fusion power (MW/m3)
    !! This subroutine numerically integrates over plasma cross-section to
    !! find the core plasma fusion power.
    !! T&amp;M/PKNIGHT/LOGBOOK24, p.6
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: echarge
    use maths_library, only: quanc8
    implicit none

    !  Arguments

    real(dp), intent(in) :: alphan, alphat, deni, fdeut, &
         fhe3, ftrit, ti
    real(dp), intent(out) :: palppv, pchargepv, pneutpv, sigvdt, &
         fusionrate, alpharate, protonrate, pdtpv, pdhe3pv, pddpv

    !  Local variables

    integer, parameter :: DT=1, DHE3=2, DD1=3, DD2=4
    integer :: ireaction,nofun
    real(dp) :: alow,arate,bhigh,epsq8,errest,etot,flag, &
         fpow,frate,pa,pc,pn,prate,sigmav

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Initialise local quantities
         alow = 0.0D0
         bhigh = 1.0D0
         epsq8 = 1.0D-9

         ! Find fusion power
         ! Integrate over plasma profiles to obtain fusion reaction rate
         palppv = 0.0D0
         pchargepv = 0.0D0
         pneutpv = 0.0D0
         fusionrate = 0.0D0
         alpharate = 0.0D0
         protonrate = 0.0D0
         pddpv = 0.0D0

         do ireaction = 1,4
             ! Fusion reaction rate (m3/s) is calculated in fint for each ireaction
             ! sigmav is the volume-averaged fusion reaction rate (m3/s)
             ! = integral(2 rho sigv(rho).ni(rho)^2 drho) / (deni**2)

             call quanc8(fint,alow,bhigh,epsq8,epsq8,sigmav,errest,nofun,flag)
             if (ireaction == DT) sigvdt = sigmav

             select case (ireaction)

                 case (DT)  ! D + T --> 4He + n reaction

                     etot = 17.59D0 * echarge  ! MJ
                     fpow = 1.0D0 * sigmav * etot * fdeut*ftrit * deni*deni  ! MW/m3
                     pa = 0.2D0 * fpow
                     pc = 0.0D0
                     pn = 0.8D0 * fpow
                     frate = fpow/etot  ! reactions/m3/second
                     arate = frate
                     prate = 0.0D0
                     pdtpv = fpow

                 case (DHE3)  ! D + 3He --> 4He + p reaction

                     etot = 18.35D0 * echarge  ! MJ
                     fpow = 1.0D0 * sigmav * etot * fdeut*fhe3 * deni*deni  ! MW/m3
                     pa = 0.2D0 * fpow
                     pc = 0.8D0 * fpow
                     pn = 0.0D0
                     frate = fpow/etot  ! reactions/m3/second
                     arate = frate
                     prate = frate      ! proton production /m3/second
                     pdhe3pv = fpow

                 case (DD1)  ! D + D --> 3He + n reaction

                     ! The 0.5 branching ratio is assumed to be included in sigmav
                     etot = 3.27D0 * echarge  ! MJ
                     fpow = 1.0D0 * sigmav * etot * 0.5D0*fdeut*fdeut * deni*deni  ! MW/m3
                     pa = 0.0D0
                     pc = 0.25D0 * fpow
                     pn = 0.75D0 * fpow
                     frate = fpow/etot  ! reactions/m3/second
                     arate = 0.0D0
                     prate = 0.0D0      ! Issue #557: No proton production
                     pddpv = pddpv + fpow

                 case (DD2)  !  D + D --> T + p reaction

                     ! The 0.5 branching ratio is assumed to be included in sigmav
                     etot = 4.03D0 * echarge  ! MJ
                     fpow = 1.0D0 * sigmav * etot * 0.5D0*fdeut*fdeut * deni*deni  ! MW/m3
                     pa = 0.0D0
                     pc = fpow
                     pn = 0.0D0
                     frate = fpow/etot  ! reactions/m3/second
                     arate = 0.0D0
                     prate = frate      ! proton production /m3/second
                     pddpv = pddpv + fpow

             end select

             palppv = palppv + pa
             pchargepv = pchargepv + pc
             pneutpv = pneutpv + pn
             fusionrate = fusionrate + frate
             alpharate = alpharate + arate
             protonrate = protonrate + prate

         end do

     contains

         ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function fint(rho)

      !! Integrand for fusion power integration
      !! author: P J Knight, CCFE, Culham Science Centre
      !! rho : input real :  Abscissa of the integration, = normalised
      !! plasma minor radius (0.0 <= rho < 1.0)
      !! This function evaluates the integrand for the fusion power
      !! integration, performed using routine
      !! <A HREF="quanc8.html">QUANC8</A>
      !! in routine <A HREF="palph.html">PALPH</A>.
      !! The fusion reaction assumed is controlled by flag
      !! <CODE>ireaction</CODE> set in <CODE>PALPH</CODE>.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use physics_variables, only: te, rhopedt, te0, teped, tesep, tbeta, &
        dene, rhopedn, ne0, neped, nesep
      use profiles_module, only: tprofile, nprofile
      implicit none

      real(dp) :: fint

      !  Arguments

      real(dp), intent(in) :: rho

      !  Local variables

      real(dp) :: nprof, nprofsq, sigv, tiofr

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Local ion temperature (keV) at r/a = rho

      tiofr = ti/te * tprofile(rho,rhopedt,te0,teped,tesep,alphat,tbeta)

      !  Fusion reaction rate (m3/s)

      sigv = bosch_hale(tiofr,ireaction)

      !  Integrand for the volume averaged fusion reaction rate sigmav:
      !  sigmav = integral(2 rho (sigv(rho) ni(rho)^2) drho),
      !  divided by the square of the volume-averaged ion density
      !  to retain the dimensions m3/s (this is multiplied back in later)

      nprof = 1.0D0/dene * nprofile(rho,rhopedn,ne0,neped,nesep,alphan)
      nprofsq = nprof*nprof

      fint = 2.0D0 * rho * sigv * nprofsq

    end function fint

  end subroutine palph

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine palph2(bt,bp,dene,deni,dnitot,falpe,falpi,palpnb, &
       ifalphap,pchargepv,pneutpv,ten,tin,vol,palpmw,pneutmw,pchargemw, &
       betaft,palppv,palpipv,palpepv,pfuscmw,powfmw)

    !! (Concluding part of) fusion power and fast alpha pressure
    !! calculations
    !! author: P J Knight, CCFE, Culham Science Centre
    !! bp       : input real :  poloidal field (T)
    !! bt       : input real :  toroidal field on axis (T)
    !! dene     : input real :  electron density (/m3)
    !! deni     : input real :  fuel ion density (/m3)
    !! dnitot   : input real :  total ion density (/m3)
    !! falpe    : input real :  fraction of alpha energy to electrons
    !! falpi    : input real :  fraction of alpha energy to ions
    !! ifalphap : input integer :  switch for fast alpha pressure method
    !! palpnb   : input real :  alpha power from hot neutral beam ions (MW)
    !! pchargepv : input real : other charged particle fusion power/volume (MW/m3)
    !! pneutpv  : input/output real : neutron fusion power per volume (MW/m3)
    !! ten      : input real :  density-weighted electron temperature (keV)
    !! tin      : input real :  density-weighted ion temperature (keV)
    !! vol      : input real :  plasma volume (m3)
    !! palpmw   : output real : alpha power (MW)
    !! pneutmw  : output real : neutron fusion power (MW)
    !! pchargemw : output real : other charged particle fusion power (MW)
    !! betaft   : output real : fast alpha beta component
    !! palppv   : input/output real : alpha power per volume (MW/m3)
    !! palpepv  : output real : alpha power per volume to electrons (MW/m3)
    !! palpipv  : output real : alpha power per volume to ions (MW/m3)
    !! pfuscmw  : output real : charged particle fusion power (MW)
    !! powfmw   : output real : fusion power (MW)
    !! This subroutine completes the calculation of the fusion power
    !! fast alpha pressure, and determines other alpha particle quantities.
    !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
    !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
    !! D J Ward, UKAEA Fusion: F/PL/PJK/PROCESS/CODE/050
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: echarge, rmu0
    use physics_variables, only: falpha, fdeut
    implicit none

    !  Arguments

    integer, intent(in) :: ifalphap
    real(dp), intent(in) :: bp, bt, dene, deni, dnitot, falpe, &
         falpi, palpnb, pchargepv, ten, tin, vol
    real(dp), intent(inout) :: palppv, pneutpv
    real(dp), intent(out) :: palpmw, pneutmw, pchargemw, betaft, palpepv, &
         palpipv, pfuscmw, powfmw

    !  Local variables

    real(dp) :: betath, fact, fact2, palppv_no_nb

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Store the palppv value without the NB alpha power
    palppv_no_nb = palppv

    !  Add neutral beam alpha power / volume
    palppv = palppv + palpnb/vol

    !  Add extra neutron power
    pneutpv = pneutpv + 4.0D0*palpnb/vol

    !  Total alpha power
    palpmw = palppv*vol

    !  Total non-alpha charged particle power
    pchargemw = pchargepv*vol

    !  Total neutron power
    pneutmw = pneutpv*vol

    !  Total fusion power
    powfmw = palpmw + pneutmw + pchargemw

    !  Charged particle fusion power
    pfuscmw = palpmw + pchargemw

    !  Alpha power to electrons and ions (used with electron
    !  and ion power balance equations only)
    !  No consideration of pchargepv here...
    palpipv = falpha * palppv*falpi
    palpepv = falpha * palppv*falpe

    !  Determine average fast alpha density
    if (fdeut < 1.0D0) then

       betath = 2.0D3*rmu0*echarge * (dene*ten + dnitot*tin)/(bt**2 + bp**2)

       ! jlion: This "fact" model is heavily flawed for smaller temperatures! It is unphysical for a stellarator (high n low T)
       ! IPDG89 fast alpha scaling
       if (ifalphap == 0) then
          fact = min( 0.30D0, &
               0.29D0*(deni/dene)**2 * ( (ten+tin)/20.0D0 - 0.37D0) )

       ! Modified scaling, D J Ward
       else
          fact = min( 0.30D0, &
               0.26D0*(deni/dene)**2 * &
               sqrt( max(0.0D0, ((ten+tin)/20.0D0 - 0.65D0)) ) )
       end if

       fact = max(fact,0.0D0)
       fact2 = palppv / palppv_no_nb
       betaft = betath * fact*fact2

    else  !  negligible alpha production, palppv = palpnb = 0
       betaft = 0.0D0
    end if

  end subroutine palph2

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function bosch_hale(t,reaction)

    !! Routine to calculate the fusion reaction rate
    !! author: R Kemp, CCFE, Culham Science Centre
    !! author: P J Knight, CCFE, Culham Science Centre
    !! t : input real : Maxwellian density-weighted ion temperature (keV)
    !! reaction : input integer : flag for fusion reaction to use:
    !! 1 : D-T reaction
    !! 2 : D-3He reaction
    !! 3 : D-D 1st reaction (50% probability)
    !! 4 : D-D 2nd reaction (50% probability)
    !! This routine calculates the volumetric fusion reaction rate
    !! <I>&lt;sigma v&gt;</I> in m3/s for one of four nuclear reactions,
    !! using the Bosch-Hale parametrization.
    !! <P>The valid range of the fit is 0.2 keV < t < 100 keV
    !! Bosch and Hale, Nuclear Fusion 32 (1992) 611-631
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(dp) :: bosch_hale

    !  Arguments

    real(dp), intent(in) :: t
    integer, intent(in) :: reaction

    !  Local variables

    integer, parameter :: DT=1, DHE3=2, DD1=3, DD2=4
    real(dp) :: theta1, theta, xi
    real(dp), dimension(4) :: bg, mrc2
    real(dp), dimension(4,7) :: cc

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if  (t == 0.0D0) then
       bosch_hale = 0.0D0
       return
    end if

    !  Gamov constant, BG

    bg(DT)   = 34.3827D0  !  D + T --> 4He + n reaction
    bg(DHE3) = 68.7508D0  !  D + 3He --> 4He + p reaction
    bg(DD1)  = 31.3970D0  !  D + D --> 3He + n reaction
    bg(DD2)  = 31.3970D0  !  D + D --> T + p reaction

    !  Reduced mass of the particles, keV

    mrc2(DT)   = 1.124656D6
    mrc2(DHE3) = 1.124572D6
    mrc2(DD1)  = 0.937814D6
    mrc2(DD2)  = 0.937814D6

    !  Parametrization coefficients

    cc(DT,1) =  1.17302D-9
    cc(DT,2) =  1.51361D-2
    cc(DT,3) =  7.51886D-2
    cc(DT,4) =  4.60643D-3
    cc(DT,5) =  1.35000D-2
    cc(DT,6) = -1.06750D-4
    cc(DT,7) =  1.36600D-5

    cc(DHE3,1) =  5.51036D-10
    cc(DHE3,2) =  6.41918D-3
    cc(DHE3,3) = -2.02896D-3
    cc(DHE3,4) = -1.91080D-5
    cc(DHE3,5) =  1.35776D-4
    cc(DHE3,6) =  0.00000D0
    cc(DHE3,7) =  0.00000D0

    cc(DD1,1) =  5.43360D-12
    cc(DD1,2) =  5.85778D-3
    cc(DD1,3) =  7.68222D-3
    cc(DD1,4) =  0.00000D0
    cc(DD1,5) = -2.96400D-6
    cc(DD1,6) =  0.00000D0
    cc(DD1,7) =  0.00000D0

    cc(DD2,1) =  5.65718D-12
    cc(DD2,2) =  3.41267D-3
    cc(DD2,3) =  1.99167D-3
    cc(DD2,4) =  0.00000D0
    cc(DD2,5) =  1.05060D-5
    cc(DD2,6) =  0.00000D0
    cc(DD2,7) =  0.00000D0

    theta1 = t*(cc(reaction,2) + t*(cc(reaction,4) + t*cc(reaction,6))) / &
         (1.0D0 + t*(cc(reaction,3) + t*(cc(reaction,5) + t*cc(reaction,7))))
    theta = t/(1.0D0 - theta1)

    xi = ((bg(reaction)**2)/(4.0D0*theta))**0.3333333333D0

    !  Volumetric reaction rate <sigma v> (m3/s)

    bosch_hale = 1.0D-6 * cc(reaction,1) * theta * &
         sqrt( xi/(mrc2(reaction)*t**3) ) * exp(-3.0D0*xi)

  end function bosch_hale

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine beamfus(beamfus0,betbm0,bp,bt,cnbeam,dene,deni,dlamie, &
       ealphadt,enbeam,fdeut,ftrit,ftritbm,sigvdt,ten,tin,vol,zeffai, &
       betanb,dnbeam2,palpnb)

    !! Routine to calculate beam slowing down properties
    !! author: P J Knight, CCFE, Culham Science Centre
    !! beamfus0: input real : multiplier for beam-background fusion calculation
    !! betbm0 : input real :  leading coefficient for neutral beam beta fraction
    !! bp     : input real :  poloidal field (T)
    !! bt     : input real :  toroidal field on axis (T)
    !! cnbeam : input real :  neutral beam current (A)
    !! dene   : input real :  electron density (/m3)
    !! deni   : input real :  fuel ion density (/m3)
    !! dlamie : input real :  ion-electron coulomb logarithm
    !! ealphadt : input real :  alpha particle birth energy (D-T) (keV)
    !! enbeam : input real :  neutral beam energy (keV)
    !! fdeut  : input real :  deuterium fraction of main plasma
    !! ftrit  : input real :  tritium fraction of main plasma
    !! ftritbm: input real :  tritium fraction of neutral beam
    !! sigvdt : input real :  profile averaged <sigma v> for D-T (m3/s)
    !! ten    : input real :  density weighted average electron temperature (keV)
    !! tin    : input real :  density weighted average ion temperature (keV)
    !! vol    : input real :  plasma volume (m3)
    !! zeffai : input real :  mass weighted plasma effective charge
    !! betanb : output real : neutral beam beta component
    !! dnbeam2: output real : hot beam ion density (/m3)
    !! palpnb : output real : alpha power from hot neutral beam ions (MW)
    !! This routine calculates the beam slowing down properties.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    real(dp), intent(in) :: beamfus0, betbm0, bp, bt, cnbeam, &
         dene, deni, dlamie, ealphadt, enbeam, fdeut, ftrit, ftritbm, &
         sigvdt, ten, tin, vol, zeffai
    real(dp), intent(out) :: betanb, dnbeam2, palpnb

    !  Local variables

    real(dp) :: denid,denit,ecritd,ecritt,ehotnb,palpdb, &
         palptb,tausl

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Velocity slowing down time

    tausl = 1.99D19 * (2.0D0*(1.0D0-ftritbm) + 3.0D0*ftritbm) * &
         ten**1.5D0 / dene / dlamie

    !  Critical energy for electron/ion slowing down of the beam ion
    !  (deuterium and tritium neutral beams, respectively) (keV)

    ecritd = 14.8D0 * ten * 2.0D0 * zeffai**0.6666D0 * &
         (dlamie+4.0D0)/dlamie
    ecritt = ecritd * 1.5D0

    !  Deuterium and tritium ion densities

    denid = deni * fdeut
    denit = deni * ftrit

    !  Perform beam calculations

    call beamcalc(denid,denit,ealphadt,enbeam,ecritd,ecritt,tausl, &
         ftritbm,cnbeam,tin,vol,sigvdt,palpdb,palptb,dnbeam2,ehotnb)

    !  Neutral beam alpha power

    palpnb = beamfus0 * (palpdb + palptb)

    !  Neutral beam beta

    betanb = betbm0 * 4.03D-22 * 0.66666D0 * dnbeam2 * ehotnb / &
         (bt**2 + bp**2)

  end subroutine beamfus

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine beamcalc(nd,nt,ealphadt,ebeam,ecritd,ecritt,tausbme, &
       ftritbm,ibeam,ti,vol,svdt,palfdb,palftb,nhot,ehot)

    !! Neutral beam alpha power and ion energy
    !! author: P J Knight, CCFE, Culham Science Centre
    !! ealphadt : input real :  alpha particle birth energy (D-T) (keV)
    !! ebeam  : input real :  beam energy (keV)
    !! ecritd : input real :  critical energy for electron/ion slowing down of
    !! the beam ion (deuterium neutral beam) (keV)
    !! ecritt : input real :  critical energy for beam slowing down
    !! (tritium neutral beam) (keV)
    !! ftritbm: input real :  beam tritium fraction (0.0 = deuterium beam)
    !! ibeam  : input real :  beam current (A)
    !! nd     : input real :  thermal deuterium density (/m3)
    !! nt     : input real :  thermal tritium density   (/m3)
    !! svdt   : input real :  profile averaged <sigma v> for D-T (m3/s)
    !! tausbme: input real :  beam ion slowing down time on electrons (s)
    !! ti     : input real :  thermal ion temperature (keV)
    !! vol    : input real :  plasma volume (m3) (95% flux surface)
    !! ehot   : output real : average hot beam ion energy (keV)
    !! nhot   : output real : hot beam ion density (/m3)
    !! palfdb : output real : alpha power from deut. beam-background fusion (MW)
    !! palftb : output real : alpha power from trit. beam-background fusion (MW)
    !! This routine calculates the neutral beam alpha power and ion energy.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        use constants, only: echarge, mproton

        implicit none

        ! Arguments
        real(kind(1.0D0)), intent(in) :: ealphadt, ebeam, ecritd, ecritt, &
            ftritbm, ibeam, nd, nt, svdt, tausbme, ti, vol
        real(kind(1.0D0)), intent(out) :: ehot, nhot, palfdb, palftb

        ! Local variables
        integer :: iabm
        real(kind(1.0D0)) :: ebmratd, ebmratt, ehotd, ehott, ifbmd, ifbmt, &
            ndhot, nhotmsd, nhotmst, nthot, presd, prest, s0d, s0t, svdhotn, &
            svthotn, tauseffd, tausefft, vcds, vcritd, vcritt, vcts, xcoefd, &
            xcoeft
        real(kind(1.0D0)) :: atmd, atmt, epsabs, epsrel

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Initialise shared variables
        atmd = 2.0D0   !  atomic mass of deuterium
        atmt = 3.0D0   !  atomic mass of tritium
        epsabs = 1.0D-7  !  absolute error
        epsrel = 1.0D-7  !  relative error

        ! D and T beam current fractions
        ifbmd = ibeam * (1.0D0 - ftritbm)
        ifbmt = ibeam * ftritbm

        ebmratd = ebeam/ecritd
        vcritd = sqrt(2.0D0*echarge*1000.0D0*ecritd/(mproton*atmd))
        tauseffd = tausbme/3.0D0 * log(1.0D0+(ebmratd)**1.5D0)
        nhotmsd = (1.0D0-ftritbm) * ibeam * tauseffd/(echarge * vol)

        ebmratt = ebeam/ecritt
        vcritt = sqrt(2.0D0*echarge*1000.0D0*ecritt/(mproton*atmt))
        tausefft = tausbme/3.0D0 * log(1.0D0+(ebmratt)**1.5D0)
        nhotmst = ftritbm * ibeam * tausefft/(echarge * vol)

        nhot = nhotmsd + nhotmst
        ndhot = nhotmsd
        nthot = nhotmst

        ! Average hot ion energy from Deng & Emmert, UWFDM-718, Jan 87
        vcds = 2.0D0 * ecritd * echarge * 1000.0D0/(2.0D0 * mproton)
        vcts = 2.0D0 * ecritt * echarge * 1000.0D0/(3.0D0 * mproton)

        s0d = ifbmd/(echarge * vol)
        s0t = ifbmt/(echarge * vol)

        xcoefd = atmd * mproton * tausbme * vcds * s0d / &
            (echarge * 1000.0D0 * 3.0D0)
        xcoeft = atmt * mproton * tausbme * vcts * s0t / &
            (echarge * 1000.0D0 * 3.0D0)

        presd = xcoefd * xbrak(ebeam,ecritd)
        prest = xcoeft * xbrak(ebeam,ecritt)

        ehotd = 1.5D0 * presd/ndhot
        ehott = 1.5D0 * prest/nthot
        ehot = (ndhot*ehotd + nthot*ehott)/nhot

        iabm = 2 ; svdhotn = 1.0D-4 * sgvhot(iabm,vcritd,ebeam)
        iabm = 3 ; svthotn = 1.0D-4 * sgvhot(iabm,vcritt,ebeam)

        palfdb = palphabm(ealphadt,ndhot,nt,svdhotn,vol,ti,svdt)
        palftb = palphabm(ealphadt,nthot,nd,svthotn,vol,ti,svdt)

    contains

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        function xbrak(e0,ec)
            !! Hot ion energy parameter
            !! author: P J Knight, CCFE, Culham Science Centre
            !! e0 : input real :  neutral beam energy (keV)
            !! ec : input real :  critical energy for electron/ion slowing down of
            !! the beam ion (keV)
            !! This routine calculates something to do with the hot ion energy...
            !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
            !
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            implicit none

            real(dp) :: xbrak

            ! Arguments
            real(dp), intent(in) :: e0, ec

      real(dp) :: ans,t1,t2,t3,t4,xarg,xc,xcs

            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            xcs = e0/ec
            xc = sqrt(xcs)

            t1 = xcs/2.0D0
            t2 = (log((xcs + 2.0D0*xc + 1.0D0)/(xcs - xc + 1.0D0)))/6.0D0

            xarg = (2.0D0*xc -1.0D0)/sqrt(3.0D0)
            t3 = (atan(xarg))/sqrt(3.0D0)
            t4 = 0.3022999D0

            ans = t1 + t2 - t3 - t4
            xbrak = ans

        end function xbrak

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function palphabm(ealphadt,nbm,nblk,sigv,vol,ti,svdt)

      !! Alpha power from beam-background fusion
      !! author: P J Knight, CCFE, Culham Science Centre
      !! ealphadt : input real :  alpha particle birth energy (D-T) (keV)
      !! nblk   : input real :  thermal ion density (/m3)
      !! nbm    : input real :  hot beam ion density (/m3)
      !! sigv   : input real :  hot beam fusion reaction rate (m3/s)
      !! svdt   : input real :  profile averaged <sigma v> for D-T (m3/s)
      !! ti     : input real :  thermal ion temperature (keV)
      !! vol    : input real :  plasma volume (m3)
      !! This routine calculates the alpha power from
      !! beam-background fusion.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: palphabm

      !  Arguments

      real(dp), intent(in) :: ealphadt,nblk,nbm,sigv,svdt,ti,vol

      !  Local variables

      real(dp) :: ratio

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ratio = svdt / bosch_hale(ti,1)

      palphabm = echarge/1000.0D0 * nbm * nblk * sigv * ealphadt * vol * ratio

    end function palphabm

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function sgvhot(iabm,vcrx,ebeam)

      !! Hot beam fusion reaction rate
      !! author: P J Knight, CCFE, Culham Science Centre
      !! ebeam  : input real :  neutral beam energy (keV)
      !! iabm   : input integer : switch denoting type of ion (2=D,3=T)
      !! vcrx   : input real :  critical velocity for electron/ion slowing down of
      !! the beam ion (m/s)
      !! This routine calculates the hot beam fusion reaction rate in m3/s.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use error_handling, only: idiags, report_error
      use maths_library, only: quanc8

      implicit none

      real(dp) :: sgvhot

      !  Arguments
      integer, intent(in) :: iabm
      real(dp), intent(in) :: ebeam, vcrx

      !  Local variables

      integer :: nofun
      real(dp) :: abm,abserr,epsabs1,flag,svint,t1,t2, &
           vbeam,vbeams,xv

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      epsabs1 = 1.0D-33

      if (iabm == 2) then
         abm = atmd
      else if (iabm == 3) then
         abm = atmt
      else
         idiags(1) = iabm ; call report_error(84)
      end if

      !  Initialise global variables

      vcritx = vcrx

      !  Beam velocity

      vbeams = ebeam * echarge * 1000.0D0 * 2.0D0/(abm * mproton)
      vbeam = sqrt(vbeams)

      xv = vbeam/vcrx
      t1 = 3.0D0 * vcrx/log(1.0D0+(xv**3))

      call quanc8(fsv,0.0D0,xv,epsabs1,epsrel,svint,abserr,nofun,flag)
      t2 = svint

      sgvhot = t1 * t2

    end function sgvhot

    end subroutine beamcalc

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function fsv(u)

    !! Integrand function for the hot beam fusion reaction rate
    !! author: P J Knight, CCFE, Culham Science Centre
    !! u : input real : abscissa of integration, = ratio of beam velocity
    !! to the critical velocity
    !! This is the integrand function for the hot beam fusion reaction rate.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: mproton, echarge

    implicit none

    real(dp) :: fsv

    !  Arguments

    real(dp), intent(in) :: u

    !  Local variables

    real(dp) :: t1,t2,xvc,xvcs

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    t1 = (u**3)/(1.0D0+u**3)

    !  vcritx : critical velocity for electron/ion slowing down of beam ion (m/s)

    xvc = vcritx*u
    xvcs = xvc * xvc * mproton/(echarge * 1000.0D0)
    t2 = sigbmfus(xvcs)

    fsv = t1 * t2

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function sigbmfus(vrelsq)

      !! Fusion reaction cross-section
      !! author: P J Knight, CCFE, Culham Science Centre
      !! vrelsq : input real :  square of the speed of the beam ion (keV/amu)
      !! This function evaluates the fusion reaction cross-section as a
      !! function of beam ion velocity (squared).
      !! The functional form of the cross-section is in terms of the equivalent
      !! deuterium energy, i.e. for a tritium beam at 500 keV the energy
      !! used in the cross-section function is 333 keV.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: sigbmfus

      !  Arguments

      real(dp), intent(in) :: vrelsq

      !  Local variables

      real(dp) :: a1,a2,a3,a4,a5,atmd,ebm,t1,t2

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      a1 = 45.95D0
      a2 = 5.02D4
      a3 = 1.368D-2
      a4 = 1.076D0
      a5 = 4.09D2

      !  Deuterium atomic mass

      atmd = 2.0D0

      !  Beam kinetic energy

      ebm = 0.5D0 * atmd * vrelsq

      !  Set limits on cross-section at low and high beam energies

      if (ebm < 10.0D0) then
         sigbmfus = 1.0D-27
      else if (ebm > 1.0D4) then
         sigbmfus = 8.0D-26
      else
         t1 = a2/(1.0D0 + (a3 * ebm - a4)**2) + a5
         t2 = ebm * (exp (a1/sqrt(ebm)) - 1.0D0)
         sigbmfus = 1.0D-24 * t1/t2
      end if

    end function sigbmfus

  end function fsv


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  function plasma_elongation_IPB()
    !! author: H Lux (UKAEA)
    !!
    !! Volume measure of plasma elongation using the IPB definition
    !!
    !! See Otto Kardaun et al 2008 Nucl. Fusion 48 099801

    ! Module variables
    use physics_variables, only : vol, rminor, rmajor
    use constants, only : pi

    real(dp) :: plasma_elongation_IPB
    !! Plasma elongation (IPB)

    plasma_elongation_IPB = vol / ( 2.0D0 * pi*pi * rminor*rminor * rmajor )
    !! \begin{equation} \kappa_{IPB} = \frac{V}{2\pi a^2 R_0} \end{equation}
    !!
    !! - \( V \) -- Plasma volume [m\(^3\)]
    !! - \( a \) -- Plasma minor radius [m]
    !! - \( R_0 \) -- Plasma major radius [m]

  end function plasma_elongation_IPB

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function total_mag_field()
    !! author: J. Morris (UKAEA)
    !!
    !! Calculates the total magnetic field

    ! Module variables
    use physics_variables, only : bt, bp

     ! Return value
     real(dp) :: total_mag_field

    total_mag_field = sqrt(bt**2 + bp**2)
    !! \begin{equation} B_{tot} = \sqrt{B_T^2 + B_p^2} \end{equation}

  end function total_mag_field

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function beta_poloidal()
    !! author: J. Morris (UKAEA)
    !!
    !! Calculates total poloidal beta

    ! Module variables
    use physics_variables, only : btot, bp, beta

     ! Return value
     real(dp) :: beta_poloidal

    beta_poloidal = beta * ( btot/bp )**2
    !! \begin{equation} \beta_p = \beta \left( \frac{B_{tot}}{B_p} \right)^2 \end{equation}
    !! See J.P. Freidberg, "Plasma physics and fusion energy", Cambridge University Press (2007)
    !! Page 270 ISBN 0521851076

  end function beta_poloidal

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function res_diff_time()
    !! author: J. Morris (UKAEA)
    !!
    !! Calculates resistive diffusion time

    ! Module variables
    use physics_variables, only : rmajor, rplas, kappa95
    use constants, only : rmu0

    ! Return value
    real(dp) :: res_diff_time

    res_diff_time = 2.0D0*rmu0*rmajor / (rplas*kappa95)
    !! Resistive diffusion time equals the current penetration time which is approximated by:
    !! \begin{equation} t_{\text{res-diff}} \sim
    !! \frac{2\mu_0.R_0}{\rho_{\text{plasma}}\kappa_{95}}\end{equation}
    !!
    !! * \( \mu_0 \) -- permittivity of free space [H/m]
    !! * \( R_0 \) -- plasma major radius [m]
    !! - \( \rho_{\text{plasma}} \) -- plasma resistivity [Ohms]
    !! - \( \kappa_{95} \) -- plasma elongation at 95% flux surface
    !!
    !! #TODO Reference needed

  end function res_diff_time

end module physics_functions_module
