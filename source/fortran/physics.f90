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

end module physics_module
