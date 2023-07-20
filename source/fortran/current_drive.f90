! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module current_drive_module
  !! Module containing current drive system routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains routines relevant for calculating the
  !! current drive parameters for a fusion power plant.
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Import modules
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cudriv(outfile,iprint)

    !! Routine to calculate the current drive power requirements
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer : output file unit
    !! iprint : input integer : switch for writing to output file (1=yes)
    !! This routine calculates the power requirements of the current
    !! drive system, using a choice of models for the current drive
    !! efficiency.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use error_handling, only: idiags, report_error
    use profiles_module, only: tprofile, nprofile
    use process_output, only: oblnkl, ocmmnt, ovarin, ovarre, ovarrf, osubhd, &
      oheadr
    use heat_transport_variables, only: pinjwpfix, pinjwp
    use current_drive_variables, only: echpwr, pnbeam, plhybd, cnbeam, porbitlossmw, &
      iefrf, iefrffix, pheat, pheatfix, pinjfixmw, irfcd, feffcd, fpion, nbshinef, &
      gamcd, gamma_ecrh, rho_ecrh, etalh, etacd, etacdfix, etaech, forbitloss, &
      pinjmw, pwpnb, etanbi, enbeam, effcd, pwplh, echwpow, pnbitot, nbshinemw, &
      pinjemw, pinjimw, bigq, bootipf, bscfmax, taubeam, pinjalw, nbshield, &
      frbeam, rtanbeam, rtanmax, diaipf, psipf, plasipf, harnum, xi_ebw
    use physics_variables, only: dene, te, rmajor, ten, zeff, dlamee, beta, &
      rhopedt, rhopedn, te0, teped, tesep, alphat, alphan, ne0, nesep, neped, &
      bt, rminor, tbeta, plascur, ipedestal, faccd, ignite, pohmmw, powfmw, &
      facoh, fvsbrnni
    use constants, only: nout, echarge, emass, pi, epsilon0
    use cost_variables, only: startupratio

    implicit none

    ! Arguments
    integer, intent(in) :: iprint, outfile

    ! Local variables !
    ! !!!!!!!!!!!!!!!!!!

    real(dp) :: dene20, effnbss, effrfss, gamnb, gamrf, power1
    real(dp) :: effcdfix, effrfssfix, effnbssfix, pinjwp1
    real(dp) :: pnbitotfix, nbshinemwfix, porbitlossmwfix, cnbeamfix
    real(dp) :: pinjimw1, pinjemw1, pinjimwfix, pinjemwfix, pinjmw1, pinjmwfix
    real(dp) :: auxiliary_cdfix, faccdfix, gamcdfix
    real(dp) :: fshift, xf, enpa,ftherm,fpp,cdeff, ampperwatt
    real(dp) :: dens_at_rho, te_at_rho
    logical :: Temperature_capped
    real(dp) :: auxiliary_cd
    real(dp) :: a, fc, fp, density_factor

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    echpwr = 0.0D0
    pnbeam = 0.0D0
    plhybd = 0.0D0
    cnbeam = 0.0D0
    cnbeamfix = 0.0D0
    porbitlossmw = 0.0D0
    porbitlossmwfix = 0.0D0

    pinjmw1 = 0.0
    pinjmwfix = 0.0
    pinjimw1 = 0.0
    pinjemw1 = 0.0
    pinjemwfix = 0.0
    pinjimwfix = 0.0
    auxiliary_cdfix = 0.0
    faccdfix = 0.0
    gamcdfix = 0.0D0

    ! To stop issues with input file we force
    ! zero secondary heating if no injection method
    if (iefrffix == 0) then
      pheatfix = 0.0
    end if

    ! check for unphysically large heating in
    ! secondary injected power source
    if (pheatfix > pinjfixmw) then
      pheatfix = pinjfixmw
    end if

    ! irfcd |  switch for current drive calculation
    ! = 0   |  turned off
    ! = 1   |  turned on
    if (irfcd /= 0) then

       ! put electron density in desired units (10^-20 m-3)
       dene20 = dene * 1.0D-20

       ! If present we must calculate second current drive
       ! efficiencies in units of Amps/Watt using the fixed
       ! values from user input
       ! iefrffix |  switch for fixed current drive efficiency model
       select case (iefrffix)

       case(0) ! second current drive not present

       case (1)  ! Fenstermacher Lower Hybrid model

          effrfssfix = (0.36D0 * (1.0D0 + (te/25.0D0)**1.16D0)) / &
               (rmajor*dene20) * feffcd
          effcdfix = effrfssfix

       case (2)  ! Ion-Cyclotron current drive

          effrfssfix = 0.63D0 * 0.1D0*ten / (2.0D0 + zeff) / &
               (rmajor*dene20) * feffcd
          effcdfix = effrfssfix

       case (3)  ! Fenstermacher Electron Cyclotron Resonance model

          effrfssfix = 0.21D0 * ten/ (rmajor * dene20 * dlamee) * feffcd
          effcdfix = effrfssfix

       case (4)  ! Ehst Lower Hybrid / Fast Wave current drive

          effrfssfix = te**0.77D0 * (0.034D0 + 0.196D0 * beta) / &
               (rmajor*dene20) * ( 32.0D0/(5.0D0+zeff) + 2.0D0 + &
               (12.0D0*(6.0D0+zeff))/(5.0D0+zeff)/(3.0D0+zeff) + &
               3.76D0/zeff) / 12.507D0 * feffcd
          effcdfix = effrfssfix

       case (5)  ! ITER Neutral Beam current drive

          call iternb(effnbss,fpion,nbshinef)
          effnbssfix = effnbss * feffcd
          effcdfix = effnbssfix

       case (6)  ! Culham Lower Hybrid current drive model

          call cullhy(effrfss)
          effrfssfix = effrfss * feffcd
          effcdfix = effrfssfix

       case (7)  ! Culham ECCD model

          call culecd(effrfss)
          effrfssfix = effrfss * feffcd
          effcdfix = effrfssfix

       case (8)  ! Culham Neutral Beam model

          call culnbi(effnbss,fpion,nbshinef)
          effnbssfix = effnbss * feffcd
          effcdfix = effnbssfix

       case (9)  ! Issue #508 Remove RFP option  Oscillating Field CD model

       case (10)  ! ECRH user input gamma

          !  Normalised current drive efficiency gamma
          gamcd = gamma_ecrh

          ! Absolute current drive efficiency
          effrfssfix = gamcd / (dene20 * rmajor)
          effcdfix = effrfssfix

       case (11)  ! ECRH Poli model "HARE", removed in issue #1811

       case (12)
       ! EBW scaling
       ! Scaling author Simon Freethy
       ! Ref : PROCESS issue 1262

          !  Normalised current drive efficiency gamma
          gamcd = (xi_ebw/32.7D0) * te

          ! Absolute current drive efficiency
          effrfssfix = gamcd / (dene20 * rmajor)
          effcdfix = effrfssfix

          ! EBWs can only couple to plasma if cyclotron harmonic is above plasma density cut-off;
          !  this behaviour is captured in the following function (ref issue #1262):
          ! harnum = cyclotron harmonic number (fundamental used as default)
          ! constant 'a' controls sharpness of transition
          a = 0.1D0

          fc = 1.0D0/(2.0D0*pi) * harnum * echarge * bt / emass
          fp = 1.0D0/(2.0D0*pi) * sqrt( dene * echarge**2 / ( emass * epsilon0 ) )

          density_factor = 0.5D0 * ( 1.0D0 + tanh( (2.0D0/a) * ( ( fp - fc )/fp - a) ) )

          effcdfix = effcdfix * density_factor

          effrfssfix = effrfssfix * density_factor


       case default
          idiags(1) = iefrffix
          call report_error(126)

       end select

       ! find the current drive wall plug and injected powers (MW)
       ! and efficiencies for secondary current drive mechanisms
       ! using the fixed injected power from the user input
       select case (iefrffix)

       case (1,2,4,6)  ! LHCD or ICCD

          !  Injected power
          pinjemwfix = pinjfixmw

          !  Wall plug power
          pinjwpfix = pinjfixmw / etalh

          !  Wall plug to injector efficiency
          etacdfix = etalh

          !  Normalised current drive efficiency gamma
          gamcdfix = effrfssfix * (dene20 * rmajor)

          ! the fixed auxiliary current
          auxiliary_cdfix = effrfssfix * ( pinjfixmw - pheatfix) * 1.0d6
          faccdfix = auxiliary_cdfix / plascur

       case (3,7,10,11,12)  ! ECCD

          !  Injected power
          pinjemwfix = pinjfixmw

          !  Wall plug power
          pinjwpfix = pinjfixmw / etaech

          !  Wall plug to injector efficiency
          etacdfix = etaech

          ! the fixed auxiliary current
          auxiliary_cdfix = effrfssfix * ( pinjfixmw - pheatfix) * 1.0d6
          faccdfix = auxiliary_cdfix / plascur

       case (5,8)  ! NBCD

          ! Account for first orbit losses
          ! (power due to particles that are ionised but not thermalised) [MW]:
          ! This includes a second order term in shinethrough*(first orbit loss)
          forbitloss = min(0.999,forbitloss) ! Should never be needed

          if(ipedestal.ne.3) then  ! When not using PLASMOD
             pnbitotfix = pinjfixmw / (1.0D0-forbitloss+forbitloss*nbshinef)
          else
             ! Netural beam power calculated by PLASMOD
             pnbitotfix = pinjmw / (1.0D0-forbitloss+forbitloss*nbshinef)
          endif

          ! Shinethrough power (atoms that are not ionised) [MW]:
          nbshinemwfix = pnbitotfix * nbshinef

          ! First orbit loss
          porbitlossmwfix = forbitloss * (pnbitotfix - nbshinemwfix)

          ! Power deposited
          pinjmwfix = pnbitotfix - nbshinemwfix - porbitlossmwfix
          pinjimwfix = pinjmwfix * fpion
          pinjemwfix = pinjmwfix * (1.0D0-fpion)

          pwpnb = pnbitotfix/etanbi ! neutral beam wall plug power
          pinjwpfix = pwpnb
          etacdfix = etanbi
          gamnb = effnbssfix * (dene20 * rmajor)
          gamcdfix = gamnb
          cnbeamfix = 1.0D-3 * (pnbitotfix*1.0D6) / enbeam !  Neutral beam current (A)
          auxiliary_cdfix = effnbssfix * ( pinjfixmw - pheatfix) * 1.0d6
          faccdfix = auxiliary_cdfix / plascur

       case (9)  ! OFCD - RFP option removed in PROCESS (issue #508)

       end select

       ! Calculate current drive efficiencies in units of Amps/Watt.
       ! iefrf |  switch for current drive efficiency model
       select case (iefrf)

       case (1)  ! Fenstermacher Lower Hybrid model

          effrfss = (0.36D0 * (1.0D0 + (te/25.0D0)**1.16D0)) / &
               (rmajor*dene20) * feffcd
          effcd = effrfss

       case (2)  ! Ion-Cyclotron current drive

          effrfss = 0.63D0 * 0.1D0*ten / (2.0D0 + zeff) / &
               (rmajor*dene20) * feffcd
          effcd = effrfss

       case (3)  ! Fenstermacher Electron Cyclotron Resonance model

          effrfss = 0.21D0 * ten/ (rmajor * dene20 * dlamee) * feffcd
          effcd = effrfss

       case (4)  ! Ehst Lower Hybrid / Fast Wave current drive

          effrfss = te**0.77D0 * (0.034D0 + 0.196D0 * beta) / &
               (rmajor*dene20) * ( 32.0D0/(5.0D0+zeff) + 2.0D0 + &
               (12.0D0*(6.0D0+zeff))/(5.0D0+zeff)/(3.0D0+zeff) + &
               3.76D0/zeff) / 12.507D0 * feffcd
          effcd = effrfss

       case (5)  ! ITER Neutral Beam current drive

          call iternb(effnbss,fpion,nbshinef)
          effnbss = effnbss * feffcd
          effcd = effnbss

       case (6)  ! Culham Lower Hybrid current drive model

          call cullhy(effrfss)
          effrfss = effrfss * feffcd
          effcd = effrfss

       case (7)  ! Culham ECCD model

          call culecd(effrfss)
          effrfss = effrfss * feffcd
          effcd = effrfss

       case (8)  ! Culham Neutral Beam model

          call culnbi(effnbss,fpion,nbshinef)
          effnbss = effnbss * feffcd
          effcd = effnbss

       case (9)  ! Issue #508 Remove RFP option  Oscillating Field CD model

       case (10)  ! ECRH user input gamma

          gamcd = gamma_ecrh

          ! Absolute current drive efficiency
          effrfss = gamcd / (dene20 * rmajor)
          effcd = effrfss

       case (11)  ! ECRH Poli model "HARE", removed in issue #1811

       case (12)
       ! EBW scaling
       ! Scaling author Simon Freethy
       ! Ref : PROCESS issue 1262

          !  Normalised current drive efficiency gamma
          gamcd = (xi_ebw/32.7D0) * te

          ! Absolute current drive efficiency
          effrfss = gamcd / (dene20 * rmajor)
          effcd = effrfss

          ! EBWs can only couple to plasma if cyclotron harmonic is above plasma density cut-off;
          !  this behaviour is captured in the following function (ref issue #1262):
          ! harnum = cyclotron harmonic number (fundamental used as default)
          ! contant 'a' controls sharpness of transition

          !TODO is the below needed?
         !  a = 0.1D0

         !  fc = 1.0D0/(2.0D0*pi) * harnum * echarge * bt / emass
         !  fp = 1.0D0/(2.0D0*pi) * sqrt( dene * echarge**2 / ( emass * epsilon0 ) )

         !  density_factor = 0.5D0 * ( 1.0D0 + tanh( (2.0D0/a) * ( ( fp - fc )/fp - a) ) )

         !  effcd = effcd * density_factor

         !  effrfss = effrfss * density_factor


       case default
          idiags(1) = iefrf
          call report_error(126)

       end select

       ! Compute current drive wall plug and injected powers (MW) and efficiencies
       auxiliary_cd = faccd * plascur
       select case (iefrf)

       case (1,2,4,6)  ! LHCD or ICCD

          !  Injected power
          plhybd = 1.0D-6 * (faccd - faccdfix) * plascur / effrfss + pheat
          pinjimw1 = 0.0D0
          pinjemw1 = plhybd

          !  Wall plug power
          pwplh = plhybd / etalh
          pinjwp1 = pwplh

          !  Wall plug to injector efficiency
          etacd = etalh

          !  Normalised current drive efficiency gamma
          gamrf = effrfss * (dene20 * rmajor)
          gamcd = gamrf

       case (3,7,10,11,12)  ! ECCD

          !  Injected power (set to close to close the Steady-state current equilibrium)
          echpwr = 1.0D-6 * (faccd - faccdfix) * plascur / effrfss + pheat
          pinjemw1 = echpwr

          !  Wall plug power
          echwpow = echpwr / etaech

          !  Wall plug to injector efficiency
          pinjwp1 = echwpow
          etacd = etaech

       case (5,8)  ! NBCD

          ! MDK. See Gitlab issue #248, and scanned note.
          power1 = 1.0D-6 * (faccd - faccdfix) * plascur / effnbss + pheat

          ! Account for first orbit losses
          ! (power due to particles that are ionised but not thermalised) [MW]:
          ! This includes a second order term in shinethrough*(first orbit loss)
          forbitloss = min(0.999,forbitloss) ! Should never be needed

          if(ipedestal.ne.3)then  ! When not using PLASMOD
             pnbitot = power1 / (1.0D0-forbitloss+forbitloss*nbshinef)
          else
             ! Neutral beam power calculated by PLASMOD
             pnbitot = pinjmw / (1.0D0-forbitloss+forbitloss*nbshinef)
          endif

          ! Shinethrough power (atoms that are not ionised) [MW]:
          nbshinemw = pnbitot * nbshinef

          ! First orbit loss
          porbitlossmw = forbitloss * (pnbitot - nbshinemw)

          ! Power deposited
          pinjmw1 = pnbitot - nbshinemw - porbitlossmw
          pinjimw1 = pinjmw1 * fpion
          pinjemw1 = pinjmw1 * (1.0D0-fpion)

          pwpnb = pnbitot/etanbi ! neutral beam wall plug power
          pinjwp1 = pwpnb
          etacd = etanbi
          gamnb = effnbss * (dene20 * rmajor)
          gamcd = gamnb
          cnbeam = 1.0D-3 * (pnbitot*1.0D6) / enbeam !  Neutral beam current (A)



       case (9)  ! OFCD - RFP option removed in PROCESS (issue #508)

       end select

       ! Total injected power
       ! sum contributions from primary and secondary systems
       pinjmw = pinjemw1 + pinjimw1 + pinjemwfix + pinjimwfix
       pinjmw1 = pinjemw1 + pinjimw1
       pinjmwfix = pinjemwfix + pinjimwfix
       pinjemw = pinjemw1 + pinjemwfix
       pinjimw = pinjimw1 + pinjimwfix
       pinjwp = pinjwp1 + pinjwpfix

       ! Reset injected power to zero for ignited plasma (fudge)
       if (ignite == 1) then
           pinjwp = 0.0D0
       end if

       ! Ratio of fusion to input (injection+ohmic) power
       if (abs(pinjmw + porbitlossmw + pohmmw) < 1.0D-6) then
          bigq = 1.0D18
       else
          bigq = powfmw / (pinjmw + porbitlossmw + pohmmw)
       end if

    end if

    ! Output !
    ! !!!!!!!!!

    if (iprint == 0) return

    call oheadr(outfile,'Current Drive System')

    if (irfcd == 0) then
       call ocmmnt(outfile,'No current drive used')
       call oblnkl(outfile)
       return
    end if

    select case (iefrf)

    case (1,4,6)
       call ocmmnt(outfile,'Lower Hybrid Current Drive')
    case (2)
       call ocmmnt(outfile,'Ion Cyclotron Current Drive')
    case (3,7)
       call ocmmnt(outfile,'Electron Cyclotron Current Drive')
    case (5,8)
       call ocmmnt(outfile,'Neutral Beam Current Drive')
    case (9)
       ! RFP option removed in PROCESS (issue #508)
    case (10)
       call ocmmnt(outfile,'Electron Cyclotron Current Drive (user input gamma_CD)')
    case (11)
       call ocmmnt(outfile,'Electron Cyclotron Current Drive (HARE). No longer implemented.')
    case (12)
       call ocmmnt(outfile,'EBW current drive')
    end select

    call ovarin(outfile,'Current drive efficiency model','(iefrf)',iefrf)

    if (iefrffix.NE.0) then
      select case (iefrffix)

      case (1,4,6)
          call ocmmnt(outfile,'Lower Hybrid Current Drive')
      case (2)
          call ocmmnt(outfile,'Ion Cyclotron Current Drive')
      case (3,7)
          call ocmmnt(outfile,'Electron Cyclotron Current Drive')
      case (5,8)
          call ocmmnt(outfile,'Neutral Beam Current Drive')
      case (9)
          ! RFP option removed in PROCESS (issue #508)
      case (10)
          call ocmmnt(outfile,'Electron Cyclotron Current Drive (user input gamma_CD)')
      case(11)
          call ocmmnt(outfile,'Electron Cyclotron Current Drive (HARE). No longer implemented.')
      case (12)
          call ocmmnt(outfile,'EBW current drive')
      end select

      call ovarin(outfile,'Secondary current drive efficiency model','(iefrffix)',iefrffix)
    end if

    if (ignite == 1) then
       call ocmmnt(outfile, &
            'Ignited plasma; injected power only used for start-up phase')
    end if

    call oblnkl(outfile)

    if (abs(facoh) > 1.0D-8) then
       call ocmmnt(outfile,'Current is driven by both inductive')
       call ocmmnt(outfile,'and non-inductive means.')
    end if

    call ovarre(outfile,'Ratio of power for flat-top to start-up (MW)', '(startupratio)', startupratio)
    call ovarre(outfile,'Auxiliary power used for plasma heating only (MW)', '(pheat)', pheat + pheatfix)
    call ovarre(outfile,'Power injected for current drive (MW)','(pcurrentdrivemw)', pinjmw - pheat - pheatfix)
    call ovarre(outfile,'Maximum Allowed Bootstrap current fraction', '(bscfmax)', bscfmax)
    if (iefrffix.NE.0) then
      call ovarre(outfile,'Power injected for main current drive (MW)','(pcurrentdrivemw1)', pinjmw1 - pheat)
      call ovarre(outfile,'Power injected for secondary current drive (MW)','(pcurrentdrivemw2)', pinjmwfix - pheatfix)
    end if
    call ovarre(outfile,'Fusion gain factor Q','(bigq)',bigq, 'OP ')
    call ovarre(outfile,'Auxiliary current drive (A)','(auxiliary_cd)',auxiliary_cd, 'OP ')
    if (iefrffix.ne.0) then
      call ovarre(outfile,'Secondary auxiliary current drive (A)','(auxiliary_cdfix)',auxiliary_cdfix, 'OP ')
    end if
    call ovarre(outfile,'Current drive efficiency (A/W)','(effcd)',effcd, 'OP ')
    call ovarre(outfile,'Normalised current drive efficiency, gamma (10^20 A/W-m2)', &
         '(gamcd)',gamcd, 'OP ')
    call ovarre(outfile,'Wall plug to injector efficiency','(etacd)',etacd)

    select case(iefrf) ! Select plasma heating choice for output in the mfile
    ! Other heating cases to be added when necessary

    case(1,2,3,4)

    case(5,8) !NBI  case(5) and Case(8) NBI is covered elsewhere

    case(6,7,9,11)

    case(10) !ECRH/ECCD Heating of Plasma

          call ovarre(outfile,'ECRH plasma heating efficiency','(gamma_ecrh)',gamma_ecrh)

    case(12) !EBW Heating of Plasma

          call ovarre(outfile,'EBW plasma heating efficiency','(xi_ebw)',xi_ebw)

    case default
          idiags(1) = iefrf
          call report_error(126)

    end select


    if (iefrffix.NE.0) then
      call ovarre(outfile,'Secondary current drive efficiency (A/W)','(effcdfix)',effcdfix, 'OP ')
      call ovarre(outfile,'Seconday wall plug to injector efficiency','(etacdfix)',etacdfix)
      call ovarre(outfile,'Normalised secondary current drive efficiency, gamma (10^20 A/W-m2)', &
         '(gamcdfix)',gamcdfix, 'OP ')
    end if

    call osubhd(outfile,'Fractions of current drive :')
    call ovarrf(outfile,'Bootstrap fraction','(bootipf)',bootipf, 'OP ')
    call ovarrf(outfile,'Diamagnetic fraction','(diaipf)',diaipf, 'OP ')
    call ovarrf(outfile,'Pfirsch-Schlueter fraction','(psipf)',psipf, 'OP ')
    call ovarrf(outfile,'Auxiliary current drive fraction','(faccd)',faccd, 'OP ')
    call ovarrf(outfile,'Inductive fraction','(facoh)',facoh, 'OP ')
    ! Add total error check.
    call ovarrf(outfile,'Total','(plasipf+faccd+facoh)',plasipf+faccd+facoh)
    if (abs(plasipf+faccd+facoh-1.0d0) > 1.0d-8) then
        call ocmmnt(outfile,'ERROR: current drive fractions do not add to 1')
    end if
    ! MDK Add fvsbrnni as it can be an iteration variable
    call ovarrf(outfile,'Fraction of the plasma current produced by non-inductive means','(fvsbrnni)',fvsbrnni)

    if (abs(bootipf-bscfmax) < 1.0D-8) then
       call ocmmnt(outfile,'Warning : bootstrap current fraction is at')
       call ocmmnt(outfile,'          its prescribed maximum.')
    end if

    call oblnkl(outfile)

    if (abs(plhybd) > 1.0D-8) then
       !if (iefrffix.NE.0) then ! needs updating
       !  call ovarre(outfile,'Secondary RF efficiency (A/W)','(effrfssfix)',effrfssfix, 'OP ')
       !end if
       call ovarre(outfile,'RF efficiency (A/W)','(effrfss)',effrfss, 'OP ')
       call ovarre(outfile,'RF gamma (10^20 A/W-m2)','(gamrf)',gamrf, 'OP ')
       call ovarre(outfile,'Lower hybrid injected power (MW)','(plhybd)',plhybd, 'OP ')
       call ovarre(outfile,'Lower hybrid wall plug efficiency','(etalh)',etalh)
       call ovarre(outfile,'Lower hybrid wall plug power (MW)','(pwplh)',pwplh, 'OP ')
    end if

    ! MDK rearranged and added nbshinemw
    !if (abs(pnbeam) > 1.0D-8) then
    if ((iefrf == 5).or.(iefrf== 8).or.(iefrffix == 5).or.(iefrffix == 8)) then
       call ovarre(outfile,'Neutral beam energy (keV)','(enbeam)',enbeam)
       if ((iefrf == 5).or.(iefrf == 8)) then
         call ovarre(outfile,'Neutral beam current (A)','(cnbeam)',cnbeam, 'OP ')
       end if
       if ((iefrffix == 5).or.(iefrffix == 8)) then
         call ovarre(outfile,'Secondary fixed neutral beam current (A)','(cnbeamfix)',cnbeamfix, 'OP ')
       end if
       if ((iefrf == 5).or.(iefrf == 8)) then
         call ovarre(outfile,'Beam efficiency (A/W)','(effnbss)',effnbss, 'OP ')
       end if
       if ((iefrffix == 5).or.(iefrffix == 8)) then
         call ovarre(outfile,'Secondary fixed beam efficiency (A/W)','(effnbssfix)',effnbssfix, 'OP ')
       end if
       call ovarre(outfile,'Beam gamma (10^20 A/W-m2)','(gamnb)',gamnb, 'OP ')
       call ovarre(outfile,'Neutral beam wall plug efficiency','(etanbi)',etanbi)
       call ovarre(outfile,'Beam decay lengths to centre','(taubeam)',taubeam, 'OP ')
       call ovarre(outfile,'Beam shine-through fraction','(nbshinef)',nbshinef, 'OP ')
       call ovarre(outfile,'Neutral beam wall plug power (MW)','(pwpnb)',pwpnb, 'OP ')

       call oblnkl(outfile)
       call ocmmnt(outfile,'Neutral beam power balance :')
       call ocmmnt(outfile,'----------------------------')
       if ((iefrf == 5).or.(iefrf == 8)) then
         call ovarrf(outfile,'Beam first orbit loss power (MW)','(porbitlossmw)', porbitlossmw, 'OP ')
         call ovarrf(outfile,'Beam shine-through power [MW]','(nbshinemw)',nbshinemw, 'OP ')
         call ovarrf(outfile,'Beam power deposited in plasma (MW)','(pinjmw)',pinjmw1, 'OP ')
         call ovarrf(outfile,'Maximum allowable beam power (MW)','(pinjalw)',pinjalw)
         call ovarrf(outfile,'Total (MW)', &
                           '(porbitlossmw+nbshinemw+pinjmw)',porbitlossmw+nbshinemw+pinjmw1)
         call oblnkl(outfile)
         call ovarrf(outfile,'Beam power entering vacuum vessel (MW)','(pnbitot)',pnbitot, 'OP ')
       end if
       if ((iefrffix == 5).or.(iefrffix == 8)) then
         call oblnkl(outfile)
         call ocmmnt(outfile,'Secondary fixed neutral beam power balance :')
         call ocmmnt(outfile,'----------------------------')
         call ovarrf(outfile,'Secondary fixed beam first orbit loss power (MW)','(porbitlossmwfix)', porbitlossmwfix, 'OP ')
         call ovarrf(outfile,'Secondary fixed beam shine-through power [MW]','(nbshinemwfix)',nbshinemwfix, 'OP ')
         call ovarrf(outfile,'Secondary fixed beam power deposited in plasma (MW)','(pinjmwfix)',pinjmwfix, 'OP ')
         call ovarrf(outfile,'Maximum allowable beam power (MW)','(pinjalw)',pinjalw)
         call ovarrf(outfile,'Secondary fixed total (MW)', &
                           '(porbitlossmwfixed+nbshinemwfix+pinjmwfix)',porbitlossmwfix+nbshinemwfix+pinjmwfix)
         call oblnkl(outfile)
         call ovarrf(outfile,'Secondary beam power entering vacuum vessel (MW)','(pnbitotfix)',pnbitotfix, 'OP ')
       end if
       call oblnkl(outfile)

       call ovarre(outfile,'Fraction of beam energy to ions','(fpion)',fpion, 'OP ')
       call ovarre(outfile,'Beam duct shielding thickness (m)','(nbshield)',nbshield)
       call ovarre(outfile,'Beam tangency radius / Plasma major radius','(frbeam)',frbeam)
       call ovarre(outfile,'Beam centreline tangency radius (m)','(rtanbeam)', rtanbeam, 'OP ')
       call ovarre(outfile,'Maximum possible tangency radius (m)','(rtanmax)', rtanmax, 'OP ')
    end if

    if (abs(echpwr) > 1.0D-8) then
       call ovarre(outfile,'Electron cyclotron injected power (MW)','(echpwr)',echpwr, 'OP ')
       call ovarrf(outfile,'Maximum allowable ECRH power (MW)','(pinjalw)',pinjalw)
       call ovarre(outfile,'ECH wall plug efficiency','(etaech)',etaech)
       call ovarre(outfile,'ECH wall plug power (MW)','(echwpow)',echwpow, 'OP ')
    end if
    if (abs(pinjfixmw) > 1.0D-8) then
       call ovarrf(outfile,'Fixed ECRH power (MW)','(pinjmwfix)',pinjmwfix)
       call ovarre(outfile,'ECH wall plug efficiency','(etaech)',etaech)
       call ovarre(outfile,'Secondary fixed ECH wall plug power (MW)','(pinjwpfix)',pinjwpfix, 'OP ')
    end if

  end subroutine cudriv

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine iternb(effnbss,fpion,fshine)

    !! Routine to calculate ITER Neutral Beam current drive parameters
    !! author: P J Knight, CCFE, Culham Science Centre
    !! effnbss : output real : neutral beam current drive efficiency (A/W)
    !! fpion   : output real : fraction of NB power given to ions
    !! fshine  : output real : shine-through fraction of beam
    !! This routine calculates the current drive parameters for a
    !! neutral beam system, based on the 1990 ITER model.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
    !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use error_handling, only: fdiags, idiags, report_error
    use current_drive_variables, only: frbeam, enbeam, taubeam, ftritbm
    use physics_variables, only: eps, rmajor, abeam, te, dene, ralpne, rncne, &
      rnone, rnfene, deni, ten, zeffai, dlamie, alphan, alphat, aspect, zeff
		use constants, only: rmu0
    implicit none

    ! Arguments !
    ! !!!!!!!!!!!!

    real(dp), intent(out) :: effnbss,fpion,fshine

    ! Local variables !
    ! !!!!!!!!!!!!!!!!!!

    real(dp) :: dend,dent,dpath,sigstop

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Check argument sanity
    if ((1.0D0+eps) < frbeam) then
       fdiags(1) = eps ; fdiags(2) = frbeam
       call report_error(15)
    end if

    ! Calculate beam path length to centre
    dpath = rmajor * sqrt( (1.0D0 + eps)**2 - frbeam**2)

    ! Calculate beam stopping cross-section
    sigstop = sigbeam(enbeam/abeam,te,dene,ralpne,rncne,rnone,rnfene)

    ! Calculate number of decay lengths to centre
    taubeam = dpath * dene * sigstop

    ! Shine-through fraction of beam
    fshine = exp(-2.0D0 * dpath*dene*sigstop)
    fshine = max(fshine, 1.0D-20)

    ! Deuterium and tritium beam densities
    dend = deni * (1.0D0-ftritbm)
    dent = deni * ftritbm

    ! Power split to ions / electrons
    call cfnbi(abeam,enbeam,ten,dene,dend,dent,zeffai,dlamie,fpion)

    ! Current drive efficiency
    effnbss = frbeam * &
         etanb(abeam,alphan,alphat,aspect,dene,enbeam,rmajor,ten,zeff)

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function etanb(abeam,alphan,alphat,aspect,dene,ebeam,rmajor,ten,zeff)

      !! Routine to find neutral beam current drive efficiency
      !! using the ITER 1990 formulation
      !! author: P J Knight, CCFE, Culham Science Centre
      !! abeam   : input real : beam ion mass (amu)
      !! alphan  : input real : density profile factor
      !! alphat  : input real : temperature profile factor
      !! aspect  : input real : aspect ratio
      !! dene    : input real : volume averaged electron density (m**-3)
      !! enbeam  : input real : neutral beam energy (keV)
      !! rmajor  : input real : plasma major radius (m)
      !! ten     : input real : density weighted average electron temp. (keV)
      !! zeff    : input real : plasma effective charge
      !! This routine calculates the current drive efficiency of
      !! a neutral beam system, based on the 1990 ITER model.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
      !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: etanb

      ! Arguments !
      ! !!!!!!!!!!!!

      real(dp), intent(in) :: abeam,alphan,alphat,aspect,dene, &
           ebeam,rmajor,ten,zeff

      ! Local variables !
      ! !!!!!!!!!!!!!!!!!!

      real(dp) :: abd,bbd,dene20,dum,epseff,ffac,gfac,rjfunc, &
           xj,xjs,yj,zbeam

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! TODO comment this subroutine.

      zbeam = 1.0D0
      bbd = 1.0D0

      dene20 = 1.0D-20*dene

      ! Ratio of E_beam/E_crit
      xjs = ebeam / (bbd*10.0D0*abeam*ten)
      xj = sqrt(xjs)

      yj = 0.8D0 * zeff/abeam

      rjfunc = xjs / (4.0D0 + 3.0D0*yj + xjs * &
           (xj + 1.39D0 + 0.61D0*yj**0.7D0))

      epseff = 0.5D0/aspect
      gfac = (1.55D0 + 0.85D0/zeff)*sqrt(epseff) - &
           (0.2D0 + 1.55D0/zeff)*epseff
      ffac = 1.0D0/zbeam - (1.0D0 - gfac)/zeff

      abd = 0.107D0 * (1.0D0 - 0.35D0*alphan + 0.14D0*alphan**2) * &
           (1.0D0 - 0.21D0*alphat) * (1.0D0 - 0.2D-3*ebeam  + &
           0.09D-6 * ebeam**2)
      dum = abd *(5.0D0/rmajor) * (0.1D0*ten/dene20) * rjfunc/0.2D0 * ffac

      etanb = dum

    end function etanb

  end subroutine iternb

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cfnbi(afast,efast,te,ne,nd,nt,zeffai,xlmbda,fpion)

    !! Routine to calculate the fraction of the fast particle energy
    !! coupled to the ions
    !! author: P J Knight, CCFE, Culham Science Centre
    !! afast   : input real : mass of fast particle (units of proton mass)
    !! efast   : input real : energy of fast particle (keV)
    !! te      : input real : density weighted average electron temp. (keV)
    !! ne      : input real : volume averaged electron density (m**-3)
    !! nd      : input real : deuterium beam density (m**-3)
    !! nt      : input real : tritium beam density (m**-3)
    !! zeffai  : input real : mass weighted plasma effective charge
    !! xlmbda  : input real : ion-electron coulomb logarithm
    !! fpion   : output real : fraction of fast particle energy coupled to ions
    !! This routine calculates the fast particle energy coupled to
    !! the ions in the neutral beam system.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: mproton, pi, echarge

    implicit none

    !  Arguments

    real(dp), intent(in) :: afast,efast,te,ne,nd,nt,zeffai,xlmbda
    real(dp), intent(out) :: fpion

    !  Local variables

    real(dp) :: ans,ecritfi,ecritfix,sum,sumln,thx,t1,t2,ve,x, &
         xlbd,xlbt,xlmbdai,xlnrat

    real(dp), parameter :: atmd = 2.0D0
    real(dp), parameter :: atmdt = 2.5D0
    real(dp), parameter :: atmt = 3.0D0
    real(dp), parameter :: c = 3.0D8
    real(dp), parameter :: me = 9.1D-31
    real(dp), parameter :: zd = 1.0D0
    real(dp), parameter :: zt = 1.0D0

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xlbd = xlmbdabi(afast,atmd,efast,te,ne)
    xlbt = xlmbdabi(afast,atmt,efast,te,ne)

    sum = nd*zd*zd * xlbd/atmd + nt*zt*zt * xlbt/atmt
    ecritfix = 16.0D0 * te * afast * (sum/(ne*xlmbda))**(2.0D0/3.0D0)

    xlmbdai = xlmbdabi(afast,atmdt,efast,te,ne)
    sumln = zeffai * xlmbdai/xlmbda
    xlnrat = (3.0D0*sqrt(pi)/4.0D0 * me/mproton * sumln)**(2.0D0/3.0D0)
    ve = c * sqrt(2.0D0*te/511.0D0)

    ecritfi = afast * mproton * ve*ve * xlnrat/(2.0D0 * echarge * 1.0D3)

    x = sqrt(efast/ecritfi)
    t1 = log( (x*x - x + 1.0D0) / ((x + 1.0D0)**2) )
    thx = (2.0D0*x - 1.0D0)/sqrt(3.0D0)
    t2 = 2.0D0*sqrt(3.0D0) *(atan(thx) + pi/6.0D0)

    ans = (t1 + t2) / (3.0D0 * x*x)
    fpion = ans

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function xlmbdabi(mb,mth,eb,t,nelec)

      !! Calculates the Coulomb logarithm for ion-ion collisions
      !! author: P J Knight, CCFE, Culham Science Centre
      !! mb     : input real : mass of fast particle (units of proton mass)
      !! mth    : input real : mass of background ions (units of proton mass)
      !! eb     : input real : energy of fast particle (keV)
      !! t      : input real : density weighted average electron temp. (keV)
      !! nelec  : input real : volume averaged electron density (m**-3)
      !! This function calculates the Coulomb logarithm for ion-ion
      !! collisions where the relative velocity may be large compared
      !! with the background ('mt') thermal velocity.
      !! Mikkelson and Singer, Nuc Tech/Fus, 4, 237 (1983)
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: xlmbdabi

      !  Arguments

      real(dp), intent(in) :: mb,mth,eb,t,nelec

      !  Local variables

      real(dp) :: ans,x1,x2

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      x1 = (t/10.0D0) * (eb/1000.0D0) * mb/(nelec/1.0D20)
      x2 = mth/(mth + mb)

      ans = 23.7D0 + log(x2 * sqrt(x1))
      xlmbdabi = ans

    end function xlmbdabi

  end subroutine cfnbi

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function sigbeam(eb,te,ne,rnhe,rnc,rno,rnfe)

    !! Calculates the stopping cross-section for a hydrogen
    !! beam in a fusion plasma
    !! author: P J Knight, CCFE, Culham Science Centre
    !! eb     : input real : beam energy (kev/amu)
    !! te     : input real : electron temperature (keV)
    !! ne     : input real : electron density (10^20m-3)
    !! rnhe   : input real : alpha density / ne
    !! rnc    : input real : carbon density /ne
    !! rno    : input real : oxygen density /ne
    !! rnfe   : input real : iron density /ne
    !! This function calculates the stopping cross-section (m^-2)
    !! for a hydrogen beam in a fusion plasma.
    !! Janev, Boley and Post, Nuclear Fusion 29 (1989) 2125
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(dp) :: sigbeam

    !  Arguments

    real(dp), intent(in) :: eb,te,ne,rnhe,rnc,rno,rnfe

    !  Local variables

    real(dp) :: ans,nen,sz,s1
    real(dp), dimension(2,3,2) :: a
    real(dp), dimension(3,2,2,4) :: b
    real(dp), dimension(4) :: nn,z

    integer :: i,is,j,k

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    data a/ 4.40D0, 2.30D-1, 7.46D-2,-2.55D-3, 3.16D-3, 1.32D-3, &
         -2.49D-2,-1.15D-2, 2.27D-3,-6.20D-4,-2.78D-5, 3.38D-5/

    data b/ &
         -2.36D0, 8.49D-1,-5.88D-2,-2.50D-1, 6.77D-2,-4.48D-3, &
         1.85D-1,-4.78D-2, 4.34D-3,-3.81D-2, 1.05D-2,-6.76D-4, &
         -1.49D0, 5.18D-1,-3.36D-2,-1.19D-1, 2.92D-2,-1.79D-3, &
         -1.54D-2, 7.18D-3, 3.41D-4,-1.50D-2, 3.66D-3,-2.04D-4, &
         -1.41D0, 4.77D-1,-3.05D-2,-1.08D-1, 2.59D-2,-1.57D-3, &
         -4.08D-4, 1.57D-3, 7.35D-4,-1.38D-2, 3.33D-3,-1.86D-4, &
         -1.03D0, 3.22D-1,-1.87D-2,-5.58D-2, 1.24D-2,-7.43D-4, &
         1.06D-1,-3.75D-2, 3.53D-3,-3.72D-3, 8.61D-4,-5.12D-5/

    z(1) = 2.0D0 ; z(2) = 6.0D0 ; z(3) = 8.0D0 ; z(4) = 26.0D0
    nn(1) = rnhe ; nn(2) = rnc ; nn(3) = rno ; nn(4) = rnfe

    nen = ne/1.0D19

    s1 = 0.0D0
    do k = 1,2
       do j = 1,3
          do i = 1,2
             s1 = s1 + a(i,j,k) * (log(eb))**(i-1) * &
                  (log(nen))**(j-1) * (log(te))**(k-1)
          end do
       end do
    end do

    !  Impurity term

    sz = 0.0D0
    do is = 1,4
       do k = 1,2
          do j = 1,2
             do i = 1,3
                sz = sz + b(i,j,k,is)* (log(eb))**(i-1) * &
                     (log(nen))**(j-1) * (log(te))**(k-1) * &
                     nn(is) * z(is) * (z(is)-1.0D0)
             end do
          end do
       end do
    end do

    ans = 1.0D-20 * ( exp(s1)/eb * (1.0D0 + sz) )

    sigbeam = max(ans,1.0D-23)

  end function sigbeam

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cullhy(effrfss)

    !! Routine to calculate Lower Hybrid current drive efficiency
    !! author: P J Knight, CCFE, Culham Science Centre
    !! effrfss : output real : lower hybrid current drive efficiency (A/W)
    !! This routine calculates the current drive parameters for a
    !! lower hybrid system, based on the AEA FUS 172 model.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !! AEA FUS 172: Physics Assessment for the European Reactor Study
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use error_handling, only: fdiags, idiags, report_error
    use profiles_module, only: nprofile, tprofile
    use physics_variables, only: rminor, rhopedn, ne0, neped, nesep, alphan, &
      rhopedt, te0, tesep, teped, alphat, tbeta, bt, rmajor, zeff

    implicit none

    !  Arguments

    real(dp), intent(out) :: effrfss

    !  Local variables

    real(dp) :: blocal,dlocal,epslh,frac,gamlh,nplacc,rpenet, &
         rratio,term01,term02,term03,term04,tlocal,x

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Calculate the penetration radius of the LH waves

    call lhrad(rratio)
    rpenet = rratio*rminor

    !  Local density, temperature, toroidal field at this minor radius

    dlocal = 1.0D-19 * nprofile(rratio,rhopedn,ne0,neped,nesep,alphan)
    tlocal = tprofile(rratio,rhopedt,te0,teped,tesep,alphat,tbeta)
    blocal = bt*rmajor/(rmajor-rpenet)  !  Calculated on inboard side

    !  Parallel refractive index needed for plasma access

    frac = sqrt(dlocal)/blocal
    nplacc = frac + sqrt(1.0D0 + frac*frac)

    !  Local inverse aspect ratio

    epslh = rpenet/rmajor

    !  LH normalised efficiency (A/W m**-2)

    x = 24.0D0 / (nplacc*sqrt(tlocal))

    term01 = 6.1D0 / (nplacc*nplacc * (zeff+5.0D0))
    term02 = 1.0D0 + (tlocal/25.0D0)**1.16D0
    term03 = epslh**0.77D0 * sqrt(12.25D0 + x*x)
    term04 = 3.5D0*epslh**0.77D0 + x

    if (term03 > term04) then
       fdiags(1) = term03 ; fdiags(2) = term04
       call report_error(129)
    end if

    gamlh = term01 * term02 * (1.0D0 - term03/term04)

    !  Current drive efficiency (A/W)

    effrfss = gamlh / ((0.1D0*dlocal)*rmajor)

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lhrad(rratio)

      !! Routine to calculate Lower Hybrid wave absorption radius
      !! author: P J Knight, CCFE, Culham Science Centre
      !! rratio  : output real : minor radius of penetration / rminor
      !! This routine determines numerically the minor radius at which the
      !! damping of Lower Hybrid waves occurs, using a Newton-Raphson method.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !! AEA FUS 172: Physics Assessment for the European Reactor Study
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !  Arguments

      real(dp), intent(out) :: rratio

      !  Local variables

      real(dp) :: dgdr,drfind,g0,g1,g2,rat0,rat1,r1,r2
      integer :: lapno
      integer, parameter :: maxlap = 100

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Correction to refractive index (kept within valid bounds)

      drfind = min(0.7D0, max(0.1D0,12.5D0/te0))

      !  Use Newton-Raphson method to establish the correct minor radius
      !  ratio. g is calculated as a function of r / r_minor, where g is
      !  the difference between the results of the two formulae for the
      !  energy E given in AEA FUS 172, p.58. The required minor radius
      !  ratio has been found when g is sufficiently close to zero.

      !  Initial guess for the minor radius ratio

      rat0 = 0.8D0

      lapno = 0
      do ; lapno = lapno+1

         !  Minor radius ratios either side of the latest guess

         r1 = rat0 - 1.0D-3*rat0
         r2 = rat0 + 1.0D-3*rat0

         !  Evaluate g at rat0, r1, r2

         call lheval(drfind,rat0,g0)
         call lheval(drfind,r1,g1)
         call lheval(drfind,r2,g2)

         !  Calculate gradient of g with respect to minor radius ratio

         dgdr = (g2-g1)/(r2-r1)

         !  New approximation

         rat1 = rat0 - g0/dgdr

         !  Force this approximation to lie within bounds

         rat1 = max(0.0001D0,rat1)
         rat1 = min(0.9999D0,rat1)

         !  Check the number of laps for convergence

         if (lapno >= maxlap) then
            idiags(1) = lapno ; call report_error(16)
            rat0 = 0.8D0
            exit
         end if

         !  Is g sufficiently close to zero?

         if (abs(g0) > 0.01D0) then
            !  No, so go around loop again
            rat0 = rat1
         else
            exit
         end if

      end do

      rratio = rat0

    end subroutine lhrad

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lheval(drfind,rratio,ediff)

      !! Routine to evaluate the difference between electron energy
      !! expressions required to find the Lower Hybrid absorption radius
      !! author: P J Knight, CCFE, Culham Science Centre
      !! drfind  : input real : correction to parallel refractive index
      !! rratio  : input real : guess for radius of penetration / rminor
      !! ediff   : output real : difference between the E values (keV)
      !! This routine evaluates the difference between the values calculated
      !! from the two equations for the electron energy E, given in
      !! AEA FUS 172, p.58. This difference is used to locate the Lower Hybrid
      !! wave absorption radius via a Newton-Raphson method, in calling
      !! routine <A HREF="lhrad.html">lhrad</A>.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !! AEA FUS 172: Physics Assessment for the European Reactor Study
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !  Arguments

      real(dp), intent(in) :: drfind,rratio
      real(dp), intent(out) :: ediff

      !  Local variables

      real(dp) :: blocal,dlocal,e1,e2,frac,nplacc,refind,tlocal

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Local electron density

      dlocal = 1.0D-19 * nprofile(rratio,rhopedn,ne0,neped,nesep,alphan)

      !  Local electron temperature

      tlocal = tprofile(rratio,rhopedt,te0,teped,tesep,alphat,tbeta)

      !  Local toroidal field (evaluated at the inboard region of the flux surface)

      blocal = bt * rmajor/(rmajor - rratio*rminor)

      !  Parallel refractive index needed for plasma access

      frac = sqrt(dlocal)/blocal
      nplacc = frac + sqrt(1.0D0 + frac*frac)

      !  Total parallel refractive index

      refind = nplacc + drfind

      !  First equation for electron energy E

      e1 = 511.0D0 * (sqrt(1.0D0 + 1.0D0/(refind*refind)) - 1.0D0)

      !  Second equation for E

      e2 = 7.0D0 * tlocal

      !  Difference

      ediff = e1-e2

    end subroutine lheval

  end subroutine cullhy

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine culecd(effrfss)

    !! Routine to calculate Electron Cyclotron current drive efficiency
    !! author: M R O'Brien, CCFE, Culham Science Centre
    !! author: P J Knight, CCFE, Culham Science Centre
    !! effrfss : output real : electron cyclotron current drive efficiency (A/W)
    !! This routine calculates the current drive parameters for a
    !! electron cyclotron system, based on the AEA FUS 172 model.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !! AEA FUS 172: Physics Assessment for the European Reactor Study
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use profiles_module, only: tprofile, nprofile
    use physics_variables, only: rhopedt, te0, teped, tesep, alphat, tbeta, &
      rhopedn, ne0, nesep, neped, alphan, rminor, rmajor, zeff

    implicit none

    real(dp), intent(out) :: effrfss

    !  Local variables

    real(dp) :: cosang,coulog,dlocal,ecgam,ecgam1,ecgam2,ecgam3,ecgam4, &
         epsloc,rrr,tlocal,zlocal

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Local plasma parameters : take r = a/3
    rrr = 1.0D0/3.0D0

    !  Temperature
    tlocal = tprofile(rrr,rhopedt,te0,teped,tesep,alphat,tbeta)

    !  Density (10**20 m**-3)
    dlocal = 1.0D-20 * nprofile(rrr,rhopedn,ne0,neped,nesep,alphan)

    !  Inverse aspect ratio
    epsloc = rrr * rminor/rmajor

    !  Effective charge (use average value)
    zlocal = zeff

    !  Coulomb logarithm for ion-electron collisions
    !  (From J. A. Wesson, 'Tokamaks', Clarendon Press, Oxford, p.293)
    coulog = 15.2D0 - 0.5D0*log(dlocal) + log(tlocal)

    !  Calculate normalised current drive efficiency at four different
    !  poloidal angles, and average.
    !  cosang = cosine of the poloidal angle at which ECCD takes place
    !         = +1 outside, -1 inside.
    cosang =  1.0D0 ; call eccdef(tlocal,epsloc,zlocal,cosang,coulog,ecgam1)
    cosang =  0.5D0 ; call eccdef(tlocal,epsloc,zlocal,cosang,coulog,ecgam2)
    cosang = -0.5D0 ; call eccdef(tlocal,epsloc,zlocal,cosang,coulog,ecgam3)
    cosang = -1.0D0 ; call eccdef(tlocal,epsloc,zlocal,cosang,coulog,ecgam4)

    !  Normalised current drive efficiency (A/W m**-2)
    ecgam = 0.25D0 * (ecgam1+ecgam2+ecgam3+ecgam4)

    !  Current drive efficiency (A/W)
    effrfss = ecgam/(dlocal*rmajor)

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine eccdef(tlocal,epsloc,zlocal,cosang,coulog,ecgam)

      !! Routine to calculate Electron Cyclotron current drive efficiency
      !! author: M R O'Brien, CCFE, Culham Science Centre
      !! author: P J Knight, CCFE, Culham Science Centre
      !! tlocal : input real : local electron temperature (keV)
      !! epsloc : input real : local inverse aspect ratio
      !! zlocal : input real : local plasma effective charge
      !! cosang : input real : cosine of the poloidal angle at which ECCD takes
      !! place (+1 outside, -1 inside)
      !! coulog : input real : local coulomb logarithm for ion-electron collisions
      !! ecgam  : output real : normalised current drive efficiency (A/W m**-2)
      !! This routine calculates the current drive parameters for a
      !! electron cyclotron system, based on the AEA FUS 172 model.
      !! It works out the ECCD efficiency using the formula due to Cohen
      !! quoted in the ITER Physics Design Guidelines : 1989
      !! (but including division by the Coulomb Logarithm omitted from
      !! IPDG89). We have assumed gamma**2-1 << 1, where gamma is the
      !! relativistic factor. The notation follows that in IPDG89.
      !! <P>The answer ECGAM is the normalised efficiency nIR/P with n the
      !! local density in 10**20 /m**3, I the driven current in MAmps,
      !! R the major radius in metres, and P the absorbed power in MWatts.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !! AEA FUS 172: Physics Assessment for the European Reactor Study
      !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
      !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use error_handling, only: report_error

      implicit none

      !  Arguments

      real(dp), intent(in) :: tlocal,epsloc,zlocal,cosang,coulog
      real(dp), intent(out) :: ecgam

      !  Local variables

      real(dp) :: f,facm,fp,h,hp,lam,lams,mcsq,palpha,palphap,palphaps, &
           palphas,y

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mcsq = 9.1095D-31 * 2.9979D8**2 /(1.0D3*1.6022D-19) !  keV
      f = 16.0D0 * (tlocal/mcsq)**2

      !  fp is the derivative of f with respect to gamma, the relativistic
      !  factor, taken equal to 1 + 2T/(m c**2)

      fp = 16.0D0 * tlocal/mcsq

      !  lam is IPDG89's lambda. LEGEND calculates the Legendre function of
      !  order alpha and argument lam, palpha, and its derivative, palphap.
      !  Here alpha satisfies alpha(alpha+1) = -8/(1+zlocal). alpha is of the
      !  form  (-1/2 + ix), with x a real number and i = sqrt(-1).

      lam = 1.0D0
      call legend(zlocal,lam,palpha,palphap)

      lams = sqrt(2.0D0*epsloc/(1.0D0+epsloc))
      call legend(zlocal,lams,palphas,palphaps)

      !  hp is the derivative of IPDG89's h function with respect to lam

      h = -4.0D0 * lam/(zlocal+5.0D0) * (1.0D0-lams*palpha/(lam*palphas))
      hp = -4.0D0 / (zlocal+5.0D0) * (1.0D0-lams*palphap/palphas)

      !  facm is IPDG89's momentum conserving factor

      facm = 1.5D0
      y = mcsq/(2.0D0*tlocal) * (1.0D0 + epsloc*cosang)

      !  We take the negative of the IPDG89 expression to get a positive
      !  number

      ecgam = -7.8D0 * facm * sqrt((1.0D0+epsloc)/(1.0D0-epsloc)) / coulog &
           * (h*fp - 0.5D0*y*f*hp)

      if (ecgam < 0.0D0) call report_error(17)

    end subroutine eccdef

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine legend(zlocal,arg,palpha,palphap)

      !! Routine to calculate Legendre function and its derivative
      !! author: M R O'Brien, CCFE, Culham Science Centre
      !! author: P J Knight, CCFE, Culham Science Centre
      !! zlocal  : input real : local plasma effective charge
      !! arg     : input real : argument of Legendre function
      !! palpha  : output real : value of Legendre function
      !! palphap : output real : derivative of Legendre function
      !! This routine calculates the Legendre function <CODE>palpha</CODE>
      !! of argument <CODE>arg</CODE> and order
      !! <CODE>alpha = -0.5 + i sqrt(xisq)</CODE>,
      !! and its derivative <CODE>palphap</CODE>.
      !! <P>This Legendre function is a conical function and we use the series
      !! in <CODE>xisq</CODE> given in Abramowitz and Stegun. The
      !! derivative is calculated from the derivative of this series.
      !! <P>The derivatives were checked by calculating <CODE>palpha</CODE> for
      !! neighbouring arguments. The calculation of <CODE>palpha</CODE> for zero
      !! argument was checked by comparison with the expression
      !! <CODE>palpha(0) = 1/sqrt(pi) * cos(pi*alpha/2) * gam1 / gam2</CODE>
      !! (Abramowitz and Stegun, eqn 8.6.1). Here <CODE>gam1</CODE> and
      !! <CODE>gam2</CODE> are the Gamma functions of arguments
      !! <CODE>0.5*(1+alpha)</CODE> and <CODE>0.5*(2+alpha)</CODE> respectively.
      !! Abramowitz and Stegun, equation 8.12.1
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use error_handling, only: fdiags, idiags, report_error

      implicit none

      real(dp), intent(in) :: zlocal,arg
      real(dp), intent(out) ::  palpha,palphap

      !  Local variables

      real(dp) :: arg2,pold,poldp,pterm,sinsq,term1,term2,xisq
      integer :: n

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Check for invalid argument

      if (abs(arg) > (1.0D0+1.0D-10)) then
         fdiags(1) = arg ; call report_error(18)
      end if

      arg2 = min(arg, (1.0D0-1.0D-10))
      sinsq = 0.5D0*(1.0D0-arg2)
      xisq = 0.25D0*(32.0D0*zlocal/(zlocal+1.0D0) - 1.0D0)
      palpha = 1.0D0
      pold = 1.0D0
      pterm = 1.0D0
      palphap = 0.0D0
      poldp = 0.0D0

      do n = 1,10001

         !  Check for convergence every 20 iterations

         if ((n > 1).and.(mod(n,20) == 1)) then
            term1 = 1.0D-10 * max(abs(pold),abs(palpha))
            term2 = 1.0D-10 * max(abs(poldp),abs(palphap))

            if ( (abs(pold-palpha) < term1) .and. &
                 (abs(poldp-palphap) < term2) ) return

            pold = palpha
            poldp = palphap
         end if

         pterm = pterm * (4.0D0*xisq+(2.0D0*n - 1.0D0)**2) / &
              (2.0D0*n)**2 * sinsq
         palpha = palpha + pterm
         palphap = palphap - n*pterm/(1.0D0-arg2)

      end do

      !  Code will only get this far if convergence has failed

      call report_error(19)

    end subroutine legend

  end subroutine culecd

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine culnbi(effnbss,fpion,fshine)

    !! Routine to calculate Neutral Beam current drive parameters
    !! author: P J Knight, CCFE, Culham Science Centre
    !! effnbss : output real : neutral beam current drive efficiency (A/W)
    !! fpion   : output real : fraction of NB power given to ions
    !! fshine  : output real : shine-through fraction of beam
    !! This routine calculates Neutral Beam current drive parameters
    !! using the corrections outlined in AEA FUS 172 to the ITER method.
    !! <P>The result cannot be guaranteed for devices with aspect ratios far
    !! from that of ITER (approx. 2.8).
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !! AEA FUS 172: Physics Assessment for the European Reactor Study
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use error_handling, only: fdiags, idiags, report_error
    use current_drive_variables, only: enbeam, frbeam, taubeam, ftritbm
    use physics_variables, only: eps, rmajor, abeam, te, dene, ralpne, rncne, &
      rnone, rnfene, dnla, deni, ten, zeffai, dlamie, alphan, alphat, aspect, &
      rminor, zeff

    implicit none

    real(dp), intent(out) :: effnbss,fpion,fshine

    !  Local variables

    real(dp) :: dend,dent,dpath,sigstop

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Check argument sanity

    if ((1.0D0+eps) < frbeam) then
       fdiags(1) = eps ; fdiags(2) = frbeam
       call report_error(20)
    end if

    !  Calculate beam path length to centre

    dpath = rmajor * sqrt( (1.0D0 + eps)**2 - frbeam**2)

    !  Calculate beam stopping cross-section

    sigstop = sigbeam(enbeam/abeam,te,dene,ralpne,rncne,rnone,rnfene)

    !  Calculate number of decay lengths to centre

    taubeam = dpath * dnla * sigstop

    !  Shine-through fraction of beam

    fshine = exp(-2.0D0 * dpath*dnla*sigstop)
    fshine = max(fshine, 1.0D-20)

    !  Deuterium and tritium beam densities

    dend = deni * (1.0D0-ftritbm)
    dent = deni * ftritbm

    !  Power split to ions / electrons

    call cfnbi(abeam,enbeam,ten,dene,dend,dent,zeffai,dlamie,fpion)

    !  Current drive efficiency

    effnbss = etanb2(abeam,alphan,alphat,aspect,dene,dnla,enbeam, &
         frbeam,fshine,rmajor,rminor,ten,zeff)

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function etanb2(abeam,alphan,alphat,aspect,dene,dnla,enbeam,frbeam, &
      fshine,rmajor,rminor,ten,zeff)
      !! Routine to find neutral beam current drive efficiency
      !! using the ITER 1990 formulation, plus correction terms
      !! outlined in Culham Report AEA FUS 172
      !! author: P J Knight, CCFE, Culham Science Centre
      !! abeam   : input real : beam ion mass (amu)
      !! alphan  : input real : density profile factor
      !! alphat  : input real : temperature profile factor
      !! aspect  : input real : aspect ratio
      !! dene    : input real : volume averaged electron density (m**-3)
      !! dnla    : input real : line averaged electron density (m**-3)
      !! enbeam  : input real : neutral beam energy (keV)
      !! frbeam  : input real : R_tangent / R_major for neutral beam injection
      !! fshine  : input real : shine-through fraction of beam
      !! rmajor  : input real : plasma major radius (m)
      !! rminor  : input real : plasma minor radius (m)
      !! ten     : input real : density weighted average electron temperature (keV)
      !! zeff    : input real : plasma effective charge
      !! This routine calculates the current drive efficiency in A/W of
      !! a neutral beam system, based on the 1990 ITER model,
      !! plus correction terms outlined in Culham Report AEA FUS 172.
      !! <P>The formulae are from AEA FUS 172, unless denoted by IPDG89.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !! AEA FUS 172: Physics Assessment for the European Reactor Study
      !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
      !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: etanb2

      !  Arguments

      real(dp), intent(in) :: abeam,alphan,alphat,aspect,dene,dnla, &
           enbeam,frbeam,fshine,rmajor,rminor,ten,zeff

      !  Local variables

      real(dp) :: abd,bbd,d,dene20,dnla20,dnorm,ebmev,ebnorm, &
           ecrit,epseff,epsitr,eps1,ffac,gamnb,gfac,j0,nnorm,r,xj, &
           xjs,yj,zbeam

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Charge of beam ions
      zbeam = 1.0D0

      !  Fitting factor (IPDG89)
      bbd = 1.0D0

      !  Volume averaged electron density (10**20 m**-3)
      dene20 = dene/1.0D20

      !  Line averaged electron density (10**20 m**-3)
      dnla20 = dnla/1.0D20

      !  Critical energy (MeV) (power to electrons = power to ions) (IPDG89)
      !  N.B. ten is in keV
      ecrit = 0.01D0 * abeam * ten

      !  Beam energy in MeV
      ebmev = enbeam/1.0D3

      !  x and y coefficients of function J0(x,y) (IPDG89)
      xjs = ebmev/(bbd*ecrit)
      xj = sqrt(xjs)

      yj = 0.8D0 * zeff/abeam

      !  Fitting function J0(x,y)
      j0 = xjs / (4.0D0 + 3.0D0*yj + xjs *(xj + 1.39D0 + &
           0.61D0 * yj**0.7D0))

      !  Effective inverse aspect ratio, with a limit on its maximum value
      epseff = min(0.2D0, (0.5D0/aspect))

      !  Reduction in the reverse electron current
      !  due to neoclassical effects
      gfac = (1.55D0 + 0.85D0/zeff)*sqrt(epseff) - &
           (0.2D0 + 1.55D0/zeff)*epseff

      !  Reduction in the net beam driven current
      !  due to the reverse electron current
      ffac = 1.0D0 - (zbeam/zeff) * (1.0D0 - gfac)

      !  Normalisation to allow results to be valid for
      !  non-ITER plasma size and density:

      !  Line averaged electron density (10**20 m**-3) normalised to ITER
      nnorm = 1.0D0

      !  Distance along beam to plasma centre
      r = max(rmajor,rmajor*frbeam)
      eps1 = rminor/r

      if ((1.0D0+eps1) < frbeam) then
         fdiags(1) = eps1 ; fdiags(2) = frbeam
         call report_error(21)
      end if

      d = rmajor * sqrt( (1.0D0+eps1)**2 - frbeam**2)

      !  Distance along beam to plasma centre for ITER
      !  assuming a tangency radius equal to the major radius
      epsitr = 2.15D0/6.0D0
      dnorm = 6.0D0 * sqrt(2.0D0*epsitr + epsitr**2)

      !  Normalisation to beam energy (assumes a simplified formula for
      !  the beam stopping cross-section)
      ebnorm = ebmev * ( (nnorm*dnorm)/(dnla20*d) )**(1.0D0/0.78D0)

      !  A_bd fitting coefficient, after normalisation with ebnorm
      abd = 0.107D0 * (1.0D0 - 0.35D0*alphan + 0.14D0*alphan**2) * &
           (1.0D0 - 0.21D0*alphat) * (1.0D0 - 0.2D0*ebnorm + 0.09D0*ebnorm**2)

      !  Normalised current drive efficiency (A/W m**-2) (IPDG89)
      gamnb = 5.0D0 * abd * 0.1D0*ten * (1.0D0-fshine) * frbeam * &
           j0/0.2D0 * ffac

      !  Current drive efficiency (A/W)
      etanb2 = gamnb / (dene20*rmajor)

    end function etanb2

  end subroutine culnbi

end module current_drive_module
