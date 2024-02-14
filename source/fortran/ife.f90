! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ife_module
  !! Module containing Inertial Fusion Energy device routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains routines for calculating the
  !! parameters of an Inertial Fusion Energy power plant.
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifephy(outfile,iprint)

    !! Routine to calculate the physics parameters of an Inertial Fusion
    !! Energy power plant
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer : output file unit
    !! iprint : input integer : switch for writing to output file (1=yes)
    !! This routine calculates the physics parameters of an Inertial Fusion
    !! Energy power plant.
    !! F/MI/PJK/LOGBOOK12, pp.68,85
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: pi
    use build_variables, only: fwarea
    use error_handling, only: idiags, report_error
    use process_output, only: oheadr, oblnkl, ocmmnt, ovarre
    use ife_variables, only: ifedrv, edrive, gainve, etave, gain, etadrv, tgain, &
      drveff, reprat, pdrive, rrin, pfusife, ifetyp, zl1, r1, zu1, flirad
    use physics_variables, only: powfmw, wallmw

		use constants, only: pi
    implicit none

    !  Arguments
    integer, intent(in) :: outfile,iprint

    !  Local variables

    real(dp) :: aaion,bmax,dpp,dtheta,emitt,etai,lf,phi,qion, &
         sang,sigma,sigma0,tauf,theta,vi
    integer :: nbeams

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Driver calculations

    select case (ifedrv)

    case (-1)
       !  Target gain and driver efficiency dependencies on
       !  driver energy are input

       call driver(edrive,gainve,etave,gain,etadrv)

    case (0)  !  Target gain and driver efficiency are input

       gain = tgain
       etadrv = drveff

    case (1)  !  Laser driver based on SOMBRERO design

       call lasdrv(edrive,gain,etadrv)

    case (2)  !  Heavy-ion beam driver based on OSIRIS design

       aaion = 131.0D0
       bmax = 10.0D0
       dpp = 1.25D-3
       dtheta = 1.3D-4
       emitt = 1.0D-5
       etai = 0.8D0
       lf = 5.0D0
       nbeams = 12
       qion = 1.0D0
       sigma = 8.0D0
       sigma0 = 80.0D0
       tauf = 1.0D-7
       theta = 30.0D-3
       vi = 3.0D6

       call iondrv(aaion,bmax,dpp,dtheta,edrive,emitt,etai,lf, &
            nbeams,qion,sigma,sigma0,tauf,theta,vi,gain,etadrv)

    case(3)
       etadrv = drveff

    case default
       idiags(1) = ifedrv
       call report_error(127)

    end select

    if (ifedrv /= 3) then
        !  Repetition rate (Hz)
        reprat = pdrive / edrive
        !  Fusion power (MW)
        powfmw = 1.0D-6 * pdrive * gain

    else
        !  Driver Power
        reprat = rrin
        pdrive = reprat * edrive
        !  Gain
        powfmw = pfusife
        gain = powfmw / (1.0D-6 * pdrive)
    end if

    !  Wall load (assume total fusion power applies)

    if (ifetyp == 1) then

       !  OSIRIS-type build: First wall subtends a solid angle of 2 pi * SANG

       phi = 0.5D0*pi + atan(zl1/r1)
       sang = 1.0D0 - cos(phi)
       wallmw = powfmw * 0.5D0*sang / fwarea

    else if (ifetyp == 4) then

       ! 2019 build only has first wall at the top which has a tube at
       ! its centre.  This calculates solid angle and removes tube.

       phi = atan(r1/zu1)
       sang = 1.0D0 - cos(phi)
       phi = atan(flirad/zu1)
       sang = sang - (1.0D0 - cos(phi))
       wallmw = powfmw * 0.5D0*sang / fwarea

    else
       wallmw = powfmw / fwarea
    end if

    if (iprint == 0) return

    !  Output section

    call oheadr(outfile,'Physics / Driver Issues')

    select case (ifedrv)

    case (-1, 0)
       call ocmmnt(outfile,'Driver type : generic')
    case (1)
       call ocmmnt(outfile,'Driver type : laser')
    case (2)
       call ocmmnt(outfile,'Driver type : heavy ion beam')
    end select
    call oblnkl(outfile)

    call ovarre(outfile,'Driver energy (J)','(edrive)',edrive)
    call ovarre(outfile,'Driver efficiency','(etadrv)',etadrv)
    call ovarre(outfile,'Driver power reaching target (W)','(pdrive)', &
         pdrive)
    call ovarre(outfile,'Driver repetition rate (Hz)','(reprat)',reprat)
    call ovarre(outfile,'Target gain','(gain)',gain)
    call ovarre(outfile,'Fusion power (MW)','(powfmw)',powfmw)
    call ovarre(outfile,'Neutron wall load (MW/m2)','(wallmw)',wallmw)

  end subroutine ifephy

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine iondrv(aaion,bmax,dpp,dtheta,edrive,emitt,etai,lf, &
       nbeams,qion,sigma,sigma0,tauf,theta,vi,gain,etadrv)

    !! Routine to calculate parameters of a heavy ion driver
    !! suitable for inertial fusion energy
    !! author: P J Knight, CCFE, Culham Science Centre
    !! aaion  : input real : Ion mass (amu)
    !! bmax   : input real : Maximum field at the superconductor (T)
    !! dpp    : input real : Beam momentum spread
    !! dtheta : input real : Pointing error (rad)
    !! edrive : input real : Driver energy (J)
    !! emitt  : input real : Normalised beam emittance (metre-rad)
    !! etai   : input real : Axial quadrupole packing fraction at injection
    !! lf     : input real : Distance from final focussing quad to target (m)
    !! nbeams : input integer : Number of beams
    !! qion   : input real : Ion charge state
    !! sigma  : input real : Depressed tune (incl. space charge effects) (deg)
    !! sigma0 : input real : Phase advance per lattice period (tune) (deg)
    !! tauf   : input real : Post-acceleration pulse duration (s)
    !! theta  : input real : Final focussing half angle (rad)
    !! vi     : input real : Injection voltage (V)
    !! gain   : output real : Target gain
    !! etadrv : output real : Driver efficiency
    !! This routine calculates the parameters of a heavy ion driver
    !! suitable for inertial fusion energy.
    !! <P>Currently, the complicated model taken from the reference is not
    !! complete, so it is recommended that the simple model is used
    !! (set <CODE>ISIMP=1</CODE>)
    !! Heavy-ion Driver Design and Scaling, R. Bieri et al.,
    !! Fusion Technology, vol.21 (1992) 1583
    !! Meier and Bieri, Fusion Technology, vol.21 (1992) 1547
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: echarge, mproton, degrad, pi
    implicit none

    !  Arguments

    real(dp), intent(in) :: aaion,bmax,dpp,dtheta,edrive,emitt,etai,lf,qion, &
         sigma,sigma0,tauf,theta,vi
    integer, intent(in) :: nbeams
    real(dp), intent(out) :: gain,etadrv

    !  Local variables

    real(dp), parameter :: c2 = 8.98755178737D16
    integer, parameter :: isimp = 1  !  Switch for simple model (1=yes)

    real(dp) :: ci,de,dgap,dlcore,drcore,e,eomc2,fins,floss, &
         ibf,ibfo,ibi,ibpc,lfocus,lpf,lpfo,lpi,lppc,lq,phif,phifo, &
         phii,phipc,rion,rs,rs1,rs2,rs3,rs4,sig,sig0,taufo,taui, &
         taupc,tbrad,vf,vfo,vpc,vs,vscore,vshe,vsle,vspc,xhe,xle,xpc
    real(dp), dimension(10) :: eve,gve
    integer :: i,ie,j,nche,ncle,ncpc,nqarrs,nqle,nqpche,nquads

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (isimp /= 1) then

       !  Complex IONDRV model (incomplete)...
       !  ====================

       !  Electron charge / (proton mass * c**2)

       eomc2 = echarge / (mproton*c2)

       !  Degrees to radians

       sig = sigma * degrad
       sig0 = sigma0 * degrad

       !  Beam current coefficient

       ci = cbeam(etai,bmax,aaion,qion,sig,sig0,emitt,vi,eomc2)

       !  Final voltage (V)

       vf = (edrive/(dble(nbeams)*tauf*ci))**(0.666666D0)

       !  Final beam current (A)

       ibf = ci*sqrt(vf)

       !  Ion range (g/cm^2)

       rion = (3.04D-5 * aaion + 349.0D0*aaion**(-2.21D0)) * &
            (qion*vf/1.0D9)**(1.72D0-0.00275D0*aaion)

       !  Convert to kg/m^2

       rion = rion*10.0D0

       !  Total spot size (rms of four terms):
       !  (1) Emittance term

       rs1 = emitt / (betgam(aaion,qion,vf)*theta)

       !  (2) Dispersion (momentum spread) term

       rs2 = 8.0D0*dpp*lf*theta

       !  (3) Alignment and jitter term

       rs3 = dtheta*lf

       !  (4) Space charge term (inaccurate at present)

       rs4 = 1.0D-3

       rs = sqrt(rs1**2 + rs2**2 + rs3**2 + rs4**2)

       !  Target gain (inaccurate at present)

       gain = 62.6D0 - 142.3D0*rs**1.25D0 * rion**0.633D0 + &
            (32.8D0+82.1D0*rs)**0.783D0 * rion**0.647D0 * &
            log(1.0D-6*edrive)

       !  Source/Injector characteristics:
       !  Transport beam radius (m)

       tbrad = ( (emitt**2 * sig0 * c2)/(sig**2 * bmax/1.5D0 * etai * &
            sqrt(vi)) )**(0.333333D0) * sqrt(2.0D0 * aaion * mproton / &
            (qion * echarge))

       !  Extractor voltage (V)

       vs = 4.0D0/3.0D0 * vi

       !  Source extraction gap width (m)

       dgap = 0.8D0

       !  Injected current per beam (A)

       ibi = 4.0D0 * pi * tbrad**2 * 5.46D-8 * &
            sqrt( qion*vs**3/(aaion*dgap**4) )

       !  Beam loss fraction between injector and target

       floss = 0.08D0

       !  Pulse duration at injection (s)

       taui = 1.0D0/(1.0D0-floss) * tauf * ibf / ibi

       !  Pulse length at injection (m)

       lpi = taui * sqrt( 2.0D0 * qion * echarge * vi / (aaion * mproton) )

       !  Initial voltage gradient (V/m)

       phii = 0.3D0 * vi / lpi

       !  End of Low Energy Transport Stage:
       !  Insulated fraction of driver length

       fins = 0.85D0

       !  Maximum voltage gradient, limited by insulator flash-over (V/m)

       phifo = 1.0D6 * fins

       !  Voltage at which this limit is reached (V)
       !  The end of the low energy stage is set by the point at which
       !  the limit is reached.

       vfo = (phifo/phii)**(2.0D0/3.0D0) * vi

       !  Beam current (A)

       ibfo = ibi * (vfo/vi)

       !  Pulse duration (s)

       taufo = taui * (vi/vfo)

       !  Pulse length (m)

       lpfo = taui * vi / vfo * sqrt( 2.0D0 * qion * echarge * vfo / &
            (aaion * mproton) )

       !  Length of the low energy transport stage (m)
       !  (rearrangement of integral of PHI(V).dl)

       xle = 2.0D0*vi**1.5D0 / phii * (1.0D0/sqrt(vi) - 1.0D0/sqrt(vfo))

       !  End of Pulse Compression Stage:
       !  Pulse length = final pulse length (m)

       lpf = tauf * sqrt( 2.0D0 * qion * echarge * vf / (aaion * mproton) )
       lppc = lpf

       !  Length of the pulse compression region (m)

       xpc = (lpfo - lpf)/0.3D0

       !  Voltage gradient (V/m) - limit already reached

       phipc = phifo

       !  Voltage (V)

       vpc = vfo + (xpc*phipc)

       !  Beam current (A)

       ibpc = ibfo * (lpfo/lppc) * sqrt(vpc/vfo)

       !  Pulse duration (s)

       taupc = taufo * (lppc/lpfo) * sqrt(vfo/vpc)

       !  End of High Energy Transport Stage:
       !  Voltage gradient (V/m) - limit already reached

       phif = phifo

       !  Length of the high energy transport stage (m)

       xhe = (vf - vfo)/phifo

       write(*,*) '                 Injector:'
       write(*,*) '   vi = ',vi
       write(*,*) '  ibi = ',ibi
       write(*,*) ' taui = ',taui
       write(*,*) ' phii = ',phii
       write(*,*) '  lpi = ',lpi
       write(*,*) '                 End of LET stage:'
       write(*,*) '  vfo = ',vfo
       write(*,*) ' ibfo = ',ibfo
       write(*,*) 'taufo = ',taufo
       write(*,*) 'phifo = ',phifo
       write(*,*) ' lpfo = ',lpfo
       write(*,*) '  xle = ',xle
       write(*,*) '                 End of PC stage:'
       write(*,*) '  vpc = ',vpc
       write(*,*) ' ibpc = ',ibpc
       write(*,*) 'taupc = ',taupc
       write(*,*) 'phipc = ',phipc
       write(*,*) ' lppc = ',lppc
       write(*,*) '  xpc = ',xpc
       write(*,*) '                 End of HET stage:'
       write(*,*) '   vf = ',vf
       write(*,*) '  ibf = ',ibf
       write(*,*) ' tauf = ',tauf
       write(*,*) ' phif = ',phif
       write(*,*) '  lpf = ',lpf
       write(*,*) '  xhe = ',xhe

       !  Volt-second requirements and number of cores

       !  Volt-seconds per core (2.4T is field swing of the core
       !  material metglas)

       drcore = 0.8D0
       dlcore = 0.2D0
       vscore = 2.4D0 * drcore * dlcore

       !  Low energy stage

       vsle = taui * vi * log(vfo/vi)
       ncle = 1 + int(vsle/vscore)

       !  Pulse compression stage

       vspc = 2.0D0 * tauf * sqrt(vf) * (sqrt(vpc)-sqrt(vfo))
       ncpc = 1 + int(vspc/vscore)

       !  High energy stage

       drcore = 0.4D0
       dlcore = 0.1D0
       vscore = 2.4D0 * DRCORE * DLCORE

       vshe = 2.0D0 * tauf * vf * (1.0D0 - sqrt(vpc/vf))
       nche = 1 + int(vshe/vscore)

       write(*,*) '                 Number of cores: '
       write(*,*) ' ncle = ',ncle
       write(*,*) ' ncpc = ',ncpc
       write(*,*) ' nche = ',nche

       !  Quadrupole requirements

       !  Effective quadrupole length (m)

       lq = ( emitt * etai * sig0**2 * sqrt(c2*vi) / &
            (sig*(bmax/1.5D0)**2) )**(0.333333D0) * &
            sqrt(2.0D0*aaion*mproton/(qion*echarge))

       write(*,*) '   lq = ',lq

       !  Actual quadrupole length (m)

       write(*,*) '  alq = ',lq/etai

       !  Length of quadrupole focussing fields: low energy stage

       lfocus = etai * vi/phii * (1.0D0 - vi/vfo)
       nqle = 1 + int(lfocus/lq)

       !  Pulse compression + high energy stages

       lfocus = 2.0D0 * etai * sqrt(vi)/phifo * (sqrt(vf) - sqrt(vfo))

       nqpche = 1 + int(lfocus/lq)

       !  Number of quadrupole arrays

       nqarrs = nqle + nqpche

       !  Number of quadrupoles

       nquads = nqarrs * nbeams

       write(*,*) 'nqarrs = ',nqarrs
       write(*,*) 'nquads = ',nquads

       etadrv = 0.28D0

    else

       !  Simple model
       !  ============

       !  gve(k): target gain at edrive = k MegaJoules

       gve(1)  = 25.0D0
       gve(2)  = 44.0D0
       gve(3)  = 62.0D0
       gve(4)  = 76.0D0
       gve(5)  = 87.0D0
       gve(6)  = 97.0D0
       gve(7)  = 107.0D0
       gve(8)  = 115.0D0
       gve(9)  = 125.0D0
       gve(10) = 132.0D0

       !  eve(k): driver efficiency at edrive = k MegaJoules

       eve(1)  = 0.232D0
       eve(2)  = 0.256D0
       eve(3)  = 0.269D0
       eve(4)  = 0.276D0
       eve(5)  = 0.282D0
       eve(6)  = 0.286D0
       eve(7)  = 0.290D0
       eve(8)  = 0.292D0
       eve(9)  = 0.294D0
       eve(10) = 0.296D0

       e = 1.0D-6 * edrive
       ie = int(e)
       de = e - dble(IE)

       !  Assume linear interpolations and extrapolations
       !  Would be better to prevent extrapolation

       if (ie <= 1) then

          gain   = gve(2) - 1.0D-6*(edrive-2.0D6)*(gve(1)-gve(2))
          etadrv = eve(2) - 1.0D-6*(edrive-2.0d6)*(eve(1)-eve(2))

       else if (ie >= 9) then

          gain   = gve(9) + 1.0D-6*(edrive-9.0D6)*(gve(10)-gve(9))
          etadrv = eve(9) + 1.0D-6*(edrive-9.0D6)*(eve(10)-eve(9))

       else

          gain   = gve(ie) + de*(gve(ie+1)-gve(ie))
          etadrv = eve(ie) + de*(eve(ie+1)-eve(ie))

       end if

       !  Ensure sensible values

       gain = max(0.01D0,gain)
       etadrv = max(0.01D0,etadrv)

    end if

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function cbeam(etai,bmax,aaion,qion,sigma,sigma0,emitt,vi,eomc2)

      !! Routine to evaluate the beam current coefficient, CI
      !! author: P J Knight, CCFE, Culham Science Centre
      !! etai   : input real : Axial quadrupole packing fraction at injection
      !! bmax   : input real : Maximum field at the superconductor (T)
      !! aaion  : input real : Ion mass (amu)
      !! qion   : input real : Ion charge state
      !! sigma  : input real : Depressed tune (incl. space charge effects) (rad)
      !! sigma0 : input real : Phase advance per lattice period (tune) (rad)
      !! emitt  : input real : Normalised beam emittance (metre-rad)
      !! vi     : input real : Injection voltage (V)
      !! eomc2  : input real : Electron charge / (proton mass * light speed**2)
      !! This routine calculates the beam current coefficient
      !! of a heavy ion driver suitable for inertial fusion energy.
      !! Heavy-ion Driver Design and Scaling, R. Bieri et al.,
      !! Fusion Technology, vol.21 (1992) 1583
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: cbeam

      !  Arguments

      real(dp), intent(in) :: etai,bmax,aaion,qion,sigma,sigma0, &
           emitt,vi,eomc2

      !  Local variables

      real(dp) :: bbe

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Maximum field at the beam edge

      bbe = bmax/1.5D0

      !  Perform calculation

      cbeam = 2.89D6*(1.0D0-(sigma/sigma0)**2) * &
           ( sigma0**4 * etai**2 * vi * aaion/qion * (emitt/sigma)**2 &
           * bbe**2 )**0.333333D0 &
           * (2.0D0 * qion * eomc2 / aaion)**(5.0D0/6.0D0)

    end function cbeam

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function betgam(aaion,qion,v)

      !! Routine to calculate the relativistic factor (beta*gamma) for the
      !! heavy ions in the beam driver
      !! author: P J Knight, CCFE, Culham Science Centre
      !! aaion  : input real : Ion mass (amu)
      !! qion   : input real : Ion charge state
      !! v      : input real : Acceleration voltage (V)
      !! This routine calculates the relativistic factor (beta*gamma) for the
      !! heavy ions in a beam driver suitable for inertial fusion energy.
      !! Heavy-ion Driver Design and Scaling, R. Bieri et al.,
      !! Fusion Technology, vol.21 (1992) 1583
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: betgam

      !  Arguments

      real(dp), intent(in) :: aaion,qion,v

      !  Local variables

      real(dp), parameter :: c2 = 8.98755178737D16

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      betgam = sqrt( 2.0D0*qion*echarge*v / (aaion*mproton*c2) )

    end function betgam

  end subroutine iondrv

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine lasdrv(edrive,gain,etadrv)

    !! Routine to calculate parameters of a laser driver
    !! suitable for inertial fusion energy
    !! author: P J Knight, CCFE, Culham Science Centre
    !! edrive : input real : Driver energy (J)
    !! gain   : output real : Target gain
    !! etadrv : output real : Driver efficiency
    !! This routine calculates the parameters of a laser driver
    !! suitable for inertial fusion energy.
    !! Gain and driver efficiency data are taken from Figures 1 and 2 of
    !! Meier and Rosenberg.
    !! Meier and Rosenberg, Fusion Technology vol.21 (1992) p.1552
    !! F/MI/PJK/LOGBOOK12, p.86
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    real(dp), intent(in) :: edrive
    real(dp), intent(out) :: etadrv,gain

    !  Local variables

    real(dp) :: e,de
    real(dp), dimension(10) :: gve,eve
    integer :: ie

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  GVE(K): target gain at EDRIVE = K MegaJoules

    gve(1)  = 63.0D0
    gve(2)  = 95.0D0
    gve(3)  = 112.0D0
    gve(4)  = 125.0D0
    gve(5)  = 136.0D0
    gve(6)  = 144.0D0
    gve(7)  = 151.0D0
    gve(8)  = 157.0D0
    gve(9)  = 162.0D0
    gve(10) = 166.0D0

    !  EVE(K): driver efficiency at EDRIVE = K MegaJoules

    eve(1)  = 0.082D0
    eve(2)  = 0.079D0
    eve(3)  = 0.076D0
    eve(4)  = 0.072D0
    eve(5)  = 0.069D0
    eve(6)  = 0.064D0
    eve(7)  = 0.059D0
    eve(8)  = 0.054D0
    eve(9)  = 0.048D0
    eve(10) = 0.042D0

    e = 1.0D-6 * edrive
    ie = int(e)
    de = e - dble(ie)

    !  Assume linear interpolations and extrapolations
    !  Would be better to prevent extrapolation

    if (ie <= 1) then

       gain   = gve(2) - 1.0D-6*(edrive-2.0D6)*(gve(1)-gve(2))
       etadrv = eve(2) - 1.0D-6*(edrive-2.0D6)*(eve(1)-eve(2))

    else if (ie >= 9) then

       gain   = gve(9) + 1.0D-6*(edrive-9.0D6)*(gve(10)-gve(9))
       etadrv = eve(9) + 1.0D-6*(edrive-9.0D6)*(eve(10)-eve(9))

    else

       gain = gve(ie)   + de*(gve(ie+1)-gve(ie))
       etadrv = eve(ie) + de*(eve(ie+1)-eve(ie))

    end if

    !  Ensure sensible values

    gain = max(0.01D0,gain)
    etadrv = max(0.01D0,etadrv)

  end subroutine lasdrv

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine driver(edrive,gainve,etave,gain,etadrv)

    !! Routine to calculate parameters of a generic driver
    !! suitable for inertial fusion energy
    !! author: P J Knight, CCFE, Culham Science Centre
    !! edrive : input real : Driver energy (J)
    !! gainve(10) : input real array : Gain vs energy data
    !! etave(10) : input real array : Driver efficiency vs energy data
    !! gain   : output real : Target gain
    !! etadrv : output real : Driver efficiency
    !! This routine calculates the parameters of a generic driver
    !! suitable for inertial fusion energy.
    !! Gain and driver efficiency data are interpolated from input data.
    !! F/MI/PJK/LOGBOOK12, p.85
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    real(dp), intent(in) :: edrive
    real(dp), dimension(10), intent(in) :: etave,gainve
    real(dp), intent(out) :: etadrv,gain

    !  Local variables

    real(dp) :: de,e
    integer :: ie

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  The arrays contain data points for EDRIVE = 1MJ, 2MJ, ... , 10MJ

    e = 1.0D-6 * edrive
    ie = int(e)
    de = e - dble(ie)

    !  Assume linear interpolations and extrapolations

    if (ie <= 1) then

       gain =  gainve(2) - 1.0D-6*(edrive-2.0D6)*(gainve(1)-gainve(2))
       etadrv = etave(2) - 1.0D-6*(edrive-2.0D6)*(etave(1)-etave(2))

    else if (ie >= 9) then

       gain =  gainve(9) + 1.0D-6*(edrive-9.0D6)*(gainve(10)-gainve(9))
       etadrv = etave(9) + 1.0D-6*(edrive-9.0D6)*(etave(10)-etave(9))

    else

       gain =  gainve(ie) + de*(gainve(ie+1)-gainve(ie))
       etadrv = etave(ie) + de*(etave(ie+1)-etave(ie))

    end if

    !  Ensure sensible values

    gain = max(0.01D0,gain)
    etadrv = max(0.01D0,etadrv)

  end subroutine driver

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifebld(outfile,iprint)
    !! Routine to create the build of an inertial fusion energy device
    !! and to calculate the material volumes for the device core
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer : output file unit
    !! iprint : input integer : switch for writing to output file (1=yes)
    !! This routine constructs the build of an inertial fusion energy device
    !! and calculates the material volumes for the device core.
    !! F/MI/PJK/LOGBOOK12, p.52
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use process_output, only: oheadr, obuild, ovarre
    use ife_variables, only: ifetyp, chrad, r1, fwdr, r2, v1dr, r3, bldr, r4, &
      v2dr, r5, shdr, r6, v3dr, r7, zl7, v3dzl, zl6, shdzl, zl5, v2dzl, zl4, &
      bldzl, zl3, v1dzl, zl2, fwdzl, zl1, chdzl, chdzu, zu1, fwdzu, zu2, &
      v1dzu, zu3, bldzu, zu4, v2dzu, zu5, shdzu, zu6, v3dzu, zu7, fwmatv, &
      v1matv, blmatv, v2matv, shmatv, v3matv, chmatv
		use constants, only: mfile
    implicit none

    !  Arguments
    integer, intent(in) :: outfile,iprint

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    select case (ifetyp)
    case (1)  !  OSIRIS-type device
       call osibld
    case (2)  !  SOMBRERO-type device
       call sombld
    case (3)  !  HYLIFE-II-type device
       call hylbld
    case (4)  !  2019 device
       call bld2019
    case default  !  Generic device
       call genbld
    end select

    if (iprint == 0) return

    !  Output section

    call oheadr(outfile,'Radial Build')
    write(outfile,20)
    20  format(T43,'Thickness (m)',T60,'Radius (m)')

    call obuild(outfile,'Device centreline',0.0D0,0.0D0)
    call obuild(outfile,'Chamber',chrad,r1)
    call ovarre(mfile,'Chamber (m)','(chrad)',chrad)
    call obuild(outfile,'First Wall',fwdr,r2)
    call ovarre(mfile,'First Wall (m)','(fwdr)',fwdr)
    call obuild(outfile,'Void 1',v1dr,r3)
    call ovarre(mfile,'Void 1 (m)','(v1dr)',v1dr)
    call obuild(outfile,'Blanket',bldr,r4)
    call ovarre(mfile,'Blanket (m)','(bldr)',bldr)
    call obuild(outfile,'Void 2',v2dr,r5)
    call ovarre(mfile,'Void 2 (m)','(v2dr)',v2dr)
    call obuild(outfile,'Shield',shdr,r6)
    call ovarre(mfile,'Shield (m)','(shdr)',shdr)
    call obuild(outfile,'Void 3',v3dr,r7)
    call ovarre(mfile,'Void 3 (m)','(v3dr)',v3dr)

    30  format(T43,'Thickness (m)',T60,'Height (m)')

    if (ifetyp /= 4) then

        call oheadr(outfile,'Vertical Build')
        write(outfile,30)

        call obuild(outfile,'Base of device',0.0D0,-zl7)
        call obuild(outfile,'Void 3',v3dzl,-zl6)
        call obuild(outfile,'Shield',shdzl,-zl5)
        call obuild(outfile,'Void 2',v2dzl,-zl4)
        call obuild(outfile,'Blanket',bldzl,-zl3)
        call obuild(outfile,'Void 1',v1dzl,-zl2)
        call obuild(outfile,'First Wall',fwdzl,-zl1)
        call obuild(outfile,'Chamber',chdzl,0.0D0)
        call obuild(outfile,'Chamber',chdzu,zu1)
        call obuild(outfile,'First Wall',fwdzu,zu2)
        call obuild(outfile,'Void 1',v1dzu,zu3)
        call obuild(outfile,'Blanket',bldzu,zu4)
        call obuild(outfile,'Void 2',v2dzu,zu5)
        call obuild(outfile,'Shield',shdzu,zu6)
        call obuild(outfile,'Void 3',v3dzu,zu7)

    else

        call oheadr(outfile,'Vertical Build - Midplane')
        write(outfile,30)

        call obuild(outfile,'Base of device',0.0D0,-zl7)
        call obuild(outfile,'Void 3',v3dzl,-zl6)
        call obuild(outfile,'Shield',shdzl,-zl5)
        call obuild(outfile,'Void 2',v2dzl,-zl4)
        call obuild(outfile,'Blanket',bldzl,-zl3)
        call obuild(outfile,'Void 1',v1dzl,-zl2)
        call obuild(outfile,'First Wall',fwdzl,-zl1)
        call obuild(outfile,'Chamber',chdzl,0.0D0)
        call obuild(outfile,'Chamber',chdzu,zu1)
        call obuild(outfile,'First Wall',fwdzu,zu2)
        call obuild(outfile,'Void 1',v1dzu,zu3)
        call obuild(outfile,'Blanket',bldzu-bldzu,zu4-bldzu)
        call obuild(outfile,'Void 2',v2dzu,zu5-bldzu)
        call obuild(outfile,'Shield',shdzu,zu6-bldzu)
        call obuild(outfile,'Void 3',v3dzu+bldzu,zu7)

        call oheadr(outfile,'Vertical Build - Edge')
        write(outfile,30)

        call obuild(outfile,'Base of device',0.0D0,-zl7)
        call obuild(outfile,'Void 3',v3dzl,-zl6)
        call ovarre(mfile,'Void 3 lower (m)','(v3dzl)',v3dzl)
        call obuild(outfile,'Shield',shdzl,-zl5)
        call ovarre(mfile,'Shield lower (m)','(shdzl)',shdzl)
        call obuild(outfile,'Void 2',v2dzl,-zl4)
        call ovarre(mfile,'Void 2 lower (m)','(v2dzl)',v2dzl)
        call obuild(outfile,'Blanket',bldzl,-zl3)
        call ovarre(mfile,'Blanket lower (m)','(bldzl)',bldzl)
        call obuild(outfile,'Void 1',v1dzl,-zl2)
        call ovarre(mfile,'Void 1 lower (m)','(v1dzl)',v1dzl)
        call obuild(outfile,'First Wall',fwdzl,-zl1)
        call ovarre(mfile,'First Wall lower (m)','(fwdzl)',fwdzl)
        call obuild(outfile,'Chamber',chdzl,0.0D0)
        call ovarre(mfile,'Chamber lower (m)','(chdzl)',chdzl)
        call obuild(outfile,'Chamber',chdzu,zu1)
        call ovarre(mfile,'Chamber upper (m)','(chdzu)',chdzu)
        call obuild(outfile,'First Wall',fwdzu,zu2)
        call ovarre(mfile,'First Wall upper (m)','(fwdzu)',fwdzu)
        call obuild(outfile,'Void 1',v1dzu,zu3)
        call ovarre(mfile,'Void 1 upper (m)','(v1dzu)',v1dzu)
        call obuild(outfile,'Blanket',bldzu,zu4)
        call ovarre(mfile,'Blanket upper (m)','(bldzu)',bldzu)
        call obuild(outfile,'Void 2',v2dzu,zu5)
        call ovarre(mfile,'Void 2 upper (m)','(v2dzu)',v2dzu)
        call obuild(outfile,'Shield',shdzu,zu6)
        call ovarre(mfile,'Shield upper (m)','(shdzu)',shdzu)
        call obuild(outfile,'Void 3',v3dzu,zu7)
        call ovarre(mfile,'Void 3 upper (m)','(v3dzu)',v3dzu)

    end if

    !  Print matrix of material volumes

    call oheadr(outfile,'Material volumes')

    write(outfile,*) '         Chamber  1st wall  Void 1  Blanket  ' &
         //' Void 2   Shield   Void 3'
    write(outfile,'(A9,7(1pe9.2))') 'void     ', &
         chmatv(0), &
         (fwmatv(1,0)+fwmatv(2,0)+fwmatv(3,0)), &
         (v1matv(1,0)+v1matv(2,0)+v1matv(3,0)), &
         (blmatv(1,0)+blmatv(2,0)+blmatv(3,0)), &
         (v2matv(1,0)+v2matv(2,0)+v2matv(3,0)), &
         (shmatv(1,0)+shmatv(2,0)+shmatv(3,0)), &
         (v3matv(1,0)+v3matv(2,0)+v3matv(3,0))
    write(outfile,'(A9,7(1pe9.2))') 'steel    ', &
         chmatv(1), &
         (fwmatv(1,1)+fwmatv(2,1)+fwmatv(3,1)), &
         (v1matv(1,1)+v1matv(2,1)+v1matv(3,1)), &
         (blmatv(1,1)+blmatv(2,1)+blmatv(3,1)), &
         (v2matv(1,1)+v2matv(2,1)+v2matv(3,1)), &
         (shmatv(1,1)+shmatv(2,1)+shmatv(3,1)), &
         (v3matv(1,1)+v3matv(2,1)+v3matv(3,1))

    write(outfile,'(A9,7(1pe9.2))') 'carbon   ', &
         chmatv(2), &
         (fwmatv(1,2)+fwmatv(2,2)+fwmatv(3,2)), &
         (v1matv(1,2)+v1matv(2,2)+v1matv(3,2)), &
         (blmatv(1,2)+blmatv(2,2)+blmatv(3,2)), &
         (v2matv(1,2)+v2matv(2,2)+v2matv(3,2)), &
         (shmatv(1,2)+shmatv(2,2)+shmatv(3,2)), &
         (v3matv(1,2)+v3matv(2,2)+v3matv(3,2))

    write(outfile,'(A9,7(1pe9.2))') 'FLiBe    ', &
         chmatv(3), &
         (fwmatv(1,3)+fwmatv(2,3)+fwmatv(3,3)), &
         (v1matv(1,3)+v1matv(2,3)+v1matv(3,3)), &
         (blmatv(1,3)+blmatv(2,3)+blmatv(3,3)), &
         (v2matv(1,3)+v2matv(2,3)+v2matv(3,3)), &
         (shmatv(1,3)+shmatv(2,3)+shmatv(3,3)), &
         (v3matv(1,3)+v3matv(2,3)+v3matv(3,3))

    write(outfile,'(A9,7(1pe9.2))') 'Li2O     ', &
         chmatv(4), &
         (fwmatv(1,4)+fwmatv(2,4)+fwmatv(3,4)), &
         (v1matv(1,4)+v1matv(2,4)+v1matv(3,4)), &
         (blmatv(1,4)+blmatv(2,4)+blmatv(3,4)), &
         (v2matv(1,4)+v2matv(2,4)+v2matv(3,4)), &
         (shmatv(1,4)+shmatv(2,4)+shmatv(3,4)), &
         (v3matv(1,4)+v3matv(2,4)+v3matv(3,4))

    write(outfile,'(A9,7(1pe9.2))') 'concrete ', &
         chmatv(5), &
         (fwmatv(1,5)+fwmatv(2,5)+fwmatv(3,5)), &
         (v1matv(1,5)+v1matv(2,5)+v1matv(3,5)), &
         (blmatv(1,5)+blmatv(2,5)+blmatv(3,5)), &
         (v2matv(1,5)+v2matv(2,5)+v2matv(3,5)), &
         (shmatv(1,5)+shmatv(2,5)+shmatv(3,5)), &
         (v3matv(1,5)+v3matv(2,5)+v3matv(3,5))

    write(outfile,'(A9,7(1pe9.2))') 'helium   ', &
         chmatv(6), &
         (fwmatv(1,6)+fwmatv(2,6)+fwmatv(3,6)), &
         (v1matv(1,6)+v1matv(2,6)+v1matv(3,6)), &
         (blmatv(1,6)+blmatv(2,6)+blmatv(3,6)), &
         (v2matv(1,6)+v2matv(2,6)+v2matv(3,6)), &
         (shmatv(1,6)+shmatv(2,6)+shmatv(3,6)), &
         (v3matv(1,6)+v3matv(2,6)+v3matv(3,6))

    write(outfile,'(A9,7(1pe9.2))') 'xenon    ', &
         chmatv(7), &
         (fwmatv(1,7)+fwmatv(2,7)+fwmatv(3,7)), &
         (v1matv(1,7)+v1matv(2,7)+v1matv(3,7)), &
         (blmatv(1,7)+blmatv(2,7)+blmatv(3,7)), &
         (v2matv(1,7)+v2matv(2,7)+v2matv(3,7)), &
         (shmatv(1,7)+shmatv(2,7)+shmatv(3,7)), &
         (v3matv(1,7)+v3matv(2,7)+v3matv(3,7))

    write(outfile,'(A9,7(1pe9.2))') 'lithium  ', &
         chmatv(8), &
         (fwmatv(1,8)+fwmatv(2,8)+fwmatv(3,8)), &
         (v1matv(1,8)+v1matv(2,8)+v1matv(3,8)), &
         (blmatv(1,8)+blmatv(2,8)+blmatv(3,8)), &
         (v2matv(1,8)+v2matv(2,8)+v2matv(3,8)), &
         (shmatv(1,8)+shmatv(2,8)+shmatv(3,8)), &
         (v3matv(1,8)+v3matv(2,8)+v3matv(3,8))

  contains

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine bld2019
     !! Routine to create the build of a 2019 inertial fusion energy
     !! device, and to calculate the material volumes for the device core
     !! author: S I Muldrew, CCFE, Culham Science Centre
     !! None
     !! This routine constructs the build of a modern inertial fusion energy
     !! device, assumed to be cylindrically-symmetric, with a pool at bottom
     !! and top corners and with a lower shield at the centre.  See diagram
     !! attached to Issue #907.
     !! Issue #907
     !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
     !
     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use constants, only: pi
    use build_variables, only: fwarea
    use error_handling, only: report_error
    use buildings_variables, only: trcl, stcl
    use fwbs_variables, only: tbr, emult
    use ife_variables, only: blmatf, bldrc, lipmw, etali, taufall, rrmax, chvol, &
      fwvol, v1vol, blvol, v2vol, shvol, flirad, v3vol, maxmat, chmatf, fwmatf, &
      v1matf, v2matf, shmatf, v3matf
    implicit none

     !  Arguments

     !  Local variables

     integer :: i,j
     real(dp), save :: g, vel, acurt, mdot, phi, sang, li_frac

     ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Check input
     if ((fwdr.gt.1.0D-9).or.(v1dr.gt.1.0d-9)) then
          call report_error(230)
     else if ((fwdzu.gt.1.0D-9).or.(v1dzu.gt.1.0D-9).or.(v2dzu.gt.1.0D-9)) then
          call report_error(231)
     else if ((fwdzl.gt.1.0D-9).or.(v1dzl.gt.1.0d-9).or.(v2dzu.gt.1.0D-9)) then
          call report_error(232)
     end if

     ! Lithium Pump

     ! Acceleration due to gravity (m/s^2)
     g = 9.81D0

     ! Velocity
     vel = sqrt(2.0D0*g*(chdzu+bldzu))

     ! Lithium Fraction
     blmatf(1,8) = 0.91*sqrt(bldzu/(chdzu+bldzu))
     blmatf(1,0) = 1.0D0 - blmatf(1,8)

     ! Spatial Thickness
     bldr = bldrc / blmatf(1,8)

     ! Area
     acurt = pi * ((chrad+bldr)**2.0D0-chrad**2.0D0)

     ! Mass Flow
     mdot = 512.0D0 * vel * blmatf(1,8) * acurt

     ! Pump Power (MW)
     lipmw = 1.0D-6 * mdot * g * (chdzl+chdzu+bldzu+bldzl) / etali

     ! Fall Time
     taufall = 2.0D0 * (chdzl+chdzu+bldzu) / vel

     rrmax = 1.0D0 / taufall

     ! TBR and Emult model was for spherical lithium
     ! Remove reactor head
     phi = atan(r1/zu1)
     sang = 1.0D0 - cos(phi)
     li_frac = 1.0D0 - 0.5D0 * sang

     ! TBR
     tbr = 3.7418D0 * (1.0D0 / (1.0D0 + exp(-2.6366D0 * bldrc)) - 0.5D0)
     tbr = tbr * li_frac

     ! Energy Multiplication
     emult = 2.2414D0 * (1.0D0 / (1.0D0 + exp(-3.0038D0 * bldrc)) - 0.5D0)
     emult = emult * li_frac

     !  Radial build

     r1 = chrad
     r2 = r1 + fwdr
     r3 = r2 + v1dr
     r4 = r3 + bldr
     r5 = r4 + v2dr
     r6 = r5 + shdr
     r7 = r6 + v3dr

     !  Vertical build (below midplane)

     zl1 = chdzl
     zl2 = zl1 + fwdzl
     zl3 = zl2 + v1dzl
     zl4 = zl3 + bldzl
     zl5 = zl4 + v2dzl
     zl6 = zl5 + shdzl
     zl7 = zl6 + v3dzl

     !  Vertical build (above midplane)

     zu1 = chdzu
     zu2 = zu1 + fwdzu
     zu3 = zu2 + v1dzu
     zu4 = zu3 + bldzu
     zu5 = zu4 + v2dzu
     zu6 = zu5 + shdzu

     v3dzu = (zu6+zl6) + trcl + stcl + 5.1D0 + 9.41D-6 * 1.0D5

     zu7 = zu6 + v3dzu

     !  Component volumes
     !  The following notation applies below:
     !  J=1 : side part
     !  J=2 : top part
     !  J=3 : bottom part

     !  Chamber

     chvol = pi * r1*r1 * (zu1 + zl1)

     !  First wall

     fwvol(1) = pi * (r2*r2 - r1*r1) * (zu1 + zl1)
     fwvol(2) = pi * r2*r2 * (zu2 - zu1)
     fwvol(3) = pi * r2*r2 * (zl2 - zl1)

     !  First void

     v1vol(1) = pi * (r3*r3 - r2*r2) * (zu2 + zl2)
     v1vol(2) = pi * r3*r3 * (zu3 - zu2)
     v1vol(3) = pi * r3*r3 * (zl3 - zl2)

     !  Blanket
     !  Radial Blanket - between void 2 and chamber
     blvol(1) = pi * (r4*r4 - r3*r3) * (zu3 + zl3)
     !  Upper Blanket - Pool radially between shield and
     !  chamber of input height.
     blvol(2) = pi * (r5*r5 - r3*r3) * bldzu
     !  Lower Blanket - Pool filling base of device
     blvol(3) = pi * r5*r5 * (zl4 - zl3)

     !  Second void

    v2vol(1) = pi * (r5*r5 - r4*r4) * (chdzl+chdzu)
    v2vol(2) = 0.0D0
    v2vol(3) = 0.0D0

    !  Shield
    shvol(1) = pi * (r6*r6 - r5*r5) * (zu5 + zl5)
    ! Top Section is in three parts to account for the dip at
    ! the centre.  The first is the horizontal top, the second is the
    ! horizontal
    shvol(2) = pi * (((r6*r6 - (chrad-shdr)*(chrad-shdr)) * shdzu) &
      + ((r1*r1 - flirad*flirad) * shdzu) &
      + ( (r1*r1 - (r1-shdzu)*(r1-shdzu)) *(bldzu-shdzu) ))
    shvol(3) = pi * r6*r6 * (zl6 - zl5)

    !  Third void

    v3vol(1) = pi * (r7*r7 - r6*r6) * (zu6 + zl6)
    v3vol(2) = pi * r7*r7 * (zu7 - zu6) + pi * ((r1-shdzu)*(r1-shdzu) - flirad*flirad) * bldzu
    v3vol(3) = pi * r7*r7 * (zl7 - zl6)

     !  Material volumes

     do i = 0,maxmat
        chmatv(i) = max(0.0D0, chvol * chmatf(i))
        do j = 1,3
           fwmatv(j,i) = max(0.0D0, fwvol(j) * fwmatf(j,i))
           v1matv(j,i) = max(0.0D0, v1vol(j) * v1matf(j,i))
           blmatv(j,i) = max(0.0D0, blvol(j) * blmatf(j,i))
           v2matv(j,i) = max(0.0D0, v2vol(j) * v2matf(j,i))
           shmatv(j,i) = max(0.0D0, shvol(j) * shmatf(j,i))
           v3matv(j,i) = max(0.0D0, v3vol(j) * v3matf(j,i))
        end do
     end do

     !  First wall area
     !  The chamber is surrounded by liquid on three sides
     !  with only the top being solid.  This is considered part
     !  of the shield. There is a target injector tube at the
     !  centre of this area.
     fwarea = pi * (r1*r1 - flirad*flirad)

   end subroutine bld2019

   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine genbld

      !! Routine to create the build of a generic inertial fusion energy
      !! device, and to calculate the material volumes for the device core
      !! author: P J Knight, CCFE, Culham Science Centre
      !! None
      !! This routine constructs the build of a generic inertial fusion energy
      !! device, assumed to be cylindrically-symmetric, and to calculate
      !! the material volumes for the device core.
      !! F/MI/PJK/LOGBOOK12, p.52
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use constants, only: pi
      use build_variables, only: fwarea
      use ife_variables, only: chvol, fwvol, v1vol, blvol, v2vol, shvol, v3vol, &
        maxmat, chmatf, fwmatf, v1matf, blmatf, v2matf, shmatf, v3matf

      implicit none

      !  Local variables
      integer :: i,j

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Radial build

      r1 = chrad
      r2 = r1 + fwdr
      r3 = r2 + v1dr
      r4 = r3 + bldr
      r5 = r4 + v2dr
      r6 = r5 + shdr
      r7 = r6 + v3dr

      !  Vertical build (below midplane)

      zl1 = chdzl
      zl2 = zl1 + fwdzl
      zl3 = zl2 + v1dzl
      zl4 = zl3 + bldzl
      zl5 = zl4 + v2dzl
      zl6 = zl5 + shdzl
      zl7 = zl6 + v3dzl

      !  Vertical build (above midplane)

      zu1 = chdzu
      zu2 = zu1 + fwdzu
      zu3 = zu2 + v1dzu
      zu4 = zu3 + bldzu
      zu5 = zu4 + v2dzu
      zu6 = zu5 + shdzu
      zu7 = zu6 + v3dzu

      !  Component volumes
      !  The following notation applies below:
      !  J=1 : side part
      !  J=2 : top part
      !  J=3 : bottom part

      !  Chamber

      chvol = pi * r1*r1 * (zu1 + zl1)

      !  First wall

      fwvol(1) = pi * (r2*r2 - r1*r1) * (zu1 + zl1)
      fwvol(2) = pi * r2*r2 * (zu2 - zu1)
      fwvol(3) = pi * r2*r2 * (zl2 - zl1)

      !  First void

      v1vol(1) = pi * (r3*r3 - r2*r2) * (zu2 + zl2)
      v1vol(2) = pi * r3*r3 * (zu3 - zu2)
      v1vol(3) = pi * r3*r3 * (zl3 - zl2)

      !  Blanket

      blvol(1) = pi * (r4*r4 - r3*r3) * (zu3 + zl3)
      blvol(2) = pi * r4*r4 * (zu4 - zu3)
      blvol(3) = pi * r4*r4 * (zl4 - zl3)

      !  Second void

      v2vol(1) = pi * (r5*r5 - r4*r4) * (zu4 + zl4)
      v2vol(2) = pi * r5*r5 * (zu5 - zu4)
      v2vol(3) = pi * r5*r5 * (zl5 - zl4)

      !  Shield

      shvol(1) = pi * (r6*r6 - r5*r5) * (zu5 + zl5)
      shvol(2) = pi * r6*r6 * (zu6 - zu5)
      shvol(3) = pi * r6*r6 * (zl6 - zl5)

      !  Third void

      v3vol(1) = pi * (r7*r7 - r6*r6) * (zu6 + zl6)
      v3vol(2) = pi * r7*r7 * (zu7 - zu6)
      v3vol(3) = pi * r7*r7 * (zl7 - zl6)

      !  Material volumes

      do i = 0,maxmat
         chmatv(i) = max(0.0D0, chvol * chmatf(i))
         do j = 1,3
            fwmatv(j,i) = max(0.0D0, fwvol(j) * fwmatf(j,i))
            v1matv(j,i) = max(0.0D0, v1vol(j) * v1matf(j,i))
            blmatv(j,i) = max(0.0D0, blvol(j) * blmatf(j,i))
            v2matv(j,i) = max(0.0D0, v2vol(j) * v2matf(j,i))
            shmatv(j,i) = max(0.0D0, shvol(j) * shmatf(j,i))
            v3matv(j,i) = max(0.0D0, v3vol(j) * v3matf(j,i))
         end do
      end do

      !  First wall area

      fwarea = 2.0D0*pi*r1 * ((zu1 + zl1) + r1)

    end subroutine genbld

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine osibld
      !! Routine to create the build of an inertial fusion energy
      !! device, based on the design of the OSIRIS study,
      !! and to calculate the material volumes for the device core
      !! author: P J Knight, CCFE, Culham Science Centre
      !! None
      !! This routine constructs the build of an inertial fusion energy
      !! device, based on the design of the OSIRIS study, and to calculate
      !! the material volumes for the device core.
      !! F/MI/PJK/LOGBOOK12, p.56
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use constants, only: pi
      use build_variables, only: fwarea

      implicit none

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Careful choice of thicknesses, and assuming that the FLiBe
      !  inlet radius is small, allows the generic build calculation
      !  to be roughly applicable.
      call genbld

      !  First wall area: no true first wall at bottom of chamber
      fwarea = 2.0D0*pi*r1*(zu1 + zl1) + pi*r1*r1

    end subroutine osibld

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine sombld
      !! Routine to create the build of an inertial fusion energy
      !! device, based on the design of the SOMBRERO study,
      !! and to calculate the material volumes for the device core
      !! author: P J Knight, CCFE, Culham Science Centre
      !! None
      !! This routine constructs the build of an inertial fusion energy
      !! device, based on the design of the SOMBRERO study, and to calculate
      !! the material volumes for the device core.
      !! Sviatoslavsky et al, Fusion Technology vol.21 (1992) 1470
      !! F/MI/PJK/LOGBOOK12, p.53
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use constants, only: pi
      use build_variables, only: fwarea
      use ife_variables, only: chvol, fwvol, v1vol, blvol, v2vol, somtdr, sombdr, &
        shvol, v3vol, maxmat, chmatf, fwmatf, blmatf, v1matf, v2matf, shmatf, v3matf

      implicit none

      !  Local variables

      real(dp), parameter :: third  = 1.0D0/3.0D0
      real(dp) :: chcylh,ddz,dvol
      integer :: i,j

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Radial build
      r1 = chrad
      r2 = r1 + fwdr
      r3 = r2 + v1dr
      r4 = r3 + bldr
      r5 = r4 + v2dr
      r6 = r5 + shdr
      r7 = r6 + v3dr

      !  Vertical build (below midplane)
      zl1 = chdzl
      zl2 = zl1 + fwdzl
      zl3 = zl2 + v1dzl
      zl4 = zl3 + bldzl
      zl5 = zl4 + v2dzl
      zl6 = zl5 + shdzl
      zl7 = zl6 + v3dzl

      !  Vertical build (above midplane)
      zu1 = chdzu
      zu2 = zu1 + fwdzu
      zu3 = zu2 + v1dzu
      zu4 = zu3 + bldzu
      zu5 = zu4 + v2dzu
      zu6 = zu5 + shdzu
      zu7 = zu6 + v3dzu

      !  The SOMBRERO chamber is made up of a cylindrical first wall/
      !  blanket, with conical regions above and below. Outside this is
      !  a cylindrical shield.

      !  Component volumes
      !  The following notation applies below:
      !  J=1 : side part
      !  J=2 : top part
      !  J=3 : bottom part

      !  Chamber : CHCYLH is the height of the cylindrical part
      chcylh = chdzu + chdzl - 2.0D0*chrad

      chvol = pi * r1*r1 * (chcylh + 2.0D0*third * chrad)

      !  First wall
      fwvol(1) = pi * (r2*r2 - r1*r1) * chcylh
      fwvol(2) = third * pi * ( r2*r2*(chrad+fwdzu) - r1*r1*chrad )
      fwvol(3) = third * pi * ( r2*r2*(chrad+fwdzl) - r1*r1*chrad )

      !  First void
      v1vol(1) = pi * (r3*r3 - r2*r2) * chcylh
      v1vol(2) = third * pi * ( r3*r3*(chrad+fwdzu+v1dzu) -  &
           r2*r2*(chrad+fwdzu) )
      v1vol(3) = third * pi * ( r3*r3*(chrad+fwdzl+v1dzl) -  &
           r2*r2*(chrad+fwdzl) )

      !  Blanket:  SOMTDR and SOMBDR are the radii of the cylindrical
      !  sections at the top/bottom of the blanket
      !  DDZ = Height of top cylindrical section (by similar triangles)
      !  DVOL = Volume of top cylindrical section, less the internal cone

      blvol(1) = pi * (r4*r4 - r3*r3) * chcylh

      blvol(2) = third * pi * ( r4*r4*(chrad+fwdzu+v1dzu+bldzu) -  &
           r3*r3*(chrad+fwdzu+v1dzu) )
      ddz = (chrad+fwdzu+v1dzu+bldzu)/(chrad+fwdr+v1dr+bldr)*somtdr
      dvol = 2.0D0*third * pi * somtdr*somtdr * ddz

      blvol(2) = blvol(2) + dvol

      !  Ditto for bottom region...

      blvol(3) = third * pi * ( r4*r4*(chrad+fwdzl+v1dzl+bldzl) -  &
           r3*r3*(chrad+fwdzl+v1dzl) )
      ddz = (chrad+fwdzl+v1dzl+bldzl)/(chrad+fwdr+v1dr+bldr)*sombdr
      dvol = 2.0D0*third * pi * sombdr*sombdr * ddz

      blvol(3) = blvol(3) + dvol

      !  Second void
      v2vol(1) = pi * (r5*r5 - r4*r4) * chcylh
      v2vol(2) = pi * r5*r5 * (zu5 - chdzu + chrad) - ( &
           fwvol(2) + v1vol(2) + blvol(2) + (third*pi*r1*r1*chrad) )
      v2vol(3) = pi * r5*r5 * (zl5 - chdzl + chrad) - ( &
           fwvol(3) + v1vol(3) + blvol(3) + (third*pi*r1*r1*chrad) )

      !  Shield
      shvol(1) = pi * (r6*r6 - r5*r5) * (zu6 + zl6)
      shvol(2) = pi * r5*r5 * (zu6 - zu5)
      shvol(3) = pi * r5*r5 * (zl6 - zl5)

      !  Third void
      v3vol(1) = pi * (r7*r7 - r6*r6) * (zu7 + zl7)
      v3vol(2) = pi * r6*r6 * (zu7 - zu6)
      v3vol(3) = pi * r6*r6 * (zl7 - zl6)

      !  Material volumes
      do i = 0,maxmat
        chmatv(i) = max(0.0D0, chvol * chmatf(i))
        do j = 1,3
          fwmatv(j,i) = max(0.0D0, fwvol(j) * fwmatf(j,i))
          v1matv(j,i) = max(0.0D0, v1vol(j) * v1matf(j,i))
          blmatv(j,i) = max(0.0D0, blvol(j) * blmatf(j,i))
          v2matv(j,i) = max(0.0D0, v2vol(j) * v2matf(j,i))
          shmatv(j,i) = max(0.0D0, shvol(j) * shmatf(j,i))
          v3matv(j,i) = max(0.0D0, v3vol(j) * v3matf(j,i))
        end do
      end do

      !  First wall area
      fwarea = 2.0D0*pi*r1*( (zu1 + zl1) + r1*sqrt(2.0D0) )

    end subroutine sombld

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine hylbld
      !! Routine to create the build of an inertial fusion energy
      !! device, based on the design of the HYLIFE-II study,
      !! and to calculate the material volumes for the device core
      !! author: P J Knight, CCFE, Culham Science Centre
      !! None
      !! This routine constructs the build of an inertial fusion energy
      !! device, based on the design of the HYLIFE-II study, and to calculate
      !! the material volumes for the device core.
      !! Moir, Fusion Technology vol.21 (1992) 1475
      !! F/MI/PJK/LOGBOOK12, p.57
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use constants, only: pi
      use build_variables, only: fwarea
      use ife_variables, only: chvol, fwvol, v1vol, blvol, v2vol, shvol, v3vol, &
        maxmat, chmatf, fwmatf, blmatf, v1matf, v2matf, shmatf, v3matf, flirad

      implicit none

      !  Local variables

      real(dp), parameter :: third = 1.0D0/3.0D0
      real(dp) :: chcylh,ddz,dvol
      integer :: i,j

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Radial build
      r1 = chrad
      r2 = r1 + fwdr
      r3 = r2 + v1dr
      r4 = r3 + bldr
      r5 = r4 + v2dr
      r6 = r5 + shdr
      r7 = r6 + v3dr

      !  Vertical build (below midplane)
      zl1 = chdzl
      zl2 = zl1 + fwdzl
      zl3 = zl2 + v1dzl
      zl4 = zl3 + bldzl
      zl5 = zl4 + v2dzl
      zl6 = zl5 + shdzl
      zl7 = zl6 + v3dzl

      !  Vertical build (above midplane)
      zu1 = chdzu
      zu2 = zu1 + fwdzu
      zu3 = zu2 + v1dzu
      zu4 = zu3 + bldzu
      zu5 = zu4 + v2dzu
      zu6 = zu5 + shdzu
      zu7 = zu6 + v3dzu

      !  The HYLIFE-II chamber is assumed to be mostly cylindrical, but
      !  with a conical region below the midplane that causes the Flibe
      !  to flow downwards and outwards towards the outlet.

      !  Component volumes
      !  The following notation applies below:
      !  J=1 : side part
      !  J=2 : top part
      !  J=3 : bottom part
      chvol = pi * r1*r1 * ( (zu1 + zl5) - third*(zl5 - zl1) )

      !  First wall
      !  FLIRAD is the radius of the Flibe inlet
      fwvol(1) = pi * (r2*r2 - r1*r1) * (zu2 + zl5)
      fwvol(2) = pi * (r1*r1 - flirad*flirad) * (zu2 - zu1)
      fwvol(3) = third * pi * ( r2*r2*(zl5-zl1) - r1*r1*(zl5-zl2) )

      !  First void
      v1vol(1) = pi * (r3*r3 - r2*r2) * (zu2 + zl3)
      v1vol(2) = pi * (r4*r4 - flirad*flirad) * (zu3 - zu2)
      v1vol(3) = third * pi * r1*r1 * (zl3 - zl2)

      !  Blanket
      blvol(1) = pi * (r4*r4 - r3*r3) * (zu2 + zl3)
      blvol(2) = pi * (r4*r4 - flirad*flirad) * (zu4 - zu3)
      blvol(3) = pi * r4*r4 * (zl4 - zl3)

      !  Second void
      v2vol(1) = pi * (r5*r5 - r4*r4) * (zu4 + zl4)
      v2vol(2) = pi * (r5*r5 - flirad*flirad) * (zu5 - zu4)
      v2vol(3) = pi * r5*r5 * (zl5 - zl4)

      !  Shield
      shvol(1) = pi * (r6*r6 - r5*r5) * (zu5 + zl5)
      shvol(2) = pi * r6*r6 * (zu6 - zu5)
      shvol(3) = pi * r6*r6 * (zl6 - zl5)

      !  Third void
      v3vol(1) = pi * (r7*r7 - r6*r6) * (zu6 + zl6)
      v3vol(2) = pi * r7*r7 * (zu7 - zu6)
      v3vol(3) = pi * r7*r7 * (zl7 - zl6)

      !  Material volumes
      do i = 0,maxmat
        chmatv(i) = max(0.0D0, chvol * chmatf(i))
        do j = 1,3
          fwmatv(j,i) = max(0.0D0, fwvol(j) * fwmatf(j,i))
          v1matv(j,i) = max(0.0D0, v1vol(j) * v1matf(j,i))
          blmatv(j,i) = max(0.0D0, blvol(j) * blmatf(j,i))
          v2matv(j,i) = max(0.0D0, v2vol(j) * v2matf(j,i))
          shmatv(j,i) = max(0.0D0, shvol(j) * shmatf(j,i))
          v3matv(j,i) = max(0.0D0, v3vol(j) * v3matf(j,i))
        end do
      end do

      !  First wall area
      fwarea = 2.0D0 * pi * r1 * (zu1 + zl5)
      fwarea = fwarea + pi * (r1*r1 - flirad*flirad)
      fwarea = fwarea + pi * r1 * sqrt(r1*r1 + (zl3-zl1)**2)

    end subroutine hylbld

  end subroutine ifebld

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifestr
    !! Routine to calculate the support structural masses for the core of
    !! an Inertial Fusion Energy power plant
    !! author: P J Knight, CCFE, Culham Science Centre
    !! This routine calculates the support structural masses for the core of
    !! an Inertial Fusion Energy power plant.
    !! <P>In fact, the output masses are all trivially zero, as they are
    !! magnetic fusion specific. F/MI/PJK/LOGBOOK12, p.87
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use structure_variables, only: aintmass, clgsmass, coldmass, fncmass, gsmass

    implicit none

    !  Set all outputs to zero, as they are magnetic fusion specific
    aintmass = 0.0D0
    clgsmass = 0.0D0
    coldmass = 0.0D0
    fncmass  = 0.0D0
    gsmass   = 0.0D0

  end subroutine ifestr

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifetgt()
    !! Routine to calculate the power requirements of the target
    !! delivery system and the target factory
    !! author: P J Knight, CCFE, Culham Science Centre
    !! This routine calculates the power requirements of the target
    !! delivery system and the target factory, for an Inertial
    !! Fusion Energy power plant.
    !! F/MI/PJK/LOGBOOK12, pp.87-88
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use ife_variables, only: tfacmw, ptargf, reprat

    implicit none

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Target delivery system power (MWe) - effectively negligible

    !  tdspmw = 1.0D-2

    !  Target factory power (MWe)
    !  Assumed to scale with repetition rate (not quite linearly)
    tfacmw = ptargf * (reprat/6.0D0)**0.7D0

  end subroutine ifetgt

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifefbs(outfile,iprint)
    !! Routine to calculate the first wall, blanket and shield volumes,
    !! masses and other parameters, for an Inertial Fusion Energy device
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer : output file unit
    !! iprint : input integer : switch for writing to output file (1=yes)
    !! This routine calculates the first wall, blanket and shield volumes,
    !! masses and other parameters, for an Inertial Fusion Energy device.
    !! F/MI/PJK/LOGBOOK12, p.86
    !! Moir et al., Fusion Technology, vol.25 (1994) p.5
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use build_variables, only: fwarea
    use error_handling, only: report_error, fdiags
    use process_output, only: ovarre, oheadr
    use cost_variables, only: tlife, abktflnc, cfactr
    use fwbs_variables, only: denstl, fwmass, whtblkt, whtshld, whtblbe, whtblvd, &
      whtblss, wtblli2o, whtblli, bktlife, fwlife
    use ife_variables, only: maxmat, chmatm, chmatv, fwmatm, fwmatv, v1matm, &
      v1matv, blmatm, blmatv, v2matm, v2matv, shmatm, shmatv, v3matm, v3matv, &
      mflibe, fbreed, ifetyp
    use physics_variables, only: wallmw

    implicit none

    !  Arguments
    integer, intent(in) :: outfile,iprint

    !  Local variables

    real(dp) :: den,life
    real(dp), dimension(0:maxmat) :: matden
    integer :: i,j

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Material densities
    !  0 = void
    !  1 = steel
    !  2 = carbon
    !  3 = FLiBe (inferred from Moir et al)
    !  4 = Li2O
    !  5 = concrete
    !  6 = helium (at typical coolant temperatures)
    !  7 = xenon (taken as ten times the normal tabulated value)
    !  8 = lithium (liquid, source Wikipedia)

    matden(0) = 0.0D0
    matden(1) = denstl
    matden(2) = 2300.0D0
    matden(3) = 2020.0D0
    matden(4) = 2010.0D0
    matden(5) = 2400.0D0
    matden(6) = 1.517D0
    matden(7) = 55.0D0
    matden(8) = 512.0D0

    !  Material masses
    do i = 0,maxmat
      den = matden(i)
      chmatm(i) = chmatv(i) * den
      do j = 1,3
        fwmatm(j,i) = fwmatv(j,i) * den
        v1matm(j,i) = v1matv(j,i) * den
        blmatm(j,i) = blmatv(j,i) * den
        v2matm(j,i) = v2matv(j,i) * den
        shmatm(j,i) = shmatv(j,i) * den
        v3matm(j,i) = v3matv(j,i) * den
      end do
    end do

    !  Total masses of components (excluding coolant)
    fwmass = 0.0D0
    whtblkt = 0.0D0
    whtshld = 0.0D0
    do i = 1,5
      do j = 1,3
        fwmass = fwmass + fwmatm(j,i)
        whtblkt = whtblkt + blmatm(j,i)
        whtshld = whtshld + shmatm(j,i)
      end do
    end do

    !  Other masses
    whtblbe = 0.0D0
    whtblvd = 0.0D0
    whtblss = 0.0D0
    wtblli2o = 0.0D0
    whtblli = 0.0D0
    do j = 1,3
      whtblss = whtblss + blmatm(j,1)
      wtblli2o = wtblli2o + blmatm(j,4)
      whtblli = whtblli + blmatm(j,8)
    end do

    !  Total mass of FLiBe
    mflibe = chmatm(3)
    do j = 1,3
      mflibe = mflibe + fwmatm(j,3) + v1matm(j,3) + blmatm(j,3) + &
        v2matm(j,3) + shmatm(j,3) + v3matm(j,3)
    end do

    !  A fraction FBREED of the total breeder inventory is outside the
    !  core region, i.e. is in the rest of the heat transport system
    if ((fbreed < 0.0D0).or.(fbreed > 0.999D0)) then
      fdiags(1) = fbreed ; call report_error(26)
    end if

    !  Following assumes that use of FLiBe and Li2O are
    !  mutually exclusive
    mflibe = mflibe / (1.0D0 - fbreed)
    wtblli2o = wtblli2o / (1.0D0 - fbreed)
    whtblli = whtblli / (1.0D0 - fbreed)

    !  Blanket and first wall lifetimes (HYLIFE-II: = plant life)
    if ((ifetyp == 3).or.(ifetyp == 4)) THEN
      life = tlife
    else
      life = min( tlife, abktflnc/(wallmw*cfactr) )
    end if

    bktlife = life
    fwlife = life

    !  Cryostat mass (=zero)
    !cryomass = 0.0D0

    if (iprint == 0) return

    !  Output section
    call oheadr(outfile,'First Wall, Blanket, Shield')
    call ovarre(outfile,'First wall area (m2)','(fwarea)',fwarea)
    call ovarre(outfile,'First wall mass (kg)','(fwmass)',fwmass)
    call ovarre(outfile,'Blanket mass (kg)','(whtblkt)',whtblkt)
    call ovarre(outfile,'Blanket lithium mass (kg)','(whtblli)',whtblli)
    call ovarre(outfile,'Total mass of FLiBe (kg)','(mflibe)',mflibe)
    call ovarre(outfile,'Shield mass (kg)','(whtshld)',whtshld)

  end subroutine ifefbs

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifepw1

    !! Routine to calculate the first part of the heat transport
    !! and plant power balance constituents, for an IFE power plant
    !! author: P J Knight, CCFE, Culham Science Centre
    !! None
    !! This routine calculates the first part of the heat transport
    !! and plant power balance constituents, for an IFE power plant.
    !! F/MI/PJK/LOGBOOK12, pp.67,89
    !! Bourque et al., Fusion Technology vol.21 (1992) 1465
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use fwbs_variables, only: emult, fhole, pnucblkt, pnucshld, pnucloss
    use heat_transport_variables, only: priheat, pthermmw, pfwdiv, nphx, pinjwp, &
      pinjht, crypmw, helpow
    use ife_variables, only: pdrive, ifetyp, etadrv, pifecr
    use physics_variables, only: powfmw

    implicit none

    !  Local variables

    real(dp) :: pdrvmw

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Driver power reaching target (MW)

    pdrvmw = 1.0D-6 * pdrive

    !  Primary nuclear heating (MW)
    !  Total thermal power removed from fusion core

    priheat = emult * powfmw

    !  Useful (high-grade) thermal power (MW)

    pthermmw = priheat * (1.0D0-fhole)

    !  Assume 0.24 of thermal power is intercepted by the first wall
    !  (Bourque et al)
    !  HYLIFE-II case: Assume FLiBe flows intercept all fusion power
    !  and provide the energy multiplication as though it were a
    !  conventional blanket

    if ((ifetyp /= 3).and.(ifetyp /= 4)) then
       pfwdiv = 0.24D0 * pthermmw
       pnucblkt = pthermmw - pfwdiv
    else
       pfwdiv = 0.0D0
       pnucblkt = pthermmw
    end if
    pnucshld = 0.0D0

    !  Lost fusion power (MW)

    pnucloss = priheat - pthermmw  !  = priheat*fhole

    !  Number of primary heat exchangers

    !rnphx = max(2.0D0, (pthermmw/400.0D0 + 0.8D0) )
    nphx = ceiling(pthermmw/1000.0D0)

    !  Secondary heat (some of it... rest calculated in IFEPW2)

    !  Wall plug driver power (MW)

    pinjwp = pdrvmw/etadrv

    !  Waste driver power (MW)

    pinjht = pinjwp - pdrvmw

    !  Cryogenic power (MW)
    !  Cryogenic temperature is assumed to be 4.5K

    crypmw = pifecr
    helpow = 1.0D6 * crypmw * (0.13D0 * 4.5D0)/(293.0D0 - 4.5D0)

  end subroutine ifepw1

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifeacp(outfile,iprint)

    !! Routine to calculate AC power requirements for an IFE power plant
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer : output file unit
    !! iprint : input integer : switch for writing to output file (1=yes)
    !! This routine calculates the AC power requirements for an IFE power plant.
    !! F/MI/PJK/LOGBOOK12, p.68
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use buildings_variables, only: efloor
    use process_output, only: oheadr, ovarre, oblnkl
    use heat_transport_variables, only: baseel, pwpm2, pacpmw, crypmw, vachtmw, &
      trithtmw, pinjwp, fcsht, tlvpmw
    use ife_variables, only: tdspmw, tfacmw, htpmw_ife, reprat, lipmw, ifetyp

    implicit none

    ! Arguments
    integer, intent(in) :: iprint,outfile

    !  Local variables

    real(dp), save :: basemw,pmwpm2

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Facility base load, MW (loads not dependent on floor area)

    basemw = baseel * 1.0D-6

    !  Power needed per floor area, MW/m2

    pmwpm2 = pwpm2 * 1.0D-6

    !  Total pulsed power system load, MW

    pacpmw = crypmw + vachtmw + tdspmw + tfacmw + &
         (htpmw_ife*reprat/6.0D0) + trithtmw + pinjwp + basemw + &
         (efloor*pmwpm2) + lipmw

    !  Total baseline power to facility loads, MW

    fcsht  = basemw + (efloor*pmwpm2)

    !  Estimate of the total low voltage power, MW

    tlvpmw = fcsht + trithtmw + (htpmw_ife*reprat/6.0D0) + vachtmw + &
         0.5D0*crypmw + tfacmw

    if (iprint == 0) return

    !  Output section
    call oheadr(outfile,'AC Power')

    call ovarre(outfile,'Facility base load (MW)','(basemw)',basemw)
    call ovarre(outfile,'Total floor space (m2)','(efloor)',efloor)
    call ovarre(outfile,'Power/floor area (MW/m2)','(pmwpm2)',pmwpm2)
    call ovarre(outfile,'Driver power supplies (MW)','(pinjwp)', &
         pinjwp)
    call ovarre(outfile,'Target delivery system (MW)','(tdspmw)', &
         tdspmw)
    call ovarre(outfile,'Target factory (MW)','(tfacmw)', &
         tfacmw)
    call ovarre(outfile,'Tritium processing plant (MW)','(trithtmw)', &
         trithtmw)
    call ovarre(outfile,'Vacuum pump motors (MW)','(vachtmw)',vachtmw)
    call ovarre(outfile,'Cryogenic comp motors (MW)','(crypmw)',crypmw)
    call ovarre(outfile,'Heat transport system pump motors (MW)', &
         '(htpmw_ife*reprat/6)',htpmw_ife*reprat/6.0D0)
    if (ifetyp.eq.4) then
        call ovarre(outfile,'Lithium Pump Power (MW)','(lipmw)',lipmw)
    end if
    call oblnkl(outfile)
    call ovarre(outfile,'Total pulsed power (MW)','(pacpmw)',pacpmw)
    call ovarre(outfile,'Total base power reqd at all times (MW)', &
         '(fcsht)',fcsht)
    call ovarre(outfile,'Total low voltage power (MW)','(tlvpmw)',tlvpmw)

  end subroutine ifeacp

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifepw2(outfile,iprint)

    !! Routine to calculate the rest of the IFE heat transport
    !! and plant power balance constituents, not already calculated in
    !! IFEPW1 or IFEACP
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer : output file unit
    !! iprint : input integer : switch for writing to output file (1=yes)
    !! This routine calculates the rest of the IFE heat transport
    !! and plant power balance constituents, not already calculated in
    !! routines <A HREF="ifepw1.html">IFEPW1</A> or
    !! <A HREF="ifeacp.html">IFEACP</A>.
    !! F/MI/PJK/LOGBOOK12, p.67
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use cost_variables, only: ireactor, ipnet
    use process_output, only: oheadr, osubhd, ovarre, oblnkl, ovarin
    use fwbs_variables, only: pnucloss, emult, tbr, pnucblkt
    use heat_transport_variables, only: fachtmw, fcsht, psechtmw, pinjht, &
      vachtmw, trithtmw, crypmw, pgrossmw, pthermmw, etath, fgrosbop, precircmw, &
      pacpmw, pnetelmw, pinjwp, pfwdiv, nphx
    use ife_variables, only: tdspmw, tfacmw, htpmw_ife, fauxbop, ifetyp, taufall

    implicit none

    !  Arguments
    integer, intent(in) :: outfile,iprint

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Facility heat removal (fcsht calculated in IFEACP)
    fachtmw = fcsht

    !  Total secondary heat
    psechtmw = pinjht + pnucloss + fachtmw + vachtmw + trithtmw + &
      tdspmw + tfacmw + crypmw + htpmw_ife

    !  Calculate powers relevant to a power-producing plant
    if (ireactor == 1) then

      !  Gross electric power
      pgrossmw = pthermmw * etath

      !  Balance of plant recirculating power fraction
      fgrosbop = min( 0.5D0, ( fauxbop/(pgrossmw/1000.0D0)**0.6D0) )

      !  Total recirculating power
      precircmw = (fgrosbop*pgrossmw) + pacpmw

      !  Net electric power
      pnetelmw = pgrossmw - precircmw

    end if

    if (iprint == 0) return

    !  Output section
    call oheadr(outfile,'Power / Heat Transport')
    call ovarre(outfile,'Fusion power escaping via holes (MW)', '(pnucloss)', pnucloss)
    call ovarre(outfile,'Power multiplication factor','(emult)',emult)
    if (ifetyp.eq.4) then
        call ovarre(outfile,'Tritium Breeding Ratio','(tbr)',tbr)
        call ovarre(outfile,'Lithium Fall Time (s)','(taufall)',taufall)
    end if
    call ovarre(outfile,'Driver wall plug power (MW)','(pinjwp)', pinjwp)
    call ovarre(outfile,'First wall nuclear heating (MW)','(pfwdiv)', pfwdiv)
    call ovarre(outfile,'Blanket nuclear heating (MW)','(pnucblkt)', pnucblkt)
    call ovarre(outfile,'Primary heat (MW)','(pthermmw)', pthermmw)
    call ovarre(outfile,'Secondary heat (MW)','(psechtmw)', psechtmw)
    call oblnkl(outfile)
    call ovarre(outfile,'Heat removal from driver power (MW)', '(pinjht)', pinjht)
    call ovarre(outfile,'Heat removal from cryogenic plant (MW)', '(crypmw)', crypmw)
    call ovarre(outfile,'Heat removal from vacuum pumps (MW)', '(vachtmw)', vachtmw)
    call ovarre(outfile,'Heat removal from target factory (MW)', '(tfacmw)', tfacmw)
    call ovarre(outfile,'Heat removal from delivery system (MW)', '(tdspmw)', tdspmw)
    call ovarre(outfile,'Heat removal from tritium plant (MW)', '(trithtmw)', trithtmw)
    call ovarre(outfile,'Heat removal from facilities (MW)', '(fachtmw)', fachtmw)
    call ovarin(outfile,'Number of primary heat exchangers', '(nphx)', nphx)

    if (ireactor /= 1) return
    call osubhd(outfile,'Reactor powers :')
    call ovarre(outfile,'Gross electric power (MW)', '(pgrossmw)', pgrossmw)
    call ovarre(outfile,'Net electric power (MW)', '(pnetelmw)', pnetelmw)
    call ovarre(outfile,'Balance of plant aux. power fraction', '(fgrosbop)', fgrosbop)

  end subroutine ifepw2

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifevac
    !! Routine to calculate parameters of the vacuum system for an
    !! Inertial Fusion Energy power plant
    !! author: P J Knight, CCFE, Culham Science Centre
    !! This routine calculates the parameters of the vacuum system for an
    !! Inertial Fusion Energy power plant.
    !! <P>The calculated values are hard-wired; they are based loosely
    !! on those for a tokamak of 6m major radius. F/MI/PJK/LOGBOOK12, p.87
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use vacuum_variables, only: dlscal, nvduct, vacdshm, vcdimax, vpumpn

    implicit none

    dlscal = 2.0D0
    nvduct = 16
    vacdshm = 0.0D0
    vcdimax = 0.3D0
    vpumpn = 32

  end subroutine ifevac

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ifebdg(outfile,iprint)
    !! Routine to calculate the volumes of the buildings required for
    !! an Inertial Fusion Energy power plant
    !! author: P J Knight, CCFE, Culham Science Centre
    !! outfile : input integer : output file unit
    !! iprint : input integer : switch for writing to output file (1=yes)
    !! This routine calculates the volumes of the buildings required for
    !! an Inertial Fusion Energy power plant. The method is based
    !! closely on that for tokamaks etc. in routine
    !! <A HREF="bldgs.html">BLDGS</A>.
    !! F/MI/PJK/LOGBOOK12, p.87
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use buildings_variables, only: wrbi, rbwt, rbrt, fndt, trcl, hcwt, hccl, &
      wgt2, stcl, pibv, efloor, triv, conv, admv, shov, admvol, convol, elevol, &
      rbvol, rmbvol, shovol, volrci, wsvol, volnucb
    use process_output, only: oheadr, ovarre
    use fwbs_variables, only: whtshld
    use heat_transport_variables, only: helpow
    use ife_variables, only: zl7, zu7, r7, zl6, zu6, r6

    implicit none

    !  Arguments
    integer, intent(in) :: outfile,iprint

    !  Local variables

    real(dp) :: cran,dcl,dcw,fac2,fac3,hcl,hcw,hrbi,rbh,rbl,rbw, &
         rmbh,rmbl,rmbw,rwl,rww,shh,tch,tcl,tcw,wgts,wsa
    real(dp), save :: cryv,elev,rbv,rmbv,vrci,wsv

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Reactor building
    !  ================

    !  Total internal height
    hrbi = zl7 + zu7

    !  Distance from centre of device to wall
    wrbi = r7

    !  Internal volume (square floor)
    vrci = (2.0D0 * wrbi)**2 * hrbi

    !  External dimensions
    !  RBWT = wall thickness
    !  RBRT = roof thickness
    !  FNDT = foundation thickness
    rbw = 2.0D0 * (r7 + rbwt)
    rbl = rbw
    rbh = hrbi + rbrt + fndt

    !  External volume
    rbv = rbw * rbl * rbh

    !  Maintenance building
    !  ====================

    !  The reactor maintenance building includes the hot cells, the
    !  decontamination chamber, the transfer corridors, and the waste
    !  treatment building.  The dimensions of these areas are scaled
    !  from a reference (tokamak) design based on the shield sector size.

    !  Shield height

    shh = zl6 + zu6

    !  Transport corridor size

    tcw = r6 + 4.0D0*trcl
    tcl = 5.0D0*tcw + 2.0D0*hcwt

    !  Decontamination cell size

    dcw = 2.0D0*tcw + 1.0D0
    dcl = 2.0D0*tcw + 1.0D0

    !  Hot cell size

    hcw = r6 + 3.0D0*hccl + 2.0D0
    hcl = 3.0D0*r6 + 4.0D0*hccl + tcw

    !  Radioactive waste treatment

    rww = dcw
    rwl = hcl-dcl-hcwt

    !  Maintenance building dimensions

    rmbw = hcw + dcw + 3.0D0*hcwt
    rmbl = hcl + 2.0D0*hcwt

    !  Height

    if (wgt2 > 1.0D0) then
       wgts = wgt2
    else
       wgts = whtshld
    end if
    cran = 9.41D-6*wgts + 5.1D0
    rmbh = 10.0D0 + (zl6+zu6) + trcl + cran + 5.1D0 + stcl + fndt
    tch = shh + stcl + fndt

    !  Volume

    fac2 = 2.8D0
    rmbv = fac2*rmbw*rmbl*rmbh + tcw*tcl*tch

    !  Warm shop and hot cell gallery

    wsa = (rmbw+7.0D0)*20.0D0 + rmbl*7.0D0
    fac3 = 1.9D0
    wsv = fac3 * wsa*rmbh

    !  Cryogenic building volume

    cryv = 55.0D0 * sqrt(helpow)

    !  Electrical building volume
    !  (set equal to power injection (i.e. driver) building volume)

    elev = pibv

    !  Calculate effective floor area for ac power module

    efloor = (rbv+rmbv+wsv+triv+elev+conv+cryv+admv+shov)/6.0D0

    !  Convert local into global variables

    admvol = admv
    convol = conv
    elevol = elev
    rbvol  = rbv
    rmbvol = rmbv
    shovol = shov
    volrci = vrci
    wsvol  = wsv

    !  Total volume of nuclear buildings

    volnucb = ( vrci + rmbv + wsv + triv + cryv )

    if (iprint == 0) return

    !  Output section

    call oheadr(outfile,'Plant Buildings System')
    call ovarre(outfile,'Internal volume of reactor building (m3)', &
         '(vrci)',vrci)
    call ovarre(outfile,'Dist from device centre to bldg wall (m)', &
         '(wrbi)',wrbi)
    call ovarre(outfile,'Effective floor area (m2)','(efloor)',efloor)
    call ovarre(outfile,'Reactor building volume (m3)','(rbv)',rbv)
    call ovarre(outfile,'Reactor maintenance building volume (m3)', &
         '(rmbv)',rmbv)
    call ovarre(outfile,'Warmshop volume (m3)','(wsv)',wsv)
    call ovarre(outfile,'Tritium building volume (m3)','(triv)',triv)
    call ovarre(outfile,'Electrical building volume (m3)','(elev)',elev)
    call ovarre(outfile,'Control building volume (m3)','(conv)',conv)
    call ovarre(outfile,'Cryogenics building volume (m3)','(cryv)',cryv)
    call ovarre(outfile,'Administration building volume (m3)','(admv)', &
         admv)
    call ovarre(outfile,'Shops volume (m3)','(shov)',shov)
    call ovarre(outfile,'Total volume of nuclear buildings (m3)', &
         '(volnucb)',volnucb)

  end subroutine ifebdg

end module ife_module
