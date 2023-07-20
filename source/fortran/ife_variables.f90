module ife_variables
    !! author: S. Muldrew (UKAEA)
    !!
    !! Module containing global variables relating to the inertial fusion energy model
    !!
    !!### References
    !!
    !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    public


    !! Default IFE builds and material volumes are those for the SOMBRERO device.
    !! The 2-dimensional arrays have indices (region, material), where 'region'
    !! is the region and maxmat is the 'material':
    !!
    !! - 'region' = 1 radially outside chamber
    !! - 'region' = 2 above chamber
    !! - 'region' = 3 below chamber

    integer, parameter ::  maxmat = 8
    !! Total number of materials in IFE device. Material numbers are as follows:
    !!
    !! - =0 void
    !! - =1 steel
    !! - =2 carbon cloth
    !! - =3 FLiBe
    !! - =4 lithium oxide Li2O
    !! - =5 concrete
    !! - =6 helium
    !! - =7 xenon
    !! - =8 lithium

    real(dp) :: bldr
    !! radial thickness of IFE blanket (m; calculated `if ifetyp=4`)

    real(dp) :: bldrc
    !! radial thickness of IFE curtain (m; `ifetyp=4`)

    real(dp) :: bldzl
    !! vertical thickness of IFE blanket below chamber (m)

    real(dp) :: bldzu
    !! vertical thickness of IFE blanket above chamber (m)

    real(dp), dimension(3,0:maxmat) :: blmatf
    !! IFE blanket material fractions

    real(dp), dimension(3,0:maxmat) :: blmatm
    !! IFE blanket material masses (kg)

    real(dp), dimension(3,0:maxmat) :: blmatv
    !! IFE blanket material volumes (m3)

    real(dp), dimension(3) :: blvol
    !! IFE blanket volume (m3)

    real(dp) :: cdriv0
    !! IFE generic/laser driver cost at edrive=0 (M$)

    real(dp) :: cdriv1
    !! IFE low energy heavy ion beam driver cost extrapolated to `edrive=0` (M$)

    real(dp) :: cdriv2
    !! IFE high energy heavy ion beam driver cost extrapolated to `edrive=0` (M$)

    real(dp) :: cdriv3
    !! IFE driver cost ($/J wall plug) (`ifedrv==3`)

    real(dp) :: chdzl
    !! vertical thickness of IFE chamber below centre (m)

    real(dp) :: chdzu
    !! vertical thickness of IFE chamber above centre (m)

    real(dp), dimension(0:maxmat) :: chmatf
    !! IFE chamber material fractions

    real(dp), dimension(0:maxmat) :: chmatm
    !! IFE chamber material masses (kg)

    real(dp), dimension(0:maxmat) :: chmatv
    !! IFE chamber material volumes (m3)

    real(dp) :: chrad
    !! radius of IFE chamber (m) (`iteration variable 84`)

    real(dp) :: chvol
    !! IFE chamber volume (m3)

    real(dp) :: dcdrv0
    !! IFE generic/laser driver cost gradient (M$/MJ)

    real(dp) :: dcdrv1
    !! HIB driver cost gradient at low energy (M$/MJ)

    real(dp) :: dcdrv2
    !! HIB driver cost gradient at high energy (M$/MJ)

    real(dp) :: drveff
    !! IFE driver wall plug to target efficiency (`ifedrv=0,3`) (`iteration variable 82`)

    real(dp) :: edrive
    !! IFE driver energy (J) (`iteration variable 81`)

    real(dp) :: etadrv
    !! IFE driver wall plug to target efficiency

    real(dp) :: etali
    !! IFE lithium pump wall plug efficiency (`ifetyp=4`)

    real(dp), dimension(10) :: etave
    !! IFE driver efficiency vs driver energy (`ifedrv=-1`)

    real(dp) :: fauxbop
    !! fraction of gross electric power to balance-of-plant (IFE)

    real(dp) :: fbreed
    !! fraction of breeder external to device core

    real(dp) :: fburn
    !! IFE burn fraction (fraction of tritium fused/target)

    real(dp) :: flirad
    !! radius of FLiBe/lithium inlet (m) (`ifetyp=3,4`)

    real(dp) :: frrmax
    !! f-value for maximum IFE repetition rate (`constraint equation 50`, `iteration variable 86`)

    real(dp) :: fwdr
    !! radial thickness of IFE first wall (m)

    real(dp) :: fwdzl
    !! vertical thickness of IFE first wall below chamber (m)

    real(dp) :: fwdzu
    !! vertical thickness of IFE first wall above chamber (m)

    real(dp), dimension(3,0:maxmat) :: fwmatf
    !! IFE first wall material fractions

    real(dp), dimension(3,0:maxmat) :: fwmatm
    !! IFE first wall material masses (kg)

    real(dp), dimension(3,0:maxmat) :: fwmatv
    !! IFE first wall material volumes (kg)

    real(dp), dimension(3) :: fwvol
    !! IFE first wall volume (m3)

    real(dp) :: gain
    !! IFE target gain

    real(dp), dimension(10) :: gainve
    !! IFE target gain vs driver energy (`ifedrv=-1`)

    real(dp) :: htpmw_ife
    !! IFE heat transport system electrical pump power (MW)

    integer :: ife
    !! Switch for IFE option:
    !!
    !! - =0 use tokamak, RFP or stellarator model
    !! - =1 use IFE model

    integer :: ifedrv
    !! Switch for type of IFE driver:
    !!
    !! - =-1 use gainve, etave for gain and driver efficiency
    !! - =0 use tgain, drveff for gain and driver efficiency
    !! - =1 use laser driver based on SOMBRERO design
    !! - =2 use heavy ion beam driver based on OSIRIS
    !! - =3 Input pfusife, rrin and drveff

    integer :: ifetyp
    !! Switch for type of IFE device build:
    !!
    !! - =0 generic (cylindrical) build
    !! - =1 OSIRIS-like build
    !! - =2 SOMBRERO-like build
    !! - =3 HYLIFE-II-like build
    !! - =4 2019 build

    real(dp) :: lipmw
    !! IFE lithium pump power (MW; `ifetyp=4`)

    real(dp) :: mcdriv
    !! IFE driver cost multiplier

    real(dp) :: mflibe
    !! total mass of FLiBe (kg)

    real(dp) :: pdrive
    !! IFE driver power reaching target (W) (`iteration variable 85`)

    real(dp) :: pfusife
    !! IFE input fusion power (MW) (`ifedrv=3 only`; `itv 155`)

    real(dp) :: pifecr
    !! IFE cryogenic power requirements (MW)

    real(dp) :: ptargf
    !! IFE target factory power at 6 Hz repetition rate (MW)

    real(dp) :: r1
    !! IFE device radial build (m)

    real(dp) :: r2
    !! IFE device radial build (m)

    real(dp) :: r3
    !! IFE device radial build (m)

    real(dp) :: r4
    !! IFE device radial build (m)

    real(dp) :: r5
    !! IFE device radial build (m)

    real(dp) :: r6
    !! IFE device radial build (m)

    real(dp) :: r7
    !! IFE device radial build (m)

    real(dp) :: reprat
    !! IFE driver repetition rate (Hz)

    real(dp) :: rrin
    !! Input IFE repetition rate (Hz) (`ifedrv=3 only`; `itv 156`)

    real(dp) :: rrmax
    !! maximum IFE repetition rate (Hz)

    real(dp) :: shdr
    !! radial thickness of IFE shield (m)

    real(dp) :: shdzl
    !! vertical thickness of IFE shield below chamber (m)

    real(dp) :: shdzu
    !! vertical thickness of IFE shield above chamber (m)

    real(dp), dimension(3,0:maxmat) :: shmatf
    !! IFE shield material fractions

    real(dp), dimension(3,0:maxmat) :: shmatm
    !! IFE shield material masses (kg)

    real(dp), dimension(3,0:maxmat) :: shmatv
    !! IFE shield material volumes (kg)

    real(dp), dimension(3) :: shvol
    !! IFE shield volume (m3)

    real(dp) :: sombdr
    !! radius of cylindrical blanket section below chamber (`ifetyp=2`)

    real(dp) :: somtdr
    !! radius of cylindrical blanket section above chamber (`ifetyp=2`)

    real(dp) :: taufall
    !! Lithium Fall Time (s)

    real(dp) :: tdspmw
    !! IFE target delivery system power (MW)

    real(dp) :: tfacmw
    !! IFE target factory power (MW)

    real(dp) :: tgain
    !! IFE target gain (if `ifedrv = 0`) (`iteration variable 83`)

    real(dp) :: uccarb
    !! cost of carbon cloth ($/kg)

    real(dp) :: ucconc
    !! cost of concrete ($/kg)

    real(dp) :: ucflib
    !! cost of FLiBe ($/kg)

    real(dp) :: uctarg
    !! cost of IFE target ($/target)

    real(dp) :: v1dr
    !! radial thickness of IFE void between first wall and blanket (m)

    real(dp) :: v1dzl
    !! vertical thickness of IFE void 1 below chamber (m)

    real(dp) :: v1dzu
    !! vertical thickness of IFE void 1 above chamber (m)

    real(dp), dimension(3,0:maxmat) :: v1matf
    !! IFE void 1 material fractions

    real(dp), dimension(3,0:maxmat) :: v1matm
    !! IFE void 1 material masses (kg)

    real(dp), dimension(3,0:maxmat) :: v1matv
    !! IFE void 1 material volumes (kg)

    real(dp), dimension(3) :: v1vol
    !! IFE void 1 volume (m3)

    real(dp) :: v2dr
    !! radial thickness of IFE void between blanket and shield (m)

    real(dp) :: v2dzl
    !! vertical thickness of IFE void 2 below chamber (m)

    real(dp) :: v2dzu
    !! vertical thickness of IFE void 2 above chamber (m)

    real(dp), dimension(3,0:maxmat) :: v2matf
    !! IFE void 2 material fractions

    real(dp), dimension(3,0:maxmat) :: v2matm
    !! IFE void 2 material masses (kg)

    real(dp), dimension(3,0:maxmat) :: v2matv
    !! IFE void 2 material volumes (kg)

    real(dp), dimension(3) :: v2vol
    !! IFE void 2 volume (m3)

    real(dp) :: v3dr
    !! radial thickness of IFE void outside shield (m)

    real(dp) :: v3dzl
    !! vertical thickness of IFE void 3 below chamber (m)

    real(dp) :: v3dzu
    !! vertical thickness of IFE void 3 above chamber (m)

    real(dp), dimension(3,0:maxmat) :: v3matf
    !! IFE void 3 material fractions

    real(dp), dimension(3,0:maxmat) :: v3matm
    !! IFE void 3 material masses (kg)

    real(dp), dimension(3,0:maxmat) :: v3matv
    !! IFE void 3 material volumes (kg)

    real(dp), dimension(3) :: v3vol
    !! IFE void 3 volume (m3)

    real(dp) :: zl1
    !! IFE vertical build below centre (m)

    real(dp) :: zl2
    !! IFE vertical build below centre (m)

    real(dp) :: zl3
    !! IFE vertical build below centre (m)

    real(dp) :: zl4
    !! IFE vertical build below centre (m)

    real(dp) :: zl5
    !! IFE vertical build below centre (m)

    real(dp) :: zl6
    !! IFE vertical build below centre (m)

    real(dp) :: zl7
    !! IFE vertical build below centre (m)

    real(dp) :: zu1
    !! IFE vertical build above centre (m)

    real(dp) :: zu2
    !! IFE vertical build above centre (m)

    real(dp) :: zu3
    !! IFE vertical build above centre (m)

    real(dp) :: zu4
    !! IFE vertical build above centre (m)

    real(dp) :: zu5
    !! IFE vertical build above centre (m)

    real(dp) :: zu6
    !! IFE vertical build above centre (m)

    real(dp) :: zu7
    !! IFE vertical build above centre (m)

    contains

    subroutine init_ife_variables
      !! Initialise module variables
      implicit none

      bldr   = 1.0D0
      bldrc   = 1.0D0
      bldzl  = 4.0D0
      bldzu  = 4.0D0
      blmatf = reshape( (/ &
        0.05D0,0.05D0,0.05D0, &
        0.0D0, 0.0D0, 0.0D0,  &
        0.45D0,0.45D0,0.45D0, &
        0.0D0, 0.0D0, 0.0D0,  &
        0.20D0,0.20D0,0.20D0, &
        0.0D0, 0.0D0, 0.0D0,  &
        0.30D0,0.30D0,0.30D0, &
        0.0D0,0.0D0,0.0D0,    &
        0.0D0, 0.0D0, 0.0D0  /), shape(blmatf))
      blmatm = 0.0D0
      blmatv = 0.0D0
      blvol = 0.0D0
      cdriv0 = 154.3D0
      cdriv1 = 163.2D0
      cdriv2 = 244.9D0
      cdriv3 = 1.463D0
      chdzl = 9.0D0
      chdzu = 9.0D0
      chmatf = &
        (/1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/)
      chmatm = 0.0D0
      chmatv = 0.0D0
      chrad = 6.5D0
      chvol = 0.0D0
      dcdrv0 = 111.4D0
      dcdrv1 = 78.0D0
      dcdrv2 = 59.9D0
      drveff = 0.28D0
      edrive = 5.0D6
      etadrv = 0.0D0
      etali = 0.4D0
      etave = (/ &
        0.082D0,0.079D0,0.076D0,0.073D0,0.069D0, &
        0.066D0,0.062D0,0.059D0,0.055D0,0.051D0 /)
      fauxbop = 0.06D0
      fbreed = 0.51D0
      fburn  = 0.3333D0
      flirad = 0.78D0
      frrmax = 1.0D0
      fwdr = 0.01D0
      fwdzl = 0.01D0
      fwdzu = 0.01D0
      fwmatf = reshape( (/ &
        0.05D0,0.05D0,0.05D0, &
        0.0D0, 0.0D0, 0.0D0,  &
        0.95D0,0.95D0,0.95D0, &
        0.0D0, 0.0D0, 0.0D0,  &
        0.0D0, 0.0D0, 0.0D0,  &
        0.0D0, 0.0D0, 0.0D0,  &
        0.0D0, 0.0D0, 0.0D0,  &
        0.0D0, 0.0D0, 0.0D0,  &
        0.0D0, 0.0D0, 0.0D0  /), shape(fwmatf))
      fwmatm = 0.0D0
      fwmatv = 0.0D0
      fwvol = 0.0D0
      gain = 0.0D0
      gainve = (/ &
        60.0D0, 95.0D0,115.0D0,125.0D0,133.0D0, &
        141.0D0,152.0D0,160.0D0,165.0D0,170.0D0 /)
      htpmw_ife = 0.0D0
      ife = 0
      ifedrv = 2
      ifetyp = 0
      lipmw = 0.0D0
      mcdriv = 1.0D0
      mflibe = 0.0D0
      pdrive = 23.0D6
      pfusife = 1000.0D0
      pifecr = 10.0D0
      ptargf = 2.0D0
      r1 = 0.0D0
      r2 = 0.0D0
      r3 = 0.0D0
      r4 = 0.0D0
      r5 = 0.0D0
      r6 = 0.0D0
      r7 = 0.0D0
      reprat = 0.0D0
      rrin = 6.0D0
      rrmax = 20.0D0
      shdr = 1.7D0
      shdzl = 5.0D0
      shdzu  = 5.0D0
      shmatf = reshape( (/ &
        0.05D0,0.05D0,0.05D0, &
        0.19D0,0.19D0,0.19D0, &
        0.0D0, 0.0D0, 0.0D0,  &
        0.0D0, 0.0D0, 0.0D0,  &
        0.0D0, 0.0D0, 0.0D0,  &
        0.665D0,0.665D0,0.665D0, &
        0.095D0,0.095D0,0.095D0, &
        0.0D0, 0.0D0, 0.0D0,  &
        0.0D0, 0.0D0, 0.0D0  /), shape(shmatf))
      shmatm = 0.0D0
      shmatv = 0.0D0
      shvol = 0.0D0
      sombdr = 2.7D0
      somtdr = 2.7D0
      taufall = 0.0D0
      tdspmw = 0.01D0
      tfacmw = 0.0D0
      tgain = 85.0D0
      uccarb = 50.0D0
      ucconc = 0.1D0
      ucflib = 84.0D0
      uctarg = 0.3D0
      v1dr = 0.0D0
      v1dzl = 0.0D0
      v1dzu = 0.0D0
      v1matf = reshape( (/ &
        1.0D0, 1.0D0, 1.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0  /), shape(v1matf))
      v1matm = 0.0D0
      v1matv = 0.0D0
      v1vol = 0.0D0
      v2dr = 2.0D0
      v2dzl = 7.0D0
      v2dzu = 7.0D0
      v2matf = reshape( (/ &
        1.0D0, 1.0D0, 1.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0  /), shape(v2matf))
      v2matm = 0.0D0
      v2matv = 0.0D0
      v2vol = 0.0D0
      v3dr   = 43.3D0
      v3dzl  = 30.0D0
      v3dzu  = 20.0D0
      v3matf = reshape( (/ &
        1.0D0, 1.0D0, 1.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0, &
        0.0D0, 0.0D0, 0.0D0  /), shape(v3matf))
      v3matm = 0.0D0
      v3matv = 0.0D0
      v3vol = 0.0D0
      zl1 = 0.0D0
      zl2 = 0.0D0
      zl3 = 0.0D0
      zl4 = 0.0D0
      zl5 = 0.0D0
      zl6 = 0.0D0
      zl7 = 0.0D0
      zu1 = 0.0D0
      zu2 = 0.0D0
      zu3 = 0.0D0
      zu4 = 0.0D0
      zu5 = 0.0D0
      zu6 = 0.0D0
      zu7 = 0.0D0
    end subroutine init_ife_variables
  end module ife_variables
