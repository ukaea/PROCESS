module ife_variables
    !! author: S. Muldrew (UKAEA)
    !!
    !! Module containing global variables relating to the inertial fusion energy model
    !!
    !!### References
    !!
    !! -
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
  end module ife_variables
