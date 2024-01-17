module build_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the machine's radial and vertical build
  !!
  !!### References
  !!
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: aplasmin
  !! minimum minor radius (m)

  real(dp) :: available_radial_space
  !! Minimal radial space between plasma and coils (m)

  real(dp) :: blarea
  !! blanket total surface area (m2)

  real(dp) :: blareaib
  !! inboard blanket surface area (m2)

  real(dp) :: blareaob
  !! outboard blanket surface area (m2)

  real(dp) :: blbmith
  !! inboard blanket box manifold thickness (m) (`blktmodel>0`)
  !#TODO: remove blktmodel and similar below

  real(dp) :: blbmoth
  !! outboard blanket box manifold thickness (m) (`blktmodel>0`)

  real(dp) :: blbpith
  !! inboard blanket base plate thickness (m) (`blktmodel>0`)

  real(dp) :: blbpoth
  !! outboard blanket base plate thickness (m) (`blktmodel>0`)

  real(dp) :: blbuith
  !! inboard blanket breeding zone thickness (m) (`blktmodel>0`) (`iteration variable 90`)

  real(dp) :: blbuoth
  !! outboard blanket breeding zone thickness (m) (`blktmodel>0`) (`iteration variable 91`)

  real(dp) :: blnkith
  !! inboard blanket thickness (m); (calculated if `blktmodel>0`) (=0.0 if `iblnkith=0`)

  real(dp) :: blnkoth
  !! outboard blanket thickness (m); calculated if `blktmodel>0`

  real(dp) :: blnktth
  !! top blanket thickness (m), = mean of inboard and outboard blanket thicknesses

  real(dp) :: bore
  !! central solenoid inboard radius (m) (`iteration variable 29`)

  real(dp) :: clhsf
  !! cryostat lid height scaling factor (tokamaks)

  real(dp) :: ddwex
  !! cryostat thickness (m)

  real(dp) :: d_vv_in
  !! vacuum vessel inboard thickness (TF coil / shield) (m)

  real(dp) :: d_vv_out
  !! vacuum vessel outboard thickness (TF coil / shield) (m)

  real(dp) :: d_vv_top
  !! vacuum vessel topside thickness (TF coil / shield) (m) (= d_vv_bot if double-null)

  real(dp) :: d_vv_bot
  !! vacuum vessel underside thickness (TF coil / shield) (m)

  real(dp) :: f_avspace
  !! F-value for stellarator radial space check (`constraint equation 83`)

  real(dp) :: fcspc
  !! Fraction of space occupied by CS pre-compression structure

  real(dp) :: fseppc
  !! Separation force in CS coil pre-compression structure

  real(dp) :: fwarea
  !! first wall total surface area (m2)

  real(dp) :: fwareaib
  !! inboard first wall surface area (m2)

  real(dp) :: fwareaob
  !! outboard first wall surface area (m2)

  real(dp) :: fwith
  !! inboard first wall thickness, initial estimate as calculated (m)

  real(dp) :: fwoth
  !! outboard first wall thickness, initial estimate as calculated (m)

  real(dp) :: gapds
  !! gap between inboard vacuum vessel and thermal shield (m) (`iteration variable 61`)

  real(dp) :: gapoh
  !! gap between central solenoid and TF coil (m) (`iteration variable 42`)

  real(dp) :: gapomin
  !! minimum gap between outboard vacuum vessel and TF coil (m) (`iteration variable 31`)

  real(dp) :: gapsto
  !! gap between outboard vacuum vessel and TF coil (m)

  real(dp) :: hmax
  !! maximum (half-)height of TF coil (inside edge) (m)

  real(dp) :: hpfdif
  !! difference in distance from midplane of upper and lower portions of TF
  !! legs (non-zero for single-null devices) (m)

  real(dp) :: hpfu
  !! height to top of (upper) TF coil leg (m)

  real(dp) :: hr1
  !! half-height of TF coil inboard leg straight section (m)

  integer :: iohcl
  !! Switch for existence of central solenoid:
  !!
  !! - =0 central solenoid not present
  !! - =1 central solenoid exists

  integer :: iprecomp
  !! Switch for existence of central solenoid pre-compression structure:
  !!
  !! - =0 no pre-compression structure
  !! - =1 calculated pre-compression structure

  integer :: tf_in_cs
  !! Switch for placing the TF coil inside the CS
  !!
  !! - = 0 TF coil is outside the CS (default)
  !! - = 1 TF coil is inside the CS

  real(dp) :: ohcth
  !! Central solenoid thickness (m) (`iteration variable 16`)

  real(dp) :: precomp
  !! CS coil precompression structure thickness (m)

  real(dp) :: rbld
  !! sum of thicknesses to the major radius (m)

  real(dp) :: required_radial_space
  !! Required space between coil and plasma for blanket shield wall etc (m)

  real(dp) :: rinboard
  !! plasma inboard radius (m) (`consistency equation 29`)

  real(dp) :: rsldi
  !! radius to inboard shield (inside point) (m)

  real(dp) :: rsldo
  !! radius to outboard shield (outside point) (m)

  real(dp) :: r_vv_inboard_out
  !! Radial plasma facing side position of inboard vacuum vessel [m]

  real(dp) :: r_sh_inboard_in
  !! Radial inner side position of inboard neutronic shield [m]

  real(dp) :: r_sh_inboard_out
  !! Radial plasma facing side position of inboard neutronic shield [m]

  real(dp) :: r_tf_inboard_in
 	!! Mid-plane inboard TF coil leg radius at the centre-machine side [m]

  real(dp) :: r_tf_inboard_mid
  !! Mid-plane inboard TF coil leg radius at middle of the coil [m]

  real(dp) :: r_tf_inboard_out
  !! Mid-plane inboard TF coil leg radius at the plasma side [m]

  real(dp) :: r_tf_outboard_mid
  !! Mid-plane outboard TF coil leg radius at the middle of the coil [m]

  integer :: i_r_cp_top
  !! Switch selecting the he parametrization of the outer radius of the top of the CP part of the TF coil
  !!  0 : `r_cp_top` is set by the plasma shape
  !!  1 : `r_cp_top` is a user input
  !!  2 : `r_cp_top` is set using the CP top and midplane CP radius ratio

  real(dp) :: r_cp_top
  !! Top outer radius of the centropost (ST only) (m)

  real(dp) :: f_r_cp
  !! Ratio between the top and the midplane TF CP outer radius [-]
  !! Not used by default (-1) must be larger than 1 otherwise

  real(dp) :: dr_tf_inner_bore
  !! TF coil horizontal inner bore (m)

  real(dp) :: dh_tf_inner_bore
  !! TF coil vertical inner bore (m)

  real(dp) :: scrapli
  !! Gap between plasma and first wall, inboard side (m) (if `iscrp=1`)
  !! Iteration variable: ixc = 73
  !! Scan variable: nsweep = 58

  real(dp) :: scraplo
  !! Gap between plasma and first wall, outboard side (m) (if `iscrp=1`)
  !! Iteration variable: ixc = 74
  !! Scan variable: nsweep = 59

  real(dp) :: sharea
  !! shield total surface area (m2)

  real(dp) :: shareaib
  !! inboard shield surface area (m2)

  real(dp) :: shareaob
  !! outboard shield surface area (m2)

  real(dp) :: shldith
  !! inboard shield thickness (m) (`iteration variable 93`)

  real(dp) :: shldlth
  !! lower (under divertor) shield thickness (m)

  real(dp) :: shldoth
  !! outboard shield thickness (m) (`iteration variable 94`)

  real(dp) :: shldtth
  !! upper/lower shield thickness (m); calculated if `blktmodel > 0` (= shldlth if double-null)

  real(dp) :: sigallpc
  !! allowable stress in CSpre-compression structure (Pa)

  !#TODO: Issue #514 Make tfcth an output not an iteration variable
  real(dp) :: tfcth
  !! inboard TF coil thickness, (centrepost for ST) (m)
  !! (input, calculated or `iteration variable 13`)

  real(dp) :: tfoffset
  !! vertical distance between centre of TF coils and centre of plasma (m)

  real(dp) :: tfootfi
  !! TF coil outboard leg / inboard leg radial thickness
  !! ratio (`i_tf_sup=0` only) (`iteration variable 75`)

  real(dp) :: tfthko
  !! Outboard TF coil thickness (m)

  real(dp) :: tftsgap
  !! Minimum metal-to-metal gap between TF coil and thermal shield (m)

  real(dp) :: thshield_ib
  !! TF-VV thermal shield thickness, inboard (m)

  real(dp) :: thshield_ob
  !! TF-VV thermal shield thickness, outboard (m)

  real(dp) :: thshield_vb
  !! TF-VV thermal shield thickness, vertical build (m)

  real(dp) :: vgap2
  !! vertical gap between vacuum vessel and thermal shields (m)

  real(dp) :: vgap
  !! vertical gap between x-point and divertor (m) (if = 0, it is calculated)

  real(dp) :: vgaptop
  !! vertical gap between top of plasma and first wall (m) (= vgap if double-null)

  real(dp) :: vvblgap
  !! gap between vacuum vessel and blanket (m)

  real(dp) :: plleni
  !! length of inboard divertor plate (m)

  real(dp) :: plleno
  !! length of outboard divertor plate (m)

  real(dp) :: plsepi
  !! poloidal length, x-point to inboard strike point (m)

  real(dp) :: plsepo
  !! poloidal length, x-point to outboard strike point (m)

  real(dp) :: rspo
  !! outboard strike point radius (m)

  contains

  subroutine init_build_variables
    !! Initialise module variables
    implicit none

    aplasmin = 0.25D0
    available_radial_space = 0.0D0
    blarea = 0.0D0
    blareaib = 0.0D0
    blareaob = 0.0D0
    blbmith = 0.17D0
    blbmoth = 0.27D0
    blbpith = 0.30D0
    blbpoth = 0.35D0
    blbuith = 0.365D0
    blbuoth = 0.465D0
    blnkith = 0.115D0
    blnkoth = 0.235D0
    blnktth = 0.0D0
    bore = 1.42D0
    clhsf = 4.268D0
    ddwex = 0.07D0
    d_vv_in = 0.07D0
    d_vv_out = 0.07D0
    d_vv_top = 0.07D0
    d_vv_bot = 0.07D0
    f_avspace = 1.0D0
    fcspc = 0.6D0
    fseppc = 3.5D8
    fwarea = 0.0D0
    fwareaib = 0.0D0
    fwareaob = 0.0D0
    fwith = 0.0D0
    fwoth = 0.0D0
    gapds = 0.155D0
    gapoh = 0.08D0
    gapomin = 0.234D0
    gapsto = 0.0D0
    hmax = 0.0D0
    hpfdif = 0.0D0
    hpfu = 0.0D0
    hr1 = 0.0D0
    iohcl = 1
    iprecomp = 1
    tf_in_cs = 0
    ohcth = 0.811D0
    precomp = 0.0D0
    rbld = 0.0D0
    required_radial_space = 0.0D0
    rinboard = 0.651D0
    rsldi = 0.0D0
    rsldo = 0.0D0
    r_vv_inboard_out = 0.0D0
    r_sh_inboard_out = 0.0D0
    r_tf_inboard_in = 0.0D0
    r_tf_inboard_mid = 0.0D0
    r_tf_inboard_out = 0.0D0
    r_tf_outboard_mid = 0.0D0
    i_r_cp_top = 0
    r_cp_top = 0.0D0
    f_r_cp = 1.4D0
    dr_tf_inner_bore = 0.0D0
    dh_tf_inner_bore = 0.0D0
    scrapli = 0.14D0
    scraplo = 0.15D0
    sharea = 0.0D0
    shareaib = 0.0D0
    shareaob = 0.0D0
    shldith = 0.69D0
    shldlth = 0.7D0
    shldoth = 1.05D0
    shldtth = 0.6D0
    sigallpc = 3.0D8
    tfcth = 0.0D0
    tfoffset = 0.0D0
    tfootfi = 1.19D0
    tfthko = 0.0D0
    tftsgap = 0.05D0
    thshield_ib = 0.05D0
    thshield_ob = 0.05D0
    thshield_vb = 0.05D0
    vgap2 = 0.163D0
    vgap= 0.0D0
    vgaptop = 0.60D0
    vvblgap = 0.05D0
    plleni = 1.0D0
    plleno = 1.0D0
    plsepi = 1.0D0
    plsepo = 1.5D0
    rspo = 0.0D0
    r_sh_inboard_in = 0.0D0
  end subroutine init_build_variables
end module build_variables
