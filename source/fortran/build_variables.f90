module build_variables
  !! author: J. Morris (UKAEA)
  !!
  !! Module containing global variables relating to the machine's radial and vertical build
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: aplasmin
  !! minimum minor radius (m)

  real(dp) :: available_radial_space
  !! Minimal radial space between plasma and coils (m)

  real(dp) :: a_blkt_total_surface
  !! blanket total surface area (m2)

  real(dp) :: a_blkt_inboard_surface
  !! inboard blanket surface area (m2)

  real(dp) :: a_blkt_outboard_surface
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

  real(dp) :: dr_blkt_inboard
  !! inboard blanket thickness (m); (calculated if `blktmodel>0`) (=0.0 if `i_blkt_inboard=0`)

  real(dp) :: dr_blkt_outboard
  !! outboard blanket thickness (m); calculated if `blktmodel>0`

  real(dp) :: dz_blkt_upper
  !! top blanket thickness (m), = mean of inboard and outboard blanket thicknesses

  real(dp) :: dz_fw_upper
  !! upper first wall thickness (m)

  real(dp) :: dr_bore
  !! central solenoid inboard radius (m) (`iteration variable 29`)

  real(dp) :: f_z_cryostat
  !! cryostat lid height scaling factor (tokamaks)

  real(dp) :: dr_cryostat
  !! cryostat thickness (m)

  real(dp) :: dr_vv_inboard
  !! vacuum vessel inboard thickness (TF coil / shield) (m)

  real(dp) :: dr_vv_outboard
  !! vacuum vessel outboard thickness (TF coil / shield) (m)

  real(dp) :: dz_vv_upper
  !! vacuum vessel topside thickness (TF coil / shield) (m) (= dz_vv_lower if double-null)

  real(dp) :: dz_vv_lower
  !! vacuum vessel underside thickness (TF coil / shield) (m)

  real(dp) :: dr_vv_shells
  !! vacuum vessel double walled shell thicknesses (m)

  real(dp) :: f_avspace
  !! F-value for stellarator radial space check (`constraint equation 83`)

  real(dp) :: fcspc
  !! Fraction of space occupied by CS pre-compression structure

  real(dp) :: fseppc
  !! Separation force in CS coil pre-compression structure

  real(dp) :: a_fw_total
  !! First wall total surface area [m^2]

  real(dp) :: a_fw_inboard
  !! Inboard first wall surface area [m^2]

  real(dp) :: a_fw_outboard
  !! Outboard first wall surface area [m^2]

  real(dp) :: dr_fw_inboard
  !! inboard first wall thickness, initial estimate as calculated (m)

  real(dp) :: dr_fw_outboard
  !! outboard first wall thickness, initial estimate as calculated (m)

  real(dp) :: dr_shld_vv_gap_inboard
  !! gap between inboard vacuum vessel and thermal shield (m) (`iteration variable 61`)

  real(dp) :: dr_cs_tf_gap
  !! gap between central solenoid and TF coil (m) (`iteration variable 42`)

  real(dp) :: gapomin
  !! minimum gap between outboard vacuum vessel and TF coil (m) (`iteration variable 31`)

  real(dp) :: dr_shld_vv_gap_outboard
  !! gap between outboard vacuum vessel and TF coil (m)

  real(dp) :: z_tf_inside_half
  !! maximum (half-)height of TF coil (inside edge) (m)

  real(dp) :: hpfdif
  !! difference in distance from midplane of upper and lower portions of TF
  !! legs (non-zero for single-null devices) (m)

  real(dp) :: z_tf_top
  !! height to top of (upper) TF coil leg (m)

  real(dp) :: hr1
  !! half-height of TF coil inboard leg straight section (m)

  integer :: iohcl
  !! Switch for existence of central solenoid:
  !!
  !! - =0 central solenoid not present
  !! - =1 central solenoid exists

  integer :: i_cs_precomp
  !! Switch for existence of central solenoid pre-compression structure:
  !!
  !! - =0 no pre-compression structure
  !! - =1 calculated pre-compression structure

  integer :: i_tf_inside_cs
  !! Switch for placing the TF coil inside the CS
  !!
  !! - = 0 TF coil is outside the CS (default)
  !! - = 1 TF coil is inside the CS

  real(dp) :: dr_cs
  !! Central solenoid thickness (m) (`iteration variable 16`)

  real(dp) :: dr_cs_precomp
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
  !! TF coil horizontal inner dr_bore (m)

  real(dp) :: dh_tf_inner_bore
  !! TF coil vertical inner dr_bore (m)

  real(dp) :: dr_fw_plasma_gap_inboard
  !! Gap between plasma and first wall, inboard side (m) (if `i_plasma_wall_gap=1`)
  !! Iteration variable: ixc = 73
  !! Scan variable: nsweep = 58

  real(dp) :: dr_fw_plasma_gap_outboard
  !! Gap between plasma and first wall, outboard side (m) (if `i_plasma_wall_gap=1`)
  !! Iteration variable: ixc = 74
  !! Scan variable: nsweep = 59

  real(dp) :: a_shld_total_surface
  !! shield total surface area (m2)

  real(dp) :: a_shld_inboard_surface
  !! inboard shield surface area (m2)

  real(dp) :: a_shld_outboard_surface
  !! outboard shield surface area (m2)

  real(dp) :: dr_shld_inboard
  !! inboard shield thickness (m) (`iteration variable 93`)

  real(dp) :: dz_shld_lower
  !! lower (under divertor) shield thickness (m)

  real(dp) :: dr_shld_outboard
  !! outboard shield thickness (m) (`iteration variable 94`)

  real(dp) :: dz_shld_upper
  !! upper/lower shield thickness (m); calculated if `blktmodel > 0` (= dz_shld_lower if double-null)

  real(dp) :: sigallpc
  !! allowable stress in CSpre-compression structure (Pa)

  !#TODO: Issue #514 Make dr_tf_inboard an output not an iteration variable
  real(dp) :: dr_tf_inboard
  !! inboard TF coil thickness, (centrepost for ST) (m)
  !! (input, calculated or `iteration variable 13`)

  real(dp) :: tfoffset
  !! vertical distance between centre of TF coils and centre of plasma (m)

  real(dp) :: tfootfi
  !! TF coil outboard leg / inboard leg radial thickness
  !! ratio (`i_tf_sup=0` only) (`iteration variable 75`)

  real(dp) :: dr_tf_outboard
  !! Outboard TF coil thickness (m)

  real(dp) :: dr_tf_shld_gap
  !! Minimum metal-to-metal gap between TF coil and thermal shield (m)

  real(dp) :: dr_shld_thermal_inboard
  !! TF-VV thermal shield thickness, inboard (m)

  real(dp) :: dr_shld_thermal_outboard
  !! TF-VV thermal shield thickness, outboard (m)

  real(dp) :: dz_shld_thermal
  !! TF-VV thermal shield thickness, vertical build (m)

  real(dp) :: dz_shld_vv_gap
  !! vertical gap between vacuum vessel and thermal shields (m)

  real(dp) :: dz_xpoint_divertor
  !! vertical gap between x-point and divertor (m) (if = 0, it is calculated)

  real(dp) :: dz_fw_plasma_gap
  !! vertical gap between top of plasma and first wall (m) (= dz_xpoint_divertor if double-null)

  real(dp) :: dr_shld_blkt_gap
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

  real(dp) :: z_plasma_xpoint_upper
  !! Vertical height of the upper plasma x-point (m)

  real(dp) :: z_plasma_xpoint_lower
  !! Vertical height of the lower plasma x-point (m)
end module build_variables
