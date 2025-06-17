module pfcoil_variables
  !! author: J. Morris, M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to the poloidal field coil systems
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  integer, parameter :: n_pf_groups_max = 10
  !! maximum number of groups of PF coils

  integer, parameter :: n_pf_coils_in_group_max = 2
  !! maximum number of PF coils in a given group

  integer, parameter :: nptsmx = 32
  !! maximum number of points across the midplane of the plasma at which the field from
  !! the PF coils is fixed

  integer, parameter :: nfixmx = 64
  !! maximum number of fixed current PF coils

  integer, parameter :: ngc = n_pf_groups_max*n_pf_coils_in_group_max
  !! maximum total number of coils across all groups

  integer, parameter :: ngc2 = ngc+2
  !! new variable to include 2 additional circuits: plasma and central solenoid

  real(dp) :: alfapf
  !! smoothing parameter used in PF coil current calculation at the beginning of pulse (BoP)

  real(dp) :: alstroh
  !! allowable hoop stress in Central Solenoid structural material (Pa)

  integer :: i_cs_stress
  !! Switch for CS stress calculation:
  !!
  !! - =0 Hoop stress only
  !! - =1 Hoop + Axial stress

  real(dp) :: a_cs_poloidal
  !! Central solenoid vertical cross-sectional area (m2)

  real(dp) :: a_cs_turn
  !! Central solenoid (OH) trun cross-sectional area (m2)

  real(dp) :: awpoh
  !! central solenoid conductor+void area with area of steel subtracted (m2)

  real(dp) :: b_cs_peak_flat_top_end
  !! maximum field in central solenoid at end of flat-top (EoF) (T)

  real(dp) :: b_cs_peak_pulse_start
  !! maximum field in central solenoid at beginning of pulse (T)

  real(dp), dimension(ngc2) :: b_pf_coil_peak
  !! peak field at coil i (T)

  real(dp), dimension(n_pf_groups_max) :: ccl0_ma
  !! PF group current array, flux-swing cancellation current (MA)
  !! Input if i_pf_current=0, computed otherwise

  real(dp), dimension(n_pf_groups_max) :: ccls_ma
  !! PF group current array, equilibrium current (MA)
  !! Input if i_pf_current=0, computed otherwise

  real(dp) :: j_cs_pulse_start
  !! Central solenoid overall current density at beginning of pulse (A/m2)

  real(dp) :: j_cs_flat_top_end
  !! Central solenoid overall current density at end of flat-top (A/m2) (`iteration variable 37`) (`sweep variable 62`)

  real(dp), dimension(ngc2,6) :: c_pf_coil_turn
  !! current per turn in coil i at time j (A)

  real(dp), dimension(ngc2) :: c_pf_coil_turn_peak_input
  !! peak current per turn input for PF coil i (A)

  real(dp), dimension(ngc2) :: c_pf_cs_coil_pulse_start_ma
  !! PF coil current array, at beginning of pulse (MA)
  !! Indexed by coil number, not group number

  real(dp), dimension(ngc2) :: c_pf_cs_coil_flat_top_ma
  !! PF coil current array, at flat top (MA)
  !! Indexed by coil number, not group number

  real(dp), dimension(ngc2) :: c_pf_cs_coil_pulse_end_ma
  !! PF coil current array, at end of pulse (MA)
  !! Indexed by coil number, not group number

  real(dp) :: etapsu
  !! Efficiency of transfer of PF stored energy into or out of storage.

  real(dp) :: f_j_cs_start_end_flat_top
  !! ratio of central solenoid overall current density at beginning of flat-top / end of flat-top

  real(dp) :: f_j_cs_start_pulse_end_flat_top
  !! ratio of central solenoid overall current density at beginning of pulse / end of flat-top
  !! (`iteration variable 41`)

  real(dp) :: fcuohsu
  !! copper fraction of strand in central solenoid

  real(dp) :: fcupfsu
  !! copper fraction of cable conductor (PF coils)

  real(dp) :: fvs_cs_pf_total_ramp
  !! F-value for `constraint equation 51`

  integer, dimension(n_pf_groups_max) :: i_pf_location
  !! Switch for location of PF coil group i:
  !!
  !! - =1 PF coil on top of central solenoid (flux ramp only)
  !! - =2 PF coil on top of TF coil (flux ramp only)
  !! - =3 PF coil outside of TF coil (equilibrium coil)
  !! - =4 PF coil, general location (equilibrium coil)

  integer :: i_pf_conductor
  !! switch for PF & CS coil conductor type:
  !!
  !! - =0 superconducting PF coils
  !! - =1 resistive PF coils
  !
  real(dp) :: itr_sum
  !! total sum of I x turns x radius for all PF coils and CS (Am)

  integer :: i_cs_superconductor
  !! switch for superconductor material in central solenoid:
  !!
  !! - =1 ITER Nb3Sn critical surface model with standard
  !!   ITER parameters
  !! - =2 Bi-2212 high temperature superconductor (range of
  !!   validity T < 20K, adjusted field b < 104 T, B > 6 T)
  !! - =3 NbTi
  !! - =4 ITER Nb3Sn model with user-specified parameters
  !! - =5 WST Nb3Sn parameterisation
  !! - =6 REBCO HTS tape in CroCo strand
  !! - =7 Durham Ginzburg-Landau critical surface model for Nb-Ti
  !! - =8 Durham Ginzburg-Landau critical surface model for REBCO
  !! - =9 Hazelton experimental data + Zhai conceptual model for REBCO

  integer :: i_pf_superconductor
  !! switch for superconductor material in PF coils:
  !!
  !! - =1 ITER Nb3Sn critical surface model with standard
  !!   ITER parameters
  !! - =2 Bi-2212 high temperature superconductor (range of
  !!   validity T < 20K, adjusted field b < 104 T, B > 6 T)
  !! - =3 NbTi
  !! - =4 ITER Nb3Sn model with user-specified parameters
  !! - =5 WST Nb3Sn parameterisation
  !! - =6 REBCO HTS tape in CroCo strand
  !! - =7 Durham Ginzburg-Landau critical surface model for Nb-Ti
  !! - =8 Durham Ginzburg-Landau critical surface model for REBCO
  !! - =9 Hazelton experimental data + Zhai conceptual model for REBCO

  real(dp) :: j_crit_str_cs
  !! superconductor strand critical current density under operating
  !! conditions in central solenoid (A/m2). Necessary for the cost calculation in $/kA m

  real(dp) :: j_crit_str_pf
  !! superconductor strand critical current density under operating
  !! conditions in PF coils (A/m2). Necessary for the cost calculation in $/kA m

  integer :: i_pf_current
  !! Switch for controlling the current of the PF coils:
  !!
  !! - =0 Input via the variables c_pf_cs_coil_pulse_start_ma, c_pf_cs_coil_flat_top_ma, c_pf_cs_coil_pulse_end_ma
  !! - =1 SVD targets zero field across midplane (flux swing
  !!   coils) and the correct vertical field at the plasma
  !!   center (equilibrium coils)

  integer :: i_sup_pf_shape
  !! Switch for the placement of Location 3 (outboard) PF coils
  !! when the TF coils are superconducting (i_tf_sup = 1)
  !!
  !! - =0 (Default) Outboard PF coils follow TF shape
  !!   in an ellipsoidal winding surface
  !! - =1 Outboard PF coils all have same radius, cylindrical
  !!   winding surface

  real(dp) :: j_cs_conductor_critical_pulse_start
  !! central solenoid superconductor critical current density (A/m2) at beginning-of-pulse

  real(dp) :: j_cs_conductor_critical_flat_top_end
  !! central solenoid superconductor critical current density (A/m2) at end-of-flattop

  real(dp) :: jcableoh_bop
  !! central solenoid cable critical current density (A/m2) at beginning-of-pulse

  real(dp) :: jcableoh_eof
  !! central solenoid cable critical current density (A/m2) at end-of-flattop

  integer :: n_pf_cs_plasma_circuits
  !! number of PF circuits (including central solenoid and plasma)

  integer, dimension(n_pf_groups_max+2) :: n_pf_coils_in_group
  !! number of PF coils in group j

  integer :: nfxfh
  !! number of filaments the top and bottom of the central solenoid should be broken
  !! into during scaling (5 - 10 is good)

  integer :: n_pf_coil_groups
  !! number of groups of PF coils. Symmetric coil pairs should all be in the same group

  integer :: n_cs_pf_coils
  !! number of PF coils (excluding the central solenoid) + 1

  real(dp) :: f_z_cs_tf_internal
  !! Central solenoid height / TF coil internal height

  real(dp) :: f_a_cs_steel
  !! central solenoid steel fraction (`iteration variable 122`)

  real(dp) :: pf_current_safety_factor
  !! Ratio of permissible PF coil conductor current density to critical conductor
  !! current density based on short-sample DC measurements

  real(dp), dimension(ngc2) :: pfcaseth
  !! steel case thickness for PF coil i (m)

  real(dp) :: rho_pf_coil
  !! PF coil resistivity (if i_pf_conductor=1) (Ohm-m)

  real(dp) :: rhopfbus
  !! Resistivity of CS and PF coil bus bars (irrespective of
  !! whether the coils themselves are superconducting or resistive) (Ohm-m)

  real(dp) :: m_pf_coil_max
  !! mass of heaviest PF coil (tonnes)

  real(dp) :: r_pf_coil_outer_max
  !! radius of largest PF coil (m)

  real(dp) :: p_pf_electric_supplies_mw
  !! Total mean wall plug power dissipated in PFC and CS power supplies (MW) (issue #713)

  real(dp) :: p_cs_resistive_flat_top
  !! central solenoid resistive power during flattop (W)

  real(dp) :: p_pf_coil_resistive_total_flat_top
  !! total PF coil resistive losses during flattop (W)

  real(dp), dimension(ngc2) :: r_pf_coil_inner
  !! inner radius of coil i (m)

  real(dp), dimension(ngc2) :: r_pf_coil_outer
  !! outer radius of coil i (m)

  real(dp), dimension(ngc2) :: c_pf_cs_coils_peak_ma
  !! peak current in coil i (MA-turns)

  real(dp), dimension(ngc2) :: j_pf_coil_wp_peak
  !! average winding pack current density of PF coil i (A/m2) at time of peak
  !! current in that coil (calculated for `i_pf_location=1` coils)

  real(dp) :: j_cs_critical_flat_top_end
  !! allowable central solenoid current density at end of flat-top (A/m2)

  real(dp) :: j_cs_critical_pulse_start
  !! allowable central solenoid current density at beginning of pulse (A/m2)

  real(dp), dimension(ngc2) :: j_pf_wp_critical
  !! allowable winding pack current density of PF coil i (A/m2)

  real(dp) :: r_cs_middle
  !! radius to the centre of the central solenoid (m)

  real(dp) :: routr
  !! radial distance (m) from outboard TF coil leg to centre of `i_pf_location=3` PF coils

  real(dp), dimension(ngc2) :: r_pf_coil_middle
  !! radius of PF coil i (m)

  real(dp) :: rpf1
  !! offset (m) of radial position of `i_pf_location=1` PF coils from being directly above
  !! the central solenoid

  real(dp) :: rpf2
  !! offset (m) of radial position of `i_pf_location=2` PF coils from being at
  !! rmajor (offset = rpf2*triang*rminor)

  real(dp), dimension(n_pf_groups_max) :: rref
  !! PF coil radial positioning adjuster:
  !!
  !! - for groups j with i_pf_location(j) = 1; rref(j) is ignored
  !! - for groups j with i_pf_location(j) = 2; rref(j) is ignored
  !! - for groups j with i_pf_location(j) = 3; rref(j) is ignored
  !! - for groups j with i_pf_location(j) = 4; rref(j) is radius of
  !!   the coil in units of minor radii from the major radius
  !!   (r = rmajor + rref*rminor)

  real(dp) :: s_shear_cs_peak
  !! Maximum shear stress (Tresca criterion) coils/central solenoid [MPa]

  real(dp) :: sigpfcalw
  !! maximum permissible tensile stress (MPa) in steel coil cases for superconducting
  !! PF coils (`i_pf_conductor=0`)

  real(dp) :: sigpfcf
  !! fraction of JxB hoop force supported by steel case for superconducting PF coils (`i_pf_conductor=0`)

  real(dp), dimension(ngc2,ngc2) :: ind_pf_cs_plasma_mutual
  !! mutual inductance matrix (H)

  real(dp) :: temp_cs_margin
  !! Central solenoid temperature margin (K)

  real(dp), dimension(ngc2) :: n_pf_coil_turns
  !! number of turns in PF coil i

  real(dp), dimension(ngc2) :: f_a_pf_coil_void
  !! winding pack void fraction of PF coil i for coolant

  real(dp) :: f_a_cs_void
  !! void fraction of central solenoid conductor for coolant

  real(dp) :: vs_cs_pf_total_burn
  !! total flux swing available for burn (Wb)

  real(dp) :: vs_pf_coils_total_burn
  !! flux swing from PF coils for burn (Wb)

  real(dp) :: vs_pf_coils_total_ramp
  !! flux swing from PF coils for startup (Wb)

  real(dp) :: vs_pf_coils_total_pulse
  !! total flux swing from PF coils (Wb)

  real(dp) :: vs_cs_total_pulse
  !! total flux swing from the central solenoid (Wb)

  real(dp) :: vs_cs_burn
  !! central solenoid flux swing for burn (Wb)

  real(dp) :: vs_cs_ramp
  !! central solenoid flux swing for startup (Wb)

  real(dp) :: vs_cs_pf_total_ramp
  !! total flux swing for startup (`constraint eqn 51` to enforce vs_cs_pf_total_ramp=vs_plasma_res_ramp+vs_plasma_ind_ramp) (Wb)

  real(dp) :: vs_cs_pf_total_pulse
  !! total flux swing for pulse (Wb)

  real(dp), dimension(ngc2,6) :: waves
  !! used in current waveform of PF coils/central solenoid

  real(dp) :: m_pf_coil_conductor_total
  !! total mass of the PF coil conductor (kg)

  real(dp) :: m_pf_coil_structure_total
  !! total mass of the PF coil structure (kg)

  real(dp), dimension(ngc2) :: m_pf_coil_conductor
  !! conductor mass for PF coil i (kg)

  real(dp), dimension(ngc2) :: m_pf_coil_structure
  !! structure mass for PF coil i (kg)

  real(dp), dimension(ngc2) :: z_pf_coil_upper
  !! upper point of PF coil i (m)

  real(dp), dimension(ngc2) :: z_pf_coil_lower
  !! lower point of PF coil i (m)

  real(dp), dimension(ngc2) :: z_pf_coil_middle
  !! z (height) location of PF coil i (m)

  real(dp), dimension(n_pf_groups_max) :: zref
  !! PF coil vertical positioning adjuster:
  !!
  !! - for groups j with i_pf_location(j) = 1; zref(j) is ignored
  !! - for groups j with i_pf_location(j) = 2 AND itart=1 (only);
  !!   zref(j) is distance of centre of PF coil from inside
  !!   edge of TF coil (remember that PF coils for STs lie
  !!   within the TF coil)
  !! - for groups j with i_pf_location(j) = 3; zref(j) = ratio of
  !!   height of coil group j to plasma minor radius</UL>
  !! - for groups j with i_pf_location(j) = 4; zref(j) = ratio of
  !!   height of coil group j to plasma minor radius</UL>

  real(dp) :: b_cs_limit_max
  !! Central solenoid max field limit [T]

  real(dp) :: fb_cs_limit_max
  !! F-value for CS mmax field (`cons. 79`, `itvar 149`)

  real(dp) :: ld_ratio_cst
  !! Ratio of CS coil turn conduit length to depth

  real(dp) :: l_cond_cst
  !! Length of CS of CS coil turn conduit

  real(dp) :: dz_cs_turn
  !! Depth/width of CS of CS coil turn conduit

  real(dp) :: r_out_cst
  !! Length of CS of CS coil turn conduit length

  real(dp) :: r_in_cst
  !! Length of CS of CS coil turn conduit length
end module pfcoil_variables
