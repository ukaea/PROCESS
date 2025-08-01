module tfcoil_variables
  !! author: J. Morris, M. Kovari, S. Kahn (UKAEA)
  !!
  !! Module containing global variables relating to the toroidal field coil systems
  !!
  !!### References
  !!
  !! - ITER Magnets design description document DDD11-2 v2 2 (2009)

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: a_tf_coil_inboard_case
  !! external case area per coil (inboard leg) (m2)

  real(dp) :: a_tf_coil_outboard_case
  !! external case area per coil (outboard leg) (m2)

  real(dp) :: a_tf_turn_steel
  !! area of the cable conduit (m2)

  real(dp) :: a_tf_wp_conductor
  !! Winding pack conductor area [m2]
  !! Does not include the area of voids and central helium channel

  real(dp) :: a_tf_turn_cable_space_no_void
  !! Cable space area (per turn)  [m2]
  !! Includes the area of voids and central helium channel

  real(dp) :: a_tf_turn_insulation
  !! single turn insulation area (m2)

  real(dp) :: a_tf_coil_wp_turn_insulation
  !! winding pack turn insulation area per coil (m2)

  real(dp) :: sig_tf_case_max
  !! Allowable maximum shear stress (Tresca criterion) in TF coil case (Pa)

  real(dp) :: sig_tf_wp_max
  !! Allowable maximum shear stress (Tresca criterion) in TF coil conduit (Pa)

  ! TODO remove below IF not needed
  ! real(dp) :: alstrtf
  !! Allowable Tresca stress in TF coil structural material (Pa)

  real(dp) :: a_tf_leg_outboard
  !! outboard TF leg area (m2)

  real(dp) :: a_tf_wp_steel
  !! winding pack structure area (m2)

  real(dp) :: a_tf_wp_extra_void
  !! winding pack void (He coolant) area (m2)

  real(dp) :: a_tf_wp_coolant_channels
  !! winding pack He coil area (m2)

  real(dp) :: bcritsc
  !! upper critical field (T) for Nb3Sn superconductor at zero temperature and
  !! strain (`i_tf_sc_mat=4, =bc20m`)

  real(dp) :: b_tf_inboard_peak
  !! mean peak field at TF coil (T)

  real(dp) :: bmaxtfrp
  !! peak field at TF conductor with ripple (T)

  real(dp) :: casestr
  !! case strain

  real(dp) :: dr_tf_plasma_case
  !! inboard TF coil case plasma side thickness (m) (calculated for stellarators)

  real(dp) :: f_dr_tf_plasma_case
  !! inboard TF coil case plasma side thickness as a fraction of dr_tf_inboard

  logical :: i_f_dr_tf_plasma_case
  !! logical switch to make dr_tf_plasma_case a fraction of TF coil thickness (`f_dr_tf_plasma_case`)

  real(dp) :: dx_tf_side_case_min
  !! inboard TF coil sidewall case thickness (m) (calculated for stellarators)

  real(dp) :: casths_fraction
  !! inboard TF coil sidewall case thickness as a fraction of dx_tf_inboard_out_toroidal

  logical :: tfc_sidewall_is_fraction
  !! logical switch to make dx_tf_side_case_min a fraction of TF coil thickness (`casths_fraction`)

  real(dp) :: t_conductor
  !! Conductor (cable + steel conduit) area averaged dimension [m]

  real(dp) :: t_turn_tf
  !! TF coil turn edge length including turn insulation [m]
  !!   If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
  !!   equivelent size is use to calculated this quantity
  !!   If the t_turn_tf is non zero, c_tf_turn is calculated

  logical :: t_turn_tf_is_input
  !! Boolean switch to activated when the user set the TF coil turn dimensions
  !! Not an input

  real(dp) :: f_t_turn_tf
  !! f-value for TF turn edge length constraint
  !!  If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
  !!  equivelent size is use for this constraint
  !!  iteration variable ixc = 175
  !!  constraint equation icc = 86

  real(dp) :: t_turn_tf_max
  !! TF turn edge length including turn insulation upper limit [m]
  !! If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
  !! equivelent size is use for this constraint
  !! constraint equation icc = 86

  real(dp) :: t_cable_tf
  !! TF coil superconducting cable squared/rounded dimensions [m]
  !!   If the turn is not a square (i_tf_turns_integer = 1) a squared cable of
  !!   equivelent size is use to calculated this quantity
  !!   If the t_cable_tf is non zero, c_tf_turn is calculated

  logical :: t_cable_tf_is_input
  !! Boolean switch to activated when the user set the TF coil cable dimensions
  !! Not an input

  real(dp) :: acs
  !! Area of space inside conductor (m2)

  real(dp) :: cdtfleg
  !! TF outboard leg current density (A/m2) (resistive coils only)

  real(dp) :: cforce
  !! centering force on inboard leg (per coil) (N/m)

  real(dp) :: cplen
  !! length of TF coil inboard leg ('centrepost') (`i_tf_sup = 1`)

  real(dp) :: c_tf_turn
  !! TF coil current per turn (A). (calculated for stellarators) (calculated for
  !! integer-turn TF coils `i_tf_turns_integer=1`) (`iteration variable 60`)

  real(dp) :: c_tf_turn_max
  !! Max TF coil current per turn [A]. (for stellarators and `i_tf_turns_integer=1`)
  !! (`constraint equation 77`)

  real(dp) :: dcase
  !! density of coil case (kg/m3)

  real(dp), dimension(9) :: dcond
  !! density of superconductor type given by i_tf_sc_mat/i_cs_superconductor/i_pf_superconductor (kg/m3)

  real(dp) :: dcondins
  !! density of conduit + ground-wall insulation (kg/m3)

  real(dp) :: dia_tf_turn_coolant_channel
  !! diameter of central helium channel in TF winding (m)

  real(dp) :: e_tf_magnetic_stored_total_gj
  !! total stored energy in the toroidal field (GJ)

  real(dp) :: b_crit_upper_nbti
  !! upper critical field of GL_nbti
  real(dp) :: t_crit_nbti
  !! critical temperature of GL_nbti
  real(dp) :: max_force_density
  !! Maximal (WP averaged) force density in TF coils at 1 point. (MN/m3)
  real(dp) :: fcutfsu
  !! copper fraction of cable conductor (TF coils)
  !! (iteration variable 59)
  real(dp) :: fhts
  !! technology adjustment factor for critical current density fit for isumat..=2
  !! Bi-2212 superconductor, to describe the level of technology assumed (i.e. to
  !! account for stress, fatigue, radiation, AC losses, joints or manufacturing
  !! variations; 1.0 would be very optimistic)

  real(dp) :: insstrain
  !! Radial strain in insulator

  integer :: i_tf_stress_model
  !! Switch for the TF coil stress model
  !!   0 : Generalized plane strain formulation, Issues #977 and #991, O(n^3)
  !!   1 : Old plane stress model (only for SC)
  !!   2 : Axisymmetric extended plane strain, Issues #1414 and #998, O(n)

  integer :: i_tf_tresca
  !! Switch for TF coil conduit Tresca stress criterion:
  !!   0 : Tresca (no adjustment);
  !!   1 : Tresca with CEA adjustment factors (radial+2%, vertical+60%) </UL>

  integer :: i_tf_wp_geom
  !! Switch for TF WP geometry selection
  !!   0 : Rectangular geometry
  !!   1 : Double rectangular geometry
  !!   2 : Trapezoidal geometry (constant lateral casing thickness)
  !! Default setting for backward compatibility
  !!   if i_tf_turns_integer = 0 : Double rectangular
  !!   if i_tf_turns_integer = 1 : Rectangular

  integer :: i_tf_case_geom
  !! Switch for TF case geometry selection
  !!   0 : Circular front case (ITER design)
  !!   1 : Straight front case

  integer :: i_tf_turns_integer
  !! Switch for TF coil integer/non-integer turns:
  !!   0 : non-integer turns
  !!   1 : integer turns

  integer :: i_tf_sc_mat
  !! Switch for superconductor material in TF coils:
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

  integer :: i_tf_sup
  !! Switch for TF coil conductor model:
  !!
  !! - =0 copper
  !! - =1 superconductor
  !! - =2 Cryogenic aluminium

  integer :: i_tf_shape
  !! Switch for TF coil toroidal shape:
  !!
  !! - =0  Default value : Picture frame coil for TART / PROCESS D-shape for non itart
  !! - =1  PROCESS D-shape : parametrise with 2 arcs
  !! - =2  Picture frame coils

  integer :: i_tf_cond_eyoung_axial
  !! Switch for the behavior of the TF coil conductor elastic axial properties
  !!
  !! - =0  Young's modulus is set to zero, and the conductor is not considered
  !!       in the stress calculation. This corresponds to the case that the
  !!       conductor is much less stiff than the conduit, or the case that the
  !!       conductor is prevented (isolated) from taking axial loads.
  !! - =1  Elastic properties are set by user input, using the variable
  !!       `eyoung_cond_axial`
  !! - =2  Elastic properties are set to reasonable defaults taking into
  !!       account the superconducting material `i_tf_sc_mat`

  integer :: i_tf_cond_eyoung_trans
  !! Switch for the behavior of the elastic properties of the TF coil
  !! conductorin the transverse direction. Only active if
  !! `i_tf_cond_eyoung_axial == 2`
  !!
  !! - =0  Cable not potted in solder. Transverse Young's modulus set to zero.
  !! - =1  Cable potted in solder. If `i_tf_cond_eyoung_axial == 2`, the
  !!       transverse Young's modulus of the conductor is equal to the axial,
  !!       which is set to a sensible material-dependent default.

  integer :: n_pancake
  !! Number of pancakes in TF coil. Only used if `i_tf_turns_integer=1`

  integer :: n_layer
  !! Number of layers in TF coil. Only used if `i_tf_turns_integer=1`

  integer :: n_rad_per_layer
  !! Size of the arrays per layers storing the radial dependent stress
  !! quantities (stresses, strain displacement etc..)

  integer :: i_tf_bucking
  !! Switch for TF inboard suport structure design:
  !!
  !! Default setting for backward compatibility
  !!     - if copper resistive TF (i_tf_sup = 0) : Free standing TF without bucking structure
  !!     - if Superconducting TF  (i_tf_sup = 1) : Free standing TF with a steel casing
  !!     - if aluminium  TF       (i_tf_sup = 2) : Free standing TF with a bucking structure
  !!     Rem : the case is a bucking structure
  !! - =0 : Free standing TF without case/bucking cyliner (only a conductor layer)
  !! - =1 : Free standing TF with a case/bucking cylinder made of
  !!     - if copper resistive     TF (i_tf_sup = 0) : used defined bucking cylinder
  !!     - if Superconducting      TF (i_tf_sup = 1) : Steel casing
  !!     - if aluminium resisitive TF (i_tf_sup = 2) : used defined bucking cylinder
  !! - =2 : The TF is in contact with the CS : "bucked and wedged design"
  !!       Fast version : thin TF-CS interface neglected in the stress calculations (3 layers)
  !!                      The CS is frictionally decoupled from the TF, does not carry axial tension
  !! - =3 : The TF is in contact with the CS : "bucked and wedged design"
  !!       Full version : thin TF-CS Kapton interface introduced in the stress calculations (4 layers)
  !!                      The CS and kaptop are frictionally decoupled from the TF, do not carry
  !!                      axial tension

  integer :: n_tf_graded_layers
  !! Number of layers of different stress properties in the WP. If `n_tf_graded_layers > 1`,
  !! a graded coil is condidered

  integer :: n_tf_stress_layers
  !! Number of layers considered for the inboard TF stress calculations
  !! set in initial.f90 from i_tf_bucking and n_tf_graded_layers

  integer :: n_tf_wp_layers
  !! Maximum number of layers that can be considered in the TF coil composited/smeared
  !! stress analysis. This is the layers of one turn, not the entire WP.
  !! Default: 5. void, conductor, copper, conduit, insulation.

  real(dp) :: j_tf_bus
  !! bussing current density (A/m2)

  real(dp) :: j_crit_str_tf
  !! j_crit_str : superconductor strand critical current density under operating
  !! conditions (A/m2). Necessary for the cost calculation in $/kAm

  real(dp), dimension(9) :: j_crit_str_0
  !! j_crit_str_pf_0 : superconductor strand critical current density at 6 T and 4.2 K (A/m2)
  !! Necessary for the cost calculation in $/kAm

  real(dp) :: j_tf_wp_critical
  !! critical current density for winding pack (A/m2)

  real(dp) :: jwdgpro
  !! allowable TF coil winding pack current density, for dump temperature rise protection (A/m2)

  real(dp) :: j_tf_wp
  !! winding pack engineering current density (A/m2)

  real(dp) :: oacdcp
  !! Overall current density in TF coil inboard legs midplane (A/m2)
  !! Rem SK : Not used in tfcoil to set the current any more. Should not be used as
  !! iteration variable 12 any more. It is now calculated.

  real(dp) :: eyoung_ins
  !! Insulator Young's modulus [Pa]. Default value (1.0D8) setup the following values
  !!  - SC TF, eyoung_ins = 20 Gpa (default value from DDD11-2 v2 2 (2009))
  !!  - Al TF, eyoung_ins = 2.5 GPa (Kapton polymer)

  real(dp) :: eyoung_steel
  !! Steel case Young's modulus (Pa) (default value from DDD11-2 v2 2 (2009))

  real(dp) :: eyoung_cond_axial
  !! SC TF coil conductor Young's modulus in the parallel (along the wire/tape)
  !! direction [Pa]
  !! Set by user input only if `i_tf_cond_eyoung_axial == 1`; otherwise
  !! set by the behavior of that switch.

  real(dp) :: eyoung_cond_trans
  !! SC TF coil conductor Young's modulus in the transverse direction [Pa]
  !! Set by user input only if `i_tf_cond_eyoung_axial == 1`; otherwise
  !! set by the behavior of that switch.

  real(dp) :: eyoung_res_tf_buck
  !! Resistive TF magnets bucking cylinder young modulus (Pa)

  real(dp) :: eyoung_copper
  !! Copper young modulus. Default value taken from wikipedia

  real(dp) :: eyoung_al
  !! Aluminium young modulus.  Default value taken from wikipedia

  real(dp) :: poisson_steel
  !! Steel Poisson's ratio, Source : https://www.engineeringtoolbox.com/metals-poissons-ratio-d_1268.html

  real(dp):: poisson_copper
  !! Copper Poisson's ratio. Source : https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html

  real(dp):: poisson_al
  !! Aluminium Poisson's ratio.
  !! Source : https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html

  real(dp):: poisson_ins
  !! Insulation Poisson's ratio. Default: Kapton.
  !! Source : DuPont™ Kapton® HN datasheet.

  real(dp) :: poisson_cond_axial
  !! SC TF coil conductor Poisson's ratio in the parallel-transverse direction

  real(dp) :: poisson_cond_trans
  !! SC TF coil conductor Poisson's ratio in the transverse-transverse direction

  real(dp) :: r_b_tf_inboard_peak
  !! Radius of maximum TF B-field (m)

  real(dp) :: res_tf_leg
  !! TF coil leg resistance (ohm)

  real(dp) :: toroidalgap
  !! Minimal distance between two toroidal coils. (m)

  real(dp) :: ftoroidalgap
  !! F-value for minimum dx_tf_inboard_out_toroidal (`constraint equation 82`)

  real(dp) :: ripmax
  !! aximum allowable toroidal field ripple amplitude at plasma edge (%)

  real(dp) :: ripple
  !! peak/average toroidal field ripple at plasma edge (%)

  real(dp) :: c_tf_total
  !! total (summed) current in TF coils (A)

  integer, parameter :: n_radial_array = 50
  !! Size of the radial distribution arrays per layers
  !! used for stress, strain and displacement distibution

  real(dp), dimension(2*n_radial_array) :: radial_array
  !! Array refining the radii of the stress calculations arrays

  real(dp), dimension(2*n_radial_array) :: sig_tf_r
  !! TF Inboard leg radial stress in steel r distribution at mid-plane [Pa]

  real(dp), dimension(2*n_radial_array) :: sig_tf_t
  !! TF Inboard leg tangential stress in steel r distribution at mid-plane [Pa]

  real(dp), dimension(2*n_radial_array) :: deflect
  !! TF coil radial deflection (displacement) radial distribution [m]

  real(dp) :: sig_tf_z
  !! TF Inboard leg vertical tensile stress in steel at mid-plane [Pa]

  real(dp), dimension(2*n_radial_array) :: sig_tf_vmises
  !! TF Inboard leg Von-Mises stress in steel r distribution at mid-plane [Pa]

  real(dp), dimension(2*n_radial_array) :: s_shear_tf
  !! TF Inboard leg maximum shear stress (Tresca criterion) in steel r distribution at mid-plane [Pa]

  real(dp) :: sig_tf_cs_bucked

  ! TODO is this needed?
  ! real(dp) :: strtf0
  !! Maximum shear stress (Tresca criterion) in CS structures at CS flux swing [Pa]:
  !!
  !!  - If superconducting CS (i_pf_conductor = 0): turn steel conduits stress
  !!  - If resistive       CS (i_pf_conductor = 1): copper conductor stress
  !!
  !! Quantity only computed for bucked and wedged design (`i_tf_bucking >= 2`)
  !! Def : CS Flux swing, instant when the current changes sign in CS (null current)

  real(dp) :: sig_tf_case
  !! Maximum shear stress (Tresca criterion) in TF casing steel structures (Pa)

  real(dp) :: sig_tf_wp

  ! TODO is this needed?
  ! real(dp) :: strtf1
  ! !! Maximum TRESCA stress in TF casing steel structures (Pa)

  ! real(dp) :: strtf2
  ! !! Maximum TRESCA stress in TF WP conduit steel structures (Pa)
  ! !! This is the TF stress condition used in the case of stellarators

  real(dp) :: str_cs_con_res
  !! Residual manufacturing strain in CS superconductor material

  real(dp) :: str_pf_con_res
  !! Residual manufacturing strain in PF superconductor material

  real(dp) :: str_tf_con_res
  !! Residual manufacturing strain in TF superconductor material
  !! If `i_str_wp == 0`, used to compute the critical surface.
  !! Otherwise, the self-consistent winding pack `str_wp` is used.

  real(dp) :: str_wp
  !! Axial (vertical) strain in the TF coil winding pack found by
  !! self-consistent stress/strain calculation.
  !! if `i_str_wp == 1`, used to compute the critical surface.
  !! Otherwise, the input value `str_tf_con_res` is used.
  !! Constrain the absolute value using `constraint equation 88`
  !! You can't have constraint 88 and i_str_wp = 0 at the same time

  real(dp) :: str_wp_max
  !! Maximum allowed absolute value of the strain in the TF coil
  !! (`Constraint equation 88`)

  integer :: i_str_wp
  !! Switch for the behavior of the TF strain used to compute
  !! the strain-dependent critical surface:
  !!
  !! - =0  str_tf_con_res is used
  !! - =1  str_wp is used

  character(len=12) :: quench_model
  !! switch for TF coil quench model (Only applies to REBCO magnet at present, issue #522):
  !!
  !! - ='exponential' exponential quench with constant discharge resistor
  !! - ='linear' quench with constant voltage

  real(dp) :: time1
  !! Time at which TF quench is detected (s)

  real(dp) :: tcritsc
  !! critical temperature (K) for superconductor at zero field and strain (`i_tf_sc_mat=4, =tc0m`)

  real(dp) :: tdmptf
  !! fast discharge time for TF coil in event of quench (s) (`iteration variable 56`)
  !!
  !! For REBCO model, meaning depends on quench_model:
  !!
  !! - exponential quench : e-folding time (s)`
  !! - linear quench : discharge time (s)

  real(dp) :: a_tf_inboard_total
  !! Area of inboard midplane TF legs (m2)

  real(dp) :: len_tf_bus
  !! TF coil bus length (m)

  real(dp) :: m_tf_bus
  !! TF coil bus mass (kg)

  real(dp) :: tfckw
  !! available DC power for charging the TF coils (kW)

  real(dp) :: tfcmw
  !! Peak power per TF power supply (MW)

  real(dp) :: p_cp_resistive_mw
  !! Peak resistive TF coil inboard leg power (MW)

  real(dp) :: p_tf_joints_resistive_mw
  !! TF joints resistive power losses (MW)

  real(dp) :: tfcryoarea
  !! surface area of toroidal shells covering TF coils (m2)

  real(dp) :: tficrn
  !! TF coil half-width - inner dr_bore (m)

  real(dp) :: ind_tf_coil
  !! TF coil inductance (H)

  real(dp) :: dx_tf_wp_insertion_gap
  !! TF coil WP insertion gap (m)

  real(dp) :: p_tf_leg_resistive_mw
  !! TF coil outboard leg resistive power (MW)

  real(dp) :: rho_cp
  !! TF coil inboard leg resistivity [Ohm-m]. If `itart=0`, this variable is the
  !! average resistivity over the whole magnet

  real(dp) :: rho_tf_leg
  !! Resistivity of a TF coil leg (Ohm-m)

  real(dp) :: rho_tf_bus
  !! Resistivity of a TF coil bus (Ohm-m). Default values is for that of GLIDCOP AL-15 (C15715) at 293K

  real(dp) :: frhocp
  !! Centrepost resistivity enhancement factor. For `itart=0`, this factor
  !! is used for the whole magnet

  real(dp) :: frholeg
  !! Ouboard legs resistivity enhancement factor. Only used for `itart=1`.

  integer :: i_cp_joints
  !! Switch for CP demoutable joints type
  !!  -= 0 : Clampled joints
  !!  -= 1 : Sliding joints
  !! Default value (-1) choses :
  !!   Sliding joints for resistive magnets (i_tf_sup = 0, 2)
  !!   Clampled joints for superconducting magents (i_tf_sup = 1)

  real(dp) :: rho_tf_joints
  !! TF joints surfacic resistivity [ohm.m]. Feldmetal joints assumed.

  integer :: n_tf_joints_contact
  !! Number of contact per turn

  integer :: n_tf_joints
  !! Number of joints
  !! Ex: n_tf_joints = 2 for top and bottom CP joints

  real(dp) :: th_joint_contact
  !! TF sliding joints contact pad width [m]

  real(dp) :: p_tf_joints_resistive
  !! Calculated TF joints resistive power losses [W]

  real(dp) :: len_tf_coil
  !! TF coil circumference (m)

  real(dp) :: eff_tf_cryo
  !! TF cryoplant efficiency (compared to pefect Carnot cycle).
  !! Using -1 set the default value depending on magnet technology:
  !!
  !!  - i_tf_sup = 1 : SC magnet, eff_tf_cryo = 0.13 (ITER design)
  !!  - i_tf_sup = 2 : Cryo-aluminium, eff_tf_cryo = 0.4

  real(dp) :: n_tf_coils
  !! Number of TF coils (default = 50 for stellarators). Number of TF coils outer legs for ST

  real(dp) :: tfocrn
  !! TF coil half-width - outer dr_bore (m)

  real(dp) :: tfsai
  !! area of the inboard TF coil legs (m2)

  real(dp) :: tfsao
  !! area of the outboard TF coil legs (m2)

  real(dp) :: tftmp
  !! peak helium coolant temperature in TF coils and PF coils (K)

  real(dp) :: dx_tf_inboard_out_toroidal
  !! TF coil toroidal thickness (m)

  real(dp) :: dx_tf_turn_insulation
  !! conduit insulation thickness (m)

  real(dp) :: layer_ins
  !! Additional insulation thickness between layers (m)

  real(dp) :: dr_tf_nose_case
  !! inboard TF coil case outer (non-plasma side) thickness (m) (`iteration variable 57`)
  !! (calculated for stellarators)

  real(dp) :: dr_tf_wp_with_insulation
  !! radial thickness of winding pack (m) (`iteration variable 140`) (issue #514)

  real(dp) :: dx_tf_turn_steel
  !! TF coil conduit case thickness (m) (`iteration variable 58`)

  real(dp) :: dx_tf_wp_insulation
  !! Thickness of the ground insulation layer surrounding (m)
  !!   - Superconductor TF (`i_tf_sup == 1`) : The TF coil Winding packs
  !!   - Resistive magnets (`i_tf_sup /= 1`) : The TF coil wedges
  !! Rem : Thickness calculated for stellarators.

  real(dp) :: tmargmin_tf
  !! minimum allowable temperature margin : TF coils (K)

  real(dp) :: tmargmin_cs
  !! minimum allowable temperature margin : CS (K)

  real(dp) :: tmargmin
  !! minimum allowable temperature margin : TFC AND CS (K)

  real(dp) :: temp_margin
  !! temperature margin (K)

  real(dp) :: tmargtf
  !! TF coil temperature margin (K)

  real(dp) :: tmaxpro
  !! maximum temp rise during a quench for protection (K)

  real(dp) :: tmax_croco
  !! CroCo strand: maximum permitted temp during a quench (K)

  real(dp) :: croco_quench_temperature
  !! CroCo strand: Actual temp reached during a quench (K)

  real(dp) :: temp_tf_cryo
  !! coil temperature for cryogenic plant power calculation (K)

  real(dp) :: n_tf_coil_turns
  !! number of turns per TF coil

  real(dp) :: vdalw
  !! max voltage across TF coil during quench (kV) (`iteration variable 52`)

  real(dp) :: vforce
  !! vertical tension on inboard leg/coil (N)

  real(dp) :: f_vforce_inboard
  !! Fraction of the total vertical force taken by the TF inboard leg tension
  !! Not used for resistive `itart=1` (sliding joints)

  real(dp) :: vforce_outboard
  !! Vertical tension on outboard leg/coil (N)

  real(dp) :: f_a_tf_turn_cable_space_extra_void
  !! coolant fraction of TFC 'cable' (`i_tf_sup=1`), or of TFC leg (`i_tf_ssup=0`)

  real(dp) :: voltfleg
  !! volume of each TF coil outboard leg (m3)

  real(dp) :: vtfkv
  !! TF coil voltage for resistive coil including bus (kV)

  real(dp) :: vtfskv
  !! voltage across a TF coil during quench (kV)

  real(dp) :: whtcas
  !! mass per coil of external case (kg)

  real(dp) :: whtcon
  !! TF coil conductor mass per coil (kg/coil).
  !! For `itart=1`, coil is return limb plus centrepost/n_tf_coils

  real(dp) :: whtconcu
  !! copper mass in TF coil conductor (kg/coil).
  !! For `itart=1`, coil is return limb plus centrepost/n_tf_coils

  real(dp) :: whtconal
  !! Aluminium mass in TF coil conductor (kg/coil).
  !! For `itart=1`, coil is return limb plus centrepost/n_tf_coils

  real(dp) :: whtconin
  !! conduit insulation mass in TF coil conductor (kg/coil)

  real(dp) :: whtconsc
  !! superconductor mass in TF coil cable (kg/coil)

  real(dp) :: m_tf_turn_steel_conduit
  !! steel conduit mass in TF coil conductor (kg/coil)

  real(dp) :: whtgw
  !! mass of ground-wall insulation layer per coil (kg/coil)

  real(dp) :: m_tf_coils_total
  !! total mass of the TF coils (kg)

  real(dp) :: dx_tf_wp_primary_toroidal
  !! width of first step of winding pack (m)

  real(dp) :: dx_tf_wp_secondary_toroidal
  !! width of second step of winding pack (m)

  ! Superconducting TF coil shape parameters;
  ! the TF inner surface top half is approximated by four circular arcs.
  ! Arc 1 goes through points 1 and 2 on the inner surface. Arc 2
  ! goes through points 2 and 3, etc.
  real(dp), dimension(4) :: dthet
  !! angle of arc i (rad)

  real(dp), dimension(4) :: radctf
  !! radius of arc i (m)

  real(dp), dimension(5) :: r_tf_arc
  !! x location of arc point i on surface (m)

  real(dp), dimension(4) :: xctfc
  !! x location of arc centre i (m)

  real(dp), dimension(5) :: z_tf_arc
  !! y location of arc point i on surface (m)

  real(dp), dimension(4) :: yctfc
  !! y location of arc centre i (m)

  ! New TF shape:  Horizontal and vertical radii of inside edge of TF coil
  ! Arcs are numbered clockwise:
  ! 1=upper inboard, 2=upper outboard, 3=lower ouboard, 4=lower inboard

  real(dp), dimension(4) :: tfa
  !! Horizontal radius of inside edge of TF coil (m)

  real(dp), dimension(4) :: tfb
  !! Vertical radius of inside edge of TF coil (m)
  ! Quantities relating to the spherical tokamak model (itart=1)
  ! (and in some cases, also to resistive TF coils, i_tf_sup=0):

  real(dp) :: drtop
  !! centrepost taper maximum radius adjustment (m)

  real(dp) :: dztop
  !! centrepost taper height adjustment (m)

  real(dp) :: etapump
  !! centrepost coolant pump efficiency

  real(dp) :: fcoolcp
  !! coolant fraction of TF coil inboard legs (`iteration variable 23`)

  real(dp) :: f_a_tf_cool_outboard
  !! coolant fraction of TF coil outboard legs

  real(dp) :: a_cp_cool
  !! Centrepost cooling area toroidal cross-section (constant over the whole CP)

  real(dp) :: ncool
  !! number of centrepost coolant tubes

  real(dp) :: p_cp_coolant_pump_elec
  !! centrepost coolant pump power (W)

  real(dp) :: p_cp_resistive
  !! resistive power in the centrepost (itart=1) [W].
  !! If `itart=0`, this variable is the ressitive power on the whole magnet

  real(dp) :: p_tf_leg_resistive
  !! Summed resistive power in the TF coil legs [W]. Remain 0 if `itart=0`.

  real(dp) :: ptempalw
  !! maximum peak centrepost temperature (K) (`constraint equation 44`)

  real(dp) :: rcool
  !! average radius of coolant channel (m) (`iteration variable 69`)

  real(dp) :: tcoolin
  !! centrepost coolant inlet temperature (K)

  real(dp) :: dtiocool
  !! inlet / outlet TF coil coolant temperature rise (K)

  real(dp) :: temp_cp_average
  !! Average temperature of centrepost called CP (K). Only used for resistive coils
  !! to compute the resisitive heating. Must be an iteration variable for
  !! ST (`itart=1`) (`iteration variable 20`)

  real(dp) :: tcpav2
  !! Computed centrepost average temperature (K) (for consistency)

  real(dp) :: temp_tf_legs_outboard
  !! Average temperature of the TF outboard legs [K]. If `temp_tf_legs_outboard=-1.0`, the ouboard
  !! legs and CP temperatures are the same. Fixed for now, should use a contraints eq like temp_cp_average

  real(dp) :: tcpmax
  !! peak centrepost temperature (K)

  real(dp) :: vcool
  !! inlet centrepost coolant flow speed at midplane (m/s) (`iteration variable 70`)

  real(dp) :: vol_cond_cp
  !! Exact conductor volume in the centrepost (m3)

  real(dp) :: whtcp
  !! mass of TF coil inboard legs (kg)

  real(dp) :: whttflgs
  !! mass of the TF coil legs (kg)

  real(dp) :: cryo_cool_req
  !! Cryo cooling requirement at helium temp 4.5K (kW)

  real(dp) :: theta1_coil
  !! The angle of the outboard arc forming the TF coil current center line [deg]

  real(dp) :: theta1_vv
  !! The angle of the outboard arc forming the Vacuum Vessel current center line [deg]

  real(dp) :: max_vv_stress
  !! The allowable peak maximum shear stress in the vacuum vessel due to quench and fast discharge of the TF coils [Pa]

  real(dp) :: t_tf_quench_detection
  !! TF coil quench detection time (s). Only used for TF coil quench protection.

  real(dp) :: rrr_tf_cu
  !! TF coil copper residual-resistance-ratio (RRR). Only used for quench protection.
end module tfcoil_variables
