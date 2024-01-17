module tfcoil_variables
  !! author: J. Morris, M. Kovari, S. Kahn (UKAEA)
  !!
  !! Module containing global variables relating to the toroidal field coil systems
  !!
  !!### References
  !!
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !! - ITER Magnets design description document DDD11-2 v2 2 (2009)

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: acasetf
  !! external case area per coil (inboard leg) (m2)

  real(dp) :: acasetfo
  !! external case area per coil (outboard leg) (m2)

  real(dp) :: acndttf
  !! area of the cable conduit (m2)

  real(dp) :: acond
  !! Winding pack conductor area [m2]
  !! Does not include the area of voids and central helium channel

  real(dp) :: acstf
  !! Cable space area (per turn)  [m2]
  !! Includes the area of voids and central helium channel

  real(dp) :: insulation_area
  !! single turn insulation area (m2)

  real(dp) :: aiwp
  !! winding pack turn insulation area per coil (m2)

  real(dp) :: sig_tf_case_max
  !! Allowable maximum shear stress (Tresca criterion) in TF coil case (Pa)

  real(dp) :: sig_tf_wp_max
  !! Allowable maximum shear stress (Tresca criterion) in TF coil conduit (Pa)

  ! TODO remove below IF not needed
  ! real(dp) :: alstrtf
  !! Allowable Tresca stress in TF coil structural material (Pa)

  real(dp) :: arealeg
  !! outboard TF leg area (m2)

  real(dp) :: aswp
  !! winding pack structure area (m2)

  real(dp) :: avwp
  !! winding pack void (He coolant) area (m2)

  real(dp) :: awphec
  !! winding pack He coil area (m2)

  real(dp) :: bcritsc
  !! upper critical field (T) for Nb3Sn superconductor at zero temperature and
  !! strain (`i_tf_sc_mat=4, =bc20m`)

  real(dp) :: bmaxtf
  !! mean peak field at TF coil (T)

  real(dp) :: bmaxtfrp
  !! peak field at TF conductor with ripple (T)

  real(dp) :: casestr
  !! case strain

  real(dp) :: casthi
  !! inboard TF coil case plasma side thickness (m) (calculated for stellarators)

  real(dp) :: casthi_fraction
  !! inboard TF coil case plasma side thickness as a fraction of tfcth

  logical :: casthi_is_fraction
  !! logical switch to make casthi a fraction of TF coil thickness (`casthi_fraction`)

  real(dp) :: casths
  !! inboard TF coil sidewall case thickness (m) (calculated for stellarators)

  real(dp) :: casths_fraction
  !! inboard TF coil sidewall case thickness as a fraction of tftort

  logical :: tfc_sidewall_is_fraction
  !! logical switch to make casths a fraction of TF coil thickness (`casths_fraction`)

  real(dp) :: t_conductor
  !! Conductor (cable + steel conduit) area averaged dimension [m]

  real(dp) :: t_turn_tf
  !! TF coil turn edge length including turn insulation [m]
  !!   If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
  !!   equivelent size is use to calculated this quantity
  !!   If the t_turn_tf is non zero, cpttf is calculated

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
  !!   If the t_cable_tf is non zero, cpttf is calculated

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

  real(dp) :: cpttf
  !! TF coil current per turn (A). (calculated for stellarators) (calculated for
  !! integer-turn TF coils `i_tf_turns_integer=1`) (`iteration variable 60`)

  real(dp) :: cpttf_max
  !! Max TF coil current per turn [A]. (for stellarators and `i_tf_turns_integer=1`)
  !! (`constraint equation 77`)

  real(dp) :: dcase
  !! density of coil case (kg/m3)

  real(dp), dimension(9) :: dcond
  !! density of superconductor type given by i_tf_sc_mat/isumatoh/isumatpf (kg/m3)

  real(dp) :: dcondins
  !! density of conduit + ground-wall insulation (kg/m3)

  real(dp) :: dhecoil
  !! diameter of central helium channel in TF winding (m)

  real(dp) :: estotftgj
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

  real(dp) :: jbus
  !! bussing current density (A/m2)

  real(dp) :: jwdgcrt
  !! critical current density for winding pack (A/m2)

  real(dp) :: jwdgpro
  !! allowable TF coil winding pack current density, for dump temperature rise protection (A/m2)

  real(dp) :: jwptf
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

  real(dp) :: rbmax
  !! Radius of maximum TF B-field (m)

  real(dp) :: tflegres
  !! TF coil leg resistance (ohm)

  real(dp) :: toroidalgap
  !! Minimal distance between two toroidal coils. (m)

  real(dp) :: ftoroidalgap
  !! F-value for minimum tftort (`constraint equation 82`)

  real(dp) :: ripmax
  !! aximum allowable toroidal field ripple amplitude at plasma edge (%)

  real(dp) :: ripple
  !! peak/average toroidal field ripple at plasma edge (%)

  real(dp) :: ritfc
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

  real(dp), dimension(2*n_radial_array) :: sig_tf_tresca
  !! TF Inboard leg maximum shear stress (Tresca criterion) in steel r distribution at mid-plane [Pa]

  real(dp) :: sig_tf_cs_bucked

  ! TODO is this needed?
  ! real(dp) :: strtf0
  !! Maximum shear stress (Tresca criterion) in CS structures at CS flux swing [Pa]:
  !!
  !!  - If superconducting CS (ipfres = 0): turn steel conduits stress
  !!  - If resistive       CS (ipfres = 1): copper conductor stress
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

  real(dp) :: tfareain
  !! Area of inboard midplane TF legs (m2)

  real(dp) :: tfbusl
  !! TF coil bus length (m)

  real(dp) :: tfbusmas
  !! TF coil bus mass (kg)

  real(dp) :: tfckw
  !! available DC power for charging the TF coils (kW)

  !#TODO: issue #781
  ! integer :: tfc_model
  ! !! tfc_model /1/ : switch for TF coil magnet stress model:<UL>
  ! !!                 <LI> = 0 simple model (solid copper coil)
  ! !!                 <LI> = 1 CCFE two-layer stress model; superconductor</UL>

  real(dp) :: tfcmw
  !! Peak power per TF power supply (MW)

  real(dp) :: tfcpmw
  !! Peak resistive TF coil inboard leg power (MW)

  real(dp) :: tfjtsmw
  !! TF joints resistive power losses (MW)

  real(dp) :: tfcryoarea
  !! surface area of toroidal shells covering TF coils (m2)

  real(dp) :: tficrn
  !! TF coil half-width - inner bore (m)

  real(dp) :: tfind
  !! TF coil inductance (H)

  real(dp) :: tfinsgap
  !! TF coil WP insertion gap (m)

  real(dp) :: tflegmw
  !! TF coil outboard leg resistive power (MW)

  real(dp) :: rhocp
  !! TF coil inboard leg resistivity [Ohm-m]. If `itart=0`, this variable is the
  !! average resistivity over the whole magnet

  real(dp) :: rhotfleg
  !! Resistivity of a TF coil leg (Ohm-m)

  real(dp) :: rhotfbus
  !! Resistivity of a TF coil bus (Ohm-m). Default value takes the same res as the leg one

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

  real(dp) :: pres_joints
  !! Calculated TF joints resistive power losses [W]

  real(dp) :: tfleng
  !! TF coil circumference (m)

  real(dp) :: eff_tf_cryo
  !! TF cryoplant efficiency (compared to pefect Carnot cycle).
  !! Using -1 set the default value depending on magnet technology:
  !!
  !!  - i_tf_sup = 1 : SC magnet, eff_tf_cryo = 0.13 (ITER design)
  !!  - i_tf_sup = 2 : Cryo-aluminium, eff_tf_cryo = 0.4

  real(dp) :: n_tf
  !! Number of TF coils (default = 50 for stellarators). Number of TF coils outer legs for ST

  real(dp) :: tfocrn
  !! TF coil half-width - outer bore (m)

  real(dp) :: tfsai
  !! area of the inboard TF coil legs (m2)

  real(dp) :: tfsao
  !! area of the outboard TF coil legs (m2)

  real(dp) :: tftmp
  !! peak helium coolant temperature in TF coils and PF coils (K)

  real(dp) :: tftort
  !! TF coil toroidal thickness (m)

  real(dp) :: thicndut
  !! conduit insulation thickness (m)

  real(dp) :: layer_ins
  !! Additional insulation thickness between layers (m)

  real(dp) :: thkcas
  !! inboard TF coil case outer (non-plasma side) thickness (m) (`iteration variable 57`)
  !! (calculated for stellarators)

  real(dp) :: dr_tf_wp
  !! radial thickness of winding pack (m) (`iteration variable 140`) (issue #514)

  real(dp) :: thwcndut
  !! TF coil conduit case thickness (m) (`iteration variable 58`)

  real(dp) :: tinstf
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

  real(dp) :: tmpcry
  !! coil temperature for cryogenic plant power calculation (K)

  real(dp) :: n_tf_turn
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

  real(dp) :: vftf
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
  !! For `itart=1`, coil is return limb plus centrepost/n_tf

  real(dp) :: whtconcu
  !! copper mass in TF coil conductor (kg/coil).
  !! For `itart=1`, coil is return limb plus centrepost/n_tf

  real(dp) :: whtconal
  !! Aluminium mass in TF coil conductor (kg/coil).
  !! For `itart=1`, coil is return limb plus centrepost/n_tf

  real(dp) :: whtconin
  !! conduit insulation mass in TF coil conductor (kg/coil)

  real(dp) :: whtconsc
  !! superconductor mass in TF coil cable (kg/coil)

  real(dp) :: whtconsh
  !! steel conduit mass in TF coil conductor (kg/coil)

  real(dp) :: whtgw
  !! mass of ground-wall insulation layer per coil (kg/coil)

  real(dp) :: whttf
  !! total mass of the TF coils (kg)

  real(dp) :: wwp1
  !! width of first step of winding pack (m)

  real(dp) :: wwp2
  !! width of second step of winding pack (m)

  ! Superconducting TF coil shape parameters;
  ! the TF inner surface top half is approximated by four circular arcs.
  ! Arc 1 goes through points 1 and 2 on the inner surface. Arc 2
  ! goes through points 2 and 3, etc.
  real(dp), dimension(4) :: dthet
  !! angle of arc i (rad)

  real(dp), dimension(4) :: radctf
  !! radius of arc i (m)

  real(dp), dimension(5) :: xarc
  !! x location of arc point i on surface (m)

  real(dp), dimension(4) :: xctfc
  !! x location of arc centre i (m)

  real(dp), dimension(5) :: yarc
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

  real(dp) :: fcoolleg
  !! coolant fraction of TF coil outboard legs

  real(dp) :: a_cp_cool
  !! Centrepost cooling area toroidal cross-section (constant over the whole CP)

  real(dp) :: ncool
  !! number of centrepost coolant tubes

  real(dp) :: ppump
  !! centrepost coolant pump power (W)

  real(dp) :: prescp
  !! resistive power in the centrepost (itart=1) [W].
  !! If `itart=0`, this variable is the ressitive power on the whole magnet

  real(dp) :: presleg
  !! Summed resistive power in the TF coil legs [W]. Remain 0 if `itart=0`.

  real(dp) :: ptempalw
  !! maximum peak centrepost temperature (K) (`constraint equation 44`)

  real(dp) :: rcool
  !! average radius of coolant channel (m) (`iteration variable 69`)

  real(dp) :: tcoolin
  !! centrepost coolant inlet temperature (K)

  real(dp) :: dtiocool
  !! inlet / outlet TF coil coolant temperature rise (K)

  real(dp) :: tcpav
  !! Average temperature of centrepost called CP (K). Only used for resistive coils
  !! to compute the resisitive heating. Must be an iteration variable for
  !! ST (`itart=1`) (`iteration variable 20`)

  real(dp) :: tcpav2
  !! Computed centrepost average temperature (K) (for consistency)

  real(dp) :: tlegav
  !! Average temperature of the TF outboard legs [K]. If `tlegav=-1.0`, the ouboard
  !! legs and CP temperatures are the same. Fixed for now, should use a contraints eq like tcpav

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

  contains

  subroutine init_tfcoil_variables
    !! Initialise module variables
    implicit none

    acasetf = 0.0D0
    acasetfo = 0.0D0
    acndttf = 0.0D0
    acond = 0.0D0
    acstf = 0.0D0
    insulation_area = 0.0D0
    aiwp = 0.0D0
    sig_tf_case_max = 6.0D8
    sig_tf_wp_max = 6.0D8
    arealeg = 0.0D0
    aswp = 0.0D0
    avwp = 0.0D0
    awphec = 0.0D0
    bcritsc = 24.0D0
    bmaxtf = 0.0D0
    bmaxtfrp = 0.0D0
    casestr = 0.0D0
    casthi = 0.0D0
    casthi_fraction = 0.05D0
    casthi_is_fraction = .false.
    casths = 0.0D0
    casths_fraction = 0.06D0
    t_conductor = 0.0D0
    t_cable_tf = 0.0D0
    t_cable_tf_is_input = .false.
    t_turn_tf = 0.0D0
    t_turn_tf_is_input = .false.
    f_t_turn_tf = 1.0D0
    t_turn_tf_max = 0.05
    acs = 0.0D0
    cdtfleg = 0.0D0
    cforce = 0.0D0
    cplen = 0.0D0
    cpttf = 7.0e4
    cpttf_max = 9.0e4
    dcase = 8000.0D0
    dcond = (/6080.0D0, 6080.0D0, 6070.0D0, 6080.0D0, 6080.0D0, 8500.0D0, &
      6070.0D0, 8500.0D0, 8500.0D0/)
    dcondins = 1800.0D0
    dhecoil = 0.005D0
    estotftgj = 0.0D0
    b_crit_upper_nbti = 14.86D0
    t_crit_nbti = 9.04D0
    max_force_density = 0.0D0
    fcutfsu = 0.69D0
    fhts = 0.5D0
    insstrain = 0.0D0
    i_tf_stress_model = 1
    i_tf_tresca = 0
    i_tf_wp_geom = -1
    i_tf_case_geom = 0
    i_tf_turns_integer = 0
    i_tf_sc_mat = 1
    i_tf_sup = 1
    i_tf_shape = 0
    i_tf_cond_eyoung_axial = 0
    i_tf_cond_eyoung_trans = 1
    n_pancake = 10
    n_layer = 20
    n_rad_per_layer = 100
    i_tf_bucking = -1
    n_tf_graded_layers = 1
    n_tf_stress_layers = 0
    n_tf_wp_layers = 5
    jbus = 1.25D6
    jwdgcrt = 0.0D0
    jwdgpro = 0.0D0
    jwptf = 0.0D0
    oacdcp = 0.0D0
    eyoung_ins = 1.0D8
    eyoung_steel = 2.05D11
    eyoung_cond_axial = 6.6D8
    eyoung_cond_trans = 0.0D0
    eyoung_res_tf_buck = 150.0D9
    eyoung_copper = 117.0D9
    eyoung_al = 69.0D9
    poisson_steel = 0.3D0
    poisson_copper = 0.35D0
    poisson_al = 0.35D0
    poisson_ins = 0.34D0
    poisson_cond_axial = 0.3
    poisson_cond_trans = 0.3
    rbmax = 0.0D0
    tflegres = 0.0D0
    toroidalgap = 1.0D0 ![m]
    ftoroidalgap = 1.0D0
    ripmax = 1.0D0
    ripple = 0.0D0
    ritfc = 0.0D0
    radial_array = 0.0D0
    sig_tf_r = 0.0D0
    sig_tf_t = 0.0D0
    deflect = 0.0D0
    sig_tf_z = 0.0D0
    sig_tf_vmises = 0.0D0
    sig_tf_tresca = 0.0D0
    sig_tf_cs_bucked = 0.0D0
    sig_tf_case = 0.0D0
    sig_tf_wp = 0.0D0
    str_cs_con_res = -0.005D0
    str_pf_con_res = -0.005D0
    str_tf_con_res = -0.005D0
    str_wp = 0.0D0
    str_wp_max = 0.7D-2
    i_str_wp = 1
    quench_model = 'exponential'
    time1 = 0D0
    tcritsc = 16.0D0
    tdmptf = 10.0D0
    tfareain = 0.0D0
    tfbusl = 0.0D0
    tfbusmas = 0.0D0
    tfckw = 0.0D0
    tfcmw = 0.0D0
    tfcpmw = 0.0D0
    tfjtsmw = 0.0D0
    tfcryoarea = 0.0D0
    tficrn = 0.0D0
    tfind = 0.0D0
    tfinsgap = 0.010D0
    tflegmw = 0.0D0
    rhocp = 0.0D0
    rhotfleg = 0.0D0
    rhotfbus = -1.0D0 ! 2.5D-8
    frhocp = 1.0D0
    frholeg = 1.0D0
    rho_tf_joints = 2.5D-10
    n_tf_joints_contact = 6
    n_tf_joints = 4
    th_joint_contact = 0.03D0
    pres_joints = 0.0D0
    tfleng = 0.0D0
    eff_tf_cryo = -1.0D0
    n_tf = 16.0D0
    tfocrn = 0.0D0
    tfsai = 0.0D0
    tfsao = 0.0D0
    tftmp = 4.5D0
    tftort = 1.0D0
    thicndut = 8.0D-4
    layer_ins = 0.0D0
    thkcas = 0.3D0
    dr_tf_wp = 0.0D0
    thwcndut = 8.0D-3
    tinstf = 0.018D0
    tmargmin_tf = 0D0
    tmargmin_cs = 0D0
    tmargmin = 0D0
    temp_margin = 0.00D0
    tmargtf = 0.0D0
    tmaxpro = 150.0D0
    tmax_croco = 200.0D0
    croco_quench_temperature = 0D0
    tmpcry = 4.5D0
    n_tf_turn = 0.0D0
    vdalw = 20.0D0
    vforce = 0.0D0
    f_vforce_inboard = 0.5D0
    vforce_outboard = 0.0D0
    vftf = 0.4D0
    voltfleg = 0.0D0
    vtfkv = 0.0D0
    vtfskv = 0.0D0
    whtcas = 0.0D0
    whtcon = 0.0D0
    whtconcu = 0.0D0
    whtconal = 0.0D0
    whtconin = 0.0D0
    whtconsc = 0.0D0
    whtconsh = 0.0D0
    whtgw = 0.0D0
    whttf = 0.0D0
    wwp1 = 0.0D0
    wwp2 = 0.0D0
    dthet = 0.0D0
    radctf = 0.0D0
    xarc = 0.0D0
    xctfc = 0.0D0
    yarc = 0.0D0
    yctfc = 0.0D0
    tfa = 0.0D0
    tfb = 0.0D0
    drtop = 0.0D0
    dztop = 0.0D0
    etapump = 0.8D0
    fcoolcp = 0.3D0
    fcoolleg = 0.2D0
    a_cp_cool = 0.0D0
    ncool = 0.0D0
    ppump = 0.0D0
    prescp = 0.0D0
    presleg = 0.0D0
    ptempalw = 473.15D0   ! 200 C
    rcool = 0.005D0
    tcoolin = 313.15D0   ! 40 C
    dtiocool = 0.0D0
    tcpav = 373.15D0     ! 100 C
    tcpav2 = 0.0D0
    tlegav = -1.0D0
    tcpmax = 0.0D0
    vcool = 20.0D0
    vol_cond_cp = 0.0D0
    whtcp = 0.0D0
    whttflgs = 0.0D0
    tfc_sidewall_is_fraction = .false.
    i_cp_joints = -1
    cryo_cool_req = 0.0D0
    theta1_coil = 45.0D0
    theta1_vv = 1.0D0
    max_vv_stress = 143.0D6
  end subroutine init_tfcoil_variables
end module tfcoil_variables
