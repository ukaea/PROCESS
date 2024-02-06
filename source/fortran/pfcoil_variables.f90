module pfcoil_variables
  !! author: J. Morris, M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to the poloidal field coil systems
  !!
  !!### References
  !!
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  integer, parameter :: ngrpmx = 10
  !! maximum number of groups of PF coils

  integer, parameter :: nclsmx = 2
  !! maximum number of PF coils in a given group

  integer, parameter :: nptsmx = 32
  !! maximum number of points across the midplane of the plasma at which the field from
  !! the PF coils is fixed

  integer, parameter :: nfixmx = 64
  !! maximum number of fixed current PF coils

  integer, parameter :: ngc = ngrpmx*nclsmx
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

  real(dp) :: areaoh
  !! Central solenoid vertical cross-sectional area (m2)

  real(dp) :: a_oh_turn
  !! Central solenoid (OH) trun cross-sectional area (m2)

  real(dp) :: awpoh
  !! central solenoid conductor+void area with area of steel subtracted (m2)

  real(dp) :: bmaxoh
  !! maximum field in central solenoid at end of flat-top (EoF) (T)

  real(dp) :: bmaxoh0
  !! maximum field in central solenoid at beginning of pulse (T)

  real(dp), dimension(ngc2) :: bpf
  !! peak field at coil i (T)

  real(dp), dimension(ngrpmx) :: ccl0_ma
  !! PF group current array, flux-swing cancellation current (MA)
  !! Input if i_pf_current=0, computed otherwise

  real(dp), dimension(ngrpmx) :: ccls_ma
  !! PF group current array, equilibrium current (MA)
  !! Input if i_pf_current=0, computed otherwise

  real(dp) :: cohbop
  !! Central solenoid overall current density at beginning of pulse (A/m2)

  real(dp) :: coheof
  !! Central solenoid overall current density at end of flat-top (A/m2) (`iteration variable 37`) (`sweep variable 62`)

  real(dp), dimension(ngc2,6) :: cpt
  !! current per turn in coil i at time j (A)

  real(dp), dimension(ngc2) :: cptdin
  !! peak current per turn input for PF coil i (A)

  real(dp), dimension(ngc2) :: curpfb
  !! PF coil current array, at beginning of pulse (MA)
  !! Indexed by coil number, not group number

  real(dp), dimension(ngc2) :: curpff
  !! PF coil current array, at flat top (MA)
  !! Indexed by coil number, not group number

  real(dp), dimension(ngc2) :: curpfs
  !! PF coil current array, at end of pulse (MA)
  !! Indexed by coil number, not group number

  real(dp) :: etapsu
  !! Efficiency of transfer of PF stored energy into or out of storage.

  real(dp) :: fcohbof
  !! ratio of central solenoid overall current density at beginning of flat-top / end of flat-top

  real(dp) :: fcohbop
  !! ratio of central solenoid overall current density at beginning of pulse / end of flat-top
  !! (`iteration variable 41`)

  real(dp) :: fcuohsu
  !! copper fraction of strand in central solenoid

  real(dp) :: fcupfsu
  !! copper fraction of cable conductor (PF coils)

  real(dp) :: fvssu
  !! F-value for `constraint equation 51`

  integer, dimension(ngrpmx) :: ipfloc
  !! Switch for location of PF coil group i:
  !!
  !! - =1 PF coil on top of central solenoid (flux ramp only)
  !! - =2 PF coil on top of TF coil (flux ramp only)
  !! - =3 PF coil outside of TF coil (equilibrium coil)
  !! - =4 PF coil, general location (equilibrium coil)

  integer :: ipfres
  !! switch for PF & CS coil conductor type:
  !!
  !! - =0 superconducting PF coils
  !! - =1 resistive PF coils
  !
  real(dp) :: itr_sum
  !! total sum of I x turns x radius for all PF coils and CS (Am)

  integer :: isumatoh
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

  integer :: isumatpf
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

  integer :: i_pf_current
  !! Switch for controlling the current of the PF coils:
  !!
  !! - =0 Input via the variables curpfb, curpff, curpfs
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

  real(dp) :: jscoh_bop
  !! central solenoid superconductor critical current density (A/m2) at beginning-of-pulse

  real(dp) :: jscoh_eof
  !! central solenoid superconductor critical current density (A/m2) at end-of-flattop

  real(dp) :: jstrandoh_bop
  !! central solenoid strand critical current density (A/m2) at beginning-of-pulse

  real(dp) :: jstrandoh_eof
  !! central solenoid strand critical current density (A/m2) at end-of-flattop

  integer :: ncirt
  !! number of PF circuits (including central solenoid and plasma)

  integer, dimension(ngrpmx+2) :: ncls
  !! number of PF coils in group j

  integer :: nfxfh
  !! number of filaments the top and bottom of the central solenoid should be broken
  !! into during scaling (5 - 10 is good)

  integer :: ngrp
  !! number of groups of PF coils. Symmetric coil pairs should all be in the same group

  integer :: nohc
  !! number of PF coils (excluding the central solenoid) + 1

  real(dp) :: ohhghf
  !! Central solenoid height / TF coil internal height

  real(dp) :: oh_steel_frac
  !! central solenoid steel fraction (`iteration variable 122`)

  real(dp) :: pf_current_safety_factor
  !! Ratio of permissible PF coil conductor current density to critical conductor
  !! current density based on short-sample DC measurements

  real(dp), dimension(ngc2) :: pfcaseth
  !! steel case thickness for PF coil i (m)

  real(dp) :: pfclres
  !! PF coil resistivity (if ipfres=1) (Ohm-m)

  real(dp) :: pfmmax
  !! mass of heaviest PF coil (tonnes)

  real(dp) :: pfrmax
  !! radius of largest PF coil (m)

  real(dp) :: pfwpmw
  !! Total mean wall plug power dissipated in PFC and CS power supplies (MW) (issue #713)

  real(dp) :: powohres
  !! central solenoid resistive power during flattop (W)

  real(dp) :: powpfres
  !! total PF coil resistive losses during flattop (W)

  real(dp), dimension(ngc2) :: ra
  !! inner radius of coil i (m)

  real(dp), dimension(ngc2) :: rb
  !! outer radius of coil i (m)

  real(dp), dimension(ngc2) :: ric
  !! peak current in coil i (MA-turns)

  real(dp), dimension(ngc2) :: rjconpf
  !! average winding pack current density of PF coil i (A/m2) at time of peak
  !! current in that coil (calculated for `ipfloc=1` coils)

  real(dp) :: rjohc
  !! allowable central solenoid current density at end of flat-top (A/m2)

  real(dp) :: rjohc0
  !! allowable central solenoid current density at beginning of pulse (A/m2)

  real(dp), dimension(ngc2) :: rjpfalw
  !! allowable winding pack current density of PF coil i (A/m2)

  real(dp) :: rohc
  !! radius to the centre of the central solenoid (m)

  real(dp) :: routr
  !! radial distance (m) from outboard TF coil leg to centre of `ipfloc=3` PF coils

  real(dp), dimension(ngc2) :: rpf
  !! radius of PF coil i (m)

  real(dp) :: rpf1
  !! offset (m) of radial position of `ipfloc=1` PF coils from being directly above
  !! the central solenoid

  real(dp) :: rpf2
  !! offset (m) of radial position of `ipfloc=2` PF coils from being at
  !! rmajor (offset = rpf2*triang*rminor)

  real(dp), dimension(ngrpmx) :: rref
  !! PF coil radial positioning adjuster:
  !!
  !! - for groups j with ipfloc(j) = 1; rref(j) is ignored
  !! - for groups j with ipfloc(j) = 2; rref(j) is ignored
  !! - for groups j with ipfloc(j) = 3; rref(j) is ignored
  !! - for groups j with ipfloc(j) = 4; rref(j) is radius of
  !!   the coil in units of minor radii from the major radius
  !!   (r = rmajor + rref*rminor)

  real(dp) :: s_tresca_oh
  !! Maximum shear stress (Tresca criterion) coils/central solenoid [MPa]

  real(dp) :: sigpfcalw
  !! maximum permissible tensile stress (MPa) in steel coil cases for superconducting
  !! PF coils (`ipfres=0`)

  real(dp) :: sigpfcf
  !! fraction of JxB hoop force supported by steel case for superconducting PF coils (`ipfres=0`)

  real(dp), dimension(ngc2,ngc2) :: sxlg
  !! mutual inductance matrix (H)

  real(dp) :: tmargoh
  !! Central solenoid temperature margin (K)

  real(dp), dimension(ngc2) :: turns
  !! number of turns in PF coil i

  real(dp), dimension(ngc2) :: vf
  !! winding pack void fraction of PF coil i for coolant

  real(dp) :: vfohc
  !! void fraction of central solenoid conductor for coolant

  real(dp) :: vsbn
  !! total flux swing available for burn (Wb)

  real(dp) :: vsefbn
  !! flux swing from PF coils for burn (Wb)

  real(dp) :: vsefsu
  !! flux swing from PF coils for startup (Wb)

  real(dp) :: vseft
  !! total flux swing from PF coils (Wb)

  real(dp) :: vsoh
  !! total flux swing from the central solenoid (Wb)

  real(dp) :: vsohbn
  !! central solenoid flux swing for burn (Wb)

  real(dp) :: vsohsu
  !! central solenoid flux swing for startup (Wb)

  real(dp) :: vssu
  !! total flux swing for startup (`constraint eqn 51` to enforce vssu=vsres+vsind) (Wb)

  real(dp) :: vstot
  !! total flux swing for pulse (Wb)

  real(dp), dimension(ngc2,6) :: waves
  !! used in current waveform of PF coils/central solenoid

  real(dp) :: whtpf
  !! total mass of the PF coil conductor (kg)

  real(dp) :: whtpfs
  !! total mass of the PF coil structure (kg)

  real(dp), dimension(ngc2) :: wtc
  !! conductor mass for PF coil i (kg)

  real(dp), dimension(ngc2) :: wts
  !! structure mass for PF coil i (kg)

  real(dp), dimension(ngc2) :: zh
  !! upper point of PF coil i (m)

  real(dp), dimension(ngc2) :: zl
  !! lower point of PF coil i (m)

  real(dp), dimension(ngc2) :: zpf
  !! z (height) location of PF coil i (m)

  real(dp), dimension(ngrpmx) :: zref
  !! PF coil vertical positioning adjuster:
  !!
  !! - for groups j with ipfloc(j) = 1; zref(j) is ignored
  !! - for groups j with ipfloc(j) = 2 AND itart=1 (only);
  !!   zref(j) is distance of centre of PF coil from inside
  !!   edge of TF coil (remember that PF coils for STs lie
  !!   within the TF coil)
  !! - for groups j with ipfloc(j) = 3; zref(j) = ratio of
  !!   height of coil group j to plasma minor radius</UL>
  !! - for groups j with ipfloc(j) = 4; zref(j) = ratio of
  !!   height of coil group j to plasma minor radius</UL>

  real(dp) :: bmaxcs_lim
  !! Central solenoid max field limit [T]

  real(dp) :: fbmaxcs
  !! F-value for CS mmax field (`cons. 79`, `itvar 149`)

  real(dp) :: ld_ratio_cst
  !! Ratio of CS coil turn conduit length to depth

  real(dp) :: l_cond_cst
  !! Length of CS of CS coil turn conduit

  real(dp) :: d_cond_cst
  !! Depth/width of CS of CS coil turn conduit

  real(dp) :: r_out_cst
  !! Length of CS of CS coil turn conduit length

  real(dp) :: r_in_cst
  !! Length of CS of CS coil turn conduit length

  contains

  subroutine init_pfcoil_variables
    !! Initialise module variables
    implicit none

    alfapf = 5.0D-10
    alstroh = 4.0D8
    i_cs_stress = 0
    areaoh = 0.0D0
    a_oh_turn = 0.0D0
    awpoh = 0.0D0
    bmaxoh = 0.0D0
    bmaxoh0 = 0.0D0
    bpf = 0.0D0
    ccl0_ma = (/0.0D0, 0.0D0, 0.0D0, &
      0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
    ccls_ma = (/0.0D0, 0.0D0, 0.0D0, &
      0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0/)
    cohbop = 0.0D0
    coheof = 1.85D7
    cpt = 0.0D0
    cptdin = 4.0D4
    curpfb = 0.0D0
    curpff = 0.0D0
    curpfs = 0.0D0
    etapsu = 0.9D0
    fcohbof = 0.0D0
    fcohbop = 0.9D0
    fcuohsu = 0.7D0
    fcupfsu = 0.69D0
    fvssu = 1.0
    ipfloc = (/2,2,3,0,0,0,0,0,0,0/)
    ipfres = 0
    itr_sum = 0.0D0
    isumatoh = 1
    isumatpf = 1
    i_pf_current = 1
    i_sup_pf_shape = 0
    jscoh_bop = 0.0D0
    jscoh_eof = 0.0D0
    jstrandoh_bop = 0.0D0
    jstrandoh_eof = 0.0D0
    ncirt = 0
    ncls = (/1,1,2,0,0,0,0,0,0,0,0,0/)
    nfxfh = 7
    ngrp = 3
    nohc = 0
    ohhghf = 0.71D0
    oh_steel_frac = 0.5D0
    pf_current_safety_factor = 1.0D0
    pfcaseth = 0.0D0
    pfclres = 2.5D-8
    pfmmax = 0.0D0
    pfrmax = 0.0D0
    pfwpmw = 0.0D0
    powohres = 0.0D0
    powpfres = 0.0D0
    ra = 0.0D0
    rb = 0.0D0
    ric = 0.0D0
    rjconpf = 3.0D7
    rjohc = 0.0D0
    rjohc0 = 0.0D0
    rjpfalw = 0.0D0
    rohc = 0.0D0
    routr = 1.5D0
    rpf = 0.0D0
    rpf1 = 0.0D0
    rpf2 = -1.63D0
    rref = (/7.0D0, 7.0D0, 7.0D0, &
      7.0D0, 7.0D0, 7.0D0, 7.0D0, 7.0D0, 7.0D0, 7.0D0/)
    s_tresca_oh = 0.0D0
    sigpfcalw = 500.0D0
    sigpfcf = 1.0D0
    sxlg = 0.0D0
    tmargoh = 0.0D0
    turns = 0.0D0
    vf = 0.3D0
    vfohc = 0.3D0
    vsbn = 0.0D0
    vsefbn = 0.0D0
    vsefsu = 0.0D0
    vseft = 0.0D0
    vsoh = 0.0D0
    vsohbn = 0.0D0
    vsohsu = 0.0D0
    vssu = 0.0D0
    vstot = 0.0D0
    waves = 0.0D0
    whtpf = 0.0D0
    whtpfs = 0.0D0
    wtc = 0.0D0
    wts = 0.0D0
    zh = 0.0D0
    zl = 0.0D0
    zpf = 0.0D0
    zref = (/3.6D0, 1.2D0, 2.5D0, &
      1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0/)
    bmaxcs_lim = 13.0
    fbmaxcs = 1.0
    ld_ratio_cst = 70/22
    l_cond_cst = 0.0D0
    d_cond_cst = 0.0D0
    r_in_cst = 0.0D0
    r_out_cst = 3.0D-3
  end subroutine init_pfcoil_variables
end module pfcoil_variables
