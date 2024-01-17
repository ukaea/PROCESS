! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module numerics
  !! Module containing callers to the main equation solvers
  !! HYBRD and VMCON
  !! author: P J Knight, CCFE, Culham Science Centre
  !! This module contains the primary numerics variables and the
  !! calling routines for the two equation solvers in the code.
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  public

  integer, parameter :: ipnvars = 175
  !!  ipnvars FIX : total number of variables available for iteration

  integer, parameter :: ipeqns = 91
  !!  ipeqns  FIX : number of constraint equations available

  integer, parameter :: ipnfoms = 19
  !!  ipnfoms FIX : number of available figures of merit

  integer, parameter :: ipvlam  = ipeqns+2*ipnvars+1
  integer, parameter :: iptnt   = (ipeqns*(3*ipeqns+13))/2
  integer, parameter :: ipvp1   = ipnvars+1

  integer :: ioptimz
  !!  ioptimz /1/ : code operation switch:<UL>
  !!           <LI> = -2 for no optimisation, no VMCOM or HYBRD;
  !!           <LI> = -1 for no optimisation, HYBRD only;
  !!           <LI> = 0  for HYBRD and VMCON (not recommended);
  !!           <LI> = 1  for optimisation, VMCON only</UL>

  !!  minmax /7/ : switch for figure-of-merit (see lablmm for descriptions)
  !!               negative => maximise, positive => minimise
  integer :: minmax
  character(len=22), dimension(ipnfoms) :: lablmm
  !!  lablmm(ipnfoms) : labels describing figures of merit:<UL>
  !!  <LI> ( 1) major radius
  !!  <LI> ( 2) not used
  !!  <LI> ( 3) neutron wall load
  !!  <LI> ( 4) P_tf + P_pf
  !!  <LI> ( 5) fusion gain Q
  !!  <LI> ( 6) cost of electricity
  !!  <LI> ( 7) capital cost (direct cost if ireactor=0,
  !!                          constructed cost otherwise)
  !!  <LI> ( 8) aspect ratio
  !!  <LI> ( 9) divertor heat load
  !!  <LI> (10) toroidal field
  !!  <LI> (11) total injected power
  !!  <LI> (12) hydrogen plant capital cost OBSOLETE
  !!  <LI> (13) hydrogen production rate OBSOLETE
  !!  <LI> (14) pulse length
  !!  <LI> (15) plant availability factor (N.B. requires
  !!            iavail=1 to be set)
  !!  <LI> (16) linear combination of major radius (minimised) and pulse length (maximised)
  !!              note: FoM should be minimised only!
  !!  <LI> (17) net electrical output
  !!  <LI> (18) Null Figure of Merit
  !!  <LI> (19) linear combination of big Q and pulse length (maximised)
  !!              note: FoM should be minimised only!</UL>

  integer :: ncalls
  !!  ncalls : number of function calls during solution
  integer :: neqns
  !!  neqns /0/ : number of equality constraints to be satisfied
  integer :: nfev1
  !!  nfev1 : number of calls to FCNHYB (HYBRD function caller) made
  integer :: nfev2
  !!  nfev2 : number of calls to FCNVMC1 (VMCON function caller) made
  integer :: nineqns
  !!  nineqns /0/ : number of inequality constraints VMCON must satisfy
  !!                (leave at zero for now)
  integer :: nvar
  !!  nvar /16/ : number of iteration variables to use
  integer :: nviter
  !!  nviter : number of VMCON iterations performed

  !!  icc(ipeqns) /0/ :
  !!           array defining which constraint equations to activate
  !!           (see lablcc for descriptions)

  ! TODO Check the dictionaries are created correctly.
  ! Issue #491 Default constraints removed.
  integer, dimension(ipeqns) :: icc

  logical, dimension(ipeqns) :: active_constraints
  !!  active_constraints(ipeqns) : Logical array showing which constraints are active

  ! #TODO Do not change the comments for lablcc: they are used to create the
  ! Python-Fortran dictionaries. This must be improved on.
  character(len=33), dimension(ipeqns) :: lablcc
  !!  lablcc(ipeqns) : labels describing constraint equations (corresponding itvs)<UL>
  !!  <LI> ( 1) Beta (consistency equation) (itv 5)
  !!  <LI> ( 2) Global power balance (consistency equation) (itv 10,1,2,3,4,6,11)
  !!  <LI> ( 3) Ion power balance DEPRECATED (itv 10,1,2,3,4,6,11)
  !!  <LI> ( 4) Electron power balance DEPRECATED (itv 10,1,2,3,4,6,11)
  !!  <LI> ( 5) Density upper limit (itv 9,1,2,3,4,5,6)
  !!  <LI> ( 6) (Epsilon x beta poloidal) upper limit (itv 8,1,2,3,4,6)
  !!  <LI> ( 7) Beam ion density (NBI) (consistency equation) (itv 7)
  !!  <LI> ( 8) Neutron wall load upper limit (itv 14,1,2,3,4,6)
  !!  <LI> ( 9) Fusion power upper limit (itv 26,1,2,3,4,6)
  !!  <LI> (10) Toroidal field 1/R (consistency equation) (itv 12,1,2,3,13 )
  !!  <LI> (11) Radial build (consistency equation) (itv 3,1,13,16,29,42,61)
  !!  <LI> (12) Volt second lower limit (STEADY STATE) (itv 15,1,2,3)
  !!  <LI> (13) Burn time lower limit (PULSE) (itv 21,1,16,17,29,42,44,61)
  !!            (itv 19,1,2,3,6)
  !!  <LI> (14) Neutral beam decay lengths to plasma centre (NBI) (consistency equation)
  !!  <LI> (15) LH power threshold limit (itv 103)
  !!  <LI> (16) Net electric power lower limit (itv 25,1,2,3)
  !!  <LI> (17) Radiation fraction upper limit (itv 28)
  !!  <LI> (18) Divertor heat load upper limit (itv 27)
  !!  <LI> (19) MVA upper limit (itv 30)
  !!  <LI> (20) Neutral beam tangency radius upper limit (NBI) (itv 33,31,3,13)
  !!  <LI> (21) Plasma minor radius lower limit (itv 32)
  !!  <LI> (22) Divertor collisionality upper limit (itv 34,43)
  !!  <LI> (23) Conducting shell to plasma minor radius ratio upper limit
  !!            (itv 104,1,74)
  !!  <LI> (24) Beta upper limit (itv 36,1,2,3,4,6,18)
  !!  <LI> (25) Peak toroidal field upper limit (itv 35,3,13,29)
  !!  <LI> (26) Central solenoid EOF current density upper limit (ipfres=0)
  !!            (itv 38,37,41,12)
  !!  <LI> (27) Central solenoid BOP current density upper limit (ipfres=0)
  !!            (itv 39,37,41,12)
  !!  <LI> (28) Fusion gain Q lower limit (itv 45,47,40)
  !!  <LI> (29) Inboard radial build consistency (itv 3,1,13,16,29,42,61)
  !!  <LI> (30) Injection power upper limit (itv 46,47,11)
  !!  <LI> (31) TF coil case stress upper limit (SCTF) (itv 48,56,57,58,59,60,24)
  !!  <LI> (32) TF coil conduit stress upper limit (SCTF) (itv 49,56,57,58,59,60,24)
  !!  <LI> (33) I_op / I_critical (TF coil) (SCTF) (itv 50,56,57,58,59,60,24)
  !!  <LI> (34) Dump voltage upper limit (SCTF) (itv 51,52,56,57,58,59,60,24)
  !!  <LI> (35) J_winding pack/J_protection upper limit (SCTF) (itv 53,56,57,58,59,60,24)
  !!  <LI> (36) TF coil temperature margin lower limit (SCTF) (itv 54,55,56,57,58,59,60,24)
  !!  <LI> (37) Current drive gamma upper limit (itv 40,47)
  !!  <LI> (38) First wall coolant temperature rise upper limit (itv 62)
  !!  <LI> (39) First wall peak temperature upper limit (itv 63)
  !!  <LI> (40) Start-up injection power lower limit (PULSE) (itv 64)
  !!  <LI> (41) Plasma current ramp-up time lower limit (PULSE) (itv  66,65)
  !!  <LI> (42) Cycle time lower limit (PULSE) (itv 17,67,65)
  !!  <LI> (43) Average centrepost temperature
  !!            (TART) (consistency equation) (itv 13,20,69,70)
  !!  <LI> (44) Peak centrepost temperature upper limit (TART) (itv 68,69,70)
  !!  <LI> (45) Edge safety factor lower limit (TART) (itv 71,1,2,3)
  !!  <LI> (46) Equation for Ip/Irod upper limit (TART) (itv 72,2,60)
  !!  <LI> (47) NOT USED
  !!  <LI> (48) Poloidal beta upper limit (itv 79,2,3,18)
  !!  <LI> (49) NOT USED
  !!  <LI> (50) IFE repetition rate upper limit (IFE)
  !!  <LI> (51) Startup volt-seconds consistency (PULSE) (itv 16,29,3,1)
  !!  <LI> (52) Tritium breeding ratio lower limit (itv 89,90,91)
  !!  <LI> (53) Neutron fluence on TF coil upper limit (itv 92,93,94)
  !!  <LI> (54) Peak TF coil nuclear heating upper limit (itv 95,93,94)
  !!  <LI> (55) Vacuum vessel helium concentration upper limit iblanket =2 (itv 96,93,94)
  !!  <LI> (56) Pseparatrix/Rmajor upper limit (itv 97,1,3)
  !!  <LI> (57) NOT USED
  !!  <LI> (58) NOT USED
  !!  <LI> (59) Neutral beam shine-through fraction upper limit (NBI) (itv 105,6,19,4 )
  !!  <LI> (60) Central solenoid temperature margin lower limit (SCTF) (itv 106)
  !!  <LI> (61) Minimum availability value (itv 107)
  !!  <LI> (62) taup/taueff the ratio of particle to energy confinement times (itv 110)
  !!  <LI> (63) The number of ITER-like vacuum pumps niterpump < tfno (itv 111)
  !!  <LI> (64) Zeff less than or equal to zeffmax (itv 112)
  !!  <LI> (65) Dump time set by VV loads (itv 56, 113)
  !!  <LI> (66) Limit on rate of change of energy in poloidal field
  !!            (Use iteration variable 65(tohs), 115)
  !!  <LI> (67) Simple Radiation Wall load limit (itv 116, 4,6)
  !!  <LI> (68) Psep * Bt / qAR upper limit (itv 117)
  !!  <LI> (69) ensure separatrix power = the value from Kallenbach divertor (itv 118)
  !!  <LI> (70) ensure that teomp = separatrix temperature in the pedestal profile,
  !!            (itv 119 (tesep))
  !!  <LI> (71) ensure that neomp = separatrix density (nesep) x neratio
  !!  <LI> (72) central solenoid shear stress limit (Tresca yield criterion) (itv 123 foh_stress)
  !!  <LI> (73) Psep >= Plh + Paux (itv 137 (fplhsep))
  !!  <LI> (74) TFC quench < tmax_croco (itv 141 (fcqt))
  !!  <LI> (75) TFC current/copper area < Maximum (itv 143 f_coppera_m2)
  !!  <LI> (76) Eich critical separatrix density
  !!  <LI> (77) TF coil current per turn upper limit
  !!  <LI> (78) Reinke criterion impurity fraction lower limit (itv  147 freinke)
  !!  <LI> (79) Peak CS field upper limit (itv  149 fbmaxcs)
  !!  <LI> (80) Divertor power lower limit pdivt (itv  153 fpdivlim)
  !!  <LI> (81) Ne(0) > ne(ped) constraint (itv  154 fne0)
  !!  <LI> (82) toroidalgap >  tftort constraint (itv  171 ftoroidalgap)
  !!  <LI> (83) Radial build consistency for stellarators (itv 172 f_avspace)
  !!  <LI> (84) Lower limit for beta (itv 173 fbetatry_lower)
  !!  <LI> (85) Constraint for CP lifetime
  !!  <LI> (86) Constraint for TF coil turn dimension
  !!  <LI> (87) Constraint for cryogenic power
  !!  <LI> (88) Constraint for TF coil strain absolute value
  !!  <LI> (89) Constraint for CS coil quench protection
  !!  <LI> (90) Lower Limit on number of stress load cycles for CS (itr 167 fncycle)
  !!  <LI> (91) Checking if the design point is ECRH ignitable (itv 168 fecrh_ignition) </UL>

  integer, dimension(ipnvars) :: ixc
  !!  ixc(ipnvars) /0/ :
  !!               array defining which iteration variables to activate
  !!               (see lablxc for descriptions)

  character*14, dimension(ipnvars) :: lablxc
  !! lablxc(ipnvars) : labels describing iteration variables<UL>
  !!  <LI> ( 1) aspect
  !!  <LI> ( 2) bt
  !!  <LI> ( 3) rmajor
  !! <LI> ( 4) te
  !! <LI> ( 5) beta
  !! <LI> ( 6) dene
  !! <LI> ( 7) rnbeam
  !! <LI> ( 8) fbeta (f-value for equation 6)
  !! <LI> ( 9) fdene (f-value for equation 5)
  !! <LI> (10) hfact
  !! <LI> (11) pheat
  !! <LI> (12) oacdcp
  !! <LI> (13) tfcth (NOT RECOMMENDED)
  !! <LI> (14) fwalld (f-value for equation 8)
  !! <LI> (15) fvs (f-value for equation 12)
  !! <LI> (16) ohcth
  !! <LI> (17) tdwell
  !! <LI> (18) q
  !! <LI> (19) enbeam
  !! <LI> (20) tcpav
  !! <LI> (21) ftburn (f-value for equation 13)
  !! <LI> (22) NOT USED
  !! <LI> (23) fcoolcp
  !! <LI> (24) NOT USED
  !! <LI> (25) fpnetel (f-value for equation 16)
  !! <LI> (26) ffuspow (f-value for equation 9)
  !! <LI> (27) fhldiv (f-value for equation 18)
  !! <LI> (28) fradpwr (f-value for equation 17), total radiation fraction
  !! <LI> (29) bore
  !! <LI> (30) fmva (f-value for equation 19)
  !! <LI> (31) gapomin
  !! <LI> (32) frminor (f-value for equation 21)
  !! <LI> (33) fportsz (f-value for equation 20)
  !! <LI> (34) fdivcol (f-value for equation 22)
  !! <LI> (35) fpeakb (f-value for equation 25)
  !! <LI> (36) fbetatry (f-value for equation 24)
  !! <LI> (37) coheof
  !! <LI> (38) fjohc (f-value for equation 26)
  !! <LI> (39) fjohc0 (f-value for equation 27)
  !! <LI> (40) fgamcd (f-value for equation 37)
  !! <LI> (41) fcohbop
  !! <LI> (42) gapoh
  !! <LI> (43) NOT USED
  !! <LI> (44) fvsbrnni
  !! <LI> (45) fqval (f-value for equation 28)
  !! <LI> (46) fpinj (f-value for equation 30)
  !! <LI> (47) feffcd
  !! <LI> (48) fstrcase (f-value for equation 31)
  !! <LI> (49) fstrcond (f-value for equation 32)
  !! <LI> (50) fiooic (f-value for equation 33)
  !! <LI> (51) fvdump (f-value for equation 34)
  !! <LI> (52) vdalw
  !! <LI> (53) fjprot (f-value for equation 35)
  !! <LI> (54) ftmargtf (f-value for equation 36)
  !! <LI> (55) NOT USED
  !! <LI> (56) tdmptf
  !! <LI> (57) thkcas
  !! <LI> (58) thwcndut
  !! <LI> (59) fcutfsu
  !! <LI> (60) cpttf
  !! <LI> (61) gapds
  !! <LI> (62) fdtmp (f-value for equation 38)
  !! <LI> (63) ftpeak (f-value for equation 39)
  !! <LI> (64) fauxmn (f-value for equation 40)
  !! <LI> (65) tohs
  !! <LI> (66) ftohs (f-value for equation 41)
  !! <LI> (67) ftcycl (f-value for equation 42)
  !! <LI> (68) fptemp (f-value for equation 44)
  !! <LI> (69) rcool
  !! <LI> (70) vcool
  !! <LI> (71) fq (f-value for equation 45)
  !! <LI> (72) fipir (f-value for equation 46)
  !! <LI> (73) scrapli
  !! <LI> (74) scraplo
  !! <LI> (75) tfootfi
  !! <LI> (76) NOT USED
  !! <LI> (77) NOT USED
  !! <LI> (78) NOT USED
  !! <LI> (79) fbetap (f-value for equation 48)
  !! <LI> (80) NOT USED
  !! <LI> (81) edrive
  !! <LI> (82) drveff
  !! <LI> (83) tgain
  !! <LI> (84) chrad
  !! <LI> (85) pdrive
  !! <LI> (86) frrmax (f-value for equation 50)
  !! <LI> (87) NOT USED
  !! <LI> (88) NOT USED
  !! <LI> (89) ftbr (f-value for equation 52)
  !! <LI> (90) blbuith
  !! <LI> (91) blbuoth
  !! <LI> (92) fflutf (f-value for equation 53)
  !! <LI> (93) shldith
  !! <LI> (94) shldoth
  !! <LI> (95) fptfnuc (f-value for equation 54)
  !! <LI> (96) fvvhe (f-value for equation 55)
  !! <LI> (97) fpsepr (f-value for equation 56)
  !! <LI> (98) li6enrich
  !! <LI> (99) NOT USED
  !! <LI> (100) NOT USED
  !! <LI> (101) NOT USED
  !! <LI> (102) fimpvar # OBSOLETE
  !! <LI> (103) flhthresh (f-value for equation 15)
  !! <LI> (104) fcwr (f-value for equation 23)
  !! <LI> (105) fnbshinef (f-value for equation 59)
  !! <LI> (106) ftmargoh (f-value for equation 60)
  !! <LI> (107) favail (f-value for equation 61)
  !! <LI> (108) breeder_f: Volume of Li4SiO4 / (Volume of Be12Ti + Li4SiO4)
  !! <LI> (109) ralpne: thermal alpha density / electron density
  !! <LI> (110) ftaulimit: Lower limit on taup/taueff the ratio of alpha
  !! <LI> (111) fniterpump: f-value for constraint that number
  !! <LI> (112) fzeffmax: f-value for max Zeff (f-value for equation 64)
  !! <LI> (113) ftaucq: f-value for minimum quench time (f-value for equation 65)
  !! <LI> (114) fw_channel_length: Length of a single first wall channel
  !! <LI> (115) fpoloidalpower: f-value for max rate of change of
  !! <LI> (116) fradwall: f-value for radiation wall load limit (eq. 67)
  !! <LI> (117) fpsepbqar: f-value for  Psep*Bt/qar upper limit (eq. 68)
  !! <LI> (118) fpsep: f-value to ensure separatrix power is less than
  !! <LI> (119) tesep:  separatrix temperature calculated by the Kallenbach divertor model
  !! <LI> (120) ttarget: Plasma temperature adjacent to divertor sheath [eV]
  !! <LI> (121) neratio: ratio of mean SOL density at OMP to separatrix density at OMP
  !! <LI> (122) oh_steel_frac : streel fraction of Central Solenoid
  !! <LI> (123) foh_stress : f-value for CS coil Tresca yield criterion (f-value for eq. 72)
  !! <LI> (124) qtargettotal : Power density on target including surface recombination [W/m2]
  !! <LI> (125) fimp(3) :  Beryllium density fraction relative to electron density
  !! <LI> (126) fimp(4) :  Carbon density fraction relative to electron density
  !! <LI> (127) fimp(5) :  Nitrogen fraction relative to electron density
  !! <LI> (128) fimp(6) :  Oxygen density fraction relative to electron density
  !! <LI> (129) fimp(7) :  Neon density fraction relative to electron density
  !! <LI> (130) fimp(8) :  Silicon density fraction relative to electron density
  !! <LI> (131) fimp(9) :  Argon density fraction relative to electron density
  !! <LI> (132) fimp(10) :  Iron density fraction relative to electron density
  !! <LI> (133) fimp(11) :  Nickel density fraction relative to electron density
  !! <LI> (134) fimp(12) :  Krypton density fraction relative to electron density
  !! <LI> (135) fimp(13) :  Xenon density fraction relative to electron density
  !! <LI> (136) fimp(14) :  Tungsten density fraction relative to electron density
  !! <LI> (137) fplhsep (f-value for equation 73)
  !! <LI> (138) rebco_thickness : thickness of REBCO layer in tape (m)
  !! <LI> (139) copper_thick : thickness of copper layer in tape (m)
  !! <LI> (140) dr_tf_wp : radial thickness of TFC winding pack (m)
  !! <LI> (141) fcqt : TF coil quench temperature < tmax_croco (f-value for equation 74)
  !! <LI> (142) nesep : electron density at separatrix [m-3]
  !! <LI> (143) f_copperA_m2 : TF coil current / copper area < Maximum value
  !! <LI> (144) fnesep : Eich critical electron density at separatrix
  !! <LI> (145) fgwped :  fraction of Greenwald density to set as pedestal-top density
  !! <LI> (146) fcpttf : F-value for TF coil current per turn limit (constraint equation 77)
  !! <LI> (147) freinke : F-value for Reinke detachment criterion (constraint equation 78)
  !! <LI> (148) fzactual : fraction of impurity at SOL with Reinke detachment criterion
  !! <LI> (149) fbmaxcs : F-value for max peak CS field (con. 79, itvar 149)
  !! <LI> (150) REMOVED
  !! <LI> (151) REMOVED
  !! <LI> (152) fbmaxcs : Ratio of separatrix density to Greenwald density
  !! <LI> (153) fpdivlim : F-value for minimum pdivt (con. 80)
  !! <LI> (154) fne0 : F-value for ne(0) > ne(ped) (con. 81)
  !! <LI> (155) pfusife : IFE input fusion power (MW) (ifedrv=3 only)
  !! <LI> (156) rrin : Input IFE repetition rate (Hz) (ifedrv=3 only)
  !! <LI> (157) fvssu : F-value for available to required start up flux (con. 51)
  !! <LI> (158) croco_thick : Thickness of CroCo copper tube (m)
  !! <LI> (159) ftoroidalgap : F-value for toroidalgap >  tftort constraint (con. 82)
  !! <LI> (160) f_avspace (f-value for equation 83)
  !! <LI> (161) fbetatry_lower (f-value for equation 84)
  !! <LI> (162) r_cp_top : Top outer radius of the centropost (ST only) (m)
  !! <LI> (163) f_t_turn_tf : f-value for TF coils WP trurn squared dimension constraint
  !! <LI> (164) f_crypmw : f-value for cryogenic plant power
  !! <LI> (165) fstr_wp : f-value for TF coil strain absolute value
  !! <LI> (166) f_copperaoh_m2 : CS coil current /copper area < Maximum value
  !! <LI> (167) fncycle : f-value for minimum CS coil stress load cycles
  !! <LI> (168) fecrh_ignition: f-value for equation 91
  !! <LI> (169) te0_ecrh_achievable: Max. achievable electron temperature at ignition point
  !! <LI> (170) beta_div : field line angle wrt divertor target plate (degrees)
  !! <LI> (171) EMPTY : Description
  !! <LI> (172) EMPTY : Description
  !! <LI> (173) EMPTY : Description
  !! <LI> (174) EMPTY : Description
  !! <LI> (175) EMPTY : Description
  ! Issue 287 iteration variables are now defined in module define_iteration_variables in iteration variables.f90

  character(len=14), dimension(:), allocatable :: name_xc

  real(dp) :: sqsumsq
  !!  sqsumsq : sqrt of the sum of the square of the constraint residuals
  real(dp) :: norm_objf
  !! Normalised objective function (figure of merit)
  real(dp) :: epsfcn
  !!  epsfcn /1.0e-3/ : finite difference step length for HYBRD/VMCON derivatives
  real(dp) :: epsvmc
  !!  epsvmc /1.0e-6/ : error tolerance for VMCON
  real(dp) :: factor
  !!  factor /0.1/ : used in HYBRD for first step size
  real(dp) :: ftol
  !!  ftol /1.0e-4/ : error tolerance for HYBRD

  real(dp), dimension(ipnvars) :: boundl
  !!  boundl(ipnvars) /../ : lower bounds used on ixc variables during
  !!                         VMCON optimisation runs

  ! Issue #287 These bounds now defined in initial.f90
  real(dp), dimension(ipnvars) :: boundu
  ! !!  boundu(ipnvars) /../ : upper bounds used on ixc variables

  real(dp), dimension(ipnvars) :: bondl
  real(dp), dimension(ipnvars) :: bondu
  real(dp), dimension(ipnvars) :: rcm
  real(dp), dimension(ipnvars) :: resdl
  real(dp), dimension(ipnvars) :: scafc
  real(dp), dimension(ipnvars) :: scale
  real(dp), dimension(ipnvars) :: xcm
  real(dp), dimension(ipnvars) :: xcs
  real(dp), dimension(ipvlam)  :: vlam

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_numerics()
    ! Initialise module variables
    ioptimz = 1
    minmax = 7
    lablmm = (/ &
      'major radius.         ', &
      'not used.             ', &
      'neutron wall load.    ', &
      'P_tf + P_pf.          ', &
      'fusion gain.          ', &
      'cost of electricity.  ', &
      'capital cost.         ', &
      'aspect ratio.         ', &
      'divertor heat load.   ', &
      'toroidal field.       ', &
      'total injected power. ', &
      'H plant capital cost. ', &
      'H production rate.    ', &
      'pulse length.         ', &
      'plant availability.   ', &
      'min R0, max tau_burn. ', &
      'net electrical output.', &
      'Null figure of merit. ', &
      'max Q, max t_burn.    '  &
      /)

    ncalls = 0
    neqns = 0
    nfev1 = 0
    nfev2 = 0
    nineqns = 0
    nvar = 16
    nviter = 0
    icc = 0
    active_constraints = .false.

    lablcc = (/ &
      'Beta consistency                 ', &
      'Global power balance consistency ', &
      'Ion power balance                ', &
      'Electron power balance           ', &
      'Density upper limit              ', &
      '(Epsilon x beta-pol) upper limit ', &
      'Beam ion density consistency     ', &
      'Neutron wall load upper limit    ', &
      'Fusion power upper limit         ', &
      'Toroidal field 1/R consistency   ', &
      'Radial build consistency         ', &
      'Volt second lower limit          ', &
      'Burn time lower limit            ', &
      'NBI decay lengths consistency    ', &
      'L-H power threshold limit        ', &
      'Net electric power lower limit   ', &
      'Radiation fraction upper limit   ', &
      'Divertor heat load upper limit   ', &
      'MVA upper limit                  ', &
      'Beam tangency radius upper limit ', &
      'Plasma minor radius lower limit  ', &
      'Divertor collisionality upper lim', &
      'Conducting shell radius upper lim', &
      'Beta upper limit                 ', &
      'Peak toroidal field upper limit  ', &
      'CS coil EOF current density limit', &
      'CS coil BOP current density limit', &
      'Fusion gain Q lower limit        ', &
      'Inboard radial build consistency ', &
      'Injection power upper limit      ', &
      'TF coil case stress upper limit  ', &
      'TF coil conduit stress upper lim ', &
      'I_op / I_critical (TF coil)      ', &
      'Dump voltage upper limit         ', &
      'J_winding pack/J_protection limit', &
      'TF coil temp. margin lower limit ', &
      'Current drive gamma limit        ', &
      '1st wall coolant temp rise limit ', &
      'First wall peak temperature limit', &
      'Start-up inj. power lower limit  ', &
      'Plasma curr. ramp time lower lim ', &
      'Cycle time lower limit           ', &
      'Average centrepost temperature   ', &
      'Peak centrepost temp. upper limit', &
      'Edge safety factor lower limit   ', &
      'Ip/Irod upper limit              ', &
      'TF coil tor. thickness upper lim ', &
      'Poloidal beta upper limit        ', &
      'RFP reversal parameter < 0       ', &
      'IFE repetition rate upper limit  ', &
      'Startup volt-seconds consistency ', &
      'Tritium breeding ratio lower lim ', &
      'Neutron fluence on TF coil limit ', &
      'Peak TF coil nucl. heating limit ', &
      'Vessel helium concentration limit', &
      'Psep / R upper limit             ', &
      'TF coil leg rad width lower limit', &
      'TF coil leg rad width lower limit', &
      'NB shine-through frac upper limit', &
      'CS temperature margin lower limit', &
      'Minimum availability value       ', &
      'taup/taueff                      ', &
      'number of ITER-like vacuum pumps ', &
      'Zeff limit                       ', &
      'Dump time set by VV stress       ', &
      'Rate of change of energy in field', &
      'Upper Lim. on Radiation Wall load', &
      'Upper Lim. on Psep * Bt / q A R  ', &
      'pdivt < psep_kallenbach divertor ', &
      'Separatrix temp consistency      ', &
      'Separatrix density consistency   ', &
      'CS Tresca yield criterion        ', &
      'Psep >= Plh + Paux               ', &
      'TFC quench < tmax_croco          ', &
      'TFC current/copper area < Max    ', &
      'Eich critical separatrix density ', &
      'TFC current per turn upper limit ', &
      'Reinke criterion fZ lower limit  ', &
      'Peak CS field upper limit        ', &
      'pdivt lower limit                ', &
      'ne0 > neped                      ', &
      'toroidalgap >  tftort            ', &
      'available_space > required_space ', &
      'beta > betalim_lower             ', &
      'CP lifetime                      ', &
      'TFC turn dimension               ', &
      'Cryogenic plant power            ', &
      'TF coil strain absolute value    ', &
      'CS current/copper area < Max     ', &
      'CS stress load cycles            ', &
      'ECRH ignitability                '  &
      /)

    ! Please note: All strings between '...' above must be exactly 33 chars long
    ! Each line of code has a comma before the ampersand, except the last one.
    ! The last ad_varc line ends with the html tag "</UL>".

    ! Issue #495.  Remove default iteration variables
    ixc = 0

    ! WARNING These labels are used as variable names by write_new_in_dat.py, and possibly
    ! other python utilities, so they cannot easily be changed.
    lablxc = ''
    ! Issue 287 iteration variables are now defined in module define_iteration_variables in iteration variables.f90
    sqsumsq = 0.0D0
    norm_objf = 0.0D0
    epsfcn = 1.0D-3
    epsvmc = 1.0D-6
    factor = 0.1D0
    ftol = 1.0D-4

    boundl = 9.d-99
    ! Issue #287 These bounds now defined in initial.f90
    boundu = 9.d99

    bondl = 0.0D0
    bondu = 0.0D0
    rcm = 0.0D0
    resdl = 0.0D0
    scafc = 0.0D0
    scale = 0.0D0
    xcm = 0.0D0
    xcs = 0.0D0
    vlam = 0.0D0
    if (allocated(name_xc)) deallocate(name_xc)
    allocate(name_xc(1))
    name_xc = ""
  end subroutine init_numerics

! eqsolv() has been temporarily commented out. Please see the comment in
! function_evaluator.fcnhyb() for an explanation.

  ! subroutine eqsolv(fcnhyb,n,x,fvec,tol,epsfcn,factor,nprint,info, &
  !      wa,lwa,resdl,nfev)

  !   !! Find the non-optimising HYBRD solution to the problem
  !   !! author: Argonne National Laboratory. Minpack Project. March 1980.
  !   !! author: Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More
  !   !! author: P J Knight, CCFE, Culham Science Centre
  !   !! fcnhyb : external routine name : see below
  !   !! n : input integer : number of functions and variables
  !   !! x(n) : input/output real array : On input X must contain
  !   !! an initial estimate of the solution vector. On output X
  !   !! contains the final estimate of the solution vector.
  !   !! fvec(n) : output real array : Functions evaluated at output X
  !   !! tol : input real : Termination occurs when the algorithm
  !   !! estimates that the relative error between X and the solution
  !   !! is at most TOL.
  !   !! epsfcn : input real : Used in determining a suitable
  !   !! step length for the forward-difference approximation
  !   !! (see <A HREF="hybrd.html">hybrd</A>)
  !   !! factor : input real : Used in determining the initial step bound
  !   !! (see <A HREF="hybrd.html">hybrd</A>)
  !   !! nprint : input integer : Number of iterations between print-outs
  !   !! info : output integer : If the user has terminated execution,
  !   !! INFO is set to the (negative) value of IFLAG, see description below.
  !   !! Otherwise, INFO is set as follows:
  !   !! <PRE>
  !   !! INFO = 0   Improper input parameters.
  !   !! INFO = 1   Algorithm estimates that the relative error
  !   !! between X and the solution is at most TOL.
  !   !! INFO = 2   Number of calls to FCNHYB has reached or exceeded
  !   !! 200*(N+1).
  !   !! INFO = 3   TOL is too small. No further improvement in
  !   !! the approximate solution X is possible.
  !   !! INFO = 4   Iteration is not making good progress.
  !   !! </PRE>
  !   !! wa(lwa) : input/output real array : work array
  !   !! lwa : input integer : work array size, not less than (N*(3*N+13))/2
  !   !! resdl(n) : output real array : residuals
  !   !! nfev : output integer : number of iterations performed
  !   !! Routine EQSOLV is the Argonne Minpack subroutine HYBRD1
  !   !! which has been modified by D.T. Blackfield FEDC/TRW.
  !   !! The routine is the same except some of the arguments are
  !   !! user supplied rather than 'hardwired'.
  !   !! <P>The purpose of EQSOLV is to find a zero of a system of
  !   !! N nonlinear functions in N variables by a modification
  !   !! of the Powell hybrid method. This is done by using the
  !   !! more general nonlinear equation solver <A HREF="hybrd.html">HYBRD</A>.
  !   !! The user must provide a subroutine which calculates the functions.
  !   !! The Jacobian is then calculated by a forward-difference
  !   !! approximation.
  !   !! <P>FCNHYB is the name of a user-supplied subroutine which
  !   !! calculates the functions. FCNHYB must be declared
  !   !! in an external statement in the user calling
  !   !! program, and should be written as follows:
  !   !! <PRE>
  !   !! subroutine fcnhyb(n,x,fvec,iflag)
  !   !! integer n,iflag
  !   !! double precision x(n),fvec(n)
  !   !! ----------
  !   !! calculate the functions at x and
  !   !! return this vector in fvec.
  !   !! ---------
  !   !! return
  !   !! end
  !   !! </PRE>
  !   !! The value of iflag should not be changed by FCNHYB unless
  !   !! the user wants to terminate execution of EQSOLV.
  !   !! In this case set IFLAG to a negative integer.
  !   !! None
  !   !
  !   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! 	use global_variables, only: icase, verbose, vlabel
  !   use constants, only: mfile, nout
  !   use maths_library, only: HYBRD
  !   implicit none

  !   ! Interface for the external fcnhyb subroutine argument
  !   ! This interface is necessary when wrapping with f2py
  !   interface
  !     subroutine fcnhyb(n, x, fvec, iflag)
  !       use, intrinsic :: iso_fortran_env, only: dp=>real64
  !       integer, intent(in) :: n
  !       real(dp), dimension(n), intent(inout) :: x
  !       real(dp), dimension(n), intent(out) :: fvec
  !       integer, intent(inout) :: iflag
  !     end subroutine fcnhyb
  !   end interface

  !   !  Arguments

  !   external :: fcnhyb
  !   integer, intent(in) :: n, nprint, lwa
  !   real(dp), dimension(n), intent(inout) :: x
  !   real(dp), dimension(n), intent(out) :: fvec, resdl
  !   real(dp), dimension(lwa), intent(out) :: wa
  !   real(dp), intent(in) :: tol, epsfcn, factor
  !   integer, intent(out) :: info, nfev

  !   !  Local variables

  !   integer :: n1,indx,lr,maxfev,ml,mode,mu
  !   real(dp), parameter :: one = 1.0D0
  !   real(dp), parameter :: zero = 0.0D0
  !   real(dp) :: xtol

  !   ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !   info = 0

  !   !  Check the input parameters for errors

  !   if ( (n == 0).or.(tol < zero).or.(lwa < ((n*(3*n + 13))/2) ) ) return

  !   !  Call HYBRD

  !   maxfev = 200*(n + 1)
  !   xtol = tol
  !   ml = n - 1
  !   mu = n - 1
  !   mode = 2

  !   wa(:) = one

  !   lr = (n*(n + 1))/2
  !   indx = 6*n + lr
  !   n1 = n

  !   !+**PJK 23/10/92 Warning produced by QA Fortran :
  !   !+**PJK 23/10/92 Arg 16 in call to HYBRD has wrong dimensions.
  !   !+**PJK 23/10/92 Code works at present, but beware of future
  !   !+**PJK 23/10/92 modifications.

  !   call hybrd(fcnhyb,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,wa(1),mode, &
  !        factor,nprint,info,nfev,wa(indx+1),n1,wa(6*n+1),lr, &
  !        wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1), &
  !        resdl)

  !   if (info == 5) info = 4

  ! end subroutine eqsolv
end module numerics
