! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constraints
  !! author: J Morris
  !!
  !! Module defining the constraint equations and the routine that evaluates the
  !! constraint equations.

  ! Import modules
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  use error_handling, only: report_error, idiags, fdiags

  implicit none

  public :: constraint_eqns

!   type constraint_args_type
!     real(dp) :: cc
!     !! Residual error in constraint equation
!     real(dp) :: con
!     !! constraint value for constraint equation in physical units
!     real(dp) :: err
!     !! residual error in constraint equation in physical units
!     character(len=1)  :: symbol
!     !! `=<`, `>`, `<` symbol for constraint equation denoting its type
!     character(len=10) :: units
!     !! constraint units in constraint equation
!   end type

contains

   subroutine constraint_eqns(m,ieqn,cc,con,err,symbol,units)
    !! Routine that formulates the constraint equations
    !!
    !! **author: P J Knight** (UKAEA)
    !!
    !! **author: J Morris** (UKAEA)
    !!
    !! if `ieqn` is zero or negative, evaluate all the constraint equations, otherwise
    !! evaluate only the `ieqn`th equation. The code attempts to make \( cc(i) = 0 \) for all
    !! \( i \in \{1,\cdots,m\} \) equations. All relevant consistency equations should be active in
    !! order to make a self-consistent machine.
    !!
    !! **References**
    !!
    !! 1. AEA FUS 251: A User's Guide to the PROCESS Systems Code

    use numerics, only: icc
    use maths_library, only: variable_error

    implicit none

    integer, intent(in) :: m
    !! Number of constraint equations to solve

    integer, intent(in) :: ieqn
    !! Switch for constraint equations to evaluate;

    real(dp), dimension(m), intent(out) :: cc
    !! Residual error in equation i

    real(dp), optional, dimension(m), intent(out) :: con
    !! constraint value for equation i in physical units

    real(dp), optional, dimension(m), intent(out) :: err
    !! residual error in equation i in physical units

    character(len=1),  optional, dimension(m), intent(out) :: symbol
    !! `=<`, `>`, `<` symbol for equation i denoting its type

    character*10, optional, dimension(m), intent(out) :: units
    !! constraint units in equation i

    ! Local variables
    integer :: i, i1, i2

    real(dp) :: tmp_cc = 0
    !! Residual error in constraint equation
    real(dp) :: tmp_con = 0
    !! constraint value for constraint equation in physical units
    real(dp) :: tmp_err = 0
    !! residual error in constraint equation in physical units
    character(len=1)  :: tmp_symbol = ' '
    !! `=<`, `>`, `<` symbol for constraint equation denoting its type
    character(len=10) :: tmp_units = ' '
    !! constraint units in constraint equation


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! If ieqn is positive, only evaluate the 'ieqn'th constraint residue,
    ! otherwise evaluate all m constraint residues
    if (ieqn > 0) then
      i1 = ieqn ; i2 = ieqn
    else
      i1 = 1 ; i2 = m
    end if

    ! Consistency (equality) constraints should converge to zero.
    do i = i1,i2

      ! The constraint value in physical units is
      ! a) for consistency equations, the quantity to be equated, or
      ! b) for limit equations, the limiting value.
      ! The symbol is = for a consistency equation, < for an upper limit
      ! or > for a lower limit.
      select case (icc(i))

	      ! Relationship between beta, temperature (keV) and density
        case (1); call constraint_eqn_001(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Global plasma power balance equation
        case (2); call constraint_eqn_002(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Global power balance equation for ions
        case (3); call constraint_eqn_003(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Global power balance equation for electrons
        case (4); call constraint_eqn_004(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for density upper limit
        case (5); call constraint_eqn_005(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for epsilon beta-poloidal upper limit
        case (6); call constraint_eqn_006(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for hot beam ion density
        case (7); call constraint_eqn_007(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for neutron wall load upper limit
        case (8); call constraint_eqn_008(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for fusion power upper limit
        case (9); call constraint_eqn_009(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Obsolete
        case (10); call constraint_eqn_010(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for radial build
        case (11); call constraint_eqn_011(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for volt-second capability lower limit
        case (12); call constraint_eqn_012(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for burn time lower limit
        case (13); call constraint_eqn_013(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation to fix number of NBI decay lengths to plasma centre
        case (14); call constraint_eqn_014(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for L-H power threshold limit
        case (15); call constraint_eqn_015(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for net electric power lower limit
        case (16); call constraint_eqn_016(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for radiation power upper limit
        case (17); call constraint_eqn_017(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for divertor heat load upper limit
        case (18); call constraint_eqn_018(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for MVA upper limit
        case (19); call constraint_eqn_019(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for neutral beam tangency radius upper limit
        case (20); call constraint_eqn_020(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for minor radius lower limit
        case (21); call constraint_eqn_021(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for divertor collision/connection length ratio upper limit
        case (22); call constraint_eqn_022(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for conducting shell radius / rminor upper limit
        case (23); call constraint_eqn_023(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for beta upper limit
        case (24); call constraint_eqn_024(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for peak toroidal field upper limit
        case (25); call constraint_eqn_025(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for Central Solenoid current density upper limit at EOF
        case (26); call constraint_eqn_026(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for Central Solenoid current density upper limit at BOP
        case (27); call constraint_eqn_027(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for fusion gain (big Q) lower limit
        case (28); call constraint_eqn_028(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for inboard major radius
        case (29); call constraint_eqn_029(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for injection power upper limit
        case (30); call constraint_eqn_030(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for TF coil case stress upper limit (SCTF)
        case (31); call constraint_eqn_031(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for TF coil conduit stress upper limit (SCTF)
        case (32); call constraint_eqn_032(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for TF coil operating/critical J upper limit (SCTF)
        case (33); call constraint_eqn_033(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for TF coil dump voltage upper limit (SCTF)
        case (34); call constraint_eqn_034(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for TF coil J_wp/J_prot upper limit (SCTF)
        case (35); call constraint_eqn_035(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for TF coil s/c temperature margin lower limit (SCTF)
        case (36); call constraint_eqn_036(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for current drive gamma upper limit
        case (37); call constraint_eqn_037(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for first wall temperature upper limit
        case (39); call constraint_eqn_039(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for auxiliary power lower limit
        case (40); call constraint_eqn_040(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for plasma current ramp-up time lower limit
        case (41); call constraint_eqn_041(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for cycle time lower limit
        case (42); call constraint_eqn_042(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for average centrepost temperature
        case (43); call constraint_eqn_043(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for centrepost temperature upper limit (TART)
        case (44); call constraint_eqn_044(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for edge safety factor lower limit (TART)
        case (45); call constraint_eqn_045(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for Ip/Irod upper limit (TART)
        case (46); call constraint_eqn_046(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for TF coil toroidal thickness upper limit
        case (47); call constraint_eqn_047(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for poloidal beta upper limit
        case (48); call constraint_eqn_048(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Issue #508 Remove RFP option Equation to ensure reversal parameter F is negative
        case (49); call constraint_eqn_049(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! IFE option: Equation for repetition rate upper limit
        case (50); call constraint_eqn_050(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation to enforce startup flux = available startup flux
        case (51); call constraint_eqn_051(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for tritium breeding ratio lower limit
        case (52); call constraint_eqn_052(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for fast neutron fluence on TF coil upper limit
        case (53); call constraint_eqn_053(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for peak TF coil nuclear heating upper limit
        case (54); call constraint_eqn_054(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for helium concentration in vacuum vessel upper limit
        case (55); call constraint_eqn_055(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for power through separatrix / major radius upper limit
        case (56); call constraint_eqn_056(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Obsolete
        case (57); call constraint_eqn_057(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Obsolete
        case (58); call constraint_eqn_058(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for neutral beam shine-through fraction upper limit
        case (59); call constraint_eqn_059(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for Central Solenoid s/c temperature margin lower limit
        case (60); call constraint_eqn_060(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Equation for availability limit
        case (61); call constraint_eqn_061(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Lower limit on taup/taueff the ratio of alpha particle to energy confinement times
        case (62); call constraint_eqn_062(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Upper limit on niterpump (vacuum_model = simple)
        case (63); call constraint_eqn_063(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Upper limit on Zeff
        case (64); call constraint_eqn_064(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Limit TF dump time to calculated quench time
        case (65); call constraint_eqn_065(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Limit on rate of change of energy in poloidal field
        case (66); call constraint_eqn_066(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Simple upper limit on radiation wall load
        case (67); call constraint_eqn_067(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! New Psep scaling (PsepB/qAR)
        case (68); call constraint_eqn_068(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Ensure separatrix power is less than value from Kallenbach divertor
        case (69); call constraint_eqn_069(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Separatrix temperature consistency
        case (70); call constraint_eqn_070(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Separatrix density consistency
        case (71); call constraint_eqn_071(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	      ! Central Solenoid Tresca yield criterion
        case (72); call constraint_eqn_072(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! ensure separatrix power is greater than the L-H power + auxiliary power
        case (73); call constraint_eqn_073(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! ensure TF coil quench temperature < tmax_croco
        case (74); call constraint_eqn_074(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	     ! ensure that TF coil current / copper area < Maximum value ONLY used for croco HTS coil
        case (75); call constraint_eqn_075(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Eich critical separatrix density model
        case (76); call constraint_eqn_076(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Equation for maximum TF current per turn upper limit
        case (77); call constraint_eqn_077(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
	  	  ! Equation for Reinke criterion, divertor impurity fraction lower limit
        case (78); call constraint_eqn_078(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Equation for maximum CS field
        case (79); call constraint_eqn_079(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Lower limit pdivt
        case (80); call constraint_eqn_080(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Constraint equation making sure that ne(0) > ne(ped)
        case (81); call constraint_eqn_081(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Constraint equation making sure that stellarator coils dont touch in toroidal direction
        case (82); call constraint_eqn_082(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Constraint ensuring radial build consistency for stellarators
        case (83); call constraint_eqn_083(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
        ! Constraint for lower limit of beta
        case (84); call constraint_eqn_084(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
         ! Constraint for CP lifetime
        case (85); call constraint_eqn_085(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
         ! Constraint for turn dimension
        case (86); call constraint_eqn_086(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
         ! Constraint for cryogenic power
        case (87); call constraint_eqn_087(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
         ! Constraint for TF coil strain
        case (88); call constraint_eqn_088(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
         ! ensure that OH coil current / copper area < Maximum value ONLY used for croco HTS coil
        case (89); call constraint_eqn_089(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
         ! Constraint for minimum CS stress load cycles
        case (90); call constraint_eqn_090(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
         ! Constraint for indication of ECRH ignitability
        case (91); call constraint_eqn_091(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
       case default

        idiags(1) = icc(i)
        call report_error(13)
        tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ' '
        tmp_units = ' '

      end select

      cc(i) = tmp_cc
      if (present(con)) then
        con(i) = tmp_con
        if (present(err))    err(i)    = tmp_err
        if (present(symbol)) symbol(i) = tmp_symbol
        if (present(units))  units(i)  = tmp_units
      end if

      ! Crude method of catching NaN errors
      !if ((abs(cc(i)) > 9.99D99).or.(cc(i) /= cc(i))) then
      if (variable_error(cc(i))) then

        ! Add debugging lines as appropriate...
        select case (icc(i))

          ! Relationship between beta, temperature (keV) and density (consistency equation)
          case (1); call constraint_err_001()
          ! Equation for net electric power lower limit
          case (16); call constraint_err_016()
          ! Equation for injection power upper limit
          case (30); call constraint_err_030()
          ! Limit on rate of change of energy in poloidal field
          case (66); call constraint_err_066()

        end select

        idiags(1) = icc(i) ; fdiags(1) = cc(i)
        call report_error(14)

      end if

    end do
    ! Issue 505 Reverse the sign so it works as an inequality constraint (cc(i) > 0)
    ! This will have no effect if it is used as an equality constraint because it will be squared.
    cc = -cc

   end subroutine constraint_eqns

   !--- Error-handling routines

   subroutine constraint_err_001()
    !! Error in: Relationship between beta, temperature (keV) and density (consistency equation)
    !! author: P B Lloyd, CCFE, Culham Science Centre
    use physics_variables, only: betaft, betanb, dene, ten, dnitot, tin, btot, beta
    write(*,*) 'betaft = ', betaft
    write(*,*) 'betanb = ', betanb
    write(*,*) 'dene = ', dene
    write(*,*) 'ten = ', ten
    write(*,*) 'dnitot = ', dnitot
    write(*,*) 'tin = ', tin
    write(*,*) 'btot = ',btot
    write(*,*) 'beta = ', beta
   end subroutine

   subroutine constraint_err_016()
    !! Error in: Equation for net electric power lower limit
    !! author: P B Lloyd, CCFE, Culham Science Centre
    use constraint_variables, only: fpnetel, pnetelin
    use heat_transport_variables, only: pnetelmw
    implicit none
    write(*,*) 'fpnetel = ', fpnetel
    write(*,*) 'pnetelmw = ', pnetelmw
    write(*,*) 'pnetelin = ', pnetelin
   end subroutine

   subroutine constraint_err_030()
    !! Error in: Equation for injection power upper limit
    !! author: P B Lloyd, CCFE, Culham Science Centre
    use current_drive_variables, only: pinjmw, pinjalw
    use constraint_variables, only: fpinj
    implicit none
    write(*,*) 'fpinj = ', fpinj
    write(*,*) 'pinjalw = ', pinjalw
    write(*,*) 'pinjmw = ', pinjmw
   end subroutine

   subroutine constraint_err_066()
    !! Error in: Limit on rate of change of energy in poloidal field
    !! author: P B Lloyd, CCFE, Culham Science Centre
    use constraint_variables, only: fpoloidalpower
    use pf_power_variables, only: maxpoloidalpower, peakpoloidalpower
    implicit none
    write(*,*) 'fpoloidalpower = ', fpoloidalpower
    write(*,*) 'maxpoloidalpower = ', maxpoloidalpower
    write(*,*) 'peakpoloidalpower = ', peakpoloidalpower
   end subroutine constraint_err_066

   !---

   subroutine constraint_eqn_001(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
    !! author: J Morris
    !! category: equality constraint
    !!
    !! Relationship between beta, temperature (keV) and density
    !!
    !! \begin{equation}
    !! c_i = 1 - \frac{1}{\beta}\left( \beta_{ft} + \beta_{NBI} + 2 \times 10^3 \mu_0 e
    !! \left( \frac{n_e T_e + n_i T_i}{B_{tot}^2} \right) \right)
    !! \end{equation}
    !!
    !! - \( \beta \) -- total plasma beta
    !! - \( \beta_{ft} \) -- fast alpha beta component
    !! - \( \beta_{NBI} \) -- neutral beam beta component
    !! - \( n_e \) -- electron density [m\(^3\)]
    !! - \( n_i \) -- total ion density [m\(^3\)]
    !! - \( T_e \) -- density weighted average electron temperature [keV]
    !! - \( T_i \) -- density weighted average ion temperature [keV]
    !! - \( B_{tot} \) -- total toroidal + poloidal field [T]

    use physics_variables, only: betaft, betanb, dene, ten, dnitot, tin, btot, beta
    use constants, only: echarge,rmu0

    implicit none

   !  type(constraint_args_type), intent(out) :: args
      real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

    !! constraint derived type

      tmp_cc = 1.0D0 - (betaft + betanb + &
        2.0D3*rmu0*echarge * (dene*ten + dnitot*tin)/btot**2 )/beta
      tmp_con = beta * (1.0D0 - tmp_cc)
      tmp_err = beta * tmp_cc
      tmp_symbol = '='
      tmp_units  = ''

   end subroutine constraint_eqn_001

   subroutine constraint_eqn_002(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
    !! author: J. Morris
    !! category: equality constraint
    !!
    !! Global plasma power balance equation
    !!
    !! \begin{equation} c_i =
    !! \end{equation}
    !!
    !! iradloss : input integer : switch for radiation loss term usage in power balance (see User Guide):<UL>
    !! <LI> = 0 total power lost is scaling power plus radiation (needed for ipedestal=2,3)
    !! <LI> = 1 total power lost is scaling power plus core radiation only
    !! <LI> = 2 total power lost is scaling power only, with no additional
    !! allowance for radiation. This is not recommended for power plant models.</UL>
    !! ignite : input integer : switch for ignition assumption:<UL>
    !! <LI> = 0 do not assume plasma ignition;
    !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
    !! ptrepv : input real : electron transport power per volume (MW/m3)
    !! ptripv : input real :  ion transport power per volume (MW/m3)
    !! pradpv : input real : total radiation power per volume (MW/m3)
    !! pcoreradpv : input real : total core radiation power per volume (MW/m3)
    !! falpha : input real : fraction of alpha power deposited in plasma
    !! palppv : input real : alpha power per volume (MW/m3)
    !! pchargepv : input real : non-alpha charged particle fusion power per volume (MW/m3)
    !! pohmpv : input real : ohmic heating power per volume (MW/m3)
    !! pinjmw : input real : total auxiliary injected power (MW)
    !! vol : input real : plasma volume (m3)

    use physics_variables, only: iradloss, ignite, ptrepv, ptripv, pradpv, &
                                  pcoreradpv, falpha, palppv, pchargepv, &
                                  pohmpv, vol
    use current_drive_variables, only: pinjmw

    implicit none

          real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
    !! constraint derived type

    ! pscaling : Local real : total transport power per volume (MW/m3)
    real(dp) :: pscaling
    real(dp) :: pnumerator, pdenom
    pscaling = ptrepv + ptripv
    ! Total power lost is scaling power plus radiation:
    if (iradloss == 0) then
        pnumerator = pscaling + pradpv
    else if (iradloss == 1) then
        pnumerator = pscaling + pcoreradpv
    else
        pnumerator = pscaling
    end if

    ! if plasma not ignited include injected power
    if (ignite == 0) then
      pdenom = falpha*palppv + pchargepv + pohmpv + pinjmw/vol
    else
      ! if plasma ignited
      pdenom = falpha*palppv + pchargepv + pohmpv
    end if

    tmp_cc = 1.0D0 - pnumerator / pdenom
    tmp_con = pdenom * (1.0D0 - tmp_cc)
    tmp_err = pdenom * tmp_cc
    tmp_symbol = '='
    tmp_units = 'MW/m3'

   end subroutine constraint_eqn_002

   subroutine constraint_eqn_003(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Global power balance equation for ions
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Global power balance equation for ions
      !! This is a consistency equation (NBI)
      !! #=# physics
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ignite : input integer : switch for ignition assumption:<UL>
      !! <LI> = 0 do not assume plasma ignition;
      !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
      !! ptripv : input real :  ion transport power per volume (MW/m3)
      !! piepv : input real : ion/electron equilibration power per volume (MW/m3)
      !! falpha : input real : fraction of alpha power deposited in plasma
      !! palpipv : input real : alpha power per volume to ions (MW/m3)
      !! pinjimw : input real : auxiliary injected power to ions (MW)
      !! vol : input real : plasma volume (m3)
      use physics_variables, only: ignite, ptripv, piepv, falpha, palpipv, vol
      use current_drive_variables, only: pinjimw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

	   ! No assume plasma ignition:
      if (ignite == 0) then
         tmp_cc     = 1.0D0 - (ptripv + piepv) / (falpha*palpipv + pinjimw/vol)
         tmp_con    = (falpha*palpipv + pinjimw/vol) * (1.0D0 - tmp_cc)
         tmp_err    = (falpha*palpipv + pinjimw/vol) * tmp_cc
         tmp_symbol = '='
         tmp_units  = 'MW/m3'
	   ! Plasma ignited:
      else
         tmp_cc     = 1.0D0 - (ptripv+piepv) / (falpha*palpipv)
         tmp_con    = (falpha*palpipv) * (1.0D0 - tmp_cc)
         tmp_err    = (falpha*palpipv) * tmp_cc
         tmp_symbol = '='
         tmp_units  = 'MW/m3'
      end if

   end subroutine constraint_eqn_003

   subroutine constraint_eqn_004(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Global power balance equation for electrons
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Global power balance equation for electrons
      !! This is a consistency equation
      !! N.B. This constraint is currently NOT RECOMMENDED for use.
      !! #=# physics
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! iradloss : input integer : switch for radiation loss term usage in power balance (see User Guide):<UL>
      !! <LI> = 0 total power lost is scaling power plus radiation (needed for ipedestal=2,3)
      !! <LI> = 1 total power lost is scaling power plus core radiation only
      !! <LI> = 2 total power lost is scaling power only, with no additional
      !! allowance for radiation. This is not recommended for power plant models.</UL>
      !! ignite : input integer : switch for ignition assumption:<UL>
      !! <LI> = 0 do not assume plasma ignition;
      !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
      !! ptrepv : input real : electron transport power per volume (MW/m3)
      !! pradpv : input real : total radiation power per volume (MW/m3)
      !! pcoreradpv : input real : total core radiation power per volume (MW/m3)
      !! falpha : input real : fraction of alpha power deposited in plasma
      !! palpepv : input real : alpha power per volume to electrons (MW/m3)
      !! piepv : input real : ion/electron equilibration power per volume (MW/m3)
      !! pinjemw : input real : auxiliary injected power to electrons (MW)
      !! vol : input real : plasma volume (m3)
      use physics_variables, only: iradloss, ignite, ptrepv, pcoreradpv, falpha, &
                                 palpepv, piepv, vol, pradpv
      use current_drive_variables, only: pinjemw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! pscaling : Local real : total transport power per volume (MW/m3)
      real(dp) :: pscaling
      real(dp) :: pnumerator, pdenom
      pscaling = ptrepv
	   ! Total power lost is scaling power plus radiation:
      if (iradloss == 0) then
         pnumerator = pscaling + pradpv
      else if (iradloss == 1) then
         pnumerator = pscaling + pcoreradpv
      else
         pnumerator = pscaling
      end if

      ! if plasma not ignited include injected power
      if (ignite == 0) then
         pdenom = falpha*palpepv + piepv + pinjemw/vol
      else
      ! if plasma ignited
         pdenom = falpha*palpepv + piepv
      end if

      tmp_cc     = 1.0D0 - pnumerator / pdenom
      tmp_con    = pdenom * (1.0D0 - tmp_cc)
      tmp_err    = pdenom * tmp_cc
      tmp_symbol = '='
      tmp_units  = 'MW/m3'

   end subroutine constraint_eqn_004

   subroutine constraint_eqn_005(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for density upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for density upper limit
      !! #=# physics
      !! #=#=# fdene, dnelimt
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! idensl : input integer : switch for density limit to enforce (constraint equation 5):<UL>
      !! <LI> = 1 old ASDEX;
      !! <LI> = 2 Borrass model for ITER (I);
      !! <LI> = 3 Borrass model for ITER (II);
      !! <LI> = 4 JET edge radiation;
      !! <LI> = 5 JET simplified;
      !! <LI> = 6 Hugill-Murakami Mq limit;
      !! <LI> = 7 Greenwald limit</UL>
      !! fdene : input real : f-value for density limit
      !! dene : input real : electron density (/m3)
      !! dnelimt : input real : density limit (/m3)
      !! dnla : input real : line averaged electron density (m-3)
      use physics_variables, only: idensl, dnelimt, dnla, dene
      use constraint_variables, only: fdene
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

	   ! Apply Greenwald limit to line-averaged density
      if (idensl == 7) then
         tmp_cc     = 1.0D0 - fdene * dnelimt/dnla
         tmp_con    = fdene * dnelimt
         tmp_err    = fdene * dnelimt - dnla
         tmp_symbol = '<'
         tmp_units  = '/m3'
      else
         tmp_cc = 1.0D0 - fdene * dnelimt/dene
         tmp_con    = dnelimt * (1.0D0 - tmp_cc)
         tmp_err    = dene * tmp_cc
         tmp_symbol = '<'
         tmp_units  = '/m3'
      end if

   end subroutine constraint_eqn_005

   subroutine constraint_eqn_006(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for epsilon beta-poloidal upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for epsilon beta-poloidal upper limit
      !! #=# physics
      !! #=#=# fbeta, epbetmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fbeta : input real : f-value for epsilon beta-poloidal
      !! epbetmax : input real : maximum (eps*beta_poloidal)
      !! eps : input real : inverse aspect ratio
      !! betap : input real : poloidal beta
      use physics_variables, only: epbetmax, eps, betap
      use constraint_variables, only: fbeta, fbeta
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fbeta * epbetmax/(eps*betap)
      tmp_con = epbetmax * (1.0D0 - tmp_cc)
      tmp_err = (eps*betap) * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_006

   subroutine constraint_eqn_007(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for hot beam ion density
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for hot beam ion density
      !! This is a consistency equation (NBI)
      !! #=# physics
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ignite : input integer : switch for ignition assumption:<UL>
      !! <LI> = 0 do not assume plasma ignition;
      !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
      !! Obviously, ignite must be zero if current drive is required.
      !! If ignite=1, any auxiliary power is assumed to be used only
      !! during plasma start-up, and is excluded from all steady-state
      !! power balance calculations.
      !! dnbeam2 : input real :  hot beam ion density from calculation (/m3)
      !! dnbeam : input real : hot beam ion density, variable (/m3)
      use physics_variables, only: ignite, dnbeam2, dnbeam
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

	   ! Do not assume plasma ignition:
      if (ignite == 0) then
         tmp_cc     = 1.0D0 - dnbeam2/dnbeam
         tmp_con    = dnbeam * (1.0D0 - tmp_cc)
         tmp_err    = dnbeam * tmp_cc
         tmp_symbol = '='
         tmp_units  = '/m3'
      else
         tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''
         call report_error(1)
      end if

   end subroutine constraint_eqn_007

   subroutine constraint_eqn_008(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for neutron wall load upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for neutron wall load upper limit
      !! #=# physics
      !! #=#=# fwalld, walalw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fwalld : input real : f-value for maximum wall load
      !! walalw : input real : allowable wall-load (MW/m2)
      !! wallmw : input real : average neutron wall load (MW/m2)
      use constraint_variables, only: fwalld, walalw
      use physics_variables, only: wallmw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fwalld * walalw/wallmw
      tmp_con = fwalld * walalw
      tmp_err = fwalld * walalw - wallmw
      tmp_symbol = '<'
      tmp_units = 'MW/m2'

   end subroutine constraint_eqn_008

   subroutine constraint_eqn_009(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for fusion power upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for fusion power upper limit
      !! #=# physics
      !! #=#=# ffuspow, powfmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ffuspow : input real : f-value for maximum fusion power
      !! powfmax : input real : maximum fusion power (MW)
      !! powfmw : input real : fusion power (MW)
      use constraint_variables, only: ffuspow, powfmax
      use physics_variables, only: powfmw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ffuspow * powfmax/powfmw
      tmp_con = powfmax * (1.0D0 - tmp_cc)
      tmp_err = powfmw * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW'

   end subroutine constraint_eqn_009

   subroutine constraint_eqn_010(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Equation for field at TF coil
      !! This is a consistency equation
      !! (do not use for stellarators)
      !! #=# tfcoil
      !! #=#=# consistency
      !! rmajor |  plasma major radius (m)
      !! bt     |  toroidal field on axis (T)
      !! rbmax  |  radius of maximum toroidal field (m)
      !! bmaxtf |  peak field at toroidal field coil (T)

      !! This constraint is depreciated

      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      !! Constraints output

      ! This constraint is depreciated
      call report_error(236)

      tmp_con = 1.0D0
      tmp_err = 0.0D0
      tmp_symbol = '='
      tmp_units = ''

   end subroutine constraint_eqn_010

   subroutine constraint_eqn_011(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for radial build
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for radial build
      !! (This is a consistency equation.)
      !! #=# build
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! rbld : input real : sum of thicknesses to the major radius (m)
      !! rmajor : input real : plasma major radius (m)
      use build_variables, only: rbld
      use physics_variables, only: rmajor
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - rbld/rmajor
      tmp_con = rmajor * (1.0D0 - tmp_cc)
      tmp_err = rmajor * tmp_cc
      tmp_symbol = '='
      tmp_units = 'm'

   end subroutine constraint_eqn_011

   subroutine constraint_eqn_012(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for volt-second capability lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for volt-second capability lower limit
      !! #=# pfcoil
      !! #=#=# fvs, vsstt
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! vsstt : input real : total V-s needed (Wb)
      !! vsstt (lower limit) is positive; vstot (available) is negative
      !! fvs : input real : f-value for flux-swing (V-s) requirement (STEADY STATE)
      !! vstot : input real :   total flux swing for pulse (Wb)
      use physics_variables, only: vsstt
      use constraint_variables, only: fvs
      use pfcoil_variables, only: vstot
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 + fvs * vstot/vsstt
      tmp_con = vsstt * (1.0D0 - tmp_cc)
      tmp_err = vsstt * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'V.sec'

   end subroutine constraint_eqn_012

   subroutine constraint_eqn_013(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for burn time lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for burn time lower limit
      !! #=# times
      !! #=#=# ftburn, tbrnmn
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftburn : input real : f-value for minimum burn time
      !! tburn : input real : burn time (s) (calculated if lpulse=1)
      !! tbrnmn : input real :  minimum burn time (s)
      use constraint_variables, only: ftburn,tbrnmn
      use times_variables, only: tburn
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ftburn * tburn/tbrnmn
      tmp_con = tbrnmn / ftburn
      tmp_err = tbrnmn / ftburn  - tburn
      tmp_symbol = '>'
      tmp_units = 'sec'

   end subroutine constraint_eqn_013

   subroutine constraint_eqn_014(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation to fix number of NBI decay lengths to plasma centre
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation to fix number of NBI decay lengths to plasma centre
      !! This is a consistency equation
      !! #=# current_drive
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! taubeam : input real : neutral beam e-decay lengths to plasma centre
      !! tbeamin : input real : permitted neutral beam e-decay lengths to plasma centre
      use current_drive_variables, only: taubeam, tbeamin
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - taubeam/tbeamin
      tmp_con = tbeamin * (1.0D0 - tmp_cc)
      tmp_err = tbeamin * tmp_cc
      tmp_symbol = '='
      tmp_units = ''

   end subroutine constraint_eqn_014

   subroutine constraint_eqn_015(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for L-H power threshold limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for L-H power threshold limit
      !! #=# physics
      !! #=#=# flhthresh, plhthresh
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! flhthresh : input real : f-value for L-H power threshold
      !! plhthresh : input real : L-H mode power threshold (MW)
      !! pdivt : input real : power to conducted to the divertor region (MW)
      use constraint_variables, only: flhthresh
      use physics_variables, only: plhthresh, pdivt
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  -(1.0D0 - flhthresh * plhthresh / pdivt)
      tmp_con = plhthresh
      tmp_err = plhthresh - pdivt / flhthresh
      if (flhthresh > 1.0D0) then
         tmp_symbol = '>'
      else
         tmp_symbol = '<'
      end if
      tmp_units = 'MW'

   end subroutine constraint_eqn_015

   subroutine constraint_eqn_016(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for net electric power lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for net electric power lower limit
      !! #=# heat_transport
      !! #=#=# fpnetel, pnetelin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpnetel : input real : f-value for net electric power
      !! pnetelmw : input real : net electric power (MW)
      !! pnetelin : input real : required net electric power (MW)
      use constraint_variables, only: fpnetel, pnetelin
      use heat_transport_variables, only: pnetelmw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fpnetel * pnetelmw / pnetelin
      tmp_con = pnetelin
      tmp_err = pnetelmw - pnetelin / fpnetel
      tmp_symbol = '>'
      tmp_units = 'MW'

   end subroutine constraint_eqn_016

   subroutine constraint_eqn_017(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for radiation power upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for radiation power upper limit
      !! #=# physics
      !! #=#=# fradpwr, pradmaxpv
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! falpha : input real : fraction of alpha power deposited in plasma
      !! pinjmw : input real : total auxiliary injected power (MW)
      !! vol : input real : plasma volume (m3)
      !! palppv : input real : alpha power per volume (MW/m3)
      !! pchargepv :  input real : non-alpha charged particle fusion power per volume (MW/m3)
      !! pohmpv : input real : ohmic heating power per volume (MW/m3)
      !! fradpwr : input real : f-value for core radiation power limit
      !! pradpv : input real : total radiation power per volume (MW/m3)
      use physics_variables, only: falpha, vol, palppv, pchargepv, pohmpv, pradpv
      use current_drive_variables, only: pinjmw
      use constraint_variables, only: fradpwr
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      real(dp) :: pradmaxpv
      !! Maximum possible power/vol that can be radiated (local)

      pradmaxpv = pinjmw/vol + palppv*falpha + pchargepv + pohmpv
      tmp_cc =  1.0D0 - fradpwr * pradmaxpv / pradpv
      tmp_con = pradmaxpv * (1.0D0 - tmp_cc)
      tmp_err = pradpv * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m3'

   end subroutine constraint_eqn_017

   subroutine constraint_eqn_018(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for divertor heat load upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for divertor heat load upper limit
      !! #=# divertor
      !! #=#=# fhldiv, hldivlim
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fhldiv : input real : peak resistive TF coil inboard leg power (MW)
      !! hldivlim : input real : heat load limit (MW/m2)
      !! hldiv : input real : divertor heat load (MW/m2)
      use constraint_variables, only: fhldiv
      use divertor_variables, only: hldivlim, hldiv
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fhldiv * hldivlim/hldiv
      tmp_con = hldivlim * (1.0D0 - tmp_cc)
      tmp_err = hldiv * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m2'

   end subroutine constraint_eqn_018

   subroutine constraint_eqn_019(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for MVA upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for MVA upper limit
      !! #=# tfcoil
      !! #=#=# fmva, mvalim
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! tfcpmw : input real : peak resistive TF coil inboard leg power (MW)
      !! tflegmw : input real : TF coil outboard leg resistive power (MW)
      !! fmva : input real : f-value for maximum MVA
      !! mvalim : input real : TF coil outboard leg resistive power (MW)
      use tfcoil_variables, only: tfcpmw, tflegmw
      use constraint_variables, only: fmva, mvalim
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      ! totmva : local real : total MVA in TF coil (MW)
      real(dp) :: totmva

      totmva = tfcpmw + tflegmw
      tmp_cc =  1.0D0 - fmva * mvalim/totmva
      tmp_con = mvalim * (1.0D0 - tmp_cc)
      tmp_err = totmva * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MVA'

   end subroutine constraint_eqn_019

   subroutine constraint_eqn_020(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for neutral beam tangency radius upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for neutral beam tangency radius upper limit
      !! #=# current_drive
      !! #=#=# fportsz, rtanmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fportsz : input real : f-value for neutral beam tangency radius limit
      !! rtanmax : input real : maximum tangency radius for centreline of beam (m)
      !! rtanbeam : input real : ratio of collision length / connection length
      use constraint_variables, only: fportsz
      use current_drive_variables, only: rtanmax, rtanbeam
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fportsz * rtanmax/rtanbeam
      tmp_con = rtanmax * (1.0D0 - tmp_cc)
      tmp_err = rtanbeam * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm'

   end subroutine constraint_eqn_020

   subroutine constraint_eqn_021(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for minor radius lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for minor radius lower limit
      !! #=# physics
      !! #=#=# frminor, aplasmin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! frminor : input real : f-value for minor radius limit
      !! rminor : input real : plasma minor radius (m)
      !! aplasmin : input real : minimum minor radius (m)
      use constraint_variables, only: frminor
      use physics_variables, only: rminor
      use build_variables, only: aplasmin
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - frminor * rminor/aplasmin
      tmp_con = aplasmin * (1.0D0 - tmp_cc)
      tmp_err = aplasmin * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_021

   subroutine constraint_eqn_022(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for divertor collision/connection length ratio upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for divertor collision/connection length ratio upper limit
      !! #=# divertor
      !! #=#=# fdivcol, rlenmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fdivcol : input real : f-value for divertor collisionality
      !! rlenmax : input real : maximum value for length ratio (rlclolcn)
      !! rlclolcn : input real : ratio of collision length / connection length
      use constraint_variables, only: fdivcol
      use divertor_variables, only: rlenmax, rlclolcn
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fdivcol * rlenmax / rlclolcn
      tmp_con = rlenmax * (1.0D0 - tmp_cc)
      tmp_err = rlclolcn * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_022

   subroutine constraint_eqn_023(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for conducting shell radius / rminor upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for conducting shell radius / rminor upper limit
      !! #=# physics
      !! #=#=# fcwr, cwrmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! rminor : input real : plasma minor radius (m)
      !! scraplo : input real : gap between plasma and first wall, outboard side (m)
      !! fwoth : input real : outboard first wall thickness, initial estimate (m)
      !! blnkoth : input real : outboard blanket thickness (m)
      !! fcwr : input real : f-value for conducting wall radius / rminor limit
      !! cwrmax : input real : maximum ratio of conducting wall distance to plasma minor radius for vertical stability
      use physics_variables, only: rminor, cwrmax
      use build_variables, only: scraplo, fwoth, blnkoth
      use constraint_variables, only: fcwr
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      ! rcw : local real : conducting shell radius (m)
      real(dp) :: rcw

      rcw = rminor + scraplo + fwoth + blnkoth
      tmp_cc =  1.0D0 - fcwr * cwrmax*rminor / rcw
      tmp_con = cwrmax*rminor * (1.0D0 - tmp_cc)
      tmp_err = rcw * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm'

   end subroutine constraint_eqn_023

   subroutine constraint_eqn_024(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for beta upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for beta upper limit
      !! #=# physics
      !! #=#=# fbetatry, betalim
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! iculbl : input integer : switch for beta limit scaling (constraint equation  24):<UL>
      !! <LI> = 0 apply limit to total beta;
      !! <LI> = 1 apply limit to thermal beta;
      !! <LI> = 2 apply limit to thermal + neutral beam beta
      !! <LI> = 3 apply limit to toroidal beta </UL>
      !! istell : input integer : switch for stellarator option (set via <CODE>device.dat</CODE>):<UL>
      !! <LI> = 0 use tokamak model;
      !! <LI> = 1 use stellarator model</UL>
      !! fbetatry : input real : f-value for beta limit
      !! betalim : input real : allowable beta
      !! beta : input real : total plasma beta (calculated if ipedestal =3)
      !! betaft : input real : fast alpha beta component
      !! betanb : input real : neutral beam beta component
      !! bt : input real : toroidal field
      !! btot : input real : total field
      use physics_variables, only: iculbl, betalim, beta, betanb, betaft, bt, btot
      use stellarator_variables, only: istell
      use constraint_variables, only: fbetatry
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! Include all beta components: relevant for both tokamaks and stellarators
      if ((iculbl == 0).or.(istell /= 0)) then
         tmp_cc =  1.0D0 - fbetatry * betalim/beta
         tmp_con = betalim
         tmp_err = betalim - beta / fbetatry
         tmp_symbol = '<'
         tmp_units = ''
      ! Here, the beta limit applies to only the thermal component, not the fast alpha or neutral beam parts
      else if (iculbl == 1) then
         tmp_cc = 1.0D0 - fbetatry * betalim/(beta-betaft-betanb)
         tmp_con = betalim
         tmp_err = betalim - (beta-betaft-betanb) / fbetatry
         tmp_symbol = '<'
         tmp_units = ''
      ! Beta limit applies to thermal + neutral beam: components of the total beta, i.e. excludes alphas
      else if (iculbl == 2) then
         tmp_cc = 1.0D0 - fbetatry * betalim/(beta-betaft)
         tmp_con = betalim * (1.0D0 - tmp_cc)
         tmp_err = (beta-betaft) * tmp_cc
         tmp_symbol = '<'
         tmp_units = ''
      ! Beta limit applies to toroidal beta
      else if (iculbl == 3) then
         tmp_cc =  1.0D0 - fbetatry * betalim/(beta*(btot/bt)**2)
         tmp_con = betalim
         tmp_err = betalim - (beta*(btot/bt)**2) / fbetatry
         tmp_symbol = '<'
         tmp_units = ''
      end if

   end subroutine constraint_eqn_024

   subroutine constraint_eqn_025(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for peak toroidal field upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for peak toroidal field upper limit
      !! #=# tfcoil
      !! #=#=# fpeakb, bmxlim
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpeakb : input real : f-value for maximum toroidal field
      !! bmxlim : input real : maximum peak toroidal field (T)
      !! bmaxtf : input real : mean peak field at TF coil (T)
      use constraint_variables, only: fpeakb, bmxlim
      use tfcoil_variables, only: bmaxtf
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fpeakb * bmxlim/bmaxtf
      tmp_con = bmxlim * (1.0D0 - tmp_cc)
      tmp_err = bmaxtf * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'T'

   end subroutine constraint_eqn_025

   subroutine constraint_eqn_026(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Central Solenoid current density upper limit at EOF
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for Central Solenoid current density upper limit at EOF
      !! #=# pfcoil
      !! #=#=# fjohc, rjohc
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fjohc : input real : f-value for central solenoid current at end-of-flattop
      !! rjohc : input real : allowable central solenoid current density at end of flat-top (A/m2)
      !! coheof : input real : central solenoid overall current density at end of flat-top (A/m2)
      use constraint_variables, only: fjohc
      use pfcoil_variables, only: rjohc, coheof
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fjohc * rjohc/coheof
      tmp_con = rjohc
      tmp_err = rjohc - coheof / fjohc
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_026

   subroutine constraint_eqn_027(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Central Solenoid current density upper limit at BOP
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for Central Solenoid current density upper limit at BOP
      !! #=# pfcoil
      !! #=#=# fjohc0, rjohc0
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fjohc0 : input real : f-value for central solenoid current at beginning of pulse
      !! rjohc0 : input real : allowable central solenoid current density at beginning of pulse (A/m2)
      !! cohbop : input real : central solenoid overall current density at beginning of pulse (A/m2)
      use constraint_variables, only: fjohc0
      use pfcoil_variables, only: rjohc0, cohbop
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fjohc0 * rjohc0/cohbop
      tmp_con = rjohc0
      tmp_err = rjohc0 - cohbop / fjohc0
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_027

   subroutine constraint_eqn_028(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for fusion gain (big Q) lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for fusion gain (big Q) lower limit
      !! #=# physics
      !! #=#=# fqval, bigqmin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fqval : input real : pf-value for Q
      !! bigq : input real : Fusion gain; P_fusion / (P_injection + P_ohmic)
      !! bigqmin : input real : minimum fusion gain Q
      !! ignite : input integer : switch for ignition assumption:<UL>
      !! <LI> = 0 do not assume plasma ignition;
      !! <LI> = 1 assume ignited (but include auxiliary power in costs)</UL>
      !! Obviously, ignite must be zero if current drive is required.
      !! If ignite=1, any auxiliary power is assumed to be used only
      !! during plasma start-up, and is excluded from all steady-state
      !! power balance calculations.
      use constraint_variables, only: fqval, bigqmin
      use current_drive_variables, only: bigq
      use physics_variables, only: ignite
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if plasma is not ignited ...
      if (ignite == 0) then
         tmp_cc =  1.0D0 - fqval * bigq/bigqmin
         tmp_con = bigqmin * (1.0D0 - tmp_cc)
         tmp_err = bigqmin * tmp_cc
         tmp_symbol = '>'
         tmp_units = ''
      ! if plasma is ignited report error
      else
         tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''
         call report_error(4)
      end if

   end subroutine constraint_eqn_028

   subroutine constraint_eqn_029(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for inboard major radius: This is a consistency equation
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for inboard major radius: This is a consistency equation
      !! #=# build
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! rmajor : input real : plasma major radius (m) (iteration variable 3)
      !! rminor : input real : plasma minor radius (m)
      !! rinboard : input real : plasma inboard radius (m)
      use physics_variables, only: rmajor, rminor
      use build_variables, only: rinboard
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - (rmajor - rminor) / rinboard
      tmp_con = rinboard * (1.0D0 - tmp_cc)
      tmp_err = rinboard * tmp_cc
      tmp_symbol = '='
      tmp_units = 'm'

   end subroutine constraint_eqn_029

   subroutine constraint_eqn_030(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for injection power upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for injection power upper limit
      !! #=# current_drive
      !! #=#=# fpinj, pinjalw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! pinjmw : input real : total auxiliary injected power (MW)
      !! fpinj : input real : f-value for injection power
      !! pinjalw : input real : Maximum allowable value for injected power (MW)
      use current_drive_variables, only: pinjmw, pinjalw
      use constraint_variables, only: fpinj
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - pinjmw / (fpinj * pinjalw)
      tmp_con = pinjalw
      tmp_err = pinjalw  - pinjmw * fpinj
      tmp_symbol = '<'
      tmp_units = 'MW'

   end subroutine constraint_eqn_030

   subroutine constraint_eqn_031(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil case stress upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil case stress upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fstrcase, sig_tf_case_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fstrcase : input real : f-value for TF coil case stress
      !! sig_tf_case_max : input real : Allowable maximum shear stress in TF coil case (Tresca criterion) (Pa)
      !! sig_tf_case : input real : Constrained stress in TF coil case (Pa)
      use constraint_variables, only: fstrcase
      use tfcoil_variables, only: sig_tf_case_max, sig_tf_case
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fstrcase * sig_tf_case_max/sig_tf_case
      tmp_con = sig_tf_case_max
      tmp_err = sig_tf_case_max - sig_tf_case / fstrcase
      tmp_symbol = '<'
      tmp_units = 'Pa'

   end subroutine constraint_eqn_031

   subroutine constraint_eqn_032(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil conduit stress upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil conduit stress upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fstrcond, sig_tf_wp_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fstrcond : input real : f-value for TF coil conduit stress
      !! sig_tf_wp_max : input real : Allowable maximum shear stress in TF coil conduit (Tresca criterion) (Pa)
      !! sig_tf_wp : input real : Constrained stress in TF conductor conduit (Pa)
      use constraint_variables, only: fstrcond
      use tfcoil_variables, only: sig_tf_wp_max, sig_tf_wp
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fstrcond * sig_tf_wp_max/sig_tf_wp
      tmp_con = sig_tf_wp_max
      tmp_err = sig_tf_wp_max - sig_tf_wp / fstrcond
      tmp_symbol = '<'
      tmp_units = 'Pa'

   end subroutine constraint_eqn_032

   subroutine constraint_eqn_033(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil operating/critical J upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil operating/critical J upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fiooic, jwdgcrt
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fiooic : input real : f-value for TF coil operating current / critical
      !! jwdgcrt : input real : critical current density for winding pack (A/m2)
      !! jwptf : input real : winding pack current density (A/m2)
      use constraint_variables, only: fiooic
      use tfcoil_variables, only: jwdgcrt, jwptf
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fiooic * jwdgcrt/jwptf
      tmp_con = jwdgcrt * (1.0D0 - tmp_cc)
      tmp_err = jwptf * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_033

   subroutine constraint_eqn_034(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil dump voltage upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil dump voltage upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fvdump, vdalw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fvdump : input real : f-value for dump voltage
      !! vdalw : input real : max voltage across TF coil during quench (kV)
      !! vtfskv : input real : voltage across a TF coil during quench (kV)
      use constraint_variables, only: fvdump
      use tfcoil_variables, only: vdalw, vtfskv
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fvdump * vdalw/vtfskv
      tmp_con = vdalw
      tmp_err = vdalw - vtfskv
      tmp_symbol = '<'
      tmp_units = 'V'

   end subroutine constraint_eqn_034

   subroutine constraint_eqn_035(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil J_wp/J_prot upper limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil J_wp/J_prot upper limit (SCTF)
      !! #=# tfcoil
      !! #=#=# fjprot, jwdgpro
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fjprot : input real : f-value for TF coil winding pack current density
      !! jwdgpro : input real : allowable TF coil winding pack current density, for dump temperature rise protection (A/m2)
      !! jwptf : input real : winding pack current density (A/m2)
      use constraint_variables, only: fjprot
      use tfcoil_variables, only: jwdgpro, jwptf
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fjprot * jwdgpro/jwptf
      tmp_con = jwdgpro
      tmp_err = jwptf - jwptf
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_035

   subroutine constraint_eqn_036(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil s/c temperature margin lower limit (SCTF)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil s/c temperature margin lower limit (SCTF)
      !! #=# tfcoil
      !! #=#=# ftmargtf, tmargmin_tf
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftmargtf : input real : f-value for TF coil temperature margin
      !! tmargtf : input real : TF coil temperature margin (K)
      !! tmargmin_tf : input real : minimum allowable temperature margin : TF coils (K)
      use constraint_variables, only: ftmargtf
      use tfcoil_variables, only: tmargtf, tmargmin_tf
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ftmargtf * tmargtf/tmargmin_tf
      tmp_con = tmargmin_tf
      tmp_err = tmargmin_tf - tmargtf
      tmp_symbol = '>'
      tmp_units = 'K'

   end subroutine constraint_eqn_036

   subroutine constraint_eqn_037(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for current drive gamma upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for current drive gamma upper limit
      !! #=# current_drive
      !! #=#=# fgamcd, gammax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fgamcd : input real : f-value for current drive gamma
      !! gammax : input real : maximum current drive gamma
      !! gamcd : input real : normalised current drive efficiency (1.0e20 A/W-m2)
      use constraint_variables, only: fgamcd, gammax
      use current_drive_variables, only: gamcd
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fgamcd * gammax/gamcd
      tmp_con = gammax * (1.0D0 - tmp_cc)
      tmp_err = gamcd * tmp_cc
      tmp_symbol = '<'
      tmp_units = '1E20 A/Wm2'

   end subroutine constraint_eqn_037

   subroutine constraint_eqn_038(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Obsolete
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Obsolete
      !! #=# empty
      !! #=#=# empty
      implicit none
	   ! Dummy formal arguments, for compliance with interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_038

   subroutine constraint_eqn_039(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for first wall temperature upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for first wall temperature upper limit
      !! #=# fwbs
      !! #=#=# ftpeak, tfwmatmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftpeak : input real : f-value for first wall peak temperature
      !! tfwmatmax : input real : maximum temperature of first wall material (K) (secondary_cycle>1)
      !! tpeak : input real : peak first wall temperature (K)
      use constraint_variables, only: ftpeak
      use fwbs_variables, only: tfwmatmax, tpeak
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! If the temperature peak == 0 then report an error
      if (tpeak < 1.0D0) call report_error(5)
      tmp_cc =  1.0D0 - ftpeak * tfwmatmax/tpeak
      tmp_con = tfwmatmax * (1.0D0 - tmp_cc)
      tmp_err = tpeak * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'K'

   end subroutine constraint_eqn_039

   subroutine constraint_eqn_040(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for auxiliary power lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for auxiliary power lower limit
      !! #=# current_drive
      !! #=#=# fauxmn, auxmin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fauxmn : input real : f-value for minimum auxiliary power
      !! pinjmw : input real : total auxiliary injected power (MW)
      !! auxmin : input real : minimum auxiliary power (MW)
      use constraint_variables, only: fauxmn, auxmin
      use current_drive_variables, only: pinjmw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fauxmn * pinjmw/auxmin
      tmp_con = auxmin * (1.0D0 - tmp_cc)
      tmp_err = auxmin * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'MW'

   end subroutine constraint_eqn_040

   subroutine constraint_eqn_041(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for plasma current ramp-up time lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for plasma current ramp-up time lower limit
      !! #=# times
      !! #=#=# ftohs, tohsmn
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftohs : input real : f-value for plasma current ramp-up time
      !! tohs : input real : plasma current ramp-up time for current initiation (s)
      !! tohsmn : input real : minimum plasma current ramp-up time (s)
      use constraint_variables, only: ftohs, tohsmn
      use times_variables, only: tohs
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ftohs * tohs/tohsmn
      tmp_con = tohsmn * (1.0D0 - tmp_cc)
      tmp_err = tohsmn * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'sec'

   end subroutine constraint_eqn_041

   subroutine constraint_eqn_042(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for cycle time lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for cycle time lower limit
      !! #=# times
      !! #=#=# ftcycl, tcycmn
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftcycl : input real : f-value for cycle time
      !! tcycle : input real : full cycle time (s)
      !! tcycmn : input real : minimum cycle time (s)
      use constraint_variables, only: ftcycl, tcycmn
      use times_variables, only: tcycle
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the minimum cycle time == 0 report an error
      if (tcycmn < 1.0D0) call report_error(6)
      tmp_cc =  1.0D0 - ftcycl * tcycle/tcycmn
      tmp_con = tcycmn * (1.0D0 - tmp_cc)
      tmp_err = tcycmn * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'sec'

   end subroutine constraint_eqn_042

   subroutine constraint_eqn_043(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for average centrepost temperature: This is a consistency equation (TART)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for average centrepost temperature: This is a consistency equation (TART)
      !! #=# tfcoil
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! tcpav : input real : average temp of TF coil inboard leg conductor (C)e
      !! tcpav2 : input real : centrepost average temperature (C) (for consistency)
      !! itart : input integer : switch for spherical tokamak (ST) models:<UL>
      !! <LI> = 0 use conventional aspect ratio models;
      !! <LI> = 1 use spherical tokamak models</UL>
      use tfcoil_variables, only: tcpav, tcpav2
      use physics_variables, only: itart
      use tfcoil_variables, only:  i_tf_sup

      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the machine isn't a ST then report error
      if (itart == 0) call report_error(7)

      ! For some reasons these lines are needed to make VMCON CONVERGE ....
      if ( i_tf_sup == 0 ) then ! Copper case
         tcpav = tcpav - 273.15D0
         tcpav2 = tcpav2 - 273.15D0
      end if

      tmp_cc =   1.0D0 - tcpav/tcpav2
      tmp_con = tcpav2 * (1.0D0 - tmp_cc)
      tmp_err = tcpav2 * tmp_cc
      tmp_symbol = '='
      tmp_units = 'deg C'

      ! For some reasons these lines are needed to make VMCON CONVERGE ....
      if ( i_tf_sup == 0 ) then ! Copper case
         tcpav = tcpav + 273.15D0
         tcpav2 = tcpav2 + 273.15D0
      end if



   end subroutine constraint_eqn_043

   subroutine constraint_eqn_044(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for centrepost temperature upper limit (TART)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for centrepost temperature upper limit (TART)
      !! #=# tfcoil
      !! #=#=# fptemp, ptempalw
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fptemp : input real : f-value for peak centrepost temperature
      !! ptempalw : input real : maximum peak centrepost temperature (K)
      !! tcpmax : input real :  peak centrepost temperature (K)
      !! itart : input integer : switch for spherical tokamak (ST) models:<UL>
      !! <LI> = 0 use conventional aspect ratio models;
      !! <LI> = 1 use spherical tokamak models</UL>
      use constraint_variables, only: fptemp
      use tfcoil_variables, only: ptempalw, tcpmax
      use physics_variables, only: itart
      use tfcoil_variables, only:  i_tf_sup
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the machine isn't a ST then report error
      if (itart == 0) call report_error(8)

      ! For some reasons these lines are needed to make VMCON CONVERGE ....
      if ( i_tf_sup == 0 ) then ! Copper case
         ptempalw = ptempalw - 273.15D0
         tcpmax = tcpmax - 273.15D0
      end if

      tmp_cc =   1.0D0 - fptemp * ptempalw / tcpmax
      tmp_con = ptempalw * (1.0D0 - tmp_cc)
      tmp_err = tcpmax * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'deg C'

      ! For some reasons these lines are needed to make VMCON CONVERGE ....
      if ( i_tf_sup == 0 ) then ! Copper case
         ptempalw = ptempalw + 273.15D0
         tcpmax = tcpmax + 273.15D0
      end if

   end subroutine constraint_eqn_044

   subroutine constraint_eqn_045(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for edge safety factor lower limit (TART)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for edge safety factor lower limit (TART)
      !! #=# tfcoil
      !! #=#=# fq, qlim
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fq : input real : f-value for edge safety factor
      !! q : safety factor 'near' plasma edge: equal to q95
      !! (unless icurr = 2 (ST current scaling), in which case q = mean edge safety factor qbar)
      !! qlim : input real :  lower limit for edge safety factor
      !! itart : input integer : switch for spherical tokamak (ST) models:<UL>
      !! <LI> = 0 use conventional aspect ratio models;
      !! <LI> = 1 use spherical tokamak models</UL>
      use constraint_variables, only: fq
      use physics_variables, only: q, qlim, itart
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the machine isn't a ST then report error
      if (itart == 0) call report_error(9)
      tmp_cc =   1.0D0 - fq * q/qlim
      tmp_con = qlim * (1.0D0 - tmp_cc)
      tmp_err = qlim * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_045

   subroutine constraint_eqn_046(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Ip/Irod upper limit (TART)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for Ip/Irod upper limit (TART)
      !! #=# tfcoil
      !! #=#=# fipir, cratmx
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! eps : input real :  inverse aspect ratio
      !! fipir : input real : f-value for Ip/Irod upper limit
      !! ritfc : input real : total (summed) current in TF coils (A)
      !! plascur : input real :  plasma current (A)
      !! itart : input integer : switch for spherical tokamak (ST) models:<UL>
      !! <LI> = 0 use conventional aspect ratio models;
      !! <LI> = 1 use spherical tokamak models</UL>
      use physics_variables, only: eps, plascur, itart
      use constraint_variables, only: fipir
      use tfcoil_variables, only: ritfc
      implicit none
      ! cratmx : local real : maximum ratio of plasma current to centrepost current
      real(dp) :: cratmx
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! if the machine isn't a ST then report error
      if (itart == 0) call report_error(10)
      cratmx = 1.0D0 + 4.91D0*(eps-0.62D0)
      tmp_cc =  1.0D0 - fipir * cratmx * ritfc/plascur
      tmp_con = cratmx * (1.0D0 - tmp_cc)
      tmp_err = plascur/ritfc * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_046

   subroutine constraint_eqn_047(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Issue #508 Remove RFP option: Relevant only to reversed field pinch devices
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Issue #508 Remove RFP option: Relevant only to reversed field pinch devices
      !! Equation for TF coil toroidal thickness upper limit
      !! #=# empty
      !! #=#=# empty
      implicit none
      ! Dummy formal arguments, just to comply with the subroutine interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_047

   subroutine constraint_eqn_048(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for poloidal beta upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for poloidal beta upper limit
      !! #=# physics
      !! #=#=# fbetap, betpmx
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fbetap : input real : rf-value for poloidal beta
      !! betpmx : input real :  maximum poloidal beta
      !! betap : input real :  poloidal beta
      use constraint_variables, only: fbetap, betpmx
      use physics_variables, only: betap
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fbetap * betpmx/betap
      tmp_con = betpmx * (1.0D0 - tmp_cc)
      tmp_err = betap * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_048

   subroutine constraint_eqn_049(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Issue #508 Remove IFE option: Equation for repetition rate upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Issue #508 Remove IFE option: Equation for repetition rate upper limit
      !! #=# empty
      !! #=#=# empty
      ! Dummy formal arguments, just to comply with the subroutine interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_049

   subroutine constraint_eqn_050(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! IFE option: Equation for repetition rate upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! author: S I Muldrew, CCFE, Culham Science Centre
      !! IFE option: Equation for repetition rate upper limit
      !! #=# IFE
      !! #=#=# frrmax, rrmax
      use ife_variables, only: frrmax, ife, rrmax, reprat
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      if (ife /= 1) then
         call report_error(12)
      end if

      tmp_cc =  1.0D0 - frrmax * rrmax/reprat
      tmp_con = rrmax * (1.0D0 - tmp_cc)
      tmp_err = reprat * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'Hz'

   end subroutine constraint_eqn_050

   subroutine constraint_eqn_051(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation to enforce startup flux = available startup flux
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation to enforce startup flux = available startup flux
      !! #=# pfcoil
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! vsres : input real : resistive losses in startup V-s (Wb)
      !! vsind : input real :  internal and external plasma inductance V-s (Wb))
      !! vssu : input real :  total flux swing for startup (Wb)
      use physics_variables, only: vsres, vsind
      use pfcoil_variables, only: vssu, fvssu
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fvssu * abs((vsres+vsind) / vssu)
      tmp_con = vssu * (1.0D0 - tmp_cc)
      tmp_err = vssu * tmp_cc
      tmp_symbol = '='
      tmp_units = 'V.s'

   end subroutine constraint_eqn_051

   subroutine constraint_eqn_052(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for tritium breeding ratio lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for tritium breeding ratio lower limit
      !! #=# fwbs
      !! #=#=# ftbr, tbrmin
      !! ? TODO should this only be for certain blanket models ?
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftbr : input real : f-value for minimum tritium breeding ratio
      !! tbr : input real :  tritium breeding ratio (iblanket=2,3 (KIT HCPB/HCLL))
      !! tbrmin : input real :  minimum tritium breeding ratio (If iblanket=1, tbrmin=minimum 5-year time-averaged tritium breeding ratio)
      use constraint_variables, only: ftbr, tbrmin
      use fwbs_variables, only: tbr
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ftbr * tbr/tbrmin
      tmp_con = tbrmin * (1.0D0 - tmp_cc)
      tmp_err = tbrmin * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_052

   subroutine constraint_eqn_053(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for fast neutron fluence on TF coil upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for fast neutron fluence on TF coil upper limit
      !! #=# fwbs
      !! #=#=# fflutf, nflutfmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fflutf : input real : f-value for maximum TF coil nuclear heating
      !! nflutfmax : input real :  max fast neutron fluence on TF coil (n/m2)
      !! nflutf : input real :  peak fast neutron fluence on TF coil superconductor (n/m2)
      use constraint_variables, only: fflutf, nflutfmax
      use fwbs_variables, only: nflutf
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fflutf * nflutfmax/nflutf
      tmp_con = nflutfmax * (1.0D0 - tmp_cc)
      tmp_err = nflutf * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'neutron/m2'

   end subroutine constraint_eqn_053

   subroutine constraint_eqn_054(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for peak TF coil nuclear heating upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for peak TF coil nuclear heating upper limit
      !! #=# fwbs
      !! #=#=# fptfnuc, ptfnucmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fptfnuc : input real : f-value for maximum TF coil nuclear heating
      !! ptfnucmax : input real :  maximum nuclear heating in TF coil (MW/m3)
      !! ptfnucpm3 : input real :  nuclear heating in the TF coil (MW/m3) (blktmodel>0)
      use constraint_variables, only: fptfnuc, ptfnucmax
      use fwbs_variables, only: ptfnucpm3
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - fptfnuc * ptfnucmax/ptfnucpm3
      tmp_con = ptfnucmax * (1.0D0 - tmp_cc)
      tmp_err = ptfnucpm3 * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m3'

   end subroutine constraint_eqn_054

   subroutine constraint_eqn_055(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! vvhemax is no longer calculated in PROCESS and this constraint is disabled
      implicit none

      real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      call report_error(173)
   end subroutine constraint_eqn_055

   subroutine constraint_eqn_056(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for power through separatrix / major radius upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for power through separatrix / major radius upper limit
      !! #=# current_drive
      !! #=#=# fnbshinef, nbshinefmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpsepr : input real : f-value for maximum Psep/R limit
      !! pseprmax : input real :  maximum ratio of power crossing the separatrix to plasma major radius (Psep/R) (MW/m)
      !! pdivt : input real :  power to be conducted to the divertor region (MW)
      !! rmajor : input real :  plasma major radius (m)
      use constraint_variables, only: fpsepr, pseprmax
      use physics_variables, only: pdivt, rmajor
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - fpsepr * pseprmax / (pdivt/rmajor)
      tmp_con = pseprmax * (1.0D0 - tmp_cc)
      tmp_err = (pdivt/rmajor) * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m'

   end subroutine constraint_eqn_056

   subroutine constraint_eqn_057(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Obsolete
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Obsolete
      !! #=# empty
      !! #=#=# empty
      ! Dummy formal arguments, just to comply with the subroutine interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_057

   subroutine constraint_eqn_058(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Obsolete
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! Obsolete
      !! #=# empty
      !! #=#=# empty
      ! Dummy formal arguments, just to comply with the subroutine interface
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 0
        tmp_con = 0
        tmp_err = 0
        tmp_symbol = ''
        tmp_units = ''

   end subroutine constraint_eqn_058

   subroutine constraint_eqn_059(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for neutral beam shine-through fraction upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for neutral beam shine-through fraction upper limit
      !! #=# current_drive
      !! #=#=# fnbshinef, nbshinefmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fnbshinef : input real : f-value for maximum neutral beam shine-through fraction
      !! nbshinefmax : input real :  maximum neutral beam shine-through fraction
      !! nbshinef : input real :  neutral beam shine-through fraction
      use constraint_variables, only: fnbshinef, nbshinefmax
      use current_drive_variables, only: nbshinef
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      tmp_cc = 1.0D0 - fnbshinef * nbshinefmax / nbshinef
      tmp_con = nbshinefmax * (1.0D0 - tmp_cc)
      tmp_err = nbshinef * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''
   end subroutine constraint_eqn_059

   subroutine constraint_eqn_060(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Central Solenoid s/c temperature margin lower limi
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for Central Solenoid s/c temperature margin lower limi
      !! #=# tfcoil
      !! #=#=# ftmargoh, tmargmin_cs
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftmargoh : input real :  f-value for central solenoid temperature margin
      !! tmargoh : input real :  Central solenoid temperature margin (K)
      !! tmargmin_cs : input real :  Minimum allowable temperature margin : CS (K)
      use constraint_variables, only: ftmargoh
      use pfcoil_variables, only: tmargoh
      use tfcoil_variables, only: tmargmin_cs
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - ftmargoh * tmargoh/tmargmin_cs
      tmp_con = tmargmin_cs
      tmp_err = tmargmin_cs - tmargoh
      tmp_symbol = '>'
      tmp_units = 'K'

   end subroutine constraint_eqn_060

   subroutine constraint_eqn_061(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for availability limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for availability limit
      !! #=# cost
      !! #=#=# favail, avail_min
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! favail : input real : F-value for minimum availability
      !! cfactr : input real : Total plant availability fraction
      !! avail_min : input real : Minimum availability
      use cost_variables, only: favail, cfactr, avail_min
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - favail * cfactr / avail_min
      tmp_con = avail_min * (1.0D0 - tmp_cc)
      tmp_err = cfactr * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_061

   subroutine constraint_eqn_062(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Lower limit on taup/taueff the ratio of alpha particle to energy confinement times
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Lower limit on taup/taueff the ratio of alpha particle to energy confinement times
      !! #=# physics
      !! #=#=# ftaulimit, taulimit
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftaulimit : input real : f-value for lower limit on taup/taueff the ratio of alpha particle to energy confinement
      !! taup : input real : alpha particle confinement time (s)
      !! taueff : input real : global thermal energy confinement time (sec)
      !! taulimit : input real : Lower limit on taup/taueff the ratio of alpha particle to energy confinement times
      use constraint_variables, only: ftaulimit, taulimit
      use physics_variables, only: taup, taueff
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - ftaulimit * (taup / taueff) / taulimit
      tmp_con = taulimit
      tmp_err = (taup / taueff) * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_062

   subroutine constraint_eqn_063(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit on niterpump (vacuum_model = simple)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Upper limit on niterpump (vacuum_model = simple)
      !! #=# vacuum
      !! #=#=# fniterpump, tfno
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fniterpump : input real : f-value for constraint that number of pumps < tfno
      !! tfno : input real : number of TF coils (default = 50 for stellarators)
      !! niterpump : input real : number of high vacuum pumps (real number), each with the throughput
      use constraint_variables, only: fniterpump
      use tfcoil_variables, only: n_tf
      use vacuum_variables, only: niterpump
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - fniterpump * n_tf / niterpump
      tmp_con = n_tf
      tmp_err = n_tf * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_063

   subroutine constraint_eqn_064(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Upper limit on Zeff
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Upper limit on Zeff
      !! #=# physics
      !! #=#=# fzeffmax, zeffmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fzeffmax : input real : f-value for maximum zeff
      !! zeffmax : input real : maximum value for Zeff
      !! zeff : input real : plasma effective charge
      use constraint_variables, only: fzeffmax, zeffmax
      use physics_variables, only: zeff
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - fzeffmax * (zeffmax/zeff)
      tmp_con = zeffmax
      tmp_err = zeffmax * tmp_cc
      tmp_symbol = '<'
      tmp_units = ''

   end subroutine constraint_eqn_064

   subroutine constraint_eqn_065(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Limit the stress of the vacuum vessel that occurs when the TF coil quenches.
      !! author: Timothy Nunn, UKAEA
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fmaxvvstress : input real : f-value for constraint on maximum VV stress
      !! max_vv_stress : input real : Maximum permitted stress of the VV (Pa)
      !! vv_stress_quench : input real : Stress of the VV (Pa)
      use constraint_variables, only: fmaxvvstress
      use tfcoil_variables, only: max_vv_stress
      use sctfcoil_module, only: vv_stress_quench
      implicit none
      real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0d0 - fmaxvvstress * max_vv_stress / vv_stress_quench
      tmp_con = max_vv_stress
      tmp_err = max_vv_stress * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'Pa'

   end subroutine constraint_eqn_065

   subroutine constraint_eqn_066(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Limit on rate of change of energy in poloidal field
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Limit on rate of change of energy in poloidal field
      !! #=# pfcoil
      !! #=#=# fpoloidalpower, maxpoloidalpower
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpoloidalpower : input real : f-value for constraint on rate of change of energy in poloidal field
      !! maxpoloidalpower : input real : Maximum permitted absolute rate of change of stored energy in poloidal field (MW)
      !! peakpoloidalpower : input real : Peak absolute rate of change of stored energy in poloidal field (MW) (11/01/16)
      use constraint_variables, only: fpoloidalpower
      use pf_power_variables, only: maxpoloidalpower, peakpoloidalpower
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0d0 - fpoloidalpower * maxpoloidalpower / peakpoloidalpower
      tmp_con = maxpoloidalpower
      tmp_err = maxpoloidalpower * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW'

   end subroutine constraint_eqn_066

   subroutine constraint_eqn_067(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Simple upper limit on radiation wall load
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Simple upper limit on radiation wall load
      !! #=# physics
      !! #=#=# fradwall, maxradwallload
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fradwall : input real : f-value for upper limit on radiation wall load
      !! maxradwallload : input real : Maximum permitted radiation wall load (MW/m^2)
      !! peakradwallload : input real : Peak radiation wall load (MW/m^2)
      use constraint_variables, only: fradwall, maxradwallload, peakradwallload
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0d0 - fradwall * maxradwallload / peakradwallload
      tmp_con = maxradwallload
      tmp_err =  maxradwallload * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW/m^2'

   end subroutine constraint_eqn_067

   subroutine constraint_eqn_068(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! New Psep scaling (PsepB/qAR)
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! New Psep scaling (PsepB/qAR)
      !! Issue #464
      !! #=# physics
      !! #=#=# fpsepbqar, psepbqarmax
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpsepbqar : input real : f-value for upper limit on psepbqar, maximum Psep*Bt/qAR limit
      !! psepbqarmax : input real : maximum permitted value of ratio of Psep*Bt/qAR (MWT/m)
      !! pdivt : input real : Power to conducted to the divertor region (MW)
      !! bt : input real : toroidal field on axis (T) (iteration variable 2)
      !! q95 : input real : safety factor q at 95% flux surface
      !! aspect : input real : aspect ratio (iteration variable 1)
      !! rmajor : input real : plasma major radius (m) (iteration variable 3)
      use constraint_variables, only: fpsepbqar, psepbqarmax
      use physics_variables, only: pdivt, bt, q95, aspect, rmajor
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0d0 - fpsepbqar * psepbqarmax / ((pdivt*bt)/(q95*aspect*rmajor))
      tmp_con = psepbqarmax
      tmp_err = (pdivt*bt)/(q95*aspect*rmajor) - psepbqarmax
      tmp_symbol = '<'
      tmp_units = 'MWT/m'

   end subroutine constraint_eqn_068

   subroutine constraint_eqn_069(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Ensure separatrix power is less than value from Kallenbach divertor
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Ensure separatrix power is less than value from Kallenbach divertor
      !! #=# divertor_kallenbach
      !! #=#=# consistency, psep_kallenbach
      !! fpsep has been removed from the equation.
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! psep_kallenbach : input real : Power conducted through the separatrix, as calculated by the divertor model [W]
      !! pdivt : input real :  power to conducted to the divertor region (MW)
      ! use div_kal_vars, only: psep_kallenbach
      use physics_variables, only: pdivt
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! From Kallenbach model, should be reserved if the model is going to be added back

   end subroutine constraint_eqn_069

   subroutine constraint_eqn_070(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Separatrix density consistency
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Separatrix density consistency
      !! #=# divertor_kallenbach
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! teomp : input real : Separatrix temperature calculated by the Kallenbach divertor model [eV]
      !! tesep : input real : Electron temperature at separatrix [keV]
      ! use div_kal_vars, only: teomp
      use  physics_variables, only: tesep
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! From Kallenbach model, should be reserved if the model is going to be added back

   end subroutine constraint_eqn_070

   subroutine constraint_eqn_071(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Separatrix density consistency
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Separatrix density consistency
      !! #=# divertor_kallenbach
      !! #=#=# consistency
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! neomp : input real : Mean SOL density at OMP calculated by the Kallenbach divertor model [m-3]
      !! nesep : input real :  electron density at separatrix [m-3] (ipedestal=1,2, calculated if 3)
      !! neratio : input real : Ratio of mean SOL density at OMP to separatrix density at OMP (iteration variable 121)
      ! use div_kal_vars, only: neomp, neratio
      use physics_variables, only: nesep
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! From Kallenbach model, should be reserved if the model is going to be added back

   end subroutine constraint_eqn_071

   subroutine constraint_eqn_072(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Central Solenoid Tresca yield criterion
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Central Solenoid Tresca yield criterion
      !! #=# pfcoil
      !! #=#=# foh_stress, alstroh
      !! In the case if the bucked and wedged option ( i_tf_bucking >= 2 ) the constrained
      !! stress is the largest the largest stress of the
      !!  - CS stress at maximum current (conservative as the TF inward pressure is not taken
      !!    into account)
      !!  - CS stress at flux swing (no current in CS) from the TF inward pressure
      !! This allow to cover the 2 worst stress scenario in the bucked and wedged design
      !! Otherwise (free standing TF), the stress limits are only set by the CS stress at max current
      !! Reverse the sign so it works as an inequality constraint (tmp_cc > 0)
      !! This will have no effect if it is used as an equality constraint because it will be squared.
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! foh_stress : input real : f-value for Tresca yield criterion in Central Solenoid
      !! alstroh : input real :  allowable hoop stress in Central Solenoid structural material (Pa)
      !! s_tresca_oh : input real : Maximum shear stress coils/central solenoid (Pa)
      !! sig_tf_cs_bucked : input real : Maximum shear stress in CS case at flux swing (no current in CS)
      !!                       can be significant for the bucked and weged design
      !! i_tf_bucking : input integer : switch for TF structure design
      use constraint_variables, only: foh_stress
      use pfcoil_variables, only: alstroh, s_tresca_oh
      use tfcoil_variables, only: sig_tf_cs_bucked, i_tf_bucking
      use build_variables, only: tf_in_cs
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! bucked and wedged desing (see subroutine comment)
      if ( i_tf_bucking >= 2 .and. tf_in_cs == 0 ) then
         tmp_cc = 1.0d0 - foh_stress * alstroh / max(s_tresca_oh, sig_tf_cs_bucked)
         tmp_err = alstroh - max(s_tresca_oh, sig_tf_cs_bucked)
      ! Free standing CS
      else
         tmp_cc = 1.0d0 - foh_stress * alstroh / s_tresca_oh
         tmp_err = alstroh - s_tresca_oh
      end if

      tmp_con = alstroh
      tmp_symbol = '<'
      tmp_units = 'Pa'

   end subroutine constraint_eqn_072

   subroutine constraint_eqn_073(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Ensure separatrix power is greater than the L-H power + auxiliary power
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Ensure separatrix power is greater than the L-H power + auxiliary power
      !! #=# physics
      !! #=#=# fplhsep, pdivt
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fplhsep : input real : F-value for Psep >= Plh + Paux : for consistency of two values of separatrix power
      !! plhthresh : input real : L-H mode power threshold (MW)
      !! pdivt : input real : power to be conducted to the divertor region (MW)
      !! pinjmw : inout real : total auxiliary injected power (MW)
      use physics_variables, only: fplhsep, plhthresh, pdivt
      use current_drive_variables, only: pinjmw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0d0 - fplhsep * pdivt / (plhthresh+pinjmw)
      tmp_con = pdivt
      tmp_err = pdivt * tmp_cc
      tmp_symbol = '>'
      tmp_units = 'MW'

   end subroutine constraint_eqn_073

   subroutine constraint_eqn_074(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Ensure TF coil quench temperature < tmax_croco ONLY used for croco HTS coil
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Ensure TF coil quench temperature < tmax_croco ONLY used for croco HTS coil
      !! #=# physics
      !! #=#=# fcqt, tmax_croco
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fcqt : input real : f-value: TF coil quench temparature remains below tmax_croco
      !! croco_quench_temperature : input real : CroCo strand: Actual temp reached during a quench (K)
      !! tmax_croco : input real : CroCo strand: maximum permitted temp during a quench (K)
      use constraint_variables, only: fcqt
      use tfcoil_variables, only: croco_quench_temperature, tmax_croco
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0d0 - fcqt * tmax_croco / croco_quench_temperature
      tmp_con = croco_quench_temperature
      tmp_err = croco_quench_temperature * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'K'

   end subroutine constraint_eqn_074

   subroutine constraint_eqn_075(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Ensure that TF coil current / copper area < Maximum value
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Ensure that TF coil current / copper area < Maximum value
      !! ONLY used for croco HTS coil
      !! #=# physics
      !! #=#=# f_coppera_m2, copperA_m2_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! copperA_m2 : input real :
      !! copperA_m2_max : input real :
      !! f_coppera_m2 : input real :
      use rebco_variables, only: copperA_m2, copperA_m2_max, f_coppera_m2
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0d0 - f_coppera_m2 * copperA_m2_max / copperA_m2
      tmp_con = copperA_m2
      tmp_err = copperA_m2 * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_075

   subroutine constraint_eqn_076(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Eich critical separatrix density model: Added for issue 558
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Eich critical separatrix density model
      !! Added for issue 558 with ref to http://iopscience.iop.org/article/10.1088/1741-4326/aaa340/pdf
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! alpha_crit : output real : critical ballooning parameter value
      !! nesep_crit : output real : critical electron density at separatrix [m-3]
      !! kappa : input real : plasma separatrix elongation (calculated if ishape = 1-5, 7 or 9)
      !! triang : input real : plasma separatrix triangularity (calculated if ishape = 1, 3-5 or 7)
      !! aspect : input real : aspect ratio (iteration variable 1)
      !! pdivt : input real : power to conducted to the divertor region (MW)
      !! dlimit(7) : input real array : density limit (/m3) as calculated using various models
      !! fnesep : input real : f-value for Eich critical separatrix density
      use physics_variables, only: alpha_crit, nesep_crit, kappa, triang, &
                                   aspect, pdivt, dlimit, nesep
      use constraint_variables, only: fnesep
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      alpha_crit = (kappa ** 1.2D0) * (1.0D0 + 1.5D0 * triang)
      nesep_crit = 5.9D0 * alpha_crit * (aspect ** (-2.0D0/7.0D0)) * &
                (((1.0D0 + (kappa ** 2.0D0)) / 2.0D0) ** (-6.0D0/7.0D0)) &
                * ((pdivt* 1.0D6) ** (-11.0D0/70.0D0)) * dlimit(7)
      tmp_cc = 1.0D0 - fnesep * nesep_crit/nesep
      tmp_con = nesep
      tmp_err = nesep * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm-3'

   end subroutine constraint_eqn_076

   subroutine constraint_eqn_077(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for maximum TF current per turn upper limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value; residual error in physical units; output string; units string
      !! Equation for maximum TF current per turn upper limit
      !! #=# tfcoil
      !! #=#=# fcpttf, cpttf, cpttf_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fcpttf : input : f-value for TF coil current per turn
      !! cpttf_max  : input : allowable TF coil current per turn [A/turn]
      !! cpttf  : input : TF coil current per turn [A/turn]
      use constraint_variables, only: fcpttf
      use tfcoil_variables, only: cpttf_max, cpttf
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0D0 - fcpttf * cpttf_max/cpttf
      tmp_con = cpttf_max
      tmp_err = cpttf_max * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'A/turn'

   end subroutine constraint_eqn_077

   subroutine constraint_eqn_078(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for Reinke criterion, divertor impurity fraction lower limit
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value; residual error in physical units; output string; units string
      !! Equation for Reinke criterion, divertor impurity fraction lower limit
      !! #=# divertor
      !! #=#=# freinke, fzactual, fzmin
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present;
      !! and con will be printed out only if present. Thesw conditions were missing.
      !! freinke : input : f-value for Reinke criterion (itv 147)
      !! fzmin : input : minimum impurity fraction from Reinke model
      !! fzactual : input : actual impurity fraction
      use constraint_variables, only: freinke
      use reinke_variables, only: fzactual, fzmin
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! write(*,*) 'freinke, fzact, fzmin = ', freinke, ', ', fzactual, ', ', fzmin
      !            1.0,    0.0,   value
      tmp_cc = 1.0D0 - freinke *  fzactual/fzmin
      !The following two pre-existing lines are not understood:
      !KE note - cc is always 1, code never enters IF statement...
      tmp_con = fzmin * (1.0D0 - tmp_cc)
      tmp_err = fzmin * tmp_cc
      tmp_symbol = '>'
      tmp_units  = ''
      ! write(*,*) 'cc, con = ', tmp_cc, ', ', tmp_con
      ! write(*,*) 'freinke, fzactual, fzmin = ', freinke, ', ', fzactual, ', ', fzmin

   end subroutine constraint_eqn_078

   subroutine constraint_eqn_079(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for maximum CS field
      !! author: P B Lloyd, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value; residual error in physical units; output string; units string
      !! Equation for maximum CS field
      !! #=# pfcoil
      !! #=#=# fbmaxcs, bmaxoh, bmaxoh0, bmaxcs_lim
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fbmaxcs : input : F-value for CS mmax field (cons. 79, itvar 149)
      !! bmaxcs_lim : input : Central solenoid max field limit [T]
      !! bmaxoh0 : input : maximum field in central solenoid at beginning of pulse (T)
      !! bmaxoh : input real : maximum field in central solenoid at end of flat-top (EoF) (T)
      !! (Note: original code has "bmaxoh/bmaxoh0 |  peak CS field [T]".)
      use pfcoil_variables, only: fbmaxcs, bmaxcs_lim, bmaxoh0, bmaxoh
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc     = 1.0D0 - fbmaxcs * bmaxcs_lim/max(bmaxoh, bmaxoh0)
      tmp_con    = bmaxcs_lim
      tmp_err    = max(bmaxoh, bmaxoh0) * tmp_cc
      tmp_symbol = '<'
      tmp_units  = 'A/turn'

   end subroutine constraint_eqn_079

   subroutine constraint_eqn_080(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for pdivt lower limit
      !! author: J Morris, Culham Science Centre
      !! args : output structure : residual error; constraint value; residual error in physical units;
      !! output string; units string
      !! Lower limit pdivt
      !! #=# physics
      !! #=#=# fpdivlim, pdivt
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fpdivlim : input : F-value for lower limit on pdivt (cons. 80, itvar 153)
      !! pdivtlim : input : Minimum power crossing separatrix pdivt [MW]
      !! pdivt : input : Power crossing separatrix [MW]
      use physics_variables, only: fpdivlim, pdivt
      use constraint_variables, only : pdivtlim
      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      tmp_cc     = 1.0D0 - fpdivlim * pdivt / pdivtlim
      tmp_con    = pdivtlim
      tmp_err    = pdivt * tmp_cc
      tmp_symbol = '>'
      tmp_units  = 'MW'

   end subroutine constraint_eqn_080

   subroutine constraint_eqn_081(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Make sure that the central density is larger that the pedestal one
      !! author: S Kahn, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Lower limit ne0 > neped
      !! !#=# physics
      !! !#=#=# ne0, neped
      !! Logic change during pre-factoring: err, symbol, units will be
      !! assigned only if present.
      !! fne0  : input : F-value for constraint on ne0 > neped
      !! ne0   : input : Central electron density [m-3]
      !! neped : input : Electron density at pedestal [m-3]
      use physics_variables, only: ne0, fne0, neped
      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      tmp_cc     = 1.0D0 - fne0 * ne0/neped
      tmp_con    = fne0
      tmp_err    = fne0 * tmp_cc
      tmp_symbol = '>'
      tmp_units  = '/m3'

   end subroutine constraint_eqn_081

   subroutine constraint_eqn_082(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for toroidal consistency of stellarator build
      !! author: J Lion, IPP Greifswald
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! toroidalgap > tftort
      !! #=# tfcoil
      !! #=#=# tftort, ftoroidalgap
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! ftoroidalgap : input real : f-value for constraint toroidalgap > tftort
      !! toroidalgap : input real :  minimal gap between two stellarator coils
      !! tftort : input real :  total toroidal width of a tf coil
      use tfcoil_variables, only: tftort,ftoroidalgap,toroidalgap
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - ftoroidalgap * toroidalgap/tftort
      tmp_con = toroidalgap
      tmp_err = toroidalgap - tftort/ftoroidalgap
      tmp_symbol = '<'
      tmp_units = 'm'

   end subroutine constraint_eqn_082

   subroutine constraint_eqn_083(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for radial consistency of stellarator build
      !! author: J Lion, IPP Greifswald
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! available_radial_space > required_radial_space
      !! #=# build
      !! #=#=# required_radial_space, f_avspace
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! f_avspace : input real : f-value for constraint available_radial_space > required_radial_space
      !! available_radial_space : input real :  avaible space in radial direction as given by each s.-configuration
      !! required_radial_space : input real :  required space in radial direction
      use build_variables, only: available_radial_space, required_radial_space, f_avspace
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - f_avspace  * available_radial_space/required_radial_space
      tmp_con = available_radial_space * (1.0D0 - tmp_cc)
      tmp_err = required_radial_space * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm'
   end subroutine constraint_eqn_083

   subroutine constraint_eqn_084(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for the lower limit of beta
      !! author: J Lion, IPP Greifswald
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !!  (beta-betaft) > betalim_lower
      !! #=# physics
      !! #=#=# betaft, beta, fbetatry_lower
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fbetatry_lower : input real : f-value for constraint beta-betaft > betalim_lower
      !! betalim_lower : input real :  Lower limit for beta
      !! beta : input real :  plasma beta
      !! betaft : input real : Alpha particle beta

      use physics_variables, only: betalim_lower, beta, betaft
      use constraint_variables, only: fbetatry_lower
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units


      tmp_cc = 1.0D0 - fbetatry_lower * (beta-betaft)/betalim_lower
      tmp_con = betalim_lower * (1.0D0 - tmp_cc)
      tmp_err = (beta-betaft) * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''


   end subroutine constraint_eqn_084

   subroutine constraint_eqn_085(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Author : S Kahn
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation constraining the centerpost (CP) lifetime
      !! Depending on the chosen option : i_cp_lifetime
      !!  - 0 : The CP full power year lifelime is set by the user (cplife_input)
      !!  - 1 : The CP lifelime is equal to the divertor one
      !!  - 2 : The CP lifetime is equal to the breeding blankets one
      !!  - 3 : The CP lifetime is equal to the plant one
      !! #=# availability
      !! #=#=# consistency
      !! Logic change during pre-factoring: err, symbol, units will be assigned
      !! only if present.
      !! cplife : input real : calculated CP full power year lifetime (years)
      !! bktlife : input real : calculated first wall/blanket power year lifetime (years)
      !! divlife : input real : calculated divertor  power year lifetime (years)
      !! i_cp_lifetime : input integer : switch chosing which plant element the CP
      !!                                 the CP lifetime must equate
      use cost_variables, only : cplife, divlife, cplife_input, &
         tlife, i_cp_lifetime
      use fwbs_variables, only : bktlife

      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      !! Constraints output


      ! The CP lifetime is equal to the the divertor one
      if  ( i_cp_lifetime == 0 ) then
         tmp_cc = 1.0D0 - cplife/cplife_input

      else if ( i_cp_lifetime == 1 ) then
         tmp_cc = 1.0D0 - cplife/divlife

      ! The CP lifetime is equal to the tritium breeding blankets / FW one
      else if ( i_cp_lifetime == 2 ) then
         tmp_cc = 1.0D0 - cplife/bktlife

      ! The CP lifetime is equal to the
      else if ( i_cp_lifetime == 3 ) then
         tmp_cc = 1.0D0 - cplife/tlife
      end if

      tmp_con = divlife * (1.0D0 - tmp_cc)
      tmp_err = divlife * tmp_cc
      tmp_symbol = '='
      tmp_units = 'years'

   end subroutine constraint_eqn_085

   subroutine constraint_eqn_086(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Author : S Kahn
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units;
      !!
      use tfcoil_variables, only : t_turn_tf, f_t_turn_tf, t_turn_tf_max

      implicit none

            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      !! Constraints output
      tmp_cc = 1.0D0 - t_turn_tf / ( f_t_turn_tf * t_turn_tf_max )
      tmp_con = t_turn_tf_max * (1.0D0 - tmp_cc)
      tmp_err = t_turn_tf_max * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'm'

   end subroutine constraint_eqn_086


   subroutine constraint_eqn_087(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! author: S. Kahn, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil cryogenic power upper limit

      use heat_transport_variables, only: crypmw, crypmw_max, f_crypmw
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - f_crypmw * crypmw_max/crypmw
      tmp_con = crypmw_max * (1.0D0 - tmp_cc)
      tmp_err = crypmw * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW'
   end subroutine constraint_eqn_087

   subroutine constraint_eqn_088(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for TF coil vertical strain upper limit (absolute value)
      !! author: CPS Swanson, PPPL, USA
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for TF coil vertical strain upper limit (absolute value)
      !! #=# tfcoil
      !! #=#=# fstr_wp, str_wp_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! fstr_wp : input real : f-value for TF coil strain
      !! str_wp_max : input real : Allowable maximum TF coil vertical strain
      !! str_wp : input real : Constrained TF coil vertical strain
      use constraint_variables, only: fstr_wp
      use tfcoil_variables, only: str_wp_max, str_wp
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc =  1.0D0 - fstr_wp * str_wp_max/abs(str_wp)
      tmp_con = str_wp_max
      tmp_err = str_wp_max - abs(str_wp) / fstr_wp
      tmp_symbol = '<'
      tmp_units = ''
   end subroutine constraint_eqn_088

   subroutine constraint_eqn_089(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Ensure that the Central Solenoid [OH] coil current / copper area < Maximum value
      !! author: G Turkington, CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! #=# physics
      !! #=#=# f_copperaoh_m2, copperaoh_m2_max
      !! and hence also optional here.
      !! Logic change during pre-factoring: err, symbol, units will be assigned only if present.
      !! copperaoh_m2 : input real : CS coil current at EOF / copper area [A/m2]
      !! copperaoh_m2_max : input real : maximum coil current / copper area [A/m2]
      !! f_copperaoh_m2 : input real : f-value for CS coil current / copper area
      use rebco_variables, only: copperaoh_m2, copperaoh_m2_max, f_copperaoh_m2
      implicit none
                  real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      tmp_cc = 1.0d0 - f_copperaoh_m2 * copperaoh_m2_max / copperaoh_m2
      tmp_con = copperaoh_m2
      tmp_err = copperaoh_m2 * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'A/m2'

   end subroutine constraint_eqn_089

   subroutine constraint_eqn_090(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! author: A. Pearce, G Turkington CCFE, Culham Science Centre
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !! Equation for minimum CS coil stress load cycles
      !! fncycle : input real : f-value for constraint n_cycle > n_cycle_min
      !! n_cycle : input real : Allowable number of cycles for CS
      !! n_cycle_min : input real :  Minimum required cycles for CS
      use CS_fatigue_variables, only: n_cycle, n_cycle_min, bkt_life_csf
      use constraint_variables, only: fncycle
      use cost_variables, only: ibkt_life, bktcycles
      implicit none
            real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units
      !! Switch to relay the calculated fw/blanket lifetime cycles as the minimum required CS stress cycles.
      !! bkt_life_cycle = 1 turns on the relay. Otherwise the models run independently.
      if (ibkt_life == 1 .and. bkt_life_csf == 1 ) then
         n_cycle_min = bktcycles
      end if

      tmp_cc =  1.0D0 - fncycle * n_cycle / n_cycle_min
      tmp_con = n_cycle_min * (1.0D0 - tmp_cc)
      tmp_err = n_cycle * tmp_cc
      tmp_symbol = '>'
      tmp_units = ''

   end subroutine constraint_eqn_090

   subroutine constraint_eqn_091(tmp_cc, tmp_con, tmp_err, tmp_symbol, tmp_units)
      !! Equation for checking if the design point is ECRH ignitable
      !! at lower values for n and B. Or if the design point is ECRH heatable (if ignite==0)
      !! stellarators only (but in principle usable also for tokamaks).
      !! author: J Lion, IPP Greifswald
      !! args : output structure : residual error; constraint value;
      !! residual error in physical units; output string; units string
      !!  powerht_local > powerscaling
      !! #=# physics
      !! #=#=# fecrh_ignition, powerht_local, powerscaling
      !! fecrh_ignition : input real : f-value for constraint powerht_local > powerscaling
      !! max_gyrotron_frequency : input real :  Max. av. gyrotron frequency
      !! te0_ecrh_achievable : input real : Max. achievable electron temperature at ignition point
      use constraint_variables, only: fecrh_ignition
      use stellarator_variables, only: max_gyrotron_frequency, te0_ecrh_achievable, powerscaling_constraint, powerht_constraint
      use physics_variables, only: ignite
      use current_drive_variables, only: pheat
      implicit none
      real(dp), intent(out) :: tmp_cc
      real(dp), intent(out) :: tmp_con
      real(dp), intent(out) :: tmp_err
      character(len=1), intent(out) :: tmp_symbol
      character(len=10), intent(out) :: tmp_units

      ! Achievable ECRH te needs to be larger than needed te for igntion
      if(ignite==0) then
         tmp_cc = 1.0D0 - fecrh_ignition* (powerht_constraint+pheat)/powerscaling_constraint
      else
         tmp_cc = 1.0D0 - fecrh_ignition* powerht_constraint/powerscaling_constraint
      endif

      tmp_con = powerscaling_constraint * (1.0D0 - tmp_cc)
      tmp_err = powerht_constraint * tmp_cc
      tmp_symbol = '<'
      tmp_units = 'MW'
   end subroutine constraint_eqn_091



end module constraints
