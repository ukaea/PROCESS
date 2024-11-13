!#######################################################################

 MODULE const_and_precisions

!########################################################################!
 IMPLICIT NONE
 PUBLIC
!------------------------------------------------------------------------
! common precisions
!------------------------------------------------------------------------
 INTEGER, PARAMETER :: sp_   = 4           ! single precision
 INTEGER, PARAMETER :: dp_   = 8           ! double precision
! INTEGER, PARAMETER :: sp_   = KIND(1.0E0) ! single precision
! INTEGER, PARAMETER :: dp_   = KIND(1.0D0) ! double precision
 INTEGER, PARAMETER :: wp_   = dp_         ! work-precision
 INTEGER, PARAMETER :: odep_ = dp_         ! ODE-solver precision
!------------------------------------------------------------------------
! precisions which are in use in CONFIG_yat
!------------------------------------------------------------------------
 INTEGER, PARAMETER :: ypi_ = 4          ! <- direct precision def.
 INTEGER, PARAMETER :: ypd_ = 8          ! <- direct precision def.
!------------------------------------------------------------------------
! length of the file names
!------------------------------------------------------------------------
 INTEGER, PARAMETER :: lfn_ = 256        ! <- requested for yat
!========================================================================
! Arithmetic constants
!========================================================================
 REAL(wp_), PARAMETER :: zero    = 0.0_wp_
 REAL(wp_), PARAMETER :: unit    = 1.0_wp_
 REAL(wp_), PARAMETER :: pi      = 3.141592653589793_wp_
 REAL(wp_), PARAMETER :: sqrt_pi = 1.772453850905516_wp_
 REAL(wp_), PARAMETER :: sqrt_2  = 1.414213562373095_wp_
 REAL(wp_), PARAMETER :: rad     = pi/180.0_wp_
!---
 REAL(wp_), PARAMETER :: ex(1:3) = (/unit,zero,zero/)
 REAL(wp_), PARAMETER :: ey(1:3) = (/zero,unit,zero/)
 REAL(wp_), PARAMETER :: ez(1:3) = (/zero,zero,unit/)
!---
 REAL(wp_), PARAMETER :: kron(3,3) = reshape((/unit,zero,zero, &
                                               zero,unit,zero, &
                                               zero,zero,unit/),(/3,3/))
 COMPLEX(wp_), PARAMETER :: im    = (0.0_wp_,1.0_wp_)
 COMPLEX(wp_), PARAMETER :: czero = (0.0_wp_,0.0_wp_)
 COMPLEX(wp_), PARAMETER :: cunit = (1.0_wp_,0.0_wp_)
 COMPLEX(wp_), PARAMETER :: ctwo  = (2.0_wp_,0.0_wp_)
!========================================================================
! Computer constants
!========================================================================
 REAL(wp_), PARAMETER :: comp_eps      = EPSILON(unit)
 REAL(wp_), PARAMETER :: comp_eps2     = comp_eps**2
 REAL(wp_), PARAMETER :: comp_tiny     = TINY(unit)
 REAL(wp_), PARAMETER :: comp_huge     = HUGE(unit)
 REAL(wp_), PARAMETER :: comp_tinylog  =-200 ! LOG10(comp_tiny)
 REAL(wp_), PARAMETER :: comp_hugelog  =+200 ! LOG10(comp_huge)
! REAL(wp_), PARAMETER :: comp_tiny1    = 1d+50*comp_tiny
! REAL(wp_), PARAMETER :: comp_huge1    = 1d-50*comp_huge
! REAL(wp_), PARAMETER :: comp_tiny1log = LOG10(comp_tiny1)
! REAL(wp_), PARAMETER :: comp_huge1log = LOG10(comp_huge1)
!------------------------------------------------------------------------
! Conventional constants
!------------------------------------------------------------------------
 REAL(wp_), PARAMETER :: output_tiny = 1.0d-66
 REAL(wp_), PARAMETER :: output_huge = 1.0d+66
!========================================================================
! Physical constants (SI)
!========================================================================
 REAL(wp_), PARAMETER :: e_     = 1.601917d-19    ! [C]
 REAL(wp_), PARAMETER :: me_    = 9.109558d-31    ! [kg]
 REAL(wp_), PARAMETER :: mp_    = 1.672614d-27    ! [kg]
 REAL(wp_), PARAMETER :: rmpe_  = mp_/me_
 REAL(wp_), PARAMETER :: c_     = 2.997925d+08    ! [m/s]
 REAL(wp_), PARAMETER :: eps0_  = 8.854188d-12    ! [F/m]
!------------------------------------------------------------------------
! Useful definitions
!------------------------------------------------------------------------
 REAL(wp_), PARAMETER :: keV_   = 1000*e_         ! [J]
 REAL(wp_), PARAMETER :: mc2_SI = me_*c_**2       ! [J]
 REAL(wp_), PARAMETER :: mc2_   = mc2_SI/keV_     ! [keV]
 REAL(wp_), PARAMETER :: mc_    = me_*c_          ! [kg*m/s]
 ! f_ce = fce1_*B (B in Tesla):                   !
 REAL(wp_), PARAMETER :: wce1_  = e_/me_          ! [rad/s]
 REAL(wp_), PARAMETER :: fce1_  = wce1_/(2*pi)    ! [1/s]
 ! f_pl = fpe1_*sqrt(Ne) (Ne in 1/m**3):          !
 REAL(wp_), PARAMETER :: wpe1_  = 56.4049201      ! [rad/s]
 REAL(wp_), PARAMETER :: fpe1_  = wpe1_/(2*pi)    ! [1/s]
 REAL(wp_), PARAMETER :: wpe12_ = wpe1_**2        !
 ! vte  = vte1_*sqrt(Te) (Te in keV):             !
 REAL(wp_), PARAMETER :: vte1_  = 1.8755328e7     ! [m/s]
 ! je   = curr1_*sqrt(Te)*Ne  (Ne in 1/m**3):     !
 REAL(wp_), PARAMETER :: curr1_ = e_*vte1_        ! [A/m**2]
!========================================================================
! Upper limit for the momentum value for integration
!========================================================================
 REAL(wp_), PARAMETER :: umax_ = 7.0d0            ! max of (p/pth)
 INTEGER,   PARAMETER :: nu_   = 700              ! size of upar-array
!========================================================================
! minimal value of Nparallel
!========================================================================
 REAL(wp_), PARAMETER :: Npar_min = 1.0d-2
!########################################################################!

 END MODULE const_and_precisions

!########################################################################!
