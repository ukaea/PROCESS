module mod_f90_kind

  implicit none

  ! Single byte integer
  integer, parameter :: int_def = selected_int_kind(9)

  ! Long byte integer
  integer, parameter :: int_long = selected_int_kind(18)

  ! Single precision
  integer, parameter :: single  = selected_real_kind(6,37)

  ! Double precision
  integer, parameter :: double  = selected_real_kind(15,300)

  ! logical
  integer, parameter :: logical_def = kind(.true.)

  !-----------------------------------------------------------------------
  ! Working precision
  !-----------------------------------------------------------------------
  integer, parameter :: rkind   = double      ! double
  integer, parameter :: r4kind  = single      ! single
  integer, parameter :: ikind   = int_def     ! integer
  integer, parameter :: i8kind   = int_long   ! integer double
  integer, parameter :: lkind   = logical_def ! integer
  integer, parameter :: skind   = 512         ! string

end module mod_f90_kind

module param
  use mod_f90_kind, only: rkind   !EXTERNAL
  implicit none
  real(rkind), parameter :: pi = 3.14159265358979e0_rkind, &
                            c  = 2.99792458e10_rkind
end module param

module freq
  use mod_f90_kind, only: rkind   !EXTERNAL
  implicit none
  real(rkind), save :: om, xk0
end module freq

module mode
  use mod_f90_kind, only: ikind   !EXTERNAL
  implicit none
  integer(ikind), save :: imod
end module mode

module machin
  use mod_f90_kind, only: rkind   !EXTERNAL
  implicit none
  real(rkind), save :: rmaj, rmin, sgnm, zsearch
  character(len=160), save :: machname
end module machin

module bsquar
  use mod_f90_kind, only: rkind   !EXTERNAL
  implicit none
  real(rkind), save :: conf, conf2, qb2
end module bsquar

module linliu
  use mod_f90_kind, only: rkind   !EXTERNAL
  implicit none
  save
  complex(rkind) :: PCEXNORMO, PCEYNORMO, PCEZNORMO
  complex(rkind) :: PCEXNORMX, PCEYNORMX, PCEZNORMX
end module linliu

module abs_cd
  use mod_f90_kind, only: ikind   !EXTERNAL
  implicit none
  integer(ikind), save :: cdroutine, absroutine, iastra, profcalc
end module abs_cd
