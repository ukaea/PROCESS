module precision_mod

  integer, parameter :: real_8_ = selected_real_kind(12,100)
  integer, parameter :: real_4_ = selected_real_kind(6,36)

  ! Precision can be adjusted below:

  integer, parameter :: p_ = real_8_ ! real_8_

end module precision_mod
