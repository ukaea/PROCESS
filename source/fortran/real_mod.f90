module real_mod

  interface real8

     function real8_real_4_(x)
     use precision_mod, only: real_4_, p_
       implicit none
       real(p_) :: real8_real_4_
       real(real_4_), intent(in) :: x
       !return real(x, real_8_)
     end function real8_real_4_

     function real8_real_8_(x)
       use precision_mod, only: real_8_, p_
       implicit none
       real(p_) :: real8_real_8_
       real(real_8_), intent(in) :: x
       !return real(x, real_8_)
     end function real8_real_8_

     function real8_i(x)
       use precision_mod, only: p_
       implicit none
       real(p_) :: real8_i
       integer, intent(in) :: x
       !return real(x, real_8_)
     end function real8_i

     function real8_c(x)
       use precision_mod, only: real_4_, p_
       implicit none
       real(p_) :: real8_c
       complex(real_4_), intent(in) :: x
       !return real(x, real_8_)
     end function real8_c

     function real8_c8(x)
       use precision_mod, only: real_8_, p_
       implicit none
       real(p_) :: real8_c8
       complex(real_8_), intent(in) :: x
       !return real(x, real_8_)
     end function real8_c8

  end interface

  interface cmplx8

     function cmplx8_real_4_(x, y)
       use precision_mod, only: real_4_, p_
       implicit none
       complex(p_) :: cmplx8_real_4_
       real(real_4_), intent(in) :: x, y
       !return cmplx(x, y, real_8_)
     end function cmplx8_real_4_

     function cmplx8_real_8_(x, y)
       use precision_mod, only: real_8_, p_
       implicit none
       complex(p_) :: cmplx8_real_8_
       real(real_8_), intent(in) :: x, y
       !return cmplx(x, y, real_8_)
     end function cmplx8_real_8_

     function cmplx8_i(x, y)
       use precision_mod, only: p_
       implicit none
       complex(p_) :: cmplx8_i
       integer, intent(in) :: x, y
       !return cmplx(x, y, real_8_)
     end function cmplx8_i

  end interface

end module real_mod
