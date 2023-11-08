! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module current_drive_module
  !! Module containing current drive system routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains routines relevant for calculating the
  !! current drive parameters for a fusion power plant.
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Import modules
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

contains
  function sigbeam(eb,te,ne,rnhe,rnc,rno,rnfe)

    !! Calculates the stopping cross-section for a hydrogen
    !! beam in a fusion plasma
    !! author: P J Knight, CCFE, Culham Science Centre
    !! eb     : input real : beam energy (kev/amu)
    !! te     : input real : electron temperature (keV)
    !! ne     : input real : electron density (10^20m-3)
    !! rnhe   : input real : alpha density / ne
    !! rnc    : input real : carbon density /ne
    !! rno    : input real : oxygen density /ne
    !! rnfe   : input real : iron density /ne
    !! This function calculates the stopping cross-section (m^-2)
    !! for a hydrogen beam in a fusion plasma.
    !! Janev, Boley and Post, Nuclear Fusion 29 (1989) 2125
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(dp) :: sigbeam

    !  Arguments

    real(dp), intent(in) :: eb,te,ne,rnhe,rnc,rno,rnfe

    !  Local variables

    real(dp) :: ans,nen,sz,s1
    real(dp), dimension(2,3,2) :: a
    real(dp), dimension(3,2,2,4) :: b
    real(dp), dimension(4) :: nn,z

    integer :: i,is,j,k

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    data a/ 4.40D0, 2.30D-1, 7.46D-2,-2.55D-3, 3.16D-3, 1.32D-3, &
         -2.49D-2,-1.15D-2, 2.27D-3,-6.20D-4,-2.78D-5, 3.38D-5/

    data b/ &
         -2.36D0, 8.49D-1,-5.88D-2,-2.50D-1, 6.77D-2,-4.48D-3, &
         1.85D-1,-4.78D-2, 4.34D-3,-3.81D-2, 1.05D-2,-6.76D-4, &
         -1.49D0, 5.18D-1,-3.36D-2,-1.19D-1, 2.92D-2,-1.79D-3, &
         -1.54D-2, 7.18D-3, 3.41D-4,-1.50D-2, 3.66D-3,-2.04D-4, &
         -1.41D0, 4.77D-1,-3.05D-2,-1.08D-1, 2.59D-2,-1.57D-3, &
         -4.08D-4, 1.57D-3, 7.35D-4,-1.38D-2, 3.33D-3,-1.86D-4, &
         -1.03D0, 3.22D-1,-1.87D-2,-5.58D-2, 1.24D-2,-7.43D-4, &
         1.06D-1,-3.75D-2, 3.53D-3,-3.72D-3, 8.61D-4,-5.12D-5/

    z(1) = 2.0D0 ; z(2) = 6.0D0 ; z(3) = 8.0D0 ; z(4) = 26.0D0
    nn(1) = rnhe ; nn(2) = rnc ; nn(3) = rno ; nn(4) = rnfe

    nen = ne/1.0D19

    s1 = 0.0D0
    do k = 1,2
       do j = 1,3
          do i = 1,2
             s1 = s1 + a(i,j,k) * (log(eb))**(i-1) * &
                  (log(nen))**(j-1) * (log(te))**(k-1)
          end do
       end do
    end do

    !  Impurity term

    sz = 0.0D0
    do is = 1,4
       do k = 1,2
          do j = 1,2
             do i = 1,3
                sz = sz + b(i,j,k,is)* (log(eb))**(i-1) * &
                     (log(nen))**(j-1) * (log(te))**(k-1) * &
                     nn(is) * z(is) * (z(is)-1.0D0)
             end do
          end do
       end do
    end do

    ans = 1.0D-20 * ( exp(s1)/eb * (1.0D0 + sz) )

    sigbeam = max(ans,1.0D-23)

  end function sigbeam

end module current_drive_module
