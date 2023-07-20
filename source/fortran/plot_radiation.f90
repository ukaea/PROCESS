module plot_radiation
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
implicit none

contains
    subroutine plot_Lz()
    !! Write loss data to file for plotting
    !! author: M Kovari, CCFE, Culham Science Centre
    !! Write loss data to file for plotting
    !! Compare to Figure 3 in Kallenbach 2016.
    !!
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use impurity_radiation_module, only: nimp, imp_label
    use read_radiation, only: read_lz
    implicit none

    ! Subroutine declarations !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=2) :: element

    real(dp) :: dummy

    integer :: i, j

    integer, parameter :: points = 27

    ! Temperature plot points
    real(dp), parameter :: te(points) = (/1., 1.2, 1.5, 2., 2.5, 3., 4., 5., 6., 7., 8., 9., &
        10., 12., 14., 16., 20., 30., 40., 50., 60., 70., 80., 90., 100., 150., 200./)

    real(dp) :: Lz_plot(nimp)

    real(dp), parameter :: netau = 0.5

    open(unit=12,file='radiative_loss_functions.txt',status='replace')
    write(12,'(30a11)')'Te (eV)', (imp_label(i), i=2,nimp)

    ! Just read data.  Exclude hydrogen by starting at 2
    do i=2,nimp
        element=imp_label(i)
        dummy=read_lz(element,30.0d0,netau, mean_z=.false., mean_qz=.false., verbose=.false.)
    enddo

    do i=1,points
        do j=2,nimp
            Lz_plot(j)=read_lz(imp_label(j),te(i),netau, mean_z=.false., mean_qz=.false., verbose=.false.)
        enddo
        write(12,'(30es11.3)')te(i), (Lz_plot(j), j=2,nimp)
    enddo
    close(unit=12)

  end subroutine plot_Lz

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine plot_z()
    !! Write z and z^2 data to file for plotting
    !! author: M Kovari, CCFE, Culham Science Centre
    !! Write z and z^2 data to file for plotting
    !!
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use impurity_radiation_module, only: nimp, imp_label
    use read_radiation, only: read_lz
    implicit none

    ! Subroutine declarations !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!

    character(len=2) :: element

    real(dp) :: dummy

    integer :: i, j, k

    integer, parameter :: points = 27

    real(dp), parameter :: te(points) = (/1., 1.2, 1.5, 2., 2.5, 3.,4., 5., 6., 7., &
        8., 9.,10.,12., 14., 16., 20., 30., 40., 50., 60., 70., 80., 90., 100., 150., 200./)

    real(dp) :: Z_plot(3,nimp)

    real(dp), parameter :: netau(3) = (/0.1, 1.0, 10.0/)

    open(unit=12,file='mean_Z.txt',status='replace')

    write(12,*)'Mean Z'

    write(12,'(a11,  42(3(a4, es7.1)))')'Te (ev)',((imp_label(i),netau(j), j=1,3),i=2,nimp)

    ! Just read data.  Exclude hydrogen by starting at 2
    do i = 2, nimp
        element = imp_label(i)
        dummy = read_lz(element,30.0d0,1.0d0, mean_z=.true., mean_qz=.false., verbose=.true.)
        dummy = read_lz(element,30.0d0,1.0d0, mean_z=.false., mean_qz=.true., verbose=.true.)
    enddo

    do i = 1, points
        do j = 1, 3
            do k = 2, nimp
                element = imp_label(k)
                Z_plot(j,k) = read_lz(element,te(i),netau(j), mean_z=.true., mean_qz=.false., verbose=.false.)
            enddo
        enddo
        write(12,'(42es11.3)')te(i), ((Z_plot(j,k), j=1,3), k=2,nimp)
    enddo

    write(12,*)
    write(12,*)'Mean Z^2'
    write(12,'(a11,  42(3(a4, es7.1)))')'Te (ev)', ((imp_label(i),netau(j), j=1,3),i=2,nimp)
    do i = 1, points
        do j = 1, 3
            do k = 2, nimp
                element = imp_label(k)
                Z_plot(j,k) = read_lz(element,te(i),netau(j), mean_z=.false., mean_qz=.true., verbose=.false.)
            enddo
        enddo
        write(12,'(42es11.3)')te(i), ((Z_plot(j,k), j=1,3), k=2,nimp)
    enddo

    close(unit=12)

  end subroutine plot_z

end module plot_radiation
