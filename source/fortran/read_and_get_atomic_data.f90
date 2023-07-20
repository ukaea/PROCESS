
module read_and_get_atomic_data
  !! Module for reading atomic data
  !! author: M Kovari, CCFE, Culham Science Centre
  !! N/A
  !! This module reads atomic data for the PROCESS Kallenbach divertor model
  !! 
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Modules to import !
  ! !!!!!!!!!!!!!!!!!!!!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  ! Var for subroutine get_h_rates requiring re-initialisation on each new run
  logical, save :: FirstCall

contains
  character(len=300) function hdatadir()
        implicit none
        character(len=200) :: process_dir
        CALL get_environment_variable("PYTHON_PROCESS_ROOT", process_dir)
        if (process_dir == "") then
          hdatadir = INSTALLDIR//'/process/data/h_data/'
        else
          hdatadir = trim(process_dir)//'/data/h_data/' 
        end if
        
  end function hdatadir

  subroutine init_read_and_get_atomic_data
    !! Initialise module variables
    implicit none

    FirstCall = .true.
  end subroutine init_read_and_get_atomic_data

  subroutine get_h_rates(density, temperature, s, al, Rcx, plt, prb, mass, verbose)
    !! 
    !! author: M Kovari, CCFE, Culham Science Centre
    !! density : input real : electron density  (m^-3)
    !! temperature : input real : electron temperature  (eV)
    !! s : output real : ionisation rate coefficient (m^3/s)
    !! al : output real : recombination rate coefficient (m^3/s)
    !! Rcx : output real : charge exchange rate coefficient (m^3/s)
    !! plt : output real : line radiation power rate coefficient (Wm^3)
    !! prb : output real : continuum radiation power rate coefficient (Wm^3)
    !! mass : input real : relative atomic mass for CX rate coefficient
    !! verbose : input logical : verbose switch
    !! 
    !! 
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use maths_library, only: interpolate
    implicit none

    ! Subroutine declarations !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!

    real(dp), intent(in) :: density, temperature, mass

    real(dp), intent(out) :: s, al, Rcx, plt, prb

    logical, intent(in) :: verbose

    integer :: m

    ! These arrays are read only once and then saved
    ! Density values: "_d"
    real(dp), save, dimension(24) :: scd_d, acd_d, ccd_d, plt_d,prb_d

    ! Temperature values: "_t"
    real(dp), save, dimension(29) :: scd_t,acd_t,ccd_t, plt_t, prb_t

    ! Rate coefficients: "_r"
    real(dp), save, dimension(24,29) :: scd_r, acd_r, ccd_r, plt_r,prb_r

    character(len=200) :: acd_file, scd_file, plt_file, prb_file, ccd_file

    real(dp) :: logdens, logtemp

    logical :: iexist

    integer :: ine, ite

    ! Maxima for log density and log temperature in each data file
    real(dp), save :: max_scd_d, max_scd_t
    real(dp), save :: max_acd_d, max_acd_t
    real(dp), save :: max_ccd_d, max_ccd_t
    real(dp), save :: max_plt_d, max_plt_t
    real(dp), save :: max_prb_d, max_prb_t

    ine = 24
    ite = 29

    logdens = log10(density/1.0D6)
    logtemp = log10(temperature)

    ! Read the data
    if(FirstCall) then

      FirstCall=.false.

      !  Add trailing / to hdatadir() if necessary
      ! if (index(hdatadir(),'/',.true.) .ne. len(trim(hdatadir()))) hdatadir() = hdatadir()//'/'
      ! if (index(hdatadir(),'\',.true.) .ne. len(trim(hdatadir()))) hdatadir() = hdatadir()//'\'

      acd_file = trim(hdatadir())//'acd12_h.dat'
      scd_file = trim(hdatadir())//'scd12_h.dat'
      plt_file = trim(hdatadir())//'plt12_h.dat'
      prb_file = trim(hdatadir())//'prb12_h.dat'

      m = floor(mass+0.5)                     ! round to an integer
      ! Select the correct atomic species for the CX rates: H, D or T
      if      (m == 1) then
          ccd_file = trim(hdatadir())//'ccd96_h.dat'
      elseif (m == 2) then
          ccd_file = trim(hdatadir())//'ccd96_d.dat'
      elseif (m == 3) then
          ccd_file = trim(hdatadir())//'ccd96_t.dat'
      else
          write(*,*) 'The atomic mass is ', m, '.  It must be in the range 1-3.'
      end if

      inquire(file=trim(acd_file), exist=iexist)
      if (.not.iexist) write(*,*) "ERROR  File "// acd_file // ' does not exist'
      inquire(file=trim(scd_file), exist=iexist)
      if (.not.iexist) write(*,*) "ERROR  File "// scd_file // ' does not exist'
      inquire(file=trim(plt_file), exist=iexist)
      if (.not.iexist) write(*,*) "ERROR  File "// plt_file // ' does not exist'
      inquire(file=trim(prb_file), exist=iexist)
      if (.not.iexist) write(*,*) "ERROR  File "// prb_file // ' does not exist'
      inquire(file=trim(ccd_file), exist=iexist)
      if (.not.iexist) write(*,*) "ERROR  File "// ccd_file // ' does not exist'

      ! ionisation data
      call read_atomdat(scd_d,scd_t,scd_r, ine=24, ite=29, filename=scd_file, verbose=verbose)
      ! recombination data
      call read_atomdat(acd_d,acd_t,acd_r, ine=24, ite=29, filename=acd_file, verbose=verbose)
      ! CX data
      call read_atomdat(ccd_d,ccd_t,ccd_r, ine=24, ite=29, filename=ccd_file, verbose=verbose)
      ! line radiation data
      call read_atomdat(plt_d,plt_t,plt_r, ine=24, ite=29, filename=plt_file, verbose=verbose)
      ! continuum radiation data
      call read_atomdat(prb_d,prb_t,prb_r, ine=24, ite=29, filename=prb_file, verbose=verbose)

      ! Store the maxima for log density and log temperature in each data file
      ! Subtract a smidgen to ensure that the values submitted for interpolation are
      ! the range.
      max_scd_d=maxval(scd_d) - 0.00001
      max_scd_t=maxval(scd_t) - 0.00001

      max_acd_d=maxval(acd_d) - 0.00001
      max_acd_t=maxval(acd_t) - 0.00001

      max_ccd_d=maxval(ccd_d) - 0.00001
      max_ccd_t=maxval(ccd_t) - 0.00001

      max_plt_d=maxval(plt_d) - 0.00001
      max_plt_t=maxval(plt_t) - 0.00001

      max_prb_d=maxval(prb_d) - 0.00001
      max_prb_t=maxval(prb_t) - 0.00001

    end if

    ! Using function interpolate(x_len, x_array, y_len, y_array, f, x, y, delta)
    ! ionisation
    s   = 1.d-6*10.0**(interpolate(ine, scd_d, ite, scd_t, scd_r,   &
                                  min(logdens,max_scd_d),           &
                                  min(logtemp,max_scd_t)))
    ! recombination
    al  = 1.d-6*10.0**(interpolate(ine, acd_d, ite, acd_t, acd_r,   &
                                  min(logdens,max_acd_d),           &
                                  min(logtemp,max_acd_t)))
    ! Rcx
    Rcx  = 1.d-6*10.0**(interpolate(ine, ccd_d, ite, ccd_t, ccd_r,   &
                                  min(logdens,max_ccd_d),           &
                                  min(logtemp,max_ccd_t)))
    ! line radiation
    plt = 1.d-6*10.0**(interpolate(ine, plt_d, ite, plt_t, plt_r,   &
                                  min(logdens,max_plt_d),           &
                                  min(logtemp,max_plt_t)))
    ! continuum radiation
    prb = 1.d-6*10.0**(interpolate(ine, prb_d, ite, prb_t, prb_r,   &
                                  min(logdens,max_prb_d),           &
                                  min(logtemp,max_prb_t)))

  end subroutine get_h_rates

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_atomdat(density,temperature,rates, ine, ite, filename, verbose)
    !! 
    !! author: M Kovari, CCFE, Culham Science Centre
    !! density : input real : electron density  (m^-3)
    !! temperature : input real : electron temperature  (eV)
    !! s : input  : ionisation rate coefficient (m^3/s)
    !! al : input  : recombination rate coefficient (m^3/s)
    !! Rcx : input  : charge exchange rate coefficient (m^3/s)
    !! plt : input  : line radiation power rate coefficient (Wm^3)
    !! prb : input  : continuum radiation power rate coefficient (Wm^3)
    !! mass : input  : relative atomic mass for CX rate coefficient
    !! verbose : output  : verbose switch
    !! 
    !! 
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    ! Subroutine declarations !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!

    real(dp), dimension(ine), intent(out) :: density
    real(dp), dimension(ite), intent(out) :: temperature
    real(dp), dimension(ine,ite), intent(out) :: rates
    character(len=200), intent(in) :: filename
    logical, intent(in), optional::verbose
    integer, intent(in)::ine,ite
    integer::i, l

    if (verbose) write(*,*) filename

    !  Open the input/output external files
    open(unit=8,file=filename,status='old')
    read(8,*)
    read(8,*)
    read(8,'(8(f10.5))')(density(i), i= 1, ine)
    read(8,'(8(f10.5))')(temperature(i), i= 1, ite)

    read(8,*)
    read(8,'(8(f10.5))')((rates(I,L),I=1,ine),L=1,ite)
    close(unit=8)

  end subroutine read_atomdat

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine unit_test_read()
        ! Radiative cooling function Lz
        ! To test the interpolation, use a point at the geometrical mean of the two
        ! first temperatures and the two first values of ne.tau
        use read_radiation, only: read_lz
        implicit none

        real(dp):: s, al, Rcx, plt, prb, density, temperature, mass
        real(dp):: te,netau,test_lz,estimate_lz
        
        te=sqrt(1.000 * 1.047)
        netau=sqrt(0.1*1.0)

        test_lz= read_lz('Ar', te, netau, mean_z=.false., mean_qz=.false., verbose=.true.)
        write(*,*) 'Compare Lz values from code with values interpolated by hand:'
        write(*,*)'test_lz = ',test_lz
        estimate_lz = sqrt(sqrt(1.772e-37))*sqrt(sqrt(2.916e-37))*sqrt(sqrt(1.772e-37))*sqrt(sqrt(2.916e-37))
        write(*,*)'estimate_lz = ',estimate_lz

        density = 10.0**((7.69897 + 8.00000)/2.0 +6.0)
        temperature = 10.0**((-0.69897 -0.52288)/2.0 )
        mass=2.0
        write(*,'(5(e12.5),2x)')density, temperature, mass

        call get_h_rates(density, temperature, &
                          s, al, Rcx, plt, prb, &
                          mass, verbose=.true.)

        write(*,*) 'Compare atomic rate values from code with values interpolated by hand:'
        write(*,'(5(e12.5),2x)')s, al, Rcx, plt, prb

        s  = 10.0**((-37.87042 -37.83901  -27.99641 -27.97238)/4.0) / 1.0e6
        al = 10.0**((-11.85004 -11.83456  -11.99890 -11.98826)/4.0) / 1.0e6
        Rcx = 10.0**((-8.59260  -8.59260    -8.50455  -8.50455)/4.0) / 1.0e6
        plt= 10.0**((-47.19666d0 -47.19666d0  -39.88720d0 -39.88720d0)/4.0d0) / 1.0d6
        prb= 10.0**((-29.50732 -29.49208  -29.65257 -29.64218)/4.0) / 1.0e6
        write(*,'(5(e12.5),2x)')s, al, Rcx, plt, prb
        return
  end subroutine unit_test_read

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine plot_rates()
    ! Reads rate coefficients for deuterium.
    ! Compare to Figure 2 in Kallenbach 2016.
    real(dp):: s(3), al(3), Rcx
    real(dp):: plt(3), prb(3), mass
    real(dp):: lz_deuterium(3)
    real(dp):: dummy1, dummy2, dummy3, dummy4, dummy5
    integer::i,j
    real(dp), parameter ::te(15)=(/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,12.,14.,16.,18.,20./)
    real(dp), parameter ::density(3)=(/1.e19,1.e20,1.e21/)

    open(unit=12,file='rate_coefficients.txt',status='replace')
    write(12,'(30a11)')'te [eV]','Rcx', 'ionis19', 'recomb19', 'line rad19', 'cont rad19', 'tot rad19',&
                                       'ionis20', 'recomb20', 'line rad20', 'cont rad20', 'tot rad20',&
                                       'ionis21', 'recomb21', 'line rad21', 'cont rad21', 'tot rad21'
    mass=2.0D0
    ! Just read data
    call get_h_rates(1.d20, 1.0d0, dummy1, dummy2, dummy3, dummy4, dummy5, &
                          mass, verbose=.true.)
    do i=1,15
        do j=1,3
            call get_h_rates(density(j), te(i), s(j), al(j), Rcx, plt(j), prb(j), &
                             mass, verbose=.false.)
            lz_deuterium(j)=plt(j)+prb(j)
        enddo
        write(12,'(30es11.3)')te(i), Rcx, (s(j), al(j), plt(j), prb(j), lz_deuterium(j),  j=1,3)
    enddo
    close(unit=12)

  end subroutine plot_rates


end module read_and_get_atomic_data
