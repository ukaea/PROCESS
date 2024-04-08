! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module process_output

  !! Module containing routines to produce a uniform output style
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains a number of routines that allow the
  !! program to write output to a file unit in a uniform style.
  !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  public



contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ocentr(file,string,width)

    !! Routine to print a centred header within a line of asterisks
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! string : input character string : Character string to be used
    !! width : input integer : Total width of header
    !! This routine writes out a centred header within a line of asterisks.
    !! It cannot cope with a zero-length string; routine
    !! <A HREF="ostars.html"><CODE>ostars</CODE></A> should be used instead.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: active_constraints, ncalls, ipnvars, ioptimz
    use global_variables, only: run_tests, verbose, output_prefix
		use constants, only: mfile
    use maths_library, only: secant_solve
    implicit none

    !  Arguments

    integer, intent(in) :: file, width
    character(len=*), intent(in) :: string

    !  Local variables

    integer :: lh, nstars, nstars2
    integer, parameter :: maxwidth = 110
    character(len=maxwidth) :: stars

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    stars = repeat('*',maxwidth)

    lh = len(string)

    if (lh == 0) then
       call ostars(file,width)
       return
    end if

    if (width > maxwidth) then
       write(*,*) 'Error in routine OCENTR :'
       write(*,*) 'Maximum width = ',maxwidth
       write(*,*) 'Requested width = ',width
       write(*,*) 'PROCESS stopping.'
       stop 1
    end if

    if (lh >= width) then
       write(*,*) 'Error in routine OCENTR :'
       write(*,*) string
       write(*,*) 'This is too long to fit into ',width,' columns.'
       write(*,*) 'PROCESS stopping.'
       stop 1
    end if

    !  Number of stars to be printed on the left

    nstars = int( (width-lh)/2 ) - 1

    !  Number of stars to be printed on the right

    nstars2 = width - (nstars+lh+2)

    !  Write the whole line

    write(file,'(t2,a)') stars(1:nstars)//' '//string//' '//stars(1:nstars2)

    write(mfile,'(t2,a)') '#'//' '//string//' '//'#'

  end subroutine ocentr

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ostars(file,width)

    !! Routine to print a line of asterisks
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! width : input integer : Total width of header
    !! This routine writes out a line of asterisks.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use global_variables, only: output_prefix, fileprefix
		use constants, only: mfile
    implicit none

    !  Arguments

    integer, intent(in) :: file, width

    !  Local variables

    integer, parameter :: maxwidth = 110
    character(len=maxwidth) :: stars

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    stars = repeat('*',maxwidth)

    write(file,'(1x,a)') stars(1:min(width,maxwidth))

  end subroutine ostars

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine oheadr(file,string)

    !! Routine to print a centred header within a line of asterisks,
    !! and between two blank lines
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! string : input character string : Character string to be used
    !! This routine writes out a centred header within a line of
    !! asterisks, and between two blank lines.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: active_constraints, sqsumsq, ioptimz
    use global_variables, only: vlabel, run_tests, verbose
		use constants, only: rmu0, pi
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: string

    !  Local variables

    integer, parameter :: width = 110

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call oblnkl(file)
    call ocentr(file,string,width)
    call oblnkl(file)

  end subroutine oheadr

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine oshead(file,string)

    !! Routine to print a short, centred header within a line of asterisks,
    !! and between two blank lines
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! string : input character string : Character string to be used
    !! This routine writes out a short, centred header within a line of
    !! asterisks, and between two blank lines.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use global_variables, only: icase
		use constants, only: pi
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: string

    !  Local variables

    integer, parameter :: width = 80

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call oblnkl(file)
    call ocentr(file,string,width)
    call oblnkl(file)

  end subroutine oshead

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine oblnkl(file)

    !! Routine to print a blank line
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! This routine writes out a simple blank line.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use constants, only: pi, nout
    implicit none

    !  Arguments

    integer, intent(in) :: file

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write(file,10)
10  format(' ')

  end subroutine oblnkl

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine osubhd(file,string)

    !! Routine to print a subheading between two blank lines
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! string : input character string : Character string to be used
    !! This routine writes out a subheading between two blank lines.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use constants, only: iotty, nout
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: string

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call oblnkl(file)
    call ocmmnt(file,string)
    call oblnkl(file)

  end subroutine osubhd

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ocmmnt(file,string)

    !! Routine to print a comment
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! string : input character string : Character string to be used
    !! This routine writes out a comment line.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: boundl, boundu, sqsumsq
		use global_variables, only: icase, vlabel, iscan_global
		use constants, only: rmu0
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: string

    !  Local variables

    integer, parameter :: maxwidth = 110
    integer :: lh
!    character(len = maxwidth) :: dummy

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lh = len(trim(string))

    if (lh == 0) then
       write(*,*) 'Error in routine OCMMNT :'
       write(*,*) 'A zero-length string is not permitted.'
       write(*,*) 'PROCESS stopping.'
       stop 1
    end if

    if (lh >= maxwidth) then
       write(*, *) 'Warning in routine OCMMNT :'
       write(*, '(A)') string
!       write(*,*) 'This is too long to fit into ',maxwidth,' columns.'
       write(*, '(A,i3,A)') 'This is longer than ',maxwidth,' columns.'  ! MK 28/10/2016 Modified previous output to reflect warning message
       !write(*,*) 'PROCESS stopping.'
       !stop 1
    end if
!    dummy = trim(string)
    write(file,'(t2,a)') trim(string)
    !write(file,'(t2,a)') dummy
  end subroutine ocmmnt

  subroutine write(file, string)
    !! Write a string to file.
    !! file : input integer : Fortran output unit identifier
    !! string : input character string : Character string to be used
    implicit none

    !  Arguments
    integer, intent(in) :: file
    character(len=*), intent(in) :: string

    write(file,*) trim(string)
  end subroutine write

  subroutine dblcol(file, desc, val1, val2)
    !! Write a description and 2 columns of values to 2dp in standard notation.
    !! file : input integer : Fortran output unit identifier
    !! desc : input character string : Character string to be used
    !! val1 : input real : Value of the left variable
    !! val2 : input real : Value of the right variable
    implicit none

    !  Arguments
    integer, intent(in) :: file
    character(len=70), intent(in) :: desc
    real(8), intent(in) :: val1, val2

    write(file,10) desc, val1, val2
    10 format(1x,a,t75,f10.2,t100,f10.2)
  end subroutine dblcol

  subroutine ovarrf(file,descr,varnam,value,output_flag)

    !! Routine to print out the details of a floating-point
    !! variable using 'F' format
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! descr : input character string : Description of the variable
    !! varnam : input character string : Name of the variable
    !! value : input real : Value of the variable
    !! output_flag : optional character
    !! This routine writes out the description, name and value of a
    !! double precision variable in F format (e.g.
    !! <CODE>-12345.000</CODE>).
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: name_xc
		use global_variables, only: verbose
		use constants, only: pi, mfile, nplot, echarge
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: descr, varnam
    real(8), intent(in) :: value
    character(len=3), intent(in), optional :: output_flag

    !  Local variables

    character(len=72) :: dum72
    character(len=20) :: dum20
    character(len=20) :: stripped
    character(len=3) :: flag
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Replace descr and varnam with dummy strings of the correct length.
    !  This counters problems that would occur if the two original strings
    !  were the wrong length.

    dum72 = descr
    dum20 = varnam
    stripped = varnam(2:len(varnam)-1)

    if (present(output_flag)) then
        flag = output_flag
    else
        flag = ''
    end if

    if (any(name_xc == stripped))  flag = 'ITV'

    if (file /= mfile) then
       !MDK add label if it is an iteration variable
       ! The ITV flag overwrites the output_flag
       if (verbose==1) then
            write(file,10) dum72, dum20, value, flag
       else
            write(file,20) dum72, dum20, value, flag
       end if
    end if

10  format(1x,a,t75,a,t100,f13.6, t115, a)
20  format(1x,a,t75,a,t100,f10.3, t112, a)

    call ovarre(mfile,descr,varnam,value)

  end subroutine ovarrf

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ovarre(file,descr,varnam,value,output_flag)

    !! Routine to print out the details of a floating-point
    !! variable using 'E' format
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! descr : input character string : Description of the variable
    !! varnam : input character string : Name of the variable
    !! value : input real : Value of the variable
    !! output_flag : optional character
    !! This routine writes out the description, name and value of a
    !! double precision variable in E format (e.g.
    !! <CODE>-1.234E+04</CODE>).
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: name_xc
		use global_variables, only: icase, vlabel
		use constants, only: mfile, nout
		use maths_library, only: variable_error
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: descr, varnam
    real(8), intent(in) :: value
    character(len=3), intent(in), optional :: output_flag

    !  Local variables

    character(len=72) :: dum72
    character(len=30) :: dum20
    character(len=20) :: stripped
    character(len=3) :: flag
    integer :: dotindex

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Replace descr and varnam with dummy strings of the correct length.
    !  This counters problems that would occur if the two original strings
    !  were the wrong length.

    dum72 = descr
    dum20 = varnam
    ! Remove the "(" and ")" from the varnam
    stripped = varnam(2:len(varnam)-1)

    ! May need to strip Python module name (e.g. the pfv. from pfv.coheof)
    ! This ensures the ITV flag is still added when required in output files
    dotindex = scan(stripped,".")
    stripped = stripped(dotindex+1:)

    if (present(output_flag)) then
        flag = output_flag
    else
        flag = ''
    end if

    if (any(name_xc == stripped))  flag = 'ITV'

    if (file /= mfile) then
       ! MDK add ITV label if it is an iteration variable
       ! The ITV flag overwrites the output_flag
       write(file,20) dum72, dum20, value, flag
    end if

    call underscore(dum72)
    call underscore(dum20)
    write(mfile,10) dum72, dum20, value, flag

  ! MFILE.DAT format
  ! Machine epsilon for double ~2.22e-16, hence require 17 sig figs in significand
  ! for full precision of a double float
10  format(1x,a,t75,a,t110,ES23.16e2," ",a,t10)
  ! OUT.DAT format
20  format(1x,a,t75,a,t100,1pe10.3, t112, a)

  end subroutine ovarre

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ovarin(file,descr,varnam,value,output_flag)

    !! Routine to print out the details of an integer variable
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! descr : input character string : Description of the variable
    !! varnam : input character string : Name of the variable
    !! value : input integer : Value of the variable
    !! output_flag : optional character
    !! This routine writes out the description, name and value of an
    !! integer variable.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: name_xc, icc, ioptimz
		use global_variables, only: xlabel_2, iscan_global
		use constants, only: mfile, nout
		use maths_library, only: variable_error
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: descr, varnam
    integer, intent(in) :: value
    character(len=3), intent(in), optional :: output_flag

    !  Local variables

    character(len=72) :: dum72
    character(len=30) :: dum20
    character(len=20) :: stripped
    character(len=3) :: flag

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Replace descr and varnam with dummy strings of the correct length.
    !  This counters problems that would occur if the two original strings
    !  were the wrong length.

    dum72 = descr
    dum20 = varnam
    stripped = varnam(2:len(varnam)-1)
    if (present(output_flag)) then
        flag = output_flag
    else
        flag = ''
    end if

    if (any(name_xc == stripped))  flag = 'ITV'

    if (file /= mfile) then
       ! MDK add ITV label if it is an iteration variable
       ! The ITV flag overwrites the output_flag
       write(file,20) dum72, dum20, value, flag
    end if

    call underscore(dum72)
    call underscore(dum20)
    write(mfile,10) dum72, dum20, value, flag

10  format(1x,a,t75,a,t110,i10," ",a,t10)
20  format(1x,a,t75,a,t100,i10,t112, a)

  end subroutine ovarin

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ovarst(file,descr,varnam,value)

    !! Routine to print out the details of a character variable
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! descr : input character string : Description of the variable
    !! varnam : input character string : Name of the variable
    !! value : input character string : Value of the variable
    !! This routine writes out the description, name and value of a
    !! character string variable.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use numerics, only: sqsumsq
		use constants, only: mfile
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: descr, varnam, value

    !  Local variables

    character(len=72) :: dum72
    character(len=30) :: dum20

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Replace descr and varnam with dummy strings of the correct length.
    !  This counters problems that would occur if the two original strings
    !  were the wrong length.

    dum72 = descr
    dum20 = varnam

    if (file /= mfile) then
       write(file,10) dum72, dum20, value
    end if

    call underscore(dum72)
    call underscore(dum20)
    write(mfile,10) dum72, dum20, value

10  format(1x,a,t75,a,t110,a)

  end subroutine ovarst

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ocosts(file,ccode,descr,value)

    !! Routine to print out the code, description and value
    !! of a cost item
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! ccode : input character string : Code number/name of the cost item
    !! descr : input character string : Description of the cost item
    !! value : input real : Value of the cost item
    !! This routine writes out the cost code, description and value
    !! of a cost item.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use constants, only: pi, mfile, nplot, twopi
		use maths_library, only: variable_error
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: ccode, descr
    real(8), intent(in) :: value

    !  Local variables

    character(len=30) :: dum20
    character(len=72) :: dum72

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Replace ccode and descr with dummy strings of the correct length.
    !  This counters problems that would occur if the two original strings
    !  were the wrong length.

    dum20 = ccode
    dum72 = descr

    write(file,10) dum20, dum72, value
10  format(1x,a,t22,a,t110,f10.2)

    call ovarrf(mfile,descr,ccode,value)

  end subroutine ocosts

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine obuild(file,descr,thick,total,variable_name)

    !! Routine to print out a description, the thickness and
    !! summed build of a component of the radial or vertical build
    !! author: P J Knight, CCFE, Culham Science Centre
    !! file : input integer : Fortran output unit identifier
    !! descr : input character string : Description of the component
    !! thick : input real : Thickness of the component (m)
    !! total : input real : Total build, including this component (m)
    !! This routine writes out a description, the thickness and
    !! summed build of a component of the radial or vertical build.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: boundl, boundu
		use constants, only: echarge
		use maths_library, only: variable_error
    implicit none

    !  Arguments

    integer, intent(in) :: file
    character(len=*), intent(in) :: descr
    character(len=*), optional :: variable_name
    real(8), intent(in) :: thick, total

    !  Local variables

    character(len=30) :: dum30
    character(len=20) :: dum20

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Replace descr with dummy string of the correct length.
    !  This counters problems that would occur if the original string
    !  was the wrong length.

    dum30 = descr
    if (present(variable_name)) then
        dum20 = variable_name
    else
        dum20 = ''
    end if

    write(file,10) dum30, thick, total, dum20
10  format(1x,a,t42,f10.3,t58,f10.3,t71,a)

  end subroutine obuild

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine underscore(string)

    !! Routine that converts spaces in a string to underscores
    !! author: P J Knight, CCFE, Culham Science Centre
    !! string : input/output string : character string of interest
    !! This routine converts any space characters in the string
    !! to underscore characters.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use numerics, only: active_constraints, boundu, boundl
		use constants, only: echarge
    implicit none

    !  Arguments

    character(len=*), intent(inout) :: string

    !  Local variables

    integer :: loop, i

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    i = index(string, ' ')
    if (i > 0) then
       do loop = i, len(string)
          if (string(loop:loop) == ' ') string(loop:loop) = '_'
       end do
    end if

  end subroutine underscore

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function int2char(i)

    !! Converts a single-digit integer into a character string
    !! author: P J Knight, CCFE, Culham Science Centre
    !! i : input integer : must be between 0 and 9
    !! This is a very simple routine that converts a single-digit
    !! integer into a character string. If the integer is outside
    !! the range 0 to 9 the program stops with an error.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use numerics, only: epsvmc, boundu
		use constants, only: rmu0
    implicit none

    character(len=1) :: int2char

    !  Arguments

    integer, intent(in) :: i

    !  Local variables

    character(len=10), parameter :: number = '0123456789'

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((i < 0).or.(i > 9)) then
       write(*,*) 'INT2CHAR: illegal argument'
       stop 1
    end if

    int2char = number(i+1:i+1)

  end function int2char

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function int_to_string2(i)

    !! Converts a positive integer into a two-digit character string
    !! author: P J Knight, CCFE, Culham Science Centre
    !! i : input integer : must be between 0 and 99
    !! This routine converts a positive integer into a two-digit
    !! character string.
    !! If the integer is negative, the routine stops with an error.
    !! If the integer is greater than 99, the routine returns a
    !! string containing its last two digits.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use numerics, only: boundu
		use constants, only: pi
    implicit none

    character(len=2) :: int_to_string2

    !  Arguments

    integer, intent(in) :: i

    !  Local variables

    character(len=1) :: a0, a1

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (i < 0) then
       write(*,*) 'INT_TO_STRING2: illegal argument'
       stop 1
    end if

    a0 = int2char(mod(i,10))
    a1 = int2char(mod(int(i/10),10))

    int_to_string2 = a1//a0

  end function int_to_string2

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function int_to_string3(i)

    !! Converts a positive integer into a 3-digit character string
    !! author: P J Knight, CCFE, Culham Science Centre
    !! i : input integer : must be between 0 and 99
    !! This routine converts a positive integer into a three-digit
    !! character string.
    !! If the integer is negative, the routine stops with an error.
    !! If the integer is greater than 999, the routine returns a
    !! string containing its last three digits.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		use numerics, only: boundu
		use constants, only: pi
    implicit none

    character(len=3) :: int_to_string3

    !  Arguments

    integer, intent(in) :: i

    !  Local variables

    character(len=1) :: a0, a1, a2

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (i < 0) then
       write(*,*) 'INT_TO_STRING3: illegal argument'
       stop 1
    end if

    a0 = int2char(mod(i,10))
    a1 = int2char(mod(int(i/10),10))
    a2 = int2char(mod(int(i/100),100))

    int_to_string3 = a2//a1//a0

  end function int_to_string3

end module process_output
