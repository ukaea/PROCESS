module maths_library

  !! Library of mathematical and numerical routines
  !! author: P J Knight, CCFE, Culham Science Centre
  !! N/A
  !! This module contains a large number of routines to enable
  !! PROCESS to perform a variety of numerical procedures, including
  !! linear algebra, zero finding, integration and minimisation.
  !! The module is an amalgamation of the contents of several
  !! different pre-existing PROCESS source files, which themselves
  !! were derived from a number of different numerical libraries
  !! including BLAS and MINPAC.
  !! http://en.wikipedia.org/wiki/Gamma_function
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  !use process_output

  implicit none

  !  Precision variable
  integer, parameter :: double = 8

contains

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function find_y_nonuniform_x(x0,x,y,n)

    !! Routine to find y0 such that y0 = y(x0) given a set of
    !! values x(1:n), y(1:n)
    !! author: P J Knight, CCFE, Culham Science Centre
    !! x0 : input real : x value at which we want to find y
    !! x(1:n) : input real array : monotonically de/increasing x values
    !! y(1:n) : input real array : y values at each x
    !! n : input integer : size of array
    !! This routine performs a simple linear interpolation method
    !! to find the y value at x = x0. If x0 lies outside the
    !! range [x(1),x(n)], the y value at the nearest 'end' of the data
    !! is returned.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(dp) :: find_y_nonuniform_x

    !  Arguments

    integer, intent(in) :: n
    real(dp), intent(in) :: x0
    real(dp), dimension(n), intent(in) :: x
    real(dp), dimension(n), intent(in) :: y

    !  Local variables
    integer :: i,j

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Step through arrays until x crosses the value of interest

    j = 0
    rough_search: do i = 1,n-1

       if (((x(i)-x0)*(x(i+1)-x0)) <= 0.0D0) then
          j = i
          exit rough_search
       end if

    end do rough_search

    if (j /= 0) then

       !  Simply do a linear interpolation between the two grid points
       !  spanning the point of interest

       find_y_nonuniform_x = y(j) + (y(j+1)-y(j))*(x0-x(j))/(x(j+1)-x(j))

    else  !  No points found, so return the 'nearest' y value

       if (x(n) > x(1)) then  !  values are monotonically increasing
          if (x0 > x(n)) then
             find_y_nonuniform_x = y(n)
          else
             find_y_nonuniform_x = y(1)
          end if

       else  !  values are monotonically decreasing
          if (x0 < x(n)) then
             find_y_nonuniform_x = y(n)
          else
             find_y_nonuniform_x = y(1)
          end if

       end if

    end if

  end function find_y_nonuniform_x

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sumup2(dx,y,inty,n)

    !! Routine to integrate a 1-D array of y values, using a process
    !! similar to Simpson's Rule, and assuming equally-spaced x values.
    !! It returns the integral at all tabulated points.
    !! author: P J Knight, CCFE, Culham Science Centre
    !! dx : input real : (constant) spacing between adjacent x values
    !! y(1:n) : input real array : y values to be integrated
    !! inty(1:n) : input/output real array : calculated integral
    !! (see below)
    !! n : input integer : length of arrays y and inty
    !! This routine uses a process similar to (but not quite the same
    !! as) Simpson's Rule to integrate an array y,
    !! returning the integral up to point i in array element inty(i).
    !! Note that the first element of inty is not calculated, and must
    !! be set to the required value on entry. Usually, but not always,
    !! this value will be zero.
    !! The original source for this algorithm is not known...
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: n
    real(dp), intent(in) :: dx
    real(dp), intent(in), dimension(n) :: y
    real(dp), intent(inout), dimension(n) :: inty

    !  Local variables

    integer :: ix
    real(dp), parameter :: third = 1.0D0/3.0D0
    real(dp) :: thirddx
    real(dp), allocatable, dimension(:) :: yhalf

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    allocate (yhalf(2:n))

    thirddx = third*dx

    do ix = 2,n-1
       yhalf(ix) = y(ix)-0.25D0*(y(ix+1)-y(ix-1))
    end do
    yhalf(n) = y(n-1)+0.25D0*(y(n)-y(n-2))

    do ix = 2,n
       inty(ix) = inty(ix-1) + thirddx*(y(ix)+yhalf(ix)+y(ix-1))
    end do

    deallocate (yhalf)

  end subroutine sumup2

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sumup3(dx,y,integral,n)

    !! Routine to integrate a 1-D array of y values using the
    !! Extended Simpson's Rule, assuming equally-spaced x values
    !! author: P J Knight, CCFE, Culham Science Centre
    !! dx : input real : (constant) spacing between adjacent x values
    !! y(1:n) : input real array : y values to be integrated
    !! integral : output real : calculated integral
    !! n : input integer : length of array y
    !! This routine uses Simpson's Rule to integrate an array y.
    !! If n is even, routine <CODE>sumup2</CODE> is called to
    !! perform the calculation.
    !! <P>Note: unlike sumup1 and sumup2, this routine returns only
    !! the complete integral, not the intermediate values as well.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: n
    real(dp), intent(in) :: dx
    real(dp), intent(in), dimension(n) :: y
    real(dp), intent(out) :: integral

    !  Local variables

    integer :: ix
    real(dp) :: sum1
    real(dp), allocatable, dimension(:) :: inty

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (mod(n,2) == 0) then

       !  Use sumup2 if the number of tabulated points is even

       allocate (inty(n))

       inty(1) = 0.0D0
       call sumup2(dx,y,inty,n)
       integral = inty(n)

       deallocate (inty)

    else

       sum1 = y(1)
       do ix = 2,n-3,2
          sum1 = sum1 + 4.0D0*y(ix) + 2.0D0*y(ix+1)
       end do
       integral = dx/3.0D0*(sum1 + 4.0D0*y(n-1) + y(n))

    end if

  end subroutine sumup3

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine tril(a,n,alower)

    !! Routine to extract the lower triangular part of a square matrix
    !! author: P J Knight, CCFE, Culham Science Centre
    !! a(n,n) : input real array : input matrix
    !! n      : input integer : number of rows and columns in A
    !! a(n,n) : output real array : lower triangular part of A
    !! This routine extracts the lower triangular part of a square matrix,
    !! excluding the diagonal, into a new square matrix. The remainder
    !! of this matrix contains zeroes on exit.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: n
    real(dp), dimension(n,n), intent(in) :: a
    real(dp), dimension(n,n), intent(out) :: alower

    !  Local variables

    integer :: row,col

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    alower = 0.0D0
    do col = 1,n-1
       do row = col+1,n
          alower(row,col) = a(row,col)
       end do
    end do

  end subroutine tril

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine ellipke(sqk,kk,ek)

    !! Routine that calculates the complete elliptic integral
    !! of the first and second kinds
    !! author: P J Knight, CCFE, Culham Science Centre
    !! sqk : input real : square of the elliptic modulus
    !! kk  : output real : complete elliptic integral of the first kind
    !! ek  : output real : complete elliptic integral of the second kind
    !! This routine calculates the complete elliptic integral
    !! of the first and second kinds.
    !! <P>The method used is that described in the reference, and
    !! the code is taken from the Culham maglib library routine FN02A.
    !! Approximations for Digital Computers, C. Hastings,
    !! Princeton University Press, 1955
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    real(dp), intent(in) :: sqk
    real(dp), intent(out) :: kk,ek

    !  Local variables

    real(dp) :: a,b,d,e

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    d = 1.0D0 - sqk
    e = log(d)

    !  Evaluate series for integral of first kind

    a = (((0.014511962D0*d + 0.037425637D0)*d + 0.035900924D0)*d &
         + 0.096663443D0)*d + 1.386294361D0
    b = (((0.004417870D0*d + 0.033283553D0)*d + 0.06880249D0)*d &
         + 0.12498594D0)*d + 0.5D0

    kk = a - b*e

    !  Evaluate series for integral of second kind

    a = (((0.017365065D0*d + 0.047573835D0)*d + 0.06260601D0)*d &
         + 0.44325141D0)*d + 1.0D0
    b = (((0.005264496D0*d + 0.040696975D0)*d + 0.09200180D0)*d &
         + 0.24998368D0)*d

    ek = a - b*e

  end subroutine ellipke

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function binomial(n,k) result(coefficient)
    ! This outputs a real approximation to the coefficient
    ! http://en.wikipedia.org/wiki/Binomial_coefficient#Multiplicative_formula
    implicit none
    integer, intent(in) :: n, k
    real(dp) :: coefficient
    integer :: numerator, i
    if (k == 0) then
        coefficient = 1
    else
        coefficient = 1.0D0
        do i = 1, k
            numerator = n + 1 -i
            coefficient = coefficient * real(numerator)/real(i)
        end do
    end if
  end function binomial
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  recursive function gamfun(x) result(gamma)

    !! Calculates the gamma function for arbitrary real x
    !! author: P J Knight, CCFE, Culham Science Centre
    !! x : input real : gamma function argument
    !! This routine evaluates the gamma function, using an
    !! asymptotic expansion based on Stirling's approximation.
    !! http://en.wikipedia.org/wiki/Gamma_function
    !! T&amp;M/PKNIGHT/LOGBOOK24, p.5
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    real(dp), intent(in) :: x
    real(dp) :: gamma

    !  Local variables

    real(dp), parameter :: sqtwopi = 2.5066282746310005D0
    real(dp), parameter :: c1 = 8.3333333333333333D-2  !  1/12
    real(dp), parameter :: c2 = 3.4722222222222222D-3  !  1/288
    real(dp), parameter :: c3 = 2.6813271604938272D-3  !  139/51840
    real(dp), parameter :: c4 = 2.2947209362139918D-4  !  571/2488320
    real(dp) :: summ, denom
    integer :: i,n

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (x > 1.0D0) then

       summ = 1.0D0 + c1/x + c2/x**2 - c3/x**3 - c4/x**4
       gamma = exp(-x) * x**(x-0.5D0) * sqtwopi * summ

    else

       !  Use recurrence formula to shift the argument to >1
       !  gamma(x) = gamma(x+n) / (x*(x+1)*(x+2)*...*(x+n-1))
       !  where n is chosen to make x+n > 1

       n = int(-x) + 2
       denom = x
       do i = 1,n-1
          denom = denom*(x+i)
       end do
       gamma = gamfun(x+n)/denom

    end if

  end function gamfun

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function binarysearch(length, array, value, delta)
    ! Given an array and a value, returns the index of the element that
    ! is closest to, but less than, the given value.
    ! Uses a binary search algorithm.
    ! "delta" is the tolerance used to determine if two values are equal
    ! if ( abs(x1 - x2) <= delta) then
    !    assume x1 = x2
    ! endif
    ! Should never return an index < 1 or > length!

    implicit none
    integer, intent(in) :: length
    real(dp), dimension(length), intent(in) :: array
    real(dp), intent(in) :: value
    real(dp), intent(in), optional :: delta
    integer :: binarysearch
    integer :: left, middle, right
    real(dp) :: d

    if (present(delta) .eqv. .true.) then
        d = delta
    else
        d = 1e-9
    endif

    left = 1
    right = length
    do
        if (left > right) then
            exit
        endif
        middle = nint((real(left+right)) / 2.0)
        if ( abs(array(middle) - value) <= d) then
            binarySearch = middle
            return
        else if (array(middle) > value) then
            right = middle - 1
        else
            left = middle + 1
        end if
    end do
    binarysearch = min(max(right,1),length)

  end function binarysearch

  real(dp) function interpolate(x_len, x_array, y_len, y_array, f, x, y)
    ! This function uses bilinear interpolation to estimate the value
    ! of a function f at point (x,y)
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by x_array and the grid y values specified by y_array
    ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation
    implicit none
    integer, intent(in) :: x_len, y_len
    real(dp), dimension(x_len), intent(in) :: x_array
    real(dp), dimension(y_len), intent(in) :: y_array
    real(dp), dimension(x_len, y_len), intent(in) :: f
    real(dp), intent(in) :: x,y
    real(dp) :: denom, x1, x2, y1, y2
    integer :: i,j

    i = binarysearch(x_len, x_array, x)
    j = binarysearch(y_len, y_array, y)

    if (i  >= x_len) then
       i = x_len -1
    end if
    if (j >= y_len) then
       j = y_len-1
    end if
    x1 = x_array(i)
    x2 = x_array(i+1)

    y1 = y_array(j)
    y2 = y_array(j+1)

    denom = (x2 - x1)*(y2 - y1)

    interpolate = (f(i,j)*(x2-x)*(y2-y) + f(i+1,j)*(x-x1)*(y2-y) + &
      f(i,j+1)*(x2-x)*(y-y1) + f(i+1, j+1)*(x-x1)*(y-y1))/denom

  end function interpolate

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine quanc8(fun,a,b,abserr,relerr,result,errest,nofun,flag)

    !! Estimate the integral of fun(x) from a to b
    !! author: P J Knight, CCFE, Culham Science Centre
    !! fun : external function : integrand function subprogram fun(x)
    !! a : input real : lower limit of integration
    !! b : input real : upper limit of integration (b may be less than a)
    !! abserr : input real : absolute error tolerance (should be non-negative)
    !! relerr : input real : relative error tolerance (should be non-negative)
    !! result : output real : approximation to the integral hopefully
    !! satisfying the least stringent of the two error tolerances
    !! errest : output real : estimate of the magnitude of the actual error
    !! nofun : output integer : number of function values used in calculation
    !! flag : output real : Reliability indicator; if flag is zero, then
    !! result probably satisfies the error tolerance.  If flag is
    !! xxx.yyy , then  xxx = the number of intervals which have
    !! not converged and  0.yyy = the fraction of the interval
    !! left to do when the limit on  nofun  was approached.
    !! This routine estimates the integral of fun(x) from a to b
    !! to a user provided tolerance. An automatic adaptive
    !! routine based on the 8-panel Newton-Cotes rule.
    !! http://www.netlib.org/fmm/index.html :
    !! Computer Methods for Mathematical Computations,
    !! G E Forsythe, M A Malcolm, and C B Moler,
    !! Prentice-Hall, Englewood Cliffs, New Jersey
    !! 1977, ISBN 0-13-165332-6
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    interface
      function fun(rho) result(fint)
#ifndef dp
        use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
        real(dp), intent(in) :: rho
        real(dp) :: fint
      end function fun
    end interface

    !  Arguments

    external :: fun
    real(dp), intent(in) :: a, b, abserr, relerr
    real(dp), intent(out) :: result, errest, flag
    integer, intent(out) :: nofun

    !  Local variables

    real(dp) :: w0,w1,w2,w3,w4,area,x0,f0,stone,step,cor11,temp
    real(dp) :: qprev,qnow,qdiff,qleft,esterr,tolerr
    real(dp), dimension(31) :: qright
    real(dp), dimension(16) :: f, x
    real(dp), dimension(8,30) :: fsave, xsave

    integer :: levmin,levmax,levout,nomax,nofin,lev,nim,i,j

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Stage 1:  general initialization
    !  Set constants

    levmin = 1
    levmax = 30
    levout = 6
    nomax = 5000
    nofin = nomax - 8*(levmax-levout+2**(levout+1))

    !  Trouble when nofun reaches nofin

    w0 =   3956.0D0 / 14175.0D0
    w1 =  23552.0D0 / 14175.0D0
    w2 =  -3712.0D0 / 14175.0D0
    w3 =  41984.0D0 / 14175.0D0
    w4 = -18160.0D0 / 14175.0D0

    !  Initialize running sums to zero

    flag = 0.0D0
    result = 0.0D0
    cor11  = 0.0D0
    errest = 0.0D0
    area  = 0.0D0
    nofun = 0

    if (a == b) return

    !  Stage 2:  initialization for first interval

    lev = 0
    nim = 1
    x0 = a
    x(16) = b
    qprev  = 0.0D0
    f0 = fun(x0)
    stone = (b - a) / 16.0D0
    x(8)  =  (x0  + x(16)) / 2.0D0
    x(4)  =  (x0  + x(8))  / 2.0D0
    x(12) =  (x(8)  + x(16)) / 2.0D0
    x(2)  =  (x0  + x(4))  / 2.0D0
    x(6)  =  (x(4)  + x(8))  / 2.0D0
    x(10) =  (x(8)  + x(12)) / 2.0D0
    x(14) =  (x(12) + x(16)) / 2.0D0
    do j = 2, 16, 2
       f(j) = fun(x(j))
    end do
    nofun = 9

    !  Stage 3:  central calculation

    !  Requires qprev,x0,x2,x4,...,x16,f0,f2,f4,...,f16
    !  Calculates x1,x3,...x15, f1,f3,...f15,qleft,qright,qnow,qdiff,area

    main_loop: do

       x(1) = (x0 + x(2)) / 2.0D0
       f(1) = fun(x(1))
       do j = 3, 15, 2
          x(j) = (x(j-1) + x(j+1)) / 2.0D0
          f(j) = fun(x(j))
       end do
       nofun = nofun + 8
       step = (x(16) - x0) / 16.0D0
       qleft = (w0*(f0 + f(8))  + w1*(f(1)+f(7))  + w2*(f(2)+f(6)) &
            + w3*(f(3)+f(5))  +  w4*f(4)) * step
       qright(lev+1) = (w0*(f(8)+f(16))+w1*(f(9)+f(15))+w2*(f(10)+f(14)) &
            + w3*(f(11)+f(13)) + w4*f(12)) * step
       qnow = qleft + qright(lev+1)
       qdiff = qnow - qprev
       area = area + qdiff

       !  Stage 4:  interval convergence test

       esterr = abs(qdiff) / 1023.0D0
       tolerr = max(abserr,relerr*abs(area)) * (step/stone)

       if ( (lev < levmin).or. &
            ((lev < levmax).and.(nofun <= nofin) &
            .and.(esterr > tolerr)) ) then

          !  Stage 5:  no convergence
          !  Locate next interval

          nim = 2*nim
          lev = lev+1

          !  Store right hand elements for future use

          do i = 1, 8
             fsave(i,lev) = f(i+8)
             xsave(i,lev) = x(i+8)
          end do

          !  Assemble left hand elements for immediate use

          qprev = qleft
          do i = 1, 8
             j = -i
             f(2*j+18) = f(j+9)
             x(2*j+18) = x(j+9)
          end do

          cycle main_loop

       else if (lev >= levmax) then

          flag = flag + 1.0D0

       else if (nofun > nofin) then

          !  Stage 6:  trouble section
          !  Number of function values is about to exceed limit

          nofin = 2*nofin
          levmax = levout
          flag = flag + (b - x0) / (b - a)

       end if

       !  Stage 7:  interval converged
       !  Add contributions into running sums

       result = result + qnow
       errest = errest + esterr
       cor11  = cor11  + qdiff / 1023.0D0

       !  Locate next interval

       do
          if (nim == (2*(nim/2))) exit
          nim = nim/2
          lev = lev-1
       end do
       nim = nim + 1
       if (lev <= 0) exit main_loop

       !  Assemble elements required for the next interval

       qprev = qright(lev)
       x0 = x(16)
       f0 = f(16)
       do i = 1, 8
          f(2*i) = fsave(i,lev)
          x(2*i) = xsave(i,lev)
       end do

    end do main_loop

    !  Stage 8:  finalize and return

    result = result + cor11

    !  Make sure errest not less than roundoff level

    if (errest == 0.0D0) return
    estimate_error: do
       temp = abs(result) + errest
       if (temp /= abs(result)) exit estimate_error
       errest = 2.0D0*errest
    end do estimate_error

  end subroutine quanc8

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine linesolv(a, ndim, b, x)

    !! Routine to solve the linear equation system Ax = b
    !! author: P J Knight, CCFE, Culham Science Centre
    !! a(ndim,ndim) : in/out real array : array A
    !! ndim         : input integer     : dimension of a
    !! b(ndim)      : in/out real array : RHS vector
    !! x(ndim)      : output real array : solution for Ax = b
    !! This routine solves the linear equation Ax = b.
    !! It calls (local copies of) the linpack routines sgefa and sgesl.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: ndim
    real(dp), dimension(ndim,ndim), intent(inout) :: a
    real(dp), dimension(ndim), intent(inout) :: b, x

    !  Local variables

    integer :: job, ndim1, info
    integer, dimension(ndim) :: ipvt

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    job = 0
    ndim1 = ndim

    call sgefa(a, ndim, ndim1, ipvt, info)
    call sgesl(a, ndim, ndim1, ipvt, b, job)

    x(:) = b(:)

  end subroutine linesolv

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine hinv(h,ih,n,ipvt,info)

    !! Matrix inversion routine
    !! author: Roger L. Crane, Kenneth E. Hillstrom, Michael Minkoff; Linpack
    !! author: P J Knight, CCFE, Culham Science Centre
    !! h(ih,ih) : input/output real array : On input, matrix to be inverted
    !! On output, the calculated inverse
    !! ih       : input integer : array size
    !! n        : input integer : order of H; n <= ih
    !! ipvt(n)  : output integer array : pivot vector
    !! info     : output integer : info flag
    !! = 1  normal return
    !! = 2  H matrix is singular
    !! This routine inverts the matrix H by use of linpack software.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use global_variables, only: verbose
    implicit none

    !  Arguments

    integer, intent(in) :: ih, n
    integer, intent(out) :: info
    integer, dimension(:), intent(out) :: ipvt
    real(dp), dimension(:,:), intent(inout) :: h

    !  Local variables

    real(dp), dimension(2) :: det

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Do LU decomposition of h

    call sgefa(h,ih,n,ipvt,info)

    if (info == 0) then  !  H is non-singular, so we can form its inverse
       call sgedi(h,ih,n,ipvt,det,1)
       info = 1
    else
       info = 2
       if (verbose == 1) then
          call sgedi(h,ih,n,ipvt,det,10)  !  Calculate determinant only
          write(*,*) 'Determinant = det(1) * 10.0**det(2)',det
          write(*,*) 'H matrix is singular in subroutine hinv'
       end if
    end if

  end subroutine hinv

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dotpmc(x,ix,y,iy,c,total,n,iflag)

    !! Calculates +/-C +/- (X.dot.Y) for arrays X, Y and scalar C
    !! author: Roger L. Crane, Kenneth E. Hillstrom, Michael Minkoff; Linpack
    !! author: P J Knight, CCFE, Culham Science Centre
    !! x(ix*n) : input real array : X array
    !! ix      : input integer : interval in storage between X array elements
    !! y(iy*n) : input real array : Y array
    !! iy      : input integer : interval in storage between Y array elements
    !! c       : input real : C value
    !! total   : output real : computed result
    !! n       : input integer : number of terms in the dot product
    !! iflag   : input integer : switch
    !! = 0    +c + (x dot y) is computed
    !! = 1    +c - (x dot y) is computed
    !! = 2    -c + (x dot y) is computed
    !! = 3    -c - (x dot y) is computed
    !! This subroutine computes
    !! total = (plus or minus c) plus or minus the dot product of x and y
    !! by invoking the basic linear algebra routine dot.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: ix,iy,n,iflag
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(:), intent(in) :: y
    real(dp), intent(in) :: c
    real(dp), intent(out) :: total

    !  Local variables

    real(dp) :: prod

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Calculate dot product

    prod = sdot(n,x,ix,y,iy)
    if (mod(iflag,2) /= 0) prod = -prod

    total = c + prod
    if (iflag > 1) total = -c + prod

  end subroutine dotpmc

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sgesl(a,lda,n,ipvt,b,job)

    !! Routine to solve the the real system  Ax = b  or  transp(A).x = b
    !! author: Cleve Moler, University of New Mexico, Argonne National Lab.
    !! author: P J Knight, CCFE, Culham Science Centre
    !! a(lda,n) : input real array : output from <A HREF="sgefa.html">sgefa</A>
    !! lda : input integer : leading dimension of the array A
    !! n : input integer : order of the matrix A
    !! ipvt(n) : input integer array : pivot vector from <CODE>sgefa</CODE>
    !! b(n) : input/output real array : RHS vector on input,
    !! solution vector x on output
    !! job : input integer : switch
    !! = 0         to solve  A*x = b ,
    !! = nonzero   to solve  transp(A)*x = b  where
    !! transp(A)  is the transpose
    !! This routine solves the real system  A*x = b  or  transp(A)*x = b
    !! using the factors computed by <A HREF="sgefa.html">sgefa</A>.
    !! <P>A division by zero will occur if the input factor contains a
    !! zero on the diagonal.  Technically this indicates singularity
    !! but it is often caused by improper arguments or improper
    !! setting of <CODE>lda</CODE>.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: lda,n,job
    integer, dimension(n), intent(in) :: ipvt
    real(dp), dimension(lda,n), intent(in) :: a
    real(dp), dimension(n), intent(inout) :: b

    !  Local variables

    integer :: k,kb,l,nm1
    real(dp) :: t

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    nm1 = n - 1

    if (job == 0) then  !  Solve  A * x = b

       !  First solve  l*y = b

       if (nm1 >= 1) then
          do k = 1, nm1
             l = ipvt(k)
             t = b(l)
             if (l /= k) then
                b(l) = b(k)
                b(k) = t
             end if
             call saxpy(n-k,t,a(k+1:n,k),1,b(k+1:n),1)
          end do
       end if

       !  Now solve  u*x = y

       do kb = 1, n
          k = n + 1 - kb
          b(k) = b(k)/a(k,k)
          t = -b(k)
          call saxpy(k-1,t,a(1:n,k),1,b(1:n),1)
       end do

    else  !  Solve  transp(A) * x = b

       !  First solve  transp(u)*y = b

       do k = 1, n
          t = sdot(k-1,a(:,k),1,b(:),1)
          b(k) = (b(k) - t)/a(k,k)
       end do

       !  Now solve transp(l)*x = y

       if (nm1 >= 1) then
          do kb = 1, nm1
             k = n - kb
             b(k) = b(k) + sdot(n-k,a(k+1:,k),1,b(k+1:),1)
             l = ipvt(k)
             if (l /= k) then
                t = b(l)
                b(l) = b(k)
                b(k) = t
             end if
          end do
       end if

    end if

  end subroutine sgesl

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sgefa(a,lda,n,ipvt,info)

    !! Routine to factor a real matrix by Gaussian elimination
    !! author: Cleve Moler, University of New Mexico, Argonne National Lab.
    !! author: P J Knight, CCFE, Culham Science Centre
    !! a(lda,n) : input/output real array : On entry, matrix to be factored.
    !! On exit, an upper triangular matrix and the multipliers
    !! which were used to obtain it.
    !! The factorization can be written  A = L*U  where
    !! L is a product of permutation and unit lower
    !! triangular matrices and U is upper triangular.
    !! lda : input integer : leading dimension of the array A
    !! n : input integer : order of the matrix A
    !! ipvt(n) : output integer array : pivot indices
    !! info : output integer : info flag
    !! = 0  normal completion
    !! = k  if  u(k,k) == 0.0
    !! This routine factors a real matrix by Gaussian elimination.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use global_variables, only: verbose
    implicit none

    !  Arguments

    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    integer, dimension(n), intent(out) :: ipvt
    real(dp), dimension(lda,n), intent(inout) :: a

    !  Local variables

    integer :: j,k,kp1,l,nm1
    real(dp) :: t

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (verbose == 1) then
       do j = 1,n
          if (all(a(j,:) == 0.0D0)) then
             write(*,*) 'Line ',j, &
                  ' in matrix a in subroutine sgefa is all zero'
          end if
       end do
    end if

    info = 0
    nm1 = n - 1

    if (nm1 >= 1) then

       do k = 1, nm1
          kp1 = k + 1

          !  Find L = pivot index

          l = isamax(n-k+1,a(k,k),1) + k - 1
          ipvt(k) = l

          !  Zero pivot implies this column already triangularized

          if (a(l,k) /= 0.0D0) then

             !  Interchange if necessary

             if (l /= k) then
                t = a(l,k)
                a(l,k) = a(k,k)
                a(k,k) = t
             end if

             !  Compute multipliers

             t = -1.0D0/a(k,k)
             call sscal(n-k,t,a(k+1:n,k),1)

             !  Row elimination with column indexing

             do j = kp1, n
                t = a(l,j)
                if (l /= k) then
                   a(l,j) = a(k,j)
                   a(k,j) = t
                end if
                call saxpy(n-k,t,a(k+1:n,k),1,a(k+1:n,j),1)
             end do

          else
             info = k
             if (verbose == 1) then
                write(*,*) 'a(l,k) = 0.0D0 in subroutine sgefa'
                write(*,*) 'info=k=',info
             end if
          end if
       end do

    end if

    ipvt(n) = n
    if (a(n,n) == 0.0D0) then
       info = n
       if (verbose == 1) then
          write(*,*) 'Error: a(n,n) == 0.0D0 in subroutine sgefa'
       end if
    end if

  end subroutine sgefa

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sgedi(a,lda,n,ipvt,det,job)

    !! Routine to compute the determinant and inverse of a matrix
    !! author: Cleve Moler, University of New Mexico, Argonne National Lab.
    !! author: P J Knight, CCFE, Culham Science Centre
    !! a(lda,n) : input/output real array :
    !! On entry, output from <A HREF="sgefa.html">sgefa</A>.
    !! On exit, the inverse if requested, otherwise unchanged
    !! lda      : input integer : leading dimension of the array A
    !! n        : input integer : order of the matrix A
    !! ipvt(n)  : input integer array : pivot vector from sgefa
    !! det(2)   : output real array : determinant of original matrix if requested,
    !! otherwise not referenced.
    !! Determinant = det(1) * 10.0**det(2)
    !! with  1.0 .le. abs(det(1)) .lt. 10.0
    !! or  det(1) .eq. 0.0 .
    !! job : input integer : switch for required outputs
    !! = 11   both determinant and inverse.
    !! = 01   inverse only.
    !! = 10   determinant only.
    !! This routine computes the determinant and inverse of a matrix
    !! using the factors computed by (SGECO or) <A HREF="sgefa.html">SGEFA</A>.
    !! <P>A division by zero will occur if the input factor contains
    !! a zero on the diagonal and the inverse is requested.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: lda,n,job
    integer, dimension(n), intent(in) :: ipvt
    real(dp), dimension(lda,n), intent(inout) :: a

    !  Local variables

    integer :: i,j,k,kk,kb,kp1,l,nm1
    real(dp), parameter :: ten = 10.0D0
    real(dp) :: t
    real(dp), dimension(2) :: det
    real(dp), dimension(n) :: work

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((job/10) /= 0) then  !  Compute determinant

       det(1) = 1.0D0
       det(2) = 0.0D0

       do i = 1, n
          if (ipvt(i) /=  i) det(1) = -det(1)
          det(1) = a(i,i)*det(1)

          if (det(1) == 0.0D0) exit

          do
             if (abs(det(1)) >= 1.0D0) exit
             det(1) = ten*det(1)
             det(2) = det(2) - 1.0D0
          end do

          do
             if (abs(det(1)) < ten) exit
             det(1) = det(1)/ten
             det(2) = det(2) + 1.0D0
          end do
       end do

    end if

    !  Compute inverse(u)

    if (mod(job,10) /= 0) then

       do k = 1, n
          a(k,k) = 1.0D0/a(k,k)
          t = -a(k,k)

          call sscal(k-1,t,a(1:n,k),1)
          kp1 = k + 1
          if (n >= kp1) then
             do j = kp1, n
                t = a(k,j)
                a(k,j) = 0.0D0
                kk = k
                call saxpy(kk,t,a(1:n,k),1,a(1:n,j),1)
             end do
          end if
       end do

       !  Form inverse(u)*inverse(l)

       nm1 = n - 1
       if (nm1 >= 1) then
          do kb = 1, nm1
             k = n - kb
             kp1 = k + 1

             do i = kp1, n
                work(i) = a(i,k)
                a(i,k) = 0.0D0
             end do

             do j = kp1, n
                t = work(j)
                call saxpy(n,t,a(1:n,j),1,a(1:n,k),1)
             end do

             l = ipvt(k)

             if (l /= k) call sswap(n,a(1:n,k),1,a(1:n,l),1)
          end do

       end if
    end if

  end subroutine sgedi

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sscal(n,sa,sx,incx)

    !! Routine to scale a vector by a constant
    !! author: Jack Dongarra, Linpack
    !! author: P J Knight, CCFE, Culham Science Centre
    !! n        : input integer : order of the matrix sx
    !! sa       : input real array : constant multiplier
    !! sx(n*incx) : input/output real array : On entry, matrix to be scaled;
    !! On exit, the scaled matrix
    !! incx     : input integer : interval in storage between sx array elements
    !! This routine scales a vector by a constant, using
    !! unrolled loops for increments equal to 1.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: n, incx
    real(dp), intent(in) :: sa
    real(dp), dimension(n*incx), intent(inout) :: sx

    !  Local variables

    integer :: i,ix,m,mp1

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (n <= 0) return

    if (incx /= 1) then

       ix = 1
       if (incx < 0) ix = (-n+1)*incx + 1
       do i = 1,n
          sx(ix) = sa*sx(ix)
          ix = ix + incx
       end do

    else

       m = mod(n,5)
       if ( m /= 0 ) then
          do i = 1,m
             sx(i) = sa*sx(i)
          end do
          if (n < 5) return
       end if

       mp1 = m + 1
       do i = mp1,n,5
          sx(i)     = sa*sx(i)
          sx(i + 1) = sa*sx(i + 1)
          sx(i + 2) = sa*sx(i + 2)
          sx(i + 3) = sa*sx(i + 3)
          sx(i + 4) = sa*sx(i + 4)
       end do

    end if

  end subroutine sscal

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine saxpy(n,sa,sx,incx,sy,incy)

    !! Routine to scale a vector by a constant, then add another vector
    !! author: Jack Dongarra, Linpack
    !! author: P J Knight, CCFE, Culham Science Centre
    !! n        : input integer : order of the matrices sx, sy
    !! sa       : input real array : constant multiplier
    !! sx(n*incx) : input real array : matrix to be scaled
    !! incx     : input integer : interval in storage between sx array elements
    !! sy(n*incy) : input/output real array : On entry, matrix being added;
    !! On exit, the final result
    !! incy     : input integer : interval in storage between sy array elements
    !! This routine calculates <CODE>sa*sx(:) + sy(:)</CODE>,
    !! using unrolled loops for increments equal to 1.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: n,incx,incy
    real(dp), intent(in) :: sa
    real(dp), dimension(n*incx), intent(in) :: sx
    real(dp), dimension(n*incy), intent(inout) :: sy

    !  Local variables

    integer :: i,ix,iy,m,mp1

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if ((n <= 0).or.(sa == 0.0D0)) return

    if ((incx /= 1).or.(incy /= 1)) then

       ix = 1 ; iy = 1
       if (incx < 0) ix = (-n+1)*incx + 1
       if (incy < 0) iy = (-n+1)*incy + 1
       do i = 1,n
          sy(iy) = sy(iy) + sa*sx(ix)
          ix = ix + incx
          iy = iy + incy
       end do

    else

       m = mod(n,4)
       if (m /= 0) then
          do i = 1,m
             sy(i) = sy(i) + sa*sx(i)
          end do
          if (n < 4) return
       end if

       mp1 = m + 1
       do i = mp1,n,4
          sy(i)     = sy(i)     + sa*sx(i)
          sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
          sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
          sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
       end do

    end if

  end subroutine saxpy

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine sswap(n,sx,incx,sy,incy)

    !! Routine to interchange two vectors
    !! author: Jack Dongarra, Linpack
    !! author: P J Knight, CCFE, Culham Science Centre
    !! n        : input integer : order of the matrices sx, sy
    !! sx(n*incx) : input/output real array : first vector
    !! incx     : input integer : interval in storage between sx array elements
    !! sy(n*incy) : input/output real array : second vector
    !! incy     : input integer : interval in storage between sy array elements
    !! This routine swaps the contents of two vectors,
    !! using unrolled loops for increments equal to 1.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: n, incx, incy
    real(dp), dimension(n*incx), intent(inout) :: sx
    real(dp), dimension(n*incy), intent(inout) :: sy

    !  Local variables

    integer :: i,ix,iy,m,mp1
    real(dp) :: stemp

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (n <= 0) return

    if ((incx /= 1).or.(incy /= 1)) then

       ix = 1 ; iy = 1
       if (incx < 0) ix = (-n+1)*incx + 1
       if (incy < 0) iy = (-n+1)*incy + 1
       do i = 1,n
          stemp = sx(ix)
          sx(ix) = sy(iy)
          sy(iy) = stemp
          ix = ix + incx
          iy = iy + incy
       end do

    else

       m = mod(n,3)
       if (m /= 0) then
          do i = 1,m
             stemp = sx(i)
             sx(i) = sy(i)
             sy(i) = stemp
          end do
          if (n < 3) return
       end if

       mp1 = m + 1
       do i = mp1,n,3
          stemp = sx(i)
          sx(i) = sy(i)
          sy(i) = stemp
          stemp = sx(i + 1)
          sx(i + 1) = sy(i + 1)
          sy(i + 1) = stemp
          stemp = sx(i + 2)
          sx(i + 2) = sy(i + 2)
          sy(i + 2) = stemp
       end do

    end if

  end subroutine sswap

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function sdot(n,sx,incx,sy,incy)

    !! Routine to compute X*Y where X and Y are vectors
    !! author: Jack Dongarra, Linpack
    !! author: P J Knight, CCFE, Culham Science Centre
    !! n        : input integer : order of the matrices sx, sy
    !! sx(n*incx) : input real array : first vector
    !! incx     : input integer : interval in storage between sx array elements
    !! sy(n*incy) : input real array : second vector
    !! incy     : input integer : interval in storage between sy array elements
    !! This routine performs the dot product of two vectors, i.e.
    !! calculates the sum from i=1 to N, of X(i)*Y(i).
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(dp) :: sdot

    !  Arguments

    integer, intent(in) :: n,incx,incy
    real(dp), dimension(:), intent(in) :: sx
    real(dp), dimension(:), intent(in) :: sy

    !  Local variables

    integer :: ix,i,iy
    real(dp) :: sw

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    sw = 0.0D0
    ix = 1
    iy = 1
    do i = 1,n
       sw = sw + (sx(ix) * sy(iy))
       ix = ix + incx
       iy = iy + incy
    end do

    sdot = sw

  end function sdot

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function isamax(n,sx,incx)

    !! Routine to finds the index of the array element having
    !! the maximum absolute value
    !! author: Jack Dongarra, Linpack
    !! author: P J Knight, CCFE, Culham Science Centre
    !! n        : input integer : order of the matrix sx
    !! sx(n*incx) : input real array : array being checked
    !! incx     : input integer : interval in storage between sx array elements
    !! This routine finds the array element with the maximum
    !! absolute value, and returns the element index.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    integer :: isamax

    !  Arguments

    integer, intent(in) :: n, incx
    real(dp), dimension(n*incx), intent(in) :: sx

    !  Local variables

    integer :: i,ix
    real(dp) :: smax

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    isamax = 0
    if (n < 1) return

    isamax = 1
    if (n == 1) return

    if (incx /= 1) then

       ix = 1
       if (incx < 0) ix = (-n+1)*incx + 1
       smax = abs(sx(ix))
       ix = ix + incx
       do i = 2,n
          if (abs(sx(ix)) > smax) then
             isamax = i
             smax = abs(sx(ix))
          end if
          ix = ix + incx
       end do

    else

       smax = abs(sx(1))
       do i = 2,n
          if (abs(sx(i)) <= smax) cycle
          isamax = i
          smax = abs(sx(i))
       end do

    end if

  end function isamax

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine svd(nm,m,n,a,w,matu,u,matv,v,ierr,rv1)

    !! Singular Value Decomposition
    !! author: P J Knight, CCFE, Culham Science Centre
    !! author: B. S. Garbow, Applied Mathematics Division, Argonne National Laboratory
    !! nm : input integer : Max number of rows of arrays a, u, v; >= m,n
    !! m : input integer : Actual number of rows of arrays a, u
    !! n : input integer : Number of columns of arrays a, u, and the order of v
    !! a(nm,n) : input/output real array : On input matrix to be decomposed;
    !! on output, either unchanged or overwritten with u or v
    !! w(n) : output real array : The n (non-negative) singular values of a
    !! (the diagonal elements of s); unordered.  If an error exit
    !! is made, the singular values should be correct for indices
    !! ierr+1,ierr+2,...,n.
    !! matu : input logical : Set to .true. if the u matrix in the
    !! decomposition is desired, and to .false. otherwise.
    !! u(nm,n) : output real array : The matrix u (orthogonal column vectors)
    !! of the decomposition if matu has been set to .true., otherwise
    !! u is used as a temporary array.  u may coincide with a.
    !! If an error exit is made, the columns of u corresponding
    !! to indices of correct singular values should be correct.
    !! matv : input logical : Set to .true. if the v matrix in the
    !! decomposition is desired, and to .false. otherwise.
    !! v(nm,n) : output real array : The matrix v (orthogonal) of the
    !! decomposition if matv has been set to .true., otherwise
    !! v is not referenced.  v may also coincide with a if u is
    !! not needed.  If an error exit is made, the columns of v
    !! corresponding to indices of correct singular values
    !! should be correct.
    !! ierr : output integer :  zero for normal return, or <I>k</I> if the
    !! k-th singular value has not been determined after 30 iterations.
    !! rv1(n) : output real array : work array
    !! This subroutine is a translation of the algol procedure SVD,
    !! Num. Math. 14, 403-420(1970) by Golub and Reinsch,
    !! Handbook for Auto. Comp., vol II - Linear Algebra, 134-151(1971).
    !! <P>It determines the singular value decomposition
    !! <I>a=usv<SUP>t</SUP></I> of a real m by n rectangular matrix.
    !! Householder bidiagonalization and a variant of the QR
    !! algorithm are used.
    !! None
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    !  Arguments

    integer, intent(in) :: nm, m, n
    logical, intent(in) :: matu, matv
    real(dp), dimension(nm,n), intent(inout) :: a
    real(dp), dimension(nm,n), intent(out) :: u, v
    real(dp), dimension(n), intent(out) :: w, rv1
    integer, intent(out) :: ierr

    !  Local variables

    integer :: i,j,k,l,ii,i1,kk,k1,ll,l1,mn,its
    real(dp) :: c,f,g,h,s,x,y,z,scale,anorm

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ierr = 0

    u = a

    !  Householder reduction to bidiagonal form

    g = 0.0D0
    scale = 0.0D0
    anorm = 0.0D0

    do i = 1, n

       l = i + 1
       rv1(i) = scale * g
       g = 0.0D0
       s = 0.0D0
       scale = 0.0D0

       if (i <= m) then

          do k = i, m
             scale = scale + abs(u(k,i))
          end do

          if (scale /= 0.0D0) then

             do k = i, m
                u(k,i) = u(k,i) / scale
                s = s + u(k,i)**2
             end do

             f = u(i,i)
             g = -sign(sqrt(s),f)
             h = f * g - s
             u(i,i) = f - g

             if (i /= n) then
                do j = l, n
                   s = 0.0D0
                   do k = i, m
                      s = s + u(k,i) * u(k,j)
                   end do
                   f = s / h
                   do k = i, m
                      u(k,j) = u(k,j) + f * u(k,i)
                   end do
                end do
             end if

             do k = i, m
                u(k,i) = scale * u(k,i)
             end do

          end if

       end if

       w(i) = scale * g
       g = 0.0D0
       s = 0.0D0
       scale = 0.0D0

       if (.not.((i > m) .or. (i == n))) then

          do k = l, n
             scale = scale + abs(u(i,k))
          end do

          if (scale /= 0.0D0) then

             do k = l, n
                u(i,k) = u(i,k) / scale
                s = s + u(i,k)**2
             end do

             f = u(i,l)
             g = -sign(sqrt(s),f)
             h = f * g - s
             u(i,l) = f - g

             do k = l, n
                rv1(k) = u(i,k) / h
             end do

             if (i /= m) then
                do j = l, m
                   s = 0.0D0
                   do k = l, n
                      s = s + u(j,k) * u(i,k)
                   end do
                   do k = l, n
                      u(j,k) = u(j,k) + s * rv1(k)
                   end do
                end do
             end if

             do k = l, n
                u(i,k) = scale * u(i,k)
             end do

          end if

       end if

       anorm = max(anorm,abs(w(i))+abs(rv1(i)))

    end do  ! i

    !  Accumulation of right-hand transformations

    if (matv) then

       !  For i=n step -1 until 1 do
       do ii = 1, n
          i = n + 1 - ii
          if (i /= n) then

             if (g /= 0.0D0) then
                do j = l, n
                   !  Double division avoids possible underflow
                   v(j,i) = (u(i,j) / u(i,l)) / g
                end do
                do j = l, n
                   s = 0.0D0
                   do k = l, n
                      s = s + u(i,k) * v(k,j)
                   end do
                   do k = l, n
                      v(k,j) = v(k,j) + s * v(k,i)
                   end do
                end do
             end if

             do j = l, n
                v(i,j) = 0.0D0
                v(j,i) = 0.0D0
             end do

          end if

          v(i,i) = 1.0D0
          g = rv1(i)
          l = i
       end do

    end if

    !  Accumulation of left-hand transformations

    if (matu) then

       !  For i=min(m,n) step -1 until 1 do
       mn = n
       if (m < n) mn = m

       do ii = 1, mn
          i = mn + 1 - ii
          l = i + 1
          g = w(i)
          if (i /= n) then
             do j = l, n
                u(i,j) = 0.0D0
             end do
          end if

          if (g /= 0.0D0) then

             if (i /= mn) then
                do j = l, n
                   s = 0.0D0
                   do k = l, m
                      s = s + u(k,i) * u(k,j)
                   end do
                   f = (s / u(i,i)) / g  !  Double division avoids possible underflow
                   do k = i, m
                      u(k,j) = u(k,j) + f * u(k,i)
                   end do
                end do
             end if

             do j = i, m
                u(j,i) = u(j,i) / g
             end do

          else
             do j = i, m
                u(j,i) = 0.0D0
             end do
          end if

          u(i,i) = u(i,i) + 1.0D0

       end do

    end if

    !  Diagonalization of the bidiagonal form
    !  For k=n step -1 until 1 do

    do kk = 1, n
       k1 = n - kk
       k = k1 + 1
       its = 0

       !  Test for splitting.
       !  For l=k step -1 until 1 do

       do
          do ll = 1, k
             l1 = k - ll
             l = l1 + 1
             if ((abs(rv1(l)) + anorm) == anorm) goto 470

             !  rv1(1) is always zero, so there is no exit
             !  through the bottom of the loop

             !+**PJK 23/05/06 Prevent problems from the code getting here with l1=0
             if (l1 == 0) then
                write(*,*) 'SVD: Shouldn''t get here...'
                goto 470
             end if

             if ((abs(w(l1)) + anorm) == anorm) exit
          end do

          !  Cancellation of rv1(l) if l greater than 1

          c = 0.0D0
          s = 1.0D0

          do i = l, k
             f = s * rv1(i)
             rv1(i) = c * rv1(i)
             if ((abs(f) + anorm) == anorm) exit
             g = w(i)
             h = sqrt(f*f+g*g)
             w(i) = h
             c = g / h
             s = -f / h
             if (.not. matu) cycle

             do j = 1, m
                y = u(j,l1)
                z = u(j,i)
                u(j,l1) = y * c + z * s
                u(j,i) = -y * s + z * c
             end do
          end do

470       continue

          !  Test for convergence

          z = w(k)
          if (l == k) exit

          !  Shift from bottom 2 by 2 minor

          if (its == 30) then
             !  Set error - no convergence to a
             !  singular value after 30 iterations
             ierr = k
             return
          end if

          its = its + 1
          x = w(l)
          y = w(k1)
          g = rv1(k1)
          h = rv1(k)
          f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.D0 * h * y)
          g = sqrt(f*f+1.D0)
          f = ((x - z) * (x + z) + h * (y / (f + sign(g,f)) - h)) / x

          !  Next QR transformation

          c = 1.0D0
          s = 1.0D0

          do i1 = l, k1
             i = i1 + 1
             g = rv1(i)
             y = w(i)
             h = s * g
             g = c * g
             z = sqrt(f*f+h*h)
             rv1(i1) = z
             c = f / z
             s = h / z
             f = x * c + g * s
             g = -x * s + g * c
             h = y * s
             y = y * c

             if (matv) then
                do j = 1, n
                   x = v(j,i1)
                   z = v(j,i)
                   v(j,i1) = x * c + z * s
                   v(j,i) = -x * s + z * c
                end do
             end if

             z = sqrt(f*f+h*h)
             w(i1) = z

             !  Rotation can be arbitrary if z is zero

             if (z /= 0.0D0) then
                c = f / z
                s = h / z
             end if

             f = c * g + s * y
             x = -s * g + c * y
             if (.not. matu) cycle

             do j = 1, m
                y = u(j,i1)
                z = u(j,i)
                u(j,i1) = y * c + z * s
                u(j,i) = -y * s + z * c
             end do

          end do

          rv1(l) = 0.0D0
          rv1(k) = f
          w(k) = x

       end do

       !  Convergence

       if (z >= 0.0D0) cycle

       !  w(k) is made non-negative
       w(k) = -z
       if (.not. matv) cycle

       do j = 1, n
          v(j,k) = -v(j,k)
       end do

    end do

  end subroutine svd

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! hybrd() has been temporarily commented out. Please see the comment in
! function_evaluator.fcnhyb() for an explanation.

!   SUBROUTINE HYBRD( &
!        fcnhyb,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,diag, &
!        mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr, &
!        qtf,wa1,wa2,wa3,wa4,resdl)

!     !  www.math.utah.edu/software/minpack/minpack/hybrd.html

!     !  The purpose of HYBRD is to find a zero of a system of
!     !  N nonlinear functions in N variables by a modification
!     !  of the Powell Hybrid method. The user must provide a
!     !  subroutine which calculates the functions. The Jacobian is
!     !  then calculated by a forward-difference approximation.
!     !
!     !  The subroutine statement is
!     !
!     !  subroutine hybrd(fcnhyb,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,
!     !                   diag,mode,factor,nprint,info,nfev,fjac,
!     !                   ldfjac,r,lr,qtf,wa1,wa2,wa3,wa4)
!     !
!     !  where
!     !
!     !  FCNHYB is the name of the user-supplied subroutine which
!     !  calculates the functions. FCNHYB must be declared
!     !  in an external statement in the user calling
!     !  program, and should be written as follows.
!     !
!     !   subroutine fcnhyb(n,x,fvec,iflag)
!     !   integer n,iflag
!     !   real x(n),fvec(n)
!     !   ----------
!     !   calculate the functions at x and
!     !   return this vector in fvec.
!     !   ---------
!     !   return
!     !   end
!     !
!     !  The value of IFLAG should not be changed by FCNHYB unless
!     !  the user wants to terminate execution of HYBRD.
!     !  In this case set IFLAG to a negative integer.
!     !
!     !  N is a positive integer input variable set to the number
!     !  of functions and variables.
!     !
!     !  X is an array of length N. On input X must contain
!     !  an initial estimate of the solution vector. On output X
!     !  contains the final estimate of the solution vector.
!     !
!     !  FVEC is an output array of length N which contains
!     !  the functions evaluated at the output X.
!     !
!     !  XTOL is a nonnegative input variable. Termination
!     !  occurs when the relative error between two consecutive
!     !  iterations is at most XTOL.
!     !
!     !  MAXFEV is a positive integer input variable. Termination
!     !  occurs when the number of calls to FCNHYB is at least MAXFEV
!     !  by the end of an iteration.
!     !
!     !  ML is a nonnegative integer input variable which specifies
!     !  the number of subdiagonals within the band of the
!     !  Jacobian matrix. If the Jacobian is not banded, set
!     !  ML to at least N - 1.
!     !
!     !  MU is a nonnegative integer input variable which specifies
!     !  the number of superdiagonals within the band of the
!     !  Jacobian matrix. If the Jacobian is not banded, set
!     !  MU to at least N - 1.
!     !
!     !  EPSFCN is an input variable used in determining a suitable
!     !  step length for the forward-difference approximation. This
!     !  approximation assumes that the relative errors in the
!     !  functions are of the order of EPSFCN. If EPSFCN is less
!     !  than the machine precision, it is assumed that the relative
!     !  errors in the functions are of the order of the machine
!     !  precision.
!     !
!     !  DIAG is an array of length N. If MODE = 1 (see
!     !  below), DIAG is internally set. If MODE = 2, DIAG
!     !  must contain positive entries that serve as
!     !  multiplicative scale factors for the variables.
!     !
!     !  MODE is an integer input variable. If MODE = 1, the
!     !  variables will be scaled internally. If MODE = 2,
!     !  the scaling is specified by the input DIAG. Other
!     !  values of MODE are equivalent to MODE = 1.
!     !
!     !  FACTOR is a positive input variable used in determining the
!     !  initial step bound. This bound is set to the product of
!     !  FACTOR and the Euclidean norm of DIAG*X if nonzero, or else
!     !  to FACTOR itself. In most cases FACTOR should lie in the
!     !  interval (.1,100.). 100. is a generally recommended value.
!     !
!     !  NPRINT is an integer input variable that enables controlled
!     !  printing of iterations if it is positive. In this case,
!     !  FCNHYB is called with IFLAG = 0 at the beginning of the first
!     !  iteration and every NPRINT iterations thereafter and
!     !  immediately prior to return, with X and FVEC available
!     !  for printing. If NPRINT is not positive, no special calls
!     !  of FCNHYB with IFLAG = 0 are made.
!     !
!     !  INFO is an integer output variable. If the user has
!     !  terminated execution, INFO is set to the (negative)
!     !  value of IFLAG. see description of FCNHYB. Otherwise,
!     !  INFO is set as follows.
!     !
!     !   INFO = 0   improper input parameters.
!     !
!     !   INFO = 1   relative error between two consecutive iterates
!     !              is at most XTOL.
!     !
!     !   INFO = 2   number of calls to FCNHYB has reached or exceeded
!     !              MAXFEV.
!     !
!     !   INFO = 3   XTOL is too small. No further improvement in
!     !              the approximate solution X is possible.
!     !
!     !   INFO = 4   iteration is not making good progress, as
!     !              measured by the improvement from the last
!     !              five Jacobian evaluations.
!     !
!     !   INFO = 5   iteration is not making good progress, as
!     !              measured by the improvement from the last
!     !              ten iterations.
!     !
!     !  NFEV is an integer output variable set to the number of
!     !  calls to FCNHYB.
!     !
!     !  FJAC is an output N by N array which contains the
!     !  orthogonal matrix Q produced by the QR factorization
!     !  of the final approximate Jacobian.
!     !
!     !  LDFJAC is a positive integer input variable not less than N
!     !  which specifies the leading dimension of the array FJAC.
!     !
!     !  R is an output array of length LR which contains the
!     !  upper triangular matrix produced by the QR factorization
!     !  of the final approximate Jacobian, stored rowwise.
!     !
!     !  LR is a positive integer input variable not less than
!     !  (N*(N+1))/2.
!     !
!     !  QTF is an output array of length N which contains
!     !  the vector (Q transpose)*FVEC.
!     !
!     !  WA1, WA2, WA3, and WA4 are work arrays of length N.
!     !
!     !  Subprograms called
!     !
!     !   user-supplied ...... fcnhyb
!     !
!     !   minpack-supplied ... dogleg,spmpar,enorm,fdjac1,
!     !                        qform,qrfac,r1mpyq,r1updt
!     !
!     !  Argonne National Laboratory. Minpack project. March 1980.
!     !  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

!     IMPLICIT NONE

!     interface
!       subroutine fcnhyb(n, x, fvec, iflag)
!         use, intrinsic :: iso_fortran_env, only: dp=>real64
!         integer, intent(in) :: n
!         real(dp), dimension(n), intent(inout) :: x
!         real(dp), dimension(n), intent(out) :: fvec
!         integer, intent(inout) :: iflag
!       end subroutine fcnhyb
!     end interface

!     INTEGER n,maxfev,ml,mu,mode,nprint,info,nfev,ldfjac,lr,irr
!     INTEGER i,iflag,iter,j,jm1,l,msum,ncfail,ncsuc,nslow1,nslow2

!     !+**PJK 08/10/92 Possible problems with the following declaration:
!     INTEGER iwa(1)

!     real(dp) xtol,epsfcn,factor
!     real(dp) x(n),fvec(n),diag(n),fjac(ldfjac,n),r(lr), &
!          qtf(n),wa1(n),wa2(n),wa3(n),wa4(n),resdl(n)
!     real(dp) actred,delta,epsmch,fnorm,fnorm1,one,pnorm, &
!          prered,p1,p5,p001,p0001,ratio,sum,temp,xnorm,zero
!     logical jeval,sing

!     EXTERNAL fcnhyb

!     one = 1.0D0
!     p1 = 0.1D0
!     p5 = 0.5D0
!     p001 = 1.0D-3
!     p0001 = 1.0D-4
!     zero = 0.0D0

!     !  Machine precision

!     epsmch = spmpar(1)

!     info = 0
!     iflag = 0
!     nfev = 0

!     !  Check the input parameters for errors.

!     if ( &
!          (n <= 0)         .or. &
!          (xtol < zero)   .or. &
!          (maxfev <= 0)    .or. &
!          (ml < 0)        .or. &
!          (mu < 0)        .or. &
!          (factor <= zero) .or. &
!          (ldfjac < n)    .or. &
!          (lr < ( ( n*(n + 1) ) /2)) &
!          ) goto 300

!     if (mode  /=  2) goto 20
!     do j = 1, n
!        if (diag(j) <= zero) goto 300
!     end do

! 20  continue

!     !  Evaluate the function at the starting point
!     !  and calculate its norm.

!     iflag = 1
!     call fcnhyb(n,x,fvec,iflag)
!     nfev = 1

!     if (iflag < 0) goto 300
!     fnorm = enorm(n,fvec)

!     !  Determine the number of calls to FCNHYB needed to compute
!     !  the Jacobian matrix.

!     msum = min(ml+mu+1,n)

!     !  Initialize iteration counter and monitors.

!     iter = 1
!     ncsuc = 0
!     ncfail = 0
!     nslow1 = 0
!     nslow2 = 0

!     !  Beginning of the outer loop.

! 30  continue
!     jeval = .true.

!     !  Calculate the Jacobian matrix.

!     iflag = 2
!     call fdjac1( &
!          fcnhyb,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)
!     nfev = nfev + msum
!     if (iflag < 0) goto 300

!     !  Compute the qr factorization of the Jacobian.

!     call qrfac(n,n,fjac,ldfjac,.false.,iwa,1,wa1,wa2,wa3)

!     !  On the first iteration and if mode is 1, scale according
!     !  to the norms of the columns of the initial Jacobian.

!     if (iter  /=  1) goto 70
!     if (mode == 2) goto 50
!     do j = 1, n
!        diag(j) = wa2(j)
!        if (wa2(j) == zero) diag(j) = one
!     end do

! 50  continue

!     !  On the first iteration, calculate the norm of the scaled x
!     !  and initialize the step bound delta.

!     do j = 1, n
!        wa3(j) = diag(j)*x(j)
!     end do
!     xnorm = enorm(n,wa3)
!     delta = factor*xnorm
!     if (delta == zero) delta = factor

! 70  continue

!     !  Form (q transpose)*fvec and store in qtf.

!     do i = 1, n
!        qtf(i) = fvec(i)
!     end do
!     do j = 1, n
!        if (fjac(j,j) == zero) goto 110
!        sum = zero
!        do i = j, n
!           sum = sum + fjac(i,j)*qtf(i)
!        end do
!        temp = -sum/fjac(j,j)
!        do i = j, n
!           qtf(i) = qtf(i) + fjac(i,j)*temp
!        end do

! 110    continue
!     end do

!     !  Copy the triangular factor of the qr factorization into r.

!     sing = .false.
!     do j = 1, n
!        l = j
!        jm1 = j - 1
!        if (jm1 < 1) goto 140
!        do i = 1, jm1
!           r(l) = fjac(i,j)
!           l = l + n - i
!        end do

! 140    continue
!        r(l) = wa1(j)
!        if (wa1(j) == zero) sing = .true.
!     end do

!     !  Accumulate the orthogonal factor in fjac.

!     call qform(n,n,fjac,ldfjac,wa1)

!     !  Rescale if necessary.

!     if (mode == 2) goto 170
!     do j = 1, n
!        diag(j) = max(diag(j),wa2(j))
!     end do

! 170 continue

!     !  Beginning of the inner loop.

! 180 continue

!     !  If requested, call FCNHYB to enable printing of iterates.

!     if (nprint <= 0) goto 190
!     iflag = 0
!     if (mod(iter-1,nprint) == 0) call fcnhyb(n,x,fvec,iflag)
!     if (iflag < 0) goto 300

! 190 continue

!     !  Determine the direction p.

!     call dogleg(n,r,lr,diag,qtf,delta,wa1,wa2,wa3)

!     !  Store the direction p and x + p. Calculate the norm of p.

!     do j = 1, n
!        wa1(j) = -wa1(j)
!        wa2(j) = x(j) + wa1(j)
!        wa3(j) = diag(j)*wa1(j)
!     end do
!     pnorm = enorm(n,wa3)

!     !  On the first iteration, adjust the initial step bound.

!     if (iter == 1) delta = min(delta,pnorm)

!     !  Evaluate the function at x + p and calculate its norm.

!     iflag = 1
!     call fcnhyb(n,wa2,wa4,iflag)
!     nfev = nfev + 1
!     if (iflag < 0) goto 300
!     fnorm1 = enorm(n,wa4)

!     !  Compute the scaled actual reduction.

!     actred = -one
!     if (fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

!     !  Compute the scaled predicted reduction.

!     l = 1
!     do i = 1, n
!        sum = zero
!        do j = i, n
!           sum = sum + r(l)*wa1(j)
!           l = l + 1
!        end do
!        wa3(i) = qtf(i) + sum
!     end do
!     temp = enorm(n,wa3)
!     prered = zero
!     if (temp < fnorm) prered = one - (temp/fnorm)**2

!     !  Compute the ratio of the actual to the predicted reduction.

!     ratio = zero
!     if (prered > zero) ratio = actred/prered

!     !  Update the step bound.

!     if (ratio >= p1) goto 230
!     ncsuc = 0
!     ncfail = ncfail + 1
!     delta = p5*delta
!     goto 240
! 230 continue
!     ncfail = 0
!     ncsuc = ncsuc + 1
!     if (ratio >= p5 .or. ncsuc > 1) &
!          delta = max(delta,pnorm/p5)
!     if (abs(ratio-one) <= p1) delta = pnorm/p5
! 240 continue

!     !  Test for successful iteration.

!     if (ratio < p0001) goto 260

!     !  Successful iteration. Update x, fvec, and their norms.

!     do j = 1, n
!        x(j) = wa2(j)
!        wa2(j) = diag(j)*x(j)
!        fvec(j) = wa4(j)
!     end do
!     xnorm = enorm(n,wa2)
!     fnorm = fnorm1
!     iter = iter + 1
! 260 continue

!     !  Determine the progress of the iteration.

!     nslow1 = nslow1 + 1
!     if (actred >= p001) nslow1 = 0
!     if (jeval) nslow2 = nslow2 + 1
!     if (actred >= p1) nslow2 = 0

!     !  Test for convergence.

!     if ((delta <= (xtol*xnorm)) .or. (fnorm == zero)) info = 1
!     if (info  /=  0) goto 300

!     !  Tests for termination and stringent tolerances.

!     if (nfev >= maxfev) info = 2
!     if ((p1*max(p1*delta,pnorm)) <= (epsmch*xnorm)) info = 3
!     if (nslow2 == 5) info = 4
!     if (nslow1 == 10) info = 5
!     if (info  /=  0) goto 300

!     !  Criterion for recalculating Jacobian approximation
!     !  by forward differences.

!     if (ncfail == 2) goto 290

!     !  Calculate the rank one modification to the Jacobian
!     !  and update qtf if necessary.

!     do j = 1, n
!        sum = zero
!        do i = 1, n
!           sum = sum + fjac(i,j)*wa4(i)
!        end do
!        wa2(j) = (sum - wa3(j))/pnorm
!        wa1(j) = diag(j)*((diag(j)*wa1(j))/pnorm)
!        if (ratio >= p0001) qtf(j) = sum
!     end do

!     !  Compute the qr factorization of the updated Jacobian.

!     call r1updt(n,n,r,lr,wa1,wa2,wa3,sing)
!     call r1mpyq(n,n,fjac,ldfjac,wa2,wa3)

!     !+**PJK 02/11/92 Warning produced by QA Fortran :
!     !+**PJK 02/11/92 Arg 3 in call to R1MPYQ has wrong dimensions.
!     !+**PJK 02/11/92 Code works at present, but beware of future
!     !+**PJK 02/11/92 modifications.

!     call r1mpyq(1,n,qtf,1,wa2,wa3)

!     !  End of the inner loop.

!     jeval = .false.
!     goto 180

! 290 continue

!     !  End of the outer loop.

!     goto 30

! 300 continue

!     !  Termination, either normal or user imposed.

!     if (iflag < 0) info = iflag
!     iflag = 0
!     if (nprint > 0) call fcnhyb(n,x,fvec,iflag)

!     do irr=1,n
!        resdl(irr)=abs(qtf(irr))
!     end do

!     return
!   end SUBROUTINE HYBRD

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE DOGLEG(n,r,lr,diag,qtb,delta,x,wa1,wa2)

    !  Given an M by N matrix A, an N by N nonsingular diagonal
    !  matrix D, an M-vector B, and a positive number DELTA, the
    !  problem is to determine the convex combination X of the
    !  Gauss-Newton and scaled gradient directions that minimizes
    !  (A*X - B) in the least squares sense, subject to the
    !  restriction that the Euclidean norm of D*X be at most DELTA.
    !
    !  This subroutine completes the solution of the problem
    !  if it is provided with the necessary information from the
    !  QR factorization of A. That is, if A = Q*R, where Q has
    !  orthogonal columns and R is an upper triangular matrix,
    !  then DOGLEG expects the full upper triangle of R and
    !  the first N components of (Q transpose)*B.
    !
    !  The subroutine statement is
    !
    !  subroutine dogleg(n,r,lr,diag,qtb,delta,x,wa1,wa2)
    !
    !  where
    !
    !  N is a positive integer input variable set to the order of R.
    !
    !  R is an input array of length LR which must contain the upper
    !  triangular matrix R stored by rows.
    !
    !  LR is a positive integer input variable not less than
    !  (N*(N+1))/2.
    !
    !  DIAG is an input array of length N which must contain the
    !  diagonal elements of the matrix D.
    !
    !  QTB is an input array of length N which must contain the first
    !  N elements of the vector (Q transpose)*B.
    !
    !  DELTA is a positive input variable which specifies an upper
    !  bound on the Euclidean norm of D*X.
    !
    !  X is an output array of length N which contains the desired
    !  convex combination of the Gauss-Newton direction and the
    !  scaled gradient direction.
    !
    !  WA1 and WA2 are work arrays of length N.
    !
    !  Argonne National Laboratory. Minpack project. March 1980.
    !  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    IMPLICIT NONE

    INTEGER n,lr,i,j,jj,jp1,k,l

    real(dp) delta
    real(dp) r(lr),diag(n),qtb(n),x(n),wa1(n),wa2(n)
    real(dp) alpha,bnorm,epsmch,gnorm,one,qnorm,sgnorm, &
         sum,temp,zero

    one = 1.0D0
    zero = 0.0D0

    !  Machine precision

    epsmch = spmpar(1)

    !  First, calculate the Gauss-Newton direction.

    jj = (n*(n + 1))/2 + 1
    do k = 1, n
       j = n - k + 1
       jp1 = j + 1
       jj = jj - k
       l = jj + 1
       sum = zero
       if (n < jp1) goto 20
       do i = jp1, n
          sum = sum + r(l)*x(i)
          l = l + 1
       end do

20     continue
       temp = r(jj)
       if (temp  /=  zero) goto 40
       l = j
       do i = 1, j
          temp = max(temp,abs(r(l)))
          l = l + n - i
       end do
       temp = epsmch*temp
       if (temp == zero) temp = epsmch

40     continue
       x(j) = (qtb(j) - sum)/temp
    end do

    !  Test whether the Gauss-Newton direction is acceptable.

    do j = 1, n
       wa1(j) = zero
       wa2(j) = diag(j)*x(j)
    end do
    qnorm = enorm(n,wa2)
    if (qnorm <= delta) goto 140

    !  The Gauss-Newton direction is not acceptable.
    !  Next, calculate the scaled gradient direction.

    l = 1
    do j = 1, n
       temp = qtb(j)
       do i = j, n
          wa1(i) = wa1(i) + r(l)*temp
          l = l + 1
       end do
       wa1(j) = wa1(j)/diag(j)
    end do

    !  Calculate the norm of the scaled gradient and test for
    !  the special case in which the scaled gradient is zero.

    gnorm = enorm(n,wa1)
    sgnorm = zero
    alpha = delta/qnorm
    if (gnorm == zero) goto 120

    !  Calculate the point along the scaled gradient
    !  at which the quadratic is minimized.

    do j = 1, n
       wa1(j) = (wa1(j)/gnorm)/diag(j)
    end do
    l = 1
    do j = 1, n
       sum = zero
       do i = j, n
          sum = sum + r(l)*wa1(i)
          l = l + 1
       end do
       wa2(j) = sum
    end do
    temp = enorm(n,wa2)
    sgnorm = (gnorm/temp)/temp

    !  Test whether the scaled gradient direction is acceptable.

    alpha = zero
    if (sgnorm >= delta) goto 120

    !  The scaled gradient direction is not acceptable.
    !  Finally, calculate the point along the dogleg
    !  at which the quadratic is minimized.

    bnorm = enorm(n,qtb)
    temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta)
    temp = temp - (delta/qnorm)*(sgnorm/delta)**2 &
         + sqrt((temp-(delta/qnorm))**2 &
         + (one-(delta/qnorm)**2)*(one-(sgnorm/delta)**2))
    alpha = ((delta/qnorm)*(one - (sgnorm/delta)**2))/temp

120 continue

    !  Form appropriate convex combination of the Gauss-Newton
    !  direction and the scaled gradient direction.

    temp = (one - alpha)*min(sgnorm,delta)
    do j = 1, n
       x(j) = temp*wa1(j) + alpha*x(j)
    end do

140 continue

    return
  end SUBROUTINE DOGLEG

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp) FUNCTION ENORM(n,x)

    !  Given an N-vector X, this function calculates the
    !  Euclidean norm of X.
    !
    !  The Euclidean norm is computed by accumulating the sum of
    !  squares in three different sums. The sums of squares for the
    !  small and large components are scaled so that no overflows
    !  occur. Non-destructive underflows are permitted. Underflows
    !  and overflows do not occur in the computation of the unscaled
    !  sum of squares for the intermediate components.
    !  The definitions of small, intermediate and large components
    !  depend on two constants, RDWARF and RGIANT. The main
    !  restrictions on these constants are that RDWARF**2 not
    !  underflow and RGIANT**2 not overflow. The constants
    !  given here are suitable for every known computer.
    !
    !  The function statement is
    !
    !  real function enorm(n,x)
    !
    !  where
    !
    !  N is a positive integer input variable.
    !
    !  X is an input array of length N.
    !
    !  Argonne National Laboratory. Minpack project. March 1980.
    !  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    IMPLICIT NONE

    INTEGER n,i

    real(dp) x(n)
    real(dp) agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs, &
         x1max,x3max,zero

    one = 1.0D0
    zero = 0.0D0
    rdwarf = 3.834d-20
    rgiant = 1.304d19

    s1 = zero
    s2 = zero
    s3 = zero
    x1max = zero
    x3max = zero
    floatn = dble(n)
    agiant = rgiant/floatn
    do i = 1, n
       xabs = abs(x(i))
       if ((xabs > rdwarf) .and. (xabs < agiant)) goto 70
       if (xabs <= rdwarf) goto 30

       !  Sum for large components.

       if (xabs <= x1max) goto 10
       s1 = one + s1*(x1max/xabs)**2
       x1max = xabs
       goto 20

10     continue
       s1 = s1 + (xabs/x1max)**2

20     continue
       goto 60

30     continue

       !  Sum for small components.

       if (xabs <= x3max) goto 40
       s3 = one + s3*(x3max/xabs)**2
       x3max = xabs
       goto 50

40     continue
       if (xabs  /=  zero) s3 = s3 + (xabs/x3max)**2

50     continue
60     continue
       goto 80

70     continue

       !  Sum for intermediate components.

       s2 = s2 + xabs**2

80     continue
    end do

    !  Calculation of norm.

    if (s1 == zero) goto 100
    enorm = x1max*sqrt(s1+(s2/x1max)/x1max)
    goto 130

100 continue
    if (s2 == zero) goto 110
    if (s2 >= x3max) &
         enorm = sqrt(s2*(one+(x3max/s2)*(x3max*s3)))
    if (s2 < x3max) &
         enorm = sqrt(x3max*((s2/x3max)+(x3max*s3)))
    goto 120
110 continue
    enorm = x3max*sqrt(s3)

120 continue
130 continue

    return
  end FUNCTION ENORM

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE FDJAC1( &
       fcnhyb,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,wa1,wa2)

    !  This subroutine computes a forward-difference approximation
    !  to the N by N Jacobian matrix associated with a specified
    !  problem of N functions in N variables. If the Jacobian has
    !  a banded form, then function evaluations are saved by only
    !  approximating the nonzero terms.
    !
    !  The subroutine statement is
    !
    !  subroutine fdjac1(fcnhyb,n,x,fvec,fjac,ldfjac,iflag,ml,mu,epsfcn,
    !                    wa1,wa2)
    !
    !  where
    !
    !  FCNHYB is the name of the user-supplied subroutine which
    !  calculates the functions. FCNHYB must be declared
    !  in an external statement in the user calling
    !  program, and should be written as follows.
    !
    !   subroutine fcnhyb(n,x,fvec,iflag)
    !   integer n,iflag
    !   real x(n),fvec(n)
    !   ----------
    !   calculate the functions at x and
    !   return this vector in fvec.
    !   ----------
    !   return
    !   end
    !
    !  The value of IFLAG should not be changed by FCNHYB unless
    !  the user wants to terminate execution of FDJAC1.
    !  In this case set IFLAG to a negative integer.
    !
    !  N is a positive integer input variable set to the number
    !  of functions and variables.
    !
    !  X is an input array of length N.
    !
    !  FVEC is an input array of length N which must contain the
    !  functions evaluated at X.
    !
    !  FJAC is an output N by N array which contains the
    !  approximation to the Jacobian matrix evaluated at X.
    !
    !  LDFJAC is a positive integer input variable not less than N
    !  which specifies the leading dimension of the array FJAC.
    !
    !  IFLAG is an integer variable which can be used to terminate
    !  the execution of FDJAC1. See description of FCNHYB.
    !
    !  ML is a nonnegative integer input variable which specifies
    !  the number of subdiagonals within the band of the
    !  Jacobian matrix. If the Jacobian is not banded, set
    !  ML to at least N - 1.
    !
    !  EPSFCN is an input variable used in determining a suitable
    !  step length for the forward-difference approximation. This
    !  approximation assumes that the relative errors in the
    !  functions are of the order of EPSFCN. If EPSFCN is less
    !  than the machine precision, it is assumed that the relative
    !  errors in the functions are of the order of the machine
    !  precision.
    !
    !  MU is a nonnegative integer input variable which specifies
    !  the number of superdiagonals within the band of the
    !  Jacobian matrix. If the Jacobian is not banded, set
    !  MU to at least N - 1.
    !
    !  WA1 and WA2 are work arrays of length N. If ML + MU + 1 is at
    !  least N, then the Jacobian is considered dense, and WA2 is
    !  not referenced.
    !
    !  Argonne National Laboratory. Minpack project. March 1980.
    !  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    IMPLICIT NONE

    INTEGER n,ldfjac,iflag,ml,mu,i,j,k,msum

    real(dp) epsfcn
    real(dp) x(n),fvec(n),fjac(ldfjac,n),wa1(n),wa2(n)
    real(dp) eps,epsmch,h,temp,zero

    EXTERNAL  fcnhyb

    zero = 0.0D0

    !  Machine precision

    epsmch = spmpar(1)

    eps = sqrt(max(epsfcn,epsmch))
    msum = ml + mu + 1
    if (msum < n) goto 40

    !  Computation of dense approximate Jacobian.

    do j = 1, n
       temp = x(j)
       h = eps*abs(temp)
       if (h == zero) h = eps
       x(j) = temp + h
       call fcnhyb(n,x,wa1,iflag)
       if (iflag < 0) goto 30
       x(j) = temp
       do i = 1, n
          fjac(i,j) = (wa1(i) - fvec(i))/h
       end do
    end do

30  continue
    goto 110

40  continue

    !  Computation of banded approximate Jacobian.

    do k = 1, msum
       do j = k, n, msum
          wa2(j) = x(j)
          h = eps*abs(wa2(j))
          if (h == zero) h = eps
          x(j) = wa2(j) + h
       end do
       call fcnhyb(n,x,wa1,iflag)
       if (iflag < 0) goto 100
       do j = k, n, msum
          x(j) = wa2(j)
          h = eps*abs(wa2(j))
          if (h == zero) h = eps
          do i = 1, n
             fjac(i,j) = zero
             if ((i >= (j - mu)).and.(i <= (j + ml))) &
                  fjac(i,j) = (wa1(i) - fvec(i))/h
          end do
       end do
    end do

100 continue
110 continue

    return
  end SUBROUTINE FDJAC1

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE QFORM(m,n,q,ldq,wa)

    !  This subroutine proceeds from the computed QR factorization of
    !  an M by N matrix A to accumulate the M by M orthogonal matrix
    !  Q from its factored form.
    !
    !  The subroutine statement is
    !
    !  subroutine qform(m,n,q,ldq,wa)
    !
    !  where
    !
    !  M is a positive integer input variable set to the number
    !  of rows of A and the order of Q.
    !
    !  N is a positive integer input variable set to the number
    !  of columns of A.
    !
    !  Q is an M by M array. On input the full lower trapezoid in
    !  the first min(M,N) columns of Q contains the factored form.
    !  On output Q has been accumulated into a square matrix.
    !
    !  LDQ is a positive integer input variable not less than M
    !  which specifies the leading dimension of the array Q.
    !
    !  WA is a work array of length M.

    !  Argonne National Laboratory. Minpack project. March 1980.
    !  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    IMPLICIT NONE

    INTEGER m,n,ldq,i,j,jm1,k,l,minmn,np1

    real(dp) q(ldq,m),wa(m)
    real(dp) one,sum,temp,zero

    one = 1.0D0
    zero = 0.0D0

    !  Zero out upper triangle of q in the first min(m,n) columns.

    minmn = min(m,n)
    if (minmn < 2) goto 30
    do j = 2, minmn
       jm1 = j - 1
       do i = 1, jm1
          q(i,j) = zero
       end do
    end do

30  continue

    !  Initialize remaining columns to those of the identity matrix.

    np1 = n + 1
    if (m < np1) goto 60
    do j = np1, m
       do i = 1, m
          q(i,j) = zero
       end do
       q(j,j) = one
    end do

60  continue

    !  Accumulate q from its factored form.

    do l = 1, minmn
       k = minmn - l + 1
       do i = k, m
          wa(i) = q(i,k)
          q(i,k) = zero
       end do
       q(k,k) = one
       if (wa(k) == zero) goto 110
       do j = k, m
          sum = zero
          do i = k, m
             sum = sum + q(i,j)*wa(i)
          end do
          temp = sum/wa(k)
          do i = k, m
             q(i,j) = q(i,j) - temp*wa(i)
          end do
       end do

110    continue
    end do

    return
  end SUBROUTINE QFORM

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE QRFAC(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)

    !  This subroutine uses householder transformations with column
    !  pivoting (optional) to compute a QR factorization of the
    !  M by N matrix A. That is, QRFAC determines an orthogonal
    !  matrix Q, a permutation matrix P, and an upper trapezoidal
    !  matrix R with diagonal elements of nonincreasing magnitude,
    !  such that A*P = Q*R. The householder transformation for
    !  column K, K = 1,2,...,min(M,N), is of the form
    !
    !  i - (1/u(k))*u*u
    !
    !  where U has zeros in the first K-1 positions. The form of
    !  this transformation and the method of pivoting first
    !  appeared in the corresponding Linpack subroutine.
    !
    !  The subroutine statement is
    !
    !  subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)
    !
    !  where
    !
    !  M is a positive integer input variable set to the number
    !  of rows of A.
    !
    !  N is a positive integer input variable set to the number
    !  of columns of A.
    !
    !  A is an M by N array. On input A contains the matrix for
    !  which the QR factorization is to be computed. On output
    !  the strict upper trapezoidal part of A contains the strict
    !  upper trapezoidal part of R, and the lower trapezoidal
    !  part of A contains a factored form of Q (the non-trivial
    !  elements of the U vectors described above).
    !
    !  LDA is a positive integer input variable not less than M
    !  which specifies the leading dimension of the array A.
    !
    !  PIVOT is a logical input variable. If PIVOT is set true,
    !  then column pivoting is enforced. If PIVOT is set false,
    !  then no column pivoting is done.
    !
    !  IPVT is an integer output array of length LIPVT. IPVT
    !  defines the permutation matrix P such that A*P = Q*R.
    !  Column J of P is column IPVT(J) of the identity matrix.
    !  If PIVOT is false, IPVT is not referenced.
    !
    !  LIPVT is a positive integer input variable. If PIVOT is false,
    !  then LIPVT may be as small as 1. If PIVOT is true, then
    !  LIPVT must be at least N.
    !
    !  RDIAG is an output array of length N which contains the
    !  diagonal elements of R.
    !
    !  ACNORM is an output array of length N which contains the
    !  norms of the corresponding columns of the input matrix A.
    !  If this information is not needed, then ACNORM can coincide
    !  with RDIAG.
    !
    !  WA is a work array of length N. If PIVOT is false, then WA
    !  can coincide with RDIAG.
    !
    !  Argonne National Laboratory. Minpack project. March 1980.
    !  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    IMPLICIT NONE

    INTEGER m,n,lda,lipvt,i,j,jp1,k,kmax,minmn
    INTEGER ipvt(lipvt)

    LOGICAL pivot

    real(dp) a(lda,n),rdiag(n),acnorm(n),wa(n)
    real(dp) ajnorm,epsmch,one,p05,sum,temp,zero

    one = 1.0D0
    p05 = 0.05D0
    zero = 0.0D0

    !  Machine precision

    epsmch = spmpar(1)

    !  Compute the initial column norms and initialize several arrays.

    do j = 1, n
       acnorm(j) = enorm(m,a(1,j))
       rdiag(j) = acnorm(j)
       wa(j) = rdiag(j)
       if (pivot) ipvt(j) = j
    end do

    !  Reduce a to r with householder transformations.

    minmn = min(m,n)
    do j = 1, minmn
       if (.not.pivot) goto 40

       !  Bring the column of largest norm into the pivot position.

       kmax = j
       do k = j, n
          if (rdiag(k) > rdiag(kmax)) kmax = k
       end do
       if (kmax == j) goto 40
       do i = 1, m
          temp = a(i,j)
          a(i,j) = a(i,kmax)
          a(i,kmax) = temp
       end do
       rdiag(kmax) = rdiag(j)
       wa(kmax) = wa(j)
       k = ipvt(j)
       ipvt(j) = ipvt(kmax)
       ipvt(kmax) = k

40     continue

       !  Compute the householder transformation to reduce the
       !  j-th column of a to a multiple of the j-th unit vector.

       ajnorm = enorm(m-j+1,a(j,j))
       if (ajnorm == zero) goto 100
       if (a(j,j) < zero) ajnorm = -ajnorm
       do i = j, m
          a(i,j) = a(i,j)/ajnorm
       end do
       a(j,j) = a(j,j) + one

       !  Apply the transformation to the remaining columns
       !  and update the norms.

       jp1 = j + 1
       if (n < jp1) goto 100
       do k = jp1, n
          sum = zero
          do i = j, m
             sum = sum + a(i,j)*a(i,k)
          end do
          temp = sum/a(j,j)
          do i = j, m
             a(i,k) = a(i,k) - temp*a(i,j)
          end do
          if ((.not.pivot).or.(rdiag(k) == zero)) goto 80
          temp = a(j,k)/rdiag(k)
          rdiag(k) = rdiag(k)*sqrt(max(zero,one-temp**2))
          if ((p05*(rdiag(k)/wa(k))**2) > epsmch) goto 80
          rdiag(k) = enorm(m-j,a(jp1,k))
          wa(k) = rdiag(k)

80        continue
       end do

100    continue
       rdiag(j) = -ajnorm
    end do

    return
  end SUBROUTINE QRFAC

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE R1MPYQ(m,n,a,lda,v,w)

    !  Given an M by N matrix A, this subroutine computes A*Q where
    !  Q is the product of 2*(N - 1) transformations
    !
    !  gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
    !
    !  and GV(I), GW(i) are Givens rotations in the (I,N) plane which
    !  eliminate elements in the I-th and N-th planes, respectively.
    !  Q itself is not given, rather the information to recover the
    !  GV, GW rotations is supplied.
    !
    !  The subroutine statement is
    !
    !  subroutine r1mpyq(m,n,a,lda,v,w)
    !
    !  where
    !
    !  M is a positive integer input variable set to the number
    !  of rows of A.
    !
    !  N is a positive integer input variable set to the number
    !  of columns of A.
    !
    !  A is an M by N array. On input A must contain the matrix
    !  to be postmultiplied by the orthogonal matrix Q
    !  described above. On output A*Q has replaced A.
    !
    !  LDA is a positive integer input variable not less than M
    !  which specifies the leading dimension of the array A.
    !
    !  V is an input array of length N. V(I) must contain the
    !  information necessary to recover the Givens rotation GV(I)
    !  described above.
    !
    !  W is an input array of length N. W(I) must contain the
    !  information necessary to recover the Givens rotation GW(I)
    !  described above.
    !
    !  Argonne National Laboratory. Minpack project. March 1980.
    !  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More

    IMPLICIT NONE

    INTEGER m,n,lda,i,j,nmj,nm1

    real(dp) a(lda,n),v(n),w(n)
    real(dp) cos1,one,sin1,temp

    one = 1.0D0

    !  Apply the first set of givens rotations to a.

    nm1 = n - 1
    if (nm1 < 1) goto 50
    do nmj = 1, nm1
       j = n - nmj
       if (abs(v(j)) > one) cos1 = one/v(j)
       if (abs(v(j)) > one) sin1 = sqrt(one-cos1**2)
       if (abs(v(j)) <= one) sin1 = v(j)
       if (abs(v(j)) <= one) cos1 = sqrt(one-sin1**2)
       do i = 1, m
          temp  = cos1*a(i,j) - sin1*a(i,n)
          a(i,n) = sin1*a(i,j) + cos1*a(i,n)
          a(i,j) = temp
       end do
    end do

    !  Apply the second set of givens rotations to a.

    do j = 1, nm1
       if (abs(w(j)) > one) cos1 = one/w(j)
       if (abs(w(j)) > one) sin1 = sqrt(one-cos1**2)
       if (abs(w(j)) <= one) sin1 = w(j)
       if (abs(w(j)) <= one) cos1 = sqrt(one-sin1**2)
       do i = 1, m
          temp   =  cos1*a(i,j) + sin1*a(i,n)
          a(i,n) = -sin1*a(i,j) + cos1*a(i,n)
          a(i,j) = temp
       end do
    end do

50  continue

    return
  end SUBROUTINE R1MPYQ

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE R1UPDT(m,n,s,ls,u,v,w,sing)

    !  Given an M by N lower trapezoidal matrix S, an M-vector U,
    !  and an N-vector V, the problem is to determine an
    !  orthogonal matrix Q such that
    !
    !          t
    !  (s + u*v )*q
    !
    !  is again lower trapezoidal.
    !
    !  This subroutine determines Q as the product of 2*(N - 1)
    !  transformations
    !
    !  gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
    !
    !  where GV(I), GW(I) are Givens rotations in the (I,N) plane
    !  which eliminate elements in the I-th and N-th planes,
    !  respectively. Q itself is not accumulated, rather the
    !  information to recover the GV, GW rotations is returned.
    !
    !  The subroutine statement is
    !
    !  subroutine r1updt(m,n,s,ls,u,v,w,sing)
    !
    !  where
    !
    !  M is a positive integer input variable set to the number
    !  of rows of S.
    !
    !  N is a positive integer input variable set to the number
    !  of columns of S. N must not exceed M.
    !
    !  S is an array of length LS. On input S must contain the lower
    !  trapezoidal matrix S stored by columns. On output S contains
    !  the lower trapezoidal matrix produced as described above.
    !
    !  LS is a positive integer input variable not less than
    !  (N*(2*M-N+1))/2.
    !
    !  U is an input array of length M which must contain the
    !  vector U.
    !
    !  V is an array of length N. On input V must contain the vector
    !  V. On output V(I) contains the information necessary to
    !  recover the Givens rotation GV(I) described above.
    !
    !  W is an output array of length M. W(I) contains information
    !  necessary to recover the Givens rotation GW(I) described
    !  above.
    !
    !  SING is a logical output variable. SING is set true if any
    !  of the diagonal elements of the output S are zero. Otherwise
    !  SING is set false.
    !
    !  Argonne National Laboratory. Minpack project. March 1980.
    !  Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More,
    !  John L. Nazareth

    IMPLICIT NONE

    INTEGER m,n,ls,i,j,jj,l,nmj,nm1

    LOGICAL sing

    real(dp) s(ls),u(m),v(n),w(m),cos1,cotan,giant,one
    real(dp) p5,p25,sin1,tan1,tau,temp,zero

    one = 1.0D0
    p5 = 0.5D0
    p25 = 0.25D0
    zero = 0.0D0

    !  giant is the largest magnitude in the computer's arithmetic range

    giant = spmpar(3)

    !  Initialize the diagonal element pointer.

    jj = (n*(2*m - n + 1))/2 - (m - n)

    !  Move the nontrivial part of the last column of s into w.

    l = jj
    do i = n, m
       w(i) = s(l)
       l = l + 1
    end do

    !  Rotate the vector v into a multiple of the n-th unit vector
    !  in such a way that a spike is introduced into w.

    nm1 = n - 1
    if (nm1 < 1) goto 70
    do nmj = 1, nm1
       j = n - nmj
       jj = jj - (m - j + 1)
       w(j) = zero
       if (v(j) == zero) goto 50

       !  Determine a givens rotation which eliminates the
       !  j-th element of v.

       if (abs(v(n)) >= abs(v(j))) goto 20
       cotan = v(n)/v(j)
       sin1 = p5/sqrt(p25+p25*cotan**2)
       cos1 = sin1*cotan
       tau = one
       if (abs(cos1)*giant > one) tau = one/cos1
       goto 30

20     continue
       tan1 = v(j)/v(n)
       cos1 = p5/sqrt(p25+p25*tan1**2)
       sin1 = cos1*tan1
       tau = sin1

30     continue

       !  Apply the transformation to v and store the information
       !  necessary to recover the givens rotation.

       v(n) = sin1*v(j) + cos1*v(n)
       v(j) = tau

       !  Apply the transformation to s and extend the spike in w.

       l = jj
       do i = j, m
          temp = cos1*s(l) - sin1*w(i)
          w(i) = sin1*s(l) + cos1*w(i)
          s(l) = temp
          l = l + 1
       end do

50     continue
    end do

70  continue

    !  Add the spike from the rank 1 update to w.

    do i = 1, m
       w(i) = w(i) + v(n)*u(i)
    end do

    !  Eliminate the spike.

    sing = .false.
    if (nm1 < 1) goto 140
    do j = 1, nm1
       if (w(j) == zero) goto 120

       !  Determine a givens rotation which eliminates the
       !  j-th element of the spike.

       if (abs(s(jj)) >= abs(w(j))) goto 90
       cotan = s(jj)/w(j)
       sin1 = p5/sqrt(p25+p25*cotan**2)
       cos1 = sin1*cotan
       tau = one
       if ((abs(cos1)*giant) > one) tau = one/cos1
       goto 100

90     continue
       tan1 = w(j)/s(jj)
       cos1 = p5/sqrt(p25+p25*tan1**2)
       sin1 = cos1*tan1
       tau = sin1

100    continue

       !  Apply the transformation to s and reduce the spike in w.

       l = jj
       do i = j, m
          temp =  cos1*s(l) + sin1*w(i)
          w(i) = -sin1*s(l) + cos1*w(i)
          s(l) = temp
          l = l + 1
       end do

       !  Store the information necessary to recover the
       !  givens rotation.

       w(j) = tau

120    continue

       !  Test for zero diagonal elements in the output s.

       if (s(jj) == zero) sing = .true.
       jj = jj + (m - j + 1)
    end do

140 continue

    !  Move w back into the last column of the output s.

    l = jj
    do i = n, m
       s(l) = w(i)
       l = l + 1
    end do
    if (s(jj) == zero) sing = .true.

    return
  end SUBROUTINE R1UPDT

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function spmpar(i)

    !! Calculates machine (computing) parameters
    !! author: P J Knight, CCFE, Culham Science Centre
    !! i : input integer : Switch for return value:
    !! i=1 : B**(1 - P), the machine precision
    !! i=2 : B**(EMIN - 1), the smallest magnitude
    !! i=3 : B**EMAX*(1 - B**(-P)), the largest magnitude
    !! where the machine being used has P base B digits, and its smallest
    !! and largest exponents are EMIN and EMAX, respectively.
    !! This routine evaluates the numerical machine parameters of the
    !! computer being used to run the program, as defined above.
    !! <P>Note that the values of these parameters can be found for a given
    !! machine if the Mark 12 or later NAg library is installed on it.
    !! <P><CODE>SPMPAR(1)</CODE> is equivalent to <CODE>X02AJF()</CODE>;
    !! <BR><CODE>SPMPAR(2)</CODE> is equivalent to <CODE>X02AKF()</CODE>;
    !! <BR><CODE>SPMPAR(3)</CODE> is equivalent to <CODE>X02ALF()</CODE>.
    !! Metcalf and Reid, Fortran 90/95 Explained, 2nd Edition (section 8.7.2)
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    implicit none

    real(dp) :: spmpar

    !  Arguments

    integer, intent(in) :: i

    !  Local variables

    real(dp), dimension(3) :: rmach

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Previously-used hardwired values shown in comments

    !rmach(1) = 1.110223024625157D-016
    rmach(1) = epsilon(0.0D0)

    !rmach(2) = 2.3D-308
    rmach(2) = tiny(0.0D0)

    !rmach(3) = 1.797693134862316D+308
    rmach(3) = huge(0.0D0)

    spmpar = rmach(i)

  end function spmpar

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eshellvol(rshell,rmini,rmino,zminor,drin,drout,dz,vin,vout,vtot)

    !! Routine to calculate the inboard, outboard and total volumes
    !! of a toroidal shell comprising two elliptical sections
    !! author: P J Knight, CCFE, Culham Science Centre
    !! rshell : input real : major radius of centre of both ellipses (m)
    !! rmini  : input real : horizontal distance from rshell to outer edge
    !! of inboard elliptical shell (m)
    !! rmino  : input real : horizontal distance from rshell to inner edge
    !! of outboard elliptical shell (m)
    !! zminor : input real : vertical internal half-height of shell (m)
    !! drin   : input real : horiz. thickness of inboard shell at midplane (m)
    !! drout  : input real : horiz. thickness of outboard shell at midplane (m)
    !! dz     : input real : vertical thickness of shell at top/bottom (m)
    !! vin    : output real : volume of inboard section (m3)
    !! vout   : output real : volume of outboard section (m3)
    !! vtot   : output real : total volume of shell (m3)
    !! This routine calculates the volume of the inboard and outboard sections
    !! of a toroidal shell defined by two co-centred semi-ellipses.
    !! Each section's internal and external surfaces are in turn defined
    !! by two semi-ellipses. The volumes of each section are calculated as
    !! the difference in those of the volumes of revolution enclosed by their
    !! inner and outer surfaces.
    !! <P>See also <A HREF="eshellarea.html"><CODE>eshellarea</CODE></A>
    !! Internal CCFE note T&amp;M/PKNIGHT/PROCESS/009, P J Knight:
    !! Surface Area and Volume Calculations for Toroidal Shells
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use constants, only: pi, twopi
    implicit none

    !  Arguments
    real(kind=double), intent(in) :: rshell, rmini, rmino, zminor, drin, drout, dz
    real(kind=double), intent(out) :: vin, vout, vtot

    !  Local variables
    real(kind=double) :: a, b, elong, v1, v2

    !  Global shared variables
    !  Input: pi,twopi
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! #TODO - Review both equations containing dz and attempt to separate
    !         top and bottom of vacuum vessel thickness
    !         See issue #433 for explanation

    !  Inboard section

    !  Volume enclosed by outer (higher R) surface of elliptical section
    !  and the vertical straight line joining its ends
    a = rmini ; b = zminor ; elong = b/a
    v1 = twopi * elong * (0.5D0*pi*rshell*a*a - 2.0D0/3.0D0*a*a*a)

    !  Volume enclosed by inner (lower R) surface of elliptical section
    !  and the vertical straight line joining its ends
    a = rmini+drin ; b = zminor+dz ; elong = b/a
    v2 = twopi * elong * (0.5D0*pi*rshell*a*a - 2.0D0/3.0D0*a*a*a)

    !  Volume of inboard section of shell
    vin = v2 - v1

    !  Outboard section

    !  Volume enclosed by inner (lower R) surface of elliptical section
    !  and the vertical straight line joining its ends
    a = rmino ; b = zminor ; elong = b/a
    v1 = twopi * elong * (0.5D0*pi*rshell*a*a + 2.0D0/3.0D0*a*a*a)

    !  Volume enclosed by outer (higher R) surface of elliptical section
    !  and the vertical straight line joining its ends
    a = rmino+drout ; b = zminor+dz ; elong = b/a
    v2 = twopi * elong * (0.5D0*pi*rshell*a*a + 2.0D0/3.0D0*a*a*a)

    !  Volume of outboard section of shell
    vout = v2 - v1

    !  Total shell volume
    vtot = vin + vout

  end subroutine eshellvol

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dshellvol(rmajor,rminor,zminor,drin,drout,dz,vin,vout,vtot)

    !! Routine to calculate the inboard, outboard and total volumes
    !! of a D-shaped toroidal shell
    !! author: P J Knight, CCFE, Culham Science Centre
    !! rmajor : input real : major radius to outer point of inboard
    !! straight section of shell (m)
    !! rminor : input real : horizontal internal width of shell (m)
    !! zminor : input real : vertical internal half-height of shell (m)
    !! drin   : input real : horiz. thickness of inboard shell at midplane (m)
    !! drout  : input real : horiz. thickness of outboard shell at midplane (m)
    !! dz     : input real : vertical thickness of shell at top/bottom (m)
    !! vin    : output real : volume of inboard straight section (m3)
    !! vout   : output real : volume of outboard curved section (m3)
    !! vtot   : output real : total volume of shell (m3)
    !! This routine calculates the volume of the inboard and outboard sections
    !! of a D-shaped toroidal shell defined by the above input parameters.
    !! The inboard section is assumed to be a cylinder of uniform thickness.
    !! The outboard section's internal and external surfaces are defined
    !! by two semi-ellipses, centred on the outer edge of the inboard section;
    !! its volume is calculated as the difference in those of the volumes of
    !! revolution enclosed by the two surfaces.
    !! <P>See also <A HREF="dshellarea.html"><CODE>dshellarea</CODE></A>
    !! Internal CCFE note T&amp;M/PKNIGHT/PROCESS/009, P J Knight:
    !! Surface Area and Volume Calculations for Toroidal Shells
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use constants, only: pi, twopi
    implicit none

    !  Arguments
    real(kind=double), intent(in) :: rmajor, rminor, zminor, drin, drout, dz
    real(kind=double), intent(out) :: vin, vout, vtot

    !  Local variables
    real(kind=double) :: a, b, elong, v1, v2

    !  Global shared variables
    !  Input: pi,twopi
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! #TODO - Review both equations containing dz and attempt to separate
    !         top and bottom of vacuum vessel thickness
    !         See issue #433 for explanation

    !  Volume of inboard cylindrical shell
    vin = 2.0D0*(zminor+dz) * pi*(rmajor**2 - (rmajor-drin)**2)

    !  Volume enclosed by inner surface of elliptical outboard section
    !  and the vertical straight line joining its ends
    a = rminor ; b = zminor ; elong = b/a
    v1 = twopi * elong * (0.5D0*pi*rmajor*a*a + 2.0D0/3.0D0*a*a*a)

    !  Volume enclosed by outer surface of elliptical outboard section
    !  and the vertical straight line joining its ends
    a = rminor+drout ; b = zminor+dz ; elong = b/a
    v2 = twopi * elong * (0.5D0*pi*rmajor*a*a + 2.0D0/3.0D0*a*a*a)

    !  Volume of elliptical outboard shell
    vout = v2 - v1

    !  Total shell volume
    vtot = vin + vout

  end subroutine dshellvol

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine dshellarea(rmajor,rminor,zminor,ain,aout,atot)

    !! Routine to calculate the inboard, outboard and total surface areas
    !! of a D-shaped toroidal shell
    !! author: P J Knight, CCFE, Culham Science Centre
    !! rmajor : input real : major radius of inboard straight section (m)
    !! rminor : input real : horizontal width of shell (m)
    !! zminor : input real : vertical half-height of shell (m)
    !! ain    : output real : surface area of inboard straight section (m3)
    !! aout   : output real : surface area of outboard curved section (m3)
    !! atot   : output real : total surface area of shell (m3)
    !! This routine calculates the surface area of the inboard and outboard
    !! sections of a D-shaped toroidal shell defined by the above input
    !! parameters.
    !! The inboard section is assumed to be a cylinder.
    !! The outboard section is defined by a semi-ellipse, centred on the
    !! major radius of the inboard section.
    !! <P>See also <A HREF="dshellvol.html"><CODE>dshellvol</CODE></A>
    !! Internal CCFE note T&amp;M/PKNIGHT/PROCESS/009, P J Knight:
    !! Surface Area and Volume Calculations for Toroidal Shells
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use constants, only: pi, twopi
    implicit none

    !  Arguments
    real(dp), intent(in) :: rmajor,rminor,zminor
    real(dp), intent(out) :: ain,aout,atot

    !  Local variables
    real(dp) :: elong

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Area of inboard cylindrical shell
    ain = 4.0D0*zminor*pi*rmajor

    !  Area of elliptical outboard section
    elong = zminor/rminor
    aout = twopi * elong * (pi*rmajor*rminor + 2.0D0*rminor*rminor)

    !  Total surface area
    atot = ain + aout

  end subroutine dshellarea

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine eshellarea(rshell,rmini,rmino,zminor,ain,aout,atot)

    !! Routine to calculate the inboard, outboard and total surface areas
    !! of a toroidal shell comprising two elliptical sections
    !! author: P J Knight, CCFE, Culham Science Centre
    !! rshell : input real : major radius of centre of both ellipses (m)
    !! rmini  : input real : horizontal distance from rshell to
    !! inboard elliptical shell (m)
    !! rmino  : input real : horizontal distance from rshell to
    !! outboard elliptical shell (m)
    !! zminor : input real : vertical internal half-height of shell (m)
    !! ain    : output real : surface area of inboard section (m3)
    !! aout   : output real : surface area of outboard section (m3)
    !! atot   : output real : total surface area of shell (m3)
    !! This routine calculates the surface area of the inboard and outboard
    !! sections of a toroidal shell defined by two co-centred semi-ellipses.
    !! <P>See also <A HREF="eshellvol.html"><CODE>eshellvol</CODE></A>
    !! Internal CCFE note T&amp;M/PKNIGHT/PROCESS/009, P J Knight:
    !! Surface Area and Volume Calculations for Toroidal Shells
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use constants, only: pi, twopi
    implicit none

    !  Arguments
    real(dp), intent(in) :: rshell,rmini,rmino,zminor
    real(dp), intent(out) :: ain,aout,atot

    !  Local variables
    real(dp) :: elong

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Inboard section
    elong = zminor/rmini
    ain = twopi * elong * (pi*rshell*rmini - 2.0D0*rmini*rmini)

    !  Outboard section
    elong = zminor/rmino
    aout = twopi * elong * (pi*rshell*rmino + 2.0D0*rmino*rmino)

    !  Total surface area
    atot = ain + aout

  end subroutine eshellarea

  ! ------------------------------------------------------------------------
  pure function variable_error(variable)
      real(dp), intent(in) ::variable
      logical::variable_error

      if((variable/=variable).or.(variable<-9.99D99).or.(variable>9.99D99))then
          variable_error = .TRUE.
      else
          variable_error = .FALSE.
      end if

  end function variable_error

  ! ------------------------------------------------------------------------
  pure function nearly_equal(variable1, variable2,tol)
      real(dp), intent(in) ::variable1, variable2
      real(dp), intent(in), optional :: tol
      real(dp) :: tolerance
      logical::nearly_equal
      if(present(tol))then
          tolerance = tol
      else
          tolerance = 1.d-5
      end if

      if(abs( (variable1 - variable2)/(variable1+variable2)) < tolerance) then
          nearly_equal = .TRUE.
      else
          nearly_equal = .FALSE.
      end if

  end function nearly_equal

  ! ------------------------------------------------------------------------
  pure function integer2string(value)
      ! Convert an integer value to a 2-digit string with leading zero if required.
      integer, intent(in) :: value
      character(len=2) integer2string
      write (integer2string,'(I2.2)') value
  end function integer2string

  pure function integer3string(value)
      ! Convert an integer value to a 3-digit string with leading zero if required.
      integer, intent(in) :: value
      character(len=3) integer3string
      write (integer3string,'(I3.3)') value
  end function integer3string

  ! ------------------------------------------------------------------------
  subroutine secant_solve(f,x1,x2,solution,error,residual,opt_tol)
      ! Solve an equation f(x)=0
      ! Only requires one function evaluation per iteration (plus two to start with)
      ! Does not require any derivatives.
      ! https://en.wikipedia.org/wiki/Secant_method
      ! Requires two initial values, x0 and x1, which should ideally be chosen to lie close to the root.
      !external:: f

      interface
        function f(x)
#ifndef dp
          use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
          real(dp), intent(in) :: x
          real(dp) :: f
        end function f
      end interface

      external :: f
      real(dp), intent(out) ::solution, residual
      real(dp), intent(in) ::x1,x2
      real(dp), intent(in), optional ::opt_tol
      real(dp),dimension(20) ::x
      integer :: i
      logical, intent(out)::error
      real(dp)::mean, tol,fximinus1, fximinus2
      error = .FALSE.
      tol=0.001d0; if (present(opt_tol)) tol=opt_tol

      x(1)=x1
      x(2)=x2
      mean = (x1+x2)/2
      ! Calculate the first two values before the loop
      fximinus1 = f(x(2))
      fximinus2 = f(x(1))
      !write(*,*)"x(1)",x(1),"x(2)",x(2),"fximinus1",fximinus1,"fximinus2",fximinus2
      do i=3,20
          x(i) = x(i-1) - fximinus1 * ((x(i-1) - x(i-2)) / (fximinus1 - fximinus2))
          ! Store values for the *next iteration
          fximinus2 = fximinus1
          fximinus1 = f(x(i))
          residual = fximinus1
          !write(*,*)"x(i)", x(i), "  residual",residual,"  fximinus1",fximinus1, "  fximinus2",fximinus2
          ! test for convergence
          if(abs(residual) < tol) then
              solution = x(i)
              return
          endif
      end do
      ! Convergence not achieved.  Return the best value and a flag.
      error = .TRUE.
      solution = x(i-1)
      write(*,*)"Secant solver not converged.  solution", solution, "  residual",residual
      !stop 1

  end subroutine secant_solve
!---------------------------------------------------------------

  subroutine test_secant_solve()
      real(dp) ::solution
      real(dp) ::residual
      logical::error
      !external:: f

      call  secant_solve(dummy,10.d0,30.0d0,solution,error,residual)
      if((abs(solution-24.7386)<0.001d0).and.(error.eqv..FALSE.).and.(abs(residual)<0.001d0))then
          write(*,*)"secant solve: PASS.  Error = ", solution-24.7386, ' residual = ', residual
      else
           write(*,*)"secant solve: FAIL. solution=", solution, 'error flag = ', error, "residual = ", residual
           write(*,*)"Correct solution should be 24.7386"
      end if

  contains
      function dummy(x)
          real(dp), intent(in) :: x
          real(dp) :: dummy
          dummy = x**2 - 612.d0
      end function
  end subroutine test_secant_solve

end module maths_library
