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

contains
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
    real(kind=dp), intent(in) :: rshell, rmini, rmino, zminor, drin, drout, dz
    real(kind=dp), intent(out) :: vin, vout, vtot

    !  Local variables
    real(kind=dp) :: a, b, elong, v1, v2

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
    real(kind=dp), intent(in) :: rmajor, rminor, zminor, drin, drout, dz
    real(kind=dp), intent(out) :: vin, vout, vtot

    !  Local variables
    real(kind=dp) :: a, b, elong, v1, v2

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

end module maths_library
