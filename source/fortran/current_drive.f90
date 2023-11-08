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
  function etanb(abeam,alphan,alphat,aspect,dene,ebeam,rmajor,ten,zeff)

      !! Routine to find neutral beam current drive efficiency
      !! using the ITER 1990 formulation
      !! author: P J Knight, CCFE, Culham Science Centre
      !! abeam   : input real : beam ion mass (amu)
      !! alphan  : input real : density profile factor
      !! alphat  : input real : temperature profile factor
      !! aspect  : input real : aspect ratio
      !! dene    : input real : volume averaged electron density (m**-3)
      !! enbeam  : input real : neutral beam energy (keV)
      !! rmajor  : input real : plasma major radius (m)
      !! ten     : input real : density weighted average electron temp. (keV)
      !! zeff    : input real : plasma effective charge
      !! This routine calculates the current drive efficiency of
      !! a neutral beam system, based on the 1990 ITER model.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
      !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: etanb

      ! Arguments !
      ! !!!!!!!!!!!!

      real(dp), intent(in) :: abeam,alphan,alphat,aspect,dene, &
           ebeam,rmajor,ten,zeff

      ! Local variables !
      ! !!!!!!!!!!!!!!!!!!

      real(dp) :: abd,bbd,dene20,dum,epseff,ffac,gfac,rjfunc, &
           xj,xjs,yj,zbeam

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! TODO comment this subroutine.

      zbeam = 1.0D0
      bbd = 1.0D0

      dene20 = 1.0D-20*dene

      ! Ratio of E_beam/E_crit
      xjs = ebeam / (bbd*10.0D0*abeam*ten)
      xj = sqrt(xjs)

      yj = 0.8D0 * zeff/abeam

      rjfunc = xjs / (4.0D0 + 3.0D0*yj + xjs * &
           (xj + 1.39D0 + 0.61D0*yj**0.7D0))

      epseff = 0.5D0/aspect
      gfac = (1.55D0 + 0.85D0/zeff)*sqrt(epseff) - &
           (0.2D0 + 1.55D0/zeff)*epseff
      ffac = 1.0D0/zbeam - (1.0D0 - gfac)/zeff

      abd = 0.107D0 * (1.0D0 - 0.35D0*alphan + 0.14D0*alphan**2) * &
           (1.0D0 - 0.21D0*alphat) * (1.0D0 - 0.2D-3*ebeam  + &
           0.09D-6 * ebeam**2)
      dum = abd *(5.0D0/rmajor) * (0.1D0*ten/dene20) * rjfunc/0.2D0 * ffac

      etanb = dum

    end function etanb

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine cfnbi(afast,efast,te,ne,nd,nt,zeffai,xlmbda,fpion)

    !! Routine to calculate the fraction of the fast particle energy
    !! coupled to the ions
    !! author: P J Knight, CCFE, Culham Science Centre
    !! afast   : input real : mass of fast particle (units of proton mass)
    !! efast   : input real : energy of fast particle (keV)
    !! te      : input real : density weighted average electron temp. (keV)
    !! ne      : input real : volume averaged electron density (m**-3)
    !! nd      : input real : deuterium beam density (m**-3)
    !! nt      : input real : tritium beam density (m**-3)
    !! zeffai  : input real : mass weighted plasma effective charge
    !! xlmbda  : input real : ion-electron coulomb logarithm
    !! fpion   : output real : fraction of fast particle energy coupled to ions
    !! This routine calculates the fast particle energy coupled to
    !! the ions in the neutral beam system.
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use constants, only: mproton, pi, echarge

    implicit none

    !  Arguments

    real(dp), intent(in) :: afast,efast,te,ne,nd,nt,zeffai,xlmbda
    real(dp), intent(out) :: fpion

    !  Local variables

    real(dp) :: ans,ecritfi,ecritfix,sum,sumln,thx,t1,t2,ve,x, &
         xlbd,xlbt,xlmbdai,xlnrat

    real(dp), parameter :: atmd = 2.0D0
    real(dp), parameter :: atmdt = 2.5D0
    real(dp), parameter :: atmt = 3.0D0
    real(dp), parameter :: c = 3.0D8
    real(dp), parameter :: me = 9.1D-31
    real(dp), parameter :: zd = 1.0D0
    real(dp), parameter :: zt = 1.0D0

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    xlbd = xlmbdabi(afast,atmd,efast,te,ne)
    xlbt = xlmbdabi(afast,atmt,efast,te,ne)

    sum = nd*zd*zd * xlbd/atmd + nt*zt*zt * xlbt/atmt
    ecritfix = 16.0D0 * te * afast * (sum/(ne*xlmbda))**(2.0D0/3.0D0)

    xlmbdai = xlmbdabi(afast,atmdt,efast,te,ne)
    sumln = zeffai * xlmbdai/xlmbda
    xlnrat = (3.0D0*sqrt(pi)/4.0D0 * me/mproton * sumln)**(2.0D0/3.0D0)
    ve = c * sqrt(2.0D0*te/511.0D0)

    ecritfi = afast * mproton * ve*ve * xlnrat/(2.0D0 * echarge * 1.0D3)

    x = sqrt(efast/ecritfi)
    t1 = log( (x*x - x + 1.0D0) / ((x + 1.0D0)**2) )
    thx = (2.0D0*x - 1.0D0)/sqrt(3.0D0)
    t2 = 2.0D0*sqrt(3.0D0) *(atan(thx) + pi/6.0D0)

    ans = (t1 + t2) / (3.0D0 * x*x)
    fpion = ans

  end subroutine cfnbi

  function xlmbdabi(mb,mth,eb,t,nelec)

      !! Calculates the Coulomb logarithm for ion-ion collisions
      !! author: P J Knight, CCFE, Culham Science Centre
      !! mb     : input real : mass of fast particle (units of proton mass)
      !! mth    : input real : mass of background ions (units of proton mass)
      !! eb     : input real : energy of fast particle (keV)
      !! t      : input real : density weighted average electron temp. (keV)
      !! nelec  : input real : volume averaged electron density (m**-3)
      !! This function calculates the Coulomb logarithm for ion-ion
      !! collisions where the relative velocity may be large compared
      !! with the background ('mt') thermal velocity.
      !! Mikkelson and Singer, Nuc Tech/Fus, 4, 237 (1983)
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      real(dp) :: xlmbdabi

      !  Arguments

      real(dp), intent(in) :: mb,mth,eb,t,nelec

      !  Local variables

      real(dp) :: ans,x1,x2

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      x1 = (t/10.0D0) * (eb/1000.0D0) * mb/(nelec/1.0D20)
      x2 = mth/(mth + mb)

      ans = 23.7D0 + log(x2 * sqrt(x1))
      xlmbdabi = ans

    end function xlmbdabi

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine legend(zlocal,arg,palpha,palphap)

      !! Routine to calculate Legendre function and its derivative
      !! author: M R O'Brien, CCFE, Culham Science Centre
      !! author: P J Knight, CCFE, Culham Science Centre
      !! zlocal  : input real : local plasma effective charge
      !! arg     : input real : argument of Legendre function
      !! palpha  : output real : value of Legendre function
      !! palphap : output real : derivative of Legendre function
      !! This routine calculates the Legendre function <CODE>palpha</CODE>
      !! of argument <CODE>arg</CODE> and order
      !! <CODE>alpha = -0.5 + i sqrt(xisq)</CODE>,
      !! and its derivative <CODE>palphap</CODE>.
      !! <P>This Legendre function is a conical function and we use the series
      !! in <CODE>xisq</CODE> given in Abramowitz and Stegun. The
      !! derivative is calculated from the derivative of this series.
      !! <P>The derivatives were checked by calculating <CODE>palpha</CODE> for
      !! neighbouring arguments. The calculation of <CODE>palpha</CODE> for zero
      !! argument was checked by comparison with the expression
      !! <CODE>palpha(0) = 1/sqrt(pi) * cos(pi*alpha/2) * gam1 / gam2</CODE>
      !! (Abramowitz and Stegun, eqn 8.6.1). Here <CODE>gam1</CODE> and
      !! <CODE>gam2</CODE> are the Gamma functions of arguments
      !! <CODE>0.5*(1+alpha)</CODE> and <CODE>0.5*(2+alpha)</CODE> respectively.
      !! Abramowitz and Stegun, equation 8.12.1
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use error_handling, only: fdiags, idiags, report_error

      implicit none

      real(dp), intent(in) :: zlocal,arg
      real(dp), intent(out) ::  palpha,palphap

      !  Local variables

      real(dp) :: arg2,pold,poldp,pterm,sinsq,term1,term2,xisq
      integer :: n

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Check for invalid argument

      if (abs(arg) > (1.0D0+1.0D-10)) then
         fdiags(1) = arg ; call report_error(18)
      end if

      arg2 = min(arg, (1.0D0-1.0D-10))
      sinsq = 0.5D0*(1.0D0-arg2)
      xisq = 0.25D0*(32.0D0*zlocal/(zlocal+1.0D0) - 1.0D0)
      palpha = 1.0D0
      pold = 1.0D0
      pterm = 1.0D0
      palphap = 0.0D0
      poldp = 0.0D0

      do n = 1,10001

         !  Check for convergence every 20 iterations

         if ((n > 1).and.(mod(n,20) == 1)) then
            term1 = 1.0D-10 * max(abs(pold),abs(palpha))
            term2 = 1.0D-10 * max(abs(poldp),abs(palphap))

            if ( (abs(pold-palpha) < term1) .and. &
                 (abs(poldp-palphap) < term2) ) return

            pold = palpha
            poldp = palphap
         end if

         pterm = pterm * (4.0D0*xisq+(2.0D0*n - 1.0D0)**2) / &
              (2.0D0*n)**2 * sinsq
         palpha = palpha + pterm
         palphap = palphap - n*pterm/(1.0D0-arg2)

      end do

      !  Code will only get this far if convergence has failed

      call report_error(19)

    end subroutine legend

end module current_drive_module
