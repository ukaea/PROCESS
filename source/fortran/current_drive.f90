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

  subroutine lhrad(rratio)

      !! Routine to calculate Lower Hybrid wave absorption radius
      !! author: P J Knight, CCFE, Culham Science Centre
      !! rratio  : output real : minor radius of penetration / rminor
      !! This routine determines numerically the minor radius at which the
      !! damping of Lower Hybrid waves occurs, using a Newton-Raphson method.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !! AEA FUS 172: Physics Assessment for the European Reactor Study
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use physics_variables, only: te0
      use error_handling, only: idiags, report_error

      implicit none

      !  Arguments

      real(dp), intent(out) :: rratio

      !  Local variables

      real(dp) :: dgdr,drfind,g0,g1,g2,rat0,rat1,r1,r2
      integer :: lapno
      integer, parameter :: maxlap = 100

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Correction to refractive index (kept within valid bounds)

      drfind = min(0.7D0, max(0.1D0,12.5D0/te0))

      !  Use Newton-Raphson method to establish the correct minor radius
      !  ratio. g is calculated as a function of r / r_minor, where g is
      !  the difference between the results of the two formulae for the
      !  energy E given in AEA FUS 172, p.58. The required minor radius
      !  ratio has been found when g is sufficiently close to zero.

      !  Initial guess for the minor radius ratio

      rat0 = 0.8D0

      lapno = 0
      do ; lapno = lapno+1

         !  Minor radius ratios either side of the latest guess

         r1 = rat0 - 1.0D-3*rat0
         r2 = rat0 + 1.0D-3*rat0

         !  Evaluate g at rat0, r1, r2

         call lheval(drfind,rat0,g0)
         call lheval(drfind,r1,g1)
         call lheval(drfind,r2,g2)

         !  Calculate gradient of g with respect to minor radius ratio

         dgdr = (g2-g1)/(r2-r1)

         !  New approximation

         rat1 = rat0 - g0/dgdr

         !  Force this approximation to lie within bounds

         rat1 = max(0.0001D0,rat1)
         rat1 = min(0.9999D0,rat1)

         !  Check the number of laps for convergence

         if (lapno >= maxlap) then
            idiags(1) = lapno ; call report_error(16)
            rat0 = 0.8D0
            exit
         end if

         !  Is g sufficiently close to zero?

         if (abs(g0) > 0.01D0) then
            !  No, so go around loop again
            rat0 = rat1
         else
            exit
         end if

      end do

      rratio = rat0

    end subroutine lhrad

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lheval(drfind,rratio,ediff)

      !! Routine to evaluate the difference between electron energy
      !! expressions required to find the Lower Hybrid absorption radius
      !! author: P J Knight, CCFE, Culham Science Centre
      !! drfind  : input real : correction to parallel refractive index
      !! rratio  : input real : guess for radius of penetration / rminor
      !! ediff   : output real : difference between the E values (keV)
      !! This routine evaluates the difference between the values calculated
      !! from the two equations for the electron energy E, given in
      !! AEA FUS 172, p.58. This difference is used to locate the Lower Hybrid
      !! wave absorption radius via a Newton-Raphson method, in calling
      !! routine <A HREF="lhrad.html">lhrad</A>.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !! AEA FUS 172: Physics Assessment for the European Reactor Study
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use profiles_module, only: nprofile, tprofile
      use physics_variables, only: rminor, rhopedn, ne0, neped, nesep, alphan, &
         rhopedt, te0, tesep, teped, alphat, tbeta, bt, rmajor

      implicit none

      !  Arguments

      real(dp), intent(in) :: drfind,rratio
      real(dp), intent(out) :: ediff

      !  Local variables

      real(dp) :: blocal,dlocal,e1,e2,frac,nplacc,refind,tlocal

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Local electron density

      dlocal = 1.0D-19 * nprofile(rratio,rhopedn,ne0,neped,nesep,alphan)

      !  Local electron temperature

      tlocal = tprofile(rratio,rhopedt,te0,teped,tesep,alphat,tbeta)

      !  Local toroidal field (evaluated at the inboard region of the flux surface)

      blocal = bt * rmajor/(rmajor - rratio*rminor)

      !  Parallel refractive index needed for plasma access

      frac = sqrt(dlocal)/blocal
      nplacc = frac + sqrt(1.0D0 + frac*frac)

      !  Total parallel refractive index

      refind = nplacc + drfind

      !  First equation for electron energy E

      e1 = 511.0D0 * (sqrt(1.0D0 + 1.0D0/(refind*refind)) - 1.0D0)

      !  Second equation for E

      e2 = 7.0D0 * tlocal

      !  Difference

      ediff = e1-e2

    end subroutine lheval

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

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine culnbi(effnbss,fpion,fshine)

    !! Routine to calculate Neutral Beam current drive parameters
    !! author: P J Knight, CCFE, Culham Science Centre
    !! effnbss : output real : neutral beam current drive efficiency (A/W)
    !! fpion   : output real : fraction of NB power given to ions
    !! fshine  : output real : shine-through fraction of beam
    !! This routine calculates Neutral Beam current drive parameters
    !! using the corrections outlined in AEA FUS 172 to the ITER method.
    !! <P>The result cannot be guaranteed for devices with aspect ratios far
    !! from that of ITER (approx. 2.8).
    !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
    !! AEA FUS 172: Physics Assessment for the European Reactor Study
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use error_handling, only: fdiags, idiags, report_error
    use current_drive_variables, only: enbeam, frbeam, taubeam, ftritbm
    use physics_variables, only: eps, rmajor, abeam, te, dene, ralpne, rncne, &
      rnone, rnfene, dnla, deni, ten, zeffai, dlamie, alphan, alphat, aspect, &
      rminor, zeff

    implicit none

    real(dp), intent(out) :: effnbss,fpion,fshine

    !  Local variables

    real(dp) :: dend,dent,dpath,sigstop

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !  Check argument sanity

    if ((1.0D0+eps) < frbeam) then
       fdiags(1) = eps ; fdiags(2) = frbeam
       call report_error(20)
    end if

    !  Calculate beam path length to centre

    dpath = rmajor * sqrt( (1.0D0 + eps)**2 - frbeam**2)

    !  Calculate beam stopping cross-section

    sigstop = sigbeam(enbeam/abeam,te,dene,ralpne,rncne,rnone,rnfene)

    !  Calculate number of decay lengths to centre

    taubeam = dpath * dnla * sigstop

    !  Shine-through fraction of beam

    fshine = exp(-2.0D0 * dpath*dnla*sigstop)
    fshine = max(fshine, 1.0D-20)

    !  Deuterium and tritium beam densities

    dend = deni * (1.0D0-ftritbm)
    dent = deni * ftritbm

    !  Power split to ions / electrons

    call cfnbi(abeam,enbeam,ten,dene,dend,dent,zeffai,dlamie,fpion)

    !  Current drive efficiency

    effnbss = etanb2(abeam,alphan,alphat,aspect,dene,dnla,enbeam, &
         frbeam,fshine,rmajor,rminor,ten,zeff)

  end subroutine culnbi

  function etanb2(abeam,alphan,alphat,aspect,dene,dnla,enbeam,frbeam, &
      fshine,rmajor,rminor,ten,zeff)
      !! Routine to find neutral beam current drive efficiency
      !! using the ITER 1990 formulation, plus correction terms
      !! outlined in Culham Report AEA FUS 172
      !! author: P J Knight, CCFE, Culham Science Centre
      !! abeam   : input real : beam ion mass (amu)
      !! alphan  : input real : density profile factor
      !! alphat  : input real : temperature profile factor
      !! aspect  : input real : aspect ratio
      !! dene    : input real : volume averaged electron density (m**-3)
      !! dnla    : input real : line averaged electron density (m**-3)
      !! enbeam  : input real : neutral beam energy (keV)
      !! frbeam  : input real : R_tangent / R_major for neutral beam injection
      !! fshine  : input real : shine-through fraction of beam
      !! rmajor  : input real : plasma major radius (m)
      !! rminor  : input real : plasma minor radius (m)
      !! ten     : input real : density weighted average electron temperature (keV)
      !! zeff    : input real : plasma effective charge
      !! This routine calculates the current drive efficiency in A/W of
      !! a neutral beam system, based on the 1990 ITER model,
      !! plus correction terms outlined in Culham Report AEA FUS 172.
      !! <P>The formulae are from AEA FUS 172, unless denoted by IPDG89.
      !! AEA FUS 251: A User's Guide to the PROCESS Systems Code
      !! AEA FUS 172: Physics Assessment for the European Reactor Study
      !! ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
      !! ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990
      !
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use error_handling, only: fdiags, idiags, report_error
      implicit none

      real(dp) :: etanb2

      !  Arguments

      real(dp), intent(in) :: abeam,alphan,alphat,aspect,dene,dnla, &
           enbeam,frbeam,fshine,rmajor,rminor,ten,zeff

      !  Local variables

      real(dp) :: abd,bbd,d,dene20,dnla20,dnorm,ebmev,ebnorm, &
           ecrit,epseff,epsitr,eps1,ffac,gamnb,gfac,j0,nnorm,r,xj, &
           xjs,yj,zbeam

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !  Charge of beam ions
      zbeam = 1.0D0

      !  Fitting factor (IPDG89)
      bbd = 1.0D0

      !  Volume averaged electron density (10**20 m**-3)
      dene20 = dene/1.0D20

      !  Line averaged electron density (10**20 m**-3)
      dnla20 = dnla/1.0D20

      !  Critical energy (MeV) (power to electrons = power to ions) (IPDG89)
      !  N.B. ten is in keV
      ecrit = 0.01D0 * abeam * ten

      !  Beam energy in MeV
      ebmev = enbeam/1.0D3

      !  x and y coefficients of function J0(x,y) (IPDG89)
      xjs = ebmev/(bbd*ecrit)
      xj = sqrt(xjs)

      yj = 0.8D0 * zeff/abeam

      !  Fitting function J0(x,y)
      j0 = xjs / (4.0D0 + 3.0D0*yj + xjs *(xj + 1.39D0 + &
           0.61D0 * yj**0.7D0))

      !  Effective inverse aspect ratio, with a limit on its maximum value
      epseff = min(0.2D0, (0.5D0/aspect))

      !  Reduction in the reverse electron current
      !  due to neoclassical effects
      gfac = (1.55D0 + 0.85D0/zeff)*sqrt(epseff) - &
           (0.2D0 + 1.55D0/zeff)*epseff

      !  Reduction in the net beam driven current
      !  due to the reverse electron current
      ffac = 1.0D0 - (zbeam/zeff) * (1.0D0 - gfac)

      !  Normalisation to allow results to be valid for
      !  non-ITER plasma size and density:

      !  Line averaged electron density (10**20 m**-3) normalised to ITER
      nnorm = 1.0D0

      !  Distance along beam to plasma centre
      r = max(rmajor,rmajor*frbeam)
      eps1 = rminor/r

      if ((1.0D0+eps1) < frbeam) then
         fdiags(1) = eps1 ; fdiags(2) = frbeam
         call report_error(21)
      end if

      d = rmajor * sqrt( (1.0D0+eps1)**2 - frbeam**2)

      !  Distance along beam to plasma centre for ITER
      !  assuming a tangency radius equal to the major radius
      epsitr = 2.15D0/6.0D0
      dnorm = 6.0D0 * sqrt(2.0D0*epsitr + epsitr**2)

      !  Normalisation to beam energy (assumes a simplified formula for
      !  the beam stopping cross-section)
      ebnorm = ebmev * ( (nnorm*dnorm)/(dnla20*d) )**(1.0D0/0.78D0)

      !  A_bd fitting coefficient, after normalisation with ebnorm
      abd = 0.107D0 * (1.0D0 - 0.35D0*alphan + 0.14D0*alphan**2) * &
           (1.0D0 - 0.21D0*alphat) * (1.0D0 - 0.2D0*ebnorm + 0.09D0*ebnorm**2)

      !  Normalised current drive efficiency (A/W m**-2) (IPDG89)
      gamnb = 5.0D0 * abd * 0.1D0*ten * (1.0D0-fshine) * frbeam * &
           j0/0.2D0 * ffac

      !  Current drive efficiency (A/W)
      etanb2 = gamnb / (dene20*rmajor)

    end function etanb2

end module current_drive_module
