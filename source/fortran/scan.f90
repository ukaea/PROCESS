! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module scan_module

  !! Module containing routines to perform a parameter scan
  !! author: P J Knight, CCFE, Culham Science Centre
  !! None
  !! This module contains routines to perform a parameter scan
  !! over a range of values of a particular scanning variable.
  !!   !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif
  implicit none

  public

  integer, parameter :: ipnscns = 1000
  !! Maximum number of scan points

  integer, parameter :: ipnscnv = 81
  !! Number of available scan variables

  integer, parameter :: noutvars = 84
  integer, parameter :: width = 110

  integer :: scan_dim
  !! 1-D or 2-D scan switch (1=1D, 2=2D)

  integer :: isweep
  !! Number of scan points to calculate

  integer :: isweep_2
  !! Number of 2D scan points to calculate

  integer :: nsweep
  !! Switch denoting quantity to scan:<UL>
  !!         <LI> 1  aspect
  !!         <LI> 2  hldivlim
  !!         <LI> 3  pnetelin
  !!         <LI> 4  hfact
  !!         <LI> 5  oacdcp
  !!         <LI> 6  walalw
  !!         <LI> 7  beamfus0
  !!         <LI> 8  fqval
  !!         <LI> 9  te
  !!         <LI> 10 boundu(15: fvs)
  !!         <LI> 11 beta_norm_max
  !!         <LI> 12 bootstrap_current_fraction_max
  !!         <LI> 13 boundu(10: hfact)
  !!         <LI> 14 fiooic
  !!         <LI> 15 fjprot
  !!         <LI> 16 rmajor
  !!         <LI> 17 bmxlim
  !!         <LI> 18 gammax
  !!         <LI> 19 boundl(16: ohcth)
  !!         <LI> 20 t_burn_min
  !!         <LI> 21 not used
  !!         <LI> 22 cfactr (N.B. requires iavail=0)
  !!         <LI> 23 boundu(72: fipir)
  !!         <LI> 24 powfmax
  !!         <LI> 25 kappa
  !!         <LI> 26 triang
  !!         <LI> 27 tbrmin (for blktmodel > 0 only)
  !!         <LI> 28 bt
  !!         <LI> 29 coreradius
  !!         <LI> 30 fimpvar # OBSOLETE
  !!         <LI> 31 taulimit
  !!         <LI> 32 epsvmc
  !!         <LI> 33 ttarget
  !!         <LI> 34 qtargettotal
  !!         <LI> 35 lambda_q_omp
  !!         <LI> 36 lambda_target
  !!         <LI> 37 lcon_factor
  !!         <LI> 38 Neon upper limit
  !!         <LI> 39 Argon upper limit
  !!         <LI> 40 Xenon upper limit
  !!         <LI> 41 blnkoth
  !!         <LI> 42 Argon fraction fimp(9)
  !!         <LI> 43 normalised minor radius at which electron cyclotron current drive is maximum
  !!         <LI> 44 Allowable maximum shear stress (Tresca) in tf coil structural material
  !!         <LI> 45 Minimum allowable temperature margin ; tf coils
  !!         <LI> 46 boundu(150) fgwsep
  !!         <LI> 47 impurity_enrichment(9) Argon impurity enrichment
  !!         <LI> 48 TF coil - n_pancake (integer turn winding pack)
  !!         <LI> 49 TF coil - n_layer (integer turn winding pack)
  !!         <LI> 50 Xenon fraction fimp(13)
  !!         <LI> 51 Power fraction to lower DN Divertor ftar
  !!         <LI> 52 SoL radiation fraction
  !!         <LI> 54 GL_nbti upper critical field at 0 Kelvin
  !!         <LI> 55 `shldith` : Inboard neutron shield thickness
  !!         <LI> 56 crypmw_max: Maximum cryogenic power (ixx=164, ixc=87)
  !!         <LI> 57 `bt` lower boundary
  !!         <LI> 58 `scrapli` : Inboard plasma-first wall gap
  !!         <LI> 59 `scraplo` : Outboard plasma-first wall gap
  !!         <LI> 60 sig_tf_wp_max: Allowable stress in TF Coil conduit (Tresca)
  !!         <LI> 61 copperaoh_m2_max : CS coil current / copper area
  !!         <LI> 62 coheof : CS coil current density at EOF
  !!         <LI> 63 ohcth : CS thickness (m)
  !!         <LI> 64 ohhghf : CS height (m)
  !!         <LI> 65 n_cycle_min : Minimum cycles for CS stress model constraint 90
  !!         <LI> 66 oh_steel_frac: Steel fraction in CS coil
  !!         <LI> 67 t_crack_vertical: Initial crack vertical dimension (m) </UL>
  !!         <LI> 68 `inlet_temp_liq' : Inlet temperature of blanket liquid metal coolant/breeder (K)
  !!         <LI> 69 `outlet_temp_liq' : Outlet temperature of blanket liquid metal coolant/breeder (K)
  !!         <LI> 70 `blpressure_liq' : Blanket liquid metal breeder/coolant pressure (Pa)
  !!         <LI> 71 `n_liq_recirc' : Selected number of liquid metal breeder recirculations per day
  !!         <LI> 72 `bz_channel_conduct_liq' : Conductance of liquid metal breeder duct walls (A V-1 m-1)
  !!         <LI> 73 `pnuc_fw_ratio_dcll' : Ratio of FW nuclear power as fraction of total (FW+BB)
  !!         <LI> 74 `f_nuc_pow_bz_struct' : Fraction of BZ power cooled by primary coolant for dual-coolant balnket
  !!         <LI> 75 pitch : pitch of first wall cooling channels (m)
  !!         <LI> 76 etath : Thermal conversion eff.
  !!         <LI> 77 startupratio : Gyrotron redundancy
  !!         <LI> 78 fkind : Multiplier for Nth of a kind costs
  !!         <LI> 79 etaech : ECH wall plug to injector efficiency

  integer :: nsweep_2
  !! nsweep_2 /3/ : switch denoting quantity to scan for 2D scan:

  real(dp), dimension(ipnscns) :: sweep
  !! sweep(ipnscns) /../: actual values to use in scan

  real(dp), dimension(ipnscns) :: sweep_2
  !! sweep_2(ipnscns) /../: actual values to use in 2D scan

  ! Vars in subroutines scan_1d and scan_2d requiring re-initialising before
  ! each new run
  logical :: first_call_1d
  logical :: first_call_2d

contains

  subroutine init_scan_module
    !! Initialise module variables
    implicit none

    scan_dim = 1
    isweep = 0
    isweep_2 = 0
    nsweep = 1
    nsweep_2 = 3
    sweep = 0.0D0
    sweep_2 = 0.0D0
    first_call_1d = .true.
    first_call_2d = .true.
  end subroutine init_scan_module

  subroutine post_optimise(ifail)
  !! Called after calling the optimising equation solver from Python.
  !! author: P J Knight, CCFE, Culham Science Centre
  !! ifail   : input integer : error flag
  !!
  use constraints
  use error_handling
  use numerics
  use process_output
  use utilities, only:upper_case
  use main_module, only:verror
  ! for ipedestal = 2 option
  use global_variables, only: convergence_parameter
  use constants, only: iotty, nout, mfile
  use physics_variables, only: ipedestal
  use define_iteration_variables, only: boundxc, loadxc
  implicit none

  !  Arguments
  integer, intent(in) :: ifail

  !  Local variables
  integer :: ii,inn,iflag
  real(dp) :: summ,xcval,xmaxx,xminn,f,xnorm
  real(dp), dimension(ipeqns) :: con1, con2, err
  character(len=1), dimension(ipeqns) :: sym
  character(len=10), dimension(ipeqns) :: lab
  character(len=60) :: string1, string2

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !  Check on accuracy of solution by summing the
  !  squares of the residuals of the equality constraints
  summ = 0.0D0
  do ii = 1,neqns
     summ = summ + rcm(ii)*rcm(ii)
  end do
  sqsumsq = sqrt(summ)

  !  Turn on error reporting
  errors_on = .true.

  !  Print out information on solution
  call oheadr(nout,'Numerics')
  call ocmmnt(nout,'PROCESS has performed a VMCON (optimisation) run.')
  if (ifail /= 1) then
     !call ocmmnt(nout,'but could not find a feasible set of parameters.')
    !  call oheadr(nout,'PROCESS COULD NOT FIND A FEASIBLE SOLUTION')
    !  call ovarin(iotty,'VMCON error flag (ifail)','',ifail)
     call ovarin(nout,'VMCON error flag','(ifail)',ifail)
     call oheadr(iotty,'PROCESS COULD NOT FIND A FEASIBLE SOLUTION')
     call oblnkl(iotty)

     idiags(1) = ifail ; call report_error(132)

  else
     call ocmmnt(nout,'and found a feasible set of parameters.')
     call oblnkl(nout)
     call ovarin(nout,'VMCON error flag','(ifail)',ifail)
     call oheadr(iotty,'PROCESS found a feasible solution')
  end if

  !call oblnkl(nout)

  !  If necessary, write out a relevant error message
  if (ifail /= 1) then
     call verror(ifail)
     call oblnkl(nout)
     call oblnkl(iotty)
  else
     !  Show a warning if the constraints appear high even if allegedly converged
     if (sqsumsq >= 1.0D-2) then
        call oblnkl(nout)
        call ocmmnt(nout,'WARNING: Constraint residues are HIGH; consider re-running')
        call ocmmnt(nout,'   with lower values of EPSVMC to confirm convergence...')
        call ocmmnt(nout,'   (should be able to get down to about 1.0E-8 okay)')
        call oblnkl(nout)
        call ocmmnt(iotty,'WARNING: Constraint residues are HIGH; consider re-running')
        call ocmmnt(iotty,'   with lower values of EPSVMC to confirm convergence...')
        call ocmmnt(iotty,'   (should be able to get down to about 1.0E-8 okay)')
        call oblnkl(iotty)

        fdiags(1) = sqsumsq ; call report_error(134)

     end if
  end if

  call ovarin(nout,'Number of iteration variables','(nvar)',nvar)
  call ovarin(nout,'Number of constraints (total)','(neqns+nineqns)',neqns+nineqns)
  call ovarin(nout,'Optimisation switch','(ioptimz)',ioptimz)
  call ovarin(nout,'Figure of merit switch','(minmax)',minmax)
!   if (ifail /= 1) then
!      call ovarin(nout,'VMCON error flag','(ifail)',ifail)
!   end if

  objf_name = '"'//trim(lablmm(abs(minmax)))//'"'
  ! Quotes required for string parsing in MFILE
  call ovarst(nout,'Objective function name','(objf_name)',objf_name)
  call ovarre(nout,'Normalised objective function','(norm_objf)',norm_objf, 'OP ')
  call ovarre(nout,'Square root of the sum of squares of the constraint residuals','(sqsumsq)',sqsumsq, 'OP ')
  call ovarre(nout,'VMCON convergence parameter','(convergence_parameter)',convergence_parameter, 'OP ')
  call ovarin(nout,'Number of VMCON iterations','(nviter)',nviter, 'OP ')
  call oblnkl(nout)

  if (ifail == 1) then
     string1 = 'PROCESS has successfully optimised the iteration variables'
  else
     string1 = 'PROCESS has tried to optimise the iteration variables'
  end if

  if (minmax > 0) then
     string2 = ' to minimise the figure of merit: '
  else
     string2 = ' to maximise the figure of merit: '
  end if

  call upper_case(objf_name)
  write(nout,10) trim(string1) // trim(string2),  trim(objf_name)
10 format(a90, t92, a22)

  call oblnkl(nout)

  !  Check which variables are at bounds
  iflag = 0
  do ii = 1,nvar
     xminn = 1.01D0*bondl(ii)
     xmaxx = 0.99D0*bondu(ii)

     if (xcm(ii) < xminn) then
        if (iflag == 0) then
           call ocmmnt(nout, &
                'Certain operating limits have been reached,')
           call ocmmnt(nout, &
                'as shown by the following iteration variables that are')
           call ocmmnt(nout, &
                'at or near to the edge of their prescribed range :')
           call oblnkl(nout)
           iflag = 1
        end if
        xcval = xcm(ii)*scafc(ii)
        !write(nout,30) ii,lablxc(ixc(ii)),xcval,bondl(ii)*scafc(ii)
        write(nout,30) lablxc(ixc(ii)),xcval,bondl(ii)*scafc(ii)
     end if

     if (xcm(ii) > xmaxx) then
        if (iflag == 0) then
           call ocmmnt(nout, &
                'Certain operating limits have been reached,')
           call ocmmnt(nout, &
                'as shown by the following iteration variables that are')
           call ocmmnt(nout, &
                'at or near to the edge of their prescribed range :')
           call oblnkl(nout)
           iflag = 1
        end if
        xcval = xcm(ii)*scafc(ii)
        write(nout,40) lablxc(ixc(ii)),xcval,bondu(ii)*scafc(ii)
     end if
  end do

!30 format(t4,'Variable ',i3,' (',a9, &
!        ',',1pe12.4,') is at or below its lower bound:',1pe12.4)
30 format(t4, a30, '=',1pe12.4,' is at or below its lower bound:',1pe12.4)
40 format(t4, a30, '=',1pe12.4,' is at or above its upper bound:',1pe12.4)
!40 format(t4,'Variable ',i3,' (',a9, &
!        ',',1pe12.4,') is at or above its upper bound:',1pe12.4)

  !  Print out information on numerics
  call osubhd(nout,'The solution vector is comprised as follows :')
!  write(nout,50)
! Remove Lagrange multipliers as no-one understands them.
! MFILE not changed
!50 format(t47,'lower',t59,'upper')

  write(nout,60)
!60 format(t23,'final',t33,'fractional',t46,'Lagrange',t58,'Lagrange')
60 format(t43,'final',t55,'final /')


  write(nout,70)
!70 format(t5,'i',t23,'value',t35,'change',t45,'multiplier', &
!        t57,'multiplier')
70 format(t5,'i',t43,'value',t55,'initial')

  call oblnkl(nout)

  do inn = 1,nvar
     xcs(inn) = xcm(inn)*scafc(inn)
!     write(nout,80) inn,lablxc(ixc(inn)),xcs(inn),xcm(inn), &
!          vlam(neqns+nineqns+inn), vlam(neqns+nineqns+1+inn+nvar)
     write(nout,80) inn,lablxc(ixc(inn)),xcs(inn),xcm(inn)
!80 format(t2,i4,t8,a9,t19,4(1pe12.4))
!80 format(t2,i4,t8,a30,t39,2(1pe12.4))
80 format(t2,i4,t8,a30,t39,1pe12.4, t52, 0pf10.4)
! MDK The 0p is needed because of a bizarre "feature"/bug in fortran:
! the 1p in the previous format continues until changed.
     call ovarre(mfile,lablxc(ixc(inn)),'(itvar'//int_to_string3(inn)//')',xcs(inn))

     !  'Range-normalised' iteration variable values for MFILE:
     !  0.0 (at lower bound) to 1.0 (at upper bound)
     if (bondl(inn) == bondu(inn)) then
        xnorm = 1.0D0
     else
        xnorm = (xcm(inn) - bondl(inn)) / (bondu(inn) - bondl(inn))
        xnorm = max(xnorm, 0.0D0)
        xnorm = min(xnorm, 1.0D0)
     end if
     ! Added ratio final/initial to MFILE
     call ovarre(mfile,trim(lablxc(ixc(inn)))//' (final value/initial value)', &
          '(xcm'//int_to_string3(inn)//')',xcm(inn))
     call ovarre(mfile,trim(lablxc(ixc(inn)))//' (range normalised)', &
          '(nitvar'//int_to_string3(inn)//')',xnorm)
  end do


  call osubhd(nout, &
       'The following equality constraint residues should be close to zero :')

  call constraint_eqns(neqns+nineqns,-1,con1,con2,err,sym,lab)
  write(nout,90)
90 format(t48,'physical',t73,'constraint',t100,'normalised')
  write(nout,100)
100 format(t47,'constraint',t74,'residue',t101,'residue')
  call oblnkl(nout)
  do inn = 1,neqns
     write(nout,110) inn,lablcc(icc(inn)),sym(inn),con2(inn), &
          lab(inn),err(inn),lab(inn),con1(inn)
     call ovarre(mfile,lablcc(icc(inn))//' normalised residue', &
          '(eq_con'//int_to_string3(icc(inn))//')',con1(inn))
  end do
110 format(t2,i4,t8,a33,t46,a1,t47,1pe12.4,t60,a10,t71,1pe12.4,t84,a10,t98,1pe12.4)

  if (nineqns > 0) then
     call osubhd(nout, &
          'The following inequality constraint residues should be greater than or approximately equal to zero :')

     do inn = neqns+1,neqns+nineqns
        !write(nout,120) inn,lablcc(icc(inn)),rcm(inn),vlam(inn)
        write(nout,110) inn,lablcc(icc(inn)),sym(inn),con2(inn), &
                        lab(inn), err(inn), lab(inn)
        call ovarre(mfile,lablcc(icc(inn)),'(ineq_con'//int_to_string3(icc(inn))//')',rcm(inn))
     end do
  end if

! 120 format(t2,i4,t8,a33,t45,1pe12.4,1pe12.4)

end subroutine post_optimise

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module scan_module
