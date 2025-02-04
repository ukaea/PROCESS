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

subroutine verror(ifail)

  !! Routine to print out relevant messages in the case of an
  !! unfeasible result from a VMCON (optimisation) run
  !! author: P J Knight, CCFE, Culham Science Centre
  !! ifail  : input integer : error flag
  !! This routine prints out relevant messages in the case of
  !! an unfeasible result from a VMCON (optimisation) run.
  !! <P>The messages are written to units NOUT and IOTTY, which are
  !! by default the output file and screen, respectively.
  !! <P>If <CODE>IFAIL=1</CODE> then a feasible solution has been
  !! found and therefore no error message is required.
  !!   !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use constants, only: nout, iotty
  use process_output, only: ocmmnt, oblnkl
  implicit none

  !  Arguments
  integer, intent(in) :: ifail

  !  Local variables

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  select case (ifail)

  case (:-1)
     call ocmmnt(nout, 'User-terminated execution of VMCON.')
     call ocmmnt(iotty,'User-terminated execution of VMCON.')

  case (0)
     call ocmmnt(nout, 'Improper input parameters to the VMCON routine.')
     call ocmmnt(nout, 'PROCESS coding must be checked.')

     call ocmmnt(iotty,'Improper input parameters to the VMCON routine.')
     call ocmmnt(iotty,'PROCESS coding must be checked.')

  case (1)
     continue

  case (2)
     call ocmmnt(nout,'The maximum number of calls has been reached without solution.')
     call ocmmnt(nout,'The code may be stuck in a minimum in the residual space that is significantly above zero.')
     call oblnkl(nout)
     call ocmmnt(nout,'There is either no solution possible, or the code')
     call ocmmnt(nout,'is failing to escape from a deep local minimum.')
     call ocmmnt(nout,'Try changing the variables in IXC, or modify their initial values.')

     call ocmmnt(iotty,'The maximum number of calls has been reached without solution.')
     call ocmmnt(iotty,'The code may be stuck in a minimum in the residual space that is significantly above zero.')
     call oblnkl(nout)
     call oblnkl(iotty)
     call ocmmnt(iotty,'There is either no solution possible, or the code')
     call ocmmnt(iotty,'is failing to escape from a deep local minimum.')
     call ocmmnt(iotty,'Try changing the variables in IXC, or modify their initial values.')

  case (3)
     call ocmmnt(nout,'The line search required the maximum of 10 calls.')
     call ocmmnt(nout,'A feasible solution may be difficult to achieve.')
     call ocmmnt(nout,'Try changing or adding variables to IXC.')

     call ocmmnt(iotty,'The line search required the maximum of 10 calls.')
     call ocmmnt(iotty,'A feasible solution may be difficult to achieve.')
     call ocmmnt(iotty,'Try changing or adding variables to IXC.')

  case (4)
     call ocmmnt(nout,'An uphill search direction was found.')
     call ocmmnt(nout,'Try changing the equations in ICC, or')
     call ocmmnt(nout,'adding new variables to IXC.')

     call ocmmnt(iotty,'An uphill search direction was found.')
     call ocmmnt(iotty,'Try changing the equations in ICC, or')
     call ocmmnt(iotty,'adding new variables to IXC.')

  case (5)
     call ocmmnt(nout, &
          'The quadratic programming technique was unable to')
     call ocmmnt(nout,'find a feasible point.')
     call oblnkl(nout)
     call ocmmnt(nout,'Try changing or adding variables to IXC, or modify')
     call ocmmnt(nout,'their initial values (especially if only 1 optimisation')
     call ocmmnt(nout,'iteration was performed).')

     call ocmmnt(iotty, &
          'The quadratic programming technique was unable to')
     call ocmmnt(iotty,'find a feasible point.')
     call oblnkl(iotty)
     call ocmmnt(iotty,'Try changing or adding variables to IXC, or modify')
     call ocmmnt(iotty,'their initial values (especially if only 1 optimisation')
     call ocmmnt(iotty,'iteration was performed).')

  case (6)
     call ocmmnt(nout, &
          'The quadratic programming technique was restricted')
     call ocmmnt(nout, &
          'by an artificial bound, or failed due to a singular')
     call ocmmnt(nout,'matrix.')
     call ocmmnt(nout,'Try changing the equations in ICC, or')
     call ocmmnt(nout,'adding new variables to IXC.')

     call ocmmnt(iotty, &
          'The quadratic programming technique was restricted')
     call ocmmnt(iotty, &
          'by an artificial bound, or failed due to a singular')
     call ocmmnt(iotty,'matrix.')
     call ocmmnt(iotty,'Try changing the equations in ICC, or')
     call ocmmnt(iotty,'adding new variables to IXC.')

  case default
     call ocmmnt(nout,'This value of IFAIL should not be possible...')
     call ocmmnt(nout,'See source code for details.')

     call ocmmnt(iotty,'This value of IFAIL should not be possible...')
     call ocmmnt(iotty,'See source code for details.')

  end select

end subroutine verror

end module scan_module
