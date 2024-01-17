module heat_transport_variables
    !! author: J. Morris, M. Kovari (UKAEA)
    !!
    !! This module contains global variables relating to the heat transport system
    !! of a fusion power plant, and also those for a hydrogen production plant.
    !!
    !!### References
    !!
    !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    public

    real(dp) :: baseel
    !! base plant electric load (W)

    real(dp) :: crypmw
    !! cryogenic plant power (MW)

    real(dp) :: crypmw_max
    !! Maximum cryogenic plant power (MW)
    !! Constraint equation icc = 87
    !! Scan variable nwseep = 56

    real(dp) :: f_crypmw
    !! f-value for maximum cryogenic plant power
    !! Iteration variable ixc = 164
    !! Constraint equation icc = 87

    real(dp) :: etatf
    !! AC to resistive power conversion for TF coils

    real(dp) :: etath
    !! thermal to electric conversion efficiency if `secondary_cycle=2`; otherwise calculated.

    real(dp) :: etath_liq

    real(dp) :: fachtmw
    !! facility heat removal (MW)

    real(dp) :: fcsht
    !! total baseline power required at all times (MW)

    real(dp) :: fgrosbop
    !! scaled fraction of gross power to balance-of-plant

    real(dp) :: fmgdmw
    !! power to mgf (motor-generator flywheel) units (MW) (ignored if `iscenr=2`)

    real(dp) :: fpumpblkt
    !! fraction of total blanket thermal power required to drive the blanket
    !! coolant pumps (default assumes water coolant) (`secondary_cycle=0`)

    real(dp) :: fpumpdiv
    !! fraction of total divertor thermal power required to drive the divertor
    !! coolant pumps (default assumes water coolant)

    real(dp) :: fpumpfw
    !! fraction of total first wall thermal power required to drive the FW coolant
    !! pumps (default assumes water coolant) (`secondary_cycle=0`)

    real(dp) :: fpumpshld
    !! fraction of total shield thermal power required to drive the shield coolant
    !! pumps (default assumes water coolant)

    real(dp) :: htpmw_min
    !! Minimum total electrical power for primary coolant pumps (MW) (NOT RECOMMENDED)

    real(dp) :: helpow
    !! Heat removal at cryogenic temperature tmpcry (W)

    real(dp) :: helpow_cryal
    !! Heat removal at cryogenic temperature tcoolin (W)

    real(dp) :: htpmw
    !! heat transport system electrical pump power (MW)

    real(dp) :: htpmw_blkt
    !! blanket primary coolant mechanical pumping power (MW)

    real(dp) :: htpmw_blkt_liq
    !! blanket secondary coolant mechanical pumping power (MW)

    real(dp) :: htpmw_blkt_tot
    !! blanket primary + secondary coolant mechanical pumping power (MW)

    real(dp) :: htpmw_div
    !! divertor coolant mechanical pumping power (MW)

    real(dp) :: htpmw_fw
    !! first wall coolant mechanical pumping power (MW)

    real(dp) :: htpmw_shld
    !! shield and vacuum vessel coolant mechanical pumping power (MW)

    real(dp) :: htpsecmw
    !! Waste power lost from primary coolant pumps (MW)

    integer :: ipowerflow
    !! switch for power flow model:
    !!
    !! - =0 pre-2014 version
    !! - =1 comprehensive 2014 model

    integer :: iprimshld
    !! Switch for shield thermal power destiny:
    !!
    !! - =0 does not contribute to energy generation cycle
    !! - =1 contributes to energy generation cycle

    integer :: nphx
    !! number of primary heat exchangers

    real(dp) :: pacpmw
    !! total pulsed power system load (MW)

    real(dp) :: peakmva
    !! peak MVA requirement

    real(dp) :: pfwdiv
    !! heat removal from first wall/divertor (MW)

    real(dp) :: pgrossmw
    !! gross electric power (MW)

    real(dp) :: pinjht
    !! power dissipated in heating and current drive system (MW)

    real(dp) :: pinjmax
    !! maximum injector power during pulse (heating and ramp-up/down phase) (MW)

    real(dp) :: pinjwp
    !! injector wall plug power (MW)

    real(dp) :: pinjwpfix
    !! secondary injector wall plug power (MW)

    real(dp) :: pnetelmw
    !! net electric power (MW)

    real(dp) :: precircmw
    !! recirculating electric power (MW)

    real(dp) :: priheat
    !! total thermal power removed from fusion core (MW)

    real(dp) :: psecdiv
    !! Low-grade heat lost in divertor (MW)

    real(dp) :: psechcd
    !! Low-grade heat lost into HCD apparatus (MW)

    real(dp) :: psechtmw
    !! Low-grade heat (MW)

    real(dp) :: pseclossmw
    !! Low-grade heat (VV + lost)(MW)

    real(dp) :: psecshld
    !! Low-grade heat deposited in shield (MW)

    real(dp) :: pthermmw
    !! High-grade heat useful for electric production (MW)

    real(dp) :: pwpm2
    !! base AC power requirement per unit floor area (W/m2)

    real(dp) :: tfacpd
    !! total steady state TF coil AC power demand (MW)

    real(dp) :: tlvpmw
    !! estimate of total low voltage power (MW)

    real(dp) :: trithtmw
    !! power required for tritium processing (MW)

    real(dp) :: tturb
    !! coolant temperature at turbine inlet (K) (`secondary_cycle = 3,4`)

    real(dp) :: vachtmw
    !! vacuum pump power (MW)

    contains

    subroutine init_heat_transport_variables
      !! Initialise module variables
      implicit none

      baseel = 5.0D6
      crypmw = 0.0D0
      crypmw_max = 50.0D0
      f_crypmw = 1.0D0
      etatf = 0.9D0
      etath = 0.35D0
      etath_liq = 0.35D0
      fachtmw = 0.0D0
      fcsht = 0.0D0
      fgrosbop = 0.0D0
      fmgdmw = 0.0D0
      fpumpblkt = 0.005D0
      fpumpdiv = 0.005D0
      fpumpfw = 0.005D0
      fpumpshld = 0.005D0
      htpmw_min = 0.0D0
      helpow = 0.0D0
      helpow_cryal = 0.0D0
      htpmw = 0.0D0
      htpmw_blkt = 0.0D0
      htpmw_blkt_liq = 0.0D0
      htpmw_blkt_tot = 0.0D0
      htpmw_div = 0.0D0
      htpmw_fw = 0.0D0
      htpmw_shld = 0.0D0
      htpsecmw = 0.0D0
      ipowerflow = 1
      iprimshld = 1
      nphx = 0
      pacpmw = 0.0D0
      peakmva = 0.0D0
      pfwdiv = 0.0D0
      pgrossmw = 0.0D0
      pinjht = 0.0D0
      pinjmax = 120.0D0
      pinjwp = 0.0D0
      pinjwpfix = 0.0D0
      pnetelmw = 0.0D0
      precircmw = 0.0D0
      priheat = 0.0D0
      psecdiv = 0.0D0
      psechcd = 0.0D0
      psechtmw = 0.0D0
      pseclossmw = 0.0D0
      psecshld = 0.0D0
      pthermmw = 0.0D0
      pwpm2 = 150.0D0
      tfacpd = 0.0D0
      tlvpmw = 0.0D0
      trithtmw = 15.0D0
      tturb = 0.0D0
      vachtmw = 0.5D0
    end subroutine init_heat_transport_variables
  end module heat_transport_variables
