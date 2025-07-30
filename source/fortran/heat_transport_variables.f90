module heat_transport_variables
    !! author: J. Morris, M. Kovari (UKAEA)
    !!
    !! This module contains global variables relating to the heat transport system
    !! of a fusion power plant, and also those for a hydrogen production plant.
    !!
    !!### References
    !!
    !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    implicit none

    public

    real(dp) :: p_plant_electric_base
    !! base plant electric load (W)

    real(dp) :: p_cryo_plant_electric_mw
    !! cryogenic plant power (MW)

    real(dp) :: p_cryo_plant_electric_max_mw
    !! Maximum cryogenic plant power (MW)
    !! Constraint equation icc = 87
    !! Scan variable nwseep = 56

    real(dp) :: f_crypmw
    !! f-value for maximum cryogenic plant power
    !! Iteration variable ixc = 164
    !! Constraint equation icc = 87

    real(dp) :: etatf
    !! AC to resistive power conversion for TF coils

    real(dp) :: eta_turbine
    !! thermal to electric conversion efficiency if `i_thermal_electric_conversion=2`; otherwise calculated.

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
    !! coolant pumps (default assumes water coolant) (`i_thermal_electric_conversion=0`)

    real(dp) :: fpumpdiv
    !! fraction of total divertor thermal power required to drive the divertor
    !! coolant pumps (default assumes water coolant)

    real(dp) :: fpumpfw
    !! fraction of total first wall thermal power required to drive the FW coolant
    !! pumps (default assumes water coolant) (`i_thermal_electric_conversion=0`)

    real(dp) :: fpumpshld
    !! fraction of total shield thermal power required to drive the shield coolant
    !! pumps (default assumes water coolant)

    real(dp) :: helpow
    !! Heat removal at cryogenic temperature temp_tf_cryo (W)

    real(dp) :: helpow_cryal
    !! Heat removal at cryogenic temperature tcoolin (W)

    real(dp) :: p_coolant_pump_elec_total_mw
    !! heat transport system electrical pump power (MW)

    real(dp) :: p_blkt_coolant_pump_mw
    !! blanket primary coolant mechanical pumping power (MW)

    real(dp) :: p_blkt_breeder_pump_mw
    !! blanket secondary coolant mechanical pumping power (MW)

    real(dp) :: htpmw_blkt_tot
    !! blanket primary + secondary coolant mechanical pumping power (MW)

    real(dp) :: p_div_coolant_pump_mw
    !! divertor coolant mechanical pumping power (MW)

    real(dp) :: p_fw_coolant_pump_mw
    !! first wall coolant mechanical pumping power (MW)

    real(dp) :: p_shld_coolant_pump_mw
    !! shield and vacuum vessel coolant mechanical pumping power (MW)

    real(dp) :: p_coolant_pump_loss_total_mw
    !! Waste power lost from primary coolant pumps (MW)

    integer :: ipowerflow
    !! switch for power flow model:
    !!
    !! - =0 pre-2014 version
    !! - =1 comprehensive 2014 model

    integer :: i_shld_primary_heat
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

    real(dp) :: p_fw_div_heat_deposited_mw
    !! heat removal from first wall/divertor (MW)

    real(dp) :: p_plant_electric_gross_mw
    !! gross electric power (MW)

    real(dp) :: pinjht
    !! power dissipated in heating and current drive system (MW)

    real(dp) :: pinjmax
    !! maximum injector power during pulse (heating and ramp-up/down phase) (MW)

    real(dp) :: p_hcd_electric_total_mw
    !! injector wall plug power (MW)

    real(dp) :: p_hcd_secondary_electric_mw
    !! Secondary HCD system injector wall plug power (MW)

    real(dp) :: p_hcd_primary_electric_mw
    !! Primary HCD system injector wall plug power (MW)

    real(dp) :: p_plant_electric_net_mw
    !! net electric power (MW)

    real(dp) :: p_plant_electric_recirc_mw
    !! recirculating electric power (MW)

    real(dp) :: priheat
    !! total thermal power removed from fusion core (MW)

    real(dp) :: p_div_secondary_heat_mw
    !! Low-grade heat lost in divertor (MW)

    real(dp) :: p_hcd_secondary_heat_mw
    !! Low-grade heat lost into HCD apparatus (MW)

    real(dp) :: p_plant_secondary_heat_mw
    !! Low-grade heat (MW)

    real(dp) :: pseclossmw
    !! Low-grade heat (VV + lost)(MW)

    real(dp) :: p_shld_secondary_heat_mw
    !! Low-grade heat deposited in shield (MW)

    real(dp) :: p_plant_primary_heat_mw
    !! High-grade heat useful for electric production (MW)

    real(dp) :: pwpm2
    !! base AC power requirement per unit floor area (W/m2)

    real(dp) :: tfacpd
    !! total steady state TF coil AC power demand (MW)

    real(dp) :: tlvpmw
    !! estimate of total low voltage power (MW)

    real(dp) :: p_tritium_plant_electric_mw
    !! power required for tritium processing (MW)

    real(dp) :: temp_turbine_coolant_in
    !! coolant temperature at turbine inlet (K) (`i_thermal_electric_conversion = 3,4`)

    real(dp) :: vachtmw
    !! vacuum pump power (MW)
  end module heat_transport_variables
