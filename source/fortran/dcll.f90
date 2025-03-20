module dcll_module

    !! This module contains the Dual Coolant Lead Lithium (DCLL) specific submods of PROCESSS.
    !!
    !! author: G. Graham, CCFE
    !!
    !! Acronyms for this module:
    !!
    !!      BB          Breeding Blanket
    !!      FW          First Wall
    !!      BZ          Breeder Zone
    !!      MF/BSS      Manifold/Back Supporting Structure
    !!      LT          Low Temperature
    !!      HT          High Temperature
    !!      MMS         Multi Module Segment
    !!      SMS         Single Modle Segment
    !!      IB          Inboard
    !!      OB          Outboard
    !!      HCD         Heating & Current Drive
    !!      FCI         Flow Channel Insert
    !!
    !! IN.DAT info for DCLL:
    !!
    !!      Select DCLL model
    !!          i_blanket_type = 5 * DCLL
    !!
    !!      Liquid Metal Breeder Material = PbLi
    !!          i_blkt_liquid_breeder_type = 0 * Liquid Metal Breeder Material = PbLi
    !!
    !!      Specify dual-coolant i.e., get mass flow required from heat extracted from liqid metal breeder
    !!          i_blkt_dual_coolant = 2
    !!
    !!      FIC switch: 0 = no FIC, Eurofer; 1 = FCIs, perfect electrical insulator, 2 = FCIs, with specified conductance
    !!          i_blkt_liquid_breeder_channel_type = 0, 1, or 2
    !!
    !!      Liquid metal duct wall conductance initilized at Eurofer value in fwbs_variables, or can input other value, used for i_blkt_liquid_breeder_channel_type = 0 or 2
    !!          (bz_channel_conduct_liq)
    !!
    !!      Choose if FW and BB structure are on the same pumping system (unless have diffent coolants), default is same coolant with flow IN->FW->BB->OUT
    !!          (i_fw_blkt_shared_coolant)
    !!
    !!      Can set inlet and oulet temperature for liquid metal breeder
    !!          (inlet_temp_liq)
    !!          (outlet_temp_liq)
    !!
    !! References:
    !!
    !!      [Nat1995]   Natesan et al. (1995), Assessment of alkali metal coolants for
    !!                  for the ITER blanket, Fusion Engineering and Design 27, 457-466
    !!
    !!      [Mal1995]   Malang and Mattas (1995), Comparison of lithium and the eutectic
    !!                  lead-lithium alloy, two candidate liquid metal breeder materials
    !!                  for self-cooled blankets, Fusion Engineering and Design 27, 399-406
    !!
    !!      [Gas2001]   Gasior and Mozer (2001), Thermodynamic study of liquid lithium–lead
    !!                  alloys using the EMF method, Journal of Nuclear Materials, 294, 77-83
    !!
    !!      [Pal2016]   Palermo et al. (2016), Neutronic analyses of the preliminary design
    !!                  of a DCLL blanket for the EUROfusion DEMO power plant,
    !!                  Fusion Engineering and Design 109-111.
    !!
    !!      [Gar2017]   Garcinuno et al. (2017), Design of a permeator against vacuum for
    !!                  tritium extraction from eutectic lithium-lead in a DCLL DEMO,
    !!                  Fusion Engineering and Design, 117, 226-231
    !!
    !!      [Fer2021]   Fernandez-Berceruelo et al. (2021), Alternatives for upgrading the
    !!                  eu dcll breeding blanket from mms to sms, Fusion Engineering and
    !!                  Design 167, 112380
    !!
    !!
    !! Note: request for when CCFE Bluemira nutronics work is added: output maximum values, as well as average values, for wall neutronics calculation if possible.
    !!
    !! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifndef dp
    use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

    use fwbs_variables, only: wht_liq, wht_liq_ib, wht_liq_ob

    implicit none

    !! DCLL Variables !!!!!!!!!!!!!!!!

    !! Radial thickness of FCIs annd backwall [m]
    real(dp) :: r_fci, r_backwall

    !! Radial BZ thickness [m]
    real(dp) :: bz_r_ib, bz_r_ob

    !! Structure/coolant compositional fractions
    real(dp) :: f_vol_stff_plates, f_vol_stl_bz_struct, f_vol_stl_back_wall, f_vol_stl_fw

    !! MF/BSS compositional fractions
    real(dp) :: f_vol_mfbss_stl, f_vol_mfbss_he, f_vol_mfbss_pbli

    !! Volume of FCIs, other BZ structure, liquid channels, backwall and MF/BSS [m^3]
    real(dp) :: vol_fci, vol_bz_struct, vol_bz_liq, vol_bz_liq_ib, vol_bz_liq_ob, vol_bw, vol_bss

    !! BZ masses by composition [kg]
    real(dp) :: wht_cer, wht_stl_struct, wht_cool_struct

    !! Backwall masses by composition [kg]
    real(dp) :: wht_bw_stl, wht_bw_cool

    !! MF/BSS masses by composition [kg]
    real(dp) :: wht_mfbss_stl, wht_mfbss_cool, wht_mfbss_pbli

    !! FW masses by composition [kg]
    real(dp) :: fwmass_stl, fwmass_cool

    !! Total masses of material in blanket [kg]
    real(dp) :: mass_cool_blanket, mass_liq_blanket, mass_stl_blanket

    !! Total mass for an inboard/outboard reactor segment [kg]
    real(dp) :: mass_segm_ib, mass_segm_ob

end module dcll_module
