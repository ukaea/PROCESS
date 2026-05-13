"""This module contains the Dual Coolant Lead Lithium (DCLL) specific submods of PROCESSS.



Acronyms for this module:

     BB          Breeding Blanket
     FW          First Wall
     BZ          Breeder Zone
     MF/BSS      Manifold/Back Supporting Structure
     LT          Low Temperature
     HT          High Temperature
     MMS         Multi Module Segment
     SMS         Single Module Segment
     IB          Inboard
     OB          Outboard
     HCD         Heating & Current Drive
     FCI         Flow Channel Insert

IN.DAT info for DCLL:

     Select DCLL model
         i_blanket_type = 5 * DCLL

     Liquid Metal Breeder Material = PbLi
         i_blkt_liquid_breeder_type = 0 * Liquid Metal Breeder Material = PbLi

     Specify dual-coolant i.e., get mass flow required from heat extracted from liquid metal breeder
         i_blkt_dual_coolant = 2

     FIC switch: 0 = no FIC, Eurofer; 1 = FCIs, perfect electrical insulator, 2 = FCIs, with specified conductance
         i_blkt_liquid_breeder_channel_type = 0, 1, or 2

     Liquid metal duct wall conductance initialized at Eurofer value in fwbs_variables, or can input other value, used for i_blkt_liquid_breeder_channel_type = 0 or 2
         (bz_channel_conduct_liq)

     Choose if FW and BB structure are on the same pumping system (unless have different coolants), default is same coolant with flow IN->FW->BB->OUT
         (i_fw_blkt_shared_coolant)

     Can set inlet and outlet temperature for liquid metal breeder
         (inlet_temp_liq)
         (outlet_temp_liq)

References
----------
     [Nat1995]   Natesan et al. (1995), Assessment of alkali metal coolants for
                 for the ITER blanket, Fusion Engineering and Design 27, 457-466

     [Mal1995]   Malang and Mattas (1995), Comparison of lithium and the eutectic
                 lead-lithium alloy, two candidate liquid metal breeder materials
                 for self-cooled blankets, Fusion Engineering and Design 27, 399-406

     [Gas2001]   Gasior and Mozer (2001), Thermodynamic study of liquid lithium-lead
                 alloys using the EMF method, Journal of Nuclear Materials, 294, 77-83

     [Pal2016]   Palermo et al. (2016), Neutronic analyses of the preliminary design
                 of a DCLL blanket for the EUROfusion DEMO power plant,
                 Fusion Engineering and Design 109-111.

     [Gar2017]   Garcinuno et al. (2017), Design of a permeator against vacuum for
                 tritium extraction from eutectic lithium-lead in a DCLL DEMO,
                 Fusion Engineering and Design, 117, 226-231

     [Fer2021]   Fernandez-Berceruelo et al. (2021), Alternatives for upgrading the
                 eu dcll breeding blanket from mms to sms, Fusion Engineering and
                 Design 167, 112380


Note: request for when CCFE Bluemira neutronics work is added: output maximum values, as well as average values, for wall neutronics calculation if possible.
"""

from dataclasses import dataclass


@dataclass
class DCLLData:
    r_fci: float = 0.0
    """Radial thickness of FCIs [m]"""

    r_backwall: float = 0.0
    """Radial thickness of backwall [m]"""

    bz_r_ib: float = 0.0
    """Radial BZ thickness [m]"""

    bz_r_ob: float = 0.0
    """Radial BZ thickness [m]"""

    f_vol_stff_plates: float = 0.0
    """Structure/coolant compositional fractions"""

    f_vol_stl_bz_struct: float = 0.0
    """Structure/coolant compositional fractions"""

    f_vol_stl_back_wall: float = 0.0
    """Structure/coolant compositional fractions"""

    f_vol_stl_fw: float = 0.0
    """Structure/coolant compositional fractions"""

    f_vol_mfbss_stl: float = 0.0
    """MF/BSS compositional fractions"""

    f_vol_mfbss_he: float = 0.0
    """MF/BSS compositional fractions"""

    f_vol_mfbss_pbli: float = 0.0
    """MF/BSS compositional fractions"""

    vol_fci: float = 0.0
    """Volume of FCIs [m^3]"""

    vol_bz_struct: float = 0.0
    """Volume of other BZ structure [m^3]"""

    vol_bz_liq: float = 0.0
    """Volume of liquid channels [m^3]"""

    vol_bz_liq_ib: float = 0.0
    """Volume of liquid channels [m^3]"""

    vol_bz_liq_ob: float = 0.0
    """Volume of liquid channels [m^3]"""

    vol_bw: float = 0.0
    """Volume of backwall [m^3]"""

    vol_bss: float = 0.0
    """Volume of MF/BSS [m^3]"""

    wht_cer: float = 0.0
    """BZ masses by composition [kg]"""

    wht_stl_struct: float = 0.0
    """BZ masses by composition [kg]"""

    wht_cool_struct: float = 0.0
    """BZ masses by composition [kg]"""

    wht_bw_stl: float = 0.0
    """Backwall masses by composition [kg]"""

    wht_bw_cool: float = 0.0
    """Backwall masses by composition [kg]"""

    wht_mfbss_stl: float = 0.0
    """MF/BSS masses by composition [kg]"""

    wht_mfbss_cool: float = 0.0
    """MF/BSS masses by composition [kg]"""

    wht_mfbss_pbli: float = 0.0
    """MF/BSS masses by composition [kg]"""

    fwmass_stl: float = 0.0
    """FW masses by composition [kg]"""

    fwmass_cool: float = 0.0
    """FW masses by composition [kg]"""

    mass_cool_blanket: float = 0.0
    """Total masses of material in blanket [kg]"""

    mass_liq_blanket: float = 0.0
    """Total masses of material in blanket [kg]"""

    mass_stl_blanket: float = 0.0
    """Total masses of material in blanket [kg]"""

    mass_segm_ib: float = 0.0
    """Total mass for an inboard/outboard reactor segment [kg]"""

    mass_segm_ob: float = 0.0
    """Total mass for an inboard/outboard reactor segment [kg]"""


CREATE_DICTS_FROM_DATACLASS = DCLLData
