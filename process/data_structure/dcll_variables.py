"""Module containing variables for the Dual Coolant Lead Lithium (DCLL) models

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
"""

from dataclasses import dataclass

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@dataclass(slots=True)
class DCLLData(PROCESSModelData):
    """Dataclass holding DCLL variables"""

    r_fci: Parameter[float] = 0.0
    """Radial thickness of FCIs [m]"""

    r_backwall: Parameter[float] = 0.0
    """Radial thickness of back wall [m]"""

    bz_r_ib: Parameter[float] = 0.0
    """Radial BZ thickness [m]"""

    bz_r_ob: Parameter[float] = 0.0
    """Radial BZ thickness [m]"""

    f_vol_stff_plates: Parameter[float] = 0.0
    """Structure/coolant compositional fractions"""

    f_vol_stl_bz_struct: Parameter[float] = 0.0
    """Structure/coolant compositional fractions"""

    f_vol_stl_back_wall: Parameter[float] = 0.0
    """Structure/coolant compositional fractions"""

    f_vol_stl_fw: Parameter[float] = 0.0
    """Structure/coolant compositional fractions"""

    f_vol_mfbss_stl: Parameter[float] = 0.0
    """MF/BSS compositional fractions"""

    f_vol_mfbss_he: Parameter[float] = 0.0
    """MF/BSS compositional fractions"""

    f_vol_mfbss_pbli: Parameter[float] = 0.0
    """MF/BSS compositional fractions"""

    vol_fci: Parameter[float] = 0.0
    """Volume of FCIs [m^3]"""

    vol_bz_struct: Parameter[float] = 0.0
    """Volume of other BZ structure [m^3]"""

    vol_bz_liq: Parameter[float] = 0.0
    """Volume of liquid channels [m^3]"""

    vol_bz_liq_ib: Parameter[float] = 0.0
    """Volume of liquid channels [m^3]"""

    vol_bz_liq_ob: Parameter[float] = 0.0
    """Volume of liquid channels [m^3]"""

    vol_bw: Parameter[float] = 0.0
    """Volume of back wall [m^3]"""

    vol_bss: Parameter[float] = 0.0
    """Volume of MF/BSS [m^3]"""

    wht_cer: Parameter[float] = 0.0
    """BZ masses by composition [kg]"""

    wht_stl_struct: Parameter[float] = 0.0
    """BZ masses by composition [kg]"""

    wht_cool_struct: Parameter[float] = 0.0
    """BZ masses by composition [kg]"""

    wht_bw_stl: Parameter[float] = 0.0
    """Back wall masses by composition [kg]"""

    wht_bw_cool: Parameter[float] = 0.0
    """Back wall masses by composition [kg]"""

    wht_mfbss_stl: Parameter[float] = 0.0
    """MF/BSS masses by composition [kg]"""

    wht_mfbss_cool: Parameter[float] = 0.0
    """MF/BSS masses by composition [kg]"""

    wht_mfbss_pbli: Parameter[float] = 0.0
    """MF/BSS masses by composition [kg]"""

    fwmass_stl: Parameter[float] = 0.0
    """FW masses by composition [kg]"""

    fwmass_cool: Parameter[float] = 0.0
    """FW masses by composition [kg]"""

    mass_cool_blanket: Parameter[float] = 0.0
    """Total masses of material in blanket [kg]"""

    mass_liq_blanket: Parameter[float] = 0.0
    """Total masses of material in blanket [kg]"""

    mass_stl_blanket: Parameter[float] = 0.0
    """Total masses of material in blanket [kg]"""

    mass_segm_ib: Parameter[float] = 0.0
    """Total mass for an inboard/outboard reactor segment [kg]"""

    mass_segm_ob: Parameter[float] = 0.0
    """Total mass for an inboard/outboard reactor segment [kg]"""


CREATE_DICTS_FROM_DATACLASS = DCLLData
