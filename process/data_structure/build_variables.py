"""Module containing variables for the build models"""

from dataclasses import dataclass
from enum import IntEnum
from types import DynamicClassAttribute

from process.core.data_structure.parameter import Parameter, PROCESSModelData


class CSPrecompressionConfiguration(IntEnum):
    """Switch for CS pre-compression structure, controlled through `BuildData.i_cs_precomp`"""

    NO_CS_PRECOMPRESSION_STRUCTURE = 0
    CS_PRECOMPRESSION_STRUCTURE_PRESENT = 1


class InboardBlanketConfiguration(IntEnum):
    """Switch for inboard blanket, controlled through `BuildData.i_blkt_inboard`"""

    NO_INBOARD_BLANKET = 0
    INBOARD_BLANKET_PRESENT = 1


class TFCSRadialConfiguration(IntEnum):
    """Switch for placing the TF coil inside the CS, controlled through `BuildData.i_tf_inside_cs`"""

    TF_OUTSIDE_CS = (0, "Inboard TF coil leg is outside the CS")
    TF_INSIDE_CS = (1, "Inboard TF coil leg is inside the CS")

    def __new__(cls, value: int, description: str):
        """Create a new TFCSRadialConfiguration enum member with value and description.

        Parameters
        ----------
        value : int
            The numeric value of the enum member.
        description : str
            The description of the TFCSRadialConfiguration enum member.

        Returns
        -------
        TFCSRadialConfiguration
            A new enum member with the specified value and description.
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._description_ = description
        return obj

    @DynamicClassAttribute
    def description(self):
        """Description of the TF inboard and CS coil radial configuration."""
        return self._description_


@dataclass(slots=True)
class BuildData(PROCESSModelData):
    """Dataclass holding build variables"""

    rminor_min: Parameter[float] = 0.25
    """Minimum allowed minor radius (m)"""

    available_radial_space: Parameter[float] = 0.0
    """Minimal radial space between plasma and coils (m)"""

    a_blkt_total_surface: Parameter[float] = 0.0
    """blanket total surface area (m2)"""

    a_blkt_total_surface_full_coverage: Parameter[float] = 0.0
    """Blanket(s) total surface area with no holes or ports (toroidally continuous) (m²)"""

    a_blkt_inboard_surface: Parameter[float] = 0.0
    """inboard blanket surface area (m2)"""

    a_blkt_inboard_surface_full_coverage: Parameter[float] = 0.0
    """Inboard blanket surface area with no holes or ports (toroidally continuous) (m²)"""

    a_blkt_outboard_surface: Parameter[float] = 0.0
    """outboard blanket surface area (m2)"""

    a_blkt_outboard_surface_full_coverage: Parameter[float] = 0.0
    """Outboard blanket surface area with no holes or ports (toroidally continuous) (m²)"""

    blbmith: Parameter[float] = 0.17
    """inboard blanket box manifold thickness (m) (`blktmodel>0`)"""

    blbmoth: Parameter[float] = 0.27
    """outboard blanket box manifold thickness (m) (`blktmodel>0`)"""

    blbpith: Parameter[float] = 0.30
    """inboard blanket base plate thickness (m) (`blktmodel>0`)"""

    blbpoth: Parameter[float] = 0.35
    """outboard blanket base plate thickness (m) (`blktmodel>0`)"""

    blbuith: Parameter[float] = 0.365
    """inboard blanket breeding zone thickness (m) (`blktmodel>0`) (`iteration variable 90`)"""

    blbuoth: Parameter[float] = 0.465
    """outboard blanket breeding zone thickness (m) (`blktmodel>0`) (`iteration variable 91`)"""

    dr_blkt_inboard: Parameter[float] = 0.115
    """inboard blanket thickness (m); (calculated if `blktmodel>0`) (=0.0 if `i_blkt_inboard=0`)"""

    dr_blkt_outboard: Parameter[float] = 0.235
    """outboard blanket thickness (m); calculated if `blktmodel>0`"""

    dz_blkt_upper: Parameter[float] = 0.0
    """top blanket thickness (m), = mean of inboard and outboard blanket thicknesses"""

    dz_fw_upper: Parameter[float] = 0.0
    """upper first wall thickness (m)"""

    dr_bore: Parameter[float] = 1.42
    """Machine centre bore gap inboard radius (m) (`iteration variable 29`)"""

    dr_cs_bore: Parameter[float] = 1.42
    """Central solenoid bore radius [m]"""

    f_z_cryostat: Parameter[float] = 4.268
    """cryostat lid height scaling factor (tokamaks)"""

    dr_cryostat: Parameter[float] = 0.07
    """cryostat thickness (m)"""

    dr_vv_inboard: Parameter[float] = 0.07
    """vacuum vessel inboard thickness (TF coil / shield) (m)"""

    dr_vv_outboard: Parameter[float] = 0.07
    """vacuum vessel outboard thickness (TF coil / shield) (m)"""

    dz_vv_upper: Parameter[float] = 0.07
    """vacuum vessel topside thickness (TF coil / shield) (m) (= dz_vv_lower if double-null)"""

    dz_vv_lower: Parameter[float] = 0.07
    """vacuum vessel underside thickness (TF coil / shield) (m)"""

    dr_vv_shells: Parameter[float] = 0.12
    """vacuum vessel double walled shell thicknesses (m)"""

    fcspc: Parameter[float] = 0.6
    """Fraction of space occupied by CS pre-compression structure"""

    fseppc: Parameter[float] = 3.5e8
    """Separation force in CS coil pre-compression structure"""

    dr_fw_inboard: Parameter[float] = 0.0
    """inboard first wall thickness, initial estimate as calculated (m)"""

    dr_fw_outboard: Parameter[float] = 0.0
    """outboard first wall thickness, initial estimate as calculated (m)"""

    dr_shld_vv_gap_inboard: Parameter[float] = 0.155
    """gap between inboard vacuum vessel and thermal shield (m) (`iteration variable 61`)"""

    dr_cs_tf_gap: Parameter[float] = 0.08
    """gap between central solenoid and TF coil (m) (`iteration variable 42`)"""

    gapomin: Parameter[float] = 0.234
    """minimum gap between outboard vacuum vessel and TF coil (m) (`iteration variable 31`)"""

    dr_shld_vv_gap_outboard: Parameter[float] = 0.0
    """gap between outboard vacuum vessel and TF coil (m)"""

    z_tf_inside_half: Parameter[float] = 0.0
    """maximum (half-)height of TF coil (inside edge) (m)"""

    dz_tf_upper_lower_midplane: Parameter[float] = 0.0
    """difference in distance from midplane of upper and lower portions of TF
    legs (non-zero for single-null devices) (m)
    """

    z_tf_top: Parameter[float] = 0.0
    """height to top of (upper) TF coil leg (m)"""

    hr1: Parameter[float] = 0.0
    """half-height of TF coil inboard leg straight section (m)"""

    iohcl: int = 1
    """Switch for existence of central solenoid:
    - =0 central solenoid not present
    - =1 central solenoid exists
    """

    i_cs_precomp: int = 1
    """Switch for existence of central solenoid pre-compression structure:
    - =0 no pre-compression structure
    - =1 calculated pre-compression structure
    """

    i_tf_inside_cs: int = 0
    """Switch for placing the TF coil inside the CS
    - = 0 TF coil is outside the CS (default)
    - = 1 TF coil is inside the CS
    """

    dr_cs: Parameter[float] = 0.811
    """Central solenoid radial thickness (m) (`iteration variable 16`)"""

    dr_cs_precomp: Parameter[float] = 0.0
    """CS coil precompression structure radial thickness (m)"""

    rbld: Parameter[float] = 0.0
    """sum of thicknesses to the major radius (m)"""

    required_radial_space: Parameter[float] = 0.0
    """Required space between coil and plasma for blanket shield wall etc (m)"""

    rinboard: Parameter[float] = 0.651
    """plasma inboard radius (m) (`consistency equation 29`)"""

    r_shld_inboard_inner: Parameter[float] = 0.0
    """radius to inboard shield (inside point) (m)"""

    r_shld_outboard_outer: Parameter[float] = 0.0
    """radius to outboard shield (outside point) (m)"""

    r_vv_inboard_out: Parameter[float] = 0.0
    """Radial plasma facing side position of inboard vacuum vessel [m]"""

    r_sh_inboard_in: Parameter[float] = 0.0
    """Radial inner side position of inboard neutronic shield [m]"""

    r_sh_inboard_out: Parameter[float] = 0.0
    """Radial plasma facing side position of inboard neutronic shield [m]"""

    r_tf_inboard_in: Parameter[float] = 0.0
    """Mid-plane inboard TF coil leg radius at the centre-machine side [m]"""

    r_tf_inboard_mid: Parameter[float] = 0.0
    """Mid-plane inboard TF coil leg radius at middle of the coil [m]"""

    r_tf_inboard_out: Parameter[float] = 0.0
    """Mid-plane inboard TF coil leg radius at the plasma side [m]"""

    r_tf_outboard_mid: Parameter[float] = 0.0
    """Mid-plane outboard TF coil leg radius at the middle of the coil [m]"""

    i_r_cp_top: int = 0
    """Switch selecting the he parametrization of the outer radius of the top of the CP part of the TF coil
    0 : `r_cp_top` is set by the plasma shape
    1 : `r_cp_top` is a user input
    2 : `r_cp_top` is set using the CP top and midplane CP radius ratio
    """

    r_cp_top: Parameter[float] = 0.0
    """Top outer radius of the centropost (ST only) (m)"""

    f_r_cp: Parameter[float] = 1.4
    """Ratio between the top and the midplane TF CP outer radius [-]
    Not used by default (-1) must be larger than 1 otherwise
    """

    dr_tf_inner_bore: Parameter[float] = 0.0
    """TF coil horizontal inner dr_bore (m)"""

    dh_tf_inner_bore: Parameter[float] = 0.0
    """TF coil vertical inner dr_bore (m)"""

    dr_fw_plasma_gap_inboard: Parameter[float] = 0.14
    """Gap between plasma and first wall, inboard side (m) (if `i_plasma_wall_gap=1`)
    Iteration variable: ixc = 73
    Scan variable: nsweep = 58
    """

    dr_fw_plasma_gap_outboard: Parameter[float] = 0.15
    """Gap between plasma and first wall, outboard side (m) (if `i_plasma_wall_gap=1`)
    Iteration variable: ixc = 74
    Scan variable: nsweep = 59
    """

    a_shld_total_surface: Parameter[float] = 0.0
    """shield total surface area (m2)"""

    a_shld_inboard_surface: Parameter[float] = 0.0
    """inboard shield surface area (m2)"""

    a_shld_outboard_surface: Parameter[float] = 0.0
    """outboard shield surface area (m2)"""

    dr_shld_inboard: Parameter[float] = 0.69
    """inboard shield thickness (m) (`iteration variable 93`)"""

    dz_shld_lower: Parameter[float] = 0.7
    """lower (under divertor) shield thickness (m)"""

    dr_shld_outboard: Parameter[float] = 1.05
    """outboard shield thickness (m) (`iteration variable 94`)"""

    dz_shld_upper: Parameter[float] = 0.6
    """upper/lower shield thickness (m); calculated if `blktmodel > 0` (= dz_shld_lower if double-null)"""

    sigallpc: Parameter[float] = 3.0e8
    """allowable stress in CSpre-compression structure (Pa)"""

    # TODO: Issue #514 Make dr_tf_inboard an output not an iteration variable
    dr_tf_inboard: Parameter[float] = 0.0
    """inboard TF coil thickness, (centrepost for ST) (m)
    (input, calculated or `iteration variable 13`)
    """

    dz_tf_plasma_centre_offset: Parameter[float] = 0.0
    """Vertical distance between centre of TF coils and centre of plasma (m)"""

    f_dr_tf_outboard_inboard: Parameter[float] = 1.19
    """TF coil outboard leg / inboard leg radial thickness
    ratio (`i_tf_sup=0` only) (`iteration variable 75`)
    """

    dr_tf_outboard: Parameter[float] = 0.0
    """Outboard TF coil thickness (m)"""

    dr_tf_shld_gap: Parameter[float] = 0.05
    """Minimum metal-to-metal gap between TF coil and thermal shield (m)"""

    dr_shld_thermal_inboard: Parameter[float] = 0.05
    """TF-VV thermal shield thickness, inboard (m)"""

    dr_shld_thermal_outboard: Parameter[float] = 0.05
    """TF-VV thermal shield thickness, outboard (m)"""

    dz_shld_thermal: Parameter[float] = 0.05
    """TF-VV thermal shield thickness, vertical build (m)"""

    dz_shld_vv_gap: Parameter[float] = 0.163
    """vertical gap between vacuum vessel and thermal shields (m)"""

    dz_xpoint_divertor: Parameter[float] = 0.0
    """vertical gap between x-point and divertor (m) (if = 0, it is calculated)"""

    dz_fw_plasma_gap: Parameter[float] = 0.60
    """vertical gap between top of plasma and first wall (m) (= dz_xpoint_divertor if double-null)"""

    dr_shld_blkt_gap: Parameter[float] = 0.05
    """gap between vacuum vessel and blanket (m)"""

    plleni: Parameter[float] = 1.0
    """length of inboard divertor plate (m)"""

    plleno: Parameter[float] = 1.0
    """length of outboard divertor plate (m)"""

    plsepi: Parameter[float] = 1.0
    """poloidal length, x-point to inboard strike point (m)"""

    plsepo: Parameter[float] = 1.5
    """poloidal length, x-point to outboard strike point (m)"""

    rspo: Parameter[float] = 0.0
    """outboard strike point radius (m)"""

    z_plasma_xpoint_upper: Parameter[float] = 0.0
    """Vertical height of the upper plasma x-point (m)"""

    z_plasma_xpoint_lower: Parameter[float] = 0.0
    """Vertical height of the lower plasma x-point (m)"""

    ripflag: int = 0
    """1 if the fitted range of applicability is exceeded for
    the ripple calculation, else 0
    """

    i_blkt_inboard: int = 1
    """Switch for inboard blanket, controlled through InboardBlanketConfiguration enum:
    - =0 No inboard blanket (dr_blkt_inboard=0.0)
    - =1 Inboard blanket present
    """


CREATE_DICTS_FROM_DATACLASS = BuildData
