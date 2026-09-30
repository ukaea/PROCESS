"""Constants shared by summary plotting modules."""

from __future__ import annotations

import numpy as np

SOLENOID_COLOUR = ["pink", "#1764ab"]
CSCOMPRESSION_COLOUR = ["maroon", "#33CCCC"]
TFC_COLOUR = ["cyan", "#084a91"]
THERMAL_SHIELD_COLOUR = ["gray", "#e3eef9"]
VESSEL_COLOUR = ["green", "#b7d4ea"]
SHIELD_COLOUR = ["green", "#94c4df"]
BLANKET_COLOUR = ["magenta", "#4a98c9"]
PLASMA_COLOUR = ["khaki", "#cc8acc"]
CRYOSTAT_COLOUR = ["red", "#2e7ebc"]
FIRSTWALL_COLOUR = ["darkblue", "darkblue"]
NBSHIELD_COLOUR = ["black", "black"]
thin = 0.0
RADIAL_BUILD = [
    "dr_bore",
    "dr_cs",
    "dr_cs_precomp",
    "dr_cs_tf_gap",
    "dr_tf_inboard",
    "dr_tf_shld_gap",
    "dr_shld_thermal_inboard",
    "dr_shld_vv_gap_inboard",
    "dr_vv_inboard",
    "dr_shld_inboard",
    "vvblgapi",
    "dr_blkt_inboard",
    "dr_fw_inboard",
    "dr_fw_plasma_gap_inboard",
    "rminori",
    "rminoro",
    "dr_fw_plasma_gap_outboard",
    "dr_fw_outboard",
    "dr_blkt_outboard",
    "vvblgapo",
    "dr_shld_outboard",
    "dr_vv_outboard",
    "dr_shld_vv_gap_outboard",
    "dr_shld_thermal_outboard",
    "dr_tf_shld_gap",
    "dr_tf_outboard",
]
vertical_lower = [
    "z_plasma_xpoint_lower",
    "dz_xpoint_divertor",
    "dz_divertor",
    "dz_shld_lower",
    "dz_vv_lower",
    "dz_shld_vv_gap",
    "dz_shld_thermal",
    "dr_tf_shld_gap",
    "dr_tf_inboard",
]
ANIMATION_INFO = [
    ("rmajor", "Major radius", "m"),
    ("rminor", "Minor radius", "m"),
    ("aspect", "Aspect ratio", ""),
]
rtangle = np.pi / 2
rtangle2 = 2 * rtangle
white_box = {"boxstyle": "round", "facecolor": "white", "alpha": 1.0}


__all__ = [
    "ANIMATION_INFO",
    "BLANKET_COLOUR",
    "CRYOSTAT_COLOUR",
    "CSCOMPRESSION_COLOUR",
    "FIRSTWALL_COLOUR",
    "NBSHIELD_COLOUR",
    "PLASMA_COLOUR",
    "RADIAL_BUILD",
    "SHIELD_COLOUR",
    "SOLENOID_COLOUR",
    "TFC_COLOUR",
    "THERMAL_SHIELD_COLOUR",
    "VESSEL_COLOUR",
    "rtangle",
    "rtangle2",
    "thin",
    "vertical_lower",
    "white_box",
]
