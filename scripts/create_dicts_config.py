# ifail value of a successful process run
IFAIL_SUCCESS = 1

# default values for making a plot file from MFILE.DAT
PARAMETER_DEFAULTS = [
    "rmajor",
    "aspect",
    "rminor",
    "bt",
    "powfmw",
    "pnetelmw",
    "te",
    "pdivt",
    "sig_tf_case",
    "sig_tf_wp",
]

# parameters that start with f, but are not f-values
NON_F_VALUES = ["fcohbop", "fvsbrnni", "feffcd", "fcutfsu"]

# PROCESS TF Coil types
DICT_TF_TYPE = {
    1: "Nb3Sn ITER",
    2: "Bi-2212",
    3: "NbTi",
    4: "Nb3Sn user",
    5: "Nb3Sn WST",
    6: "REBCO Croco",
    7: "NbTi Ginzburg-Landau",
    8: "REBCO Ginzburg-Landau",
    9: "REBCO Hazelton-Zhai",
}

# FIMP Values
DICT_FIMP = {
    "fimp(1)": "Hydrogen (fraction calculated by code)",
    "fimp(2)": "Helium",
    "fimp(3)": "Beryllium",
    "fimp(4)": "Carbon",
    "fimp(5)": "Nitrogen",
    "fimp(6)": "Oxygen",
    "fimp(7)": "Neon",
    "fimp(8)": "Silicon",
    "fimp(9)": "Argon",
    "fimp(10)": "Iron",
    "fimp(11)": "Nickel",
    "fimp(12)": "Krypton",
    "fimp(13)": "Xenon",
    "fimp(14)": "Tungsten",
}

# Optimisation variable dictionary
DICT_OPTIMISATION_VARS = {
    1: "Plasma major radius",
    2: "ratio fusion power:input power",
    3: "neutron wall load",
    4: "total TF + PF coil power",
    5: "ratio fusion power:injection power",
    6: "cost of electricity",
    7: "constructed cost",
    8: "aspect ratio",
    9: "divertor heat load",
    10: "toroidal field on axis",
    11: "injection power",
    12: "hydrogen production capital cost",
    13: "hydrogen production rate",
    14: "pulse length",
    15: "plant availability factor",
    16: "linear combination of major radius (minimised) and pulse length (maximised)",
    17: "net electrical output",
}
