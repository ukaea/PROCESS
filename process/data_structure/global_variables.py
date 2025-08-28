icase: str = None
"""Power plant type"""

runtitle: str
"""A short descriptive title for the run"""

run_tests: int = None

verbose: int = None

maxcal: int = None
"""Maximum number of solver iterations"""

fileprefix: str = None
"""Input file path prefix"""

output_prefix: str = None
"""Output file path prefix"""

xlabel: str = None
"""Scan parameter description label"""

vlabel: str = None
"""Scan value name label"""

xlabel_2: str = None
"""Scan parameter description label (2nd dimension)"""

vlabel_2: str = None
"""Scan value name label (2nd dimension)"""

iscan_global: int = None

convergence_parameter: int = None
"""VMCON convergence parameter 'sum'"""


def init_global_variables():
    global icase
    global runtitle
    global run_tests
    global verbose
    global maxcal
    global fileprefix
    global output_prefix
    global xlabel
    global vlabel
    global xlabel_2
    global vlabel_2
    global iscan_global
    global convergence_parameter

    icase = "Steady-state tokamak model"
    runtitle = "Run Title (change this line using input variable 'runtitle')"
    run_tests = 0
    verbose = 0
    maxcal = 200
    fileprefix = ""
    output_prefix = ""
    xlabel = ""
    vlabel = ""
    xlabel_2 = ""
    vlabel_2 = ""
    iscan_global = 0
    convergence_parameter = 0.0
