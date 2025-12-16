impvardiv: int = None
"""Index of impurity to be iterated for Reinke divertor detachment criterion"""

lhat: float = None
"""Connection length factor L|| = lhat qstar R for Reinke criterion, default value from
Post et al. 1995 J. Nucl. Mat.  220-2 1014
"""

fzmin: float = None
"""Minimum impurity fraction necessary for detachment. This is the impurity at the SOL/Div."""

fzactual: float = None
"""Actual impurity fraction of divertor impurity (impvardiv) in the SoL (taking
impurity_enrichment into account) (`iteration variable 148`)
"""

reinke_mode: int = None
"""Switch for Reinke criterion H/I mode:
- =0 H-mode
- =1 I-mode
"""


def init_reinke_variables():
    """Initialise Reinke criterion variables"""
    global impvardiv, lhat, fzmin, fzactual, reinke_mode

    impvardiv = 9
    lhat = 4.33
    fzmin = 0.0
    fzactual = 0.001
    reinke_mode = 0
