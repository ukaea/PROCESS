from dataclasses import dataclass


@dataclass
class ReinkeData:
    impvardiv: int = 9
    """Index of impurity to be iterated for Reinke divertor detachment criterion"""

    lhat: float = 4.33
    """Connection length factor L|| = lhat qstar R for Reinke criterion, default value from
    Post et al. 1995 J. Nucl. Mat.  220-2 1014
    """

    fzmin: float = 0.0
    """Minimum impurity fraction necessary for detachment. This is the impurity at the SOL/Div."""

    fzactual: float = 0.001
    """Actual impurity fraction of divertor impurity (impvardiv) in the SoL (taking
    impurity_enrichment into account) (`iteration variable 148`)
    """

    reinke_mode: int = 0
    """Switch for Reinke criterion H/I mode:
    - =0 H-mode
    - =1 I-mode
    """


CREATE_DICTS_FROM_DATACLASS = ReinkeData
