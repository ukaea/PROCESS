"""Module containing global variables relating to the inertial fusion energy model



Default IFE builds and material volumes are those for the SOMBRERO device.
The 2-dimensional arrays have indices (region, material), where 'region'
is the region and maxmat is the 'material':

- 'region' = 1 radially outside chamber
- 'region' = 2 above chamber
- 'region' = 3 below chamber
"""

from dataclasses import dataclass, field

import numpy as np

MAXMAT = 8
"""Total number of materials in IFE device. Material numbers are as follows:
- =0 void
- =1 steel
- =2 carbon cloth
- =3 FLiBe
- =4 lithium oxide Li2O
- =5 concrete
- =6 helium
- =7 xenon
- =8 lithium
"""


@dataclass
class IFEData:
    bldr: float = 1.0
    """radial thickness of IFE blanket (m; calculated `if ifetyp=4`)"""

    bldrc: float = 1.0
    """radial thickness of IFE curtain (m; `ifetyp=4`)"""

    bldzl: float = 4.0
    """vertical thickness of IFE blanket below chamber (m)"""

    bldzu: float = 4.0
    """vertical thickness of IFE blanket above chamber (m)"""

    blmatf: list[float] = field(
        default_factory=lambda: np.reshape(
            [
                0.05,
                0.05,
                0.05,
                0.0,
                0.0,
                0.0,
                0.45,
                0.45,
                0.45,
                0.0,
                0.0,
                0.0,
                0.20,
                0.20,
                0.20,
                0.0,
                0.0,
                0.0,
                0.30,
                0.30,
                0.30,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            (3, MAXMAT + 1),
        )
    )
    """IFE blanket material fractions"""

    blmatm: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE blanket material masses (kg)"""

    blmatv: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE blanket material volumes (m3)"""

    blvol: list[float] = field(default_factory=lambda: np.zeros(3))
    """IFE blanket volume (m3)"""

    cdriv0: float = 154.3
    """IFE generic/laser driver cost at edrive=0 (M$)"""

    cdriv1: float = 163.2
    """IFE low energy heavy ion beam driver cost extrapolated to `edrive=0` (M$)"""

    cdriv2: float = 244.9
    """IFE high energy heavy ion beam driver cost extrapolated to `edrive=0` (M$)"""

    cdriv3: float = 1.463
    """IFE driver cost ($/J wall plug) (`ifedrv==3`)"""

    chdzl: float = 9.0
    """vertical thickness of IFE chamber below centre (m)"""

    chdzu: float = 9.0
    """vertical thickness of IFE chamber above centre (m)"""

    chmatf: list[float] = field(
        default_factory=lambda: np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    )
    """IFE chamber material fractions"""

    chmatm: list[float] = field(default_factory=lambda: np.zeros(MAXMAT + 1))
    """IFE chamber material masses (kg)"""

    chmatv: list[float] = field(default_factory=lambda: np.zeros(MAXMAT + 1))
    """IFE chamber material volumes (m3)"""

    chrad: float = 6.5
    """radius of IFE chamber (m) (`iteration variable 84`)"""

    chvol: float = 0.0
    """IFE chamber volume (m3)"""

    dcdrv0: float = 111.4
    """IFE generic/laser driver cost gradient (M$/MJ)"""

    dcdrv1: float = 78.0
    """HIB driver cost gradient at low energy (M$/MJ)"""

    dcdrv2: float = 59.9
    """HIB driver cost gradient at high energy (M$/MJ)"""

    drveff: float = 0.28
    """IFE driver wall plug to target efficiency (`ifedrv=0,3`) (`iteration variable 82`)"""

    edrive: float = 5.0e6
    """IFE driver energy (J) (`iteration variable 81`)"""

    etadrv: float = 0.0
    """IFE driver wall plug to target efficiency"""

    etali: float = 0.4
    """IFE lithium pump wall plug efficiency (`ifetyp=4`)"""

    etave: list[float] = field(
        default_factory=lambda: np.array([
            0.082,
            0.079,
            0.076,
            0.073,
            0.069,
            0.066,
            0.062,
            0.059,
            0.055,
            0.051,
        ])
    )
    """IFE driver efficiency vs driver energy (`ifedrv=-1`)"""

    fauxbop: float = 0.06
    """fraction of gross electric power to balance-of-plant (IFE)"""

    fbreed: float = 0.51
    """fraction of breeder external to device core"""

    fburn: float = 0.3333
    """IFE burn fraction (fraction of tritium fused/target)"""

    flirad: float = 0.78
    """radius of FLiBe/lithium inlet (m) (`ifetyp=3,4`)"""

    fwdr: float = 0.01
    """radial thickness of IFE first wall (m)"""

    fwdzl: float = 0.01
    """vertical thickness of IFE first wall below chamber (m)"""

    fwdzu: float = 0.01
    """vertical thickness of IFE first wall above chamber (m)"""

    fwmatf: list[float] = field(
        default_factory=lambda: np.reshape(
            [
                0.05,
                0.05,
                0.05,
                0.0,
                0.0,
                0.0,
                0.95,
                0.95,
                0.95,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            (3, MAXMAT + 1),
        )
    )
    """IFE first wall material fractions"""

    fwmatm: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE first wall material masses (kg)"""

    fwmatv: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE first wall material volumes (kg)"""

    fwvol: list[float] = field(default_factory=lambda: np.zeros(3))
    """IFE first wall volume (m3)"""

    gain: float = 0.0
    """IFE target gain"""

    gainve: list[float] = field(
        default_factory=lambda: np.array([
            60.0,
            95.0,
            115.0,
            125.0,
            133.0,
            141.0,
            152.0,
            160.0,
            165.0,
            170.0,
        ])
    )
    """IFE target gain vs driver energy (`ifedrv=-1`)"""

    htpmw_ife: float = 0.0
    """IFE heat transport system electrical pump power (MW)"""

    ife: int = 0
    """Switch for IFE option:
    - =0 use tokamak, RFP or stellarator model
    - =1 use IFE model
    """

    ifedrv: int = 2
    """Switch for type of IFE driver:
    - =-1 use gainve, etave for gain and driver efficiency
    - =0 use tgain, drveff for gain and driver efficiency
    - =1 use laser driver based on SOMBRERO design
    - =2 use heavy ion beam driver based on OSIRIS
    - =3 Input pfusife, rrin and drveff
    """

    ifetyp: int = 0
    """Switch for type of IFE device build:
    - =0 generic (cylindrical) build
    - =1 OSIRIS-like build
    - =2 SOMBRERO-like build
    - =3 HYLIFE-II-like build
    - =4 2019 build
    """

    lipmw: float = 0.0
    """IFE lithium pump power (MW; `ifetyp=4`)"""

    mcdriv: float = 1.0
    """IFE driver cost multiplier"""

    mflibe: float = 0.0
    """total mass of FLiBe (kg)"""

    pdrive: float = 23.0e6
    """IFE driver power reaching target (W) (`iteration variable 85`)"""

    pfusife: float = 1000.0
    """IFE input fusion power (MW) (`ifedrv=3 only`; `itv 155`)"""

    pifecr: float = 10.0
    """IFE cryogenic power requirements (MW)"""

    ptargf: float = 2.0
    """IFE target factory power at 6 Hz repetition rate (MW)"""

    r1: float = 0.0
    """IFE device radial build (m)"""

    r2: float = 0.0
    """IFE device radial build (m)"""

    r3: float = 0.0
    """IFE device radial build (m)"""

    r4: float = 0.0
    """IFE device radial build (m)"""

    r5: float = 0.0
    """IFE device radial build (m)"""

    r6: float = 0.0
    """IFE device radial build (m)"""

    r7: float = 0.0
    """IFE device radial build (m)"""

    reprat: float = 0.0
    """IFE driver repetition rate (Hz)"""

    rrin: float = 6.0
    """Input IFE repetition rate (Hz) (`ifedrv=3 only`; `itv 156`)"""

    rrmax: float = 20.0
    """maximum IFE repetition rate (Hz)"""

    shdr: float = 1.7
    """radial thickness of IFE shield (m)"""

    shdzl: float = 5.0
    """vertical thickness of IFE shield below chamber (m)"""

    shdzu: float = 5.0
    """vertical thickness of IFE shield above chamber (m)"""

    shmatf: list[float] = field(
        default_factory=lambda: np.reshape(
            [
                0.05,
                0.05,
                0.05,
                0.19,
                0.19,
                0.19,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.665,
                0.665,
                0.665,
                0.095,
                0.095,
                0.095,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            (3, MAXMAT + 1),
        )
    )
    """IFE shield material fractions"""

    shmatm: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE shield material masses (kg)"""

    shmatv: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE shield material volumes (kg)"""

    shvol: list[float] = field(default_factory=lambda: np.zeros(3))
    """IFE shield volume (m3)"""

    sombdr: float = 2.7
    """radius of cylindrical blanket section below chamber (`ifetyp=2`)"""

    somtdr: float = 2.7
    """radius of cylindrical blanket section above chamber (`ifetyp=2`)"""

    taufall: float = 0.0
    """Lithium Fall Time (s)"""

    tdspmw: float = 0.01
    """IFE target delivery system power (MW)"""

    tfacmw: float = 0.0
    """IFE target factory power (MW)"""

    tgain: float = 85.0
    """IFE target gain (if `ifedrv = 0`) (`iteration variable 83`)"""

    uccarb: float = 50.0
    """cost of carbon cloth ($/kg)"""

    ucconc: float = 0.1
    """cost of concrete ($/kg)"""

    ucflib: float = 84.0
    """cost of FLiBe ($/kg)"""

    uctarg: float = 0.3
    """cost of IFE target ($/target)"""

    v1dr: float = 0.0
    """radial thickness of IFE void between first wall and blanket (m)"""

    v1dzl: float = 0.0
    """vertical thickness of IFE void 1 below chamber (m)"""

    v1dzu: float = 0.0
    """vertical thickness of IFE void 1 above chamber (m)"""

    v1matf: list[float] = field(
        default_factory=lambda: np.reshape(
            [
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            (3, MAXMAT + 1),
        )
    )
    """IFE void 1 material fractions"""

    v1matm: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE void 1 material masses (kg)"""

    v1matv: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE void 1 material volumes (kg)"""

    v1vol: list[float] = field(default_factory=lambda: np.zeros(3))
    """IFE void 1 volume (m3)"""

    v2dr: float = 2.0
    """radial thickness of IFE void between blanket and shield (m)"""

    v2dzl: float = 7.0
    """vertical thickness of IFE void 2 below chamber (m)"""

    v2dzu: float = 7.0
    """vertical thickness of IFE void 2 above chamber (m)"""

    v2matf: list[float] = field(
        default_factory=lambda: np.reshape(
            [
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            (3, MAXMAT + 1),
        )
    )
    """IFE void 2 material fractions"""

    v2matm: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE void 2 material masses (kg)"""

    v2matv: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE void 2 material volumes (kg)"""

    v2vol: list[float] = field(default_factory=lambda: np.zeros(3))
    """IFE void 2 volume (m3)"""

    v3dr: float = 43.3
    """radial thickness of IFE void outside shield (m)"""

    v3dzl: float = 30.0
    """vertical thickness of IFE void 3 below chamber (m)"""

    v3dzu: float = 20.0
    """vertical thickness of IFE void 3 above chamber (m)"""

    v3matf: list[float] = field(
        default_factory=lambda: np.reshape(
            [
                1.0,
                1.0,
                1.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
            (3, MAXMAT + 1),
        )
    )
    """IFE void 3 material fractions"""

    v3matm: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE void 3 material masses (kg)"""

    v3matv: list[float] = field(default_factory=lambda: np.zeros((3, MAXMAT + 1)))
    """IFE void 3 material volumes (kg)"""

    v3vol: list[float] = field(default_factory=lambda: np.zeros(3))
    """IFE void 3 volume (m3)"""

    zl1: float = 0.0
    """IFE vertical build below centre (m)"""

    zl2: float = 0.0
    """IFE vertical build below centre (m)"""

    zl3: float = 0.0
    """IFE vertical build below centre (m)"""

    zl4: float = 0.0
    """IFE vertical build below centre (m)"""

    zl5: float = 0.0
    """IFE vertical build below centre (m)"""

    zl6: float = 0.0
    """IFE vertical build below centre (m)"""

    zl7: float = 0.0
    """IFE vertical build below centre (m)"""

    zu1: float = 0.0
    """IFE vertical build above centre (m)"""

    zu2: float = 0.0
    """IFE vertical build above centre (m)"""

    zu3: float = 0.0
    """IFE vertical build above centre (m)"""

    zu4: float = 0.0
    """IFE vertical build above centre (m)"""

    zu5: float = 0.0
    """IFE vertical build above centre (m)"""

    zu6: float = 0.0
    """IFE vertical build above centre (m)"""

    zu7: float = 0.0
    """IFE vertical build above centre (m)"""


CREATE_DICTS_FROM_DATACLASS = IFEData
