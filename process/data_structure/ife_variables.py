"""Module containing global variables for the Inertial Fusion Energy (IFE) models



Default IFE builds and material volumes are those for the SOMBRERO device.
The 2-dimensional arrays have indices (region, material), where 'region'
is the region and maxmat is the 'material':

- 'region' = 1 radially outside chamber
- 'region' = 2 above chamber
- 'region' = 3 below chamber
"""

from dataclasses import dataclass, field

import numpy as np

from process.core.data_structure.parameter import Parameter, PROCESSModelData

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


@dataclass(slots=True)
class IFEData(PROCESSModelData):
    """Dataclass holding IFE variables"""

    bldr: Parameter[float] = 1.0
    """radial thickness of IFE blanket (m; calculated `if ifetyp=4`)"""

    bldrc: Parameter[float] = 1.0
    """radial thickness of IFE curtain (m; `ifetyp=4`)"""

    bldzl: Parameter[float] = 4.0
    """vertical thickness of IFE blanket below chamber (m)"""

    bldzu: Parameter[float] = 4.0
    """vertical thickness of IFE blanket above chamber (m)"""

    blmatf: Parameter[list[float]] = field(
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

    blmatm: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE blanket material masses (kg)"""

    blmatv: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE blanket material volumes (m3)"""

    blvol: Parameter[list[float]] = field(default_factory=lambda: np.zeros(3))
    """IFE blanket volume (m3)"""

    cdriv0: Parameter[float] = 154.3
    """IFE generic/laser driver cost at edrive=0 (M$)"""

    cdriv1: Parameter[float] = 163.2
    """IFE low energy heavy ion beam driver cost extrapolated to `edrive=0` (M$)"""

    cdriv2: Parameter[float] = 244.9
    """IFE high energy heavy ion beam driver cost extrapolated to `edrive=0` (M$)"""

    cdriv3: Parameter[float] = 1.463
    """IFE driver cost ($/J wall plug) (`ifedrv==3`)"""

    chdzl: Parameter[float] = 9.0
    """vertical thickness of IFE chamber below centre (m)"""

    chdzu: Parameter[float] = 9.0
    """vertical thickness of IFE chamber above centre (m)"""

    chmatf: Parameter[list[float]] = field(
        default_factory=lambda: np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    )
    """IFE chamber material fractions"""

    chmatm: Parameter[list[float]] = field(default_factory=lambda: np.zeros(MAXMAT + 1))
    """IFE chamber material masses (kg)"""

    chmatv: Parameter[list[float]] = field(default_factory=lambda: np.zeros(MAXMAT + 1))
    """IFE chamber material volumes (m3)"""

    chrad: Parameter[float] = 6.5
    """radius of IFE chamber (m) (`iteration variable 84`)"""

    chvol: Parameter[float] = 0.0
    """IFE chamber volume (m3)"""

    dcdrv0: Parameter[float] = 111.4
    """IFE generic/laser driver cost gradient (M$/MJ)"""

    dcdrv1: Parameter[float] = 78.0
    """HIB driver cost gradient at low energy (M$/MJ)"""

    dcdrv2: Parameter[float] = 59.9
    """HIB driver cost gradient at high energy (M$/MJ)"""

    drveff: Parameter[float] = 0.28
    """IFE driver wall plug to target efficiency (`ifedrv=0,3`) (`iteration variable 82`)"""

    edrive: Parameter[float] = 5.0e6
    """IFE driver energy (J) (`iteration variable 81`)"""

    etadrv: Parameter[float] = 0.0
    """IFE driver wall plug to target efficiency"""

    etali: Parameter[float] = 0.4
    """IFE lithium pump wall plug efficiency (`ifetyp=4`)"""

    etave: Parameter[list[float]] = field(
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

    fauxbop: Parameter[float] = 0.06
    """fraction of gross electric power to balance-of-plant (IFE)"""

    fbreed: Parameter[float] = 0.51
    """fraction of breeder external to device core"""

    fburn: Parameter[float] = 0.3333
    """IFE burn fraction (fraction of tritium fused/target)"""

    flirad: Parameter[float] = 0.78
    """radius of FLiBe/lithium inlet (m) (`ifetyp=3,4`)"""

    fwdr: Parameter[float] = 0.01
    """radial thickness of IFE first wall (m)"""

    fwdzl: Parameter[float] = 0.01
    """vertical thickness of IFE first wall below chamber (m)"""

    fwdzu: Parameter[float] = 0.01
    """vertical thickness of IFE first wall above chamber (m)"""

    fwmatf: Parameter[list[float]] = field(
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

    fwmatm: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE first wall material masses (kg)"""

    fwmatv: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE first wall material volumes (kg)"""

    fwvol: Parameter[list[float]] = field(default_factory=lambda: np.zeros(3))
    """IFE first wall volume (m3)"""

    gain: Parameter[float] = 0.0
    """IFE target gain"""

    gainve: Parameter[list[float]] = field(
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

    htpmw_ife: Parameter[float] = 0.0
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

    lipmw: Parameter[float] = 0.0
    """IFE lithium pump power (MW; `ifetyp=4`)"""

    mcdriv: Parameter[float] = 1.0
    """IFE driver cost multiplier"""

    mflibe: Parameter[float] = 0.0
    """total mass of FLiBe (kg)"""

    pdrive: Parameter[float] = 23.0e6
    """IFE driver power reaching target (W) (`iteration variable 85`)"""

    pfusife: Parameter[float] = 1000.0
    """IFE input fusion power (MW) (`ifedrv=3 only`; `itv 155`)"""

    pifecr: Parameter[float] = 10.0
    """IFE cryogenic power requirements (MW)"""

    ptargf: Parameter[float] = 2.0
    """IFE target factory power at 6 Hz repetition rate (MW)"""

    r1: Parameter[float] = 0.0
    """IFE device radial build (m)"""

    r2: Parameter[float] = 0.0
    """IFE device radial build (m)"""

    r3: Parameter[float] = 0.0
    """IFE device radial build (m)"""

    r4: Parameter[float] = 0.0
    """IFE device radial build (m)"""

    r5: Parameter[float] = 0.0
    """IFE device radial build (m)"""

    r6: Parameter[float] = 0.0
    """IFE device radial build (m)"""

    r7: Parameter[float] = 0.0
    """IFE device radial build (m)"""

    reprat: Parameter[float] = 0.0
    """IFE driver repetition rate (Hz)"""

    rrin: Parameter[float] = 6.0
    """Input IFE repetition rate (Hz) (`ifedrv=3 only`; `itv 156`)"""

    rrmax: Parameter[float] = 20.0
    """maximum IFE repetition rate (Hz)"""

    shdr: Parameter[float] = 1.7
    """radial thickness of IFE shield (m)"""

    shdzl: Parameter[float] = 5.0
    """vertical thickness of IFE shield below chamber (m)"""

    shdzu: Parameter[float] = 5.0
    """vertical thickness of IFE shield above chamber (m)"""

    shmatf: Parameter[list[float]] = field(
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

    shmatm: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE shield material masses (kg)"""

    shmatv: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE shield material volumes (kg)"""

    shvol: Parameter[list[float]] = field(default_factory=lambda: np.zeros(3))
    """IFE shield volume (m3)"""

    sombdr: Parameter[float] = 2.7
    """radius of cylindrical blanket section below chamber (`ifetyp=2`)"""

    somtdr: Parameter[float] = 2.7
    """radius of cylindrical blanket section above chamber (`ifetyp=2`)"""

    taufall: Parameter[float] = 0.0
    """Lithium Fall Time (s)"""

    tdspmw: Parameter[float] = 0.01
    """IFE target delivery system power (MW)"""

    tfacmw: Parameter[float] = 0.0
    """IFE target factory power (MW)"""

    tgain: Parameter[float] = 85.0
    """IFE target gain (if `ifedrv = 0`) (`iteration variable 83`)"""

    uccarb: Parameter[float] = 50.0
    """cost of carbon cloth ($/kg)"""

    ucconc: Parameter[float] = 0.1
    """cost of concrete ($/kg)"""

    ucflib: Parameter[float] = 84.0
    """cost of FLiBe ($/kg)"""

    uctarg: Parameter[float] = 0.3
    """cost of IFE target ($/target)"""

    v1dr: Parameter[float] = 0.0
    """radial thickness of IFE void between first wall and blanket (m)"""

    v1dzl: Parameter[float] = 0.0
    """vertical thickness of IFE void 1 below chamber (m)"""

    v1dzu: Parameter[float] = 0.0
    """vertical thickness of IFE void 1 above chamber (m)"""

    v1matf: Parameter[list[float]] = field(
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

    v1matm: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE void 1 material masses (kg)"""

    v1matv: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE void 1 material volumes (kg)"""

    v1vol: Parameter[list[float]] = field(default_factory=lambda: np.zeros(3))
    """IFE void 1 volume (m3)"""

    v2dr: Parameter[float] = 2.0
    """radial thickness of IFE void between blanket and shield (m)"""

    v2dzl: Parameter[float] = 7.0
    """vertical thickness of IFE void 2 below chamber (m)"""

    v2dzu: Parameter[float] = 7.0
    """vertical thickness of IFE void 2 above chamber (m)"""

    v2matf: Parameter[list[float]] = field(
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

    v2matm: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE void 2 material masses (kg)"""

    v2matv: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE void 2 material volumes (kg)"""

    v2vol: Parameter[list[float]] = field(default_factory=lambda: np.zeros(3))
    """IFE void 2 volume (m3)"""

    v3dr: Parameter[float] = 43.3
    """radial thickness of IFE void outside shield (m)"""

    v3dzl: Parameter[float] = 30.0
    """vertical thickness of IFE void 3 below chamber (m)"""

    v3dzu: Parameter[float] = 20.0
    """vertical thickness of IFE void 3 above chamber (m)"""

    v3matf: Parameter[list[float]] = field(
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

    v3matm: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE void 3 material masses (kg)"""

    v3matv: Parameter[list[float]] = field(
        default_factory=lambda: np.zeros((3, MAXMAT + 1))
    )
    """IFE void 3 material volumes (kg)"""

    v3vol: Parameter[list[float]] = field(default_factory=lambda: np.zeros(3))
    """IFE void 3 volume (m3)"""

    zl1: Parameter[float] = 0.0
    """IFE vertical build below centre (m)"""

    zl2: Parameter[float] = 0.0
    """IFE vertical build below centre (m)"""

    zl3: Parameter[float] = 0.0
    """IFE vertical build below centre (m)"""

    zl4: Parameter[float] = 0.0
    """IFE vertical build below centre (m)"""

    zl5: Parameter[float] = 0.0
    """IFE vertical build below centre (m)"""

    zl6: Parameter[float] = 0.0
    """IFE vertical build below centre (m)"""

    zl7: Parameter[float] = 0.0
    """IFE vertical build below centre (m)"""

    zu1: Parameter[float] = 0.0
    """IFE vertical build above centre (m)"""

    zu2: Parameter[float] = 0.0
    """IFE vertical build above centre (m)"""

    zu3: Parameter[float] = 0.0
    """IFE vertical build above centre (m)"""

    zu4: Parameter[float] = 0.0
    """IFE vertical build above centre (m)"""

    zu5: Parameter[float] = 0.0
    """IFE vertical build above centre (m)"""

    zu6: Parameter[float] = 0.0
    """IFE vertical build above centre (m)"""

    zu7: Parameter[float] = 0.0
    """IFE vertical build above centre (m)"""


CREATE_DICTS_FROM_DATACLASS = IFEData
