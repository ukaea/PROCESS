"""Module containing global variables relating to the inertial fusion energy model

author: S. Muldrew (UKAEA)

Default IFE builds and material volumes are those for the SOMBRERO device.
The 2-dimensional arrays have indices (region, material), where 'region'
is the region and maxmat is the 'material':

- 'region' = 1 radially outside chamber
- 'region' = 2 above chamber
- 'region' = 3 below chamber
"""

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

bldr: float = None
"""radial thickness of IFE blanket (m; calculated `if ifetyp=4`)"""


bldrc: float = None
"""radial thickness of IFE curtain (m; `ifetyp=4`)"""


bldzl: float = None
"""vertical thickness of IFE blanket below chamber (m)"""


bldzu: float = None
"""vertical thickness of IFE blanket above chamber (m)"""


blmatf: list[float] = None
"""IFE blanket material fractions"""


blmatm: list[float] = None
"""IFE blanket material masses (kg)"""

blmatv: list[float] = None
"""IFE blanket material volumes (m3)"""

blvol: list[float] = None
"""IFE blanket volume (m3)"""


cdriv0: float = None
"""IFE generic/laser driver cost at edrive=0 (M$)"""


cdriv1: float = None
"""IFE low energy heavy ion beam driver cost extrapolated to `edrive=0` (M$)"""


cdriv2: float = None
"""IFE high energy heavy ion beam driver cost extrapolated to `edrive=0` (M$)"""


cdriv3: float = None
"""IFE driver cost ($/J wall plug) (`ifedrv==3`)"""


chdzl: float = None
"""vertical thickness of IFE chamber below centre (m)"""


chdzu: float = None
"""vertical thickness of IFE chamber above centre (m)"""


chmatf: list[float] = None
"""IFE chamber material fractions"""


chmatm: list[float] = None
"""IFE chamber material masses (kg)"""


chmatv: list[float] = None
"""IFE chamber material volumes (m3)"""


chrad: float = None
"""radius of IFE chamber (m) (`iteration variable 84`)"""


chvol: float = None
"""IFE chamber volume (m3)"""


dcdrv0: float = None
"""IFE generic/laser driver cost gradient (M$/MJ)"""


dcdrv1: float = None
"""HIB driver cost gradient at low energy (M$/MJ)"""


dcdrv2: float = None
"""HIB driver cost gradient at high energy (M$/MJ)"""


drveff: float = None
"""IFE driver wall plug to target efficiency (`ifedrv=0,3`) (`iteration variable 82`)"""


edrive: float = None
"""IFE driver energy (J) (`iteration variable 81`)"""


etadrv: float = None
"""IFE driver wall plug to target efficiency"""


etali: float = None
"""IFE lithium pump wall plug efficiency (`ifetyp=4`)"""


etave: list[float] = None
"""IFE driver efficiency vs driver energy (`ifedrv=-1`)"""


fauxbop: float = None
"""fraction of gross electric power to balance-of-plant (IFE)"""


fbreed: float = None
"""fraction of breeder external to device core"""


fburn: float = None
"""IFE burn fraction (fraction of tritium fused/target)"""


flirad: float = None
"""radius of FLiBe/lithium inlet (m) (`ifetyp=3,4`)"""


frrmax: float = None
"""f-value for maximum IFE repetition rate (`constraint equation 50`, `iteration variable 86`)"""


fwdr: float = None
"""radial thickness of IFE first wall (m)"""


fwdzl: float = None
"""vertical thickness of IFE first wall below chamber (m)"""


fwdzu: float = None
"""vertical thickness of IFE first wall above chamber (m)"""


fwmatf: list[float] = None
"""IFE first wall material fractions"""


fwmatm: list[float] = None
"""IFE first wall material masses (kg)"""


fwmatv: list[float] = None
"""IFE first wall material volumes (kg)"""


fwvol: list[float] = None
"""IFE first wall volume (m3)"""


gain: float = None
"""IFE target gain"""


gainve: list[float] = None
"""IFE target gain vs driver energy (`ifedrv=-1`)"""


htpmw_ife: float = None
"""IFE heat transport system electrical pump power (MW)"""


ife: int = None
"""Switch for IFE option:
- =0 use tokamak, RFP or stellarator model
- =1 use IFE model
"""


ifedrv: int = None
"""Switch for type of IFE driver:
- =-1 use gainve, etave for gain and driver efficiency
- =0 use tgain, drveff for gain and driver efficiency
- =1 use laser driver based on SOMBRERO design
- =2 use heavy ion beam driver based on OSIRIS
- =3 Input pfusife, rrin and drveff
"""


ifetyp: int = None
"""Switch for type of IFE device build:
- =0 generic (cylindrical) build
- =1 OSIRIS-like build
- =2 SOMBRERO-like build
- =3 HYLIFE-II-like build
- =4 2019 build
"""


lipmw: float = None
"""IFE lithium pump power (MW; `ifetyp=4`)"""


mcdriv: float = None
"""IFE driver cost multiplier"""


mflibe: float = None
"""total mass of FLiBe (kg)"""


pdrive: float = None
"""IFE driver power reaching target (W) (`iteration variable 85`)"""


pfusife: float = None
"""IFE input fusion power (MW) (`ifedrv=3 only`; `itv 155`)"""


pifecr: float = None
"""IFE cryogenic power requirements (MW)"""


ptargf: float = None
"""IFE target factory power at 6 Hz repetition rate (MW)"""


r1: float = None
"""IFE device radial build (m)"""


r2: float = None
"""IFE device radial build (m)"""


r3: float = None
"""IFE device radial build (m)"""


r4: float = None
"""IFE device radial build (m)"""


r5: float = None
"""IFE device radial build (m)"""


r6: float = None
"""IFE device radial build (m)"""


r7: float = None
"""IFE device radial build (m)"""


reprat: float = None
"""IFE driver repetition rate (Hz)"""


rrin: float = None
"""Input IFE repetition rate (Hz) (`ifedrv=3 only`; `itv 156`)"""


rrmax: float = None
"""maximum IFE repetition rate (Hz)"""


shdr: float = None
"""radial thickness of IFE shield (m)"""


shdzl: float = None
"""vertical thickness of IFE shield below chamber (m)"""


shdzu: float = None
"""vertical thickness of IFE shield above chamber (m)"""


shmatf: list[float] = None
"""IFE shield material fractions"""


shmatm: list[float] = None
"""IFE shield material masses (kg)"""


shmatv: list[float] = None
"""IFE shield material volumes (kg)"""


shvol: list[float] = None
"""IFE shield volume (m3)"""


sombdr: float = None
"""radius of cylindrical blanket section below chamber (`ifetyp=2`)"""


somtdr: float = None
"""radius of cylindrical blanket section above chamber (`ifetyp=2`)"""


taufall: float = None
"""Lithium Fall Time (s)"""


tdspmw: float = None
"""IFE target delivery system power (MW)"""


tfacmw: float = None
"""IFE target factory power (MW)"""


tgain: float = None
"""IFE target gain (if `ifedrv = 0`) (`iteration variable 83`)"""


uccarb: float = None
"""cost of carbon cloth ($/kg)"""


ucconc: float = None
"""cost of concrete ($/kg)"""


ucflib: float = None
"""cost of FLiBe ($/kg)"""


uctarg: float = None
"""cost of IFE target ($/target)"""


v1dr: float = None
"""radial thickness of IFE void between first wall and blanket (m)"""


v1dzl: float = None
"""vertical thickness of IFE void 1 below chamber (m)"""


v1dzu: float = None
"""vertical thickness of IFE void 1 above chamber (m)"""


v1matf: list[float] = None
"""IFE void 1 material fractions"""


v1matm: list[float] = None
"""IFE void 1 material masses (kg)"""


v1matv: list[float] = None
"""IFE void 1 material volumes (kg)"""


v1vol: list[float] = None
"""IFE void 1 volume (m3)"""


v2dr: float = None
"""radial thickness of IFE void between blanket and shield (m)"""


v2dzl: float = None
"""vertical thickness of IFE void 2 below chamber (m)"""


v2dzu: float = None
"""vertical thickness of IFE void 2 above chamber (m)"""


v2matf: list[float] = None
"""IFE void 2 material fractions"""


v2matm: list[float] = None
"""IFE void 2 material masses (kg)"""


v2matv: list[float] = None
"""IFE void 2 material volumes (kg)"""


v2vol: list[float] = None
"""IFE void 2 volume (m3)"""


v3dr: float = None
"""radial thickness of IFE void outside shield (m)"""


v3dzl: float = None
"""vertical thickness of IFE void 3 below chamber (m)"""


v3dzu: float = None
"""vertical thickness of IFE void 3 above chamber (m)"""


v3matf: list[float] = None
"""IFE void 3 material fractions"""


v3matm: list[float] = None
"""IFE void 3 material masses (kg)"""


v3matv: list[float] = None
"""IFE void 3 material volumes (kg)"""


v3vol: list[float] = None
"""IFE void 3 volume (m3)"""


zl1: float = None
"""IFE vertical build below centre (m)"""


zl2: float = None
"""IFE vertical build below centre (m)"""


zl3: float = None
"""IFE vertical build below centre (m)"""


zl4: float = None
"""IFE vertical build below centre (m)"""


zl5: float = None
"""IFE vertical build below centre (m)"""


zl6: float = None
"""IFE vertical build below centre (m)"""


zl7: float = None
"""IFE vertical build below centre (m)"""


zu1: float = None
"""IFE vertical build above centre (m)"""


zu2: float = None
"""IFE vertical build above centre (m)"""


zu3: float = None
"""IFE vertical build above centre (m)"""


zu4: float = None
"""IFE vertical build above centre (m)"""


zu5: float = None
"""IFE vertical build above centre (m)"""


zu6: float = None
"""IFE vertical build above centre (m)"""


zu7: float = None
"""IFE vertical build above centre (m)"""


def init_ife_variables():
    global bldr
    global bldrc
    global bldzl
    global bldzu
    global blmatf
    global blmatm
    global blmatv
    global blvol
    global cdriv0
    global cdriv1
    global cdriv2
    global cdriv3
    global chdzl
    global chdzu
    global chmatf
    global chmatm
    global chmatv
    global chrad
    global chvol
    global dcdrv0
    global dcdrv1
    global dcdrv2
    global drveff
    global edrive
    global etadrv
    global etali
    global etave
    global fauxbop
    global fbreed
    global fburn
    global flirad
    global frrmax
    global fwdr
    global fwdzl
    global fwdzu
    global fwmatf
    global fwmatm
    global fwmatv
    global fwvol
    global gain
    global gainve
    global htpmw_ife
    global ife
    global ifedrv
    global ifetyp
    global lipmw
    global mcdriv
    global mflibe
    global pdrive
    global pfusife
    global pifecr
    global ptargf
    global r1
    global r2
    global r3
    global r4
    global r5
    global r6
    global r7
    global reprat
    global rrin
    global rrmax
    global shdr
    global shdzl
    global shdzu
    global shmatf
    global shmatm
    global shmatv
    global shvol
    global sombdr
    global somtdr
    global taufall
    global tdspmw
    global tfacmw
    global tgain
    global uccarb
    global ucconc
    global ucflib
    global uctarg
    global v1dr
    global v1dzl
    global v1dzu
    global v1matf
    global v1matm
    global v1matv
    global v1vol
    global v2dr
    global v2dzl
    global v2dzu
    global v2matf
    global v2matm
    global v2matv
    global v2vol
    global v3dr
    global v3dzl
    global v3dzu
    global v3matf
    global v3matm
    global v3matv
    global v3vol
    global zl1
    global zl2
    global zl3
    global zl4
    global zl5
    global zl6
    global zl7
    global zu1
    global zu2
    global zu3
    global zu4
    global zu5
    global zu6
    global zu7

    """Initialise IFE variables"""
    bldr = 1.0
    bldrc = 1.0
    bldzl = 4.0
    bldzu = 4.0
    blmatf = np.reshape(
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
    blmatm = np.zeros((3, MAXMAT + 1))
    blmatv = np.zeros((3, MAXMAT + 1))
    blvol = np.zeros(3)
    cdriv0 = 154.3
    cdriv1 = 163.2
    cdriv2 = 244.9
    cdriv3 = 1.463
    chdzl = 9.0
    chdzu = 9.0
    chmatf = np.array([1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    chmatm = np.zeros(MAXMAT + 1)
    chmatv = np.zeros(MAXMAT + 1)
    chrad = 6.5
    chvol = 0.0
    dcdrv0 = 111.4
    dcdrv1 = 78.0
    dcdrv2 = 59.9
    drveff = 0.28
    edrive = 5.0e6
    etadrv = 0.0
    etali = 0.4
    etave = np.array([
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
    fauxbop = 0.06
    fbreed = 0.51
    fburn = 0.3333
    flirad = 0.78
    frrmax = 1.0
    fwdr = 0.01
    fwdzl = 0.01
    fwdzu = 0.01
    fwmatf = np.reshape(
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
    fwmatm = np.zeros((3, MAXMAT + 1))
    fwmatv = np.zeros((3, MAXMAT + 1))
    fwvol = np.zeros(3)
    gain = 0.0
    gainve = np.array([
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
    htpmw_ife = 0.0
    ife = 0
    ifedrv = 2
    ifetyp = 0
    lipmw = 0.0
    mcdriv = 1.0
    mflibe = 0.0
    pdrive = 23.0e6
    pfusife = 1000.0
    pifecr = 10.0
    ptargf = 2.0
    r1 = 0.0
    r2 = 0.0
    r3 = 0.0
    r4 = 0.0
    r5 = 0.0
    r6 = 0.0
    r7 = 0.0
    reprat = 0.0
    rrin = 6.0
    rrmax = 20.0
    shdr = 1.7
    shdzl = 5.0
    shdzu = 5.0
    shmatf = np.reshape(
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
    shmatm = np.zeros((3, MAXMAT + 1))
    shmatv = np.zeros((3, MAXMAT + 1))
    shvol = np.zeros(3)
    sombdr = 2.7
    somtdr = 2.7
    taufall = 0.0
    tdspmw = 0.01
    tfacmw = 0.0
    tgain = 85.0
    uccarb = 50.0
    ucconc = 0.1
    ucflib = 84.0
    uctarg = 0.3
    v1dr = 0.0
    v1dzl = 0.0
    v1dzu = 0.0
    v1matf = np.reshape(
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
    v1matm = np.zeros((3, MAXMAT + 1))
    v1matv = np.zeros((3, MAXMAT + 1))
    v1vol = np.zeros(3)
    v2dr = 2.0
    v2dzl = 7.0
    v2dzu = 7.0
    v2matf = np.reshape(
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
    v2matm = np.zeros((3, MAXMAT + 1))
    v2matv = np.zeros((3, MAXMAT + 1))
    v2vol = np.zeros(3)
    v3dr = 43.3
    v3dzl = 30.0
    v3dzu = 20.0
    v3matf = np.reshape(
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
    v3matm = np.zeros((3, MAXMAT + 1))
    v3matv = np.zeros((3, MAXMAT + 1))
    v3vol = np.zeros(3)
    zl1 = 0.0
    zl2 = 0.0
    zl3 = 0.0
    zl4 = 0.0
    zl5 = 0.0
    zl6 = 0.0
    zl7 = 0.0
    zu1 = 0.0
    zu2 = 0.0
    zu3 = 0.0
    zu4 = 0.0
    zu5 = 0.0
    zu6 = 0.0
    zu7 = 0.0
