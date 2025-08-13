import numpy as np

nef: int = None

nfxf: int = None

ricpf: float = None

ssq0: float = None

sig_axial: float = None

sig_hoop: float = None

axial_force: float = None

rfxf: list[float] = None

zfxf: list[float] = None

cfxf: list[float] = None

xind: list[float] = None

rcls: list[float] = None

zcls: list[float] = None

ccls: list[float] = None

ccl0: list[float] = None

bpf2: list[float] = None

vsdum: list[float] = None

first_call: bool = None

cslimit: bool = None


def init_pfcoil_module():
    global first_call
    global cslimit
    global nef
    global nfxf
    global ricpf
    global ssq0
    global sig_axial
    global sig_hoop
    global axial_force
    global rfxf
    global zfxf

    global cfxf
    global xind
    global rcls
    global zcls
    global ccls
    global ccl0
    global bpf2
    global vsdum

    first_call = True
    cslimit = False
    nef = 0
    nfxf = 0
    ricpf = 0.0
    ssq0 = 0.0
    sig_axial = 0.0
    sig_hoop = 0.0
    axial_force = 0

    rfxf = np.zeros(64)
    zfxf = np.zeros(64)
    cfxf = np.zeros(64)
    xind = np.zeros(64)
    rcls = np.zeros((10, 2))
    zcls = np.zeros((10, 2))
    ccls = np.zeros(10)
    ccl0 = np.zeros(10)
    bpf2 = np.zeros(22)
    vsdum = np.zeros((22, 3))
