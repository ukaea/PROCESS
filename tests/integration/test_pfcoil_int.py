"""Integration tests for the PFCoil module.

Vaguely realistic mocked values are taken from baseline2019 output, init values
in the pfcoil_variables module, or where necessary, guesses.

Many of these subroutines are long and perform multiple gets/sets on many "use"
dependencies. As a result, many mocks are required to isolate the tests. There
are also many variables that could be asserted, so a few key variables central
to the testing of the subroutine have been chosen.
"""

import numpy as np
from numpy.testing import assert_array_almost_equal
import pytest
from process.fortran import pfcoil_module as pf
from process.fortran import build_variables as bv
from process.fortran import pfcoil_variables as pfv
from process.fortran import physics_variables as pv
from process.fortran import error_handling as eh
from process.fortran import fwbs_variables as fwbsv
from process.fortran import tfcoil_variables as tfv
from process.fortran import times_variables as tv
from process.fortran import constants
from process.pfcoil import PFCoil, mtrx, fixb
from process.cs_fatigue import CsFatigue


@pytest.fixture
def pfcoil():
    """Fixture to create PFCoil object.

    :return: a PFCoil instance
    :rtype: process.pfcoil.PFCoil
    """
    return PFCoil(cs_fatigue=CsFatigue())


def test_pfcoil(monkeypatch, pfcoil):
    """Test pfcoil subroutine.

    :param monkeypatch: mocking fixture
    :type monkeypatch: MonkeyPatch
    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    """

    monkeypatch.setattr(bv, "iohcl", 1)
    monkeypatch.setattr(bv, "hpfdif", 0.0)
    monkeypatch.setattr(bv, "hpfu", 4.0)  # guess
    monkeypatch.setattr(bv, "hmax", 8.8)
    monkeypatch.setattr(bv, "ohcth", 0.65)
    monkeypatch.setattr(bv, "tfthko", 1.4)
    monkeypatch.setattr(bv, "tfcth", 1.4)
    monkeypatch.setattr(bv, "r_tf_outboard_mid", 1.66e1)
    monkeypatch.setattr(bv, "bore", 2.15)
    monkeypatch.setattr(eh, "idiags", np.full(8, -999999))
    monkeypatch.setattr(fwbsv, "denstl", 7.8e3)
    monkeypatch.setattr(pfv, "rpf1", 0.0)
    monkeypatch.setattr(pfv, "whtpfs", 0.0)
    monkeypatch.setattr(pfv, "curpff", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "nohc", 0)
    monkeypatch.setattr(pfv, "pfrmax", 0.0)
    monkeypatch.setattr(pfv, "fcohbop", 1.0)
    monkeypatch.setattr(pfv, "rjconpf", np.full(22, 1.1e7))
    monkeypatch.setattr(pfv, "ngrp", 4)
    monkeypatch.setattr(pfv, "rohc", 0.0)
    monkeypatch.setattr(pfv, "ncls", np.array([1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
    monkeypatch.setattr(pfv, "zpf", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "cptdin", np.full(22, 4.22e4))
    monkeypatch.setattr(pfv, "pfcaseth", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "itr_sum", 0.0)
    monkeypatch.setattr(pfv, "sigpfcf", 6.66e-1)
    monkeypatch.setattr(pfv, "ohhghf", 9.0e-1)
    monkeypatch.setattr(pfv, "ipfloc", np.array([2, 2, 3, 3, 0, 0, 0, 0, 0, 0]))
    monkeypatch.setattr(pfv, "wts", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "powpfres", 0.0)
    monkeypatch.setattr(pfv, "curpfb", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "routr", 1.5)
    monkeypatch.setattr(pfv, "ric", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "fcohbof", 2.654e-1)
    monkeypatch.setattr(pfv, "rpf2", -1.825)
    monkeypatch.setattr(pfv, "nfxfh", 7)
    monkeypatch.setattr(pfv, "bpf", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "zl", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "wtc", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "vf", np.full(22, 3.0e-1))
    monkeypatch.setattr(pfv, "turns", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "curpfs", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "rpf", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "zref", [3.6, 1.2, 2.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    monkeypatch.setattr(pfv, "pfmmax", 0.0)
    monkeypatch.setattr(pfv, "ipfres", 0)
    monkeypatch.setattr(pfv, "alfapf", 5.0e-10)
    monkeypatch.setattr(pfv, "ncirt", 8)
    monkeypatch.setattr(pfv, "pfclres", 2.5e-8)
    monkeypatch.setattr(pfv, "cpt", np.full([22, 6], 0.0))
    monkeypatch.setattr(pfv, "waves", np.full([22, 6], 0.0))
    monkeypatch.setattr(pfv, "sxlg", np.full([22, 22], 0.0))
    monkeypatch.setattr(pfv, "sigpfcalw", 5.0e2)
    monkeypatch.setattr(pfv, "coheof", 1.6932e7)
    monkeypatch.setattr(pfv, "zh", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "fcohbof", 2.654e-1)
    monkeypatch.setattr(pfv, "ra", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "rb", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "isumatpf", 3)
    monkeypatch.setattr(pfv, "isumatoh", 1)
    monkeypatch.setattr(pfv, "whtpf", 0.0)
    monkeypatch.setattr(pfv, "fcupfsu", 6.900e-1)
    monkeypatch.setattr(pfv, "cohbop", 1.693e7)
    monkeypatch.setattr(pfv, "rjpfalw", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "i_sup_pf_shape", 0)
    monkeypatch.setattr(pfv, "rref", np.full(10, 7.0))
    monkeypatch.setattr(pfv, "i_pf_current", 1)
    monkeypatch.setattr(pfv, "ccl0_ma", np.full(10, 0.0))
    monkeypatch.setattr(pfv, "ccls_ma", np.full(10, 0.0))
    monkeypatch.setattr(pv, "bvert", -6.51e-1)
    monkeypatch.setattr(pv, "kappa", 1.727)
    monkeypatch.setattr(pv, "rli", 1.693)
    monkeypatch.setattr(pv, "itartpf", 0)
    monkeypatch.setattr(pv, "vsres", 6.151e1)
    monkeypatch.setattr(pv, "plascur", 1.8254e7)
    monkeypatch.setattr(pv, "triang", 0.413)
    monkeypatch.setattr(pv, "rminor", 2.883)
    monkeypatch.setattr(pv, "rmajor", 8.938)
    monkeypatch.setattr(pv, "vsind", 3.497e2)
    monkeypatch.setattr(pv, "aspect", 3.1)
    monkeypatch.setattr(pv, "itart", 0)
    monkeypatch.setattr(pv, "betap", 6.313e-1)
    monkeypatch.setattr(tfv, "tftmp", 4.750)
    monkeypatch.setattr(tfv, "dcond", np.full(9, 9.0e3))
    monkeypatch.setattr(tfv, "i_tf_sup", 1)
    monkeypatch.setattr(tfv, "fhts", 0.5)
    monkeypatch.setattr(tfv, "tcritsc", 1.6e1)
    monkeypatch.setattr(tfv, "str_pf_con_res", -5.0e-3)
    monkeypatch.setattr(tfv, "bcritsc", 2.4e1)
    monkeypatch.setattr(tfv, "b_crit_upper_nbti", 1.486e1)
    monkeypatch.setattr(tfv, "t_crit_nbti", 9.04)
    monkeypatch.setattr(tv, "tim", np.full(6, 0.0))
    monkeypatch.setattr(tv, "tramp", 5.0e2)
    monkeypatch.setattr(tv, "tburn", 7.1263e-1)
    monkeypatch.setattr(tv, "tohs", 1.82538e2)
    monkeypatch.setattr(tv, "tqnch", 1.82538e2)
    monkeypatch.setattr(tv, "t_fusion_ramp", 1.0e1)
    monkeypatch.setattr(constants, "dcopper", 8.9e3)
    monkeypatch.setattr(pf, "first_call", True)

    pfcoil.pfcoil()

    assert pytest.approx(pv.bvert) == -0.65121393
    assert pytest.approx(pfv.zpf) == np.array(
        [
            4.86,
            -4.86,
            7.2075,
            -7.2075,
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
        ]
    )


def test_ohcalc(monkeypatch, reinitialise_error_module, pfcoil):
    """Test ohcalc subroutine.

    :param monkeypatch: mocking fixture
    :type monkeypatch: MonkeyPatch
    :param reinitialise_error_module: teardown any error side-effects
    :type reinitialise_error_module: None
    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    """
    # Mocks for ohcalc()
    monkeypatch.setattr(bv, "hmax", 8.864)
    monkeypatch.setattr(bv, "ohcth", 6.510e-1)
    monkeypatch.setattr(fwbsv, "denstl", 7.8e3)
    monkeypatch.setattr(eh, "idiags", np.full(8, 0))
    monkeypatch.setattr(pfv, "nohc", 5)
    monkeypatch.setattr(pfv, "bmaxoh", 1.4e1)
    monkeypatch.setattr(pfv, "i_cs_stress", 0)
    monkeypatch.setattr(pfv, "coheof", 1.693e7)
    monkeypatch.setattr(pfv, "rohc", 0.0)
    monkeypatch.setattr(pfv, "vfohc", 3.0e-1)
    monkeypatch.setattr(pfv, "jstrandoh_bop", 1.069e8)
    monkeypatch.setattr(pfv, "fcuohsu", 7.000e-1)
    monkeypatch.setattr(pfv, "isumatoh", 5)
    monkeypatch.setattr(pfv, "ohhghf", 0.9)
    monkeypatch.setattr(pfv, "areaoh", 1.039e1)
    monkeypatch.setattr(pfv, "powpfres", 0.0)
    monkeypatch.setattr(pfv, "jstrandoh_eof", 1.427e8)
    monkeypatch.setattr(pfv, "powohres", 0.0)
    monkeypatch.setattr(pfv, "rjohc0", 3.048e7)
    monkeypatch.setattr(pfv, "s_tresca_oh", 5.718e8)
    monkeypatch.setattr(pfv, "awpoh", 4.232)
    monkeypatch.setattr(pfv, "oh_steel_frac", 5.926e-1)
    monkeypatch.setattr(pfv, "bmaxoh0", 1.4e1)
    monkeypatch.setattr(pfv, "rjohc", 4.070e7)
    monkeypatch.setattr(pfv, "tmargoh", 1.5)
    monkeypatch.setattr(pfv, "ipfres", 0)
    monkeypatch.setattr(pfv, "rjpfalw", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "pfclres", 2.8e-8)
    monkeypatch.setattr(pfv, "vf", np.full(22, 0.3))
    monkeypatch.setattr(pfv, "ric", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "bpf", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "jscoh_eof", 4.758e8)
    monkeypatch.setattr(pfv, "zpf", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "rb", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "ra", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "jscoh_bop", 3.562e8)
    monkeypatch.setattr(pfv, "cptdin", np.full(22, 4.22e4))
    monkeypatch.setattr(pfv, "pfcaseth", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "rpf", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "cohbop", 1.693e7)
    monkeypatch.setattr(pfv, "zh", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "wtc", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "zl", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "turns", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "wts", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "a_oh_turn", 0.0)
    monkeypatch.setattr(tfv, "dcond", np.full(9, 9.0e3))
    monkeypatch.setattr(tfv, "tftmp", 4.750)
    monkeypatch.setattr(tfv, "tcritsc", 1.6e1)
    monkeypatch.setattr(tfv, "str_cs_con_res", -5.000e-3)
    monkeypatch.setattr(tfv, "fhts", 0.5)
    monkeypatch.setattr(tfv, "bcritsc", 2.4e1)
    monkeypatch.setattr(tfv, "b_crit_upper_nbti", 1.486e1)
    monkeypatch.setattr(tfv, "t_crit_nbti", 9.04)
    monkeypatch.setattr(constants, "dcopper", 8.9e3)

    # Mocks for peakb()
    monkeypatch.setattr(bv, "iohcl", 1)
    monkeypatch.setattr(pfv, "waves", np.full([22, 6], 0.0))
    monkeypatch.setattr(pfv, "ngrp", 4)
    monkeypatch.setattr(pfv, "ncls", np.array([1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]))
    monkeypatch.setattr(pfv, "curpfb", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "curpff", np.full(22, 0.0))
    monkeypatch.setattr(pfv, "curpfs", np.full(22, -175.84911993600002))
    monkeypatch.setattr(pv, "rmajor", 8.938)
    monkeypatch.setattr(pv, "plascur", 1.8254e7)

    # Mocks for hoop_stress()
    monkeypatch.setattr(tfv, "poisson_steel", 3.0e-1)

    # Mocks for superconpf()
    monkeypatch.setattr(eh, "fdiags", np.full(8, -9.99999e5))
    monkeypatch.setattr(tfv, "tmargmin_cs", 1.5)
    monkeypatch.setattr(tfv, "temp_margin", 0.0)
    monkeypatch.setattr(tfv, "b_crit_upper_nbti", 1.486e1)
    monkeypatch.setattr(tfv, "b_crit_upper_nbti", 9.04)

    pfcoil.ohcalc()

    assert pytest.approx(pfv.bpf[4]) == 9.299805e2
    assert pytest.approx(pfv.rjohc) == -7.728453e9


def test_efc(pfcoil, monkeypatch):
    """Test efc subroutine.

    efc() requires specific arguments in order to work; these were discovered
    using gdb to break on the first call of efc() when running the baseline 2019
    IN.DAT.

    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    :param monkeypatch: mocking fixture
    :type monkeypatch: MonkeyPatch
    """
    ngrpmx = 10
    nclsmx = 2
    nptsmx = 32
    nfixmx = 64
    lrow1 = 2 * nptsmx + ngrpmx
    lcol1 = ngrpmx
    npts = 32
    rpts = np.array(
        [
            6.0547741935483881,
            6.2407887617065567,
            6.4268033298647254,
            6.612817898022894,
            6.7988324661810626,
            6.9848470343392313,
            7.1708616024973999,
            7.3568761706555676,
            7.5428907388137372,
            7.7289053069719049,
            7.9149198751300744,
            8.1009344432882422,
            8.2869490114464099,
            8.4729635796045795,
            8.658978147762749,
            8.8449927159209167,
            9.0310072840790845,
            9.217021852237254,
            9.4030364203954235,
            9.5890509885535913,
            9.775065556711759,
            9.9610801248699286,
            10.147094693028098,
            10.333109261186266,
            10.519123829344434,
            10.705138397502601,
            10.891152965660771,
            11.07716753381894,
            11.263182101977108,
            11.449196670135276,
            11.635211238293445,
            11.821225806451615,
        ]
    )
    zpts = np.full(nptsmx, 0.0)
    brin = np.full(nptsmx, 0.0)
    bzin = np.full(nptsmx, 0.0)
    nfix = 14
    rfix = np.full(nfixmx, 0.0)
    rfix[0:14] = 2.3936999999999999
    zfix = np.full(nfixmx, 0.0)
    zfix[0:14] = [
        0.56988544739721259,
        1.7096563421916378,
        2.8494272369860631,
        3.9891981317804879,
        5.1289690265749135,
        6.2687399213693382,
        7.4085108161637638,
        -0.56988544739721259,
        -1.7096563421916378,
        -2.8494272369860631,
        -3.9891981317804879,
        -5.1289690265749135,
        -6.2687399213693382,
        -7.4085108161637638,
    ]
    cfix = np.full(nfixmx, 0.0)
    cfix[0:14] = 12547065.315963898
    ngrp = 4
    ncls = np.array([1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0])

    # This 2D array argument discovered via gdb prints as a 1D array, therefore
    # needs to be reshaped into its original 2D. Fortran ordering is essential
    # when passing greater-than-1D arrays from Python to Fortran
    rcls = np.reshape(
        [
            6.7651653417201345,
            6.7651653417201345,
            18.597693381136555,
            17.000392357531304,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            18.597693381136555,
            17.000392357531304,
            0,
            0,
            0,
            0,
            0,
            0,
        ],
        (10, 2),
        order="F",
    )
    zcls = np.reshape(
        [
            9.8904697261474404,
            -11.124884737289973,
            2.883225806451613,
            8.0730322580645151,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            -2.883225806451613,
            -8.0730322580645151,
            0,
            0,
            0,
            0,
            0,
            0,
        ],
        (10, 2),
        order="F",
    )
    alfa = 5.0e-10
    bfix = np.full(lrow1, 0.0)
    gmat = np.full([lrow1, lcol1], 0.0, order="F")
    bvec = np.full(lrow1, 0.0)
    rc = np.full(nclsmx, 0.0)
    zc = np.full(nclsmx, 0.0)
    cc = np.full(nclsmx, 0.0)
    xc = np.full(nclsmx, 0.0)
    umat = np.full([lrow1, lcol1], 0.0, order="F")
    vmat = np.full([lrow1, lcol1], 0.0, order="F")
    sigma = np.full(ngrpmx, 0.0)
    work2 = np.full(ngrpmx, 0.0)

    ssq, ccls = pfcoil.efc(
        npts,
        rpts,
        zpts,
        brin,
        bzin,
        nfix,
        rfix,
        zfix,
        cfix,
        ngrp,
        ncls,
        rcls,
        zcls,
        alfa,
        bfix,
        gmat,
        bvec,
        rc,
        zc,
        cc,
        xc,
        umat,
        vmat,
        sigma,
        work2,
    )

    assert pytest.approx(ssq) == 4.208729e-4
    assert pytest.approx(ccls[0:4]) == np.array(
        [
            12846165.42893886,
            16377261.02000236,
            579111.6216917,
            20660782.82356247,
        ]
    )


def test_mtrx(pfcoil):
    """Test mtrx subroutine.

    mtrx() requires specific arguments in order to work; these were discovered
    using gdb to break on the first call of mtrx() when running the baseline 2019
    IN.DAT.

    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    """
    lrow1 = 74
    lcol1 = 10
    nptsmx = 32
    npts = 32
    rpts = np.array(
        [
            6.0547741935483881,
            6.2407887617065567,
            6.4268033298647254,
            6.612817898022894,
            6.7988324661810626,
            6.9848470343392313,
            7.1708616024973999,
            7.3568761706555676,
            7.5428907388137372,
            7.7289053069719049,
            7.9149198751300744,
            8.1009344432882422,
            8.2869490114464099,
            8.4729635796045795,
            8.658978147762749,
            8.8449927159209167,
            9.0310072840790845,
            9.217021852237254,
            9.4030364203954235,
            9.5890509885535913,
            9.775065556711759,
            9.9610801248699286,
            10.147094693028098,
            10.333109261186266,
            10.519123829344434,
            10.705138397502601,
            10.891152965660771,
            11.07716753381894,
            11.263182101977108,
            11.449196670135276,
            11.635211238293445,
            11.821225806451615,
        ]
    )
    zpts = np.zeros(nptsmx)
    brin = np.zeros(nptsmx)
    bzin = np.zeros(nptsmx)
    ngrp = 4
    ncls = np.array([1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0])
    rcls = np.reshape(
        [
            0,
            0,
            0,
            18.597693381136555,
            17.000392357531304,
            0,
            0,
            0,
            0,
            0,
            0,
            6.9169190417774516e-323,
            4.9406564584124654e-324,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ],
        (10, 2),
        order="F",
    )
    zcls = np.reshape(
        [
            0,
            0,
            0,
            -2.883225806451613,
            -8.0730322580645151,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ],
        (10, 2),
        order="F",
    )
    alfa = 5.0e-10
    bfix = np.array(
        [
            -4.163336342344337e-17,
            0,
            1.3877787807814457e-17,
            -3.4694469519536142e-17,
            -6.9388939039072284e-18,
            2.0816681711721685e-17,
            -1.3877787807814457e-17,
            3.1225022567582528e-17,
            1.3877787807814457e-17,
            2.7755575615628914e-17,
            3.1225022567582528e-17,
            1.0408340855860843e-17,
            -4.163336342344337e-17,
            -3.4694469519536142e-18,
            -1.7347234759768071e-17,
            1.3877787807814457e-17,
            -2.0816681711721685e-17,
            1.7347234759768071e-17,
            -2.4286128663675299e-17,
            0,
            -3.4694469519536142e-18,
            -1.0408340855860843e-17,
            0,
            0,
            1.7347234759768071e-18,
            -1.0408340855860843e-17,
            1.7347234759768071e-18,
            -6.9388939039072284e-18,
            6.9388939039072284e-18,
            -6.9388939039072284e-18,
            3.4694469519536142e-18,
            5.2041704279304213e-18,
            -0.3130693525427572,
            -0.30317412503141067,
            -0.29349361903088056,
            -0.28403698539500122,
            -0.27481156279537977,
            -0.26582301409269415,
            -0.25707546259584763,
            -0.24857162691582502,
            -0.24031295298976921,
            -0.23229974199262643,
            -0.22453127307719892,
            -0.21700592011907616,
            -0.20972126186523041,
            -0.20267418508484636,
            -0.19586098049525957,
            -0.18927743138451436,
            -0.18291889497602498,
            -0.17678037668189961,
            -0.17085659747167894,
            -0.16514205464473081,
            -0.15963107633981927,
            -0.15431787014642362,
            -0.14919656620141122,
            -0.1442612551639135,
            -0.13950602146198512,
            -0.13492497219911381,
            -0.13051226209776626,
            -0.12626211484245656,
            -0.12216884116737772,
            -0.11822685401402291,
            -0.11443068106371825,
            -0.11077497492862906,
            1.244050972187992e-316,
            9.655957481515668e-97,
            0,
            6.9533558071276999e-310,
            6.9533558069559627e-310,
            6.9533474562307984e-310,
            0,
            0,
            1.2440161899665248e-316,
            9.655957481515668e-97,
        ]
    )

    nrws, gmat, bvec = mtrx(
        lrow1,
        lcol1,
        npts,
        rpts,
        zpts,
        brin,
        bzin,
        ngrp,
        ncls,
        rcls,
        zcls,
        alfa,
        bfix,
        int(pfv.nclsmx),
    )

    gmat_exp = np.array(
        [
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                2.92158969e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                3.04960151e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                3.18186728e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                3.31865171e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                3.46023678e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                3.60692320e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                3.75903204e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                3.91690643e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                4.08091348e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                4.25144640e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                4.42892676e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                4.61380701e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                4.80657328e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                5.00774837e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                5.21789513e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                5.43762013e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                5.66757773e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                5.90847452e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                6.16107431e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                6.42620349e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                6.70475709e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                6.99770541e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                7.30610130e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                7.63108823e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                7.97390922e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                8.33591662e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                8.71858292e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                9.12351264e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                9.55245534e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                1.00073199e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                1.04901901e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                1.10033415e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -3.30317878e-19,
                -3.30317878e-19,
                -6.60635757e-19,
                3.50145552e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -3.20472338e-19,
                -3.20472338e-19,
                -6.40944676e-19,
                3.51776252e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -3.11196728e-19,
                -3.11196728e-19,
                -6.22393456e-19,
                3.53472730e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -3.02442952e-19,
                -3.02442952e-19,
                -6.04885904e-19,
                3.55236630e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.94168179e-19,
                -2.94168179e-19,
                -5.88336358e-19,
                3.57069672e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.86334140e-19,
                -2.86334140e-19,
                -5.72668279e-19,
                3.58973648e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.78906536e-19,
                -2.78906536e-19,
                -5.57813071e-19,
                3.60950425e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.71854537e-19,
                -2.71854537e-19,
                -5.43709074e-19,
                3.63001946e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.65150356e-19,
                -2.65150356e-19,
                -5.30300712e-19,
                3.65130230e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.58768880e-19,
                -2.58768880e-19,
                -5.17537759e-19,
                3.67337373e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.52687355e-19,
                -2.52687355e-19,
                -5.05374710e-19,
                3.69625546e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.46885119e-19,
                -2.46885119e-19,
                -4.93770239e-19,
                3.71996994e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.41343366e-19,
                -2.41343366e-19,
                -4.82686732e-19,
                3.74454037e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.36044938e-19,
                -2.36044938e-19,
                -4.72089877e-19,
                3.76999062e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.30974156e-19,
                -2.30974156e-19,
                -4.61948311e-19,
                3.79634522e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.26116655e-19,
                -2.26116655e-19,
                -4.52233310e-19,
                3.82362931e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.21459257e-19,
                -2.21459257e-19,
                -4.42918515e-19,
                3.85186853e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.16989848e-19,
                -2.16989848e-19,
                -4.33979695e-19,
                3.88108893e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.12697269e-19,
                -2.12697269e-19,
                -4.25394538e-19,
                3.91131688e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.08571231e-19,
                -2.08571231e-19,
                -4.17142461e-19,
                3.94257887e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.04602225e-19,
                -2.04602225e-19,
                -4.09204451e-19,
                3.97490133e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -2.00781456e-19,
                -2.00781456e-19,
                -4.01562911e-19,
                4.00831037e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.97100769e-19,
                -1.97100769e-19,
                -3.94201538e-19,
                4.04283150e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.93552600e-19,
                -1.93552600e-19,
                -3.87105201e-19,
                4.07848927e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.90129919e-19,
                -1.90129919e-19,
                -3.80259839e-19,
                4.11530678e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.86826185e-19,
                -1.86826185e-19,
                -3.73652370e-19,
                4.15330516e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.83635302e-19,
                -1.83635302e-19,
                -3.67270604e-19,
                4.19250290e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.80551586e-19,
                -1.80551586e-19,
                -3.61103172e-19,
                4.23291505e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.77569727e-19,
                -1.77569727e-19,
                -3.55139453e-19,
                4.27455221e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.74684759e-19,
                -1.74684759e-19,
                -3.49369519e-19,
                4.31741941e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.71892037e-19,
                -1.71892037e-19,
                -3.43784075e-19,
                4.36151462e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                -1.69187206e-19,
                -1.69187206e-19,
                -3.38374412e-19,
                4.40682706e-08,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                5.00000000e-10,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                5.00000000e-10,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                1.00000000e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                1.00000000e-09,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
        ]
    )

    bvec_exp = np.array(
        [
            4.16333634e-17,
            0.00000000e00,
            -1.38777878e-17,
            3.46944695e-17,
            6.93889390e-18,
            -2.08166817e-17,
            1.38777878e-17,
            -3.12250226e-17,
            -1.38777878e-17,
            -2.77555756e-17,
            -3.12250226e-17,
            -1.04083409e-17,
            4.16333634e-17,
            3.46944695e-18,
            1.73472348e-17,
            -1.38777878e-17,
            2.08166817e-17,
            -1.73472348e-17,
            2.42861287e-17,
            0.00000000e00,
            3.46944695e-18,
            1.04083409e-17,
            0.00000000e00,
            0.00000000e00,
            -1.73472348e-18,
            1.04083409e-17,
            -1.73472348e-18,
            6.93889390e-18,
            -6.93889390e-18,
            6.93889390e-18,
            -3.46944695e-18,
            -5.20417043e-18,
            3.13069353e-01,
            3.03174125e-01,
            2.93493619e-01,
            2.84036985e-01,
            2.74811563e-01,
            2.65823014e-01,
            2.57075463e-01,
            2.48571627e-01,
            2.40312953e-01,
            2.32299742e-01,
            2.24531273e-01,
            2.17005920e-01,
            2.09721262e-01,
            2.02674185e-01,
            1.95860980e-01,
            1.89277431e-01,
            1.82918895e-01,
            1.76780377e-01,
            1.70856597e-01,
            1.65142055e-01,
            1.59631076e-01,
            1.54317870e-01,
            1.49196566e-01,
            1.44261255e-01,
            1.39506021e-01,
            1.34924972e-01,
            1.30512262e-01,
            1.26262115e-01,
            1.22168841e-01,
            1.18226854e-01,
            1.14430681e-01,
            1.10774975e-01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
        ]
    )

    assert nrws == 68
    assert_array_almost_equal(gmat, gmat_exp)
    assert_array_almost_equal(bvec, bvec_exp)


def test_solv(pfcoil):
    """Test solv() with simple arguments.

    Running baseline_2019 results in 2D array args with 740 elements: unfeasible
    for a unit test. Made-up simplified args are used instead.

    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    """
    ngrpmx = 3
    ngrp = 3
    nrws = 3
    gmat = np.full((3, 3), 2.0, order="F")
    bvec = np.full(3, 1.0)

    ccls, umat, vmat, sigma, work2 = pfcoil.solv(ngrpmx, ngrp, nrws, gmat, bvec)

    assert_array_almost_equal(ccls, np.array([0.16666667, 0.37079081, -0.03745748]))
    assert_array_almost_equal(
        umat,
        np.array(
            [
                [-0.81649658, -0.57735027, 0.0],
                [0.40824829, -0.57735027, -0.70710678],
                [0.40824829, -0.57735027, 0.70710678],
            ]
        ),
    )
    assert_array_almost_equal(
        vmat,
        np.array(
            [
                [-0.81649658, -0.57735027, 0.0],
                [0.40824829, -0.57735027, -0.70710678],
                [0.40824829, -0.57735027, 0.70710678],
            ]
        ),
    )
    assert_array_almost_equal(
        sigma, np.array([5.1279005e-16, 6.0000000e00, 0.0000000e00])
    )
    assert_array_almost_equal(
        work2, np.array([-2.22044605e-16, -1.73205081e00, 0.00000000e00])
    )


def test_fixb(pfcoil):
    """Test fixb subroutine.

    fixb() requires specific arguments in order to work; these were discovered
    using gdb to break on the first subroutine call when running the baseline
    2018 IN.DAT.

    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    """
    nptsmx = 32
    lrow1 = 74
    npts = 32
    rpts = np.array(
        [
            6.0223258064516134,
            6.2073434963579608,
            6.3923611862643082,
            6.5773788761706555,
            6.7623965660770038,
            6.9474142559833512,
            7.1324319458896985,
            7.3174496357960459,
            7.5024673257023942,
            7.6874850156087415,
            7.8725027055150889,
            8.0575203954214363,
            8.2425380853277836,
            8.427555775234131,
            8.6125734651404784,
            8.7975911550468275,
            8.9826088449531731,
            9.1676265348595223,
            9.3526442247658697,
            9.537661914672217,
            9.7226796045785644,
            9.9076972944849118,
            10.092714984391259,
            10.277732674297607,
            10.462750364203954,
            10.647768054110301,
            10.83278574401665,
            11.017803433922996,
            11.202821123829345,
            11.387838813735691,
            11.57285650364204,
            11.757874193548387,
        ]
    )
    zpts = np.zeros(nptsmx)
    nfix = 14
    rfix = np.array(
        [
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            2.6084100000000001,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
    )
    zfix = np.array(
        [
            0.58327007281470211,
            1.7498102184441064,
            2.9163503640735104,
            4.0828905097029144,
            5.2494306553323193,
            6.4159708009617233,
            7.5825109465911273,
            -0.58327007281470211,
            -1.7498102184441064,
            -2.9163503640735104,
            -4.0828905097029144,
            -5.2494306553323193,
            -6.4159708009617233,
            -7.5825109465911273,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
    )
    cfix = np.array(
        [
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            12444820.564847374,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ]
    )

    bfix_exp = np.array(
        [
            -2.77555756e-17,
            -2.77555756e-17,
            -2.08166817e-17,
            -6.24500451e-17,
            4.85722573e-17,
            8.32667268e-17,
            6.93889390e-18,
            4.16333634e-17,
            2.08166817e-17,
            1.38777878e-17,
            -2.77555756e-17,
            -2.77555756e-17,
            1.38777878e-17,
            3.12250226e-17,
            -1.38777878e-17,
            -3.46944695e-18,
            6.93889390e-18,
            -2.77555756e-17,
            3.46944695e-18,
            2.42861287e-17,
            4.51028104e-17,
            1.38777878e-17,
            1.04083409e-17,
            -6.93889390e-18,
            1.38777878e-17,
            1.73472348e-17,
            2.77555756e-17,
            -6.93889390e-18,
            -1.73472348e-18,
            2.25514052e-17,
            1.04083409e-17,
            1.73472348e-18,
            -3.53728301e-01,
            -3.43046326e-01,
            -3.32568315e-01,
            -3.22305794e-01,
            -3.12268396e-01,
            -3.02463988e-01,
            -2.92898791e-01,
            -2.83577512e-01,
            -2.74503461e-01,
            -2.65678678e-01,
            -2.57104045e-01,
            -2.48779412e-01,
            -2.40703694e-01,
            -2.32874990e-01,
            -2.25290668e-01,
            -2.17947471e-01,
            -2.10841592e-01,
            -2.03968763e-01,
            -1.97324322e-01,
            -1.90903284e-01,
            -1.84700401e-01,
            -1.78710219e-01,
            -1.72927125e-01,
            -1.67345392e-01,
            -1.61959222e-01,
            -1.56762779e-01,
            -1.51750219e-01,
            -1.46915716e-01,
            -1.42253492e-01,
            -1.37757829e-01,
            -1.33423089e-01,
            -1.29243733e-01,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
        ]
    )

    bfix = fixb(lrow1, npts, rpts, zpts, nfix, rfix, zfix, cfix)

    assert_array_almost_equal(bfix, bfix_exp)


def test_peakb(monkeypatch, pfcoil):
    """Test peakb subroutine.

    peakb() requires specific arguments in order to work; these were discovered
    using gdb to break on the first subroutine call when running the baseline
    2018 IN.DAT.
    :param monkeypatch: mocking fixture
    :type monkeypatch: MonkeyPatch
    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    """
    # Mock module vars
    monkeypatch.setattr(pf, "nfxf", 14)
    monkeypatch.setattr(
        pf,
        "rfxf",
        np.array(
            (
                6.2732560483870969,
                6.2732560483870969,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                2.6084100000000001,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            )
        ),
    )
    monkeypatch.setattr(
        pf,
        "zfxf",
        np.array(
            (
                9.606146709677418,
                -11.141090021562032,
                2.9163503640735104,
                4.0828905097029144,
                5.2494306553323193,
                6.4159708009617233,
                7.5825109465911273,
                -0.58327007281470211,
                -1.7498102184441064,
                -2.9163503640735104,
                -4.0828905097029144,
                -5.2494306553323193,
                -6.4159708009617233,
                -7.5825109465911273,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            )
        ),
    )
    monkeypatch.setattr(
        pf,
        "cfxf",
        np.array(
            (
                15889161.548344759,
                18583102.707854092,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                12444820.564847374,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            )
        ),
    )
    monkeypatch.setattr(pf, "xind", np.zeros(64))
    monkeypatch.setattr(pf, "bpf2", np.zeros(22))

    monkeypatch.setattr(bv, "iohcl", 1)
    monkeypatch.setattr(bv, "hmax", 9.0730900215620327)
    monkeypatch.setattr(bv, "ohcth", 0.55242000000000002)
    monkeypatch.setattr(
        eh,
        "idiags",
        np.array(
            [-999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "ra",
        np.array(
            [
                5.6944236847973242,
                5.5985055619292972,
                17.819978201682968,
                17.819978201682968,
                16.385123084628962,
                16.385123084628962,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(pfv, "nohc", 7)
    monkeypatch.setattr(
        pfv,
        "ric",
        np.array(
            [
                14.742063826112622,
                20.032681634901664,
                -8.1098913365453491,
                -8.1098913365453491,
                -5.5984385047179153,
                -5.5984385047179153,
                -186.98751599968148,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(pfv, "ohhghf", 0.90000000000000002)
    monkeypatch.setattr(
        pfv,
        "rb",
        np.array(
            [
                6.8520884119768697,
                6.9480065348448967,
                18.98258241468535,
                18.98258241468535,
                17.22166645654087,
                17.22166645654087,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "rpf",
        np.array(
            [
                6.2732560483870969,
                6.2732560483870969,
                18.401280308184159,
                18.401280308184159,
                16.803394770584916,
                16.803394770584916,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "waves",
        np.array(
            [
                [0.0, 1.0, 0.00457346, 0.00457346, 0.00457346, 0.0],
                [0.0, 1.0, -0.14559845, -0.14559845, -0.14559845, 0.0],
                [0.0, -0.07156774, 1.0, 1.0, 1.0, 0.0],
                [0.0, -0.07156774, 1.0, 1.0, 1.0, 0.0],
                [0.0, -0.07676189, 1.0, 1.0, 1.0, 0.0],
                [0.0, -0.07676189, 1.0, 1.0, 1.0, 0.0],
                [0.0, -0.93176, 1.0, 1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            ],
            order="F",
        ),
    )
    monkeypatch.setattr(pfv, "coheof", 20726000)
    monkeypatch.setattr(pfv, "ngrp", 4)
    monkeypatch.setattr(pfv, "bpf", np.zeros(22, dtype=int))  # maybe
    monkeypatch.setattr(
        pfv,
        "zpf",
        np.array(
            [
                9.606146709677418,
                -11.141090021562032,
                2.8677741935483869,
                -2.8677741935483869,
                8.0297677419354834,
                -8.0297677419354834,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(pfv, "ncls", np.array([1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0]))
    monkeypatch.setattr(
        pfv,
        "zl",
        np.array(
            [
                9.0273143460876444,
                -10.466339535104233,
                2.2864720870471951,
                -2.2864720870471951,
                7.6114960559795311,
                -7.6114960559795311,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "curpfb",
        np.array(
            [
                0.067422231232391661,
                -2.9167273287450968,
                -8.1098913365453491,
                -8.1098913365453491,
                -5.5984385047179153,
                -5.5984385047179153,
                -186.98751599968148,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "curpff",
        np.array(
            [
                0.067422231232391661,
                -2.9167273287450968,
                -8.1098913365453491,
                -8.1098913365453491,
                -5.5984385047179153,
                -5.5984385047179153,
                -186.98751599968148,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "curpfs",
        np.array(
            [
                14.742063826112622,
                20.032681634901664,
                0.58040662653667285,
                0.58040662653667285,
                0.42974674788703021,
                0.42974674788703021,
                174.22748790786324,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(pfv, "cohbop", 19311657.760000002)
    monkeypatch.setattr(
        pfv,
        "zh",
        np.array(
            [
                10.184979073267192,
                -11.815840508019832,
                3.4490763000495788,
                -3.4490763000495788,
                8.4480394278914357,
                -8.4480394278914357,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(pv, "rmajor", 8.8901000000000003)
    monkeypatch.setattr(pv, "plascur", 17721306.969367817)

    i = 1
    ii = 1
    it = 32767

    bri_exp = 6.250031e-1
    bro_exp = 4.339524e-1
    bzi_exp = 1.049564e1
    bzo_exp = -6.438987

    bri, bro, bzi, bzo = pfcoil.peakb(i, ii, it)

    assert pytest.approx(bri) == bri_exp
    assert pytest.approx(bro) == bro_exp
    assert pytest.approx(bzi) == bzi_exp
    assert pytest.approx(bzo) == bzo_exp


def test_superconpf(monkeypatch, pfcoil):
    """Test superconpf subroutine.

    superconpf() requires specific arguments in order to work; these were
    discovered using gdb to break on the first subroutine call when running the
    baseline 2018 IN.DAT.
    :param monkeypatch: mocking fixture
    :type monkeypatch: _pytest.monkeypatch.MonkeyPatch
    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    """
    # TODO This test would benefit from parameterisation for different SC
    # materials (isumat)
    monkeypatch.setattr(eh, "fdiags", np.zeros(8))
    monkeypatch.setattr(eh, "idiags", np.zeros(8))
    monkeypatch.setattr(tfv, "tmargmin_cs", 0.0)
    monkeypatch.setattr(tfv, "temp_margin", 0.0)
    monkeypatch.setattr(tfv, "b_crit_upper_nbti", 0.0)
    monkeypatch.setattr(tfv, "t_crit_nbti", 0.0)

    bmax = 10.514241695080285
    fhe = 0.29999999999999999
    fcu = 0.68999999999999995
    jwp = 11000000
    isumat = 3
    fhts = 0.5
    strain = -0.0050000000000000001
    thelium = 4.75
    bcritsc = 24
    tcritsc = 16

    jcritwp_exp = -2.6716372e7
    jcritstr_exp = -3.8166246e7
    jcritsc_exp = -1.23116924e8
    tmarg_exp = -2.651537e-1

    jcritwp, jcritstr, jcritsc, tmarg = pfcoil.superconpf(
        bmax, fhe, fcu, jwp, isumat, fhts, strain, thelium, bcritsc, tcritsc
    )

    assert pytest.approx(jcritwp) == jcritwp_exp
    assert pytest.approx(jcritstr) == jcritstr_exp
    assert pytest.approx(jcritsc) == jcritsc_exp
    assert pytest.approx(tmarg) == tmarg_exp


def test_axial_stress(pfcoil, monkeypatch):
    """Test axial_stress subroutine.

    axial_stress() requires specific mocked vars in order to work; these were
    discovered using gdb to break on the first subroutine call when running the
    baseline 2018 IN.DAT.

    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    :param monkeypatch: mocking fixture
    :type monkeypatch: _pytest.monkeypatch.MonkeyPatch
    """
    monkeypatch.setattr(pfv, "oh_steel_frac", 0.57874999999999999)
    monkeypatch.setattr(pfv, "nohc", 7)
    monkeypatch.setattr(
        pfv,
        "rb",
        np.array(
            [
                6.8520884119768697,
                6.9480065348448967,
                18.98258241468535,
                18.98258241468535,
                17.22166645654087,
                17.22166645654087,
                2.88462,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "ra",
        np.array(
            [
                5.6944236847973242,
                5.5985055619292972,
                17.819978201682968,
                17.819978201682968,
                16.385123084628962,
                16.385123084628962,
                2.3321999999999998,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "ric",
        np.array(
            [
                14.742063826112622,
                20.032681634901664,
                -8.1098913365453491,
                -8.1098913365453491,
                -5.5984385047179153,
                -5.5984385047179153,
                -186.98751599968145,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "zh",
        np.array(
            [
                10.184979073267192,
                -11.815840508019832,
                3.4490763000495788,
                -3.4490763000495788,
                8.4480394278914357,
                -8.4480394278914357,
                8.1657810194058289,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )

    s_axial_exp = -7.468967e8
    axial_force_exp = -1.956801e9
    s_axial, axial_force = pfcoil.axial_stress()

    assert pytest.approx(s_axial) == s_axial_exp
    assert pytest.approx(axial_force) == axial_force_exp


def test_induct(pfcoil, monkeypatch):
    """Test induct subroutine.

    induct() requires specific mocked vars in order to work; these were
    discovered using gdb to break on the first subroutine call when running the
    baseline 2018 IN.DAT.

    :param pfcoil: a PFCoil instance
    :type pfcoil: process.pfcoil.PFCoil
    :param monkeypatch: mocking fixture
    :type monkeypatch: _pytest.monkeypatch.MonkeyPatch
    """
    monkeypatch.setattr(bv, "iohcl", 1)
    monkeypatch.setattr(bv, "ohcth", 0.55242000000000002)
    monkeypatch.setattr(
        eh,
        "fdiags",
        np.array(
            [-999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999]
        ),
    )
    monkeypatch.setattr(
        eh,
        "idiags",
        np.array(
            [-999999, -999999, -999999, -999999, -999999, -999999, -999999, -999999]
        ),
    )
    monkeypatch.setattr(pfv, "nohc", 7)
    monkeypatch.setattr(
        pfv,
        "turns",
        np.array(
            [
                349.33800535811901,
                474.70809561378354,
                192.17751982334951,
                192.17751982334951,
                130.19624429576547,
                130.19624429576547,
                4348.5468837135222,
                1,
                100,
                100,
                100,
                100,
                100,
                100,
                100,
                100,
                100,
                100,
                100,
                100,
                100,
                100,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "zpf",
        np.array(
            [
                9.606146709677418,
                -11.141090021562032,
                2.8677741935483869,
                -2.8677741935483869,
                8.0297677419354834,
                -8.0297677419354834,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "rpf",
        np.array(
            [
                6.2732560483870969,
                6.2732560483870969,
                18.401280308184159,
                18.401280308184159,
                16.803394770584916,
                16.803394770584916,
                2.6084100000000001,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(pfv, "sxlg", np.ones((22, 22), dtype=int))
    monkeypatch.setattr(pfv, "rohc", 2.6084100000000001)
    monkeypatch.setattr(pfv, "ngrp", 4)
    monkeypatch.setattr(pfv, "ncls", np.array([1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0]))
    monkeypatch.setattr(
        pfv,
        "zl",
        np.array(
            [
                9.0273143460876444,
                -10.466339535104233,
                2.2864720870471951,
                -2.2864720870471951,
                7.6114960559795311,
                -7.6114960559795311,
                -8.1657810194058289,
                -5.2996467096774191,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(pfv, "ncirt", 8)
    monkeypatch.setattr(
        pfv,
        "ra",
        np.array(
            [
                5.6944236847973242,
                5.5985055619292972,
                17.819978201682968,
                17.819978201682968,
                16.385123084628962,
                16.385123084628962,
                2.3321999999999998,
                6.0223258064516134,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "zh",
        np.array(
            [
                10.184979073267192,
                -11.815840508019832,
                3.4490763000495788,
                -3.4490763000495788,
                8.4480394278914357,
                -8.4480394278914357,
                8.1657810194058289,
                5.2996467096774191,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(
        pfv,
        "rb",
        np.array(
            [
                6.8520884119768697,
                6.9480065348448967,
                18.98258241468535,
                18.98258241468535,
                17.22166645654087,
                17.22166645654087,
                2.88462,
                11.757874193548387,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
            ]
        ),
    )
    monkeypatch.setattr(pv, "rmajor", 8.8901000000000003)
    monkeypatch.setattr(pv, "rlp", 1.6039223939491056e-05)

    sxlg_exp = np.array(
        [
            [
                2.49332453e00,
                4.46286166e-02,
                2.38094100e-01,
                1.57653632e-01,
                2.18695928e-01,
                6.62002005e-02,
                8.81068392e-01,
                8.15132226e-04,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                4.46286166e-02,
                4.33166888e00,
                1.89207090e-01,
                2.93333228e-01,
                7.84212470e-02,
                2.83752898e-01,
                8.54403195e-01,
                8.60878436e-04,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                2.38094100e-01,
                1.89207090e-01,
                3.12872133e00,
                1.10843611e00,
                7.24769254e-01,
                3.90823361e-01,
                5.46263544e-01,
                1.70440906e-03,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                1.57653632e-01,
                2.93333228e-01,
                1.10843611e00,
                3.12872133e00,
                3.90823361e-01,
                7.24769254e-01,
                5.46263544e-01,
                1.70440906e-03,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                2.18695928e-01,
                7.84212470e-02,
                7.24769254e-01,
                3.90823361e-01,
                1.39661265e00,
                1.50164883e-01,
                3.27696035e-01,
                8.81560519e-04,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                6.62002005e-02,
                2.83752898e-01,
                3.90823361e-01,
                7.24769254e-01,
                1.50164883e-01,
                1.39661265e00,
                3.27696035e-01,
                8.81560519e-04,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                8.81068392e-01,
                8.54403195e-01,
                5.46263544e-01,
                5.46263544e-01,
                3.27696035e-01,
                3.27696035e-01,
                2.50139308e01,
                4.90307127e-03,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                8.15132226e-04,
                8.60878436e-04,
                1.70440906e-03,
                1.70440906e-03,
                8.81560519e-04,
                8.81560519e-04,
                4.90307127e-03,
                1.60392239e-05,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
                0.00000000e00,
            ],
        ]
    )
    pfcoil.induct(False)
    assert_array_almost_equal(pfv.sxlg, sxlg_exp)
