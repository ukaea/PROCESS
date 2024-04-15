import pytest
import numpy
from typing import NamedTuple, Any

from process.fortran import sctfcoil_module
from process.fortran import tfcoil_variables
from process.fortran import global_variables
from process.fortran import physics_variables
from process.fortran import build_variables
from process.fortran import fwbs_variables
from process.fortran import divertor_variables
from process.sctfcoil import Sctfcoil
from process import sctfcoil as sctf


@pytest.fixture
def sctfcoil():
    """Provides Sctfcoil object for testing.

    :returns: initialised Sctfcoil object
    :rtype: process.sctfcoil.Sctfcoil
    """
    return Sctfcoil()


class ProtectParam(NamedTuple):
    aio: Any = None

    tfes: Any = None

    acs: Any = None

    aturn: Any = None

    tdump: Any = None

    fcond: Any = None

    fcu: Any = None

    tba: Any = None

    tmax: Any = None

    expected_ajwpro: Any = None

    expected_vd: Any = None


@pytest.mark.parametrize(
    "protectparam",
    (
        ProtectParam(
            aio=74026.751437500003,
            tfes=9561415368.8360519,
            acs=0.001293323051622732,
            aturn=0.0032012300777680192,
            tdump=25.829000000000001,
            fcond=0.63927285511442711,
            fcu=0.80884,
            tba=4.75,
            tmax=150,
            expected_ajwpro=17475706.393616617,
            expected_vd=10001.287165953383,
        ),
    ),
)
def test_protect(protectparam, sctfcoil):
    """
    Automatically generated Regression Unit Test for protect.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param protectparam: the data used to mock and assert in this test.
    :type protectparam: protectparam

    :param sctfcoil: initialised Sctfcoil object
    :type sctfcoil: process.sctfcoil.Sctfcoil
    """

    ajwpro, vd = sctfcoil.protect(
        aio=protectparam.aio,
        tfes=protectparam.tfes,
        acs=protectparam.acs,
        aturn=protectparam.aturn,
        tdump=protectparam.tdump,
        fcond=protectparam.fcond,
        fcu=protectparam.fcu,
        tba=protectparam.tba,
        tmax=protectparam.tmax,
    )

    assert ajwpro == pytest.approx(protectparam.expected_ajwpro)

    assert vd == pytest.approx(protectparam.expected_vd)


class SuperconParam(NamedTuple):
    tmargmin_tf: Any = None

    n_tf: Any = None

    temp_margin: Any = None

    jwdgpro: Any = None

    dhecoil: Any = None

    cpttf: Any = None

    bmaxtfrp: Any = None

    str_tf_con_res: Any = None

    b_crit_upper_nbti: Any = None

    i_str_wp: Any = None

    str_wp: Any = None

    t_crit_nbti: Any = None

    tf_fit_t: Any = None

    tf_fit_z: Any = None

    tf_fit_y: Any = None

    run_tests: Any = None

    isumat: Any = None

    iprint: Any = None

    outfile: Any = None

    acs: Any = None

    aturn: Any = None

    bmax: Any = None

    fcu: Any = None

    fhe: Any = None

    fhts: Any = None

    iop: Any = None

    jwp: Any = None

    tdmptf: Any = None

    tfes: Any = None

    thelium: Any = None

    tmax: Any = None

    bcritsc: Any = None

    tcritsc: Any = None

    expected_temp_margin: Any = None

    expected_jwdgpro: Any = None

    expected_jwdgcrt: Any = None

    expected_vd: Any = None

    expected_tmarg: Any = None


@pytest.mark.parametrize(
    "superconparam",
    (
        SuperconParam(
            tmargmin_tf=1.5,
            n_tf=16,
            temp_margin=0,
            jwdgpro=0,
            dhecoil=0.010000000000000002,
            cpttf=74026.751437500003,
            bmaxtfrp=12.48976756562082,
            str_tf_con_res=-0.0050000000000000001,
            b_crit_upper_nbti=14.859999999999999,
            i_str_wp=1,
            str_wp=0.0015619754370069119,
            t_crit_nbti=9.0399999999999991,
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            tf_fit_y=1.0658869305062604,
            run_tests=0,
            isumat=5,
            iprint=0,
            outfile=11,
            acs=0.001293323051622732,
            aturn=0.0032012300777680192,
            bmax=12.48976756562082,
            fcu=0.80884,
            fhe=0.30000000000000004,
            fhts=0.5,
            iop=74026.751437500003,
            jwp=23124470.793774806,
            tdmptf=25.829000000000001,
            tfes=9548964780.4287167,
            thelium=4.75,
            tmax=150,
            bcritsc=24,
            tcritsc=16,
            expected_temp_margin=2.3431632224075836,
            expected_jwdgpro=17475706.393616617,
            expected_jwdgcrt=41107234.360397324,
            expected_vd=9988.2637896807955,
            expected_tmarg=2.3431632224075836,
        ),
        SuperconParam(
            tmargmin_tf=1.5,
            n_tf=16,
            temp_margin=2.3431632224075836,
            jwdgpro=17475706.393616617,
            dhecoil=0.010000000000000002,
            cpttf=74026.751437500003,
            bmaxtfrp=12.48976756562082,
            str_tf_con_res=-0.0050000000000000001,
            b_crit_upper_nbti=14.859999999999999,
            i_str_wp=1,
            str_wp=0.0015619754370069119,
            t_crit_nbti=9.0399999999999991,
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            tf_fit_y=1.0658869305062604,
            run_tests=0,
            isumat=5,
            iprint=0,
            outfile=11,
            acs=0.001293323051622732,
            aturn=0.0032012300777680192,
            bmax=12.48976756562082,
            fcu=0.80884,
            fhe=0.30000000000000004,
            fhts=0.5,
            iop=74026.751437500003,
            jwp=23124470.793774806,
            tdmptf=25.829000000000001,
            tfes=9561415368.8360519,
            thelium=4.75,
            tmax=150,
            bcritsc=24,
            tcritsc=16,
            expected_temp_margin=2.3431632224075836,
            expected_jwdgpro=17475706.393616617,
            expected_jwdgcrt=41107234.360397324,
            expected_vd=10001.287165953383,
            expected_tmarg=2.3431632224075836,
        ),
        SuperconParam(
            tmargmin_tf=1.5,
            n_tf=16,
            temp_margin=2.3431632224075836,
            jwdgpro=17475706.393616617,
            dhecoil=0.010000000000000002,
            cpttf=74026.751437500003,
            bmaxtfrp=12.48976756562082,
            str_tf_con_res=-0.0050000000000000001,
            b_crit_upper_nbti=14.859999999999999,
            i_str_wp=1,
            str_wp=0.0015619754370069119,
            t_crit_nbti=9.0399999999999991,
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            tf_fit_y=1.0658869305062604,
            run_tests=0,
            isumat=5,
            iprint=0,
            outfile=11,
            acs=0.001293323051622732,
            aturn=0.0032012300777680192,
            bmax=12.48976756562082,
            fcu=0.80884,
            fhe=0.30000000000000004,
            fhts=0.5,
            iop=74026.751437500003,
            jwp=23124470.793774806,
            tdmptf=25.829000000000001,
            tfes=9561415368.8360519,
            thelium=4.75,
            tmax=150,
            bcritsc=24,
            tcritsc=16,
            expected_temp_margin=2.3431632224075836,
            expected_jwdgpro=17475706.393616617,
            expected_jwdgcrt=41107234.360397324,
            expected_vd=10001.287165953383,
            expected_tmarg=2.3431632224075836,
        ),
    ),
)
def test_supercon(superconparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for supercon.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param superconparam: the data used to mock and assert in this test.
    :type superconparam: superconparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param sctfcoil: initialised Sctfcoil object
    :type sctfcoil: process.sctfcoil.Sctfcoil
    """

    monkeypatch.setattr(tfcoil_variables, "tmargmin_tf", superconparam.tmargmin_tf)

    monkeypatch.setattr(tfcoil_variables, "n_tf", superconparam.n_tf)

    monkeypatch.setattr(tfcoil_variables, "temp_margin", superconparam.temp_margin)

    monkeypatch.setattr(tfcoil_variables, "jwdgpro", superconparam.jwdgpro)

    monkeypatch.setattr(tfcoil_variables, "dhecoil", superconparam.dhecoil)

    monkeypatch.setattr(tfcoil_variables, "cpttf", superconparam.cpttf)

    monkeypatch.setattr(tfcoil_variables, "bmaxtfrp", superconparam.bmaxtfrp)

    monkeypatch.setattr(
        tfcoil_variables, "str_tf_con_res", superconparam.str_tf_con_res
    )

    monkeypatch.setattr(
        tfcoil_variables, "b_crit_upper_nbti", superconparam.b_crit_upper_nbti
    )

    monkeypatch.setattr(tfcoil_variables, "i_str_wp", superconparam.i_str_wp)

    monkeypatch.setattr(tfcoil_variables, "str_wp", superconparam.str_wp)

    monkeypatch.setattr(tfcoil_variables, "t_crit_nbti", superconparam.t_crit_nbti)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_t", superconparam.tf_fit_t)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_z", superconparam.tf_fit_z)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_y", superconparam.tf_fit_y)

    monkeypatch.setattr(global_variables, "run_tests", superconparam.run_tests)

    jwdgcrt, vd, tmarg = sctfcoil.supercon(
        isumat=superconparam.isumat,
        acs=superconparam.acs,
        aturn=superconparam.aturn,
        bmax=superconparam.bmax,
        fcu=superconparam.fcu,
        fhe=superconparam.fhe,
        fhts=superconparam.fhts,
        iop=superconparam.iop,
        jwp=superconparam.jwp,
        tdmptf=superconparam.tdmptf,
        tfes=superconparam.tfes,
        thelium=superconparam.thelium,
        tmax=superconparam.tmax,
        bcritsc=superconparam.bcritsc,
        tcritsc=superconparam.tcritsc,
        output=False,
    )

    assert tfcoil_variables.temp_margin == pytest.approx(
        superconparam.expected_temp_margin
    )

    assert tfcoil_variables.jwdgpro == pytest.approx(superconparam.expected_jwdgpro)

    assert jwdgcrt == pytest.approx(superconparam.expected_jwdgcrt)

    assert vd == pytest.approx(superconparam.expected_vd)

    assert tmarg == pytest.approx(superconparam.expected_tmarg)


class TfCurrentParam(NamedTuple):
    casthi: Any = None

    ritfc: Any = None

    rbmax: Any = None

    i_tf_sup: Any = None

    casths_fraction: Any = None

    tinstf: Any = None

    tftort: Any = None

    bmaxtf: Any = None

    tfinsgap: Any = None

    tfc_sidewall_is_fraction: Any = None

    casths: Any = None

    casthi_is_fraction: Any = None

    casthi_fraction: Any = None

    n_tf: Any = None

    thicndut: Any = None

    thkcas: Any = None

    oacdcp: Any = None

    tfareain: Any = None

    r_tf_inboard_out: Any = None

    r_tf_inboard_in: Any = None

    tfcth: Any = None

    bt: Any = None

    rmajor: Any = None

    tfc_current: Any = None

    theta_coil: Any = None

    expected_ritfc: Any = None

    expected_rbmax: Any = None

    expected_bmaxtf: Any = None

    expected_oacdcp: Any = None

    expected_tfc_current: Any = None


@pytest.mark.parametrize(
    "tfcurrentparam",
    (
        TfCurrentParam(
            casthi=0.060000000000000012,
            ritfc=0,
            rbmax=0,
            i_tf_sup=1,
            casths_fraction=0.059999999999999998,
            tinstf=0.0080000000000000019,
            tftort=1.6395161177915356,
            bmaxtf=0,
            tfinsgap=0.01,
            tfc_sidewall_is_fraction=False,
            casths=0.05000000000000001,
            casthi_is_fraction=False,
            casthi_fraction=0.050000000000000003,
            n_tf=16,
            thicndut=0.002,
            thkcas=0.52465000000000006,
            oacdcp=8673900,
            tfareain=27.308689677971632,
            r_tf_inboard_out=4.20194118510911,
            r_tf_inboard_in=2.9939411851091102,
            tfcth=1.208,
            bt=5.3292000000000002,
            rmajor=8.8901000000000003,
            tfc_current=0,
            theta_coil=0.19634954084936207,
            expected_ritfc=236885604.60000002,
            expected_rbmax=4.0432020634751211,
            expected_bmaxtf=11.717722779177526,
            expected_oacdcp=8674367.2945641987,
            expected_tfc_current=14805350.287500001,
        ),
    ),
)
def test_tf_current(tfcurrentparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_current.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcurrentparam: the data used to mock and assert in this test.
    :type tfcurrentparam: tfcurrentparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param sctfcoil: initialised Sctfcoil object
    :type sctfcoil: process.sctfcoil.Sctfcoil
    """

    monkeypatch.setattr(tfcoil_variables, "casthi", tfcurrentparam.casthi)

    monkeypatch.setattr(tfcoil_variables, "ritfc", tfcurrentparam.ritfc)

    monkeypatch.setattr(tfcoil_variables, "rbmax", tfcurrentparam.rbmax)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tfcurrentparam.i_tf_sup)

    monkeypatch.setattr(
        tfcoil_variables, "casths_fraction", tfcurrentparam.casths_fraction
    )

    monkeypatch.setattr(tfcoil_variables, "tinstf", tfcurrentparam.tinstf)

    monkeypatch.setattr(tfcoil_variables, "tftort", tfcurrentparam.tftort)

    monkeypatch.setattr(tfcoil_variables, "bmaxtf", tfcurrentparam.bmaxtf)

    monkeypatch.setattr(tfcoil_variables, "tfinsgap", tfcurrentparam.tfinsgap)

    monkeypatch.setattr(
        tfcoil_variables,
        "tfc_sidewall_is_fraction",
        tfcurrentparam.tfc_sidewall_is_fraction,
    )

    monkeypatch.setattr(tfcoil_variables, "casths", tfcurrentparam.casths)

    monkeypatch.setattr(
        tfcoil_variables, "casthi_is_fraction", tfcurrentparam.casthi_is_fraction
    )

    monkeypatch.setattr(
        tfcoil_variables, "casthi_fraction", tfcurrentparam.casthi_fraction
    )

    monkeypatch.setattr(tfcoil_variables, "n_tf", tfcurrentparam.n_tf)

    monkeypatch.setattr(tfcoil_variables, "thicndut", tfcurrentparam.thicndut)

    monkeypatch.setattr(tfcoil_variables, "thkcas", tfcurrentparam.thkcas)

    monkeypatch.setattr(tfcoil_variables, "oacdcp", tfcurrentparam.oacdcp)

    monkeypatch.setattr(tfcoil_variables, "tfareain", tfcurrentparam.tfareain)

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfcurrentparam.r_tf_inboard_out
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfcurrentparam.r_tf_inboard_in
    )

    monkeypatch.setattr(build_variables, "tfcth", tfcurrentparam.tfcth)

    monkeypatch.setattr(physics_variables, "bt", tfcurrentparam.bt)

    monkeypatch.setattr(physics_variables, "rmajor", tfcurrentparam.rmajor)

    monkeypatch.setattr(sctfcoil_module, "tfc_current", tfcurrentparam.tfc_current)

    monkeypatch.setattr(sctfcoil_module, "theta_coil", tfcurrentparam.theta_coil)

    sctfcoil.tf_current()

    assert tfcoil_variables.ritfc == pytest.approx(tfcurrentparam.expected_ritfc)

    assert tfcoil_variables.rbmax == pytest.approx(tfcurrentparam.expected_rbmax)

    assert tfcoil_variables.bmaxtf == pytest.approx(tfcurrentparam.expected_bmaxtf)

    assert tfcoil_variables.oacdcp == pytest.approx(tfcurrentparam.expected_oacdcp)

    assert sctfcoil_module.tfc_current == pytest.approx(
        tfcurrentparam.expected_tfc_current
    )


class TfGlobalGeometryParam(NamedTuple):
    r_tf_outboard_mid: Any = None

    r_cp_top: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    tfthko: Any = None

    tfareain: Any = None

    ritfc: Any = None

    tftort: Any = None

    n_tf: Any = None

    arealeg: Any = None

    i_tf_sup: Any = None

    dztop: Any = None

    i_tf_case_geom: Any = None

    itart: Any = None

    kappa: Any = None

    rminor: Any = None

    h_cp_top: Any = None

    r_tf_outboard_in: Any = None

    r_tf_outboard_out: Any = None

    theta_coil: Any = None

    tan_theta_coil: Any = None

    expected_tfareain: Any = None

    expected_tftort: Any = None

    expected_arealeg: Any = None

    expected_r_tf_outboard_in: Any = None

    expected_r_tf_outboard_out: Any = None

    expected_theta_coil: Any = None

    expected_tan_theta_coil: Any = None


@pytest.mark.parametrize(
    "tfglobalgeometryparam",
    (
        TfGlobalGeometryParam(
            r_tf_outboard_mid=16.519405859443332,
            r_cp_top=4.20194118510911,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            tfthko=1.208,
            tfareain=0,
            ritfc=0,
            tftort=1,
            n_tf=16,
            arealeg=0,
            i_tf_sup=1,
            dztop=0,
            i_tf_case_geom=0,
            itart=0,
            kappa=1.8480000000000001,
            rminor=2.8677741935483869,
            h_cp_top=0,
            r_tf_outboard_in=0,
            r_tf_outboard_out=0,
            theta_coil=0,
            tan_theta_coil=0,
            expected_tfareain=27.308689677971632,
            expected_tftort=1.6395161177915356,
            expected_arealeg=1.9805354702921749,
            expected_r_tf_outboard_in=15.915405859443332,
            expected_r_tf_outboard_out=17.123405859443331,
            expected_theta_coil=0.19634954084936207,
            expected_tan_theta_coil=0.19891236737965801,
        ),
        TfGlobalGeometryParam(
            r_tf_outboard_mid=17.063351291812893,
            r_cp_top=4.4822055399518357,
            r_tf_inboard_in=2.9538679176819831,
            r_tf_inboard_out=4.4822055399518357,
            tfthko=1.5283376222698528,
            tfareain=35.703669036223495,
            ritfc=241812532.66279837,
            tftort=1.7488698442633552,
            n_tf=16,
            arealeg=2.6728635794409041,
            i_tf_sup=1,
            dztop=0,
            i_tf_case_geom=0,
            itart=0,
            kappa=1.8480000000000001,
            rminor=2.9620024998595755,
            h_cp_top=0,
            r_tf_outboard_in=16.299182480677967,
            r_tf_outboard_out=17.827520102947819,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            expected_tfareain=35.703669036223495,
            expected_tftort=1.7488698442633552,
            expected_arealeg=2.6728635794409041,
            expected_r_tf_outboard_in=16.299182480677967,
            expected_r_tf_outboard_out=17.827520102947819,
            expected_theta_coil=0.19634954084936207,
            expected_tan_theta_coil=0.19891236737965801,
        ),
    ),
)
def test_tf_global_geometry(tfglobalgeometryparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_global_geometry.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfglobalgeometryparam: the data used to mock and assert in this test.
    :type tfglobalgeometryparam: tfglobalgeometryparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param sctfcoil: initialised Sctfcoil object
    :type sctfcoil: process.sctfcoil.Sctfcoil
    """

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", tfglobalgeometryparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(build_variables, "r_cp_top", tfglobalgeometryparam.r_cp_top)

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfglobalgeometryparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfglobalgeometryparam.r_tf_inboard_out
    )

    monkeypatch.setattr(build_variables, "tfthko", tfglobalgeometryparam.tfthko)

    monkeypatch.setattr(tfcoil_variables, "tfareain", tfglobalgeometryparam.tfareain)

    monkeypatch.setattr(tfcoil_variables, "ritfc", tfglobalgeometryparam.ritfc)

    monkeypatch.setattr(tfcoil_variables, "tftort", tfglobalgeometryparam.tftort)

    monkeypatch.setattr(tfcoil_variables, "n_tf", tfglobalgeometryparam.n_tf)

    monkeypatch.setattr(tfcoil_variables, "arealeg", tfglobalgeometryparam.arealeg)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tfglobalgeometryparam.i_tf_sup)

    monkeypatch.setattr(tfcoil_variables, "dztop", tfglobalgeometryparam.dztop)

    monkeypatch.setattr(
        tfcoil_variables, "i_tf_case_geom", tfglobalgeometryparam.i_tf_case_geom
    )

    monkeypatch.setattr(physics_variables, "itart", tfglobalgeometryparam.itart)

    monkeypatch.setattr(physics_variables, "kappa", tfglobalgeometryparam.kappa)

    monkeypatch.setattr(physics_variables, "rminor", tfglobalgeometryparam.rminor)

    monkeypatch.setattr(sctfcoil_module, "h_cp_top", tfglobalgeometryparam.h_cp_top)

    monkeypatch.setattr(
        sctfcoil_module, "r_tf_outboard_in", tfglobalgeometryparam.r_tf_outboard_in
    )

    monkeypatch.setattr(
        sctfcoil_module, "r_tf_outboard_out", tfglobalgeometryparam.r_tf_outboard_out
    )

    monkeypatch.setattr(sctfcoil_module, "theta_coil", tfglobalgeometryparam.theta_coil)

    monkeypatch.setattr(
        sctfcoil_module, "tan_theta_coil", tfglobalgeometryparam.tan_theta_coil
    )

    sctfcoil.tf_global_geometry()

    assert tfcoil_variables.tfareain == pytest.approx(
        tfglobalgeometryparam.expected_tfareain
    )

    assert tfcoil_variables.tftort == pytest.approx(
        tfglobalgeometryparam.expected_tftort
    )

    assert tfcoil_variables.arealeg == pytest.approx(
        tfglobalgeometryparam.expected_arealeg
    )

    assert sctfcoil_module.r_tf_outboard_in == pytest.approx(
        tfglobalgeometryparam.expected_r_tf_outboard_in
    )

    assert sctfcoil_module.r_tf_outboard_out == pytest.approx(
        tfglobalgeometryparam.expected_r_tf_outboard_out
    )

    assert sctfcoil_module.theta_coil == pytest.approx(
        tfglobalgeometryparam.expected_theta_coil
    )

    assert sctfcoil_module.tan_theta_coil == pytest.approx(
        tfglobalgeometryparam.expected_tan_theta_coil
    )


@pytest.mark.parametrize(
    "a, b, expected_circumference",
    (
        (2.667950e9, 6.782819e8, 11464316399.111176),
        (4.7186039761812131, 3.6192586838709673, 26.308134540723429),
    ),
    ids=["johndcook", "baseline_2018"],
)
def test_circumference(a, b, expected_circumference, sctfcoil):
    """Unit test for the sctfcoil circumference routine.

    This unit test uses values from an external blog referenced in the
    routine header (https://www.johndcook.com/blog/2013/05/05/ramanujan-circumference-ellipse/)
    as well as results obtained from baseline 2018.

    :param a: the value of a (x/a)^2...  in the formula of an ellipse
    :type a: float

    :param b: the value of b ...(y/b)^2...  in the formula of an ellipse
    :type b: float

    :param expected_circumference: the expected result of the routine given inputs a and b
    :type expected_circumference: float

    :param sctfcoil: initialised Sctfcoil object
    :type sctfcoil: process.sctfcoil.Sctfcoil
    """
    assert sctfcoil.circumference(a, b) == pytest.approx(expected_circumference)


class ResTfInternalGeomParam(NamedTuple):
    n_tf_turn: Any = None

    thicndut: Any = None

    thkcas: Any = None

    dr_tf_wp: Any = None

    tftort: Any = None

    tfareain: Any = None

    ritfc: Any = None

    fcoolcp: Any = None

    cpttf: Any = None

    cdtfleg: Any = None

    casthi: Any = None

    aiwp: Any = None

    acasetf: Any = None

    tinstf: Any = None

    n_tf: Any = None

    tfthko: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    r_cp_top: Any = None

    itart: Any = None

    expected_n_tf_turn: Any = None

    expected_cpttf: Any = None

    expected_cdtfleg: Any = None

    expected_aiwp: Any = None

    expected_acasetf: Any = None


class TfResHeatingParam(NamedTuple):
    rhocp: Any = None
    tlegav: Any = None
    thicndut: Any = None
    th_joint_contact: Any = None
    rhotfleg: Any = None
    vol_cond_cp: Any = None
    n_tf_turn: Any = None
    thkcas: Any = None
    tftort: Any = None
    tfleng: Any = None
    tflegres: Any = None
    tcpav: Any = None
    arealeg: Any = None
    ritfc: Any = None
    rho_tf_joints: Any = None
    presleg: Any = None
    prescp: Any = None
    pres_joints: Any = None
    n_tf_joints_contact: Any = None
    n_tf_joints: Any = None
    n_tf: Any = None
    i_tf_sup: Any = None
    frholeg: Any = None
    frhocp: Any = None
    fcoolcp: Any = None
    casthi: Any = None
    a_cp_cool: Any = None
    fcoolleg: Any = None
    i_cp_joints: Any = None
    tinstf: Any = None
    tfthko: Any = None
    tfcth: Any = None
    r_cp_top: Any = None
    hmax: Any = None
    r_tf_inboard_in: Any = None
    r_tf_inboard_out: Any = None
    itart: Any = None
    h_cp_top: Any = None
    is_leg_cp_temp_same: Any = None
    expected_rhocp: Any = None
    expected_rhotfleg: Any = None
    expected_vol_cond_cp: Any = None
    expected_tflegres: Any = None
    expected_presleg: Any = None
    expected_prescp: Any = None
    expected_pres_joints: Any = None
    expected_a_cp_cool: Any = None
    expected_is_leg_cp_temp_same: Any = None


@pytest.mark.parametrize(
    "restfinternalgeomparam",
    (
        ResTfInternalGeomParam(
            n_tf_turn=0,
            thicndut=0.00080000000000000004,
            thkcas=0,
            dr_tf_wp=0.15483000000000002,
            tftort=0.45367650933034859,
            tfareain=0.0753112923616783,
            ritfc=25500000,
            fcoolcp=0.12725,
            cpttf=70000,
            cdtfleg=0,
            casthi=0.0077415000000000019,
            aiwp=0,
            acasetf=0,
            tinstf=0,
            n_tf=12,
            tfthko=0.15483000000000002,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.87643571428571443,
            itart=1,
            expected_n_tf_turn=1,
            expected_cpttf=2125000,
            expected_cdtfleg=421788350.27812088,
            expected_aiwp=0.00030678028680367151,
            expected_acasetf=0.00061190425043863676,
        ),
        ResTfInternalGeomParam(
            n_tf_turn=1,
            thicndut=0.00080000000000000004,
            thkcas=0,
            dr_tf_wp=0.14708850000000001,
            tftort=0.44435902370665786,
            tfareain=0.0753112923616783,
            ritfc=25500000,
            fcoolcp=0.12725,
            cpttf=2125000,
            cdtfleg=421788350.27812088,
            casthi=0.0077415000000000019,
            aiwp=0.00030678028680367151,
            acasetf=0.00061190425043863676,
            tinstf=0,
            n_tf=12,
            tfthko=0.15483000000000002,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.85843571428571441,
            itart=1,
            expected_n_tf_turn=1,
            expected_cpttf=2125000,
            expected_cdtfleg=430664525.98439038,
            expected_aiwp=0.00029439388680367086,
            expected_acasetf=0.00061190425043863676,
        ),
    ),
)
def test_res_tf_internal_geom(restfinternalgeomparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for res_tf_internal_geom.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param restfinternalgeomparam: the data used to mock and assert in this test.
    :type restfinternalgeomparam: restfinternalgeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch

    :param sctfcoil: initialised Sctfcoil object
    :type sctfcoil: process.sctfcoil.Sctfcoil
    """

    monkeypatch.setattr(tfcoil_variables, "n_tf_turn", restfinternalgeomparam.n_tf_turn)

    monkeypatch.setattr(tfcoil_variables, "thicndut", restfinternalgeomparam.thicndut)

    monkeypatch.setattr(tfcoil_variables, "thkcas", restfinternalgeomparam.thkcas)

    monkeypatch.setattr(tfcoil_variables, "dr_tf_wp", restfinternalgeomparam.dr_tf_wp)

    monkeypatch.setattr(tfcoil_variables, "tftort", restfinternalgeomparam.tftort)

    monkeypatch.setattr(tfcoil_variables, "tfareain", restfinternalgeomparam.tfareain)

    monkeypatch.setattr(tfcoil_variables, "ritfc", restfinternalgeomparam.ritfc)

    monkeypatch.setattr(tfcoil_variables, "fcoolcp", restfinternalgeomparam.fcoolcp)

    monkeypatch.setattr(tfcoil_variables, "cpttf", restfinternalgeomparam.cpttf)

    monkeypatch.setattr(tfcoil_variables, "cdtfleg", restfinternalgeomparam.cdtfleg)

    monkeypatch.setattr(tfcoil_variables, "casthi", restfinternalgeomparam.casthi)

    monkeypatch.setattr(tfcoil_variables, "aiwp", restfinternalgeomparam.aiwp)

    monkeypatch.setattr(tfcoil_variables, "acasetf", restfinternalgeomparam.acasetf)

    monkeypatch.setattr(tfcoil_variables, "tinstf", restfinternalgeomparam.tinstf)

    monkeypatch.setattr(tfcoil_variables, "n_tf", restfinternalgeomparam.n_tf)

    monkeypatch.setattr(build_variables, "tfthko", restfinternalgeomparam.tfthko)

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", restfinternalgeomparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", restfinternalgeomparam.r_tf_inboard_out
    )

    monkeypatch.setattr(build_variables, "r_cp_top", restfinternalgeomparam.r_cp_top)

    monkeypatch.setattr(physics_variables, "itart", restfinternalgeomparam.itart)

    sctfcoil.res_tf_internal_geom()

    assert tfcoil_variables.n_tf_turn == pytest.approx(
        restfinternalgeomparam.expected_n_tf_turn
    )

    assert tfcoil_variables.cpttf == pytest.approx(
        restfinternalgeomparam.expected_cpttf
    )

    assert tfcoil_variables.cdtfleg == pytest.approx(
        restfinternalgeomparam.expected_cdtfleg
    )

    assert tfcoil_variables.aiwp == pytest.approx(restfinternalgeomparam.expected_aiwp)

    assert tfcoil_variables.acasetf == pytest.approx(
        restfinternalgeomparam.expected_acasetf
    )

    assert tfcoil_variables.n_tf_turn == pytest.approx(
        restfinternalgeomparam.expected_n_tf_turn
    )

    assert tfcoil_variables.cpttf == pytest.approx(
        restfinternalgeomparam.expected_cpttf
    )

    assert tfcoil_variables.cdtfleg == pytest.approx(
        restfinternalgeomparam.expected_cdtfleg
    )

    assert tfcoil_variables.aiwp == pytest.approx(restfinternalgeomparam.expected_aiwp)

    assert tfcoil_variables.acasetf == pytest.approx(
        restfinternalgeomparam.expected_acasetf
    )


@pytest.mark.parametrize(
    "tfresheatingparam",
    (
        TfResHeatingParam(
            rhocp=0,
            tlegav=-1,
            thicndut=0.00080000000000000004,
            th_joint_contact=0.029999999999999999,
            rhotfleg=0,
            vol_cond_cp=0,
            n_tf_turn=1,
            thkcas=0,
            tftort=0.45367650933034859,
            tfleng=15.582502857142856,
            tflegres=0,
            tcpav=347.13,
            arealeg=0.070242733939617885,
            ritfc=25500000,
            rho_tf_joints=2.5000000000000002e-10,
            presleg=0,
            prescp=0,
            pres_joints=0,
            n_tf_joints_contact=6,
            n_tf_joints=4,
            n_tf=12,
            i_tf_sup=0,
            frholeg=1,
            frhocp=1,
            fcoolcp=0.12725,
            casthi=0.0077415000000000019,
            a_cp_cool=0,
            fcoolleg=0.20000000000000001,
            i_cp_joints=1,
            tinstf=0,
            tfthko=0.15483000000000002,
            tfcth=0.15483000000000002,
            r_cp_top=0.87643571428571443,
            hmax=4.4214285714285717,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            itart=1,
            h_cp_top=2.6714285714285717,
            is_leg_cp_temp_same=0,
            expected_rhocp=2.1831760869565221e-08,
            expected_rhotfleg=2.1831760869565221e-08,
            expected_vol_cond_cp=12.020160732580297,
            expected_tflegres=6.1387543007600344e-06,
            expected_presleg=332643748.67243439,
            expected_prescp=470083798.99090022,
            expected_pres_joints=1944336.7995005273,
            expected_a_cp_cool=0.00068328705812121333,
            expected_is_leg_cp_temp_same=1,
        ),
        TfResHeatingParam(
            rhocp=2.1831760869565221e-08,
            tlegav=-1,
            thicndut=0.00080000000000000004,
            th_joint_contact=0.029999999999999999,
            rhotfleg=2.1831760869565221e-08,
            vol_cond_cp=12.020160732580297,
            n_tf_turn=1,
            thkcas=0,
            tftort=0.44435902370665786,
            tfleng=15.654502857142857,
            tflegres=6.1387543007600344e-06,
            tcpav=347.13,
            arealeg=0.068800107640501845,
            ritfc=25500000,
            rho_tf_joints=2.5000000000000002e-10,
            presleg=332643748.67243439,
            prescp=470083798.99090022,
            pres_joints=1944336.7995005273,
            n_tf_joints_contact=6,
            n_tf_joints=4,
            n_tf=12,
            i_tf_sup=0,
            frholeg=1,
            frhocp=1,
            fcoolcp=0.12725,
            casthi=0.0077415000000000019,
            a_cp_cool=0.00068328705812121333,
            fcoolleg=0.20000000000000001,
            i_cp_joints=1,
            tinstf=0,
            tfthko=0.15483000000000002,
            tfcth=0.15483000000000002,
            r_cp_top=0.85843571428571441,
            hmax=4.4214285714285717,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            itart=1,
            h_cp_top=2.6714285714285717,
            is_leg_cp_temp_same=1,
            expected_rhocp=2.1831760869565221e-08,
            expected_rhotfleg=2.1831760869565221e-08,
            expected_vol_cond_cp=11.545770024935592,
            expected_tflegres=6.2969005770928158e-06,
            expected_presleg=341213300.02121693,
            expected_prescp=475710489.56122422,
            expected_pres_joints=1944336.7995005273,
            expected_a_cp_cool=0.00068328705812121333,
            expected_is_leg_cp_temp_same=1,
        ),
    ),
)
def test_tf_res_heating(tfresheatingparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_res_heating.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param tfresheatingparam: the data used to mock and assert in this test.
    :type tfresheatingparam: tfresheatingparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(tfcoil_variables, "rhocp", tfresheatingparam.rhocp)

    monkeypatch.setattr(tfcoil_variables, "tlegav", tfresheatingparam.tlegav)

    monkeypatch.setattr(tfcoil_variables, "thicndut", tfresheatingparam.thicndut)

    monkeypatch.setattr(
        tfcoil_variables, "th_joint_contact", tfresheatingparam.th_joint_contact
    )

    monkeypatch.setattr(tfcoil_variables, "rhotfleg", tfresheatingparam.rhotfleg)

    monkeypatch.setattr(tfcoil_variables, "vol_cond_cp", tfresheatingparam.vol_cond_cp)

    monkeypatch.setattr(tfcoil_variables, "n_tf_turn", tfresheatingparam.n_tf_turn)

    monkeypatch.setattr(tfcoil_variables, "thkcas", tfresheatingparam.thkcas)

    monkeypatch.setattr(tfcoil_variables, "tftort", tfresheatingparam.tftort)

    monkeypatch.setattr(tfcoil_variables, "tfleng", tfresheatingparam.tfleng)

    monkeypatch.setattr(tfcoil_variables, "tflegres", tfresheatingparam.tflegres)

    monkeypatch.setattr(tfcoil_variables, "tcpav", tfresheatingparam.tcpav)

    monkeypatch.setattr(tfcoil_variables, "arealeg", tfresheatingparam.arealeg)

    monkeypatch.setattr(tfcoil_variables, "ritfc", tfresheatingparam.ritfc)

    monkeypatch.setattr(
        tfcoil_variables, "rho_tf_joints", tfresheatingparam.rho_tf_joints
    )

    monkeypatch.setattr(tfcoil_variables, "presleg", tfresheatingparam.presleg)

    monkeypatch.setattr(tfcoil_variables, "prescp", tfresheatingparam.prescp)

    monkeypatch.setattr(tfcoil_variables, "pres_joints", tfresheatingparam.pres_joints)

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_joints_contact", tfresheatingparam.n_tf_joints_contact
    )

    monkeypatch.setattr(tfcoil_variables, "n_tf_joints", tfresheatingparam.n_tf_joints)

    monkeypatch.setattr(tfcoil_variables, "n_tf", tfresheatingparam.n_tf)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tfresheatingparam.i_tf_sup)

    monkeypatch.setattr(tfcoil_variables, "frholeg", tfresheatingparam.frholeg)

    monkeypatch.setattr(tfcoil_variables, "frhocp", tfresheatingparam.frhocp)

    monkeypatch.setattr(tfcoil_variables, "fcoolcp", tfresheatingparam.fcoolcp)

    monkeypatch.setattr(tfcoil_variables, "casthi", tfresheatingparam.casthi)

    monkeypatch.setattr(tfcoil_variables, "a_cp_cool", tfresheatingparam.a_cp_cool)

    monkeypatch.setattr(tfcoil_variables, "fcoolleg", tfresheatingparam.fcoolleg)

    monkeypatch.setattr(tfcoil_variables, "i_cp_joints", tfresheatingparam.i_cp_joints)

    monkeypatch.setattr(tfcoil_variables, "tinstf", tfresheatingparam.tinstf)

    monkeypatch.setattr(build_variables, "tfthko", tfresheatingparam.tfthko)

    monkeypatch.setattr(build_variables, "tfcth", tfresheatingparam.tfcth)

    monkeypatch.setattr(build_variables, "r_cp_top", tfresheatingparam.r_cp_top)

    monkeypatch.setattr(build_variables, "hmax", tfresheatingparam.hmax)

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfresheatingparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfresheatingparam.r_tf_inboard_out
    )

    monkeypatch.setattr(physics_variables, "itart", tfresheatingparam.itart)

    monkeypatch.setattr(sctfcoil_module, "h_cp_top", tfresheatingparam.h_cp_top)

    monkeypatch.setattr(
        sctfcoil_module, "is_leg_cp_temp_same", tfresheatingparam.is_leg_cp_temp_same
    )

    sctfcoil.tf_res_heating()

    assert tfcoil_variables.rhocp == pytest.approx(tfresheatingparam.expected_rhocp)

    assert tfcoil_variables.rhotfleg == pytest.approx(
        tfresheatingparam.expected_rhotfleg
    )

    assert tfcoil_variables.vol_cond_cp == pytest.approx(
        tfresheatingparam.expected_vol_cond_cp
    )

    assert tfcoil_variables.tflegres == pytest.approx(
        tfresheatingparam.expected_tflegres
    )

    assert tfcoil_variables.presleg == pytest.approx(tfresheatingparam.expected_presleg)

    assert tfcoil_variables.prescp == pytest.approx(tfresheatingparam.expected_prescp)

    assert tfcoil_variables.pres_joints == pytest.approx(
        tfresheatingparam.expected_pres_joints
    )

    assert tfcoil_variables.a_cp_cool == pytest.approx(
        tfresheatingparam.expected_a_cp_cool
    )

    assert sctfcoil_module.is_leg_cp_temp_same == pytest.approx(
        tfresheatingparam.expected_is_leg_cp_temp_same
    )


class CpostParam(NamedTuple):
    n_tf: Any = None

    hmax: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    r_cp_top: Any = None

    ztop: Any = None

    hmaxi: Any = None

    cas_in_th: Any = None

    cas_out_th: Any = None

    gr_ins_th: Any = None

    ins_th: Any = None

    n_tf_turn: Any = None

    curr: Any = None

    rho: Any = None

    fcool: Any = None

    expected_vol_ins_cp: Any = None

    expected_vol_gr_ins_cp: Any = None

    expected_vol_case_cp: Any = None

    expected_respow: Any = None

    expected_vol_cond_cp: Any = None

    expected_a_cp_cool: Any = None


@pytest.mark.parametrize(
    "cpostparam",
    (
        CpostParam(
            n_tf=12,
            hmax=4.4214285714285717,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.87643571428571443,
            ztop=2.6714285714285717,
            hmaxi=4.5762585714285713,
            cas_in_th=0,
            cas_out_th=0.0077415000000000019,
            gr_ins_th=0,
            ins_th=0.00080000000000000004,
            n_tf_turn=1,
            curr=25500000,
            rho=2.1831760869565221e-08,
            fcool=0.12725,
            expected_vol_ins_cp=0.12917075053120922,
            expected_vol_gr_ins_cp=0,
            expected_vol_case_cp=0.12791418544773489,
            expected_respow=470083798.99090022,
            expected_vol_cond_cp=12.020160732580297,
            expected_a_cp_cool=0.00068328705812121333,
        ),
        CpostParam(
            n_tf=12,
            hmax=4.4214285714285717,
            r_tf_inboard_in=0,
            r_tf_inboard_out=0.15483000000000002,
            r_cp_top=0.85843571428571441,
            ztop=2.6714285714285717,
            hmaxi=4.5762585714285713,
            cas_in_th=0,
            cas_out_th=0.0077415000000000019,
            gr_ins_th=0,
            ins_th=0.00080000000000000004,
            n_tf_turn=1,
            curr=25500000,
            rho=2.1831760869565221e-08,
            fcool=0.12725,
            expected_vol_ins_cp=0.12679799009998483,
            expected_vol_gr_ins_cp=0,
            expected_vol_case_cp=0.12648575512245444,
            expected_respow=475710489.56122422,
            expected_vol_cond_cp=11.545770024935592,
            expected_a_cp_cool=0.00068328705812121333,
        ),
    ),
)
def test_cpost(cpostparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for cpost.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param cpostparam: the data used to mock and assert in this test.
    :type cpostparam: cpostparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "hmax", cpostparam.hmax)

    (
        a_cp_cool,
        vol_cond_cp,
        respow,
        vol_ins_cp,
        vol_case_cp,
        vol_gr_ins_cp,
    ) = sctfcoil.cpost(
        r_tf_inboard_in=cpostparam.r_tf_inboard_in,
        r_tf_inboard_out=cpostparam.r_tf_inboard_out,
        r_cp_top=cpostparam.r_cp_top,
        ztop=cpostparam.ztop,
        hmaxi=cpostparam.hmaxi,
        cas_in_th=cpostparam.cas_in_th,
        cas_out_th=cpostparam.cas_out_th,
        gr_ins_th=cpostparam.gr_ins_th,
        ins_th=cpostparam.ins_th,
        n_tf_turn=cpostparam.n_tf_turn,
        curr=cpostparam.curr,
        rho=cpostparam.rho,
        fcool=cpostparam.fcool,
        n_tf=cpostparam.n_tf,
    )

    assert vol_ins_cp == pytest.approx(cpostparam.expected_vol_ins_cp)

    assert vol_gr_ins_cp == pytest.approx(cpostparam.expected_vol_gr_ins_cp)

    assert vol_case_cp == pytest.approx(cpostparam.expected_vol_case_cp)

    assert respow == pytest.approx(cpostparam.expected_respow)

    assert vol_cond_cp == pytest.approx(cpostparam.expected_vol_cond_cp)

    assert a_cp_cool == pytest.approx(cpostparam.expected_a_cp_cool)


class TfFieldAndForceParam(NamedTuple):
    rminor: Any = None

    rmajor: Any = None

    bt: Any = None

    itart: Any = None

    r_tf_outboard_mid: Any = None

    r_vv_inboard_out: Any = None

    r_tf_inboard_mid: Any = None

    r_cp_top: Any = None

    vforce: Any = None

    n_tf: Any = None

    taucq: Any = None

    sigvvall: Any = None

    cforce: Any = None

    ritfc: Any = None

    bmaxtf: Any = None

    i_tf_sup: Any = None

    f_vforce_inboard: Any = None

    vforce_outboard: Any = None

    tinstf: Any = None

    thicndut: Any = None

    dr_tf_wp: Any = None

    tfinsgap: Any = None

    i_cp_joints: Any = None

    casthi: Any = None

    r_tf_outboard_in: Any = None

    r_wp_inner: Any = None

    r_wp_outer: Any = None

    vforce_inboard_tot: Any = None

    expected_vforce: Any = None

    expected_cforce: Any = None

    expected_f_vforce_inboard: Any = None

    expected_vforce_outboard: Any = None

    expected_vforce_inboard_tot: Any = None


@pytest.mark.parametrize(
    "tffieldandforceparam",
    (
        TfFieldAndForceParam(
            rminor=0.97142857142857153,
            rmajor=1.7000000000000002,
            bt=3,
            itart=1,
            r_tf_outboard_mid=4.1688435714285719,
            r_vv_inboard_out=0.20483000000000001,
            r_tf_inboard_mid=0.077415000000000012,
            r_cp_top=0.87643571428571443,
            vforce=0,
            n_tf=12,
            taucq=30,
            sigvvall=93000000,
            cforce=0,
            ritfc=25500000,
            bmaxtf=34.862617362267024,
            i_tf_sup=0,
            f_vforce_inboard=0.5,
            vforce_outboard=0,
            tinstf=0,
            thicndut=0.00080000000000000004,
            dr_tf_wp=0.15483000000000002,
            tfinsgap=0.01,
            i_cp_joints=1,
            casthi=0.0077415000000000019,
            r_tf_outboard_in=4.0914285714285716,
            r_wp_inner=0,
            r_wp_outer=0.14708850000000001,
            vforce_inboard_tot=0,
            expected_vforce=12380916.66459452,
            expected_cforce=37041530.947408713,
            expected_f_vforce_inboard=0.59539634897566385,
            expected_vforce_outboard=8413494.7991220243,
            expected_vforce_inboard_tot=148570999.97513425,
        ),
        TfFieldAndForceParam(
            rminor=0.97142857142857153,
            rmajor=1.7000000000000002,
            bt=3,
            itart=1,
            r_tf_outboard_mid=4.1868435714285717,
            r_vv_inboard_out=0.20483000000000001,
            r_tf_inboard_mid=0.077415000000000012,
            r_cp_top=0.85843571428571441,
            vforce=12380916.66459452,
            n_tf=12,
            taucq=30,
            sigvvall=93000000,
            cforce=37041530.947408713,
            ritfc=25500000,
            bmaxtf=34.862617362267024,
            i_tf_sup=0,
            f_vforce_inboard=0.59539634897566385,
            vforce_outboard=8413494.7991220243,
            tinstf=0,
            thicndut=0.00080000000000000004,
            dr_tf_wp=0.14708850000000001,
            tfinsgap=0.01,
            i_cp_joints=1,
            casthi=0.0077415000000000019,
            r_tf_outboard_in=4.1094285714285714,
            r_wp_inner=0,
            r_wp_outer=0.14708850000000001,
            vforce_inboard_tot=148570999.97513425,
            expected_vforce=12268469.138442248,
            expected_cforce=37041530.947408713,
            expected_f_vforce_inboard=0.58932254522566518,
            expected_vforce_outboard=8549450.0771621168,
            expected_vforce_inboard_tot=147221629.66130698,
        ),
    ),
)
def test_tf_field_and_force(tffieldandforceparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_field_and_force.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param tffieldandforceparam: the data used to mock and assert in this test.
    :type tffieldandforceparam: tffieldandforceparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "rminor", tffieldandforceparam.rminor)

    monkeypatch.setattr(physics_variables, "rmajor", tffieldandforceparam.rmajor)

    monkeypatch.setattr(physics_variables, "bt", tffieldandforceparam.bt)

    monkeypatch.setattr(physics_variables, "itart", tffieldandforceparam.itart)

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", tffieldandforceparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(
        build_variables, "r_vv_inboard_out", tffieldandforceparam.r_vv_inboard_out
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_mid", tffieldandforceparam.r_tf_inboard_mid
    )

    monkeypatch.setattr(build_variables, "r_cp_top", tffieldandforceparam.r_cp_top)

    monkeypatch.setattr(tfcoil_variables, "vforce", tffieldandforceparam.vforce)

    monkeypatch.setattr(tfcoil_variables, "n_tf", tffieldandforceparam.n_tf)

    monkeypatch.setattr(tfcoil_variables, "cforce", tffieldandforceparam.cforce)

    monkeypatch.setattr(tfcoil_variables, "ritfc", tffieldandforceparam.ritfc)

    monkeypatch.setattr(tfcoil_variables, "bmaxtf", tffieldandforceparam.bmaxtf)

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tffieldandforceparam.i_tf_sup)

    monkeypatch.setattr(
        tfcoil_variables, "f_vforce_inboard", tffieldandforceparam.f_vforce_inboard
    )

    monkeypatch.setattr(
        tfcoil_variables, "vforce_outboard", tffieldandforceparam.vforce_outboard
    )

    monkeypatch.setattr(tfcoil_variables, "tinstf", tffieldandforceparam.tinstf)

    monkeypatch.setattr(tfcoil_variables, "thicndut", tffieldandforceparam.thicndut)

    monkeypatch.setattr(tfcoil_variables, "dr_tf_wp", tffieldandforceparam.dr_tf_wp)

    monkeypatch.setattr(tfcoil_variables, "tfinsgap", tffieldandforceparam.tfinsgap)

    monkeypatch.setattr(
        tfcoil_variables, "i_cp_joints", tffieldandforceparam.i_cp_joints
    )

    monkeypatch.setattr(tfcoil_variables, "casthi", tffieldandforceparam.casthi)

    monkeypatch.setattr(
        sctfcoil_module, "r_tf_outboard_in", tffieldandforceparam.r_tf_outboard_in
    )

    monkeypatch.setattr(sctfcoil_module, "r_wp_inner", tffieldandforceparam.r_wp_inner)

    monkeypatch.setattr(sctfcoil_module, "r_wp_outer", tffieldandforceparam.r_wp_outer)

    monkeypatch.setattr(
        sctfcoil_module, "vforce_inboard_tot", tffieldandforceparam.vforce_inboard_tot
    )

    sctfcoil.tf_field_and_force()

    assert tfcoil_variables.vforce == pytest.approx(
        tffieldandforceparam.expected_vforce
    )

    assert tfcoil_variables.cforce == pytest.approx(
        tffieldandforceparam.expected_cforce
    )

    assert tfcoil_variables.f_vforce_inboard == pytest.approx(
        tffieldandforceparam.expected_f_vforce_inboard
    )

    assert tfcoil_variables.vforce_outboard == pytest.approx(
        tffieldandforceparam.expected_vforce_outboard
    )

    assert sctfcoil_module.vforce_inboard_tot == pytest.approx(
        tffieldandforceparam.expected_vforce_inboard_tot
    )


class TfcindParam(NamedTuple):
    yarc: Any = None

    xarc: Any = None

    tfind: Any = None

    tfthk: Any = None

    expected_yarc: Any = None

    expected_tfind: Any = None


@pytest.mark.parametrize(
    "tfcindparam",
    (
        TfcindParam(
            yarc=numpy.array(
                (
                    4.5228880258064512,
                    7.5381467096774184,
                    0,
                    -9.0730900215620327,
                    -5.4438540129372193,
                ),
                order="F",
            ),
            xarc=numpy.array(
                (
                    4.20194118510911,
                    8.316545161290323,
                    15.915405859443332,
                    8.316545161290323,
                    4.20194118510911,
                ),
                order="F",
            ),
            tfind=0,
            tfthk=1.208,
            expected_tfind=5.4453892599192845e-06,
        ),
        TfcindParam(
            yarc=numpy.array(
                (
                    4.5336880258064509,
                    7.5561467096774191,
                    0,
                    -9.0730900215620327,
                    -5.4438540129372193,
                ),
                order="F",
            ),
            xarc=numpy.array(
                (
                    4.20194118510911,
                    8.316545161290323,
                    15.915405859443332,
                    8.316545161290323,
                    4.20194118510911,
                ),
                order="F",
            ),
            tfind=5.4524893280368181e-06,
            tfthk=1.208,
            expected_tfind=5.4524893280368181e-06,
        ),
    ),
)
def test_tfcind(tfcindparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tfcind.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcindparam: the data used to mock and assert in this test.
    :type tfcindparam: tfcindparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(tfcoil_variables, "tfind", tfcindparam.tfind)

    tfind = sctfcoil.tfcind(
        tfthk=tfcindparam.tfthk, xarc=tfcindparam.xarc, yarc=tfcindparam.yarc
    )

    assert tfind == pytest.approx(tfcindparam.expected_tfind)


class TfCoilAreaAndMassesParam(NamedTuple):
    hr1: Any = None

    r_tf_outboard_mid: Any = None

    tfcth: Any = None

    r_tf_inboard_mid: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    hmax: Any = None

    denstl: Any = None

    whtconsh: Any = None

    whttf: Any = None

    whtcas: Any = None

    tficrn: Any = None

    tfcryoarea: Any = None

    tfsao: Any = None

    whtgw: Any = None

    tfocrn: Any = None

    whtconsc: Any = None

    whtconcu: Any = None

    whtcon: Any = None

    whtconin: Any = None

    tfsai: Any = None

    vftf: Any = None

    dcond: Any = None

    dcondins: Any = None

    tfleng: Any = None

    dcase: Any = None

    acndttf: Any = None

    n_tf_turn: Any = None

    n_tf: Any = None

    aiwp: Any = None

    acasetfo: Any = None

    acasetf: Any = None

    fcutfsu: Any = None

    awphec: Any = None

    acstf: Any = None

    whttflgs: Any = None

    whtcp: Any = None

    whtconal: Any = None

    vol_cond_cp: Any = None

    i_tf_sup: Any = None

    i_tf_sc_mat: Any = None

    arealeg: Any = None

    thkcas: Any = None

    voltfleg: Any = None

    cplen: Any = None

    itart: Any = None

    awpc: Any = None

    awptf: Any = None

    vol_ins_cp: Any = None

    vol_gr_ins_cp: Any = None

    vol_case_cp: Any = None

    a_leg_ins: Any = None

    a_leg_gr_ins: Any = None

    a_leg_cond: Any = None

    theta_coil: Any = None

    tan_theta_coil: Any = None

    expected_whtconsh: Any = None

    expected_whtcas: Any = None

    expected_tficrn: Any = None

    expected_tfcryoarea: Any = None

    expected_tfsao: Any = None

    expected_whtgw: Any = None

    expected_tfocrn: Any = None

    expected_whtconsc: Any = None

    expected_whtconcu: Any = None

    expected_whtcon: Any = None

    expected_whtconin: Any = None

    expected_cplen: Any = None


@pytest.mark.parametrize(
    "tfcoilareaandmassesparam",
    (
        TfCoilAreaAndMassesParam(
            hr1=0,
            r_tf_outboard_mid=16.519405859443332,
            tfcth=1.208,
            r_tf_inboard_mid=3.5979411851091103,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            hmax=9.0730900215620327,
            denstl=7800,
            whtconsh=0,
            whttf=0,
            whtcas=0,
            tficrn=0,
            tfcryoarea=0,
            tfsao=0,
            whtgw=0,
            tfocrn=0,
            whtconsc=0,
            whtconcu=0,
            whtcon=0,
            whtconin=0,
            tfsai=0,
            vftf=0.30000000000000004,
            dcond=numpy.array(
                numpy.array(
                    (6080, 6080, 6070, 6080, 6080, 8500, 6070, 8500, 8500), order="F"
                ),
                order="F",
            ).transpose(),
            dcondins=1800,
            tfleng=50.483843027201402,
            dcase=8000,
            acndttf=0.0014685061538103825,
            n_tf_turn=200,
            n_tf=16,
            aiwp=0.087880174466980876,
            acasetfo=1.2752592893394648,
            acasetf=1.0015169239205168,
            fcutfsu=0.80884,
            awphec=0.015707963267948974,
            acstf=0.001293323051622732,
            whttflgs=0,
            whtcp=0,
            whtconal=0,
            vol_cond_cp=0,
            i_tf_sup=1,
            i_tf_sc_mat=5,
            arealeg=1.9805354702921749,
            thkcas=0.52465000000000006,
            voltfleg=0,
            cplen=0,
            itart=0,
            awpc=0.70527618095271016,
            awptf=0.64024601555360383,
            vol_ins_cp=0,
            vol_gr_ins_cp=0,
            vol_case_cp=0,
            a_leg_ins=0,
            a_leg_gr_ins=0,
            a_leg_cond=0,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            expected_whtconsh=115651.90127937049,
            expected_whtcas=1034021.9996272125,
            expected_tficrn=0.8197580588957678,
            expected_tfcryoarea=6381.2092203414386,
            expected_tfsao=1324.3051892984724,
            expected_whtgw=5909.3507916745702,
            expected_tfocrn=0.59553192892551199,
            expected_whtconsc=5802.5700395134345,
            expected_whtconcu=58744.465423173802,
            expected_whtcon=188184.68882144717,
            expected_whtconin=7985.7520793894437,
            expected_cplen=20.562180043124066,
        ),
        TfCoilAreaAndMassesParam(
            hr1=0,
            r_tf_outboard_mid=16.519405859443332,
            tfcth=1.208,
            r_tf_inboard_mid=3.5979411851091103,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            hmax=9.0730900215620327,
            denstl=7800,
            whtconsh=115651.90127937049,
            whttf=19649856.627845347,
            whtcas=1034021.9996272125,
            tficrn=0.8197580588957678,
            tfcryoarea=6381.2092203414386,
            tfsao=1324.3051892984724,
            whtgw=5909.3507916745702,
            tfocrn=0.59553192892551199,
            whtconsc=5802.5700395134345,
            whtconcu=58744.465423173802,
            whtcon=0,
            whtconin=0,
            tfsai=0,
            vftf=0.30000000000000004,
            dcond=numpy.array(
                numpy.array(
                    (6080, 6080, 6070, 6080, 6080, 8500, 6070, 8500, 8500), order="F"
                ),
                order="F",
            ).transpose(),
            dcondins=1800,
            tfleng=50.514015976170839,
            dcase=8000,
            acndttf=0.0014685061538103825,
            n_tf_turn=200,
            n_tf=16,
            aiwp=0.087880174466980876,
            acasetfo=1.2752592893394648,
            acasetf=1.0015169239205168,
            fcutfsu=0.80884,
            awphec=0.015707963267948974,
            acstf=0.001293323051622732,
            whttflgs=0,
            whtcp=0,
            whtconal=0,
            vol_cond_cp=0,
            i_tf_sup=1,
            i_tf_sc_mat=5,
            arealeg=1.9805354702921749,
            thkcas=0.52465000000000006,
            voltfleg=0,
            cplen=20.562180043124066,
            itart=0,
            awpc=0.70527618095271016,
            awptf=0.64024601555360383,
            vol_ins_cp=0,
            vol_gr_ins_cp=0,
            vol_case_cp=0,
            a_leg_ins=0,
            a_leg_gr_ins=0,
            a_leg_cond=0,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            expected_whtconsh=115721.02357090525,
            expected_whtcas=1034699.2182961091,
            expected_tficrn=0.8197580588957678,
            expected_tfcryoarea=6385.0231118485681,
            expected_tfsao=1325.0966938769795,
            expected_whtgw=5912.8826650262808,
            expected_tfocrn=0.59553192892551199,
            expected_whtconsc=5806.038092640837,
            expected_whtconcu=58779.575542593491,
            expected_whtcon=188297.16217276,
            expected_whtconin=7990.5249666247555,
            expected_cplen=20.562180043124066,
        ),
    ),
)
def test_tf_coil_area_and_masses(tfcoilareaandmassesparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_coil_area_and_masses.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcoilareaandmassesparam: the data used to mock and assert in this test.
    :type tfcoilareaandmassesparam: tfcoilareaandmassesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "hr1", tfcoilareaandmassesparam.hr1)

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", tfcoilareaandmassesparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(build_variables, "tfcth", tfcoilareaandmassesparam.tfcth)

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_mid", tfcoilareaandmassesparam.r_tf_inboard_mid
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfcoilareaandmassesparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfcoilareaandmassesparam.r_tf_inboard_out
    )

    monkeypatch.setattr(build_variables, "hmax", tfcoilareaandmassesparam.hmax)

    monkeypatch.setattr(fwbs_variables, "denstl", tfcoilareaandmassesparam.denstl)

    monkeypatch.setattr(tfcoil_variables, "whtconsh", tfcoilareaandmassesparam.whtconsh)

    monkeypatch.setattr(tfcoil_variables, "whttf", tfcoilareaandmassesparam.whttf)

    monkeypatch.setattr(tfcoil_variables, "whtcas", tfcoilareaandmassesparam.whtcas)

    monkeypatch.setattr(tfcoil_variables, "tficrn", tfcoilareaandmassesparam.tficrn)

    monkeypatch.setattr(
        tfcoil_variables, "tfcryoarea", tfcoilareaandmassesparam.tfcryoarea
    )

    monkeypatch.setattr(tfcoil_variables, "tfsao", tfcoilareaandmassesparam.tfsao)

    monkeypatch.setattr(tfcoil_variables, "whtgw", tfcoilareaandmassesparam.whtgw)

    monkeypatch.setattr(tfcoil_variables, "tfocrn", tfcoilareaandmassesparam.tfocrn)

    monkeypatch.setattr(tfcoil_variables, "whtconsc", tfcoilareaandmassesparam.whtconsc)

    monkeypatch.setattr(tfcoil_variables, "whtconcu", tfcoilareaandmassesparam.whtconcu)

    monkeypatch.setattr(tfcoil_variables, "whtcon", tfcoilareaandmassesparam.whtcon)

    monkeypatch.setattr(tfcoil_variables, "whtconin", tfcoilareaandmassesparam.whtconin)

    monkeypatch.setattr(tfcoil_variables, "tfsai", tfcoilareaandmassesparam.tfsai)

    monkeypatch.setattr(tfcoil_variables, "vftf", tfcoilareaandmassesparam.vftf)

    monkeypatch.setattr(tfcoil_variables, "dcond", tfcoilareaandmassesparam.dcond)

    monkeypatch.setattr(tfcoil_variables, "dcondins", tfcoilareaandmassesparam.dcondins)

    monkeypatch.setattr(tfcoil_variables, "tfleng", tfcoilareaandmassesparam.tfleng)

    monkeypatch.setattr(tfcoil_variables, "dcase", tfcoilareaandmassesparam.dcase)

    monkeypatch.setattr(tfcoil_variables, "acndttf", tfcoilareaandmassesparam.acndttf)

    monkeypatch.setattr(
        tfcoil_variables, "n_tf_turn", tfcoilareaandmassesparam.n_tf_turn
    )

    monkeypatch.setattr(tfcoil_variables, "n_tf", tfcoilareaandmassesparam.n_tf)

    monkeypatch.setattr(tfcoil_variables, "aiwp", tfcoilareaandmassesparam.aiwp)

    monkeypatch.setattr(tfcoil_variables, "acasetfo", tfcoilareaandmassesparam.acasetfo)

    monkeypatch.setattr(tfcoil_variables, "acasetf", tfcoilareaandmassesparam.acasetf)

    monkeypatch.setattr(tfcoil_variables, "fcutfsu", tfcoilareaandmassesparam.fcutfsu)

    monkeypatch.setattr(tfcoil_variables, "awphec", tfcoilareaandmassesparam.awphec)

    monkeypatch.setattr(tfcoil_variables, "acstf", tfcoilareaandmassesparam.acstf)

    monkeypatch.setattr(tfcoil_variables, "whttflgs", tfcoilareaandmassesparam.whttflgs)

    monkeypatch.setattr(tfcoil_variables, "whtcp", tfcoilareaandmassesparam.whtcp)

    monkeypatch.setattr(tfcoil_variables, "whtconal", tfcoilareaandmassesparam.whtconal)

    monkeypatch.setattr(
        tfcoil_variables, "vol_cond_cp", tfcoilareaandmassesparam.vol_cond_cp
    )

    monkeypatch.setattr(tfcoil_variables, "i_tf_sup", tfcoilareaandmassesparam.i_tf_sup)

    monkeypatch.setattr(
        tfcoil_variables, "i_tf_sc_mat", tfcoilareaandmassesparam.i_tf_sc_mat
    )

    monkeypatch.setattr(tfcoil_variables, "arealeg", tfcoilareaandmassesparam.arealeg)

    monkeypatch.setattr(tfcoil_variables, "thkcas", tfcoilareaandmassesparam.thkcas)

    monkeypatch.setattr(tfcoil_variables, "voltfleg", tfcoilareaandmassesparam.voltfleg)

    monkeypatch.setattr(tfcoil_variables, "cplen", tfcoilareaandmassesparam.cplen)

    monkeypatch.setattr(physics_variables, "itart", tfcoilareaandmassesparam.itart)

    monkeypatch.setattr(sctfcoil_module, "awpc", tfcoilareaandmassesparam.awpc)

    monkeypatch.setattr(sctfcoil_module, "awptf", tfcoilareaandmassesparam.awptf)

    monkeypatch.setattr(
        sctfcoil_module, "vol_ins_cp", tfcoilareaandmassesparam.vol_ins_cp
    )

    monkeypatch.setattr(
        sctfcoil_module, "vol_gr_ins_cp", tfcoilareaandmassesparam.vol_gr_ins_cp
    )

    monkeypatch.setattr(
        sctfcoil_module, "vol_case_cp", tfcoilareaandmassesparam.vol_case_cp
    )

    monkeypatch.setattr(
        sctfcoil_module, "a_leg_ins", tfcoilareaandmassesparam.a_leg_ins
    )

    monkeypatch.setattr(
        sctfcoil_module, "a_leg_gr_ins", tfcoilareaandmassesparam.a_leg_gr_ins
    )

    monkeypatch.setattr(
        sctfcoil_module, "a_leg_cond", tfcoilareaandmassesparam.a_leg_cond
    )

    monkeypatch.setattr(
        sctfcoil_module, "theta_coil", tfcoilareaandmassesparam.theta_coil
    )

    monkeypatch.setattr(
        sctfcoil_module, "tan_theta_coil", tfcoilareaandmassesparam.tan_theta_coil
    )

    sctfcoil.tf_coil_area_and_masses()

    assert tfcoil_variables.whtconsh == pytest.approx(
        tfcoilareaandmassesparam.expected_whtconsh
    )

    assert tfcoil_variables.whtcas == pytest.approx(
        tfcoilareaandmassesparam.expected_whtcas
    )

    assert tfcoil_variables.tficrn == pytest.approx(
        tfcoilareaandmassesparam.expected_tficrn
    )

    assert tfcoil_variables.tfcryoarea == pytest.approx(
        tfcoilareaandmassesparam.expected_tfcryoarea
    )

    assert tfcoil_variables.tfsao == pytest.approx(
        tfcoilareaandmassesparam.expected_tfsao
    )

    assert tfcoil_variables.whtgw == pytest.approx(
        tfcoilareaandmassesparam.expected_whtgw
    )

    assert tfcoil_variables.tfocrn == pytest.approx(
        tfcoilareaandmassesparam.expected_tfocrn
    )

    assert tfcoil_variables.whtconsc == pytest.approx(
        tfcoilareaandmassesparam.expected_whtconsc
    )

    assert tfcoil_variables.whtconcu == pytest.approx(
        tfcoilareaandmassesparam.expected_whtconcu
    )

    assert tfcoil_variables.whtcon == pytest.approx(
        tfcoilareaandmassesparam.expected_whtcon
    )

    assert tfcoil_variables.whtconin == pytest.approx(
        tfcoilareaandmassesparam.expected_whtconin
    )

    assert tfcoil_variables.cplen == pytest.approx(
        tfcoilareaandmassesparam.expected_cplen
    )


class PeakTfWithRippleParam(NamedTuple):
    tf_fit_t: Any = None

    tf_fit_z: Any = None

    tf_fit_y: Any = None

    n_tf: Any = None

    wwp1: Any = None

    dr_tf_wp: Any = None

    tfin: Any = None

    bmaxtf: Any = None

    expected_tf_fit_t: Any = None

    expected_tf_fit_z: Any = None

    expected_tf_fit_y: Any = None

    expected_bmaxtfrp: Any = None

    expected_flag: Any = None


@pytest.mark.parametrize(
    "peaktfwithrippleparam",
    (
        PeakTfWithRippleParam(
            tf_fit_t=0,
            tf_fit_z=0,
            tf_fit_y=0,
            n_tf=16,
            wwp1=1.299782604942499,
            dr_tf_wp=0.50661087836601015,
            tfin=3.789896624292115,
            bmaxtf=11.717722779177526,
            expected_tf_fit_t=0.80807838916035957,
            expected_tf_fit_z=0.3149613642807837,
            expected_tf_fit_y=1.0658869305062604,
            expected_bmaxtfrp=12.48976756562082,
            expected_flag=0,
        ),
        PeakTfWithRippleParam(
            tf_fit_t=0.80807838916035957,
            tf_fit_z=0.3149613642807837,
            tf_fit_y=1.0658869305062604,
            n_tf=16,
            wwp1=1.299782604942499,
            dr_tf_wp=0.50661087836601015,
            tfin=3.789896624292115,
            bmaxtf=11.717722779177526,
            expected_tf_fit_t=0.80807838916035957,
            expected_tf_fit_z=0.3149613642807837,
            expected_tf_fit_y=1.0658869305062604,
            expected_bmaxtfrp=12.48976756562082,
            expected_flag=0,
        ),
    ),
)
def test_peak_tf_with_ripple(peaktfwithrippleparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for peak_tf_with_ripple.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param peaktfwithrippleparam: the data used to mock and assert in this test.
    :type peaktfwithrippleparam: peaktfwithrippleparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(sctfcoil_module, "tf_fit_t", peaktfwithrippleparam.tf_fit_t)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_z", peaktfwithrippleparam.tf_fit_z)

    monkeypatch.setattr(sctfcoil_module, "tf_fit_y", peaktfwithrippleparam.tf_fit_y)

    bmaxtfrp, flag = sctfcoil.peak_tf_with_ripple(
        n_tf=peaktfwithrippleparam.n_tf,
        wwp1=peaktfwithrippleparam.wwp1,
        dr_tf_wp=peaktfwithrippleparam.dr_tf_wp,
        tfin=peaktfwithrippleparam.tfin,
        bmaxtf=peaktfwithrippleparam.bmaxtf,
    )

    assert sctfcoil_module.tf_fit_t == pytest.approx(
        peaktfwithrippleparam.expected_tf_fit_t
    )

    assert sctfcoil_module.tf_fit_z == pytest.approx(
        peaktfwithrippleparam.expected_tf_fit_z
    )

    assert sctfcoil_module.tf_fit_y == pytest.approx(
        peaktfwithrippleparam.expected_tf_fit_y
    )

    assert bmaxtfrp == pytest.approx(peaktfwithrippleparam.expected_bmaxtfrp)

    assert flag == pytest.approx(peaktfwithrippleparam.expected_flag)


class TfWpGeomParam(NamedTuple):
    tfcth: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    dr_tf_wp: Any = None

    casthi: Any = None

    thkcas: Any = None

    casths: Any = None

    wwp1: Any = None

    wwp2: Any = None

    tinstf: Any = None

    tfinsgap: Any = None

    awpc: Any = None

    awptf: Any = None

    r_wp_inner: Any = None

    r_wp_outer: Any = None

    r_wp_centre: Any = None

    t_wp_toroidal: Any = None

    t_wp_toroidal_av: Any = None

    a_ground_ins: Any = None

    theta_coil: Any = None

    tan_theta_coil: Any = None

    i_tf_wp_geom: Any = None

    expected_wwp1: Any = None

    expected_awpc: Any = None

    expected_awptf: Any = None

    expected_r_wp_inner: Any = None

    expected_r_wp_outer: Any = None

    expected_r_wp_centre: Any = None

    expected_t_wp_toroidal: Any = None

    expected_t_wp_toroidal_av: Any = None

    expected_a_ground_ins: Any = None


@pytest.mark.parametrize(
    "tfwpgeomparam",
    (
        TfWpGeomParam(
            tfcth=1.208,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            dr_tf_wp=0.54261087836601019,
            casthi=0.060000000000000012,
            thkcas=0.52465000000000006,
            casths=0.05000000000000001,
            wwp1=0,
            wwp2=0,
            tinstf=0.0080000000000000019,
            tfinsgap=0.01,
            awpc=0,
            awptf=0,
            r_wp_inner=0,
            r_wp_outer=0,
            r_wp_centre=0,
            t_wp_toroidal=0,
            t_wp_toroidal_av=0,
            a_ground_ins=0,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=0,
            expected_wwp1=1.299782604942499,
            expected_awpc=0.70527618095271016,
            expected_awptf=0.64024601555360383,
            expected_r_wp_inner=3.5185911851091101,
            expected_r_wp_outer=4.06120206347512,
            expected_r_wp_centre=3.789896624292115,
            expected_t_wp_toroidal=1.299782604942499,
            expected_t_wp_toroidal_av=1.299782604942499,
            expected_a_ground_ins=0.028582295732936136,
        ),
        TfWpGeomParam(
            tfcth=1.208,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            dr_tf_wp=0.54261087836601019,
            casthi=0.060000000000000012,
            thkcas=0.52465000000000006,
            casths=0.05000000000000001,
            wwp1=1.299782604942499,
            wwp2=0,
            tinstf=0.0080000000000000019,
            tfinsgap=0.01,
            awpc=0.70527618095271016,
            awptf=0.64024601555360383,
            r_wp_inner=3.5185911851091101,
            r_wp_outer=4.06120206347512,
            r_wp_centre=3.789896624292115,
            t_wp_toroidal=1.299782604942499,
            t_wp_toroidal_av=1.299782604942499,
            a_ground_ins=0.028582295732936136,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=0,
            expected_wwp1=1.299782604942499,
            expected_awpc=0.70527618095271016,
            expected_awptf=0.64024601555360383,
            expected_r_wp_inner=3.5185911851091101,
            expected_r_wp_outer=4.06120206347512,
            expected_r_wp_centre=3.789896624292115,
            expected_t_wp_toroidal=1.299782604942499,
            expected_t_wp_toroidal_av=1.299782604942499,
            expected_a_ground_ins=0.028582295732936136,
        ),
    ),
)
def test_tf_wp_geom(tfwpgeomparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_wp_geom.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfwpgeomparam: the data used to mock and assert in this test.
    :type tfwpgeomparam: tfwpgeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "tfcth", tfwpgeomparam.tfcth)

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfwpgeomparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfwpgeomparam.r_tf_inboard_out
    )

    monkeypatch.setattr(tfcoil_variables, "dr_tf_wp", tfwpgeomparam.dr_tf_wp)

    monkeypatch.setattr(tfcoil_variables, "casthi", tfwpgeomparam.casthi)

    monkeypatch.setattr(tfcoil_variables, "thkcas", tfwpgeomparam.thkcas)

    monkeypatch.setattr(tfcoil_variables, "casths", tfwpgeomparam.casths)

    monkeypatch.setattr(tfcoil_variables, "wwp1", tfwpgeomparam.wwp1)

    monkeypatch.setattr(tfcoil_variables, "wwp2", tfwpgeomparam.wwp2)

    monkeypatch.setattr(tfcoil_variables, "tinstf", tfwpgeomparam.tinstf)

    monkeypatch.setattr(tfcoil_variables, "tfinsgap", tfwpgeomparam.tfinsgap)

    monkeypatch.setattr(sctfcoil_module, "awpc", tfwpgeomparam.awpc)

    monkeypatch.setattr(sctfcoil_module, "awptf", tfwpgeomparam.awptf)

    monkeypatch.setattr(sctfcoil_module, "r_wp_inner", tfwpgeomparam.r_wp_inner)

    monkeypatch.setattr(sctfcoil_module, "r_wp_outer", tfwpgeomparam.r_wp_outer)

    monkeypatch.setattr(sctfcoil_module, "r_wp_centre", tfwpgeomparam.r_wp_centre)

    monkeypatch.setattr(sctfcoil_module, "t_wp_toroidal", tfwpgeomparam.t_wp_toroidal)

    monkeypatch.setattr(
        sctfcoil_module, "t_wp_toroidal_av", tfwpgeomparam.t_wp_toroidal_av
    )

    monkeypatch.setattr(sctfcoil_module, "a_ground_ins", tfwpgeomparam.a_ground_ins)

    monkeypatch.setattr(sctfcoil_module, "theta_coil", tfwpgeomparam.theta_coil)

    monkeypatch.setattr(sctfcoil_module, "tan_theta_coil", tfwpgeomparam.tan_theta_coil)

    sctfcoil.tf_wp_geom(i_tf_wp_geom=tfwpgeomparam.i_tf_wp_geom)

    assert tfcoil_variables.wwp1 == pytest.approx(tfwpgeomparam.expected_wwp1)

    assert sctfcoil_module.awpc == pytest.approx(tfwpgeomparam.expected_awpc)

    assert sctfcoil_module.awptf == pytest.approx(tfwpgeomparam.expected_awptf)

    assert sctfcoil_module.r_wp_inner == pytest.approx(
        tfwpgeomparam.expected_r_wp_inner
    )

    assert sctfcoil_module.r_wp_outer == pytest.approx(
        tfwpgeomparam.expected_r_wp_outer
    )

    assert sctfcoil_module.r_wp_centre == pytest.approx(
        tfwpgeomparam.expected_r_wp_centre
    )

    assert sctfcoil_module.t_wp_toroidal == pytest.approx(
        tfwpgeomparam.expected_t_wp_toroidal
    )

    assert sctfcoil_module.t_wp_toroidal_av == pytest.approx(
        tfwpgeomparam.expected_t_wp_toroidal_av
    )

    assert sctfcoil_module.a_ground_ins == pytest.approx(
        tfwpgeomparam.expected_a_ground_ins
    )


class TfCaseGeomParam(NamedTuple):
    acasetf: Any = None

    acasetfo: Any = None

    arealeg: Any = None

    tfareain: Any = None

    n_tf: Any = None

    casths: Any = None

    casthi: Any = None

    dr_tf_wp: Any = None

    r_tf_inboard_in: Any = None

    r_tf_inboard_out: Any = None

    awpc: Any = None

    r_wp_inner: Any = None

    r_wp_outer: Any = None

    t_lat_case_av: Any = None

    a_case_front: Any = None

    a_case_nose: Any = None

    theta_coil: Any = None

    tan_theta_coil: Any = None

    i_tf_wp_geom: Any = None

    i_tf_case_geom: Any = None

    expected_acasetf: Any = None

    expected_acasetfo: Any = None

    expected_t_lat_case_av: Any = None

    expected_a_case_front: Any = None

    expected_a_case_nose: Any = None


@pytest.mark.parametrize(
    "tfcasegeomparam",
    (
        TfCaseGeomParam(
            acasetf=0,
            acasetfo=0,
            arealeg=1.9805354702921749,
            tfareain=27.308689677971632,
            n_tf=16,
            casths=0.05000000000000001,
            casthi=0.060000000000000012,
            dr_tf_wp=0.54261087836601019,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            awpc=0.70527618095271016,
            r_wp_inner=3.5185911851091101,
            r_wp_outer=4.06120206347512,
            t_lat_case_av=0,
            a_case_front=0,
            a_case_nose=0,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=0,
            i_tf_case_geom=0,
            expected_acasetf=1.0015169239205168,
            expected_acasetfo=1.2752592893394648,
            expected_t_lat_case_av=0.10396600719086938,
            expected_a_case_front=0.18607458590131154,
            expected_a_case_nose=0.70261616505511615,
        ),
        TfCaseGeomParam(
            acasetf=1.0015169239205168,
            acasetfo=1.2752592893394648,
            arealeg=1.9805354702921749,
            tfareain=27.308689677971632,
            n_tf=16,
            casths=0.05000000000000001,
            casthi=0.060000000000000012,
            dr_tf_wp=0.54261087836601019,
            r_tf_inboard_in=2.9939411851091102,
            r_tf_inboard_out=4.20194118510911,
            awpc=0.70527618095271016,
            r_wp_inner=3.5185911851091101,
            r_wp_outer=4.06120206347512,
            t_lat_case_av=0.10396600719086938,
            a_case_front=0.18607458590131154,
            a_case_nose=0.70261616505511615,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            i_tf_wp_geom=0,
            i_tf_case_geom=0,
            expected_acasetf=1.0015169239205168,
            expected_acasetfo=1.2752592893394648,
            expected_t_lat_case_av=0.10396600719086938,
            expected_a_case_front=0.18607458590131154,
            expected_a_case_nose=0.70261616505511615,
        ),
    ),
)
def test_tf_case_geom(tfcasegeomparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_case_geom.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfcasegeomparam: the data used to mock and assert in this test.
    :type tfcasegeomparam: tfcasegeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(tfcoil_variables, "acasetf", tfcasegeomparam.acasetf)

    monkeypatch.setattr(tfcoil_variables, "acasetfo", tfcasegeomparam.acasetfo)

    monkeypatch.setattr(tfcoil_variables, "arealeg", tfcasegeomparam.arealeg)

    monkeypatch.setattr(tfcoil_variables, "tfareain", tfcasegeomparam.tfareain)

    monkeypatch.setattr(tfcoil_variables, "n_tf", tfcasegeomparam.n_tf)

    monkeypatch.setattr(tfcoil_variables, "casths", tfcasegeomparam.casths)

    monkeypatch.setattr(tfcoil_variables, "casthi", tfcasegeomparam.casthi)

    monkeypatch.setattr(tfcoil_variables, "dr_tf_wp", tfcasegeomparam.dr_tf_wp)

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_in", tfcasegeomparam.r_tf_inboard_in
    )

    monkeypatch.setattr(
        build_variables, "r_tf_inboard_out", tfcasegeomparam.r_tf_inboard_out
    )

    monkeypatch.setattr(sctfcoil_module, "awpc", tfcasegeomparam.awpc)

    monkeypatch.setattr(sctfcoil_module, "r_wp_inner", tfcasegeomparam.r_wp_inner)

    monkeypatch.setattr(sctfcoil_module, "r_wp_outer", tfcasegeomparam.r_wp_outer)

    monkeypatch.setattr(sctfcoil_module, "t_lat_case_av", tfcasegeomparam.t_lat_case_av)

    monkeypatch.setattr(sctfcoil_module, "a_case_front", tfcasegeomparam.a_case_front)

    monkeypatch.setattr(sctfcoil_module, "a_case_nose", tfcasegeomparam.a_case_nose)

    monkeypatch.setattr(sctfcoil_module, "theta_coil", tfcasegeomparam.theta_coil)

    monkeypatch.setattr(
        sctfcoil_module, "tan_theta_coil", tfcasegeomparam.tan_theta_coil
    )

    sctfcoil.tf_case_geom(
        i_tf_wp_geom=tfcasegeomparam.i_tf_wp_geom,
        i_tf_case_geom=tfcasegeomparam.i_tf_case_geom,
    )

    assert tfcoil_variables.acasetf == pytest.approx(tfcasegeomparam.expected_acasetf)

    assert tfcoil_variables.acasetfo == pytest.approx(tfcasegeomparam.expected_acasetfo)

    assert sctfcoil_module.t_lat_case_av == pytest.approx(
        tfcasegeomparam.expected_t_lat_case_av
    )

    assert sctfcoil_module.a_case_front == pytest.approx(
        tfcasegeomparam.expected_a_case_front
    )

    assert sctfcoil_module.a_case_nose == pytest.approx(
        tfcasegeomparam.expected_a_case_nose
    )


class TfIntegerTurnGeomParam(NamedTuple):
    dr_tf_wp: Any = None

    tinstf: Any = None

    tfinsgap: Any = None

    t_conductor: Any = None

    t_turn_tf: Any = None

    tfc_current: Any = None

    t_wp_toroidal: Any = None

    t_conductor_radial: Any = None

    t_conductor_toroidal: Any = None

    t_cable_radial: Any = None

    t_cable_toroidal: Any = None

    t_turn_radial: Any = None

    t_turn_toroidal: Any = None

    t_cable: Any = None

    n_layer: Any = None

    n_pancake: Any = None

    thwcndut: Any = None

    thicndut: Any = None

    expected_t_conductor: Any = None

    expected_t_turn_tf: Any = None

    expected_t_conductor_radial: Any = None

    expected_t_conductor_toroidal: Any = None

    expected_t_cable_radial: Any = None

    expected_t_cable_toroidal: Any = None

    expected_t_turn_radial: Any = None

    expected_t_turn_toroidal: Any = None

    expected_t_cable: Any = None

    expected_acstf: Any = None

    expected_acndttf: Any = None

    expected_insulation_area: Any = None

    expected_cpttf: Any = None

    expected_n_tf_turn: Any = None


@pytest.mark.parametrize(
    "tfintegerturngeomparam",
    (
        TfIntegerTurnGeomParam(
            dr_tf_wp=0.54261087836601019,
            tinstf=0.0080000000000000019,
            tfinsgap=0.01,
            t_conductor=0,
            t_turn_tf=0,
            tfc_current=14805350.287500001,
            t_wp_toroidal=1.299782604942499,
            t_conductor_radial=0,
            t_conductor_toroidal=0,
            t_cable_radial=0,
            t_cable_toroidal=0,
            t_turn_radial=0,
            t_turn_toroidal=0,
            t_cable=0,
            n_layer=10,
            n_pancake=20,
            thwcndut=0.0080000000000000002,
            thicndut=0.002,
            expected_t_conductor=0.052553108427885735,
            expected_t_turn_tf=0.056579413904423038,
            expected_t_conductor_radial=0.046661087836601015,
            expected_t_conductor_toroidal=0.059189130247124938,
            expected_t_cable_radial=0.030661087836601014,
            expected_t_cable_toroidal=0.043189130247124938,
            expected_t_turn_radial=0.050661087836601018,
            expected_t_turn_toroidal=0.063189130247124942,
            expected_t_cable=0.036389912284773368,
            expected_acstf=0.001293323051622732,
            expected_acndttf=0.0014685061538103825,
            expected_insulation_area=0.00043940087233490435,
            expected_cpttf=74026.751437500003,
            expected_n_tf_turn=200,
        ),
        TfIntegerTurnGeomParam(
            dr_tf_wp=0.54261087836601019,
            tinstf=0.0080000000000000019,
            tfinsgap=0.01,
            t_conductor=0.052553108427885735,
            t_turn_tf=0.056579413904423038,
            tfc_current=14805350.287500001,
            t_wp_toroidal=1.299782604942499,
            t_conductor_radial=0.046661087836601015,
            t_conductor_toroidal=0.059189130247124938,
            t_cable_radial=0.030661087836601014,
            t_cable_toroidal=0.043189130247124938,
            t_turn_radial=0.050661087836601018,
            t_turn_toroidal=0.063189130247124942,
            t_cable=0.036389912284773368,
            n_layer=10,
            n_pancake=20,
            thwcndut=0.0080000000000000002,
            thicndut=0.002,
            expected_t_conductor=0.052553108427885735,
            expected_t_turn_tf=0.056579413904423038,
            expected_t_conductor_radial=0.046661087836601015,
            expected_t_conductor_toroidal=0.059189130247124938,
            expected_t_cable_radial=0.030661087836601014,
            expected_t_cable_toroidal=0.043189130247124938,
            expected_t_turn_radial=0.050661087836601018,
            expected_t_turn_toroidal=0.063189130247124942,
            expected_t_cable=0.036389912284773368,
            expected_acstf=0.001293323051622732,
            expected_acndttf=0.0014685061538103825,
            expected_insulation_area=0.00043940087233490435,
            expected_cpttf=74026.751437500003,
            expected_n_tf_turn=200,
        ),
    ),
)
def test_tf_integer_turn_geom(tfintegerturngeomparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_integer_turn_geom.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param tfintegerturngeomparam: the data used to mock and assert in this test.
    :type tfintegerturngeomparam: tfintegerturngeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(tfcoil_variables, "dr_tf_wp", tfintegerturngeomparam.dr_tf_wp)

    monkeypatch.setattr(tfcoil_variables, "tinstf", tfintegerturngeomparam.tinstf)

    monkeypatch.setattr(tfcoil_variables, "tfinsgap", tfintegerturngeomparam.tfinsgap)

    monkeypatch.setattr(
        tfcoil_variables, "t_conductor", tfintegerturngeomparam.t_conductor
    )

    monkeypatch.setattr(tfcoil_variables, "t_turn_tf", tfintegerturngeomparam.t_turn_tf)

    monkeypatch.setattr(
        sctfcoil_module, "tfc_current", tfintegerturngeomparam.tfc_current
    )

    monkeypatch.setattr(
        sctfcoil_module, "t_wp_toroidal", tfintegerturngeomparam.t_wp_toroidal
    )

    monkeypatch.setattr(
        sctfcoil_module, "t_conductor_radial", tfintegerturngeomparam.t_conductor_radial
    )

    monkeypatch.setattr(
        sctfcoil_module,
        "t_conductor_toroidal",
        tfintegerturngeomparam.t_conductor_toroidal,
    )

    monkeypatch.setattr(
        sctfcoil_module, "t_cable_radial", tfintegerturngeomparam.t_cable_radial
    )

    monkeypatch.setattr(
        sctfcoil_module, "t_cable_toroidal", tfintegerturngeomparam.t_cable_toroidal
    )

    monkeypatch.setattr(
        sctfcoil_module, "t_turn_radial", tfintegerturngeomparam.t_turn_radial
    )

    monkeypatch.setattr(
        sctfcoil_module, "t_turn_toroidal", tfintegerturngeomparam.t_turn_toroidal
    )

    monkeypatch.setattr(sctfcoil_module, "t_cable", tfintegerturngeomparam.t_cable)

    (
        acstf,
        acndttf,
        insulation_area,
        cpttf,
        n_tf_turn,
    ) = sctfcoil.tf_integer_turn_geom(
        n_layer=tfintegerturngeomparam.n_layer,
        n_pancake=tfintegerturngeomparam.n_pancake,
        thwcndut=tfintegerturngeomparam.thwcndut,
        thicndut=tfintegerturngeomparam.thicndut,
    )

    assert tfcoil_variables.t_conductor == pytest.approx(
        tfintegerturngeomparam.expected_t_conductor
    )

    assert tfcoil_variables.t_turn_tf == pytest.approx(
        tfintegerturngeomparam.expected_t_turn_tf
    )

    assert sctfcoil_module.t_conductor_radial == pytest.approx(
        tfintegerturngeomparam.expected_t_conductor_radial
    )

    assert sctfcoil_module.t_conductor_toroidal == pytest.approx(
        tfintegerturngeomparam.expected_t_conductor_toroidal
    )

    assert sctfcoil_module.t_cable_radial == pytest.approx(
        tfintegerturngeomparam.expected_t_cable_radial
    )

    assert sctfcoil_module.t_cable_toroidal == pytest.approx(
        tfintegerturngeomparam.expected_t_cable_toroidal
    )

    assert sctfcoil_module.t_turn_radial == pytest.approx(
        tfintegerturngeomparam.expected_t_turn_radial
    )

    assert sctfcoil_module.t_turn_toroidal == pytest.approx(
        tfintegerturngeomparam.expected_t_turn_toroidal
    )

    assert sctfcoil_module.t_cable == pytest.approx(
        tfintegerturngeomparam.expected_t_cable
    )

    assert acstf == pytest.approx(tfintegerturngeomparam.expected_acstf)

    assert acndttf == pytest.approx(tfintegerturngeomparam.expected_acndttf)

    assert insulation_area == pytest.approx(
        tfintegerturngeomparam.expected_insulation_area
    )

    assert cpttf == pytest.approx(tfintegerturngeomparam.expected_cpttf)

    assert n_tf_turn == pytest.approx(tfintegerturngeomparam.expected_n_tf_turn)


class TfAveragedTurnGeomParam(NamedTuple):
    layer_ins: Any = None

    t_conductor: Any = None

    t_turn_tf: Any = None

    t_turn_tf_is_input: Any = None

    cpttf: Any = None

    t_cable_tf: Any = None

    t_cable_tf_is_input: Any = None

    awptf: Any = None

    t_turn_radial: Any = None

    t_turn_toroidal: Any = None

    t_cable: Any = None

    i_tf_sc_mat: Any = None

    jwptf: Any = None

    thwcndut: Any = None

    thicndut: Any = None

    expected_t_conductor: Any = None

    expected_t_turn_tf: Any = None

    expected_t_turn_radial: Any = None

    expected_t_turn_toroidal: Any = None

    expected_t_cable: Any = None

    expected_acstf: Any = None

    expected_acndttf: Any = None

    expected_insulation_area: Any = None

    expected_n_tf_turn: Any = None


@pytest.mark.parametrize(
    "tfaveragedturngeomparam",
    (
        TfAveragedTurnGeomParam(
            layer_ins=0,
            t_conductor=0,
            t_turn_tf=0,
            t_turn_tf_is_input=False,
            cpttf=65000,
            t_cable_tf=0,
            t_cable_tf_is_input=False,
            awptf=0.60510952642236249,
            t_turn_radial=0,
            t_turn_toroidal=0,
            t_cable=0,
            i_tf_sc_mat=5,
            jwptf=26493137.688284047,
            thwcndut=0.0080000000000000019,
            thicndut=0.00080000000000000004,
            expected_t_conductor=0.047932469413859431,
            expected_t_turn_tf=0.049532469413859428,
            expected_t_turn_radial=0.049532469413859428,
            expected_t_turn_toroidal=0.049532469413859428,
            expected_t_cable=0.031932469413859424,
            expected_acstf=0.00098877993839630008,
            expected_acndttf=0.0013087416857142699,
            expected_insulation_area=0.00015594390212434958,
            expected_n_tf_turn=246.63461538461544,
        ),
        TfAveragedTurnGeomParam(
            layer_ins=0,
            t_conductor=0.047932469413859431,
            t_turn_tf=0.049532469413859428,
            t_turn_tf_is_input=False,
            cpttf=65000,
            t_cable_tf=0,
            t_cable_tf_is_input=False,
            awptf=0.60510952642236249,
            t_turn_radial=0.049532469413859428,
            t_turn_toroidal=0.049532469413859428,
            t_cable=0.031932469413859424,
            i_tf_sc_mat=5,
            jwptf=26493137.688284047,
            thwcndut=0.0080000000000000019,
            thicndut=0.00080000000000000004,
            expected_t_conductor=0.047932469413859431,
            expected_t_turn_tf=0.049532469413859428,
            expected_t_turn_radial=0.049532469413859428,
            expected_t_turn_toroidal=0.049532469413859428,
            expected_t_cable=0.031932469413859424,
            expected_acstf=0.00098877993839630008,
            expected_acndttf=0.0013087416857142699,
            expected_insulation_area=0.00015594390212434958,
            expected_n_tf_turn=246.63461538461544,
        ),
        TfAveragedTurnGeomParam(
            layer_ins=0,
            t_conductor=5.712e-02,
            t_turn_tf=0.05872,
            t_turn_tf_is_input=True,
            cpttf=0,
            t_cable_tf=0,
            t_cable_tf_is_input=False,
            awptf=0.60510952642236249,
            t_turn_radial=0.05872,
            t_turn_toroidal=0.05872,
            t_cable=0.04109,
            i_tf_sc_mat=1,
            jwptf=2.301e07,
            thwcndut=8.015e-03,
            thicndut=8.0e-4,
            expected_t_conductor=5.712e-02,
            expected_t_turn_tf=0.05872,
            expected_t_turn_radial=0.05872,
            expected_t_turn_toroidal=0.05872,
            expected_t_cable=0.04109,
            expected_acstf=0.001657369442,
            expected_acndttf=0.001605324958,
            expected_insulation_area=0.000185344,
            expected_n_tf_turn=175.49384787,
        ),
        TfAveragedTurnGeomParam(
            layer_ins=0,
            t_conductor=0.058296,
            t_turn_tf=0,
            t_turn_tf_is_input=False,
            cpttf=0,
            t_cable_tf=0.042,
            t_cable_tf_is_input=True,
            awptf=0.60510952642236249,
            t_turn_radial=0.05872,
            t_turn_toroidal=0.05872,
            t_cable=0.04109,
            i_tf_sc_mat=1,
            jwptf=2.673e07,
            thwcndut=8.148e-03,
            thicndut=8.0e-4,
            expected_t_conductor=0.058296,
            expected_t_turn_tf=0.059896,
            expected_t_turn_radial=0.059896,
            expected_t_turn_toroidal=0.059896,
            expected_t_cable=0.042,
            expected_acstf=0.001731943361,
            expected_acndttf=0.001666480255,
            expected_insulation_area=0.00018910719999999962,
            expected_n_tf_turn=168.6701961481806,
        ),
    ),
)
def test_tf_averaged_turn_geom(tfaveragedturngeomparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_averaged_turn_geom.

    This test was generated using data from tests/regression/scenarios/i_mode/IN.DAT.

    :param tfaveragedturngeomparam: the data used to mock and assert in this test.
    :type tfaveragedturngeomparam: tfaveragedturngeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(
        tfcoil_variables, "layer_ins", tfaveragedturngeomparam.layer_ins
    )

    monkeypatch.setattr(
        tfcoil_variables, "t_conductor", tfaveragedturngeomparam.t_conductor
    )

    monkeypatch.setattr(
        tfcoil_variables, "t_turn_tf", tfaveragedturngeomparam.t_turn_tf
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "t_turn_tf_is_input",
        tfaveragedturngeomparam.t_turn_tf_is_input,
    )

    monkeypatch.setattr(tfcoil_variables, "cpttf", tfaveragedturngeomparam.cpttf)

    monkeypatch.setattr(
        tfcoil_variables, "t_cable_tf", tfaveragedturngeomparam.t_cable_tf
    )

    monkeypatch.setattr(
        tfcoil_variables,
        "t_cable_tf_is_input",
        tfaveragedturngeomparam.t_cable_tf_is_input,
    )

    monkeypatch.setattr(sctfcoil_module, "awptf", tfaveragedturngeomparam.awptf)

    monkeypatch.setattr(
        sctfcoil_module, "t_turn_radial", tfaveragedturngeomparam.t_turn_radial
    )

    monkeypatch.setattr(
        sctfcoil_module, "t_turn_toroidal", tfaveragedturngeomparam.t_turn_toroidal
    )

    monkeypatch.setattr(sctfcoil_module, "t_cable", tfaveragedturngeomparam.t_cable)

    acstf, acndttf, insulation_area, n_tf_turn = sctfcoil.tf_averaged_turn_geom(
        i_tf_sc_mat=tfaveragedturngeomparam.i_tf_sc_mat,
        jwptf=tfaveragedturngeomparam.jwptf,
        thwcndut=tfaveragedturngeomparam.thwcndut,
        thicndut=tfaveragedturngeomparam.thicndut,
    )

    assert tfcoil_variables.t_conductor == pytest.approx(
        tfaveragedturngeomparam.expected_t_conductor
    )

    assert tfcoil_variables.t_turn_tf == pytest.approx(
        tfaveragedturngeomparam.expected_t_turn_tf
    )

    assert sctfcoil_module.t_turn_radial == pytest.approx(
        tfaveragedturngeomparam.expected_t_turn_radial
    )

    assert sctfcoil_module.t_turn_toroidal == pytest.approx(
        tfaveragedturngeomparam.expected_t_turn_toroidal
    )

    assert sctfcoil_module.t_cable == pytest.approx(
        tfaveragedturngeomparam.expected_t_cable
    )

    assert acstf == pytest.approx(tfaveragedturngeomparam.expected_acstf)

    assert acndttf == pytest.approx(tfaveragedturngeomparam.expected_acndttf)

    assert insulation_area == pytest.approx(
        tfaveragedturngeomparam.expected_insulation_area
    )

    assert n_tf_turn == pytest.approx(tfaveragedturngeomparam.expected_n_tf_turn)


class TfWpCurrentsParam(NamedTuple):
    ritfc: Any = None

    n_tf: Any = None

    jwptf: Any = None

    awptf: Any = None

    expected_jwptf: Any = None


@pytest.mark.parametrize(
    "tfwpcurrentsparam",
    (
        TfWpCurrentsParam(
            ritfc=256500000.00000003,
            n_tf=16,
            jwptf=0,
            awptf=0.60510952642236249,
            expected_jwptf=26493137.688284047,
        ),
        TfWpCurrentsParam(
            ritfc=256500000.00000003,
            n_tf=16,
            jwptf=26493137.688284047,
            awptf=0.60510952642236249,
            expected_jwptf=26493137.688284047,
        ),
    ),
)
def test_tf_wp_currents(tfwpcurrentsparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for tf_wp_currents.

    This test was generated using data from tests/regression/scenarios/i_mode/IN.DAT.

    :param tfwpcurrentsparam: the data used to mock and assert in this test.
    :type tfwpcurrentsparam: tfwpcurrentsparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(tfcoil_variables, "ritfc", tfwpcurrentsparam.ritfc)

    monkeypatch.setattr(tfcoil_variables, "n_tf", tfwpcurrentsparam.n_tf)

    monkeypatch.setattr(tfcoil_variables, "jwptf", tfwpcurrentsparam.jwptf)

    monkeypatch.setattr(sctfcoil_module, "awptf", tfwpcurrentsparam.awptf)

    sctfcoil.tf_wp_currents()

    assert tfcoil_variables.jwptf == pytest.approx(tfwpcurrentsparam.expected_jwptf)


class StressclParam(NamedTuple):
    tfcth: Any = None

    r_tf_inboard_mid: Any = None

    bore: Any = None

    ohcth: Any = None

    tf_in_cs: Any = None

    gapoh: Any = None

    hmax: Any = None

    r_tf_inboard_in: Any = None

    casestr: Any = None

    n_tf_turn: Any = None

    dr_tf_wp: Any = None

    i_tf_tresca: Any = None

    acstf: Any = None

    vforce: Any = None

    ritfc: Any = None

    jwptf: Any = None

    sig_tf_cs_bucked: Any = None

    sig_tf_case: Any = None

    sig_tf_wp: Any = None

    thwcndut: Any = None

    insstrain: Any = None

    tinstf: Any = None

    thicndut: Any = None

    acndttf: Any = None

    tfinsgap: Any = None

    acasetf: Any = None

    sig_tf_case_max: Any = None

    poisson_steel: Any = None

    poisson_copper: Any = None

    poisson_al: Any = None

    n_tf_graded_layers: Any = None

    i_tf_sup: Any = None

    i_tf_bucking: Any = None

    fcoolcp: Any = None

    eyoung_cond_axial: Any = None

    eyoung_steel: Any = None

    eyoung_res_tf_buck: Any = None

    eyoung_ins: Any = None

    eyoung_al: Any = None

    eyoung_copper: Any = None

    aiwp: Any = None

    aswp: Any = None

    cpttf: Any = None

    n_tf: Any = None

    i_tf_stress_model: Any = None

    sig_tf_wp_max: Any = None

    i_tf_turns_integer: Any = None

    casthi: Any = None

    acond: Any = None

    avwp: Any = None

    awphec: Any = None

    poisson_ins: Any = None

    eyoung_cond_trans: Any = None

    poisson_cond_axial: Any = None

    poisson_cond_trans: Any = None

    dhecoil: Any = None

    fcutfsu: Any = None

    str_wp: Any = None

    n_tf_wp_layers: Any = None

    ipfres: Any = None

    oh_steel_frac: Any = None

    ohhghf: Any = None

    coheof: Any = None

    cohbop: Any = None

    ncls: Any = None

    cptdin: Any = None

    awpc: Any = None

    a_tf_steel: Any = None

    a_tf_ins: Any = None

    r_wp_inner: Any = None

    r_wp_outer: Any = None

    t_wp_toroidal: Any = None

    t_wp_toroidal_av: Any = None

    t_lat_case_av: Any = None

    a_case_front: Any = None

    a_case_nose: Any = None

    theta_coil: Any = None

    tan_theta_coil: Any = None

    t_cable_radial: Any = None

    t_cable: Any = None

    vforce_inboard_tot: Any = None

    iprint: Any = None

    outfile: Any = None

    n_radial_array: Any = None

    n_tf_layer: Any = None

    expected_casestr: Any = None

    expected_sig_tf_case: Any = None

    expected_sig_tf_wp: Any = None

    expected_insstrain: Any = None

    expected_str_wp: Any = None


@pytest.mark.parametrize(
    "stressclparam",
    (
        StressclParam(
            tfcth=1.208,
            tf_in_cs=0,
            gapoh=0.01,
            r_tf_inboard_mid=3.5979411851091103,
            bore=2.3322000000000003,
            ohcth=0.55242000000000002,
            hmax=9.0730900215620327,
            r_tf_inboard_in=2.9939411851091102,
            casestr=0,
            n_tf_turn=200,
            dr_tf_wp=0.54261087836601019,
            i_tf_tresca=0,
            acstf=0.001293323051622732,
            vforce=250545611.13801825,
            ritfc=236885604.60000002,
            jwptf=23124470.793774806,
            sig_tf_cs_bucked=0,
            sig_tf_case=0,
            sig_tf_wp=0,
            thwcndut=0.0080000000000000002,
            insstrain=0,
            tinstf=0.0080000000000000019,
            thicndut=0.002,
            acndttf=0.0014685061538103825,
            tfinsgap=0.01,
            acasetf=1.0015169239205168,
            sig_tf_case_max=580000000,
            poisson_steel=0.29999999999999999,
            poisson_copper=0.34999999999999998,
            poisson_al=0.34999999999999998,
            n_tf_graded_layers=1,
            i_tf_sup=1,
            i_tf_bucking=1,
            fcoolcp=0.29999999999999999,
            eyoung_cond_axial=0,
            eyoung_steel=205000000000,
            eyoung_res_tf_buck=150000000000,
            eyoung_ins=20000000000,
            eyoung_al=numpy.array(69000000000.0),
            eyoung_copper=numpy.array(117000000000.0),
            aiwp=0.087880174466980876,
            aswp=0.29370123076207649,
            cpttf=74026.751437500003,
            n_tf=16,
            i_tf_stress_model=1,
            sig_tf_wp_max=580000000,
            i_tf_turns_integer=1,
            casthi=0.060000000000000012,
            acond=0.1653572639592335,
            avwp=0.07759938309736393,
            awphec=0.015707963267948974,
            poisson_ins=0.34000000000000002,
            eyoung_cond_trans=0,
            poisson_cond_axial=0.30000001192092896,
            poisson_cond_trans=0.30000001192092896,
            dhecoil=0.010000000000000002,
            fcutfsu=0.80884,
            str_wp=0,
            n_tf_wp_layers=5,
            ipfres=0,
            oh_steel_frac=0.57874999999999999,
            ohhghf=0.90000000000000002,
            coheof=20726000,
            cohbop=0,
            ncls=numpy.array(
                numpy.array((1, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            cptdin=numpy.array(
                numpy.array(
                    (
                        42200,
                        42200,
                        42200,
                        42200,
                        43000,
                        43000,
                        43000,
                        43000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            awpc=0.70527618095271016,
            a_tf_steel=1.2952181546825934,
            a_tf_ins=0.11646247019991701,
            r_wp_inner=3.5185911851091101,
            r_wp_outer=4.06120206347512,
            t_wp_toroidal=1.299782604942499,
            t_wp_toroidal_av=1.299782604942499,
            t_lat_case_av=0.10396600719086938,
            a_case_front=0.18607458590131154,
            a_case_nose=0.70261616505511615,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            t_cable_radial=0.030661087836601014,
            t_cable=0.036389912284773368,
            vforce_inboard_tot=4008729778.208292,
            iprint=0,
            outfile=11,
            n_radial_array=100,
            n_tf_layer=3,
            expected_casestr=0.00094360452596334093,
            expected_sig_tf_case=543381805.25001633,
            expected_sig_tf_wp=397005702.35272157,
            expected_insstrain=-0.006687152422925652,
            expected_str_wp=0.0015619754370069119,
        ),
        StressclParam(
            tfcth=1.208,
            tf_in_cs=0,
            gapoh=0.01,
            r_tf_inboard_mid=3.5979411851091103,
            bore=2.3322000000000003,
            ohcth=0.55242000000000002,
            hmax=9.0730900215620327,
            r_tf_inboard_in=2.9939411851091102,
            casestr=0.00094360452596334093,
            n_tf_turn=200,
            dr_tf_wp=0.54261087836601019,
            i_tf_tresca=0,
            acstf=0.001293323051622732,
            vforce=250545611.13801825,
            ritfc=236885604.60000002,
            jwptf=23124470.793774806,
            sig_tf_cs_bucked=0,
            sig_tf_case=543381805.25001633,
            sig_tf_wp=397005702.35272157,
            thwcndut=0.0080000000000000002,
            insstrain=0,
            tinstf=0.0080000000000000019,
            thicndut=0.002,
            acndttf=0.0014685061538103825,
            tfinsgap=0.01,
            acasetf=1.0015169239205168,
            sig_tf_case_max=580000000,
            poisson_steel=0.29999999999999999,
            poisson_copper=0.34999999999999998,
            poisson_al=0.34999999999999998,
            n_tf_graded_layers=1,
            i_tf_sup=1,
            i_tf_bucking=1,
            fcoolcp=0.29999999999999999,
            eyoung_cond_axial=0,
            eyoung_steel=205000000000,
            eyoung_res_tf_buck=150000000000,
            eyoung_ins=20000000000,
            eyoung_al=numpy.array(69000000000.0),
            eyoung_copper=numpy.array(117000000000.0),
            aiwp=0.087880174466980876,
            aswp=0.29370123076207649,
            cpttf=74026.751437500003,
            n_tf=16,
            i_tf_stress_model=1,
            sig_tf_wp_max=580000000,
            i_tf_turns_integer=1,
            casthi=0.060000000000000012,
            acond=0.1653572639592335,
            avwp=0.07759938309736393,
            awphec=0.015707963267948974,
            poisson_ins=0.34000000000000002,
            eyoung_cond_trans=0,
            poisson_cond_axial=0.30000001192092896,
            poisson_cond_trans=0.30000001192092896,
            dhecoil=0.010000000000000002,
            fcutfsu=0.80884,
            str_wp=0.0015619754370069119,
            n_tf_wp_layers=5,
            ipfres=0,
            oh_steel_frac=0.57874999999999999,
            ohhghf=0.90000000000000002,
            coheof=20726000,
            cohbop=19311657.760000002,
            ncls=numpy.array(
                numpy.array((1, 1, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            cptdin=numpy.array(
                numpy.array(
                    (
                        42200,
                        42200,
                        42200,
                        42200,
                        43000,
                        43000,
                        43000,
                        43000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                        40000,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            awpc=0.70527618095271016,
            a_tf_steel=1.2952181546825934,
            a_tf_ins=0.11646247019991701,
            r_wp_inner=3.5185911851091101,
            r_wp_outer=4.06120206347512,
            t_wp_toroidal=1.299782604942499,
            t_wp_toroidal_av=1.299782604942499,
            t_lat_case_av=0.10396600719086938,
            a_case_front=0.18607458590131154,
            a_case_nose=0.70261616505511615,
            theta_coil=0.19634954084936207,
            tan_theta_coil=0.19891236737965801,
            t_cable_radial=0.030661087836601014,
            t_cable=0.036389912284773368,
            vforce_inboard_tot=4008729778.208292,
            iprint=0,
            outfile=11,
            n_radial_array=100,
            n_tf_layer=3,
            expected_casestr=0.00094360452596334093,
            expected_sig_tf_case=543381805.25001633,
            expected_sig_tf_wp=397005702.35272157,
            expected_insstrain=-0.006687152422925652,
            expected_str_wp=0.0015619754370069119,
        ),
    ),
)
def test_stresscl(stressclparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for stresscl.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param stressclparam: the data used to mock and assert in this test.
    :type stressclparam: stressclparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    (
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        _,
        sig_tf_wp,
        sig_tf_case,
        _,
        str_wp,
        casestr,
        insstrain,
        _,
    ) = sctfcoil.stresscl(
        stressclparam.n_tf_layer,
        stressclparam.n_radial_array,
        stressclparam.n_tf_wp_layers,
        stressclparam.i_tf_bucking,
        stressclparam.r_tf_inboard_in,
        stressclparam.bore,
        stressclparam.hmax,
        stressclparam.ohhghf,
        stressclparam.ohcth,
        stressclparam.tf_in_cs,
        stressclparam.tfcth,
        stressclparam.gapoh,
        stressclparam.ipfres,
        stressclparam.coheof,
        stressclparam.cohbop,
        stressclparam.cptdin,
        stressclparam.ncls,
        70 / 22,
        3e-3,
        stressclparam.oh_steel_frac,
        stressclparam.eyoung_steel,
        stressclparam.poisson_steel,
        stressclparam.eyoung_cond_axial,
        stressclparam.poisson_cond_axial,
        stressclparam.eyoung_cond_trans,
        stressclparam.poisson_cond_trans,
        stressclparam.eyoung_ins,
        stressclparam.poisson_ins,
        stressclparam.thicndut,
        stressclparam.eyoung_copper,
        stressclparam.poisson_copper,
        stressclparam.i_tf_sup,
        stressclparam.eyoung_res_tf_buck,
        stressclparam.r_wp_inner,
        stressclparam.tan_theta_coil,
        stressclparam.theta_coil,
        stressclparam.r_wp_outer,
        stressclparam.a_tf_steel,
        stressclparam.a_case_front,
        stressclparam.a_case_nose,
        stressclparam.tfinsgap,
        stressclparam.tinstf,
        stressclparam.n_tf_turn,
        stressclparam.i_tf_turns_integer,
        stressclparam.t_cable,
        stressclparam.t_cable_radial,
        stressclparam.dhecoil,
        stressclparam.fcutfsu,
        stressclparam.thwcndut,
        stressclparam.t_lat_case_av,
        stressclparam.t_wp_toroidal_av,
        stressclparam.a_tf_ins,
        stressclparam.aswp,
        stressclparam.acond,
        stressclparam.awpc,
        stressclparam.eyoung_al,
        stressclparam.poisson_al,
        stressclparam.fcoolcp,
        stressclparam.n_tf_graded_layers,
        stressclparam.ritfc,
        stressclparam.casthi,
        stressclparam.i_tf_stress_model,
        stressclparam.vforce_inboard_tot,
        stressclparam.i_tf_tresca,
        stressclparam.acasetf,
        stressclparam.vforce,
        stressclparam.acndttf,
    )

    assert casestr == pytest.approx(stressclparam.expected_casestr, rel=0.01)

    assert sig_tf_case == pytest.approx(stressclparam.expected_sig_tf_case, rel=0.01)

    assert sig_tf_wp == pytest.approx(stressclparam.expected_sig_tf_wp, rel=0.01)

    assert insstrain == pytest.approx(stressclparam.expected_insstrain, rel=0.01)

    assert str_wp == pytest.approx(stressclparam.expected_str_wp, rel=0.01)


class PlaneStressParam(NamedTuple):
    linesolv: Any = None

    n_radial_array: Any = None

    nlayers: Any = None

    nu: Any = None

    rad: Any = None

    ey: Any = None

    j: Any = None

    expected_sigr: Any = None

    expected_sigt: Any = None

    expected_r_deflect: Any = None

    expected_rradius: Any = None


@pytest.mark.parametrize(
    "planestressparam",
    (
        PlaneStressParam(
            n_radial_array=100,
            nlayers=3,
            nu=numpy.array(
                numpy.array(
                    (0.29999999999999999, 0.30904421667064924, 0.29999999999999999),
                    order="F",
                ),
                order="F",
            ).transpose(),
            rad=numpy.array(
                numpy.array(
                    (
                        2.9939411851091102,
                        3.5414797139565706,
                        4.0876202904571599,
                        4.1476202904571595,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            ey=numpy.array(
                numpy.array(
                    (205000000000, 43126670035.025253, 205000000000), order="F"
                ),
                order="F",
            ).transpose(),
            j=numpy.array(
                numpy.array((0, 18097185.781970859, 0), order="F"), order="F"
            ).transpose(),
            expected_sigr=[
                0.0,
                -638231.855766482,
                -1272978.0149700674,
                -1904263.8142719064,
                -2532114.3605434406,
                -3156554.533362813,
                -3777608.9874793817,
                -4395302.155247777,
                -5009658.2490304215,
                -5620701.263570679,
                -6228454.97833536,
                -6832942.959828614,
                -7434188.563876412,
                -8032214.937882571,
                -8627045.023056716,
                -9218701.556614282,
                -9807207.07394974,
                -10392583.910782054,
                -10974854.205274347,
                -11554039.900126804,
                -12130162.744643984,
                -12703244.296776803,
                -13273305.925139055,
                -13840368.810999524,
                -14404453.950249525,
                -14965582.155346265,
                -15523774.05723289,
                -16079050.107234456,
                -16631430.578931632,
                -17180935.57001094,
                -17727585.00409306,
                -18271398.632538624,
                -18812396.036232498,
                -19350596.627346307,
                -19886019.651079975,
                -20418684.18738201,
                -20948609.1526495,
                -21475813.301407494,
                -22000315.22796835,
                -22522133.368071266,
                -23041286.00050247,
                -23557791.24869583,
                -24071667.082314458,
                -24582931.31881398,
                -25091601.624986775,
                -25597695.518488254,
                -26101230.36934514,
                -26602223.401445802,
                -27100691.6940132,
                -27596652.183060553,
                -28090121.6628299,
                -28581116.78721381,
                -29069654.071160544,
                -29555749.892062735,
                -30039420.491130088,
                -30520681.974745892,
                -30999550.31580798,
                -31476041.355053987,
                -31950170.802371547,
                -32421954.23809305,
                -32891407.114275593,
                -33358544.755966377,
                -33823382.3624533,
                -34285935.00850141,
                -34746217.64557501,
                -35204245.10304597,
                -35660032.08938813,
                -36113593.19335821,
                -36564942.88516308,
                -37014095.51761406,
                -37461065.327267796,
                -37905866.4355546,
                -38348512.84989362,
                -38789018.46479585,
                -39227397.06295431,
                -39663662.31632225,
                -40097827.78717919,
                -40529906.92918485,
                -40959913.08842154,
                -41387859.50442459,
                -41813759.31120163,
                -42237625.5382402,
                -42659471.11150431,
                -43079308.85441998,
                -43497151.48884965,
                -43913011.6360561,
                -44326901.817655385,
                -44738834.45655951,
                -45148821.87790862,
                -45556876.30999292,
                -45963009.88516445,
                -46367234.64073896,
                -46769562.51988788,
                -47170005.372520424,
                -47568574.956156105,
                -47965282.936788,
                -48360140.88973615,
                -48753160.30049217,
                -49144352.56555428,
                -49533728.99325363,
                -49921300.804571584,
                -49954359.75079646,
                -49974994.421951994,
                -49983237.62599215,
                -49979121.95724278,
                -49962679.79819553,
                -49933943.321276575,
                -49892944.49059213,
                -49839715.063674085,
                -49774286.59318251,
                -49696690.42861684,
                -49606957.71798797,
                -49505119.40949372,
                -49391206.253157824,
                -49265248.802470304,
                -49127277.41599873,
                -48977322.25898817,
                -48815413.304947995,
                -48641580.33721941,
                -48455852.950529814,
                -48258260.55252331,
                -48048832.365294494,
                -47827597.42688793,
                -47594584.59279222,
                -47349822.53742223,
                -47093339.75557484,
                -46825164.56388664,
                -46545325.10226525,
                -46253849.335310586,
                -45950765.05372232,
                -45636099.875691466,
                -45309881.24828925,
                -44972136.44881964,
                -44622892.586189084,
                -44262176.60223319,
                -43890015.27305536,
                -43506435.21033433,
                -43111462.86262834,
                -42705124.516666815,
                -42287446.29862875,
                -41858454.17540391,
                -41418173.95584675,
                -40966631.29202375,
                -40503851.68043686,
                -40029860.46324105,
                -39544682.82945457,
                -39048343.816151306,
                -38540868.309644595,
                -38022281.04666105,
                -37492606.61550055,
                -36951869.45718825,
                -36400093.8666107,
                -35837303.99364805,
                -35263523.8442963,
                -34678777.28176432,
                -34083088.027581,
                -33476479.662681654,
                -32858975.628480688,
                -32230599.22794289,
                -31591373.626639504,
                -30941321.853795685,
                -30280466.803329386,
                -29608831.23487964,
                -28926437.774824906,
                -28233308.91729218,
                -27529467.025158796,
                -26814934.33104299,
                -26089732.938284557,
                -25353884.821918193,
                -24607411.829640254,
                -23850335.682763584,
                -23082677.977158453,
                -22304460.184201445,
                -21515703.651693877,
                -20716429.60478688,
                -19906659.14690248,
                -19086413.26062489,
                -18255712.808604427,
                -17414578.534444164,
                -16563031.063583972,
                -15701090.904165361,
                -14828778.44790294,
                -13946113.970938675,
                -13053117.634692835,
                -12149809.48670104,
                -11236209.461454624,
                -10312337.381225117,
                -9378212.956880113,
                -8433855.788701806,
                -7479285.367186312,
                -6514521.073844995,
                -5539582.181995865,
                -4554487.85754239,
                -3559257.1597569585,
                -2553909.0420464603,
                -1538462.3527192033,
                -512935.83573593997,
                522651.86853248376,
                1568282.2225666705,
                2623936.7909615375,
                3689597.239691208,
                3651895.231238659,
                3614209.818926577,
                3576540.9930157103,
                3538888.743774036,
                3501253.061476492,
                3463633.9364053686,
                3426031.358849868,
                3388445.3191064713,
                3350875.807478765,
                3313322.8142772503,
                3275786.3298197784,
                3238266.3444311135,
                3200762.8484432492,
                3163275.832195042,
                3125805.2860326516,
                3088351.200309296,
                3050913.5653850324,
                3013492.3716272456,
                2976087.6094101095,
                2938699.2691151258,
                2901327.3411305123,
                2863971.8158518155,
                2826632.683681444,
                2789309.93502884,
                2752003.560310531,
                2714713.5499499296,
                2677439.8943775576,
                2640182.5840308485,
                2602941.609354294,
                2565716.9607993956,
                2528508.6288244454,
                2491316.6038949895,
                2454140.8764832406,
                2416981.4370685695,
                2379838.2761372332,
                2342711.3841824015,
                2305600.7517043273,
                2268506.36920998,
                2231428.2272134367,
                2194366.316235539,
                2157320.6268042377,
                2120291.149454175,
                2083277.8747270512,
                2046280.7931714796,
                2009299.8953427651,
                1972335.1718033939,
                1935386.6131224218,
                1898454.2098760372,
                1861537.9526470467,
                1824637.8320254136,
                1787753.838607672,
                1750885.962997364,
                1714034.1958049217,
                1677198.5276473942,
                1640378.949148915,
                1603575.450940211,
                1566788.0236590446,
                1530016.6579497962,
                1493261.3444638075,
                1456522.0738591857,
                1419798.8368006812,
                1383091.6239601525,
                1346400.4260158325,
                1309725.23365311,
                1273066.0375638467,
                1236422.8284469144,
                1199795.5970078032,
                1163184.33395872,
                1126589.0300187848,
                1090009.6759136862,
                1053446.262375977,
                1016898.7801448042,
                980367.2199661784,
                943851.5725928267,
                907351.8287839973,
                870867.979305924,
                834400.0149312643,
                797947.9264396608,
                761511.704617082,
                725091.3402565308,
                688686.8241575557,
                652298.1471262753,
                615925.2999756474,
                579568.2735251015,
                543227.0586009554,
                506901.646035926,
                470592.02666951995,
                434298.1913479606,
                398020.1309238218,
                361757.8362566382,
                325511.29821222165,
                289280.5076633691,
                253065.4554891793,
                216866.13257544299,
                180682.52981469216,
                144514.63810580914,
                108362.44835451504,
                72225.95147280802,
                36105.13837954948,
            ],
            expected_sigt=[
                -349942877.4275314,
                -349304645.57176495,
                -348669899.41256136,
                -348038613.61325955,
                -347410763.066988,
                -346786322.8941686,
                -346165268.44005203,
                -345547575.2722837,
                -344933219.178501,
                -344322176.1639607,
                -343714422.44919604,
                -343109934.4677028,
                -342508688.863655,
                -341910662.4896488,
                -341315832.4044747,
                -340724175.87091714,
                -340135670.35358167,
                -339550293.5167493,
                -338968023.222257,
                -338388837.5274046,
                -337812714.68288743,
                -337239633.1307546,
                -336669571.5023923,
                -336102508.61653185,
                -335538423.4772818,
                -334977295.2721851,
                -334419103.37029845,
                -333863827.3202969,
                -333311446.84859973,
                -332761941.85752034,
                -332215292.42343825,
                -331671478.7949927,
                -331130481.39129883,
                -330592280.800185,
                -330056857.7764513,
                -329524193.2401493,
                -328994268.2748818,
                -328467064.1261238,
                -327942562.1995629,
                -327420744.05946004,
                -326901591.42702883,
                -326385086.17883545,
                -325871210.3452168,
                -325359946.10871726,
                -324851275.8025445,
                -324345181.90904295,
                -323841647.0581861,
                -323340654.02608544,
                -322842185.73351806,
                -322346225.24447066,
                -321852755.7647013,
                -321361760.64031744,
                -320873223.3563707,
                -320387127.5354685,
                -319903456.9364011,
                -319422195.4527854,
                -318943327.11172324,
                -318466836.0724772,
                -317992706.6251597,
                -317520923.18943816,
                -317051470.3132556,
                -316584332.6715648,
                -316119495.0650779,
                -315656942.4190298,
                -315196659.7819562,
                -314738632.32448524,
                -314282845.33814305,
                -313829284.23417294,
                -313377934.54236805,
                -312928781.9099172,
                -312481812.10026336,
                -312037010.99197656,
                -311594364.57763755,
                -311153858.9627353,
                -310715480.3645769,
                -310279215.1112089,
                -309845049.640352,
                -309412970.4983463,
                -308982964.3391096,
                -308555017.92310655,
                -308129118.11632955,
                -307705251.8892909,
                -307283406.3160268,
                -306863568.5731112,
                -306445725.9386815,
                -306029865.791475,
                -305615975.60987574,
                -305204042.97097164,
                -304794055.5496225,
                -304386001.1175382,
                -303979867.5423667,
                -303575642.78679216,
                -303173314.9076432,
                -302772872.0550107,
                -302374302.471375,
                -301977594.4907431,
                -301582736.53779495,
                -301189717.12703896,
                -300798524.8619768,
                -300409148.4342775,
                -75393985.39087032,
                -75352895.95330013,
                -75308186.27146304,
                -75259848.14709243,
                -75207873.48230305,
                -75152254.27849749,
                -75092982.635284,
                -75030050.74940914,
                -74963450.91370073,
                -74893175.51602733,
                -74819217.03826267,
                -74741568.05526796,
                -74660221.23388477,
                -74575169.33193818,
                -74486405.1972532,
                -74393921.76667804,
                -74297712.06512664,
                -74197769.20462456,
                -74094086.38337062,
                -73986656.88480622,
                -73875474.07669795,
                -73760531.41022801,
                -73641822.41909933,
                -73519340.71864453,
                -73393080.00495018,
                -73263034.05399033,
                -73129196.72076792,
                -72991561.93846779,
                -72850123.71761769,
                -72704876.14526016,
                -72555813.3841336,
                -72402929.67186344,
                -72246219.32015613,
                -72085676.71401668,
                -71921296.31095727,
                -71753072.64022708,
                -71581000.30204895,
                -71405073.9668568,
                -71225288.37455449,
                -71041638.33377115,
                -70854118.72113262,
                -70662724.48053451,
                -70467450.62243325,
                -70268292.22313243,
                -70065244.4240877,
                -69858302.43121327,
                -69647461.51419875,
                -69432717.00583516,
                -69214064.30134167,
                -68991498.85770856,
                -68765016.19304176,
                -68534611.8859144,
                -68300281.57473324,
                -68062020.95709684,
                -67819825.78918111,
                -67573691.88511074,
                -67323615.11635438,
                -67069591.41111759,
                -66811616.75374359,
                -66549687.18412141,
                -66283798.797103256,
                -66013947.74192255,
                -65740130.22162387,
                -65462342.49249586,
                -65180580.86350979,
                -64894841.695768714,
                -64605121.401957214,
                -64311416.44580017,
                -64013723.341527745,
                -63712038.653342426,
                -63406358.99489821,
                -63096681.02878188,
                -62783001.46599489,
                -62465317.06545305,
                -62143624.63347747,
                -61817921.02330379,
                -61488203.13458894,
                -61154467.91292381,
                -60816712.34935414,
                -60474933.47990544,
                -60129128.38510988,
                -59779294.189546056,
                -59425428.061372176,
                -59067527.21187808,
                -58705588.895027705,
                -58339610.40701998,
                -57969589.08584027,
                -57595522.31083187,
                -57217407.502260946,
                -56835242.12088977,
                -56449023.667553514,
                -56058749.68274543,
                -55664417.74619955,
                -55266025.4764855,
                -54863570.530602135,
                -54457050.60357594,
                -54046463.4280663,
                -53631806.773974754,
                -53213078.44805462,
                -52790276.293530844,
                -253219448.39249134,
                -253181746.38403878,
                -253144060.97172672,
                -253106392.14581582,
                -253068739.89657417,
                -253031104.21427664,
                -252993485.0892055,
                -252955882.51165,
                -252918296.4719066,
                -252880726.9602789,
                -252843173.96707734,
                -252805637.48261994,
                -252768117.49723122,
                -252730614.00124338,
                -252693126.98499516,
                -252655656.43883282,
                -252618202.35310942,
                -252580764.71818522,
                -252543343.52442738,
                -252505938.76221025,
                -252468550.42191526,
                -252431178.4939306,
                -252393822.96865195,
                -252356483.8364816,
                -252319161.08782896,
                -252281854.7131107,
                -252244564.70275006,
                -252207291.04717773,
                -252170033.736831,
                -252132792.76215443,
                -252095568.1135995,
                -252058359.7816246,
                -252021167.75669512,
                -251983992.0292834,
                -251946832.5898687,
                -251909689.42893735,
                -251872562.53698257,
                -251835451.90450448,
                -251798357.5220101,
                -251761279.38001359,
                -251724217.46903566,
                -251687171.77960438,
                -251650142.3022543,
                -251613129.02752718,
                -251576131.94597158,
                -251539151.04814288,
                -251502186.32460356,
                -251465237.76592252,
                -251428305.36267614,
                -251391389.1054472,
                -251354488.98482555,
                -251317604.99140784,
                -251280737.1157975,
                -251243885.34860504,
                -251207049.68044755,
                -251170230.10194904,
                -251133426.60374033,
                -251096639.1764592,
                -251059867.81074992,
                -251023112.4972639,
                -250986373.22665933,
                -250949649.98960084,
                -250912942.77676025,
                -250876251.57881597,
                -250839576.38645324,
                -250802917.19036397,
                -250766273.98124704,
                -250729646.7498079,
                -250693035.48675886,
                -250656440.18281895,
                -250619860.82871377,
                -250583297.4151761,
                -250546749.9329449,
                -250510218.37276635,
                -250473702.72539294,
                -250437202.9815841,
                -250400719.13210604,
                -250364251.1677314,
                -250327799.0792398,
                -250291362.8574172,
                -250254942.49305668,
                -250218537.97695765,
                -250182149.29992643,
                -250145776.45277578,
                -250109419.42632523,
                -250073078.21140105,
                -250036752.79883608,
                -250000443.17946965,
                -249964149.34414807,
                -249927871.28372398,
                -249891608.98905677,
                -249855362.45101237,
                -249819131.66046348,
                -249782916.60828927,
                -249746717.28537557,
                -249710533.6826148,
                -249674365.79090595,
                -249638213.60115463,
                -249602077.10427296,
                -249565956.29117972,
            ],
            expected_r_deflect=[
                -0.005110772649589637,
                -0.005107979732115243,
                -0.005105208914710189,
                -0.005102460076784812,
                -0.00509973309862519,
                -0.0050970278613852085,
                -0.005094344247078711,
                -0.0050916821385717315,
                -0.005089041419574819,
                -0.005086421974635437,
                -0.005083823689130451,
                -0.00508124644925869,
                -0.005078690142033596,
                -0.005076154655275944,
                -0.005073639877606648,
                -0.005071145698439637,
                -0.005068672007974808,
                -0.0050662186971910635,
                -0.005063785657839406,
                -0.005061372782436119,
                -0.005058979964256019,
                -0.005056607097325768,
                -0.005054254076417275,
                -0.005051920797041146,
                -0.005049607155440222,
                -0.005047313048583168,
                -0.005045038374158145,
                -0.0050427830305665375,
                -0.005040546916916744,
                -0.005038329933018046,
                -0.005036131979374529,
                -0.005033952957179067,
                -0.00503179276830738,
                -0.005029651315312134,
                -0.005027528501417125,
                -0.0050254242305115045,
                -0.005023338407144072,
                -0.005021270936517628,
                -0.005019221724483378,
                -0.005017190677535402,
                -0.005015177702805176,
                -0.00501318270805615,
                -0.0050112056016783785,
                -0.005009246292683213,
                -0.00500730469069804,
                -0.0050053807059610815,
                -0.005003474249316234,
                -0.00500158523220798,
                -0.004999713566676328,
                -0.004997859165351823,
                -0.004996021941450596,
                -0.004994201808769465,
                -0.004992398681681088,
                -0.004990612475129163,
                -0.004988843104623673,
                -0.004987090486236186,
                -0.004985354536595191,
                -0.004983635172881489,
                -0.0049819323128236295,
                -0.004980245874693384,
                -0.004978575777301276,
                -0.004976921939992143,
                -0.004975284282640753,
                -0.004973662725647458,
                -0.00497205718993389,
                -0.004970467596938702,
                -0.00496889386861335,
                -0.004967335927417915,
                -0.004965793696316967,
                -0.0049642670987754675,
                -0.004962756058754715,
                -0.004961260500708327,
                -0.004959780349578261,
                -0.004958315530790879,
                -0.004956865970253039,
                -0.004955431594348239,
                -0.004954012329932786,
                -0.004952608104332005,
                -0.00495121884533649,
                -0.004949844481198384,
                -0.004948484940627697,
                -0.00494714015278866,
                -0.004945810047296113,
                -0.0049444945542119305,
                -0.00494319360404147,
                -0.00494190712773007,
                -0.004940635056659575,
                -0.004939377322644884,
                -0.004938133857930543,
                -0.0049369045951873705,
                -0.004935689467509104,
                -0.004934488408409089,
                -0.004933301351816993,
                -0.0049321282320755515,
                -0.004930968983937346,
                -0.004929823542561608,
                -0.0049286918435110655,
                -0.004927573822748794,
                -0.004926469416635126,
                -0.004925378561924566,
                -0.0049243011957627175,
                -0.00492767542842093,
                -0.00493105473941545,
                -0.0049344376993040034,
                -0.004937822883104445,
                -0.004941208870263711,
                -0.0049445942446273705,
                -0.004947977594409306,
                -0.004951357512161553,
                -0.004954732594744646,
                -0.00495810144329803,
                -0.004961462663210653,
                -0.004964814864092157,
                -0.004968156659744025,
                -0.004971486668130923,
                -0.004974803511352557,
                -0.004978105815615613,
                -0.00498139221120579,
                -0.004984661332460477,
                -0.004987911817741164,
                -0.004991142309406588,
                -0.004994351453785614,
                -0.004997537901150734,
                -0.005000700305691663,
                -0.0050038373254889384,
                -0.005006947622488253,
                -0.005010029862474369,
                -0.005013082715045802,
                -0.005016104853589359,
                -0.0050190949552551545,
                -0.005022051700931451,
                -0.005024973775220265,
                -0.005027859866412632,
                -0.0050307086664643785,
                -0.005033518870972115,
                -0.005036289179149139,
                -0.0050390182938020445,
                -0.005041704921306817,
                -0.005044347771585994,
                -0.0050469455580851585,
                -0.005049496997750119,
                -0.005052000811004348,
                -0.0050544557217261366,
                -0.005056860457226531,
                -0.005059213748226904,
                -0.005061514328837002,
                -0.00506376093653324,
                -0.00506595231213694,
                -0.005068087199792903,
                -0.0050701643469480695,
                -0.005072182504330391,
                -0.005074140425927964,
                -0.005076036868968126,
                -0.0050778705938968605,
                -0.005079640364358401,
                -0.0050813449471749095,
                -0.0050829831123263,
                -0.005084553632930311,
                -0.005086055285222824,
                -0.005087486848537798,
                -0.005088847105288369,
                -0.005090134840946975,
                -0.005091348844026539,
                -0.005092487906061177,
                -0.00509355082158755,
                -0.005094536388125848,
                -0.005095443406161526,
                -0.005096270679126519,
                -0.0050970170133812776,
                -0.005097681218196398,
                -0.005098262105734608,
                -0.005098758491033112,
                -0.005099169191985525,
                -0.005099493029324692,
                -0.005099728826604727,
                -0.00509987541018414,
                -0.005099931609208519,
                -0.005099896255593428,
                -0.005099768184007369,
                -0.005099546231855212,
                -0.005099229239261538,
                -0.005098816049053878,
                -0.0050983055067464755,
                -0.005097696460524048,
                -0.005096987761225497,
                -0.005096178262328027,
                -0.005095266819931027,
                -0.005094252292740603,
                -0.005093133542053485,
                -0.0050919094317419256,
                -0.005090578828238018,
                -0.005089140600518294,
                -0.005087593620088843,
                -0.005085936760970022,
                -0.005084168899681357,
                -0.0050822889152271344,
                -0.005080295689081138,
                -0.005078188105172354,
                -0.005075965049870346,
                -0.005073625411970845,
                -0.0050711680826814156,
                -0.0050709350165634055,
                -0.005070702093865363,
                -0.005070469314524161,
                -0.005070236678476709,
                -0.005070004185659953,
                -0.005069771836010877,
                -0.005069539629466499,
                -0.0050693075659638785,
                -0.005069075645440109,
                -0.00506884386783232,
                -0.0050686122330776805,
                -0.005068380741113396,
                -0.005068149391876707,
                -0.005067918185304892,
                -0.005067687121335266,
                -0.005067456199905181,
                -0.005067225420952024,
                -0.005066994784413221,
                -0.005066764290226234,
                -0.0050665339383285605,
                -0.005066303728657736,
                -0.005066073661151332,
                -0.005065843735746956,
                -0.005065613952382253,
                -0.0050653843109949035,
                -0.005065154811522624,
                -0.00506492545390317,
                -0.0050646962380743316,
                -0.005064467163973935,
                -0.0050642382315398415,
                -0.005064009440709953,
                -0.005063780791422204,
                -0.005063552283614566,
                -0.0050633239172250466,
                -0.005063095692191691,
                -0.0050628676084525795,
                -0.005062639665945829,
                -0.0050624118646095916,
                -0.0050621842043820555,
                -0.005061956685201447,
                -0.0050617293070060266,
                -0.005061502069734092,
                -0.005061274973323974,
                -0.005061048017714043,
                -0.0050608212028427045,
                -0.005060594528648399,
                -0.005060367995069603,
                -0.0050601416020448296,
                -0.005059915349512625,
                -0.005059689237411577,
                -0.005059463265680303,
                -0.005059237434257461,
                -0.005059011743081741,
                -0.005058786192091869,
                -0.005058560781226611,
                -0.005058335510424764,
                -0.005058110379625163,
                -0.005057885388766677,
                -0.005057660537788212,
                -0.005057435826628708,
                -0.005057211255227144,
                -0.005056986823522529,
                -0.005056762531453913,
                -0.005056538378960378,
                -0.0050563143659810425,
                -0.005056090492455059,
                -0.00505586675832162,
                -0.005055643163519947,
                -0.005055419707989301,
                -0.005055196391668977,
                -0.005054973214498306,
                -0.005054750176416653,
                -0.005054527277363419,
                -0.005054304517278041,
                -0.0050540818960999895,
                -0.005053859413768771,
                -0.005053637070223927,
                -0.005053414865405034,
                -0.0050531927992517045,
                -0.005052970871703585,
                -0.005052749082700356,
                -0.005052527432181736,
                -0.005052305920087477,
                -0.005052084546357364,
                -0.005051863310931219,
                -0.005051642213748901,
                -0.005051421254750298,
                -0.005051200433875337,
                -0.005050979751063981,
                -0.005050759206256223,
                -0.005050538799392096,
                -0.005050318530411664,
                -0.005050098399255028,
                -0.005049878405862322,
                -0.005049658550173715,
                -0.005049438832129412,
                -0.005049219251669652,
                -0.005048999808734708,
                -0.0050487805032648865,
            ],
            expected_rradius=[
                2.9939411851091102,
                2.999416570397585,
                3.0048919556860594,
                3.010367340974534,
                3.0158427262630085,
                3.021318111551483,
                3.0267934968399577,
                3.0322688821284323,
                3.037744267416907,
                3.043219652705382,
                3.0486950379938564,
                3.054170423282331,
                3.0596458085708056,
                3.06512119385928,
                3.0705965791477547,
                3.0760719644362293,
                3.081547349724704,
                3.0870227350131785,
                3.092498120301653,
                3.0979735055901276,
                3.103448890878602,
                3.108924276167077,
                3.1143996614555514,
                3.119875046744026,
                3.1253504320325005,
                3.130825817320975,
                3.13630120260945,
                3.1417765878979247,
                3.1472519731863993,
                3.152727358474874,
                3.1582027437633484,
                3.163678129051823,
                3.1691535143402976,
                3.174628899628772,
                3.1801042849172467,
                3.1855796702057213,
                3.191055055494196,
                3.1965304407826705,
                3.202005826071145,
                3.2074812113596196,
                3.212956596648094,
                3.218431981936569,
                3.223907367225044,
                3.2293827525135184,
                3.234858137801993,
                3.2403335230904675,
                3.245808908378942,
                3.2512842936674167,
                3.2567596789558912,
                3.262235064244366,
                3.2677104495328404,
                3.273185834821315,
                3.2786612201097896,
                3.284136605398264,
                3.2896119906867387,
                3.2950873759752133,
                3.300562761263688,
                3.3060381465521624,
                3.311513531840637,
                3.316988917129112,
                3.3224643024175866,
                3.327939687706061,
                3.3334150729945358,
                3.3388904582830103,
                3.344365843571485,
                3.3498412288599595,
                3.355316614148434,
                3.3607919994369087,
                3.3662673847253832,
                3.371742770013858,
                3.3772181553023324,
                3.382693540590807,
                3.3881689258792815,
                3.393644311167756,
                3.399119696456231,
                3.4045950817447057,
                3.4100704670331803,
                3.415545852321655,
                3.4210212376101294,
                3.426496622898604,
                3.4319720081870786,
                3.437447393475553,
                3.4429227787640277,
                3.4483981640525023,
                3.453873549340977,
                3.4593489346294515,
                3.464824319917926,
                3.4702997052064006,
                3.475775090494875,
                3.48125047578335,
                3.4867258610718244,
                3.492201246360299,
                3.4976766316487735,
                3.5031520169372485,
                3.508627402225723,
                3.5141027875141977,
                3.5195781728026723,
                3.525053558091147,
                3.5305289433796214,
                3.536004328668096,
                3.5414797139565706,
                3.5469411197215766,
                3.5524025254865825,
                3.557863931251588,
                3.563325337016594,
                3.5687867427816,
                3.574248148546606,
                3.579709554311612,
                3.5851709600766175,
                3.5906323658416235,
                3.5960937716066295,
                3.6015551773716354,
                3.6070165831366414,
                3.6124779889016474,
                3.617939394666653,
                3.623400800431659,
                3.628862206196665,
                3.634323611961671,
                3.639785017726677,
                3.6452464234916824,
                3.6507078292566884,
                3.6561692350216943,
                3.6616306407867003,
                3.6670920465517063,
                3.672553452316712,
                3.678014858081718,
                3.6834762638467238,
                3.6889376696117298,
                3.6943990753767357,
                3.6998604811417417,
                3.7053218869067472,
                3.710783292671753,
                3.716244698436759,
                3.721706104201765,
                3.727167509966771,
                3.7326289157317767,
                3.7380903214967827,
                3.7435517272617886,
                3.7490131330267946,
                3.7544745387918006,
                3.759935944556806,
                3.765397350321812,
                3.770858756086818,
                3.776320161851824,
                3.78178156761683,
                3.7872429733818356,
                3.7927043791468416,
                3.7981657849118475,
                3.8036271906768535,
                3.8090885964418595,
                3.8145500022068655,
                3.820011407971871,
                3.825472813736877,
                3.830934219501883,
                3.836395625266889,
                3.841857031031895,
                3.8473184367969004,
                3.8527798425619064,
                3.8582412483269124,
                3.8637026540919184,
                3.8691640598569244,
                3.8746254656219303,
                3.880086871386936,
                3.885548277151942,
                3.891009682916948,
                3.896471088681954,
                3.9019324944469593,
                3.9073939002119653,
                3.9128553059769713,
                3.9183167117419773,
                3.9237781175069832,
                3.929239523271989,
                3.9347009290369948,
                3.9401623348020007,
                3.9456237405670067,
                3.9510851463320127,
                3.9565465520970187,
                3.962007957862024,
                3.96746936362703,
                3.972930769392036,
                3.978392175157042,
                3.983853580922048,
                3.989314986687054,
                3.9947763924520596,
                4.000237798217066,
                4.005699203982071,
                4.011160609747077,
                4.016622015512083,
                4.022083421277089,
                4.027544827042095,
                4.033006232807101,
                4.038467638572107,
                4.043929044337113,
                4.049390450102119,
                4.054851855867125,
                4.06031326163213,
                4.065774667397136,
                4.071236073162142,
                4.076697478927148,
                4.082158884692154,
                4.08762029045716,
                4.08822029045716,
                4.08882029045716,
                4.08942029045716,
                4.09002029045716,
                4.09062029045716,
                4.0912202904571595,
                4.09182029045716,
                4.09242029045716,
                4.09302029045716,
                4.09362029045716,
                4.09422029045716,
                4.09482029045716,
                4.0954202904571595,
                4.09602029045716,
                4.09662029045716,
                4.09722029045716,
                4.09782029045716,
                4.09842029045716,
                4.09902029045716,
                4.0996202904571595,
                4.10022029045716,
                4.10082029045716,
                4.10142029045716,
                4.10202029045716,
                4.10262029045716,
                4.10322029045716,
                4.1038202904571595,
                4.10442029045716,
                4.10502029045716,
                4.10562029045716,
                4.10622029045716,
                4.10682029045716,
                4.10742029045716,
                4.108020290457159,
                4.10862029045716,
                4.10922029045716,
                4.10982029045716,
                4.11042029045716,
                4.1110202904571596,
                4.11162029045716,
                4.112220290457159,
                4.11282029045716,
                4.11342029045716,
                4.11402029045716,
                4.11462029045716,
                4.1152202904571595,
                4.11582029045716,
                4.116420290457159,
                4.11702029045716,
                4.117620290457159,
                4.11822029045716,
                4.11882029045716,
                4.1194202904571595,
                4.12002029045716,
                4.120620290457159,
                4.12122029045716,
                4.121820290457159,
                4.12242029045716,
                4.12302029045716,
                4.1236202904571595,
                4.12422029045716,
                4.124820290457159,
                4.12542029045716,
                4.126020290457159,
                4.12662029045716,
                4.12722029045716,
                4.1278202904571595,
                4.12842029045716,
                4.129020290457159,
                4.12962029045716,
                4.130220290457159,
                4.13082029045716,
                4.13142029045716,
                4.1320202904571595,
                4.13262029045716,
                4.133220290457159,
                4.13382029045716,
                4.134420290457159,
                4.13502029045716,
                4.13562029045716,
                4.136220290457159,
                4.13682029045716,
                4.137420290457159,
                4.13802029045716,
                4.138620290457159,
                4.13922029045716,
                4.13982029045716,
                4.140420290457159,
                4.14102029045716,
                4.141620290457159,
                4.14222029045716,
                4.142820290457159,
                4.1434202904571595,
                4.14402029045716,
                4.144620290457159,
                4.14522029045716,
                4.145820290457159,
                4.14642029045716,
                4.147020290457159,
            ],
        ),
        PlaneStressParam(
            n_radial_array=100,
            nlayers=3,
            nu=numpy.array(
                numpy.array(
                    (0.29999999999999999, 0.30904421667064924, 0.29999999999999999),
                    order="F",
                ),
                order="F",
            ).transpose(),
            rad=numpy.array(
                numpy.array(
                    (
                        2.9939411851091102,
                        3.5414797139565706,
                        4.0876202904571599,
                        4.1476202904571595,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            ey=numpy.array(
                numpy.array(
                    (205000000000, 43126670035.025253, 205000000000), order="F"
                ),
                order="F",
            ).transpose(),
            j=numpy.array(
                numpy.array((0, 18097185.781970859, 0), order="F"), order="F"
            ).transpose(),
            expected_sigr=[
                0.0,
                -638231.855766482,
                -1272978.0149700674,
                -1904263.8142719064,
                -2532114.3605434406,
                -3156554.533362813,
                -3777608.9874793817,
                -4395302.155247777,
                -5009658.2490304215,
                -5620701.263570679,
                -6228454.97833536,
                -6832942.959828614,
                -7434188.563876412,
                -8032214.937882571,
                -8627045.023056716,
                -9218701.556614282,
                -9807207.07394974,
                -10392583.910782054,
                -10974854.205274347,
                -11554039.900126804,
                -12130162.744643984,
                -12703244.296776803,
                -13273305.925139055,
                -13840368.810999524,
                -14404453.950249525,
                -14965582.155346265,
                -15523774.05723289,
                -16079050.107234456,
                -16631430.578931632,
                -17180935.57001094,
                -17727585.00409306,
                -18271398.632538624,
                -18812396.036232498,
                -19350596.627346307,
                -19886019.651079975,
                -20418684.18738201,
                -20948609.1526495,
                -21475813.301407494,
                -22000315.22796835,
                -22522133.368071266,
                -23041286.00050247,
                -23557791.24869583,
                -24071667.082314458,
                -24582931.31881398,
                -25091601.624986775,
                -25597695.518488254,
                -26101230.36934514,
                -26602223.401445802,
                -27100691.6940132,
                -27596652.183060553,
                -28090121.6628299,
                -28581116.78721381,
                -29069654.071160544,
                -29555749.892062735,
                -30039420.491130088,
                -30520681.974745892,
                -30999550.31580798,
                -31476041.355053987,
                -31950170.802371547,
                -32421954.23809305,
                -32891407.114275593,
                -33358544.755966377,
                -33823382.3624533,
                -34285935.00850141,
                -34746217.64557501,
                -35204245.10304597,
                -35660032.08938813,
                -36113593.19335821,
                -36564942.88516308,
                -37014095.51761406,
                -37461065.327267796,
                -37905866.4355546,
                -38348512.84989362,
                -38789018.46479585,
                -39227397.06295431,
                -39663662.31632225,
                -40097827.78717919,
                -40529906.92918485,
                -40959913.08842154,
                -41387859.50442459,
                -41813759.31120163,
                -42237625.5382402,
                -42659471.11150431,
                -43079308.85441998,
                -43497151.48884965,
                -43913011.6360561,
                -44326901.817655385,
                -44738834.45655951,
                -45148821.87790862,
                -45556876.30999292,
                -45963009.88516445,
                -46367234.64073896,
                -46769562.51988788,
                -47170005.372520424,
                -47568574.956156105,
                -47965282.936788,
                -48360140.88973615,
                -48753160.30049217,
                -49144352.56555428,
                -49533728.99325363,
                -49921300.804571584,
                -49954359.75079646,
                -49974994.421951994,
                -49983237.62599215,
                -49979121.95724278,
                -49962679.79819553,
                -49933943.321276575,
                -49892944.49059213,
                -49839715.063674085,
                -49774286.59318251,
                -49696690.42861684,
                -49606957.71798797,
                -49505119.40949372,
                -49391206.253157824,
                -49265248.802470304,
                -49127277.41599873,
                -48977322.25898817,
                -48815413.304947995,
                -48641580.33721941,
                -48455852.950529814,
                -48258260.55252331,
                -48048832.365294494,
                -47827597.42688793,
                -47594584.59279222,
                -47349822.53742223,
                -47093339.75557484,
                -46825164.56388664,
                -46545325.10226525,
                -46253849.335310586,
                -45950765.05372232,
                -45636099.875691466,
                -45309881.24828925,
                -44972136.44881964,
                -44622892.586189084,
                -44262176.60223319,
                -43890015.27305536,
                -43506435.21033433,
                -43111462.86262834,
                -42705124.516666815,
                -42287446.29862875,
                -41858454.17540391,
                -41418173.95584675,
                -40966631.29202375,
                -40503851.68043686,
                -40029860.46324105,
                -39544682.82945457,
                -39048343.816151306,
                -38540868.309644595,
                -38022281.04666105,
                -37492606.61550055,
                -36951869.45718825,
                -36400093.8666107,
                -35837303.99364805,
                -35263523.8442963,
                -34678777.28176432,
                -34083088.027581,
                -33476479.662681654,
                -32858975.628480688,
                -32230599.22794289,
                -31591373.626639504,
                -30941321.853795685,
                -30280466.803329386,
                -29608831.23487964,
                -28926437.774824906,
                -28233308.91729218,
                -27529467.025158796,
                -26814934.33104299,
                -26089732.938284557,
                -25353884.821918193,
                -24607411.829640254,
                -23850335.682763584,
                -23082677.977158453,
                -22304460.184201445,
                -21515703.651693877,
                -20716429.60478688,
                -19906659.14690248,
                -19086413.26062489,
                -18255712.808604427,
                -17414578.534444164,
                -16563031.063583972,
                -15701090.904165361,
                -14828778.44790294,
                -13946113.970938675,
                -13053117.634692835,
                -12149809.48670104,
                -11236209.461454624,
                -10312337.381225117,
                -9378212.956880113,
                -8433855.788701806,
                -7479285.367186312,
                -6514521.073844995,
                -5539582.181995865,
                -4554487.85754239,
                -3559257.1597569585,
                -2553909.0420464603,
                -1538462.3527192033,
                -512935.83573593997,
                522651.86853248376,
                1568282.2225666705,
                2623936.7909615375,
                3689597.239691208,
                3651895.231238659,
                3614209.818926577,
                3576540.9930157103,
                3538888.743774036,
                3501253.061476492,
                3463633.9364053686,
                3426031.358849868,
                3388445.3191064713,
                3350875.807478765,
                3313322.8142772503,
                3275786.3298197784,
                3238266.3444311135,
                3200762.8484432492,
                3163275.832195042,
                3125805.2860326516,
                3088351.200309296,
                3050913.5653850324,
                3013492.3716272456,
                2976087.6094101095,
                2938699.2691151258,
                2901327.3411305123,
                2863971.8158518155,
                2826632.683681444,
                2789309.93502884,
                2752003.560310531,
                2714713.5499499296,
                2677439.8943775576,
                2640182.5840308485,
                2602941.609354294,
                2565716.9607993956,
                2528508.6288244454,
                2491316.6038949895,
                2454140.8764832406,
                2416981.4370685695,
                2379838.2761372332,
                2342711.3841824015,
                2305600.7517043273,
                2268506.36920998,
                2231428.2272134367,
                2194366.316235539,
                2157320.6268042377,
                2120291.149454175,
                2083277.8747270512,
                2046280.7931714796,
                2009299.8953427651,
                1972335.1718033939,
                1935386.6131224218,
                1898454.2098760372,
                1861537.9526470467,
                1824637.8320254136,
                1787753.838607672,
                1750885.962997364,
                1714034.1958049217,
                1677198.5276473942,
                1640378.949148915,
                1603575.450940211,
                1566788.0236590446,
                1530016.6579497962,
                1493261.3444638075,
                1456522.0738591857,
                1419798.8368006812,
                1383091.6239601525,
                1346400.4260158325,
                1309725.23365311,
                1273066.0375638467,
                1236422.8284469144,
                1199795.5970078032,
                1163184.33395872,
                1126589.0300187848,
                1090009.6759136862,
                1053446.262375977,
                1016898.7801448042,
                980367.2199661784,
                943851.5725928267,
                907351.8287839973,
                870867.979305924,
                834400.0149312643,
                797947.9264396608,
                761511.704617082,
                725091.3402565308,
                688686.8241575557,
                652298.1471262753,
                615925.2999756474,
                579568.2735251015,
                543227.0586009554,
                506901.646035926,
                470592.02666951995,
                434298.1913479606,
                398020.1309238218,
                361757.8362566382,
                325511.29821222165,
                289280.5076633691,
                253065.4554891793,
                216866.13257544299,
                180682.52981469216,
                144514.63810580914,
                108362.44835451504,
                72225.95147280802,
                36105.13837954948,
            ],
            expected_sigt=[
                -349942877.4275314,
                -349304645.57176495,
                -348669899.41256136,
                -348038613.61325955,
                -347410763.066988,
                -346786322.8941686,
                -346165268.44005203,
                -345547575.2722837,
                -344933219.178501,
                -344322176.1639607,
                -343714422.44919604,
                -343109934.4677028,
                -342508688.863655,
                -341910662.4896488,
                -341315832.4044747,
                -340724175.87091714,
                -340135670.35358167,
                -339550293.5167493,
                -338968023.222257,
                -338388837.5274046,
                -337812714.68288743,
                -337239633.1307546,
                -336669571.5023923,
                -336102508.61653185,
                -335538423.4772818,
                -334977295.2721851,
                -334419103.37029845,
                -333863827.3202969,
                -333311446.84859973,
                -332761941.85752034,
                -332215292.42343825,
                -331671478.7949927,
                -331130481.39129883,
                -330592280.800185,
                -330056857.7764513,
                -329524193.2401493,
                -328994268.2748818,
                -328467064.1261238,
                -327942562.1995629,
                -327420744.05946004,
                -326901591.42702883,
                -326385086.17883545,
                -325871210.3452168,
                -325359946.10871726,
                -324851275.8025445,
                -324345181.90904295,
                -323841647.0581861,
                -323340654.02608544,
                -322842185.73351806,
                -322346225.24447066,
                -321852755.7647013,
                -321361760.64031744,
                -320873223.3563707,
                -320387127.5354685,
                -319903456.9364011,
                -319422195.4527854,
                -318943327.11172324,
                -318466836.0724772,
                -317992706.6251597,
                -317520923.18943816,
                -317051470.3132556,
                -316584332.6715648,
                -316119495.0650779,
                -315656942.4190298,
                -315196659.7819562,
                -314738632.32448524,
                -314282845.33814305,
                -313829284.23417294,
                -313377934.54236805,
                -312928781.9099172,
                -312481812.10026336,
                -312037010.99197656,
                -311594364.57763755,
                -311153858.9627353,
                -310715480.3645769,
                -310279215.1112089,
                -309845049.640352,
                -309412970.4983463,
                -308982964.3391096,
                -308555017.92310655,
                -308129118.11632955,
                -307705251.8892909,
                -307283406.3160268,
                -306863568.5731112,
                -306445725.9386815,
                -306029865.791475,
                -305615975.60987574,
                -305204042.97097164,
                -304794055.5496225,
                -304386001.1175382,
                -303979867.5423667,
                -303575642.78679216,
                -303173314.9076432,
                -302772872.0550107,
                -302374302.471375,
                -301977594.4907431,
                -301582736.53779495,
                -301189717.12703896,
                -300798524.8619768,
                -300409148.4342775,
                -75393985.39087032,
                -75352895.95330013,
                -75308186.27146304,
                -75259848.14709243,
                -75207873.48230305,
                -75152254.27849749,
                -75092982.635284,
                -75030050.74940914,
                -74963450.91370073,
                -74893175.51602733,
                -74819217.03826267,
                -74741568.05526796,
                -74660221.23388477,
                -74575169.33193818,
                -74486405.1972532,
                -74393921.76667804,
                -74297712.06512664,
                -74197769.20462456,
                -74094086.38337062,
                -73986656.88480622,
                -73875474.07669795,
                -73760531.41022801,
                -73641822.41909933,
                -73519340.71864453,
                -73393080.00495018,
                -73263034.05399033,
                -73129196.72076792,
                -72991561.93846779,
                -72850123.71761769,
                -72704876.14526016,
                -72555813.3841336,
                -72402929.67186344,
                -72246219.32015613,
                -72085676.71401668,
                -71921296.31095727,
                -71753072.64022708,
                -71581000.30204895,
                -71405073.9668568,
                -71225288.37455449,
                -71041638.33377115,
                -70854118.72113262,
                -70662724.48053451,
                -70467450.62243325,
                -70268292.22313243,
                -70065244.4240877,
                -69858302.43121327,
                -69647461.51419875,
                -69432717.00583516,
                -69214064.30134167,
                -68991498.85770856,
                -68765016.19304176,
                -68534611.8859144,
                -68300281.57473324,
                -68062020.95709684,
                -67819825.78918111,
                -67573691.88511074,
                -67323615.11635438,
                -67069591.41111759,
                -66811616.75374359,
                -66549687.18412141,
                -66283798.797103256,
                -66013947.74192255,
                -65740130.22162387,
                -65462342.49249586,
                -65180580.86350979,
                -64894841.695768714,
                -64605121.401957214,
                -64311416.44580017,
                -64013723.341527745,
                -63712038.653342426,
                -63406358.99489821,
                -63096681.02878188,
                -62783001.46599489,
                -62465317.06545305,
                -62143624.63347747,
                -61817921.02330379,
                -61488203.13458894,
                -61154467.91292381,
                -60816712.34935414,
                -60474933.47990544,
                -60129128.38510988,
                -59779294.189546056,
                -59425428.061372176,
                -59067527.21187808,
                -58705588.895027705,
                -58339610.40701998,
                -57969589.08584027,
                -57595522.31083187,
                -57217407.502260946,
                -56835242.12088977,
                -56449023.667553514,
                -56058749.68274543,
                -55664417.74619955,
                -55266025.4764855,
                -54863570.530602135,
                -54457050.60357594,
                -54046463.4280663,
                -53631806.773974754,
                -53213078.44805462,
                -52790276.293530844,
                -253219448.39249134,
                -253181746.38403878,
                -253144060.97172672,
                -253106392.14581582,
                -253068739.89657417,
                -253031104.21427664,
                -252993485.0892055,
                -252955882.51165,
                -252918296.4719066,
                -252880726.9602789,
                -252843173.96707734,
                -252805637.48261994,
                -252768117.49723122,
                -252730614.00124338,
                -252693126.98499516,
                -252655656.43883282,
                -252618202.35310942,
                -252580764.71818522,
                -252543343.52442738,
                -252505938.76221025,
                -252468550.42191526,
                -252431178.4939306,
                -252393822.96865195,
                -252356483.8364816,
                -252319161.08782896,
                -252281854.7131107,
                -252244564.70275006,
                -252207291.04717773,
                -252170033.736831,
                -252132792.76215443,
                -252095568.1135995,
                -252058359.7816246,
                -252021167.75669512,
                -251983992.0292834,
                -251946832.5898687,
                -251909689.42893735,
                -251872562.53698257,
                -251835451.90450448,
                -251798357.5220101,
                -251761279.38001359,
                -251724217.46903566,
                -251687171.77960438,
                -251650142.3022543,
                -251613129.02752718,
                -251576131.94597158,
                -251539151.04814288,
                -251502186.32460356,
                -251465237.76592252,
                -251428305.36267614,
                -251391389.1054472,
                -251354488.98482555,
                -251317604.99140784,
                -251280737.1157975,
                -251243885.34860504,
                -251207049.68044755,
                -251170230.10194904,
                -251133426.60374033,
                -251096639.1764592,
                -251059867.81074992,
                -251023112.4972639,
                -250986373.22665933,
                -250949649.98960084,
                -250912942.77676025,
                -250876251.57881597,
                -250839576.38645324,
                -250802917.19036397,
                -250766273.98124704,
                -250729646.7498079,
                -250693035.48675886,
                -250656440.18281895,
                -250619860.82871377,
                -250583297.4151761,
                -250546749.9329449,
                -250510218.37276635,
                -250473702.72539294,
                -250437202.9815841,
                -250400719.13210604,
                -250364251.1677314,
                -250327799.0792398,
                -250291362.8574172,
                -250254942.49305668,
                -250218537.97695765,
                -250182149.29992643,
                -250145776.45277578,
                -250109419.42632523,
                -250073078.21140105,
                -250036752.79883608,
                -250000443.17946965,
                -249964149.34414807,
                -249927871.28372398,
                -249891608.98905677,
                -249855362.45101237,
                -249819131.66046348,
                -249782916.60828927,
                -249746717.28537557,
                -249710533.6826148,
                -249674365.79090595,
                -249638213.60115463,
                -249602077.10427296,
                -249565956.29117972,
            ],
            expected_r_deflect=[
                -0.005110772649589637,
                -0.005107979732115243,
                -0.005105208914710189,
                -0.005102460076784812,
                -0.00509973309862519,
                -0.0050970278613852085,
                -0.005094344247078711,
                -0.0050916821385717315,
                -0.005089041419574819,
                -0.005086421974635437,
                -0.005083823689130451,
                -0.00508124644925869,
                -0.005078690142033596,
                -0.005076154655275944,
                -0.005073639877606648,
                -0.005071145698439637,
                -0.005068672007974808,
                -0.0050662186971910635,
                -0.005063785657839406,
                -0.005061372782436119,
                -0.005058979964256019,
                -0.005056607097325768,
                -0.005054254076417275,
                -0.005051920797041146,
                -0.005049607155440222,
                -0.005047313048583168,
                -0.005045038374158145,
                -0.0050427830305665375,
                -0.005040546916916744,
                -0.005038329933018046,
                -0.005036131979374529,
                -0.005033952957179067,
                -0.00503179276830738,
                -0.005029651315312134,
                -0.005027528501417125,
                -0.0050254242305115045,
                -0.005023338407144072,
                -0.005021270936517628,
                -0.005019221724483378,
                -0.005017190677535402,
                -0.005015177702805176,
                -0.00501318270805615,
                -0.0050112056016783785,
                -0.005009246292683213,
                -0.00500730469069804,
                -0.0050053807059610815,
                -0.005003474249316234,
                -0.00500158523220798,
                -0.004999713566676328,
                -0.004997859165351823,
                -0.004996021941450596,
                -0.004994201808769465,
                -0.004992398681681088,
                -0.004990612475129163,
                -0.004988843104623673,
                -0.004987090486236186,
                -0.004985354536595191,
                -0.004983635172881489,
                -0.0049819323128236295,
                -0.004980245874693384,
                -0.004978575777301276,
                -0.004976921939992143,
                -0.004975284282640753,
                -0.004973662725647458,
                -0.00497205718993389,
                -0.004970467596938702,
                -0.00496889386861335,
                -0.004967335927417915,
                -0.004965793696316967,
                -0.0049642670987754675,
                -0.004962756058754715,
                -0.004961260500708327,
                -0.004959780349578261,
                -0.004958315530790879,
                -0.004956865970253039,
                -0.004955431594348239,
                -0.004954012329932786,
                -0.004952608104332005,
                -0.00495121884533649,
                -0.004949844481198384,
                -0.004948484940627697,
                -0.00494714015278866,
                -0.004945810047296113,
                -0.0049444945542119305,
                -0.00494319360404147,
                -0.00494190712773007,
                -0.004940635056659575,
                -0.004939377322644884,
                -0.004938133857930543,
                -0.0049369045951873705,
                -0.004935689467509104,
                -0.004934488408409089,
                -0.004933301351816993,
                -0.0049321282320755515,
                -0.004930968983937346,
                -0.004929823542561608,
                -0.0049286918435110655,
                -0.004927573822748794,
                -0.004926469416635126,
                -0.004925378561924566,
                -0.0049243011957627175,
                -0.00492767542842093,
                -0.00493105473941545,
                -0.0049344376993040034,
                -0.004937822883104445,
                -0.004941208870263711,
                -0.0049445942446273705,
                -0.004947977594409306,
                -0.004951357512161553,
                -0.004954732594744646,
                -0.00495810144329803,
                -0.004961462663210653,
                -0.004964814864092157,
                -0.004968156659744025,
                -0.004971486668130923,
                -0.004974803511352557,
                -0.004978105815615613,
                -0.00498139221120579,
                -0.004984661332460477,
                -0.004987911817741164,
                -0.004991142309406588,
                -0.004994351453785614,
                -0.004997537901150734,
                -0.005000700305691663,
                -0.0050038373254889384,
                -0.005006947622488253,
                -0.005010029862474369,
                -0.005013082715045802,
                -0.005016104853589359,
                -0.0050190949552551545,
                -0.005022051700931451,
                -0.005024973775220265,
                -0.005027859866412632,
                -0.0050307086664643785,
                -0.005033518870972115,
                -0.005036289179149139,
                -0.0050390182938020445,
                -0.005041704921306817,
                -0.005044347771585994,
                -0.0050469455580851585,
                -0.005049496997750119,
                -0.005052000811004348,
                -0.0050544557217261366,
                -0.005056860457226531,
                -0.005059213748226904,
                -0.005061514328837002,
                -0.00506376093653324,
                -0.00506595231213694,
                -0.005068087199792903,
                -0.0050701643469480695,
                -0.005072182504330391,
                -0.005074140425927964,
                -0.005076036868968126,
                -0.0050778705938968605,
                -0.005079640364358401,
                -0.0050813449471749095,
                -0.0050829831123263,
                -0.005084553632930311,
                -0.005086055285222824,
                -0.005087486848537798,
                -0.005088847105288369,
                -0.005090134840946975,
                -0.005091348844026539,
                -0.005092487906061177,
                -0.00509355082158755,
                -0.005094536388125848,
                -0.005095443406161526,
                -0.005096270679126519,
                -0.0050970170133812776,
                -0.005097681218196398,
                -0.005098262105734608,
                -0.005098758491033112,
                -0.005099169191985525,
                -0.005099493029324692,
                -0.005099728826604727,
                -0.00509987541018414,
                -0.005099931609208519,
                -0.005099896255593428,
                -0.005099768184007369,
                -0.005099546231855212,
                -0.005099229239261538,
                -0.005098816049053878,
                -0.0050983055067464755,
                -0.005097696460524048,
                -0.005096987761225497,
                -0.005096178262328027,
                -0.005095266819931027,
                -0.005094252292740603,
                -0.005093133542053485,
                -0.0050919094317419256,
                -0.005090578828238018,
                -0.005089140600518294,
                -0.005087593620088843,
                -0.005085936760970022,
                -0.005084168899681357,
                -0.0050822889152271344,
                -0.005080295689081138,
                -0.005078188105172354,
                -0.005075965049870346,
                -0.005073625411970845,
                -0.0050711680826814156,
                -0.0050709350165634055,
                -0.005070702093865363,
                -0.005070469314524161,
                -0.005070236678476709,
                -0.005070004185659953,
                -0.005069771836010877,
                -0.005069539629466499,
                -0.0050693075659638785,
                -0.005069075645440109,
                -0.00506884386783232,
                -0.0050686122330776805,
                -0.005068380741113396,
                -0.005068149391876707,
                -0.005067918185304892,
                -0.005067687121335266,
                -0.005067456199905181,
                -0.005067225420952024,
                -0.005066994784413221,
                -0.005066764290226234,
                -0.0050665339383285605,
                -0.005066303728657736,
                -0.005066073661151332,
                -0.005065843735746956,
                -0.005065613952382253,
                -0.0050653843109949035,
                -0.005065154811522624,
                -0.00506492545390317,
                -0.0050646962380743316,
                -0.005064467163973935,
                -0.0050642382315398415,
                -0.005064009440709953,
                -0.005063780791422204,
                -0.005063552283614566,
                -0.0050633239172250466,
                -0.005063095692191691,
                -0.0050628676084525795,
                -0.005062639665945829,
                -0.0050624118646095916,
                -0.0050621842043820555,
                -0.005061956685201447,
                -0.0050617293070060266,
                -0.005061502069734092,
                -0.005061274973323974,
                -0.005061048017714043,
                -0.0050608212028427045,
                -0.005060594528648399,
                -0.005060367995069603,
                -0.0050601416020448296,
                -0.005059915349512625,
                -0.005059689237411577,
                -0.005059463265680303,
                -0.005059237434257461,
                -0.005059011743081741,
                -0.005058786192091869,
                -0.005058560781226611,
                -0.005058335510424764,
                -0.005058110379625163,
                -0.005057885388766677,
                -0.005057660537788212,
                -0.005057435826628708,
                -0.005057211255227144,
                -0.005056986823522529,
                -0.005056762531453913,
                -0.005056538378960378,
                -0.0050563143659810425,
                -0.005056090492455059,
                -0.00505586675832162,
                -0.005055643163519947,
                -0.005055419707989301,
                -0.005055196391668977,
                -0.005054973214498306,
                -0.005054750176416653,
                -0.005054527277363419,
                -0.005054304517278041,
                -0.0050540818960999895,
                -0.005053859413768771,
                -0.005053637070223927,
                -0.005053414865405034,
                -0.0050531927992517045,
                -0.005052970871703585,
                -0.005052749082700356,
                -0.005052527432181736,
                -0.005052305920087477,
                -0.005052084546357364,
                -0.005051863310931219,
                -0.005051642213748901,
                -0.005051421254750298,
                -0.005051200433875337,
                -0.005050979751063981,
                -0.005050759206256223,
                -0.005050538799392096,
                -0.005050318530411664,
                -0.005050098399255028,
                -0.005049878405862322,
                -0.005049658550173715,
                -0.005049438832129412,
                -0.005049219251669652,
                -0.005048999808734708,
                -0.0050487805032648865,
            ],
            expected_rradius=[
                2.9939411851091102,
                2.999416570397585,
                3.0048919556860594,
                3.010367340974534,
                3.0158427262630085,
                3.021318111551483,
                3.0267934968399577,
                3.0322688821284323,
                3.037744267416907,
                3.043219652705382,
                3.0486950379938564,
                3.054170423282331,
                3.0596458085708056,
                3.06512119385928,
                3.0705965791477547,
                3.0760719644362293,
                3.081547349724704,
                3.0870227350131785,
                3.092498120301653,
                3.0979735055901276,
                3.103448890878602,
                3.108924276167077,
                3.1143996614555514,
                3.119875046744026,
                3.1253504320325005,
                3.130825817320975,
                3.13630120260945,
                3.1417765878979247,
                3.1472519731863993,
                3.152727358474874,
                3.1582027437633484,
                3.163678129051823,
                3.1691535143402976,
                3.174628899628772,
                3.1801042849172467,
                3.1855796702057213,
                3.191055055494196,
                3.1965304407826705,
                3.202005826071145,
                3.2074812113596196,
                3.212956596648094,
                3.218431981936569,
                3.223907367225044,
                3.2293827525135184,
                3.234858137801993,
                3.2403335230904675,
                3.245808908378942,
                3.2512842936674167,
                3.2567596789558912,
                3.262235064244366,
                3.2677104495328404,
                3.273185834821315,
                3.2786612201097896,
                3.284136605398264,
                3.2896119906867387,
                3.2950873759752133,
                3.300562761263688,
                3.3060381465521624,
                3.311513531840637,
                3.316988917129112,
                3.3224643024175866,
                3.327939687706061,
                3.3334150729945358,
                3.3388904582830103,
                3.344365843571485,
                3.3498412288599595,
                3.355316614148434,
                3.3607919994369087,
                3.3662673847253832,
                3.371742770013858,
                3.3772181553023324,
                3.382693540590807,
                3.3881689258792815,
                3.393644311167756,
                3.399119696456231,
                3.4045950817447057,
                3.4100704670331803,
                3.415545852321655,
                3.4210212376101294,
                3.426496622898604,
                3.4319720081870786,
                3.437447393475553,
                3.4429227787640277,
                3.4483981640525023,
                3.453873549340977,
                3.4593489346294515,
                3.464824319917926,
                3.4702997052064006,
                3.475775090494875,
                3.48125047578335,
                3.4867258610718244,
                3.492201246360299,
                3.4976766316487735,
                3.5031520169372485,
                3.508627402225723,
                3.5141027875141977,
                3.5195781728026723,
                3.525053558091147,
                3.5305289433796214,
                3.536004328668096,
                3.5414797139565706,
                3.5469411197215766,
                3.5524025254865825,
                3.557863931251588,
                3.563325337016594,
                3.5687867427816,
                3.574248148546606,
                3.579709554311612,
                3.5851709600766175,
                3.5906323658416235,
                3.5960937716066295,
                3.6015551773716354,
                3.6070165831366414,
                3.6124779889016474,
                3.617939394666653,
                3.623400800431659,
                3.628862206196665,
                3.634323611961671,
                3.639785017726677,
                3.6452464234916824,
                3.6507078292566884,
                3.6561692350216943,
                3.6616306407867003,
                3.6670920465517063,
                3.672553452316712,
                3.678014858081718,
                3.6834762638467238,
                3.6889376696117298,
                3.6943990753767357,
                3.6998604811417417,
                3.7053218869067472,
                3.710783292671753,
                3.716244698436759,
                3.721706104201765,
                3.727167509966771,
                3.7326289157317767,
                3.7380903214967827,
                3.7435517272617886,
                3.7490131330267946,
                3.7544745387918006,
                3.759935944556806,
                3.765397350321812,
                3.770858756086818,
                3.776320161851824,
                3.78178156761683,
                3.7872429733818356,
                3.7927043791468416,
                3.7981657849118475,
                3.8036271906768535,
                3.8090885964418595,
                3.8145500022068655,
                3.820011407971871,
                3.825472813736877,
                3.830934219501883,
                3.836395625266889,
                3.841857031031895,
                3.8473184367969004,
                3.8527798425619064,
                3.8582412483269124,
                3.8637026540919184,
                3.8691640598569244,
                3.8746254656219303,
                3.880086871386936,
                3.885548277151942,
                3.891009682916948,
                3.896471088681954,
                3.9019324944469593,
                3.9073939002119653,
                3.9128553059769713,
                3.9183167117419773,
                3.9237781175069832,
                3.929239523271989,
                3.9347009290369948,
                3.9401623348020007,
                3.9456237405670067,
                3.9510851463320127,
                3.9565465520970187,
                3.962007957862024,
                3.96746936362703,
                3.972930769392036,
                3.978392175157042,
                3.983853580922048,
                3.989314986687054,
                3.9947763924520596,
                4.000237798217066,
                4.005699203982071,
                4.011160609747077,
                4.016622015512083,
                4.022083421277089,
                4.027544827042095,
                4.033006232807101,
                4.038467638572107,
                4.043929044337113,
                4.049390450102119,
                4.054851855867125,
                4.06031326163213,
                4.065774667397136,
                4.071236073162142,
                4.076697478927148,
                4.082158884692154,
                4.08762029045716,
                4.08822029045716,
                4.08882029045716,
                4.08942029045716,
                4.09002029045716,
                4.09062029045716,
                4.0912202904571595,
                4.09182029045716,
                4.09242029045716,
                4.09302029045716,
                4.09362029045716,
                4.09422029045716,
                4.09482029045716,
                4.0954202904571595,
                4.09602029045716,
                4.09662029045716,
                4.09722029045716,
                4.09782029045716,
                4.09842029045716,
                4.09902029045716,
                4.0996202904571595,
                4.10022029045716,
                4.10082029045716,
                4.10142029045716,
                4.10202029045716,
                4.10262029045716,
                4.10322029045716,
                4.1038202904571595,
                4.10442029045716,
                4.10502029045716,
                4.10562029045716,
                4.10622029045716,
                4.10682029045716,
                4.10742029045716,
                4.108020290457159,
                4.10862029045716,
                4.10922029045716,
                4.10982029045716,
                4.11042029045716,
                4.1110202904571596,
                4.11162029045716,
                4.112220290457159,
                4.11282029045716,
                4.11342029045716,
                4.11402029045716,
                4.11462029045716,
                4.1152202904571595,
                4.11582029045716,
                4.116420290457159,
                4.11702029045716,
                4.117620290457159,
                4.11822029045716,
                4.11882029045716,
                4.1194202904571595,
                4.12002029045716,
                4.120620290457159,
                4.12122029045716,
                4.121820290457159,
                4.12242029045716,
                4.12302029045716,
                4.1236202904571595,
                4.12422029045716,
                4.124820290457159,
                4.12542029045716,
                4.126020290457159,
                4.12662029045716,
                4.12722029045716,
                4.1278202904571595,
                4.12842029045716,
                4.129020290457159,
                4.12962029045716,
                4.130220290457159,
                4.13082029045716,
                4.13142029045716,
                4.1320202904571595,
                4.13262029045716,
                4.133220290457159,
                4.13382029045716,
                4.134420290457159,
                4.13502029045716,
                4.13562029045716,
                4.136220290457159,
                4.13682029045716,
                4.137420290457159,
                4.13802029045716,
                4.138620290457159,
                4.13922029045716,
                4.13982029045716,
                4.140420290457159,
                4.14102029045716,
                4.141620290457159,
                4.14222029045716,
                4.142820290457159,
                4.1434202904571595,
                4.14402029045716,
                4.144620290457159,
                4.14522029045716,
                4.145820290457159,
                4.14642029045716,
                4.147020290457159,
            ],
        ),
        PlaneStressParam(
            linesolv=None,
            n_radial_array=100,
            nlayers=3,
            nu=numpy.array([0.3, 0.34006912702297704, 0.3]),
            rad=numpy.array(
                [
                    3.6732023601326333,
                    3.7688101124061717,
                    3.7649909451102674,
                    3.8249909451102675,
                ]
            ),
            ey=numpy.array([2.05000000e11, 21085960915.80571, 2.05000000e11]),
            j=numpy.array([0.00000000e00, -2245759961.294637, 0.00000000e00]),
            expected_sigr=[
                9.769733861957293e-08,
                -299296.0160031287,
                -598358.477675115,
                -897187.6279567069,
                -1195783.7094728944,
                -1494146.9645336915,
                -1792277.6351339405,
                -2090175.9629544844,
                -2387842.1893622647,
                -2685276.5554106147,
                -2982479.301840138,
                -3279450.669079198,
                -3576190.8972443077,
                -3872700.226140228,
                -4168978.895261434,
                -4465027.143791233,
                -4760845.210603429,
                -5056433.334262221,
                -5351791.7530228915,
                -5646920.704832389,
                -5941820.427329133,
                -6236491.157844482,
                -6530933.133402437,
                -6825146.590720718,
                -7119131.76621057,
                -7412888.89597764,
                -7706418.215822659,
                -7999719.961241548,
                -8292794.367425895,
                -8585641.669263843,
                -8878262.101339698,
                -9170655.897935387,
                -9462823.293030469,
                -9754764.52030252,
                -10046479.813127715,
                -10337969.404581523,
                -10629233.5274385,
                -10920272.414173374,
                -11211086.296961622,
                -11501675.407679182,
                -11792039.977903819,
                -12082180.238914639,
                -12372096.421693357,
                -12661788.756924588,
                -12951257.474995853,
                -13240502.80599855,
                -13529524.979728054,
                -13818324.225684403,
                -14106900.773072781,
                -14395254.850803819,
                -14683386.687493877,
                -14971296.511466429,
                -15258984.550750671,
                -15546451.033083988,
                -15833696.185911354,
                -16120720.23638582,
                -16407523.411369402,
                -16694105.937432678,
                -16980468.04085617,
                -17266609.94763033,
                -17552531.883456137,
                -17838234.073745485,
                -18123716.74362137,
                -18408980.11791859,
                -18694024.421184506,
                -18978849.877678875,
                -19263456.71137479,
                -19547845.145958323,
                -19832015.404830262,
                -20115967.711105447,
                -20399702.287613522,
                -20683219.356899556,
                -20966519.141224694,
                -21249601.862565402,
                -21532467.742615502,
                -21815117.002785694,
                -22097549.864204034,
                -22379766.547716625,
                -22661767.27388752,
                -22943552.26300009,
                -23225121.735056523,
                -23506475.90977892,
                -23787615.006609373,
                -24068539.244710173,
                -24349248.842964787,
                -24629744.019978136,
                -24910024.994076516,
                -25190091.98330896,
                -25469945.205446444,
                -25749584.87798367,
                -26029011.218138255,
                -26308224.44285203,
                -26587224.76879102,
                -26866012.412346132,
                -27144587.589632966,
                -27422950.51649318,
                -27701101.408494093,
                -27979040.48092967,
                -28256767.948820613,
                -28534284.026914664,
                -28811588.927430246,
                -28805999.455120668,
                -28791165.464775197,
                -28767086.858809426,
                -28733763.46373996,
                -28691195.10608341,
                -28639381.64488451,
                -28578322.939188015,
                -28508018.804667804,
                -28428469.111211345,
                -28339673.66364981,
                -28241632.364398792,
                -28134345.018289477,
                -28017811.46268118,
                -27892031.56746136,
                -27757005.191674754,
                -27612732.129309833,
                -27459212.261096764,
                -27296445.424394865,
                -27124431.456563458,
                -26943170.194961857,
                -26752661.50947752,
                -26552905.23746976,
                -26343901.2162979,
                -26125649.27247854,
                -25898149.26505643,
                -25661401.063919023,
                -25415404.4738975,
                -25160159.37572203,
                -24895665.574223787,
                -24621922.93929023,
                -24338931.30828068,
                -24046690.540239878,
                -23745200.450841717,
                -23434460.88828823,
                -23114471.73330958,
                -22785232.790736947,
                -22446743.930457793,
                -22099004.97898872,
                -21742015.79537447,
                -21375776.206131652,
                -21000286.070305,
                -20615545.225253846,
                -20221553.530022923,
                -19818310.789443415,
                -19405816.884245485,
                -18984071.640945747,
                -18553074.92943165,
                -18112826.55453438,
                -17663326.396984097,
                -17204574.294140123,
                -16736570.061676355,
                -16259313.558637533,
                -15772804.65491111,
                -15277043.166170983,
                -14772028.929776467,
                -14257761.81561502,
                -13734241.661045957,
                -13201468.3034286,
                -12659441.590964975,
                -12108161.361014403,
                -11547627.450936202,
                -10977839.730617827,
                -10398798.037418595,
                -9810502.208697824,
                -9212952.070972122,
                -8606147.494128942,
                -7990088.326370316,
                -7364774.372527423,
                -6730205.513330433,
                -6086381.586138662,
                -5433302.417468717,
                -4770967.877208053,
                -4099377.8027159884,
                -3418532.0313518397,
                -2728430.389632215,
                -2029072.7582872796,
                -1320458.9529909282,
                -602588.8436306156,
                124537.7541197646,
                860920.9812154698,
                1606561.0111398937,
                2361457.974005582,
                3125612.021610504,
                3899023.349123478,
                4681692.064971625,
                5473618.364323763,
                6274802.377292438,
                7085244.277361044,
                7904944.216327549,
                8733902.345989924,
                9572118.818146138,
                10419593.806279587,
                11276327.451345524,
                12142319.905141924,
                13017571.319466753,
                13902081.88948883,
                14795851.734477988,
                15698881.02791762,
                16611169.910762986,
                17532718.54152625,
                17353232.12053368,
                17173831.47581403,
                16994516.55271951,
                16815287.296645615,
                16636143.653031997,
                16457085.567361193,
                16278112.985159123,
                16099225.851995375,
                15920424.113482526,
                15741707.71527633,
                15563076.603076415,
                15384530.722625,
                15206070.019707484,
                15027694.440152455,
                14849403.929831775,
                14671198.434660003,
                14493077.900594682,
                14315042.273636637,
                14137091.499829292,
                13959225.525259051,
                13781444.29605492,
                13603747.758388791,
                13426135.858475542,
                13248608.542572359,
                13071165.756979309,
                12893807.448039062,
                12716533.562136687,
                12539344.045699945,
                12362238.845199,
                12185217.907146512,
                12008281.178097446,
                11831428.60464926,
                11654660.133441824,
                11477975.711157104,
                11301375.28451918,
                11124858.800294822,
                10948426.205292325,
                10772077.44636287,
                10595812.470398966,
                10419631.224335717,
                10243533.655149944,
                10067519.70986048,
                9891589.335528163,
                9715742.479255551,
                9539979.088187408,
                9364299.109509919,
                9188702.490451185,
                9013189.178281117,
                8837759.120311053,
                8662412.26389434,
                8487148.556425454,
                8311967.945341077,
                8136870.378118926,
                7961855.802278334,
                7786924.165380256,
                7612075.415026678,
                7437309.498861204,
                7262626.3645690605,
                7088025.959876113,
                6913508.232550042,
                6739073.130399465,
                6564720.601274324,
                6390450.593065496,
                6216263.053704988,
                6042157.931166332,
                5868135.173463404,
                5694194.728651507,
                5520336.544826585,
                5346560.57012581,
                5172866.752726997,
                4999255.040848799,
                4825725.382750801,
                4652277.726733331,
                4478912.021137162,
                4305628.214344004,
                4132426.2547761104,
                3959306.0908961794,
                3786267.671208043,
                3613310.9442552943,
                3440435.858622461,
                3267642.3629344213,
                3094930.4058564017,
                2922299.9360940754,
                2749750.902393269,
                2577283.2535404516,
                2404896.9383620503,
                2232591.9057247434,
                2060368.1045355585,
                1888225.483741481,
                1716163.9923295528,
                1544183.579326871,
                1372284.1938010778,
                1200465.7848589912,
                1028728.3016480722,
                857071.693355251,
                685495.9092075138,
                514000.8984715119,
                342586.610454148,
                171252.99450159894,
            ],
            expected_sigt=[
                -1150329392.874935,
                -1150030096.8589315,
                -1149731034.3972595,
                -1149432205.2469778,
                -1149133609.1654618,
                -1148835245.9104009,
                -1148537115.2398007,
                -1148239216.9119802,
                -1147941550.6855724,
                -1147644116.3195238,
                -1147346913.5730944,
                -1147049942.2058554,
                -1146753201.9776905,
                -1146456692.6487944,
                -1146160413.9796731,
                -1145864365.7311432,
                -1145568547.6643312,
                -1145272959.5406723,
                -1144977601.1219118,
                -1144682472.1701021,
                -1144387572.4476054,
                -1144092901.7170901,
                -1143798459.741532,
                -1143504246.284214,
                -1143210261.108724,
                -1142916503.978957,
                -1142622974.659112,
                -1142329672.9136932,
                -1142036598.5075085,
                -1141743751.2056708,
                -1141451130.7735949,
                -1141158736.976999,
                -1140866569.5819042,
                -1140574628.3546321,
                -1140282913.061807,
                -1139991423.4703531,
                -1139700159.347496,
                -1139409120.4607613,
                -1139118306.577973,
                -1138827717.4672556,
                -1138537352.8970308,
                -1138247212.6360202,
                -1137957296.4532413,
                -1137667604.11801,
                -1137378135.3999386,
                -1137088890.068936,
                -1136799867.8952067,
                -1136511068.6492503,
                -1136222492.1018617,
                -1135934138.0241308,
                -1135646006.1874406,
                -1135358096.3634682,
                -1135070408.324184,
                -1134782941.8418505,
                -1134495696.6890233,
                -1134208672.6385489,
                -1133921869.463565,
                -1133635286.937502,
                -1133348924.8340786,
                -1133062782.9273043,
                -1132776860.9914784,
                -1132491158.8011892,
                -1132205676.131313,
                -1131920412.757016,
                -1131635368.4537501,
                -1131350542.9972558,
                -1131065936.16356,
                -1130781547.7289762,
                -1130497377.4701045,
                -1130213425.163829,
                -1129929690.587321,
                -1129646173.5180352,
                -1129362873.73371,
                -1129079791.0123694,
                -1128796925.132319,
                -1128514275.8721488,
                -1128231843.0107305,
                -1127949626.3272178,
                -1127667625.601047,
                -1127385840.6119344,
                -1127104271.139878,
                -1126822916.9651556,
                -1126541777.8683252,
                -1126260853.6302245,
                -1125980144.0319698,
                -1125699648.8549564,
                -1125419367.8808582,
                -1125139300.8916256,
                -1124859447.6694882,
                -1124579807.996951,
                -1124300381.6567965,
                -1124021168.4320827,
                -1123742168.1061435,
                -1123463380.4625885,
                -1123184805.2853017,
                -1122906442.3584416,
                -1122628291.4664407,
                -1122350352.3940048,
                -1122072624.926114,
                -1121795108.8480198,
                -124266344.44614087,
                -124265739.92156953,
                -124261991.7532102,
                -124255099.97359103,
                -124245064.60981879,
                -124231885.69984297,
                -124215563.26534902,
                -124196097.34970774,
                -124173487.95834053,
                -124147735.14546093,
                -124118838.92733301,
                -124086799.34190626,
                -124051616.41628748,
                -124013290.17216207,
                -123971820.65290089,
                -123927207.88018936,
                -123879451.89739834,
                -123828552.71537052,
                -123774510.37205541,
                -123717324.89455979,
                -123656996.32083315,
                -123593524.67798226,
                -123526909.9931139,
                -123457152.29875623,
                -123384251.61659464,
                -123308207.98457865,
                -123229021.4352364,
                -123146692.00651737,
                -123061219.69842154,
                -122972604.56516251,
                -122880846.63926838,
                -122785945.9424246,
                -122687902.49631658,
                -122586716.34431519,
                -122482387.50810583,
                -122374916.02563798,
                -122264301.92401847,
                -122150545.23035404,
                -122033645.97717284,
                -121913604.20242435,
                -121790419.91152994,
                -121664093.14786047,
                -121534623.9602081,
                -121402012.36483695,
                -121266258.37801106,
                -121127362.05394399,
                -120985323.3980571,
                -120840142.45914261,
                -120691819.25346456,
                -120540353.83523655,
                -120385746.19903718,
                -120227996.41534409,
                -120067104.46789323,
                -119903070.42716222,
                -119735894.30399376,
                -119565576.13633737,
                -119392115.9512998,
                -119215513.7705665,
                -119035769.62666558,
                -118852883.56296791,
                -118666855.61200161,
                -118477685.7737667,
                -118285374.10789809,
                -118089920.63065982,
                -117891325.38000143,
                -117689588.3776083,
                -117484709.67227267,
                -117276689.26941587,
                -117065527.20698741,
                -116851223.52293678,
                -116633778.23894939,
                -116413191.39297475,
                -116189463.00127691,
                -115962593.11806946,
                -115732581.74877371,
                -115499428.93676054,
                -115263134.71455811,
                -115023699.09300908,
                -114781122.12632704,
                -114535403.83077605,
                -114286544.24972697,
                -114034543.3940225,
                -113779401.30703348,
                -113521118.02670942,
                -113259693.55305032,
                -112995127.94569108,
                -112727421.23173851,
                -112456573.42745665,
                -112182584.57079498,
                -111905454.67801762,
                -111625183.80875945,
                -111341771.97928454,
                -111055219.20585696,
                -110765525.53184757,
                -110472690.99520583,
                -110176715.61219586,
                -109877599.42618847,
                -109575342.43718366,
                -109269944.70481637,
                -108961406.26703608,
                -1109011101.25244,
                -1108831614.8314474,
                -1108652214.1867275,
                -1108472899.263633,
                -1108293670.0075593,
                -1108114526.3639457,
                -1107935468.2782748,
                -1107756495.6960728,
                -1107577608.5629091,
                -1107398806.8243961,
                -1107220090.42619,
                -1107041459.31399,
                -1106862913.4335387,
                -1106684452.730621,
                -1106506077.151066,
                -1106327786.6407454,
                -1106149581.1455736,
                -1105971460.6115084,
                -1105793424.9845502,
                -1105615474.210743,
                -1105437608.2361727,
                -1105259827.0069685,
                -1105082130.4693024,
                -1104904518.569389,
                -1104726991.253486,
                -1104549548.467893,
                -1104372190.1589527,
                -1104194916.2730505,
                -1104017726.7566135,
                -1103840621.5561128,
                -1103663600.61806,
                -1103486663.8890111,
                -1103309811.315563,
                -1103133042.8443556,
                -1102956358.4220707,
                -1102779757.9954326,
                -1102603241.5112083,
                -1102426808.916206,
                -1102250460.1572766,
                -1102074195.1813126,
                -1101898013.9352493,
                -1101721916.3660636,
                -1101545902.420774,
                -1101369972.0464418,
                -1101194125.190169,
                -1101018361.7991009,
                -1100842681.8204236,
                -1100667085.2013648,
                -1100491571.8891947,
                -1100316141.8312247,
                -1100140794.974808,
                -1099965531.267339,
                -1099790350.6562548,
                -1099615253.0890326,
                -1099440238.5131922,
                -1099265306.876294,
                -1099090458.1259403,
                -1098915692.209775,
                -1098741009.0754828,
                -1098566408.6707897,
                -1098391890.9434638,
                -1098217455.8413131,
                -1098043103.312188,
                -1097868833.3039792,
                -1097694645.7646186,
                -1097520540.64208,
                -1097346517.8843772,
                -1097172577.4395652,
                -1096998719.2557402,
                -1096824943.2810395,
                -1096651249.4636407,
                -1096477637.7517624,
                -1096304108.0936644,
                -1096130660.4376469,
                -1095957294.732051,
                -1095784010.9252577,
                -1095610808.9656897,
                -1095437688.8018098,
                -1095264650.3821218,
                -1095091693.6551688,
                -1094918818.5695362,
                -1094746025.0738482,
                -1094573313.11677,
                -1094400682.6470077,
                -1094228133.613307,
                -1094055665.9644542,
                -1093883279.6492758,
                -1093710974.6166384,
                -1093538750.8154492,
                -1093366608.1946552,
                -1093194546.7032433,
                -1093022566.2902405,
                -1092850666.9047146,
                -1092678848.4957728,
                -1092507111.0125618,
                -1092335454.404269,
                -1092163878.620121,
                -1091992383.609385,
                -1091820969.3213677,
                -1091649635.7054152,
            ],
            expected_r_deflect=[
                -0.02061167141872268,
                -0.02061006285676633,
                -0.020608456108712828,
                -0.020606851173146885,
                -0.02060524804865469,
                -0.0206036467338239,
                -0.020602047227243632,
                -0.020600449527504477,
                -0.02059885363319849,
                -0.02059725954291919,
                -0.020595667255261546,
                -0.020594076768822,
                -0.020592488082198438,
                -0.020590901193990215,
                -0.020589316102798128,
                -0.020587732807224424,
                -0.020586151305872807,
                -0.020584571597348424,
                -0.020582993680257868,
                -0.020581417553209174,
                -0.02057984321481182,
                -0.02057827066367673,
                -0.020576699898416255,
                -0.020575130917644187,
                -0.020573563719975752,
                -0.020571998304027617,
                -0.02057043466841786,
                -0.020568872811766007,
                -0.020567312732693,
                -0.020565754429821212,
                -0.020564197901774432,
                -0.02056264314717788,
                -0.020561090164658184,
                -0.0205595389528434,
                -0.020557989510363,
                -0.020556441835847858,
                -0.02055489592793027,
                -0.020553351785243944,
                -0.020551809406423985,
                -0.020550268790106924,
                -0.02054872993493068,
                -0.02054719283953458,
                -0.020545657502559354,
                -0.02054412392264713,
                -0.020542592098441434,
                -0.02054106202858719,
                -0.020539533711730715,
                -0.020538007146519717,
                -0.020536482331603295,
                -0.02053495926563194,
                -0.020533437947257517,
                -0.020531918375133296,
                -0.02053040054791392,
                -0.020528884464255416,
                -0.02052737012281518,
                -0.020525857522252007,
                -0.020524346661226045,
                -0.020522837538398835,
                -0.020521330152433286,
                -0.02051982450199367,
                -0.02051832058574564,
                -0.020516818402356207,
                -0.02051531795049375,
                -0.02051381922882802,
                -0.020512322236030122,
                -0.02051082697077252,
                -0.020509333431729046,
                -0.020507841617574876,
                -0.020506351526986554,
                -0.02050486315864197,
                -0.02050337651122037,
                -0.02050189158340234,
                -0.02050040837386984,
                -0.020498926881306147,
                -0.020497447104395897,
                -0.020495969041825067,
                -0.02049449269228098,
                -0.02049301805445229,
                -0.020491545127029,
                -0.020490073908702437,
                -0.020488604398165272,
                -0.02048713659411151,
                -0.020485670495236478,
                -0.02048420610023684,
                -0.020482743407810592,
                -0.02048128241665704,
                -0.02047982312547683,
                -0.02047836553297193,
                -0.020476909637845615,
                -0.020475455438802498,
                -0.020474002934548496,
                -0.02047255212379085,
                -0.020471103005238116,
                -0.02046965557760015,
                -0.02046820983958814,
                -0.020466765789914566,
                -0.020465323427293222,
                -0.02046388275043921,
                -0.02046244375806893,
                -0.020461006448900098,
                -0.020459570823732065,
                -0.020459595184547652,
                -0.020459619563553133,
                -0.020459643971662445,
                -0.02045966842888447,
                -0.020459692946133146,
                -0.02045971754068887,
                -0.020459742226194066,
                -0.020459767017200647,
                -0.020459791929170024,
                -0.02045981697665411,
                -0.020459842176023812,
                -0.02045986753910256,
                -0.020459893083170755,
                -0.02045991882096132,
                -0.020459944768845162,
                -0.020459970941374195,
                -0.02045999735400983,
                -0.02046002401948499,
                -0.020460050955080078,
                -0.020460078173528018,
                -0.02046010569029022,
                -0.020460133522647084,
                -0.020460161681512545,
                -0.02046019018325751,
                -0.020460219041524397,
                -0.0204602482745031,
                -0.020460277894017054,
                -0.020460307915527665,
                -0.020460338355405838,
                -0.0204603692263845,
                -0.020460400543015567,
                -0.020460432322579436,
                -0.02046046457871853,
                -0.02046049732416577,
                -0.020460530577111058,
                -0.020460564350287314,
                -0.020460598660065443,
                -0.020460633518268878,
                -0.0204606689449065,
                -0.020460704948163766,
                -0.020460741548049555,
                -0.020460778757296794,
                -0.020460816590457398,
                -0.020460855064811767,
                -0.020460894191273837,
                -0.020460933988033503,
                -0.020460974468733184,
                -0.0204610156470153,
                -0.020461057538341265,
                -0.02046110015908198,
                -0.020461143522879865,
                -0.020461187644286838,
                -0.0204612325396738,
                -0.020461278219954693,
                -0.0204613247042289,
                -0.020461372006138845,
                -0.020461420139326947,
                -0.020461469119254616,
                -0.020461518962292757,
                -0.020461569679355307,
                -0.020461621289541654,
                -0.020461673805584724,
                -0.02046172724203643,
                -0.020461781615267682,
                -0.020461836938011402,
                -0.020461893226638495,
                -0.020461950497519865,
                -0.020462008760659955,
                -0.020462068035158154,
                -0.020462128333747387,
                -0.02046218967097957,
                -0.020462252065044595,
                -0.0204623155268564,
                -0.020462380073695385,
                -0.020462445718294475,
                -0.020462512477024575,
                -0.020462580365347094,
                -0.020462649395085464,
                -0.02046271958442958,
                -0.02046279094702186,
                -0.02046286349741422,
                -0.020462937250158575,
                -0.020463012219806842,
                -0.020463088422729925,
                -0.020463165872570244,
                -0.02046324458478921,
                -0.02046332457302924,
                -0.020463405854570738,
                -0.020463488441237132,
                -0.020463572350308823,
                -0.02046365759542823,
                -0.020463744190237776,
                -0.020463832152927353,
                -0.020463921495320392,
                -0.0204640122337878,
                -0.020464104381971993,
                -0.020464197955334384,
                -0.02046429296842689,
                -0.020464389438529906,
                -0.020464487376623448,
                -0.02046346263699384,
                -0.020462438580124982,
                -0.020461415205690593,
                -0.02046039251336463,
                -0.020459370502821237,
                -0.020458349173734777,
                -0.020457328525779812,
                -0.020456308558631123,
                -0.02045528927196368,
                -0.020454270665452676,
                -0.020453252738773502,
                -0.020452235491601762,
                -0.020451218923613262,
                -0.020450203034484013,
                -0.020449187823890234,
                -0.020448173291508352,
                -0.020447159437014997,
                -0.020446146260087002,
                -0.020445133760401413,
                -0.020444121937635475,
                -0.020443110791466636,
                -0.020442100321572558,
                -0.0204410905276311,
                -0.02044008140932033,
                -0.020439072966318514,
                -0.020438065198304132,
                -0.020437058104955858,
                -0.020436051685952576,
                -0.020435045940973375,
                -0.020434040869697537,
                -0.020433036471804562,
                -0.020432032746974152,
                -0.020431029694886197,
                -0.0204300273152208,
                -0.020429025607658274,
                -0.02042802457187912,
                -0.02042702420756405,
                -0.02042602451439398,
                -0.02042502549205002,
                -0.020424027140213494,
                -0.020423029458565917,
                -0.020422032446789005,
                -0.020421036104564687,
                -0.02042004043157508,
                -0.02041904542750252,
                -0.02041805109202952,
                -0.020417057424838812,
                -0.02041606442561332,
                -0.020415072094036174,
                -0.020414080429790704,
                -0.020413089432560438,
                -0.0204120991020291,
                -0.020411109437880624,
                -0.02041012043979913,
                -0.020409132107468954,
                -0.02040814444057462,
                -0.02040715743880086,
                -0.020406171101832585,
                -0.020405185429354933,
                -0.020404200421053224,
                -0.020403216076612975,
                -0.020402232395719914,
                -0.020401249378059955,
                -0.020400267023319218,
                -0.020399285331184014,
                -0.02039830430134086,
                -0.02039732393347646,
                -0.020396344227277735,
                -0.020395365182431777,
                -0.0203943867986259,
                -0.020393409075547594,
                -0.020392432012884563,
                -0.0203914556103247,
                -0.020390479867556092,
                -0.02038950478426703,
                -0.02038853036014599,
                -0.020387556594881662,
                -0.020386583488162913,
                -0.020385611039678817,
                -0.020384639249118642,
                -0.020383668116171847,
                -0.020382697640528094,
                -0.020381727821877235,
                -0.020380758659909314,
                -0.02037979015431458,
                -0.02037882230478347,
                -0.020377855111006614,
                -0.02037688857267484,
                -0.020375922689479173,
                -0.02037495746111082,
                -0.020373992887261196,
                -0.020373028967621908,
                -0.02037206570188475,
                -0.02037110308974171,
                -0.020370141130884974,
                -0.020369179825006918,
                -0.020368219171800116,
                -0.02036725917095733,
                -0.020366299822171516,
            ],
            expected_rradius=[
                3.6732023601326333,
                3.6741584376553686,
                3.675114515178104,
                3.6760705927008397,
                3.677026670223575,
                3.6779827477463103,
                3.6789388252690456,
                3.679894902791781,
                3.6808509803145166,
                3.681807057837252,
                3.6827631353599872,
                3.6837192128827225,
                3.684675290405458,
                3.685631367928193,
                3.686587445450929,
                3.687543522973664,
                3.6884996004963995,
                3.6894556780191348,
                3.69041175554187,
                3.691367833064606,
                3.692323910587341,
                3.6932799881100764,
                3.6942360656328117,
                3.695192143155547,
                3.6961482206782827,
                3.697104298201018,
                3.6980603757237533,
                3.6990164532464886,
                3.699972530769224,
                3.7009286082919597,
                3.701884685814695,
                3.7028407633374303,
                3.7037968408601656,
                3.704752918382901,
                3.705708995905636,
                3.706665073428372,
                3.7076211509511072,
                3.7085772284738425,
                3.709533305996578,
                3.710489383519313,
                3.711445461042049,
                3.712401538564784,
                3.7133576160875195,
                3.7143136936102548,
                3.71526977113299,
                3.716225848655726,
                3.717181926178461,
                3.7181380037011964,
                3.7190940812239317,
                3.720050158746667,
                3.7210062362694023,
                3.721962313792138,
                3.7229183913148733,
                3.7238744688376086,
                3.724830546360344,
                3.7257866238830792,
                3.726742701405815,
                3.7276987789285503,
                3.7286548564512856,
                3.729610933974021,
                3.730567011496756,
                3.731523089019492,
                3.7324791665422272,
                3.7334352440649625,
                3.734391321587698,
                3.735347399110433,
                3.736303476633169,
                3.737259554155904,
                3.7382156316786395,
                3.7391717092013748,
                3.74012778672411,
                3.7410838642468454,
                3.742039941769581,
                3.7429960192923164,
                3.7439520968150517,
                3.744908174337787,
                3.7458642518605223,
                3.746820329383258,
                3.7477764069059933,
                3.7487324844287286,
                3.749688561951464,
                3.7506446394741992,
                3.751600716996935,
                3.7525567945196703,
                3.7535128720424056,
                3.754468949565141,
                3.755425027087876,
                3.756381104610612,
                3.757337182133347,
                3.7582932596560825,
                3.759249337178818,
                3.760205414701553,
                3.7611614922242884,
                3.762117569747024,
                3.7630736472697595,
                3.7640297247924948,
                3.76498580231523,
                3.7659418798379654,
                3.766897957360701,
                3.7678540348834364,
                3.7688101124061717,
                3.7687719207332124,
                3.7687337290602536,
                3.768695537387295,
                3.7686573457143355,
                3.7686191540413763,
                3.7685809623684174,
                3.7685427706954586,
                3.7685045790224994,
                3.76846638734954,
                3.7684281956765813,
                3.7683900040036225,
                3.768351812330663,
                3.768313620657704,
                3.768275428984745,
                3.7682372373117863,
                3.768199045638827,
                3.7681608539658678,
                3.768122662292909,
                3.76808447061995,
                3.768046278946991,
                3.7680080872740316,
                3.7679698956010728,
                3.767931703928114,
                3.7678935122551547,
                3.7678553205821954,
                3.7678171289092366,
                3.7677789372362778,
                3.7677407455633185,
                3.7677025538903592,
                3.7676643622174004,
                3.7676261705444416,
                3.7675879788714823,
                3.767549787198523,
                3.7675115955255642,
                3.7674734038526054,
                3.767435212179646,
                3.767397020506687,
                3.767358828833728,
                3.7673206371607693,
                3.76728244548781,
                3.7672442538148507,
                3.767206062141892,
                3.767167870468933,
                3.767129678795974,
                3.7670914871230146,
                3.7670532954500557,
                3.767015103777097,
                3.7669769121041377,
                3.7669387204311784,
                3.7669005287582196,
                3.7668623370852607,
                3.7668241454123015,
                3.766785953739342,
                3.7667477620663834,
                3.7667095703934246,
                3.7666713787204653,
                3.766633187047506,
                3.7665949953745472,
                3.7665568037015884,
                3.766518612028629,
                3.76648042035567,
                3.766442228682711,
                3.7664040370097522,
                3.766365845336793,
                3.7663276536638337,
                3.766289461990875,
                3.766251270317916,
                3.766213078644957,
                3.7661748869719975,
                3.7661366952990387,
                3.76609850362608,
                3.7660603119531206,
                3.7660221202801614,
                3.7659839286072025,
                3.7659457369342437,
                3.7659075452612845,
                3.765869353588325,
                3.7658311619153664,
                3.7657929702424076,
                3.7657547785694483,
                3.765716586896489,
                3.76567839522353,
                3.7656402035505714,
                3.765602011877612,
                3.765563820204653,
                3.765525628531694,
                3.765487436858735,
                3.765449245185776,
                3.7654110535128167,
                3.765372861839858,
                3.765334670166899,
                3.76529647849394,
                3.7652582868209805,
                3.7652200951480217,
                3.765181903475063,
                3.7651437118021036,
                3.7651055201291443,
                3.7650673284561855,
                3.7650291367832267,
                3.7649909451102674,
                3.7655909451102674,
                3.7661909451102673,
                3.7667909451102672,
                3.7673909451102676,
                3.7679909451102676,
                3.7685909451102675,
                3.7691909451102674,
                3.7697909451102674,
                3.7703909451102673,
                3.7709909451102677,
                3.7715909451102676,
                3.7721909451102675,
                3.7727909451102675,
                3.7733909451102674,
                3.7739909451102673,
                3.7745909451102673,
                3.7751909451102676,
                3.7757909451102676,
                3.7763909451102675,
                3.7769909451102675,
                3.7775909451102674,
                3.7781909451102673,
                3.7787909451102673,
                3.7793909451102676,
                3.7799909451102676,
                3.7805909451102675,
                3.7811909451102674,
                3.7817909451102674,
                3.7823909451102673,
                3.7829909451102672,
                3.7835909451102676,
                3.7841909451102675,
                3.7847909451102675,
                3.7853909451102674,
                3.7859909451102673,
                3.7865909451102673,
                3.7871909451102677,
                3.7877909451102676,
                3.7883909451102675,
                3.7889909451102675,
                3.7895909451102674,
                3.7901909451102673,
                3.7907909451102673,
                3.7913909451102676,
                3.7919909451102676,
                3.7925909451102675,
                3.7931909451102674,
                3.7937909451102674,
                3.7943909451102673,
                3.7949909451102672,
                3.7955909451102676,
                3.7961909451102676,
                3.7967909451102675,
                3.7973909451102674,
                3.7979909451102674,
                3.7985909451102673,
                3.7991909451102677,
                3.7997909451102676,
                3.8003909451102675,
                3.8009909451102675,
                3.8015909451102674,
                3.8021909451102673,
                3.8027909451102673,
                3.8033909451102677,
                3.8039909451102676,
                3.8045909451102675,
                3.8051909451102675,
                3.8057909451102674,
                3.8063909451102673,
                3.8069909451102673,
                3.8075909451102676,
                3.8081909451102676,
                3.8087909451102675,
                3.8093909451102674,
                3.8099909451102674,
                3.8105909451102673,
                3.8111909451102677,
                3.8117909451102676,
                3.8123909451102675,
                3.8129909451102675,
                3.8135909451102674,
                3.8141909451102674,
                3.8147909451102673,
                3.8153909451102677,
                3.8159909451102676,
                3.8165909451102675,
                3.8171909451102675,
                3.8177909451102674,
                3.8183909451102673,
                3.8189909451102673,
                3.8195909451102676,
                3.8201909451102676,
                3.8207909451102675,
                3.8213909451102674,
                3.8219909451102674,
                3.8225909451102673,
                3.8231909451102677,
                3.8237909451102676,
                3.8243909451102676,
            ],
        ),
    ),
)
def test_plane_stress(planestressparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for plane_stress.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param planestressparam: the data used to mock and assert in this test.
    :type planestressparam: planestressparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    sigr, sigt, r_deflect, rradius = sctf.plane_stress(
        n_radial_array=planestressparam.n_radial_array,
        nlayers=planestressparam.nlayers,
        nu=planestressparam.nu,
        rad=planestressparam.rad,
        ey=planestressparam.ey,
        j=planestressparam.j,
    )

    assert sigr == pytest.approx(planestressparam.expected_sigr)
    assert sigt == pytest.approx(planestressparam.expected_sigt)
    assert r_deflect == pytest.approx(planestressparam.expected_r_deflect)
    assert rradius == pytest.approx(planestressparam.expected_rradius)


class ExtendedPlaneStrainParam(NamedTuple):
    n_radial_array: Any = None

    nlayers: Any = None

    i_tf_bucking: Any = None

    nu_t: Any = None

    nu_zt: Any = None

    ey_t: Any = None

    ey_z: Any = None

    d_curr: Any = None

    rad: Any = None

    v_force: Any = None

    expected_sigr: Any = None

    expected_sigt: Any = None

    expected_sigz: Any = None

    expected_str_r: Any = None

    expected_str_t: Any = None

    expected_str_z: Any = None

    expected_r_deflect: Any = None

    expected_rradius: Any = None


@pytest.mark.parametrize(
    "extendedplanestrainparam",
    (
        ExtendedPlaneStrainParam(
            n_radial_array=100,
            nlayers=2,
            i_tf_bucking=0,
            nu_t=numpy.array(
                numpy.array((0.34999999999999998, 0.29999999999999999), order="F"),
                order="F",
            ).transpose(),
            nu_zt=numpy.array(
                numpy.array((0.34948024015688117, 0.29999999999999999), order="F"),
                order="F",
            ).transpose(),
            ey_t=numpy.array(
                numpy.array((117000000000, 205000000000), order="F"), order="F"
            ).transpose(),
            ey_z=numpy.array(
                numpy.array((97843910970.178864, 205000000000.00003), order="F"),
                order="F",
            ).transpose(),
            d_curr=numpy.array(
                numpy.array((375174117.44492483, 0), order="F"), order="F"
            ).transpose(),
            rad=numpy.array(
                numpy.array((0, 0.14708850000000001, 0.15483000000000002), order="F"),
                order="F",
            ).transpose(),
            v_force=147221629.66130698,
            expected_sigr=[
                numpy.nan,
                -829140939.1329328,
                -828878786.3701528,
                -828441865.0984682,
                -827830175.3181431,
                -827043717.0291449,
                -826082490.231469,
                -824946494.925121,
                -823635731.1101007,
                -822150198.7864125,
                -820489897.9540532,
                -818654828.6130259,
                -816644990.7633287,
                -814460384.4049616,
                -812101009.5379249,
                -809566866.1622195,
                -806857954.2778445,
                -803974273.8847997,
                -800915824.983086,
                -797682607.5727025,
                -794274621.6536493,
                -790691867.2259276,
                -786934344.2895358,
                -783002052.8444749,
                -778894992.890745,
                -774613164.4283451,
                -770156567.4572762,
                -765525201.9775376,
                -760719067.9891297,
                -755738165.4920528,
                -750582494.4863062,
                -745252054.9718906,
                -739746846.9488052,
                -734066870.4170501,
                -728212125.3766268,
                -722182611.8275331,
                -715978329.7697704,
                -709599279.2033384,
                -703045460.1282363,
                -696316872.5444654,
                -689413516.4520254,
                -682335391.8509156,
                -675082498.7411368,
                -667654837.1226883,
                -660052406.9955704,
                -652275208.3597835,
                -644323241.215327,
                -636196505.5622008,
                -627895001.4004059,
                -619418728.729941,
                -610767687.550807,
                -601941877.8630037,
                -592941299.6665307,
                -583765952.9613888,
                -574415837.7475771,
                -564890954.025096,
                -555191301.793946,
                -545316881.054126,
                -535267691.80563694,
                -525043734.0484789,
                -514645007.7826511,
                -504071513.0081537,
                -493323249.72498757,
                -482400217.9331514,
                -471302417.63264626,
                -460029848.8234718,
                -448582511.5056272,
                -436960405.6791144,
                -425163531.34393185,
                -413191888.5000792,
                -401045477.1475577,
                -388724297.2863667,
                -376228348.91650677,
                -363557632.03797716,
                -350712146.6507778,
                -337691892.7549095,
                -324496870.3503715,
                -311127079.4371644,
                -297582520.0152877,
                -283863192.084742,
                -269969095.6455273,
                -255900230.69764256,
                -241656597.24108845,
                -227238195.27586463,
                -212645024.80197182,
                -197877085.81941003,
                -182934378.32817852,
                -167816902.32827768,
                -152524657.81970748,
                -137057644.80246758,
                -121415863.2765587,
                -105599313.2419801,
                -89607994.69873217,
                -73441907.64681524,
                -57101052.08622896,
                -40585428.01697298,
                -23895035.439047653,
                -7029874.352452978,
                10010055.242811045,
                27224753.346744414,
                27224753.346743274,
                26928096.53957756,
                26631912.241218276,
                26336199.448728982,
                26040957.16183203,
                25746184.3829009,
                25451880.116952077,
                25158043.371635936,
                24864673.157229107,
                24571768.48662679,
                24279328.375331763,
                23987351.841448147,
                23695837.905674234,
                23404785.59129195,
                23114193.9241592,
                22824061.932703156,
                22534388.647910696,
                22245173.10331978,
                21956414.335014176,
                21668111.381612953,
                21380263.284262314,
                21092869.08662891,
                20805927.834891703,
                20519438.577733815,
                20233400.366334405,
                19947812.254360523,
                19662673.29796089,
                19377982.555756312,
                19093739.08883156,
                18809941.960731525,
                18526590.23744783,
                18243682.98741451,
                17961219.281500816,
                17679198.19300216,
                17397618.797632422,
                17116480.173517253,
                16835781.40118595,
                16555521.56356474,
                16275699.745967694,
                15996315.03609098,
                15717366.524004238,
                15438853.302144371,
                15160774.465306437,
                14883129.110638388,
                14605916.337631976,
                14329135.248115571,
                14052784.94624889,
                13776864.538512478,
                13501373.133703861,
                13226309.842928467,
                12951673.779591953,
                12677464.05939495,
                12403679.800323,
                12130320.122644654,
                11857384.148897571,
                11584871.0038871,
                11312779.814676696,
                11041109.71058122,
                10769859.823160235,
                10499029.286210349,
                10228617.23576043,
                9958622.810062023,
                9689045.14958362,
                9419883.397004897,
                9151136.697207155,
                8882804.197269963,
                8614885.04646254,
                8347378.396235619,
                8080283.400217615,
                7813599.214206012,
                7547324.996161614,
                7281459.906201365,
                7016003.106592132,
                6750953.761742556,
                6486311.038200666,
                6222074.104640953,
                5958242.131863889,
                5694814.292785397,
                5431789.762433497,
                5169167.717938738,
                4906947.338528923,
                4645127.80552433,
                4383708.30232957,
                4122688.0144268856,
                3862066.129371364,
                3601841.836783756,
                3342014.3283437733,
                3082582.797784345,
                2823546.4408853957,
                2564904.455468096,
                2306656.041387208,
                2048800.4005258176,
                1791336.7367895895,
                1534264.2560995868,
                1277582.1663874835,
                1021289.6775879061,
                765386.0016326876,
                509870.35244656017,
                254741.94593901612,
                -4.787169592359073e-07,
            ],
            expected_sigt=[
                numpy.nan,
                -829161394.3774397,
                -828960607.3461856,
                -828625962.2944794,
                -828157459.2220584,
                -827555098.1289548,
                -826818879.0151731,
                -825948801.8807079,
                -824944866.7255567,
                -823807073.5497193,
                -822535422.3531967,
                -821129913.1359856,
                -819590545.89809,
                -817917320.6395062,
                -816110237.3602372,
                -814169296.0602809,
                -812094496.7396379,
                -809885839.3983089,
                -807543324.0362918,
                -805066950.6535902,
                -802456719.2502012,
                -799712629.8261249,
                -796834682.3813626,
                -793822876.9159135,
                -790677213.4297779,
                -787397691.9229555,
                -783984312.3954465,
                -780437074.8472508,
                -776755979.2783681,
                -772941025.6887994,
                -768992214.0785445,
                -764909544.4476022,
                -760693016.7959738,
                -756342631.1236581,
                -751858387.4306564,
                -747240285.7169673,
                -742488325.9825923,
                -737602508.2275307,
                -732582832.451782,
                -727429298.655347,
                -722141906.8382254,
                -716720657.0004174,
                -711165549.141922,
                -705476583.2627406,
                -699653759.3628722,
                -693697077.442318,
                -687606537.5010763,
                -681382139.539148,
                -675023883.5565336,
                -668531769.5532327,
                -661905797.5292447,
                -655145967.4845699,
                -648252279.4192085,
                -641224733.3331609,
                -634063329.2264262,
                -626768067.0990052,
                -619338946.9508975,
                -611775968.7821032,
                -604079132.5926222,
                -596248438.3824548,
                -588283886.1516008,
                -580185475.9000597,
                -571953207.6278324,
                -563587081.3349184,
                -555087097.0213174,
                -546453254.6870303,
                -537685554.3320568,
                -528783995.956396,
                -519748579.560049,
                -510579305.14301497,
                -501276172.7052946,
                -491839182.24688727,
                -482268333.76779395,
                -472563627.2680136,
                -462725062.74754655,
                -452752640.2063928,
                -442646359.6445528,
                -432406221.06202614,
                -422032224.4588125,
                -411524369.8349124,
                -400882657.19032645,
                -390107086.5250528,
                -379197657.8390931,
                -368154371.13244605,
                -356977226.40511304,
                -345666223.65709335,
                -334221362.888387,
                -322642644.09899396,
                -310930067.28891426,
                -299083632.45814794,
                -287103339.6066949,
                -274989188.7345552,
                -262741179.84172884,
                -250359312.92821652,
                -237843587.99401715,
                -225194005.03913078,
                -212410564.0635581,
                -199493265.0672991,
                -186442108.0503527,
                -173257093.01272035,
                -531231725.5608198,
                -530935068.75365406,
                -530638884.4552948,
                -530343171.6628055,
                -530047929.37590855,
                -529753156.5969774,
                -529458852.3310286,
                -529165015.58571297,
                -528871645.3713066,
                -528578740.7007033,
                -528286300.5894083,
                -527994324.0555251,
                -527702810.1197512,
                -527411757.80536896,
                -527121166.1382362,
                -526831034.14677966,
                -526541360.86198723,
                -526252145.31739676,
                -525963386.54909116,
                -525675083.5956895,
                -525387235.4983388,
                -525099841.3007059,
                -524812900.04896873,
                -524526410.7918108,
                -524240372.58041143,
                -523954784.4684371,
                -523669645.5120379,
                -523384954.76983285,
                -523100711.302909,
                -522816914.1748085,
                -522533562.4515244,
                -522250655.20149153,
                -521968191.4955778,
                -521686170.40707916,
                -521404591.011709,
                -521123452.3875943,
                -520842753.61526245,
                -520562493.77764124,
                -520282671.9600442,
                -520003287.2501675,
                -519724338.73808074,
                -519445825.51622134,
                -519167746.67938393,
                -518890101.3247149,
                -518612888.5517085,
                -518336107.4621926,
                -518059757.1603254,
                -517783836.752589,
                -517508345.3477809,
                -517233282.05700547,
                -516958645.993669,
                -516684436.273471,
                -516410652.0144,
                -516137292.3367212,
                -515864356.3629741,
                -515591843.2179641,
                -515319752.0287532,
                -515048081.92465824,
                -514776832.03723675,
                -514506001.50028735,
                -514235589.44983745,
                -513965595.02413857,
                -513696017.36366063,
                -513426855.6110814,
                -513158108.9112837,
                -512889776.41134745,
                -512621857.26053953,
                -512354350.61031216,
                -512087255.6142941,
                -511820571.42828256,
                -511554297.21023816,
                -511288432.12027836,
                -511022975.32066864,
                -510757925.9758195,
                -510493283.2522772,
                -510229046.3187175,
                -509965214.3459409,
                -509701786.5068629,
                -509438761.9765105,
                -509176139.93201524,
                -508913919.55260545,
                -508652100.01960135,
                -508390680.5164061,
                -508129660.2285034,
                -507869038.34344786,
                -507608814.0508603,
                -507348986.54242027,
                -507089555.01186085,
                -506830518.65496194,
                -506571876.66954464,
                -506313628.2554637,
                -506055772.6146023,
                -505798308.9508661,
                -505541236.4701766,
                -505284554.380464,
                -505028261.89166397,
                -504772358.2157097,
                -504516842.56652355,
                -504261714.160016,
                -504006972.21407604,
            ],
            expected_sigz=[
                numpy.nan,
                1457137239.0435781,
                1457299027.3539753,
                1457568674.5379708,
                1457946180.5955644,
                1458431545.526756,
                1459024769.3315454,
                1459725852.0099328,
                1460534793.5619192,
                1461451593.9875033,
                1462476253.2866855,
                1463608771.4594655,
                1464849148.5058439,
                1466197384.4258204,
                1467653479.219395,
                1469217432.8865676,
                1470889245.4273384,
                1472668916.8417072,
                1474556447.1296747,
                1476551836.2912397,
                1478655084.326403,
                1480866191.2351644,
                1483185157.017524,
                1485611981.673482,
                1488146665.2030375,
                1490789207.6061914,
                1493539608.8829434,
                1496397869.0332937,
                1499363988.0572422,
                1502437965.9547884,
                1505619802.7259328,
                1508909498.3706753,
                1512307052.8890164,
                1515812466.2809553,
                1519425738.5464919,
                1523146869.6856275,
                1526975859.6983604,
                1530912708.5846918,
                1534957416.3446217,
                1539109982.978149,
                1543370408.4852748,
                1547738692.8659983,
                1552214836.1203206,
                1556798838.2482407,
                1561490699.249759,
                1566290419.1248748,
                1571197997.8735893,
                1576213435.4959023,
                1581336731.9918127,
                1586567887.3613212,
                1591906901.6044283,
                1597353774.7211332,
                1602908506.7114365,
                1608571097.5753376,
                1614341547.3128374,
                1620219855.9239347,
                1626206023.4086304,
                1632300049.7669241,
                1638501934.998816,
                1644811679.1043057,
                1651229282.0833943,
                1657754743.9360802,
                1664388064.6623647,
                1671129244.262247,
                1677978282.7357278,
                1684935180.082806,
                1691999936.3034832,
                1699172551.397758,
                1706453025.365631,
                1713841358.2071023,
                1721337549.9221714,
                1728941600.5108392,
                1736653509.9731042,
                1744473278.308968,
                1752400905.5184302,
                1760436391.6014898,
                1768579736.558148,
                1776830940.3884041,
                1785190003.0922582,
                1793656924.6697106,
                1802231705.1207607,
                1810914344.4454095,
                1819704842.6436563,
                1828603199.7155013,
                1837609415.6609445,
                1846723490.4799852,
                1855945424.1726246,
                1865275216.738862,
                1874712868.178697,
                1884258378.492131,
                1893911747.6791627,
                1903672975.7397926,
                1913542062.6740208,
                1923519008.481846,
                1933603813.1632702,
                1943796476.7182927,
                1954096999.1469133,
                1964505380.4491317,
                1975021620.6249483,
                1985645719.6743631,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
                4115998891.068471,
            ],
            expected_str_r=[
                numpy.nan,
                -0.009810900388423915,
                -0.009809838289219257,
                -0.009808068123873725,
                -0.009805589892390361,
                -0.009802403594768786,
                -0.009798509231008954,
                -0.009793906801110926,
                -0.009788596305074707,
                -0.009782577742900339,
                -0.009775851114587788,
                -0.00976841642013709,
                -0.009760273659548216,
                -0.009751422832821184,
                -0.009741863939955984,
                -0.009731596980952627,
                -0.009720621955811103,
                -0.009708938864531409,
                -0.009696547707113564,
                -0.009683448483557545,
                -0.009669641193863366,
                -0.009655125838031029,
                -0.009639902416060523,
                -0.009623970927951857,
                -0.009607331373705029,
                -0.009589983753320035,
                -0.009571928066796881,
                -0.00955316431413556,
                -0.009533692495336081,
                -0.009513512610398435,
                -0.009492624659322628,
                -0.009471028642108659,
                -0.009448724558756526,
                -0.009425712409266228,
                -0.009401992193637769,
                -0.009377563911871147,
                -0.009352427563966361,
                -0.009326583149923413,
                -0.009300030669742299,
                -0.009272770123423027,
                -0.009244801510965588,
                -0.009216124832369988,
                -0.009186740087636225,
                -0.009156647276764298,
                -0.009125846399754209,
                -0.009094337456605957,
                -0.009062120447319542,
                -0.009029195371894963,
                -0.008995562230332222,
                -0.008961221022631317,
                -0.00892617174879225,
                -0.008890414408815021,
                -0.008853949002699627,
                -0.008816775530446071,
                -0.008778893992054351,
                -0.00874030438752447,
                -0.008701006716856426,
                -0.008661000980050217,
                -0.008620287177105845,
                -0.008578865308023313,
                -0.008536735372802614,
                -0.008493897371443754,
                -0.008450351303946732,
                -0.008406097170311546,
                -0.008361134970538198,
                -0.008315464704626685,
                -0.00826908637257701,
                -0.008221999974389172,
                -0.008174205510063173,
                -0.008125702979599008,
                -0.00807649238299668,
                -0.00802657372025619,
                -0.007975946991377538,
                -0.007924612196360724,
                -0.007872569335205745,
                -0.007819818407912602,
                -0.007766359414481298,
                -0.0077121923549118305,
                -0.007657317229204199,
                -0.007601734037358406,
                -0.007545442779374451,
                -0.007488443455252331,
                -0.0074307360649920485,
                -0.007372320608593602,
                -0.007313197086056994,
                -0.0072533654973822225,
                -0.0071928258425692886,
                -0.007131578121618191,
                -0.007069622334528931,
                -0.007006958481301506,
                -0.00694358656193592,
                -0.006879506576432171,
                -0.006814718524790258,
                -0.006749222407010184,
                -0.006683018223091946,
                -0.006616105973035545,
                -0.006548485656840981,
                -0.006480157274508253,
                -0.006411120826037364,
                -0.00634137631142831,
                -0.005113197055149033,
                -0.005115078293438378,
                -0.005116956535330411,
                -0.005118831787185222,
                -0.005120704055346032,
                -0.005122573346139254,
                -0.005124439665874539,
                -0.005126303020844835,
                -0.005128163417326438,
                -0.00513002086157904,
                -0.00513187535984579,
                -0.005133726918353343,
                -0.005135575543311909,
                -0.00513742124091531,
                -0.005139264017341029,
                -0.005141103878750263,
                -0.005142940831287973,
                -0.005144774881082937,
                -0.005146606034247801,
                -0.0051484342968791275,
                -0.005150259675057449,
                -0.0051520821748473195,
                -0.005153901802297359,
                -0.005155718563440312,
                -0.0051575324642930895,
                -0.005159343510856825,
                -0.005161151709116921,
                -0.005162957065043098,
                -0.005164759584589447,
                -0.005166559273694473,
                -0.005168356138281152,
                -0.005170150184256971,
                -0.0051719414175139845,
                -0.005173729843928855,
                -0.0051755154693629065,
                -0.005177298299662174,
                -0.005179078340657444,
                -0.005180855598164311,
                -0.00518263007798322,
                -0.005184401785899512,
                -0.0051861707276834755,
                -0.005187936909090392,
                -0.005189700335860582,
                -0.005191461013719453,
                -0.005193218948377544,
                -0.0051949741455305735,
                -0.005196726610859487,
                -0.005198476350030496,
                -0.0052002233686951355,
                -0.005201967672490298,
                -0.0052037092670382855,
                -0.005205448157946854,
                -0.005207184350809258,
                -0.005208917851204294,
                -0.0052106486646963475,
                -0.005212376796835437,
                -0.00521410225315726,
                -0.005215825039183232,
                -0.005217545160420537,
                -0.0052192626223621676,
                -0.005220977430486972,
                -0.005222689590259694,
                -0.00522439910713102,
                -0.005226105986537619,
                -0.005227810233902188,
                -0.005229511854633496,
                -0.005231210854126423,
                -0.005232907237762009,
                -0.0052346010109074895,
                -0.005236292178916344,
                -0.005237980747128331,
                -0.0052396667208695425,
                -0.005241350105452432,
                -0.005243030906175865,
                -0.005244709128325158,
                -0.005246384777172121,
                -0.005248057857975098,
                -0.005249728375979008,
                -0.005251396336415387,
                -0.005253061744502428,
                -0.005254724605445026,
                -0.005256384924434811,
                -0.005258042706650194,
                -0.005259697957256405,
                -0.005261350681405538,
                -0.005263000884236582,
                -0.00526464857087547,
                -0.005266293746435115,
                -0.005267936416015449,
                -0.005269576584703461,
                -0.005271214257573243,
                -0.0052728494396860216,
                -0.005274482136090204,
                -0.0052761123518214065,
                -0.005277740091902509,
                -0.005279365361343679,
                -0.0052809881651424175,
                -0.005282608508283598,
                -0.005284226395739498,
                -0.005285841832469847,
            ],
            expected_str_t=[
                numpy.nan,
                -0.009811136410475919,
                -0.00981078237740425,
                -0.009810192322289238,
                -0.00980936624512784,
                -0.009808304145920438,
                -0.009807006024667076,
                -0.009805471881367694,
                -0.00980370171602228,
                -0.009801695528630801,
                -0.009799453319193288,
                -0.009796975087709704,
                -0.00979426083418008,
                -0.009791310558604395,
                -0.009788124260982661,
                -0.00978470194131487,
                -0.009781043599601025,
                -0.009777149235841134,
                -0.009773018850035174,
                -0.009768652442183175,
                -0.009764050012285117,
                -0.009759211560340998,
                -0.00975413708635083,
                -0.009748826590314608,
                -0.009743280072232331,
                -0.009737497532104001,
                -0.009731478969929614,
                -0.009725224385709174,
                -0.009718733779442677,
                -0.009712007151130131,
                -0.00970504450077153,
                -0.00969784582836687,
                -0.009690411133916159,
                -0.009682740417419393,
                -0.009674833678876575,
                -0.009666690918287698,
                -0.00965831213565277,
                -0.009649697330971787,
                -0.00964084650424475,
                -0.009631759655471658,
                -0.009622436784652513,
                -0.009612877891787313,
                -0.009603082976876056,
                -0.009593052039918748,
                -0.009582785080915385,
                -0.00957228209986597,
                -0.009561543096770497,
                -0.009550568071628969,
                -0.00953935702444139,
                -0.009527909955207756,
                -0.009516226863928067,
                -0.009504307750602322,
                -0.009492152615230525,
                -0.009479761457812673,
                -0.009467134278348766,
                -0.009454271076838805,
                -0.00944117185328279,
                -0.00942783660768072,
                -0.009414265340032597,
                -0.009400458050338418,
                -0.009386414738598187,
                -0.0093721354048119,
                -0.009357620048979559,
                -0.009342868671101163,
                -0.009327881271176713,
                -0.00931265784920621,
                -0.009297198405189652,
                -0.009281502939127039,
                -0.009265571451018372,
                -0.00924940394086365,
                -0.009233000408662876,
                -0.009216360854416045,
                -0.009199485278123162,
                -0.009182373679784223,
                -0.009165026059399228,
                -0.009147442416968181,
                -0.00912962275249108,
                -0.009111567065967925,
                -0.009093275357398714,
                -0.009074747626783449,
                -0.009055983874122132,
                -0.009036984099414759,
                -0.00901774830266133,
                -0.008998276483861848,
                -0.008978568643016311,
                -0.008958624780124722,
                -0.008938444895187077,
                -0.008918028988203378,
                -0.008897377059173626,
                -0.008876489108097817,
                -0.008855365134975954,
                -0.008834005139808037,
                -0.008812409122594067,
                -0.008790577083334042,
                -0.008768509022027964,
                -0.008746204938675829,
                -0.00872366483327764,
                -0.008700888705833399,
                -0.008677876556343101,
                -0.00865462838480675,
                -0.008654628384806752,
                -0.008652747146517406,
                -0.008650868904625374,
                -0.008648993652770562,
                -0.008647121384609753,
                -0.00864525209381653,
                -0.008643385774081245,
                -0.00864152241911095,
                -0.008639662022629347,
                -0.008637804578376745,
                -0.008635950080109995,
                -0.008634098521602442,
                -0.008632249896643875,
                -0.008630404199040475,
                -0.008628561422614756,
                -0.008626721561205522,
                -0.008624884608667812,
                -0.008623050558872847,
                -0.008621219405707984,
                -0.008619391143076657,
                -0.008617565764898335,
                -0.008615743265108465,
                -0.008613923637658425,
                -0.008612106876515473,
                -0.008610292975662695,
                -0.008608481929098959,
                -0.008606673730838864,
                -0.008604868374912686,
                -0.008603065855366339,
                -0.008601266166261312,
                -0.008599469301674632,
                -0.008597675255698814,
                -0.0085958840224418,
                -0.00859409559602693,
                -0.008592309970592878,
                -0.00859052714029361,
                -0.00858874709929834,
                -0.008586969841791473,
                -0.008585195361972563,
                -0.008583423654056272,
                -0.008581654712272308,
                -0.008579888530865392,
                -0.008578125104095204,
                -0.00857636442623633,
                -0.00857460649157824,
                -0.008572851294425211,
                -0.008571098829096298,
                -0.008569349089925288,
                -0.008567602071260649,
                -0.008565857767465487,
                -0.0085641161729175,
                -0.00856237728200893,
                -0.008560641089146528,
                -0.00855890758875149,
                -0.008557176775259436,
                -0.008555448643120347,
                -0.008553723186798524,
                -0.008552000400772553,
                -0.008550280279535248,
                -0.008548562817593617,
                -0.008546848009468812,
                -0.00854513584969609,
                -0.008543426332824764,
                -0.008541719453418166,
                -0.008540015206053596,
                -0.00853831358532229,
                -0.008536614585829362,
                -0.008534918202193775,
                -0.008533224429048294,
                -0.008531533261039441,
                -0.008529844692827453,
                -0.008528158719086243,
                -0.008526475334503352,
                -0.00852479453377992,
                -0.008523116311630627,
                -0.008521440662783663,
                -0.008519767581980687,
                -0.008518097063976778,
                -0.008516429103540397,
                -0.008514763695453356,
                -0.008513100834510758,
                -0.008511440515520973,
                -0.00850978273330559,
                -0.00850812748269938,
                -0.008506474758550247,
                -0.008504824555719203,
                -0.008503176869080314,
                -0.00850153169352067,
                -0.008499889023940336,
                -0.008498248855252324,
                -0.00849661118238254,
                -0.008494976000269763,
                -0.00849334330386558,
                -0.008491713088134378,
                -0.008490085348053275,
                -0.008488460078612104,
                -0.008486837274813367,
                -0.008485216931672188,
                -0.008483599044216287,
                -0.008481983607485937,
            ],
            expected_str_z=[
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
                0.020815614549915578,
            ],
            expected_r_deflect=[
                numpy.nan,
                -1.457682159507361e-05,
                -2.9152591186238894e-05,
                -4.3726256769607287e-05,
                -5.829676634127218e-05,
                -7.286306789733426e-05,
                -8.742410943389353e-05,
                -0.00010197883894704914,
                -0.00011652620443290045,
                -0.00013106515388754654,
                -0.00014559463530708708,
                -0.00016011359668762097,
                -0.00017462098602524814,
                -0.00018911575131606742,
                -0.00020359684055617854,
                -0.00021806320174168066,
                -0.0002325137828686732,
                -0.00024694753193325576,
                -0.000261363396931527,
                -0.0002757603258595873,
                -0.00029013726671353525,
                -0.00030449316748947024,
                -0.0003188269761834921,
                -0.0003331376407916998,
                -0.0003474241093101928,
                -0.00036168532973507056,
                -0.0003759202500624323,
                -0.0003901278182883774,
                -0.0004043069824090053,
                -0.00041845669042041544,
                -0.00043257589031870703,
                -0.00044666353009997935,
                -0.00046071855776033194,
                -0.0004747399212958642,
                -0.0004887265687026754,
                -0.0005026774479768648,
                -0.000516591507114532,
                -0.000530467694111776,
                -0.0005443049569646965,
                -0.0005581022436693927,
                -0.0005718585022219641,
                -0.00058557268061851,
                -0.0005992437268551295,
                -0.0006128705889279224,
                -0.000626452214832988,
                -0.0006399875525664253,
                -0.0006534755501243339,
                -0.000666915155502813,
                -0.0006803053166979625,
                -0.0006936449817058811,
                -0.0007069330985226685,
                -0.0007201686151444238,
                -0.0007333504795672467,
                -0.0007464776397872365,
                -0.0007595490438004923,
                -0.0007725636396031136,
                -0.0007855203751912001,
                -0.0007984181985608506,
                -0.0008112560577081647,
                -0.0008240329006292419,
                -0.0008367476753201815,
                -0.0008493993297770826,
                -0.0008619868119960451,
                -0.0008745090699731677,
                -0.0008869650517045502,
                -0.000899353705186292,
                -0.0009116739784144921,
                -0.0009239248193852503,
                -0.0009361051760946656,
                -0.0009482139965388374,
                -0.0009602502287138653,
                -0.0009722128206158484,
                -0.0009841007202408864,
                -0.0009959128755850785,
                -0.0010076482346445235,
                -0.0010193057454153216,
                -0.0010308843558935717,
                -0.0010423830140753737,
                -0.001053800667956826,
                -0.0010651362655340292,
                -0.0010763887548030813,
                -0.001087557083760083,
                -0.0010986402004011322,
                -0.0011096370527223297,
                -0.001120546588719774,
                -0.0011313677563895644,
                -0.0011420995037278012,
                -0.0011527407787305826,
                -0.0011632905293940089,
                -0.0011737477037141785,
                -0.0011841112496871916,
                -0.0011943801153091473,
                -0.0012045532485761446,
                -0.0012146295974842835,
                -0.001224608110029663,
                -0.0012344877342083827,
                -0.0012442674180165416,
                -0.0012539461094502394,
                -0.0012635227565055748,
                -0.0012729963071786477,
                -0.001272996307178648,
                -0.0012733962172669375,
                -0.0012737962743451642,
                -0.001274196478179268,
                -0.001274596828535685,
                -0.001274997325181347,
                -0.0012753979678836796,
                -0.0012757987564106011,
                -0.0012761996905305216,
                -0.001276600770012341,
                -0.001277001994625448,
                -0.0012774033641397191,
                -0.0012778048783255171,
                -0.0012782065369536895,
                -0.0012786083397955676,
                -0.001279010286622965,
                -0.0012794123772081772,
                -0.0012798146113239782,
                -0.0012802169887436224,
                -0.0012806195092408394,
                -0.0012810221725898369,
                -0.0012814249785652963,
                -0.001281827926942373,
                -0.0012822310174966948,
                -0.0012826342500043605,
                -0.0012830376242419388,
                -0.001283441139986467,
                -0.00128384479701545,
                -0.0012842485951068592,
                -0.0012846525340391296,
                -0.0012850566135911615,
                -0.0012854608335423172,
                -0.0012858651936724203,
                -0.001286269693761755,
                -0.001286674333591063,
                -0.0012870791129415454,
                -0.0012874840315948595,
                -0.0012878890893331173,
                -0.0012882942859388856,
                -0.0012886996211951833,
                -0.0012891050948854823,
                -0.0012895107067937043,
                -0.0012899164567042208,
                -0.0012903223444018517,
                -0.0012907283696718629,
                -0.0012911345322999684,
                -0.0012915408320723254,
                -0.0012919472687755351,
                -0.0012923538421966412,
                -0.001292760552123129,
                -0.0012931673983429235,
                -0.0012935743806443897,
                -0.0012939814988163294,
                -0.0012943887526479824,
                -0.001294796141929023,
                -0.0012952036664495604,
                -0.0012956113260001385,
                -0.0012960191203717317,
                -0.0012964270493557467,
                -0.0012968351127440195,
                -0.0012972433103288158,
                -0.001297651641902829,
                -0.0012980601072591796,
                -0.0012984687061914123,
                -0.0012988774384934984,
                -0.0012992863039598315,
                -0.0012996953023852273,
                -0.0013001044335649245,
                -0.00130051369729458,
                -0.0013009230933702712,
                -0.001301332621588493,
                -0.0013017422817461575,
                -0.0013021520736405932,
                -0.0013025619970695422,
                -0.0013029720518311619,
                -0.0013033822377240215,
                -0.0013037925545471023,
                -0.0013042030020997966,
                -0.0013046135801819052,
                -0.0013050242885936387,
                -0.0013054351271356147,
                -0.0013058460956088568,
                -0.0013062571938147953,
                -0.0013066684215552635,
                -0.0013070797786324992,
                -0.0013074912648491418,
                -0.0013079028800082324,
                -0.0013083146239132126,
                -0.001308726496367923,
                -0.0013091384971766027,
                -0.0013095506261438876,
                -0.0013099628830748108,
                -0.0013103752677747999,
                -0.0013107877800496767,
                -0.0013112004197056567,
                -0.0013116131865493477,
                -0.0013120260803877485,
                -0.0013124391010282486,
                -0.001312852248278626,
                -0.0013132655219470477,
            ],
            expected_rradius=[
                0.0,
                0.0014857424242424244,
                0.0029714848484848487,
                0.004457227272727273,
                0.0059429696969696974,
                0.007428712121212122,
                0.008914454545454547,
                0.01040019696969697,
                0.011885939393939395,
                0.013371681818181819,
                0.014857424242424243,
                0.01634316666666667,
                0.017828909090909093,
                0.019314651515151517,
                0.02080039393939394,
                0.022286136363636366,
                0.02377187878787879,
                0.025257621212121214,
                0.026743363636363638,
                0.028229106060606062,
                0.029714848484848486,
                0.03120059090909091,
                0.03268633333333334,
                0.03417207575757576,
                0.035657818181818186,
                0.03714356060606061,
                0.038629303030303035,
                0.04011504545454546,
                0.04160078787878788,
                0.04308653030303031,
                0.04457227272727273,
                0.046058015151515155,
                0.04754375757575758,
                0.049029500000000004,
                0.05051524242424243,
                0.05200098484848485,
                0.053486727272727276,
                0.0549724696969697,
                0.056458212121212124,
                0.05794395454545455,
                0.05942969696969697,
                0.0609154393939394,
                0.06240118181818182,
                0.06388692424242425,
                0.06537266666666668,
                0.0668584090909091,
                0.06834415151515152,
                0.06982989393939394,
                0.07131563636363637,
                0.07280137878787879,
                0.07428712121212122,
                0.07577286363636364,
                0.07725860606060607,
                0.07874434848484849,
                0.08023009090909092,
                0.08171583333333333,
                0.08320157575757577,
                0.08468731818181818,
                0.08617306060606061,
                0.08765880303030303,
                0.08914454545454546,
                0.09063028787878788,
                0.09211603030303031,
                0.09360177272727273,
                0.09508751515151516,
                0.09657325757575759,
                0.09805900000000001,
                0.09954474242424244,
                0.10103048484848486,
                0.10251622727272729,
                0.1040019696969697,
                0.10548771212121213,
                0.10697345454545455,
                0.10845919696969698,
                0.1099449393939394,
                0.11143068181818183,
                0.11291642424242425,
                0.11440216666666668,
                0.1158879090909091,
                0.11737365151515153,
                0.11885939393939395,
                0.12034513636363638,
                0.1218308787878788,
                0.12331662121212122,
                0.12480236363636364,
                0.12628810606060606,
                0.1277738484848485,
                0.12925959090909092,
                0.13074533333333335,
                0.13223107575757576,
                0.1337168181818182,
                0.13520256060606062,
                0.13668830303030305,
                0.13817404545454545,
                0.13965978787878788,
                0.14114553030303031,
                0.14263127272727275,
                0.14411701515151518,
                0.14560275757575758,
                0.1470885,
                0.1470885,
                0.14716669696969698,
                0.14724489393939394,
                0.14732309090909093,
                0.1474012878787879,
                0.14747948484848486,
                0.14755768181818182,
                0.1476358787878788,
                0.14771407575757578,
                0.14779227272727274,
                0.1478704696969697,
                0.14794866666666667,
                0.14802686363636364,
                0.14810506060606063,
                0.1481832575757576,
                0.14826145454545456,
                0.14833965151515152,
                0.14841784848484849,
                0.14849604545454548,
                0.14857424242424244,
                0.1486524393939394,
                0.14873063636363637,
                0.14880883333333333,
                0.14888703030303033,
                0.1489652272727273,
                0.14904342424242426,
                0.14912162121212122,
                0.14919981818181818,
                0.14927801515151518,
                0.14935621212121214,
                0.1494344090909091,
                0.14951260606060607,
                0.14959080303030303,
                0.14966900000000002,
                0.149747196969697,
                0.14982539393939395,
                0.14990359090909092,
                0.14998178787878788,
                0.15005998484848487,
                0.15013818181818184,
                0.1502163787878788,
                0.15029457575757577,
                0.15037277272727273,
                0.15045096969696972,
                0.15052916666666669,
                0.15060736363636365,
                0.15068556060606061,
                0.1507637575757576,
                0.15084195454545457,
                0.15092015151515153,
                0.1509983484848485,
                0.15107654545454546,
                0.15115474242424243,
                0.15123293939393942,
                0.15131113636363638,
                0.15138933333333335,
                0.1514675303030303,
                0.1515457272727273,
                0.15162392424242427,
                0.15170212121212123,
                0.1517803181818182,
                0.15185851515151516,
                0.15193671212121215,
                0.15201490909090912,
                0.15209310606060608,
                0.15217130303030305,
                0.1522495,
                0.152327696969697,
                0.15240589393939397,
                0.15248409090909093,
                0.1525622878787879,
                0.15264048484848486,
                0.15271868181818185,
                0.15279687878787881,
                0.15287507575757578,
                0.15295327272727274,
                0.1530314696969697,
                0.1531096666666667,
                0.15318786363636366,
                0.15326606060606063,
                0.1533442575757576,
                0.15342245454545456,
                0.15350065151515155,
                0.1535788484848485,
                0.15365704545454548,
                0.15373524242424244,
                0.1538134393939394,
                0.1538916363636364,
                0.15396983333333336,
                0.15404803030303033,
                0.1541262272727273,
                0.15420442424242425,
                0.15428262121212125,
                0.1543608181818182,
                0.15443901515151517,
                0.15451721212121214,
                0.1545954090909091,
                0.1546736060606061,
                0.15475180303030306,
                0.15483000000000002,
            ],
        ),
        ExtendedPlaneStrainParam(
            n_radial_array=100,
            nlayers=5,
            i_tf_bucking=3,
            nu_t=numpy.array(
                numpy.array(
                    (
                        0.30000000502169133,
                        0.34000000000000002,
                        0.29999999999999999,
                        0.30901178507421895,
                        0.29999999999999999,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            nu_zt=numpy.array(
                numpy.array(
                    (
                        0.31163570564277626,
                        0.34000000000000002,
                        0.29999999999999999,
                        0.31377709779186291,
                        0.29999999999999999,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            ey_t=numpy.array(
                numpy.array(
                    (
                        118643750000,
                        2500000000,
                        205000000000,
                        43163597776.087654,
                        205000000000,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            ey_z=numpy.array(
                numpy.array(
                    (
                        48005309351.198608,
                        2500000000,
                        205000000000,
                        124208626934.75433,
                        390554854819.81116,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            d_curr=numpy.array(
                numpy.array((0, 0, 0, 18343613.563061949, 0), order="F"), order="F"
            ).transpose(),
            rad=numpy.array(
                numpy.array(
                    (
                        2.3322000000000003,
                        2.8846200000000004,
                        2.9346200000000002,
                        3.4817726429672304,
                        4.0290604740948242,
                        4.0890604740948238,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            v_force=4051971733.3410816,
            expected_sigr=[
                -8.205145940973878e-07,
                -415231.5173357017,
                -827500.3245084467,
                -1236834.540195528,
                -1643261.9502749545,
                -2046810.012542044,
                -2447505.8613465936,
                -2845376.3121548747,
                -3240447.8660354516,
                -3632746.7140727337,
                -4022298.741707528,
                -4409129.533007751,
                -4793264.374868267,
                -5174728.261144401,
                -5553545.896717039,
                -5929741.701492815,
                -6303339.814339905,
                -6674364.096960006,
                -7042838.1376985,
                -7408785.255294099,
                -7772228.5025681555,
                -8133190.670054784,
                -8491694.289574316,
                -8847761.637749428,
                -9201414.739465648,
                -9552675.371277362,
                -9901565.064760374,
                -10248105.109811509,
                -10592316.557896374,
                -10934220.225246798,
                -11273836.696007827,
                -11611186.325336073,
                -11946289.242449744,
                -12279165.353631394,
                -12609834.345184395,
                -12938315.686343344,
                -13264628.632140255,
                -13588792.226226069,
                -13910825.303649498,
                -14230746.493592657,
                -14548574.2220656,
                -14864326.714559477,
                -15178021.99865938,
                -15489677.906617459,
                -15799312.077887261,
                -16106941.961619942,
                -16412584.819122143,
                -16716257.726278093,
                -17017977.57593395,
                -17317761.080247745,
                -17615624.77300322,
                -17911585.011889998,
                -18205657.980748963,
                -18497859.691784922,
                -18788205.987746373,
                -19076712.544073083,
                -19363394.87101139,
                -19648268.31569944,
                -19931348.06422071,
                -20212649.143628374,
                -20492186.42393912,
                -20769974.620098636,
                -21046028.293917917,
                -21320361.85598206,
                -21592989.56753048,
                -21863925.542310633,
                -22133183.748404544,
                -22400778.01002948,
                -22666722.009311974,
                -22931029.288036954,
                -23193713.24937192,
                -23454787.15956555,
                -23714264.149623252,
                -23972157.216958147,
                -24228479.227018487,
                -24483242.914892554,
                -24736460.8868904,
                -24988145.622103136,
                -25238309.47394085,
                -25486964.67164802,
                -25734123.32179869,
                -25979797.409769565,
                -26223998.80119323,
                -26466739.243391115,
                -26708030.366786003,
                -26947883.68629526,
                -27186310.602704655,
                -27423322.404023346,
                -27658930.26681952,
                -27893145.2575384,
                -28125978.333801653,
                -28357440.345688872,
                -28587542.037001874,
                -28816294.046511214,
                -29043706.90918619,
                -29269791.057407647,
                -29494556.82216431,
                -29718014.43423333,
                -29940174.02534434,
                -30161045.629327673,
                -30161045.629327886,
                -30158858.1982904,
                -30156671.915702056,
                -30154486.78075904,
                -30152302.792658247,
                -30150119.950597282,
                -30147938.25377445,
                -30145757.701388743,
                -30143578.292639863,
                -30141400.026728213,
                -30139222.902854886,
                -30137046.920221683,
                -30134872.078031097,
                -30132698.375486318,
                -30130525.81179124,
                -30128354.38615044,
                -30126184.097769186,
                -30124014.945853464,
                -30121846.92960993,
                -30119680.048245925,
                -30117514.300969515,
                -30115349.686989423,
                -30113186.205515087,
                -30111023.855756607,
                -30108862.636924807,
                -30106702.54823116,
                -30104543.58888785,
                -30102385.758107755,
                -30100229.055104423,
                -30098073.47909208,
                -30095919.02928566,
                -30093765.70490075,
                -30091613.505153637,
                -30089462.429261312,
                -30087312.47644141,
                -30085163.645912264,
                -30083015.936892882,
                -30080869.348602958,
                -30078723.880262855,
                -30076579.53109363,
                -30074436.300316997,
                -30072294.187155362,
                -30070153.1908318,
                -30068013.31057005,
                -30065874.54559455,
                -30063736.8951304,
                -30061600.35840336,
                -30059464.93463988,
                -30057330.623067077,
                -30055197.42291273,
                -30053065.333405305,
                -30050934.353773918,
                -30048804.483248368,
                -30046675.721059114,
                -30044548.066437297,
                -30042421.518614687,
                -30040296.07682377,
                -30038171.740297664,
                -30036048.508270167,
                -30033926.379975718,
                -30031805.354649447,
                -30029685.431527138,
                -30027566.60984523,
                -30025448.888840828,
                -30023332.267751694,
                -30021216.74581625,
                -30019102.322273582,
                -30016988.99636343,
                -30014876.767326206,
                -30012765.634402946,
                -30010655.596835367,
                -30008546.653865848,
                -30006438.804737404,
                -30004332.048693698,
                -30002226.384979088,
                -30000121.812838543,
                -29998018.331517696,
                -29995915.94026284,
                -29993814.638320915,
                -29991714.424939506,
                -29989615.29936684,
                -29987517.26085183,
                -29985420.30864398,
                -29983324.441993505,
                -29981229.660151202,
                -29979135.962368567,
                -29977043.347897716,
                -29974951.81599142,
                -29972861.365903076,
                -29970771.996886745,
                -29968683.708197117,
                -29966596.499089547,
                -29964510.36882,
                -29962425.316645097,
                -29960341.34182211,
                -29958258.443608917,
                -29956176.62126409,
                -29954095.874046776,
                -29952016.20121681,
                -29949937.602034632,
                -29949937.60203456,
                -30300932.11275363,
                -30649952.82628681,
                -30997014.51431862,
                -31342131.81060446,
                -31685319.212513402,
                -32026591.082550496,
                -32365961.649859756,
                -32703445.011707082,
                -33039055.13494445,
                -33372805.857454974,
                -33704710.88957932,
                -34034783.81552361,
                -34363038.09474966,
                -34689487.063346975,
                -35014143.93538717,
                -35337021.80426137,
                -35658133.64400036,
                -35977492.31057784,
                -36295110.543197155,
                -36611000.9655619,
                -36925176.087130085,
                -37237648.30435279,
                -37548429.90189682,
                -37857533.05385223,
                -38164969.82492459,
                -38470752.17161217,
                -38774891.94336859,
                -39077400.88375059,
                -39378290.63155183,
                -39677572.721922435,
                -39975258.58747437,
                -40271359.55937333,
                -40565886.8684172,
                -40858851.646100774,
                -41150264.92566782,
                -41440137.64314948,
                -41728480.63839061,
                -42015304.656063296,
                -42300620.34666741,
                -42584438.26751965,
                -42866768.883730076,
                -43147622.569166854,
                -43427009.607409135,
                -43704940.19268834,
                -43981424.4308183,
                -44256472.34011384,
                -44530093.85229813,
                -44802298.81339952,
                -45073096.984637186,
                -45342498.04329618,
                -45610511.583592065,
                -45877147.117525175,
                -46142414.075724564,
                -46406321.80828193,
                -46668879.58557528,
                -46930096.599083625,
                -47189981.96219079,
                -47448544.71098072,
                -47705793.80502272,
                -47961738.12814784,
                -48216386.48921557,
                -48469747.62287199,
                -48721830.19029882,
                -48972642.77995357,
                -49222193.90830097,
                -49470492.02053617,
                -49717545.49129898,
                -49963362.62538021,
                -50207951.658419676,
                -50451320.757595725,
                -50693478.022307225,
                -50934431.484847434,
                -51174189.111070134,
                -51412758.80104776,
                -51650148.38972252,
                -51886365.64754948,
                -52121418.28113245,
                -52355313.933852635,
                -52588060.186489895,
                -52819664.55783726,
                -53050134.50530772,
                -53279477.42553502,
                -53507700.654966846,
                -53734811.470451556,
                -53960817.08981835,
                -54185724.67245068,
                -54409541.319853276,
                -54632274.07621272,
                -54853929.92895165,
                -55074515.80927691,
                -55294038.592721224,
                -55512505.09967921,
                -55729922.0959369,
                -55946296.29319571,
                -56161634.34959044,
                -56375942.87020137,
                -56589228.40756096,
                -56801497.462154366,
                -57012756.48291527,
                -57012756.4829153,
                -56998898.77366743,
                -56972187.46576252,
                -56932655.747310705,
                -56880336.58824644,
                -56815262.742183134,
                -56737466.74824842,
                -56646980.93290269,
                -56543837.41173834,
                -56428068.091256954,
                -56299704.67063484,
                -56158778.64346349,
                -56005321.29947959,
                -55839363.72626889,
                -55660936.810961984,
                -55470071.24190524,
                -55266797.51031879,
                -55051145.91193752,
                -54823146.54863426,
                -54582829.33002596,
                -54330223.97506812,
                -54065360.01362707,
                -53788266.788039714,
                -53498973.4546594,
                -53197508.985382356,
                -52883902.16916223,
                -52558181.613506556,
                -52220375.74596186,
                -51870512.81558069,
                -51508620.89437677,
                -51134727.87876348,
                -50748861.49097967,
                -50351049.28050165,
                -49941318.625439376,
                -49519696.733920544,
                -49086210.645462126,
                -48640887.23232431,
                -48183753.20085678,
                -47714835.09282659,
                -47234159.28673844,
                -46741751.999135256,
                -46237639.285894305,
                -45721847.04350362,
                -45194401.01033023,
                -44655326.76787479,
                -44104649.742013805,
                -43542395.20423049,
                -42968588.27283408,
                -42383253.9141667,
                -41786416.94379792,
                -41178102.02771103,
                -40558333.683473185,
                -39927136.28140024,
                -39284534.045703985,
                -38630551.05563459,
                -37965211.246607006,
                -37288538.41132125,
                -36600556.20086938,
                -35901288.125830434,
                -35190757.55736128,
                -34468987.728269,
                -33736001.734080136,
                -32991822.53409599,
                -32236472.952438682,
                -31469975.67908989,
                -30692353.270914197,
                -29903628.15268061,
                -29103822.61806818,
                -28292958.83066472,
                -27471058.82495762,
                -26638144.507311735,
                -25794237.656943187,
                -24939359.926879253,
                -24073532.844912816,
                -23196777.8145469,
                -22309116.115929943,
                -21410568.90678382,
                -20501157.223321468,
                -19580901.9811581,
                -18649823.97621408,
                -17707943.885606322,
                -16755282.268536685,
                -15791859.56716701,
                -14817696.107491374,
                -13832812.100195508,
                -12837227.641513025,
                -11830962.714069963,
                -10814037.18772364,
                -9786470.820395064,
                -8748283.258890826,
                -7699494.03972132,
                -6640122.589907372,
                -5570188.2277836045,
                -4489710.163792794,
                -3398707.5012716674,
                -2297199.237235021,
                -1185204.2631434032,
                -62741.36567673539,
                1070170.772511864,
                2213513.5720397322,
                2213513.5720397546,
                2190659.3327056696,
                2167815.402857395,
                2144981.776294998,
                2122158.4468234032,
                2099345.4082520274,
                2076542.6543950075,
                2053750.179071122,
                2030967.9761037766,
                2008196.039321003,
                1985434.3625555413,
                1962682.939644665,
                1939941.764430341,
                1917210.8307590964,
                1894490.1324821794,
                1871779.663455358,
                1849079.417539027,
                1826389.3885982349,
                1803709.5705025634,
                1781039.9571262214,
                1758380.5423480715,
                1735731.320051415,
                1713092.2841242745,
                1690463.4284593,
                1667844.746953472,
                1645236.2335086272,
                1622637.88203096,
                1600049.686431346,
                1577471.640625207,
                1554903.7385324438,
                1532345.9740776382,
                1509798.341189797,
                1487260.833802528,
                1464733.445854011,
                1442216.1712868656,
                1419709.0040483377,
                1397211.9380901663,
                1374724.967368583,
                1352248.0858444187,
                1329781.2874829038,
                1307324.5662538684,
                1284877.9161316748,
                1262441.3310950447,
                1240014.805127393,
                1217598.3322164528,
                1195191.9063545435,
                1172795.5215385184,
                1150409.1717695605,
                1128032.8510534272,
                1105666.553400435,
                1083310.2728251917,
                1060964.0033469184,
                1038627.7389891675,
                1016301.4737800913,
                993985.2017522273,
                971678.9169425247,
                949382.6133924657,
                927096.2851479175,
                904819.9262591992,
                882553.5307810823,
                860297.0927726965,
                838050.6062976907,
                815814.0654241928,
                793587.4642245143,
                771370.7967756338,
                749164.057158767,
                726967.2394596086,
                704780.3377683188,
                682603.3461793356,
                660436.2587915623,
                638279.0697083006,
                616131.7730371836,
                593994.3628903095,
                571866.8333840411,
                549749.1786392601,
                527641.3927811258,
                505543.46993916895,
                483455.404247332,
                461377.1898438619,
                439308.82087143057,
                417250.2914769876,
                395201.5958118269,
                373162.7280317351,
                351133.682296655,
                329114.4527709279,
                307105.0336233069,
                285105.4190267686,
                263115.60315863416,
                241135.58020066356,
                219165.34433870582,
                197204.8897631966,
                175254.2106686337,
                153313.3012539806,
                131382.15572247806,
                109460.76828160389,
                87549.13314318044,
                65647.24452332086,
                43755.09664238875,
                21872.68372509228,
                3.3626858094974444e-07,
            ],
            expected_sigt=[
                -174172010.1748907,
                -173756778.65755582,
                -173344509.85038307,
                -172935175.63469604,
                -172528748.2246165,
                -172125200.1623495,
                -171724504.31354484,
                -171326633.86273664,
                -170931562.30885598,
                -170539263.46081883,
                -170149711.4331839,
                -169762880.6418838,
                -169378745.80002326,
                -168997281.91374716,
                -168618464.27817452,
                -168242268.47339863,
                -167868670.36055154,
                -167497646.07793155,
                -167129172.03719303,
                -166763224.91959733,
                -166399781.6723233,
                -166038819.50483674,
                -165680315.88531715,
                -165324248.53714204,
                -164970595.4354259,
                -164619334.80361417,
                -164270445.11013114,
                -163923905.06508,
                -163579693.61699513,
                -163237789.94964468,
                -162898173.4788837,
                -162560823.84955546,
                -162225720.93244174,
                -161892844.82126004,
                -161562175.82970718,
                -161233694.48854813,
                -160907381.54275116,
                -160583217.94866538,
                -160261184.871242,
                -159941263.6812989,
                -159623435.9528259,
                -159307683.460332,
                -158993988.17623213,
                -158682332.26827404,
                -158372698.0970042,
                -158065068.21327165,
                -157759425.35576928,
                -157455752.44861346,
                -157154032.59895757,
                -156854249.09464377,
                -156556385.40188825,
                -156260425.1630015,
                -155966352.19414258,
                -155674150.48310658,
                -155383804.18714508,
                -155095297.63081843,
                -154808615.30388018,
                -154523741.85919216,
                -154240662.11067075,
                -153959361.03126308,
                -153679823.75095242,
                -153402035.55479294,
                -153125981.88097355,
                -152851648.3189094,
                -152579020.607361,
                -152308084.63258088,
                -152038826.42648688,
                -151771232.164862,
                -151505288.16557953,
                -151240980.8868545,
                -150978296.92551962,
                -150717223.015326,
                -150457746.0252682,
                -150199852.9579334,
                -149943530.94787297,
                -149688767.25999886,
                -149435549.2880011,
                -149183864.5527884,
                -148933700.70095074,
                -148685045.50324348,
                -148437886.85309276,
                -148192212.76512197,
                -147948011.37369823,
                -147705270.93150038,
                -147463979.80810553,
                -147224126.48859626,
                -146985699.5721869,
                -146748687.77086818,
                -146513079.90807205,
                -146278864.91735312,
                -146046031.84108987,
                -145814569.82920262,
                -145584468.13788968,
                -145355716.12838024,
                -145128303.26570526,
                -144902219.11748385,
                -144677453.35272714,
                -144453995.74065813,
                -144231836.1495472,
                -144010964.54556382,
                -17664147.894719966,
                -17666335.325757444,
                -17668521.608345795,
                -17670706.743288808,
                -17672890.731389597,
                -17675073.573450554,
                -17677255.27027339,
                -17679435.822659105,
                -17681615.23140799,
                -17683793.497319635,
                -17685970.621192966,
                -17688146.60382617,
                -17690321.44601675,
                -17692495.148561526,
                -17694667.712256603,
                -17696839.137897402,
                -17699009.42627865,
                -17701178.578194376,
                -17703346.59443792,
                -17705513.47580192,
                -17707679.22307833,
                -17709843.83705842,
                -17712007.31853276,
                -17714169.668291237,
                -17716330.887123045,
                -17718490.975816686,
                -17720649.935159985,
                -17722807.76594008,
                -17724964.46894342,
                -17727120.044955764,
                -17729274.49476219,
                -17731427.8191471,
                -17733580.018894203,
                -17735731.094786532,
                -17737881.047606435,
                -17740029.878135584,
                -17742177.587154962,
                -17744324.175444886,
                -17746469.643784985,
                -17748613.99295421,
                -17750757.223730847,
                -17752899.336892482,
                -17755040.33321605,
                -17757180.21347779,
                -17759318.978453286,
                -17761456.628917444,
                -17763593.16564449,
                -17765728.58940796,
                -17767862.900980763,
                -17769996.10113511,
                -17772128.190642536,
                -17774259.170273922,
                -17776389.040799472,
                -17778517.80298872,
                -17780645.45761055,
                -17782772.005433153,
                -17784897.44722407,
                -17787021.783750176,
                -17789145.015777677,
                -17791267.144072123,
                -17793388.169398393,
                -17795508.092520703,
                -17797626.914202612,
                -17799744.635207012,
                -17801861.256296147,
                -17803976.77823159,
                -17806091.201774254,
                -17808204.5276844,
                -17810316.756721642,
                -17812427.889644895,
                -17814537.927212473,
                -17816646.870181993,
                -17818754.719310444,
                -17820861.47535414,
                -17822967.139068753,
                -17825071.711209297,
                -17827175.192530144,
                -17829277.583784997,
                -17831378.885726925,
                -17833479.099108342,
                -17835578.224680994,
                -17837676.26319601,
                -17839773.215403855,
                -17841869.082054336,
                -17843963.863896634,
                -17846057.561679266,
                -17848150.17615012,
                -17850241.708056424,
                -17852332.15814477,
                -17854421.527161095,
                -17856509.81585072,
                -17858597.024958294,
                -17860683.15522784,
                -17862768.207402743,
                -17864852.182225738,
                -17866935.08043892,
                -17869016.90278376,
                -17871097.650001064,
                -17873177.322831035,
                -17875255.92201321,
                -216847832.74234107,
                -216496838.23162204,
                -216147817.5180889,
                -215800755.83005714,
                -215455638.53377122,
                -215112451.13186234,
                -214771179.2618252,
                -214431808.69451594,
                -214094325.33266863,
                -213758715.20943123,
                -213424964.48692068,
                -213093059.4547964,
                -212762986.52885205,
                -212434732.249626,
                -212108283.28102878,
                -211783626.40898857,
                -211460748.54011437,
                -211139636.70037535,
                -210820278.0337979,
                -210502659.80117857,
                -210186769.37881377,
                -209872594.25724557,
                -209560122.0400229,
                -209249340.44247892,
                -208940237.29052353,
                -208632800.51945105,
                -208327018.17276353,
                -208022878.4010071,
                -207720369.46062514,
                -207419479.71282384,
                -207120197.6224532,
                -206822511.7569014,
                -206526410.78500238,
                -206231883.47595856,
                -205938918.6982749,
                -205647505.4187079,
                -205357632.7012263,
                -205069289.70598507,
                -204782465.68831244,
                -204497149.99770835,
                -204213332.07685608,
                -203931001.4606456,
                -203650147.7752089,
                -203370760.73696658,
                -203092830.15168735,
                -202816345.91355738,
                -202541298.00426194,
                -202267676.4920776,
                -201995471.53097615,
                -201724673.35973856,
                -201455272.30107954,
                -201187258.7607837,
                -200920623.2268505,
                -200655356.26865107,
                -200391448.53609383,
                -200128890.7588004,
                -199867673.7452921,
                -199607788.38218486,
                -199349225.633395,
                -199091976.53935292,
                -198836032.2162279,
                -198581383.8551602,
                -198328022.7215037,
                -198075940.15407687,
                -197825127.56442213,
                -197575576.43607476,
                -197327278.32383955,
                -197080224.85307673,
                -196834407.71899545,
                -196589818.68595603,
                -196346449.58677998,
                -196104292.32206848,
                -195863338.85952827,
                -195623581.23330563,
                -195385011.543328,
                -195147621.9546532,
                -194911404.69682622,
                -194676352.0632433,
                -194442456.41052315,
                -194209710.1578858,
                -193978105.78653848,
                -193747635.83906797,
                -193518292.91884065,
                -193290069.68940887,
                -193062958.87392414,
                -192836953.2545574,
                -192612045.67192504,
                -192388229.02452245,
                -192165496.268163,
                -191943840.41542408,
                -191723254.5350988,
                -191503731.75165448,
                -191285265.24469644,
                -191067848.24843884,
                -190851474.05117998,
                -190636135.9947853,
                -190421827.47417432,
                -190208541.93681476,
                -189996272.88222134,
                -189785013.86146045,
                -52339586.701563686,
                -52344691.67092526,
                -52345163.23240329,
                -52340995.80914617,
                -52332183.911686264,
                -52318722.13691016,
                -52300605.1670393,
                -52277827.76862575,
                -52250384.79156041,
                -52218271.16808853,
                -52181481.911844626,
                -52140012.11689179,
                -52093856.956779905,
                -52043011.68360702,
                -51987471.62710146,
                -51927232.193706565,
                -51862288.865682095,
                -51792637.20021451,
                -51718272.82853808,
                -51639191.455065206,
                -51555388.85653211,
                -51466860.88114854,
                -51373603.44776061,
                -51275612.545025885,
                -51172884.230594546,
                -51065414.63030283,
                -50953199.937374085,
                -50836236.41163265,
                -50714520.378723435,
                -50588048.2293438,
                -50456816.41848262,
                -50320821.46466913,
                -50180059.949231826,
                -50034528.51556382,
                -49884223.868398644,
                -49729142.77309577,
                -49569282.05493002,
                -49404638.59839572,
                -49235209.34651348,
                -49060991.30014996,
                -48881981.517339766,
                -48698177.11262308,
                -48509575.25638367,
                -48316173.17419965,
                -48117968.1461997,
                -47914957.50642701,
                -47707138.64221125,
                -47494508.99354834,
                -47277066.05248651,
                -47054807.362518564,
                -46827730.51798531,
                -46595833.16348018,
                -46359112.99326728,
                -46117567.750699386,
                -45871195.22764921,
                -45619993.26394161,
                -45363959.74679763,
                -45103092.61028124,
                -44837389.83475229,
                -44566849.44633101,
                -44291469.51636085,
                -44011248.16088445,
                -43726183.54012203,
                -43436273.85795625,
                -43141517.361425206,
                -42841912.34021622,
                -42537457.12617256,
                -42228150.0928,
                -41913989.654780984,
                -41594974.26749641,
                -41271102.426548205,
                -40942372.66729347,
                -40608783.564378165,
                -40270333.731280066,
                -39927021.81985544,
                -39578846.51989073,
                -39225806.55866059,
                -38867900.70048869,
                -38505127.74631624,
                -38137486.53327459,
                -37764975.934260234,
                -37387594.85751937,
                -37005342.24623179,
                -36618217.07810502,
                -36226218.364967994,
                -35829345.15237361,
                -35427596.51920241,
                -35020971.577272564,
                -34609469.47095455,
                -34193089.37678883,
                -33771830.503109515,
                -33345692.089670055,
                -32914673.40727543,
                -32478773.757417858,
                -32037992.4719143,
                -31592328.9125537,
                -31141782.47073961,
                -30686352.567147497,
                -30226038.651375063,
                -29760840.20160454,
                -149754604.05101642,
                -149731749.81168228,
                -149708905.881834,
                -149686072.25527167,
                -149663248.92580009,
                -149640435.88722867,
                -149617633.13337162,
                -149594840.65804774,
                -149572058.45508042,
                -149549286.51829767,
                -149526524.84153217,
                -149503773.41862133,
                -149481032.243407,
                -149458301.3097358,
                -149435580.6114588,
                -149412870.14243203,
                -149390169.89651567,
                -149367479.86757484,
                -149344800.04947922,
                -149322130.43610284,
                -149299471.0213247,
                -149276821.79902807,
                -149254182.76310098,
                -149231553.90743592,
                -149208935.22593015,
                -149186326.7124853,
                -149163728.3610076,
                -149141140.165408,
                -149118562.11960188,
                -149095994.21750915,
                -149073436.45305434,
                -149050888.82016647,
                -149028351.31277913,
                -149005823.92483068,
                -148983306.65026352,
                -148960799.48302498,
                -148938302.41706684,
                -148915815.44634524,
                -148893338.56482106,
                -148870871.76645955,
                -148848415.04523054,
                -148825968.39510828,
                -148803531.81007168,
                -148781105.28410405,
                -148758688.81119308,
                -148736282.38533118,
                -148713886.0005152,
                -148691499.65074623,
                -148669123.33003008,
                -148646757.03237712,
                -148624400.75180185,
                -148602054.4823236,
                -148579718.21796584,
                -148557391.95275673,
                -148535075.68072885,
                -148512769.39591923,
                -148490473.09236908,
                -148468186.7641246,
                -148445910.40523586,
                -148423644.0097577,
                -148401387.57174936,
                -148379141.08527437,
                -148356904.54440084,
                -148334677.94320115,
                -148312461.27575228,
                -148290254.53613538,
                -148268057.71843624,
                -148245870.816745,
                -148223693.825156,
                -148201526.73776823,
                -148179369.54868498,
                -148157222.2520138,
                -148135084.8418669,
                -148112957.3123607,
                -148090839.6576159,
                -148068731.8717578,
                -148046633.9489158,
                -148024545.883224,
                -148002467.66882053,
                -147980399.29984805,
                -147958340.77045363,
                -147936292.0747885,
                -147914253.20700833,
                -147892224.16127324,
                -147870204.93174756,
                -147848195.5126,
                -147826195.8980034,
                -147804206.08213535,
                -147782226.0591773,
                -147760255.82331538,
                -147738295.36873987,
                -147716344.6896453,
                -147694403.7802306,
                -147672472.63469914,
                -147650551.24725825,
                -147628639.61211982,
                -147606737.7235,
                -147584845.575619,
                -147562963.16270173,
                -147541090.47897696,
            ],
            expected_sigz=[
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                1349227.729224144,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633991,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633986,
                -13363623.471633991,
                -13363623.471633988,
                -13363623.471633991,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633986,
                -13363623.471633986,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633986,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                -13363623.471633988,
                167753504.69823205,
                167753504.69823202,
                167753504.69823202,
                167753504.698232,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823205,
                167753504.69823202,
                167753504.69823205,
                167753504.69823205,
                167753504.698232,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.698232,
                167753504.69823202,
                167753504.69823202,
                167753504.69823205,
                167753504.69823202,
                167753504.69823202,
                167753504.698232,
                167753504.69823205,
                167753504.69823202,
                167753504.69823202,
                167753504.698232,
                167753504.69823202,
                167753504.69823205,
                167753504.698232,
                167753504.69823202,
                167753504.698232,
                167753504.69823205,
                167753504.698232,
                167753504.698232,
                167753504.69823205,
                167753504.698232,
                167753504.698232,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823205,
                167753504.69823202,
                167753504.698232,
                167753504.69823202,
                167753504.69823205,
                167753504.698232,
                167753504.69823202,
                167753504.698232,
                167753504.69823205,
                167753504.69823205,
                167753504.69823202,
                167753504.69823205,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.698232,
                167753504.698232,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.698232,
                167753504.69823202,
                167753504.69823205,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.698232,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823205,
                167753504.698232,
                167753504.69823205,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.698232,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823205,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                167753504.69823202,
                112188988.57039726,
                112191734.9797165,
                112199968.41119711,
                112213680.20105855,
                112232861.72655964,
                112257504.40573992,
                112287599.69716345,
                112323139.09966344,
                112364114.15208875,
                112410516.4330553,
                112462337.56069489,
                112519569.19241038,
                112582203.02462852,
                112650230.79255977,
                112723644.26995425,
                112802435.26886454,
                112886595.6394075,
                112976117.26952855,
                113070992.084768,
                113171212.04803002,
                113276769.1593505,
                113387655.45567037,
                113503863.01060896,
                113625383.93423767,
                113752210.37285809,
                113884334.50877978,
                114021748.56010136,
                114164444.78049125,
                114312415.4589722,
                114465652.91970573,
                114624149.52177952,
                114787897.65899567,
                114956889.75966033,
                115131118.28637597,
                115310575.73583439,
                115495254.63861054,
                115685147.55896066,
                115880247.09461819,
                116080545.8765946,
                116286036.56897861,
                116496711.86874084,
                116712564.50553493,
                116933587.24150482,
                117159772.87109038,
                117391114.22083598,
                117627604.14920005,
                117869235.54636602,
                118116001.3340544,
                118367894.4653366,
                118624907.92445098,
                118887034.72661762,
                119154267.91785805,
                119426600.57481246,
                119704025.80456218,
                119986536.74444973,
                120274126.5619036,
                120566788.45426103,
                120864515.64859442,
                121167301.40153939,
                121475138.99912035,
                121788021.75658281,
                122105943.01822218,
                122428896.15721622,
                122756874.5754587,
                123089871.70339227,
                123427880.99984638,
                123770895.95187148,
                124118910.07457829,
                124471916.91097689,
                124829910.0318159,
                125192883.03542532,
                125560829.54755741,
                125933743.22123154,
                126311617.73657821,
                126694446.80068503,
                127082224.14744382,
                127474943.53739801,
                127872598.75759272,
                128275183.6214241,
                128682691.96849015,
                129095117.66444443,
                129512454.60084756,
                129934696.69502331,
                130361837.88991223,
                130793872.15392943,
                131230793.4808207,
                131672595.88952184,
                132119273.42401777,
                132570820.15320235,
                133027230.17074038,
                133488497.59492895,
                133954616.56856179,
                134425581.25879264,
                134901385.85700023,
                135382024.5786556,
                135867491.66318652,
                136357781.373849,
                136852887.99759188,
                137352805.84493017,
                137857529.24981353,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
                416388238.13895464,
            ],
            expected_str_r=[
                0.0004316487841185448,
                0.0004270990208598772,
                0.0004225817205267804,
                0.0004180965750180727,
                0.0004136432798791019,
                0.0004092215342500769,
                0.0004048310408152506,
                0.00040047150575294246,
                0.00039614263868637813,
                0.0003918441526353354,
                0.0003875757639685784,
                0.00038333719235707133,
                0.000379128160727949,
                0.00037494839521923435,
                0.0003707976251352919,
                0.000366675582902998,
                0.00036258200402861714,
                0.0003585166270553742,
                0.00035447919352170386,
                0.00035046944792016933,
                0.0003464871376570386,
                0.0003425320130125021,
                0.00033860382710152294,
                0.0003347023358353104,
                0.00033082729788340025,
                0.0003269784746363335,
                0.00032315563016892604,
                0.0003193585312041118,
                0.0003155869470773554,
                0.0003118406497016194,
                0.0003081194135328806,
                0.00030442301553618243,
                0.00030075123515221317,
                0.0002971038542644058,
                0.00029348065716654464,
                0.000289881430530871,
                0.0002863059633766834,
                0.0002827540470394181,
                0.0002792254751402025,
                0.00027572004355587855,
                0.00027223755038948013,
                0.00026877779594116193,
                0.00026534058267956863,
                0.0002619257152136421,
                0.0002585330002648508,
                0.000255162246639841,
                0.0002518132652035008,
                0.00024848586885242794,
                0.00024517987248879496,
                0.00024189509299460876,
                0.00023863134920635242,
                0.0002353884618900039,
                0.0002321662537164308,
                0.00022896454923714634,
                0.0002257831748604264,
                0.00022262195882777856,
                0.0002194807311907608,
                0.00021635932378813742,
                0.00021325757022337525,
                0.00021017530584246556,
                0.00020711236771207365,
                0.0002040685945980048,
                0.00020104382694398663,
                0.0001980379068507581,
                0.00019505067805546128,
                0.00019208198591133476,
                0.0001891316773676967,
                0.00018619960095021858,
                0.00018328560674148262,
                0.0001803895463618183,
                0.00017751127295041215,
                0.0001746506411466909,
                0.0001718075070719665,
                0.00016898172831134603,
                0.00016617316389589516,
                0.0001633816742850591,
                0.00016060712134932718,
                0.00015784936835314773,
                0.0001551082799380797,
                0.0001523837221061834,
                0.00014967556220364295,
                0.0001469836689046212,
                0.00014430791219533737,
                0.0001416481633583677,
                0.0001390042949571664,
                0.0001363761808208003,
                0.0001337636960288956,
                0.00013116671689679252,
                0.00012858512096090749,
                0.00012601878696429374,
                0.00012346759484240422,
                0.00012093142570904713,
                0.00011841016184253601,
                0.00011590368667202847,
                0.00011341188476405196,
                0.00011093464180921301,
                0.00010847184460908738,
                0.00010602338106328847,
                0.00010358914015671009,
                0.00010116901194694335,
                -0.007844641345907013,
                -0.007843468882870921,
                -0.007842297035403567,
                -0.007841125803074112,
                -0.007839955185452089,
                -0.007838785182107413,
                -0.007837615792610372,
                -0.007836447016531634,
                -0.007835278853442233,
                -0.007834111302913588,
                -0.007832944364517484,
                -0.007831778037826087,
                -0.007830612322411935,
                -0.007829447217847934,
                -0.007828282723707372,
                -0.007827118839563902,
                -0.007825955564991553,
                -0.007824792899564725,
                -0.007823630842858188,
                -0.007822469394447085,
                -0.007821308553906928,
                -0.007820148320813599,
                -0.007818988694743352,
                -0.00781782967527281,
                -0.007816671261978962,
                -0.007815513454439168,
                -0.007814356252231158,
                -0.007813199654933026,
                -0.0078120436621232385,
                -0.007810888273380622,
                -0.0078097334882843775,
                -0.007808579306414067,
                -0.007807425727349619,
                -0.007806272750671331,
                -0.007805120375959864,
                -0.00780396860279624,
                -0.007802817430761852,
                -0.007801666859438453,
                -0.00780051688840816,
                -0.007799367517253454,
                -0.0077982187455571785,
                -0.007797070572902541,
                -0.00779592299887311,
                -0.007794776023052815,
                -0.007793629645025947,
                -0.007792483864377161,
                -0.0077913386806914686,
                -0.0077901940935542445,
                -0.0077890501025512215,
                -0.007787906707268493,
                -0.007786763907292512,
                -0.007785621702210089,
                -0.007784480091608394,
                -0.007783339075074955,
                -0.007782198652197658,
                -0.007781058822564742,
                -0.0077799195857648105,
                -0.0077787809413868165,
                -0.007777642889020076,
                -0.007776505428254254,
                -0.007775368558679373,
                -0.007774232279885816,
                -0.007773096591464312,
                -0.007771961493005953,
                -0.007770826984102176,
                -0.007769693064344779,
                -0.007768559733325911,
                -0.007767426990638071,
                -0.007766294835874115,
                -0.0077651632686272475,
                -0.007764032288491027,
                -0.007762901895059363,
                -0.007761772087926515,
                -0.007760642866687093,
                -0.0077595142309360605,
                -0.007758386180268728,
                -0.007757258714280755,
                -0.007756131832568153,
                -0.00775500553472728,
                -0.007753879820354843,
                -0.0077527546890478975,
                -0.007751630140403849,
                -0.0077505061740204455,
                -0.007749382789495786,
                -0.007748259986428314,
                -0.007747137764416823,
                -0.007746016123060447,
                -0.007744895061958668,
                -0.007743774580711317,
                -0.007742654678918564,
                -0.007741535356180926,
                -0.007740416612099265,
                -0.007739298446274788,
                -0.007738180858309041,
                -0.007737063847803918,
                -0.007735947414361651,
                -0.00773483155758482,
                -0.007733716277076342,
                -0.007732601572439478,
                -0.0077314874432778315,
                -7.425189848196039e-05,
                -7.647771733042263e-05,
                -7.869101941624279e-05,
                -8.089189841351755e-05,
                -8.308044712167185e-05,
                -8.525675747524065e-05,
                -8.742092055352463e-05,
                -8.957302659011999e-05,
                -9.17131649823225e-05,
                -9.384142430041327e-05,
                -9.595789229682142e-05,
                -9.806265591517081e-05,
                -0.00010015580129920786,
                -0.00010223741380161703,
                -0.00010430757799272182,
                -0.00010636637766907428,
                -0.00010841389586193511,
                -0.00011045021484564583,
                -0.00011247541614589315,
                -0.00011448958054786928,
                -0.0001164927881043287,
                -0.0001184851181435417,
                -0.00012046664927714903,
                -0.00012243745940791599,
                -0.00012439762573738936,
                -0.00012634722477345795,
                -0.00012828633233781831,
                -0.00013021502357334677,
                -0.00013213337295137886,
                -0.0001340414542788991,
                -0.00013593934070563944,
                -0.00013782710473109053,
                -0.00013970481821142552,
                -0.00014157255236633786,
                -0.0001434303777857949,
                -0.0001452783644367077,
                -0.00014711658166951809,
                -0.00014894509822470597,
                -0.00015076398223921557,
                -0.00015257330125280257,
                -0.00015437312221430462,
                -0.00015616351148783425,
                -0.00015794453485889668,
                -0.00015971625754043306,
                -0.00016147874417878907,
                -0.00016323205885961328,
                -0.00016497626511368247,
                -0.00016671142592265603,
                -0.00016843760372476258,
                -0.00017015486042041585,
                -0.00017186325737776555,
                -0.00017356285543817845,
                -0.00017525371492165692,
                -0.00017693589563218964,
                -0.00017860945686304087,
                -0.0001802744574019747,
                -0.0001819309555364177,
                -0.00018357900905856083,
                -0.00018521867527039934,
                -0.00018685001098871457,
                -0.00018847307254999562,
                -0.00019008791581530307,
                -0.00019169459617507557,
                -0.00019329316855387998,
                -0.00019488368741510511,
                -0.0001964662067656008,
                -0.00019804078016026307,
                -0.0001996074607065638,
                -0.00020116630106903035,
                -0.00020271735347367073,
                -0.00020426066971234802,
                -0.00020579630114710384,
                -0.0002073242987144321,
                -0.00020884471292950273,
                -0.00021035759389033653,
                -0.00021186299128193262,
                -0.00021336095438034743,
                -0.0002148515320567273,
                -0.0002163347727812942,
                -0.0002178107246272868,
                -0.0002192794352748553,
                -0.00022074095201491198,
                -0.0002221953217529388,
                -0.0002236425910127502,
                -0.00022508280594021423,
                -0.0002265160123069305,
                -0.0002279422555138673,
                -0.0002293615805949569,
                -0.0002307740322206509,
                -0.00023217965470143445,
                -0.00023357849199130183,
                -0.00023497058769119278,
                -0.00023635598505238975,
                -0.00023773472697987747,
                -0.0002391068560356651,
                -0.0002404724144420706,
                -0.0002418314440849694,
                -0.00024318398651700557,
                -0.00024453008296076886,
                -0.00024586977431193544,
                -0.0012295622544535054,
                -0.001229211594722353,
                -0.0012286101793795104,
                -0.0012277587953547834,
                -0.00122665822400144,
                -0.0012253092411458973,
                -0.0012237126171368911,
                -0.0012218691168941506,
                -0.0012197794999565437,
                -0.0012174445205296687,
                -0.0012148649275330343,
                -0.0012120414646466387,
                -0.001208974870357171,
                -0.0012056658780035953,
                -0.0012021152158223952,
                -0.0011983236069922382,
                -0.0011942917696782175,
                -0.0011900204170756455,
                -0.001185510257453361,
                -0.0011807619941965818,
                -0.001175776325849375,
                -0.0011705539461565966,
                -0.0011650955441054492,
                -0.0011594018039666287,
                -0.0011534734053349905,
                -0.0011473110231698506,
                -0.0011409153278348228,
                -0.0011342869851372975,
                -0.001127426656367473,
                -0.0011203349983370232,
                -0.0011130126634173335,
                -0.0011054602995773757,
                -0.001097678550421196,
                -0.0010896680552249952,
                -0.0010814294489738575,
                -0.0010729633623981182,
                -0.0010642704220093098,
                -0.0010553512501358187,
                -0.0010462064649581074,
                -0.0010368366805436659,
                -0.001027242506881503,
                -0.0010174245499164152,
                -0.00100738341158281,
                -0.0009971196898382452,
                -0.0009866339786966207,
                -0.0009759268682610304,
                -0.000964998944756297,
                -0.000953850790561187,
                -0.0009424829842402942,
                -0.0009308961005756023,
                -0.0009190907105977701,
                -0.0009070673816170494,
                -0.0008948266772539748,
                -0.0008823691574696606,
                -0.0008696953785958863,
                -0.0008568058933648142,
                -0.0008437012509384773,
                -0.0008303819969379465,
                -0.0008168486734721799,
                -0.0008031018191666899,
                -0.0007891419691918068,
                -0.0007749696552907727,
                -0.0007605854058074998,
                -0.0007459897457140734,
                -0.0007311831966380244,
                -0.0007161662768892525,
                -0.0007009395014868042,
                -0.0006855033821852989,
                -0.0006698584275011368,
                -0.0006540051427384719,
                -0.0006379440300148893,
                -0.0006216755882868963,
                -0.0006052003133751101,
                -0.0005885186979892466,
                -0.0005716312317528597,
                -0.0005545384012278325,
                -0.0005372406899386568,
                -0.0005197385783964569,
                -0.0005020325441228089,
                -0.00048412306167334144,
                -0.0004660106026610602,
                -0.00044769563577952854,
                -0.00042917862682575683,
                -0.0004104600387229445,
                -0.0003915403315429419,
                -0.00037241996252856827,
                -0.0003530993861156598,
                -0.0003335790539549403,
                -0.0003138594149337106,
                -0.00029394091519727855,
                -0.00027382399817025877,
                -0.00025350910457759806,
                -0.0002329966724654745,
                -0.00021228713722196905,
                -0.00019138093159752948,
                -0.00017027848572532026,
                -0.0001489802271412412,
                -0.00012748658080396238,
                -0.00010579796911457139,
                -8.391481193618206e-05,
                -8.989289545790189e-05,
                -9.003782478050838e-05,
                -9.01826887258877e-05,
                -9.032748733335647e-05,
                -9.047222064220066e-05,
                -9.061688869167772e-05,
                -9.076149152101499e-05,
                -9.09060291694103e-05,
                -9.105050167603255e-05,
                -9.119490908002085e-05,
                -9.133925142048476e-05,
                -9.148352873650494e-05,
                -9.16277410671323e-05,
                -9.177188845138894e-05,
                -9.191597092826707e-05,
                -9.205998853672976e-05,
                -9.220394131571142e-05,
                -9.234782930411643e-05,
                -9.249165254082069e-05,
                -9.263541106467068e-05,
                -9.277910491448338e-05,
                -9.292273412904749e-05,
                -9.306629874712192e-05,
                -9.32097988074365e-05,
                -9.335323434869292e-05,
                -9.349660540956264e-05,
                -9.36399120286894e-05,
                -9.378315424468696e-05,
                -9.39263320961405e-05,
                -9.406944562160674e-05,
                -9.42124948596128e-05,
                -9.435547984865773e-05,
                -9.449840062721119e-05,
                -9.46412572337139e-05,
                -9.47840497065788e-05,
                -9.492677808418899e-05,
                -9.506944240489932e-05,
                -9.521204270703616e-05,
                -9.535457902889668e-05,
                -9.549705140875025e-05,
                -9.563945988483682e-05,
                -9.578180449536786e-05,
                -9.592408527852691e-05,
                -9.606630227246808e-05,
                -9.620845551531804e-05,
                -9.635054504517396e-05,
                -9.649257090010486e-05,
                -9.663453311815196e-05,
                -9.677643173732739e-05,
                -9.691826679561461e-05,
                -9.706003833096985e-05,
                -9.720174638131992e-05,
                -9.734339098456412e-05,
                -9.748497217857295e-05,
                -9.762649000118868e-05,
                -9.776794449022579e-05,
                -9.790933568347009e-05,
                -9.805066361867935e-05,
                -9.819192833358348e-05,
                -9.833312986588376e-05,
                -9.847426825325395e-05,
                -9.861534353333928e-05,
                -9.875635574375669e-05,
                -9.889730492209612e-05,
                -9.903819110591828e-05,
                -9.917901433275695e-05,
                -9.931977464011746e-05,
                -9.946047206547682e-05,
                -9.960110664628506e-05,
                -9.974167841996363e-05,
                -9.988218742390623e-05,
                -0.00010002263369547924,
                -0.00010016301727202041,
                -0.00010030333819084056,
                -0.00010044359648922214,
                -0.00010058379220441997,
                -0.00010072392537366175,
                -0.00010086399603414652,
                -0.00010100400422304658,
                -0.00010114394997750599,
                -0.00010128383333464144,
                -0.00010142365433154235,
                -0.00010156341300526989,
                -0.00010170310939285824,
                -0.00010184274353131398,
                -0.0001019823154576159,
                -0.00010212182520871597,
                -0.00010226127282153823,
                -0.00010240065833297911,
                -0.00010253998177990853,
                -0.00010267924319916789,
                -0.00010281844262757242,
                -0.00010295758010190925,
                -0.00010309665565893822,
                -0.00010323566933539258,
                -0.00010337462116797778,
                -0.00010351351119337193,
                -0.00010365233944822672,
                -0.00010379110596916568,
                -0.0001039298107927861,
            ],
            expected_str_t=[
                -0.0014767839323287615,
                -0.0014722341690700939,
                -0.001467716868736997,
                -0.0014632317232282895,
                -0.0014587784280893185,
                -0.0014543566824602935,
                -0.0014499661890254672,
                -0.0014456066539631591,
                -0.0014412777868965947,
                -0.0014369793008455522,
                -0.001432710912178795,
                -0.001428472340567288,
                -0.0014242633089381657,
                -0.0014200835434294511,
                -0.0014159327733455087,
                -0.0014118107311132146,
                -0.0014077171522388337,
                -0.001403651775265591,
                -0.0013996143417319205,
                -0.0013956045961303859,
                -0.0013916222858672551,
                -0.0013876671612227187,
                -0.0013837389753117396,
                -0.001379837484045527,
                -0.001375962446093617,
                -0.0013721136228465503,
                -0.0013682907783791427,
                -0.0013644936794143284,
                -0.001360722095287572,
                -0.001356975797911836,
                -0.0013532545617430973,
                -0.001349558163746399,
                -0.0013458863833624298,
                -0.0013422390024746223,
                -0.0013386158053767614,
                -0.0013350165787410877,
                -0.0013314411115869,
                -0.0013278891952496347,
                -0.0013243606233504191,
                -0.0013208551917660953,
                -0.0013173726985996968,
                -0.0013139129441513786,
                -0.0013104757308897853,
                -0.0013070608634238587,
                -0.0013036681484750673,
                -0.0013002973948500578,
                -0.0012969484134137173,
                -0.0012936210170626447,
                -0.0012903150206990117,
                -0.0012870302412048254,
                -0.001283766497416569,
                -0.0012805236101002206,
                -0.0012773014019266476,
                -0.001274099697447363,
                -0.001270918323070643,
                -0.0012677571070379952,
                -0.0012646158794009776,
                -0.0012614944719983542,
                -0.0012583927184335918,
                -0.001255310454052682,
                -0.0012522475159222903,
                -0.0012492037428082216,
                -0.0012461789751542032,
                -0.0012431730550609746,
                -0.0012401858262656778,
                -0.0012372171341215514,
                -0.0012342668255779132,
                -0.0012313347491604352,
                -0.0012284207549516993,
                -0.0012255246945720348,
                -0.001222646421160629,
                -0.0012197857893569077,
                -0.001216942655282183,
                -0.0012141168765215627,
                -0.0012113083121061117,
                -0.0012085168224952756,
                -0.0012057422695595438,
                -0.0012029845165633644,
                -0.0012002434281482965,
                -0.0011975188703164,
                -0.0011948107104138596,
                -0.0011921188171148379,
                -0.001189443060405554,
                -0.0011867833115685844,
                -0.001184139443167383,
                -0.001181511329031017,
                -0.0011788988442391123,
                -0.0011763018651070093,
                -0.0011737202691711243,
                -0.0011711539351745104,
                -0.0011686027430526209,
                -0.0011660665739192638,
                -0.0011635453100527528,
                -0.001161038834882245,
                -0.0011585470329742685,
                -0.0011560697900194297,
                -0.001153606992819304,
                -0.001151158529273505,
                -0.0011487242883669267,
                -0.00114630416015716,
                -0.00114630416015717,
                -0.001147476623193261,
                -0.0011486484706606155,
                -0.0011498197029900714,
                -0.001150990320612093,
                -0.0011521603239567683,
                -0.0011533297134538089,
                -0.0011544984895325494,
                -0.0011556666526219496,
                -0.0011568342031505951,
                -0.001158001141546699,
                -0.0011591674682380954,
                -0.0011603331836522479,
                -0.0011614982882162481,
                -0.00116266278235681,
                -0.0011638266665002792,
                -0.0011649899410726286,
                -0.0011661526064994567,
                -0.0011673146632059945,
                -0.0011684761116170978,
                -0.0011696369521572545,
                -0.0011707971852505839,
                -0.00117195681132083,
                -0.0011731158307913725,
                -0.0011742742440852207,
                -0.0011754320516250134,
                -0.0011765892538330236,
                -0.0011777458511311556,
                -0.0011789018439409436,
                -0.00118005723268356,
                -0.0011812120177798046,
                -0.0011823661996501152,
                -0.0011835197787145634,
                -0.001184672755392851,
                -0.001185825130104318,
                -0.0011869769032679424,
                -0.00118812807530233,
                -0.0011892786466257284,
                -0.0011904286176560221,
                -0.0011915779888107278,
                -0.0011927267605070036,
                -0.0011938749331616406,
                -0.0011950225071910717,
                -0.0011961694830113678,
                -0.0011973158610382342,
                -0.0011984616416870209,
                -0.001199606825372714,
                -0.001200751412509938,
                -0.0012018954035129606,
                -0.001203038798795689,
                -0.0012041815987716697,
                -0.001205323803854093,
                -0.0012064654144557873,
                -0.0012076064309892265,
                -0.001208746853866525,
                -0.00120988668349944,
                -0.0012110259202993716,
                -0.0012121645646773651,
                -0.0012133026170441062,
                -0.0012144400778099284,
                -0.0012155769473848085,
                -0.0012167132261783663,
                -0.00121784891459987,
                -0.0012189840130582296,
                -0.0012201185219620054,
                -0.0012212524417194027,
                -0.0012223857727382713,
                -0.001223518515426111,
                -0.0012246506701900678,
                -0.0012257822374369342,
                -0.0012269132175731549,
                -0.00122804361100482,
                -0.0012291734181376674,
                -0.0012303026393770894,
                -0.0012314312751281215,
                -0.0012325593257954542,
                -0.0012336867917834273,
                -0.0012348136734960297,
                -0.001235939971336902,
                -0.0012370656857093394,
                -0.0012381908170162842,
                -0.0012393153656603327,
                -0.0012404393320437366,
                -0.0012415627165683956,
                -0.0012426855196358676,
                -0.001243807741647359,
                -0.0012449293830037352,
                -0.001246050444105514,
                -0.0012471709253528652,
                -0.0012482908271456185,
                -0.0012494101498832563,
                -0.0012505288939649163,
                -0.0012516470597893938,
                -0.0012527646477551416,
                -0.0012538816582602645,
                -0.0012549980917025312,
                -0.0012561139484793629,
                -0.0012572292289878403,
                -0.0012583439336247043,
                -0.00125945806278635,
                -0.0012594580627863432,
                -0.001257232243937881,
                -0.0012550189418520609,
                -0.0012528180628547862,
                -0.0012506295141466318,
                -0.001248453203793063,
                -0.001246289040714779,
                -0.0012441369346781837,
                -0.0012419967962859811,
                -0.0012398685369678904,
                -0.0012377520689714821,
                -0.0012356473053531328,
                -0.0012335541599690957,
                -0.0012314725474666865,
                -0.001229402383275582,
                -0.0012273435835992294,
                -0.0012252960654063685,
                -0.0012232597464226578,
                -0.0012212345451224106,
                -0.0012192203807204344,
                -0.001217217173163975,
                -0.0012152248431247618,
                -0.0012132433119911546,
                -0.0012112725018603877,
                -0.0012093123355309144,
                -0.0012073627364948456,
                -0.0012054236289304853,
                -0.0012034949376949569,
                -0.001201576588316925,
                -0.0011996685069894045,
                -0.001197770620562664,
                -0.0011958828565372132,
                -0.0011940051430568781,
                -0.001192137408901966,
                -0.0011902795834825087,
                -0.001188431596831596,
                -0.0011865933795987857,
                -0.0011847648630435976,
                -0.0011829459790290882,
                -0.0011811366600155012,
                -0.001179336839053999,
                -0.0011775464497804694,
                -0.001175765426409407,
                -0.0011739937037278706,
                -0.0011722312170895145,
                -0.0011704779024086904,
                -0.0011687336961546213,
                -0.0011669985353456476,
                -0.001165272357543541,
                -0.001163555100847888,
                -0.001161846703890538,
                -0.0011601471058301253,
                -0.0011584562463466466,
                -0.001156774065636114,
                -0.0011551005044052628,
                -0.0011534355038663288,
                -0.001151779005731886,
                -0.0011501309522097428,
                -0.0011484912859979044,
                -0.001146859950279589,
                -0.0011452368887183081,
                -0.0011436220454530007,
                -0.0011420153650932281,
                -0.0011404167927144236,
                -0.0011388262738531985,
                -0.001137243754502703,
                -0.0011356691811080406,
                -0.0011341025005617397,
                -0.0011325436601992732,
                -0.001130992607794633,
                -0.0011294492915559556,
                -0.0011279136601211998,
                -0.0011263856625538715,
                -0.001124865248338801,
                -0.0011233523673779672,
                -0.001121846969986371,
                -0.0011203490068879562,
                -0.0011188584292115764,
                -0.0011173751884870095,
                -0.0011158992366410167,
                -0.0011144305259934485,
                -0.0011129690092533916,
                -0.0011115146395153648,
                -0.0011100673702555534,
                -0.0011086271553280894,
                -0.0011071939489613733,
                -0.0011057677057544364,
                -0.0011043483806733468,
                -0.0011029359290476528,
                -0.0011015303065668692,
                -0.0011001314692770019,
                -0.001098739373577111,
                -0.0010973539762159138,
                -0.0010959752342884262,
                -0.0010946031052326386,
                -0.001093237546826233,
                -0.0010918785171833343,
                -0.001090525974751298,
                -0.0010891798783075348,
                -0.0010878401869563682,
                -0.0010878401869563682,
                -0.0010880646037029714,
                -0.0010882875565362711,
                -0.0010885086567706374,
                -0.0010887275194105366,
                -0.001088943763112741,
                -0.001089157010148928,
                -0.0010893668863687501,
                -0.001089573021163336,
                -0.0010897750474291312,
                -0.001089972601532264,
                -0.0010901653232732132,
                -0.001090352855851976,
                -0.001090534845833519,
                -0.0010907109431137567,
                -0.001090880800885783,
                -0.001091044075606585,
                -0.0010912004269640833,
                -0.0010913495178445616,
                -0.0010914910143004356,
                -0.0010916245855184584,
                -0.0010917499037881935,
                -0.0010918666444708802,
                -0.0010919744859686847,
                -0.0010920731096942394,
                -0.001092162200040554,
                -0.0010922414443512388,
                -0.0010923105328911069,
                -0.001092369158817046,
                -0.0010924170181492614,
                -0.001092453809742811,
                -0.0010924792352594614,
                -0.0010924929991398763,
                -0.001092494808576079,
                -0.0010924843734842402,
                -0.0010924614064777875,
                -0.0010924256228407384,
                -0.0010923767405014262,
                -0.0010923144800064294,
                -0.0010922385644948534,
                -0.0010921487196728166,
                -0.0010920446737883232,
                -0.0010919261576063004,
                -0.0010917929043839888,
                -0.0010916446498465588,
                -0.0010914811321630043,
                -0.0010913020919222998,
                -0.0010911072721098288,
                -0.0010908964180840507,
                -0.001090669277553402,
                -0.0010904256005535168,
                -0.0010901651394246,
                -0.0010898876487891427,
                -0.0010895928855297792,
                -0.0010892806087674752,
                -0.0010889505798398682,
                -0.0010886025622799183,
                -0.0010882363217947283,
                -0.001087851626244595,
                -0.0010874482456223538,
                -0.0010870259520328344,
                -0.00108658451967263,
                -0.0010861237248100373,
                -0.0010856433457652104,
                -0.001085143162890558,
                -0.0010846229585512797,
                -0.0010840825171062052,
                -0.0010835216248887573,
                -0.0010829400701881406,
                -0.001082337643230761,
                -0.0010817141361617697,
                -0.0010810693430268936,
                -0.0010804030597543711,
                -0.0010797150841371352,
                -0.0010790052158151588,
                -0.0010782732562579886,
                -0.0010775190087474744,
                -0.0010767422783606476,
                -0.0010759428719528197,
                -0.00107512059814084,
                -0.0010742752672864985,
                -0.0010734066914801666,
                -0.0010725146845245266,
                -0.0010715990619185574,
                -0.0010706596408416068,
                -0.0010696962401376914,
                -0.0010687086802999128,
                -0.0010676967834550612,
                -0.0010666603733483738,
                -0.0010655992753284453,
                -0.0010645133163323022,
                -0.001063402324870611,
                -0.0010622661310130633,
                -0.0010611045663739059,
                -0.0010599174640975825,
                -0.0010587046588446051,
                -0.0010574659867774545,
                -0.0010562012855467593,
                -0.0010549103942774788,
                -0.0010535931535553314,
                -0.0010535931535553312,
                -0.0010534482242327246,
                -0.0010533033602873453,
                -0.0010531585616798767,
                -0.0010530138283710326,
                -0.0010528691603215554,
                -0.001052724557492218,
                -0.0010525800198438228,
                -0.0010524355473372005,
                -0.0010522911399332123,
                -0.0010521467975927483,
                -0.0010520025202767282,
                -0.0010518583079461008,
                -0.0010517141605618442,
                -0.001051570078084966,
                -0.0010514260604765034,
                -0.0010512821076975216,
                -0.0010511382197091166,
                -0.0010509943964724124,
                -0.0010508506379485623,
                -0.0010507069440987497,
                -0.0010505633148841856,
                -0.0010504197502661113,
                -0.0010502762502057966,
                -0.0010501328146645403,
                -0.0010499894436036706,
                -0.0010498461369845436,
                -0.0010497028947685462,
                -0.0010495597169170926,
                -0.0010494166033916264,
                -0.0010492735541536204,
                -0.0010491305691645754,
                -0.0010489876483860219,
                -0.0010488447917795193,
                -0.0010487019993066543,
                -0.001048559270929044,
                -0.0010484166066083338,
                -0.001048274006306197,
                -0.0010481314699843365,
                -0.0010479889976044829,
                -0.0010478465891283963,
                -0.0010477042445178652,
                -0.0010475619637347062,
                -0.001047419746740765,
                -0.001047277593497915,
                -0.0010471355039680592,
                -0.0010469934781131283,
                -0.0010468515158950811,
                -0.0010467096172759058,
                -0.0010465677822176186,
                -0.0010464260106822633,
                -0.0010462843026319132,
                -0.001046142658028669,
                -0.0010460010768346601,
                -0.0010458595590120444,
                -0.0010457181045230074,
                -0.001045576713329763,
                -0.0010454353853945539,
                -0.0010452941206796496,
                -0.0010451529191473493,
                -0.0010450117807599792,
                -0.001044870705479894,
                -0.0010447296932694764,
                -0.001044588744091137,
                -0.0010444478579073148,
                -0.0010443070346804761,
                -0.0010441662743731156,
                -0.0010440255769477564,
                -0.001043884942366948,
                -0.0010437443705932695,
                -0.001043603861589327,
                -0.0010434634153177538,
                -0.0010433230317412126,
                -0.0010431827108223926,
                -0.001043042452524011,
                -0.0010429022568088132,
                -0.0010427621236395713,
                -0.0010426220529790866,
                -0.0010424820447901865,
                -0.001042342099035727,
                -0.0010422022156785916,
                -0.0010420623946816908,
                -0.0010419226360079631,
                -0.0010417829396203748,
                -0.0010416433054819191,
                -0.0010415037335556173,
                -0.0010413642238045172,
                -0.001041224776191695,
                -0.001041085390680254,
                -0.0010409460672333247,
                -0.0010408068058140652,
                -0.0010406676063856608,
                -0.0010405284689113238,
                -0.001040389393354295,
                -0.0010402503796778405,
                -0.0010401114278452553,
                -0.0010399725378198612,
                -0.0010398337095650064,
                -0.0010396949430440674,
                -0.001039556238220447,
            ],
            expected_str_z=[
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011587769306169113,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
                0.0011794772478124136,
            ],
            expected_r_deflect=[
                -0.0034441554869771383,
                -0.003441759595768685,
                -0.0034393890015235297,
                -0.0034370435239598584,
                -0.0034347229845048633,
                -0.003432427206274539,
                -0.0034301560140537674,
                -0.0034279092342766816,
                -0.0034256866950073033,
                -0.00342348822592046,
                -0.0034213136582829633,
                -0.00341916282493505,
                -0.0034170355602720897,
                -0.0034149317002265386,
                -0.0034128510822501465,
                -0.0034107935452964153,
                -0.0034087589298032917,
                -0.0034067470776761055,
                -0.0034047578322707393,
                -0.0034027910383770303,
                -0.003400846542202399,
                -0.003398924191355703,
                -0.0033970238348313086,
                -0.0033951453229933816,
                -0.003393288507560391,
                -0.0033914532415898185,
                -0.0033896393794630828,
                -0.00338784677687066,
                -0.0033860752907974062,
                -0.003384324779508078,
                -0.0033825951025330467,
                -0.003380886120654204,
                -0.0033791976958910544,
                -0.0033775296914869916,
                -0.003375881971895762,
                -0.0033742544027680994,
                -0.0033726468509385453,
                -0.0033710591844124376,
                -0.0033694912723530707,
                -0.003367942985069025,
                -0.0033664141940016657,
                -0.003364904771712798,
                -0.003363414591872488,
                -0.0033619435292470447,
                -0.0033604914596871514,
                -0.003359058260116154,
                -0.0033576438085185056,
                -0.003356247983928349,
                -0.0033548706664182585,
                -0.003353511737088118,
                -0.003352171078054145,
                -0.0033508485724380553,
                -0.0033495441043563637,
                -0.003348257558909824,
                -0.003346988822173,
                -0.003345737781183974,
                -0.003344504323934177,
                -0.003343288339358358,
                -0.0033420897173246714,
                -0.0033409083486248903,
                -0.0033397441249647486,
                -0.003338596938954397,
                -0.003337466684098982,
                -0.003336353254789341,
                -0.0033352565462928135,
                -0.0033341764547441697,
                -0.0033331128771366438,
                -0.003332065711313088,
                -0.0033310348559572266,
                -0.003330020210585025,
                -0.0033290216755361606,
                -0.0033280391519656,
                -0.0033270725418352784,
                -0.003326121747905882,
                -0.0033251866737287303,
                -0.003324267223637755,
                -0.00332336330274158,
                -0.003322474816915694,
                -0.0033216016727947223,
                -0.003320743777764784,
                -0.003319901039955951,
                -0.0033190733682347904,
                -0.0033182606721969985,
                -0.003317462862160127,
                -0.003316679849156387,
                -0.00331591154492555,
                -0.003315157861907923,
                -0.0033144187132374157,
                -0.0033136940127346847,
                -0.003312983674900359,
                -0.003312287614908349,
                -0.0033116057485992312,
                -0.0033109379924737116,
                -0.0033102842636861646,
                -0.0033096444800382527,
                -0.003309018559972614,
                -0.0033084064225666265,
                -0.003307807987526247,
                -0.0033072231751799176,
                -0.0033066519064725475,
                -0.0033066519064725765,
                -0.003310613550443823,
                -0.0033145746024174916,
                -0.003318535062704355,
                -0.0033224949316149888,
                -0.0033264542094597425,
                -0.003330412896548743,
                -0.0033343709931919104,
                -0.0033383284996989417,
                -0.0033422854163793186,
                -0.0033462417435423126,
                -0.0033501974814969648,
                -0.003354152630552113,
                -0.0033581071910163806,
                -0.0033620611631981594,
                -0.003366014547405644,
                -0.0033699673439468066,
                -0.0033739195531293977,
                -0.0033778711752609665,
                -0.003381822210648835,
                -0.0033857726596001145,
                -0.0033897225224217113,
                -0.0033936717994203017,
                -0.003397620490902359,
                -0.003401568597174144,
                -0.003405516118541691,
                -0.0034094630553108355,
                -0.0034134094077871973,
                -0.003417355176276171,
                -0.0034213003610829573,
                -0.0034252449625125237,
                -0.003429188980869642,
                -0.003433132416458871,
                -0.0034370752695845382,
                -0.0034410175405507858,
                -0.0034449592296615286,
                -0.003448900337220468,
                -0.003452840863531101,
                -0.0034567808088967183,
                -0.0034607201736203824,
                -0.0034646589580049642,
                -0.0034685971623531087,
                -0.003472534786967258,
                -0.0034764718321496496,
                -0.0034804082982022953,
                -0.003484344185427009,
                -0.003488279494125399,
                -0.003492214224598845,
                -0.0034961483771485374,
                -0.00350008195207545,
                -0.003504014949680343,
                -0.0035079473702637753,
                -0.003511879214126091,
                -0.0035158104815674295,
                -0.0035197411728877265,
                -0.0035236712883866955,
                -0.0035276008283638555,
                -0.0035315297931185156,
                -0.0035354581829497717,
                -0.003539385998156513,
                -0.003543313239037435,
                -0.0035472399058910017,
                -0.0035511659990154973,
                -0.003555091518708973,
                -0.003559016465269297,
                -0.0035629408389941195,
                -0.0035668646401808817,
                -0.003570787869126825,
                -0.0035747105261289876,
                -0.003578632611484193,
                -0.0035825541254890653,
                -0.0035864750684400284,
                -0.0035903954406332848,
                -0.003594315242364853,
                -0.0035982344739305273,
                -0.00360215313562591,
                -0.0036060712277464022,
                -0.003609988750587184,
                -0.003613905704443249,
                -0.0036178220896093786,
                -0.0036217379063801514,
                -0.003625653155049941,
                -0.003629567835912926,
                -0.00363348194926307,
                -0.0036373954953941457,
                -0.003641308474599708,
                -0.003645220887173126,
                -0.003649132733407558,
                -0.0036530440135959553,
                -0.003656954728031075,
                -0.00366086487700548,
                -0.0036647744608115044,
                -0.00366868347974131,
                -0.0036725919340868437,
                -0.003676499824139848,
                -0.0036804071501918764,
                -0.0036843139125342667,
                -0.0036882201114581663,
                -0.003692125747254525,
                -0.0036960308202140784,
                -0.0036960308202140593,
                -0.003696447351796846,
                -0.0036968761503905876,
                -0.003697317147076975,
                -0.003697770273452995,
                -0.003698235461626122,
                -0.0036987126442095675,
                -0.003699201754317577,
                -0.003699702725560783,
                -0.0037002154920416103,
                -0.003700739988349728,
                -0.003701276149557556,
                -0.003701823911215819,
                -0.003702383209349151,
                -0.003702953980451746,
                -0.00370353616148306,
                -0.003704129689863555,
                -0.0037047345034704952,
                -0.0037053505406337848,
                -0.0037059777401318543,
                -0.00370661604118759,
                -0.003707265383464306,
                -0.003707925707061766,
                -0.0037085969525122394,
                -0.003709279060776608,
                -0.0037099719732405097,
                -0.0037106756317105277,
                -0.003711389978410414,
                -0.003712114955977365,
                -0.0037128505074583208,
                -0.0037135965763063217,
                -0.0037143531063768896,
                -0.0037151200419244567,
                -0.003715897327598826,
                -0.003716684908441675,
                -0.003717482729883091,
                -0.00371829073773815,
                -0.0037191088782035225,
                -0.0037199370978541254,
                -0.003720775343639799,
                -0.003721623562882027,
                -0.003722481703270686,
                -0.0037233497128608337,
                -0.003724227540069525,
                -0.0037251151336726634,
                -0.003726012442801889,
                -0.003726919416941493,
                -0.0037278360059253687,
                -0.0037287621599339923,
                -0.0037296978294914336,
                -0.0037306429654623994,
                -0.0037315975190493084,
                -0.0037325614417893947,
                -0.0037335346855518372,
                -0.0037345172025349286,
                -0.0037355089452632624,
                -0.003736509866584955,
                -0.0037375199196688948,
                -0.003738539058002019,
                -0.0037395672353866182,
                -0.003740604405937667,
                -0.0037416505240801855,
                -0.003742705544546623,
                -0.003743769422374271,
                -0.003744842112902702,
                -0.003745923571771231,
                -0.003747013754916408,
                -0.003748112618569531,
                -0.0037492201192541865,
                -0.003750336213783812,
                -0.003751460859259287,
                -0.0037525940130665432,
                -0.003753735632874203,
                -0.0037548856766312366,
                -0.0037560441025646482,
                -0.0037572108691771786,
                -0.0037583859352450363,
                -0.003759569259815648,
                -0.0037607608022054307,
                -0.0037619605219975875,
                -0.003763168379039925,
                -0.00376438433344269,
                -0.0037656083455764305,
                -0.003766840376069875,
                -0.0037680803858078323,
                -0.003769328335929115,
                -0.0037705841878244795,
                -0.003771847903134589,
                -0.0037731194437479925,
                -0.003774398771799128,
                -0.0037756858496663425,
                -0.0037769806399699287,
                -0.0037782831055701853,
                -0.003779593209565495,
                -0.003780910915290417,
                -0.0037822361863138016,
                -0.003783568986436922,
                -0.003784909279691625,
                -0.003786257030338494,
                -0.00378761220286504,
                -0.00378761220286504,
                -0.0037944085660760814,
                -0.0038012022973116606,
                -0.0038079920125323945,
                -0.00381477633203371,
                -0.0038215538804152196,
                -0.0038283232865502057,
                -0.003835083183555538,
                -0.003841832208761871,
                -0.003848569003683834,
                -0.0038552922139908877,
                -0.003862000489478057,
                -0.0038686924840372694,
                -0.00387536685562855,
                -0.0038820222662519404,
                -0.003888657381919198,
                -0.0038952708726260883,
                -0.003901861412324752,
                -0.003908427678896362,
                -0.003914968354123923,
                -0.003921482123665571,
                -0.003927967677027777,
                -0.003934423707538923,
                -0.003940848912323231,
                -0.003947241992274724,
                -0.003953601652031513,
                -0.0039599265999502205,
                -0.003966215548080795,
                -0.00397246721214133,
                -0.003978680311493265,
                -0.003984853569116701,
                -0.003990985711585948,
                -0.003997075469045343,
                -0.004003121575185166,
                -0.004009122767217826,
                -0.004015077785854328,
                -0.004020985375280672,
                -0.004026844283134816,
                -0.0040326532604835165,
                -0.0040384110617996145,
                -0.004044116444939213,
                -0.004049768171119451,
                -0.004055365004896054,
                -0.004060905714141372,
                -0.00406638907002243,
                -0.004071813846979222,
                -0.00407717882270318,
                -0.004082482778115867,
                -0.004087724497347787,
                -0.0040929027677172915,
                -0.004098016379709932,
                -0.00410306412695762,
                -0.004108044806218331,
                -0.004112957217355607,
                -0.004117800163318575,
                -0.004122572450121809,
                -0.00412727288682567,
                -0.004131900285516557,
                -0.0041364534612873555,
                -0.0041409312322182945,
                -0.0041453324193575485,
                -0.004149655846702353,
                -0.0041539003411801,
                -0.004158064732629598,
                -0.0041621478537825895,
                -0.004166148540245168,
                -0.004170065630479705,
                -0.004173897965786626,
                -0.004177644390286417,
                -0.00418130375090191,
                -0.004184874897340417,
                -0.004188356682076385,
                -0.004191747960333831,
                -0.00419504759006916,
                -0.004198254431954011,
                -0.004201367349358261,
                -0.004204385208333197,
                -0.0042073068775947174,
                -0.004210131228506818,
                -0.0042128571350651305,
                -0.004215483473880495,
                -0.004218009124162905,
                -0.004220432967705254,
                -0.0042227538888675725,
                -0.004224970774561047,
                -0.004227082514232434,
                -0.004229087999848406,
                -0.004230986125880116,
                -0.004232775789287874,
                -0.004234455889505924,
                -0.0042360253284273565,
                -0.004237483010389088,
                -0.004238827842157023,
                -0.004240058732911337,
                -0.004241174594231711,
                -0.004242174340083015,
                -0.004243056886800643,
                -0.004243821153076487,
                -0.004244466059944489,
                -0.004244990530766704,
                -0.004244990530766704,
                -0.004245045055230685,
                -0.004245099667510806,
                -0.004245154367567456,
                -0.0042452091553610485,
                -0.00424526403085202,
                -0.004245318994000833,
                -0.0042453740447679685,
                -0.004245429183113938,
                -0.004245484408999273,
                -0.004245539722384528,
                -0.004245595123230282,
                -0.00424565061149714,
                -0.004245706187145728,
                -0.004245761850136696,
                -0.004245817600430718,
                -0.0042458734379884925,
                -0.0042459293627707395,
                -0.004245985374738204,
                -0.004246041473851655,
                -0.004246097660071883,
                -0.004246153933359705,
                -0.00424621029367596,
                -0.004246266740981509,
                -0.004246323275237238,
                -0.004246379896404057,
                -0.004246436604442898,
                -0.004246493399314717,
                -0.004246550280980494,
                -0.004246607249401232,
                -0.004246664304537957,
                -0.004246721446351719,
                -0.004246778674803589,
                -0.004246835989854666,
                -0.004246893391466068,
                -0.004246950879598938,
                -0.004247008454214442,
                -0.004247066115273769,
                -0.004247123862738131,
                -0.004247181696568765,
                -0.004247239616726929,
                -0.0042472976231739045,
                -0.004247355715870998,
                -0.004247413894779537,
                -0.004247472159860874,
                -0.004247530511076382,
                -0.004247588948387459,
                -0.0042476474717555255,
                -0.004247706081142027,
                -0.004247764776508428,
                -0.004247823557816219,
                -0.004247882425026913,
                -0.004247941378102045,
                -0.004248000417003175,
                -0.004248059541691883,
                -0.0042481187521297755,
                -0.0042481780482784785,
                -0.004248237430099643,
                -0.004248296897554941,
                -0.0042483564506060716,
                -0.00424841608921475,
                -0.00424847581334272,
                -0.004248535622951746,
                -0.004248595518003615,
                -0.004248655498460137,
                -0.0042487155642831455,
                -0.004248775715434496,
                -0.004248835951876067,
                -0.004248896273569759,
                -0.0042489566804774954,
                -0.004249017172561223,
                -0.0042490777497829115,
                -0.004249138412104552,
                -0.004249199159488159,
                -0.0042492599918957695,
                -0.004249320909289442,
                -0.00424938191163126,
                -0.004249442998883328,
                -0.004249504171007772,
                -0.004249565427966743,
                -0.004249626769722412,
                -0.004249688196236974,
                -0.004249749707472647,
                -0.00424981130339167,
                -0.004249872983956305,
                -0.0042499347491288374,
                -0.004249996598871573,
                -0.00425005853314684,
                -0.004250120551916993,
                -0.004250182655144402,
                -0.004250244842791467,
                -0.004250307114820605,
                -0.004250369471194258,
                -0.004250431911874887,
                -0.0042504944368249795,
                -0.004250557046007042,
                -0.004250619739383606,
                -0.004250682516917222,
                -0.004250745378570465,
                -0.004250808324305933,
            ],
            expected_rradius=[
                2.3322000000000003,
                2.3377800000000004,
                2.34336,
                2.3489400000000002,
                2.3545200000000004,
                2.3601,
                2.3656800000000002,
                2.3712600000000004,
                2.3768400000000005,
                2.38242,
                2.3880000000000003,
                2.3935800000000005,
                2.39916,
                2.4047400000000003,
                2.4103200000000005,
                2.4159,
                2.4214800000000003,
                2.4270600000000004,
                2.43264,
                2.4382200000000003,
                2.4438000000000004,
                2.44938,
                2.4549600000000003,
                2.4605400000000004,
                2.46612,
                2.4717000000000002,
                2.4772800000000004,
                2.4828600000000005,
                2.48844,
                2.4940200000000003,
                2.4996000000000005,
                2.50518,
                2.5107600000000003,
                2.5163400000000005,
                2.52192,
                2.5275000000000003,
                2.5330800000000004,
                2.53866,
                2.5442400000000003,
                2.5498200000000004,
                2.5554000000000006,
                2.5609800000000003,
                2.5665600000000004,
                2.5721400000000005,
                2.5777200000000002,
                2.5833000000000004,
                2.5888800000000005,
                2.59446,
                2.6000400000000004,
                2.6056200000000005,
                2.6112,
                2.6167800000000003,
                2.6223600000000005,
                2.6279400000000006,
                2.6335200000000003,
                2.6391000000000004,
                2.64468,
                2.6502600000000003,
                2.6558400000000004,
                2.6614200000000006,
                2.6670000000000003,
                2.6725800000000004,
                2.6781600000000005,
                2.6837400000000002,
                2.6893200000000004,
                2.6949000000000005,
                2.70048,
                2.7060600000000004,
                2.7116400000000005,
                2.71722,
                2.7228000000000003,
                2.7283800000000005,
                2.7339600000000006,
                2.7395400000000003,
                2.7451200000000004,
                2.7507,
                2.7562800000000003,
                2.7618600000000004,
                2.7674400000000006,
                2.7730200000000003,
                2.7786000000000004,
                2.7841800000000005,
                2.7897600000000002,
                2.7953400000000004,
                2.8009200000000005,
                2.8065000000000007,
                2.8120800000000004,
                2.8176600000000005,
                2.82324,
                2.8288200000000003,
                2.8344000000000005,
                2.8399800000000006,
                2.8455600000000003,
                2.8511400000000005,
                2.85672,
                2.8623000000000003,
                2.8678800000000004,
                2.8734600000000006,
                2.8790400000000007,
                2.8846200000000004,
                2.8846200000000004,
                2.8851250505050507,
                2.8856301010101015,
                2.886135151515152,
                2.886640202020202,
                2.887145252525253,
                2.8876503030303033,
                2.888155353535354,
                2.8886604040404045,
                2.889165454545455,
                2.8896705050505056,
                2.890175555555556,
                2.8906806060606063,
                2.891185656565657,
                2.8916907070707074,
                2.8921957575757578,
                2.8927008080808085,
                2.893205858585859,
                2.8937109090909097,
                2.89421595959596,
                2.8947210101010104,
                2.895226060606061,
                2.8957311111111115,
                2.896236161616162,
                2.8967412121212126,
                2.897246262626263,
                2.8977513131313133,
                2.898256363636364,
                2.8987614141414144,
                2.899266464646465,
                2.8997715151515155,
                2.900276565656566,
                2.9007816161616167,
                2.901286666666667,
                2.9017917171717174,
                2.902296767676768,
                2.9028018181818185,
                2.903306868686869,
                2.9038119191919196,
                2.90431696969697,
                2.9048220202020207,
                2.905327070707071,
                2.9058321212121214,
                2.906337171717172,
                2.9068422222222225,
                2.907347272727273,
                2.9078523232323237,
                2.908357373737374,
                2.9088624242424244,
                2.909367474747475,
                2.9098725252525255,
                2.9103775757575763,
                2.9108826262626266,
                2.911387676767677,
                2.9118927272727277,
                2.912397777777778,
                2.9129028282828284,
                2.913407878787879,
                2.9139129292929296,
                2.91441797979798,
                2.9149230303030307,
                2.915428080808081,
                2.915933131313132,
                2.916438181818182,
                2.9169432323232325,
                2.9174482828282833,
                2.9179533333333336,
                2.918458383838384,
                2.9189634343434347,
                2.919468484848485,
                2.9199735353535354,
                2.920478585858586,
                2.9209836363636366,
                2.9214886868686873,
                2.9219937373737377,
                2.922498787878788,
                2.923003838383839,
                2.923508888888889,
                2.9240139393939395,
                2.9245189898989903,
                2.9250240404040406,
                2.925529090909091,
                2.9260341414141418,
                2.926539191919192,
                2.927044242424243,
                2.9275492929292932,
                2.9280543434343436,
                2.9285593939393944,
                2.9290644444444447,
                2.929569494949495,
                2.930074545454546,
                2.930579595959596,
                2.9310846464646465,
                2.9315896969696973,
                2.9320947474747476,
                2.9325997979797984,
                2.9331048484848488,
                2.933609898989899,
                2.93411494949495,
                2.9346200000000002,
                2.9346200000000002,
                2.9401467943734065,
                2.945673588746813,
                2.951200383120219,
                2.956727177493626,
                2.962253971867032,
                2.9677807662404385,
                2.973307560613845,
                2.978834354987251,
                2.9843611493606574,
                2.9898879437340637,
                2.99541473810747,
                3.000941532480877,
                3.006468326854283,
                3.0119951212276894,
                3.0175219156010957,
                3.023048709974502,
                3.0285755043479083,
                3.0341022987213147,
                3.039629093094721,
                3.0451558874681277,
                3.050682681841534,
                3.0562094762149403,
                3.0617362705883466,
                3.067263064961753,
                3.0727898593351592,
                3.0783166537085656,
                3.0838434480819723,
                3.0893702424553786,
                3.094897036828785,
                3.1004238312021912,
                3.1059506255755975,
                3.111477419949004,
                3.11700421432241,
                3.1225310086958165,
                3.128057803069223,
                3.1335845974426295,
                3.139111391816036,
                3.144638186189442,
                3.1501649805628484,
                3.1556917749362547,
                3.161218569309661,
                3.1667453636830674,
                3.172272158056474,
                3.1777989524298804,
                3.1833257468032867,
                3.188852541176693,
                3.1943793355500993,
                3.1999061299235056,
                3.205432924296912,
                3.2109597186703187,
                3.216486513043725,
                3.2220133074171313,
                3.2275401017905376,
                3.233066896163944,
                3.2385936905373502,
                3.2441204849107566,
                3.2496472792841633,
                3.2551740736575696,
                3.260700868030976,
                3.2662276624043822,
                3.2717544567777885,
                3.277281251151195,
                3.282808045524601,
                3.2883348398980075,
                3.2938616342714138,
                3.2993884286448205,
                3.304915223018227,
                3.310442017391633,
                3.3159688117650394,
                3.3214956061384457,
                3.327022400511852,
                3.3325491948852584,
                3.338075989258665,
                3.3436027836320714,
                3.3491295780054777,
                3.354656372378884,
                3.3601831667522903,
                3.3657099611256966,
                3.371236755499103,
                3.3767635498725097,
                3.382290344245916,
                3.3878171386193223,
                3.3933439329927286,
                3.398870727366135,
                3.4043975217395412,
                3.4099243161129476,
                3.415451110486354,
                3.42097790485976,
                3.426504699233167,
                3.4320314936065732,
                3.4375582879799795,
                3.443085082353386,
                3.448611876726792,
                3.4541386711001985,
                3.4596654654736048,
                3.4651922598470115,
                3.470719054220418,
                3.476245848593824,
                3.4817726429672304,
                3.4817726429672304,
                3.48730080287761,
                3.49282896278799,
                3.4983571226983696,
                3.503885282608749,
                3.5094134425191292,
                3.514941602429509,
                3.5204697623398884,
                3.5259979222502684,
                3.531526082160648,
                3.5370542420710276,
                3.5425824019814076,
                3.548110561891787,
                3.5536387218021668,
                3.559166881712547,
                3.5646950416229264,
                3.570223201533306,
                3.575751361443686,
                3.5812795213540656,
                3.5868076812644456,
                3.592335841174825,
                3.5978640010852048,
                3.603392160995585,
                3.6089203209059644,
                3.614448480816344,
                3.619976640726724,
                3.6255048006371036,
                3.631032960547483,
                3.636561120457863,
                3.6420892803682428,
                3.6476174402786223,
                3.6531456001890024,
                3.658673760099382,
                3.6642019200097615,
                3.6697300799201416,
                3.675258239830521,
                3.6807863997409007,
                3.6863145596512807,
                3.6918427195616603,
                3.69737087947204,
                3.70289903938242,
                3.7084271992927995,
                3.713955359203179,
                3.719483519113559,
                3.7250116790239387,
                3.7305398389343187,
                3.7360679988446983,
                3.741596158755078,
                3.747124318665458,
                3.7526524785758375,
                3.758180638486217,
                3.763708798396597,
                3.7692369583069767,
                3.7747651182173563,
                3.7802932781277363,
                3.785821438038116,
                3.7913495979484955,
                3.7968777578588755,
                3.802405917769255,
                3.8079340776796347,
                3.8134622375900147,
                3.8189903975003943,
                3.824518557410774,
                3.830046717321154,
                3.8355748772315335,
                3.841103037141913,
                3.846631197052293,
                3.8521593569626726,
                3.8576875168730522,
                3.8632156767834323,
                3.868743836693812,
                3.8742719966041914,
                3.8798001565145714,
                3.885328316424951,
                3.8908564763353306,
                3.8963846362457106,
                3.90191279615609,
                3.90744095606647,
                3.91296911597685,
                3.9184972758872294,
                3.924025435797609,
                3.929553595707989,
                3.9350817556183686,
                3.940609915528748,
                3.946138075439128,
                3.951666235349508,
                3.9571943952598874,
                3.9627225551702674,
                3.968250715080647,
                3.973778874991027,
                3.9793070349014066,
                3.984835194811786,
                3.9903633547221657,
                3.9958915146325458,
                4.001419674542926,
                4.006947834453305,
                4.012475994363685,
                4.018004154274065,
                4.023532314184444,
                4.029060474094824,
                4.029060474094824,
                4.029666534700885,
                4.030272595306945,
                4.030878655913006,
                4.031484716519066,
                4.032090777125127,
                4.032696837731188,
                4.033302898337248,
                4.033908958943309,
                4.03451501954937,
                4.03512108015543,
                4.035727140761491,
                4.036333201367551,
                4.036939261973612,
                4.037545322579673,
                4.038151383185733,
                4.038757443791794,
                4.039363504397854,
                4.039969565003915,
                4.040575625609976,
                4.041181686216036,
                4.041787746822097,
                4.042393807428158,
                4.042999868034218,
                4.043605928640279,
                4.044211989246339,
                4.0448180498524,
                4.045424110458461,
                4.046030171064521,
                4.046636231670582,
                4.047242292276642,
                4.047848352882703,
                4.048454413488764,
                4.049060474094824,
                4.049666534700885,
                4.050272595306946,
                4.050878655913006,
                4.051484716519067,
                4.052090777125127,
                4.052696837731188,
                4.0533028983372485,
                4.053908958943309,
                4.0545150195493695,
                4.05512108015543,
                4.0557271407614905,
                4.0563332013675515,
                4.0569392619736115,
                4.0575453225796725,
                4.058151383185733,
                4.0587574437917935,
                4.059363504397854,
                4.0599695650039145,
                4.060575625609975,
                4.061181686216036,
                4.061787746822096,
                4.062393807428157,
                4.062999868034218,
                4.063605928640278,
                4.064211989246339,
                4.064818049852399,
                4.06542411045846,
                4.066030171064521,
                4.066636231670581,
                4.067242292276642,
                4.067848352882702,
                4.068454413488763,
                4.069060474094824,
                4.069666534700884,
                4.070272595306945,
                4.070878655913006,
                4.071484716519066,
                4.072090777125127,
                4.072696837731187,
                4.073302898337248,
                4.073908958943309,
                4.074515019549369,
                4.07512108015543,
                4.07572714076149,
                4.076333201367551,
                4.076939261973612,
                4.077545322579672,
                4.078151383185733,
                4.078757443791794,
                4.079363504397854,
                4.079969565003915,
                4.080575625609975,
                4.081181686216036,
                4.081787746822097,
                4.082393807428157,
                4.082999868034218,
                4.083605928640278,
                4.084211989246339,
                4.0848180498524,
                4.08542411045846,
                4.086030171064521,
                4.086636231670582,
                4.087242292276642,
                4.087848352882703,
                4.088454413488763,
                4.089060474094824,
            ],
        ),
    ),
)
def test_extended_plane_strain(extendedplanestrainparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for extended_plane_strain.

    This test was generated using data from tests/regression/scenarios/FNSF/IN.DAT.

    :param extendedplanestrainparam: the data used to mock and assert in this test.
    :type extendedplanestrainparam: extendedplanestrainparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    (
        rradius,
        sigr,
        sigt,
        sigz,
        str_r,
        str_t,
        str_z,
        r_deflect,
    ) = sctf.extended_plane_strain(
        n_radial_array=extendedplanestrainparam.n_radial_array,
        nlayers=extendedplanestrainparam.nlayers,
        i_tf_bucking=extendedplanestrainparam.i_tf_bucking,
        nu_t=extendedplanestrainparam.nu_t,
        nu_zt=extendedplanestrainparam.nu_zt,
        ey_t=extendedplanestrainparam.ey_t,
        ey_z=extendedplanestrainparam.ey_z,
        d_curr=extendedplanestrainparam.d_curr,
        rad=extendedplanestrainparam.rad,
        v_force=extendedplanestrainparam.v_force,
    )

    # assert sigr == pytest.approx(extendedplanestrainparam.expected_sigr, rel=0.01)

    numpy.testing.assert_array_almost_equal(
        sigr, numpy.array(extendedplanestrainparam.expected_sigr), decimal=3
    )

    numpy.testing.assert_array_almost_equal(
        sigt, numpy.array(extendedplanestrainparam.expected_sigt), decimal=3
    )

    numpy.testing.assert_array_almost_equal(
        sigz, numpy.array(extendedplanestrainparam.expected_sigz), decimal=3
    )

    numpy.testing.assert_array_almost_equal(
        str_r, numpy.array(extendedplanestrainparam.expected_str_r)
    )

    numpy.testing.assert_array_almost_equal(
        str_t, numpy.array(extendedplanestrainparam.expected_str_t)
    )

    numpy.testing.assert_array_almost_equal(
        str_z, numpy.array(extendedplanestrainparam.expected_str_z)
    )

    numpy.testing.assert_array_almost_equal(
        r_deflect, numpy.array(extendedplanestrainparam.expected_r_deflect)
    )

    numpy.testing.assert_array_almost_equal(
        rradius, numpy.array(extendedplanestrainparam.expected_rradius)
    )


class EyoungParallelParam(NamedTuple):
    eyoung_j_1: Any = None

    eyoung_j_2: Any = None

    a_1: Any = None

    a_2: Any = None

    poisson_j_perp_1: Any = None

    poisson_j_perp_2: Any = None

    expected_eyoung_j_3: Any = None

    expected_a_3: Any = None

    expected_poisson_j_perp_3: Any = None


@pytest.mark.parametrize(
    "eyoungparallelparam",
    (
        EyoungParallelParam(
            eyoung_j_1=0,
            eyoung_j_2=0,
            a_1=0.010000000000000002,
            a_2=0,
            poisson_j_perp_1=0.30000001192092896,
            poisson_j_perp_2=0,
            expected_eyoung_j_3=0,
            expected_a_3=0.010000000000000002,
            expected_poisson_j_perp_3=0.30000001192092896,
        ),
        EyoungParallelParam(
            eyoung_j_1=0,
            eyoung_j_2=0,
            a_1=0.020661087836601012,
            a_2=0.010000000000000002,
            poisson_j_perp_1=0.30000001192092896,
            poisson_j_perp_2=0.30000001192092896,
            expected_eyoung_j_3=0,
            expected_a_3=0.030661087836601014,
            expected_poisson_j_perp_3=0.30000001192092896,
        ),
    ),
)
def test_eyoung_parallel(eyoungparallelparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for eyoung_parallel.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param eyoungparallelparam: the data used to mock and assert in this test.
    :type eyoungparallelparam: eyoungparallelparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    eyoung_j_3, a_3, poisson_j_perp_3 = sctf.eyoung_parallel(
        eyoung_j_1=eyoungparallelparam.eyoung_j_1,
        eyoung_j_2=eyoungparallelparam.eyoung_j_2,
        a_1=eyoungparallelparam.a_1,
        a_2=eyoungparallelparam.a_2,
        poisson_j_perp_1=eyoungparallelparam.poisson_j_perp_1,
        poisson_j_perp_2=eyoungparallelparam.poisson_j_perp_2,
    )

    assert eyoung_j_3 == pytest.approx(eyoungparallelparam.expected_eyoung_j_3)

    assert a_3 == pytest.approx(eyoungparallelparam.expected_a_3)

    assert poisson_j_perp_3 == pytest.approx(
        eyoungparallelparam.expected_poisson_j_perp_3
    )


class EyoungTNestedSquaresParam(NamedTuple):
    n: Any = None

    eyoung_j_in: Any = None

    l_in: Any = None

    poisson_j_perp_in: Any = None

    expected_eyoung_j_out: Any = None

    expected_l_out: Any = None

    expected_poisson_j_perp_out: Any = None

    expected_eyoung_stiffest: Any = None


class EyoungSeriesParam(NamedTuple):
    eyoung_j_1: Any = None

    eyoung_j_2: Any = None

    l_1: Any = None

    l_2: Any = None

    poisson_j_perp_1: Any = None

    poisson_j_perp_2: Any = None

    expected_eyoung_j_3: Any = None

    expected_l_3: Any = None

    expected_poisson_j_perp_3: Any = None


class EyoungParallelArrayParam(NamedTuple):
    n: Any = None

    eyoung_j_in: Any = None

    a_in: Any = None

    poisson_j_perp_in: Any = None

    expected_eyoung_j_out: Any = None

    expected_a_out: Any = None

    expected_poisson_j_perp_out: Any = None


@pytest.mark.parametrize(
    "eyoungtnestedsquaresparam",
    (
        EyoungTNestedSquaresParam(
            n=4,
            eyoung_j_in=numpy.array(
                numpy.array((0, 0, 205000000000, 20000000000, 0), order="F"), order="F"
            ).transpose(),
            l_in=numpy.array(
                numpy.array(
                    (
                        0.010000000000000002,
                        0.020661087836601012,
                        0.016,
                        0.0041799999999999997,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            poisson_j_perp_in=numpy.array(
                numpy.array(
                    (
                        0.29999999999999999,
                        0.30000001192092896,
                        0.29999999999999999,
                        0.34000000000000002,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_eyoung_j_out=38289891367.115105,
            expected_l_out=0.050841087836601018,
            expected_poisson_j_perp_out=0.30931445806415137,
            expected_eyoung_stiffest=116443733140.5881,
        ),
        EyoungTNestedSquaresParam(
            n=4,
            eyoung_j_in=numpy.array(
                numpy.array((0, 0, 205000000000, 20000000000, 0), order="F"), order="F"
            ).transpose(),
            l_in=numpy.array(
                numpy.array(
                    (
                        0.010000000000000002,
                        0.020661087836601012,
                        0.016,
                        0.0041799999999999997,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            poisson_j_perp_in=numpy.array(
                numpy.array(
                    (
                        0.29999999999999999,
                        0.30000001192092896,
                        0.29999999999999999,
                        0.34000000000000002,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_eyoung_j_out=38289891367.115105,
            expected_l_out=0.050841087836601018,
            expected_poisson_j_perp_out=0.30931445806415137,
            expected_eyoung_stiffest=116443733140.5881,
        ),
    ),
)
def test_eyoung_t_nested_squares(eyoungtnestedsquaresparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for eyoung_t_nested_squares.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param eyoungtnestedsquaresparam: the data used to mock and assert in this test.
    :type eyoungtnestedsquaresparam: eyoungtnestedsquaresparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    (
        eyoung_j_out,
        l_out,
        poisson_j_perp_out,
        eyoung_stiffest,
    ) = sctf.eyoung_t_nested_squares(
        n=eyoungtnestedsquaresparam.n,
        eyoung_j_in=eyoungtnestedsquaresparam.eyoung_j_in,
        l_in=eyoungtnestedsquaresparam.l_in,
        poisson_j_perp_in=eyoungtnestedsquaresparam.poisson_j_perp_in,
    )

    # assert eyoung_j_out == pytest.approx(
    #     eyoungtnestedsquaresparam.expected_eyoung_j_out
    # )

    # assert l_out == pytest.approx(eyoungtnestedsquaresparam.expected_l_out)

    # assert poisson_j_perp_out == pytest.approx(
    #     eyoungtnestedsquaresparam.expected_poisson_j_perp_out
    # )

    assert eyoung_stiffest == pytest.approx(
        eyoungtnestedsquaresparam.expected_eyoung_stiffest
    )


@pytest.mark.parametrize(
    "eyoungseriesparam",
    (
        EyoungSeriesParam(
            eyoung_j_1=0,
            eyoung_j_2=117000000000,
            l_1=0.003949573550844649,
            l_2=0.016711514285756363,
            poisson_j_perp_1=0.30000001192092896,
            poisson_j_perp_2=0.34999999999999998,
            expected_eyoung_j_3=0,
            expected_l_3=0.020661087836601012,
            expected_poisson_j_perp_3=0.30000001192092896,
        ),
        EyoungSeriesParam(
            eyoung_j_1=0,
            eyoung_j_2=0,
            l_1=0.020661087836601012,
            l_2=0.010000000000000002,
            poisson_j_perp_1=0.30000001192092896,
            poisson_j_perp_2=0.29999999999999999,
            expected_eyoung_j_3=0,
            expected_l_3=0.030661087836601014,
            expected_poisson_j_perp_3=0.30000001192092896,
        ),
    ),
)
def test_eyoung_series(eyoungseriesparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for eyoung_series.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param eyoungseriesparam: the data used to mock and assert in this test.
    :type eyoungseriesparam: eyoungseriesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    eyoung_j_3, l_3, poisson_j_perp_3 = sctf.eyoung_series(
        eyoung_j_1=eyoungseriesparam.eyoung_j_1,
        eyoung_j_2=eyoungseriesparam.eyoung_j_2,
        l_1=eyoungseriesparam.l_1,
        l_2=eyoungseriesparam.l_2,
        poisson_j_perp_1=eyoungseriesparam.poisson_j_perp_1,
        poisson_j_perp_2=eyoungseriesparam.poisson_j_perp_2,
    )

    assert eyoung_j_3 == pytest.approx(eyoungseriesparam.expected_eyoung_j_3)

    assert l_3 == pytest.approx(eyoungseriesparam.expected_l_3)

    assert poisson_j_perp_3 == pytest.approx(
        eyoungseriesparam.expected_poisson_j_perp_3
    )


@pytest.mark.parametrize(
    "eyoungparallelarrayparam",
    (
        EyoungParallelArrayParam(
            n=5,
            eyoung_j_in=numpy.array(
                numpy.array((205000000000, 20000000000, 117000000000, 0, 0), order="F"),
                order="F",
            ).transpose(),
            a_in=numpy.array(
                numpy.array(
                    (
                        0.29370123076207649,
                        0.11646247019991701,
                        0.13374756938078641,
                        0.031609694578447076,
                        0.1297552160314831,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            poisson_j_perp_in=numpy.array(
                numpy.array(
                    (
                        0.29999999999999999,
                        0.34000000000000002,
                        0.34999999999999998,
                        0.30000001192092896,
                        0.29999999999999999,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_eyoung_j_out=110859361820.72557,
            expected_a_out=0.70527618095271016,
            expected_poisson_j_perp_out=0.31608714140682664,
        ),
        EyoungParallelArrayParam(
            n=5,
            eyoung_j_in=numpy.array(
                numpy.array((205000000000, 20000000000, 117000000000, 0, 0), order="F"),
                order="F",
            ).transpose(),
            a_in=numpy.array(
                numpy.array(
                    (
                        0.29370123076207649,
                        0.11646247019991701,
                        0.13374756938078641,
                        0.031609694578447076,
                        0.1297552160314831,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            poisson_j_perp_in=numpy.array(
                numpy.array(
                    (
                        0.29999999999999999,
                        0.34000000000000002,
                        0.34999999999999998,
                        0.30000001192092896,
                        0.29999999999999999,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_eyoung_j_out=110859361820.72557,
            expected_a_out=0.70527618095271016,
            expected_poisson_j_perp_out=0.31608714140682664,
        ),
    ),
)
def test_eyoung_parallel_array(eyoungparallelarrayparam, monkeypatch, sctfcoil):
    """
    Automatically generated Regression Unit Test for eyoung_parallel_array.

    This test was generated using data from tracking/baseline_2018/baseline_2018_IN.DAT.

    :param eyoungparallelarrayparam: the data used to mock and assert in this test.
    :type eyoungparallelarrayparam: eyoungparallelarrayparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    eyoung_j_out, a_out, poisson_j_perp_out = sctf.eyoung_parallel_array(
        n=eyoungparallelarrayparam.n,
        eyoung_j_in=eyoungparallelarrayparam.eyoung_j_in,
        a_in=eyoungparallelarrayparam.a_in,
        poisson_j_perp_in=eyoungparallelarrayparam.poisson_j_perp_in,
    )

    assert eyoung_j_out == pytest.approx(eyoungparallelarrayparam.expected_eyoung_j_out)

    assert a_out == pytest.approx(eyoungparallelarrayparam.expected_a_out)

    assert poisson_j_perp_out == pytest.approx(
        eyoungparallelarrayparam.expected_poisson_j_perp_out
    )


@pytest.mark.parametrize(
    "sx, sy, sz, expected",
    (
        (0, -3.2e8, 2.4e8, 486621002.42385757),
        (-2.8e8, 0, 2.4e8, 450777106.7833858),
    ),
)
def test_sigvm(sx, sy, sz, expected, sctfcoil):
    # could not find an example of a use in PROCESS where
    # tx, ty, or tz were anything other than 0
    ret = sctf.sigvm(sx, sy, sz, 0, 0, 0)

    assert ret == pytest.approx(expected)


def test_vv_stress_on_quench():
    """Tests the VV stress on TF quench model presented in Itoh et al using the
    values they use to test the model for JA DEMO concept.
    """
    assert (
        pytest.approx(
            sctf.vv_stress_on_quench(
                # TF shape
                H_coil=9.5,
                Ri_coil=3.55,
                Ro_coil=15.62,
                Rm_coil=7.66,
                ccl_length_coil=51.1,
                theta1_coil=48,
                # VV shape
                H_vv=7.9,
                Ri_vv=4.45,
                Ro_vv=13.09,
                Rm_vv=7.88,
                theta1_vv=1,
                # TF properties
                n_tf=18,
                n_tf_turn=192,
                S_rp=0.55,
                S_cc=0.94,
                taud=30,
                I_op=83200,
                # VV properties
                d_vv=0.12,  # for 6.6 restistance -> lambda2 = 2.1
            )
        )
        == 56835032.21809308
    )


def test_vv_stress_on_quench_integration(sctfcoil, monkeypatch):
    """Tests the VV stress on TF quench model presented in Itoh et al using the
    values they use to test the model for JA DEMO concept. Includes the assumptions
    and approximations in the models integration with PROCESS.
    """
    monkeypatch.setattr(build_variables, "tfcth", 1.4)  # Baseline 2018 value
    monkeypatch.setattr(build_variables, "hmax", 8.8)  # Table 2
    monkeypatch.setattr(build_variables, "r_tf_inboard_mid", 3.55)  # Table 2
    monkeypatch.setattr(build_variables, "r_tf_outboard_mid", 15.62)  # Table 2
    monkeypatch.setattr(tfcoil_variables, "theta1_coil", 48)  # Table 2
    monkeypatch.setattr(tfcoil_variables, "theta1_vv", 1)  # Table 2
    monkeypatch.setattr(
        build_variables,
        "r_tf_inboard_out",
        build_variables.r_tf_inboard_mid + (build_variables.tfcth / 2),
    )
    monkeypatch.setattr(build_variables, "tfthko", 0)  # simplifies

    monkeypatch.setattr(physics_variables, "rminor", 2.96)  # Baseline 2018
    monkeypatch.setattr(physics_variables, "kappa", 1.848)  # Baseline 2018

    monkeypatch.setattr(sctfcoil_module, "a_tf_steel", 0.55)  # Section 3

    # Sum from Section 3
    monkeypatch.setattr(sctfcoil_module, "a_case_front", 0.47)
    monkeypatch.setattr(sctfcoil_module, "a_case_nose", 0.47)

    monkeypatch.setattr(build_variables, "vgap", 0.05)  # Baseline 2018
    monkeypatch.setattr(build_variables, "shldtth", 0.3)  # Baseline 2018
    monkeypatch.setattr(
        divertor_variables, "divfix", 2.05
    )  # chosen to achieve H_vv in Table 2

    monkeypatch.setattr(build_variables, "tftsgap", 0.05)  # Baseline 2018
    monkeypatch.setattr(build_variables, "thshield_ob", 0.05)  # Baseline 2018
    monkeypatch.setattr(build_variables, "tfthko", 1.4)  # Baseline 2018
    monkeypatch.setattr(
        build_variables, "gapsto", 1.7
    )  # chosen to achieve Ro_vv in Table 2

    monkeypatch.setattr(build_variables, "d_vv_out", 0.06)  # Section 3
    monkeypatch.setattr(build_variables, "d_vv_in", 0.06)  # Section 3
    monkeypatch.setattr(build_variables, "d_vv_top", 0.06)  # Section 3

    monkeypatch.setattr(tfcoil_variables, "tfleng", 51.1)  # Table 2
    monkeypatch.setattr(
        tfcoil_variables, "tfa", [3.41, 7.77, 7.77, 3.41]
    )  # chosen to achieve Rm_coil in Table 2
    monkeypatch.setattr(tfcoil_variables, "n_tf", 18)  # Section 3
    monkeypatch.setattr(tfcoil_variables, "n_tf_turn", 192)  # Section 3
    monkeypatch.setattr(tfcoil_variables, "tdmptf", 30)  # Figure 6
    monkeypatch.setattr(sctfcoil_module, "tfc_current", 83200 * 192)  # Section 3

    monkeypatch.setattr(
        build_variables, "r_vv_inboard_out", 4.45 + (build_variables.d_vv_in / 2)
    )  # Table 2

    sctfcoil.vv_stress_on_quench()

    assert pytest.approx(sctfcoil_module.vv_stress_quench) == 56834395.24352395
