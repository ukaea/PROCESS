import pytest
from typing import NamedTuple, Any
import numpy

from process.fortran import (
    physics_variables,
    stellarator_configuration,
    stellarator_module,
    build_variables,
    fwbs_variables,
    heat_transport_variables,
    structure_variables,
    tfcoil_variables,
    impurity_radiation_module,
)
from process.power import Power
from process.stellarator import Stellarator
from process.vacuum import Vacuum
from process.availability import Availability
from process.buildings import Buildings
from process.costs import Costs
from process.plasma_profiles import PlasmaProfile
from process.hcpb import CCFE_HCPB
from process.blanket_library import BlanketLibrary
from process.fw import Fw
from process.current_drive import CurrentDrive


@pytest.fixture
def stellarator():
    """Provides Stellarator object for testing.

    :returns: initialised Stellarator object
    :rtype: process.stellerator.Stellarator
    """
    return Stellarator(
        Availability(),
        Vacuum(),
        Buildings(),
        Costs(),
        Power(),
        PlasmaProfile(),
        CCFE_HCPB(BlanketLibrary(Fw())),
        CurrentDrive(),
    )


class StgeomParam(NamedTuple):
    aspect: Any = None

    rmajor: Any = None

    rminor: Any = None

    sarea: Any = None

    sareao: Any = None

    vol: Any = None

    xarea: Any = None

    bt: Any = None

    stella_config_plasma_volume: Any = None

    stella_config_plasma_surface: Any = None

    f_r: Any = None

    f_a: Any = None

    expected_sarea: Any = None

    expected_sareao: Any = None

    expected_vol: Any = None

    expected_xarea: Any = None


@pytest.mark.parametrize(
    "stgeomparam",
    (
        StgeomParam(
            aspect=12.33,
            rmajor=22,
            rminor=1.7842660178426601,
            sarea=0,
            sareao=0,
            vol=0,
            xarea=0,
            bt=5.5,
            stella_config_plasma_volume=1422.6300000000001,
            stella_config_plasma_surface=1960,
            f_r=0.99099099099099097,
            f_a=0.99125889880147788,
            expected_sarea=1925.3641313657533,
            expected_sareao=962.68206568287667,
            expected_vol=1385.2745877380669,
            expected_xarea=10.001590778710231,
        ),
        StgeomParam(
            aspect=12.33,
            rmajor=22,
            rminor=1.7842660178426601,
            sarea=1925.3641313657533,
            sareao=962.68206568287667,
            vol=1385.2745877380669,
            xarea=10.001590778710231,
            bt=5.5,
            stella_config_plasma_volume=1422.6300000000001,
            stella_config_plasma_surface=1960,
            f_r=0.99099099099099097,
            f_a=0.99125889880147788,
            expected_sarea=1925.3641313657533,
            expected_sareao=962.68206568287667,
            expected_vol=1385.2745877380669,
            expected_xarea=10.001590778710231,
        ),
    ),
)
def test_stgeom(stgeomparam, monkeypatch, stellarator):
    """
    Automatically generated Regression Unit Test for stgeom.

    This test was generated using data from tests/regression/scenarios/stellarator/IN.DAT.

    :param stgeomparam: the data used to mock and assert in this test.
    :type stgeomparam: stgeomparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "aspect", stgeomparam.aspect)

    monkeypatch.setattr(physics_variables, "rmajor", stgeomparam.rmajor)

    monkeypatch.setattr(physics_variables, "rminor", stgeomparam.rminor)

    monkeypatch.setattr(physics_variables, "sarea", stgeomparam.sarea)

    monkeypatch.setattr(physics_variables, "sareao", stgeomparam.sareao)

    monkeypatch.setattr(physics_variables, "vol", stgeomparam.vol)

    monkeypatch.setattr(physics_variables, "xarea", stgeomparam.xarea)

    monkeypatch.setattr(physics_variables, "bt", stgeomparam.bt)

    monkeypatch.setattr(
        stellarator_configuration,
        "stella_config_plasma_volume",
        stgeomparam.stella_config_plasma_volume,
    )

    monkeypatch.setattr(
        stellarator_configuration,
        "stella_config_plasma_surface",
        stgeomparam.stella_config_plasma_surface,
    )

    monkeypatch.setattr(stellarator_module, "f_r", stgeomparam.f_r)

    monkeypatch.setattr(stellarator_module, "f_a", stgeomparam.f_a)

    stellarator.stgeom()

    assert physics_variables.sarea == pytest.approx(stgeomparam.expected_sarea)

    assert physics_variables.sareao == pytest.approx(stgeomparam.expected_sareao)

    assert physics_variables.vol == pytest.approx(stgeomparam.expected_vol)

    assert physics_variables.xarea == pytest.approx(stgeomparam.expected_xarea)


class StbildParam(NamedTuple):
    blbmith: Any = None

    blbmoth: Any = None

    blbpith: Any = None

    blbpoth: Any = None

    blbuith: Any = None

    blbuoth: Any = None

    blnkith: Any = None

    blnkoth: Any = None

    blnktth: Any = None

    bore: Any = None

    d_vv_in: Any = None

    d_vv_out: Any = None

    fwarea: Any = None

    fwith: Any = None

    fwoth: Any = None

    gapds: Any = None

    gapoh: Any = None

    gapomin: Any = None

    gapsto: Any = None

    hmax: Any = None

    ohcth: Any = None

    r_tf_outboard_mid: Any = None

    rbld: Any = None

    rsldi: Any = None

    rsldo: Any = None

    rspo: Any = None

    scrapli: Any = None

    scraplo: Any = None

    shldith: Any = None

    shldoth: Any = None

    shldtth: Any = None

    tfcth: Any = None

    tfthko: Any = None

    available_radial_space: Any = None

    f_avspace: Any = None

    required_radial_space: Any = None

    afw: Any = None

    blktmodel: Any = None

    fdiv: Any = None

    fhcd: Any = None

    fhole: Any = None

    fw_wall: Any = None

    ipowerflow: Any = None

    rmajor: Any = None

    rminor: Any = None

    sarea: Any = None

    stella_config_derivative_min_lcfs_coils_dist: Any = None

    stella_config_rminor_ref: Any = None

    stella_config_min_plasma_coil_distance: Any = None

    f_r: Any = None

    f_aspect: Any = None

    f_a: Any = None

    iprint: Any = None

    outfile: Any = None

    expected_blnktth: Any = None

    expected_bore: Any = None

    expected_fwarea: Any = None

    expected_fwith: Any = None

    expected_fwoth: Any = None

    expected_gapsto: Any = None

    expected_hmax: Any = None

    expected_r_tf_outboard_mid: Any = None

    expected_rbld: Any = None

    expected_rsldi: Any = None

    expected_rsldo: Any = None

    expected_rspo: Any = None

    expected_available_radial_space: Any = None

    expected_required_radial_space: Any = None


@pytest.mark.parametrize(
    "stbildparam",
    (
        StbildParam(
            blbmith=0.17000000000000001,
            blbmoth=0.27000000000000002,
            blbpith=0.29999999999999999,
            blbpoth=0.34999999999999998,
            blbuith=0.36499999999999999,
            blbuoth=0.46500000000000002,
            blnkith=0.70000000000000007,
            blnkoth=0.80000000000000004,
            blnktth=0,
            bore=1.4199999999999999,
            d_vv_in=0.35000000000000003,
            d_vv_out=0.35000000000000003,
            fwarea=0,
            fwith=0,
            fwoth=0,
            gapds=0.025000000000000005,
            gapoh=0,
            gapomin=0.025000000000000005,
            gapsto=0,
            hmax=6.2927927927927927,
            ohcth=0,
            r_tf_outboard_mid=0,
            rbld=0,
            rsldi=0,
            rsldo=0,
            rspo=0,
            scrapli=0.15000000000000002,
            scraplo=0.30000000000000004,
            shldith=0.40000000000000002,
            shldoth=0.70000000000000007,
            shldtth=0.70000000000000007,
            tfcth=0.78058448071757114,
            tfthko=0.78058448071757114,
            available_radial_space=0,
            f_avspace=1,
            required_radial_space=0,
            afw=0.0060000000000000001,
            blktmodel=0,
            fdiv=0.115,
            fhcd=0,
            fhole=0,
            fw_wall=0.0030000000000000001,
            ipowerflow=1,
            rmajor=22,
            rminor=1.7842660178426601,
            sarea=1925.3641313657533,
            stella_config_derivative_min_lcfs_coils_dist=-1,
            stella_config_rminor_ref=1.8,
            stella_config_min_plasma_coil_distance=1.8999999999999999,
            f_r=0.99099099099099097,
            f_aspect=1,
            f_a=0.99125889880147788,
            iprint=0,
            outfile=11,
            expected_blnktth=0.75,
            expected_bore=17.79214950143977,
            expected_fwarea=1918.8188778803135,
            expected_fwith=0.018000000000000002,
            expected_fwoth=0.018000000000000002,
            expected_gapsto=0.025000000000000005,
            expected_hmax=3.7022660178426601,
            expected_r_tf_outboard_mid=26.367558258201448,
            expected_rbld=22,
            expected_rsldi=18.947733982157342,
            expected_rsldo=25.602266017842663,
            expected_rspo=22,
            expected_available_radial_space=1.8828828828828827,
            expected_required_radial_space=2.0332922403587861,
        ),
        StbildParam(
            blbmith=0.17000000000000001,
            blbmoth=0.27000000000000002,
            blbpith=0.29999999999999999,
            blbpoth=0.34999999999999998,
            blbuith=0.36499999999999999,
            blbuoth=0.46500000000000002,
            blnkith=0.70000000000000007,
            blnkoth=0.80000000000000004,
            blnktth=0.75,
            bore=17.79214950143977,
            d_vv_in=0.35000000000000003,
            d_vv_out=0.35000000000000003,
            fwarea=1918.8188778803135,
            fwith=0.018000000000000002,
            fwoth=0.018000000000000002,
            gapds=0.025000000000000005,
            gapoh=0,
            gapomin=0.025000000000000005,
            gapsto=0.025000000000000005,
            hmax=6.2927927927927927,
            ohcth=0,
            r_tf_outboard_mid=26.367558258201448,
            rbld=22,
            rsldi=18.947733982157342,
            rsldo=25.602266017842663,
            rspo=22,
            scrapli=0.15000000000000002,
            scraplo=0.30000000000000004,
            shldith=0.40000000000000002,
            shldoth=0.70000000000000007,
            shldtth=0.70000000000000007,
            tfcth=0.78058448071757114,
            tfthko=0.78058448071757114,
            available_radial_space=1.8828828828828827,
            f_avspace=1,
            required_radial_space=2.0332922403587861,
            afw=0.0060000000000000001,
            blktmodel=0,
            fdiv=0.021924555536480182,
            fhcd=0,
            fhole=0,
            fw_wall=0.0030000000000000001,
            ipowerflow=1,
            rmajor=22,
            rminor=1.7842660178426601,
            sarea=1925.3641313657533,
            stella_config_derivative_min_lcfs_coils_dist=-1,
            stella_config_rminor_ref=1.8,
            stella_config_min_plasma_coil_distance=1.8999999999999999,
            f_r=0.99099099099099097,
            f_aspect=1,
            f_a=0.99125889880147788,
            iprint=0,
            outfile=11,
            expected_blnktth=0.75,
            expected_bore=17.79214950143977,
            expected_fwarea=2120.6210472630282,
            expected_fwith=0.018000000000000002,
            expected_fwoth=0.018000000000000002,
            expected_gapsto=0.025000000000000005,
            expected_hmax=3.7022660178426601,
            expected_r_tf_outboard_mid=26.367558258201448,
            expected_rbld=22,
            expected_rsldi=18.947733982157342,
            expected_rsldo=25.602266017842663,
            expected_rspo=22,
            expected_available_radial_space=1.8828828828828827,
            expected_required_radial_space=2.0332922403587861,
        ),
    ),
)
def test_stbild(stbildparam, monkeypatch, stellarator):
    """
    Automatically generated Regression Unit Test for stbild.

    This test was generated using data from tests/regression/scenarios/stellarator/IN.DAT.

    :param stbildparam: the data used to mock and assert in this test.
    :type stbildparam: stbildparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(build_variables, "blbmith", stbildparam.blbmith)

    monkeypatch.setattr(build_variables, "blbmoth", stbildparam.blbmoth)

    monkeypatch.setattr(build_variables, "blbpith", stbildparam.blbpith)

    monkeypatch.setattr(build_variables, "blbpoth", stbildparam.blbpoth)

    monkeypatch.setattr(build_variables, "blbuith", stbildparam.blbuith)

    monkeypatch.setattr(build_variables, "blbuoth", stbildparam.blbuoth)

    monkeypatch.setattr(build_variables, "blnkith", stbildparam.blnkith)

    monkeypatch.setattr(build_variables, "blnkoth", stbildparam.blnkoth)

    monkeypatch.setattr(build_variables, "blnktth", stbildparam.blnktth)

    monkeypatch.setattr(build_variables, "bore", stbildparam.bore)

    monkeypatch.setattr(build_variables, "d_vv_in", stbildparam.d_vv_in)

    monkeypatch.setattr(build_variables, "d_vv_out", stbildparam.d_vv_out)

    monkeypatch.setattr(build_variables, "fwarea", stbildparam.fwarea)

    monkeypatch.setattr(build_variables, "fwith", stbildparam.fwith)

    monkeypatch.setattr(build_variables, "fwoth", stbildparam.fwoth)

    monkeypatch.setattr(build_variables, "gapds", stbildparam.gapds)

    monkeypatch.setattr(build_variables, "gapoh", stbildparam.gapoh)

    monkeypatch.setattr(build_variables, "gapomin", stbildparam.gapomin)

    monkeypatch.setattr(build_variables, "gapsto", stbildparam.gapsto)

    monkeypatch.setattr(build_variables, "hmax", stbildparam.hmax)

    monkeypatch.setattr(build_variables, "ohcth", stbildparam.ohcth)

    monkeypatch.setattr(
        build_variables, "r_tf_outboard_mid", stbildparam.r_tf_outboard_mid
    )

    monkeypatch.setattr(build_variables, "rbld", stbildparam.rbld)

    monkeypatch.setattr(build_variables, "rsldi", stbildparam.rsldi)

    monkeypatch.setattr(build_variables, "rsldo", stbildparam.rsldo)

    monkeypatch.setattr(build_variables, "rspo", stbildparam.rspo)

    monkeypatch.setattr(build_variables, "scrapli", stbildparam.scrapli)

    monkeypatch.setattr(build_variables, "scraplo", stbildparam.scraplo)

    monkeypatch.setattr(build_variables, "shldith", stbildparam.shldith)

    monkeypatch.setattr(build_variables, "shldoth", stbildparam.shldoth)

    monkeypatch.setattr(build_variables, "shldtth", stbildparam.shldtth)

    monkeypatch.setattr(build_variables, "tfcth", stbildparam.tfcth)

    monkeypatch.setattr(build_variables, "tfthko", stbildparam.tfthko)

    monkeypatch.setattr(
        build_variables, "available_radial_space", stbildparam.available_radial_space
    )

    monkeypatch.setattr(build_variables, "f_avspace", stbildparam.f_avspace)

    monkeypatch.setattr(
        build_variables, "required_radial_space", stbildparam.required_radial_space
    )

    monkeypatch.setattr(fwbs_variables, "afw", stbildparam.afw)

    monkeypatch.setattr(fwbs_variables, "blktmodel", stbildparam.blktmodel)

    monkeypatch.setattr(fwbs_variables, "fdiv", stbildparam.fdiv)

    monkeypatch.setattr(fwbs_variables, "fhcd", stbildparam.fhcd)

    monkeypatch.setattr(fwbs_variables, "fhole", stbildparam.fhole)

    monkeypatch.setattr(fwbs_variables, "fw_wall", stbildparam.fw_wall)

    monkeypatch.setattr(heat_transport_variables, "ipowerflow", stbildparam.ipowerflow)

    monkeypatch.setattr(physics_variables, "rmajor", stbildparam.rmajor)

    monkeypatch.setattr(physics_variables, "rminor", stbildparam.rminor)

    monkeypatch.setattr(physics_variables, "sarea", stbildparam.sarea)

    monkeypatch.setattr(
        stellarator_configuration,
        "stella_config_derivative_min_lcfs_coils_dist",
        stbildparam.stella_config_derivative_min_lcfs_coils_dist,
    )

    monkeypatch.setattr(
        stellarator_configuration,
        "stella_config_rminor_ref",
        stbildparam.stella_config_rminor_ref,
    )

    monkeypatch.setattr(
        stellarator_configuration,
        "stella_config_min_plasma_coil_distance",
        stbildparam.stella_config_min_plasma_coil_distance,
    )

    monkeypatch.setattr(stellarator_module, "f_r", stbildparam.f_r)

    monkeypatch.setattr(stellarator_module, "f_aspect", stbildparam.f_aspect)

    monkeypatch.setattr(stellarator_module, "f_a", stbildparam.f_a)

    stellarator.stbild(False)
    assert build_variables.blnktth == pytest.approx(stbildparam.expected_blnktth)

    assert build_variables.bore == pytest.approx(stbildparam.expected_bore)

    assert build_variables.fwarea == pytest.approx(stbildparam.expected_fwarea)

    assert build_variables.fwith == pytest.approx(stbildparam.expected_fwith)

    assert build_variables.fwoth == pytest.approx(stbildparam.expected_fwoth)

    assert build_variables.gapsto == pytest.approx(stbildparam.expected_gapsto)

    assert build_variables.hmax == pytest.approx(stbildparam.expected_hmax)

    assert build_variables.r_tf_outboard_mid == pytest.approx(
        stbildparam.expected_r_tf_outboard_mid
    )

    assert build_variables.rbld == pytest.approx(stbildparam.expected_rbld)

    assert build_variables.rsldi == pytest.approx(stbildparam.expected_rsldi)

    assert build_variables.rsldo == pytest.approx(stbildparam.expected_rsldo)

    assert build_variables.rspo == pytest.approx(stbildparam.expected_rspo)

    assert build_variables.available_radial_space == pytest.approx(
        stbildparam.expected_available_radial_space
    )

    assert build_variables.required_radial_space == pytest.approx(
        stbildparam.expected_required_radial_space
    )


class StstrcParam(NamedTuple):
    dewmkg: Any = None

    denstl: Any = None

    aintmass: Any = None

    clgsmass: Any = None

    coldmass: Any = None

    fncmass: Any = None

    gsmass: Any = None

    whttf: Any = None

    tcritsc: Any = None

    estotftgj: Any = None

    vtfskv: Any = None

    tftort: Any = None

    stella_config_coilsurface: Any = None

    stella_config_coillength: Any = None

    f_n: Any = None

    f_r: Any = None

    f_b: Any = None

    iprint: Any = None

    outfile: Any = None

    expected_aintmass: Any = None

    expected_clgsmass: Any = None

    expected_coldmass: Any = None


@pytest.mark.parametrize(
    "ststrcparam",
    (
        StstrcParam(
            dewmkg=0,
            denstl=7800,
            aintmass=0,
            clgsmass=0,
            coldmass=0,
            fncmass=0,
            gsmass=0,
            whttf=5204872.8206625767,
            tcritsc=16,
            estotftgj=132.55990646265246,
            vtfskv=4.3242392290600487,
            tftort=0.67648706726464258,
            stella_config_coilsurface=4817.6999999999998,
            stella_config_coillength=1680,
            f_n=1,
            f_r=0.99099099099099097,
            f_b=0.98214285714285721,
            iprint=0,
            outfile=11,
            expected_aintmass=4882304.266547408,
            expected_clgsmass=976460.85330948164,
            expected_coldmass=10087177.087209985,
        ),
        StstrcParam(
            dewmkg=22397931.480129492,
            denstl=7800,
            aintmass=4882304.266547408,
            clgsmass=976460.85330948164,
            coldmass=10087177.087209985,
            fncmass=0,
            gsmass=0,
            whttf=5204872.8206625767,
            tcritsc=16,
            estotftgj=132.55990646265246,
            vtfskv=4.3242392290600487,
            tftort=0.67648706726464258,
            stella_config_coilsurface=4817.6999999999998,
            stella_config_coillength=1680,
            f_n=1,
            f_r=0.99099099099099097,
            f_b=0.98214285714285721,
            iprint=0,
            outfile=11,
            expected_aintmass=4882304.266547408,
            expected_clgsmass=976460.85330948164,
            expected_coldmass=32485108.567339476,
        ),
    ),
)
def test_ststrc(ststrcparam, monkeypatch, stellarator):
    """
    Automatically generated Regression Unit Test for ststrc.

    This test was generated using data from tests/regression/scenarios/stellarator/IN.DAT.

    :param ststrcparam: the data used to mock and assert in this test.
    :type ststrcparam: ststrcparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(fwbs_variables, "dewmkg", ststrcparam.dewmkg)

    monkeypatch.setattr(fwbs_variables, "denstl", ststrcparam.denstl)

    monkeypatch.setattr(structure_variables, "aintmass", ststrcparam.aintmass)

    monkeypatch.setattr(structure_variables, "clgsmass", ststrcparam.clgsmass)

    monkeypatch.setattr(structure_variables, "coldmass", ststrcparam.coldmass)

    monkeypatch.setattr(structure_variables, "fncmass", ststrcparam.fncmass)

    monkeypatch.setattr(structure_variables, "gsmass", ststrcparam.gsmass)

    monkeypatch.setattr(tfcoil_variables, "whttf", ststrcparam.whttf)

    monkeypatch.setattr(tfcoil_variables, "tcritsc", ststrcparam.tcritsc)

    monkeypatch.setattr(tfcoil_variables, "estotftgj", ststrcparam.estotftgj)

    monkeypatch.setattr(tfcoil_variables, "vtfskv", ststrcparam.vtfskv)

    monkeypatch.setattr(tfcoil_variables, "tftort", ststrcparam.tftort)

    monkeypatch.setattr(
        stellarator_configuration,
        "stella_config_coilsurface",
        ststrcparam.stella_config_coilsurface,
    )

    monkeypatch.setattr(
        stellarator_configuration,
        "stella_config_coillength",
        ststrcparam.stella_config_coillength,
    )

    monkeypatch.setattr(stellarator_module, "f_n", ststrcparam.f_n)

    monkeypatch.setattr(stellarator_module, "f_r", ststrcparam.f_r)

    monkeypatch.setattr(stellarator_module, "f_b", ststrcparam.f_b)

    stellarator.ststrc(False)

    assert structure_variables.aintmass == pytest.approx(ststrcparam.expected_aintmass)

    assert structure_variables.clgsmass == pytest.approx(ststrcparam.expected_clgsmass)

    assert structure_variables.coldmass == pytest.approx(ststrcparam.expected_coldmass)


def test_u_max_protect_v(stellarator):
    assert stellarator.u_max_protect_v(
        tfes=2651198129.2530489, tdump=10, aio=122620.32643505408
    ) == pytest.approx(4324.2392290600483)


def test_j_max_protect_am2(stellarator):
    assert stellarator.j_max_protect_am2(
        tau_quench=10,
        t_detect=0,
        fcu=0.69000000000000017,
        fcond=0.69999999999999996,
        temp=4.2000000000000002,
        acs=0.0022141440000000008,
        aturn=0.0031360000000000008,
    ) == pytest.approx(54919989.379449144)


def test_bmax_from_awp(stellarator, monkeypatch):
    monkeypatch.setattr(stellarator_configuration, "stella_config_a1", 0.688)
    monkeypatch.setattr(stellarator_configuration, "stella_config_a2", 0.025)

    assert stellarator.bmax_from_awp(
        wp_width_radial=0.11792792792792792,
        current=12.711229086229087,
        n_tf=50,
        r_coil_major=22.237837837837837,
        r_coil_minor=4.7171171171171169,
    ) == pytest.approx(39.193416982177489)


class IntersectParam(NamedTuple):
    x1: Any = None

    y1: Any = None

    x2: Any = None

    y2: Any = None

    xin: Any = None

    expected_x: Any = None


@pytest.mark.parametrize(
    "intersectparam",
    (
        IntersectParam(
            x1=numpy.array(
                numpy.array(
                    (
                        0.11792792792792792,
                        0.14103943139119018,
                        0.16415093485445242,
                        0.18726243831771466,
                        0.21037394178097696,
                        0.2334854452442392,
                        0.25659694870750144,
                        0.27970845217076368,
                        0.30281995563402597,
                        0.32593145909728821,
                        0.34904296256055045,
                        0.37215446602381275,
                        0.39526596948707493,
                        0.41837747295033723,
                        0.44148897641359952,
                        0.46460047987686171,
                        0.487711983340124,
                        0.51082348680338618,
                        0.53393499026664848,
                        0.55704649372991077,
                        0.58015799719317307,
                        0.60326950065643525,
                        0.62638100411969755,
                        0.64949250758295984,
                        0.67260401104622203,
                        0.69571551450948421,
                        0.71882701797274662,
                        0.7419385214360088,
                        0.7650500248992711,
                        0.78816152836253328,
                        0.81127303182579558,
                        0.83438453528905787,
                        0.85749603875232017,
                        0.88060754221558235,
                        0.90371904567884453,
                        0.92683054914210694,
                        0.94994205260536912,
                        0.97305355606863142,
                        0.9961650595318936,
                        1.0192765629951557,
                        1.0423880664584182,
                        1.0654995699216803,
                        1.0886110733849426,
                        1.1117225768482049,
                        1.1348340803114672,
                        1.1579455837747294,
                        1.1810570872379917,
                        1.2041685907012538,
                        1.2272800941645161,
                        1.2503915976277784,
                        1.2735031010910405,
                        1.296614604554303,
                        1.3197261080175653,
                        1.3428376114808274,
                        1.3659491149440897,
                        1.389060618407352,
                        1.4121721218706142,
                        1.4352836253338765,
                        1.4583951287971386,
                        1.4815066322604009,
                        1.5046181357236632,
                        1.5277296391869255,
                        1.5508411426501878,
                        1.5739526461134499,
                        1.5970641495767124,
                        1.6201756530399747,
                        1.6432871565032368,
                        1.6663986599664991,
                        1.6895101634297611,
                        1.7126216668930236,
                        1.7357331703562859,
                        1.758844673819548,
                        1.7819561772828103,
                        1.8050676807460724,
                        1.8281791842093349,
                        1.8512906876725972,
                        1.8744021911358593,
                        1.8975136945991216,
                        1.9206251980623836,
                        1.9437367015256461,
                        1.9668482049889084,
                        1.9899597084521705,
                        2.0130712119154328,
                        2.0361827153786947,
                        2.0592942188419574,
                        2.0824057223052197,
                        2.1055172257684815,
                        2.1286287292317438,
                        2.1517402326950061,
                        2.1748517361582684,
                        2.1979632396215307,
                        2.221074743084793,
                        2.2441862465480553,
                        2.2672977500113176,
                        2.2904092534745795,
                        2.3135207569378422,
                        2.3366322604011041,
                        2.3597437638643664,
                        2.3828552673276286,
                        2.4059667707908909,
                        2.4290782742541528,
                        2.4521897777174155,
                        2.4753012811806778,
                        2.4984127846439401,
                        2.5215242881072024,
                        2.5446357915704643,
                        2.5677472950337266,
                        2.5908587984969889,
                        2.6139703019602512,
                        2.6370818054235139,
                        2.6601933088867757,
                        2.683304812350038,
                        2.7064163158133003,
                        2.7295278192765626,
                        2.7526393227398249,
                        2.7757508262030868,
                        2.7988623296663491,
                        2.8219738331296114,
                        2.8450853365928737,
                        2.8681968400561364,
                        2.8913083435193982,
                        2.9144198469826605,
                        2.9375313504459228,
                        2.9606428539091851,
                        2.9837543573724474,
                        3.0068658608357097,
                        3.0299773642989716,
                        3.0530888677622339,
                        3.0762003712254966,
                        3.0993118746887589,
                        3.1224233781520212,
                        3.1455348816152831,
                        3.1686463850785453,
                        3.1917578885418076,
                        3.2148693920050699,
                        3.2379808954683322,
                        3.2610923989315941,
                        3.2842039023948564,
                        3.3073154058581191,
                        3.3304269093213814,
                        3.3535384127846437,
                        3.3766499162479056,
                        3.3997614197111679,
                        3.4228729231744301,
                        3.4459844266376924,
                        3.4690959301009547,
                        3.4922074335642166,
                        3.5153189370274789,
                        3.5384304404907416,
                        3.5615419439540039,
                        3.5846534474172662,
                        3.6077649508805281,
                        3.6308764543437904,
                        3.6539879578070527,
                        3.677099461270315,
                        3.7002109647335772,
                        3.7233224681968391,
                        3.7464339716601018,
                        3.7695454751233641,
                        3.7926569785866264,
                        3.8157684820498887,
                        3.8388799855131506,
                        3.8619914889764129,
                        3.8851029924396752,
                        3.9082144959029375,
                        3.9313259993661998,
                        3.9544375028294616,
                        3.9775490062927243,
                        4.0006605097559866,
                        4.0237720132192489,
                        4.0468835166825112,
                        4.0699950201457735,
                        4.0931065236090358,
                        4.1162180270722981,
                        4.1393295305355604,
                        4.1624410339988227,
                        4.185552537462085,
                        4.2086640409253473,
                        4.2317755443886096,
                        4.2548870478518719,
                        4.2779985513151342,
                        4.3011100547783965,
                        4.3242215582416588,
                        4.3473330617049211,
                        4.3704445651681834,
                        4.3935560686314457,
                        4.4166675720947079,
                        4.4397790755579694,
                        4.4628905790212317,
                        4.4860020824844939,
                        4.5091135859477571,
                        4.5322250894110194,
                        4.5553365928742808,
                        4.5784480963375431,
                        4.6015595998008054,
                        4.6246711032640677,
                        4.64778260672733,
                        4.6708941101905923,
                        4.6940056136538546,
                        4.7171171171171169,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            y1=(
                5.0000000000000004e-16,
                5.0000000000000004e-16,
                2.86227506204928e-06,
                3.1894470530728806e-06,
                3.501769913275426e-06,
                3.8002321316979702e-06,
                1.2911174666617371,
                8.5056437504875237,
                20.310566119714178,
                35.08860867240201,
                51.70918610389198,
                69.396748622632828,
                87.620394571705575,
                106.01674628541932,
                124.33769120996891,
                142.41502100154381,
                160.13633420240603,
                177.42845527181157,
                194.24590993730428,
                210.56283761600093,
                226.36726601523148,
                241.65702673605452,
                256.43682253220271,
                270.71611044750449,
                284.50756794064279,
                297.82597881823477,
                310.68742355249833,
                323.10869161811542,
                335.10685659501451,
                346.69897109328275,
                357.90185016614276,
                368.73192020856345,
                379.20511636351574,
                389.33681584460641,
                399.14179779875741,
                408.63422270274015,
                417.82762604428638,
                426.73492234723324,
                435.36841657939482,
                443.73982071745854,
                451.86027379765881,
                459.74036420016017,
                467.39015323269564,
                474.81919932015148,
                482.03658229003128,
                489.05092738310611,
                495.87042872446028,
                502.50287207059142,
                508.95565670912913,
                515.23581643372211,
                521.3500395511744,
                527.30468790362568,
                533.10581490756374,
                538.75918262525659,
                544.2702778940095,
                549.64432754547875,
                554.88631275174816,
                560.0009825376078,
                564.99286649991768,
                569.86628677533872,
                574.62536929745033,
                579.2740543834276,
                583.81610668926919,
                588.25512457110517,
                592.59454888849484,
                596.83767128391446,
                600.98764197086132,
                605.04747706123726,
                609.02006546092423,
                612.90817536075372,
                616.71446034842711,
                620.44146516535295,
                624.09163113086288,
                627.66730125482525,
                631.17072505831391,
                634.60406312072689,
                637.96939137051618,
                641.26870513558185,
                644.50392296830967,
                647.67689025924653,
                650.7893826524753,
                653.84310927488525,
                656.83971579073147,
                659.78078729211506,
                662.66785103531993,
                665.50237903228447,
                668.2857905058803,
                671.01945421709422,
                673.70469067169677,
                676.34277421347099,
                678.93493501062949,
                681.48236094161882,
                683.98619938610727,
                686.4475589265902,
                688.86751096569674,
                691.24709126396783,
                693.58730140256318,
                695.88911017509122,
                698.15345491248911,
                700.3812427446336,
                702.57335180214466,
                704.73063236162898,
                706.85390793741499,
                708.94397632264361,
                711.0016105824119,
                713.0275600015043,
                715.02255098909779,
                716.98728794267436,
                718.92245407327232,
                720.82871219404842,
                722.70670547403643,
                724.55705815886756,
                726.38037626011544,
                728.17724821484387,
                729.94824551683462,
                731.69392332090069,
                733.41482102160455,
                735.11146280763103,
                736.78435819299762,
                738.43400252621257,
                740.06087747844401,
                741.66545151169328,
                743.24818032791723,
                744.80950729999677,
                746.34986388539869,
                747.86967002332915,
                749.36933451614414,
                750.84925539573658,
                752.30982027558071,
                753.75140668908614,
                755.17438241487491,
                756.57910578956728,
                757.96592600862562,
                759.33518341579338,
                760.6872097816165,
                762.02232857153524,
                763.34085520398844,
                764.64309729896991,
                765.92935491743606,
                767.19992079195924,
                768.45508054900256,
                769.69511292315394,
                770.92028996367185,
                772.13087723365243,
                773.32713400212515,
                774.5093134293686,
                775.67766274572489,
                776.83242342417338,
                777.97383134691734,
                779.10211696622468,
                780.21750545975488,
                781.32021688058262,
                782.41046630213918,
                783.48846395826172,
                784.55441537854597,
                785.60852151918687,
                786.65097888947412,
                787.68197967412061,
                788.70171185157085,
                789.71035930845301,
                790.70810195031345,
                791.6951158087769,
                792.67157314526526,
                793.63764255140313,
                794.5934890462313,
                795.53927417035038,
                796.47515607710181,
                797.40128962089454,
                798.31782644278928,
                799.22491505342578,
                800.12270091340042,
                801.01132651117723,
                801.89093143862772,
                802.76165246427263,
                803.62362360431916,
                804.47697619156065,
                805.32183894221828,
                806.15833802079419,
                806.98659710300694,
                807.80673743687316,
                808.61887790199683,
                809.42313506713379,
                810.21962324608114,
                811.00845455195702,
                811.78973894991509,
                812.5635843083561,
                813.33009644867627,
                814.08937919361006,
                814.84153441420722,
                815.58666207549129,
                816.32486028084327,
                817.05622531514791,
                817.78085168675068,
                818.49883216825299,
                819.21025783619109,
                819.91521810963047,
                820.61380078771219,
                821.30609208618284,
                821.99217667293931,
                822.67213770261992,
            ),
            x2=numpy.array(
                numpy.array(
                    (
                        0.11792792792792792,
                        0.14103943139119018,
                        0.16415093485445242,
                        0.18726243831771466,
                        0.21037394178097696,
                        0.2334854452442392,
                        0.25659694870750144,
                        0.27970845217076368,
                        0.30281995563402597,
                        0.32593145909728821,
                        0.34904296256055045,
                        0.37215446602381275,
                        0.39526596948707493,
                        0.41837747295033723,
                        0.44148897641359952,
                        0.46460047987686171,
                        0.487711983340124,
                        0.51082348680338618,
                        0.53393499026664848,
                        0.55704649372991077,
                        0.58015799719317307,
                        0.60326950065643525,
                        0.62638100411969755,
                        0.64949250758295984,
                        0.67260401104622203,
                        0.69571551450948421,
                        0.71882701797274662,
                        0.7419385214360088,
                        0.7650500248992711,
                        0.78816152836253328,
                        0.81127303182579558,
                        0.83438453528905787,
                        0.85749603875232017,
                        0.88060754221558235,
                        0.90371904567884453,
                        0.92683054914210694,
                        0.94994205260536912,
                        0.97305355606863142,
                        0.9961650595318936,
                        1.0192765629951557,
                        1.0423880664584182,
                        1.0654995699216803,
                        1.0886110733849426,
                        1.1117225768482049,
                        1.1348340803114672,
                        1.1579455837747294,
                        1.1810570872379917,
                        1.2041685907012538,
                        1.2272800941645161,
                        1.2503915976277784,
                        1.2735031010910405,
                        1.296614604554303,
                        1.3197261080175653,
                        1.3428376114808274,
                        1.3659491149440897,
                        1.389060618407352,
                        1.4121721218706142,
                        1.4352836253338765,
                        1.4583951287971386,
                        1.4815066322604009,
                        1.5046181357236632,
                        1.5277296391869255,
                        1.5508411426501878,
                        1.5739526461134499,
                        1.5970641495767124,
                        1.6201756530399747,
                        1.6432871565032368,
                        1.6663986599664991,
                        1.6895101634297611,
                        1.7126216668930236,
                        1.7357331703562859,
                        1.758844673819548,
                        1.7819561772828103,
                        1.8050676807460724,
                        1.8281791842093349,
                        1.8512906876725972,
                        1.8744021911358593,
                        1.8975136945991216,
                        1.9206251980623836,
                        1.9437367015256461,
                        1.9668482049889084,
                        1.9899597084521705,
                        2.0130712119154328,
                        2.0361827153786947,
                        2.0592942188419574,
                        2.0824057223052197,
                        2.1055172257684815,
                        2.1286287292317438,
                        2.1517402326950061,
                        2.1748517361582684,
                        2.1979632396215307,
                        2.221074743084793,
                        2.2441862465480553,
                        2.2672977500113176,
                        2.2904092534745795,
                        2.3135207569378422,
                        2.3366322604011041,
                        2.3597437638643664,
                        2.3828552673276286,
                        2.4059667707908909,
                        2.4290782742541528,
                        2.4521897777174155,
                        2.4753012811806778,
                        2.4984127846439401,
                        2.5215242881072024,
                        2.5446357915704643,
                        2.5677472950337266,
                        2.5908587984969889,
                        2.6139703019602512,
                        2.6370818054235139,
                        2.6601933088867757,
                        2.683304812350038,
                        2.7064163158133003,
                        2.7295278192765626,
                        2.7526393227398249,
                        2.7757508262030868,
                        2.7988623296663491,
                        2.8219738331296114,
                        2.8450853365928737,
                        2.8681968400561364,
                        2.8913083435193982,
                        2.9144198469826605,
                        2.9375313504459228,
                        2.9606428539091851,
                        2.9837543573724474,
                        3.0068658608357097,
                        3.0299773642989716,
                        3.0530888677622339,
                        3.0762003712254966,
                        3.0993118746887589,
                        3.1224233781520212,
                        3.1455348816152831,
                        3.1686463850785453,
                        3.1917578885418076,
                        3.2148693920050699,
                        3.2379808954683322,
                        3.2610923989315941,
                        3.2842039023948564,
                        3.3073154058581191,
                        3.3304269093213814,
                        3.3535384127846437,
                        3.3766499162479056,
                        3.3997614197111679,
                        3.4228729231744301,
                        3.4459844266376924,
                        3.4690959301009547,
                        3.4922074335642166,
                        3.5153189370274789,
                        3.5384304404907416,
                        3.5615419439540039,
                        3.5846534474172662,
                        3.6077649508805281,
                        3.6308764543437904,
                        3.6539879578070527,
                        3.677099461270315,
                        3.7002109647335772,
                        3.7233224681968391,
                        3.7464339716601018,
                        3.7695454751233641,
                        3.7926569785866264,
                        3.8157684820498887,
                        3.8388799855131506,
                        3.8619914889764129,
                        3.8851029924396752,
                        3.9082144959029375,
                        3.9313259993661998,
                        3.9544375028294616,
                        3.9775490062927243,
                        4.0006605097559866,
                        4.0237720132192489,
                        4.0468835166825112,
                        4.0699950201457735,
                        4.0931065236090358,
                        4.1162180270722981,
                        4.1393295305355604,
                        4.1624410339988227,
                        4.185552537462085,
                        4.2086640409253473,
                        4.2317755443886096,
                        4.2548870478518719,
                        4.2779985513151342,
                        4.3011100547783965,
                        4.3242215582416588,
                        4.3473330617049211,
                        4.3704445651681834,
                        4.3935560686314457,
                        4.4166675720947079,
                        4.4397790755579694,
                        4.4628905790212317,
                        4.4860020824844939,
                        4.5091135859477571,
                        4.5322250894110194,
                        4.5553365928742808,
                        4.5784480963375431,
                        4.6015595998008054,
                        4.6246711032640677,
                        4.64778260672733,
                        4.6708941101905923,
                        4.6940056136538546,
                        4.7171171171171169,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            y2=numpy.array(
                numpy.array(
                    (
                        7158.8937047628706,
                        5004.9316715329842,
                        3694.8135594405553,
                        2839.0817737773841,
                        2249.5484991256844,
                        1826.2474529253161,
                        1512.0852402131027,
                        1272.530117074451,
                        1085.7010719257139,
                        937.18793256963431,
                        817.18705296685539,
                        718.84090024522891,
                        637.23614115501732,
                        568.77783627650172,
                        510.786630516309,
                        461.23254243400288,
                        418.55486713593308,
                        381.53776842598086,
                        349.2227154220239,
                        320.84580081746759,
                        295.79217667265311,
                        273.56246656674404,
                        253.74768704408271,
                        236.01030089701965,
                        220.06974682397441,
                        205.69127634981271,
                        192.67726151788867,
                        180.86036756672672,
                        170.09814691732512,
                        160.26872610227511,
                        151.26734021373531,
                        143.00352975000524,
                        135.39885901833182,
                        128.3850480674335,
                        121.90243465834536,
                        115.8987012784103,
                        110.32781625621425,
                        105.14914879151058,
                        100.32672600489369,
                        95.828606544860236,
                        91.626350312843897,
                        87.694567812436034,
                        84.010535746139567,
                        80.553867959068427,
                        77.306232806090279,
                        74.251109605447908,
                        71.373578121094255,
                        68.660136052082635,
                        66.098540350235496,
                        63.677668875882041,
                        61.387399466223862,
                        59.218503955912581,
                        57.162555073740741,
                        55.211844458103641,
                        53.35931029918256,
                        51.598473337326723,
                        49.923380132688088,
                        48.328552677103559,
                        46.808943550648046,
                        45.359895936370307,
                        43.977107900883965,
                        42.656600428510146,
                        41.39468876485914,
                        40.187956683990443,
                        39.033233343175098,
                        37.927572432105393,
                        36.868233360237589,
                        35.852664257722239,
                        34.878486592828409,
                        33.943481232542041,
                        33.045575793648069,
                        32.182833149542383,
                        31.353440973645949,
                        30.555702213931628,
                        29.788026404998764,
                        29.048921734577537,
                        28.336987790509955,
                        27.650908922311828,
                        26.989448158512342,
                        26.351441627222858,
                        25.73579343291075,
                        25.141470947240073,
                        24.567500476169894,
                        24.012963269341064,
                        23.476991841193676,
                        22.958766576292955,
                        22.457512594044672,
                        21.972496850393419,
                        21.503025456251216,
                        21.048441194330234,
                        20.608121217778795,
                        20.181474915566099,
                        19.767941930949156,
                        19.366990320602941,
                        18.978114843116536,
                        18.600835366568806,
                        18.234695385808415,
                        17.879260640884635,
                        17.53411782881939,
                        17.198873401581842,
                        16.873152443735915,
                        16.55659762378199,
                        16.248868213714225,
                        15.949639171768748,
                        15.658600283750749,
                        15.375455358703663,
                        15.099921475025205,
                        14.831728273446556,
                        14.570617293574582,
                        14.316341350956414,
                        14.068663951862252,
                        13.827358743198703,
                        13.592208995163054,
                        13.363007114430141,
                        13.139554185829615,
                        12.921659540623899,
                        12.709140349636851,
                        12.501821239611756,
                        12.299533931295143,
                        12.102116897851683,
                        11.909415042315253,
                        11.721279392873315,
                        11.537566814866656,
                        11.358139738464688,
                        11.182865901048743,
                        11.011618103402494,
                        10.844273978870199,
                        10.680715774700444,
                        10.520830144845792,
                        10.364507953537444,
                        10.211644088999289,
                        10.062137286707557,
                        9.91588996164114,
                        9.7728080490036398,
                        9.6328008529317479,
                        9.4957809027354614,
                        9.3616638162447625,
                        9.2303681698639899,
                        9.101815374960319,
                        8.9759295602358886,
                        8.8526374597548152,
                        8.7318683063165654,
                        8.6135537298858829,
                        8.4976276608070904,
                        8.3840262375469035,
                        8.2726877187252423,
                        8.1635523992077417,
                        8.0565625300470547,
                        7.9516622420725014,
                        7.848797472939288,
                        7.7479158974594009,
                        7.648966861046592,
                        7.5519013161173332,
                        7.456671761298697,
                        7.3632321833024221,
                        7.2715380013323854,
                        7.1815460139000722,
                        7.0932143479295533,
                        7.0065024100400572,
                        6.9213708399002529,
                        6.8377814655542144,
                        6.7556972606243528,
                        6.6750823033017541,
                        6.5959017370391155,
                        6.5181217328659393,
                        6.4417094532499517,
                        6.3666330174326555,
                        6.2928614681706652,
                        6.2203647398180992,
                        6.1491136276885232,
                        6.0790797586382022,
                        6.0102355628153186,
                        5.9425542465226266,
                        5.8760097661436923,
                        5.8105768030853397,
                        5.7462307396912662,
                        5.6829476360840934,
                        5.6207042078951268,
                        5.5594778048431799,
                        5.4992463901256654,
                        5.439988520586919,
                        5.3816833276304639,
                        5.3243104988434826,
                        5.2678502603032751,
                        5.2122833595369489,
                        5.157591049106899,
                        5.1037550707959838,
                        5.0507576403674586,
                        4.9985814328759144,
                        4.9472095685066089,
                        4.8966255989215437,
                        4.8468134940916663,
                        4.797757629595548,
                        4.7494427743657095,
                        4.7018540788646925,
                        4.6549770636737042,
                        4.6087976084774764,
                        4.5633019414297058,
                        4.5184766288841098,
                        4.4743085654767931,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            xin=0.22251193896599297,
            expected_x=0.624584480717571,
        ),
        IntersectParam(
            x1=numpy.array(
                numpy.array(
                    (
                        0.11792792792792792,
                        0.14103943139119018,
                        0.16415093485445242,
                        0.18726243831771466,
                        0.21037394178097696,
                        0.2334854452442392,
                        0.25659694870750144,
                        0.27970845217076368,
                        0.30281995563402597,
                        0.32593145909728821,
                        0.34904296256055045,
                        0.37215446602381275,
                        0.39526596948707493,
                        0.41837747295033723,
                        0.44148897641359952,
                        0.46460047987686171,
                        0.487711983340124,
                        0.51082348680338618,
                        0.53393499026664848,
                        0.55704649372991077,
                        0.58015799719317307,
                        0.60326950065643525,
                        0.62638100411969755,
                        0.64949250758295984,
                        0.67260401104622203,
                        0.69571551450948421,
                        0.71882701797274662,
                        0.7419385214360088,
                        0.7650500248992711,
                        0.78816152836253328,
                        0.81127303182579558,
                        0.83438453528905787,
                        0.85749603875232017,
                        0.88060754221558235,
                        0.90371904567884453,
                        0.92683054914210694,
                        0.94994205260536912,
                        0.97305355606863142,
                        0.9961650595318936,
                        1.0192765629951557,
                        1.0423880664584182,
                        1.0654995699216803,
                        1.0886110733849426,
                        1.1117225768482049,
                        1.1348340803114672,
                        1.1579455837747294,
                        1.1810570872379917,
                        1.2041685907012538,
                        1.2272800941645161,
                        1.2503915976277784,
                        1.2735031010910405,
                        1.296614604554303,
                        1.3197261080175653,
                        1.3428376114808274,
                        1.3659491149440897,
                        1.389060618407352,
                        1.4121721218706142,
                        1.4352836253338765,
                        1.4583951287971386,
                        1.4815066322604009,
                        1.5046181357236632,
                        1.5277296391869255,
                        1.5508411426501878,
                        1.5739526461134499,
                        1.5970641495767124,
                        1.6201756530399747,
                        1.6432871565032368,
                        1.6663986599664991,
                        1.6895101634297611,
                        1.7126216668930236,
                        1.7357331703562859,
                        1.758844673819548,
                        1.7819561772828103,
                        1.8050676807460724,
                        1.8281791842093349,
                        1.8512906876725972,
                        1.8744021911358593,
                        1.8975136945991216,
                        1.9206251980623836,
                        1.9437367015256461,
                        1.9668482049889084,
                        1.9899597084521705,
                        2.0130712119154328,
                        2.0361827153786947,
                        2.0592942188419574,
                        2.0824057223052197,
                        2.1055172257684815,
                        2.1286287292317438,
                        2.1517402326950061,
                        2.1748517361582684,
                        2.1979632396215307,
                        2.221074743084793,
                        2.2441862465480553,
                        2.2672977500113176,
                        2.2904092534745795,
                        2.3135207569378422,
                        2.3366322604011041,
                        2.3597437638643664,
                        2.3828552673276286,
                        2.4059667707908909,
                        2.4290782742541528,
                        2.4521897777174155,
                        2.4753012811806778,
                        2.4984127846439401,
                        2.5215242881072024,
                        2.5446357915704643,
                        2.5677472950337266,
                        2.5908587984969889,
                        2.6139703019602512,
                        2.6370818054235139,
                        2.6601933088867757,
                        2.683304812350038,
                        2.7064163158133003,
                        2.7295278192765626,
                        2.7526393227398249,
                        2.7757508262030868,
                        2.7988623296663491,
                        2.8219738331296114,
                        2.8450853365928737,
                        2.8681968400561364,
                        2.8913083435193982,
                        2.9144198469826605,
                        2.9375313504459228,
                        2.9606428539091851,
                        2.9837543573724474,
                        3.0068658608357097,
                        3.0299773642989716,
                        3.0530888677622339,
                        3.0762003712254966,
                        3.0993118746887589,
                        3.1224233781520212,
                        3.1455348816152831,
                        3.1686463850785453,
                        3.1917578885418076,
                        3.2148693920050699,
                        3.2379808954683322,
                        3.2610923989315941,
                        3.2842039023948564,
                        3.3073154058581191,
                        3.3304269093213814,
                        3.3535384127846437,
                        3.3766499162479056,
                        3.3997614197111679,
                        3.4228729231744301,
                        3.4459844266376924,
                        3.4690959301009547,
                        3.4922074335642166,
                        3.5153189370274789,
                        3.5384304404907416,
                        3.5615419439540039,
                        3.5846534474172662,
                        3.6077649508805281,
                        3.6308764543437904,
                        3.6539879578070527,
                        3.677099461270315,
                        3.7002109647335772,
                        3.7233224681968391,
                        3.7464339716601018,
                        3.7695454751233641,
                        3.7926569785866264,
                        3.8157684820498887,
                        3.8388799855131506,
                        3.8619914889764129,
                        3.8851029924396752,
                        3.9082144959029375,
                        3.9313259993661998,
                        3.9544375028294616,
                        3.9775490062927243,
                        4.0006605097559866,
                        4.0237720132192489,
                        4.0468835166825112,
                        4.0699950201457735,
                        4.0931065236090358,
                        4.1162180270722981,
                        4.1393295305355604,
                        4.1624410339988227,
                        4.185552537462085,
                        4.2086640409253473,
                        4.2317755443886096,
                        4.2548870478518719,
                        4.2779985513151342,
                        4.3011100547783965,
                        4.3242215582416588,
                        4.3473330617049211,
                        4.3704445651681834,
                        4.3935560686314457,
                        4.4166675720947079,
                        4.4397790755579694,
                        4.4628905790212317,
                        4.4860020824844939,
                        4.5091135859477571,
                        4.5322250894110194,
                        4.5553365928742808,
                        4.5784480963375431,
                        4.6015595998008054,
                        4.6246711032640677,
                        4.64778260672733,
                        4.6708941101905923,
                        4.6940056136538546,
                        4.7171171171171169,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            y1=(
                5.0000000000000004e-16,
                5.0000000000000004e-16,
                2.86227506204928e-06,
                3.1894470530728806e-06,
                3.501769913275426e-06,
                3.8002321316979702e-06,
                1.2911174666617371,
                8.5056437504875237,
                20.310566119714178,
                35.08860867240201,
                51.70918610389198,
                69.396748622632828,
                87.620394571705575,
                106.01674628541932,
                124.33769120996891,
                142.41502100154381,
                160.13633420240603,
                177.42845527181157,
                194.24590993730428,
                210.56283761600093,
                226.36726601523148,
                241.65702673605452,
                256.43682253220271,
                270.71611044750449,
                284.50756794064279,
                297.82597881823477,
                310.68742355249833,
                323.10869161811542,
                335.10685659501451,
                346.69897109328275,
                357.90185016614276,
                368.73192020856345,
                379.20511636351574,
                389.33681584460641,
                399.14179779875741,
                408.63422270274015,
                417.82762604428638,
                426.73492234723324,
                435.36841657939482,
                443.73982071745854,
                451.86027379765881,
                459.74036420016017,
                467.39015323269564,
                474.81919932015148,
                482.03658229003128,
                489.05092738310611,
                495.87042872446028,
                502.50287207059142,
                508.95565670912913,
                515.23581643372211,
                521.3500395511744,
                527.30468790362568,
                533.10581490756374,
                538.75918262525659,
                544.2702778940095,
                549.64432754547875,
                554.88631275174816,
                560.0009825376078,
                564.99286649991768,
                569.86628677533872,
                574.62536929745033,
                579.2740543834276,
                583.81610668926919,
                588.25512457110517,
                592.59454888849484,
                596.83767128391446,
                600.98764197086132,
                605.04747706123726,
                609.02006546092423,
                612.90817536075372,
                616.71446034842711,
                620.44146516535295,
                624.09163113086288,
                627.66730125482525,
                631.17072505831391,
                634.60406312072689,
                637.96939137051618,
                641.26870513558185,
                644.50392296830967,
                647.67689025924653,
                650.7893826524753,
                653.84310927488525,
                656.83971579073147,
                659.78078729211506,
                662.66785103531993,
                665.50237903228447,
                668.2857905058803,
                671.01945421709422,
                673.70469067169677,
                676.34277421347099,
                678.93493501062949,
                681.48236094161882,
                683.98619938610727,
                686.4475589265902,
                688.86751096569674,
                691.24709126396783,
                693.58730140256318,
                695.88911017509122,
                698.15345491248911,
                700.3812427446336,
                702.57335180214466,
                704.73063236162898,
                706.85390793741499,
                708.94397632264361,
                711.0016105824119,
                713.0275600015043,
                715.02255098909779,
                716.98728794267436,
                718.92245407327232,
                720.82871219404842,
                722.70670547403643,
                724.55705815886756,
                726.38037626011544,
                728.17724821484387,
                729.94824551683462,
                731.69392332090069,
                733.41482102160455,
                735.11146280763103,
                736.78435819299762,
                738.43400252621257,
                740.06087747844401,
                741.66545151169328,
                743.24818032791723,
                744.80950729999677,
                746.34986388539869,
                747.86967002332915,
                749.36933451614414,
                750.84925539573658,
                752.30982027558071,
                753.75140668908614,
                755.17438241487491,
                756.57910578956728,
                757.96592600862562,
                759.33518341579338,
                760.6872097816165,
                762.02232857153524,
                763.34085520398844,
                764.64309729896991,
                765.92935491743606,
                767.19992079195924,
                768.45508054900256,
                769.69511292315394,
                770.92028996367185,
                772.13087723365243,
                773.32713400212515,
                774.5093134293686,
                775.67766274572489,
                776.83242342417338,
                777.97383134691734,
                779.10211696622468,
                780.21750545975488,
                781.32021688058262,
                782.41046630213918,
                783.48846395826172,
                784.55441537854597,
                785.60852151918687,
                786.65097888947412,
                787.68197967412061,
                788.70171185157085,
                789.71035930845301,
                790.70810195031345,
                791.6951158087769,
                792.67157314526526,
                793.63764255140313,
                794.5934890462313,
                795.53927417035038,
                796.47515607710181,
                797.40128962089454,
                798.31782644278928,
                799.22491505342578,
                800.12270091340042,
                801.01132651117723,
                801.89093143862772,
                802.76165246427263,
                803.62362360431916,
                804.47697619156065,
                805.32183894221828,
                806.15833802079419,
                806.98659710300694,
                807.80673743687316,
                808.61887790199683,
                809.42313506713379,
                810.21962324608114,
                811.00845455195702,
                811.78973894991509,
                812.5635843083561,
                813.33009644867627,
                814.08937919361006,
                814.84153441420722,
                815.58666207549129,
                816.32486028084327,
                817.05622531514791,
                817.78085168675068,
                818.49883216825299,
                819.21025783619109,
                819.91521810963047,
                820.61380078771219,
                821.30609208618284,
                821.99217667293931,
                822.67213770261992,
            ),
            x2=numpy.array(
                numpy.array(
                    (
                        0.11792792792792792,
                        0.14103943139119018,
                        0.16415093485445242,
                        0.18726243831771466,
                        0.21037394178097696,
                        0.2334854452442392,
                        0.25659694870750144,
                        0.27970845217076368,
                        0.30281995563402597,
                        0.32593145909728821,
                        0.34904296256055045,
                        0.37215446602381275,
                        0.39526596948707493,
                        0.41837747295033723,
                        0.44148897641359952,
                        0.46460047987686171,
                        0.487711983340124,
                        0.51082348680338618,
                        0.53393499026664848,
                        0.55704649372991077,
                        0.58015799719317307,
                        0.60326950065643525,
                        0.62638100411969755,
                        0.64949250758295984,
                        0.67260401104622203,
                        0.69571551450948421,
                        0.71882701797274662,
                        0.7419385214360088,
                        0.7650500248992711,
                        0.78816152836253328,
                        0.81127303182579558,
                        0.83438453528905787,
                        0.85749603875232017,
                        0.88060754221558235,
                        0.90371904567884453,
                        0.92683054914210694,
                        0.94994205260536912,
                        0.97305355606863142,
                        0.9961650595318936,
                        1.0192765629951557,
                        1.0423880664584182,
                        1.0654995699216803,
                        1.0886110733849426,
                        1.1117225768482049,
                        1.1348340803114672,
                        1.1579455837747294,
                        1.1810570872379917,
                        1.2041685907012538,
                        1.2272800941645161,
                        1.2503915976277784,
                        1.2735031010910405,
                        1.296614604554303,
                        1.3197261080175653,
                        1.3428376114808274,
                        1.3659491149440897,
                        1.389060618407352,
                        1.4121721218706142,
                        1.4352836253338765,
                        1.4583951287971386,
                        1.4815066322604009,
                        1.5046181357236632,
                        1.5277296391869255,
                        1.5508411426501878,
                        1.5739526461134499,
                        1.5970641495767124,
                        1.6201756530399747,
                        1.6432871565032368,
                        1.6663986599664991,
                        1.6895101634297611,
                        1.7126216668930236,
                        1.7357331703562859,
                        1.758844673819548,
                        1.7819561772828103,
                        1.8050676807460724,
                        1.8281791842093349,
                        1.8512906876725972,
                        1.8744021911358593,
                        1.8975136945991216,
                        1.9206251980623836,
                        1.9437367015256461,
                        1.9668482049889084,
                        1.9899597084521705,
                        2.0130712119154328,
                        2.0361827153786947,
                        2.0592942188419574,
                        2.0824057223052197,
                        2.1055172257684815,
                        2.1286287292317438,
                        2.1517402326950061,
                        2.1748517361582684,
                        2.1979632396215307,
                        2.221074743084793,
                        2.2441862465480553,
                        2.2672977500113176,
                        2.2904092534745795,
                        2.3135207569378422,
                        2.3366322604011041,
                        2.3597437638643664,
                        2.3828552673276286,
                        2.4059667707908909,
                        2.4290782742541528,
                        2.4521897777174155,
                        2.4753012811806778,
                        2.4984127846439401,
                        2.5215242881072024,
                        2.5446357915704643,
                        2.5677472950337266,
                        2.5908587984969889,
                        2.6139703019602512,
                        2.6370818054235139,
                        2.6601933088867757,
                        2.683304812350038,
                        2.7064163158133003,
                        2.7295278192765626,
                        2.7526393227398249,
                        2.7757508262030868,
                        2.7988623296663491,
                        2.8219738331296114,
                        2.8450853365928737,
                        2.8681968400561364,
                        2.8913083435193982,
                        2.9144198469826605,
                        2.9375313504459228,
                        2.9606428539091851,
                        2.9837543573724474,
                        3.0068658608357097,
                        3.0299773642989716,
                        3.0530888677622339,
                        3.0762003712254966,
                        3.0993118746887589,
                        3.1224233781520212,
                        3.1455348816152831,
                        3.1686463850785453,
                        3.1917578885418076,
                        3.2148693920050699,
                        3.2379808954683322,
                        3.2610923989315941,
                        3.2842039023948564,
                        3.3073154058581191,
                        3.3304269093213814,
                        3.3535384127846437,
                        3.3766499162479056,
                        3.3997614197111679,
                        3.4228729231744301,
                        3.4459844266376924,
                        3.4690959301009547,
                        3.4922074335642166,
                        3.5153189370274789,
                        3.5384304404907416,
                        3.5615419439540039,
                        3.5846534474172662,
                        3.6077649508805281,
                        3.6308764543437904,
                        3.6539879578070527,
                        3.677099461270315,
                        3.7002109647335772,
                        3.7233224681968391,
                        3.7464339716601018,
                        3.7695454751233641,
                        3.7926569785866264,
                        3.8157684820498887,
                        3.8388799855131506,
                        3.8619914889764129,
                        3.8851029924396752,
                        3.9082144959029375,
                        3.9313259993661998,
                        3.9544375028294616,
                        3.9775490062927243,
                        4.0006605097559866,
                        4.0237720132192489,
                        4.0468835166825112,
                        4.0699950201457735,
                        4.0931065236090358,
                        4.1162180270722981,
                        4.1393295305355604,
                        4.1624410339988227,
                        4.185552537462085,
                        4.2086640409253473,
                        4.2317755443886096,
                        4.2548870478518719,
                        4.2779985513151342,
                        4.3011100547783965,
                        4.3242215582416588,
                        4.3473330617049211,
                        4.3704445651681834,
                        4.3935560686314457,
                        4.4166675720947079,
                        4.4397790755579694,
                        4.4628905790212317,
                        4.4860020824844939,
                        4.5091135859477571,
                        4.5322250894110194,
                        4.5553365928742808,
                        4.5784480963375431,
                        4.6015595998008054,
                        4.6246711032640677,
                        4.64778260672733,
                        4.6708941101905923,
                        4.6940056136538546,
                        4.7171171171171169,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            y2=numpy.array(
                numpy.array(
                    (
                        7158.8937047628706,
                        5004.9316715329842,
                        3694.8135594405553,
                        2839.0817737773841,
                        2249.5484991256844,
                        1826.2474529253161,
                        1512.0852402131027,
                        1272.530117074451,
                        1085.7010719257139,
                        937.18793256963431,
                        817.18705296685539,
                        718.84090024522891,
                        637.23614115501732,
                        568.77783627650172,
                        510.786630516309,
                        461.23254243400288,
                        418.55486713593308,
                        381.53776842598086,
                        349.2227154220239,
                        320.84580081746759,
                        295.79217667265311,
                        273.56246656674404,
                        253.74768704408271,
                        236.01030089701965,
                        220.06974682397441,
                        205.69127634981271,
                        192.67726151788867,
                        180.86036756672672,
                        170.09814691732512,
                        160.26872610227511,
                        151.26734021373531,
                        143.00352975000524,
                        135.39885901833182,
                        128.3850480674335,
                        121.90243465834536,
                        115.8987012784103,
                        110.32781625621425,
                        105.14914879151058,
                        100.32672600489369,
                        95.828606544860236,
                        91.626350312843897,
                        87.694567812436034,
                        84.010535746139567,
                        80.553867959068427,
                        77.306232806090279,
                        74.251109605447908,
                        71.373578121094255,
                        68.660136052082635,
                        66.098540350235496,
                        63.677668875882041,
                        61.387399466223862,
                        59.218503955912581,
                        57.162555073740741,
                        55.211844458103641,
                        53.35931029918256,
                        51.598473337326723,
                        49.923380132688088,
                        48.328552677103559,
                        46.808943550648046,
                        45.359895936370307,
                        43.977107900883965,
                        42.656600428510146,
                        41.39468876485914,
                        40.187956683990443,
                        39.033233343175098,
                        37.927572432105393,
                        36.868233360237589,
                        35.852664257722239,
                        34.878486592828409,
                        33.943481232542041,
                        33.045575793648069,
                        32.182833149542383,
                        31.353440973645949,
                        30.555702213931628,
                        29.788026404998764,
                        29.048921734577537,
                        28.336987790509955,
                        27.650908922311828,
                        26.989448158512342,
                        26.351441627222858,
                        25.73579343291075,
                        25.141470947240073,
                        24.567500476169894,
                        24.012963269341064,
                        23.476991841193676,
                        22.958766576292955,
                        22.457512594044672,
                        21.972496850393419,
                        21.503025456251216,
                        21.048441194330234,
                        20.608121217778795,
                        20.181474915566099,
                        19.767941930949156,
                        19.366990320602941,
                        18.978114843116536,
                        18.600835366568806,
                        18.234695385808415,
                        17.879260640884635,
                        17.53411782881939,
                        17.198873401581842,
                        16.873152443735915,
                        16.55659762378199,
                        16.248868213714225,
                        15.949639171768748,
                        15.658600283750749,
                        15.375455358703663,
                        15.099921475025205,
                        14.831728273446556,
                        14.570617293574582,
                        14.316341350956414,
                        14.068663951862252,
                        13.827358743198703,
                        13.592208995163054,
                        13.363007114430141,
                        13.139554185829615,
                        12.921659540623899,
                        12.709140349636851,
                        12.501821239611756,
                        12.299533931295143,
                        12.102116897851683,
                        11.909415042315253,
                        11.721279392873315,
                        11.537566814866656,
                        11.358139738464688,
                        11.182865901048743,
                        11.011618103402494,
                        10.844273978870199,
                        10.680715774700444,
                        10.520830144845792,
                        10.364507953537444,
                        10.211644088999289,
                        10.062137286707557,
                        9.91588996164114,
                        9.7728080490036398,
                        9.6328008529317479,
                        9.4957809027354614,
                        9.3616638162447625,
                        9.2303681698639899,
                        9.101815374960319,
                        8.9759295602358886,
                        8.8526374597548152,
                        8.7318683063165654,
                        8.6135537298858829,
                        8.4976276608070904,
                        8.3840262375469035,
                        8.2726877187252423,
                        8.1635523992077417,
                        8.0565625300470547,
                        7.9516622420725014,
                        7.848797472939288,
                        7.7479158974594009,
                        7.648966861046592,
                        7.5519013161173332,
                        7.456671761298697,
                        7.3632321833024221,
                        7.2715380013323854,
                        7.1815460139000722,
                        7.0932143479295533,
                        7.0065024100400572,
                        6.9213708399002529,
                        6.8377814655542144,
                        6.7556972606243528,
                        6.6750823033017541,
                        6.5959017370391155,
                        6.5181217328659393,
                        6.4417094532499517,
                        6.3666330174326555,
                        6.2928614681706652,
                        6.2203647398180992,
                        6.1491136276885232,
                        6.0790797586382022,
                        6.0102355628153186,
                        5.9425542465226266,
                        5.8760097661436923,
                        5.8105768030853397,
                        5.7462307396912662,
                        5.6829476360840934,
                        5.6207042078951268,
                        5.5594778048431799,
                        5.4992463901256654,
                        5.439988520586919,
                        5.3816833276304639,
                        5.3243104988434826,
                        5.2678502603032751,
                        5.2122833595369489,
                        5.157591049106899,
                        5.1037550707959838,
                        5.0507576403674586,
                        4.9985814328759144,
                        4.9472095685066089,
                        4.8966255989215437,
                        4.8468134940916663,
                        4.797757629595548,
                        4.7494427743657095,
                        4.7018540788646925,
                        4.6549770636737042,
                        4.6087976084774764,
                        4.5633019414297058,
                        4.5184766288841098,
                        4.4743085654767931,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            xin=0.22251193896599297,
            expected_x=0.624584480717571,
        ),
    ),
)
def test_intersect(intersectparam, stellarator):
    """
    Automatically generated Regression Unit Test for intersect.

    This test was generated using data from tests/regression/scenarios/stellarator/IN.DAT.

    :param intersectparam: the data used to mock and assert in this test.
    :type intersectparam: intersectparam
    """

    x = stellarator.intersect(
        x1=intersectparam.x1,
        y1=intersectparam.y1,
        x2=intersectparam.x2,
        y2=intersectparam.y2,
        xin=intersectparam.xin,
    )

    assert x == pytest.approx(intersectparam.expected_x)


class StdlimParam(NamedTuple):
    dene: Any = None

    dnla: Any = None

    dnelimt: Any = None

    bt: Any = None

    powht: Any = None

    rmajor: Any = None

    rminor: Any = None

    expected_dnelimt: Any = None

    expected_dlimit: Any = None


@pytest.mark.parametrize(
    "stdlimparam",
    (
        StdlimParam(
            dene=2.0914e20,
            dnla=2.357822619799476e20,
            dnelimt=0,
            bt=5.5,
            powht=432.20449197454559,
            rmajor=22,
            rminor=1.7842660178426601,
            expected_dnelimt=1.2918765671497731e20,
            expected_dlimit=1.2918765671497731e20,
        ),
        StdlimParam(
            dene=2.0914e20,
            dnla=2.357822619799476e20,
            dnelimt=1.2918765671497731e20,
            bt=5.5,
            powht=431.98698920075435,
            rmajor=22,
            rminor=1.7842660178426601,
            expected_dnelimt=1.2915514639846759e20,
            expected_dlimit=1.2915514639846759e20,
        ),
    ),
)
def test_stdlim(stdlimparam, monkeypatch, stellarator):
    """
    Automatically generated Regression Unit Test for stdlim.

    This test was generated using data from tests/regression/scenarios/stellarator/IN.DAT.

    :param stdlimparam: the data used to mock and assert in this test.
    :type stdlimparam: stdlimparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "dene", stdlimparam.dene)

    monkeypatch.setattr(physics_variables, "dnla", stdlimparam.dnla)

    monkeypatch.setattr(physics_variables, "dnelimt", stdlimparam.dnelimt)

    dlimit = stellarator.stdlim(
        bt=stdlimparam.bt,
        powht=stdlimparam.powht,
        rmajor=stdlimparam.rmajor,
        rminor=stdlimparam.rminor,
    )

    assert physics_variables.dnelimt == pytest.approx(stdlimparam.expected_dnelimt)

    assert dlimit == pytest.approx(stdlimparam.expected_dlimit)


class StdlimEcrhParam(NamedTuple):
    ipedestal: Any = None

    bt_input: Any = None

    gyro_frequency_max: Any = None

    expected_dlimit_ecrh: Any = None

    expected_bt_max: Any = None


@pytest.mark.parametrize(
    "stdlimecrhparam",
    (
        StdlimEcrhParam(
            ipedestal=0,
            bt_input=6.9100000000000001,
            gyro_frequency_max=400000000000,
            expected_dlimit_ecrh=4.6472737339514113e20,
            expected_bt_max=14.279966607226331,
        ),
    ),
)
def test_stdlim_ecrh(stdlimecrhparam, monkeypatch, stellarator):
    """
    Automatically generated Regression Unit Test for stdlim_ecrh.

    This test was generated using data from tests/regression/scenarios/stellarator_config/IN.DAT.

    :param stdlimecrhparam: the data used to mock and assert in this test.
    :type stdlimecrhparam: stdlimecrhparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "ipedestal", stdlimecrhparam.ipedestal)

    dlimit_ecrh, bt_max = stellarator.stdlim_ecrh(
        bt_input=stdlimecrhparam.bt_input,
        gyro_frequency_max=stdlimecrhparam.gyro_frequency_max,
    )

    assert dlimit_ecrh == pytest.approx(stdlimecrhparam.expected_dlimit_ecrh)
    assert bt_max == pytest.approx(stdlimecrhparam.expected_bt_max)


class StCalcEffChiParam(NamedTuple):
    te0: Any = None

    ne0: Any = None

    falpha: Any = None

    palppv: Any = None

    pcoreradpv: Any = None

    alphan: Any = None

    alphat: Any = None

    vol: Any = None

    sarea: Any = None

    rminor: Any = None

    coreradius: Any = None

    stella_config_rminor_ref: Any = None

    f_r: Any = None

    expected_output: Any = None


@pytest.mark.parametrize(
    "stcalceffchiparam",
    (
        StCalcEffChiParam(
            te0=19.108573496973477,
            ne0=3.4479000000000007e20,
            falpha=0.95000000000000007,
            palppv=1.2629524018077414,
            pcoreradpv=0.10762698429338043,
            alphan=0.35000000000000003,
            alphat=1.2,
            vol=1385.8142655379029,
            sarea=1926.0551116585129,
            rminor=1.7863900994187722,
            coreradius=0.60000000000000009,
            stella_config_rminor_ref=1.80206932,
            f_r=0.99129932482229,
            expected_output=0.2620230359599852,
            # expected_output=0.26206561772729992, used old e_
        ),
        StCalcEffChiParam(
            te0=17.5,
            ne0=3.4479000000000007e20,
            falpha=0.95000000000000007,
            palppv=1.0570658694225301,
            pcoreradpv=0.1002475669217598,
            alphan=0.35000000000000003,
            alphat=1.2,
            vol=1385.8142655379029,
            sarea=1926.0551116585129,
            rminor=1.7863900994187722,
            coreradius=0.60000000000000009,
            stella_config_rminor_ref=1.80206932,
            f_r=0.99129932482229,
            expected_output=0.2368034193234161,
            # expected_output=0.23684190261197124, used old e_
        ),
    ),
)
def test_st_calc_eff_chi(stcalceffchiparam, monkeypatch, stellarator):
    """
    Automatically generated Regression Unit Test for st_calc_eff_chi.

    This test was generated using data from tests/regression/scenarios/stellarator_config/IN.DAT.

    :param stcalceffchiparam: the data used to mock and assert in this test.
    :type stcalceffchiparam: stcalceffchiparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(physics_variables, "te0", stcalceffchiparam.te0)

    monkeypatch.setattr(physics_variables, "ne0", stcalceffchiparam.ne0)

    monkeypatch.setattr(physics_variables, "falpha", stcalceffchiparam.falpha)

    monkeypatch.setattr(physics_variables, "palppv", stcalceffchiparam.palppv)

    monkeypatch.setattr(physics_variables, "pcoreradpv", stcalceffchiparam.pcoreradpv)

    monkeypatch.setattr(physics_variables, "alphan", stcalceffchiparam.alphan)

    monkeypatch.setattr(physics_variables, "alphat", stcalceffchiparam.alphat)

    monkeypatch.setattr(physics_variables, "vol", stcalceffchiparam.vol)

    monkeypatch.setattr(physics_variables, "sarea", stcalceffchiparam.sarea)

    monkeypatch.setattr(physics_variables, "rminor", stcalceffchiparam.rminor)

    monkeypatch.setattr(
        impurity_radiation_module, "coreradius", stcalceffchiparam.coreradius
    )

    monkeypatch.setattr(
        stellarator_configuration,
        "stella_config_rminor_ref",
        stcalceffchiparam.stella_config_rminor_ref,
    )

    monkeypatch.setattr(stellarator_module, "f_r", stcalceffchiparam.f_r)

    output = stellarator.st_calc_eff_chi()

    assert output == pytest.approx(stcalceffchiparam.expected_output)
