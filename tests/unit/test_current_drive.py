import pytest
from typing import NamedTuple, Any


from process.fortran import (
    current_drive_variables,
    cost_variables,
    physics_variables,
    heat_transport_variables,
)
from process.current_drive import CurrentDrive


@pytest.fixture
def current_drive():
    """Provides CurrentDrive object for testing.

    :returns current_drive: initialised CurrentDrive object
    :rtype: process.current_drive.CurrentDrive
    """
    return CurrentDrive()


class CudrivParam(NamedTuple):
    pinjwpfix: Any = None

    pinjwp: Any = None

    echpwr: Any = None

    pnbeam: Any = None

    plhybd: Any = None

    cnbeam: Any = None

    porbitlossmw: Any = None

    iefrf: Any = None

    iefrffix: Any = None

    pheat: Any = None

    pheatfix: Any = None

    pinjfixmw: Any = None

    irfcd: Any = None

    feffcd: Any = None

    fpion: Any = None

    nbshinef: Any = None

    gamcd: Any = None

    gamma_ecrh: Any = None

    etalh: Any = None

    etacd: Any = None

    etacdfix: Any = None

    etaech: Any = None

    forbitloss: Any = None

    pinjmw: Any = None

    pwpnb: Any = None

    etanbi: Any = None

    enbeam: Any = None

    effcd: Any = None

    pwplh: Any = None

    echwpow: Any = None

    pnbitot: Any = None

    nbshinemw: Any = None

    pinjemw: Any = None

    pinjimw: Any = None

    bigq: Any = None

    bootipf: Any = None

    bscfmax: Any = None

    taubeam: Any = None

    pinjalw: Any = None

    nbshield: Any = None

    frbeam: Any = None

    rtanbeam: Any = None

    rtanmax: Any = None

    diaipf: Any = None

    psipf: Any = None

    plasipf: Any = None

    harnum: Any = None

    xi_ebw: Any = None

    dene: Any = None

    te: Any = None

    rmajor: Any = None

    ten: Any = None

    zeff: Any = None

    dlamee: Any = None

    beta: Any = None

    rhopedt: Any = None

    rhopedn: Any = None

    te0: Any = None

    teped: Any = None

    tesep: Any = None

    alphat: Any = None

    alphan: Any = None

    ne0: Any = None

    nesep: Any = None

    neped: Any = None

    bt: Any = None

    rminor: Any = None

    tbeta: Any = None

    plascur: Any = None

    ipedestal: Any = None

    faccd: Any = None

    ignite: Any = None

    pohmmw: Any = None

    powfmw: Any = None

    facoh: Any = None

    fvsbrnni: Any = None

    startupratio: Any = None

    iprint: Any = None

    outfile: Any = None

    expected_pinjwp: Any = None

    expected_echpwr: Any = None

    expected_gamcd: Any = None

    expected_etacd: Any = None

    expected_pinjmw: Any = None

    expected_effcd: Any = None

    expected_echwpow: Any = None

    expected_pinjemw: Any = None

    expected_bigq: Any = None


@pytest.mark.parametrize(
    "cudrivparam",
    (
        CudrivParam(
            pinjwpfix=0,
            pinjwp=0,
            echpwr=0,
            pnbeam=0,
            plhybd=0,
            cnbeam=0,
            porbitlossmw=0,
            iefrf=10,
            iefrffix=0,
            pheat=75,
            pheatfix=0,
            pinjfixmw=0,
            irfcd=1,
            feffcd=1,
            fpion=0.5,
            nbshinef=0,
            gamcd=0,
            gamma_ecrh=0.30000000000000004,
            etalh=0.29999999999999999,
            etacd=0,
            etacdfix=0,
            etaech=0.5,
            forbitloss=0,
            pinjmw=0,
            pwpnb=0,
            etanbi=0.29999999999999999,
            enbeam=1000,
            effcd=0,
            pwplh=0,
            echwpow=0,
            pnbitot=0,
            nbshinemw=0,
            pinjemw=0,
            pinjimw=0,
            bigq=0,
            bootipf=0.27635918746616817,
            bscfmax=0.95000000000000007,
            taubeam=0,
            pinjalw=200,
            nbshield=0.5,
            frbeam=1.05,
            rtanbeam=0,
            rtanmax=0,
            diaipf=0,
            psipf=0,
            plasipf=0.27635918746616817,
            harnum=1,
            xi_ebw=0.80000000000000004,
            dene=7.5e19,
            te=12,
            rmajor=8,
            ten=12.626131115905864,
            zeff=2.0909945616489103,
            dlamee=17.510652035055571,
            beta=0.030000000000000006,
            rhopedt=0.94000000000000006,
            rhopedn=0.94000000000000006,
            te0=24.402321098330372,
            teped=5.5,
            tesep=0.10000000000000001,
            alphat=1.45,
            alphan=1,
            ne0=8.515060981068918e19,
            nesep=4.1177885154594193e19,
            neped=7.000240476281013e19,
            bt=5.7000000000000002,
            rminor=2.6666666666666665,
            tbeta=2,
            plascur=18398455.678867526,
            ipedestal=1,
            faccd=0.12364081253383186,
            ignite=0,
            pohmmw=0,
            powfmw=0,
            facoh=0.59999999999999998,
            fvsbrnni=0.40000000000000002,
            startupratio=1,
            iprint=0,
            outfile=11,
            expected_pinjwp=240.99200038011492,
            expected_echpwr=120.49600019005746,
            expected_gamcd=0.30000000000000004,
            expected_etacd=0.5,
            expected_pinjmw=120.49600019005746,
            expected_effcd=0.05000000000000001,
            expected_echwpow=240.99200038011492,
            expected_pinjemw=120.49600019005746,
            expected_bigq=0,
        ),
        CudrivParam(
            pinjwpfix=0,
            pinjwp=240.99200038011492,
            echpwr=120.49600019005746,
            pnbeam=0,
            plhybd=0,
            cnbeam=0,
            porbitlossmw=0,
            iefrf=10,
            iefrffix=0,
            pheat=75,
            pheatfix=0,
            pinjfixmw=0,
            irfcd=1,
            feffcd=1,
            fpion=0.5,
            nbshinef=0,
            gamcd=0.30000000000000004,
            gamma_ecrh=0.30000000000000004,
            etalh=0.29999999999999999,
            etacd=0.5,
            etacdfix=0,
            etaech=0.5,
            forbitloss=0,
            pinjmw=120.49600019005746,
            pwpnb=0,
            etanbi=0.29999999999999999,
            enbeam=1000,
            effcd=0.05000000000000001,
            pwplh=0,
            echwpow=240.99200038011492,
            pnbitot=0,
            nbshinemw=0,
            pinjemw=120.49600019005746,
            pinjimw=0,
            bigq=0,
            bootipf=0.27635918746616817,
            bscfmax=0.95000000000000007,
            taubeam=0,
            pinjalw=200,
            nbshield=0.5,
            frbeam=1.05,
            rtanbeam=8.4000000000000004,
            rtanmax=13.179564451855533,
            diaipf=0,
            psipf=0,
            plasipf=0.27635918746616817,
            harnum=1,
            xi_ebw=0.80000000000000004,
            dene=7.5e19,
            te=12,
            rmajor=8,
            ten=12.626131115905864,
            zeff=2.0909945616489103,
            dlamee=17.510652035055571,
            beta=0.030000000000000006,
            rhopedt=0.94000000000000006,
            rhopedn=0.94000000000000006,
            te0=24.402321098330372,
            teped=5.5,
            tesep=0.10000000000000001,
            alphat=1.45,
            alphan=1,
            ne0=8.515060981068918e19,
            nesep=4.1177885154594193e19,
            neped=7.000240476281013e19,
            bt=5.7000000000000002,
            rminor=2.6666666666666665,
            tbeta=2,
            plascur=18398455.678867526,
            ipedestal=1,
            faccd=0.12364081253383186,
            ignite=0,
            pohmmw=0.76707314489379119,
            powfmw=1051.6562748933977,
            facoh=0.59999999999999998,
            fvsbrnni=0.40000000000000002,
            startupratio=1,
            iprint=0,
            outfile=11,
            expected_pinjwp=240.99200038011492,
            expected_echpwr=120.49600019005746,
            expected_gamcd=0.30000000000000004,
            expected_etacd=0.5,
            expected_pinjmw=120.49600019005746,
            expected_effcd=0.05000000000000001,
            expected_echwpow=240.99200038011492,
            expected_pinjemw=120.49600019005746,
            expected_bigq=8.6725187311435423,
        ),
    ),
)
def test_cudriv(cudrivparam, monkeypatch, current_drive):
    """
    Automatically generated Regression Unit Test for cudriv.

    This test was generated using data from tests/regression/scenarios/large-tokamak/IN.DAT.

    :param cudrivparam: the data used to mock and assert in this test.
    :type cudrivparam: cudrivparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(heat_transport_variables, "pinjwpfix", cudrivparam.pinjwpfix)

    monkeypatch.setattr(heat_transport_variables, "pinjwp", cudrivparam.pinjwp)

    monkeypatch.setattr(current_drive_variables, "echpwr", cudrivparam.echpwr)

    monkeypatch.setattr(current_drive_variables, "pnbeam", cudrivparam.pnbeam)

    monkeypatch.setattr(current_drive_variables, "plhybd", cudrivparam.plhybd)

    monkeypatch.setattr(current_drive_variables, "cnbeam", cudrivparam.cnbeam)

    monkeypatch.setattr(
        current_drive_variables, "porbitlossmw", cudrivparam.porbitlossmw
    )

    monkeypatch.setattr(current_drive_variables, "iefrf", cudrivparam.iefrf)

    monkeypatch.setattr(current_drive_variables, "iefrffix", cudrivparam.iefrffix)

    monkeypatch.setattr(current_drive_variables, "pheat", cudrivparam.pheat)

    monkeypatch.setattr(current_drive_variables, "pheatfix", cudrivparam.pheatfix)

    monkeypatch.setattr(current_drive_variables, "pinjfixmw", cudrivparam.pinjfixmw)

    monkeypatch.setattr(current_drive_variables, "irfcd", cudrivparam.irfcd)

    monkeypatch.setattr(current_drive_variables, "feffcd", cudrivparam.feffcd)

    monkeypatch.setattr(current_drive_variables, "fpion", cudrivparam.fpion)

    monkeypatch.setattr(current_drive_variables, "nbshinef", cudrivparam.nbshinef)

    monkeypatch.setattr(current_drive_variables, "gamcd", cudrivparam.gamcd)

    monkeypatch.setattr(current_drive_variables, "gamma_ecrh", cudrivparam.gamma_ecrh)

    monkeypatch.setattr(current_drive_variables, "etalh", cudrivparam.etalh)

    monkeypatch.setattr(current_drive_variables, "etacd", cudrivparam.etacd)

    monkeypatch.setattr(current_drive_variables, "etacdfix", cudrivparam.etacdfix)

    monkeypatch.setattr(current_drive_variables, "etaech", cudrivparam.etaech)

    monkeypatch.setattr(current_drive_variables, "forbitloss", cudrivparam.forbitloss)

    monkeypatch.setattr(current_drive_variables, "pinjmw", cudrivparam.pinjmw)

    monkeypatch.setattr(current_drive_variables, "pwpnb", cudrivparam.pwpnb)

    monkeypatch.setattr(current_drive_variables, "etanbi", cudrivparam.etanbi)

    monkeypatch.setattr(current_drive_variables, "enbeam", cudrivparam.enbeam)

    monkeypatch.setattr(current_drive_variables, "effcd", cudrivparam.effcd)

    monkeypatch.setattr(current_drive_variables, "pwplh", cudrivparam.pwplh)

    monkeypatch.setattr(current_drive_variables, "echwpow", cudrivparam.echwpow)

    monkeypatch.setattr(current_drive_variables, "pnbitot", cudrivparam.pnbitot)

    monkeypatch.setattr(current_drive_variables, "nbshinemw", cudrivparam.nbshinemw)

    monkeypatch.setattr(current_drive_variables, "pinjemw", cudrivparam.pinjemw)

    monkeypatch.setattr(current_drive_variables, "pinjimw", cudrivparam.pinjimw)

    monkeypatch.setattr(current_drive_variables, "bigq", cudrivparam.bigq)

    monkeypatch.setattr(current_drive_variables, "bootipf", cudrivparam.bootipf)

    monkeypatch.setattr(current_drive_variables, "bscfmax", cudrivparam.bscfmax)

    monkeypatch.setattr(current_drive_variables, "taubeam", cudrivparam.taubeam)

    monkeypatch.setattr(current_drive_variables, "pinjalw", cudrivparam.pinjalw)

    monkeypatch.setattr(current_drive_variables, "nbshield", cudrivparam.nbshield)

    monkeypatch.setattr(current_drive_variables, "frbeam", cudrivparam.frbeam)

    monkeypatch.setattr(current_drive_variables, "rtanbeam", cudrivparam.rtanbeam)

    monkeypatch.setattr(current_drive_variables, "rtanmax", cudrivparam.rtanmax)

    monkeypatch.setattr(current_drive_variables, "diaipf", cudrivparam.diaipf)

    monkeypatch.setattr(current_drive_variables, "psipf", cudrivparam.psipf)

    monkeypatch.setattr(current_drive_variables, "plasipf", cudrivparam.plasipf)

    monkeypatch.setattr(current_drive_variables, "harnum", cudrivparam.harnum)

    monkeypatch.setattr(current_drive_variables, "xi_ebw", cudrivparam.xi_ebw)

    monkeypatch.setattr(physics_variables, "dene", cudrivparam.dene)

    monkeypatch.setattr(physics_variables, "te", cudrivparam.te)

    monkeypatch.setattr(physics_variables, "rmajor", cudrivparam.rmajor)

    monkeypatch.setattr(physics_variables, "ten", cudrivparam.ten)

    monkeypatch.setattr(physics_variables, "zeff", cudrivparam.zeff)

    monkeypatch.setattr(physics_variables, "dlamee", cudrivparam.dlamee)

    monkeypatch.setattr(physics_variables, "beta", cudrivparam.beta)

    monkeypatch.setattr(physics_variables, "rhopedt", cudrivparam.rhopedt)

    monkeypatch.setattr(physics_variables, "rhopedn", cudrivparam.rhopedn)

    monkeypatch.setattr(physics_variables, "te0", cudrivparam.te0)

    monkeypatch.setattr(physics_variables, "teped", cudrivparam.teped)

    monkeypatch.setattr(physics_variables, "tesep", cudrivparam.tesep)

    monkeypatch.setattr(physics_variables, "alphat", cudrivparam.alphat)

    monkeypatch.setattr(physics_variables, "alphan", cudrivparam.alphan)

    monkeypatch.setattr(physics_variables, "ne0", cudrivparam.ne0)

    monkeypatch.setattr(physics_variables, "nesep", cudrivparam.nesep)

    monkeypatch.setattr(physics_variables, "neped", cudrivparam.neped)

    monkeypatch.setattr(physics_variables, "bt", cudrivparam.bt)

    monkeypatch.setattr(physics_variables, "rminor", cudrivparam.rminor)

    monkeypatch.setattr(physics_variables, "tbeta", cudrivparam.tbeta)

    monkeypatch.setattr(physics_variables, "plascur", cudrivparam.plascur)

    monkeypatch.setattr(physics_variables, "ipedestal", cudrivparam.ipedestal)

    monkeypatch.setattr(physics_variables, "faccd", cudrivparam.faccd)

    monkeypatch.setattr(physics_variables, "ignite", cudrivparam.ignite)

    monkeypatch.setattr(physics_variables, "pohmmw", cudrivparam.pohmmw)

    monkeypatch.setattr(physics_variables, "powfmw", cudrivparam.powfmw)

    monkeypatch.setattr(physics_variables, "facoh", cudrivparam.facoh)

    monkeypatch.setattr(physics_variables, "fvsbrnni", cudrivparam.fvsbrnni)

    monkeypatch.setattr(cost_variables, "startupratio", cudrivparam.startupratio)

    current_drive.cudriv(output=False)

    assert heat_transport_variables.pinjwp == pytest.approx(cudrivparam.expected_pinjwp)

    assert current_drive_variables.echpwr == pytest.approx(cudrivparam.expected_echpwr)

    assert current_drive_variables.gamcd == pytest.approx(cudrivparam.expected_gamcd)

    assert current_drive_variables.etacd == pytest.approx(cudrivparam.expected_etacd)

    assert current_drive_variables.pinjmw == pytest.approx(cudrivparam.expected_pinjmw)

    assert current_drive_variables.effcd == pytest.approx(cudrivparam.expected_effcd)

    assert current_drive_variables.echwpow == pytest.approx(
        cudrivparam.expected_echwpow
    )

    assert current_drive_variables.pinjemw == pytest.approx(
        cudrivparam.expected_pinjemw
    )

    assert current_drive_variables.bigq == pytest.approx(cudrivparam.expected_bigq)


def test_sigbeam(current_drive):
    assert current_drive.sigbeam(
        1e3, 13.07, 8.0e-1, 0.1, 1e-4, 1e-4, 1e-4
    ) == pytest.approx(2.013589662302492e-11)
