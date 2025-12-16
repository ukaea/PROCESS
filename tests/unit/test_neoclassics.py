from typing import Any, NamedTuple

import numpy as np
import pytest

from process.data_structure import neoclassics_variables, physics_variables
from process.stellarator import Neoclassics


@pytest.fixture
def neoclassics():
    """Provides Neoclassics object for testing.

    :returns: initialised Neoclassics object
    :rtype: process.stellerator.Neoclassics
    """
    return Neoclassics()


class InitNeoclassicsParam(NamedTuple):
    nd_plasma_electron_on_axis: Any = None
    temp_plasma_electron_on_axis_kev: Any = None
    alphan: Any = None
    alphat: Any = None
    temp_plasma_ion_on_axis_kev: Any = None
    nd_plasma_ions_on_axis: Any = None
    f_plasma_fuel_deuterium: Any = None
    nd_plasma_alphas_vol_avg: Any = None
    rminor: Any = None
    rmajor: Any = None
    b_plasma_toroidal_on_axis: Any = None
    te: Any = None
    ti: Any = None
    nd_plasma_electrons_vol_avg: Any = None
    nd_plasma_fuel_ions_vol_avg: Any = None
    densities: Any = None
    temperatures: Any = None
    dr_densities: Any = None
    dr_temperatures: Any = None
    roots: Any = None
    weights: Any = None
    nu: Any = None
    nu_star: Any = None
    nu_star_averaged: Any = None
    vd: Any = None
    iota: Any = None
    q_flux: Any = None
    eps_eff: Any = None
    r_eff: Any = None
    r_effin: Any = None
    eps_effin: Any = None
    iotain: Any = None
    expected_densities: Any = None
    expected_temperatures: Any = None
    expected_dr_densities: Any = None
    expected_dr_temperatures: Any = None
    expected_roots: Any = None
    expected_weights: Any = None
    expected_nu: Any = None
    expected_nu_star: Any = None
    expected_nu_star_averaged: Any = None
    expected_vd: Any = None
    expected_q_flux: Any = None


@pytest.mark.parametrize(
    "initneoclassicsparam",
    (
        InitNeoclassicsParam(
            nd_plasma_electron_on_axis=2.7956610000000002e20,
            temp_plasma_electron_on_axis_kev=13.241800000000001,
            alphan=0.35000000000000003,
            alphat=1.2,
            temp_plasma_ion_on_axis_kev=12.579710000000002,
            nd_plasma_ions_on_axis=2.3930858160000005e20,
            f_plasma_fuel_deuterium=0.5,
            nd_plasma_alphas_vol_avg=2.9820384000000004e19,
            rminor=1.7993820274145451,
            rmajor=22.16,
            b_plasma_toroidal_on_axis=5.2400000000000002,
            te=6.0190000000000001,
            ti=5.7180500000000007,
            nd_plasma_electrons_vol_avg=2.07086e20,
            nd_plasma_fuel_ions_vol_avg=1.47415411616e20,
            densities=np.array(np.array((0, 0, 0, 0), order="F"), order="F").transpose(),
            temperatures=np.array(
                np.array((0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            dr_densities=np.array(
                np.array((0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            dr_temperatures=np.array(
                np.array((0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            roots=np.array(
                np.array(
                    (
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            weights=np.array(
                np.array(
                    (
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                        0,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            nu=np.array(
                (
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                ),
                order="F",
            ).transpose(),
            nu_star=np.array(
                (
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                ),
                order="F",
            ).transpose(),
            nu_star_averaged=np.array(
                np.array((0, 0, 0, 0), order="F"), order="F"
            ).transpose(),
            vd=np.array(
                (
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                    (0, 0, 0, 0),
                ),
                order="F",
            ).transpose(),
            iota=1,
            q_flux=np.array(np.array((0, 0, 0, 0), order="F"), order="F").transpose(),
            eps_eff=1.0000000000000001e-05,
            r_eff=0,
            r_effin=0.59999999999999998,
            eps_effin=0.01464553,
            iotain=0.90000000000000002,
            expected_densities=(
                2.391373976836772e20,
                1.0235080620861384e20,
                1.0235080620861384e20,
                3.4435785266449519e19,
            ),
            expected_temperatures=(
                1.2416608937807636e-15,
                1.1795778490917256e-15,
                1.1795778490917256e-15,
                1.1795778490917256e-15,
            ),
            expected_dr_densities=(
                -8.7215452215783653e19,
                -3.7328213548355412e19,
                -3.7328213548355412e19,
                -1.2559025119072848e19,
            ),
            expected_dr_temperatures=(
                -1.5526091560561597e-15,
                -1.4749786982533515e-15,
                -1.4749786982533515e-15,
                -1.4749786982533515e-15,
            ),
            expected_roots=np.array(
                np.array(
                    (
                        0.047407180540805262,
                        0.24992391675315939,
                        0.61483345439276837,
                        1.1431958256661015,
                        1.8364545546225723,
                        2.6965218745572161,
                        3.7258145077795093,
                        4.9272937658498819,
                        6.3045155909650736,
                        7.8616932933702603,
                        9.6037759854792633,
                        11.536546597956139,
                        13.666744693064235,
                        16.002221188981068,
                        18.55213484014315,
                        21.327204321783128,
                        24.340035764532693,
                        27.605554796780961,
                        31.141586701111237,
                        34.969652008249071,
                        39.11608494906789,
                        43.613652908484831,
                        48.5039861638042,
                        53.841385406507506,
                        59.699121859235497,
                        66.18061779443849,
                        73.441238595559881,
                        81.736810506727679,
                        91.556466522536837,
                        104.15752443105889,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_weights=(
                0.11604408602043889,
                0.22085112475067714,
                0.24139982758785372,
                0.19463676844641709,
                0.12372841596687649,
                0.063678780368986609,
                0.026860475273379727,
                0.0093380708816039257,
                0.0026806968913368197,
                0.00063512912194085564,
                0.00012390745990688301,
                1.9828788438952331e-05,
                2.5893509291313925e-06,
                2.7409428405360132e-07,
                2.3328311650257382e-08,
                1.580745574778328e-09,
                8.4274791230567164e-11,
                3.4851612349078554e-12,
                1.0990180597534515e-13,
                2.5883126649590802e-15,
                4.437838059840029e-17,
                5.3659183082120453e-19,
                4.3939468922916045e-21,
                2.3114097943885432e-23,
                7.2745884982922481e-26,
                1.2391497014482679e-28,
                9.8323750831058875e-32,
                2.8423235534027009e-35,
                1.8786080317495154e-39,
                8.7459804404650116e-45,
            ),
            expected_nu=np.array(
                (
                    (
                        3695464.938733798,
                        11862.660224787611,
                        7933.188679624609,
                        23814.418394017364,
                    ),
                    (
                        343316.72274464089,
                        2123.530294833155,
                        1450.1975058356002,
                        4403.93523774,
                    ),
                    (
                        97325.601943625719,
                        782.90996572377162,
                        550.46150951687332,
                        1698.5476568456329,
                    ),
                    (
                        40967.943564057088,
                        373.3569294227749,
                        270.83070812853452,
                        851.34936648572614,
                    ),
                    (
                        21050.27896218509,
                        204.31685155512426,
                        152.39427303068328,
                        487.91926000009528,
                    ),
                    (
                        12196.141744781391,
                        122.4627549323195,
                        93.37151291439298,
                        303.82292065987531,
                    ),
                    (
                        7662.855304124023,
                        78.497186099853351,
                        60.822352763230512,
                        200.532111745934,
                    ),
                    (
                        5107.6181124745644,
                        53.008588809359139,
                        41.543171451498225,
                        138.36568682166154,
                    ),
                    (
                        3562.1896935110958,
                        37.316304271771997,
                        29.478466134989471,
                        98.931218624509825,
                    ),
                    (
                        2575.1493019063892,
                        27.17114569450117,
                        21.584157832183529,
                        72.845381503203299,
                    ),
                    (
                        1916.5626815978521,
                        20.341237798241099,
                        16.222662407259669,
                        54.977961633660357,
                    ),
                    (
                        1461.0373149373809,
                        15.58437516643423,
                        12.464365949260028,
                        42.371363838150558,
                    ),
                    (
                        1136.3293377214109,
                        12.174519399214853,
                        9.7573846405961202,
                        33.245772325916597,
                    ),
                    (
                        898.87370477934678,
                        9.6691588841578611,
                        7.7612618185478333,
                        26.490662053885323,
                    ),
                    (
                        721.36428722004121,
                        7.7886363224542317,
                        6.2588248584263582,
                        21.391037801763538,
                    ),
                    (
                        586.10784951621235,
                        6.3505475050965412,
                        5.1074008781905587,
                        17.473656809416916,
                    ),
                    (
                        481.30424894598372,
                        5.2325813464582778,
                        4.2107797561372902,
                        14.417477347614785,
                    ),
                    (
                        398.88189106951904,
                        4.3506921903280844,
                        3.5025544153548962,
                        11.999840321777146,
                    ),
                    (
                        333.19454764267971,
                        3.6458582076791037,
                        2.9359214172811088,
                        10.063206020290355,
                    ),
                    (
                        280.21278073497348,
                        3.0758139209321369,
                        2.4772660955148873,
                        8.494068483233054,
                    ),
                    (
                        237.00858514975135,
                        2.6097529846718648,
                        2.1020279060352189,
                        7.2092769760436468,
                    ),
                    (
                        201.41883814132996,
                        2.2248531750213925,
                        1.7919748195144374,
                        6.1469708670487027,
                    ),
                    (
                        171.82051852815496,
                        1.9039461312634718,
                        1.5333681490810633,
                        5.260448904005937,
                    ),
                    (
                        146.97724533175426,
                        1.6339205441745102,
                        1.3157006442816281,
                        4.5139376970122962,
                    ),
                    (
                        125.93193667489498,
                        1.4046010159515565,
                        1.130807807982801,
                        3.8796018446013152,
                    ),
                    (
                        107.92910284917983,
                        1.2079326736732123,
                        0.97221975350348488,
                        3.335356496099501,
                    ),
                    (
                        92.354617547302752,
                        1.0373444948341899,
                        0.83465363159517136,
                        2.8631484655280879,
                    ),
                    (
                        78.68026674668576,
                        0.88715494590298993,
                        0.71353811705010739,
                        2.447338816503998,
                    ),
                    (
                        66.385967521689778,
                        0.75172062458023658,
                        0.60432964829100722,
                        2.0723651470862698,
                    ),
                    (
                        54.725799270711548,
                        0.62283558834560981,
                        0.50041870619678697,
                        1.7155597347573166,
                    ),
                ),
                order="F",
            ).transpose(),
            expected_nu_star=np.array(
                (
                    (
                        10.191074937923535,
                        2.0331966316122507,
                        1.6652931278769849,
                        5.7723445649189209,
                    ),
                    (
                        0.41298040929536733,
                        0.15851635531154254,
                        0.13258314893305531,
                        0.4649127239401013,
                    ),
                    (
                        0.074848077541108754,
                        0.037260863337660155,
                        0.03208584666735572,
                        0.11432298738955823,
                    ),
                    (
                        0.023197077266323505,
                        0.01303119334609255,
                        0.011577190579052517,
                        0.042022540114017584,
                    ),
                    (
                        0.0094525677652312987,
                        0.0056264565943291555,
                        0.0051397801137085886,
                        0.019001722106245209,
                    ),
                    (
                        0.0045482132226884351,
                        0.0027830670017342913,
                        0.0025988361335270451,
                        0.0097645901733738206,
                    ),
                    (
                        0.0024492474025462241,
                        0.001517630205259042,
                        0.0014401895657356347,
                        0.0054828859338448806,
                    ),
                    (
                        0.0014317922467676544,
                        0.00089118019395046364,
                        0.00085538801619877453,
                        0.003289732501027202,
                    ),
                    (
                        0.00089132611440290003,
                        0.00055462245458156732,
                        0.00053659613981289387,
                        0.0020794309035960076,
                    ),
                    (
                        0.00058320422633617944,
                        0.00036163954478676596,
                        0.00035184101734703275,
                        0.0013711394734907242,
                    ),
                    (
                        0.0003973238458688693,
                        0.00024495361972729277,
                        0.00023926078180568864,
                        0.00093628060138880047,
                    ),
                    (
                        0.00027986715624797473,
                        0.00017123030253161586,
                        0.00016772730703600486,
                        0.00065837564428247015,
                    ),
                    (
                        0.00020271700094666678,
                        0.0001228995355882499,
                        0.00012063507255191103,
                        0.00047461779631969821,
                    ),
                    (
                        0.00015035093166913791,
                        9.0205186378467152e-05,
                        8.8678069715659373e-05,
                        0.0003494972842656118,
                    ),
                    (
                        0.00011379147025090292,
                        6.7483810604850391e-05,
                        6.6415788403404255e-05,
                        0.00026210609671657975,
                    ),
                    (
                        8.7635514886206191e-05,
                        5.1319413760932046e-05,
                        5.0548705442571374e-05,
                        0.00019969180051899852,
                    ),
                    (
                        6.8517128585695404e-05,
                        3.9581790112080437e-05,
                        3.9010425484079875e-05,
                        0.00015423159537601957,
                    ),
                    (
                        5.4275274472921955e-05,
                        3.0903164171056251e-05,
                        3.0469636469663725e-05,
                        0.00012053781020124761,
                    ),
                    (
                        4.3485101647524803e-05,
                        2.4382332471219346e-05,
                        2.4046766615879637e-05,
                        9.5172958985582652e-05,
                    ),
                    (
                        3.5184659036602886e-05,
                        1.9411697131159546e-05,
                        1.914748805113672e-05,
                        7.5808745485372899e-05,
                    ),
                    (
                        2.8710573684483027e-05,
                        1.5573074040919529e-05,
                        1.5362006532190414e-05,
                        6.0836608331866169e-05,
                    ),
                    (
                        2.3596485999114507e-05,
                        1.2573222214622194e-05,
                        1.2402532788994784e-05,
                        4.9125028134176529e-05,
                    ),
                    (
                        1.9508685715836919e-05,
                        1.0202969147411138e-05,
                        1.0063526905956131e-05,
                        3.9864751628194759e-05,
                    ),
                    (
                        1.6204365827687939e-05,
                        8.310707531476441e-06,
                        8.1958593174412469e-06,
                        3.246793683964068e-05,
                    ),
                    (
                        1.3503921608329314e-05,
                        6.7848331921788916e-06,
                        6.6896542854973157e-06,
                        2.6501040497028117e-05,
                    ),
                    (
                        1.1272082960471908e-05,
                        5.5418275702450897e-06,
                        5.462627168620374e-06,
                        2.1639107879962109e-05,
                    ),
                    (
                        9.4045207292483327e-06,
                        4.5178819909572341e-06,
                        4.4518747990969652e-06,
                        1.7633538292126709e-05,
                    ),
                    (
                        7.8173435227433508e-06,
                        3.6625163533856096e-06,
                        3.6076109486006602e-06,
                        1.4287438740418093e-05,
                    ),
                    (
                        6.4359614189096637e-06,
                        2.93230577098992e-06,
                        2.8869969639588791e-06,
                        1.1431285631297126e-05,
                    ),
                    (
                        5.1694249624199883e-06,
                        2.2779077252622501e-06,
                        2.2413635447347678e-06,
                        8.8723608152712553e-06,
                    ),
                ),
                order="F",
            ).transpose(),
            expected_nu_star_averaged=np.array(
                np.array(
                    (
                        0.030792808467814459,
                        0.01661476041777428,
                        0.014667916022065946,
                        0.053033247562414364,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
            expected_vd=np.array(
                (
                    (
                        3.1645071193768142,
                        3.0062817634079737,
                        3.0062817634079737,
                        1.5031408817039869,
                    ),
                    (
                        16.682831690173259,
                        15.848690105664597,
                        15.848690105664597,
                        7.9243450528322983,
                    ),
                    (
                        41.04114232193708,
                        38.989085205840233,
                        38.989085205840233,
                        19.494542602920117,
                    ),
                    (
                        76.310197904479367,
                        72.494688009255412,
                        72.494688009255412,
                        36.247344004627706,
                    ),
                    (
                        122.58635603762467,
                        116.45703823574345,
                        116.45703823574345,
                        58.228519117871727,
                    ),
                    (
                        179.99726143272298,
                        170.99739836108685,
                        170.99739836108685,
                        85.498699180543426,
                    ),
                    (
                        248.70423427095045,
                        236.26902255740296,
                        236.26902255740296,
                        118.13451127870148,
                    ),
                    (
                        328.90494695992078,
                        312.45969961192475,
                        312.45969961192475,
                        156.22984980596237,
                    ),
                    (
                        420.8367644782997,
                        399.79492625438479,
                        399.79492625438479,
                        199.8974631271924,
                    ),
                    (
                        524.78093220104745,
                        498.54188559099515,
                        498.54188559099515,
                        249.27094279549758,
                    ),
                    (
                        641.06781150569088,
                        609.01442093040646,
                        609.01442093040646,
                        304.50721046520323,
                    ),
                    (
                        770.08342250666237,
                        731.57925138132941,
                        731.57925138132941,
                        365.78962569066471,
                    ),
                    (
                        912.27764204794528,
                        866.66375994554812,
                        866.66375994554812,
                        433.33187997277406,
                    ),
                    (
                        1068.17453180507,
                        1014.7658052148167,
                        1014.7658052148167,
                        507.38290260740837,
                    ),
                    (
                        1238.3854536706519,
                        1176.4661809871195,
                        1176.4661809871195,
                        588.23309049355976,
                    ),
                    (
                        1423.6258968110376,
                        1352.4446019704858,
                        1352.4446019704858,
                        676.2223009852429,
                    ),
                    (
                        1624.7373411386941,
                        1543.5004740817596,
                        1543.5004740817596,
                        771.75023704087982,
                    ),
                    (
                        1842.7160968488222,
                        1750.5802920063811,
                        1750.5802920063811,
                        875.29014600319056,
                    ),
                    (
                        2078.7520308138296,
                        1974.8144292731388,
                        1974.8144292731388,
                        987.40721463656939,
                    ),
                    (
                        2334.2816737853227,
                        2217.5675900960568,
                        2217.5675900960568,
                        1108.7837950480284,
                    ),
                    (
                        2611.0628789014031,
                        2480.5097349563334,
                        2480.5097349563334,
                        1240.2548674781667,
                    ),
                    (
                        2911.2829228925325,
                        2765.7187767479063,
                        2765.7187767479063,
                        1382.8593883739532,
                    ),
                    (
                        3237.7206951043381,
                        3075.8346603491214,
                        3075.8346603491214,
                        1537.9173301745607,
                    ),
                    (
                        3594.0008558270988,
                        3414.3008130357439,
                        3414.3008130357439,
                        1707.1504065178719,
                    ),
                    (
                        3985.0143794458554,
                        3785.763660473563,
                        3785.763660473563,
                        1892.8818302367815,
                    ),
                    (
                        4417.6648724129982,
                        4196.7816287923488,
                        4196.7816287923488,
                        2098.3908143961744,
                    ),
                    (
                        4902.3232291036566,
                        4657.2070676484745,
                        4657.2070676484745,
                        2328.6035338242373,
                    ),
                    (
                        5456.0662712488684,
                        5183.2629576864265,
                        5183.2629576864265,
                        2591.6314788432132,
                    ),
                    (
                        6111.5444291433778,
                        5805.9672076862098,
                        5805.9672076862098,
                        2902.9836038431049,
                    ),
                    (
                        6952.6857290119615,
                        6605.0514425613637,
                        6605.0514425613637,
                        3302.5257212806819,
                    ),
                ),
                order="F",
            ).transpose(),
            expected_q_flux=np.array(
                np.array(
                    (
                        13499.211058929142,
                        472549.72072198224,
                        598376.09017838724,
                        14996.485377861178,
                    ),
                    order="F",
                ),
                order="F",
            ).transpose(),
        ),
    ),
)
def test_init_neoclassics(initneoclassicsparam, monkeypatch, neoclassics):
    """
    Automatically generated Regression Unit Test for init_neoclassics.

    This test was generated using data from stellarator/stellarator.IN.DAT.

    :param initneoclassicsparam: the data used to mock and assert in this test.
    :type initneoclassicsparam: initneoclassicsparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """
    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_electron_on_axis",
        initneoclassicsparam.nd_plasma_electron_on_axis,
    )
    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_electron_on_axis_kev",
        initneoclassicsparam.temp_plasma_electron_on_axis_kev,
    )
    monkeypatch.setattr(physics_variables, "alphan", initneoclassicsparam.alphan)
    monkeypatch.setattr(physics_variables, "alphat", initneoclassicsparam.alphat)
    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_ion_on_axis_kev",
        initneoclassicsparam.temp_plasma_ion_on_axis_kev,
    )
    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_ions_on_axis",
        initneoclassicsparam.nd_plasma_ions_on_axis,
    )
    monkeypatch.setattr(
        physics_variables,
        "f_plasma_fuel_deuterium",
        initneoclassicsparam.f_plasma_fuel_deuterium,
    )
    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_alphas_vol_avg",
        initneoclassicsparam.nd_plasma_alphas_vol_avg,
    )
    monkeypatch.setattr(physics_variables, "rminor", initneoclassicsparam.rminor)
    monkeypatch.setattr(physics_variables, "rmajor", initneoclassicsparam.rmajor)
    monkeypatch.setattr(
        physics_variables,
        "b_plasma_toroidal_on_axis",
        initneoclassicsparam.b_plasma_toroidal_on_axis,
    )
    monkeypatch.setattr(
        physics_variables, "temp_plasma_electron_vol_avg_kev", initneoclassicsparam.te
    )
    monkeypatch.setattr(
        physics_variables, "temp_plasma_ion_vol_avg_kev", initneoclassicsparam.ti
    )
    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_electrons_vol_avg",
        initneoclassicsparam.nd_plasma_electrons_vol_avg,
    )
    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_fuel_ions_vol_avg",
        initneoclassicsparam.nd_plasma_fuel_ions_vol_avg,
    )
    monkeypatch.setattr(
        neoclassics_variables, "densities", initneoclassicsparam.densities
    )
    monkeypatch.setattr(
        neoclassics_variables, "temperatures", initneoclassicsparam.temperatures
    )
    monkeypatch.setattr(
        neoclassics_variables, "dr_densities", initneoclassicsparam.dr_densities
    )
    monkeypatch.setattr(
        neoclassics_variables, "dr_temperatures", initneoclassicsparam.dr_temperatures
    )
    monkeypatch.setattr(neoclassics_variables, "roots", initneoclassicsparam.roots)
    monkeypatch.setattr(neoclassics_variables, "weights", initneoclassicsparam.weights)
    monkeypatch.setattr(neoclassics_variables, "nu", initneoclassicsparam.nu)
    monkeypatch.setattr(neoclassics_variables, "nu_star", initneoclassicsparam.nu_star)
    monkeypatch.setattr(
        neoclassics_variables, "nu_star_averaged", initneoclassicsparam.nu_star_averaged
    )
    monkeypatch.setattr(neoclassics_variables, "vd", initneoclassicsparam.vd)
    monkeypatch.setattr(neoclassics_variables, "iota", initneoclassicsparam.iota)
    monkeypatch.setattr(neoclassics_variables, "q_flux", initneoclassicsparam.q_flux)
    monkeypatch.setattr(neoclassics_variables, "eps_eff", initneoclassicsparam.eps_eff)
    monkeypatch.setattr(neoclassics_variables, "r_eff", initneoclassicsparam.r_eff)

    neoclassics.init_neoclassics(
        r_effin=initneoclassicsparam.r_effin,
        eps_effin=initneoclassicsparam.eps_effin,
        iotain=initneoclassicsparam.iotain,
    )

    assert neoclassics_variables.densities == pytest.approx(
        initneoclassicsparam.expected_densities, rel=0.001
    )
    assert neoclassics_variables.temperatures == pytest.approx(
        initneoclassicsparam.expected_temperatures, rel=0.001
    )
    assert neoclassics_variables.dr_densities == pytest.approx(
        initneoclassicsparam.expected_dr_densities, rel=0.001
    )
    assert neoclassics_variables.dr_temperatures == pytest.approx(
        initneoclassicsparam.expected_dr_temperatures, rel=0.001
    )
    assert neoclassics_variables.roots == pytest.approx(
        initneoclassicsparam.expected_roots, rel=0.001
    )
    assert neoclassics_variables.weights == pytest.approx(
        initneoclassicsparam.expected_weights, rel=0.001
    )
    assert neoclassics_variables.nu == pytest.approx(
        initneoclassicsparam.expected_nu, rel=0.001
    )
    assert neoclassics_variables.nu_star == pytest.approx(
        initneoclassicsparam.expected_nu_star, rel=0.001
    )
    assert neoclassics_variables.nu_star_averaged == pytest.approx(
        initneoclassicsparam.expected_nu_star_averaged, rel=0.001
    )
    assert neoclassics_variables.vd == pytest.approx(
        initneoclassicsparam.expected_vd, rel=0.001
    )
    assert neoclassics_variables.q_flux == pytest.approx(
        initneoclassicsparam.expected_q_flux, rel=0.001
    )
