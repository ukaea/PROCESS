from typing import Any, NamedTuple

import pytest

from process.data_structure import physics_variables as pv
from process.data_structure import tfcoil_variables as tfv
from process.data_structure import times_variables as tv
from process.data_structure import vacuum_variables as vacv
from process.vacuum import Vacuum, VacuumVessel


@pytest.fixture
def vacuum():
    """Provides Vacuum object for testing.

    :return vacuum: initialised Vacuum object
    :type vacuum: process.vacuum.Vacuum
    """
    return Vacuum()


@pytest.fixture
def vacuum_vessel():
    """Provides Vacuum object for testing.

    :return vacuum: initialised Vacuum object
    :type vacuum: process.vacuum.Vacuum
    """
    return VacuumVessel()


class TestVacuum:
    def test_simple_model(self, monkeypatch, vacuum):
        """Tests `vacuum_simple` subroutine.

        Values taken from first calling of the model in i_vacuum_pumping regression test.

        :param monkeypatch: Mock fixture
        :type monkeypatch: object

        :param tfcoil: fixture containing an initialised `TFCoil` object
        :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
        """
        monkeypatch.setattr(
            pv, "molflow_plasma_fuelling_required", 7.5745668997694112e22
        )
        monkeypatch.setattr(pv, "a_plasma_surface", 1500.3146527709359)
        monkeypatch.setattr(tfv, "n_tf_coils", 18)
        monkeypatch.setattr(tv, "t_plant_pulse_dwell", 500)
        monkeypatch.setattr(vacv, "outgasfactor", 0.0235)
        monkeypatch.setattr(vacv, "outgasindex", 1)
        monkeypatch.setattr(vacv, "pres_vv_chamber_base", 0.0005)
        monkeypatch.setattr(vacv, "f_a_vac_pump_port_plasma_surface", 0.0203)
        monkeypatch.setattr(vacv, "f_volflow_vac_pumps_impedance", 0.4)
        monkeypatch.setattr(vacv, "volflow_vac_pumps_max", 27.3)
        monkeypatch.setattr(vacv, "molflow_vac_pumps", 1.2155e22)

        n_iter_vacuum_pumps = vacuum.vacuum_simple(output=False)

        assert n_iter_vacuum_pumps == pytest.approx(14.082585474801862)

    def test_old_model(self, monkeypatch, vacuum):
        """Test `vacuum` subroutine.

        Values taken from first calling of the model in G-L_Nb-Ti regression test.
        """
        monkeypatch.setattr(pv, "p_fusion_total_mw", 2115.3899563651776)
        monkeypatch.setattr(pv, "temp_plasma_electron_vol_avg_kev", 15.872999999999999)
        monkeypatch.setattr(tv, "t_plant_pulse_coil_precharge", 30)
        monkeypatch.setattr(vacv, "i_vac_pump_dwell", 0)
        monkeypatch.setattr(vacv, "i_vacuum_pump_type", 1)
        monkeypatch.setattr(vacv, "pres_vv_chamber_base", 0.00050000000000000001)
        monkeypatch.setattr(vacv, "pres_div_chamber_burn", 0.35999999999999999)
        monkeypatch.setattr(vacv, "outgrat_fw", 1.3000000000000001e-08)
        monkeypatch.setattr(vacv, "temp_vv_chamber_gas_burn_end", 300)

        ndiv = 1
        pfusmw = 2115.3899563651776
        r0 = 8.1386000000000003
        aw = 3.2664151549205331
        dsol = 0.22500000000000003
        plasma_a_plasma_surface = 1468.3151179059994
        plasma_vol = 2907.2299918381777
        thshldo = 0.40000000000000002
        thshldi = 0.12000000000000001
        thtf = 0.63812000000000002
        ritf = 3.6371848450794664
        n_tf_coils = 18
        t_plant_pulse_dwell = 1800
        nplasma = 7.2834e19
        qtorus = 0
        gasld = 2.7947500651998464e-05
        nduct = 0
        pumpn = 0
        dlscalc = 0
        mvdsh = 0
        dimax = 0

        pumpn, nduct, dlscalc, mvdsh, dimax = vacuum.vacuum(
            pfusmw,
            r0,
            aw,
            dsol,
            plasma_a_plasma_surface,
            plasma_vol,
            thshldo,
            thshldi,
            thtf,
            ritf,
            n_tf_coils,
            t_plant_pulse_dwell,
            nplasma,
            ndiv,
            qtorus,
            gasld,
            output=False,
        )

        assert pumpn == 36.0
        assert nduct == 18
        assert dlscalc == pytest.approx(2.798765707267961)
        assert mvdsh == 0.0
        assert dimax == pytest.approx(0.42414752916950604)


class EllipticalVesselVolumes(NamedTuple):
    rmajor: Any = None
    rminor: Any = None
    triang: Any = None
    r_shld_inboard_inner: Any = None
    r_shld_outboard_outer: Any = None
    dz_vv_half: Any = None
    dr_vv_inboard: Any = None
    dr_vv_outboard: Any = None
    dz_vv_upper: Any = None
    dz_vv_lower: Any = None


@pytest.mark.parametrize(
    "elliptical_vessel_volumes, expected",
    [
        (
            EllipticalVesselVolumes(
                rmajor=8,
                rminor=2.6666666666666665,
                triang=0.5,
                r_shld_inboard_inner=4.083333333333334,
                r_shld_outboard_outer=12.716666666666667,
                dz_vv_half=7.5032752487304135,
                dr_vv_inboard=0.30000000000000004,
                dr_vv_outboard=0.30000000000000004,
                dz_vv_upper=0.30000000000000004,
                dz_vv_lower=0.30000000000000004,
            ),
            (
                pytest.approx(143.03162449152501),
                pytest.approx(441.04172325889158),
                pytest.approx(584.07334775041659),
            ),
        )
    ],
)
def test_elliptical_vessel_volumes(vacuum_vessel, elliptical_vessel_volumes, expected):
    """Tests `elliptical_vessel_volumes` function.

    :param elliptical_vessel_volumes: input parameters for the function
    :type elliptical_vessel_volumes: EllipticalVesselVolumes


    """

    vol_vv_inboard, vol_vv_outboard, vol_vv = (
        vacuum_vessel.calculate_elliptical_vessel_volumes(
            rmajor=elliptical_vessel_volumes.rmajor,
            rminor=elliptical_vessel_volumes.rminor,
            triang=elliptical_vessel_volumes.triang,
            r_shld_inboard_inner=elliptical_vessel_volumes.r_shld_inboard_inner,
            r_shld_outboard_outer=elliptical_vessel_volumes.r_shld_outboard_outer,
            dz_vv_half=elliptical_vessel_volumes.dz_vv_half,
            dr_vv_inboard=elliptical_vessel_volumes.dr_vv_inboard,
            dr_vv_outboard=elliptical_vessel_volumes.dr_vv_outboard,
            dz_vv_upper=elliptical_vessel_volumes.dz_vv_upper,
            dz_vv_lower=elliptical_vessel_volumes.dz_vv_lower,
        )
    )

    assert vol_vv_inboard == expected[0]
    assert vol_vv_outboard == expected[1]
    assert vol_vv == expected[2]


class DShapedVesselVolumes(NamedTuple):
    r_shld_inboard_inner: Any = None
    r_shld_outboard_outer: Any = None
    dz_vv_half: Any = None
    dr_vv_inboard: Any = None
    dr_vv_outboard: Any = None
    dz_vv_upper: Any = None
    dz_vv_lower: Any = None


@pytest.mark.parametrize(
    "dshaped_vessel_volumes, expected",
    [
        (
            DShapedVesselVolumes(
                r_shld_inboard_inner=1.5,
                r_shld_outboard_outer=8.4000000000000004,
                dz_vv_half=9.4349999999999987,
                dr_vv_inboard=0.20000000000000001,
                dr_vv_outboard=0.30000000000000004,
                dz_vv_upper=0.30000000000000004,
                dz_vv_lower=0.30000000000000004,
            ),
            (
                pytest.approx(34.253413020620215),
                pytest.approx(306.20028292282814),
                pytest.approx(340.45369594344834),
            ),
        )
    ],
)
def test_dshaped_vessel_volumes(vacuum_vessel, dshaped_vessel_volumes, expected):
    """Tests `dshaped_vessel_volumes` function.

    :param dshaped_vessel_volumes: input parameters for the function
    :type dshaped_vessel_volumes: DShapedVesselVolumes

    """

    vol_vv_inboard, vol_vv_outboard, vol_vv = (
        vacuum_vessel.calculate_dshaped_vessel_volumes(
            r_shld_inboard_inner=dshaped_vessel_volumes.r_shld_inboard_inner,
            r_shld_outboard_outer=dshaped_vessel_volumes.r_shld_outboard_outer,
            dz_vv_half=dshaped_vessel_volumes.dz_vv_half,
            dr_vv_inboard=dshaped_vessel_volumes.dr_vv_inboard,
            dr_vv_outboard=dshaped_vessel_volumes.dr_vv_outboard,
            dz_vv_upper=dshaped_vessel_volumes.dz_vv_upper,
            dz_vv_lower=dshaped_vessel_volumes.dz_vv_lower,
        )
    )

    assert vol_vv_inboard == expected[0]
    assert vol_vv_outboard == expected[1]
    assert vol_vv == expected[2]
