import pytest

from process.data_structure import physics_variables as pv
from process.data_structure import tfcoil_variables as tfv
from process.data_structure import times_variables as tv
from process.data_structure import vacuum_variables as vacv
from process.vacuum import Vacuum


@pytest.fixture
def vacuum():
    """Provides Vacuum object for testing.

    :return vacuum: initialised Vacuum object
    :type vacuum: process.vacuum.Vacuum
    """
    return Vacuum()


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
