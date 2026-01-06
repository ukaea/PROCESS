from typing import Any, NamedTuple

import numpy as np
import pytest

from process.data_structure import divertor_variables, physics_variables
from process.plasma_profiles import PlasmaProfile
from process.profiles import NeProfile, TeProfile


class ProfileParam(NamedTuple):
    rho: int
    expected_profile: np.array


class NeProfileParam(NamedTuple):
    nd_plasma_separatrix_electron: float = 0.0
    i_plasma_pedestal: float = 0.0
    nd_plasma_electron_on_axis: float = 0.0
    nd_plasma_pedestal_electron: float = 0.0
    radius_plasma_pedestal_density_norm: float = 0.0
    alphan: float = 0.0
    expected_neprofile: list | None = None


@pytest.mark.parametrize(
    "neprofileparam",
    (
        NeProfileParam(
            nd_plasma_separatrix_electron=3.6421334486704804e19,
            i_plasma_pedestal=1,
            nd_plasma_electron_on_axis=0.0,
            nd_plasma_pedestal_electron=6.1916268627398164e19,
            radius_plasma_pedestal_density_norm=0.94000000000000006,
            alphan=1,
            expected_neprofile=[
                1.12500000e20,
                1.12275191e20,
                1.11587869e20,
                1.10396307e20,
                1.08619112e20,
                1.06109483e20,
                1.02592152e20,
                9.74783439e19,
                8.90714603e19,
                4.00000000e19,
            ],
        ),
    ),
    ids=["baseline_2018"],
)
def test_neprofile(neprofileparam: ProfileParam, monkeypatch):
    neprofile = NeProfile(10)
    neprofile.run()
    assert neprofile.profile_y == pytest.approx(neprofileparam.expected_neprofile)


def test_ncore():
    neprofile = NeProfile(10)
    radius_plasma_pedestal_density_norm = 0.94
    nped = 5.8300851381352219e19
    nsep = 3.4294618459618943e19
    nav = 7.4321e19
    alphan = 1.0
    assert neprofile.ncore(
        radius_plasma_pedestal_density_norm, nped, nsep, nav, alphan
    ) == pytest.approx(9.7756974320342041e19)


class TeProfileParam(NamedTuple):
    radius_plasma_pedestal_temp_norm: float = 0.0
    tbeta: float = 0.0
    temp_plasma_separatrix_kev: float = 0.0
    i_plasma_pedestal: float = 0.0
    alphat: float = 0.0
    temp_plasma_pedestal_kev: float = 0.0
    expected_teprofile: Any = np.array


@pytest.mark.parametrize(
    "teprofileparam",
    (
        TeProfileParam(
            radius_plasma_pedestal_temp_norm=0.94000000000000006,
            tbeta=2,
            temp_plasma_separatrix_kev=0.10000000000000001,
            i_plasma_pedestal=1,
            alphat=1.45,
            temp_plasma_pedestal_kev=5.5,
            expected_teprofile=[
                18.85,
                18.739472621498333,
                18.403679368327712,
                17.829141392239833,
                16.990144534125456,
                15.841907634203302,
                14.30460446612375,
                12.219427594826554,
                9.177492824141694,
                1.0,
            ],
        ),
    ),
    ids=["baseline_2018"],
)
def test_teprofile(teprofileparam: ProfileParam, monkeypatch):
    monkeypatch.setattr(
        physics_variables, "i_plasma_pedestal", teprofileparam.i_plasma_pedestal
    )
    teprofile = TeProfile(10)
    teprofile.run()
    assert teprofile.profile_y == pytest.approx(teprofileparam.expected_teprofile)


def test_tcore():
    teprofile = TeProfile(10)
    radius_plasma_pedestal_temp_norm = 0.94
    tped = 3.7775374842470044
    tsep = 0.1
    tav = 12.33
    alphat = 1.45
    tbeta = 2.0

    assert teprofile.tcore(
        radius_plasma_pedestal_temp_norm, tped, tsep, tav, alphat, tbeta
    ) == pytest.approx(28.09093632260765)


class PlasmaProfilesParam(NamedTuple):
    prn1: float = 0.0

    radius_plasma_pedestal_temp_norm: float = 0.0

    temp_plasma_electron_density_weighted_kev: float = 0.0

    temp_plasma_ion_density_weighted_kev: float = 0.0

    alphap: float = 0.0

    tbeta: float = 0.0

    temp_plasma_electron_on_axis_kev: float = 0.0

    pres_plasma_thermal_on_axis: float = 0.0

    nd_plasma_separatrix_electron: float = 0.0

    temp_plasma_separatrix_kev: float = 0.0

    f_temp_plasma_electron_density_vol_avg: float = 0.0

    i_plasma_pedestal: float = 0.0

    nd_plasma_ions_on_axis: float = 0.0

    nd_plasma_electron_on_axis: float = 0.0

    temp_plasma_ion_on_axis_kev: float = 0.0

    f_temp_plasma_ion_electron: float = 0.0

    nd_plasma_electron_line: float = 0.0

    alphat: float = 0.0

    nd_plasma_ions_total_vol_avg: float = 0.0

    nd_plasma_pedestal_electron: float = 0.0

    temp_plasma_ion_vol_avg_kev: float = 0.0

    radius_plasma_pedestal_density_norm: float = 0.0

    nd_plasma_electrons_vol_avg: float = 0.0

    temp_plasma_pedestal_kev: float = 0.0

    alphan: float = 0.0

    te: float = 0.0

    rho_ne_max: float = 0.0

    rho_te_max: float = 0.0

    gradient_length_ne: float = 0.0

    gradient_length_te: float = 0.0

    rminor: float = 0.0

    expected_prn1: float = 0.0

    expected_ten: float = 0.0

    expected_tin: float = 0.0

    expected_alphap: float = 0.0

    expected_te0: float = 0.0

    expected_p0: float = 0.0

    expected_pcoef: float = 0.0

    expected_ni0: float = 0.0

    expected_ne0: float = 0.0

    expected_ti0: float = 0.0

    expected_nd_electron_line: float = 0.0

    expected_ti: float = 0.0


@pytest.mark.parametrize(
    "plasmaprofilesparam",
    (
        PlasmaProfilesParam(
            prn1=0.40000000000000002,
            radius_plasma_pedestal_temp_norm=0.94000000000000006,
            temp_plasma_electron_density_weighted_kev=0.0,
            temp_plasma_ion_density_weighted_kev=0.0,
            alphap=0.0,
            tbeta=2,
            temp_plasma_electron_on_axis_kev=0.0,
            pres_plasma_thermal_on_axis=0.0,
            nd_plasma_separatrix_electron=3.6421334486704804e19,
            temp_plasma_separatrix_kev=0.10000000000000001,
            f_temp_plasma_electron_density_vol_avg=0.0,
            i_plasma_pedestal=1,
            nd_plasma_ions_on_axis=0.0,
            nd_plasma_electron_on_axis=0.0,
            temp_plasma_ion_on_axis_kev=0.0,
            f_temp_plasma_ion_electron=1,
            nd_plasma_electron_line=0.0,
            alphat=1.45,
            nd_plasma_ions_total_vol_avg=6.9461125748017857e19,
            nd_plasma_pedestal_electron=6.1916268627398164e19,
            temp_plasma_ion_vol_avg_kev=12.9,
            radius_plasma_pedestal_density_norm=0.94000000000000006,
            nd_plasma_electrons_vol_avg=7.983e19,
            temp_plasma_pedestal_kev=5.5,
            alphan=1,
            te=13.07,
            rho_ne_max=0.0,
            rho_te_max=0.0,
            gradient_length_ne=0.0,
            gradient_length_te=0.0,
            rminor=2.9264516129032256,
            expected_prn1=0.45623618297262686,
            expected_ten=14.52233022043558,
            expected_tin=14.52233022043558,
            expected_alphap=2.4500000000000002,
            expected_te0=27.370104119511087,
            expected_p0=868106.0658743214,
            expected_pcoef=1.111119374172577,
            expected_ni0=9.210720071916929e19,
            expected_ne0=1.0585658890823703e20,
            expected_ti0=27.370104119511087,
            expected_nd_electron_line=8.8687354645836431e19,
            expected_ti=13.07,
        ),
        PlasmaProfilesParam(
            prn1=0.45623618297262686,
            radius_plasma_pedestal_temp_norm=0.94000000000000006,
            temp_plasma_electron_density_weighted_kev=14.521871327399182,
            temp_plasma_ion_density_weighted_kev=14.521871327399182,
            alphap=2.4500000000000002,
            tbeta=2,
            temp_plasma_electron_on_axis_kev=27.369013322953624,
            pres_plasma_thermal_on_axis=868071.46874220832,
            nd_plasma_separatrix_electron=3.6421334486704804e19,
            temp_plasma_separatrix_kev=0.10000000000000001,
            f_temp_plasma_electron_density_vol_avg=1.1110842637642833,
            i_plasma_pedestal=1,
            nd_plasma_ions_on_axis=9.210720071916929e19,
            nd_plasma_electron_on_axis=1.0585658890823703e20,
            temp_plasma_ion_on_axis_kev=27.369013322953624,
            f_temp_plasma_ion_electron=1,
            nd_plasma_electron_line=8.8687354645836431e19,
            alphat=1.45,
            nd_plasma_ions_total_vol_avg=6.9461125748017857e19,
            nd_plasma_pedestal_electron=6.1916268627398164e19,
            temp_plasma_ion_vol_avg_kev=13.07,
            radius_plasma_pedestal_density_norm=0.94000000000000006,
            nd_plasma_electrons_vol_avg=7.983e19,
            temp_plasma_pedestal_kev=5.5,
            alphan=1,
            te=13.07,
            rho_ne_max=0.0,
            rho_te_max=0.0,
            gradient_length_ne=0.0,
            gradient_length_te=0.0,
            rminor=2.9264516129032256,
            expected_prn1=0.45623618297262686,
            expected_ten=14.52233022043558,
            expected_tin=14.52233022043558,
            expected_alphap=2.4500000000000002,
            expected_te0=27.370104119511087,
            expected_p0=868106.0658743214,
            expected_pcoef=1.111119374172577,
            expected_ni0=9.210720071916929e19,
            expected_ne0=1.0585658890823703e20,
            expected_ti0=27.370104119511087,
            expected_nd_electron_line=8.8687354645836431e19,
            expected_ti=13.07,
        ),
    ),
)
def test_plasma_profiles(plasmaprofilesparam, monkeypatch):
    """
    Automatically generated Regression Unit Test for plasma_profiles.

    This test was generated using data from the HARE regression test, which
    has since been removed in preparation for open-sourcing (#1889).

    :param plasmaprofilesparam: the data used to mock and assert in this test.
    :type plasmaprofilesparam: plasmaprofilesparam

    :param monkeypatch: pytest fixture used to mock module/class variables
    :type monkeypatch: _pytest.monkeypatch.monkeypatch
    """

    monkeypatch.setattr(divertor_variables, "prn1", plasmaprofilesparam.prn1)

    monkeypatch.setattr(
        physics_variables,
        "radius_plasma_pedestal_temp_norm",
        plasmaprofilesparam.radius_plasma_pedestal_temp_norm,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_electron_density_weighted_kev",
        plasmaprofilesparam.temp_plasma_electron_density_weighted_kev,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_ion_density_weighted_kev",
        plasmaprofilesparam.temp_plasma_ion_density_weighted_kev,
    )

    monkeypatch.setattr(physics_variables, "alphap", plasmaprofilesparam.alphap)

    monkeypatch.setattr(physics_variables, "tbeta", plasmaprofilesparam.tbeta)

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_electron_on_axis_kev",
        plasmaprofilesparam.temp_plasma_electron_on_axis_kev,
    )

    monkeypatch.setattr(
        physics_variables,
        "pres_plasma_thermal_on_axis",
        plasmaprofilesparam.pres_plasma_thermal_on_axis,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_separatrix_electron",
        plasmaprofilesparam.nd_plasma_separatrix_electron,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_separatrix_kev",
        plasmaprofilesparam.temp_plasma_separatrix_kev,
    )

    monkeypatch.setattr(
        physics_variables,
        "f_temp_plasma_electron_density_vol_avg",
        plasmaprofilesparam.f_temp_plasma_electron_density_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables, "i_plasma_pedestal", plasmaprofilesparam.i_plasma_pedestal
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_ions_on_axis",
        plasmaprofilesparam.nd_plasma_ions_on_axis,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_electron_on_axis",
        plasmaprofilesparam.nd_plasma_electron_on_axis,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_ion_on_axis_kev",
        plasmaprofilesparam.temp_plasma_ion_on_axis_kev,
    )

    monkeypatch.setattr(
        physics_variables,
        "f_temp_plasma_ion_electron",
        plasmaprofilesparam.f_temp_plasma_ion_electron,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_electron_line",
        plasmaprofilesparam.nd_plasma_electron_line,
    )

    monkeypatch.setattr(physics_variables, "alphat", plasmaprofilesparam.alphat)

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_ions_total_vol_avg",
        plasmaprofilesparam.nd_plasma_ions_total_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_pedestal_electron",
        plasmaprofilesparam.nd_plasma_pedestal_electron,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_ion_vol_avg_kev",
        plasmaprofilesparam.temp_plasma_ion_vol_avg_kev,
    )

    monkeypatch.setattr(
        physics_variables,
        "radius_plasma_pedestal_density_norm",
        plasmaprofilesparam.radius_plasma_pedestal_density_norm,
    )

    monkeypatch.setattr(
        physics_variables,
        "nd_plasma_electrons_vol_avg",
        plasmaprofilesparam.nd_plasma_electrons_vol_avg,
    )

    monkeypatch.setattr(
        physics_variables,
        "temp_plasma_pedestal_kev",
        plasmaprofilesparam.temp_plasma_pedestal_kev,
    )

    monkeypatch.setattr(physics_variables, "alphan", plasmaprofilesparam.alphan)

    monkeypatch.setattr(
        physics_variables, "temp_plasma_electron_vol_avg_kev", plasmaprofilesparam.te
    )

    monkeypatch.setattr(physics_variables, "rho_ne_max", plasmaprofilesparam.rho_ne_max)

    monkeypatch.setattr(physics_variables, "rho_te_max", plasmaprofilesparam.rho_te_max)

    monkeypatch.setattr(
        physics_variables, "gradient_length_ne", plasmaprofilesparam.gradient_length_ne
    )

    monkeypatch.setattr(
        physics_variables, "gradient_length_te", plasmaprofilesparam.gradient_length_te
    )

    monkeypatch.setattr(physics_variables, "rminor", plasmaprofilesparam.rminor)

    monkeypatch.setattr(physics_variables, "a_plasma_poloidal", 4.0)

    plasmaprofile = PlasmaProfile()
    plasmaprofile.run()

    assert divertor_variables.prn1 == pytest.approx(plasmaprofilesparam.expected_prn1)

    assert physics_variables.temp_plasma_electron_density_weighted_kev == pytest.approx(
        plasmaprofilesparam.expected_ten
    )

    assert physics_variables.temp_plasma_ion_density_weighted_kev == pytest.approx(
        plasmaprofilesparam.expected_tin
    )

    assert physics_variables.alphap == pytest.approx(plasmaprofilesparam.expected_alphap)

    assert physics_variables.temp_plasma_electron_on_axis_kev == pytest.approx(
        plasmaprofilesparam.expected_te0
    )

    assert physics_variables.pres_plasma_thermal_on_axis == pytest.approx(
        plasmaprofilesparam.expected_p0
    )

    assert physics_variables.f_temp_plasma_electron_density_vol_avg == pytest.approx(
        plasmaprofilesparam.expected_pcoef
    )

    assert physics_variables.nd_plasma_ions_on_axis == pytest.approx(
        plasmaprofilesparam.expected_ni0
    )

    assert physics_variables.nd_plasma_electron_on_axis == pytest.approx(
        plasmaprofilesparam.expected_ne0
    )

    assert physics_variables.temp_plasma_ion_on_axis_kev == pytest.approx(
        plasmaprofilesparam.expected_ti0
    )

    assert physics_variables.nd_plasma_electron_line == pytest.approx(
        plasmaprofilesparam.expected_nd_electron_line
    )

    assert physics_variables.temp_plasma_ion_vol_avg_kev == pytest.approx(
        plasmaprofilesparam.expected_ti
    )
