import pytest

from process.models.engineering.pumping import (
    calculate_reynolds_number,
    darcy_friction_haaland,
    gnielinski_heat_transfer_coefficient,
)


def test_darcy_friction_haaland():
    assert darcy_friction_haaland(
        reynolds=5500, roughness_channel=1e-6, radius_channel=0.1
    ) == pytest.approx(0.0366668931278784)


def test_gnielinski_heat_transfer_coefficient():
    assert gnielinski_heat_transfer_coefficient(
        mflux_coolant=112.19853108876258,
        den_coolant=8.8673250601290707,
        radius_channel=0.0060000000000000001,
        heatcap_coolant=5184.9330299967578,
        visc_coolant=4.0416219836935569e-05,
        thermcond_coolant=0.3211653052986152,
        roughness_channel=6e-8,
    ) == pytest.approx(1929.2042015869506)


def test_calculate_reynolds_number():

    assert calculate_reynolds_number(
        den_coolant=8.8673250601290707,
        vel_coolant=12.649110769896881,
        radius_channel=0.0060000000000000001,
        visc_coolant=4.0416219836935569e-05,
    ) == pytest.approx(33302.602975971815)
