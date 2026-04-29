import pytest

from process.models.engineering.pumping import darcy_friction_haaland


def test_darcy_friction_haaland():
    assert darcy_friction_haaland(
        reynolds=5500, roughness_channel=1e-6, radius_channel=0.1
    ) == pytest.approx(0.0366668931278784)
