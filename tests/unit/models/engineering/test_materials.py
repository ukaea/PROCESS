import pytest

from process.data_structure import fwbs_variables
from process.models.engineering.materials import eurofer97_thermal_conductivity


def test_eurofer97_thermal_conductivity(monkeypatch):
    monkeypatch.setattr(fwbs_variables, "fw_th_conductivity", 28.9)

    assert eurofer97_thermal_conductivity(1900.0) == pytest.approx(326.70406785462256)
