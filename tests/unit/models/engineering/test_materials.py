import pytest

from process.models.engineering.materials import eurofer97_thermal_conductivity


def test_eurofer97_thermal_conductivity(monkeypatch, process_models):
    monkeypatch.setattr(process_models.data.fwbs, "fw_th_conductivity", 28.9)

    assert eurofer97_thermal_conductivity(1900.0, process_models.data) == pytest.approx(
        326.70406785462256
    )
