"""Unit tests for ife.f90."""
from pytest import approx
from process import fortran


def test_ifetgt(monkeypatch):
    """Test ifetgt.

    :param monkeypatch: Mock fixture
    :type monkeypatch: object
    """
    # Mock module variables
    # Repetition Rate (Hz)
    monkeypatch.setattr(fortran.ife_variables, "reprat", 4.0)
    # IFE target factory power at 6 Hz repetition rate
    monkeypatch.setattr(fortran.ife_variables, "ptargf", 2.0)
    monkeypatch.setattr(fortran.ife_variables, "tfacmw", 0)

    fortran.ife_module.ifetgt()
    assert fortran.ife_variables.tfacmw == approx(1.506, abs=0.001)
