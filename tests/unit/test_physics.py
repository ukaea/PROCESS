"""Unit tests for physics.f90."""
from process.fortran import physics_module as pm
from pytest import approx


def test_diamagnetic_fraction_hender():
    """Test diamagnetic_fraction_hender()."""
    beta = 0.14
    diacf = pm.diamagnetic_fraction_hender(beta)
    assert diacf == approx(0.05, abs=0.0001)


def test_diamagnetic_fraction_scene():
    """Test diamagnetic_fraction_scene."""
    beta = 0.15
    q95 = 3.0
    q0 = 1.0
    diacf = pm.diamagnetic_fraction_scene(beta, q95, q0)
    assert diacf == approx(0.0460, abs=0.0001)


def test_ps_fraction_scene():
    """Test ps_fraction_scene."""
    beta = 0.15
    pscf = pm.ps_fraction_scene(beta)
    assert pscf == approx(-0.0135, abs=0.0001)
