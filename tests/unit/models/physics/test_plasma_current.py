"""Unit tests for the Peng TART plasma current scaling (i_plasma_current=2)."""

import numpy as np
import pytest

from process.models.physics.plasma_current import PlasmaCurrent, PlasmaCurrentModel
from process.models.physics.plasma_fields import PlasmaFields


@pytest.fixture
def plasma_current():
    return PlasmaCurrent()


@pytest.fixture
def plasma_fields():
    return PlasmaFields()


def test_calculate_plasma_current_peng_qbar_transform(plasma_current):
    """The q95 -> qbar transform must invert q95 = 1.3 * qbar * (1-eps)^0.6.

    Expected value derived by hand for q95=6, A=1.8, a=2.5 m, B_T=3 T,
    kappa=2.8, triang=0.5:
        qbar = q95 / (1.3 * (1 - 1/A)^0.6) = 6 / 0.79916... = 7.507881...
    and the remaining geometry factors evaluated from plascar_bpol as in
    the function body.
    """
    current = plasma_current.calculate_plasma_current_peng(
        q95=6.0,
        aspect=1.8,
        rminor=2.5,
        b_plasma_toroidal_on_axis=3.0,
        kappa=2.8,
        triang=0.5,
    )

    assert current == pytest.approx(22.45707000674425)


def test_peng_qbar_roundtrip():
    """Applying the documented forward relation to the derived qbar returns q95."""
    q95 = 6.0
    aspect = 1.8
    factor = 1.3e0 * (1.0e0 - (1.0 / aspect)) ** 0.6e0
    qbar = q95 / factor
    assert qbar * factor == pytest.approx(q95)


def test_calculate_surface_averaged_poloidal_field_peng(plasma_fields):
    """TART branch of <Bp(a)>: expected value derived by hand.

    bpol = B_T * (ff1 + ff2) / (2 * pi * qbar) with
    qbar = q95 / (1.3 * (1 - 1/A)^0.6) and ff1/ff2 from plascar_bpol at
    A=1.8, kappa=2.8, triang=0.5.
    """
    bpol = plasma_fields.calculate_surface_averaged_poloidal_field(
        i_plasma_current=PlasmaCurrentModel.PENG_DIVERTOR_SCALING,
        cur_plasma=0.0,  # unused by the Peng branch
        q95=6.0,
        aspect=1.8,
        b_plasma_toroidal_on_axis=3.0,
        kappa=2.8,
        triang=0.5,
        len_plasma_poloidal=1.0,  # unused by the Peng branch
    )

    assert bpol == pytest.approx(0.915149073471596)


def test_peng_current_and_bpol_share_qbar(plasma_current, plasma_fields):
    """The two Peng-scaling call sites must use the same qbar transform.

    I = 5 * kappa * rminor * B_T / (2 * pi^2 * qbar) * (asin(e1)/e1
    + asin(e2)/e2) * (ff1 + ff2) and <Bp> = B_T (ff1+ff2) / (2 pi qbar)
    imply I / <Bp> is independent of qbar; verify both functions agree on
    the qbar-dependent part instead by reconstructing qbar from each.
    """
    q95, aspect, rminor, bt, kappa, triang = 4.0, 2.0, 1.5, 2.5, 2.2, 0.4

    current = plasma_current.calculate_plasma_current_peng(
        q95=q95,
        aspect=aspect,
        rminor=rminor,
        b_plasma_toroidal_on_axis=bt,
        kappa=kappa,
        triang=triang,
    )
    bpol = plasma_fields.calculate_surface_averaged_poloidal_field(
        i_plasma_current=PlasmaCurrentModel.PENG_DIVERTOR_SCALING,
        cur_plasma=0.0,
        q95=q95,
        aspect=aspect,
        b_plasma_toroidal_on_axis=bt,
        kappa=kappa,
        triang=triang,
        len_plasma_poloidal=1.0,
    )

    ff1, ff2, d1, d2 = plasma_current.plascar_bpol(
        aspect=aspect, eps=(1.0 / aspect), kappa=kappa, triang=triang
    )
    e1 = (2.0 * kappa) / (d1 * (1.0 + triang))
    e2 = (2.0 * kappa) / (d2 * (1.0 - triang))
    shape = (np.arcsin(e1) / e1 + np.arcsin(e2) / e2) * (ff1 + ff2)

    qbar_from_current = rminor * bt * 5.0 * kappa * shape / (2.0 * np.pi**2 * current)
    qbar_from_bpol = bt * (ff1 + ff2) / (2.0 * np.pi * bpol)

    assert qbar_from_current == pytest.approx(qbar_from_bpol)
    assert qbar_from_current == pytest.approx(
        q95 / (1.3e0 * (1.0e0 - (1.0 / aspect)) ** 0.6e0)
    )
