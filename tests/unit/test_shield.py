from typing import Any, NamedTuple

import pytest

from process.shield import Shield


@pytest.fixture
def shield():
    """Provides Shield object for testing.

    :return shield: initialised Shield object
    :type shield: process.shield.Shield
    """
    return Shield()


class EllipticalShieldVolumes(NamedTuple):
    rmajor: Any = None
    rminor: Any = None
    triang: Any = None
    rsldi: Any = None
    rsldo: Any = None
    dz_shield_half: Any = None
    dr_shield_inboard: Any = None
    dr_shield_outboard: Any = None
    dz_shield_upper: Any = None
    dz_shield_lower: Any = None


@pytest.mark.parametrize(
    "elliptical_shield_volumes, expected",
    [
        (
            EllipticalShieldVolumes(
                rmajor=8,
                rminor=2.6666666666666665,
                triang=0.5,
                rsldi=4.083333333333334,
                rsldo=12.716666666666667,
                dz_shield_half=6.8032752487304133,
                dr_shield_inboard=0.30000000000000004,
                dr_shield_outboard=0.80000000000000004,
                dz_shield_upper=0.59999999999999998,
                dz_shield_lower=0.59999999999999998,
            ),
            (
                pytest.approx(177.89822933168091),
                pytest.approx(946.56393192782434),
                pytest.approx(1124.4621612595051),
            ),
        )
    ],
)
def test_elliptical_shield_volumes(shield, elliptical_shield_volumes, expected):
    """Tests `elliptical_shield_volumes` function.

    :param elliptical_shield_volumes: input parameters for the function
    :type elliptical_shield_volumes: EllipticalShieldVolumes


    """

    vol_shield_inboard, vol_shield_outboard, vol_shield = (
        shield.calculate_elliptical_shield_volumes(
            rsldi=elliptical_shield_volumes.rsldi,
            rsldo=elliptical_shield_volumes.rsldo,
            rmajor=elliptical_shield_volumes.rmajor,
            triang=elliptical_shield_volumes.triang,
            dr_shld_inboard=elliptical_shield_volumes.dr_shield_inboard,
            rminor=elliptical_shield_volumes.rminor,
            dz_shld_half=elliptical_shield_volumes.dz_shield_half,
            dr_shld_outboard=elliptical_shield_volumes.dr_shield_outboard,
            dz_shld_upper=elliptical_shield_volumes.dz_shield_upper,
        )
    )

    assert vol_shield_inboard == expected[0]
    assert vol_shield_outboard == expected[1]
    assert vol_shield == expected[2]
