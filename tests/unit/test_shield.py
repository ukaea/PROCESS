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


class EllipticalShieldAreas(NamedTuple):
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
    "elliptical_shield_areas, expected",
    [
        (
            EllipticalShieldAreas(
                rmajor=8,
                rminor=2.6666666666666665,
                triang=0.5,
                rsldi=4.083333333333334,
                rsldo=12.716666666666667,
                dz_shield_half=6.8032752487304133,
                dr_shield_inboard=0.30000000000000004,
                dr_shield_outboard=0.80000000000000004,
                dz_shield_upper=0.59999999999999998,
            ),
            (
                pytest.approx(700.06731267447844),
                pytest.approx(1344.1106481995357),
                pytest.approx(2044.1779608740142),
            ),
        )
    ],
)
def test_elliptical_shield_areas(shield, elliptical_shield_areas, expected):
    """Tests `elliptical_shield_areas` function.
    :param elliptical_shield_areas: input parameters for the function
    :type elliptical_shield_areas: EllipticalShieldAreas


    """

    a_shield_inboard, a_shield_outboard, a_shield = (
        shield.calculate_elliptical_shield_areas(
            rsldi=elliptical_shield_areas.rsldi,
            rsldo=elliptical_shield_areas.rsldo,
            rmajor=elliptical_shield_areas.rmajor,
            triang=elliptical_shield_areas.triang,
            dr_shld_inboard=elliptical_shield_areas.dr_shield_inboard,
            rminor=elliptical_shield_areas.rminor,
            dz_shld_half=elliptical_shield_areas.dz_shield_half,
            dr_shld_outboard=elliptical_shield_areas.dr_shield_outboard,
        )
    )

    assert a_shield_inboard == expected[0]
    assert a_shield_outboard == expected[1]
    assert a_shield == expected[2]


class DShapedShieldVolumes(NamedTuple):
    rminor: Any = None
    rsldi: Any = None
    dr_shld_inboard: Any = None
    dr_fw_plasma_gap_inboard: Any = None
    dr_fw_plasma_gap_outboard: Any = None
    dr_fw_inboard: Any = None
    dr_fw_outboard: Any = None
    dr_blkt_inboard: Any = None
    dr_blkt_outboard: Any = None
    dz_shld_half: Any = None
    dr_shld_outboard: Any = None
    dz_shld_upper: Any = None


@pytest.mark.parametrize(
    "dshaped_shield_volumes, expected",
    [
        (
            DShapedShieldVolumes(
                rminor=2.5,
                dr_shld_inboard=0.40000000000000002,
                rsldi=1.5,
                dr_fw_plasma_gap_inboard=0.10000000000000001,
                dr_fw_plasma_gap_outboard=0.10000000000000001,
                dr_fw_inboard=0.018000000000000002,
                dr_fw_outboard=0.018000000000000002,
                dr_blkt_inboard=0.0,
                dr_blkt_outboard=1.0,
                dz_shld_half=8.75,
                dr_shld_outboard=0.30000000000000004,
                dz_shld_upper=0.6000000000000009,
            ),
            (
                pytest.approx(79.896984366095609),
                pytest.approx(370.5642451119993),
                pytest.approx(450.46122947809488),
            ),
        )
    ],
)
def test_dshaped_shield_volumes(shield, dshaped_shield_volumes, expected):
    """Tests `dshaped_shield_volumes` function.

    :param dshaped_shield_volumes: input parameters for the function
    :type dshaped_shield_volumes: DShapedShieldVolumes


    """

    vol_shield_inboard, vol_shield_outboard, vol_shield = (
        shield.calculate_dshaped_shield_volumes(
            rsldi=dshaped_shield_volumes.rsldi,
            dr_shld_inboard=dshaped_shield_volumes.dr_shld_inboard,
            dr_fw_inboard=dshaped_shield_volumes.dr_fw_inboard,
            dr_fw_plasma_gap_inboard=dshaped_shield_volumes.dr_fw_plasma_gap_inboard,
            rminor=dshaped_shield_volumes.rminor,
            dr_fw_plasma_gap_outboard=dshaped_shield_volumes.dr_fw_plasma_gap_outboard,
            dr_fw_outboard=dshaped_shield_volumes.dr_fw_outboard,
            dr_blkt_inboard=dshaped_shield_volumes.dr_blkt_inboard,
            dr_blkt_outboard=dshaped_shield_volumes.dr_blkt_outboard,
            dz_shld_half=dshaped_shield_volumes.dz_shld_half,
            dr_shld_outboard=dshaped_shield_volumes.dr_shld_outboard,
            dz_shld_upper=dshaped_shield_volumes.dz_shld_upper,
        )
    )

    assert vol_shield_inboard == expected[0]
    assert vol_shield_outboard == expected[1]
    assert vol_shield == expected[2]
