import numpy as np
import pytest

from process.models.physics.scrape_off_layer import ScrapeOffLayer


@pytest.mark.parametrize(
    ("p_plasma_separatrix_mw", "rmajor", "b_plasma_surface_poloidal_average", "aspect"),
    [
        (100.0, 3.0, 0.5, 3.0),
        (10.0, 3.0, 0.5, 3.0),
        (500.0, 3.0, 0.5, 3.0),
        (100.0, 10.0, 0.5, 3.0),
        (100.0, 1.0, 0.5, 3.0),
    ],
)
def test_calculate_eich2013_sol_power_decay_length(
    p_plasma_separatrix_mw, rmajor, b_plasma_surface_poloidal_average, aspect
):
    """Test Eich 2013 SOL power decay length with various parameters."""
    result = ScrapeOffLayer.calculate_eich2013_sol_power_decay_length(
        p_plasma_separatrix_mw=p_plasma_separatrix_mw,
        rmajor=rmajor,
        b_plasma_surface_poloidal_average=b_plasma_surface_poloidal_average,
        aspect=aspect,
    )
    assert isinstance(result, float)
    assert result > 0


def test_calculate_eich2013_sol_power_decay_length_exact():
    """Test Eich 2013 SOL power decay length with exact value check."""
    result = ScrapeOffLayer.calculate_eich2013_sol_power_decay_length(
        p_plasma_separatrix_mw=100.0,
        rmajor=3.0,
        b_plasma_surface_poloidal_average=0.5,
        aspect=3.0,
    )
    assert isinstance(result, float)
    assert pytest.approx(result) == 0.0015345296622855315


@pytest.mark.parametrize(
    ("p_plasma_separatrix_mw", "b_plasma_surface_poloidal_average"),
    [
        (100.0, 0.5),
        (10.0, 0.5),
        (500.0, 0.5),
        (100.0, 2.0),
    ],
)
def test_calculate_mast2014_sol_power_decay_length_1(
    p_plasma_separatrix_mw, b_plasma_surface_poloidal_average
):
    """Test MAST 2014 SOL power decay length 1 with various parameters."""
    result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_1(
        p_plasma_separatrix_mw=p_plasma_separatrix_mw,
        b_plasma_surface_poloidal_average=b_plasma_surface_poloidal_average,
    )
    assert isinstance(result, float)
    assert result > 0


def test_calculate_mast2014_sol_power_decay_length_1_exact():
    """Test MAST 2014 SOL power decay length 1 with exact value check."""
    result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_1(
        p_plasma_separatrix_mw=100.0,
        b_plasma_surface_poloidal_average=0.5,
    )
    assert isinstance(result, float)
    assert pytest.approx(result) == 0.006753333858250275


@pytest.mark.parametrize(
    ("p_plasma_separatrix_mw", "cur_plasma_ma"),
    [
        (100.0, 1.0),
        (10.0, 1.0),
        (500.0, 1.0),
        (100.0, 3.0),
    ],
)
def test_calculate_mast2014_sol_power_decay_length_2(
    p_plasma_separatrix_mw, cur_plasma_ma
):
    """Test MAST 2014 SOL power decay length 2 with various parameters."""
    result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_2(
        p_plasma_separatrix_mw=p_plasma_separatrix_mw,
        cur_plasma_ma=cur_plasma_ma,
    )
    assert isinstance(result, float)
    assert result > 0


def test_calculate_mast2014_sol_power_decay_length_2_exact():
    """Test MAST 2014 SOL power decay length 2 with exact value check."""
    result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_2(
        p_plasma_separatrix_mw=100.0,
        cur_plasma_ma=1.0,
    )
    assert isinstance(result, float)
    assert pytest.approx(result) == 0.01258682517425542


def test_calculate_upstream_sol_outboard_parallel_area_exact():
    """Test upstream SOL outboard parallel area with exact value check."""
    result = ScrapeOffLayer.calculate_upstream_sol_outboard_parallel_area(
        rmajor=6.0,
        rminor=2.0,
        len_plasma_sol_power_decay=0.001,
        b_plasma_outboard_total=4.0,
        b_plasma_surface_poloidal_average=0.5,
    )
    assert isinstance(result, float)
    assert pytest.approx(result) == 0.006283185307179587


@pytest.mark.parametrize(
    "r",
    [
        8.001,
        8.01,
        8.1,
    ],
)
def test_calculate_outboard_midplane_near_sol_radial_profile(r):
    """Test outboard midplane near SOL radial profile with various parameters."""
    result = ScrapeOffLayer.calculate_outboard_midplane_near_sol_radial_profile(
        rmajor=6.0,
        rminor=2.0,
        len_plasma_sol_power_decay=0.001,
        pflux_plasma_outboard_sol_parallel_mw=10.0,
        r=r,
    )
    assert isinstance(result, float)
    assert result > 0


def test_calculate_outboard_midplane_near_sol_radial_profile_exact():
    """Test outboard midplane near SOL radial profile with exact value check."""
    result = ScrapeOffLayer.calculate_outboard_midplane_near_sol_radial_profile(
        rmajor=6.0,
        rminor=2.0,
        len_plasma_sol_power_decay=0.001,
        pflux_plasma_outboard_sol_parallel_mw=10.0,
        r=8.001,
    )
    assert isinstance(result, float)
    assert pytest.approx(result) == 3.678794411714423


def test_calculate_outboard_midplane_near_sol_radial_profile_array():
    """Test outboard midplane near SOL radial profile with array input."""
    r = np.array([8.001, 8.002, 8.003])
    result = ScrapeOffLayer.calculate_outboard_midplane_near_sol_radial_profile(
        rmajor=6.0,
        rminor=2.0,
        len_plasma_sol_power_decay=0.001,
        pflux_plasma_outboard_sol_parallel_mw=10.0,
        r=r,
    )
    assert isinstance(result, np.ndarray)
    assert np.all(result > 0)


def test_calculate_outboard_midplane_near_sol_radial_profile_invalid_r():
    """Test outboard midplane near SOL radial profile raises for r inside plasma edge."""
    with pytest.raises(ValueError, match=r"inside plasma edge|outside plasma"):
        ScrapeOffLayer.calculate_outboard_midplane_near_sol_radial_profile(
            rmajor=6.0,
            rminor=2.0,
            len_plasma_sol_power_decay=0.001,
            pflux_plasma_outboard_sol_parallel_mw=10.0,
            r=7.0,
        )


@pytest.mark.parametrize(
    "r",
    [
        8.001,
        8.01,
        8.1,
    ],
)
def test_calculate_eich_target_heat_flux_profile(r):
    """Test Eich target heat flux profile with various parameters."""
    result = ScrapeOffLayer.calculate_eich_target_heat_flux_profile(
        rmajor=6.0,
        rminor=2.0,
        pflux_plasma_sol_parallel_mw=10.0,
        len_plasma_sol_power_decay=0.001,
        f_b_div_flux_expansion=2.0,
        len_plasma_sol_power_spreading=0.001,
        pflux_target_background_heat_flux_mw=0.01,
        r=r,
    )
    assert isinstance(result, float)
    assert result > 0


def test_calculate_eich_target_heat_flux_profile_exact():
    """Test Eich target heat flux profile with exact value check."""
    result = ScrapeOffLayer.calculate_eich_target_heat_flux_profile(
        rmajor=6.0,
        rminor=2.0,
        pflux_plasma_sol_parallel_mw=10.0,
        len_plasma_sol_power_decay=0.001,
        f_b_div_flux_expansion=2.0,
        len_plasma_sol_power_spreading=0.001,
        pflux_target_background_heat_flux_mw=0.01,
        r=8.001,
    )
    assert isinstance(result, float)
    assert pytest.approx(result) == 3.8999590240461988
