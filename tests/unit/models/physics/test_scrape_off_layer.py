from process.models.physics.scrape_off_layer import ScrapeOffLayer


class TestScrapeOffLayer:
    """Test suite for ScrapeOffLayer class."""

    @staticmethod
    def test_calculate_eich2013_sol_power_decay_length_nominal():
        """Test Eich 2013 SOL power decay length with nominal values."""
        result = ScrapeOffLayer.calculate_eich2013_sol_power_decay_length(
            p_plasma_separatrix_mw=100.0,
            rmajor=3.0,
            b_plasma_surface_poloidal_average=0.5,
            aspect=3.0,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_eich2013_sol_power_decay_length_low_power():
        """Test with low plasma separatrix power."""
        result = ScrapeOffLayer.calculate_eich2013_sol_power_decay_length(
            p_plasma_separatrix_mw=10.0,
            rmajor=3.0,
            b_plasma_surface_poloidal_average=0.5,
            aspect=3.0,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_eich2013_sol_power_decay_length_high_power():
        """Test with high plasma separatrix power."""
        result = ScrapeOffLayer.calculate_eich2013_sol_power_decay_length(
            p_plasma_separatrix_mw=500.0,
            rmajor=3.0,
            b_plasma_surface_poloidal_average=0.5,
            aspect=3.0,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_eich2013_sol_power_decay_length_large_rmajor():
        """Test with large major radius."""
        result = ScrapeOffLayer.calculate_eich2013_sol_power_decay_length(
            p_plasma_separatrix_mw=100.0,
            rmajor=10.0,
            b_plasma_surface_poloidal_average=0.5,
            aspect=3.0,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_eich2013_sol_power_decay_length_small_rmajor():
        """Test with small major radius."""
        result = ScrapeOffLayer.calculate_eich2013_sol_power_decay_length(
            p_plasma_separatrix_mw=100.0,
            rmajor=1.0,
            b_plasma_surface_poloidal_average=0.5,
            aspect=3.0,
        )
        assert isinstance(result, float)
        assert result > 0
