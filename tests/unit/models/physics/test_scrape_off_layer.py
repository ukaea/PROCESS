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

    @staticmethod
    def test_calculate_mast2014_sol_power_decay_length_1_nominal():
        """Test MAST 2014 SOL power decay length 1 with nominal values."""
        result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_1(
            p_plasma_separatrix_mw=100.0,
            b_plasma_surface_poloidal_average=0.5,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_mast2014_sol_power_decay_length_1_low_power():
        """Test MAST 2014 SOL power decay length 1 with low power."""
        result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_1(
            p_plasma_separatrix_mw=10.0,
            b_plasma_surface_poloidal_average=0.5,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_mast2014_sol_power_decay_length_1_high_power():
        """Test MAST 2014 SOL power decay length 1 with high power."""
        result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_1(
            p_plasma_separatrix_mw=500.0,
            b_plasma_surface_poloidal_average=0.5,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_mast2014_sol_power_decay_length_1_high_bpol():
        """Test MAST 2014 SOL power decay length 1 with high poloidal field."""
        result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_1(
            p_plasma_separatrix_mw=100.0,
            b_plasma_surface_poloidal_average=2.0,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_mast2014_sol_power_decay_length_2_nominal():
        """Test MAST 2014 SOL power decay length 2 with nominal values."""
        result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_2(
            p_plasma_separatrix_mw=100.0,
            cur_plasma_ma=1.0,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_mast2014_sol_power_decay_length_2_low_power():
        """Test MAST 2014 SOL power decay length 2 with low power."""
        result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_2(
            p_plasma_separatrix_mw=10.0,
            cur_plasma_ma=1.0,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_mast2014_sol_power_decay_length_2_high_power():
        """Test MAST 2014 SOL power decay length 2 with high power."""
        result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_2(
            p_plasma_separatrix_mw=500.0,
            cur_plasma_ma=1.0,
        )
        assert isinstance(result, float)
        assert result > 0

    @staticmethod
    def test_calculate_mast2014_sol_power_decay_length_2_high_current():
        """Test MAST 2014 SOL power decay length 2 with high plasma current."""
        result = ScrapeOffLayer.calculate_mast2014_sol_power_decay_length_2(
            p_plasma_separatrix_mw=100.0,
            cur_plasma_ma=3.0,
        )
        assert isinstance(result, float)
        assert result > 0
