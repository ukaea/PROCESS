"""Unit and Integration tests for tfcoil.f90."""

from typing import NamedTuple
from process.sctfcoil import Sctfcoil
import pytest

from process.tfcoil import TFcoil

from process.fortran import tfcoil_variables as tfv
from process.fortran import build_variables as bv
from process.fortran import fwbs_variables as fwbsv
from process.build import Build


@pytest.fixture
def tfcoil():
    """Provides TFcoil object for testing.

    :param monkeypatch: pytest mocking fixture
    :type monkeypatch: MonkeyPatch

    :return tfcoil: initialised TFcoil object
    :type tfcoil: process.tfcoil.TFcoil
    """
    return TFcoil(build=Build(), sctfcoil=Sctfcoil())


@pytest.mark.parametrize(
    "temperature, expected_density",
    [(24.6, 130.02434313053487), (30.2, 113.09539723009078), (43.6, 85.26924709595201)],
)
def test_he_density(temperature, expected_density, tfcoil):
    """Tests `he_density` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param expected_density: expected result of the routine.
    :type expected_density: float

    :param tfcoil: fixture containing an initialised `TFcoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    density = tfcoil.he_density(temperature)

    assert pytest.approx(density) == expected_density


@pytest.mark.parametrize(
    "temperature, expected_cp",
    [(24.6, 5674.909063980127), (30.2, 5798.42049712345), (43.6, 5673.218322000001)],
)
def test_he_cp(temperature, expected_cp, tfcoil):
    """Tests `he_cp` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param expected_cp: expected result of the routine.
    :type expected_cp: float

    :param tfcoil: fixture containing an initialised `TFcoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """

    cp = tfcoil.he_cp(temperature)

    assert pytest.approx(cp) == expected_cp


@pytest.mark.parametrize(
    "temperature, expected_visco",
    [
        (20.6, 6.889108080243641e-06),
        (26.2, 6.859929884441028e-06),
        (43.6, 7.717393982e-06),
    ],
)
def test_he_visco(temperature, expected_visco, tfcoil):
    """Tests `he_visco` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param expected_visco: expected result of the routine.
    :type expected_visco: float

    :param tfcoil: fixture containing an initialised `TFcoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    visco = tfcoil.he_visco(temperature)

    assert pytest.approx(visco) == expected_visco


# error module needs to be initialised here because the temperature ranges are (in some cases) greater than the error conditions
@pytest.mark.parametrize(
    "temperature, expected_th_cond",
    [
        (20.6, 0.0585183573711527),
        (24.2, 0.05720100686027678),
        (43.6, 0.061189437089717184),
        (50.6, 0.06409264503),
        (54.4, 0.065706872),
    ],
)
def test_he_th_cond(temperature, expected_th_cond, reinitialise_error_module, tfcoil):
    """Tests `he_th_cond` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param expected_th_cond: expected result of the routine.
    :type expected_th_cond: float

    :param reinitialise_error_module: teardown any error side-effects
    :type reinitialise_error_module: None

    :param tfcoil: fixture containing an initialised `TFcoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    th_cond = tfcoil.he_th_cond(temperature)

    assert pytest.approx(th_cond) == expected_th_cond


@pytest.mark.parametrize(
    "temperature, expected_th_cond",
    [
        (54.4, 844.9049012800042),
        (66.9, 571.151543384937),
        (109.5, 233.66333020125035),
        (151, 250.4911087866094),
    ],
)
def test_al_th_cond(temperature, expected_th_cond, tfcoil):
    """Tests `he_th_cond` subroutine.

    :param temperature: test asset passed to the routine representing the temperature, in Kelvin.
    :type temperature: float

    :param al_th_cond: expected result of the routine.
    :type al_th_cond: float

    :param tfcoil: fixture containing an initialised `TFcoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    th_cond = tfcoil.al_th_cond(temperature)

    assert pytest.approx(th_cond) == expected_th_cond


class CntrpstTestAsset(NamedTuple):
    """Test asset for a test case of cntrpst

    :i_tf_sup: value for tfcoil_variables.i_tf_sup to be mocked with (0=Copper, 2=Cryogenic aluminium)
    :type i_tf_sup: integer
    :tcoolin: value for tfcoil_variables.tcoolin to be mocked with (centrepost coolant inlet temperature)
    :type tcoolin: float

    :expected_dtiocool: expected value of tfcoil_variables.dtiocool after tfcoil.cntrpst routine has run
    :type expected_dtiocool: float
    :expected_tcpav2: expected value of tfcoil_variables.tcpav2 after tfcoil.cntrpst routine has run
    :type expected_tcpav2: float
    :expected_tcpmax: expected value of tfcoil_variables.tcpmax after tfcoil.cntrpst routine has run
    :type expected_tcpmax: float
    :expected_ppump: expected value of tfcoil_variables.ppump after tfcoil.cntrpst routine has run
    :type expected_ppump: float
    """

    i_tf_sup: int
    tcoolin: float
    expected_dtiocool: float
    expected_tcpav2: float
    expected_tcpmax: float
    expected_ppump: float


@pytest.mark.parametrize(
    "cntrpst_asset",
    [
        CntrpstTestAsset(
            0, 100.0, 0.00075899, 100.00109611, 100.00147829, 7.05905966e08
        ),
        # CntrpstTestAsset(1, 45.0),
        CntrpstTestAsset(
            2, 43.6, 0.00645998, 43.60678774, 43.61001841, 80926408.5501315
        ),
    ],
)
def test_cntrpst(cntrpst_asset, monkeypatch, reinitialise_error_module, tfcoil):
    """Integration test for cntrpst

    Testing tfcoil module variables being set:
        - dtiocool
        - tcpav2
        - tcpmax
        - ppump

    :param cntrpst_asset: test asset containing values to mock and expected results for the represented test case
    :type cntrpst_asset: CntrpstTestAsset

    :param monkeypatch: Mock fixture
    :type monkeypatch: object

    :param reinitialise_error_module: teardown any error side-effects
    :type reinitialise_error_module: None

    :param tfcoil: fixture containing an initialised `TFcoil` object
    :type tfcoil: tests.unit.test_tfcoil.tfcoil (functional fixture)
    """
    monkeypatch.setattr(tfv, "a_cp_cool", 1)
    monkeypatch.setattr(tfv, "n_tf", 16)
    monkeypatch.setattr(tfv, "rcool", 0.005)
    monkeypatch.setattr(tfv, "vcool", 20.0)
    monkeypatch.setattr(tfv, "vol_cond_cp", 2)
    monkeypatch.setattr(tfv, "prescp", 1)
    monkeypatch.setattr(tfv, "i_tf_sup", cntrpst_asset.i_tf_sup)
    monkeypatch.setattr(tfv, "tcoolin", cntrpst_asset.tcoolin)
    monkeypatch.setattr(fwbsv, "pnuc_cp_tf", 1)
    monkeypatch.setattr(bv, "hmax", 1)
    monkeypatch.setattr(bv, "tfthko", 0.5)

    tfcoil.cntrpst()

    # appears to be the same for all cases?
    assert pytest.approx(tfv.ncool) == 203718.3271576

    assert pytest.approx(tfv.dtiocool, abs=1e-8) == cntrpst_asset.expected_dtiocool
    assert pytest.approx(tfv.tcpav2) == cntrpst_asset.expected_tcpav2
    assert pytest.approx(tfv.tcpmax) == cntrpst_asset.expected_tcpmax
    assert pytest.approx(tfv.ppump) == cntrpst_asset.expected_ppump
