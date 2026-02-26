import logging
from enum import IntEnum

import numpy as np

from process import constants
from process import process_output as po
from process.data_structure import (
    current_drive_variables,
    divertor_variables,
    physics_variables,
)
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)


class DensityLimitModel(IntEnum):
    """Electron density model types"""

    ASDEX = 1
    BORRASS_ITER_I = 2
    BORRASS_ITER_II = 3
    JET_EDGE_RADIATION = 4
    JET_SIMPLE = 5
    HUGILL_MURAKAMI = 6
    GREENWALD = 7
    ASDEX_NEW = 8


class PlasmaDensityLimit:
    """Class to hold plasma density limit calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self):
        nd_plasma_electron_max_array = np.empty((8,))

        p_perp = (
            physics_variables.p_plasma_separatrix_mw / physics_variables.a_plasma_surface
        )

        # Old ASDEX density limit formula
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        nd_plasma_electron_max_array[0] = self.calculate_asdex_density_limit(
            p_perp=p_perp,
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            q95=physics_variables.q95,
            rmajor=physics_variables.rmajor,
            prn1=divertor_variables.prn1,
        )

        # Borrass density limit model for ITER (I)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # Borrass et al, ITER-TN-PH-9-6 (1989)

        nd_plasma_electron_max_array[1] = self.calculate_borrass_iter_i_density_limit(
            p_perp=p_perp,
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            q95=physics_variables.q95,
            rmajor=physics_variables.rmajor,
            prn1=divertor_variables.prn1,
        )

        # Borrass density limit model for ITER (II)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # This formula is (almost) identical to that in the original routine
        # denlim (now deleted).

        nd_plasma_electron_max_array[2] = self.calculate_borrass_iter_ii_density_limit(
            p_perp=p_perp,
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            q95=physics_variables.q95,
            rmajor=physics_variables.rmajor,
            prn1=divertor_variables.prn1,
        )

        # JET edge radiation density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # qcyl=qstar here, but literature is not clear.

        nd_plasma_electron_max_array[3] = (
            self.calculate_jet_edge_radiation_density_limit(
                zeff=physics_variables.n_charge_plasma_effective_vol_avg,
                p_hcd_injected_total_mw=current_drive_variables.p_hcd_injected_total_mw,
                prn1=divertor_variables.prn1,
                qcyl=physics_variables.qstar,
            )
        )

        # JET simplified density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.

        nd_plasma_electron_max_array[4] = self.calculate_jet_simple_density_limit(
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            p_plasma_separatrix_mw=physics_variables.p_plasma_separatrix_mw,
            rmajor=physics_variables.rmajor,
            prn1=divertor_variables.prn1,
        )

        # Hugill-Murakami M.q limit
        # qcyl=qstar here, which is okay according to the literature

        nd_plasma_electron_max_array[5] = self.calculate_hugill_murakami_density_limit(
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            rmajor=physics_variables.rmajor,
            qcyl=physics_variables.qstar,
        )

        # Greenwald limit

        nd_plasma_electron_max_array[6] = self.calculate_greenwald_density_limit(
            c_plasma=physics_variables.plasma_current, rminor=physics_variables.rminor
        )

        nd_plasma_electron_max_array[7] = self.calculate_asdex_new_density_limit(
            p_hcd_injected_total_mw=current_drive_variables.p_hcd_injected_total_mw,
            c_plasma=physics_variables.plasma_current,
            q95=physics_variables.q95,
            prn1=divertor_variables.prn1,
        )

        physics_variables.nd_plasma_electron_max_array = nd_plasma_electron_max_array

        # Calculate beta_norm_max based on i_beta_norm_max
        try:
            model = DensityLimitModel(int(physics_variables.i_density_limit))
            physics_variables.nd_plasma_electrons_max = self.get_density_limit_value(
                model
            )
        except ValueError:
            raise ProcessValueError(
                "Illegal value of i_density_limit",
                i_density_limit=physics_variables.i_density_limit,
            ) from None

    def get_density_limit_value(self, model: DensityLimitModel) -> float:
        """Get the density limit value (n_e_max) for the specified model."""
        model_map = {
            DensityLimitModel.ASDEX: physics_variables.nd_plasma_electron_max_array[0],
            DensityLimitModel.BORRASS_ITER_I: physics_variables.nd_plasma_electron_max_array[
                1
            ],
            DensityLimitModel.BORRASS_ITER_II: physics_variables.nd_plasma_electron_max_array[
                2
            ],
            DensityLimitModel.JET_EDGE_RADIATION: physics_variables.nd_plasma_electron_max_array[
                3
            ],
            DensityLimitModel.JET_SIMPLE: physics_variables.nd_plasma_electron_max_array[
                4
            ],
            DensityLimitModel.HUGILL_MURAKAMI: physics_variables.nd_plasma_electron_max_array[
                5
            ],
            DensityLimitModel.GREENWALD: physics_variables.nd_plasma_electron_max_array[
                6
            ],
            DensityLimitModel.ASDEX_NEW: physics_variables.nd_plasma_electron_max_array[
                7
            ],
        }
        return model_map[model]

    @staticmethod
    def calculate_asdex_density_limit(
        p_perp: float,
        b_plasma_toroidal_on_axis: float,
        q95: float,
        rmajor: float,
        prn1: float,
    ) -> float:
        """
        Calculate the ASDEX density limit.

        :param p_perp: Perpendicular power density (MW/m²).
        :type p_perp: float
        :param b_plasma_toroidal_on_axis: Toroidal field on axis (T).
        :type b_plasma_toroidal_on_axis: float
        :param q95: Safety factor at 95% of the plasma poloidal flux.
        :type q95: float
        :param rmajor: Plasma major radius (m).
        :type rmajor: float
        :param prn1: Edge density / average plasma density.
        :type prn1: float
        :return: The ASDEX density limit (m⁻³).
        :rtype: float

        :references:
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

        """
        return (
            1.54e20
            * p_perp**0.43
            * b_plasma_toroidal_on_axis**0.31
            / (q95 * rmajor) ** 0.45
        ) / prn1

    @staticmethod
    def calculate_borrass_iter_i_density_limit(
        p_perp: float,
        b_plasma_toroidal_on_axis: float,
        q95: float,
        rmajor: float,
        prn1: float,
    ) -> float:
        """
        Calculate the Borrass ITER I density limit.

        :param p_perp: Perpendicular power density (MW/m²).
        :type p_perp: float
        :param b_plasma_toroidal_on_axis: Toroidal field on axis (T).
        :type b_plasma_toroidal_on_axis: float
        :param q95: Safety factor at 95% of the plasma poloidal flux.
        :type q95: float
        :param rmajor: Plasma major radius (m).
        :type rmajor: float
        :param prn1: Edge density / average plasma density.
        :type prn1: float
        :return: The Borrass ITER I density limit (m⁻³).
        :rtype: float


        :references:
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

        """
        return (
            1.8e20
            * p_perp**0.53
            * b_plasma_toroidal_on_axis**0.31
            / (q95 * rmajor) ** 0.22
        ) / prn1

    @staticmethod
    def calculate_borrass_iter_ii_density_limit(
        p_perp: float,
        b_plasma_toroidal_on_axis: float,
        q95: float,
        rmajor: float,
        prn1: float,
    ) -> float:
        """
        Calculate the Borrass ITER II density limit.

        :param p_perp: Perpendicular power density (MW/m²).
        :type p_perp: float
        :param b_plasma_toroidal_on_axis: Toroidal field on axis (T).
        :type b_plasma_toroidal_on_axis: float
        :param q95: Safety factor at 95% of the plasma poloidal flux.
        :type q95: float
        :param rmajor: Plasma major radius (m).
        :type rmajor: float
        :param prn1: Edge density / average plasma density.
        :type prn1: float
        :return: The Borrass ITER II density limit (m⁻³).
        :rtype: float


        :references:
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

        """
        return (
            0.5e20
            * p_perp**0.57
            * b_plasma_toroidal_on_axis**0.31
            / (q95 * rmajor) ** 0.09
        ) / prn1

    @staticmethod
    def calculate_jet_edge_radiation_density_limit(
        zeff: float, p_hcd_injected_total_mw: float, prn1: float, qcyl: float
    ) -> float:
        """
        Calculate the JET edge radiation density limit.

        :param zeff: Effective charge (Z_eff).
        :type zeff: float
        :param p_hcd_injected_total_mw: Power injected into the plasma (MW).
        :type p_hcd_injected_total_mw: float
        :param prn1: Edge density / average plasma density.
        :type prn1: float
        :param qcyl: Equivalent cylindrical safety factor (qstar).
        :type qcyl: float
        :return: The JET edge radiation density limit (m⁻³).
        :rtype: float

        :references:
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992
        """

        denom = (zeff - 1.0) * (1.0 - 4.0 / (3.0 * qcyl))
        if denom <= 0.0:
            return 0.0
        return (1.0e20 * np.sqrt(p_hcd_injected_total_mw / denom)) / prn1

    @staticmethod
    def calculate_jet_simple_density_limit(
        b_plasma_toroidal_on_axis: float,
        p_plasma_separatrix_mw: float,
        rmajor: float,
        prn1: float,
    ) -> float:
        """
        Calculate the JET simple density limit.

        :param b_plasma_toroidal_on_axis: Toroidal field on axis (T).
        :type b_plasma_toroidal_on_axis: float
        :param p_plasma_separatrix_mw: Power crossing the separatrix (MW).
        :type p_plasma_separatrix_mw: float
        :param rmajor: Plasma major radius (m).
        :type rmajor: float
        :param prn1: Edge density / average plasma density.
        :type prn1: float
        :return: The JET simple density limit (m⁻³).
        :rtype: float

        :references:
            - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study', AEA FUS 172, 1992

        """
        return (
            0.237e20
            * b_plasma_toroidal_on_axis
            * np.sqrt(p_plasma_separatrix_mw)
            / rmajor
        ) / prn1

    @staticmethod
    def calculate_hugill_murakami_density_limit(
        b_plasma_toroidal_on_axis: float, rmajor: float, qcyl: float
    ) -> float:
        """
        Calculate the Hugill-Murakami density limit.

        :param b_plasma_toroidal_on_axis: Toroidal field on axis (T).
        :type b_plasma_toroidal_on_axis: float
        :param rmajor: Plasma major radius (m).
        :type rmajor: float
        :param qcyl: Equivalent cylindrical safety factor (qstar).
        :type qcyl: float
        :return: The Hugill-Murakami density limit (m⁻³).
        :rtype: float

        :references:
            - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',
        """

        return 3.0e20 * b_plasma_toroidal_on_axis / (rmajor * qcyl)

    @staticmethod
    def calculate_greenwald_density_limit(c_plasma: float, rminor: float) -> float:
        """
        Calculate the Greenwald density limit (n_GW).

        :param c_plasma: Plasma current (A).
        :type c_plasma: float
        :param rminor: Plasma minor radius (m).
        :type rminor: float
        :return: The Greenwald density limit (m⁻³).
        :rtype: float

        :notes: The Greenwald limit is typically applied to the line averaged electron density

        :references:
            - M. Greenwald et al., “A new look at density limits in tokamaks,”
            Nuclear Fusion, vol. 28, no. 12, pp. 2199-2207, Dec. 1988,
            doi: https://doi.org/10.1088/0029-5515/28/12/009.

            - M. Greenwald, “Density limits in toroidal plasmas,”
            Plasma Physics and Controlled Fusion, vol. 44, no. 8, pp. R27-R53, Jul. 2002,
            doi: https://doi.org/10.1088/0741-3335/44/8/201.
        """

        return 1.0e14 * c_plasma / (np.pi * rminor**2)

    @staticmethod
    def calculate_asdex_new_density_limit(
        p_hcd_injected_total_mw: float, c_plasma: float, q95: float, prn1: float
    ) -> float:
        """
        Calculate the ASDEX Upgrade new density limit.

        :param p_hcd_injected_total_mw: Power injected into the plasma (MW).
        :type p_hcd_injected_total_mw: float
        :param plasma_current: Plasma current (A).
        :type plasma_current: float
        :param q95: Safety factor at 95% surface.
        :type q95: float
        :param prn1: Edge density / average plasma density.
        :type prn1: float
        :return: The ASDEX Upgrade new density limit (m⁻³).
        :rtype: float

        :notes: This limit is for the separatrix density so wee scale by `prn1` to get it as a volume average

        :references:

            - J. W. Berkery et al., “Density limits as disruption forecasters for spherical tokamaks,”
            Plasma Physics and Controlled Fusion, vol. 65, no. 9, pp. 095003-095003, Jul. 2023,
            doi: https://doi.org/10.1088/1361-6587/ace476.

            - M. Bernert et al., “The H-mode density limit in the full tungsten ASDEX Upgrade tokamak,” vol. 57, no. 1, pp. 014038-014038, Nov. 2014,
            doi: https://doi.org/10.1088/0741-3335/57/1/014038.
        """
        return (
            1.0e20
            * 0.506
            * (p_hcd_injected_total_mw**0.396 * (c_plasma / 1.0e6) ** 0.265)
            / (q95**0.323)
        ) / prn1

    def calculate_density_limit(
        self,
        b_plasma_toroidal_on_axis: float,
        i_density_limit: int,
        p_plasma_separatrix_mw: float,
        p_hcd_injected_total_mw: float,
        plasma_current: float,
        prn1: float,
        qcyl: float,
        q95: float,
        rmajor: float,
        rminor: float,
        a_plasma_surface: float,
        zeff: float,
    ) -> tuple[np.ndarray, float]:
        """
        Calculate the density limit using various models.

        Args:
            b_plasma_toroidal_on_axis (float): Toroidal field on axis (T).
            i_density_limit (int): Switch denoting which formula to enforce.
            p_plasma_separatrix_mw (float): Power flowing to the edge plasma via charged particles (MW).
            p_hcd_injected_total_mw (float): Power injected into the plasma (MW).
            plasma_current (float): Plasma current (A).
            prn1 (float): Edge density / average plasma density.
            qcyl (float): Equivalent cylindrical safety factor (qstar).
            q95 (float): Safety factor at 95% surface.
            rmajor (float): Plasma major radius (m).
            rminor (float): Plasma minor radius (m).
            a_plasma_surface (float): Plasma surface area (m^2).
            zeff (float): Plasma effective charge.

        Returns:
            Tuple[np.ndarray, float]: A tuple containing:
                - nd_plasma_electron_max_array (np.ndarray): Average plasma density limit using seven different models (m^-3).
                - nd_plasma_electrons_max (float): Enforced average plasma density limit (m^-3).

        Raises:
            ValueError: If i_density_limit is not between 1 and 7.

        Notes:
            This routine calculates several different formulae for the density limit and enforces the one chosen by the user.
            For i_density_limit = 1-5, 8, we scale the sepatrix density limit output by the ratio of the separatrix to volume averaged density

        References:
            - AEA FUS 172: Physics Assessment for the European Reactor Study

            - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989

            - M. Bernert et al., “The H-mode density limit in the full tungsten ASDEX Upgrade tokamak,”
              vol. 57, no. 1, pp. 014038-014038, Nov. 2014, doi: https://doi.org/10.1088/0741-3335/57/1/014038. ‌
        """

        if i_density_limit < 1 or i_density_limit > 7:
            raise ProcessValueError(
                "Illegal value for i_density_limit", i_density_limit=i_density_limit
            )

        nd_plasma_electron_max_array = np.empty((8,))

        # Power per unit area crossing the plasma edge
        # (excludes radiation and neutrons)

        p_perp = p_plasma_separatrix_mw / a_plasma_surface

        # Old ASDEX density limit formula
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.

        nd_plasma_electron_max_array[0] = self.calculate_asdex_density_limit(
            p_perp=p_perp,
            b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
            q95=q95,
            rmajor=rmajor,
            prn1=prn1,
        )

        # Borrass density limit model for ITER (I)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # Borrass et al, ITER-TN-PH-9-6 (1989)

        nd_plasma_electron_max_array[1] = self.calculate_borrass_iter_i_density_limit(
            p_perp=p_perp,
            b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
            q95=q95,
            rmajor=rmajor,
            prn1=prn1,
        )

        # Borrass density limit model for ITER (II)
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # This formula is (almost) identical to that in the original routine
        # denlim (now deleted).

        nd_plasma_electron_max_array[2] = self.calculate_borrass_iter_ii_density_limit(
            p_perp=p_perp,
            b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
            q95=q95,
            rmajor=rmajor,
            prn1=prn1,
        )

        # JET edge radiation density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.
        # qcyl=qstar here, but literature is not clear.

        nd_plasma_electron_max_array[3] = (
            self.calculate_jet_edge_radiation_density_limit(
                zeff=zeff,
                p_hcd_injected_total_mw=p_hcd_injected_total_mw,
                prn1=prn1,
                qcyl=qcyl,
            )
        )

        # JET simplified density limit model
        # This applies to the density at the plasma edge, so must be scaled
        # to give the density limit applying to the average plasma density.

        nd_plasma_electron_max_array[4] = self.calculate_jet_simple_density_limit(
            b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
            p_plasma_separatrix_mw=p_plasma_separatrix_mw,
            rmajor=rmajor,
            prn1=prn1,
        )

        # Hugill-Murakami M.q limit
        # qcyl=qstar here, which is okay according to the literature

        nd_plasma_electron_max_array[5] = self.calculate_hugill_murakami_density_limit(
            b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis, rmajor=rmajor, qcyl=qcyl
        )

        # Greenwald limit

        nd_plasma_electron_max_array[6] = self.calculate_greenwald_density_limit(
            c_plasma=plasma_current, rminor=rminor
        )

        nd_plasma_electron_max_array[7] = self.calculate_asdex_new_density_limit(
            p_hcd_injected_total_mw=p_hcd_injected_total_mw,
            c_plasma=plasma_current,
            q95=q95,
            prn1=prn1,
        )

        # Enforce the chosen density limit

        return nd_plasma_electron_max_array, nd_plasma_electron_max_array[
            i_density_limit - 1
        ]

    def output_density_limit_information(self):
        """Output density limit information to file."""

        po.osubhd(self.outfile, "Density Limit using different models :")
        po.ovarre(
            self.outfile,
            "Old ASDEX model",
            "(nd_plasma_electron_max_array(1))",
            physics_variables.nd_plasma_electron_max_array[0],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Borrass ITER model I",
            "(nd_plasma_electron_max_array(2))",
            physics_variables.nd_plasma_electron_max_array[1],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Borrass ITER model II",
            "(nd_plasma_electron_max_array(3))",
            physics_variables.nd_plasma_electron_max_array[2],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "JET edge radiation model",
            "(nd_plasma_electron_max_array(4))",
            physics_variables.nd_plasma_electron_max_array[3],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "JET simplified model",
            "(nd_plasma_electron_max_array(5))",
            physics_variables.nd_plasma_electron_max_array[4],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Hugill-Murakami Mq model",
            "(nd_plasma_electron_max_array(6))",
            physics_variables.nd_plasma_electron_max_array[5],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Greenwald model",
            "(nd_plasma_electron_max_array(7))",
            physics_variables.nd_plasma_electron_max_array[6],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "ASDEX New",
            "(nd_plasma_electron_max_array(8))",
            physics_variables.nd_plasma_electron_max_array[7],
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Density limit from scaling (/m3)",
            "(nd_plasma_electrons_max)",
            physics_variables.nd_plasma_electrons_max,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ostars(self.outfile, 110)
        po.oblnkl(self.outfile)
