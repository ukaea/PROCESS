"""Plasma current scaling models and calculations.

This module provides various plasma current scaling models used in fusion
physics calculations, including empirical and analytical fits based on
different tokamak designs and experimental data.
"""

import logging
from enum import IntEnum
from types import DynamicClassAttribute

import numba as nb
import numpy as np

from process.core import constants
from process.core import process_output as po
from process.core.exceptions import ProcessValueError
from process.core.model import Model
from process.data_structure import (
    current_drive_variables,
    physics_variables,
    stellarator_variables,
)
from process.models.physics.plasma_geometry import PlasmaGeometryModelType

logger = logging.getLogger(__name__)


class PlasmaCurrentModel(IntEnum):
    """Enumeration of plasma current scaling models available for calculations.

    Each model represents a different scaling law used to calculate plasma
    current based on various plasma and machine parameters.
    """

    PENG_ANALYTIC_FIT = (1, "Peng analytic fit")
    PENG_DIVERTOR_SCALING = (2, "Peng divertor scaling")
    ITER_SCALING = (3, "Simple ITER scaling (cylindrical case)")
    IPDG89_SCALING = (4, "IPDG89 scaling")
    TODD_EMPIRICAL_SCALING_I = (5, "Todd empirical scaling I")
    TODD_EMPIRICAL_SCALING_II = (6, "Todd empirical scaling II")
    CONNOR_HASTIE_MODEL = (7, "Connor-Hastie model")
    SAUTER_SCALING = (8, "Sauter scaling")
    FIESTA_ST_SCALING = (9, "FIESTA ST scaling")

    def __new__(cls, value: int, full_name: str):
        """Create a new PlasmaCurrentModel enum member with value and full_name.

        Parameters
        ----------
        value : int
            The numeric value of the enum member.
        full_name : str
            The full name description of the plasma current model.

        Returns
        -------
        PlasmaCurrentModel
            A new enum member with the specified value and full_name.
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._full_name_ = full_name
        return obj

    @DynamicClassAttribute
    def full_name(self):
        """Return the full name of the plasma current model."""
        return self._full_name_


class PlasmaCurrent:
    """Class to hold plasma current calculations for plasma processing."""

    def __init__(self):
        """Initialize PlasmaCurrent with output file constants."""
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def output(self) -> None:
        """Output plasma current and safety factor information."""
        po.oheadr(self.outfile, "Plasma Current and Safety Factor")

        if stellarator_variables.istell == 0:
            po.oblnkl(self.outfile)
            po.ovarin(
                self.outfile,
                "Plasma current scaling law used",
                "(i_plasma_current)",
                physics_variables.i_plasma_current,
            )
            full_model_name = PlasmaCurrentModel(
                physics_variables.i_plasma_current
            ).full_name
            po.ocmmnt(
                self.outfile,
                f"Plasma current model selected: {full_model_name} ",
            )

            po.ovarrf(
                self.outfile,
                "Plasma current (Iₚ) (MA)",
                "(plasma_current_MA)",
                physics_variables.plasma_current / 1.0e6,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Plasma current (Iₚ) (A)",
                "(plasma_current)",
                physics_variables.plasma_current,
                "OP ",
            )

            self.output_plasma_current_models()
            po.oblnkl(self.outfile)

            if physics_variables.i_alphaj == 1:
                po.ovarrf(
                    self.outfile,
                    "Current density profile factor (αⱼ)",
                    "(alphaj)",
                    physics_variables.alphaj,
                    "OP ",
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "Current density profile factor (αⱼ)",
                    "(alphaj)",
                    physics_variables.alphaj,
                )
            po.ocmmnt(self.outfile, "Current profile index scalings:")
            po.oblnkl(self.outfile)

            po.ovarrf(
                self.outfile,
                "J. Wesson plasma current profile index",
                "(alphaj_wesson)",
                physics_variables.alphaj_wesson,
                "OP ",
            )
            po.oblnkl(self.outfile)
            po.ovarrf(
                self.outfile,
                "On-axis plasma current density (j₀) (A/m²)",
                "(j_plasma_on_axis)",
                physics_variables.j_plasma_on_axis,
                "OP ",
            )

        if stellarator_variables.istell == 0:
            po.ovarrf(
                self.outfile, "Safety factor on axis (q₀)", "(q0)", physics_variables.q0
            )

            if physics_variables.i_plasma_current == 2:
                po.ovarrf(
                    self.outfile,
                    "Mean edge safety factor (q₉₅)",
                    "(q95)",
                    physics_variables.q95,
                )

            po.ovarrf(
                self.outfile,
                "Safety factor at 95% flux surface (q₉₅)",
                "(q95)",
                physics_variables.q95,
            )

            po.ovarrf(
                self.outfile,
                "Cylindrical safety factor (qcyl)",
                "(qstar)",
                physics_variables.qstar,
                "OP ",
            )

            if (
                physics_variables.i_plasma_geometry
                == PlasmaGeometryModelType.STAR_FIESTA
            ):
                po.ovarrf(
                    self.outfile,
                    "Lower limit for edge safety factor q95",
                    "(q95_min)",
                    physics_variables.q95_min,
                    "OP ",
                )
            po.oblnkl(self.outfile)

        else:
            po.ovarrf(
                self.outfile,
                "Rotational transform",
                "(iotabar)",
                stellarator_variables.iotabar,
            )

    def calculate_plasma_current(
        self,
        alphaj: float,
        alphap: float,
        b_plasma_toroidal_on_axis: float,
        eps: float,
        i_plasma_current: int,
        kappa: float,
        kappa95: float,
        pres_plasma_on_axis: float,
        len_plasma_poloidal: float,
        q95: float,
        rmajor: float,
        rminor: float,
        triang: float,
        triang95: float,
    ) -> tuple[float, float, float, float, float]:
        """Calculate the plasma current.

        Parameters
        ----------
        alphaj : float
            Current profile index.
        alphap : float
            Pressure profile index.
        b_plasma_toroidal_on_axis : float
            Toroidal field on axis (T).
        eps : float
            Inverse aspect ratio.
        i_plasma_current : int
            Current scaling model to use.
            1 = Peng analytic fit
            2 = Peng divertor scaling (TART,STAR)
            3 = Simple ITER scaling
            4 = IPDG89 scaling
            5 = Todd empirical scaling I
            6 = Todd empirical scaling II
            7 = Connor-Hastie model
            8 = Sauter scaling (allowing negative triangularity)
            9 = FIESTA ST scaling
        kappa : float
            Plasma elongation.
        kappa95 : float
            Plasma elongation at 95% surface.
        pres_plasma_on_axis : float
            Central plasma pressure (Pa).
        len_plasma_poloidal : float
            Plasma perimeter length (m).
        q95 : float
            Plasma safety factor at 95% flux (= q-bar for i_plasma_current=2).
        rmajor : float
            Major radius (m).
        rminor : float
            Minor radius (m).
        triang : float
            Plasma triangularity.
        triang95 : float
            Plasma triangularity at 95% surface.

        Returns
        -------
        tuple[float, float, float, float, float]
            Tuple containing (b_plasma_poloidal_average, qstar, plasma_current,
            betap, li).

        Raises
        ------
        ValueError
            If invalid value for `i_plasma_current` is provided.
        ProcessValueError
            If triangularity is negative without `i_plasma_current = 8` selected.

        Notes
        -----
        This routine calculates the plasma current based on the edge safety factor
        q95. It will also make the current profile parameters consistent with the
        q-profile if required.

        References
        ----------
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
              unpublished internal Oak Ridge document

            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
              'Small Tokamaks for Fusion Technology Testing'. Fusion Technology,
              21(3P2A), 1729-1738. https://doi.org/10.13182/FST92-A29971

            - ITER Physics Design Guidelines: 1989 [IPDG89], N. A. Uckan et al,
              ITER Documentation Series No.10, IAEA/ITER/DS/10, IAEA, Vienna, 1990

            - M. Kovari et al, 2014, "PROCESS": A systems code for fusion power plants
              - Part 1: Physics

            - H. Zohm et al, 2013, On the Physics Guidelines for a Tokamak DEMO

            - T. Hartmann, 2013, Development of a modular systems code to analyse the
              implications of physics assumptions on the design of a demonstration
              fusion power plant

            - Sauter, Geometric formulas for systems codes..., FED 2016
        """
        # Aspect ratio
        aspect_ratio = 1.0 / eps

        # Only the Sauter scaling (i_plasma_current=8) is suitable for negative
        # triangularity:
        if i_plasma_current != 8 and triang < 0.0:
            raise ProcessValueError(
                f"Triangularity is negative without i_plasma_current = 8 selected:"
                f" {triang=}, {i_plasma_current=}"
            )
        try:
            model = PlasmaCurrentModel(int(i_plasma_current))
            # Calculate the function Fq that scales the edge q from the
            # circular cross-section cylindrical case

            # Peng analytical fit
            if model == PlasmaCurrentModel.PENG_ANALYTIC_FIT:
                fq = self.calculate_current_coefficient_peng(
                    eps, len_plasma_poloidal, rminor
                )

            # Peng scaling for double null divertor; TARTs [STAR Code]
            elif model == PlasmaCurrentModel.PENG_DIVERTOR_SCALING:
                plasma_current = 1.0e6 * self.calculate_plasma_current_peng(
                    q95,
                    aspect_ratio,
                    eps,
                    rminor,
                    b_plasma_toroidal_on_axis,
                    kappa,
                    triang,
                )

            # Simple ITER scaling (simply the cylindrical case)
            elif model == PlasmaCurrentModel.ITER_SCALING:
                fq = 1.0

            # ITER formula (IPDG89)
            elif model == PlasmaCurrentModel.IPDG89_SCALING:
                fq = self.calculate_current_coefficient_ipdg89(eps, kappa95, triang95)

            # Todd empirical scalings
            # D.C.Robinson and T.N.Todd, Plasma and Contr Fusion 28 (1986) 1181
            elif model in {
                PlasmaCurrentModel.TODD_EMPIRICAL_SCALING_I,
                PlasmaCurrentModel.TODD_EMPIRICAL_SCALING_II,
            }:
                fq = self.calculate_current_coefficient_todd(
                    eps, kappa95, triang95, model=1
                )

                if model == PlasmaCurrentModel.TODD_EMPIRICAL_SCALING_II:
                    fq = self.calculate_current_coefficient_todd(
                        eps, kappa95, triang95, model=2
                    )

            # Connor-Hastie asymptotically-correct expression
            elif model == PlasmaCurrentModel.CONNOR_HASTIE_MODEL:
                fq = self.calculate_current_coefficient_hastie(
                    alphaj,
                    alphap,
                    b_plasma_toroidal_on_axis,
                    triang95,
                    eps,
                    kappa95,
                    pres_plasma_on_axis,
                    constants.RMU0,
                )

            # Sauter scaling allowing negative triangularity [FED May 2016]
            # https://doi.org/10.1016/j.fusengdes.2016.04.033.
            elif model == PlasmaCurrentModel.SAUTER_SCALING:
                # Assumes zero squareness, note takes kappa, delta at separatrix not _95
                fq = self.calculate_current_coefficient_sauter(eps, kappa, triang)

            # FIESTA ST scaling
            # https://doi.org/10.1016/j.fusengdes.2020.111530.
            elif model == PlasmaCurrentModel.FIESTA_ST_SCALING:
                fq = self.calculate_current_coefficient_fiesta(eps, kappa, triang)

        except ValueError as e:
            raise ProcessValueError(
                "Illegal value of i_plasma_current",
                i_plasma_current=physics_variables.i_plasma_current,
            ) from e

        # Main plasma current calculation using the fq value from the different settings
        if model != PlasmaCurrentModel.PENG_DIVERTOR_SCALING:
            plasma_current = (
                self.calculate_cyclindrical_plasma_current(
                    rminor=rminor,
                    rmajor=rmajor,
                    q95=q95,
                    b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
                )
                * fq
            )

        return plasma_current

    def calculate_all_plasma_current_models(
        self,
        alphaj: float,
        alphap: float,
        b_plasma_toroidal_on_axis: float,
        eps: float,
        kappa: float,
        kappa95: float,
        pres_plasma_on_axis: float,
        len_plasma_poloidal: float,
        q95: float,
        rmajor: float,
        rminor: float,
        triang: float,
        triang95: float,
    ) -> dict[PlasmaCurrentModel, float]:
        """Calculate the plasma current for all models.

        This function calculates the plasma current for all models and returns a
        dictionary of the results.

        Parameters
        ----------
        alphaj :
            current profile index
        alphap :
            pressure profile index
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        eps :
            inverse aspect ratio
        kappa :
            plasma elongation
        kappa95 :
            plasma elongation at 95% surface
        pres_plasma_on_axis :
            central plasma pressure (Pa)
        len_plasma_poloidal :
            plasma perimeter length (m)
        q95 :
            plasma safety factor at 95% flux
        rmajor :
            major radius (m)
        rminor :
            minor radius (m)
        triang :
            plasma triangularity
        triang95 :
            plasma triangularity at 95% surface

        Returns
        -------
        dict[PlasmaCurrentModel, float]
            Dictionary containing the plasma current for each model.
        """
        results = {}
        for model in PlasmaCurrentModel:
            results[model] = self.calculate_plasma_current(
                alphaj=alphaj,
                alphap=alphap,
                b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
                eps=eps,
                i_plasma_current=model.value,
                kappa=kappa,
                kappa95=kappa95,
                pres_plasma_on_axis=pres_plasma_on_axis,
                len_plasma_poloidal=len_plasma_poloidal,
                q95=q95,
                rmajor=rmajor,
                rminor=rminor,
                triang=triang,
                triang95=triang95,
            )
        return results

    def output_plasma_current_models(self) -> None:
        """Output the plasma current for all models.

        This function outputs the plasma current for all models to the output file.
        """
        plasma_currents = self.calculate_all_plasma_current_models(
            alphaj=physics_variables.alphaj,
            alphap=physics_variables.alphap,
            b_plasma_toroidal_on_axis=physics_variables.b_plasma_toroidal_on_axis,
            eps=physics_variables.eps,
            kappa=physics_variables.kappa,
            kappa95=physics_variables.kappa95,
            pres_plasma_on_axis=physics_variables.pres_plasma_thermal_on_axis,
            len_plasma_poloidal=physics_variables.len_plasma_poloidal,
            q95=physics_variables.q95,
            rmajor=physics_variables.rmajor,
            rminor=physics_variables.rminor,
            triang=physics_variables.triang,
            triang95=physics_variables.triang95,
        )

        physics_variables.c_plasma_peng_analytic = plasma_currents.get(
            PlasmaCurrentModel.PENG_ANALYTIC_FIT
        )
        physics_variables.c_plasma_peng_double_null = plasma_currents.get(
            PlasmaCurrentModel.PENG_DIVERTOR_SCALING
        )
        physics_variables.c_plasma_cyclindrical = plasma_currents.get(
            PlasmaCurrentModel.ITER_SCALING
        )
        physics_variables.c_plasma_ipdg89 = plasma_currents.get(
            PlasmaCurrentModel.IPDG89_SCALING
        )
        physics_variables.c_plasma_todd_empirical_i = plasma_currents.get(
            PlasmaCurrentModel.TODD_EMPIRICAL_SCALING_I
        )
        physics_variables.c_plasma_todd_empirical_ii = plasma_currents.get(
            PlasmaCurrentModel.TODD_EMPIRICAL_SCALING_II
        )
        physics_variables.c_plasma_connor_hastie = plasma_currents.get(
            PlasmaCurrentModel.CONNOR_HASTIE_MODEL
        )
        physics_variables.c_plasma_sauter = plasma_currents.get(
            PlasmaCurrentModel.SAUTER_SCALING
        )
        physics_variables.c_plasma_fiesta_st = plasma_currents.get(
            PlasmaCurrentModel.FIESTA_ST_SCALING
        )

        po.osubhd(self.outfile, "Plasma Currents using different models :")

        po.ovarre(
            self.outfile,
            f"{PlasmaCurrentModel.PENG_ANALYTIC_FIT.full_name}",
            "(c_plasma_peng_analytic)",
            physics_variables.c_plasma_peng_analytic,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"{PlasmaCurrentModel.PENG_DIVERTOR_SCALING.full_name}",
            "(c_plasma_peng_double_null)",
            physics_variables.c_plasma_peng_double_null,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"{PlasmaCurrentModel.ITER_SCALING.full_name}",
            "(c_plasma_cyclindrical)",
            physics_variables.c_plasma_cyclindrical,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"{PlasmaCurrentModel.IPDG89_SCALING.full_name}",
            "(c_plasma_ipdg89)",
            physics_variables.c_plasma_ipdg89,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"{PlasmaCurrentModel.TODD_EMPIRICAL_SCALING_I.full_name}",
            "(c_plasma_todd_empirical_i)",
            physics_variables.c_plasma_todd_empirical_i,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"{PlasmaCurrentModel.TODD_EMPIRICAL_SCALING_II.full_name}",
            "(c_plasma_todd_empirical_ii)",
            physics_variables.c_plasma_todd_empirical_ii,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"{PlasmaCurrentModel.CONNOR_HASTIE_MODEL.full_name}",
            "(c_plasma_connor_hastie)",
            physics_variables.c_plasma_connor_hastie,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"{PlasmaCurrentModel.SAUTER_SCALING.full_name}",
            "(c_plasma_sauter)",
            physics_variables.c_plasma_sauter,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            f"{PlasmaCurrentModel.FIESTA_ST_SCALING.full_name}",
            "(c_plasma_fiesta_st)",
            physics_variables.c_plasma_fiesta_st,
            "OP ",
        )

    @staticmethod
    def calculate_cyclindrical_plasma_current(
        rminor: float, rmajor: float, q95: float, b_plasma_toroidal_on_axis: float
    ) -> float:
        """Calculate the plasma current for a cylindrical plasma.

        Parameters
        ----------
        rminor :
            plasma minor radius (m)
        rmajor :
            plasma major radius (m)
        q95 :
            plasma safety factor at 95% flux
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)

        Returns
        -------
        float
            plasma current (A)

        """
        return (
            (2.0 * np.pi / constants.RMU0)
            * rminor**2
            / (rmajor * q95)
            * b_plasma_toroidal_on_axis
        )

    @staticmethod
    def plascar_bpol(
        aspect: float, eps: float, kappa: float, delta: float
    ) -> tuple[float, float, float, float]:
        """Calculate the poloidal field coefficients for determining the plasma current
        and poloidal field.


        This internal function calculates the poloidal field coefficients,
        which is used to calculate the poloidal field and the plasma current.

        Parameters
        ----------
        aspect :
            plasma aspect ratio
        eps :
            inverse aspect ratio
        kappa :
            plasma elongation
        delta :
            plasma triangularity

        Returns
        -------
        :
            coefficients ff1, ff2, d1, d2

        References
        ----------
            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
            'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
            1729-1738. https://doi.org/10.13182/FST92-A29971
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
            unpublished internal Oak Ridge document

        """
        # Original coding, only suitable for TARTs [STAR Code]

        c1 = (kappa**2 / (1.0 + delta)) + delta
        c2 = (kappa**2 / (1.0 - delta)) - delta

        d1 = (kappa / (1.0 + delta)) ** 2 + 1.0
        d2 = (kappa / (1.0 - delta)) ** 2 + 1.0

        c1_aspect = ((c1 * eps) - 1.0) if aspect < c1 else (1.0 - (c1 * eps))

        y1 = np.sqrt(c1_aspect / (1.0 + eps)) * ((1.0 + delta) / kappa)
        y2 = np.sqrt((c2 * eps + 1.0) / (1.0 - eps)) * ((1.0 - delta) / kappa)

        h2 = (1.0 + (c2 - 1.0) * (eps / 2.0)) / np.sqrt((1.0 - eps) * (c2 * eps + 1.0))
        f2 = (d2 * (1.0 - delta) * eps) / ((1.0 - eps) * ((c2 * eps) + 1.0))
        g = (eps * kappa) / (1.0 - (eps * delta))
        ff2 = f2 * (g + 2.0 * h2 * np.arctan(y2))

        h1 = (1.0 + (1.0 - c1) * (eps / 2.0)) / np.sqrt((1.0 + eps) * c1_aspect)
        f1 = (d1 * (1.0 + delta) * eps) / ((1.0 + eps) * (c1 * eps - 1.0))

        if aspect < c1:
            ff1 = f1 * (g - h1 * np.log((1.0 + y1) / (1.0 - y1)))
        else:
            ff1 = -f1 * (-g + 2.0 * h1 * np.arctan(y1))

        return ff1, ff2, d1, d2

    def calculate_poloidal_field(
        self,
        i_plasma_current: int,
        ip: float,
        q95: float,
        aspect: float,
        eps: float,
        b_plasma_toroidal_on_axis: float,
        kappa: float,
        delta: float,
        perim: float,
    ) -> float:
        """Function to calculate poloidal field from the plasma current

        This function calculates the poloidal field from the plasma current in Tesla,
        using a simple calculation using Ampere's law for conventional
        tokamaks, or for TARTs, a scaling from Peng, Galambos and
        Shipe (1992).

        Parameters
        ----------
        i_plasma_current :
            current scaling model to use
        ip :
            plasma current (A)
        q95 :
            95% flux surface safety factor
        aspect :
            plasma aspect ratio
        eps :
            inverse aspect ratio
        b_plasma_toroidal_on_axis :
            toroidal field on axis (T)
        kappa :
            plasma elongation
        delta :
            plasma triangularity
        perim :
            plasma perimeter (m)

        Returns
        -------
        :
            poloidal field in Tesla


        References
        ----------
            - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
            unpublished internal Oak Ridge document
            - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
            'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
            1729-1738. https://doi.org/10.13182/FST92-A29971

        """
        # Use Ampere's law using the plasma poloidal cross-section
        if i_plasma_current != 2:
            return constants.RMU0 * ip / perim
        # Use the relation from Peng, Galambos and Shipe (1992) [STAR code] otherwise
        ff1, ff2, _, _ = self.plascar_bpol(aspect, eps, kappa, delta)

        # Transform q95 to qbar
        qbar = q95 * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0

        return b_plasma_toroidal_on_axis * (ff1 + ff2) / (2.0 * np.pi * qbar)

    @staticmethod
    def calculate_current_coefficient_peng(
        eps: float, len_plasma_poloidal: float, rminor: float
    ) -> float:
        """
        Calculate the plasma current scaling coefficient for the Peng scaling from
        the STAR code.

        Parameters
        ----------
        eps:
            Plasma inverse aspect ratio.
        len_plasma_poloidal:
            Plasma poloidal perimeter length (m).
        rminor:
            Plasma minor radius (m).

        Returns
        -------
        :
            The plasma current scaling coefficient.

        References
        ----------
            None
        """
        return (
            (1.22 - 0.68 * eps)
            / ((1.0 - eps * eps) ** 2)
            * (len_plasma_poloidal / (2.0 * np.pi * rminor)) ** 2
        )

    def calculate_plasma_current_peng(
        self,
        q95: float,
        aspect: float,
        eps: float,
        rminor: float,
        b_plasma_toroidal_on_axis: float,
        kappa: float,
        delta: float,
    ) -> float:
        """
        Function to calculate plasma current (Peng scaling from the STAR code)

        Parameters
        ----------
        - q95: float, 95% flux surface safety factor
        - aspect: float, plasma aspect ratio
        - eps: float, inverse aspect ratio
        - rminor: float, plasma minor radius (m)
        - b_plasma_toroidal_on_axis: float, toroidal field on axis (T)
        - kappa: float, plasma elongation
        - delta: float, plasma triangularity

        Returns
        -------
        - float, plasma current in MA

        This function calculates the plasma current in MA,
        using a scaling from Peng, Galambos and Shipe (1992).
        It is primarily used for Tight Aspect Ratio Tokamaks and is
        selected via i_plasma_current=2.

        References
        ----------
        - J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document
        - Peng, Y. K. M., Galambos, J. D., & Shipe, P. C. (1992).
        'Small Tokamaks for Fusion Technology Testing'. Fusion Technology, 21(3P2A),
        1729-1738. https://doi.org/10.13182/FST92-A29971
        """
        # Transform q95 to qbar
        qbar = q95 * 1.3e0 * (1.0e0 - physics_variables.eps) ** 0.6e0

        ff1, ff2, d1, d2 = self.plascar_bpol(aspect, eps, kappa, delta)

        e1 = (2.0 * kappa) / (d1 * (1.0 + delta))
        e2 = (2.0 * kappa) / (d2 * (1.0 - delta))

        return (
            rminor
            * b_plasma_toroidal_on_axis
            / qbar
            * 5.0
            * kappa
            / (2.0 * np.pi**2)
            * (np.arcsin(e1) / e1 + np.arcsin(e2) / e2)
            * (ff1 + ff2)
        )

    @staticmethod
    def calculate_current_coefficient_ipdg89(
        eps: float, kappa95: float, triang95: float
    ) -> float:
        """
        Calculate the fq coefficient from the IPDG89 guidlines used in the plasma
        current scaling.

        Parameters
        ----------
        - eps: float, plasma inverse aspect ratio
        - kappa95: float, plasma elongation 95%
        - triang95: float, plasma triangularity 95%

        Returns
        -------
        - float, the fq plasma current coefficient

        This function calculates the fq coefficient used in the IPDG89 plasma current
        scaling, based on the given plasma parameters.

        References
        ----------
        - N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989'
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study',
          AEA FUS 172, 1992
        """
        return (
            0.5
            * (1.17 - 0.65 * eps)
            / ((1.0 - eps * eps) ** 2)
            * (1.0 + kappa95**2 * (1.0 + 2.0 * triang95**2 - 1.2 * triang95**3))
        )

    @staticmethod
    def calculate_current_coefficient_todd(
        eps: float, kappa95: float, triang95: float, model: int
    ) -> float:
        """
        Calculate the fq coefficient used in the two Todd plasma current scalings.

        Parameters
        ----------
        - eps: float, plasma inverse aspect ratio
        - kappa95: float, plasma elongation 95%
        - triang95: float, plasma triangularity 95%

        Returns
        -------
        - float, the fq plasma current coefficient

        Raises
        ------
        ProcessValueError
            If model is not 1 or 2

        This function calculates the fq coefficient based on the given plasma parameters
        for the two Todd scalings.

        References
        ----------
        - D.C.Robinson and T.N.Todd, Plasma and Contr Fusion 28 (1986) 1181
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study',
          AEA FUS 172, 1992
        """
        # Calculate the Todd scaling based on the model
        base_scaling = (
            (1.0 + 2.0 * eps**2)
            * ((1.0 + kappa95**2) / 2)
            * (
                1.24
                - 0.54 * kappa95
                + 0.3 * (kappa95**2 + triang95**2)
                + 0.125 * triang95
            )
        )
        if model == 1:
            return base_scaling
        if model == 2:
            return base_scaling * (1.0 + (abs(kappa95 - 1.2)) ** 3)
        raise ProcessValueError(f"model = {model} is an invalid option")

    @staticmethod
    def calculate_current_coefficient_hastie(
        alphaj: float,
        alphap: float,
        b_plasma_toroidal_on_axis: float,
        delta95: float,
        eps: float,
        kappa95: float,
        pres_plasma_on_axis: float,
        rmu0: float,
    ) -> float:
        """
        Routine to calculate the f_q coefficient for the Connor-Hastie model used for
        scaling the plasma current.

        Parameters
        ----------
        - alphaj: float, the current profile index
        - alphap: float, the pressure profile index
        - b_plasma_toroidal_on_axis: float, the toroidal field on axis (T)
        - delta95: float, the plasma triangularity 95%
        - eps: float, the inverse aspect ratio
        - kappa95: float, the plasma elongation 95%
        - pres_plasma_on_axis: float, the central plasma pressure (Pa)
        - rmu0: float, the vacuum permeability (H/m)

        Returns
        -------
        - float, the F coefficient

        This routine calculates the f_q coefficient used for scaling the plasma current,
        using the Connor-Hastie scaling

        Reference:
        - J.W.Connor and R.J.Hastie, Culham Lab Report CLM-M106 (1985).
          https://scientific-publications.ukaea.uk/wp-content/uploads/CLM-M106-1.pdf
        - T.C.Hender et.al., 'Physics Assesment of the European Reactor Study',
          AEA FUS 172, 1992
        """
        # Exponent in Connor-Hastie current profile
        lamda = alphaj

        # Exponent in Connor-Hastie pressure profile
        nu = alphap

        # Central plasma beta
        beta0 = 2.0 * rmu0 * pres_plasma_on_axis / (b_plasma_toroidal_on_axis**2)

        # Plasma internal inductance
        lamp1 = 1.0 + lamda
        li = lamp1 / lamda * (lamp1 / lamda * np.log(lamp1) - 1.0)

        # T/r in AEA FUS 172
        kap1 = kappa95 + 1.0
        tr = kappa95 * delta95 / kap1**2

        # E/r in AEA FUS 172
        er = (kappa95 - 1.0) / kap1

        # T primed in AEA FUS 172
        tprime = 2.0 * tr * lamp1 / (1.0 + 0.5 * lamda)

        # E primed in AEA FUS 172
        eprime = er * lamp1 / (1.0 + lamda / 3.0)

        # Delta primed in AEA FUS 172
        deltap = (0.5 * kap1 * eps * 0.5 * li) + (
            beta0 / (0.5 * kap1 * eps)
        ) * lamp1**2 / (1.0 + nu)

        # Delta/R0 in AEA FUS 172
        deltar = beta0 / 6.0 * (1.0 + 5.0 * lamda / 6.0 + 0.25 * lamda**2) + (
            0.5 * kap1 * eps
        ) ** 2 * 0.125 * (1.0 - (lamda**2) / 3.0)

        # F coefficient
        return (0.5 * kap1) ** 2 * (
            1.0
            + eps**2 * (0.5 * kap1) ** 2
            + 0.5 * deltap**2
            + 2.0 * deltar
            + 0.5 * (eprime**2 + er**2)
            + 0.5 * (tprime**2 + 4.0 * tr**2)
        )

    @staticmethod
    def calculate_current_coefficient_sauter(
        eps: float,
        kappa: float,
        triang: float,
    ) -> float:
        """
        Routine to calculate the f_q coefficient for the Sauter model used for scaling
        the plasma current.

        Parameters
        ----------
        - eps: float, inverse aspect ratio
        - kappa: float, plasma elongation at the separatrix
        - triang: float, plasma triangularity at the separatrix

        Returns
        -------
        - float, the fq coefficient

        Reference:
        - O. Sauter, Geometric formulas for system codes including the effect of
        negative triangularity, Fusion Engineering and Design, Volume 112, 2016,
        Pages 633-645, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2016.04.033.
        """
        w07 = 1.0  # zero squareness - can be modified later if required

        return (
            (4.1e6 / 5.0e6)
            * (1.0 + 1.2 * (kappa - 1.0) + 0.56 * (kappa - 1.0) ** 2)
            * (1.0 + 0.09 * triang + 0.16 * triang**2)
            * (1.0 + 0.45 * triang * eps)
            / (1.0 - 0.74 * eps)
            * (1.0 + 0.55 * (w07 - 1.0))
        )

    @staticmethod
    def calculate_current_coefficient_fiesta(
        eps: float, kappa: float, triang: float
    ) -> float:
        """
        Calculate the fq coefficient used in the FIESTA plasma current scaling.

        Parameters
        ----------
        - eps: float, plasma inverse aspect ratio
        - kappa: float, plasma elongation at the separatrix
        - triang: float, plasma triangularity at the separatrix

        Returns
        -------
        - float, the fq plasma current coefficient

        This function calculates the fq coefficient based on the given plasma parameters
        for the FIESTA scaling.

        References
        ----------
        - S.Muldrew et.al,“PROCESS”: Systems studies of spherical tokamaks,
          Fusion Engineering and Design, Volume 154, 2020, 111530, ISSN 0920-3796,
          https://doi.org/10.1016/j.fusengdes.2020.111530.
        """
        return 0.538 * (1.0 + 2.440 * eps**2.736) * kappa**2.154 * triang**0.060


class PlasmaDiamagneticCurrentModel(IntEnum):
    """Enum for plasma diamagnetic current method types"""

    NONE = (0, "None")
    HENDER_ST_FIT = (1, "Hender ST fit")
    SCENE_FIT = (2, "SCENE fit")

    def __new__(cls, value: int, full_name: str):
        """Create a new enum member with value and full name.

        Parameters
        ----------
        value :
            The enum value
        full_name :
            The full name description of the enum member

        Returns
        -------
        PlasmaDiamagneticCurrentModel
            The new enum member
        """
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj._full_name_ = full_name
        return obj

    @DynamicClassAttribute
    def full_name(self):
        """Get the full name of the enum member.

        Returns
        -------
        str
            The full name description
        """
        return self._full_name_


class PlasmaDiamagneticCurrent(Model):
    """Class to hold plasma diamagnetic current calculations for plasma processing."""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def run(self):
        """Calculate plasma diamagnetic current fractions using scalings."""
        # Hender scaling for diamagnetic current at tight physics_variables.aspect ratio
        current_drive_variables.f_c_plasma_diamagnetic_hender = (
            self.diamagnetic_fraction_hender(physics_variables.beta_total_vol_avg)
        )

        # SCENE scaling for diamagnetic current
        current_drive_variables.f_c_plasma_diamagnetic_scene = (
            self.diamagnetic_fraction_scene(
                physics_variables.beta_total_vol_avg,
                physics_variables.q95,
                physics_variables.q0,
            )
        )

        if (
            physics_variables.i_diamagnetic_current
            == PlasmaDiamagneticCurrentModel.HENDER_ST_FIT
        ):
            current_drive_variables.f_c_plasma_diamagnetic = (
                current_drive_variables.f_c_plasma_diamagnetic_hender
            )
        elif (
            physics_variables.i_diamagnetic_current
            == PlasmaDiamagneticCurrentModel.SCENE_FIT
        ):
            current_drive_variables.f_c_plasma_diamagnetic = (
                current_drive_variables.f_c_plasma_diamagnetic_scene
            )

    def output(self):
        """Output the plasma diamagnetic current model results."""
        po.oblnkl(self.outfile)
        po.ovarin(
            self.outfile,
            "Plasma diamagnetic current fraction scaling law used",
            "(i_diamagnetic_current)",
            physics_variables.i_diamagnetic_current,
        )
        full_model_name = PlasmaDiamagneticCurrentModel(
            physics_variables.i_diamagnetic_current
        ).full_name
        po.ocmmnt(
            self.outfile,
            f"Diamagnetic current fraction model selected: {full_model_name} ",
        )
        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile,
            "Diamagnetic fraction (Hender)",
            "(f_c_plasma_diamagnetic_hender)",
            current_drive_variables.f_c_plasma_diamagnetic_hender,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Diamagnetic fraction (SCENE)",
            "(f_c_plasma_diamagnetic_scene)",
            current_drive_variables.f_c_plasma_diamagnetic_scene,
            "OP ",
        )
        po.oblnkl(self.outfile)
        if (
            physics_variables.i_diamagnetic_current == PlasmaDiamagneticCurrentModel.NONE
            and current_drive_variables.f_c_plasma_diamagnetic_scene > 0.01e0
        ):
            # Error to show if diamagnetic current is above 1% but not used
            logger.error(
                "Diamagnetic fraction is more than 1%, but not calculated. "
                "Consider using i_diamagnetic_current=2 and i_pfirsch_schluter_current=1"
            )

    @staticmethod
    @nb.njit(cache=True)
    def diamagnetic_fraction_hender(beta: float) -> float:
        """Calculate the diamagnetic fraction based on the Hender fit.

        Parameters
        ----------
        beta :
            the plasma beta value

        Returns
        -------
        float
            the diamagnetic fraction

        """
        return beta / 2.8

    @staticmethod
    @nb.njit(cache=True)
    def diamagnetic_fraction_scene(beta: float, q95: float, q0: float) -> float:
        """Calculate the diamagnetic fraction based on the SCENE fit by Tim Hender.

        Parameters
        ----------
        beta :
            the plasma beta value
        q95 :
            the normalized safety factor at 95% of the plasma radius
        q0 :
            the normalized safety factor at the magnetic axis

        Returns
        -------
        float
            the diamagnetic fraction
        """
        return beta * (0.1 * q95 / q0 + 0.44) * 0.414
