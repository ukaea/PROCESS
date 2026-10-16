"""Module for calculating radiation power densities in tokamak plasmas.

This module provides functions to compute synchrotron and impurity radiation
power densities for plasma systems.
"""

import logging
from dataclasses import dataclass

import numpy as np

import process.models.physics.impurity_radiation as impurity
from process.core import constants
from process.core import process_output as po
from process.core.data_structure.base import DataStructure
from process.core.model import Model
from process.data_structure.stellarator_variables import StellaratorModel
from process.models.physics.plasma_profiles import PlasmaProfile

logger = logging.getLogger(__name__)


@dataclass
class RadpwrData:
    """DataClass which holds the output of the function radpwr"""

    pden_plasma_sync_vol_avg_mw: float
    "Volume-averaged synchrotron radiation power density [MW/m³]"
    pden_plasma_core_rad_vol_avg_mw: float
    "Volume-averaged core radiation power density [MW/m³]"
    pden_plasma_core_rad_reduced_vol_avg_mw: float
    "Volume-averaged reduced core radiation power density [MW/m³]"
    pden_plasma_edge_rad_vol_avg_mw: float
    "Volume-averaged edge radiation power density [MW/m³]"
    pden_plasma_rad_vol_avg_mw: float
    "Volume-averaged total radiation power density [MW/m³]"


class PlasmaRadiation(Model):
    """Class representing the radiation properties of a plasma."""

    def __init__(
        self,
        plasma_profile: PlasmaProfile,
        data_structure: DataStructure,
    ):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE
        self.plasma_profile = plasma_profile
        self.data = data_structure
        self.impurity_rad = impurity.ImpurityRadiation(
            plasma_profile=self.plasma_profile, data_structure=self.data
        )

    def run(self) -> None:
        """Run the plasma radiation calculation and update the data structure."""
        radpwrdata = self.calculate_radiation_powers(
            nd_plasma_electron_on_axis=self.data.physics.nd_plasma_electron_on_axis,
            rminor=self.data.physics.rminor,
            b_plasma_toroidal_on_axis=self.data.physics.b_plasma_toroidal_on_axis,
            aspect=self.data.physics.aspect,
            alphan=self.data.physics.alphan,
            alphat=self.data.physics.alphat,
            tbeta=self.data.physics.tbeta,
            temp_plasma_electron_on_axis_kev=self.data.physics.temp_plasma_electron_on_axis_kev,
            f_sync_reflect=self.data.physics.f_sync_reflect,
            rmajor=self.data.physics.rmajor,
            kappa=self.data.physics.kappa,
            vol_plasma=self.data.physics.vol_plasma,
        )
        self.data.physics.pden_plasma_sync_vol_avg_mw = (
            radpwrdata.pden_plasma_sync_vol_avg_mw
        )
        self.data.physics.pden_plasma_core_rad_reduced_vol_avg_mw = (
            radpwrdata.pden_plasma_core_rad_reduced_vol_avg_mw
        )
        self.data.physics.pden_plasma_core_rad_vol_avg_mw = (
            radpwrdata.pden_plasma_core_rad_vol_avg_mw
        )
        self.data.physics.pden_plasma_edge_rad_vol_avg_mw = (
            radpwrdata.pden_plasma_edge_rad_vol_avg_mw
        )
        self.data.physics.pden_plasma_rad_vol_avg_mw = (
            radpwrdata.pden_plasma_rad_vol_avg_mw
        )

        self.data.physics.p_plasma_sync_mw = (
            self.data.physics.pden_plasma_sync_vol_avg_mw * self.data.physics.vol_plasma
        )
        self.data.physics.p_plasma_core_rad_reduced_mw = (
            self.data.physics.pden_plasma_core_rad_reduced_vol_avg_mw
            * self.data.physics.vol_plasma
        )

        self.data.physics.p_plasma_core_rad_mw = (
            self.data.physics.pden_plasma_core_rad_vol_avg_mw
            * self.data.physics.vol_plasma
        )
        self.data.physics.p_plasma_edge_rad_mw = (
            self.data.physics.pden_plasma_edge_rad_vol_avg_mw
            * self.data.physics.vol_plasma
        )
        self.data.physics.p_plasma_rad_mw = (
            self.data.physics.pden_plasma_rad_vol_avg_mw * self.data.physics.vol_plasma
        )

    def output(self) -> None:
        """Subroutine to output the plasma radiation information"""
        po.oheadr(self.outfile, "Plasma Radiation (excluding SOL):")

        po.osubhd(self.outfile, "Synchrotron Radiation:")

        po.ovarre(
            self.outfile,
            "Plasma total synchrotron radiation power (Pₛₙ) [MW]",
            "(p_plasma_sync_mw)",
            self.data.physics.p_plasma_sync_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma volume-averaged synchrotron radiation power density [MW/m³]",
            "(pden_plasma_sync_vol_avg_mw)",
            self.data.physics.pden_plasma_sync_vol_avg_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Synchrotron wall reflectivity factor",
            "(f_sync_reflect)",
            self.data.physics.f_sync_reflect,
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Core region:")

        po.ocmmnt(
            self.outfile, "Note: Only the core region contains the synchrotron radiation"
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma normalised minor radius defining 'core' region (ρᵢₙₙₑᵣ)",
            "(radius_plasma_core_norm)",
            self.data.impurity_radiation.radius_plasma_core_norm,
        )
        po.ovarre(
            self.outfile,
            "Fractional scaling of radiation power along core profile",
            "(f_p_plasma_core_rad_reduction)",
            self.data.impurity_radiation.f_p_plasma_core_rad_reduction,
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma reduced core radiation power (Pᵧ,ᵢₙₙₑᵣ) [MW]",
            "(p_plasma_core_rad_reduced_mw)",
            self.data.physics.p_plasma_core_rad_reduced_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Plasma core radiation power (Pᵧ,ᵢₙₙₑᵣ) [MW]",
            "(p_plasma_core_rad_mw)",
            self.data.physics.p_plasma_core_rad_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(
            self.outfile,
            "Note: Volume-averaged quantities are calculated over the whole plasma",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Volume-averaged total reduced core radiation power per volume [MW/m³]",
            "(pden_plasma_core_rad_reduced_vol_avg_mw)",
            self.data.physics.pden_plasma_core_rad_reduced_vol_avg_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Volume-averaged total core radiation power per volume [MW/m³]",
            "(pden_plasma_core_rad_vol_avg_mw)",
            self.data.physics.pden_plasma_core_rad_vol_avg_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Edge region:")
        po.ovarre(
            self.outfile,
            "Plasma total radiation power from edge region (Pᵧ,ₒᵤₜₑᵣ) [MW] ",
            "(p_plasma_edge_rad_mw)",
            self.data.physics.p_plasma_edge_rad_mw,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ocmmnt(
            self.outfile,
            "Note: Volume-averaged quantities are calculated over the whole plasma",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Volume-averaged total edge radiation power per volume [MW/m³]",
            "(pden_plasma_edge_rad_vol_avg_mw)",
            self.data.physics.pden_plasma_edge_rad_vol_avg_mw,
            "OP ",
        )

        if self.data.stellarator.istell != StellaratorModel.DISABLED:
            po.ovarre(
                self.outfile,
                "SOL radiation power as imposed by f_rad (MW)",
                "(psolradmw)",
                self.data.physics.psolradmw,
                "OP ",
            )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")
        po.osubhd(self.outfile, "Totals:")

        po.ovarre(
            self.outfile,
            "Plasma total radiation power from inside last closed flux surface "
            "(Pᵧ) [MW]",
            "(p_plasma_rad_mw)",
            self.data.physics.p_plasma_rad_mw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Separatrix radiation fraction (fᵧ)",
            "(f_p_plasma_separatrix_rad)",
            self.data.physics.f_p_plasma_separatrix_rad,
            "OP ",
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Plasma total volume-averaged radiation power per unit volume [MW/m³]",
            "(pden_plasma_rad_vol_avg_mw)",
            self.data.physics.pden_plasma_rad_vol_avg_mw,
            "OP ",
        )

        po.oblnkl(self.outfile)
        po.ocmmnt(self.outfile, "----------------------------")

    def calculate_radiation_powers(
        self,
        nd_plasma_electron_on_axis: float,
        rminor: float,
        b_plasma_toroidal_on_axis: float,
        aspect: float,
        alphan: float,
        alphat: float,
        tbeta: float,
        temp_plasma_electron_on_axis_kev: float,
        f_sync_reflect: float,
        rmajor: float,
        kappa: float,
        vol_plasma: float,
    ) -> RadpwrData:
        """Calculate the radiation powers in [MW/m³] by calling relevant routines.

        This function computes the radiation power densities for the plasma, including
        impurity radiation and synchrotron radiation. It returns a dataclass containing
        the calculated radiation power densities.

        Parameters
        ----------
        nd_plasma_electron_on_axis : float
            Central electron density [m⁻³].
        rminor : float
            Minor radius of the plasma [m].
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T].
        aspect : float
            Aspect ratio of the plasma [-].
        alphan : float
            Alpha parameter for density profile [-].
        alphat : float
            Alpha parameter for temperature profile [-].
        tbeta : float
            Beta parameter for temperature profile [-].
        temp_plasma_electron_on_axis_kev : float
            Central electron temperature [keV].
        f_sync_reflect : float
            Fraction of synchrotron radiation reflected.
        rmajor : float
            Major radius of the plasma [m].
        kappa : float
            Elongation of the plasma [-].
        vol_plasma : float
            Plasma volume [m³].

        Returns
        -------
        RadpwrData
            A dataclass containing the radiation power densities in [MW/m³].

        References
        ----------
        [1] F. Albajar, J. Johner, and G. Granata, “Improved calculation of synchrotron
        radiation losses in realistic tokamak plasmas,” Nuclear Fusion, vol. 41,
        no. 6, pp. 665-678, Jun. 2001,
        doi: https://doi.org/10.1088/0029-5515/41/6/301.

        [2] I. Fidone, G Giruzzi, and G. Granata, “Synchrotron radiation loss in tokamaks
        of arbitrary geometry,” Nuclear Fusion, vol. 41, no. 12, pp. 1755-1758,
        Dec. 2001, doi: https://doi.org/10.1088/0029-5515/41/12/102.
        """
        imp_rad = self.impurity_rad
        imp_rad.calculate_imprad()

        # Synchrotron radiation power/volume; assumed to be from core only.
        pden_plasma_sync_vol_avg_mw = self.psync_albajar_fidone(
            nd_plasma_electron_on_axis=nd_plasma_electron_on_axis,
            rminor=rminor,
            b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
            aspect=aspect,
            alphan=alphan,
            alphat=alphat,
            tbeta=tbeta,
            temp_plasma_electron_on_axis_kev=temp_plasma_electron_on_axis_kev,
            f_sync_reflect=f_sync_reflect,
            rmajor=rmajor,
            kappa=kappa,
            vol_plasma=vol_plasma,
        )

        # Total core radiation power/volume.
        pden_plasma_core_rad_reduced_vol_avg_mw = (
            imp_rad.pden_impurity_core_rad_reduced_vol_avg_mw
            + pden_plasma_sync_vol_avg_mw
        )

        # Total core radiation power/volume.
        pden_plasma_core_rad_vol_avg_mw = (
            imp_rad.pden_impurity_core_rad_vol_avg_mw + pden_plasma_sync_vol_avg_mw
        )

        # Total radiation power/volume.
        pden_plasma_rad_vol_avg_mw = (
            imp_rad.pden_impurity_rad_total_vol_avg_mw + pden_plasma_sync_vol_avg_mw
        )

        return RadpwrData(
            pden_plasma_sync_vol_avg_mw=pden_plasma_sync_vol_avg_mw,
            pden_plasma_core_rad_vol_avg_mw=pden_plasma_core_rad_vol_avg_mw,
            pden_plasma_core_rad_reduced_vol_avg_mw=pden_plasma_core_rad_reduced_vol_avg_mw,
            pden_plasma_edge_rad_vol_avg_mw=imp_rad.pden_impurity_rad_edge_vol_avg_mw,
            pden_plasma_rad_vol_avg_mw=pden_plasma_rad_vol_avg_mw,
        )

    @staticmethod
    def psync_albajar_fidone(
        nd_plasma_electron_on_axis: float,
        rminor: float,
        b_plasma_toroidal_on_axis: float,
        aspect: float,
        alphan: float,
        alphat: float,
        tbeta: float,
        temp_plasma_electron_on_axis_kev: float,
        f_sync_reflect: float,
        rmajor: float,
        kappa: float,
        vol_plasma: float,
    ) -> float:
        """Calculate the synchrotron radiation power in [MW/m³].

        This function computes the synchrotron radiation power density for the plasma
        based on the plasma shape, major and minor radii, electron density, and
        temperature profiles.

        Parameters
        ----------
        nd_plasma_electron_on_axis : float
            Central electron density [m⁻³].
        rminor : float
            Minor radius of the plasma [m].
        b_plasma_toroidal_on_axis : float
            Toroidal magnetic field [T].
        aspect : float
            Aspect ratio of the plasma [-].
        alphan : float
            Alpha parameter for density profile [-].
        alphat : float
            Alpha parameter for temperature profile [-].
        tbeta : float
            Beta parameter for temperature profile [-].
        temp_plasma_electron_on_axis_kev : float
            Central electron temperature [keV].
        f_sync_reflect : float
            Fraction of synchrotron radiation reflected [-].
        rmajor : float
            Major radius of the plasma [m].
        kappa : float
            Elongation of the plasma [-].
        vol_plasma : float
            Plasma volume [m³].

        Returns
        -------
        float
            Synchrotron radiation power per unit volume [MW/m³].

        References
        ----------
        [1] F. Albajar, J. Johner, and G. Granata, “Improved calculation of synchrotron
        radiation losses in realistic tokamak plasmas,” Nuclear Fusion, vol. 41, no. 6,
        pp. 665-678, Jun. 2001, doi: https://doi.org/10.1088/0029-5515/41/6/301.

        [2] I. Fidone, G Giruzzi, and G. Granata, “Synchrotron radiation loss in
        tokamaks of arbitrary geometry,” Nuclear Fusion, vol. 41, no. 12, pp. 1755-1758,
        Dec. 2001, doi: https://doi.org/10.1088/0029-5515/41/12/102.
        """
        # Variable names are created to closely match those from the reference papers.

        ne0_20 = 1.0e-20 * nd_plasma_electron_on_axis

        p_a0 = 6.04e3 * (rminor * ne0_20) / b_plasma_toroidal_on_axis

        g_function = 0.93 * (1.0 + 0.85 * np.exp(-0.82 * aspect))

        k_function = (
            (alphan + 3.87 * alphat + 1.46) ** -0.79
            * (1.98 + alphat) ** 1.36
            * tbeta**2.14
            * (tbeta**1.53 + 1.87 * alphat - 0.16) ** -1.33
        )

        dum = (
            1.0
            + 0.12
            * (temp_plasma_electron_on_axis_kev / p_a0**0.41)
            * (1.0 - f_sync_reflect) ** 0.41
        ) ** -1.51

        p_sync_mw = (
            3.84e-8
            * (1.0 - f_sync_reflect) ** 0.62
            * rmajor
            * rminor**1.38
            * kappa**0.79
            * b_plasma_toroidal_on_axis**2.62
            * ne0_20**0.38
            * temp_plasma_electron_on_axis_kev
            * (16.0 + temp_plasma_electron_on_axis_kev) ** 2.61
            * dum
            * g_function
            * k_function
        )

        # pden_plasma_sync_vol_avg_mw should be per unit volume; Albajar gives
        # it as total

        return p_sync_mw / vol_plasma
