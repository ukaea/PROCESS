import numpy as np

from process.core import constants
from process.core import process_output as po
from process.core.model import Model
from process.data_structure import (
    pfcoil_variables,
)


class Cryostat(Model):
    def __init__(self):
        self.outfile = constants.NOUT

    def run(self):
        """Run the cryostat calculations.

        This method runs the cryostat calculations, including the calculation of the cryostat geometry.
        """
        # Calculate cryostat geometry
        self.external_cryo_geometry()

    def external_cryo_geometry(self):
        """Calculate cryostat geometry.

        This method calculates the geometry of the cryostat, including the inboard radius,
        the vertical clearance between the uppermost PF coil and the cryostat lid, the half-height
        of the cryostat, the vertical clearance between the TF coil and the cryostat, the cryostat volume,
        the vacuum vessel mass, and the sum of internal vacuum vessel and cryostat masses.
        """
        # Cryostat radius [m]
        # Take radius of furthest PF coil and add clearance
        self.data.fwbs.r_cryostat_inboard = (
            np.max(pfcoil_variables.r_pf_coil_outer) + self.data.fwbs.dr_pf_cryostat
        )

        # Clearance between uppermost PF coil and cryostat lid [m].
        # Scaling from ITER by M. Kovari
        self.data.blanket.dz_pf_cryostat = (
            self.data.build.f_z_cryostat
            * (2.0 * self.data.fwbs.r_cryostat_inboard)
            / 28.440
        )

        # Half-height of cryostat [m]
        # Take height of furthest PF coil and add clearance
        self.data.fwbs.z_cryostat_half_inside = (
            np.max(pfcoil_variables.z_pf_coil_upper) + self.data.blanket.dz_pf_cryostat
        )

        # Vertical clearance between TF coil and cryostat (m)
        self.data.buildings.dz_tf_cryostat = self.data.fwbs.z_cryostat_half_inside - (
            self.data.build.z_tf_inside_half + self.data.build.dr_tf_inboard
        )

        # Internal cryostat space volume [m^3]
        self.data.fwbs.vol_cryostat_internal = (
            np.pi
            * (self.data.fwbs.r_cryostat_inboard) ** 2
            * 2
            * self.data.fwbs.z_cryostat_half_inside
        )

        # Cryostat structure volume [m^3]
        # Calculate by taking the volume of the outer cryostat and subtracting the volume of the inner cryostat
        self.data.fwbs.vol_cryostat = (
            (
                np.pi
                * (self.data.fwbs.r_cryostat_inboard + self.data.build.dr_cryostat) ** 2
            )
            * 2
            * (self.data.build.dr_cryostat + self.data.fwbs.z_cryostat_half_inside)
        ) - (self.data.fwbs.vol_cryostat_internal)

        # Sum of internal vacuum vessel and cryostat masses (kg)
        self.data.fwbs.dewmkg = (
            self.data.fwbs.vol_vv + self.data.fwbs.vol_cryostat
        ) * self.data.fwbs.den_steel

    def output(self):
        """Outputs the cryostat geometry details to the output file."""
        po.oheadr(self.outfile, "Cryostat build")

        po.ovarrf(
            self.outfile,
            "Cryostat thickness (m)",
            "(dr_cryostat)",
            self.data.build.dr_cryostat,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Cryostat internal radius (m)",
            "(r_cryostat_inboard)",
            self.data.fwbs.r_cryostat_inboard,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Cryostat internal half height (m)",
            "(z_cryostat_half_inside)",
            self.data.fwbs.z_cryostat_half_inside,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Vertical clearance from highest PF coil to cryostat (m)",
            "(dz_pf_cryostat)",
            self.data.blanket.dz_pf_cryostat,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Cryostat structure volume (m³)",
            "(vol_cryostat)",
            self.data.fwbs.vol_cryostat,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Cryostat internal volume (m³)",
            "(vol_cryostat_internal)",
            self.data.fwbs.vol_cryostat_internal,
            "OP ",
        )
