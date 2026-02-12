import numpy as np

from process import constants
from process import process_output as po
from process.data_structure import (
    blanket_library,
    build_variables,
    buildings_variables,
    fwbs_variables,
    pfcoil_variables,
)


class Cryostat:
    def __init__(self) -> None:
        self.outfile = constants.NOUT

    def run(self) -> None:
        """Run the cryostat calculations.

        This method runs the cryostat calculations, including the calculation of the cryostat geometry.
        """

        # Calculate cryostat geometry
        self.external_cryo_geometry()

    @staticmethod
    def external_cryo_geometry() -> None:
        """Calculate cryostat geometry.

        This method calculates the geometry of the cryostat, including the inboard radius,
        the vertical clearance between the uppermost PF coil and the cryostat lid, the half-height
        of the cryostat, the vertical clearance between the TF coil and the cryostat, the cryostat volume,
        the vacuum vessel mass, and the sum of internal vacuum vessel and cryostat masses.
        """

        # Cryostat radius [m]
        # Take radius of furthest PF coil and add clearance
        fwbs_variables.r_cryostat_inboard = (
            np.max(pfcoil_variables.r_pf_coil_outer) + fwbs_variables.dr_pf_cryostat
        )

        # Clearance between uppermost PF coil and cryostat lid [m].
        # Scaling from ITER by M. Kovari
        blanket_library.dz_pf_cryostat = (
            build_variables.f_z_cryostat
            * (2.0 * fwbs_variables.r_cryostat_inboard)
            / 28.440
        )

        # Half-height of cryostat [m]
        # Take height of furthest PF coil and add clearance
        fwbs_variables.z_cryostat_half_inside = (
            np.max(pfcoil_variables.z_pf_coil_upper) + blanket_library.dz_pf_cryostat
        )

        # Vertical clearance between TF coil and cryostat (m)
        buildings_variables.dz_tf_cryostat = fwbs_variables.z_cryostat_half_inside - (
            build_variables.z_tf_inside_half + build_variables.dr_tf_inboard
        )

        # Internal cryostat space volume [m^3]
        fwbs_variables.vol_cryostat_internal = (
            np.pi
            * (fwbs_variables.r_cryostat_inboard) ** 2
            * 2
            * fwbs_variables.z_cryostat_half_inside
        )

        # Cryostat structure volume [m^3]
        # Calculate by taking the volume of the outer cryostat and subtracting the volume of the inner cryostat
        fwbs_variables.vol_cryostat = (
            (
                np.pi
                * (fwbs_variables.r_cryostat_inboard + build_variables.dr_cryostat) ** 2
            )
            * 2
            * (build_variables.dr_cryostat + fwbs_variables.z_cryostat_half_inside)
        ) - (fwbs_variables.vol_cryostat_internal)

        # Sum of internal vacuum vessel and cryostat masses (kg)
        fwbs_variables.dewmkg = (
            fwbs_variables.vol_vv + fwbs_variables.vol_cryostat
        ) * fwbs_variables.den_steel

    def cryostat_output(self):
        """Outputs the cryostat geometry details to the output file."""
        po.oheadr(self.outfile, "Cryostat build")

        po.ovarrf(
            self.outfile,
            "Cryostat thickness (m)",
            "(dr_cryostat)",
            build_variables.dr_cryostat,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Cryostat internal radius (m)",
            "(r_cryostat_inboard)",
            fwbs_variables.r_cryostat_inboard,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Cryostat intenral half height (m)",
            "(z_cryostat_half_inside)",
            fwbs_variables.z_cryostat_half_inside,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Vertical clearance from highest PF coil to cryostat (m)",
            "(dz_pf_cryostat)",
            blanket_library.dz_pf_cryostat,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Cryostat structure volume (m^3)",
            "(vol_cryostat)",
            fwbs_variables.vol_cryostat,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Cryostat internal volume (m^3)",
            "(vol_cryostat_internal)",
            fwbs_variables.vol_cryostat_internal,
            "OP ",
        )
