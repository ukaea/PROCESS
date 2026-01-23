import logging

from process.blanket_library import dshellarea, dshellvol, eshellarea, eshellvol
from process.data_structure import blanket_library as blanket_library
from process.data_structure import build_variables as build_variables
from process.data_structure import ccfe_hcpb_module as ccfe_hcpb_module
from process.data_structure import divertor_variables as divertor_variables
from process.data_structure import fwbs_variables as fwbs_variables
from process.data_structure import physics_variables as physics_variables

logger = logging.getLogger(__name__)


class Shield:
    def __init__(self) -> None:
        pass

    def run(self) -> None:
        blanket_library.dz_shld_half = self.calculate_shield_half_height(
            z_plasma_xpoint_lower=build_variables.z_plasma_xpoint_lower,
            dz_xpoint_divertor=build_variables.dz_xpoint_divertor,
            dz_divertor=divertor_variables.dz_divertor,
            n_divertors=physics_variables.n_divertors,
            z_plasma_xpoint_upper=build_variables.z_plasma_xpoint_upper,
            dr_fw_plasma_gap_inboard=build_variables.dr_fw_plasma_gap_inboard,
            dr_fw_plasma_gap_outboard=build_variables.dr_fw_plasma_gap_outboard,
            dr_fw_inboard=build_variables.dr_fw_inboard,
            dr_fw_outboard=build_variables.dr_fw_outboard,
            dz_blkt_upper=build_variables.dz_blkt_upper,
        )
        # D-shaped blanket and shield
        if physics_variables.itart == 1 or fwbs_variables.i_fw_blkt_vv_shape == 1:
            (
                blanket_library.a_shld_inboard_surface,
                blanket_library.a_shld_outboard_surface,
                blanket_library.a_shld_total_surface,
            ) = self.calculate_dshaped_shield_areas(
                rsldi=build_variables.rsldi,
                dr_shld_inboard=build_variables.dr_shld_inboard,
                dr_fw_inboard=build_variables.dr_fw_inboard,
                dr_fw_plasma_gap_inboard=build_variables.dr_fw_plasma_gap_inboard,
                rminor=physics_variables.rminor,
                dr_fw_plasma_gap_outboard=build_variables.dr_fw_plasma_gap_outboard,
                dr_fw_outboard=build_variables.dr_fw_outboard,
                dr_blkt_inboard=build_variables.dr_blkt_inboard,
                dr_blkt_outboard=build_variables.dr_blkt_outboard,
                dz_shld_half=blanket_library.dz_shld_half,
            )

            (
                blanket_library.vol_shld_inboard,
                blanket_library.vol_shld_outboard,
                blanket_library.vol_shld_total,
            ) = self.calculate_dshaped_shield_volumes(
                rsldi=build_variables.rsldi,
                dr_shld_inboard=build_variables.dr_shld_inboard,
                dr_fw_inboard=build_variables.dr_fw_inboard,
                dr_fw_plasma_gap_inboard=build_variables.dr_fw_plasma_gap_inboard,
                rminor=physics_variables.rminor,
                dr_fw_plasma_gap_outboard=build_variables.dr_fw_plasma_gap_outboard,
                dr_fw_outboard=build_variables.dr_fw_outboard,
                dr_blkt_inboard=build_variables.dr_blkt_inboard,
                dr_blkt_outboard=build_variables.dr_blkt_outboard,
                dz_shld_half=blanket_library.dz_shld_half,
                dr_shld_outboard=build_variables.dr_shld_outboard,
                dz_shld_upper=build_variables.dz_shld_upper,
            )

        else:
            (
                build_variables.a_shld_inboard_surface,
                build_variables.a_shld_outboard_surface,
                build_variables.a_shld_total_surface,
            ) = self.calculate_elliptical_shield_areas(
                rsldi=build_variables.rsldi,
                rsldo=build_variables.rsldo,
                rmajor=physics_variables.rmajor,
                triang=physics_variables.triang,
                dr_shld_inboard=build_variables.dr_shld_inboard,
                rminor=physics_variables.rminor,
                dz_shld_half=blanket_library.dz_shld_half,
                dr_shld_outboard=build_variables.dr_shld_outboard,
            )

            (
                build_variables.vol_shld_inboard,
                build_variables.vol_shld_outboard,
                build_variables.vol_shld_total,
            ) = self.calculate_elliptical_shield_volumes(
                rsldi=build_variables.rsldi,
                rsldo=build_variables.rsldo,
                rmajor=physics_variables.rmajor,
                triang=physics_variables.triang,
                dr_shld_inboard=build_variables.dr_shld_inboard,
                rminor=physics_variables.rminor,
                dz_shld_half=blanket_library.dz_shld_half,
                dr_shld_outboard=build_variables.dr_shld_outboard,
                dz_shld_upper=build_variables.dz_shld_upper,
            )

    @staticmethod
    def calculate_shield_half_height(
        z_plasma_xpoint_lower: float,
        dz_xpoint_divertor: float,
        dz_divertor: float,
        n_divertors: int,
        z_plasma_xpoint_upper: float,
        dr_fw_plasma_gap_inboard: float,
        dr_fw_plasma_gap_outboard: float,
        dr_fw_inboard: float,
        dr_fw_outboard: float,
        dz_blkt_upper: float,
    ) -> float:
        """Calculate shield half-height."""

        z_bottom = z_plasma_xpoint_lower + dz_xpoint_divertor + dz_divertor

        # Calculate component internal upper half-height (m)
        # If a double null machine then symmetric
        if n_divertors == 2:
            z_top = z_bottom
        else:
            z_top = z_plasma_xpoint_upper + 0.5 * (
                dr_fw_plasma_gap_inboard
                + dr_fw_plasma_gap_outboard
                + dr_fw_inboard
                + dr_fw_outboard
            )

            z_top = z_top + dz_blkt_upper

        # Average of top and bottom (m)
        return 0.5 * (z_top + z_bottom)

    @staticmethod
    def calculate_dshaped_shield_volumes(
        rsldi: float,
        dr_shld_inboard: float,
        dr_fw_inboard: float,
        dr_fw_plasma_gap_inboard: float,
        rminor: float,
        dr_fw_plasma_gap_outboard: float,
        dr_fw_outboard: float,
        dr_blkt_inboard: float,
        dr_blkt_outboard: float,
        dz_shld_half: float,
        dr_shld_outboard: float,
        dz_shld_upper: float,
    ) -> tuple[float, float, float]:
        """Calculate volumes of D-shaped shield segments."""

        r_1 = rsldi + dr_shld_inboard
        r_2 = (
            dr_fw_inboard
            + dr_fw_plasma_gap_inboard
            + 2.0 * rminor
            + dr_fw_plasma_gap_outboard
            + dr_fw_outboard
        )

        r_2 = dr_blkt_inboard + r_2 + dr_blkt_outboard

        (
            vol_shld_inboard,
            vol_shld_outboard,
            vol_shld_total,
        ) = dshellvol(
            rmajor=r_1,
            rminor=r_2,
            zminor=dz_shld_half,
            drin=dr_shld_inboard,
            drout=dr_shld_outboard,
            dz=dz_shld_upper,
        )

        return vol_shld_inboard, vol_shld_outboard, vol_shld_total

    @staticmethod
    def calculate_dshaped_shield_areas(
        rsldi: float,
        dr_shld_inboard: float,
        dr_fw_inboard: float,
        dr_fw_plasma_gap_inboard: float,
        rminor: float,
        dr_fw_plasma_gap_outboard: float,
        dr_fw_outboard: float,
        dr_blkt_inboard: float,
        dr_blkt_outboard: float,
        dz_shld_half: float,
    ) -> tuple[float, float, float]:
        """Calculate areas of D-shaped shield segments."""

        r_1 = rsldi + dr_shld_inboard
        r_2 = (
            dr_fw_inboard
            + dr_fw_plasma_gap_inboard
            + 2.0 * rminor
            + dr_fw_plasma_gap_outboard
            + dr_fw_outboard
        )

        r_2 = dr_blkt_inboard + r_2 + dr_blkt_outboard

        (
            a_shld_inboard_surface,
            a_shld_outboard_surface,
            a_shld_total_surface,
        ) = dshellarea(rmajor=r_1, rminor=r_2, zminor=dz_shld_half)

        return a_shld_inboard_surface, a_shld_outboard_surface, a_shld_total_surface

    @staticmethod
    def calculate_elliptical_shield_volumes(
        rsldi: float,
        rsldo: float,
        rmajor: float,
        triang: float,
        dr_shld_inboard: float,
        rminor: float,
        dz_shld_half: float,
        dr_shld_outboard: float,
        dz_shld_upper: float,
    ) -> tuple[float, float, float]:
        """Calculate volumes of elliptical shield segments."""

        # Major radius to centre of inboard and outboard ellipses (m)
        # (coincident in radius with top of plasma)
        r_1 = rmajor - rminor * triang
        r_2 = r_1 - rsldi

        r_2 = r_2 - dr_shld_inboard

        r_3 = rsldo - r_1
        r_3 = r_3 - dr_shld_outboard

        (
            vol_shld_inboard,
            vol_shld_outboard,
            vol_shld_total,
        ) = eshellvol(
            rshell=r_1,
            rmini=r_2,
            rmino=r_3,
            zminor=dz_shld_half,
            drin=dr_shld_inboard,
            drout=dr_shld_outboard,
            dz=dz_shld_upper,
        )

        return vol_shld_inboard, vol_shld_outboard, vol_shld_total

    @staticmethod
    def calculate_elliptical_shield_areas(
        rsldi: float,
        rsldo: float,
        rmajor: float,
        triang: float,
        dr_shld_inboard: float,
        rminor: float,
        dz_shld_half: float,
        dr_shld_outboard: float,
    ) -> tuple[float, float, float]:
        """Calculate areas of elliptical shield segments."""

        # Major radius to centre of inboard and outboard ellipses (m)
        # (coincident in radius with top of plasma)
        r_1 = rmajor - rminor * triang
        r_2 = r_1 - rsldi

        r_2 = r_2 - dr_shld_inboard

        r_3 = rsldo - r_1
        r_3 = r_3 - dr_shld_outboard

        (
            a_shld_inboard_surface,
            a_shld_outboard_surface,
            a_shld_total_surface,
        ) = eshellarea(rshell=r_1, rmini=r_2, rmino=r_3, zminor=dz_shld_half)

        return a_shld_inboard_surface, a_shld_outboard_surface, a_shld_total_surface
