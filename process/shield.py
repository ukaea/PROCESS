import logging

from process.blanket_library import dshellarea, dshellvol

logger = logging.getLogger(__name__)


class Shield:
    def __init__(self) -> None:
        pass

    def run(self) -> None:
        pass

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
