import logging

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
