"""This library contains routines that can be shared by the blanket modules used in PROCESS.

author: G Graham, CCFE, Culham Science Centre
"""

import logging

import numpy as np

from process import constants
from process import process_output as po
from process.coolprop_interface import FluidProperties
from process.data_structure import (
    blanket_library,
    build_variables,
    divertor_variables,
    fwbs_variables,
    heat_transport_variables,
    physics_variables,
    primary_pumping_variables,
)
from process.exceptions import ProcessValueError

logger = logging.getLogger(__name__)

# Acronyms for this module:
# BB          Breeding Blanket
# FW          First Wall
# BZ          Breeder Zone
# MF/BSS      Manifold/Back Supporting Structure
# LT          Low Temperature
# HT          High Temperature
# MMS         Multi Module Segment
# SMS         Single Modle Segment
# IB          Inboard
# OB          Outboard
# HCD         Heating & Current Drive
# FCI         Flow Channel Insert


class BlanketLibrary:
    def __init__(self, fw) -> None:
        self.outfile = constants.NOUT

        self.fw = fw

    def component_volumes(self):
        """Calculate the blanket, shield, vacuum vessel and cryostat volumes
        author: J. Morris, CCFE, Culham Science Centre
        Calculate the blanket, shield, vacuum vessel and cryostat volumes
        """
        # N.B. icomponent is a switch used to specify selected component: blanket=0, sheild=1, vacuum vessel=2
        # Replaced seperate subroutines for blnkt, shld and vv with fuction/subroutine with icomponent switch.

        # Calculate half-height
        # Blanket
        blanket_library.dz_blkt_half = self.component_half_height(icomponent=0)
        # Shield
        blanket_library.dz_shld_half = self.component_half_height(icomponent=1)

        # D-shaped blanket and shield
        if physics_variables.itart == 1 or fwbs_variables.i_fw_blkt_vv_shape == 1:
            for icomponent in range(2):
                self.dshaped_component(icomponent)

        # Elliptical blanket and shield
        else:
            for icomponent in range(2):
                self.elliptical_component(icomponent)

        # Apply coverage factors to volumes and surface areas
        self.apply_coverage_factors()

    def component_half_height(self, icomponent: int):
        """Calculate the blanket, shield or vacuum vessel half-height
        Based on blanket_half_height, shield_half_height, vv_half_height
        original author: J. Morris, CCFE, Culham Science Centre
        author: G. Graham, CCFE, Culham Science Centre
        """
        # Calculate component internal lower half-height (m)
        # Blanket
        if icomponent == 0:
            hbot = (
                build_variables.z_plasma_xpoint_lower
                + build_variables.dz_xpoint_divertor
                + divertor_variables.dz_divertor
                - build_variables.dz_blkt_upper
            )
        # Sheild
        elif icomponent == 1:
            hbot = (
                build_variables.z_plasma_xpoint_lower
                + build_variables.dz_xpoint_divertor
                + divertor_variables.dz_divertor
            )
        # Vacuum vessel
        elif icomponent == 2:
            hbot = (
                build_variables.z_tf_inside_half
                - build_variables.dz_shld_vv_gap
                - build_variables.dz_vv_lower
            )
        else:
            raise ProcessValueError(f"{icomponent=} is invalid, it must be either 0,1")

        # Calculate component internal upper half-height (m)
        # If a double null machine then symmetric
        if divertor_variables.n_divertors == 2:
            htop = hbot
        else:
            # Blanket
            htop = build_variables.z_plasma_xpoint_upper + 0.5 * (
                build_variables.dr_fw_plasma_gap_inboard
                + build_variables.dr_fw_plasma_gap_outboard
                + build_variables.dr_fw_inboard
                + build_variables.dr_fw_outboard
            )
            # Shield
            if icomponent == 1:
                htop = htop + build_variables.dz_blkt_upper

        # Average of top and bottom (m)
        return 0.5 * (htop + hbot)

    def dshaped_component(self, icomponent: int):
        """Calculate component surface area and volume using dshaped scheme
        Based on dshaped_blanket, dshaped_shield, dshaped_vv
        original author: J. Morris, CCFE, Culham Science Centre
        author: G. Graham, CCFE, Culham Science Centre
        """
        # Calculate major radius to outer edge of inboard ...
        # ... section (m)
        r1 = build_variables.rsldi
        # ... shield (m)
        if icomponent == 1:
            r1 = r1 + build_variables.dr_shld_inboard
        # ... blanket (m)
        elif icomponent == 0:
            r1 = r1 + build_variables.dr_shld_inboard + build_variables.dr_blkt_inboard

        # Horizontal distance between inside edges (m)
        # i.e. outer radius of inboard part to inner radius of outboard part
        # Blanket
        r2 = (
            build_variables.dr_fw_inboard
            + build_variables.dr_fw_plasma_gap_inboard
            + 2.0 * physics_variables.rminor
            + build_variables.dr_fw_plasma_gap_outboard
            + build_variables.dr_fw_outboard
        )
        # Sheild
        if icomponent == 1:
            r2 = build_variables.dr_blkt_inboard + r2 + build_variables.dr_blkt_outboard
        # Vaccum Vessel
        if icomponent == 2:
            r2 = build_variables.rsldo - r1

        # Calculate surface area, assuming 100% coverage
        if icomponent == 0:
            (
                build_variables.a_blkt_inboard_surface,
                build_variables.a_blkt_outboard_surface,
                build_variables.a_blkt_total_surface,
            ) = dshellarea(r1, r2, blanket_library.dz_blkt_half)
        if icomponent == 1:
            (
                build_variables.a_shld_inboard_surface,
                build_variables.a_shld_outboard_surface,
                build_variables.a_shld_total_surface,
            ) = dshellarea(r1, r2, blanket_library.dz_shld_half)

        # Calculate volumes, assuming 100% coverage
        if icomponent == 0:
            (
                fwbs_variables.vol_blkt_inboard,
                fwbs_variables.vol_blkt_outboard,
                fwbs_variables.vol_blkt_total,
            ) = dshellvol(
                r1,
                r2,
                blanket_library.dz_blkt_half,
                build_variables.dr_blkt_inboard,
                build_variables.dr_blkt_outboard,
                build_variables.dz_blkt_upper,
            )
        elif icomponent == 1:
            (
                blanket_library.vol_shld_inboard,
                blanket_library.vol_shld_outboard,
                fwbs_variables.vol_shld_total,
            ) = dshellvol(
                r1,
                r2,
                blanket_library.dz_shld_half,
                build_variables.dr_shld_inboard,
                build_variables.dr_shld_outboard,
                build_variables.dz_shld_upper,
            )

    def elliptical_component(self, icomponent: int):
        """Calculate component surface area and volume using elliptical scheme
        Based on elliptical_blanket, elliptical_shield, elliptical_vv
        original author: J. Morris, CCFE, Culham Science Centre
        author: G. Graham, CCFE, Culham Science Centre
        """
        # Major radius to centre of inboard and outboard ellipses (m)
        # (coincident in radius with top of plasma)
        r1 = (
            physics_variables.rmajor
            - physics_variables.rminor * physics_variables.triang
        )

        # Calculate distance between r1 and outer edge of inboard ...
        # ... section (m)
        r2 = r1 - build_variables.rsldi
        # ... shield (m)
        if icomponent == 1:
            r2 = r2 - build_variables.dr_shld_inboard
        # ... blanket (m)
        if icomponent == 0:
            r2 = r2 - build_variables.dr_shld_inboard - build_variables.dr_blkt_inboard

        # Calculate distance between r1 and inner edge of outboard ...
        # ... section (m)
        r3 = build_variables.rsldo - r1
        # ... shield (m)
        if icomponent == 1:
            r3 = r3 - build_variables.dr_shld_outboard
        # ... blanket (m)
        if icomponent == 0:
            r3 = r3 - build_variables.dr_shld_outboard - build_variables.dr_blkt_outboard

        # Calculate surface area, assuming 100% coverage
        if icomponent == 0:
            (
                build_variables.a_blkt_inboard_surface,
                build_variables.a_blkt_outboard_surface,
                build_variables.a_blkt_total_surface,
            ) = eshellarea(r1, r2, r3, blanket_library.dz_blkt_half)
        if icomponent == 1:
            (
                build_variables.a_shld_inboard_surface,
                build_variables.a_shld_outboard_surface,
                build_variables.a_shld_total_surface,
            ) = eshellarea(r1, r2, r3, blanket_library.dz_shld_half)

        # Calculate volumes, assuming 100% coverage
        if icomponent == 0:
            (
                fwbs_variables.vol_blkt_inboard,
                fwbs_variables.vol_blkt_outboard,
                fwbs_variables.vol_blkt_total,
            ) = eshellvol(
                r1,
                r2,
                r3,
                blanket_library.dz_blkt_half,
                build_variables.dr_blkt_inboard,
                build_variables.dr_blkt_outboard,
                build_variables.dz_blkt_upper,
            )
        if icomponent == 1:
            (
                blanket_library.vol_shld_inboard,
                blanket_library.vol_shld_outboard,
                fwbs_variables.vol_shld_total,
            ) = eshellvol(
                r1,
                r2,
                r3,
                blanket_library.dz_shld_half,
                build_variables.dr_shld_inboard,
                build_variables.dr_shld_outboard,
                build_variables.dz_shld_upper,
            )

    def apply_coverage_factors(self):
        """Apply coverage factors to volumes
        author: J. Morris, CCFE, Culham Science Centre
        Apply coverage factors to volumes
        """
        # Apply blanket coverage factors
        if divertor_variables.n_divertors == 2:
            # double null configuration
            build_variables.a_blkt_outboard_surface = (
                build_variables.a_blkt_total_surface
                * (
                    1.0
                    - 2.0 * fwbs_variables.f_ster_div_single
                    - fwbs_variables.f_a_fw_outboard_hcd
                )
                - build_variables.a_blkt_inboard_surface
            )
        else:
            # single null configuration
            build_variables.a_blkt_outboard_surface = (
                build_variables.a_blkt_total_surface
                * (
                    1.0
                    - fwbs_variables.f_ster_div_single
                    - fwbs_variables.f_a_fw_outboard_hcd
                )
                - build_variables.a_blkt_inboard_surface
            )

        build_variables.a_blkt_total_surface = (
            build_variables.a_blkt_inboard_surface
            + build_variables.a_blkt_outboard_surface
        )

        fwbs_variables.vol_blkt_outboard = (
            fwbs_variables.vol_blkt_total
            * (
                1.0
                - fwbs_variables.f_ster_div_single
                - fwbs_variables.f_a_fw_outboard_hcd
            )
            - fwbs_variables.vol_blkt_inboard
        )
        fwbs_variables.vol_blkt_total = (
            fwbs_variables.vol_blkt_inboard + fwbs_variables.vol_blkt_outboard
        )

        # Apply shield coverage factors
        build_variables.a_shld_inboard_surface = (
            fwbs_variables.fvolsi * build_variables.a_shld_inboard_surface
        )
        build_variables.a_shld_outboard_surface = (
            fwbs_variables.fvolso * build_variables.a_shld_outboard_surface
        )
        build_variables.a_shld_total_surface = (
            build_variables.a_shld_inboard_surface
            + build_variables.a_shld_outboard_surface
        )

        blanket_library.vol_shld_inboard = (
            fwbs_variables.fvolsi * blanket_library.vol_shld_inboard
        )
        blanket_library.vol_shld_outboard = (
            fwbs_variables.fvolso * blanket_library.vol_shld_outboard
        )
        fwbs_variables.vol_shld_total = (
            blanket_library.vol_shld_inboard + blanket_library.vol_shld_outboard
        )

    def primary_coolant_properties(self, output: bool):
        """Calculates the fluid properties of the Primary Coolant in the FW and BZ.
        Uses middle value of input and output temperatures of coolant.
        Curently have H20 and He options.

        original author: P. J. Knight, CCFE
        adapted from previous version of pumppower function by: G Graham, CCFE
        References: see pumppower function description
        """

        # Make sure that, if the inputs for the FW and blanket inputs are different,
        # the i_fw_blkt_shared_coolant variable is appropriately set for seperate coolants
        if (
            fwbs_variables.i_fw_coolant_type == "Helium"
            and fwbs_variables.i_blkt_coolant_type == 2
        ):
            fwbs_variables.i_fw_blkt_shared_coolant = 1
        if (
            fwbs_variables.i_fw_coolant_type == "Water"
            and fwbs_variables.i_blkt_coolant_type == 1
        ):
            fwbs_variables.i_fw_blkt_shared_coolant = 1

        # If FW and BB have same coolant...
        if fwbs_variables.i_fw_blkt_shared_coolant == 0:
            # Use FW inlet temp and BB outlet temp
            mid_temp = (
                fwbs_variables.temp_fw_coolant_in + fwbs_variables.temp_blkt_coolant_out
            ) * 0.5
            # FW/BB
            fw_bb_fluid_properties = FluidProperties.of(
                fwbs_variables.i_fw_coolant_type,
                temperature=mid_temp,
                pressure=fwbs_variables.pres_fw_coolant,
            )
            fwbs_variables.den_fw_coolant = fw_bb_fluid_properties.density
            fwbs_variables.cp_fw = fw_bb_fluid_properties.specific_heat_const_p
            fwbs_variables.cv_fw = fw_bb_fluid_properties.specific_heat_const_v
            fwbs_variables.visc_fw_coolant = fw_bb_fluid_properties.viscosity

            fwbs_variables.den_blkt_coolant = fwbs_variables.den_fw_coolant
            fwbs_variables.visc_blkt_coolant = fwbs_variables.visc_fw_coolant
            fwbs_variables.cp_bl = fwbs_variables.cp_fw
            fwbs_variables.cv_bl = fwbs_variables.cv_fw

        # If FW and BB have different coolants...
        else:
            # FW
            mid_temp_fw = (
                fwbs_variables.temp_fw_coolant_in + fwbs_variables.temp_fw_coolant_out
            ) * 0.5
            fw_fluid_properties = FluidProperties.of(
                fwbs_variables.i_fw_coolant_type,
                temperature=mid_temp_fw,
                pressure=fwbs_variables.pres_fw_coolant,
            )
            fwbs_variables.den_fw_coolant = fw_fluid_properties.density
            fwbs_variables.cp_fw = fw_fluid_properties.specific_heat_const_p
            fwbs_variables.cv_fw = fw_fluid_properties.specific_heat_const_v
            fwbs_variables.visc_fw_coolant = fw_fluid_properties.viscosity

            # BB
            mid_temp_bl = (
                fwbs_variables.temp_blkt_coolant_in
                + fwbs_variables.temp_blkt_coolant_out
            ) * 0.5
            bb_fluid_properties = FluidProperties.of(
                "Helium" if fwbs_variables.i_blkt_coolant_type == 1 else "Water",
                temperature=mid_temp_bl,
                pressure=fwbs_variables.pres_blkt_coolant,
            )
            fwbs_variables.den_blkt_coolant = bb_fluid_properties.density
            fwbs_variables.cp_bl = bb_fluid_properties.specific_heat_const_p
            fwbs_variables.cv_bl = bb_fluid_properties.specific_heat_const_v
            fwbs_variables.visc_blkt_coolant = bb_fluid_properties.viscosity

        if (
            fwbs_variables.den_fw_coolant > 1e9
            or fwbs_variables.den_fw_coolant <= 0
            or np.isnan(fwbs_variables.den_fw_coolant)
        ):
            raise ProcessValueError(
                f"Error in primary_coolant_properties. {fwbs_variables.den_fw_coolant = }"
            )
        if (
            fwbs_variables.den_blkt_coolant > 1e9
            or fwbs_variables.den_blkt_coolant <= 0
            or np.isnan(fwbs_variables.den_blkt_coolant)
        ):
            raise ProcessValueError(
                f"Error in primary_coolant_properties. {fwbs_variables.den_blkt_coolant = }"
            )

        if output:
            po.oheadr(
                self.outfile, "First wall and blanket : (Primary) Coolant Properties"
            )
            po.ocmmnt(
                self.outfile,
                "Calculated using mid temp(s) of system (or systems if use different collant types).",
            )

            # FW (or FW/BB)
            if fwbs_variables.i_fw_blkt_shared_coolant == 1:
                po.osubhd(self.outfile, "First Wall :")

            po.ovarst(
                self.outfile,
                "Coolant type",
                "(i_fw_coolant_type)",
                f'"{fwbs_variables.i_fw_coolant_type}"',
            )
            po.ovarrf(
                self.outfile,
                "Density (kg m-3)",
                "(den_fw_coolant)",
                fwbs_variables.den_fw_coolant,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Viscosity (Pa s)",
                "(visc_fw_coolant)",
                fwbs_variables.visc_fw_coolant,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Inlet Temperature (Celcius)",
                "(temp_fw_coolant_in)",
                fwbs_variables.temp_fw_coolant_in,
                "OP ",
            )

            if fwbs_variables.i_fw_blkt_shared_coolant == 0:
                po.ovarre(
                    self.outfile,
                    "Outlet Temperature (Celcius)",
                    "(temp_blkt_coolant_out)",
                    fwbs_variables.temp_blkt_coolant_out,
                    "OP ",
                )

            else:
                po.ovarre(
                    self.outfile,
                    "Outlet Temperature (Celcius)",
                    "(temp_fw_coolant_out)",
                    fwbs_variables.temp_fw_coolant_out,
                    "OP ",
                )

            # BB
            if fwbs_variables.i_fw_blkt_shared_coolant == 1:
                po.osubhd(self.outfile, "Breeding Blanket :")

                if fwbs_variables.i_blkt_coolant_type == 1:
                    po.ocmmnt(
                        self.outfile, "Coolant type (i_blkt_coolant_type=1), Helium"
                    )
                if fwbs_variables.i_blkt_coolant_type == 2:
                    po.ocmmnt(
                        self.outfile, "Coolant type (i_blkt_coolant_type=2), Water"
                    )
                po.ovarrf(
                    self.outfile,
                    "Density (kg m-3)",
                    "(den_blkt_coolant)",
                    fwbs_variables.den_blkt_coolant,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Viscosity (Pa s)",
                    "(visc_blkt_coolant)",
                    fwbs_variables.visc_blkt_coolant,
                    "OP ",
                )

                po.ovarre(
                    self.outfile,
                    "Inlet Temperature (Celcius)",
                    "(temp_blkt_coolant_in)",
                    fwbs_variables.temp_blkt_coolant_in,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Outlet Temperature (Celcius)",
                    "(temp_blkt_coolant_out)",
                    fwbs_variables.temp_blkt_coolant_out,
                    "OP ",
                )

    def set_blanket_module_geometry(self):
        """
        Sets the geometry parameters for blanket modules, including coolant channel dimensions,
        module segmentation, and flow lengths, based on the current configuration and input variables.

        The method performs the following steps:
        - Determines inboard and outboard coolant channel radial lengths based on blanket type.
        - Segments the blanket modules poloidally and toroidally according to input segmentation settings.
        - Calculates the toroidal segment lengths for inboard and outboard blanket modules.
        - Computes the poloidal height of blanket modules.
        - For dual coolant blankets, calculates the minimum available space for liquid breeder pipes
          in radial, toroidal, and poloidal directions, and checks for geometric constraints.
        - Calculates total flow lengths for primary coolant channels, used in pressure drop calculations.

        Raises:
            Error: If the poloidal segment length is less than three times the minimum liquid breeder pipe width.
        """

        if fwbs_variables.i_blanket_type == 5:
            # Unless DCLL then we will use BZ
            blanket_library.len_blkt_inboard_coolant_channel_radial = (
                build_variables.blbuith
            )
            blanket_library.len_blkt_outboard_coolant_channel_radial = (
                build_variables.blbuoth
            )
        else:
            blanket_library.len_blkt_inboard_coolant_channel_radial = (
                0.8e0 * build_variables.dr_blkt_inboard
            )
            blanket_library.len_blkt_outboard_coolant_channel_radial = (
                0.8e0 * build_variables.dr_blkt_outboard
            )

        # Using the total perimeter of the machine, segment the outboard
        # blanket into nblktmodp*nblktmodt modules, all assumed to be the same size

        # If SMS blanket then do not have seperate poloidal modules....
        # Should not need this as n_blkt_inboard_modules_poloidal is input but make sure here.
        if fwbs_variables.i_blkt_module_segmentation == 1:
            fwbs_variables.n_blkt_inboard_modules_poloidal = 1
            fwbs_variables.n_blkt_outboard_modules_poloidal = 1

        # Calculate poloidal height of blanket modules
        self.blanket_module_poloidal_height()

        # If liquid breeder or dual coolant blanket then calculate
        if fwbs_variables.i_blkt_dual_coolant > 0:
            # Use smallest space available to pipes for pipe sizes in pumping calculations (worst case)
            if fwbs_variables.i_blkt_inboard == 1:
                # Radial direction
                fwbs_variables.b_bz_liq = (
                    min(
                        (
                            blanket_library.len_blkt_inboard_coolant_channel_radial
                            * fwbs_variables.r_f_liq_ib
                        ),
                        (
                            blanket_library.len_blkt_outboard_coolant_channel_radial
                            * fwbs_variables.r_f_liq_ob
                        ),
                    )
                    / fwbs_variables.nopol
                )
                # Toroidal direction
                fwbs_variables.a_bz_liq = (
                    min(
                        (
                            blanket_library.len_blkt_inboard_segment_toroidal
                            * fwbs_variables.w_f_liq_ib
                        ),
                        (
                            blanket_library.len_blkt_outboard_segment_toroidal
                            * fwbs_variables.w_f_liq_ob
                        ),
                    )
                    / fwbs_variables.nopipes
                )
                # Poloidal
                if (
                    blanket_library.len_blkt_inboard_segment_poloidal
                    < (fwbs_variables.b_bz_liq * 3)
                ) or (
                    blanket_library.len_blkt_outboard_segment_poloidal
                    < (fwbs_variables.b_bz_liq * 3)
                ):
                    logger.error(
                        "Your blanket modules are too small for the Liquid Metal pipes"
                    )

            # Unless there is no IB blanket...
            else:
                # Radial direction
                fwbs_variables.b_bz_liq = (
                    blanket_library.len_blkt_outboard_coolant_channel_radial
                    * fwbs_variables.r_f_liq_ob
                ) / fwbs_variables.nopol
                # Toroidal direction
                fwbs_variables.a_bz_liq = (
                    blanket_library.len_blkt_outboard_segment_toroidal
                    * fwbs_variables.w_f_liq_ob
                ) / fwbs_variables.nopipes
                # Poloidal
                if blanket_library.len_blkt_outboard_segment_poloidal < (
                    fwbs_variables.b_bz_liq * 3
                ):
                    logger.error(
                        "Your blanket modules are too small for the Liquid Metal pipes"
                    )

        # Calculate total flow lengths, used for pressure drop calculation
        # Blanket primary coolant flow
        blanket_library.len_blkt_inboard_channel_total = (
            fwbs_variables.n_blkt_inboard_module_coolant_sections_radial
            * blanket_library.len_blkt_inboard_coolant_channel_radial
            + fwbs_variables.n_blkt_inboard_module_coolant_sections_poloidal
            * blanket_library.len_blkt_inboard_segment_poloidal
        )
        blanket_library.len_blkt_outboard_channel_total = (
            fwbs_variables.n_blkt_outboard_module_coolant_sections_radial
            * blanket_library.len_blkt_outboard_coolant_channel_radial
            + fwbs_variables.n_blkt_outboard_module_coolant_sections_poloidal
            * blanket_library.len_blkt_outboard_segment_poloidal
        )

    def thermo_hydraulic_model_pressure_drop_calculations(self, output: bool):
        """
        Function that calculates the pressure drops for the thermo-hydraulic model
        when i_p_coolant_pumping = 2.

        Within are calculations necessary for the deltap_tot function but not required
        for other calculations within the thermo-hydraulic model as then they are just
        included there.

        Returns the pressure drops as a list with the number of entries dependent upon
        the switches i_blkt_dual_coolant and i_blkt_inboard.
        """
        npoltoti = 0
        npoltoto = 0
        npblkti_liq = 0
        npblkto_liq = 0

        # Blanket secondary coolant/breeder flow
        pollengi = blanket_library.len_blkt_inboard_segment_poloidal
        pollengo = blanket_library.len_blkt_outboard_segment_poloidal
        fwbs_variables.nopol = 2
        fwbs_variables.nopipes = 4
        bzfllengi_liq = (
            fwbs_variables.bzfllengi_n_rad_liq
            * blanket_library.len_blkt_inboard_coolant_channel_radial
            + fwbs_variables.bzfllengi_n_pol_liq
            * blanket_library.len_blkt_inboard_segment_poloidal
        )
        bzfllengo_liq = (
            fwbs_variables.bzfllengo_n_rad_liq
            * blanket_library.len_blkt_outboard_coolant_channel_radial
            + fwbs_variables.bzfllengo_n_pol_liq
            * blanket_library.len_blkt_outboard_segment_poloidal
        )

        # ======================================================================

        # Coolant channel bends

        # Number of angle turns in FW and blanket flow channels, n.b. these are the
        # same for CCFE HCPB and KIT DCLL. FW is also be the same for DCLL MMS ans SMS.

        N_FW_PIPE_90_DEG_BENDS = 2
        N_FW_PIPE_180_DEG_BENDS = 0

        # N.B. This is for BZ only, does not include MF/BSS.
        if fwbs_variables.i_blkt_dual_coolant in (1, 2):
            N_BLKT_PIPE_90_DEG_BENDS = 4
            N_BLKT_PIPE_180_DEG_BENDS = 1
            no90bz_liq = 2
            no180bz_liq = 1
        else:
            N_BLKT_PIPE_90_DEG_BENDS = 4
            N_BLKT_PIPE_180_DEG_BENDS = 1

        # ======================================================================

        # FW Pipe Flow and Velocity

        # Mass flow rate per FW coolant pipe (kg/s):
        blanket_library.mflow_fw_inboard_coolant_channel = (
            blanket_library.mflow_fw_inboard_coolant_total
            / blanket_library.n_fw_inboard_channels
        )
        blanket_library.mflow_fw_outboard_coolant_channel = (
            blanket_library.mflow_fw_outboard_coolant_total
            / blanket_library.n_fw_outboard_channels
        )

        # Coolant velocity in FW (m/s)
        vel_fw_inboard_coolant = self.flow_velocity(
            i_channel_shape=1,
            mass_flow_rate=blanket_library.mflow_fw_inboard_coolant_channel,
            flow_density=fwbs_variables.den_fw_coolant,
        )
        vel_fw_outboard_coolant = self.flow_velocity(
            i_channel_shape=1,
            mass_flow_rate=blanket_library.mflow_fw_outboard_coolant_channel,
            flow_density=fwbs_variables.den_fw_coolant,
        )

        # If the blanket is dual-coolant...
        if fwbs_variables.i_blkt_dual_coolant == 2:
            # Calc total num of pipes (in all inboard modules) from
            # coolant frac and channel dimensions
            # Assumes up/down flow, two 90 deg bends per length
            blanket_library.n_blkt_outboard_channels = (
                fwbs_variables.f_a_blkt_cooling_channels
                * fwbs_variables.vol_blkt_outboard
            ) / (
                np.pi
                * fwbs_variables.radius_fw_channel
                * fwbs_variables.radius_fw_channel
                * blanket_library.len_blkt_outboard_channel_total
            )
            npblkto_liq = (
                fwbs_variables.nopipes
                * fwbs_variables.n_blkt_outboard_modules_toroidal
                * fwbs_variables.n_blkt_outboard_modules_poloidal
            )

            # Mass flow rate per coolant pipe
            blanket_library.mfblktpo = (
                blanket_library.mflow_blkt_outboard_coolant
                / blanket_library.n_blkt_outboard_channels
            )
            mfblktpo_liq = blanket_library.mfblkto_liq / npblkto_liq
            # Coolant velocites in blanket (m/s)
            # Assume BZ structure has same channel width as FW
            blanket_library.vel_blkt_outboard_coolant = self.flow_velocity(
                i_channel_shape=1,
                mass_flow_rate=blanket_library.mfblktpo,
                flow_density=fwbs_variables.den_blkt_coolant,
            )
            velblkto_liq = self.flow_velocity(
                i_channel_shape=2,
                mass_flow_rate=mfblktpo_liq,
                flow_density=fwbs_variables.den_liq,
            )

            if fwbs_variables.i_blkt_inboard == 1:
                # Calc total num of pipes (in all inboard modules) from
                # coolant frac and channel dimensions
                # Assumes up/down flow, two 90 deg bends per length
                blanket_library.n_blkt_inboard_channels = (
                    fwbs_variables.f_a_blkt_cooling_channels
                    * fwbs_variables.vol_blkt_inboard
                ) / (
                    np.pi
                    * fwbs_variables.radius_fw_channel
                    * fwbs_variables.radius_fw_channel
                    * blanket_library.len_blkt_inboard_channel_total
                )
                # Have DEMO DCLL set here for now
                npblkti_liq = (
                    fwbs_variables.nopipes
                    * fwbs_variables.n_blkt_inboard_modules_toroidal
                    * fwbs_variables.n_blkt_inboard_modules_poloidal
                )

                # Mass flow rate per coolant pipe
                blanket_library.mfblktpi = (
                    blanket_library.mflow_blkt_inboard_coolant
                    / blanket_library.n_blkt_inboard_channels
                )
                blanket_library.mfblktpi_liq = blanket_library.mfblkti_liq / npblkti_liq

                # Coolant velocites in blanket (m/s)
                # Assume BZ structure has same channel width as FW
                blanket_library.vel_blkt_inboard_coolant = self.flow_velocity(
                    i_channel_shape=1,
                    mass_flow_rate=blanket_library.mfblktpi,
                    flow_density=fwbs_variables.den_blkt_coolant,
                )
                velblkti_liq = self.flow_velocity(
                    i_channel_shape=2,
                    mass_flow_rate=blanket_library.mfblktpi_liq,
                    flow_density=fwbs_variables.den_liq,
                )

        # If the blanket is single-coolant with liquid metal breeder...
        elif fwbs_variables.i_blkt_dual_coolant == 1:
            # Calc total num of pipes (in all inboard modules) from
            # coolant frac and channel dimensions
            # Assumes up/down flow, two 90 deg bends per length
            blanket_library.n_blkt_outboard_channels = (
                fwbs_variables.f_a_blkt_cooling_channels
                * fwbs_variables.vol_blkt_outboard
            ) / (
                np.pi
                * fwbs_variables.radius_fw_channel
                * fwbs_variables.radius_fw_channel
                * blanket_library.len_blkt_outboard_channel_total
            )
            npblkto_liq = (
                fwbs_variables.nopipes
                * fwbs_variables.n_blkt_outboard_modules_toroidal
                * fwbs_variables.n_blkt_outboard_modules_poloidal
            )

            # Mass flow rate per coolant pipe
            blanket_library.mfblktpo = (
                blanket_library.mflow_blkt_outboard_coolant
                / blanket_library.n_blkt_outboard_channels
            )

            # Coolant velocity in blanket (m/s)
            # Assume BZ structure has same channel width as FW
            blanket_library.vel_blkt_outboard_coolant = self.flow_velocity(
                i_channel_shape=1,
                mass_flow_rate=blanket_library.mfblktpo,
                flow_density=fwbs_variables.den_blkt_coolant,
            )

            # Get mass flow rate etc. for inboard blanket breeder flow for tritium extraction
            # Use the number of desired recirculations ([Aub2013]=10) and mass from dcll_masses
            # N.B. wht_liq is BZ mass, does not include manifold.
            blanket_library.mfblkto_liq = (
                fwbs_variables.n_liq_recirc * fwbs_variables.wht_liq_ob
            ) / (24 * 3600)
            blanket_library.mfblktpo_liq = blanket_library.mfblkto_liq / npblkto_liq
            velblkto_liq = self.flow_velocity(
                i_channel_shape=2,
                mass_flow_rate=blanket_library.mfblktpo_liq,
                flow_density=fwbs_variables.den_liq,
            )

            if fwbs_variables.i_blkt_inboard == 1:
                # Calc total num of pipes (in all inboard modules) from
                # coolant frac and channel dimensions
                # Assumes up/down flow, two 90 deg bends per length
                blanket_library.n_blkt_inboard_channels = (
                    fwbs_variables.f_a_blkt_cooling_channels
                    * fwbs_variables.vol_blkt_inboard
                ) / (
                    np.pi
                    * fwbs_variables.radius_fw_channel
                    * fwbs_variables.radius_fw_channel
                    * blanket_library.len_blkt_inboard_channel_total
                )
                # Have DEMO DCLL set here for now
                npblkti_liq = (
                    fwbs_variables.nopipes
                    * fwbs_variables.n_blkt_inboard_modules_toroidal
                    * fwbs_variables.n_blkt_inboard_modules_poloidal
                )

                # Mass flow rate per coolant pipe
                blanket_library.mfblktpi = (
                    blanket_library.mflow_blkt_inboard_coolant
                    / blanket_library.n_blkt_inboard_channels
                )

                # Coolant velocity in blanket (m/s)
                # Assume BZ structure has same channel width as FW
                blanket_library.vel_blkt_inboard_coolant = self.flow_velocity(
                    i_channel_shape=1,
                    mass_flow_rate=blanket_library.mfblktpi,
                    flow_density=fwbs_variables.den_blkt_coolant,
                )

                # Get mass flow rate etc. for inboard blanket breeder flow for tritium extraction
                # Use the number of desired recirculations ([Aub2013]=10) and mass from dcll_masses
                # N.B. wht_liq is BZ mass, does not include manifold.
                blanket_library.mfblkti_liq = (
                    fwbs_variables.n_liq_recirc * fwbs_variables.wht_liq_ib
                ) / (24 * 3600)
                blanket_library.mfblktpi_liq = blanket_library.mfblkti_liq / npblkti_liq
                velblkti_liq = self.flow_velocity(
                    i_channel_shape=2,
                    mass_flow_rate=blanket_library.mfblktpi_liq,
                    flow_density=fwbs_variables.den_liq,
                )

        # If the blanket is single-coolant with solid breeder...
        else:
            # Calculate total number of pipes (in all outboard modules) from coolant fraction and
            # channel dimensions (assumes up/down flow, two 90 deg bends per length)
            blanket_library.n_blkt_outboard_channels = (
                fwbs_variables.f_a_blkt_cooling_channels
                * fwbs_variables.vol_blkt_outboard
            ) / (
                np.pi
                * fwbs_variables.radius_fw_channel
                * fwbs_variables.radius_fw_channel
                * blanket_library.len_blkt_outboard_channel_total
            )

            # Mass flow rate per coolant pipe
            blanket_library.mfblktpo = (
                blanket_library.mflow_blkt_outboard_coolant
                / blanket_library.n_blkt_outboard_channels
            )

            # Coolant velocity in blanket (m/s)
            # Assume BZ structure has same channel width as FW
            blanket_library.vel_blkt_outboard_coolant = self.flow_velocity(
                i_channel_shape=1,
                mass_flow_rate=blanket_library.mfblktpo,
                flow_density=fwbs_variables.den_blkt_coolant,
            )

            if fwbs_variables.i_blkt_inboard == 1:
                # Calc total num of pipes (in all inboard modules) from
                # coolant frac and channel dimensions
                # Assumes up/down flow, two 90 deg bends per length
                blanket_library.n_blkt_inboard_channels = (
                    fwbs_variables.f_a_blkt_cooling_channels
                    * fwbs_variables.vol_blkt_inboard
                ) / (
                    np.pi
                    * fwbs_variables.radius_fw_channel
                    * fwbs_variables.radius_fw_channel
                    * blanket_library.len_blkt_inboard_channel_total
                )

                # Mass flow rate per coolant pipe
                blanket_library.mfblktpi = (
                    blanket_library.mflow_blkt_inboard_coolant
                    / blanket_library.n_blkt_inboard_channels
                )

                # Coolant velocity in blanket (m/s)
                # Assume BZ structure has same channel width as FW
                blanket_library.vel_blkt_inboard_coolant = self.flow_velocity(
                    i_channel_shape=1,
                    mass_flow_rate=blanket_library.mfblktpi,
                    flow_density=fwbs_variables.den_blkt_coolant,
                )

        # FW Presure Drops ###############

        (
            fwbs_variables.radius_blkt_channel_90_bend,
            fwbs_variables.radius_blkt_channel_180_bend,
        ) = self.calculate_pipe_bend_radius(i_ps=1)

        dpres_fw_inboard_coolant = self.total_pressure_drop(
            output,
            icoolpump=1,
            vel_coolant=vel_fw_inboard_coolant,
            len_pipe=fwbs_variables.len_fw_channel,
            n_pipe_90_deg_bends=N_FW_PIPE_90_DEG_BENDS,
            n_pipe_180_deg_bends=N_FW_PIPE_180_DEG_BENDS,
            den_coolant=fwbs_variables.den_fw_coolant,
            visc_coolant_dynamic=fwbs_variables.visc_fw_coolant,
            coolant_electrical_conductivity=0.0e0,
            pol_channel_length=pollengi,
            nopolchan=npoltoti,
            label="Inboard first wall",
        )

        dpres_fw_outboard_coolant = self.total_pressure_drop(
            output,
            icoolpump=1,
            vel_coolant=vel_fw_outboard_coolant,
            len_pipe=fwbs_variables.len_fw_channel,
            n_pipe_90_deg_bends=N_FW_PIPE_90_DEG_BENDS,
            n_pipe_180_deg_bends=N_FW_PIPE_180_DEG_BENDS,
            den_coolant=fwbs_variables.den_fw_coolant,
            visc_coolant_dynamic=fwbs_variables.visc_fw_coolant,
            coolant_electrical_conductivity=0.0e0,
            pol_channel_length=pollengo,
            nopolchan=npoltoto,
            label="Outboard first wall",
        )

        # BB Presure Drops ###############
        (
            fwbs_variables.radius_blkt_channel_90_bend,
            fwbs_variables.radius_blkt_channel_180_bend,
        ) = self.calculate_pipe_bend_radius(i_ps=1)

        # Long polodal flows
        if fwbs_variables.i_blkt_inboard == 1:
            npoltoti = fwbs_variables.nopol * npblkti_liq
        npoltoto = fwbs_variables.nopol * npblkto_liq

        dpres_blkt_outboard_coolant = self.total_pressure_drop(
            output,
            icoolpump=1,
            vel_coolant=blanket_library.vel_blkt_outboard_coolant,
            len_pipe=blanket_library.len_blkt_outboard_channel_total,
            n_pipe_90_deg_bends=N_BLKT_PIPE_90_DEG_BENDS,
            n_pipe_180_deg_bends=N_BLKT_PIPE_180_DEG_BENDS,
            den_coolant=fwbs_variables.den_blkt_coolant,
            visc_coolant_dynamic=fwbs_variables.visc_blkt_coolant,
            coolant_electrical_conductivity=0.0e0,
            pol_channel_length=pollengo,
            nopolchan=npoltoto,
            label="Outboard blanket",
        )

        if fwbs_variables.i_blkt_inboard == 1:
            dpres_blkt_inboard_coolant = self.total_pressure_drop(
                output,
                icoolpump=1,
                vel_coolant=blanket_library.vel_blkt_inboard_coolant,
                len_pipe=blanket_library.len_blkt_inboard_channel_total,
                n_pipe_90_deg_bends=N_BLKT_PIPE_90_DEG_BENDS,
                n_pipe_180_deg_bends=N_BLKT_PIPE_180_DEG_BENDS,
                den_coolant=fwbs_variables.den_blkt_coolant,
                visc_coolant_dynamic=fwbs_variables.visc_blkt_coolant,
                coolant_electrical_conductivity=0.0e0,
                pol_channel_length=pollengi,
                nopolchan=npoltoti,
                label="Inboard blanket",
            )

        # If the blanket has a liquid metal breeder...
        if fwbs_variables.i_blkt_dual_coolant > 0:
            deltap_blo_liq = self.total_pressure_drop(
                output,
                icoolpump=2,
                vel_coolant=velblkto_liq,
                len_pipe=bzfllengo_liq,
                n_pipe_90_deg_bends=no90bz_liq,
                n_pipe_180_deg_bends=no180bz_liq,
                den_coolant=fwbs_variables.den_liq,
                visc_coolant_dynamic=fwbs_variables.dynamic_viscosity_liq,
                coolant_electrical_conductivity=fwbs_variables.electrical_conductivity_liq,
                pol_channel_length=pollengo,
                nopolchan=npoltoto,
                label="Outboard blanket breeder liquid",
            )
            if fwbs_variables.i_blkt_inboard == 1:
                deltap_bli_liq = self.total_pressure_drop(
                    output,
                    icoolpump=2,
                    vel_coolant=velblkti_liq,
                    len_pipe=bzfllengi_liq,
                    n_pipe_90_deg_bends=no90bz_liq,
                    n_pipe_180_deg_bends=no180bz_liq,
                    den_coolant=fwbs_variables.den_liq,
                    visc_coolant_dynamic=fwbs_variables.dynamic_viscosity_liq,
                    coolant_electrical_conductivity=fwbs_variables.electrical_conductivity_liq,
                    pol_channel_length=pollengi,
                    nopolchan=npoltoti,
                    label="Inboard blanket breeder liquid",
                )

                return [
                    dpres_fw_inboard_coolant,
                    dpres_fw_outboard_coolant,
                    dpres_blkt_outboard_coolant,
                    dpres_blkt_inboard_coolant,
                    deltap_blo_liq,
                    deltap_bli_liq,
                ]
            return [
                dpres_fw_inboard_coolant,
                dpres_fw_outboard_coolant,
                dpres_blkt_outboard_coolant,
                deltap_blo_liq,
            ]

        if fwbs_variables.i_blkt_inboard == 1:
            return [
                dpres_fw_inboard_coolant,
                dpres_fw_outboard_coolant,
                dpres_blkt_outboard_coolant,
                dpres_blkt_inboard_coolant,
            ]
        return [
            dpres_fw_inboard_coolant,
            dpres_fw_outboard_coolant,
            dpres_blkt_outboard_coolant,
        ]

    def blanket_module_poloidal_height(self):
        """Calculations for blanket module poloidal height
        author: J. Morris, CCFE, Culham Science Centre
        Calculations for blanket module poloidal height for D shaped and elliptical machines
        """
        if (
            physics_variables.itart == 1 or fwbs_variables.i_fw_blkt_vv_shape == 1
        ):  # D-shaped machine
            # Segment vertical inboard surface (m)
            blanket_library.len_blkt_inboard_segment_poloidal = (
                2.0 * blanket_library.dz_blkt_half
            ) / fwbs_variables.n_blkt_inboard_modules_poloidal

            # Calculate perimeter of ellipse that defines the internal
            # surface of the outboard first wall / blanket

            # Mid-plane distance from inboard to outboard side (m)
            a = (
                build_variables.dr_fw_plasma_gap_inboard
                + 2.0 * physics_variables.rminor
                + build_variables.dr_fw_plasma_gap_outboard
            )

            # Internal half-height of blanket (m)
            b = blanket_library.dz_blkt_half

            # Calculate ellipse circumference using Ramanujan approximation (m)
            ptor = np.pi * (3.0 * (a + b) - np.sqrt((3.0 * a + b) * (a + 3.0 * b)))

            # Calculate blanket poloidal length and segment, subtracting divertor length (m)
            # kit hcll version only had the single null option
            if divertor_variables.n_divertors == 2:
                # Double null configuration
                blanket_library.len_blkt_outboard_segment_poloidal = (
                    0.5
                    * ptor
                    * (1.0 - 2.0 * fwbs_variables.f_ster_div_single)
                    / fwbs_variables.n_blkt_outboard_modules_poloidal
                )
            else:
                # single null configuration
                blanket_library.len_blkt_outboard_segment_poloidal = (
                    0.5
                    * ptor
                    * (1.0 - fwbs_variables.f_ster_div_single)
                    / fwbs_variables.n_blkt_outboard_modules_poloidal
                )

        # shape defined by two half-ellipses
        else:
            # Major radius where half-ellipses 'meet' (m)
            r1 = (
                physics_variables.rmajor
                - physics_variables.rminor * physics_variables.triang
            )

            # Internal half-height of blanket (m)
            b = blanket_library.dz_blkt_half

            # Distance between r1 and nearest edge of inboard first wall / blanket (m)
            a = r1 - (
                physics_variables.rmajor
                - physics_variables.rminor
                - build_variables.dr_fw_plasma_gap_inboard
            )

            # Calculate ellipse circumference using Ramanujan approximation (m)
            ptor = np.pi * (3.0 * (a + b) - np.sqrt((3.0 * a + b) * (a + 3.0 * b)))

            # Calculate inboard blanket poloidal length and segment, subtracting divertor length (m)
            # Assume divertor lies between the two ellipses, so fraction f_ster_div_single still applies

            # kit hcll version only had the single null option
            if divertor_variables.n_divertors == 2:
                # Double null configuration
                blanket_library.len_blkt_inboard_segment_poloidal = (
                    0.5
                    * ptor
                    * (1.0 - 2.0 * fwbs_variables.f_ster_div_single)
                    / fwbs_variables.n_blkt_inboard_modules_poloidal
                )
            else:
                # single null configuration
                blanket_library.len_blkt_inboard_segment_poloidal = (
                    0.5
                    * ptor
                    * (1.0 - fwbs_variables.f_ster_div_single)
                    / fwbs_variables.n_blkt_inboard_modules_poloidal
                )

            # Distance between r1 and inner edge of outboard first wall / blanket (m)
            a = (
                physics_variables.rmajor
                + physics_variables.rminor
                + build_variables.dr_fw_plasma_gap_outboard
                - r1
            )

            # Calculate ellipse circumference using Ramanujan approximation (m)
            ptor = np.pi * (3.0 * (a + b) - np.sqrt((3.0 * a + b) * (a + 3.0 * b)))

            # kit hcll version only had the single null option
            # Calculate outboard blanket poloidal length and segment, subtracting divertor length (m)
            if divertor_variables.n_divertors == 2:
                # Double null configuration
                blanket_library.len_blkt_outboard_segment_poloidal = (
                    0.5
                    * ptor
                    * (1.0 - 2.0 * fwbs_variables.f_ster_div_single)
                    / fwbs_variables.n_blkt_outboard_modules_poloidal
                )
            else:
                # single null configuration
                blanket_library.len_blkt_outboard_segment_poloidal = (
                    0.5
                    * ptor
                    * (1.0 - fwbs_variables.f_ster_div_single)
                    / fwbs_variables.n_blkt_outboard_modules_poloidal
                )

    def liquid_breeder_properties(self, output: bool = False):
        """Calculates the fluid properties of the Liquid Metal Breeder/Coolant in the Blanket BZ
        Uses middle value of input and output temperatures of Liquid Metal Breeder/Coolant
        Curently have PbLi but can expand with e.g., Lithium

        author: G Graham, CCFE

        References:

             [Mal1995]   Malang and Mattas (1995), Comparison of lithium and the eutectic
                         lead-lithium alloy, two candidate liquid metal breeder materials
                         for self-cooled blankets, Fusion Engineering and Design 27, 399-406.

             [Mas2008]   Mas de les Valles et al. (2008), Lead-lithium material database for
                         nuclear fusion technology, Journal of Nuclear Materials, Vol. 376(6).

             [Mar2019]   Martelli et al. (2019), Literature review of lead-lithium
                         thermophysical properties, Fusion Engineering and Design, 138, 183-195.
        """

        # Use mid temp
        if fwbs_variables.inlet_temp_liq == fwbs_variables.outlet_temp_liq:
            mid_temp_liq = fwbs_variables.outlet_temp_liq
        else:
            mid_temp_liq = (
                fwbs_variables.inlet_temp_liq + fwbs_variables.outlet_temp_liq
            ) * 0.5

        # If the liquid metal is PbLi...
        if fwbs_variables.i_blkt_liquid_breeder_type == 0:
            # PbLi from [Mar2019]
            # Constant pressure ~ 17 atmospheres ~ 1.7D6 Pa
            # Li content is ~ 17%
            #
            # density                      kg m-3          T in Kelvin     range = 508-880 K
            #
            # specific_heat                J kg-1 K-1      T in Kelvin     range = 508-880 K
            #
            # thermal_conductivity         W m-1 K-1       T in Celcius    range = 508-773 K
            #
            # dynamic_viscosity            Pa s            T in Celcius    range = 508-873 K
            #
            # electrical_conductivity      A V-1 m-1       T in Kelvin     range = 600-800 K

            # Caculate properties
            fwbs_variables.den_liq = 1.052e4 * (1 - mid_temp_liq * 1.13e-4)

            fwbs_variables.specific_heat_liq = 1.95e2 - mid_temp_liq * 9.116e-3

            fwbs_variables.thermal_conductivity_liq = (
                1.95 + (mid_temp_liq - 273.15) * 1.96e-2
            )

            fwbs_variables.dynamic_viscosity_liq = (
                6.11e-3
                - (2.257e-5 * (mid_temp_liq - 273.15))
                + (3.766e-8 * (mid_temp_liq - 273.15) ** 2)
                - (2.289e-11 * (mid_temp_liq - 273.15) ** 3)
            )

            t_ranges = np.zeros((5, 2))

            t_ranges[:4, 0] = 508.0
            t_ranges[:4, 1] = 880.0

            fwbs_variables.electrical_conductivity_liq = 1.0 / (
                1.03e-6 - (6.75e-11 * mid_temp_liq) + (4.18e-13 * mid_temp_liq**2)
            )

            t_ranges[4, 0] = 600.0
            t_ranges[4, 1] = 800.0

        # If the liquid metal is Li...
        elif fwbs_variables.i_blkt_liquid_breeder_type == 1:
            # Temporary - should be updated with information from Li reviews conducted at CCFE once completed
            # Li Properties from [Mal1995] at 300 Celcius
            # den_liq = 505                            kg/m3
            # specific_heat_liq = 4260                 J kg-1 K-1
            # thermal_conductivity_liq = 46            W m-1 K-1
            # dynamic_viscosity_liq = 1.0D-6           m2 s-1
            # electrical_conductivity_liq = 3.03D6     A V-1 m-1

            # New from 'Application of lithium in systems of fusion reactors. 1. Physical and chemical properties of lithium'
            # Lyublinski et al., 2009, Plasma Devicec and Operations
            fwbs_variables.den_liq = (
                504.43
                - (0.2729 * mid_temp_liq)
                - (8.0035e-5 * mid_temp_liq**2)
                + (3.799e-8 * mid_temp_liq**3)
            )
            fwbs_variables.specific_heat_liq = (
                31.227
                + (0.205e6 * mid_temp_liq ** (-2))
                - (5.265e-3 * mid_temp_liq)
                + (2.628e6 * mid_temp_liq ** (-2))
            )
            # thermal_conductivity_liq also in paper
            fwbs_variables.dynamic_viscosity_liq = np.exp(
                -4.16e0 - (0.64 * np.log(mid_temp_liq)) + (262.1 / mid_temp_liq)
            )
            fwbs_variables.electrical_conductivity_liq = (
                (0.9249e9 * mid_temp_liq) + 2.3167e6 - (0.7131e3 * mid_temp_liq)
            )

        # Magnetic feild strength in T for Hartmann calculation
        # IB
        if fwbs_variables.i_blkt_inboard == 1:
            fwbs_variables.b_mag_blkt[0] = (
                physics_variables.b_plasma_toroidal_on_axis
                * physics_variables.rmajor
                / (
                    physics_variables.rmajor
                    - (physics_variables.rmajor / physics_variables.aspect)
                    - (build_variables.dr_blkt_inboard / 2)
                )
            )
        # We do not use this if there is no IB blanket, but will use edge as fill value
        if fwbs_variables.i_blkt_inboard == 0:
            fwbs_variables.b_mag_blkt[0] = (
                physics_variables.b_plasma_toroidal_on_axis
                * physics_variables.rmajor
                / (
                    physics_variables.rmajor
                    - (physics_variables.rmajor / physics_variables.aspect)
                )
            )
        # OB
        fwbs_variables.b_mag_blkt[1] = (
            physics_variables.b_plasma_toroidal_on_axis
            * physics_variables.rmajor
            / (
                physics_variables.rmajor
                + (physics_variables.rmajor / physics_variables.aspect)
                + (build_variables.dr_blkt_outboard / 2)
            )
        )

        # Calculate Hartmann number
        con_vsc_rat = (
            fwbs_variables.electrical_conductivity_liq
            / fwbs_variables.dynamic_viscosity_liq
        )
        # Use toroidal width of the rectangular cooling channel as characteristic length scale
        fwbs_variables.hartmann_liq = (
            np.asarray(fwbs_variables.b_mag_blkt)
            * fwbs_variables.a_bz_liq
            / 2.0
            * np.sqrt(con_vsc_rat)
        )

        # Error for temperature range of breeder property realtions
        if fwbs_variables.i_blkt_liquid_breeder_type == 0 and (
            (t_ranges[:, 0] > mid_temp_liq).any()
            or (t_ranges[:, 1] < mid_temp_liq).any()
        ):
            logger.error(
                "Outside temperature limit for one or more liquid metal breeder properties"
            )

            if output:
                po.ocmmnt(
                    self.outfile,
                    "Outside temperature limit for one or more liquid metal breeder properties.",
                )
                po.ovarrf(
                    self.outfile,
                    "Liquid metal temperature (K)",
                    "(mid_temp_liq)",
                    mid_temp_liq,
                    "OP ",
                )
                po.ocmmnt(self.outfile, "Density: Max T = 880 K, Min T = 508 K")
                po.ocmmnt(self.outfile, "Specific heat: Max T = 880 K, Min T = 508 K")
                po.ocmmnt(
                    self.outfile, "Thermal conductivity: Max T = 880 K, Min T = 508 K"
                )
                po.ocmmnt(
                    self.outfile, "Dynamic viscosity : Max T = 880 K, Min T = 508 K"
                )
                po.ocmmnt(
                    self.outfile,
                    "Electrical conductivity: Max T = 800 K, Min T = 600 K",
                )

        if not output:
            return

        po.oheadr(self.outfile, "Blanket : Liquid Breeder Properties")

        if fwbs_variables.i_blkt_dual_coolant == 1:
            po.ocmmnt(
                self.outfile,
                "Single coolant: liquid metal circulted for tritium extraction.",
            )
        if fwbs_variables.i_blkt_dual_coolant == 2:
            po.ocmmnt(self.outfile, "Dual coolant: self-cooled liquid metal breeder.")

        if fwbs_variables.i_blkt_liquid_breeder_type == 0:
            po.ocmmnt(
                self.outfile,
                "Blanket breeder type (i_blkt_liquid_breeder_type=0), PbLi (~ 17% Li)",
            )
        if fwbs_variables.i_blkt_liquid_breeder_type == 1:
            po.ocmmnt(
                self.outfile, "Blanket breeder type (i_blkt_liquid_breeder_type=1), Li"
            )

        po.ovarrf(
            self.outfile, "Density (kg m-3)", "(den_liq)", fwbs_variables.den_liq, "OP "
        )
        po.ovarrf(
            self.outfile,
            "Viscosity (Pa s)",
            "(dynamic_viscosity_liq)",
            fwbs_variables.dynamic_viscosity_liq,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Electrical Conductivity (A V-1 m-1)",
            "(electrical_conductivity_liq)",
            fwbs_variables.electrical_conductivity_liq,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Hartmann Number IB",
            "(hartmann_liq)",
            fwbs_variables.hartmann_liq[0],
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Hartmann Number OB",
            "(hartmann_liq)",
            fwbs_variables.hartmann_liq[0],
            "OP ",
        )

        po.ovarre(
            self.outfile,
            "Inlet Temperature (Celcius)",
            "(inlet_temp_liq)",
            fwbs_variables.inlet_temp_liq,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Outlet Temperature (Celcius)",
            "(outlet_temp_liq)",
            fwbs_variables.outlet_temp_liq,
            "OP ",
        )

    def flow_velocity(self, i_channel_shape, mass_flow_rate, flow_density):
        """Calculate the coolant flow velocity (m/s) for given pipe mass flow rate and pipe size/shape.
        N.B. Assumed that primary BB and FW coolants have same pipe radius (= radius_fw_channel).
        author: G. Graham, CCFE

        :param i_channel_shape: Switch for circular or rectangular channel crossection.
            Shape depends on whether primary or secondary coolant.
            1: circle (primary)
            2: rectangle (secondary)
        :param mass_flow_rate: Coolant mass flow rate per pipe (kg/s)
        :param flow_density: Coolant density
        """

        if i_channel_shape == 1:
            return mass_flow_rate / (
                flow_density
                * np.pi
                * fwbs_variables.radius_fw_channel
                * fwbs_variables.radius_fw_channel
            )

        # If secondary coolant then rectangular channels assumed
        if i_channel_shape == 2:
            return mass_flow_rate / (
                flow_density * fwbs_variables.a_bz_liq * fwbs_variables.b_bz_liq
            )

        raise ProcessValueError(
            f"i_channel_shape ={i_channel_shape} is an invalid option."
        )

    def thermo_hydraulic_model(self, output: bool):
        """
        Thermo-hydraulic model for first wall and blanket
        ONLY CALLED if i_p_coolant_pumping = 2 or 3

        Calculations for detailed powerflow model i_thermal_electric_conversion > 1

        original author: J. Morris, CCFE, Culham Science Centre
        Dual-coolant modifications and generalisation refactor: G. Graham, CCFE

        Three options:
        1.   Solid breeder - nuclear heating in the blanket is exctrated by the primary coolant.
        2.   Liquid metal breeder, single-coolant
                 - nuclear heating in the blanket is exctrated by the primary coolant.
                 - liquid metal is circulated for tritium extraction, specified by number of circulations/day.
        3.   Liquid metal breeder, dual-coolant -
                 - nuclear heating in the liquid breeder/coolant is extracted by the liquid breeder/coolant.
                 - nuclear heating in the blanket structure is extracted by the primary coolant

        Flow Channel and Coolant Input Info:

            N.B. Primary coolant applies to single-coolant BB, or structural cooling of dual-coolant BB.
            Secondary coolant applies to self-cooled breeder material.

            Coolant Channels            FW                      BB primary          BB Liquid Breeder/Coolant

            length (m)                  len_fw_channel
            width (m)                   radius_fw_channel (radius, cicular)   radius_fw_channel                 a_bz_liq, b_bz_liq (rectangular)
            wall thickness (m)          dr_fw_wall                 dr_fw_wall             th_wall_secondary
            dx_fw_module (m)                   dx_fw_module
            roughness epsilon           roughness_fw_channel
            peak FW temp (K)            temp_fw_peak
            maximum temp (K)            temp_fw_max
            FCI switch                  ---                     ---                 i_blkt_liquid_breeder_channel_type

            Coolant                     FW                      BB primary          BB secondary

            primary coolant switch      i_fw_coolant_type               i_blkt_coolant_type              ---
            secondary coolant switch    ---                     ---                 i_blkt_liquid_breeder_type
            inlet temp (K)              temp_fw_coolant_in                 temp_blkt_coolant_in          inlet_temp_liq
            outlet temp (K)             temp_fw_coolant_out                temp_blkt_coolant_out         outlet_temp_liq
            pressure (Pa)               pres_fw_coolant              pres_blkt_coolant          blpressure_liq
        """
        ######################################################
        # Pre calculations needed for thermo-hydraulic model #
        ######################################################
        # IB/OB FW (MW)
        blanket_library.p_fw_inboard_nuclear_heat_mw = (
            fwbs_variables.p_fw_nuclear_heat_total_mw
            * build_variables.a_fw_inboard
            / build_variables.a_fw_total
        )
        blanket_library.p_fw_outboard_nuclear_heat_mw = (
            fwbs_variables.p_fw_nuclear_heat_total_mw
            * build_variables.a_fw_outboard
            / build_variables.a_fw_total
        )

        # IB/OB Blanket (MW)

        # Neutron power deposited in inboard blanket (MW)
        if fwbs_variables.i_blkt_inboard == 1:
            blanket_library.p_blkt_nuclear_heat_inboard_mw = (
                fwbs_variables.p_blkt_nuclear_heat_total_mw
                * fwbs_variables.vol_blkt_inboard
                / fwbs_variables.vol_blkt_total
            )

        # Neutron power deposited in outboard blanket (MW)
        blanket_library.p_blkt_nuclear_heat_outboard_mw = (
            fwbs_variables.p_blkt_nuclear_heat_total_mw
            * fwbs_variables.vol_blkt_outboard
            / fwbs_variables.vol_blkt_total
        )

        # For a dual-coolant blanket, some fraction of the power goes into the
        # structure of the BZ and is cooled by the primary coolant, and some fraction
        # goes into the liquid breeder to be cooled by itself.

        # If the blanket is dual-coolant...
        if fwbs_variables.i_blkt_dual_coolant == 2:
            f_nuc_pow_bz_liq = 1 - fwbs_variables.f_nuc_pow_bz_struct

            # Inboard blanket calc. Will return 0 if no inboard dr_shld_inboard thickness
            pnucblkti_struct = (
                fwbs_variables.p_blkt_nuclear_heat_total_mw
                * fwbs_variables.f_nuc_pow_bz_struct
            ) * (fwbs_variables.vol_blkt_inboard / fwbs_variables.vol_blkt_total)
            pnucblkti_liq = (
                fwbs_variables.p_blkt_nuclear_heat_total_mw * f_nuc_pow_bz_liq
            ) * (fwbs_variables.vol_blkt_inboard / fwbs_variables.vol_blkt_total)
            pnucblkto_struct = (
                fwbs_variables.p_blkt_nuclear_heat_total_mw
                * fwbs_variables.f_nuc_pow_bz_struct
            ) * (fwbs_variables.vol_blkt_outboard / fwbs_variables.vol_blkt_total)
            pnucblkto_liq = (
                fwbs_variables.p_blkt_nuclear_heat_total_mw * f_nuc_pow_bz_liq
            ) * (fwbs_variables.vol_blkt_outboard / fwbs_variables.vol_blkt_total)

        # FW and BB Mass Flow ###########

        # Make sure that, if the inputs for the FW and blanket inputs are different,
        # the i_fw_blkt_shared_coolant variable is appropriately set for seperate coolants
        if (
            fwbs_variables.i_fw_coolant_type == "Helium"
            and fwbs_variables.i_blkt_coolant_type == 2
        ):
            fwbs_variables.i_fw_blkt_shared_coolant = 1
        if (
            fwbs_variables.i_fw_coolant_type == "Water"
            and fwbs_variables.i_blkt_coolant_type == 1
        ):
            fwbs_variables.i_fw_blkt_shared_coolant = 1

        # If FW and BB have the same coolant...
        if fwbs_variables.i_fw_blkt_shared_coolant == 0:
            # Fraction of heat to be removed by IB/OB FW
            if fwbs_variables.i_blkt_dual_coolant == 2:
                f_nuc_fwi = (
                    blanket_library.p_fw_inboard_nuclear_heat_mw
                    + fwbs_variables.psurffwi
                ) / (
                    blanket_library.p_fw_inboard_nuclear_heat_mw
                    + fwbs_variables.psurffwi
                    + pnucblkti_struct
                )
                f_nuc_fwo = (
                    blanket_library.p_fw_outboard_nuclear_heat_mw
                    + fwbs_variables.psurffwo
                ) / (
                    blanket_library.p_fw_outboard_nuclear_heat_mw
                    + fwbs_variables.psurffwo
                    + pnucblkto_struct
                )
            else:
                f_nuc_fwi = (
                    blanket_library.p_fw_inboard_nuclear_heat_mw
                    + fwbs_variables.psurffwi
                ) / (
                    blanket_library.p_fw_inboard_nuclear_heat_mw
                    + fwbs_variables.psurffwi
                    + blanket_library.p_blkt_nuclear_heat_inboard_mw
                )
                f_nuc_fwo = (
                    blanket_library.p_fw_outboard_nuclear_heat_mw
                    + fwbs_variables.psurffwo
                ) / (
                    blanket_library.p_fw_outboard_nuclear_heat_mw
                    + fwbs_variables.psurffwo
                    + blanket_library.p_blkt_nuclear_heat_outboard_mw
                )

            # Outlet FW/inlet BB temp (mass flow FW = mass flow BB)
            if fwbs_variables.i_blkt_inboard == 1:
                fwoutleti = (f_nuc_fwi * fwbs_variables.temp_blkt_coolant_out) + (
                    1 - f_nuc_fwi
                ) * fwbs_variables.temp_fw_coolant_in
                inlet_tempi = fwoutleti
            else:
                fwoutleti = fwbs_variables.temp_fw_coolant_out

            fwoutleto = (f_nuc_fwo * fwbs_variables.temp_blkt_coolant_out) + (
                1 - f_nuc_fwo
            ) * fwbs_variables.temp_fw_coolant_in
            inlet_tempo = fwoutleto

        elif fwbs_variables.i_fw_blkt_shared_coolant == 1:
            fwoutleti = fwbs_variables.temp_fw_coolant_out
            inlet_tempi = fwbs_variables.temp_blkt_coolant_in
            fwoutleto = fwbs_variables.temp_fw_coolant_out
            inlet_tempo = fwbs_variables.temp_blkt_coolant_in

        # Maximum FW temperature. (27/11/2015) Issue #348
        # First wall flow is just along the first wall, with no allowance for radial
        # pipes, manifolds etc. The outputs are mid quantities of inlet and outlet.
        # This subroutine recalculates cp and rhof.
        (
            blanket_library.temp_fw_inboard_peak,
            _,
            _,
            blanket_library.mflow_fw_inboard_coolant_channel,
        ) = self.fw.fw_temp(
            output,
            fwbs_variables.radius_fw_channel,
            build_variables.dr_fw_inboard,
            build_variables.a_fw_inboard,
            fwbs_variables.psurffwi,
            blanket_library.p_fw_inboard_nuclear_heat_mw,
            "Inboard first wall",
        )
        (
            blanket_library.temp_fw_outboard_peak,
            _cf,
            _rhof,
            blanket_library.mflow_fw_outboard_coolant_channel,
        ) = self.fw.fw_temp(
            output,
            fwbs_variables.radius_fw_channel,
            build_variables.dr_fw_outboard,
            build_variables.a_fw_outboard,
            fwbs_variables.psurffwo,
            blanket_library.p_fw_outboard_nuclear_heat_mw,
            "Outboard first wall",
        )

        # Peak first wall temperature (K)
        fwbs_variables.temp_fw_peak = max(
            blanket_library.temp_fw_inboard_peak, blanket_library.temp_fw_outboard_peak
        )

        # Total mass flow rate to remove inboard FW power (kg/s)
        blanket_library.mflow_fw_inboard_coolant_total = (
            1.0e6
            * (blanket_library.p_fw_inboard_nuclear_heat_mw + fwbs_variables.psurffwi)
            / (fwbs_variables.cp_fw * (fwoutleti - fwbs_variables.temp_fw_coolant_in))
        )
        # Total mass flow rate to remove outboard FW power (kg/s)
        blanket_library.mflow_fw_outboard_coolant_total = (
            1.0e6
            * (blanket_library.p_fw_outboard_nuclear_heat_mw + fwbs_variables.psurffwo)
            / (fwbs_variables.cp_fw * (fwoutleto - fwbs_variables.temp_fw_coolant_in))
        )

        # If the blanket is dual-coolant...
        if fwbs_variables.i_blkt_dual_coolant == 2:
            # Mass flow rates for outboard blanket coolants (kg/s)
            blanket_library.mflow_blkt_outboard_coolant = (
                1.0e6
                * (pnucblkto_struct)
                / (
                    fwbs_variables.cp_bl
                    * (fwbs_variables.temp_blkt_coolant_out - inlet_tempo)
                )
            )
            blanket_library.mfblkto_liq = (
                1.0e6
                * (pnucblkto_liq)
                / (
                    fwbs_variables.specific_heat_liq
                    * (fwbs_variables.outlet_temp_liq - fwbs_variables.inlet_temp_liq)
                )
            )

            # If there is an IB blanket...
            if fwbs_variables.i_blkt_inboard == 1:
                # Mass flow rates for inboard blanket coolants (kg/s)
                blanket_library.mflow_blkt_inboard_coolant = (
                    1.0e6
                    * (pnucblkti_struct)
                    / (
                        fwbs_variables.cp_bl
                        * (fwbs_variables.temp_blkt_coolant_out - inlet_tempi)
                    )
                )
                blanket_library.mfblkti_liq = (
                    1.0e6
                    * (pnucblkti_liq)
                    / (
                        fwbs_variables.specific_heat_liq
                        * (
                            fwbs_variables.outlet_temp_liq
                            - fwbs_variables.inlet_temp_liq
                        )
                    )
                )

        # If the blanket is single-coolant with liquid metal breeder...
        elif fwbs_variables.i_blkt_dual_coolant == 1:
            # Mass flow rate for outboard blanket coolant (kg/s)
            blanket_library.mflow_blkt_outboard_coolant = (
                1.0e6
                * (blanket_library.p_blkt_nuclear_heat_outboard_mw)
                / (
                    fwbs_variables.cp_bl
                    * (fwbs_variables.temp_blkt_coolant_out - inlet_tempo)
                )
            )

            # Get mass flow rate etc. for inboard blanket breeder flow for tritium extraction
            # Use the number of desired recirculations ([Aub2013]=10) and mass from dcll_masses
            # N.B. wht_liq is BZ mass, does not include manifold.
            blanket_library.mfblkto_liq = (
                fwbs_variables.n_liq_recirc * fwbs_variables.wht_liq_ob
            ) / (24 * 3600)

            # If there is an IB blanket...
            if fwbs_variables.i_blkt_inboard == 1:
                # Mass flow rate for inboard blanket coolant (kg/s)
                blanket_library.mflow_blkt_inboard_coolant = (
                    1.0e6
                    * (blanket_library.p_blkt_nuclear_heat_inboard_mw)
                    / (
                        fwbs_variables.cp_bl
                        * (fwbs_variables.temp_blkt_coolant_out - inlet_tempi)
                    )
                )
                # Mass flow rate for inboard breeder flow (kg/s)
                fwbs_variables.mfblkti_liq = (
                    fwbs_variables.n_liq_recirc * fwbs_variables.wht_liq_ib
                ) / (24 * 3600)

        # If the blanket is single-coolant with solid breeder...
        else:
            # Mass flow rate for inboard blanket coolant (kg/s)
            blanket_library.mflow_blkt_outboard_coolant = (
                1.0e6
                * (blanket_library.p_blkt_nuclear_heat_outboard_mw)
                / (
                    fwbs_variables.cp_bl
                    * (fwbs_variables.temp_blkt_coolant_out - inlet_tempo)
                )
            )

            # If there is an IB blanket...
            # Mass flow rate for inboard blanket coolant (kg/s)
            if fwbs_variables.i_blkt_inboard == 1:
                blanket_library.mflow_blkt_inboard_coolant = (
                    1.0e6
                    * (blanket_library.p_blkt_nuclear_heat_inboard_mw)
                    / (
                        fwbs_variables.cp_bl
                        * (fwbs_variables.temp_blkt_coolant_out - inlet_tempi)
                    )
                )

        ########################################################
        # Handling of pressure drops and coolant pumping power #
        ########################################################

        # load in pressures if primary pumping == 2
        if fwbs_variables.i_p_coolant_pumping == 2:
            deltap = self.thermo_hydraulic_model_pressure_drop_calculations(
                output=output
            )
            deltap_fwi = deltap[0]
            deltap_fwo = deltap[1]
            deltap_blo = deltap[2]
            if fwbs_variables.i_blkt_dual_coolant > 0:
                if fwbs_variables.i_blkt_inboard == 1:
                    deltap_bli = deltap[3]
                    deltap_blo_liq = deltap[4]
                    deltap_bli_liq = deltap[5]
                else:
                    deltap_blo_liq = deltap[3]
            else:
                if fwbs_variables.i_blkt_inboard == 1:
                    deltap_bli = deltap[3]

        # Pumping Power
        # If FW and BB have the same coolant...
        if fwbs_variables.i_fw_blkt_shared_coolant == 0:
            # Total pressure drop in the first wall/blanket  (Pa)
            if fwbs_variables.i_p_coolant_pumping == 2:
                if fwbs_variables.i_blkt_inboard == 1:
                    deltap_fw_blkt = deltap_fwi + deltap_bli + deltap_fwo + deltap_blo
                if fwbs_variables.i_blkt_inboard == 0:
                    deltap_fw_blkt = deltap_fwi + deltap_fwo + deltap_blo
            elif fwbs_variables.i_p_coolant_pumping == 3:
                deltap_fw_blkt = primary_pumping_variables.dp_fw_blkt
            # Total coolant mass flow rate in the first wall/blanket (kg/s)
            blanket_library.mftotal = (
                blanket_library.mflow_fw_inboard_coolant_total
                + blanket_library.mflow_fw_outboard_coolant_total
            )

            # Total mechanical pumping power (MW)
            primary_pumping_variables.p_fw_blkt_coolant_pump_mw = (
                self.coolant_pumping_power(
                    output=output,
                    i_liquid_breeder=1,
                    temp_coolant_pump_outlet=fwbs_variables.temp_fw_coolant_in,
                    temp_coolant_pump_inlet=fwbs_variables.temp_blkt_coolant_out,
                    pres_coolant_pump_inlet=fwbs_variables.pres_fw_coolant,
                    dpres_coolant=deltap_fw_blkt,
                    mflow_coolant_total=blanket_library.mftotal,
                    primary_coolant_switch=fwbs_variables.i_fw_coolant_type,
                    den_coolant=fwbs_variables.den_fw_coolant,
                    label="First Wall and Blanket",
                )
            )

        # If FW and BB have different coolants...
        elif fwbs_variables.i_fw_blkt_shared_coolant == 1:
            if fwbs_variables.i_p_coolant_pumping == 2:
                # Total pressure drop in the first wall (Pa)
                deltap_fw = deltap_fwi + deltap_fwo

                # Total pressure drop in the blanket (Pa)
                if fwbs_variables.i_blkt_inboard == 1:
                    deltap_blkt = deltap_bli + deltap_blo
                if fwbs_variables.i_blkt_inboard == 0:
                    deltap_blkt = deltap_blo
            elif fwbs_variables.i_p_coolant_pumping == 3:
                deltap_fw = primary_pumping_variables.dp_fw
                deltap_blkt = primary_pumping_variables.dp_blkt

            # Total coolant mass flow rate in the first wall (kg/s)
            blanket_library.mflow_fw_coolant_total = (
                blanket_library.mflow_fw_inboard_coolant_total
                + blanket_library.mflow_fw_outboard_coolant_total
            )
            # Total coolant mass flow rate in the blanket (kg/s)
            blanket_library.mflow_blkt_coolant_total = (
                blanket_library.mflow_blkt_inboard_coolant
                + blanket_library.mflow_blkt_outboard_coolant
            )

            # Mechanical pumping power for the first wall (MW)
            heat_transport_variables.p_fw_coolant_pump_mw = self.coolant_pumping_power(
                output=output,
                i_liquid_breeder=1,
                temp_coolant_pump_outlet=fwbs_variables.temp_fw_coolant_in,
                temp_coolant_pump_inlet=fwbs_variables.temp_fw_coolant_out,
                pres_coolant_pump_inlet=fwbs_variables.pres_fw_coolant,
                dpres_coolant=deltap_fw,
                mflow_coolant_total=blanket_library.mflow_fw_coolant_total,
                primary_coolant_switch=fwbs_variables.i_fw_coolant_type,
                den_coolant=fwbs_variables.den_fw_coolant,
                label="First Wall",
            )

            # Mechanical pumping power for the blanket (MW)
            heat_transport_variables.p_blkt_coolant_pump_mw = self.coolant_pumping_power(
                output=output,
                i_liquid_breeder=1,
                temp_coolant_pump_outlet=fwbs_variables.temp_blkt_coolant_in,
                temp_coolant_pump_inlet=fwbs_variables.temp_blkt_coolant_out,
                pres_coolant_pump_inlet=fwbs_variables.pres_blkt_coolant,
                dpres_coolant=deltap_blkt,
                mflow_coolant_total=blanket_library.mflow_blkt_coolant_total,
                primary_coolant_switch=(
                    "Helium" if fwbs_variables.i_blkt_coolant_type == 1 else "Water"
                ),
                den_coolant=fwbs_variables.den_blkt_coolant,
                label="Blanket",
            )

            # Total mechanical pumping power (MW)
            primary_pumping_variables.p_fw_blkt_coolant_pump_mw = (
                heat_transport_variables.p_fw_coolant_pump_mw
                + heat_transport_variables.p_blkt_coolant_pump_mw
            )

        # If the blanket has a liquid metal breeder...
        if fwbs_variables.i_blkt_dual_coolant > 0:
            # Total pressure drop in the blanket (Pa)
            if fwbs_variables.i_p_coolant_pumping == 2:
                if fwbs_variables.i_blkt_inboard == 1:
                    deltap_bl_liq = deltap_bli_liq + deltap_blo_liq
                if fwbs_variables.i_blkt_inboard == 0:
                    deltap_bl_liq = deltap_blo_liq
            elif fwbs_variables.i_p_coolant_pumping == 3:
                deltap_bl_liq = primary_pumping_variables.dp_liq
            # Total liquid metal breeder/coolant mass flow rate in the blanket (kg/s)
            blanket_library.mfblkt_liq = (
                blanket_library.mfblkti_liq + blanket_library.mfblkto_liq
            )

            # Mechanical pumping power for the blanket (MW)
            heat_transport_variables.p_blkt_breeder_pump_mw = self.coolant_pumping_power(
                output=output,
                i_liquid_breeder=2,
                temp_coolant_pump_outlet=fwbs_variables.inlet_temp_liq,
                temp_coolant_pump_inlet=fwbs_variables.outlet_temp_liq,
                pres_coolant_pump_inlet=fwbs_variables.blpressure_liq,
                dpres_coolant=deltap_bl_liq,
                mflow_coolant_total=blanket_library.mfblkt_liq,
                primary_coolant_switch=(
                    "Helium" if fwbs_variables.i_blkt_coolant_type == 1 else "Water"
                ),
                den_coolant=fwbs_variables.den_liq,
                label="Liquid Metal Breeder/Coolant",
            )

            heat_transport_variables.htpmw_blkt_tot = (
                primary_pumping_variables.p_fw_blkt_coolant_pump_mw
                + heat_transport_variables.p_blkt_breeder_pump_mw
            )

        if output:
            po.oheadr(self.outfile, "Summary of first wall and blanket thermohydraulics")

            # FW
            po.osubhd(self.outfile, "First wall: ")

            po.ovarst(
                self.outfile,
                "First wall coolant type",
                "(i_fw_coolant_type)",
                fwbs_variables.i_fw_coolant_type,
            )
            po.ovarre(
                self.outfile,
                "Wall thickness of first wall cooling channels (m)",
                "(dr_fw_wall)",
                fwbs_variables.dr_fw_wall,
            )
            po.ovarre(
                self.outfile,
                "Radius of first wall cooling channels (m)",
                "(radius_fw_channel)",
                fwbs_variables.radius_fw_channel,
            )
            po.ovarre(
                self.outfile,
                "Radius of blanket cooling channels (m)",
                "(radius_blkt_channel)",
                fwbs_variables.radius_blkt_channel,
            )
            po.ovarre(
                self.outfile,
                "Roughness of first wall cooling channels (m)",
                "(roughness_fw_channel)",
                fwbs_variables.roughness_fw_channel,
            )
            po.ovarrf(
                self.outfile,
                "Inlet temperature of first wall coolant (K)",
                "(temp_fw_coolant_in)",
                fwbs_variables.temp_fw_coolant_in,
            )
            po.ovarrf(
                self.outfile,
                "Outlet temperature of first wall coolant (K)",
                "(temp_fw_coolant_out)",
                fwbs_variables.temp_fw_coolant_out,
            )
            po.ovarre(
                self.outfile,
                "First wall coolant pressure (Pa)",
                "(pres_fw_coolant)",
                fwbs_variables.pres_fw_coolant,
            )
            if fwbs_variables.i_fw_blkt_shared_coolant == 1:
                po.ovarre(
                    self.outfile,
                    "First wall coolant mass flow rate (kg/s)",
                    "(mflow_fw_coolant_total)",
                    blanket_library.mflow_fw_coolant_total,
                    "OP ",
                )
            po.ovarrf(
                self.outfile,
                "Allowable temperature of first wall material, excluding armour (K)",
                "(temp_fw_max)",
                fwbs_variables.temp_fw_max,
            )
            po.ovarrf(
                self.outfile,
                "Actual peak temperature of first wall material (K)",
                "(temp_fw_peak)",
                fwbs_variables.temp_fw_peak,
                "OP ",
            )

            # BB
            po.osubhd(self.outfile, "Breeding Blanket (primary): ")
            po.ovarre(
                self.outfile,
                "Blanket half height (m)",
                "(dz_blkt_half)",
                blanket_library.dz_blkt_half,
            )
            po.ovarin(
                self.outfile,
                "Blanket coolant type (1=He, 2=H20)",
                "(i_blkt_coolant_type)",
                fwbs_variables.i_blkt_coolant_type,
            )
            po.ovarrf(
                self.outfile,
                "Inlet temperature of blanket coolant (K)",
                "(temp_blkt_coolant_in)",
                fwbs_variables.temp_blkt_coolant_in,
            )
            po.ovarrf(
                self.outfile,
                "Outlet temperature of blanket coolant (K)",
                "(temp_blkt_coolant_out)",
                fwbs_variables.temp_blkt_coolant_out,
            )
            po.ovarre(
                self.outfile,
                "Blanket (primary) coolant pressure (Pa)",
                "(pres_blkt_coolant)",
                fwbs_variables.pres_blkt_coolant,
            )
            if fwbs_variables.i_fw_blkt_shared_coolant == 1:
                po.ovarre(
                    self.outfile,
                    "Blanket coolant mass flow rate (kg/s)",
                    "(mflow_blkt_coolant_total)",
                    blanket_library.mflow_blkt_coolant_total,
                    "OP ",
                )

            # Total primary coolant mass flow rate (if they are the same coolant)
            if fwbs_variables.i_fw_blkt_shared_coolant == 0:
                po.ovarre(
                    self.outfile,
                    "Total (FW+BB) primary coolant mass flow rate(kg/s)",
                    "(mftotal)",
                    blanket_library.mftotal,
                    "OP ",
                )

            # BB Liquid Metal Breeder !
            if fwbs_variables.i_blkt_dual_coolant > 0:
                po.osubhd(self.outfile, "Breeding Blanket (breeder): ")

                po.ovarin(
                    self.outfile,
                    "Blanket liquid breeder type (0=PbLi, 1=Li)",
                    "(i_blkt_liquid_breeder_type)",
                    fwbs_variables.i_blkt_liquid_breeder_type,
                )
                if fwbs_variables.i_blkt_dual_coolant == 2:
                    po.ocmmnt(self.outfile, "Dual-coolant BB, i.e. self-cooled breeder.")
                    po.ovarrf(
                        self.outfile,
                        "Inlet temperature of blanket liquid breeder (K)",
                        "(inlet_temp_liq)",
                        fwbs_variables.inlet_temp_liq,
                    )
                    po.ovarrf(
                        self.outfile,
                        "Outlet temperature of blanket liquid breeder (K)",
                        "(outlet_temp_liq)",
                        fwbs_variables.outlet_temp_liq,
                    )
                    po.ovarre(
                        self.outfile,
                        "Blanket liquid breeder pressure (Pa)",
                        "(blpressure_liq)",
                        fwbs_variables.blpressure_liq,
                    )
                else:
                    po.ocmmnt(
                        self.outfile,
                        "single-coolant BB, breeder circulated for tritium extraction.",
                    )

                po.ovarre(
                    self.outfile,
                    "Blanket liquid breeder mass flow rate (kg/s)",
                    "(mfblkt_liq)",
                    blanket_library.mfblkt_liq,
                    "OP ",
                )

            # Pumping Power
            po.osubhd(self.outfile, "Mechanical pumping power: ")

            if fwbs_variables.i_fw_blkt_shared_coolant == 1:
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for FW (MW)",
                    "(p_fw_coolant_pump_mw)",
                    heat_transport_variables.p_fw_coolant_pump_mw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for blanket (primary) coolant (MW)",
                    "(p_blkt_coolant_pump_mw)",
                    heat_transport_variables.p_blkt_coolant_pump_mw,
                    "OP ",
                )
            if fwbs_variables.i_blkt_dual_coolant > 0:
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for blanket liquid breeder (MW)",
                    "(p_blkt_breeder_pump_mw)",
                    heat_transport_variables.p_blkt_breeder_pump_mw,
                    "OP ",
                )
            po.ovarre(
                self.outfile,
                "Total mechanical pumping power for FW and blanket (MW)",
                "(p_fw_blkt_coolant_pump_mw)",
                primary_pumping_variables.p_fw_blkt_coolant_pump_mw,
                "OP ",
            )
            if fwbs_variables.i_blkt_dual_coolant > 0:
                po.ovarre(
                    self.outfile,
                    "Total mechanical pumping power for FW, blanket and liquid metal breeder(MW)",
                    "(htpmw_blkt_tot)",
                    heat_transport_variables.htpmw_blkt_tot,
                    "OP ",
                )
            po.ovarre(
                self.outfile,
                "Pumping power for divertor (MW)",
                "(p_div_coolant_pump_mw)",
                heat_transport_variables.p_div_coolant_pump_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Pumping power for shield and vacuum vessel (MW)",
                "(p_shld_coolant_pump_mw)",
                heat_transport_variables.p_shld_coolant_pump_mw,
                "OP ",
            )

    def total_pressure_drop(
        self,
        output: bool,
        icoolpump: int,
        vel_coolant: float,
        len_pipe: float,
        n_pipe_90_deg_bends: int,
        n_pipe_180_deg_bends: int,
        den_coolant: float,
        visc_coolant_dynamic: float,
        coolant_electrical_conductivity: float,
        pol_channel_length: float,
        nopolchan: int,
        label: str,
    ) -> float:
        """
        Calculate the total pressure drop (Pa) for coolant flow in the first wall (FW) and breeding blanket (BZ).

        This includes frictional losses and, for liquid breeder coolants, magnetohydrodynamic (MHD) losses.

        :param output: Whether to write output to file.
        :type output: bool
        :param icoolpump: Switch for coolant type (1=primary He/H2O, 2=secondary PbLi/Li).
        :type icoolpump: int
        :param flow_velocity: Coolant flow velocity (m/s).
        :type flow_velocity: float
        :param len_pipe: Total flow length along pipe (m).
        :type len_pipe: float
        :param n_pipe_90_deg_bends: Number of 90 degree bends in pipe.
        :type n_pipe_90_deg_bends: int
        :param n_pipe_180_deg_bends: Number of 180 degree bends in pipe.
        :type n_pipe_180_deg_bends: int
        :param den_coolant: Coolant density (kg/m³).
        :type den_coolant: float
        :param visc_coolant_dynamic: Coolant dynamic viscosity (Pa s).
        :type visc_coolant_dynamic: float
        :param coolant_electrical_conductivity: Coolant electrical conductivity (A V⁻¹ m⁻¹).
        :type coolant_electrical_conductivity: float
        :param pol_channel_length: Length of poloidal channel section (m).
        :type pol_channel_length: float
        :param nopolchan: Number of poloidal channel sections.
        :type nopolchan: int
        :param label: Description label for output.
        :type label: str
        :return: Total pressure drop (Pa).
        :rtype: float
        """

        radius_pipe_90_deg_bend, radius_pipe_180_deg_bend = (
            self.calculate_pipe_bend_radius(i_ps=icoolpump)
        )

        # Friction - for all coolants
        dpres_friction = self.coolant_friction_pressure_drop(
            i_ps=icoolpump,
            radius_pipe_90_deg_bend=radius_pipe_90_deg_bend,
            radius_pipe_180_deg_bend=radius_pipe_180_deg_bend,
            n_pipe_90_deg_bends=n_pipe_90_deg_bends,
            n_pipe_180_deg_bends=n_pipe_180_deg_bends,
            len_pipe=len_pipe,
            den_coolant=den_coolant,
            visc_coolant=visc_coolant_dynamic,
            vel_coolant=vel_coolant,
            label=label,
            output=output,
        )

        if icoolpump == 2:
            dpres_mhd = self.liquid_breeder_mhd_pressure_drop(
                vel_coolant,
                visc_coolant_dynamic,
                coolant_electrical_conductivity,
                pol_channel_length,
                nopolchan,
                label,
                output=output,
            )
        else:
            dpres_mhd = 0

        # Total pressure drop (Pa)
        dpres_total = dpres_friction + dpres_mhd

        if output:
            po.osubhd(self.outfile, f"Total pressure drop for {label}")

            po.ocmmnt(self.outfile, "Friction drops plus MHD drops if applicaple")
            po.ovarre(
                self.outfile, "Total pressure drop (Pa)", "(deltap)", dpres_total, "OP "
            )
            po.ovarre(
                self.outfile,
                "Coolant flow velocity (m/s)",
                "(flow_velocity, formerly vv)",
                vel_coolant,
                "OP ",
            )

        return dpres_total

    def liquid_breeder_mhd_pressure_drop(
        self,
        vel: float,
        vsc: float,
        conduct_liq: float,
        l_channel: float,
        num_pol: int,
        label: str,
        output: bool = False,
    ):
        """Calculates the pressure drop in a liquid metal flow channel due to MHD effects. The total pressure
        drop is the sum of contributions. This is only used for secondary coolant/breeder so rectangular flow
        channels are assumed.

        author: G Graham, CCFE

        :param vel: liquid metal coolant/breeder flow velocity (m/s)
        :param vsc: liquid metal visosity
        :param conduct_liq: liquid metal conductivity
        :param l_channel: length long poloidal sections of channel
        :param num_pol: number long poloidal sections of channel
        :param label: description of calculation

        References:

             [Miy1986]   Miyazaki et al. (1986), Magneto-Hydro-Dynamic Pressure Drop of Lithium
                         Flow in Rectangular Ducts, Fusion Technology, 10:3P2A, 830-836, DOI: 10.13182/FST10-830

             [Mal1995]   Malang and Mattas (1995), Comparison of lithium and the eutectic
                         lead-lithium alloy, two candidate liquid metal breeder materials
                         for self-cooled blankets, Fusion Engineering and Design 27, 399-406

             [Iba2013]   Ibano et al (2013), Nutronics and pumping power analysis on the
                         Tokamak reactor for the fusion-biomass hybrid concept,
                         Fusion Engineering and Design, 88

             [Sho2018]   Shoki et al (2018), MHD pressure drop measurement of PbLi flow
                         in double-bended pipe, Fusion Engineering and Design, 136, 17-23

             [Klu2019]   Kluber et al. (2019), Numerical simulations of 3D magnetohydrodynamic
                         flows in dual-coolant lead lithium blankets, Fusion Engineering and Design,
                         146, 684-687

             [Sua2021]   MHD effects in geometrical sigularities on high velocity breeding
                         blanket designs. Part II, ENR-PRD.BB-T007-D002, EFDA_D_2PDT9U.
                         Also, see asssociated paper: Suarez et al. (2021), On the use of CFD
                         to obtain head loss coefficients in hydraulic systems and it's appliaction
                         to liquid metal flows in nuclear fusion reactor blankets, Plasma. Phys.
                         Control fusion, 63, 124002
        """
        # Magnetic feild strength in IB or OB blanket
        if label == "Inboard blanket breeder liquid":
            b_mag = fwbs_variables.b_mag_blkt[0]  # IB
        if label == "Outboard blanket breeder liquid":
            b_mag = fwbs_variables.b_mag_blkt[1]  # OB

        # Half-widths
        # N.B. a_bz_liq (width in the toroidal direction) is in B direction
        half_wth_a = fwbs_variables.a_bz_liq * 0.5
        half_wth_b = fwbs_variables.b_bz_liq * 0.5

        # If have thin conducting walls...
        if fwbs_variables.i_blkt_liquid_breeder_channel_type != 1:
            # Caculate resistances of fluid and walls
            r_i = half_wth_b / (conduct_liq * half_wth_a)
            r_w = half_wth_b / (
                fwbs_variables.bz_channel_conduct_liq * fwbs_variables.th_wall_secondary
            )
            big_c = r_i / r_w
            #  Calculate pressure drop for conducting wall [Miy1986]
            kp = big_c / (1 + half_wth_a / (3 * half_wth_b) + big_c)
            mhd_pressure_drop = kp * conduct_liq * vel * (b_mag**2) * l_channel

        # If have perfcetly insulating FCIs...
        else:
            # Calculate pressure drop for (perfectly) insulating FCI [Mal1995]
            mhd_pressure_drop = (
                vel * b_mag * l_channel * np.sqrt(conduct_liq * vsc / half_wth_a)
            )

        # Total (Pa)
        liquid_breeder_pressure_drop_mhd = num_pol * mhd_pressure_drop

        if output:
            po.osubhd(
                self.outfile,
                f"Liquid metal breeder/coolant MHD pressure drop for {label}",
            )

            if fwbs_variables.i_blkt_liquid_breeder_channel_type == 0:
                po.ocmmnt(
                    self.outfile,
                    "Flow channels have thin conducting walls (i_blkt_liquid_breeder_channel_type==0)",
                )
                po.ovarre(
                    self.outfile,
                    "Wall conductance (A V-1 m-1)",
                    "(bz_channel_conduct_liq)",
                    fwbs_variables.bz_channel_conduct_liq,
                    "OP ",
                )
            elif fwbs_variables.i_blkt_liquid_breeder_channel_type == 2:
                po.ocmmnt(
                    self.outfile,
                    "Flow Channel Inserts (FCIs) used (i_blkt_liquid_breeder_channel_type==2)",
                )
                po.ovarre(
                    self.outfile,
                    "FCI conductance (A V-1 m-1)",
                    "(bz_channel_conduct_liq)",
                    fwbs_variables.bz_channel_conduct_liq,
                    "OP ",
                )
            else:
                po.ocmmnt(
                    self.outfile,
                    "Flow Channel Inserts - assumed perfect insulator (i_blkt_liquid_breeder_channel_type==1)",
                )

            po.ovarre(
                self.outfile,
                "Length of long poloidal secion of channel (m)",
                "(l_channel)",
                l_channel,
                "OP ",
            )
            po.ovarin(
                self.outfile,
                "Number of long poloidal secions of channel",
                "(num_pol)",
                num_pol,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "MHD pressure drop (Pa)",
                "(liquid_breeder_pressure_drop_mhd)",
                liquid_breeder_pressure_drop_mhd,
                "OP ",
            )

        return liquid_breeder_pressure_drop_mhd

    def calculate_pipe_bend_radius(self, i_ps: int):
        """Set the pipe bend radius based on the coolant type.

        :param i_ps: switch for primary or secondary coolant
        """
        # If primary coolant or secondary coolant (See DCLL)
        radius_pipe_90_deg_bend = (
            (3 * fwbs_variables.radius_fw_channel)
            if (i_ps == 1)
            else fwbs_variables.b_bz_liq
        )
        radius_pipe_180_deg_bend = radius_pipe_90_deg_bend / 2

        return radius_pipe_90_deg_bend, radius_pipe_180_deg_bend

    def coolant_friction_pressure_drop(
        self,
        i_ps: int,
        radius_pipe_90_deg_bend: float,
        radius_pipe_180_deg_bend: float,
        n_pipe_90_deg_bends: float,
        n_pipe_180_deg_bends: float,
        len_pipe: float,
        den_coolant: float,
        visc_coolant: float,
        vel_coolant: float,
        label: str,
        output: bool = False,
    ):
        """
        Pressure drops are calculated for a pipe with a number of 90
        and 180 degree bends. The pressure drop due to frictional forces along
        the total straight length of the pipe is calculated, then the pressure
        drop due to the bends is calculated. The total pressure drop is the sum
        of all contributions.

        :param i_ps: switch for primary or secondary coolant
        :param radius_pipe_90_deg_bend: radius of 90 degree bend in pipe (m)
        :param radius_pipe_180_deg_bend: radius of 180 degree bend in pipe (m)
        :param n_pipe_90_deg_bends: number of 90 degree bends in the pipe
        :param n_pipe_180_deg_bends: number of 180 degree bends in the pipe
        :param len_pipe: total flow length along pipe (m)
        :param den_coolant: coolant density (kg/m³)
        :param visc_coolant: coolant viscosity (Pa s)
        :param vel_coolant: coolant flow velocity (m/s)
        :param label: component name
        :param output: boolean of whether to write data to output file

        :Notes:
            Darcy-Weisbach Equation (straight pipe):

            ΔP = λ * L/D * (p 〈v〉²) / 2

            λ - Darcy friction factor, L - pipe length, D - hydraulic diameter,
            p - fluid density, 〈v〉 - fluid flow average velocity

            This function also calculates pressure drop equations for elbow bends,
            with modified coefficients.

            N.B. Darcy friction factor is estimated from the Haaland approximation.
        """

        # Calculate hydraulic dimater for round or retancular pipe (m)
        dia_pipe = self.pipe_hydraulic_diameter(i_ps)

        # Reynolds number
        reynolds_number = den_coolant * vel_coolant * dia_pipe / visc_coolant

        # Calculate Darcy friction factor
        # N.B. friction function Uses Haaland approx. which assumes a filled circular pipe.
        # Use dh which allows us to do fluid calculations for non-cicular tubes
        # (dh is estimate appropriate for fully developed flow).

        darcy_friction_factor = self.fw.darcy_friction_haaland(
            reynolds_number,
            fwbs_variables.roughness_fw_channel,
            fwbs_variables.radius_fw_channel,
        )

        # Pressure drop coefficient

        # Straight section
        f_straight = darcy_friction_factor * len_pipe / dia_pipe

        # 90 degree elbow pressure drop coefficient
        f_elbow_90 = self.elbow_coeff(
            radius_pipe_elbow=radius_pipe_90_deg_bend,
            deg_pipe_elbow=90.0,
            darcy_friction=darcy_friction_factor,
            dia_pipe=dia_pipe,
        )

        # 180 degree elbow pressure drop coefficient
        f_elbow_180 = self.elbow_coeff(
            radius_pipe_elbow=radius_pipe_180_deg_bend,
            deg_pipe_elbow=180.0,
            darcy_friction=darcy_friction_factor,
            dia_pipe=dia_pipe,
        )

        # Pressure drop due to friction in straight sections
        dpres_straight = f_straight * 0.5 * den_coolant * vel_coolant**2

        # Pressure drop due to 90 and 180 degree bends
        dpres_90 = n_pipe_90_deg_bends * f_elbow_90 * 0.5 * den_coolant * vel_coolant**2
        dpres_180 = (
            n_pipe_180_deg_bends * f_elbow_180 * 0.5 * den_coolant * vel_coolant**2
        )

        # Total pressure drop (Pa)
        dpres_total = dpres_straight + dpres_90 + dpres_180

        if output:
            po.osubhd(self.outfile, f"Pressure drop (friction) for {label}")
            po.ovarre(self.outfile, "Reynolds number", "(reyn)", reynolds_number, "OP ")
            po.ovarre(
                self.outfile,
                "Darcy friction factor",
                "(lambda)",
                darcy_friction_factor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Pressure drop (Pa)",
                "(pressure_drop)",
                dpres_total,
                "OP ",
            )
            po.ocmmnt(self.outfile, "This is the sum of the following:")
            po.ovarre(
                self.outfile,
                "            Straight sections (Pa)",
                "(pdropstraight)",
                dpres_straight,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "            90 degree bends (Pa)",
                "(pdrop90)",
                dpres_90,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "            180 degree bends (Pa)",
                "(pdrop180)",
                dpres_180,
                "OP ",
            )

            # TN: always write verbose stuff, it has no harm
            po.ovarre(
                self.outfile,
                "Straight section pressure drop coefficient",
                "(kstrght)",
                f_straight,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "90 degree elbow coefficient",
                "(kelbwn)",
                f_elbow_90,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "180 degree elbow coefficient coefficient",
                "(kelbwt)",
                f_elbow_180,
                "OP ",
            )

        return dpres_total

    def pipe_hydraulic_diameter(self, i_channel_shape):
        """Caculate the hydraulic diameter (m) for a given coolant pipe size/shape.
        author: G. Graham

        :param i_channel_shape: switch for circular or rectangular channel crossection.
            Shape depends on whether primary or secondary coolant
        """
        # If primary coolant then circular channels assumed
        if i_channel_shape == 1:
            return 2.0 * fwbs_variables.radius_fw_channel

        # If secondary coolant then rectangular channels assumed
        if i_channel_shape == 2:
            return (
                2
                * fwbs_variables.a_bz_liq
                * fwbs_variables.b_bz_liq
                / (fwbs_variables.a_bz_liq + fwbs_variables.b_bz_liq)
            )

        raise ProcessValueError(
            f"i_channel_shape ={i_channel_shape} is an invalid option."
        )

    def elbow_coeff(
        self,
        radius_pipe_elbow: float,
        deg_pipe_elbow: float,
        darcy_friction: float,
        dia_pipe: float,
    ) -> float:
        """
        Calculates elbow bend coefficients for pressure drop calculations.

        :param radius_pipe_elbow: Pipe elbow radius (m)
        :type radius_pipe_elbow: float
        :param deg_pipe_elbow: Pipe elbow angle (degrees)
        :type deg_pipe_elbow: float
        :param darcy_friction: Darcy friction factor
        :type darcy_friction: float
        :param dia_pipe: Pipe diameter (m)
        :type dia_pipe: float
        :return: Elbow coefficient for pressure drop calculation
        :rtype: float

        :References:
        - [Ide1969] Idel'Cik, I. E. (1969), Memento des pertes de charge,
          Collection de la Direction des Etudes et Recherches d'Electricité de France.

        """

        if deg_pipe_elbow == 90:
            a = 1.0
        elif deg_pipe_elbow < 70:
            a = 0.9 * np.sin(deg_pipe_elbow * np.pi / 180.0)
        elif deg_pipe_elbow > 100:
            a = 0.7 + (0.35 * np.sin((deg_pipe_elbow / 90.0) * (np.pi / 180.0)))
        else:
            raise ProcessValueError(
                "No formula for 70 <= elbow angle(deg) <= 100, only 90 deg option available in this range."
            )

        r_ratio = radius_pipe_elbow / dia_pipe

        if r_ratio > 1:
            b = 0.21 / r_ratio**0.5
        elif r_ratio < 1:
            b = 0.21 / r_ratio**2.5
        else:
            b = 0.21

        # Singularity
        ximt = a * b

        # Friction
        xift = (
            (np.pi / 180.0)
            * darcy_friction
            * (radius_pipe_elbow / dia_pipe)
            * deg_pipe_elbow
        )

        # Elbow Coefficient
        return ximt + xift

    def coolant_pumping_power(
        self,
        output: bool,
        i_liquid_breeder: int,
        temp_coolant_pump_outlet: float,
        temp_coolant_pump_inlet: float,
        pres_coolant_pump_inlet: float,
        dpres_coolant: float,
        mflow_coolant_total: float,
        primary_coolant_switch: str,
        den_coolant: float,
        label: str,
    ) -> float:
        """
        Calculate the coolant pumping power in MW for the first wall (FW) or breeding blanket (BZ) coolant.

        :param output: Whether to write data to output file.
        :type output: bool
        :param i_liquid_breeder: Switch for primary coolant or secondary coolant/breeder (1=primary He/H2O, 2=secondary PbLi/Li).
        :type i_liquid_breeder: int
        :param temp_coolant_pump_outlet: Pump outlet temperature (K).
        :type temp_coolant_pump_outlet: float
        :param temp_coolant_pump_inlet: Pump inlet temperature (K).
        :type temp_coolant_pump_inlet: float
        :param pressure: Outlet (pump inlet) coolant pressure (Pa).
        :type pressure: float
        :param dpres_coolant: Coolant pressure drop (Pa).
        :type dpres_coolant: float
        :param mflow_coolant_total: Total coolant mass flow rate in (kg/s).
        :type mflow_coolant_total: float
        :param primary_coolant_switch: Name of FW/blanket coolant (e.g., "Helium" or "Water") if icoolpump=1.
        :type primary_coolant_switch: str
        :param den_coolant: Density of coolant or liquid breeder (kg/m³).
        :type den_coolant: float
        :param label: Description label for output.
        :type label: str
        :return: Pumping power in MW.
        :rtype: float

        :references:
            - Idel'Cik, I. E. (1969), Memento des pertes de charge
            - S.P. Sukhatme (2005), A Textbook on Heat Transfer
        """

        # Pump outlet pressure (Pa)
        # The pump adds the pressure lost going through the coolant channels back
        pres_coolant_pump_outlet = pres_coolant_pump_inlet + dpres_coolant

        # Adiabatic index for helium or water
        gamma = (5 / 3) if fwbs_variables.i_blkt_coolant_type == 1 else (4 / 3)

        # If calculating for primary coolant
        if i_liquid_breeder == 1:
            # The pumping power is be calculated in the most general way,
            # using enthalpies before and after the pump.

            pump_outlet_fluid_properties = FluidProperties.of(
                fluid_name=primary_coolant_switch,
                temperature=temp_coolant_pump_outlet,
                pressure=pres_coolant_pump_outlet,
            )

            # Assume isentropic pump so that s1 = s2
            s1 = pump_outlet_fluid_properties.entropy

            # Get specific enthalpy at the outlet (J/kg) before pump using pressure and entropy s1
            pump_inlet_fluid_properties = FluidProperties.of(
                fluid_name=primary_coolant_switch,
                pressure=pres_coolant_pump_inlet,
                entropy=s1,
            )

            # Pumping power (MW) is given by enthalpy change, with a correction for
            # the isentropic efficiency of the pump.
            fp = (
                temp_coolant_pump_outlet
                * (
                    1
                    - (pres_coolant_pump_outlet / pres_coolant_pump_inlet)
                    ** -((gamma - 1) / gamma)
                )
                / (
                    fwbs_variables.etaiso
                    * (temp_coolant_pump_inlet - temp_coolant_pump_outlet)
                )
            )
            pumppower = (
                1e-6
                * mflow_coolant_total
                * (
                    pump_outlet_fluid_properties.enthalpy
                    - pump_inlet_fluid_properties.enthalpy
                )
                / fwbs_variables.etaiso
            ) / (1 - fp)

        # If calculating for secondary coolant/breeder...
        else:
            # Calculate specific volume
            spec_vol = 1 / den_coolant

            # Pumping power (MW) is given by pressure change, with a correction for
            # the isentropic efficiency of the pump.
            fp = (
                temp_coolant_pump_outlet
                * (
                    1
                    - (pres_coolant_pump_outlet / pres_coolant_pump_inlet)
                    ** -((gamma - 1) / gamma)
                )
                / (
                    fwbs_variables.etaiso_liq
                    * (temp_coolant_pump_inlet - temp_coolant_pump_outlet)
                )
            )
            pumppower = (
                1e-6
                * mflow_coolant_total
                * spec_vol
                * dpres_coolant
                / fwbs_variables.etaiso_liq
            ) / (1 - fp)

        # Error for dpres_coolant too large
        if fp >= 1:
            raise ProcessValueError(
                "Pressure drops in coolant are too large to be feasible"
            )

        if output:
            po.oheadr(self.outfile, "Mechanical Pumping Power for " + label)
            po.osubhd(self.outfile, "Pumping power for " + label)

            po.ovarre(
                self.outfile, "Pumping power (MW)", "(pumppower)", pumppower, "OP "
            )
            po.ovarre(
                self.outfile,
                "FW or Blanket inlet (pump oulet) pressure (Pa)",
                "(coolpin)",
                pres_coolant_pump_outlet,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "FW or Blanket oulet (pump inlet) pressure (Pa)",
                "(pres_coolant_pump_inlet)",
                pres_coolant_pump_inlet,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "FW or Blanket total pressure drop (Pa)",
                "(dpres_coolant)",
                dpres_coolant,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Mass flow rate in (kg/s) = ",
                "(mf)",
                mflow_coolant_total,
                "OP ",
            )

        return pumppower


def set_pumping_powers_as_fractions(
    f_p_fw_coolant_pump_total_heat: float,
    f_p_blkt_coolant_pump_total_heat: float,
    f_p_shld_coolant_pump_total_heat: float,
    f_p_div_coolant_pump_total_heat: float,
    p_fw_nuclear_heat_total_mw: float,
    psurffwi: float,
    psurffwo: float,
    p_blkt_nuclear_heat_total_mw: float,
    p_shld_nuclear_heat_mw: float,
    p_cp_shield_nuclear_heat_mw: float,
    p_plasma_separatrix_mw: float,
    p_div_nuclear_heat_total_mw: float,
    p_div_rad_total_mw: float,
) -> tuple[float, float, float, float]:
    """
    Calculate mechanical pumping powers as fractions of thermal power in each component.

    :param f_p_fw_coolant_pump_total_heat: Fraction for FW coolant pump.
    :type f_p_fw_coolant_pump_total_heat: float
    :param f_p_blkt_coolant_pump_total_heat: Fraction for blanket coolant pump.
    :type f_p_blkt_coolant_pump_total_heat: float
    :param f_p_shld_coolant_pump_total_heat: Fraction for shield coolant pump.
    :type f_p_shld_coolant_pump_total_heat: float
    :param f_p_div_coolant_pump_total_heat: Fraction for divertor coolant pump.
    :type f_p_div_coolant_pump_total_heat: float
    :param p_fw_nuclear_heat_total_mw: Total FW nuclear heating (MW).
    :type p_fw_nuclear_heat_total_mw: float
    :param psurffwi: Inboard FW surface heating (MW).
    :type psurffwi: float
    :param psurffwo: Outboard FW surface heating (MW).
    :type psurffwo: float
    :param p_blkt_nuclear_heat_total_mw: Total blanket nuclear heating (MW).
    :type p_blkt_nuclear_heat_total_mw: float
    :param p_shld_nuclear_heat_mw: Shield nuclear heating (MW).
    :type p_shld_nuclear_heat_mw: float
    :param p_cp_shield_nuclear_heat_mw: CP shield nuclear heating (MW).
    :type p_cp_shield_nuclear_heat_mw: float
    :param p_plasma_separatrix_mw: Plasma separatrix power (MW).
    :type p_plasma_separatrix_mw: float
    :param p_div_nuclear_heat_total_mw: Divertor nuclear heating (MW).
    :type p_div_nuclear_heat_total_mw: float
    :param p_div_rad_total_mw: Divertor radiative power (MW).
    :type p_div_rad_total_mw: float

    :return: Tuple of pumping powers (MW) for FW, blanket, shield, and divertor.
    :rtype: tuple[float, float, float, float]
    """
    p_fw_coolant_pump_mw = f_p_fw_coolant_pump_total_heat * (
        p_fw_nuclear_heat_total_mw + psurffwi + psurffwo
    )
    p_blkt_coolant_pump_mw = (
        f_p_blkt_coolant_pump_total_heat * p_blkt_nuclear_heat_total_mw
    )
    p_shld_coolant_pump_mw = f_p_shld_coolant_pump_total_heat * (
        p_shld_nuclear_heat_mw + p_cp_shield_nuclear_heat_mw
    )
    p_div_coolant_pump_mw = f_p_div_coolant_pump_total_heat * (
        p_plasma_separatrix_mw + p_div_nuclear_heat_total_mw + p_div_rad_total_mw
    )
    return (
        p_fw_coolant_pump_mw,
        p_blkt_coolant_pump_mw,
        p_shld_coolant_pump_mw,
        p_div_coolant_pump_mw,
    )


def eshellarea(rshell, rmini, rmino, zminor):
    """Routine to calculate the inboard, outboard and total surface areas
    of a toroidal shell comprising two elliptical sections
    author: P J Knight, CCFE, Culham Science Centre
    rshell : input real : major radius of centre of both ellipses (m)
    rmini  : input real : horizontal distance from rshell to
    inboard elliptical shell (m)
    rmino  : input real : horizontal distance from rshell to
    outboard elliptical shell (m)
    zminor : input real : vertical internal half-height of shell (m)
    ain    : output real : surface area of inboard section (m3)
    aout   : output real : surface area of outboard section (m3)
    atot   : output real : total surface area of shell (m3)
    This routine calculates the surface area of the inboard and outboard
    sections of a toroidal shell defined by two co-centred semi-ellipses.
    """

    # Inboard section
    elong = zminor / rmini
    ain = 2.0 * np.pi * elong * (np.pi * rshell * rmini - 2.0 * rmini * rmini)

    # Outboard section
    elong = zminor / rmino
    aout = 2.0 * np.pi * elong * (np.pi * rshell * rmino + 2.0 * rmino * rmino)

    return ain, aout, ain + aout


def dshellarea(
    rmajor: float, rminor: float, zminor: float
) -> tuple[float, float, float]:
    """
    Calculate the inboard, outboard, and total surface areas of a D-shaped toroidal shell.

    :param rmajor: Major radius of inboard straight section (m)
    :type rmajor: float
    :param rminor: Horizontal width of shell (m)
    :type rminor: float
    :param zminor: Vertical half-height of shell (m)
    :type zminor: float

    :return: Tuple containing:
        - ain: Surface area of inboard straight section (m²)
        - aout: Surface area of outboard curved section (m²)
        - atot: Total surface area of shell (m²)
    :rtype: tuple[float, float, float]

    The inboard section is assumed to be a cylinder.
    The outboard section is defined by a semi-ellipse, centred on the major radius of the inboard section.
    """
    # Area of inboard cylindrical shell
    ain = 4.0 * zminor * np.pi * rmajor

    # Area of elliptical outboard section
    elong = zminor / rminor
    aout = 2.0 * np.pi * elong * (np.pi * rmajor * rminor + 2.0 * rminor * rminor)

    return ain, aout, ain + aout


def eshellvol(rshell, rmini, rmino, zminor, drin, drout, dz):
    """Routine to calculate the inboard, outboard and total volumes
    of a toroidal shell comprising two elliptical sections
    author: P J Knight, CCFE, Culham Science Centre
    rshell : input real : major radius of centre of both ellipses (m)
    rmini  : input real : horizontal distance from rshell to outer edge
    of inboard elliptical shell (m)
    rmino  : input real : horizontal distance from rshell to inner edge
    of outboard elliptical shell (m)
    zminor : input real : vertical internal half-height of shell (m)
    drin   : input real : horiz. thickness of inboard shell at midplane (m)
    drout  : input real : horiz. thickness of outboard shell at midplane (m)
    dz     : input real : vertical thickness of shell at top/bottom (m)
    vin    : output real : volume of inboard section (m3)
    vout   : output real : volume of outboard section (m3)
    vtot   : output real : total volume of shell (m3)
    This routine calculates the volume of the inboard and outboard sections
    of a toroidal shell defined by two co-centred semi-ellipses.
    Each section's internal and external surfaces are in turn defined
    by two semi-ellipses. The volumes of each section are calculated as
    the difference in those of the volumes of revolution enclosed by their
    inner and outer surfaces.
    """
    # Inboard section

    # Volume enclosed by outer (higher R) surface of elliptical section
    # and the vertical straight line joining its ends
    a = rmini
    b = zminor
    elong = b / a
    v1 = 2.0 * np.pi * elong * (0.5 * np.pi * rshell * a**2 - (2.0 / 3.0) * a**3)

    # Volume enclosed by inner (lower R) surface of elliptical section
    # and the vertical straight line joining its ends
    a = rmini + drin
    b = zminor + dz
    elong = b / a
    v2 = 2.0 * np.pi * elong * (0.5 * np.pi * rshell * a**2 - (2.0 / 3.0) * a**3)

    # Volume of inboard section of shell
    vin = v2 - v1

    # Outboard section

    # Volume enclosed by inner (lower R) surface of elliptical section
    # and the vertical straight line joining its ends
    a = rmino
    b = zminor
    elong = b / a
    v1 = 2.0 * np.pi * elong * (0.5 * np.pi * rshell * a**2 + (2.0 / 3.0) * a**3)

    # Volume enclosed by outer (higher R) surface of elliptical section
    # and the vertical straight line joining its ends
    a = rmino + drout
    b = zminor + dz
    elong = b / a
    v2 = 2.0 * np.pi * elong * (0.5 * np.pi * rshell * a**2 + (2.0 / 3.0) * a**3)

    # Volume of outboard section of shell
    vout = v2 - v1

    return vin, vout, vin + vout


def dshellvol(rmajor, rminor, zminor, drin, drout, dz):
    """Routine to calculate the inboard, outboard and total volumes
    of a D-shaped toroidal shell
    author: P J Knight, CCFE, Culham Science Centre
    rmajor : input real : major radius to outer point of inboard
    straight section of shell (m)
    rminor : input real : horizontal internal width of shell (m)
    zminor : input real : vertical internal half-height of shell (m)
    drin   : input real : horiz. thickness of inboard shell at midplane (m)
    drout  : input real : horiz. thickness of outboard shell at midplane (m)
    dz     : input real : vertical thickness of shell at top/bottom (m)
    vin    : output real : volume of inboard straight section (m3)
    vout   : output real : volume of outboard curved section (m3)
    vtot   : output real : total volume of shell (m3)
    This routine calculates the volume of the inboard and outboard sections
    of a D-shaped toroidal shell defined by the above input parameters.
    The inboard section is assumed to be a cylinder of uniform thickness.
    The outboard section's internal and external surfaces are defined
    by two semi-ellipses, centred on the outer edge of the inboard section;
    its volume is calculated as the difference in those of the volumes of
    revolution enclosed by the two surfaces.
    """
    # Volume of inboard cylindrical shell
    vin = 2.0 * (zminor + dz) * np.pi * (rmajor**2 - (rmajor - drin) ** 2)

    # Volume enclosed by inner surface of elliptical outboard section
    # and the vertical straight line joining its ends
    a = rminor
    b = zminor
    elong = b / a
    v1 = 2.0 * np.pi * elong * (0.5 * np.pi * rmajor * a**2 + (2.0 / 3.0) * a**3)

    # Volume enclosed by outer surface of elliptical outboard section
    # and the vertical straight line joining its ends
    a = rminor + drout
    b = zminor + dz
    elong = b / a
    v2 = 2.0 * np.pi * elong * (0.5 * np.pi * rmajor * a**2 + (2.0 / 3.0) * a**3)

    # Volume of elliptical outboard shell
    vout = v2 - v1

    return vin, vout, vin + vout


class OutboardBlanket(BlanketLibrary):
    def calculate_basic_geometry(self):
        self.component_volumes()

        dia_blkt_channel = self.pipe_hydraulic_diameter(i_channel_shape=1)
        fwbs_variables.radius_blkt_channel = dia_blkt_channel / 2
        (
            fwbs_variables.radius_blkt_channel_90_bend,
            fwbs_variables.radius_blkt_channel_180_bend,
        ) = self.calculate_pipe_bend_radius(i_ps=1)

    def calculate_blanket_outboard_module_geometry(
        self,
        n_blkt_outboard_modules_toroidal: int,
        rmajor: float,
        rminor: float,
        dr_fw_plasma_gap_outboard: float,
    ) -> float:
        """
        Calculate the mid-plane toroidal circumference and segment length of the outboard blanket.

        :param n_blkt_outboard_modules_toroidal: Number of outboard blanket modules in the toroidal direction.
        :type n_blkt_outboard_modules_toroidal: int
        :param rmajor: Major radius (m).
        :type rmajor: float
        :param rminor: Minor radius (m).
        :type rminor: float
        :param dr_fw_plasma_gap_outboard: Outboard first wall to plasma gap (m).
        :type dr_fw_plasma_gap_outboard: float
        :return: Length of outboard blanket segment in the toroidal direction (m).
        :rtype: float
        """
        return (
            2.0 * np.pi * (rmajor + rminor + dr_fw_plasma_gap_outboard)
        ) / n_blkt_outboard_modules_toroidal


class InboardBlanket(BlanketLibrary):
    def calculate_basic_geometry(self):
        self.component_volumes()

        dia_blkt_channel = self.pipe_hydraulic_diameter(i_channel_shape=1)
        fwbs_variables.radius_blkt_channel = dia_blkt_channel / 2
        (
            fwbs_variables.radius_blkt_channel_90_bend,
            fwbs_variables.radius_blkt_channel_180_bend,
        ) = self.calculate_pipe_bend_radius(i_ps=1)

        self.set_blanket_module_geometry()

    def calculate_blanket_inboard_module_geometry(
        self,
        n_blkt_inboard_modules_toroidal: int,
        rmajor: float,
        rminor: float,
        dr_fw_plasma_gap_inboard: float,
    ) -> float:
        """
        Calculate the mid-plane toroidal circumference and segment length of the inboard blanket.

        :param n_blkt_inboard_modules_toroidal: Number of inboard blanket modules in the toroidal direction.
        :type n_blkt_inboard_modules_toroidal: int
        :param rmajor: Major radius (m).
        :type rmajor: float
        :param rminor: Minor radius (m).
        :type rminor: float
        :param dr_fw_plasma_gap_inboard: Inboard first wall to plasma gap (m).
        :type dr_fw_plasma_gap_inboard: float
        :return: Length of inboard blanket segment in the toroidal direction (m).
        :rtype: float
        """
        return (
            2.0 * np.pi * (rmajor + rminor + dr_fw_plasma_gap_inboard)
        ) / n_blkt_inboard_modules_toroidal
