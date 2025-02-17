"""This library contains routines that can be shared by the blanket modules used in PROCESS.

author: G Graham, CCFE, Culham Science Centre
"""

import numpy as np

from process import (
    process_output as po,
)
from process.coolprop_interface import FluidProperties
from process.fortran import (
    blanket_library,
    build_variables,
    buildings_variables,
    constants,
    divertor_variables,
    error_handling,
    fwbs_variables,
    heat_transport_variables,
    pfcoil_variables,
    physics_variables,
    primary_pumping_variables,
)
from process.fortran import (
    error_handling as eh,
)
from process.utilities.f2py_string_patch import f2py_compatible_to_string

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
        self.outfile = constants.nout

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
        blanket_library.hblnkt = self.component_half_height(icomponent=0)
        # Shield
        blanket_library.hshld = self.component_half_height(icomponent=1)
        # Vacuum Vessel
        blanket_library.hvv = self.component_half_height(icomponent=2)

        # D-shaped blanket and shield
        if physics_variables.itart == 1 or fwbs_variables.i_fw_blkt_vv_shape == 1:
            for icomponent in range(3):
                self.dshaped_component(icomponent)

        # Elliptical blanket and shield
        else:
            for icomponent in range(3):
                self.elliptical_component(icomponent)

            # This will fail the hts_REBCO and 2D_scan regression tests,
            # the number of VMCON iterations (nviter) is different.
            # Seems to be because in the blanket calculations (icomponent=0):
            # r2 = 1.3836567143743970 rather than old value of r2 = 1.3836567143743972,
            # r3 = 3.7009701431231936 rather than r3 = 3.7009701431231923.

        # Apply coverage factors to volumes and surface areas
        self.apply_coverage_factors()

        # Calculate cryostat geometry
        self.external_cryo_geometry()

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
                build_variables.hmax
                - build_variables.dz_shld_vv_gap
                - build_variables.dz_vv_lower
            )
        else:
            raise ValueError(f"{icomponent=} is invalid, it must be either 0,1,2")

        # Calculate component internal upper half-height (m)
        # If a double null machine then symmetric
        if physics_variables.idivrt == 2:
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
            # Vacuum Vessel
            if icomponent == 2:
                htop = (
                    htop + build_variables.dz_blkt_upper + build_variables.dz_shld_upper
                )

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
                build_variables.blareaib,
                build_variables.blareaob,
                build_variables.blarea,
            ) = dshellarea(r1, r2, blanket_library.hblnkt)
        if icomponent == 1:
            (
                build_variables.shareaib,
                build_variables.shareaob,
                build_variables.sharea,
            ) = dshellarea(r1, r2, blanket_library.hshld)

        # Calculate volumes, assuming 100% coverage
        if icomponent == 0:
            (
                fwbs_variables.vol_blkt_inboard,
                fwbs_variables.vol_blkt_outboard,
                fwbs_variables.vol_blkt_total,
            ) = dshellvol(
                r1,
                r2,
                blanket_library.hblnkt,
                build_variables.dr_blkt_inboard,
                build_variables.dr_blkt_outboard,
                build_variables.dz_blkt_upper,
            )
        elif icomponent == 1:
            (
                blanket_library.volshldi,
                blanket_library.volshldo,
                fwbs_variables.volshld,
            ) = dshellvol(
                r1,
                r2,
                blanket_library.hshld,
                build_variables.dr_shld_inboard,
                build_variables.dr_shld_outboard,
                build_variables.dz_shld_upper,
            )
        elif icomponent == 2:
            (
                blanket_library.volvvi,
                blanket_library.volvvo,
                fwbs_variables.vol_vv,
            ) = dshellvol(
                r1,
                r2,
                blanket_library.hvv,
                build_variables.dr_vv_inboard,
                build_variables.dr_vv_outboard,
                (build_variables.dz_vv_upper + build_variables.dz_vv_lower) / 2,
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
            r3 = (
                r3 - build_variables.dr_shld_outboard - build_variables.dr_blkt_outboard
            )

        # Calculate surface area, assuming 100% coverage
        if icomponent == 0:
            (
                build_variables.blareaib,
                build_variables.blareaob,
                build_variables.blarea,
            ) = eshellarea(r1, r2, r3, blanket_library.hblnkt)
        if icomponent == 1:
            (
                build_variables.shareaib,
                build_variables.shareaob,
                build_variables.sharea,
            ) = eshellarea(r1, r2, r3, blanket_library.hshld)

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
                blanket_library.hblnkt,
                build_variables.dr_blkt_inboard,
                build_variables.dr_blkt_outboard,
                build_variables.dz_blkt_upper,
            )
        if icomponent == 1:
            (
                blanket_library.volshldi,
                blanket_library.volshldo,
                fwbs_variables.volshld,
            ) = eshellvol(
                r1,
                r2,
                r3,
                blanket_library.hshld,
                build_variables.dr_shld_inboard,
                build_variables.dr_shld_outboard,
                build_variables.dz_shld_upper,
            )
        if icomponent == 2:
            (
                blanket_library.volvvi,
                blanket_library.volvvo,
                fwbs_variables.vol_vv,
            ) = eshellvol(
                r1,
                r2,
                r3,
                blanket_library.hvv,
                build_variables.dr_vv_inboard,
                build_variables.dr_vv_outboard,
                (build_variables.dz_vv_upper + build_variables.dz_vv_lower) / 2,
            )

    def apply_coverage_factors(self):
        """Apply coverage factors to volumes
        author: J. Morris, CCFE, Culham Science Centre
        Apply coverage factors to volumes
        """
        # Apply blanket coverage factors
        if physics_variables.idivrt == 2:
            # double null configuration
            build_variables.blareaob = (
                build_variables.blarea
                * (1.0 - 2.0 * fwbs_variables.fdiv - fwbs_variables.fhcd)
                - build_variables.blareaib
            )
        else:
            # single null configuration
            build_variables.blareaob = (
                build_variables.blarea
                * (1.0 - fwbs_variables.fdiv - fwbs_variables.fhcd)
                - build_variables.blareaib
            )

        build_variables.blarea = build_variables.blareaib + build_variables.blareaob

        fwbs_variables.vol_blkt_outboard = (
            fwbs_variables.vol_blkt_total
            * (1.0 - fwbs_variables.fdiv - fwbs_variables.fhcd)
            - fwbs_variables.vol_blkt_inboard
        )
        fwbs_variables.vol_blkt_total = (
            fwbs_variables.vol_blkt_inboard + fwbs_variables.vol_blkt_outboard
        )

        # Apply shield coverage factors
        build_variables.shareaib = fwbs_variables.fvolsi * build_variables.shareaib
        build_variables.shareaob = fwbs_variables.fvolso * build_variables.shareaob
        build_variables.sharea = build_variables.shareaib + build_variables.shareaob

        blanket_library.volshldi = fwbs_variables.fvolsi * blanket_library.volshldi
        blanket_library.volshldo = fwbs_variables.fvolso * blanket_library.volshldo
        fwbs_variables.volshld = blanket_library.volshldi + blanket_library.volshldo

        # Apply vacuum vessel coverage factor
        # moved from dshaped_* and elliptical_* to keep coverage factor
        # changes in the same location.
        fwbs_variables.vol_vv = fwbs_variables.fvoldw * fwbs_variables.vol_vv

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
            build_variables.hmax + build_variables.dr_tf_inboard
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

        # Vacuum vessel mass (kg)
        fwbs_variables.vvmass = fwbs_variables.vol_vv * fwbs_variables.denstl

        # Sum of internal vacuum vessel and cryostat masses (kg)
        fwbs_variables.dewmkg = (
            fwbs_variables.vol_vv + fwbs_variables.vol_cryostat
        ) * fwbs_variables.denstl

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
            f2py_compatible_to_string(fwbs_variables.i_fw_coolant_type).title()
            == "Helium"
            and fwbs_variables.i_blkt_coolant_type == 2
        ):
            fwbs_variables.i_fw_blkt_shared_coolant = 1
        if (
            f2py_compatible_to_string(fwbs_variables.i_fw_coolant_type).title()
            == "Water"
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
                f2py_compatible_to_string(fwbs_variables.i_fw_coolant_type),
                temperature=mid_temp,
                pressure=fwbs_variables.pres_fw_coolant.item(),
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
                f2py_compatible_to_string(fwbs_variables.i_fw_coolant_type),
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
            raise RuntimeError(
                f"Error in primary_coolant_properties. {fwbs_variables.den_fw_coolant = }"
            )
        if (
            fwbs_variables.den_blkt_coolant > 1e9
            or fwbs_variables.den_blkt_coolant <= 0
            or np.isnan(fwbs_variables.den_blkt_coolant)
        ):
            raise RuntimeError(
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

    def thermo_hydraulic_model_pressure_drop_calculations(self, output: bool):
        """
        Function that calculates the pressure drops for the thermo-hydraulic model
        when i_coolant_pumping = 2.

        Within are calculations necessary for the deltap_tot function but not required
        for other calculations within the thermo-hydraulic model as then they are just
        included there.

        Returns the pressure drops as a list with the number of entries dependent upon
        the switches icooldual and i_blkt_inboard.
        """
        npoltoti = 0
        npoltoto = 0
        npblkti_liq = 0
        npblkto_liq = 0

        if fwbs_variables.i_blanket_type == 5:
            # Unless DCLL then we will use BZ
            blanket_library.bldepti = build_variables.blbuith
            blanket_library.bldepto = build_variables.blbuoth
        else:
            blanket_library.bldepti = 0.8e0 * build_variables.dr_blkt_inboard
            blanket_library.bldepto = 0.8e0 * build_variables.dr_blkt_outboard

        # Using the total perimeter of the machine, segment the outboard
        # blanket into nblktmodp*nblktmodt modules, all assumed to be the same size

        # If SMS blanket then do not have seperate poloidal modules....
        # Should not need this as n_blkt_inboard_modules_poloidal is input but make sure here.
        if fwbs_variables.ims == 1:
            fwbs_variables.n_blkt_inboard_modules_poloidal = 1
            fwbs_variables.n_blkt_outboard_modules_poloidal = 1

        # Calculate mid-plane toroidal circumference and segment
        blanket_library.blwidti = (
            2.0e0
            * np.pi
            * (
                physics_variables.rmajor
                - physics_variables.rminor
                - build_variables.dr_fw_plasma_gap_inboard
            )
        ) / fwbs_variables.n_blkt_inboard_modules_toroidal
        blanket_library.blwidto = (
            2.0e0
            * np.pi
            * (
                physics_variables.rmajor
                + physics_variables.rminor
                + build_variables.dr_fw_plasma_gap_outboard
            )
        ) / fwbs_variables.n_blkt_outboard_modules_toroidal

        # Calculate poloidal height of blanket modules
        self.blanket_mod_pol_height()

        if fwbs_variables.icooldual > 0:
            # Use smallest space available to pipes for pipe sizes in pumping calculations (worst case)
            if fwbs_variables.i_blkt_inboard == 1:
                # Radial direction
                fwbs_variables.b_bz_liq = (
                    min(
                        (blanket_library.bldepti * fwbs_variables.r_f_liq_ib),
                        (blanket_library.bldepto * fwbs_variables.r_f_liq_ob),
                    )
                    / fwbs_variables.nopol
                )
                # Toroidal direction
                fwbs_variables.a_bz_liq = (
                    min(
                        (blanket_library.blwidti * fwbs_variables.w_f_liq_ib),
                        (blanket_library.blwidto * fwbs_variables.w_f_liq_ob),
                    )
                    / fwbs_variables.nopipes
                )
                # Poloidal
                if (blanket_library.bllengi < (fwbs_variables.b_bz_liq * 3)) or (
                    blanket_library.bllengo < (fwbs_variables.b_bz_liq * 3)
                ):
                    eh.report_error(278)

            # Unless there is no IB blanket...
            else:
                # Radial direction
                fwbs_variables.b_bz_liq = (
                    blanket_library.bldepto * fwbs_variables.r_f_liq_ob
                ) / fwbs_variables.nopol
                # Toroidal direction
                fwbs_variables.a_bz_liq = (
                    blanket_library.blwidto * fwbs_variables.w_f_liq_ob
                ) / fwbs_variables.nopipes
                # Poloidal
                if blanket_library.bllengo < (fwbs_variables.b_bz_liq * 3):
                    eh.report_error(278)

        # Calculate total flow lengths, used for pressure drop calculation
        # Blanket primary coolant flow
        blanket_library.bzfllengi = (
            fwbs_variables.bzfllengi_n_rad * blanket_library.bldepti
            + fwbs_variables.bzfllengi_n_pol * blanket_library.bllengi
        )
        blanket_library.bzfllengo = (
            fwbs_variables.bzfllengo_n_rad * blanket_library.bldepto
            + fwbs_variables.bzfllengo_n_pol * blanket_library.bllengo
        )
        # Blanket secondary coolant/breeder flow
        pollengi = blanket_library.bllengi
        pollengo = blanket_library.bllengo
        fwbs_variables.nopol = 2
        fwbs_variables.nopipes = 4
        bzfllengi_liq = (
            fwbs_variables.bzfllengi_n_rad_liq * blanket_library.bldepti
            + fwbs_variables.bzfllengi_n_pol_liq * blanket_library.bllengi
        )
        bzfllengo_liq = (
            fwbs_variables.bzfllengo_n_rad_liq * blanket_library.bldepto
            + fwbs_variables.bzfllengo_n_pol_liq * blanket_library.bllengo
        )

        # Coolant channel bends #########

        # Number of angle turns in FW and blanket flow channels, n.b. these are the
        # same for ccfe hcpb and kit hcll. FW is also be the same for DCLL MMS ans SMS.
        no90fw = 2
        no180fw = 0

        # N.B. This is for BZ only, does not include MF/BSS.
        if fwbs_variables.icooldual == 2 or fwbs_variables.icooldual == 1:
            no90bz = 4
            no180bz = 1
            no90bz_liq = 2
            no180bz_liq = 1
        else:
            no90bz = 4
            no180bz = 1

        # FW Pipe Flow and Velocity ######

        # Total number of first wall pipes from channel length and dx_fw_module (02/12/2015)
        blanket_library.npfwi = build_variables.a_fw_inboard / (
            fwbs_variables.len_fw_channel * fwbs_variables.dx_fw_module
        )
        blanket_library.npfwo = build_variables.a_fw_outboard / (
            fwbs_variables.len_fw_channel * fwbs_variables.dx_fw_module
        )

        # Mass flow rate per FW coolant pipe (kg/s):
        blanket_library.mffwpi = blanket_library.mffwi / blanket_library.npfwi
        blanket_library.mffwpo = blanket_library.mffwo / blanket_library.npfwo

        # Coolant velocite in FW (m/s)
        velfwi = self.flow_velocity(
            i_channel_shape=1,
            mass_flow_rate=blanket_library.mffwpi,
            flow_density=fwbs_variables.den_fw_coolant,
        )
        velfwo = self.flow_velocity(
            i_channel_shape=1,
            mass_flow_rate=blanket_library.mffwpo,
            flow_density=fwbs_variables.den_fw_coolant,
        )

        # If the blanket is dual-coolant...
        if fwbs_variables.icooldual == 2:
            # Calc total num of pipes (in all inboard modules) from
            # coolant frac and channel dimensions
            # Assumes up/down flow, two 90 deg bends per length
            blanket_library.npblkto = (
                fwbs_variables.vfblkt * fwbs_variables.vol_blkt_outboard
            ) / (
                np.pi
                * fwbs_variables.radius_fw_channel
                * fwbs_variables.radius_fw_channel
                * blanket_library.bzfllengo
            )
            npblkto_liq = (
                fwbs_variables.nopipes
                * fwbs_variables.n_blkt_outboard_modules_toroidal
                * fwbs_variables.n_blkt_outboard_modules_poloidal
            )

            # Mass flow rate per coolant pipe
            blanket_library.mfblktpo = blanket_library.mfblkto / blanket_library.npblkto
            mfblktpo_liq = blanket_library.mfblkto_liq / npblkto_liq
            # Coolant velocites in blanket (m/s)
            # Assume BZ structure has same channel width as FW
            blanket_library.velblkto = self.flow_velocity(
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
                blanket_library.npblkti = (
                    fwbs_variables.vfblkt * fwbs_variables.vol_blkt_inboard
                ) / (
                    np.pi
                    * fwbs_variables.radius_fw_channel
                    * fwbs_variables.radius_fw_channel
                    * blanket_library.bzfllengi
                )
                # Have DEMO DCLL set here for now
                npblkti_liq = (
                    fwbs_variables.nopipes
                    * fwbs_variables.n_blkt_inboard_modules_toroidal
                    * fwbs_variables.n_blkt_inboard_modules_poloidal
                )

                # Mass flow rate per coolant pipe
                blanket_library.mfblktpi = (
                    blanket_library.mfblkti / blanket_library.npblkti
                )
                blanket_library.mfblktpi_liq = blanket_library.mfblkti_liq / npblkti_liq

                # Coolant velocites in blanket (m/s)
                # Assume BZ structure has same channel width as FW
                blanket_library.velblkti = self.flow_velocity(
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
        elif fwbs_variables.icooldual == 1:
            # Calc total num of pipes (in all inboard modules) from
            # coolant frac and channel dimensions
            # Assumes up/down flow, two 90 deg bends per length
            blanket_library.npblkto = (
                fwbs_variables.vfblkt * fwbs_variables.vol_blkt_outboard
            ) / (
                np.pi
                * fwbs_variables.radius_fw_channel
                * fwbs_variables.radius_fw_channel
                * blanket_library.bzfllengo
            )
            npblkto_liq = (
                fwbs_variables.nopipes
                * fwbs_variables.n_blkt_outboard_modules_toroidal
                * fwbs_variables.n_blkt_outboard_modules_poloidal
            )

            # Mass flow rate per coolant pipe
            blanket_library.mfblktpo = blanket_library.mfblkto / blanket_library.npblkto

            # Coolant velocity in blanket (m/s)
            # Assume BZ structure has same channel width as FW
            blanket_library.velblkto = self.flow_velocity(
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
                blanket_library.npblkti = (
                    fwbs_variables.vfblkt * fwbs_variables.vol_blkt_inboard
                ) / (
                    np.pi
                    * fwbs_variables.radius_fw_channel
                    * fwbs_variables.radius_fw_channel
                    * blanket_library.bzfllengi
                )
                # Have DEMO DCLL set here for now
                npblkti_liq = (
                    fwbs_variables.nopipes
                    * fwbs_variables.n_blkt_inboard_modules_toroidal
                    * fwbs_variables.n_blkt_inboard_modules_poloidal
                )

                # Mass flow rate per coolant pipe
                blanket_library.mfblktpi = (
                    blanket_library.mfblkti / blanket_library.npblkti
                )

                # Coolant velocity in blanket (m/s)
                # Assume BZ structure has same channel width as FW
                blanket_library.velblkti = self.flow_velocity(
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
                blanket_library.mfblktpi_liq = fwbs_variables.mfblkti_liq / npblkti_liq
                velblkti_liq = self.flow_velocity(
                    i_channel_shape=2,
                    mass_flow_rate=blanket_library.mfblktpi_liq,
                    flow_density=fwbs_variables.den_liq,
                )

        # If the blanket is single-coolant with solid breeder...
        else:
            # Calculate total number of pipes (in all outboard modules) from coolant fraction and
            # channel dimensions (assumes up/down flow, two 90 deg bends per length)
            blanket_library.npblkto = (
                fwbs_variables.vfblkt * fwbs_variables.vol_blkt_outboard
            ) / (
                np.pi
                * fwbs_variables.radius_fw_channel
                * fwbs_variables.radius_fw_channel
                * blanket_library.bzfllengo
            )

            # Mass flow rate per coolant pipe
            blanket_library.mfblktpo = blanket_library.mfblkto / blanket_library.npblkto

            # Coolant velocity in blanket (m/s)
            # Assume BZ structure has same channel width as FW
            blanket_library.velblkto = self.flow_velocity(
                i_channel_shape=1,
                mass_flow_rate=blanket_library.mfblktpo,
                flow_density=fwbs_variables.den_blkt_coolant,
            )

            if fwbs_variables.i_blkt_inboard == 1:
                # Calc total num of pipes (in all inboard modules) from
                # coolant frac and channel dimensions
                # Assumes up/down flow, two 90 deg bends per length
                blanket_library.npblkti = (
                    fwbs_variables.vfblkt * fwbs_variables.vol_blkt_inboard
                ) / (
                    np.pi
                    * fwbs_variables.radius_fw_channel
                    * fwbs_variables.radius_fw_channel
                    * blanket_library.bzfllengi
                )

                # Mass flow rate per coolant pipe
                blanket_library.mfblktpi = (
                    blanket_library.mfblkti / blanket_library.npblkti
                )

                # Coolant velocity in blanket (m/s)
                # Assume BZ structure has same channel width as FW
                blanket_library.velblkti = self.flow_velocity(
                    i_channel_shape=1,
                    mass_flow_rate=blanket_library.mfblktpi,
                    flow_density=fwbs_variables.den_blkt_coolant,
                )

        # FW Presure Drops ###############

        deltap_fwi = self.deltap_tot(
            output,
            icoolpump=1,
            flow_velocity=velfwi,
            flleng=fwbs_variables.len_fw_channel,
            no90=no90fw,
            no180=no180fw,
            coolant_density=fwbs_variables.den_fw_coolant,
            coolant_dynamic_viscosity=fwbs_variables.visc_fw_coolant,
            coolant_electrical_conductivity=0.0e0,
            pol_channel_length=pollengi,
            nopolchan=npoltoti,
            label="Inboard first wall",
        )

        deltap_fwo = self.deltap_tot(
            output,
            icoolpump=1,
            flow_velocity=velfwo,
            flleng=fwbs_variables.len_fw_channel,
            no90=no90fw,
            no180=no180fw,
            coolant_density=fwbs_variables.den_fw_coolant,
            coolant_dynamic_viscosity=fwbs_variables.visc_fw_coolant,
            coolant_electrical_conductivity=0.0e0,
            pol_channel_length=pollengo,
            nopolchan=npoltoto,
            label="Outboard first wall",
        )

        # BB Presure Drops ###############

        # Long polodal flows
        if fwbs_variables.i_blkt_inboard == 1:
            npoltoti = fwbs_variables.nopol * npblkti_liq
        npoltoto = fwbs_variables.nopol * npblkto_liq

        deltap_blo = self.deltap_tot(
            output,
            icoolpump=1,
            flow_velocity=blanket_library.velblkto,
            flleng=blanket_library.bzfllengo,
            no90=no90bz,
            no180=no180bz,
            coolant_density=fwbs_variables.den_blkt_coolant,
            coolant_dynamic_viscosity=fwbs_variables.visc_blkt_coolant,
            coolant_electrical_conductivity=0.0e0,
            pol_channel_length=pollengo,
            nopolchan=npoltoto,
            label="Outboard blanket",
        )

        if fwbs_variables.i_blkt_inboard == 1:
            deltap_bli = self.deltap_tot(
                output,
                icoolpump=1,
                flow_velocity=blanket_library.velblkti,
                flleng=blanket_library.bzfllengi,
                no90=no90bz,
                no180=no180bz,
                coolant_density=fwbs_variables.den_blkt_coolant,
                coolant_dynamic_viscosity=fwbs_variables.visc_blkt_coolant,
                coolant_electrical_conductivity=0.0e0,
                pol_channel_length=pollengi,
                nopolchan=npoltoti,
                label="Inboard blanket",
            )

        # If the blanket has a liquid metal breeder...
        if fwbs_variables.icooldual > 0:
            deltap_blo_liq = self.deltap_tot(
                output,
                icoolpump=2,
                flow_velocity=velblkto_liq,
                flleng=bzfllengo_liq,
                no90=no90bz_liq,
                no180=no180bz_liq,
                coolant_density=fwbs_variables.den_liq,
                coolant_dynamic_viscosity=fwbs_variables.dynamic_viscosity_liq,
                coolant_electrical_conductivity=fwbs_variables.electrical_conductivity_liq,
                pol_channel_length=pollengo,
                nopolchan=npoltoto,
                label="Outboard blanket breeder liquid",
            )
            if fwbs_variables.i_blkt_inboard == 1:
                deltap_bli_liq = self.deltap_tot(
                    output,
                    icoolpump=2,
                    flow_velocity=velblkti_liq,
                    flleng=bzfllengi_liq,
                    no90=no90bz_liq,
                    no180=no180bz_liq,
                    coolant_density=fwbs_variables.den_liq,
                    coolant_dynamic_viscosity=fwbs_variables.dynamic_viscosity_liq,
                    coolant_electrical_conductivity=fwbs_variables.electrical_conductivity_liq,
                    pol_channel_length=pollengi,
                    nopolchan=npoltoti,
                    label="Inboard blanket breeder liquid",
                )

                return [
                    deltap_fwi,
                    deltap_fwo,
                    deltap_blo,
                    deltap_bli,
                    deltap_blo_liq,
                    deltap_bli_liq,
                ]
            return [deltap_fwi, deltap_fwo, deltap_blo, deltap_blo_liq]

        if fwbs_variables.i_blkt_inboard == 1:
            return [deltap_fwi, deltap_fwo, deltap_blo, deltap_bli]
        return [deltap_fwi, deltap_fwo, deltap_blo]

    def blanket_mod_pol_height(self):
        """Calculations for blanket module poloidal height
        author: J. Morris, CCFE, Culham Science Centre
        Calculations for blanket module poloidal height for D shaped and elliptical machines
        """
        if (
            physics_variables.itart == 1 or fwbs_variables.i_fw_blkt_vv_shape == 1
        ):  # D-shaped machine
            # Segment vertical inboard surface (m)
            blanket_library.bllengi = (
                2.0 * blanket_library.hblnkt
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
            b = blanket_library.hblnkt

            # Calculate ellipse circumference using Ramanujan approximation (m)
            ptor = np.pi * (3.0 * (a + b) - np.sqrt((3.0 * a + b) * (a + 3.0 * b)))

            # Calculate blanket poloidal length and segment, subtracting divertor length (m)
            # kit hcll version only had the single null option
            if physics_variables.idivrt == 2:
                # Double null configuration
                blanket_library.bllengo = (
                    0.5
                    * ptor
                    * (1.0 - 2.0 * fwbs_variables.fdiv)
                    / fwbs_variables.n_blkt_outboard_modules_poloidal
                )
            else:
                # single null configuration
                blanket_library.bllengo = (
                    0.5
                    * ptor
                    * (1.0 - fwbs_variables.fdiv)
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
            b = blanket_library.hblnkt

            # Distance between r1 and nearest edge of inboard first wall / blanket (m)
            a = r1 - (
                physics_variables.rmajor
                - physics_variables.rminor
                - build_variables.dr_fw_plasma_gap_inboard
            )

            # Calculate ellipse circumference using Ramanujan approximation (m)
            ptor = np.pi * (3.0 * (a + b) - np.sqrt((3.0 * a + b) * (a + 3.0 * b)))

            # Calculate inboard blanket poloidal length and segment, subtracting divertor length (m)
            # Assume divertor lies between the two ellipses, so fraction fdiv still applies

            # kit hcll version only had the single null option
            if physics_variables.idivrt == 2:
                # Double null configuration
                blanket_library.bllengi = (
                    0.5
                    * ptor
                    * (1.0 - 2.0 * fwbs_variables.fdiv)
                    / fwbs_variables.n_blkt_inboard_modules_poloidal
                )
            else:
                # single null configuration
                blanket_library.bllengi = (
                    0.5
                    * ptor
                    * (1.0 - fwbs_variables.fdiv)
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
            if physics_variables.idivrt == 2:
                # Double null configuration
                blanket_library.bllengo = (
                    0.5
                    * ptor
                    * (1.0 - 2.0 * fwbs_variables.fdiv)
                    / fwbs_variables.n_blkt_outboard_modules_poloidal
                )
            else:
                # single null configuration
                blanket_library.bllengo = (
                    0.5
                    * ptor
                    * (1.0 - fwbs_variables.fdiv)
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
        if fwbs_variables.i_bb_liq == 0:
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
        elif fwbs_variables.i_bb_liq == 1:
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
                physics_variables.bt
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
                physics_variables.bt
                * physics_variables.rmajor
                / (
                    physics_variables.rmajor
                    - (physics_variables.rmajor / physics_variables.aspect)
                )
            )
        # OB
        fwbs_variables.b_mag_blkt[1] = (
            physics_variables.bt
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
            fwbs_variables.b_mag_blkt
            * fwbs_variables.a_bz_liq
            / 2.0
            * np.sqrt(con_vsc_rat)
        )

        # Error for temperature range of breeder property realtions
        if fwbs_variables.i_bb_liq == 0 and (
            (t_ranges[:, 0] > mid_temp_liq).any()
            or (t_ranges[:, 1] < mid_temp_liq).any()
        ):
            error_handling.report_error(280)

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

        if fwbs_variables.icooldual == 1:
            po.ocmmnt(
                self.outfile,
                "Single coolant: liquid metal circulted for tritium extraction.",
            )
        if fwbs_variables.icooldual == 2:
            po.ocmmnt(self.outfile, "Dual coolant: self-cooled liquid metal breeder.")

        if fwbs_variables.i_bb_liq == 0:
            po.ocmmnt(
                self.outfile, "Blanket breeder type (i_bb_liq=0), PbLi (~ 17% Li)"
            )
        if fwbs_variables.i_bb_liq == 1:
            po.ocmmnt(self.outfile, "Blanket breeder type (i_bb_liq=1), Li")

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

        raise ValueError(f"i_channel_shape ={i_channel_shape} is an invalid option.")

    def thermo_hydraulic_model(self, output: bool):
        """
        Thermo-hydraulic model for first wall and blanket
        ONLY CALLED if i_coolant_pumping = 2 or 3

        Calculations for detailed powerflow model secondary_cycle > 1

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
            roughness epsilon           roughness
            peak FW temp (K)            temp_fw_peak
            maximum temp (K)            temp_fw_max
            FCI switch                  ---                     ---                 ifci

            Coolant                     FW                      BB primary          BB secondary

            primary coolant switch      i_fw_coolant_type               i_blkt_coolant_type              ---
            secondary coolant switch    ---                     ---                 i_bb_liq
            inlet temp (K)              temp_fw_coolant_in                 temp_blkt_coolant_in          inlet_temp_liq
            outlet temp (K)             temp_fw_coolant_out                temp_blkt_coolant_out         outlet_temp_liq
            pressure (Pa)               pres_fw_coolant              pres_blkt_coolant          blpressure_liq
        """
        ######################################################
        # Pre calculations needed for thermo-hydraulic model #
        ######################################################
        # IB/OB FW (MW)
        blanket_library.pnucfwi = (
            fwbs_variables.p_fw_nuclear_heat_total_mw
            * build_variables.a_fw_inboard
            / build_variables.a_fw_total
        )
        blanket_library.pnucfwo = (
            fwbs_variables.p_fw_nuclear_heat_total_mw
            * build_variables.a_fw_outboard
            / build_variables.a_fw_total
        )

        # IB/OB Blanket (MW)

        # Neutron power deposited in inboard blanket (MW)
        if fwbs_variables.i_blkt_inboard == 1:
            blanket_library.pnucblkti = (
                fwbs_variables.p_blkt_nuclear_heat_total_mw
                * fwbs_variables.vol_blkt_inboard
                / fwbs_variables.vol_blkt_total
            )

        # Neutron power deposited in outboard blanket (MW)
        blanket_library.pnucblkto = (
            fwbs_variables.p_blkt_nuclear_heat_total_mw
            * fwbs_variables.vol_blkt_outboard
            / fwbs_variables.vol_blkt_total
        )

        # For a dual-coolant blanket, some fraction of the power goes into the
        # structure of the BZ and is cooled by the primary coolant, and some fraction
        # goes into the liquid breeder to be cooled by itself.

        # If the blanket is dual-coolant...
        if fwbs_variables.icooldual == 2:
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
            f2py_compatible_to_string(fwbs_variables.i_fw_coolant_type).title()
            == "Helium"
            and fwbs_variables.i_blkt_coolant_type == 2
        ):
            fwbs_variables.i_fw_blkt_shared_coolant = 1
        if (
            f2py_compatible_to_string(fwbs_variables.i_fw_coolant_type).title()
            == "Water"
            and fwbs_variables.i_blkt_coolant_type == 1
        ):
            fwbs_variables.i_fw_blkt_shared_coolant = 1

        # If FW and BB have the same coolant...
        if fwbs_variables.i_fw_blkt_shared_coolant == 0:
            # Fraction of heat to be removed by IB/OB FW
            if fwbs_variables.icooldual == 2:
                f_nuc_fwi = (blanket_library.pnucfwi + fwbs_variables.psurffwi) / (
                    blanket_library.pnucfwi + fwbs_variables.psurffwi + pnucblkti_struct
                )
                f_nuc_fwo = (blanket_library.pnucfwo + fwbs_variables.psurffwo) / (
                    blanket_library.pnucfwo + fwbs_variables.psurffwo + pnucblkto_struct
                )
            else:
                f_nuc_fwi = (blanket_library.pnucfwi + fwbs_variables.psurffwi) / (
                    blanket_library.pnucfwi
                    + fwbs_variables.psurffwi
                    + blanket_library.pnucblkti
                )
                f_nuc_fwo = (blanket_library.pnucfwo + fwbs_variables.psurffwo) / (
                    blanket_library.pnucfwo
                    + fwbs_variables.psurffwo
                    + blanket_library.pnucblkto
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
        (blanket_library.tpeakfwi, _, _, blanket_library.mffwpi) = self.fw.fw_temp(
            output,
            fwbs_variables.radius_fw_channel,
            build_variables.dr_fw_inboard,
            build_variables.a_fw_inboard,
            fwbs_variables.psurffwi,
            blanket_library.pnucfwi,
            "Inboard first wall",
        )
        # (
        #     blanket_library.tpeakfwi,
        #     cf,
        #     rhof,
        #     blanket_library.mffwpi,
        # ) = fw_module.fw_temp(
        #     int(output),
        #     self.outfile,
        #     fwbs_variables.radius_fw_channel,
        #     build_variables.dr_fw_inboard,
        #     build_variables.a_fw_inboard,
        #     fwbs_variables.psurffwi,
        #     blanket_library.pnucfwi,
        #     "Inboard first wall",
        # )
        (fwbs_variables.tpeakfwo, cf, rhof, fwbs_variables.mffwpo) = self.fw.fw_temp(
            output,
            fwbs_variables.radius_fw_channel,
            build_variables.dr_fw_outboard,
            build_variables.a_fw_outboard,
            fwbs_variables.psurffwo,
            blanket_library.pnucfwo,
            "Outboard first wall",
        )
        # (fwbs_variables.tpeakfwo, cf, rhof, fwbs_variables.mffwpo) = fw_module.fw_temp(
        #     int(output),
        #     self.outfile,
        #     fwbs_variables.radius_fw_channel,
        #     build_variables.dr_fw_outboard,
        #     build_variables.a_fw_outboard,
        #     fwbs_variables.psurffwo,
        #     blanket_library.pnucfwo,
        #     "Outboard first wall",
        # )

        # Peak first wall temperature (K)
        fwbs_variables.temp_fw_peak = max(
            blanket_library.tpeakfwi, blanket_library.tpeakfwo
        )

        # Total mass flow rate to remove inboard FW power (kg/s)
        blanket_library.mffwi = (
            1.0e6
            * (blanket_library.pnucfwi + fwbs_variables.psurffwi)
            / (fwbs_variables.cp_fw * (fwoutleti - fwbs_variables.temp_fw_coolant_in))
        )
        # Total mass flow rate to remove outboard FW power (kg/s)
        blanket_library.mffwo = (
            1.0e6
            * (blanket_library.pnucfwo + fwbs_variables.psurffwo)
            / (fwbs_variables.cp_fw * (fwoutleto - fwbs_variables.temp_fw_coolant_in))
        )

        # If the blanket is dual-coolant...
        if fwbs_variables.icooldual == 2:
            # Mass flow rates for outboard blanket coolants (kg/s)
            blanket_library.mfblkto = (
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
                blanket_library.mfblkti = (
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
        elif fwbs_variables.icooldual == 1:
            # Mass flow rate for outboard blanket coolant (kg/s)
            blanket_library.mfblkto = (
                1.0e6
                * (blanket_library.pnucblkto)
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
                blanket_library.mfblkti = (
                    1.0e6
                    * (blanket_library.pnucblkti)
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
            blanket_library.mfblkto = (
                1.0e6
                * (blanket_library.pnucblkto)
                / (
                    fwbs_variables.cp_bl
                    * (fwbs_variables.temp_blkt_coolant_out - inlet_tempo)
                )
            )

            # If there is an IB blanket...
            # Mass flow rate for inboard blanket coolant (kg/s)
            if fwbs_variables.i_blkt_inboard == 1:
                blanket_library.mfblkti = (
                    1.0e6
                    * (blanket_library.pnucblkti)
                    / (
                        fwbs_variables.cp_bl
                        * (fwbs_variables.temp_blkt_coolant_out - inlet_tempi)
                    )
                )

        ########################################################
        # Handling of pressure drops and coolant pumping power #
        ########################################################

        # load in pressures if primary pumping == 2
        if fwbs_variables.i_coolant_pumping == 2:
            deltap = self.thermo_hydraulic_model_pressure_drop_calculations(
                output=output
            )
            deltap_fwi = deltap[0]
            deltap_fwo = deltap[1]
            deltap_blo = deltap[2]
            if fwbs_variables.icooldual > 0:
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
            if fwbs_variables.i_coolant_pumping == 2:
                if fwbs_variables.i_blkt_inboard == 1:
                    deltap_fw_blkt = deltap_fwi + deltap_bli + deltap_fwo + deltap_blo
                if fwbs_variables.i_blkt_inboard == 0:
                    deltap_fw_blkt = deltap_fwi + deltap_fwo + deltap_blo
            elif fwbs_variables.i_coolant_pumping == 3:
                deltap_fw_blkt = primary_pumping_variables.dp_fw_blkt
            # Total coolant mass flow rate in the first wall/blanket (kg/s)
            blanket_library.mftotal = blanket_library.mffwi + blanket_library.mffwo

            # Total mechanical pumping power (MW)
            primary_pumping_variables.htpmw_fw_blkt = self.pumppower(
                output=output,
                icoolpump=1,
                temp_in=fwbs_variables.temp_fw_coolant_in.item(),
                temp_out=fwbs_variables.temp_blkt_coolant_out.item(),
                pressure=fwbs_variables.pres_fw_coolant.item(),
                pdrop=deltap_fw_blkt,
                mf=blanket_library.mftotal,
                primary_coolant_switch=f2py_compatible_to_string(
                    fwbs_variables.i_fw_coolant_type
                ),
                coolant_density=fwbs_variables.den_fw_coolant,
                label="First Wall and Blanket",
            )

        # If FW and BB have different coolants...
        elif fwbs_variables.i_fw_blkt_shared_coolant == 1:
            if fwbs_variables.i_coolant_pumping == 2:
                # Total pressure drop in the first wall (Pa)
                deltap_fw = deltap_fwi + deltap_fwo

                # Total pressure drop in the blanket (Pa)
                if fwbs_variables.i_blkt_inboard == 1:
                    deltap_blkt = deltap_bli + deltap_blo
                if fwbs_variables.i_blkt_inboard == 0:
                    deltap_blkt = deltap_blo
            elif fwbs_variables.i_coolant_pumping == 3:
                deltap_fw = primary_pumping_variables.dp_fw
                deltap_blkt = primary_pumping_variables.dp_blkt

            # Total coolant mass flow rate in the first wall (kg/s)
            blanket_library.mffw = blanket_library.mffwi + blanket_library.mffwo
            # Total coolant mass flow rate in the blanket (kg/s)
            blanket_library.mfblkt = blanket_library.mfblkti + blanket_library.mfblkto

            # Mechanical pumping power for the first wall (MW)
            heat_transport_variables.htpmw_fw = self.pumppower(
                output=output,
                icoolpump=1,
                temp_in=fwbs_variables.temp_fw_coolant_in.item(),
                temp_out=fwbs_variables.temp_fw_coolant_out.item(),
                pressure=fwbs_variables.pres_fw_coolant.item(),
                pdrop=deltap_fw.item(),
                mf=blanket_library.mffw,
                primary_coolant_switch=f2py_compatible_to_string(
                    fwbs_variables.i_fw_coolant_type
                ),
                coolant_density=fwbs_variables.den_fw_coolant,
                label="First Wall",
            )

            # Mechanical pumping power for the blanket (MW)
            heat_transport_variables.htpmw_blkt = self.pumppower(
                output=output,
                icoolpump=1,
                temp_in=fwbs_variables.temp_blkt_coolant_in.item(),
                temp_out=fwbs_variables.temp_blkt_coolant_out.item(),
                pressure=fwbs_variables.pres_blkt_coolant.item(),
                pdrop=deltap_blkt.item(),
                mf=blanket_library.mfblkt,
                primary_coolant_switch=(
                    "Helium" if fwbs_variables.i_blkt_coolant_type == 1 else "Water"
                ),
                coolant_density=blanket_library.den_blkt_coolant,
                label="Blanket",
            )

            # Total mechanical pumping power (MW)
            primary_pumping_variables.htpmw_fw_blkt = (
                heat_transport_variables.htpmw_fw + heat_transport_variables.htpmw_blkt
            )

        # If the blanket has a liquid metal breeder...
        if fwbs_variables.icooldual > 0:
            # Total pressure drop in the blanket (Pa)
            if fwbs_variables.i_coolant_pumping == 2:
                if fwbs_variables.i_blkt_inboard == 1:
                    deltap_bl_liq = deltap_bli_liq + deltap_blo_liq
                if fwbs_variables.i_blkt_inboard == 0:
                    deltap_bl_liq = deltap_blo_liq
            elif fwbs_variables.i_coolant_pumping == 3:
                deltap_bl_liq = primary_pumping_variables.dp_liq
            # Total liquid metal breeder/coolant mass flow rate in the blanket (kg/s)
            fwbs_variables.mfblkt_liq = (
                blanket_library.mfblkti_liq + blanket_library.mfblkto_liq
            )

            # Mechanical pumping power for the blanket (MW)
            heat_transport_variables.htpmw_blkt_liq = self.pumppower(
                output=output,
                icoolpump=2,
                temp_in=fwbs_variables.inlet_temp_liq.item(),
                temp_out=fwbs_variables.outlet_temp_liq.item(),
                pressure=fwbs_variables.blpressure_liq.item(),
                pdrop=deltap_bl_liq,
                mf=fwbs_variables.mfblkt_liq,
                primary_coolant_switch=(
                    "Helium" if fwbs_variables.i_blkt_coolant_type == 1 else "Water"
                ),
                coolant_density=fwbs_variables.den_liq,
                label="Liquid Metal Breeder/Coolant",
            )

            heat_transport_variables.htpmw_blkt_tot = (
                primary_pumping_variables.htpmw_fw_blkt
                + heat_transport_variables.htpmw_blkt_liq
            )

        if output:
            po.oheadr(
                self.outfile, "Summary of first wall and blanket thermohydraulics"
            )

            # FW
            po.osubhd(self.outfile, "First wall: ")

            po.ovarst(
                self.outfile,
                "First wall coolant type",
                "(i_fw_coolant_type)",
                f'"{fwbs_variables.i_fw_coolant_type}"',
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
                "Roughness of first wall cooling channels (m)",
                "(roughness)",
                fwbs_variables.roughness,
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
                    "(mffw)",
                    fwbs_variables.mffw,
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
                    "(mfblkt)",
                    fwbs_variables.mfblkt,
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
            if fwbs_variables.icooldual > 0:
                po.osubhd(self.outfile, "Breeding Blanket (breeder): ")

                po.ovarin(
                    self.outfile,
                    "Blanket liquid breeder type (0=PbLi, 1=Li)",
                    "(i_bb_liq)",
                    fwbs_variables.i_bb_liq,
                )
                if fwbs_variables.icooldual == 2:
                    po.ocmmnt(
                        self.outfile, "Dual-coolant BB, i.e. self-cooled breeder."
                    )
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
                    fwbs_variables.mfblkt_liq,
                    "OP ",
                )

            # Pumping Power
            po.osubhd(self.outfile, "Mechanical pumping power: ")

            if fwbs_variables.i_fw_blkt_shared_coolant == 1:
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for FW (MW)",
                    "(htpmw_fw)",
                    fwbs_variables.htpmw_fw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for blanket (primary) coolant (MW)",
                    "(htpmw_blkt)",
                    fwbs_variables.htpmw_blkt,
                    "OP ",
                )
            if fwbs_variables.icooldual > 0:
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for blanket liquid breeder (MW)",
                    "(htpmw_blkt_liq)",
                    heat_transport_variables.htpmw_blkt_liq,
                    "OP ",
                )
            po.ovarre(
                self.outfile,
                "Total mechanical pumping power for FW and blanket (MW)",
                "(htpmw_fw_blkt)",
                primary_pumping_variables.htpmw_fw_blkt,
                "OP ",
            )
            if fwbs_variables.icooldual > 0:
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
                "(htpmw_div)",
                heat_transport_variables.htpmw_div,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Pumping power for shield and vacuum vessel (MW)",
                "(htpmw_shld)",
                heat_transport_variables.htpmw_shld,
                "OP ",
            )

    def deltap_tot(
        self,
        output: bool,
        icoolpump,
        flow_velocity,
        flleng,
        no90,
        no180,
        coolant_density,
        coolant_dynamic_viscosity,
        coolant_electrical_conductivity,
        pol_channel_length,
        nopolchan,
        label,
    ):
        """Routine to calculate the coolant pumping power in MW in the FW and BZ.
        Adapted from previous pumppower function.

        original author: P. J. Knight, CCFE
        original references: Idel'Cik, I. E. (1969), Memento des pertes de charge;
        A Textbook on Heat Transfer, S.P. Sukhatme, 2005

        author: G. Graham

        :param icoolpump: Switch for primary coolant or secondary coolant/breeder (1=primary He/H2O, 2=secondary PbLi/Li)
        :param flow_velocity: Coolant flow velocity (m/s)
        :param flleng: Total flow length along pipe (m)
        :param no90: Number of 90 degree bends in pipe
        :param no180: Number of 180 degree bends in pipe
        """
        # Friction - for all coolants
        frict_drop = self.pressure_drop(
            icoolpump,
            no90,
            no180,
            flleng,
            coolant_density,
            coolant_dynamic_viscosity,
            flow_velocity,
            label,
            output=output,
        )

        if icoolpump == 2:
            mhd_drop = self.liquid_breeder_pressure_drop_mhd(
                flow_velocity,
                coolant_dynamic_viscosity,
                coolant_electrical_conductivity,
                pol_channel_length,
                nopolchan,
                label,
                output=output,
            )
        else:
            mhd_drop = 0

        # Total pressure drop (Pa)
        deltap_tot = frict_drop + mhd_drop

        if output:
            po.osubhd(self.outfile, f"Total pressure drop for {label}")

            po.ocmmnt(self.outfile, "Friction drops plus MHD drops if applicaple")
            po.ovarre(
                self.outfile, "Total pressure drop (Pa)", "(deltap)", deltap_tot, "OP "
            )
            po.ovarre(
                self.outfile,
                "Coolant flow velocity (m/s)",
                "(flow_velocity, formerly vv)",
                flow_velocity,
                "OP ",
            )

        return deltap_tot

    def liquid_breeder_pressure_drop_mhd(
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
        if fwbs_variables.ifci != 1:
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

            if fwbs_variables.ifci == 0:
                po.ocmmnt(
                    self.outfile, "Flow channels have thin conducting walls (ifci==0)"
                )
                po.ovarre(
                    self.outfile,
                    "Wall conductance (A V-1 m-1)",
                    "(bz_channel_conduct_liq)",
                    fwbs_variables.bz_channel_conduct_liq,
                    "OP ",
                )
            elif fwbs_variables.ifci == 2:
                po.ocmmnt(self.outfile, "Flow Channel Inserts (FCIs) used (ifci==2)")
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
                    "Flow Channel Inserts - assumed perfect insulator (ifci==1)",
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

    def pressure_drop(
        self,
        i_ps: int,
        num_90: float,
        num_180: float,
        l_pipe: float,
        den: float,
        vsc: float,
        vv: float,
        label: str,
        output: bool = False,
    ):
        """Pressure drops are calculated for a pipe with a number of 90
        and 180 degree bends. The pressure drop due to frictional forces along
        the total straight length of the pipe is calculated, then the pressure
        drop due to the bends is calculated. The total pressure drop is the sum
        of all contributions.

        original author: P. J. Knight, CCFE
        moved from previous version of pumppower function by: G Graham, CCFE

        :param i_ps: switch for primary or secondary coolant
        :param num_90: number of 90 degree bends in the pipe
        :param num_180: number of 180 degree bends in the pipe
        :param l_pipe: total flow length along pipe (m)
        :param den: coolant density (kg/m3)
        :param vsc: coolant viscosity (Pa s)
        :param vv: coolant flow velocity (m/s)
        :param label: component name
        :param output: boolean of whether to write data to output file

        N.B Darcy-Weisbach Equation (straight pipe):

         kstrght = lambda * L/D

         pressure drop = kstrght * (rho*V^2)/2

         lambda - Darcy friction factor, L - pipe length, D - hydraulic diameter,
         rho - fluid density, V - fluid flow average velocity

        This function also calculates pressure drop equations for elbow bends,
        with modified coefficients.

        N.B. Darcy friction factor is estimated from the Haaland approximation.
        """

        # Calculate hydraulic dimater for round or retancular pipe (m)
        dh = self.hydraulic_diameter(i_ps)

        # Reynolds number
        reyn = den * vv * dh / vsc

        # Calculate Darcy friction factor
        # N.B. friction function Uses Haaland approx. which assumes a filled circular pipe.
        # Use dh which allows us to do fluid calculations for non-cicular tubes
        # (dh is estimate appropriate for fully developed flow).
        lamda = self.fw.friction(reyn)

        # Pressure drop coefficient

        # Straight section
        kstrght = lamda * l_pipe / dh

        # In preveious version of pumppower:
        # - elbow radius assumed = 0.018m for 90 degree elbow, from WCLL
        # - elbow radius assumed half that of 90 deg case for 180 deg elbow
        # Intialised value for radius_fw_channel is 0.006m, so elbow radius = 3 * radius_fw_channel,
        # aka 1.5 * pipe diameter, which seems to be engineering standard for
        # a steel pipe long-radius elbow (short-radius elbow = 2 * radius_fw_channel).

        # If primary coolant or secondary coolant (See DCLL)
        elbow_radius = (
            (3 * fwbs_variables.radius_fw_channel)
            if (i_ps == 1)
            else fwbs_variables.b_bz_liq
        )

        # 90 degree elbow pressure drop coefficient
        kelbwn = self.elbow_coeff(elbow_radius, 90.0, lamda, dh)

        # 180 degree elbow pressure drop coefficient
        kelbwt = self.elbow_coeff(elbow_radius / 2, 180.0, lamda, dh)

        # Total (Pa)
        pdropstraight = kstrght * 0.5 * den * vv * vv
        pdrop90 = num_90 * kelbwn * 0.5 * den * vv * vv
        pdrop180 = num_180 * kelbwt * 0.5 * den * vv * vv

        pressure_drop = pdropstraight + pdrop90 + pdrop180

        if output:
            po.osubhd(self.outfile, f"Pressure drop (friction) for {label}")
            po.ovarre(self.outfile, "Reynolds number", "(reyn)", reyn, "OP ")
            po.ovarre(self.outfile, "Darcy friction factor", "(lambda)", lamda, "OP ")
            po.ovarre(
                self.outfile,
                "Pressure drop (Pa)",
                "(pressure_drop)",
                pressure_drop,
                "OP ",
            )
            po.ocmmnt(self.outfile, "This is the sum of the following:")
            po.ovarre(
                self.outfile,
                "            Straight sections (Pa)",
                "(pdropstraight)",
                pdropstraight,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "            90 degree bends (Pa)",
                "(pdrop90)",
                pdrop90,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "            180 degree bends (Pa)",
                "(pdrop180)",
                pdrop180,
                "OP ",
            )

            # TN: always write verbose stuff, it has no harm
            po.ovarre(
                self.outfile,
                "Straight section pressure drop coefficient",
                "(kstrght)",
                kstrght,
                "OP ",
            )
            po.ovarre(
                self.outfile, "90 degree elbow coefficient", "(kelbwn)", kelbwn, "OP "
            )
            po.ovarre(
                self.outfile,
                "180 degree elbow coefficient coefficient",
                "(kelbwt)",
                kelbwt,
                "OP ",
            )

        return pressure_drop

    def hydraulic_diameter(self, i_channel_shape):
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

        raise ValueError(f"i_channel_shape ={i_channel_shape} is an invalid option.")

    def elbow_coeff(self, r_elbow, ang_elbow, lamda, dh):
        """Function calculates elbow bends coefficients for pressure drop
        calculations.

        author: G. Graham, CCFE

        :param r_elbow: pipe elbow radius (m)
        :param ang_elbow: pipe elbow angle (degrees)
        :param lamda: darcy Friction Factor
        :param dh: hydraulic diameter (m)

        References:

             [Ide1969]   Idel'Cik, I. E. (1969), Memento des pertes de charge,
                         Collection de la Direction des Etudes et Recherches d'Electricité de France.
        """

        if ang_elbow == 90:
            a = 1.0
        elif ang_elbow < 70:
            a = 0.9 * np.sin(ang_elbow * np.pi / 180.0)
        elif ang_elbow > 100:
            a = 0.7 + (0.35 * np.sin((ang_elbow / 90.0) * (np.pi / 180.0)))
        else:
            raise ValueError(
                "No formula for 70 <= elbow angle(deg) <= 100, only 90 deg option available in this range."
            )

        r_ratio = r_elbow / dh
        if r_ratio > 1:
            b = 0.21 / r_ratio**0.5
        elif r_ratio < 1:
            b = 0.21 / r_ratio**2.5
        else:
            b = 0.21

        # Singularity
        ximt = a * b

        # Friction
        xift = 0.0175 * lamda * (r_elbow / dh) * ang_elbow

        # Elbow Coefficient
        return ximt + xift

    def pumppower(
        self,
        output,
        icoolpump,
        temp_in,
        temp_out,
        pressure,
        pdrop,
        mf,
        primary_coolant_switch,
        coolant_density,
        label,
    ):
        """Routine to calculate the coolant pumping power in MW in the FW and BZ.
        Adapted from previous pumppower function.

        original author: P. J. Knight, CCFE
        original references: Idel'Cik, I. E. (1969), Memento des pertes de charge;
        A Textbook on Heat Transfer, S.P. Sukhatme, 2005

        author: G. Graham

        :param icoolpump: Switch for primary coolant or secondary coolant/breeder (1=primary He/H2O, 2=secondary PbLi/Li)
        :param temp_in: Inlet (pump oulet) temperature (K)
        :param temp_out: Oulet (pump inlet) temperature (K)
        :param pressure: Outlet (pump inlet) coolant pressure (Pa)
        :param pdrop: Pressure drop (Pa)
        :param mf: Total coolant mass flow rate in (kg/s)
        :param primary_coolant_switch: Switch for FW/blanket coolant, (1=He or 2=H2O) if icoolpump=1
        :param coolant_density: Density of coolant or liquid breeder
        """
        # Pumping power !

        # Outlet pressure is 'pressure'
        # Inlet pressure (Pa)
        coolpin = pressure + pdrop

        # Adiabatic index for helium or water
        gamma = (5 / 3) if fwbs_variables.i_blkt_coolant_type == 1 else (4 / 3)

        # If caculating for primary coolant...
        if icoolpump == 1:
            # Comments from original pumppower function:
            # The pumping power is be calculated in the most general way,
            # using enthalpies before and after the pump.

            fluid_properties = FluidProperties.of(
                fluid_name=primary_coolant_switch,
                temperature=temp_in,
                pressure=coolpin,
            )

            # Assume isentropic pump so that s1 = s2
            s1 = fluid_properties.entropy

            # Get specific enthalpy at the outlet (J/kg) before pump using pressure and entropy s1
            outlet_fluid_properties = FluidProperties.of(
                fluid_name=primary_coolant_switch,
                pressure=pressure,
                entropy=s1,
            )

            # Pumping power (MW) is given by enthalpy change, with a correction for
            # the isentropic efficiency of the pump.
            fp = (
                temp_in
                * (1 - (coolpin / pressure) ** -((gamma - 1) / gamma))
                / (fwbs_variables.etaiso * (temp_out - temp_in))
            )
            pumppower = (
                1e-6
                * mf
                * (fluid_properties.enthalpy - outlet_fluid_properties.enthalpy)
                / fwbs_variables.etaiso
            ) / (1 - fp)

        # If calculating for secondary coolant/breeder...
        else:
            # Calculate specific volume
            spec_vol = 1 / coolant_density

            # Pumping power (MW) is given by pressure change, with a correction for
            # the isentropic efficiency of the pump.
            fp = (
                temp_in
                * (1 - (coolpin / pressure) ** -((gamma - 1) / gamma))
                / (fwbs_variables.etaiso_liq * (temp_out - temp_in))
            )
            pumppower = (1e-6 * mf * spec_vol * pdrop / fwbs_variables.etaiso_liq) / (
                1 - fp
            )

        # Error for pdrop too large
        if fp >= 1:
            eh.report_error(279)

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
                coolpin,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "FW or Blanket oulet (pump inlet) pressure (Pa)",
                "(pressure)",
                pressure,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "FW or Blanket total pressure drop (Pa)",
                "(pdrop)",
                pdrop,
                "OP ",
            )
            po.ovarre(self.outfile, "Mass flow rate in (kg/s) = ", "(mf)", mf, "OP ")

        return pumppower


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


def dshellarea(rmajor, rminor, zminor):
    """Routine to calculate the inboard, outboard and total surface areas
    of a D-shaped toroidal shell
    author: P J Knight, CCFE, Culham Science Centre
    rmajor : input real : major radius of inboard straight section (m)
    rminor : input real : horizontal width of shell (m)
    zminor : input real : vertical half-height of shell (m)
    ain    : output real : surface area of inboard straight section (m3)
    aout   : output real : surface area of outboard curved section (m3)
    atot   : output real : total surface area of shell (m3)
    This routine calculates the surface area of the inboard and outboard
    sections of a D-shaped toroidal shell defined by the above input
    parameters.
    The inboard section is assumed to be a cylinder.
    The outboard section is defined by a semi-ellipse, centred on the
    major radius of the inboard section.
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
