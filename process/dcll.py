import process.blanket_library as blanket_library
from process import constants
from process import (
    process_output as po,
)
from process.blanket_library import InboardBlanket, OutboardBlanket
from process.data_structure import (
    build_variables,
    current_drive_variables,
    dcll_variables,
    divertor_variables,
    first_wall_variables,
    fwbs_variables,
    heat_transport_variables,
    physics_variables,
    primary_pumping_variables,
)


class DCLL(InboardBlanket, OutboardBlanket):
    """This module contains the Dual Coolant Lead Lithium (DCLL) specific submods of PROCESSS.



    Acronyms for this module:

         BB          Breeding Blanket
         FW          First Wall
         BZ          Breeder Zone
         MF/BSS      Manifold/Back Supporting Structure
         LT          Low Temperature
         HT          High Temperature
         MMS         Multi Module Segment
         SMS         Single Modle Segment
         IB          Inboard
         OB          Outboard
         HCD         Heating & Current Drive
         FCI         Flow Channel Insert

    IN.DAT info for DCLL:

         Select DCLL model
             i_blanket_type = 5 * DCLL

         Liquid Metal Breeder Material = PbLi
             i_blkt_liquid_breeder_type = 0 * Liquid Metal Breeder Material = PbLi

         Specify dual-coolant i.e., get mass flow required from heat extracted from liqid metal breeder
             i_blkt_dual_coolant = 2

         FIC switch: 0 = no FIC, Eurofer; 1 = FCIs, perfect electrical insulator, 2 = FCIs, with specified conductance
             i_blkt_liquid_breeder_channel_type = 0, 1, or 2

         Liquid metal duct wall conductance initilized at Eurofer value in fwbs_variables, or can input other value, used for i_blkt_liquid_breeder_channel_type = 0 or 2
             (bz_channel_conduct_liq)

         Choose if FW and BB structure are on the same pumping system (unless have diffent coolants), default is same coolant with flow IN->FW->BB->OUT
             (i_fw_blkt_shared_coolant)

         Can set inlet and oulet temperature for liquid metal breeder
             (inlet_temp_liq)
             (outlet_temp_liq)

    References:

         [Nat1995]   Natesan et al. (1995), Assessment of alkali metal coolants for
                     for the ITER blanket, Fusion Engineering and Design 27, 457-466

         [Mal1995]   Malang and Mattas (1995), Comparison of lithium and the eutectic
                     lead-lithium alloy, two candidate liquid metal breeder materials
                     for self-cooled blankets, Fusion Engineering and Design 27, 399-406

         [Gas2001]   Gasior and Mozer (2001), Thermodynamic study of liquid lithium-lead
                     alloys using the EMF method, Journal of Nuclear Materials, 294, 77-83

         [Pal2016]   Palermo et al. (2016), Neutronic analyses of the preliminary design
                     of a DCLL blanket for the EUROfusion DEMO power plant,
                     Fusion Engineering and Design 109-111.

         [Gar2017]   Garcinuno et al. (2017), Design of a permeator against vacuum for
                     tritium extraction from eutectic lithium-lead in a DCLL DEMO,
                     Fusion Engineering and Design, 117, 226-231

         [Fer2021]   Fernandez-Berceruelo et al. (2021), Alternatives for upgrading the
                     eu dcll breeding blanket from mms to sms, Fusion Engineering and
                     Design 167, 112380


    Note: request for when CCFE Bluemira nutronics work is added: output maximum values, as well as average values, for wall neutronics calculation if possible.
    """

    def run(self, output: bool):
        self.component_volumes()
        dia_blkt_channel = self.pipe_hydraulic_diameter(i_channel_shape=1)
        fwbs_variables.radius_blkt_channel = dia_blkt_channel / 2
        (
            fwbs_variables.radius_blkt_channel_90_bend,
            fwbs_variables.radius_blkt_channel_180_bend,
        ) = self.calculate_pipe_bend_radius(i_ps=1)

        self.set_blanket_module_geometry()

        blanket_library.len_blkt_inboard_segment_toroidal = self.calculate_blanket_inboard_module_geometry(
            n_blkt_inboard_modules_toroidal=fwbs_variables.n_blkt_inboard_modules_toroidal,
            rmajor=physics_variables.rmajor,
            rminor=physics_variables.rminor,
            dr_fw_plasma_gap_inboard=build_variables.dr_fw_plasma_gap_inboard,
        )
        blanket_library.len_blkt_outboard_segment_toroidal = self.calculate_blanket_outboard_module_geometry(
            n_blkt_outboard_modules_toroidal=fwbs_variables.n_blkt_outboard_modules_toroidal,
            rmajor=physics_variables.rmajor,
            rminor=physics_variables.rminor,
            dr_fw_plasma_gap_outboard=build_variables.dr_fw_plasma_gap_outboard,
        )

        self.primary_coolant_properties(output=output)
        self.liquid_breeder_properties(output=output)
        self.dcll_neutronics_and_power(output=output)
        self.dcll_masses(output=output)
        self.dcll_power_and_heating(output=output)

        if output:
            self.write_output()

    def dcll_neutronics_and_power(self, output: bool):
        """This is a tempory module that will use results from CCFE Bluemira nutronics work (once completed).
        Database will provide values for power deposition in FW & BB, BB TBR, and nuron fluence at TF coil for
        different thicknesses of BB and meterial fractions.

        For now we use the same method as KIT HCLL and the user can select approprite fractional
        values from DCLL nutronics studies as inputs.
        See fwbs_variables:
             - pnuc_fw_ratio_dcll
             - pnuc_blkt_ratio_dcll
             - f_nuc_pow_bz_struct
             - f_nuc_pow_bz_liq

        Parameters
        ----------
        output: bool

        """

        if divertor_variables.n_divertors == 2:
            # Double null configuration
            covf = (
                1
                - (2 * fwbs_variables.f_ster_div_single)
                - fwbs_variables.f_a_fw_outboard_hcd
            )
        else:
            # Single null configuration
            covf = (
                1 - fwbs_variables.f_ster_div_single - fwbs_variables.f_a_fw_outboard_hcd
            )

        # Nuclear heating in the first wall (MW)
        fwbs_variables.p_fw_nuclear_heat_total_mw = (
            physics_variables.p_neutron_total_mw
            * fwbs_variables.pnuc_fw_ratio_dcll
            * covf
        )

        # Nuclear heating in the blanket with energy multiplication (MW)
        fwbs_variables.pnuc_blkt_ratio_dcll = 1 - fwbs_variables.pnuc_fw_ratio_dcll
        fwbs_variables.p_blkt_nuclear_heat_total_mw = (
            physics_variables.p_neutron_total_mw
            * fwbs_variables.pnuc_blkt_ratio_dcll
            * fwbs_variables.f_p_blkt_multiplication
            * covf
        )

        # Energy multiplication energy (MW)
        fwbs_variables.p_blkt_multiplication_mw = (
            (physics_variables.p_neutron_total_mw * fwbs_variables.pnuc_blkt_ratio_dcll)
            * (fwbs_variables.f_p_blkt_multiplication - 1)
            * covf
        )

        # HCD Apperatus

        # No nuclear heating of the H & CD
        fwbs_variables.p_fw_hcd_nuclear_heat_mw = 0
        # Radiation power incident on HCD apparatus (MW)
        fwbs_variables.p_fw_hcd_rad_total_mw = (
            physics_variables.p_plasma_rad_mw * fwbs_variables.f_a_fw_outboard_hcd
        )

        # FW

        # Radiation power incident on first wall (MW)
        fwbs_variables.p_fw_rad_total_mw = (
            physics_variables.p_plasma_rad_mw
            - fwbs_variables.p_div_rad_total_mw
            - fwbs_variables.p_fw_hcd_rad_total_mw
        )

        # Surface heat flux on first wall (MW)
        # All of the fast particle losses go to the outer wall.
        fwbs_variables.psurffwo = (
            fwbs_variables.p_fw_rad_total_mw
            * first_wall_variables.a_fw_outboard
            / first_wall_variables.a_fw_total
            + current_drive_variables.p_beam_orbit_loss_mw
            + physics_variables.p_fw_alpha_mw
        )
        fwbs_variables.psurffwi = fwbs_variables.p_fw_rad_total_mw * (
            1 - first_wall_variables.a_fw_outboard / first_wall_variables.a_fw_total
        )

        if output:
            po.osubhd(
                self.outfile, "DCLL model: Nuclear and Radiation Heating of Components"
            )

            po.osubhd(self.outfile, "Component Coverage :")

            po.ovarre(
                self.outfile,
                "Solid angle fraction taken by on divertor",
                "(f_ster_div_single)",
                fwbs_variables.f_ster_div_single,
            )
            po.ovarre(
                self.outfile,
                "Fraction of outboard first wall area covered by HCD and diagnostics",
                "(f_a_fw_outboard_hcd)",
                fwbs_variables.f_a_fw_outboard_hcd,
            )
            po.ovarre(self.outfile, "Blanket coverage factor", "(covf)", covf)

            po.osubhd(self.outfile, "Nuclear heating :")

            po.ovarre(
                self.outfile,
                "Total nuclear heating in FW (MW)",
                "(p_fw_nuclear_heat_total_mw)",
                fwbs_variables.p_fw_nuclear_heat_total_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Energy multiplication in the blanket",
                "(f_p_blkt_multiplication)",
                fwbs_variables.f_p_blkt_multiplication,
                "OP ",
            )
            po.ocmmnt(
                self.outfile,
                "(Note: f_p_blkt_multiplication is fixed for this model inside the code)",
            )
            po.ovarre(
                self.outfile,
                "Total nuclear heating in the blanket (including f_p_blkt_multiplication) (MW)",
                "(p_blkt_nuclear_heat_total_mw)",
                fwbs_variables.p_blkt_nuclear_heat_total_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total nuclear heating in the shield (MW)",
                "(p_shld_nuclear_heat_mw)",
                fwbs_variables.p_shld_nuclear_heat_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total nuclear heating in the divertor (MW)",
                "(p_div_nuclear_heat_total_mw)",
                fwbs_variables.p_div_nuclear_heat_total_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total nuclear heating in TF+PF coils (CS is negligible) (MW)",
                "(p_tf_nuclear_heat_mw)",
                fwbs_variables.p_tf_nuclear_heat_mw,
                "OP ",
            )

            po.osubhd(self.outfile, "Radiation heating :")

            po.ovarrf(
                self.outfile,
                "Radiation heating power into the divertor (MW)",
                "(p_div_rad_total_mw)",
                fwbs_variables.p_div_rad_total_mw,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Radiation heating power into the first wall (MW)",
                "(p_fw_rad_total_mw)",
                fwbs_variables.p_fw_rad_total_mw,
                "OP ",
            )

    def dcll_power_and_heating(self, output: bool):
        # Mechanical Pumping

        # For i_p_coolant_pumping == 0:
        # User sets mechanical pumping power directly (primary_pumping_power)
        # Values of p_blkt_coolant_pump_mw, p_div_coolant_pump_mw, p_fw_coolant_pump_mw, p_shld_coolant_pump_mw set in input file

        if fwbs_variables.i_p_coolant_pumping == 1:
            # User sets mechanical pumping power directly
            (
                heat_transport_variables.p_fw_coolant_pump_mw,
                heat_transport_variables.p_blkt_coolant_pump_mw,
                heat_transport_variables.p_shld_coolant_pump_mw,
                heat_transport_variables.p_div_coolant_pump_mw,
            ) = blanket_library.set_pumping_powers_as_fractions(
                f_p_fw_coolant_pump_total_heat=heat_transport_variables.f_p_fw_coolant_pump_total_heat,
                f_p_blkt_coolant_pump_total_heat=heat_transport_variables.f_p_blkt_coolant_pump_total_heat,
                f_p_shld_coolant_pump_total_heat=heat_transport_variables.f_p_shld_coolant_pump_total_heat,
                f_p_div_coolant_pump_total_heat=heat_transport_variables.f_p_div_coolant_pump_total_heat,
                p_fw_nuclear_heat_total_mw=fwbs_variables.p_fw_nuclear_heat_total_mw,
                psurffwi=fwbs_variables.psurffwi,
                psurffwo=fwbs_variables.psurffwo,
                p_blkt_nuclear_heat_total_mw=fwbs_variables.p_blkt_nuclear_heat_total_mw,
                p_shld_nuclear_heat_mw=fwbs_variables.p_shld_nuclear_heat_mw,
                p_cp_shield_nuclear_heat_mw=fwbs_variables.p_cp_shield_nuclear_heat_mw,
                p_plasma_separatrix_mw=physics_variables.p_plasma_separatrix_mw,
                p_div_nuclear_heat_total_mw=fwbs_variables.p_div_nuclear_heat_total_mw,
                p_div_rad_total_mw=fwbs_variables.p_div_rad_total_mw,
            )

        elif fwbs_variables.i_p_coolant_pumping in [2, 3]:
            # Mechanical pumping power is calculated for first wall and blanket
            self.thermo_hydraulic_model(output=output)
            # For divertor,mechanical pumping power is a fraction of thermal power removed by coolant
            heat_transport_variables.p_div_coolant_pump_mw = (
                heat_transport_variables.f_p_div_coolant_pump_total_heat
                * (
                    physics_variables.p_plasma_separatrix_mw
                    + fwbs_variables.p_div_nuclear_heat_total_mw
                    + fwbs_variables.p_div_rad_total_mw
                )
            )

            # Shield power is negligible and this model doesn't have nuclear heating to the shield
            heat_transport_variables.p_shld_coolant_pump_mw = (
                heat_transport_variables.f_p_shld_coolant_pump_total_heat * 0.0
            )

        if output:
            po.osubhd(self.outfile, "DCLL model: Thermal-hydraulics Component Totals")

            if (fwbs_variables.i_p_coolant_pumping != 2) and (
                fwbs_variables.i_p_coolant_pumping != 3
            ):
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for first wall (MW)",
                    "(p_fw_coolant_pump_mw)",
                    heat_transport_variables.p_fw_coolant_pump_mw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for blanket (MW)",
                    "(p_blkt_coolant_pump_mw)",
                    heat_transport_variables.p_blkt_coolant_pump_mw,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for FW and blanket cooling loop including heat exchanger (MW)",
                    "(p_fw_blkt_coolant_pump_mw)",
                    primary_pumping_variables.p_fw_blkt_coolant_pump_mw,
                    "OP ",
                )

            if fwbs_variables.i_blkt_dual_coolant > 0:
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for liquid metal breeder (MW)",
                    "(p_blkt_breeder_pump_mw)",
                    heat_transport_variables.p_blkt_breeder_pump_mw,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Mechanical pumping power for divertor (MW)",
                "(p_div_coolant_pump_mw)",
                heat_transport_variables.p_div_coolant_pump_mw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Mechanical pumping power for shield and vacuum vessel (MW)",
                "(p_shld_coolant_pump_mw)",
                heat_transport_variables.p_shld_coolant_pump_mw,
                "OP ",
            )

            po.ovarin(
                self.outfile,
                "Switch for plant secondary cycle ",
                "(i_thermal_electric_conversion)",
                fwbs_variables.i_thermal_electric_conversion,
            )
            po.ovarin(
                self.outfile,
                "Switch for plant secondary cycle (liquid metal breeder) ",
                "(secondary_cycle_liq)",
                fwbs_variables.secondary_cycle_liq,
            )
            po.ovarre(
                self.outfile,
                "First wall coolant pressure (Pa)",
                "(pres_fw_coolant)",
                fwbs_variables.pres_fw_coolant,
            )
            po.ovarre(
                self.outfile,
                "Blanket coolant pressure (Pa)",
                "(pres_blkt_coolant)",
                fwbs_variables.pres_blkt_coolant,
            )
            if fwbs_variables.i_blkt_dual_coolant > 0:
                po.ovarre(
                    self.outfile,
                    "Blanket liquid metal breeder pressure (Pa)",
                    "(blpressure_liq)",
                    fwbs_variables.blpressure_liq,
                )

    def dcll_masses(self, output: bool):
        """Material Density Info !

        FW Armour
             - Tungsten
             - Use den_tungsten form constants.f90
        FW and BB Structure Coolant
             - Helium
             - See primary_coolant_properties for denisty etc.
        BB Breeder
             - PbLi
             - See submodule liquid_breeder_properties for density etc.
        Structure
             - EUROFER
             - den_steel in fwbs_variables
        Ceramic FCIs
             - SiC
             - den_ceramic in fwbs_variables


        LT MMS DCLL model [Pal2016]

             Radial Build (m, % vol):

                  FW
                 IB/OB armour = 2.0D-3 m, 100% W
                 IB/OB FW = 1.98D-2 m, 85.54% EUROfer, 14.46% He

                  BZ
                 IB/OB BZ radial stiffening plates total = 6.0D-2 m, 91.33% EUROfer, 8.67% He
                 IB PbLi Channels = 3.0D-1 m, 100% PbLi
                 OB PbLi Channels = 6.4D-1 m, 100% PbLi
                 IB He plena EUROfer walls = 1.0D-1 m, 53% EUROfer, 47% He
                 OB He plena EUROfer walls = 1.7D-1 m, 53% EUROfer, 47% He

                 Back wall = 2.0D-2 m, 85.54% EUROfer, 14.46% He

                 MF/BSS = variable thickness, 51.29% EUROfer, 4.35% He, 44.36% PbLi

             Other info (m, % vol):

                 Side walls = 2.0D-2 m, 85.54% EUROfer, 14.46% He
                 Top walls = 2.0D-2 m, 85.54% EUROfer, 14.46% He
                 Bottom walls = 2.0D-2 m, 85.54% EUROfer, 14.46% He

        Parameters
        ----------
        output: bool

        """
        # If there are FCIs then how much of the radial build is FCI?
        if fwbs_variables.i_blkt_liquid_breeder_channel_type > 0:
            dcll_variables.r_fci = (
                2 * fwbs_variables.nopol * fwbs_variables.th_wall_secondary
            )
        else:
            dcll_variables.r_fci = 0.0

        # Back wall set 0.02m thickness but will vary BZ (structure and breeder) thickness
        dcll_variables.bz_r_ib = build_variables.blbuith - dcll_variables.r_fci
        dcll_variables.bz_r_ob = build_variables.blbuoth - dcll_variables.r_fci
        # Back wall thickness (m)
        dcll_variables.r_backwall = 2.0e-2

        # Manifold/BSS (m) also vars from elsewhere in process but set here
        build_variables.blbmith = (
            build_variables.dr_blkt_inboard
            - dcll_variables.r_backwall
            - build_variables.blbuith
        )
        build_variables.blbmoth = (
            build_variables.dr_blkt_outboard
            - dcll_variables.r_backwall
            - build_variables.blbuoth
        )

        # Fraction of EUROfer (volume composition for EURO + He structures)
        dcll_variables.f_vol_stff_plates = 0.91
        dcll_variables.f_vol_stl_bz_struct = 0.53
        dcll_variables.f_vol_stl_back_wall = 0.86
        dcll_variables.f_vol_stl_fw = 0.86

        # Radial Fraction of BZ Liquid Breeder/Coolant (from DEMO)
        fwbs_variables.r_f_liq_ib = 0.75
        fwbs_variables.r_f_liq_ib = 0.79
        fwbs_variables.w_f_liq_ib = fwbs_variables.r_f_liq_ib
        fwbs_variables.w_f_liq_ob = fwbs_variables.r_f_liq_ib

        # Manifold/BSS Fractions
        dcll_variables.f_vol_mfbss_stl = 0.5129
        dcll_variables.f_vol_mfbss_he = 0.0435
        dcll_variables.f_vol_mfbss_pbli = 0.4436

        # Calculate Volumes
        if fwbs_variables.i_blkt_inboard == 1:
            # IB and OB blanket

            # BZ
            dcll_variables.vol_bz_struct = (
                fwbs_variables.vol_blkt_inboard
                * dcll_variables.bz_r_ib
                * (1 - fwbs_variables.r_f_liq_ib)
                / build_variables.dr_blkt_inboard
            ) + (
                fwbs_variables.vol_blkt_outboard
                * (dcll_variables.bz_r_ob * (1 - fwbs_variables.r_f_liq_ob))
                / build_variables.dr_blkt_outboard
            )
            if fwbs_variables.i_blkt_dual_coolant > 0:
                fwbs_variables.f_a_blkt_cooling_channels = (
                    (1 - dcll_variables.f_vol_stl_bz_struct)
                    * dcll_variables.vol_bz_struct
                ) / fwbs_variables.vol_blkt_total

            dcll_variables.vol_bz_liq = (
                fwbs_variables.vol_blkt_inboard
                * dcll_variables.bz_r_ib
                * fwbs_variables.r_f_liq_ib
                / build_variables.dr_blkt_inboard
            ) + (
                fwbs_variables.vol_blkt_outboard
                * dcll_variables.bz_r_ob
                * fwbs_variables.r_f_liq_ob
                / build_variables.dr_blkt_outboard
            )
            dcll_variables.vol_bz_liq_ib = (
                fwbs_variables.vol_blkt_inboard
                * dcll_variables.bz_r_ib
                * fwbs_variables.r_f_liq_ib
                / build_variables.dr_blkt_inboard
            )
            dcll_variables.vol_bz_liq_ob = (
                fwbs_variables.vol_blkt_outboard
                * dcll_variables.bz_r_ob
                * fwbs_variables.r_f_liq_ob
                / build_variables.dr_blkt_outboard
            )

            if fwbs_variables.i_blkt_liquid_breeder_channel_type > 0:
                dcll_variables.vol_fci = (
                    fwbs_variables.vol_blkt_inboard
                    * dcll_variables.r_fci
                    / build_variables.dr_blkt_inboard
                ) + (
                    fwbs_variables.vol_blkt_outboard
                    * dcll_variables.r_fci
                    / build_variables.dr_blkt_outboard
                )

            # Back Wall
            dcll_variables.vol_bw = (
                fwbs_variables.vol_blkt_inboard
                * dcll_variables.r_backwall
                / build_variables.dr_blkt_inboard
            ) + (
                fwbs_variables.vol_blkt_outboard
                * dcll_variables.r_backwall
                / build_variables.dr_blkt_outboard
            )

            # Manifold/BSS
            dcll_variables.vol_bss = (
                fwbs_variables.vol_blkt_inboard
                * build_variables.blbmith
                / build_variables.dr_blkt_inboard
            ) + (
                fwbs_variables.vol_blkt_outboard
                * build_variables.blbmoth
                / build_variables.dr_blkt_outboard
            )

        else:
            # Only OB blanket

            # BZ
            dcll_variables.vol_bz_struct = (
                fwbs_variables.vol_blkt_outboard
                * dcll_variables.bz_r_ob
                * (1 - fwbs_variables.r_f_liq_ob)
                / build_variables.dr_blkt_outboard
            )
            if fwbs_variables.i_blkt_dual_coolant > 0:
                fwbs_variables.f_a_blkt_cooling_channels = (
                    (1 - dcll_variables.f_vol_stl_bz_struct)
                    * dcll_variables.vol_bz_struct
                ) / fwbs_variables.vol_blkt_total

            dcll_variables.vol_bz_liq = (
                fwbs_variables.vol_blkt_outboard
                * dcll_variables.bz_r_ob
                * fwbs_variables.r_f_liq_ob
                / build_variables.dr_blkt_outboard
            )
            dcll_variables.vol_bz_liq_ob = (
                fwbs_variables.vol_blkt_outboard
                * dcll_variables.bz_r_ob
                * fwbs_variables.r_f_liq_ob
                / build_variables.dr_blkt_outboard
            )
            if fwbs_variables.i_blkt_liquid_breeder_channel_type > 0:
                dcll_variables.vol_fci = (
                    fwbs_variables.vol_blkt_outboard
                    * dcll_variables.r_fci
                    / build_variables.dr_blkt_outboard
                )

            # Back Wall
            dcll_variables.vol_bw = (
                fwbs_variables.vol_blkt_outboard
                * dcll_variables.r_backwall
                / build_variables.dr_blkt_outboard
            )

            # Manifold/BSS
            dcll_variables.vol_bss = (
                fwbs_variables.vol_blkt_outboard
                * build_variables.blbmoth
                / build_variables.dr_blkt_outboard
            )

        # Calculate masses
        # BZ
        dcll_variables.wht_stl_struct = (
            fwbs_variables.den_steel
            * dcll_variables.f_vol_stl_bz_struct
            * dcll_variables.vol_bz_struct
        )
        dcll_variables.wht_cool_struct = (
            fwbs_variables.den_blkt_coolant
            * (1 - dcll_variables.f_vol_stl_bz_struct)
            * dcll_variables.vol_bz_struct
        )
        fwbs_variables.wht_liq = fwbs_variables.den_liq * dcll_variables.vol_bz_liq
        fwbs_variables.wht_liq_ib = fwbs_variables.den_liq * dcll_variables.vol_bz_liq_ib
        fwbs_variables.wht_liq_ob = fwbs_variables.den_liq * dcll_variables.vol_bz_liq_ob
        dcll_variables.wht_cer = fwbs_variables.den_ceramic * dcll_variables.vol_fci
        # Back Wall
        dcll_variables.wht_bw_stl = (
            fwbs_variables.den_steel
            * dcll_variables.f_vol_stl_back_wall
            * dcll_variables.vol_bw
        )
        dcll_variables.wht_bw_cool = (
            fwbs_variables.den_blkt_coolant
            * (1 - dcll_variables.f_vol_stl_back_wall)
            * dcll_variables.vol_bw
        )

        # Manifold/BSS
        dcll_variables.wht_mfbss_stl = (
            fwbs_variables.den_steel
            * dcll_variables.f_vol_mfbss_stl
            * dcll_variables.vol_bss
        )
        dcll_variables.wht_mfbss_cool = (
            fwbs_variables.den_blkt_coolant
            * dcll_variables.f_vol_mfbss_he
            * dcll_variables.vol_bss
        )
        dcll_variables.wht_mfbss_pbli = (
            fwbs_variables.den_liq
            * dcll_variables.f_vol_mfbss_pbli
            * dcll_variables.vol_bss
        )

        # FW
        # First wall volume (m^3)
        fwbs_variables.vol_fw_total = (
            first_wall_variables.a_fw_inboard * build_variables.dr_fw_inboard
            + first_wall_variables.a_fw_outboard * build_variables.dr_fw_outboard
        )
        # First wall mass, excluding armour (kg)
        dcll_variables.fwmass_stl = (
            fwbs_variables.den_steel
            * dcll_variables.f_vol_stl_fw
            * fwbs_variables.vol_fw_total
        )
        dcll_variables.fwmass_cool = (
            fwbs_variables.den_fw_coolant
            * (1 - dcll_variables.f_vol_stl_fw)
            * fwbs_variables.vol_fw_total
        )
        fwbs_variables.m_fw_total = (
            dcll_variables.fwmass_stl + dcll_variables.fwmass_cool
        )
        # First wall armour volume (m^3)
        fwbs_variables.fw_armour_vol = (
            physics_variables.a_plasma_surface * fwbs_variables.fw_armour_thickness
        )
        # First wall armour mass (kg)
        fwbs_variables.fw_armour_mass = (
            fwbs_variables.fw_armour_vol * constants.DEN_TUNGSTEN
        )

        # Total mass of blanket
        fwbs_variables.m_blkt_total = (
            dcll_variables.wht_stl_struct
            + dcll_variables.wht_cool_struct
            + fwbs_variables.wht_liq
            + dcll_variables.wht_bw_stl
            + dcll_variables.wht_bw_cool
            + dcll_variables.wht_mfbss_stl
            + dcll_variables.wht_mfbss_cool
            + dcll_variables.wht_mfbss_pbli
            + dcll_variables.wht_cer
        )

        # Total mass of first wall and blanket
        fwbs_variables.armour_fw_bl_mass = (
            fwbs_variables.fw_armour_mass
            + fwbs_variables.m_fw_total
            + fwbs_variables.m_blkt_total
        )

        # Total mass of IB/OB segment
        if fwbs_variables.i_blkt_inboard == 1:
            dcll_variables.mass_segm_ib = (
                fwbs_variables.m_blkt_total
                * (fwbs_variables.vol_blkt_inboard / fwbs_variables.vol_blkt_total)
                + fwbs_variables.m_fw_total
                * (
                    first_wall_variables.a_fw_inboard
                    * build_variables.dr_fw_inboard
                    / fwbs_variables.vol_fw_total
                )
                + fwbs_variables.fw_armour_mass
                * (
                    (
                        physics_variables.a_plasma_surface
                        - physics_variables.a_plasma_surface_outboard
                    )
                    * fwbs_variables.fw_armour_thickness
                    / fwbs_variables.fw_armour_vol
                )
            ) / fwbs_variables.n_blkt_inboard_modules_toroidal

        dcll_variables.mass_segm_ob = (
            fwbs_variables.m_blkt_total
            * (fwbs_variables.vol_blkt_outboard / fwbs_variables.vol_blkt_total)
            + fwbs_variables.m_fw_total
            * (
                first_wall_variables.a_fw_outboard
                * build_variables.dr_fw_outboard
                / fwbs_variables.vol_fw_total
            )
            + fwbs_variables.fw_armour_mass
            * (
                physics_variables.a_plasma_surface_outboard
                * fwbs_variables.fw_armour_thickness
                / fwbs_variables.fw_armour_vol
            )
        ) / fwbs_variables.n_blkt_outboard_modules_toroidal

        # Total FW/Structure Coolant Mass
        dcll_variables.mass_cool_blanket = (
            dcll_variables.fwmass_cool
            + dcll_variables.wht_cool_struct
            + dcll_variables.wht_bw_cool
            + dcll_variables.wht_mfbss_cool
        )
        # Total Liquid Breeder/Coolant Mass
        dcll_variables.mass_liq_blanket = (
            fwbs_variables.wht_liq + dcll_variables.wht_mfbss_pbli
        )
        # Total Steel Mass
        dcll_variables.mass_stl_blanket = (
            dcll_variables.fwmass_stl
            + dcll_variables.wht_stl_struct
            + dcll_variables.wht_bw_stl
            + dcll_variables.wht_mfbss_stl
        )

        # Mass of material =   density of material * fraction of material by volume * (
        #                      (volume OB blanket * blanket OB zone thickness/ total OB blanket thickness) +
        #                      (volume IB blanket * blanket IB zone thickness/ total IB blanket thickness)

        if output:
            po.osubhd(self.outfile, "DCLL model: Masses")

            po.osubhd(self.outfile, "Component Masses: ")

            po.ovarre(
                self.outfile,
                "First Wall Armour Mass (kg)",
                "(fw_armour_mass)",
                fwbs_variables.fw_armour_mass,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "First Wall Mass, excluding armour (kg)",
                "(m_fw_total)",
                fwbs_variables.m_fw_total,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total Blanket Mass (kg)",
                "(m_blkt_total)",
                fwbs_variables.m_blkt_total,
                "OP ",
            )
            if fwbs_variables.i_blkt_liquid_breeder_channel_type == 1:
                po.ovarre(
                    self.outfile,
                    "Blanket FCI Mass (kg)",
                    "(wht_cer)",
                    dcll_variables.wht_cer,
                    "OP ",
                )
            po.ovarre(
                self.outfile,
                "Total mass of armour, first wall and blanket (kg)",
                "(armour_fw_bl_mass)",
                fwbs_variables.armour_fw_bl_mass,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Total mass for an inboard blanket segment (kg)",
                "(mass_segm_ib)",
                dcll_variables.mass_segm_ib,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total mass for an outboard blanket segment (kg)",
                "(mass_segm_ob)",
                dcll_variables.mass_segm_ob,
                "OP ",
            )

            po.osubhd(self.outfile, "Compositional Masses: ")

            po.ovarre(
                self.outfile,
                "Total FW/Structure Coolant Mass (kg)",
                "(mass_cool_blanket)",
                dcll_variables.mass_cool_blanket,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total Liquid Breeder/Coolant mass (kg)",
                "(mass_liq_blanket)",
                dcll_variables.mass_liq_blanket,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total Steel Mass (FW + Structure) (kg)",
                "(mass_stl_blanket)",
                dcll_variables.mass_stl_blanket,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total W mass (kg)",
                "(fw_armour_mass)",
                fwbs_variables.fw_armour_mass,
                "OP ",
            )

            po.osubhd(self.outfile, "Radial Thickness: ")

            po.ovarrf(
                self.outfile,
                "Inboard radial first wall thickness (m)",
                "(dr_fw_inboard)",
                build_variables.dr_fw_inboard,
            )
            po.ovarrf(
                self.outfile,
                "Outboard radial first wall thickness (m)",
                "(dr_fw_outboard)",
                build_variables.dr_fw_outboard,
            )
            po.ovarrf(
                self.outfile,
                "Inboard radial breeder zone thickness (m)",
                "(blbuith)",
                build_variables.blbuith,
            )
            po.ovarrf(
                self.outfile,
                "Outboard radial breeder zone thickness (m)",
                "(blbuoth)",
                build_variables.blbuoth,
            )

    def write_output(self):
        # Component Volumes
        po.osubhd(self.outfile, "Component Volumes :")

        po.ovarrf(
            self.outfile,
            "First Wall Armour Volume (m3)",
            "(fw_armour_vol)",
            fwbs_variables.fw_armour_vol,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "First Wall Volume (m3)",
            "(vol_fw_total)",
            fwbs_variables.vol_fw_total,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Blanket Volume (m3)",
            "(vol_blkt_total)",
            fwbs_variables.vol_blkt_total,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Shield Volume (m3)",
            "(vol_shld_total)",
            fwbs_variables.vol_shld_total,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Vacuum vessel volume (m3)",
            "(vol_vv)",
            fwbs_variables.vol_vv,
            "OP ",
        )
