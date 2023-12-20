from process.fortran import (
    constants,
    build_variables,
    fwbs_variables,
    dcll_module,
    blanket_library,
    physics_variables,
    current_drive_variables,
    process_output as po,
    primary_pumping_variables,
    heat_transport_variables,
)


class DCLL:
    """This module contains the Dual Coolant Lead Lithium (DCLL) specific submods of PROCESSS.

    author: G. Graham, CCFE

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
             iblanket = 5 * DCLL

         Liquid Metal Breeder Material = PbLi
             i_bb_liq = 0 * Liquid Metal Breeder Material = PbLi

         Specify dual-coolant i.e., get mass flow required from heat extracted from liqid metal breeder
             icooldual = 2

         FIC switch: 0 = no FIC, Eurofer; 1 = FCIs, perfect electrical insulator, 2 = FCIs, with specified conductance
             ifci = 0, 1, or 2

         Liquid metal duct wall conductance initilized at Eurofer value in fwbs_variables, or can input other value, used for ifci = 0 or 2
             (bz_channel_conduct_liq)

         Choose if FW and BB structure are on the same pumping system (unless have diffent coolants), default is same coolant with flow IN->FW->BB->OUT
             (ipump)

         Can set inlet and oulet temperature for liquid metal breeder
             (inlet_temp_liq)
             (outlet_temp_liq)

    References:

         [Nat1995]   Natesan et al. (1995), Assessment of alkali metal coolants for
                     for the ITER blanket, Fusion Engineering and Design 27, 457-466

         [Mal1995]   Malang and Mattas (1995), Comparison of lithium and the eutectic
                     lead-lithium alloy, two candidate liquid metal breeder materials
                     for self-cooled blankets, Fusion Engineering and Design 27, 399-406

         [Gas2001]   Gasior and Mozer (2001), Thermodynamic study of liquid lithium–lead
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

    def __init__(self, blanket_library) -> None:
        self.outfile = constants.nout

        self.blanket_library = blanket_library

    def run(self, output: bool):
        # MDK (27/11/2015)
        build_variables.fwith = 2 * fwbs_variables.afw + 2 * fwbs_variables.fw_wall
        build_variables.fwoth = build_variables.fwith

        blanket_library.component_volumes()
        self.blanket_library.primary_coolant_properties(output=output)
        blanket_library.liquid_breeder_properties(int(output), self.outfile)
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
        """

        if physics_variables.idivrt == 2:
            # Double null configuration
            covf = 1 - (2 * fwbs_variables.fdiv) - fwbs_variables.fhcd
        else:
            # Single null configuration
            covf = 1 - fwbs_variables.fdiv - fwbs_variables.fhcd

        # Nuclear heating in the first wall (MW)
        fwbs_variables.pnucfw = (
            physics_variables.pneutmw * fwbs_variables.pnuc_fw_ratio_dcll * covf
        )

        # Nuclear heating in the blanket with energy multiplication (MW)
        fwbs_variables.pnuc_blkt_ratio_dcll = 1 - fwbs_variables.pnuc_fw_ratio_dcll
        fwbs_variables.pnucblkt = (
            physics_variables.pneutmw
            * fwbs_variables.pnuc_blkt_ratio_dcll
            * fwbs_variables.emult
            * covf
        )

        # Energy multiplication energy (MW)
        fwbs_variables.emultmw = (
            (physics_variables.pneutmw * fwbs_variables.pnuc_blkt_ratio_dcll)
            * (fwbs_variables.emult - 1)
            * covf
        )

        # Divertor

        if physics_variables.idivrt == 2:
            # Double null configuration
            # Nuclear heating in the divertor (MW), neutron power times fdiv
            fwbs_variables.pnucdiv = physics_variables.pneutmw * 2 * fwbs_variables.fdiv
            # Radiation power incident on divertor (MW)
            fwbs_variables.praddiv = physics_variables.pradmw * 2 * fwbs_variables.fdiv
        else:
            # Single null configuration
            # Nuclear heating in the divertor (MW), neutron power times fdiv
            fwbs_variables.pnucdiv = physics_variables.pneutmw * fwbs_variables.fdiv
            # Radiation power incident on divertor (MW)
            fwbs_variables.praddiv = physics_variables.pradmw * fwbs_variables.fdiv

        # HCD Apperatus

        # No nuclear heating of the H & CD
        fwbs_variables.pnuchcd = 0
        # Radiation power incident on HCD apparatus (MW)
        fwbs_variables.pradhcd = physics_variables.pradmw * fwbs_variables.fhcd

        # FW

        # Radiation power incident on first wall (MW)
        fwbs_variables.pradfw = (
            physics_variables.pradmw - fwbs_variables.praddiv - fwbs_variables.pradhcd
        )

        # Surface heat flux on first wall (MW)
        # All of the fast particle losses go to the outer wall.
        fwbs_variables.psurffwo = (
            fwbs_variables.pradfw * build_variables.fwareaob / build_variables.fwarea
            + current_drive_variables.porbitlossmw
            + physics_variables.palpfwmw
        )
        fwbs_variables.psurffwi = fwbs_variables.pradfw * (
            1 - build_variables.fwareaob / build_variables.fwarea
        )

        if output:
            po.osubhd(
                self.outfile, "DCLL model: Nuclear and Radiation Heating of Components"
            )

            po.osubhd(self.outfile, "Component Coverage :")

            po.ovarre(
                self.outfile,
                "Solid angle fraction taken by on divertor",
                "(fdiv)",
                fwbs_variables.fdiv,
            )
            po.ovarre(self.outfile, "Blanket coverage factor", "(covf)", covf)

            po.osubhd(self.outfile, "Nuclear heating :")

            po.ovarre(
                self.outfile,
                "Total nuclear heating in FW (MW)",
                "(pnucfw)",
                fwbs_variables.pnucfw,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Energy multiplication in the blanket",
                "(emult)",
                fwbs_variables.emult,
                "OP ",
            )
            po.ocmmnt(
                self.outfile, "(Note: emult is fixed for this model inside the code)"
            )
            po.ovarre(
                self.outfile,
                "Total nuclear heating in the blanket (including emult) (MW)",
                "(pnucblkt)",
                fwbs_variables.pnucblkt,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total nuclear heating in the shield (MW)",
                "(pnucshld)",
                fwbs_variables.pnucshld,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total nuclear heating in the divertor (MW)",
                "(pnucdiv)",
                fwbs_variables.pnucdiv,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total nuclear heating in TF+PF coils (CS is negligible) (MW)",
                "(ptfnuc)",
                fwbs_variables.ptfnuc,
                "OP ",
            )

            po.osubhd(self.outfile, "Radiation heating :")

            po.ovarrf(
                self.outfile,
                "Radiation heating power into the divertor (MW)",
                "(praddiv)",
                fwbs_variables.praddiv,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Radiation heating power into the first wall (MW)",
                "(pradfw)",
                fwbs_variables.pradfw,
                "OP ",
            )

    def dcll_power_and_heating(self, output: bool):
        # Mechanical Pumping

        # For primary_pumping == 0:
        # User sets mechanical pumping power directly (primary_pumping_power)
        # Values of htpmw_blkt, htpmw_div, htpmw_fw, htpmw_shld set in input file

        if fwbs_variables.primary_pumping == 1:
            # User sets mechanical pumping power as a fraction of thermal power
            # removed by coolant
            heat_transport_variables.htpmw_fw = heat_transport_variables.fpumpfw * (
                fwbs_variables.pnucfw
                + fwbs_variables.psurffwi
                + fwbs_variables.psurffwo
            )
            primary_pumping_variables.htpmw_blkt = (
                heat_transport_variables.fpumpblkt * fwbs_variables.pnucblkt
            )
            # For CCFE HCPB: htpmw_shld = fpumpshld * ( pnucshld + pnuc_cp_sh )
            # Use same as KIT HCLL for now "pnucshld is not available and is very small
            # compared to other powers so set to zero."
            heat_transport_variables.htpmw_shld = (
                heat_transport_variables.fpumpshld * 0.0
            )
            heat_transport_variables.htpmw_div = heat_transport_variables.fpumpdiv * (
                physics_variables.pdivt
                + fwbs_variables.pnucdiv
                + fwbs_variables.praddiv
            )

        elif fwbs_variables.primary_pumping == 2:
            # Mechanical pumping power is calculated for first wall and blanket
            self.blanket_library.thermo_hydraulic_model(output=output)

            # For divertor,mechanical pumping power is a fraction of thermal power removed by coolant
            heat_transport_variables.htpmw_div = heat_transport_variables.fpumpdiv * (
                physics_variables.pdivt
                + fwbs_variables.pnucdiv
                + fwbs_variables.praddiv
            )

            # Shield power is negligible and this model doesn't have nuclear heating to the shield
            heat_transport_variables.htpmw_shld = heat_transport_variables.fpumpshld * 0

        if output:
            po.osubhd(self.outfile, "DCLL model: Thermal-hydraulics Component Totals")

            if (fwbs_variables.primary_pumping != 2) and (
                fwbs_variables.primary_pumping != 3
            ):
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for first wall (MW)",
                    "(htpmw_fw)",
                    heat_transport_variables.htpmw_fw,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for blanket (MW)",
                    "(htpmw_blkt)",
                    heat_transport_variables.htpmw_blkt,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for FW and blanket cooling loop including heat exchanger (MW)",
                    "(htpmw_fw_blkt)",
                    primary_pumping_variables.htpmw_fw_blkt,
                    "OP ",
                )

            if fwbs_variables.icooldual > 0:
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for liquid metal breeder (MW)",
                    "(htpmw_blkt_liq)",
                    heat_transport_variables.htpmw_blkt_liq,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Mechanical pumping power for divertor (MW)",
                "(htpmw_div)",
                heat_transport_variables.htpmw_div,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Mechanical pumping power for shield and vacuum vessel (MW)",
                "(htpmw_shld)",
                heat_transport_variables.htpmw_shld,
                "OP ",
            )

            po.ovarin(
                self.outfile,
                "Switch for plant secondary cycle ",
                "(secondary_cycle)",
                fwbs_variables.secondary_cycle,
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
                "(fwpressure)",
                fwbs_variables.fwpressure,
            )
            po.ovarre(
                self.outfile,
                "Blanket coolant pressure (Pa)",
                "(blpressure)",
                fwbs_variables.blpressure,
            )
            if fwbs_variables.icooldual > 0:
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
             - Use denw from fwbs_variables
        FW and BB Structure Coolant
             - Helium
             - See primary_coolant_properties for denisty etc.
        BB Breeder
             - PbLi
             - See submodule liquid_breeder_properties for density etc.
        Structure
             - EUROFER
             - denstl in fwbs_variables
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
        """
        # If there are FCIs then how much of the radial build is FCI?
        if fwbs_variables.ifci > 0:
            dcll_module.r_fci = (
                2 * fwbs_variables.nopol * fwbs_variables.th_wall_secondary
            )
        else:
            dcll_module.r_fci = 0.0

        # Back wall set 0.02m thickness but will vary BZ (structure and breeder) thickness
        dcll_module.bz_r_ib = build_variables.blbuith - dcll_module.r_fci
        dcll_module.bz_r_ob = build_variables.blbuoth - dcll_module.r_fci
        # Back wall thickness (m)
        dcll_module.r_backwall = 2.0e-2

        # Manifold/BSS (m) also vars from elsewhere in process but set here
        build_variables.blbmith = (
            build_variables.blnkith - dcll_module.r_backwall - build_variables.blbuith
        )
        build_variables.blbmoth = (
            build_variables.blnkoth - dcll_module.r_backwall - build_variables.blbuoth
        )

        # Fraction of EUROfer (volume composition for EURO + He structures)
        dcll_module.f_vol_stff_plates = 0.91
        dcll_module.f_vol_stl_bz_struct = 0.53
        dcll_module.f_vol_stl_back_wall = 0.86
        dcll_module.f_vol_stl_fw = 0.86

        # Radial Fraction of BZ Liquid Breeder/Coolant (from DEMO)
        fwbs_variables.r_f_liq_ib = 0.75
        fwbs_variables.r_f_liq_ib = 0.79
        fwbs_variables.w_f_liq_ib = fwbs_variables.r_f_liq_ib
        fwbs_variables.w_f_liq_ob = fwbs_variables.r_f_liq_ib

        # Manifold/BSS Fractions
        dcll_module.f_vol_mfbss_stl = 0.5129
        dcll_module.f_vol_mfbss_he = 0.0435
        dcll_module.f_vol_mfbss_pbli = 0.4436

        # Calculate Volumes
        if fwbs_variables.iblnkith == 1:
            # IB and OB blanket

            # BZ
            dcll_module.vol_bz_struct = (
                fwbs_variables.volblkti
                * dcll_module.bz_r_ib
                * (1 - fwbs_variables.r_f_liq_ib)
                / build_variables.blnkith
            ) + (
                fwbs_variables.volblkto
                * (dcll_module.bz_r_ob * (1 - fwbs_variables.r_f_liq_ob))
                / build_variables.blnkoth
            )
            if fwbs_variables.icooldual > 0:
                fwbs_variables.vfblkt = (
                    (1 - dcll_module.f_vol_stl_bz_struct) * dcll_module.vol_bz_struct
                ) / fwbs_variables.volblkt

            dcll_module.vol_bz_liq = (
                fwbs_variables.volblkti
                * dcll_module.bz_r_ib
                * fwbs_variables.r_f_liq_ib
                / build_variables.blnkith
            ) + (
                fwbs_variables.volblkto
                * dcll_module.bz_r_ob
                * fwbs_variables.r_f_liq_ob
                / build_variables.blnkoth
            )
            dcll_module.vol_bz_liq_ib = (
                fwbs_variables.volblkti
                * dcll_module.bz_r_ib
                * fwbs_variables.r_f_liq_ib
                / build_variables.blnkith
            )
            dcll_module.vol_bz_liq_ob = (
                fwbs_variables.volblkto
                * dcll_module.bz_r_ob
                * fwbs_variables.r_f_liq_ob
                / build_variables.blnkoth
            )

            if fwbs_variables.ifci > 0:
                dcll_module.vol_fci = (
                    fwbs_variables.volblkti
                    * dcll_module.r_fci
                    / build_variables.blnkith
                ) + (
                    fwbs_variables.volblkto
                    * dcll_module.r_fci
                    / build_variables.blnkoth
                )

            # Back Wall
            dcll_module.vol_bw = (
                fwbs_variables.volblkti
                * dcll_module.r_backwall
                / build_variables.blnkith
            ) + (
                fwbs_variables.volblkto
                * dcll_module.r_backwall
                / build_variables.blnkoth
            )

            # Manifold/BSS
            dcll_module.vol_bss = (
                fwbs_variables.volblkti
                * build_variables.blbmith
                / build_variables.blnkith
            ) + (
                fwbs_variables.volblkto
                * build_variables.blbmoth
                / build_variables.blnkoth
            )

        else:
            # Only OB blanket

            # BZ
            dcll_module.vol_bz_struct = (
                fwbs_variables.volblkto
                * dcll_module.bz_r_ob
                * (1 - fwbs_variables.r_f_liq_ob)
                / build_variables.blnkoth
            )
            if fwbs_variables.icooldual > 0:
                fwbs_variables.vfblkt = (
                    (1 - dcll_module.f_vol_stl_bz_struct) * dcll_module.vol_bz_struct
                ) / fwbs_variables.volblkt

            dcll_module.vol_bz_liq = (
                fwbs_variables.volblkto
                * dcll_module.bz_r_ob
                * fwbs_variables.r_f_liq_ob
                / build_variables.blnkoth
            )
            dcll_module.vol_bz_liq_ob = (
                fwbs_variables.volblkto
                * dcll_module.bz_r_ob
                * fwbs_variables.r_f_liq_ob
                / build_variables.blnkoth
            )
            if fwbs_variables.ifci > 0:
                dcll_module.vol_fci = (
                    fwbs_variables.volblkto
                    * dcll_module.r_fci
                    / build_variables.blnkoth
                )

            # Back Wall
            dcll_module.vol_bw = (
                fwbs_variables.volblkto
                * dcll_module.r_backwall
                / build_variables.blnkoth
            )

            # Manifold/BSS
            dcll_module.vol_bss = (
                fwbs_variables.volblkto
                * build_variables.blbmoth
                / build_variables.blnkoth
            )

        # Calculate masses
        # BZ
        dcll_module.wht_stl_struct = (
            fwbs_variables.denstl
            * dcll_module.f_vol_stl_bz_struct
            * dcll_module.vol_bz_struct
        )
        dcll_module.wht_cool_struct = (
            fwbs_variables.rhof_bl
            * (1 - dcll_module.f_vol_stl_bz_struct)
            * dcll_module.vol_bz_struct
        )
        fwbs_variables.wht_liq = fwbs_variables.den_liq * dcll_module.vol_bz_liq
        fwbs_variables.wht_liq_ib = fwbs_variables.den_liq * dcll_module.vol_bz_liq_ib
        fwbs_variables.wht_liq_ob = fwbs_variables.den_liq * dcll_module.vol_bz_liq_ob
        dcll_module.wht_cer = fwbs_variables.den_ceramic * dcll_module.vol_fci
        # Back Wall
        dcll_module.wht_bw_stl = (
            fwbs_variables.denstl * dcll_module.f_vol_stl_back_wall * dcll_module.vol_bw
        )
        dcll_module.wht_bw_cool = (
            fwbs_variables.rhof_bl
            * (1 - dcll_module.f_vol_stl_back_wall)
            * dcll_module.vol_bw
        )

        # Manifold/BSS
        dcll_module.wht_mfbss_stl = (
            fwbs_variables.denstl * dcll_module.f_vol_mfbss_stl * dcll_module.vol_bss
        )
        dcll_module.wht_mfbss_cool = (
            fwbs_variables.rhof_bl * dcll_module.f_vol_mfbss_he * dcll_module.vol_bss
        )
        dcll_module.wht_mfbss_pbli = (
            fwbs_variables.den_liq * dcll_module.f_vol_mfbss_pbli * dcll_module.vol_bss
        )

        # FW
        # First wall volume (m^3)
        fwbs_variables.volfw = (
            build_variables.fwareaib * build_variables.fwith
            + build_variables.fwareaob * build_variables.fwoth
        )
        # First wall mass, excluding armour (kg)
        dcll_module.fwmass_stl = (
            fwbs_variables.denstl * dcll_module.f_vol_stl_fw * fwbs_variables.volfw
        )
        dcll_module.fwmass_cool = (
            fwbs_variables.rhof_fw
            * (1 - dcll_module.f_vol_stl_fw)
            * fwbs_variables.volfw
        )
        fwbs_variables.fwmass = dcll_module.fwmass_stl + dcll_module.fwmass_cool
        # First wall armour volume (m^3)
        fwbs_variables.fw_armour_vol = (
            physics_variables.sarea * fwbs_variables.fw_armour_thickness
        )
        # First wall armour mass (kg)
        fwbs_variables.fw_armour_mass = (
            fwbs_variables.fw_armour_vol * fwbs_variables.denw
        )

        # Total mass of blanket
        fwbs_variables.whtblkt = (
            dcll_module.wht_stl_struct
            + dcll_module.wht_cool_struct
            + fwbs_variables.wht_liq
            + dcll_module.wht_bw_stl
            + dcll_module.wht_bw_cool
            + dcll_module.wht_mfbss_stl
            + dcll_module.wht_mfbss_cool
            + dcll_module.wht_mfbss_pbli
            + dcll_module.wht_cer
        )

        # Total mass of first wall and blanket
        fwbs_variables.armour_fw_bl_mass = (
            fwbs_variables.fw_armour_mass
            + fwbs_variables.fwmass
            + fwbs_variables.whtblkt
        )

        # Total mass of IB/OB segment
        if fwbs_variables.iblnkith == 1:
            dcll_module.mass_segm_ib = (
                fwbs_variables.whtblkt
                * (fwbs_variables.volblkti / fwbs_variables.volblkt)
                + fwbs_variables.fwmass
                * (
                    build_variables.fwareaib
                    * build_variables.fwith
                    / fwbs_variables.volfw
                )
                + fwbs_variables.fw_armour_mass
                * (
                    (physics_variables.sarea - physics_variables.sareao)
                    * fwbs_variables.fw_armour_thickness
                    / fwbs_variables.fw_armour_vol
                )
            ) / fwbs_variables.nblktmodti

        dcll_module.mass_segm_ob = (
            fwbs_variables.whtblkt * (fwbs_variables.volblkto / fwbs_variables.volblkt)
            + fwbs_variables.fwmass
            * (build_variables.fwareaob * build_variables.fwoth / fwbs_variables.volfw)
            + fwbs_variables.fw_armour_mass
            * (
                physics_variables.sareao
                * fwbs_variables.fw_armour_thickness
                / fwbs_variables.fw_armour_vol
            )
        ) / fwbs_variables.nblktmodto

        # Total FW/Structure Coolant Mass
        dcll_module.mass_cool_blanket = (
            dcll_module.fwmass_cool
            + dcll_module.wht_cool_struct
            + dcll_module.wht_bw_cool
            + dcll_module.wht_mfbss_cool
        )
        # Total Liquid Breeder/Coolant Mass
        dcll_module.mass_liq_blanket = (
            fwbs_variables.wht_liq + dcll_module.wht_mfbss_pbli
        )
        # Total Steel Mass
        dcll_module.mass_stl_blanket = (
            dcll_module.fwmass_stl
            + dcll_module.wht_stl_struct
            + dcll_module.wht_bw_stl
            + dcll_module.wht_mfbss_stl
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
                "(fwmass)",
                fwbs_variables.fwmass,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total Blanket Mass (kg)",
                "(whtblkt)",
                fwbs_variables.whtblkt,
                "OP ",
            )
            if fwbs_variables.ifci == 1:
                po.ovarre(
                    self.outfile,
                    "Blanket FCI Mass (kg)",
                    "(wht_cer)",
                    dcll_module.wht_cer,
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
                dcll_module.mass_segm_ib,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total mass for an outboard blanket segment (kg)",
                "(mass_segm_ob)",
                dcll_module.mass_segm_ob,
                "OP ",
            )

            po.osubhd(self.outfile, "Compositional Masses: ")

            po.ovarre(
                self.outfile,
                "Total FW/Structure Coolant Mass (kg)",
                "(mass_cool_blanket)",
                dcll_module.mass_cool_blanket,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total Liquid Breeder/Coolant mass (kg)",
                "(mass_liq_blanket)",
                dcll_module.mass_liq_blanket,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total Steel Mass (FW + Structure) (kg)",
                "(mass_stl_blanket)",
                dcll_module.mass_stl_blanket,
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
                "(fwith)",
                build_variables.fwith,
            )
            po.ovarrf(
                self.outfile,
                "Outboard radial first wall thickness (m)",
                "(fwoth)",
                build_variables.fwoth,
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
            "(volfw)",
            fwbs_variables.volfw,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Blanket Volume (m3)",
            "(volblkt)",
            fwbs_variables.volblkt,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Shield Volume (m3)",
            "(volshld)",
            fwbs_variables.volshld,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Vacuum vessel volume (m3)",
            "(vdewin)",
            fwbs_variables.vdewin,
            "OP ",
        )
