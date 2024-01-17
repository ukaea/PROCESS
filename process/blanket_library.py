import numpy as np

from process.fortran import (
    constants,
    fwbs_variables,
    process_output as po,
    blanket_library,
    build_variables,
    physics_variables,
    primary_pumping_variables,
    error_handling as eh,
    heat_transport_variables,
)
from process.utilities.f2py_string_patch import f2py_compatible_to_string
from process.coolprop_interface import FluidProperties


class BlanketLibrary:
    def __init__(self, fw) -> None:
        self.outfile = constants.nout

        self.fw = fw

    def primary_coolant_properties(self, output: bool):
        """Calculates the fluid properties of the Primary Coolant in the FW and BZ.
        Uses middle value of input and output temperatures of coolant.
        Curently have H20 and He options.

        original author: P. J. Knight, CCFE
        adapted from previous version of pumppower function by: G Graham, CCFE
        References: see pumppower function description
        """

        # Make sure that, if the inputs for the FW and blanket inputs are different,
        # the ipump variable is appropriately set for seperate coolants
        if (
            f2py_compatible_to_string(fwbs_variables.fwcoolant).title() == "Helium"
            and fwbs_variables.coolwh == 2
        ):
            fwbs_variables.ipump = 1
        if (
            f2py_compatible_to_string(fwbs_variables.fwcoolant).title() == "Water"
            and fwbs_variables.coolwh == 1
        ):
            fwbs_variables.ipump = 1

        # If FW and BB have same coolant...
        if fwbs_variables.ipump == 0:
            # Use FW inlet temp and BB outlet temp
            mid_temp = (fwbs_variables.fwinlet + fwbs_variables.outlet_temp) * 0.5
            # FW/BB
            fw_bb_fluid_properties = FluidProperties.of(
                f2py_compatible_to_string(fwbs_variables.fwcoolant),
                temperature=mid_temp,
                pressure=fwbs_variables.fwpressure.item(),
            )
            fwbs_variables.rhof_fw = fw_bb_fluid_properties.density
            fwbs_variables.cp_fw = fw_bb_fluid_properties.specific_heat_const_p
            fwbs_variables.cv_fw = fw_bb_fluid_properties.specific_heat_const_v
            fwbs_variables.visc_fw = fw_bb_fluid_properties.viscosity

            fwbs_variables.rhof_bl = fwbs_variables.rhof_fw
            fwbs_variables.visc_bl = fwbs_variables.visc_fw
            fwbs_variables.cp_bl = fwbs_variables.cp_fw
            fwbs_variables.cv_bl = fwbs_variables.cv_fw

        # If FW and BB have different coolants...
        else:
            # FW
            mid_temp_fw = (fwbs_variables.fwinlet + fwbs_variables.fwoutlet) * 0.5
            fw_fluid_properties = FluidProperties.of(
                f2py_compatible_to_string(fwbs_variables.fwcoolant),
                temperature=mid_temp_fw,
                pressure=fwbs_variables.fwpressure,
            )
            fwbs_variables.rhof_fw = fw_fluid_properties.density
            fwbs_variables.cp_fw = fw_fluid_properties.specific_heat_const_p
            fwbs_variables.cv_fw = fw_fluid_properties.specific_heat_const_v
            fwbs_variables.visc_fw = fw_fluid_properties.viscosity

            # BB
            mid_temp_bl = (fwbs_variables.inlet_temp + fwbs_variables.outlet_temp) * 0.5
            bb_fluid_properties = FluidProperties.of(
                "Helium" if fwbs_variables.coolwh == 1 else "Water",
                temperature=mid_temp_bl,
                pressure=fwbs_variables.blpressure,
            )
            fwbs_variables.rhof_bl = bb_fluid_properties.density
            fwbs_variables.cp_bl = bb_fluid_properties.specific_heat_const_p
            fwbs_variables.cv_bl = bb_fluid_properties.specific_heat_const_v
            fwbs_variables.visc_bl = bb_fluid_properties.viscosity

        if (
            fwbs_variables.rhof_fw > 1e9
            or fwbs_variables.rhof_fw <= 0
            or np.isnan(fwbs_variables.rhof_fw)
        ):
            raise RuntimeError(
                f"Error in primary_coolant_properties. {fwbs_variables.rhof_fw = }"
            )
        if (
            fwbs_variables.rhof_bl > 1e9
            or fwbs_variables.rhof_bl <= 0
            or np.isnan(fwbs_variables.rhof_bl)
        ):
            raise RuntimeError(
                f"Error in primary_coolant_properties. {fwbs_variables.rhof_bl = }"
            )

        if output:
            po.oheadr(
                self.outfile, "First wall and blanket : (Primary) Coolant Properties"
            )
            po.ocmmnt(
                self.outfile,
                "Calculated using mid temp(s) of system (or systems if use different collant types).",
            )

            # FW (or FW/BB)!!!!!!!!!!!!!!!!!!
            if fwbs_variables.ipump == 1:
                po.osubhd(self.outfile, "First Wall :")

            po.ovarst(
                self.outfile,
                "Coolant type",
                "(fwcoolant)",
                f'"{fwbs_variables.fwcoolant}"',
            )
            po.ovarrf(
                self.outfile,
                "Density (kg m-3)",
                "(rhof_fw)",
                fwbs_variables.rhof_fw,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Viscosity (Pa s)",
                "(visc_fw)",
                fwbs_variables.visc_fw,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Inlet Temperature (Celcius)",
                "(fwinlet)",
                fwbs_variables.fwinlet,
                "OP ",
            )

            if fwbs_variables.ipump == 0:
                po.ovarre(
                    self.outfile,
                    "Outlet Temperature (Celcius)",
                    "(outlet_temp)",
                    fwbs_variables.outlet_temp,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Outlet Temperature (Celcius)",
                    "(fwoutlet)",
                    fwbs_variables.fwoutlet,
                    "OP ",
                )

            # BB !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if fwbs_variables.ipump == 1:
                po.osubhd(self.outfile, "Breeding Blanket :")

                if fwbs_variables.coolwh == 1:
                    po.ocmmnt(self.outfile, "Coolant type (coolwh=1), Helium")
                if fwbs_variables.coolwh == 2:
                    po.ocmmnt(self.outfile, "Coolant type (coolwh=2), Water")
                po.ovarrf(
                    self.outfile,
                    "Density (kg m-3)",
                    "(rhof_bl)",
                    fwbs_variables.rhof_bl,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Viscosity (Pa s)",
                    "(visc_bl)",
                    fwbs_variables.visc_bl,
                    "OP ",
                )

                po.ovarre(
                    self.outfile,
                    "Inlet Temperature (Celcius)",
                    "(inlet_temp)",
                    fwbs_variables.inlet_temp,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Outlet Temperature (Celcius)",
                    "(outlet_temp)",
                    fwbs_variables.outlet_temp,
                    "OP ",
                )

    def thermo_hydraulic_model(self, output: bool):
        """Thermo-hydraulic model for first wall and blanket
        ONLY CALLED if primary_pumping = 2

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

            length (m)                  fw_channel_length
            width (m)                   afw (radius, cicular)   afw                 a_bz_liq, b_bz_liq (rectangular)
            wall thickness (m)          fw_wall                 fw_wall             th_wall_secondary
            pitch (m)                   pitch
            roughness epsilon           roughness
            peak FW temp (K)            tpeak
            maximum temp (K)            tfwmatmax
            FCI switch                  ---                     ---                 ifci

            Coolant                     FW                      BB primary          BB secondary

            primary coolant switch      fwcoolant               coolwh              ---
            secondary coolant switch    ---                     ---                 i_bb_liq
            inlet temp (K)              fwinlet                 inlet_temp          inlet_temp_liq
            outlet temp (K)             fwoutlet                outlet_temp         outlet_temp_liq
            pressure (Pa)               fwpressure              blpressure          blpressure_liq
        """
        npoltoti = 0
        npoltoto = 0
        npblkti_liq = 0
        npblkto_liq = 0

        if fwbs_variables.iblanket == 5:
            # Unless DCLL then we will use BZ
            blanket_library.bldepti = build_variables.blbuith
            blanket_library.bldepto = build_variables.blbuoth
        else:
            blanket_library.bldepti = 0.8e0 * build_variables.blnkith
            blanket_library.bldepto = 0.8e0 * build_variables.blnkoth

        # Using the total perimeter of the machine, segment the outboard
        # blanket into nblktmodp*nblktmodt modules, all assumed to be the same size

        # If SMS blanket then do not have seperate poloidal modules....
        # Should not need this as nblktmodpi is input but make sure here.
        if fwbs_variables.ims == 1:
            fwbs_variables.nblktmodpi = 1
            fwbs_variables.nblktmodpo = 1

        # Calculate mid-plane toroidal circumference and segment
        blanket_library.blwidti = (
            2.0e0
            * np.pi
            * (
                physics_variables.rmajor
                - physics_variables.rminor
                - build_variables.scrapli
            )
        ) / fwbs_variables.nblktmodti
        blanket_library.blwidto = (
            2.0e0
            * np.pi
            * (
                physics_variables.rmajor
                + physics_variables.rminor
                + build_variables.scraplo
            )
        ) / fwbs_variables.nblktmodto

        # Calculate poloidal height of blanket modules
        blanket_library.blanket_mod_pol_height()

        if fwbs_variables.icooldual > 0:
            # Use smallest space available to pipes for pipe sizes in pumping calculations (worst case)
            if fwbs_variables.iblnkith == 1:
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
        if fwbs_variables.icooldual == 2:
            no90bz = 4
            no180bz = 1
            no90bz_liq = 2
            no180bz_liq = 1
        elif fwbs_variables.icooldual == 1:
            no90bz = 4
            no180bz = 1
            no90bz_liq = 2
            no180bz_liq = 1
        else:
            no90bz = 4
            no180bz = 1

        # Nuclear Power Deposited #######

        # IB/OB FW (MW)
        blanket_library.pnucfwi = (
            fwbs_variables.pnucfw * build_variables.fwareaib / build_variables.fwarea
        )
        blanket_library.pnucfwo = (
            fwbs_variables.pnucfw * build_variables.fwareaob / build_variables.fwarea
        )

        # IB/OB Blanket (MW)

        # Neutron power deposited in inboard blanket (MW)
        if fwbs_variables.iblnkith == 1:
            blanket_library.pnucblkti = (
                fwbs_variables.pnucblkt
                * fwbs_variables.volblkti
                / fwbs_variables.volblkt
            )

        # Neutron power deposited in outboard blanket (MW)
        blanket_library.pnucblkto = (
            fwbs_variables.pnucblkt * fwbs_variables.volblkto / fwbs_variables.volblkt
        )

        # For a dual-coolant blanket, some fraction of the power goes into the
        # structure of the BZ and is cooled by the primary coolant, and some fraction
        # goes into the liquid breeder to be cooled by itself.

        # If the blanket is dual-coolant...
        if fwbs_variables.icooldual == 2:
            f_nuc_pow_bz_liq = 1 - fwbs_variables.f_nuc_pow_bz_struct

            # Inboard blanket calc. Will return 0 if no inboard shldith thickness
            pnucblkti_struct = (
                fwbs_variables.pnucblkt * fwbs_variables.f_nuc_pow_bz_struct
            ) * (fwbs_variables.volblkti / fwbs_variables.volblkt)
            pnucblkti_liq = (fwbs_variables.pnucblkt * f_nuc_pow_bz_liq) * (
                fwbs_variables.volblkti / fwbs_variables.volblkt
            )
            pnucblkto_struct = (
                fwbs_variables.pnucblkt * fwbs_variables.f_nuc_pow_bz_struct
            ) * (fwbs_variables.volblkto / fwbs_variables.volblkt)
            pnucblkto_liq = (fwbs_variables.pnucblkt * f_nuc_pow_bz_liq) * (
                fwbs_variables.volblkto / fwbs_variables.volblkt
            )

        # FW and BB Mass Flow ###########

        # Make sure that, if the inputs for the FW and blanket inputs are different,
        # the ipump variable is appropriately set for seperate coolants
        if (
            f2py_compatible_to_string(fwbs_variables.fwcoolant).title() == "Helium"
            and fwbs_variables.coolwh == 2
        ):
            fwbs_variables.ipump = 1
        if (
            f2py_compatible_to_string(fwbs_variables.fwcoolant).title() == "Water"
            and fwbs_variables.coolwh == 1
        ):
            fwbs_variables.ipump = 1

        # If FW and BB have the same coolant...
        if fwbs_variables.ipump == 0:
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
            if fwbs_variables.iblnkith == 1:
                fwoutleti = (f_nuc_fwi * fwbs_variables.outlet_temp) + (
                    1 - f_nuc_fwi
                ) * fwbs_variables.fwinlet
                inlet_tempi = fwoutleti
            else:
                fwoutleti = fwbs_variables.fwoutlet

            fwoutleto = (f_nuc_fwo * fwbs_variables.outlet_temp) + (
                1 - f_nuc_fwo
            ) * fwbs_variables.fwinlet
            inlet_tempo = fwoutleto

        elif fwbs_variables.ipump == 1:
            fwoutleti = fwbs_variables.fwoutlet
            inlet_tempi = fwbs_variables.inlet_temp
            fwoutleto = fwbs_variables.fwoutlet
            inlet_tempo = fwbs_variables.inlet_temp

        # Maximum FW temperature. (27/11/2015) Issue #348
        # First wall flow is just along the first wall, with no allowance for radial
        # pipes, manifolds etc. The outputs are mid quantities of inlet and outlet.
        # This subroutine recalculates cp and rhof.
        (blanket_library.tpeakfwi, cf, rhof, blanket_library.mffwpi,) = self.fw.fw_temp(
            output,
            fwbs_variables.afw,
            build_variables.fwith,
            build_variables.fwareaib,
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
        #     fwbs_variables.afw,
        #     build_variables.fwith,
        #     build_variables.fwareaib,
        #     fwbs_variables.psurffwi,
        #     blanket_library.pnucfwi,
        #     "Inboard first wall",
        # )
        (fwbs_variables.tpeakfwo, cf, rhof, fwbs_variables.mffwpo) = self.fw.fw_temp(
            output,
            fwbs_variables.afw,
            build_variables.fwoth,
            build_variables.fwareaob,
            fwbs_variables.psurffwo,
            blanket_library.pnucfwo,
            "Outboard first wall",
        )
        # (fwbs_variables.tpeakfwo, cf, rhof, fwbs_variables.mffwpo) = fw_module.fw_temp(
        #     int(output),
        #     self.outfile,
        #     fwbs_variables.afw,
        #     build_variables.fwoth,
        #     build_variables.fwareaob,
        #     fwbs_variables.psurffwo,
        #     blanket_library.pnucfwo,
        #     "Outboard first wall",
        # )

        # Peak first wall temperature (K)
        fwbs_variables.tpeak = max(blanket_library.tpeakfwi, blanket_library.tpeakfwo)

        # Total mass flow rate to remove inboard FW power (kg/s)
        blanket_library.mffwi = (
            1.0e6
            * (blanket_library.pnucfwi + fwbs_variables.psurffwi)
            / (fwbs_variables.cp_fw * (fwoutleti - fwbs_variables.fwinlet))
        )
        # Total mass flow rate to remove outboard FW power (kg/s)
        blanket_library.mffwo = (
            1.0e6
            * (blanket_library.pnucfwo + fwbs_variables.psurffwo)
            / (fwbs_variables.cp_fw * (fwoutleto - fwbs_variables.fwinlet))
        )

        # If the blanket is dual-coolant...
        if fwbs_variables.icooldual == 2:
            # Mass flow rates for outboard blanket coolants (kg/s)
            blanket_library.mfblkto = (
                1.0e6
                * (pnucblkto_struct)
                / (fwbs_variables.cp_bl * (fwbs_variables.outlet_temp - inlet_tempo))
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
            if fwbs_variables.iblnkith == 1:
                # Mass flow rates for inboard blanket coolants (kg/s)
                blanket_library.mfblkti = (
                    1.0e6
                    * (pnucblkti_struct)
                    / (
                        fwbs_variables.cp_bl
                        * (fwbs_variables.outlet_temp - inlet_tempi)
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
                / (fwbs_variables.cp_bl * (fwbs_variables.outlet_temp - inlet_tempo))
            )

            # Get mass flow rate etc. for inboard blanket breeder flow for tritium extraction
            # Use the number of desired recirculations ([Aub2013]=10) and mass from dcll_masses
            # N.B. wht_liq is BZ mass, does not include manifold.
            blanket_library.mfblkto_liq = (
                fwbs_variables.n_liq_recirc * fwbs_variables.wht_liq_ob
            ) / (24 * 3600)

            # If there is an IB blanket...
            if fwbs_variables.iblnkith == 1:
                # Mass flow rate for inboard blanket coolant (kg/s)
                blanket_library.mfblkti = (
                    1.0e6
                    * (blanket_library.pnucblkti)
                    / (
                        fwbs_variables.cp_bl
                        * (fwbs_variables.outlet_temp - inlet_tempi)
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
                / (fwbs_variables.cp_bl * (fwbs_variables.outlet_temp - inlet_tempo))
            )

            # If there is an IB blanket...
            # Mass flow rate for inboard blanket coolant (kg/s)
            if fwbs_variables.iblnkith == 1:
                blanket_library.mfblkti = (
                    1.0e6
                    * (blanket_library.pnucblkti)
                    / (
                        fwbs_variables.cp_bl
                        * (fwbs_variables.outlet_temp - inlet_tempi)
                    )
                )

        # FW Pipe Flow and Velocity ######

        # Total number of first wall pipes from channel length and pitch (02/12/2015)
        blanket_library.npfwi = build_variables.fwareaib / (
            fwbs_variables.fw_channel_length * fwbs_variables.pitch
        )
        blanket_library.npfwo = build_variables.fwareaob / (
            fwbs_variables.fw_channel_length * fwbs_variables.pitch
        )

        # Mass flow rate per FW coolant pipe (kg/s):
        blanket_library.mffwpi = blanket_library.mffwi / blanket_library.npfwi
        blanket_library.mffwpo = blanket_library.mffwo / blanket_library.npfwo

        # Coolant velocite in FW (m/s)
        velfwi = blanket_library.flow_velocity(
            i_channel_shape=1,
            mass_flow_rate=blanket_library.mffwpi,
            flow_density=fwbs_variables.rhof_fw,
        )
        velfwo = blanket_library.flow_velocity(
            i_channel_shape=1,
            mass_flow_rate=blanket_library.mffwpo,
            flow_density=fwbs_variables.rhof_fw,
        )

        # If the blanket is dual-coolant...
        if fwbs_variables.icooldual == 2:
            # Calc total num of pipes (in all inboard modules) from
            # coolant frac and channel dimensions
            # Assumes up/down flow, two 90 deg bends per length
            blanket_library.npblkto = (
                fwbs_variables.vfblkt * fwbs_variables.volblkto
            ) / (
                np.pi
                * fwbs_variables.afw
                * fwbs_variables.afw
                * blanket_library.bzfllengo
            )
            npblkto_liq = (
                fwbs_variables.nopipes
                * fwbs_variables.nblktmodto
                * fwbs_variables.nblktmodpo
            )

            # Mass flow rate per coolant pipe
            blanket_library.mfblktpo = blanket_library.mfblkto / blanket_library.npblkto
            mfblktpo_liq = blanket_library.mfblkto_liq / npblkto_liq
            # Coolant velocites in blanket (m/s)
            # Assume BZ structure has same channel width as FW
            blanket_library.velblkto = blanket_library.flow_velocity(
                i_channel_shape=1,
                mass_flow_rate=blanket_library.mfblktpo,
                flow_density=fwbs_variables.rhof_bl,
            )
            velblkto_liq = blanket_library.flow_velocity(
                i_channel_shape=2,
                mass_flow_rate=mfblktpo_liq,
                flow_density=fwbs_variables.den_liq,
            )

            if fwbs_variables.iblnkith == 1:
                # Calc total num of pipes (in all inboard modules) from
                # coolant frac and channel dimensions
                # Assumes up/down flow, two 90 deg bends per length
                blanket_library.npblkti = (
                    fwbs_variables.vfblkt * fwbs_variables.volblkti
                ) / (
                    np.pi
                    * fwbs_variables.afw
                    * fwbs_variables.afw
                    * blanket_library.bzfllengi
                )
                # Have DEMO DCLL set here for now
                npblkti_liq = (
                    fwbs_variables.nopipes
                    * fwbs_variables.nblktmodti
                    * fwbs_variables.nblktmodpi
                )

                # Mass flow rate per coolant pipe
                blanket_library.mfblktpi = (
                    blanket_library.mfblkti / blanket_library.npblkti
                )
                blanket_library.mfblktpi_liq = blanket_library.mfblkti_liq / npblkti_liq

                # Coolant velocites in blanket (m/s)
                # Assume BZ structure has same channel width as FW
                blanket_library.velblkti = blanket_library.flow_velocity(
                    i_channel_shape=1,
                    mass_flow_rate=blanket_library.mfblktpi,
                    flow_density=fwbs_variables.rhof_bl,
                )
                velblkti_liq = blanket_library.flow_velocity(
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
                fwbs_variables.vfblkt * fwbs_variables.volblkto
            ) / (
                np.pi
                * fwbs_variables.afw
                * fwbs_variables.afw
                * blanket_library.bzfllengo
            )
            npblkto_liq = (
                fwbs_variables.nopipes
                * fwbs_variables.nblktmodto
                * fwbs_variables.nblktmodpo
            )

            # Mass flow rate per coolant pipe
            blanket_library.mfblktpo = blanket_library.mfblkto / blanket_library.npblkto

            # Coolant velocity in blanket (m/s)
            # Assume BZ structure has same channel width as FW
            blanket_library.velblkto = blanket_library.flow_velocity(
                i_channel_shape=1,
                mass_flow_rate=blanket_library.mfblktpo,
                flow_density=fwbs_variables.rhof_bl,
            )

            # Get mass flow rate etc. for inboard blanket breeder flow for tritium extraction
            # Use the number of desired recirculations ([Aub2013]=10) and mass from dcll_masses
            # N.B. wht_liq is BZ mass, does not include manifold.
            blanket_library.mfblkto_liq = (
                fwbs_variables.n_liq_recirc * fwbs_variables.wht_liq_ob
            ) / (24 * 3600)
            blanket_library.mfblktpo_liq = blanket_library.mfblkto_liq / npblkto_liq
            velblkto_liq = blanket_library.flow_velocity(
                i_channel_shape=2,
                mass_flow_rate=blanket_library.mfblktpo_liq,
                flow_density=fwbs_variables.den_liq,
            )

            if fwbs_variables.iblnkith == 1:
                # Calc total num of pipes (in all inboard modules) from
                # coolant frac and channel dimensions
                # Assumes up/down flow, two 90 deg bends per length
                blanket_library.npblkti = (
                    fwbs_variables.vfblkt * fwbs_variables.volblkti
                ) / (
                    np.pi
                    * fwbs_variables.afw
                    * fwbs_variables.afw
                    * blanket_library.bzfllengi
                )
                # Have DEMO DCLL set here for now
                npblkti_liq = (
                    fwbs_variables.nopipes
                    * fwbs_variables.nblktmodti
                    * fwbs_variables.nblktmodpi
                )

                # Mass flow rate per coolant pipe
                blanket_library.mfblktpi = (
                    blanket_library.mfblkti / blanket_library.npblkti
                )

                # Coolant velocity in blanket (m/s)
                # Assume BZ structure has same channel width as FW
                blanket_library.velblkti = blanket_library.flow_velocity(
                    i_channel_shape=1,
                    mass_flow_rate=fwbs_variables.mfblktpi,
                    flow_density=fwbs_variables.rhof_bl,
                )

                # Get mass flow rate etc. for inboard blanket breeder flow for tritium extraction
                # Use the number of desired recirculations ([Aub2013]=10) and mass from dcll_masses
                # N.B. wht_liq is BZ mass, does not include manifold.
                blanket_library.mfblkti_liq = (
                    fwbs_variables.n_liq_recirc * fwbs_variables.wht_liq_ib
                ) / (24 * 3600)
                blanket_library.mfblktpi_liq = fwbs_variables.mfblkti_liq / npblkti_liq
                velblkti_liq = blanket_library.flow_velocity(
                    i_channel_shape=2,
                    mass_flow_rate=blanket_library.mfblktpi_liq,
                    flow_density=fwbs_variables.den_liq,
                )

        # If the blanket is single-coolant with solid breeder...
        else:
            # Calculate total number of pipes (in all outboard modules) from coolant fraction and
            # channel dimensions (assumes up/down flow, two 90 deg bends per length)
            blanket_library.npblkto = (
                fwbs_variables.vfblkt * fwbs_variables.volblkto
            ) / (
                np.pi
                * fwbs_variables.afw
                * fwbs_variables.afw
                * blanket_library.bzfllengo
            )

            # Mass flow rate per coolant pipe
            blanket_library.mfblktpo = blanket_library.mfblkto / blanket_library.npblkto

            # Coolant velocity in blanket (m/s)
            # Assume BZ structure has same channel width as FW
            blanket_library.velblkto = blanket_library.flow_velocity(
                i_channel_shape=1,
                mass_flow_rate=blanket_library.mfblktpo,
                flow_density=fwbs_variables.rhof_bl,
            )

            if fwbs_variables.iblnkith == 1:
                # Calc total num of pipes (in all inboard modules) from
                # coolant frac and channel dimensions
                # Assumes up/down flow, two 90 deg bends per length
                blanket_library.npblkti = (
                    fwbs_variables.vfblkt * fwbs_variables.volblkti
                ) / (
                    np.pi
                    * fwbs_variables.afw
                    * fwbs_variables.afw
                    * blanket_library.bzfllengi
                )

                # Mass flow rate per coolant pipe
                blanket_library.mfblktpi = (
                    blanket_library.mfblkti / blanket_library.npblkti
                )

                # Coolant velocity in blanket (m/s)
                # Assume BZ structure has same channel width as FW
                blanket_library.velblkti = blanket_library.flow_velocity(
                    i_channel_shape=1,
                    mass_flow_rate=blanket_library.mfblktpi,
                    flow_density=fwbs_variables.rhof_bl,
                )

        # FW Presure Drops ###############

        deltap_fwi = self.deltap_tot(
            output,
            icoolpump=1,
            flow_velocity=velfwi,
            flleng=fwbs_variables.fw_channel_length,
            no90=no90fw,
            no180=no180fw,
            coolant_density=fwbs_variables.rhof_fw,
            coolant_dynamic_viscosity=fwbs_variables.visc_fw,
            coolant_electrical_conductivity=0.0e0,
            pol_channel_length=pollengi,
            nopolchan=npoltoti,
            label="Inboard first wall",
        )

        deltap_fwo = self.deltap_tot(
            output,
            icoolpump=1,
            flow_velocity=velfwo,
            flleng=fwbs_variables.fw_channel_length,
            no90=no90fw,
            no180=no180fw,
            coolant_density=fwbs_variables.rhof_fw,
            coolant_dynamic_viscosity=fwbs_variables.visc_fw,
            coolant_electrical_conductivity=0.0e0,
            pol_channel_length=pollengo,
            nopolchan=npoltoto,
            label="Outboard first wall",
        )

        # BB Presure Drops ###############

        # Long polodal flows
        if fwbs_variables.iblnkith == 1:
            npoltoti = fwbs_variables.nopol * npblkti_liq
        npoltoto = fwbs_variables.nopol * npblkto_liq

        deltap_blo = self.deltap_tot(
            output,
            icoolpump=1,
            flow_velocity=blanket_library.velblkto,
            flleng=blanket_library.bzfllengo,
            no90=no90bz,
            no180=no180bz,
            coolant_density=fwbs_variables.rhof_bl,
            coolant_dynamic_viscosity=fwbs_variables.visc_bl,
            coolant_electrical_conductivity=0.0e0,
            pol_channel_length=pollengo,
            nopolchan=npoltoto,
            label="Outboard blanket",
        )

        if fwbs_variables.iblnkith == 1:
            deltap_bli = self.deltap_tot(
                output,
                icoolpump=1,
                flow_velocity=blanket_library.velblkti,
                flleng=blanket_library.bzfllengi,
                no90=no90bz,
                no180=no180bz,
                coolant_density=fwbs_variables.rhof_bl,
                coolant_dynamic_viscosity=fwbs_variables.visc_bl,
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
            if fwbs_variables.iblnkith == 1:
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

        # Pumping Power

        # If FW and BB have the same coolant...
        if fwbs_variables.ipump == 0:
            # Total pressure drop in the first wall/blanket  (Pa)
            if fwbs_variables.iblnkith == 1:
                deltap_fw_blkt = deltap_fwi + deltap_bli + deltap_fwo + deltap_blo
            if fwbs_variables.iblnkith == 0:
                deltap_fw_blkt = deltap_fwi + deltap_fwo + deltap_blo

            # Total coolant mass flow rate in the first wall/blanket (kg/s)
            blanket_library.mftotal = blanket_library.mffwi + blanket_library.mffwo

            # Total mechanical pumping power (MW)
            primary_pumping_variables.htpmw_fw_blkt = self.pumppower(
                output=output,
                icoolpump=1,
                temp_in=fwbs_variables.fwinlet.item(),
                temp_out=fwbs_variables.outlet_temp.item(),
                pressure=fwbs_variables.fwpressure.item(),
                pdrop=deltap_fw_blkt,
                mf=blanket_library.mftotal,
                primary_coolant_switch=f2py_compatible_to_string(
                    fwbs_variables.fwcoolant
                ),
                coolant_density=fwbs_variables.rhof_fw,
                label="First Wall and Blanket",
            )

        # If FW and BB have different coolants...
        elif fwbs_variables.ipump == 1:
            # Total pressure drop in the first wall (Pa)
            deltap_fw = deltap_fwi + deltap_fwo

            # Total pressure drop in the blanket (Pa)
            if fwbs_variables.iblnkith == 1:
                deltap_blkt = deltap_bli + deltap_blo
            if fwbs_variables.iblnkith == 0:
                deltap_blkt = deltap_blo

            # Total coolant mass flow rate in the first wall (kg/s)
            blanket_library.mffw = blanket_library.mffwi + blanket_library.mffwo
            # Total coolant mass flow rate in the blanket (kg/s)
            blanket_library.mfblkt = blanket_library.mfblkti + blanket_library.mfblkto

            # Mechanical pumping power for the first wall (MW)
            heat_transport_variables.htpmw_fw = self.pumppower(
                output=output,
                icoolpump=1,
                temp_in=fwbs_variables.fwinlet.item(),
                temp_out=fwbs_variables.fwoutlet.item(),
                pressure=fwbs_variables.fwpressure.item(),
                pdrop=deltap_fw.item(),
                mf=blanket_library.mffw,
                primary_coolant_switch=f2py_compatible_to_string(
                    fwbs_variables.fwcoolant
                ),
                coolant_density=fwbs_variables.rhof_fw,
                label="First Wall",
            )

            # Mechanical pumping power for the blanket (MW)
            heat_transport_variables.htpmw_blkt = self.pumppower(
                output=output,
                icoolpump=1,
                temp_in=fwbs_variables.inlet_temp.item(),
                temp_out=fwbs_variables.outlet_temp.item(),
                pressure=fwbs_variables.blpressure.item(),
                pdrop=deltap_blkt.item(),
                mf=blanket_library.mfblkt,
                primary_coolant_switch="Helium"
                if fwbs_variables.coolwh == 1
                else "Water",
                coolant_density=blanket_library.rhof_bl,
                label="Blanket",
            )

            # Total mechanical pumping power (MW)
            primary_pumping_variables.htpmw_fw_blkt = (
                heat_transport_variables.htpmw_fw + heat_transport_variables.htpmw_blkt
            )

        # If the blanket has a liquid metal breeder...
        if fwbs_variables.icooldual > 0:
            # Total pressure drop in the blanket (Pa)
            if fwbs_variables.iblnkith == 1:
                deltap_bl_liq = deltap_bli_liq + deltap_blo_liq
            if fwbs_variables.iblnkith == 0:
                deltap_bl_liq = deltap_blo_liq

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
                primary_coolant_switch="Helium"
                if fwbs_variables.coolwh == 1
                else "Water",
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

            # FW !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            po.osubhd(self.outfile, "First wall: ")

            po.ovarst(
                self.outfile,
                "First wall coolant type",
                "(fwcoolant)",
                f'"{fwbs_variables. fwcoolant}"',
            )
            po.ovarre(
                self.outfile,
                "Wall thickness of first wall cooling channels (m)",
                "(fw_wall)",
                fwbs_variables.fw_wall,
            )
            po.ovarre(
                self.outfile,
                "Radius of first wall cooling channels (m)",
                "(afw)",
                fwbs_variables.afw,
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
                "(fwinlet)",
                fwbs_variables.fwinlet,
            )
            po.ovarrf(
                self.outfile,
                "Outlet temperature of first wall coolant (K)",
                "(fwoutlet)",
                fwbs_variables.fwoutlet,
            )
            po.ovarre(
                self.outfile,
                "First wall coolant pressure (Pa)",
                "(fwpressure)",
                fwbs_variables.fwpressure,
            )
            if fwbs_variables.ipump == 1:
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
                "(tfwmatmax)",
                fwbs_variables.tfwmatmax,
            )
            po.ovarrf(
                self.outfile,
                "Actual peak temperature of first wall material (K)",
                "(tpeak)",
                fwbs_variables.tpeak,
                "OP ",
            )

            # BB !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            po.osubhd(self.outfile, "Breeding Blanket (primary): ")

            po.ovarin(
                self.outfile,
                "Blanket coolant type (1=He, 2=H20)",
                "(coolwh)",
                fwbs_variables.coolwh,
            )
            po.ovarrf(
                self.outfile,
                "Inlet temperature of blanket coolant (K)",
                "(inlet_temp)",
                fwbs_variables.inlet_temp,
            )
            po.ovarrf(
                self.outfile,
                "Outlet temperature of blanket coolant (K)",
                "(outlet_temp)",
                fwbs_variables.outlet_temp,
            )
            po.ovarre(
                self.outfile,
                "Blanket (primary) coolant pressure (Pa)",
                "(blpressure)",
                fwbs_variables.blpressure,
            )
            if fwbs_variables.ipump == 1:
                po.ovarre(
                    self.outfile,
                    "Blanket coolant mass flow rate (kg/s)",
                    "(mfblkt)",
                    fwbs_variables.mfblkt,
                    "OP ",
                )

            # Total primary coolant mass flow rate (if they are the same coolant)
            if fwbs_variables.ipump == 0:
                po.ovarre(
                    self.outfile,
                    "Total (FW+BB) primary coolant mass flow rate(kg/s)",
                    "(mftotal)",
                    blanket_library.mftotal,
                    "OP ",
                )

            # BB Liquid Metal Breeder !!!!!!!
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

            # Pumping Power !!!!!!!!!!!!!!!!!!!!!!!!
            po.osubhd(self.outfile, "Mechanical pumping power: ")

            if fwbs_variables.ipump == 1:
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
        frict_drop = blanket_library.pressure_drop(
            int(output),
            self.outfile,
            icoolpump,
            no90,
            no180,
            flleng,
            coolant_density,
            coolant_dynamic_viscosity,
            flow_velocity,
            label,
        )

        if icoolpump == 2:
            mhd_drop = blanket_library.liquid_breeder_pressure_drop_mhd(
                int(output),
                self.outfile,
                flow_velocity,
                coolant_dynamic_viscosity,
                coolant_electrical_conductivity,
                pol_channel_length,
                nopolchan,
                label,
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
        # Pumping power !!!!!!!!!!!!!!!!!

        # Outlet pressure is 'pressure'
        # Inlet pressure (Pa)
        coolpin = pressure + pdrop

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
                * (1 - (coolpin / pressure) ** -0.4)
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
                * (1 - (coolpin / pressure) ** -0.4)
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
