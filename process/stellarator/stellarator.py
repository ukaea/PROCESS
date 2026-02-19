import logging
from pathlib import Path

import numpy as np

import process.fusion_reactions as reactions
import process.physics_functions as physics_funcs
from process import constants
from process import process_output as po
from process.coolprop_interface import FluidProperties
from process.data_structure import (
    build_variables,
    constraint_variables,
    cost_variables,
    current_drive_variables,
    divertor_variables,
    first_wall_variables,
    fwbs_variables,
    global_variables,
    heat_transport_variables,
    numerics,
    physics_variables,
    stellarator_configuration,
    stellarator_variables,
    structure_variables,
    tfcoil_variables,
)
from process.exceptions import ProcessValueError
from process.physics import Physics, rether
from process.stellarator.build import st_build
from process.stellarator.coils.calculate import st_coil
from process.stellarator.denisty_limits import power_at_ignition_point, st_denisty_limits
from process.stellarator.divertor import st_div
from process.stellarator.heating import st_heat
from process.stellarator.preset_config import load_stellarator_config

logger = logging.getLogger(__name__)

# NOTE: a different value of electron_charge was used in the original implementation
# making the post-Python results slightly different. As a result, there is a
# relative tolerance on the neoclassics tests of 1e-3
KEV = 1e3 * constants.ELECTRON_CHARGE  # Kiloelectron-volt (keV)


class Stellarator:
    """Module containing stellarator routines
    author: P J Knight, CCFE, Culham Science Centre
    N/A
    This module contains routines for calculating the
    parameters of the first wall, blanket and shield components
    of a fusion power plant.

    """

    def __init__(
        self,
        availability,
        vacuum,
        buildings,
        costs,
        power,
        plasma_profile,
        hcpb,
        current_drive,
        physics: Physics,
        neoclassics,
        plasma_beta,
    ) -> None:
        """Initialises the Stellarator model's variables

        :param availability: a pointer to the availability model, allowing use of availability's variables/methods
        :type availability: process.availability.Availability
        :param buildings: a pointer to the buildings model, allowing use of buildings's variables/methods
        :type buildings: process.buildings.Buildings
        :param Vacuum: a pointer to the vacuum model, allowing use of vacuum's variables/methods
        :type Vacuum: process.vacuum.Vacuum
        :param costs: a pointer to the costs model, allowing use of costs' variables/methods
        :type costs: process.costs.Costs
        :param plasma_profile: a pointer to the plasma_profile model, allowing use of plasma_profile's variables/methods
        :type plasma_profile: process.plasma_profile.PlasmaProfile
        :param hcpb: a pointer to the ccfe_hcpb model, allowing use of ccfe_hcpb's variables/methods
        :type hcpb: process.hcpb.CCFE_HCPB
        :param current_drive: a pointer to the CurrentDrive model, allowing use of CurrentDrives's variables/methods
        :type current_drive: process.current_drive.CurrentDrive
        :param physics: a pointer to the Physics model, allowing use of Physics's variables/methods
        :type physics: process.physics.Physics
        :param neoclassics: a pointer to the Neoclassics model, allowing use of neoclassics's variables/methods
        :type neoclassics: process.stellarator.Neoclassics
        """

        self.outfile: int = constants.NOUT
        self.first_call_stfwbs = True

        self.availability = availability
        self.buildings = buildings
        self.vacuum = vacuum
        self.costs = costs
        self.power = power
        self.plasma_profile = plasma_profile
        self.hcpb = hcpb
        self.current_drive = current_drive
        self.physics = physics
        self.neoclassics = neoclassics
        self.beta = plasma_beta

    def run(self, output: bool):
        """Routine to call the physics and engineering modules
        relevant to stellarators
        author: P J Knight, CCFE, Culham Science Centre
        author: F Warmer, IPP Greifswald

        This routine is the caller for the stellarator models.

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """

        if output:
            self.costs.run()
            self.costs.output()
            self.availability.run(output=True)
            self.physics.calculate_effective_charge_ionisation_profiles()
            self.physics.outplas()
            st_heat(self, True)
            self.st_phys(True)
            st_denisty_limits(self, True)

            # Change in density limit can result in changed dene?
            # A second call of st_phys is used to make sure it is consitent.
            # st_phys and denisty limits should be integarted to avoid this double call.
            # Problem was probably bigger in the older version

            self.st_phys(False)

            st_div(self, True)
            st_build(self, True)
            st_coil(self, True)
            self.st_strc(True)
            self.st_fwbs(True)

            self.power.tfpwr(output=True)
            self.buildings.run(output=True)
            self.vacuum.run(output=True)
            self.power.acpow(output=True)
            self.power.output_plant_electric_powers()

            return

        self.st_new_config()
        self.st_geom()
        self.st_phys(False)
        st_denisty_limits(self, False)
        st_coil(self, False)
        st_build(self, False)
        self.st_strc(False)
        self.st_fwbs(False)
        st_div(self, False)

        self.power.tfpwr(output=False)
        self.power.component_thermal_powers()
        self.power.calculate_cryo_loads()
        self.buildings.run(output=False)
        self.vacuum.run(output=False)
        self.power.acpow(output=False)
        self.power.plant_electric_production()
        # TODO: should availability.run be called
        # rather than availability.avail?
        self.availability.avail(output=False)
        self.costs.run()

        if 91 in numerics.icc:
            # This call is comparably time consuming..
            # If the respective constraint equation is not called, do not set the values
            (
                stellarator_variables.powerht_constraint,
                stellarator_variables.powerscaling_constraint,
            ) = power_at_ignition_point(
                stellarator_variables.max_gyrotron_frequency,
                stellarator_variables.te0_ecrh_achievable,
            )

        stellarator_variables.first_call = False

    def st_new_config(self):
        """author: J Lion, IPP Greifswald
        Routine to initialise the stellarator configuration

        Routine to initialise the stellarator configuration.
        This routine is called right before the calculation and could
        in principle overwrite variables from the input file.
        It overwrites rminor with rmajor and aspect ratio e.g.

        To clarify the coils scaling factor:
        Coil aspect ratio factor can be described with the reversed equation (so if we would know r_coil_minor)
        stellarator_variables.f_coil_aspect = (
            (physics_variables.rmajor / stellarator_variables.r_coil_minor) /
            (stellarator_configuration.stella_config_rmajor_ref /
             stellarator_configuration.stella_config_coil_rminor)
        )
        """

        load_stellarator_config(
            stellarator_variables.istell,
            Path(f"{global_variables.output_prefix}stella_conf.json"),
        )

        # If physics_variables.aspect ratio is not in numerics.ixc set it to default value
        # Or when you call it the first time
        if 1 not in numerics.ixc:
            physics_variables.aspect = stellarator_configuration.stella_config_aspect_ref

        # Set the physics_variables.rminor radius as result here.
        physics_variables.rminor = physics_variables.rmajor / physics_variables.aspect
        physics_variables.eps = 1.0e0 / physics_variables.aspect

        tfcoil_variables.n_tf_coils = (
            stellarator_configuration.stella_config_coilspermodule
            * stellarator_configuration.stella_config_symmetry
        )  # This overwrites tfcoil_variables.n_tf_coils in input file.

        stellarator_variables.f_st_rmajor = (
            physics_variables.rmajor / stellarator_configuration.stella_config_rmajor_ref
        )  # Size scaling factor with respect to the reference calculation
        stellarator_variables.f_st_rminor = (
            physics_variables.rminor / stellarator_configuration.stella_config_rminor_ref
        )  # Size scaling factor with respect to the reference calculation

        stellarator_variables.f_st_aspect = (
            physics_variables.aspect / stellarator_configuration.stella_config_aspect_ref
        )
        stellarator_variables.f_st_n_coils = tfcoil_variables.n_tf_coils / (
            stellarator_configuration.stella_config_coilspermodule
            * stellarator_configuration.stella_config_symmetry
        )  # Coil number factor
        stellarator_variables.f_st_b = (
            physics_variables.b_plasma_toroidal_on_axis
            / stellarator_configuration.stella_config_bt_ref
        )  # B-field scaling factor

        # Coil aspect ratio factor to the reference calculation (we use it to scale the coil minor radius)
        stellarator_variables.f_coil_aspect = stellarator_variables.f_st_coil_aspect

        # Coil major radius, scaled with respect to the reference calculation
        stellarator_variables.r_coil_major = (
            stellarator_configuration.stella_config_coil_rmajor
            * stellarator_variables.f_st_rmajor
        )
        # Coil minor radius, scaled with respect to the reference calculation
        stellarator_variables.r_coil_minor = (
            stellarator_configuration.stella_config_coil_rminor
            * stellarator_variables.f_st_rmajor
            / stellarator_variables.f_coil_aspect
        )

        stellarator_variables.f_coil_shape = (
            stellarator_configuration.stella_config_min_plasma_coil_distance
            + stellarator_configuration.stella_config_rminor_ref
        ) / stellarator_configuration.stella_config_coil_rminor

    def st_geom(self):
        """
                author: J Lion, IPP Greifswald
        Routine to calculate the plasma volume and surface area for
        a stellarator using precalculated effective values

        This routine calculates the plasma volume and surface area for
        a stellarator configuration.
        It is simple scaling based on a Fourier representation based on
        that described in Geiger documentation.

        J. Geiger, IPP Greifswald internal document:  'Darstellung von
        ineinandergeschachtelten toroidal geschlossenen Flaechen mit
        Fourierkoeffizienten' ('Representation of nested, closed
        surfaces with Fourier coefficients')
        """
        physics_variables.vol_plasma = (
            stellarator_variables.f_st_rmajor
            * stellarator_variables.f_st_rminor**2
            * stellarator_configuration.stella_config_vol_plasma
        )

        # Plasma surface scaled from effective parameter:
        physics_variables.a_plasma_surface = (
            stellarator_variables.f_st_rmajor
            * stellarator_variables.f_st_rminor
            * stellarator_configuration.stella_config_plasma_surface
        )

        # Plasma cross section area. Approximated
        physics_variables.a_plasma_poloidal = (
            np.pi * physics_variables.rminor * physics_variables.rminor
        )  # average, could be calculated for every toroidal angle if desired

        #  physics_variables.a_plasma_surface_outboard is retained only for obsolescent fispact calculation...

        #  Cross-sectional area, averaged over toroidal angle
        physics_variables.a_plasma_surface_outboard = (
            0.5e0 * physics_variables.a_plasma_surface
        )  # Used only in the divertor model; approximate as for tokamaks

    def st_strc(self, output):
        """
                Routine to calculate the structural masses for a stellarator
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This routine calculates the structural masses for a stellarator.
        This is the stellarator version of routine
        <A HREF="struct.html">STRUCT</A>. In practice, many of the masses
        are simply set to zero to avoid double-counting of structural
        components that are specified differently for tokamaks.
        """
        structure_variables.fncmass = 0.0e0

        #  Reactor core gravity support mass
        structure_variables.gsmass = 0.0e0  # ? Not sure about this.

        # This is the previous scaling law for intercoil structure
        # We keep is here as a reference to the new model, which
        # we do not really trust yet.
        #  Mass of support structure (includes casing) (tonnes)
        #  Scaling for required structure mass (Steel) from:
        #  F.C. Moon, J. Appl. Phys. 53(12) (1982) 9112
        #
        #  Values based on regression analysis by Greifswald, March 2014
        m_struc = (
            1.3483e0
            * (1000.0e0 * tfcoil_variables.e_tf_magnetic_stored_total_gj) ** 0.7821e0
        )
        msupstr = 1000.0e0 * m_struc  # kg

        ################################################################
        # Intercoil support structure calculation:
        # Calculate the intercoil bolted plates structure from the coil surface

        intercoil_surface = (
            stellarator_configuration.stella_config_coilsurface
            * stellarator_variables.f_st_rmajor
            * (
                stellarator_variables.r_coil_minor
                / stellarator_configuration.stella_config_coil_rminor
            )
            - tfcoil_variables.dx_tf_inboard_out_toroidal
            * tfcoil_variables.len_tf_coil
            * tfcoil_variables.n_tf_coils
        )

        # This 0.18 m is an effective thickness which is scaled with empirial 1.5 law. 5.6 T is reference point of Helias
        # The thickness 0.18m was obtained as a measured value from Schauer, F. and Bykov, V. design of Helias 5-B. (Nucl Fus. 2013)
        structure_variables.aintmass = (
            0.18e0
            * (physics_variables.b_plasma_toroidal_on_axis / 5.6) ** 2
            * intercoil_surface
            * fwbs_variables.den_steel
        )

        structure_variables.clgsmass = (
            0.2e0 * structure_variables.aintmass
        )  # Very simple approximation for the gravity support.
        # This fits for the Helias 5b reactor design point ( F. and Bykov, V. design of Helias 5-B. (nucl Fus. 2013)).

        #  Total mass of cooled components
        structure_variables.coldmass = (
            tfcoil_variables.m_tf_coils_total
            + structure_variables.aintmass
            + fwbs_variables.dewmkg
        )

        #  Output section

        if output:
            po.oheadr(self.outfile, "Support Structure")
            po.ovarre(
                self.outfile,
                "Intercoil support structure mass (from intercoil calculation) (kg)",
                "(aintmass)",
                structure_variables.aintmass,
            )
            po.ovarre(
                self.outfile,
                "Intercoil support structure mass (scaling, for comparison) (kg)",
                "(empiricalmass)",
                msupstr,
            )
            po.ovarre(
                self.outfile,
                "Gravity support structure mass (kg)",
                "(clgsmass)",
                structure_variables.clgsmass,
            )
            po.ovarre(
                self.outfile,
                "Mass of cooled components (kg)",
                "(coldmass)",
                structure_variables.coldmass,
            )

    def blanket_neutronics(self):
        """Routine to calculate neutronic properties for a stellarator"""

        # heating of the blanket
        if fwbs_variables.breedmat == 1:
            fwbs_variables.breeder = "Orthosilicate"
            fwbs_variables.densbreed = 1.50e3
        elif fwbs_variables.breedmat == 2:
            fwbs_variables.breeder = "Metatitanate"
            fwbs_variables.densbreed = 1.78e3
        else:
            fwbs_variables.breeder = (
                "Zirconate"  # (In reality, rarely used - activation problems)
            )
            fwbs_variables.densbreed = 2.12e3

        fwbs_variables.m_blkt_total = (
            fwbs_variables.vol_blkt_total * fwbs_variables.densbreed
        )
        self.hcpb.nuclear_heating_blanket()

        # Heating of the magnets
        self.hcpb.nuclear_heating_magnets(False)

        # Rough estimate of TF coil volume used, assuming 25% of the total
        # TF coil perimeter is inboard, 75% outboard
        tf_volume = (
            0.25 * tfcoil_variables.len_tf_coil * tfcoil_variables.a_tf_inboard_total
            + 0.75
            * tfcoil_variables.len_tf_coil
            * tfcoil_variables.a_tf_leg_outboard
            * tfcoil_variables.n_tf_coils
        )

        fwbs_variables.ptfnucpm3 = fwbs_variables.p_tf_nuclear_heat_mw / tf_volume

        # heating of the shield
        self.hcpb.nuclear_heating_shield()

        # Energy multiplication factor
        fwbs_variables.f_p_blkt_multiplication = 1.269

        # Use older model to calculate neutron fluence since it
        # is not calculated in the CCFE blanket model
        (
            _,
            _,
            _,
            fwbs_variables.nflutf,
            _,
            _,
            _,
            _,
            _,
            _,
        ) = self.sc_tf_coil_nuclear_heating_iter90()

        # blktlife calculation left entierly to availability
        # Cannot find calculation for vvhemax in CCFE blanket

    def st_fwbs(self, output: bool):
        """Routine to calculate first wall, blanket and shield properties
        for a stellarator
        author: P J Knight, CCFE, Culham Science Centre
        author: F Warmer, IPP Greifswald
        outfile : input integer : Fortran output unit identifier
        iprint : input integer : Switch to write output to file (1=yes)
        This routine calculates a stellarator's first wall, blanket and
        shield properties.
        It calculates the nuclear heating in the blanket / shield, and
        estimates the volume and masses of the first wall,
        blanket, shield and vacuum vessel.
        <P>The arrays <CODE>coef(i,j)</CODE> and <CODE>decay(i,j)</CODE>
        are used for exponential decay approximations of the
        (superconducting) TF coil nuclear parameters.
        <UL><P><LI><CODE>j = 1</CODE> : stainless steel shield (assumed)
        <P><LI><CODE>j = 2</CODE> : tungsten shield (not used)</UL>
        Note: Costing and mass calculations elsewhere assume
        stainless steel only.
        <P>The method is the same as for tokamaks (as performed via
        <A HREF="fwbs.html">fwbs</A>), except for the volume calculations,
        which scale the surface area of the components from that
        of the plasma.
        """
        fwbs_variables.life_fw_fpy = min(
            cost_variables.abktflnc / physics_variables.pflux_fw_neutron_mw,
            cost_variables.life_plant,
        )

        #  First wall inboard, outboard areas (assume 50% of total each)
        first_wall_variables.a_fw_inboard = 0.5e0 * first_wall_variables.a_fw_total
        first_wall_variables.a_fw_outboard = 0.5e0 * first_wall_variables.a_fw_total

        #  Blanket volume; assume that its surface area is scaled directly from the
        #  plasma surface area.
        #  Uses fwbs_variables.fhole etc. to take account of gaps due to ports etc.

        r1 = physics_variables.rminor + 0.5e0 * (
            build_variables.dr_fw_plasma_gap_inboard
            + build_variables.dr_fw_inboard
            + build_variables.dr_fw_plasma_gap_outboard
            + build_variables.dr_fw_outboard
        )
        if heat_transport_variables.ipowerflow == 0:
            build_variables.a_blkt_total_surface = (
                physics_variables.a_plasma_surface
                * r1
                / physics_variables.rminor
                * (1.0e0 - fwbs_variables.fhole)
            )
        else:
            build_variables.a_blkt_total_surface = (
                physics_variables.a_plasma_surface
                * r1
                / physics_variables.rminor
                * (
                    1.0e0
                    - fwbs_variables.fhole
                    - fwbs_variables.f_ster_div_single
                    - fwbs_variables.f_a_fw_outboard_hcd
                )
            )

        build_variables.a_blkt_inboard_surface = (
            0.5e0 * build_variables.a_blkt_total_surface
        )
        build_variables.a_blkt_outboard_surface = (
            0.5e0 * build_variables.a_blkt_total_surface
        )

        fwbs_variables.vol_blkt_inboard = (
            build_variables.a_blkt_inboard_surface * build_variables.dr_blkt_inboard
        )
        fwbs_variables.vol_blkt_outboard = (
            build_variables.a_blkt_outboard_surface * build_variables.dr_blkt_outboard
        )
        fwbs_variables.vol_blkt_total = (
            fwbs_variables.vol_blkt_inboard + fwbs_variables.vol_blkt_outboard
        )

        #  Shield volume
        #  Uses fvolsi, fwbs_variables.fvolso as area coverage factors

        r1 = r1 + 0.5e0 * (
            build_variables.dr_blkt_inboard + build_variables.dr_blkt_outboard
        )
        build_variables.a_shld_total_surface = (
            physics_variables.a_plasma_surface * r1 / physics_variables.rminor
        )
        build_variables.a_shld_inboard_surface = (
            0.5e0 * build_variables.a_shld_total_surface * fwbs_variables.fvolsi
        )
        build_variables.a_shld_outboard_surface = (
            0.5e0 * build_variables.a_shld_total_surface * fwbs_variables.fvolso
        )

        vol_shld_inboard = (
            build_variables.a_shld_inboard_surface * build_variables.dr_shld_inboard
        )
        vol_shld_outboard = (
            build_variables.a_shld_outboard_surface * build_variables.dr_shld_outboard
        )
        fwbs_variables.vol_shld_total = vol_shld_inboard + vol_shld_outboard

        #  Neutron power lost through holes in first wall (eventually absorbed by
        #  shield)

        fwbs_variables.pnucloss = (
            physics_variables.p_neutron_total_mw * fwbs_variables.fhole
        )

        # The peaking factor, obtained as precalculated parameter
        fwbs_variables.wallpf = (
            stellarator_configuration.stella_config_neutron_peakfactor
        )

        #  Blanket neutronics calculations
        if fwbs_variables.blktmodel == 1:
            self.blanket_neutronics()

            if heat_transport_variables.ipowerflow == 1:
                fwbs_variables.p_div_nuclear_heat_total_mw = (
                    physics_variables.p_neutron_total_mw
                    * fwbs_variables.f_ster_div_single
                )
                fwbs_variables.p_fw_hcd_nuclear_heat_mw = (
                    physics_variables.p_neutron_total_mw
                    * fwbs_variables.f_a_fw_outboard_hcd
                )
                fwbs_variables.p_fw_nuclear_heat_total_mw = (
                    physics_variables.p_neutron_total_mw
                    - fwbs_variables.p_div_nuclear_heat_total_mw
                    - fwbs_variables.pnucloss
                    - fwbs_variables.p_fw_hcd_nuclear_heat_mw
                )

                fwbs_variables.pradloss = (
                    physics_variables.p_plasma_rad_mw * fwbs_variables.fhole
                )
                fwbs_variables.p_div_rad_total_mw = (
                    physics_variables.p_plasma_rad_mw * fwbs_variables.f_ster_div_single
                )
                fwbs_variables.p_fw_hcd_rad_total_mw = (
                    physics_variables.p_plasma_rad_mw
                    * fwbs_variables.f_a_fw_outboard_hcd
                )
                fwbs_variables.p_fw_rad_total_mw = (
                    physics_variables.p_plasma_rad_mw
                    - fwbs_variables.p_div_rad_total_mw
                    - fwbs_variables.pradloss
                    - fwbs_variables.p_fw_hcd_rad_total_mw
                )

                heat_transport_variables.p_fw_coolant_pump_mw = (
                    heat_transport_variables.f_p_fw_coolant_pump_total_heat
                    * (
                        fwbs_variables.p_fw_nuclear_heat_total_mw
                        + fwbs_variables.p_fw_rad_total_mw
                        + current_drive_variables.p_beam_orbit_loss_mw
                    )
                )
                heat_transport_variables.p_blkt_coolant_pump_mw = (
                    heat_transport_variables.f_p_blkt_coolant_pump_total_heat
                    * fwbs_variables.p_blkt_nuclear_heat_total_mw
                )
                heat_transport_variables.p_shld_coolant_pump_mw = (
                    heat_transport_variables.f_p_shld_coolant_pump_total_heat
                    * fwbs_variables.p_shld_nuclear_heat_mw
                )
                heat_transport_variables.p_div_coolant_pump_mw = (
                    heat_transport_variables.f_p_div_coolant_pump_total_heat
                    * (
                        physics_variables.p_plasma_separatrix_mw
                        + fwbs_variables.p_div_nuclear_heat_total_mw
                        + fwbs_variables.p_div_rad_total_mw
                    )
                )

                #  Void fraction in first wall / breeding zone,
                #  for use in fwbs_variables.m_fw_total and coolvol calculation below

                f_a_fw_coolant_inboard = (
                    1.0e0
                    - fwbs_variables.fblbe
                    - fwbs_variables.fblbreed
                    - fwbs_variables.fblss
                )
                f_a_fw_coolant_outboard = f_a_fw_coolant_inboard

        else:
            fwbs_variables.pnuc_cp = 0.0e0

            if heat_transport_variables.ipowerflow == 0:
                #  Energy-multiplied neutron power

                pneut2 = (
                    physics_variables.p_neutron_total_mw
                    - fwbs_variables.pnucloss
                    - fwbs_variables.pnuc_cp
                ) * fwbs_variables.f_p_blkt_multiplication

                fwbs_variables.p_blkt_multiplication_mw = pneut2 - (
                    physics_variables.p_neutron_total_mw
                    - fwbs_variables.pnucloss
                    - fwbs_variables.pnuc_cp
                )

                #  Nuclear heating in the blanket

                decaybl = 0.075e0 / (
                    1.0e0
                    - fwbs_variables.f_a_blkt_cooling_channels
                    - fwbs_variables.fblli2o
                    - fwbs_variables.fblbe
                )

                fwbs_variables.p_blkt_nuclear_heat_total_mw = pneut2 * (
                    1.0e0 - np.exp(-build_variables.dr_blkt_outboard / decaybl)
                )

                #  Nuclear heating in the shield
                fwbs_variables.p_shld_nuclear_heat_mw = (
                    pneut2 - fwbs_variables.p_blkt_nuclear_heat_total_mw
                )

                #  Superconducting coil shielding calculations
                (
                    coilhtmx,
                    dpacop,
                    htheci,
                    fwbs_variables.nflutf,
                    pheci,
                    pheco,
                    ptfiwp,
                    ptfowp,
                    raddose,
                    fwbs_variables.p_tf_nuclear_heat_mw,
                ) = self.sc_tf_coil_nuclear_heating_iter90()

            else:  # heat_transport_variables.ipowerflow == 1
                #  Neutron power incident on divertor (MW)

                fwbs_variables.p_div_nuclear_heat_total_mw = (
                    physics_variables.p_neutron_total_mw
                    * fwbs_variables.f_ster_div_single
                )

                #  Neutron power incident on HCD apparatus (MW)

                fwbs_variables.p_fw_hcd_nuclear_heat_mw = (
                    physics_variables.p_neutron_total_mw
                    * fwbs_variables.f_a_fw_outboard_hcd
                )

                #  Neutron power deposited in first wall, blanket and shield (MW)

                pnucfwbs = (
                    physics_variables.p_neutron_total_mw
                    - fwbs_variables.p_div_nuclear_heat_total_mw
                    - fwbs_variables.pnucloss
                    - fwbs_variables.pnuc_cp
                    - fwbs_variables.p_fw_hcd_nuclear_heat_mw
                )

                #  Split between inboard and outboard by first wall area fractions

                pnucfwbsi = (
                    pnucfwbs
                    * first_wall_variables.a_fw_inboard
                    / first_wall_variables.a_fw_total
                )
                pnucfwbso = (
                    pnucfwbs
                    * first_wall_variables.a_fw_outboard
                    / first_wall_variables.a_fw_total
                )

                #  Radiation power incident on divertor (MW)

                fwbs_variables.p_fw_hcd_rad_total_mw = (
                    physics_variables.p_plasma_rad_mw
                    * fwbs_variables.f_a_fw_outboard_hcd
                )

                #  Radiation power incident on HCD apparatus (MW)

                fwbs_variables.p_fw_hcd_rad_total_mw = (
                    physics_variables.p_plasma_rad_mw
                    * fwbs_variables.f_a_fw_outboard_hcd
                )

                #  Radiation power lost through holes (eventually hits shield) (MW)

                fwbs_variables.pradloss = (
                    physics_variables.p_plasma_rad_mw * fwbs_variables.fhole
                )

                #  Radiation power incident on first wall (MW)

                fwbs_variables.p_fw_rad_total_mw = (
                    physics_variables.p_plasma_rad_mw
                    - fwbs_variables.p_div_rad_total_mw
                    - fwbs_variables.pradloss
                    - fwbs_variables.p_fw_hcd_rad_total_mw
                )

                #  Calculate the power deposited in the first wall, blanket and shield,
                #  and the required coolant pumping power

                #  If we have chosen pressurised water as the coolant, set the
                #  coolant outlet temperature as 20 deg C below the boiling point

                if fwbs_variables.i_blkt_coolant_type == 2:
                    if fwbs_variables.irefprop:
                        fwbs_variables.temp_blkt_coolant_out = (
                            FluidProperties.of(
                                "Water",
                                pressure=fwbs_variables.coolp,
                                vapor_quality=0,
                            )
                            - 20
                        )
                    else:
                        fwbs_variables.temp_blkt_coolant_out = (
                            273.15
                            + 168.396
                            + 0.314653 / fwbs_variables.coolp
                            + -0.000728 / fwbs_variables.coolp**2
                            + 31.588979 * np.log(fwbs_variables.coolp)
                            + 11.473141 * fwbs_variables.coolp
                            + -0.575335 * fwbs_variables.coolp**2
                            + 0.013165 * fwbs_variables.coolp**3
                        ) - 20

                bfwi = 0.5e0 * build_variables.dr_fw_inboard
                bfwo = 0.5e0 * build_variables.dr_fw_outboard

                f_a_fw_coolant_inboard = (
                    fwbs_variables.radius_fw_channel
                    * fwbs_variables.radius_fw_channel
                    / (bfwi * bfwi)
                )  # inboard FW coolant void fraction
                f_a_fw_coolant_outboard = (
                    fwbs_variables.radius_fw_channel
                    * fwbs_variables.radius_fw_channel
                    / (bfwo * bfwo)
                )  # outboard FW coolant void fraction

                #  First wall decay length (m) - improved calculation required

                decayfwi = fwbs_variables.declfw
                decayfwo = fwbs_variables.declfw

                #  Surface heat flux on first wall (MW) (sum = fwbs_variables.p_fw_rad_total_mw)

                psurffwi = (
                    fwbs_variables.p_fw_rad_total_mw
                    * first_wall_variables.a_fw_inboard
                    / first_wall_variables.a_fw_total
                )
                psurffwo = (
                    fwbs_variables.p_fw_rad_total_mw
                    * first_wall_variables.a_fw_outboard
                    / first_wall_variables.a_fw_total
                )

                #  Simple blanket model (fwbs_variables.i_p_coolant_pumping = 0 or 1) is assumed for stellarators

                #  The power deposited in the first wall, breeder zone and shield is
                #  calculated according to their dimensions and materials assuming
                #  an exponential attenuation of nuclear heating with increasing
                #  radial distance.  The pumping power for the coolant is calculated
                #  as a fraction of the total thermal power deposited in the
                #  coolant.

                p_fw_inboard_nuclear_heat_mw = pnucfwbsi * (
                    1.0e0 - np.exp(-2.0e0 * bfwi / decayfwi)
                )
                p_fw_outboard_nuclear_heat_mw = pnucfwbso * (
                    1.0e0 - np.exp(-2.0e0 * bfwo / decayfwo)
                )

                #  Neutron power reaching blanket and shield (MW)

                pnucbsi = pnucfwbsi - p_fw_inboard_nuclear_heat_mw
                pnucbso = pnucfwbso - p_fw_outboard_nuclear_heat_mw

                #  Blanket decay length (m) - improved calculation required

                decaybzi = fwbs_variables.declblkt
                decaybzo = fwbs_variables.declblkt

                #  Neutron power deposited in breeder zone (MW)

                pnucbzi = pnucbsi * (
                    1.0e0 - np.exp(-build_variables.dr_blkt_inboard / decaybzi)
                )
                pnucbzo = pnucbso * (
                    1.0e0 - np.exp(-build_variables.dr_blkt_outboard / decaybzo)
                )

                #  Calculate coolant pumping powers from input fraction.
                #  The pumping power is assumed to be a fraction, fpump, of the
                #  incident thermal power to each component so that
                #  htpmw_i = fpump_i*C, where C is the non-pumping thermal power
                #  deposited in the coolant

                #  First wall and Blanket pumping power (MW)

                if fwbs_variables.i_p_coolant_pumping == 0:
                    #    Use input
                    pass
                elif fwbs_variables.i_p_coolant_pumping == 1:
                    heat_transport_variables.p_fw_coolant_pump_mw = (
                        heat_transport_variables.f_p_fw_coolant_pump_total_heat
                        * (
                            p_fw_inboard_nuclear_heat_mw
                            + p_fw_outboard_nuclear_heat_mw
                            + psurffwi
                            + psurffwo
                            + current_drive_variables.p_beam_orbit_loss_mw
                        )
                    )
                    heat_transport_variables.p_blkt_coolant_pump_mw = (
                        heat_transport_variables.f_p_blkt_coolant_pump_total_heat
                        * (
                            pnucbzi * fwbs_variables.f_p_blkt_multiplication
                            + pnucbzo * fwbs_variables.f_p_blkt_multiplication
                        )
                    )
                else:
                    raise ProcessValueError(
                        "i_p_coolant_pumping = 0 or 1 only for stellarator"
                    )

                fwbs_variables.p_blkt_multiplication_mw = (
                    heat_transport_variables.f_p_blkt_coolant_pump_total_heat
                    * (pnucbzi * fwbs_variables.f_p_blkt_multiplication + pnucbzo)
                    * (fwbs_variables.f_p_blkt_multiplication - 1.0e0)
                )

                #  Total nuclear heating of first wall (MW)

                fwbs_variables.p_fw_nuclear_heat_total_mw = (
                    p_fw_inboard_nuclear_heat_mw + p_fw_outboard_nuclear_heat_mw
                )

                #  Total nuclear heating of blanket (MW)

                fwbs_variables.p_blkt_nuclear_heat_total_mw = (
                    pnucbzi + pnucbzo
                ) * fwbs_variables.f_p_blkt_multiplication

                fwbs_variables.p_blkt_multiplication_mw = (
                    fwbs_variables.p_blkt_multiplication_mw
                    + (pnucbzi + pnucbzo)
                    * (fwbs_variables.f_p_blkt_multiplication - 1.0e0)
                )

                #  Calculation of shield and divertor powers
                #  Shield and divertor powers and pumping powers are calculated using the same
                #  simplified method as the first wall and breeder zone when fwbs_variables.i_p_coolant_pumping = 1.
                #  i.e. the pumping power is a fraction of the total thermal power deposited in the
                #  coolant.

                #  Neutron power reaching the shield (MW)
                #  The power lost from the fwbs_variables.fhole area fraction is assumed to be incident upon the shield

                pnucsi = (
                    pnucbsi
                    - pnucbzi
                    + (fwbs_variables.pnucloss + fwbs_variables.pradloss)
                    * first_wall_variables.a_fw_inboard
                    / first_wall_variables.a_fw_total
                )
                pnucso = (
                    pnucbso
                    - pnucbzo
                    + (fwbs_variables.pnucloss + fwbs_variables.pradloss)
                    * first_wall_variables.a_fw_outboard
                    / first_wall_variables.a_fw_total
                )

                #  Improved calculation of shield power decay lengths required

                decayshldi = fwbs_variables.declshld
                decayshldo = fwbs_variables.declshld

                #  Neutron power deposited in the shield (MW)

                pnucshldi = pnucsi * (
                    1.0e0 - np.exp(-build_variables.dr_shld_inboard / decayshldi)
                )
                pnucshldo = pnucso * (
                    1.0e0 - np.exp(-build_variables.dr_shld_outboard / decayshldo)
                )

                fwbs_variables.p_shld_nuclear_heat_mw = pnucshldi + pnucshldo

                #  Calculate coolant pumping powers from input fraction.
                #  The pumping power is assumed to be a fraction, fpump, of the incident
                #  thermal power to each component so that,
                #     htpmw_i = fpump_i*C
                #  where C is the non-pumping thermal power deposited in the coolant

                if fwbs_variables.i_p_coolant_pumping == 1:
                    #  Shield pumping power (MW)
                    heat_transport_variables.p_shld_coolant_pump_mw = (
                        heat_transport_variables.f_p_shld_coolant_pump_total_heat
                        * (pnucshldi + pnucshldo)
                    )

                    #  Divertor pumping power (MW)
                    heat_transport_variables.p_div_coolant_pump_mw = (
                        heat_transport_variables.f_p_div_coolant_pump_total_heat
                        * (
                            physics_variables.p_plasma_separatrix_mw
                            + fwbs_variables.p_div_nuclear_heat_total_mw
                            + fwbs_variables.p_div_rad_total_mw
                        )
                    )

                #  Remaining neutron power to coils and else:where. This is assumed
                #  (for superconducting coils at least) to be absorbed by the
                #  coils, and so contributes to the cryogenic load

                if tfcoil_variables.i_tf_sup == 1:
                    fwbs_variables.p_tf_nuclear_heat_mw = (
                        pnucsi + pnucso - pnucshldi - pnucshldo
                    )
                else:  # resistive coils
                    fwbs_variables.p_tf_nuclear_heat_mw = 0.0e0

        #  Divertor mass
        #  N.B. divertor_variables.a_div_surface_total is calculated in stdiv after this point, so will
        #  be zero on first lap, hence the initial approximation

        if self.first_call_stfwbs:
            divertor_variables.a_div_surface_total = 50.0e0
            self.first_call_stfwbs = False

        divertor_variables.m_div_plate = (
            divertor_variables.a_div_surface_total
            * divertor_variables.den_div_structure
            * (1.0e0 - divertor_variables.f_vol_div_coolant)
            * divertor_variables.dx_div_plate
        )

        #  Start adding components of the coolant mass:
        #  Divertor coolant volume (m3)

        coolvol = (
            divertor_variables.a_div_surface_total
            * divertor_variables.f_vol_div_coolant
            * divertor_variables.dx_div_plate
        )

        #  Blanket mass, excluding coolant

        if fwbs_variables.blktmodel == 0:
            if (fwbs_variables.blkttype == 1) or (
                fwbs_variables.blkttype == 2
            ):  # liquid breeder (WCLL or HCLL)
                fwbs_variables.wtbllipb = (
                    fwbs_variables.vol_blkt_total * fwbs_variables.fbllipb * 9400.0e0
                )
                fwbs_variables.m_blkt_lithium = (
                    fwbs_variables.vol_blkt_total * fwbs_variables.fblli * 534.0e0
                )
                fwbs_variables.m_blkt_total = (
                    fwbs_variables.wtbllipb + fwbs_variables.m_blkt_lithium
                )
            else:  # solid breeder (HCPB); always for ipowerflow=0
                fwbs_variables.m_blkt_li2o = (
                    fwbs_variables.vol_blkt_total * fwbs_variables.fblli2o * 2010.0e0
                )
                fwbs_variables.m_blkt_beryllium = (
                    fwbs_variables.vol_blkt_total * fwbs_variables.fblbe * 1850.0e0
                )
                fwbs_variables.m_blkt_total = (
                    fwbs_variables.m_blkt_li2o + fwbs_variables.m_blkt_beryllium
                )

            fwbs_variables.m_blkt_steel_total = (
                fwbs_variables.vol_blkt_total
                * fwbs_variables.den_steel
                * fwbs_variables.fblss
            )
            fwbs_variables.m_blkt_vanadium = (
                fwbs_variables.vol_blkt_total * 5870.0e0 * fwbs_variables.fblvd
            )

            fwbs_variables.m_blkt_total = (
                fwbs_variables.m_blkt_total
                + fwbs_variables.m_blkt_steel_total
                + fwbs_variables.m_blkt_vanadium
            )

        else:  # volume fractions proportional to sub-assembly thicknesses
            fwbs_variables.m_blkt_steel_total = fwbs_variables.den_steel * (
                fwbs_variables.vol_blkt_inboard
                / build_variables.dr_blkt_inboard
                * (
                    build_variables.blbuith * fwbs_variables.fblss
                    + build_variables.blbmith * (1.0e0 - fwbs_variables.fblhebmi)
                    + build_variables.blbpith * (1.0e0 - fwbs_variables.fblhebpi)
                )
                + fwbs_variables.vol_blkt_outboard
                / build_variables.dr_blkt_outboard
                * (
                    build_variables.blbuoth * fwbs_variables.fblss
                    + build_variables.blbmoth * (1.0e0 - fwbs_variables.fblhebmo)
                    + build_variables.blbpoth * (1.0e0 - fwbs_variables.fblhebpo)
                )
            )
            fwbs_variables.m_blkt_beryllium = (
                1850.0e0
                * fwbs_variables.fblbe
                * (
                    (
                        fwbs_variables.vol_blkt_inboard
                        * build_variables.blbuith
                        / build_variables.dr_blkt_inboard
                    )
                    + (
                        fwbs_variables.vol_blkt_outboard
                        * build_variables.blbuoth
                        / build_variables.dr_blkt_outboard
                    )
                )
            )
            fwbs_variables.whtblbreed = (
                fwbs_variables.densbreed
                * fwbs_variables.fblbreed
                * (
                    (
                        fwbs_variables.vol_blkt_inboard
                        * build_variables.blbuith
                        / build_variables.dr_blkt_inboard
                    )
                    + (
                        fwbs_variables.vol_blkt_outboard
                        * build_variables.blbuoth
                        / build_variables.dr_blkt_outboard
                    )
                )
            )
            fwbs_variables.m_blkt_total = (
                fwbs_variables.m_blkt_steel_total
                + fwbs_variables.m_blkt_beryllium
                + fwbs_variables.whtblbreed
            )

            fwbs_variables.f_a_blkt_cooling_channels = (
                fwbs_variables.vol_blkt_inboard
                / fwbs_variables.vol_blkt_total
                * (  # inboard portion
                    (build_variables.blbuith / build_variables.dr_blkt_inboard)
                    * (
                        1.0e0
                        - fwbs_variables.fblbe
                        - fwbs_variables.fblbreed
                        - fwbs_variables.fblss
                    )
                    + (build_variables.blbmith / build_variables.dr_blkt_inboard)
                    * fwbs_variables.fblhebmi
                    + (build_variables.blbpith / build_variables.dr_blkt_inboard)
                    * fwbs_variables.fblhebpi
                )
            )
            fwbs_variables.f_a_blkt_cooling_channels = (
                fwbs_variables.f_a_blkt_cooling_channels
                + fwbs_variables.vol_blkt_outboard
                / fwbs_variables.vol_blkt_total
                * (  # outboard portion
                    (build_variables.blbuoth / build_variables.dr_blkt_outboard)
                    * (
                        1.0e0
                        - fwbs_variables.fblbe
                        - fwbs_variables.fblbreed
                        - fwbs_variables.fblss
                    )
                    + (build_variables.blbmoth / build_variables.dr_blkt_outboard)
                    * fwbs_variables.fblhebmo
                    + (build_variables.blbpoth / build_variables.dr_blkt_outboard)
                    * fwbs_variables.fblhebpo
                )
            )

        #  When fwbs_variables.blktmodel > 0, although the blanket is by definition helium-cooled
        #  in this case, the shield etc. are assumed to be water-cooled, and since
        #  water is heavier the calculation for fwbs_variables.m_fw_blkt_div_coolant_total is better done with
        #  i_blkt_coolant_type=2 if fwbs_variables.blktmodel > 0; thus we can ignore the helium coolant mass
        #  in the blanket.

        if fwbs_variables.blktmodel == 0:
            coolvol = (
                coolvol
                + fwbs_variables.vol_blkt_total
                * fwbs_variables.f_a_blkt_cooling_channels
            )

        # Shield mass
        fwbs_variables.whtshld = (
            fwbs_variables.vol_shld_total
            * fwbs_variables.den_steel
            * (1.0e0 - fwbs_variables.vfshld)
        )

        coolvol = coolvol + fwbs_variables.vol_shld_total * fwbs_variables.vfshld

        #  Penetration shield (set = internal shield)

        fwbs_variables.wpenshld = fwbs_variables.whtshld

        if heat_transport_variables.ipowerflow == 0:
            #  First wall mass
            #  (first wall area is calculated else:where)

            fwbs_variables.m_fw_total = (
                first_wall_variables.a_fw_total
                * (build_variables.dr_fw_inboard + build_variables.dr_fw_outboard)
                / 2.0e0
                * fwbs_variables.den_steel
                * (1.0e0 - fwbs_variables.fwclfr)
            )

            #  Surface areas adjacent to plasma

            coolvol = (
                coolvol
                + first_wall_variables.a_fw_total
                * (build_variables.dr_fw_inboard + build_variables.dr_fw_outboard)
                / 2.0e0
                * fwbs_variables.fwclfr
            )

        else:
            fwbs_variables.m_fw_total = fwbs_variables.den_steel * (
                first_wall_variables.a_fw_inboard
                * build_variables.dr_fw_inboard
                * (1.0e0 - f_a_fw_coolant_inboard)
                + first_wall_variables.a_fw_outboard
                * build_variables.dr_fw_outboard
                * (1.0e0 - f_a_fw_coolant_outboard)
            )
            coolvol = (
                coolvol
                + first_wall_variables.a_fw_inboard
                * build_variables.dr_fw_inboard
                * f_a_fw_coolant_inboard
                + first_wall_variables.a_fw_outboard
                * build_variables.dr_fw_outboard
                * f_a_fw_coolant_outboard
            )

            #  Average first wall coolant fraction, only used by old routines
            #  in fispact.f90, safety.f90

            fwbs_variables.fwclfr = (
                first_wall_variables.a_fw_inboard
                * build_variables.dr_fw_inboard
                * f_a_fw_coolant_inboard
                + first_wall_variables.a_fw_outboard
                * build_variables.dr_fw_outboard
                * f_a_fw_coolant_outboard
            ) / (
                first_wall_variables.a_fw_total
                * 0.5e0
                * (build_variables.dr_fw_inboard + build_variables.dr_fw_outboard)
            )

        #  Mass of coolant = volume * density at typical coolant
        #  temperatures and pressures
        #  N.B. for fwbs_variables.blktmodel > 0, mass of *water* coolant in the non-blanket
        #  structures is used (see comment above)

        if (fwbs_variables.blktmodel > 0) or (
            fwbs_variables.i_blkt_coolant_type == 2
        ):  # pressurised water coolant
            fwbs_variables.m_fw_blkt_div_coolant_total = coolvol * 806.719e0
        else:  # gaseous helium coolant
            fwbs_variables.m_fw_blkt_div_coolant_total = coolvol * 1.517e0

        #  Assume external cryostat is a torus with circular cross-section,
        #  centred on plasma major radius.
        #  N.B. No check made to see if coils etc. lie wholly within cryostat...

        #  External cryostat outboard major radius (m)

        fwbs_variables.r_cryostat_inboard = (
            build_variables.r_tf_outboard_mid
            + 0.5e0 * build_variables.dr_tf_outboard
            + fwbs_variables.dr_pf_cryostat
        )
        adewex = fwbs_variables.r_cryostat_inboard - physics_variables.rmajor

        #  External cryostat volume

        fwbs_variables.vol_cryostat = (
            4.0e0
            * (np.pi**2)
            * physics_variables.rmajor
            * adewex
            * build_variables.dr_cryostat
        )

        #  Internal vacuum vessel volume
        #  fwbs_variables.fvoldw accounts for ports, support, etc. additions

        r1 = physics_variables.rminor + 0.5e0 * (
            build_variables.dr_fw_plasma_gap_inboard
            + build_variables.dr_fw_inboard
            + build_variables.dr_blkt_inboard
            + build_variables.dr_shld_inboard
            + build_variables.dr_fw_plasma_gap_outboard
            + build_variables.dr_fw_outboard
            + build_variables.dr_blkt_outboard
            + build_variables.dr_shld_outboard
        )
        fwbs_variables.vol_vv = (
            (build_variables.dr_vv_inboard + build_variables.dr_vv_outboard)
            / 2.0e0
            * physics_variables.a_plasma_surface
            * r1
            / physics_variables.rminor
            * fwbs_variables.fvoldw
        )

        #  Vacuum vessel mass

        fwbs_variables.m_vv = fwbs_variables.vol_vv * fwbs_variables.den_steel

        #  Sum of internal vacuum vessel and external cryostat masses

        fwbs_variables.dewmkg = (
            fwbs_variables.vol_vv + fwbs_variables.vol_cryostat
        ) * fwbs_variables.den_steel

        if output:
            #  Output section

            po.oheadr(self.outfile, "First Wall / Blanket / Shield")
            po.ovarre(
                self.outfile,
                "Average neutron wall load (MW/m2)",
                "(pflux_fw_neutron_mw)",
                physics_variables.pflux_fw_neutron_mw,
            )
            if fwbs_variables.blktmodel > 0:
                po.ovarre(
                    self.outfile,
                    "Neutron wall load peaking factor",
                    "(wallpf)",
                    fwbs_variables.wallpf,
                )

            po.ovarre(
                self.outfile,
                "First wall full-power lifetime (years)",
                "(life_fw_fpy)",
                fwbs_variables.life_fw_fpy,
            )

            po.ovarre(
                self.outfile,
                "Inboard shield thickness (m)",
                "(dr_shld_inboard)",
                build_variables.dr_shld_inboard,
            )
            po.ovarre(
                self.outfile,
                "Outboard shield thickness (m)",
                "(dr_shld_outboard)",
                build_variables.dr_shld_outboard,
            )
            po.ovarre(
                self.outfile,
                "Top shield thickness (m)",
                "(dz_shld_upper)",
                build_variables.dz_shld_upper,
            )

            if fwbs_variables.blktmodel > 0:
                po.ovarre(
                    self.outfile,
                    "Inboard breeding zone thickness (m)",
                    "(blbuith)",
                    build_variables.blbuith,
                )
                po.ovarre(
                    self.outfile,
                    "Inboard box manifold thickness (m)",
                    "(blbmith)",
                    build_variables.blbmith,
                )
                po.ovarre(
                    self.outfile,
                    "Inboard back plate thickness (m)",
                    "(blbpith)",
                    build_variables.blbpith,
                )

            po.ovarre(
                self.outfile,
                "Inboard blanket thickness (m)",
                "(dr_blkt_inboard)",
                build_variables.dr_blkt_inboard,
            )
            if fwbs_variables.blktmodel > 0:
                po.ovarre(
                    self.outfile,
                    "Outboard breeding zone thickness (m)",
                    "(blbuoth)",
                    build_variables.blbuoth,
                )
                po.ovarre(
                    self.outfile,
                    "Outboard box manifold thickness (m)",
                    "(blbmoth)",
                    build_variables.blbmoth,
                )
                po.ovarre(
                    self.outfile,
                    "Outboard back plate thickness (m)",
                    "(blbpoth)",
                    build_variables.blbpoth,
                )

            po.ovarre(
                self.outfile,
                "Outboard blanket thickness (m)",
                "(dr_blkt_outboard)",
                build_variables.dr_blkt_outboard,
            )
            po.ovarre(
                self.outfile,
                "Top blanket thickness (m)",
                "(dz_blkt_upper)",
                build_variables.dz_blkt_upper,
            )

            if (heat_transport_variables.ipowerflow == 0) and (
                fwbs_variables.blktmodel == 0
            ):
                po.osubhd(self.outfile, "Coil nuclear parameters :")
                po.ovarre(
                    self.outfile, "Peak magnet heating (MW/m3)", "(coilhtmx)", coilhtmx
                )
                po.ovarre(
                    self.outfile,
                    "Inboard coil winding pack heating (MW)",
                    "(ptfiwp)",
                    ptfiwp,
                )
                po.ovarre(
                    self.outfile,
                    "Outboard coil winding pack heating (MW)",
                    "(ptfowp)",
                    ptfowp,
                )
                po.ovarre(
                    self.outfile, "Peak coil case heating (MW/m3)", "(htheci)", htheci
                )
                po.ovarre(
                    self.outfile, "Inboard coil case heating (MW)", "(pheci)", pheci
                )
                po.ovarre(
                    self.outfile, "Outboard coil case heating (MW)", "(pheco)", pheco
                )
                po.ovarre(self.outfile, "Insulator dose (rad)", "(raddose)", raddose)
                po.ovarre(
                    self.outfile,
                    "Maximum neutron fluence (n/m2)",
                    "(nflutf)",
                    fwbs_variables.nflutf,
                )
                po.ovarre(
                    self.outfile,
                    "Copper stabiliser displacements/atom",
                    "(dpacop)",
                    dpacop,
                )

            if fwbs_variables.blktmodel == 0:
                po.osubhd(self.outfile, "Nuclear heating :")
                po.ovarre(
                    self.outfile,
                    "Blanket heating (including energy multiplication) (MW)",
                    "(p_blkt_nuclear_heat_total_mw)",
                    fwbs_variables.p_blkt_nuclear_heat_total_mw,
                )
                po.ovarre(
                    self.outfile,
                    "Shield nuclear heating (MW)",
                    "(p_shld_nuclear_heat_mw)",
                    fwbs_variables.p_shld_nuclear_heat_mw,
                )
                po.ovarre(
                    self.outfile,
                    "Coil nuclear heating (MW)",
                    "(p_tf_nuclear_heat_mw)",
                    fwbs_variables.p_tf_nuclear_heat_mw,
                )
            else:
                po.osubhd(self.outfile, "Blanket neutronics :")
                po.ovarre(
                    self.outfile,
                    "Blanket heating (including energy multiplication) (MW)",
                    "(p_blkt_nuclear_heat_total_mw)",
                    fwbs_variables.p_blkt_nuclear_heat_total_mw,
                )
                po.ovarre(
                    self.outfile,
                    "Shield heating (MW)",
                    "(p_shld_nuclear_heat_mw)",
                    fwbs_variables.p_shld_nuclear_heat_mw,
                )
                po.ovarre(
                    self.outfile,
                    "Energy multiplication in blanket",
                    "(f_p_blkt_multiplication)",
                    fwbs_variables.f_p_blkt_multiplication,
                )
                po.ovarin(
                    self.outfile,
                    "Number of divertor ports assumed",
                    "(npdiv)",
                    fwbs_variables.npdiv,
                )
                po.ovarin(
                    self.outfile,
                    "Number of inboard H/CD ports assumed",
                    "(nphcdin)",
                    fwbs_variables.nphcdin,
                )
                po.ovarin(
                    self.outfile,
                    "Number of outboard H/CD ports assumed",
                    "(nphcdout)",
                    fwbs_variables.nphcdout,
                )
                if fwbs_variables.hcdportsize == 1:
                    po.ocmmnt(
                        self.outfile, "     (small heating/current drive ports assumed)"
                    )
                else:
                    po.ocmmnt(
                        self.outfile, "     (large heating/current drive ports assumed)"
                    )

                if fwbs_variables.breedmat == 1:
                    po.ocmmnt(
                        self.outfile,
                        "Breeder material: Lithium orthosilicate (Li4Si04)",
                    )
                elif fwbs_variables.breedmat == 2:
                    po.ocmmnt(
                        self.outfile,
                        "Breeder material: Lithium methatitanate (Li2TiO3)",
                    )
                elif fwbs_variables.breedmat == 3:
                    po.ocmmnt(
                        self.outfile, "Breeder material: Lithium zirconate (Li2ZrO3)"
                    )
                else:  # shouldn't get here...
                    po.ocmmnt(self.outfile, "Unknown breeder material...")

                po.ovarre(
                    self.outfile,
                    "Lithium-6 enrichment (%)",
                    "(f_blkt_li6_enrichment)",
                    fwbs_variables.f_blkt_li6_enrichment,
                )
                po.ovarre(
                    self.outfile,
                    "Tritium production rate (g/day)",
                    "(tritprate)",
                    fwbs_variables.tritprate,
                )
                po.ovarre(
                    self.outfile,
                    "Nuclear heating on i/b coil (MW/m3)",
                    "(pnuctfi)",
                    fwbs_variables.pnuctfi,
                )
                po.ovarre(
                    self.outfile,
                    "Nuclear heating on o/b coil (MW/m3)",
                    "(pnuctfo)",
                    fwbs_variables.pnuctfo,
                )
                po.ovarre(
                    self.outfile,
                    "Total nuclear heating on coil (MW)",
                    "(p_tf_nuclear_heat_mw)",
                    fwbs_variables.p_tf_nuclear_heat_mw,
                )
                po.ovarre(
                    self.outfile,
                    "Fast neut. fluence on i/b coil (n/m2)",
                    "(nflutfi)",
                    fwbs_variables.nflutfi * 1.0e4,
                )
                po.ovarre(
                    self.outfile,
                    "Fast neut. fluence on o/b coil (n/m2)",
                    "(nflutfo)",
                    fwbs_variables.nflutfo * 1.0e4,
                )
                po.ovarre(
                    self.outfile,
                    "Minimum final He conc. in IB VV (appm)",
                    "(vvhemini)",
                    fwbs_variables.vvhemini,
                )
                po.ovarre(
                    self.outfile,
                    "Minimum final He conc. in OB VV (appm)",
                    "(vvhemino)",
                    fwbs_variables.vvhemino,
                )
                po.ovarre(
                    self.outfile,
                    "Maximum final He conc. in IB VV (appm)",
                    "(vvhemaxi)",
                    fwbs_variables.vvhemaxi,
                )
                po.ovarre(
                    self.outfile,
                    "Maximum final He conc. in OB VV (appm)",
                    "(vvhemaxo)",
                    fwbs_variables.vvhemaxo,
                )
                po.ovarre(
                    self.outfile,
                    "Blanket lifetime (full power years)",
                    "(life_blkt_fpy)",
                    fwbs_variables.life_blkt_fpy,
                )
                po.ovarre(
                    self.outfile,
                    "Blanket lifetime (calendar years)",
                    "(t_bl_y)",
                    fwbs_variables.t_bl_y,
                )

            if (heat_transport_variables.ipowerflow == 1) and (
                fwbs_variables.blktmodel == 0
            ):
                po.oblnkl(self.outfile)
                po.ovarin(
                    self.outfile,
                    "First wall / blanket thermodynamic model",
                    "(i_thermal_electric_conversion)",
                    fwbs_variables.i_thermal_electric_conversion,
                )
                if fwbs_variables.i_thermal_electric_conversion == 0:
                    po.ocmmnt(self.outfile, "   (Simple calculation)")

            po.osubhd(self.outfile, "Blanket / shield volumes and weights :")

            po.osubhd(self.outfile, "Other volumes, masses and areas :")
            po.ovarre(
                self.outfile,
                "First wall area (m2)",
                "(a_fw_total)",
                first_wall_variables.a_fw_total,
            )
            po.ovarre(
                self.outfile,
                "First wall mass (kg)",
                "(m_fw_total)",
                fwbs_variables.m_fw_total,
            )
            po.ovarre(
                self.outfile,
                "External cryostat inner radius (m)",
                "",
                fwbs_variables.r_cryostat_inboard - 2.0e0 * adewex,
            )
            po.ovarre(
                self.outfile,
                "External cryostat outer radius (m)",
                "(r_cryostat_inboard)",
                fwbs_variables.r_cryostat_inboard,
            )
            po.ovarre(
                self.outfile, "External cryostat minor radius (m)", "(adewex)", adewex
            )
            po.ovarre(
                self.outfile,
                "External cryostat shell volume (m^3)",
                "(vol_cryostat)",
                fwbs_variables.vol_cryostat,
            )
            po.ovarre(
                self.outfile,
                "Internal volume of the cryostat structure (m^3)",
                "(vol_cryostat_internal)",
                fwbs_variables.vol_cryostat_internal,
            )
            po.ovarre(
                self.outfile,
                "External cryostat mass (kg)",
                "",
                fwbs_variables.dewmkg - fwbs_variables.m_vv,
            )
            po.ovarre(
                self.outfile,
                "Internal vacuum vessel shell volume (m3)",
                "(vol_vv)",
                fwbs_variables.vol_vv,
            )
            po.ovarre(
                self.outfile,
                "Vacuum vessel mass (kg)",
                "(m_vv)",
                fwbs_variables.m_vv,
            )
            po.ovarre(
                self.outfile,
                "Total cryostat + vacuum vessel mass (kg)",
                "(dewmkg)",
                fwbs_variables.dewmkg,
            )
            po.ovarre(
                self.outfile,
                "Divertor area (m2)",
                "(a_div_surface_total)",
                divertor_variables.a_div_surface_total,
            )
            po.ovarre(
                self.outfile,
                "Divertor mass (kg)",
                "(m_div_plate)",
                divertor_variables.m_div_plate,
            )

    def sc_tf_coil_nuclear_heating_iter90(self):
        """Superconducting TF coil nuclear heating estimate
        author: P J Knight, CCFE, Culham Science Centre
        coilhtmx : output real : peak magnet heating (MW/m3)
        dpacop : output real : copper stabiliser displacements/atom
        htheci : output real : peak TF coil case heating (MW/m3)
        nflutf : output real : maximum neutron fluence (n/m2)
        pheci : output real : inboard coil case heating (MW)
        pheco : output real : outboard coil case heating (MW)
        ptfiwp : output real : inboard TF coil winding pack heating (MW)
        ptfowp : output real : outboard TF coil winding pack heating (MW)
        raddose : output real : insulator dose (rad)
        p_tf_nuclear_heat_mw : output real : TF coil nuclear heating (MW)
        This subroutine calculates the nuclear heating in the
        superconducting TF coils, assuming an exponential neutron
        attenuation through the blanket and shield materials.
        The estimates are based on 1990 ITER data.
        <P>The arrays <CODE>coef(i,j)</CODE> and <CODE>decay(i,j)</CODE>
        are used for exponential decay approximations of the
        (superconducting) TF coil nuclear parameters.
        <UL><P><LI><CODE>j = 1</CODE> : stainless steel shield (assumed)
        <P><LI><CODE>j = 2</CODE> : tungsten shield (not used)</UL>
        Note: Costing and mass calculations elsewhere assume
        stainless steel only.
        """

        ishmat = 0  # stainless steel coil casing is assumed

        if tfcoil_variables.i_tf_sup != 1:  # Resistive coils
            coilhtmx = 0.0
            ptfiwp = 0.0
            ptfowp = 0.0
            htheci = 0.0
            pheci = 0.0
            pheco = 0.0
            raddose = 0.0
            nflutf = 0.0
            dpacop = 0.0
            p_tf_nuclear_heat_mw = 0.0

        else:
            # TF coil nuclear heating coefficients in region i (first element),
            # assuming shield material j (second element where present)

            fact = np.array([8.0, 8.0, 6.0, 4.0, 4.0])
            coef = np.array([
                [10.3, 11.6, 7.08e5, 2.19e18, 3.33e-7],
                [8.32, 10.6, 7.16e5, 2.39e18, 3.84e-7],
            ]).T

            decay = np.array([
                [10.05, 17.61, 13.82, 13.24, 14.31, 13.26, 13.25],
                [10.02, 3.33, 15.45, 14.47, 15.87, 15.25, 17.25],
            ]).T

            # N.B. The vacuum vessel appears to be ignored

            dshieq = (
                build_variables.dr_shld_inboard
                + build_variables.dr_fw_inboard
                + build_variables.dr_blkt_inboard
            )
            dshoeq = (
                build_variables.dr_shld_outboard
                + build_variables.dr_fw_outboard
                + build_variables.dr_blkt_outboard
            )

            # Winding pack radial thickness, including groundwall insulation

            wpthk = (
                tfcoil_variables.dr_tf_wp_with_insulation
                + 2.0 * tfcoil_variables.dx_tf_wp_insulation
            )

            # Nuclear heating rate in inboard TF coil (MW/m**3)

            coilhtmx = (
                fact[0]
                * physics_variables.pflux_fw_neutron_mw
                * coef[0, ishmat]
                * np.exp(
                    -decay[5, ishmat] * (dshieq + tfcoil_variables.dr_tf_plasma_case)
                )
            )

            # Total nuclear heating (MW)

            ptfiwp = (
                coilhtmx
                * tfcoil_variables.tfsai
                * (1.0 - np.exp(-decay[0, ishmat] * wpthk))
                / decay[0, ishmat]
            )
            ptfowp = (
                fact[0]
                * physics_variables.pflux_fw_neutron_mw
                * coef[0, ishmat]
                * np.exp(
                    -decay[5, ishmat] * (dshoeq + tfcoil_variables.dr_tf_plasma_case)
                )
                * tfcoil_variables.tfsao
                * (1.0 - np.exp(-decay[0, ishmat] * wpthk))
                / decay[0, ishmat]
            )

            # Nuclear heating in plasma-side TF coil case (MW)

            htheci = (
                fact[1]
                * physics_variables.pflux_fw_neutron_mw
                * coef[1, ishmat]
                * np.exp(-decay[6, ishmat] * dshieq)
            )
            pheci = (
                htheci
                * tfcoil_variables.tfsai
                * (1.0 - np.exp(-decay[1, ishmat] * tfcoil_variables.dr_tf_plasma_case))
                / decay[1, ishmat]
            )
            pheco = (
                fact[1]
                * physics_variables.pflux_fw_neutron_mw
                * coef[1, ishmat]
                * np.exp(-decay[6, ishmat] * dshoeq)
                * tfcoil_variables.tfsao
                * (1.0 - np.exp(-decay[1, ishmat] * tfcoil_variables.dr_tf_plasma_case))
                / decay[1, ishmat]
            )
            ptfi = ptfiwp + pheci
            ptfo = ptfowp + pheco

            p_tf_nuclear_heat_mw = ptfi + ptfo

            # Full power DT operation years for replacement of TF Coil
            # (or plant life)

            fpydt = cost_variables.f_t_plant_available * cost_variables.life_plant
            fpsdt = fpydt * 3.154e7  # seconds

            # Insulator dose (rad)

            raddose = (
                coef[2, ishmat]
                * fpsdt
                * fact[2]
                * physics_variables.pflux_fw_neutron_mw
                * np.exp(
                    -decay[2, ishmat] * (dshieq + tfcoil_variables.dr_tf_plasma_case)
                )
            )

            # Maximum neutron fluence in superconductor (n/m**2)

            nflutf = (
                fpsdt
                * fact[3]
                * physics_variables.pflux_fw_neutron_mw
                * coef[3, ishmat]
                * np.exp(
                    -decay[3, ishmat] * (dshieq + tfcoil_variables.dr_tf_plasma_case)
                )
            )

            # Atomic displacement in copper stabilizer

            dpacop = (
                fpsdt
                * fact[4]
                * physics_variables.pflux_fw_neutron_mw
                * coef[4, ishmat]
                * np.exp(
                    -decay[4, ishmat] * (dshieq + tfcoil_variables.dr_tf_plasma_case)
                )
            )

        return (
            coilhtmx,
            dpacop,
            htheci,
            nflutf,
            pheci,
            pheco,
            ptfiwp,
            ptfowp,
            raddose,
            p_tf_nuclear_heat_mw,
        )

    def st_phys(self, output):
        """Routine to calculate stellarator plasma physics information
        author: P J Knight, CCFE, Culham Science Centre
        author: F Warmer, IPP Greifswald
        None
        This routine calculates the physics quantities relevant to
        a stellarator device.
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        # ###############################################
        #  Calculate plasma composition
        # Issue #261 Remove old radiation model

        self.physics.plasma_composition()

        # Calculate density and temperature profile quantities
        self.plasma_profile.run()

        #  Total field
        physics_variables.b_plasma_total = np.sqrt(
            physics_variables.b_plasma_toroidal_on_axis**2
            + physics_variables.b_plasma_poloidal_average**2
        )

        # Check if physics_variables.beta (iteration variable 5) is an iteration variable
        if 5 in numerics.ixc:
            raise ProcessValueError(
                "Beta should not be in ixc if istell>0. Use Constraints 24 and 84 instead"
            )

        #  Set physics_variables.beta as a consequence:
        #  This replaces constraint equation 1 as it is just an equality.
        physics_variables.beta_total_vol_avg = (
            physics_variables.beta_fast_alpha
            + physics_variables.beta_beam
            + 2.0e3
            * constants.RMU0
            * constants.ELECTRON_CHARGE
            * (
                physics_variables.nd_plasma_electrons_vol_avg
                * physics_variables.temp_plasma_electron_density_weighted_kev
                + physics_variables.nd_plasma_ions_total_vol_avg
                * physics_variables.temp_plasma_ion_density_weighted_kev
            )
            / physics_variables.b_plasma_total**2
        )
        physics_variables.e_plasma_beta = (
            1.5e0
            * physics_variables.beta_total_vol_avg
            * physics_variables.b_plasma_total
            * physics_variables.b_plasma_total
            / (2.0e0 * constants.RMU0)
            * physics_variables.vol_plasma
        )

        physics_variables.rho_star = np.sqrt(
            2.0e0
            * constants.PROTON_MASS
            * physics_variables.m_ions_total_amu
            * physics_variables.e_plasma_beta
            / (
                3.0e0
                * physics_variables.vol_plasma
                * physics_variables.nd_plasma_electron_line
            )
        ) / (
            constants.ELECTRON_CHARGE
            * physics_variables.b_plasma_toroidal_on_axis
            * physics_variables.eps
            * physics_variables.rmajor
        )

        #  Calculate poloidal field using rotation transform
        physics_variables.b_plasma_poloidal_average = (
            physics_variables.rminor
            * physics_variables.b_plasma_toroidal_on_axis
            / physics_variables.rmajor
            * stellarator_variables.iotabar
        )

        #  Perform auxiliary power calculations

        st_heat(self, False)

        #  Calculate fusion power

        fusion_reactions = reactions.FusionReactionRate(self.plasma_profile)
        fusion_reactions.deuterium_branching(
            physics_variables.temp_plasma_ion_vol_avg_kev
        )
        fusion_reactions.calculate_fusion_rates()
        fusion_reactions.set_physics_variables()

        # D-T power density is named differently to differentiate it from the beam given component
        physics_variables.p_plasma_dt_mw = (
            physics_variables.dt_power_density_plasma * physics_variables.vol_plasma
        )
        physics_variables.p_dhe3_total_mw = (
            physics_variables.dhe3_power_density * physics_variables.vol_plasma
        )
        physics_variables.p_dd_total_mw = (
            physics_variables.dd_power_density * physics_variables.vol_plasma
        )

        #  Calculate neutral beam slowing down effects
        #  If ignited, then ignore beam fusion effects

        if (current_drive_variables.p_hcd_beam_injected_total_mw != 0.0e0) and (
            physics_variables.i_plasma_ignited == 0
        ):
            (
                physics_variables.beta_beam,
                physics_variables.nd_beam_ions_out,
                physics_variables.p_beam_alpha_mw,
            ) = reactions.beam_fusion(
                physics_variables.beamfus0,
                physics_variables.betbm0,
                physics_variables.b_plasma_poloidal_average,
                physics_variables.b_plasma_toroidal_on_axis,
                current_drive_variables.c_beam_total,
                physics_variables.nd_plasma_electrons_vol_avg,
                physics_variables.nd_plasma_fuel_ions_vol_avg,
                physics_variables.dlamie,
                current_drive_variables.e_beam_kev,
                physics_variables.f_plasma_fuel_deuterium,
                physics_variables.f_plasma_fuel_tritium,
                current_drive_variables.f_beam_tritium,
                physics_variables.sigmav_dt_average,
                physics_variables.temp_plasma_electron_density_weighted_kev,
                physics_variables.temp_plasma_ion_density_weighted_kev,
                physics_variables.vol_plasma,
                physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
            )
            physics_variables.fusden_total = (
                physics_variables.fusden_plasma
                + 1.0e6
                * physics_variables.p_beam_alpha_mw
                / (constants.DT_ALPHA_ENERGY)
                / physics_variables.vol_plasma
            )
            physics_variables.fusden_alpha_total = (
                physics_variables.fusden_plasma_alpha
                + 1.0e6
                * physics_variables.p_beam_alpha_mw
                / (constants.DT_ALPHA_ENERGY)
                / physics_variables.vol_plasma
            )
            physics_variables.p_dt_total_mw = (
                physics_variables.p_plasma_dt_mw
                + 5.0e0 * physics_variables.p_beam_alpha_mw
            )
        else:
            # If no beams present then the total alpha rates and power are the same as the plasma values
            physics_variables.fusden_total = physics_variables.fusden_plasma
            physics_variables.fusden_alpha_total = physics_variables.fusden_plasma_alpha
            physics_variables.p_dt_total_mw = physics_variables.p_plasma_dt_mw

        # Create some derived values and add beam contribution to fusion power
        (
            physics_variables.pden_neutron_total_mw,
            physics_variables.p_plasma_alpha_mw,
            physics_variables.p_alpha_total_mw,
            physics_variables.p_plasma_neutron_mw,
            physics_variables.p_neutron_total_mw,
            physics_variables.p_non_alpha_charged_mw,
            physics_variables.pden_alpha_total_mw,
            physics_variables.f_pden_alpha_electron_mw,
            physics_variables.f_pden_alpha_ions_mw,
            physics_variables.p_charged_particle_mw,
            physics_variables.p_fusion_total_mw,
        ) = reactions.set_fusion_powers(
            physics_variables.f_alpha_electron,
            physics_variables.f_alpha_ion,
            physics_variables.p_beam_alpha_mw,
            physics_variables.pden_non_alpha_charged_mw,
            physics_variables.pden_plasma_neutron_mw,
            physics_variables.vol_plasma,
            physics_variables.pden_plasma_alpha_mw,
        )

        physics_variables.beta_fast_alpha = self.beta.fast_alpha_beta(
            physics_variables.b_plasma_poloidal_average,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.nd_plasma_fuel_ions_vol_avg,
            physics_variables.nd_plasma_ions_total_vol_avg,
            physics_variables.temp_plasma_electron_density_weighted_kev,
            physics_variables.temp_plasma_ion_density_weighted_kev,
            physics_variables.pden_alpha_total_mw,
            physics_variables.pden_plasma_alpha_mw,
            physics_variables.i_beta_fast_alpha,
        )

        #  Neutron wall load

        if physics_variables.i_pflux_fw_neutron == 1:
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.ffwal
                * physics_variables.p_neutron_total_mw
                / physics_variables.a_plasma_surface
            )
        else:
            if heat_transport_variables.ipowerflow == 0:
                physics_variables.pflux_fw_neutron_mw = (
                    (1.0e0 - fwbs_variables.fhole)
                    * physics_variables.p_neutron_total_mw
                    / first_wall_variables.a_fw_total
                )
            else:
                physics_variables.pflux_fw_neutron_mw = (
                    (
                        1.0e0
                        - fwbs_variables.fhole
                        - fwbs_variables.f_a_fw_outboard_hcd
                        - fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_neutron_total_mw
                    / first_wall_variables.a_fw_total
                )

        #  Calculate ion/electron equilibration power

        physics_variables.pden_ion_electron_equilibration_mw = rether(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.dlamie,
            physics_variables.temp_plasma_electron_vol_avg_kev,
            physics_variables.temp_plasma_ion_vol_avg_kev,
            physics_variables.n_charge_plasma_effective_mass_weighted_vol_avg,
        )

        #  Calculate radiation power
        radpwr_data = physics_funcs.calculate_radiation_powers(
            self.plasma_profile,
            physics_variables.nd_plasma_electron_on_axis,
            physics_variables.rminor,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.aspect,
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.tbeta,
            physics_variables.temp_plasma_electron_on_axis_kev,
            physics_variables.f_sync_reflect,
            physics_variables.rmajor,
            physics_variables.kappa,
            physics_variables.vol_plasma,
        )
        physics_variables.pden_plasma_sync_mw = radpwr_data.pden_plasma_sync_mw
        physics_variables.pden_plasma_core_rad_mw = radpwr_data.pden_plasma_core_rad_mw
        physics_variables.pden_plasma_outer_rad_mw = radpwr_data.pden_plasma_outer_rad_mw
        physics_variables.pden_plasma_rad_mw = radpwr_data.pden_plasma_rad_mw

        physics_variables.pden_plasma_core_rad_mw = max(
            physics_variables.pden_plasma_core_rad_mw, 0.0e0
        )
        physics_variables.pden_plasma_outer_rad_mw = max(
            physics_variables.pden_plasma_outer_rad_mw, 0.0e0
        )

        physics_variables.p_plasma_inner_rad_mw = (
            physics_variables.pden_plasma_core_rad_mw * physics_variables.vol_plasma
        )  # Should probably be vol_core
        physics_variables.p_plasma_outer_rad_mw = (
            physics_variables.pden_plasma_outer_rad_mw * physics_variables.vol_plasma
        )

        physics_variables.p_plasma_rad_mw = (
            physics_variables.pden_plasma_rad_mw * physics_variables.vol_plasma
        )

        #  Heating power to plasma (= Psol in divertor model)
        #  Ohmic power is zero in a stellarator
        #  physics_variables.p_plasma_rad_mw here is core + edge (no SOL)

        powht = (
            physics_variables.f_p_alpha_plasma_deposited
            * physics_variables.p_alpha_total_mw
            + physics_variables.p_non_alpha_charged_mw
            + physics_variables.p_plasma_ohmic_mw
            - physics_variables.pden_plasma_rad_mw * physics_variables.vol_plasma
        )
        powht = max(
            0.00001e0, powht
        )  # To avoid negative heating power. This line is VERY important

        if physics_variables.i_plasma_ignited == 0:
            powht = (
                powht + current_drive_variables.p_hcd_injected_total_mw
            )  # if not ignited add the auxiliary power

        # Here the implementation sometimes leaves the accessible regime when p_plasma_rad_mw> powht which is unphysical and
        # is not taken care of by the rad module. We restrict the radiation power here by the heating power:
        physics_variables.p_plasma_rad_mw = max(0.0e0, physics_variables.p_plasma_rad_mw)

        #  Power to divertor, = (1-stellarator_variables.f_rad)*Psol

        # The SOL radiation needs to be smaller than the physics_variables.p_plasma_rad_mw
        physics_variables.psolradmw = stellarator_variables.f_rad * powht
        physics_variables.p_plasma_separatrix_mw = powht - physics_variables.psolradmw

        # Add SOL Radiation to total
        physics_variables.p_plasma_rad_mw = (
            physics_variables.p_plasma_rad_mw + physics_variables.psolradmw
        )

        #  The following line is unphysical, but prevents -ve sqrt argument
        #  Should be obsolete if constraint eqn 17 is turned on (but beware -
        #  this may not be quite correct for stellarators)
        physics_variables.p_plasma_separatrix_mw = max(
            0.001e0, physics_variables.p_plasma_separatrix_mw
        )

        #  Power transported to the first wall by escaped alpha particles

        physics_variables.p_fw_alpha_mw = physics_variables.p_alpha_total_mw * (
            1.0e0 - physics_variables.f_p_alpha_plasma_deposited
        )

        # Nominal mean photon wall load
        if physics_variables.i_pflux_fw_neutron == 1:
            physics_variables.pflux_fw_rad_mw = (
                physics_variables.ffwal
                * physics_variables.p_plasma_rad_mw
                / physics_variables.a_plasma_surface
            )
        else:
            if heat_transport_variables.ipowerflow == 0:
                physics_variables.pflux_fw_rad_mw = (
                    (1.0e0 - fwbs_variables.fhole)
                    * physics_variables.p_plasma_rad_mw
                    / first_wall_variables.a_fw_total
                )
            else:
                physics_variables.pflux_fw_rad_mw = (
                    (
                        1.0e0
                        - fwbs_variables.fhole
                        - fwbs_variables.f_a_fw_outboard_hcd
                        - fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_plasma_rad_mw
                    / first_wall_variables.a_fw_total
                )

        constraint_variables.pflux_fw_rad_max_mw = (
            physics_variables.pflux_fw_rad_mw * constraint_variables.f_fw_rad_max
        )

        physics_variables.rad_fraction_total = physics_variables.p_plasma_rad_mw / (
            physics_variables.f_p_alpha_plasma_deposited
            * physics_variables.p_alpha_total_mw
            + physics_variables.p_non_alpha_charged_mw
            + physics_variables.p_plasma_ohmic_mw
            + current_drive_variables.p_hcd_injected_total_mw
        )

        #  Calculate transport losses and energy confinement time using the
        #  chosen scaling law
        #  N.B. stellarator_variables.iotabar replaces tokamak physics_variables.q95 in argument list

        (
            physics_variables.pden_electron_transport_loss_mw,
            physics_variables.pden_ion_transport_loss_mw,
            physics_variables.t_electron_energy_confinement,
            physics_variables.t_ion_energy_confinement,
            physics_variables.t_energy_confinement,
            physics_variables.p_plasma_loss_mw,
        ) = self.physics.calculate_confinement_time(
            physics_variables.m_fuel_amu,
            physics_variables.p_alpha_total_mw,
            physics_variables.aspect,
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.nd_plasma_ions_total_vol_avg,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.nd_plasma_electron_line,
            physics_variables.eps,
            physics_variables.hfact,
            physics_variables.i_confinement_time,
            physics_variables.i_plasma_ignited,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.p_non_alpha_charged_mw,
            current_drive_variables.p_hcd_injected_total_mw,
            physics_variables.plasma_current,
            physics_variables.pden_plasma_core_rad_mw,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.temp_plasma_electron_density_weighted_kev,
            physics_variables.temp_plasma_ion_density_weighted_kev,
            stellarator_variables.iotabar,
            physics_variables.qstar,
            physics_variables.vol_plasma,
            physics_variables.n_charge_plasma_effective_vol_avg,
        )

        physics_variables.p_electron_transport_loss_mw = (
            physics_variables.pden_electron_transport_loss_mw
            * physics_variables.vol_plasma
        )
        physics_variables.p_ion_transport_loss_mw = (
            physics_variables.pden_ion_transport_loss_mw * physics_variables.vol_plasma
        )

        physics_variables.pscalingmw = (
            physics_variables.p_electron_transport_loss_mw
            + physics_variables.p_ion_transport_loss_mw
        )

        #  Calculate auxiliary physics related information
        #  for the rest of the code

        sbar = 1.0e0
        (
            physics_variables.burnup,
            physics_variables.ntau,
            physics_variables.nTtau,
            physics_variables.figmer,
            _fusrat,
            physics_variables.molflow_plasma_fuelling_required,
            physics_variables.rndfuel,
            physics_variables.t_alpha_confinement,
            physics_variables.f_alpha_energy_confinement,
        ) = self.physics.phyaux(
            physics_variables.aspect,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.temp_plasma_electron_vol_avg_kev,
            physics_variables.nd_plasma_fuel_ions_vol_avg,
            physics_variables.fusden_total,
            physics_variables.fusden_alpha_total,
            physics_variables.plasma_current,
            sbar,
            physics_variables.nd_plasma_alphas_vol_avg,
            physics_variables.t_energy_confinement,
            physics_variables.vol_plasma,
        )

        # Calculate the neoclassical sanity check with PROCESS parameters
        (
            q_PROCESS,
            q_PROCESS_r1,
            _q_neo,
            _gamma_neo,
            _total_q_neo,
            total_q_neo_e,
            q_neo_e,
            _q_neo_D,
            _q_neo_a,
            _q_neo_T,
            g_neo_e,
            _g_neo_D,
            _g_neo_a,
            _g_neo_T,
            dndt_neo_e,
            _dndt_neo_D,
            _dndt_neo_a,
            _dndt_neo_T,
            _dndt_neo_fuel,
            _dmdt_neo_fuel,
            dmdt_neo_fuel_from_e,
            chi_neo_e,
            chi_PROCESS_e,
            nu_star_e,
            nu_star_d,
            nu_star_T,
            nu_star_He,
        ) = self.neoclassics.calc_neoclassics()

        if output:
            self.st_phys_output(
                q_PROCESS,
                total_q_neo_e,
                dmdt_neo_fuel_from_e,
                q_PROCESS_r1,
                chi_PROCESS_e,
                chi_neo_e,
                q_neo_e,
                g_neo_e,
                dndt_neo_e,
                physics_variables.rho_ne_max,
                physics_variables.rho_te_max,
                physics_variables.gradient_length_ne,
                physics_variables.gradient_length_te,
                physics_variables.rho_star,
                nu_star_e,
                nu_star_d,
                nu_star_T,
                nu_star_He,
                physics_variables.nd_plasma_electron_line,
                physics_variables.nd_plasma_electrons_max,
            )

    def st_phys_output(
        self,
        q_PROCESS,
        total_q_neo_e,
        dmdt_neo_fuel_from_e,
        q_PROCESS_r1,
        chi_PROCESS_e,
        chi_neo_e,
        q_neo_e,
        g_neo_e,
        dndt_neo_e,
        rho_ne_max,
        rho_te_max,
        gradient_length_ne,
        gradient_length_te,
        rho_star,
        nu_star_e,
        nu_star_D,
        nu_star_T,
        nu_star_He,
        nd_plasma_electron_line,
        nd_plasma_electrons_max,
    ):
        po.oheadr(self.outfile, "Stellarator Specific Physics:")

        po.ovarre(
            self.outfile,
            "Total 0D heat flux (r=rhocore) (MW/m2)",
            "(q_PROCESS)",
            q_PROCESS,
        )
        po.ovarre(
            self.outfile,
            "Total neoclassical flux from 4*q_e (r=rhocore) (MW/m2)",
            "(total_q_neo_e)",
            total_q_neo_e,
        )

        po.ovarre(
            self.outfile,
            "Total fuel (DT) mass flux by using 4 * neoclassical e transport (mg/s): ",
            "(dmdt_neo_fuel_from_e)",
            dmdt_neo_fuel_from_e,
        )
        po.ovarre(
            self.outfile,
            "Considered Heatflux by LCFS heat flux ratio (1)",
            "(q_PROCESS/q_PROCESS_r1)",
            q_PROCESS / q_PROCESS_r1,
        )

        po.ovarre(
            self.outfile,
            "Resulting electron effective chi (0D) (r=rhocore): ",
            "(chi_PROCESS_e)",
            chi_PROCESS_e,
        )
        po.ovarre(
            self.outfile,
            "Neoclassical electron effective chi (r=rhocore): ",
            "(chi_neo_e)",
            chi_neo_e,
        )

        po.ovarre(
            self.outfile,
            "Heat flux due to neoclassical energy transport (e) (MW/m2): ",
            "(q_neo_e)",
            q_neo_e,
        )
        po.ovarre(
            self.outfile,
            "Heat flux due to neoclassical particle transport (e) (MW/m2): ",
            "(g_neo_e)",
            g_neo_e,
        )
        po.ovarre(
            self.outfile,
            "Particle flux due to neoclassical particle transport (e) (1/m2/s): ",
            "(dndt_neo_e)",
            dndt_neo_e,
        )

        po.ovarre(
            self.outfile, "r/a of maximum ne gradient (m)", "(rho_ne_max)", rho_ne_max
        )
        po.ovarre(
            self.outfile, "r/a of maximum te gradient (m)", "(rho_te_max)", rho_te_max
        )
        po.ovarre(
            self.outfile,
            "Maxium ne gradient length (1)",
            "(gradient_length_ne)",
            gradient_length_ne,
        )
        po.ovarre(
            self.outfile,
            "Maxium te gradient length (1)",
            "(gradient_length_te)",
            gradient_length_te,
        )
        po.ovarre(
            self.outfile,
            "Gradient Length Ratio (T/n) (1)",
            "(gradient_length_ratio)",
            gradient_length_te / gradient_length_ne,
        )

        po.ovarre(self.outfile, "Normalized ion Larmor radius", "(rho_star)", rho_star)
        po.ovarre(
            self.outfile,
            "Normalized collisionality (electrons)",
            "(nu_star_e)",
            nu_star_e,
        )
        po.ovarre(
            self.outfile, "Normalized collisionality (D)", "(nu_star_D)", nu_star_D
        )
        po.ovarre(
            self.outfile, "Normalized collisionality (T)", "(nu_star_T)", nu_star_T
        )
        po.ovarre(
            self.outfile, "Normalized collisionality (He)", "(nu_star_He)", nu_star_He
        )

        po.ovarre(
            self.outfile,
            "Obtained line averaged density at op. point (/m3)",
            "(nd_plasma_electron_line)",
            nd_plasma_electron_line,
        )
        po.ovarre(
            self.outfile,
            "Sudo density limit (/m3)",
            "(nd_plasma_electrons_max)",
            nd_plasma_electrons_max,
        )
        po.ovarre(
            self.outfile,
            "Ratio density to sudo limit (1)",
            "(nd_plasma_electron_line/nd_plasma_electrons_max)",
            nd_plasma_electron_line / nd_plasma_electrons_max,
        )
