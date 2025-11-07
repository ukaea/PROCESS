import logging
from copy import copy
from pathlib import Path

import numpy as np

import process.fusion_reactions as reactions
import process.physics_functions as physics_funcs
import process.superconductors as superconductors
from process import constants
from process import process_output as po
from process.coolprop_interface import FluidProperties
from process.data_structure import (
    build_variables,
    constraint_variables,
    cost_variables,
    current_drive_variables,
    divertor_variables,
    fwbs_variables,
    global_variables,
    heat_transport_variables,
    impurity_radiation_module,
    neoclassics_variables,
    numerics,
    pfcoil_variables,
    physics_variables,
    rebco_variables,
    stellarator_configuration,
    stellarator_variables,
    structure_variables,
    superconducting_tf_coil_variables,
    tfcoil_variables,
    times_variables,
)
from process.exceptions import ProcessValueError
from process.physics import rether
from process.stellarator_config import load_stellarator_config

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
        physics,
        neoclassics,
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
            self.physics.outplas()
            self.stheat(True)
            self.stphys(True)
            self.stopt(True)

            # As stopt changes nd_plasma_electrons_vol_avg, te and b_plasma_toroidal_on_axis, stphys needs two calls
            # to correct for larger changes (it is only consistent after
            # two or three fix point iterations) call stphys here again, just to be sure.
            # This can be removed once the bad practice in stopt is removed!
            self.stphys(False)

            self.stdiv(True)
            self.stbild(True)
            self.stcoil(True)
            self.ststrc(True)
            self.stfwbs(True)

            self.power.tfpwr(output=True)
            self.buildings.run(output=True)
            self.vacuum.run(output=True)
            self.power.acpow(output=True)
            self.power.output_plant_electric_powers()

            return

        self.stnewconfig()
        self.stgeom()
        self.stphys(False)
        self.stopt(False)
        self.stcoil(False)
        self.stbild(False)
        self.ststrc(False)
        self.stfwbs(False)
        self.stdiv(False)

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
            ) = self.power_at_ignition_point(
                stellarator_variables.max_gyrotron_frequency,
                stellarator_variables.te0_ecrh_achievable,
            )

        stellarator_variables.first_call = False

    def stnewconfig(self):
        """author: J Lion, IPP Greifswald
        Routine to initialise the stellarator configuration

        Routine to initialise the stellarator configuration.
        This routine is called right before the calculation and could
        in principle overwrite variables from the input file.
        It overwrites rminor with rmajor and aspect ratio e.g.
        """

        load_stellarator_config(
            stellarator_variables.istell,
            Path(f"{global_variables.output_prefix}stella_conf.json"),
        )

        # If physics_variables.aspect ratio is not in numerics.ixc set it to default value
        # Or when you call it the first time
        if 1 not in numerics.ixc:
            physics_variables.aspect = (
                stellarator_configuration.stella_config_aspect_ref
            )

        # Set the physics_variables.rminor radius as result here.
        physics_variables.rminor = physics_variables.rmajor / physics_variables.aspect
        physics_variables.eps = 1.0e0 / physics_variables.aspect

        tfcoil_variables.n_tf_coils = (
            stellarator_configuration.stella_config_coilspermodule
            * stellarator_configuration.stella_config_symmetry
        )  # This overwrites tfcoil_variables.n_tf_coils in input file.

        #  Factors used to scale the reference point.
        stellarator_variables.f_r = (
            physics_variables.rmajor
            / stellarator_configuration.stella_config_rmajor_ref
        )  # Size scaling factor with respect to the reference calculation
        stellarator_variables.f_a = (
            physics_variables.rminor
            / stellarator_configuration.stella_config_rminor_ref
        )  # Size scaling factor with respect to the reference calculation

        stellarator_variables.f_aspect = (
            physics_variables.aspect
            / stellarator_configuration.stella_config_aspect_ref
        )
        stellarator_variables.f_n = tfcoil_variables.n_tf_coils / (
            stellarator_configuration.stella_config_coilspermodule
            * stellarator_configuration.stella_config_symmetry
        )  # Coil number factor
        stellarator_variables.f_b = (
            physics_variables.b_plasma_toroidal_on_axis
            / stellarator_configuration.stella_config_bt_ref
        )  # B-field scaling factor

    def stgeom(self):
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
            stellarator_variables.f_r
            * stellarator_variables.f_a**2
            * stellarator_configuration.stella_config_vol_plasma
        )

        # Plasma surface scaled from effective parameter:
        physics_variables.a_plasma_surface = (
            stellarator_variables.f_r
            * stellarator_variables.f_a
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

    def stopt(self, output: bool):
        """Routine to reiterate the physics loop
        author: J Lion, IPP Greifswald
        None
        This routine reiterates some physics modules.
        """

        physics_variables.nd_plasma_electrons_max = self.stdlim(
            physics_variables.b_plasma_toroidal_on_axis,
            physics_variables.p_plasma_loss_mw,
            physics_variables.rmajor,
            physics_variables.rminor,
        )

        # Calculates the ECRH parameters

        ne0_max_ECRH, bt_ecrh = self.stdlim_ecrh(
            stellarator_variables.max_gyrotron_frequency,
            physics_variables.b_plasma_toroidal_on_axis,
        )

        ne0_max_ECRH = min(physics_variables.nd_plasma_electron_on_axis, ne0_max_ECRH)
        bt_ecrh = min(physics_variables.b_plasma_toroidal_on_axis, bt_ecrh)

        if output:
            self.stopt_output(
                stellarator_variables.max_gyrotron_frequency,
                physics_variables.b_plasma_toroidal_on_axis,
                bt_ecrh,
                ne0_max_ECRH,
                stellarator_variables.te0_ecrh_achievable,
            )

    def stbild(self, output: bool):
        """
                Routine to determine the build of a stellarator machine
        author: P J Knight, CCFE, Culham Science Centre
        author: F Warmer, IPP Greifswald
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This routine determines the build of the stellarator machine.
        The values calculated are based on the mean minor radius, etc.,
        as the actual radial and vertical build thicknesses vary with
        toroidal angle.
        """
        if fwbs_variables.blktmodel > 0:
            build_variables.dr_blkt_inboard = (
                build_variables.blbuith
                + build_variables.blbmith
                + build_variables.blbpith
            )
            build_variables.dr_blkt_outboard = (
                build_variables.blbuoth
                + build_variables.blbmoth
                + build_variables.blbpoth
            )
            build_variables.dz_shld_upper = 0.5e0 * (
                build_variables.dr_shld_inboard + build_variables.dr_shld_outboard
            )

        #  Top/bottom blanket thickness

        build_variables.dz_blkt_upper = 0.5e0 * (
            build_variables.dr_blkt_inboard + build_variables.dr_blkt_outboard
        )

        # First Wall
        build_variables.dr_fw_inboard = (
            2.0e0 * fwbs_variables.radius_fw_channel + 2.0e0 * fwbs_variables.dr_fw_wall
        )
        build_variables.dr_fw_outboard = build_variables.dr_fw_inboard

        build_variables.dr_bore = physics_variables.rmajor - (
            build_variables.dr_cs
            + build_variables.dr_cs_tf_gap
            + build_variables.dr_tf_inboard
            + build_variables.dr_shld_vv_gap_inboard
            + build_variables.dr_vv_inboard
            + build_variables.dr_shld_inboard
            + build_variables.dr_blkt_inboard
            + build_variables.dr_fw_inboard
            + build_variables.dr_fw_plasma_gap_inboard
            + physics_variables.rminor
        )

        #  Radial build to centre of plasma (should be equal to physics_variables.rmajor)
        build_variables.rbld = (
            build_variables.dr_bore
            + build_variables.dr_cs
            + build_variables.dr_cs_tf_gap
            + build_variables.dr_tf_inboard
            + build_variables.dr_shld_vv_gap_inboard
            + build_variables.dr_vv_inboard
            + build_variables.dr_shld_inboard
            + build_variables.dr_blkt_inboard
            + build_variables.dr_fw_inboard
            + build_variables.dr_fw_plasma_gap_inboard
            + physics_variables.rminor
        )

        # Bc stellarators cannot scale physics_variables.rminor reasonably well an additional constraint equation is required,
        # that ensures that there is enough space between coils and plasma.
        build_variables.required_radial_space = (
            build_variables.dr_tf_inboard / 2.0e0
            + build_variables.dr_shld_vv_gap_inboard
            + build_variables.dr_vv_inboard
            + build_variables.dr_shld_inboard
            + build_variables.dr_blkt_inboard
            + build_variables.dr_fw_inboard
            + build_variables.dr_fw_plasma_gap_inboard
        )

        # derivative_min_LCFS_coils_dist  for how strong the stellarator shape changes wrt to aspect ratio
        build_variables.available_radial_space = stellarator_variables.f_r * (
            stellarator_configuration.stella_config_derivative_min_lcfs_coils_dist
            * stellarator_configuration.stella_config_rminor_ref
            * (1 / stellarator_variables.f_aspect - 1)
            + stellarator_configuration.stella_config_min_plasma_coil_distance
        )

        #  Radius to inner edge of inboard shield
        build_variables.rsldi = (
            physics_variables.rmajor
            - physics_variables.rminor
            - build_variables.dr_fw_plasma_gap_inboard
            - build_variables.dr_fw_inboard
            - build_variables.dr_blkt_inboard
            - build_variables.dr_shld_inboard
        )

        #  Radius to outer edge of outboard shield
        build_variables.rsldo = (
            physics_variables.rmajor
            + physics_variables.rminor
            + build_variables.dr_fw_plasma_gap_outboard
            + build_variables.dr_fw_outboard
            + build_variables.dr_blkt_outboard
            + build_variables.dr_shld_outboard
        )

        #  Thickness of outboard TF coil legs
        build_variables.dr_tf_outboard = build_variables.dr_tf_inboard

        #  Radius to centre of outboard TF coil legs

        build_variables.dr_shld_vv_gap_outboard = build_variables.gapomin
        build_variables.r_tf_outboard_mid = (
            build_variables.rsldo
            + build_variables.dr_vv_outboard
            + build_variables.dr_shld_vv_gap_outboard
            + 0.5e0 * build_variables.dr_tf_outboard
        )

        #  Height to inside edge of TF coil
        #  Roughly equal to average of (inboard build from TF coil to plasma
        #  centre) and (outboard build from plasma centre to TF coil)

        build_variables.z_tf_inside_half = 0.5e0 * (
            (
                build_variables.dr_shld_vv_gap_inboard
                + build_variables.dr_vv_inboard
                + build_variables.dr_shld_inboard
                + build_variables.dr_blkt_inboard
                + build_variables.dr_fw_inboard
                + build_variables.dr_fw_plasma_gap_inboard
                + physics_variables.rminor
            )
            + (
                physics_variables.rminor
                + build_variables.dr_fw_plasma_gap_outboard
                + build_variables.dr_fw_outboard
                + build_variables.dr_blkt_outboard
                + build_variables.dr_shld_outboard
                + build_variables.dr_vv_outboard
                + build_variables.dr_shld_vv_gap_outboard
            )
        )

        #  Outer divertor strike point radius, set equal to major radius

        build_variables.rspo = physics_variables.rmajor

        #  First wall area: scales with minor radius

        # Average minor radius of the first wall
        awall = physics_variables.rminor + 0.5e0 * (
            build_variables.dr_fw_plasma_gap_inboard
            + build_variables.dr_fw_plasma_gap_outboard
        )
        build_variables.a_fw_total = (
            physics_variables.a_plasma_surface * awall / physics_variables.rminor
        )

        if heat_transport_variables.ipowerflow == 0:
            build_variables.a_fw_total = (
                1.0e0 - fwbs_variables.fhole
            ) * build_variables.a_fw_total
        else:
            build_variables.a_fw_total = (
                1.0e0
                - fwbs_variables.fhole
                - fwbs_variables.f_ster_div_single
                - fwbs_variables.f_a_fw_outboard_hcd
            ) * build_variables.a_fw_total

        if output:
            #  Print out device build

            po.oheadr(self.outfile, "Radial Build")

            po.ovarre(
                self.outfile,
                "Avail. Space (m)",
                "(available_radial_space)",
                build_variables.available_radial_space,
            )
            po.ovarre(
                self.outfile,
                "Req. Space (m)",
                "(required_radial_space)",
                build_variables.required_radial_space,
            )
            po.ovarre(
                self.outfile, "f value: ", "(f_avspace)", build_variables.f_avspace
            )

            #     po.write(self.outfile,10)
            # 10  format(t43,'Thickness (m)',t60,'Radius (m)')

            radius = 0.0e0
            po.obuild(self.outfile, "Device centreline", 0.0e0, radius)

            drbild = (
                build_variables.dr_bore
                + build_variables.dr_cs
                + build_variables.dr_cs_tf_gap
            )
            radius = radius + drbild
            po.obuild(self.outfile, "Machine dr_bore", drbild, radius, "(dr_bore)")
            po.ovarre(
                self.outfile, "Machine build_variables.dr_bore (m)", "(dr_bore)", drbild
            )

            radius = radius + build_variables.dr_tf_inboard
            po.obuild(
                self.outfile,
                "Coil inboard leg",
                build_variables.dr_tf_inboard,
                radius,
                "(dr_tf_inboard)",
            )
            po.ovarre(
                self.outfile,
                "Coil inboard leg (m)",
                "(deltf)",
                build_variables.dr_tf_inboard,
            )

            radius = radius + build_variables.dr_shld_vv_gap_inboard
            po.obuild(
                self.outfile,
                "Gap",
                build_variables.dr_shld_vv_gap_inboard,
                radius,
                "(dr_shld_vv_gap_inboard)",
            )
            po.ovarre(
                self.outfile,
                "Gap (m)",
                "(dr_shld_vv_gap_inboard)",
                build_variables.dr_shld_vv_gap_inboard,
            )

            radius = radius + build_variables.dr_vv_inboard
            po.obuild(
                self.outfile,
                "Vacuum vessel",
                build_variables.dr_vv_inboard,
                radius,
                "(dr_vv_inboard)",
            )
            po.ovarre(
                self.outfile,
                "Vacuum vessel radial thickness (m)",
                "(dr_vv_inboard)",
                build_variables.dr_vv_inboard,
            )

            radius = radius + build_variables.dr_shld_inboard
            po.obuild(
                self.outfile,
                "Inboard shield",
                build_variables.dr_shld_inboard,
                radius,
                "(dr_shld_inboard)",
            )
            po.ovarre(
                self.outfile,
                "Inner radiation shield radial thickness (m)",
                "(dr_shld_inboard)",
                build_variables.dr_shld_inboard,
            )

            radius = radius + build_variables.dr_blkt_inboard
            po.obuild(
                self.outfile,
                "Inboard blanket",
                build_variables.dr_blkt_inboard,
                radius,
                "(dr_blkt_inboard)",
            )
            po.ovarre(
                self.outfile,
                "Inboard blanket radial thickness (m)",
                "(dr_blkt_inboard)",
                build_variables.dr_blkt_inboard,
            )

            radius = radius + build_variables.dr_fw_inboard
            po.obuild(
                self.outfile,
                "Inboard first wall",
                build_variables.dr_fw_inboard,
                radius,
                "(dr_fw_inboard)",
            )
            po.ovarre(
                self.outfile,
                "Inboard first wall radial thickness (m)",
                "(dr_fw_inboard)",
                build_variables.dr_fw_inboard,
            )

            radius = radius + build_variables.dr_fw_plasma_gap_inboard
            po.obuild(
                self.outfile,
                "Inboard scrape-off",
                build_variables.dr_fw_plasma_gap_inboard,
                radius,
                "(dr_fw_plasma_gap_inboard)",
            )
            po.ovarre(
                self.outfile,
                "Inboard scrape-off radial thickness (m)",
                "(dr_fw_plasma_gap_inboard)",
                build_variables.dr_fw_plasma_gap_inboard,
            )

            radius = radius + physics_variables.rminor
            po.obuild(
                self.outfile,
                "Plasma geometric centre",
                physics_variables.rminor,
                radius,
                "(rminor)",
            )

            radius = radius + physics_variables.rminor
            po.obuild(
                self.outfile,
                "Plasma outboard edge",
                physics_variables.rminor,
                radius,
                "(rminor)",
            )

            radius = radius + build_variables.dr_fw_plasma_gap_outboard
            po.obuild(
                self.outfile,
                "Outboard scrape-off",
                build_variables.dr_fw_plasma_gap_outboard,
                radius,
                "(dr_fw_plasma_gap_outboard)",
            )
            po.ovarre(
                self.outfile,
                "Outboard scrape-off radial thickness (m)",
                "(dr_fw_plasma_gap_outboard)",
                build_variables.dr_fw_plasma_gap_outboard,
            )

            radius = radius + build_variables.dr_fw_outboard
            po.obuild(
                self.outfile,
                "Outboard first wall",
                build_variables.dr_fw_outboard,
                radius,
                "(dr_fw_outboard)",
            )
            po.ovarre(
                self.outfile,
                "Outboard first wall radial thickness (m)",
                "(dr_fw_outboard)",
                build_variables.dr_fw_outboard,
            )

            radius = radius + build_variables.dr_blkt_outboard
            po.obuild(
                self.outfile,
                "Outboard blanket",
                build_variables.dr_blkt_outboard,
                radius,
                "(dr_blkt_outboard)",
            )
            po.ovarre(
                self.outfile,
                "Outboard blanket radial thickness (m)",
                "(dr_blkt_outboard)",
                build_variables.dr_blkt_outboard,
            )

            radius = radius + build_variables.dr_shld_outboard
            po.obuild(
                self.outfile,
                "Outboard shield",
                build_variables.dr_shld_outboard,
                radius,
                "(dr_shld_outboard)",
            )
            po.ovarre(
                self.outfile,
                "Outer radiation shield radial thickness (m)",
                "(dr_shld_outboard)",
                build_variables.dr_shld_outboard,
            )

            radius = radius + build_variables.dr_vv_outboard
            po.obuild(
                self.outfile,
                "Vacuum vessel",
                build_variables.dr_vv_outboard,
                radius,
                "(dr_vv_outboard)",
            )

            radius = radius + build_variables.dr_shld_vv_gap_outboard
            po.obuild(
                self.outfile,
                "Gap",
                build_variables.dr_shld_vv_gap_outboard,
                radius,
                "(dr_shld_vv_gap_outboard)",
            )
            po.ovarre(
                self.outfile,
                "Gap (m)",
                "(dr_shld_vv_gap_outboard)",
                build_variables.dr_shld_vv_gap_outboard,
            )

            radius = radius + build_variables.dr_tf_outboard
            po.obuild(
                self.outfile,
                "Coil outboard leg",
                build_variables.dr_tf_outboard,
                radius,
                "(dr_tf_outboard)",
            )
            po.ovarre(
                self.outfile,
                "Coil outboard leg radial thickness (m)",
                "(dr_tf_outboard)",
                build_variables.dr_tf_outboard,
            )

    def ststrc(self, output):
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
        # which needs to be precalculated (or calculated in PROCESS but this not done here)
        # The coil width is subtracted from that:
        # total_coil_width = b + 2* d_ic + 2* case_thickness_constant
        # total_coil_thickness = h + 2* d_ic + 2* case_thickness_constant

        # The following line is correct AS LONG AS we do not scale the coil sizes
        intercoil_surface = (
            stellarator_configuration.stella_config_coilsurface
            * stellarator_variables.f_r**2
            - tfcoil_variables.dx_tf_inboard_out_toroidal
            * stellarator_configuration.stella_config_coillength
            * stellarator_variables.f_r
            * stellarator_variables.f_n
        )

        # This 0.18 m is an effective thickness which is scaled with empirial 1.5 law. 5.6 T is reference point of Helias
        # The thickness 0.18m was obtained as a measured value from Schauer, F. and Bykov, V. design of Helias 5-B. (Nucl Fus. 2013)
        structure_variables.aintmass = (
            0.18e0
            * stellarator_variables.f_b**2
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

    def stdiv(self, output: bool):
        """Routine to call the stellarator divertor model
        author: P J Knight, CCFE, Culham Science Centre
        author: F Warmer, IPP Greifswald
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This routine calls the divertor model for a stellarator,
        developed by Felix Warmer.
        Stellarator Divertor Model for the Systems
        Code PROCESS, F. Warmer, 21/06/2013
        """
        Theta = stellarator_variables.flpitch  # ~bmn [rad] field line pitch
        r = physics_variables.rmajor
        p_div = physics_variables.p_plasma_separatrix_mw
        alpha = divertor_variables.anginc
        xi_p = divertor_variables.xpertin
        T_scrape = divertor_variables.tdiv

        #  Scrape-off temperature in Joules

        e = T_scrape * constants.ELECTRON_CHARGE

        #  Sound speed of particles (m/s)

        c_s = np.sqrt(e / (physics_variables.m_fuel_amu * constants.UMASS))

        #  Island size (m)

        w_r = 4.0e0 * np.sqrt(
            stellarator_variables.bmn
            * r
            / (stellarator_variables.shear * stellarator_variables.n_res)
        )

        #  Perpendicular (to plate) distance from X-point to divertor plate (m)

        Delta = stellarator_variables.f_w * w_r

        #  Length 'along' plasma (m)

        l_p = (
            2 * np.pi * r * (stellarator_variables.m_res) / stellarator_variables.n_res
        )

        #  Connection length from X-point to divertor plate (m)

        l_x_t = Delta / Theta

        #  Power decay length (m)

        l_q = np.sqrt(xi_p * (l_x_t / c_s))

        #  Channel broadening length (m)

        l_b = np.sqrt(xi_p * l_p / (c_s))

        #  Channel broadening factor

        f_x = 1.0e0 + (l_b / (l_p * Theta))

        #  Length of a single divertor plate (m)

        l_d = f_x * l_p * (Theta / alpha)

        #  Total length of divertor plates (m)

        l_t = 2.0e0 * stellarator_variables.n_res * l_d

        #  Wetted area (m2)

        a_eff = l_t * l_q

        #  Divertor plate width (m): assume total area is wetted area/stellarator_variables.fdivwet

        darea = a_eff / stellarator_variables.fdivwet
        l_w = darea / l_t

        #  Divertor heat load (MW/m2)

        q_div = stellarator_variables.f_asym * (p_div / a_eff)

        #  Transfer to global variables

        divertor_variables.pflux_div_heat_load_mw = q_div
        divertor_variables.a_div_surface_total = darea

        fwbs_variables.f_ster_div_single = darea / build_variables.a_fw_total

        if output:
            po.oheadr(self.outfile, "Divertor")

            po.ovarre(
                self.outfile,
                "Power to divertor (MW)",
                "(p_plasma_separatrix_mw.)",
                physics_variables.p_plasma_separatrix_mw,
            )
            po.ovarre(
                self.outfile,
                "Angle of incidence (deg)",
                "(anginc)",
                divertor_variables.anginc * 180.0e0 / np.pi,
            )
            po.ovarre(
                self.outfile,
                "Perp. heat transport coefficient (m2/s)",
                "(xpertin)",
                divertor_variables.xpertin,
            )
            po.ovarre(
                self.outfile,
                "Divertor plasma temperature (eV)",
                "(tdiv)",
                divertor_variables.tdiv,
            )
            po.ovarre(
                self.outfile,
                "Radiated power fraction in SOL",
                "(f_rad)",
                stellarator_variables.f_rad,
            )
            po.ovarre(
                self.outfile,
                "Heat load peaking factor",
                "(f_asym)",
                stellarator_variables.f_asym,
            )
            po.ovarin(
                self.outfile,
                "Poloidal resonance number",
                "(m_res)",
                stellarator_variables.m_res,
            )
            po.ovarin(
                self.outfile,
                "Toroidal resonance number",
                "(n_res)",
                stellarator_variables.n_res,
            )
            po.ovarre(
                self.outfile,
                "Relative radial field perturbation",
                "(bmn)",
                stellarator_variables.bmn,
            )
            po.ovarre(
                self.outfile,
                "Field line pitch (rad)",
                "(flpitch)",
                stellarator_variables.flpitch,
            )
            po.ovarre(
                self.outfile,
                "Island size fraction factor",
                "(f_w)",
                stellarator_variables.f_w,
            )
            po.ovarre(
                self.outfile,
                "Magnetic stellarator_variables.shear (/m)",
                "(shear)",
                stellarator_variables.shear,
            )
            po.ovarre(self.outfile, "Divertor wetted area (m2)", "(A_eff)", a_eff)
            po.ovarre(
                self.outfile,
                "Wetted area fraction of total plate area",
                "(fdivwet)",
                stellarator_variables.fdivwet,
            )
            po.ovarre(self.outfile, "Divertor plate length (m)", "(L_d)", l_d)
            po.ovarre(self.outfile, "Divertor plate width (m)", "(L_w)", l_w)
            po.ovarre(self.outfile, "Flux channel broadening factor", "(F_x)", f_x)
            po.ovarre(
                self.outfile, "Power decay width (cm)", "(100*l_q)", 100.0e0 * l_q
            )
            po.ovarre(self.outfile, "Island width (m)", "(w_r)", w_r)
            po.ovarre(
                self.outfile,
                "Perp. distance from X-point to plate (m)",
                "(Delta)",
                Delta,
            )
            po.ovarre(
                self.outfile,
                "Peak heat load (MW/m2)",
                "(pflux_div_heat_load_mw)",
                divertor_variables.pflux_div_heat_load_mw,
            )

    def blanket_neutronics(self):
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
        ) = self.sctfcoil_nuclear_heating_iter90()

        # blktlife calculation left entierly to availability
        # Cannot find calculation for vvhemax in CCFE blanket

    def stfwbs(self, output: bool):
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
            cost_variables.tlife,
        )

        #  First wall inboard, outboard areas (assume 50% of total each)
        build_variables.a_fw_inboard = 0.5e0 * build_variables.a_fw_total
        build_variables.a_fw_outboard = 0.5e0 * build_variables.a_fw_total

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
                ) = self.sctfcoil_nuclear_heating_iter90()

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
                    pnucfwbs * build_variables.a_fw_inboard / build_variables.a_fw_total
                )
                pnucfwbso = (
                    pnucfwbs
                    * build_variables.a_fw_outboard
                    / build_variables.a_fw_total
                )

                #  Radiation power incident on divertor (MW)

                fwbs_variables.p_div_rad_total_mw = (
                    physics_variables.p_plasma_rad_mw * fwbs_variables.f_ster_div_single
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
                    * build_variables.a_fw_inboard
                    / build_variables.a_fw_total
                )
                psurffwo = (
                    fwbs_variables.p_fw_rad_total_mw
                    * build_variables.a_fw_outboard
                    / build_variables.a_fw_total
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
                    * build_variables.a_fw_inboard
                    / build_variables.a_fw_total
                )
                pnucso = (
                    pnucbso
                    - pnucbzo
                    + (fwbs_variables.pnucloss + fwbs_variables.pradloss)
                    * build_variables.a_fw_outboard
                    / build_variables.a_fw_total
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

        #  heat_transport_variables.ipowerflow

        #  fwbs_variables.blktmodel

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
                build_variables.a_fw_total
                * (build_variables.dr_fw_inboard + build_variables.dr_fw_outboard)
                / 2.0e0
                * fwbs_variables.den_steel
                * (1.0e0 - fwbs_variables.fwclfr)
            )

            #  Surface areas adjacent to plasma

            coolvol = (
                coolvol
                + build_variables.a_fw_total
                * (build_variables.dr_fw_inboard + build_variables.dr_fw_outboard)
                / 2.0e0
                * fwbs_variables.fwclfr
            )

        else:
            fwbs_variables.m_fw_total = fwbs_variables.den_steel * (
                build_variables.a_fw_inboard
                * build_variables.dr_fw_inboard
                * (1.0e0 - f_a_fw_coolant_inboard)
                + build_variables.a_fw_outboard
                * build_variables.dr_fw_outboard
                * (1.0e0 - f_a_fw_coolant_outboard)
            )
            coolvol = (
                coolvol
                + build_variables.a_fw_inboard
                * build_variables.dr_fw_inboard
                * f_a_fw_coolant_inboard
                + build_variables.a_fw_outboard
                * build_variables.dr_fw_outboard
                * f_a_fw_coolant_outboard
            )

            #  Average first wall coolant fraction, only used by old routines
            #  in fispact.f90, safety.f90

            fwbs_variables.fwclfr = (
                build_variables.a_fw_inboard
                * build_variables.dr_fw_inboard
                * f_a_fw_coolant_inboard
                + build_variables.a_fw_outboard
                * build_variables.dr_fw_outboard
                * f_a_fw_coolant_outboard
            ) / (
                build_variables.a_fw_total
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
                    self.outfile, "Tritium breeding ratio", "(tbr)", fwbs_variables.tbr
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

            #     if (fwbs_variables.blktmodel == 0) :
            #         if ((fwbs_variables.blkttype == 1)or(fwbs_variables.blkttype == 2)) :
            #             po.write(self.outfile,601) vol_blkt_inboard, vol_blkt_outboard, vol_blkt_total,                m_blkt_total, f_a_blkt_cooling_channels, fbllipb, wtbllipb, fblli, m_blkt_lithium,                fblss, m_blkt_steel_total, fblvd, m_blkt_vanadium, vol_shld_inboard, vol_shld_outboard,                vol_shld_total, whtshld, vfshld, fwbs_variables.wpenshld
            #         else:  #  (also if ipowerflow=0)
            #             po.write(self.outfile,600) vol_blkt_inboard, vol_blkt_outboard, vol_blkt_total,                m_blkt_total, f_a_blkt_cooling_channels, fblbe, m_blkt_beryllium, fblli2o, m_blkt_li2o,                fblss, m_blkt_steel_total, fblvd, m_blkt_vanadium, vol_shld_inboard, vol_shld_outboard,                vol_shld_total, whtshld, vfshld, fwbs_variables.wpenshld

            #     else:
            #         po.write(self.outfile,602) vol_blkt_inboard, vol_blkt_outboard, vol_blkt_total, m_blkt_total, f_a_blkt_cooling_channels,             (fwbs_variables.vol_blkt_inboard/fwbs_variables.vol_blkt_total * build_variables.blbuith/build_variables.dr_blkt_inboard +             fwbs_variables.vol_blkt_outboard/fwbs_variables.vol_blkt_total * build_variables.blbuoth/build_variables.dr_blkt_outboard) * fblbe, m_blkt_beryllium,             (fwbs_variables.vol_blkt_inboard/fwbs_variables.vol_blkt_total * build_variables.blbuith/build_variables.dr_blkt_inboard +             fwbs_variables.vol_blkt_outboard/fwbs_variables.vol_blkt_total * build_variables.blbuoth/build_variables.dr_blkt_outboard) * fblbreed, whtblbreed,             fwbs_variables.vol_blkt_inboard/fwbs_variables.vol_blkt_total/build_variables.dr_blkt_inboard * (build_variables.blbuith * fwbs_variables.fblss             + build_variables.blbmith * (1.0e0-fwbs_variables.fblhebmi) + build_variables.blbpith * (1.0e0-fwbs_variables.fblhebpi)) +             fwbs_variables.vol_blkt_outboard/fwbs_variables.vol_blkt_total/build_variables.dr_blkt_outboard * (build_variables.blbuoth * fwbs_variables.fblss             + build_variables.blbmoth * (1.0e0-fwbs_variables.fblhebmo) + build_variables.blbpoth * (1.0e0-fwbs_variables.fblhebpo)),             m_blkt_steel_total,             vol_shld_inboard, vol_shld_outboard, vol_shld_total, whtshld, vfshld, fwbs_variables.wpenshld

            # 600 format(          t32,'volume (m3)',t45,'vol fraction',t62,'weight (kg)'/          t32,'-----------',t45,'------------',t62,'-----------'/          '    Inboard blanket' ,t32,1pe10.3,/          '    Outboard blanket' ,t32,1pe10.3,/          '    Total blanket' ,t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '       Blanket Be   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket Li2O ',t45,1pe10.3,t62,1pe10.3/          '       Blanket ss   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket Vd   ',t45,1pe10.3,t62,1pe10.3/          '    Inboard shield'  ,t32,1pe10.3,/          '    Outboard shield'  ,t32,1pe10.3,/          '    Primary shield',t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '    Penetration shield'        ,t62,1pe10.3)

            # 601 format(          t32,'volume (m3)',t45,'vol fraction',t62,'weight (kg)'/          t32,'-----------',t45,'------------',t62,'-----------'/          '    Inboard blanket' ,t32,1pe10.3,/          '    Outboard blanket' ,t32,1pe10.3,/          '    Total blanket' ,t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '       Blanket LiPb ',t45,1pe10.3,t62,1pe10.3/          '       Blanket Li   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket ss   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket Vd   ',t45,1pe10.3,t62,1pe10.3/          '    Inboard shield'  ,t32,1pe10.3,/          '    Outboard shield'  ,t32,1pe10.3,/          '    Primary shield',t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '    Penetration shield'        ,t62,1pe10.3)

            # 602 format(          t32,'volume (m3)',t45,'vol fraction',t62,'weight (kg)'/          t32,'-----------',t45,'------------',t62,'-----------'/          '    Inboard blanket' ,t32,1pe10.3,/          '    Outboard blanket' ,t32,1pe10.3,/          '    Total blanket' ,t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '       Blanket Be   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket breeder',t45,1pe10.3,t62,1pe10.3/          '       Blanket steel',t45,1pe10.3,t62,1pe10.3/          '    Inboard shield'  ,t32,1pe10.3,/          '    Outboard shield'  ,t32,1pe10.3,/          '    Primary shield',t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '    Penetration shield'        ,t62,1pe10.3)

            po.osubhd(self.outfile, "Other volumes, masses and areas :")
            po.ovarre(
                self.outfile,
                "First wall area (m2)",
                "(a_fw_total)",
                build_variables.a_fw_total,
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

    def sctfcoil_nuclear_heating_iter90(self):
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

            fpydt = cost_variables.cfactr * cost_variables.tlife
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

    def stcoil(self, output: bool):
        """Routine that performs the calculations for stellarator coils
        author: J Lion, IPP Greifswald
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This routine calculates the properties of the coils for
        a stellarator device.
        <P>Some precalculated effective parameters for a stellarator power
        plant design are used as the basis for the calculations. The coils
        are assumed to be a fixed shape, but are scaled in size
        appropriately for the machine being modelled.
        """
        r_coil_major = (
            stellarator_configuration.stella_config_coil_rmajor
            * stellarator_variables.f_r
        )
        r_coil_minor = (
            stellarator_configuration.stella_config_coil_rminor
            * stellarator_variables.f_r
        )

        ########################################################################################
        # Winding Pack Geometry: for one conductor
        #
        # This one conductor will just be multiplied later to fit the winding pack size.
        #
        # [m] Dimension of square cable space inside insulation
        #     and case of the conduit of each turn
        dx_tf_turn_cable_space_average = tfcoil_variables.dx_tf_turn_general - 2.0e0 * (
            tfcoil_variables.dx_tf_turn_steel + tfcoil_variables.dx_tf_turn_insulation
        )  # dx_tf_turn_cable_space_average = t_w
        if dx_tf_turn_cable_space_average < 0:
            print(
                "dx_tf_turn_cable_space_average is negative. Check t_turn, tfcoil_variables.dx_tf_turn_steel and dx_tf_turn_insulation."
            )
        # [m^2] Cross-sectional area of cable space per turn
        tfcoil_variables.a_tf_turn_cable_space_no_void = (
            0.9e0 * dx_tf_turn_cable_space_average**2
        )  # 0.9 to include some rounded corners. (tfcoil_variables.a_tf_turn_cable_space_no_void = pi (dx_tf_turn_cable_space_average/2)**2 = pi/4 *dx_tf_turn_cable_space_average**2 for perfect round conductor). This factor depends on how round the corners are.
        # [m^2] Cross-sectional area of conduit case per turn
        tfcoil_variables.a_tf_turn_steel = (
            dx_tf_turn_cable_space_average + 2.0e0 * tfcoil_variables.dx_tf_turn_steel
        ) ** 2 - tfcoil_variables.a_tf_turn_cable_space_no_void
        #######################################################################################

        #######################################################################################
        # Winding Pack total size:
        #
        # Total coil current (MA)
        coilcurrent = (
            stellarator_variables.f_b
            * stellarator_configuration.stella_config_i0
            * stellarator_variables.f_r
            / stellarator_variables.f_n
        )
        stellarator_variables.f_i = (
            coilcurrent / stellarator_configuration.stella_config_i0
        )

        n_it = 200  # number of iterations

        rhs = np.zeros((n_it,))
        lhs = np.zeros((n_it,))
        jcrit_vector = np.zeros((n_it,))
        wp_width_r = np.zeros((n_it,))
        b_max_k = np.zeros((n_it,))

        for k in range(n_it):
            # Sample coil winding pack
            wp_width_r[k] = (r_coil_minor / 40.0e0) + (k / (n_it - 1e0)) * (
                r_coil_minor / 1.0e0 - r_coil_minor / 40.0e0
            )
            if tfcoil_variables.i_tf_sc_mat == 6:
                wp_width_r[k] = (r_coil_minor / 150.0e0) + (k / (n_it - 1e0)) * (
                    r_coil_minor / 1.0e0 - r_coil_minor / 150.0e0
                )

            #  B-field calculation
            b_max_k[k] = self.bmax_from_awp(
                wp_width_r[k],
                coilcurrent,
                tfcoil_variables.n_tf_coils,
                r_coil_major,
                r_coil_minor,
            )

            # jcrit for this bmax:
            jcrit_vector[k] = self.jcrit_frommaterial(
                b_max_k[k],
                tfcoil_variables.tftmp + 1.5,
                tfcoil_variables.i_tf_sc_mat,
                tfcoil_variables.b_crit_upper_nbti,
                tfcoil_variables.bcritsc,
                tfcoil_variables.f_a_tf_turn_cable_copper,
                tfcoil_variables.fhts,
                tfcoil_variables.t_crit_nbti,
                tfcoil_variables.tcritsc,
                tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
                tfcoil_variables.j_tf_wp,
            )  # Get here a temperature margin of 1.5K.

        # The operation current density weighted with the global iop/icrit fraction
        lhs[:] = constraint_variables.fiooic * jcrit_vector

        # Conduct fraction of conduit * Superconductor fraction in conductor
        f_scu = (
            (
                tfcoil_variables.a_tf_turn_cable_space_no_void
                * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
            )
            / (tfcoil_variables.dx_tf_turn_general**2)
            * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_copper)
        )  # fraction that is SC of wp.
        # print *, "f_scu. ",f_scu,"Awp min: ",Awp(1)

        rhs[:] = coilcurrent / (
            wp_width_r**2 / stellarator_configuration.stella_config_wp_ratio * f_scu
        )  # f_scu should be the fraction of the sc that is in the winding pack.

        wp_width_r_min = (
            r_coil_minor / 10.0e0
        ) ** 2  # Initial guess for intersection routine
        if tfcoil_variables.i_tf_sc_mat == 6:
            wp_width_r_min = (
                r_coil_minor / 20.0e0
            ) ** 2  # If REBCO, : start at smaller winding pack ratios

        # Find the intersection between LHS and RHS (or: how much awp do I need to get to the desired coil current)
        wp_width_r_min = self.intersect(
            wp_width_r, lhs, wp_width_r, rhs, wp_width_r_min
        )

        # Maximum field at superconductor surface (T)
        wp_width_r_min = max(tfcoil_variables.dx_tf_turn_general**2, wp_width_r_min)

        # Recalculate tfcoil_variables.b_tf_inboard_peak_symmetric at the found awp_min:
        tfcoil_variables.b_tf_inboard_peak_symmetric = self.bmax_from_awp(
            wp_width_r_min,
            coilcurrent,
            tfcoil_variables.n_tf_coils,
            r_coil_major,
            r_coil_minor,
        )

        # Winding pack toroidal, radial cross-sections (m)
        awp_tor = (
            wp_width_r_min / stellarator_configuration.stella_config_wp_ratio
        )  # Toroidal dimension
        awp_rad = wp_width_r_min  # Radial dimension

        tfcoil_variables.dx_tf_wp_primary_toroidal = (
            awp_tor  # [m] toroidal thickness of winding pack
        )
        tfcoil_variables.dx_tf_wp_secondary_toroidal = (
            awp_tor  # [m] toroidal thickness of winding pack (region in front)
        )
        tfcoil_variables.dr_tf_wp_with_insulation = (
            awp_rad  # [m] radial thickness of winding pack
        )

        #  [m^2] winding-pack cross sectional area including insulation (not global)
        a_tf_wp_with_insulation = (
            tfcoil_variables.dr_tf_wp_with_insulation
            + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
        ) * (
            tfcoil_variables.dx_tf_wp_primary_toroidal
            + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
        )

        a_tf_wp_no_insulation = (
            awp_tor * awp_rad
        )  # [m^2] winding-pack cross sectional area
        tfcoil_variables.j_tf_wp = (
            coilcurrent * 1.0e6 / a_tf_wp_no_insulation
        )  # [A/m^2] winding pack current density
        tfcoil_variables.n_tf_coil_turns = (
            a_tf_wp_no_insulation / (tfcoil_variables.dx_tf_turn_general**2)
        )  # estimated number of turns for a given turn size (not global). Take at least 1.
        tfcoil_variables.c_tf_turn = (
            coilcurrent * 1.0e6 / tfcoil_variables.n_tf_coil_turns
        )  # [A] current per turn - estimation
        # [m^2] Total conductor cross-sectional area, taking account of void area
        tfcoil_variables.a_tf_wp_conductor = (
            tfcoil_variables.a_tf_turn_cable_space_no_void
            * tfcoil_variables.n_tf_coil_turns
            * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
        )
        # [m^2] Void area in cable, for He
        tfcoil_variables.a_tf_wp_extra_void = (
            tfcoil_variables.a_tf_turn_cable_space_no_void
            * tfcoil_variables.n_tf_coil_turns
            * tfcoil_variables.f_a_tf_turn_cable_space_extra_void
        )
        # [m^2] Insulation area (not including ground-wall)
        tfcoil_variables.a_tf_coil_wp_turn_insulation = (
            tfcoil_variables.n_tf_coil_turns
            * (
                tfcoil_variables.dx_tf_turn_general**2
                - tfcoil_variables.a_tf_turn_steel
                - tfcoil_variables.a_tf_turn_cable_space_no_void
            )
        )
        # [m^2] Structure area for cable
        tfcoil_variables.a_tf_wp_steel = (
            tfcoil_variables.n_tf_coil_turns * tfcoil_variables.a_tf_turn_steel
        )
        # End of winding pack calculations
        #######################################################################################

        #######################################################################################
        #  Casing calculations
        #
        # Coil case thickness (m). Here assumed to be constant
        # until something better comes up.
        # case_thickness_constant = tfcoil_variables.dr_tf_nose_case #0.2e0 # #? Leave this constant for now... Check this## Should be scaled with forces I think.
        #  For now assumed to be constant in a bolted plate model.
        #
        tfcoil_variables.dr_tf_plasma_case = (
            tfcoil_variables.dr_tf_nose_case
        )  # [m] coil case thickness outboard distance (radial)
        # dr_tf_nose_case = case_thickness_constant/2.0e0 # [m] coil case thickness inboard distance  (radial).
        tfcoil_variables.dx_tf_side_case_min = (
            tfcoil_variables.dr_tf_nose_case
        )  # [m] coil case thickness toroidal distance (toroidal)

        # End of casing calculations
        #######################################################################################

        #######################################################################################
        #  Port calculations
        #
        #  Maximal toroidal port size (vertical ports) (m)
        #  The maximal distance is correct but the vertical extension of this port is not clear#
        #  This is simplified for now and can be made more accurate in the future#
        stellarator_variables.vporttmax = (
            0.4e0
            * stellarator_configuration.stella_config_max_portsize_width
            * stellarator_variables.f_r
            / stellarator_variables.f_n
        )  # This is not accurate yet. Needs more insight#

        #  Maximal poloidal port size (vertical ports) (m)
        stellarator_variables.vportpmax = (
            2.0 * stellarator_variables.vporttmax
        )  # Simple approximation

        #  Maximal vertical port clearance area (m2)
        stellarator_variables.vportamax = (
            stellarator_variables.vporttmax * stellarator_variables.vportpmax
        )

        #  Horizontal ports
        #  Maximal toroidal port size (horizontal ports) (m)
        stellarator_variables.hporttmax = (
            0.8e0
            * stellarator_configuration.stella_config_max_portsize_width
            * stellarator_variables.f_r
            / stellarator_variables.f_n
        )  # Factor 0.8 to take the variation with height into account

        #  Maximal poloidal port size (horizontal ports) (m)
        stellarator_variables.hportpmax = (
            2.0e0 * stellarator_variables.hporttmax
        )  # Simple approximation

        #  Maximal horizontal port clearance area (m2)
        stellarator_variables.hportamax = (
            stellarator_variables.hporttmax * stellarator_variables.hportpmax
        )
        # End of port calculations
        #######################################################################################

        #######################################################################################
        #  General Coil Geometry values
        #
        tfcoil_variables.dx_tf_inboard_out_toroidal = (
            tfcoil_variables.dx_tf_wp_primary_toroidal
            + 2.0e0 * tfcoil_variables.dx_tf_side_case_min
            + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
        )  # [m] Thickness of inboard leg in toroidal direction

        build_variables.dr_tf_inboard = (
            tfcoil_variables.dr_tf_nose_case
            + tfcoil_variables.dr_tf_wp_with_insulation
            + tfcoil_variables.dr_tf_plasma_case
            + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
        )  # [m] Thickness of inboard leg in radial direction
        build_variables.dr_tf_outboard = (
            tfcoil_variables.dr_tf_nose_case
            + tfcoil_variables.dr_tf_wp_with_insulation
            + tfcoil_variables.dr_tf_plasma_case
            + 2.0e0 * tfcoil_variables.dx_tf_wp_insulation
        )  # [m] Thickness of outboard leg in radial direction (same as inboard)
        tfcoil_variables.a_tf_leg_outboard = (
            build_variables.dr_tf_inboard * tfcoil_variables.dx_tf_inboard_out_toroidal
        )  # [m^2] overall coil cross-sectional area (assuming inboard and
        #       outboard leg are the same)
        tfcoil_variables.a_tf_coil_inboard_case = (
            build_variables.dr_tf_inboard * tfcoil_variables.dx_tf_inboard_out_toroidal
        ) - a_tf_wp_with_insulation  # [m^2] Cross-sectional area of surrounding case

        tfcoil_variables.tfocrn = (
            0.5e0 * tfcoil_variables.dx_tf_inboard_out_toroidal
        )  # [m] Half-width of side of coil nearest torus centreline
        tfcoil_variables.tficrn = (
            0.5e0 * tfcoil_variables.dx_tf_inboard_out_toroidal
        )  # [m] Half-width of side of coil nearest plasma

        # [m^2] Total surface area of coil side facing plasma: inboard region
        tfcoil_variables.tfsai = (
            tfcoil_variables.n_tf_coils
            * tfcoil_variables.dx_tf_inboard_out_toroidal
            * 0.5e0
            * tfcoil_variables.len_tf_coil
        )
        # [m^2] Total surface area of coil side facing plasma: outboard region
        tfcoil_variables.tfsao = (
            tfcoil_variables.tfsai
        )  # depends, how 'inboard' and 'outboard' are defined

        # [m] Minimal distance in toroidal direction between two stellarator coils (from mid to mid)
        # Consistency with coil width is checked in constraint equation 82
        tfcoil_variables.toroidalgap = (
            stellarator_configuration.stella_config_dmin
            * (r_coil_major - r_coil_minor)
            / (
                stellarator_configuration.stella_config_coil_rmajor
                - stellarator_configuration.stella_config_coil_rminor
            )
        )
        # Left-Over coil gap between two coils (m)
        coilcoilgap = (
            tfcoil_variables.toroidalgap - tfcoil_variables.dx_tf_inboard_out_toroidal
        )

        #  Variables for ALL coils.
        tfcoil_variables.a_tf_inboard_total = (
            tfcoil_variables.n_tf_coils * tfcoil_variables.a_tf_leg_outboard
        )  # [m^2] Total area of all coil legs (midplane)
        tfcoil_variables.c_tf_total = (
            tfcoil_variables.n_tf_coils * coilcurrent * 1.0e6
        )  # [A] Total current in ALL coils
        tfcoil_variables.oacdcp = (
            tfcoil_variables.c_tf_total / tfcoil_variables.a_tf_inboard_total
        )  # [A / m^2] overall current density
        tfcoil_variables.r_b_tf_inboard_peak = (
            r_coil_major - r_coil_minor + awp_rad
        )  # [m] radius of peak field occurrence, average
        # jlion: not sure what this will be used for. Not very
        # useful for stellarators

        # This uses the reference value for the inductance and scales it with a^2/R (toroid inductance scaling)
        inductance = (
            stellarator_configuration.stella_config_inductance
            / stellarator_variables.f_r
            * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor) ** 2
            * stellarator_variables.f_n**2
        )
        tfcoil_variables.e_tf_magnetic_stored_total_gj = (
            0.5e0
            * (
                stellarator_configuration.stella_config_inductance
                / stellarator_variables.f_r
                * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
                ** 2
                * stellarator_variables.f_n**2
            )
            * (tfcoil_variables.c_tf_total / tfcoil_variables.n_tf_coils) ** 2
            * 1.0e-9
        )  # [GJ] Total magnetic energy

        #  Coil dimensions
        build_variables.z_tf_inside_half = (
            0.5e0
            * stellarator_configuration.stella_config_maximal_coil_height
            * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
        )  # [m] maximum half-height of coil
        r_tf_inleg_mid = (
            r_coil_major - r_coil_minor
        )  # This is not very well defined for a stellarator.
        # Though, this is taken as an average value.
        tf_total_h_width = (
            r_coil_minor  # ? not really sure what this is supposed to be. Estimated as
        )
        # the average minor coil radius

        tfborev = (
            2.0e0 * build_variables.z_tf_inside_half
        )  # [m] estimated vertical coil dr_bore

        tfcoil_variables.len_tf_coil = (
            stellarator_configuration.stella_config_coillength
            * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
            / tfcoil_variables.n_tf_coils
        )  # [m] estimated average length of a coil

        # [m^2] Total surface area of toroidal shells covering coils
        tfcoil_variables.tfcryoarea = (
            stellarator_configuration.stella_config_coilsurface
            * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor) ** 2
            * 1.1e0
        )  # 1.1 to scale it out a bit. Should be coupled to winding pack maybe.

        # Minimal bending radius:
        min_bending_radius = (
            stellarator_configuration.stella_config_min_bend_radius
            * stellarator_variables.f_r
            * 1.0
            / (1.0 - tfcoil_variables.dr_tf_wp_with_insulation / (2.0 * r_coil_minor))
        )

        # End of general coil geometry values
        #######################################################################################

        #######################################################################################
        #  Masses of conductor constituents
        #
        # [kg] Mass of case
        #  (no need for correction factors as is the case for tokamaks)
        # This is only correct if the winding pack is 'thin' (len_tf_coil>>sqrt(tfcoil_variables.a_tf_coil_inboard_case)).
        tfcoil_variables.m_tf_coil_case = (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.a_tf_coil_inboard_case
            * tfcoil_variables.den_tf_coil_case
        )
        # Mass of ground-wall insulation [kg]
        # (assumed to be same density/material as conduit insulation)
        tfcoil_variables.m_tf_coil_wp_insulation = (
            tfcoil_variables.len_tf_coil
            * (a_tf_wp_with_insulation - a_tf_wp_no_insulation)
            * tfcoil_variables.den_tf_wp_turn_insulation
        )
        # [kg] mass of Superconductor
        tfcoil_variables.m_tf_coil_superconductor = (
            (
                tfcoil_variables.len_tf_coil
                * tfcoil_variables.n_tf_coil_turns
                * tfcoil_variables.a_tf_turn_cable_space_no_void
                * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
                * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_copper)
                - tfcoil_variables.len_tf_coil
                * tfcoil_variables.a_tf_wp_coolant_channels
            )
            * tfcoil_variables.dcond[tfcoil_variables.i_tf_sc_mat - 1]
        )  # a_tf_wp_coolant_channels is 0 for a stellarator. but keep this term for now.
        # [kg] mass of Copper in conductor
        tfcoil_variables.m_tf_coil_copper = (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.n_tf_coil_turns
            * tfcoil_variables.a_tf_turn_cable_space_no_void
            * (1.0e0 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void)
            * tfcoil_variables.f_a_tf_turn_cable_copper
            - tfcoil_variables.len_tf_coil * tfcoil_variables.a_tf_wp_coolant_channels
        ) * constants.den_copper
        # [kg] mass of Steel conduit (sheath)
        tfcoil_variables.m_tf_wp_steel_conduit = (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.n_tf_coil_turns
            * tfcoil_variables.a_tf_turn_steel
            * fwbs_variables.den_steel
        )
        # if (i_tf_sc_mat==6)   tfcoil_variables.m_tf_wp_steel_conduit = fcondsteel * a_tf_wp_no_insulation *tfcoil_variables.len_tf_coil* fwbs_variables.den_steel
        # Conduit insulation mass [kg]
        # (tfcoil_variables.a_tf_coil_wp_turn_insulation already contains tfcoil_variables.n_tf_coil_turns)
        tfcoil_variables.m_tf_coil_wp_turn_insulation = (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.a_tf_coil_wp_turn_insulation
            * tfcoil_variables.den_tf_wp_turn_insulation
        )
        # [kg] Total conductor mass
        tfcoil_variables.m_tf_coil_conductor = (
            tfcoil_variables.m_tf_coil_superconductor
            + tfcoil_variables.m_tf_coil_copper
            + tfcoil_variables.m_tf_wp_steel_conduit
            + tfcoil_variables.m_tf_coil_wp_turn_insulation
        )
        # [kg] Total coil mass
        tfcoil_variables.m_tf_coils_total = (
            tfcoil_variables.m_tf_coil_case
            + tfcoil_variables.m_tf_coil_conductor
            + tfcoil_variables.m_tf_coil_wp_insulation
        ) * tfcoil_variables.n_tf_coils
        # End of general coil geometry values
        #######################################################################################

        #######################################################################################
        # Quench protection:
        #
        # This copied from the tokamak module:
        # Radial position of vacuum vessel [m]
        radvv = (
            physics_variables.rmajor
            - physics_variables.rminor
            - build_variables.dr_fw_plasma_gap_inboard
            - build_variables.dr_fw_inboard
            - build_variables.dr_blkt_inboard
            - build_variables.dr_shld_blkt_gap
            - build_variables.dr_shld_inboard
        )

        # Actual VV force density
        # Based on reference values from W-7X:
        # Bref = 3;
        # Iref = 1.3*50;
        # aref = 0.92;
        # \[Tau]ref = 1.;
        # Rref = 5.2;
        # dref = 14*10^-3;

        # NOTE: original implementation used taucq which used a EUROfusion
        # constant in the calculation. This was the minimum allowed quench time.
        # Replacing with the actual quench time.
        # MN/m^3
        f_vv_actual = (
            2.54e6
            * (3e0 * 1.3e0 * 50e0 * 0.92e0**2e0)
            / (1e0 * 5.2e0 * 0.014e0)
            * (
                physics_variables.b_plasma_toroidal_on_axis
                * tfcoil_variables.c_tf_total
                * physics_variables.rminor**2
                / (
                    (build_variables.dr_vv_inboard + build_variables.dr_vv_outboard)
                    / 2
                    * tfcoil_variables.t_tf_superconductor_quench
                    * radvv
                )
            )
            ** (-1)
        )

        # N/m^2
        # is the vv width the correct length to multiply by to turn the
        # force density into a stress?
        superconducting_tf_coil_variables.vv_stress_quench = (
            f_vv_actual
            * 1e6
            * ((build_variables.dr_vv_inboard + build_variables.dr_vv_outboard) / 2)
        )

        # the conductor fraction is meant of the cable space#

        vd = self.u_max_protect_v(
            tfcoil_variables.e_tf_magnetic_stored_total_gj
            / tfcoil_variables.n_tf_coils
            * 1.0e9,
            tfcoil_variables.t_tf_superconductor_quench,
            tfcoil_variables.c_tf_turn,
        )

        # comparison
        # the new quench protection routine, see #1047
        tfcoil_variables.j_tf_wp_quench_heat_max = self.j_max_protect_am2(
            tfcoil_variables.t_tf_superconductor_quench,
            0.0e0,
            tfcoil_variables.f_a_tf_turn_cable_copper,
            1 - tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
            tfcoil_variables.tftmp,
            tfcoil_variables.a_tf_turn_cable_space_no_void,
            tfcoil_variables.dx_tf_turn_general**2,
        )

        # print *, "Jmax, comparison: ", j_tf_wp_quench_heat_max, "  ", jwdgpro2,"  ",j_tf_wp/j_tf_wp_quench_heat_max, "   , tfcoil_variables.t_tf_superconductor_quench: ",t_tf_superconductor_quench, " tfcoil_variables.f_a_tf_turn_cable_copper: ",f_a_tf_turn_cable_copper
        # print *, "a_tf_turn_cable_space_no_void: ", tfcoil_variables.a_tf_turn_cable_space_no_void
        # Also give the copper area for REBCO quench calculations:
        rebco_variables.coppera_m2 = (
            coilcurrent
            * 1.0e6
            / (
                tfcoil_variables.a_tf_wp_conductor
                * tfcoil_variables.f_a_tf_turn_cable_copper
            )
        )
        tfcoil_variables.v_tf_coil_dump_quench_kv = vd / 1.0e3  # Dump voltage
        #
        #######################################################################################

        # Forces scaling #
        tfcoil_variables.max_force_density = (
            stellarator_configuration.stella_config_max_force_density
            * stellarator_variables.f_i
            / stellarator_variables.f_n
            * tfcoil_variables.b_tf_inboard_peak_symmetric
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_wp_area
            / a_tf_wp_no_insulation
        )

        # Approximate, very simple maxiumum stress: (needed for limitation of icc 32)
        tfcoil_variables.sig_tf_wp = (
            tfcoil_variables.max_force_density
            * tfcoil_variables.dr_tf_wp_with_insulation
            * 1.0e6
        )  # in Pa

        # Units: MN/m
        max_force_density_mnm = (
            stellarator_configuration.stella_config_max_force_density_mnm
            * stellarator_variables.f_i
            / stellarator_variables.f_n
            * tfcoil_variables.b_tf_inboard_peak_symmetric
            / stellarator_configuration.stella_config_wp_bmax
        )
        #
        max_lateral_force_density = (
            stellarator_configuration.stella_config_max_lateral_force_density
            * stellarator_variables.f_i
            / stellarator_variables.f_n
            * tfcoil_variables.b_tf_inboard_peak_symmetric
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_wp_area
            / a_tf_wp_no_insulation
        )
        max_radial_force_density = (
            stellarator_configuration.stella_config_max_radial_force_density
            * stellarator_variables.f_i
            / stellarator_variables.f_n
            * tfcoil_variables.b_tf_inboard_peak_symmetric
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_wp_area
            / a_tf_wp_no_insulation
        )
        #
        # F = f*V = B*j*V \propto B/B0 * I/I0 * A0/A * A/A0 * len/len0
        centering_force_max_mn = (
            stellarator_configuration.stella_config_centering_force_max_mn
            * stellarator_variables.f_i
            / stellarator_variables.f_n
            * tfcoil_variables.b_tf_inboard_peak_symmetric
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_coillength
            / tfcoil_variables.n_tf_coils
            / tfcoil_variables.len_tf_coil
        )
        centering_force_min_mn = (
            stellarator_configuration.stella_config_centering_force_min_mn
            * stellarator_variables.f_i
            / stellarator_variables.f_n
            * tfcoil_variables.b_tf_inboard_peak_symmetric
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_coillength
            / tfcoil_variables.n_tf_coils
            / tfcoil_variables.len_tf_coil
        )
        centering_force_avg_mn = (
            stellarator_configuration.stella_config_centering_force_avg_mn
            * stellarator_variables.f_i
            / stellarator_variables.f_n
            * tfcoil_variables.b_tf_inboard_peak_symmetric
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_coillength
            / tfcoil_variables.n_tf_coils
            / tfcoil_variables.len_tf_coil
        )
        #
        ####################################

        if output:
            self.stcoil_output(
                a_tf_wp_no_insulation,
                centering_force_avg_mn,
                centering_force_max_mn,
                centering_force_min_mn,
                coilcoilgap,
                rebco_variables.coppera_m2,
                rebco_variables.coppera_m2_max,
                f_scu,
                f_vv_actual,
                constraint_variables.fiooic,
                inductance,
                tfcoil_variables.max_force_density,
                max_force_density_mnm,
                max_lateral_force_density,
                max_radial_force_density,
                min_bending_radius,
                r_coil_major,
                r_coil_minor,
                r_tf_inleg_mid,
                tfcoil_variables.sig_tf_wp,
                tfcoil_variables.dx_tf_turn_general,
                tfcoil_variables.t_tf_superconductor_quench,
                tf_total_h_width,
                tfborev,
                tfcoil_variables.toroidalgap,
                tfcoil_variables.v_tf_coil_dump_quench_max_kv,
                tfcoil_variables.v_tf_coil_dump_quench_kv,
            )

    def u_max_protect_v(self, tfes, tdump, aio):
        """tfes : input real : Energy stored in one TF coil (J)
        tdump : input real : Dump time (sec)
        aio : input real : Operating current (A)
        """
        return 2 * tfes / (tdump * aio)

    def j_max_protect_am2(self, tau_quench, t_detect, fcu, fcond, temp, acs, aturn):
        temp_k = [4, 14, 24, 34, 44, 54, 64, 74, 84, 94, 104, 114, 124]
        q_cu_array_sa2m4 = [
            1.08514e17,
            1.12043e17,
            1.12406e17,
            1.05940e17,
            9.49741e16,
            8.43757e16,
            7.56346e16,
            6.85924e16,
            6.28575e16,
            5.81004e16,
            5.40838e16,
            5.06414e16,
            4.76531e16,
        ]
        q_he_array_sa2m4 = [
            3.44562e16,
            9.92398e15,
            4.90462e15,
            2.41524e15,
            1.26368e15,
            7.51617e14,
            5.01632e14,
            3.63641e14,
            2.79164e14,
            2.23193e14,
            1.83832e14,
            1.54863e14,
            1.32773e14,
        ]

        q_he = np.interp(temp, temp_k, q_he_array_sa2m4)
        q_cu = np.interp(temp, temp_k, q_cu_array_sa2m4)

        # This leaves out the contribution from the superconductor fraction for now
        return (acs / aturn) * np.sqrt(
            1
            / (0.5 * tau_quench + t_detect)
            * (fcu**2 * fcond**2 * q_cu + fcu * fcond * (1 - fcond) * q_he)
        )

    def jcrit_frommaterial(
        self,
        bmax,
        thelium,
        i_tf_sc_mat,
        b_crit_upper_nbti,
        bcritsc,
        f_a_tf_turn_cable_copper,
        fhts,
        t_crit_nbti,
        tcritsc,
        f_a_tf_turn_cable_space_extra_void,
        jwp,
    ):
        strain = -0.005  # for now a small value
        fhe = f_a_tf_turn_cable_space_extra_void  # this is helium fraction in the superconductor (set it to the fixed global variable here)

        fcu = f_a_tf_turn_cable_copper  # f_a_tf_turn_cable_copper is a global variable. Is the copper fraction
        # of a cable conductor.

        if i_tf_sc_mat == 1:  # ITER Nb3Sn critical surface parameterization
            bc20m = 32.97  # these are values taken from sctfcoil.f90
            tc0m = 16.06

            #  j_crit_sc returned by itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper
            if bmax > bc20m:
                j_crit_sc = 1.0e-9  # Set to a small nonzero value
            else:
                (
                    j_crit_sc,
                    bcrit,
                    tcrit,
                ) = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)

            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1.0 - fcu) * (1.0e0 - fhe)

            # This is needed right now. Can we change it later?
            j_crit_sc = max(1.0e-9, j_crit_sc)
            j_crit_cable = max(1.0e-9, j_crit_cable)

        elif i_tf_sc_mat == 2:
            # Bi-2212 high temperature superconductor parameterization
            #  Current density in a strand of Bi-2212 conductor
            #  N.B. jcrit returned by bi2212 is the critical current density
            #  in the strand, not just the superconducting portion.
            #  The parameterization for j_crit_cable assumes a particular strand
            #  composition that does not require a user-defined copper fraction,
            #  so this is irrelevant in this model

            jstrand = jwp / (1 - fhe)
            #  jstrand = 0  # as far as I can tell this will always be 0
            #  because jwp was never set in fortran (so 0)

            j_crit_cable, tmarg = superconductors.bi2212(
                bmax, jstrand, thelium, fhts
            )  # bi2212 outputs j_crit_cable
            j_crit_sc = j_crit_cable / (1 - fcu)
            tcrit = thelium + tmarg
        elif i_tf_sc_mat == 3:  # NbTi data
            bc20m = 15.0
            tc0m = 9.3
            c0 = 1.0

            if bmax > bc20m:
                j_crit_sc = 1.0e-9  # Set to a small nonzero value
            else:
                j_crit_sc, tcrit = superconductors.jcrit_nbti(
                    thelium,
                    bmax,
                    c0,
                    bc20m,
                    tc0m,
                )
                # I dont need tcrit here so dont use it.

            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1 - fcu) * (1 - fhe)

            # This is needed right now. Can we change it later?
            j_crit_sc = max(1.0e-9, j_crit_sc)
            j_crit_cable = max(1.0e-9, j_crit_cable)
        elif i_tf_sc_mat == 4:  # As (1), but user-defined parameters
            bc20m = bcritsc
            tc0m = tcritsc
            j_crit_sc, bcrit, tcrit = superconductors.itersc(
                thelium, bmax, strain, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1 - fcu) * (1 - fhe)
        elif i_tf_sc_mat == 5:  # WST Nb3Sn parameterisation
            bc20m = 32.97
            tc0m = 16.06

            #  j_crit_sc returned by itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper

            j_crit_sc, bcrit, tcrit = superconductors.western_superconducting_nb3sn(
                thelium,
                bmax,
                strain,
                bc20m,
                tc0m,
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1 - fcu) * (1 - fhe)
        elif (
            i_tf_sc_mat == 6
        ):  # ! "REBCO" 2nd generation HTS superconductor in CrCo strand
            j_crit_sc, validity = superconductors.jcrit_rebco(thelium, bmax, 0)
            j_crit_sc = max(1.0e-9, j_crit_sc)
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1 - fcu) * (1 - fhe)

        elif i_tf_sc_mat == 7:  # Durham Ginzburg-Landau Nb-Ti parameterisation
            bc20m = b_crit_upper_nbti
            tc0m = t_crit_nbti
            j_crit_sc, bcrit, tcrit = superconductors.gl_nbti(
                thelium, bmax, strain, bc20m, tc0m
            )
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1 - fcu) * (1 - fhe)
        elif i_tf_sc_mat == 8:
            bc20m = 429
            tc0m = 185
            j_crit_sc, bcrit, tcrit = superconductors.gl_rebco(
                thelium, bmax, strain, bc20m, tc0m
            )
            # A0 calculated for tape cross section already
            # j_crit_cable = j_crit_sc * non-copper fraction of conductor * conductor fraction of cable
            j_crit_cable = j_crit_sc * (1 - fcu) * (1 - fhe)
        else:
            raise ProcessValueError(
                "Illegal value for i_pf_superconductor", i_tf_sc_mat=i_tf_sc_mat
            )

        return j_crit_sc * 1e-6

    def bmax_from_awp(
        self, wp_width_radial, current, n_tf_coils, r_coil_major, r_coil_minor
    ):
        """Returns a fitted function for bmax for stellarators

        author: J Lion, IPP Greifswald
        Returns a fitted function for bmax in dependece
        of the winding pack. The stellarator type config
        is taken from the parent scope.
        """

        return (
            2e-1
            * current
            * n_tf_coils
            / (r_coil_major - r_coil_minor)
            * (
                stellarator_configuration.stella_config_a1
                + stellarator_configuration.stella_config_a2
                * r_coil_major
                / wp_width_radial
            )
        )

    def intersect(self, x1, y1, x2, y2, xin):
        """Routine to find the x (abscissa) intersection point of two curves
        each defined by tabulated (x,y) values
        author: P J Knight, CCFE, Culham Science Centre
        x1(1:n1) : input real array : x values for first curve
        y1(1:n1) : input real array : y values for first curve
        n1       : input integer : length of arrays x1, y1
        x2(1:n2) : input real array : x values for first curve
        y2(1:n2) : input real array : y values for first curve
        n2       : input integer : length of arrays x2, y2
        x        : input/output real : initial x value guess on entry;
        x value at point of intersection on exit
        This routine estimates the x point (abscissa) at which two curves
        defined by tabulated (x,y) values intersect, using simple
        linear interpolation and the Newton-Raphson method.
        The routine will stop with an error message if no crossing point
        is found within the x ranges of the two curves.
        None
        """
        x = xin
        n1 = len(x1)
        n2 = len(x2)

        xmin = max(np.amin(x1), np.amin(x2))
        xmax = min(np.max(x1), np.amax(x2))

        if xmin >= xmax:
            logger.error(
                f"X ranges not overlapping. {np.amin(x1)=} {np.amin(x2)=} "
                f"{np.amax(x1)=} {np.amax(x2)=}"
            )

        #  Ensure input guess for x is within this range

        if x < xmin:
            x = xmin
        elif x > xmax:
            x = xmax

        #  Find overall y range, and set tolerance
        #  in final difference in y values

        ymin = min(np.amin(y1), np.amin(y2))
        ymax = max(np.max(y1), np.max(y2))

        epsy = 1.0e-6 * (ymax - ymin)

        #  Finite difference dx

        dx = 0.01e0 / max(n1, n2) * (xmax - xmin)

        for _i in range(100):
            #  Find difference in y values at x

            y01 = np.interp(x, x1, y1)
            y02 = np.interp(x, x2, y2)
            y = y01 - y02

            if abs(y) < epsy:
                break

            #  Find difference in y values at x+dx

            y01 = np.interp(x + dx, x1, y1)
            y02 = np.interp(x + dx, x2, y2)
            yright = y01 - y02

            #  Find difference in y values at x-dx

            y01 = np.interp(x - dx, x1, y1)
            y02 = np.interp(x - dx, x2, y2)
            yleft = y01 - y02

            #  Adjust x using Newton-Raphson method

            x = x - 2.0e0 * dx * y / (yright - yleft)

            if x < xmin:
                logger.error(
                    f"X has dropped below Xmin; X={x} has been set equal to Xmin={xmin}"
                )
                x = xmin
                break

            if x > xmax:
                logger.error(
                    f"X has risen above Xmax; X={x} has been set equal to Xmax={xmin}"
                )
                x = xmax
                break
        else:
            logger.error("Convergence too slow; X may be wrong...")

        return x

    def stopt_output(
        self,
        max_gyrotron_frequency,
        b_plasma_toroidal_on_axis,
        bt_ecrh,
        ne0_max_ECRH,
        te0_ecrh_achievable,
    ):
        po.oheadr(self.outfile, "ECRH Ignition at lower values. Information:")

        po.ovarre(
            self.outfile,
            "Maximal available gyrotron freq (input)",
            "(max_gyro_frequency)",
            max_gyrotron_frequency,
        )

        po.ovarre(
            self.outfile,
            "Operating point: bfield",
            "(b_plasma_toroidal_on_axis)",
            b_plasma_toroidal_on_axis,
        )
        po.ovarre(
            self.outfile,
            "Operating point: Peak density",
            "(nd_plasma_electron_on_axis)",
            physics_variables.nd_plasma_electron_on_axis,
        )
        po.ovarre(
            self.outfile,
            "Operating point: Peak temperature",
            "(temp_plasma_electron_on_axis_kev)",
            physics_variables.temp_plasma_electron_on_axis_kev,
        )

        po.ovarre(self.outfile, "Ignition point: bfield (T)", "(bt_ecrh)", bt_ecrh)
        po.ovarre(
            self.outfile,
            "Ignition point: density (/m3)",
            "(ne0_max_ECRH)",
            ne0_max_ECRH,
        )
        po.ovarre(
            self.outfile,
            "Maximum reachable ECRH temperature (pseudo) (KEV)",
            "(te0_ecrh_achievable)",
            te0_ecrh_achievable,
        )

        powerht_local, pscalingmw_local = self.power_at_ignition_point(
            max_gyrotron_frequency, te0_ecrh_achievable
        )
        po.ovarre(
            self.outfile,
            "Ignition point: Heating Power (MW)",
            "(powerht_ecrh)",
            powerht_local,
        )
        po.ovarre(
            self.outfile,
            "Ignition point: Loss Power (MW)",
            "(pscalingmw_ecrh)",
            pscalingmw_local,
        )

        if powerht_local >= pscalingmw_local:
            po.ovarin(self.outfile, "Operation point ECRH ignitable?", "(ecrh_bool)", 1)
        else:
            po.ovarin(self.outfile, "Operation point ECRH ignitable?", "(ecrh_bool)", 0)

    def power_at_ignition_point(self, gyro_frequency_max, te0_available):
        """Routine to calculate if the plasma is ignitable with the current values for the B field. Assumes
        current ECRH achievable peak temperature (which is inaccurate as the cordey pass should be calculated)
        author: J Lion, IPP Greifswald
        gyro_frequency_max : input real : Maximal available Gyrotron frequency (1/s) NOT (rad/s)
        te0_available : input real : Reachable peak electron temperature, reached by ECRH (KEV)
        powerht_out : output real: Heating Power at ignition point (MW)
        pscalingmw_out : output real: Heating Power loss at ignition point (MW)
        This routine calculates the density limit due to an ECRH heating scheme on axis
        Assumes current peak temperature (which is inaccurate as the cordey pass should be calculated)
        Maybe use this: https://doi.org/10.1088/0029-5515/49/8/085026
        """
        te_old = copy(physics_variables.temp_plasma_electron_vol_avg_kev)
        # Volume averaged physics_variables.te from te0_achievable
        physics_variables.temp_plasma_electron_vol_avg_kev = te0_available / (
            1.0e0 + physics_variables.alphat
        )
        ne0_max, bt_ecrh_max = self.stdlim_ecrh(
            gyro_frequency_max, physics_variables.b_plasma_toroidal_on_axis
        )
        # Now go to point where ECRH is still available
        # In density..
        dene_old = copy(physics_variables.nd_plasma_electrons_vol_avg)
        physics_variables.nd_plasma_electrons_vol_avg = min(
            dene_old, ne0_max / (1.0e0 + physics_variables.alphan)
        )

        # And B-field..
        bt_old = copy(physics_variables.b_plasma_toroidal_on_axis)
        physics_variables.b_plasma_toroidal_on_axis = min(
            bt_ecrh_max, physics_variables.b_plasma_toroidal_on_axis
        )

        self.stphys(False)
        self.stphys(
            False
        )  # The second call seems to be necessary for all values to "converge" (and is sufficient)

        powerht_out = max(
            copy(physics_variables.p_plasma_loss_mw), 0.00001e0
        )  # the radiation module sometimes returns negative heating power
        pscalingmw_out = copy(physics_variables.pscalingmw)

        # Reverse it and do it again because anything more efficiently isn't suitable with the current implementation
        # This is bad practice but seems to be necessary as of now:
        physics_variables.temp_plasma_electron_vol_avg_kev = te_old
        physics_variables.nd_plasma_electrons_vol_avg = dene_old
        physics_variables.b_plasma_toroidal_on_axis = bt_old

        self.stphys(False)
        self.stphys(False)

        return powerht_out, pscalingmw_out

    def stdlim(self, b_plasma_toroidal_on_axis, powht, rmajor, rminor):
        """Routine to calculate the Sudo density limit in a stellarator
        author: P J Knight, CCFE, Culham Science Centre
        b_plasma_toroidal_on_axis     : input real : Toroidal field on axis (T)
        powht  : input real : Absorbed heating power (MW)
        rmajor : input real : Plasma major radius (m)
        rminor : input real : Plasma minor radius (m)
        nd_plasma_electron_max_array : output real : Maximum volume-averaged plasma density (/m3)
        This routine calculates the density limit for a stellarator.
        S.Sudo, Y.Takeiri, H.Zushi et al., Scalings of Energy Confinement
        and Density Limit in Stellarator/Heliotron Devices, Nuclear Fusion
        vol.30, 11 (1990).
        """
        arg = powht * b_plasma_toroidal_on_axis / (rmajor * rminor * rminor)

        if arg <= 0.0e0:
            raise ProcessValueError(
                "Negative square root imminent",
                arg=arg,
                powht=powht,
                b_plasma_toroidal_on_axis=b_plasma_toroidal_on_axis,
                rmajor=rmajor,
                rminor=rminor,
            )

        #  Maximum line-averaged electron density

        dnlamx = 0.25e20 * np.sqrt(arg)

        #  Scale the result so that it applies to the volume-averaged
        #  electron density

        nd_plasma_electron_max_array = (
            dnlamx
            * physics_variables.nd_plasma_electrons_vol_avg
            / physics_variables.nd_plasma_electron_line
        )

        #  Set the required value for icc=5

        physics_variables.nd_plasma_electrons_max = nd_plasma_electron_max_array

        return nd_plasma_electron_max_array

    def stcoil_output(
        self,
        a_tf_wp_no_insulation,
        centering_force_avg_mn,
        centering_force_max_mn,
        centering_force_min_mn,
        coilcoilgap,
        coppera_m2,
        coppera_m2_max,
        f_scu,
        f_vv_actual,
        fiooic,
        inductance,
        max_force_density,
        max_force_density_mnm,
        max_lateral_force_density,
        max_radial_force_density,
        min_bending_radius,
        r_coil_major,
        r_coil_minor,
        r_tf_inleg_mid,
        sig_tf_wp,
        dx_tf_turn_general,
        t_tf_superconductor_quench,
        tf_total_h_width,
        tfborev,
        toroidalgap,
        v_tf_coil_dump_quench_max_kv,
        v_tf_coil_dump_quench_kv,
    ):
        """Writes stellarator modular coil output to file
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        This routine writes the stellarator modular coil results
        to the output file.
        None
        """
        po.oheadr(self.outfile, "Modular Coils")

        po.osubhd(self.outfile, "General Coil Parameters :")

        po.ovarre(
            self.outfile,
            "Number of modular coils",
            "(n_tf_coils)",
            tfcoil_variables.n_tf_coils,
        )
        po.ovarre(self.outfile, "Av. coil major radius", "(coil_r)", r_coil_major)
        po.ovarre(self.outfile, "Av. coil minor radius", "(coil_a)", r_coil_minor)
        po.ovarre(
            self.outfile,
            "Av. coil aspect ratio",
            "(coil_aspect)",
            r_coil_major / r_coil_minor,
        )

        po.ovarre(
            self.outfile,
            "Cross-sectional area per coil (m2)",
            "(tfarea/n_tf_coils)",
            tfcoil_variables.a_tf_inboard_total / tfcoil_variables.n_tf_coils,
        )
        po.ovarre(
            self.outfile,
            "Total inboard leg radial thickness (m)",
            "(dr_tf_inboard)",
            build_variables.dr_tf_inboard,
        )
        po.ovarre(
            self.outfile,
            "Total outboard leg radial thickness (m)",
            "(dr_tf_outboard)",
            build_variables.dr_tf_outboard,
        )
        po.ovarre(
            self.outfile,
            "Inboard leg outboard half-width (m)",
            "(tficrn)",
            tfcoil_variables.tficrn,
        )
        po.ovarre(
            self.outfile,
            "Inboard leg inboard half-width (m)",
            "(tfocrn)",
            tfcoil_variables.tfocrn,
        )
        po.ovarre(
            self.outfile,
            "Outboard leg toroidal thickness (m)",
            "(dx_tf_inboard_out_toroidal)",
            tfcoil_variables.dx_tf_inboard_out_toroidal,
        )
        po.ovarre(
            self.outfile, "Minimum coil distance (m)", "(toroidalgap)", toroidalgap
        )
        po.ovarre(
            self.outfile,
            "Minimal left gap between coils (m)",
            "(coilcoilgap)",
            coilcoilgap,
        )
        po.ovarre(
            self.outfile,
            "Minimum coil bending radius (m)",
            "(min_bend_radius)",
            min_bending_radius,
        )
        po.ovarre(
            self.outfile,
            "Mean coil circumference (m)",
            "(len_tf_coil)",
            tfcoil_variables.len_tf_coil,
        )
        po.ovarre(
            self.outfile,
            "Total current (MA)",
            "(c_tf_total)",
            1.0e-6 * tfcoil_variables.c_tf_total,
        )
        po.ovarre(
            self.outfile,
            "Current per coil(MA)",
            "(c_tf_total/n_tf_coils)",
            1.0e-6 * tfcoil_variables.c_tf_total / tfcoil_variables.n_tf_coils,
        )
        po.ovarre(
            self.outfile,
            "Winding pack current density (A/m2)",
            "(j_tf_wp)",
            tfcoil_variables.j_tf_wp,
        )
        po.ovarre(
            self.outfile,
            "Max allowable current density as restricted by quench (A/m2)",
            "(j_tf_wp_quench_heat_max)",
            tfcoil_variables.j_tf_wp_quench_heat_max,
        )
        po.ovarre(
            self.outfile,
            "Overall current density (A/m2)",
            "(oacdcp)",
            tfcoil_variables.oacdcp,
        )
        po.ovarre(
            self.outfile,
            "Maximum field on superconductor (T)",
            "(b_tf_inboard_peak_symmetric)",
            tfcoil_variables.b_tf_inboard_peak_symmetric,
        )
        po.ovarre(
            self.outfile,
            "Total Stored energy (GJ)",
            "(e_tf_magnetic_stored_total_gj)",
            tfcoil_variables.e_tf_magnetic_stored_total_gj,
        )
        po.ovarre(
            self.outfile, "Inductance of TF Coils (H)", "(inductance)", inductance
        )
        po.ovarre(
            self.outfile,
            "Total mass of coils (kg)",
            "(m_tf_coils_total)",
            tfcoil_variables.m_tf_coils_total,
        )

        po.osubhd(self.outfile, "Coil Geometry :")
        po.ovarre(
            self.outfile,
            "Inboard leg centre radius (m)",
            "(r_tf_inleg_mid)",
            r_tf_inleg_mid,
        )
        po.ovarre(
            self.outfile,
            "Outboard leg centre radius (m)",
            "(r_tf_outboard_mid)",
            build_variables.r_tf_outboard_mid,
        )
        po.ovarre(
            self.outfile,
            "Maximum inboard edge height (m)",
            "(z_tf_inside_half)",
            build_variables.z_tf_inside_half,
        )
        po.ovarre(
            self.outfile,
            "Clear horizontal dr_bore (m)",
            "(tf_total_h_width)",
            tf_total_h_width,
        )
        po.ovarre(self.outfile, "Clear vertical dr_bore (m)", "(tfborev)", tfborev)

        po.osubhd(self.outfile, "Conductor Information :")
        po.ovarre(
            self.outfile,
            "Superconductor mass per coil (kg)",
            "(m_tf_coil_superconductor)",
            tfcoil_variables.m_tf_coil_superconductor,
        )
        po.ovarre(
            self.outfile,
            "Copper mass per coil (kg)",
            "(m_tf_coil_copper)",
            tfcoil_variables.m_tf_coil_copper,
        )
        po.ovarre(
            self.outfile,
            "Steel conduit mass per coil (kg)",
            "(m_tf_wp_steel_conduit)",
            tfcoil_variables.m_tf_wp_steel_conduit,
        )
        po.ovarre(
            self.outfile,
            "Total conductor cable mass per coil (kg)",
            "(m_tf_coil_conductor)",
            tfcoil_variables.m_tf_coil_conductor,
        )
        po.ovarre(
            self.outfile,
            "Cable conductor + void area (m2)",
            "(a_tf_turn_cable_space_no_void)",
            tfcoil_variables.a_tf_turn_cable_space_no_void,
        )
        po.ovarre(
            self.outfile,
            "Cable space coolant fraction",
            "(f_a_tf_turn_cable_space_extra_void)",
            tfcoil_variables.f_a_tf_turn_cable_space_extra_void,
        )
        po.ovarre(
            self.outfile,
            "Conduit case thickness (m)",
            "(dx_tf_turn_steel)",
            tfcoil_variables.dx_tf_turn_steel,
        )
        po.ovarre(
            self.outfile,
            "Cable insulation thickness (m)",
            "(dx_tf_turn_insulation)",
            tfcoil_variables.dx_tf_turn_insulation,
        )

        ap = a_tf_wp_no_insulation
        po.osubhd(self.outfile, "Winding Pack Information :")
        po.ovarre(self.outfile, "Winding pack area", "(ap)", ap)
        po.ovarre(
            self.outfile,
            "Conductor fraction of winding pack",
            "(a_tf_wp_conductor/ap)",
            tfcoil_variables.a_tf_wp_conductor / ap,
        )
        po.ovarre(
            self.outfile,
            "Copper fraction of conductor",
            "(f_a_tf_turn_cable_copper)",
            tfcoil_variables.f_a_tf_turn_cable_copper,
        )
        po.ovarre(
            self.outfile,
            "Structure fraction of winding pack",
            "(a_tf_wp_steel/ap)",
            tfcoil_variables.a_tf_wp_steel / ap,
        )
        po.ovarre(
            self.outfile,
            "Insulator fraction of winding pack",
            "(a_tf_coil_wp_turn_insulation/ap)",
            tfcoil_variables.a_tf_coil_wp_turn_insulation / ap,
        )
        po.ovarre(
            self.outfile,
            "Helium fraction of winding pack",
            "(a_tf_wp_extra_void/ap)",
            tfcoil_variables.a_tf_wp_extra_void / ap,
        )
        po.ovarre(
            self.outfile,
            "Winding radial thickness (m)",
            "(dr_tf_wp_with_insulation)",
            tfcoil_variables.dr_tf_wp_with_insulation,
        )
        po.ovarre(
            self.outfile,
            "Winding toroidal thickness (m)",
            "(dx_tf_wp_primary_toroidal)",
            tfcoil_variables.dx_tf_wp_primary_toroidal,
        )
        po.ovarre(
            self.outfile,
            "Ground wall insulation thickness (m)",
            "(dx_tf_wp_insulation)",
            tfcoil_variables.dx_tf_wp_insulation,
        )
        po.ovarre(
            self.outfile,
            "Number of turns per coil",
            "(n_tf_coil_turns)",
            tfcoil_variables.n_tf_coil_turns,
        )
        po.ovarre(
            self.outfile,
            "Width of each turn (incl. insulation) (m)",
            "(dx_tf_turn_general)",
            dx_tf_turn_general,
        )
        po.ovarre(
            self.outfile,
            "Current per turn (A)",
            "(c_tf_turn)",
            tfcoil_variables.c_tf_turn,
        )
        po.ovarre(self.outfile, "jop/jcrit", "(fiooic)", fiooic)
        po.ovarre(
            self.outfile,
            "Current density in conductor area (A/m2)",
            "(c_tf_total/a_tf_wp_conductor)",
            1.0e-6
            * tfcoil_variables.c_tf_total
            / tfcoil_variables.n_tf_coils
            / tfcoil_variables.a_tf_wp_conductor,
        )
        po.ovarre(
            self.outfile,
            "Current density in SC area (A/m2)",
            "(c_tf_total/a_tf_wp_conductor/f_scu)",
            1.0e-6
            * tfcoil_variables.c_tf_total
            / tfcoil_variables.n_tf_coils
            / ap
            / f_scu,
        )
        po.ovarre(self.outfile, "Superconductor faction of WP (1)", "(f_scu)", f_scu)

        po.osubhd(self.outfile, "Forces and Stress :")
        po.ovarre(
            self.outfile,
            "Maximal toroidally and radially av. force density (MN/m3)",
            "(max_force_density)",
            max_force_density,
        )
        po.ovarre(
            self.outfile,
            "Maximal force density (MN/m)",
            "(max_force_density_Mnm)",
            max_force_density_mnm,
        )
        po.ovarre(
            self.outfile,
            "Maximal stress (approx.) (MPa)",
            "(sig_tf_wp)",
            sig_tf_wp * 1.0e-6,
        )

        po.ovarre(
            self.outfile,
            "Maximal lateral force density (MN/m3)",
            "(max_lateral_force_density)",
            max_lateral_force_density,
        )
        po.ovarre(
            self.outfile,
            "Maximal radial force density (MN/m3)",
            "(max_radial_force_density)",
            max_radial_force_density,
        )

        po.ovarre(
            self.outfile,
            "Max. centering force (coil) (MN)",
            "(centering_force_max_MN)",
            centering_force_max_mn,
        )
        po.ovarre(
            self.outfile,
            "Min. centering force (coil) (MN)",
            "(centering_force_min_MN)",
            centering_force_min_mn,
        )
        po.ovarre(
            self.outfile,
            "Avg. centering force per coil (MN)",
            "(centering_force_avg_MN)",
            centering_force_avg_mn,
        )

        po.osubhd(self.outfile, "Quench Restrictions :")
        po.ovarre(
            self.outfile,
            "Actual quench time (or time constant) (s)",
            "(t_tf_superconductor_quench)",
            t_tf_superconductor_quench,
        )
        po.ovarre(
            self.outfile,
            "Actual quench vaccuum vessel force density (MN/m^3)",
            "(f_vv_actual)",
            f_vv_actual,
        )
        po.ovarre(
            self.outfile,
            "Maximum allowed voltage during quench due to insulation (kV)",
            "(v_tf_coil_dump_quench_max_kv)",
            v_tf_coil_dump_quench_max_kv,
        )
        po.ovarre(
            self.outfile,
            "Actual quench voltage (kV)",
            "(v_tf_coil_dump_quench_kv)",
            v_tf_coil_dump_quench_kv,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Current (A) per mm^2 copper (A/mm2)",
            "(coppera_m2)",
            coppera_m2 * 1.0e-6,
        )
        po.ovarre(
            self.outfile,
            "Max Copper current fraction:",
            "(coppera_m2/coppera_m2_max)",
            coppera_m2 / coppera_m2_max,
        )

        po.osubhd(self.outfile, "External Case Information :")

        po.ovarre(
            self.outfile,
            "Case thickness, plasma side (m)",
            "(dr_tf_plasma_case)",
            tfcoil_variables.dr_tf_plasma_case,
        )
        po.ovarre(
            self.outfile,
            "Case thickness, outer side (m)",
            "(dr_tf_nose_case)",
            tfcoil_variables.dr_tf_nose_case,
        )
        po.ovarre(
            self.outfile,
            "Case toroidal thickness (m)",
            "(dx_tf_side_case_min)",
            tfcoil_variables.dx_tf_side_case_min,
        )
        po.ovarre(
            self.outfile,
            "Case area per coil (m2)",
            "(a_tf_coil_inboard_case)",
            tfcoil_variables.a_tf_coil_inboard_case,
        )
        po.ovarre(
            self.outfile,
            "External case mass per coil (kg)",
            "(m_tf_coil_case)",
            tfcoil_variables.m_tf_coil_case,
        )

        po.osubhd(self.outfile, "Available Space for Ports :")

        po.ovarre(
            self.outfile,
            "Max toroidal size of vertical ports (m)",
            "(vporttmax)",
            stellarator_variables.vporttmax,
        )
        po.ovarre(
            self.outfile,
            "Max poloidal size of vertical ports (m)",
            "(vportpmax)",
            stellarator_variables.vportpmax,
        )
        po.ovarre(
            self.outfile,
            "Max area of vertical ports (m2)",
            "(vportamax)",
            stellarator_variables.vportamax,
        )
        po.ovarre(
            self.outfile,
            "Max toroidal size of horizontal ports (m)",
            "(hporttmax)",
            stellarator_variables.hporttmax,
        )
        po.ovarre(
            self.outfile,
            "Max poloidal size of horizontal ports (m)",
            "(hportpmax)",
            stellarator_variables.hportpmax,
        )
        po.ovarre(
            self.outfile,
            "Max area of horizontal ports (m2)",
            "(hportamax)",
            stellarator_variables.hportamax,
        )

    def stdlim_ecrh(self, gyro_frequency_max, bt_input):
        """Routine to calculate the density limit due to an ECRH heating scheme on axis
        depending on an assumed maximal available gyrotron frequency.
        author: J Lion, IPP Greifswald
        gyro_frequency_max     : input real : Maximal available Gyrotron frequency (1/s) NOT (rad/s)
        b_plasma_toroidal_on_axis  : input real : Maximal magnetic field on axis (T)
        dlimit_ecrh : output real : Maximum peak plasma density by ECRH constraints (/m3)
        bt_max : output real : Maximum allowable b field for ecrh heating (T)
        This routine calculates the density limit due to an ECRH heating scheme on axis
        """
        gyro_frequency = min(1.76e11 * bt_input, gyro_frequency_max * 2.0e0 * np.pi)

        # Restrict b field to the maximal available gyrotron frequency
        bt_max = (gyro_frequency_max * 2.0e0 * np.pi) / 1.76e11

        #                      me*e0/e^2       * w^2
        ne0_max = max(0.0e0, 3.142077e-4 * gyro_frequency**2)

        # Check if parabolic profiles are used:
        if physics_variables.i_plasma_pedestal == 0:
            # Parabolic profiles used, use analytical formula:
            dlimit_ecrh = ne0_max
        else:
            logger.error(
                "It was used physics_variables.i_plasma_pedestal = 1 in a stellarator routine. PROCESS will pretend it got parabolic profiles (physics_variables.i_plasma_pedestal = 0)."
            )
            dlimit_ecrh = ne0_max

        return dlimit_ecrh, bt_max

    def stphys(self, output):
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

        #  Poloidal physics_variables.beta

        #  Perform auxiliary power calculations

        self.stheat(False)

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
                physics_variables.zeffai,
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

        physics_variables.beta_fast_alpha = physics_funcs.fast_alpha_beta(
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
                    / build_variables.a_fw_total
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
                    / build_variables.a_fw_total
                )

        #  Calculate ion/electron equilibration power

        physics_variables.pden_ion_electron_equilibration_mw = rether(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.dlamie,
            physics_variables.temp_plasma_electron_vol_avg_kev,
            physics_variables.temp_plasma_ion_vol_avg_kev,
            physics_variables.zeffai,
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
        physics_variables.pden_plasma_outer_rad_mw = (
            radpwr_data.pden_plasma_outer_rad_mw
        )
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
        physics_variables.p_plasma_rad_mw = max(
            0.0e0, physics_variables.p_plasma_rad_mw
        )

        #  Power to divertor, = (1-stellarator_variables.f_rad)*Psol

        # The SOL radiation needs to be smaller than the physics_variables.p_plasma_rad_mw
        physics_variables.psolradmw = stellarator_variables.f_rad * powht
        physics_variables.p_plasma_separatrix_mw = powht - physics_variables.psolradmw

        # Add SOL Radiation to total
        physics_variables.p_plasma_rad_mw = (
            physics_variables.p_plasma_rad_mw + physics_variables.psolradmw
        )
        # pden_plasma_rad_mw = physics_variables.p_plasma_rad_mw / physics_variables.vol_plasma # this line OVERWRITES the original definition of pden_plasma_rad_mw, probably shouldn't be defined like that as the core does not lose SOL power.

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
                    / build_variables.a_fw_total
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
                    / build_variables.a_fw_total
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
            fusrat,
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

        # Calculate physics_variables.beta limit. Does nothing atm so commented out
        # call stblim(physics_variables.beta_vol_avg_max)

        # Calculate the neoclassical sanity check with PROCESS parameters
        (
            q_PROCESS,
            q_PROCESS_r1,
            q_neo,
            gamma_neo,
            total_q_neo,
            total_q_neo_e,
            q_neo_e,
            q_neo_D,
            q_neo_a,
            q_neo_T,
            g_neo_e,
            g_neo_D,
            g_neo_a,
            g_neo_T,
            dndt_neo_e,
            dndt_neo_D,
            dndt_neo_a,
            dndt_neo_T,
            dndt_neo_fuel,
            dmdt_neo_fuel,
            dmdt_neo_fuel_from_e,
            chi_neo_e,
            chi_PROCESS_e,
            nu_star_e,
            nu_star_d,
            nu_star_T,
            nu_star_He,
        ) = self.calc_neoclassics()

        if output:
            self.stphys_output(
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

    def stphys_output(
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

    def calc_neoclassics(self):
        self.neoclassics.init_neoclassics(
            0.6,
            stellarator_configuration.stella_config_epseff,
            stellarator_variables.iotabar,
        )

        q_PROCESS = (
            (
                physics_variables.f_p_alpha_plasma_deposited
                * physics_variables.pden_alpha_total_mw
                - physics_variables.pden_plasma_core_rad_mw
            )
            * physics_variables.vol_plasma
            / physics_variables.a_plasma_surface
            * impurity_radiation_module.radius_plasma_core_norm
        )
        q_PROCESS_r1 = (
            (
                physics_variables.f_p_alpha_plasma_deposited
                * physics_variables.pden_alpha_total_mw
                - physics_variables.pden_plasma_core_rad_mw
            )
            * physics_variables.vol_plasma
            / physics_variables.a_plasma_surface
        )

        q_neo = sum(neoclassics_variables.q_flux * 1e-6)
        gamma_neo = sum(
            neoclassics_variables.gamma_flux * neoclassics_variables.temperatures * 1e-6
        )

        total_q_neo = sum(
            neoclassics_variables.q_flux * 1e-6
            + neoclassics_variables.gamma_flux
            * neoclassics_variables.temperatures
            * 1e-6
        )

        total_q_neo_e = (
            2
            * 2
            * (
                neoclassics_variables.q_flux[0] * 1e-6
                + neoclassics_variables.gamma_flux[0]
                * neoclassics_variables.temperatures[0]
                * 1e-6
            )
        )

        q_neo_e = neoclassics_variables.q_flux[0] * 1e-6
        q_neo_D = neoclassics_variables.q_flux[1] * 1e-6
        q_neo_a = neoclassics_variables.q_flux[3] * 1e-6
        q_neo_T = neoclassics_variables.q_flux[2] * 1e-6

        g_neo_e = (
            neoclassics_variables.gamma_flux[0]
            * 1e-6
            * neoclassics_variables.temperatures[0]
        )
        g_neo_D = (
            neoclassics_variables.gamma_flux[1]
            * 1e-6
            * neoclassics_variables.temperatures[1]
        )
        g_neo_a = (
            neoclassics_variables.gamma_flux[3]
            * 1e-6
            * neoclassics_variables.temperatures[3]
        )
        g_neo_T = (
            neoclassics_variables.gamma_flux[2]
            * 1e-6
            * neoclassics_variables.temperatures[2]
        )

        dndt_neo_e = neoclassics_variables.gamma_flux[0]
        dndt_neo_D = neoclassics_variables.gamma_flux[1]
        dndt_neo_a = neoclassics_variables.gamma_flux[3]
        dndt_neo_T = neoclassics_variables.gamma_flux[2]

        dndt_neo_fuel = (
            (dndt_neo_D + dndt_neo_T)
            * physics_variables.a_plasma_surface
            * impurity_radiation_module.radius_plasma_core_norm
        )
        dmdt_neo_fuel = (
            dndt_neo_fuel * physics_variables.m_fuel_amu * constants.PROTON_MASS * 1.0e6
        )  # mg
        dmdt_neo_fuel_from_e = (
            4
            * dndt_neo_e
            * physics_variables.a_plasma_surface
            * impurity_radiation_module.radius_plasma_core_norm
            * physics_variables.m_fuel_amu
            * constants.PROTON_MASS
            * 1.0e6
        )  # kg

        chi_neo_e = -(
            neoclassics_variables.q_flux[0]
            + neoclassics_variables.gamma_flux[0]
            * neoclassics_variables.temperatures[0]
        ) / (
            neoclassics_variables.densities[0]
            * neoclassics_variables.dr_temperatures[0]
            + neoclassics_variables.temperatures[0]
            * neoclassics_variables.dr_densities[0]
        )

        chi_PROCESS_e = self.st_calc_eff_chi()

        nu_star_e = neoclassics_variables.nu_star_averaged[0]
        nu_star_d = neoclassics_variables.nu_star_averaged[1]
        nu_star_T = neoclassics_variables.nu_star_averaged[2]
        nu_star_He = neoclassics_variables.nu_star_averaged[3]

        return (
            q_PROCESS,
            q_PROCESS_r1,
            q_neo,
            gamma_neo,
            total_q_neo,
            total_q_neo_e,
            q_neo_e,
            q_neo_D,
            q_neo_a,
            q_neo_T,
            g_neo_e,
            g_neo_D,
            g_neo_a,
            g_neo_T,
            dndt_neo_e,
            dndt_neo_D,
            dndt_neo_a,
            dndt_neo_T,
            dndt_neo_fuel,
            dmdt_neo_fuel,
            dmdt_neo_fuel_from_e,
            chi_neo_e,
            chi_PROCESS_e,
            nu_star_e,
            nu_star_d,
            nu_star_T,
            nu_star_He,
        )

    def st_calc_eff_chi(self):
        volscaling = (
            physics_variables.vol_plasma
            * stellarator_variables.f_r
            * (
                impurity_radiation_module.radius_plasma_core_norm
                * physics_variables.rminor
                / stellarator_configuration.stella_config_rminor_ref
            )
            ** 2
        )
        surfacescaling = (
            physics_variables.a_plasma_surface
            * stellarator_variables.f_r
            * (
                impurity_radiation_module.radius_plasma_core_norm
                * physics_variables.rminor
                / stellarator_configuration.stella_config_rminor_ref
            )
        )

        nominator = (
            physics_variables.f_p_alpha_plasma_deposited
            * physics_variables.pden_alpha_total_mw
            - physics_variables.pden_plasma_core_rad_mw
        ) * volscaling

        # in fortran there was a 0*alphan term which I have removed for obvious reasons
        # the following comment seems to describe this?
        # "include alphan if chi should be incorporate density gradients too"
        # but the history can be consulted if required (23/11/22 TN)
        denominator = (
            (
                3
                * physics_variables.nd_plasma_electron_on_axis
                * constants.ELECTRON_CHARGE
                * physics_variables.temp_plasma_electron_on_axis_kev
                * 1e3
                * physics_variables.alphat
                * impurity_radiation_module.radius_plasma_core_norm
                * (1 - impurity_radiation_module.radius_plasma_core_norm**2)
                ** (physics_variables.alphan + physics_variables.alphat - 1)
            )
            * surfacescaling
            * 1e-6
        )

        return nominator / denominator

    def stheat(self, output: bool):
        """Routine to calculate the auxiliary heating power
        in a stellarator
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This routine calculates the auxiliary heating power for
        a stellarator device.
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        if stellarator_variables.isthtr == 1:
            current_drive_variables.p_hcd_ecrh_injected_total_mw = (
                current_drive_variables.p_hcd_primary_extra_heat_mw
            )
            current_drive_variables.p_hcd_injected_ions_mw = 0
            current_drive_variables.p_hcd_injected_electrons_mw = (
                current_drive_variables.p_hcd_ecrh_injected_total_mw
            )
            current_drive_variables.eta_hcd_primary_injector_wall_plug = (
                current_drive_variables.eta_ecrh_injector_wall_plug
            )
            current_drive_variables.p_hcd_electric_total_mw = (
                current_drive_variables.p_hcd_injected_ions_mw
                + current_drive_variables.p_hcd_injected_electrons_mw
            ) / current_drive_variables.eta_hcd_primary_injector_wall_plug
        elif stellarator_variables.isthtr == 2:
            current_drive_variables.p_hcd_lowhyb_injected_total_mw = (
                current_drive_variables.p_hcd_primary_extra_heat_mw
            )
            current_drive_variables.p_hcd_injected_ions_mw = 0
            current_drive_variables.p_hcd_injected_electrons_mw = (
                current_drive_variables.p_hcd_lowhyb_injected_total_mw
            )
            current_drive_variables.eta_hcd_primary_injector_wall_plug = (
                current_drive_variables.eta_lowhyb_injector_wall_plug
            )
            current_drive_variables.p_hcd_electric_total_mw = (
                current_drive_variables.p_hcd_injected_ions_mw
                + current_drive_variables.p_hcd_injected_electrons_mw
            ) / current_drive_variables.eta_hcd_primary_injector_wall_plug
        elif stellarator_variables.isthtr == 3:
            (
                effnbss,
                f_p_beam_injected_ions,
                current_drive_variables.f_p_beam_shine_through,
            ) = self.current_drive.culnbi()
            current_drive_variables.p_hcd_beam_injected_total_mw = (
                current_drive_variables.p_hcd_primary_extra_heat_mw
                * (1 - current_drive_variables.f_p_beam_orbit_loss)
            )
            current_drive_variables.p_beam_orbit_loss_mw = (
                current_drive_variables.p_hcd_primary_extra_heat_mw
                * current_drive_variables.f_p_beam_orbit_loss
            )
            current_drive_variables.p_hcd_injected_ions_mw = (
                current_drive_variables.p_hcd_beam_injected_total_mw
                * f_p_beam_injected_ions
            )
            current_drive_variables.p_hcd_injected_electrons_mw = (
                current_drive_variables.p_hcd_beam_injected_total_mw
                * (1 - f_p_beam_injected_ions)
            )
            current_drive_variables.eta_hcd_primary_injector_wall_plug = (
                current_drive_variables.eta_beam_injector_wall_plug
            )
            current_drive_variables.p_hcd_electric_total_mw = (
                current_drive_variables.p_hcd_injected_ions_mw
                + current_drive_variables.p_hcd_injected_electrons_mw
            ) / current_drive_variables.eta_hcd_primary_injector_wall_plug
        else:
            raise ProcessValueError(
                "Illegal value for isthtr", isthtr=stellarator_variables.isthtr
            )

        #  Total injected power

        current_drive_variables.p_hcd_injected_total_mw = (
            current_drive_variables.p_hcd_injected_electrons_mw
            + current_drive_variables.p_hcd_injected_ions_mw
        )

        #  Calculate neutral beam current

        if abs(current_drive_variables.p_hcd_beam_injected_total_mw) > 1e-8:
            current_drive_variables.c_beam_total = (
                1e-3
                * (current_drive_variables.p_hcd_beam_injected_total_mw * 1e6)
                / current_drive_variables.e_beam_kev
            )
        else:
            current_drive_variables.c_beam_total = 0

        #  Ratio of fusion to input (injection+ohmic) power

        if (
            abs(
                current_drive_variables.p_hcd_injected_total_mw
                + current_drive_variables.p_beam_orbit_loss_mw
                + physics_variables.p_plasma_ohmic_mw
            )
            < 1e-6
        ):
            current_drive_variables.big_q_plasma = 1e18
        else:
            current_drive_variables.big_q_plasma = (
                physics_variables.p_fusion_total_mw
                / (
                    current_drive_variables.p_hcd_injected_total_mw
                    + current_drive_variables.p_beam_orbit_loss_mw
                    + physics_variables.p_plasma_ohmic_mw
                )
            )

        if output:
            po.oheadr(self.outfile, "Auxiliary Heating System")

            if stellarator_variables.isthtr == 1:
                po.ocmmnt(self.outfile, "Electron Cyclotron Resonance Heating")
            elif stellarator_variables.isthtr == 2:
                po.ocmmnt(self.outfile, "Lower Hybrid Heating")
            elif stellarator_variables.isthtr == 3:
                po.ocmmnt(self.outfile, "Neutral Beam Injection Heating")

            if physics_variables.i_plasma_ignited == 1:
                po.ocmmnt(
                    self.outfile,
                    "Ignited plasma; injected power only used for start-up phase",
                )

            po.oblnkl(self.outfile)

            po.ovarre(
                self.outfile,
                "Auxiliary power supplied to plasma (MW)",
                "(p_hcd_primary_extra_heat_mw)",
                current_drive_variables.p_hcd_primary_extra_heat_mw,
            )
            po.ovarre(
                self.outfile,
                "Fusion gain factor Q",
                "(big_q_plasma)",
                current_drive_variables.big_q_plasma,
            )

            if abs(current_drive_variables.p_hcd_beam_injected_total_mw) > 1e-8:
                po.ovarre(
                    self.outfile,
                    "Neutral beam energy (KEV)",
                    "(enbeam)",
                    current_drive_variables.enbeam,
                )
                po.ovarre(
                    self.outfile,
                    "Neutral beam current (A)",
                    "(c_beam_total)",
                    current_drive_variables.c_beam_total,
                )
                po.ovarre(
                    self.outfile,
                    "Fraction of beam energy to ions",
                    "(f_p_beam_injected_ions)",
                    f_p_beam_injected_ions,
                )
                po.ovarre(
                    self.outfile,
                    "Neutral beam shine-through fraction",
                    "(f_p_beam_shine_through)",
                    current_drive_variables.f_p_beam_shine_through,
                )
                po.ovarre(
                    self.outfile,
                    "Neutral beam orbit loss power (MW)",
                    "(p_beam_orbit_loss_mw)",
                    current_drive_variables.p_beam_orbit_loss_mw,
                )
                po.ovarre(
                    self.outfile,
                    "Beam duct shielding thickness (m)",
                    "(dx_beam_shield)",
                    current_drive_variables.dx_beam_shield,
                )
                po.ovarre(
                    self.outfile,
                    "R injection tangent / R-major",
                    "(f_radius_beam_tangency_rmajor)",
                    current_drive_variables.f_radius_beam_tangency_rmajor,
                )
                po.ovarre(
                    self.outfile,
                    "Beam centreline tangency radius (m)",
                    "(radius_beam_tangency)",
                    current_drive_variables.radius_beam_tangency,
                )
                po.ovarre(
                    self.outfile,
                    "Maximum possible tangency radius (m)",
                    "(radius_beam_tangency_max)",
                    current_drive_variables.radius_beam_tangency_max,
                )
                po.ovarre(
                    self.outfile,
                    "Beam decay lengths to centre",
                    "(n_beam_decay_lengths_core)",
                    current_drive_variables.n_beam_decay_lengths_core,
                )


class Neoclassics:
    @property
    def no_roots(self):
        return neoclassics_variables.roots.shape[0]

    def init_neoclassics(self, r_effin, eps_effin, iotain):
        """Constructor of the neoclassics object from the effective radius,
        epsilon effective and iota only.
        """
        (
            neoclassics_variables.densities,
            neoclassics_variables.temperatures,
            neoclassics_variables.dr_densities,
            neoclassics_variables.dr_temperatures,
        ) = self.init_profile_values_from_PROCESS(r_effin)
        neoclassics_variables.roots = np.array([
            4.740718054080526184e-2,
            2.499239167531593919e-1,
            6.148334543927683749e-1,
            1.143195825666101451,
            1.836454554622572344,
            2.696521874557216147,
            3.725814507779509288,
            4.927293765849881879,
            6.304515590965073635,
            7.861693293370260349,
            9.603775985479263255,
            1.153654659795613924e1,
            1.366674469306423489e1,
            1.600222118898106771e1,
            1.855213484014315029e1,
            2.132720432178312819e1,
            2.434003576453269346e1,
            2.760555479678096091e1,
            3.114158670111123683e1,
            3.496965200824907072e1,
            3.911608494906788991e1,
            4.361365290848483056e1,
            4.850398616380419980e1,
            5.384138540650750571e1,
            5.969912185923549686e1,
            6.618061779443848991e1,
            7.344123859555988076e1,
            8.173681050672767867e1,
            9.155646652253683726e1,
            1.041575244310588886e2,
        ])
        neoclassics_variables.weights = np.array([
            1.160440860204388913e-1,
            2.208511247506771413e-1,
            2.413998275878537214e-1,
            1.946367684464170855e-1,
            1.237284159668764899e-1,
            6.367878036898660943e-2,
            2.686047527337972682e-2,
            9.338070881603925677e-3,
            2.680696891336819664e-3,
            6.351291219408556439e-4,
            1.239074599068830081e-4,
            1.982878843895233056e-5,
            2.589350929131392509e-6,
            2.740942840536013206e-7,
            2.332831165025738197e-8,
            1.580745574778327984e-9,
            8.427479123056716393e-11,
            3.485161234907855443e-12,
            1.099018059753451500e-13,
            2.588312664959080167e-15,
            4.437838059840028968e-17,
            5.365918308212045344e-19,
            4.393946892291604451e-21,
            2.311409794388543236e-23,
            7.274588498292248063e-26,
            1.239149701448267877e-28,
            9.832375083105887477e-32,
            2.842323553402700938e-35,
            1.878608031749515392e-39,
            8.745980440465011553e-45,
        ])

        neoclassics_variables.kt = self.neoclassics_calc_KT()
        neoclassics_variables.nu = self.neoclassics_calc_nu()
        neoclassics_variables.nu_star = self.neoclassics_calc_nu_star()
        neoclassics_variables.nu_star_averaged = self.neoclassics_calc_nu_star_fromT(
            iotain
        )
        neoclassics_variables.vd = self.neoclassics_calc_vd()

        neoclassics_variables.d11_plateau = self.neoclassics_calc_D11_plateau()

        neoclassics_variables.d11_mono = self.neoclassics_calc_d11_mono(
            eps_effin
        )  # for using epseff

        neoclassics_variables.d111 = self.calc_integrated_radial_transport_coeffs(
            index=1
        )
        neoclassics_variables.d112 = self.calc_integrated_radial_transport_coeffs(
            index=2
        )
        neoclassics_variables.d113 = self.calc_integrated_radial_transport_coeffs(
            index=3
        )

        neoclassics_variables.gamma_flux = self.neoclassics_calc_gamma_flux(
            neoclassics_variables.densities,
            neoclassics_variables.temperatures,
            neoclassics_variables.dr_densities,
            neoclassics_variables.dr_temperatures,
        )
        neoclassics_variables.q_flux = self.neoclassics_calc_q_flux()

    def init_profile_values_from_PROCESS(self, rho):
        """Initializes the profile_values object from PROCESS' parabolic profiles"""
        tempe = (
            physics_variables.temp_plasma_electron_on_axis_kev
            * (1 - rho**2) ** physics_variables.alphat
            * KEV
        )
        tempT = (
            physics_variables.temp_plasma_ion_on_axis_kev
            * (1 - rho**2) ** physics_variables.alphat
            * KEV
        )
        tempD = (
            physics_variables.temp_plasma_ion_on_axis_kev
            * (1 - rho**2) ** physics_variables.alphat
            * KEV
        )
        tempa = (
            physics_variables.temp_plasma_ion_on_axis_kev
            * (1 - rho**2) ** physics_variables.alphat
            * KEV
        )

        dense = (
            physics_variables.nd_plasma_electron_on_axis
            * (1 - rho**2) ** physics_variables.alphan
        )
        densT = (
            (1 - physics_variables.f_plasma_fuel_deuterium)
            * physics_variables.nd_plasma_ions_on_axis
            * (1 - rho**2) ** physics_variables.alphan
        )
        densD = (
            physics_variables.f_plasma_fuel_deuterium
            * physics_variables.nd_plasma_ions_on_axis
            * (1 - rho**2) ** physics_variables.alphan
        )
        densa = (
            physics_variables.nd_plasma_alphas_vol_avg
            * (1 + physics_variables.alphan)
            * (1 - rho**2) ** physics_variables.alphan
        )

        # Derivatives in real space
        dr_tempe = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.temp_plasma_electron_on_axis_kev
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )
        dr_tempT = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.temp_plasma_ion_on_axis_kev
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )
        dr_tempD = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.temp_plasma_ion_on_axis_kev
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )
        dr_tempa = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.temp_plasma_ion_on_axis_kev
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )

        dr_dense = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * physics_variables.nd_plasma_electron_on_axis
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )
        dr_densT = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * (1 - physics_variables.f_plasma_fuel_deuterium)
            * physics_variables.nd_plasma_ions_on_axis
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )
        dr_densD = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * physics_variables.f_plasma_fuel_deuterium
            * physics_variables.nd_plasma_ions_on_axis
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )
        dr_densa = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * physics_variables.nd_plasma_alphas_vol_avg
            * (1 + physics_variables.alphan)
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )

        dens = np.array([dense, densD, densT, densa])
        temp = np.array([tempe, tempD, tempT, tempa])
        dr_dens = np.array([dr_dense, dr_densD, dr_densT, dr_densa])
        dr_temp = np.array([dr_tempe, dr_tempD, dr_tempT, dr_tempa])

        return dens, temp, dr_dens, dr_temp

    def neoclassics_calc_KT(self):
        """Calculates the energy on the given grid
        which is given by the gauss laguerre roots.
        """
        k = np.repeat((neoclassics_variables.roots / KEV)[:, np.newaxis], 4, axis=1)

        return (k * neoclassics_variables.temperatures).T

    def neoclassics_calc_nu(self):
        """Calculates the collision frequency"""
        mass = np.array([
            constants.ELECTRON_MASS,
            constants.PROTON_MASS * 2.0,
            constants.PROTON_MASS * 3.0,
            constants.PROTON_MASS * 4.0,
        ])
        z = np.array([-1.0, 1.0, 1.0, 2.0]) * constants.ELECTRON_CHARGE

        # transform the temperature back in eV
        # Formula from L. Spitzer.Physics of fully ionized gases.  Interscience, New York, 1962
        lnlambda = (
            32.2
            - 1.15 * np.log10(neoclassics_variables.densities[0])
            + 2.3
            * np.log10(
                neoclassics_variables.temperatures[0] / constants.ELECTRON_CHARGE
            )
        )

        neoclassics_calc_nu = np.zeros((4, self.no_roots), order="F")

        for j in range(4):
            for i in range(self.no_roots):
                x = neoclassics_variables.roots[i]
                for k in range(4):
                    xk = (
                        (mass[k] / mass[j])
                        * (
                            neoclassics_variables.temperatures[j]
                            / neoclassics_variables.temperatures[k]
                        )
                        * x
                    )
                    expxk = np.exp(-xk)
                    t = 1.0 / (1.0 + 0.3275911 * np.sqrt(xk))
                    erfn = (
                        1.0
                        - t
                        * (
                            0.254829592
                            + t
                            * (
                                -0.284496736
                                + t
                                * (1.421413741 + t * (-1.453152027 + t * 1.061405429))
                            )
                        )
                        * expxk
                    )
                    phixmgx = (1.0 - 0.5 / xk) * erfn + expxk / np.sqrt(np.pi * xk)
                    v = np.sqrt(
                        2.0 * x * neoclassics_variables.temperatures[j] / mass[j]
                    )
                    neoclassics_calc_nu[j, i] = neoclassics_calc_nu[
                        j, i
                    ] + neoclassics_variables.densities[k] * (
                        z[j] * z[k]
                    ) ** 2 * lnlambda * phixmgx / (
                        4.0 * np.pi * constants.EPSILON0**2 * mass[j] ** 2 * v**3
                    )

        return neoclassics_calc_nu

    def neoclassics_calc_nu_star(self):
        """Calculates the normalized collision frequency"""
        k = np.repeat(neoclassics_variables.roots[:, np.newaxis], 4, axis=1)
        kk = (k * neoclassics_variables.temperatures).T

        mass = np.array([
            constants.ELECTRON_MASS,
            constants.PROTON_MASS * 2.0,
            constants.PROTON_MASS * 3.0,
            constants.PROTON_MASS * 4.0,
        ])

        v = np.empty((4, self.no_roots))
        v[0, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0 - (kk[0, :] / (mass[0] * constants.SPEED_LIGHT**2) + 1) ** (-1)
        )
        v[1, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0 - (kk[1, :] / (mass[1] * constants.SPEED_LIGHT**2) + 1) ** (-1)
        )
        v[2, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0 - (kk[2, :] / (mass[2] * constants.SPEED_LIGHT**2) + 1) ** (-1)
        )
        v[3, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0 - (kk[3, :] / (mass[3] * constants.SPEED_LIGHT**2) + 1) ** (-1)
        )

        return (
            physics_variables.rmajor
            * neoclassics_variables.nu
            / (neoclassics_variables.iota * v)
        )

    def neoclassics_calc_nu_star_fromT(self, iota):
        """Calculates the collision frequency"""
        temp = (
            np.array([
                physics_variables.temp_plasma_electron_vol_avg_kev,
                physics_variables.temp_plasma_ion_vol_avg_kev,
                physics_variables.temp_plasma_ion_vol_avg_kev,
                physics_variables.temp_plasma_ion_vol_avg_kev,
            ])
            * KEV
        )
        density = np.array([
            physics_variables.nd_plasma_electrons_vol_avg,
            physics_variables.nd_plasma_fuel_ions_vol_avg
            * physics_variables.f_plasma_fuel_deuterium,
            physics_variables.nd_plasma_fuel_ions_vol_avg
            * (1 - physics_variables.f_plasma_fuel_deuterium),
            physics_variables.nd_plasma_alphas_vol_avg,
        ])

        mass = np.array([
            constants.ELECTRON_MASS,
            constants.PROTON_MASS * 2.0,
            constants.PROTON_MASS * 3.0,
            constants.PROTON_MASS * 4.0,
        ])
        z = np.array([-1.0, 1.0, 1.0, 2.0]) * constants.ELECTRON_CHARGE

        # transform the temperature back in eV
        # Formula from L. Spitzer.Physics of fully ionized gases.  Interscience, New York, 1962
        lnlambda = (
            32.2
            - 1.15 * np.log10(density[0])
            + 2.3 * np.log10(temp[0] / constants.ELECTRON_CHARGE)
        )

        neoclassics_calc_nu_star_fromT = np.zeros((4,))

        for j in range(4):
            v = np.sqrt(2.0 * temp[j] / mass[j])
            for k in range(4):
                xk = (mass[k] / mass[j]) * (temp[j] / temp[k])

                expxk = 0.0
                if xk < 200.0:
                    expxk = np.exp(-xk)

                t = 1.0 / (1.0 + 0.3275911 * np.sqrt(xk))
                erfn = (
                    1.0
                    - t
                    * (
                        0.254829592
                        + t
                        * (
                            -0.284496736
                            + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))
                        )
                    )
                    * expxk
                )
                phixmgx = (1.0 - 0.5 / xk) * erfn + expxk / np.sqrt(np.pi * xk)
                neoclassics_calc_nu_star_fromT[j] = (
                    neoclassics_calc_nu_star_fromT[j]
                    + density[k]
                    * (z[j] * z[k]) ** 2
                    * lnlambda
                    * phixmgx
                    / (4.0 * np.pi * constants.EPSILON0**2 * mass[j] ** 2 * v**4)
                    * physics_variables.rmajor
                    / iota
                )
        return neoclassics_calc_nu_star_fromT

    def neoclassics_calc_vd(self):
        vde = (
            neoclassics_variables.roots
            * neoclassics_variables.temperatures[0]
            / (
                constants.ELECTRON_CHARGE
                * physics_variables.rmajor
                * physics_variables.b_plasma_toroidal_on_axis
            )
        )
        vdD = (
            neoclassics_variables.roots
            * neoclassics_variables.temperatures[1]
            / (
                constants.ELECTRON_CHARGE
                * physics_variables.rmajor
                * physics_variables.b_plasma_toroidal_on_axis
            )
        )
        vdT = (
            neoclassics_variables.roots
            * neoclassics_variables.temperatures[2]
            / (
                constants.ELECTRON_CHARGE
                * physics_variables.rmajor
                * physics_variables.b_plasma_toroidal_on_axis
            )
        )
        vda = (
            neoclassics_variables.roots
            * neoclassics_variables.temperatures[3]
            / (
                2.0
                * constants.ELECTRON_CHARGE
                * physics_variables.rmajor
                * physics_variables.b_plasma_toroidal_on_axis
            )
        )

        vd = np.empty((4, self.no_roots))

        vd[0, :] = vde
        vd[1, :] = vdD
        vd[2, :] = vdT
        vd[3, :] = vda

        return vd

    def neoclassics_calc_D11_plateau(self):
        """Calculates the plateau transport coefficients (D11_star sometimes)"""
        mass = np.array([
            constants.ELECTRON_MASS,
            constants.PROTON_MASS * 2.0,
            constants.PROTON_MASS * 3.0,
            constants.PROTON_MASS * 4.0,
        ])

        v = np.empty((4, self.no_roots))
        v[0, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (
                neoclassics_variables.kt[0, :] / (mass[0] * constants.SPEED_LIGHT**2)
                + 1
            )
            ** (-1)
        )
        v[1, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (
                neoclassics_variables.kt[1, :] / (mass[1] * constants.SPEED_LIGHT**2)
                + 1
            )
            ** (-1)
        )
        v[2, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (
                neoclassics_variables.kt[2, :] / (mass[2] * constants.SPEED_LIGHT**2)
                + 1
            )
            ** (-1)
        )
        v[3, :] = constants.SPEED_LIGHT * np.sqrt(
            1.0
            - (
                neoclassics_variables.kt[3, :] / (mass[3] * constants.SPEED_LIGHT**2)
                + 1
            )
            ** (-1)
        )

        return (
            np.pi
            / 4.0
            * neoclassics_variables.vd**2
            * physics_variables.rmajor
            / neoclassics_variables.iota
            / v
        )

    def neoclassics_calc_d11_mono(self, eps_eff):
        """Calculates the monoenergetic radial transport coefficients
        using epsilon effective
        """
        return (
            4.0
            / (9.0 * np.pi)
            * (2.0 * eps_eff) ** (3.0 / 2.0)
            * neoclassics_variables.vd**2
            / neoclassics_variables.nu
        )

    def calc_integrated_radial_transport_coeffs(self, index: int):
        """Calculates the integrated radial transport coefficients (index `index`)
        It uses Gauss laguerre integration
        https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature
        """
        return np.sum(
            2.0
            / np.sqrt(np.pi)
            * neoclassics_variables.d11_mono
            * neoclassics_variables.roots ** (index - 0.5)
            * neoclassics_variables.weights,
            axis=1,
        )

    def neoclassics_calc_gamma_flux(
        self, densities, temperatures, dr_densities, dr_temperatures
    ):
        """Calculates the Energy flux by neoclassical particle transport"""

        z = np.array([-1.0, 1.0, 1.0, 2.0])

        return (
            -densities
            * neoclassics_variables.d111
            * (
                (dr_densities / densities - z * neoclassics_variables.er / temperatures)
                + (neoclassics_variables.d112 / neoclassics_variables.d111 - 3.0 / 2.0)
                * dr_temperatures
                / temperatures
            )
        )

    def neoclassics_calc_q_flux(self):
        """Calculates the Energy flux by neoclassicsal energy transport"""

        z = np.array([-1.0, 1.0, 1.0, 2.0])

        return (
            -neoclassics_variables.densities
            * neoclassics_variables.temperatures
            * neoclassics_variables.d112
            * (
                (
                    neoclassics_variables.dr_densities / neoclassics_variables.densities
                    - z * neoclassics_variables.er / neoclassics_variables.temperatures
                )
                + (neoclassics_variables.d113 / neoclassics_variables.d112 - 3.0 / 2.0)
                * neoclassics_variables.dr_temperatures
                / neoclassics_variables.temperatures
            )
        )


def stinit():
    """Routine to initialise the variables relevant to stellarators
    author: P J Knight, CCFE, Culham Science Centre
    author: F Warmer, IPP Greifswald

    This routine initialises the variables relevant to stellarators.
    Many of these may override the values set in routine
    """
    if stellarator_variables.istell == 0:
        return

    numerics.boundu[0] = 40.0  # allow higher aspect ratio

    # These lines switch off tokamak specifics (solenoid, pf coils, pulses etc.).
    # Are they still up to date? (26/07/22 JL)

    # Build quantities

    build_variables.dr_cs = 0.0
    build_variables.iohcl = 0
    pfcoil_variables.f_z_cs_tf_internal = 0.0
    build_variables.dr_cs_tf_gap = 0.0
    build_variables.f_dr_tf_outboard_inboard = 1.0

    #  Physics quantities

    physics_variables.beta_norm_max = 0.0
    physics_variables.kappa95 = 1.0
    physics_variables.triang = 0.0
    physics_variables.q95 = 1.03

    #  Turn off current drive

    current_drive_variables.i_hcd_calculations = 0

    #  Times for different phases

    times_variables.t_plant_pulse_coil_precharge = 0.0
    times_variables.t_plant_pulse_plasma_current_ramp_up = 0.0
    times_variables.t_plant_pulse_burn = 3.15576e7  # one year
    times_variables.t_plant_pulse_plasma_current_ramp_down = 0.0
    times_variables.t_plant_pulse_plasma_present = (
        times_variables.t_plant_pulse_plasma_current_ramp_up
        + times_variables.t_plant_pulse_fusion_ramp
        + times_variables.t_plant_pulse_burn
        + times_variables.t_plant_pulse_plasma_current_ramp_down
    )
    times_variables.t_plant_pulse_no_burn = (
        times_variables.t_plant_pulse_coil_precharge
        + times_variables.t_plant_pulse_plasma_current_ramp_up
        + times_variables.t_plant_pulse_plasma_current_ramp_down
        + times_variables.t_plant_pulse_dwell
        + times_variables.t_plant_pulse_fusion_ramp
    )
    times_variables.t_plant_pulse_total = (
        times_variables.t_plant_pulse_coil_precharge
        + times_variables.t_plant_pulse_plasma_current_ramp_up
        + times_variables.t_plant_pulse_fusion_ramp
        + times_variables.t_plant_pulse_burn
        + times_variables.t_plant_pulse_plasma_current_ramp_down
        + times_variables.t_plant_pulse_dwell
    )
