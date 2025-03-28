import logging
from copy import copy
from pathlib import Path

import numpy as np

import process.fusion_reactions as reactions
import process.physics_functions as physics_funcs
import process.superconductors as superconductors
from process import (
    process_output as po,
)
from process.coolprop_interface import FluidProperties
from process.fortran import (
    build_variables,
    constants,
    constraint_variables,
    cost_variables,
    current_drive_variables,
    divertor_variables,
    error_handling,
    fwbs_variables,
    global_variables,
    heat_transport_variables,
    impurity_radiation_module,
    neoclassics_module,
    numerics,
    physics_module,
    physics_variables,
    rebco_variables,
    sctfcoil_module,
    stellarator_configuration,
    stellarator_variables,
    structure_variables,
    tfcoil_variables,
)
from process.fortran import (
    stellarator_module as st,
)
from process.physics import rether
from process.stellarator_config import load_stellarator_config
from process.utilities.f2py_string_patch import f2py_compatible_to_string

logger = logging.getLogger(__name__)
# Logging handler for console output
s_handler = logging.StreamHandler()
s_handler.setLevel(logging.ERROR)
logger.addHandler(s_handler)

# NOTE: a different value of electron_charge was used in the original implementation
# making the post-Python results slightly different. As a result, there is a
# relative tolerance on the neoclassics tests of 1e-3
KEV = 1e3 * constants.electron_charge  # Kiloelectron-volt (keV)


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

        self.outfile: int = constants.nout
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

            # As stopt changes dene, te and bt, stphys needs two calls
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
            self.power.power2(output=True)

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
        self.power.power1()
        self.buildings.run(output=False)
        self.vacuum.run(output=False)
        self.power.acpow(output=False)
        self.power.power2(output=False)
        # TODO: should availability.run be called
        # rather than availability.avail?
        self.availability.avail(output=False)
        self.costs.run()

        if any(numerics.icc == 91):
            # This call is comparably time consuming..
            # If the respective constraint equation is not called, do not set the values
            (
                stellarator_variables.powerht_constraint,
                stellarator_variables.powerscaling_constraint,
            ) = self.power_at_ignition_point(
                stellarator_variables.max_gyrotron_frequency,
                stellarator_variables.te0_ecrh_achievable,
            )

        st.first_call = False

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
            Path(
                f"{f2py_compatible_to_string(global_variables.output_prefix)}stella_conf.json"
            ),
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
        st.f_r = (
            physics_variables.rmajor
            / stellarator_configuration.stella_config_rmajor_ref
        )  # Size scaling factor with respect to the reference calculation
        st.f_a = (
            physics_variables.rminor
            / stellarator_configuration.stella_config_rminor_ref
        )  # Size scaling factor with respect to the reference calculation

        st.f_aspect = (
            physics_variables.aspect
            / stellarator_configuration.stella_config_aspect_ref
        )
        st.f_n = tfcoil_variables.n_tf_coils / (
            stellarator_configuration.stella_config_coilspermodule
            * stellarator_configuration.stella_config_symmetry
        )  # Coil number factor
        st.f_b = (
            physics_variables.bt / stellarator_configuration.stella_config_bt_ref
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
            st.f_r * st.f_a**2 * stellarator_configuration.stella_config_vol_plasma
        )

        # Plasma surface scaled from effective parameter:
        physics_variables.a_plasma_surface = (
            st.f_r * st.f_a * stellarator_configuration.stella_config_plasma_surface
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

        physics_variables.dnelimt = self.stdlim(
            physics_variables.bt,
            physics_variables.p_plasma_loss_mw,
            physics_variables.rmajor,
            physics_variables.rminor,
        )

        # Calculates the ECRH parameters

        ne0_max_ECRH, bt_ecrh = self.stdlim_ecrh(
            stellarator_variables.max_gyrotron_frequency, physics_variables.bt
        )

        ne0_max_ECRH = min(physics_variables.ne0, ne0_max_ECRH)
        bt_ecrh = min(physics_variables.bt, bt_ecrh)

        if output:
            self.stopt_output(
                stellarator_variables.max_gyrotron_frequency,
                physics_variables.bt,
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
        build_variables.available_radial_space = st.f_r * (
            stellarator_configuration.stella_config_derivative_min_lcfs_coils_dist
            * stellarator_configuration.stella_config_rminor_ref
            * (1 / st.f_aspect - 1)
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

        build_variables.hmax = 0.5e0 * (
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
                - fwbs_variables.f_a_fw_hcd
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
        m_struc = 1.3483e0 * (1000.0e0 * tfcoil_variables.estotftgj) ** 0.7821e0
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
            stellarator_configuration.stella_config_coilsurface * st.f_r**2
            - tfcoil_variables.tftort
            * stellarator_configuration.stella_config_coillength
            * st.f_r
            * st.f_n
        )

        # This 0.18 m is an effective thickness which is scaled with empirial 1.5 law. 5.6 T is reference point of Helias
        # The thickness 0.18m was obtained as a measured value from Schauer, F. and Bykov, V. design of Helias 5-B. (Nucl Fus. 2013)
        structure_variables.aintmass = (
            0.18e0 * st.f_b**2 * intercoil_surface * fwbs_variables.denstl
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
        p_div = physics_variables.pdivt
        alpha = divertor_variables.anginc
        xi_p = divertor_variables.xpertin
        T_scrape = divertor_variables.tdiv

        #  Scrape-off temperature in Joules

        e = T_scrape * constants.electron_charge

        #  Sound speed of particles (m/s)

        c_s = np.sqrt(e / (physics_variables.m_fuel_amu * constants.umass))

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

        divertor_variables.hldiv = q_div
        divertor_variables.divsur = darea

        fwbs_variables.f_ster_div_single = darea / build_variables.a_fw_total

        if output:
            po.oheadr(self.outfile, "Divertor")

            po.ovarre(
                self.outfile,
                "Power to divertor (MW)",
                "(pdivt.)",
                physics_variables.pdivt,
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
                "(hldiv)",
                divertor_variables.hldiv,
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
            0.25 * tfcoil_variables.len_tf_coil * tfcoil_variables.tfareain
            + 0.75
            * tfcoil_variables.len_tf_coil
            * tfcoil_variables.a_tf_leg_outboard
            * tfcoil_variables.n_tf_coils
        )

        fwbs_variables.ptfnucpm3 = fwbs_variables.ptfnuc / tf_volume

        # heating of the shield
        self.hcpb.nuclear_heating_shield()

        # Energy multiplication factor
        fwbs_variables.emult = 1.269

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
            build_variables.blarea = (
                physics_variables.a_plasma_surface
                * r1
                / physics_variables.rminor
                * (1.0e0 - fwbs_variables.fhole)
            )
        else:
            build_variables.blarea = (
                physics_variables.a_plasma_surface
                * r1
                / physics_variables.rminor
                * (
                    1.0e0
                    - fwbs_variables.fhole
                    - fwbs_variables.f_ster_div_single
                    - fwbs_variables.f_a_fw_hcd
                )
            )

        build_variables.blareaib = 0.5e0 * build_variables.blarea
        build_variables.blareaob = 0.5e0 * build_variables.blarea

        fwbs_variables.vol_blkt_inboard = (
            build_variables.blareaib * build_variables.dr_blkt_inboard
        )
        fwbs_variables.vol_blkt_outboard = (
            build_variables.blareaob * build_variables.dr_blkt_outboard
        )
        fwbs_variables.vol_blkt_total = (
            fwbs_variables.vol_blkt_inboard + fwbs_variables.vol_blkt_outboard
        )

        #  Shield volume
        #  Uses fvolsi, fwbs_variables.fvolso as area coverage factors

        r1 = r1 + 0.5e0 * (
            build_variables.dr_blkt_inboard + build_variables.dr_blkt_outboard
        )
        build_variables.sharea = (
            physics_variables.a_plasma_surface * r1 / physics_variables.rminor
        )
        build_variables.shareaib = (
            0.5e0 * build_variables.sharea * fwbs_variables.fvolsi
        )
        build_variables.shareaob = (
            0.5e0 * build_variables.sharea * fwbs_variables.fvolso
        )

        volshldi = build_variables.shareaib * build_variables.dr_shld_inboard
        volshldo = build_variables.shareaob * build_variables.dr_shld_outboard
        fwbs_variables.volshld = volshldi + volshldo

        #  Neutron power lost through holes in first wall (eventually absorbed by
        #  shield)

        fwbs_variables.pnucloss = (
            physics_variables.neutron_power_total * fwbs_variables.fhole
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
                    physics_variables.neutron_power_total
                    * fwbs_variables.f_ster_div_single
                )
                fwbs_variables.p_fw_hcd_nuclear_heat_mw = (
                    physics_variables.neutron_power_total * fwbs_variables.f_a_fw_hcd
                )
                fwbs_variables.p_fw_nuclear_heat_total_mw = (
                    physics_variables.neutron_power_total
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
                    physics_variables.p_plasma_rad_mw * fwbs_variables.f_a_fw_hcd
                )
                fwbs_variables.p_fw_rad_total_mw = (
                    physics_variables.p_plasma_rad_mw
                    - fwbs_variables.p_div_rad_total_mw
                    - fwbs_variables.pradloss
                    - fwbs_variables.p_fw_hcd_rad_total_mw
                )

                heat_transport_variables.htpmw_fw = heat_transport_variables.fpumpfw * (
                    fwbs_variables.p_fw_nuclear_heat_total_mw
                    + fwbs_variables.p_fw_rad_total_mw
                    + current_drive_variables.porbitlossmw
                )
                heat_transport_variables.htpmw_blkt = (
                    heat_transport_variables.fpumpblkt
                    * fwbs_variables.p_blkt_nuclear_heat_total_mw
                )
                heat_transport_variables.htpmw_shld = (
                    heat_transport_variables.fpumpshld * fwbs_variables.pnucshld
                )
                heat_transport_variables.htpmw_div = (
                    heat_transport_variables.fpumpdiv
                    * (
                        physics_variables.pdivt
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
                    physics_variables.neutron_power_total
                    - fwbs_variables.pnucloss
                    - fwbs_variables.pnuc_cp
                ) * fwbs_variables.emult

                fwbs_variables.emultmw = pneut2 - (
                    physics_variables.neutron_power_total
                    - fwbs_variables.pnucloss
                    - fwbs_variables.pnuc_cp
                )

                #  Nuclear heating in the blanket

                decaybl = 0.075e0 / (
                    1.0e0
                    - fwbs_variables.vfblkt
                    - fwbs_variables.fblli2o
                    - fwbs_variables.fblbe
                )

                fwbs_variables.p_blkt_nuclear_heat_total_mw = pneut2 * (
                    1.0e0 - np.exp(-build_variables.dr_blkt_outboard / decaybl)
                )

                #  Nuclear heating in the shield
                fwbs_variables.pnucshld = (
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
                    fwbs_variables.ptfnuc,
                ) = self.sctfcoil_nuclear_heating_iter90()

            else:  # heat_transport_variables.ipowerflow == 1
                #  Neutron power incident on divertor (MW)

                fwbs_variables.p_div_nuclear_heat_total_mw = (
                    physics_variables.neutron_power_total
                    * fwbs_variables.f_ster_div_single
                )

                #  Neutron power incident on HCD apparatus (MW)

                fwbs_variables.p_fw_hcd_nuclear_heat_mw = (
                    physics_variables.neutron_power_total * fwbs_variables.f_a_fw_hcd
                )

                #  Neutron power deposited in first wall, blanket and shield (MW)

                pnucfwbs = (
                    physics_variables.neutron_power_total
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
                    physics_variables.p_plasma_rad_mw * fwbs_variables.f_a_fw_hcd
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

                #  Simple blanket model (fwbs_variables.i_coolant_pumping = 0 or 1) is assumed for stellarators

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

                if fwbs_variables.i_coolant_pumping == 0:
                    #    Use input
                    pass
                elif fwbs_variables.i_coolant_pumping == 1:
                    heat_transport_variables.htpmw_fw = (
                        heat_transport_variables.fpumpfw
                        * (
                            p_fw_inboard_nuclear_heat_mw
                            + p_fw_outboard_nuclear_heat_mw
                            + psurffwi
                            + psurffwo
                            + current_drive_variables.porbitlossmw
                        )
                    )
                    heat_transport_variables.htpmw_blkt = (
                        heat_transport_variables.fpumpblkt
                        * (
                            pnucbzi * fwbs_variables.emult
                            + pnucbzo * fwbs_variables.emult
                        )
                    )
                else:
                    error_handling.report_error(215)

                fwbs_variables.emultmw = (
                    heat_transport_variables.fpumpblkt
                    * (pnucbzi * fwbs_variables.emult + pnucbzo)
                    * (fwbs_variables.emult - 1.0e0)
                )

                #  Total nuclear heating of first wall (MW)

                fwbs_variables.p_fw_nuclear_heat_total_mw = (
                    p_fw_inboard_nuclear_heat_mw + p_fw_outboard_nuclear_heat_mw
                )

                #  Total nuclear heating of blanket (MW)

                fwbs_variables.p_blkt_nuclear_heat_total_mw = (
                    pnucbzi + pnucbzo
                ) * fwbs_variables.emult

                fwbs_variables.emultmw = fwbs_variables.emultmw + (
                    pnucbzi + pnucbzo
                ) * (fwbs_variables.emult - 1.0e0)

                #  Calculation of shield and divertor powers
                #  Shield and divertor powers and pumping powers are calculated using the same
                #  simplified method as the first wall and breeder zone when fwbs_variables.i_coolant_pumping = 1.
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

                fwbs_variables.pnucshld = pnucshldi + pnucshldo

                #  Calculate coolant pumping powers from input fraction.
                #  The pumping power is assumed to be a fraction, fpump, of the incident
                #  thermal power to each component so that,
                #     htpmw_i = fpump_i*C
                #  where C is the non-pumping thermal power deposited in the coolant

                if fwbs_variables.i_coolant_pumping == 1:
                    #  Shield pumping power (MW)
                    heat_transport_variables.htpmw_shld = (
                        heat_transport_variables.fpumpshld * (pnucshldi + pnucshldo)
                    )

                    #  Divertor pumping power (MW)
                    heat_transport_variables.htpmw_div = (
                        heat_transport_variables.fpumpdiv
                        * (
                            physics_variables.pdivt
                            + fwbs_variables.p_div_nuclear_heat_total_mw
                            + fwbs_variables.p_div_rad_total_mw
                        )
                    )

                #  Remaining neutron power to coils and else:where. This is assumed
                #  (for superconducting coils at least) to be absorbed by the
                #  coils, and so contributes to the cryogenic load

                if tfcoil_variables.i_tf_sup == 1:
                    fwbs_variables.ptfnuc = pnucsi + pnucso - pnucshldi - pnucshldo
                else:  # resistive coils
                    fwbs_variables.ptfnuc = 0.0e0

        #  heat_transport_variables.ipowerflow

        #  fwbs_variables.blktmodel

        #  Divertor mass
        #  N.B. divertor_variables.divsur is calculated in stdiv after this point, so will
        #  be zero on first lap, hence the initial approximation

        if self.first_call_stfwbs:
            divertor_variables.divsur = 50.0e0
            self.first_call_stfwbs = False

        divertor_variables.divmas = (
            divertor_variables.divsur
            * divertor_variables.divdens
            * (1.0e0 - divertor_variables.divclfr)
            * divertor_variables.divplt
        )

        #  Start adding components of the coolant mass:
        #  Divertor coolant volume (m3)

        coolvol = (
            divertor_variables.divsur
            * divertor_variables.divclfr
            * divertor_variables.divplt
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
                * fwbs_variables.denstl
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
            fwbs_variables.m_blkt_steel_total = fwbs_variables.denstl * (
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

            fwbs_variables.vfblkt = (
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
            fwbs_variables.vfblkt = (
                fwbs_variables.vfblkt
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
            coolvol = coolvol + fwbs_variables.vol_blkt_total * fwbs_variables.vfblkt

        # Shield mass
        fwbs_variables.whtshld = (
            fwbs_variables.volshld
            * fwbs_variables.denstl
            * (1.0e0 - fwbs_variables.vfshld)
        )

        coolvol = coolvol + fwbs_variables.volshld * fwbs_variables.vfshld

        #  Penetration shield (set = internal shield)

        fwbs_variables.wpenshld = fwbs_variables.whtshld

        if heat_transport_variables.ipowerflow == 0:
            #  First wall mass
            #  (first wall area is calculated else:where)

            fwbs_variables.m_fw_total = (
                build_variables.a_fw_total
                * (build_variables.dr_fw_inboard + build_variables.dr_fw_outboard)
                / 2.0e0
                * fwbs_variables.denstl
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
            fwbs_variables.m_fw_total = fwbs_variables.denstl * (
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

        fwbs_variables.m_vv = fwbs_variables.vol_vv * fwbs_variables.denstl

        #  Sum of internal vacuum vessel and external cryostat masses

        fwbs_variables.dewmkg = (
            fwbs_variables.vol_vv + fwbs_variables.vol_cryostat
        ) * fwbs_variables.denstl

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
                    "(pnucshld)",
                    fwbs_variables.pnucshld,
                )
                po.ovarre(
                    self.outfile,
                    "Coil nuclear heating (MW)",
                    "(ptfnuc)",
                    fwbs_variables.ptfnuc,
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
                    "(pnucshld)",
                    fwbs_variables.pnucshld,
                )
                po.ovarre(
                    self.outfile,
                    "Energy multiplication in blanket",
                    "(emult)",
                    fwbs_variables.emult,
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
                    "(ptfnuc)",
                    fwbs_variables.ptfnuc,
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
            #             po.write(self.outfile,601) vol_blkt_inboard, vol_blkt_outboard, vol_blkt_total,                m_blkt_total, vfblkt, fbllipb, wtbllipb, fblli, m_blkt_lithium,                fblss, m_blkt_steel_total, fblvd, m_blkt_vanadium, volshldi, volshldo,                volshld, whtshld, vfshld, fwbs_variables.wpenshld
            #         else:  #  (also if ipowerflow=0)
            #             po.write(self.outfile,600) vol_blkt_inboard, vol_blkt_outboard, vol_blkt_total,                m_blkt_total, vfblkt, fblbe, m_blkt_beryllium, fblli2o, m_blkt_li2o,                fblss, m_blkt_steel_total, fblvd, m_blkt_vanadium, volshldi, volshldo,                volshld, whtshld, vfshld, fwbs_variables.wpenshld

            #     else:
            #         po.write(self.outfile,602) vol_blkt_inboard, vol_blkt_outboard, vol_blkt_total, m_blkt_total, vfblkt,             (fwbs_variables.vol_blkt_inboard/fwbs_variables.vol_blkt_total * build_variables.blbuith/build_variables.dr_blkt_inboard +             fwbs_variables.vol_blkt_outboard/fwbs_variables.vol_blkt_total * build_variables.blbuoth/build_variables.dr_blkt_outboard) * fblbe, m_blkt_beryllium,             (fwbs_variables.vol_blkt_inboard/fwbs_variables.vol_blkt_total * build_variables.blbuith/build_variables.dr_blkt_inboard +             fwbs_variables.vol_blkt_outboard/fwbs_variables.vol_blkt_total * build_variables.blbuoth/build_variables.dr_blkt_outboard) * fblbreed, whtblbreed,             fwbs_variables.vol_blkt_inboard/fwbs_variables.vol_blkt_total/build_variables.dr_blkt_inboard * (build_variables.blbuith * fwbs_variables.fblss             + build_variables.blbmith * (1.0e0-fwbs_variables.fblhebmi) + build_variables.blbpith * (1.0e0-fwbs_variables.fblhebpi)) +             fwbs_variables.vol_blkt_outboard/fwbs_variables.vol_blkt_total/build_variables.dr_blkt_outboard * (build_variables.blbuoth * fwbs_variables.fblss             + build_variables.blbmoth * (1.0e0-fwbs_variables.fblhebmo) + build_variables.blbpoth * (1.0e0-fwbs_variables.fblhebpo)),             m_blkt_steel_total,             volshldi, volshldo, volshld, whtshld, vfshld, fwbs_variables.wpenshld

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
                "(divsur)",
                divertor_variables.divsur,
            )
            po.ovarre(
                self.outfile,
                "Divertor mass (kg)",
                "(divmas)",
                divertor_variables.divmas,
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
        ptfnuc : output real : TF coil nuclear heating (MW)
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
            ptfnuc = 0.0

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

            wpthk = tfcoil_variables.dr_tf_wp + 2.0 * tfcoil_variables.tinstf

            # Nuclear heating rate in inboard TF coil (MW/m**3)

            coilhtmx = (
                fact[0]
                * physics_variables.pflux_fw_neutron_mw
                * coef[0, ishmat]
                * np.exp(-decay[5, ishmat] * (dshieq + tfcoil_variables.casthi))
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
                * np.exp(-decay[5, ishmat] * (dshoeq + tfcoil_variables.casthi))
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
                * (1.0 - np.exp(-decay[1, ishmat] * tfcoil_variables.casthi))
                / decay[1, ishmat]
            )
            pheco = (
                fact[1]
                * physics_variables.pflux_fw_neutron_mw
                * coef[1, ishmat]
                * np.exp(-decay[6, ishmat] * dshoeq)
                * tfcoil_variables.tfsao
                * (1.0 - np.exp(-decay[1, ishmat] * tfcoil_variables.casthi))
                / decay[1, ishmat]
            )
            ptfi = ptfiwp + pheci
            ptfo = ptfowp + pheco

            ptfnuc = ptfi + ptfo

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
                * np.exp(-decay[2, ishmat] * (dshieq + tfcoil_variables.casthi))
            )

            # Maximum neutron fluence in superconductor (n/m**2)

            nflutf = (
                fpsdt
                * fact[3]
                * physics_variables.pflux_fw_neutron_mw
                * coef[3, ishmat]
                * np.exp(-decay[3, ishmat] * (dshieq + tfcoil_variables.casthi))
            )

            # Atomic displacement in copper stabilizer

            dpacop = (
                fpsdt
                * fact[4]
                * physics_variables.pflux_fw_neutron_mw
                * coef[4, ishmat]
                * np.exp(-decay[4, ishmat] * (dshieq + tfcoil_variables.casthi))
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
            ptfnuc,
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
        r_coil_major = stellarator_configuration.stella_config_coil_rmajor * st.f_r
        r_coil_minor = stellarator_configuration.stella_config_coil_rminor * st.f_r

        ########################################################################################
        # Winding Pack Geometry: for one conductor
        #
        # This one conductor will just be multiplied later to fit the winding pack size.
        #
        # [m] Dimension of square cable space inside insulation
        #     and case of the conduit of each turn
        t_cable = tfcoil_variables.t_turn_tf - 2.0e0 * (
            tfcoil_variables.thwcndut + tfcoil_variables.thicndut
        )  # t_cable = t_w
        if t_cable < 0:
            print(
                "t_cable is negative. Check t_turn, tfcoil_variables.thwcndut and thicndut."
            )
        # [m^2] Cross-sectional area of cable space per turn
        tfcoil_variables.acstf = (
            0.9e0 * t_cable**2
        )  # 0.9 to include some rounded corners. (tfcoil_variables.acstf = pi (t_cable/2)**2 = pi/4 *t_cable**2 for perfect round conductor). This factor depends on how round the corners are.
        # [m^2] Cross-sectional area of conduit case per turn
        tfcoil_variables.acndttf = (
            t_cable + 2.0e0 * tfcoil_variables.thwcndut
        ) ** 2 - tfcoil_variables.acstf
        #######################################################################################

        #######################################################################################
        # Winding Pack total size:
        #
        # Total coil current (MA)
        coilcurrent = (
            st.f_b * stellarator_configuration.stella_config_i0 * st.f_r / st.f_n
        )
        st.f_i = coilcurrent / stellarator_configuration.stella_config_i0

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
                tfcoil_variables.fcutfsu,
                tfcoil_variables.fhts,
                tfcoil_variables.t_crit_nbti,
                tfcoil_variables.tcritsc,
                tfcoil_variables.vftf,
                tfcoil_variables.jwptf,
            )  # Get here a temperature margin of 1.5K.

        # The operation current density weighted with the global iop/icrit fraction
        lhs[:] = constraint_variables.fiooic * jcrit_vector

        # Conduct fraction of conduit * Superconductor fraction in conductor
        f_scu = (
            (tfcoil_variables.acstf * (1.0e0 - tfcoil_variables.vftf))
            / (tfcoil_variables.t_turn_tf**2)
            * (1.0e0 - tfcoil_variables.fcutfsu)
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
        wp_width_r_min = max(tfcoil_variables.t_turn_tf**2, wp_width_r_min)

        # Recalculate tfcoil_variables.bmaxtf at the found awp_min:
        tfcoil_variables.bmaxtf = self.bmax_from_awp(
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

        tfcoil_variables.wwp1 = awp_tor  # [m] toroidal thickness of winding pack
        tfcoil_variables.wwp2 = (
            awp_tor  # [m] toroidal thickness of winding pack (region in front)
        )
        tfcoil_variables.dr_tf_wp = awp_rad  # [m] radial thickness of winding pack

        #  [m^2] winding-pack cross sectional area including insulation (not global)
        awpc = (tfcoil_variables.dr_tf_wp + 2.0e0 * tfcoil_variables.tinstf) * (
            tfcoil_variables.wwp1 + 2.0e0 * tfcoil_variables.tinstf
        )

        awptf = awp_tor * awp_rad  # [m^2] winding-pack cross sectional area
        tfcoil_variables.jwptf = (
            coilcurrent * 1.0e6 / awptf
        )  # [A/m^2] winding pack current density
        tfcoil_variables.n_tf_turn = (
            awptf / (tfcoil_variables.t_turn_tf**2)
        )  # estimated number of turns for a given turn size (not global). Take at least 1.
        tfcoil_variables.cpttf = (
            coilcurrent * 1.0e6 / tfcoil_variables.n_tf_turn
        )  # [A] current per turn - estimation
        # [m^2] Total conductor cross-sectional area, taking account of void area
        tfcoil_variables.acond = (
            tfcoil_variables.acstf
            * tfcoil_variables.n_tf_turn
            * (1.0e0 - tfcoil_variables.vftf)
        )
        # [m^2] Void area in cable, for He
        tfcoil_variables.avwp = (
            tfcoil_variables.acstf * tfcoil_variables.n_tf_turn * tfcoil_variables.vftf
        )
        # [m^2] Insulation area (not including ground-wall)
        tfcoil_variables.aiwp = tfcoil_variables.n_tf_turn * (
            tfcoil_variables.t_turn_tf**2
            - tfcoil_variables.acndttf
            - tfcoil_variables.acstf
        )
        # [m^2] Structure area for cable
        tfcoil_variables.aswp = tfcoil_variables.n_tf_turn * tfcoil_variables.acndttf
        # End of winding pack calculations
        #######################################################################################

        #######################################################################################
        #  Casing calculations
        #
        # Coil case thickness (m). Here assumed to be constant
        # until something better comes up.
        # case_thickness_constant = tfcoil_variables.thkcas #0.2e0 # #? Leave this constant for now... Check this## Should be scaled with forces I think.
        #  For now assumed to be constant in a bolted plate model.
        #
        tfcoil_variables.casthi = (
            tfcoil_variables.thkcas
        )  # [m] coil case thickness outboard distance (radial)
        # thkcas = case_thickness_constant/2.0e0 # [m] coil case thickness inboard distance  (radial).
        tfcoil_variables.casths = (
            tfcoil_variables.thkcas
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
            * st.f_r
            / st.f_n
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
            * st.f_r
            / st.f_n
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
        tfcoil_variables.tftort = (
            tfcoil_variables.wwp1
            + 2.0e0 * tfcoil_variables.casths
            + 2.0e0 * tfcoil_variables.tinstf
        )  # [m] Thickness of inboard leg in toroidal direction

        build_variables.dr_tf_inboard = (
            tfcoil_variables.thkcas
            + tfcoil_variables.dr_tf_wp
            + tfcoil_variables.casthi
            + 2.0e0 * tfcoil_variables.tinstf
        )  # [m] Thickness of inboard leg in radial direction
        build_variables.dr_tf_outboard = (
            tfcoil_variables.thkcas
            + tfcoil_variables.dr_tf_wp
            + tfcoil_variables.casthi
            + 2.0e0 * tfcoil_variables.tinstf
        )  # [m] Thickness of outboard leg in radial direction (same as inboard)
        tfcoil_variables.a_tf_leg_outboard = (
            build_variables.dr_tf_inboard * tfcoil_variables.tftort
        )  # [m^2] overall coil cross-sectional area (assuming inboard and
        #       outboard leg are the same)
        tfcoil_variables.acasetf = (
            build_variables.dr_tf_inboard * tfcoil_variables.tftort
        ) - awpc  # [m^2] Cross-sectional area of surrounding case

        tfcoil_variables.tfocrn = (
            0.5e0 * tfcoil_variables.tftort
        )  # [m] Half-width of side of coil nearest torus centreline
        tfcoil_variables.tficrn = (
            0.5e0 * tfcoil_variables.tftort
        )  # [m] Half-width of side of coil nearest plasma

        # [m^2] Total surface area of coil side facing plasma: inboard region
        tfcoil_variables.tfsai = (
            tfcoil_variables.n_tf_coils
            * tfcoil_variables.tftort
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
        coilcoilgap = tfcoil_variables.toroidalgap - tfcoil_variables.tftort

        #  Variables for ALL coils.
        tfcoil_variables.tfareain = (
            tfcoil_variables.n_tf_coils * tfcoil_variables.a_tf_leg_outboard
        )  # [m^2] Total area of all coil legs (midplane)
        tfcoil_variables.c_tf_total = (
            tfcoil_variables.n_tf_coils * coilcurrent * 1.0e6
        )  # [A] Total current in ALL coils
        tfcoil_variables.oacdcp = (
            tfcoil_variables.c_tf_total / tfcoil_variables.tfareain
        )  # [A / m^2] overall current density
        tfcoil_variables.rbmax = (
            r_coil_major - r_coil_minor + awp_rad
        )  # [m] radius of peak field occurrence, average
        # jlion: not sure what this will be used for. Not very
        # useful for stellarators

        # This uses the reference value for the inductance and scales it with a^2/R (toroid inductance scaling)
        inductance = (
            stellarator_configuration.stella_config_inductance
            / st.f_r
            * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor) ** 2
            * st.f_n**2
        )
        tfcoil_variables.estotftgj = (
            0.5e0
            * (
                stellarator_configuration.stella_config_inductance
                / st.f_r
                * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
                ** 2
                * st.f_n**2
            )
            * (tfcoil_variables.c_tf_total / tfcoil_variables.n_tf_coils) ** 2
            * 1.0e-9
        )  # [GJ] Total magnetic energy

        #  Coil dimensions
        build_variables.hmax = (
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

        tfborev = 2.0e0 * build_variables.hmax  # [m] estimated vertical coil dr_bore

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
            * st.f_r
            * 1.0
            / (1.0 - tfcoil_variables.dr_tf_wp / (2.0 * r_coil_minor))
        )

        # End of general coil geometry values
        #######################################################################################

        #######################################################################################
        #  Masses of conductor constituents
        #
        # [kg] Mass of case
        #  (no need for correction factors as is the case for tokamaks)
        # This is only correct if the winding pack is 'thin' (len_tf_coil>>sqrt(tfcoil_variables.acasetf)).
        tfcoil_variables.whtcas = (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.acasetf
            * tfcoil_variables.dcase
        )
        # Mass of ground-wall insulation [kg]
        # (assumed to be same density/material as conduit insulation)
        tfcoil_variables.whtgw = (
            tfcoil_variables.len_tf_coil * (awpc - awptf) * tfcoil_variables.dcondins
        )
        # [kg] mass of Superconductor
        tfcoil_variables.whtconsc = (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.n_tf_turn
            * tfcoil_variables.acstf
            * (1.0e0 - tfcoil_variables.vftf)
            * (1.0e0 - tfcoil_variables.fcutfsu)
            - tfcoil_variables.len_tf_coil * tfcoil_variables.awphec
        ) * tfcoil_variables.dcond[
            tfcoil_variables.i_tf_sc_mat - 1
        ]  # awphec is 0 for a stellarator. but keep this term for now.
        # [kg] mass of Copper in conductor
        tfcoil_variables.whtconcu = (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.n_tf_turn
            * tfcoil_variables.acstf
            * (1.0e0 - tfcoil_variables.vftf)
            * tfcoil_variables.fcutfsu
            - tfcoil_variables.len_tf_coil * tfcoil_variables.awphec
        ) * constants.dcopper
        # [kg] mass of Steel conduit (sheath)
        tfcoil_variables.whtconsh = (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.n_tf_turn
            * tfcoil_variables.acndttf
            * fwbs_variables.denstl
        )
        # if (i_tf_sc_mat==6)   tfcoil_variables.whtconsh = fcondsteel * awptf *tfcoil_variables.len_tf_coil* fwbs_variables.denstl
        # Conduit insulation mass [kg]
        # (tfcoil_variables.aiwp already contains tfcoil_variables.n_tf_turn)
        tfcoil_variables.whtconin = (
            tfcoil_variables.len_tf_coil
            * tfcoil_variables.aiwp
            * tfcoil_variables.dcondins
        )
        # [kg] Total conductor mass
        tfcoil_variables.whtcon = (
            tfcoil_variables.whtconsc
            + tfcoil_variables.whtconcu
            + tfcoil_variables.whtconsh
            + tfcoil_variables.whtconin
        )
        # [kg] Total coil mass
        tfcoil_variables.m_tf_coils_total = (
            tfcoil_variables.whtcas + tfcoil_variables.whtcon + tfcoil_variables.whtgw
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
                physics_variables.bt
                * tfcoil_variables.c_tf_total
                * physics_variables.rminor**2
                / (
                    (build_variables.dr_vv_inboard + build_variables.dr_vv_outboard)
                    / 2
                    * tfcoil_variables.tdmptf
                    * radvv
                )
            )
            ** (-1)
        )

        # N/m^2
        # is the vv width the correct length to multiply by to turn the
        # force density into a stress?
        sctfcoil_module.vv_stress_quench = (
            f_vv_actual
            * 1e6
            * ((build_variables.dr_vv_inboard + build_variables.dr_vv_outboard) / 2)
        )

        # the conductor fraction is meant of the cable space#
        # This is the old routine which is being replaced for now by the new one below
        #    protect(aio,  tfes,               acs,       aturn,   tdump,  fcond,  fcu,   tba,  tmax   ,ajwpro, vd)
        # call protect(cpttf,estotftgj/tfcoil_variables.n_tf_coils*1.0e9,acstf,   tfcoil_variables.t_turn_tf**2   ,tdmptf,1-vftf,fcutfsu,tftmp,tmaxpro,jwdgpro2,vd)

        vd = self.u_max_protect_v(
            tfcoil_variables.estotftgj / tfcoil_variables.n_tf_coils * 1.0e9,
            tfcoil_variables.tdmptf,
            tfcoil_variables.cpttf,
        )

        # comparison
        # the new quench protection routine, see #1047
        tfcoil_variables.jwdgpro = self.j_max_protect_am2(
            tfcoil_variables.tdmptf,
            0.0e0,
            tfcoil_variables.fcutfsu,
            1 - tfcoil_variables.vftf,
            tfcoil_variables.tftmp,
            tfcoil_variables.acstf,
            tfcoil_variables.t_turn_tf**2,
        )

        # print *, "Jmax, comparison: ", jwdgpro, "  ", jwdgpro2,"  ",jwptf/jwdgpro, "   , tfcoil_variables.tdmptf: ",tdmptf, " tfcoil_variables.fcutfsu: ",fcutfsu
        # print *, "acstf: ", tfcoil_variables.acstf
        # Also give the copper area for REBCO quench calculations:
        rebco_variables.coppera_m2 = (
            coilcurrent * 1.0e6 / (tfcoil_variables.acond * tfcoil_variables.fcutfsu)
        )
        tfcoil_variables.vtfskv = vd / 1.0e3  # Dump voltage
        #
        #######################################################################################

        # Forces scaling #
        tfcoil_variables.max_force_density = (
            stellarator_configuration.stella_config_max_force_density
            * st.f_i
            / st.f_n
            * tfcoil_variables.bmaxtf
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_wp_area
            / awptf
        )

        # Approximate, very simple maxiumum stress: (needed for limitation of icc 32)
        tfcoil_variables.sig_tf_wp = (
            tfcoil_variables.max_force_density * tfcoil_variables.dr_tf_wp * 1.0e6
        )  # in Pa

        # Units: MN/m
        max_force_density_mnm = (
            stellarator_configuration.stella_config_max_force_density_mnm
            * st.f_i
            / st.f_n
            * tfcoil_variables.bmaxtf
            / stellarator_configuration.stella_config_wp_bmax
        )
        #
        max_lateral_force_density = (
            stellarator_configuration.stella_config_max_lateral_force_density
            * st.f_i
            / st.f_n
            * tfcoil_variables.bmaxtf
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_wp_area
            / awptf
        )
        max_radial_force_density = (
            stellarator_configuration.stella_config_max_radial_force_density
            * st.f_i
            / st.f_n
            * tfcoil_variables.bmaxtf
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_wp_area
            / awptf
        )
        #
        # F = f*V = B*j*V \propto B/B0 * I/I0 * A0/A * A/A0 * len/len0
        centering_force_max_mn = (
            stellarator_configuration.stella_config_centering_force_max_mn
            * st.f_i
            / st.f_n
            * tfcoil_variables.bmaxtf
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_coillength
            / tfcoil_variables.n_tf_coils
            / tfcoil_variables.len_tf_coil
        )
        centering_force_min_mn = (
            stellarator_configuration.stella_config_centering_force_min_mn
            * st.f_i
            / st.f_n
            * tfcoil_variables.bmaxtf
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_coillength
            / tfcoil_variables.n_tf_coils
            / tfcoil_variables.len_tf_coil
        )
        centering_force_avg_mn = (
            stellarator_configuration.stella_config_centering_force_avg_mn
            * st.f_i
            / st.f_n
            * tfcoil_variables.bmaxtf
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_coillength
            / tfcoil_variables.n_tf_coils
            / tfcoil_variables.len_tf_coil
        )
        #
        ####################################

        if output:
            self.stcoil_output(
                awptf,
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
                tfcoil_variables.t_turn_tf,
                tfcoil_variables.tdmptf,
                tf_total_h_width,
                tfborev,
                tfcoil_variables.toroidalgap,
                tfcoil_variables.vdalw,
                tfcoil_variables.vtfskv,
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
        fcutfsu,
        fhts,
        t_crit_nbti,
        tcritsc,
        vftf,
        jwp,
    ):
        strain = -0.005  # for now a small value
        fhe = vftf  # this is helium fraction in the superconductor (set it to the fixed global variable here)

        fcu = fcutfsu  # fcutfsu is a global variable. Is the copper fraction
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

            j_crit_sc, bcrit, tcrit = superconductors.wstsc(
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
            error_handling.idiags[0] = i_tf_sc_mat
            error_handling.report_error(156)

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
            error_handling.fdiags[0] = np.amin(x1)
            error_handling.fdiags[1] = np.amin(x2)
            error_handling.fdiags[2] = np.amax(x1)
            error_handling.fdiags[3] = np.amax(x2)
            error_handling.report_error(111)

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
                error_handling.fdiags[0] = x
                error_handling.fdiags[1] = xmin
                error_handling.report_error(112)
                x = xmin
                break

            if x > xmax:
                error_handling.fdiags[0] = x
                error_handling.fdiags[1] = xmax
                error_handling.report_error(113)
                x = xmax
                break
        else:
            error_handling.report_error(114)

        return x

    def stopt_output(
        self, max_gyrotron_frequency, bt, bt_ecrh, ne0_max_ECRH, te0_ecrh_achievable
    ):
        po.oheadr(self.outfile, "ECRH Ignition at lower values. Information:")

        po.ovarre(
            self.outfile,
            "Maximal available gyrotron freq (input)",
            "(max_gyro_frequency)",
            max_gyrotron_frequency,
        )

        po.ovarre(self.outfile, "Operating point: bfield", "(bt)", bt)
        po.ovarre(
            self.outfile,
            "Operating point: Peak density",
            "(ne0)",
            physics_variables.ne0,
        )
        po.ovarre(
            self.outfile,
            "Operating point: Peak temperature",
            "(te0)",
            physics_variables.te0,
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
        te_old = copy(physics_variables.te)
        # Volume averaged physics_variables.te from te0_achievable
        physics_variables.te = te0_available / (1.0e0 + physics_variables.alphat)
        ne0_max, bt_ecrh_max = self.stdlim_ecrh(
            gyro_frequency_max, physics_variables.bt
        )
        # Now go to point where ECRH is still available
        # In density..
        dene_old = copy(physics_variables.dene)
        physics_variables.dene = min(
            dene_old, ne0_max / (1.0e0 + physics_variables.alphan)
        )

        # And B-field..
        bt_old = copy(physics_variables.bt)
        physics_variables.bt = min(bt_ecrh_max, physics_variables.bt)

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
        physics_variables.te = te_old
        physics_variables.dene = dene_old
        physics_variables.bt = bt_old

        self.stphys(False)
        self.stphys(False)

        return powerht_out, pscalingmw_out

    def stdlim(self, bt, powht, rmajor, rminor):
        """Routine to calculate the Sudo density limit in a stellarator
        author: P J Knight, CCFE, Culham Science Centre
        bt     : input real : Toroidal field on axis (T)
        powht  : input real : Absorbed heating power (MW)
        rmajor : input real : Plasma major radius (m)
        rminor : input real : Plasma minor radius (m)
        dlimit : output real : Maximum volume-averaged plasma density (/m3)
        This routine calculates the density limit for a stellarator.
        S.Sudo, Y.Takeiri, H.Zushi et al., Scalings of Energy Confinement
        and Density Limit in Stellarator/Heliotron Devices, Nuclear Fusion
        vol.30, 11 (1990).
        """
        arg = powht * bt / (rmajor * rminor * rminor)

        if arg <= 0.0e0:
            error_handling.fdiags[0] = arg
            error_handling.fdiags[1] = powht
            error_handling.fdiags[2] = bt
            error_handling.fdiags[3] = rmajor
            error_handling.fdiags[4] = rminor
            error_handling.report_error(108)

        #  Maximum line-averaged electron density

        dnlamx = 0.25e20 * np.sqrt(arg)

        #  Scale the result so that it applies to the volume-averaged
        #  electron density

        dlimit = dnlamx * physics_variables.dene / physics_variables.dnla

        #  Set the required value for icc=5

        physics_variables.dnelimt = dlimit

        return dlimit

    def stcoil_output(
        self,
        awptf,
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
        t_turn_tf,
        tdmptf,
        tf_total_h_width,
        tfborev,
        toroidalgap,
        vdalw,
        vtfskv,
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
            tfcoil_variables.tfareain / tfcoil_variables.n_tf_coils,
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
            "(tftort)",
            tfcoil_variables.tftort,
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
            "(jwptf)",
            tfcoil_variables.jwptf,
        )
        po.ovarre(
            self.outfile,
            "Max allowable current density as restricted by quench (A/m2)",
            "(jwdgpro)",
            tfcoil_variables.jwdgpro,
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
            "(bmaxtf)",
            tfcoil_variables.bmaxtf,
        )
        po.ovarre(
            self.outfile,
            "Total Stored energy (GJ)",
            "(estotftgj)",
            tfcoil_variables.estotftgj,
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
            "(hmax)",
            build_variables.hmax,
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
            "(whtconsc)",
            tfcoil_variables.whtconsc,
        )
        po.ovarre(
            self.outfile,
            "Copper mass per coil (kg)",
            "(whtconcu)",
            tfcoil_variables.whtconcu,
        )
        po.ovarre(
            self.outfile,
            "Steel conduit mass per coil (kg)",
            "(whtconsh)",
            tfcoil_variables.whtconsh,
        )
        po.ovarre(
            self.outfile,
            "Total conductor cable mass per coil (kg)",
            "(whtcon)",
            tfcoil_variables.whtcon,
        )
        po.ovarre(
            self.outfile,
            "Cable conductor + void area (m2)",
            "(acstf)",
            tfcoil_variables.acstf,
        )
        po.ovarre(
            self.outfile,
            "Cable space coolant fraction",
            "(vftf)",
            tfcoil_variables.vftf,
        )
        po.ovarre(
            self.outfile,
            "Conduit case thickness (m)",
            "(thwcndut)",
            tfcoil_variables.thwcndut,
        )
        po.ovarre(
            self.outfile,
            "Cable insulation thickness (m)",
            "(thicndut)",
            tfcoil_variables.thicndut,
        )

        ap = awptf
        po.osubhd(self.outfile, "Winding Pack Information :")
        po.ovarre(self.outfile, "Winding pack area", "(ap)", ap)
        po.ovarre(
            self.outfile,
            "Conductor fraction of winding pack",
            "(acond/ap)",
            tfcoil_variables.acond / ap,
        )
        po.ovarre(
            self.outfile,
            "Copper fraction of conductor",
            "(fcutfsu)",
            tfcoil_variables.fcutfsu,
        )
        po.ovarre(
            self.outfile,
            "Structure fraction of winding pack",
            "(aswp/ap)",
            tfcoil_variables.aswp / ap,
        )
        po.ovarre(
            self.outfile,
            "Insulator fraction of winding pack",
            "(aiwp/ap)",
            tfcoil_variables.aiwp / ap,
        )
        po.ovarre(
            self.outfile,
            "Helium fraction of winding pack",
            "(avwp/ap)",
            tfcoil_variables.avwp / ap,
        )
        po.ovarre(
            self.outfile,
            "Winding radial thickness (m)",
            "(dr_tf_wp)",
            tfcoil_variables.dr_tf_wp,
        )
        po.ovarre(
            self.outfile,
            "Winding toroidal thickness (m)",
            "(wwp1)",
            tfcoil_variables.wwp1,
        )
        po.ovarre(
            self.outfile,
            "Ground wall insulation thickness (m)",
            "(tinstf)",
            tfcoil_variables.tinstf,
        )
        po.ovarre(
            self.outfile,
            "Number of turns per coil",
            "(n_tf_turn)",
            tfcoil_variables.n_tf_turn,
        )
        po.ovarre(
            self.outfile,
            "Width of each turn (incl. insulation) (m)",
            "(t_turn_tf)",
            t_turn_tf,
        )
        po.ovarre(
            self.outfile, "Current per turn (A)", "(cpttf)", tfcoil_variables.cpttf
        )
        po.ovarre(self.outfile, "jop/jcrit", "(fiooic)", fiooic)
        po.ovarre(
            self.outfile,
            "Current density in conductor area (A/m2)",
            "(c_tf_total/acond)",
            1.0e-6
            * tfcoil_variables.c_tf_total
            / tfcoil_variables.n_tf_coils
            / tfcoil_variables.acond,
        )
        po.ovarre(
            self.outfile,
            "Current density in SC area (A/m2)",
            "(c_tf_total/acond/f_scu)",
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
            "(tdmptf)",
            tdmptf,
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
            "(vdalw)",
            vdalw,
        )
        po.ovarre(self.outfile, "Actual quench voltage (kV)", "(vtfskv)", vtfskv, "OP ")
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
            "(casthi)",
            tfcoil_variables.casthi,
        )
        po.ovarre(
            self.outfile,
            "Case thickness, outer side (m)",
            "(thkcas)",
            tfcoil_variables.thkcas,
        )
        po.ovarre(
            self.outfile,
            "Case toroidal thickness (m)",
            "(casths)",
            tfcoil_variables.casths,
        )
        po.ovarre(
            self.outfile,
            "Case area per coil (m2)",
            "(acasetf)",
            tfcoil_variables.acasetf,
        )
        po.ovarre(
            self.outfile,
            "External case mass per coil (kg)",
            "(whtcas)",
            tfcoil_variables.whtcas,
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
        bt  : input real : Maximal magnetic field on axis (T)
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
        if physics_variables.ipedestal == 0:
            # Parabolic profiles used, use analytical formula:
            dlimit_ecrh = ne0_max
        else:
            logger.warning(
                "It was used physics_variables.ipedestal = 1 in a stellarator routine. PROCESS will pretend it got parabolic profiles (physics_variables.ipedestal = 0)."
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
        physics_variables.btot = np.sqrt(
            physics_variables.bt**2 + physics_variables.bp**2
        )

        if (
            5 in numerics.ixc
        ):  # Check if physics_variables.beta (iteration variable 5) is an iteration variable
            error_handling.report_error(251)

        #  Set physics_variables.beta as a consequence:
        #  This replaces constraint equation 1 as it is just an equality.
        physics_variables.beta = (
            physics_variables.beta_fast_alpha
            + physics_variables.beta_beam
            + 2.0e3
            * constants.rmu0
            * constants.electron_charge
            * (
                physics_variables.dene * physics_variables.ten
                + physics_variables.nd_ions_total * physics_variables.tin
            )
            / physics_variables.btot**2
        )
        physics_module.e_plasma_beta = (
            1.5e0
            * physics_variables.beta
            * physics_variables.btot
            * physics_variables.btot
            / (2.0e0 * constants.rmu0)
            * physics_variables.vol_plasma
        )

        physics_module.rho_star = np.sqrt(
            2.0e0
            * constants.proton_mass
            * physics_variables.m_ions_total_amu
            * physics_module.e_plasma_beta
            / (3.0e0 * physics_variables.vol_plasma * physics_variables.dnla)
        ) / (
            constants.electron_charge
            * physics_variables.bt
            * physics_variables.eps
            * physics_variables.rmajor
        )

        #  Calculate poloidal field using rotation transform
        physics_variables.bp = (
            physics_variables.rminor
            * physics_variables.bt
            / physics_variables.rmajor
            * stellarator_variables.iotabar
        )

        #  Poloidal physics_variables.beta

        # beta_poloidal = physics_variables.beta * ( physics_variables.btot/physics_variables.bp )**2 # Dont need this I think.

        #  Perform auxiliary power calculations

        self.stheat(False)

        #  Calculate fusion power

        fusion_reactions = reactions.FusionReactionRate(self.plasma_profile)
        fusion_reactions.deuterium_branching(physics_variables.ti)
        fusion_reactions.calculate_fusion_rates()
        fusion_reactions.set_physics_variables()

        # D-T power density is named differently to differentiate it from the beam given component
        physics_variables.dt_power_plasma = (
            physics_module.dt_power_density_plasma * physics_variables.vol_plasma
        )
        physics_variables.dhe3_power = (
            physics_module.dhe3_power_density * physics_variables.vol_plasma
        )
        physics_variables.dd_power = (
            physics_module.dd_power_density * physics_variables.vol_plasma
        )

        #  Calculate neutral beam slowing down effects
        #  If ignited, then ignore beam fusion effects

        if (current_drive_variables.pnbeam != 0.0e0) and (
            physics_variables.ignite == 0
        ):
            (
                physics_variables.beta_beam,
                physics_variables.beam_density_out,
                physics_variables.alpha_power_beams,
            ) = reactions.beam_fusion(
                physics_variables.beamfus0,
                physics_variables.betbm0,
                physics_variables.bp,
                physics_variables.bt,
                current_drive_variables.beam_current,
                physics_variables.dene,
                physics_variables.nd_fuel_ions,
                physics_variables.dlamie,
                current_drive_variables.beam_energy,
                physics_variables.f_deuterium,
                physics_variables.f_tritium,
                current_drive_variables.f_tritium_beam,
                physics_module.sigmav_dt_average,
                physics_variables.ten,
                physics_variables.tin,
                physics_variables.vol_plasma,
                physics_variables.zeffai,
            )
            physics_variables.fusion_rate_density_total = (
                physics_variables.fusion_rate_density_plasma
                + 1.0e6
                * physics_variables.alpha_power_beams
                / (constants.dt_alpha_energy)
                / physics_variables.vol_plasma
            )
            physics_variables.alpha_rate_density_total = (
                physics_variables.alpha_rate_density_plasma
                + 1.0e6
                * physics_variables.alpha_power_beams
                / (constants.dt_alpha_energy)
                / physics_variables.vol_plasma
            )
            physics_variables.dt_power_total = (
                physics_variables.dt_power_plasma
                + 5.0e0 * physics_variables.alpha_power_beams
            )
        else:
            # If no beams present then the total alpha rates and power are the same as the plasma values
            physics_variables.fusion_rate_density_total = (
                physics_variables.fusion_rate_density_plasma
            )
            physics_variables.alpha_rate_density_total = (
                physics_variables.alpha_rate_density_plasma
            )
            physics_variables.dt_power_total = physics_variables.dt_power_plasma

        # Create some derived values and add beam contribution to fusion power
        (
            physics_variables.neutron_power_density_total,
            physics_variables.alpha_power_plasma,
            physics_variables.alpha_power_total,
            physics_variables.neutron_power_plasma,
            physics_variables.neutron_power_total,
            physics_variables.non_alpha_charged_power,
            physics_variables.alpha_power_density_total,
            physics_variables.alpha_power_electron_density,
            physics_variables.alpha_power_ions_density,
            physics_variables.charged_particle_power,
            physics_variables.fusion_power,
        ) = reactions.set_fusion_powers(
            physics_variables.f_alpha_electron,
            physics_variables.f_alpha_ion,
            physics_variables.alpha_power_beams,
            physics_variables.charged_power_density,
            physics_variables.neutron_power_density_plasma,
            physics_variables.vol_plasma,
            physics_variables.alpha_power_density_plasma,
        )

        physics_variables.beta_fast_alpha = physics_funcs.fast_alpha_beta(
            physics_variables.bp,
            physics_variables.bt,
            physics_variables.dene,
            physics_variables.nd_fuel_ions,
            physics_variables.nd_ions_total,
            physics_variables.ten,
            physics_variables.tin,
            physics_variables.alpha_power_density_total,
            physics_variables.alpha_power_density_plasma,
            physics_variables.i_beta_fast_alpha,
        )

        #  Neutron wall load

        if physics_variables.iwalld == 1:
            physics_variables.pflux_fw_neutron_mw = (
                physics_variables.ffwal
                * physics_variables.neutron_power_total
                / physics_variables.a_plasma_surface
            )
        else:
            if heat_transport_variables.ipowerflow == 0:
                physics_variables.pflux_fw_neutron_mw = (
                    (1.0e0 - fwbs_variables.fhole)
                    * physics_variables.neutron_power_total
                    / build_variables.a_fw_total
                )
            else:
                physics_variables.pflux_fw_neutron_mw = (
                    (
                        1.0e0
                        - fwbs_variables.fhole
                        - fwbs_variables.f_a_fw_hcd
                        - fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.neutron_power_total
                    / build_variables.a_fw_total
                )

        #  Calculate ion/electron equilibration power

        physics_variables.piepv = rether(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.dene,
            physics_variables.dlamie,
            physics_variables.te,
            physics_variables.ti,
            physics_variables.zeffai,
        )

        #  Calculate radiation power
        radpwr_data = physics_funcs.calculate_radiation_powers(
            self.plasma_profile,
            physics_variables.ne0,
            physics_variables.rminor,
            physics_variables.bt,
            physics_variables.aspect,
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.tbeta,
            physics_variables.te0,
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
            physics_variables.f_alpha_plasma * physics_variables.alpha_power_total
            + physics_variables.non_alpha_charged_power
            + physics_variables.p_plasma_ohmic_mw
            - physics_variables.pden_plasma_rad_mw * physics_variables.vol_plasma
        )
        powht = max(
            0.00001e0, powht
        )  # To avoid negative heating power. This line is VERY important

        if physics_variables.ignite == 0:
            powht = (
                powht + current_drive_variables.pinjmw
            )  # if not ignited add the auxiliary power

        # Here the implementation sometimes leaves the accessible regime when p_plasma_rad_mw> powht which is unphysical and
        # is not taken care of by the rad module. We restrict the radiation power here by the heating power:
        physics_variables.p_plasma_rad_mw = max(
            0.0e0, physics_variables.p_plasma_rad_mw
        )

        #  Power to divertor, = (1-stellarator_variables.f_rad)*Psol

        # The SOL radiation needs to be smaller than the physics_variables.p_plasma_rad_mw
        physics_variables.psolradmw = stellarator_variables.f_rad * powht
        physics_variables.pdivt = powht - physics_variables.psolradmw

        # Add SOL Radiation to total
        physics_variables.p_plasma_rad_mw = (
            physics_variables.p_plasma_rad_mw + physics_variables.psolradmw
        )
        # pden_plasma_rad_mw = physics_variables.p_plasma_rad_mw / physics_variables.vol_plasma # this line OVERWRITES the original definition of pden_plasma_rad_mw, probably shouldn't be defined like that as the core does not lose SOL power.

        #  The following line is unphysical, but prevents -ve sqrt argument
        #  Should be obsolete if constraint eqn 17 is turned on (but beware -
        #  this may not be quite correct for stellarators)
        physics_variables.pdivt = max(0.001e0, physics_variables.pdivt)

        #  Power transported to the first wall by escaped alpha particles

        physics_variables.p_fw_alpha_mw = physics_variables.alpha_power_total * (
            1.0e0 - physics_variables.f_alpha_plasma
        )

        # Nominal mean photon wall load
        if physics_variables.iwalld == 1:
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
                        - fwbs_variables.f_a_fw_hcd
                        - fwbs_variables.f_ster_div_single
                    )
                    * physics_variables.p_plasma_rad_mw
                    / build_variables.a_fw_total
                )

        constraint_variables.pflux_fw_rad_max_mw = (
            physics_variables.pflux_fw_rad_mw * constraint_variables.f_fw_rad_max
        )

        physics_variables.rad_fraction_total = physics_variables.p_plasma_rad_mw / (
            physics_variables.f_alpha_plasma * physics_variables.alpha_power_total
            + physics_variables.non_alpha_charged_power
            + physics_variables.p_plasma_ohmic_mw
            + current_drive_variables.pinjmw
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
            physics_variables.alpha_power_total,
            physics_variables.aspect,
            physics_variables.bt,
            physics_variables.nd_ions_total,
            physics_variables.dene,
            physics_variables.dnla,
            physics_variables.eps,
            physics_variables.hfact,
            physics_variables.i_confinement_time,
            physics_variables.ignite,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.non_alpha_charged_power,
            current_drive_variables.pinjmw,
            physics_variables.plasma_current,
            physics_variables.pden_plasma_core_rad_mw,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.ten,
            physics_variables.tin,
            stellarator_variables.iotabar,
            physics_variables.qstar,
            physics_variables.vol_plasma,
            physics_variables.zeff,
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
            physics_variables.qfuel,
            physics_variables.rndfuel,
            physics_variables.t_alpha_confinement,
            physics_variables.f_alpha_energy_confinement,
        ) = self.physics.phyaux(
            physics_variables.aspect,
            physics_variables.dene,
            physics_variables.te,
            physics_variables.nd_fuel_ions,
            physics_variables.fusion_rate_density_total,
            physics_variables.alpha_rate_density_total,
            physics_variables.plasma_current,
            sbar,
            physics_variables.nd_alphas,
            physics_variables.t_energy_confinement,
            physics_variables.vol_plasma,
        )

        # Calculate physics_variables.beta limit. Does nothing atm so commented out
        # call stblim(physics_variables.beta_max)

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
                physics_module.rho_star,
                nu_star_e,
                nu_star_d,
                nu_star_T,
                nu_star_He,
                physics_variables.dnla,
                physics_variables.dnelimt,
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
        dnla,
        dnelimt,
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
            "(dnla)",
            dnla,
        )
        po.ovarre(self.outfile, "Sudo density limit (/m3)", "(dnelimt)", dnelimt)
        po.ovarre(
            self.outfile,
            "Ratio density to sudo limit (1)",
            "(dnla/dnelimt)",
            dnla / dnelimt,
        )

    def calc_neoclassics(self):
        self.neoclassics.init_neoclassics(
            0.6,
            stellarator_configuration.stella_config_epseff,
            stellarator_variables.iotabar,
        )

        q_PROCESS = (
            (
                physics_variables.f_alpha_plasma
                * physics_variables.alpha_power_density_total
                - physics_variables.pden_plasma_core_rad_mw
            )
            * physics_variables.vol_plasma
            / physics_variables.a_plasma_surface
            * impurity_radiation_module.radius_plasma_core_norm
        )
        q_PROCESS_r1 = (
            (
                physics_variables.f_alpha_plasma
                * physics_variables.alpha_power_density_total
                - physics_variables.pden_plasma_core_rad_mw
            )
            * physics_variables.vol_plasma
            / physics_variables.a_plasma_surface
        )

        q_neo = sum(neoclassics_module.q_flux * 1e-6)
        gamma_neo = sum(
            neoclassics_module.gamma_flux * neoclassics_module.temperatures * 1e-6
        )

        total_q_neo = sum(
            neoclassics_module.q_flux * 1e-6
            + neoclassics_module.gamma_flux * neoclassics_module.temperatures * 1e-6
        )

        total_q_neo_e = (
            2
            * 2
            * (
                neoclassics_module.q_flux[0] * 1e-6
                + neoclassics_module.gamma_flux[0]
                * neoclassics_module.temperatures[0]
                * 1e-6
            )
        )

        q_neo_e = neoclassics_module.q_flux[0] * 1e-6
        q_neo_D = neoclassics_module.q_flux[1] * 1e-6
        q_neo_a = neoclassics_module.q_flux[3] * 1e-6
        q_neo_T = neoclassics_module.q_flux[2] * 1e-6

        g_neo_e = (
            neoclassics_module.gamma_flux[0] * 1e-6 * neoclassics_module.temperatures[0]
        )
        g_neo_D = (
            neoclassics_module.gamma_flux[1] * 1e-6 * neoclassics_module.temperatures[1]
        )
        g_neo_a = (
            neoclassics_module.gamma_flux[3] * 1e-6 * neoclassics_module.temperatures[3]
        )
        g_neo_T = (
            neoclassics_module.gamma_flux[2] * 1e-6 * neoclassics_module.temperatures[2]
        )

        dndt_neo_e = neoclassics_module.gamma_flux[0]
        dndt_neo_D = neoclassics_module.gamma_flux[1]
        dndt_neo_a = neoclassics_module.gamma_flux[3]
        dndt_neo_T = neoclassics_module.gamma_flux[2]

        dndt_neo_fuel = (
            (dndt_neo_D + dndt_neo_T)
            * physics_variables.a_plasma_surface
            * impurity_radiation_module.radius_plasma_core_norm
        )
        dmdt_neo_fuel = (
            dndt_neo_fuel * physics_variables.m_fuel_amu * constants.proton_mass * 1.0e6
        )  # mg
        dmdt_neo_fuel_from_e = (
            4
            * dndt_neo_e
            * physics_variables.a_plasma_surface
            * impurity_radiation_module.radius_plasma_core_norm
            * physics_variables.m_fuel_amu
            * constants.proton_mass
            * 1.0e6
        )  # kg

        chi_neo_e = -(
            neoclassics_module.q_flux[0]
            + neoclassics_module.gamma_flux[0] * neoclassics_module.temperatures[0]
        ) / (
            neoclassics_module.densities[0] * neoclassics_module.dr_temperatures[0]
            + neoclassics_module.temperatures[0] * neoclassics_module.dr_densities[0]
        )

        chi_PROCESS_e = self.st_calc_eff_chi()

        nu_star_e = neoclassics_module.nu_star_averaged[0]
        nu_star_d = neoclassics_module.nu_star_averaged[1]
        nu_star_T = neoclassics_module.nu_star_averaged[2]
        nu_star_He = neoclassics_module.nu_star_averaged[3]

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
            * st.f_r
            * (
                impurity_radiation_module.radius_plasma_core_norm
                * physics_variables.rminor
                / stellarator_configuration.stella_config_rminor_ref
            )
            ** 2
        )
        surfacescaling = (
            physics_variables.a_plasma_surface
            * st.f_r
            * (
                impurity_radiation_module.radius_plasma_core_norm
                * physics_variables.rminor
                / stellarator_configuration.stella_config_rminor_ref
            )
        )

        nominator = (
            physics_variables.f_alpha_plasma
            * physics_variables.alpha_power_density_total
            - physics_variables.pden_plasma_core_rad_mw
        ) * volscaling

        # in fortran there was a 0*alphan term which I have removed for obvious reasons
        # the following comment seems to describe this?
        # "include alphan if chi should be incorporate density gradients too"
        # but the history can be consulted if required (23/11/22 TN)
        denominator = (
            (
                3
                * physics_variables.ne0
                * constants.electron_charge
                * physics_variables.te0
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
            current_drive_variables.echpwr = current_drive_variables.pheat
            current_drive_variables.pinjimw = 0
            current_drive_variables.pinjemw = current_drive_variables.echpwr
            current_drive_variables.etacd = current_drive_variables.etaech
            current_drive_variables.pinjwp = (
                current_drive_variables.pinjimw + current_drive_variables.pinjemw
            ) / current_drive_variables.etacd
        elif stellarator_variables.isthtr == 2:
            current_drive_variables.plhybd = current_drive_variables.pheat
            current_drive_variables.pinjimw = 0
            current_drive_variables.pinjemw = current_drive_variables.plhybd
            current_drive_variables.etacd = current_drive_variables.etalh
            current_drive_variables.pinjwp = (
                current_drive_variables.pinjimw + current_drive_variables.pinjemw
            ) / current_drive_variables.etacd
        elif stellarator_variables.isthtr == 3:
            (
                effnbss,
                fpion,
                current_drive_variables.nbshinef,
            ) = self.current_drive.culnbi()
            current_drive_variables.pnbeam = current_drive_variables.pheat * (
                1 - current_drive_variables.forbitloss
            )
            current_drive_variables.porbitlossmw = (
                current_drive_variables.pheat * current_drive_variables.forbitloss
            )
            current_drive_variables.pinjimw = current_drive_variables.pnbeam * fpion
            current_drive_variables.pinjemw = current_drive_variables.pnbeam * (
                1 - fpion
            )
            current_drive_variables.etacd = current_drive_variables.etanbi
            current_drive_variables.pinjwp = (
                current_drive_variables.pinjimw + current_drive_variables.pinjemw
            ) / current_drive_variables.etacd
        else:
            error_handling.idiags[0] = stellarator_variables.isthtr
            error_handling.report_error(107)

        #  Total injected power

        current_drive_variables.pinjmw = (
            current_drive_variables.pinjemw + current_drive_variables.pinjimw
        )

        #  Calculate neutral beam current

        if abs(current_drive_variables.pnbeam) > 1e-8:
            current_drive_variables.beam_current = (
                1e-3
                * (current_drive_variables.pnbeam * 1e6)
                / current_drive_variables.beam_energy
            )
        else:
            current_drive_variables.beam_current = 0

        #  Ratio of fusion to input (injection+ohmic) power

        if (
            abs(
                current_drive_variables.pinjmw
                + current_drive_variables.porbitlossmw
                + physics_variables.p_plasma_ohmic_mw
            )
            < 1e-6
        ):
            current_drive_variables.bigq = 1e18
        else:
            current_drive_variables.bigq = physics_variables.fusion_power / (
                current_drive_variables.pinjmw
                + current_drive_variables.porbitlossmw
                + physics_variables.p_plasma_ohmic_mw
            )

        if output:
            po.oheadr(self.outfile, "Auxiliary Heating System")

            if stellarator_variables.isthtr == 1:
                po.ocmmnt(self.outfile, "Electron Cyclotron Resonance Heating")
            elif stellarator_variables.isthtr == 2:
                po.ocmmnt(self.outfile, "Lower Hybrid Heating")
            elif stellarator_variables.isthtr == 3:
                po.ocmmnt(self.outfile, "Neutral Beam Injection Heating")

            if physics_variables.ignite == 1:
                po.ocmmnt(
                    self.outfile,
                    "Ignited plasma; injected power only used for start-up phase",
                )

            po.oblnkl(self.outfile)

            po.ovarre(
                self.outfile,
                "Auxiliary power supplied to plasma (MW)",
                "(pheat)",
                current_drive_variables.pheat,
            )
            po.ovarre(
                self.outfile,
                "Fusion gain factor Q",
                "(bigq)",
                current_drive_variables.bigq,
            )

            if abs(current_drive_variables.pnbeam) > 1e-8:
                po.ovarre(
                    self.outfile,
                    "Neutral beam energy (KEV)",
                    "(enbeam)",
                    current_drive_variables.enbeam,
                )
                po.ovarre(
                    self.outfile,
                    "Neutral beam current (A)",
                    "(beam_current)",
                    current_drive_variables.beam_current,
                )
                po.ovarre(
                    self.outfile, "Fraction of beam energy to ions", "(fpion)", fpion
                )
                po.ovarre(
                    self.outfile,
                    "Neutral beam shine-through fraction",
                    "(nbshinef)",
                    current_drive_variables.nbshinef,
                )
                po.ovarre(
                    self.outfile,
                    "Neutral beam orbit loss power (MW)",
                    "(porbitlossmw)",
                    current_drive_variables.porbitlossmw,
                )
                po.ovarre(
                    self.outfile,
                    "Beam duct shielding thickness (m)",
                    "(nbshield)",
                    current_drive_variables.nbshield,
                )
                po.ovarre(
                    self.outfile,
                    "R injection tangent / R-major",
                    "(frbeam)",
                    current_drive_variables.frbeam,
                )
                po.ovarre(
                    self.outfile,
                    "Beam centreline tangency radius (m)",
                    "(rtanbeam)",
                    current_drive_variables.rtanbeam,
                )
                po.ovarre(
                    self.outfile,
                    "Maximum possible tangency radius (m)",
                    "(rtanmax)",
                    current_drive_variables.rtanmax,
                )
                po.ovarre(
                    self.outfile,
                    "Beam decay lengths to centre",
                    "(taubeam)",
                    current_drive_variables.taubeam,
                )


class Neoclassics:
    @property
    def no_roots(self):
        return neoclassics_module.roots.shape[0]

    def init_neoclassics(self, r_effin, eps_effin, iotain):
        """Constructor of the neoclassics object from the effective radius,
        epsilon effective and iota only.
        """
        (
            neoclassics_module.densities,
            neoclassics_module.temperatures,
            neoclassics_module.dr_densities,
            neoclassics_module.dr_temperatures,
        ) = self.init_profile_values_from_PROCESS(r_effin)
        neoclassics_module.roots = np.array([
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
        neoclassics_module.weights = np.array([
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

        neoclassics_module.kt = self.neoclassics_calc_KT()
        neoclassics_module.nu = self.neoclassics_calc_nu()
        neoclassics_module.nu_star = self.neoclassics_calc_nu_star()
        neoclassics_module.nu_star_averaged = self.neoclassics_calc_nu_star_fromT(
            iotain
        )
        neoclassics_module.vd = self.neoclassics_calc_vd()

        neoclassics_module.d11_plateau = self.neoclassics_calc_D11_plateau()

        neoclassics_module.d11_mono = self.neoclassics_calc_d11_mono(
            eps_effin
        )  # for using epseff

        neoclassics_module.d111 = self.calc_integrated_radial_transport_coeffs(index=1)
        neoclassics_module.d112 = self.calc_integrated_radial_transport_coeffs(index=2)
        neoclassics_module.d113 = self.calc_integrated_radial_transport_coeffs(index=3)

        neoclassics_module.gamma_flux = self.neoclassics_calc_gamma_flux(
            neoclassics_module.densities,
            neoclassics_module.temperatures,
            neoclassics_module.dr_densities,
            neoclassics_module.dr_temperatures,
        )
        neoclassics_module.q_flux = self.neoclassics_calc_q_flux()

    def init_profile_values_from_PROCESS(self, rho):
        """Initializes the profile_values object from PROCESS' parabolic profiles"""
        tempe = physics_variables.te0 * (1 - rho**2) ** physics_variables.alphat * KEV
        tempT = physics_variables.ti0 * (1 - rho**2) ** physics_variables.alphat * KEV
        tempD = physics_variables.ti0 * (1 - rho**2) ** physics_variables.alphat * KEV
        tempa = physics_variables.ti0 * (1 - rho**2) ** physics_variables.alphat * KEV

        dense = physics_variables.ne0 * (1 - rho**2) ** physics_variables.alphan
        densT = (
            (1 - physics_variables.f_deuterium)
            * physics_variables.ni0
            * (1 - rho**2) ** physics_variables.alphan
        )
        densD = (
            physics_variables.f_deuterium
            * physics_variables.ni0
            * (1 - rho**2) ** physics_variables.alphan
        )
        densa = (
            physics_variables.nd_alphas
            * (1 + physics_variables.alphan)
            * (1 - rho**2) ** physics_variables.alphan
        )

        # Derivatives in real space
        dr_tempe = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.te0
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )
        dr_tempT = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.ti0
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )
        dr_tempD = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.ti0
            * rho
            * (1.0 - rho**2) ** (physics_variables.alphat - 1.0)
            * physics_variables.alphat
            * KEV
        )
        dr_tempa = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * physics_variables.ti0
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
            * physics_variables.ne0
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )
        dr_densT = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * (1 - physics_variables.f_deuterium)
            * physics_variables.ni0
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )
        dr_densD = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * physics_variables.f_deuterium
            * physics_variables.ni0
            * (1.0 - rho**2) ** (physics_variables.alphan - 1.0)
            * physics_variables.alphan
        )
        dr_densa = (
            -2.0
            * 1.0
            / physics_variables.rminor
            * rho
            * physics_variables.nd_alphas
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
        k = np.repeat((neoclassics_module.roots / KEV)[:, np.newaxis], 4, axis=1)

        return (k * neoclassics_module.temperatures).T

    def neoclassics_calc_nu(self):
        """Calculates the collision frequency"""
        mass = np.array([
            constants.electron_mass,
            constants.proton_mass * 2.0,
            constants.proton_mass * 3.0,
            constants.proton_mass * 4.0,
        ])
        z = np.array([-1.0, 1.0, 1.0, 2.0]) * constants.electron_charge

        # transform the temperature back in eV
        # Formula from L. Spitzer.Physics of fully ionized gases.  Interscience, New York, 1962
        lnlambda = (
            32.2
            - 1.15 * np.log10(neoclassics_module.densities[0])
            + 2.3
            * np.log10(neoclassics_module.temperatures[0] / constants.electron_charge)
        )

        neoclassics_calc_nu = np.zeros((4, self.no_roots), order="F")

        for j in range(4):
            for i in range(self.no_roots):
                x = neoclassics_module.roots[i]
                for k in range(4):
                    xk = (
                        (mass[k] / mass[j])
                        * (
                            neoclassics_module.temperatures[j]
                            / neoclassics_module.temperatures[k]
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
                    v = np.sqrt(2.0 * x * neoclassics_module.temperatures[j] / mass[j])
                    neoclassics_calc_nu[j, i] = neoclassics_calc_nu[
                        j, i
                    ] + neoclassics_module.densities[k] * (
                        z[j] * z[k]
                    ) ** 2 * lnlambda * phixmgx / (
                        4.0 * np.pi * constants.epsilon0**2 * mass[j] ** 2 * v**3
                    )

        return neoclassics_calc_nu

    def neoclassics_calc_nu_star(self):
        """Calculates the normalized collision frequency"""
        k = np.repeat(neoclassics_module.roots[:, np.newaxis], 4, axis=1)
        kk = (k * neoclassics_module.temperatures).T

        mass = np.array([
            constants.electron_mass,
            constants.proton_mass * 2.0,
            constants.proton_mass * 3.0,
            constants.proton_mass * 4.0,
        ])

        v = np.empty((4, self.no_roots))
        v[0, :] = constants.speed_light * np.sqrt(
            1.0 - (kk[0, :] / (mass[0] * constants.speed_light**2) + 1) ** (-1)
        )
        v[1, :] = constants.speed_light * np.sqrt(
            1.0 - (kk[1, :] / (mass[1] * constants.speed_light**2) + 1) ** (-1)
        )
        v[2, :] = constants.speed_light * np.sqrt(
            1.0 - (kk[2, :] / (mass[2] * constants.speed_light**2) + 1) ** (-1)
        )
        v[3, :] = constants.speed_light * np.sqrt(
            1.0 - (kk[3, :] / (mass[3] * constants.speed_light**2) + 1) ** (-1)
        )

        return (
            physics_variables.rmajor
            * neoclassics_module.nu
            / (neoclassics_module.iota * v)
        )

    def neoclassics_calc_nu_star_fromT(self, iota):
        """Calculates the collision frequency"""
        temp = (
            np.array([
                physics_variables.te,
                physics_variables.ti,
                physics_variables.ti,
                physics_variables.ti,
            ])
            * KEV
        )
        density = np.array([
            physics_variables.dene,
            physics_variables.nd_fuel_ions * physics_variables.f_deuterium,
            physics_variables.nd_fuel_ions * (1 - physics_variables.f_deuterium),
            physics_variables.nd_alphas,
        ])

        mass = np.array([
            constants.electron_mass,
            constants.proton_mass * 2.0,
            constants.proton_mass * 3.0,
            constants.proton_mass * 4.0,
        ])
        z = np.array([-1.0, 1.0, 1.0, 2.0]) * constants.electron_charge

        # transform the temperature back in eV
        # Formula from L. Spitzer.Physics of fully ionized gases.  Interscience, New York, 1962
        lnlambda = (
            32.2
            - 1.15 * np.log10(density[0])
            + 2.3 * np.log10(temp[0] / constants.electron_charge)
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
                    / (4.0 * np.pi * constants.epsilon0**2 * mass[j] ** 2 * v**4)
                    * physics_variables.rmajor
                    / iota
                )
        return neoclassics_calc_nu_star_fromT

    def neoclassics_calc_vd(self):
        vde = (
            neoclassics_module.roots
            * neoclassics_module.temperatures[0]
            / (
                constants.electron_charge
                * physics_variables.rmajor
                * physics_variables.bt
            )
        )
        vdD = (
            neoclassics_module.roots
            * neoclassics_module.temperatures[1]
            / (
                constants.electron_charge
                * physics_variables.rmajor
                * physics_variables.bt
            )
        )
        vdT = (
            neoclassics_module.roots
            * neoclassics_module.temperatures[2]
            / (
                constants.electron_charge
                * physics_variables.rmajor
                * physics_variables.bt
            )
        )
        vda = (
            neoclassics_module.roots
            * neoclassics_module.temperatures[3]
            / (
                2.0
                * constants.electron_charge
                * physics_variables.rmajor
                * physics_variables.bt
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
            constants.electron_mass,
            constants.proton_mass * 2.0,
            constants.proton_mass * 3.0,
            constants.proton_mass * 4.0,
        ])

        v = np.empty((4, self.no_roots))
        v[0, :] = constants.speed_light * np.sqrt(
            1.0
            - (neoclassics_module.kt[0, :] / (mass[0] * constants.speed_light**2) + 1)
            ** (-1)
        )
        v[1, :] = constants.speed_light * np.sqrt(
            1.0
            - (neoclassics_module.kt[1, :] / (mass[1] * constants.speed_light**2) + 1)
            ** (-1)
        )
        v[2, :] = constants.speed_light * np.sqrt(
            1.0
            - (neoclassics_module.kt[2, :] / (mass[2] * constants.speed_light**2) + 1)
            ** (-1)
        )
        v[3, :] = constants.speed_light * np.sqrt(
            1.0
            - (neoclassics_module.kt[3, :] / (mass[3] * constants.speed_light**2) + 1)
            ** (-1)
        )

        return (
            np.pi
            / 4.0
            * neoclassics_module.vd**2
            * physics_variables.rmajor
            / neoclassics_module.iota
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
            * neoclassics_module.vd**2
            / neoclassics_module.nu
        )

    def calc_integrated_radial_transport_coeffs(self, index: int):
        """Calculates the integrated radial transport coefficients (index `index`)
        It uses Gauss laguerre integration
        https://en.wikipedia.org/wiki/Gauss%E2%80%93Laguerre_quadrature
        """
        return np.sum(
            2.0
            / np.sqrt(np.pi)
            * neoclassics_module.d11_mono
            * neoclassics_module.roots ** (index - 0.5)
            * neoclassics_module.weights,
            axis=1,
        )

    def neoclassics_calc_gamma_flux(
        self, densities, temperatures, dr_densities, dr_temperatures
    ):
        """Calculates the Energy flux by neoclassical particle transport"""

        z = np.array([-1.0, 1.0, 1.0, 2.0])

        return (
            -densities
            * neoclassics_module.d111
            * (
                (dr_densities / densities - z * neoclassics_module.er / temperatures)
                + (neoclassics_module.d112 / neoclassics_module.d111 - 3.0 / 2.0)
                * dr_temperatures
                / temperatures
            )
        )

    def neoclassics_calc_q_flux(self):
        """Calculates the Energy flux by neoclassicsal energy transport"""

        z = np.array([-1.0, 1.0, 1.0, 2.0])

        return (
            -neoclassics_module.densities
            * neoclassics_module.temperatures
            * neoclassics_module.d112
            * (
                (
                    neoclassics_module.dr_densities / neoclassics_module.densities
                    - z * neoclassics_module.er / neoclassics_module.temperatures
                )
                + (neoclassics_module.d113 / neoclassics_module.d112 - 3.0 / 2.0)
                * neoclassics_module.dr_temperatures
                / neoclassics_module.temperatures
            )
        )


def init_stellarator_variables():
    stellarator_variables.istell = 0
    stellarator_variables.bmn = 1e-3
    stellarator_variables.f_asym = 1.0
    stellarator_variables.f_rad = 0.85
    stellarator_variables.f_w = 0.5
    stellarator_variables.fdivwet = 0.333333333333333
    stellarator_variables.flpitch = 1e-3
    stellarator_variables.hportamax = 0.0
    stellarator_variables.hportpmax = 0.0
    stellarator_variables.hporttmax = 0.0
    stellarator_variables.iotabar = 1.0
    stellarator_variables.isthtr = 3
    stellarator_variables.m_res = 5
    stellarator_variables.n_res = 5
    stellarator_variables.shear = 0.5
    stellarator_variables.vportamax = 0.0
    stellarator_variables.vportpmax = 0.0
    stellarator_variables.vporttmax = 0.0
    stellarator_variables.max_gyrotron_frequency = 1.0e9
    stellarator_variables.te0_ecrh_achievable = 1.0e2


def init_stellarator_module():
    st.first_call = True
    st.first_call_stfwbs = True
    st.f_n = 0.0
    st.f_r = 0.0
    st.f_a = 0.0
    st.f_b = 0.0
    st.f_i = 0.0
