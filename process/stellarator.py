import logging
from copy import copy
import numpy as np

from process.fortran import (
    constants,
    physics_module as ph,
    stellarator_module as st,
    process_output as po,
    physics_variables,
    physics_module,
    current_drive_variables,
    tfcoil_variables,
    stellarator_configuration,
    stellarator_variables,
    numerics,
    build_variables,
    fwbs_variables,
    heat_transport_variables,
    structure_variables,
    divertor_variables,
    cost_variables,
    fwbs_module,
    error_handling,
    constraint_variables,
    rebco_variables,
    maths_library,
    profiles_module,
    physics_functions_module,
    neoclassics_module,
    impurity_radiation_module,
)
import process.superconductors as superconductors
import process.physics_functions as physics_funcs
from process.coolprop_interface import FluidProperties

logger = logging.getLogger(__name__)
# Logging handler for console output
s_handler = logging.StreamHandler()
s_handler.setLevel(logging.ERROR)
logger.addHandler(s_handler)


class Stellarator:
    """Module containing stellarator routines
    author: P J Knight, CCFE, Culham Science Centre
    N/A
    This module contains routines for calculating the
    parameters of the first wall, blanket and shield components
    of a fusion power plant.

    AEA FUS 251: A User's Guide to the PROCESS Systems Code
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

    def run(self, output: bool):
        """Routine to call the physics and engineering modules
        relevant to stellarators
        author: P J Knight, CCFE, Culham Science Centre
        author: F Warmer, IPP Greifswald

        This routine is the caller for the stellarator models.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code

        :param output: indicate whether output should be written to the output file, or not
        :type output: boolean
        """

        if output:
            self.costs.costs(output=True)
            self.availability.run(output=True)
            ph.outplas(self.outfile)
            self.stigma()
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
        self.costs.costs(output=False)

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

    def stigma(self):
        """Routine to calculate ignition margin at the final point
        with different stellarator confinement time scaling laws
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        This routine calculates the ignition margin at the final
        point with different stellarator confinement time scaling laws
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        po.osubhd(self.outfile, "Confinement times, and required H-factors :")

        po.write(
            self.outfile,
            f"{' '*5}scaling law{' '*30}confinement time (s){' '*55}H-factor for",
        )
        po.write(self.outfile, f"{' '*34}for H = 2{' '*54}power balance")

        #  Label stellarator scaling laws (update if more are added)

        istlaw = [21, 22, 23, 37, 38]

        #  Calculate power balances for all stellarator scaling laws
        #  assuming H = 2

        for iisc, i in enumerate(istlaw):
            (
                physics_variables.kappaa,
                physics_variables.ptrepv,
                physics_variables.ptripv,
                physics_variables.tauee,
                physics_variables.tauei,
                physics_variables.taueff,
                physics_variables.powerht,
            ) = physics_module.pcond(
                physics_variables.afuel,
                physics_variables.palpmw,
                physics_variables.aspect,
                physics_variables.bt,
                physics_variables.dnitot,
                physics_variables.dene,
                physics_variables.dnla,
                physics_variables.eps,
                2.0,
                physics_variables.iinvqd,
                physics_variables.isc,
                physics_variables.ignite,
                physics_variables.kappa,
                physics_variables.kappa95,
                physics_variables.pchargemw,
                current_drive_variables.pinjmw,
                physics_variables.plascur,
                physics_variables.pcoreradpv,
                physics_variables.rmajor,
                physics_variables.rminor,
                physics_variables.te,
                physics_variables.ten,
                physics_variables.tin,
                physics_variables.q,
                physics_variables.qstar,
                physics_variables.vol,
                physics_variables.xarea,
                physics_variables.zeff,
            )

            physics_variables.hfac[iisc] = physics_module.fhfac(i)

    def stnewconfig(self):
        """author: J Lion, IPP Greifswald
        Routine to initialise the stellarator configuration

        Routine to initialise the stellarator configuration.
        This routine is called right before the calculation and could
        in principle overwrite variables from the input file.
        It overwrites rminor with rmajor and aspect ratio e.g.
        """

        stellarator_configuration.new_stella_config(stellarator_variables.istell)

        # If physics_variables.aspect ratio is not in numerics.ixc set it to default value
        # Or when you call it the first time
        if 1 not in numerics.ixc:
            physics_variables.aspect = (
                stellarator_configuration.stella_config_aspect_ref
            )

        # Set the physics_variables.rminor radius as result here.
        physics_variables.rminor = physics_variables.rmajor / physics_variables.aspect
        physics_variables.eps = 1.0e0 / physics_variables.aspect

        tfcoil_variables.n_tf = (
            stellarator_configuration.stella_config_coilspermodule
            * stellarator_configuration.stella_config_symmetry
        )  # This overwrites tfcoil_variables.n_tf in input file.

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
        st.f_n = tfcoil_variables.n_tf / (
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
        physics_variables.vol = (
            st.f_r * st.f_a**2 * stellarator_configuration.stella_config_plasma_volume
        )

        # Plasma surface scaled from effective parameter:
        physics_variables.sarea = (
            st.f_r * st.f_a * stellarator_configuration.stella_config_plasma_surface
        )

        # Plasma cross section area. Approximated
        physics_variables.xarea = (
            np.pi * physics_variables.rminor * physics_variables.rminor
        )  # average, could be calculated for every toroidal angle if desired

        #  physics_variables.sareao is retained only for obsolescent fispact calculation...

        #  Cross-sectional area, averaged over toroidal angle
        physics_variables.sareao = (
            0.5e0 * physics_variables.sarea
        )  # Used only in the divertor model; approximate as for tokamaks

    def stopt(self, output: bool):
        """Routine to reiterate the physics loop
        author: J Lion, IPP Greifswald
        None
        This routine reiterates some physics modules.
        """

        physics_variables.dnelimt = self.stdlim(
            physics_variables.bt,
            physics_variables.powerht,
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
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        if fwbs_variables.blktmodel > 0:
            build_variables.blnkith = (
                build_variables.blbuith
                + build_variables.blbmith
                + build_variables.blbpith
            )
            build_variables.blnkoth = (
                build_variables.blbuoth
                + build_variables.blbmoth
                + build_variables.blbpoth
            )
            build_variables.shldtth = 0.5e0 * (
                build_variables.shldith + build_variables.shldoth
            )

        #  Top/bottom blanket thickness

        build_variables.blnktth = 0.5e0 * (
            build_variables.blnkith + build_variables.blnkoth
        )

        # First Wall
        build_variables.fwith = (
            2.0e0 * fwbs_variables.afw + 2.0e0 * fwbs_variables.fw_wall
        )
        build_variables.fwoth = build_variables.fwith

        build_variables.bore = physics_variables.rmajor - (
            build_variables.ohcth
            + build_variables.gapoh
            + build_variables.tfcth
            + build_variables.gapds
            + build_variables.d_vv_in
            + build_variables.shldith
            + build_variables.blnkith
            + build_variables.fwith
            + build_variables.scrapli
            + physics_variables.rminor
        )

        #  Radial build to centre of plasma (should be equal to physics_variables.rmajor)
        build_variables.rbld = (
            build_variables.bore
            + build_variables.ohcth
            + build_variables.gapoh
            + build_variables.tfcth
            + build_variables.gapds
            + build_variables.d_vv_in
            + build_variables.shldith
            + build_variables.blnkith
            + build_variables.fwith
            + build_variables.scrapli
            + physics_variables.rminor
        )

        # Bc stellarators cannot scale physics_variables.rminor reasonably well an additional constraint equation is required,
        # that ensures that there is enough space between coils and plasma.
        build_variables.required_radial_space = (
            build_variables.tfcth / 2.0e0
            + build_variables.gapds
            + build_variables.d_vv_in
            + build_variables.shldith
            + build_variables.blnkith
            + build_variables.fwith
            + build_variables.scrapli
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
            - build_variables.scrapli
            - build_variables.fwith
            - build_variables.blnkith
            - build_variables.shldith
        )

        #  Radius to outer edge of outboard shield
        build_variables.rsldo = (
            physics_variables.rmajor
            + physics_variables.rminor
            + build_variables.scraplo
            + build_variables.fwoth
            + build_variables.blnkoth
            + build_variables.shldoth
        )

        #  Thickness of outboard TF coil legs
        build_variables.tfthko = build_variables.tfcth

        #  Radius to centre of outboard TF coil legs

        build_variables.gapsto = build_variables.gapomin
        build_variables.r_tf_outboard_mid = (
            build_variables.rsldo
            + build_variables.d_vv_out
            + build_variables.gapsto
            + 0.5e0 * build_variables.tfthko
        )

        #  Height to inside edge of TF coil
        #  Roughly equal to average of (inboard build from TF coil to plasma
        #  centre) and (outboard build from plasma centre to TF coil)

        build_variables.hmax = 0.5e0 * (
            (
                build_variables.gapds
                + build_variables.d_vv_in
                + build_variables.shldith
                + build_variables.blnkith
                + build_variables.fwith
                + build_variables.scrapli
                + physics_variables.rminor
            )
            + (
                physics_variables.rminor
                + build_variables.scraplo
                + build_variables.fwoth
                + build_variables.blnkoth
                + build_variables.shldoth
                + build_variables.d_vv_out
                + build_variables.gapsto
            )
        )

        #  Outer divertor strike point radius, set equal to major radius

        build_variables.rspo = physics_variables.rmajor

        #  First wall area: scales with minor radius

        # Average minor radius of the first wall
        awall = physics_variables.rminor + 0.5e0 * (
            build_variables.scrapli + build_variables.scraplo
        )
        build_variables.fwarea = (
            physics_variables.sarea * awall / physics_variables.rminor
        )

        if heat_transport_variables.ipowerflow == 0:
            build_variables.fwarea = (
                1.0e0 - fwbs_variables.fhole
            ) * build_variables.fwarea
        else:
            build_variables.fwarea = (
                1.0e0 - fwbs_variables.fhole - fwbs_variables.fdiv - fwbs_variables.fhcd
            ) * build_variables.fwarea

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
                build_variables.bore + build_variables.ohcth + build_variables.gapoh
            )
            radius = radius + drbild
            po.obuild(self.outfile, "Machine bore", drbild, radius, "(bore)")
            po.ovarre(
                self.outfile, "Machine build_variables.bore (m)", "(bore)", drbild
            )

            radius = radius + build_variables.tfcth
            po.obuild(
                self.outfile,
                "Coil inboard leg",
                build_variables.tfcth,
                radius,
                "(tfcth)",
            )
            po.ovarre(
                self.outfile, "Coil inboard leg (m)", "(deltf)", build_variables.tfcth
            )

            radius = radius + build_variables.gapds
            po.obuild(self.outfile, "Gap", build_variables.gapds, radius, "(gapds)")
            po.ovarre(self.outfile, "Gap (m)", "(gapds)", build_variables.gapds)

            radius = radius + build_variables.d_vv_in
            po.obuild(
                self.outfile,
                "Vacuum vessel",
                build_variables.d_vv_in,
                radius,
                "(d_vv_in)",
            )
            po.ovarre(
                self.outfile,
                "Vacuum vessel radial thickness (m)",
                "(d_vv_in)",
                build_variables.d_vv_in,
            )

            radius = radius + build_variables.shldith
            po.obuild(
                self.outfile,
                "Inboard shield",
                build_variables.shldith,
                radius,
                "(shldith)",
            )
            po.ovarre(
                self.outfile,
                "Inner radiation shield radial thickness (m)",
                "(shldith)",
                build_variables.shldith,
            )

            radius = radius + build_variables.blnkith
            po.obuild(
                self.outfile,
                "Inboard blanket",
                build_variables.blnkith,
                radius,
                "(blnkith)",
            )
            po.ovarre(
                self.outfile,
                "Inboard blanket radial thickness (m)",
                "(blnkith)",
                build_variables.blnkith,
            )

            radius = radius + build_variables.fwith
            po.obuild(
                self.outfile,
                "Inboard first wall",
                build_variables.fwith,
                radius,
                "(fwith)",
            )
            po.ovarre(
                self.outfile,
                "Inboard first wall radial thickness (m)",
                "(fwith)",
                build_variables.fwith,
            )

            radius = radius + build_variables.scrapli
            po.obuild(
                self.outfile,
                "Inboard scrape-off",
                build_variables.scrapli,
                radius,
                "(scrapli)",
            )
            po.ovarre(
                self.outfile,
                "Inboard scrape-off radial thickness (m)",
                "(scrapli)",
                build_variables.scrapli,
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

            radius = radius + build_variables.scraplo
            po.obuild(
                self.outfile,
                "Outboard scrape-off",
                build_variables.scraplo,
                radius,
                "(scraplo)",
            )
            po.ovarre(
                self.outfile,
                "Outboard scrape-off radial thickness (m)",
                "(scraplo)",
                build_variables.scraplo,
            )

            radius = radius + build_variables.fwoth
            po.obuild(
                self.outfile,
                "Outboard first wall",
                build_variables.fwoth,
                radius,
                "(fwoth)",
            )
            po.ovarre(
                self.outfile,
                "Outboard first wall radial thickness (m)",
                "(fwoth)",
                build_variables.fwoth,
            )

            radius = radius + build_variables.blnkoth
            po.obuild(
                self.outfile,
                "Outboard blanket",
                build_variables.blnkoth,
                radius,
                "(blnkoth)",
            )
            po.ovarre(
                self.outfile,
                "Outboard blanket radial thickness (m)",
                "(blnkoth)",
                build_variables.blnkoth,
            )

            radius = radius + build_variables.shldoth
            po.obuild(
                self.outfile,
                "Outboard shield",
                build_variables.shldoth,
                radius,
                "(shldoth)",
            )
            po.ovarre(
                self.outfile,
                "Outer radiation shield radial thickness (m)",
                "(shldoth)",
                build_variables.shldoth,
            )

            radius = radius + build_variables.d_vv_out
            po.obuild(
                self.outfile,
                "Vacuum vessel",
                build_variables.d_vv_out,
                radius,
                "(d_vv_out)",
            )

            radius = radius + build_variables.gapsto
            po.obuild(self.outfile, "Gap", build_variables.gapsto, radius, "(gapsto)")
            po.ovarre(self.outfile, "Gap (m)", "(gapsto)", build_variables.gapsto)

            radius = radius + build_variables.tfthko
            po.obuild(
                self.outfile,
                "Coil outboard leg",
                build_variables.tfthko,
                radius,
                "(tfthko)",
            )
            po.ovarre(
                self.outfile,
                "Coil outboard leg radial thickness (m)",
                "(tfthko)",
                build_variables.tfthko,
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
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
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
        M_struc = 1.3483e0 * (1000.0e0 * tfcoil_variables.estotftgj) ** 0.7821e0
        msupstr = 1000.0e0 * M_struc  # kg

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
            tfcoil_variables.whttf
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
        R = physics_variables.rmajor
        P_div = physics_variables.pdivt
        alpha = divertor_variables.anginc
        xi_p = divertor_variables.xpertin
        T_scrape = divertor_variables.tdiv

        #  Scrape-off temperature in Joules

        E = T_scrape * constants.echarge

        #  Sound speed of particles (m/s)

        c_s = np.sqrt(E / (physics_variables.afuel * constants.umass))

        #  Island size (m)

        w_r = 4.0e0 * np.sqrt(
            stellarator_variables.bmn
            * R
            / (stellarator_variables.shear * stellarator_variables.n_res)
        )

        #  Perpendicular (to plate) distance from X-point to divertor plate (m)

        Delta = stellarator_variables.f_w * w_r

        #  Length 'along' plasma (m)

        l_p = (
            2 * np.pi * R * (stellarator_variables.m_res) / stellarator_variables.n_res
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

        q_div = stellarator_variables.f_asym * (P_div / a_eff)

        #  Transfer to global variables

        divertor_variables.hldiv = q_div
        divertor_variables.divsur = darea

        fwbs_variables.fdiv = darea / build_variables.fwarea

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

        fwbs_variables.whtblkt = fwbs_variables.volblkt * fwbs_variables.densbreed
        self.hcpb.nuclear_heating_blanket()

        # Heating of the magnets
        self.hcpb.nuclear_heating_magnets(False)

        # Rough estimate of TF coil volume used, assuming 25% of the total
        # TF coil perimeter is inboard, 75% outboard
        tf_volume = (
            0.25 * tfcoil_variables.tfleng * tfcoil_variables.tfareain
            + 0.75
            * tfcoil_variables.tfleng
            * tfcoil_variables.arealeg
            * tfcoil_variables.n_tf
        )

        fwbs_variables.ptfnucpm3 = fwbs_variables.ptfnuc / tf_volume

        # heating of the shield
        self.hcpb.nuclear_heating_shield()

        # Energy multiplication factor
        fwbs_variables.emult = 1.269

        # Tritium breeding ratio
        fwbs_variables.tbr = self.hcpb.tbr_shimwell(
            fwbs_variables.volblkt, fwbs_variables.li6enrich, 1
        )

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
        ) = fwbs_module.sctfcoil_nuclear_heating_iter90()

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
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        fwbs_variables.fwlife = min(
            cost_variables.abktflnc / physics_variables.wallmw, cost_variables.tlife
        )

        #  First wall inboard, outboard areas (assume 50% of total each)
        build_variables.fwareaib = 0.5e0 * build_variables.fwarea
        build_variables.fwareaob = 0.5e0 * build_variables.fwarea

        #  Blanket volume; assume that its surface area is scaled directly from the
        #  plasma surface area.
        #  Uses fwbs_variables.fhole etc. to take account of gaps due to ports etc.

        r1 = physics_variables.rminor + 0.5e0 * (
            build_variables.scrapli
            + build_variables.fwith
            + build_variables.scraplo
            + build_variables.fwoth
        )
        if heat_transport_variables.ipowerflow == 0:
            build_variables.blarea = (
                physics_variables.sarea
                * r1
                / physics_variables.rminor
                * (1.0e0 - fwbs_variables.fhole)
            )
        else:
            build_variables.blarea = (
                physics_variables.sarea
                * r1
                / physics_variables.rminor
                * (
                    1.0e0
                    - fwbs_variables.fhole
                    - fwbs_variables.fdiv
                    - fwbs_variables.fhcd
                )
            )

        build_variables.blareaib = 0.5e0 * build_variables.blarea
        build_variables.blareaob = 0.5e0 * build_variables.blarea

        fwbs_variables.volblkti = build_variables.blareaib * build_variables.blnkith
        fwbs_variables.volblkto = build_variables.blareaob * build_variables.blnkoth
        fwbs_variables.volblkt = fwbs_variables.volblkti + fwbs_variables.volblkto

        #  Shield volume
        #  Uses fvolsi, fwbs_variables.fvolso as area coverage factors

        r1 = r1 + 0.5e0 * (build_variables.blnkith + build_variables.blnkoth)
        build_variables.sharea = physics_variables.sarea * r1 / physics_variables.rminor
        build_variables.shareaib = (
            0.5e0 * build_variables.sharea * fwbs_variables.fvolsi
        )
        build_variables.shareaob = (
            0.5e0 * build_variables.sharea * fwbs_variables.fvolso
        )

        volshldi = build_variables.shareaib * build_variables.shldith
        volshldo = build_variables.shareaob * build_variables.shldoth
        fwbs_variables.volshld = volshldi + volshldo

        fwbs_variables.whtshld = (
            fwbs_variables.volshld
            * fwbs_variables.denstl
            * (1.0e0 - fwbs_variables.vfshld)
        )

        #  Neutron power lost through holes in first wall (eventually absorbed by
        #  shield)

        fwbs_variables.pnucloss = physics_variables.pneutmw * fwbs_variables.fhole

        # The peaking factor, obtained as precalculated parameter
        fwbs_variables.wallpf = (
            stellarator_configuration.stella_config_neutron_peakfactor
        )

        #  Blanket neutronics calculations
        if fwbs_variables.blktmodel == 1:
            self.blanket_neutronics()

            if heat_transport_variables.ipowerflow == 1:
                fwbs_variables.pnucdiv = physics_variables.pneutmw * fwbs_variables.fdiv
                fwbs_variables.pnuchcd = physics_variables.pneutmw * fwbs_variables.fhcd
                fwbs_variables.pnucfw = (
                    physics_variables.pneutmw
                    - fwbs_variables.pnucdiv
                    - fwbs_variables.pnucloss
                    - fwbs_variables.pnuchcd
                )

                fwbs_variables.pradloss = (
                    physics_variables.pradmw * fwbs_variables.fhole
                )
                fwbs_variables.praddiv = physics_variables.pradmw * fwbs_variables.fdiv
                fwbs_variables.pradhcd = physics_variables.pradmw * fwbs_variables.fhcd
                fwbs_variables.pradfw = (
                    physics_variables.pradmw
                    - fwbs_variables.praddiv
                    - fwbs_variables.pradloss
                    - fwbs_variables.pradhcd
                )

                heat_transport_variables.htpmw_fw = heat_transport_variables.fpumpfw * (
                    fwbs_variables.pnucfw
                    + fwbs_variables.pradfw
                    + current_drive_variables.porbitlossmw
                )
                heat_transport_variables.htpmw_blkt = (
                    heat_transport_variables.fpumpblkt * fwbs_variables.pnucblkt
                )
                heat_transport_variables.htpmw_shld = (
                    heat_transport_variables.fpumpshld * fwbs_variables.pnucshld
                )
                heat_transport_variables.htpmw_div = (
                    heat_transport_variables.fpumpdiv
                    * (
                        physics_variables.pdivt
                        + fwbs_variables.pnucdiv
                        + fwbs_variables.praddiv
                    )
                )

                #  Void fraction in first wall / breeding zone,
                #  for use in fwbs_variables.fwmass and coolvol calculation below

                vffwi = (
                    1.0e0
                    - fwbs_variables.fblbe
                    - fwbs_variables.fblbreed
                    - fwbs_variables.fblss
                )
                vffwo = vffwi

        else:
            fwbs_variables.pnuc_cp = 0.0e0

            if heat_transport_variables.ipowerflow == 0:
                #  Energy-multiplied neutron power

                pneut2 = (
                    physics_variables.pneutmw
                    - fwbs_variables.pnucloss
                    - fwbs_variables.pnuc_cp
                ) * fwbs_variables.emult

                fwbs_variables.emultmw = pneut2 - (
                    physics_variables.pneutmw
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

                fwbs_variables.pnucblkt = pneut2 * (
                    1.0e0 - np.exp(-build_variables.blnkoth / decaybl)
                )

                #  Nuclear heating in the shield
                fwbs_variables.pnucshld = pneut2 - fwbs_variables.pnucblkt

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
                ) = fwbs_module.sctfcoil_nuclear_heating_iter90()

            else:  # heat_transport_variables.ipowerflow == 1
                #  Neutron power incident on divertor (MW)

                fwbs_variables.pnucdiv = physics_variables.pneutmw * fwbs_variables.fdiv

                #  Neutron power incident on HCD apparatus (MW)

                fwbs_variables.pnuchcd = physics_variables.pneutmw * fwbs_variables.fhcd

                #  Neutron power deposited in first wall, blanket and shield (MW)

                pnucfwbs = (
                    physics_variables.pneutmw
                    - fwbs_variables.pnucdiv
                    - fwbs_variables.pnucloss
                    - fwbs_variables.pnuc_cp
                    - fwbs_variables.pnuchcd
                )

                #  Split between inboard and outboard by first wall area fractions

                pnucfwbsi = pnucfwbs * build_variables.fwareaib / build_variables.fwarea
                pnucfwbso = pnucfwbs * build_variables.fwareaob / build_variables.fwarea

                #  Radiation power incident on divertor (MW)

                fwbs_variables.praddiv = physics_variables.pradmw * fwbs_variables.fdiv

                #  Radiation power incident on HCD apparatus (MW)

                fwbs_variables.pradhcd = physics_variables.pradmw * fwbs_variables.fhcd

                #  Radiation power lost through holes (eventually hits shield) (MW)

                fwbs_variables.pradloss = (
                    physics_variables.pradmw * fwbs_variables.fhole
                )

                #  Radiation power incident on first wall (MW)

                fwbs_variables.pradfw = (
                    physics_variables.pradmw
                    - fwbs_variables.praddiv
                    - fwbs_variables.pradloss
                    - fwbs_variables.pradhcd
                )

                #  Calculate the power deposited in the first wall, blanket and shield,
                #  and the required coolant pumping power

                #  If we have chosen pressurised water as the coolant, set the
                #  coolant outlet temperature as 20 deg C below the boiling point

                if fwbs_variables.coolwh == 2:
                    if fwbs_variables.irefprop:
                        fwbs_variables.outlet_temp = (
                            FluidProperties.of(
                                "Water",
                                pressure=fwbs_variables.coolp,
                                vapor_quality=0,
                            )
                            - 20
                        )
                    else:
                        fwbs_variables.outlet_temp = (
                            273.15
                            + 168.396
                            + 0.314653 / fwbs_variables.coolp
                            + -0.000728 / fwbs_variables.coolp**2
                            + 31.588979 * np.log(fwbs_variables.coolp)
                            + 11.473141 * fwbs_variables.coolp
                            + -0.575335 * fwbs_variables.coolp**2
                            + 0.013165 * fwbs_variables.coolp**3
                        ) - 20

                bfwi = 0.5e0 * build_variables.fwith
                bfwo = 0.5e0 * build_variables.fwoth

                vffwi = (
                    fwbs_variables.afwi * fwbs_variables.afwi / (bfwi * bfwi)
                )  # inboard FW coolant void fraction
                vffwo = (
                    fwbs_variables.afwo * fwbs_variables.afwo / (bfwo * bfwo)
                )  # outboard FW coolant void fraction

                #  First wall decay length (m) - improved calculation required

                decayfwi = fwbs_variables.declfw
                decayfwo = fwbs_variables.declfw

                #  Surface heat flux on first wall (MW) (sum = fwbs_variables.pradfw)

                psurffwi = (
                    fwbs_variables.pradfw
                    * build_variables.fwareaib
                    / build_variables.fwarea
                )
                psurffwo = (
                    fwbs_variables.pradfw
                    * build_variables.fwareaob
                    / build_variables.fwarea
                )

                #  Simple blanket model (fwbs_variables.primary_pumping = 0 or 1) is assumed for stellarators

                #  The power deposited in the first wall, breeder zone and shield is
                #  calculated according to their dimensions and materials assuming
                #  an exponential attenuation of nuclear heating with increasing
                #  radial distance.  The pumping power for the coolant is calculated
                #  as a fraction of the total thermal power deposited in the
                #  coolant.

                pnucfwi = pnucfwbsi * (1.0e0 - np.exp(-2.0e0 * bfwi / decayfwi))
                pnucfwo = pnucfwbso * (1.0e0 - np.exp(-2.0e0 * bfwo / decayfwo))

                #  Neutron power reaching blanket and shield (MW)

                pnucbsi = pnucfwbsi - pnucfwi
                pnucbso = pnucfwbso - pnucfwo

                #  Blanket decay length (m) - improved calculation required

                decaybzi = fwbs_variables.declblkt
                decaybzo = fwbs_variables.declblkt

                #  Neutron power deposited in breeder zone (MW)

                pnucbzi = pnucbsi * (
                    1.0e0 - np.exp(-build_variables.blnkith / decaybzi)
                )
                pnucbzo = pnucbso * (
                    1.0e0 - np.exp(-build_variables.blnkoth / decaybzo)
                )

                #  Calculate coolant pumping powers from input fraction.
                #  The pumping power is assumed to be a fraction, fpump, of the
                #  incident thermal power to each component so that
                #  htpmw_i = fpump_i*C, where C is the non-pumping thermal power
                #  deposited in the coolant

                #  First wall and Blanket pumping power (MW)

                if fwbs_variables.primary_pumping == 0:
                    #    Use input
                    pass
                elif fwbs_variables.primary_pumping == 1:
                    heat_transport_variables.htpmw_fw = (
                        heat_transport_variables.fpumpfw
                        * (
                            pnucfwi
                            + pnucfwo
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

                fwbs_variables.pnucfw = pnucfwi + pnucfwo

                #  Total nuclear heating of blanket (MW)

                fwbs_variables.pnucblkt = (pnucbzi + pnucbzo) * fwbs_variables.emult

                fwbs_variables.emultmw = fwbs_variables.emultmw + (
                    pnucbzi + pnucbzo
                ) * (fwbs_variables.emult - 1.0e0)

                #  Calculation of shield and divertor powers
                #  Shield and divertor powers and pumping powers are calculated using the same
                #  simplified method as the first wall and breeder zone when fwbs_variables.primary_pumping = 1.
                #  i.e. the pumping power is a fraction of the total thermal power deposited in the
                #  coolant.

                #  Neutron power reaching the shield (MW)
                #  The power lost from the fwbs_variables.fhole area fraction is assumed to be incident upon the shield

                pnucsi = (
                    pnucbsi
                    - pnucbzi
                    + (fwbs_variables.pnucloss + fwbs_variables.pradloss)
                    * build_variables.fwareaib
                    / build_variables.fwarea
                )
                pnucso = (
                    pnucbso
                    - pnucbzo
                    + (fwbs_variables.pnucloss + fwbs_variables.pradloss)
                    * build_variables.fwareaob
                    / build_variables.fwarea
                )

                #  Improved calculation of shield power decay lengths required

                decayshldi = fwbs_variables.declshld
                decayshldo = fwbs_variables.declshld

                #  Neutron power deposited in the shield (MW)

                pnucshldi = pnucsi * (
                    1.0e0 - np.exp(-build_variables.shldith / decayshldi)
                )
                pnucshldo = pnucso * (
                    1.0e0 - np.exp(-build_variables.shldoth / decayshldo)
                )

                fwbs_variables.pnucshld = pnucshldi + pnucshldo

                #  Calculate coolant pumping powers from input fraction.
                #  The pumping power is assumed to be a fraction, fpump, of the incident
                #  thermal power to each component so that,
                #     htpmw_i = fpump_i*C
                #  where C is the non-pumping thermal power deposited in the coolant

                if fwbs_variables.primary_pumping == 1:
                    #  Shield pumping power (MW)
                    heat_transport_variables.htpmw_shld = (
                        heat_transport_variables.fpumpshld * (pnucshldi + pnucshldo)
                    )

                    #  Divertor pumping power (MW)
                    heat_transport_variables.htpmw_div = (
                        heat_transport_variables.fpumpdiv
                        * (
                            physics_variables.pdivt
                            + fwbs_variables.pnucdiv
                            + fwbs_variables.praddiv
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
                    fwbs_variables.volblkt * fwbs_variables.fbllipb * 9400.0e0
                )
                fwbs_variables.whtblli = (
                    fwbs_variables.volblkt * fwbs_variables.fblli * 534.0e0
                )
                fwbs_variables.whtblkt = (
                    fwbs_variables.wtbllipb + fwbs_variables.whtblli
                )
            else:  # solid breeder (HCPB); always for ipowerflow=0
                fwbs_variables.wtblli2o = (
                    fwbs_variables.volblkt * fwbs_variables.fblli2o * 2010.0e0
                )
                fwbs_variables.whtblbe = (
                    fwbs_variables.volblkt * fwbs_variables.fblbe * 1850.0e0
                )
                fwbs_variables.whtblkt = (
                    fwbs_variables.wtblli2o + fwbs_variables.whtblbe
                )

            fwbs_variables.whtblss = (
                fwbs_variables.volblkt * fwbs_variables.denstl * fwbs_variables.fblss
            )
            fwbs_variables.whtblvd = (
                fwbs_variables.volblkt * 5870.0e0 * fwbs_variables.fblvd
            )

            fwbs_variables.whtblkt = (
                fwbs_variables.whtblkt + fwbs_variables.whtblss + fwbs_variables.whtblvd
            )

        else:  # volume fractions proportional to sub-assembly thicknesses
            fwbs_variables.whtblss = fwbs_variables.denstl * (
                fwbs_variables.volblkti
                / build_variables.blnkith
                * (
                    build_variables.blbuith * fwbs_variables.fblss
                    + build_variables.blbmith * (1.0e0 - fwbs_variables.fblhebmi)
                    + build_variables.blbpith * (1.0e0 - fwbs_variables.fblhebpi)
                )
                + fwbs_variables.volblkto
                / build_variables.blnkoth
                * (
                    build_variables.blbuoth * fwbs_variables.fblss
                    + build_variables.blbmoth * (1.0e0 - fwbs_variables.fblhebmo)
                    + build_variables.blbpoth * (1.0e0 - fwbs_variables.fblhebpo)
                )
            )
            fwbs_variables.whtblbe = (
                1850.0e0
                * fwbs_variables.fblbe
                * (
                    (
                        fwbs_variables.volblkti
                        * build_variables.blbuith
                        / build_variables.blnkith
                    )
                    + (
                        fwbs_variables.volblkto
                        * build_variables.blbuoth
                        / build_variables.blnkoth
                    )
                )
            )
            fwbs_variables.whtblbreed = (
                fwbs_variables.densbreed
                * fwbs_variables.fblbreed
                * (
                    (
                        fwbs_variables.volblkti
                        * build_variables.blbuith
                        / build_variables.blnkith
                    )
                    + (
                        fwbs_variables.volblkto
                        * build_variables.blbuoth
                        / build_variables.blnkoth
                    )
                )
            )
            fwbs_variables.whtblkt = (
                fwbs_variables.whtblss
                + fwbs_variables.whtblbe
                + fwbs_variables.whtblbreed
            )

            fwbs_variables.vfblkt = (
                fwbs_variables.volblkti
                / fwbs_variables.volblkt
                * (  # inboard portion
                    (build_variables.blbuith / build_variables.blnkith)
                    * (
                        1.0e0
                        - fwbs_variables.fblbe
                        - fwbs_variables.fblbreed
                        - fwbs_variables.fblss
                    )
                    + (build_variables.blbmith / build_variables.blnkith)
                    * fwbs_variables.fblhebmi
                    + (build_variables.blbpith / build_variables.blnkith)
                    * fwbs_variables.fblhebpi
                )
            )
            fwbs_variables.vfblkt = (
                fwbs_variables.vfblkt
                + fwbs_variables.volblkto
                / fwbs_variables.volblkt
                * (  # outboard portion
                    (build_variables.blbuoth / build_variables.blnkoth)
                    * (
                        1.0e0
                        - fwbs_variables.fblbe
                        - fwbs_variables.fblbreed
                        - fwbs_variables.fblss
                    )
                    + (build_variables.blbmoth / build_variables.blnkoth)
                    * fwbs_variables.fblhebmo
                    + (build_variables.blbpoth / build_variables.blnkoth)
                    * fwbs_variables.fblhebpo
                )
            )

        #  When fwbs_variables.blktmodel > 0, although the blanket is by definition helium-cooled
        #  in this case, the shield etc. are assumed to be water-cooled, and since
        #  water is heavier the calculation for fwbs_variables.coolmass is better done with
        #  coolwh=2 if fwbs_variables.blktmodel > 0; thus we can ignore the helium coolant mass
        #  in the blanket.

        if fwbs_variables.blktmodel == 0:
            coolvol = coolvol + fwbs_variables.volblkt * fwbs_variables.vfblkt

        coolvol = coolvol + fwbs_variables.volshld * fwbs_variables.vfshld

        #  Penetration shield (set = internal shield)

        fwbs_variables.wpenshld = fwbs_variables.whtshld

        if heat_transport_variables.ipowerflow == 0:
            #  First wall mass
            #  (first wall area is calculated else:where)

            fwbs_variables.fwmass = (
                build_variables.fwarea
                * (build_variables.fwith + build_variables.fwoth)
                / 2.0e0
                * fwbs_variables.denstl
                * (1.0e0 - fwbs_variables.fwclfr)
            )

            #  Surface areas adjacent to plasma

            coolvol = (
                coolvol
                + build_variables.fwarea
                * (build_variables.fwith + build_variables.fwoth)
                / 2.0e0
                * fwbs_variables.fwclfr
            )

        else:
            fwbs_variables.fwmass = fwbs_variables.denstl * (
                build_variables.fwareaib * build_variables.fwith * (1.0e0 - vffwi)
                + build_variables.fwareaob * build_variables.fwoth * (1.0e0 - vffwo)
            )
            coolvol = (
                coolvol
                + build_variables.fwareaib * build_variables.fwith * vffwi
                + build_variables.fwareaob * build_variables.fwoth * vffwo
            )

            #  Average first wall coolant fraction, only used by old routines
            #  in fispact.f90, safety.f90

            fwbs_variables.fwclfr = (
                build_variables.fwareaib * build_variables.fwith * vffwi
                + build_variables.fwareaob * build_variables.fwoth * vffwo
            ) / (
                build_variables.fwarea
                * 0.5e0
                * (build_variables.fwith + build_variables.fwoth)
            )

        #  Mass of coolant = volume * density at typical coolant
        #  temperatures and pressures
        #  N.B. for fwbs_variables.blktmodel > 0, mass of *water* coolant in the non-blanket
        #  structures is used (see comment above)

        if (fwbs_variables.blktmodel > 0) or (
            fwbs_variables.coolwh == 2
        ):  # pressurised water coolant
            fwbs_variables.coolmass = coolvol * 806.719e0
        else:  # gaseous helium coolant
            fwbs_variables.coolmass = coolvol * 1.517e0

        #  Assume external cryostat is a torus with circular cross-section,
        #  centred on plasma major radius.
        #  N.B. No check made to see if coils etc. lie wholly within cryostat...

        #  External cryostat outboard major radius (m)

        fwbs_variables.rdewex = (
            build_variables.r_tf_outboard_mid
            + 0.5e0 * build_variables.tfthko
            + fwbs_variables.rpf2dewar
        )
        adewex = fwbs_variables.rdewex - physics_variables.rmajor

        #  External cryostat volume

        fwbs_variables.vdewex = (
            4.0e0
            * (np.pi**2)
            * physics_variables.rmajor
            * adewex
            * build_variables.ddwex
        )

        #  Internal vacuum vessel volume
        #  fwbs_variables.fvoldw accounts for ports, support, etc. additions

        r1 = physics_variables.rminor + 0.5e0 * (
            build_variables.scrapli
            + build_variables.fwith
            + build_variables.blnkith
            + build_variables.shldith
            + build_variables.scraplo
            + build_variables.fwoth
            + build_variables.blnkoth
            + build_variables.shldoth
        )
        fwbs_variables.vdewin = (
            (build_variables.d_vv_in + build_variables.d_vv_out)
            / 2.0e0
            * physics_variables.sarea
            * r1
            / physics_variables.rminor
            * fwbs_variables.fvoldw
        )

        #  Vacuum vessel mass

        fwbs_variables.vvmass = fwbs_variables.vdewin * fwbs_variables.denstl

        #  Sum of internal vacuum vessel and external cryostat masses

        fwbs_variables.dewmkg = (
            fwbs_variables.vdewin + fwbs_variables.vdewex
        ) * fwbs_variables.denstl

        if output:
            #  Output section

            po.oheadr(self.outfile, "First Wall / Blanket / Shield")
            po.ovarre(
                self.outfile,
                "Average neutron wall load (MW/m2)",
                "(wallmw)",
                physics_variables.wallmw,
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
                "(fwlife)",
                fwbs_variables.fwlife,
            )

            po.ovarre(
                self.outfile,
                "Inboard shield thickness (m)",
                "(shldith)",
                build_variables.shldith,
            )
            po.ovarre(
                self.outfile,
                "Outboard shield thickness (m)",
                "(shldoth)",
                build_variables.shldoth,
            )
            po.ovarre(
                self.outfile,
                "Top shield thickness (m)",
                "(shldtth)",
                build_variables.shldtth,
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
                "(blnkith)",
                build_variables.blnkith,
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
                "(blnkoth)",
                build_variables.blnkoth,
            )
            po.ovarre(
                self.outfile,
                "Top blanket thickness (m)",
                "(blnktth)",
                build_variables.blnktth,
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
                    "(pnucblkt)",
                    fwbs_variables.pnucblkt,
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
                    "(pnucblkt)",
                    fwbs_variables.pnucblkt,
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
                    "(li6enrich)",
                    fwbs_variables.li6enrich,
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
                    "(bktlife)",
                    fwbs_variables.bktlife,
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
                    "(secondary_cycle)",
                    fwbs_variables.secondary_cycle,
                )
                if fwbs_variables.secondary_cycle == 0:
                    po.ocmmnt(self.outfile, "   (Simple calculation)")

            po.osubhd(self.outfile, "Blanket / shield volumes and weights :")

            #     if (fwbs_variables.blktmodel == 0) :
            #         if ((fwbs_variables.blkttype == 1)or(fwbs_variables.blkttype == 2)) :
            #             po.write(self.outfile,601) volblkti, volblkto, volblkt,                whtblkt, vfblkt, fbllipb, wtbllipb, fblli, whtblli,                fblss, whtblss, fblvd, whtblvd, volshldi, volshldo,                volshld, whtshld, vfshld, fwbs_variables.wpenshld
            #         else:  #  (also if ipowerflow=0)
            #             po.write(self.outfile,600) volblkti, volblkto, volblkt,                whtblkt, vfblkt, fblbe, whtblbe, fblli2o, wtblli2o,                fblss, whtblss, fblvd, whtblvd, volshldi, volshldo,                volshld, whtshld, vfshld, fwbs_variables.wpenshld

            #     else:
            #         po.write(self.outfile,602) volblkti, volblkto, volblkt, whtblkt, vfblkt,             (fwbs_variables.volblkti/fwbs_variables.volblkt * build_variables.blbuith/build_variables.blnkith +             fwbs_variables.volblkto/fwbs_variables.volblkt * build_variables.blbuoth/build_variables.blnkoth) * fblbe, whtblbe,             (fwbs_variables.volblkti/fwbs_variables.volblkt * build_variables.blbuith/build_variables.blnkith +             fwbs_variables.volblkto/fwbs_variables.volblkt * build_variables.blbuoth/build_variables.blnkoth) * fblbreed, whtblbreed,             fwbs_variables.volblkti/fwbs_variables.volblkt/build_variables.blnkith * (build_variables.blbuith * fwbs_variables.fblss             + build_variables.blbmith * (1.0e0-fwbs_variables.fblhebmi) + build_variables.blbpith * (1.0e0-fwbs_variables.fblhebpi)) +             fwbs_variables.volblkto/fwbs_variables.volblkt/build_variables.blnkoth * (build_variables.blbuoth * fwbs_variables.fblss             + build_variables.blbmoth * (1.0e0-fwbs_variables.fblhebmo) + build_variables.blbpoth * (1.0e0-fwbs_variables.fblhebpo)),             whtblss,             volshldi, volshldo, volshld, whtshld, vfshld, fwbs_variables.wpenshld

            # 600 format(          t32,'volume (m3)',t45,'vol fraction',t62,'weight (kg)'/          t32,'-----------',t45,'------------',t62,'-----------'/          '    Inboard blanket' ,t32,1pe10.3,/          '    Outboard blanket' ,t32,1pe10.3,/          '    Total blanket' ,t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '       Blanket Be   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket Li2O ',t45,1pe10.3,t62,1pe10.3/          '       Blanket ss   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket Vd   ',t45,1pe10.3,t62,1pe10.3/          '    Inboard shield'  ,t32,1pe10.3,/          '    Outboard shield'  ,t32,1pe10.3,/          '    Primary shield',t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '    Penetration shield'        ,t62,1pe10.3)

            # 601 format(          t32,'volume (m3)',t45,'vol fraction',t62,'weight (kg)'/          t32,'-----------',t45,'------------',t62,'-----------'/          '    Inboard blanket' ,t32,1pe10.3,/          '    Outboard blanket' ,t32,1pe10.3,/          '    Total blanket' ,t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '       Blanket LiPb ',t45,1pe10.3,t62,1pe10.3/          '       Blanket Li   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket ss   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket Vd   ',t45,1pe10.3,t62,1pe10.3/          '    Inboard shield'  ,t32,1pe10.3,/          '    Outboard shield'  ,t32,1pe10.3,/          '    Primary shield',t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '    Penetration shield'        ,t62,1pe10.3)

            # 602 format(          t32,'volume (m3)',t45,'vol fraction',t62,'weight (kg)'/          t32,'-----------',t45,'------------',t62,'-----------'/          '    Inboard blanket' ,t32,1pe10.3,/          '    Outboard blanket' ,t32,1pe10.3,/          '    Total blanket' ,t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '       Blanket Be   ',t45,1pe10.3,t62,1pe10.3/          '       Blanket breeder',t45,1pe10.3,t62,1pe10.3/          '       Blanket steel',t45,1pe10.3,t62,1pe10.3/          '    Inboard shield'  ,t32,1pe10.3,/          '    Outboard shield'  ,t32,1pe10.3,/          '    Primary shield',t32,1pe10.3,t62,1pe10.3/          '       Void fraction' ,t45,1pe10.3,/          '    Penetration shield'        ,t62,1pe10.3)

            po.osubhd(self.outfile, "Other volumes, masses and areas :")
            po.ovarre(
                self.outfile, "First wall area (m2)", "(fwarea)", build_variables.fwarea
            )
            po.ovarre(
                self.outfile, "First wall mass (kg)", "(fwmass)", fwbs_variables.fwmass
            )
            po.ovarre(
                self.outfile,
                "External cryostat inner radius (m)",
                "",
                fwbs_variables.rdewex - 2.0e0 * adewex,
            )
            po.ovarre(
                self.outfile,
                "External cryostat outer radius (m)",
                "(rdewex)",
                fwbs_variables.rdewex,
            )
            po.ovarre(
                self.outfile, "External cryostat minor radius (m)", "(adewex)", adewex
            )
            po.ovarre(
                self.outfile,
                "External cryostat shell volume (m3)",
                "(vdewex)",
                fwbs_variables.vdewex,
            )
            po.ovarre(
                self.outfile,
                "External cryostat mass (kg)",
                "",
                fwbs_variables.dewmkg - fwbs_variables.vvmass,
            )
            po.ovarre(
                self.outfile,
                "Internal vacuum vessel shell volume (m3)",
                "(vdewin)",
                fwbs_variables.vdewin,
            )
            po.ovarre(
                self.outfile,
                "Vacuum vessel mass (kg)",
                "(vvmass)",
                fwbs_variables.vvmass,
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
                tfcoil_variables.n_tf,
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
            tfcoil_variables.n_tf,
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
        tfcoil_variables.n_tf_turn = awptf / (
            tfcoil_variables.t_turn_tf**2
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

        build_variables.tfcth = (
            tfcoil_variables.thkcas
            + tfcoil_variables.dr_tf_wp
            + tfcoil_variables.casthi
            + 2.0e0 * tfcoil_variables.tinstf
        )  # [m] Thickness of inboard leg in radial direction
        build_variables.tfthko = (
            tfcoil_variables.thkcas
            + tfcoil_variables.dr_tf_wp
            + tfcoil_variables.casthi
            + 2.0e0 * tfcoil_variables.tinstf
        )  # [m] Thickness of outboard leg in radial direction (same as inboard)
        tfcoil_variables.arealeg = (
            build_variables.tfcth * tfcoil_variables.tftort
        )  # [m^2] overall coil cross-sectional area (assuming inboard and
        #       outboard leg are the same)
        tfcoil_variables.acasetf = (
            build_variables.tfcth * tfcoil_variables.tftort
        ) - awpc  # [m^2] Cross-sectional area of surrounding case

        tfcoil_variables.tfocrn = (
            0.5e0 * tfcoil_variables.tftort
        )  # [m] Half-width of side of coil nearest torus centreline
        tfcoil_variables.tficrn = (
            0.5e0 * tfcoil_variables.tftort
        )  # [m] Half-width of side of coil nearest plasma

        # [m^2] Total surface area of coil side facing plasma: inboard region
        tfcoil_variables.tfsai = (
            tfcoil_variables.n_tf
            * tfcoil_variables.tftort
            * 0.5e0
            * tfcoil_variables.tfleng
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
            tfcoil_variables.n_tf * tfcoil_variables.arealeg
        )  # [m^2] Total area of all coil legs (midplane)
        tfcoil_variables.ritfc = (
            tfcoil_variables.n_tf * coilcurrent * 1.0e6
        )  # [A] Total current in ALL coils
        tfcoil_variables.oacdcp = (
            tfcoil_variables.ritfc / tfcoil_variables.tfareain
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
            * (tfcoil_variables.ritfc / tfcoil_variables.n_tf) ** 2
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

        tfborev = 2.0e0 * build_variables.hmax  # [m] estimated vertical coil bore

        tfcoil_variables.tfleng = (
            stellarator_configuration.stella_config_coillength
            * (r_coil_minor / stellarator_configuration.stella_config_coil_rminor)
            / tfcoil_variables.n_tf
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
        # This is only correct if the winding pack is 'thin' (tfleng>>sqrt(tfcoil_variables.acasetf)).
        tfcoil_variables.whtcas = (
            tfcoil_variables.tfleng * tfcoil_variables.acasetf * tfcoil_variables.dcase
        )
        # Mass of ground-wall insulation [kg]
        # (assumed to be same density/material as conduit insulation)
        tfcoil_variables.whtgw = (
            tfcoil_variables.tfleng * (awpc - awptf) * tfcoil_variables.dcondins
        )
        # [kg] mass of Superconductor
        tfcoil_variables.whtconsc = (
            tfcoil_variables.tfleng
            * tfcoil_variables.n_tf_turn
            * tfcoil_variables.acstf
            * (1.0e0 - tfcoil_variables.vftf)
            * (1.0e0 - tfcoil_variables.fcutfsu)
            - tfcoil_variables.tfleng * tfcoil_variables.awphec
        ) * tfcoil_variables.dcond[
            tfcoil_variables.i_tf_sc_mat - 1
        ]  # awphec is 0 for a stellarator. but keep this term for now.
        # [kg] mass of Copper in conductor
        tfcoil_variables.whtconcu = (
            tfcoil_variables.tfleng
            * tfcoil_variables.n_tf_turn
            * tfcoil_variables.acstf
            * (1.0e0 - tfcoil_variables.vftf)
            * tfcoil_variables.fcutfsu
            - tfcoil_variables.tfleng * tfcoil_variables.awphec
        ) * constants.dcopper
        # [kg] mass of Steel conduit (sheath)
        tfcoil_variables.whtconsh = (
            tfcoil_variables.tfleng
            * tfcoil_variables.n_tf_turn
            * tfcoil_variables.acndttf
            * fwbs_variables.denstl
        )
        # if (i_tf_sc_mat==6)   tfcoil_variables.whtconsh = fcondsteel * awptf *tfcoil_variables.tfleng* fwbs_variables.denstl
        # Conduit insulation mass [kg]
        # (tfcoil_variables.aiwp already contains tfcoil_variables.n_tf_turn)
        tfcoil_variables.whtconin = (
            tfcoil_variables.tfleng * tfcoil_variables.aiwp * tfcoil_variables.dcondins
        )
        # [kg] Total conductor mass
        tfcoil_variables.whtcon = (
            tfcoil_variables.whtconsc
            + tfcoil_variables.whtconcu
            + tfcoil_variables.whtconsh
            + tfcoil_variables.whtconin
        )
        # [kg] Total coil mass
        tfcoil_variables.whttf = (
            tfcoil_variables.whtcas + tfcoil_variables.whtcon + tfcoil_variables.whtgw
        ) * tfcoil_variables.n_tf
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
            - build_variables.scrapli
            - build_variables.fwith
            - build_variables.blnkith
            - build_variables.vvblgap
            - build_variables.shldith
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
        f_vv_actual = (
            2.54e6
            * (3e0 * 1.3e0 * 50e0 * 0.92e0**2e0)
            / (1e0 * 5.2e0 * 0.014e0)
            * (
                physics_variables.bt
                * tfcoil_variables.ritfc
                * physics_variables.rminor**2
                / (
                    (build_variables.d_vv_in + build_variables.d_vv_out)
                    / 2
                    * tfcoil_variables.tdmptf
                    * radvv
                )
            )
            ** (-1)
        )

        # the conductor fraction is meant of the cable space#
        # This is the old routine which is being replaced for now by the new one below
        #    protect(aio,  tfes,               acs,       aturn,   tdump,  fcond,  fcu,   tba,  tmax   ,ajwpro, vd)
        # call protect(cpttf,estotftgj/tfcoil_variables.n_tf*1.0e9,acstf,   tfcoil_variables.t_turn_tf**2   ,tdmptf,1-vftf,fcutfsu,tftmp,tmaxpro,jwdgpro2,vd)

        vd = self.u_max_protect_v(
            tfcoil_variables.estotftgj / tfcoil_variables.n_tf * 1.0e9,
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
            / tfcoil_variables.n_tf
            / tfcoil_variables.tfleng
        )
        centering_force_min_mn = (
            stellarator_configuration.stella_config_centering_force_min_mn
            * st.f_i
            / st.f_n
            * tfcoil_variables.bmaxtf
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_coillength
            / tfcoil_variables.n_tf
            / tfcoil_variables.tfleng
        )
        centering_force_avg_mn = (
            stellarator_configuration.stella_config_centering_force_avg_mn
            * st.f_i
            / st.f_n
            * tfcoil_variables.bmaxtf
            / stellarator_configuration.stella_config_wp_bmax
            * stellarator_configuration.stella_config_coillength
            / tfcoil_variables.n_tf
            / tfcoil_variables.tfleng
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
        q_cu_array_sA2m4 = [
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
        q_he_array_sA2m4 = [
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

        q_he = maths_library.find_y_nonuniform_x(temp, temp_k, q_he_array_sA2m4, 13)
        q_cu = maths_library.find_y_nonuniform_x(temp, temp_k, q_cu_array_sA2m4, 13)

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
    ):
        strain = -0.005  # for now a small value
        # fhe = vftf  # this is helium fraction in the superconductor (set it to the fixed global variable here)

        fcu = fcutfsu  # fcutfsu is a global variable. Is the copper fraction
        # of a cable conductor.

        if i_tf_sc_mat == 1:  # ITER Nb3Sn critical surface parameterization
            bc20m = 32.97  # these are values taken from sctfcoil.f90
            tc0m = 16.06

            #  jcritsc returned by itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper
            if bmax > bc20m:
                jcritsc = 1.0e-9  # Set to a small nonzero value
            else:
                (
                    jcritsc,
                    bcrit,
                    tcrit,
                ) = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)

            jcritstr = jcritsc * (1.0 - fcu)

            # This is needed right now. Can we change it later?
            jcritsc = max(1.0e-9, jcritsc)
            jcritstr = max(1.0e-9, jcritstr)

        elif (
            i_tf_sc_mat == 2
        ):  # Bi-2212 high temperature superconductor parameterization
            #  Current density in a strand of Bi-2212 conductor
            #  N.B. jcrit returned by bi2212 is the critical current density
            #  in the strand, not just the superconducting portion.
            #  The parameterization for jcritstr assumes a particular strand
            #  composition that does not require a user-defined copper fraction,
            #  so this is irrelevant in this model

            # jstrand = jwp / (1 - fhe)
            jstrand = 0  # as far as I can tell this will always be 0
            # because jwp was never set in fortran (so 0)

            jcritstr, tmarg = superconductors.bi2212(
                bmax, jstrand, thelium, fhts
            )  # bi2212 outputs jcritstr
            jcritsc = jcritstr / (1 - fcu)
            tcrit = thelium + tmarg
        elif i_tf_sc_mat == 3:  # NbTi data
            bc20m = 15.0
            tc0m = 9.3
            c0 = 1.0

            if bmax > bc20m:
                jcritsc = 1.0e-9  # Set to a small nonzero value
            else:
                jcritsc, tcrit = superconductors.jcrit_nbti(
                    thelium,
                    bmax,
                    c0,
                    bc20m,
                    tc0m,
                )
                # I dont need tcrit here so dont use it.

            jcritstr = jcritsc * (1 - fcu)

            # This is needed right now. Can we change it later?
            jcritsc = max(1.0e-9, jcritsc)
            jcritstr = max(1.0e-9, jcritstr)
        elif i_tf_sc_mat == 4:  # As (1), but user-defined parameters
            bc20m = bcritsc
            tc0m = tcritsc
            jcritsc, bcrit, tcrit = superconductors.itersc(
                thelium, bmax, strain, bc20m, tc0m
            )
            jcritstr = jcritsc * (1 - fcu)
        elif i_tf_sc_mat == 5:  # WST Nb3Sn parameterisation
            bc20m = 32.97
            tc0m = 16.06

            #  jcritsc returned by itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper

            jcritsc, bcrit, tcrit = superconductors.wstsc(
                thelium,
                bmax,
                strain,
                bc20m,
                tc0m,
            )
            jcritstr = jcritsc * (1 - fcu)
        elif (
            i_tf_sc_mat == 6
        ):  # ! "REBCO" 2nd generation HTS superconductor in CrCo strand
            jcritsc, validity = superconductors.jcrit_rebco(thelium, bmax, 0)
            jcritsc = max(1.0e-9, jcritsc)
            jcritstr = jcritsc * (1 - fcu)

        elif i_tf_sc_mat == 7:  # Durham Ginzburg-Landau Nb-Ti parameterisation
            bc20m = b_crit_upper_nbti
            tc0m = t_crit_nbti
            jcritsc, bcrit, tcrit = superconductors.gl_nbti(
                thelium, bmax, strain, bc20m, tc0m
            )
            jcritstr = jcritsc * (1 - fcu)
        elif i_tf_sc_mat == 8:
            bc20m = 429
            tc0m = 185
            jcritsc, bcrit, tcrit = superconductors.gl_rebco(
                thelium, bmax, strain, bc20m, tc0m
            )
            # A0 calculated for tape cross section already
            jcritstr = jcritsc * (1 - fcu)
        else:
            error_handling.idiags[0] = i_tf_sc_mat
            error_handling.report_error(156)

        return jcritsc * 1e-6

    def bmax_from_awp(self, wp_width_radial, current, n_tf, r_coil_major, r_coil_minor):
        """Returns a fitted function for bmax for stellarators

        author: J Lion, IPP Greifswald
        Returns a fitted function for bmax in dependece
        of the winding pack. The stellarator type config
        is taken from the parent scope.
        """

        return (
            2e-1
            * current
            * n_tf
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

        for i in range(100):
            #  Find difference in y values at x

            y01 = maths_library.find_y_nonuniform_x(x, x1, y1, n1)
            y02 = maths_library.find_y_nonuniform_x(x, x2, y2, n2)
            y = y01 - y02

            if abs(y) < epsy:
                break

            #  Find difference in y values at x+dx

            y01 = maths_library.find_y_nonuniform_x(x + dx, x1, y1, n1)
            y02 = maths_library.find_y_nonuniform_x(x + dx, x2, y2, n2)
            yright = y01 - y02

            #  Find difference in y values at x-dx

            y01 = maths_library.find_y_nonuniform_x(x - dx, x1, y1, n1)
            y02 = maths_library.find_y_nonuniform_x(x - dx, x2, y2, n2)
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
            "Maximum reachable ECRH temperature (pseudo) (keV)",
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
        te0_available : input real : Reachable peak electron temperature, reached by ECRH (keV)
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
            copy(physics_variables.powerht), 0.00001e0
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
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
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
            self.outfile, "Number of modular coils", "(n_tf)", tfcoil_variables.n_tf
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
            "(tfarea/n_tf)",
            tfcoil_variables.tfareain / tfcoil_variables.n_tf,
        )
        po.ovarre(
            self.outfile,
            "Total inboard leg radial thickness (m)",
            "(tfcth)",
            build_variables.tfcth,
        )
        po.ovarre(
            self.outfile,
            "Total outboard leg radial thickness (m)",
            "(tfthko)",
            build_variables.tfthko,
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
            "(tfleng)",
            tfcoil_variables.tfleng,
        )
        po.ovarre(
            self.outfile,
            "Total current (MA)",
            "(ritfc)",
            1.0e-6 * tfcoil_variables.ritfc,
        )
        po.ovarre(
            self.outfile,
            "Current per coil(MA)",
            "(ritfc/n_tf)",
            1.0e-6 * tfcoil_variables.ritfc / tfcoil_variables.n_tf,
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
            self.outfile, "Total mass of coils (kg)", "(whttf)", tfcoil_variables.whttf
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
            "Clear horizontal bore (m)",
            "(tf_total_h_width)",
            tf_total_h_width,
        )
        po.ovarre(self.outfile, "Clear vertical bore (m)", "(tfborev)", tfborev)

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
            "(ritfc/acond)",
            1.0e-6
            * tfcoil_variables.ritfc
            / tfcoil_variables.n_tf
            / tfcoil_variables.acond,
        )
        po.ovarre(
            self.outfile,
            "Current density in SC area (A/m2)",
            "(ritfc/acond/f_scu)",
            1.0e-6 * tfcoil_variables.ritfc / tfcoil_variables.n_tf / ap / f_scu,
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
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        # ###############################################
        #  Calculate plasma composition
        # Issue #261 Remove old radiation model

        physics_module.plasma_composition()

        # Calculate density and temperature profile quantities
        profiles_module.plasma_profiles()

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
            physics_variables.betaft
            + physics_variables.betanb
            + 2.0e3
            * constants.rmu0
            * constants.echarge
            * (
                physics_variables.dene * physics_variables.ten
                + physics_variables.dnitot * physics_variables.tin
            )
            / physics_variables.btot**2
        )
        physics_module.total_plasma_internal_energy = (
            1.5e0
            * physics_variables.beta
            * physics_variables.btot
            * physics_variables.btot
            / (2.0e0 * constants.rmu0)
            * physics_variables.vol
        )

        physics_module.rho_star = np.sqrt(
            2.0e0
            * constants.mproton
            * physics_variables.aion
            * physics_module.total_plasma_internal_energy
            / (3.0e0 * physics_variables.vol * physics_variables.dnla)
        ) / (
            constants.echarge
            * physics_variables.bt
            * physics_variables.eps
            * physics_variables.rmajor
        )

        physics_variables.q95 = physics_variables.q

        #  Calculate poloidal field using rotation transform
        physics_variables.bp = (
            physics_variables.rminor
            * physics_variables.bt
            / physics_variables.rmajor
            * stellarator_variables.iotabar
        )

        #  Poloidal physics_variables.beta

        # betap = physics_variables.beta * ( physics_variables.btot/physics_variables.bp )**2 # Dont need this I think.

        #  Perform auxiliary power calculations

        self.stheat(False)

        #  Calculate fusion power

        (
            physics_variables.palppv,
            physics_variables.pchargepv,
            physics_variables.pneutpv,
            sigvdt,
            physics_variables.fusionrate,
            physics_variables.alpharate,
            physics_variables.protonrate,
            pdtpv,
            pdhe3pv,
            pddpv,
        ) = physics_functions_module.palph(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.deni,
            physics_variables.fdeut,
            physics_variables.fhe3,
            physics_variables.ftrit,
            physics_variables.ti,
        )

        physics_variables.pdt = pdtpv * physics_variables.vol
        physics_variables.pdhe3 = pdhe3pv * physics_variables.vol
        physics_variables.pdd = pddpv * physics_variables.vol

        #  Calculate neutral beam slowing down effects
        #  If ignited, then ignore beam fusion effects

        if (current_drive_variables.pnbeam != 0.0e0) and (
            physics_variables.ignite == 0
        ):
            (
                physics_variables.betanb,
                physics_variables.dnbeam2,
                physics_variables.palpnb,
            ) = physics_functions_module.beamfus(
                physics_variables.beamfus0,
                physics_variables.betbm0,
                physics_variables.bp,
                physics_variables.bt,
                current_drive_variables.cnbeam,
                physics_variables.dene,
                physics_variables.deni,
                physics_variables.dlamie,
                physics_variables.ealphadt,
                current_drive_variables.enbeam,
                physics_variables.fdeut,
                physics_variables.ftrit,
                current_drive_variables.ftritbm,
                sigvdt,
                physics_variables.ten,
                physics_variables.tin,
                physics_variables.vol,
                physics_variables.zeffai,
            )
            physics_variables.fusionrate = (
                physics_variables.fusionrate
                + 1.0e6
                * physics_variables.palpnb
                / (1.0e3 * physics_variables.ealphadt * constants.echarge)
                / physics_variables.vol
            )
            physics_variables.alpharate = (
                physics_variables.alpharate
                + 1.0e6
                * physics_variables.palpnb
                / (1.0e3 * physics_variables.ealphadt * constants.echarge)
                / physics_variables.vol
            )

        physics_variables.pdt = physics_variables.pdt + 5.0e0 * physics_variables.palpnb

        (
            physics_variables.palpmw,
            physics_variables.pneutmw,
            physics_variables.pchargemw,
            physics_variables.betaft,
            physics_variables.palpipv,
            physics_variables.palpepv,
            physics_variables.pfuscmw,
            physics_variables.powfmw,
        ) = physics_functions_module.palph2(
            physics_variables.bt,
            physics_variables.bp,
            physics_variables.dene,
            physics_variables.deni,
            physics_variables.dnitot,
            physics_variables.falpe,
            physics_variables.falpi,
            physics_variables.palpnb,
            physics_variables.ifalphap,
            physics_variables.pchargepv,
            physics_variables.pneutpv,
            physics_variables.ten,
            physics_variables.tin,
            physics_variables.vol,
            physics_variables.palppv,
        )

        #  Neutron wall load

        if physics_variables.iwalld == 1:
            physics_variables.wallmw = (
                physics_variables.ffwal
                * physics_variables.pneutmw
                / physics_variables.sarea
            )
        else:
            if heat_transport_variables.ipowerflow == 0:
                physics_variables.wallmw = (
                    (1.0e0 - fwbs_variables.fhole)
                    * physics_variables.pneutmw
                    / build_variables.fwarea
                )
            else:
                physics_variables.wallmw = (
                    (
                        1.0e0
                        - fwbs_variables.fhole
                        - fwbs_variables.fhcd
                        - fwbs_variables.fdiv
                    )
                    * physics_variables.pneutmw
                    / build_variables.fwarea
                )

        #  Calculate ion/electron equilibration power

        physics_variables.piepv = physics_module.rether(
            physics_variables.alphan,
            physics_variables.alphat,
            physics_variables.dene,
            physics_variables.dlamie,
            physics_variables.te,
            physics_variables.ti,
            physics_variables.zeffai,
        )

        #  Calculate radiation power
        radpwr_data = physics_funcs.radpwr(self.plasma_profile)
        physics_variables.pbrempv = radpwr_data.pbrempv
        physics_variables.plinepv = radpwr_data.plinepv
        physics_variables.psyncpv = radpwr_data.psyncpv
        physics_variables.pcoreradpv = radpwr_data.pcoreradpv
        physics_variables.pedgeradpv = radpwr_data.pedgeradpv
        physics_variables.pradpv = radpwr_data.pradpv

        physics_variables.pcoreradpv = max(physics_variables.pcoreradpv, 0.0e0)
        physics_variables.pedgeradpv = max(physics_variables.pedgeradpv, 0.0e0)

        physics_variables.pinnerzoneradmw = (
            physics_variables.pcoreradpv * physics_variables.vol
        )  # Should probably be vol_core
        physics_variables.pouterzoneradmw = (
            physics_variables.pedgeradpv * physics_variables.vol
        )

        physics_variables.pradmw = physics_variables.pradpv * physics_variables.vol

        #  Heating power to plasma (= Psol in divertor model)
        #  Ohmic power is zero in a stellarator
        #  physics_variables.pradmw here is core + edge (no SOL)

        powht = (
            physics_variables.falpha * physics_variables.palpmw
            + physics_variables.pchargemw
            + physics_variables.pohmmw
            - physics_variables.pradpv * physics_variables.vol
        )
        powht = max(
            0.00001e0, powht
        )  # To avoid negative heating power. This line is VERY important

        if physics_variables.ignite == 0:
            powht = (
                powht + current_drive_variables.pinjmw
            )  # if not ignited add the auxiliary power

        # Here the implementation sometimes leaves the accessible regime when pradmw> powht which is unphysical and
        # is not taken care of by the rad module. We restrict the radiation power here by the heating power:
        physics_variables.pradmw = max(0.0e0, physics_variables.pradmw)

        #  Power to divertor, = (1-stellarator_variables.f_rad)*Psol

        # The SOL radiation needs to be smaller than the physics_variables.pradmw
        physics_variables.psolradmw = stellarator_variables.f_rad * powht
        physics_variables.pdivt = powht - physics_variables.psolradmw

        # Add SOL Radiation to total
        physics_variables.pradmw = (
            physics_variables.pradmw + physics_variables.psolradmw
        )
        # pradpv = physics_variables.pradmw / physics_variables.vol # this line OVERWRITES the original definition of pradpv, probably shouldn't be defined like that as the core does not lose SOL power.

        #  The following line is unphysical, but prevents -ve sqrt argument
        #  Should be obsolete if constraint eqn 17 is turned on (but beware -
        #  this may not be quite correct for stellarators)
        physics_variables.pdivt = max(0.001e0, physics_variables.pdivt)

        #  Power transported to the first wall by escaped alpha particles

        physics_variables.palpfwmw = physics_variables.palpmw * (
            1.0e0 - physics_variables.falpha
        )

        # Nominal mean photon wall load
        if physics_variables.iwalld == 1:
            physics_variables.photon_wall = (
                physics_variables.ffwal
                * physics_variables.pradmw
                / physics_variables.sarea
            )
        else:
            if heat_transport_variables.ipowerflow == 0:
                physics_variables.photon_wall = (
                    (1.0e0 - fwbs_variables.fhole)
                    * physics_variables.pradmw
                    / build_variables.fwarea
                )
            else:
                physics_variables.photon_wall = (
                    (
                        1.0e0
                        - fwbs_variables.fhole
                        - fwbs_variables.fhcd
                        - fwbs_variables.fdiv
                    )
                    * physics_variables.pradmw
                    / build_variables.fwarea
                )

        constraint_variables.peakradwallload = (
            physics_variables.photon_wall * constraint_variables.peakfactrad
        )

        physics_variables.rad_fraction_total = physics_variables.pradmw / (
            physics_variables.falpha * physics_variables.palpmw
            + physics_variables.pchargemw
            + physics_variables.pohmmw
            + current_drive_variables.pinjmw
        )

        #  Calculate transport losses and energy confinement time using the
        #  chosen scaling law
        #  N.B. stellarator_variables.iotabar replaces tokamak physics_variables.q95 in argument list

        (
            physics_variables.kappaa,
            physics_variables.ptrepv,
            physics_variables.ptripv,
            physics_variables.tauee,
            physics_variables.tauei,
            physics_variables.taueff,
            physics_variables.powerht,
        ) = physics_module.pcond(
            physics_variables.afuel,
            physics_variables.palpmw,
            physics_variables.aspect,
            physics_variables.bt,
            physics_variables.dnitot,
            physics_variables.dene,
            physics_variables.dnla,
            physics_variables.eps,
            physics_variables.hfact,
            physics_variables.iinvqd,
            physics_variables.isc,
            physics_variables.ignite,
            physics_variables.kappa,
            physics_variables.kappa95,
            physics_variables.pchargemw,
            current_drive_variables.pinjmw,
            physics_variables.plascur,
            physics_variables.pcoreradpv,
            physics_variables.rmajor,
            physics_variables.rminor,
            physics_variables.te,
            physics_variables.ten,
            physics_variables.tin,
            stellarator_variables.iotabar,
            physics_variables.qstar,
            physics_variables.vol,
            physics_variables.xarea,
            physics_variables.zeff,
        )

        physics_variables.ptremw = physics_variables.ptrepv * physics_variables.vol
        physics_variables.ptrimw = physics_variables.ptripv * physics_variables.vol

        physics_variables.pscalingmw = (
            physics_variables.ptremw + physics_variables.ptrimw
        )

        #  Calculate auxiliary physics related information
        #  for the rest of the code

        sbar = 1.0e0
        (
            physics_variables.burnup,
            physics_variables.dntau,
            physics_variables.figmer,
            fusrat,
            physics_variables.qfuel,
            physics_variables.rndfuel,
            physics_variables.taup,
        ) = physics_module.phyaux(
            physics_variables.aspect,
            physics_variables.dene,
            physics_variables.deni,
            physics_variables.fusionrate,
            physics_variables.alpharate,
            physics_variables.plascur,
            sbar,
            physics_variables.dnalp,
            physics_variables.taueff,
            physics_variables.vol,
        )

        # Calculate physics_variables.beta limit. Does nothing atm so commented out
        # call stblim(physics_variables.betalim)

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
        neoclassics_module.init_neoclassics(
            0.6,
            stellarator_configuration.stella_config_epseff,
            stellarator_variables.iotabar,
        )

        q_PROCESS = (
            (
                physics_variables.falpha * physics_variables.palppv
                - physics_variables.pcoreradpv
            )
            * physics_variables.vol
            / physics_variables.sarea
            * impurity_radiation_module.coreradius
        )
        q_PROCESS_r1 = (
            (
                physics_variables.falpha * physics_variables.palppv
                - physics_variables.pcoreradpv
            )
            * physics_variables.vol
            / physics_variables.sarea
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
            * physics_variables.sarea
            * impurity_radiation_module.coreradius
        )
        dmdt_neo_fuel = (
            dndt_neo_fuel * physics_variables.afuel * constants.mproton * 1.0e6
        )  # mg
        dmdt_neo_fuel_from_e = (
            4
            * dndt_neo_e
            * physics_variables.sarea
            * impurity_radiation_module.coreradius
            * physics_variables.afuel
            * constants.mproton
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
            physics_variables.vol
            * st.f_r
            * (
                impurity_radiation_module.coreradius
                * physics_variables.rminor
                / stellarator_configuration.stella_config_rminor_ref
            )
            ** 2
        )
        surfacescaling = (
            physics_variables.sarea
            * st.f_r
            * (
                impurity_radiation_module.coreradius
                * physics_variables.rminor
                / stellarator_configuration.stella_config_rminor_ref
            )
        )

        nominator = (
            physics_variables.falpha * physics_variables.palppv
            - physics_variables.pcoreradpv
        ) * volscaling

        # in fortran there was a 0d0*alphan term which I have removed for obvious reasons
        # the following comment seems to describe this?
        # "include alphan if chi should be incorporate density gradients too"
        # but the history can be consulted if required (23/11/22 TN)
        denominator = (
            (
                3
                * physics_variables.ne0
                * constants.echarge
                * physics_variables.te0
                * 1e3
                * physics_variables.alphat
                * impurity_radiation_module.coreradius
                * (1 - impurity_radiation_module.coreradius**2)
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
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
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
            current_drive_variables.cnbeam = (
                1e-3
                * (current_drive_variables.pnbeam * 1e6)
                / current_drive_variables.enbeam
            )
        else:
            current_drive_variables.cnbeam = 0

        #  Ratio of fusion to input (injection+ohmic) power

        if (
            abs(
                current_drive_variables.pinjmw
                + current_drive_variables.porbitlossmw
                + physics_variables.pohmmw
            )
            < 1e-6
        ):
            current_drive_variables.bigq = 1e18
        else:
            current_drive_variables.bigq = physics_variables.powfmw / (
                current_drive_variables.pinjmw
                + current_drive_variables.porbitlossmw
                + physics_variables.pohmmw
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
                    "Neutral beam energy (keV)",
                    "(enbeam)",
                    current_drive_variables.enbeam,
                )
                po.ovarre(
                    self.outfile,
                    "Neutral beam current (A)",
                    "(cnbeam)",
                    current_drive_variables.cnbeam,
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
