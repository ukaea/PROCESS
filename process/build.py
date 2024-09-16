from process.fortran import tfcoil_variables
from process.fortran import divertor_variables
from process.fortran import current_drive_variables
from process.fortran import physics_variables
from process.fortran import maths_library
from process.fortran import constants
from process.fortran import pfcoil_variables
from process.fortran import build_variables
from process.fortran import numerics
from process.fortran import fwbs_variables
from process.fortran import error_handling
from process.fortran import process_output as po
from process.fortran import buildings_variables
from process.variables import AnnotatedVariable
import numpy as np
import logging

logger = logging.getLogger(__name__)


class Build:
    def __init__(self):
        self.outfile = constants.nout
        self.mfile = constants.mfile
        self.ripflag = AnnotatedVariable(int, 0)

    def portsz(self):
        """Port size calculation
        author: P J Knight, CCFE, Culham Science Centre
        author: M D Kovari, CCFE, Culham Science Centre
        None
        This subroutine finds the maximum possible tangency radius
        for adequate beam access.
        <P>The outputs from the routine are
        <UL> <P><LI>rtanbeam : Beam tangency radius (m)
        <P><LI>rtanmax : Maximum possible tangency radius (m) </UL>
        A User's Guide to the PROCESS Systems Code
        """
        current_drive_variables.rtanbeam = (
            current_drive_variables.frbeam * physics_variables.rmajor
        )

        #  Toroidal angle between adjacent TF coils

        omega = constants.twopi / tfcoil_variables.n_tf

        #  Half-width of outboard TF coil in toroidal direction (m)
        a = 0.5e0 * tfcoil_variables.tftort  # (previously used inboard leg width)
        try:
            assert a < np.inf
        except AssertionError:
            logger.exception("a is inf. Kludging to 1e10.")
            a = 1e10

        #  Radial thickness of outboard TF coil leg (m)
        b = build_variables.tfthko
        try:
            assert b < np.inf
        except AssertionError:
            logger.exception("b is inf. Kludging to 1e10.")
            b = 1e10

        #  Width of beam duct, including shielding on both sides (m)
        c = current_drive_variables.beamwd + 2.0e0 * current_drive_variables.nbshield

        #  Major radius of inner edge of outboard TF coil (m)
        d = build_variables.r_tf_outboard_mid - 0.5e0 * b
        try:
            assert d < np.inf
        except AssertionError:
            logger.exception("d is inf. Kludging to 1e10.")
            d = 1e10

        #  Refer to figure in User Guide for remaining geometric calculations
        e = np.sqrt(a * a + (d + b) * (d + b))
        f = np.sqrt(a * a + d * d)

        theta = omega - np.arctan(a / d)
        phi = theta - np.arcsin(a / e)

        g = np.sqrt(e * e + f * f - 2.0e0 * e * f * np.cos(phi))  # cosine rule

        if g > c:

            h = np.sqrt(g * g - c * c)

            alpha = np.arctan(h / c)
            eps = np.arcsin(e * np.sin(phi) / g) - alpha  # from sine rule

            #  Maximum tangency radius for centreline of beam (m)

            current_drive_variables.rtanmax = f * np.cos(eps) - 0.5e0 * c

        else:  # coil separation is too narrow for beam...

            error_handling.fdiags[0] = g
            error_handling.fdiags[1] = c
            error_handling.report_error(63)

            current_drive_variables.rtanmax = 0.0e0

    def calculate_vertical_build(self, output: bool) -> None:
        """
        This method determines the vertical build of the machine.
        It calculates various parameters related to the build of the machine,
        such as thicknesses, radii, and areas.
        Results can be outputted with the `output` flag.

        Args:
            output (bool): Flag indicating whether to output results

        Returns:
            None

        """
        if output:
            po.oheadr(self.outfile, "Vertical Build")

            po.ovarin(
                self.mfile,
                "Divertor null switch",
                "(i_single_null)",
                physics_variables.i_single_null,
            )

            if physics_variables.i_single_null == 0:
                po.ocmmnt(self.outfile, "Double null case")

                # Start at the top and work down.

                vbuild = (
                    buildings_variables.clh1
                    + build_variables.tfcth
                    + build_variables.tftsgap
                    + build_variables.thshield_vb
                    + build_variables.vgap2
                    + build_variables.d_vv_top
                    + build_variables.shldtth
                    + divertor_variables.divfix
                    + build_variables.vgaptop
                    + physics_variables.rminor * physics_variables.kappa
                )

                # To calculate vertical offset between TF coil centre and plasma centre
                vbuile1 = vbuild

                po.obuild(
                    self.outfile,
                    "Cryostat roof structure*",
                    buildings_variables.clh1,
                    vbuild,
                    "(clh1)",
                )
                po.ovarre(
                    self.mfile,
                    "Cryostat roof structure*",
                    "(clh1)",
                    buildings_variables.clh1,
                )
                vbuild = vbuild - buildings_variables.clh1

                # Top of TF coil
                tf_top = vbuild

                po.obuild(
                    self.outfile,
                    "TF coil",
                    build_variables.tfcth,
                    vbuild,
                    "(tfcth)",
                )
                vbuild = vbuild - build_variables.tfcth

                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.tftsgap,
                    vbuild,
                    "(tftsgap)",
                )
                vbuild = vbuild - build_variables.tftsgap

                po.obuild(
                    self.outfile,
                    "Thermal shield, vertical",
                    build_variables.thshield_vb,
                    vbuild,
                    "(thshield_vb)",
                )

                po.ovarre(
                    self.mfile,
                    "Thermal shield, vertical (m)",
                    "(thshield_vb)",
                    build_variables.thshield_vb,
                )
                vbuild = vbuild - build_variables.thshield_vb

                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.vgap2,
                    vbuild,
                    "(vgap2)",
                )
                po.ovarre(
                    self.mfile,
                    "Vessel - TF coil vertical gap (m)",
                    "(vgap2)",
                    build_variables.vgap2,
                )
                vbuild = vbuild - build_variables.vgap2

                po.obuild(
                    self.outfile,
                    "Vacuum vessel (and shielding)",
                    build_variables.d_vv_top + build_variables.shldtth,
                    vbuild,
                    "(d_vv_top+shldtth)",
                )
                vbuild = vbuild - build_variables.d_vv_top - build_variables.shldtth
                po.ovarre(
                    self.mfile,
                    "Topside vacuum vessel radial thickness (m)",
                    "(d_vv_top)",
                    build_variables.d_vv_top,
                )
                po.ovarre(
                    self.mfile,
                    "Top radiation shield thickness (m)",
                    "(shldtth)",
                    build_variables.shldtth,
                )

                po.obuild(
                    self.outfile,
                    "Divertor structure",
                    divertor_variables.divfix,
                    vbuild,
                    "(divfix)",
                )
                po.ovarre(
                    self.mfile,
                    "Divertor structure vertical thickness (m)",
                    "(divfix)",
                    divertor_variables.divfix,
                )
                vbuild = vbuild - divertor_variables.divfix

                po.obuild(
                    self.outfile,
                    "Top scrape-off",
                    build_variables.vgaptop,
                    vbuild,
                    "(vgaptop)",
                )
                po.ovarre(
                    self.mfile,
                    "Top scrape-off vertical thickness (m)",
                    "(vgaptop)",
                    build_variables.vgaptop,
                )
                vbuild = vbuild - build_variables.vgaptop

                po.obuild(
                    self.outfile,
                    "Plasma top",
                    physics_variables.rminor * physics_variables.kappa,
                    vbuild,
                    "(rminor*kappa)",
                )
                po.ovarre(
                    self.mfile,
                    "Plasma half-height (m)",
                    "(rminor*kappa)",
                    physics_variables.rminor * physics_variables.kappa,
                )
                vbuild = vbuild - physics_variables.rminor * physics_variables.kappa

                po.obuild(self.outfile, "Midplane", 0.0e0, vbuild)

                vbuild = vbuild - physics_variables.rminor * physics_variables.kappa
                po.obuild(
                    self.outfile,
                    "Plasma bottom",
                    physics_variables.rminor * physics_variables.kappa,
                    vbuild,
                    "(rminor*kappa)",
                )

                vbuild = vbuild - build_variables.vgap
                po.obuild(
                    self.outfile,
                    "Lower scrape-off",
                    build_variables.vgap,
                    vbuild,
                    "(vgap)",
                )
                po.ovarre(
                    self.mfile,
                    "Bottom scrape-off vertical thickness (m)",
                    "(vgap)",
                    build_variables.vgap,
                )

                vbuild = vbuild - divertor_variables.divfix
                po.obuild(
                    self.outfile,
                    "Divertor structure",
                    divertor_variables.divfix,
                    vbuild,
                    "(divfix)",
                )
                po.ovarre(
                    self.mfile,
                    "Divertor structure vertical thickness (m)",
                    "(divfix)",
                    divertor_variables.divfix,
                )

                vbuild = vbuild - build_variables.shldlth

                vbuild = vbuild - build_variables.d_vv_bot
                po.obuild(
                    self.outfile,
                    "Vacuum vessel (and shielding)",
                    build_variables.d_vv_bot + build_variables.shldlth,
                    vbuild,
                    "(d_vv_bot+shldlth)",
                )
                po.ovarre(
                    self.mfile,
                    "Bottom radiation shield thickness (m)",
                    "(shldlth)",
                    build_variables.shldlth,
                )
                po.ovarre(
                    self.mfile,
                    "Underside vacuum vessel radial thickness (m)",
                    "(d_vv_bot)",
                    build_variables.d_vv_bot,
                )

                vbuild = vbuild - build_variables.vgap2
                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.vgap2,
                    vbuild,
                    "(vgap2)",
                )

                vbuild = vbuild - build_variables.thshield_vb
                po.obuild(
                    self.outfile,
                    "Thermal shield, vertical",
                    build_variables.thshield_vb,
                    vbuild,
                    "(thshield_vb)",
                )

                vbuild = vbuild - build_variables.tftsgap
                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.tftsgap,
                    vbuild,
                    "(tftsgap)",
                )

                vbuild = vbuild - build_variables.tfcth
                po.obuild(
                    self.outfile,
                    "TF coil",
                    build_variables.tfcth,
                    vbuild,
                    "(tfcth)",
                )

                # Total height of TF coil
                tf_height = tf_top - vbuild
                # Inner vertical dimension of TF coil
                build_variables.dh_tf_inner_bore = tf_height - 2 * build_variables.tfcth

                vbuild = vbuild - buildings_variables.clh1
                po.obuild(
                    self.outfile,
                    "Cryostat floor structure**",
                    buildings_variables.clh1,
                    vbuild,
                    "(clh1)",
                )

                # To calculate vertical offset between TF coil centre and plasma centre
                build_variables.tfoffset = (vbuile1 + vbuild) / 2.0e0

                # End of Double null case
            else:
                #  po.ocmmnt(self.outfile, "Single null case")
                #  write(self.outfile, 20)

                vbuild = (
                    buildings_variables.clh1
                    + build_variables.tfcth
                    + build_variables.tftsgap
                    + build_variables.thshield_vb
                    + build_variables.vgap2
                    + 0.5e0 * (build_variables.d_vv_top + build_variables.d_vv_bot)
                    + build_variables.vvblgap
                    + build_variables.shldtth
                    + build_variables.blnktth
                    + 0.5e0 * (build_variables.fwith + build_variables.fwoth)
                    + build_variables.vgaptop
                    + physics_variables.rminor * physics_variables.kappa
                )

                # To calculate vertical offset between TF coil centre and plasma centre
                vbuile1 = vbuild

                po.obuild(
                    self.outfile,
                    "Cryostat roof structure*",
                    buildings_variables.clh1,
                    vbuild,
                    "(clh1)",
                )
                po.ovarre(
                    self.mfile,
                    "Cryostat roof structure*",
                    "(clh1)",
                    buildings_variables.clh1,
                )
                vbuild = vbuild - buildings_variables.clh1

                # Top of TF coil
                tf_top = vbuild

                po.obuild(
                    self.outfile,
                    "TF coil",
                    build_variables.tfcth,
                    vbuild,
                    "(tfcth)",
                )
                vbuild = vbuild - build_variables.tfcth

                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.tftsgap,
                    vbuild,
                    "(tftsgap)",
                )
                vbuild = vbuild - build_variables.tftsgap

                po.obuild(
                    self.outfile,
                    "Thermal shield, vertical",
                    build_variables.thshield_vb,
                    vbuild,
                    "(thshield_vb)",
                )
                po.ovarre(
                    self.mfile,
                    "Thermal shield, vertical (m)",
                    "(thshield_vb)",
                    build_variables.thshield_vb,
                )
                vbuild = vbuild - build_variables.thshield_vb

                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.vgap2,
                    vbuild,
                    "(vgap2)",
                )
                po.ovarre(
                    self.mfile,
                    "Vessel - TF coil vertical gap (m)",
                    "(vgap2)",
                    build_variables.vgap2,
                )
                vbuild = vbuild - build_variables.vgap2

                po.obuild(
                    self.outfile,
                    "Vacuum vessel (and shielding)",
                    build_variables.d_vv_top + build_variables.shldtth,
                    vbuild,
                    "(d_vv_top+shldtth)",
                )
                vbuild = vbuild - build_variables.d_vv_top - build_variables.shldtth
                po.ovarre(
                    self.mfile,
                    "Topside vacuum vessel radial thickness (m)",
                    "(d_vv_top)",
                    build_variables.d_vv_top,
                )
                po.ovarre(
                    self.mfile,
                    "Top radiation shield thickness (m)",
                    "(shldtth)",
                    build_variables.shldtth,
                )

                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.vvblgap,
                    vbuild,
                    "(vvblgap)",
                )
                vbuild = vbuild - build_variables.vvblgap

                po.obuild(
                    self.outfile,
                    "Top blanket",
                    build_variables.blnktth,
                    vbuild,
                    "(blnktth)",
                )
                po.ovarre(
                    self.mfile,
                    "Top blanket vertical thickness (m)",
                    "(blnktth)",
                    build_variables.blnktth,
                )
                vbuild = vbuild - build_variables.blnktth

                fwtth = 0.5e0 * (build_variables.fwith + build_variables.fwoth)
                po.obuild(self.outfile, "Top first wall", fwtth, vbuild, "(fwtth)")
                po.ovarre(
                    self.mfile,
                    "Top first wall vertical thickness (m)",
                    "(fwtth)",
                    fwtth,
                )
                vbuild = vbuild - fwtth

                po.obuild(
                    self.outfile,
                    "Top scrape-off",
                    build_variables.vgaptop,
                    vbuild,
                    "(vgaptop)",
                )
                po.ovarre(
                    self.mfile,
                    "Top scrape-off vertical thickness (m)",
                    "(vgaptop)",
                    build_variables.vgaptop,
                )
                vbuild = vbuild - build_variables.vgaptop

                po.obuild(
                    self.outfile,
                    "Plasma top",
                    physics_variables.rminor * physics_variables.kappa,
                    vbuild,
                    "(rminor*kappa)",
                )
                po.ovarre(
                    self.mfile,
                    "Plasma half-height (m)",
                    "(rminor*kappa)",
                    physics_variables.rminor * physics_variables.kappa,
                )
                vbuild = vbuild - physics_variables.rminor * physics_variables.kappa

                po.obuild(self.outfile, "Midplane", 0.0e0, vbuild)

                vbuild = vbuild - physics_variables.rminor * physics_variables.kappa
                po.obuild(
                    self.outfile,
                    "Plasma bottom",
                    physics_variables.rminor * physics_variables.kappa,
                    vbuild,
                    "(rminor*kappa)",
                )

                vbuild = vbuild - build_variables.vgap
                po.obuild(
                    self.outfile,
                    "Lower scrape-off",
                    build_variables.vgap,
                    vbuild,
                    "(vgap)",
                )
                po.ovarre(
                    self.mfile,
                    "Bottom scrape-off vertical thickness (m)",
                    "(vgap)",
                    build_variables.vgap,
                )

                vbuild = vbuild - divertor_variables.divfix
                po.obuild(
                    self.outfile,
                    "Divertor structure",
                    divertor_variables.divfix,
                    vbuild,
                    "(divfix)",
                )
                po.ovarre(
                    self.mfile,
                    "Divertor structure vertical thickness (m)",
                    "(divfix)",
                    divertor_variables.divfix,
                )

                vbuild = vbuild - build_variables.shldlth

                vbuild = vbuild - build_variables.d_vv_bot
                po.obuild(
                    self.outfile,
                    "Vacuum vessel (and shielding)",
                    build_variables.d_vv_bot + build_variables.shldlth,
                    vbuild,
                    "(d_vv_bot+shldlth)",
                )
                po.ovarre(
                    self.mfile,
                    "Bottom radiation shield thickness (m)",
                    "(shldlth)",
                    build_variables.shldlth,
                )
                po.ovarre(
                    self.mfile,
                    "Underside vacuum vessel radial thickness (m)",
                    "(d_vv_bot)",
                    build_variables.d_vv_bot,
                )

                vbuild = vbuild - build_variables.vgap2
                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.vgap2,
                    vbuild,
                    "(vgap2)",
                )

                vbuild = vbuild - build_variables.thshield_vb
                po.obuild(
                    self.outfile,
                    "Thermal shield, vertical",
                    build_variables.thshield_vb,
                    vbuild,
                    "(thshield_vb)",
                )

                vbuild = vbuild - build_variables.tftsgap
                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.tftsgap,
                    vbuild,
                    "(tftsgap)",
                )

                vbuild = vbuild - build_variables.tfcth
                po.obuild(
                    self.outfile,
                    "TF coil",
                    build_variables.tfcth,
                    vbuild,
                    "(tfcth)",
                )

                # Total height of TF coil
                tf_height = tf_top - vbuild
                # Inner vertical dimension of TF coil
                build_variables.dh_tf_inner_bore = tf_height - 2 * build_variables.tfcth

                vbuild = vbuild - buildings_variables.clh1

                po.obuild(
                    self.outfile,
                    "Cryostat floor structure**",
                    buildings_variables.clh1,
                    vbuild,
                    "(clh1)",
                )

                # To calculate vertical offset between TF coil centre and plasma centre
                build_variables.tfoffset = (vbuile1 + vbuild) / 2.0e0

                # end of Single null case

            po.ovarre(
                self.mfile,
                "Ratio of Central Solenoid height to TF coil internal height",
                "(ohhghf)",
                pfcoil_variables.ohhghf,
            )
            po.ocmmnt(
                self.outfile,
                "\n*Cryostat roof allowance includes uppermost PF coil and outer thermal shield.\n*Cryostat floor allowance includes lowermost PF coil, outer thermal shield and gravity support.",
            )

        #  Other build quantities

        divht = self.divgeom(output)
        # Issue #481 Remove build_variables.vgaptf
        if build_variables.vgap < 0.00001e0:
            build_variables.vgap = divht

        # If build_variables.vgap /= 0 use the value set by the user.

        # Height to inside edge of TF coil. TF coils are assumed to be symmetrical.
        # Therefore this applies to single and double null cases.
        build_variables.hmax = (
            physics_variables.rminor * physics_variables.kappa
            + build_variables.vgap
            + divertor_variables.divfix
            + build_variables.shldlth
            + build_variables.d_vv_bot
            + build_variables.vgap2
            + build_variables.thshield_vb
            + build_variables.tftsgap
        )

        #  Vertical locations of divertor coils
        if physics_variables.i_single_null == 0:
            build_variables.hpfu = build_variables.hmax + build_variables.tfcth
            build_variables.hpfdif = 0.0e0
        else:
            build_variables.hpfu = (
                build_variables.tfcth
                + build_variables.tftsgap
                + build_variables.thshield_vb
                + build_variables.vgap2
                + build_variables.d_vv_top
                + build_variables.shldtth
                + build_variables.vvblgap
                + build_variables.blnktth
                + 0.5e0 * (build_variables.fwith + build_variables.fwoth)
                + build_variables.vgaptop
                + physics_variables.rminor * physics_variables.kappa
            )
            build_variables.hpfdif = (
                build_variables.hpfu - (build_variables.hmax + build_variables.tfcth)
            ) / 2.0e0

    def divgeom(self, output: bool):
        """
                Divertor geometry calculation
        author: J Galambos, ORNL
        author: P J Knight, CCFE, Culham Science Centre
        divht : output real : divertor height (m)
        self.outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This subroutine determines the divertor geometry.
        The inboard (i) and outboard (o) plasma surfaces
        are approximated by arcs, and followed past the X-point to
        determine the maximum height.
        TART option: Peng SOFT paper
        """
        if physics_variables.itart == 1:
            return 1.75e0 * physics_variables.rminor
        #  Conventional tokamak divertor model
        #  options for seperate upper and lower physics_variables.triangularity

        kap = physics_variables.kappa
        triu = physics_variables.triang
        tril = physics_variables.triang

        # Old method: assumes that divertor arms are continuations of arcs
        #
        # Outboard side
        # build_variables.plsepo = poloidal length along the separatrix from null to
        # strike point on outboard [default 1.5 m]
        # thetao = arc angle between the strike point and the null point
        #
        # xpointo = physics_variables.rmajor + 0.5e0*physics_variables.rminor*(kap**2 + tri**2 - 1.0e0) /     #     (1.0e0 - tri)
        # rprimeo = (xpointo - physics_variables.rmajor + physics_variables.rminor)
        # phio = asin(kap*physics_variables.rminor/rprimeo)
        # thetao = build_variables.plsepo/rprimeo
        #
        # Initial strike point
        #
        # yspointo = rprimeo * sin(thetao + phio)
        # xspointo = xpointo - rprimeo * cos(thetao + phio)
        #
        # Outboard strike point radius - normalized to ITER
        #
        # rstrko = xspointo + 0.14e0
        #
        # Uppermost divertor strike point (end of power decay)
        # anginc = angle of incidence of scrape-off field lines on the
        # divertor (rad)
        #
        # +**PJK 25/07/11 Changed sign of anginc contribution
        # yprimeb = soleno * cos(thetao + phio - anginc)
        #
        # divht = yprimeb + yspointo - kap*physics_variables.rminor

        # New method, assuming straight legs -- superceded by new method 26/5/2016
        # Assumed 90 degrees at X-pt -- wrong#
        #
        #  Find half-angle of outboard arc
        # denomo = (tril**2 + kap**2 - 1.0e0)/( 2.0e0*(1.0e0+tril) ) - tril
        # thetao = atan(kap/denomo)
        # Angle between horizontal and inner divertor leg
        # alphad = (pi/2.0e0) - thetao

        # Method 26/05/2016
        # Find radius of inner and outer plasma arcs

        rco = 0.5 * np.sqrt(
            (physics_variables.rminor**2 * ((tril + 1.0e0) ** 2 + kap**2) ** 2)
            / ((tril + 1.0e0) ** 2)
        )
        rci = 0.5 * np.sqrt(
            (physics_variables.rminor**2 * ((tril - 1.0e0) ** 2 + kap**2) ** 2)
            / ((tril - 1.0e0) ** 2)
        )

        # Find angles between vertical and legs
        # Inboard arc angle = outboard leg angle

        thetao = np.arcsin(1.0e0 - (physics_variables.rminor * (1.0e0 - tril)) / rci)

        # Outboard arc angle = inboard leg angle

        thetai = np.arcsin(1.0e0 - (physics_variables.rminor * (1.0e0 + tril)) / rco)

        #  Position of lower x-pt
        rxpt = physics_variables.rmajor - tril * physics_variables.rminor
        zxpt = -1.0e0 * kap * physics_variables.rminor

        # Position of inner strike point
        # rspi = rxpt - build_variables.plsepi*cos(alphad)
        # zspi = zxpt - build_variables.plsepi*sin(alphad)
        rspi = rxpt - build_variables.plsepi * np.cos(thetai)
        zspi = zxpt - build_variables.plsepi * np.sin(thetai)

        # Position of outer strike point
        # build_variables.rspo = rxpt + build_variables.plsepo*cos((pi/2.0e0)-alphad)
        # zspo = zxpt - build_variables.plsepo*sin((pi/2.0e0)-alphad)
        build_variables.rspo = rxpt + build_variables.plsepo * np.cos(thetao)
        zspo = zxpt - build_variables.plsepo * np.sin(thetao)

        # Position of inner plate ends
        # rplti = rspi - (build_variables.plleni/2.0e0)*sin(divertor_variables.betai + alphad - pi/2.0e0)
        # zplti = zspi + (build_variables.plleni/2.0e0)*cos(divertor_variables.betai + alphad - pi/2.0e0)
        # rplbi = rspi + (build_variables.plleni/2.0e0)*sin(divertor_variables.betai + alphad - pi/2.0e0)
        # zplbi = zspi - (build_variables.plleni/2.0e0)*cos(divertor_variables.betai + alphad - pi/2.0e0)
        rplti = rspi + (build_variables.plleni / 2.0e0) * np.cos(
            thetai + divertor_variables.betai
        )
        zplti = zspi + (build_variables.plleni / 2.0e0) * np.sin(
            thetai + divertor_variables.betai
        )
        rplbi = rspi - (build_variables.plleni / 2.0e0) * np.cos(
            thetai + divertor_variables.betai
        )
        zplbi = zspi - (build_variables.plleni / 2.0e0) * np.sin(
            thetai + divertor_variables.betai
        )

        # Position of outer plate ends
        # rplto = build_variables.rspo + (build_variables.plleno/2.0e0)*sin(divertor_variables.betao - alphad)
        # zplto = zspo + (build_variables.plleno/2.0e0)*cos(divertor_variables.betao - alphad)
        # rplbo = build_variables.rspo - (build_variables.plleno/2.0e0)*sin(divertor_variables.betao - alphad)
        # zplbo = zspo - (build_variables.plleno/2.0e0)*cos(divertor_variables.betao - alphad)
        rplto = build_variables.rspo - (build_variables.plleno / 2.0e0) * np.cos(
            thetao + divertor_variables.betao
        )
        zplto = zspo + (build_variables.plleno / 2.0e0) * np.sin(
            thetao + divertor_variables.betao
        )
        rplbo = build_variables.rspo + (build_variables.plleno / 2.0e0) * np.cos(
            thetao + divertor_variables.betao
        )
        zplbo = zspo - (build_variables.plleno / 2.0e0) * np.sin(
            thetao + divertor_variables.betao
        )

        divht = max(zplti, zplto) - min(zplbo, zplbi)

        if output:
            if physics_variables.idivrt == 1:
                po.oheadr(self.outfile, "Divertor build and plasma position")
                po.ocmmnt(self.outfile, "Divertor Configuration = Single Null Divertor")
                po.oblnkl(self.outfile)
                ptop_radial = physics_variables.rmajor - triu * physics_variables.rminor
                ptop_vertical = kap * physics_variables.rminor
                po.ovarrf(
                    self.outfile,
                    "Plasma top position, radial (m)",
                    "(ptop_radial)",
                    ptop_radial,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma top position, vertical (m)",
                    "(ptop_vertical)",
                    ptop_vertical,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma geometric centre, radial (m)",
                    "(physics_variables.rmajor.)",
                    physics_variables.rmajor,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma geometric centre, vertical (m)",
                    "(0.0)",
                    0.0e0,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma lower physics_variables.triangularity",
                    "(tril)",
                    tril,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma elongation",
                    "(physics_variables.kappa.)",
                    kap,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "TF coil vertical offset (m)",
                    "(tfoffset)",
                    build_variables.tfoffset,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma outer arc radius of curvature (m)",
                    "(rco)",
                    rco,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma inner arc radius of curvature (m)",
                    "(rci)",
                    rci,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile, "Plasma lower X-pt, radial (m)", "(rxpt)", rxpt, "OP "
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma lower X-pt, vertical (m)",
                    "(zxpt)",
                    zxpt,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Poloidal plane angle between vertical and inner leg (rad)",
                    "(thetai)",
                    thetai,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Poloidal plane angle between vertical and outer leg (rad)",
                    "(thetao)",
                    thetao,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Poloidal plane angle between inner leg and plate (rad)",
                    "(betai)",
                    divertor_variables.betai,
                )
                po.ovarrf(
                    self.outfile,
                    "Poloidal plane angle between outer leg and plate (rad)",
                    "(betao)",
                    divertor_variables.betao,
                )
                po.ovarrf(
                    self.outfile,
                    "Inner divertor leg poloidal length (m)",
                    "(plsepi)",
                    build_variables.plsepi,
                )
                po.ovarrf(
                    self.outfile,
                    "Outer divertor leg poloidal length (m)",
                    "(plsepo)",
                    build_variables.plsepo,
                )
                po.ovarrf(
                    self.outfile,
                    "Inner divertor plate length (m)",
                    "(plleni)",
                    build_variables.plleni,
                )
                po.ovarrf(
                    self.outfile,
                    "Outer divertor plate length (m)",
                    "(plleno)",
                    build_variables.plleno,
                )
                po.ovarrf(
                    self.outfile,
                    "Inner strike point, radial (m)",
                    "(rspi)",
                    rspi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Inner strike point, vertical (m)",
                    "(zspi)",
                    zspi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile, "Inner plate top, radial (m)", "(rplti)", rplti, "OP "
                )
                po.ovarrf(
                    self.outfile,
                    "Inner plate top, vertical (m)",
                    "(zplti)",
                    zplti,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Inner plate bottom, radial (m)",
                    "(rplbi)",
                    rplbi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Inner plate bottom, vertical (m)",
                    "(zplbi)",
                    zplbi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Outer strike point, radial (m)",
                    "(rspo)",
                    build_variables.rspo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Outer strike point, vertical (m)",
                    "(zspo)",
                    zspo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile, "Outer plate top, radial (m)", "(rplto)", rplto, "OP "
                )
                po.ovarrf(
                    self.outfile,
                    "Outer plate top, vertical (m)",
                    "(zplto)",
                    zplto,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Outer plate bottom, radial (m)",
                    "(rplbo)",
                    rplbo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Outer plate bottom, vertical (m)",
                    "(zplbo)",
                    zplbo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Calculated maximum divertor height (m)",
                    "(divht)",
                    divht,
                    "OP ",
                )

            elif physics_variables.idivrt == 2:
                po.oheadr(self.outfile, "Divertor build and plasma position")
                po.ocmmnt(self.outfile, "Divertor Configuration = Double Null Divertor")
                po.oblnkl(self.outfile)
                # Assume upper and lower divertors geometries are symmetric.
                ptop_radial = physics_variables.rmajor - triu * physics_variables.rminor
                ptop_vertical = kap * physics_variables.rminor
                po.ovarrf(
                    self.outfile,
                    "Plasma top position, radial (m)",
                    "(ptop_radial)",
                    ptop_radial,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma top position, vertical (m)",
                    "(ptop_vertical)",
                    ptop_vertical,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma geometric centre, radial (m)",
                    "(rmajor.)",
                    physics_variables.rmajor,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma geometric centre, vertical (m)",
                    "(0.0)",
                    0.0e0,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma physics_variables.triangularity",
                    "(tril)",
                    tril,
                    "OP ",
                )
                po.ovarrf(self.outfile, "Plasma elongation", "(kappa.)", kap, "OP ")
                po.ovarrf(
                    self.outfile,
                    "TF coil vertical offset (m)",
                    "(tfoffset)",
                    build_variables.tfoffset,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile, "Plasma upper X-pt, radial (m)", "(rxpt)", rxpt, "OP "
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma upper X-pt, vertical (m)",
                    "(-zxpt)",
                    -zxpt,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma outer arc radius of curvature (m)",
                    "(rco)",
                    rco,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma inner arc radius of curvature (m)",
                    "(rci)",
                    rci,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile, "Plasma lower X-pt, radial (m)", "(rxpt)", rxpt, "OP "
                )
                po.ovarrf(
                    self.outfile,
                    "Plasma lower X-pt, vertical (m)",
                    "(zxpt)",
                    zxpt,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Poloidal plane angle between vertical and inner leg (rad)",
                    "(thetai)",
                    thetai,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Poloidal plane angle between vertical and outer leg (rad)",
                    "(thetao)",
                    thetao,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Poloidal plane angle between inner leg and plate (rad)",
                    "(betai)",
                    divertor_variables.betai,
                )
                po.ovarrf(
                    self.outfile,
                    "Poloidal plane angle between outer leg and plate (rad)",
                    "(betao)",
                    divertor_variables.betao,
                )
                po.ovarrf(
                    self.outfile,
                    "Inner divertor leg poloidal length (m)",
                    "(plsepi)",
                    build_variables.plsepi,
                )
                po.ovarrf(
                    self.outfile,
                    "Outer divertor leg poloidal length (m)",
                    "(plsepo)",
                    build_variables.plsepo,
                )
                po.ovarrf(
                    self.outfile,
                    "Inner divertor plate length (m)",
                    "(lleni)",
                    build_variables.plleni,
                )
                po.ovarrf(
                    self.outfile,
                    "Outer divertor plate length (m)",
                    "(plleno)",
                    build_variables.plleno,
                )
                po.ovarrf(
                    self.outfile,
                    "Upper inner strike point, radial (m)",
                    "(rspi)",
                    rspi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper inner strike point, vertical (m)",
                    "(-zspi)",
                    -zspi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper inner plate top, radial (m)",
                    "(rplti)",
                    rplti,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper inner plate top, vertical (m)",
                    "(-zplti)",
                    -zplti,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper inner plate bottom, radial (m)",
                    "(rplbi)",
                    rplbi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper inner plate bottom, vertical (m)",
                    "(-zplbi)",
                    -zplbi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper outer strike point, radial (m)",
                    "(rspo)",
                    build_variables.rspo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper outer strike point, vertical (m)",
                    "(-zspo)",
                    -zspo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper outer plate top, radial (m)",
                    "(rplto)",
                    rplto,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper outer plate top, vertical (m)",
                    "(-zplto)",
                    -zplto,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper outer plate bottom, radial (m)",
                    "(rplbo)",
                    rplbo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Upper outer plate bottom, vertical (m)",
                    "(-zplbo)",
                    -zplbo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower inner strike point, radial (m)",
                    "(rspi)",
                    rspi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower inner strike point, vertical (m)",
                    "(zspi)",
                    zspi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower inner plate top, radial (m)",
                    "(rplti)",
                    rplti,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower inner plate top, vertical (m)",
                    "(zplti)",
                    zplti,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower inner plate bottom, radial (m)",
                    "(rplbi)",
                    rplbi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower inner plate bottom, vertical (m)",
                    "(zplbi)",
                    zplbi,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower outer strike point, radial (m)",
                    "(rspo)",
                    build_variables.rspo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower outer strike point, vertical (m)",
                    "(zspo)",
                    zspo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower outer plate top, radial (m)",
                    "(rplto)",
                    rplto,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower outer plate top, vertical (m)",
                    "(zplto)",
                    zplto,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower outer plate bottom, radial (m)",
                    "(rplbo)",
                    rplbo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Lower outer plate bottom, vertical (m)",
                    "(zplbo)",
                    zplbo,
                    "OP ",
                )
                po.ovarrf(
                    self.outfile,
                    "Calculated maximum divertor height (m)",
                    "(divht)",
                    divht,
                    "OP ",
                )
            else:
                po.oheadr(self.outfile, "Divertor build and plasma position")
                po.ocmmnt(
                    self.outfile,
                    "ERROR: null value not supported, check i_single_null value.",
                )
        return divht

    def ripple_amplitude(self, ripmax: float, r_tf_outboard_mid: float) -> float:
        """
        TF ripple calculation
        author: P J Knight and C W Ashe, CCFE, Culham Science Centre
        ripmax : input real  : maximum allowed ripple at plasma edge (%)
        ripple : output real : actual ripple at plasma edge (%)
        rtot   : input real  : radius to the centre of the outboard
        TF coil leg (m)
        rtotmin : output real : radius to the centre of the outboard
        TF coil leg which would produce
        a ripple of amplitude ripmax (m)
        flag : output integer : on exit, =1 if the fitted
        range of applicability is exceeded
        This routine calculates the toroidal field ripple amplitude
        at the midplane outboard plasma edge. The fitted coefficients
        were produced from MATLAB runs by M. Kovari using the CCFE
        MAGINT code to model the coils and fields.
        <P>The minimum radius of the centre of the TF coil legs
        to produce the maximum allowed ripple is also calculated.
        M. Kovari, Toroidal Field Coils - Maximum Field and Ripple -
        Parametric Calculation, July 2014
        ##############################################################

        Picture frame coil model by Ken McClements 2022 gives analytical
        solutions within 10% agreement with numerical models.
        Activated when i_tf_shape == 2 (picture frame)

        """
        n = float(tfcoil_variables.n_tf)
        if tfcoil_variables.i_tf_sup == 1:

            # Minimal inboard WP radius [m]
            r_wp_min = build_variables.r_tf_inboard_in + tfcoil_variables.thkcas

            # Rectangular WP
            if tfcoil_variables.i_tf_wp_geom == 0:
                r_wp_max = r_wp_min

            # Double rectangle WP
            elif tfcoil_variables.i_tf_wp_geom == 1:
                r_wp_max = r_wp_min + 0.5e0 * tfcoil_variables.dr_tf_wp

            # Trapezoidal WP
            elif tfcoil_variables.i_tf_wp_geom == 2:
                r_wp_max = r_wp_min + tfcoil_variables.dr_tf_wp

            # Calculated maximum toroidal WP toroidal thickness [m]
            if tfcoil_variables.tfc_sidewall_is_fraction:
                t_wp_max = 2.0e0 * (
                    (r_wp_max - tfcoil_variables.casths_fraction * r_wp_min)
                    * np.tan(np.pi / n)
                    - tfcoil_variables.tinstf
                    - tfcoil_variables.tfinsgap
                )
            else:
                t_wp_max = 2.0e0 * (
                    r_wp_max * np.tan(np.pi / n)
                    - tfcoil_variables.casths
                    - tfcoil_variables.tinstf
                    - tfcoil_variables.tfinsgap
                )

        # Resistive magnet case
        else:
            # Radius used to define the t_wp_max [m]
            r_wp_max = (
                build_variables.r_tf_inboard_in
                + tfcoil_variables.thkcas
                + tfcoil_variables.dr_tf_wp
            )

            # Calculated maximum toroidal WP toroidal thickness [m]
            t_wp_max = 2.0e0 * r_wp_max * np.tan(np.pi / n)

        flag = 0
        if tfcoil_variables.i_tf_shape == 2:
            # Ken McClements ST picture frame coil analytical ripple calc
            # Calculated ripple for coil at r_tf_outboard_mid (%)
            ripple = 100.0e0 * (
                (physics_variables.rmajor + physics_variables.rminor)
                / r_tf_outboard_mid
            ) ** (n)
            #  Calculated r_tf_outboard_mid to produce a ripple of amplitude ripmax
            r_tf_outboard_midmin = (
                physics_variables.rmajor + physics_variables.rminor
            ) / ((0.01e0 * ripmax) ** (1.0e0 / n))
        else:

            # Winding pack to iter-coil at plasma centre toroidal lenth ratio
            x = t_wp_max * n / physics_variables.rmajor

            # Fitting parameters
            c1 = 0.875e0 - 0.0557e0 * x
            c2 = 1.617e0 + 0.0832e0 * x

            #  Calculated ripple for coil at r_tf_outboard_mid (%)
            ripple = (
                100.0e0
                * c1
                * (
                    (physics_variables.rmajor + physics_variables.rminor)
                    / r_tf_outboard_mid
                )
                ** (n - c2)
            )

            #  Calculated r_tf_outboard_mid to produce a ripple of amplitude ripmax
            base = 0.01 * ripmax / c1
            # Avoid potential negative or complex result: kludge base to be
            # small and positive if required
            try:
                assert base > 1e-6
            except AssertionError:
                logger.exception("base is <= 1e-6. Kludging to 1e-6.")
                base = 1e-6

            r_tf_outboard_midmin = (
                physics_variables.rmajor + physics_variables.rminor
            ) / (base ** (1.0 / (n - c2)))

            try:
                assert r_tf_outboard_midmin < np.inf
            except AssertionError:
                logger.exception(
                    "r_tf_outboard_midmin is inf. Kludging to a large value instead."
                )
                r_tf_outboard_midmin = (
                    physics_variables.rmajor + physics_variables.rminor
                ) * 3

            #  Notify via flag if a range of applicability is violated
            flag = 0
            if (x < 0.737e0) or (x > 2.95e0):
                flag = 1
            if (tfcoil_variables.n_tf < 16) or (tfcoil_variables.n_tf > 20):
                flag = 2
            if (
                (physics_variables.rmajor + physics_variables.rminor)
                / r_tf_outboard_mid
                < 0.7e0
            ) or (
                (physics_variables.rmajor + physics_variables.rminor)
                / r_tf_outboard_mid
                > 0.8e0
            ):
                flag = 3

        return ripple, r_tf_outboard_midmin, flag

    def tf_in_cs_bore_calc(self):
        build_variables.bore += build_variables.tfcth + build_variables.gapoh

    def calculate_radial_build(self, output: bool) -> None:
        """
        This method determines the radial build of the machine.
        It calculates various parameters related to the build of the machine,
        such as thicknesses, radii, and areas.
        Results can be outputted with the `output` flag.

        Args:
            output (bool): Flag indicating whether to output the results

        Returns:
            None

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

        if physics_variables.i_single_null == 1:
            #  Check if build_variables.vgaptop has been set too small
            build_variables.vgaptop = max(
                0.5e0 * (build_variables.scrapli + build_variables.scraplo),
                build_variables.vgaptop,
            )

        # Calculate pre-compression structure thickness is build_variables.iprecomp=1
        if build_variables.iprecomp == 1:
            build_variables.precomp = build_variables.fseppc / (
                2.0e0
                * np.pi
                * build_variables.fcspc
                * build_variables.sigallpc
                * (build_variables.bore + build_variables.bore + build_variables.ohcth)
            )
        else:
            build_variables.precomp = 0.0e0

        if build_variables.tf_in_cs == 1:
            build_variables.r_tf_inboard_in = (
                build_variables.bore - build_variables.tfcth - build_variables.gapoh
            )
        else:
            # Inboard side inner radius [m]
            build_variables.r_tf_inboard_in = (
                build_variables.bore
                + build_variables.ohcth
                + build_variables.precomp
                + build_variables.gapoh
            )

        # Issue #514 Radial dimensions of inboard leg
        # Calculate build_variables.tfcth if tfcoil_variables.dr_tf_wp is an iteration variable (140)
        if any(numerics.ixc[0 : numerics.nvar] == 140):

            # SC TF coil thickness defined using its maximum (diagonal)
            if tfcoil_variables.i_tf_sup == 1:
                build_variables.tfcth = (
                    build_variables.r_tf_inboard_in
                    + tfcoil_variables.dr_tf_wp
                    + tfcoil_variables.casthi
                    + tfcoil_variables.thkcas
                ) / np.cos(
                    np.pi / tfcoil_variables.n_tf
                ) - build_variables.r_tf_inboard_in

            # Rounded resistive TF geometry
            else:
                build_variables.tfcth = (
                    tfcoil_variables.dr_tf_wp
                    + tfcoil_variables.casthi
                    + tfcoil_variables.thkcas
                )

        # Radial build to tfcoil middle [m]
        build_variables.r_tf_inboard_mid = (
            build_variables.r_tf_inboard_in + 0.5e0 * build_variables.tfcth
        )

        # Radial build to tfcoil plasma facing side [m]
        build_variables.r_tf_inboard_out = (
            build_variables.r_tf_inboard_in + build_variables.tfcth
        )

        # WP radial thickness [m]
        # Calculated only if not used as an iteration variable
        if not any(numerics.ixc[0 : numerics.nvar] == 140):

            # SC magnets
            if tfcoil_variables.i_tf_sup == 1:
                tfcoil_variables.dr_tf_wp = (
                    np.cos(np.pi / tfcoil_variables.n_tf)
                    * build_variables.r_tf_inboard_out
                    - build_variables.r_tf_inboard_in
                    - tfcoil_variables.casthi
                    - tfcoil_variables.thkcas
                )

            # Resistive magnets
            else:
                tfcoil_variables.dr_tf_wp = (
                    build_variables.tfcth
                    - tfcoil_variables.casthi
                    - tfcoil_variables.thkcas
                )

        # Radius of the centrepost at the top of the machine
        if physics_variables.itart == 1 and tfcoil_variables.i_tf_sup != 1:

            # build_variables.r_cp_top is set using the plasma shape
            if build_variables.i_r_cp_top == 0:
                build_variables.r_cp_top = (
                    physics_variables.rmajor
                    - physics_variables.rminor * physics_variables.triang
                    - (
                        build_variables.tftsgap
                        + build_variables.thshield_ib
                        + build_variables.shldith
                        + build_variables.vvblgap
                        + build_variables.blnkith
                        + build_variables.fwith
                        + 3.0e0 * build_variables.scrapli
                    )
                    + tfcoil_variables.drtop
                )

                # Notify user that build_variables.r_cp_top has been set to 1.01*build_variables.r_tf_inboard_out (lvl 2 error)
                if build_variables.r_cp_top < 1.01e0 * build_variables.r_tf_inboard_out:
                    error_handling.fdiags[0] = build_variables.r_cp_top
                    error_handling.fdiags[1] = build_variables.r_tf_inboard_out
                    error_handling.report_error(268)

                    # build_variables.r_cp_top correction
                    build_variables.r_cp_top = build_variables.r_tf_inboard_out * 1.01e0

                # Top and mid-plane TF coil CP radius ratio
                build_variables.f_r_cp = (
                    build_variables.r_cp_top / build_variables.r_tf_inboard_out
                )

            # User defined build_variables.r_cp_top
            elif build_variables.i_r_cp_top == 1:

                # Notify user that build_variables.r_cp_top has been set to 1.01*build_variables.r_tf_inboard_out (lvl 2 error)
                if build_variables.r_cp_top < 1.01e0 * build_variables.r_tf_inboard_out:
                    error_handling.fdiags[0] = build_variables.r_cp_top
                    error_handling.fdiags[1] = build_variables.r_tf_inboard_out
                    error_handling.report_error(268)

                    # build_variables.r_cp_top correction
                    build_variables.r_cp_top = build_variables.r_tf_inboard_out * 1.01e0

                # Top / mid-plane TF CP radius ratio
                build_variables.f_r_cp = (
                    build_variables.r_cp_top / build_variables.r_tf_inboard_out
                )

            # build_variables.r_cp_top set as a fraction of the outer TF midplane radius
            elif build_variables.i_r_cp_top == 2:
                build_variables.r_cp_top = (
                    build_variables.f_r_cp * build_variables.r_tf_inboard_out
                )

        else:  # End of physics_variables.itart == 1 .and. tfcoil_variables.i_tf_sup /= 1
            build_variables.r_cp_top = build_variables.r_tf_inboard_out

        if build_variables.i_r_cp_top != 0 and (
            build_variables.r_cp_top
            > physics_variables.rmajor
            - physics_variables.rminor * physics_variables.triang
            - (
                build_variables.tftsgap
                + build_variables.thshield_ib
                + build_variables.shldith
                + build_variables.vvblgap
                + build_variables.blnkith
                + build_variables.fwith
                + 3.0e0 * build_variables.scrapli
            )
            + tfcoil_variables.drtop
        ):

            error_handling.fdiags[0] = build_variables.r_cp_top
            error_handling.report_error(256)
        if build_variables.tf_in_cs == 1:
            #  Radial position of vacuum vessel [m]
            build_variables.r_vv_inboard_out = (
                build_variables.r_tf_inboard_out
                + build_variables.ohcth
                + build_variables.gapoh
                + build_variables.precomp
                + build_variables.tftsgap
                + build_variables.thshield_ib
                + build_variables.gapds
                + build_variables.d_vv_in
            )
        else:
            build_variables.r_vv_inboard_out = (
                build_variables.r_tf_inboard_out
                + build_variables.tftsgap
                + build_variables.thshield_ib
                + build_variables.gapds
                + build_variables.d_vv_in
            )
        # Radial position of the inner side of inboard neutronic shield [m]
        build_variables.r_sh_inboard_in = build_variables.r_vv_inboard_out

        # Radial position of the plasma facing side of inboard neutronic shield [m]
        build_variables.r_sh_inboard_out = (
            build_variables.r_sh_inboard_in + build_variables.shldith
        )

        #  Radial build to centre of plasma (should be equal to physics_variables.rmajor)
        build_variables.rbld = (
            build_variables.r_sh_inboard_out
            + build_variables.vvblgap
            + build_variables.blnkith
            + build_variables.fwith
            + build_variables.scrapli
            + physics_variables.rminor
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
        if tfcoil_variables.i_tf_sup != 1:
            build_variables.tfthko = build_variables.tfootfi * build_variables.tfcth
        else:
            build_variables.tfthko = build_variables.tfcth

        #  Radius to centre of outboard TF coil legs
        build_variables.r_tf_outboard_mid = (
            build_variables.rsldo
            + build_variables.vvblgap
            + build_variables.d_vv_out
            + build_variables.gapomin
            + build_variables.thshield_ob
            + build_variables.tftsgap
            + 0.5e0 * build_variables.tfthko
        )

        # TF coil horizontal build_variables.bore [m]
        build_variables.dr_tf_inner_bore = (
            build_variables.r_tf_outboard_mid - 0.5e0 * build_variables.tfthko
        ) - (build_variables.r_tf_inboard_mid - 0.5e0 * build_variables.tfcth)

        (
            tfcoil_variables.ripple,
            r_tf_outboard_midl,
            self.ripflag,
        ) = self.ripple_amplitude(
            tfcoil_variables.ripmax,
            build_variables.r_tf_outboard_mid,
        )

        #  If the tfcoil_variables.ripple is too large then move the outboard TF coil leg
        if r_tf_outboard_midl > build_variables.r_tf_outboard_mid:
            build_variables.r_tf_outboard_mid = r_tf_outboard_midl
            build_variables.gapsto = (
                build_variables.r_tf_outboard_mid
                - 0.5e0 * build_variables.tfthko
                - build_variables.d_vv_out
                - build_variables.rsldo
                - build_variables.thshield_ob
                - build_variables.tftsgap
                - build_variables.vvblgap
            )
            build_variables.dr_tf_inner_bore = (
                build_variables.r_tf_outboard_mid - 0.5e0 * build_variables.tfthko
            ) - (build_variables.r_tf_inboard_mid - 0.5e0 * build_variables.tfcth)
        else:
            build_variables.gapsto = build_variables.gapomin

        #  Call tfcoil_variables.ripple calculation again with new build_variables.r_tf_outboard_mid/build_variables.gapsto value
        #  call rippl(tfcoil_variables.ripmax,rmajor,rminor,r_tf_outboard_mid,n_tf,ripple,r_tf_outboard_midl)
        (
            tfcoil_variables.ripple,
            r_tf_outboard_midl,
            self.ripflag,
        ) = self.ripple_amplitude(
            tfcoil_variables.ripmax,
            build_variables.r_tf_outboard_mid,
        )

        #  Calculate first wall area
        #  Old calculation... includes a mysterious factor 0.875
        # fwarea = 0.875e0 *     #     ( 4.0e0*pi**2*sf*physics_variables.rmajor*(physics_variables.rminor+0.5e0*(build_variables.scrapli+build_variables.scraplo)) )

        #  Half-height of first wall (internal surface)
        hbot = (
            physics_variables.rminor * physics_variables.kappa
            + build_variables.vgap
            + divertor_variables.divfix
            - build_variables.blnktth
            - 0.5e0 * (build_variables.fwith + build_variables.fwoth)
        )
        if physics_variables.idivrt == 2:  # (i.e. physics_variables.i_single_null=0)
            htop = hbot
        else:
            htop = (
                physics_variables.rminor * physics_variables.kappa
                + build_variables.vgaptop
            )

        hfw = 0.5e0 * (htop + hbot)

        if (physics_variables.itart == 1) or (
            fwbs_variables.fwbsshape == 1
        ):  # D-shaped

            #  Major radius to outer edge of inboard section
            r1 = (
                physics_variables.rmajor
                - physics_variables.rminor
                - build_variables.scrapli
            )

            #  Horizontal distance between inside edges,
            #  i.e. outer radius of inboard part to inner radius of outboard part

            r2 = (
                physics_variables.rmajor
                + physics_variables.rminor
                + build_variables.scraplo
            ) - r1
            #  Calculate surface area, assuming 100% coverage
            # maths_library.eshellarea was not working across
            # the interface so has been reimplemented here
            # as a test

            (
                build_variables.fwareaib,
                build_variables.fwareaob,
                build_variables.fwarea,
            ) = maths_library.dshellarea(r1, r2, hfw)

        else:  # Cross-section is assumed to be defined by two ellipses

            #  Major radius to centre of inboard and outboard ellipses
            #  (coincident in radius with top of plasma)

            r1 = (
                physics_variables.rmajor
                - physics_variables.rminor * physics_variables.triang
            )

            #  Distance between r1 and outer edge of inboard section

            r2 = r1 - (
                physics_variables.rmajor
                - physics_variables.rminor
                - build_variables.scrapli
            )

            #  Distance between r1 and inner edge of outboard section

            r3 = (
                physics_variables.rmajor
                + physics_variables.rminor
                + build_variables.scraplo
            ) - r1

            #  Calculate surface area, assuming 100% coverage

            # maths_library.eshellarea was not working across
            # the interface so has been reimplemented here
            # as a test

            (
                build_variables.fwareaib,
                build_variables.fwareaob,
                build_variables.fwarea,
            ) = maths_library.eshellarea(r1, r2, r3, hfw)

        #  Apply area coverage factor

        if physics_variables.idivrt == 2:
            # Double null configuration
            build_variables.fwareaob = build_variables.fwareaob * (
                1.0e0 - 2.0e0 * fwbs_variables.fdiv - fwbs_variables.fhcd
            )
            build_variables.fwareaib = build_variables.fwareaib * (
                1.0e0 - 2.0e0 * fwbs_variables.fdiv - fwbs_variables.fhcd
            )
        else:
            # Single null configuration
            build_variables.fwareaob = build_variables.fwareaob * (
                1.0e0 - fwbs_variables.fdiv - fwbs_variables.fhcd
            )
            build_variables.fwareaib = build_variables.fwareaib * (
                1.0e0 - fwbs_variables.fdiv - fwbs_variables.fhcd
            )

        build_variables.fwarea = build_variables.fwareaib + build_variables.fwareaob

        if build_variables.fwareaob <= 0.0e0:
            error_handling.fdiags[0] = fwbs_variables.fdiv
            error_handling.fdiags[1] = fwbs_variables.fhcd
            error_handling.report_error(61)

        #

        if output:

            #  Print out device build

            po.oheadr(self.outfile, "Radial Build")

            if self.ripflag != 0:
                po.ocmmnt(
                    self.outfile,
                    "(Ripple result may not be accurate, as the fit was outside",
                )
                po.ocmmnt(self.outfile, " its range of applicability.)")
                po.oblnkl(self.outfile)
                error_handling.report_error(62)

                if self.ripflag == 1:
                    error_handling.fdiags[0] = (
                        tfcoil_variables.wwp1
                        * tfcoil_variables.n_tf
                        / physics_variables.rmajor
                    )
                    error_handling.report_error(141)
                elif self.ripflag == 2:
                    # Convert to integer as idiags is integer array
                    error_handling.idiags[0] = int(tfcoil_variables.n_tf)
                    error_handling.report_error(142)
                else:
                    error_handling.fdiags[0] = (
                        physics_variables.rmajor + physics_variables.rminor
                    ) / build_variables.r_tf_outboard_mid
                    error_handling.report_error(143)

            po.ovarin(
                self.outfile,
                "TF coil radial placement switch",
                "(tf_in_cs)",
                build_variables.tf_in_cs,
            )
            po.ovarrf(
                self.outfile,
                "Inboard build thickness (m)",
                "(inboard_build)",
                physics_variables.rmajor - physics_variables.rminor,
                "OP ",
            )

            if build_variables.tf_in_cs == 1:
                po.ocmmnt(
                    self.outfile,
                    ("\n (The stated machine bore size is just for the hollow space, "),
                )
                po.ocmmnt(
                    self.outfile,
                    (
                        "the true bore size used for calculations is bore + tfcth + gapoh)\n"
                    ),
                )
            if build_variables.tf_in_cs == 1 and tfcoil_variables.i_tf_bucking >= 2:
                po.ocmmnt(
                    self.outfile,
                    "(Bore hollow space has been filled with a solid metal cyclinder to act as wedge support)\n",
                )

            radius = 0.0e0
            po.obuild(self.outfile, "Device centreline", 0.0e0, radius)
            if build_variables.tf_in_cs == 1 and tfcoil_variables.i_tf_bucking >= 2:
                radius = (
                    radius
                    + build_variables.bore
                    - build_variables.tfcth
                    - build_variables.gapoh
                )
                po.obuild(
                    self.outfile,
                    "Machine bore wedge support",
                    build_variables.bore
                    - build_variables.tfcth
                    - build_variables.gapoh,
                    radius,
                    "(bore)",
                )
                po.ovarre(
                    self.mfile,
                    "Machine bore wedge support cylinder (m)",
                    "(bore)",
                    build_variables.bore
                    - build_variables.tfcth
                    - build_variables.gapoh,
                )
            elif build_variables.tf_in_cs == 1 and tfcoil_variables.i_tf_bucking < 2:
                radius = (
                    radius
                    + build_variables.bore
                    - build_variables.tfcth
                    - build_variables.gapoh
                )
                po.obuild(
                    self.outfile,
                    "Machine bore hole",
                    build_variables.bore
                    - build_variables.tfcth
                    - build_variables.gapoh,
                    radius,
                    "(bore)",
                )
                po.ovarre(
                    self.mfile,
                    "Machine bore hole (m)",
                    "(bore)",
                    build_variables.bore
                    - build_variables.tfcth
                    - build_variables.gapoh,
                )
            else:
                radius = radius + build_variables.bore
                po.obuild(
                    self.outfile,
                    "Machine bore",
                    build_variables.bore,
                    radius,
                    "(bore)",
                )
                po.ovarre(
                    self.mfile,
                    "Machine bore (m)",
                    "(bore)",
                    build_variables.bore,
                )
            if build_variables.tf_in_cs == 1:
                radius += build_variables.tfcth
                po.obuild(
                    self.outfile,
                    "TF coil inboard leg,in bore",
                    build_variables.tfcth,
                    radius,
                    "(tfcth)",
                )
                po.ovarre(
                    self.mfile,
                    "TF coil inboard leg (m)",
                    "(tfcth)",
                    build_variables.tfcth,
                )
            if build_variables.tf_in_cs == 1:
                radius += build_variables.gapoh
                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.gapoh,
                    radius,
                    "(gapoh)",
                )
                po.ovarre(
                    self.mfile,
                    "CS precompresion to TF coil radial gap (m)",
                    "(gapoh)",
                    build_variables.gapoh,
                )

            radius = radius + build_variables.ohcth
            po.obuild(
                self.outfile,
                "Central solenoid",
                build_variables.ohcth,
                radius,
                "(ohcth)",
            )

            po.ovarre(
                self.mfile,
                "CS radial thickness (m)",
                "(ohcth)",
                build_variables.ohcth,
            )

            radius = radius + build_variables.precomp
            po.obuild(
                self.outfile,
                "CS precompression",
                build_variables.precomp,
                radius,
                "(precomp)",
            )
            po.ovarre(
                self.mfile,
                "CS precompression (m)",
                "(precomp)",
                build_variables.precomp,
            )
            if build_variables.tf_in_cs == 0:
                radius = radius + build_variables.gapoh
                po.obuild(
                    self.outfile,
                    "Gap",
                    build_variables.gapoh,
                    radius,
                    "(gapoh)",
                )
                po.ovarre(
                    self.mfile,
                    "CS precompresion to TF coil radial gap (m)",
                    "(gapoh)",
                    build_variables.gapoh,
                )
            if build_variables.tf_in_cs == 0:
                radius = radius + build_variables.tfcth
                po.obuild(
                    self.outfile,
                    "TF coil inboard leg",
                    build_variables.tfcth,
                    radius,
                    "(tfcth)",
                )
                po.ovarre(
                    self.mfile,
                    "TF coil inboard leg (m)",
                    "(tfcth)",
                    build_variables.tfcth,
                )

            radius = radius + build_variables.tftsgap
            po.obuild(
                self.outfile,
                "Gap",
                build_variables.tftsgap,
                radius,
                "(tftsgap)",
            )
            po.ovarre(
                self.mfile,
                "TF coil inboard leg insulation gap (m)",
                "(tftsgap)",
                build_variables.tftsgap,
            )

            radius = radius + build_variables.thshield_ib
            po.obuild(
                self.outfile,
                "Thermal shield, inboard",
                build_variables.thshield_ib,
                radius,
                "(thshield_ib)",
            )
            po.ovarre(
                self.mfile,
                "Thermal shield, inboard (m)",
                "(thshield_ib)",
                build_variables.thshield_ib,
            )

            radius = radius + build_variables.gapds
            po.obuild(
                self.outfile,
                "Gap",
                build_variables.gapds,
                radius,
                "(gapds)",
            )
            po.ovarre(
                self.mfile,
                "thermal shield to vessel radial gap (m)",
                "(gapds)",
                build_variables.gapds,
            )

            radius = radius + build_variables.d_vv_in + build_variables.shldith
            po.obuild(
                self.outfile,
                "Vacuum vessel (and shielding)",
                build_variables.d_vv_in + build_variables.shldith,
                radius,
                "(d_vv_in + shldith)",
            )
            po.ovarre(
                self.mfile,
                "Inboard vacuum vessel radial thickness (m)",
                "(d_vv_in)",
                build_variables.d_vv_in,
            )
            po.ovarre(
                self.mfile,
                "Inner radiation shield radial thickness (m)",
                "(shldith)",
                build_variables.shldith,
            )

            radius = radius + build_variables.vvblgap
            po.obuild(
                self.outfile,
                "Gap",
                build_variables.vvblgap,
                radius,
                "(vvblgap)",
            )
            po.ovarre(
                self.mfile,
                "Gap (m)",
                "(vvblgap)",
                build_variables.vvblgap,
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
                self.mfile,
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
                self.mfile,
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
                self.mfile,
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
                self.mfile,
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
                self.mfile,
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
                self.mfile,
                "Outboard blanket radial thickness (m)",
                "(blnkoth)",
                build_variables.blnkoth,
            )

            radius = radius + build_variables.vvblgap
            po.obuild(
                self.outfile,
                "Gap",
                build_variables.vvblgap,
                radius,
                "(vvblgap)",
            )

            radius = radius + build_variables.d_vv_out + build_variables.shldoth
            po.obuild(
                self.outfile,
                "Vacuum vessel (and shielding)",
                build_variables.d_vv_out + build_variables.shldoth,
                radius,
                "(d_vv_out+shldoth)",
            )
            po.ovarre(
                self.mfile,
                "Outer radiation shield radial thickness (m)",
                "(shldoth)",
                build_variables.shldoth,
            )
            po.ovarre(
                self.mfile,
                "Outboard vacuum vessel radial thickness (m)",
                "(d_vv_out)",
                build_variables.d_vv_out,
            )

            radius = radius + build_variables.gapsto
            po.obuild(
                self.outfile,
                "Gap",
                build_variables.gapsto,
                radius,
                "(gapsto)",
            )
            po.ovarre(
                self.mfile,
                "Vessel to TF radial gap (m)",
                "(gapsto)",
                build_variables.gapsto,
            )

            radius = radius + build_variables.thshield_ob
            po.obuild(
                self.outfile,
                "Thermal shield, outboard",
                build_variables.thshield_ob,
                radius,
                "(thshield_ob)",
            )
            po.ovarre(
                self.mfile,
                "Thermal shield, outboard (m)",
                "(thshield_ob)",
                build_variables.thshield_ob,
            )

            radius = radius + build_variables.tftsgap
            po.obuild(
                self.outfile,
                "Gap",
                build_variables.tftsgap,
                radius,
                "(tftsgap)",
            )
            po.ovarre(
                self.mfile,
                "Gap (m)",
                "(tftsgap)",
                build_variables.tftsgap,
            )

            radius = radius + build_variables.tfthko
            po.obuild(
                self.outfile,
                "TF coil outboard leg",
                build_variables.tfthko,
                radius,
                "(tfthko)",
            )
            po.ovarre(
                self.mfile,
                "TF coil outboard leg radial thickness (m)",
                "(tfthko)",
                build_variables.tfthko,
            )

            po.ovarre(
                self.mfile,
                "External cryostat thickness (excludes structure) (m)",
                "(ddwex)",
                build_variables.ddwex,
            )

            if (current_drive_variables.iefrf in [5, 8]) or (
                current_drive_variables.iefrffix in [5, 8]
            ):
                po.ovarre(
                    self.mfile,
                    "Width of neutral beam duct where it passes between the TF coils (m)",
                    "(beamwd)",
                    current_drive_variables.beamwd,
                )
