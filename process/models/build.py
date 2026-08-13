"""Module containing routines for build calculations"""

import logging
from enum import IntEnum

import numpy as np

from process.core import constants
from process.core import process_output as po
from process.core.model import Model
from process.data_structure.build_variables import (
    CSPrecompressionConfiguration,
    TFCSRadialConfiguration,
)
from process.data_structure.physics_variables import DivertorNumberModels
from process.models.physics.current_drive import (
    CurrentDriveMethodType,
    CurrentDriveModel,
)
from process.models.tfcoil.base import TFCoilShapeModel
from process.models.tfcoil.superconducting import SuperconductingTFWPShapeType

logger = logging.getLogger(__name__)


class FwBlktVVShape(IntEnum):
    """Enum for first wall, blanket, and vacuum vessel shape options."""

    D_SHAPED = 1
    ELLIPTICAL_SHAPED = 2


class Build(Model):
    """Routines for build calculations"""

    def __init__(self):
        self.outfile = constants.NOUT
        self.mfile = constants.MFILE

    def output(self):
        """Output the build information"""
        # Radial build
        self.calculate_radial_build(output=True)

        # Vertical build
        self.calculate_vertical_build(output=True)

    def run(self):
        """Run the build routines"""
        self.calculate_radial_build(output=False)
        self.calculate_vertical_build(output=False)

        (
            self.data.current_drive.radius_beam_tangency,
            self.data.current_drive.radius_beam_tangency_max,
        ) = self.calculate_beam_port_size(
            f_radius_beam_tangency_rmajor=self.data.current_drive.f_radius_beam_tangency_rmajor,
            rmajor=self.data.physics.rmajor,
            n_tf_coils=self.data.tfcoil.n_tf_coils,
            dx_tf_inboard_out_toroidal=self.data.tfcoil.dx_tf_inboard_out_toroidal,
            dr_tf_outboard=self.data.build.dr_tf_outboard,
            r_tf_outboard_mid=self.data.build.r_tf_outboard_mid,
            dx_beam_duct=self.data.current_drive.dx_beam_duct,
            dx_beam_shield=self.data.current_drive.dx_beam_shield,
        )

    @staticmethod
    def calculate_beam_port_size(
        f_radius_beam_tangency_rmajor: float,
        rmajor: float,
        n_tf_coils: int,
        dx_tf_inboard_out_toroidal: float,
        dr_tf_outboard: float,
        r_tf_outboard_mid: float,
        dx_beam_duct: float,
        dx_beam_shield: float,
    ) -> tuple[float, float]:
        """Calculates the maximum possible tangency radius for adequate beam access.

        Parameters
        ----------
        f_radius_beam_tangency_rmajor : float
            Fraction of rmajor for beam tangency
        rmajor : float
            Major radius
        n_tf_coils : int
            Number of TF coils
        dx_tf_inboard_out_toroidal : float
            Toroidal width of outboard TF coil
        dr_tf_outboard : float
            Radial thickness of outboard TF coil leg
        r_tf_outboard_mid : float
            Major radius of centre of outboard TF coil
        dx_beam_duct : float
            Width of beam duct
        dx_beam_shield : float
            Shielding width on both sides of beam duct

        Returns
        -------
        tuple[float, float]
            Tuple containing (radius_beam_tangency, radius_beam_tangency_max)
        """
        # Have kept the single letter variable names to match the original code and
        # documentation diagram.
        radius_beam_tangency = f_radius_beam_tangency_rmajor * rmajor

        omega = 2.0 * np.pi / n_tf_coils

        a = 0.5e0 * dx_tf_inboard_out_toroidal

        if np.isinf(a):
            logger.error("a is inf. Kludging to 1e10.")
            a = 1e10

        b = dr_tf_outboard
        if np.isinf(b):
            logger.error("b is inf. Kludging to 1e10.")
            b = 1e10

        c = dx_beam_duct + 2.0e0 * dx_beam_shield

        d = r_tf_outboard_mid - 0.5e0 * b
        if np.isinf(d):
            logger.error("d is inf. Kludging to 1e10.")
            d = 1e10

        e = np.sqrt(a**2 + (d + b) ** 2)
        f = np.sqrt(a**2 + d**2)

        theta = omega - np.arctan(a / d)
        phi = theta - np.arcsin(a / e)

        g = np.sqrt(e**2 + f**2 - 2.0e0 * e * f * np.cos(phi))

        if g > c:
            h = np.sqrt(g**2 - c**2)
            alpha = np.arctan(h / c)
            eps = np.arcsin(e * np.sin(phi) / g) - alpha
            radius_beam_tangency_max = f * np.cos(eps) - 0.5e0 * c
        else:
            logger.error(
                "Max beam tangency radius set =0 temporarily; "
                "change dx_beam_duct. %s %s",
                g,
                c,
            )
            radius_beam_tangency_max = 0.0e0

        return radius_beam_tangency, radius_beam_tangency_max

    def _vertical_build_out(self, i_single_null: DivertorNumberModels):
        bld = self.data.build
        sn = i_single_null == i_single_null.SINGLE_NULL

        po.ocmmnt(self.outfile, f"{'Single' if sn else 'Double'} null case")
        bld.dz_vv_upper = 0.5 * (bld.dz_vv_upper + bld.dz_vv_lower)
        bld.dz_fw_upper = 0.5 * (bld.dr_fw_inboard + bld.dr_fw_outboard)

        vbuild = (
            self.data.buildings.dz_tf_cryostat
            + bld.dr_tf_inboard
            + bld.dr_tf_shld_gap
            + bld.dz_shld_thermal
            + bld.dz_shld_vv_gap
            + bld.dz_vv_upper
            + (bld.dr_shld_blkt_gap if sn else 0)
            + bld.dz_shld_upper
            + (bld.dz_blkt_upper if sn else 0)
            + (
                (bld.dz_fw_upper + bld.dz_fw_plasma_gap)
                if sn
                else (
                    self.data.divertor.dz_divertor + self.data.build.dz_xpoint_divertor
                )
            )
            + bld.z_plasma_xpoint_upper
        )

        # To calculate vertical offset between TF coil centre and plasma centre
        vbuile1 = vbuild

        dz_fw_upper = 0.5e0 * (bld.dr_fw_inboard + bld.dr_fw_outboard)

        # Top of TF coil
        tf_top = vbuild - self.data.buildings.dz_tf_cryostat
        vbuild = self.write_obuild(
            vbuild,
            [
                (
                    "Cryostat roof structure*",
                    self.data.buildings.dz_tf_cryostat,
                    "(dz_tf_cryostat)",
                ),
                ("TF coil", bld.dr_tf_inboard, "(dr_tf_inboard)"),
                ("Gap", bld.dr_tf_shld_gap, "(dr_tf_shld_gap)"),
                ("Thermal shield, vertical", bld.dz_shld_thermal, "(dz_shld_thermal)"),
                ("Gap", bld.dz_shld_vv_gap, "(dz_shld_vv_gap)"),
                (
                    "Vacuum vessel (and shielding)",
                    bld.dz_vv_upper + bld.dz_shld_upper,
                    "(dz_vv_upper+dz_shld_upper)",
                ),
                ("Gap", bld.dr_shld_blkt_gap, "(dr_shld_blkt_gap)"),
                ("Top blanket", bld.dz_blkt_upper, "(dz_blkt_upper)"),
                ("Top first wall", dz_fw_upper, "(dz_fw_upper)")
                if sn
                else (
                    "Divertor structure",
                    self.data.divertor.dz_divertor,
                    "(dz_divertor)",
                ),
                ("Top scrape-off", bld.dz_fw_plasma_gap, "(dz_fw_plasma_gap)"),
                (
                    "Plasma upper X-point height (m)",
                    bld.z_plasma_xpoint_upper,
                    "(z_plasma_xpoint_upper)",
                ),
            ],
        )

        for desc, name, val in [
            (
                "Cryostat roof structure*",
                "(dz_tf_cryostat)",
                self.data.buildings.dz_tf_cryostat,
            ),
            ("Thermal shield, vertical (m)", "(dz_shld_thermal)", bld.dz_shld_thermal),
            (
                "Vessel - TF coil vertical gap (m)",
                "(dz_shld_vv_gap)",
                bld.dz_shld_vv_gap,
            ),
            (
                "Topside vacuum vessel radial thickness (m)",
                "(dz_vv_upper)",
                bld.dz_vv_upper,
            ),
            ("Top radiation shield thickness (m)", "(dz_shld_upper)", bld.dz_shld_upper),
            ("Top blanket vertical thickness (m)", "(dz_blkt_upper)", bld.dz_blkt_upper),
            ("Top first wall vertical thickness (m)", "(dz_fw_upper)", dz_fw_upper)
            if sn
            else (
                "Divertor structure vertical thickness (m)",
                "(dz_divertor)",
                self.data.divertor.dz_divertor,
            ),
            (
                "Top scrape-off vertical thickness (m)",
                "(dz_fw_plasma_gap)",
                bld.dz_fw_plasma_gap,
            ),
            (
                "Plasma upper X-point height (m)",
                "(z_plasma_xpoint_upper)",
                bld.z_plasma_xpoint_upper,
            ),
        ]:
            po.ovarre(self.mfile, desc, name, val)

        po.obuild(self.outfile, "Midplane", 0.0e0, vbuild)

        vbuild = self.write_obuild(
            vbuild,
            [
                (
                    "Plasma lower X-point height (m)",
                    bld.z_plasma_xpoint_lower,
                    "(z_plasma_xpoint_lower)",
                ),
                ("Lower scrape-off", bld.dz_xpoint_divertor, "(dz_xpoint_divertor)"),
                ("Divertor structure", self.data.divertor.dz_divertor, "(dz_divertor)"),
                (
                    "Vacuum vessel (and shielding)",
                    bld.dz_vv_lower + bld.dz_shld_lower,
                    "(dz_vv_lower+dz_shld_lower)",
                ),
                ("Gap", bld.dz_shld_vv_gap, "(dz_shld_vv_gap)"),
                ("Thermal shield, vertical", bld.dz_shld_thermal, "(dz_shld_thermal)"),
                ("Gap", bld.dr_tf_shld_gap, "(dr_tf_shld_gap)"),
                ("TF coil", bld.dr_tf_inboard, "(dr_tf_inboard)"),
                (
                    "Cryostat floor structure**",
                    self.data.buildings.dz_tf_cryostat,
                    "(dz_tf_cryostat)",
                ),
            ],
            before=True,
        )
        for desc, name, val in [
            (
                "Plasma lower X-point height (m)",
                "(z_plasma_xpoint_lower)",
                bld.z_plasma_xpoint_lower,
            ),
            (
                "Bottom scrape-off vertical thickness (m)",
                "(dz_xpoint_divertor)",
                bld.dz_xpoint_divertor,
            ),
            (
                "Divertor structure vertical thickness (m)",
                "(dz_divertor)",
                self.data.divertor.dz_divertor,
            ),
            (
                "Bottom radiation shield thickness (m)",
                "(dz_shld_lower)",
                bld.dz_shld_lower,
            ),
            (
                "Underside vacuum vessel radial thickness (m)",
                "(dz_vv_lower)",
                bld.dz_vv_lower,
            ),
        ]:
            po.ovarre(self.mfile, desc, name, val)

        # Total height of TF coil
        tf_height = tf_top - vbuild + self.data.buildings.dz_tf_cryostat
        # Inner vertical dimension of TF coil
        bld.dh_tf_inner_bore = tf_height - 2 * bld.dr_tf_inboard

        # To calculate vertical offset between TF coil centre and plasma centre
        bld.dz_tf_plasma_centre_offset = (vbuile1 + vbuild) / 2.0e0

        # end of Single null case

    def write_obuild(self, vbuild, entry: tuple | list[tuple], *, before=False):
        """Write obuild entry"""
        if not isinstance(entry, list):
            entry = [entry]
        for desc, var, name in entry:
            if before:
                vbuild -= var
            po.obuild(self.outfile, desc, var, vbuild, name)
            if not before:
                vbuild -= var

        return vbuild

    def calculate_vertical_build(self, output: bool):
        """Determines the vertical build of the machine.

        This method calculates various parameters related
        to the vertical build of the machine,
        such as thicknesses, radii, and areas.
        Results can be outputted with the `output` flag.

        Parameters
        ----------
        output : bool
            Flag indicating whether to output results
        """
        # Set the X-point heights for the top and bottom of the plasma
        # Assumes top-down plasma symmetry
        self.data.build.z_plasma_xpoint_upper = (
            self.data.physics.rminor * self.data.physics.kappa
        )
        self.data.build.z_plasma_xpoint_lower = (
            self.data.physics.rminor * self.data.physics.kappa
        )
        i_single_null = DivertorNumberModels(self.data.physics.i_single_null)

        if output:
            po.oheadr(self.outfile, "Vertical Build")

            po.ovarre(
                self.mfile,
                "Divertor null switch",
                "(i_single_null)",
                self.data.physics.i_single_null,
            )

            self._vertical_build_out(i_single_null)

            po.ovarre(
                self.mfile,
                "Ratio of Central Solenoid height to TF coil internal height",
                "(f_z_cs_tf_internal)",
                self.data.pf_coil.f_z_cs_tf_internal,
            )
            po.ocmmnt(
                self.outfile,
                "\n*Cryostat roof allowance includes uppermost PF coil and "
                "outer thermal shield.\n*Cryostat floor allowance "
                "includes lowermost PF coil, outer thermal shield and gravity support.",
            )

        # Output the cdivertor geometry
        divht = self.divgeom(output)
        # Issue #481 Remove self.data.build.vgaptf
        if self.data.build.dz_xpoint_divertor < 0.00001e0:
            self.data.build.dz_xpoint_divertor = divht

        # If self.data.build.dz_xpoint_divertor /= 0 use the value set by the user.

        # Height to inside edge of TF coil. TF coils are assumed to be symmetrical.
        # Therefore this applies to single and double null cases.
        self.data.build.z_tf_inside_half = (
            self.data.build.z_plasma_xpoint_upper
            + self.data.build.dz_xpoint_divertor
            + self.data.divertor.dz_divertor
            + self.data.build.dz_shld_lower
            + self.data.build.dz_vv_lower
            + self.data.build.dz_shld_vv_gap
            + self.data.build.dz_shld_thermal
            + self.data.build.dr_tf_shld_gap
        )

        #  Vertical locations of divertor coils
        if i_single_null == DivertorNumberModels.DOUBLE_NULL:
            self.data.build.z_tf_top = (
                self.data.build.z_tf_inside_half + self.data.build.dr_tf_inboard
            )
            self.data.build.dz_tf_upper_lower_midplane = 0.0e0
        else:
            self.data.build.z_tf_top = (
                self.data.build.dr_tf_inboard
                + self.data.build.dr_tf_shld_gap
                + self.data.build.dz_shld_thermal
                + self.data.build.dz_shld_vv_gap
                + self.data.build.dz_vv_upper
                + self.data.build.dz_shld_upper
                + self.data.build.dr_shld_blkt_gap
                + self.data.build.dz_blkt_upper
                + 0.5e0
                * (self.data.build.dr_fw_inboard + self.data.build.dr_fw_outboard)
                + self.data.build.dz_fw_plasma_gap
                + self.data.build.z_plasma_xpoint_upper
            )
            self.data.build.dz_tf_upper_lower_midplane = self.data.build.z_tf_top - (
                self.data.build.z_tf_inside_half + self.data.build.dr_tf_inboard
            )

    def divgeom(self, output: bool):
        """Divertor geometry calculation

        This subroutine determines the divertor geometry.
        The inboard (i) and outboard (o) plasma surfaces
        are approximated by arcs, and followed past the X-point to
        determine the maximum height.
        TART option: Peng SOFT paper

        Parameters
        ----------
        output: bool

        Returns
        -------
        divht:
            divertor height (m)
        """
        if self.data.physics.itart == 1:
            return 1.75e0 * self.data.physics.rminor
        #  Conventional tokamak divertor model
        #  options for separate upper and lower self.data.physics.triangularity

        kap = self.data.physics.kappa
        triu = tril = self.data.physics.triang

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
        def rc(trilio):
            return 0.5 * np.sqrt(
                (self.data.physics.rminor**2 * (trilio**2 + kap**2) ** 2)
                / ((tril + 1.0e0) ** 2)
            )

        rco = rc(tril + 1)
        rci = rc(tril - 1)

        # Find angles between vertical and legs
        def theta(trilio, rcio):
            return np.arcsin(
                1.0e0 - (self.data.physics.rminor * (1.0e0 - trilio)) / rcio
            )

        # Inboard arc angle = outboard leg angle
        thetai = theta(1 - tril, rci)
        # Outboard arc angle = inboard leg angle
        thetao = theta(1 + tril, rco)

        #  Position of lower x-pt
        rxpt = self.data.physics.rmajor - tril * self.data.physics.rminor
        zxpt = -1.0e0 * kap * self.data.physics.rminor

        # Position of outer strike point
        self.data.build.rspo = rxpt + self.data.build.plsepo * np.cos(thetao)
        zspo = zxpt - self.data.build.plsepo * np.sin(thetao)

        # Position of inner plate ends
        def plate_cos(l_div_plate, theta, beta):
            return (l_div_plate / 2) * np.cos(theta + beta)

        def plate_sin(l_div_plate, theta, beta):
            return (l_div_plate / 2) * np.sin(theta + beta)

        inner_plte_cos = plate_cos(
            self.data.build.plleni, thetai, self.data.divertor.betai
        )
        inner_plte_sin = plate_sin(
            self.data.build.plleni, thetai, self.data.divertor.betai
        )
        outer_plte_cos = plate_cos(
            self.data.build.plleno, thetao, self.data.divertor.betao
        )
        outer_plte_sin = plate_sin(
            self.data.build.plleno, thetao, self.data.divertor.betao
        )

        # Position of inner strike point
        # rspi = rxpt - self.data.build.plsepi*cos(alphad)
        # zspi = zxpt - self.data.build.plsepi*sin(alphad)
        rspi = rxpt - self.data.build.plsepi * np.cos(thetai)
        zspi = zxpt - self.data.build.plsepi * np.sin(thetai)
        zplti = zspi + inner_plte_sin
        zplbi = zspi - inner_plte_sin
        zplto = zspo + outer_plte_sin
        zplbo = zspo - outer_plte_sin

        divht = max(zplti, zplto) - min(zplbo, zplbi)

        if output:
            self.divertor_geom_output(
                kap,
                divht,
                tril,
                triu,
                thetai,
                thetao,
                rco,
                rci,
                rxpt,
                zxpt,
                rspi,
                zspi,
                zspo,
                # Position of inner strike points
                rplti=rspi + inner_plte_cos,
                rplbi=rspi - inner_plte_cos,
                zplti=zplti,
                zplbi=zplbi,
                # Position of outer plate ends
                rplto=self.data.build.rspo - outer_plte_cos,
                rplbo=self.data.build.rspo + outer_plte_cos,
                zplto=zplto,
                zplbo=zplbo,
            )

        return divht

    def divertor_geom_output(
        self,
        kap,
        divht,
        tril,
        triu,
        thetai,
        thetao,
        rco,
        rci,
        rxpt,
        zxpt,
        rspi,
        zspi,
        zspo,
        rplti,
        rplbi,
        zplti,
        zplbi,
        rplto,
        rplbo,
        zplto,
        zplbo,
    ):
        """Divertor geometry output"""
        po.oheadr(self.outfile, "Divertor build and plasma position")

        if self.data.divertor.n_divertors == 1:
            po.ocmmnt(self.outfile, "Divertor Configuration = Single Null Divertor")
            po.oblnkl(self.outfile)

            for desc, name, var in [
                (
                    "Plasma top position, radial (m)",
                    "(ptop_radial)",
                    self.data.physics.rmajor - triu * self.data.physics.rminor,
                ),
                (
                    "Plasma top position, vertical (m)",
                    "(ptop_vertical)",
                    kap * self.data.physics.rminor,
                ),
                (
                    "Plasma geometric centre, radial (m)",
                    "(rmajor.)",
                    self.data.physics.rmajor,
                ),
                ("Plasma geometric centre, vertical (m)", "(0.0)", 0.0e0),
                ("Plasma lower triangularity", "(tril)", tril),
                ("Plasma elongation", "(kappa.)", kap),
                (
                    "TF coil vertical offset (m)",
                    "(dz_tf_plasma_centre_offset)",
                    self.data.build.dz_tf_plasma_centre_offset,
                ),
                ("Plasma outer arc radius of curvature (m)", "(rco)", rco),
                ("Plasma inner arc radius of curvature (m)", "(rci)", rci),
                ("Plasma lower X-pt, radial (m)", "(rxpt)", rxpt),
                ("Plasma lower X-pt, vertical (m)", "(zxpt)", zxpt),
                (
                    "Poloidal plane angle between vertical and inner leg (rad)",
                    "(thetai)",
                    thetai,
                ),
                (
                    "Poloidal plane angle between vertical and outer leg (rad)",
                    "(thetao)",
                    thetao,
                ),
                (
                    "Poloidal plane angle between inner leg and plate (rad)",
                    "(betai)",
                    self.data.divertor.betai,
                ),
                (
                    "Poloidal plane angle between outer leg and plate (rad)",
                    "(betao)",
                    self.data.divertor.betao,
                ),
                (
                    "Inner divertor leg poloidal length (m)",
                    "(plsepi)",
                    self.data.build.plsepi,
                ),
                (
                    "Outer divertor leg poloidal length (m)",
                    "(plsepo)",
                    self.data.build.plsepo,
                ),
                ("Inner divertor plate length (m)", "(plleni)", self.data.build.plleni),
                ("Outer divertor plate length (m)", "(plleno)", self.data.build.plleno),
                ("Inner strike point, radial (m)", "(rspi)", rspi),
                ("Inner strike point, vertical (m)", "(zspi)", zspi),
                ("Inner plate top, radial (m)", "(rplti)", rplti),
                ("Inner plate top, vertical (m)", "(zplti)", zplti),
                ("Inner plate bottom, radial (m)", "(rplbi)", rplbi),
                ("Inner plate bottom, vertical (m)", "(zplbi)", zplbi),
                ("Outer strike point, radial (m)", "(rspo)", self.data.build.rspo),
                ("Outer strike point, vertical (m)", "(zspo)", zspo),
                ("Outer plate top, radial (m)", "(rplto)", rplto),
                ("Outer plate top, vertical (m)", "(zplto)", zplto),
                ("Outer plate bottom, radial (m)", "(rplbo)", rplbo),
                ("Outer plate bottom, vertical (m)", "(zplbo)", zplbo),
                ("Calculated maximum divertor height (m)", "(divht)", divht),
            ]:
                po.ovarre(self.outfile, desc, name, var, "OP ")

        elif self.data.divertor.n_divertors == 2:
            po.ocmmnt(self.outfile, "Divertor Configuration = Double Null Divertor")
            po.oblnkl(self.outfile)
            # Assume upper and lower divertors geometries are symmetric.
            for desc, name, var in [
                (
                    "Plasma top position, radial (m)",
                    "(ptop_radial)",
                    self.data.physics.rmajor - triu * self.data.physics.rminor,
                ),
                (
                    "Plasma top position, vertical (m)",
                    "(ptop_vertical)",
                    kap * self.data.physics.rminor,
                ),
                (
                    "Plasma geometric centre, radial (m)",
                    "(rmajor.)",
                    self.data.physics.rmajor,
                ),
                ("Plasma geometric centre, vertical (m)", "(0.0)", 0.0e0),
                ("Plasma data.physics.triangularity", "(tril)", tril),
                ("Plasma elongation", "(kappa.)", kap),
                (
                    "TF coil vertical offset (m)",
                    "(dz_tf_plasma_centre_offset)",
                    self.data.build.dz_tf_plasma_centre_offset,
                ),
                ("Plasma upper X-pt, radial (m)", "(rxpt)", rxpt),
                ("Plasma upper X-pt, vertical (m)", "(-zxpt)", -zxpt),
                ("Plasma outer arc radius of curvature (m)", "(rco)", rco),
                ("Plasma inner arc radius of curvature (m)", "(rci)", rci),
                ("Plasma lower X-pt, radial (m)", "(rxpt)", rxpt),
                ("Plasma lower X-pt, vertical (m)", "(zxpt)", zxpt),
                (
                    "Poloidal plane angle between vertical and inner leg (rad)",
                    "(thetai)",
                    thetai,
                ),
                (
                    "Poloidal plane angle between vertical and outer leg (rad)",
                    "(thetao)",
                    thetao,
                ),
                (
                    "Poloidal plane angle between inner leg and plate (rad)",
                    "(betai)",
                    self.data.divertor.betai,
                ),
                (
                    "Poloidal plane angle between outer leg and plate (rad)",
                    "(betao)",
                    self.data.divertor.betao,
                ),
                (
                    "Inner divertor leg poloidal length (m)",
                    "(plsepi)",
                    self.data.build.plsepi,
                ),
                (
                    "Outer divertor leg poloidal length (m)",
                    "(plsepo)",
                    self.data.build.plsepo,
                ),
                ("Inner divertor plate length (m)", "(plleni)", self.data.build.plleni),
                ("Outer divertor plate length (m)", "(plleno)", self.data.build.plleno),
                ("Upper inner strike point, radial (m)", "(rspi)", rspi),
                ("Upper inner strike point, vertical (m)", "(-zspi)", -zspi),
                ("Upper inner plate top, radial (m)", "(rplti)", rplti),
                ("Upper inner plate top, vertical (m)", "(-zplti)", -zplti),
                ("Upper inner plate bottom, radial (m)", "(rplbi)", rplbi),
                ("Upper inner plate bottom, vertical (m)", "(-zplbi)", -zplbi),
                ("Upper outer strike point, radial (m)", "(rspo)", self.data.build.rspo),
                ("Upper outer strike point, vertical (m)", "(-zspo)", -zspo),
                ("Upper outer plate top, radial (m)", "(rplto)", rplto),
                ("Upper outer plate top, vertical (m)", "(-zplto)", -zplto),
                ("Upper outer plate bottom, radial (m)", "(rplbo)", rplbo),
                ("Upper outer plate bottom, vertical (m)", "(-zplbo)", -zplbo),
                ("Lower inner strike point, radial (m)", "(rspi)", rspi),
                ("Lower inner strike point, vertical (m)", "(zspi)", zspi),
                ("Lower inner plate top, radial (m)", "(rplti)", rplti),
                ("Lower inner plate top, vertical (m)", "(zplti)", zplti),
                ("Lower inner plate bottom, radial (m)", "(rplbi)", rplbi),
                ("Lower inner plate bottom, vertical (m)", "(zplbi)", zplbi),
                ("Lower outer strike point, radial (m)", "(rspo)", self.data.build.rspo),
                ("Lower outer strike point, vertical (m)", "(zspo)", zspo),
                ("Lower outer plate top, radial (m)", "(rplto)", rplto),
                ("Lower outer plate top, vertical (m)", "(zplto)", zplto),
                ("Lower outer plate bottom, radial (m)", "(rplbo)", rplbo),
                ("Lower outer plate bottom, vertical (m)", "(zplbo)", zplbo),
                ("Calculated maximum divertor height (m)", "(divht)", divht),
            ]:
                po.ovarre(self.outfile, desc, name, var, "OP ")

        else:
            po.ocmmnt(
                self.outfile,
                "ERROR: null value not supported, check i_single_null value.",
            )

    @staticmethod
    def plasma_outboard_edge_toroidal_ripple(
        ripple_b_tf_plasma_edge_max: float,
        r_tf_outboard_mid: float,
        n_tf_coils: int,
        rmajor: float,
        rminor: float,
        r_tf_wp_inboard_inner,
        r_tf_wp_inboard_centre: float,
        r_tf_wp_inboard_outer: float,
        dx_tf_wp_primary_toroidal: float,
        i_tf_shape: int,
        i_tf_sup: int,
        dx_tf_wp_insulation: float,
        dx_tf_wp_insertion_gap: float,
        i_tf_wp_geom: int,
    ) -> tuple[float, float, int]:
        """Plasma outboard toroidal field (TF) ripple calculation.

        This routine computes the TF ripple amplitude at the midplane outboard
        plasma edge and the minimum radius of the TF coil centre that would
        produce a specified maximum allowed ripple. The calculation uses
        fitted coefficients derived from numerical modelling (MAGINT) and
        includes a simplified analytical picture-frame coil model for
        i_tf_shape == 2.

        Parameters
        ----------
        ripple_b_tf_plasma_edge_max : float
            Maximum allowed ripple at plasma edge (percent)
        r_tf_outboard_mid : float
            Radius to the centre of the outboard TF coil leg (m)
        n_tf_coils : int
            Number of TF coils
        rmajor : float
            Plasma major radius (m)
        rminor : float
            Plasma minor radius (m)
        r_tf_wp_inboard_inner : float
            Inner winding-pack inboard radius (m)
        r_tf_wp_inboard_centre : float
            Centre winding-pack inboard radius (m)
        r_tf_wp_inboard_outer : float
            Outer winding-pack inboard radius (m)
        dx_tf_wp_primary_toroidal : float
            Primary toroidal winding-pack thickness (m)
        i_tf_shape : int
            TF coil shape switch (2 => picture-frame analytical model)
        i_tf_sup : int
            TF coil support flag (1 => superconducting)
        dx_tf_wp_insulation : float
            Winding-pack insulation thickness (m)
        dx_tf_wp_insertion_gap : float
            Winding-pack insertion gap (m)
        i_tf_wp_geom: int
            Specifies type of shape for the winding pack
            See SuperconductingTFWPShapeType for more information

        Returns
        -------
        tuple[float, float, int]
            Tuple containing:
            - ripple: Calculated ripple at plasma edge (percent)
            - r_tf_outboard_midmin: Minimum r_tf_outboard_mid that yields the specified
                maximum ripple (m)
            - flag: Applicability flag (0 = OK, non-zero = fitted-range concern)

        Notes
        -----
        - Fitted coefficients originate from parametric MAGINT runs (M. Kovari, 2014).
        - Picture-frame coil analytical model (Ken McClements, 2022) is used when
        `i_tf_shape == 2` and gives approximate results (within ~10% of numerical).
        - The routine sets an applicability flag when fitted-range assumptions are
        exceeded.
        """
        if i_tf_sup == 1:
            # Minimal inboard WP radius [m]
            r_wp_min = r_tf_wp_inboard_inner

            i_tf_wp_geom = SuperconductingTFWPShapeType(i_tf_wp_geom)

            # Rectangular WP
            if i_tf_wp_geom == SuperconductingTFWPShapeType.RECTANGULAR:
                r_wp_max = r_wp_min

            # Double rectangle WP
            elif i_tf_wp_geom == SuperconductingTFWPShapeType.DOUBLE_RECTANGULAR:
                r_wp_max = r_tf_wp_inboard_centre

            # Trapezoidal WP
            elif i_tf_wp_geom == SuperconductingTFWPShapeType.TRAPEZOIDAL:
                r_wp_max = r_tf_wp_inboard_outer

            # Calculated maximum toroidal WP toroidal thickness [m]
            dx_tf_wp_conductor_max = dx_tf_wp_primary_toroidal - 2.0 * (
                dx_tf_wp_insulation + dx_tf_wp_insertion_gap
            )

        # Resistive magnet case
        else:
            # Radius used to define the dx_tf_wp_conductor_max [m]
            r_wp_max = r_tf_wp_inboard_outer
            # Calculated maximum toroidal WP toroidal thickness [m]
            dx_tf_wp_conductor_max = 2.0e0 * r_wp_max * np.tan(np.pi / n_tf_coils)

        flag = 0
        if i_tf_shape == TFCoilShapeModel.PICTURE_FRAME:
            # Ken McClements ST picture frame coil analytical ripple calc
            # Calculated ripple for coil at r_tf_outboard_mid (%)
            ripple = 100.0e0 * ((rmajor + rminor) / r_tf_outboard_mid) ** (n_tf_coils)
            #  Calculated r_tf_outboard_mid to produce a ripple of amplitude
            # ripple_b_tf_plasma_edge_max
            r_tf_outboard_midmin = (rmajor + rminor) / (
                (0.01e0 * ripple_b_tf_plasma_edge_max) ** (1.0e0 / n_tf_coils)
            )
        else:
            # Winding pack to iter-coil at plasma centre toroidal lenth ratio
            x = dx_tf_wp_conductor_max * n_tf_coils / rmajor

            # Fitting parameters
            c1 = 0.875e0 - 0.0557e0 * x
            c2 = 1.617e0 + 0.0832e0 * x

            #  Calculated ripple for coil at r_tf_outboard_mid (%)
            ripple = (
                100.0e0
                * c1
                * ((rmajor + rminor) / r_tf_outboard_mid) ** (n_tf_coils - c2)
            )

            #  Calculated r_tf_outboard_mid to produce a ripple of amplitude
            # ripple_b_tf_plasma_edge_max
            base = 0.01 * ripple_b_tf_plasma_edge_max / c1
            # Avoid potential negative or complex result: kludge base to be
            # small and positive if required
            if base <= 1e-6:
                logger.error("base is <= 1e-6. Kludging to 1e-6.")
                base = 1e-6

            r_tf_outboard_midmin = (rmajor + rminor) / (
                base ** (1.0 / (n_tf_coils - c2))
            )

            if np.isinf(r_tf_outboard_midmin):
                logger.error(
                    "r_tf_outboard_midmin is inf. Kludging to a large value instead."
                )
                r_tf_outboard_midmin = (rmajor + rminor) * 3

            #  Notify via flag if a range of applicability is violated
            flag = 0
            if (x < 0.737e0) or (x > 2.95e0):
                flag = 1
            if (n_tf_coils < 16) or (n_tf_coils > 20):
                flag = 2
            if ((rmajor + rminor) / r_tf_outboard_mid < 0.7e0) or (
                (rmajor + rminor) / r_tf_outboard_mid > 0.8e0
            ):
                flag = 3

        return ripple, r_tf_outboard_midmin, flag

    def calculate_radial_build(self, output: bool):
        """Method determining the radial build of the machine.
        It calculates various parameters related to the build of the machine,
        such as thicknesses, radii, and areas.
        Results can be outputted with the `output` flag.

        Parameters
        ----------
        output : bool
            Flag indicating whether to output the results
        """
        bld = self.data.build
        if self.data.fwbs.blktmodel > 0:
            bld.dr_blkt_inboard = bld.blbuith + bld.blbmith + bld.blbpith
            bld.dr_blkt_outboard = bld.blbuoth + bld.blbmoth + bld.blbpoth
            bld.dz_shld_upper = 0.5e0 * (bld.dr_shld_inboard + bld.dr_shld_outboard)

        #  Top/bottom blanket thickness
        bld.dz_blkt_upper = 0.5e0 * (bld.dr_blkt_inboard + bld.dr_blkt_outboard)

        i_single_null = DivertorNumberModels(self.data.physics.i_single_null)
        if i_single_null == DivertorNumberModels.SINGLE_NULL:
            #  Check if bld.dz_fw_plasma_gap has been set too small
            bld.dz_fw_plasma_gap = max(
                0.5e0 * (bld.dr_fw_plasma_gap_inboard + bld.dr_fw_plasma_gap_outboard),
                bld.dz_fw_plasma_gap,
            )

        # Issue #514 Radial dimensions of inboard leg
        # Calculate bld.dr_tf_inboard if
        # self.data.tfcoil.dr_tf_wp_with_insulation is an iteration variable (140)
        if 140 in self.data.numerics.ixc[0 : self.data.numerics.nvar]:
            bld.dr_tf_inboard = (
                self.data.tfcoil.dr_tf_wp_with_insulation
                + self.data.tfcoil.dr_tf_plasma_case
                + self.data.tfcoil.dr_tf_nose_case
            )

        if bld.i_tf_inside_cs == TFCSRadialConfiguration.TF_INSIDE_CS:
            bld.r_tf_inboard_in = bld.dr_bore
            # CS bore radius [m]
            bld.dr_cs_bore = bld.dr_bore + bld.dr_tf_inboard + bld.dr_cs_tf_gap
        else:
            bld.dr_cs_bore = bld.dr_bore

        # Calculate pre-compression structure thickness
        if (
            CSPrecompressionConfiguration(bld.i_cs_precomp)
            == CSPrecompressionConfiguration.CS_PRECOMPRESSION_STRUCTURE_PRESENT
        ):
            bld.dr_cs_precomp = bld.fseppc / (
                2.0e0
                * np.pi
                * bld.fcspc
                * bld.sigallpc
                * (2.0 * bld.dr_cs_bore + bld.dr_cs)
            )
        else:
            bld.dr_cs_precomp = 0.0e0

        if bld.i_tf_inside_cs != TFCSRadialConfiguration.TF_INSIDE_CS:
            # Inboard side inner radius [m]
            # This is not calculated above because it requires the dr_cs_precomp
            bld.r_tf_inboard_in = (
                bld.dr_bore + bld.dr_cs + bld.dr_cs_precomp + bld.dr_cs_tf_gap
            )

        # Radial build to tfcoil middle [m]
        bld.r_tf_inboard_mid = bld.r_tf_inboard_in + 0.5e0 * bld.dr_tf_inboard

        # Radial build to tfcoil plasma facing side [m]
        bld.r_tf_inboard_out = bld.r_tf_inboard_in + bld.dr_tf_inboard

        # WP radial thickness [m]
        # Calculated only if not used as an iteration variable
        if 140 not in self.data.numerics.ixc[0 : self.data.numerics.nvar]:
            self.data.tfcoil.dr_tf_wp_with_insulation = (
                bld.dr_tf_inboard
                - self.data.tfcoil.dr_tf_plasma_case
                - self.data.tfcoil.dr_tf_nose_case
            )

        # Radius of the centrepost at the top of the machine
        if self.data.physics.itart == 1 and self.data.tfcoil.i_tf_sup != 1:
            # bld.r_cp_top is set using the plasma shape
            if bld.i_r_cp_top == 0:
                bld.r_cp_top = (
                    self.data.physics.rmajor
                    - self.data.physics.rminor * self.data.physics.triang
                    - (
                        bld.dr_tf_shld_gap
                        + bld.dr_shld_thermal_inboard
                        + bld.dr_shld_inboard
                        + bld.dr_shld_blkt_gap
                        + bld.dr_blkt_inboard
                        + bld.dr_fw_inboard
                        + 3.0e0 * bld.dr_fw_plasma_gap_inboard
                    )
                    + self.data.tfcoil.drtop
                )

                # Notify user that bld.r_cp_top has been set to
                # 1.01*bld.r_tf_inboard_out (lvl 2 error)
                if bld.r_cp_top < 1.01e0 * bld.r_tf_inboard_out:
                    logger.error(
                        "TF CP top radius (r_cp_top) replaced by 1.01*r_tf_inboard_out "
                        "-> potential top rbuild issue"
                        f"{bld.r_cp_top=} "
                        f"{bld.r_tf_inboard_out=}"
                    )

                    # bld.r_cp_top correction
                    bld.r_cp_top = bld.r_tf_inboard_out * 1.01e0

                # Top and mid-plane TF coil CP radius ratio
                bld.f_r_cp = bld.r_cp_top / bld.r_tf_inboard_out

            # User defined bld.r_cp_top
            elif bld.i_r_cp_top == 1:
                # Notify user that bld.r_cp_top has been set to
                # 1.01*bld.r_tf_inboard_out (lvl 2 error)
                if bld.r_cp_top < 1.01e0 * bld.r_tf_inboard_out:
                    logger.error(
                        "TF CP top radius (r_cp_top) replaced by 1.01*r_tf_inboard_out "
                        "-> potential top rbuild issue"
                        f"{bld.r_cp_top=} "
                        f"{bld.r_tf_inboard_out=}"
                    )

                    # bld.r_cp_top correction
                    bld.r_cp_top = bld.r_tf_inboard_out * 1.01e0

                # Top / mid-plane TF CP radius ratio
                bld.f_r_cp = bld.r_cp_top / bld.r_tf_inboard_out

            # bld.r_cp_top set as a fraction of the outer TF midplane radius
            elif bld.i_r_cp_top == 2:
                bld.r_cp_top = bld.f_r_cp * bld.r_tf_inboard_out

        else:  # End of self.data.physics.itart == 1 .and. self.data.tfcoil.i_tf_sup /= 1
            bld.r_cp_top = bld.r_tf_inboard_out

        if bld.i_r_cp_top != 0 and (
            bld.r_cp_top
            > self.data.physics.rmajor
            - self.data.physics.rminor * self.data.physics.triang
            - (
                bld.dr_tf_shld_gap
                + bld.dr_shld_thermal_inboard
                + bld.dr_shld_inboard
                + bld.dr_shld_blkt_gap
                + bld.dr_blkt_inboard
                + bld.dr_fw_inboard
                + 3.0e0 * bld.dr_fw_plasma_gap_inboard
            )
            + self.data.tfcoil.drtop
        ):
            logger.error(
                "Top CP radius larger that its value determined with plasma shape "
                f"{bld.r_cp_top=}"
            )
        if bld.i_tf_inside_cs == TFCSRadialConfiguration.TF_INSIDE_CS:
            #  Radial position of vacuum vessel [m]
            bld.r_vv_inboard_out = (
                bld.r_tf_inboard_out
                + bld.dr_cs
                + bld.dr_cs_tf_gap
                + bld.dr_cs_precomp
                + bld.dr_tf_shld_gap
                + bld.dr_shld_thermal_inboard
                + bld.dr_shld_vv_gap_inboard
                + bld.dr_vv_inboard
            )
        else:
            bld.r_vv_inboard_out = (
                bld.r_tf_inboard_out
                + bld.dr_tf_shld_gap
                + bld.dr_shld_thermal_inboard
                + bld.dr_shld_vv_gap_inboard
                + bld.dr_vv_inboard
            )
        # Radial position of the inner side of inboard neutronic shield [m]
        bld.r_sh_inboard_in = bld.r_vv_inboard_out

        # Radial position of the plasma facing side of inboard neutronic shield [m]
        bld.r_sh_inboard_out = bld.r_sh_inboard_in + bld.dr_shld_inboard

        #  Radial build to centre of plasma (should be equal to self.data.physics.rmajor)
        bld.rbld = (
            bld.r_sh_inboard_out
            + bld.dr_shld_blkt_gap
            + bld.dr_blkt_inboard
            + bld.dr_fw_inboard
            + bld.dr_fw_plasma_gap_inboard
            + self.data.physics.rminor
        )

        #  Radius to inner edge of inboard shield
        bld.r_shld_inboard_inner = (
            self.data.physics.rmajor
            - self.data.physics.rminor
            - bld.dr_fw_plasma_gap_inboard
            - bld.dr_fw_inboard
            - bld.dr_blkt_inboard
            - bld.dr_shld_inboard
        )

        #  Radius to outer edge of outboard shield
        bld.r_shld_outboard_outer = (
            self.data.physics.rmajor
            + self.data.physics.rminor
            + bld.dr_fw_plasma_gap_outboard
            + bld.dr_fw_outboard
            + bld.dr_blkt_outboard
            + bld.dr_shld_outboard
        )

        #  Thickness of outboard TF coil legs
        if self.data.tfcoil.i_tf_sup != 1:
            bld.dr_tf_outboard = bld.f_dr_tf_outboard_inboard * bld.dr_tf_inboard
        else:
            bld.dr_tf_outboard = bld.dr_tf_inboard

        #  Radius to centre of outboard TF coil legs
        bld.r_tf_outboard_mid = (
            bld.r_shld_outboard_outer
            + bld.dr_shld_blkt_gap
            + bld.dr_vv_outboard
            + bld.gapomin
            + bld.dr_shld_thermal_outboard
            + bld.dr_tf_shld_gap
            + 0.5e0 * bld.dr_tf_outboard
        )

        # TF coil horizontal bore at mid-plane [m]
        bld.dr_tf_inner_bore = (bld.r_tf_outboard_mid - 0.5e0 * bld.dr_tf_outboard) - (
            bld.r_tf_inboard_mid - 0.5e0 * bld.dr_tf_inboard
        )

        (
            self.data.tfcoil.ripple_b_tf_plasma_edge,
            r_tf_outboard_midl,
            bld.ripflag,
        ) = self.plasma_outboard_edge_toroidal_ripple(
            ripple_b_tf_plasma_edge_max=self.data.tfcoil.ripple_b_tf_plasma_edge_max,
            r_tf_outboard_mid=bld.r_tf_outboard_mid,
            n_tf_coils=self.data.tfcoil.n_tf_coils,
            rmajor=self.data.physics.rmajor,
            rminor=self.data.physics.rminor,
            r_tf_wp_inboard_inner=self.data.superconducting_tfcoil.r_tf_wp_inboard_inner,
            r_tf_wp_inboard_centre=self.data.superconducting_tfcoil.r_tf_wp_inboard_centre,
            r_tf_wp_inboard_outer=self.data.superconducting_tfcoil.r_tf_wp_inboard_outer,
            dx_tf_wp_primary_toroidal=self.data.tfcoil.dx_tf_wp_primary_toroidal,
            i_tf_shape=self.data.tfcoil.i_tf_shape,
            i_tf_sup=self.data.tfcoil.i_tf_sup,
            dx_tf_wp_insulation=self.data.tfcoil.dx_tf_wp_insulation,
            dx_tf_wp_insertion_gap=self.data.tfcoil.dx_tf_wp_insertion_gap,
            i_tf_wp_geom=self.data.tfcoil.i_tf_wp_geom,
        )

        #  If the self.data.tfcoil.ripple is too large then move the outboard TF coil leg
        if r_tf_outboard_midl > bld.r_tf_outboard_mid:
            bld.r_tf_outboard_mid = r_tf_outboard_midl
            bld.dr_shld_vv_gap_outboard = (
                bld.r_tf_outboard_mid
                - 0.5e0 * bld.dr_tf_outboard
                - bld.dr_vv_outboard
                - bld.r_shld_outboard_outer
                - bld.dr_shld_thermal_outboard
                - bld.dr_tf_shld_gap
                - bld.dr_shld_blkt_gap
            )
            bld.dr_tf_inner_bore = (
                bld.r_tf_outboard_mid - 0.5e0 * bld.dr_tf_outboard
            ) - (bld.r_tf_inboard_mid - 0.5e0 * bld.dr_tf_inboard)
        else:
            bld.dr_shld_vv_gap_outboard = bld.gapomin

        (
            self.data.tfcoil.ripple_b_tf_plasma_edge,
            _r_tf_outboard_midl,
            self.data.build.ripflag,
        ) = self.plasma_outboard_edge_toroidal_ripple(
            ripple_b_tf_plasma_edge_max=self.data.tfcoil.ripple_b_tf_plasma_edge_max,
            r_tf_outboard_mid=self.data.build.r_tf_outboard_mid,
            n_tf_coils=self.data.tfcoil.n_tf_coils,
            rmajor=self.data.physics.rmajor,
            rminor=self.data.physics.rminor,
            r_tf_wp_inboard_inner=self.data.superconducting_tfcoil.r_tf_wp_inboard_inner,
            r_tf_wp_inboard_centre=self.data.superconducting_tfcoil.r_tf_wp_inboard_centre,
            r_tf_wp_inboard_outer=self.data.superconducting_tfcoil.r_tf_wp_inboard_outer,
            dx_tf_wp_primary_toroidal=self.data.tfcoil.dx_tf_wp_primary_toroidal,
            i_tf_shape=self.data.tfcoil.i_tf_shape,
            i_tf_sup=self.data.tfcoil.i_tf_sup,
            dx_tf_wp_insulation=self.data.tfcoil.dx_tf_wp_insulation,
            dx_tf_wp_insertion_gap=self.data.tfcoil.dx_tf_wp_insertion_gap,
            i_tf_wp_geom=self.data.tfcoil.i_tf_wp_geom,
        )

        if output:
            self.radial_build_output()

    def _ripple_flag_validation(self):
        po.ocmmnt(
            self.outfile,
            "(Ripple result may not be accurate, as the fit was outside",
        )
        po.ocmmnt(self.outfile, " its range of applicability.)")
        po.oblnkl(self.outfile)
        logger.warning(
            "Ripple result may be inaccurate, as the fit has been extrapolated"
        )

        if self.data.build.ripflag == 1:
            warning_str = (
                "(TF coil ripple calculation) "
                "Dimensionless coil width X out of fitted range. %s"
            )
            diagnostic = (
                self.data.tfcoil.dx_tf_wp_primary_toroidal
                * self.data.tfcoil.n_tf_coils
                / self.data.physics.rmajor
            )
        elif self.data.build.ripflag == 2:
            warning_str = (
                "(TF coil ripple calculation) "
                "No of TF coils not between 16 and 20 inclusive "
            )
            diagnostic = f"{self.data.tfcoil.n_tf_coils=}"
        else:
            diagnostic = (
                self.data.physics.rmajor + self.data.physics.rminor
            ) / self.data.build.r_tf_outboard_mid
            warning_str = (
                "(TF coil ripple calculation) (R+a)/rtot=%s out of fitted range.",
            )

        logger.warning(warning_str, diagnostic)

    def radial_build_output(self):
        """Print out device build"""
        bld = self.data.build
        po.oheadr(self.outfile, "Radial Build")

        if self.data.build.ripflag != 0:
            self._ripple_flag_validation()

        po.ovarre(
            self.outfile,
            "TF coil radial placement switch",
            "(i_tf_inside_cs)",
            self.data.build.i_tf_inside_cs,
        )
        po.ocmmnt(
            self.outfile,
            (
                "  -> "
                f"{TFCSRadialConfiguration(self.data.build.i_tf_inside_cs).description}"
            ),
        )
        po.oblnkl(self.outfile)
        po.ovarre(
            self.outfile,
            "Inboard build thickness (m)",
            "(dr_inboard_build)",
            self.data.physics.rmajor - self.data.physics.rminor,
            "OP ",
        )

        if (
            self.data.build.i_tf_inside_cs == TFCSRadialConfiguration.TF_INSIDE_CS
            and self.data.tfcoil.i_tf_bucking >= 2
        ):
            po.ocmmnt(
                self.outfile,
                "(Bore hollow space has been filled with a solid metal cyclinder to"
                " act as wedge support)\n",
            )

        # an array that holds the following information
        # description, variable name, thickness, radius
        radial_build_data = []

        radius = 0.0e0
        radial_build_data.append(["Device centreline", None, 0.0, radius])

        bore_descr = "Machine dr_bore"
        if bld.i_tf_inside_cs == TFCSRadialConfiguration.TF_INSIDE_CS:
            bore_descr += (
                "wedge support cylinder"
                if self.data.tfcoil.i_tf_bucking >= 2
                else "hole"
            )

        radial_build_data.append([bore_descr, "dr_bore", bld.dr_bore])

        if bld.i_tf_inside_cs == TFCSRadialConfiguration.TF_INSIDE_CS:
            radial_build_data.extend((
                [
                    "TF coil inboard leg (in dr_bore)",
                    "dr_tf_inboard",
                    bld.dr_tf_inboard,
                ],
                [
                    "CS precompresion to TF coil radial gap",
                    "dr_cs_tf_gap",
                    bld.dr_cs_tf_gap,
                ],
            ))

        radial_build_data.extend((
            ["Central solenoid", "dr_cs", bld.dr_cs],
            ["CS precompression", "dr_cs_precomp", bld.dr_cs_precomp],
        ))
        if bld.i_tf_inside_cs == TFCSRadialConfiguration.TF_OUTSIDE_CS:
            radial_build_data.extend((
                [
                    "CS precompresion to TF coil radial gap",
                    "dr_cs_tf_gap",
                    bld.dr_cs_tf_gap,
                ],
                [
                    "TF coil inboard leg",
                    "dr_tf_inboard",
                    bld.dr_tf_inboard,
                ],
            ))

        radial_build_data.extend((
            ["TF coil inboard leg insulation gap", "dr_tf_shld_gap", bld.dr_tf_shld_gap],
            [
                "Thermal shield, inboard",
                "dr_shld_thermal_inboard",
                bld.dr_shld_thermal_inboard,
            ],
            [
                "Thermal shield to vessel radial gap",
                "dr_shld_vv_gap_inboard",
                bld.dr_shld_vv_gap_inboard,
            ],
            ["Inboard vacuum vessel", "dr_vv_inboard", bld.dr_vv_inboard],
            ["Inner radiation shield", "dr_shld_inboard", bld.dr_shld_inboard],
            ["Gap", "dr_shld_blkt_gap", bld.dr_shld_blkt_gap],
            ["Inboard blanket", "dr_blkt_inboard", bld.dr_blkt_inboard],
            ["Inboard first wall", "dr_fw_inboard", bld.dr_fw_inboard],
            [
                "Inboard scrape-off",
                "dr_fw_plasma_gap_inboard",
                bld.dr_fw_plasma_gap_inboard,
            ],
            ["Plasma geometric centre", "rminor", self.data.physics.rminor],
            ["Plasma outboard edge", "rminor", self.data.physics.rminor],
            [
                "Outboard scrape-off",
                "dr_fw_plasma_gap_outboard",
                bld.dr_fw_plasma_gap_outboard,
            ],
            ["Outboard first wall", "dr_fw_outboard", bld.dr_fw_outboard],
            ["Outboard blanket", "dr_blkt_outboard", bld.dr_blkt_outboard],
            ["Gap", "dr_shld_blkt_gap", bld.dr_shld_blkt_gap],
            ["Outer radiation shield", "dr_shld_outboard", bld.dr_shld_outboard],
            ["Outboard vacuum vessel", "dr_vv_outboard", bld.dr_vv_outboard],
            ["Vessel to TF gap", "dr_shld_vv_gap_outboard", bld.dr_shld_vv_gap_outboard],
            [
                "Outboard thermal shield",
                "dr_shld_thermal_outboard",
                bld.dr_shld_thermal_outboard,
            ],
            ["Gap", "dr_tf_shld_gap", bld.dr_tf_shld_gap],
            ["TF coil outboard leg", "dr_tf_outboard", bld.dr_tf_outboard],
        ))

        for description, variable, thickness in radial_build_data:
            radius += thickness
            var = f"({variable})" if variable else ""
            po.obuild(self.outfile, description, thickness, radius, var)

        # use manual index to ensure count is contiguous in the event
        # of a `None` variable component
        index = 0
        for description, variable, thickness, radius in radial_build_data:
            if variable is None:
                continue

            index += 1

            po.ovarre(
                self.mfile,
                f"{description} radial thickness (m)",
                f"({variable})",
                thickness,
            )

            po.ovarre(
                self.mfile,
                f"Radial build component {index}",
                f"(radial_label({index}))",
                f'"{variable}"',
            )
            po.ovarre(
                self.mfile,
                f"Radial build cumulative radius {index}",
                f"(radial_cum({index}))",
                radius,
            )

        if (
            CurrentDriveModel(
                self.data.current_drive.i_hcd_primary
                or self.data.current_drive.i_hcd_secondary
            ).method
            == CurrentDriveMethodType.NEUTRAL_BEAM
        ):
            po.ovarre(
                self.mfile,
                "Width of neutral beam duct where it passes between the TF coils (m)",
                "(dx_beam_duct)",
                self.data.current_drive.dx_beam_duct,
            )
