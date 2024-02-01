import json
import numpy
import logging
import copy
import numba

from process.fortran import rebco_variables
from process.fortran import global_variables
from process.fortran import tfcoil_variables
from process.fortran import physics_variables
from process.fortran import build_variables
from process.fortran import constants
from process.fortran import sctfcoil_module
from process.fortran import process_output as po
from process.fortran import error_handling
from process.fortran import fwbs_variables
from process.fortran import pfcoil_variables
from process.fortran import numerics
from process.fortran import divertor_variables

import process.superconductors as superconductors

from process.utilities.f2py_string_patch import f2py_compatible_to_string


logger = logging.getLogger(__name__)


RMU0 = constants.rmu0
EPS = numpy.finfo(1.0).eps


class Sctfcoil:
    def __init__(self):
        self.outfile = constants.nout

    def run(self, output: bool):
        """
        Routine to call the superconductor module for the TF coils
        """
        tfes = sctfcoil_module.estotft / tfcoil_variables.n_tf
        # Cross-sectional area per turn
        aturn = tfcoil_variables.ritfc / (
            tfcoil_variables.jwptf * tfcoil_variables.n_tf * tfcoil_variables.n_tf_turn
        )

        if tfcoil_variables.i_tf_sc_mat == 6:
            (tfcoil_variables.jwdgcrt, tfcoil_variables.tmargtf,) = self.supercon_croco(
                aturn,
                tfcoil_variables.bmaxtfrp,
                tfcoil_variables.cpttf,
                tfcoil_variables.tftmp,
                output=output,
            )

            tfcoil_variables.vtfskv = (
                self.croco_voltage() / 1.0e3
            )  # TFC Quench voltage in kV

        else:
            (
                tfcoil_variables.jwdgcrt,
                vdump,
                tfcoil_variables.tmargtf,
            ) = self.supercon(
                tfcoil_variables.acstf,
                aturn,
                tfcoil_variables.bmaxtfrp,
                tfcoil_variables.vftf,
                tfcoil_variables.fcutfsu,
                tfcoil_variables.cpttf,
                tfcoil_variables.jwptf,
                tfcoil_variables.i_tf_sc_mat,
                tfcoil_variables.fhts,
                tfcoil_variables.tdmptf,
                tfes,
                tfcoil_variables.tftmp,
                tfcoil_variables.tmaxpro,
                tfcoil_variables.bcritsc,
                tfcoil_variables.tcritsc,
                output=output,
            )

            tfcoil_variables.vtfskv = vdump / 1.0e3  # TFC Quench voltage in kV

    def croco_voltage(self) -> float:
        if f2py_compatible_to_string(tfcoil_variables.quench_model) == "linear":
            sctfcoil_module.time2 = tfcoil_variables.tdmptf
            croco_voltage = (
                2.0e0
                / sctfcoil_module.time2
                * (sctfcoil_module.estotft / tfcoil_variables.n_tf)
                / tfcoil_variables.cpttf
            )
        elif f2py_compatible_to_string(tfcoil_variables.quench_model) == "exponential":
            sctfcoil_module.tau2 = tfcoil_variables.tdmptf
            croco_voltage = (
                2.0e0
                / sctfcoil_module.tau2
                * (sctfcoil_module.estotft / tfcoil_variables.n_tf)
                / tfcoil_variables.cpttf
            )
        else:
            return 0.0

        return croco_voltage

    def supercon_croco(self, aturn, bmax, iop, thelium, output: bool):
        """TF superconducting CroCo conductor using REBCO tape
        author: M Kovari, CCFE, Culham Science Centre
        bmax : input real : Peak field at conductor (T)
        iop : input real : Operating current per turn (A)
        thelium : input real : He temperature at peak field point (K)
        iprint : input integer : Switch for printing (1 = yes, 0 = no)
        outfile : input integer : Fortran output unit identifier
        jwdgcrt : output real : Critical winding pack current density (A/m2)
        tmarg : output real : Temperature margin (K)
        """

        jcritsc: float = 0.0
        #  Find critical current density in superconducting strand, jcritstr
        jcritsc, _ = superconductors.jcrit_rebco(thelium, bmax)
        # tfcoil_variables.acstf : Cable space - inside area (m2)
        # Set new rebco_variables.croco_od
        # allowing for scaling of rebco_variables.croco_od
        rebco_variables.croco_od = (
            tfcoil_variables.t_conductor / 3.0e0
            - tfcoil_variables.thwcndut * (2.0e0 / 3.0e0)
        )
        sctfcoil_module.conductor_acs = (
            9.0e0 / 4.0e0 * numpy.pi * rebco_variables.croco_od**2
        )
        tfcoil_variables.acstf = sctfcoil_module.conductor_acs
        sctfcoil_module.conductor_area = (
            tfcoil_variables.t_conductor**2
        )  # does this not assume it's a sqaure???

        sctfcoil_module.conductor_jacket_area = (
            sctfcoil_module.conductor_area - sctfcoil_module.conductor_acs
        )
        tfcoil_variables.acndttf = sctfcoil_module.conductor_jacket_area

        sctfcoil_module.conductor_jacket_fraction = (
            sctfcoil_module.conductor_jacket_area / sctfcoil_module.conductor_area
        )
        (
            sctfcoil_module.croco_strand_area,
            sctfcoil_module.croco_strand_critical_current,
            sctfcoil_module.conductor_copper_area,
            sctfcoil_module.conductor_copper_fraction,
            sctfcoil_module.conductor_copper_bar_area,
            sctfcoil_module.conductor_hastelloy_area,
            sctfcoil_module.conductor_hastelloy_fraction,
            sctfcoil_module.conductor_helium_area,
            sctfcoil_module.conductor_helium_fraction,
            sctfcoil_module.conductor_solder_area,
            sctfcoil_module.conductor_solder_fraction,
            sctfcoil_module.conductor_rebco_area,
            sctfcoil_module.conductor_rebco_fraction,
            sctfcoil_module.conductor_critical_current,
        ) = superconductors.croco(
            jcritsc,
            sctfcoil_module.conductor_area,
            rebco_variables.croco_od,
            rebco_variables.croco_thick,
        )

        rebco_variables.coppera_m2 = iop / sctfcoil_module.conductor_copper_area

        icrit = sctfcoil_module.conductor_critical_current
        jcritstr = (
            sctfcoil_module.croco_strand_critical_current
            / sctfcoil_module.croco_strand_area
        )

        # Critical current density in winding pack
        # aturn : Area per turn (i.e. entire jacketed conductor with insulation) (m2)
        jwdgcrt = icrit / aturn
        #  Ratio of operating / critical current
        iooic = iop / icrit
        #  Operating current density
        jwdgop = iop / aturn
        #  Actual current density in superconductor,
        # which should be equal to jcrit(thelium+tmarg)

        #  when we have found the desired value of tmarg
        jsc = iooic * jcritsc

        # Temperature margin using secant solver
        current_sharing_t = superconductors.current_sharing_rebco(bmax, jsc)
        tmarg = current_sharing_t - thelium
        tfcoil_variables.temp_margin = (
            tmarg  # Only used in the availabilty routine - see comment to Issue #526
        )

        if output:  # Output ----------------------------------
            total = (
                sctfcoil_module.conductor_copper_area
                + sctfcoil_module.conductor_hastelloy_area
                + sctfcoil_module.conductor_solder_area
                + sctfcoil_module.conductor_jacket_area
                + sctfcoil_module.conductor_helium_area
                + sctfcoil_module.conductor_rebco_area
            )

            if tfcoil_variables.temp_margin <= 0.0e0:
                logger.warning(
                    f"""Negative TFC temperature margin
                temp_margin: {tfcoil_variables.temp_margin}
                bmax: {bmax}"""
                )

            po.oheadr(self.outfile, "Superconducting TF Coils")
            po.ovarin(self.outfile, "Superconductor switch", "(isumat)", 6)
            po.ocmmnt(
                self.outfile, "Superconductor used: REBCO HTS tape in CroCo strand"
            )

            po.ovarre(
                self.outfile,
                "Thickness of REBCO layer in tape (m)",
                "(rebco_thickness)",
                rebco_variables.rebco_thickness,
            )
            po.ovarre(
                self.outfile,
                "Thickness of copper layer in tape (m)",
                "(copper_thick  )",
                rebco_variables.copper_thick,
            )
            po.ovarre(
                self.outfile,
                "Thickness of Hastelloy layer in tape (m) ",
                "(hastelloy_thickness)",
                rebco_variables.hastelloy_thickness,
            )

            po.ovarre(
                self.outfile,
                "Mean width of tape (m)",
                "(tape_width)",
                rebco_variables.tape_width,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Outer diameter of CroCo copper tube (m) ",
                "(croco_od)",
                rebco_variables.croco_od,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Inner diameter of CroCo copper tube (m) ",
                "(croco_id)",
                rebco_variables.croco_id,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Thickness of CroCo copper tube (m) ",
                "(croco_thick)",
                rebco_variables.croco_thick,
            )

            po.ovarre(
                self.outfile,
                "Thickness of each HTS tape ",
                "(tape_thickness)",
                rebco_variables.tape_thickness,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Thickness of stack of rebco_variables.tapes (m) ",
                "(stack_thickness)",
                rebco_variables.stack_thickness,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Number of rebco_variables.tapes in strand",
                "(tapes)",
                rebco_variables.tapes,
                "OP ",
            )
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Area of REBCO in strand (m2)",
                "(rebco_area)",
                rebco_variables.rebco_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of copper in strand (m2)",
                "(copper_area)",
                rebco_variables.copper_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of hastelloy substrate in strand (m2) ",
                "(hastelloy_area)",
                rebco_variables.hastelloy_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of solder in strand (m2)  ",
                "(solder_area)",
                rebco_variables.solder_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total: area of CroCo strand (m2)  ",
                "(croco_strand_area)",
                sctfcoil_module.croco_strand_area,
                "OP ",
            )
            if (
                abs(
                    sctfcoil_module.croco_strand_area
                    - (
                        rebco_variables.rebco_area
                        + rebco_variables.copper_area
                        + rebco_variables.hastelloy_area
                        + rebco_variables.solder_area
                    )
                )
                > 1e-6
            ):
                po.ocmmnt(self.outfile, "ERROR: Areas in CroCo strand do not add up")
                logger.warning("Areas in CroCo strand do not add up - see OUT.DAT")

            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Cable information")
            po.ovarin(
                self.outfile,
                "Number of CroCo strands in the cable (fixed) ",
                "",
                6,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total area of cable space (m2)",
                "(acstf)",
                tfcoil_variables.acstf,
                "OP ",
            )

            po.oblnkl(self.outfile)
            po.ocmmnt(
                self.outfile,
                "Conductor information (includes jacket, not including insulation)",
            )
            po.ovarre(
                self.outfile,
                "Width of square conductor (cable + steel jacket) (m)",
                "(t_conductor)",
                tfcoil_variables.t_conductor,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of conductor (m2)",
                "(area)",
                sctfcoil_module.conductor_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "REBCO area of conductor (mm2)",
                "(rebco_area)",
                sctfcoil_module.conductor_rebco_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Area of central copper bar (mm2)",
                "(copper_bar_area)",
                sctfcoil_module.conductor_copper_bar_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Total copper area of conductor, total (mm2)",
                "(copper_area)",
                sctfcoil_module.conductor_copper_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Hastelloy area of conductor (mm2)",
                "(hastelloy_area)",
                sctfcoil_module.conductor_hastelloy_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Solder area of conductor (mm2)",
                "(solder_area)",
                sctfcoil_module.conductor_solder_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Jacket area of conductor (mm2)",
                "(jacket_area)",
                sctfcoil_module.conductor_jacket_area,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Helium area of conductor (mm2)",
                "(helium_area)",
                sctfcoil_module.conductor_helium_area,
                "OP ",
            )
            if abs(total - sctfcoil_module.conductor_area) > 1e-8:
                po.ovarre(
                    self.outfile,
                    "ERROR: conductor areas do not add up:",
                    "(total)",
                    total,
                    "OP ",
                )
                logger.warning(f"conductor areas do not add up. total: {total}")

            po.ovarre(
                self.outfile,
                "Critical current of CroCo strand (A)",
                "(croco_strand_critical_current)",
                sctfcoil_module.croco_strand_critical_current,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current of conductor (A) ",
                "(conductor_critical_current)",
                sctfcoil_module.conductor_critical_current,
                "OP ",
            )

            if global_variables.run_tests == 1:
                po.oblnkl(self.outfile)
                po.ocmmnt(
                    self.outfile,
                    "PROCESS TF Coil peak field fit. Values for t, z and y:",
                )
                po.oblnkl(self.outfile)
                po.ovarre(
                    self.outfile,
                    "Dimensionless winding pack width",
                    "(tf_fit_t)",
                    sctfcoil_module.tf_fit_t,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Dimensionless winding pack radial thickness",
                    "(tf_fit_z)",
                    sctfcoil_module.tf_fit_z,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Ratio of actual peak field to nominal axisymmetric peak field",
                    "(tf_fit_y)",
                    sctfcoil_module.tf_fit_y,
                    "OP ",
                )

            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Helium temperature at peak field (= superconductor temperature) (K)",
                "(thelium)",
                thelium,
            )
            po.ovarre(
                self.outfile,
                "Critical current density in superconductor (A/m2)",
                "(jcritsc)",
                jcritsc,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current density in strand (A/m2)",
                "(jcritstr)",
                jcritstr,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current density in winding pack (A/m2)",
                "(jwdgcrt)",
                jwdgcrt,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current density in winding pack (A/m2)",
                "(jwdgop)",
                jwdgop,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Minimum allowed temperature margin in superconductor (K)",
                "(tmargmin_tf)",
                tfcoil_variables.tmargmin_tf,
            )

            po.ovarre(
                self.outfile,
                "Actual temperature margin in superconductor (K)",
                "(tmarg)",
                tmarg,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Current sharing temperature (K)",
                "(current_sharing_t)",
                current_sharing_t,
                "OP ",
            )
            po.ovarre(self.outfile, "Critical current (A)", "(icrit)", icrit, "OP ")
            po.ovarre(
                self.outfile,
                "Actual current (A)",
                "(cpttf)",
                tfcoil_variables.cpttf,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current / critical current",
                "(iooic)",
                iooic,
                "OP ",
            )

        return jwdgcrt, tmarg

    def supercon(
        self,
        acs,
        aturn,
        bmax,
        fhe,
        fcu,
        iop,
        jwp,
        isumat,
        fhts,
        tdmptf,
        tfes,
        thelium,
        tmax,
        bcritsc,
        tcritsc,
        output: bool,
    ):
        """Routine to calculate the TF superconducting conductor  properties
        author: P J Knight, CCFE, Culham Science Centre
        author: J Galambos, ORNL
        author: R Kemp, CCFE, Culham Science Centre
        author: M Kovari, CCFE, Culham Science Centre
        author: J Miller, ORNL
        acs : input real : Cable space - inside area (m2)
        aturn : input real : Area per turn (i.e. entire jacketed conductor) (m2)
        bmax : input real : Peak field at conductor (T)
        fhe : input real : Fraction of cable space that is for He cooling
        fcu : input real : Fraction of conductor that is copper
        iop : input real : Operating current per turn (A)
        jwp : input real : Actual winding pack current density (A/m2)
        isumat : input integer : Switch for conductor type:
        1 = ITER Nb3Sn, standard parameters,
        2 = Bi-2212 High Temperature Superconductor,
        3 = NbTi,
        4 = ITER Nb3Sn, user-defined parameters
        5 = WST Nb3Sn parameterisation
        7 = Durham Ginzburg-Landau Nb-Ti parameterisation
        fhts    : input real : Adjustment factor (<= 1) to account for strain,
        radiation damage, fatigue or AC losses
        tdmptf : input real : Dump time (sec)
        tfes : input real : Energy stored in one TF coil (J)
        thelium : input real : He temperature at peak field point (K)
        tmax : input real : Max conductor temperature during quench (K)
        bcritsc : input real : Critical field at zero temperature and strain (T) (isumat=4 only)
        tcritsc : input real : Critical temperature at zero field and strain (K) (isumat=4 only)
        iprint : input integer : Switch for printing (1 = yes, 0 = no)
        outfile : input integer : Fortran output unit identifier
        jwdgpro : output real : Winding pack current density from temperature
        rise protection (A/m2)
        jwdgcrt : output real : Critical winding pack current density (A/m2)
        vd : output real : Discharge voltage imposed on a TF coil (V)
        tmarg : output real : Temperature margin (K)
        This routine calculates the superconductor properties for the TF coils.
        It was originally programmed by J. Galambos 1991, from algorithms provided
        by J. Miller.
        <P>The routine calculates the critical current density (winding pack)
        and also the protection information (for a quench).
        NOT used for the Croco conductor
        """
        tdump = tdmptf

        # Helium channel
        fhetot = (
            fhe
            + (numpy.pi / 4.0e0)
            * tfcoil_variables.dhecoil
            * tfcoil_variables.dhecoil
            / acs
        )

        # Guard against negative conductor fraction fcond
        # Kludge to allow solver to continue and hopefully be constrained away
        # from this point
        if fhetot > 0.99:
            fhetot = 0.99

        #  Conductor fraction (including central helium channel)
        fcond = 1.0e0 - fhetot

        if tfcoil_variables.i_str_wp == 0:
            strain = tfcoil_variables.str_tf_con_res
        else:
            strain = tfcoil_variables.str_wp

        #  Find critical current density in superconducting strand, jcritstr

        if isumat == 1:  # ITER Nb3Sn critical surface parameterization
            bc20m = 32.97e0
            tc0m = 16.06e0
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.5e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = numpy.sign(strain) * 0.5e-2

            #  jcritsc returned by superconductors.itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper
            jcritsc, _, _ = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)
            #  Critical current in cable
            icrit = jcritstr * acs * fcond

        elif isumat == 2:  # Bi-2212 high temperature superconductor parameterization
            #  Current density in a strand of Bi-2212 conductor
            #  N.B. jcrit returned by superconductors.bi2212 is the critical current density
            #  in the strand, not just the superconducting portion.
            #  The parameterization for jcritstr assumes a particular strand
            #  composition that does not require a user-defined copper fraction,
            #  so this is irrelevant in this model
            jstrand = jwp * aturn / (acs * fcond)

            jcritstr, tmarg = superconductors.bi2212(bmax, jstrand, thelium, fhts)
            jcritsc = jcritstr / (1.0e0 - fcu)
            #  Critical current in cable
            icrit = jcritstr * acs * fcond

        elif isumat == 3:  # NbTi data
            bc20m = 15.0e0
            tc0m = 9.3e0
            c0 = 1.0e10
            jcritsc, _ = superconductors.jcrit_nbti(thelium, bmax, c0, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)
            #  Critical current in cable
            icrit = jcritstr * acs * fcond

        elif isumat == 4:  # ITER Nb3Sn parameterization, but user-defined parameters
            bc20m = bcritsc
            tc0m = tcritsc
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.5e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = numpy.sign(strain) * 0.5e-2

            jcritsc, _, _ = superconductors.itersc(thelium, bmax, strain, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)
            #  Critical current in cable
            icrit = jcritstr * acs * fcond

        elif isumat == 5:  # WST Nb3Sn parameterisation
            bc20m = 32.97e0
            tc0m = 16.06e0
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.5e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = numpy.sign(strain) * 0.5e-2

            #  jcritsc returned by superconductors.itersc is the critical current density in the
            #  superconductor - not the whole strand, which contains copper
            jcritsc, _, _ = superconductors.wstsc(thelium, bmax, strain, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)
            #  Critical current in cable
            icrit = jcritstr * acs * fcond

        elif isumat == 6:  # "REBCO" 2nd generation HTS superconductor in CrCo strand
            raise ValueError(
                "sctfcoil.supercon has been called but tfcoil_variables.i_tf_sc_mat=6"
            )

        elif isumat == 7:  # Durham Ginzburg-Landau Nb-Ti parameterisation
            bc20m = tfcoil_variables.b_crit_upper_nbti
            tc0m = tfcoil_variables.t_crit_nbti
            jcritsc, _, _ = superconductors.gl_nbti(thelium, bmax, strain, bc20m, tc0m)
            jcritstr = jcritsc * (1.0e0 - fcu)
            #  Critical current in cable
            icrit = jcritstr * acs * fcond

        elif isumat == 8:  # Durham Ginzburg-Landau critical surface model for REBCO
            bc20m = 430
            tc0m = 185
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.7e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = numpy.sign(strain) * 0.7e-2

            jcritsc, _, _ = superconductors.gl_rebco(thelium, bmax, strain, bc20m, tc0m)
            # A0 calculated for tape cross section already
            jcritstr = jcritsc * (1.0e0 - fcu)
            #  Critical current in cable (copper added at this stage in HTS cables)
            icrit = jcritstr * acs * fcond

        elif (
            isumat == 9
        ):  # Hazelton experimental data + Zhai conceptual model for REBCO
            bc20m = 138
            tc0m = 92
            # If strain limit achieved, throw a warning and use the lower strain
            if abs(strain) > 0.7e-2:
                error_handling.fdiags[0] = strain
                error_handling.report_error(261)
                strain = numpy.sign(strain) * 0.7e-2

            # 'high current density' as per parameterisation described in Wolf,
            #  and based on Hazelton experimental data and Zhai conceptual model;
            #  see subroutine for full references
            jcritsc, _, _ = superconductors.hijc_rebco(
                thelium, bmax, strain, bc20m, tc0m
            )

            jcritstr = jcritsc * (1.0e0 - fcu)
            #  Critical current in cable (copper added at this stage in HTS cables)
            icrit = jcritstr * acs * fcond

        else:
            error_handling.idiags[0] = isumat
            error_handling.report_error(105)

        # Critical current density in winding pack
        # aturn : Area per turn (i.e. entire jacketed conductor with insulation) (m2)
        jwdgcrt = icrit / aturn
        #  Ratio of operating / critical current
        iooic = iop / icrit
        #  Operating current density
        jwdgop = iop / aturn
        #  Actual current density in superconductor, which should be equal to jcrit(thelium+tmarg)
        #  when we have found the desired value of tmarg
        jsc = iooic * jcritsc

        if iooic <= 0e0:
            logger.warning(
                f"""Negative Iop/Icrit for TF coil
            jsc: {jsc}
            iooic: {iooic}
            jcritsc: {jcritsc}
            Check conductor dimensions. fcond likely gone negative. fcond: {fcond}
            """
            )

        # REBCO measurements from 2 T to 14 T, extrapolating outside this
        if (isumat == 8) and (tfcoil_variables.bmaxtfrp >= 14):
            error_handling.report_error(266)

        #  Temperature margin (already calculated in superconductors.bi2212 for isumat=2)
        if (
            (isumat == 1)
            or (isumat == 3)
            or (isumat == 4)
            or (isumat == 5)
            or (isumat == 7)
            or (isumat == 8)
            or (isumat == 9)
        ):
            #  Newton-Raphson method; start approx at requested minimum temperature margin
            ttest = thelium + tfcoil_variables.tmargmin_tf + 0.001e0
            delt = 0.01e0
            jtol = 1.0e4

            for lap in range(100):
                if ttest <= 0.0e0:
                    error_handling.idiags[0] = lap
                    error_handling.fdiags[0] = ttest
                    error_handling.report_error(157)
                    break

                # Calculate derivative numerically
                ttestm = ttest - delt
                ttestp = ttest + delt

                # Issue #483 to be on the safe side, check the fractional as well as the absolute error
                if isumat in (1, 4):
                    jcrit0, _, _ = superconductors.itersc(
                        ttest, bmax, strain, bc20m, tc0m
                    )
                    if (abs(jsc - jcrit0) <= jtol) and (
                        abs((jsc - jcrit0) / jsc) <= 0.01
                    ):
                        break
                    jcritm, _, _ = superconductors.itersc(
                        ttestm, bmax, strain, bc20m, tc0m
                    )
                    jcritp, _, _ = superconductors.itersc(
                        ttestp, bmax, strain, bc20m, tc0m
                    )
                elif isumat == 3:
                    jcrit0, _ = superconductors.jcrit_nbti(ttest, bmax, c0, bc20m, tc0m)
                    if (abs(jsc - jcrit0) <= jtol) and (
                        abs((jsc - jcrit0) / jsc) <= 0.01
                    ):
                        break
                    jcritm, _ = super.jcrit_nbti(ttestm, bmax, c0, bc20m, tc0m)
                    jcritp, _ = superconductors.jcrit_nbti(
                        ttestp, bmax, c0, bc20m, tc0m
                    )
                elif isumat == 5:
                    jcrit0, _, _ = superconductors.wstsc(
                        ttest, bmax, strain, bc20m, tc0m
                    )
                    if (abs(jsc - jcrit0) <= jtol) and (
                        abs((jsc - jcrit0) / jsc) <= 0.01
                    ):
                        break
                    jcritm, _, _ = superconductors.wstsc(
                        ttestm, bmax, strain, bc20m, tc0m
                    )
                    jcritp, _, _ = superconductors.wstsc(
                        ttestp, bmax, strain, bc20m, tc0m
                    )
                elif isumat == 7:
                    jcrit0, _, _ = superconductors.gl_nbti(
                        ttest, bmax, strain, bc20m, tc0m
                    )
                    if (abs(jsc - jcrit0) <= jtol) and (
                        abs((jsc - jcrit0) / jsc) <= 0.01
                    ):
                        break
                    jcritm, _, _ = superconductors.gl_nbti(
                        ttestm, bmax, strain, bc20m, tc0m
                    )
                    jcritp, _, _ = superconductors.gl_nbti(
                        ttestp, bmax, strain, bc20m, tc0m
                    )
                elif isumat == 8:
                    jcrit0, _, _ = superconductors.gl_rebco(
                        ttest, bmax, strain, bc20m, tc0m
                    )
                    if (abs(jsc - jcrit0) <= jtol) and (
                        abs((jsc - jcrit0) / jsc) <= 0.01
                    ):
                        break
                    jcritm, _, _ = superconductors.gl_rebco(
                        ttestm, bmax, strain, bc20m, tc0m
                    )
                    jcritp, _, _ = superconductors.gl_rebco(
                        ttestp, bmax, strain, bc20m, tc0m
                    )
                elif isumat == 9:
                    jcrit0, _, _ = superconductors.hijc_rebco(
                        ttest, bmax, strain, bc20m, tc0m
                    )
                    if (abs(jsc - jcrit0) <= jtol) and (
                        abs((jsc - jcrit0) / jsc) <= 0.01
                    ):
                        break
                    jcritm, _, _ = superconductors.hijc_rebco(
                        ttestm, bmax, strain, bc20m, tc0m
                    )
                    jcritp, _, _ = superconductors.hijc_rebco(
                        ttestp, bmax, strain, bc20m, tc0m
                    )

                ttest = ttest - 2.0e0 * delt * (jcrit0 - jsc) / (jcritp - jcritm)
            else:
                error_handling.idiags[0] = lap
                error_handling.fdiags[0] = ttest
                error_handling.report_error(157)

            tmarg = ttest - thelium
            tfcoil_variables.temp_margin = tmarg

        #  Find the current density limited by the protection limit
        #  (N.B. Unclear of this routine's relevance for Bi-2212 (isumat=2), due
        #  to presence of fcu argument, which is not used for this model above)

        tfcoil_variables.jwdgpro, vd = self.protect(
            iop, tfes, acs, aturn, tdump, fcond, fcu, thelium, tmax
        )

        if output:  # Output --------------------------
            if ttest <= 0.0e0:
                logger.warning(
                    """Negative TFC temperature margin
                ttest: {ttest}
                bmax: {bmax}
                jcrit0: {jcrit0}
                jsc: {jsc}
                ttestp: {ttestp}
                ttestm: {ttestm}
                jcritp: {jcritp}
                jcritm: {jcritm}
                """
                )

            po.oheadr(self.outfile, "Superconducting TF Coils")
            po.ovarin(self.outfile, "Superconductor switch", "(isumat)", isumat)

            if isumat == 1:
                po.ocmmnt(self.outfile, "Superconductor used: Nb3Sn")
                po.ocmmnt(self.outfile, "  (ITER Jcrit model, standard parameters)")
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if isumat == 2:
                po.ocmmnt(self.outfile, "Superconductor used: Bi-2212 HTS")
            if isumat == 3:
                po.ocmmnt(self.outfile, "Superconductor used: NbTi")
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if isumat == 4:
                po.ocmmnt(self.outfile, "Superconductor used: Nb3Sn")
                po.ocmmnt(self.outfile, "  (ITER Jcrit model, user-defined parameters)")
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if isumat == 5:
                po.ocmmnt(self.outfile, "Superconductor used: Nb3Sn")
                po.ocmmnt(self.outfile, " (WST Nb3Sn critical surface model)")
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if isumat == 7:
                po.ocmmnt(self.outfile, "Superconductor used: Nb-Ti")
                po.ocmmnt(
                    self.outfile, " (Durham Ginzburg-Landau critical surface model)"
                )
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if isumat == 8:
                po.ocmmnt(self.outfile, "Superconductor used: REBCO")
                po.ocmmnt(
                    self.outfile, " (Durham Ginzburg-Landau critical surface model)"
                )
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )
            if isumat == 9:
                po.ocmmnt(self.outfile, "Superconductor used: REBCO")
                po.ocmmnt(
                    self.outfile,
                    " (Hazelton experimental data + Zhai conceptual model)",
                )
                po.ovarre(
                    self.outfile,
                    "Critical field at zero temperature and strain (T)",
                    "(bc20m)",
                    bc20m,
                )
                po.ovarre(
                    self.outfile,
                    "Critical temperature at zero field and strain (K)",
                    "(tc0m)",
                    tc0m,
                )

            if global_variables.run_tests == 1:
                po.oblnkl(self.outfile)
                po.ocmmnt(
                    self.outfile,
                    "PROCESS TF Coil peak field fit. Values for t, z and y:",
                )
                po.oblnkl(self.outfile)
                po.ovarre(
                    self.outfile,
                    "Dimensionless winding pack width",
                    "(tf_fit_t)",
                    sctfcoil_module.tf_fit_t,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Dimensionless winding pack radial thickness",
                    "(tf_fit_z)",
                    sctfcoil_module.tf_fit_z,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Ratio of peak field with ripple to nominal axisymmetric peak field",
                    "(tf_fit_y)",
                    sctfcoil_module.tf_fit_y,
                    "OP ",
                )

            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Helium temperature at peak field (= superconductor temperature) (K)",
                "(thelium)",
                thelium,
            )
            po.ovarre(
                self.outfile,
                "Total helium fraction inside cable space",
                "(fhetot)",
                fhetot,
                "OP ",
            )
            po.ovarre(self.outfile, "Copper fraction of conductor", "(fcutfsu)", fcu)
            po.ovarre(
                self.outfile,
                "Residual manufacturing strain on superconductor",
                "(str_tf_con_res)",
                tfcoil_variables.str_tf_con_res,
            )
            po.ovarre(
                self.outfile,
                "Self-consistent strain on superconductor",
                "(str_wp)",
                tfcoil_variables.str_wp,
            )
            po.ovarre(
                self.outfile,
                "Critical current density in superconductor (A/m2)",
                "(jcritsc)",
                jcritsc,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current density in strand (A/m2)",
                "(jcritstr)",
                jcritstr,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Critical current density in winding pack (A/m2)",
                "(jwdgcrt)",
                jwdgcrt,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current density in winding pack (A/m2)",
                "(jwdgop)",
                jwdgop,
                "OP ",
            )

            po.ovarre(
                self.outfile,
                "Minimum allowed temperature margin in superconductor (K)",
                "(tmargmin_tf)",
                tfcoil_variables.tmargmin_tf,
            )
            po.ovarre(
                self.outfile,
                "Actual temperature margin in superconductor (K)",
                "(tmarg)",
                tmarg,
                "OP ",
            )
            po.ovarre(self.outfile, "Critical current (A)", "(icrit)", icrit, "OP ")
            po.ovarre(
                self.outfile,
                "Actual current (A)",
                "(cpttf)",
                tfcoil_variables.cpttf,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Actual current / critical current",
                "(iooic)",
                iooic,
                "OP ",
            )

        return jwdgcrt, vd, tmarg

    def sctfcoil(self, output: bool):
        """TF coil module
        author: P J Knight, CCFE, Culham Science Centre
        author: J Galambos, FEDC/ORNL
        author: R Kemp, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        author: S Kahn, CCFE, Culham Science Centre
        This subroutine calculates various parameters for a TF coil set.
        The primary outputs are coil size, shape, weight, stress and and fields.
        It is a variant from the original FEDC/Tokamak systems code.
        """
        peaktfflag = 0

        self.tf_global_geometry()

        # Calculation of the TF current from bt
        self.tf_current()

        # Conductor section internal geometry
        # ---
        # Superconducting magnets
        if tfcoil_variables.i_tf_sup == 1:
            self.sc_tf_internal_geom(
                tfcoil_variables.i_tf_wp_geom,
                tfcoil_variables.i_tf_case_geom,
                tfcoil_variables.i_tf_turns_integer,
            )

        # Resitive magnets
        else:
            self.res_tf_internal_geom()

        # ---

        # Coil vertical geometry
        self.coilshap()

        # TF resistive heating (res TF only)
        if tfcoil_variables.i_tf_sup != 1:
            self.tf_res_heating()

        # Vertical force
        self.tf_field_and_force()

        # TF coil inductance
        # ---
        if physics_variables.itart == 0 and tfcoil_variables.i_tf_shape == 1:
            tfcoil_variables.tfind = self.tfcind(
                build_variables.tfcth, tfcoil_variables.xarc, tfcoil_variables.yarc
            )
        else:
            tfcoil_variables.tfind = (
                (build_variables.hmax + build_variables.tfthko)
                * RMU0
                / constants.pi
                * numpy.log(
                    build_variables.r_tf_outboard_mid / build_variables.r_tf_inboard_mid
                )
            )

        # Total TF coil stored magnetic energy [J]
        sctfcoil_module.estotft = (
            0.5e0 * tfcoil_variables.tfind * tfcoil_variables.ritfc**2
        )

        # Total TF coil stored magnetic energy [Gigajoule]
        tfcoil_variables.estotftgj = 1.0e-9 * sctfcoil_module.estotft
        # ---

        # Calculate TF coil areas and masses
        self.tf_coil_area_and_masses()

        # Peak field including ripple
        # Rem : as resistive magnets are axisymmetric, no inboard ripple is present
        if tfcoil_variables.i_tf_sup == 1:
            tfcoil_variables.bmaxtfrp, peaktfflag = self.peak_tf_with_ripple(
                tfcoil_variables.n_tf,
                tfcoil_variables.wwp1,
                tfcoil_variables.dr_tf_wp
                - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap),
                sctfcoil_module.r_wp_centre,
                tfcoil_variables.bmaxtf,
            )
        else:
            tfcoil_variables.bmaxtfrp = tfcoil_variables.bmaxtf

        # Do stress calculations (writes the stress output)
        if output:
            tfcoil_variables.n_rad_per_layer = 500

        try:
            (
                sig_tf_r_max,
                sig_tf_t_max,
                sig_tf_z_max,
                sig_tf_vmises_max,
                sig_tf_tresca_max,
                deflect,
                eyoung_axial,
                eyoung_trans,
                eyoung_wp_axial,
                eyoung_wp_trans,
                poisson_wp_trans,
                radial_array,
                s_tresca_cond_cea,
                poisson_wp_axial,
                sig_tf_r,
                sig_tf_smeared_r,
                sig_tf_smeared_t,
                sig_tf_smeared_z,
                sig_tf_t,
                sig_tf_tresca,
                sig_tf_vmises,
                sig_tf_z,
                str_tf_r,
                str_tf_t,
                str_tf_z,
                n_radial_array,
                n_tf_bucking,
                tfcoil_variables.sig_tf_wp,
                sig_tf_case,
                sig_tf_cs_bucked,
                str_wp,
                casestr,
                insstrain,
                sig_tf_wp_av_z,
            ) = self.stresscl(
                int(tfcoil_variables.n_tf_stress_layers),
                int(tfcoil_variables.n_rad_per_layer),
                int(tfcoil_variables.n_tf_wp_layers),
                int(tfcoil_variables.i_tf_bucking),
                float(build_variables.r_tf_inboard_in),
                build_variables.bore,
                build_variables.hmax,
                pfcoil_variables.ohhghf,
                build_variables.ohcth,
                build_variables.tf_in_cs,
                build_variables.tfcth,
                build_variables.gapoh,
                pfcoil_variables.ipfres,
                pfcoil_variables.coheof,
                pfcoil_variables.cohbop,
                pfcoil_variables.cptdin,
                pfcoil_variables.ncls,
                pfcoil_variables.ld_ratio_cst,
                pfcoil_variables.r_out_cst,
                pfcoil_variables.oh_steel_frac,
                tfcoil_variables.eyoung_steel,
                tfcoil_variables.poisson_steel,
                tfcoil_variables.eyoung_cond_axial,
                tfcoil_variables.poisson_cond_axial,
                tfcoil_variables.eyoung_cond_trans,
                tfcoil_variables.poisson_cond_trans,
                tfcoil_variables.eyoung_ins,
                tfcoil_variables.poisson_ins,
                tfcoil_variables.thicndut,
                tfcoil_variables.eyoung_copper,
                tfcoil_variables.poisson_copper,
                tfcoil_variables.i_tf_sup,
                tfcoil_variables.eyoung_res_tf_buck,
                sctfcoil_module.r_wp_inner,
                sctfcoil_module.tan_theta_coil,
                sctfcoil_module.theta_coil,
                sctfcoil_module.r_wp_outer,
                sctfcoil_module.a_tf_steel,
                sctfcoil_module.a_case_front,
                sctfcoil_module.a_case_nose,
                tfcoil_variables.tfinsgap,
                tfcoil_variables.tinstf,
                tfcoil_variables.n_tf_turn,
                int(tfcoil_variables.i_tf_turns_integer),
                sctfcoil_module.t_cable,
                sctfcoil_module.t_cable_radial,
                tfcoil_variables.dhecoil,
                tfcoil_variables.fcutfsu,
                tfcoil_variables.thwcndut,
                sctfcoil_module.t_lat_case_av,
                sctfcoil_module.t_wp_toroidal_av,
                sctfcoil_module.a_tf_ins,
                tfcoil_variables.aswp,
                tfcoil_variables.acond,
                sctfcoil_module.awpc,
                tfcoil_variables.eyoung_al,
                tfcoil_variables.poisson_al,
                tfcoil_variables.fcoolcp,
                tfcoil_variables.n_tf_graded_layers,
                tfcoil_variables.ritfc,
                tfcoil_variables.casthi,
                tfcoil_variables.i_tf_stress_model,
                sctfcoil_module.vforce_inboard_tot,
                tfcoil_variables.i_tf_tresca,
                tfcoil_variables.acasetf,
                tfcoil_variables.vforce,
                tfcoil_variables.acndttf,
            )

            tfcoil_variables.sig_tf_case = (
                tfcoil_variables.sig_tf_case
                if tfcoil_variables.sig_tf_case is None
                else sig_tf_case
            )

            tfcoil_variables.sig_tf_cs_bucked = (
                tfcoil_variables.sig_tf_cs_bucked
                if tfcoil_variables.sig_tf_cs_bucked is None
                else sig_tf_cs_bucked
            )

            tfcoil_variables.str_wp = (
                tfcoil_variables.str_wp if tfcoil_variables.str_wp is None else str_wp
            )

            tfcoil_variables.casestr = (
                tfcoil_variables.casestr
                if tfcoil_variables.casestr is None
                else casestr
            )

            tfcoil_variables.insstrain = (
                tfcoil_variables.insstrain
                if tfcoil_variables.insstrain is None
                else insstrain
            )

            if output:
                self.out_stress(
                    sig_tf_r_max,
                    sig_tf_t_max,
                    sig_tf_z_max,
                    sig_tf_vmises_max,
                    sig_tf_tresca_max,
                    deflect,
                    eyoung_axial,
                    eyoung_trans,
                    eyoung_wp_axial,
                    eyoung_wp_trans,
                    poisson_wp_trans,
                    radial_array,
                    s_tresca_cond_cea,
                    poisson_wp_axial,
                    sig_tf_r,
                    sig_tf_smeared_r,
                    sig_tf_smeared_t,
                    sig_tf_smeared_z,
                    sig_tf_t,
                    sig_tf_tresca,
                    sig_tf_vmises,
                    sig_tf_z,
                    str_tf_r,
                    str_tf_t,
                    str_tf_z,
                    n_radial_array,
                    n_tf_bucking,
                    sig_tf_wp_av_z,
                )
        except ValueError as e:
            if e.args[1] == 245 and e.args[2] == 0:
                error_handling.report_error(245)
                tfcoil_variables.sig_tf_case = 0.0e0
                tfcoil_variables.sig_tf_wp = 0.0e0

        if output:
            self.outtf(peaktfflag)

    def protect(self, aio, tfes, acs, aturn, tdump, fcond, fcu, tba, tmax):
        """Finds the current density limited by the protection limit
        author: P J Knight, CCFE, Culham Science Centre
        author: J Miller, ORNL
        aio : input real : Operating current (A)
        tfes : input real : Energy stored in one TF coil (J)
        acs : input real : Cable space - inside area (m2)
        aturn : input real : Area per turn (i.e.  entire cable) (m2)
        tdump : input real : Dump time (sec)
        fcond : input real : Fraction of cable space containing conductor
        fcu : input real : Fraction of conductor that is copper
        tba : input real : He temperature at peak field point (K)
        tmax : input real : Max conductor temperature during quench (K)
        ajwpro : output real :  Winding pack current density from temperature
        rise protection (A/m2)
        vd : output real :  Discharge voltage imposed on a TF coil (V)
        This routine calculates maximum conductor current density which
        limits the peak temperature in the winding to a given limit (tmax).
        It also finds the dump voltage.
        <P>These calculations are based on Miller's formulations.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        # Integration coefficients p1,p2,p3
        p1 = (
            0.0e0,
            0.8e0,
            1.75e0,
            2.4e0,
            2.7e0,
            2.95e0,
            3.1e0,
            3.2e0,
            3.3e0,
            3.4e0,
            3.5e0,
        )
        p2 = (
            0.0e0,
            0.05e0,
            0.5e0,
            1.4e0,
            2.6e0,
            3.7e0,
            4.6e0,
            5.3e0,
            5.95e0,
            6.55e0,
            7.1e0,
        )
        p3 = (
            0.0e0,
            0.05e0,
            0.5e0,
            1.4e0,
            2.6e0,
            3.7e0,
            4.6e0,
            5.4e0,
            6.05e0,
            6.8e0,
            7.2e0,
        )

        #  Dump voltage

        vd = 2.0e0 * tfes / (tdump * aio)

        #  Current density limited by temperature rise during quench

        tav = 1.0e0 + (tmax - tba) / 20.0e0
        no = int(tav)
        np = no + 1
        np = min(np, 11)

        ai1 = 1.0e16 * (p1[no - 1] + (p1[np - 1] - p1[no - 1]) * (tav - no))
        ai2 = 1.0e16 * (p2[no - 1] + (p2[np - 1] - p2[no - 1]) * (tav - no))
        ai3 = 1.0e16 * (p3[no - 1] + (p3[np - 1] - p3[no - 1]) * (tav - no))

        aa = vd * aio / tfes
        bb = (1.0e0 - fcond) * fcond * fcu * ai1
        cc = (fcu * fcond) ** 2 * ai2
        dd = (1.0e0 - fcu) * fcu * fcond**2 * ai3
        ajcp = numpy.sqrt(aa * (bb + cc + dd))
        ajwpro = ajcp * (acs / aturn)

        return ajwpro, vd

    def tf_global_geometry(self):
        """Subroutine for calculating the TF coil geometry
        This includes:
        - Overall geometry of coil (radii and toroidal planes area)
        - Winding Pack NOT included
        """

        sctfcoil_module.theta_coil = numpy.pi / tfcoil_variables.n_tf
        sctfcoil_module.tan_theta_coil = numpy.tan(sctfcoil_module.theta_coil)

        # TF coil inboard legs mid-plane cross-section area (WP + casing ) [m2]
        if tfcoil_variables.i_tf_case_geom == 0:
            # Circular front case
            tfcoil_variables.tfareain = numpy.pi * (
                build_variables.r_tf_inboard_out**2
                - build_variables.r_tf_inboard_in**2
            )
        else:
            # Straight front case
            tfcoil_variables.tfareain = (
                tfcoil_variables.n_tf
                * numpy.sin(sctfcoil_module.theta_coil)
                * numpy.cos(sctfcoil_module.theta_coil)
                * build_variables.r_tf_inboard_out**2
                - numpy.pi * build_variables.r_tf_inboard_in**2
            )

        # Vertical distance from the midplane to the top of the tapered section [m]
        if physics_variables.itart == 1:
            sctfcoil_module.h_cp_top = (
                physics_variables.rminor * physics_variables.kappa
                + tfcoil_variables.dztop
            )
        # ---

        # Outer leg geometry
        # ---
        # Mid-plane inner/out radial position of the TF coil outer leg [m]

        sctfcoil_module.r_tf_outboard_in = (
            build_variables.r_tf_outboard_mid - build_variables.tfthko * 0.5e0
        )
        sctfcoil_module.r_tf_outboard_out = (
            build_variables.r_tf_outboard_mid + build_variables.tfthko * 0.5e0
        )

        # TF coil width in toroidal direction at inboard leg outer edge [m]
        # ***
        # Sliding joints geometry
        if physics_variables.itart == 1 and tfcoil_variables.i_tf_sup != 1:
            tfcoil_variables.tftort = (
                2.0e0 * build_variables.r_cp_top * numpy.sin(sctfcoil_module.theta_coil)
            )

        # Default thickness, initially written for DEMO SC magnets
        elif physics_variables.itart == 1 and tfcoil_variables.i_tf_sup == 1:
            tfcoil_variables.tftort = (
                2.0e0
                * build_variables.r_tf_inboard_out
                * numpy.sin(sctfcoil_module.theta_coil)
            )
        else:
            tfcoil_variables.tftort = (
                2.0e0
                * build_variables.r_tf_inboard_out
                * numpy.sin(sctfcoil_module.theta_coil)
            )

        # Area of rectangular cross-section TF outboard leg [m2]
        tfcoil_variables.arealeg = tfcoil_variables.tftort * build_variables.tfthko
        # ---

    def tf_current(self):
        """
        Calculation of the maximum B field and the corresponding TF current
        """
        if tfcoil_variables.casthi_is_fraction:
            tfcoil_variables.casthi = (
                tfcoil_variables.casthi_fraction * build_variables.tfcth
            )

        # Case thickness of side wall [m]
        if tfcoil_variables.tfc_sidewall_is_fraction:
            tfcoil_variables.casths = (
                tfcoil_variables.casths_fraction
                * (build_variables.r_tf_inboard_in + tfcoil_variables.thkcas)
                * numpy.tan(numpy.pi / tfcoil_variables.n_tf)
            )

        # Radial position of peak toroidal field [m]
        if tfcoil_variables.i_tf_sup == 1:
            # SC : conservative assumption as the radius is calculated with the
            # WP radial distances defined at the TF middle (cos)
            tfcoil_variables.rbmax = (
                build_variables.r_tf_inboard_out * numpy.cos(sctfcoil_module.theta_coil)
                - tfcoil_variables.casthi
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap
            )
        else:
            # Resistive coils : No approx necessary as the symmetry is cylindrical
            # The turn insulation th (tfcoil_variables.thicndut) is also subtracted too here
            tfcoil_variables.rbmax = (
                build_variables.r_tf_inboard_out
                - tfcoil_variables.casthi
                - tfcoil_variables.thicndut
                - tfcoil_variables.tinstf
            )

        # Calculation of the maximum B field on the magnet [T]
        tfcoil_variables.bmaxtf = (
            physics_variables.bt * physics_variables.rmajor / tfcoil_variables.rbmax
        )

        # Total current in TF coils [A]
        # rem SK : ritcf is no longer an input
        tfcoil_variables.ritfc = (
            tfcoil_variables.bmaxtf * tfcoil_variables.rbmax * 5.0e6
        )

        # Current per TF coil [A]
        sctfcoil_module.tfc_current = tfcoil_variables.ritfc / tfcoil_variables.n_tf

        # Global inboard leg average current in TF coils [A/m2]
        tfcoil_variables.oacdcp = tfcoil_variables.ritfc / tfcoil_variables.tfareain

    def coilshap(self):
        """Calculates the TF coil shape
        Calculates the shape of the INSIDE of the TF coil. The coil is
        approximated by a straight inboard section and four elliptical arcs
        This is a totally ad hoc model, with no physics or engineering basis.

        The referenced equations can be found in draft/unpublished document
        attached in GitLab to issue #1328.
        """
        FSTRAIGHT = 0.6
        if tfcoil_variables.i_tf_shape == 1 and physics_variables.itart == 0:
            # PROCESS D-shape parametrisation

            # X position of the arcs, eq(21)
            # The tfcoil_variables.xarc/tfcoil_variables.yarc are defined in the INSIDE part of the TF
            tfcoil_variables.xarc[0] = build_variables.r_tf_inboard_out
            tfcoil_variables.xarc[1] = (
                physics_variables.rmajor - 0.2e0 * physics_variables.rminor
            )
            tfcoil_variables.xarc[2] = sctfcoil_module.r_tf_outboard_in
            tfcoil_variables.xarc[3] = tfcoil_variables.xarc[1]
            tfcoil_variables.xarc[4] = tfcoil_variables.xarc[0]

            # Height of straight section as a fraction of the coil inner height
            if physics_variables.i_single_null == 0:
                # Double null
                tfcoil_variables.yarc[0] = FSTRAIGHT * build_variables.hmax
                tfcoil_variables.yarc[1] = build_variables.hmax
                tfcoil_variables.yarc[2] = 0
                tfcoil_variables.yarc[3] = -build_variables.hmax
                tfcoil_variables.yarc[4] = -FSTRAIGHT * build_variables.hmax
            else:
                # Single null
                tfcoil_variables.yarc[0] = FSTRAIGHT * (
                    build_variables.hpfu - build_variables.tfcth
                )
                tfcoil_variables.yarc[1] = build_variables.hpfu - build_variables.tfcth
                tfcoil_variables.yarc[2] = 0
                tfcoil_variables.yarc[3] = -build_variables.hmax
                tfcoil_variables.yarc[4] = -FSTRAIGHT * build_variables.hmax

            # Horizontal and vertical radii of inside edge of TF coil
            # Arcs are numbered clockwise:
            # 1=upper inboard, 2=upper outboard, 3=lower ouboard, 4=lower inboard
            # 'tfleng' is the length of the coil midline.
            tfcoil_variables.tfleng = (
                tfcoil_variables.yarc[0] - tfcoil_variables.yarc[4]
            )
            for ii in range(4):
                tfcoil_variables.tfa[ii] = abs(
                    tfcoil_variables.xarc[ii + 1] - tfcoil_variables.xarc[ii]
                )
                tfcoil_variables.tfb[ii] = abs(
                    tfcoil_variables.yarc[ii + 1] - tfcoil_variables.yarc[ii]
                )
                # Radii and length of midline of coil segments
                aa = tfcoil_variables.tfa[ii] + 0.5e0 * build_variables.tfcth
                bb = tfcoil_variables.tfb[ii] + 0.5e0 * build_variables.tfcth
                tfcoil_variables.tfleng = (
                    tfcoil_variables.tfleng + 0.25e0 * self.circumference(aa, bb)
                )
                # note: final tfcoil_variables.tfleng includes inboard leg length; eq(22)

        # Centrepost with D-shaped
        # ---
        elif tfcoil_variables.i_tf_shape == 1 and physics_variables.itart == 1:
            # X position of the arcs, eq(23) and text before it
            tfcoil_variables.yarc[0] = build_variables.r_cp_top
            tfcoil_variables.yarc[1] = (
                physics_variables.rmajor - 0.2e0 * physics_variables.rminor
            )
            tfcoil_variables.yarc[2] = sctfcoil_module.r_tf_outboard_in
            tfcoil_variables.yarc[3] = tfcoil_variables.xarc[1]
            tfcoil_variables.yarc[4] = tfcoil_variables.xarc[0]

            # Double null, eq(23) and text before it
            tfcoil_variables.yarc[0] = build_variables.hpfu - build_variables.tfcth
            tfcoil_variables.yarc[1] = build_variables.hpfu - build_variables.tfcth
            tfcoil_variables.yarc[2] = 0
            tfcoil_variables.yarc[3] = -build_variables.hmax
            tfcoil_variables.yarc[4] = -build_variables.hmax

            # TF middle circumference
            tfcoil_variables.tfleng = 2 * (
                tfcoil_variables.xarc[1] - tfcoil_variables.xarc[0]
            )

            for ii in range(1, 3):
                tfcoil_variables.tfa[ii] = abs(
                    tfcoil_variables.xarc[ii + 1] - tfcoil_variables.xarc[ii]
                )
                tfcoil_variables.tfb[ii] = abs(
                    tfcoil_variables.yarc[ii + 1] - tfcoil_variables.yarc[ii]
                )

                # Radii and length of midline of coil segments
                aa = tfcoil_variables.tfa[ii] + 0.5e0 * build_variables.tfthko
                bb = tfcoil_variables.tfb[ii] + 0.5e0 * build_variables.tfthko
                tfcoil_variables.tfleng = (
                    tfcoil_variables.tfleng + 0.25e0 * self.circumference(aa, bb)
                )
                # IMPORTANT : THE CENTREPOST LENGTH IS NOT INCLUDED IN TFLENG FOR TART; eq(24)
        # ---

        # Picture frame coil
        # ---
        elif tfcoil_variables.i_tf_shape == 2:
            # X position of the arcs
            if physics_variables.itart == 0:
                tfcoil_variables.xarc[0] = build_variables.r_tf_inboard_out
            if physics_variables.itart == 1:
                tfcoil_variables.xarc[0] = build_variables.r_cp_top
            tfcoil_variables.xarc[1] = sctfcoil_module.r_tf_outboard_in
            tfcoil_variables.xarc[2] = tfcoil_variables.xarc[1]
            tfcoil_variables.xarc[3] = tfcoil_variables.xarc[1]
            tfcoil_variables.xarc[4] = tfcoil_variables.xarc[0]

            # Y position of the arcs
            tfcoil_variables.yarc[0] = build_variables.hpfu - build_variables.tfcth
            tfcoil_variables.yarc[1] = build_variables.hpfu - build_variables.tfcth
            tfcoil_variables.yarc[2] = 0
            tfcoil_variables.yarc[3] = -build_variables.hmax
            tfcoil_variables.yarc[4] = -build_variables.hmax

            # TF middle circumference
            # IMPORTANT : THE CENTREPOST LENGTH IS NOT INCLUDED IN TFLENG FOR TART
            if physics_variables.itart == 0:
                tfcoil_variables.tfleng = 2.0e0 * (
                    2.0e0 * build_variables.hmax
                    + build_variables.tfcth
                    + build_variables.r_tf_outboard_mid
                    - build_variables.r_tf_inboard_mid
                )  # eq(25)
            elif physics_variables.itart == 1:
                tfcoil_variables.tfleng = (
                    build_variables.hmax
                    + build_variables.hpfu
                    + 2.0e0
                    * (build_variables.r_tf_outboard_mid - build_variables.r_cp_top)
                )  # eq(26)

        # ---

    @staticmethod
    def circumference(aaa, bbb):
        """Calculate ellipse arc circumference using Ramanujan approximation (m)
        See https://www.johndcook.com/blog/2013/05/05/ramanujan-circumference-ellipse/
        for a discussion of the precision of the formula

        An ellipse has the following formula: (x/a)² + (y/b)² = 1

        :param aaa: the value of a in the formula of the ellipse.
        :type aaa: float

        :param bbb: the value of b in the formula of the ellipse.
        :type bbb: float

        :returns: an approximation of the circumference of the ellipse
        :rtype: float
        """
        hh = (aaa - bbb) ** 2 / (aaa + bbb) ** 2
        return (
            numpy.pi
            * (aaa + bbb)
            * (1.0e0 + (3.0e0 * hh) / (10.0e0 + numpy.sqrt(4.0e0 - 3.0e0 * hh)))
        )

    def tf_res_heating(self):
        """
        Resitive magnet resitive heating calculations
        Rem SK : Clamped joined superconductors might have resistive power losses on the joints
        Rem SK : Sliding joints might have a region of high resistivity
        """
        if tfcoil_variables.i_tf_sup == 0:
            tfcoil_variables.rhocp = (
                (tfcoil_variables.frhocp / 0.92e0)
                * (1.72e0 + 0.0039e0 * (tfcoil_variables.tcpav - 273.15e0))
                * 1.0e-8
            )

        # Aluminium
        if tfcoil_variables.i_tf_sup == 2:
            tfcoil_variables.rhocp = tfcoil_variables.frhocp * (
                2.00016e-14 * tfcoil_variables.tcpav**3
                - 6.75384e-13 * tfcoil_variables.tcpav**2
                + 8.89159e-12 * tfcoil_variables.tcpav
            )

        # Calculations dedicated for configurations with CP
        if physics_variables.itart == 1:
            # Tricky trick to make the leg / CP tempearture the same
            if (
                abs(tfcoil_variables.tlegav + 1.0e0)
                < numpy.finfo(float(tfcoil_variables.tlegav)).eps
            ):
                sctfcoil_module.is_leg_cp_temp_same = 1
                tfcoil_variables.tlegav = tfcoil_variables.tcpav

            # Leg resistivity (different leg temperature as separate cooling channels)
            if tfcoil_variables.i_tf_sup == 0:
                tfcoil_variables.rhotfleg = (
                    (tfcoil_variables.frholeg / 0.92e0)
                    * (1.72e0 + 0.0039e0 * (tfcoil_variables.tlegav - 273.15e0))
                    * 1.0e-8
                )
            elif tfcoil_variables.i_tf_sup == 2:
                tfcoil_variables.rhotfleg = tfcoil_variables.frholeg * (
                    2.00016e-14 * tfcoil_variables.tlegav**3
                    - 6.75384e-13 * tfcoil_variables.tlegav**2
                    + 8.89159e-12 * tfcoil_variables.tlegav
                )

            # Tricky trick to make the leg / CP tempearture the same
            if sctfcoil_module.is_leg_cp_temp_same == 1:
                tfcoil_variables.tlegav = -1.0e0

            # Centrepost resisitivity and conductor/insulation volume

            (
                tfcoil_variables.a_cp_cool,
                tfcoil_variables.vol_cond_cp,
                tfcoil_variables.prescp,
                sctfcoil_module.vol_ins_cp,
                sctfcoil_module.vol_case_cp,
                sctfcoil_module.vol_gr_ins_cp,
            ) = self.cpost(
                build_variables.r_tf_inboard_in,
                build_variables.r_tf_inboard_out,
                build_variables.r_cp_top,
                sctfcoil_module.h_cp_top,
                build_variables.hmax + build_variables.tfthko,
                tfcoil_variables.thkcas,
                tfcoil_variables.casthi,
                tfcoil_variables.tinstf,
                tfcoil_variables.thicndut,
                tfcoil_variables.n_tf_turn,
                tfcoil_variables.ritfc,
                tfcoil_variables.rhocp,
                tfcoil_variables.fcoolcp,
                tfcoil_variables.n_tf,
            )

        # Leg cross-section areas
        # Rem : For physics_variables.itart = 1, these quantitire corresponds to the outer leg only
        # ---
        # Leg ground insulation area per coil [m2]
        sctfcoil_module.a_leg_gr_ins = tfcoil_variables.arealeg - (
            tfcoil_variables.tftort - 2.0e0 * tfcoil_variables.tinstf
        ) * (build_variables.tfthko - 2.0e0 * tfcoil_variables.tinstf)

        # Outboard leg turns insulation area per coil [m2]
        sctfcoil_module.a_leg_ins = 2.0e0 * tfcoil_variables.thicndut * (
            tfcoil_variables.tftort - 2.0e0 * tfcoil_variables.tinstf
        ) + 2.0e0 * tfcoil_variables.thicndut * tfcoil_variables.n_tf_turn * (
            build_variables.tfthko
            - 2.0e0 * (tfcoil_variables.thicndut + tfcoil_variables.tinstf)
        )  # toroidal direction + radial direction

        # Exact TF outboard leg conductor area per coil [m2]
        sctfcoil_module.a_leg_cond = (1.0e0 - tfcoil_variables.fcoolleg) * (
            tfcoil_variables.arealeg
            - sctfcoil_module.a_leg_gr_ins
            - sctfcoil_module.a_leg_ins
        )
        # ---

        if physics_variables.itart == 1:
            # Outer leg resistive power loss
            # ---
            # TF outboard leg's resistance calculation (per leg) [ohm]
            tfcoil_variables.tflegres = (
                tfcoil_variables.rhotfleg
                * tfcoil_variables.tfleng
                / sctfcoil_module.a_leg_cond
            )

            # TF outer leg resistive power (TOTAL) [W]
            tfcoil_variables.presleg = (
                tfcoil_variables.tflegres
                * tfcoil_variables.ritfc**2
                / tfcoil_variables.n_tf
            )
            # ---

            # Sliding joints resistive heating
            # ---
            if tfcoil_variables.i_cp_joints != 0:
                # Number of contact area per joint (all legs)
                n_contact_tot = (
                    tfcoil_variables.n_tf_joints_contact
                    * numpy.round(tfcoil_variables.n_tf_turn)
                    * numpy.round(tfcoil_variables.n_tf)
                )

                # Area of joint contact (all legs)
                a_joints = (
                    build_variables.tfthko
                    * tfcoil_variables.th_joint_contact
                    * n_contact_tot
                )

                # Total joints resistive power losses
                tfcoil_variables.pres_joints = (
                    tfcoil_variables.n_tf_joints
                    * tfcoil_variables.rho_tf_joints
                    * tfcoil_variables.ritfc**2
                    / a_joints
                )
            else:
                # Joints resistance to be evaluated for SC
                tfcoil_variables.pres_joints = 0.0e0

            # ---

        # Case of a resistive magnet without joints
        # ***
        else:
            # TF resistive powers
            tfcoil_variables.prescp = (
                tfcoil_variables.rhocp
                * tfcoil_variables.ritfc**2
                * tfcoil_variables.tfleng
                / (sctfcoil_module.a_leg_cond * tfcoil_variables.n_tf)
            )

            # tfcoil_variables.prescp containts the the total resistive power losses
            tfcoil_variables.presleg = 0.0e0

            # No joints if physics_variables.itart = 0
            tfcoil_variables.pres_joints = 0.0e0

    @staticmethod
    @numba.njit(cache=True)
    def cpost(
        r_tf_inboard_in,
        r_tf_inboard_out,
        r_cp_top,
        ztop,
        hmaxi,
        cas_in_th,
        cas_out_th,
        gr_ins_th,
        ins_th,
        n_tf_turn,
        curr,
        rho,
        fcool,
        n_tf,
    ):
        """
        author: P J Knight, CCFE, Culham Science Centre
        Calculates the volume and resistive power losses of a TART centrepost
        This routine calculates the volume and resistive power losses
        of a TART centrepost. It is assumed to be tapered - narrowest at
        the midplane and reaching maximum thickness at the height of the
        plasma. Above/below the plasma, the centrepost is cylindrical.
        The shape of the taper is assumed to be an arc of a circle.
        P J Knight, CCFE, Culham Science Centre
        21/10/96 PJK Initial version
        08/05/12 PJK Initial F90 version
        16/10/12 PJK Added constants; removed argument pi
        26/06/14 PJK Added error handling
        12/11/19 SK Using fixed cooling cross-section area along the CP
        26/11/19 SK added the coolant area, the conuctor/isulator/outer casing volume
        30/11/20 SK added the ground outer ground insulation volume
        F/MI/PJK/LOGBOOK12, pp.33,34
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        yy_ins = numpy.zeros((101,))  # Exact conductor area (to be integrated)
        yy_cond = numpy.zeros((101,))  # Turn insulation area (to be integrated)
        yy_gr_ins = numpy.zeros(
            (101,)
        )  # Outter ground insulation area (to be integrated)
        yy_casout = numpy.zeros((101,))  # Outter case area (to be integrated)

        rtop = r_cp_top - cas_out_th - gr_ins_th

        # Conductor outer radius at CP mid-plane [m]
        rmid = r_tf_inboard_out - cas_out_th - gr_ins_th

        # Conductor inner radius [m]
        r_tfin_inleg = r_tf_inboard_in + cas_in_th + gr_ins_th
        # -#

        #  Error traps
        # ------------
        # if rtop <= 0.0e0:
        #     error_handling.fdiags[0] = rtop
        #     error_handling.report_error(115)

        # if ztop <= 0.0e0:
        #     error_handling.fdiags[0] = ztop
        #     error_handling.report_error(116)

        # if rmid <= 0.0e0:
        #     error_handling.fdiags[0] = rmid
        #     error_handling.report_error(117)

        # if build_variables.hmax <= 0.0e0:
        #     error_handling.fdiags[0] = build_variables.hmax
        #     error_handling.report_error(118)

        # if (fcool < 0.0e0) or (fcool > 1.0e0):
        #     error_handling.fdiags[0] = fcool
        #     error_handling.report_error(119)

        # if rtop < rmid:
        #     error_handling.fdiags[0] = rtop
        #     error_handling.fdiags[1] = rmid
        #     error_handling.report_error(120)

        # if build_variables.hmax < ztop:
        #     error_handling.fdiags[0] = build_variables.hmax
        #     error_handling.fdiags[1] = ztop
        #     error_handling.report_error(121)

        # ------------

        # Mid-plane area calculations
        # ---------------------------
        # Total number of CP turns
        n_turns_tot = n_tf * n_tf_turn

        # Area of the innner TF central hole [m2]
        a_tfin_hole = numpy.pi * r_tfin_inleg**2

        # Mid-plane outer casing cross-section area [m2]
        a_casout = numpy.pi * (
            (rmid + gr_ins_th + cas_out_th) ** 2 - (rmid + gr_ins_th) ** 2
        )

        # Mid-plane outter ground insulation thickness [m2]
        a_cp_gr_ins = (
            numpy.pi * ((rmid + gr_ins_th) ** 2 - rmid**2)
            + 2.0e0 * gr_ins_th * (rmid - r_tfin_inleg) * n_tf
        )

        # Mid-plane turn layer cross-section area [m2]
        a_cp_ins = (
            numpy.pi
            * ((r_tfin_inleg + ins_th) ** 2 - r_tfin_inleg**2)  # Inner layer volume
            + numpy.pi * (rmid**2 - (rmid - ins_th) ** 2)  # Outter layer volume
            + 2.0e0 * n_turns_tot * ins_th * (rmid - r_tfin_inleg - 2.0e0 * ins_th)
        )  # inter turn separtion

        # Cooling pipes cross-section per coil [m2]
        a_cp_cool = fcool * (
            (numpy.pi * rmid**2 - a_tfin_hole - a_cp_ins) / n_tf
            - 2.0e0 * gr_ins_th * (rmid - r_tfin_inleg)
        )  # Wedge ground insulation
        # ---------------------------

        #  Trivial solutions
        # ------------------
        if numpy.abs(fcool) < EPS:
            vol_cond_cp = 0.0e0
            respow = 0.0e0
            vol_case_cp = 0.0e0
            vol_gr_ins_cp = 0.0e0
            vol_ins_cp = 0.0e0
            # error_handling.report_error(122)
            return (
                a_cp_cool,
                vol_cond_cp,
                respow,
                vol_ins_cp,
                vol_case_cp,
                vol_gr_ins_cp,
            )

        if numpy.abs(rmid - rtop) < EPS:
            # Exact conductor cross-section
            a_cond_midplane = (
                numpy.pi * rmid**2 - a_tfin_hole - n_tf * a_cp_cool - a_cp_ins
            )

            # Volumes and resisitive losses calculations
            vol_cond_cp = 2.0e0 * hmaxi * a_cond_midplane
            vol_ins_cp = 2.0e0 * hmaxi * a_cp_ins
            vol_gr_ins_cp = 2.0e0 * hmaxi * a_cp_gr_ins
            respow = 2.0e0 * hmaxi * curr**2 * rho / a_cond_midplane
            vol_case_cp = 2.0e0 * hmaxi * a_casout

            return (
                a_cp_cool,
                vol_cond_cp,
                respow,
                vol_ins_cp,
                vol_case_cp,
                vol_gr_ins_cp,
            )

        # ------------------

        # Find centre of circle (RC,0) defining the taper's arc
        # (r1,z1) is midpoint of line joining (rmid,0) and (rtop,ztop)
        # Rem : The taper arc is defined using the outer radius of the
        #       conductor including turn unsulation
        # -------------------------------------------------------------
        r1 = 0.5e0 * (rmid + rtop)
        z1 = 0.5e0 * ztop

        x = (r1 - rmid) ** 2 + z1**2
        y = ztop**2 / ((rtop - rmid) ** 2 + ztop**2)

        rc = rmid + numpy.sqrt(x / (1.0e0 - y))
        # -------------------------------------------------------------

        #  Find volume of tapered section of centrepost, and the resistive
        #  power losses, by integrating along the centrepost from the midplane
        # --------------------------------------------------------------------
        #  Calculate centrepost radius and cross-sectional areas at each Z
        dz = 0.01e0 * ztop

        for ii in range(101):
            z = ii * dz
            z = numpy.fmin(numpy.array(z), ztop)

            r = rc - numpy.sqrt((rc - rmid) ** 2 - z * z)

            # if r <= 0.0e0:
            #     error_handling.fdiags[0] = r
            #     error_handling.fdiags[1] = rc
            #     error_handling.fdiags[2] = rmid
            #     error_handling.fdiags[3] = z

            #     error_handling.report_error(123)

            # Insulation cross-sectional area at z
            yy_ins[ii] = (
                numpy.pi * ((r_tfin_inleg + ins_th) ** 2 - r_tfin_inleg**2)
                + numpy.pi * (r**2 - (r - ins_th) ** 2)  # Inner layer volume
                + 2.0e0  # Outter layer volume
                * ins_th
                * (r - r_tfin_inleg - 2.0e0 * ins_th)
                * n_turns_tot
            )  # inter turn layers

            #  Conductor cross-sectional area at z
            yy_cond[ii] = (
                numpy.pi * r**2
                - a_tfin_hole
                - n_tf * a_cp_cool
                - yy_ins[ii]
                - 2.0e0 * n_tf * gr_ins_th * (r - r_tfin_inleg)
            )  # Wedge ground insulation

            #  Outer ground insulation area at z
            yy_gr_ins[ii] = numpy.pi * (
                (r + gr_ins_th) ** 2 - r**2
            ) + 2.0e0 * n_tf * gr_ins_th * (r - r_tfin_inleg)

            #  Outer casing Cross-sectional area at z
            yy_casout[ii] = numpy.pi * (
                (r + gr_ins_th + cas_out_th) ** 2 - (r + gr_ins_th) ** 2
            )

        #  Perform integrals using trapezium rule
        sum1 = 0.0e0
        sum2 = 0.0e0
        sum3 = 0.0e0
        sum4 = 0.0e0
        sum5 = 0.0e0
        for ii in range(1, 100):
            sum1 = sum1 + yy_cond[ii]
            sum2 = sum2 + 1.0e0 / yy_cond[ii]
            sum3 = sum3 + yy_ins[ii]
            sum4 = sum4 + yy_casout[ii]
            sum5 = sum5 + yy_gr_ins[ii]

        sum1 = 0.5e0 * dz * (yy_cond[0] + yy_cond[100] + 2.0e0 * sum1)
        sum2 = 0.5e0 * dz * (1.0e0 / yy_cond[0] + 1.0e0 / yy_cond[100] + 2.0e0 * sum2)
        sum3 = 0.5e0 * dz * (yy_ins[0] + yy_ins[100] + 2.0e0 * sum3)
        sum4 = 0.5e0 * dz * (yy_casout[0] + yy_casout[100] + 2.0e0 * sum4)
        sum5 = 0.5e0 * dz * (yy_gr_ins[0] + yy_gr_ins[100] + 2.0e0 * sum5)

        # Turn insulation layer cross section at CP top  [m2]
        a_cp_ins = (
            numpy.pi * ((r_tfin_inleg + ins_th) ** 2 - r_tfin_inleg**2)
            + numpy.pi * (rtop**2 - (rtop - ins_th) ** 2)  # Inner layer volume
            + 2.0e0  # Outter layer volume
            * ins_th
            * (rtop - r_tfin_inleg - 2.0e0 * ins_th)
            * n_turns_tot
        )  # turn separtion layers

        # Ground insulation layer cross-section at CP top [m2]
        a_cp_gr_ins = (
            numpy.pi * ((rtop + gr_ins_th) ** 2 - rtop**2)
            + 2.0e0 * gr_ins_th * (rtop - r_tfin_inleg) * n_tf
        )

        # Outer casing cross-section area at CP top [m2]
        a_casout = numpy.pi * (
            (rmid + gr_ins_th + cas_out_th) ** 2 - (rmid + gr_ins_th) ** 2
        )

        # Centrepost volume (ignoring coolant fraction) [m3]
        vol_cond_cp = 2.0e0 * sum1 + 2.0e0 * (  # Tapered section
            hmaxi - ztop
        ) * (  # Straight section vertical height
            numpy.pi * rtop**2
            - a_tfin_hole
            - a_cp_ins
            - n_tf * a_cp_cool
            - 2.0e0 * n_tf * gr_ins_th * (rtop - r_tfin_inleg)
        )  # subtracting ground insulation wedge separation

        # Resistive power losses in taped section (variable radius section) [W]
        res_taped = rho * curr**2 * sum2

        # Centrepost insulator volume [m3]
        vol_ins_cp = 2.0e0 * (sum3 + (hmaxi - ztop) * a_cp_ins)

        # Ground insulation volume [m3]
        vol_gr_ins_cp = 2.0e0 * (
            sum5
            + (hmaxi - ztop) * a_cp_gr_ins
            + hmaxi * numpy.pi * (r_tfin_inleg**2 - (r_tfin_inleg - gr_ins_th) ** 2)
        )

        # CP casing volume [m3]
        vol_case_cp = 2.0e0 * (
            sum4
            + (hmaxi - ztop) * a_casout
            + hmaxi
            * numpy.pi
            * (
                (r_tfin_inleg - gr_ins_th) ** 2
                - (r_tfin_inleg - gr_ins_th - cas_in_th) ** 2
            )
        )

        # Resistive power losses in cylindrical section (constant radius) [W]
        res_cyl = (
            rho
            * curr**2
            * (
                (hmaxi - ztop)
                / (
                    numpy.pi * rtop**2
                    - a_tfin_hole
                    - a_cp_ins
                    - n_tf * a_cp_cool
                    - 2.0e0 * n_tf * gr_ins_th * (rtop - r_tfin_inleg)
                )
            )
        )  # ground insulation separation

        # Total CP resistive power [W]
        respow = 2.0e0 * (res_cyl + res_taped)

        return (
            a_cp_cool,
            vol_cond_cp,
            respow,
            vol_ins_cp,
            vol_case_cp,
            vol_gr_ins_cp,
        )

    def vv_stress_on_quench(self):
        """Calculate the Tresca stress [Pa] of the Vacuum Vessel (VV)
        experienced when the TF coil quenches.

        Author: Timothy Nunn, UKAEA

        Assumes the current center line (CCL) of the TF coil is the
        middle of the coil.

        We assume vertical symmetry which is only true for double null
        machines.
        """
        H_coil = build_variables.hmax + (build_variables.tfcth / 2)
        Ri_coil = build_variables.r_tf_inboard_mid
        Ro_coil = build_variables.r_tf_outboard_mid
        # NOTE: Rm is measured from the outside edge of the coil because thats where
        # the radius of the first ellipse is measured from
        Rm_coil = build_variables.r_tf_inboard_out + tfcoil_variables.tfa[0]

        H_vv = (
            physics_variables.rminor * physics_variables.kappa
            + build_variables.vgap
            + divertor_variables.divfix
            + build_variables.shldtth
            + (build_variables.d_vv_top / 2)
        )
        # Ri and Ro for VV dont consider the shield widths
        # because it is assumed the shield is on the plasma side
        # of the VV
        Ri_vv = build_variables.r_vv_inboard_out - (build_variables.d_vv_out / 2)
        Ro_vv = (
            build_variables.r_tf_outboard_mid
            - (build_variables.tfthko / 2)
            - build_variables.tftsgap
            - build_variables.thshield_ob
            - build_variables.gapsto
            - (build_variables.d_vv_out / 2)
        )

        # Assume the radius of the first ellipse of the VV is in the same proportion to
        # that of the plasma facing radii of the two structures
        tf_vv_frac = build_variables.r_tf_inboard_out / build_variables.r_vv_inboard_out
        Rm_vv = build_variables.r_vv_inboard_out + (
            tfcoil_variables.tfa[0] * tf_vv_frac
        )

        sctfcoil_module.vv_stress_quench = vv_stress_on_quench(
            # TF shape
            H_coil=H_coil,
            Ri_coil=Ri_coil,
            Ro_coil=Ro_coil,
            Rm_coil=Rm_coil,
            ccl_length_coil=tfcoil_variables.tfleng,
            theta1_coil=tfcoil_variables.theta1_coil,
            # VV shape
            H_vv=H_vv,
            Ri_vv=Ri_vv,
            Ro_vv=Ro_vv,
            Rm_vv=Rm_vv,
            theta1_vv=tfcoil_variables.theta1_vv,
            # TF properties
            n_tf=tfcoil_variables.n_tf,
            n_tf_turn=tfcoil_variables.n_tf_turn,
            # Area of the radial plate taken to be the area of steel in the WP
            # TODO: value clipped due to #1883
            S_rp=numpy.clip(sctfcoil_module.a_tf_steel, 0, None),
            # TODO: Does this calculation of Scc exclude the area of the case down the side?
            S_cc=sctfcoil_module.a_case_front + sctfcoil_module.a_case_nose,
            taud=tfcoil_variables.tdmptf,
            # TODO: is this the correct current?
            I_op=sctfcoil_module.tfc_current / tfcoil_variables.n_tf_turn,
            # VV properties
            d_vv=build_variables.d_vv_in,
        )

    def tf_field_and_force(self):
        """
        Calculate the TF coil field, force and VV quench consideration, and the resistive magnets resistance/volume
        """
        if tfcoil_variables.i_tf_sup == 1:
            self.vv_stress_on_quench()

        # Outer/inner WP radius removing the ground insulation layer and the insertion gap [m]
        if tfcoil_variables.i_tf_sup == 1:
            r_out_wp = (
                sctfcoil_module.r_wp_outer
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap
            )
            r_in_wp = (
                sctfcoil_module.r_wp_inner
                + tfcoil_variables.tinstf
                + tfcoil_variables.tfinsgap
            )
        else:
            r_out_wp = sctfcoil_module.r_wp_outer - tfcoil_variables.tinstf
            r_in_wp = sctfcoil_module.r_wp_inner + tfcoil_variables.tinstf

        # Associated WP thickness
        dr_wp = r_out_wp - r_in_wp

        # In plane forces
        # ---
        # Centering force = net inwards radial force per meters per TF coil [N/m]
        tfcoil_variables.cforce = (
            0.5e0
            * tfcoil_variables.bmaxtf
            * tfcoil_variables.ritfc
            / tfcoil_variables.n_tf
        )

        # Vertical force per coil [N]
        # ***
        # Rem : this force does not depends on the TF shape or the presence of
        #        sliding joints, the in/outboard vertical tension repartition is
        # -#
        # Ouboard leg WP plasma side radius without ground insulation/insertion gat [m]
        if tfcoil_variables.i_tf_sup == 1:
            r_in_outwp = (
                sctfcoil_module.r_tf_outboard_in
                + tfcoil_variables.casthi
                + tfcoil_variables.tinstf
                + tfcoil_variables.tfinsgap
            )
        else:
            r_in_outwp = sctfcoil_module.r_tf_outboard_in + tfcoil_variables.tinstf

        # If the TF coil has no bore it would induce division by 0.
        # In this situation, the bore radius is set to a very small value : 1.0e-9 m
        if abs(r_in_wp) < numpy.finfo(float(r_in_wp)).eps:
            r_in_wp = 1.0e-9

        # May the force be with you
        vforce_tot = (
            0.5e0
            * (physics_variables.bt * physics_variables.rmajor * tfcoil_variables.ritfc)
            / (tfcoil_variables.n_tf * dr_wp**2)
            * (
                r_out_wp**2 * numpy.log(r_out_wp / r_in_wp)
                + r_in_outwp**2 * numpy.log((r_in_outwp + dr_wp) / r_in_outwp)
                + dr_wp**2 * numpy.log((r_in_outwp + dr_wp) / r_in_wp)
                - dr_wp * (r_out_wp + r_in_outwp)
                + 2.0e0
                * dr_wp
                * (
                    r_out_wp * numpy.log(r_in_wp / r_out_wp)
                    + r_in_outwp * numpy.log((r_in_outwp + dr_wp) / r_in_outwp)
                )
            )
        )

        # Case of a centrepost (physics_variables.itart == 1) with sliding joints (the CP vertical are separated from the leg ones)
        # Rem SK : casing/insulation thickness not subtracted as part of the CP is genuinely connected to the legs..
        if physics_variables.itart == 1 and tfcoil_variables.i_cp_joints == 1:
            # CP vertical tension [N]
            tfcoil_variables.vforce = (
                0.25e0
                * (
                    physics_variables.bt
                    * physics_variables.rmajor
                    * tfcoil_variables.ritfc
                )
                / (tfcoil_variables.n_tf * dr_wp**2)
                * (
                    2.0e0 * r_out_wp**2 * numpy.log(r_out_wp / r_in_wp)
                    + 2.0e0 * dr_wp**2 * numpy.log(build_variables.r_cp_top / r_in_wp)
                    + 3.0e0 * dr_wp**2
                    - 2.0e0 * dr_wp * r_out_wp
                    + 4.0e0 * dr_wp * r_out_wp * numpy.log(r_in_wp / r_out_wp)
                )
            )

            # Vertical tension applied on the outer leg [N]
            tfcoil_variables.vforce_outboard = vforce_tot - tfcoil_variables.vforce

            # Inboard vertical tension fraction
            tfcoil_variables.f_vforce_inboard = tfcoil_variables.vforce / vforce_tot

        # Case of TF without joints or with clamped joints vertical tension
        else:
            # Inboard vertical tension [N]
            tfcoil_variables.vforce = tfcoil_variables.f_vforce_inboard * vforce_tot

            # Ouboard vertical tension [N]
            tfcoil_variables.vforce_outboard = tfcoil_variables.vforce * (
                (1.0e0 / tfcoil_variables.f_vforce_inboard) - 1.0e0
            )

        # ***

        # Total vertical force
        sctfcoil_module.vforce_inboard_tot = (
            tfcoil_variables.vforce * tfcoil_variables.n_tf
        )

    @staticmethod
    @numba.njit(cache=True)
    def tfcind(tfthk, xarc, yarc):
        """Calculates the self inductance of a TF coil
        This routine calculates the self inductance of a TF coil
        approximated by a straight inboard section and two elliptical arcs.
        The inductance of the TFC (considered as a single axisymmetric turn)
        is calculated by numerical integration over the cross-sectional area.
        The contribution from the cross-sectional area of the
        coil itself is calculated by taking the field as B(r)/2.
        The field in the bore is calculated for unit current.
        Top/bottom symmetry is assumed.

        :param tfthk: TF coil thickness (m)
        :type tfthk: float
        """
        NINTERVALS = 100

        # Integrate over the whole TF area, including the coil thickness.
        x0 = xarc[1]
        y0 = yarc[1]

        # Minor and major radii of the inside and outside perimeters of the the
        # Inboard leg and arc.
        # Average the upper and lower halves, which are different in the
        # single null case
        ai = xarc[1] - xarc[0]
        bi = (yarc[1] - yarc[3]) / 2.0e0 - yarc[0]
        ao = ai + tfthk
        bo = bi + tfthk
        # Interval used for integration
        dr = ao / NINTERVALS
        # Start both integrals from the centre-point where the arcs join.
        # Initialise major radius
        r = x0 - dr / 2.0e0

        tfind = 0

        for _ in range(NINTERVALS):
            # Field in the bore for unit current
            b = RMU0 / (2.0e0 * numpy.pi * r)
            # Find out if there is a bore
            if x0 - r < ai:
                h_bore = y0 + bi * numpy.sqrt(1 - ((r - x0) / ai) ** 2)
                h_thick = bo * numpy.sqrt(1 - ((r - x0) / ao) ** 2) - h_bore
            else:
                h_bore = 0.0e0
                # Include the contribution from the straight section
                h_thick = bo * numpy.sqrt(1 - ((r - x0) / ao) ** 2) + yarc[0]

            # Assume B in TF coil = 1/2  B in bore
            # Multiply by 2 for upper and lower halves of coil
            tfind += b * dr * (2.0e0 * h_bore + h_thick)
            r = r - dr

        # Outboard arc
        ai = xarc[2] - xarc[1]
        bi = (yarc[1] - yarc[3]) / 2.0e0
        ao = ai + tfthk
        bo = bi + tfthk
        dr = ao / NINTERVALS
        # Initialise major radius
        r = x0 + dr / 2.0e0

        for _ in range(NINTERVALS):
            # Field in the bore for unit current
            b = RMU0 / (2.0e0 * numpy.pi * r)
            # Find out if there is a bore
            if r - x0 < ai:
                h_bore = y0 + bi * numpy.sqrt(1 - ((r - x0) / ai) ** 2)
                h_thick = bo * numpy.sqrt(1 - ((r - x0) / ao) ** 2) - h_bore
            else:
                h_bore = 0.0e0
                h_thick = bo * numpy.sqrt(1 - ((r - x0) / ao) ** 2)

            # Assume B in TF coil = 1/2  B in bore
            # Multiply by 2 for upper and lower halves of coil
            tfind += b * dr * (2.0e0 * h_bore + h_thick)
            r = r + dr

        return tfind

    def tf_coil_area_and_masses(self):
        """Subroutine to calculate the TF coil areas and masses"""
        vol_case = 0.0e0  # Total TF case volume [m3]
        vol_ins = 0.0e0  # Total leg turn insulation volume [m3]
        vol_gr_ins = 0.0e0  # Total leg turn insulation volume [m3]
        vol_cond = 0.0e0  # Total conductor insulator volume [m3]
        vol_ins_leg = 0.0e0  # Outboard leg turn isulation volume [m3]
        vol_gr_ins_leg = 0.0e0  # Outboard leg turn isulation volume [m3]
        vol_cond_leg = 0.0e0  # Outboard leg conductor insulator volume [m3]
        # ---

        # Surface areas (for cryo system) [m2]
        # tfsai, tfcoil_variables.tfsao are retained for the (obsolescent) TF coil nuclear heating calculation
        wbtf = (
            build_variables.r_tf_inboard_out * numpy.sin(sctfcoil_module.theta_coil)
            - build_variables.r_tf_inboard_in * sctfcoil_module.tan_theta_coil
        )
        tfcoil_variables.tfocrn = (
            build_variables.r_tf_inboard_in * sctfcoil_module.tan_theta_coil
        )
        tfcoil_variables.tficrn = tfcoil_variables.tfocrn + wbtf
        tfcoil_variables.tfsai = (
            4.0e0
            * tfcoil_variables.n_tf
            * tfcoil_variables.tficrn
            * build_variables.hr1
        )
        tfcoil_variables.tfsao = (
            2.0e0
            * tfcoil_variables.n_tf
            * tfcoil_variables.tficrn
            * (tfcoil_variables.tfleng - 2.0e0 * build_variables.hr1)
        )

        # Total surface area of two toroidal shells covering the TF coils [m2]
        # (inside and outside surfaces)
        # = 2 * centroid coil length * 2 pi R, where R is average of i/b and o/b centres
        # (This will possibly be used to replace 2*tfcoil_variables.tfsai in the calculation of qss
        # in subroutine cryo - not done at present.)
        tfcoil_variables.tfcryoarea = (
            2.0e0
            * tfcoil_variables.tfleng
            * constants.twopi
            * 0.5e0
            * (build_variables.r_tf_inboard_mid + build_variables.r_tf_outboard_mid)
        )

        # Superconductor coil design specific calculation
        # ---
        if tfcoil_variables.i_tf_sup == 1:
            # Mass of case [kg]
            # ***

            # Mass of ground-wall insulation [kg]
            # (assumed to be same density/material as turn insulation)
            tfcoil_variables.whtgw = (
                tfcoil_variables.tfleng
                * (sctfcoil_module.awpc - sctfcoil_module.awptf)
                * tfcoil_variables.dcondins
            )

            # The length of the vertical section is that of the first (inboard) segment
            # = height of TF coil inner edge + (2 * coil thickness)
            tfcoil_variables.cplen = (2.0e0 * build_variables.hmax) + (
                2.0e0 * build_variables.tfcth
            )

            # The 2.2 factor is used as a scaling factor to fit
            # to the ITER-FDR value of 450 tonnes; see CCFE note T&M/PKNIGHT/PROCESS/026
            if physics_variables.itart == 1:
                # tfcoil_variables.tfleng does not include inboard leg ('centrepost') length in TART
                tfcoil_variables.whtcas = (
                    2.2e0
                    * tfcoil_variables.dcase
                    * (
                        tfcoil_variables.cplen * tfcoil_variables.acasetf
                        + tfcoil_variables.tfleng * tfcoil_variables.acasetfo
                    )
                )
            else:
                tfcoil_variables.whtcas = (
                    2.2e0
                    * tfcoil_variables.dcase
                    * (
                        tfcoil_variables.cplen * tfcoil_variables.acasetf
                        + (tfcoil_variables.tfleng - tfcoil_variables.cplen)
                        * tfcoil_variables.acasetfo
                    )
                )

            # ***

            # Masses of conductor constituents
            # ---------------------------------
            # Superconductor mass [kg]
            # Includes space allowance for central helium channel, area tfcoil_variables.awphec
            tfcoil_variables.whtconsc = (
                tfcoil_variables.tfleng
                * tfcoil_variables.n_tf_turn
                * tfcoil_variables.acstf
                * (1.0e0 - tfcoil_variables.vftf)
                * (1.0e0 - tfcoil_variables.fcutfsu)
                - tfcoil_variables.tfleng * tfcoil_variables.awphec
            ) * tfcoil_variables.dcond[tfcoil_variables.i_tf_sc_mat - 1]

            # Copper mass [kg]
            tfcoil_variables.whtconcu = (
                tfcoil_variables.tfleng
                * tfcoil_variables.n_tf_turn
                * tfcoil_variables.acstf
                * (1.0e0 - tfcoil_variables.vftf)
                * tfcoil_variables.fcutfsu
                - tfcoil_variables.tfleng * tfcoil_variables.awphec
            ) * constants.dcopper
            if tfcoil_variables.whtconcu <= 0.0e0:
                tfcoil_variables.whtconcu = 0.0e0

            # Steel conduit (sheath) mass [kg]
            tfcoil_variables.whtconsh = (
                tfcoil_variables.tfleng
                * tfcoil_variables.n_tf_turn
                * tfcoil_variables.acndttf
                * fwbs_variables.denstl
            )

            # Conduit insulation mass [kg]
            # (tfcoil_variables.aiwp already contains tfcoil_variables.n_tf_turn)
            tfcoil_variables.whtconin = (
                tfcoil_variables.tfleng
                * tfcoil_variables.aiwp
                * tfcoil_variables.dcondins
            )

            # Total conductor mass [kg]
            tfcoil_variables.whtcon = (
                tfcoil_variables.whtconsc
                + tfcoil_variables.whtconcu
                + tfcoil_variables.whtconsh
                + tfcoil_variables.whtconin
            )
            # ---------------------------------

            # Total TF coil mass [kg] (all coils)
            tfcoil_variables.whttf = (
                tfcoil_variables.whtcas
                + tfcoil_variables.whtcon
                + tfcoil_variables.whtgw
            ) * tfcoil_variables.n_tf

            # If spherical tokamak, distribute between centrepost and outboard legs
            # (in this case, total TF coil length = inboard `cplen` + outboard `tfleng`)
            if physics_variables.itart == 1:
                tfleng_sph = tfcoil_variables.cplen + tfcoil_variables.tfleng
                tfcoil_variables.whtcp = tfcoil_variables.whttf * (
                    tfcoil_variables.cplen / tfleng_sph
                )
                tfcoil_variables.whttflgs = tfcoil_variables.whttf * (
                    tfcoil_variables.tfleng / tfleng_sph
                )

        # Resitivive magnets weights
        # ---
        # Rem SK : No casing for the outboard leg is considered for now #
        else:
            # Volumes
            # -------
            # CP with joints
            # ---
            if physics_variables.itart == 1:
                # Total volume of one outerleg [m3]
                tfcoil_variables.voltfleg = (
                    tfcoil_variables.tfleng * tfcoil_variables.arealeg
                )

                # Outboard leg TF conductor volume [m3]
                vol_cond_leg = tfcoil_variables.tfleng * sctfcoil_module.a_leg_cond

                # Total TF conductor volume [m3]
                vol_cond = (
                    tfcoil_variables.vol_cond_cp + tfcoil_variables.n_tf * vol_cond_leg
                )

                # Outboard leg TF turn insulation layer volume (per leg) [m3]
                vol_ins_leg = tfcoil_variables.tfleng * sctfcoil_module.a_leg_ins

                # Total turn insulation layer volume [m3]
                vol_ins = (
                    sctfcoil_module.vol_ins_cp + tfcoil_variables.n_tf * vol_ins_leg
                )

                # Ouboard leg TF ground insulation layer volume (per leg) [m3]
                vol_gr_ins_leg = tfcoil_variables.tfleng * sctfcoil_module.a_leg_gr_ins

                # Total ground insulation layer volume [m3]
                vol_gr_ins = (
                    sctfcoil_module.vol_gr_ins_cp
                    + tfcoil_variables.n_tf * vol_gr_ins_leg
                )

                # Total volume of the CP casing [m3]
                # Rem : no outer leg case
                vol_case = sctfcoil_module.vol_case_cp

            # No joints
            # ---
            else:
                # Total TF outer leg conductor volume [m3]
                vol_cond = (
                    tfcoil_variables.tfleng
                    * sctfcoil_module.a_leg_cond
                    * tfcoil_variables.n_tf
                )

                # Total turn insulation layer volume [m3]
                vol_ins = (
                    tfcoil_variables.tfleng
                    * sctfcoil_module.a_leg_ins
                    * tfcoil_variables.n_tf
                )

                # Total ground insulation volume [m3]
                vol_gr_ins = (
                    tfcoil_variables.tfleng
                    * sctfcoil_module.a_leg_gr_ins
                    * tfcoil_variables.n_tf
                )

                # Total case volume [m3]
                vol_case = (
                    tfcoil_variables.tfleng
                    * tfcoil_variables.acasetf
                    * tfcoil_variables.n_tf
                )

            # ---
            # -------

            # Copper magnets casing/conductor weights per coil [kg]
            if tfcoil_variables.i_tf_sup == 0:
                tfcoil_variables.whtcas = (
                    fwbs_variables.denstl * vol_case / tfcoil_variables.n_tf
                )  # Per TF leg, no casing for outer leg
                tfcoil_variables.whtconcu = (
                    constants.dcopper * vol_cond / tfcoil_variables.n_tf
                )
                tfcoil_variables.whtconal = 0.0e0

                # Outer legs/CP weights
                if physics_variables.itart == 1:
                    # Weight of all the TF legs
                    tfcoil_variables.whttflgs = tfcoil_variables.n_tf * (
                        constants.dcopper * vol_cond_leg
                        + tfcoil_variables.dcondins * (vol_ins_leg + vol_gr_ins_leg)
                    )

                    # CP weight
                    tfcoil_variables.whtcp = (
                        constants.dcopper * tfcoil_variables.vol_cond_cp
                        + tfcoil_variables.dcondins
                        * (sctfcoil_module.vol_ins_cp + sctfcoil_module.vol_gr_ins_cp)
                        + sctfcoil_module.vol_case_cp * fwbs_variables.denstl
                    )

            # Cryo-aluminium conductor weights
            # Casing made of re-inforced aluminium alloy
            elif tfcoil_variables.i_tf_sup == 2:
                # Casing weight (CP only if physics_variables.itart = 1)bper leg/coil
                tfcoil_variables.whtcas = (
                    constants.dalu * vol_case / tfcoil_variables.n_tf
                )
                tfcoil_variables.whtconcu = 0.0e0
                tfcoil_variables.whtconal = (
                    constants.dalu * vol_cond / tfcoil_variables.n_tf
                )

                # Outer legs/CP weights
                if physics_variables.itart == 1:
                    # Weight of all the TF legs
                    tfcoil_variables.whttflgs = tfcoil_variables.n_tf * (
                        constants.dalu * vol_cond_leg
                        + tfcoil_variables.dcondins * (vol_ins_leg + vol_gr_ins_leg)
                    )

                    # CP weight
                    tfcoil_variables.whtcp = (
                        constants.dalu * tfcoil_variables.vol_cond_cp
                        + tfcoil_variables.dcondins
                        * (sctfcoil_module.vol_ins_cp + sctfcoil_module.vol_gr_ins_cp)
                        + sctfcoil_module.vol_case_cp * fwbs_variables.denstl
                    )

            # Turn insulation mass [kg]
            tfcoil_variables.whtconin = (
                tfcoil_variables.dcondins * vol_ins / tfcoil_variables.n_tf
            )

            # Ground wall insulation layer weight
            tfcoil_variables.whtgw = (
                tfcoil_variables.dcondins * vol_gr_ins / tfcoil_variables.n_tf
            )

            # Total weight
            tfcoil_variables.whttf = (
                tfcoil_variables.whtcas
                + tfcoil_variables.whtconcu
                + tfcoil_variables.whtconal
                + tfcoil_variables.whtconin
                + tfcoil_variables.whtgw
            ) * tfcoil_variables.n_tf

    def peak_tf_with_ripple(self, n_tf, wwp1, dr_tf_wp, tfin, bmaxtf):
        """Peak toroidal field on the conductor
        author: P J Knight, CCFE, Culham Science Centre
        This subroutine calculates the peak toroidal field at the
        outboard edge of the inboard TF coil winding pack, including
        the effects of ripple.
        <P>For 16, 18 or 20 coils, the calculation uses fitting formulae
        derived by M. Kovari using MAGINT calculations on coil sets based
        on a DEMO1 case.
        <P>For other numbers of coils, the original estimate using a 9%
        increase due to ripple from the axisymmetric calculation is used.
        M. Kovari, Toroidal Field Coils - Maximum Field and Ripple -
        Parametric Calculation, July 2014

        :param n_tf: number of TF coils
        :type n_tf: float
        :param wwp1: width of plasma-facing face of winding pack (m)
        :type wwp1: float
        :param dr_tf_wp: radial thickness of winding pack (m)
        :type dr_tf_wp: float
        :param tfin: major radius of centre of winding pack (m)
        :type tfin: float
        :param bmaxtf: nominal (axisymmetric) peak toroidal field (T)
        :type bmaxtf: float

        :returns: (bmaxtfrp, flag)
        * bmaxtfrp: peak toroidal field including ripple (T)
        * flag: flag warning of applicability problems

        :rtype: Tuple[float, int]

        """
        a = numpy.zeros((4,))
        flag = 0

        #  Set fitting coefficients for different numbers of TF coils

        int_n_tf = numpy.round(n_tf)

        if int_n_tf == 16:
            a[0] = 0.28101e0
            a[1] = 1.8481e0
            a[2] = -0.88159e0
            a[3] = 0.93834e0
        elif int_n_tf == 18:
            a[0] = 0.29153e0
            a[1] = 1.81600e0
            a[2] = -0.84178e0
            a[3] = 0.90426e0
        elif int_n_tf == 20:
            a[0] = 0.29853e0
            a[1] = 1.82130e0
            a[2] = -0.85031e0
            a[3] = 0.89808e0

        else:
            bmaxtfrp = 1.09e0 * bmaxtf
            return bmaxtfrp, flag

        #  Maximum winding pack width before adjacent packs touch
        #  (ignoring the external case and ground wall thicknesses)

        wmax = (2.0e0 * tfin + dr_tf_wp) * numpy.tan(numpy.pi / n_tf)

        #  Dimensionless winding pack width

        sctfcoil_module.tf_fit_t = wwp1 / wmax
        if (sctfcoil_module.tf_fit_t < 0.3e0) or (sctfcoil_module.tf_fit_t > 1.1e0):
            # write(*,*) 'PEAK_TF_WITH_RIPPLE: fitting problem; t = ',t
            flag = 1

        #  Dimensionless winding pack radial thickness

        sctfcoil_module.tf_fit_z = dr_tf_wp / wmax
        if (sctfcoil_module.tf_fit_z < 0.26e0) or (sctfcoil_module.tf_fit_z > 0.7e0):
            # write(*,*) 'PEAK_TF_WITH_RIPPLE: fitting problem; z = ',z
            flag = 2

        #  Ratio of peak field with ripple to nominal axisymmetric peak field

        sctfcoil_module.tf_fit_y = (
            a[0]
            + a[1] * numpy.exp(-sctfcoil_module.tf_fit_t)
            + a[2] * sctfcoil_module.tf_fit_z
            + a[3] * sctfcoil_module.tf_fit_z * sctfcoil_module.tf_fit_t
        )

        bmaxtfrp = sctfcoil_module.tf_fit_y * bmaxtf

        return bmaxtfrp, flag

    def res_tf_internal_geom(self):
        """
        Author : S. Kahn
        Resisitve TF turn geometry, equivalent to winding_pack subroutines

        """
        sctfcoil_module.r_wp_inner = (
            build_variables.r_tf_inboard_in + tfcoil_variables.thkcas
        )
        sctfcoil_module.r_wp_outer = (
            build_variables.r_tf_inboard_out - tfcoil_variables.casthi
        )

        # Conductor layer radial thickness at centercollumn top [m]
        if physics_variables.itart == 1:
            sctfcoil_module.dr_tf_wp_top = (
                build_variables.r_cp_top
                - tfcoil_variables.casthi
                - tfcoil_variables.thkcas
                - build_variables.r_tf_inboard_in
            )

        # Number of turns
        # Set by user (no turn structure by default, i.e. tfcoil_variables.n_tf_turn = 1 )
        if (
            abs(tfcoil_variables.n_tf_turn)
            < numpy.finfo(float(tfcoil_variables.n_tf_turn)).eps
        ):
            tfcoil_variables.n_tf_turn = 1.0e0

        # Total mid-plane cross-sectional area of winding pack, [m2]
        # including the surrounding ground-wall insulation layer
        sctfcoil_module.awpc = (
            numpy.pi
            * (sctfcoil_module.r_wp_outer**2 - sctfcoil_module.r_wp_inner**2)
            / tfcoil_variables.n_tf
        )

        # Area of the front case, the plasma-facing case of the inner TF coil [m2]
        sctfcoil_module.a_case_front = (
            numpy.pi
            * (
                (sctfcoil_module.r_wp_outer + tfcoil_variables.casthi) ** 2
                - sctfcoil_module.r_wp_outer**2
            )
            / tfcoil_variables.n_tf
        )

        # WP mid-plane cross-section excluding ground insulation per coil [m2]
        sctfcoil_module.awptf = numpy.pi * (
            (sctfcoil_module.r_wp_outer - tfcoil_variables.tinstf) ** 2
            - (sctfcoil_module.r_wp_inner + tfcoil_variables.tinstf) ** 2
        ) / tfcoil_variables.n_tf - 2.0e0 * tfcoil_variables.tinstf * (
            tfcoil_variables.dr_tf_wp - 2.0e0 * tfcoil_variables.tinstf
        )

        # Ground insulation cross-section area per coil [m2]
        sctfcoil_module.a_ground_ins = sctfcoil_module.awpc - sctfcoil_module.awptf

        # Exact mid-plane cross-section area of the conductor per TF coil [m2]
        a_tf_cond = numpy.pi * (
            (
                sctfcoil_module.r_wp_outer
                - tfcoil_variables.tinstf
                - tfcoil_variables.thicndut
            )
            ** 2
            - (
                sctfcoil_module.r_wp_inner
                + tfcoil_variables.tinstf
                + tfcoil_variables.thicndut
            )
            ** 2
        ) / tfcoil_variables.n_tf - (
            tfcoil_variables.dr_tf_wp
            - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.thicndut)
        ) * 2.0e0 * (
            tfcoil_variables.tinstf
            + tfcoil_variables.thicndut * tfcoil_variables.n_tf_turn
        )
        a_tf_cond = a_tf_cond * (1.0e0 - tfcoil_variables.fcoolcp)

        # Inter turn insulation area per coil [m2]
        tfcoil_variables.aiwp = sctfcoil_module.awptf - a_tf_cond / (
            1.0e0 - tfcoil_variables.fcoolcp
        )

        # Total insulation cross-section per coil [m2]
        sctfcoil_module.a_tf_ins = tfcoil_variables.aiwp + sctfcoil_module.a_ground_ins

        # Insulation fraction [-]
        sctfcoil_module.f_tf_ins = (
            tfcoil_variables.n_tf * sctfcoil_module.a_tf_ins / tfcoil_variables.tfareain
        )

        # Total cross-sectional area of the bucking cylindre and the outer support
        # support structure per coil [m2]
        # physics_variables.itart = 1 : Only valid at mid-plane
        tfcoil_variables.acasetf = (
            tfcoil_variables.tfareain / tfcoil_variables.n_tf
        ) - sctfcoil_module.awpc

        # Current per turn
        tfcoil_variables.cpttf = tfcoil_variables.ritfc / (
            tfcoil_variables.n_tf_turn * tfcoil_variables.n_tf
        )

        # Exact current density on TF oubard legs
        tfcoil_variables.cdtfleg = tfcoil_variables.ritfc / (
            (1.0e0 - tfcoil_variables.fcoolcp)
            * (
                tfcoil_variables.tftort
                - 2.0e0
                * (
                    tfcoil_variables.n_tf_turn * tfcoil_variables.thicndut
                    + tfcoil_variables.tinstf
                )
            )
            * (
                build_variables.tfthko
                - 2.0e0 * (tfcoil_variables.thicndut + tfcoil_variables.tinstf)
            )
        )

        # Reporting negative WP areas issues
        if sctfcoil_module.awpc < 0.0e0:
            error_handling.fdiags[0] = sctfcoil_module.awpc
            error_handling.fdiags[0] = tfcoil_variables.dr_tf_wp
            error_handling.report_error(99)

        elif sctfcoil_module.awptf < 0.0e0:
            error_handling.fdiags[0] = sctfcoil_module.awptf
            error_handling.report_error(101)

    def sc_tf_internal_geom(self, i_tf_wp_geom, i_tf_case_geom, i_tf_turns_integer):
        """
        Author : S. Kahn, CCFE
        Seting the WP, case and turns geometry for SC magnets
        """

        # Calculating the WP / ground insulation areas
        self.tf_wp_geom(i_tf_wp_geom)

        # Calculating the TF steel casing areas
        self.tf_case_geom(i_tf_wp_geom, i_tf_case_geom)

        # WP/trun currents
        self.tf_wp_currents()

        # sctfcoil_module.sc_tf_internal_geom(
        #     i_tf_wp_geom, i_tf_case_geom, i_tf_turns_integer
        # )

        # Setting the WP turn geometry / areas
        if i_tf_turns_integer == 0:
            # Non-ingeger number of turns
            (
                tfcoil_variables.acstf,
                tfcoil_variables.acndttf,
                tfcoil_variables.insulation_area,
                tfcoil_variables.n_tf_turn,
            ) = self.tf_averaged_turn_geom(
                tfcoil_variables.jwptf,
                tfcoil_variables.thwcndut,
                tfcoil_variables.thicndut,
                tfcoil_variables.i_tf_sc_mat,
            )

        else:
            # Integer number of turns
            (
                tfcoil_variables.acstf,
                tfcoil_variables.acndttf,
                tfcoil_variables.insulation_area,
                tfcoil_variables.cpttf,
                tfcoil_variables.n_tf_turn,
            ) = self.tf_integer_turn_geom(
                tfcoil_variables.n_layer,
                tfcoil_variables.n_pancake,
                tfcoil_variables.thwcndut,
                tfcoil_variables.thicndut,
            )

        # Areas and fractions
        # -------------------
        # Central helium channel down the conductor core [m2]
        tfcoil_variables.awphec = (
            0.25e0
            * tfcoil_variables.n_tf_turn
            * numpy.pi
            * tfcoil_variables.dhecoil**2
        )

        # Total conductor cross-sectional area, taking account of void area
        # and central helium channel [m2]
        tfcoil_variables.acond = (
            tfcoil_variables.acstf
            * tfcoil_variables.n_tf_turn
            * (1.0e0 - tfcoil_variables.vftf)
            - tfcoil_variables.awphec
        )

        # Void area in conductor for He, not including central channel [m2]
        tfcoil_variables.avwp = (
            tfcoil_variables.acstf * tfcoil_variables.n_tf_turn * tfcoil_variables.vftf
        )

        # Area of inter-turn insulation: total [m2]
        tfcoil_variables.aiwp = (
            tfcoil_variables.n_tf_turn * tfcoil_variables.insulation_area
        )

        # Area of steel structure in winding pack [m2]
        tfcoil_variables.aswp = tfcoil_variables.n_tf_turn * tfcoil_variables.acndttf

        # Inboard coil steel area [m2]
        sctfcoil_module.a_tf_steel = tfcoil_variables.acasetf + tfcoil_variables.aswp

        # Inboard coil steel fraction [-]
        sctfcoil_module.f_tf_steel = (
            tfcoil_variables.n_tf
            * sctfcoil_module.a_tf_steel
            / tfcoil_variables.tfareain
        )

        # Inboard coil insulation cross-section [m2]
        sctfcoil_module.a_tf_ins = tfcoil_variables.aiwp + sctfcoil_module.a_ground_ins

        #  Inboard coil insulation fraction [-]
        sctfcoil_module.f_tf_ins = (
            tfcoil_variables.n_tf * sctfcoil_module.a_tf_ins / tfcoil_variables.tfareain
        )

        # Negative areas or fractions error reporting
        if (
            tfcoil_variables.acond <= 0.0e0
            or tfcoil_variables.avwp <= 0.0e0
            or tfcoil_variables.aiwp <= 0.0e0
            or tfcoil_variables.aswp <= 0.0e0
            or sctfcoil_module.a_tf_steel <= 0.0e0
            or sctfcoil_module.f_tf_steel <= 0.0e0
            or sctfcoil_module.a_tf_ins <= 0.0e0
            or sctfcoil_module.f_tf_ins <= 0.0e0
        ):
            error_handling.fdiags[0] = tfcoil_variables.acond
            error_handling.fdiags[1] = tfcoil_variables.avwp
            error_handling.fdiags[2] = tfcoil_variables.aiwp
            error_handling.fdiags[3] = tfcoil_variables.aswp
            error_handling.fdiags[4] = sctfcoil_module.a_tf_steel
            error_handling.fdiags[5] = sctfcoil_module.f_tf_steel
            error_handling.fdiags[6] = sctfcoil_module.a_tf_ins
            error_handling.fdiags[7] = sctfcoil_module.f_tf_ins
            error_handling.report_error(276)

    def tf_wp_geom(self, i_tf_wp_geom):
        """
         Author : S. Kahn, CCFE
        Seting the WP geometry and area for SC magnets
        """
        sctfcoil_module.r_wp_inner = (
            build_variables.r_tf_inboard_in + tfcoil_variables.thkcas
        )

        # Radial position of outer edge of winding pack [m]
        sctfcoil_module.r_wp_outer = (
            sctfcoil_module.r_wp_inner + tfcoil_variables.dr_tf_wp
        )

        # Radius of geometrical centre of winding pack [m]
        sctfcoil_module.r_wp_centre = 0.5e0 * (
            sctfcoil_module.r_wp_inner + sctfcoil_module.r_wp_outer
        )

        # TF toroidal thickness at the WP inner radius [m]
        t_tf_at_wp = 2.0e0 * sctfcoil_module.r_wp_inner * sctfcoil_module.tan_theta_coil

        # Minimal toroidal thickness of winding pack [m]
        sctfcoil_module.t_wp_toroidal = t_tf_at_wp - 2.0e0 * tfcoil_variables.casths

        # Rectangular WP
        # --------------
        if i_tf_wp_geom == 0:
            # Outer WP layer toroidal thickness [m]
            tfcoil_variables.wwp1 = sctfcoil_module.t_wp_toroidal

            # Averaged toroidal thickness of of winding pack [m]
            sctfcoil_module.t_wp_toroidal_av = sctfcoil_module.t_wp_toroidal

            # Total cross-sectional area of winding pack [m2]
            sctfcoil_module.awpc = (
                tfcoil_variables.dr_tf_wp * sctfcoil_module.t_wp_toroidal
            )

            # WP cross-section without insertion gap and ground insulation [m2]
            sctfcoil_module.awptf = (
                tfcoil_variables.dr_tf_wp
                - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
            ) * (
                sctfcoil_module.t_wp_toroidal
                - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
            )

            # Cross-section area of the WP ground insulation [m2]
            sctfcoil_module.a_ground_ins = (
                tfcoil_variables.dr_tf_wp - 2.0e0 * tfcoil_variables.tfinsgap
            ) * (
                sctfcoil_module.t_wp_toroidal - 2.0e0 * tfcoil_variables.tfinsgap
            ) - sctfcoil_module.awptf

        # Double rectangular WP
        # ---------------------
        elif i_tf_wp_geom == 1:
            # Thickness of winding pack section at R > sctfcoil_module.r_wp_centre [m]
            tfcoil_variables.wwp1 = 2.0e0 * (
                sctfcoil_module.r_wp_centre * sctfcoil_module.tan_theta_coil
                - tfcoil_variables.casths
            )

            # Thickness of winding pack section at R < sctfcoil_module.r_wp_centre [m]
            tfcoil_variables.wwp2 = 2.0e0 * (
                sctfcoil_module.r_wp_inner * sctfcoil_module.tan_theta_coil
                - tfcoil_variables.casths
            )

            # Averaged toroidal thickness of of winding pack [m]
            sctfcoil_module.t_wp_toroidal_av = 0.5e0 * (
                tfcoil_variables.wwp1 + tfcoil_variables.wwp2
            )

            # Total cross-sectional area of winding pack [m2]
            # Including ground insulation and insertion gap
            sctfcoil_module.awpc = (
                tfcoil_variables.dr_tf_wp * sctfcoil_module.t_wp_toroidal_av
            )

            # WP cross-section without insertion gap and ground insulation [m2]
            sctfcoil_module.awptf = (
                0.5e0
                * (
                    tfcoil_variables.dr_tf_wp
                    - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
                )
                * (
                    tfcoil_variables.wwp1
                    + tfcoil_variables.wwp2
                    - 4.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
                )
            )

            # Cross-section area of the WP ground insulation [m2]
            sctfcoil_module.a_ground_ins = (
                0.5e0
                * (tfcoil_variables.dr_tf_wp - 2.0e0 * tfcoil_variables.tfinsgap)
                * (
                    tfcoil_variables.wwp1
                    + tfcoil_variables.wwp2
                    - 4.0e0 * tfcoil_variables.tfinsgap
                )
                - sctfcoil_module.awptf
            )

        # Trapezoidal WP
        # --------------
        else:
            # Thickness of winding pack section at sctfcoil_module.r_wp_outer [m]
            tfcoil_variables.wwp1 = 2.0e0 * (
                sctfcoil_module.r_wp_outer * sctfcoil_module.tan_theta_coil
                - tfcoil_variables.casths
            )

            # Thickness of winding pack section at sctfcoil_module.r_wp_inner [m]
            tfcoil_variables.wwp2 = 2.0e0 * (
                sctfcoil_module.r_wp_inner * sctfcoil_module.tan_theta_coil
                - tfcoil_variables.casths
            )

            # Averaged toroidal thickness of of winding pack [m]
            sctfcoil_module.t_wp_toroidal_av = 0.5e0 * (
                tfcoil_variables.wwp1 + tfcoil_variables.wwp2
            )

            # Total cross-sectional area of winding pack [m2]
            # Including ground insulation and insertion gap
            sctfcoil_module.awpc = tfcoil_variables.dr_tf_wp * (
                tfcoil_variables.wwp2
                + 0.5e0 * (tfcoil_variables.wwp1 - tfcoil_variables.wwp2)
            )

            # WP cross-section without insertion gap and ground insulation [m2]
            sctfcoil_module.awptf = (
                tfcoil_variables.dr_tf_wp
                - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
            ) * (
                tfcoil_variables.wwp2
                - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
                + 0.5e0 * (tfcoil_variables.wwp1 - tfcoil_variables.wwp2)
            )

            # Cross-section area of the WP ground insulation [m2]
            sctfcoil_module.a_ground_ins = (
                tfcoil_variables.dr_tf_wp - 2.0e0 * tfcoil_variables.tfinsgap
            ) * (
                tfcoil_variables.wwp2
                - 2.0e0 * tfcoil_variables.tfinsgap
                + 0.5e0 * (tfcoil_variables.wwp1 - tfcoil_variables.wwp2)
            ) - sctfcoil_module.awptf

        # --------------

        # Negative WP area error reporting
        if sctfcoil_module.awptf <= 0.0e0 or sctfcoil_module.awpc <= 0.0e0:
            error_handling.fdiags[0] = sctfcoil_module.awptf
            error_handling.fdiags[1] = sctfcoil_module.awpc
            error_handling.report_error(99)

    def tf_case_geom(self, i_tf_wp_geom, i_tf_case_geom):
        """
        Author : S. Kahn, CCFE
        Setting the case geometry and area for SC magnets
        """
        tfcoil_variables.acasetf = (
            tfcoil_variables.tfareain / tfcoil_variables.n_tf
        ) - sctfcoil_module.awpc

        # Outboard leg cross-sectional area of surrounding case [m2]
        tfcoil_variables.acasetfo = tfcoil_variables.arealeg - sctfcoil_module.awpc

        # Front casing area [m2]
        if i_tf_case_geom == 0:
            # Circular front case
            sctfcoil_module.a_case_front = (
                sctfcoil_module.theta_coil * build_variables.r_tf_inboard_out**2
                - sctfcoil_module.tan_theta_coil * sctfcoil_module.r_wp_outer**2
            )
        else:
            # Straight front case
            sctfcoil_module.a_case_front = (
                (sctfcoil_module.r_wp_outer + tfcoil_variables.casthi) ** 2
                - sctfcoil_module.r_wp_outer**2
            ) * sctfcoil_module.tan_theta_coil

        # Nose casing area [m2]
        sctfcoil_module.a_case_nose = (
            sctfcoil_module.tan_theta_coil * sctfcoil_module.r_wp_inner**2
            - sctfcoil_module.theta_coil * build_variables.r_tf_inboard_in**2
        )

        # Report error if the casing area is negative
        if tfcoil_variables.acasetf <= 0.0e0 or tfcoil_variables.acasetfo <= 0.0e0:
            error_handling.fdiags[0] = tfcoil_variables.acasetf
            error_handling.fdiags[1] = tfcoil_variables.acasetfo
            error_handling.report_error(99)

        # Average lateral casing thickness
        # --------------
        # Rectangular casing
        if i_tf_wp_geom == 0:
            sctfcoil_module.t_lat_case_av = (
                tfcoil_variables.casths
                + 0.5e0 * sctfcoil_module.tan_theta_coil * tfcoil_variables.dr_tf_wp
            )

        # Double rectangular WP
        elif i_tf_wp_geom == 1:
            sctfcoil_module.t_lat_case_av = (
                tfcoil_variables.casths
                + 0.25e0 * sctfcoil_module.tan_theta_coil * tfcoil_variables.dr_tf_wp
            )

        # Trapezoidal WP
        else:
            sctfcoil_module.t_lat_case_av = tfcoil_variables.casths

        # --------------

    def tf_integer_turn_geom(self, n_layer, n_pancake, thwcndut, thicndut):
        """
        Authors: J Morris & S Khan
        Setting the TF WP turn geometry for SC magnets from the number
        of turns rows in the radial direction. The turns can have any
        rectangular shapes.
        This calculation has two purposes, first to check if a turn can exist
        (positive cable space) and the second to provide its dimenesions,
        areas and the its associated current

        """
        sctfcoil_module.rbcndut = thwcndut * 0.75e0

        # Radial turn dimension [m]
        sctfcoil_module.t_turn_radial = (
            tfcoil_variables.dr_tf_wp
            - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
        ) / n_layer

        if sctfcoil_module.t_turn_radial <= (2.0e0 * thicndut + 2.0e0 * thwcndut):
            error_handling.fdiags[0] = sctfcoil_module.t_turn_radial
            error_handling.fdiags[1] = thicndut
            error_handling.fdiags[2] = thwcndut
            error_handling.report_error(100)

        # Toroidal turn dimension [m]
        sctfcoil_module.t_turn_toroidal = (
            sctfcoil_module.t_wp_toroidal
            - 2.0e0 * (tfcoil_variables.tinstf + tfcoil_variables.tfinsgap)
        ) / n_pancake

        if sctfcoil_module.t_turn_toroidal <= (2.0e0 * thicndut + 2.0e0 * thwcndut):
            error_handling.fdiags[0] = sctfcoil_module.t_turn_toroidal
            error_handling.fdiags[1] = thicndut
            error_handling.fdiags[2] = thwcndut
            error_handling.report_error(100)

        tfcoil_variables.t_turn_tf = numpy.sqrt(
            sctfcoil_module.t_turn_radial * sctfcoil_module.t_turn_toroidal
        )

        # Number of TF turns
        n_tf_turn = numpy.double(n_layer * n_pancake)

        # Current per turn [A/turn]
        cpttf = sctfcoil_module.tfc_current / n_tf_turn

        # Radial and toroidal dimension of conductor [m]
        sctfcoil_module.t_conductor_radial = (
            sctfcoil_module.t_turn_radial - 2.0e0 * thicndut
        )
        sctfcoil_module.t_conductor_toroidal = (
            sctfcoil_module.t_turn_toroidal - 2.0e0 * thicndut
        )
        tfcoil_variables.t_conductor = numpy.sqrt(
            sctfcoil_module.t_conductor_radial * sctfcoil_module.t_conductor_toroidal
        )

        # Dimension of square cable space inside conduit [m]
        sctfcoil_module.t_cable_radial = (
            sctfcoil_module.t_conductor_radial - 2.0e0 * thwcndut
        )
        sctfcoil_module.t_cable_toroidal = (
            sctfcoil_module.t_conductor_toroidal - 2.0e0 * thwcndut
        )
        sctfcoil_module.t_cable = numpy.sqrt(
            sctfcoil_module.t_cable_radial * sctfcoil_module.t_cable_toroidal
        )

        # Cross-sectional area of cable space per turn
        # taking account of rounded inside corners [m2]
        acstf = (sctfcoil_module.t_cable_radial * sctfcoil_module.t_cable_toroidal) - (
            4.0e0 - numpy.pi
        ) * sctfcoil_module.rbcndut**2

        if acstf <= 0.0e0:
            if (sctfcoil_module.t_cable_radial < 0.0e0) or (
                sctfcoil_module.t_cable_toroidal < 0.0e0
            ):
                error_handling.fdiags[0] = acstf
                error_handling.fdiags[1] = sctfcoil_module.t_cable_radial
                error_handling.fdiags[2] = sctfcoil_module.t_cable_toroidal
                error_handling.report_error(101)
            else:
                error_handling.fdiags[0] = acstf
                error_handling.fdiags[1] = sctfcoil_module.t_cable_radial
                error_handling.fdiags[1] = sctfcoil_module.t_cable_toroidal
                error_handling.report_error(102)
                sctfcoil_module.rbcndut = 0.0e0
                acstf = (
                    sctfcoil_module.t_cable_radial * sctfcoil_module.t_cable_toroidal
                )

        # Cross-sectional area of conduit jacket per turn [m2]
        acndttf = (
            sctfcoil_module.t_conductor_radial * sctfcoil_module.t_conductor_toroidal
            - acstf
        )

        # Area of inter-turn insulation: single turn [m2]
        insulation_area = (
            sctfcoil_module.t_turn_radial * sctfcoil_module.t_turn_toroidal
            - acndttf
            - acstf
        )
        return acstf, acndttf, insulation_area, cpttf, n_tf_turn

        # -------------

    def tf_wp_currents(self):
        """
        Author : S. Kahn, CCFE
        Turn engineering turn currents/densities
        """
        tfcoil_variables.jwptf = max(
            1.0e0,
            tfcoil_variables.ritfc / (tfcoil_variables.n_tf * sctfcoil_module.awptf),
        )

    @staticmethod
    @numba.njit(cache=True)
    def stresscl(
        n_tf_layer,
        n_radial_array,
        n_tf_wp_layers,
        i_tf_bucking,
        r_tf_inboard_in,
        bore,
        hmax,
        ohhghf,
        ohcth,
        tf_in_cs,
        tfcth,
        gapoh,
        ipfres,
        coheof,
        cohbop,
        cptdin,
        ncls,
        ld_ratio_cst,
        r_out_cst,
        oh_steel_frac,
        eyoung_steel,
        poisson_steel,
        eyoung_cond_axial,
        poisson_cond_axial,
        eyoung_cond_trans,
        poisson_cond_trans,
        eyoung_ins,
        poisson_ins,
        thicndut,
        eyoung_copper,
        poisson_copper,
        i_tf_sup,
        eyoung_res_tf_buck,
        r_wp_inner,
        tan_theta_coil,
        theta_coil,
        r_wp_outer,
        a_tf_steel,
        a_case_front,
        a_case_nose,
        tfinsgap,
        tinstf,
        n_tf_turn,
        i_tf_turns_integer,
        t_cable,
        t_cable_radial,
        dhecoil,
        fcutfsu,
        thwcndut,
        t_lat_case_av,
        t_wp_toroidal_av,
        a_tf_ins,
        aswp,
        acond,
        awpc,
        eyoung_al,
        poisson_al,
        fcoolcp,
        n_tf_graded_layers,
        ritfc,
        casthi,
        i_tf_stress_model,
        vforce_inboard_tot,
        i_tf_tresca,
        acasetf,
        vforce,
        acndttf,
    ):
        """TF coil stress routine


        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        author: S Kahn, CCFE, Culham Science Centre
        author: J Galambos, FEDC/ORNL

        This subroutine sets up the stress calculations for the
        TF coil set.
        PROCESS Superconducting TF Coil Model, J. Morris, CCFE, 1st May 2014
        """
        jeff = numpy.zeros((n_tf_layer,))
        # Effective current density [A/m2]

        radtf = numpy.zeros((n_tf_layer + 1,))
        # Radii used to define the layers used in the stress models [m]
        # Layers are labelled from inboard to outbard

        eyoung_trans = numpy.zeros((n_tf_layer,))
        # Young's moduli (one per layer) of the TF coil in the
        # transverse (radial/toroidal) direction. Used in the stress
        # models [Pa]

        poisson_trans = numpy.zeros(
            (n_tf_layer,),
        )
        # Poisson's ratios (one per layer) of the TF coil between the
        # two transverse directions (radial and toroidal). Used in the
        # stress models.

        eyoung_member_array = numpy.zeros((n_tf_wp_layers,))
        # Array to store the Young's moduli of the members to composite into smeared
        # properties [Pa]

        poisson_member_array = numpy.zeros((n_tf_wp_layers,))
        # Array to store the Poisson's ratios of the members to composite into smeared
        # properti

        l_member_array = numpy.zeros((n_tf_wp_layers,))
        # Array to store the linear dimension (thickness) of the members to composite into smeared
        # properties [m]

        eyoung_axial = numpy.zeros((n_tf_layer,))
        # Young's moduli (one per layer) of the TF coil in the vertical
        # direction used in the stress models [Pa]

        poisson_axial = numpy.zeros((n_tf_layer,))
        # Poisson's ratios (one per layer) of the TF coil between the
        # vertical and transverse directions (in that order). Used in the
        # stress models. d(transverse strain)/d(vertical strain) with
        # only vertical stress.

        sig_tf_wp_av_z = numpy.zeros(((n_tf_layer - i_tf_bucking) * n_radial_array,))
        # TF Inboard leg WP smeared vertical stress r distribution at mid-plane [Pa]

        sig_tf_r_max = numpy.zeros((n_tf_layer,))
        # Radial stress of the point of maximum shear stress of each layer [Pa]

        sig_tf_t_max = numpy.zeros((n_tf_layer,))
        # Toroidal stress of the point of maximum shear stress of each layer [Pa]

        sig_tf_z_max = numpy.zeros((n_tf_layer,))
        # Vertical stress of the point of maximum shear stress of each layer [Pa]
        # Rem : Currently constant but will be r dependent in the future

        sig_tf_vmises_max = numpy.zeros((n_tf_layer,))
        # Von-Mises stress of the point of maximum shear stress of each layer [Pa]

        sig_tf_tresca_max = numpy.zeros((n_tf_layer,))
        # Maximum shear stress, for the Tresca yield criterion of each layer [Pa]
        # If the CEA correction is addopted, the CEA corrected value is used

        sig_tf_z = numpy.zeros((n_tf_layer * n_radial_array,))
        # TF Inboard leg vertical tensile stress at mid-plane [Pa]

        sig_tf_smeared_r = numpy.zeros((n_tf_layer * n_radial_array,))
        # TF Inboard leg radial smeared stress r distribution at mid-plane [Pa]

        sig_tf_smeared_t = numpy.zeros((n_tf_layer * n_radial_array,))
        # TF Inboard leg tangential smeared stress r distribution at mid-plane [Pa]

        sig_tf_smeared_z = numpy.zeros((n_tf_layer * n_radial_array,))
        # TF Inboard leg vertical smeared stress r distribution at mid-plane [Pa]

        fac_sig_z_wp_av = 0.0
        # WP averaged vertical stress unsmearing factor

        str_tf_r = numpy.zeros((n_radial_array * n_tf_layer,))
        str_tf_t = numpy.zeros((n_radial_array * n_tf_layer,))
        str_tf_z = numpy.zeros((n_radial_array * n_tf_layer,))

        sig_tf_case = None
        sig_tf_cs_bucked = None
        str_wp = None
        casestr = None
        insstrain = None

        if abs(r_tf_inboard_in) < numpy.finfo(float(r_tf_inboard_in)).eps:
            # New extended plane strain model can handle it
            if i_tf_stress_model != 2:
                raise ValueError("r_tf_inboard_in is ~= 0", 245)

        # TODO: following is no longer used/needed?
        # if tfcoil_variables.acstf >= 0.0e0:
        #     tcbs = numpy.sqrt(tfcoil_variables.acstf)
        # else:
        #     tcbs = 0.0e0

        # LAYER ELASTIC PROPERTIES
        # ------------------------
        # Number of bucking layers
        n_tf_bucking = i_tf_bucking

        # CS properties (bucked and wedged)
        # ---
        if i_tf_bucking >= 2:
            # Calculation performed at CS flux swing (no current on CS)
            jeff[0] = 0.0e0

            # Inner radius of the CS
            if tf_in_cs == 1:
                # CS not used as wedge support tf_in_cs = 1
                radtf[0] = 0.001
            else:
                radtf[0] = bore

            # Superconducting CS
            if ipfres == 0:
                # Getting the turn dimention from scratch
                # as the TF is called before CS in caller.f90
                # -#

                # CS vertical cross-section area [m2]
                if tf_in_cs == 1:
                    a_oh = 2.0e0 * hmax * ohhghf * (bore - tfcth)
                else:
                    a_oh = 2.0e0 * hmax * ohhghf * ohcth

                # Maximum current in Central Solenoid, at either BOP or EOF [MA-turns]
                # Absolute value
                curr_oh_max = 1.0e-6 * numpy.maximum(coheof, cohbop) * a_oh

                #  Number of turns
                n_oh_turns = 1.0e6 * curr_oh_max / cptdin[sum(ncls)]

                # CS Turn vertical cross-sectionnal area
                a_oh_turn = a_oh / n_oh_turns

                # CS coil turn geometry calculation - stadium shape
                # Literature: https://doi.org/10.1016/j.fusengdes.2017.04.052
                d_cond_cst = (
                    a_oh_turn / ld_ratio_cst
                ) ** 0.5  # width of cs turn conduit
                l_cond_cst = ld_ratio_cst * d_cond_cst  # length of cs turn conduit
                # Radius of turn space = r_in_cst
                # Radius of curved outer corrner r_out_cst = 3mm from literature
                # ld_ratio_cst = 70 / 22 from literature
                p1 = ((l_cond_cst - d_cond_cst) / numpy.pi) ** 2
                p2 = (
                    (l_cond_cst * d_cond_cst)
                    - (4 - numpy.pi) * (r_out_cst**2)
                    - (a_oh_turn * oh_steel_frac)
                ) / numpy.pi
                r_in_cst = -((l_cond_cst - d_cond_cst) / numpy.pi) + numpy.sqrt(p1 + p2)
                t_cond_oh = (
                    d_cond_cst / 2
                ) - r_in_cst  # thickness of steel conduit in cs turn

                # OH/CS conduit thickness calculated assuming square conduit [m]
                # The CS insulation layer is assumed to the same as the TF one

                # CS turn cable space thickness
                t_cable_oh = r_in_cst * 2
                # -#

                # Smeared elastic properties of the CS
                # These smearing functions were written assuming transverse-
                # isotropic materials; that is not true of the CS, where the
                # stiffest dimension is toroidal and the radial and vertical
                # dimension are less stiff. Nevertheless this attempts to
                # hit the mark.
                # [EDIT: eyoung_cond is for the TF coil, not the CS coil]

                # Get transverse properties
                (eyoung_trans[0], a_working, poisson_trans[0],) = eyoung_parallel(
                    eyoung_steel,
                    oh_steel_frac,
                    poisson_steel,
                    eyoung_cond_axial,
                    1e0 - oh_steel_frac,
                    poisson_cond_axial,
                )

                # Get vertical properties
                # Split up into "members", concentric squares in cross section
                # (described in Figure 10 of the TF coil documentation)
                # Conductor
                eyoung_member_array[0] = eyoung_cond_trans
                poisson_member_array[0] = poisson_cond_trans
                l_member_array[0] = t_cable_oh
                # Steel conduit
                eyoung_member_array[1] = eyoung_steel
                poisson_member_array[1] = poisson_steel
                l_member_array[1] = 2 * t_cond_oh
                # Insulation
                eyoung_member_array[2] = eyoung_ins
                poisson_member_array[2] = poisson_ins
                l_member_array[2] = 2 * thicndut
                # [EDIT: Add central cooling channel? Would be new member #1]

                # Compute the composited (smeared) properties
                (
                    eyoung_axial[0],
                    a_working,
                    poisson_axial[0],
                    eyoung_cs_stiffest_leg,
                ) = eyoung_t_nested_squares(
                    3, eyoung_member_array, l_member_array, poisson_member_array
                )

            # resistive CS (copper)
            else:
                # Here is a rough approximation
                eyoung_trans[0] = eyoung_copper
                eyoung_axial[0] = eyoung_copper
                poisson_trans[0] = poisson_copper
                poisson_axial[0] = poisson_copper

        # ---

        # CS-TF interlayer properties
        # ---
        if i_tf_bucking == 3:
            # No current in this layer
            jeff[1] = 0.0e0

            # Outer radius of the CS
            if tf_in_cs == 1:
                radtf[1] = bore - tfcth - gapoh
            else:
                radtf[1] = bore + ohcth

            # Assumed to be Kapton for the moment
            # Ref : https://www.dupont.com/content/dam/dupont/products-and-services/membranes-and-films/polyimde-films/documents/DEC-Kapton-summary-of-properties.pdf
            eyoung_trans[1] = 2.5e9
            eyoung_axial[1] = 2.5e9
            poisson_trans[1] = 0.34e0  # Default value for young modulus
            poisson_axial[1] = 0.34e0  # Default value for young modulus

        # ---

        # bucking cylinder/casing properties
        # ---
        if i_tf_bucking >= 1:
            # No current in bucking cylinder/casing
            jeff[n_tf_bucking - 1] = 0.0e0

            if i_tf_sup == 1:
                eyoung_trans[n_tf_bucking - 1] = eyoung_steel
                eyoung_axial[n_tf_bucking - 1] = eyoung_steel
                poisson_trans[n_tf_bucking - 1] = poisson_steel
                poisson_axial[n_tf_bucking - 1] = poisson_steel

            # Bucking cylinder properties
            else:
                eyoung_trans[n_tf_bucking - 1] = eyoung_res_tf_buck
                eyoung_axial[n_tf_bucking - 1] = eyoung_res_tf_buck
                poisson_trans[n_tf_bucking - 1] = poisson_steel  # Seek better value #
                poisson_axial[n_tf_bucking - 1] = poisson_steel  # Seek better value #

            # Innernost TF casing radius
            radtf[n_tf_bucking - 1] = r_tf_inboard_in

        # ---

        # (Super)conductor layer properties
        # ---
        # SC coil
        if i_tf_sup == 1:
            # Inner/outer radii of the layer representing the WP in stress calculations [m]
            # These radii are chosen to preserve the true WP area; see Issue #1048
            r_wp_inner_eff = r_wp_inner * numpy.sqrt(tan_theta_coil / theta_coil)
            r_wp_outer_eff = r_wp_outer * numpy.sqrt(tan_theta_coil / theta_coil)

            # Area of the cylinder representing the WP in stress calculations [m2]
            a_wp_eff = (r_wp_outer_eff**2 - r_wp_inner_eff**2) * theta_coil

            # Steel cross-section under the area representing the WP in stress calculations [m2]
            a_wp_steel_eff = a_tf_steel - a_case_front - a_case_nose

            # WP effective insulation thickness (SC only) [m]
            # include groundwall insulation + insertion gap in tfcoil_variables.thicndut
            # inertion gap is tfcoil_variables.tfinsgap on 4 sides
            t_ins_eff = thicndut + (tfinsgap + tinstf) / n_tf_turn

            # Effective WP young modulus in the toroidal direction [Pa]
            # The toroidal property drives the stress calculation (J. Last report no 4)
            # Hence, the radial direction is relevant for the property smearing
            # Rem : This assumption might be re-defined for bucked and wedged design
            if i_tf_turns_integer == 0:
                # Non-integer number of turns
                t_cable_eyng = t_cable
            else:
                # Integer number of turns
                t_cable_eyng = t_cable_radial

            # Average WP Young's modulus in the transverse
            # (radial and toroidal) direction
            # Split up into "members", concentric squares in cross section
            # (described in Figure 10 of the TF coil documentation)
            # Helium
            eyoung_member_array[0] = 0e0
            poisson_member_array[0] = poisson_steel
            l_member_array[0] = dhecoil
            # Conductor and co-wound copper
            (
                eyoung_member_array[1],
                l_member_array[1],
                poisson_member_array[1],
            ) = eyoung_series(
                numpy.double(eyoung_cond_trans),
                (t_cable_eyng - dhecoil) * (1.0e0 - fcutfsu),
                numpy.double(poisson_cond_trans),
                numpy.double(eyoung_copper),
                (t_cable_eyng - dhecoil) * fcutfsu,
                numpy.double(poisson_copper),
            )
            # Steel conduit
            eyoung_member_array[2] = eyoung_steel
            poisson_member_array[2] = poisson_steel
            l_member_array[2] = 2 * thwcndut
            # Insulation
            eyoung_member_array[3] = eyoung_ins
            poisson_member_array[3] = poisson_ins
            l_member_array[3] = 2 * t_ins_eff

            # Compute the composited (smeared) properties
            (
                eyoung_wp_trans,
                a_working,
                poisson_wp_trans,
                eyoung_wp_stiffest_leg,
            ) = eyoung_t_nested_squares(
                4,
                eyoung_member_array,
                l_member_array,
                poisson_member_array,
            )

            # Lateral casing correction (series-composition)
            (eyoung_wp_trans_eff, a_working, poisson_wp_trans_eff,) = eyoung_series(
                eyoung_wp_trans,
                numpy.double(t_wp_toroidal_av),
                poisson_wp_trans,
                numpy.double(eyoung_steel),
                2.0e0 * t_lat_case_av,
                numpy.double(poisson_steel),
            )

            poisson_wp_trans_eff = numpy.double(poisson_wp_trans_eff)

            # Average WP Young's modulus in the vertical direction
            # Split up into "members", concentric squares in cross section
            # (described in Figure 10 of the TF coil documentation)
            # Steel conduit
            eyoung_member_array[0] = eyoung_steel
            poisson_member_array[0] = poisson_steel
            l_member_array[0] = aswp
            # Insulation
            eyoung_member_array[1] = eyoung_ins
            poisson_member_array[1] = poisson_ins
            l_member_array[1] = a_tf_ins
            # Copper
            eyoung_member_array[2] = eyoung_copper
            poisson_member_array[2] = poisson_copper
            l_member_array[2] = acond * fcutfsu
            # Conductor
            eyoung_member_array[3] = eyoung_cond_axial
            poisson_member_array[3] = poisson_cond_axial
            l_member_array[3] = acond * (1.0e0 - fcutfsu)
            # Helium and void
            eyoung_member_array[4] = 0e0
            poisson_member_array[4] = poisson_steel
            l_member_array[4] = awpc - acond - a_tf_ins - aswp
            # Compute the composite / smeared properties:
            (eyoung_wp_axial, a_working, poisson_wp_axial,) = eyoung_parallel_array(
                5,
                eyoung_member_array,
                l_member_array,
                poisson_member_array,
            )

            # Average WP Young's modulus in the vertical direction, now including the lateral case
            # Parallel-composite the steel and insulation, now including the lateral case (sidewalls)
            (eyoung_wp_axial_eff, a_working, poisson_wp_axial_eff,) = eyoung_parallel(
                eyoung_steel,
                a_wp_steel_eff - aswp,
                poisson_steel,
                eyoung_wp_axial,
                a_working,
                poisson_wp_axial,
            )

        # Resistive coil
        else:
            # Picking the conductor material Young's modulus
            if i_tf_sup == 0:
                eyoung_cond = eyoung_copper
                poisson_cond = poisson_copper
            elif i_tf_sup == 2:
                eyoung_cond = eyoung_al
                poisson_cond = poisson_al

            # Effective WP young modulus in the toroidal direction [Pa]
            # Rem : effect of cooling pipes and insulation not taken into account
            #       for now as it needs a radially dependent Young modulus
            eyoung_wp_trans_eff = eyoung_cond
            eyoung_wp_trans = numpy.double(eyoung_cond)
            poisson_wp_trans_eff = numpy.double(poisson_cond)
            poisson_wp_trans = numpy.double(poisson_cond)

            # WP area using the stress model circular geometry (per coil) [m2]
            a_wp_eff = (r_wp_outer**2 - r_wp_inner**2) * theta_coil

            # Effective conductor region young modulus in the vertical direction [Pa]
            # Parallel-composite conductor and insulator
            (eyoung_wp_axial, a_working, poisson_wp_axial,) = eyoung_parallel(
                eyoung_cond,
                (a_wp_eff - a_tf_ins) * (1.0e0 - fcoolcp),
                poisson_cond,
                eyoung_ins,
                a_tf_ins,
                poisson_ins,
            )
            # Parallel-composite cooling pipes into that
            (eyoung_wp_axial, a_working, poisson_wp_axial,) = eyoung_parallel(
                0e0,
                (a_wp_eff - a_tf_ins) * fcoolcp,
                poisson_cond,
                eyoung_wp_axial,
                a_working,
                poisson_wp_axial,
            )

            # Effective young modulus used in stress calculations
            eyoung_wp_axial_eff = eyoung_wp_axial
            poisson_wp_axial_eff = poisson_wp_axial

            # Effect conductor layer inner/outer radius
            r_wp_inner_eff = numpy.double(r_wp_inner)
            r_wp_outer_eff = numpy.double(r_wp_outer)

        # Thickness of the layer representing the WP in stress calcualtions [m]
        dr_tf_wp_eff = r_wp_outer_eff - r_wp_outer_eff

        # Thickness of WP with homogeneous stress property [m]
        dr_wp_layer = dr_tf_wp_eff / n_tf_graded_layers

        for ii in range(numpy.intc(n_tf_graded_layers)):
            # Homogeneous current in (super)conductor
            jeff[n_tf_bucking + ii] = ritfc / (
                numpy.pi * (r_wp_outer_eff**2 - r_wp_inner_eff**2)
            )

            # Same thickness for all WP layers in stress calculation
            radtf[n_tf_bucking + ii] = r_wp_inner_eff + ii * dr_wp_layer

            # Young modulus
            eyoung_trans[n_tf_bucking + ii] = eyoung_wp_trans_eff
            eyoung_axial[n_tf_bucking + ii] = eyoung_wp_axial_eff

            # Poisson's ratio
            poisson_trans[n_tf_bucking + ii] = poisson_wp_trans_eff
            poisson_axial[n_tf_bucking + ii] = poisson_wp_axial_eff

        # Steel case on the plasma side of the inboard TF coil
        # As per Issue #1509
        jeff[n_tf_layer - 1] = 0.0e0
        radtf[n_tf_layer - 1] = r_wp_outer_eff
        eyoung_trans[n_tf_layer - 1] = eyoung_steel
        eyoung_axial[n_tf_layer - 1] = eyoung_steel
        poisson_trans[n_tf_layer - 1] = poisson_steel
        poisson_axial[n_tf_layer - 1] = poisson_steel

        # last layer radius
        radtf[n_tf_layer] = r_wp_outer_eff + casthi

        # The ratio between the true cross sectional area of the
        # front case, and that considered by the plane strain solver
        f_tf_stress_front_case = (
            a_case_front
            / theta_coil
            / (radtf[n_tf_layer] ** 2 - radtf[n_tf_layer - 1] ** 2)
        )

        # Correct for the missing axial stiffness from the missing
        # outer case steel as per the updated description of
        # Issue #1509
        eyoung_axial[n_tf_layer - 1] = (
            eyoung_axial[n_tf_layer - 1] * f_tf_stress_front_case
        )

        # ---
        # ------------------------

        # RADIAL STRESS SUBROUTINES CALL
        # ------------------------------
        # Stress model not valid the TF does not contain any hole
        # (Except if i_tf_plane_stress == 2; see Issue 1414)
        # Current action : trigger and error and add a little hole
        #                  to allow stress calculations
        # Rem SK : Can be easily ameneded playing around the boundary conditions
        if abs(radtf[0]) < numpy.finfo(float(radtf[0])).eps:
            # New extended plane strain model can handle it
            if i_tf_stress_model != 2:
                # error_handling.report_error(245)
                radtf[0] = 1.0e-9
            # elif abs(radtf[1]) < numpy.finfo(float(radtf[0])).eps:
            #     logger.error(
            #         """ERROR: First TF layer has zero thickness.
            #     Perhaps you meant to have thkcas nonzero or tfcoil_variables.i_tf_bucking = 0?
            #     """
            #     )

        # ---

        # Old generalized plane stress model
        # ---
        if i_tf_stress_model == 1:
            # Plane stress calculation (SC) [Pa]

            (sig_tf_r, sig_tf_t, deflect, radial_array,) = plane_stress(
                nu=poisson_trans,
                rad=radtf,
                ey=eyoung_trans,
                j=jeff,
                nlayers=int(n_tf_layer),
                n_radial_array=int(n_radial_array),
            )

            # Vertical stress [Pa]
            sig_tf_z[:] = vforce / (
                acasetf + acndttf * n_tf_turn
            )  # Array equation [EDIT: Are you sure? It doesn't look like one to me]

            # Strain in vertical direction on WP
            str_wp = sig_tf_z[n_tf_bucking] / eyoung_wp_axial_eff

            # Case strain
            casestr = sig_tf_z[n_tf_bucking - 1] / eyoung_steel

            # Radial strain in insulator
            insstrain = (
                sig_tf_r[n_radial_array - 1]
                * eyoung_wp_stiffest_leg
                / eyoung_wp_trans_eff
                / eyoung_ins
            )
        # ---

        # New generalized plane strain formulation
        # ---
        # See #1670 for why i_tf_stress_model=0 was removed
        elif i_tf_stress_model in (0, 2):
            # Extended plane strain calculation [Pa]
            # Issues #1414 and #998
            # Permits build_variables.bore >= 0, O(n) in layers
            # If build_variables.bore > 0, same result as generalized plane strain calculation

            (
                radial_array,
                sig_tf_r,
                sig_tf_t,
                sig_tf_z,
                str_tf_r,
                str_tf_t,
                str_tf_z,
                deflect,
            ) = extended_plane_strain(
                poisson_trans,
                poisson_axial,
                eyoung_trans,
                eyoung_axial,
                radtf,
                jeff,
                vforce_inboard_tot,
                int(n_tf_layer),
                int(n_radial_array),
                int(n_tf_bucking),
            )

            # Strain in TF conductor material
            str_wp = str_tf_z[n_tf_bucking * n_radial_array]

        # ---

        # Storing the smeared properties for output

        sig_tf_smeared_r[:] = sig_tf_r  # Array equation
        sig_tf_smeared_t[:] = sig_tf_t  # Array equation
        sig_tf_smeared_z[:] = sig_tf_z  # Array equation

        # ------------------------------

        # STRESS DISTRIBUTIONS CORRECTIONS
        # --------------------------------
        # SC central solenoid coil stress unsmearing (bucked and wedged only)
        # ---
        if i_tf_bucking >= 2 and ipfres == 0:
            # Central Solenoid (OH) steel conduit stress unsmearing factors
            for ii in range(n_radial_array):
                sig_tf_r[ii] = sig_tf_r[ii] * eyoung_cs_stiffest_leg / eyoung_axial[0]
                sig_tf_t[ii] = sig_tf_t[ii] * eyoung_steel / eyoung_trans[0]
                sig_tf_z[ii] = sig_tf_z[ii] * eyoung_cs_stiffest_leg / eyoung_axial[0]

        # ---

        # No TF vertical forces on CS and CS-TF layer (bucked and wedged only)
        # ---
        # This correction is only applied if the plane stress model is used
        # as the generalized plane strain calculates the vertical stress properly
        if i_tf_bucking >= 2 and i_tf_stress_model == 1:
            for ii in range((n_tf_bucking - 1) * n_radial_array):
                sig_tf_z[ii] = 0.0e0

        # ---

        # Toroidal coil unsmearing
        # ---
        # Copper magnets
        if i_tf_sup == 0:
            # Vertical force unsmearing factor
            fac_sig_z = eyoung_copper / eyoung_wp_axial_eff

            # Toroidal WP steel stress unsmearing factor
            fac_sig_t = 1.0e0
            fac_sig_r = 1.0e0

        elif i_tf_sup == 1:
            # Vertical WP steel stress unsmearing factor
            if i_tf_stress_model != 1:
                fac_sig_z = eyoung_steel / eyoung_wp_axial_eff
                fac_sig_z_wp_av = eyoung_wp_axial / eyoung_wp_axial_eff
            else:
                fac_sig_z = 1.0e0

            # Toroidal WP steel conduit stress unsmearing factor
            fac_sig_t = eyoung_wp_stiffest_leg / eyoung_wp_trans_eff

            # Radial WP steel conduit stress unsmearing factor
            fac_sig_r = eyoung_wp_stiffest_leg / eyoung_wp_trans_eff

        elif i_tf_sup == 2:
            # Vertical WP steel stress unsmearing factor
            fac_sig_z = eyoung_al / eyoung_wp_axial_eff

            # Toroidal WP steel stress unsmearing factor
            # NO CALCULTED FOR THE MOMENT (to be done later)
            fac_sig_t = 1.0e0
            fac_sig_r = 1.0e0

        # Application of the unsmearing to the WP layers
        # For each point within the winding pack / conductor, unsmear the
        # stress. This is n_radial_array test points within tfcoil_variables.n_tf_graded_layers
        # layers starting at n_tf_bucking + 1
        # GRADED MODIF : add another do loop to allow the graded properties
        #                to be taken into account
        for ii in range(
            n_tf_bucking * n_radial_array,
            ((n_tf_bucking + n_tf_graded_layers) * n_radial_array),
        ):
            sig_tf_wp_av_z[ii - n_tf_bucking * n_radial_array] = (
                sig_tf_z[ii] * fac_sig_z_wp_av
            )
            sig_tf_r[ii] = sig_tf_r[ii] * fac_sig_r
            sig_tf_t[ii] = sig_tf_t[ii] * fac_sig_t
            sig_tf_z[ii] = sig_tf_z[ii] * fac_sig_z

        # For each point within the front case,
        # remove the correction for the missing axial
        # stiffness as per the updated description of
        # Issue #1509
        for ii in range(
            (n_tf_bucking + n_tf_graded_layers) * n_radial_array,
            (n_tf_layer * n_radial_array),
        ):
            sig_tf_z[ii] = sig_tf_z[ii] / f_tf_stress_front_case
        # ---

        # Tresca / Von Mises yield criteria calculations
        # -----------------------------
        # Array equation
        sig_tf_tresca_tmp1 = numpy.maximum(
            numpy.absolute(sig_tf_r - sig_tf_t), numpy.absolute(sig_tf_r - sig_tf_z)
        )
        sig_tf_tresca = numpy.maximum(
            sig_tf_tresca_tmp1, numpy.absolute(sig_tf_z - sig_tf_t)
        )

        # Array equation

        sig_tf_vmises = numpy.sqrt(
            0.5e0
            * (
                (sig_tf_r - sig_tf_t) ** 2
                + (sig_tf_r - sig_tf_z) ** 2
                + (sig_tf_z - sig_tf_t) ** 2
            )
        )

        # Array equation
        s_tresca_cond_cea = sig_tf_tresca.copy()

        # SC conducting layer stress distribution corrections
        # ---
        if i_tf_sup == 1:
            # GRADED MODIF : add another do loop to allow the graded properties
            #                to be taken into account
            for ii in range(
                (n_tf_bucking * n_radial_array), n_tf_layer * n_radial_array
            ):
                # Addaped Von-mises stress calculation to WP strucure [Pa]

                svmxz = sigvm(0.0e0, sig_tf_t[ii], sig_tf_z[ii], 0.0e0, 0.0e0, 0.0e0)

                svmyz = sigvm(sig_tf_r[ii], 0.0e0, sig_tf_z[ii], 0.0e0, 0.0e0, 0.0e0)
                sig_tf_vmises[ii] = max(svmxz, svmyz)

                # Maximum shear stress for the Tresca yield criterion using CEA calculation [Pa]
                s_tresca_cond_cea[ii] = (
                    1.02e0 * abs(sig_tf_r[ii]) + 1.6e0 * sig_tf_z[ii]
                )

        # ---
        # -----------------------------

        # Output formating : Maximum shear stress of each layer for the Tresca yield criterion
        # ----------------
        for ii in range(n_tf_layer):
            sig_max = 0.0e0
            ii_max = 0

            for jj in range(ii * n_radial_array, (ii + 1) * n_radial_array):
                # CEA out of plane approximation
                if i_tf_tresca == 1 and i_tf_sup == 1 and ii >= i_tf_bucking + 1:
                    if sig_max < s_tresca_cond_cea[jj]:
                        sig_max = s_tresca_cond_cea[jj]
                        ii_max = jj

                # Conventional Tresca
                else:
                    if sig_max < sig_tf_tresca[jj]:
                        sig_max = sig_tf_tresca[jj]
                        ii_max = jj

            # OUT.DAT output

            sig_tf_r_max[ii] = sig_tf_r[ii_max]
            sig_tf_t_max[ii] = sig_tf_t[ii_max]
            sig_tf_z_max[ii] = sig_tf_z[ii_max]
            sig_tf_vmises_max[ii] = sig_tf_vmises[ii_max]

            # Maximum shear stress for the Tresca yield criterion (or CEA OOP correction)

            if i_tf_tresca == 1 and i_tf_sup == 1 and ii >= i_tf_bucking + 1:
                sig_tf_tresca_max[ii] = s_tresca_cond_cea[ii_max]
            else:
                sig_tf_tresca_max[ii] = sig_tf_tresca[ii_max]

        # Constraint equation for the Tresca yield criterion

        sig_tf_wp = sig_tf_tresca_max[n_tf_bucking]
        # Maximum assumed in the first graded layer

        if i_tf_bucking >= 1:
            sig_tf_case = sig_tf_tresca_max[n_tf_bucking - 1]
        if i_tf_bucking >= 2:
            sig_tf_cs_bucked = sig_tf_tresca_max[0]
        # ----------------

        return (
            sig_tf_r_max,
            sig_tf_t_max,
            sig_tf_z_max,
            sig_tf_vmises_max,
            sig_tf_tresca_max,
            deflect,
            eyoung_axial,
            eyoung_trans,
            eyoung_wp_axial,
            eyoung_wp_trans,
            poisson_wp_trans,
            radial_array,
            s_tresca_cond_cea,
            poisson_wp_axial,
            sig_tf_r,
            sig_tf_smeared_r,
            sig_tf_smeared_t,
            sig_tf_smeared_z,
            sig_tf_t,
            sig_tf_tresca,
            sig_tf_vmises,
            sig_tf_z,
            str_tf_r,
            str_tf_t,
            str_tf_z,
            n_radial_array,
            n_tf_bucking,
            sig_tf_wp,
            sig_tf_case,
            sig_tf_cs_bucked,
            str_wp,
            casestr,
            insstrain,
            sig_tf_wp_av_z,
        )

    def outtf(self, peaktfflag):
        """Writes superconducting TF coil output to file
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        peaktfflag : input integer : warning flag from peak TF calculation
        This routine writes the superconducting TF coil results
        to the output file.
        PROCESS Superconducting TF Coil Model, J. Morris, CCFE, 1st May 2014
        """

        # General coil parameters
        po.osubhd(self.outfile, "TF design")
        po.ovarin(
            self.outfile,
            "Conductor technology",
            "(i_tf_sup)",
            tfcoil_variables.i_tf_sup,
        )

        if tfcoil_variables.i_tf_sup == 0:
            po.ocmmnt(
                self.outfile, "  -> Resitive coil : Water cooled copper (GLIDCOP AL-15)"
            )
        elif tfcoil_variables.i_tf_sup == 1:
            po.ocmmnt(self.outfile, "  -> Superconducting coil (SC)")
        elif tfcoil_variables.i_tf_sup == 2:
            po.ocmmnt(self.outfile, "  -> Reisitive coil : Helium cooled aluminium")

        # SC material scaling
        if tfcoil_variables.i_tf_sup == 1:
            po.ovarin(
                self.outfile,
                "Superconductor material",
                "(i_tf_sc_mat)",
                tfcoil_variables.i_tf_sc_mat,
            )

            if tfcoil_variables.i_tf_sc_mat == 1:
                po.ocmmnt(self.outfile, "  -> ITER Nb3Sn critical surface model")
            elif tfcoil_variables.i_tf_sc_mat == 2:
                po.ocmmnt(self.outfile, "  -> Bi-2212 high temperature superconductor")
            elif tfcoil_variables.i_tf_sc_mat == 3:
                po.ocmmnt(self.outfile, "  -> NbTi")
            elif tfcoil_variables.i_tf_sc_mat == 4:
                po.ocmmnt(
                    self.outfile,
                    "  -> ITER Nb3Sn critical surface model, user-defined parameters",
                )
            elif tfcoil_variables.i_tf_sc_mat == 5:
                po.ocmmnt(self.outfile, "  -> WST Nb3Sn")
            elif tfcoil_variables.i_tf_sc_mat == 6:
                po.ocmmnt(
                    self.outfile,
                    "  -> High temperature superconductor: REBCO HTS tape in CroCo strand",
                )
            elif tfcoil_variables.i_tf_sc_mat == 7:
                po.ocmmnt(
                    self.outfile,
                    "  ->  Durham Ginzburg-Landau critical surface model for Nb-Ti",
                )
            elif tfcoil_variables.i_tf_sc_mat == 8:
                po.ocmmnt(
                    self.outfile,
                    "  ->  Durham Ginzburg-Landau critical surface model for REBCO",
                )
            elif tfcoil_variables.i_tf_sc_mat == 9:
                po.ocmmnt(
                    self.outfile,
                    "  ->  Hazelton experimental data + Zhai conceptual model for REBCO",
                )

        # Joints strategy
        po.ovarin(
            self.outfile,
            "Presence of TF demountable joints",
            "(itart)",
            physics_variables.itart,
        )
        if physics_variables.itart == 1:
            po.ocmmnt(
                self.outfile, "  -> TF coil made of a Centerpost (CP) and outer legs"
            )
            po.ocmmnt(self.outfile, "     interfaced with demountable joints")
        else:
            po.ocmmnt(self.outfile, "  -> Coils without demountable joints")

        # Centring forces support strategy
        po.ovarin(
            self.outfile,
            "TF inboard leg support strategy",
            "(i_tf_bucking)",
            tfcoil_variables.i_tf_bucking,
        )

        if tfcoil_variables.i_tf_bucking == 0:
            po.ocmmnt(self.outfile, "  -> No support structure")
        elif tfcoil_variables.i_tf_bucking == 1:
            if tfcoil_variables.i_tf_sup == 1:
                po.ocmmnt(self.outfile, "  -> Steel casing")
            elif (
                abs(tfcoil_variables.eyoung_res_tf_buck - 205.0e9)
                < numpy.finfo(float(tfcoil_variables.eyoung_res_tf_buck)).eps
            ):
                po.ocmmnt(self.outfile, "  -> Steel bucking cylinder")
            else:
                po.ocmmnt(self.outfile, "  -> Bucking cylinder")

        elif tfcoil_variables.i_tf_bucking in (2, 3) and build_variables.tf_in_cs == 1:
            po.ocmmnt(
                self.outfile,
                "  -> TF in contact with bore filler support (bucked and weged design)",
            )

        elif tfcoil_variables.i_tf_bucking in (2, 3) and build_variables.tf_in_cs == 0:
            po.ocmmnt(
                self.outfile, "  -> TF in contact with CS (bucked and weged design)"
            )

        # TF coil geometry
        po.osubhd(self.outfile, "TF coil Geometry :")
        po.ovarin(
            self.outfile, "Number of TF coils", "(n_tf)", int(tfcoil_variables.n_tf)
        )
        po.ovarre(
            self.outfile,
            "Inboard leg centre radius (m)",
            "(r_tf_inboard_mid)",
            build_variables.r_tf_inboard_mid,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Outboard leg centre radius (m)",
            "(r_tf_outboard_mid)",
            build_variables.r_tf_outboard_mid,
            "OP ",
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
            "Outboard leg toroidal thickness (m)",
            "(tftort)",
            tfcoil_variables.tftort,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Maximum inboard edge height (m)",
            "(hmax)",
            build_variables.hmax,
            "OP ",
        )
        if physics_variables.itart == 1:
            po.ovarre(
                self.outfile,
                "Mean coil circumference (inboard leg not included) (m)",
                "(tfleng)",
                tfcoil_variables.tfleng,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Length of the inboard segment (m)",
                "(cplen)",
                tfcoil_variables.cplen,
                "OP ",
            )
        else:
            po.ovarre(
                self.outfile,
                "Mean coil circumference (including inboard leg length) (m)",
                "(tfleng)",
                tfcoil_variables.tfleng,
                "OP ",
            )

        # Vertical shape
        po.ovarin(
            self.outfile,
            "Vertical TF shape",
            "(i_tf_shape)",
            tfcoil_variables.i_tf_shape,
        )
        if tfcoil_variables.i_tf_shape == 1:
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "D-shape coil, inner surface shape approximated by")
            po.ocmmnt(
                self.outfile,
                "by a straight segment and elliptical arcs between the following points:",
            )
            po.oblnkl(self.outfile)
        elif tfcoil_variables.i_tf_shape == 2:
            po.oblnkl(self.outfile)
            po.ocmmnt(self.outfile, "Picture frame coil, inner surface approximated by")
            po.ocmmnt(
                self.outfile, "by a straight segment between the following points:"
            )
            po.oblnkl(self.outfile)

        po.write(self.outfile, "  point              x(m)              y(m)")
        for ii in range(5):
            po.write(
                self.outfile,
                f"  {ii}              {tfcoil_variables.xarc[ii]}              {tfcoil_variables.yarc[ii]}",
            )
            po.ovarre(
                constants.mfile,
                f"TF coil arc point {ii} R (m)",
                f"(xarc({ii+1}))",
                tfcoil_variables.xarc[ii],
            )
            po.ovarre(
                constants.mfile,
                f"TF coil arc point {ii} Z (m)",
                f"(yarc({ii+1}))",
                tfcoil_variables.yarc[ii],
            )

        # CP tapering geometry
        if physics_variables.itart == 1 and tfcoil_variables.i_tf_sup != 1:
            po.osubhd(self.outfile, "Tapered Centrepost TF coil Dimensions:")
            po.ovarre(
                self.outfile,
                "TF coil centrepost outer radius at midplane (m)",
                "(r_tf_inboard_out)",
                build_variables.r_tf_inboard_out,
            )
            po.ovarre(
                self.outfile,
                "TF coil centrepost outer radius at its top (m)",
                "(r_cp_top)",
                build_variables.r_cp_top,
            )
            po.ovarre(
                self.outfile,
                "Top/miplane TF CP radius ratio (-)",
                "(f_r_cp)",
                build_variables.f_r_cp,
            )
            po.ovarre(
                self.outfile,
                "Distance from the midplane to the top of the tapered section (m)",
                "(h_cp_top)",
                sctfcoil_module.h_cp_top,
            )
            po.ovarre(
                self.outfile,
                "Distance from the midplane to the top of the centrepost (m)",
                "(hmax + tfthko)",
                build_variables.hmax + build_variables.tfthko,
            )

        # Turn/WP gemoetry
        if tfcoil_variables.i_tf_sup == 1:
            # Total material fraction
            po.osubhd(self.outfile, "Global material area/fractions:")
            po.ovarre(
                self.outfile,
                "TF cross-section (total) (m2)",
                "(tfareain)",
                tfcoil_variables.tfareain,
            )
            po.ovarre(
                self.outfile,
                "Total steel cross-section (m2)",
                "(a_tf_steel*n_tf)",
                sctfcoil_module.a_tf_steel * tfcoil_variables.n_tf,
            )
            po.ovarre(
                self.outfile,
                "Total steel TF fraction",
                "(f_tf_steel)",
                sctfcoil_module.f_tf_steel,
            )
            po.ovarre(
                self.outfile,
                "Total Insulation cross-section (total) (m2)",
                "(a_tf_ins*n_tf)",
                sctfcoil_module.a_tf_ins * tfcoil_variables.n_tf,
            )
            po.ovarre(
                self.outfile,
                "Total Insulation fraction",
                "(f_tf_ins)",
                sctfcoil_module.f_tf_ins,
            )

            # External casing
            po.osubhd(self.outfile, "External steel Case Information :")
            po.ovarre(
                self.outfile,
                "Casing cross section area (per leg) (m2)",
                "(acasetf)",
                tfcoil_variables.acasetf,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg case plasma side wall thickness (m)",
                "(casthi)",
                tfcoil_variables.casthi,
            )
            po.ovarre(
                self.outfile,
                'Inboard leg case inboard "nose" thickness (m)',
                "(thkcas)",
                tfcoil_variables.thkcas,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg case sidewall thickness at its narrowest point (m)",
                "(casths)",
                tfcoil_variables.casths,
            )
            po.ovarre(
                self.outfile,
                "External case mass per coil (kg)",
                "(whtcas)",
                tfcoil_variables.whtcas,
                "OP ",
            )

            # Winding pack structure
            po.osubhd(self.outfile, "TF winding pack (WP) geometry:")
            po.ovarre(
                self.outfile,
                "WP cross section area with insulation and insertion (per coil) (m2)",
                "(awpc)",
                sctfcoil_module.awpc,
            )
            po.ovarre(
                self.outfile,
                "WP cross section area (per coil) (m2)",
                "(aswp)",
                sctfcoil_module.awptf,
            )
            po.ovarre(
                self.outfile,
                "Winding pack radial thickness (m)",
                "(dr_tf_wp)",
                tfcoil_variables.dr_tf_wp,
                "OP ",
            )
            if tfcoil_variables.i_tf_turns_integer == 1:
                po.ovarre(
                    self.outfile,
                    "Winding pack toroidal width (m)",
                    "(wwp1)",
                    tfcoil_variables.wwp1,
                    "OP ",
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Winding pack toroidal width 1 (m)",
                    "(wwp1)",
                    tfcoil_variables.wwp1,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Winding pack toroidal width 2 (m)",
                    "(wwp2)",
                    tfcoil_variables.wwp2,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Ground wall insulation thickness (m)",
                "(tinstf)",
                tfcoil_variables.tinstf,
            )
            po.ovarre(
                self.outfile,
                "Winding pack insertion gap (m)",
                "(tfinsgap)",
                tfcoil_variables.tfinsgap,
            )

            # WP material fraction
            po.osubhd(self.outfile, "TF winding pack (WP) material area/fractions:")
            po.ovarre(
                self.outfile,
                "Steel WP cross-section (total) (m2)",
                "(aswp*n_tf)",
                tfcoil_variables.aswp * tfcoil_variables.n_tf,
            )
            po.ovarre(
                self.outfile,
                "Steel WP fraction",
                "(aswp/awpc)",
                tfcoil_variables.aswp / sctfcoil_module.awpc,
            )
            po.ovarre(
                self.outfile,
                "Insulation WP fraction",
                "(aiwp/awpc)",
                tfcoil_variables.aiwp / sctfcoil_module.awpc,
            )
            po.ovarre(
                self.outfile,
                "Cable WP fraction",
                "((awpc-aswp-aiwp)/awpc)",
                (sctfcoil_module.awpc - tfcoil_variables.aswp - tfcoil_variables.aiwp)
                / sctfcoil_module.awpc,
            )

            # Number of turns
            po.osubhd(self.outfile, "WP turn information:")
            po.ovarin(
                self.outfile,
                "Turn parametrisation",
                "(i_tf_turns_integer)",
                tfcoil_variables.i_tf_turns_integer,
            )
            if tfcoil_variables.i_tf_turns_integer == 0:
                po.ocmmnt(self.outfile, "  Non-integer number of turns")
            else:
                po.ocmmnt(self.outfile, "  Integer number of turns")

            po.ovarre(
                self.outfile,
                "Number of turns per TF coil",
                "(n_tf_turn)",
                tfcoil_variables.n_tf_turn,
                "OP ",
            )
            if tfcoil_variables.i_tf_turns_integer == 1:
                po.ovarin(
                    self.outfile,
                    "Number of TF pancakes",
                    "(n_pancake)",
                    tfcoil_variables.n_pancake,
                )
                po.ovarin(
                    self.outfile,
                    "Number of TF layers",
                    "(n_layer)",
                    tfcoil_variables.n_layer,
                )

            po.oblnkl(self.outfile)

            if tfcoil_variables.i_tf_turns_integer == 1:
                po.ovarre(
                    self.outfile,
                    "Radial width of turn (m)",
                    "(t_turn_radial)",
                    sctfcoil_module.t_turn_radial,
                )
                po.ovarre(
                    self.outfile,
                    "Toroidal width of turn (m)",
                    "(t_turn_toroidal)",
                    sctfcoil_module.t_turn_toroidal,
                )
                po.ovarre(
                    self.outfile,
                    "Radial width of conductor (m)",
                    "(elonductor_radial)",
                    sctfcoil_module.t_conductor_radial,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Toroidal width of conductor (m)",
                    "(t_conductor_toroidal)",
                    sctfcoil_module.t_conductor_toroidal,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Radial width of cable space",
                    "(t_cable_radial)",
                    sctfcoil_module.t_cable_radial,
                )
                po.ovarre(
                    self.outfile,
                    "Toroidal width of cable space",
                    "(t_cable_toroidal)",
                    sctfcoil_module.t_cable_toroidal,
                )
            else:
                po.ovarre(
                    self.outfile,
                    "Width of turn including inter-turn insulation (m)",
                    "(t_turn_tf)",
                    tfcoil_variables.t_turn_tf,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Width of conductor (square) (m)",
                    "(t_conductor)",
                    tfcoil_variables.t_conductor,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Width of space inside conductor (m)",
                    "(t_cable)",
                    sctfcoil_module.t_cable,
                    "OP ",
                )

            po.ovarre(
                self.outfile,
                "Steel conduit thickness (m)",
                "(thwcndut)",
                tfcoil_variables.thwcndut,
            )
            po.ovarre(
                self.outfile,
                "Inter-turn insulation thickness (m)",
                "(thicndut)",
                tfcoil_variables.thicndut,
            )

            if tfcoil_variables.i_tf_sc_mat in (1, 2, 3, 4, 5, 7, 8, 9):
                po.osubhd(self.outfile, "Conductor information:")
                po.ovarre(
                    self.outfile,
                    "Diameter of central helium channel in cable",
                    "(dhecoil)",
                    tfcoil_variables.dhecoil,
                )
                po.ocmmnt(self.outfile, "Fractions by area")
                po.ovarre(
                    self.outfile,
                    "internal area of the cable space",
                    "(acstf)",
                    tfcoil_variables.acstf,
                )
                po.ovarre(
                    self.outfile,
                    "Coolant fraction in conductor excluding central channel",
                    "(vftf)",
                    tfcoil_variables.vftf,
                )
                po.ovarre(
                    self.outfile,
                    "Copper fraction of conductor",
                    "(fcutfsu)",
                    tfcoil_variables.fcutfsu,
                )
                po.ovarre(
                    self.outfile,
                    "Superconductor fraction of conductor",
                    "(1-fcutfsu)",
                    1 - tfcoil_variables.fcutfsu,
                )
                # TODO
                # po.ovarre(self.outfile,'Conductor fraction of winding pack','(tfcoil_variables.acond/ap)',acond/ap, 'OP ')
                # po.ovarre(self.outfile,'Conduit fraction of winding pack','(tfcoil_variables.n_tf_turn*tfcoil_variables.acndttf/ap)',n_tf_turn*tfcoil_variables.acndttf/ap, 'OP ')
                # po.ovarre(self.outfile,'Insulator fraction of winding pack','(tfcoil_variables.aiwp/ap)',aiwp/ap, 'OP ')
                # po.ovarre(self.outfile,'Helium area fraction of winding pack excluding central channel','(tfcoil_variables.avwp/ap)',avwp/ap, 'OP ')
                # po.ovarre(self.outfile,'Central helium channel area as fraction of winding pack','(tfcoil_variables.awphec/ap)',awphec/ap, 'OP ')
                ap = (
                    tfcoil_variables.acond
                    + tfcoil_variables.n_tf_turn * tfcoil_variables.acndttf
                    + tfcoil_variables.aiwp
                    + tfcoil_variables.avwp
                    + tfcoil_variables.awphec
                )
                po.ovarrf(
                    self.outfile,
                    "Check total area fractions in winding pack = 1",
                    "",
                    (
                        tfcoil_variables.acond
                        + tfcoil_variables.n_tf_turn * tfcoil_variables.acndttf
                        + tfcoil_variables.aiwp
                        + tfcoil_variables.avwp
                        + tfcoil_variables.awphec
                    )
                    / ap,
                )
                po.ovarrf(
                    self.outfile,
                    "minimum TF conductor temperature margin  (K)",
                    "(tmargmin_tf)",
                    tfcoil_variables.tmargmin_tf,
                )
                po.ovarrf(
                    self.outfile,
                    "TF conductor temperature margin (K)",
                    "(tmargtf)",
                    tfcoil_variables.tmargtf,
                )

                po.ovarin(
                    self.outfile,
                    "Elastic properties behavior",
                    "(i_tf_cond_eyoung_axial)",
                    tfcoil_variables.i_tf_cond_eyoung_axial,
                )
                if tfcoil_variables.i_tf_cond_eyoung_axial == 0:
                    po.ocmmnt(self.outfile, "  Conductor stiffness neglected")
                elif tfcoil_variables.i_tf_cond_eyoung_axial == 1:
                    po.ocmmnt(self.outfile, "  Conductor stiffness is user-input")
                elif tfcoil_variables.i_tf_cond_eyoung_axial == 2:
                    po.ocmmnt(
                        self.outfile,
                        "  Conductor stiffness is set by material-specific default",
                    )

                po.ovarre(
                    self.outfile,
                    "Conductor axial Young" "s modulus",
                    "(eyoung_cond_axial)",
                    tfcoil_variables.eyoung_cond_axial,
                )
                po.ovarre(
                    self.outfile,
                    "Conductor transverse Young" "s modulus",
                    "(eyoung_cond_trans)",
                    tfcoil_variables.eyoung_cond_trans,
                )
        else:
            # External casing
            po.osubhd(self.outfile, "Bucking cylinder information:")
            po.ovarre(
                self.outfile,
                "Casing cross section area (per leg) (m2)",
                "(acasetf)",
                tfcoil_variables.acasetf,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg case plasma side wall thickness (m)",
                "(casthi)",
                tfcoil_variables.casthi,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg bucking cylinder thickness (m)",
                "(thkcas)",
                tfcoil_variables.thkcas,
            )

            # Conductor layer geometry
            po.osubhd(self.outfile, "Inboard TFC conductor sector geometry:")
            po.ovarre(
                self.outfile,
                "Inboard TFC conductor sector area with gr insulation (per leg) (m2)",
                "(awpc)",
                sctfcoil_module.awpc,
            )
            po.ovarre(
                self.outfile,
                "Inboard TFC conductor sector area, NO ground & gap (per leg) (m2)",
                "(awptf)",
                sctfcoil_module.awptf,
            )
            po.ovarre(
                self.outfile,
                "Inboard conductor sector radial thickness (m)",
                "(dr_tf_wp)",
                tfcoil_variables.dr_tf_wp,
            )
            if physics_variables.itart == 1:
                po.ovarre(
                    self.outfile,
                    "Central collumn top conductor sector radial thickness (m)",
                    "(dr_tf_wp_top)",
                    sctfcoil_module.dr_tf_wp_top,
                )

            po.ovarre(
                self.outfile,
                "Ground wall insulation thickness (m)",
                "(tinstf)",
                tfcoil_variables.tinstf,
            )
            # Turn info
            po.osubhd(self.outfile, "Coil turn information:")
            po.ovarre(
                self.outfile,
                "Number of turns per TF leg",
                "(n_tf_turn)",
                tfcoil_variables.n_tf_turn,
            )
            po.ovarre(
                self.outfile,
                "Turn insulation thickness",
                "(thicndut)",
                tfcoil_variables.thicndut,
            )
            po.ovarre(
                self.outfile,
                "Mid-plane CP cooling fraction",
                "(fcoolcp)",
                tfcoil_variables.fcoolcp,
            )
            po.ovarre(
                self.outfile,
                "Outboard leg current per turn (A)",
                "(cpttf)",
                tfcoil_variables.cpttf,
            )
            po.ovarre(
                self.outfile,
                "Inboard leg conductor volume (m3)",
                "(vol_cond_cp)",
                tfcoil_variables.vol_cond_cp,
            )
            po.ovarre(
                self.outfile,
                "Outboard leg volume per coil (m3)",
                "(voltfleg)",
                tfcoil_variables.voltfleg,
            )

        # Coil masses
        po.osubhd(self.outfile, "TF coil mass:")
        po.ovarre(
            self.outfile,
            "Superconductor mass per coil (kg)",
            "(whtconsc)",
            tfcoil_variables.whtconsc,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Copper mass per coil (kg)",
            "(whtconcu)",
            tfcoil_variables.whtconcu,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Steel conduit mass per coil (kg)",
            "(whtconsh)",
            tfcoil_variables.whtconsh,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Conduit insulation mass per coil (kg)",
            "(whtconin)",
            tfcoil_variables.whtconin,
            "OP ",
        )
        if tfcoil_variables.i_tf_sup == 1:
            po.ovarre(
                self.outfile,
                "Total conduit mass per coil (kg)",
                "(whtcon)",
                tfcoil_variables.whtcon,
                "OP ",
            )

        if physics_variables.itart == 1:
            po.ovarre(
                self.outfile,
                "Mass of inboard legs (kg)",
                "(whtcp)",
                tfcoil_variables.whtcp,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Mass of outboard legs (kg)",
                "(whttflgs)",
                tfcoil_variables.whttflgs,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Mass of each TF coil (kg)",
            "(whttf/n_tf)",
            tfcoil_variables.whttf / tfcoil_variables.n_tf,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total TF coil mass (kg)",
            "(whttf)",
            tfcoil_variables.whttf,
            "OP ",
        )

        # TF current and field
        po.osubhd(self.outfile, "Maximum B field and currents:")
        po.ovarre(
            self.outfile,
            "Nominal peak field assuming toroidal symmetry (T)",
            "(bmaxtf)",
            tfcoil_variables.bmaxtf,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total current in all TF coils (MA)",
            "(ritfc/1.D6)",
            1.0e-6 * tfcoil_variables.ritfc,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "TF coil current (summed over all coils) (A)",
            "(ritfc)",
            tfcoil_variables.ritfc,
        )
        if tfcoil_variables.i_tf_sup == 1:
            po.ovarre(
                self.outfile,
                "Actual peak field at discrete conductor (T)",
                "(bmaxtfrp)",
                tfcoil_variables.bmaxtfrp,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Winding pack current density (A/m2)",
                "(jwptf)",
                tfcoil_variables.jwptf,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Inboard leg mid-plane conductor current density (A/m2)",
            "(oacdcp)",
            tfcoil_variables.oacdcp,
        )
        if physics_variables.itart == 1:
            po.ovarre(
                self.outfile,
                "Outboard leg conductor current density (A/m2)",
                "(cdtfleg)",
                tfcoil_variables.cdtfleg,
            )

        po.ovarre(
            self.outfile,
            "Total stored energy in TF coils (GJ)",
            "(estotftgj)",
            tfcoil_variables.estotftgj,
            "OP ",
        )

        # TF forces
        po.osubhd(self.outfile, "TF Forces:")
        po.ovarre(
            self.outfile,
            "Inboard vertical tension per coil (N)",
            "(vforce)",
            tfcoil_variables.vforce,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Outboard vertical tension per coil (N)",
            "(vforce_outboard)",
            tfcoil_variables.vforce_outboard,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "inboard vertical tension fraction (-)",
            "(f_vforce_inboard)",
            tfcoil_variables.f_vforce_inboard,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Centring force per coil (N/m)",
            "(cforce)",
            tfcoil_variables.cforce,
            "OP ",
        )

        # Resistive coil parameters
        if tfcoil_variables.i_tf_sup != 1:
            po.osubhd(self.outfile, "Resitive loss parameters:")
            if tfcoil_variables.i_tf_sup == 0:
                po.ocmmnt(
                    self.outfile,
                    "Resistive Material : GLIDCOP AL-15 - Dispersion Strengthened Copper",
                )
            elif tfcoil_variables.i_tf_sup == 2:
                po.ocmmnt(
                    self.outfile, "Resistive Material : Pure Aluminium (99.999+ %)"
                )

            if physics_variables.itart == 1:
                po.ovarre(
                    self.outfile,
                    "CP resistivity (ohm.m)",
                    "(rhocp)",
                    tfcoil_variables.rhocp,
                )
                po.ovarre(
                    self.outfile,
                    "Leg resistivity (ohm.m)",
                    "(rhotfleg)",
                    tfcoil_variables.rhotfleg,
                )
                po.ovarre(
                    self.outfile,
                    "CP resistive power loss (W)",
                    "(prescp)",
                    tfcoil_variables.prescp,
                )
                po.ovarre(
                    self.outfile,
                    "Leg resitive power loss, (per leg) (W)",
                    "(presleg)",
                    tfcoil_variables.presleg,
                )
                po.ovarre(
                    self.outfile,
                    "joints resistive power loss (W)",
                    "(pres_joints)",
                    tfcoil_variables.pres_joints,
                )
                po.ovarre(
                    self.outfile,
                    "Outboard leg resistance per coil (ohm)",
                    "(tflegres)",
                    tfcoil_variables.tflegres,
                )
                po.ovarre(
                    self.outfile,
                    "Average CP temperature (K)",
                    "(tcpav)",
                    tfcoil_variables.tcpav,
                )
                po.ovarre(
                    self.outfile,
                    "Average leg temperature (K)",
                    "(tlegav)",
                    tfcoil_variables.tlegav,
                )

            else:
                po.ovarre(
                    self.outfile,
                    "TF resistivity (ohm.m)",
                    "(prescp)",
                    tfcoil_variables.rhocp,
                )
                po.ovarre(
                    self.outfile,
                    "TF coil resistive power less (total) (ohm.m)",
                    "(prescp)",
                    tfcoil_variables.prescp,
                )
                po.ovarre(
                    self.outfile,
                    "Average coil temperature (K)",
                    "(tcpav)",
                    tfcoil_variables.tcpav,
                )

        # Ripple calculations
        po.osubhd(self.outfile, "Ripple information:")
        if tfcoil_variables.i_tf_shape == 1:
            if peaktfflag == 1:
                error_handling.report_error(144)
            elif peaktfflag == 2:
                error_handling.report_error(145)

            po.ovarre(
                self.outfile,
                "Max allowed tfcoil_variables.ripple amplitude at plasma outboard midplane (%)",
                "(ripmax)",
                tfcoil_variables.ripmax,
            )
            po.ovarre(
                self.outfile,
                "Ripple amplitude at plasma outboard midplane (%)",
                "(ripple)",
                tfcoil_variables.ripple,
                "OP ",
            )
        else:
            po.ovarre(
                self.outfile,
                "Max allowed tfcoil_variables.ripple amplitude at plasma (%)",
                "(ripmax)",
                tfcoil_variables.ripmax,
            )
            po.ovarre(
                self.outfile,
                "Ripple at plasma edge (%)",
                "(ripple)",
                tfcoil_variables.ripple,
            )
            po.ocmmnt(
                self.outfile,
                "  Ripple calculation to be re-defined for picure frame coils",
            )

        # Quench information
        if tfcoil_variables.i_tf_sup == 1:
            po.osubhd(self.outfile, "Quench information :")
            po.ovarre(
                self.outfile,
                "Actual quench time (or time constant) (s)",
                "(tdmptf)",
                tfcoil_variables.tdmptf,
            )
            po.ovarre(
                self.outfile,
                "Vacuum Vessel stress on quench (Pa)",
                "(vv_stress_quench)",
                sctfcoil_module.vv_stress_quench,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Maximum allowed voltage during quench due to insulation (kV)",
                "(vdalw)",
                tfcoil_variables.vdalw,
            )
            po.ovarre(
                self.outfile,
                "Actual quench voltage (kV)",
                "(vtfskv)",
                tfcoil_variables.vtfskv,
                "OP ",
            )

            if tfcoil_variables.i_tf_sc_mat in (1, 2, 3, 4, 5):
                po.ovarre(
                    self.outfile,
                    "Maximum allowed temp rise during a quench (K)",
                    "(tmaxpro)",
                    tfcoil_variables.tmaxpro,
                )
            elif tfcoil_variables == 6:
                po.ocmmnt(self.outfile, "CroCo cable with jacket: ")

                if 75 in numerics.icc:
                    po.ovarre(
                        self.outfile,
                        "Maximum permitted TF coil current / copper area (A/m2)",
                        "(copperA_m2_max)",
                        rebco_variables.copperA_m2_max,
                    )

                po.ovarre(
                    self.outfile,
                    "Actual TF coil current / copper area (A/m2)",
                    "(copperA_m2)",
                    rebco_variables.copperA_m2,
                )

        # TF coil radial build
        po.osubhd(self.outfile, "Radial build of TF coil centre-line :")
        # po.write(self.outfile,5)
        # 5   format(t43,'Thickness (m)',t60,'Outer radius (m)')

        radius = build_variables.r_tf_inboard_in
        po.obuild(self.outfile, "Innermost edge of TF coil", radius, radius)

        # Radial build for SC TF coils
        if tfcoil_variables.i_tf_sup == 1:
            radius = radius + tfcoil_variables.thkcas
            po.obuild(
                self.outfile,
                'Coil case ("nose")',
                tfcoil_variables.thkcas,
                radius,
                "(thkcas)",
            )

            radius = radius + tfcoil_variables.tfinsgap
            po.obuild(
                self.outfile,
                "Insertion gap for winding pack",
                tfcoil_variables.tfinsgap,
                radius,
                "(tfinsgap)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Winding pack ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = (
                radius
                + 0.5e0 * tfcoil_variables.dr_tf_wp
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap
            )
            po.obuild(
                self.outfile,
                "Winding - first half",
                tfcoil_variables.dr_tf_wp / 2e0
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap,
                radius,
                "(dr_tf_wp/2-tinstf-tfinsgap)",
            )

            radius = (
                radius
                + 0.5e0 * tfcoil_variables.dr_tf_wp
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap
            )
            po.obuild(
                self.outfile,
                "Winding - second half",
                tfcoil_variables.dr_tf_wp / 2e0
                - tfcoil_variables.tinstf
                - tfcoil_variables.tfinsgap,
                radius,
                "(dr_tf_wp/2-tinstf-tfinsgap)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Winding pack insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = radius + tfcoil_variables.tfinsgap
            po.obuild(
                self.outfile,
                "Insertion gap for winding pack",
                tfcoil_variables.tfinsgap,
                radius,
                "(tfinsgap)",
            )

            radius = radius + tfcoil_variables.casthi
            po.obuild(
                self.outfile,
                "Plasma side case min radius",
                tfcoil_variables.casthi,
                radius,
                "(casthi)",
            )

            radius = radius / numpy.cos(numpy.pi / tfcoil_variables.n_tf)
            po.obuild(
                self.outfile,
                "Plasma side case max radius",
                build_variables.r_tf_inboard_out,
                radius,
                "(r_tf_inboard_out)",
            )

        # Radial build for restive coil
        else:
            radius = radius + tfcoil_variables.thkcas
            po.obuild(
                self.outfile,
                "Coil bucking cylindre",
                tfcoil_variables.thkcas,
                radius,
                "(thkcas)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = (
                radius + 0.5e0 * tfcoil_variables.dr_tf_wp - tfcoil_variables.tinstf
            )
            po.obuild(
                self.outfile,
                "Conductor - first half",
                tfcoil_variables.dr_tf_wp / 2e0 - tfcoil_variables.tinstf,
                radius,
                "(tfcoil_variables.dr_tf_wp/2-tfcoil_variables.tinstf)",
            )

            radius = (
                radius + 0.5e0 * tfcoil_variables.dr_tf_wp - tfcoil_variables.tinstf
            )
            po.obuild(
                self.outfile,
                "Conductor - second half",
                tfcoil_variables.dr_tf_wp / 2e0 - tfcoil_variables.tinstf,
                radius,
                "(tfcoil_variables.dr_tf_wp/2-tfcoil_variables.tinstf)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = radius + tfcoil_variables.casthi
            po.obuild(
                self.outfile,
                "Plasma side TF coil support",
                tfcoil_variables.casthi,
                radius,
                "(casthi)",
            )

        # Radial build consistency check
        if (
            abs(radius - build_variables.r_tf_inboard_in - build_variables.tfcth)
            < 10.0e0 * numpy.finfo(float(radius)).eps
        ):
            po.ocmmnt(self.outfile, "TF coil dimensions are consistent")
        else:
            po.ocmmnt(self.outfile, "ERROR: TF coil dimensions are NOT consistent:")
            po.ovarre(
                self.outfile,
                "Radius of plasma-facing side of inner leg SHOULD BE [m]",
                "",
                build_variables.r_tf_inboard_in + build_variables.tfcth,
            )
            po.ovarre(
                self.outfile,
                "Inboard TF coil radial thickness [m]",
                "(tfcth)",
                build_variables.tfcth,
            )
            po.oblnkl(self.outfile)

        # Top section TF coil radial build (physics_variables.itart = 1 only)
        if physics_variables.itart == 1 and tfcoil_variables.i_tf_sup != 1:
            po.osubhd(self.outfile, "Radial build of TF coil at central collumn top :")
            # write(self.outfile,5)

            # Restart the radial build at bucking cylindre inner radius
            radius = build_variables.r_tf_inboard_in
            po.obuild(self.outfile, "Innermost edge of TF coil", radius, radius)

            radius = radius + tfcoil_variables.thkcas
            po.obuild(
                self.outfile,
                "Coil bucking cylindre",
                tfcoil_variables.thkcas,
                radius,
                "(thkcas)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = (
                radius + 0.5e0 * sctfcoil_module.dr_tf_wp_top - tfcoil_variables.tinstf
            )
            po.obuild(
                self.outfile,
                "Conductor - first half",
                0.5e0 * sctfcoil_module.dr_tf_wp_top - tfcoil_variables.tinstf,
                radius,
                "(dr_tf_wp_top/2-tinstf)",
            )

            radius = (
                radius + 0.5e0 * sctfcoil_module.dr_tf_wp_top - tfcoil_variables.tinstf
            )
            po.obuild(
                self.outfile,
                "Conductor - second half",
                0.5e0 * sctfcoil_module.dr_tf_wp_top - tfcoil_variables.tinstf,
                radius,
                "(dr_tf_wp_top/2-tinstf)",
            )

            radius = radius + tfcoil_variables.tinstf
            po.obuild(
                self.outfile,
                "Conductor ground insulation",
                tfcoil_variables.tinstf,
                radius,
                "(tinstf)",
            )

            radius = radius + tfcoil_variables.casthi
            po.obuild(
                self.outfile,
                "Plasma side TF coil support",
                tfcoil_variables.casthi,
                radius,
                "(casthi)",
            )

            # Consistency check
            if abs(radius - build_variables.r_cp_top) < numpy.finfo(float(radius)).eps:
                po.ocmmnt(self.outfile, "Top TF coil dimensions are consistent")
            else:
                po.ocmmnt(self.outfile, "ERROR: TF coil dimensions are NOT consistent:")
                po.ovarre(
                    self.outfile,
                    "Radius of plasma-facing side of inner leg SHOULD BE [m]",
                    "",
                    build_variables.r_cp_top,
                )
                po.oblnkl(self.outfile)

    def out_stress(
        self,
        sig_tf_r_max,
        sig_tf_t_max,
        sig_tf_z_max,
        sig_tf_vmises_max,
        sig_tf_tresca_max,
        deflect,
        eyoung_axial,
        eyoung_trans,
        eyoung_wp_axial,
        eyoung_wp_trans,
        poisson_wp_trans,
        radial_array,
        s_tresca_cond_cea,
        poisson_wp_axial,
        sig_tf_r,
        sig_tf_smeared_r,
        sig_tf_smeared_t,
        sig_tf_smeared_z,
        sig_tf_t,
        sig_tf_tresca,
        sig_tf_vmises,
        sig_tf_z,
        str_tf_r,
        str_tf_t,
        str_tf_z,
        n_radial_array,
        n_tf_bucking,
        sig_tf_wp_av_z,
    ):
        """Subroutine showing the writing the TF midplane stress analysis
        in the output file and the stress distribution in the SIG_TF.json
        file used to plot stress distributions
        Author : S. Kahn
        """

        def table_format_arrays(a, mult=1, delim="\t\t"):
            return delim.join([str(i * mult) for i in a])

        # Stress output section
        po.oheadr(self.outfile, "TF coils ")
        po.osubhd(self.outfile, "TF Coil Stresses (CCFE model) :")

        if tfcoil_variables.i_tf_stress_model == 1:
            po.ocmmnt(self.outfile, "Plane stress model with smeared properties")
        else:
            po.ocmmnt(self.outfile, "Generalized plane strain model")

        po.ovarre(
            self.outfile,
            "Allowable maximum shear stress in TF coil case (Tresca criterion) (Pa)",
            "(sig_tf_case_max)",
            tfcoil_variables.sig_tf_case_max,
        )

        po.ovarre(
            self.outfile,
            "Allowable maximum shear stress in TF coil conduit (Tresca criterion) (Pa)",
            "(sig_tf_wp_max)",
            tfcoil_variables.sig_tf_wp_max,
        )
        if tfcoil_variables.i_tf_tresca == 1 and tfcoil_variables.i_tf_sup == 1:
            po.ocmmnt(
                self.outfile,
                "WP conduit Tresca criterion corrected using CEA formula (i_tf_tresca = 1)",
            )

        if tfcoil_variables.i_tf_bucking >= 3:
            po.ocmmnt(
                self.outfile, "No stress limit imposed on the TF-CS interface layer"
            )
            po.ocmmnt(
                self.outfile, "  -> Too much unknow on it material choice/properties"
            )

        # OUT.DAT data on maximum shear stress values for the Tresca criterion
        po.ocmmnt(
            self.outfile,
            "Materal stress of the point of maximum shear stress (Tresca criterion) for each layer",
        )
        po.ocmmnt(
            self.outfile,
            "Please use utilities/plot_stress_tf.py for radial plots plots summary",
        )

        if tfcoil_variables.i_tf_bucking == 0:
            if tfcoil_variables.i_tf_sup == 1:
                po.write(self.outfile, "  Layers \t\t\t\t WP \t\t Outer case")
            else:
                po.write(self.outfile, "  Layers \t\t\t\t conductor \t\t Outer case")

        elif tfcoil_variables.i_tf_bucking == 1:
            if tfcoil_variables.i_tf_sup == 1:
                po.write(
                    self.outfile, "  Layers \t\t\t\t Steel case \t\t WP \t\t Outer case"
                )
            else:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t bucking \t\t conductor \t\t Outer case",
                )

        elif tfcoil_variables.i_tf_bucking == 2:
            if tfcoil_variables.i_tf_sup == 1:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t CS \t\t Steel case \t\t WP \t\t Outer case",
                )
            else:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t CS \t\t bucking \t\t conductor \t\t Outer case",
                )

        elif tfcoil_variables.i_tf_bucking == 3:
            if tfcoil_variables.i_tf_sup == 1:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t CS \t\t interface \t\t Steel case \t\t WP \t\t Outer case",
                )
            else:
                po.write(
                    self.outfile,
                    "  Layers \t\t\t\t CS \t\t interface \t\t bucking \t\t conductor \t\t Outer case",
                )

        po.write(
            self.outfile,
            f"  Radial stress \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_r_max, 1e-6)}",
        )
        po.write(
            self.outfile,
            f"  Toroidal stress \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_t_max, 1e-6)}",
        )
        po.write(
            self.outfile,
            f"  Vertical stress \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_z_max, 1e-6)}",
        )
        po.write(
            self.outfile,
            f"  Von-Mises stress \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_vmises_max, 1e-6)}",
        )

        if tfcoil_variables.i_tf_tresca == 1 and tfcoil_variables.i_tf_sup == 1:
            po.write(
                self.outfile,
                f"  Shear (CEA Tresca) \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_tresca_max, 1e-6)}",
            )
        else:
            po.write(
                self.outfile,
                f"  Shear (Tresca) \t\t\t (MPa) \t\t {table_format_arrays(sig_tf_tresca_max, 1e-6)}",
            )

        po.write(self.outfile, "")
        po.write(
            self.outfile,
            f"  Toroidal modulus \t\t\t (GPa) \t\t {table_format_arrays(eyoung_trans, 1e-9)}",
        )
        po.write(
            self.outfile,
            f"  Vertical modulus \t\t\t (GPa) \t\t {table_format_arrays(eyoung_axial, 1e-9)}",
        )
        po.write(self.outfile, "")
        po.ovarre(
            self.outfile,
            "WP transverse modulus (GPa)",
            "(eyoung_wp_trans*1.0d-9)",
            eyoung_wp_trans * 1.0e-9,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "WP vertical modulus (GPa)",
            "(eyoung_wp_axial*1.0d-9)",
            eyoung_wp_axial * 1.0e-9,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "WP transverse Poisson" "s ratio",
            "(poisson_wp_trans)",
            poisson_wp_trans,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "WP vertical-transverse Pois. rat.",
            "(poisson_wp_axial)",
            poisson_wp_axial,
            "OP ",
        )

        # MFILE.DAT data
        for ii in range(n_tf_bucking + 2):
            po.ovarre(
                constants.mfile,
                f"Radial    stress at maximum shear of layer {ii+1} (Pa)",
                f"(sig_tf_r_max({ii+1}))",
                sig_tf_r_max[ii],
            )
            po.ovarre(
                constants.mfile,
                f"toroidal  stress at maximum shear of layer {ii+1} (Pa)",
                f"(sig_tf_t_max({ii+1}))",
                sig_tf_t_max[ii],
            )
            po.ovarre(
                constants.mfile,
                f"Vertical  stress at maximum shear of layer {ii+1} (Pa)",
                f"(sig_tf_z_max({ii+1}))",
                sig_tf_z_max[ii],
            )
            po.ovarre(
                constants.mfile,
                f"Von-Mises stress at maximum shear of layer {ii+1} (Pa)",
                f"(sig_tf_vmises_max({ii+1}))",
                sig_tf_vmises_max[ii],
            )
            if tfcoil_variables.i_tf_tresca == 1 and tfcoil_variables.i_tf_sup == 1:
                po.ovarre(
                    constants.mfile,
                    f"Maximum shear stress for CEA Tresca yield criterion {ii+1} (Pa)",
                    f"(sig_tf_tresca_max({ii+1}))",
                    sig_tf_tresca_max[ii],
                )
            else:
                po.ovarre(
                    constants.mfile,
                    f"Maximum shear stress for the Tresca yield criterion {ii+1} (Pa)",
                    f"(sig_tf_tresca_max({ii+1}))",
                    sig_tf_tresca_max[ii],
                )

        # SIG_TF.json storage
        sig_file_data = {
            "Points per layers": n_radial_array,
            "Radius (m)": radial_array,
            "Radial stress (MPa)": sig_tf_r * 1e-6,
            "Toroidal stress (MPa)": sig_tf_t * 1e-6,
            "Vertical stress (MPa)": sig_tf_z * 1e-6,
            "Radial smear stress (MPa)": sig_tf_smeared_r * 1e-6,
            "Toroidal smear stress (MPa)": sig_tf_smeared_t * 1e-6,
            "Vertical smear stress (MPa)": sig_tf_smeared_z * 1e-6,
            "Von-Mises stress (MPa)": sig_tf_vmises * 1e-6,
            "CEA Tresca stress (MPa)": s_tresca_cond_cea * 1e-6
            if tfcoil_variables.i_tf_sup == 1
            else sig_tf_tresca * 1e-6,
            "rad. displacement (mm)": deflect * 1e3,
        }
        if tfcoil_variables.i_tf_stress_model != 1:
            sig_file_data = {
                **sig_file_data,
                "Radial strain": str_tf_r,
                "Toroidal strain": str_tf_t,
                "Vertical strain": str_tf_z,
            }

        if tfcoil_variables.i_tf_sup == 1:
            sig_file_data = {
                **sig_file_data,
                "WP smeared stress (MPa)": sig_tf_wp_av_z * 1.0e-6,
            }

        sig_file_data = {
            k: list(v) if isinstance(v, numpy.ndarray) else v
            for k, v in sig_file_data.items()
        }

        sig_filename = (
            f2py_compatible_to_string(global_variables.output_prefix) + "SIG_TF.json"
        )
        with open(sig_filename, "w") as f:
            json.dump(sig_file_data, f)

        # TODO sig_tf_wp_av_z is always undefined here. This needs correcting or removing
        # if ( tfcoil_variables.i_tf_sup == 1 ) :
        # write(constants.sig_file,'(t2, "WP"    ," smeared stress", t20, "(MPa)",t26, *(F11.3,3x))') sig_tf_wp_av_z*1.0e-6
        #

        # Quantities from the plane stress stress formulation (no resitive coil output)
        if tfcoil_variables.i_tf_stress_model == 1 and tfcoil_variables.i_tf_sup == 1:
            # Other quantities (displacement strain, etc..)
            po.ovarre(
                self.outfile,
                "Maximum radial deflection at midplane (m)",
                "(deflect)",
                deflect[n_radial_array - 1],
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Vertical strain on casing",
                "(casestr)",
                tfcoil_variables.casestr,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "Radial strain on insulator",
                "(insstrain)",
                tfcoil_variables.insstrain,
                "OP ",
            )

    def tf_averaged_turn_geom(self, jwptf, thwcndut, thicndut, i_tf_sc_mat):
        """
        subroutine straight from Python, see comments in tf_averaged_turn_geom_wrapper
        Authors : J. Morris, CCFE
        Authors : S. Kahn, CCFE
        Setting the TF WP turn geometry for SC magnets from the number
        the current per turn.
        This calculation has two purposes, first to check if a turn can exist
        (positive cable space) and the second to provide its dimensions,
        areas and the (float) number of turns
        """

        # Turn dimension is a an input
        if tfcoil_variables.t_turn_tf_is_input:
            # Turn area [m2]
            a_turn = tfcoil_variables.t_turn_tf**2

            # Current per turn [A]
            tfcoil_variables.cpttf = a_turn * jwptf

        # Turn cable dimension is an input
        elif tfcoil_variables.t_cable_tf_is_input:
            # Turn squared dimension [m]
            tfcoil_variables.t_turn_tf = tfcoil_variables.t_cable_tf + 2.0e0 * (
                thicndut + thwcndut
            )

            # Turn area [m2]
            a_turn = tfcoil_variables.t_turn_tf**2

            # Current per turn [A]
            tfcoil_variables.cpttf = a_turn * jwptf

        # Current per turn is an input
        else:
            # Turn area [m2]
            # Allow for additional inter-layer insulation MDK 13/11/18
            # Area of turn including conduit and inter-layer insulation
            a_turn = tfcoil_variables.cpttf / jwptf

            # Dimension of square cross-section of each turn including inter-turn insulation [m]
            tfcoil_variables.t_turn_tf = numpy.sqrt(a_turn)

        # Square turn assumption
        sctfcoil_module.t_turn_radial = tfcoil_variables.t_turn_tf
        sctfcoil_module.t_turn_toroidal = tfcoil_variables.t_turn_tf

        # See derivation in the following document
        # k:\power plant physics and technology\process\hts\hts coil module for process.docx
        tfcoil_variables.t_conductor = (
            -tfcoil_variables.layer_ins
            + numpy.sqrt(tfcoil_variables.layer_ins**2 + 4.0e00 * a_turn)
        ) / 2 - 2.0e0 * thicndut

        # Total number of turns per TF coil (not required to be an integer)
        n_tf_turn = sctfcoil_module.awptf / a_turn

        # Area of inter-turn insulation: single turn [m2]
        insulation_area = a_turn - tfcoil_variables.t_conductor**2

        # NOTE: Fortran has acstf as an intent(out) variable that was outputting
        # into tfcoil_variables.acstf. The local variable, however, appears to
        # initially hold the value of tfcoil_variables.acstf despite not being
        # intent(in). I have replicated this behaviour in Python for now.
        acstf = copy.copy(tfcoil_variables.acstf)

        # ITER like turn structure
        if i_tf_sc_mat != 6:
            # Radius of rounded corners of cable space inside conduit [m]
            rbcndut = thwcndut * 0.75e0

            # Dimension of square cable space inside conduit [m]
            sctfcoil_module.t_cable = tfcoil_variables.t_conductor - 2.0e0 * thwcndut

            # Cross-sectional area of cable space per turn
            # taking account of rounded inside corners [m2]
            acstf = sctfcoil_module.t_cable**2 - (4.0e0 - numpy.pi) * rbcndut**2

            if acstf <= 0.0e0:
                if tfcoil_variables.t_conductor < 0.0e0:
                    error_handling.fdiags[0] = acstf
                    error_handling.fdiags[1] = sctfcoil_module.t_cable
                    error_handling.report_error(101)
                else:
                    error_handling.fdiags[0] = acstf
                    error_handling.fdiags[1] = sctfcoil_module.t_cable
                    error_handling.report_error(102)
                    rbcndut = 0.0e0
                    acstf = sctfcoil_module.t_cable**2

            # Cross-sectional area of conduit jacket per turn [m2]
            acndttf = tfcoil_variables.t_conductor**2 - acstf

        # REBCO turn structure
        elif i_tf_sc_mat == 6:
            # Diameter of circular cable space inside conduit [m]
            sctfcoil_module.t_cable = tfcoil_variables.t_conductor - 2.0e0 * thwcndut

            # Cross-sectional area of conduit jacket per turn [m2]
            acndttf = tfcoil_variables.t_conductor**2 - acstf

        return acstf, acndttf, insulation_area, n_tf_turn


@numba.njit(cache=True)
def eyoung_t_nested_squares(n, eyoung_j_in, l_in, poisson_j_perp_in):
    """
    Author : C. Swanson, PPPL
    January 2022
    This subroutine gives the smeared transverse elastic
    properties of n members whose cross sectional areas are
    nested squares. It uses the subroutines eyoung_series and
    eyoung_parallel, above, so please be aware of the assumptions
    inherent in those subroutines.

    It assumes that each "leg" of the square cross section
    (vertical slice, as described in Figure 10 of the TF coil
    documentation) is composed of several layers under stress in
    series, and each leg is under stress in parallel with every
    other leg.
    """
    eyoung_j_working = numpy.zeros((n,))
    l_working = numpy.zeros((n,))
    poisson_j_perp_working = numpy.zeros((n,))

    # First member
    eyoung_j_working[0] = eyoung_j_in[0]
    l_working[0] = l_in[0]
    poisson_j_perp_working[0] = poisson_j_perp_in[0]

    for ii in range(1, n):
        # Initialize the leg of which this is the new member
        eyoung_j_working[ii] = eyoung_j_in[ii]
        l_working[ii] = l_working[ii - 1] + l_in[ii]
        poisson_j_perp_working[ii] = poisson_j_perp_in[ii]

        # Serial-composite the new layer of this member into the previous legs
        # changed from range(ii-1) because range(0) == []
        for jj in range(ii):
            (
                eyoung_j_working[jj],
                l_working[jj],
                poisson_j_perp_working[jj],
            ) = eyoung_series(
                eyoung_j_working[ii],
                l_in[ii],
                poisson_j_perp_working[ii],
                eyoung_j_working[jj],
                l_working[jj],
                poisson_j_perp_working[jj],
            )

    # Find stiffest leg
    eyoung_stiffest = max(eyoung_j_working)

    eyoung_j_out = 0
    l_out = 0
    poisson_j_perp_out = 0

    # Parallel-composite them all together
    for ii in range(n):
        eyoung_j_out, l_out, poisson_j_perp_out = eyoung_parallel(
            eyoung_j_working[ii],
            l_in[ii],
            poisson_j_perp_working[ii],
            eyoung_j_out,
            l_out,
            poisson_j_perp_out,
        )

    return eyoung_j_out, l_out, poisson_j_perp_out, eyoung_stiffest


@numba.njit(cache=True)
def eyoung_series(eyoung_j_1, l_1, poisson_j_perp_1, eyoung_j_2, l_2, poisson_j_perp_2):
    """
    Author : C. Swanson, PPPL
    January 2022
    See Issue #1205 for derivation PDF
    This subroutine gives the smeared elastic properties of two
    members that are carrying a force in series with each other.
    The force goes in direction j.
    The assumption is that the stresses in j are equal.
    The smeared Young's modulus is the inverse of the average of
    the inverse of the Young's moduli, weighted by the length
    of the members in j.
    Members 1 and 2 are the individual members to be smeared.
    Member 3 is the effective smeared member (output).
    The smeared Poisson's ratio is the averaged of the Poisson's
    ratios, weighted by the quantity (Young's modulus / length of
    the members in j).

    If you're dealing with anisotropy, please pay attention to the
    fact that the specific Young's Modulus used here is that in
    the j direction, and the specific Poisson's ratio used here is
    that between the j and transverse directions in that order.
    (transverse strain / j strain, under j stress)
    The smeared Poisson's ratio is computed assuming the transverse
    dynamics are isotropic, and that the two members are free to
    shrink/expand under Poisson effects without interference from
    each other. This may not be true of your case.

    To build up a composite smeared member of any number of
    individual members, you can pass the same properties for
    members 2 and 3, and call it successively, using the properties
    of each member as the first triplet of arguments. This way, the
    last triplet acts as a "working sum":
    call eyoung_series(triplet1, triplet2, tripletOUT)
    call eyoung_series(triplet3, tripletOUT, tripletOUT)
    call eyoung_series(triplet4, tripletOUT, tripletOUT)
    ... etc.
    So that tripletOUT would eventually have the smeared properties
    of the total composite member.
    """

    if eyoung_j_1 * eyoung_j_2 == 0:
        # poisson_j_perp_3 = 0
        if eyoung_j_1 == 0:
            poisson_j_perp_3 = poisson_j_perp_1
        else:
            poisson_j_perp_3 = poisson_j_perp_2

        eyoung_j_3 = 0.0
        l_3 = l_1 + l_2
    else:
        poisson_j_perp_3 = (
            poisson_j_perp_1 * l_1 / eyoung_j_1 + poisson_j_perp_2 * l_2 / eyoung_j_2
        ) / (l_1 / eyoung_j_1 + l_2 / eyoung_j_2)
        eyoung_j_3 = (l_1 + l_2) / (l_1 / eyoung_j_1 + l_2 / eyoung_j_2)
        l_3 = l_1 + l_2

    eyoung_j_3 = numpy.array(eyoung_j_3)
    poisson_j_perp_3 = numpy.array(poisson_j_perp_3)

    return eyoung_j_3, l_3, poisson_j_perp_3


@numba.njit(cache=True)
def eyoung_parallel(
    eyoung_j_1, a_1, poisson_j_perp_1, eyoung_j_2, a_2, poisson_j_perp_2
):
    """
    Author : C. Swanson, PPPL
    January 2022
    See Issue #1205 for derivation PDF
    This subroutine gives the smeared elastic properties of two
    members that are carrying a force in parallel with each other.
    The force goes in direction j.
    Members 1 and 2 are the individual members to be smeared.
    Member 3 is the effective smeared member (output).
    This is pretty easy because the smeared properties are simply
    the average weighted by the cross-sectional areas perpendicular
    to j.
    The assumption is that the strains in j are equal.
    If you're dealing with anisotropy, please pay attention to the
    fact that the specific Young's Modulus used here is that in
    the j direction, and the specific Poisson's ratio used here is
    that between the j and transverse directions in that order.
    (transverse strain / j strain, under j stress)
    The smeared Poisson's ratio is computed assuming the transverse
    dynamics are isotropic, and that the two members are free to
    shrink/expand under Poisson effects without interference from
    each other. This may not be true of your case.

    To build up a composite smeared member of any number of
    individual members, you can pass the same properties for
    members 2 and 3, and call it successively, using the properties
    of each member as the first triplet of arguments. This way, the
    last triplet acts as a "working sum":
    call eyoung_parallel(triplet1, triplet2, tripletOUT)
    call eyoung_parallel(triplet3, tripletOUT, tripletOUT)
    call eyoung_parallel(triplet4, tripletOUT, tripletOUT)
    ... etc.
    So that tripletOUT would eventually have the smeared properties
    of the total composite member.
    """
    poisson_j_perp_3 = (poisson_j_perp_1 * a_1 + poisson_j_perp_2 * a_2) / (a_1 + a_2)
    eyoung_j_3 = (eyoung_j_1 * a_1 + eyoung_j_2 * a_2) / (a_1 + a_2)
    a_3 = a_1 + a_2

    return eyoung_j_3, a_3, poisson_j_perp_3


@numba.njit(cache=True)
def sigvm(sx: float, sy: float, sz: float, txy: float, txz: float, tyz: float) -> float:
    """Calculates Von Mises stress in a TF coil
    author: P J Knight, CCFE, Culham Science Centre
    author: B Reimer, FEDC
    This routine calculates the Von Mises combination of
    stresses (Pa) in a TF coil.
    AEA FUS 251: A User's Guide to the PROCESS Systems Code

    :param sx: In-plane stress in X direction [Pa]
    :param sy: In-plane stress in Y direction [Pa]
    :param sz: In-plane stress in Z direction [Pa]
    :param txy: Out-of-plane stress in X-Y plane [Pa]
    :param txz: Out-of-plane stress in X-Z plane [Pa]
    :param tyz: Out-of-plane stress in Y-Z plane [Pa]

    :returns: Von Mises combination of stresses (Pa) in a TF coil.
    """

    return numpy.sqrt(
        0.5
        * (
            (sx - sy) ** 2
            + (sx - sz) ** 2
            + (sz - sy) ** 2
            + 6 * (txy**2 + txz**2 + tyz**2)
        )
    )


@numba.njit(cache=True, error_model="numpy")
def extended_plane_strain(
    nu_t,
    nu_zt,
    ey_t,
    ey_z,
    rad,
    d_curr,
    v_force,
    nlayers,
    n_radial_array,
    i_tf_bucking,
):
    """Author : C. Swanson, PPPL and S. Kahn, CCFE
    September 2021
    There is a writeup of the derivation of this model on the gitlab server.
    https://git.ccfe.ac.uk/process/process/-/issues/1414
    This surboutine estimates the radial displacement, stresses, and strains of
    the inboard midplane of the TF. It assumes that structures are axisymmetric
    and long in the axial (z) direction, the "axisymmetric extended plane strain"
    problem. The TF is assumed to be constructed from some number of layers,
    within which materials properties and current densities are uniform.
    The 1D radially-resolved solution is reduced to a 0D matrix inversion problem
    using analytic solutions to Lame's thick cylinder problem. Materials may be
    transverse-isotropic in Young's modulus and Poisson's ratio. The consraints
    are: Either zero radial stress or zero radial displacement at the inner
    surface (depending on whether the inner radius is zero), zero radial stress
    at the outer surface, total axial force (tension) is equal to the input,
    and optionally the axial force of an inner subset of cylinders is zero (slip
    conditions between the inner and outer subset). The matrix inversion / linear
    solve is always 4x4, no matter how many layers there are.
    The problem is formulated around a solution vector (A,B,eps_z,1.0,eps_z_slip)
    where A and B are the parameters in Lame's solution where u = A*r + B/r, u
    is the radial displacement. eps_z is the axial strain on the outer, force-
    carrying layers. eps_z_slip is the axial strain on the inner, non-force-
    carrying layers (optionally). The solution vector is A,B at the outermost
    radius, and is transformed via matrix multiplication into those A,B
    values at other radii. The constraints are inner products with this vector,
    and so when stacked form a matrix to invert.
    """
    # outputs
    sigr = numpy.zeros((n_radial_array * nlayers,))
    # Stress distribution in the radial direction (r) [Pa]

    sigt = numpy.zeros((n_radial_array * nlayers,))
    # Stress distribution in the toroidal direction (t) [Pa]

    sigz = numpy.zeros((n_radial_array * nlayers,))
    # Stress distribution in the vertical direction (z)

    str_r = numpy.zeros((n_radial_array * nlayers,))
    # Strain distribution in the radial direction (r)

    str_t = numpy.zeros((n_radial_array * nlayers,))
    # Strain distribution in the toroidal direction (t)

    str_z = numpy.zeros((n_radial_array * nlayers,))
    # Uniform strain in the vertical direction (z)

    r_deflect = numpy.zeros((n_radial_array * nlayers,))
    # Radial displacement radial distribution [m]

    rradius = numpy.zeros((n_radial_array * nlayers,))
    # Radius array [m]

    # local arrays
    # Stiffness form of compliance tensor
    nu_tz = numpy.zeros((nlayers,))
    # Transverse-axial Poisson's ratio
    # (ratio of axial strain to transverse strain upon transverse stress)
    ey_bar_z = numpy.zeros((nlayers,))
    # Axial effective Young's modulus (zero cross-strains, not stresses) [Pa]
    ey_bar_t = numpy.zeros((nlayers,))
    # Transverse effective Young's modulus [Pa]
    nu_bar_t = numpy.zeros((nlayers,))
    # Transverse effective Poisson's ratio
    nu_bar_tz = numpy.zeros((nlayers,))
    # Transverse-axial effective Poisson's ratio
    nu_bar_zt = numpy.zeros((nlayers,))
    # Axial-transverse effective Poisson's ratio

    # Lorentz force parameters
    currents = numpy.zeros((nlayers,))
    # Currents in each layer [A]
    currents_enclosed = numpy.zeros((nlayers,))
    # Currents enclosed by inner radius of each layer [A]
    f_lin_fac = numpy.zeros((nlayers,))
    # Factor that multiplies r linearly in the force density
    f_rec_fac = numpy.zeros((nlayers,))
    # Factor that multiplies r reciprocally in the force density
    f_int_A = numpy.zeros((nlayers,))
    # Force density integral that adds to Lame parameter A
    f_int_B = numpy.zeros((nlayers,))
    # Force density integral that adds to Lame parameter B

    # Layer transfer matrices
    M_int = numpy.zeros(
        (
            5,
            5,
            nlayers,
        ),
    )
    # Matrix that transforms the Lame parmeter vector from the
    # outer radius to the inner radius of each layer
    M_ext = numpy.zeros(
        (
            5,
            5,
            nlayers,
        ),
    )
    # Matrix that transforms the Lame parmeter vector from the
    # inner radius of one layer to the outer radius of the
    # next inner.
    M_tot = numpy.zeros(
        (
            5,
            5,
            nlayers,
        ),
    )
    # Matrix that transforms the Lame parmeter vector from the
    # outer radius of the outer layer to the inner radius of
    # each.

    # Axial force inner product
    v_force_row = numpy.zeros(
        (
            1,
            5,
        ),
    )
    # Row vector (matrix multiplication is inner product) to
    # obtain the axial force from the force-carrying layers
    v_force_row_slip = numpy.zeros(
        (
            1,
            5,
        ),
    )
    # Row vector (matrix multiplication is inner product) to
    # obtain the axial force inner slip layers (no net force)
    rad_row_helper = numpy.zeros(
        (
            1,
            5,
        ),
    )
    # A helper variable to store [radius, 1, 0, 0, 0] in row

    # Boundary condition matrix
    M_bc = numpy.zeros(
        (
            4,
            5,
        ),
    )
    # Boundary condition matrix. Multiply this with the
    # outermost solution vector, (A,B,eps_z,1.0,eps_z_slip),
    # to obtain a zero vector.
    M_toinv = numpy.zeros(
        (
            4,
            4,
        ),
    )
    # Matrix to invert to get the solution
    RHS_vec = numpy.zeros((4,))
    # Right-hand-side vector to divide M_toinv
    A_vec_solution = numpy.zeros((5,))
    # Solution vector, Lame parameters at outer radius, strain
    # of force-carrying layers, and strain of slip layers
    # (A,B,eps_z,1,eps_z_slip)

    # Constructing the solution everywhere
    A_vec_layer = numpy.zeros((5,))
    # Lame parameters and strains vector at outer radius
    # of each layer

    # The stress calcualtion differential equations is analytically sloved
    # The final solution is given by the layer boundary conditions on
    # radial stress and displacement between layers solved
    # The problem is set as aa.cc = bb, cc being the constant we search

    # Inner slip layers parameters
    # Section 15 in the writeup
    # Innermost layer that takes axial force. Layers inner of this
    # have zero net axial force, to include CS decoupling.
    # This configuration JUST HAPPENS to work out because of
    # the specific i_tf_bucking options; if those are changed,
    # will need a switch here.
    nonslip_layer = i_tf_bucking

    if nonslip_layer < 1:
        nonslip_layer = 1

    # Stiffness tensor factors
    # Section 3 in the writeup
    # With Section 12 anisotropic materials properties
    # ***
    # Dependent Poisson's ratio: nu-transverse-axial
    # from nu-axial-transverse and the Young's moduli
    nu_tz[:] = nu_zt * ey_t / ey_z

    # Effective Young's Moduli and Poisson's ratios
    # holding strain, not stress, cross-terms constant
    ey_bar_z[:] = ey_z * (1 - nu_t) / (1 - nu_t - 2 * nu_tz * nu_zt)
    ey_bar_t[:] = (
        ey_t * (1 - nu_tz * nu_zt) / (1 - nu_t - 2 * nu_tz * nu_zt) / (1 + nu_t)
    )

    nu_bar_t[:] = (nu_t + nu_tz * nu_zt) / (1 - nu_tz * nu_zt)
    nu_bar_tz[:] = nu_tz / (1 - nu_t)
    nu_bar_zt[:] = nu_zt * (1 + nu_t) / (1 - nu_tz * nu_zt)

    # Lorentz force parameters
    # Section 13 in the writeup
    # ***
    # Currents in each layer [A]
    currents[:] = numpy.pi * d_curr * (rad[1 : nlayers + 1] ** 2 - rad[:nlayers] ** 2)
    # Currents enclosed by inner radius of each layer [A]
    currents_enclosed[0] = 0.0e0

    for ii in range(1, nlayers):
        currents_enclosed[ii] = currents_enclosed[ii - 1] + currents[ii - 1]
    # Factor that multiplies r linearly in the force density
    f_lin_fac[:] = RMU0 / 2.0e0 * d_curr**2
    # Factor that multiplies r reciprocally in the force density
    f_rec_fac[:] = (
        RMU0
        / 2.0e0
        * (d_curr * currents_enclosed / numpy.pi - d_curr**2 * rad[:nlayers] ** 2)
    )
    # Force density integral that adds to Lame parameter A
    f_int_A[:] = 0.5e0 * f_lin_fac * (
        rad[1 : nlayers + 1] ** 2 - rad[:nlayers] ** 2
    ) + f_rec_fac * numpy.log(rad[1 : nlayers + 1] / rad[:nlayers])
    if f_rec_fac[0] == 0e0:
        f_int_A[0] = 0.5e0 * f_lin_fac[0] * (rad[1] ** 2 - rad[0] ** 2)

    # Force density integral that adds to Lame parameter B
    f_int_B[:] = 0.25e0 * f_lin_fac * (
        rad[1 : nlayers + 1] ** 4 - rad[:nlayers] ** 4
    ) + 0.5e0 * f_rec_fac * (rad[1 : nlayers + 1] ** 2 - rad[:nlayers] ** 2)

    # Transformation matrix from outer to inner Lame parameters
    # Section 5 in the writeup
    # With Section 12 anisotropic materials properties
    # ***
    # M_int[kk] multiplies Lame parameter vector of layer kk (A,B,eps_z,1.0,eps_z_slip)
    # and transforms the values at the outer radius to the values at the inner radius
    for kk in range(nlayers):
        M_int[0, 0, kk] = 1.0e0
        M_int[1, 1, kk] = 1.0e0
        M_int[2, 2, kk] = 1.0e0
        M_int[3, 3, kk] = 1.0e0
        M_int[4, 4, kk] = 1.0e0

        M_int[0, 3, kk] = -0.5e0 / ey_bar_t[kk] * f_int_A[kk]
        M_int[1, 3, kk] = 0.5e0 / ey_bar_t[kk] * f_int_B[kk]

    # Transformation matrix between layers
    # Section 6 in the writeup
    # With Section 12 anisotropic materials properties
    # With Section 15 inner slip-decoupled layers
    # ***
    # M_ext[kk] multiplies Lame parameter vector of layer kk (A,B,eps_z,1.0,eps_z_slip)
    # and transforms the values at the inner radius to the values at the outer radius
    # of layer kk-1
    for kk in range(1, nonslip_layer - 1):
        ey_fac = ey_bar_t[kk] / ey_bar_t[kk - 1]
        M_ext[0, 2, kk] = 0.0e0
        M_ext[0, 4, kk] = 0.5e0 * (ey_fac * nu_bar_zt[kk] - nu_bar_zt[kk - 1])

    if nonslip_layer > 1:
        ey_fac = ey_bar_t[nonslip_layer - 1] / ey_bar_t[nonslip_layer - 2]
        M_ext[0, 2, nonslip_layer - 1] = 0.5e0 * ey_fac * nu_bar_zt[nonslip_layer - 1]
        M_ext[0, 4, nonslip_layer - 1] = 0.5e0 * (-nu_bar_zt[nonslip_layer - 2])

    for kk in range(nonslip_layer, nlayers):
        ey_fac = ey_bar_t[kk] / ey_bar_t[kk - 1]
        M_ext[0, 2, kk] = 0.5e0 * (ey_fac * nu_bar_zt[kk] - nu_bar_zt[kk - 1])
        M_ext[0, 4, kk] = 0.0e0

    for kk in range(1, nlayers):
        ey_fac = ey_bar_t[kk] / ey_bar_t[kk - 1]
        M_ext[0, 0, kk] = 0.5e0 * (ey_fac * (1 + nu_bar_t[kk]) + 1 - nu_bar_t[kk - 1])
        if rad[kk] > 0e0:
            M_ext[0, 1, kk] = (
                0.5e0
                / rad[kk] ** 2
                * (1 - nu_bar_t[kk - 1] - ey_fac * (1 - nu_bar_t[kk]))
            )

        M_ext[1, 0, kk] = rad[kk] ** 2 * (1 - M_ext[0, 0, kk])
        M_ext[1, 1, kk] = 1 - rad[kk] ** 2 * M_ext[0, 1, kk]
        M_ext[1, 2, kk] = -rad[kk] ** 2 * M_ext[0, 2, kk]
        M_ext[1, 4, kk] = -rad[kk] ** 2 * M_ext[0, 4, kk]
        M_ext[2, 2, kk] = 1.0e0
        M_ext[3, 3, kk] = 1.0e0
        M_ext[4, 4, kk] = 1.0e0

    # Total transformation matrix, from Lame parmeters at outside to
    # Lame parameters at inside of each layer
    # Section 7 in the writeup
    # ***
    M_tot[:, :, nlayers - 1] = M_int[:, :, nlayers - 1]

    for kk in range(nlayers - 2, -1, -1):
        M_tot[:, :, kk] = numpy.ascontiguousarray(M_int[:, :, kk]) @ (
            numpy.ascontiguousarray(M_ext[:, :, kk + 1])
            @ numpy.ascontiguousarray(M_tot[:, :, kk + 1])
        )

    # Axial force inner product. Dot-product this with the
    # outermost solution vector, (A,B,eps_z,1.0,eps_z_slip),
    # to obtain the axial force.
    # Section 8 in the writeup
    # ***
    # Axial stiffness products
    ey_bar_z_area = numpy.pi * sum(
        ey_bar_z[nonslip_layer - 1 : nlayers]
        * (
            rad[nonslip_layer : nlayers + 1] ** 2
            - rad[nonslip_layer - 1 : nlayers] ** 2
        )
    )
    ey_bar_z_area_slip = numpy.pi * sum(
        ey_bar_z[: nonslip_layer - 1]
        * (rad[1:nonslip_layer] ** 2 - rad[: nonslip_layer - 1] ** 2)
    )

    # Axial stiffness inner product, for layers which carry axial force
    rad_row_helper[0, :] = [rad[nlayers] ** 2, 1e0, 0e0, 0e0, 0e0]
    v_force_row[:, :] = (
        2e0 * numpy.pi * ey_bar_z[nlayers - 1] * nu_bar_tz[nlayers - 1] * rad_row_helper
    )
    rad_row_helper[0, :] = [rad[nonslip_layer - 1] ** 2, 1e0, 0e0, 0e0, 0e0]
    v_force_row[:, :] = v_force_row - 2e0 * numpy.pi * ey_bar_z[
        nonslip_layer - 1
    ] * nu_bar_tz[nonslip_layer - 1] * (
        rad_row_helper @ numpy.ascontiguousarray(M_tot[:, :, nonslip_layer - 1])
    )
    for kk in range(nonslip_layer, nlayers):
        rad_row_helper[0, :] = [rad[kk] ** 2, 1e0, 0e0, 0e0, 0e0]
        v_force_row[:, :] = v_force_row + 2e0 * numpy.pi * (
            ey_bar_z[kk - 1] * nu_bar_tz[kk - 1] - ey_bar_z[kk] * nu_bar_tz[kk]
        ) * (rad_row_helper @ numpy.ascontiguousarray(M_tot[:, :, kk]))

    # Include the effect of axial stiffness
    v_force_row[0, 2] += ey_bar_z_area

    # Axial stiffness inner product, for layers which DON'T carry force
    if nonslip_layer > 1:
        rad_row_helper[0, :] = [rad[nonslip_layer - 1] ** 2, 1e0, 0e0, 0e0, 0e0]
        v_force_row_slip[:, :] = (
            2e0
            * numpy.pi
            * ey_bar_z[nonslip_layer - 2]
            * nu_bar_tz[nonslip_layer - 2]
            * (rad_row_helper @ numpy.ascontiguousarray(M_tot[:, :, nonslip_layer - 1]))
        )
        rad_row_helper[0, :] = [rad[0] ** 2, 1e0, 0e0, 0e0, 0e0]
        v_force_row_slip[:, :] -= (
            2e0
            * numpy.pi
            * ey_bar_z[0]
            * nu_bar_tz[0]
            * (rad_row_helper @ numpy.ascontiguousarray(M_tot[:, :, 0]))
        )
        for kk in range(1, nonslip_layer - 1):
            rad_row_helper[0, :] = [rad[kk] ** 2, 1e0, 0e0, 0e0, 0e0]
            v_force_row_slip[:, :] += (
                2e0
                * numpy.pi
                * (ey_bar_z[kk - 1] * nu_bar_tz[kk - 1] - ey_bar_z[kk] * nu_bar_tz[kk])
                * (rad_row_helper @ numpy.ascontiguousarray(M_tot[:, :, kk]))
            )
        # Include the effect of axial stiffness
        v_force_row_slip[0, 4] += ey_bar_z_area_slip
    else:
        # If there's no inner slip layer, still need a finite 5th
        # element to ensure no singular matrix
        v_force_row_slip[0, :] = [0e0, 0e0, 0e0, 0e0, 1e0]

    # Boundary condition matrix. Multiply this with the
    # outermost solution vector, (A,B,eps_z,1.0,eps_z_slip),
    # to obtain a zero vector.
    # Solved to get the Lame parameters.
    # Section 9 in the writeup
    # ***
    # Outer boundary condition row, zero radial stress
    M_bc[0, :] = [
        (1e0 + nu_bar_t[nlayers - 1]) * rad[nlayers] ** 2,
        -1e0 + nu_bar_t[nlayers - 1],
        nu_bar_zt[nlayers - 1] * rad[nlayers] ** 2,
        0e0,
        0e0,
    ]
    # Inner boundary condition row, zero radial stress
    # or zero displacement if rad(1)=0
    if nonslip_layer > 1:
        M_bc[1, :] = [
            (1e0 + nu_bar_t[0]) * rad[0] ** 2,
            -1e0 + nu_bar_t[0],
            0e0,
            0e0,
            nu_bar_zt[0] * rad[0] ** 2,
        ]
    else:
        M_bc[1, :] = [
            (1e0 + nu_bar_t[0]) * rad[0] ** 2,
            -1e0 + nu_bar_t[0],
            nu_bar_zt[0] * rad[0] ** 2,
            0e0,
            0e0,
        ]

    M_bc[1, :] = numpy.ascontiguousarray(M_bc[1, :]) @ numpy.ascontiguousarray(
        M_tot[:, :, 0]
    )
    # Axial force boundary condition
    M_bc[2, :] = v_force_row[0, :]
    M_bc[2, 3] = M_bc[2, 3] - v_force
    # Axial force boundary condition of slip layers
    M_bc[3, :] = v_force_row_slip[0, :]

    # The solution, the outermost Lame parameters A,B
    # and the axial strains of the force-carrying and
    # slip layers eps_z and eps_z_slip.
    # Section 10 in the writeup
    # ***
    M_toinv[:, :3] = M_bc[
        :,
        :3,
    ]
    M_toinv[:, 3] = M_bc[:, 4]
    RHS_vec[:] = -M_bc[:, 3]

    A_vec_solution[:4] = numpy.linalg.solve(M_toinv, RHS_vec)

    # maths_library.linesolv(M_toinv, RHS_vec, A_vec_solution[:4])
    A_vec_solution[4] = A_vec_solution[3]
    A_vec_solution[3] = 1

    # Radial/toroidal/vertical stress radial distribution
    # ------
    # Radial displacement, stress and strain distributions

    A_vec_layer[:] = A_vec_solution[:]
    for ii in range(nlayers - 1, -1, -1):
        A_layer = A_vec_layer[0]
        B_layer = A_vec_layer[1]

        dradius = (rad[ii + 1] - rad[ii]) / (n_radial_array - 1)

        for jj in range(ii * n_radial_array, (ii + 1) * n_radial_array):
            rradius[jj] = rad[ii] + dradius * (jj - (n_radial_array * ii))

            f_int_A_plot = 0.5e0 * f_lin_fac[ii] * (
                rad[ii + 1] ** 2 - rradius[jj] ** 2
            ) + f_rec_fac[ii] * numpy.log(rad[ii + 1] / (rradius[jj]))
            f_int_B_plot = 0.25e0 * f_lin_fac[ii] * (
                rad[ii + 1] ** 4 - rradius[jj] ** 4
            ) + 0.5e0 * f_rec_fac[ii] * (rad[ii + 1] ** 2 - rradius[jj] ** 2)
            A_plot = A_layer - 0.5e0 / ey_bar_t[ii] * f_int_A_plot
            B_plot = B_layer + 0.5e0 / ey_bar_t[ii] * f_int_B_plot

            # Radial displacement
            r_deflect[jj] = A_plot * rradius[jj] + B_plot / rradius[jj]

            # Radial strain
            str_r[jj] = A_plot - B_plot / rradius[jj] ** 2
            # Azimuthal strain
            str_t[jj] = A_plot + B_plot / rradius[jj] ** 2
            # Axial strain
            if ii < nonslip_layer - 1:
                str_z[jj] = A_vec_solution[4]
            else:
                str_z[jj] = A_vec_solution[2]

            # Radial stress
            sigr[jj] = ey_bar_t[ii] * (
                str_r[jj] + (nu_bar_t[ii] * str_t[jj]) + (nu_bar_zt[ii] * str_z[jj])
            )
            # Aximuthal stress
            sigt[jj] = ey_bar_t[ii] * (
                str_t[jj] + (nu_bar_t[ii] * str_r[jj]) + (nu_bar_zt[ii] * str_z[jj])
            )
            # Axial stress
            sigz[jj] = ey_bar_z[ii] * (
                str_z[jj] + (nu_bar_tz[ii] * (str_r[jj] + str_t[jj]))
            )

        A_vec_layer = numpy.ascontiguousarray(M_tot[:, :, ii]) @ A_vec_solution
        A_vec_layer = numpy.ascontiguousarray(M_ext[:, :, ii]) @ A_vec_layer
    # ------

    return rradius, sigr, sigt, sigz, str_r, str_t, str_z, r_deflect


@numba.njit(cache=True)
def plane_stress(nu, rad, ey, j, nlayers, n_radial_array):
    """Calculates the stresses in a superconductor TF coil
    inboard leg at the midplane using the plain stress approximation
    author: P J Knight, CCFE, Culham Science Centre
    author: J Morris, CCFE, Culham Science Centre
    author: S Kahn, CCFE, Culham Science Centre
    This routine calculates the stresses in a superconductor TF coil
    inboard leg at midplane.
    <P>A 2 layer plane stress model developed by CCFE is used. The first layer
    is the steel case inboard of the winding pack, and the second
    layer is the winding pack itself.
    PROCESS Superconducting TF Coil Model, J. Morris, CCFE, 1st May 2014
    """
    alpha = numpy.zeros((nlayers,))
    beta = numpy.zeros((nlayers,))
    # Lorentz body force parametres

    area = numpy.zeros((nlayers,))
    # Layer area

    aa = numpy.zeros(
        (
            2 * nlayers,
            2 * nlayers,
        )
    )
    # Matrix encoding the integration constant cc coeficients

    bb = numpy.zeros((2 * nlayers,))
    # Vector encoding the alpha/beta (lorentz forces) contribution

    cc = numpy.zeros((2 * nlayers,))
    c1 = numpy.zeros((nlayers,))
    c2 = numpy.zeros((nlayers,))
    # Integration constants vector (solution)

    rradius = numpy.zeros((nlayers * n_radial_array,))
    # Radius array [m]

    sigr = numpy.zeros((nlayers * n_radial_array,))
    # Radial stress radial distribution [Pa]

    sigt = numpy.zeros((nlayers * n_radial_array,))
    # Toroidal stress radial distribution [Pa]

    r_deflect = numpy.zeros((nlayers * n_radial_array,))
    # Radial deflection (displacement) radial distribution [m]

    kk = ey / (1 - nu**2)

    # Lorentz forces parametrisation coeficients (array equation)
    alpha = 0.5e0 * RMU0 * j**2 / kk

    inner_layer_curr = 0.0e0
    for ii in range(nlayers):
        beta[ii] = (
            0.5e0
            * RMU0
            * j[ii]
            * (inner_layer_curr - numpy.pi * j[ii] * rad[ii] ** 2)
            / (numpy.pi * kk[ii])
        )

        # Layer area
        area[ii] = numpy.pi * (rad[ii + 1] ** 2 - rad[ii] ** 2)

        # Total current carried by the inners layers
        inner_layer_curr = inner_layer_curr + area[ii] * j[ii]
    # ***

    # Null radial stress at R(1)
    aa[0, 0] = kk[0] * (1.0e0 + nu[0])
    aa[0, 1] = -kk[0] * (1.0e0 - nu[0]) / (rad[0] ** 2)

    # Inter-layer boundary conditions
    if nlayers != 1:
        for ii in range(nlayers - 1):
            # Continuous radial normal stress at R(ii+1)
            aa[2 * ii + 1, 2 * ii] = kk[ii] * (1.0e0 + nu[ii])
            aa[2 * ii + 1, 2 * ii + 1] = -kk[ii] * (1.0e0 - nu[ii]) / rad[ii + 1] ** 2
            aa[2 * ii + 1, 2 * ii + 2] = -kk[ii + 1] * (1.0e0 + nu[ii + 1])
            aa[2 * ii + 1, 2 * ii + 3] = (
                kk[ii + 1] * (1.0e0 - nu[ii + 1]) / rad[ii + 1] ** 2
            )

            # Continuous displacement at R(ii+1)
            aa[2 * ii + 2, 2 * ii] = rad[ii + 1]
            aa[2 * ii + 2, 2 * ii + 1] = 1.0e0 / rad[ii + 1]
            aa[2 * ii + 2, 2 * ii + 2] = -rad[ii + 1]
            aa[2 * ii + 2, 2 * ii + 3] = -1.0e0 / rad[ii + 1]

    # Radial stress = 0
    aa[2 * (nlayers - 1) + 1, 2 * (nlayers - 1)] = kk[nlayers - 1] * (
        1.0e0 + nu[nlayers - 1]
    )
    aa[2 * (nlayers - 1) + 1, 2 * (nlayers - 1) + 1] = (
        -kk[nlayers - 1] * (1.0e0 - nu[nlayers - 1]) / rad[nlayers] ** 2
    )
    # ***

    # Right hand side vector bb
    # ***
    # Null radial stress at R(1)
    bb[0] = -kk[0] * (
        0.125e0 * alpha[0] * (3.0e0 + nu[0]) * rad[0] ** 2
        + 0.5e0 * beta[0] * (1.0e0 + (1.0e0 + nu[0]) * numpy.log(rad[0]))
    )

    # Inter-layer boundary conditions
    if nlayers != 1:
        for ii in range(nlayers - 1):
            # Continuous radial normal stress at R[ii+1]
            bb[2 * ii + 1] = -kk[ii] * (
                0.125e0 * alpha[ii] * (3.0e0 + nu[ii]) * rad[ii + 1] ** 2
                + 0.5e0 * beta[ii] * (1.0e0 + (1.0e0 + nu[ii]) * numpy.log(rad[ii + 1]))
            ) + kk[ii + 1] * (
                0.125e0 * alpha[ii + 1] * (3.0e0 + nu[ii + 1]) * rad[ii + 1] ** 2
                + 0.5e0
                * beta[ii + 1]
                * (1.0e0 + (1.0e0 + nu[ii + 1]) * numpy.log(rad[ii + 1]))
            )

            # Continuous displacement at R(ii+1)
            bb[2 * ii + 2] = (
                -0.125e0 * alpha[ii] * rad[ii + 1] ** 3
                - 0.5e0 * beta[ii] * rad[ii + 1] * numpy.log(rad[ii + 1])
                + 0.125e0 * alpha[ii + 1] * rad[ii + 1] ** 3
                + 0.5e0 * beta[ii + 1] * rad[ii + 1] * numpy.log(rad[ii + 1])
            )

    # Null radial stress at R(nlayers+1)
    bb[2 * (nlayers - 1) + 1] = -kk[nlayers - 1] * (
        0.125e0 * alpha[nlayers - 1] * (3.0e0 + nu[nlayers - 1]) * rad[nlayers] ** 2
        + 0.5e0
        * beta[nlayers - 1]
        * (1.0e0 + (1.0e0 + nu[nlayers - 1]) * numpy.log(rad[nlayers]))
    )
    # ***

    #  Find solution vector cc
    # ***
    aa = numpy.asfortranarray(aa)
    cc = numpy.linalg.solve(aa, bb)

    # maths_library.linesolv(aa, bb, cc)

    #  Multiply c by (-1) (John Last, internal CCFE memorandum, 21/05/2013)
    for ii in range(nlayers):
        c1[ii] = cc[2 * ii]
        c2[ii] = cc[2 * ii + 1]
    # ***
    # ------

    # Radial/toroidal/vertical stress radial distribution
    # ------

    for ii in range(nlayers):
        dradius = (rad[ii + 1] - rad[ii]) / n_radial_array
        for jj in range(ii * n_radial_array, (ii + 1) * n_radial_array):
            rad_c = rad[ii] + dradius * (jj - n_radial_array * ii)
            rradius[jj] = rad_c

            # Radial stress radial distribution [Pa]
            sigr[jj] = kk[ii] * (
                (1.0e0 + nu[ii]) * c1[ii]
                - ((1.0e0 - nu[ii]) * c2[ii]) / rad_c**2
                + 0.125e0 * (3.0e0 + nu[ii]) * alpha[ii] * rad_c**2
                + 0.5e0 * beta[ii] * (1.0e0 + (1.0e0 + nu[ii]) * numpy.log(rad_c))
            )

            # Radial stress radial distribution [Pa]
            sigt[jj] = kk[ii] * (
                (1.0e0 + nu[ii]) * c1[ii]
                + (1.0e0 - nu[ii]) * c2[ii] / rad_c**2
                + 0.125e0 * (1.0e0 + 3.0e0 * nu[ii]) * alpha[ii] * rad_c**2
                + 0.5e0 * beta[ii] * (nu[ii] + (1.0e0 + nu[ii]) * numpy.log(rad_c))
            )

            #  Deflection [m]
            r_deflect[jj] = (
                c1[ii] * rad_c
                + c2[ii] / rad_c
                + 0.125e0 * alpha[ii] * rad_c**3
                + 0.5e0 * beta[ii] * rad_c * numpy.log(rad_c)
            )

    return sigr, sigt, r_deflect, rradius


@numba.njit(cache=True)
def eyoung_parallel_array(n, eyoung_j_in, a_in, poisson_j_perp_in):
    """
    Author : C. Swanson, PPPL
    January 2022
    See Issue #1205 for derivation PDF
    This subroutine gives the smeared elastic properties of two
    members that are carrying a force in parallel with each other.
    The force goes in direction j.
    Members 1 and 2 are the individual members to be smeared.
    Member 3 is the effective smeared member (output).
    This is pretty easy because the smeared properties are simply
    the average weighted by the cross-sectional areas perpendicular
    to j.
    The assumption is that the strains in j are equal.
    If you're dealing with anisotropy, please pay attention to the
    fact that the specific Young's Modulus used here is that in
    the j direction, and the specific Poisson's ratio used here is
    that between the j and transverse directions in that order.
    (transverse strain / j strain, under j stress)
    The smeared Poisson's ratio is computed assuming the transverse
    dynamics are isotropic, and that the two members are free to
    shrink/expand under Poisson effects without interference from
    each other. This may not be true of your case.

    To build up a composite smeared member of any number of
    individual members, you can pass the same properties for
    members 2 and 3, and call it successively, using the properties
    of each member as the first triplet of arguments. This way, the
    last triplet acts as a "working sum":
    call eyoung_parallel(triplet1, triplet2, tripletOUT)
    call eyoung_parallel(triplet3, tripletOUT, tripletOUT)
    call eyoung_parallel(triplet4, tripletOUT, tripletOUT)
    ... etc.
    So that tripletOUT would eventually have the smeared properties
    of the total composite member.
    """
    eyoung_j_out = 0
    a_out = 0
    poisson_j_perp_out = 0

    # Parallel-composite them all together
    for ii in range(n):
        eyoung_j_out, a_out, poisson_j_perp_out = eyoung_parallel(
            eyoung_j_in[ii],
            a_in[ii],
            poisson_j_perp_in[ii],
            eyoung_j_out,
            a_out,
            poisson_j_perp_out,
        )

    return eyoung_j_out, a_out, poisson_j_perp_out


def _inductance_factor(H, Ri, Ro, Rm, theta1):
    """Calculates the inductance factor for a toroidal structure
    using surrogate 2 #1866.

    Author: Timothy Nunn, UKAEA

    :param H: the maximum height of the structure
    :param Ri: the radius of the inboard side of the structure CCL
    :param Ro: the radius of the outboard side of the structure CCL
    :param Rm: the radius at which `H` occurs
    :param theta1: the polar angle of the point at which one circular arc is
    joined to another circular arc in the approximation to the structure CCL,
    using an arbitrary origin of coordinates (Rc2, Zc2).
    """
    # NOTE: the following properties are not those of the plasma but of
    # the VV/coil structures
    major_radius = (Ro + Ri) / 2
    minor_radius = (Ro - Ri) / 2

    aspect_ratio = major_radius / minor_radius
    triangularity = (major_radius - Rm) / minor_radius
    elongation = H / minor_radius

    return (
        4.933
        + 0.03728 * elongation
        + 0.06980 * triangularity
        - 3.551 * aspect_ratio
        + 0.7629 * aspect_ratio**2
        - 0.06298 * (theta1 / 90)
    )


def vv_stress_on_quench(
    # TF shape
    H_coil,
    Ri_coil,
    Ro_coil,
    Rm_coil,
    ccl_length_coil,
    theta1_coil,
    # VV shape
    H_vv,
    Ri_vv,
    Ro_vv,
    Rm_vv,
    theta1_vv,
    # TF properties
    n_tf,
    n_tf_turn,
    S_rp,
    S_cc,
    taud,
    I_op,
    # VV properties
    d_vv,
):
    """Generic model to calculate the Tresca stress of the
    Vacuum Vessel (VV), experienced when the TF coil quenches.

    Author: Timothy Nunn, UKAEA

    The current center line (CCL) of a structure is an appoximation
    of where the poloidal current of a structure acts. This model
    considers how the current (self-)induced in the following structures
    affects the Tresca stress on the VV (corresponding to the structure
    index used in [1]):

    0. the TF coil conductors
    1. the TF coil steel structure
    2. the VV

    :param H_coil: the maximum height of the TF coil CCL
    :param Ri_coil: the radius of the inboard edge of the TF coil CCL
    :param Ro_coil: the radius of the outboard edge of the TF coil CCL
    :param Rm_coil: the radius where the maximum height of the TF coil CCL occurs
    :param ccl_length_coil: the length of the TF coil CCL
    :param theta1_coil: the polar angle of the point at which one circular arc is
    joined to another circular arc in the approximation to the coil CCL,
    using an arbitrary origin of coordinates (Rc2, Zc2).

    :param H_vv: the maximum height of the VV CCL
    :param Ri_vv: the radius of the inboard edge of the VV CCL
    :param Ro_vv: the radius of the outboard edge of the VV CCL
    :param Rm_vv: the radius where the maximum height of the VV CCL occurs
    :param theta1_vv: the polar angle of the point at which one circular arc is
    joined to another circular arc in the approximation to the VV CCL,
    using an arbitrary origin of coordinates (Rc2, Zc2).

    :param n_tf: the number of TF coils
    :param n_tf_turn: the number of turns per TF coil
    :param S_rp: the cross-sectional area of the radial plates of the TF coil
    :param S_cc: the cross-sectional area of the TF coil case
    :param taud: the discharge time of the TF coil when quench occurs
    :param I_op: the 'normal' operating current of the TF coil

    :param d_vv: the thickness of the vacuum vessel

    :returns: the maximum stress experienced by the vacuum vessel

    Notes
    -----
    The theta1 quantity for the TF coil and VV is not very meaningful. The
    impact of it of the inductance is rather small. Generally, the paper seems to
    suggest the TF coil is between 40 and 60, as this is the range they calculate
    the surrogates over. No range is provided for the VV but the example using
    JA DEMO is 1 degree suggesting the quantity will be very small.

    References
    ----------
    1. ITOH, Yasuyuki & Utoh, Hiroyasu & SAKAMOTO, Yoshiteru & Hiwatari, Ryoji. (2020).
    Empirical Formulas for Estimating Self and Mutual Inductances of Toroidal Field Coils and Structures.
    Plasma and Fusion Research. 15. 1405078-1405078. 10.1585/pfr.15.1405078.
    """
    # Poloidal loop resistance (PLR) in ohms
    plr_coil = ((0.5 * ccl_length_coil) / (n_tf * (S_cc + S_rp))) * 1e-6
    plr_vv = ((0.84 / d_vv) * 0.94) * 1e-6

    # relevant self-inductances in henry (H)
    coil_structure_self_inductance = (
        (constants.rmu0 / numpy.pi)
        * H_coil
        * _inductance_factor(H_coil, Ri_coil, Ro_coil, Rm_coil, theta1_coil)
    )
    vv_self_inductance = (
        (constants.rmu0 / numpy.pi)
        * H_vv
        * _inductance_factor(H_vv, Ri_vv, Ro_vv, Rm_vv, theta1_vv)
    )

    # s^-1
    lambda0 = 1 / taud
    lambda1 = (plr_coil) / coil_structure_self_inductance
    lambda2 = (plr_vv) / vv_self_inductance

    # approximate time at which the maximum force (and stress) will occur on the VV
    tmaxforce = numpy.log((lambda0 + lambda1) / (2 * lambda0)) / (lambda1 - lambda0)

    I0 = I_op * numpy.exp(-lambda0 * tmaxforce)
    I1 = (
        lambda0
        * n_tf
        * n_tf_turn
        * I_op
        * (
            (numpy.exp(-lambda1 * tmaxforce) - numpy.exp(-lambda0 * tmaxforce))
            / (lambda0 - lambda1)
        )
    )
    I2 = (lambda1 / lambda2) * I1

    A_vv = (Ro_vv + Ri_vv) / (Ro_vv - Ri_vv)
    B_vvi = (constants.rmu0 * (n_tf * n_tf_turn * I0 + I1 + (I2 / 2))) / (
        2 * numpy.pi * Ri_vv
    )
    J_vvi = I2 / (2 * numpy.pi * d_vv * Ri_vv)

    zeta = 1 + ((A_vv - 1) * numpy.log((A_vv + 1) / (A_vv - 1)) / (2 * A_vv))

    tresca_stress_vv = zeta * B_vvi * J_vvi * Ri_vv

    return tresca_stress_vv
