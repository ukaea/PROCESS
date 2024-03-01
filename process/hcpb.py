import numpy as np

from process.fortran import (
    constants,
    blanket_library,
    ccfe_hcpb_module,
    build_variables,
    fwbs_variables,
    physics_variables,
    process_output as po,
    tfcoil_variables,
    heat_transport_variables,
    cost_variables,
    divertor_variables,
    buildings_variables,
    error_handling as eh,
    current_drive_variables,
    primary_pumping_variables,
    constraint_variables,
)
from process.coolprop_interface import FluidProperties


class CCFE_HCPB:
    """author: J Morris (UKAEA)

    This module contains the PROCESS CCFE HCPB blanket model
    based on CCFE HCPB model from the PROCESS engineering paper
    PROCESS Engineering paper (M. Kovari et al.)

    ### References
    - Kovari et al., Fusion Engineering and Design 104 (2016) 9-20
    """

    def __init__(self, blanket_library) -> None:
        self.outfile = constants.nout

        self.blanket_library = blanket_library

    def run(self, output: bool):
        ccfe_hcpb_module.ip = int(output)
        ccfe_hcpb_module.ofile = self.outfile

        # MDK (27/11/2015)
        build_variables.fwith = 2 * fwbs_variables.afw + 2 * fwbs_variables.fw_wall
        build_variables.fwoth = build_variables.fwith

        # Coolant type
        fwbs_variables.coolwh = 1
        # Note that the first wall coolant is now input separately.

        # Calculate blanket, shield, vacuum vessel and cryostat volumes
        blanket_library.component_volumes()

        # Centrepost neutronics
        if physics_variables.itart == 1:
            # CP radius at the point of maximum sield radius [m]
            # The maximum shield radius is assumed to be at the X-point
            r_sh_inboard_out_top = (
                physics_variables.rmajor
                - physics_variables.rminor * physics_variables.triang
                - 3 * build_variables.scrapli
            )

            # Half height of the CP at the largest shield radius [m]
            h_sh_max_r = physics_variables.rminor * physics_variables.kappa

            # Solid angle fraction of neutrons that hit the centrepost shield [-]
            # Calculating the CP solid angle coverage fraction
            # Rem : This calculation considered the shield flaring
            #       while the MCNP based neutronincs considers a
            #       cylindre
            f_geom_cp = self.st_cp_angle_fraction(
                h_sh_max_r,
                build_variables.r_sh_inboard_out,
                r_sh_inboard_out_top,
                physics_variables.rmajor,
            )

            # TF fast neutron flux (E > 0.1 MeV) [m^{-2}.s^{-1}]
            fwbs_variables.neut_flux_cp = self.st_tf_centrepost_fast_neut_flux(
                physics_variables.pneutmw,
                build_variables.shldith,
                physics_variables.rmajor,
            )

            # TF, shield and total CP nuclear heating [MW]
            (
                fwbs_variables.pnuc_cp_tf,
                fwbs_variables.pnuc_cp_sh,
                fwbs_variables.pnuc_cp,
            ) = self.st_centrepost_nuclear_heating(
                physics_variables.pneutmw, build_variables.shldith
            )

        else:  # No CP
            f_geom_cp = 0
            fwbs_variables.pnuc_cp_tf = 0
            fwbs_variables.pnuc_cp_sh = 0
            fwbs_variables.pnuc_cp = 0
            fwbs_variables.neut_flux_cp = 0

        blanket_library.component_masses()

        # Calculate the nuclear heating
        # Rem : The heating power will be normalized to the neutron power using
        #       the divertor and the centrepost (for itart == 1),
        self.nuclear_heating_magnets(output=output)
        self.nuclear_heating_fw()
        self.nuclear_heating_blanket()
        self.nuclear_heating_shield()
        self.nuclear_heating_divertor()

        # Normalisation of the nuclear heating
        # The nuclear heating are noramalized assuming no energy multiplication
        # in the divertor and the centrepost
        # Assume that all the neutrons are absorbed. (Not applicable for very thin blankets)
        # Rem SK : This calculation effectively only uses the angular fractions to get
        #          the energy multiplication and hence the power balance ...

        # Split neutron power to main wall between fw, bkt, shld and TF with same
        # fractions as before.
        # Total nuclear power deposited in the blancket sector (MW)
        ccfe_hcpb_module.pnuc_tot_blk_sector = (
            fwbs_variables.pnucfw
            + fwbs_variables.pnucblkt
            + fwbs_variables.pnucshld
            + fwbs_variables.ptfnuc
        )

        # Total nuclear power deposited in the
        # if ( pnuc_tot_blk_sector < 1.0d0 .or. pnuc_tot_blk_sector /= pnuc_tot_blk_sector ) then
        # #TODO This can flood the terminal, and should be logged once in Python
        # write(*,*)'pnucfw =', pnucfw, ' and ', 'pnucblkt =', pnucblkt
        # write(*,*)'pnucshld =', pnucshld, ' ptfnuc =', ptfnuc
        # end if

        # Solid angle fraction taken by the breeding blankets/shields
        f_geom_blanket = 1 - physics_variables.idivrt * fwbs_variables.fdiv - f_geom_cp

        # Power to the first wall (MW)
        fwbs_variables.pnucfw = (
            (fwbs_variables.pnucfw / ccfe_hcpb_module.pnuc_tot_blk_sector)
            * fwbs_variables.emult
            * f_geom_blanket
            * physics_variables.pneutmw
        )

        # Power to the blanket (MW)
        fwbs_variables.pnucblkt = (
            (fwbs_variables.pnucblkt / ccfe_hcpb_module.pnuc_tot_blk_sector)
            * fwbs_variables.emult
            * f_geom_blanket
            * physics_variables.pneutmw
        )

        # Power to the shield(MW)
        # The power deposited in the CP shield is added back in powerflow_calc
        fwbs_variables.pnucshld = (
            (fwbs_variables.pnucshld / ccfe_hcpb_module.pnuc_tot_blk_sector)
            * fwbs_variables.emult
            * f_geom_blanket
            * physics_variables.pneutmw
        )

        # Power to the TF coils (MW)
        # The power deposited in the CP conductor is added back here
        fwbs_variables.ptfnuc = (
            (fwbs_variables.ptfnuc / ccfe_hcpb_module.pnuc_tot_blk_sector)
            * fwbs_variables.emult
            * f_geom_blanket
            * physics_variables.pneutmw
            + fwbs_variables.pnuc_cp_tf
        )

        # Power deposited in the CP
        fwbs_variables.pnuc_cp_sh = (
            f_geom_cp * physics_variables.pneutmw - fwbs_variables.pnuc_cp_tf
        )

        # Old code kept for backward compatibility
        # ---
        # pnucdiv is not changed.
        # The energy due to multiplication, by subtraction:
        # emultmw = pnucfw + pnucblkt + pnucshld + ptfnuc + pnucdiv - pneutmw
        # ---

        # New code, a bit simpler
        fwbs_variables.emultmw = (
            (fwbs_variables.emult - 1) * f_geom_blanket * physics_variables.pneutmw
        )

        # powerflow calculation for pumping power
        self.powerflow_calc(output=output)

        # output
        if output:
            self.write_output()

    def nuclear_heating_magnets(self, output: bool):
        """Nuclear heating in the magnets for CCFE HCPB model
        author: Michael Kovari, CCFE, Culham Science Centre
        author: J. Morris, CCFE, Culham Science Centre
        This subroutine calculates the nuclear heating in the
        coils.
        PROCESS Engineering paper (M. Kovari et al.)
        """

        # Model factors and coefficients
        a = 2.830  # Exponential factor (m2/tonne)
        b = 0.583  # Exponential factor (m2/tonne)
        e = 9.062  # Pre-factor (1/kg). Corrected see issue #272

        # First wall void fractions

        # inboard FW coolant void fraction
        fwbs_variables.vffwi = (
            np.pi
            * fwbs_variables.afw**2
            / (fwbs_variables.pitch * build_variables.fwith)
        )

        # outboard FW coolant void fraction
        fwbs_variables.vffwo = fwbs_variables.vffwi

        # mean FW coolant void fraction
        vffwm = fwbs_variables.vffwi

        # Calculate smeared densities of blanket sections
        # gaseous He coolant in armour, FW & blanket: He mass is neglected
        ccfe_hcpb_module.armour_density = fwbs_variables.denw * (1.0 - vffwm)
        ccfe_hcpb_module.fw_density = fwbs_variables.denstl * (1.0 - vffwm)
        ccfe_hcpb_module.blanket_density = (
            fwbs_variables.whtblkt / fwbs_variables.volblkt
        )
        ccfe_hcpb_module.shield_density = (
            fwbs_variables.whtshld / fwbs_variables.volshld
        )
        # Picking the largest value for VV thickness
        d_vv_all = build_variables.d_vv_in
        if build_variables.d_vv_out > d_vv_all:
            d_vv_all = build_variables.d_vv_out

        if d_vv_all > 1.0e-6:
            ccfe_hcpb_module.vv_density = fwbs_variables.vvmass / fwbs_variables.vdewin
        else:
            ccfe_hcpb_module.vv_density = 0.0

        # Calculation of average blanket/shield thickness [m]
        if physics_variables.itart == 1:
            # There is no inner blanket for TART design [m]
            th_blanket_av = build_variables.blnkoth

            # The CP shield in considered in a separate calcualtion [m]
            th_shield_av = build_variables.shldoth

        else:
            # Average breeding blanket thickness [m]
            th_blanket_av = 0.5 * (build_variables.blnkoth + build_variables.blnkith)

            # Average neutronic shield thickness [m]
            th_shield_av = 0.5 * (build_variables.shldoth + build_variables.shldith)

        # Exponents (tonne/m2)
        # Blanket exponent (/1000 for kg -> tonnes)
        ccfe_hcpb_module.x_blanket = (
            ccfe_hcpb_module.armour_density * fwbs_variables.fw_armour_thickness
            + ccfe_hcpb_module.fw_density
            * (build_variables.fwith + build_variables.fwoth)
            / 2.0
            + ccfe_hcpb_module.blanket_density * th_blanket_av
        ) / 1000.0

        # Shield exponent(/1000 for kg -> tonnes)
        ccfe_hcpb_module.x_shield = (
            ccfe_hcpb_module.shield_density * th_shield_av
            + ccfe_hcpb_module.vv_density
            * (build_variables.d_vv_in + build_variables.d_vv_out)
            / 2.0
        ) / 1000.0

        # If spherical tokamak, this is outboard only. pnuc_cp_tf is evaluated separately
        if physics_variables.itart == 1:
            # Nuclear heating in outobard TF coil legs (whttflgs)
            # Unit heating (W/kg/GW of fusion power) x legs mass only (kg)
            ccfe_hcpb_module.tfc_nuc_heating = (
                e
                * np.exp(-a * ccfe_hcpb_module.x_blanket)
                * np.exp(-b * ccfe_hcpb_module.x_shield)
                * tfcoil_variables.whttflgs
            )
        else:
            # Nuclear heating in TF coil
            # Unit heating (W/kg/GW of fusion power) x total mass (kg)
            ccfe_hcpb_module.tfc_nuc_heating = (
                e
                * np.exp(-a * ccfe_hcpb_module.x_blanket)
                * np.exp(-b * ccfe_hcpb_module.x_shield)
                * tfcoil_variables.whttf
            )

        # Total heating (MW)
        fwbs_variables.ptfnuc = (
            ccfe_hcpb_module.tfc_nuc_heating
            * (physics_variables.powfmw / 1000.0)
            / 1.0e6
        )

        if output:
            po.oheadr(self.outfile, "Nuclear Heating Magnets Before Renormalisation")
            po.ovarre(
                self.outfile,
                "Shield line density (tonne/m2)",
                "(x_shield)",
                ccfe_hcpb_module.x_shield,
            )
            po.ovarre(
                self.outfile,
                "Blanket line density (tonne/m2)",
                "(x_blanket)",
                ccfe_hcpb_module.x_blanket,
            )
            po.ovarre(
                self.outfile,
                "Unit nuclear heating in TF coil (W/GW)",
                "(tfc_nuc_heating)",
                ccfe_hcpb_module.tfc_nuc_heating,
            )
            po.ovarre(
                self.outfile,
                "Total nuclear heating in TF coil (MW)",
                "(ptfnuc.)",
                fwbs_variables.ptfnuc,
            )
            po.ovarre(self.outfile, "powfmw", "(powfmw.)", physics_variables.powfmw)
            po.ovarre(
                self.outfile,
                "total mass of the TF coils (kg)",
                "(whttf)",
                tfcoil_variables.whttf,
            )

    def nuclear_heating_fw(self):
        """Nuclear heating in the FW for CCFE HCPB model
        author: J. Morris, CCFE, Culham Science Centre

        This subroutine calculates the nuclear heating in the FW
        """
        # Unit heating of FW and armour (W/kg per W of fusion power)
        ccfe_hcpb_module.fw_armour_u_nuc_heating = 6.25e-7

        # Total nuclear heating in FW (MW)
        fwbs_variables.pnucfw = (
            fwbs_variables.fwmass
            * ccfe_hcpb_module.fw_armour_u_nuc_heating
            * physics_variables.powfmw
        )

        if fwbs_variables.pnucfw < 0:
            raise RuntimeError(
                f"""Error in nuclear_heating_fw. {fwbs_variables.pnucfw = },
                {physics_variables.powfmw = }, {fwbs_variables.fwmass = }"""
            )

    def nuclear_heating_blanket(self):
        """Nuclear heating in the blanket for CCFE HCPB model
        author: J. Morris, CCFE, Culham Science Centre
        This subroutine calculates the nuclear heating in the blanket
        """
        # Blanket nuclear heating coefficient and exponent
        a = 0.764
        b = 2.476e-3  # 1/tonne

        # Mass of the blanket in tonnes
        mass = fwbs_variables.whtblkt / 1000

        # Total blanket nuclear heating (MW)
        ccfe_hcpb_module.exp_blanket = 1 - np.exp(-b * mass)
        fwbs_variables.pnucblkt = (
            physics_variables.powfmw * a * ccfe_hcpb_module.exp_blanket
        )

        if fwbs_variables.pnucblkt < 1:
            eh.fdiags[0] = fwbs_variables.pnucblkt
            eh.fdiags[1] = ccfe_hcpb_module.exp_blanket
            eh.fdiags[2] = physics_variables.powfmw
            eh.fdiags[3] = mass
            eh.report_error(274)

    def nuclear_heating_shield(self):
        """Nuclear heating in the shield for CCFE HCPB model
        author: J. Morris, CCFE, Culham Science Centre
        This subroutine calculates the nuclear heating in the shield
        """

        # Shield nuclear heating coefficients and exponents
        f = 6.88e2  # Shield nuclear heating coefficient (W/kg/W)
        g = 2.723  # Shield nuclear heating exponent m2/tonne
        h = 0.798  # Shield nuclear heating exponent m2/tonne

        # Calculation of average blanket/shield thickness [m]
        if physics_variables.itart == 1:
            # The CP shield in considered in a separate calcualtion
            th_shield_av = build_variables.shldoth
        else:
            # Average neutronic shield thickness [m]
            th_shield_av = 0.5 * (build_variables.shldoth + build_variables.shldith)

        # Decay length [m-2]
        y = (ccfe_hcpb_module.shield_density / 1000) * th_shield_av

        # Unit nuclear heating of shield (W/kg/GW of fusion power) x mass
        ccfe_hcpb_module.exp_shield1 = np.exp(-g * ccfe_hcpb_module.x_blanket)
        ccfe_hcpb_module.exp_shield2 = np.exp(-h * y)
        ccfe_hcpb_module.shld_u_nuc_heating = (
            fwbs_variables.whtshld
            * f
            * ccfe_hcpb_module.exp_shield1
            * ccfe_hcpb_module.exp_shield2
        )

        # Total nuclear heating in shield (MW)
        fwbs_variables.pnucshld = (
            ccfe_hcpb_module.shld_u_nuc_heating
            * (physics_variables.powfmw / 1000)
            / 1.0e6
        )

    def nuclear_heating_divertor(self):
        """Nuclear heating in the divertor for CCFE HCPB model
        author: J. Morris, CCFE, Culham Science Centre
        This subroutine calculates the nuclear heating in the divertor
        """
        # Unfortunately the divertor heating was not tallied in the neutronics calcs
        # Assume that all the neutron energy + energy multiplication is absorbed in the reactor +
        # coils. It turns out that emult is also approx constant, but this is not used. No energy
        # multiplication in the divertor

        # Overwrite global variable for fdiv 07/11/18 SIM: Removed having spoken to JM
        # fdiv = 0.115D0

        # Nuclear heating in the divertor just the neutron power times fdiv
        if physics_variables.idivrt == 2:
            # Double null configuration
            fwbs_variables.pnucdiv = (
                0.8 * physics_variables.powfmw * 2 * fwbs_variables.fdiv
            )
        else:
            # single null configuration
            fwbs_variables.pnucdiv = (
                0.8 * physics_variables.powfmw * fwbs_variables.fdiv
            )

        # No heating of the H & CD
        fwbs_variables.pnuchcd = 0.0

    def powerflow_calc(self, output: bool):
        """Calculations for powerflow
        author: J. Morris, CCFE, Culham Science Centre
        Calculations for powerflow
        """
        # Radiation power incident on divertor (MW)
        if physics_variables.idivrt == 2:
            # Double null configuration
            fwbs_variables.praddiv = (
                physics_variables.pradmw * 2.0 * fwbs_variables.fdiv
            )
        else:
            # single null configuration
            fwbs_variables.praddiv = physics_variables.pradmw * fwbs_variables.fdiv

        # Radiation power incident on HCD apparatus (MW)
        fwbs_variables.pradhcd = physics_variables.pradmw * fwbs_variables.fhcd

        # Radiation power incident on first wall (MW)
        fwbs_variables.pradfw = (
            physics_variables.pradmw - fwbs_variables.praddiv - fwbs_variables.pradhcd
        )

        # If we have chosen pressurised water as the blanket coolant, set the
        # coolant outlet temperature as 20 deg C below the boiling point
        if fwbs_variables.coolwh == 2:
            outlet_saturated_fluid_properties = FluidProperties.of(
                "Water", pressure=fwbs_variables.blpressure * 1.0e6, vapor_quality=0
            )
            fwbs_variables.outlet_temp = (
                outlet_saturated_fluid_properties.temperature - 20.0
            )  # in K

        # Surface heat flux on first wall (outboard and inboard) (MW)
        # All of the fast particle losses go to the outer wall.
        fwbs_variables.psurffwo = (
            fwbs_variables.pradfw * build_variables.fwareaob / build_variables.fwarea
            + current_drive_variables.porbitlossmw
            + physics_variables.palpfwmw
        )
        fwbs_variables.psurffwi = fwbs_variables.pradfw * (
            1 - build_variables.fwareaob / build_variables.fwarea
        )

        # primary_pumping == 0
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
            heat_transport_variables.htpmw_blkt = (
                heat_transport_variables.fpumpblkt * fwbs_variables.pnucblkt
            )
            heat_transport_variables.htpmw_shld = heat_transport_variables.fpumpshld * (
                fwbs_variables.pnucshld + fwbs_variables.pnuc_cp_sh
            )
            heat_transport_variables.htpmw_div = heat_transport_variables.fpumpdiv * (
                physics_variables.pdivt
                + fwbs_variables.pnucdiv
                + fwbs_variables.praddiv
            )

        elif fwbs_variables.primary_pumping == 2:
            # Calculate the required material properties of the FW and BB coolant.
            self.blanket_library.primary_coolant_properties(output=output)
            # Mechanical pumping power is calculated for first wall and blanket
            self.blanket_library.thermo_hydraulic_model(output)

            # For divertor and shield, mechanical pumping power is a fraction of thermal
            # power removed by coolant
            heat_transport_variables.htpmw_shld = heat_transport_variables.fpumpshld * (
                fwbs_variables.pnucshld + fwbs_variables.pnuc_cp_sh
            )
            heat_transport_variables.htpmw_div = heat_transport_variables.fpumpdiv * (
                physics_variables.pdivt
                + fwbs_variables.pnucdiv
                + fwbs_variables.praddiv
            )

        elif fwbs_variables.primary_pumping == 3:
            # Issue #503
            # Mechanical pumping power is calculated using specified pressure drop for
            # first wall and blanket circuit, including heat exchanger and pipes
            pfactor = (
                primary_pumping_variables.p_he
                / (primary_pumping_variables.p_he - primary_pumping_variables.dp_he)
            ) ** (
                (primary_pumping_variables.gamma_he - 1)
                / primary_pumping_variables.gamma_he
            )
            # N.B. Currenlty primary_pumping==3 uses seperate variables found in
            # primary_pumping_variables rather than fwbs_variables.
            # The pressure (p_he) is assumed to be the pressure at the
            # blanket inlet/pump oulet.
            # The pressures (found in fwbs_variables) for coolants using
            # primary_pumping==2 are assumed to be the pressure at the
            # blanket oulet/pump inlet. The equation below is used for primary_pumping==2:
            # pfactor = ((pressure+deltap)/pressure)**((gamma-1.0d0)/gamma)
            t_in_compressor = primary_pumping_variables.t_in_bb / pfactor
            dt_he = (
                primary_pumping_variables.t_out_bb - primary_pumping_variables.t_in_bb
            )
            fpump = t_in_compressor / (fwbs_variables.etaiso * dt_he) * (pfactor - 1)
            p_plasma = (
                fwbs_variables.pnucfw
                + fwbs_variables.psurffwi
                + fwbs_variables.psurffwo
                + fwbs_variables.pnucblkt
            )
            primary_pumping_variables.htpmw_fw_blkt = fpump / (1 - fpump) * p_plasma

            # For divertor and shield, mechanical pumping power is a fraction of thermal
            # power removed by coolant
            heat_transport_variables.htpmw_shld = heat_transport_variables.fpumpshld * (
                fwbs_variables.pnucshld + fwbs_variables.pnuc_cp_sh
            )
            heat_transport_variables.htpmw_div = heat_transport_variables.fpumpdiv * (
                physics_variables.pdivt
                + fwbs_variables.pnucdiv
                + fwbs_variables.praddiv
            )
            if output:
                po.oheadr(self.outfile, "Pumping for primary coolant (helium)")
                po.ovarre(
                    self.outfile,
                    "Pressure drop in FW and blanket coolant incl. hx and pipes (Pa)",
                    "(dp_he)",
                    primary_pumping_variables.dp_he,
                )
                po.ovarre(
                    self.outfile,
                    "Fraction of FW and blanket thermal power required for pumping",
                    "(fpump)",
                    fpump,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Total power absorbed by FW & blanket (MW)",
                    "(p_plasma)",
                    p_plasma,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Inlet temperature of FW & blanket coolant pump (K)",
                    "(t_in_compressor)",
                    t_in_compressor,
                    "OP ",
                )
                po.ovarre(
                    self.outfile,
                    "Coolant pump outlet/Inlet temperature of FW & blanket (K)",
                    "(t_in_bb)",
                    primary_pumping_variables.t_in_bb,
                )
                po.ovarre(
                    self.outfile,
                    "Outlet temperature of FW & blanket (K)",
                    "(t_out_bb)",
                    primary_pumping_variables.t_out_bb,
                )
                po.ovarre(
                    self.outfile,
                    "Mechanical pumping power for FW and blanket cooling loop including heat exchanger (MW)",
                    "(htpmw_fw_blkt)",
                    primary_pumping_variables.htpmw_fw_blkt,
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

    def st_cp_angle_fraction(self, h_cp_top, r_cp_mid, r_cp_top, rmajor):
        """Author : S. Kahn, CCFE, Culham science centre
        Estimates the CP angular solid angle coverage fration
        Equation (1-3) from
        ref : P. Guest THE-REVIEW OF SCIENTIFIC INSTRUMENTS, vol 32, n 2 (1960)
        Initial, but undocumented calculation kept as commented section
        without any talor expansion approximation

        :param h_cp_top: Centrepost shield half height [m]
        :param r_cp_top: Centrepost top radius [m]
        :param r_cp_mid: Centrepost mid-plane radius [m]
        :param rmajor: Plasma major radius [m]

        :returns: Solid angle fraction covered by the CP [-]
        """

        n_integral = 10

        # Initial calculation
        # -------------------
        # Kept as commented section for documentation
        # ! Fraction of neutrons that hit the centre post neutronic shield
        # f_geom_cp = cphalflen / sqrt(cphalflen**2 + (rmajor-r_cp_outer)**2 ) &
        #           * atan(r_cp_outer/(rmajor-r_cp_outer) )/pi
        # -------------------

        # Major radius normalized to the CP average radius [-]
        rho_maj = 2.0 * rmajor / (r_cp_mid + r_cp_top)

        # Average CP extent in the toroidal plane [rad]
        phy_cp = np.arcsin(1.0 / rho_maj)

        # toroidal plane infinitesimal angle used in the integral [rad]
        d_phy_cp = phy_cp / n_integral

        # CP solid angle integral using trapezoidal method
        phy_cp_calc = 0.0
        cp_sol_angle = 0.0

        for _ in range(n_integral):
            # Little tricks to avoild NaNs due to rounding
            int_calc_3 = 1.0 - rho_maj**2 * np.sin(phy_cp_calc) ** 2
            if int_calc_3 < 0.0:
                int_calc_3 = 0.0

            int_calc_1 = 1.0 / np.sqrt(
                h_cp_top**2
                + (rho_maj * np.cos(phy_cp_calc) - np.sqrt(int_calc_3)) ** 2
            )

            phy_cp_calc = phy_cp_calc + d_phy_cp

            # Little tricks to avoild NaNs due to rounding
            int_calc_3 = 1.0 - rho_maj**2 * np.sin(phy_cp_calc) ** 2
            if int_calc_3 < 0.0:
                int_calc_3 = 0.0

            int_calc_2 = 1.0 / np.sqrt(
                h_cp_top**2
                + (rho_maj * np.cos(phy_cp_calc) - np.sqrt(int_calc_3)) ** 2
            )

            cp_sol_angle = cp_sol_angle + d_phy_cp * 0.5 * (int_calc_1 + int_calc_2)

        cp_sol_angle = cp_sol_angle * 4.0 * h_cp_top

        # Solid angle fraction covered by the CP (OUTPUT) [-]
        return 0.25 * cp_sol_angle / np.pi

    def st_tf_centrepost_fast_neut_flux(self, pneutmw, sh_width, rmajor):
        """Author S Kahn
        Routine calculating the fast neutron (E > 0.1 MeV) flux reaching the TF
        at the centerpost. These calcualtion are made from a CP only MCNP fit
        with a variable tungsten carbyde shield with 13% water cooling. The
        TF size is kept constant in the MCNP runs in such a way tha it increases
        size.
        This subroutine uses an shielding length per decade (/10 drop in flux)
        of 16.6 cm, close to the "15 - 16 cm" of Menard et al. 2016.
        (This is an e-folding lenth of 7.22 cm.)

        :param pneutmw: neutron fusion power [MW]
        :param sh_width: Neutron shield width [m]
        :param rmajor: Plasma major radius [m]
        """
        # Fraction of fast neutrons originating from the outer wall reflection [-]
        f_neut_flux_out_wall = 1

        # Tungsten density may vary with different manufacturing processes.
        f_wc_density = 2

        # Fraction of steel structures
        f_steel_struct = 0.1

        # CP fast neutron flux (E > 0.1 MeV) [m^{-2}.s^}{-1}]
        neut_flux_cp = 0

        if tfcoil_variables.i_tf_sup == 1:
            # Effecting shield width, removing steel structures
            sh_width_eff = sh_width * (1.0 - f_steel_struct)

            # Fit [10^{-13}.cm^{-2}]
            neut_flux_cp = 5.835 * np.exp(-15.392 * sh_width_eff) + 39.70 * (
                sh_width_eff / rmajor
            ) * np.exp(-24.722 * sh_width_eff)

            # Units conversion [10^{-13}.cm^{-2}] -> [m^{-2}]
            neut_flux_cp = neut_flux_cp * 1.0e17

            # Scaling to the actual plasma neutron power
            neut_flux_cp = (
                f_wc_density * f_neut_flux_out_wall * neut_flux_cp * (pneutmw / 800)
            )

        return neut_flux_cp

    def st_centrepost_nuclear_heating(self, pneut, sh_width):
        """Author : P J Knight, CCFE, Culham Science Centre
        Author : S Kahn, CCFE, Culham Science Centre
        Estimates the nuclear power absorbed by the ST centrepost magnet
        This routine calculates the neutron power absorbed by a
        copper spherical tokamak centrepost.
        The calculation estimates the fraction of neutrons hitting
        the centrepost from a point source at the plasma centre,
        and assumes an average path length of 2*r_cp_outer, and an
        e-folding decay length of 0.08m (copper-water mixture).
        J D Galambos, STAR Code : Spherical Tokamak Analysis and Reactor Code,
        unpublished internal Oak Ridge document

        :param pneut: 14 MeV plasma neutron power generated by the plasma [MW]
        :param sh_width: Thickeness of the centrepost neutron shield [m]

        :returns: (Nuclear nuclear heat deposited in the centrepost TF coil [MW],
        Nuclear nuclear heat deposited in the centrepost shield [MW],
        Total nuclear heat deposited in the centrepost shield [MW])
        """
        # Outer wall reflection TF nuclear heating enhancement factor [-]
        f_pnuc_cp_tf = 1

        # Outer wall reflection shield nuclear heating enhancement factor [-]
        f_pnuc_cp_sh = 1.7

        # Tungsten density may vary with different manufacturing processes.
        f_wc_density = 2

        # Fraction of steel structures
        f_steel_struct = 0.1

        # Former nuclear heating calculations for Copper magnets
        # Commented out as no nuclear shielding was included
        # ---
        # ! Fraction of the nuclear power absorbed by the copper centrepost
        # ! (0.08 m e-folding decay length)
        # f_neut_absorb_cp = 1.0D0 - exp( -2.0D0*tfcth / 0.08D0)
        #
        # ! Nuclear power
        # pnuc_cp = pneut * f_geom_cp * f_neut_absorb_cp
        # ---

        # Steel support structure effective WC shield thickness reduction
        sh_width_eff = sh_width * (1 - f_steel_struct)

        # Aluminium CP
        # ------------
        # From Pfus = 1 GW ST MCNP neutronic calculations assuming
        # Tungsten carbyde with 13% water cooling fraction
        if tfcoil_variables.i_tf_sup == 2:
            pnuc_cp_tf = (pneut / 800) * np.exp(3.882 - 16.69 * sh_width_eff)

            # WARINING, this is an extraoilation from TF heat ...
            # DO NOT TRUST THIS VALUE !!
            pnuc_cp_sh = (pneut / 800.0) * np.exp(3.882) - pnuc_cp_tf
            # ------------

        # Superconducting / copper CP
        # ---------------------------
        # MCNP calculations made with a TF magnet model with very large WP
        # so the TF is mostly copper, making the calculation also valid for
        # Copper TF centrepost
        else:
            # This subroutine uses an shielding length per decade (/10 drop in neutron heating)
            # of 15.5 cm, within to the "15 - 16 cm" of Menard et al. 2016.
            # (This is an e-folding lenth of 6.72 cm.)

            # Nuclear powers fits for a 800 MW plasma neutron source
            # ***
            # Nuclear power deposited in the CP winding pack by gammas [MW]
            pnuc_cp_wp_gam = 16.3 * np.exp(
                -14.63 * sh_width_eff
            ) + 143.08 * sh_width_eff * (sh_width / physics_variables.rmajor) * np.exp(
                -21.747 * sh_width_eff
            )

            # Nuclear power deposited in the CP winding pack by neutrons [MW]
            pnuc_cp_wp_n = 1.403 * np.exp(
                -16.535 * sh_width_eff
            ) + 3.812 * sh_width_eff * (sh_width / physics_variables.rmajor) * np.exp(
                -23.631 * sh_width_eff
            )

            # Nuclear power deposited in the CP steel case by gammas [MW]
            pnuc_cp_case_gam = 1.802 * np.exp(
                -13.993 * sh_width_eff
            ) + 38.592 * sh_width * (sh_width_eff / physics_variables.rmajor) * np.exp(
                -27.051 * sh_width_eff
            )

            # Nuclear power deposited in the CP steel case by neutrons [MW]
            pnuc_cp_case_n = 0.158 * np.exp(
                -55.046 * sh_width_eff
            ) + 2.0742 * sh_width_eff * (sh_width / physics_variables.rmajor) * np.exp(
                -24.401 * sh_width_eff
            )

            # Nuclear power density deposited in the tungsten carbyde shield by photons [MW]
            pnuc_cp_sh_gam = sh_width_eff * (
                596 * np.exp(-4.130 * sh_width_eff)
                + 90.586 * np.exp(0.6837 * sh_width_eff)
            )

            # Nuclear power density deposited in the tungsten carbyde shield by neutrons [MW]
            pnuc_cp_sh_n = sh_width_eff * (
                202.10 * np.exp(-10.533 * sh_width_eff)
                + 80.510 * np.exp(-0.9801 * sh_width_eff)
            )

            # Fit generalisation
            # ***
            # Correction for the actual 14 MeV plasma neutron power
            pnuc_cp_wp_gam = (pneut / 800) * pnuc_cp_wp_gam
            pnuc_cp_wp_n = (pneut / 800) * pnuc_cp_wp_n
            pnuc_cp_case_gam = (pneut / 800) * pnuc_cp_case_gam
            pnuc_cp_case_n = (pneut / 800) * pnuc_cp_case_n
            pnuc_cp_sh_gam = (pneut / 800) * pnuc_cp_sh_gam
            pnuc_cp_sh_n = (pneut / 800) * pnuc_cp_sh_n

            # Correction for neutron reflected by the outer wall hitting the CP
            pnuc_cp_wp_gam = f_pnuc_cp_tf * pnuc_cp_wp_gam
            pnuc_cp_wp_n = f_pnuc_cp_tf * pnuc_cp_wp_n
            pnuc_cp_case_gam = f_pnuc_cp_tf * pnuc_cp_case_gam
            pnuc_cp_case_n = f_pnuc_cp_tf * pnuc_cp_case_n
            pnuc_cp_sh_gam = f_pnuc_cp_sh * pnuc_cp_sh_gam
            pnuc_cp_sh_n = f_pnuc_cp_sh * pnuc_cp_sh_n

            # TF nuclear heat [MW]
            pnuc_cp_tf = (
                pnuc_cp_wp_gam + pnuc_cp_wp_n + pnuc_cp_case_gam + pnuc_cp_case_n
            )

            # Tungsten density correction
            pnuc_cp_tf = pnuc_cp_tf * f_wc_density

            # Shield nuclear heat [MW]
            pnuc_cp_sh = pnuc_cp_sh_gam + pnuc_cp_sh_n

        # Total CP nuclear heat [MW]
        pnuc_cp = pnuc_cp_tf + pnuc_cp_sh

        return pnuc_cp_tf, pnuc_cp_sh, pnuc_cp

    def tbr_shimwell(self, breeder_f, li6enrich, iblanket_thickness, output: bool):
        """Calculates TBR
        author: Michael Kovari
        breeder_f   : input real : Volume of Li4SiO4 / (Volume of Be12Ti + Li4SiO4)
        li6enrich   : input real : lithium-6 enrichment (%)
        iblanket_thickness   : input integer : blanket thickness switch
        tbr         : output real : 5-year time-averaged tritium breeding ratio
        """
        # for the v array of expansion terms:
        # the first element is for a thin blanket,
        # the second element is for a medium blanket
        # the third element is for a thick blanket
        v1 = [1.93920586301, 1.96122608615, 1.95893103797]
        v2 = [-0.948494854004, -0.860855681012, -0.809792727863]
        v3 = [-0.0186700302911, 0.0193393390622, 0.016958778333]
        v4 = [0.483417432982, 0.279977226537, -0.120230857418]
        v5 = [0.785901227724, 0.659918133027, 0.461211316443]
        v6 = [-0.0120169189644, 0.013070435947, -0.0478789050674]
        v7 = [-3.45723121388, -3.48450356973, -2.1978304461]
        v8 = [-2.05212472576, -2.3360647329, -1.38785787744]
        v9 = [6.45375263346, 7.38314099334, 4.93883798388]
        v10 = [-0.436421277881, -0.365511595682, -0.223668963335]
        v11 = [0.0129809166177, -0.0181287662329, 0.0178181886132]
        v12 = [2.26116309299, 2.30397890094, 1.42583418972]
        v13 = [-3.87538808736, -4.37481611533, -2.80720698559]
        v14 = [1.05778783291, 1.30804004777, 0.814691647096]
        v15 = [-3.12644013943, -3.71450110227, -2.48568193656]
        v16 = [1.86242247177, 2.1588023402, 1.37932384899]
        v17 = [0.253324925437, 0.253324925437, 0.253355839249]
        v18 = [0.18795823903, 0.198976219881, 0.190845918447]
        v19 = [-0.0256707269253, -0.0192924115968, -0.0257699008284]

        y = li6enrich / 100
        tbr = (
            v1[iblanket_thickness - 1]
            + v2[iblanket_thickness - 1] * breeder_f
            + v3[iblanket_thickness - 1] * y
            + v4[iblanket_thickness - 1] * y * breeder_f
            + v5[iblanket_thickness - 1] * breeder_f**2
            + v6[iblanket_thickness - 1] * y**2
            + v7[iblanket_thickness - 1] * breeder_f**2 * y
            + v8[iblanket_thickness - 1] * breeder_f * y**2
            + v9[iblanket_thickness - 1] * breeder_f**2 * y**2
            + v10[iblanket_thickness - 1] * breeder_f**3
            + v11[iblanket_thickness - 1] * y**3
            + v12[iblanket_thickness - 1] * y * breeder_f**3
            + v13[iblanket_thickness - 1] * y**2 * breeder_f**3
            + v14[iblanket_thickness - 1] * y**3 * breeder_f
            + v15[iblanket_thickness - 1] * y**3 * breeder_f**2
            + v16[iblanket_thickness - 1] * y**3 * breeder_f**3
            + v17[iblanket_thickness - 1] * np.log(breeder_f)
            + v18[iblanket_thickness - 1] * np.log(y)
            + v19[iblanket_thickness - 1] * np.log(breeder_f) * np.log(y)
        )

        if output:
            po.ovarrf(
                self.outfile, "Lithium-6 enrichment (%)", "(li6enrich)", li6enrich
            )
            po.ovarrf(
                self.outfile,
                "Breeder fraction by volume: Li4SiO4/(Be12Ti+Li4SiO4)",
                "(breeder_f)",
                breeder_f,
            )
            if iblanket_thickness == 1:
                po.ovarin(
                    self.outfile,
                    "Blanket thickness choice: THIN (0.53 m inboard, 0.91 m outboard)",
                    "[iblanket_thickness-1]",
                    iblanket_thickness,
                )
            elif iblanket_thickness == 2:
                po.ovarin(
                    self.outfile,
                    "Blanket thickness choice: MEDIUM (0.64 m inboard, 1.11 m outboard)",
                    "[iblanket_thickness-1]",
                    iblanket_thickness,
                )
            elif iblanket_thickness == 3:
                po.ovarin(
                    self.outfile,
                    "Blanket thickness choice: THICK (0.75 m inboard, 1.30 m outboard)",
                    "[iblanket_thickness-1]",
                    iblanket_thickness,
                )
            po.ovarrf(
                self.outfile,
                "Tritium breeding ratio (5-year time-averaged)",
                "(tbr)",
                tbr,
                "OP ",
            )
            po.ovarrf(
                self.outfile,
                "Minimum Tritium breeding ratio",
                "(tbrmin)",
                constraint_variables.tbrmin,
            )

            po.ocmmnt(
                self.outfile,
                '(See "A parameter study of time-varying tritium production in solid-type breeder blankets,',
            )
            po.ocmmnt(self.outfile, "J. Shimwell et al, Fusion Engineering and Design")
            po.ovarre(
                self.outfile,
                "For consistency, inboard first wall thicknesses should be 0.03 (m)",
                "(fwith)",
                build_variables.fwith,
            )
            po.ovarre(
                self.outfile,
                "For consistency, outboard first wall thicknesses should be 0.03 (m)",
                "(fwoth)",
                build_variables.fwoth,
            )
            po.ovarre(
                self.outfile,
                "For consistency, first wall armour thickness should be 0.003 (m)",
                "(fw_armour_thickness)",
                fwbs_variables.fw_armour_thickness,
            )

        return tbr

    def write_output(self):
        po.oheadr(self.outfile, "First wall and blanket : CCFE HCPB model")
        po.osubhd(self.outfile, "Blanket Composition by volume :")

        po.ovarrf(
            self.outfile,
            "Titanium beryllide fraction",
            "(fbltibe12)",
            fwbs_variables.fbltibe12,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Lithium orthosilicate fraction",
            "(fblli2sio4)",
            fwbs_variables.fblli2sio4,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Steel fraction",
            "(fblss_ccfe)",
            fwbs_variables.fblss_ccfe,
            "OP ",
        )
        po.ovarrf(self.outfile, "Coolant fraction", "(vfcblkt)", fwbs_variables.vfcblkt)
        po.ovarrf(
            self.outfile, "Purge gas fraction", "(vfpblkt)", fwbs_variables.vfpblkt
        )

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

        po.osubhd(self.outfile, "Component Masses :")

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
            "Blanket Mass - Total(kg)",
            "(whtblkt)",
            fwbs_variables.whtblkt,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "    Blanket Mass - TiBe12 (kg)",
            "(whtbltibe12)",
            fwbs_variables.whtbltibe12,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "    Blanket Mass - Li4SiO4 (kg)",
            "(whtblli4sio4)",
            fwbs_variables.whtblli4sio4,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "    Blanket Mass - Steel (kg)",
            "(whtblss)",
            fwbs_variables.whtblss,
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
            self.outfile, "Shield Mass (kg)", "(whtshld)", fwbs_variables.whtshld, "OP "
        )
        po.ovarre(
            self.outfile,
            "Vacuum vessel mass (kg)",
            "(vvmass)",
            fwbs_variables.vvmass,
            "OP ",
        )

        #  Nuclear heating section
        po.osubhd(self.outfile, "Nuclear heating :")

        #  ST centre post
        if physics_variables.itart == 1:
            if tfcoil_variables.i_tf_sup == 0:
                po.osubhd(self.outfile, "(Copper resistive centrepost used)")
            elif tfcoil_variables.i_tf_sup == 1:
                po.osubhd(self.outfile, "(Superdonducting magnet centrepost used)")
                po.ovarre(
                    self.outfile,
                    "ST centrepost TF fast neutron fllux (E > 0.1 MeV) (m^(-2).s^(-1))",
                    "(neut_flux_cp)",
                    fwbs_variables.neut_flux_cp,
                    "OP ",
                )
            elif tfcoil_variables.i_tf_sup == 2:
                po.osubhd(self.outfile, "(Aluminium magnet centrepost used)")

            po.ovarre(
                self.outfile,
                "ST centrepost TF heating (MW)",
                "(pnuc_cp_tf)",
                fwbs_variables.pnuc_cp_tf,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ST centrepost shield heating (MW)",
                "(pnuc_cp_sh)",
                fwbs_variables.pnuc_cp_sh,
                "OP ",
            )
            po.ovarre(
                self.outfile,
                "ST centrepost total heating (MW)",
                "(pnuc_cp)",
                fwbs_variables.pnuc_cp,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Total nuclear heating in TF+PF coils (CS is negligible) (MW)",
            "(ptfnuc)",
            fwbs_variables.ptfnuc,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total nuclear heating in FW (MW)",
            "(pnucfw)",
            fwbs_variables.pnucfw,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Total nuclear heating in the blanket (including emult) (MW)",
            "(pnucblkt)",
            fwbs_variables.pnucblkt,
            "OP ",
        )
        po.ocmmnt(self.outfile, "(Note: emult is fixed for this model inside the code)")
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
        po.osubhd(self.outfile, " Diagostic output for nuclear heating :")
        po.ovarre(
            self.outfile,
            "Blanket exponential factor",
            "(exp_blanket)",
            ccfe_hcpb_module.exp_blanket,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Shield: first exponential",
            "(exp_shield1)",
            ccfe_hcpb_module.exp_shield1,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Shield: second exponential",
            "(exp_shield2)",
            ccfe_hcpb_module.exp_shield2,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Solid angle fraction taken by on divertor",
            "(fdiv)",
            fwbs_variables.fdiv,
        )

        po.ovarin(
            self.outfile,
            "Switch for plant secondary cycle ",
            "(secondary_cycle)",
            fwbs_variables.secondary_cycle,
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

        if fwbs_variables.primary_pumping != 3:
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
            po.ovarre(
                self.outfile,
                "Total electrical coolant pumping power: first wall, blanket, shield and divertor (MW)",
                "(htpmw)",
                heat_transport_variables.htpmw,
                "OP ",
            )

        po.ovarre(
            self.outfile,
            "Allowable nominal neutron fluence at first wall (MW.year/m2)",
            "(abktflnc)",
            cost_variables.abktflnc,
        )
        po.ovarin(
            self.outfile,
            "No of inboard blanket modules poloidally",
            "(nblktmodpi)",
            fwbs_variables.nblktmodpi,
        )
        po.ovarin(
            self.outfile,
            "No of inboard blanket modules toroidally",
            "(nblktmodti)",
            fwbs_variables.nblktmodti,
        )
        po.ovarin(
            self.outfile,
            "No of outboard blanket modules poloidally",
            "(nblktmodpo)",
            fwbs_variables.nblktmodpo,
        )
        po.ovarin(
            self.outfile,
            "No of outboard blanket modules toroidally",
            "(nblktmodto)",
            fwbs_variables.nblktmodto,
        )
        po.ovarre(
            self.outfile,
            "Isentropic efficiency of first wall / blanket coolant pumps",
            "(etaiso)",
            fwbs_variables.etaiso,
        )

        po.osubhd(self.outfile, "Other volumes, masses and areas :")
        po.ovarre(
            self.outfile,
            "First wall area (m2)",
            "(fwarea)",
            build_variables.fwarea,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Cryostat internal radius (m)",
            "(rdewex)",
            fwbs_variables.rdewex,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Cryostat internal half-height (m)",
            "(zdewex)",
            fwbs_variables.zdewex,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Vertical clearance from TF coil to cryostat (m)",
            "(clh1)",
            buildings_variables.clh1,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Divertor area (m2)",
            "(divsur)",
            divertor_variables.divsur,
            "OP ",
        )
        po.ovarre(
            self.outfile,
            "Divertor mass (kg)",
            "(divmas)",
            divertor_variables.divmas,
            "OP ",
        )
