import logging

import numpy as np

from process import constants
from process import process_output as po
from process.data_structure import (
    build_variables,
    buildings_variables,
    cost_variables,
    current_drive_variables,
    divertor_variables,
    fwbs_variables,
    heat_transport_variables,
    pfcoil_variables,
    physics_variables,
    tfcoil_variables,
)

logger = logging.getLogger(__name__)


class Buildings:
    """

    This module contains routines for calculating the
    """

    def __init__(self):
        """

        This routine calls the buildings calculations.
        """
        self.outfile = constants.NOUT  # output file unit

    def run(self, output: bool = False):
        # Find TF coil radial positions
        # outboard edge: outboard mid-leg radial position + half-thickness of outboard leg
        tfro = build_variables.r_tf_outboard_mid + (
            build_variables.dr_tf_outboard * 0.5e0
        )
        # inboard edge: inboard mid-leg radial position - half-thickness of inboard leg
        tfri = build_variables.r_tf_inboard_mid - (build_variables.dr_tf_inboard * 0.5e0)

        # Find width, in radial dimension, of TF coil (m)
        tf_radial_dim = tfro - tfri

        # Find full height of TF coil (m)
        #  = 2 * (mid-plane to TF coil inside edge + thickness of coil)
        tf_vertical_dim = 2.0e0 * (
            build_variables.z_tf_inside_half + build_variables.dr_tf_outboard
        )

        # Find mass of each TF coil, in tonnes
        tfmtn = 1.0e-3 * tfcoil_variables.m_tf_coils_total / tfcoil_variables.n_tf_coils

        # Calculate building areas and volumes

        if buildings_variables.i_bldgs_size == 1:
            # Updated building estimates
            self.bldgs_sizes(output, tf_radial_dim, tf_vertical_dim)

        else:
            # Previous estimation work
            (
                buildings_variables.cryvol,
                buildings_variables.volrci,
                buildings_variables.rbvol,
                buildings_variables.rmbvol,
                buildings_variables.wsvol,
                buildings_variables.elevol,
            ) = self.bldgs(
                output,
                pfcoil_variables.r_pf_coil_outer_max,
                pfcoil_variables.m_pf_coil_max,
                tfro,
                tfri,
                tf_vertical_dim,
                tfmtn,
                tfcoil_variables.n_tf_coils,
                build_variables.r_shld_outboard_outer,
                build_variables.r_shld_inboard_inner,
                2.0e0
                * (build_variables.z_tf_inside_half - build_variables.dz_shld_vv_gap)
                - build_variables.dz_vv_upper
                - build_variables.dz_vv_lower,
                fwbs_variables.whtshld,
                fwbs_variables.r_cryostat_inboard,
                heat_transport_variables.helpow,
            )

    def bldgs(
        self,
        output: bool,
        pfr,
        pfm,
        tfro,
        tfri,
        tfh,
        tfm,
        n_tf_coils,
        shro,
        shri,
        shh,
        shm,
        crr,
        helpow,
    ):
        """Determines the sizes of the plant buildings


        pfr : input/output real :  largest PF coil outer radius, m
        pfm : : input real : largest PF coil mass, tonne
        tfro : input real : outer radius of TF coil, m
        tfri : input real : inner radius of TF coil, m
        tfh : input real : full height of TF coil, m
        tfm : input real : mass of one TF coil, tonne
        tfno : input real : number of tf coils
        shro : input real : outer radius of attached shield, m
        shri : input real : inner radius of attached shield, m
        shh : input real : height of attached shield, m
        shm : input real : total mass of attached shield, kg
        crr : input real : outer radius of common cryostat, m
        helpow : input real : total cryogenic load, W
        iprint : input integer : switch for writing to output file (1=yes)
        outfile : input integer : output file unit
        cryv : output real : volume of cryogenic building, m3
        vrci : output real : inner volume of reactor building, m3
        rbv : output real : outer volume of reactor building, m3
        rmbv : output real : volume of reactor maintenance building, m3
        wsv : output real : volume of warm shop, m3
        elev : output real : volume of electrical buildings, m3
        This routine determines the size of the plant buildings.
        The reactor building and maintenance building are sized
        based on the tokamak dimensions. The cryogenic building volume is
        scaled based on the total cryogenic load. The other building
        sizes are input from other modules or by the user.
        This routine was modified to include fudge factors (fac1,2,...)
        to fit the ITER design, September 1990 (J. Galambos).
        This routine was included in PROCESS in January 1992 by
        P. C. Shipe.

        Parameters
        ----------
        output: bool :

        pfr :

        pfm :

        tfro :

        tfri :

        tfh :

        tfm :

        n_tf_coils :

        shro :

        shri :

        shh :

        shm :

        crr :

        helpow :

        """

        # Reactor building

        # Determine basic machine radius (m)
        # crr  :  cryostat radius (m)
        # pfr  :  radius of largest PF coil (m)
        # tfro :  outer radius of TF coil (m)
        bmr = max(crr, pfr, tfro)

        # Determine largest transported piece
        sectl = shro - shri  # Shield thicknes (m)
        coill = tfro - tfri  # TF coil thickness (m)
        sectl = max(coill, sectl)

        # Calculate half width of building (m)
        # rxcl : clearance around reactor, m
        # trcl : transportation clearance between components, m
        # row  : clearance to building wall for crane operation, m
        # 19.48258241468535 + 4 + max(13.764874193548387 - 4.7423258064516141, 17.123405859443331 - 2.9939411851091102) + 1 + 4 = 42.61204708901957
        buildings_variables.wrbi = (
            bmr
            + buildings_variables.rxcl
            + sectl
            + buildings_variables.trcl
            + buildings_variables.row
        )

        # Calculate length to allow PF or cryostat laydown (m)

        # Laydown length (m)
        layl = max(crr, pfr)

        # Diagnoal length (m)
        hy = bmr + buildings_variables.rxcl + sectl + buildings_variables.trcl + layl

        # Angle between diagnoal length and floor (m)
        ang = (buildings_variables.wrbi - buildings_variables.trcl - layl) / hy

        # Cap angle at 1
        if abs(ang) > 1.0e0:
            ang = abs(ang) / ang

        # Length to allow laydown (m)
        drbi = (
            buildings_variables.trcl
            + layl
            + hy * np.sin(np.arccos(ang))
            + buildings_variables.wrbi
        )

        # Crane height based on maximum lift (m)
        # wgt : reactor building crane capacity (kg)
        #       Calculated if 0 is input
        # shmf : fraction of shield mass per TF coil to be moved in
        #        the maximum shield lift
        if buildings_variables.wgt > 1.0e0:
            wt = buildings_variables.wgt
        else:
            wt = buildings_variables.shmf * shm / n_tf_coils
            wt = max(wt, 1.0e3 * pfm, 1.0e3 * tfm)

        # Crane height (m)
        crcl = 9.41e-6 * wt + 5.1e0

        # Building height (m)
        # dz_tf_cryostat : clearance from TF coil to cryostat top, m
        # clh2 : clearance beneath TF coil to foundation, including basement, m
        # stcl : clearance above crane to roof, m
        # Additional tfh allows TF coil to be lifted right out
        hrbi = (
            buildings_variables.clh2
            + 2.0e0 * tfh
            + buildings_variables.dz_tf_cryostat
            + buildings_variables.trcl
            + crcl
            + buildings_variables.stcl
        )

        # Internal volume (m3)
        vrci = (
            buildings_variables.rbvfac * 2.0e0 * buildings_variables.wrbi * drbi * hrbi
        )
        try:
            assert vrci < np.inf
        except AssertionError:
            logger.exception("vrci is inf. Kludging to 1e10.")
            vrci = 1e10

        # External dimensions of reactor building (m)
        # rbwt : reactor building wall thickness, m
        # rbrt : reactor building roof thickness, m
        # fndt : foundation thickness, m
        rbw = 2.0e0 * buildings_variables.wrbi + 2.0e0 * buildings_variables.rbwt
        rbl = drbi + 2.0e0 * buildings_variables.rbwt
        rbh = hrbi + buildings_variables.rbrt + buildings_variables.fndt
        rbv = buildings_variables.rbvfac * rbw * rbl * rbh

        # Maintenance building
        # The reactor maintenance building includes the hot cells, the
        # decontamination chamber, the transfer corridors, and the waste
        # treatment building.  The dimensions of these areas are scaled
        # from a reference design based on the shield sector size.

        # Transport corridor size
        # hcwt : hot cell wall thickness, m
        tcw = shro - shri + 4.0e0 * buildings_variables.trcl
        tcl = 5.0e0 * tcw + 2.0e0 * buildings_variables.hcwt

        # Decontamination cell size
        dcw = 2.0e0 * tcw + 1.0e0

        # Hot cell size
        # hccl : clearance around components in hot cell, m
        hcw = shro - shri + 3.0e0 * buildings_variables.hccl + 2.0e0
        hcl = 3.0e0 * (shro - shri) + 4.0e0 * buildings_variables.hccl + tcw

        # Maintenance building dimensions
        rmbw = hcw + dcw + 3.0e0 * buildings_variables.hcwt
        rmbl = hcl + 2.0e0 * buildings_variables.hcwt

        # Height
        # wgt2 : hot cell crane capacity (kg)
        #        Calculated if 0 is input
        if buildings_variables.wgt2 > 1.0e0:
            wgts = buildings_variables.wgt2
        else:
            wgts = buildings_variables.shmf * shm / n_tf_coils

        cran = 9.41e-6 * wgts + 5.1e0
        rmbh = (
            10.0e0
            + shh
            + buildings_variables.trcl
            + cran
            + buildings_variables.stcl
            + buildings_variables.fndt
        )
        tch = shh + buildings_variables.stcl + buildings_variables.fndt

        # Volume
        rmbv = buildings_variables.mbvfac * rmbw * rmbl * rmbh + tcw * tcl * tch

        # Warm shop and hot cell gallery
        wsa = (rmbw + 7.0e0) * 20.0e0 + rmbl * 7.0e0
        wsv = buildings_variables.wsvfac * wsa * rmbh

        # Cryogenic building volume
        cryv = 55.0e0 * helpow**0.5
        # Other building volumes
        # pibv : power injection building volume, m3
        # esbldgm3 is forced to be zero if no energy storage is required (i_pulsed_plant=0)
        elev = (
            buildings_variables.tfcbv
            + buildings_variables.pfbldgm3
            + buildings_variables.esbldgm3
            + buildings_variables.pibv
        )

        # Calculate effective floor area for ac power module
        buildings_variables.a_plant_floor_effective = (
            rbv
            + rmbv
            + wsv
            + buildings_variables.triv
            + elev
            + buildings_variables.conv
            + cryv
            + buildings_variables.admv
            + buildings_variables.shov
        ) / 6.0e0
        buildings_variables.admvol = buildings_variables.admv
        buildings_variables.shovol = buildings_variables.shov
        buildings_variables.convol = buildings_variables.conv

        # Total volume of nuclear buildings
        buildings_variables.volnucb = vrci + rmbv + wsv + buildings_variables.triv + cryv

        # Output !
        # !!!!!!!!!

        if output:
            po.oheadr(self.outfile, "Plant Buildings System")
            po.ovarre(
                self.outfile, "Internal volume of reactor building (m3)", "(vrci)", vrci
            )
            po.ovarre(
                self.outfile,
                "Dist from centre of torus to bldg wall (m)",
                "(wrbi)",
                buildings_variables.wrbi,
            )
            po.ovarre(
                self.outfile,
                "Effective floor area (m2)",
                "(a_plant_floor_effective)",
                buildings_variables.a_plant_floor_effective,
            )
            po.ovarre(self.outfile, "Reactor building volume (m3)", "(rbv)", rbv)
            po.ovarre(
                self.outfile, "Reactor maintenance building volume (m3)", "(rmbv)", rmbv
            )
            po.ovarre(self.outfile, "Warmshop volume (m3)", "(wsv)", wsv)
            po.ovarre(
                self.outfile,
                "Tritium building volume (m3)",
                "(triv)",
                buildings_variables.triv,
            )
            po.ovarre(self.outfile, "Electrical building volume (m3)", "(elev)", elev)
            po.ovarre(
                self.outfile,
                "Control building volume (m3)",
                "(conv)",
                buildings_variables.conv,
            )
            po.ovarre(self.outfile, "Cryogenics building volume (m3)", "(cryv)", cryv)
            po.ovarre(
                self.outfile,
                "Administration building volume (m3)",
                "(admv)",
                buildings_variables.admv,
            )
            po.ovarre(
                self.outfile, "Shops volume (m3)", "(shov)", buildings_variables.shov
            )
            po.ovarre(
                self.outfile,
                "Total volume of nuclear buildings (m3)",
                "(volnucb)",
                buildings_variables.volnucb,
            )

        return cryv, vrci, rbv, rmbv, wsv, elev

    def bldgs_sizes(self, output, tf_radial_dim, tf_vertical_dim):
        """Subroutine that estimates the sizes (footprints and volumes) of
        buildings within a fusion power plant.
        Some estimates are scaled with parameters of the fusion plant,
        some are based on engineering/specialist assumptions,
        some are derived from footprints/volumes based on
        assessment of other power plants and/or similar facilities.
        !!


        Parameters
        ----------
        output :

        tf_radial_dim :

        tf_vertical_dim :

        """

        buildings_total_vol = 0.0e0

        # Reactor building

        # Lateral size driven by radial width of largest component, from:
        #  PF coil max radius, cryostat radius, TF coil outer radius
        width_reactor_piece = max(
            pfcoil_variables.r_pf_coil_outer_max,
            fwbs_variables.r_cryostat_inboard,
            tf_radial_dim,
        )
        # Allow for biological shielding around reactor
        width_reactor_piece = width_reactor_piece + buildings_variables.bioshld_thk

        # Calculate key-width of building (m)
        # include radial width of largest component *twice*, to allow for construction;
        # include clearance around reactor, transportation clearance between components,
        # clearance to building wall for crane operation
        key_width = (
            (2.0e0 * width_reactor_piece)
            + buildings_variables.reactor_clrnc
            + buildings_variables.transp_clrnc
            + buildings_variables.crane_clrnc_h
        )

        # Width of reactor building
        # allows for laydown of large components during construction
        buildings_variables.reactor_hall_w = 3.0e0 * key_width

        # Length of reactor building
        buildings_variables.reactor_hall_l = 3.0e0 * key_width

        # Calculate vertical clearance required (above and below reactor):
        # include clearance around reactor, transportation clearance between components,
        # clearance from TF coil to cryostat, clearance beneath TF coil,
        # clearance to ceiling for crane operation, crane arm height
        height_clrnc = (
            buildings_variables.reactor_clrnc
            + buildings_variables.transp_clrnc
            + buildings_variables.cryostat_clrnc
            + buildings_variables.ground_clrnc
            + buildings_variables.crane_clrnc_h
            + buildings_variables.crane_arm_h
        )

        # Height of reactor building
        # include height of TF coil *twice*, to allow for construction/maintenance
        buildings_variables.reactor_hall_h = (2.0e0 * tf_vertical_dim) + height_clrnc

        # Heating and Current Drive facility
        # Dimensions based upon estimates from M. Henderson, HCD Development Group
        # current_drive_variables.i_hcd_primary = switch for current drive model
        if (current_drive_variables.i_hcd_primary == 5) or (
            current_drive_variables.i_hcd_primary == 8
        ):
            # NBI technology will be situated within the reactor building
            buildings_variables.reactor_hall_l = (
                buildings_variables.reactor_hall_l
                + buildings_variables.nbi_sys_l
                + buildings_variables.reactor_clrnc
                + buildings_variables.transp_clrnc
            )
            buildings_variables.reactor_hall_w = (
                buildings_variables.reactor_hall_w
                + buildings_variables.nbi_sys_w
                + buildings_variables.reactor_clrnc
                + buildings_variables.transp_clrnc
            )
            hcd_building_area = 0.0e0
            hcd_building_vol = 0.0e0
        else:
            # Assume external building designed for EC or EBW is appropriate
            hcd_building_area = (
                buildings_variables.hcd_building_l * buildings_variables.hcd_building_w
            )
            hcd_building_vol = hcd_building_area * buildings_variables.hcd_building_h

        # Fuel Cycle facilities: include within reactor building
        # Dimensions based upon estimates from W. Smith
        buildings_variables.reactor_hall_l = (
            buildings_variables.reactor_hall_l + buildings_variables.fc_building_l
        )
        buildings_variables.reactor_hall_w = (
            buildings_variables.reactor_hall_w + buildings_variables.fc_building_w
        )

        # Reactor hall internal footprint and volume
        reactor_hall_area = (
            buildings_variables.reactor_hall_l * buildings_variables.reactor_hall_w
        )
        reactor_hall_vol = reactor_hall_area * buildings_variables.reactor_hall_h

        # Reactor building external footprint and volume
        reactor_building_l = (
            buildings_variables.reactor_hall_l
            + 2.0e0 * buildings_variables.reactor_wall_thk
        )
        reactor_building_w = (
            buildings_variables.reactor_hall_w
            + 2.0e0 * buildings_variables.reactor_wall_thk
        )
        reactor_building_h = (
            buildings_variables.reactor_hall_h
            + buildings_variables.reactor_roof_thk
            + buildings_variables.reactor_fndtn_thk
        )

        reactor_building_area = reactor_building_l * reactor_building_w

        reactor_building_vol = reactor_building_area * reactor_building_h

        # Reactor maintenance basement and tunnel
        # Architecture proposed here is a basement directly beneath the reactor enabling the
        # downwards extraction of hot components. The footprint estimated here is oversized to
        # include allowance for a tunnel to the hot cell storage/maintenance building.
        reactor_basement_l = buildings_variables.reactor_hall_w
        reactor_basement_w = buildings_variables.reactor_hall_w
        reactor_basement_area = reactor_basement_l * reactor_basement_w

        # basement height still includes some clearances
        reactor_basement_h = (
            tf_vertical_dim
            + buildings_variables.transp_clrnc
            + buildings_variables.crane_clrnc_h
            + buildings_variables.crane_arm_h
        )

        reactor_basement_vol = reactor_basement_area * reactor_basement_h

        reactor_build_totvol = reactor_building_vol + reactor_basement_vol

        buildings_total_vol = reactor_hall_vol + reactor_basement_vol

        # Hot Cell Facility
        # Provides hot cell facilities to maintain or dismantle highly radioactive components.
        # These are simplifications of R. Gowland's estimates of Operational Active Waste Storage,
        # which assumes all in-vessel components used through the life of the plant will need storage.
        # The storage area required is derived from the sizes and number of components, allowing
        # for a margin in component numbers as set by the quantity safety factor (buildings_variables.qnty_sfty_fac).
        # Footprints and volumes required for storage include hot separation distance (buildings_variables.hot_sepdist).

        # Assumptions:
        # tokomak is toroidally segmented based on number of TF coils (tfcoil_variables.n_tf_coils);
        # component will be stored with the largest dimension oriented horizontally;
        # height is the largest dimension;
        # if a component lifetime == 0, that component is not in the current machine build.

        # Inboard 'component': shield, blanket, first wall:
        # find height, maximum radial dimension, maximum toroidal dimension
        if cost_variables.life_plant != 0.0e0:
            hcomp_height = 2 * (
                build_variables.z_tf_inside_half
                - (
                    build_variables.dr_tf_inboard
                    + build_variables.dr_tf_shld_gap
                    + build_variables.dz_shld_thermal
                    + build_variables.dz_shld_vv_gap
                )
            )
            hcomp_rad_thk = (
                build_variables.dr_shld_inboard
                + build_variables.dr_blkt_inboard
                + build_variables.dr_fw_inboard
            )
            hcomp_tor_thk = (
                2
                * np.pi
                * (
                    physics_variables.rmajor
                    - (
                        physics_variables.rminor
                        + build_variables.dr_fw_plasma_gap_inboard
                        + build_variables.dr_fw_inboard
                        + build_variables.dr_blkt_inboard
                        + build_variables.dr_shld_inboard
                    )
                )
            ) / tfcoil_variables.n_tf_coils
            # find footprint and volume for storing component
            hcomp_footprint = (hcomp_height + buildings_variables.hot_sepdist) * (
                max(hcomp_rad_thk, hcomp_tor_thk) + buildings_variables.hot_sepdist
            )
            hcomp_vol = hcomp_footprint * (
                min(hcomp_rad_thk, hcomp_tor_thk) + buildings_variables.hot_sepdist
            )
            # required lifetime supply of components =
            #   ( number in build * (plant lifetime / component lifetime) ) * quantity safety factor
            hcomp_req_supply = (
                tfcoil_variables.n_tf_coils
                * (cost_variables.life_plant / cost_variables.life_plant)
            ) * buildings_variables.qnty_sfty_fac
            # total storage space for required supply of inboard shield-blanket-wall
            ib_hotcell_vol = hcomp_req_supply * hcomp_vol

            # Outboard 'component': first wall, blanket, shield
            hcomp_height = 2 * (
                build_variables.z_tf_inside_half
                - (
                    build_variables.dr_tf_inboard
                    + build_variables.dr_tf_shld_gap
                    + build_variables.dz_shld_thermal
                    + build_variables.dz_shld_vv_gap
                )
            )
            hcomp_rad_thk = (
                build_variables.dr_fw_outboard
                + build_variables.dr_blkt_outboard
                + build_variables.dr_shld_outboard
            )
            hcomp_tor_thk = (
                2
                * np.pi
                * (
                    physics_variables.rmajor
                    + physics_variables.rminor
                    + build_variables.dr_fw_plasma_gap_outboard
                    + build_variables.dr_fw_outboard
                    + build_variables.dr_blkt_outboard
                    + build_variables.dr_shld_outboard
                )
            ) / tfcoil_variables.n_tf_coils
            hcomp_footprint = (hcomp_height + buildings_variables.hot_sepdist) * (
                max(hcomp_rad_thk, hcomp_tor_thk) + buildings_variables.hot_sepdist
            )
            hcomp_vol = hcomp_footprint * (
                min(hcomp_rad_thk, hcomp_tor_thk) + buildings_variables.hot_sepdist
            )
            hcomp_req_supply = (
                tfcoil_variables.n_tf_coils
                * (cost_variables.life_plant / cost_variables.life_plant)
            ) * buildings_variables.qnty_sfty_fac
            # total storage space for required supply of outboard wall-blanket-shield
            ob_hotcell_vol = hcomp_req_supply * hcomp_vol
        else:
            ib_hotcell_vol = 0.0e0
            ob_hotcell_vol = 0.0e0

        # Divertor
        # Note: this estimation developed before the divertor design has been finalised
        if cost_variables.life_div_fpy != 0.0e0:
            hcomp_height = divertor_variables.dz_divertor
            hcomp_rad_thk = 2 * physics_variables.rminor
            hcomp_tor_thk = physics_variables.rmajor + physics_variables.rminor
            hcomp_footprint = (hcomp_height + buildings_variables.hot_sepdist) * (
                max(hcomp_rad_thk, hcomp_tor_thk) + buildings_variables.hot_sepdist
            )
            hcomp_vol = hcomp_footprint * (
                min(hcomp_rad_thk, hcomp_tor_thk) + buildings_variables.hot_sepdist
            )
            hcomp_req_supply = (
                tfcoil_variables.n_tf_coils
                * (cost_variables.life_plant / cost_variables.life_div_fpy)
            ) * buildings_variables.qnty_sfty_fac
            # total storage space for required supply of divertor segments
            div_hotcell_vol = hcomp_req_supply * hcomp_vol
        else:
            div_hotcell_vol = 0.0e0

        # Centre post
        if cost_variables.cplife != 0.0e0:
            hcomp_height = 2 * build_variables.z_tf_inside_half
            if tfcoil_variables.i_tf_sup != 1:
                hcomp_rad_thk = build_variables.r_cp_top
            else:
                hcomp_rad_thk = build_variables.dr_tf_inboard

            hcomp_footprint = (hcomp_height + buildings_variables.hot_sepdist) * (
                hcomp_rad_thk + buildings_variables.hot_sepdist
            )
            hcomp_vol = hcomp_footprint * (
                hcomp_rad_thk + buildings_variables.hot_sepdist
            )
            hcomp_req_supply = (
                cost_variables.life_plant / cost_variables.cplife
            ) * buildings_variables.qnty_sfty_fac
            # total storage space for required supply of centre posts
            cp_hotcell_vol = hcomp_req_supply * hcomp_vol
        else:
            cp_hotcell_vol = 0.0e0

        # building required internal volume and footprint
        hotcell_vol = ib_hotcell_vol + ob_hotcell_vol + div_hotcell_vol + cp_hotcell_vol
        # assumed building height based on R Gowland's estimates
        hotcell_area = hotcell_vol / buildings_variables.hotcell_h

        # derive estimates for length and width by assuming a square building
        hotcell_l = hotcell_area**0.5
        hotcell_w = hotcell_l

        # external dimensions include same wall and roof thicknesses as reactor building
        hotcell_area_ext = (hotcell_l + 2.0e0 * buildings_variables.reactor_wall_thk) * (
            hotcell_w + 2.0e0 * buildings_variables.reactor_wall_thk
        )
        hotcell_vol_ext = hotcell_area_ext * (
            buildings_variables.hotcell_h
            + buildings_variables.reactor_roof_thk
            + buildings_variables.reactor_fndtn_thk
        )

        buildings_total_vol = buildings_total_vol + hotcell_vol

        # Reactor Auxiliary Buildings
        # Derived from W. Smith's estimates of necessary facilities and their sizes;
        # these values amalgamate multiple individual buildings.

        # Chemistry labs: includes RA, non-RA and environmental labs,
        # and chemical treatment facilities for coolant circuits
        chemlab_area = buildings_variables.chemlab_l * buildings_variables.chemlab_w
        chemlab_vol = chemlab_area * buildings_variables.chemlab_h

        # Heat sink facilities, includes aux heat sink at heat energy island,
        # low temp and emergency heat sink facilities, ultimate heat sink facility
        # to sea/river/cooling towers, including pumping, chemical dosing and heat exchangers
        heat_sink_area = (
            buildings_variables.heat_sink_l * buildings_variables.heat_sink_w
        )
        heat_sink_vol = heat_sink_area * buildings_variables.heat_sink_h

        # auxiliary buildings supporting tokamak processes & systems, includes non-RA
        # interfacing services such as, hydraulics, compressed air, chilled water...
        aux_build_area = (
            buildings_variables.aux_build_l * buildings_variables.aux_build_w
        )
        aux_build_vol = aux_build_area * buildings_variables.aux_build_h

        # Total auxiliary buildings supporting reactor processes & systems
        reactor_aux_area = chemlab_area + heat_sink_area + aux_build_area
        reactor_aux_vol = chemlab_vol + heat_sink_vol + aux_build_vol

        buildings_total_vol = buildings_total_vol + reactor_aux_vol

        # Magnet power facilities
        # Providing specific electrical supplies for reactor magnets;
        # based upon dimensions of comparable equipment at ITER site.
        # Steady state power trains:
        magnet_trains_area = (
            buildings_variables.magnet_trains_l * buildings_variables.magnet_trains_w
        )
        magnet_trains_vol = magnet_trains_area * buildings_variables.magnet_trains_h
        # Pulsed power for central solenoid
        magnet_pulse_area = (
            buildings_variables.magnet_pulse_l * buildings_variables.magnet_pulse_w
        )
        magnet_pulse_vol = magnet_pulse_area * buildings_variables.magnet_pulse_h

        # Total power buildings areas and volumes
        power_buildings_area = hcd_building_area + magnet_trains_area + magnet_pulse_area
        power_buildings_vol = hcd_building_vol + magnet_trains_vol + magnet_pulse_vol

        buildings_total_vol = buildings_total_vol + power_buildings_vol

        # Control
        # Derived from W. Smith's estimates of necessary facilities and their sizes:
        # includes Main Control Room, Back-up Control Room,
        # Signal Processing and Distribution Centres [Safety Train A, Safety Train B],
        # HP offices & Data Logging centre, Data Storage centre;
        # these values amalgamate multiple individual buildings.
        control_buildings_area = (
            buildings_variables.control_buildings_l
            * buildings_variables.control_buildings_w
        )
        control_buildings_vol = (
            control_buildings_area * buildings_variables.control_buildings_h
        )

        buildings_total_vol = buildings_total_vol + control_buildings_vol

        # Warm Shop
        # Values taken from W. Smith's estimates of necessary facility size:
        # 'hands on maintenance workshops for low RA dose equipment'
        warm_shop_area = (
            buildings_variables.warm_shop_l * buildings_variables.warm_shop_w
        )
        warm_shop_vol = warm_shop_area * buildings_variables.warm_shop_h

        buildings_total_vol = buildings_total_vol + warm_shop_vol

        # Maintenance
        # Derived from W. Smith's estimates of necessary facilities and their sizes;
        # these values amalgamate multiple individual buildings.

        # Maintenance workshops and clean rooms for components with *no* radiation
        # inventory; should include allowance for overhead gantry and crane access
        workshop_area = buildings_variables.workshop_l * buildings_variables.workshop_w
        workshop_vol = workshop_area * buildings_variables.workshop_h

        # Robot construction, testing, mock-up facilities
        # To allow robots to be fully assembled, commissioned and tested
        # in mock-ups of the real environment. Height should allow for mock-up of
        # central column, but building also houses offices and classrooms.
        robotics_area = buildings_variables.robotics_l * buildings_variables.robotics_w
        robotics_vol = robotics_area * buildings_variables.robotics_h

        # Maintenance control and inspection facilities: includes operations centre,
        # inbound inspection and QA storage facilities.
        maint_cont_area = (
            buildings_variables.maint_cont_l * buildings_variables.maint_cont_w
        )
        maint_cont_vol = maint_cont_area * buildings_variables.maint_cont_h

        maintenance_area = workshop_area + robotics_area + maint_cont_area
        maintenance_vol = workshop_vol + robotics_vol + maint_cont_vol

        buildings_total_vol = buildings_total_vol + maintenance_vol

        # Cryogenic & cooling facilities
        # Derived from W. Smith's estimates of necessary facilities and their sizes.

        # Cryogenic Buildings for Magnet and Fuel Cycle
        cryomag_area = buildings_variables.cryomag_l * buildings_variables.cryomag_w
        cryomag_vol = cryomag_area * buildings_variables.cryomag_h

        # Magnet Cryo Storage Tanks
        cryostore_area = (
            buildings_variables.cryostore_l * buildings_variables.cryostore_w
        )
        cryostore_vol = cryostore_area * buildings_variables.cryostore_h

        # Site-Wide Auxiliary Cooling Water facility, including pumping,
        # chemical dosing, filtration and heat exchangers.
        auxcool_area = buildings_variables.auxcool_l * buildings_variables.auxcool_w
        auxcool_vol = auxcool_area * buildings_variables.auxcool_h

        cryocool_area = cryomag_area + cryostore_area + auxcool_area
        cryocool_vol = cryomag_vol + cryostore_vol + auxcool_vol

        buildings_total_vol = buildings_total_vol + cryocool_vol

        # Electrical
        # Derived from W. Smith's estimates of necessary facilities and their sizes;
        # these values amalgamate multiple individual buildings.

        # Transformers and electrical distribution facilities; includes
        # main step down & step up transformers and substation, reactive power buildings
        elecdist_area = buildings_variables.elecdist_l * buildings_variables.elecdist_w
        elecdist_vol = elecdist_area * buildings_variables.elecdist_h

        # Load centres (essential and non-essential supplies)
        elecload_area = buildings_variables.elecload_l * buildings_variables.elecload_w
        elecload_vol = elecload_area * buildings_variables.elecload_h

        # Energy Storage Systems (batteries & flywheels) and back-up generators
        elecstore_area = (
            buildings_variables.elecstore_l * buildings_variables.elecstore_w
        )
        elecstore_vol = elecstore_area * buildings_variables.elecstore_h

        elec_buildings_area = elecdist_area + elecload_area + elecstore_area
        elec_buildings_vol = elecdist_vol + elecload_vol + elecstore_vol

        buildings_total_vol = buildings_total_vol + elec_buildings_vol

        # Turbine Hall
        # As proposed by R. Gowland, based on assessment of 18 existing fission power plants:
        # turbine hall size is largely independent of plant output power.
        # The default footprint used here represents a weighted mean of those plants
        # and the design of a Steam Rankine cycle turbine building,
        # produced by Morsons as part of the Year 1 work.
        turbine_hall_area = (
            buildings_variables.turbine_hall_l * buildings_variables.turbine_hall_w
        )
        turbine_hall_vol = turbine_hall_area * buildings_variables.turbine_hall_h

        buildings_total_vol = buildings_total_vol + turbine_hall_vol

        # Waste
        # Derived from W. Smith's estimates of necessary facilities and their sizes.

        # Intermediate Level Waste
        # Radioactive waste melt, separation and size reduction facility
        ilw_smelter_area = (
            buildings_variables.ilw_smelter_l * buildings_variables.ilw_smelter_w
        )
        ilw_smelter_vol = ilw_smelter_area * buildings_variables.ilw_smelter_h
        # ILW process and storage, amalgamated buildings
        ilw_storage_area = (
            buildings_variables.ilw_storage_l * buildings_variables.ilw_storage_w
        )
        ilw_storage_vol = ilw_storage_area * buildings_variables.ilw_storage_h

        # Low Level Waste process and storage, amalgamated buildings
        llw_storage_area = (
            buildings_variables.llw_storage_l * buildings_variables.llw_storage_w
        )
        llw_storage_vol = llw_storage_area * buildings_variables.llw_storage_h

        # Hazardous Waste process and storage, amalgamated buildings
        hw_storage_area = (
            buildings_variables.hw_storage_l * buildings_variables.hw_storage_w
        )
        hw_storage_vol = hw_storage_area * buildings_variables.hw_storage_h

        # Tritiated Waste Store
        tw_storage_area = (
            buildings_variables.tw_storage_l * buildings_variables.tw_storage_w
        )
        tw_storage_vol = tw_storage_area * buildings_variables.tw_storage_h

        # Total waste buildings areas and volumes
        waste_buildings_area = (
            ilw_smelter_area
            + ilw_storage_area
            + llw_storage_area
            + hw_storage_area
            + tw_storage_area
        )
        waste_buildings_vol = (
            ilw_smelter_vol
            + ilw_storage_vol
            + llw_storage_vol
            + hw_storage_vol
            + tw_storage_vol
        )

        buildings_total_vol = buildings_total_vol + waste_buildings_vol

        # Site Services
        # Derived from W. Smith's estimates of necessary facilities and their sizes;
        # buildings grouped by function.

        # Air & Gas supplies
        # Includes compressed air facility, common gas systems facility, bottled gas
        # storage compounds; these values amalgamate multiple individual buildings.
        gas_buildings_area = (
            buildings_variables.gas_buildings_l * buildings_variables.gas_buildings_w
        )
        gas_buildings_vol = gas_buildings_area * buildings_variables.gas_buildings_h

        # Water, Laundry & Drainage
        # Includes facilities for potable water, firewater, chilled water; PPE laundry &
        # Respiratory Protective Equipment cleaning; industrial drains & sewage
        # process and discharge; these values amalgamate multiple individual buildings.
        water_buildings_area = (
            buildings_variables.water_buildings_l * buildings_variables.water_buildings_w
        )
        water_buildings_vol = (
            water_buildings_area * buildings_variables.water_buildings_h
        )

        # Site Security & Safety
        # Includes Security Control Centre and Fire and Ambulance Garages;
        # these values amalgamate multiple individual buildings.
        sec_buildings_area = (
            buildings_variables.sec_buildings_l * buildings_variables.sec_buildings_w
        )
        sec_buildings_vol = sec_buildings_area * buildings_variables.sec_buildings_h

        buildings_total_vol = (
            buildings_total_vol
            + gas_buildings_vol
            + water_buildings_vol
            + sec_buildings_vol
        )

        # Staff Services
        # Derived from W. Smith's estimates of necessary facilities and their sizes;
        # includes main office buildings, contractor offices, staff restaurant and cafe,
        # staff induction and training facilities, main gate and reception, access control
        # and site pass office, occupational health centre.
        # Amalgamates estimates of floor area for all individual buildings, uses average height.
        staff_buildings_vol = (
            buildings_variables.staff_buildings_area
            * buildings_variables.staff_buildings_h
        )

        buildings_total_vol = buildings_total_vol + staff_buildings_vol

        # Calculate 'effective floor area for AC power module'
        # This is the total floor area (m2) across the site, allowing for multiple floors
        # within buildings by assuming an average storey height of 6m:
        buildings_variables.a_plant_floor_effective = buildings_total_vol / 6.0e0

        # Total volume of nuclear buildings
        buildings_variables.volnucb = reactor_build_totvol + hotcell_vol_ext

        # Output
        if output:
            po.oheadr(self.outfile, "Power Plant Buildings")
            po.ovarre(
                self.outfile,
                "Reactor hall (internal) footprint (m2)",
                "(reactor_hall_area)",
                reactor_hall_area,
            )
            po.ovarre(
                self.outfile,
                "Reactor hall (internal) volume (m3)",
                "(reactor_hall_vol)",
                reactor_hall_vol,
            )
            po.ovarre(
                self.outfile,
                "   Reactor hall length (m)",
                "(reactor_hall_l)",
                buildings_variables.reactor_hall_l,
            )
            po.ovarre(
                self.outfile,
                "   Reactor hall width (m)",
                "(reactor_hall_w)",
                buildings_variables.reactor_hall_w,
            )
            po.ovarre(
                self.outfile,
                "   Reactor hall height (m)",
                "(reactor_hall_h)",
                buildings_variables.reactor_hall_h,
            )
            if (current_drive_variables.i_hcd_primary == 5) or (
                current_drive_variables.i_hcd_primary == 8
            ):
                po.ocmmnt(
                    self.outfile,
                    "   NBI HCD facility included within reactor building:",
                )
                po.ovarre(
                    self.outfile,
                    "      NBI system length (m)",
                    "(nbi_sys_l)",
                    buildings_variables.nbi_sys_l,
                )
                po.ovarre(
                    self.outfile,
                    "      NBI system width (m)",
                    "(nbi_sys_w)",
                    buildings_variables.nbi_sys_w,
                )

            po.ovarre(
                self.outfile,
                "Reactor building external footprint (m2)",
                "(reactor_building_area)",
                reactor_building_area,
            )
            po.ovarre(
                self.outfile,
                "Reactor building external volume (m3)",
                "(reactor_building_vol)",
                reactor_building_vol,
            )
            po.ovarre(
                self.outfile,
                "   Reactor building length (m)",
                "(reactor_building_l)",
                reactor_building_l,
            )
            po.ovarre(
                self.outfile,
                "   Reactor building width (m)",
                "(reactor_building_w)",
                reactor_building_w,
            )
            po.ovarre(
                self.outfile,
                "   Reactor building height (m)",
                "(reactor_building_h)",
                reactor_building_h,
            )
            po.ovarre(
                self.outfile,
                "Reactor basement footprint (m2)",
                "(reactor_basement_area)",
                reactor_basement_area,
            )
            po.ovarre(
                self.outfile,
                "Reactor basement volume (m3)",
                "(reactor_basement_vol)",
                reactor_basement_vol,
            )
            po.ovarre(
                self.outfile,
                "Reactor building + basement volume (m3)",
                "(reactor_build_totvol)",
                reactor_build_totvol,
            )
            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Hot cell facility internal footprint (m2)",
                "(hotcell_area)",
                hotcell_area,
            )
            po.ovarre(
                self.outfile,
                "Hot cell facility internal volume (m3)",
                "(hotcell_vol)",
                hotcell_vol,
            )
            po.ovarre(
                self.outfile,
                "Hot cell facility external footprint (m2)",
                "(hotcell_area_ext)",
                hotcell_area_ext,
            )
            po.ovarre(
                self.outfile,
                "Hot cell facility external volume (m3)",
                "(hotcell_vol_ext)",
                hotcell_vol_ext,
            )
            po.oblnkl(self.outfile)
            if (current_drive_variables.i_hcd_primary != 5) and (
                current_drive_variables.i_hcd_primary != 8
            ):
                po.ovarre(
                    self.outfile,
                    "HCD (EC/EBW) building footprint (m2)",
                    "(hcd_building_area)",
                    hcd_building_area,
                )
                po.ovarre(
                    self.outfile,
                    "HCD (EC/EBW) building volume (m3)",
                    "(hcd_building_vol)",
                    hcd_building_vol,
                )
                if buildings_variables.i_bldgs_v == 1:
                    po.ovarre(
                        self.outfile,
                        "   HCD (EC/EBW) building length (m)",
                        "(hcd_building_l)",
                        buildings_variables.hcd_building_l,
                    )
                    po.ovarre(
                        self.outfile,
                        "   HCD (EC/EBW) building width (m)",
                        "(hcd_building_w)",
                        buildings_variables.hcd_building_w,
                    )
                    po.ovarre(
                        self.outfile,
                        "   HCD (EC/EBW) building height (m)",
                        "(hcd_building_h)",
                        buildings_variables.hcd_building_h,
                    )
                po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Turbine hall footprint (m2)",
                "(turbine_hall_area)",
                turbine_hall_area,
            )
            po.ovarre(
                self.outfile,
                "Turbine hall volume (m3)",
                "(turbine_hall_vol)",
                turbine_hall_vol,
            )
            if buildings_variables.i_bldgs_v == 1:
                po.ovarre(
                    self.outfile,
                    "   Turbine hall length (m)",
                    "(turbine_hall_l)",
                    buildings_variables.turbine_hall_l,
                )
                po.ovarre(
                    self.outfile,
                    "   Turbine hall width (m)",
                    "(turbine_hall_w)",
                    buildings_variables.turbine_hall_w,
                )
                po.ovarre(
                    self.outfile,
                    "   Turbine hall height (m)",
                    "(turbine_hall_h)",
                    buildings_variables.turbine_hall_h,
                )

            po.oblnkl(self.outfile)
            po.ovarre(
                self.outfile,
                "Effective floor area (m2)",
                "(a_plant_floor_effective)",
                buildings_variables.a_plant_floor_effective,
            )
            po.ovarre(
                self.outfile,
                "Total volume of nuclear buildings (m3)",
                "(volnucb)",
                buildings_variables.volnucb,
            )

            if buildings_variables.i_bldgs_v == 1:
                po.oblnkl(self.outfile)
                # verbose output of building sizes, areas and volumes
                po.ovarre(
                    self.outfile,
                    "Chemistry labs and facilities footprint (m2)",
                    "(chemlab_area)",
                    chemlab_area,
                )
                po.ovarre(
                    self.outfile,
                    "Chemistry labs and facilities volume (m3)",
                    "(chemlab_vol)",
                    chemlab_vol,
                )

                po.ovarre(
                    self.outfile,
                    "Reactor support buildings footprint (m2)",
                    "(aux_build_area)",
                    aux_build_area,
                )
                po.ovarre(
                    self.outfile,
                    "Reactor support buildings volume (m3)",
                    "(aux_build_vol)",
                    aux_build_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Heat sinks footprint (m2)",
                    "(heat_sink_area)",
                    heat_sink_area,
                )
                po.ovarre(
                    self.outfile,
                    "Heat sinks volume (m3)",
                    "(heat_sink_vol)",
                    heat_sink_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Reactor auxiliary buildings footprint (m2)",
                    "(reactor_aux_area)",
                    reactor_aux_area,
                )
                po.ovarre(
                    self.outfile,
                    "Reactor auxiliary buildings volume (m3)",
                    "(reactor_aux_vol)",
                    reactor_aux_vol,
                )

                po.ovarre(
                    self.outfile,
                    "Magnet trains footprint (m2)",
                    "(magnet_trains_area)",
                    magnet_trains_area,
                )
                po.ovarre(
                    self.outfile,
                    "Magnet trains volume (m3)",
                    "(magnet_trains_vol)",
                    magnet_trains_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Magnet pulse footprint (m2)",
                    "(magnet_pulse_area)",
                    magnet_pulse_area,
                )
                po.ovarre(
                    self.outfile,
                    "Magnet pulse volume (m3)",
                    "(magnet_pulse_vol)",
                    magnet_pulse_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Power buildings footprint (m2)",
                    "(power_buildings_area)",
                    power_buildings_area,
                )
                po.ovarre(
                    self.outfile,
                    "Power buildings volume (m3)",
                    "(power_buildings_vol)",
                    power_buildings_vol,
                )

                po.ovarre(
                    self.outfile,
                    "Control buildings area (m2)",
                    "(control_buildings_area)",
                    control_buildings_area,
                )
                po.ovarre(
                    self.outfile,
                    "Control buildings volume (m3)",
                    "(control_buildings_vol)",
                    control_buildings_vol,
                )

                po.ovarre(
                    self.outfile,
                    "Warm shop footprint (m2)",
                    "(warm_shop_area)",
                    warm_shop_area,
                )
                po.ovarre(
                    self.outfile,
                    "Warm shop volume (m3)",
                    "(warm_shop_vol)",
                    warm_shop_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Workshop footprint (m2)",
                    "(workshop_area)",
                    workshop_area,
                )
                po.ovarre(
                    self.outfile, "Workshop volume (m3)", "(workshop_vol)", workshop_vol
                )
                po.ovarre(
                    self.outfile,
                    "Robotics building footprint (m2)",
                    "(robotics_area)",
                    robotics_area,
                )
                po.ovarre(
                    self.outfile,
                    "Robotics building volume (m3)",
                    "(robotics_vol)",
                    robotics_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Maintenance control footprint (m2)",
                    "(maint_cont_area)",
                    maint_cont_area,
                )
                po.ovarre(
                    self.outfile,
                    "Maintenance control volume (m3)",
                    "(maint_cont_vol)",
                    maint_cont_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Maintenance buildings footprint (m2)",
                    "(maintenance_area)",
                    maintenance_area,
                )
                po.ovarre(
                    self.outfile,
                    "Maintenance buildings volume (m3)",
                    "(maintenance_vol)",
                    maintenance_vol,
                )

                po.ovarre(
                    self.outfile,
                    "Cryogenic buildings footprint (m2)",
                    "(cryomag_area)",
                    cryomag_area,
                )
                po.ovarre(
                    self.outfile,
                    "Cryogenic buildings volume (m3)",
                    "(cryomag_vol)",
                    cryomag_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Magnet cryo storage tanks footprint (m2)",
                    "(cryostore_area)",
                    cryostore_area,
                )
                po.ovarre(
                    self.outfile,
                    "Magnet cryo storage tanks volume (m3)",
                    "(cryostore_vol)",
                    cryostore_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Auxiliary cooling footprint (m2)",
                    "(auxcool_area)",
                    auxcool_area,
                )
                po.ovarre(
                    self.outfile,
                    "Auxiliary cooling volume (m3)",
                    "(auxcool_vol)",
                    auxcool_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Cryogenic & cooling total footprint (m2)",
                    "(cryocool_area)",
                    cryocool_area,
                )
                po.ovarre(
                    self.outfile,
                    "Cryogenic & cooling total volume (m3)",
                    "(cryocool_vol)",
                    cryocool_vol,
                )

                po.ovarre(
                    self.outfile,
                    "Electrical transformers footprint (m2)",
                    "(elecdist_area)",
                    elecdist_area,
                )
                po.ovarre(
                    self.outfile,
                    "Electrical transformers volume (m3)",
                    "(elecdist_vol)",
                    elecdist_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Electrical load centres footprint (m2)",
                    "(elecload_area)",
                    elecload_area,
                )
                po.ovarre(
                    self.outfile,
                    "Electrical load centres volume (m3)",
                    "(elecload_vol)",
                    elecload_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Energy storage systems footprint (m2)",
                    "(elecstore_area)",
                    elecstore_area,
                )
                po.ovarre(
                    self.outfile,
                    "Energy storage systems volume (m3)",
                    "(elecstore_vol)",
                    elecstore_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Electrical buildings total footprint (m2)",
                    "(elec_buildings_area)",
                    elec_buildings_area,
                )
                po.ovarre(
                    self.outfile,
                    "Electrical buildings total volume (m3)",
                    "(elec_buildings_vol)",
                    elec_buildings_vol,
                )

                po.ovarre(
                    self.outfile,
                    "Waste buildings footprint (m2)",
                    "(waste_buildings_area)",
                    waste_buildings_area,
                )
                po.ovarre(
                    self.outfile,
                    "Waste buildings volume (m3)",
                    "(waste_buildings_vol)",
                    waste_buildings_vol,
                )

                po.ovarre(
                    self.outfile,
                    "Air & gas supplies footprint (m2)",
                    "(gas_buildings_area)",
                    gas_buildings_area,
                )
                po.ovarre(
                    self.outfile,
                    "Air & gas supplies volume (m3)",
                    "(gas_buildings_vol)",
                    gas_buildings_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Water supplies footprint (m2)",
                    "(water_buildings_area)",
                    water_buildings_area,
                )
                po.ovarre(
                    self.outfile,
                    "Water supplies volume (m3)",
                    "(water_buildings_vol)",
                    water_buildings_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Security & Safety buildings footprint (m2)",
                    "(sec_buildings_area)",
                    sec_buildings_area,
                )
                po.ovarre(
                    self.outfile,
                    "Security & Safety buildings volume (m3)",
                    "(sec_buildings_vol)",
                    sec_buildings_vol,
                )
                po.ovarre(
                    self.outfile,
                    "Staff buildings footprint (m2)",
                    "(staff_buildings_area)",
                    buildings_variables.staff_buildings_area,
                )
                po.ovarre(
                    self.outfile,
                    "Staff buildings volume (m3)",
                    "(staff_buildings_vol)",
                    staff_buildings_vol,
                )
