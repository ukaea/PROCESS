import logging

import numpy as np

from process import process_output as po
from process.data_structure import cost_2015_variables, cost_variables
from process.fortran import (
    build_variables,
    constants,
    current_drive_variables,
    fwbs_variables,
    global_variables,
    heat_transport_variables,
    pf_power_variables,
    pfcoil_variables,
    physics_variables,
    tfcoil_variables,
)

logger = logging.getLogger(__name__)


class Costs2015:
    def __init__(self):
        self.outfile = constants.nout

    def run(self):
        """
        Cost accounting for a fusion power plant
        author: J Morris, CCFE, Culham Science Centre
        outfile : input integer : output file unit
        iprint : input integer : switch for writing to output file (1=yes)
        This routine performs the cost accounting for a fusion power plant.
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        self.outfile = self.outfile

        # ###############################################

        # Calculate building costs
        self.calc_building_costs()

        # Calculate land costs
        self.calc_land_costs()

        # Calculate tf coil costs
        self.calc_tf_coil_costs()

        # Calculate fwbs costs
        self.calc_fwbs_costs()

        # Calculate remote handling costs
        self.calc_remote_handling_costs()

        # Calculate N plant and vacuum vessel costs
        self.calc_n_plant_and_vv_costs()

        # Calculate energy conversion system costs
        self.calc_energy_conversion_system()

        # Calculate remaining subsystems costs
        self.calc_remaining_subsystems()

        # Calculate total capital cost
        cost_2015_variables.total_costs = (
            cost_2015_variables.s_cost[8]
            + cost_2015_variables.s_cost[12]
            + cost_2015_variables.s_cost[20]
            + cost_2015_variables.s_cost[26]
            + cost_2015_variables.s_cost[30]
            + cost_2015_variables.s_cost[33]
            + cost_2015_variables.s_cost[34]
            + cost_2015_variables.s_cost[60]
        )

        # Save as concost, the variable used as a Figure of Merit (M$)
        cost_variables.concost = cost_2015_variables.total_costs / 1.0e6

        # Electrical output (given availability) for a whole year
        cost_2015_variables.mean_electric_output = (
            heat_transport_variables.p_plant_electric_net_mw * cost_variables.cpfact
        )
        cost_2015_variables.annual_electric_output = (
            cost_2015_variables.mean_electric_output * 24.0e0 * 365.25e0
        )

        # Annual maintenance cost.
        cost_2015_variables.maintenance = (
            cost_2015_variables.s_cost[26] + cost_2015_variables.s_cost[37]
        ) * cost_variables.maintenance_fwbs + (
            cost_2015_variables.s_cost[8]
            + cost_2015_variables.s_cost[30]
            + cost_2015_variables.s_cost[33]
            + cost_2015_variables.s_cost[34]
            + cost_2015_variables.s_cost[40]
            + cost_2015_variables.s_cost[42]
            + cost_2015_variables.s_cost[44]
            + cost_2015_variables.s_cost[46]
            + cost_2015_variables.s_cost[47]
            + cost_2015_variables.s_cost[48]
            + cost_2015_variables.s_cost[49]
            + cost_2015_variables.s_cost[50]
            + cost_2015_variables.s_cost[51]
            + cost_2015_variables.s_cost[52]
            + cost_2015_variables.s_cost[53]
            + cost_2015_variables.s_cost[57]
        ) * cost_variables.maintenance_gen

        # Levelized cost of electricity (LCOE) ($/MWh)
        if cost_2015_variables.annual_electric_output > 0.00001:
            cost_variables.coe = (
                1.0e0 / cost_2015_variables.annual_electric_output
            ) * (
                cost_2015_variables.total_costs / cost_variables.amortization
                + cost_2015_variables.maintenance
            )

        # Switch on output if there is a NaN error
        if (abs(cost_variables.concost) > 9.99e99) or (
            cost_variables.concost != cost_variables.concost
        ):
            self.output()

            for i in range(100):
                nan_diags = [
                    cost_2015_variables.s_label[i],
                    cost_2015_variables.s_kref[i],
                    cost_2015_variables.s_k[i],
                    cost_2015_variables.s_cref[i],
                    cost_2015_variables.s_cost[i],
                    cost_2015_variables.s_cost_factor[i],
                ]

                nan_diags_str = ",".join(str(x) for x in nan_diags)

                logger.info(nan_diags_str)
                po.ocmmnt(self.outfile, nan_diags_str)

            return

    def calc_fwbs_costs(self):
        """
        Function to calculate the cost of the first wall, blanket and shield
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine calculates the cost of the first wall, blanket and shield
        coils for a fusion power plant based on the costings in the PROCESS costs paper.
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """

        for i in range(21, 27):
            cost_2015_variables.s_cost_factor[i] = cost_variables.cost_factor_fwbs

        # Enrichment
        # Costs based on the number of separative work units (SWU) required
        #
        # SWU = P V(x_p) + T V(x_t) - F V(x_f)
        #
        # where V(x) is the value function
        #
        # V(x) = (1 - 2x)ln((1-x)/x)

        # Percentage of lithium 6 in the feed (natural abundance)
        feed_li6 = 0.0742e0
        # Percentage of lithium 6 in the tail (waste) (75% natural abundance)
        tail_li6 = feed_li6 * 0.75e0

        # Built-in test
        if global_variables.run_tests == 1:
            product_li6 = 0.3
            feed_to_product_mass_ratio = (product_li6 - tail_li6) / (
                feed_li6 - tail_li6
            )
            tail_to_product_mass_ratio = (product_li6 - feed_li6) / (
                feed_li6 - tail_li6
            )
            p_v = self.value_function(product_li6)
            t_v = self.value_function(tail_li6)
            f_v = self.value_function(feed_li6)
            swu = (
                p_v
                + tail_to_product_mass_ratio * t_v
                - feed_to_product_mass_ratio * f_v
            )
            if abs(swu - 2.66e0) < 2.0e-2:
                po.ocmmnt(
                    self.outfile,
                    "SWU for default 30% enrichment.  Should = 2.66. CORRECT",
                )
            else:
                po.ocmmnt(
                    self.outfile,
                    "SWU for default 30% enrichment.  Should = 2.66. ERROR",
                )

            # Reference cost
            cost_2015_variables.s_label[21] = "Lithium enrichment"
            cost_2015_variables.s_cref[21] = 0.1e6
            cost_2015_variables.s_k[21] = 64.7e0
            cost_2015_variables.s_kref[21] = 64.7e0
            cost_2015_variables.s_cost[21] = (
                cost_2015_variables.s_cost_factor[21]
                * cost_2015_variables.s_cref[21]
                * (cost_2015_variables.s_k[21] / cost_2015_variables.s_kref[21])
                ** cost_variables.costexp
            )
            if abs(cost_2015_variables.s_cost[21] - 0.1e6) / 0.1e6 < 1.0e-3:
                po.ocmmnt(self.outfile, "Reference cost for enrichment CORRECT")
            else:
                po.ocmmnt(self.outfile, "Reference cost for enrichment ERROR")

        # Lithium 6 enrichment cost ($)
        cost_2015_variables.s_label[21] = "Lithium enrichment"

        # Zero cost for natural enrichment
        if fwbs_variables.f_blkt_li6_enrichment <= 7.42e0:
            cost_2015_variables.s_cost[21] = 0.0e0
        else:
            # Percentage of lithium 6 in the product
            product_li6 = min(fwbs_variables.f_blkt_li6_enrichment, 99.99e0) / 100.0e0
            # SWU will be calculated for a unit mass of product (P=1)

            # Feed to product mass ratio
            feed_to_product_mass_ratio = (product_li6 - tail_li6) / (
                feed_li6 - tail_li6
            )

            # Tail to product mass ratio
            tail_to_product_mass_ratio = (product_li6 - feed_li6) / (
                feed_li6 - tail_li6
            )

            # Calculate value functions
            p_v = self.value_function(product_li6)
            t_v = self.value_function(tail_li6)
            f_v = self.value_function(feed_li6)

            # Calculate separative work units per kg
            swu = (
                p_v
                + tail_to_product_mass_ratio * t_v
                - feed_to_product_mass_ratio * f_v
            )

            # Mass of lithium (kg).  Lithium orthosilicate is 22% lithium by mass.
            mass_li = fwbs_variables.m_blkt_li2o * 0.22

            # Total swu for lithium in blanket
            total_swu = swu * mass_li

            # Reference cost for lithium enrichment (2014 $)
            cost_2015_variables.s_cref[21] = 0.1e6
            # Reference case of lithium SWU
            cost_2015_variables.s_k[21] = total_swu
            cost_2015_variables.s_kref[21] = 64.7e0
            cost_2015_variables.s_cost[21] = (
                cost_2015_variables.s_cost_factor[21]
                * cost_2015_variables.s_cref[21]
                * (cost_2015_variables.s_k[21] / cost_2015_variables.s_kref[21])
                ** cost_variables.costexp
            )

        cost_2015_variables.s_label[22] = "Lithium orthosilicate pebble manufacturing"
        # Reference cost of lithium pebble manufacture (2014 $)
        cost_2015_variables.s_cref[22] = 6.5e4
        # Scale with mass of pebbles (kg)
        cost_2015_variables.s_k[22] = fwbs_variables.m_blkt_li2o
        cost_2015_variables.s_kref[22] = 10.0e0
        cost_2015_variables.s_cost[22] = (
            cost_2015_variables.s_cost_factor[22]
            * cost_2015_variables.s_cref[22]
            * (cost_2015_variables.s_k[22] / cost_2015_variables.s_kref[22])
            ** cost_variables.costexp_pebbles
        )

        cost_2015_variables.s_label[23] = "Titanium beryllide pebble manufacturing"
        #  Reference cost of titanium beryllide pebble manufacture (2014 $)
        cost_2015_variables.s_cref[23] = 450.0e6
        #  Scale with mass of titanium beryllide pebbles (kg)
        cost_2015_variables.s_k[23] = fwbs_variables.m_blkt_beryllium
        cost_2015_variables.s_kref[23] = 1.0e5
        cost_2015_variables.s_cost[23] = (
            cost_2015_variables.s_cost_factor[23]
            * cost_2015_variables.s_cref[23]
            * (cost_2015_variables.s_k[23] / cost_2015_variables.s_kref[23])
            ** cost_variables.costexp_pebbles
        )

        cost_2015_variables.s_label[24] = "First wall W coating manufacturing"
        #  Reference (PPCS A) first wall W coating cost (2014 $)
        cost_2015_variables.s_cref[24] = 25.0e6
        #  First wall W coating mass (kg)
        cost_2015_variables.s_k[24] = (
            build_variables.a_fw_total
            * fwbs_variables.fw_armour_thickness
            * constants.den_tungsten
        )
        cost_2015_variables.s_kref[24] = 29000.0e0
        cost_2015_variables.s_cost[24] = (
            cost_2015_variables.s_cost_factor[24]
            * cost_2015_variables.s_cref[24]
            * (cost_2015_variables.s_k[24] / cost_2015_variables.s_kref[24])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[25] = (
            "Blanket and shield materials and manufacturing"
        )
        # The cost of making the blanket was estimated for PPCS A.
        # This cost includes only manufacturing - not R&D, transport, or assembly in the reactor.
        # It includes the first wall, blanket and shield, but excludes the breeder and multiplier materials.
        cost_2015_variables.s_cref[25] = 317.0e6
        #  Scale with steel mass in blanket + shield mass
        cost_2015_variables.s_k[25] = (
            fwbs_variables.m_blkt_steel_total + fwbs_variables.whtshld
        )
        cost_2015_variables.s_kref[25] = 4.07e6
        cost_2015_variables.s_cost[25] = (
            cost_2015_variables.s_cost_factor[25]
            * cost_2015_variables.s_cref[25]
            * (cost_2015_variables.s_k[25] / cost_2015_variables.s_kref[25])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[26] = "Total first wall and blanket cost"
        cost_2015_variables.s_cost[26] = 0.0e0
        for j in range(21, 26):
            cost_2015_variables.s_cost[26] = (
                cost_2015_variables.s_cost[26] + cost_2015_variables.s_cost[j]
            )

    def output(self):
        """
        Function to output the costs calculations
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine outputs the costs to output file
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        po.oheadr(
            self.outfile,
            'Estimate of "overnight" capital cost for a first of kind power plant (2014 M$)',
        )

        po.oshead(self.outfile, "Buildings (M$)")
        for i in range(9):
            self.ocost(
                self.outfile,
                cost_2015_variables.s_label[i],
                i + 1,
                cost_2015_variables.s_cost[i] / 1.0e6,
            )

        po.oshead(self.outfile, "Land (M$)")

        for j in range(9, 13):
            self.ocost(
                self.outfile,
                cost_2015_variables.s_label[j],
                j + 1,
                cost_2015_variables.s_cost[j] / 1.0e6,
            )

        po.oshead(self.outfile, "TF Coils (M$)")

        for k in range(13, 21):
            self.ocost(
                self.outfile,
                cost_2015_variables.s_label[k],
                k + 1,
                cost_2015_variables.s_cost[k] / 1.0e6,
            )

        po.oshead(self.outfile, "First wall and blanket (M$)")
        for l in range(21, 27):  # noqa: E741
            self.ocost(
                self.outfile,
                cost_2015_variables.s_label[l],
                l + 1,
                cost_2015_variables.s_cost[l] / 1.0e6,
            )

        po.oshead(self.outfile, "Active maintenance and remote handling (M$)")
        self.ocost(
            self.outfile,
            cost_2015_variables.s_label[27],
            28,
            cost_2015_variables.s_cost[27] / 1.0e6,
        )
        self.ocost(
            self.outfile,
            cost_2015_variables.s_label[28],
            29,
            cost_2015_variables.s_cost[28] / 1.0e6,
        )
        self.ocost(
            self.outfile,
            cost_2015_variables.s_label[30],
            31,
            cost_2015_variables.s_cost[30] / 1.0e6,
        )

        po.oshead(self.outfile, "Vacuum vessel and liquid nitrogen plant (M$)")
        for n in range(31, 34):
            self.ocost(
                self.outfile,
                cost_2015_variables.s_label[n],
                n + 1,
                cost_2015_variables.s_cost[n] / 1.0e6,
            )

        po.oshead(self.outfile, "System for converting heat to electricity (M$)")
        self.ocost(
            self.outfile,
            cost_2015_variables.s_label[34],
            35,
            cost_2015_variables.s_cost[34] / 1.0e6,
        )

        po.oshead(self.outfile, "Remaining subsystems (M$)")
        for q in range(35, 61):
            self.ocost(
                self.outfile,
                cost_2015_variables.s_label[q],
                q + 1,
                cost_2015_variables.s_cost[q] / 1.0e6,
            )

        po.oblnkl(self.outfile)
        self.ocost(
            self.outfile,
            "TOTAL OVERNIGHT CAPITAL COST (M$)",
            "(total_costs)",
            cost_2015_variables.total_costs / 1.0e6,
        )
        self.ocost(
            self.outfile,
            "Annual maintenance cost (M$)",
            "(maintenance)",
            cost_2015_variables.maintenance / 1.0e6,
        )
        po.oblnkl(self.outfile)
        po.ovarrf(
            self.outfile,
            "Net electric output (MW)",
            "(p_plant_electric_net_mw)",
            heat_transport_variables.p_plant_electric_net_mw,
            "OP ",
        )
        po.ovarrf(
            self.outfile, "Capacity factor", "(cpfact)", cost_variables.cpfact, "OP "
        )
        po.ovarrf(
            self.outfile,
            "Mean electric output (MW)",
            "(mean_electric_output)",
            cost_2015_variables.mean_electric_output,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Capital cost / mean electric output ($/W)",
            "",
            cost_2015_variables.total_costs
            / cost_2015_variables.mean_electric_output
            / 1.0e6,
            "OP ",
        )
        po.ovarrf(
            self.outfile,
            "Levelized cost of electricity ($/MWh)",
            "(coe)",
            cost_variables.coe,
            "OP ",
        )

    def calc_building_costs(self):
        """
        Function to calculate the cost of all buildings.
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine calculates the building costs for a fusion power plant
        based on the costings in the PROCESS costs Paper.
        Buildings have a different scaling law, with fixed cost per unit volume.
        Cref is therefore now f.Viter.unit_cost
        The costs for individual buildings must not be output,
        as the same mean cost per unit volume has been used both for light
        and for shielded buildings
        The exponent =1
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        for i in range(9):
            cost_2015_variables.s_cost_factor[i] = cost_variables.cost_factor_buildings

        # Power plant admin buildings cost ($)
        cost_2015_variables.s_label[0] = "Admin Buildings"
        cost_2015_variables.s_cref[0] = (
            129000.0e0 * cost_variables.light_build_cost_per_vol
        )
        cost_2015_variables.s_cost[0] = (
            cost_2015_variables.s_cost_factor[0] * cost_2015_variables.s_cref[0]
        )

        # Tokamak complex excluding hot cell cost ($)
        cost_2015_variables.s_label[1] = "Tokamak Complex (excluding hot cell)"
        cost_2015_variables.s_cref[1] = (
            1100000.0e0 * cost_variables.tok_build_cost_per_vol
        )
        # ITER cryostat volume (m^3)
        cost_2015_variables.s_k[1] = (
            (np.pi * fwbs_variables.r_cryostat_inboard**2)
            * 2.0e0
            * fwbs_variables.z_cryostat_half_inside
        )
        cost_2015_variables.s_kref[1] = 18712.0e0
        cost_2015_variables.s_cost[1] = (
            cost_2015_variables.s_cost_factor[1]
            * cost_2015_variables.s_cref[1]
            * (cost_2015_variables.s_k[1] / cost_2015_variables.s_kref[1])
        )

        # Neutral beam buildings cost ($)
        cost_2015_variables.s_label[2] = "Neutral beam buildings"
        cost_2015_variables.s_cref[2] = (
            28000.0e0 * cost_variables.light_build_cost_per_vol
        )
        # Scale with neutral beam wall plug power (MW)
        cost_2015_variables.s_k[2] = current_drive_variables.pwpnb
        cost_2015_variables.s_kref[2] = 120.0e0
        cost_2015_variables.s_cost[2] = (
            cost_2015_variables.s_cost_factor[2]
            * cost_2015_variables.s_cref[2]
            * (cost_2015_variables.s_k[2] / cost_2015_variables.s_kref[2])
        )

        # Cryoplant buildings cost ($)
        cost_2015_variables.s_label[3] = "Cryoplant buildings"
        cost_2015_variables.s_cref[3] = (
            130000.0e0 * cost_variables.light_build_cost_per_vol
        )
        # Scale with the total heat load on the cryoplant at ~4.5K (kW)
        cost_2015_variables.s_k[3] = heat_transport_variables.helpow / 1.0e3
        cost_2015_variables.s_kref[3] = 61.0e0
        cost_2015_variables.s_cost[3] = (
            cost_2015_variables.s_cost_factor[3]
            * cost_2015_variables.s_cref[3]
            * (cost_2015_variables.s_k[3] / cost_2015_variables.s_kref[3])
        )

        # PF Coil winding building cost ($)
        cost_2015_variables.s_label[4] = "PF Coil winding building"
        cost_2015_variables.s_cref[4] = (
            190000.0e0 * cost_variables.light_build_cost_per_vol
        )
        # Scale with the radius of the largest PF coil squared (m^2)
        cost_2015_variables.s_k[4] = pfcoil_variables.r_pf_coil_outer_max**2
        cost_2015_variables.s_kref[4] = 12.4e0**2
        cost_2015_variables.s_cost[4] = (
            cost_2015_variables.s_cost_factor[4]
            * cost_2015_variables.s_cref[4]
            * (cost_2015_variables.s_k[4] / cost_2015_variables.s_kref[4])
        )

        # Magnet power supplies and related buildings cost ($)
        cost_2015_variables.s_label[5] = "Magnet power supplies and related buildings"
        cost_2015_variables.s_cref[5] = (
            110000.0e0 * cost_variables.light_build_cost_per_vol
        )
        # Scale with TF current per coil (MA)
        cost_2015_variables.s_k[5] = (
            tfcoil_variables.c_tf_total / tfcoil_variables.n_tf_coils
        ) / 1.0e6
        cost_2015_variables.s_kref[5] = 9.1e0
        cost_2015_variables.s_cost[5] = (
            cost_2015_variables.s_cost_factor[5]
            * cost_2015_variables.s_cref[5]
            * (cost_2015_variables.s_k[5] / cost_2015_variables.s_kref[5])
        )

        # Magnet discharge buildings cost ($)
        cost_2015_variables.s_label[6] = "Magnet discharge buildings"
        cost_2015_variables.s_cref[6] = (
            35000.0e0 * cost_variables.light_build_cost_per_vol
        )
        # Scale with total stored energy in TF coils (GJ)
        cost_2015_variables.s_k[6] = tfcoil_variables.e_tf_magnetic_stored_total_gj
        cost_2015_variables.s_kref[6] = 41.0e0
        cost_2015_variables.s_cost[6] = (
            cost_2015_variables.s_cost_factor[6]
            * cost_2015_variables.s_cref[6]
            * (cost_2015_variables.s_k[6] / cost_2015_variables.s_kref[6])
        )

        # Heat removal system buildings cost ($)
        cost_2015_variables.s_label[7] = "Heat removal system buildings"
        # ITER volume of cooling water buildings (m^3)
        cost_2015_variables.s_cref[7] = (
            51000.0e0 * cost_variables.light_build_cost_per_vol
        )
        # Scale with total thermal power removed from the core (MW)
        cost_2015_variables.s_k[7] = (
            heat_transport_variables.p_plant_primary_heat_mw
            + heat_transport_variables.p_plant_secondary_heat_mw
        )
        cost_2015_variables.s_kref[7] = 880.0e0
        cost_2015_variables.s_cost[7] = (
            cost_2015_variables.s_cost_factor[7]
            * cost_2015_variables.s_cref[7]
            * (cost_2015_variables.s_k[7] / cost_2015_variables.s_kref[7])
        )

        # Total cost of buildings ($)
        cost_2015_variables.s_label[8] = "Total cost of buildings"
        cost_2015_variables.s_cost[8] = 0.0e0
        for j in range(8):
            cost_2015_variables.s_cost[8] = (
                cost_2015_variables.s_cost[8] + cost_2015_variables.s_cost[j]
            )

    def calc_land_costs(self):
        """
        Function to calculate the cost of land for the power plant
        author: J Morris, CCFE, Culham Science Centre
        None
        Land also uses a unit cost, but area is scaled.
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        for i in range(9, 13):
            cost_2015_variables.s_cost_factor[i] = cost_variables.cost_factor_land

        # Land purchasing cost ($)
        cost_2015_variables.s_label[9] = "Land purchasing"
        # ITER Land area (hectares)
        ITER_total_land_area = 180.0e0
        # ITER Land area for key buildings (hectares)
        ITER_key_buildings_land_area = 42.0e0
        # ITER buffer land (hectares)
        ITER_buffer_land_area = ITER_total_land_area - ITER_key_buildings_land_area

        # Scale with area of cryostat (m)
        cost_2015_variables.s_k[9] = np.pi * fwbs_variables.r_cryostat_inboard**2
        cost_2015_variables.s_kref[9] = 638.0e0
        # Cost of land per hectare (2014 $ / ha)
        cost_2015_variables.s_cref[9] = 318000.0e0
        # Cost of power plant land (2014 $)
        cost_2015_variables.s_cost[9] = (
            cost_2015_variables.s_cost_factor[9]
            * cost_2015_variables.s_cref[9]
            * (
                ITER_key_buildings_land_area
                * (cost_2015_variables.s_k[9] / cost_2015_variables.s_kref[9])
                ** cost_variables.costexp
                + ITER_buffer_land_area
            )
        )

        # Land improvement costs ($)
        cost_2015_variables.s_label[10] = "Land improvement"
        # Cost of clearing ITER land
        cost_2015_variables.s_cref[10] = 214.0e6
        # Scale with area of cryostat (m)
        cost_2015_variables.s_k[10] = np.pi * fwbs_variables.r_cryostat_inboard**2
        cost_2015_variables.s_kref[10] = 638.0e0
        cost_2015_variables.s_cost[10] = (
            cost_2015_variables.s_cost_factor[10]
            * (cost_2015_variables.s_k[10] / cost_2015_variables.s_kref[10])
            ** cost_variables.costexp
            * cost_2015_variables.s_cref[10]
        )

        # Road improvements cost ($)
        cost_2015_variables.s_label[11] = "Road improvements"
        # Cost of ITER road improvements
        cost_2015_variables.s_cref[11] = 150.0e6
        # Scale with TF coil longest dimension
        cost_2015_variables.s_k[11] = (
            max(build_variables.dh_tf_inner_bore, build_variables.dr_tf_inner_bore)
            + 2.0e0 * build_variables.dr_tf_inboard
        )
        cost_2015_variables.s_kref[11] = 14.0e0
        cost_2015_variables.s_cost[11] = (
            cost_2015_variables.s_cost_factor[11]
            * cost_2015_variables.s_cref[11]
            * (cost_2015_variables.s_k[11] / cost_2015_variables.s_kref[11])
            ** cost_variables.costexp
        )

        # Total land costs ($)
        cost_2015_variables.s_label[12] = "Total land costs"
        cost_2015_variables.s_cost[12] = 0.0e0
        for j in range(9, 12):
            cost_2015_variables.s_cost[12] = (
                cost_2015_variables.s_cost[12] + cost_2015_variables.s_cost[j]
            )

    def calc_tf_coil_costs(self):
        """
                Function to calculate the cost of the TF coils for the power plant
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine calculates the cost of the TF coils for a fusion power
        plant based on the costings in the PROCESS costs Paper.
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        for i in range(13, 20):
            cost_2015_variables.s_cost_factor[i] = cost_variables.cost_factor_tf_coils

        # TF coil insertion and welding costs ($)
        cost_2015_variables.s_label[13] = "TF Coil insertion and welding"
        # ITER coil insertion and welding cost (2014 $)
        cost_2015_variables.s_cref[13] = 258.0e6
        # Scale with total TF coil length (m)
        cost_2015_variables.s_k[13] = (
            tfcoil_variables.n_tf_coils * tfcoil_variables.len_tf_coil
        )
        cost_2015_variables.s_kref[13] = 18.0e0 * 34.1e0
        cost_2015_variables.s_cost[13] = (
            cost_2015_variables.s_cost_factor[13]
            * cost_2015_variables.s_cref[13]
            * (cost_2015_variables.s_k[13] / cost_2015_variables.s_kref[13])
            ** cost_variables.costexp
        )

        # TF coil winding costs ($)
        cost_2015_variables.s_label[15] = "TF coil winding"
        # ITER winding cost (2014 $)
        cost_2015_variables.s_cref[15] = 414.0e6
        # Scale with the total turn length (m)
        cost_2015_variables.s_k[15] = (
            tfcoil_variables.n_tf_coils
            * tfcoil_variables.len_tf_coil
            * tfcoil_variables.n_tf_coil_turns
        )
        cost_2015_variables.s_kref[15] = 82249.0e0
        cost_2015_variables.s_cost[15] = (
            cost_2015_variables.s_cost_factor[15]
            * cost_2015_variables.s_cref[15]
            * (cost_2015_variables.s_k[15] / cost_2015_variables.s_kref[15])
            ** cost_variables.costexp
        )

        # Copper stand cost for TF coil ($)
        cost_2015_variables.s_label[16] = "Copper strand for TF coil"
        # ITER Chromium plated Cu strand for TF SC cost (2014 $)
        cost_2015_variables.s_cref[16] = 21.0e6
        # Scale with total copper mass (kg)
        cost_2015_variables.s_k[16] = (
            tfcoil_variables.whtconcu * tfcoil_variables.n_tf_coils
        )
        cost_2015_variables.s_kref[16] = 244.0e3
        cost_2015_variables.s_cost[16] = (
            cost_2015_variables.s_cost_factor[16]
            * cost_2015_variables.s_cref[16]
            * (cost_2015_variables.s_k[16] / cost_2015_variables.s_kref[16])
            ** cost_variables.costexp
        )

        # superconductor strand cost ($)
        cost_2015_variables.s_label[17] = (
            "Strands with Nb3Sn superconductor and copper stabiliser"
        )
        # ITER Nb3Sn SC strands cost (2014 $)
        cost_2015_variables.s_cref[17] = 526.0e6
        # Scale with the total mass of Nb3Sn (kg)
        cost_2015_variables.s_k[17] = (
            tfcoil_variables.whtconsc * tfcoil_variables.n_tf_coils
        )
        cost_2015_variables.s_kref[17] = 210.0e3
        cost_2015_variables.s_cost[17] = (
            cost_2015_variables.s_cost_factor[17]
            * cost_2015_variables.s_cref[17]
            * (cost_2015_variables.s_k[17] / cost_2015_variables.s_kref[17])
            ** cost_variables.costexp
        )

        # Superconductor testing cost ($)
        cost_2015_variables.s_label[18] = "Testing of superconducting strands"
        # ITER Nb3Sn strand test costs (2014 $)
        cost_2015_variables.s_cref[18] = 4.0e6
        cost_2015_variables.s_cost[18] = (
            cost_2015_variables.s_cost_factor[18] * cost_2015_variables.s_cref[18]
        )

        # Superconductor cabling and jacketing cost ($)
        cost_2015_variables.s_label[19] = "Cabling and jacketing"
        # ITER cabling and jacketing costs (2014 $)
        cost_2015_variables.s_cref[19] = 81.0e6
        # Scale with total turn length.
        cost_2015_variables.s_k[19] = (
            tfcoil_variables.n_tf_coils
            * tfcoil_variables.len_tf_coil
            * tfcoil_variables.n_tf_coil_turns
        )
        cost_2015_variables.s_kref[19] = 82249.0e0
        cost_2015_variables.s_cost[19] = (
            cost_2015_variables.s_cost_factor[19]
            * cost_2015_variables.s_cref[19]
            * (cost_2015_variables.s_k[19] / cost_2015_variables.s_kref[19])
            ** cost_variables.costexp
        )

        # Total TF coil costs ($)
        cost_2015_variables.s_label[20] = "Total TF coil costs"
        cost_2015_variables.s_cost[20] = 0.0e0
        for j in range(13, 20):
            cost_2015_variables.s_cost[20] = (
                cost_2015_variables.s_cost[20] + cost_2015_variables.s_cost[j]
            )

    def calc_remote_handling_costs(self):
        """
        Function to calculate the cost of the remote handling facilities
        author: J Morris, CCFE, Culham Science Centre
        None
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        for i in range(27, 31):
            cost_2015_variables.s_cost_factor[i] = cost_variables.cost_factor_rh

        # K:\Power Plant Physics and Technology\Costs\Remote handling
        # From Sam Ha.

        cost_2015_variables.s_label[27] = "Moveable equipment"
        cost_2015_variables.s_cref[27] = 1.0e6 * (
            139.0e0 * cost_variables.num_rh_systems + 410.0e0
        )
        #  Scale with total mass of armour, first wall and blanket (kg)
        cost_2015_variables.s_kref[27] = 4.35e6
        cost_2015_variables.s_k[27] = fwbs_variables.armour_fw_bl_mass
        cost_2015_variables.s_cost[27] = (
            cost_2015_variables.s_cost_factor[27]
            * cost_2015_variables.s_cref[27]
            * (cost_2015_variables.s_k[27] / cost_2015_variables.s_kref[27])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[28] = (
            "Active maintenance facility with fixed equipment"
        )
        cost_2015_variables.s_cref[28] = 1.0e6 * (
            95.0e0 * cost_variables.num_rh_systems + 2562.0e0
        )
        #  Scale with total mass of armour, first wall and blanket (kg)
        cost_2015_variables.s_kref[28] = 4.35e6
        cost_2015_variables.s_k[28] = fwbs_variables.armour_fw_bl_mass
        cost_2015_variables.s_cost[28] = (
            cost_2015_variables.s_cost_factor[28]
            * cost_2015_variables.s_cref[28]
            * (cost_2015_variables.s_k[28] / cost_2015_variables.s_kref[28])
            ** cost_variables.costexp
        )

        # s(30) is not in use

        cost_2015_variables.s_label[30] = "Total remote handling costs"
        cost_2015_variables.s_cost[30] = (
            cost_2015_variables.s_cost[27] + cost_2015_variables.s_cost[28]
        )

    def calc_n_plant_and_vv_costs(self):
        """
        Function to calculate the cost of the nitrogen plant and vacuum vessel
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine calculates the cost of the nitrogen plant and vacuum vessel
        for a fusion power plant based on the costings in the PROCESS costs paper.
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        for i in range(31, 34):
            cost_2015_variables.s_cost_factor[i] = cost_variables.cost_factor_vv

        #  Vacuum vessel
        cost_2015_variables.s_label[31] = "Vacuum vessel"
        #  ITER reference vacuum vessel cost (2014 $)
        cost_2015_variables.s_cref[31] = 537.0e6
        #  Scale with outermost midplane radius of vacuum vessel squared (m2)
        cost_2015_variables.s_k[31] = (
            build_variables.rsldo + build_variables.dr_vv_outboard
        ) ** 2
        cost_2015_variables.s_kref[31] = 94.09e0
        cost_2015_variables.s_cost[31] = (
            cost_2015_variables.s_cost_factor[31]
            * cost_2015_variables.s_cref[31]
            * (cost_2015_variables.s_k[31] / cost_2015_variables.s_kref[31])
            ** cost_variables.costexp
        )

        #  Nitrogen plant
        cost_2015_variables.s_label[32] = "Liquid nitrogen plant"
        #  ITER reference cost (2014 $)
        cost_2015_variables.s_cref[32] = 86.0e6
        #  Scale with 4.5K cryopower (W)
        cost_2015_variables.s_k[32] = heat_transport_variables.helpow
        cost_2015_variables.s_kref[32] = 50.0e3
        cost_2015_variables.s_cost[32] = (
            cost_2015_variables.s_cost_factor[32]
            * cost_2015_variables.s_cref[32]
            * (cost_2015_variables.s_k[32] / cost_2015_variables.s_kref[32])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[33] = (
            "Total liquid nitrogen plant and vacuum vessel"
        )
        cost_2015_variables.s_cost[33] = 0.0e0
        for j in range(31, 33):
            cost_2015_variables.s_cost[33] = (
                cost_2015_variables.s_cost[33] + cost_2015_variables.s_cost[j]
            )

    def calc_energy_conversion_system(self):
        """
        Function to calculate the cost of the energy conversion system
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine calculates the cost of the energy conversion system
        for a fusion power plant based on the costings in the PROCESS costs paper.
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        cost_2015_variables.s_label[34] = "Energy conversion system"
        #  Set cost factor for energy conversion system
        cost_2015_variables.s_cost_factor[34] = cost_variables.cost_factor_bop
        #  Cost of reference energy conversion system (Rolls Royce)
        cost_2015_variables.s_cref[34] = 511.0e6
        #  Scale with gross electric power (MWe)
        cost_2015_variables.s_k[34] = heat_transport_variables.p_plant_electric_gross_mw
        cost_2015_variables.s_kref[34] = 692.0e0
        cost_2015_variables.s_cost[34] = (
            cost_2015_variables.s_cost_factor[34]
            * cost_2015_variables.s_cref[34]
            * (cost_2015_variables.s_k[34] / cost_2015_variables.s_kref[34])
            ** cost_variables.costexp
        )

    def calc_remaining_subsystems(self):
        """
        Function to calculate the cost of the remaining subsystems
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine calculates the cost of the remaining subsystems
        for a fusion power plant based on the costings in the PROCESS costs paper.
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        for i in range(35, 60):
            cost_2015_variables.s_cost_factor[i] = cost_variables.cost_factor_misc

        cost_2015_variables.s_label[35] = "CS and PF coils"
        # #  Cost of ITER CS and PF magnets
        cost_2015_variables.s_cref[35] = 1538.0e6
        #  Scale with sum of (A x turns x radius) of CS and all PF coils
        cost_2015_variables.s_k[35] = pfcoil_variables.itr_sum
        cost_2015_variables.s_kref[35] = 7.4e8
        cost_2015_variables.s_cost[35] = (
            cost_2015_variables.s_cost_factor[35]
            * cost_2015_variables.s_cref[35]
            * (cost_2015_variables.s_k[35] / cost_2015_variables.s_kref[35])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[36] = (
            "Vacuum vessel in-wall shielding, ports and in-vessel coils"
        )
        #  Cost of ITER VV in-wall shielding, ports and in-vessel coils
        cost_2015_variables.s_cref[36] = 211.0e6
        #  Scale with vacuum vessel mass (kg)
        cost_2015_variables.s_k[36] = fwbs_variables.m_vv
        cost_2015_variables.s_kref[36] = 5.2360e6
        cost_2015_variables.s_cost[36] = (
            cost_2015_variables.s_cost_factor[36]
            * cost_2015_variables.s_cref[36]
            * (cost_2015_variables.s_k[36] / cost_2015_variables.s_kref[36])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[37] = "Divertor"
        #  Cost of ITER divertor
        cost_2015_variables.s_cref[37] = 381.0e6
        #  Scale with max power to SOL (MW)
        cost_2015_variables.s_k[37] = physics_variables.p_plasma_separatrix_mw
        cost_2015_variables.s_kref[37] = 140.0e0
        cost_2015_variables.s_cost[37] = (
            cost_2015_variables.s_cost_factor[37]
            * cost_2015_variables.s_cref[37]
            * (cost_2015_variables.s_k[37] / cost_2015_variables.s_kref[37])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[38] = "not used"
        cost_2015_variables.s_label[39] = "not used"

        cost_2015_variables.s_label[40] = (
            "Ex-vessel neutral beam remote handling equipment"
        )
        #  Cost of ITER Ex-vessel NBI RH equipment
        # Increased to 90 Mdollar because of press release
        cost_2015_variables.s_cref[40] = 90.0e6
        #  Scale with total aux injected power (MW)
        cost_2015_variables.s_k[40] = current_drive_variables.p_hcd_injected_total_mw
        cost_2015_variables.s_kref[40] = 50.0e0
        cost_2015_variables.s_cost[40] = (
            cost_2015_variables.s_cost_factor[40]
            * cost_2015_variables.s_cref[40]
            * (cost_2015_variables.s_k[40] / cost_2015_variables.s_kref[40])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[41] = "not used"

        cost_2015_variables.s_label[42] = "Vacuum vessel pressure suppression system"
        #  Cost of ITER Vacuum vessel pressure suppression system
        cost_2015_variables.s_cref[42] = 40.0e6
        #  Scale with total thermal power removed from fusion core (MW)
        cost_2015_variables.s_k[42] = (
            heat_transport_variables.p_plant_primary_heat_mw
            + heat_transport_variables.p_plant_secondary_heat_mw
        )
        cost_2015_variables.s_kref[42] = 550.0e0
        cost_2015_variables.s_cost[42] = (
            cost_2015_variables.s_cost_factor[42]
            * cost_2015_variables.s_cref[42]
            * (cost_2015_variables.s_k[42] / cost_2015_variables.s_kref[42])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[43] = "Cryostat"
        #  Cost of ITER cryostat
        cost_2015_variables.s_cref[43] = 351.0e6
        #  Scale with cryostat external volume (m3)
        cost_2015_variables.s_k[43] = (
            (np.pi * fwbs_variables.r_cryostat_inboard**2.0e0)
            * 2.0e0
            * fwbs_variables.z_cryostat_half_inside
        )
        cost_2015_variables.s_kref[43] = 18700.0e0
        cost_2015_variables.s_cost[43] = (
            cost_2015_variables.s_cost_factor[43]
            * cost_2015_variables.s_cref[43]
            * (cost_2015_variables.s_k[43] / cost_2015_variables.s_kref[43])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[44] = "Heat removal system"
        #  Cost of ITER cooling water system
        cost_2015_variables.s_cref[44] = 724.0e6
        #  Scale with total thermal power removed from fusion core (MW)
        cost_2015_variables.s_k[44] = (
            heat_transport_variables.p_plant_primary_heat_mw
            + heat_transport_variables.p_plant_secondary_heat_mw
        )
        cost_2015_variables.s_kref[44] = 550.0e0
        cost_2015_variables.s_cost[44] = (
            cost_2015_variables.s_cost_factor[44]
            * cost_2015_variables.s_cref[44]
            * (cost_2015_variables.s_k[44] / cost_2015_variables.s_kref[44])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[45] = "Thermal shields"
        #  Cost of ITER thermal shields
        cost_2015_variables.s_cref[45] = 126.0e6
        #  Scale with cryostat surface area (m2)
        cost_2015_variables.s_k[45] = (
            2.0e0
            * np.pi
            * fwbs_variables.r_cryostat_inboard
            * 2.0e0
            * fwbs_variables.z_cryostat_half_inside
            + 2 * (np.pi * fwbs_variables.r_cryostat_inboard**2)
        )
        cost_2015_variables.s_kref[45] = 3902.0e0
        cost_2015_variables.s_cost[45] = (
            cost_2015_variables.s_cost_factor[45]
            * cost_2015_variables.s_cref[45]
            * (cost_2015_variables.s_k[45] / cost_2015_variables.s_kref[45])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[46] = "Pellet injection system"
        #  Cost of ITER pellet injector and pellet injection system
        cost_2015_variables.s_cref[46] = 25.0e6
        #  Scale with fusion power (MW)
        cost_2015_variables.s_k[46] = physics_variables.p_fusion_total_mw
        cost_2015_variables.s_kref[46] = 500.0e0
        cost_2015_variables.s_cost[46] = (
            cost_2015_variables.s_cost_factor[46]
            * cost_2015_variables.s_cref[46]
            * (cost_2015_variables.s_k[46] / cost_2015_variables.s_kref[46])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[47] = "Gas injection and wall conditioning system"
        # #  Cost of ITER gas injection system, GDC, Gi valve boxes
        cost_2015_variables.s_cref[47] = 32.0e6
        #  Scale with fusion power (MW)
        cost_2015_variables.s_k[47] = physics_variables.p_fusion_total_mw
        cost_2015_variables.s_kref[47] = 500.0e0
        cost_2015_variables.s_cost[47] = (
            cost_2015_variables.s_cost_factor[47]
            * cost_2015_variables.s_cref[47]
            * (cost_2015_variables.s_k[47] / cost_2015_variables.s_kref[47])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[48] = "Vacuum pumping"
        #  Cost of ITER vacuum pumping
        cost_2015_variables.s_cref[48] = 201.0e6
        #  Scale with fusion power (MW)
        cost_2015_variables.s_k[48] = physics_variables.p_fusion_total_mw
        cost_2015_variables.s_kref[48] = 500.0e0
        cost_2015_variables.s_cost[48] = (
            cost_2015_variables.s_cost_factor[48]
            * cost_2015_variables.s_cref[48]
            * (cost_2015_variables.s_k[48] / cost_2015_variables.s_kref[48])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[49] = "Tritium plant"
        #  Cost of ITER tritium plant
        cost_2015_variables.s_cref[49] = 226.0e6
        #  Scale with fusion power (MW)
        cost_2015_variables.s_k[49] = physics_variables.p_fusion_total_mw
        cost_2015_variables.s_kref[49] = 500.0e0
        cost_2015_variables.s_cost[49] = (
            cost_2015_variables.s_cost_factor[49]
            * cost_2015_variables.s_cref[49]
            * (cost_2015_variables.s_k[49] / cost_2015_variables.s_kref[49])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[50] = "Cryoplant and distribution"
        #  Cost of ITER Cryoplant and distribution
        cost_2015_variables.s_cref[50] = 397.0e6
        #  Scale with heat removal at 4.5 K approx (W)
        cost_2015_variables.s_k[50] = heat_transport_variables.helpow
        cost_2015_variables.s_kref[50] = 50000.0e0
        cost_2015_variables.s_cost[50] = (
            cost_2015_variables.s_cost_factor[50]
            * cost_2015_variables.s_cref[50]
            * (cost_2015_variables.s_k[50] / cost_2015_variables.s_kref[50])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[51] = "Electrical power supply and distribution"
        #  Cost of ITER electrical power supply and distribution
        cost_2015_variables.s_cref[51] = 1188.0e6
        #  Scale with total magnetic energy in the poloidal field / resistive diffusion time (W)
        #  For ITER value see
        #  K:\Power Plant Physics and Technology\PROCESS\PROCESS documentation papers\resistive diffusion time.xmcd or pdf
        cost_2015_variables.s_k[51] = (
            pf_power_variables.ensxpfm
            * 1.0e6
            / physics_variables.t_plasma_res_diffusion
        )
        cost_2015_variables.s_kref[51] = 8.0e9 / 953.0e0
        cost_2015_variables.s_cost[51] = (
            cost_2015_variables.s_cost_factor[51]
            * cost_2015_variables.s_cref[51]
            * (cost_2015_variables.s_k[51] / cost_2015_variables.s_kref[51])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[52] = (
            "Neutral beam heating and current drive system"
        )
        #  Cost of ITER NB H and CD
        cost_2015_variables.s_cref[52] = 814.0e6
        #  Scale with total auxiliary injected power (MW)
        cost_2015_variables.s_k[52] = current_drive_variables.p_hcd_injected_total_mw
        cost_2015_variables.s_kref[52] = 50.0e0
        cost_2015_variables.s_cost[52] = (
            cost_2015_variables.s_cost_factor[52]
            * cost_2015_variables.s_cref[52]
            * (cost_2015_variables.s_k[52] / cost_2015_variables.s_kref[52])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[53] = "Diagnostics systems"
        #  Cost of ITER diagnostic systems
        cost_2015_variables.s_cref[53] = 640.0e6
        # No scaling
        cost_2015_variables.s_cost[53] = (
            cost_2015_variables.s_cost_factor[53] * cost_2015_variables.s_cref[53]
        )

        cost_2015_variables.s_label[54] = "Radiological protection"
        #  Cost of ITER radiological protection
        cost_2015_variables.s_cref[54] = 19.0e6
        #  Scale with fusion power (MW)
        cost_2015_variables.s_k[54] = physics_variables.p_fusion_total_mw
        cost_2015_variables.s_kref[54] = 500.0e0
        cost_2015_variables.s_cost[54] = (
            cost_2015_variables.s_cost_factor[54]
            * cost_2015_variables.s_cref[54]
            * (cost_2015_variables.s_k[54] / cost_2015_variables.s_kref[54])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[55] = "Access control and security systems"
        #  Cost of ITER access control and security systems
        #  Scale with area of cryostat (m2)
        cost_2015_variables.s_k[55] = np.pi * fwbs_variables.r_cryostat_inboard**2
        cost_2015_variables.s_kref[55] = 640.0e0
        cost_2015_variables.s_cref[55] = 42.0e6
        cost_2015_variables.s_cost[55] = (
            cost_2015_variables.s_cost_factor[55]
            * cost_2015_variables.s_cref[55]
            * (cost_2015_variables.s_k[55] / cost_2015_variables.s_kref[55])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[56] = "Assembly"
        #  Cost of ITER assembly
        cost_2015_variables.s_cref[56] = 732.0e6
        #  Scale with total cost of reactor items (cryostat and everything inside it)
        cost_2015_variables.s_k[56] = (
            cost_2015_variables.s_cost[20]
            + cost_2015_variables.s_cost[26]
            + cost_2015_variables.s_cost[31]
            + cost_2015_variables.s_cost[35]
            + cost_2015_variables.s_cost[36]
            + cost_2015_variables.s_cost[37]
            + cost_2015_variables.s_cost[43]
            + cost_2015_variables.s_cost[45]
            + cost_2015_variables.s_cost[48]
        )
        cost_2015_variables.s_kref[56] = (
            cost_2015_variables.s_cref[20]
            + cost_2015_variables.s_cref[26]
            + cost_2015_variables.s_cref[31]
            + cost_2015_variables.s_cref[35]
            + cost_2015_variables.s_cref[36]
            + cost_2015_variables.s_cref[37]
            + cost_2015_variables.s_cref[43]
            + cost_2015_variables.s_cref[45]
            + cost_2015_variables.s_cref[48]
        )
        cost_2015_variables.s_cost[56] = (
            cost_2015_variables.s_cost_factor[56]
            * cost_2015_variables.s_cref[56]
            * (cost_2015_variables.s_k[56] / cost_2015_variables.s_kref[56])
        )

        cost_2015_variables.s_label[57] = "Control and communication"
        #  Cost of ITER control and data access and communication
        cost_2015_variables.s_cref[57] = 219.0e6
        #  Scale with total cost of reactor items (cryostat and everythign inside it)
        cost_2015_variables.s_k[57] = (
            cost_2015_variables.s_cost[20]
            + cost_2015_variables.s_cost[26]
            + cost_2015_variables.s_cost[31]
            + cost_2015_variables.s_cost[35]
            + cost_2015_variables.s_cost[36]
            + cost_2015_variables.s_cost[37]
            + cost_2015_variables.s_cost[43]
            + cost_2015_variables.s_cost[45]
            + cost_2015_variables.s_cost[48]
        )
        cost_2015_variables.s_kref[57] = (
            cost_2015_variables.s_cref[20]
            + cost_2015_variables.s_cref[26]
            + cost_2015_variables.s_cref[31]
            + cost_2015_variables.s_cref[35]
            + cost_2015_variables.s_cref[36]
            + cost_2015_variables.s_cref[37]
            + cost_2015_variables.s_cref[43]
            + cost_2015_variables.s_cref[45]
            + cost_2015_variables.s_cref[48]
        )
        cost_2015_variables.s_cost[57] = (
            cost_2015_variables.s_cost_factor[57]
            * cost_2015_variables.s_cref[57]
            * (cost_2015_variables.s_k[57] / cost_2015_variables.s_kref[57])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[58] = "Additional project expenditure"
        #  Cost of ITER additional ITER IO expenditure
        cost_2015_variables.s_cref[58] = 1624.0e6
        cost_2015_variables.s_cost[58] = (
            cost_2015_variables.s_cost_factor[58] * cost_2015_variables.s_cref[58]
        )

        # Calculate miscellaneous costs
        cost_2015_variables.s_label[59] = "Logistics"
        cost_2015_variables.s_cref[59] = 129.0e6
        #  Scale with cryostat external volume (m)
        cost_2015_variables.s_k[59] = (
            np.pi
            * fwbs_variables.r_cryostat_inboard**2
            * 2.0e0
            * fwbs_variables.z_cryostat_half_inside
        )
        cost_2015_variables.s_kref[59] = 18700.0e0
        cost_2015_variables.s_cost[59] = (
            cost_2015_variables.s_cost_factor[59]
            * cost_2015_variables.s_cref[59]
            * (cost_2015_variables.s_k[59] / cost_2015_variables.s_kref[59])
            ** cost_variables.costexp
        )

        cost_2015_variables.s_label[60] = "Total remaining subsystem costs"
        cost_2015_variables.s_cost[60] = 0.0e0
        for j in range(35, 60):
            cost_2015_variables.s_cost[60] = (
                cost_2015_variables.s_cost[60] + cost_2015_variables.s_cost[j]
            )

    def value_function(self, x):
        """
        Value function
        author: J Morris, CCFE, Culham Science Centre
        None
        Function for separative work unit calculation for enrichment cost
        PROCESS Costs Paper (M. Kovari, J. Morris)
        """
        return (1.0e0 - 2.0e0 * x) * np.log((1.0e0 - x) / x)

    def ocost(self, file, descr, vname, value):
        """
        Routine to print out the code, description and value
        of a cost item from array s in costs_2015
        """

        #  Local variables
        # character(len=70) :: dum70

        if descr == "not used":
            return

        # !TODO: Convert this

        # Replace descr with dummy string of the correct length.
        #       dum70 = descr
        #       write(file,10) dum70, value, ' '
        # 10    format(1x,a,t73,f10.0, tl1, a)

        # Create variable name of format s + array entry

        po.ovarrf(file, descr, vname, value)

    def ocost_vname(self, file, descr, vname, value):
        """
        Routine to print out the code, description and value
        of a cost item not in the array s in costs_2015
        """

        # character(len=70) :: dum70

        if descr == "not used":
            return

        # !TODO: Convert this

        # Replace descr with dummy string of the correct length.
        #       dum70 = descr
        #       write(file,10) dum70, value, ' '
        # 10    format(1x,a,t73,f10.0, tl1, a)

        po.ovarrf(file, descr, vname, value)
