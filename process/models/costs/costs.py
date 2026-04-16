import logging

import numpy as np

from process.core import constants
from process.core import process_output as po
from process.core.exceptions import ProcessValueError
from process.core.model import Model
from process.data_structure import (
    build_variables,
    buildings_variables,
    current_drive_variables,
    divertor_variables,
    fwbs_variables,
    heat_transport_variables,
    ife_variables,
    pf_power_variables,
    pfcoil_variables,
    physics_variables,
    pulse_variables,
    structure_variables,
    tfcoil_variables,
    times_variables,
)
from process.models.tfcoil.base import TFConductorModel

logger = logging.getLogger(__name__)


class Costs(Model):
    def __init__(self):
        self.outfile = constants.NOUT

    def run(self):
        """Cost accounting for a fusion power plant


        This routine performs the cost accounting for a fusion power plant.
        The direct costs are calculated based on parameters input
        from other sections of the code.
        <P>Costs are in 1990 $, and assume first-of-a-kind components
        unless otherwise stated. Account 22 costs include a multiplier
        to account for Nth-of-a-kind cost reductions.
        <P>The code is arranged in the order of the standard accounts.
        """
        # Convert FPY component lifetimes to calendar years
        # for replacement components
        self.convert_fpy_to_calendar()

        self.acc21()

        #  Account 22 : Fusion power island
        self.acc22()

        #  Account 23 : Turbine plant equipment
        self.acc23()

        #  Account 24 : Electric plant equipment
        self.acc241()  # Account 241 : Switchyard
        self.acc242()  # Account 242 : Transformers
        self.acc243()  # Account 243 : Low voltage
        self.acc244()  # Account 244 : Diesel generators
        self.acc245()  # Account 245 : Auxiliary facility power equipment
        self.acc24()  # Account 24  : Total

        #  Account 25 : Miscellaneous plant equipment
        self.acc25()

        #  Account 26 : Heat rejection system
        self.acc26()

        #  Total plant direct cost
        # cdirt = c21 + c22 + self.data.costs.c23 + self.data.costs.c24 + self.data.costs.c25 + self.data.costs.c26
        self.data.costs.cdirt = (
            self.data.costs.c21
            + self.data.costs.c22
            + self.data.costs.c23
            + self.data.costs.c24
            + self.data.costs.c25
            + self.data.costs.c26
        )

        #  Account 9 : Indirect cost and project contingency
        self.acc9()

        #  Constructed cost
        self.data.costs.concost = (
            self.data.costs.cdirt + self.data.costs.cindrt + self.data.costs.ccont
        )

        #  Cost of electricity
        if (self.data.costs.ireactor == 1) and (self.data.costs.ipnet == 0):
            self.coelc()

    def output(self):
        self.run()
        if self.data.costs.output_costs == 0:
            return

        po.oheadr(self.outfile, "Power Reactor Costs (1990 US$)")

        po.ovarrf(
            self.outfile,
            "First wall / blanket life (years)",
            "(life_blkt)",
            fwbs_variables.life_blkt,
        )

        if ife_variables.ife != 1:
            po.ovarrf(
                self.outfile,
                "Divertor life (years)",
                "(life_div)",
                self.data.costs.life_div,
            )
            if physics_variables.itart == 1:
                po.ovarrf(
                    self.outfile,
                    "Centrepost life (years)",
                    "(cplife_cal)",
                    self.data.costs.cplife_cal,
                )

        po.ovarrf(
            self.outfile, "Cost of electricity (m$/kWh)", "(coe)", self.data.costs.coe
        )

        po.osubhd(self.outfile, "Power Generation Costs :")
        # TODO: Convert fortran format to Python
        # if ((annfwbl != annfwbl) or (annfwbl > 1.0e10) or (annfwbl < 0.0e0)) :
        #     write(outfile,*)'Problem with annfwbl'
        #     write(outfile,*)'fwallcst=', fwallcst, '  blkcst=', self.data.costs.blkcst
        #     write(outfile,*)'crffwbl=', crffwbl,   '  fcap0cp=', self.data.costs.fcap0cp
        #     write(outfile,*)'feffwbl=', feffwbl,   '  fwbllife=', fwbllife

        #       write(outfile,200) #          anncap,coecap, #          annoam,coeoam, #          anndecom,coedecom, #          annfwbl,coefwbl, #          anndiv,coediv, #          anncp,coecp, #          anncdr,coecdr, #          annfuel,coefuel, #          annwst,coewst, #          annfuelt,coefuelt, #          anntot,coe

        # 200   format( #          t76,'Annual Costs, M$       COE, m$/kWh'// #          1x,'Capital Investment',t80,f10.2,10x,f10.2/ #          1x,'Operation & Maintenance',t80,f10.2,10x,f10.2/ #          1x,'Decommissioning Fund',t80,f10.2,10x,f10.2/ #          1x,'Fuel Charge Breakdown'// #          5x,'Blanket & first wall',t72,f10.2,10x,f10.2/ #          5x,'Divertors',t72,f10.2,10x,f10.2/ #          5x,'Centrepost (TART only)',t72,f10.2,10x,f10.2/ #          5x,'Auxiliary Heating',t72,f10.2,10x,f10.2/ #          5x,'Actual Fuel',t72,f10.2,10x,f10.2/ #          5x,'Waste Disposal',t72,f10.2,10x,f10.2/ #          1x,'Total Fuel Cost',t80,f10.2,10x,f10.2// #          1x,'Total Cost',t80,f10.2,10x,f10.2 )

        if self.data.costs.ifueltyp == 1:
            po.oshead(self.outfile, "Replaceable Components Direct Capital Cost")
            po.ovarrf(
                self.outfile,
                "First wall direct capital cost (M$)",
                "(fwallcst)",
                self.data.costs.fwallcst,
            )
            po.ovarrf(
                self.outfile,
                "Blanket direct capital cost (M$)",
                "(blkcst)",
                self.data.costs.blkcst,
            )
            if ife_variables.ife != 1:
                po.ovarrf(
                    self.outfile,
                    "Divertor direct capital cost (M$)",
                    "(divcst)",
                    self.data.costs.divcst,
                )
                if physics_variables.itart == 1:
                    po.ovarrf(
                        self.outfile,
                        "Centrepost direct capital cost (M$)",
                        "(cpstcst)",
                        self.data.costs.cpstcst,
                    )

                po.ovarrf(
                    self.outfile,
                    "Plasma heating/CD system cap cost (M$)",
                    "",
                    self.data.costs.cdcost
                    * self.data.costs.fcdfuel
                    / (1.0e0 - self.data.costs.fcdfuel),
                )
                po.ovarrf(
                    self.outfile,
                    "Fraction of CD cost --> fuel cost",
                    "(fcdfuel)",
                    self.data.costs.fcdfuel,
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "IFE driver system direct cap cost (M$)",
                    "",
                    self.data.costs.cdcost
                    * self.data.costs.fcdfuel
                    / (1.0e0 - self.data.costs.fcdfuel),
                )
                po.ovarrf(
                    self.outfile,
                    "Fraction of driver cost --> fuel cost",
                    "(fcdfuel)",
                    self.data.costs.fcdfuel,
                )

        po.oheadr(self.outfile, "Detailed Costings (1990 US$)")
        po.ovarre(
            self.outfile,
            "Acc.22 multiplier for Nth of a kind",
            "(fkind)",
            self.data.costs.fkind,
        )
        po.ovarin(
            self.outfile, "Level of Safety Assurance", "(lsa)", self.data.costs.lsa
        )
        po.oblnkl(self.outfile)
        po.oshead(self.outfile, "Structures and Site Facilities")
        po.ocosts(
            self.outfile,
            "(c211)",
            "Site improvements, facilities, land (M$)",
            self.data.costs.c211,
        )
        po.ocosts(
            self.outfile,
            "(c212)",
            "Reactor building cost (M$)",
            self.data.costs.c212,
        )
        po.ocosts(
            self.outfile,
            "(c213)",
            "Turbine building cost (M$)",
            self.data.costs.c213,
        )
        po.ocosts(
            self.outfile,
            "(c2141)",
            "Reactor maintenance building cost (M$)",
            self.data.costs.c2141,
        )
        po.ocosts(self.outfile, "(c2142)", "Warm shop cost (M$)", self.data.costs.c2142)
        po.ocosts(
            self.outfile,
            "(c215)",
            "Tritium building cost (M$)",
            self.data.costs.c215,
        )
        po.ocosts(
            self.outfile,
            "(c216)",
            "Electrical equipment building cost (M$)",
            self.data.costs.c216,
        )
        po.ocosts(
            self.outfile,
            "(c2171)",
            "Additional buildings cost (M$)",
            self.data.costs.c2171,
        )
        po.ocosts(
            self.outfile,
            "(c2172)",
            "Control room buildings cost (M$)",
            self.data.costs.c2172,
        )
        po.ocosts(
            self.outfile,
            "(c2173)",
            "Shop and warehouses cost (M$)",
            self.data.costs.c2173,
        )
        po.ocosts(
            self.outfile,
            "(c2174)",
            "Cryogenic building cost (M$)",
            self.data.costs.c2174,
        )
        po.oblnkl(self.outfile)
        po.ocosts(
            self.outfile,
            "(c21)",
            "Total account 21 cost (M$)",
            self.data.costs.c21,
        )

        po.oshead(self.outfile, "Reactor Systems")
        po.ocosts(self.outfile, "(c2211)", "First wall cost (M$)", self.data.costs.c2211)
        if ife_variables.ife != 1:
            po.ocosts(
                self.outfile,
                "(c22121)",
                "Blanket beryllium cost (M$)",
                self.data.costs.c22121,
            )
            po.ocosts(
                self.outfile,
                "(c22122)",
                "Blanket breeder material cost (M$)",
                self.data.costs.c22122,
            )

            po.ocosts(
                self.outfile,
                "(c22123)",
                "Blanket stainless steel cost (M$)",
                self.data.costs.c22123,
            )
            po.ocosts(
                self.outfile,
                "(c22124)",
                "Blanket vanadium cost (M$)",
                self.data.costs.c22124,
            )
        else:  # IFE
            po.ocosts(
                self.outfile,
                "(c22121)",
                "Blanket beryllium cost (M$)",
                self.data.costs.c22121,
            )
            po.ocosts(
                self.outfile,
                "(c22122)",
                "Blanket lithium oxide cost (M$)",
                self.data.costs.c22122,
            )
            po.ocosts(
                self.outfile,
                "(c22123)",
                "Blanket stainless steel cost (M$)",
                self.data.costs.c22123,
            )
            po.ocosts(
                self.outfile,
                "(c22124)",
                "Blanket vanadium cost (M$)",
                self.data.costs.c22124,
            )
            po.ocosts(
                self.outfile,
                "(c22125)",
                "Blanket carbon cloth cost (M$)",
                self.data.costs.c22125,
            )
            po.ocosts(
                self.outfile,
                "(c22126)",
                "Blanket concrete cost (M$)",
                self.data.costs.c22126,
            )
            po.ocosts(
                self.outfile,
                "(c22127)",
                "Blanket FLiBe cost (M$)",
                self.data.costs.c22127,
            )
            po.ocosts(
                self.outfile,
                "(c22128)",
                "Blanket lithium cost (M$)",
                self.data.costs.c22128,
            )

        po.ocosts(
            self.outfile,
            "(c2212)",
            "Blanket total cost (M$)",
            self.data.costs.c2212,
        )
        po.ocosts(
            self.outfile,
            "(c22131)",
            "Bulk shield cost (M$)",
            self.data.costs.c22131,
        )
        po.ocosts(
            self.outfile,
            "(c22132)",
            "Penetration shielding cost (M$)",
            self.data.costs.c22132,
        )
        po.ocosts(
            self.outfile,
            "(c2213)",
            "Total shield cost (M$)",
            self.data.costs.c2213,
        )
        po.ocosts(
            self.outfile,
            "(c2214)",
            "Total support structure cost (M$)",
            self.data.costs.c2214,
        )
        po.ocosts(self.outfile, "(c2215)", "Divertor cost (M$)", self.data.costs.c2215)
        # TODO: Convert fortran format to Python
        #     if (self.data.costs.ifueltyp == 1) :
        #         po.oblnkl(self.outfile)
        #         write(self.outfile,20)
        #     20     format(t2,             'First wall, total blanket and divertor direct costs',/,             t2,'are zero as they are assumed to be fuel costs.')
        #     elif  (self.data.costs.ifueltyp == 2) :
        #         po.oblnkl(self.outfile)
        #         write(self.outfile,31)
        # 21     format(t2,             'Initial First wall, total blanket and divertor direct costs',/,             t2,'are in capital and replacemnet are in cost of electricity')

        po.oblnkl(self.outfile)
        po.ocosts(
            self.outfile,
            "(c221)",
            "Total account 221 cost (M$)",
            self.data.costs.c221,
        )

        if ife_variables.ife != 1:
            po.oshead(self.outfile, "Magnets")

            if (
                tfcoil_variables.i_tf_sup != TFConductorModel.SUPERCONDUCTING
            ):  # Resistive TF coils
                if physics_variables.itart == 1:
                    po.ocosts(
                        self.outfile,
                        "(c22211)",
                        "Centrepost costs (M$)",
                        self.data.costs.c22211,
                    )
                else:
                    po.ocosts(
                        self.outfile,
                        "(c22211)",
                        "Inboard leg cost (M$)",
                        self.data.costs.c22211,
                    )

                po.ocosts(
                    self.outfile,
                    "(c22212)",
                    "Outboard leg cost (M$)",
                    self.data.costs.c22212,
                )
                po.ocosts(
                    self.outfile,
                    "(c2221)",
                    "TF magnet assemblies cost (M$)",
                    self.data.costs.c2221,
                )
            else:  # Superconducting TF coils
                po.ocosts(
                    self.outfile,
                    "(c22211)",
                    "TF coil conductor cost (M$)",
                    self.data.costs.c22211,
                )
                po.ocosts(
                    self.outfile,
                    "(c22212)",
                    "TF coil winding cost (M$)",
                    self.data.costs.c22212,
                )
                po.ocosts(
                    self.outfile,
                    "(c22213)",
                    "TF coil case cost (M$)",
                    self.data.costs.c22213,
                )
                po.ocosts(
                    self.outfile,
                    "(c22214)",
                    "TF intercoil structure cost (M$)",
                    self.data.costs.c22214,
                )
                po.ocosts(
                    self.outfile,
                    "(c22215)",
                    "TF coil gravity support structure (M$)",
                    self.data.costs.c22215,
                )
                po.ocosts(
                    self.outfile,
                    "(c2221)",
                    "TF magnet assemblies cost (M$)",
                    self.data.costs.c2221,
                )

            po.ocosts(
                self.outfile,
                "(c22221)",
                "PF coil conductor cost (M$)",
                self.data.costs.c22221,
            )
            po.ocosts(
                self.outfile,
                "(c22222)",
                "PF coil winding cost (M$)",
                self.data.costs.c22222,
            )
            po.ocosts(
                self.outfile,
                "(c22223)",
                "PF coil case cost (M$)",
                self.data.costs.c22223,
            )
            po.ocosts(
                self.outfile,
                "(c22224)",
                "PF coil support structure cost (M$)",
                self.data.costs.c22224,
            )
            po.ocosts(
                self.outfile,
                "(c2222)",
                "PF magnet assemblies cost (M$)",
                self.data.costs.c2222,
            )
            po.ocosts(
                self.outfile,
                "(c2223)",
                "Vacuum vessel assembly cost (M$)",
                self.data.costs.c2223,
            )
            # TODO: Convert fortran format to Python
            #     if ((physics_variables.itart == 1)and(self.data.costs.ifueltyp == 1)) :
            #         po.oblnkl(self.outfile)
            #         write(self.outfile,30)
            # 30        format(t2,                'Centrepost direct cost is zero, as it ',                'is assumed to be a fuel cost.')
            #     elif  ((physics_variables.itart == 1)and(self.data.costs.ifueltyp == 2)) :
            #         po.oblnkl(self.outfile)
            #         write(self.outfile,31)
            # 31        format(t2,                'Initial centrepost direct cost in included in capital ',                'cost and replacements are assumed to be a fuel cost.')

            po.oblnkl(self.outfile)
            po.ocosts(
                self.outfile,
                "(c222)",
                "Total account 222 cost (M$)",
                self.data.costs.c222,
            )

        po.oshead(self.outfile, "Power Injection")

        if ife_variables.ife == 1:
            po.ocosts(
                self.outfile,
                "(c2231)",
                "IFE driver system cost (M$)",
                self.data.costs.c2231,
            )
        else:
            po.ocosts(
                self.outfile,
                "(c2231)",
                "ECH system cost (M$)",
                self.data.costs.c2231,
            )
            po.ocosts(
                self.outfile,
                "(c2232)",
                "Lower hybrid system cost (M$)",
                self.data.costs.c2232,
            )
            po.ocosts(
                self.outfile,
                "(c2233)",
                "Neutral beam system cost (M$)",
                self.data.costs.c2233,
            )

        po.oblnkl(self.outfile)
        po.ocosts(
            self.outfile,
            "(c223)",
            "Total account 223 cost (M$)",
            self.data.costs.c223,
        )

        po.oshead(self.outfile, "Vacuum Systems")
        po.ocosts(
            self.outfile,
            "(c2241)",
            "High vacuum pumps cost (M$)",
            self.data.costs.c2241,
        )
        po.ocosts(
            self.outfile,
            "(c2242)",
            "Backing pumps cost (M$)",
            self.data.costs.c2242,
        )
        po.ocosts(
            self.outfile,
            "(c2243)",
            "Vacuum duct cost (M$)",
            self.data.costs.c2243,
        )
        po.ocosts(self.outfile, "(c2244)", "Valves cost (M$)", self.data.costs.c2244)
        po.ocosts(
            self.outfile,
            "(c2245)",
            "Duct shielding cost (M$)",
            self.data.costs.c2245,
        )
        po.ocosts(
            self.outfile,
            "(c2246)",
            "Instrumentation cost (M$)",
            self.data.costs.c2246,
        )
        po.oblnkl(self.outfile)
        po.ocosts(
            self.outfile,
            "(c224)",
            "Total account 224 cost (M$)",
            self.data.costs.c224,
        )

        if ife_variables.ife != 1:
            po.oshead(self.outfile, "Power Conditioning")
            po.ocosts(
                self.outfile,
                "(c22511)",
                "TF coil power supplies cost (M$)",
                self.data.costs.c22511,
            )
            po.ocosts(
                self.outfile,
                "(c22512)",
                "TF coil breakers cost (M$)",
                self.data.costs.c22512,
            )
            po.ocosts(
                self.outfile,
                "(c22513)",
                "TF coil dump resistors cost (M$)",
                self.data.costs.c22513,
            )
            po.ocosts(
                self.outfile,
                "(c22514)",
                "TF coil instrumentation and control (M$)",
                self.data.costs.c22514,
            )
            po.ocosts(
                self.outfile,
                "(c22515)",
                "TF coil bussing cost (M$)",
                self.data.costs.c22515,
            )
            po.ocosts(
                self.outfile,
                "(c2251)",
                "Total, TF coil power costs (M$)",
                self.data.costs.c2251,
            )
            po.ocosts(
                self.outfile,
                "(c22521)",
                "PF coil power supplies cost (M$)",
                self.data.costs.c22521,
            )
            po.ocosts(
                self.outfile,
                "(c22522)",
                "PF coil instrumentation and control (M$)",
                self.data.costs.c22522,
            )
            po.ocosts(
                self.outfile,
                "(c22523)",
                "PF coil bussing cost (M$)",
                self.data.costs.c22523,
            )
            po.ocosts(
                self.outfile,
                "(c22524)",
                "PF coil burn power supplies cost (M$)",
                self.data.costs.c22524,
            )
            po.ocosts(
                self.outfile,
                "(c22525)",
                "PF coil breakers cost (M$)",
                self.data.costs.c22525,
            )
            po.ocosts(
                self.outfile,
                "(c22526)",
                "PF coil dump resistors cost (M$)",
                self.data.costs.c22526,
            )
            po.ocosts(
                self.outfile,
                "(c22527)",
                "PF coil ac breakers cost (M$)",
                self.data.costs.c22527,
            )
            po.ocosts(
                self.outfile,
                "(c2252)",
                "Total, PF coil power costs (M$)",
                self.data.costs.c2252,
            )
            po.ocosts(
                self.outfile,
                "(c2253)",
                "Total, energy storage cost (M$)",
                self.data.costs.c2253,
            )
            po.oblnkl(self.outfile)
            po.ocosts(
                self.outfile,
                "(c225)",
                "Total account 225 cost (M$)",
                self.data.costs.c225,
            )

        po.oshead(self.outfile, "Heat Transport System")
        po.ocosts(
            self.outfile,
            "(cpp)",
            "Pumps and piping system cost (M$)",
            self.data.costs.cpp,
        )
        po.ocosts(
            self.outfile,
            "(chx)",
            "Primary heat exchanger cost (M$)",
            self.data.costs.chx,
        )
        po.ocosts(
            self.outfile,
            "(c2261)",
            "Total, reactor cooling system cost (M$)",
            self.data.costs.c2261,
        )
        po.ocosts(
            self.outfile,
            "(cppa)",
            "Pumps, piping cost (M$)",
            self.data.costs.cppa,
        )
        po.ocosts(
            self.outfile,
            "(c2262)",
            "Total, auxiliary cooling system cost (M$)",
            self.data.costs.c2262,
        )
        po.ocosts(
            self.outfile,
            "(c2263)",
            "Total, cryogenic system cost (M$)",
            self.data.costs.c2263,
        )
        po.oblnkl(self.outfile)
        po.ocosts(
            self.outfile,
            "(c226)",
            "Total account 226 cost (M$)",
            self.data.costs.c226,
        )

        po.oshead(self.outfile, "Fuel Handling System")
        po.ocosts(
            self.outfile,
            "(c2271)",
            "Fuelling system cost (M$)",
            self.data.costs.c2271,
        )
        po.ocosts(
            self.outfile,
            "(c2272)",
            "Fuel processing and purification cost (M$)",
            self.data.costs.c2272,
        )
        po.ocosts(
            self.outfile,
            "(c2273)",
            "Atmospheric recovery systems cost (M$)",
            self.data.costs.c2273,
        )
        po.ocosts(
            self.outfile,
            "(c2274)",
            "Nuclear building ventilation cost (M$)",
            self.data.costs.c2274,
        )
        po.oblnkl(self.outfile)
        po.ocosts(
            self.outfile,
            "(c227)",
            "Total account 227 cost (M$)",
            self.data.costs.c227,
        )

        po.oshead(self.outfile, "Instrumentation and Control")
        po.ocosts(
            self.outfile,
            "(c228)",
            "Instrumentation and control cost (M$)",
            self.data.costs.c228,
        )

        po.oshead(self.outfile, "Maintenance Equipment")
        po.ocosts(
            self.outfile,
            "(c229)",
            "Maintenance equipment cost (M$)",
            self.data.costs.c229,
        )

        po.oshead(self.outfile, "Total Account 22 Cost")
        po.ocosts(
            self.outfile,
            "(c22)",
            "Total account 22 cost (M$)",
            self.data.costs.c22,
        )

        po.oshead(self.outfile, "Turbine Plant Equipment")
        po.ocosts(
            self.outfile,
            "(c23)",
            "Turbine plant equipment cost (M$)",
            self.data.costs.c23,
        )

        po.oshead(self.outfile, "Electric Plant Equipment")
        po.ocosts(
            self.outfile,
            "(c241)",
            "Switchyard equipment cost (M$)",
            self.data.costs.c241,
        )
        po.ocosts(self.outfile, "(c242)", "Transformers cost (M$)", self.data.costs.c242)
        po.ocosts(
            self.outfile,
            "(c243)",
            "Low voltage equipment cost (M$)",
            self.data.costs.c243,
        )
        po.ocosts(
            self.outfile,
            "(c244)",
            "Diesel backup equipment cost (M$)",
            self.data.costs.c244,
        )
        po.ocosts(
            self.outfile,
            "(c245)",
            "Auxiliary facilities cost (M$)",
            self.data.costs.c245,
        )
        po.oblnkl(self.outfile)
        po.ocosts(
            self.outfile,
            "(c24)",
            "Total account 24 cost (M$)",
            self.data.costs.c24,
        )

        po.oshead(self.outfile, "Miscellaneous Plant Equipment")
        po.ocosts(
            self.outfile,
            "(c25)",
            "Miscellaneous plant equipment cost (M$)",
            self.data.costs.c25,
        )

        po.oshead(self.outfile, "Heat Rejection System")
        po.ocosts(
            self.outfile,
            "(c26)",
            "Heat rejection system cost (M$)",
            self.data.costs.c26,
        )

        po.oshead(self.outfile, "Plant Direct Cost")
        po.ocosts(
            self.outfile, "(cdirt)", "Plant direct cost (M$)", self.data.costs.cdirt
        )

        po.oshead(self.outfile, "Reactor Core Cost")
        po.ocosts(
            self.outfile,
            "(crctcore)",
            "Reactor core cost (M$)",
            self.data.costs.crctcore,
        )

        po.oshead(self.outfile, "Indirect Cost")
        po.ocosts(self.outfile, "(c9)", "Indirect cost (M$)", self.data.costs.cindrt)

        po.oshead(self.outfile, "Total Contingency")
        po.ocosts(
            self.outfile,
            "(ccont)",
            "Total contingency (M$)",
            self.data.costs.ccont,
        )

        po.oshead(self.outfile, "Constructed Cost")
        po.ocosts(
            self.outfile,
            "(concost)",
            "Constructed cost (M$)",
            self.data.costs.concost,
        )

        if self.data.costs.ireactor == 1:
            po.oshead(self.outfile, "Interest during Construction")
            po.ocosts(
                self.outfile,
                "(moneyint)",
                "Interest during construction (M$)",
                self.data.costs.moneyint,
            )

            po.oshead(self.outfile, "Total Capital Investment")
            po.ocosts(
                self.outfile,
                "(capcost)",
                "Total capital investment (M$)",
                self.data.costs.capcost,
            )

    def acc22(self):
        """Account 22 : Fusion power island

        This routine evaluates the Account 22 (fusion power island
        - the tokamak itself plus auxiliary power systems, etc.) costs.
        """
        self.acc221()

        #  Account 222 : Magnets
        self.acc222()

        #  Account 223 : Power injection
        self.acc223()

        #  Account 224 : Vacuum system
        self.acc224()

        #  Account 225 : Power conditioning
        self.acc225()

        #  Account 226 : Heat transport system
        self.acc2261()  # Account 2261 : Reactor cooling system
        self.acc2262()  # Account 2262 : Auxiliary component coolin
        self.acc2263()  # Account 2263 : Cryogenic system
        self.acc226()  # ccount 226  : Total

        #  Account 227 : Fuel handling
        self.acc2271()  # Account 2271 : Fuelling system
        self.acc2272()  # Account 2272 : Fuel processing and purification
        self.acc2273()  # Account 2273 : Atmospheric recovery systems
        self.acc2274()  # Account 2274 : Nuclear building ventilation
        self.acc227()  # Account 227  : Total

        #  Account 228 : Instrumentation and control
        self.acc228()

        #  Account 229 : Maintenance equipment
        self.acc229()

        #  Reactor core costs
        self.data.costs.crctcore = (
            self.data.costs.c221 + self.data.costs.c222 + self.data.costs.c223
        )

        #  Total account 22
        self.data.costs.c22 = (
            self.data.costs.c221
            + self.data.costs.c222
            + self.data.costs.c223
            + self.data.costs.c224
            + self.data.costs.c225
            + self.data.costs.c226
            + self.data.costs.c227
            + self.data.costs.c228
            + self.data.costs.c229
        )

    def acc221(self):
        """Account 221 : Reactor
        This routine evaluates the Account 221 (reactor) costs.
        These include the first wall, blanket, shield, support structure
        and divertor plates.
        <P>If ifueltyp = 1, the first wall, blanket and divertor costs are
        treated as fuel costs, rather than as capital costs.
        <P>If ifueltyp = 2, the initial first wall, blanket and divertor costs are
        treated as capital costs, and replacemnts are included as fuel costs.
        """
        self.acc2211()

        #  Account 221.2 : Blanket

        self.acc2212()

        #  Account 221.3 : Shield

        self.acc2213()

        #  Account 221.4 : Reactor structure

        self.acc2214()

        #  Account 221.5 : Divertor

        self.acc2215()

        #  Total account 221

        self.data.costs.c221 = (
            self.data.costs.c2211
            + self.data.costs.c2212
            + self.data.costs.c2213
            + self.data.costs.c2214
            + self.data.costs.c2215
        )

    def acc222(self):
        """Account 222 : Magnets, including cryostat
        This routine evaluates the Account 222 (magnet) costs,
        including the costs of associated cryostats.
        """
        if ife_variables.ife == 1:
            self.data.costs.c222 = 0.0e0
            return

        #  Account 222.1 : TF magnet assemblies

        self.acc2221()

        #  Account 222.2 : PF magnet assemblies
        self.acc2222()

        #  Account 222.3 : Cryostat

        self.acc2223()

        #  Total account 222

        self.data.costs.c222 = (
            self.data.costs.c2221 + self.data.costs.c2222 + self.data.costs.c2223
        )

    def acc225(self):
        """Account 225 : Power conditioning
        This routine evaluates the Account 225 (power conditioning) costs.
        """
        if ife_variables.ife == 1:
            self.data.costs.c225 = 0.0e0
        else:
            #  Account 225.1 : TF coil power conditioning

            self.acc2251()

            #  Account 225.2 : PF coil power conditioning

            self.acc2252()

            #  Account 225.3 : Energy storage

            self.acc2253()

            #  Total account 225

            self.data.costs.c225 = (
                self.data.costs.c2251 + self.data.costs.c2252 + self.data.costs.c2253
            )

    def acc21(self):
        """Account 21 : Structures and site facilities
        This routine evaluates the Account 21 (structures and site
        facilities) costs.
        Building costs are scaled with volume according to algorithms
        developed from TFCX, TFTR, and commercial power plant buildings.
        Costs include equipment, materials and installation labour, but
        no engineering or construction management.
        <P>The general form of the cost algorithm is cost=ucxx*volume**expxx.
        Allowances are used for site improvements and for miscellaneous
        buildings and land costs.
        """
        cmlsa = [0.6800e0, 0.8400e0, 0.9200e0, 1.0000e0]
        exprb = 1.0e0
        #  Account 211 : Site improvements, facilities and land
        #  N.B. Land unaffected by LSA

        self.data.costs.c211 = (
            self.data.costs.csi * cmlsa[self.data.costs.lsa - 1] + self.data.costs.cland
        )

        #  Account 212 : Reactor building

        self.data.costs.c212 = (
            1.0e-6
            * self.data.costs.ucrb
            * buildings_variables.rbvol**exprb
            * cmlsa[self.data.costs.lsa - 1]
        )

        #  Account 213 : Turbine building

        if self.data.costs.ireactor == 1:
            self.data.costs.c213 = (
                self.data.costs.cturbb * cmlsa[self.data.costs.lsa - 1]
            )
        else:
            self.data.costs.c213 = 0.0e0

        #  Account 214 : Reactor maintenance and warm shops buildings

        self.data.costs.c2141 = (
            1.0e-6
            * self.data.costs.UCMB
            * buildings_variables.rmbvol**exprb
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c2142 = (
            1.0e-6
            * self.data.costs.UCWS
            * buildings_variables.wsvol**exprb
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c214 = self.data.costs.c2141 + self.data.costs.c2142

        #  Account 215 : Tritium building

        self.data.costs.c215 = (
            1.0e-6
            * self.data.costs.UCTR
            * buildings_variables.triv**exprb
            * cmlsa[self.data.costs.lsa - 1]
        )

        #  Account 216 : Electrical equipment building

        self.data.costs.c216 = (
            1.0e-6
            * self.data.costs.UCEL
            * buildings_variables.elevol**exprb
            * cmlsa[self.data.costs.lsa - 1]
        )

        #  Account 217 : Other buildings
        #  Includes administration, control, shops, cryogenic
        #  plant and an allowance for miscellaneous structures

        self.data.costs.c2171 = (
            1.0e-6
            * self.data.costs.UCAD
            * buildings_variables.admvol**exprb
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c2172 = (
            1.0e-6
            * self.data.costs.UCCO
            * buildings_variables.convol**exprb
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c2173 = (
            1.0e-6
            * self.data.costs.UCSH
            * buildings_variables.shovol**exprb
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c2174 = (
            1.0e-6
            * self.data.costs.UCCR
            * buildings_variables.cryvol**exprb
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c217 = (
            self.data.costs.c2171
            + self.data.costs.c2172
            + self.data.costs.c2173
            + self.data.costs.c2174
        )

        #  Total for Account 21

        self.data.costs.c21 = (
            self.data.costs.c211
            + self.data.costs.c212
            + self.data.costs.c213
            + self.data.costs.c214
            + self.data.costs.c215
            + self.data.costs.c216
            + self.data.costs.c217
        )

    def acc2211(self):
        """Account 221.1 : First wall
        This routine evaluates the Account 221.1 (first wall) costs.
        The first wall cost is scaled linearly with surface area from TFCX.
        If ifueltyp = 1, the first wall cost is treated as a fuel cost,
        rather than as a capital cost.
        If ifueltyp = 2, inital first wall is included as a capital cost,
        and the replacement first wall cost is treated as a fuel costs.
        """
        cmlsa = [0.5000e0, 0.7500e0, 0.8750e0, 1.0000e0]

        if ife_variables.ife != 1:
            self.data.costs.c2211 = (
                1.0e-6
                * cmlsa[self.data.costs.lsa - 1]
                * (
                    (self.data.costs.UCFWA + self.data.costs.UCFWS)
                    * self.data.first_wall.a_fw_total
                    + self.data.costs.UCFWPS
                )
            )
        else:
            self.data.costs.c2211 = (
                1.0e-6
                * cmlsa[self.data.costs.lsa - 1]
                * (
                    self.data.costs.ucblss
                    * (
                        ife_variables.fwmatm[0, 0]
                        + ife_variables.fwmatm[1, 0]
                        + ife_variables.fwmatm[2, 0]
                    )
                    + ife_variables.uccarb
                    * (
                        ife_variables.fwmatm[0, 1]
                        + ife_variables.fwmatm[1, 1]
                        + ife_variables.fwmatm[2, 1]
                    )
                    + self.data.costs.ucblli2o
                    * (
                        ife_variables.fwmatm[0, 3]
                        + ife_variables.fwmatm[1, 3]
                        + ife_variables.fwmatm[2, 3]
                    )
                    + ife_variables.ucconc
                    * (
                        ife_variables.fwmatm[0, 4]
                        + ife_variables.fwmatm[1, 4]
                        + ife_variables.fwmatm[2, 4]
                    )
                )
            )

        self.data.costs.c2211 = self.data.costs.fkind * self.data.costs.c2211

        if self.data.costs.ifueltyp == 1:
            self.data.costs.fwallcst = self.data.costs.c2211
            self.data.costs.c2211 = 0.0e0
        elif self.data.costs.ifueltyp == 2:
            self.data.costs.fwallcst = self.data.costs.c2211
        else:
            self.data.costs.fwallcst = 0.0e0

    def acc2212(self):
        """Account 221.2 : Blanket
        This routine evaluates the Account 221.2 (blanket) costs.
        If ifueltyp = 1, the blanket cost is treated as a fuel cost,
        rather than as a capital cost.
        If ifueltyp = 2, the initial blanket is included as a capital cost
        and the replacement blanket costs are treated as a fuel cost.
        """
        cmlsa = [0.5000e0, 0.7500e0, 0.8750e0, 1.0000e0]

        if ife_variables.ife != 1:
            #  Solid blanket (Li2O + Be)
            self.data.costs.c22121 = (
                1.0e-6 * fwbs_variables.m_blkt_beryllium * self.data.costs.ucblbe
            )

            # CCFE model
            self.data.costs.c22122 = (
                1.0e-6 * fwbs_variables.m_blkt_li2o * self.data.costs.ucblli2o
            )

            self.data.costs.c22123 = (
                1.0e-6 * fwbs_variables.m_blkt_steel_total * self.data.costs.ucblss
            )
            self.data.costs.c22124 = (
                1.0e-6 * fwbs_variables.m_blkt_vanadium * self.data.costs.ucblvd
            )
            self.data.costs.c22125 = 0.0e0
            self.data.costs.c22126 = 0.0e0
            self.data.costs.c22127 = 0.0e0

        else:
            #  IFE blanket; materials present are Li2O, steel, carbon, concrete,
            #  FLiBe and lithium

            self.data.costs.c22121 = 0.0e0
            self.data.costs.c22122 = (
                1.0e-6 * fwbs_variables.m_blkt_li2o * self.data.costs.ucblli2o
            )
            self.data.costs.c22123 = (
                1.0e-6 * fwbs_variables.m_blkt_steel_total * self.data.costs.ucblss
            )
            self.data.costs.c22124 = 0.0e0
            self.data.costs.c22125 = (
                1.0e-6
                * ife_variables.uccarb
                * (
                    ife_variables.blmatm[0, 1]
                    + ife_variables.blmatm[1, 1]
                    + ife_variables.blmatm[2, 1]
                )
            )
            self.data.costs.c22126 = (
                1.0e-6
                * ife_variables.ucconc
                * (
                    ife_variables.blmatm[0, 4]
                    + ife_variables.blmatm[1, 4]
                    + ife_variables.blmatm[2, 4]
                )
            )
            self.data.costs.c22127 = 1.0e-6 * ife_variables.ucflib * ife_variables.mflibe
            self.data.costs.c22128 = (
                1.0e-6 * self.data.costs.ucblli * fwbs_variables.m_blkt_lithium
            )

        self.data.costs.c22121 = (
            self.data.costs.fkind
            * self.data.costs.c22121
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c22122 = (
            self.data.costs.fkind
            * self.data.costs.c22122
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c22123 = (
            self.data.costs.fkind
            * self.data.costs.c22123
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c22124 = (
            self.data.costs.fkind
            * self.data.costs.c22124
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c22125 = (
            self.data.costs.fkind
            * self.data.costs.c22125
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c22126 = (
            self.data.costs.fkind
            * self.data.costs.c22126
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c22127 = (
            self.data.costs.fkind
            * self.data.costs.c22127
            * cmlsa[self.data.costs.lsa - 1]
        )

        self.data.costs.c2212 = (
            self.data.costs.c22121
            + self.data.costs.c22122
            + self.data.costs.c22123
            + self.data.costs.c22124
            + self.data.costs.c22125
            + self.data.costs.c22126
            + self.data.costs.c22127
        )

        if self.data.costs.ifueltyp == 1:
            self.data.costs.blkcst = self.data.costs.c2212
            self.data.costs.c2212 = 0.0e0
        elif self.data.costs.ifueltyp == 2:
            self.data.costs.blkcst = self.data.costs.c2212
        else:
            self.data.costs.blkcst = 0.0e0

    def acc2213(self):
        """Account 221.3 : Shield
        This routine evaluates the Account 221.3 (shield) costs.
        """
        cmlsa = [0.5000e0, 0.7500e0, 0.8750e0, 1.0000e0]

        if ife_variables.ife != 1:
            self.data.costs.c22131 = (
                1.0e-6
                * fwbs_variables.whtshld
                * self.data.costs.ucshld
                * cmlsa[self.data.costs.lsa - 1]
            )
        else:
            self.data.costs.c22131 = (
                1.0e-6
                * cmlsa[self.data.costs.lsa - 1]
                * (
                    self.data.costs.ucshld
                    * (
                        ife_variables.shmatm[0, 0]
                        + ife_variables.shmatm[1, 0]
                        + ife_variables.shmatm[2, 0]
                    )
                    + ife_variables.uccarb
                    * (
                        ife_variables.shmatm[0, 1]
                        + ife_variables.shmatm[1, 1]
                        + ife_variables.shmatm[2, 1]
                    )
                    + self.data.costs.ucblli2o
                    * (
                        ife_variables.shmatm[0, 3]
                        + ife_variables.shmatm[1, 3]
                        + ife_variables.shmatm[2, 3]
                    )
                    + ife_variables.ucconc
                    * (
                        ife_variables.shmatm[0, 4]
                        + ife_variables.shmatm[1, 4]
                        + ife_variables.shmatm[2, 4]
                    )
                )
            )

        self.data.costs.c22131 = self.data.costs.fkind * self.data.costs.c22131

        #  Penetration shield assumed to be typical steel plate
        if ife_variables.ife != 1:
            self.data.costs.c22132 = (
                1.0e-6
                * fwbs_variables.wpenshld
                * self.data.costs.ucpens
                * cmlsa[self.data.costs.lsa - 1]
            )
        else:
            self.data.costs.c22132 = 0.0e0

        self.data.costs.c22132 = self.data.costs.fkind * self.data.costs.c22132

        self.data.costs.c2213 = self.data.costs.c22131 + self.data.costs.c22132

    def acc2214(self):
        """Account 221.4 : Reactor structure
        This routine evaluates the Account 221.4 (reactor structure) costs.
        The structural items are costed as standard steel elements.
        """
        cmlsa = [0.6700e0, 0.8350e0, 0.9175e0, 1.0000e0]

        self.data.costs.c2214 = (
            1.0e-6
            * structure_variables.gsmass
            * self.data.costs.UCGSS
            * cmlsa[self.data.costs.lsa - 1]
        )
        self.data.costs.c2214 = self.data.costs.fkind * self.data.costs.c2214

    def acc2215(self):
        """Account 221.5 : Divertor
        This routine evaluates the Account 221.5 (divertor) costs.
        The cost of the divertor blade is scaled linearly with
        surface area from TFCX. The graphite armour is assumed to
        be brazed to water-cooled machined copper substrate.
        Tenth-of-a-kind engineering and installation is assumed.
        <P>If ifueltyp = 1, the divertor cost is treated as a fuel cost,
        rather than as a capital cost.
        <P>If ifueltyp = 2, the initial divertor is included as a capital cost
        and the replacement divertor costs ae treated as a fuel cost,
        """
        if ife_variables.ife != 1:
            self.data.costs.c2215 = (
                1.0e-6 * divertor_variables.a_div_surface_total * self.data.costs.ucdiv
            )
            self.data.costs.c2215 = self.data.costs.fkind * self.data.costs.c2215

            if self.data.costs.ifueltyp == 1:
                self.data.costs.divcst = self.data.costs.c2215
                self.data.costs.c2215 = 0.0e0
            elif self.data.costs.ifueltyp == 2:
                self.data.costs.divcst = self.data.costs.c2215
            else:
                self.data.costs.divcst = 0.0e0

        else:
            self.data.costs.c2215 = 0.0e0
            self.data.costs.divcst = 0.0e0

    def acc2221(self):
        """Account 222.1 : TF magnet assemblies
        This routine evaluates the Account 222.1 (TF magnet) costs.
        Copper magnets are costed from the TFCX data base ($/kg).
        Superconductor magnets are costed using a new method devised
        by R. Hancox under contract to Culham Laboratory, Jan/Feb 1994.
        If ifueltyp = 1, the TART centrepost cost is treated as a fuel
        cost, rather than as a capital cost.
        If ifueltyp = 2, the  initial centrepost is included as a capital cost
        and the replacement TART centrepost costs are treated as a fuel
        """
        cmlsa = [0.6900e0, 0.8450e0, 0.9225e0, 1.0000e0]

        if (
            tfcoil_variables.i_tf_sup != TFConductorModel.SUPERCONDUCTING
        ):  # Resistive TF coils
            #  Account 222.1.1 : Inboard TF coil legs

            self.data.costs.c22211 = (
                1.0e-6
                * tfcoil_variables.whtcp
                * self.data.costs.uccpcl1
                * cmlsa[self.data.costs.lsa - 1]
            )
            self.data.costs.c22211 = self.data.costs.fkind * self.data.costs.c22211

            self.data.costs.cpstcst = 0.0e0  # TART centrepost
            if (physics_variables.itart == 1) and (self.data.costs.ifueltyp == 1):
                self.data.costs.cpstcst = self.data.costs.c22211
                self.data.costs.c22211 = 0.0e0
            elif (physics_variables.itart == 1) and (self.data.costs.ifueltyp == 2):
                self.data.costs.cpstcst = self.data.costs.c22211

            #  Account 222.1.2 : Outboard TF coil legs

            self.data.costs.c22212 = (
                1.0e-6
                * tfcoil_variables.whttflgs
                * self.data.costs.uccpclb
                * cmlsa[self.data.costs.lsa - 1]
            )
            self.data.costs.c22212 = self.data.costs.fkind * self.data.costs.c22212

            #  Total (copper) TF coil costs

            self.data.costs.c2221 = self.data.costs.c22211 + self.data.costs.c22212

        else:  # Superconducting TF coils
            #  Account 222.1.1 : Conductor

            #  Superconductor ($/m)

            if self.data.costs.supercond_cost_model == 0:
                costtfsc = (
                    self.data.costs.ucsc[tfcoil_variables.i_tf_sc_mat - 1]
                    * tfcoil_variables.m_tf_coil_superconductor
                    / (tfcoil_variables.len_tf_coil * tfcoil_variables.n_tf_coil_turns)
                )
            else:
                costtfsc = (
                    self.data.costs.sc_mat_cost_0[tfcoil_variables.i_tf_sc_mat - 1]
                    * tfcoil_variables.j_crit_str_0[tfcoil_variables.i_tf_sc_mat - 1]
                    / tfcoil_variables.j_crit_str_tf
                )

            #  Copper ($/m)

            costtfcu = (
                self.data.costs.uccu
                * tfcoil_variables.m_tf_coil_copper
                / (tfcoil_variables.len_tf_coil * tfcoil_variables.n_tf_coil_turns)
            )

            #  Total cost/metre of superconductor and copper wire

            costwire = costtfsc + costtfcu

            #  Total cost/metre of conductor (including sheath and fixed costs)

            ctfconpm = costwire + self.data.costs.cconshtf + self.data.costs.cconfix

            #  Total conductor costs

            self.data.costs.c22211 = (
                1.0e-6
                * ctfconpm
                * tfcoil_variables.n_tf_coils
                * tfcoil_variables.len_tf_coil
                * tfcoil_variables.n_tf_coil_turns
            )
            self.data.costs.c22211 = (
                self.data.costs.fkind
                * self.data.costs.c22211
                * cmlsa[self.data.costs.lsa - 1]
            )

            #  Account 222.1.2 : Winding

            self.data.costs.c22212 = (
                1.0e-6
                * self.data.costs.ucwindtf
                * tfcoil_variables.n_tf_coils
                * tfcoil_variables.len_tf_coil
                * tfcoil_variables.n_tf_coil_turns
            )
            self.data.costs.c22212 = (
                self.data.costs.fkind
                * self.data.costs.c22212
                * cmlsa[self.data.costs.lsa - 1]
            )

            #  Account 222.1.3 : Case

            self.data.costs.c22213 = (
                1.0e-6
                * (tfcoil_variables.m_tf_coil_case * self.data.costs.uccase)
                * tfcoil_variables.n_tf_coils
            )
            self.data.costs.c22213 = (
                self.data.costs.fkind
                * self.data.costs.c22213
                * cmlsa[self.data.costs.lsa - 1]
            )

            #  Account 222.1.4 : Intercoil structure

            self.data.costs.c22214 = (
                1.0e-6 * structure_variables.aintmass * self.data.costs.UCINT
            )
            self.data.costs.c22214 = (
                self.data.costs.fkind
                * self.data.costs.c22214
                * cmlsa[self.data.costs.lsa - 1]
            )

            #  Account 222.1.5 : Gravity support structure

            self.data.costs.c22215 = (
                1.0e-6 * structure_variables.clgsmass * self.data.costs.UCGSS
            )
            self.data.costs.c22215 = (
                self.data.costs.fkind
                * self.data.costs.c22215
                * cmlsa[self.data.costs.lsa - 1]
            )

            #  Total (superconducting) TF coil costs

            self.data.costs.c2221 = (
                self.data.costs.c22211
                + self.data.costs.c22212
                + self.data.costs.c22213
                + self.data.costs.c22214
                + self.data.costs.c22215
            )

    def acc2222(self):
        """Account 222.2 : PF magnet assemblies
        This routine evaluates the Account 222.2 (PF magnet) costs.
        Conductor costs previously used an algorithm devised by R. Hancox,
        January 1994, under contract to Culham, which took into
        account the fact that the superconductor/copper ratio in
        the conductor is proportional to the maximum field that
        each coil will experience. Now, the input copper fractions
        are used instead.
        Maximum values for current, current density and field
        are used.
        """
        cmlsa = [0.6900e0, 0.8450e0, 0.9225e0, 1.0000e0]

        #  Total length of PF coil windings (m)

        pfwndl = 0.0e0
        for i in range(pfcoil_variables.n_cs_pf_coils):
            pfwndl += (
                2.0
                * np.pi
                * pfcoil_variables.r_pf_coil_middle[i]
                * pfcoil_variables.n_pf_coil_turns[i]
            )

        #  Account 222.2.1 : Conductor

        #  The following lines take care of resistive coils.
        #  costpfsh is the cost per metre of the steel conduit/sheath around
        #  each superconducting cable (so is zero for resistive coils)

        costpfsh = (
            0.0 if pfcoil_variables.i_pf_conductor == 1 else self.data.costs.cconshpf
        )

        #  Non-Central Solenoid coils

        if build_variables.iohcl == 1:
            npf = pfcoil_variables.n_cs_pf_coils - 1
        else:
            npf = pfcoil_variables.n_cs_pf_coils

        self.data.costs.c22221 = 0.0e0

        for i in range(npf):
            #  Superconductor ($/m)
            if self.data.costs.supercond_cost_model == 0:
                if pfcoil_variables.i_pf_conductor == 0:
                    costpfsc = (
                        self.data.costs.ucsc[pfcoil_variables.i_pf_superconductor - 1]
                        * (1.0e0 - pfcoil_variables.fcupfsu)
                        * (1.0e0 - pfcoil_variables.f_a_pf_coil_void[i])
                        * abs(
                            pfcoil_variables.c_pf_cs_coils_peak_ma[i]
                            / pfcoil_variables.n_pf_coil_turns[i]
                        )
                        * 1.0e6
                        / pfcoil_variables.j_pf_coil_wp_peak[i]
                        * tfcoil_variables.dcond[
                            pfcoil_variables.i_pf_superconductor - 1
                        ]
                    )
                else:
                    costpfsc = 0.0e0
            elif pfcoil_variables.i_pf_conductor == 0:
                costpfsc = (
                    self.data.costs.sc_mat_cost_0[
                        pfcoil_variables.i_pf_superconductor - 1
                    ]
                    * tfcoil_variables.j_crit_str_0[
                        pfcoil_variables.i_pf_superconductor - 1
                    ]
                    / pfcoil_variables.j_crit_str_pf
                )
            else:
                costpfsc = 0.0

            #  Copper ($/m)
            if pfcoil_variables.i_pf_conductor == 0:
                costpfcu = (
                    self.data.costs.uccu
                    * pfcoil_variables.fcupfsu
                    * (1.0e0 - pfcoil_variables.f_a_pf_coil_void[i])
                    * abs(
                        pfcoil_variables.c_pf_cs_coils_peak_ma[i]
                        / pfcoil_variables.n_pf_coil_turns[i]
                    )
                    * 1.0e6
                    / pfcoil_variables.j_pf_coil_wp_peak[i]
                    * constants.den_copper
                )
            else:
                costpfcu = (
                    self.data.costs.uccu
                    * (1.0e0 - pfcoil_variables.f_a_pf_coil_void[i])
                    * abs(
                        pfcoil_variables.c_pf_cs_coils_peak_ma[i]
                        / pfcoil_variables.n_pf_coil_turns[i]
                    )
                    * 1.0e6
                    / pfcoil_variables.j_pf_coil_wp_peak[i]
                    * constants.den_copper
                )

            #  Total cost/metre of superconductor and copper wire

            costwire = costpfsc + costpfcu

            #  Total cost/metre of conductor (including sheath and fixed costs)

            cpfconpm = costwire + costpfsh + self.data.costs.cconfix

            #  Total account 222.2.1 (PF coils excluding Central Solenoid)

            self.data.costs.c22221 += (
                1.0e-6
                * 2.0
                * np.pi
                * pfcoil_variables.r_pf_coil_middle[i]
                * pfcoil_variables.n_pf_coil_turns[i]
                * cpfconpm
            )

        #  Central Solenoid

        if build_variables.iohcl == 1:
            #  Superconductor ($/m)
            if self.data.costs.supercond_cost_model == 0:
                #  Issue #328  Use CS conductor cross-sectional area (m2)
                if pfcoil_variables.i_pf_conductor == 0:
                    costpfsc = (
                        self.data.costs.ucsc[pfcoil_variables.i_cs_superconductor - 1]
                        * pfcoil_variables.awpoh
                        * (1 - pfcoil_variables.f_a_cs_void)
                        * (1 - pfcoil_variables.fcuohsu)
                        / pfcoil_variables.n_pf_coil_turns[
                            pfcoil_variables.n_cs_pf_coils - 1
                        ]
                        * tfcoil_variables.dcond[
                            pfcoil_variables.i_cs_superconductor - 1
                        ]
                    )
                else:
                    costpfsc = 0.0e0
            elif pfcoil_variables.i_pf_conductor == 0:
                costpfsc = (
                    self.data.costs.sc_mat_cost_0[
                        pfcoil_variables.i_cs_superconductor - 1
                    ]
                    * tfcoil_variables.j_crit_str_0[
                        pfcoil_variables.i_cs_superconductor - 1
                    ]
                    / pfcoil_variables.j_crit_str_cs
                )
            else:
                costpfsc = 0.0e0

            #  Copper ($/m)

            if pfcoil_variables.i_pf_conductor == 0:
                costpfcu = (
                    self.data.costs.uccu
                    * pfcoil_variables.awpoh
                    * (1 - pfcoil_variables.f_a_cs_void)
                    * pfcoil_variables.fcuohsu
                    / pfcoil_variables.n_pf_coil_turns[
                        pfcoil_variables.n_cs_pf_coils - 1
                    ]
                    * constants.den_copper
                )
            else:
                # MDK I don't know if this is ccorrect as we never use the resistive model
                costpfcu = (
                    self.data.costs.uccu
                    * pfcoil_variables.awpoh
                    * (1 - pfcoil_variables.f_a_cs_void)
                    / pfcoil_variables.n_pf_coil_turns[
                        pfcoil_variables.n_cs_pf_coils - 1
                    ]
                    * constants.den_copper
                )

            #  Total cost/metre of superconductor and copper wire (Central Solenoid)

            costwire = costpfsc + costpfcu

            #  Total cost/metre of conductor (including sheath and fixed costs)

            cpfconpm = costwire + costpfsh + self.data.costs.cconfix

            #  Total account 222.2.1 (PF+Central Solenoid coils)

            self.data.costs.c22221 += (
                1.0e-6
                * 2.0
                * np.pi
                * pfcoil_variables.r_pf_coil_middle[pfcoil_variables.n_cs_pf_coils - 1]
                * pfcoil_variables.n_pf_coil_turns[pfcoil_variables.n_cs_pf_coils - 1]
                * cpfconpm
            )

        self.data.costs.c22221 = (
            self.data.costs.fkind
            * self.data.costs.c22221
            * cmlsa[self.data.costs.lsa - 1]
        )

        #  Account 222.2.2 : Winding

        self.data.costs.c22222 = 1.0e-6 * self.data.costs.ucwindpf * pfwndl
        self.data.costs.c22222 = (
            self.data.costs.fkind
            * self.data.costs.c22222
            * cmlsa[self.data.costs.lsa - 1]
        )

        #  Account 222.2.3 : Steel case - will be zero for resistive coils

        self.data.costs.c22223 = (
            1.0e-6 * self.data.costs.uccase * pfcoil_variables.m_pf_coil_structure_total
        )
        self.data.costs.c22223 = (
            self.data.costs.fkind
            * self.data.costs.c22223
            * cmlsa[self.data.costs.lsa - 1]
        )

        #  Account 222.2.4 : Support structure

        self.data.costs.c22224 = (
            1.0e-6 * self.data.costs.ucfnc * structure_variables.fncmass
        )
        self.data.costs.c22224 = (
            self.data.costs.fkind
            * self.data.costs.c22224
            * cmlsa[self.data.costs.lsa - 1]
        )

        #  Total account 222.2

        self.data.costs.c2222 = (
            self.data.costs.c22221
            + self.data.costs.c22222
            + self.data.costs.c22223
            + self.data.costs.c22224
        )

    def acc2223(self):
        """Account 222.3 : Vacuum vessel
        This routine evaluates the Account 222.3 (vacuum vessel) costs.
        """
        cmlsa = [0.6900e0, 0.8450e0, 0.9225e0, 1.0000e0]

        self.data.costs.c2223 = 1.0e-6 * fwbs_variables.m_vv * self.data.costs.uccryo
        self.data.costs.c2223 = (
            self.data.costs.fkind
            * self.data.costs.c2223
            * cmlsa[self.data.costs.lsa - 1]
        )

    def acc223(self):
        """Account 223 : Power injection
        This routine evaluates the Account 223 (power injection) costs.
        The costs are from TETRA, updated to 1990$.
        Nominal TIBER values are used pending system designs. Costs are
        scaled linearly with power injected into the plasma and include
        the power supplies.
        <P>If ifueltyp=1, the fraction (1-fcdfuel) of the cost of the
        current drive system is considered as capital cost, and the
        fraction (fcdfuel) is considered a recurring fuel cost due
        to the system's short life.
        """
        exprf = 1.0e0
        if ife_variables.ife != 1:
            #  Account 223.1 : ECH

            self.data.costs.c2231 = (
                1.0e-6
                * self.data.costs.ucech
                * (1.0e6 * current_drive_variables.p_hcd_ecrh_injected_total_mw) ** exprf
            )

            if self.data.costs.ifueltyp == 1:
                self.data.costs.c2231 = (
                    1.0e0 - self.data.costs.fcdfuel
                ) * self.data.costs.c2231
                self.data.costs.c2231 = self.data.costs.fkind * self.data.costs.c2231

            #  Account 223.2 : Lower Hybrid or ICH

            if current_drive_variables.i_hcd_primary != 2:
                self.data.costs.c2232 = (
                    1.0e-6
                    * self.data.costs.uclh
                    * (1.0e6 * current_drive_variables.p_hcd_lowhyb_injected_total_mw)
                    ** exprf
                )
            else:
                self.data.costs.c2232 = (
                    1.0e-6
                    * self.data.costs.ucich
                    * (1.0e6 * current_drive_variables.p_hcd_lowhyb_injected_total_mw)
                    ** exprf
                )

            if self.data.costs.ifueltyp == 1:
                self.data.costs.c2232 = (
                    1.0e0 - self.data.costs.fcdfuel
                ) * self.data.costs.c2232
                self.data.costs.c2232 = self.data.costs.fkind * self.data.costs.c2232

                #  Account 223.3 : Neutral Beam

                # self.data.costs.c2233 = 1.0e-6 * self.data.costs.ucnbi * (1.0e6*p_hcd_beam_injected_total_mw)**exprf
                # #327

                self.data.costs.c2233 = (
                    1.0e-6
                    * self.data.costs.ucnbi
                    * (1.0e6 * current_drive_variables.p_beam_injected_mw) ** exprf
                )

            if self.data.costs.ifueltyp == 1:
                self.data.costs.c2233 = (
                    1.0e0 - self.data.costs.fcdfuel
                ) * self.data.costs.c2233
                self.data.costs.c2233 = self.data.costs.fkind * self.data.costs.c2233

        else:
            #  IFE driver costs (depends on driver type)
            #  Assume offset linear form for generic and SOMBRERO types,
            #  or one of two offset linear forms for OSIRIS type

            if ife_variables.ifedrv == 2:
                if ife_variables.dcdrv1 <= ife_variables.dcdrv2:
                    switch = 0.0e0
                else:
                    switch = (ife_variables.cdriv2 - ife_variables.cdriv1) / (
                        ife_variables.dcdrv1 - ife_variables.dcdrv2
                    )

                if ife_variables.edrive <= switch:
                    self.data.costs.c2231 = ife_variables.mcdriv * (
                        ife_variables.cdriv1
                        + ife_variables.dcdrv1 * 1.0e-6 * ife_variables.edrive
                    )
                else:
                    self.data.costs.c2231 = ife_variables.mcdriv * (
                        ife_variables.cdriv2
                        + ife_variables.dcdrv2 * 1.0e-6 * ife_variables.edrive
                    )

            elif ife_variables.ifedrv == 3:
                self.data.costs.c2231 = (
                    ife_variables.mcdriv
                    * 1.0e-6
                    * ife_variables.cdriv3
                    * (ife_variables.edrive / ife_variables.etadrv)
                )
            else:
                self.data.costs.c2231 = ife_variables.mcdriv * (
                    ife_variables.cdriv0
                    + ife_variables.dcdrv0 * 1.0e-6 * ife_variables.edrive
                )

            if self.data.costs.ifueltyp == 1:
                self.data.costs.c2231 = (
                    1.0e0 - self.data.costs.fcdfuel
                ) * self.data.costs.c2231
                self.data.costs.c2231 = self.data.costs.fkind * self.data.costs.c2231
                self.data.costs.c2232 = 0.0e0
                self.data.costs.c2233 = 0.0e0
                self.data.costs.c2234 = 0.0e0

        #  Total account 223

        self.data.costs.c223 = (
            self.data.costs.c2231
            + self.data.costs.c2232
            + self.data.costs.c2233
            + self.data.costs.c2234
        )
        self.data.costs.cdcost = self.data.costs.c223

    def acc224(self):
        """Account 224 : Vacuum system
        This routine evaluates the Account 224 (vacuum system) costs.
        The costs are scaled from TETRA reactor code runs.
        """
        if self.data.vacuum.i_vacuum_pump_type == 1:
            self.data.costs.c2241 = (
                1.0e-6 * self.data.vacuum.n_vac_pumps_high * self.data.costs.UCCPMP
            )
        else:
            self.data.costs.c2241 = (
                1.0e-6 * self.data.vacuum.n_vac_pumps_high * self.data.costs.UCTPMP
            )

        self.data.costs.c2241 = self.data.costs.fkind * self.data.costs.c2241

        #  Account 224.2 : Backing pumps

        self.data.costs.c2242 = (
            1.0e-6 * self.data.vacuum.n_vv_vacuum_ducts * self.data.costs.UCBPMP
        )
        self.data.costs.c2242 = self.data.costs.fkind * self.data.costs.c2242

        #  Account 224.3 : Vacuum duct

        self.data.costs.c2243 = (
            1.0e-6
            * self.data.vacuum.n_vv_vacuum_ducts
            * self.data.vacuum.dlscal
            * self.data.costs.UCDUCT
        )
        self.data.costs.c2243 = self.data.costs.fkind * self.data.costs.c2243

        #  Account 224.4 : Valves

        self.data.costs.c2244 = (
            1.0e-6
            * 2.0e0
            * self.data.vacuum.n_vv_vacuum_ducts
            * (self.data.vacuum.dia_vv_vacuum_ducts * 1.2e0) ** 1.4e0
            * self.data.costs.UCVALV
        )
        self.data.costs.c2244 = self.data.costs.fkind * self.data.costs.c2244

        #  Account 224.5 : Duct shielding

        self.data.costs.c2245 = (
            1.0e-6
            * self.data.vacuum.n_vv_vacuum_ducts
            * self.data.vacuum.m_vv_vacuum_duct_shield
            * self.data.costs.UCVDSH
        )
        self.data.costs.c2245 = self.data.costs.fkind * self.data.costs.c2245

        #  Account 224.6 : Instrumentation

        self.data.costs.c2246 = 1.0e-6 * self.data.costs.UCVIAC
        self.data.costs.c2246 = self.data.costs.fkind * self.data.costs.c2246

        #  Total account 224

        self.data.costs.c224 = (
            self.data.costs.c2241
            + self.data.costs.c2242
            + self.data.costs.c2243
            + self.data.costs.c2244
            + self.data.costs.c2245
            + self.data.costs.c2246
        )

    def acc2251(self):
        """Account 225.1 : TF coil power conditioning
        This routine evaluates the Account 225.1 (TF coil power
        conditioning) costs.
        Costs are developed based on the major equipment specification
        of the tfcpwr module.  A multiplier is used to account for bulk
        materials and installation.
        """
        expel = 0.7e0
        self.data.costs.c22511 = (
            1.0e-6
            * self.data.costs.uctfps
            * (tfcoil_variables.tfckw * 1.0e3 + tfcoil_variables.tfcmw * 1.0e6) ** expel
        )
        self.data.costs.c22511 = self.data.costs.fkind * self.data.costs.c22511

        #  Account 225.1.2 : TF coil breakers (zero cost for copper coils)

        if tfcoil_variables.i_tf_sup == TFConductorModel.SUPERCONDUCTING:
            self.data.costs.c22512 = 1.0e-6 * (
                self.data.costs.uctfbr
                * tfcoil_variables.n_tf_coils
                * (
                    tfcoil_variables.c_tf_turn
                    * tfcoil_variables.v_tf_coil_dump_quench_kv
                    * 1.0e3
                )
                ** expel
                + self.data.costs.uctfsw * tfcoil_variables.c_tf_turn
            )
        else:
            self.data.costs.c22512 = 0.0e0

        self.data.costs.c22512 = self.data.costs.fkind * self.data.costs.c22512

        #  Account 225.1.3 : TF coil dump resistors

        self.data.costs.c22513 = 1.0e-6 * (
            1.0e9
            * self.data.costs.UCTFDR
            * tfcoil_variables.e_tf_magnetic_stored_total_gj
            + self.data.costs.UCTFGR * 0.5e0 * tfcoil_variables.n_tf_coils
        )
        self.data.costs.c22513 = self.data.costs.fkind * self.data.costs.c22513

        #  Account 225.1.4 : TF coil instrumentation and control

        self.data.costs.c22514 = (
            1.0e-6 * self.data.costs.UCTFIC * (30.0e0 * tfcoil_variables.n_tf_coils)
        )
        self.data.costs.c22514 = self.data.costs.fkind * self.data.costs.c22514

        #  Account 225.1.5 : TF coil bussing

        if tfcoil_variables.i_tf_sup != TFConductorModel.SUPERCONDUCTING:
            self.data.costs.c22515 = (
                1.0e-6 * self.data.costs.uctfbus * tfcoil_variables.m_tf_bus
            )
        else:
            self.data.costs.c22515 = (
                1.0e-6
                * self.data.costs.ucbus
                * tfcoil_variables.c_tf_turn
                * tfcoil_variables.len_tf_bus
            )

        self.data.costs.c22515 = self.data.costs.fkind * self.data.costs.c22515

        #  Total account 225.1

        self.data.costs.c2251 = (
            self.data.costs.c22511
            + self.data.costs.c22512
            + self.data.costs.c22513
            + self.data.costs.c22514
            + self.data.costs.c22515
        )

    def acc2252(self):
        """Account 225.2 : PF coil power conditioning
        This routine evaluates the Account 225.2 (PF coil power
        conditioning) costs.
        Costs are taken from the equipment specification of the
        <A HREF="pfpwr.html">pfpwr</A> routine from the plant power module.
        """
        self.data.costs.c22521 = (
            1.0e-6 * self.data.costs.ucpfps * heat_transport_variables.peakmva
        )
        self.data.costs.c22521 = self.data.costs.fkind * self.data.costs.c22521

        #  Account 225.2.2 : PF coil instrumentation and control

        self.data.costs.c22522 = (
            1.0e-6 * self.data.costs.ucpfic * pf_power_variables.pfckts * 30.0e0
        )
        self.data.costs.c22522 = self.data.costs.fkind * self.data.costs.c22522

        #  Account 225.2.3 : PF coil bussing

        self.data.costs.c22523 = (
            1.0e-6
            * self.data.costs.ucpfb
            * pf_power_variables.spfbusl
            * pf_power_variables.acptmax
        )
        self.data.costs.c22523 = self.data.costs.fkind * self.data.costs.c22523

        #  Account 225.2.4 : PF coil burn power supplies

        if pf_power_variables.pfckts != 0.0e0:  # noqa: RUF069
            self.data.costs.c22524 = (
                1.0e-6
                * self.data.costs.ucpfbs
                * pf_power_variables.pfckts
                * (pf_power_variables.srcktpm / pf_power_variables.pfckts) ** 0.7e0
            )
        else:
            self.data.costs.c22524 = 0.0e0

        self.data.costs.c22524 = self.data.costs.fkind * self.data.costs.c22524

        #  Account 225.2.5 : PF coil breakers

        self.data.costs.c22525 = (
            1.0e-6
            * self.data.costs.ucpfbk
            * pf_power_variables.pfckts
            * (pf_power_variables.acptmax * pf_power_variables.vpfskv) ** 0.7e0
        )
        self.data.costs.c22525 = self.data.costs.fkind * self.data.costs.c22525

        #  Account 225.2.6 : PF coil dump resistors

        self.data.costs.c22526 = (
            1.0e-6 * self.data.costs.ucpfdr1 * pf_power_variables.ensxpfm
        )
        self.data.costs.c22526 = self.data.costs.fkind * self.data.costs.c22526

        #  Account 225.2.7 : PF coil AC breakers

        self.data.costs.c22527 = (
            1.0e-6 * self.data.costs.ucpfcb * pf_power_variables.pfckts
        )
        self.data.costs.c22527 = self.data.costs.fkind * self.data.costs.c22527

        #  Total account 225.2

        self.data.costs.c2252 = (
            self.data.costs.c22521
            + self.data.costs.c22522
            + self.data.costs.c22523
            + self.data.costs.c22524
            + self.data.costs.c22525
            + self.data.costs.c22526
            + self.data.costs.c22527
        )

    def acc226(self):
        """Account 226 : Heat transport system
        This routine evaluates the Account 226 (heat transport system) costs.
        Costs are estimated from major equipment and heat transport
        system loops developed in the heatpwr module of the code.
        """
        self.data.costs.c226 = (
            self.data.costs.c2261 + self.data.costs.c2262 + self.data.costs.c2263
        )

    def acc2261(self):
        """Account 2261 : Reactor cooling system
        This routine evaluates the Account 2261 -
        """
        cmlsa = [0.4000e0, 0.7000e0, 0.8500e0, 1.0000e0]
        exphts = 0.7e0

        #  Pumps and piping system
        #  N.B. with blktmodel > 0, the blanket is assumed to be helium-cooled,
        #  but the shield etc. is water-cooled (i_blkt_coolant_type=2). Therefore, a slight
        #  inconsistency exists here...
        self.data.costs.cpp = (
            1.0e-6
            * self.data.costs.uchts[fwbs_variables.i_blkt_coolant_type - 1]
            * (
                (1.0e6 * heat_transport_variables.p_fw_div_heat_deposited_mw) ** exphts
                + (1.0e6 * fwbs_variables.p_blkt_nuclear_heat_total_mw) ** exphts
                + (1.0e6 * fwbs_variables.p_shld_nuclear_heat_mw) ** exphts
            )
        )

        self.data.costs.cpp = (
            self.data.costs.fkind * self.data.costs.cpp * cmlsa[self.data.costs.lsa - 1]
        )

        #  Primary heat exchangers
        self.data.costs.chx = (
            1.0e-6
            * self.data.costs.UCPHX
            * heat_transport_variables.n_primary_heat_exchangers
            * (
                1.0e6
                * heat_transport_variables.p_plant_primary_heat_mw
                / heat_transport_variables.n_primary_heat_exchangers
            )
            ** exphts
        )
        self.data.costs.chx = (
            self.data.costs.fkind * self.data.costs.chx * cmlsa[self.data.costs.lsa - 1]
        )

        self.data.costs.c2261 = self.data.costs.chx + self.data.costs.cpp

    def acc2262(self):
        """Account 2262 : Auxiliary component cooling
        This routine evaluates the Account 2262 - Auxiliary component cooling
        """
        cmlsa = 0.4000e0, 0.7000e0, 0.8500e0, 1.0000e0
        exphts = 0.7e0

        #  Pumps and piping system
        self.data.costs.cppa = (
            1.0e-6
            * self.data.costs.UCAHTS
            * (
                (1.0e6 * heat_transport_variables.p_hcd_electric_loss_mw) ** exphts
                + (1.0e6 * heat_transport_variables.p_cryo_plant_electric_mw) ** exphts
                + (1.0e6 * heat_transport_variables.vachtmw) ** exphts
                + (1.0e6 * heat_transport_variables.p_tritium_plant_electric_mw)
                ** exphts
                + (1.0e6 * heat_transport_variables.fachtmw) ** exphts
            )
        )

        if ife_variables.ife == 1:
            self.data.costs.cppa += (
                1.0e-6
                * self.data.costs.UCAHTS
                * (
                    (1.0e6 * ife_variables.tdspmw) ** exphts
                    + (1.0e6 * ife_variables.tfacmw) ** exphts
                )
            )

        #  Apply Nth kind and safety assurance factors
        self.data.costs.cppa = (
            self.data.costs.fkind * self.data.costs.cppa * cmlsa[self.data.costs.lsa - 1]
        )

        self.data.costs.c2262 = self.data.costs.cppa

    def acc2263(self):
        """Account 2263 : Cryogenic system
        This routine evaluates the Account 2263 - Cryogenic system
        """
        cmlsa = 0.4000e0, 0.7000e0, 0.8500e0, 1.0000e0
        expcry = 0.67e0

        self.data.costs.c2263 = (
            1.0e-6
            * self.data.costs.uccry
            * 4.5e0
            / tfcoil_variables.temp_tf_cryo
            * heat_transport_variables.helpow**expcry
        )

        #  Apply Nth kind and safety factors
        self.data.costs.c2263 = (
            self.data.costs.fkind
            * self.data.costs.c2263
            * cmlsa[self.data.costs.lsa - 1]
        )

    def acc227(self):
        """Account 227 : Fuel handling
        This routine evaluates the Account 227 (fuel handling) costs.
        Costs are scaled from TETRA reactor code runs.
        """
        self.data.costs.c227 = (
            self.data.costs.c2271
            + self.data.costs.c2272
            + self.data.costs.c2273
            + self.data.costs.c2274
        )

    def acc2271(self):
        """Account 2271 : Fuelling system
        This routine evaluates the Account 2271 - Fuelling system
        """
        self.data.costs.c2271 = 1.0e-6 * self.data.costs.ucf1

        #  Apply Nth kind factor
        self.data.costs.c2271 = self.data.costs.fkind * self.data.costs.c2271

    def acc2272(self):
        """Account 2272 : Fuel processing and purification
        This routine evaluates the Account 2272 - Fuel processing
        """
        if ife_variables.ife != 1:
            #  Previous calculation, using molflow_plasma_fuelling_required in Amps:
            #  1.3 should have been physics_variables.m_fuel_amu*umass/electron_charge*1000*s/day = 2.2
            # wtgpd = burnup * molflow_plasma_fuelling_required * 1.3e0

            #  New calculation: 2 nuclei * reactions/sec * kg/nucleus * g/kg * sec/day
            physics_variables.wtgpd = (
                2.0e0
                * physics_variables.rndfuel
                * physics_variables.m_fuel_amu
                * constants.UMASS
                * 1000.0e0
                * 86400.0e0
            )
        else:
            targtm = (
                ife_variables.gain
                * ife_variables.edrive
                * 3.0e0
                * 1.67e-27
                * 1.0e3
                / (constants.ELECTRON_VOLT * 17.6e6 * ife_variables.fburn)
            )
            physics_variables.wtgpd = targtm * ife_variables.reprat * 86400.0e0

        #  Assumes that He3 costs same as tritium to process...
        self.data.costs.c2272 = (
            1.0e-6
            * self.data.costs.UCFPR
            * (0.5e0 + 0.5e0 * (physics_variables.wtgpd / 60.0e0) ** 0.67e0)
        )

        self.data.costs.c2272 = self.data.costs.fkind * self.data.costs.c2272

    def acc2273(self):
        """Account 2273 : Atmospheric recovery systems
        This routine evaluates the Account 2273 - Atmospheric recovery systems
        """
        cfrht = 1.0e5

        #  No detritiation needed if purely D-He3 reaction
        if physics_variables.f_plasma_fuel_tritium > 1.0e-3:
            self.data.costs.c2273 = (
                1.0e-6
                * self.data.costs.UCDTC
                * (
                    (cfrht / 1.0e4) ** 0.6e0
                    * (buildings_variables.volrci + buildings_variables.wsvol)
                )
            )
        else:
            self.data.costs.c2273 = 0.0e0

        self.data.costs.c2273 = self.data.costs.fkind * self.data.costs.c2273

    def acc2274(self):
        """Account 2274 : Nuclear building ventilation
        This routine evaluates the Account 2274 - Nuclear building ventilation
        """
        self.data.costs.c2274 = (
            1.0e-6
            * self.data.costs.UCNBV
            * (buildings_variables.volrci + buildings_variables.wsvol) ** 0.8e0
        )

        #  Apply Nth kind factor
        self.data.costs.c2274 = self.data.costs.fkind * self.data.costs.c2274

    def acc228(self):
        """Account 228 : Instrumentation and control

        This routine evaluates the Account 228 (instrumentation and
        control) costs.
        Costs are based on TFCX and INTOR.
        """
        self.data.costs.c228 = 1.0e-6 * self.data.costs.uciac
        self.data.costs.c228 = self.data.costs.fkind * self.data.costs.c228

    def acc229(self):
        """Account 229 : Maintenance equipment

        This routine evaluates the Account 229 (maintenance equipment) costs.
        """
        self.data.costs.c229 = 1.0e-6 * self.data.costs.ucme
        self.data.costs.c229 = self.data.costs.fkind * self.data.costs.c229

    def acc23(self):
        """Account 23 : Turbine plant equipment

        This routine evaluates the Account 23 (turbine plant equipment) costs.
        """
        exptpe = 0.83e0
        if self.data.costs.ireactor == 1:
            self.data.costs.c23 = (
                1.0e-6
                * self.data.costs.ucturb[fwbs_variables.i_blkt_coolant_type - 1]
                * (heat_transport_variables.p_plant_electric_gross_mw / 1200.0e0)
                ** exptpe
            )

    def acc24(self):
        """Account 24 : Electric plant equipment

        This routine evaluates the Account 24 (electric plant equipment) costs.
        """
        self.data.costs.c24 = (
            self.data.costs.c241
            + self.data.costs.c242
            + self.data.costs.c243
            + self.data.costs.c244
            + self.data.costs.c245
        )

    def acc241(self):
        """Account 241 : Electric plant equipment - switchyard
        This routine evaluates the Account 241 - switchyard
        """
        cmlsa = 0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0

        #  Account 241 : Switchyard
        self.data.costs.c241 = (
            1.0e-6 * self.data.costs.UCSWYD * cmlsa[self.data.costs.lsa - 1]
        )

    def acc242(self):
        """Account 242 : Electric plant equipment - Transformers
        This routine evaluates the Account 242 - Transformers
        """
        cmlsa = 0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0
        expepe = 0.9e0

        #  Account 242 : Transformers
        self.data.costs.c242 = 1.0e-6 * (
            self.data.costs.UCPP * (heat_transport_variables.pacpmw * 1.0e3) ** expepe
            + self.data.costs.UCAP
            * (heat_transport_variables.p_plant_electric_base_total_mw * 1.0e3)
        )

        #  Apply safety assurance factor
        self.data.costs.c242 *= cmlsa[self.data.costs.lsa - 1]

    def acc243(self):
        """Account 243 : Electric plant equipment - Low voltage
        This routine evaluates the Account 243 - Low voltage
        """
        cmlsa = 0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0

        #  Account 243 : Low voltage
        #  (include 0.8 factor for transformer efficiency)
        self.data.costs.c243 = (
            1.0e-6
            * self.data.costs.UCLV
            * heat_transport_variables.tlvpmw
            * 1.0e3
            / 0.8e0
            * cmlsa[self.data.costs.lsa - 1]
        )

    def acc244(self):
        """Account 244 : Electric plant equipment - Diesel generators
        This routine evaluates the Account 244 - Diesel generators
        """
        cmlsa = [0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0]

        #  Account 244 : Diesel generator (8 MW per generator,  assume 4 )
        self.data.costs.c244 = (
            1.0e-6 * self.data.costs.UCDGEN * 4.0e0 * cmlsa[self.data.costs.lsa - 1]
        )

    def acc245(self):
        """Account 245 : Electric plant equipment - Aux facility power
        This routine evaluates the Account 245 - Aux facility power
        """
        cmlsa = 0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0

        #  Account 245 : Auxiliary facility power needs
        self.data.costs.c245 = (
            1.0e-6 * self.data.costs.UCAF * cmlsa[self.data.costs.lsa - 1]
        )

    def acc25(self):
        """Account 25 : Miscellaneous plant equipment

        This routine evaluates the Account 25 (miscellaneous plant
        equipment) costs, such as waste treatment.
        """
        cmlsa = 0.7700e0, 0.8850e0, 0.9425e0, 1.0000e0

        self.data.costs.c25 = (
            1.0e-6 * self.data.costs.ucmisc * cmlsa[self.data.costs.lsa - 1]
        )

    def acc26(self):
        """Account 26 : Heat rejection system

        This routine evaluates the Account 26 (heat rejection system) costs.
        Costs are scaled with the total plant heat rejection based on
        commercial systems.
        J. Delene, private communication, ORNL, June 1990
        """
        cmlsa = [0.8000e0, 0.9000e0, 0.9500e0, 1.0000e0]

        # Calculate rejected heat for non-reactor (==0) and reactor (==1)
        if self.data.costs.ireactor == 0:
            pwrrej = (
                physics_variables.p_fusion_total_mw
                + heat_transport_variables.p_hcd_electric_total_mw
                + tfcoil_variables.tfcmw
            )
        else:
            pwrrej = (
                heat_transport_variables.p_plant_primary_heat_mw
                - heat_transport_variables.p_plant_electric_gross_mw
            )

        # self.data.costs.uchrs - reference cost of heat rejection system [$]
        self.data.costs.c26 = (
            1.0e-6
            * self.data.costs.uchrs
            * pwrrej
            / 2300.0e0
            * cmlsa[self.data.costs.lsa - 1]
        )

    def acc9(self):
        """Account 9 : Indirect cost and contingency allowances

        This routine evaluates the Account 9 (indirect cost and
        contingency allowances) costs.
        The cost modelling is based on the commercial plant model of a
        single contractor performing all plant engineering and construction
        management, using commercially purchased equipment and materials.
        <P>The project contingency is an allowance for incomplete design
        specification and unforeseen events during the plant construction.
        <P>The factors used are estimated from commercial plant experience.
        J. Delene, private communication, ORNL, June 1990
        """
        self.data.costs.cindrt = (
            self.data.costs.cfind[self.data.costs.lsa - 1]
            * self.data.costs.cdirt
            * (1.0e0 + self.data.costs.cowner)
        )

        #  Contingency costs

        self.data.costs.ccont = self.data.costs.fcontng * (
            self.data.costs.cdirt + self.data.costs.cindrt
        )

    def acc2253(self):
        """Account 225.3 : Energy storage
        This routine evaluates the Account 225.3 (energy storage) costs.
        """
        self.data.costs.c2253 = 0.0e0

        #  Thermal storage options for a pulsed reactor
        #  See F/MPE/MOD/CAG/PROCESS/PULSE/0008 and 0014

        if pulse_variables.i_pulsed_plant == 1:
            if pulse_variables.istore == 1:
                #  Option 1 from ELECTROWATT report
                #  Pulsed Fusion Reactor Study : AEA FUS 205

                #  Increased condensate tank capacity
                self.data.costs.c2253 = 0.1e0

                #  Additional electrically-driven feedpump (50 per cent duty)
                self.data.costs.c2253 += 0.8e0

                #  Increased turbine-generator duty (5 per cent duty)
                self.data.costs.c2253 += 4.0e0

                #  Additional auxiliary transformer capacity and ancillaries
                self.data.costs.c2253 += 0.5e0

                #  Increased drum capacity
                self.data.costs.c2253 += 2.8e0

                #  Externally fired superheater
                self.data.costs.c2253 += 29.0e0

            elif pulse_variables.istore == 2:
                #  Option 2 from ELECTROWATT report
                #  Pulsed Fusion Reactor Study : AEA FUS 205

                #  Increased condensate tank capacity
                self.data.costs.c2253 = 0.1e0

                #  Additional electrically-driven feedpump (50 per cent duty)
                self.data.costs.c2253 += 0.8e0

                #  Increased drum capacity
                self.data.costs.c2253 += 2.8e0

                #  Increased turbine-generator duty (5 per cent duty)
                self.data.costs.c2253 += 4.0e0

                #  Additional fired boiler (1 x 100 per cent duty)
                self.data.costs.c2253 += 330.0e0

                #  HP/LP steam bypass system for auxiliary boiler
                #  (30 per cent boiler capacity)
                self.data.costs.c2253 += 1.0e0

                #  Dump condenser
                self.data.costs.c2253 += 2.0e0

                #  Increased cooling water system capacity
                self.data.costs.c2253 += 18.0e0

            elif pulse_variables.istore == 3:
                #  Simplistic approach that assumes that a large stainless steel
                #  block acts as the thermal storage medium. No account is taken
                #  of the cost of the piping within the block, etc.
                #
                #  shcss is the specific heat capacity of stainless steel (J/kg/K)
                #  pulse_variables.dtstor is the maximum allowable temperature change in the
                #  stainless steel block (input)

                shcss = 520.0e0
                self.data.costs.c2253 = (
                    self.data.costs.ucblss
                    * (heat_transport_variables.p_plant_primary_heat_mw * 1.0e6)
                    * times_variables.t_plant_pulse_no_burn
                    / (shcss * pulse_variables.dtstor)
                )

            else:
                raise ProcessValueError(
                    "Illegal value for istore", istore=pulse_variables.istore
                )

        if pulse_variables.istore < 3:
            #  Scale self.data.costs.c2253 with net electric power

            self.data.costs.c2253 = (
                self.data.costs.c2253
                * heat_transport_variables.p_plant_electric_net_mw
                / 1200.0e0
            )

            #  It is necessary to convert from 1992 pounds to 1990 dollars
            #  Reasonable guess for the exchange rate + inflation factor
            #  inflation = 5% per annum; exchange rate = 1.5 dollars per pound

            self.data.costs.c2253 *= 1.36e0

        self.data.costs.c2253 = self.data.costs.fkind * self.data.costs.c2253

    def coelc(self):
        """Routine to calculate the cost of electricity for a fusion power plant

        outfile : input integer : output file unit

        This routine performs the calculation of the cost of electricity
        for a fusion power plant.
        <P>Annual costs are in megadollars/year, electricity costs are in
        millidollars/kWh, while other costs are in megadollars.
        All values are based on 1990 dollars.
        """
        if ife_variables.ife == 1:
            kwhpy = (
                1.0e3
                * heat_transport_variables.p_plant_electric_net_mw
                * (24.0e0 * constants.N_DAY_YEAR)
                * self.data.costs.f_t_plant_available
            )
        else:
            kwhpy = (
                1.0e3
                * heat_transport_variables.p_plant_electric_net_mw
                * (24.0e0 * constants.N_DAY_YEAR)
                * self.data.costs.f_t_plant_available
                * times_variables.t_plant_pulse_burn
                / times_variables.t_plant_pulse_total
            )

        #  Costs due to reactor plant
        #  ==========================

        #  Interest on construction costs

        self.data.costs.moneyint = self.data.costs.concost * (
            self.data.costs.fcap0 - 1.0e0
        )

        #  Capital costs

        self.data.costs.capcost = self.data.costs.concost + self.data.costs.moneyint

        #  Annual cost of plant capital cost

        anncap = self.data.costs.capcost * self.data.costs.fcr0

        # SJP Issue #836
        # Check for the condition when kwhpy=0

        kwhpy = max(kwhpy, 1.0e-10)

        #  Cost of electricity due to plant capital cost

        self.data.costs.coecap = 1.0e9 * anncap / kwhpy

        #  Costs due to first wall and blanket renewal
        #  ===========================================

        #  Compound interest factor

        feffwbl = (1.0e0 + self.data.costs.discount_rate) ** fwbs_variables.life_blkt

        #  Capital recovery factor

        crffwbl = (feffwbl * self.data.costs.discount_rate) / (feffwbl - 1.0e0)

        #  Annual cost of replacements

        annfwbl = (
            (self.data.costs.fwallcst + self.data.costs.blkcst)
            * (1.0e0 + self.data.costs.cfind[self.data.costs.lsa - 1])
            * self.data.costs.fcap0cp
            * crffwbl
        )

        if self.data.costs.ifueltyp == 2:
            annfwbl *= 1.0e0 - fwbs_variables.life_blkt_fpy / self.data.costs.life_plant

        #  Cost of electricity due to first wall/blanket replacements

        coefwbl = 1.0e9 * annfwbl / kwhpy

        #  Costs due to divertor renewal
        #  =============================

        if ife_variables.ife == 1:
            anndiv = 0.0e0
            coediv = 0.0e0
        else:
            #  Compound interest factor

            fefdiv = (1.0e0 + self.data.costs.discount_rate) ** self.data.costs.life_div

            #  Capital recovery factor

            crfdiv = (fefdiv * self.data.costs.discount_rate) / (fefdiv - 1.0e0)

            #  Annual cost of replacements

            anndiv = (
                self.data.costs.divcst
                * (1.0e0 + self.data.costs.cfind[self.data.costs.lsa - 1])
                * self.data.costs.fcap0cp
                * crfdiv
            )

            #  Cost of electricity due to divertor replacements

            if self.data.costs.ifueltyp == 2:
                anndiv *= (
                    1.0e0 - self.data.costs.life_div_fpy / self.data.costs.life_plant
                )

            coediv = 1.0e9 * anndiv / kwhpy

        #  Costs due to centrepost renewal
        #  ===============================

        if (physics_variables.itart == 1) and (ife_variables.ife != 1):
            #  Compound interest factor

            fefcp = (1.0e0 + self.data.costs.discount_rate) ** self.data.costs.cplife_cal

            #  Capital recovery factor

            crfcp = (fefcp * self.data.costs.discount_rate) / (fefcp - 1.0e0)

            #  Annual cost of replacements

            anncp = (
                self.data.costs.cpstcst
                * (1.0e0 + self.data.costs.cfind[self.data.costs.lsa - 1])
                * self.data.costs.fcap0cp
                * crfcp
            )

            #  Cost of electricity due to centrepost replacements
            if self.data.costs.ifueltyp == 2:
                anncp *= 1.0e0 - self.data.costs.cplife / self.data.costs.life_plant

            coecp = 1.0e9 * anncp / kwhpy

        else:
            anncp = 0.0e0
            coecp = 0.0e0

        #  Costs due to partial current drive system renewal
        #  =================================================

        #  Compound interest factor

        fefcdr = (1.0e0 + self.data.costs.discount_rate) ** self.data.costs.cdrlife_cal

        #  Capital recovery factor

        crfcdr = (fefcdr * self.data.costs.discount_rate) / (fefcdr - 1.0e0)

        #  Annual cost of replacements

        if self.data.costs.ifueltyp == 0:
            anncdr = 0.0e0
        else:
            anncdr = (
                self.data.costs.cdcost
                * self.data.costs.fcdfuel
                / (1.0e0 - self.data.costs.fcdfuel)
                * (1.0e0 + self.data.costs.cfind[self.data.costs.lsa - 1])
                * self.data.costs.fcap0cp
                * crfcdr
            )

        #  Cost of electricity due to current drive system replacements

        coecdr = 1.0e9 * anncdr / kwhpy

        #  Costs due to operation and maintenance
        #  ======================================

        #  Annual cost of operation and maintenance

        if heat_transport_variables.p_plant_electric_net_mw < 0:
            sqrt_p_plant_electric_net_mw_1200 = 0.0
            logger.warning(
                "p_plant_electric_net_mw has gone negative! Clamping it to 0 for the calculation of annoam and annwst (cost of maintenance and cost of waste)."
            )
        else:
            sqrt_p_plant_electric_net_mw_1200 = np.sqrt(
                heat_transport_variables.p_plant_electric_net_mw / 1200.0e0
            )
        annoam = (
            self.data.costs.ucoam[self.data.costs.lsa - 1]
            * sqrt_p_plant_electric_net_mw_1200
        )

        #  Additional cost due to pulsed reactor thermal storage
        #  See F/MPE/MOD/CAG/PROCESS/PULSE/0008
        #
        #      if (i_pulsed_plant.eq.1) :
        #         if (istore.eq.1) :
        #            annoam1 = 51.0e0
        #         elif  (istore.eq.2) :
        #            annoam1 = 22.2e0
        #         else:
        #            continue
        #
        #
        #  Scale with net electric power
        #
        #         annoam1 = annoam1 * heat_transport_variables.p_plant_electric_net_mw/1200.0e0
        #
        #  It is necessary to convert from 1992 pounds to 1990 dollars
        #  Reasonable guess for the exchange rate + inflation factor
        #  inflation = 5% per annum; exchange rate = 1.5 dollars per pound
        #
        #         annoam1 = annoam1 * 1.36e0
        #
        #         annoam = annoam + annoam1
        #
        #

        #  Cost of electricity due to operation and maintenance

        self.data.costs.coeoam = 1.0e9 * annoam / kwhpy

        #  Costs due to reactor fuel
        #  =========================

        #  Annual cost of fuel

        if ife_variables.ife != 1:
            #  Sum D-T fuel cost and He3 fuel cost
            annfuel = (
                self.data.costs.ucfuel
                * heat_transport_variables.p_plant_electric_net_mw
                / 1200.0e0
                + 1.0e-6
                * physics_variables.f_plasma_fuel_helium3
                * physics_variables.wtgpd
                * 1.0e-3
                * self.data.costs.uche3
                * constants.N_DAY_YEAR
                * self.data.costs.f_t_plant_available
            )
        else:
            annfuel = (
                1.0e-6
                * ife_variables.uctarg
                * ife_variables.reprat
                * 3.1536e7
                * self.data.costs.f_t_plant_available
            )

        #  Cost of electricity due to reactor fuel

        coefuel = 1.0e9 * annfuel / kwhpy

        #  Costs due to waste disposal
        #  ===========================

        #  Annual cost of waste disposal

        annwst = (
            self.data.costs.ucwst[self.data.costs.lsa - 1]
            * sqrt_p_plant_electric_net_mw_1200
        )

        #  Cost of electricity due to waste disposal

        coewst = 1.0e9 * annwst / kwhpy

        #  Costs due to decommissioning fund
        #  =================================

        #  Annual contributions to fund for decommissioning
        #  A fraction self.data.costs.decomf of the construction cost is set aside for
        #  this purpose at the start of the plant life.
        #  Final factor takes into account inflation over the plant lifetime
        #  (suggested by Tim Hender 07/03/96)
        #  Difference (self.data.costs.dintrt) between borrowing and saving interest rates is
        #  included, along with the possibility of completing the fund self.data.costs.dtlife
        #  years before the end of the plant's lifetime

        anndecom = (
            self.data.costs.decomf
            * self.data.costs.concost
            * self.data.costs.fcr0
            / (1.0e0 + self.data.costs.discount_rate - self.data.costs.dintrt)
            ** (self.data.costs.life_plant - self.data.costs.dtlife)
        )

        #  Cost of electricity due to decommissioning fund

        coedecom = 1.0e9 * anndecom / kwhpy

        #  Total costs
        #  ===========

        #  Annual costs due to 'fuel-like' components

        # annfuelt = annfwbl + anndiv + anncdr + anncp + annfuel + annwst

        #  Total cost of electricity due to 'fuel-like' components

        self.data.costs.coefuelt = coefwbl + coediv + coecdr + coecp + coefuel + coewst

        #  Total annual costs

        # anntot = anncap + annfuelt + annoam + anndecom

        #  Total cost of electricity

        self.data.costs.coe = (
            self.data.costs.coecap
            + self.data.costs.coefuelt
            + self.data.costs.coeoam
            + coedecom
        )

    def convert_fpy_to_calendar(self):
        """Routine to convert component lifetimes in FPY to calendar years.
        Required for replacement component costs.
        """
        # FW/Blanket and HCD
        if fwbs_variables.life_blkt_fpy < self.data.costs.life_plant:
            fwbs_variables.life_blkt = (
                fwbs_variables.life_blkt_fpy * self.data.costs.f_t_plant_available
            )
            # Current drive system lifetime (assumed equal to first wall and blanket lifetime)
            self.data.costs.cdrlife_cal = fwbs_variables.life_blkt
        else:
            fwbs_variables.life_blkt = fwbs_variables.life_blkt_fpy

        # Divertor
        if self.data.costs.life_div_fpy < self.data.costs.life_plant:
            self.data.costs.life_div = (
                self.data.costs.life_div_fpy * self.data.costs.f_t_plant_available
            )
        else:
            self.data.costs.life_div = self.data.costs.life_div_fpy

        # Centrepost
        if physics_variables.itart == 1:
            if self.data.costs.cplife < self.data.costs.life_plant:
                self.data.costs.cplife_cal = (
                    self.data.costs.cplife * self.data.costs.f_t_plant_available
                )
            else:
                self.data.costs.cplife_cal = self.data.costs.cplife
