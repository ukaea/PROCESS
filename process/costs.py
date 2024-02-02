from process.fortran import constants, cost_variables
from process.fortran import process_output as po
from process.fortran import ife_variables, fwbs_variables
from process.fortran import tfcoil_variables
from process.fortran import physics_variables
from process.fortran import buildings_variables
from process.fortran import build_variables
from process.fortran import structure_variables
from process.fortran import divertor_variables
from process.fortran import pfcoil_variables
from process.fortran import current_drive_variables
from process.fortran import vacuum_variables
from process.fortran import heat_transport_variables
from process.fortran import pf_power_variables
from process.fortran import pulse_variables
from process.fortran import times_variables
from process.fortran import error_handling
from process.variables import AnnotatedVariable
import numpy


class Costs:
    def __init__(self):
        self.outfile = constants.nout

        # Various cost account values (M$)
        self.c228 = AnnotatedVariable(
            float, 0.0, docstring="c228 account cost", units="M$"
        )
        self.c229 = AnnotatedVariable(
            float, 0.0, docstring="c229 account cost", units="M$"
        )
        self.c23 = AnnotatedVariable(
            float, 0.0, docstring="c23 account cost", units="M$"
        )
        self.c25 = AnnotatedVariable(
            float, 0.0, docstring="c25 account cost", units="M$"
        )
        self.c26 = AnnotatedVariable(
            float, 0.0, docstring="c26 account cost", units="M$"
        )
        self.cindrt = AnnotatedVariable(
            float, 0.0, docstring="cindrt account cost", units="M$"
        )
        self.ccont = AnnotatedVariable(
            float, 0.0, docstring="ccont account cost", units="M$"
        )

        for names in (
            (f"c226{no}" for no in ("", 1, 2, 3)),  # Accnt 226: Heat transport system
            (f"c227{no}" for no in ("", 1, 2, 3, 4)),  # Accnt 227: Fuel handling
            (f"c24{no}" for no in ("", 1, 2, 3, 4, 5)),  # Accnt 24: elec equipment
            (f"c21{no}" for no in ("", 1, 2, 3, 4, 41, 42, 5, 6, 7, 71, 72, 73, 74)),
            ("c22",),
            (f"c221{no}" for no in (1, 2, 21, 22, 23, 24, 25, 26, 27, 3, 31, 32, 4, 5)),
            (f"c222{no}" for no in (1, 11, 12, 13, 14, 15, 2, 21, 22, 23, 24, 3)),
            (f"c223{no}" for no in ("", 1, 2, 3, 4)),
            (f"c224{no}" for no in ("", 1, 2, 3, 4, 5, 6)),
            (
                f"c225{no}"
                for no in ("", 1, 11, 12, 13, 14, 15, 2, 21, 22, 23, 24, 25, 26, 27, 3)
            ),
            ("chx", "cpp", "cppa", "c22128"),
        ):
            for i in names:
                setattr(self, i, AnnotatedVariable(float, 0.0, docstring="", units=""))

    def run(self):
        """
        Cost accounting for a fusion power plant
        author: P J Knight, CCFE, Culham Science Centre

        This routine performs the cost accounting for a fusion power plant.
        The direct costs are calculated based on parameters input
        from other sections of the code.
        <P>Costs are in 1990 $, and assume first-of-a-kind components
        unless otherwise stated. Account 22 costs include a multiplier
        to account for Nth-of-a-kind cost reductions.
        <P>The code is arranged in the order of the standard accounts.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
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
        # cdirt = c21 + c22 + self.c23 + self.c24 + self.c25 + self.c26 + chplant
        cost_variables.cdirt = (
            self.c21 + self.c22 + self.c23 + self.c24 + self.c25 + self.c26
        )

        #  Account 9 : Indirect cost and project contingency
        self.acc9()

        #  Constructed cost
        cost_variables.concost = cost_variables.cdirt + self.cindrt + self.ccont

        #  Cost of electricity
        if (cost_variables.ireactor == 1) and (cost_variables.ipnet == 0):
            self.coelc()

    def output(self):
        if cost_variables.output_costs == 0:
            return

        po.oheadr(self.outfile, "Power Reactor Costs (1990 US$)")

        po.ovarrf(
            self.outfile,
            "First wall / blanket life (years)",
            "(fwbllife)",
            fwbs_variables.bktlife,
        )

        if ife_variables.ife != 1:
            po.ovarrf(
                self.outfile,
                "Divertor life (years)",
                "(divlife.)",
                cost_variables.divlife,
            )
            if physics_variables.itart == 1:
                po.ovarrf(
                    self.outfile,
                    "Centrepost life (years)",
                    "(cplife.)",
                    cost_variables.cplife,
                )

        po.ovarrf(
            self.outfile, "Cost of electricity (m$/kWh)", "(coe)", cost_variables.coe
        )

        po.osubhd(self.outfile, "Power Generation Costs :")
        # TODO: Convert fortran format to Python
        # if ((annfwbl != annfwbl) or (annfwbl > 1.0e10) or (annfwbl < 0.0e0)) :
        #     write(outfile,*)'Problem with annfwbl'
        #     write(outfile,*)'fwallcst=', fwallcst, '  blkcst=', cost_variables.blkcst
        #     write(outfile,*)'crffwbl=', crffwbl,   '  fcap0cp=', cost_variables.fcap0cp
        #     write(outfile,*)'feffwbl=', feffwbl,   '  fwbllife=', fwbllife

        #       write(outfile,200) #          anncap,coecap, #          annoam,coeoam, #          anndecom,coedecom, #          annfwbl,coefwbl, #          anndiv,coediv, #          anncp,coecp, #          anncdr,coecdr, #          annfuel,coefuel, #          annwst,coewst, #          annfuelt,coefuelt, #          anntot,coe

        # 200   format( #          t76,'Annual Costs, M$       COE, m$/kWh'// #          1x,'Capital Investment',t80,f10.2,10x,f10.2/ #          1x,'Operation & Maintenance',t80,f10.2,10x,f10.2/ #          1x,'Decommissioning Fund',t80,f10.2,10x,f10.2/ #          1x,'Fuel Charge Breakdown'// #          5x,'Blanket & first wall',t72,f10.2,10x,f10.2/ #          5x,'Divertors',t72,f10.2,10x,f10.2/ #          5x,'Centrepost (TART only)',t72,f10.2,10x,f10.2/ #          5x,'Auxiliary Heating',t72,f10.2,10x,f10.2/ #          5x,'Actual Fuel',t72,f10.2,10x,f10.2/ #          5x,'Waste Disposal',t72,f10.2,10x,f10.2/ #          1x,'Total Fuel Cost',t80,f10.2,10x,f10.2// #          1x,'Total Cost',t80,f10.2,10x,f10.2 )

        if cost_variables.ifueltyp == 1:
            po.oshead(self.outfile, "Replaceable Components Direct Capital Cost")
            po.ovarrf(
                self.outfile,
                "First wall direct capital cost (M$)",
                "(fwallcst)",
                cost_variables.fwallcst,
            )
            po.ovarrf(
                self.outfile,
                "Blanket direct capital cost (M$)",
                "(blkcst)",
                cost_variables.blkcst,
            )
            if ife_variables.ife != 1:
                po.ovarrf(
                    self.outfile,
                    "Divertor direct capital cost (M$)",
                    "(divcst)",
                    cost_variables.divcst,
                )
                if physics_variables.itart == 1:
                    po.ovarrf(
                        self.outfile,
                        "Centrepost direct capital cost (M$)",
                        "(cpstcst)",
                        cost_variables.cpstcst,
                    )

                po.ovarrf(
                    self.outfile,
                    "Plasma heating/CD system cap cost (M$)",
                    "",
                    cost_variables.cdcost
                    * cost_variables.fcdfuel
                    / (1.0e0 - cost_variables.fcdfuel),
                )
                po.ovarrf(
                    self.outfile,
                    "Fraction of CD cost --> fuel cost",
                    "(fcdfuel)",
                    cost_variables.fcdfuel,
                )
            else:
                po.ovarrf(
                    self.outfile,
                    "IFE driver system direct cap cost (M$)",
                    "",
                    cost_variables.cdcost
                    * cost_variables.fcdfuel
                    / (1.0e0 - cost_variables.fcdfuel),
                )
                po.ovarrf(
                    self.outfile,
                    "Fraction of driver cost --> fuel cost",
                    "(fcdfuel)",
                    cost_variables.fcdfuel,
                )

        po.oheadr(self.outfile, "Detailed Costings (1990 US$)")
        po.ovarre(
            self.outfile,
            "Acc.22 multiplier for Nth of a kind",
            "(fkind)",
            cost_variables.fkind,
        )
        po.ovarin(
            self.outfile, "Level of Safety Assurance", "(lsa)", cost_variables.lsa
        )
        po.oblnkl(self.outfile)
        po.oshead(self.outfile, "Structures and Site Facilities")
        po.ocosts(
            self.outfile,
            "(c211)",
            "Site improvements, facilities, land (M$)",
            self.c211,
        )
        po.ocosts(self.outfile, "(c212)", "Reactor building cost (M$)", self.c212)
        po.ocosts(self.outfile, "(c213)", "Turbine building cost (M$)", self.c213)
        po.ocosts(
            self.outfile,
            "(c2141)",
            "Reactor maintenance building cost (M$)",
            self.c2141,
        )
        po.ocosts(self.outfile, "(c2142)", "Warm shop cost (M$)", self.c2142)
        po.ocosts(self.outfile, "(c215)", "Tritium building cost (M$)", self.c215)
        po.ocosts(
            self.outfile,
            "(c216)",
            "Electrical equipment building cost (M$)",
            self.c216,
        )
        po.ocosts(
            self.outfile,
            "(c2171)",
            "Additional buildings cost (M$)",
            self.c2171,
        )
        po.ocosts(
            self.outfile,
            "(c2172)",
            "Control room buildings cost (M$)",
            self.c2172,
        )
        po.ocosts(
            self.outfile,
            "(c2173)",
            "Shop and warehouses cost (M$)",
            self.c2173,
        )
        po.ocosts(
            self.outfile,
            "(c2174)",
            "Cryogenic building cost (M$)",
            self.c2174,
        )
        po.oblnkl(self.outfile)
        po.ocosts(self.outfile, "(c21)", "Total account 21 cost (M$)", self.c21)

        po.oshead(self.outfile, "Reactor Systems")
        po.ocosts(self.outfile, "(c2211)", "First wall cost (M$)", self.c2211)
        if ife_variables.ife != 1:
            if fwbs_variables.iblanket == 4:
                po.ocosts(
                    self.outfile,
                    "(c22121)",
                    "Blanket lithium-lead cost (M$)",
                    self.c22121,
                )
                po.ocosts(
                    self.outfile,
                    "(c22122)",
                    "Blanket lithium cost (M$)",
                    self.c22122,
                )
            else:
                po.ocosts(
                    self.outfile,
                    "(c22121)",
                    "Blanket beryllium cost (M$)",
                    self.c22121,
                )
                po.ocosts(
                    self.outfile,
                    "(c22122)",
                    "Blanket breeder material cost (M$)",
                    self.c22122,
                )

            po.ocosts(
                self.outfile,
                "(c22123)",
                "Blanket stainless steel cost (M$)",
                self.c22123,
            )
            po.ocosts(
                self.outfile,
                "(c22124)",
                "Blanket vanadium cost (M$)",
                self.c22124,
            )
        else:  # IFE
            po.ocosts(
                self.outfile,
                "(c22121)",
                "Blanket beryllium cost (M$)",
                self.c22121,
            )
            po.ocosts(
                self.outfile,
                "(c22122)",
                "Blanket lithium oxide cost (M$)",
                self.c22122,
            )
            po.ocosts(
                self.outfile,
                "(c22123)",
                "Blanket stainless steel cost (M$)",
                self.c22123,
            )
            po.ocosts(
                self.outfile,
                "(c22124)",
                "Blanket vanadium cost (M$)",
                self.c22124,
            )
            po.ocosts(
                self.outfile,
                "(c22125)",
                "Blanket carbon cloth cost (M$)",
                self.c22125,
            )
            po.ocosts(
                self.outfile,
                "(c22126)",
                "Blanket concrete cost (M$)",
                self.c22126,
            )
            po.ocosts(
                self.outfile,
                "(c22127)",
                "Blanket FLiBe cost (M$)",
                self.c22127,
            )
            po.ocosts(
                self.outfile,
                "(c22128)",
                "Blanket lithium cost (M$)",
                self.c22128,
            )

        po.ocosts(self.outfile, "(c2212)", "Blanket total cost (M$)", self.c2212)
        po.ocosts(self.outfile, "(c22131)", "Bulk shield cost (M$)", self.c22131)
        po.ocosts(
            self.outfile,
            "(c22132)",
            "Penetration shielding cost (M$)",
            self.c22132,
        )
        po.ocosts(self.outfile, "(c2213)", "Total shield cost (M$)", self.c2213)
        po.ocosts(
            self.outfile,
            "(c2214)",
            "Total support structure cost (M$)",
            self.c2214,
        )
        po.ocosts(self.outfile, "(c2215)", "Divertor cost (M$)", self.c2215)
        # TODO: Convert fortran format to Python
        #     if (cost_variables.ifueltyp == 1) :
        #         po.oblnkl(self.outfile)
        #         write(self.outfile,20)
        #     20     format(t2,             'First wall, total blanket and divertor direct costs',/,             t2,'are zero as they are assumed to be fuel costs.')
        #     elif  (cost_variables.ifueltyp == 2) :
        #         po.oblnkl(self.outfile)
        #         write(self.outfile,31)
        # 21     format(t2,             'Initial First wall, total blanket and divertor direct costs',/,             t2,'are in capital and replacemnet are in cost of electricity')

        po.oblnkl(self.outfile)
        po.ocosts(
            self.outfile,
            "(c221)",
            "Total account 221 cost (M$)",
            cost_variables.c221,
        )

        if ife_variables.ife != 1:
            po.oshead(self.outfile, "Magnets")

            if tfcoil_variables.i_tf_sup != 1:  # Resistive TF coils
                if physics_variables.itart == 1:
                    po.ocosts(
                        self.outfile,
                        "(c22211)",
                        "Centrepost costs (M$)",
                        self.c22211,
                    )
                else:
                    po.ocosts(
                        self.outfile,
                        "(c22211)",
                        "Inboard leg cost (M$)",
                        self.c22211,
                    )

                po.ocosts(
                    self.outfile,
                    "(c22212)",
                    "Outboard leg cost (M$)",
                    self.c22212,
                )
                po.ocosts(
                    self.outfile,
                    "(c2221)",
                    "TF magnet assemblies cost (M$)",
                    self.c2221,
                )
            else:  # Superconducting TF coils
                po.ocosts(
                    self.outfile,
                    "(c22211)",
                    "TF coil conductor cost (M$)",
                    self.c22211,
                )
                po.ocosts(
                    self.outfile,
                    "(c22212)",
                    "TF coil winding cost (M$)",
                    self.c22212,
                )
                po.ocosts(
                    self.outfile,
                    "(c22213)",
                    "TF coil case cost (M$)",
                    self.c22213,
                )
                po.ocosts(
                    self.outfile,
                    "(c22214)",
                    "TF intercoil structure cost (M$)",
                    self.c22214,
                )
                po.ocosts(
                    self.outfile,
                    "(c22215)",
                    "TF coil gravity support structure (M$)",
                    self.c22215,
                )
                po.ocosts(
                    self.outfile,
                    "(c2221)",
                    "TF magnet assemblies cost (M$)",
                    self.c2221,
                )

            po.ocosts(
                self.outfile,
                "(c22221)",
                "PF coil conductor cost (M$)",
                self.c22221,
            )
            po.ocosts(
                self.outfile,
                "(c22222)",
                "PF coil winding cost (M$)",
                self.c22222,
            )
            po.ocosts(
                self.outfile,
                "(c22223)",
                "PF coil case cost (M$)",
                self.c22223,
            )
            po.ocosts(
                self.outfile,
                "(c22224)",
                "PF coil support structure cost (M$)",
                self.c22224,
            )
            po.ocosts(
                self.outfile,
                "(c2222)",
                "PF magnet assemblies cost (M$)",
                self.c2222,
            )
            po.ocosts(
                self.outfile,
                "(c2223)",
                "Vacuum vessel assembly cost (M$)",
                self.c2223,
            )
            # TODO: Convert fortran format to Python
            #     if ((physics_variables.itart == 1)and(cost_variables.ifueltyp == 1)) :
            #         po.oblnkl(self.outfile)
            #         write(self.outfile,30)
            # 30        format(t2,                'Centrepost direct cost is zero, as it ',                'is assumed to be a fuel cost.')
            #     elif  ((physics_variables.itart == 1)and(cost_variables.ifueltyp == 2)) :
            #         po.oblnkl(self.outfile)
            #         write(self.outfile,31)
            # 31        format(t2,                'Initial centrepost direct cost in included in capital ',                'cost and replacements are assumed to be a fuel cost.')

            po.oblnkl(self.outfile)
            po.ocosts(
                self.outfile,
                "(c222)",
                "Total account 222 cost (M$)",
                cost_variables.c222,
            )

        po.oshead(self.outfile, "Power Injection")

        if ife_variables.ife == 1:
            po.ocosts(
                self.outfile,
                "(c2231)",
                "IFE driver system cost (M$)",
                self.c2231,
            )
        else:
            po.ocosts(self.outfile, "(c2231)", "ECH system cost (M$)", self.c2231)
            po.ocosts(
                self.outfile,
                "(c2232)",
                "Lower hybrid system cost (M$)",
                self.c2232,
            )
            po.ocosts(
                self.outfile,
                "(c2233)",
                "Neutral beam system cost (M$)",
                self.c2233,
            )

        po.oblnkl(self.outfile)
        po.ocosts(self.outfile, "(c223)", "Total account 223 cost (M$)", self.c223)

        po.oshead(self.outfile, "Vacuum Systems")
        po.ocosts(
            self.outfile,
            "(c2241)",
            "High vacuum pumps cost (M$)",
            self.c2241,
        )
        po.ocosts(self.outfile, "(c2242)", "Backing pumps cost (M$)", self.c2242)
        po.ocosts(self.outfile, "(c2243)", "Vacuum duct cost (M$)", self.c2243)
        po.ocosts(self.outfile, "(c2244)", "Valves cost (M$)", self.c2244)
        po.ocosts(self.outfile, "(c2245)", "Duct shielding cost (M$)", self.c2245)
        po.ocosts(self.outfile, "(c2246)", "Instrumentation cost (M$)", self.c2246)
        po.oblnkl(self.outfile)
        po.ocosts(self.outfile, "(c224)", "Total account 224 cost (M$)", self.c224)

        if ife_variables.ife != 1:
            po.oshead(self.outfile, "Power Conditioning")
            po.ocosts(
                self.outfile,
                "(c22511)",
                "TF coil power supplies cost (M$)",
                self.c22511,
            )
            po.ocosts(
                self.outfile,
                "(c22512)",
                "TF coil breakers cost (M$)",
                self.c22512,
            )
            po.ocosts(
                self.outfile,
                "(c22513)",
                "TF coil dump resistors cost (M$)",
                self.c22513,
            )
            po.ocosts(
                self.outfile,
                "(c22514)",
                "TF coil instrumentation and control (M$)",
                self.c22514,
            )
            po.ocosts(
                self.outfile,
                "(c22515)",
                "TF coil bussing cost (M$)",
                self.c22515,
            )
            po.ocosts(
                self.outfile,
                "(c2251)",
                "Total, TF coil power costs (M$)",
                self.c2251,
            )
            po.ocosts(
                self.outfile,
                "(c22521)",
                "PF coil power supplies cost (M$)",
                self.c22521,
            )
            po.ocosts(
                self.outfile,
                "(c22522)",
                "PF coil instrumentation and control (M$)",
                self.c22522,
            )
            po.ocosts(
                self.outfile,
                "(c22523)",
                "PF coil bussing cost (M$)",
                self.c22523,
            )
            po.ocosts(
                self.outfile,
                "(c22524)",
                "PF coil burn power supplies cost (M$)",
                self.c22524,
            )
            po.ocosts(
                self.outfile,
                "(c22525)",
                "PF coil breakers cost (M$)",
                self.c22525,
            )
            po.ocosts(
                self.outfile,
                "(c22526)",
                "PF coil dump resistors cost (M$)",
                self.c22526,
            )
            po.ocosts(
                self.outfile,
                "(c22527)",
                "PF coil ac breakers cost (M$)",
                self.c22527,
            )
            po.ocosts(
                self.outfile,
                "(c2252)",
                "Total, PF coil power costs (M$)",
                self.c2252,
            )
            po.ocosts(
                self.outfile,
                "(c2253)",
                "Total, energy storage cost (M$)",
                self.c2253,
            )
            po.oblnkl(self.outfile)
            po.ocosts(
                self.outfile,
                "(c225)",
                "Total account 225 cost (M$)",
                self.c225,
            )

        po.oshead(self.outfile, "Heat Transport System")
        po.ocosts(
            self.outfile,
            "(cpp)",
            "Pumps and piping system cost (M$)",
            self.cpp,
        )
        po.ocosts(
            self.outfile,
            "(chx)",
            "Primary heat exchanger cost (M$)",
            self.chx,
        )
        po.ocosts(
            self.outfile,
            "(c2261)",
            "Total, reactor cooling system cost (M$)",
            self.c2261,
        )
        po.ocosts(self.outfile, "(cppa)", "Pumps, piping cost (M$)", self.cppa)
        po.ocosts(
            self.outfile,
            "(c2262)",
            "Total, auxiliary cooling system cost (M$)",
            self.c2262,
        )
        po.ocosts(
            self.outfile,
            "(c2263)",
            "Total, cryogenic system cost (M$)",
            self.c2263,
        )
        po.oblnkl(self.outfile)
        po.ocosts(self.outfile, "(c226)", "Total account 226 cost (M$)", self.c226)

        po.oshead(self.outfile, "Fuel Handling System")
        po.ocosts(self.outfile, "(c2271)", "Fuelling system cost (M$)", self.c2271)
        po.ocosts(
            self.outfile,
            "(c2272)",
            "Fuel processing and purification cost (M$)",
            self.c2272,
        )
        po.ocosts(
            self.outfile,
            "(c2273)",
            "Atmospheric recovery systems cost (M$)",
            self.c2273,
        )
        po.ocosts(
            self.outfile,
            "(c2274)",
            "Nuclear building ventilation cost (M$)",
            self.c2274,
        )
        po.oblnkl(self.outfile)
        po.ocosts(self.outfile, "(c227)", "Total account 227 cost (M$)", self.c227)

        po.oshead(self.outfile, "Instrumentation and Control")
        po.ocosts(
            self.outfile,
            "(c228)",
            "Instrumentation and control cost (M$)",
            self.c228,
        )

        po.oshead(self.outfile, "Maintenance Equipment")
        po.ocosts(
            self.outfile,
            "(c229)",
            "Maintenance equipment cost (M$)",
            self.c229,
        )

        po.oshead(self.outfile, "Total Account 22 Cost")
        po.ocosts(self.outfile, "(c22)", "Total account 22 cost (M$)", self.c22)

        po.oshead(self.outfile, "Turbine Plant Equipment")
        po.ocosts(
            self.outfile,
            "(c23)",
            "Turbine plant equipment cost (M$)",
            self.c23,
        )

        po.oshead(self.outfile, "Electric Plant Equipment")
        po.ocosts(
            self.outfile,
            "(c241)",
            "Switchyard equipment cost (M$)",
            self.c241,
        )
        po.ocosts(self.outfile, "(c242)", "Transformers cost (M$)", self.c242)
        po.ocosts(
            self.outfile,
            "(c243)",
            "Low voltage equipment cost (M$)",
            self.c243,
        )
        po.ocosts(
            self.outfile,
            "(c244)",
            "Diesel backup equipment cost (M$)",
            self.c244,
        )
        po.ocosts(
            self.outfile,
            "(c245)",
            "Auxiliary facilities cost (M$)",
            self.c245,
        )
        po.oblnkl(self.outfile)
        po.ocosts(self.outfile, "(c24)", "Total account 24 cost (M$)", self.c24)

        po.oshead(self.outfile, "Miscellaneous Plant Equipment")
        po.ocosts(
            self.outfile,
            "(c25)",
            "Miscellaneous plant equipment cost (M$)",
            self.c25,
        )

        po.oshead(self.outfile, "Heat Rejection System")
        po.ocosts(
            self.outfile,
            "(c26)",
            "Heat rejection system cost (M$)",
            self.c26,
        )

        po.oshead(self.outfile, "Plant Direct Cost")
        po.ocosts(
            self.outfile, "(cdirt)", "Plant direct cost (M$)", cost_variables.cdirt
        )

        po.oshead(self.outfile, "Reactor Core Cost")
        po.ocosts(
            self.outfile,
            "(crctcore)",
            "Reactor core cost (M$)",
            cost_variables.crctcore,
        )

        po.oshead(self.outfile, "Indirect Cost")
        po.ocosts(self.outfile, "(c9)", "Indirect cost (M$)", self.cindrt)

        po.oshead(self.outfile, "Total Contingency")
        po.ocosts(self.outfile, "(ccont)", "Total contingency (M$)", self.ccont)

        po.oshead(self.outfile, "Constructed Cost")
        po.ocosts(
            self.outfile,
            "(concost)",
            "Constructed cost (M$)",
            cost_variables.concost,
        )

        if cost_variables.ireactor == 1:
            po.oshead(self.outfile, "Interest during Construction")
            po.ocosts(
                self.outfile,
                "(moneyint)",
                "Interest during construction (M$)",
                cost_variables.moneyint,
            )

            po.oshead(self.outfile, "Total Capital Investment")
            po.ocosts(
                self.outfile,
                "(capcost)",
                "Total capital investment (M$)",
                cost_variables.capcost,
            )

    def acc22(self):
        """
        Account 22 : Fusion power island
        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 22 (fusion power island
        - the tokamak itself plus auxiliary power systems, etc.) costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
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
        cost_variables.crctcore = cost_variables.c221 + cost_variables.c222 + self.c223

        #  Total account 22
        self.c22 = (
            cost_variables.c221
            + cost_variables.c222
            + self.c223
            + self.c224
            + self.c225
            + self.c226
            + self.c227
            + self.c228
            + self.c229
        )

    def acc221(self):
        """
        Account 221 : Reactor
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 221 (reactor) costs.
        These include the first wall, blanket, shield, support structure
        and divertor plates.
        <P>If ifueltyp = 1, the first wall, blanket and divertor costs are
        treated as fuel costs, rather than as capital costs.
        <P>If ifueltyp = 2, the initial first wall, blanket and divertor costs are
        treated as capital costs, and replacemnts are included as fuel costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
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

        cost_variables.c221 = (
            self.c2211 + self.c2212 + self.c2213 + self.c2214 + self.c2215
        )

    def acc222(self):
        """
        Account 222 : Magnets, including cryostat
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 222 (magnet) costs,
        including the costs of associated cryostats.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        if ife_variables.ife == 1:
            cost_variables.c222 = 0.0e0
            return

        #  Account 222.1 : TF magnet assemblies

        self.acc2221()

        #  Account 222.2 : PF magnet assemblies
        self.acc2222()

        #  Account 222.3 : Cryostat

        self.acc2223()

        #  Total account 222

        cost_variables.c222 = self.c2221 + self.c2222 + self.c2223

    def acc225(self):
        """
        Account 225 : Power conditioning
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 225 (power conditioning) costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        if ife_variables.ife == 1:
            self.c225 = 0.0e0
        else:
            #  Account 225.1 : TF coil power conditioning

            self.acc2251()

            #  Account 225.2 : PF coil power conditioning

            self.acc2252()

            #  Account 225.3 : Energy storage

            self.acc2253()

            #  Total account 225

            self.c225 = self.c2251 + self.c2252 + self.c2253

    def acc21(self):
        """
        Account 21 : Structures and site facilities
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 21 (structures and site
        facilities) costs.
        Building costs are scaled with volume according to algorithms
        developed from TFCX, TFTR, and commercial power plant buildings.
        Costs include equipment, materials and installation labour, but
        no engineering or construction management.
        <P>The general form of the cost algorithm is cost=ucxx*volume**expxx.
        Allowances are used for site improvements and for miscellaneous
        buildings and land costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.6800e0, 0.8400e0, 0.9200e0, 1.0000e0]
        exprb = 1.0e0
        #  Account 211 : Site improvements, facilities and land
        #  N.B. Land unaffected by LSA

        self.c211 = (
            cost_variables.csi * cmlsa[cost_variables.lsa - 1] + cost_variables.cland
        )

        #  Account 212 : Reactor building

        self.c212 = (
            1.0e-6
            * cost_variables.ucrb
            * buildings_variables.rbvol**exprb
            * cmlsa[cost_variables.lsa - 1]
        )

        #  Account 213 : Turbine building

        if cost_variables.ireactor == 1:
            self.c213 = cost_variables.cturbb * cmlsa[cost_variables.lsa - 1]
        else:
            self.c213 = 0.0e0

        #  Account 214 : Reactor maintenance and warm shops buildings

        self.c2141 = (
            1.0e-6
            * cost_variables.ucmb
            * buildings_variables.rmbvol**exprb
            * cmlsa[cost_variables.lsa - 1]
        )
        self.c2142 = (
            1.0e-6
            * cost_variables.ucws
            * buildings_variables.wsvol**exprb
            * cmlsa[cost_variables.lsa - 1]
        )
        self.c214 = self.c2141 + self.c2142

        #  Account 215 : Tritium building

        self.c215 = (
            1.0e-6
            * cost_variables.uctr
            * buildings_variables.triv**exprb
            * cmlsa[cost_variables.lsa - 1]
        )

        #  Account 216 : Electrical equipment building

        self.c216 = (
            1.0e-6
            * cost_variables.ucel
            * buildings_variables.elevol**exprb
            * cmlsa[cost_variables.lsa - 1]
        )

        #  Account 217 : Other buildings
        #  Includes administration, control, shops, cryogenic
        #  plant and an allowance for miscellaneous structures

        self.c2171 = (
            1.0e-6
            * cost_variables.ucad
            * buildings_variables.admvol**exprb
            * cmlsa[cost_variables.lsa - 1]
        )
        self.c2172 = (
            1.0e-6
            * cost_variables.ucco
            * buildings_variables.convol**exprb
            * cmlsa[cost_variables.lsa - 1]
        )
        self.c2173 = (
            1.0e-6
            * cost_variables.ucsh
            * buildings_variables.shovol**exprb
            * cmlsa[cost_variables.lsa - 1]
        )
        self.c2174 = (
            1.0e-6
            * cost_variables.uccr
            * buildings_variables.cryvol**exprb
            * cmlsa[cost_variables.lsa - 1]
        )
        self.c217 = self.c2171 + self.c2172 + self.c2173 + self.c2174

        #  Total for Account 21

        self.c21 = (
            self.c211
            + self.c212
            + self.c213
            + self.c214
            + self.c215
            + self.c216
            + self.c217
        )

    def acc2211(self):
        """
        Account 221.1 : First wall
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 221.1 (first wall) costs.
        The first wall cost is scaled linearly with surface area from TFCX.
        If ifueltyp = 1, the first wall cost is treated as a fuel cost,
        rather than as a capital cost.
        If ifueltyp = 2, inital first wall is included as a capital cost,
        and the replacement first wall cost is treated as a fuel costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.5000e0, 0.7500e0, 0.8750e0, 1.0000e0]

        if ife_variables.ife != 1:
            self.c2211 = (
                1.0e-6
                * cmlsa[cost_variables.lsa - 1]
                * (
                    (cost_variables.ucfwa + cost_variables.ucfws)
                    * build_variables.fwarea
                    + cost_variables.ucfwps
                )
            )
        else:
            self.c2211 = (
                1.0e-6
                * cmlsa[cost_variables.lsa - 1]
                * (
                    cost_variables.ucblss
                    * (
                        ife_variables.fwmatm(1, 1)
                        + ife_variables.fwmatm(2, 1)
                        + ife_variables.fwmatm(3, 1)
                    )
                    + ife_variables.uccarb
                    * (
                        ife_variables.fwmatm(1, 2)
                        + ife_variables.fwmatm(2, 2)
                        + ife_variables.fwmatm(3, 2)
                    )
                    + cost_variables.ucblli2o
                    * (
                        ife_variables.fwmatm(1, 4)
                        + ife_variables.fwmatm(2, 4)
                        + ife_variables.fwmatm(3, 4)
                    )
                    + ife_variables.ucconc
                    * (
                        ife_variables.fwmatm(1, 5)
                        + ife_variables.fwmatm(2, 5)
                        + ife_variables.fwmatm(3, 5)
                    )
                )
            )

        self.c2211 = cost_variables.fkind * self.c2211

        if cost_variables.ifueltyp == 1:
            cost_variables.fwallcst = self.c2211
            self.c2211 = 0.0e0
        elif cost_variables.ifueltyp == 2:
            cost_variables.fwallcst = self.c2211
        else:
            cost_variables.fwallcst = 0.0e0

    def acc2212(self):
        """
        Account 221.2 : Blanket
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 221.2 (blanket) costs.
        If ifueltyp = 1, the blanket cost is treated as a fuel cost,
        rather than as a capital cost.
        If ifueltyp = 2, the initial blanket is included as a capital cost
        and the replacement blanket costs are treated as a fuel cost.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.5000e0, 0.7500e0, 0.8750e0, 1.0000e0]

        if ife_variables.ife != 1:
            # iblanket=4 is used for KIT HCLL model. iblanket<4 are all
            # HCPB (CCFE, KIT and CCFE + Shimwell TBR calculation).

            if fwbs_variables.iblanket == 4:
                #  Liquid blanket (LiPb + Li)
                self.c22121 = 1.0e-6 * fwbs_variables.wtbllipb * cost_variables.ucbllipb
                self.c22122 = 1.0e-6 * fwbs_variables.whtblli * cost_variables.ucblli
            else:
                #  Solid blanket (Li2O + Be)
                self.c22121 = 1.0e-6 * fwbs_variables.whtblbe * cost_variables.ucblbe
                if fwbs_variables.iblanket == 2:
                    # KIT model
                    self.c22122 = (
                        1.0e-6 * fwbs_variables.whtblbreed * cost_variables.ucblbreed
                    )
                else:
                    # CCFE model
                    self.c22122 = (
                        1.0e-6 * fwbs_variables.wtblli2o * cost_variables.ucblli2o
                    )

            self.c22123 = 1.0e-6 * fwbs_variables.whtblss * cost_variables.ucblss
            self.c22124 = 1.0e-6 * fwbs_variables.whtblvd * cost_variables.ucblvd
            self.c22125 = 0.0e0
            self.c22126 = 0.0e0
            self.c22127 = 0.0e0

        else:
            #  IFE blanket; materials present are Li2O, steel, carbon, concrete,
            #  FLiBe and lithium

            self.c22121 = 0.0e0
            self.c22122 = 1.0e-6 * fwbs_variables.wtblli2o * cost_variables.ucblli2o
            self.c22123 = 1.0e-6 * fwbs_variables.whtblss * cost_variables.ucblss
            self.c22124 = 0.0e0
            self.c22125 = (
                1.0e-6
                * ife_variables.uccarb
                * (
                    ife_variables.blmatm(1, 2)
                    + ife_variables.blmatm(2, 2)
                    + ife_variables.blmatm(3, 2)
                )
            )
            self.c22126 = (
                1.0e-6
                * ife_variables.ucconc
                * (
                    ife_variables.blmatm(1, 5)
                    + ife_variables.blmatm(2, 5)
                    + ife_variables.blmatm(3, 5)
                )
            )
            self.c22127 = 1.0e-6 * ife_variables.ucflib * ife_variables.mflibe
            self.c22128 = 1.0e-6 * cost_variables.ucblli * fwbs_variables.whtblli

        self.c22121 = cost_variables.fkind * self.c22121 * cmlsa[cost_variables.lsa - 1]
        self.c22122 = cost_variables.fkind * self.c22122 * cmlsa[cost_variables.lsa - 1]
        self.c22123 = cost_variables.fkind * self.c22123 * cmlsa[cost_variables.lsa - 1]
        self.c22124 = cost_variables.fkind * self.c22124 * cmlsa[cost_variables.lsa - 1]
        self.c22125 = cost_variables.fkind * self.c22125 * cmlsa[cost_variables.lsa - 1]
        self.c22126 = cost_variables.fkind * self.c22126 * cmlsa[cost_variables.lsa - 1]
        self.c22127 = cost_variables.fkind * self.c22127 * cmlsa[cost_variables.lsa - 1]

        self.c2212 = (
            self.c22121
            + self.c22122
            + self.c22123
            + self.c22124
            + self.c22125
            + self.c22126
            + self.c22127
        )

        if cost_variables.ifueltyp == 1:
            cost_variables.blkcst = self.c2212
            self.c2212 = 0.0e0
        elif cost_variables.ifueltyp == 2:
            cost_variables.blkcst = self.c2212
        else:
            cost_variables.blkcst = 0.0e0

    def acc2213(self):
        """
        Account 221.3 : Shield
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 221.3 (shield) costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.5000e0, 0.7500e0, 0.8750e0, 1.0000e0]

        if ife_variables.ife != 1:
            self.c22131 = (
                1.0e-6
                * fwbs_variables.whtshld
                * cost_variables.ucshld
                * cmlsa[cost_variables.lsa - 1]
            )
        else:
            self.c22131 = (
                1.0e-6
                * cmlsa[cost_variables.lsa - 1]
                * (
                    cost_variables.ucshld
                    * (
                        ife_variables.shmatm(1, 1)
                        + ife_variables.shmatm(2, 1)
                        + ife_variables.shmatm(3, 1)
                    )
                    + ife_variables.uccarb
                    * (
                        ife_variables.shmatm(1, 2)
                        + ife_variables.shmatm(2, 2)
                        + ife_variables.shmatm(3, 2)
                    )
                    + cost_variables.ucblli2o
                    * (
                        ife_variables.shmatm(1, 4)
                        + ife_variables.shmatm(2, 4)
                        + ife_variables.shmatm(3, 4)
                    )
                    + ife_variables.ucconc
                    * (
                        ife_variables.shmatm(1, 5)
                        + ife_variables.shmatm(2, 5)
                        + ife_variables.shmatm(3, 5)
                    )
                )
            )

        self.c22131 = cost_variables.fkind * self.c22131

        #  Penetration shield assumed to be typical steel plate
        if ife_variables.ife != 1:
            self.c22132 = (
                1.0e-6
                * fwbs_variables.wpenshld
                * cost_variables.ucpens
                * cmlsa[cost_variables.lsa - 1]
            )
        else:
            self.c22132 = 0.0e0

        self.c22132 = cost_variables.fkind * self.c22132

        self.c2213 = self.c22131 + self.c22132

    def acc2214(self):
        """
        Account 221.4 : Reactor structure
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 221.4 (reactor structure) costs.
        The structural items are costed as standard steel elements.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.6700e0, 0.8350e0, 0.9175e0, 1.0000e0]

        self.c2214 = (
            1.0e-6
            * structure_variables.gsmass
            * cost_variables.ucgss
            * cmlsa[cost_variables.lsa - 1]
        )
        self.c2214 = cost_variables.fkind * self.c2214

    def acc2215(self):
        """
        Account 221.5 : Divertor
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 221.5 (divertor) costs.
        The cost of the divertor blade is scaled linearly with
        surface area from TFCX. The graphite armour is assumed to
        be brazed to water-cooled machined copper substrate.
        Tenth-of-a-kind engineering and installation is assumed.
        <P>If ifueltyp = 1, the divertor cost is treated as a fuel cost,
        rather than as a capital cost.
        <P>If ifueltyp = 2, the initial divertor is included as a capital cost
        and the replacement divertor costs ae treated as a fuel cost,
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        if ife_variables.ife != 1:
            self.c2215 = 1.0e-6 * divertor_variables.divsur * cost_variables.ucdiv
            self.c2215 = cost_variables.fkind * self.c2215

            if cost_variables.ifueltyp == 1:
                cost_variables.divcst = self.c2215
                self.c2215 = 0.0e0
            elif cost_variables.ifueltyp == 2:
                cost_variables.divcst = self.c2215
            else:
                cost_variables.divcst = 0.0e0

        else:
            self.c2215 = 0.0e0
            cost_variables.divcst = 0.0e0

    def acc2221(self):
        """
        Account 222.1 : TF magnet assemblies
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 222.1 (TF magnet) costs.
        Copper magnets are costed from the TFCX data base ($/kg).
        Superconductor magnets are costed using a new method devised
        by R. Hancox under contract to Culham Laboratory, Jan/Feb 1994.
        If ifueltyp = 1, the TART centrepost cost is treated as a fuel
        cost, rather than as a capital cost.
        If ifueltyp = 2, the  initial centrepost is included as a capital cost
        and the replacement TART centrepost costs are treated as a fuel
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.6900e0, 0.8450e0, 0.9225e0, 1.0000e0]

        if tfcoil_variables.i_tf_sup != 1:  # Resistive TF coils
            #  Account 222.1.1 : Inboard TF coil legs

            self.c22211 = (
                1.0e-6
                * tfcoil_variables.whtcp
                * cost_variables.uccpcl1
                * cmlsa[cost_variables.lsa - 1]
            )
            self.c22211 = cost_variables.fkind * self.c22211

            cost_variables.cpstcst = 0.0e0  # TART centrepost
            if (physics_variables.itart == 1) and (cost_variables.ifueltyp == 1):
                cost_variables.cpstcst = self.c22211
                self.c22211 = 0.0e0
            elif (physics_variables.itart == 1) and (cost_variables.ifueltyp == 2):
                cost_variables.cpstcst = self.c22211

            #  Account 222.1.2 : Outboard TF coil legs

            self.c22212 = (
                1.0e-6
                * tfcoil_variables.whttflgs
                * cost_variables.uccpclb
                * cmlsa[cost_variables.lsa - 1]
            )
            self.c22212 = cost_variables.fkind * self.c22212

            #  Total (copper) TF coil costs

            self.c2221 = self.c22211 + self.c22212

        else:  # Superconducting TF coils
            #  Account 222.1.1 : Conductor

            #  Superconductor ($/m)

            costtfsc = (
                cost_variables.ucsc[tfcoil_variables.i_tf_sc_mat - 1]
                * tfcoil_variables.whtconsc
                / (tfcoil_variables.tfleng * tfcoil_variables.n_tf_turn)
            )

            #  Copper ($/m)

            costtfcu = (
                cost_variables.uccu
                * tfcoil_variables.whtconcu
                / (tfcoil_variables.tfleng * tfcoil_variables.n_tf_turn)
            )

            #  Total cost/metre of superconductor and copper wire

            costwire = costtfsc + costtfcu

            #  Total cost/metre of conductor (including sheath and fixed costs)

            ctfconpm = costwire + cost_variables.cconshtf + cost_variables.cconfix

            #  Total conductor costs

            self.c22211 = (
                1.0e-6
                * ctfconpm
                * tfcoil_variables.n_tf
                * tfcoil_variables.tfleng
                * tfcoil_variables.n_tf_turn
            )
            self.c22211 = (
                cost_variables.fkind * self.c22211 * cmlsa[cost_variables.lsa - 1]
            )

            #  Account 222.1.2 : Winding

            self.c22212 = (
                1.0e-6
                * cost_variables.ucwindtf
                * tfcoil_variables.n_tf
                * tfcoil_variables.tfleng
                * tfcoil_variables.n_tf_turn
            )
            self.c22212 = (
                cost_variables.fkind * self.c22212 * cmlsa[cost_variables.lsa - 1]
            )

            #  Account 222.1.3 : Case

            self.c22213 = (
                1.0e-6
                * (tfcoil_variables.whtcas * cost_variables.uccase)
                * tfcoil_variables.n_tf
            )
            self.c22213 = (
                cost_variables.fkind * self.c22213 * cmlsa[cost_variables.lsa - 1]
            )

            #  Account 222.1.4 : Intercoil structure

            self.c22214 = 1.0e-6 * structure_variables.aintmass * cost_variables.ucint
            self.c22214 = (
                cost_variables.fkind * self.c22214 * cmlsa[cost_variables.lsa - 1]
            )

            #  Account 222.1.5 : Gravity support structure

            self.c22215 = 1.0e-6 * structure_variables.clgsmass * cost_variables.ucgss
            self.c22215 = (
                cost_variables.fkind * self.c22215 * cmlsa[cost_variables.lsa - 1]
            )

            #  Total (superconducting) TF coil costs

            self.c2221 = (
                self.c22211 + self.c22212 + self.c22213 + self.c22214 + self.c22215
            )

    def acc2222(self):
        """
        Account 222.2 : PF magnet assemblies
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 222.2 (PF magnet) costs.
        Conductor costs previously used an algorithm devised by R. Hancox,
        January 1994, under contract to Culham, which took into
        account the fact that the superconductor/copper ratio in
        the conductor is proportional to the maximum field that
        each coil will experience. Now, the input copper fractions
        are used instead.
        Maximum values for current, current density and field
        are used.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.6900e0, 0.8450e0, 0.9225e0, 1.0000e0]

        #  Total length of PF coil windings (m)

        pfwndl = 0.0e0
        for i in range(0, pfcoil_variables.nohc):
            pfwndl = (
                pfwndl
                + constants.twopi * pfcoil_variables.rpf[i] * pfcoil_variables.turns[i]
            )

        #  Account 222.2.1 : Conductor

        #  The following lines take care of resistive coils.
        #  costpfsh is the cost per metre of the steel conduit/sheath around
        #  each superconducting cable (so is zero for resistive coils)

        if pfcoil_variables.ipfres == 1:
            costpfsh = 0.0e0
        else:
            costpfsh = cost_variables.cconshpf

        #  Non-Central Solenoid coils

        if build_variables.iohcl == 1:
            npf = pfcoil_variables.nohc - 1
        else:
            npf = pfcoil_variables.nohc

        self.c22221 = 0.0e0

        for i in range(0, npf):
            #  Superconductor ($/m)
            if pfcoil_variables.ipfres == 0:
                costpfsc = (
                    cost_variables.ucsc[pfcoil_variables.isumatpf - 1]
                    * (1.0e0 - pfcoil_variables.fcupfsu)
                    * (1.0e0 - pfcoil_variables.vf[i])
                    * abs(pfcoil_variables.ric[i] / pfcoil_variables.turns[i])
                    * 1.0e6
                    / pfcoil_variables.rjconpf[i]
                    * tfcoil_variables.dcond[pfcoil_variables.isumatpf - 1]
                )
            else:
                costpfsc = 0.0e0

            #  Copper ($/m)
            if pfcoil_variables.ipfres == 0:
                costpfcu = (
                    cost_variables.uccu
                    * pfcoil_variables.fcupfsu
                    * (1.0e0 - pfcoil_variables.vf[i])
                    * abs(pfcoil_variables.ric[i] / pfcoil_variables.turns[i])
                    * 1.0e6
                    / pfcoil_variables.rjconpf[i]
                    * constants.dcopper
                )
            else:
                costpfcu = (
                    cost_variables.uccu
                    * (1.0e0 - pfcoil_variables.vf[i])
                    * abs(pfcoil_variables.ric[i] / pfcoil_variables.turns[i])
                    * 1.0e6
                    / pfcoil_variables.rjconpf[i]
                    * constants.dcopper
                )

            #  Total cost/metre of superconductor and copper wire

            costwire = costpfsc + costpfcu

            #  Total cost/metre of conductor (including sheath and fixed costs)

            cpfconpm = costwire + costpfsh + cost_variables.cconfix

            #  Total account 222.2.1 (PF coils excluding Central Solenoid)

            self.c22221 = self.c22221 + (
                1.0e-6
                * constants.twopi
                * pfcoil_variables.rpf[i]
                * pfcoil_variables.turns[i]
                * cpfconpm
            )

        #  Central Solenoid

        if build_variables.iohcl == 1:
            #  Superconductor ($/m)
            #  Issue #328  Use CS conductor cross-sectional area (m2)
            if pfcoil_variables.ipfres == 0:
                costpfsc = (
                    cost_variables.ucsc[pfcoil_variables.isumatoh - 1]
                    * pfcoil_variables.awpoh
                    * (1 - pfcoil_variables.vfohc)
                    * (1 - pfcoil_variables.fcuohsu)
                    / pfcoil_variables.turns[pfcoil_variables.nohc - 1]
                    * tfcoil_variables.dcond[pfcoil_variables.isumatoh - 1]
                )
            else:
                costpfsc = 0.0e0

            #  Copper ($/m)

            if pfcoil_variables.ipfres == 0:
                costpfcu = (
                    cost_variables.uccu
                    * pfcoil_variables.awpoh
                    * (1 - pfcoil_variables.vfohc)
                    * pfcoil_variables.fcuohsu
                    / pfcoil_variables.turns[pfcoil_variables.nohc - 1]
                    * constants.dcopper
                )
            else:
                # MDK I don't know if this is ccorrect as we never use the resistive model
                costpfcu = (
                    cost_variables.uccu
                    * pfcoil_variables.awpoh
                    * (1 - pfcoil_variables.vfohc)
                    / pfcoil_variables.turns[pfcoil_variables.nohc - 1]
                    * constants.dcopper
                )

            #  Total cost/metre of superconductor and copper wire (Central Solenoid)

            costwire = costpfsc + costpfcu

            #  Total cost/metre of conductor (including sheath and fixed costs)

            cpfconpm = costwire + costpfsh + cost_variables.cconfix

            #  Total account 222.2.1 (PF+Central Solenoid coils)

            self.c22221 = self.c22221 + (
                1.0e-6
                * constants.twopi
                * pfcoil_variables.rpf[pfcoil_variables.nohc - 1]
                * pfcoil_variables.turns[pfcoil_variables.nohc - 1]
                * cpfconpm
            )

        self.c22221 = cost_variables.fkind * self.c22221 * cmlsa[cost_variables.lsa - 1]

        #  Account 222.2.2 : Winding

        self.c22222 = 1.0e-6 * cost_variables.ucwindpf * pfwndl
        self.c22222 = cost_variables.fkind * self.c22222 * cmlsa[cost_variables.lsa - 1]

        #  Account 222.2.3 : Steel case - will be zero for resistive coils

        self.c22223 = 1.0e-6 * cost_variables.uccase * pfcoil_variables.whtpfs
        self.c22223 = cost_variables.fkind * self.c22223 * cmlsa[cost_variables.lsa - 1]

        #  Account 222.2.4 : Support structure

        self.c22224 = 1.0e-6 * cost_variables.ucfnc * structure_variables.fncmass
        self.c22224 = cost_variables.fkind * self.c22224 * cmlsa[cost_variables.lsa - 1]

        #  Total account 222.2

        self.c2222 = self.c22221 + self.c22222 + self.c22223 + self.c22224

    def acc2223(self):
        """
        Account 222.3 : Vacuum vessel
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 222.3 (vacuum vessel) costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.6900e0, 0.8450e0, 0.9225e0, 1.0000e0]

        self.c2223 = 1.0e-6 * fwbs_variables.vvmass * cost_variables.uccryo
        self.c2223 = cost_variables.fkind * self.c2223 * cmlsa[cost_variables.lsa - 1]

    def acc223(self):
        """
        Account 223 : Power injection
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 223 (power injection) costs.
        The costs are from TETRA, updated to 1990$.
        Nominal TIBER values are used pending system designs. Costs are
        scaled linearly with power injected into the plasma and include
        the power supplies.
        <P>If ifueltyp=1, the fraction (1-fcdfuel) of the cost of the
        current drive system is considered as capital cost, and the
        fraction (fcdfuel) is considered a recurring fuel cost due
        to the system's short life.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """

        exprf = 1.0e0
        if ife_variables.ife != 1:
            #  Account 223.1 : ECH

            self.c2231 = (
                1.0e-6
                * cost_variables.ucech
                * (1.0e6 * current_drive_variables.echpwr) ** exprf
            )

            if cost_variables.ifueltyp == 1:
                self.c2231 = (1.0e0 - cost_variables.fcdfuel) * self.c2231
                self.c2231 = cost_variables.fkind * self.c2231

            #  Account 223.2 : Lower Hybrid or ICH

            if current_drive_variables.iefrf != 2:
                self.c2232 = (
                    1.0e-6
                    * cost_variables.uclh
                    * (1.0e6 * current_drive_variables.plhybd) ** exprf
                )
            else:
                self.c2232 = (
                    1.0e-6
                    * cost_variables.ucich
                    * (1.0e6 * current_drive_variables.plhybd) ** exprf
                )

            if cost_variables.ifueltyp == 1:
                self.c2232 = (1.0e0 - cost_variables.fcdfuel) * self.c2232
                self.c2232 = cost_variables.fkind * self.c2232

                #  Account 223.3 : Neutral Beam

                # self.c2233 = 1.0e-6 * cost_variables.ucnbi * (1.0e6*pnbeam)**exprf
                # #327

                self.c2233 = (
                    1.0e-6
                    * cost_variables.ucnbi
                    * (1.0e6 * current_drive_variables.pnbitot) ** exprf
                )

            if cost_variables.ifueltyp == 1:
                self.c2233 = (1.0e0 - cost_variables.fcdfuel) * self.c2233
                self.c2233 = cost_variables.fkind * self.c2233

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
                    self.c2231 = ife_variables.mcdriv * (
                        ife_variables.cdriv1
                        + ife_variables.dcdrv1 * 1.0e-6 * ife_variables.edrive
                    )
                else:
                    self.c2231 = ife_variables.mcdriv * (
                        ife_variables.cdriv2
                        + ife_variables.dcdrv2 * 1.0e-6 * ife_variables.edrive
                    )

            elif ife_variables.ifedrv == 3:
                self.c2231 = (
                    ife_variables.mcdriv
                    * 1.0e-6
                    * ife_variables.cdriv3
                    * (ife_variables.edrive / ife_variables.etadrv)
                )
            else:
                self.c2231 = ife_variables.mcdriv * (
                    ife_variables.cdriv0
                    + ife_variables.dcdrv0 * 1.0e-6 * ife_variables.edrive
                )

            if cost_variables.ifueltyp == 1:
                self.c2231 = (1.0e0 - cost_variables.fcdfuel) * self.c2231
                self.c2231 = cost_variables.fkind * self.c2231
                self.c2232 = 0.0e0
                self.c2233 = 0.0e0
                self.c2234 = 0.0e0

        #  Total account 223

        self.c223 = self.c2231 + self.c2232 + self.c2233 + self.c2234
        cost_variables.cdcost = self.c223

    def acc224(self):
        """
        Account 224 : Vacuum system
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 224 (vacuum system) costs.
        The costs are scaled from TETRA reactor code runs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        if vacuum_variables.ntype == 1:
            self.c2241 = 1.0e-6 * vacuum_variables.vpumpn * cost_variables.uccpmp
        else:
            self.c2241 = 1.0e-6 * vacuum_variables.vpumpn * cost_variables.uctpmp

        self.c2241 = cost_variables.fkind * self.c2241

        #  Account 224.2 : Backing pumps

        self.c2242 = 1.0e-6 * vacuum_variables.nvduct * cost_variables.ucbpmp
        self.c2242 = cost_variables.fkind * self.c2242

        #  Account 224.3 : Vacuum duct

        self.c2243 = (
            1.0e-6
            * vacuum_variables.nvduct
            * vacuum_variables.dlscal
            * cost_variables.ucduct
        )
        self.c2243 = cost_variables.fkind * self.c2243

        #  Account 224.4 : Valves

        self.c2244 = (
            1.0e-6
            * 2.0e0
            * vacuum_variables.nvduct
            * (vacuum_variables.vcdimax * 1.2e0) ** 1.4e0
            * cost_variables.ucvalv
        )
        self.c2244 = cost_variables.fkind * self.c2244

        #  Account 224.5 : Duct shielding

        self.c2245 = (
            1.0e-6
            * vacuum_variables.nvduct
            * vacuum_variables.vacdshm
            * cost_variables.ucvdsh
        )
        self.c2245 = cost_variables.fkind * self.c2245

        #  Account 224.6 : Instrumentation

        self.c2246 = 1.0e-6 * cost_variables.ucviac
        self.c2246 = cost_variables.fkind * self.c2246

        #  Total account 224

        self.c224 = (
            self.c2241 + self.c2242 + self.c2243 + self.c2244 + self.c2245 + self.c2246
        )

    def acc2251(self):
        """
        Account 225.1 : TF coil power conditioning
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 225.1 (TF coil power
        conditioning) costs.
        Costs are developed based on the major equipment specification
        of the tfcpwr module.  A multiplier is used to account for bulk
        materials and installation.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """

        expel = 0.7e0
        self.c22511 = (
            1.0e-6
            * cost_variables.uctfps
            * (tfcoil_variables.tfckw * 1.0e3 + tfcoil_variables.tfcmw * 1.0e6) ** expel
        )
        self.c22511 = cost_variables.fkind * self.c22511

        #  Account 225.1.2 : TF coil breakers (zero cost for copper coils)

        if tfcoil_variables.i_tf_sup == 1:
            self.c22512 = 1.0e-6 * (
                cost_variables.uctfbr
                * tfcoil_variables.n_tf
                * (tfcoil_variables.cpttf * tfcoil_variables.vtfskv * 1.0e3) ** expel
                + cost_variables.uctfsw * tfcoil_variables.cpttf
            )
        else:
            self.c22512 = 0.0e0

        self.c22512 = cost_variables.fkind * self.c22512

        #  Account 225.1.3 : TF coil dump resistors

        self.c22513 = 1.0e-6 * (
            1.0e9 * cost_variables.uctfdr * tfcoil_variables.estotftgj
            + cost_variables.uctfgr * 0.5e0 * tfcoil_variables.n_tf
        )
        self.c22513 = cost_variables.fkind * self.c22513

        #  Account 225.1.4 : TF coil instrumentation and control

        self.c22514 = 1.0e-6 * cost_variables.uctfic * (30.0e0 * tfcoil_variables.n_tf)
        self.c22514 = cost_variables.fkind * self.c22514

        #  Account 225.1.5 : TF coil bussing

        if tfcoil_variables.i_tf_sup != 1:
            self.c22515 = 1.0e-6 * cost_variables.uctfbus * tfcoil_variables.tfbusmas
        else:
            self.c22515 = (
                1.0e-6
                * cost_variables.ucbus
                * tfcoil_variables.cpttf
                * tfcoil_variables.tfbusl
            )

        self.c22515 = cost_variables.fkind * self.c22515

        #  Total account 225.1

        self.c2251 = self.c22511 + self.c22512 + self.c22513 + self.c22514 + self.c22515

    def acc2252(self):
        """
        Account 225.2 : PF coil power conditioning
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 225.2 (PF coil power
        conditioning) costs.
        Costs are taken from the equipment specification of the
        <A HREF="pfpwr.html">pfpwr</A> routine from the plant power module.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.c22521 = 1.0e-6 * cost_variables.ucpfps * heat_transport_variables.peakmva
        self.c22521 = cost_variables.fkind * self.c22521

        #  Account 225.2.2 : PF coil instrumentation and control

        self.c22522 = (
            1.0e-6 * cost_variables.ucpfic * pf_power_variables.pfckts * 30.0e0
        )
        self.c22522 = cost_variables.fkind * self.c22522

        #  Account 225.2.3 : PF coil bussing

        self.c22523 = (
            1.0e-6
            * cost_variables.ucpfb
            * pf_power_variables.spfbusl
            * pf_power_variables.acptmax
        )
        self.c22523 = cost_variables.fkind * self.c22523

        #  Account 225.2.4 : PF coil burn power supplies

        if pf_power_variables.pfckts != 0.0e0:
            self.c22524 = (
                1.0e-6
                * cost_variables.ucpfbs
                * pf_power_variables.pfckts
                * (pf_power_variables.srcktpm / pf_power_variables.pfckts) ** 0.7e0
            )
        else:
            self.c22524 = 0.0e0

        self.c22524 = cost_variables.fkind * self.c22524

        #  Account 225.2.5 : PF coil breakers

        self.c22525 = (
            1.0e-6
            * cost_variables.ucpfbk
            * pf_power_variables.pfckts
            * (pf_power_variables.acptmax * pf_power_variables.vpfskv) ** 0.7e0
        )
        self.c22525 = cost_variables.fkind * self.c22525

        #  Account 225.2.6 : PF coil dump resistors

        self.c22526 = 1.0e-6 * cost_variables.ucpfdr1 * pf_power_variables.ensxpfm
        self.c22526 = cost_variables.fkind * self.c22526

        #  Account 225.2.7 : PF coil AC breakers

        self.c22527 = 1.0e-6 * cost_variables.ucpfcb * pf_power_variables.pfckts
        self.c22527 = cost_variables.fkind * self.c22527

        #  Total account 225.2

        self.c2252 = (
            self.c22521
            + self.c22522
            + self.c22523
            + self.c22524
            + self.c22525
            + self.c22526
            + self.c22527
        )

    def acc226(self):
        """
        Account 226 : Heat transport system
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 226 (heat transport system) costs.
        Costs are estimated from major equipment and heat transport
        system loops developed in the heatpwr module of the code.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.c226 = self.c2261 + self.c2262 + self.c2263

    def acc2261(self):
        """
        Account 2261 : Reactor cooling system
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 2261 -
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.4000e0, 0.7000e0, 0.8500e0, 1.0000e0]
        exphts = 0.7e0

        #  Pumps and piping system
        #  N.B. with blktmodel > 0, the blanket is assumed to be helium-cooled,
        #  but the shield etc. is water-cooled (coolwh=2). Therefore, a slight
        #  inconsistency exists here...
        self.cpp = (
            1.0e-6
            * cost_variables.uchts[fwbs_variables.coolwh - 1]
            * (
                (1.0e6 * heat_transport_variables.pfwdiv) ** exphts
                + (1.0e6 * fwbs_variables.pnucblkt) ** exphts
                + (1.0e6 * fwbs_variables.pnucshld) ** exphts
            )
        )

        self.cpp = cost_variables.fkind * self.cpp * cmlsa[cost_variables.lsa - 1]

        #  Primary heat exchangers
        self.chx = (
            1.0e-6
            * cost_variables.ucphx
            * heat_transport_variables.nphx
            * (
                1.0e6
                * heat_transport_variables.pthermmw
                / heat_transport_variables.nphx
            )
            ** exphts
        )
        self.chx = cost_variables.fkind * self.chx * cmlsa[cost_variables.lsa - 1]

        self.c2261 = self.chx + self.cpp

    def acc2262(self):
        """
        Account 2262 : Auxiliary component cooling
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 2262 - Auxiliary component cooling
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = 0.4000e0, 0.7000e0, 0.8500e0, 1.0000e0
        exphts = 0.7e0

        #  Pumps and piping system
        self.cppa = (
            1.0e-6
            * cost_variables.ucahts
            * (
                (1.0e6 * heat_transport_variables.pinjht) ** exphts
                + (1.0e6 * heat_transport_variables.crypmw) ** exphts
                + (1.0e6 * heat_transport_variables.vachtmw) ** exphts
                + (1.0e6 * heat_transport_variables.trithtmw) ** exphts
                + (1.0e6 * heat_transport_variables.fachtmw) ** exphts
            )
        )

        if ife_variables.ife == 1:
            self.cppa = self.cppa + 1.0e-6 * cost_variables.ucahts * (
                (1.0e6 * ife_variables.tdspmw) ** exphts
                + (1.0e6 * ife_variables.tfacmw) ** exphts
            )

        #  Apply Nth kind and safety assurance factors
        self.cppa = cost_variables.fkind * self.cppa * cmlsa[cost_variables.lsa - 1]

        self.c2262 = self.cppa

    def acc2263(self):
        """
        Account 2263 : Cryogenic system
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 2263 - Cryogenic system
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = 0.4000e0, 0.7000e0, 0.8500e0, 1.0000e0
        expcry = 0.67e0

        self.c2263 = (
            1.0e-6
            * cost_variables.uccry
            * 4.5e0
            / tfcoil_variables.tmpcry
            * heat_transport_variables.helpow**expcry
        )

        #  Apply Nth kind and safety factors
        self.c2263 = cost_variables.fkind * self.c2263 * cmlsa[cost_variables.lsa - 1]

    def acc227(self):
        """
        Account 227 : Fuel handling
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 227 (fuel handling) costs.
        Costs are scaled from TETRA reactor code runs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.c227 = self.c2271 + self.c2272 + self.c2273 + self.c2274

    def acc2271(self):
        """
        Account 2271 : Fuelling system
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 2271 - Fuelling system
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.c2271 = 1.0e-6 * cost_variables.ucf1

        #  Apply Nth kind factor
        self.c2271 = cost_variables.fkind * self.c2271

    def acc2272(self):
        """
        Account 2272 : Fuel processing and purification
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 2272 - Fuel processing
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        if ife_variables.ife != 1:
            #  Previous calculation, using qfuel in Amps:
            #  1.3 should have been physics_variables.afuel*umass/echarge*1000*s/day = 2.2
            # wtgpd = burnup * qfuel * 1.3e0

            #  New calculation: 2 nuclei * reactions/sec * kg/nucleus * g/kg * sec/day
            physics_variables.wtgpd = (
                2.0e0
                * physics_variables.rndfuel
                * physics_variables.afuel
                * constants.umass
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
                / (1.602e-19 * 17.6e6 * ife_variables.fburn)
            )
            physics_variables.wtgpd = targtm * ife_variables.reprat * 86400.0e0

        #  Assumes that He3 costs same as tritium to process...
        self.c2272 = (
            1.0e-6
            * cost_variables.ucfpr
            * (0.5e0 + 0.5e0 * (physics_variables.wtgpd / 60.0e0) ** 0.67e0)
        )

        self.c2272 = cost_variables.fkind * self.c2272

    def acc2273(self):
        """
        Account 2273 : Atmospheric recovery systems
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 2273 - Atmospheric recovery systems
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cfrht = 1.0e5

        #  No detritiation needed if purely D-He3 reaction
        if physics_variables.ftrit > 1.0e-3:
            self.c2273 = (
                1.0e-6
                * cost_variables.ucdtc
                * (
                    (cfrht / 1.0e4) ** 0.6e0
                    * (buildings_variables.volrci + buildings_variables.wsvol)
                )
            )
        else:
            self.c2273 = 0.0e0

        self.c2273 = cost_variables.fkind * self.c2273

    def acc2274(self):
        """
        Account 2274 : Nuclear building ventilation
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 2274 - Nuclear building ventilation
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.c2274 = (
            1.0e-6
            * cost_variables.ucnbv
            * (buildings_variables.volrci + buildings_variables.wsvol) ** 0.8e0
        )

        #  Apply Nth kind factor
        self.c2274 = cost_variables.fkind * self.c2274

    def acc228(self):
        """
        Account 228 : Instrumentation and control
        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 228 (instrumentation and
        control) costs.
        Costs are based on TFCX and INTOR.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.c228 = 1.0e-6 * cost_variables.uciac
        self.c228 = cost_variables.fkind * self.c228

    def acc229(self):
        """
        Account 229 : Maintenance equipment
        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 229 (maintenance equipment) costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.c229 = 1.0e-6 * cost_variables.ucme
        self.c229 = cost_variables.fkind * self.c229

    def acc23(self):
        """
        Account 23 : Turbine plant equipment
        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 23 (turbine plant equipment) costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """

        exptpe = 0.83e0
        if cost_variables.ireactor == 1:
            self.c23 = (
                1.0e-6
                * cost_variables.ucturb[fwbs_variables.coolwh - 1]
                * (heat_transport_variables.pgrossmw / 1200.0e0) ** exptpe
            )

    def acc24(self):
        """
        Account 24 : Electric plant equipment
        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 24 (electric plant equipment) costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.c24 = self.c241 + self.c242 + self.c243 + self.c244 + self.c245

    def acc241(self):
        """
        Account 241 : Electric plant equipment - switchyard
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 241 - switchyard
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = 0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0

        #  Account 241 : Switchyard
        self.c241 = 1.0e-6 * cost_variables.ucswyd * cmlsa[cost_variables.lsa - 1]

    def acc242(self):
        """
        Account 242 : Electric plant equipment - Transformers
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 242 - Transformers
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = 0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0
        expepe = 0.9e0

        #  Account 242 : Transformers
        self.c242 = 1.0e-6 * (
            cost_variables.ucpp * (heat_transport_variables.pacpmw * 1.0e3) ** expepe
            + cost_variables.ucap * (heat_transport_variables.fcsht * 1.0e3)
        )

        #  Apply safety assurance factor
        self.c242 = self.c242 * cmlsa[cost_variables.lsa - 1]

    def acc243(self):
        """
        Account 243 : Electric plant equipment - Low voltage
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 243 - Low voltage
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = 0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0

        #  Account 243 : Low voltage
        #  (include 0.8 factor for transformer efficiency)
        self.c243 = (
            1.0e-6
            * cost_variables.uclv
            * heat_transport_variables.tlvpmw
            * 1.0e3
            / 0.8e0
            * cmlsa[cost_variables.lsa - 1]
        )

    def acc244(self):
        """
        Account 244 : Electric plant equipment - Diesel generators
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 244 - Diesel generators
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0]

        #  Account 244 : Diesel generator (8 MW per generator,  assume 4 )
        self.c244 = (
            1.0e-6 * cost_variables.ucdgen * 4.0e0 * cmlsa[cost_variables.lsa - 1]
        )

    def acc245(self):
        """
        Account 245 : Electric plant equipment - Aux facility power
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 245 - Aux facility power
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = 0.5700e0, 0.7850e0, 0.8925e0, 1.0000e0

        #  Account 245 : Auxiliary facility power needs
        self.c245 = 1.0e-6 * cost_variables.ucaf * cmlsa[cost_variables.lsa - 1]

    def acc25(self):
        """
        Account 25 : Miscellaneous plant equipment
        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 25 (miscellaneous plant
        equipment) costs, such as waste treatment.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = 0.7700e0, 0.8850e0, 0.9425e0, 1.0000e0

        self.c25 = 1.0e-6 * cost_variables.ucmisc * cmlsa[cost_variables.lsa - 1]

    def acc26(self):
        """
        Account 26 : Heat rejection system
        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 26 (heat rejection system) costs.
        Costs are scaled with the total plant heat rejection based on
        commercial systems.
        J. Delene, private communication, ORNL, June 1990
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        cmlsa = [0.8000e0, 0.9000e0, 0.9500e0, 1.0000e0]

        # Calculate rejected heat for non-reactor (==0) and reactor (==1)
        if cost_variables.ireactor == 0:
            pwrrej = (
                physics_variables.powfmw
                + heat_transport_variables.pinjwp
                + tfcoil_variables.tfcmw
            )
        else:
            pwrrej = (
                heat_transport_variables.pthermmw - heat_transport_variables.pgrossmw
            )

        # cost_variables.uchrs - reference cost of heat rejection system [$]
        self.c26 = (
            1.0e-6
            * cost_variables.uchrs
            * pwrrej
            / 2300.0e0
            * cmlsa[cost_variables.lsa - 1]
        )

    def acc9(self):
        """
        Account 9 : Indirect cost and contingency allowances
        author: P J Knight, CCFE, Culham Science Centre
        author: J Morris, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 9 (indirect cost and
        contingency allowances) costs.
        The cost modelling is based on the commercial plant model of a
        single contractor performing all plant engineering and construction
        management, using commercially purchased equipment and materials.
        <P>The project contingency is an allowance for incomplete design
        specification and unforeseen events during the plant construction.
        <P>The factors used are estimated from commercial plant experience.
        J. Delene, private communication, ORNL, June 1990
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.cindrt = (
            cost_variables.cfind[cost_variables.lsa - 1]
            * cost_variables.cdirt
            * (1.0e0 + cost_variables.cowner)
        )

        #  Contingency costs

        self.ccont = cost_variables.fcontng * (cost_variables.cdirt + self.cindrt)

    def acc2253(self):
        """
        Account 225.3 : Energy storage
        author: P J Knight, CCFE, Culham Science Centre
        None
        This routine evaluates the Account 225.3 (energy storage) costs.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        self.c2253 = 0.0e0

        #  Thermal storage options for a pulsed reactor
        #  See F/MPE/MOD/CAG/PROCESS/PULSE/0008 and 0014

        if pulse_variables.lpulse == 1:
            if pulse_variables.istore == 1:
                #  Option 1 from ELECTROWATT report
                #  Pulsed Fusion Reactor Study : AEA FUS 205

                #  Increased condensate tank capacity
                self.c2253 = 0.1e0

                #  Additional electrically-driven feedpump (50 per cent duty)
                self.c2253 = self.c2253 + 0.8e0

                #  Increased turbine-generator duty (5 per cent duty)
                self.c2253 = self.c2253 + 4.0e0

                #  Additional auxiliary transformer capacity and ancillaries
                self.c2253 = self.c2253 + 0.5e0

                #  Increased drum capacity
                self.c2253 = self.c2253 + 2.8e0

                #  Externally fired superheater
                self.c2253 = self.c2253 + 29.0e0

            elif pulse_variables.istore == 2:
                #  Option 2 from ELECTROWATT report
                #  Pulsed Fusion Reactor Study : AEA FUS 205

                #  Increased condensate tank capacity
                self.c2253 = 0.1e0

                #  Additional electrically-driven feedpump (50 per cent duty)
                self.c2253 = self.c2253 + 0.8e0

                #  Increased drum capacity
                self.c2253 = self.c2253 + 2.8e0

                #  Increased turbine-generator duty (5 per cent duty)
                self.c2253 = self.c2253 + 4.0e0

                #  Additional fired boiler (1 x 100 per cent duty)
                self.c2253 = self.c2253 + 330.0e0

                #  HP/LP steam bypass system for auxiliary boiler
                #  (30 per cent boiler capacity)
                self.c2253 = self.c2253 + 1.0e0

                #  Dump condenser
                self.c2253 = self.c2253 + 2.0e0

                #  Increased cooling water system capacity
                self.c2253 = self.c2253 + 18.0e0

            elif pulse_variables.istore == 3:
                #  Simplistic approach that assumes that a large stainless steel
                #  block acts as the thermal storage medium. No account is taken
                #  of the cost of the piping within the block, etc.
                #
                #  shcss is the specific heat capacity of stainless steel (J/kg/K)
                #  pulse_variables.dtstor is the maximum allowable temperature change in the
                #  stainless steel block (input)

                shcss = 520.0e0
                self.c2253 = (
                    cost_variables.ucblss
                    * (heat_transport_variables.pthermmw * 1.0e6)
                    * times_variables.tdown
                    / (shcss * pulse_variables.dtstor)
                )

            else:
                error_handling.idiags[0] = pulse_variables.istore

                error_handling.report_error(125)

        if pulse_variables.istore < 3:
            #  Scale self.c2253 with net electric power

            self.c2253 = self.c2253 * heat_transport_variables.pnetelmw / 1200.0e0

            #  It is necessary to convert from 1992 pounds to 1990 dollars
            #  Reasonable guess for the exchange rate + inflation factor
            #  inflation = 5% per annum; exchange rate = 1.5 dollars per pound

            self.c2253 = self.c2253 * 1.36e0

        self.c2253 = cost_variables.fkind * self.c2253

    def coelc(self):
        """
        Routine to calculate the cost of electricity for a fusion power plant
        author: P J Knight, CCFE, Culham Science Centre
        outfile : input integer : output file unit

        This routine performs the calculation of the cost of electricity
        for a fusion power plant.
        <P>Annual costs are in megadollars/year, electricity costs are in
        millidollars/kWh, while other costs are in megadollars.
        All values are based on 1990 dollars.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        """
        if ife_variables.ife == 1:
            kwhpy = (
                1.0e3
                * heat_transport_variables.pnetelmw
                * (24.0e0 * constants.n_day_year)
                * cost_variables.cfactr
            )
        else:
            kwhpy = (
                1.0e3
                * heat_transport_variables.pnetelmw
                * (24.0e0 * constants.n_day_year)
                * cost_variables.cfactr
                * times_variables.tburn
                / times_variables.tcycle
            )

        #  Costs due to reactor plant
        #  ==========================

        #  Interest on construction costs

        cost_variables.moneyint = cost_variables.concost * (
            cost_variables.fcap0 - 1.0e0
        )

        #  Capital costs

        cost_variables.capcost = cost_variables.concost + cost_variables.moneyint

        #  Annual cost of plant capital cost

        anncap = cost_variables.capcost * cost_variables.fcr0

        # SJP Issue #836
        # Check for the condition when kwhpy=0

        if kwhpy < 1.0e-10:
            kwhpy = 1.0e-10

        #  Cost of electricity due to plant capital cost

        cost_variables.coecap = 1.0e9 * anncap / kwhpy

        #  Costs due to first wall and blanket renewal
        #  ===========================================

        #  Operational life

        fwbllife = fwbs_variables.bktlife

        #  Compound interest factor

        feffwbl = (1.0e0 + cost_variables.discount_rate) ** fwbllife

        #  Capital recovery factor

        crffwbl = (feffwbl * cost_variables.discount_rate) / (feffwbl - 1.0e0)

        #  Annual cost of replacements

        annfwbl = (
            (cost_variables.fwallcst + cost_variables.blkcst)
            * (1.0e0 + cost_variables.cfind[cost_variables.lsa - 1])
            * cost_variables.fcap0cp
            * crffwbl
        )

        if cost_variables.ifueltyp == 2:
            annfwbl = annfwbl * (1.0e0 - fwbllife / cost_variables.tlife)

        #  Cost of electricity due to first wall/blanket replacements

        coefwbl = 1.0e9 * annfwbl / kwhpy

        #  Costs due to divertor renewal
        #  =============================

        if ife_variables.ife == 1:
            anndiv = 0.0e0
            coediv = 0.0e0
        else:
            #  Compound interest factor

            fefdiv = (1.0e0 + cost_variables.discount_rate) ** cost_variables.divlife

            #  Capital recovery factor

            crfdiv = (fefdiv * cost_variables.discount_rate) / (fefdiv - 1.0e0)

            #  Annual cost of replacements

            anndiv = (
                cost_variables.divcst
                * (1.0e0 + cost_variables.cfind[cost_variables.lsa - 1])
                * cost_variables.fcap0cp
                * crfdiv
            )

            #  Cost of electricity due to divertor replacements

            if cost_variables.ifueltyp == 2:
                anndiv = anndiv * (
                    1.0e0 - cost_variables.divlife / cost_variables.tlife
                )

            coediv = 1.0e9 * anndiv / kwhpy

        #  Costs due to centrepost renewal
        #  ===============================

        if (physics_variables.itart == 1) and (ife_variables.ife != 1):
            #  Compound interest factor

            fefcp = (1.0e0 + cost_variables.discount_rate) ** cost_variables.cplife

            #  Capital recovery factor

            crfcp = (fefcp * cost_variables.discount_rate) / (fefcp - 1.0e0)

            #  Annual cost of replacements

            anncp = (
                cost_variables.cpstcst
                * (1.0e0 + cost_variables.cfind[cost_variables.lsa - 1])
                * cost_variables.fcap0cp
                * crfcp
            )

            #  Cost of electricity due to centrepost replacements
            if cost_variables.ifueltyp == 2:
                anncp = anncp * (1.0e0 - cost_variables.cplife / cost_variables.tlife)

            coecp = 1.0e9 * anncp / kwhpy

        else:
            anncp = 0.0e0
            coecp = 0.0e0

        #  Costs due to partial current drive system renewal
        #  =================================================

        #  Compound interest factor

        fefcdr = (1.0e0 + cost_variables.discount_rate) ** cost_variables.cdrlife

        #  Capital recovery factor

        crfcdr = (fefcdr * cost_variables.discount_rate) / (fefcdr - 1.0e0)

        #  Annual cost of replacements

        if cost_variables.ifueltyp == 0:
            anncdr = 0.0e0
        else:
            anncdr = (
                cost_variables.cdcost
                * cost_variables.fcdfuel
                / (1.0e0 - cost_variables.fcdfuel)
                * (1.0e0 + cost_variables.cfind[cost_variables.lsa - 1])
                * cost_variables.fcap0cp
                * crfcdr
            )

        #  Cost of electricity due to current drive system replacements

        coecdr = 1.0e9 * anncdr / kwhpy

        #  Costs due to operation and maintenance
        #  ======================================

        #  Annual cost of operation and maintenance

        annoam = cost_variables.ucoam[cost_variables.lsa - 1] * numpy.sqrt(
            heat_transport_variables.pnetelmw / 1200.0e0
        )

        #  Additional cost due to pulsed reactor thermal storage
        #  See F/MPE/MOD/CAG/PROCESS/PULSE/0008
        #
        #      if (lpulse.eq.1) :
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
        #         annoam1 = annoam1 * heat_transport_variables.pnetelmw/1200.0e0
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

        cost_variables.coeoam = 1.0e9 * annoam / kwhpy

        #  Costs due to reactor fuel
        #  =========================

        #  Annual cost of fuel

        if ife_variables.ife != 1:
            #  Sum D-T fuel cost and He3 fuel cost
            annfuel = (
                cost_variables.ucfuel * heat_transport_variables.pnetelmw / 1200.0e0
                + 1.0e-6
                * physics_variables.fhe3
                * physics_variables.wtgpd
                * 1.0e-3
                * cost_variables.uche3
                * constants.n_day_year
                * cost_variables.cfactr
            )
        else:
            annfuel = (
                1.0e-6
                * ife_variables.uctarg
                * ife_variables.reprat
                * 3.1536e7
                * cost_variables.cfactr
            )

        #  Cost of electricity due to reactor fuel

        coefuel = 1.0e9 * annfuel / kwhpy

        #  Costs due to waste disposal
        #  ===========================

        #  Annual cost of waste disposal

        annwst = cost_variables.ucwst[cost_variables.lsa - 1] * numpy.sqrt(
            heat_transport_variables.pnetelmw / 1200.0e0
        )

        #  Cost of electricity due to waste disposal

        coewst = 1.0e9 * annwst / kwhpy

        #  Costs due to decommissioning fund
        #  =================================

        #  Annual contributions to fund for decommissioning
        #  A fraction cost_variables.decomf of the construction cost is set aside for
        #  this purpose at the start of the plant life.
        #  Final factor takes into account inflation over the plant lifetime
        #  (suggested by Tim Hender 07/03/96)
        #  Difference (cost_variables.dintrt) between borrowing and saving interest rates is
        #  included, along with the possibility of completing the fund cost_variables.dtlife
        #  years before the end of the plant's lifetime

        anndecom = (
            cost_variables.decomf
            * cost_variables.concost
            * cost_variables.fcr0
            / (1.0e0 + cost_variables.discount_rate - cost_variables.dintrt)
            ** (cost_variables.tlife - cost_variables.dtlife)
        )

        #  Cost of electricity due to decommissioning fund

        coedecom = 1.0e9 * anndecom / kwhpy

        #  Total costs
        #  ===========

        #  Annual costs due to 'fuel-like' components

        # annfuelt = annfwbl + anndiv + anncdr + anncp + annfuel + annwst

        #  Total cost of electricity due to 'fuel-like' components

        cost_variables.coefuelt = coefwbl + coediv + coecdr + coecp + coefuel + coewst

        #  Total annual costs

        # anntot = anncap + annfuelt + annoam + anndecom

        #  Total cost of electricity

        cost_variables.coe = (
            cost_variables.coecap
            + cost_variables.coefuelt
            + cost_variables.coeoam
            + coedecom
        )
