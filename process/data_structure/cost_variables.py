"""Module containing variables for the costs models"""

from dataclasses import dataclass, field
from enum import IntEnum, unique

from process.core.data_structure.parameter import Parameter, PROCESSModelData


@unique
class CostModels(IntEnum):
    """Enum for cost model selection"""

    PROCESS_1990 = 0
    """Use $ 1990 PROCESS model"""

    KOVARI_2014 = 1
    """Use $ 2014 Kovari model"""

    USER_PROVIDED = 2
    """Use user-provided model"""


@dataclass(slots=True)
class CostData(PROCESSModelData):
    """Dataclass holding cost variables"""

    c228: Parameter[float] = 0.0
    """c228 account cost"""

    c229: Parameter[float] = 0.0
    """c229 account cost"""

    c23: Parameter[float] = 0.0
    """c23 account cost"""

    c25: Parameter[float] = 0.0
    """c25 account cost"""

    c26: Parameter[float] = 0.0
    """c26 account cost"""

    cindrt: Parameter[float] = 0.0
    """cindrt account cost"""

    ccont: Parameter[float] = 0.0
    """ccont account cost"""

    c226: Parameter[float] = 0.0

    c2261: Parameter[float] = 0.0

    c2262: Parameter[float] = 0.0

    c2263: Parameter[float] = 0.0

    c227: Parameter[float] = 0.0

    c2271: Parameter[float] = 0.0

    c2272: Parameter[float] = 0.0

    c2273: Parameter[float] = 0.0

    c2274: Parameter[float] = 0.0

    c24: Parameter[float] = 0.0

    c241: Parameter[float] = 0.0

    c242: Parameter[float] = 0.0

    c243: Parameter[float] = 0.0

    c244: Parameter[float] = 0.0

    c245: Parameter[float] = 0.0

    c21: Parameter[float] = 0.0

    c211: Parameter[float] = 0.0

    c212: Parameter[float] = 0.0

    c213: Parameter[float] = 0.0

    c214: Parameter[float] = 0.0

    c2141: Parameter[float] = 0.0

    c2142: Parameter[float] = 0.0

    c215: Parameter[float] = 0.0

    c216: Parameter[float] = 0.0

    c217: Parameter[float] = 0.0

    c2171: Parameter[float] = 0.0

    c2172: Parameter[float] = 0.0

    c2173: Parameter[float] = 0.0

    c2174: Parameter[float] = 0.0

    c22: Parameter[float] = 0.0

    c2211: Parameter[float] = 0.0

    c2212: Parameter[float] = 0.0

    c22121: Parameter[float] = 0.0

    c22122: Parameter[float] = 0.0

    c22123: Parameter[float] = 0.0

    c22124: Parameter[float] = 0.0

    c22125: Parameter[float] = 0.0

    c22126: Parameter[float] = 0.0

    c22127: Parameter[float] = 0.0

    c2213: Parameter[float] = 0.0

    c22131: Parameter[float] = 0.0

    c22132: Parameter[float] = 0.0

    c2214: Parameter[float] = 0.0

    c2215: Parameter[float] = 0.0

    c2221: Parameter[float] = 0.0

    c22211: Parameter[float] = 0.0

    c22212: Parameter[float] = 0.0

    c22213: Parameter[float] = 0.0

    c22214: Parameter[float] = 0.0

    c22215: Parameter[float] = 0.0

    c2222: Parameter[float] = 0.0

    c22221: Parameter[float] = 0.0

    c22222: Parameter[float] = 0.0

    c22223: Parameter[float] = 0.0

    c22224: Parameter[float] = 0.0

    c2223: Parameter[float] = 0.0

    c223: Parameter[float] = 0.0

    c2231: Parameter[float] = 0.0

    c2232: Parameter[float] = 0.0

    c2233: Parameter[float] = 0.0

    c2234: Parameter[float] = 0.0

    c224: Parameter[float] = 0.0

    c2241: Parameter[float] = 0.0

    c2242: Parameter[float] = 0.0

    c2243: Parameter[float] = 0.0

    c2244: Parameter[float] = 0.0

    c2245: Parameter[float] = 0.0

    c2246: Parameter[float] = 0.0

    c225: Parameter[float] = 0.0

    c2251: Parameter[float] = 0.0

    c22511: Parameter[float] = 0.0

    c22512: Parameter[float] = 0.0

    c22513: Parameter[float] = 0.0

    c22514: Parameter[float] = 0.0

    c22515: Parameter[float] = 0.0

    c2252: Parameter[float] = 0.0

    c22521: Parameter[float] = 0.0

    c22522: Parameter[float] = 0.0

    c22523: Parameter[float] = 0.0

    c22524: Parameter[float] = 0.0

    c22525: Parameter[float] = 0.0

    c22526: Parameter[float] = 0.0

    c22527: Parameter[float] = 0.0

    c2253: Parameter[float] = 0.0

    chx: Parameter[float] = 0.0

    cpp: Parameter[float] = 0.0

    cppa: Parameter[float] = 0.0

    c22128: Parameter[float] = 0.0

    abktflnc: Parameter[float] = 5.0
    """allowable first wall/blanket neutron fluence (MW-yr/m2) (`blktmodel=0`)"""

    adivflnc: Parameter[float] = 7.0
    """allowable divertor heat fluence (MW-yr/m2)"""

    blkcst: Parameter[float] = 0.0
    """blanket direct cost (M$)"""

    c221: Parameter[float] = 0.0
    """total account 221 cost (M$) - first wall, blanket, shield, support structure and div plates"""

    c222: Parameter[float] = 0.0
    """total account 222 cost (M$) - TF coils + PF coils"""

    capcost: Parameter[float] = 0.0
    """total capital cost including interest (M$)"""

    cconfix: Parameter[float] = 80.0
    """fixed cost of superconducting cable ($/m)"""

    cconshpf: Parameter[float] = 70.0
    """cost of PF coil steel conduit/sheath ($/m)"""

    cconshtf: Parameter[float] = 75.0
    """cost of TF coil steel conduit/sheath ($/m)"""

    cdcost: Parameter[float] = 0.0
    """current drive direct costs (M$)"""

    cdirt: Parameter[float] = 0.0
    """total plant direct cost (M$)"""

    life_hcd_fpy: Parameter[float] = 0.0
    """Full power year lifetime of heating/current drive system (y)"""

    cdrlife_cal: Parameter[float] = 0.0
    """Calendar year lifetime of heating/current drive system (y)"""

    f_t_plant_available: Parameter[float] = 0.75
    """Total plant availability fraction; input if `i_plant_availability=0`"""

    cpfact: Parameter[float] = 0.0
    """Total plant capacity factor"""

    cfind: Parameter[list[float]] = field(
        default_factory=lambda: [0.244, 0.244, 0.244, 0.29]
    )
    """indirect cost factor (func of lsa) (cost model = 0)"""

    cland: Parameter[float] = 19.2
    """cost of land (M$)"""

    coe: Parameter[float] = 0.0
    """cost of electricity ($/MW-hr)"""

    coecap: Parameter[float] = 0.0
    """capital cost of electricity (m$/kW-hr)"""

    coefuelt: Parameter[float] = 0.0
    """'fuel' (including replaceable components) contribution to cost of electricity (m$/kW-hr)"""

    coeoam: Parameter[float] = 0.0
    """operation and maintenance contribution to cost of electricity (m$/kW-hr)"""

    concost: Parameter[float] = 0.0
    """plant construction cost (M$)"""

    costexp: Parameter[float] = 0.8
    """cost exponent for scaling in 2015 costs model"""

    costexp_pebbles: Parameter[float] = 0.6
    """cost exponent for pebbles in 2015 costs model"""

    cost_factor_buildings: Parameter[float] = 1.0
    """cost scaling factor for buildings"""

    cost_factor_land: Parameter[float] = 1.0
    """cost scaling factor for land"""

    cost_factor_tf_coils: Parameter[float] = 1.0
    """cost scaling factor for TF coils"""

    cost_factor_fwbs: Parameter[float] = 1.0
    """cost scaling factor for fwbs"""

    cost_factor_rh: Parameter[float] = 1.0
    """cost scaling factor for remote handling"""

    cost_factor_vv: Parameter[float] = 1.0
    """cost scaling factor for vacuum vessel"""

    cost_factor_bop: Parameter[float] = 1.0
    """cost scaling factor for energy conversion system"""

    cost_factor_misc: Parameter[float] = 1.0
    """cost scaling factor for remaining subsystems"""

    maintenance_fwbs: Parameter[float] = 0.2
    """Maintenance cost factor: first wall, blanket, shield, divertor"""

    maintenance_gen: Parameter[float] = 0.05
    """Maintenance cost factor: All other components except coils, vacuum vessel,
    thermal shield, cryostat, land"""

    amortization: Parameter[float] = 13.6
    """amortization factor (fixed charge factor) "A" (years)"""

    i_cost_model: int = 1
    """Switch for cost model:
    - =0 use $ 1990 PROCESS model
    - =1 use $ 2014 Kovari model
    - =2 use user-provided model
    """

    i_cp_lifetime: int = 0
    """Switch for the centrepost lifetime constraint
    0 : The CP full power year lifetime is set by the user via cplife_input
    1 : The CP lifetime is equal to the divertor lifetime
    2 : The CP lifetime is equal to the breeding blankets lifetime
    3 : The CP lifetime is equal to the plant lifetime
    """

    cowner: Parameter[float] = 0.15
    """owner cost factor"""

    cplife_input: Parameter[float] = 2.0
    """User input full power year lifetime of the centrepost (years) (i_cp_lifetime = 0)"""

    cplife: Parameter[float] = 0.0
    """Calculated full power year lifetime of centrepost (years)"""

    cplife_cal: Parameter[float] = 0.0
    """Calculated calendar year lifetime of centrepost (years)"""

    cpstcst: Parameter[float] = 0.0
    """ST centrepost direct cost (M$)"""

    cpstflnc: Parameter[float] = 10.0
    """allowable ST centrepost neutron fluence (MW-yr/m2)"""

    crctcore: Parameter[float] = 0.0
    """reactor core costs (categories 221, 222 and 223)"""

    csi: Parameter[float] = 16.0
    """allowance for site costs (M$)"""

    cturbb: Parameter[float] = 38.0
    """cost of turbine building (M$)"""

    decomf: Parameter[float] = 0.1
    """proportion of constructed cost required for decommissioning fund"""

    dintrt: Parameter[float] = 0.0
    """diff between borrowing and saving interest rates"""

    divcst: Parameter[float] = 0.0
    """divertor direct cost (M$)"""

    life_div_fpy: Parameter[float] = 0.0
    """Full power year lifetime of divertor (fpy)"""

    life_div: Parameter[float] = 0.0
    """Calendar year lifetime of divertor (y)"""

    dtlife: Parameter[float] = 0.0
    """period prior to the end of the plant life that the decommissioning fund is used (years)"""

    fcap0: Parameter[float] = 1.165
    """average cost of money for construction of plant assuming design/construction time of six years"""

    fcap0cp: Parameter[float] = 1.08
    """average cost of money for replaceable components assuming lead time for these of two years"""

    fcdfuel: Parameter[float] = 0.1
    """fraction of current drive cost treated as fuel (if `ifueltyp = 1`)"""

    fcontng: Parameter[float] = 0.195
    """project contingency factor"""

    fcr0: Parameter[float] = 0.0966
    """fixed charge rate during construction"""

    fkind: Parameter[float] = 1.0
    """multiplier for Nth of a kind costs"""

    fwallcst: Parameter[float] = 0.0
    """first wall cost (M$)"""

    i_plant_availability: int = 2
    """Switch for plant availability model:
    - =0 use input value for f_t_plant_available
    - =1 calculate f_t_plant_available using Taylor and Ward 1999 model
    - =2 calculate f_t_plant_available using new (2015) model
    - =3 calculate f_t_plant_available using ST model
    """

    ibkt_life: int = 0
    """Switch for fw/blanket lifetime calculation in availability module:
    - =0 use neutron fluence model
    - =1 use fusion power model (DEMO only)"""

    life_dpa: Parameter[float] = 50
    """Allowable DPA from DEMO fw/blanket lifetime calculation in availability module"""

    bktcycles: Parameter[float] = 1.0e3
    """Number of fusion cycles to reach allowable DPA from DEMO fw/blanket lifetime calculation"""

    f_t_plant_available_min: Parameter[float] = 0.75
    """Minimum allowed plant availability (`constraint equation 61`)"""

    tok_build_cost_per_vol: Parameter[float] = 1283.0
    """Unit cost for tokamak complex buildings, including building and site services ($/m3)"""

    light_build_cost_per_vol: Parameter[float] = 270.0
    """Unit cost for unshielded non-active buildings ($/m3)"""

    num_rh_systems: int = 4
    """Number of remote handling systems (1-10)"""

    conf_mag: Parameter[float] = 0.99
    """c parameter, which determines the temperature margin at which magnet lifetime starts to decline"""

    div_prob_fail: Parameter[float] = 0.0002
    """Divertor probability of failure (per op day)"""

    div_umain_time: Parameter[float] = 0.25
    """Divertor unplanned maintenance time (years)"""

    div_nref: Parameter[float] = 7000.0
    """Reference value for cycle cycle life of divertor"""

    div_nu: Parameter[float] = 14000.0
    """The cycle when the divertor fails with 100% probability"""

    fwbs_nref: Parameter[float] = 20000.0
    """Reference value for cycle life of blanket"""

    fwbs_nu: Parameter[float] = 40000.0
    """The cycle when the blanket fails with 100% probability"""

    fwbs_prob_fail: Parameter[float] = 0.0002
    """Fwbs probability of failure (per op day)"""

    fwbs_umain_time: Parameter[float] = 0.25
    """Fwbs unplanned maintenance time (years)"""

    redun_vacp: Parameter[float] = 25.0
    """Vacuum system pump redundancy level (%)"""

    redun_vac: int = 0
    """Number of redundant vacuum pumps"""

    t_plant_operational_total_yrs: Parameter[float] = 0.0
    """Operational time (yrs)"""

    t_blkt_replace_yrs: Parameter[float] = 0.5
    """time taken to replace blanket (y) (`i_plant_availability=1`)"""

    tcomrepl: Parameter[float] = 0.5
    """time taken to replace both blanket and divertor (y) (`i_plant_availability=1`)"""

    t_div_replace_yrs: Parameter[float] = 0.25
    """Time taken to replace divertor (y) (`i_plant_availability=1`)"""

    uubop: Parameter[float] = 0.02
    """unplanned unavailability factor for balance of plant (`i_plant_availability=1`)"""

    uucd: Parameter[float] = 0.02
    """unplanned unavailability factor for current drive (`i_plant_availability=1`)"""

    uudiv: Parameter[float] = 0.04
    """unplanned unavailability factor for divertor (`i_plant_availability=1`)"""

    uufuel: Parameter[float] = 0.02
    """unplanned unavailability factor for fuel system (`i_plant_availability=1`)"""

    uufw: Parameter[float] = 0.04
    """unplanned unavailability factor for first wall (`i_plant_availability=1`)"""

    uumag: Parameter[float] = 0.02
    """unplanned unavailability factor for magnets (`i_plant_availability=1`)"""

    uuves: Parameter[float] = 0.04
    """unplanned unavailability factor for vessel (`i_plant_availability=1`)"""

    ifueltyp: int = 0
    """Switch for fuel type:
    - =2 treat initial blanket, divertor, first wall
    as capital costs. Treat all later items and
    fraction fcdfuel of CD equipment as fuel costs
    - =1 treat blanket divertor, first wall and
    fraction fcdfuel of CD equipment as fuel cost
    - =0 treat these as capital cost
    """

    ipnet: int = 0
    """Switch for net electric power calculation:
    - =0 scale so that always > 0
    - =1 let go < 0 (no c-o-e)
    """

    ireactor: int = 1
    """Switch for net electric power and cost of electricity calculations:
    - =0 do not calculate MW(electric) or c-o-e
    - =1 calculate MW(electric) and c-o-e
    """

    lsa: int = 4
    """Level of safety assurance switch (generally, use 3 or 4):
    - =1 truly passively safe plant
    - =2,3 in-between
    - =4 like current fission plant
    """

    moneyint: Parameter[float] = 0.0
    """interest portion of capital cost (M$)"""

    output_costs: int = 1
    """Switch for costs output:
    - =0 do not write cost-related outputs to file
    - =1 write cost-related outputs to file
    """

    discount_rate: Parameter[float] = 0.0435
    """effective cost of money in constant dollars"""

    startupratio: Parameter[float] = 1.0
    """ratio of additional HCD power for start-up to flat-top operational requirements"""

    startuppwr: Parameter[float] = 0.0
    """cost associated with additional HCD system power required on start-up ($)"""

    supercond_cost_model: int = 0
    """Switch for superconductor cost model:
    - =0 use $/kg
    - =1 use $/kAm"""

    life_plant: Parameter[float] = 30.0
    """Full power year plant lifetime (years)"""

    tmain: Parameter[float] = None
    """Maintenance time for replacing CP (years) (i_plant_availability = 3)"""

    u_unplanned_cp: Parameter[float] = 0.0
    """User-input CP unplanned unavailability (i_plant_availability = 3)"""

    UCAD: Parameter[float] = 180.0
    """unit cost for administration buildings (M$/m3)"""

    UCAF: Parameter[float] = 1.5e6
    """unit cost for aux facility power equipment ($)"""

    UCAHTS: Parameter[float] = 31.0
    """unit cost for aux heat transport equipment ($/W**exphts)"""

    UCAP: Parameter[float] = 17.0
    """unit cost of auxiliary transformer ($/kVA)"""

    ucblbe: Parameter[float] = 260.0
    """unit cost for blanket beryllium ($/kg)"""

    ucblbreed: Parameter[float] = 875.0
    """unit cost for breeder material ($/kg) (`blktmodel>0`)"""

    ucblli: Parameter[float] = 875.0
    """unit cost for blanket lithium ($/kg) (30% Li6)"""

    ucblli2o: Parameter[float] = 600.0
    """unit cost for blanket Li_2O ($/kg)"""

    ucbllipb: Parameter[float] = 10.3
    """unit cost for blanket Li-Pb ($/kg) (30% Li6)"""

    ucblss: Parameter[float] = 90.0
    """unit cost for blanket stainless steel ($/kg)"""

    ucblvd: Parameter[float] = 200.0
    """unit cost for blanket vanadium ($/kg)"""

    UCBPMP: Parameter[float] = 2.925e5
    """vacuum system backing pump cost ($)"""

    ucbus: Parameter[float] = 0.123
    """cost of aluminium bus for TF coil ($/A-m)"""

    uccase: Parameter[float] = 50.0
    """cost of superconductor case ($/kg)"""

    UCCO: Parameter[float] = 350.0
    """unit cost for control buildings (M$/m3)"""

    uccpcl1: Parameter[float] = 250.0
    """cost of high strength tapered copper ($/kg)"""

    uccpclb: Parameter[float] = 150.0
    """cost of TF outboard leg plate coils ($/kg)"""

    UCCPMP: Parameter[float] = 3.9e5
    """vacuum system cryopump cost ($)"""

    UCCR: Parameter[float] = 460.0
    """unit cost for cryogenic building (M$/vol)"""

    uccry: Parameter[float] = 9.3e4
    """heat transport system cryoplant costs ($/W**expcry)"""

    uccryo: Parameter[float] = 32.0
    """unit cost for vacuum vessel ($/kg)"""

    uccu: Parameter[float] = 75.0
    """unit cost for copper in superconducting cable ($/kg)"""

    UCDGEN: Parameter[float] = 1.7e6
    """cost per 8 MW diesel generator ($)"""

    ucdiv: Parameter[float] = 2.8e5
    """cost of divertor blade ($)"""

    UCDTC: Parameter[float] = 13.0
    """detritiation, air cleanup cost ($/10000m3/hr)"""

    UCDUCT: Parameter[float] = 4.225e4
    """vacuum system duct cost ($/m)"""

    ucech: Parameter[float] = 3.0
    """ECH system cost ($/W)"""

    UCEL: Parameter[float] = 380.0
    """unit cost for electrical equipment building (M$/m3)"""

    UCES1: Parameter[float] = 3.2e4
    """MGF (motor-generator flywheel) cost factor ($/MVA**0.8)"""

    UCES2: Parameter[float] = 8.8e3
    """MGF (motor-generator flywheel) cost factor ($/MJ**0.8)"""

    ucf1: Parameter[float] = 2.23e7
    """cost of fuelling system ($)"""

    ucfnc: Parameter[float] = 35.0
    """outer PF coil fence support cost ($/kg)"""

    UCFPR: Parameter[float] = 4.4e7
    """cost of 60g/day tritium processing unit ($)"""

    ucfuel: Parameter[float] = 3.45
    """unit cost of D-T fuel (M$/year/1200MW)"""

    UCFWA: Parameter[float] = 6.0e4
    """first wall armour cost ($/m2)"""

    UCFWPS: Parameter[float] = 1.0e7
    """first wall passive stabiliser cost ($)"""

    UCFWS: Parameter[float] = 5.3e4
    """first wall structure cost ($/m2)"""

    UCGSS: Parameter[float] = 35.0
    """cost of reactor structure ($/kg)"""

    uche3: Parameter[float] = 1.0e6
    """cost of helium-3 ($/kg)"""

    uchrs: Parameter[float] = 87.9e6
    """cost of heat rejection system ($)"""

    uchts: Parameter[list[float]] = field(default_factory=lambda: [15.3, 19.1])
    """cost of heat transport system equipment per loop ($/W); dependent on coolant type (i_blkt_coolant_type)"""

    uciac: Parameter[float] = 1.5e8
    """cost of instrumentation, control & diagnostics ($)"""

    ucich: Parameter[float] = 3.0
    """ICH system cost ($/W)"""

    UCINT: Parameter[float] = 35.0
    """superconductor intercoil structure cost ($/kg)"""

    uclh: Parameter[float] = 3.3
    """lower hybrid system cost ($/W)"""

    UCLV: Parameter[float] = 16.0
    """low voltage system cost ($/kVA)"""

    UCMB: Parameter[float] = 260.0
    """unit cost for reactor maintenance building (M$/m3)"""

    ucme: Parameter[float] = 1.25e8
    """cost of maintenance equipment ($)"""

    ucmisc: Parameter[float] = 2.5e7
    """miscellaneous plant allowance ($)"""

    ucnbi: Parameter[float] = 3.3
    """NBI system cost ($/W)"""

    UCNBV: Parameter[float] = 1000.0
    """cost of nuclear building ventilation ($/m3)"""

    ucoam: Parameter[list[float]] = field(
        default_factory=lambda: [68.8, 68.8, 68.8, 74.4]
    )
    """annual cost of operation and maintenance (M$/year/1200MW**0.5)"""

    ucpens: Parameter[float] = 32.0
    """penetration shield cost ($/kg)"""

    ucpfb: Parameter[float] = 210.0
    """cost of PF coil buses ($/kA-m)"""

    ucpfbk: Parameter[float] = 1.66e4
    """cost of PF coil DC breakers ($/MVA**0.7)"""

    ucpfbs: Parameter[float] = 4.9e3
    """cost of PF burn power supplies ($/kW**0.7)"""

    ucpfcb: Parameter[float] = 7.5e4
    """cost of PF coil AC breakers ($/circuit)"""

    ucpfdr1: Parameter[float] = 150.0
    """cost factor for dump resistors ($/MJ)"""

    ucpfic: Parameter[float] = 1.0e4
    """cost of PF instrumentation and control ($/channel)"""

    ucpfps: Parameter[float] = 3.5e4
    """cost of PF coil pulsed power supplies ($/MVA)"""

    UCPHX: Parameter[float] = 15.0
    """primary heat transport cost ($/W**exphts)"""

    UCPP: Parameter[float] = 48.0
    """cost of primary power transformers ($/kVA**0.9)"""

    ucrb: Parameter[float] = 400.0
    """cost of reactor building (M$/m3)"""

    ucsc: Parameter[list[float]] = field(
        default_factory=lambda: [
            600.0,
            600.0,
            300.0,
            600.0,
            600.0,
            600.0,
            300.0,
            1200.0,
            1200.0,
        ]
    )
    """cost of superconductor ($/kg)"""

    sc_mat_cost_0: Parameter[list[float]] = field(
        default_factory=lambda: [4.8, 2.0, 1.0, 4.8, 4.8, 47.4, 1.0, 47.4, 47.4]
    )
    """cost of superconductor ($/kA m) at 6.4 T, 4.2 K"""

    UCSH: Parameter[float] = 115.0
    """cost of shops and warehouses (M$/m3)"""

    ucshld: Parameter[float] = 32.0
    """cost of shield structural steel ($/kg)"""

    UCSWYD: Parameter[float] = 1.84e7
    """switchyard equipment costs ($)"""

    uctfbr: Parameter[float] = 1.22
    """cost of TF coil breakers ($/W**0.7)"""

    uctfbus: Parameter[float] = 100.0
    """cost of TF coil bus ($/kg)"""

    UCTFDR: Parameter[float] = 1.75e-4
    """cost of TF coil dump resistors ($/J)"""

    UCTFGR: Parameter[float] = 5000.0
    """additional cost of TF coil dump resistors ($/coil)"""

    UCTFIC: Parameter[float] = 1.0e4
    """cost of TF coil instrumentation and control ($/coil/30)"""

    uctfps: Parameter[float] = 24.0
    """cost of TF coil power supplies ($/W**0.7)"""

    uctfsw: Parameter[float] = 1.0
    """cost of TF coil slow dump switches ($/A)"""

    UCTPMP: Parameter[float] = 1.105e5
    """cost of turbomolecular pump ($)"""

    UCTR: Parameter[float] = 370.0
    """cost of tritium building ($/m3)"""

    ucturb: Parameter[list[float]] = field(default_factory=lambda: [230.0e6, 245.0e6])
    """cost of turbine plant equipment ($) (dependent on coolant type i_blkt_coolant_type)"""

    UCVALV: Parameter[float] = 3.9e5
    """vacuum system valve cost ($)"""

    UCVDSH: Parameter[float] = 26.0
    """vacuum duct shield cost ($/kg)"""

    UCVIAC: Parameter[float] = 1.3e6
    """vacuum system instrumentation and control cost ($)"""

    ucwindpf: Parameter[float] = 465.0
    """cost of PF coil superconductor windings ($/m)"""

    ucwindtf: Parameter[float] = 480.0
    """cost of TF coil superconductor windings ($/m)"""

    UCWS: Parameter[float] = 460.0
    """cost of active assembly shop ($/m3)"""

    ucwst: Parameter[list[float]] = field(
        default_factory=lambda: [0.0, 3.94, 5.91, 7.88]
    )
    """cost of waste disposal (M$/y/1200MW)"""


CREATE_DICTS_FROM_DATACLASS = CostData
