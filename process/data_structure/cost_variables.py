from dataclasses import dataclass, field


@dataclass(slots=True)
class CostData:
    c228: float = 0.0
    """c228 account cost - Instrumentation and control cost (M$)"""

    c229: float = 0.0
    """c229 account cost - Maintenance equipment cost (M$)"""

    c23: float = 0.0
    """c23 account cost - Turbine plant equipment cost (M$)"""

    c25: float = 0.0
    """c25 account cost - Miscellaneous plant equipment cost (M$)"""

    c26: float = 0.0
    """c26 account cost - Heat rejection system cost (M$)"""

    cindrt: float = 0.0
    """Indirect cost (M$)"""

    ccont: float = 0.0
    """Total contingency (M$)"""

    c226: float = 0.0
    """c226 account cost - Heat transport system (M$)"""

    c2261: float = 0.0
    """c2261 account cost - Reactor cooling system (M$)"""

    c2262: float = 0.0
    """c2262 account cost - Auxiliary component cooling (M$)"""

    c2263: float = 0.0
    """c2263 account cost - Cryogenic system (M$)"""

    c227: float = 0.0
    """c227 account cost - Fuel handling (M$)"""

    c2271: float = 0.0
    """c2271 account cost - Fuelling system (M$)"""

    c2272: float = 0.0
    """c2272 account cost - Fuel processing and purification (M$)"""

    c2273: float = 0.0
    """c2273 account cost - Atmospheric recovery systems (M$)"""

    c2274: float = 0.0
    """c2274 account cost - Nuclear building ventilation (M$)"""

    c24: float = 0.0
    """c24 account cost - Electric plant equipment (M$)"""

    c241: float = 0.0
    """c241 account cost - Electric plant equipment - switchyard (M$)"""

    c242: float = 0.0
    """c242 account cost - Electric plant equipment - Transformers (M$)"""

    c243: float = 0.0
    """c243 account cost - Electric plant equipment - Low voltage (M$)"""

    c244: float = 0.0
    """c244 account cost - Electric plant equipment - Diesel generators (M$)"""

    c245: float = 0.0
    """c245 account cost - Electric plant equipment - Aux facility power (M$)"""

    c21: float = 0.0
    """c21 account cost - Structures and site facilities (M$)"""

    c211: float = 0.0
    """c211 account cost - Site improvements, facilities and land (M$)"""

    c212: float = 0.0
    """c212 account cost - Reactor building (M$)"""

    c213: float = 0.0
    """c213 account cost - Turbine building (M$)"""

    c214: float = 0.0
    """c214 account cost - Reactor maintenance and warm shops buildings (M$)"""

    c2141: float = 0.0
    """c2141 account cost - Reactor maintenance building cost (M$)"""

    c2142: float = 0.0
    """c2142 account cost - Warm shop cost (M$) (M$)"""

    c215: float = 0.0
    """c215 account cost - Tritium building cost (M$)"""

    c216: float = 0.0
    """c216 account cost - Electrical equipment building cost (M$)"""

    c217: float = 0.0
    """c217 account cost - Other buildings (M$)"""

    c2171: float = 0.0
    """c2171 account cost - Additional buildings cost (M$)"""

    c2172: float = 0.0
    """c2172 account cost - Control room buildings cost (M$)"""

    c2173: float = 0.0
    """c2173 account cost - Shop and warehouses cost (M$)"""

    c2174: float = 0.0
    """c2174 account cost - Cryogenic building cost (M$)"""

    c22: float = 0.0
    """c22 account cost - Fusion power island (M$)"""

    c2211: float = 0.0
    """c2211 account cost - First wall (M$)"""

    c2212: float = 0.0
    """c2212 account cost - Blanket (M$)"""

    c22121: float = 0.0
    """c22121 account cost - Blanket beryllium cost (M$)"""

    c22122: float = 0.0
    """c22122 account cost - Blanket breeder material cost (M$)"""

    c22123: float = 0.0
    """c22123 account cost - Blanket stainless steel cost (M$)"""

    c22124: float = 0.0
    """c22124 account cost - Blanket vanadium cost (M$)"""

    c22125: float = 0.0
    """c22125 account cost - Blanket carbon cloth cost (M$)"""

    c22126: float = 0.0
    """c22126 account cost - Blanket concrete cost (M$)"""

    c22127: float = 0.0
    """c22127 account cost - Blanket FLiBe cost (M$)"""

    c22128: float = 0.0
    """c22128 account cost - Blanket lithium cost (M$)"""

    c2213: float = 0.0
    """c2213 account cost - Shield (M$)"""

    c22131: float = 0.0
    """c22131 account cost - Bulk shield cost (M$)"""

    c22132: float = 0.0
    """c22132 account cost - Penetration shielding cost (M$)"""

    c2214: float = 0.0
    """c2214 account cost - Total support structure cost (M$)"""

    c2215: float = 0.0
    """c2215 account cost - TF coil gravity support structure (M$)"""

    c2221: float = 0.0
    """c2221 account cost - TF magnet assemblies cost (M$)"""

    c22211: float = 0.0
    """c22211 account cost - TF coil conductor cost (M$)"""

    c22212: float = 0.0
    """c22212 account cost - TF coil winding cost (M$)"""

    c22213: float = 0.0
    """c22213 account cost - TF coil case cost (M$)"""

    c22214: float = 0.0
    """c22214 account cost - TF intercoil structure cost (M$)"""

    c22215: float = 0.0
    """c22215 account cost - TF coil gravity support structure (M$)"""

    c2222: float = 0.0
    """c2222 account cost - PF magnet assemblies (M$)"""

    c22221: float = 0.0
    """c22221 account cost - PF coil conductor cost (M$)"""

    c22222: float = 0.0
    """c22222 account cost - PF coil winding cost (M$)"""

    c22223: float = 0.0
    """c22223 account cost - PF coil case cost (M$)"""

    c22224: float = 0.0
    """c22224 account cost - PF coil support structure cost (M$)"""

    c2223: float = 0.0
    """c2223 account cost - Vacuum vessel assembly cost (M$)"""

    c223: float = 0.0
    """c223 account cost - Power injection (M$)"""

    c2231: float = 0.0
    """c2231 account cost - ECH system cost (M$)"""

    c2232: float = 0.0
    """c2232 account cost - Lower hybrid system cost (M$)"""

    c2233: float = 0.0
    """c2233 account cost - Neutral beam system cost (M$)"""

    c2234: float = 0.0
    """c2234 account cost - ? (M$)"""

    c224: float = 0.0
    """c224 account cost - Vacuum system (M$)"""

    c2241: float = 0.0
    """c2241 account cost - High vacuum pumps cost (M$)"""

    c2242: float = 0.0
    """c2242 account cost - Backing pumps cost (M$)"""

    c2243: float = 0.0
    """c2243 account cost - Vacuum duct cost (M$)"""

    c2244: float = 0.0
    """c2244 account cost - Valves cost (M$)"""

    c2245: float = 0.0
    """c2245 account cost - Duct shielding cost (M$)"""

    c2246: float = 0.0
    """c2246 account cost - Instrumentation cost (M$)"""

    c225: float = 0.0
    """c225 account cost - Power conditioning (M$)"""

    c2251: float = 0.0
    """c2251 account cost - TF coil power conditioning (M$)"""

    c22511: float = 0.0
    """c22511 account cost - TF coil power supplies cost (M$)"""

    c22512: float = 0.0
    """c22512 account cost - TF coil breakers cost (M$)"""

    c22513: float = 0.0
    """c22513 account cost - TF coil dump resistors cost (M$)"""

    c22514: float = 0.0
    """c22514 account cost - TF coil instrumentation and control (M$)"""

    c22515: float = 0.0
    """c22515 account cost - TF coil bussing cost (M$)"""

    c2252: float = 0.0
    """c2252 account cost - PF coil power costs (M$)"""

    c22521: float = 0.0
    """c22521 account cost - PF coil power supplies cost (M$)"""

    c22522: float = 0.0
    """c22522 account cost - PF coil instrumentation and control (M$)"""

    c22523: float = 0.0
    """c22523 account cost - PF coil bussing cost (M$)"""

    c22524: float = 0.0
    """c22524 account cost - PF coil burn power supplies cost (M$)"""

    c22525: float = 0.0
    """c22525 account cost - PF coil breakers cost (M$)"""

    c22526: float = 0.0
    """c22526 account cost - PF coil dump resistors cost (M$)"""

    c22527: float = 0.0
    """c account cost - PF coil ac breakers cost (M$)"""

    c2253: float = 0.0
    """c2253 account cost - Energy storage cost (M$)"""

    chx: float = 0.0
    """Primary heat exchanger cost (M$)"""

    cpp: float = 0.0
    """Pumps and piping system cost (M$)"""

    cppa: float = 0.0
    """Pumps, piping cost (M$)"""

    abktflnc: float = 5.0
    """allowable first wall/blanket neutron fluence (MW-yr/m2) (`blktmodel=0`)"""

    adivflnc: float = 7.0
    """allowable divertor heat fluence (MW-yr/m2)"""

    blkcst: float = 0.0
    """blanket direct cost (M$)"""

    c221: float = 0.0
    """total account 221 cost (M$) - first wall, blanket, shield, support structure and div plates"""

    c222: float = 0.0
    """total account 222 cost (M$) - TF coils + PF coils"""

    capcost: float = 0.0
    """total capital cost including interest (M$)"""

    cconfix: float = 80.0
    """fixed cost of superconducting cable ($/m)"""

    cconshpf: float = 70.0
    """cost of PF coil steel conduit/sheath ($/m)"""

    cconshtf: float = 75.0
    """cost of TF coil steel conduit/sheath ($/m)"""

    cdcost: float = 0.0
    """current drive direct costs (M$)"""

    cdirt: float = 0.0
    """total plant direct cost (M$)"""

    life_hcd_fpy: float = 0.0
    """Full power year lifetime of heating/current drive system (y)"""

    cdrlife_cal: float = 0.0
    """Calendar year lifetime of heating/current drive system (y)"""

    f_t_plant_available: float = 0.75
    """Total plant availability fraction; input if `i_plant_availability=0`"""

    cpfact: float = 0.0
    """Total plant capacity factor"""

    cfind: list[float] = field(default_factory=lambda: [0.244, 0.244, 0.244, 0.29])
    """indirect cost factor (func of lsa) (cost model = 0)"""

    cland: float = 19.2
    """cost of land (M$)"""

    coe: float = 0.0
    """cost of electricity ($/MW-hr)"""

    coecap: float = 0.0
    """capital cost of electricity (m$/kW-hr)"""

    coefuelt: float = 0.0
    """'fuel' (including replaceable components) contribution to cost of electricity (m$/kW-hr)"""

    coeoam: float = 0.0
    """operation and maintenance contribution to cost of electricity (m$/kW-hr)"""

    concost: float = 0.0
    """plant construction cost (M$)"""

    costexp: float = 0.8
    """cost exponent for scaling in 2015 costs model"""

    costexp_pebbles: float = 0.6
    """cost exponent for pebbles in 2015 costs model"""

    cost_factor_buildings: float = 1.0
    """cost scaling factor for buildings"""

    cost_factor_land: float = 1.0
    """cost scaling factor for land"""

    cost_factor_tf_coils: float = 1.0
    """cost scaling factor for TF coils"""

    cost_factor_fwbs: float = 1.0
    """cost scaling factor for fwbs"""

    cost_factor_rh: float = 1.0
    """cost scaling factor for remote handling"""

    cost_factor_vv: float = 1.0
    """cost scaling factor for vacuum vessel"""

    cost_factor_bop: float = 1.0
    """cost scaling factor for energy conversion system"""

    cost_factor_misc: float = 1.0
    """cost scaling factor for remaining subsystems"""

    maintenance_fwbs: float = 0.2
    """Maintenance cost factor: first wall, blanket, shield, divertor"""

    maintenance_gen: float = 0.05
    """Maintenance cost factor: All other components except coils, vacuum vessel,
    thermal shield, cryostat, land"""

    amortization: float = 13.6
    """amortization factor (fixed charge factor) "A" (years)"""

    cost_model: int = 1
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

    cowner: float = 0.15
    """owner cost factor"""

    cplife_input: float = 2.0
    """User input full power year lifetime of the centrepost (years) (i_cp_lifetime = 0)"""

    cplife: float = 0.0
    """Calculated full power year lifetime of centrepost (years)"""

    cplife_cal: float = 0.0
    """Calculated calendar year lifetime of centrepost (years)"""

    cpstcst: float = 0.0
    """ST centrepost direct cost (M$)"""

    cpstflnc: float = 10.0
    """allowable ST centrepost neutron fluence (MW-yr/m2)"""

    crctcore: float = 0.0
    """reactor core costs (categories 221, 222 and 223)"""

    csi: float = 16.0
    """allowance for site costs (M$)"""

    cturbb: float = 38.0
    """cost of turbine building (M$)"""

    decomf: float = 0.1
    """proportion of constructed cost required for decommissioning fund"""

    dintrt: float = 0.0
    """diff between borrowing and saving interest rates"""

    divcst: float = 0.0
    """divertor direct cost (M$)"""

    life_div_fpy: float = 0.0
    """Full power year lifetime of divertor (fpy)"""

    life_div: float = 0.0
    """Calendar year lifetime of divertor (y)"""

    dtlife: float = 0.0
    """period prior to the end of the plant life that the decommissioning fund is used (years)"""

    fcap0: float = 1.165
    """average cost of money for construction of plant assuming design/construction time of six years"""

    fcap0cp: float = 1.08
    """average cost of money for replaceable components assuming lead time for these of two years"""

    fcdfuel: float = 0.1
    """fraction of current drive cost treated as fuel (if `ifueltyp = 1`)"""

    fcontng: float = 0.195
    """project contingency factor"""

    fcr0: float = 0.0966
    """fixed charge rate during construction"""

    fkind: float = 1.0
    """multiplier for Nth of a kind costs"""

    fwallcst: float = 0.0
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

    life_dpa: float = 50
    """Allowable DPA from DEMO fw/blanket lifetime calculation in availability module"""

    bktcycles: float = 1.0e3
    """Number of fusion cycles to reach allowable DPA from DEMO fw/blanket lifetime calculation"""

    f_t_plant_available_min: float = 0.75
    """Minimum allowed plant availability (`constraint equation 61`)"""

    tok_build_cost_per_vol: float = 1283.0
    """Unit cost for tokamak complex buildings, including building and site services ($/m3)"""

    light_build_cost_per_vol: float = 270.0
    """Unit cost for unshielded non-active buildings ($/m3)"""

    num_rh_systems: int = 4
    """Number of remote handling systems (1-10)"""

    conf_mag: float = 0.99
    """c parameter, which determines the temperature margin at which magnet lifetime starts to decline"""

    div_prob_fail: float = 0.0002
    """Divertor probability of failure (per op day)"""

    div_umain_time: float = 0.25
    """Divertor unplanned maintenance time (years)"""

    div_nref: float = 7000.0
    """Reference value for cycle cycle life of divertor"""

    div_nu: float = 14000.0
    """The cycle when the divertor fails with 100% probability"""

    fwbs_nref: float = 20000.0
    """Reference value for cycle life of blanket"""

    fwbs_nu: float = 40000.0
    """The cycle when the blanket fails with 100% probability"""

    fwbs_prob_fail: float = 0.0002
    """Fwbs probability of failure (per op day)"""

    fwbs_umain_time: float = 0.25
    """Fwbs unplanned maintenance time (years)"""

    redun_vacp: float = 25.0
    """Vacuum system pump redundancy level (%)"""

    redun_vac: int = 0
    """Number of redundant vacuum pumps"""

    t_plant_operational_total_yrs: float = 0.0
    """Operational time (yrs)"""

    t_blkt_replace_yrs: float = 0.5
    """time taken to replace blanket (y) (`i_plant_availability=1`)"""

    tcomrepl: float = 0.5
    """time taken to replace both blanket and divertor (y) (`i_plant_availability=1`)"""

    t_div_replace_yrs: float = 0.25
    """Time taken to replace divertor (y) (`i_plant_availability=1`)"""

    uubop: float = 0.02
    """unplanned unavailability factor for balance of plant (`i_plant_availability=1`)"""

    uucd: float = 0.02
    """unplanned unavailability factor for current drive (`i_plant_availability=1`)"""

    uudiv: float = 0.04
    """unplanned unavailability factor for divertor (`i_plant_availability=1`)"""

    uufuel: float = 0.02
    """unplanned unavailability factor for fuel system (`i_plant_availability=1`)"""

    uufw: float = 0.04
    """unplanned unavailability factor for first wall (`i_plant_availability=1`)"""

    uumag: float = 0.02
    """unplanned unavailability factor for magnets (`i_plant_availability=1`)"""

    uuves: float = 0.04
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

    moneyint: float = 0.0
    """interest portion of capital cost (M$)"""

    output_costs: int = 1
    """Switch for costs output:
    - =0 do not write cost-related outputs to file
    - =1 write cost-related outputs to file
    """

    discount_rate: float = 0.0435
    """effective cost of money in constant dollars"""

    startupratio: float = 1.0
    """ratio of additional HCD power for start-up to flat-top operational requirements"""

    startuppwr: float = 0.0
    """cost associated with additional HCD system power required on start-up ($)"""

    supercond_cost_model: int = 0
    """Switch for superconductor cost model:
    - =0 use $/kg
    - =1 use $/kAm"""

    life_plant: float = 30.0
    """Full power year plant lifetime (years)"""

    tmain: float = None
    """Maintenance time for replacing CP (years) (i_plant_availability = 3)"""

    u_unplanned_cp: float = 0.0
    """User-input CP unplanned unavailability (i_plant_availability = 3)"""

    UCAD: float = 180.0
    """unit cost for administration buildings (M$/m3)"""

    UCAF: float = 1.5e6
    """unit cost for aux facility power equipment ($)"""

    UCAHTS: float = 31.0
    """unit cost for aux heat transport equipment ($/W**exphts)"""

    UCAP: float = 17.0
    """unit cost of auxiliary transformer ($/kVA)"""

    ucblbe: float = 260.0
    """unit cost for blanket beryllium ($/kg)"""

    ucblbreed: float = 875.0
    """unit cost for breeder material ($/kg) (`blktmodel>0`)"""

    ucblli: float = 875.0
    """unit cost for blanket lithium ($/kg) (30% Li6)"""

    ucblli2o: float = 600.0
    """unit cost for blanket Li_2O ($/kg)"""

    ucbllipb: float = 10.3
    """unit cost for blanket Li-Pb ($/kg) (30% Li6)"""

    ucblss: float = 90.0
    """unit cost for blanket stainless steel ($/kg)"""

    ucblvd: float = 200.0
    """unit cost for blanket vanadium ($/kg)"""

    UCBPMP: float = 2.925e5
    """vacuum system backing pump cost ($)"""

    ucbus: float = 0.123
    """cost of aluminium bus for TF coil ($/A-m)"""

    uccase: float = 50.0
    """cost of superconductor case ($/kg)"""

    UCCO: float = 350.0
    """unit cost for control buildings (M$/m3)"""

    uccpcl1: float = 250.0
    """cost of high strength tapered copper ($/kg)"""

    uccpclb: float = 150.0
    """cost of TF outboard leg plate coils ($/kg)"""

    UCCPMP: float = 3.9e5
    """vacuum system cryopump cost ($)"""

    UCCR: float = 460.0
    """unit cost for cryogenic building (M$/vol)"""

    uccry: float = 9.3e4
    """heat transport system cryoplant costs ($/W**expcry)"""

    uccryo: float = 32.0
    """unit cost for vacuum vessel ($/kg)"""

    uccu: float = 75.0
    """unit cost for copper in superconducting cable ($/kg)"""

    UCDGEN: float = 1.7e6
    """cost per 8 MW diesel generator ($)"""

    ucdiv: float = 2.8e5
    """cost of divertor blade ($)"""

    UCDTC: float = 13.0
    """detritiation, air cleanup cost ($/10000m3/hr)"""

    UCDUCT: float = 4.225e4
    """vacuum system duct cost ($/m)"""

    ucech: float = 3.0
    """ECH system cost ($/W)"""

    UCEL: float = 380.0
    """unit cost for electrical equipment building (M$/m3)"""

    UCES1: float = 3.2e4
    """MGF (motor-generator flywheel) cost factor ($/MVA**0.8)"""

    UCES2: float = 8.8e3
    """MGF (motor-generator flywheel) cost factor ($/MJ**0.8)"""

    ucf1: float = 2.23e7
    """cost of fuelling system ($)"""

    ucfnc: float = 35.0
    """outer PF coil fence support cost ($/kg)"""

    UCFPR: float = 4.4e7
    """cost of 60g/day tritium processing unit ($)"""

    ucfuel: float = 3.45
    """unit cost of D-T fuel (M$/year/1200MW)"""

    UCFWA: float = 6.0e4
    """first wall armour cost ($/m2)"""

    UCFWPS: float = 1.0e7
    """first wall passive stabiliser cost ($)"""

    UCFWS: float = 5.3e4
    """first wall structure cost ($/m2)"""

    UCGSS: float = 35.0
    """cost of reactor structure ($/kg)"""

    uche3: float = 1.0e6
    """cost of helium-3 ($/kg)"""

    uchrs: float = 87.9e6
    """cost of heat rejection system ($)"""

    uchts: list[float] = field(default_factory=lambda: [15.3, 19.1])
    """cost of heat transport system equipment per loop ($/W); dependent on coolant type (i_blkt_coolant_type)"""

    uciac: float = 1.5e8
    """cost of instrumentation, control & diagnostics ($)"""

    ucich: float = 3.0
    """ICH system cost ($/W)"""

    UCINT: float = 35.0
    """superconductor intercoil structure cost ($/kg)"""

    uclh: float = 3.3
    """lower hybrid system cost ($/W)"""

    UCLV: float = 16.0
    """low voltage system cost ($/kVA)"""

    UCMB: float = 260.0
    """unit cost for reactor maintenance building (M$/m3)"""

    ucme: float = 1.25e8
    """cost of maintenance equipment ($)"""

    ucmisc: float = 2.5e7
    """miscellaneous plant allowance ($)"""

    ucnbi: float = 3.3
    """NBI system cost ($/W)"""

    UCNBV: float = 1000.0
    """cost of nuclear building ventilation ($/m3)"""

    ucoam: list[float] = field(default_factory=lambda: [68.8, 68.8, 68.8, 74.4])
    """annual cost of operation and maintenance (M$/year/1200MW**0.5)"""

    ucpens: float = 32.0
    """penetration shield cost ($/kg)"""

    ucpfb: float = 210.0
    """cost of PF coil buses ($/kA-m)"""

    ucpfbk: float = 1.66e4
    """cost of PF coil DC breakers ($/MVA**0.7)"""

    ucpfbs: float = 4.9e3
    """cost of PF burn power supplies ($/kW**0.7)"""

    ucpfcb: float = 7.5e4
    """cost of PF coil AC breakers ($/circuit)"""

    ucpfdr1: float = 150.0
    """cost factor for dump resistors ($/MJ)"""

    ucpfic: float = 1.0e4
    """cost of PF instrumentation and control ($/channel)"""

    ucpfps: float = 3.5e4
    """cost of PF coil pulsed power supplies ($/MVA)"""

    UCPHX: float = 15.0
    """primary heat transport cost ($/W**exphts)"""

    UCPP: float = 48.0
    """cost of primary power transformers ($/kVA**0.9)"""

    ucrb: float = 400.0
    """cost of reactor building (M$/m3)"""

    ucsc: list[float] = field(
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

    sc_mat_cost_0: list[float] = field(
        default_factory=lambda: [4.8, 2.0, 1.0, 4.8, 4.8, 47.4, 1.0, 47.4, 47.4]
    )
    """cost of superconductor ($/kA m) at 6.4 T, 4.2 K"""

    UCSH: float = 115.0
    """cost of shops and warehouses (M$/m3)"""

    ucshld: float = 32.0
    """cost of shield structural steel ($/kg)"""

    UCSWYD: float = 1.84e7
    """switchyard equipment costs ($)"""

    uctfbr: float = 1.22
    """cost of TF coil breakers ($/W**0.7)"""

    uctfbus: float = 100.0
    """cost of TF coil bus ($/kg)"""

    UCTFDR: float = 1.75e-4
    """cost of TF coil dump resistors ($/J)"""

    UCTFGR: float = 5000.0
    """additional cost of TF coil dump resistors ($/coil)"""

    UCTFIC: float = 1.0e4
    """cost of TF coil instrumentation and control ($/coil/30)"""

    uctfps: float = 24.0
    """cost of TF coil power supplies ($/W**0.7)"""

    uctfsw: float = 1.0
    """cost of TF coil slow dump switches ($/A)"""

    UCTPMP: float = 1.105e5
    """cost of turbomolecular pump ($)"""

    UCTR: float = 370.0
    """cost of tritium building ($/m3)"""

    ucturb: list[float] = field(default_factory=lambda: [230.0e6, 245.0e6])
    """cost of turbine plant equipment ($) (dependent on coolant type i_blkt_coolant_type)"""

    UCVALV: float = 3.9e5
    """vacuum system valve cost ($)"""

    UCVDSH: float = 26.0
    """vacuum duct shield cost ($/kg)"""

    UCVIAC: float = 1.3e6
    """vacuum system instrumentation and control cost ($)"""

    ucwindpf: float = 465.0
    """cost of PF coil superconductor windings ($/m)"""

    ucwindtf: float = 480.0
    """cost of TF coil superconductor windings ($/m)"""

    UCWS: float = 460.0
    """cost of active assembly shop ($/m3)"""

    ucwst: list[float] = field(default_factory=lambda: [0.0, 3.94, 5.91, 7.88])
    """cost of waste disposal (M$/y/1200MW)"""


CREATE_DICTS_FROM_DATACLASS = CostData
