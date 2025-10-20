c228: float = None
"""c228 account cost"""

c229: float = None
"""c229 account cost"""

c23: float = None
"""c23 account cost"""

c25: float = None
"""c25 account cost"""

c26: float = None
"""c26 account cost"""

cindrt: float = None
"""cindrt account cost"""

ccont: float = None
"""ccont account cost"""

c226: float = None

c2261: float = None

c2262: float = None

c2263: float = None

c227: float = None

c2271: float = None

c2272: float = None

c2273: float = None

c2274: float = None

c24: float = None

c241: float = None

c242: float = None

c243: float = None

c244: float = None

c245: float = None

c21: float = None

c211: float = None

c212: float = None

c213: float = None

c214: float = None

c2141: float = None

c2142: float = None

c215: float = None

c216: float = None

c217: float = None

c2171: float = None

c2172: float = None

c2173: float = None

c2174: float = None

c22: float = None

c2211: float = None

c2212: float = None

c22121: float = None

c22122: float = None

c22123: float = None

c22124: float = None

c22125: float = None

c22126: float = None

c22127: float = None

c2213: float = None

c22131: float = None

c22132: float = None

c2214: float = None

c2215: float = None

c2221: float = None

c22211: float = None

c22212: float = None

c22213: float = None

c22214: float = None

c22215: float = None

c2222: float = None

c22221: float = None

c22222: float = None

c22223: float = None

c22224: float = None

c2223: float = None

c223: float = None

c2231: float = None

c2232: float = None

c2233: float = None

c2234: float = None

c224: float = None

c2241: float = None

c2242: float = None

c2243: float = None

c2244: float = None

c2245: float = None

c2246: float = None

c225: float = None

c2251: float = None

c22511: float = None

c22512: float = None

c22513: float = None

c22514: float = None

c22515: float = None

c2252: float = None

c22521: float = None

c22522: float = None

c22523: float = None

c22524: float = None

c22525: float = None

c22526: float = None

c22527: float = None

c2253: float = None

chx: float = None

cpp: float = None

cppa: float = None

c22128: float = None

abktflnc: float = None
"""allowable first wall/blanket neutron fluence (MW-yr/m2) (`blktmodel=0`)"""


adivflnc: float = None
"""allowable divertor heat fluence (MW-yr/m2)"""


blkcst: float = None
"""blanket direct cost (M$)"""


c221: float = None
"""total account 221 cost (M$) - first wall, blanket, shield, support structure and div plates"""


c222: float = None
"""total account 222 cost (M$) - TF coils + PF coils"""


capcost: float = None
"""total capital cost including interest (M$)"""


cconfix: float = None
"""fixed cost of superconducting cable ($/m)"""


cconshpf: float = None
"""cost of PF coil steel conduit/sheath ($/m)"""


cconshtf: float = None
"""cost of TF coil steel conduit/sheath ($/m)"""


cdcost: float = None
"""current drive direct costs (M$)"""


cdirt: float = None
"""total plant direct cost (M$)"""


cdrlife: float = None
"""Full power year lifetime of heating/current drive system (y)"""


cdrlife_cal: float = None
"""Calendar year lifetime of heating/current drive system (y)"""


cfactr: float = None
"""Total plant availability fraction; input if `i_plant_availability=0`"""


cpfact: float = None
"""Total plant capacity factor"""


cfind: list[float] = None
"""indirect cost factor (func of lsa) (cost model = 0)"""


cland: float = None
"""cost of land (M$)"""


coe: float = None
"""cost of electricity ($/MW-hr)"""


coecap: float = None
"""capital cost of electricity (m$/kW-hr)"""


coefuelt: float = None
"""'fuel' (including replaceable components) contribution to cost of electricity (m$/kW-hr)"""


coeoam: float = None
"""operation and maintenance contribution to cost of electricity (m$/kW-hr)"""


concost: float = None
"""plant construction cost (M$)"""


costexp: float = None
"""cost exponent for scaling in 2015 costs model"""


costexp_pebbles: float = None
"""cost exponent for pebbles in 2015 costs model"""


cost_factor_buildings: float = None
"""cost scaling factor for buildings"""


cost_factor_land: float = None
"""cost scaling factor for land"""


cost_factor_tf_coils: float = None
"""cost scaling factor for TF coils"""


cost_factor_tf_coils: float = None
"""cost scaling factor for fwbs"""


cost_factor_rh: float = None
"""cost scaling factor for remote handling"""


cost_factor_vv: float = None
"""cost scaling factor for vacuum vessel"""


cost_factor_bop: float = None
"""cost scaling factor for energy conversion system"""


cost_factor_misc: float = None
"""cost scaling factor for remaining subsystems"""


maintenance_fwbs: float = None
"""Maintenance cost factor: first wall, blanket, shield, divertor"""


maintenance_gen: float = None
"""Maintenance cost factor: All other components except coils, vacuum vessel,
thermal shield, cryostat, land"""


amortization: float = None
"""amortization factor (fixed charge factor) "A" (years)"""


cost_model: int = None
"""Switch for cost model:
- =0 use $ 1990 PROCESS model
- =1 use $ 2014 Kovari model
- =2 use user-provided model
"""


i_cp_lifetime: int = None
"""Switch for the centrepost lifetime constraint
0 : The CP full power year lifetime is set by the user via cplife_input
1 : The CP lifetime is equal to the divertor lifetime
2 : The CP lifetime is equal to the breeding blankets lifetime
3 : The CP lifetime is equal to the plant lifetime
"""


cowner: float = None
"""owner cost factor"""


cplife_input: float = None
"""User input full power year lifetime of the centrepost (years) (i_cp_lifetime = 0)"""


cplife: float = None
"""Calculated full power year lifetime of centrepost (years)"""


cplife_cal: float = None
"""Calculated calendar year lifetime of centrepost (years)"""


cpstcst: float = None
"""ST centrepost direct cost (M$)"""


cpstflnc: float = None
"""allowable ST centrepost neutron fluence (MW-yr/m2)"""


crctcore: float = None
"""reactor core costs (categories 221, 222 and 223)"""


csi: float = None
"""allowance for site costs (M$)"""


cturbb: float = None
"""cost of turbine building (M$)"""


decomf: float = None
"""proportion of constructed cost required for decommissioning fund"""


dintrt: float = None
"""diff between borrowing and saving interest rates"""


divcst: float = None
"""divertor direct cost (M$)"""


divlife: float = None
"""Full power lifetime of divertor (y)"""


divlife_cal: float = None
"""Calendar year lifetime of divertor (y)"""


dtlife: float = None
"""period prior to the end of the plant life that the decommissioning fund is used (years)"""


fcap0: float = None
"""average cost of money for construction of plant assuming design/construction time of six years"""


fcap0cp: float = None
"""average cost of money for replaceable components assuming lead time for these of two years"""


fcdfuel: float = None
"""fraction of current drive cost treated as fuel (if `ifueltyp = 1`)"""


fcontng: float = None
"""project contingency factor"""


fcr0: float = None
"""fixed charge rate during construction"""


fkind: float = None
"""multiplier for Nth of a kind costs"""


fwallcst: float = None
"""first wall cost (M$)"""


i_plant_availability: int = None
"""Switch for plant availability model:
- =0 use input value for cfactr
- =1 calculate cfactr using Taylor and Ward 1999 model
- =2 calculate cfactr using new (2015) model
- =3 calculate cfactr using ST model
"""


ibkt_life: int = None
"""Switch for fw/blanket lifetime calculation in availability module:
- =0 use neutron fluence model
- =1 use fusion power model (DEMO only)"""


life_dpa: float = None
"""Allowable DPA from DEMO fw/blanket lifetime calculation in availability module"""


bktcycles: float = None
"""Number of fusion cycles to reach allowable DPA from DEMO fw/blanket lifetime calculation"""


avail_min: float = None
"""Minimum availability (`constraint equation 61`)"""


tok_build_cost_per_vol: float = None
"""Unit cost for tokamak complex buildings, including building and site services ($/m3)"""


light_build_cost_per_vol: float = None
"""Unit cost for unshielded non-active buildings ($/m3)"""


favail: float = None
"""F-value for minimum availability (`constraint equation 61`)"""


num_rh_systems: int = None
"""Number of remote handling systems (1-10)"""


conf_mag: float = None
"""c parameter, which determines the temperature margin at which magnet lifetime starts to decline"""


div_prob_fail: float = None
"""Divertor probability of failure (per op day)"""


div_umain_time: float = None
"""Divertor unplanned maintenance time (years)"""


div_nref: float = None
"""Reference value for cycle cycle life of divertor"""


div_nu: float = None
"""The cycle when the divertor fails with 100% probability"""


fwbs_nref: float = None
"""Reference value for cycle life of blanket"""


fwbs_nu: float = None
"""The cycle when the blanket fails with 100% probability"""


fwbs_prob_fail: float = None
"""Fwbs probability of failure (per op day)"""


fwbs_umain_time: float = None
"""Fwbs unplanned maintenance time (years)"""


redun_vacp: float = None
"""Vacuum system pump redundancy level (%)"""


redun_vac: int = None
"""Number of redundant vacuum pumps"""


t_plant_operational_total_yrs: float = None
"""Operational time (yrs)"""


t_blkt_replace_yrs: float = None
"""time taken to replace blanket (y) (`i_plant_availability=1`)"""


tcomrepl: float = None
"""time taken to replace both blanket and divertor (y) (`i_plant_availability=1`)"""


t_div_replace_yrs: float = None
"""Time taken to replace divertor (y) (`i_plant_availability=1`)"""


uubop: float = None
"""unplanned unavailability factor for balance of plant (`i_plant_availability=1`)"""


uucd: float = None
"""unplanned unavailability factor for current drive (`i_plant_availability=1`)"""


uudiv: float = None
"""unplanned unavailability factor for divertor (`i_plant_availability=1`)"""


uufuel: float = None
"""unplanned unavailability factor for fuel system (`i_plant_availability=1`)"""


uufw: float = None
"""unplanned unavailability factor for first wall (`i_plant_availability=1`)"""


uumag: float = None
"""unplanned unavailability factor for magnets (`i_plant_availability=1`)"""


uuves: float = None
"""unplanned unavailability factor for vessel (`i_plant_availability=1`)"""


ifueltyp: int = None
"""Switch for fuel type:
- =2 treat initial blanket, divertor, first wall
as capital costs. Treat all later items and
fraction fcdfuel of CD equipment as fuel costs
- =1 treat blanket divertor, first wall and
fraction fcdfuel of CD equipment as fuel cost
- =0 treat these as capital cost
"""


ipnet: int = None
"""Switch for net electric power calculation:
- =0 scale so that always > 0
- =1 let go < 0 (no c-o-e)
"""


ireactor: int = None
"""Switch for net electric power and cost of electricity calculations:
- =0 do not calculate MW(electric) or c-o-e
- =1 calculate MW(electric) and c-o-e
"""


lsa: int = None
"""Level of safety assurance switch (generally, use 3 or 4):
- =1 truly passively safe plant
- =2,3 in-between
- =4 like current fission plant
"""


moneyint: float = None
"""interest portion of capital cost (M$)"""


output_costs: int = None
"""Switch for costs output:
- =0 do not write cost-related outputs to file
- =1 write cost-related outputs to file
"""


discount_rate: float = None
"""effective cost of money in constant dollars"""


startupratio: float = None
"""ratio of additional HCD power for start-up to flat-top operational requirements"""


startuppwr: float = None
"""cost associated with additional HCD system power required on start-up ($)"""


supercond_cost_model: int = None
"""Switch for superconductor cost model:
- =0 use $/kg
- =1 use $/kAm"""


tlife: float = None
"""Full power year plant lifetime (years)"""


tmain: float = None
"""Maintenance time for replacing CP (years) (i_plant_availability = 3)"""


u_unplanned_cp: float = None
"""User-input CP unplanned unavailability (i_plant_availability = 3)"""


UCAD: float = 180.0
"""unit cost for administration buildings (M$/m3)"""


UCAF: float = 1.5e6
"""unit cost for aux facility power equipment ($)"""


UCAHTS: float = 31.0
"""unit cost for aux heat transport equipment ($/W**exphts)"""


UCAP: float = 17.0
"""unit cost of auxiliary transformer ($/kVA)"""


ucblbe: float = None
"""unit cost for blanket beryllium ($/kg)"""


ucblbreed: float = None
"""unit cost for breeder material ($/kg) (`blktmodel>0`)"""


ucblli: float = None
"""unit cost for blanket lithium ($/kg) (30% Li6)"""


ucblli2o: float = None
"""unit cost for blanket Li_2O ($/kg)"""


ucbllipb: float = None
"""unit cost for blanket Li-Pb ($/kg) (30% Li6)"""


ucblss: float = None
"""unit cost for blanket stainless steel ($/kg)"""


ucblvd: float = None
"""unit cost for blanket vanadium ($/kg)"""


UCBPMP: float = 2.925e5
"""vacuum system backing pump cost ($)"""


ucbus: float = None
"""cost of aluminium bus for TF coil ($/A-m)"""


uccase: float = None
"""cost of superconductor case ($/kg)"""


UCCO: float = 350.0

"""unit cost for control buildings (M$/m3)"""


uccpcl1: float = None
"""cost of high strength tapered copper ($/kg)"""


uccpclb: float = None
"""cost of TF outboard leg plate coils ($/kg)"""


UCCPMP: float = 3.9e5
"""vacuum system cryopump cost ($)"""


UCCR: float = 460.0
"""unit cost for cryogenic building (M$/vol)"""


uccry: float = None
"""heat transport system cryoplant costs ($/W**expcry)"""


uccryo: float = None
"""unit cost for vacuum vessel ($/kg)"""


uccu: float = None
"""unit cost for copper in superconducting cable ($/kg)"""


UCDGEN: float = 1.7e6
"""cost per 8 MW diesel generator ($)"""


ucdiv: float = None
"""cost of divertor blade ($)"""


UCDTC: float = 13.0
"""detritiation, air cleanup cost ($/10000m3/hr)"""


UCDUCT: float = 4.225e4
"""vacuum system duct cost ($/m)"""


ucech: float = None
"""ECH system cost ($/W)"""


UCEL: float = 380.0
"""unit cost for electrical equipment building (M$/m3)"""


UCES1: float = 3.2e4
"""MGF (motor-generator flywheel) cost factor ($/MVA**0.8)"""


UCES2: float = 8.8e3
"""MGF (motor-generator flywheel) cost factor ($/MJ**0.8)"""


ucf1: float = None
"""cost of fuelling system ($)"""


ucfnc: float = None
"""outer PF coil fence support cost ($/kg)"""


UCFPR: float = 4.4e7

"""cost of 60g/day tritium processing unit ($)"""


ucfuel: float = None
"""unit cost of D-T fuel (M$/year/1200MW)"""


UCFWA: float = 6.0e4
"""first wall armour cost ($/m2)"""


UCFWPS: float = 1.0e7
"""first wall passive stabiliser cost ($)"""


UCFWS: float = 5.3e4
"""first wall structure cost ($/m2)"""


UCGSS: float = 35.0
"""cost of reactor structure ($/kg)"""


uche3: float = None
"""cost of helium-3 ($/kg)"""


uchrs: float = None
"""cost of heat rejection system ($)"""


uchts: list[float] = None
"""cost of heat transport system equipment per loop ($/W); dependent on coolant type (i_blkt_coolant_type)"""


uciac: float = None
"""cost of instrumentation, control & diagnostics ($)"""


ucich: float = None
"""ICH system cost ($/W)"""


UCINT: float = 35.0
"""superconductor intercoil structure cost ($/kg)"""


uclh: float = None
"""lower hybrid system cost ($/W)"""


UCLV: float = 16.0
"""low voltage system cost ($/kVA)"""


UCMB: float = 260.0
"""unit cost for reactor maintenance building (M$/m3)"""


ucme: float = None
"""cost of maintenance equipment ($)"""


ucmisc: float = None
"""miscellaneous plant allowance ($)"""


ucnbi: float = None
"""NBI system cost ($/W)"""


UCNBV: float = 1000.0
"""cost of nuclear building ventilation ($/m3)"""


ucoam: list[float] = None
"""annual cost of operation and maintenance (M$/year/1200MW**0.5)"""


ucpens: float = None
"""penetration shield cost ($/kg)"""


ucpfb: float = None
"""cost of PF coil buses ($/kA-m)"""


ucpfbk: float = None
"""cost of PF coil DC breakers ($/MVA**0.7)"""


ucpfbs: float = None
"""cost of PF burn power supplies ($/kW**0.7)"""


ucpfcb: float = None
"""cost of PF coil AC breakers ($/circuit)"""


ucpfdr1: float = None
"""cost factor for dump resistors ($/MJ)"""


ucpfic: float = None
"""cost of PF instrumentation and control ($/channel)"""


ucpfps: float = None
"""cost of PF coil pulsed power supplies ($/MVA)"""


UCPHX: float = 15.0
"""primary heat transport cost ($/W**exphts)"""


UCPP: float = 48.0
"""cost of primary power transformers ($/kVA**0.9)"""


ucrb: float = None
"""cost of reactor building (M$/m3)"""


ucsc: list[float] = None
"""cost of superconductor ($/kg)"""


sc_mat_cost_0: list[float] = None
"""cost of superconductor ($/kA m) at 6.4 T, 4.2 K"""


UCSH: float = 115.0
"""cost of shops and warehouses (M$/m3)"""


ucshld: float = None
"""cost of shield structural steel ($/kg)"""


UCSWYD: float = 1.84e7
"""switchyard equipment costs ($)"""


uctfbr: float = None
"""cost of TF coil breakers ($/W**0.7)"""


uctfbus: float = None
"""cost of TF coil bus ($/kg)"""


UCTFDR: float = 1.75e-4
"""cost of TF coil dump resistors ($/J)"""


UCTFGR: float = 5000.0
"""additional cost of TF coil dump resistors ($/coil)"""


UCTFIC: float = 1.0e4
"""cost of TF coil instrumentation and control ($/coil/30)"""


uctfps: float = None
"""cost of TF coil power supplies ($/W**0.7)"""


uctfsw: float = None
"""cost of TF coil slow dump switches ($/A)"""


UCTPMP: float = 1.105e5
"""cost of turbomolecular pump ($)"""


UCTR: float = 370.0
"""cost of tritium building ($/m3)"""


ucturb: list[float] = None
"""cost of turbine plant equipment ($) (dependent on coolant type i_blkt_coolant_type)"""


UCVALV: float = 3.9e5
"""vacuum system valve cost ($)"""


UCVDSH: float = 26.0
"""vacuum duct shield cost ($/kg)"""


UCVIAC: float = 1.3e6
"""vacuum system instrumentation and control cost ($)"""


ucwindpf: float = None
"""cost of PF coil superconductor windings ($/m)"""


ucwindtf: float = None
"""cost of TF coil superconductor windings ($/m)"""


UCWS: float = 460.0
"""cost of active assembly shop ($/m3)"""


ucwst: list[float] = None
"""cost of waste disposal (M$/y/1200MW)"""


def init_cost_variables():
    global c228
    c228 = 0.0

    global c229
    c229 = 0.0

    global c23
    c23 = 0.0

    global c25
    c25 = 0.0

    global c26
    c26 = 0.0

    global cindrt
    cindrt = 0.0

    global ccont
    ccont = 0.0

    global c226
    c226 = 0.0

    global c2261
    c2261 = 0.0

    global c2262
    c2262 = 0.0

    global c2263
    c2263 = 0.0

    global c227
    c227 = 0.0

    global c2271
    c2271 = 0.0

    global c2272
    c2272 = 0.0

    global c2273
    c2273 = 0.0

    global c2274
    c2274 = 0.0

    global c24
    c24 = 0.0

    global c241
    c241 = 0.0

    global c242
    c242 = 0.0

    global c243
    c243 = 0.0

    global c244
    c244 = 0.0

    global c245
    c245 = 0.0

    global c21
    c21 = 0.0

    global c211
    c211 = 0.0

    global c212
    c212 = 0.0

    global c213
    c213 = 0.0

    global c214
    c214 = 0.0

    global c2141
    c2141 = 0.0

    global c2142
    c2142 = 0.0

    global c215
    c215 = 0.0

    global c216
    c216 = 0.0

    global c217
    c217 = 0.0

    global c2171
    c2171 = 0.0

    global c2172
    c2172 = 0.0

    global c2173
    c2173 = 0.0

    global c2174
    c2174 = 0.0

    global c22
    c22 = 0.0

    global c2211
    c2211 = 0.0

    global c2212
    c2212 = 0.0

    global c22121
    c22121 = 0.0

    global c22122
    c22122 = 0.0

    global c22123
    c22123 = 0.0

    global c22124
    c22124 = 0.0

    global c22125
    c22125 = 0.0

    global c22126
    c22126 = 0.0

    global c22127
    c22127 = 0.0

    global c2213
    c2213 = 0.0

    global c22131
    c22131 = 0.0

    global c22132
    c22132 = 0.0

    global c2214
    c2214 = 0.0

    global c2215
    c2215 = 0.0

    global c2221
    c2221 = 0.0

    global c22211
    c22211 = 0.0

    global c22212
    c22212 = 0.0

    global c22213
    c22213 = 0.0

    global c22214
    c22214 = 0.0

    global c22215
    c22215 = 0.0

    global c2222
    c2222 = 0.0

    global c22221
    c22221 = 0.0

    global c22222
    c22222 = 0.0

    global c22223
    c22223 = 0.0

    global c22224
    c22224 = 0.0

    global c2223
    c2223 = 0.0

    global c223
    c223 = 0.0

    global c2231
    c2231 = 0.0

    global c2232
    c2232 = 0.0

    global c2233
    c2233 = 0.0

    global c2234
    c2234 = 0.0

    global c224
    c224 = 0.0

    global c2241
    c2241 = 0.0

    global c2242
    c2242 = 0.0

    global c2243
    c2243 = 0.0

    global c2244
    c2244 = 0.0

    global c2245
    c2245 = 0.0

    global c2246
    c2246 = 0.0

    global c225
    c225 = 0.0

    global c2251
    c2251 = 0.0

    global c22511
    c22511 = 0.0

    global c22512
    c22512 = 0.0

    global c22513
    c22513 = 0.0

    global c22514
    c22514 = 0.0

    global c22515
    c22515 = 0.0

    global c2252
    c2252 = 0.0

    global c22521
    c22521 = 0.0

    global c22522
    c22522 = 0.0

    global c22523
    c22523 = 0.0

    global c22524
    c22524 = 0.0

    global c22525
    c22525 = 0.0

    global c22526
    c22526 = 0.0

    global c22527
    c22527 = 0.0

    global c2253
    c2253 = 0.0

    global chx
    chx = 0.0

    global cpp
    cpp = 0.0

    global cppa
    cppa = 0.0

    global c22128
    c22128 = 0.0

    global abktflnc
    global adivflnc
    global blkcst
    global c221
    global c222
    global capcost
    global cconfix
    global cconshpf
    global cconshtf
    global cdcost
    global cdirt
    global cdrlife
    global cdrlife_cal
    global cfactr
    global cpfact
    global cfind
    global cland
    global coe
    global coecap
    global coefuelt
    global coeoam
    global concost
    global costexp
    global costexp_pebbles
    global cost_factor_buildings
    global cost_factor_land
    global cost_factor_tf_coils
    global cost_factor_fwbs
    global cost_factor_rh
    global cost_factor_vv
    global cost_factor_bop
    global cost_factor_misc
    global maintenance_fwbs
    global maintenance_gen
    global amortization
    global cost_model
    global i_cp_lifetime
    global cowner
    global cplife_input
    global cplife
    global cplife_cal
    global cpstcst
    global cpstflnc
    global crctcore
    global csi
    global cturbb
    global decomf
    global dintrt
    global divcst
    global divlife
    global divlife_cal
    global dtlife
    global fcap0
    global fcap0cp
    global fcdfuel
    global fcontng
    global fcr0
    global fkind
    global fwallcst
    global i_plant_availability
    global ibkt_life
    global life_dpa
    global bktcycles
    global avail_min
    global tok_build_cost_per_vol
    global light_build_cost_per_vol
    global favail
    global num_rh_systems
    global conf_mag
    global div_prob_fail
    global div_umain_time
    global div_nref
    global div_nu
    global fwbs_nref
    global fwbs_nu
    global fwbs_prob_fail
    global fwbs_umain_time
    global redun_vacp
    global redun_vac
    global t_plant_operational_total_yrs
    global t_blkt_replace_yrs
    global tcomrepl
    global t_div_replace_yrs
    global uubop
    global uucd
    global uudiv
    global uufuel
    global uufw
    global uumag
    global uuves
    global ifueltyp
    global ipnet
    global ireactor
    global lsa
    global moneyint
    global output_costs
    global discount_rate
    global startupratio
    global startuppwr
    global supercond_cost_model
    global tlife
    global tmain
    global u_unplanned_cp
    global ucblbe
    global ucblbreed
    global ucblli
    global ucblli2o
    global ucbllipb
    global ucblss
    global ucblvd
    global ucbus
    global uccase
    global uccpcl1
    global uccpclb
    global uccry
    global uccryo
    global uccu
    global ucdiv
    global ucech
    global uces1
    global uces2
    global ucf1
    global ucfnc
    global ucfuel
    global uche3
    global uchrs
    global uchts
    global uciac
    global ucich
    global uclh
    global ucme
    global ucmisc
    global ucnbi
    global ucoam
    global ucpens
    global ucpfb
    global ucpfbk
    global ucpfbs
    global ucpfcb
    global ucpfdr1
    global ucpfps
    global ucrb
    global ucsc
    global sc_mat_cost_0
    global ucshld
    global uctfbr
    global uctfbus
    global uctfps
    global uctfsw
    global uctpmp
    global ucturb
    global ucwindpf
    global ucwindtf
    global ucws
    global ucwst
    global ucpfic

    abktflnc = 5.0
    adivflnc = 7.0
    blkcst = 0.0
    c221 = 0.0
    c222 = 0.0
    capcost = 0.0
    cconfix = 80.0
    cconshpf = 70.0
    cconshtf = 75.0
    cdcost = 0.0
    cdirt = 0.0
    cdrlife = 0.0
    cdrlife_cal = 0.0
    cfactr = 0.75
    cpfact = 0.0
    cfind = [0.244, 0.244, 0.244, 0.29]
    cland = 19.2
    coe = 0.0
    coecap = 0.0
    coefuelt = 0.0
    coeoam = 0.0
    concost = 0.0
    costexp = 0.8
    costexp_pebbles = 0.6
    cost_factor_buildings = 1.0
    cost_factor_land = 1.0
    cost_factor_tf_coils = 1.0
    cost_factor_fwbs = 1.0
    cost_factor_rh = 1.0
    cost_factor_vv = 1.0
    cost_factor_bop = 1.0
    cost_factor_misc = 1.0
    maintenance_fwbs = 0.2
    maintenance_gen = 0.05
    amortization = 13.6
    cost_model = 1
    cowner = 0.15
    cplife = 0.0
    cplife_cal = 0.0
    cpstcst = 0.0
    cpstflnc = 10.0
    crctcore = 0.0
    csi = 16.0
    cturbb = 38.0
    decomf = 0.1
    dintrt = 0.0
    divcst = 0.0
    divlife = 0.0
    divlife_cal = 0.0
    dtlife = 0.0
    fcap0 = 1.165
    fcap0cp = 1.08
    fcdfuel = 0.1
    fcontng = 0.195
    fcr0 = 0.0966
    fkind = 1.0
    fwallcst = 0.0
    i_plant_availability = 2
    ibkt_life = 0
    life_dpa = 50
    bktcycles = 1.0e3
    avail_min = 0.75
    tok_build_cost_per_vol = 1283.0
    light_build_cost_per_vol = 270.0
    favail = 1.0
    num_rh_systems = 4
    conf_mag = 0.99
    div_prob_fail = 0.0002
    div_umain_time = 0.25
    div_nref = 7000.0
    div_nu = 14000.0
    fwbs_nref = 20000.0
    fwbs_nu = 40000.0
    fwbs_prob_fail = 0.0002
    fwbs_umain_time = 0.25
    redun_vacp = 25.0
    redun_vac = 0
    t_plant_operational_total_yrs = 0.0
    t_blkt_replace_yrs = 0.5
    tcomrepl = 0.5
    t_div_replace_yrs = 0.25
    uubop = 0.02
    uucd = 0.02
    uudiv = 0.04
    uufuel = 0.02
    uufw = 0.04
    uumag = 0.02
    uuves = 0.04
    ifueltyp = 0
    ipnet = 0
    ireactor = 1
    lsa = 4
    moneyint = 0.0
    output_costs = 1
    discount_rate = 0.0435
    startupratio = 1.0
    startuppwr = 0.0
    tlife = 30.0
    ucblbe = 260.0
    ucblbreed = 875.0
    ucblli = 875.0
    ucblli2o = 600.0
    ucbllipb = 10.3
    ucblss = 90.0
    ucblvd = 200.0
    ucbus = 0.123
    uccase = 50.0
    uccpcl1 = 250.0
    uccpclb = 150.0
    uccry = 9.3e4
    uccryo = 32.0
    uccu = 75.0
    ucdiv = 2.8e5
    ucech = 3.0
    ucf1 = 2.23e7
    ucfnc = 35.0
    ucfuel = 3.45
    uche3 = 1.0e6
    uchrs = 87.9e6
    uchts = [15.3, 19.1]
    uciac = 1.5e8
    ucich = 3.0
    uclh = 3.3
    ucme = 1.25e8
    ucmisc = 2.5e7
    ucnbi = 3.3
    ucoam = [68.8, 68.8, 68.8, 74.4]
    ucpens = 32.0
    ucpfb = 210.0
    ucpfbk = 1.66e4
    ucpfbs = 4.9e3
    ucpfcb = 7.5e4
    ucpfdr1 = 150.0
    ucpfic = 1.0e4
    ucpfps = 3.5e4
    ucrb = 400.0
    ucsc = [
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
    sc_mat_cost_0 = [4.8, 2.0, 1.0, 4.8, 4.8, 47.4, 1.0, 47.4, 47.4]
    supercond_cost_model = 0
    ucshld = 32.0
    uctfbr = 1.22
    uctfbus = 100.0
    uctfps = 24.0
    uctfsw = 1.0
    ucturb = [230.0e6, 245.0e6]
    ucwindpf = 465.0
    ucwindtf = 480.0
    ucwst = [0.0, 3.94, 5.91, 7.88]
    u_unplanned_cp = 0.0
    i_cp_lifetime = 0
    cplife_input = 2.0
