module cost_variables
  !! author: J. Morris, S. Muldrew, M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to the costing algorithms of a fusion power plant.
  !!
  !!### References
  !!
  !! - AEA FUS 251: A User's Guide to the PROCESS Systems Code

#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: abktflnc
  !! allowable first wall/blanket neutron fluence (MW-yr/m2) (`blktmodel=0`)

  real(dp) :: adivflnc
  !! allowable divertor heat fluence (MW-yr/m2)

  real(dp) :: blkcst
  !! blanket direct cost (M$)

  real(dp) :: c221
  !! total account 221 cost (M$) - first wall, blanket, shield, support structure and div plates

  real(dp) :: c222
  !! total account 222 cost (M$) - TF coils + PF coils

  real(dp) :: capcost
  !! total capital cost including interest (M$)

  real(dp) :: cconfix
  !! fixed cost of superconducting cable ($/m)

  real(dp) :: cconshpf
  !! cost of PF coil steel conduit/sheath ($/m)

  real(dp) :: cconshtf
  !! cost of TF coil steel conduit/sheath ($/m)

  real(dp) :: cdcost
  !! current drive direct costs (M$)

  real(dp) :: cdirt
  !! total plant direct cost (M$)

  real(dp) :: cdrlife
  !! lifetime of heating/current drive system (y)

  real(dp) :: cfactr
  !! Total plant availability fraction; input if `iavail=0`

  real(dp) :: cpfact
  !! Total plant capacity factor

  real(dp), dimension(4) :: cfind
  !! indirect cost factor (func of lsa) (cost model = 0)

  real(dp) :: cland
  !! cost of land (M$)

  real(dp) :: coe
  !! cost of electricity ($/MW-hr)

  real(dp) :: coecap
  !! capital cost of electricity (m$/kW-hr)

  real(dp) :: coefuelt
  !! 'fuel' (including replaceable components) contribution to cost of electricity (m$/kW-hr)

  real(dp) :: coeoam
  !! operation and maintenance contribution to cost of electricity (m$/kW-hr)

  real(dp) :: concost
  !! plant construction cost (M$)

  real(dp) :: costexp
  !! cost exponent for scaling in 2015 costs model

  real(dp) :: costexp_pebbles
  !! cost exponent for pebbles in 2015 costs model

  real(dp) :: cost_factor_buildings
  !! cost scaling factor for buildings

  real(dp) :: cost_factor_land
  !! cost scaling factor for land

  real(dp) :: cost_factor_tf_coils
  !! cost scaling factor for TF coils

  real(dp) :: cost_factor_fwbs
  !! cost scaling factor for fwbs

  real(dp) :: cost_factor_rh
  !! cost scaling factor for remote handling

  real(dp) :: cost_factor_vv
  !! cost scaling factor for vacuum vessel

  real(dp) :: cost_factor_bop
  !! cost scaling factor for energy conversion system

  real(dp) :: cost_factor_misc
  !! cost scaling factor for remaining subsystems

  real(dp) :: maintenance_fwbs
  !! Maintenance cost factor: first wall, blanket, shield, divertor

  real(dp) :: maintenance_gen
  !! Maintenance cost factor: All other components except coils, vacuum vessel,
  !! thermal shield, cryostat, land

  real(dp) :: amortization
  !! amortization factor (fixed charge factor) "A" (years)

  integer :: cost_model
  !! Switch for cost model:
  !!
  !! - =0 use $ 1990 PROCESS model
  !! - =1 use $ 2014 Kovari model
  !! - =2 use user-provided model

  integer :: i_cp_lifetime
  !! Switch for the centrepost lifetime constraint
  !!  0 : The CP full power year lifetime is set by the user via cplife_input
  !!  1 : The CP lifetime is equal to the divertor lifetime
  !!  2 : The CP lifetime is equal to the breeding blankets lifetime
  !!  3 : The CP lifetime is equal to the plant lifetime

  real(dp) :: cowner
  !! owner cost factor

  real(dp) :: cplife_input
  !! User input full power year lifetime of the centrepost (years) (i_cp_lifetime = 0)

  real(dp) :: cplife
  !! Calculated full power year lifetime of centrepost (years)

  real(dp) :: cpstcst
  !! ST centrepost direct cost (M$)

  real(dp) :: cpstflnc
  !! allowable ST centrepost neutron fluence (MW-yr/m2)

  real(dp) :: crctcore
  !! reactor core costs (categories 221, 222 and 223)

  real(dp) :: csi
  !! allowance for site costs (M$)

  real(dp) :: cturbb
  !! cost of turbine building (M$)

  real(dp) :: decomf
  !! proportion of constructed cost required for decommissioning fund

  real(dp) :: dintrt
  !! diff between borrowing and saving interest rates

  real(dp) :: divcst
  !! divertor direct cost (M$)

  real(dp) :: divlife
  !! Full power lifetime of divertor (y)

  real(dp) :: dtlife
  !! period prior to the end of the plant life that the decommissioning fund is used (years)

  real(dp) :: fcap0
  !! average cost of money for construction of plant assuming design/construction time of six years

  real(dp) :: fcap0cp
  !! average cost of money for replaceable components assuming lead time for these of two years

  real(dp) :: fcdfuel
  !! fraction of current drive cost treated as fuel (if `ifueltyp = 1`)

  real(dp) :: fcontng
  !! project contingency factor

  real(dp) :: fcr0
  !! fixed charge rate during construction

  real(dp) :: fkind
  !! multiplier for Nth of a kind costs

  real(dp) :: fwallcst
  !! first wall cost (M$)

  integer :: iavail
  !! Switch for plant availability model:
  !!
  !! - =0 use input value for cfactr
  !! - =1 calculate cfactr using Taylor and Ward 1999 model
  !! - =2 calculate cfactr using new (2015) model
  !! - =3 calculate cfactr using ST model

  integer :: ibkt_life
  !! Switch for fw/blanket lifetime calculation in availability module:
  !!
  !! - =0 use neutron fluence model
  !! - =1 use fusion power model (DEMO only)

  real(dp) :: life_dpa
  !! Allowable DPA from DEMO fw/blanket lifetime calculation in availability module

  real(dp) :: bktcycles
  !! Number of fusion cycles to reach allowable DPA from DEMO fw/blanket lifetime calculation

  real(dp) :: avail_min
  !! Minimum availability (`constraint equation 61`)

  real(dp) :: tok_build_cost_per_vol
  !! Unit cost for tokamak complex buildings, including building and site services ($/m3)

  real(dp) :: light_build_cost_per_vol
  !! Unit cost for unshielded non-active buildings ($/m3)

  real(dp) :: favail
  !! F-value for minimum availability (`constraint equation 61`)

  integer :: num_rh_systems
  !! Number of remote handling systems (1-10)

  real(dp) :: conf_mag
  !! c parameter, which determines the temperature margin at which magnet lifetime starts to decline

  real(dp) :: div_prob_fail
  !! Divertor probability of failure (per op day)

  real(dp) :: div_umain_time
  !! Divertor unplanned maintenance time (years)

  real(dp) :: div_nref
  !! Reference value for cycle cycle life of divertor

  real(dp) :: div_nu
  !! The cycle when the divertor fails with 100% probability

  real(dp) :: fwbs_nref
  !! Reference value for cycle life of blanket

  real(dp) :: fwbs_nu
  !! The cycle when the blanket fails with 100% probability

  real(dp) :: fwbs_prob_fail
  !! Fwbs probability of failure (per op day)

  real(dp) :: fwbs_umain_time
  !! Fwbs unplanned maintenance time (years)

  real(dp) :: redun_vacp
  !! Vacuum system pump redundancy level (%)

  integer :: redun_vac
  !! Number of redundant vacuum pumps

  real(dp) :: t_operation
  !! Operational time (yrs)

  real(dp) :: tbktrepl
  !! time taken to replace blanket (y) (`iavail=1`)

  real(dp) :: tcomrepl
  !! time taken to replace both blanket and divertor (y) (`iavail=1`)

  real(dp) :: tdivrepl
  !! time taken to replace divertor (y) (`iavail=1`)

  real(dp) :: uubop
  !! unplanned unavailability factor for balance of plant (`iavail=1`)

  real(dp) :: uucd
  !! unplanned unavailability factor for current drive (`iavail=1`)

  real(dp) :: uudiv
  !! unplanned unavailability factor for divertor (`iavail=1`)

  real(dp) :: uufuel
  !! unplanned unavailability factor for fuel system (`iavail=1`)

  real(dp) :: uufw
  !! unplanned unavailability factor for first wall (`iavail=1`)

  real(dp) :: uumag
  !! unplanned unavailability factor for magnets (`iavail=1`)

  real(dp) :: uuves
  !! unplanned unavailability factor for vessel (`iavail=1`)

  integer :: ifueltyp
  !! Switch for fuel type:
  !!
  !! - =2 treat initial blanket, divertor, first wall
  !!   as capital costs. Treat all later items and
  !!   fraction fcdfuel of CD equipment as fuel costs
  !! - =1 treat blanket divertor, first wall and
  !!   fraction fcdfuel of CD equipment as fuel cost
  !! - =0 treat these as capital cost

  integer :: ipnet
  !! Switch for net electric power calculation:
  !!
  !! - =0 scale so that always > 0
  !! - =1 let go < 0 (no c-o-e)

  integer :: ireactor
  !! Switch for net electric power and cost of electricity calculations:
  !!
  !! - =0 do not calculate MW(electric) or c-o-e
  !! - =1 calculate MW(electric) and c-o-e

  integer :: lsa
  !! Level of safety assurance switch (generally, use 3 or 4):
  !!
  !! - =1 truly passively safe plant
  !! - =2,3 in-between
  !! - =4 like current fission plant

  real(dp) :: moneyint
  !! interest portion of capital cost (M$)

  integer :: output_costs
  !! Switch for costs output:
  !!
  !! - =0 do not write cost-related outputs to file
  !! - =1 write cost-related outputs to file

  real(dp) :: discount_rate
  !! effective cost of money in constant dollars

  real(dp) :: startupratio
  !! ratio of additional HCD power for start-up to flat-top operational requirements

  real(dp) :: startuppwr
  !! cost associated with additional HCD system power required on start-up ($)

  real(dp) :: tlife
  !! Full power year plant lifetime (years)

  real(dp) :: tmain
  !! Maintenance time for replacing CP (years) (iavail = 3)

  real(dp) :: u_unplanned
  !! User-input CP unplanned unavailability (iavail = 3)

  real(dp), parameter :: ucad = 180.0D0
  !! unit cost for administration buildings (M$/m3)

  real(dp), parameter :: ucaf = 1.5D6
  !! unit cost for aux facility power equipment ($)

  real(dp), parameter :: ucahts = 31.0D0
  !! unit cost for aux heat transport equipment ($/W**exphts)

  real(dp), parameter :: ucap = 17.0D0
  !! unit cost of auxiliary transformer ($/kVA)

  real(dp) :: ucblbe
  !! unit cost for blanket beryllium ($/kg)

  real(dp) :: ucblbreed
  !! unit cost for breeder material ($/kg) (`blktmodel>0`)

  real(dp) :: ucblli
  !! unit cost for blanket lithium ($/kg) (30% Li6)

  real(dp) :: ucblli2o
  !! unit cost for blanket Li_2O ($/kg)

  real(dp) :: ucbllipb
  !! unit cost for blanket Li-Pb ($/kg) (30% Li6)

  real(dp) :: ucblss
  !! unit cost for blanket stainless steel ($/kg)

  real(dp) :: ucblvd
  !! unit cost for blanket vanadium ($/kg)

  real(dp), parameter :: ucbpmp = 2.925D5
  !! vacuum system backing pump cost ($)

  real(dp) :: ucbus
  !! cost of aluminium bus for TF coil ($/A-m)

  real(dp) :: uccase
  !! cost of superconductor case ($/kg)

  real(dp), parameter :: ucco = 350.0D0
  !! unit cost for control buildings (M$/m3)

  real(dp) :: uccpcl1
  !! cost of high strength tapered copper ($/kg)

  real(dp) :: uccpclb
  !! cost of TF outboard leg plate coils ($/kg)

  real(dp), parameter :: uccpmp = 3.9D5
  !! vacuum system cryopump cost ($)

  real(dp), parameter :: uccr = 460.0D0
  !! unit cost for cryogenic building (M$/vol)

  real(dp) :: uccry
  !! heat transport system cryoplant costs ($/W**expcry)

  real(dp) :: uccryo
  !! unit cost for vacuum vessel ($/kg)

  real(dp) :: uccu
  !! unit cost for copper in superconducting cable ($/kg)

  real(dp), parameter :: ucdgen = 1.7D6
  !! cost per 8 MW diesel generator ($)

  real(dp) :: ucdiv
  !! cost of divertor blade ($)

  real(dp), parameter :: ucdtc = 13.0D0
  !! detritiation, air cleanup cost ($/10000m3/hr)

  real(dp), parameter :: ucduct = 4.225D4
  !! vacuum system duct cost ($/m)

  real(dp) :: ucech
  !! ECH system cost ($/W)

  real(dp), parameter :: ucel = 380.0D0
  !! unit cost for electrical equipment building (M$/m3)

  real(dp), parameter :: uces1 = 3.2D4
  !! MGF (motor-generator flywheel) cost factor ($/MVA**0.8)

  real(dp), parameter :: uces2 = 8.8D3
  !! MGF (motor-generator flywheel) cost factor ($/MJ**0.8)

  real(dp) :: ucf1
  !! cost of fuelling system ($)

  real(dp) :: ucfnc
  !! outer PF coil fence support cost ($/kg)

  real(dp), parameter :: ucfpr = 4.4D7
  !! cost of 60g/day tritium processing unit ($)

  real(dp) :: ucfuel
  !! unit cost of D-T fuel (M$/year/1200MW)

  real(dp), parameter :: ucfwa = 6.0D4
  !! first wall armour cost ($/m2)

  real(dp), parameter :: ucfwps = 1.0D7
  !! first wall passive stabiliser cost ($)

  real(dp), parameter :: ucfws = 5.3D4
  !! first wall structure cost ($/m2)

  real(dp), parameter :: ucgss = 35.0D0
  !! cost of reactor structure ($/kg)

  real(dp) :: uche3
  !! cost of helium-3 ($/kg)

  real(dp) :: uchrs
  !! cost of heat rejection system ($)

  real(dp), dimension(2) :: uchts
  !! cost of heat transport system equipment per loop ($/W); dependent on coolant type (coolwh)

  real(dp) :: uciac
  !! cost of instrumentation, control & diagnostics ($)

  real(dp) :: ucich
  !! ICH system cost ($/W)

  real(dp), parameter :: ucint = 35.0D0
  !! superconductor intercoil structure cost ($/kg)

  real(dp) :: uclh
  !! lower hybrid system cost ($/W)

  real(dp), parameter :: uclv = 16.0D0
  !! low voltage system cost ($/kVA)

  real(dp), parameter :: ucmb = 260.0D0
  !! unit cost for reactor maintenance building (M$/m3)

  real(dp) :: ucme
  !! cost of maintenance equipment ($)

  real(dp) :: ucmisc
  !! miscellaneous plant allowance ($)

  real(dp) :: ucnbi
  !! NBI system cost ($/W)

  real(dp), parameter :: ucnbv = 1000.0D0
  !! cost of nuclear building ventilation ($/m3)

  real(dp), dimension(4) :: ucoam
  !! annual cost of operation and maintenance (M$/year/1200MW**0.5)

  real(dp) :: ucpens
  !! penetration shield cost ($/kg)

  real(dp) :: ucpfb
  !! cost of PF coil buses ($/kA-m)

  real(dp) :: ucpfbk
  !! cost of PF coil DC breakers ($/MVA**0.7)

  real(dp) :: ucpfbs
  !! cost of PF burn power supplies ($/kW**0.7)

  real(dp) :: ucpfcb
  !! cost of PF coil AC breakers ($/circuit)

  real(dp) :: ucpfdr1
  !! cost factor for dump resistors ($/MJ)

  real(dp) :: ucpfic
  !! cost of PF instrumentation and control ($/channel)

  real(dp) :: ucpfps
  !! cost of PF coil pulsed power supplies ($/MVA)

  real(dp), parameter :: ucphx = 15.0D0
  !! primary heat transport cost ($/W**exphts)

  real(dp), parameter :: ucpp = 48.0D0
  !! cost of primary power transformers ($/kVA**0.9)

  real(dp) :: ucrb
  !! cost of reactor building (M$/m3)

  real(dp), dimension(9) :: ucsc
  !! cost of superconductor ($/kg)

  real(dp), parameter :: ucsh = 115.0D0
  !! cost of shops and warehouses (M$/m3)

  real(dp) :: ucshld
  !! cost of shield structural steel ($/kg)

  real(dp), parameter :: ucswyd = 1.84D7
  !! switchyard equipment costs ($)

  real(dp) :: uctfbr
  !! cost of TF coil breakers ($/W**0.7)

  real(dp) :: uctfbus
  !! cost of TF coil bus ($/kg)

  real(dp), parameter :: uctfdr = 1.75D-4
  !! cost of TF coil dump resistors ($/J)

  real(dp), parameter :: uctfgr = 5000.0D0
  !! additional cost of TF coil dump resistors ($/coil)

  real(dp), parameter :: uctfic = 1.0D4
  !! cost of TF coil instrumentation and control ($/coil/30)

  real(dp) :: uctfps
  !! cost of TF coil power supplies ($/W**0.7)

  real(dp) :: uctfsw
  !! cost of TF coil slow dump switches ($/A)

  real(dp), parameter :: uctpmp = 1.105D5
  !! cost of turbomolecular pump ($)

  real(dp), parameter :: uctr = 370.0D0
  !! cost of tritium building ($/m3)

  real(dp), dimension(2) :: ucturb
  !! cost of turbine plant equipment ($) (dependent on coolant type coolwh)

  real(dp), parameter :: ucvalv = 3.9D5
  !! vacuum system valve cost ($)

  real(dp), parameter :: ucvdsh = 26.0D0
  !! vacuum duct shield cost ($/kg)

  real(dp), parameter :: ucviac = 1.3D6
  !! vacuum system instrumentation and control cost ($)

  real(dp) :: ucwindpf
  !! cost of PF coil superconductor windings ($/m)

  real(dp) :: ucwindtf
  !! cost of TF coil superconductor windings ($/m)

  real(dp), parameter :: ucws = 460.0D0
  !! cost of active assembly shop ($/m3)

  real(dp), dimension(4) :: ucwst
  !! cost of waste disposal (M$/y/1200MW)

  contains

  subroutine init_cost_variables
    !! Initialise cost variables
    implicit none

    abktflnc = 5.0D0
    adivflnc = 7.0D0
    blkcst = 0.0D0
    c221 = 0.0D0
    c222 = 0.0D0
    capcost = 0.0D0
    cconfix = 80.0D0
    cconshpf = 70.0D0
    cconshtf = 75.0D0
    cdcost = 0.0D0
    cdirt = 0.0D0
    cdrlife = 0.0D0
    cfactr = 0.75D0
    cpfact = 0.0D0
    cfind = (/0.244D0, 0.244D0, 0.244D0, 0.29D0/)
    cland = 19.2D0
    coe = 0.0D0
    coecap = 0.0D0
    coefuelt = 0.0D0
    coeoam = 0.0D0
    concost = 0.0D0
    costexp = 0.8D0
    costexp_pebbles = 0.6D0
    cost_factor_buildings = 1.0D0
    cost_factor_land = 1.0D0
    cost_factor_tf_coils = 1.0D0
    cost_factor_fwbs = 1.0D0
    cost_factor_rh = 1.0D0
    cost_factor_vv = 1.0D0
    cost_factor_bop = 1.0D0
    cost_factor_misc = 1.0D0
    maintenance_fwbs = 0.2D0
    maintenance_gen = 0.05D0
    amortization = 13.6D0
    cost_model = 1
    cowner = 0.15D0
    cplife = 0.0D0
    cpstcst = 0.0D0
    cpstflnc = 10.0D0
    crctcore = 0.0D0
    csi = 16.0D0
    cturbb = 38.0D0
    decomf = 0.1D0
    dintrt = 0.0D0
    divcst = 0.0D0
    divlife = 0.0D0
    dtlife = 0.0D0
    fcap0 = 1.165D0
    fcap0cp = 1.08D0
    fcdfuel = 0.1D0
    fcontng = 0.195D0
    fcr0 = 0.0966D0
    fkind = 1.0D0
    fwallcst = 0.0D0
    iavail= 2
    ibkt_life = 0
    life_dpa = 50
    bktcycles = 1.0D3
    avail_min = 0.75D0
    tok_build_cost_per_vol = 1283.0D0
    light_build_cost_per_vol = 270.0D0
    favail = 1.0D0
    num_rh_systems = 4
    conf_mag = 0.99D0
    div_prob_fail = 0.0002D0
    div_umain_time = 0.25D0
    div_nref = 7000.0D0
    div_nu = 14000.0D0
    fwbs_nref = 20000.0D0
    fwbs_nu = 40000.0D0
    fwbs_prob_fail = 0.0002D0
    fwbs_umain_time = 0.25D0
    redun_vacp = 25.0D0
    redun_vac = 0
    t_operation = 0.0D0
    tbktrepl = 0.5D0
    tcomrepl = 0.5D0
    tdivrepl = 0.25D0
    uubop = 0.02D0
    uucd = 0.02D0
    uudiv = 0.04D0
    uufuel = 0.02D0
    uufw = 0.04D0
    uumag = 0.02D0
    uuves = 0.04D0
    ifueltyp = 0
    ipnet = 0
    ireactor = 1
    lsa = 4
    moneyint = 0.0D0
    output_costs = 1
    discount_rate = 0.0435D0
    startupratio = 1.0
    startuppwr = 0.0
    tlife = 30.0D0
    ucblbe = 260.0D0
    ucblbreed = 875.0D0
    ucblli = 875.0D0
    ucblli2o = 600.0D0
    ucbllipb = 10.3D0
    ucblss = 90.0D0
    ucblvd = 200.0D0
    ucbus = 0.123D0
    uccase = 50.0D0
    uccpcl1 = 250.0D0
    uccpclb = 150.0D0
    uccry = 9.3D4
    uccryo = 32.0D0
    uccu = 75.0D0
    ucdiv = 2.8D5
    ucech = 3.0D0
    ucf1 = 2.23D7
    ucfnc = 35.0D0
    ucfuel = 3.45D0
    uche3 = 1.0D6
    uchrs = 87.9D6
    uchts = (/15.3D0, 19.1D0/)
    uciac = 1.5D8
    ucich = 3.0D0
    uclh = 3.3D0
    ucme = 1.25D8
    ucmisc = 2.5D7
    ucnbi = 3.3D0
    ucoam = (/68.8D0, 68.8D0, 68.8D0, 74.4D0/)
    ucpens = 32.0D0
    ucpfb = 210.0D0
    ucpfbk = 1.66D4
    ucpfbs = 4.9D3
    ucpfcb = 7.5D4
    ucpfdr1 = 150.0D0
    ucpfic = 1.0D4
    ucpfps = 3.5D4
    ucrb = 400.0D0
    ucsc = &
      (/600.0D0, 600.0D0, 300.0D0, 600.0D0, 600.0D0, 600.0D0, 300.0D0, 1200.0D0, &
      1200.0D0/)
    ucshld = 32.0D0
    uctfbr = 1.22D0
    uctfbus = 100.0D0
    uctfps = 24.0D0
    uctfsw = 1.0D0
    ucturb = (/230.0D6, 245.0D6/)
    ucwindpf = 465.0D0
    ucwindtf = 480.0D0
    ucwst = (/0.0D0, 3.94D0, 5.91D0, 7.88D0/)
    i_cp_lifetime = 0
    cplife_input = 2.0D0

  end subroutine init_cost_variables
end module cost_variables
