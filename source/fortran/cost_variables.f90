module cost_variables
  !! author: J. Morris, S. Muldrew, M. Kovari (UKAEA)
  !!
  !! Module containing global variables relating to the costing algorithms of a fusion power plant.
  !!
  !!### References
  !!
  !! -
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
  !! Full power year lifetime of heating/current drive system (y)

  real(dp) :: cdrlife_cal
  !! Calendar year lifetime of heating/current drive system (y)

  real(dp) :: cfactr
  !! Total plant availability fraction; input if `i_plant_availability=0`

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

  real(dp) :: cplife_cal
  !! Calculated calendar year lifetime of centrepost (years)

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

  real(dp) :: life_div_fpy
  !! Full power lifetime of divertor (fpy)

  real(dp) :: divlife_cal
  !! Calendar year lifetime of divertor (y)

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

  integer :: i_plant_availability
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

  real(dp) :: n_blkt_pulse_cycles
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

  real(dp) :: t_blkt_replace_years
  !! time taken to replace blanket (y) (`i_plant_availability=1`)

  real(dp) :: t_div_blkt_replace_years
  !! time taken to replace both blanket and divertor (y) (`i_plant_availability=1`)

  real(dp) :: t_div_replace_years
  !! time taken to replace divertor (y) (`i_plant_availability=1`)

  real(dp) :: uubop
  !! unplanned unavailability factor for balance of plant (`i_plant_availability=1`)

  real(dp) :: uucd
  !! unplanned unavailability factor for current drive (`i_plant_availability=1`)

  real(dp) :: uudiv
  !! unplanned unavailability factor for divertor (`i_plant_availability=1`)

  real(dp) :: uufuel
  !! unplanned unavailability factor for fuel system (`i_plant_availability=1`)

  real(dp) :: uufw
  !! unplanned unavailability factor for first wall (`i_plant_availability=1`)

  real(dp) :: uumag
  !! unplanned unavailability factor for magnets (`i_plant_availability=1`)

  real(dp) :: uuves
  !! unplanned unavailability factor for vessel (`i_plant_availability=1`)

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

  integer :: supercond_cost_model
  !! Switch for superconductor cost model:
  !!
  !! - =0 use $/kg
  !! - =1 use $/kAm

  real(dp) :: tlife
  !! Full power year plant lifetime (years)

  real(dp) :: tmain
  !! Maintenance time for replacing CP (years) (i_plant_availability = 3)

  real(dp) :: u_unplanned_cp
  !! User-input CP unplanned unavailability (i_plant_availability = 3)

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
  !! cost of heat transport system equipment per loop ($/W); dependent on coolant type (i_blkt_coolant_type)

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

  real(dp), dimension(9) :: sc_mat_cost_0
  !!cost of superconductor ($/kA m) at 6.4 T, 4.2 K

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
  !! cost of turbine plant equipment ($) (dependent on coolant type i_blkt_coolant_type)

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
end module cost_variables
