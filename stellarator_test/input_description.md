# Input for SQuID
Variables names update is needed

### Equaltities

- Global power balance   
`icc = 2`   
*(consistency equation, specific value not needed)*

- Net electric power lower limit [MW]  
*Necessary for cost optimization to avoid minimal working radius*  
`icc = 16`  
pnetelin = 1000  
*In our scan electric power was used as varying parameter*

### Inequaltities

#### Beta limits
- ~~Lower beta limit~~   
~~`icc = 84`~~   
~~beta_min = 0.01~~  
*In the old Helias input, there was a minimum values set for beta. In practice it doesn't change anything, as we want to maximize beta anyway. There is no physical limit restricting it from the bottom. If it goes belowe 1%, it probably means taht something is wrong with the setup, as for some reason the configuration found is very inefficient.*  
 **Not used**

- Upper beta limit  
*Physical limit. This is the ration of plasma pressure to magnetic pressure. It should be maximized for the best utilization of the magnetic field. MHD stability and edge stochastisation limits this values. Limit should be known from magnetic configuration calculations done outside PROCESS*  
`icc = 24`  
beta_max = 0.04  
*Depends on specific configuration, 4% seems to be optimistic guess for SQuID (for economic reasons, this will often be a stongly limiting factor)*

#### Radiation limits
 - Mean neutron load on vessel first wall upper limit (MW/M2)   
 *Physical limit. Defines lifetime of FW, blanket, VV, coils*  
`icc = 8`   
walalw = 1.5  
*Depends on blanket concept and thickness, can be a strongly limiting factor for the design space. TODO: needs full 3D neutronics simulations to ensure VV/coil lifetime. According to DEMO 2014. Max load for this version is 1.35 MW/m2, so it should be corrected according to peaking factor. Proxima assumed 4.05 as peak value for Stellaris, at a cost of only 4 years of full-power operation.*

- Radiation power density upper limit  
`icc = 17`  
*Ensures that the calculated total radiation power density does not exceed the total heating power to the plasma. It is recommended to have this constraint on as it is a plasma stability model. Without it, plasma can radiate more power than it produces, which would mean non-steady-state solution.* 

 - Divertor heat load upper limit (MW/M2)
`icc = 18`   
pflux_div_heat_load_max_mw = 10.0
*The exact limit depends on the divertor concept. In stellarators, divertor heat load is not limiting in most cases, unless very agressive assumptions are used for other constrains. 10 MW/M2 is considered avichable for reactor divertors.*

 - Radiation Wall load limit (MW/M2)   
 *Physical limit, cooling capability of FW. Should not limit the design in most cases (radiation load is much smaller than neutron load)*  
`icc = 67`   
maxradwallload = 1.2  
*For HCPB. This is a limit for total FW heat load, but radiation should make 95% of heat load. At increased pressure loss higher values could be possible. The problem with remaining 5% from charged particles is that they will be localized, so local overheating is possible (this would require 3D study in case of stellartors and is outside the PROCESS capabilities)*

#### Build limits
- Toroidal build consistency  
*Physical limit. Checks if coils don overlap toroidally.*  
`icc = 82`  
toroidalgap > tftort constraint, "tftort": "TF Coil Toroidal Thickness (M)",
*Calculated coil size compared vs current filament distance*

 - Radial build consistency   
 *Physical limit. Checks if the radial thickness of the components (blanked, VV, etc.) will fit into the machine dimensions.*  
`icc = 83`   
*Thickness of the blanked + shielding + vaccume vessel + plasma-wall distance*

#### Quench limits 
 - TF coil conduit stress upper limit (SCTF)  
 *Physical limit, Maximal Coil Stress on Ground insulation (approx.) This is the stress limit for the radial force which acts on the coil during operation. (Centering force is not limiting for stellartor build, as there is enough space for the supprot structure in the machine center)*  
`icc = 32`  
sig_tf_wp_max = 4.0e8  
*Maximal allowable Tresca stress. We keept the conservative assumption for stellarator, as we want to keep it safe in generic case. This values is averaged for winding-pack, so it depends on the exact structure (like the steel content in the winding pack). A 3D FEM simulation is needed to get the exact limit for specific case. For comparison: Default for tokamke is 6.0e8. Stellaris assumes 800MPa as a limit for steel, which is limited to 650MPa for WP load (based on FEM ananlysis). For DEMO, the structural material (TF nose and SC conduit) (Tresca) stress was not permitted to exceed 660 MPa, being the lower of 2/3 the yield stress and ½ the ultimate tensile stress for the cryogenic steel of choice.*


- ~~I_op / I_critical (TF coil) (SCTF)~~  
~~`icc = 33`~~  
  *Limit on the critical current density in relation to maximum current denisty allowed for superconductor is hardcoded into stellarator modeule, so this option should not be used!*  
 **Not used**
 
 - Dump voltage upper limit (SCTF) [kV]  
 *Physical limit - the maximum voltage which can be induced during the quench discharge. From the economic perspective, it is good to keep it low, as higher values require more insulation.*  
`icc = 34`  
 vdalw = 12.0  
 *PROCESS default is 20kV. 12.64kV was proposed for the old stellarator input. W7-X have 8kV, ITER 10kV. 12kV seems feasible with current technology, some proposals go to around 20kV. This is usualy not a limiting factor, so 12kV seems fine as a starting point. (for HELIAS, the ITER coils were taken as reference and technical specifications essentially copied)*  
 **TODO review after coil specification check**

  - J_winding pack/J_protection upper limit (SCTF)  
  *Physical limit - To ensure that J_op does not exceed the current density protection limit. It is not related to superconductor limits, but to the quench protection. Simple 0D model of the ohmic heating is used to define tmperature rise of the copper during quench. (copper resistance must be low enoughn to prevent high temperature rise during quench, which is a function of initial current, quench time and copper fraction in the coil).
`icc = 35`  

  - Dump time set by VV loads  
  *Physical limit - during the quench a current is induced in the vaccume vessel, which results in the stress from electromagnetic force acting on the VV. To limit induced current, a dump time can't be too short.*  
`icc = 65`  
 max_vv_stress = 9.3e7  
 The Allowable Peak Max Shear Stress In The Vacuum Vessel Due To Quench And Fast Discharge Of The TF Coils [Pa]  
 *This limit is inherited from the old Helias design. Depends on the VV geometry, steel material and norms used. Excat value should be checked in 3D FEM analysis. Note that this limit is calculated a bit diffrent than for the tokamaks. In stellarator, the model is basicly an extrapolation of the results obtained for W7-X. Such detailed studies of the other stellartor designs are not available at the moment.*

#### Other limits 
- ratio of particle to energy confinement times  
`icc = 62`  
 f_alpha_energy_confinement_min = 4
 *Proxima assumes 8 and a factor of 0.5 to supress helium denisty: "This assumption is informed by two considerations: first, by suppression of helium ash due to a positive ambipolar radial electric field, and second, by the effect that fast particles do not slow down directly in the core, but are radially displaced due to their finite drift orbit.*

 - ~~ECRH ignitability (denisty limit)~~  
~~`icc = 91`~~  
 *We don't use ECHR limits, becouse they are not strict and can be overcome with proper technology. Note that the model used is computationaly demanding and can increas execution time.*  
 **Not used**


### Iteration Variables

ixc = 2                 *bt
bt       = 5.5          *Toroidal field on axis (T)
boundl(2) = 4.0
boundu(2) = 10.0

- Plasma major radius (m)  
`ixc = 3 `  
rmajor   = 21.0   
boundl(3) = 2.  
boundu(3) = 25.  

- Volume averaged electron temperature [keV]  
`ixc = 4`  
te       = 7.0          
boundl(4) = 4.  
boundu(4) = 25.

- Electron density (/m3)  
`ixc = 6`               
dene = 2.0E20    
boundl(6) = 3.005E19  
boundu(6) = 3.005E20
*We need some constrain on the denisty. PROCESS can manage Sudo denisty, ECHR limit or limit by impurity radiation (https://ukaea.github.io/PROCESS/fusion-devices/stellarator/?h=sudo#density-limit). Sudo is not really relaible, it was more a heusristic limit for small stellerators. ECHR limit is tricki in PROCESS, at some point maybe I will write a script to check it automatically in pre/post processing in Python. Impurity radiation should work, but it is uncertain if it will limit denisty to realistic values. For now, the upper bound is set manually and it will be evaluated based on the first results. [in some HTS studies, very high values ~6x10^20 are used]*

- scaling factor on energy confinement times  
*This factor is used to scale the confinement law. It reflects the fact, that we expect better confinemnt in future machines, as in the experiments and simulations it is possible to reach better confinemnt scaling in the righ conditions. In principle, it is specific for the configuration and can be up to some degree taken into account in the design phase.*  
`ixc = 10`         
hfact = 1.0        
boundu(10) = 1.00  
*scenarios upper bound limit: conservative 1; advanced: 1.3*

- f-value (scaling factor) for net electric power  
*This is not strictly required, but it gives some flexibility to the power, so it is not necessary equal to the max value. We used it to easen convergence, but it is probably not necessary in most cases.*  
`ixc = 25`               
fpnetel = 1.0000        
boundl(25) = 0.99  
boundu(25) = 1.0

- Fraction TF coil critical current to operation current  
*In principle it can be made an iteration variable, but we don't reccoment this. This is a safty margin between operational current denisty and superconductors maximum current denisty coming from the material specification (maximum current density as a function of magnetic field and temperature is deined in the PROCESS for diffrent SC types). Margin defined as fraction of maximum current is can be defined alternatively as to temperature margin between operation temperature and temperatue for which maximum current denisty is defined*  
~~`ixc = 50`  ~~      
fiooic = 0.8  
~~boundu(50) = 0.9~~  
~~boundl(50) = 0.001~~  
**Fixed to 0.8**  
*For NbSn magnets it should correspond to 1,5- 2K of temperature safety margin. 0.8 was used by the Proxima for HTS magnets.*

- thermal alpha density / electron   
*This iteration variable is needed to vary helium content in plasma*
`ixc = 109`  
f_nd_alpha_electron: density  
boundl(109) = 0.0001  
boundu(109) = 0.4


- ~~Achievable Temperature of the ECRH at the ignition point~~  
~~`ixc = 169`~~  
~~te0_ecrh_achievable = 17.5 * keV~~   
~~boundl(169) = 4.~~  
~~boundu(169) = 35.~~  
 *We don't use ECHR limits, becouse they are not strict and can be overcome with proper technology. Note that the model used is computationaly demanding and can increas execution time.*  
**Not used**

- Copper Fraction Of Cable Conductor (TF Coils)  
*Copper fraction changes are necessary for optimization of WP dimensions vs quench protection requirements*  
`ixc = 59 `   
fcutfsu = 0.69  
boundu(59) = 0.90  
boundl(59) = 0.2

- Fast Discharge Time For TF Coil In Event Of Quench (S)  
*For self-consistancy of quench calculation*  
`ixc = 56`  
tdmptf  
boundl(56) = 1  
boundu(56) = 100.
*Usually quench duration is on the order of seconds, but as it should be self-consistantly defined, we apply wide margines to not hinder convergence*


### Physics Variables

* Density profile index    
`alphan   = 0.35`
*Values defined in W7-X. It can change in specific stellartor design, but is a reasonable guess for generic machine*

* Temperature profile index  
`alphat   = 1.20`    
*Values defined in W7-X. It can change in specific stellartor design, but is a reasonable guess for generic machine*

* Aspect ratio   
`aspect   = 11.1`   
*11.1 for SQuID, 12.3 for Helias* 

* Switch for ignition assumption (1: Ignited)  
`ignite   = 1`     
*This is necessary for burning plasma reactor*

* Switch for pedestal profiles (0: Parabolic Profiles)        
`ipedestal = 0`     
*We don't use pedestal model with stellartor. Pedestals as understood in tokamaks are not observed in W7-X.*

* Switch for radiation loss term usage in power balance (1: Total power lost is scaling power plus core radiation only)     
`i_rad_loss = 1`     
*Confinement time is power depandant, but the exact definition of the power which should be taken into equation is not clear in the case of a reactor. This option is some kind of a compromise between two limits, as we expect that we should take into account part of the radiated power, but taking all would lead to overestimation of the confinement time (at the radical case, radiated power would cancel out the heating power)* 

* Switch for energy confinement time scaling law (38: ISS04)    
`i_confinement_time = 38`  
*This is the confinement scaling which we consider most relaible*

* Synchrotron wall reflectivity factor  
`f_sync_reflect    = 0.6`  
*This value depned on the geometry and plasma facing components. In principle, higer value results in lower synchrotron looses. In some recent studies more optimistic values ver ptroposed, like 0.8 in the recent Stellaris analysis - but it also involved diffrend synchrotron rqadiation model. We consider 0.6 reasonably conservative, but the exact value would require detailed plasma + raytracing analysis in 3D geometry. As stellarators in general operate at higher denisty and lower temperature than tokamaks, synchrotron radiation is a minor issue in our case.*

* Ion temperature / electron temperature  
`tratio    = 0.95`          

* F-Value For Lower Limit On Taup/Taueff The Ratio Of Alpha Particle To Energy  
`falpha_energy_confinement = 1` 
*This can be used as alpha energy confinement safety factor. (not used in our analysis)*

* F-Value For Core Radiation Power Limit  
`fradpwr = 1`        
*This can be used as core radiation safety factor. (not used in our analysis)*

### Stellarator variables
* Switch for stellarator option (6: Use stella_config file)  
`istell   = 6`  
*This option mean we use external stellarator configuration file, which contains precalculated parameters for our magnetic configuration*

* Switch for stellarator auxiliary heating method (1: Electron cyclotron resonance heating)  
`isthtr   = 1`
*ECHR is considered the most reasonable option for stellartor reactor - it is efficient, operate in focused beams (low TBR impact) and can be easly transmitted, which means that the source can be outside the torus hall.*

### Build Variables

*From neutronic calculation of the HCPB blanket. Blanket is here understood as breeding zone, and manifold is treated as shielding. Inboard and outbourd sizes are the same as in DEMO, as we expect similar non-uniformity in stellarator. Inboard side value is used for radial build constraint. Outboard and inboard side define average value, which is used in mass and shielding calculations*

* Inboard blanket thickness (m)  
`dr_blkt_inboard  = 0.41`  

* Outboard blanket thickness (m)  
`dr_blkt_outboard  = 0.63`

* Cryostat thickness (m)  
`dr_cryostat    = 0.05`

* Inboard vacuum vessel thickness (tf coil / shield) (m)  
`dr_vv_inboard  = 0.6`

* Outboard vacuum vessel thickness (tf coil / shield) (m)  
`dr_vv_outboard = 0.6`

* Topside vacuum vessel thickness (tf coil / shield) (m)  
`d_vv_top = 0.6`

* Underside vacuum vessel thickness (tf coil / shield) (m)  
`d_vv_bot = 0.6`

* Gap between inboard vacuum vessel and tf coil (m)  
`dr_shld_vv_gap_inboard    = 0.025`

* Minimum gap between outboard vacuum vessel and TF coil (m)  
`gapomin  = 0.025`

*Distance from plasma to fisrt wall cannot be too small, as magnetic islands must fit into this space.*
* Gap between plasma and first wall; inboard side (m)  
`dr_fw_plasma_gap_inboard  = 0.3`

* Gap between plasma and first wall; outboard side (m)  
`dr_fw_plasma_gap_outboard  = 0.3`

* Inboard shield thickness (m)  
`dr_shld_inboard  = 0.3`  

* Outboard shield thickness (m)  
`dr_shld_outboard  = 0.3`

* Upper/lower shield thickness (m)  
`shldtth  = 0.3`

* Vertical gap between x-point and divertor (m)  
`vgap_xpoint_divertor     = 0.0`
*This is not applicable in stellarator model*


### Current Drive Variables
* ECH wall plug to injector efficiency  
`etaech   = 0.5`    
*Nowdays we can reach 50%. 60% is expected ECRH efficiency in the (not so far) future (according to Proxima)*

* Heating power not used for current drive (MW)     
`pheat    = 0.0`          
*We don't need curretn in stellarstor.*


### Divertor Variables
* Angle of incidence of field line on plate (rad)  
`anginc   = 0.03`        

* Temperature at divertor (eV)        
`tdiv     = 5.0`          
*It is used to calculate ion speed in SOL. It should be around 5-10 eV in detached scenario (Infinity Two assumes 2-10eV). We use 5eV, as this is reasonable whn considering plasma closer to divertor plates (in case if this will used also by some other models)*

* Perpendicular heat transport coefficient (m2/s)  
`xpertin  = 1.5`        

* Wetted Fraction Of The Divertor Area    
`fdivwet	= 0.333 `  
*As first approximation we can assume that divertor will be 3x as wide as strike line (The designs of the divertors for stellarators are still so immature that this is very rough estimate)*

* Relative radial field perturbation  
`bmn      = 0.001`  
*This should be approximetly equal to the flpitch value*

* Divertor heat load peaking factor  
`f_asym   = 1.1`

* Radiated power fraction in sol  
`f_rad    = 0.85`  
*In W7-X it was slightly above the limit. Proxima used 0.9 (at 0.9 in W7-X detachment front breached the LCFS). Keep in mind, that change from 0.9 to 0.85 results in 50% increase of the divertor direct heat load.*

* Island size fraction factor  
`f_w      = 0.5`  
*This parameter describes where in the island the divertor is placed (and so what fraction of island 'deepth' is effectively used). This value can change in specific divertor designs, but 0.5 represents good value for generic stellarator machine.*

* Field line pitch (rad)  
`flpitch  = 0.001`  
*Describes the radial displacement of a field line in the SOL along its arc-length and depends on the specific magnetic configuration. In principle it can vary between 1e-3 - 1e-4.*

* Rotational transform (reciprocal of tokamak q)  
`iotabar  = 1.0`
*This value is specific for magnetic configuration*

* Magnetic shear, derivative of iotabar  
`shear    = 0.5`
*This value is specific for magnetic configuration*

### FWBs Variables

*The values give are for HCPB blanket*

* Density of steel (kg/m3)  
`denstl   = 7800.0`  
*EUROFER97 steel density*

* Energy multiplication in blanket and shield  
`emult    = 1.35`  
*For beryllium breeder the nuclear analysis of the HCPB blanked gives 1.35 multiplication.*

* Electrical efficiency of primary coolant pumps  
`etahtp   = 1.0`  
*We don't calculate mechanical power for pumps, it is in fact electric power. So we use 1 conversion factor. This way we loose some recirculating heat, but this is a conservative assumption.*

* Beryllium fraction of blanket by volume  
`fblbe    = 0.3663`

* Lithium oxide fraction of blanket by volume  
`fblli2o  = 0.1491`

* Lithium lead fraction of blanket by volume  
`fbllipb  = 0.00`

* Stainless steel fraction of blanket by volume  
`fblss    = 0.0985`

* Vanadium fraction of blanket by volume  
`fblvd    = 0.00`

* Area fraction taken up by other holes (not used)  
`fhole    = 0.0`

* First wall coolant fraction  
`fwclfr   = 0.35`

* Coolant void fraction in blanket (blktmodel=0)  
`vfblkt   = 0.386`

* Coolant void fraction in shield  
`vfshld   = 0.40`

* Neutron Power Deposition Decay Length Of Blanket Structural Material [M] (Stellarators Only)  
`declblkt = 0.1`

* Neutron Power Deposition Decay Length Of First Wall Structural Material [M] (Stellarators Only)  
`declfw = 0.1`

* Neutron Power Deposition Decay Length Of Shield Structural Material [M] (Stellarators Only)  
`declshld = 0.056`

### Heat Transport Variables

* Switch for pumping power (0: User sets pump power directly)  
`i_coolant_pumping = 1`   
*This sets mechanical pumping power as a fraction of thermal power removed by coolant.*

* fraction of total blanket thermal power required to drive the blanket coolant pumps 
`fpumpblkt    = 0.033`   
*Estimated from reference pump power 80MW, for 2,4GW of thermal power. There is no distinction between blanket and FW heat in the blanket cooling concept.*

* fraction of total first wall thermal power required to drive the FW coolant pumps
`fpumpfw    = 0.033`   
*Estimated from reference pump power 80MW, for 2,4GW of thermal power. There is no distinction between blanket and FW heat in the blanket cooling concept.*

* fraction of total divertor thermal power required to drive the divertor coolant pumps
`fpumpdiv    = 0.107`   
*Estimated from reference pump power 14.5MW, 135MW of cooling for the divertor. (135 MW heat rejected from the divertor was calculated in Helias5b). Values are neglecting the divertor caseete cooling, but the highest pump power/cooling power is at the divertor plasma facing units.*

* Fraction Of Total Shield Thermal Power Required To Drive The Shield Coolant 
`fpumpshld    = 0.0`  
*Heat extracted from shield was less than 1MW, so we can neglect that. In the old version it was not inclueded in primary coolant pumps.* 

* Switch for power conversion cycle (2: user input thermal-electric efficiency)  
`i_thermal_electric_conversion = 2`

* Thermal to electric conversion efficiency; if seconday_cycle=2  
`eta_turbine    = 0.375` 

### Impurity Radiation Module

* Normalised radius defining the 'core' region  
`coreradius = 0.6`

* Hydrogen (fraction calculated by code)  
`fimp(1)  = 1.0`

* Helium (fraction calculated by code)  
`fimp(2)  = 0.1`

* Beryllium
`fimp(3)  = 0.0`

* Carbon
`fimp(4)  = 0.0`

* Nitrogen
`fimp(5)  = 0.0`

* Oxygen
`fimp(6)  = 0.0`

* Neon
`fimp(7)  = 0.0`

* Silicon
`fimp(8)  = 0.0`

* Argon
`fimp(9)  = 0.0`

* Iron
`fimp(10) = 0.0`

* Nickel
`fimp(11) = 0.0`

* Krypton
`fimp(12) = 0.0`

* Xenon
`fimp(13) = 0.0`

* Tungsten
`fimp(14) = 1.0E-5`

### Numerics

* Code operation switch (1: Optimisation, VMCON only)  
`ioptimz  = 1`

* Maximum number of VMCON iterations
*Number of iterations should depend on case specification. For well defined constrains, 10-30 iterations should be enough, but sime cases requires much more iterations and increasing the limit to 100 is not uncommon.*
`maxcal   = 50`

* Switch for figure-of-merit (6: cost of electricity). Positive number looks for minimum, negative for maximum.
`minmax   = 6`

* Name of the run (written in output file)  
`runtitle = SQuID`

### Tfcoil Variables

* Switch for superconductor material in tf coils
`i_tf_sc_mat = 1`
*We used 1 for LTS and 8 for HTS*

* Peak helium coolant temperature in TF coils and PF coils (K)  
`tftmp     = 4.5`  
*Critial temperature of TF coils is ≈6K. The safty margin is estimated as 0.5K, but we don't know the exact termperature distribution inside the coil. We assume that the peak temperature is up to 1K higher than inlet temperature (4,5K at inlet). The total margin is 1.5K, which we take into account by lowering the operating current with fiooic variable. So the peak temperature given at tftmp is the inlet temperature, otherwise we would double-count the safty margin. When the thermal analysis of the coil will be avaliable, this should be changed.*

* Coil temperature for cryogenic plant power calculation (K)  
`tmpcry    = 4.5`  
*This should be temperature achived in the cryoplant, which is the inlet temperature of the coolant in TFC. 4.5 K is the value from DEMO.*

* Dimension conductor area including steel and insulation.  
`t_turn_tf = 0.068`

* Conduit insulation thickness (m)  
`thicndut  = 0.0015`

* Thickness of steel around each conductor  
`thwcndut  = 0.006`

* Coolant fraction of TF coil leg (i_tf_sup=0)  
`vftf      = 0.3`

* Case thickness  
`thkcas    = 0.06`

* Insulation on top of winding pack  
`tinstf    = 0.01`


### Cost Variables

* 0: 1990 cost module, the 2015 does not work yet for stellarators  
`cost_model = 0`

* Allowable first wall/blanket neutron (MW-yr/m2)  
`abktflnc = 15.0`

* Allowable divertor heat fluence (MW-yr/m2)  
`adivflnc = 25.0`

* Total plant capacity fraction  
`cfactr   = 0.75`

* Diff between borrowing and saving interest rates  
`dintrt   = 0.00`

* Average cost of money for construction of plant  
`fcap0    = 1.15`

* Average cost of money for replaceable components  
`fcap0cp  = 1.06`

* Project contingency factor  
`fcontng  = 0.15`

* Fixed charge rate during construction  
`fcr0     = 0.065`

* Multiplier for nth of a kind costs  
`fkind    = 1.0`

* Switch for plant availability model (0: Use input value for cfactr)  
`iavail   = 0`

* Switch (0: treat blanket divertor first wall and fraction fcdfuel of CD equipment as capital cost)  
`ifueltyp = 0`

* Switch for net electric power calculations (1: Calculate MW electric and c-o-e)  
`ireactor = 1`

* Level of safety assurance switch (2: In-between)  
`lsa      = 2`

* Effective cost of money in constant dollars  
`discount_rate = 0.06`

* Plant life (years)  
`tlife    = 40.0`

* Unit cost for blanket vanadium ($/kg)  
`ucblvd   = 280.0`

* Cost of divertor blade ($)  
`ucdiv    = 5.0E5`

* Unit cost of maintenance equipment ($)  
`ucme     = 3.0E8`