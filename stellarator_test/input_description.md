# Input for SQuID, LTS

### Equaltities

- Global power balance   
`icc = 2`   
(consistency equation, specific value not needed)

- Net electric power lower limit [MW]  
*Necessary for cost optimization to avoid minimal working radius*  
`icc = 16`  
pnetelin = 1000  
*Should be defined by CAPEX/LCOE comparison (historically often 1 GW was taken as base value)*



### Inequaltities

#### Beta limits
- ~~Lower beta limit~~   
~~`icc = 84`~~   
~~beta_min = 0.01~~  
*There is probably no need for this, we want to maximize beta anyway. There should be no physical hard limit restricting it from the bottom.*  
 **Not used**

- Upper beta limit  
*Physical limit. This is the ration of plasma pressure to magnetic pressure, it should be maximized for the best utilization of the magnetic field. MHD stability and  edge stochastisation limits this values.*  
`icc = 24`  
beta_max = 0.04  
*Depends on specific configuration, 4% seems to be optimistic guess for SQuID (for economic reasons, this will often be a stongly limiting factor)*

#### Radiation limits
 - Neutron wall load upper limit (MW/M2)   
 *Physical limit. Defines lifetime of FW, blanket, VV, coils [It lacks precise description in documentation, I assume it is an average (there is usually mention for peak values)]*  
`icc = 8`   
walalw = 1.07  
*Depends on blanket concept and thickness, can be a strongly limiting factor for the design space. TODO: needs full 3D neutronics simulations to ensure VV/coil lifetime. According to DEMO 2014. Max load for this version is 1.35 MW/m2, so it should be corrected according to peaking factor. Proxima assumed 4.05 as peak value for Stellaris, at a cost of only 4 years of full-power operation.*

- Radiation power density upper limit  
`icc = 17`  
*Ensures that the calculated total radiation power density does not exceed the total heating power to the plasma. It is recommended to have this constraint on as it is a plasma stability model* 

 - Radiation Wall load limit (MW/M2)   
`icc = 18`   
maxradwallload = 10.0

 - Radiation Wall load limit (MW/M2)   
 *Physical limit, cooling capability of FW. Should not limit the design in most cases (radiation load is much smaller than neutron load)*  
`icc = 67`   
maxradwallload = 1.2  
*For HCPB. This is a limit for total FW heat load, but radiation should make 95% of heat load. At increased pressure loss higher values could be possible. The problem with remaining 5% from charged particles is that they will be localized, so local overheating is possible (to be considered in the future)*

#### Build limits
- Toroidal build consistency  
*Physical limit. Checks if coils don overlap toroidally.*  
`icc = 82`  
toroidalgap > tftort constraint, "tftort": "TF Coil Toroidal Thickness (M)",
Calculated coil size compared vs current filament distance

 - Radial build consistency   
 *Physical limit. Checks if the radial thickness of the components (blanked, VV, etc.) will fit into the machine dimensions.*  
`icc = 83`   
Thickness of the blanked + shielding + vaccume vessel + plasma-wall distance

#### Quench limits 
 - TF coil conduit stress upper limit (SCTF)  
 *Physical limit, Maximal Coil Stress on Ground insulation (approx.) Up to my understanding, this is the stress limit for the lorenz force which acts on the coil during operation.*  
`icc = 32`  
sig_tf_wp_max = 4.0e8  
Maximal allowable Tresca stress  
*Inherited from Process documentation. Default for tokamke is 6.0e8, so this seems conservative. Stellaris papaer assumes 800MPa as a limit for steel, which is limited to 650MPa for WP load (based on FEM ananlysis). For DEMO, the structural material (TF nose and SC conduit) (Tresca) stress was not permitted to exceed 660 MPa, being the lower of 2/3 the yield stress and ½ the ultimate tensile stress for the cryogenic steel of choice. 400MPa seems a conservative choice.*


- ~~I_op / I_critical (TF coil) (SCTF)~~  
  *Physical limit - Jop must not exceed the critical value Jcrit. Iteration variable 50 must be active (fiooic). The current density margin can be set using the upper bound of fiooic. This seems to be typical PROCESS mess. It was not used in the old inputs. Seems redundant with icc=35. In current setup it cause simulation to fail on fisrt iteration*  
~~`icc = 33`~~  
 **Not used**
 
 - Dump voltage upper limit (SCTF) [kV]  
 *Physical limit - the maximum voltage which can be induced during the quench discharge. From the economic perspective, it is good to keep it low, as higher values require more insulation.*  
`icc = 34`  
 vdalw = 12.0  
 *PROCESS default is 20kV. 12.64kV is proposed for stellarator input. W7-X have 8kV, ITER 10kV. 12kV seems feasible with current technology, some proposals go to around 20kV. This is usualy not a limiting factor, so 12kV seems fine as a starting point. (for HELIAS, the ITER coils were taken as reference and technical specifications essentially copied)*

  - J_winding pack/J_protection upper limit (SCTF)  
  *Physical limit - To ensure that J_op does not exceed the current density protection limit (documentation suggests that the limit is defined by the superconducting material selection. See the comment for icc=33*  
`icc = 35`  


  - Dump time set by VV loads  
  *Physical limit - during the quench a current is induced in the vaccume vesse, which results in the stress from electromagnetic force acting on the VV. To limit induced current, a dump time can't be too short.*  
`icc = 65`  
 max_vv_stress = 9.3e7  
 The Allowable Peak Max Shear Stress In The Vacuum Vessel Due To Quench And Fast Discharge Of The TF Coils [Pa]  
 *This limit is inherited from PROCESS. It is usualy not limiting factor for a design. Depends on the steel material and norms used, but seems realistic. Some safty factor due to more complex stellarator geometry would be desired. DEMO propose 91 MPa in magnets design guide.*

#### Other limits 
- ratio of particle to energy confinement times  
`icc = 62`  
 f_alpha_energy_confinement_min = 4
 *Proxima assumes 8 and a factor of 0.5 to supress helium denisty: "This assumption is informed by two considerations: first, by suppression of helium ash due to a positive ambipolar radial electric field, and second, by the effect that fast particles do not slow down directly in the core, but are radially displaced due to their finite drift orbit.*

 - ~~ECRH ignitability (denisty limit)~~  
 *Rejected, becouse ECHR are not strict limits and can be overcome with proper technology*  
~~`icc = 91`~~  
 **Not used**


### Iteration Variables

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
*In the base case it can be switched off, but the impact on the convergence is uncertain for now. It is a strongly limiting factor and we expect possible inprovement in the confinement, so it is worth to take into account an increased value to show the 'advanced' options.*  
`ixc = 10`         
hfact = 1.0        
boundu(10) = 1.00  
*scenarios: conservative 1; advanced: 1.3*

- f-value (scaling factor) for net electric power  
*I am not sure if this is really necessary, but it puts some flexibility to the power, so it is not necessary equal to the max value. If in calculations I will get always 1.0, I will turn it off.*  
`ixc = 25`               
fpnetel = 1.0000        
boundl(25) = 0.98  
boundu(25) = 1.0

- Fraction TF coil critical current to operation current  
*Following the approach of Jorrit, I keep it fixed (makes calculation easier with less iteration variables)*  
`ixc = 50`        
fiooic = 0.8  
~~boundu(50) = 0.9~~  
~~boundl(50) = 0.001~~  
**Fixed to 0.8 in the scan**  
*0,8 was used by the Proxima. For NbSn magnets it should correspond to 1,5- 2K of temperature margin.*

- thermal alpha density / electron   
`ixc = 109`  
f_nd_alpha_electron: density  
boundl(109) = 0.0001  
boundu(109) = 0.4

- ~~Achievable Temperature of the ECRH at the ignition point~~  
~~`ixc = 169`~~  
~~te0_ecrh_achievable = 17.5 * keV~~   
~~boundl(169) = 4.~~  
~~boundu(169) = 35.~~  
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



### Physics Variables

* Density profile index    
`alphan   = 0.35`

* Temperature profile index  
`alphat   = 1.20`    

* Aspect ratio   
`aspect   = 11.1`   
*11.1 for SQuID, 12.3 for Helias* 

* Toroidal field on axis (T)  
`bt       = 7.0`     

* Switch for ignition assumption (1: Ignited)  
`ignite   = 1`     

* Switch for pedestal profiles (0: Parabolic Profiles)        
`ipedestal = 0`     

* Switch for radiation loss term usage in power balance (1: Total power lost is scaling power plus core radiation only)     
`i_rad_loss = 1`     
*Confinement time is power depandant, but the exact definition of the power which should be taken into equation is not clear in the case of a reactor. This option is some kind of a compromise between two limits, as we expect that we should take into account part of the radiated power, but taking all would lead to overestimation of the confinement time (at the radical case, radiated power would cancel out the heating power)* 

* Switch for energy confinement time scaling law (38: ISS04)    
`i_confinement_time = 38`  

* ~~Plasma separatrix elongation~~  
~~`kappa    = 1.001`~~        
*Shouldn't be needed. I don't see any effect after removing this.*

* Synchrotron wall reflectivity factor  
`f_sync_reflect    = 0.6`  
*Assumption inherited from PROCESS. The lower the value, the higer the radiation loss. Proxima used 0.8 and a bit diffrent radiation model.*

* Ion temperature / electron temperature  
`tratio    = 0.95`        

* ~~High-Z impurity switch (0: Iron)~~  **(not used)**  
~~`*zfear    = 0`~~          

* F-Value For Lower Limit On Taup/Taueff The Ratio Of Alpha Particle To Energy  
`falpha_energy_confinement = 1` 

* F-Value For Core Radiation Power Limit  
`fradpwr = 1`        

### Stellarator variables
* Switch for stellarator option (6: Use stella_config file)  
`istell   = 6`  


* Switch for stellarator auxiliary heating method (1: Electron cyclotron resonance heating)  
`isthtr   = 1`


### Build Variables

* Inboard blanket thickness (m)  
`dr_blkt_inboard  = 0.41`  
*From neutronic calculation of the HCPB blanket. Blanket is here understood as breeding zone, and manifold is treated as shielding. Inboard and outbourd sizes are the same as in DEMO, to as we expect similar non-uniformity in stellarator.*

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

* Gap between plasma and first wall; inboard side (m)  
`dr_fw_plasma_gap_inboard  = 0.3`

* Gap between plasma and first wall; outboard side (m)  
`dr_fw_plasma_gap_outboard  = 0.3`

* Inboard shield thickness (m)  
`dr_shld_inboard  = 0.3`  
*From neutronic calculation of the HCPB blanket. Blanket is here understood as breeding zone, and manifold is treated as shielding.*

* Outboard shield thickness (m)  
`dr_shld_outboard  = 0.3`

* Upper/lower shield thickness (m)  
`shldtth  = 0.3`

* Vertical gap between x-point and divertor (m)  
`vgap_xpoint_divertor     = 0.0`


### Current Drive Variables
* ECH wall plug to injector efficiency  
`etaech   = 0.5`    
*Nowdays we can reach 50%. In the future possible to 60%: This is expected ECRH efficiency in the (not so far) future (according to Proxima)*

* Heating power not used for current drive (MW)     
`pheat    = 0.0`          


### Divertor Variables
* Angle of incidence of field line on plate (rad)  
`anginc   = 0.03`        

* Switch for divertor zeff model (1: input)  
`divdum   = 1`      
*We use prescribed Z_eff in divertor (we don't have any model for stellarator right now)*

* Temperature at divertor (eV)        
`tdiv     = 5.0`          
*I think it is used only to calculate ion speed in SOL. It should be around 5-10 eV in detached scenario (Infinity Two assumes 2-10eV). In case if this is used also by some other models, 5eV closer to divertor plates is also reasonable.*

* Perpendicular heat transport coefficient (m2/s)  
`xpertin  = 1.5`        

* Zeff in the divertor region (if divdum /= 0)    
`zeffdiv  = 2.5 `  

* Wetted Fraction Of The Divertor Area    
`fdivwet	= 0.333 `  
*As first approximation we can assume that divertor will be 3x as wide as strike line (The designs of the divertors are still so immature that this is currently only a guess.)*

* Relative radial field perturbation  
`bmn      = 0.001`  
*This should be approximetly equal to the flpitch value*

* Divertor heat load peaking factor  
`f_asym   = 1.1`

* Radiated power fraction in sol  
`f_rad    = 0.85`  
*In W7-X it was slightly above the limit. Proxima used 0.9 (at 0.9 in W7-X detachment front breached the LCFS)*

* Island size fraction factor  
`f_w      = 0.5`  
*According to sprocess paper, this should be around 0.5. In the old input I found 0.6, but 0.5 should be more conservative. This parameter describes where in the island the divertor is placed (and so what fraction of island 'deepth' is effectively used)*

* Field line pitch (rad)  
`flpitch  = 0.001`  
*Describes the radial displacement of a field line in the SOL along its arc-length and depends on the specific magnetic configuration. Process paper suggest 1e-3 - 1e-4.

* Rotational transform (reciprocal of tokamak q)  
`iotabar  = 1.0`

* Magnetic shear, derivative of iotabar  
`shear    = 0.5`

### FWBs Variables

* Density of steel (kg/m3)  
`denstl   = 7800.0`  
*EUROFER97 steel density*

* Energy multiplication in blanket and shield  
`emult    = 1.35`  
*For beryllium breeder the nuclear analysis of the HCPB blanked gives 1,35 multiplication.*

* Electrical efficiency of primary coolant pumps  
`etahtp   = 1.0`  
*We don't calculate mechanical power for pumps, it is in fact electric power. So we use 1 conversion factor. This way we loose some recirculating heat, but this is a conservative assumption. Electrica efficiency in the blanket analysis was 90%*

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

* Switch for pumping power (0: User sets pump power directly)  
`primary_pumping = 1`   
*Before prescribed pump powers were used. This is not the best approach, especially if I will change net power. i_coolant_pumping=1 seems better, it sets mechanical pumping power as a fraction of thermal power removed by coolant.*

* ~~Blanket coolant mechanical pumping power (MW)~~  
~~`htpmw_blkt = 120.0`~~

* ~~First wall coolant mechanical pumping power (MW)~~  
~~`htpmw_fw = 56.0`~~

* ~~Divertor coolant mechanical pumping power (MW)~~  
~~`htpmw_div = 24.0`~~

* fraction of total blanket thermal power required to drive the blanket coolant pumps 
`fpumpblkt    = 0.033`   
*Estimated from reference pump power 80MW, for 2,4GW of thermal power. There is no distinction between blanket and FW heat in the blanket cooling concept. A test is needed to check consistency with previous approach*

* fraction of total first wall thermal power required to drive the FW coolant pumps
`fpumpfw    = 0.033`   
*Estimated from reference pump power 80MW, for 2,4GW of thermal power. There is no distinction between blanket and FW heat in the blanket cooling concept. A test is needed to check consistency with previous approach*

* fraction of total divertor thermal power required to drive the divertor coolant pumps
`fpumpdiv    = 0.107`   
*Estimated from reference pump power 14.5MW, 135MW of cooling for the divertor. (135 MW heat rejected from the divertor was calculated in Helias5b). Values are neglecting the divertor caseete cooling, but the highest pump power/cooling power is at the divertor plasma facing units.*

* Fraction Of Total Shield Thermal Power Required To Drive The Shield Coolant 
`fpumpshld    = 0.0`  
*Heat extracted from shield in Helias5b example was less than 1MW, so we can probably neglect that. In the old version it was not inclueded in primary coolant pumps.* 

* Switch for power conversion cycle (2: user input thermal-electric efficiency)  
`secondary_cycle = 2`

* Thermal to electric conversion efficiency; if seconday_cycle=2  
`etath    = 0.375` 

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


### Impurity Radiation Module

* ~~Switch for impurity radiation model~~  
~~`imprad_model = 1`~~

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

* Maximum number of VMCON iterations (for well defined constrains, 10-30 iterations should be enough)  
`maxcal   = 50`

* Switch for figure-of-merit (1: Min Major Radius, 7: Min Capital Cost). Positive number looks for minimum, negative for maximum.  Full list of variables at /PROCESS/documentation/figure_of_merit.md  
`minmax   = 1`

* Convergence important to start with 10^-4, 10^-3 is NOT sufficient to match the constraint resonably well!!  
`epsfcn   = 0.0001`

* Name of the run (written in output file)  
`runtitle = SQuID`

### Tfcoil Variables

* Switch for superconductor material in tf coils (1: ITER Nb3Sn)  
`i_tf_sc_mat = 1`

* Peak helium coolant temperature in TF coils and PF coils (K)  
`tftmp     = 4.5`  
*Critial temperature of TF coils is ≈6K. The safty margin is estimated as 0.5K, but we don't know the exact termperature distribution inside the coil. We assume that the peak temperature is up to 1K higher than inlet temperature (4,5K at inlet). The total margin is 1.5K, which we take into account by lowering the operating current with fiooic variable. So the peak temperature given at tftmp is the inlet temperature, otherwise we would double-count the safty margin. When the thermal analysis of the coil will be done, this should be changed.*

* Coil temperature for cryogenic plant power calculation (K)  
`tmpcry    = 4.5`  
*This should be temperature achived in the cryoplant, which is the inlet temperature of the coolant in TFC. 4,5 K is the value from DEMO.*

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

### ~~Pfcoil Variables~~

* ~~PF coil vertical positioning adjuster~~  
~~zref(1)  = 3.6~~  
~~zref(2)  = 1.2~~  
~~zref(3)  = 2.5~~  
~~zref(4)  = 1.0~~  
~~zref(5)  = 1.0~~  
~~zref(6)  = 1.0~~  
~~zref(7)  = 1.0~~  
~~zref(8)  = 1.0~~  

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