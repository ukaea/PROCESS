# Power Requirements

The main power flow is controlled by `power.py`. The main class `Power` controls this.
## Power requirements | `Power`

### TF Coils

#### Resistive TF coil power requirements | `tfpwr()`

---

#### Superconducting TF coil power requirements | `superconducting_tf_power_iter_1988()`

---

#### TF coil power conversion system parameters | `superconducting_tf_power_iter_1988()`

---

### PF Coils

#### PF power | `pfpwr()`

---

#### AC power | `acpow()`

---

### Power Plant

#### Main power conversion | `component_thermal_powers()`

1: Calculate the electric wall plug power for the different coolant systems by dividing by the pump wall plug efficiencies

$$
P_{\text{pump,electric}} = \frac{P_{\text{pump}}}{\underbrace{\eta_{\text{pump}}}_{\texttt{eta_coolant_pump_electric}}}
$$

This is done for the first wall, blanket, vacuum vessel shield, divertor and blanket secondary breeder coolant (if present).

2: The total mechanical pumping required for all the coolants is calculated:

$$
\underbrace{P_{\text{pump,total}}}_{\texttt{p_coolant_pump_total_mw}} = P_{\text{pump,FW-Blkt}} + P_{\text{pump,blkt-secondary}} \\
  + P_{\text{pump,shield}} + P_{\text{pump,div}}
$$

3: The total electric power fo the coolant pumps is calculated:

$$
\underbrace{P_{\text{pump total,electric}}}_{\texttt{p_coolant_pump_elec_total_mw}} = \sum_{n}{P_{\text{pump,electric},n}}
$$

4: The electrical heat loss due to pump inefficiencies is calculated:

$$
\underbrace{P_{\text{pump loss,total}}}_{\texttt{p_coolant_pump_loss_total_mw}} = P_{\text{pump total, electric}} - P_{\text{pump,total}}
$$

5: The electrical heat loss in the heating and current drive power supplies is calculated:

$$
\overbrace{P_{\text{HCD, electric-loss}}}^{\texttt{p_hcd_electric_loss_mw}} = P_{\text{HCD, electric}} + P_{\text{HCD, injected}}
$$

5: If their is a secondary breeder/coolant loop in the blanket (`i_blkt_dual_coolant == 2`) then the heat deposited in the secondary coolant is the fraction of the total blanket nuclear heat and the mechanical pump power.

$$
P_{\text{blkt-breeder-heat}} = \left(P_{\text{blkt,nuclear}} \times \texttt{f_nuc_pow_bz_liq}\right) + P_{\text{pump,blkt-secondary}}
$$

6: If `i_blkt_dual_coolant == 1` the secondary breeder is only pumped for tritium extraction and not cooling so: 

$$
P_{\text{blkt-breeder-heat}} =  P_{\text{pump,blkt-secondary}}
$$

7: The total heat deposited in the first wall and blanket is simply:

$$
\overbrace{P_{\text{FW-Blkt, heat}}}^{\texttt{p_fw_blkt_heat_deposited_mw}} = \\
 \underbrace{\left[P_{\text{FW, nuclear}} + P_{\text{FW,}\gamma} + P_{\text{FW, pump}} + \underbrace{P_{\text{NB, orbit-loss}} + P_{\text{NB, shine-through}}}_{\text{Neutral beam effects (if present)}} + P_{\alpha,\text{loss}}\right]}_{\texttt{p_fw_heat_deposited_mw}} \\
 + \underbrace{\left[P_{\text{Blkt, nuclear}} + P_{\text{Blkt, pump}}\right]}_{\texttt{p_blkt_heat_deposited_mw}}
$$

  - $P_{\text{FW, nuclear}}$ & $P_{\text{Blkt, nuclear}}$ is the nuclear heating from neutron interaction (which includes the energy multiplication (`f_p_blkt_multiplication`) for the blanket.)
  - $P_{\text{FW,}\gamma}$ is the photon radiation incident on the FW (`p_fw_rad_total_mw`).
  - $P_{\alpha,\text{loss}}$ is the plasma lost alpha power (`p_fw_alpha_mw`)


8: The thermal power deposited in the shields is calculated:

$$
\overbrace{P_{\text{Shield, heat}}}^{\texttt{p_shld_heat_deposited_mw}} = \underbrace{P_{\text{CP shield, nuclear}}}_{\text{ST Centrepost nuclear heating (if present)}} \\
 + P_{\text{VV-Shield, nuclear}} + P_{\text{VV-Shield, pump}}
$$

9: The thermal power deposited in the divertor(s) is calculated:

$$
\overbrace{P_{\text{div, heat}}}^{\texttt{p_div_heat_deposited_mw}} = \underbrace{P_{\text{plasma,sep}}}_{\texttt{p_plasma_separatrix_mw}} + P_{\text{div, nuclear}} + P_{\text{div,}\gamma} + P_{\text{div, pump}}
$$

10: The thermal to electric conversion efficiency for the turbine is calculated with the [`plant_thermal_efficiency()`](#plant-thermal-efficiency--plant_thermal_efficiency)

The same is done for the liquid breeder but with plant_thermal_efficiency_2

11: The primary thermal power used to generate electricity is calculated:

$$
\overbrace{P_{\text{plant,primary-heat}}}^{\texttt{p_plant_primary_heat_mw}} = P_{\text{FW-Blkt, heat}} + \underbrace{P_{\text{Shield, heat}}}_{\text{If } \texttt{i_shld_primary_heat = 1}} + P_{\text{div, heat}}
$$

12: The secondary thermal powers for the shield and H&CD are calculated:

$$
\overbrace{P_{\text{Shield,secondary-heat}}}^{\texttt{p_shld_secondary_heat_mw}} = \underbrace{P_{\text{Shield, heat}}}_{\texttt{i_shld_primary_heat = 0}}
$$

$$
\overbrace{P_{\text{HCD,secondary-heat}}}^{\texttt{p_hcd_secondary_heat_mw}} = P_{\text{HCD, nuclear}} + P_{\text{HCD,}\gamma}
$$

13: The number of primary heat exchangers is calculated as follows:

$$
N_{\text{PHX}} = \left\lceil \dfrac{P_{\text{plant,primary-heat}}}{1000} \right\rceil
$$

-----------

#### Electric power production | `plant_electric_production()`

1: Calculate the total base plant facility load as a combination of base load and the power needed per $\text{m}^2$ of floor space

$$
\overbrace{P_{\text{base,total}}}^\texttt{p_plant_electric_base_total_mw} = \overbrace{P_{\text{base}}}^\texttt{p_plant_electric_base} \\
+ \underbrace{A_{\text{floor,effective}}}_{\texttt{a_plant_floor_effective}} \times \underbrace{q_{\text{floor,effective}}}_{\texttt{pflux_plant_floor_electric}}
$$

2: The electric demands of the plant core systems are calculated

$$
\overbrace{P_{\text{core systems}}}^\texttt{p_plant_core_systems_elec_mw} = P_{\text{base,total}} 
\\ + \overbrace{P_{\text{cryo plant}}}^\texttt{p_cryo_plant_electric_mw} + \overbrace{P_{\text{tritium plant}}}^\texttt{p_tritium_plant_electric_mw}
\\ + \overbrace{P_{\text{TF}}}^\texttt{p_tf_electric_supplies_mw} + \overbrace{P_{\text{PF}}}^\texttt{p_pf_electric_supplies_mw}
\\ + \overbrace{P_{\text{vacuum pumps}}}^\texttt{vachtmw} + \underbrace{\overbrace{P_{\text{CP pumps}}}}_{\text{If present}}^\texttt{p_cp_coolant_pump_elec_mw}
$$

3: The total secondary heat, which is the low-grade thermal power not used for power production is calculated

$$
\overbrace{P_{\text{secondary heat}}}^\texttt{p_plant_secondary_heat_mw} = P_{\text{core systems}}
\\ +  \overbrace{P_{\text{HCD, electric loss}}}^\texttt{p_hcd_electric_loss_mw} + \overbrace{P_{\text{pump, electric loss}}}^\texttt{p_coolant_pump_loss_total_mw}
\\ + \overbrace{P_{\text{div, secondary heat}}}^\texttt{p_div_secondary_heat_mw } + \underbrace{\overbrace{P_{\text{shield, secondary heat}}}}_{\text{If} \ \texttt{i_shld_primary_heat = 0}}^\texttt{p_shld_secondary_heat_mw}
\\ + \overbrace{P_{\text{HCD, secondary heat}}}^\texttt{p_hcd_secondary_heat_mw} + \overbrace{P_{\text{TF, nuclear}}}^\texttt{p_tf_nuclear_heat_mw}
$$

4: The gross-electric power produced is calculated

$$
\overbrace{P_{\text{gross, electric}}}^\texttt{p_plant_electric_gross_mw} = \overbrace{P_{\text{plant, primary-heat}}}^\texttt{p_plant_primary_heat_mw} \times \overbrace{\eta_{\text{turbine}}}^\texttt{eta_turbine}
\\ + \underbrace{\overbrace{P_{\text{liquid breeder heat}}}^\texttt{p_blkt_liquid_breeder_heat_deposited_mw} \times \overbrace{\eta_{\text{liquid,turbine}}}}_{\text{If present}}^\texttt{etath_liq}
$$

5: The thermal power lost in the turbine is calculated

$$
\overbrace{P_{\text{turbine,loss}}}^\texttt{p_turbine_loss_mw} = \overbrace{P_{\text{plant, primary-heat}}}^\texttt{p_plant_primary_heat_mw} \times (1-\overbrace{\eta_{\text{turbine}}}^\texttt{eta_turbine})
$$

6: The required recirculated electric power is calculated as the core system power and then the HCD and coolant pump electric powers that would also be on during a pulse.

$$
\overbrace{P_{\text{recirc, electric}}}^\texttt{p_plant_electric_recirc_mw} = P_{\text{core systems}}
\\ + \overbrace{P_{\text{HCD electric}}}^\texttt{p_hcd_electric_total_mw} + \overbrace{P_{\text{coolant pumps, electric}}}^\texttt{p_coolant_pump_elec_total_mw}
$$


7: The net-electric power is found by the different of gross and net

$$
\overbrace{P_{\text{net, electric}}}^\texttt{p_plant_electric_net_mw} = P_{\text{gross, electric}} - P_{\text{recirc, electric}}
$$


8: The recirculated power fraction is then quickly found as

$$
\overbrace{f_{\text{recirc}}}^\texttt{f_p_plant_electric_recirc} = \frac{P_{\text{gross, electric}} - P_{\text{net, electric}}}{P_{\text{gross, electric}}}
$$


---

#### Plant thermal efficiency | `plant_thermal_efficiency()`

`i_thermal_electric_conversion` : This switch controls the calculation of the thermal to electric conversion 
efficiency in the secondary cycle.

----------------

##### Use CCFE HCPB Model Value

This model is set with `i_thermal_electric_conversion = 0`.

It can be used with water or helium primary coolants.

The efficiency of the power generation cycle is set to a single value 
obtained from previous cycle modelling studies. The heat deposited in the Toroidal Field coils 
divertor coolant is assumed to be at such low temperature that it cannot be used for power 
generation and is dumped to the environment.

The resulting thermal efficiencies used are taken from studies that modelled Rankine cycles for
the different options of a helium-cooled primary circuit with a top temperature of 500℃,(For historical reasons in both cases the divertor was cooled by water with a top temperature of 150℃).
Hence, no variation of efficiency with primary coolant temperature
is possible using the simplified model; indeed, no temperatures
are even considered in the model.

$$
\eta_{\text{turbine}} = 0.411
$$


-------------------

##### Use CCFE HCPB Model Value with divertor

This model is set with `i_thermal_electric_conversion = 1`.

It can be used with water or helium primary coolants.

This model is the same as above but a penalty is applied as the coolant in the divertor has to operate at much lower temperature than the blanket, which may be the case because of the greater heat flux that has to be removed.

$$
\eta_{\text{turbine}} = 0.411 -0.339f
$$

where $f$ is the fraction of heat to the divertor.

--------------

##### User input

This model is set with `i_thermal_electric_conversion = 3`.

The efficiency of the power generation cycle is input by the user via `eta_turbine`.
  
$$
\eta_{\text{turbine}} = \texttt{eta_turbine}
$$

-----------------

##### Steam Rankine Cycle

This model is set with `i_thermal_electric_conversion = 3`.

It can be used with helium primary coolant.

$$
T_{\text{turbine,inlet}} = T_{\text{blkt,outlet}} - 20.0
$$
 
$$
\eta_{\text{turbine}} = 0.1802 \ln{(T_{\text{turbine,inlet}})}-0.7823 - \Delta \eta \\
\text{for} \quad  657.15 \le T_{\text{turbine,inlet}} \le 915.15 \text{K}
$$

The peak divertor coolant temperature ($T_{\text{div,outlet}}$) is assumed to be $423.15 \text{K}$.

If the Rankine cycle is chosen and the primary coolant is water, it is assumed that the cycle is similar to that of pressurised water reactors currently in operation.


------------

##### Supercritical CO2 Brayton Cycle
  
This model is set with `i_thermal_electric_conversion = 4`.

It can be used with water or helium primary coolants.

A supercritical CO$_2$ Brayton cycle is assumed. The secondary cycle 
efficiency (`eta_turbine`) is calculated from the coolant outlet temperature using simple relations 
between temperature and efficiency from :

$$
T_{\text{turbine,inlet}} = T_{\text{blkt,outlet}} - 20.0
$$
 
$$
\eta_{\text{turbine}} = 0.4347 \ln{(T_{\text{turbine,inlet}})}-2.5043 \\
\text{for} \quad  408.15 \le T_{\text{turbine,inlet}} \le 1023.15 \text{K}
$$

The peak divertor coolant temperature ($T_{\text{div,outlet}}$) is assumed to be the same as the turbine inlet ($T_{\text{turbine,inlet}}$).


The correlation of efficiency with temperature is derived from results of cycle modelling carried out by CCFE in collaboration with industry. The divertor heat is used in the main heat exchanger. The divertor heat is counted as primary heat, and is included in the calculation of the efficiency.


---

#### Liquid metal breeder plant thermal efficiency | `plant_thermal_efficiency_2()`

---




### Cryogenic power requirements | `cryo()`

---





Figure 1 shows a simplified description of the power flow. 

<figure>
    <center>
    <img src="../../images/Overall-power-flow.png" alt="Overall power flow" 
    title="Power flows" 
    width="650" height="100" />
    <br><br>
    <figcaption><i>Figure 1: Power flows
    </i></figcaption>
    <br>
    </center>
</figure>

Some details of the auxiliary systems are as follows.

`tfcpwr` calculates the TF coil power conversion system parameters.  Only the steady-state power consumption for a superconducting TFC system is described here.

The TF current is carried from the power supplies to the reactor by room-temperature aluminium busbars, organised in $N_{circuit}$ circuits.  The total length of the busbars is (somehwat arbitrarily) given by

$$
L_bus = 8 \pi R_0 + (1 + N_{circuit}) (12 R_0 + 80) 
$$

The resistivity of the busbar is 2.62e-8 ohm.m (0.0262 ohm.mm²/m) (hard-coded).

"TF coil resistive power" (`rpower`) includes the dissipation of the cryogenic current leads (assumed to be resistive).

The AC power required is determined by the efficiency of the coil power supply: `etatf` (default = 90%).






