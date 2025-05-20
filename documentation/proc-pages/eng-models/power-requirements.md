# Power Requirements

The main power flow is controlled by `power.py`. The main class `Power` controls this.
## Power requirements | `Power`

---

### PF power | `pfpwr()`

---

### AC power | `acpow()`

---

### Main power conversion | `power1()`

---

### Remainder power conversion | `power2()`

---

### Time dependent power requirements | `power3()`

---

### Cryogenic power requirements | `cryo()`

---

### Plant thermal efficiency | `plant_thermal_efficiency()`

`i_thermal_electric_conversion` : This switch controls the calculation of the thermal to electric conversion 
efficiency in the secondary cycle.

----------------

#### Use CCFE HCPB Model Value

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
\eta_{\text{thermal → electric}} = 0.411
$$


-------------------

#### Use CCFE HCPB Model Value with divertor

This model is set with `i_thermal_electric_conversion = 1`.

It can be used with water or helium primary coolants.

This model is the same as above but a penalty is applied as the coolant in the divertor has to operate at much lower temperature than the blanket, which may be the case because of the greater heat flux that has to be removed.

$$
\eta_{\text{thermal → electric}} = 0.411 -0.339f
$$

where $f$ is the fraction of heat to the divertor.

--------------

#### User input

This model is set with `i_thermal_electric_conversion = 3`.

The efficiency of the power generation cycle is input by the user via `etath`.
  
$$
\eta_{\text{thermal → electric}} = \texttt{etath}
$$

-----------------

#### Steam Rankine Cycle

This model is set with `i_thermal_electric_conversion = 3`.

It can be used with helium primary coolant.

$$
T_{\text{turbine,inlet}} = T_{\text{blkt,outlet}} - 20.0
$$
 
$$
\eta_{\text{thermal → electric}} = 0.1802 \ln{(T_{\text{turbine,inlet}})}-0.7823 - \Delta \eta \\
\text{for} \quad  657.15 \le T_{\text{turbine,inlet}} \le 915.15 \text{K}
$$

The peak divertor coolant temperature ($T_{\text{div,outlet}}$) is assumed to be $423.15 \text{K}$.

If the Rankine cycle is chosen and the primary coolant is water, it is assumed that the cycle is similar to that of pressurised water reactors currently in operation.


------------

#### Supercritical CO2 Brayton Cycle
  
This model is set with `i_thermal_electric_conversion = 4`.

It can be used with water or helium primary coolants.

A supercritical CO$_2$ Brayton cycle is assumed. The secondary cycle 
efficiency (`etath`) is calculated from the coolant outlet temperature using simple relations 
between temperature and efficiency from :

$$
T_{\text{turbine,inlet}} = T_{\text{blkt,outlet}} - 20.0
$$
 
$$
\eta_{\text{thermal → electric}} = 0.4347 \ln{(T_{\text{turbine,inlet}})}-2.5043 \\
\text{for} \quad  408.15 \le T_{\text{turbine,inlet}} \le 1023.15 \text{K}
$$

The peak divertor coolant temperature ($T_{\text{div,outlet}}$) is assumed to be the same as the turbine inlet ($T_{\text{turbine,inlet}}$).


The correlation of efficiency with temperature is derived from results of cycle modelling carried out by CCFE in collaboration with industry. The divertor heat is used in the main heat exchanger. The divertor heat is counted as primary heat, and is included in the calculation of the efficiency.


---

### Liquid metal breeder plant thermal efficiency | `plant_thermal_efficiency_2()`

---

### Resistive TF coil power requirements | `tfpwr()`

---

### Superconducting TF coil power requirements | `tfpwr()`

---

### TF coil power conversion system parameters | `tfcpwr()`

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






