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






