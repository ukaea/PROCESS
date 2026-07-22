# Plasma Fuelling | `PlasmaFuelling()`

## Particle balance

The control of fuelling is governed by 4 key particle flux equations for each of the primary fuel species and the helium ash, $\alpha$.

$$
\frac{dn_{\text{T}}}{dt} = f_{\text{fuelling,T}}\eta_{\text{fuelling}}\Gamma_{\text{fuel}} + \Gamma_{\text{D+D} \rightarrow \text{T}} - \Gamma_{\text{D+T}} - \frac{N_{\text{T}}}{\tau_{\text{T}}^*}
$$

$$
\frac{dn_{\text{D}}}{dt} = f_{\text{fuelling,D}}\eta_{\text{fuelling}}\Gamma_{\text{fuel}} -2 \Gamma_{\text{D+D}}- \Gamma_{\text{D+3He}} - \Gamma_{\text{D+T}} - \frac{N_{\text{T}}}{\tau_{\text{D}}^*}
$$

$$
\frac{dn_{\text{3He}}}{dt} = f_{\text{fuelling,3He}}\eta_{\text{fuelling}}\Gamma_{\text{fuel}} + \Gamma_{\text{D+D} \rightarrow \text{3He}} - \frac{N_{\text{T}}}{\tau_{\text{3He}}^*}
$$

$$
\frac{dn_{\alpha}}{dt} = \Gamma_{\text{D+3He}} + \Gamma_{\text{D+T}} - \frac{N_{\alpha}}{\tau_{\alpha}^*}
$$

In a steady state equilibrium all 4 of these equations should balance, therefore:

$$
\frac{dn_{\text{D}}}{dt} = \frac{dn_{\text{T}}}{dt} = \frac{dn_{\text{3He}}}{dt} = \frac{dn_{\alpha}}{dt} = 0
$$

Here, $\eta_{\text{fuelling}}$ is the fuelling efficiency, which quantifies how effectively fuel injected into the vacuum vessel reaches the plasma core. The value of $\eta_{\text{fuelling}}$ depends on the fuelling method. Typical values are around 0.01--0.1 for low-field-side gas puffing, 0.1--0.2 for supersonic gas injection, and 0.5--0.9 for pellet injection, which can approach unity under favourable conditions. $\Gamma_{\text{fuelling}}$ is the fuel injection rate into the vacuum vessel. The product $\eta_{\text{fuelling}} \Gamma_{\text{fuelling}}$ therefore represents the effective fuelling rate, i.e. the fraction of the injected fuel that successfully penetrates into the plasma and becomes available for fusion reactions.


The fuelling fractional compositions is given by $f$

 - $N$ is the total amount of ions in the plasma.

 - $\tau_{\text{fuel}}^*$ is the recycling corrected fuel particle confinement time given by:

    $\tau_{\text{fuel}}^* = (\tau_p) /  (1-R)$

The (effective) exhaust efficiency is is given by , $\eta_{\text{eff}} = 1- R$ 

The factor $\frac{R}{1-R}$ is the mean number of recycling events back into the burning region experieced by a particle before it is pumped away.

The definition of the recycling coefficient $R = 1- \frac{\Gamma_{\text{pumps}}}{\Gamma_{\text{out}}}$, where $\Gamma_{\text{pumps}}$ is the number of particles exhausted by the pumps per second and $\Gamma_{\text{out}}$ is the number of particles per second transported radially outwards across the separatrix.



Where $\tau_p$ is the particle confinement time which we can assume is approximately equal to the energy confinement time ($\tau_p = \tau_E$). 

!!! warning "Relation between $\tau_p$ and $\tau_E$"

    In this model the "raw" particle confinement time ($\tau_p$) is set to always match the energy confinement time ($\tau_E$). The only variation of this in terms of the "effective" particle confinement time $\tau_{p}^*$ is via the recycling coefficient ($R$). If needed new coeffcients may be added to scale ($\tau_p$) before recyling corrections if needed.

!!! note "Quantifying $R$"

    The recycling coefficient $R$, defined as the fraction of particles crossing the LCFS that return to the plasma, can depend on numerous factors—including vessel pumping speed, neutral pressure in the private‑divertor region, impurity seeding levels, and the detailed properties of the SOL. Among these parameters, $R$ is the least certain and the most difficult to quantify. In next‑step devices, the SOL temperature is expected to be high, so particles reflected from the vessel walls are mostly ionized within the SOL and are removed by pumping before they can effectively refuel the burning plasma. As a result, the recycling coefficient is anticipated to be lower than in present‑day tokamaks, where $R$ can often approach unity. An additional uncertainty is the extent of neutral penetration at the plasma edge, which influences both the pedestal density and the density profile, and therefore also affects $R$[^1].




### METIS Alpha Confinement

$$
\tau_{\alpha} = f_{\alpha}\tau_{\text{E}}\frac{R}{1-R}\tau_{\text{ne}}
$$

This is the model currently in METIS[^2] and is found in [^3]


--------------

### Tritium Flow Rate | `calculate_plasma_tritium_flow_rate()`

$$
\frac{dn_{\text{T}}}{dt} = f_{\text{fuelling,T}}\eta_{\text{fuelling}}\Gamma_{\text{fuel}} + \Gamma_{\text{D+D} \rightarrow \text{T}} - \Gamma_{\text{D+T}} - \frac{N_{\text{T}}}{\tau_{\text{T}}^*}
$$

---------------

### Deuterium Flow Rate | `calculate_plasma_deuterium_flow_rate()`

$$
\frac{dn_{\text{D}}}{dt} = f_{\text{fuelling,D}}\eta_{\text{fuelling}}\Gamma_{\text{fuel}} -2 \Gamma_{\text{D+D}}- \Gamma_{\text{D+3He}} - \Gamma_{\text{D+T}} - \frac{N_{\text{T}}}{\tau_{\text{D}}^*}
$$

---------------

### Helium-3 Flow Rate | `calculate_plasma_helium3_flow_rate()`

$$
\frac{dn_{\text{3He}}}{dt} = f_{\text{fuelling,3He}}\eta_{\text{fuelling}}\Gamma_{\text{fuel}} + \Gamma_{\text{D+D} \rightarrow \text{3He}} - \frac{N_{\text{T}}}{\tau_{\text{3He}}^*}
$$

---------------

### Alpha Particle Flow Rate | `calculate_plasma_alphas_flow_rate()`

$$
\frac{dn_{\alpha}}{dt} = \Gamma_{\text{D+3He}} + \Gamma_{\text{D+T}} - \frac{N_{\alpha}}{\tau_{\alpha}^*}
$$

-----------------

## Fuel Burnup Fration

In a steady state tokamak, the burnup fraction ($f_b$) is explicitly defined by the following rate equation:

$$
f_b = \frac{\Gamma_{\text{fusion}}}{\Gamma_{\text{fuel}}}
$$

where $\Gamma_{\text{fusion}}$ is the fusion reaction rate per second $[\text{s}^{-1}]$, and $\Gamma_{\text{fuel}}$ is the fuel injection rate into the vacuum vessel in particles per second $[\text{s}^{-1}]$.

-----------------

### Total Fuel Burnup Fraction | `calculate_fuel_burnup_fraction()`

For the total burnup fraction we state:

$$
\overbrace{f_b}^{\texttt{f_plasma_fuel_burnup}} = \frac{2\left(\Gamma_{\text{D+D}}+\Gamma_{\text{D+T}}+\Gamma_{\text{D+3He}}\right)}{\Gamma_{\text{fuel}}}
$$

Here the factor of 2 is included as each fusion reaction removes 2 particles but our fuelling rate looks at indivudal particles injected.

-------------------

### Tritium Burnup Fraction | `calculate_tritium_burnup_fraction()`

For just the tritium burnup fraction we state:

$$
\overbrace{f_b}^{\texttt{f_plasma_tritium_burnup}} = \frac{\left(\Gamma_{\text{D+T}}\right)}{\Gamma_{\text{fuel}}f_{\text{fuelling,T}}} 
$$

------------------

### Deuterium Burnup Fraction | `calculate_deuterium_burnup_fraction()`

For just the deuterium burnup fraction we state:

$$
\overbrace{f_b}^{\texttt{f_plasma_deuterium_burnup}} = \frac{2\left(\Gamma_{\text{D+D}}+\Gamma_{\text{D+3He}}\right)}{\Gamma_{\text{fuel}}f_{\text{fuelling,D}}}
$$

------------------

## Key Constraints

### Tritium Flow Consistency

This constraint can be activated by stating `icc = 93` in the input file.

This constraint ensures that the change in tritium particles as a function of time is zero. It ensures the output of `calculate_plasma_tritium_flow_rate()` is zero

**It is recommended to have this constraint on as it is a plasma consistency model**

-----------------

### Deuterium Flow Consistency

This constraint can be activated by stating `icc = 94` in the input file.

This constraint ensures that the change in deuterium particles as a function of time is zero. It ensures the output of `calculate_plasma_deuterium_flow_rate()` is zero

**It is recommended to have this constraint on as it is a plasma consistency model**

----------------

### Helium-3 Flow Consistency

This constraint can be activated by stating `icc = 95` in the input file.

This constraint ensures that the change in helium-3 particles as a function of time is zero. It ensures the output of `calculate_plasma_helium3_flow_rate()` is zero

**It is recommended to have this constraint on as it is a plasma consistency model**

----------------

### Alpha Particle Flow Consistency

This constraint can be activated by stating `icc = 96` in the input file.

This constraint ensures that the change in alpha particles as a function of time is zero. It ensures the output of `calculate_plasma_alphas_flow_rate(()` is zero

**It is recommended to have this constraint on as it is a plasma consistency model**

------------------

### Fuelling Proportion Consistency

This constraint can be activated by stating `icc = 97` in the input file.

This ensures that all 3 injected fuelling fractions sum up to 1:

$$
f_{\text{fuelling,D}} + f_{\text{fuelling,T}} + f_{\text{fuelling,3He}} = 1.0
$$

**It is recommended to have this constraint on as it is a plasma consistency model**

-----------------


[^1]: G. L. Jackson, V. S. Chan, and R. D. Stambaugh, “An Analytic Expression for the Tritium Burnup Fraction in Burning-Plasma Devices,” Fusion Science and Technology, vol. 64, no. 1, pp. 8–12, Jul. 2013, doi: https://doi.org/10.13182/fst13-a17042.

[^2]: J. F. Artaud et al., “Metis: a fast integrated tokamak modelling tool for scenario design,” Nuclear Fusion, vol. 58, no. 10, pp. 105001–105001, Aug. 2018, doi: https://doi.org/10.1088/1741-4326/aad5b1.

[^3]: D. Reiter, H. Kever, G. H. Wolf, M. Baelmans, R. Behrisch, and R. Schneider, “Helium removal from tokamaks,” Plasma Physics and Controlled Fusion, vol. 33, no. 13, pp. 1579–1600, Nov. 1991, doi: https://doi.org/10.1088/0741-3335/33/13/008.