# Plasma Fuelling

The control of fuelling is governed by 4 key particle flux equations for each of the primary fuel species and the helium ash, $\alpha$.

$$
\frac{dn_{\text{T}}}{dt} = \eta_{\text{fuelling}}\Gamma_{\text{T}} + \Gamma_{\text{D+D} \rightarrow \text{T}} - \Gamma_{\text{D+T}} - \frac{N_{\text{T}}}{\tau_{\text{T}}^*}
$$

$$
\frac{dn_{\text{D}}}{dt} = \eta_{\text{fuelling}}\Gamma_{\text{T}} -2 \Gamma_{\text{D+D}}- \Gamma_{\text{D+3He}} - \Gamma_{\text{D+T}} - \frac{N_{\text{T}}}{\tau_{\text{D}}^*}
$$

$$
\frac{dn_{\text{3He}}}{dt} = \eta_{\text{fuelling}}\Gamma_{\text{3He}} + \Gamma_{\text{D+D} \rightarrow \text{3He}} - \frac{N_{\text{T}}}{\tau_{\text{3He}}^*}
$$

$$
\frac{dn_{\alpha}}{dt} = \Gamma_{\text{D+3He}} + \Gamma_{\text{D+T}} - \frac{N_{\alpha}}{\tau_{\alpha}^*}
$$

In a steady state equilibrium all 4 of these equations should balance, therefore:

$$
\frac{dn_{\text{D}}}{dt} = \frac{dn_{\text{T}}}{dt} = \frac{dn_{\text{3He}}}{dt} = \frac{dn_{\alpha}}{dt} = 0
$$

Here $\eta_{\text{fuelling}}$ is the fuelling efficiecny which represents the method of injecting fuel into the plasma. Gas puffing on the low field side is probably around 0.01-0.1, supersonic gas is 0.1 and 0.2 and using pellets can get you close to unity with 0.5-0.9. $\Gamma_{\text{fuelling}}$ is the fuel injection rate into the vacuum vessel, so $\eta_{\text{fuelling}} \Gamma_{\text{fuelling}}$ together presents the fraction of injected fuel that actually makes it into the plasma core to fuse. 

 - $N$ is the total amount of ions in the plasma.

 - $\tau_{\text{fuel}}^*$ is the recycling corrected fuel particle confinement time given by:

    $\tau_{\text{fuel}}^* = (\tau_p) /  (1-R)$


Where $\tau_p$ is the particle confinement time which we can assume is approximately equal to the energy confinement time ($\tau_p = \tau_E$). 

The recycling coefficient $R$, defined as the fraction of particles crossing the LCFS that return to the plasma, can depend on numerous factors—including vessel pumping speed, neutral pressure in the private‑divertor region, impurity seeding levels, and the detailed properties of the SOL. Among these parameters, $R$ is the least certain and the most difficult to quantify. In next‑step devices, the SOL temperature is expected to be high, so particles reflected from the vessel walls are mostly ionized within the SOL and are removed by pumping before they can effectively refuel the burning plasma. As a result, the recycling coefficient is anticipated to be lower than in present‑day tokamaks, where $R$ can often approach unity. An additional uncertainty is the extent of neutral penetration at the plasma edge, which influences both the pedestal density and the density profile, and therefore also affects $R$[^1].




[^1]: G. L. Jackson, V. S. Chan, and R. D. Stambaugh, “An Analytic Expression for the Tritium Burnup Fraction in Burning-Plasma Devices,” Fusion Science and Technology, vol. 64, no. 1, pp. 8–12, Jul. 2013, doi: https://doi.org/10.13182/fst13-a17042.
‌