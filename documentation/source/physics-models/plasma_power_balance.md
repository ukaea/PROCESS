# Plasma Power Balance

There are several equations that need to be satisfied in order to ensure that the plasma is in 
power equilibrium. This includes the total energy and power leaving the plasma as a whole and then also the rate of energy transfer between the ion and electron species themselves.

## Global plasma power balance

This constraint can be activated by stating `icc = 2` in the input file.

This constraint ensures self consistency between the the transport loss power used for the confinement scalings and the calculated confinement time in relation to the plasmas total thermal energy:

$$
P_{\text{L}} = \frac{W}{\tau_{\text{E}}}
$$

$$
\underbrace{\frac{3}{2}\frac{n_{\text{i}} \langle T_{\text{i}} \rangle_{\text{n}}}{\tau_{\text{i}}} + \frac{3}{2}\frac{n_{\text{e}} \langle T_{\text{e}} \rangle_{\text{n}}}{\tau_{\text{e}}}}_{\frac{W}{\tau_{\text{E}}}}  = \underbrace{\frac{f_{\alpha}P_{\alpha} + P_{\text{c}} + P_{\text{OH}} + P_{\text{HCD}}}{V_{\text{P}}} - \frac{P_{\text{rad}}}{V_{\text{p}}}}_{P_{\text{L}}}
$$

The $\frac{3}{2}n_{\text{i}} \langle T_{\text{i}} \rangle_{\text{n}}$ value is simply the volume averaged ion thermal energy density where $\langle T_{\text{i}} \rangle_{\text{n}}$ is the density weighted temperature. The same goes for the $\frac{3}{2}n_{\text{e}} \langle T_{\text{e}} \rangle_{\text{e}}$ electron thermal energy density term. $\tau_{\text{E}}$ is the confinement time calculated from the chosen confinement scaling via `i_confinement_time`. 

The constraint uses the loss power and thermal densities hence the inclusion of the $V_{\text{p}}$ plasma volume term. The constraint is adapted depending on the condition of `i_rad_loss` which governs the radiation contribution to the loss power definition, see the [radiation and energy confinement section](#effect-of-radiation-on-energy-confinement) for more info. The injected heating and current drive contribution $P_{\text{HCD}}$ is also included or excluded depending if the plasma is deemed to be ignited with the `i_plasma_ignited` switch.

**It is highly recommended to always have this constraint on as it is a global consistency checker**

----------

## Plasma Species Power Balance

While the plasma may be in global power equilibrium, the electrons and ions are continually transferring energy to each other via Coulomb collisions and certain heating and current drive systems only heat particular species. Therefore in order to be in true equilibrium the net power transfer between all species in the plasma must also be equal to zero.

!!! tip "Solving for particle power balance"

    It is highly recommended to have `f_temp_plasma_ion_electron` which represent the ion to electron temperature ratio $\frac{T_i}{T_e}$, as an iteration variable. This allows the solver to have another degree of freedom upon which to change the fusion power and the equilibration power.

-------------------

## Electron Power Balance

This constraint can be activated by stating `icc = 3` in the input file.

For the electrons this internal power blance must be satisifed:

$$
\underbrace{P_{\Omega} + P_{\text{HCD,e}}+ (f_{e}\times(f_{\alpha}P_{\alpha}))}_{\text{Power gain}} = \underbrace{P_{\text{rad}} - P_{\text{ei}} - \frac{3}{2}\frac{n_{\text{e}} \langle T_{\text{e}} \rangle_{\text{n}}}{\tau_{\text{e}}}}_{\text{Power loss}}
$$

where $P_{\Omega}$ is the plasma ohmic heating power, $P_{\text{HCD,e}}$ is the external heating and current drive power that only goes to electrons,$(f_{e}\times(f_{\alpha}P_{\alpha}))$ is the fraction of the coupled alpha particle power that goes to the electrons, $P_{\text{rad}}$ is the total radiation power given off by the plasma, $P_{\text{ei}}$ is the [electron-ion equilibration power](#ion-electron-equilibration-power-density--calculate_ion_electron_equilibration_power).

**It is highly recommended to always have this constraint on as it is a plasma equilibrium checker**

------------

## Ion Power Balance

This constraint can be activated by stating `icc = 4` in the input file.

For the ions this internal power blance must be satisifed:

$$
\underbrace{P_{\text{HCD,i}} + (f_{i}\times(f_{\alpha}P_{\alpha})) + P_{\text{ei}}}_{\text{Power gain}} = \underbrace{\frac{3}{2}\frac{n_{\text{i}} \langle T_{\text{i}} \rangle_{\text{n}}}{\tau_{\text{i}}}}_{\text{Power Loss}}
$$

where $P_{\text{HCD,i}}$ is the external heating and current drive power that only goes to ions,$(f_{e}\times(f_{\alpha}P_{\alpha}))$ is the fraction of the coupled alpha particle power that goes to the ions, $P_{\text{ei}}$ is the [electron-ion equilibration power](#ion-electron-equilibration-power-density--calculate_ion_electron_equilibration_power).

**It is highly recommended to always have this constraint on as it is a plasma equilibrium checker**

-----------

## Ion-electron equilibration power density | `calculate_ion_electron_equilibration_power()`

The equilibration power $P_{\text{ei}}$ refers to the rate of energy transfer between plasma species—typically ions and electrons—due to elastic Coulomb collisions. It is driven by the temperature difference between the species. The equation below represents the rate of energy transfer from the electrons to the ions due to these elastic Coulomb collisions caused by the temperature difference between the electrons and the ions:

$$
P_{\text{ei}} = \frac{3}{2} \frac{n_e \text{e} (T_e - T_i)}{\tau_{\text{eq}}}
$$

where $n_e$ is the electron density, $T_i$ and $T_e$ is the ion and electron temperatures in $[\text{eV}]$ and $\tau_{\text{eq}}$ is the ion-electron equilibration time.