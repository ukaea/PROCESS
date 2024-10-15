# Fusion Reactions

## Overview

The most likely fusion reaction to be utilised in a power plant is the
deuterium-tritium reaction:

$$
\mathrm{D + T} \Longrightarrow \mathrm{^{4}He + n + 17.6 \,MeV}
$$

Roughly 20% of the energy produced is given to the alpha particles (\(^4\)He). The remaining 80% is carried
away by the neutrons, which deposit their energy within the blanket and shield and other reactor components.
The fraction of the alpha energy deposited in the plasma is `falpha`.

PROCESS can also model D-\(^3\)He power plants, which utilise the following
primary fusion reaction:

$$
\mathrm{D + \text{$^3$He}} \Longrightarrow \mathrm{^{4}He + p + 18.3 \,MeV}
$$

The fusion reaction rate is significantly different to that for D-T fusion,
and the power flow from the plasma is modified since charged particles are
produced rather than neutrons. Because only charged particles (which remain in
the plasma) are produced by this reaction, the whole of the fusion power is
used to heat the plasma. Useful energy is extracted from the plasma since the
radiation power produced is very high, and this, in theory, can be converted to
electricity without using a thermal cycle.

Since the temperature required to ignite the D-\(^3\)He reaction is considerably
higher than that for D-T, it is necessary to take into account the following
D-D reactions, which have significant reaction rates at such temperatures:

$$\begin{aligned}
\mathrm{D + D}  & \Longrightarrow \mathrm{^{3}He + n + 3.27 \,MeV} \\
\mathrm{D + D}  & \Longrightarrow \mathrm{T + p + 4.03 \,MeV}
\end{aligned}$$

Also, as tritium is produced by the latter reaction, D-T fusion also occurs. 
As a result, there is still a small amount of neutron power
extracted from the plasma.

Pure D-\(^3\)He tokamak power plants do not include breeding blankets, because 
no tritium needs to be produced for fuel.

The contributions from all four of the above fusion reactions are included in
the total fusion power production calculation. The fusion reaction rates are
calculated using the parameterizations in [^1], integrated over the plasma 
profiles (correctly, with or without pedestals).

The fractional composition of the 'fuel' ions (D, T and \(^3\)He) is
controlled using the three variables `fdeut`, `ftrit` and `fhe3`, respectively:

$$\begin{aligned}
n_{\mbox{fuel}}  & = n_D + n_T + n_{\mathrm{^{3}He}}  \;\;\; \mbox{particles/m$^3$} \\
n_D  & = \mathtt{fdeut} \, n_{\mbox{fuel}} \\
n_T  & = \mathtt{ftrit} \, n_{\mbox{fuel}} \\
n_{\mathrm{^{3}He}} & = \mathtt{fhe3} \, n_{\mbox{fuel}}
\end{aligned}$$

PROCESS checks that $\mathtt{fdeut} + \mathtt{ftrit} + \mathtt{fhe3} = 1.0$, and stops with an error message otherwise. Shouldnt this account for impurities also?

-------------------------

## Fusion Reaction Class | `FusionReactionRate`

### Initialization | `__init__()`

Initialize the FusionReactionRate class with the given plasma profile.

#### Parameters:
- `plasma_profile (PlasmaProfile)`: The parameterized temperature and density profiles of the plasma. Taken from the plasma_profile object.

#### Attributes:
- `plasma_profile (PlasmaProfile)`: The parameterized temperature and density profiles of the plasma.
- `sigmav_dt_average (float)`: Average fusion reaction rate $<\sigma v>$ for D-T.
- `dhe3_power_density (float)`: Fusion power density produced by the D-3He reaction.
- `dd_power_density (float)`: Fusion power density produced by the D-D reactions.
- `dt_power_density (float)`: Fusion power density produced by the D-T reaction.
- `alpha_power_density (float)`: Power density of alpha particles produced.
- `charged_power_density (float)`: Power density of charged particles produced.
- `neutron_power_density (float)`: Power density of neutrons produced.
- `fusion_rate_density (float)`: Fusion reaction rate density.
- `alpha_rate_density (float)`: Alpha particle production rate density.
- `proton_rate_density (float)`: Proton production rate density.

All variables above are initialized to be 0.0.

---------------------

### Deuterium Branching fraction | `deuterium_branching()`

Calculates the relative rate of tritium producing D-D reactions to 3He ones based on the volume averaged ion temperature[^1].
Valid for ion temperatures between 0.5 keV and 200 keV.
The deviation of the fit from the R-matrix branching ratio is always smaller than 0.5%.

$$
 = 1.02934 - 8.3264\times 10^{-3}\langle T_{\text{e}} \rangle + 1.7631\times 10^{-4}\langle T_{\text{e}} \rangle^2 \\
-1.8201\times 10^{-6}\langle T_{\text{e}} \rangle^3 + 6.9855\times 10^{-9}\langle T_{\text{e}} \rangle^4
$$

-----------------------

### Calculate fusion reactions

There are 4 key functions for calculating the fusion reaction for the plasma. They are `dt_reaction()`, `dhe3_reaction()`, `dd_helion_reaction()` and `dd_triton_reaction()`. They all perform the same key calculations below but with their own specific values for their reactions.

#### Detailed Steps
1. **Initialize Bosch-Hale Constants**: Initializes the Bosch-Hale constants for the required reaction using predefined reaction constants stored in the BoschHaleConstants dataclass.
2. **Calculate Fusion Reaction Rate**: Uses Simpson's rule to integrate the fusion reaction rate over the plasma profile.
3. **Calculate Fusion Power Density**: Compute the fusion power density produced by the given reaction. Using the reaction energy calculated and stored in `constants.f90`. The reactant density is is given by $\mathtt{fdeut, ftrit}$ or $\mathtt{fhe3}$ multiplied by the volume averaged ion density.
4. **Calculate Fusion Power Densities**: Compute the fusion power density for alpha particles, neutrons and other charged particles, depending on the reaction. Energy branching fractions used are calculated and called from `constants.f90`
5. **Calculate Fusion Rate Densities**: Compute the fusion rate density and for alpha particles, neutrons and other charged particles, depending on the reaction.
6. **Update Reaction Power Density**: Updates the object attribute for the specific reaction power density.
7. **Sum Fusion Rates**: Call the [`sum_fusion_rates()`](#sum-the-fusion-rates--sum_fusion_rates) function to add the reaction to the global plasma power balance.

#### Attributes Updated
The method updates the following attributes:

- `self.sigmav_dt_average`: Average fusion reaction rate `<sigma v>` for D-T.
- `self.dt_power_density`: Fusion power density produced by the D-T reaction.
- `self.alpha_power_density`: Power density of alpha particles produced.
- `self.charged_power_density`: Power density of charged particles produced.
- `self.neutron_power_density`: Power density of neutrons produced.
- `self.fusion_rate_density`: Fusion reaction rate density.
- `self.alpha_rate_density`: Alpha particle production rate density.
- `self.proton_rate_density`: Proton production rate density.

-----------------------

### Sum the fusion rates | `sum_fusion_rates()`

This method updates the cumulative class fusion power densities and reaction rates for alpha particles, charged particles, neutrons, and protons.

#### Parameters:
- `alpha_power_add` (float): Alpha particle fusion power per unit volume [MW/m³].
- `charged_power_add` (float): Other charged particle fusion power per unit volume [MW/m³].
- `neutron_power_add` (float): Neutron fusion power per unit volume [MW/m³].
- `fusion_rate_add` (float): Fusion reaction rate per unit volume [reactions/m³/s].
- `alpha_rate_add` (float): Alpha particle production rate per unit volume [/m³/s].
- `proton_rate_add` (float): Proton production rate per unit volume [/m³/s].

The above input values are added to the current object values for their respective power density.

-----------------------

### Calculate fusion rates | `calculate_fusion_rates()`

This runner function is called to run all 4 fusion reaction calculation functions.

This method sequentially calculates the fusion reaction rates and power densities for the following reactions:
- Deuterium-Tritium (D-T)
- Deuterium-Helium-3 (D-3He)
- Deuterium-Deuterium (D-D) first branch
- Deuterium-Deuterium (D-D) second branch

It updates the instance attributes for the cumulative power densities and reaction rates for alpha particles, charged particles, neutrons, and protons.

-----------------------

### Set global physics variables | `set_physics_variables()`

This method sets the required physics variables in the `physics_variables` and `physics_module` modules. It updates the global physics variables and module variables with the current instance's fusion power densities and reaction rates.

#### Updates:
- `physics_variables.alpha_power_density`: Updated with `self.alpha_power_density`
- `physics_variables.charged_power_density`: Updated with `self.charged_power_density`
- `physics_variables.neutron_power_density`: Updated with `self.neutron_power_density`
- `physics_variables.fusion_rate_density`: Updated with `self.fusion_rate_density`
- `physics_variables.alpha_rate_density`: Updated with `self.alpha_rate_density`
- `physics_variables.proton_rate_density`: Updated with `self.proton_rate_density`
- `physics_module.sigmav_dt_average`: Updated with `self.sigmav_dt_average`
- `physics_module.dt_power_density`: Updated with `self.dt_power_density`
- `physics_module.dhe3_power_density`: Updated with `self.dhe3_power_density`
- `physics_module.dd_power_density`: Updated with `self.dd_power_density`

-----------------------

## Key Constraints

### Fusion Power Upper limit

This constraint can be activated by stating `icc = 9` in the input file.

The value of `powfmax` can be set to the desired maximum fusion power. The scaling value `ffuspow` can be varied also.

---------------------------------

### Global power balance for electrons

This constraint can be activated by stating `icc = 41` in the input file.

-------------------------------

### Global power balance for ions

This constraint can be activated by stating `icc = 41` in the input file.

-----------------------------

### Q value lower limit

This constraint can be activated by stating `icc = 28` in the input file.

The value of `bigqmin` can be set to the minimum desired $Q_{\text{plasma}}$ value. The scaling value `fqval` can be varied also.

-----------------------

[^1]: H.-S. Bosch and G. M. Hale, “Improved formulas for fusion cross-sections and thermal reactivities,” Nuclear Fusion, vol. 32, no. 4, pp. 611–631, Apr. 1992, doi: https://doi.org/10.1088/0029-5515/32/4/i07.