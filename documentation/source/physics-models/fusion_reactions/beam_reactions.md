# Beam Fusion Model

This page describes the neutral beam fusion model currently implemented in `beam_fusion()`.

The model is a reduced beam-target treatment for neutral beam ions injected into the plasma. `PROCESS` first determines an effective beam current that already accounts for beam transport effects such as shine-through losses. The beam-fusion model then uses this net beam current to build a steady population of fast ions and estimates how many of those ions fuse with the background plasma. `beam_fusion()` does not explicitly model beam attenuation, orbit effects, or spatial beam evolution, and instead treats the fast-ion population in a volume-averaged (0D) sense. The model is calculated in the following steps:

1. computes a beam slowing-down time and a critical energy for deuterium and tritium beam ions,
2. estimates the steady-state hot beam ion densities from the beam source rate and slowing-down residence time,
3. evaluates a fast-ion pressure and a pressure-equivalent deposited beam energy,
4. computes effective beam-target fusion rate coefficients for hot deuterium and hot tritium beam ions,
5. forms total beam-target reaction rates and converts them to alpha power,
6. computes the neutral beam beta contribution from the hot beam population.

The main functions involved are:

- `beam_fusion()`  
Top-level routine. Computes beam slowing-down properties, beam-target fusion alpha power, and neutral beam beta.

- `beam_slowing_down_state()`  
Splits the beam into deuterium and tritium components, computes steady-state hot beam densities, critical speeds, fast-ion pressures, and the density-weighted deposited beam energy.

- `fast_ion_pressure_integral()`  
Returns the closed-form dimensionless integral factor used in the fast-ion pressure expression.

- `beam_reaction_rate_coefficient()`  
Evaluates an effective beam-target fusion rate coefficient for a slowing-down fast-ion population.

- `_hot_beam_fusion_reaction_rate_integrand()`  
Internal integrand used in the beam-target rate coefficient calculation.

- `_beam_fusion_cross_section()`  
Internal beam fusion cross-section fit used by the beam-target rate coefficient calculation.

- `beam_target_reaction_rate()`  
Converts beam density, target density, and beam-target reactivity into a total reaction rate.

- `alpha_power_beam()`  
Converts the total beam-target reaction rate into alpha power.

Due to the small contribution of fusion power from the neutral beams, only D-T beam-target reactions are included. D-D beam contributions are neglected in this model.

The beam fusion calculations will only run if the calculated beam current is greater than 0, i.e. when an NBI heating and current drive configuration is active. Currently, the model only runs when the plasma is not ignited `i_plasma_ignited = 0`.

The NBI parameters taken from the current drive and plasma state are the total beam current (`c_beam_total`), the beam energy (`e_beam_kev`), the tritium fraction in the beam (`f_beam_tritium`), the bulk deuterium and tritium fractions (`f_deuterium_plasma`, `f_tritium_plasma`), and the plasma state quantities required for slowing down.

Please see the [H&CD section](../../eng-models/heating_and_current_drive/heating-and-current-drive.md) of the docs for more info.

---

## Beam slowing down properties | `beam_fusion()`

This section explains the stages of calculation in the function `beam_fusion()`

### Calculate the beam ion slowing down time

The beam slowing down time used in the model is implemented as[^deng1987]:

$$
\tau_{\text{slow}} =
1.99\times10^{19}
\left[
\text{A}_{\text{D}}\left(1-\text{f}_{\text{beam,T}}\right)
+
\text{A}_{\text{T}}\text{f}_{\text{beam,T}}
\right]
\frac{\langle \text{T}_{\text{e}}\rangle^{3/2}}
{\langle \text{n}_{\text{e}}\rangle \ln \Lambda_{\text{ie}}}
$$

where:

- $A_{\text{D}}$ is the deuteron mass in amu,
- $A_{\text{T}}$ is the triton mass in amu,
- $f_{\text{beam,T}}$ is the tritium fraction in the injected beam,
- $\langle T_{\text{e}}\rangle$ is the density-weighted electron temperature in keV,
- $\langle n_{\text{e}}\rangle$ is the volume-averaged electron density in $\text{m}^{-3}$
- $\ln\Lambda_{\text{ie}}$ is the ion-electron Coulomb logarithm.
- The numerical prefactor $1.99\times10^{19}$ arises from rewriting the classical slowing-down time expression in practical units.

---

### Calculate the beam critical energy

For an energetic ion slowing down in a plasma, there is a critical energy $E_{\text{crit}}$ at which the rate of energy loss to ions equals the rate of energy loss to electrons. Above this energy electron drag dominates, while below it ion drag dominates.

In the current implementation, the deuterium critical energy [^sheffield1994] is

$$
\text{E}_{\text{crit,D}} =
14.8
\text{A}_{\text{D}}
\text{T}_{\text{e}}
\text{Z}_{\text{eff,mw}}^{2/3}
\frac{\ln\Lambda_{\text{ie}}+4.0}{\ln\Lambda_{\text{ie}}}
\qquad [\text{keV}]
$$

where $Z_{\text{eff,mw}}$ is the mass-weighted effective plasma charge.

The tritium critical energy is then scaled by the beam ion mass ratio:

$$
\text{E}_{\text{crit,T}} =
\text{E}_{\text{crit,D}}
\left(\frac{\text{A}_{\text{T}}}{\text{A}_{\text{D}}}\right)
$$

#### Derivation of beam slowing down rate and critical energy

The rate of slowing down of a test particle of mass $M$, charge $Ze$ and energy $E$, due to Coulomb collisions with a background species of mass $m_j$, charge $Z_j e$, density $n_j$ and temperature $T_j$, is given by[^deng1987]

$$
\frac{\text{dE}}{\text{dt}}
=
\left[
-\Phi(\text{x}_{\text{j}})
+
\text{x}_{\text{j}}\left(1+\frac{\text{m}_{\text{j}}}{\text{M}}\,\Phi^{\prime}(\text{x}_{\text{j}})\right)
\right]
\frac{4\pi \text{n}_{\text{j}}}{\text{m}_{\text{j}} \text{V}}
\left(\frac{\text{Z} \text{Z}_{\text{j}} \text{e}^2}{4\pi\varepsilon_0}\right)^2
\ln \Lambda_{\text{j}}
$$

where

$$
\Phi^{\prime}(\text{x}) = \frac{\text{d}\Phi}{\text{dx}}
$$

$$
\text{V} = \sqrt{\frac{2\text{E}}{\text{M}}}, \qquad
\text{V}_{\text{j}} = \sqrt{\frac{2\text{kT}_{\text{j}}}{\text{m}_{\text{j}}}}, \qquad
\text{x}_{\text{j}} = \frac{\text{V}}{\text{V}_{\text{j}}}
$$

For fast ions in fusion plasmas, the ion contribution and electron contribution may be approximated separately, leading to the standard slowing-down form[^wesson2011]:

$$
\frac{\mathrm{d}\text{E}}{\mathrm{d}\text{t}}
=
-\frac{\text{A} \text{Z}^2 \sqrt{\text{M}}}{\sqrt{\text{E}}}
-
\frac{\text{B} \text{Z}^2 \text{E}}{\text{M}}
$$

with coefficients

$$
\begin{aligned}
\text{A} &=
\frac{4\pi}{\sqrt{2}}
\left(\frac{\text{e}^2}{4\pi\varepsilon_0}\right)^2
\sum_{\text{j}}
\left(
\frac{\text{n}_{\text{j}} \text{Z}_{\text{j}}^2}{\text{m}_{\text{j}}}\ln\Lambda_{\text{j}}
\right) \\
\text{B} &=
\frac{16\sqrt{\pi}}{3\text{kT}_{\text{e}}}
\sqrt{\frac{\text{m}_{\text{e}}}{2\text{kT}_{\text{e}}}}
\left(\frac{\text{e}^2}{4\pi\varepsilon_0}\right)^2
\text{n}_{\text{e}} \ln\Lambda_{\text{e}}
\end{aligned}
$$

This can be rewritten in the form

$$
\frac{\mathrm{d}\text{E}}{\mathrm{d}\text{t}}
=
-\frac{2\text{E}}{\tau_{\text{slow}}}
\left[
1+\left(\frac{\text{E}_{\text{c}}}{\text{E}}\right)^{3/2}
\right]
$$

where

$$
\text{E}_{\text{c}}
=
\left[
\frac{3\sqrt{\pi}}{4}
\frac{\text{M}^{3/2}}{\text{n}_{\text{e}}\sqrt{\text{m}_{\text{e}}}}
\sum_{\text{j}}
\left(
\frac{\text{n}_{\text{j}} \text{Z}_{\text{j}}^2}{\text{m}_{\text{j}}}\ln\Lambda_{\text{j}}
\right)
\frac{1}{\ln\Lambda_{\text{e}}}
\right]^{2/3}
\text{kT}_{\text{e}}
$$

and

$$
\tau_{\text{slow}}
=
\frac{3(\text{kT}_{\text{e}})^{3/2}}{4\sqrt{2\pi \text{m}_{\text{e}}} \text{Z}^2}
\left(\frac{4\pi\varepsilon_0}{\text{e}^2}\right)^2
\frac{\text{M}}{\text{n}_{\text{e}} \ln\Lambda_{\text{e}}}
$$

In this regime, $\tau_{\text{slow}}$ is the characteristic electron-drag slowing-down timescale.

$\blacksquare$

#### Set the plasma deuterium and tritium ion densities

The bulk target ion densities used in the beam-target reactions are

$$
\text{n}_{\text{D,plasma}}
=
\text{n}_{\text{fuel}}
\text{f}_{\text{D,plasma}}
$$

$$
\text{n}_{\text{T,plasma}}
=
\text{n}_{\text{fuel}}
\text{f}_{\text{T,plasma}}
$$

where $n_{\text{fuel}}$ is the volume-averaged total fuel ion density.

---

### Calculate the beam alpha powers, beam densities and deposited energy

[`beam_slowing_down_state()`](#neutral-beam-alpha-power-beam-densities-and-deposited-energy--beam_slowing_down_state) is run to calculate the hot beam densities, critical speeds, and deposited beam energy.

[`beam_reaction_rate_coefficient()`](#beam-fusion-reaction-rate-coefficient--beam_reaction_rate_coefficient) is then run for the deuterium and tritium beam components to obtain effective beam-target rate coefficients.

[`beam_target_reaction_rate()`](#beam-target-fusion-reaction-rate--beam_target_reaction_rate) and [`alpha_power_beam()`](#beam-fusion-alpha-power--alpha_power_beam) are used to convert these into alpha power.

### Set the returned alpha power

The total neutral beam alpha power is

$$
\text{P}_{\alpha,\text{beam}}
=
\mathtt{beamfus0}
\left(
\text{P}_{\alpha,\text{D-beam}}
+
\text{P}_{\alpha,\text{T-beam}}
\right)
$$

### Calculate the neutral beam beta

The neutral beam beta is computed from the hot beam density and the deposited beam energy:

$$
\beta_{\text{beam}}
=
\mathtt{betbm0}
\times
4.03\times10^{-22}
\times
\frac{2}{3}
\frac{\text{n}_{\text{beam,hot}} \text{E}_{\text{beam,deposited}}}
{\text{B}_{\phi}^2 + \text{B}_{\theta}^2}
$$

where:

- $n_{\text{beam,hot}}$ is the total hot beam ion density,
- $E_{\text{beam,deposited}}$ is the density-weighted deposited beam ion energy in keV,
- $B_{\phi}$ is the toroidal magnetic field on axis,
- $B_{\theta}$ is the poloidal magnetic field.

The value of $E_{\text{beam,deposited}}$ is the pressure-equivalent deposited energy of the hot beam ions, **not** the initial beam injection energy.

---

## Neutral beam alpha power, beam densities and deposited energy | `beam_slowing_down_state()`

### Calculate the beam current fractions

The beam current is split into deuterium and tritium components:

$$
\text{I}_{\text{beam,D}}
=
\text{I}_{\text{beam}}
\left(1-\text{f}_{\text{beam,T}}\right)
$$

$$
\text{I}_{\text{beam,T}}
=
\text{I}_{\text{beam}}
\text{f}_{\text{beam,T}}
$$

### Calculate the characteristic slowing-down time to the thermal range

Using the classical slowing-down model, the characteristic time for the beam ion energy to slow from the birth energy to the thermal range is implemented as:

$$
\tau_{\text{slow,D}}^{*}
=
\frac{\tau_{\text{slow}}}{3}
\ln\left[
1+
\left(
\frac{\text{E}_{\text{beam}}}{\text{E}_{\text{crit,D}}}
\right)^{3/2}
\right]
$$

$$
\tau_{\text{slow,T}}^{*}
=
\frac{\tau_{\text{slow}}}{3}
\ln\left[
1+
\left(
\frac{\text{E}_{\text{beam}}}{\text{E}_{\text{crit,T}}}
\right)^{3/2}
\right]
$$

### Set the fast beam ion densities

The steady-state hot beam ion densities are set from source rate times residence time:

$$
\langle \text{n}_{\text{beam}} \rangle_{\text{D}}
=
\frac{\text{I}_{\text{beam,D}}\tau_{\text{slow,D}}^{*}}
{\text{e} \text{V}_{\text{plasma}}}
$$

$$
\langle \text{n}_{\text{beam}} \rangle_{\text{T}}
=
\frac{\text{I}_{\text{beam,T}}\tau_{\text{slow,T}}^{*}}
{\text{e} \text{V}_{\text{plasma}}}
$$

The total hot beam ion density is then

$$
\langle \text{n}_{\text{beam}} \rangle_{\text{hot}}
=
\langle \text{n}_{\text{beam}} \rangle_{\text{D}}
+
\langle \text{n}_{\text{beam}} \rangle_{\text{T}}
$$

### Calculate the speeds of ions at the critical energy

Assuming non-relativistic energies, the beam ion speeds at the critical energy are

$$
\text{v}_{\text{crit,D}}
=
\sqrt{
\frac{2 \text{e}_{\text{keV}} \text{E}_{\text{crit,D}}}
{\text{m}_{\text{u}} \text{A}_{\text{D}}}
}
$$

$$
\text{v}_{\text{crit,T}}
=
\sqrt{
\frac{2 \text{e}_{\text{keV}} \text{E}_{\text{crit,T}}}
{\text{m}_{\text{u}} \text{A}_{\text{T}}}
}
$$

where:

- $e_{\text{keV}}$ is the conversion from keV to joules,
- $m_u$ is the atomic mass unit.

### Calculate the fast ion pressures

First define the source rates per unit volume:

$$
\text{S}_{\text{D}}
=
\frac{\text{I}_{\text{beam,D}}}{\text{e} \text{V}_{\text{plasma}}}
$$

$$
\text{S}_{\text{T}}
=
\frac{\text{I}_{\text{beam,T}}}{\text{e} \text{V}_{\text{plasma}}}
$$

The implemented pressure coefficients are then

$$
\text{C}_{\text{p},\text{D}}
=
\frac{
\text{A}_{\text{D}} \text{m}_{\text{u}} \tau_{\text{slow}} \text{v}_{\text{crit,D}}^2 \text{S}_{\text{D}}
}
{3 \text{e}_{\text{keV}}}
$$

$$
\text{C}_{\text{p},\text{T}}
=
\frac{
\text{A}_{\text{T}} \text{m}_{\text{u}} \tau_{\text{slow}} \text{v}_{\text{crit,T}}^2 \text{S}_{\text{T}}
}
{3 \text{e}_{\text{keV}}}
$$

The fast ion pressures are

$$
\text{p}_{\text{D}}
=
\text{C}_{\text{p},\text{D}}
\times
\underbrace{
\left[
\frac{\text{x}_{\text{c}}^2}{2}
+
\frac{1}{6}\ln\left(\frac{\text{x}_{\text{c}}^2+2\text{x}_{\text{c}}+1}{\text{x}_{\text{c}}^2-\text{x}_{\text{c}}+1}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{2\text{x}_{\text{c}}-1}{\sqrt{3}}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{1}{\sqrt{3}}\right)
\right]
}_{\mathtt{fast\_ion\_pressure\_integral}(\text{E}_{\text{beam}},\text{E}_{\text{crit,D}})}
$$

$$
\text{p}_{\text{T}}
=
\text{C}_{\text{p},\text{T}}
\times
\underbrace{
\left[
\frac{\text{x}_{\text{c}}^2}{2}
+
\frac{1}{6}\ln\left(\frac{\text{x}_{\text{c}}^2+2\text{x}_{\text{c}}+1}{\text{x}_{\text{c}}^2-\text{x}_{\text{c}}+1}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{2\text{x}_{\text{c}}-1}{\sqrt{3}}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{1}{\sqrt{3}}\right)
\right]
}_{\mathtt{fast\_ion\_pressure\_integral}(\text{E}_{\text{beam}},\text{E}_{\text{crit,T}})}
$$

with

$$
\text{x}_{\text{c}} = \sqrt{\frac{\text{E}_{\text{beam}}}{\text{E}_{\text{crit}}}}
$$

### Beam fast ion pressure integral | `fast_ion_pressure_integral()`

This internal function returns the dimensionless pressure integral factor used in the fast-ion pressure expression. It takes the beam birth energy $E_{\text{beam}}$ and the critical energy $E_{\text{crit}}$ for the required beam ion species.

#### Derivation

The fast ions, because of their finite slowing down time, develop a finite pressure which must be supported by the magnetic field.

To calculate this pressure, introduce the distribution function $g(E)$ of fast ions, where $g(E)$ is the number of ions per unit energy per unit volume, satisfying the steady-state kinetic equation

$$
\frac{\partial}{\partial \text{E}}
\left(
\text{g} \frac{\mathrm{d}\text{E}}{\mathrm{d}\text{t}}
\right)
=
\text{S}(\text{E})
$$

with slowing-down law

$$
\frac{\mathrm{d}\text{E}}{\mathrm{d}\text{t}}
=
-\frac{2\text{E}}{\tau_{\text{s}}}
\left[
1+\left(\frac{\text{E}_{\text{c}}}{\text{E}}\right)^{3/2}
\right]
$$

If the ions are born monoenergetically,

$$
\text{S}(\text{E})=\text{S}_0\delta(\text{E}-\text{E}_0)
$$

and with the boundary condition $g(E)=0$ for $E>E_0$, the steady-state distribution becomes

$$
\text{g}(\text{E})=
\begin{cases}
\dfrac{\text{S}_0\tau_{\text{s}}}{2\text{E}\left[1+\left(\text{E}_{\text{c}}/\text{E}\right)^{3/2}\right]}, & \text{E}<\text{E}_0 \\
0, & \text{E}>\text{E}_0
\end{cases}
$$

The fast ion pressure is then

$$
\text{p}
=
\frac{2}{3}
\int_0^{\text{E}_0} \text{g}(\text{E}) \text{E} \, \mathrm{d}\text{E}
$$

which can be written as

$$
\text{p}
=
\frac{\text{M} \text{S}_0 \tau_{\text{s}} \text{V}_{\text{c}}^2}{3}
\int_0^{\text{x}_{\text{c}}}
\frac{\text{x}^4}{1+\text{x}^3}\,\mathrm{d}\text{x}
$$

and evaluated analytically as

$$
\text{p}
=
\frac{\text{M} \text{S}_0 \tau_{\text{s}} \text{V}_{\text{c}}^2}{3}
\left[
\frac{\text{x}_{\text{c}}^2}{2}
+
\frac{1}{6}\ln\left(\frac{\text{x}_{\text{c}}^2+2\text{x}_{\text{c}}+1}{\text{x}_{\text{c}}^2-\text{x}_{\text{c}}+1}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{2\text{x}_{\text{c}}-1}{\sqrt{3}}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{1}{\sqrt{3}}\right)
\right]
$$

`fast_ion_pressure_integral()` returns the bracketed dimensionless factor only.

$\blacksquare$

### Calculate the deposited fast ion energy from the pressure

The deposited beam ion energy is inferred from the pressure relation

$$
\text{p} = \frac{1}{3}\text{nmv}_{\text{rms}}^2 = \frac{2}{3}\text{n}\langle \text{E}\rangle
$$

so that the pressure-equivalent deposited energies are

$$
\text{E}_{\text{hot,D}}
=
\frac{3}{2}
\frac{\text{p}_{\text{D}}}{\langle \text{n}_{\text{beam}} \rangle_{\text{D}}}
$$

$$
\text{E}_{\text{hot,T}}
=
\frac{3}{2}
\frac{\text{p}_{\text{T}}}{\langle \text{n}_{\text{beam}} \rangle_{\text{T}}}
$$

with the convention that the deposited energy is set to zero if the corresponding hot beam density is zero.

The total deposited beam ion energy returned to `beam_fusion()` is the density-weighted average:

$$
\text{E}_{\text{hot,total}}
=
\frac{
\text{E}_{\text{hot,D}} \langle \text{n}_{\text{beam}} \rangle_{\text{D}}
+
\text{E}_{\text{hot,T}} \langle \text{n}_{\text{beam}} \rangle_{\text{T}}
}
{\langle \text{n}_{\text{beam}} \rangle_{\text{hot}}}
$$

with the convention that this is set to zero if $\langle n_{\text{beam}} \rangle_{\text{hot}}=0$.

1. **Return the slowing-down state**

`beam_slowing_down_state()` returns:

- deuterium beam ion density,
- tritium beam ion density,
- deuterium critical speed,
- tritium critical speed,
- total hot beam ion density,
- density-weighted deposited beam energy.

---

## Beam fusion reaction rate coefficient | `beam_reaction_rate_coefficient()`

### Calculate the beam velocity

The beam velocity is calculated from the input beam energy and beam ion mass:

$$
\text{v}_{\text{beam}}
=
\sqrt{
\frac{2 \text{e}_{\text{keV}} \text{E}_{\text{beam}}}
{\text{m}_{\text{u}} \text{A}}
}
$$

### Define the integral coefficient

The implemented coefficient is

$$
\frac{3\text{v}_{\text{crit}}}
{\ln\left(1+\left(\frac{\text{v}_{\text{beam}}}{\text{v}_{\text{crit}}}\right)^3\right)}
$$

### Perform the fusion rate integral

The slowing-down-weighted integral is evaluated as

$$
\int_0^{\text{v}_{\text{beam}}/\text{v}_{\text{crit}}}
\frac{\text{u}^3}{1+\text{u}^3}
\sigma_{\text{bmfus}}(\text{E}_{\text{amu}}(\text{u}))
\,\mathrm{d}\text{u}
$$

where $u = v/v_{\text{crit}}$.

The quantity returned by `_beam_fusion_cross_section()` is in cm$^2$, so the integrated area is converted to m$^2$ before forming the final rate coefficient.

### Hot beam fusion reaction rate integrand | `_hot_beam_fusion_reaction_rate_integrand()`

This internal function evaluates the integrand used in the beam-target fusion rate coefficient.

The implemented integrand is

$$
\frac{\text{u}^3}{1+\text{u}^3}\sigma_{\text{bmfus}}(\text{E}_{\text{amu}}(\text{u}))
$$

where:

- $u$ is the ratio of the instantaneous beam speed to the critical speed,
- $E_{\text{amu}}(u)$ is the instantaneous beam kinetic energy per amu,
- $\sigma_{\text{bmfus}}$ is returned by [`_beam_fusion_cross_section()`](#beam-fusion-cross-section--_beam_fusion_cross_section).

The beam kinetic energy per amu used in the code is constructed from

$$
\text{E}_{\text{amu}}
=
\frac{(\text{u} \text{v}_{\text{crit}})^2 \text{m}_{\text{u}}}{\text{e}_{\text{keV}}}
$$

#### Beam fusion cross section | `_beam_fusion_cross_section()`

This internal function returns the beam fusion cross-section fit used by the beam-target rate coefficient calculation. Note: the provenance of this fit is currently unverified in the documentation.

The plasma ions are assumed to be stationary.

The implementation is

$$
\sigma_{\text{bm}}(\text{E}) =
\begin{cases}
1.0\times10^{-27}\ \text{cm}^2, & \text{E} < 10.0\ \text{keV} \\
8.0\times10^{-26}\ \text{cm}^2, & \text{E} > 10^4\ \text{keV} \\
\dfrac{
1.0\times10^{-24}
\left[
\dfrac{\text{a}_2}{1.0+(\text{a}_3 \text{E}-\text{a}_4)^2}+\text{a}_5
\right]
}{
\text{E}\left[\exp\left(\dfrac{\text{a}_1}{\sqrt{\text{E}}}\right)-1.0\right]
}
\ \text{cm}^2, & \text{otherwise}
\end{cases}
$$

where the beam energy used inside the fit is

$$
\text{E}
=
\frac{1}{2} \text{A}_{\text{D}} \text{E}_{\text{amu}}
$$

and the constants are

- $a_1 = 45.95$
- $a_2 = 5.02 \times 10^4$
- $a_3 = 1.368 \times 10^{-2}$
- $a_4 = 1.076$
- $a_5 = 4.09 \times 10^2$

The exact provenance of this particular fit and its coefficients has not yet been fully traced in the present documentation. It is retained here to document the current implementation.

### Multiply by the coefficient to get the full rate coefficient

The final effective beam-target fusion rate coefficient is

$$
\langle \sigma \text{v} \rangle_{\text{beam}}
=
\frac{3\text{v}_{\text{crit}}}
{\ln\left(1+\left(\frac{\text{v}_{\text{beam}}}{\text{v}_{\text{crit}}}\right)^3\right)}
\int_0^{\text{v}_{\text{beam}}/\text{v}_{\text{crit}}}
\frac{\text{u}^3}{1+\text{u}^3}
\sigma_{\text{bmfus}}(\text{E}_{\text{amu}}(\text{u}))
\,\mathrm{d}\text{u}
$$

---

## Beam-target fusion reaction rate | `beam_target_reaction_rate()`

The present implementation evaluates an effective beam-target fusion rate coefficient for a slowing-down fast-ion population. This is consistent with simple beam-plasma reaction-rate models in which fast ions interact with a Maxwellian background plasma during the slowing-down process[^niikura1990]. The total beam-target fusion reaction rate is calculated as

$$
\text{R}_{\text{beam-target}}
=
\text{n}_{\text{beam}}
\text{n}_{\text{target}}
\langle \sigma \text{v} \rangle_{\text{beam}}
\text{V}_{\text{plasma}}
$$

where:

- $n_{\text{beam}}$ is the hot beam ion density,
- $n_{\text{target}}$ is the thermal target ion density,
- $\langle \sigma v \rangle_{\text{beam}}$ is the effective beam-target rate coefficient,
- $V_{\text{plasma}}$ is the plasma volume.

This returns the total reaction rate in s$^{-1}$.

---

## Beam fusion alpha power | `alpha_power_beam()`

The beam-target alpha power is obtained from the total beam-target reaction rate by

$$
\text{P}_{\alpha,\text{beam}}
=
\text{R}_{\text{beam-target}} \text{E}_{\alpha}
$$

and is converted from W to MW in the implementation.

This returns the alpha power in MW.

---

## Full beam-target alpha power calculation in `beam_fusion()`

For the deuterium beam component reacting with thermal tritium:

$$
\text{R}_{\text{D-beam,DT}}
=
\langle \text{n}_{\text{beam}} \rangle_{\text{D}}
\text{n}_{\text{T,plasma}}
\langle \sigma \text{v} \rangle_{\text{beam,D}}
\text{V}_{\text{plasma}}
$$

$$
\text{P}_{\alpha,\text{D-beam}}
=
\text{R}_{\text{D-beam,DT}} \text{E}_{\alpha}
$$

For the tritium beam component reacting with thermal deuterium:

$$
\text{R}_{\text{T-beam,DT}}
=
\langle \text{n}_{\text{beam}} \rangle_{\text{T}}
\text{n}_{\text{D,plasma}}
\langle \sigma \text{v} \rangle_{\text{beam,T}}
\text{V}_{\text{plasma}}
$$

$$
\text{P}_{\alpha,\text{T-beam}}
=
\text{R}_{\text{T-beam,DT}} \text{E}_{\alpha}
$$

The returned beam alpha power is then

$$
\text{P}_{\alpha,\text{beam}}
=
\mathtt{beamfus0}
\left(
\text{P}_{\alpha,\text{D-beam}}
+
\text{P}_{\alpha,\text{T-beam}}
\right)
$$

---

## Key Constraints

### Hot beam ion density limit

This constraint can be activated by stating `icc = 7` in the input file.

The desired value of the hot ion beam density calculated from the code (`nd_beam_ions_out`) can be constrained using the input variable `f_nd_beam_electron`, which is the ratio of the beam density to the plasma electron density. It can be set as an iteration variable by setting `ixc = 7`.

## References

[^sheffield1994]: J. W. Sheffield, “The physics of magnetic fusion reactors,” *Rev. Mod. Phys.*, vol. 66, no. 3, pp. 1015–1103, Jul. 1994. <https://doi.org/10.1103/RevModPhys.66.1015>

[^deng1987]: B. Deng and G. A. Emmert, “Fast ion pressure in fusion plasma,” *Nuclear Fusion and Plasma Physics*, vol. 9, no. 3, pp. 136–141, 1987. Available: <https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf>

[^wesson2011]: J. Wesson, *Tokamaks*, 4th ed., Oxford Science Publications, 2011.

[^niikura1990]: S. Niikura and M. Nagami, “Improvement of fusion reactivity and fusion power multiplication factor in the presence of fast ions,” *Fusion Engineering and Design*, vol. 12, pp. 467–480, 1990.
