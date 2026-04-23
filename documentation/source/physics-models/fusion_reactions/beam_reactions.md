# Beam Fusion Model

This page describes the neutral beam fusion model currently implemented in `beam_fusion()`.

The model is a reduced beam–target treatment for neutral beam ions injected into the plasma. PROCESS first determines an effective beam current that already accounts for beam transport effects such as shine-through losses. The beam-fusion model then uses this net beam current to build a steady population of fast ions and estimates how many of those ions fuse with the background plasma. `beam_fusion()` does not explicitly model beam attenuation, orbit effects, or spatial beam evolution, and instead treats the fast-ion population in a volume-averaged (0D) sense.

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

The beam slowing down time used in the model is implemented as [2]:

```math
\tau_{\text{slow}} =
1.99\times10^{19}
\left[
A_{\text{D}}\left(1-f_{\text{beam,T}}\right)
+
A_{\text{T}}f_{\text{beam,T}}
\right]
\frac{\langle T_{\text{e}}\rangle^{3/2}}
{\langle n_{\text{e}}\rangle \ln \Lambda_{\text{ie}}}
```

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

In the current implementation, the deuterium critical energy [^1] is

```math
E_{\text{crit,D}} =
14.8
A_{\text{D}}
T_{\text{e}}
Z_{\text{eff,mw}}^{2/3}
\frac{\ln\Lambda_{\text{ie}}+4.0}{\ln\Lambda_{\text{ie}}}
\qquad [\text{keV}]
```

where $Z_{\text{eff,mw}}$ is the mass-weighted effective plasma charge.

The tritium critical energy is then scaled by the beam ion mass ratio:

```math
E_{\text{crit,T}} =
E_{\text{crit,D}}
\left(\frac{A_{\text{T}}}{A_{\text{D}}}\right)
```

#### Derivation of beam slowing down rate and critical energy

The rate of slowing down of a test particle of mass $M$, charge $Ze$ and energy $E$, due to Coulomb collisions with a background species of mass $m_j$, charge $Z_j e$, density $n_j$ and temperature $T_j$, is given by[^2]

```math
\frac{dE}{dt}
=
\left[
-\Phi(x_j)
+
x_j\left(1+\frac{m_j}{M}\,\Phi^{\prime}(x_j)\right)
\right]
\frac{4\pi n_j}{m_j V}
\left(\frac{Z Z_j e^2}{4\pi\varepsilon_0}\right)^2
\ln \Lambda_j
```

where

```math
\Phi^{\prime}(x) = \frac{d\Phi}{dx}
```

```math
V = \sqrt{\frac{2E}{M}}, \qquad
V_j = \sqrt{\frac{2kT_j}{m_j}}, \qquad
x_j = \frac{V}{V_j}
```

For fast ions in fusion plasmas, the ion contribution and electron contribution may be approximated separately, leading to the standard slowing-down form[^3]:

```math
\frac{\mathrm{d}E}{\mathrm{d}t}
=
-\frac{A Z^2 \sqrt{M}}{\sqrt{E}}
-
\frac{B Z^2 E}{M}
```

with coefficients

```math
\begin{aligned}
A &=
\frac{4\pi}{\sqrt{2}}
\left(\frac{e^2}{4\pi\varepsilon_0}\right)^2
\sum_j
\left(
\frac{n_j Z_j^2}{m_j}\ln\Lambda_j
\right) \\
B &=
\frac{16\sqrt{\pi}}{3kT_e}
\sqrt{\frac{m_e}{2kT_e}}
\left(\frac{e^2}{4\pi\varepsilon_0}\right)^2
n_e \ln\Lambda_e
\end{aligned}
```

This can be rewritten in the form

```math
\frac{\mathrm{d}E}{\mathrm{d}t}
=
-\frac{2E}{\tau_{\text{slow}}}
\left[
1+\left(\frac{E_c}{E}\right)^{3/2}
\right]
```

where

```math
E_c
=
\left[
\frac{3\sqrt{\pi}}{4}
\frac{M^{3/2}}{n_e\sqrt{m_e}}
\sum_j
\left(
\frac{n_j Z_j^2}{m_j}\ln\Lambda_j
\right)
\frac{1}{\ln\Lambda_e}
\right]^{2/3}
kT_e
```

and

```math
\tau_{\text{slow}}
=
\frac{3(kT_e)^{3/2}}{4\sqrt{2\pi m_e} Z^2}
\left(\frac{4\pi\varepsilon_0}{e^2}\right)^2
\frac{M}{n_e \ln\Lambda_e}
```

In this regime, $\tau_{\text{slow}}$ is the characteristic electron-drag slowing-down timescale.

$\blacksquare$

#### **Set the plasma deuterium and tritium ion densities**

The bulk target ion densities used in the beam-target reactions are

```math
n_{\text{D,plasma}}
=
n_{\text{fuel}}
f_{\text{D,plasma}}
```

```math
n_{\text{T,plasma}}
=
n_{\text{fuel}}
f_{\text{T,plasma}}
```

where $n_{\text{fuel}}$ is the volume-averaged total fuel ion density.

---

### Calculate the beam alpha powers, beam densities and deposited energy

[`beam_slowing_down_state()`](#neutral-beam-alpha-power-beam-densities-and-deposited-energy--beam_slowing_down_state) is run to calculate the hot beam densities, critical speeds, and deposited beam energy.

[`beam_reaction_rate_coefficient()`](#beam-fusion-reaction-rate-coefficient--beam_reaction_rate_coefficient) is then run for the deuterium and tritium beam components to obtain effective beam-target rate coefficients.

[`beam_target_reaction_rate()`](#beam-target-fusion-reaction-rate--beam_target_reaction_rate) and [`alpha_power_beam()`](#beam-fusion-alpha-power--alpha_power_beam) are used to convert these into alpha power.

### Set the returned alpha power

The total neutral beam alpha power is

```math
P_{\alpha,\text{beam}}
=
\mathtt{beamfus0}
\left(
P_{\alpha,\text{D-beam}}
+
P_{\alpha,\text{T-beam}}
\right)
```

### Calculate the neutral beam beta

The neutral beam beta is computed from the hot beam density and the deposited beam energy:

```math
\beta_{\text{beam}}
=
\mathtt{betbm0}
\times
4.03\times10^{-22}
\times
\frac{2}{3}
\frac{n_{\text{beam,hot}} E_{\text{beam,deposited}}}
{B_{\phi}^2 + B_{\theta}^2}
```

where:

- $n_{\text{beam,hot}}$ is the total hot beam ion density,
- $E_{\text{beam,deposited}}$ is the density-weighted deposited beam ion energy in keV,
- $B_{\phi}$ is the toroidal magnetic field on axis,
- $B_{\theta}$ is the poloidal magnetic field.

The value of $E_{\text{beam,deposited}}$ is the pressure-equivalent deposited energy of the hot beam ions, **not** the initial beam injection energy.

## Neutral beam alpha power, beam densities and deposited energy | `beam_slowing_down_state()`

1. **Calculate the beam current fractions**

The beam current is split into deuterium and tritium components:

```math
I_{\text{beam,D}}
=
I_{\text{beam}}
\left(1-f_{\text{beam,T}}\right)
```

```math
I_{\text{beam,T}}
=
I_{\text{beam}}
f_{\text{beam,T}}
```

1. **Calculate the characteristic slowing-down time to the thermal range**

Using the classical slowing-down model, the characteristic time for the beam ion energy to slow from the birth energy to the thermal range is implemented as:

```math
\tau_{\text{slow,D}}^{*}
=
\frac{\tau_{\text{slow}}}{3}
\ln\left[
1+
\left(
\frac{E_{\text{beam}}}{E_{\text{crit,D}}}
\right)^{3/2}
\right]
```

```math
\tau_{\text{slow,T}}^{*}
=
\frac{\tau_{\text{slow}}}{3}
\ln\left[
1+
\left(
\frac{E_{\text{beam}}}{E_{\text{crit,T}}}
\right)^{3/2}
\right]
```

1. **Set the fast beam ion densities**

The steady-state hot beam ion densities are set from source rate times residence time:

```math
\langle n_{\text{beam}} \rangle_{\text{D}}
=
\frac{I_{\text{beam,D}}\tau_{\text{slow,D}}^{*}}
{e V_{\text{plasma}}}
```

```math
\langle n_{\text{beam}} \rangle_{\text{T}}
=
\frac{I_{\text{beam,T}}\tau_{\text{slow,T}}^{*}}
{e V_{\text{plasma}}}
```

The total hot beam ion density is then

```math
\langle n_{\text{beam}} \rangle_{\text{hot}}
=
\langle n_{\text{beam}} \rangle_{\text{D}}
+
\langle n_{\text{beam}} \rangle_{\text{T}}
```

1. **Calculate the speeds of ions at the critical energy**

Assuming non-relativistic energies, the beam ion speeds at the critical energy are

```math
v_{\text{crit,D}}
=
\sqrt{
\frac{2 e_{\text{keV}} E_{\text{crit,D}}}
{m_u A_{\text{D}}}
}
```

```math
v_{\text{crit,T}}
=
\sqrt{
\frac{2 e_{\text{keV}} E_{\text{crit,T}}}
{m_u A_{\text{T}}}
}
```

where:

- $e_{\text{keV}}$ is the conversion from keV to joules,
- $m_u$ is the atomic mass unit.

1. **Calculate the fast ion pressures**

First define the source rates per unit volume:

```math
S_{\text{D}}
=
\frac{I_{\text{beam,D}}}{e V_{\text{plasma}}}
```

```math
S_{\text{T}}
=
\frac{I_{\text{beam,T}}}{e V_{\text{plasma}}}
```

The implemented pressure coefficients are then

```math
C_{p,\text{D}}
=
\frac{
A_{\text{D}} m_u \tau_{\text{slow}} v_{\text{crit,D}}^2 S_{\text{D}}
}
{3 e_{\text{keV}}}
```

```math
C_{p,\text{T}}
=
\frac{
A_{\text{T}} m_u \tau_{\text{slow}} v_{\text{crit,T}}^2 S_{\text{T}}
}
{3 e_{\text{keV}}}
```

The fast ion pressures are

```math
p_{\text{D}}
=
C_{p,\text{D}}
\times
\underbrace{
\left[
\frac{x_c^2}{2}
+
\frac{1}{6}\ln\left(\frac{x_c^2+2x_c+1}{x_c^2-x_c+1}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{2x_c-1}{\sqrt{3}}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{1}{\sqrt{3}}\right)
\right]
}_{\mathtt{fast\_ion\_pressure\_integral}(E_{\text{beam}},E_{\text{crit,D}})}
```

```math
p_{\text{T}}
=
C_{p,\text{T}}
\times
\underbrace{
\left[
\frac{x_c^2}{2}
+
\frac{1}{6}\ln\left(\frac{x_c^2+2x_c+1}{x_c^2-x_c+1}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{2x_c-1}{\sqrt{3}}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{1}{\sqrt{3}}\right)
\right]
}_{\mathtt{fast\_ion\_pressure\_integral}(E_{\text{beam}},E_{\text{crit,T}})}
```

with

```math
x_c = \sqrt{\frac{E_{\text{beam}}}{E_{\text{crit}}}}
```

### Beam fast ion pressure integral | `fast_ion_pressure_integral()`

This internal function returns the dimensionless pressure integral factor used in the fast-ion pressure expression. It takes the beam birth energy $E_{\text{beam}}$ and the critical energy $E_{\text{crit}}$ for the required beam ion species.

#### Derivation

The fast ions, because of their finite slowing down time, develop a finite pressure which must be supported by the magnetic field.

To calculate this pressure, introduce the distribution function $g(E)$ of fast ions, where $g(E)$ is the number of ions per unit energy per unit volume, satisfying the steady-state kinetic equation

```math
\frac{\partial}{\partial E}
\left(
g \frac{\mathrm{d}E}{\mathrm{d}t}
\right)
=
S(E)
```

with slowing-down law

```math
\frac{\mathrm{d}E}{\mathrm{d}t}
=
-\frac{2E}{\tau_s}
\left[
1+\left(\frac{E_c}{E}\right)^{3/2}
\right]
```

If the ions are born monoenergetically,

```math
S(E)=S_0\delta(E-E_0)
```

and with the boundary condition $g(E)=0$ for $E>E_0$, the steady-state distribution becomes

```math
g(E)=
\begin{cases}
\dfrac{S_0\tau_s}{2E\left[1+\left(E_c/E\right)^{3/2}\right]}, & E<E_0 \\
0, & E>E_0
\end{cases}
```

The fast ion pressure is then

```math
p
=
\frac{2}{3}
\int_0^{E_0} g(E) E \, \mathrm{d}E
```

which can be written as

```math
p
=
\frac{M S_0 \tau_s V_c^2}{3}
\int_0^{x_c}
\frac{x^4}{1+x^3}\,\mathrm{d}x
```

and evaluated analytically as

```math
p
=
\frac{M S_0 \tau_s V_c^2}{3}
\left[
\frac{x_c^2}{2}
+
\frac{1}{6}\ln\left(\frac{x_c^2+2x_c+1}{x_c^2-x_c+1}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{2x_c-1}{\sqrt{3}}\right)
-
\frac{1}{\sqrt{3}}\tan^{-1}\left(\frac{1}{\sqrt{3}}\right)
\right]
```

`fast_ion_pressure_integral()` returns the bracketed dimensionless factor only.

$\blacksquare$

--
6. **Calculate the deposited fast ion energy from the pressure**

The deposited beam ion energy is inferred from the pressure relation

```math
p = \frac{1}{3}nmv_{\text{rms}}^2 = \frac{2}{3}n\langle E\rangle
```

so that the pressure-equivalent deposited energies are

```math
E_{\text{hot,D}}
=
\frac{3}{2}
\frac{p_{\text{D}}}{\langle n_{\text{beam}} \rangle_{\text{D}}}
```

```math
E_{\text{hot,T}}
=
\frac{3}{2}
\frac{p_{\text{T}}}{\langle n_{\text{beam}} \rangle_{\text{T}}}
```

with the convention that the deposited energy is set to zero if the corresponding hot beam density is zero.

The total deposited beam ion energy returned to `beam_fusion()` is the density-weighted average:

```math
E_{\text{hot,total}}
=
\frac{
E_{\text{hot,D}} \langle n_{\text{beam}} \rangle_{\text{D}}
+
E_{\text{hot,T}} \langle n_{\text{beam}} \rangle_{\text{T}}
}
{\langle n_{\text{beam}} \rangle_{\text{hot}}}
```

with the convention that this is set to zero if $\langle n_{\text{beam}} \rangle_{\text{hot}}=0$.

1. **Return the slowing-down state**

`beam_slowing_down_state()` returns:

- deuterium beam ion density,
- tritium beam ion density,
- deuterium critical speed,
- tritium critical speed,
- total hot beam ion density,
- density-weighted deposited beam energy.

## Beam fusion reaction rate coefficient | `beam_reaction_rate_coefficient()`

1. **Calculate the beam velocity**

The beam velocity is calculated from the input beam energy and beam ion mass:

```math
v_{\text{beam}}
=
\sqrt{
\frac{2 e_{\text{keV}} E_{\text{beam}}}
{m_u A}
}
```

1. **Define the integral coefficient**

The implemented coefficient is

```math
\frac{3v_{\text{crit}}}
{\ln\left(1+\left(\frac{v_{\text{beam}}}{v_{\text{crit}}}\right)^3\right)}
```

1. **Perform the fusion rate integral**

The slowing-down-weighted integral is evaluated as

```math
\int_0^{v_{\text{beam}}/v_{\text{crit}}}
\frac{u^3}{1+u^3}
\sigma_{\text{bmfus}}(E_{\text{amu}}(u))
\,\mathrm{d}u
```

where $u = v/v_{\text{crit}}$.

The quantity returned by `_beam_fusion_cross_section()` is in cm$^2$, so the integrated area is converted to m$^2$ before forming the final rate coefficient.

### Hot beam fusion reaction rate integrand | `_hot_beam_fusion_reaction_rate_integrand()`

This internal function evaluates the integrand used in the beam-target fusion rate coefficient.

The implemented integrand is

```math
\frac{u^3}{1+u^3}\sigma_{\text{bmfus}}(E_{\text{amu}}(u))
```

where:

- $u$ is the ratio of the instantaneous beam speed to the critical speed,
- $E_{\text{amu}}(u)$ is the instantaneous beam kinetic energy per amu,
- $\sigma_{\text{bmfus}}$ is returned by [`_beam_fusion_cross_section()`](#beam-fusion-cross-section--_beam_fusion_cross_section).

The beam kinetic energy per amu used in the code is constructed from

```math
E_{\text{amu}}
=
\frac{(u v_{\text{crit}})^2 m_u}{e_{\text{keV}}}
```

#### Beam fusion cross section | `_beam_fusion_cross_section()`

This internal function returns the beam fusion cross-section fit used by the beam-target rate coefficient calculation. Note: the provenance of this fit is currently unverified in the documentation.

The plasma ions are assumed to be stationary.

The implementation is

```math
\sigma_{\text{bm}}(E) =
\begin{cases}
1.0\times10^{-27}\ \text{cm}^2, & E < 10.0\ \text{keV} \\
8.0\times10^{-26}\ \text{cm}^2, & E > 10^4\ \text{keV} \\
\dfrac{
1.0\times10^{-24}
\left[
\dfrac{a_2}{1.0+(a_3 E-a_4)^2}+a_5
\right]
}{
E\left[\exp\left(\dfrac{a_1}{\sqrt{E}}\right)-1.0\right]
}
\ \text{cm}^2, & \text{otherwise}
\end{cases}
```

where the beam energy used inside the fit is

```math
E
=
\frac{1}{2} A_{\text{D}} E_{\text{amu}}
```

and the constants are

- $a_1 = 45.95$
- $a_2 = 5.02 \times 10^4$
- $a_3 = 1.368 \times 10^{-2}$
- $a_4 = 1.076$
- $a_5 = 4.09 \times 10^2$

The exact provenance of this particular fit and its coefficients has not yet been fully traced in the present documentation. It is retained here to document the current implementation.

1. **Multiply by the coefficient to get the full rate coefficient**

The final effective beam-target fusion rate coefficient is

```math
\langle \sigma v \rangle_{\text{beam}}
=
\frac{3v_{\text{crit}}}
{\ln\left(1+\left(\frac{v_{\text{beam}}}{v_{\text{crit}}}\right)^3\right)}
\int_0^{v_{\text{beam}}/v_{\text{crit}}}
\frac{u^3}{1+u^3}
\sigma_{\text{bmfus}}(E_{\text{amu}}(u))
\,\mathrm{d}u
```

## Beam-target fusion reaction rate | `beam_target_reaction_rate()`

The present implementation evaluates an effective beam-target fusion rate coefficient for a slowing-down fast-ion population. This is consistent with simple beam-plasma reaction-rate models in which fast ions interact with a Maxwellian background plasma during the slowing-down process[^4]. The total beam-target fusion reaction rate is calculated as

```math
R_{\text{beam-target}}
=
n_{\text{beam}}
n_{\text{target}}
\langle \sigma v \rangle_{\text{beam}}
V_{\text{plasma}}
```

where:

- $n_{\text{beam}}$ is the hot beam ion density,
- $n_{\text{target}}$ is the thermal target ion density,
- $\langle \sigma v \rangle_{\text{beam}}$ is the effective beam-target rate coefficient,
- $V_{\text{plasma}}$ is the plasma volume.

This returns the total reaction rate in s$^{-1}$.

## Beam fusion alpha power | `alpha_power_beam()`

The beam-target alpha power is obtained from the total beam-target reaction rate by

```math
P_{\alpha,\text{beam}}
=
R_{\text{beam-target}} E_{\alpha}
```

and is converted from W to MW in the implementation.

This returns the alpha power in MW.

## Full beam-target alpha power calculation in `beam_fusion()`

For the deuterium beam component reacting with thermal tritium:

```math
R_{\text{D-beam,DT}}
=
\langle n_{\text{beam}} \rangle_{\text{D}}
n_{\text{T,plasma}}
\langle \sigma v \rangle_{\text{beam,D}}
V_{\text{plasma}}
```

```math
P_{\alpha,\text{D-beam}}
=
R_{\text{D-beam,DT}} E_{\alpha}
```

For the tritium beam component reacting with thermal deuterium:

```math
R_{\text{T-beam,DT}}
=
\langle n_{\text{beam}} \rangle_{\text{T}}
n_{\text{D,plasma}}
\langle \sigma v \rangle_{\text{beam,T}}
V_{\text{plasma}}
```

```math
P_{\alpha,\text{T-beam}}
=
R_{\text{T-beam,DT}} E_{\alpha}
```

The returned beam alpha power is then

```math
P_{\alpha,\text{beam}}
=
\mathtt{beamfus0}
\left(
P_{\alpha,\text{D-beam}}
+
P_{\alpha,\text{T-beam}}
\right)
```

## Key Constraints

### Hot beam ion density limit

This constraint can be activated by stating `icc = 7` in the input file.

The desired value of the hot ion beam density calculated from the code (`nd_beam_ions_out`) can be constrained using the input variable `f_nd_beam_electron`, which is the ratio of the beam density to the plasma electron density. It can be set as an iteration variable by setting `ixc = 7`.

## References

1. J. W. Sheffield, “The physics of magnetic fusion reactors,” *Rev. Mod. Phys.*, vol. 66, no. 3, pp. 1015–1103, Jul. 1994. <https://doi.org/10.1103/RevModPhys.66.1015>

2. B. Deng and G. A. Emmert, “Fast ion pressure in fusion plasma,” *Nuclear Fusion and Plasma Physics*, vol. 9, no. 3, pp. 136–141, 1987. Available: <https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf>

3. J. Wesson, *Tokamaks*, 4th ed., Oxford Science Publications, 2011.

4. S. Niikura and M. Nagami, “Improvement of fusion reactivity and fusion power multiplication factor in the presence of fast ions,” *Fusion Engineering and Design*, vol. 12, pp. 467–480, 1990.
