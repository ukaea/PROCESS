## NBI Fusion

The main function called for calculating the fusion reactions produced by neutral beam injection is `beam_fusion()`

### Beam slowing down properties | `beam_fusion()`

Only take into account D-T beam fusion

1. Calculate the beam ion slowing down time given by:

    $$
    \tau_{\text{beam-slow}} = 1.99 \times 10^{19}\left(A_{\text{D}}\left(1.0-f_{\text{tritium-beam}}\right)+(A_{\text{T}}f_{\text{tritium-beam}})\right)\frac{\langle T_{\text{e}}\rangle^{3/2}}{\langle n_{\text{e}} \rangle \Lambda_{\text{ie}}}
    $$

2. The alpha particles are born with an energy of 3.5 MeV and initially slow down mainly by collisions with electrons. At a critical energy $E_{\text{crit}}$ the rate of loss to the ions becomes equal to that to the electrons, and at lower energies the loss to the ions predominates.[^1]

    $$
    E_{\text{crit}} = 14.8 A T_{\text{e}}\left[\frac{1}{n_{\text{e} \ln{ \Lambda_{\text{e}}}}}\left[\Sigma \frac{n_j Z_j^2\ln{\Lambda_{\text{i}}}}{A_j}\right]\right]^{2/3} \ [\text{eV}]
    $$

    This can be approximated to:

    $$
    E_{\text{crit}} \approx 0.1AT_{10} \ \ [\text{MeV}]
    $$

    ---------------------------

    #### Derivation of beam slowing down rate and critical energy

    The rate of slowing down of a test particle of mass $M$, charge $Z\text{e}$ and energy $E$, due to Coulomb collisions with a background species off mass $m_{\text{j}}$, charge $Z_{\text{j}}\text{e}$, density $n_{\text{j}}$ and temperature $T_{\text{j}}$ is given by:

    $$
    \frac{\mathrm{dE}}{\mathrm{dt}}=\left[-\Phi\left(x_{\text{j}}\right)+x_{\text{j}}\left(1+\frac{m_{\text{j}}}{M} \Phi^{\prime}\left(x_{\text{j}}\right)\right)\right] \frac{4 \pi n_{\text{j}}}{m_{\text{j}} V}\left(\frac{Z z_{\text{j}} Z \text{e}^2}{4 \pi \varepsilon_0}\right)^2 \ln \Lambda_{\text{j}},
    $$

    where $\Phi(x)$ is the error function,

    $$
    \Phi^{\prime}(x) = \frac{\mathrm{d\Phi}}{\mathrm{dx}}, \ \ V = \sqrt{\frac{2E}{M}}, \ \ V_{\text{j}} = \sqrt{\frac{2kT_{\text{j}}}{m_{\text{j}}}}, \ \ x_{\text{j}} = \frac{V}{V_{\text{j}}}
    $$

    and $\ln \Lambda_{\text{j}}$ is the usual Coulomb logarithm for the test particle and background species interactions.

    For the contribution to the slowing down due to interaction with the background ions, we can use a large argument expansion for the error function. This is because the fusion born ions have a velocity $V$, much greater than the thermal velocity $V_{\text{j}}$, of the background ions. The velocity of the fast ions is much less than the thermal velocity of electrons, however. For the electron contribution to the slowing down we use the small argument expansion for the error function. The net slowing down rate is then given by
    $$
    \frac{\mathrm{d E}}{\mathrm{d t}}=-\frac{A Z^2 \sqrt{M}}{\sqrt{E}}-\frac{B Z^2 E}{M}
    $$
    where the coefficients $A$ and $B$ are given by
    $$
    \begin{aligned}
    & A=\frac{4 \pi}{\sqrt{2}}\left(\frac{\text{e}^2}{4 \pi \varepsilon_0}\right)^2 \sum_j\left(\frac{n_{\text{j}} Z_{\text{j}}^2}{m_{\text{j}}} \ln \Lambda_j\right) \\
    & B=\frac{16 \sqrt{\pi}}{3 k T_{\text{e}}} \sqrt{\frac{m_{\text{e}}}{2 k T_{\text{e}}}}\left(\frac{\text{e}^2}{4 \pi \varepsilon_0}\right)^2 n_{\text{e}} \ln \Lambda_{\text{e}}
    \end{aligned}
    $$

    The sum over $\text{j}$ in $A$ is over the various ionic species. The quantities with subscript $\text{e}$ refer to the electrons.

    $\frac{\mathrm{d E}}{\mathrm{d t}}$ can be rewritten in the form,

    $$
    \frac{d E}{d t}=-\frac{2 E}{\tau_s}\left[1+\left(\frac{E}{E}\right)^{3 / 2}\right]
    $$

    The critical energy $E_{\text{c}}$ is given by

    $$
    E_c=\left[\frac{3 \sqrt{\pi}}{4} \frac{M^{3 / 2}}{n_{\text{e}} \sqrt{m}_{\text{e}}} \sum_j\left(\frac{n_{\text{j}} z_{\text{j}}^2}{m_{\text{j}}} \ln \Lambda_{\text{j}}\right) \frac{1}{\ln \Lambda_{\text{e}}}\right]^{2 / 3} k T_{\text{e}}
    $$

    Some authors take $\ln \Lambda_i=\ln \Lambda_{\text{e}}$ in the expression for $E_C$, but this is not a good approximation for fast ion slowing down in fusion plasmas since $\ln \Lambda_{\text{e}} \approx 17$, while $\ln \Lambda_i \approx 22$. The slowing down time, $\tau_s$ is given by

    $$
    \tau_s=\frac{3(k T_{\text{e}})^{3 / 2}}{4 \sqrt{2 \pi m_{\text{e}}} Z^2}\left(\frac{4 \pi \varepsilon_0}{\text{e}^2}\right)^2 \frac{M}{n_{\text{e}} \ln \Lambda_{\text{e}}}
    $$

    When the particle energy $E$ is above $E_C$ the contribution of the electrons to the slowing down is larger than that of the ions. The slowing down time $\tau_s$ is actually the time scale for $V$ to decrease due to electron drag, i.e. $\tau_{s}=-V /(d V / d t) e^*$

    $\blacksquare$

    -------------------------

3. Set the plasma tritium and ion densities.

    $$
    \mathtt{deuterium\_density = deni * f\_deuterium\_plasma} \\
    \mathtt{tritium\_density = deni * f\_tritium\_plasma}
    $$

4. Calculate the beam alpha powers

    [`beamcalc()`](#neutral-beam-alpha-power-and-ion-energy--beamcalc) is ran to find the alpha power from the beams and the beam densities.

5. Set the returned alpha power

    $$
    \mathtt{alpha_power_beams = beamfus0 * (palpdb + palptb)}
    $$

6. Calculate the neutral beam beta

    $$
    \mathtt{betanb = betbm0 * 4.03e-22 * 0.66666 * beam_density_out * ehotnb / (bt**2 + bp**2)}
    $$

------------------------

### Neutral beam alpha power and ion energy | `beamcalc()`

Set the fraction of the beam current for deuterium and tritium

calculate the velocity at which the deuterium ions are going at the critical speed

The characteristic time take for the beam energy to comparable to that of the thermal energy is:

$$
E_{\text{beam}} = E_{\text{beam,0}} \left[e^{-\frac{3t}{\tau_s}}-\left(\frac{E_{\text{crit}}}{E_{\text{beam,0}}}\right)^{\frac{3}{2}}\left(1-e^{-\frac{3t}{\tau_s}}\right)\right]^{\frac{2}{3}}
$$

$$
\tau = \frac{\tau_s}{3}\ln\left(1+\left(\frac{E_{\text{beam,0}}}{E_{\text{crit}}}\right)^{\frac{3}{2}}\right)
$$

[^3]

Calculate the alpha heating and energy from the fast ion pressure 
### Beam Current Fractions

1. **Calculate the fractions of beam current for deuterium and tritium:**

    ```python
    f_beam_current_deuterium = beam_current * (1.0 - f_tritium_beam)
    f_beam_current_tritium = beam_current * f_tritium_beam
    ```

2. **Determine the ratio of beam energy to critical energy for deuterium:**

    ```python
    beam_energy_ratio_deuterium = beam_energy / critical_energy_deuterium
    ```

3. **Calculate the characteristic time for the deuterium ions to slow down to the thermal energy:**

    ```python
    characteristic_deuterium_beam_slow_time = beam_slow_time / 3.0 * np.log(1.0 + (beam_energy_ratio_deuterium) ** 1.5)
    ```

4. **Compute the deuterium beam density:**

    ```python
    deuterium_beam_desnity = (
        (1.0 - f_tritium_beam) * beam_current * characteristic_deuterium_beam_slow_time / (constants.electron_charge * plasma_volume)
    )
    ```

5. **Determine the ratio of beam energy to critical energy for tritium:**

    ```python
    beam_energy_ratio_tritium = beam_energy / critical_energy_tritium
    ```

6. **Calculate the characteristic time for the tritium ions to slow down to the thermal energy:**

    ```python
    characteristic_tritium_beam_slow_time = beam_slow_time / 3.0 * np.log(1.0 + (beam_energy_ratio_tritium) ** 1.5)
    ```

7. **Compute the tritium beam density:**

    ```python
    tritium_beam_desnity = f_tritium_beam * beam_current * characteristic_tritium_beam_slow_time / (constants.electron_charge * plasma_volume)
    ```

8. **Calculate the total hot beam density:**

    ```python
    hot_beam_density = deuterium_beam_desnity + tritium_beam_desnity
    ```

### Critical Energy Speed

9. **Find the speed of the deuterium particle at critical energy:**

    ```python
    deuterium_critical_energy_speed = np.sqrt(
        2.0
        * constants.kiloelectron_volt
        * critical_energy_deuterium
        / (constants.atomic_mass_unit * ATOMIC_MASS_DEUTERIUM)
    )
    ```

10. **Find the speed of the tritium particle at critical energy:**

    ```python
    tritium_critical_energy_speed = np.sqrt(
        2.0
        * constants.kiloelectron_volt
        * critical_energy_tritium
        / (constants.atomic_mass_unit * ATOMIC_MASS_TRITIUM)
    )
    ```

### Source Term and Pressure Coefficients

11. **Calculate the source term for deuterium:**

    ```python
    source_deuterium = f_beam_current_deuterium / (constants.electron_charge * plasma_volume)
    ```

12. **Calculate the source term for tritium:**

    ```python
    source_tritium = f_beam_current_tritium / (constants.electron_charge * plasma_volume)
    ```

13. **Compute the pressure coefficient for deuterium:**

    ```python
    pressure_coeff_deuterium = (
        ATOMIC_MASS_DEUTERIUM
        * constants.atomic_mass_unit
        * beam_slow_time
        * deuterium_critical_energy_speed**2
        * source_deuterium
        / (constants.kiloelectron_volt * 3.0)
    )
    ```

14. **Compute the pressure coefficient for tritium:**

    ```python
    pressure_coeff_tritium = (
        ATOMIC_MASS_TRITIUM
        * constants.atomic_mass_unit
        * beam_slow_time
        * tritium_critical_energy_speed**2
        * source_tritium
        / (constants.kiloelectron_volt * 3.0)
    )
    ```

### Fast Ion Pressure and Deposited Energy

15. **Calculate the fast ion pressure for deuterium:**

    ```python
    deuterium_pressure = pressure_coeff_deuterium * _fast_ion_pressure_integral(beam_energy, critical_energy_deuterium)
    ```

16. **Calculate the fast ion pressure for tritium:**

    ```python
    tritium_pressure = pressure_coeff_tritium * _fast_ion_pressure_integral(beam_energy, critical_energy_tritium)
    ```

17. **Compute the beam deposited energy for deuterium:**

    ```python
    deuterium_depsoited_energy = 1.5 * deuterium_pressure / deuterium_beam_desnity
    ```

18. **Compute the beam deposited energy for tritium:**

    ```python
    tritium_depsoited_energy = 1.5 * tritium_pressure / tritium_beam_desnity
    ```

19. **Calculate the total deposited energy:**

    ```python
    total_depsoited_energy = ((deuterium_beam_desnity * deuterium_depsoited_energy) + (tritium_beam_desnity * tritium_depsoited_energy)) / hot_beam_density
    ```

### Reaction Rates and Alpha Power

20. **Calculate the hot deuterium reaction rate:**

    ```python
    hot_deuterium_rate = 1e-4 * beam_reaction_rate(ATOMIC_MASS_DEUTERIUM, deuterium_critical_energy_speed, beam_energy)
    ```

21. **Calculate the hot tritium reaction rate:**

    ```python
    hot_tritium_rate = 1e-4 * beam_reaction_rate(ATOMIC_MASS_TRITIUM, tritium_critical_energy_speed, beam_energy)
    ```

22. **Compute the deuterium beam alpha power:**

    ```python
    deuterium_beam_alpha_power = alpha_power_beam(deuterium_beam_desnity, nt, hot_deuterium_rate, plasma_volume, ti, svdt)
    ```

23. **Compute the tritium beam alpha power:**

    ```python
    tritium_beam_alpha_power = alpha_power_beam(tritium_beam_desnity, nd, hot_tritium_rate, plasma_volume, ti, svdt)
    ```

24. **Return the calculated values:**

    ```python
    return deuterium_beam_alpha_power, tritium_beam_alpha_power, hot_beam_density, total_depsoited_energy
    ```

------------------------------

### Beam fusion reaction rate | `beam_reaction_rate()`

1. Calculate beam velocity

    The beam velocity ($v_{\text{beam}}$) is calculated from the inputted beam energy and relative ion mass of the beam by simply re-arranging the kinetic energy equation

2. Define the integral coefficient

    $$
    \frac{3v_{\text{critical}}}{\ln\left(1+\frac{v_{\text{beam}}}{v_{\text{critical}}}\right)^{3}}
    $$

3. Perform the fusion rate integral

    $$
    \int_0^{v_{\text{relative}}} \frac{u^3}{1+u^3}\sigma_{\text{bmfus}}(E_{\text{amu}})
    $$

3. Multiply by the coefficient to get the full fusion rate

    $$
    \frac{3v_{\text{critical}}}{\ln\left(1+\frac{v_{\text{beam}}}{v_{\text{critical}}}\right)^{3}}\int_0^{v_{\text{relative}}} \frac{u^3}{1+u^3}\sigma_{\text{bmfus}}(E_{\text{amu}})
    $$


#### Hot Beam Fusion Reaction Rate Integrand | `_hot_beam_fusion_reaction_rate_integrand()`

This function computes the integrand for the hot beam fusion reaction rate based on the ratio of beam velocity to the critical velocity and the critical velocity for electron/ion slowing down of the beam ion.

The integrand function is:

$$
\int \frac{u^3}{1+u^3}\sigma_{\text{bmfus}}(E_{\text{amu}})
$$

Where $u$ is the inputted ratio of the beam to the critical velocity.
$E_{\text{amu}}$ represents the beam kinetic energy per atomic mass unit.

The calculated beam fusion cross section $\sigma_{\text{bmfus}}$ is calculated from [`_beam_fusion_cross_section()`](#beam-fusion-cross-section--_beam_fusion_cross_section)

-------------------------------

#### Beam fusion cross section | `_beam_fusion_cross_section()`

This internal function is used to find the beam cross section.
It sets limits on cross-section at low and high beam energies. The plasma ions are assumed to be stationary:

$$
\sigma_{\text{bm}}(E) = 
\begin{cases} 
1.0 \times 10^{-27} \ \text{cm}^2 & \text{if } E < 10.0 \ \text{keV/amu} \\
8.0 \times 10^{-26} \ \text{cm}^2 & \text{if } E > 10^4 \ \text{keV/amu} \\
\frac{1.0 \times 10^{-24} \cdot \left( \frac{a_2}{1.0 + (a_3 E - a_4)^2} + a_5 \right)}{E \left( \exp\left(\frac{a_1}{\sqrt{E}}\right) - 1.0 \right)} \ \text{cm}^2 & \text{otherwise}
\end{cases}
$$

where:

- \( E \) is the beam energy in, $\text{keV/amu}$

- \( a_1, a_2, a_3, a_4, a_5 \) are constants.

The constants are defined as:

- \( a_1 = 45.95 \)
- \( a_2 = 5.02 \times 10^4 \)
- \( a_3 = 1.368 \times 10^{-2} \)
- \( a_4 = 1.076 \)
- \( a_5 = 4.09 \times 10^2 \)

----------------------

### Beam fusion alpha power | `alpha_power_beam()`

### Beam fast ion pressure integral | `_fast_ion_pressure_integral()`

#### Derivation

#### Fast ion pressure derivation

[^2]

The fast ions, because of their finite slowing down time, develop a certain amount of pressure which has to be supported by the magnetic field.

This is in addition to the pressure of accumulated thermal "ash" in the plasma. To calculate this pressure we introduce a kinetic equation for the slowing down particles and solve it for their distribution function. Integration of the distribution function then determines the fast ion pressure.

We introduce the distribution function $g$ , of fast ions; $g$ is defined as the number of ions per unit energy per unit spatial volume and satisfies the steady-state kinetic equation,

$$
\frac{\partial}{\partial E}\left(g \frac{d E}{d t}\right)=S(E),
$$

where $S(E)$ is the source function in this "phase" space and $d E / d t$ is given by:

$$
\frac{dE}{dT} = -\frac{2E}{\tau_s}\left[1+\left(\frac{E_c}{E}^{\frac{3}{2}}\right)\right]
$$

This is the same form as above for deriving the critical energy and beam slow time in [`beam_fusion()`](#beam-slowing-down-properties--beam_fusion)

 This equation assumes the ions have a confinement time much longer than the slowing down time. For cases in which this is not true, an additional loss term would have to be introduced in the equation above. If we assume the ions are born monoenergetically, then:

$$
S(E)=S_0 \delta\left(E-E_0\right),
$$

where $\mathrm{S}_0$ is the number of ions born per unit time per unit volume and $\delta$ is the [Dirac delta function](https://en.wikipedia.org/wiki/Dirac_delta_function).

We impose the boundary condition that $\mathrm{g}=0$ for $\mathrm{E}>\mathrm{E}_0$. Our steady state kinetic distribution function can then be integrated to yield:

$$
g(E)= \begin{cases}\frac{S_0 T_s}{2 E\left(1+\left(E_c / E\right)^{3 / 2}\right)}, & E<E_0 \\ 0, & E>E_0\end{cases}
$$

If we consider the appropriate limiting cases.

The pressure $p$ of the fast ions is

$$
p=\frac{2}{3} \int_0^{E_0}  g(E) E  \ \ \mathrm{dE}
$$

which can be rewritten as

$$
p=\frac{M S_0 \tau_s V_c^2}{3} \int_0^{X_c}  \frac{x^4}{1+x^3} \quad \mathrm{dx}
$$

The integral above can also be evaluated analytically  to yield:

$$
p=\frac{M \tau_s V_c^2 S_0}{3}\left[\frac{x_c^2}{2}+\frac{1}{6} \ln \left(\frac{x_c^2+2 x_c+1}{x_c^2-x_c+1}\right) \\
-\frac{1}{\sqrt{3}} \tan ^{-1}\left(\frac{2 x_c-1}{\sqrt{3}}\right)-\frac{1}{\sqrt{3}} \tan \left(\frac{1}{\sqrt{3}}\right)\right]
$$

For most applications to fusion born particles, the dominant term is the first term in the square brackets.
This function returns the terms in the square brackets.

[^1]: J. W. Sheffield, “The physics of magnetic fusion reactors,” vol. 66, no. 3, pp. 1015–1103,Jul. 1994, doi: https://doi.org/10.1103/revmodphys.66.1015.
[^2]: Deng Baiquan and G. A. Emmert, “Fast ion pressure in fusion plasma,” Nuclear Fusion and Plasma Physics,vol. 9, no. 3, pp. 136–141, 2022, Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf  
[^3]: Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,International Series of Monographs on Physics, Volume 149.