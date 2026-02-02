# Neutral Beam Injection Fusion

The main function called for calculating the fusion reactions produced by neutral beam injection is `beam_fusion()`
Due to the small contribution of fusion power from the neutral beams only D-T reactions are taken into account, as D-D additions to fusion power are deemed to be negligible.
The beam fusion calculations will only run if the calculated beam current is greater than 0. This is done by having a NBI heating and current drive configuration. 

The NBI parameters taken from the current drive module to be used in the beam fusion calculations are the beam current (`c_beam_total`), beam energy (`e_beam_kev`) and the tritium component of the beam (`f_beam_tritium`).

Please see the [H&CD section](../../eng-models/heating_and_current_drive/heating-and-current-drive.md) of the docs for more info.

------------------------

## Beam slowing down properties | `beam_fusion()`

1. **Calculate the beam ion slowing down time given by**:

    $$
    \tau_{\text{slow}} = 1.99 \times 10^{19}\left(A_{\text{D}}\left(1.0-f_{\text{tritium-beam}}\right)+(A_{\text{T}}f_{\text{tritium-beam}})\right)\frac{\langle T_{\text{e}}\rangle^{3/2}}{\langle n_{\text{e}} \rangle \Lambda_{\text{ie}}}
    $$

2. **Calculate the the beam critical energy**

    The alpha particles are born with an energy of 3.5 MeV and initially slow down mainly by collisions with electrons. At a critical energy $E_{\text{crit}}$ the rate of loss to the ions becomes equal to that to the electrons, and at lower energies the loss to the ions predominates.[^1]

    $$
    E_{\text{crit}} = 14.8 A T_{\text{e}}\left[\frac{1}{n_{\text{e} \ln{ \Lambda_{\text{e}}}}}\left[\Sigma \frac{n_j Z_j^2\ln{\Lambda_{\text{i}}}}{A_j}\right]\right]^{2/3} \ [\text{eV}]
    $$

    This can be approximated to:

    $$
    E_{\text{crit}} \approx 0.1AT_{10} \ \ [\text{MeV}]
    $$

    Though is currently implemented for deuterium as:

    $$
    E_{\text{crit,D}} = 14.8A_{\text{D}}T_{\text{e}}Z_{\text{eff}}^{2/3}\frac{\ln \Lambda_{\text{ie}}+4.0}{\Lambda_{\text{ie}}}
    $$

    The tritium critical energy is simply just scaled with the ratio of atomic mass numbers

    $$
    E_{\text{crit,T}} = E_{\text{crit,D}}\left(\frac{A_{\text{T}}}{A_{\text{D}}}\right)
    $$

    ---------------------------

    ### Derivation of beam slowing down rate and critical energy

    The rate of slowing down of a test particle of mass $M$, charge $Z\text{e}$ and energy $E$, due to Coulomb collisions with a background species off mass $m_{\text{j}}$, charge $Z_{\text{j}}\text{e}$, density $n_{\text{j}}$ and temperature $T_{\text{j}}$ is given by[^2]:

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
    \frac{d E}{d t}=-\frac{2 E}{\tau_{\text{slow}}}\left[1+\left(\frac{E}{E}\right)^{3 / 2}\right]
    $$

    The critical energy $E_{\text{c}}$ is given by

    $$
    E_c=\left[\frac{3 \sqrt{\pi}}{4} \frac{M^{3 / 2}}{n_{\text{e}} \sqrt{m}_{\text{e}}} \sum_j\left(\frac{n_{\text{j}} z_{\text{j}}^2}{m_{\text{j}}} \ln \Lambda_{\text{j}}\right) \frac{1}{\ln \Lambda_{\text{e}}}\right]^{2 / 3} k T_{\text{e}}
    $$

    Some authors take $\ln \Lambda_i=\ln \Lambda_{\text{e}}$ in the expression for $E_C$, but this is not a good approximation for fast ion slowing down in fusion plasmas since $\ln \Lambda_{\text{e}} \approx 17$, while $\ln \Lambda_i \approx 22$. The slowing down time, $\tau_{\text{slow}}$ is given by

    $$
    \tau_{\text{slow}}=\frac{3(k T_{\text{e}})^{3 / 2}}{4 \sqrt{2 \pi m_{\text{e}}} Z^2}\left(\frac{4 \pi \varepsilon_0}{\text{e}^2}\right)^2 \frac{M}{n_{\text{e}} \ln \Lambda_{\text{e}}}
    $$

    When the particle energy $E$ is above $E_C$ the contribution of the electrons to the slowing down is larger than that of the ions. The slowing down time $\tau_s$ is actually the time scale for $V$ to decrease due to electron drag, i.e. $\tau_{\text{slow}}=-V /(d V / d t) e^*$

    $\blacksquare$

    -------------------------

3. **Set the plasma tritium and ion densities**

    $$
    \mathtt{deuterium\_density = nd_plasma_fuel_ions_vol_avg * f\_deuterium\_plasma} \\
    \mathtt{tritium\_density = nd_plasma_fuel_ions_vol_avg * f\_tritium\_plasma}
    $$

4. **Calculate the beam alpha powers, density and deposited energy**

    [`beamcalc()`](#neutral-beam-alpha-power-and-ion-energy--beamcalc) is ran to find the alpha power from the beams, the beam densities and the total energy deposited into the plasma.

5. **Set the returned alpha power**

    $$
    P_{\alpha,\text{beam}} = \mathtt{beamfus0} \times \left(P_{\alpha,\text{D-beam}} + P_{\alpha,\text{T-beam}}\right)
    $$

6. **Calculate the neutral beam beta**

    $$
    \beta_{\text{beam}} = \mathtt{betbm0}\times \frac{2}{3}4.03\times10^{-22} \frac{n_{\text{beam}}E_{\text{hot,beam}}}{B_{\text{tot}}^2}
    $$

    The value of $E_{\text{hot,beam}}$ is the energy deposited by the fast beam ions into the plasma, NOT the initial energy of the beam.

------------------------

## Neutral beam alpha powers, density and deposited energy | `beamcalc()`

1. **Calculate the beam current densities**

    $$
    I_{\text{beam,D}} = I_{\text{beam}} \times \left(1-\mathtt{f\_tritium\_beam}\right) \\
    I_{\text{beam,T}} = I_{\text{beam}} \times \mathtt{f\_tritium\_beam}
    $$

2. **Calculate the characteristic time taken for the beam energy to comparable to that of the thermal energy**

    The attenuation of the beam energy as it penetrates into the plasma is given by[^3]:

    $$
    E_{\text{beam}} = E_{\text{beam,0}} \left[e^{-\frac{3t}{\tau_{\text{slow}}}}-\left(\frac{E_{\text{crit}}}{E_{\text{beam,0}}}\right)^{\frac{3}{2}}\left(1-e^{-\frac{3t}{\tau_{\text{slow}}}}\right)\right]^{\frac{2}{3}}
    $$

    Where the characteristic time taken for the beam energy to fall to that of the thermal energy:

    $$
    \tau_{\text{slow,D*}} = \frac{\tau_{\text{slow}}}{3}\ln\left(1+\left(\frac{E_{\text{beam,0}}}{E_{\text{crit,D}}}\right)^{\frac{3}{2}}\right) \\
    \tau_{\text{slow,T*}} = \frac{\tau_{\text{slow}}}{3}\ln\left(1+\left(\frac{E_{\text{beam,0}}}{E_{\text{crit,T}}}\right)^{\frac{3}{2}}\right)
    $$

    This is calculated for the deuterium and the tritium beam components

3. **Set the fast beam ion densities**

    $$
    \langle n_{\text{beam}} \rangle_{\text{D}} = \frac{I_{\text{beam,D}}\tau_{\text{slow,D*}}}{V_{\text{plasma}}} \\
    \langle n_{\text{beam}} \rangle_{\text{T}} = \frac{I_{\text{beam,T}}\tau_{\text{slow,T*}}}{V_{\text{plasma}}}
    $$

    We can also set the total hot ion beam density as:

    $$
    \langle n_{\text{beam}} \rangle_{\text{total}} = \langle n_{\text{beam}} \rangle_{\text{D}} + \langle n_{\text{beam}} \rangle_{\text{T}}
    $$

4. **Calculate the speeds of ions when at the critical energy**

    Assuming non-relativistic energies, we set the velocities of the deuterium and tritium particles when they have the critical energy:

    $$
    v_{\text{crit,D}} = \sqrt{\left(2m_{\text{D}}E_{\text{crit,D}}\right)} \quad
    v_{\text{crit,T}} = \sqrt{\left(2m_{\text{T}}E_{\text{crit,T}}\right)}
    $$

5. **Calculate the fast ion pressures**
    
    The fast ion pressure is set as:

    $$
    p_{\text{D}}=\frac{m_{\text{D}} \tau_{\text{slow}} v_{\text{crit,D}}^2 I_{\text{beam,D}}}{3 V_{\text{plasma}}} \times \\
    \underbrace{\left[\frac{x_c^2}{2}+\frac{1}{6} \ln \left(\frac{x_c^2+2 x_c+1}{x_c^2-x_c+1}\right)-\frac{1}{\sqrt{3}} \tan ^{-1}\left(\frac{2 x_c-1}{\sqrt{3}}\right)-\frac{1}{\sqrt{3}} \tan \left(\frac{1}{\sqrt{3}}\right)\right]}_{\mathtt{\_fast\_ion\_pressure\_integral(E_{\text{beam}},E_{\text{crit,D}})}} \\
    p_{\text{T}}=\frac{m_{\text{T}} \tau_{\text{slow}} v_{\text{crit,T}}^2 I_{\text{beam,T}}}{3 V_{\text{plasma}}} \times \\ \quad
    \underbrace{\left[\frac{x_c^2}{2}+\frac{1}{6} \ln \left(\frac{x_c^2+2 x_c+1}{x_c^2-x_c+1}\right)-\frac{1}{\sqrt{3}} \tan ^{-1}\left(\frac{2 x_c-1}{\sqrt{3}}\right)-\frac{1}{\sqrt{3}} \tan \left(\frac{1}{\sqrt{3}}\right)\right]}_{\mathtt{\_fast\_ion\_pressure\_integral(E_{\text{beam}},E_{\text{crit,T}})}}
    $$

    ---------------------------------

    ### Beam fast ion pressure integral | `fast_ion_pressure_integral()`

    This internal function is for returning the main integral value for the fast ion pressure. It takes the initial beam energy ($E_{\text{beam}}$) and the critical energy ($E_{\text{crit}}$) for the required ion species.

    This integral derivation is originally of the form derived by D.Baiquan et.al.[^2]

    #### Derivation

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

    $\blacksquare$

    --------------------------

6. **Calculate the deposited fast ion energy from the pressure**

    This can be simply done by applying the Ideal gas law approximation that, $P = \frac{1}{3}nmv_{\text{rms}^2} = \frac{2}{3}n \langle E \rangle$

    $$
    E_{\text{hot,D}} = \frac{3}{2}p_{\text{D}} \langle n_{\text{beam}} \rangle_{\text{D}} \\
    E_{\text{hot,T}} = \frac{3}{2}p_{\text{T}} \langle n_{\text{beam}} \rangle_{\text{T}}
    $$

    We can thus also define the total hot ion depoisted energy as:

    $$
    E_{\text{hot,total}} = \frac{\left( E_{\text{hot,D}}  \langle n_{\text{beam}} \rangle_{\text{D}}\right) + \left( E_{\text{hot,T}}  \langle n_{\text{beam}} \rangle_{\text{T}}\right)}{\langle n_{\text{beam}} \rangle_{\text{total}}}
    $$

7. **Calculate the hot ion species fusion rates**

    The D-T fusion reaction rate is calculated from the [`beam_reaction_rate()`](#beam-fusion-reaction-rate--beam_reaction_rate) function.

    ------------------------

    ### Beam fusion reaction rate | `beam_reaction_rate()`

    1. **Calculate beam velocity**

        The beam velocity ($v_{\text{beam}}$) is calculated from the inputted beam energy and relative ion mass of the beam by simply re-arranging the kinetic energy equation

    2. **Define the integral coefficient**

        $$
        \frac{3v_{\text{critical}}}{\ln\left(1+\frac{v_{\text{beam}}}{v_{\text{critical}}}\right)^{3}}
        $$

    3. **Perform the fusion rate integral**

        $$
        \int_0^{v_{\text{relative}}} \frac{u^3}{1+u^3}\sigma_{\text{bmfus}}(E_{\text{amu}})
        $$
        
        --------------------

        #### Hot Beam Fusion Reaction Rate Integrand | `_hot_beam_fusion_reaction_rate_integrand()`

        This function computes the integrand for the hot beam fusion reaction rate based on the ratio of beam velocity to the critical velocity and the critical velocity for electron/ion slowing down of the beam ion.

        The integrand function is:

        $$
        \int \frac{u^3}{1+u^3}\sigma_{\text{bmfus}}(E_{\text{amu}})
        $$

        Where $u$ is the inputted ratio of the beam to the critical velocity.
        $E_{\text{amu}}$ represents the beam kinetic energy per atomic mass unit.

        The calculated beam fusion cross section $\sigma_{\text{bmfus}}$ is calculated from [`_beam_fusion_cross_section()`](#beam-fusion-cross-section--_beam_fusion_cross_section)

        ------------------------

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

    4. **Multiply by the coefficient to get the full fusion rate**

        $$
        \frac{3v_{\text{critical}}}{\ln\left(1+\frac{v_{\text{beam}}}{v_{\text{critical}}}\right)^{3}}\int_0^{v_{\text{relative}}} \frac{u^3}{1+u^3}\sigma_{\text{bmfus}}(E_{\text{amu}})
        $$

    -------------------------    

8. **Calculate the alpha power produced by the hot ion species**

    The function [`alpha_power_beam()`](#beam-fusion-alpha-power--alpha_power_beam) is ran to calculate the alpha power produced by the deuterium and tritium fast ions.

    -----------------------

    ### Beam fusion alpha power | `alpha_power_beam()`

    1. **Calculate reactivity ratio**

        The ratio between the profile averaged reactivity for D-T reactions and the reactivity for the D-T reactions if the plasma is assumed to be homogeneously at the volume averaged ion temperature ($T_{\text{i}}$) is calculated.

        $$
        f_{\text{DT}} = \frac{\langle\langle \sigma v \rangle\rangle_{\text{DT}}}{\langle \sigma v \rangle_{\text{DT}}}
        $$

    2. **Calculate the alpha fusion power**

        The alpha powers from the deuterium and tritium beam components are calculated:

        $$
        P_{\alpha,\text{D}} = \langle n_{\text{beam}} \rangle_{\text{D}} n_{\text{i}} \langle \sigma v \rangle_{\text{beam}} E_{\alpha} V_{\text{plasma}} f_{\text{DT}} \\
        P_{\alpha,\text{T}} = \langle n_{\text{beam}} \rangle_{\text{T}} n_{\text{i}} \langle \sigma v \rangle_{\text{beam}} E_{\alpha} V_{\text{plasma}} f_{\text{DT}}
        $$

------------------------------

## Key Constraints

### Hot beam ion density limit

This constraint can be activated by stating `icc = 7` in the input file.

The desired value of the hot ion beam density calculated from the code (`nd_beam_ions_out`) can be constrained using the input variable, `f_nd_beam_electron`. Which is the ratio of the beam density to the plasma electron density. It can be set as an iteration variable by setting `ixc = 7`.

[^1]: J. W. Sheffield, “The physics of magnetic fusion reactors,” vol. 66, no. 3, pp. 1015–1103,Jul. 1994, doi: https://doi.org/10.1103/revmodphys.66.1015.
[^2]: Deng Baiquan and G. A. Emmert, “Fast ion pressure in fusion plasma,” Nuclear Fusion and Plasma Physics,vol. 9, no. 3, pp. 136–141, 2022, Available: https://fti.neep.wisc.edu/fti.neep.wisc.edu/pdf/fdm718.pdf  
[^3]: Wesson, J. (2011) Tokamaks. 4th Edition, 2011 Oxford Science Publications,International Series of Monographs on Physics, Volume 149.