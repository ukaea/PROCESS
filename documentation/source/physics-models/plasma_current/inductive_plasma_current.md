# Inductive Current

Currently in `PROCESS` the inductive current fraction from the CS is not calculated directly but is just equal to ($1 - \texttt{f_c_plasma_non_inductive}$). Where $\texttt{f_c_plasma_non_inductive}$ is the sum of the fractions of current driven by non inductive means.

This calculated fraction (`f_c_plasma_inductive`) is then used in the `calculate_volt_second_requirements()` and `burn()` functions to calculate the volt-second requirements and the burn time for a pulsed machine.

!!! info "Inductive plasma current fraction refactor"

    It is hoped for the near future to have a more engineering based calculation of the fraction of the plasma current driven by the central solenoid. This would hopefully allow the setting of required ramp and flat-top times for which a given inductive current fraction can be given based on the operational performance and margin in the central solenoid.


--------------------

## Volt-second requirements | `calculate_volt_second_requirements()`

The plasma requires a constants magnetic flux change in order to keep inductively driving a current in itself.

By Faraday's law of induction, any change in flux through a circuit induces an electromotive force (EMF, $V$),
in the circuit, proportional to the rate of change of flux.

$$
V(t) = - \frac{\mathrm{d}}{\mathrm{d}t} \Phi(t)
$$

Inductance, $L$ is the ratio between the induced voltage from the flux change and the rate of change of current that produced the flux change.

$$
V(t) = L \ \frac{\mathrm{d}I}{\mathrm{d}t} 
$$


The flux requirements are defined by the sum of the pulse ramp up and flat top / burn phases.

-----------

### Current ramp phase

#### Resistive Component

In the ramp up phase we need to take the plasma current from 0 Amps to the plasma current value, which is normally ten's of Mega Amps. Since the plasma has a non zero resistance, the current induced in the plasma will be dissipated due to resistive losses. This can be tricky to calculate as the plasma resistance decreases as the plasma temperature increases towards our flat-top full plasma scenario state. The inductance of the plasma also varies during this ramp phase and so the amount of current driven for the same change in flux will vary too. 

Thankfully a formulation of the resistive flux consumption during ramp phase is given by Ejima et.al [^1].

$$
\overbrace{\Phi_{\text{res,ramp}}}^{\texttt{vs_res_ramp}} = C_{\text{eji}}\mu_0I_{\text{p}}R
$$

where $C_{\text{ejima}}$ is the empirical Ejima coefficient defined by the user or at its default of 0.4.

This relation is is done by analyzing a wide range of cross-sections,
ranging from circular to doublet produced in the Doublet III machine. The plasma cross-section is
slightly elongated, with $\kappa=1.2$. The toroidal field is $2.4 \ \text{T}$. The current swing of the Ohmic-heating transformer is nominally from $-25 \  \text{kA}$ to $+80 \ \text{kA}$, corresponding to a flux swing of $\approx 2.6 \ \text{Vs}$.

The initial one-turn loop voltage is around $50 \ \text{V}$, causing the plasma current to rise to $300 \ \text{kA}$ in about $80 \  \text{ms}$. The plasma current is then increased to higher flat-top currents at a steady rate of $2 \text{MA} \text{s}^{-1}$.

The calculation of $\Phi_{\text{res,ramp}}$ is based on the fact that the ramp up takes of the order of the resistive current penetration time $\frac{\mu_0 a^2}{\rho_{\text{p}}}$, where $\rho_{\text{p}}$ is the resistivity of the plasma. The
flux consumption is therefore independent of the resistivity and the minor radius.

-------------

#### Self-Inductance Component


The internal inductance is defined as the part of the inductance obtained by integrating over the plasma volume.

The internal component of the plasma self inductance flux consumption during the ramp up phase is given by:

$$
\overbrace{\Phi_{\text{ind,internal,ramp}}}^{\texttt{vs_plasma_internal}} =  I_{\text{p}} \times \underbrace{\left[\frac{\mu_0 R l_i}{2}\right]}_{\texttt{ind_plasma_internal}}
$$

Here $l_i$ is known as the normalised internal inductance defined for circular cross section plasmas with minor radius $a$:

$$
l_i = \frac{\langle B_{\text{p}}^2 \rangle }{B_{\text{p}}^2(a)}
$$

You may also see the following approximation term used, $l_i(3)$[^2].

$$
l_i(3) = \frac{2V\langle B_{\text{p}}^2 \rangle }{\mu_0^2I^2R}
$$

which is equal to $l_i$ if the plasma has a perfect circular cross-section.

The external inductance needs to be accounted for also as even though we assume the toroidal current density vanishes at the plasma edge, there still exists a vacuum poloidal magnetic field around the plasma.

Hirshman et.al[^3] gives a formula for the external plasma inductance in the form:

$$
\overbrace{L_{\text{ext}}}^{\texttt{ind_plasma_external}} = \mu_0 R\frac{a(\epsilon)(1-\epsilon)}{1-\epsilon + b(\epsilon)\kappa}
$$

$$
a(\epsilon)  = (1+1.18\sqrt{\epsilon}+2.05\epsilon)\ln{\left(\frac{8}{\epsilon}\right)} \\
- (2.0+9.25\sqrt{\epsilon}-1.21\epsilon)
$$

$$
b(\epsilon) = 0.73\sqrt{\epsilon}(1+2\epsilon^4 - 6\epsilon^5+3.7\epsilon^6)
$$

where $\epsilon$ is the plasma inverse aspect ratio and $\kappa$ is the separatrix elongation.


The total plasma inductance is then calculated as:

$$
\texttt{ind_plasma_total} = \texttt{ind_plasma_external} + \texttt{ind_plasma_internal}
$$

Therefore the total inductive flux consumption during ramp up is given by:

$$
\overbrace{\Phi_{\text{ind,ramp}}}^{\texttt{vs_self_ind_ramp}} =  \texttt{ind_plasma_total} \times I_{\text{p}}
$$

So the total resisitive and inductive flux consumption at current ramp up which is the total flux requires is given by:

$$
\overbrace{\Phi_{\text{tot,ramp}}}^{\texttt{vs_plasma_ramp_required}} =  \overbrace{\Phi_{\text{res,ramp}}}^{\texttt{vs_res_ramp}} + \overbrace{\Phi_{\text{ind,ramp}}}^{\texttt{vs_self_ind_ramp}}
$$

------------------

### Steady current burn phase

At plasma current flat-top there is no self inductance contribution as the plasma current is expected to be at a constant value. However there is still a resistive contribution.



#### Resistive Component


For the flat top resistive component we can just take the loop voltage value based on the plasmas resistivity and the known fraction of current driven inductively:

$$
\overbrace{V_{\text{loop}}}^{\texttt{v_burn_resistive}} = I_{\text{p}}  \rho_\text{p} f_{\text{ind}}
$$

where $\rho_\text{p}$ is the calculated [plasma resistivity](./plasma_resistive_heating.md) and $f_{\text{ind}}$ is the inductive current fraction.

The total flux required is then simply found by multiplying the loop voltage above by the required duration of the burn phase:

$$
\overbrace{\Phi_{\text{res,burn}}}^{\texttt{vs_burn_required}} = \overbrace{I_{\text{p}}  \rho_\text{p} f_{\text{ind}}}^{\texttt{v_burn_resistive}} \times \overbrace{T_{\text{burn}}}^{\texttt{t_plant_pulse_burn}}
$$

----------------

Finally we can now find the minimum flux required for the full duration of the pulse:

$$
\overbrace{\Phi_{\text{tot}}}^{\texttt{vs_total_required}} = \overbrace{\Phi_{\text{tot,ramp}}}^{\texttt{vs_plasma_ramp_required}} + \overbrace{\Phi_{\text{res,burn}}}^{\texttt{vs_burn_required}}
$$



[^1]: S. Ejima, R. W. Callis, J. L. Luxon, R. D. Stambaugh, T. S. Taylor, and J. C. Wesley, “Volt-second analysis and consumption in Doublet III plasmas,” Nuclear Fusion, vol. 22, no. 10, pp. 1313-1319, Oct. 1982, doi: https://doi.org/10.1088/0029-5515/22/10/006.
[^2]: G. L. Jackson et al., “ITER startup studies in the DIII-D tokamak,” Nuclear Fusion, vol. 48, no. 12, p. 125002, Nov. 2008, doi: https://doi.org/10.1088/0029-5515/48/12/125002.
[^3]: S. P. Hirshman and G. H. Neilson, “External inductance of an axisymmetric plasma,” The Physics of Fluids, vol. 29, no. 3, pp. 790-793, Mar. 1986, doi: https://doi.org/10.1063/1.865934.
‌