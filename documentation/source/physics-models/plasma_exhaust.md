# Plasma Exhaust | `PlasmaExhaust`

## Separatrix Power | `calculate_separatrix_power()`

The separatrix power is defined as the thermal and conducted power flowing outward from the core plasma across the separatrix into the scrape-off layer.

$$
\overbrace{P_{\text{sep}}}^{\texttt{p_plasma_separatrix_mw}} = \underbrace{f_{\alpha}P_{\alpha} + P_{\text{c}} + P_{\text{OH}} + P_{\text{HCD}}}_{\text{Plasma heating}} - \underbrace{P_{\text{rad}}}_{\text{Plasma Loss}}
$$

where $f_{\alpha}$ is the [fraction of alpha power that is coupled to the plasma](../physics-models/fusion_reactions/plasma_reactions.md#coupled-alpha-particle-power), $P_{\alpha}$ is the alpha power, $P_{\text{c}}$ is the charged particle power, $P_{\text{OH}}$ is the ohmic heating power, $P_{\text{HCD}}$ is the plasma heating done by the external heating & current drive systems and $P_{\text{rad}}$ is the total radiation power given off by the plasma.

--------------

## Divertor Protection Metrics

### `calculate_psep_over_r_metric()`

A simple metric for the divertor heat flux challenege is to simply divide by total separatrix power by the plasma major radius. This works as a stand in metric for the total available divertor area in which to spread the heat


$$
\texttt{p_plasma_separatrix_rmajor_mw}=\frac{P_{\text{sep}}}{R_0}
$$

------------

### EU-DEMO Protection Re-Attachment Metric | `calculate_eu_demo_re_attachment_metric()`

The divertor protection metric is obtained by starting from the target heat flux formula, inserting the Eich heat-flux-width scaling and $q_{95}-B_{\text{p}}$​ relationship, and reducing the resulting dependence to the engineering quantity whose value is constrained to remain below a reference EU-DEMO limit to ensure the divertor can withstand temporary reattachment without[^1][^2].

$$
\texttt{p_div_bt_q_aspect_rmajor_mw}=\frac{P_{\text{sep}}B_{\text{T}}}{q_{95}AR}
$$


-----------

## Divertor Separatrix Power Heat Splits

### Brunner Model | `calculate_brunner_divertor_power_splits()`

The Brunner model[^brunner_model] was made to estimate the power-sharing between four divertors in a dynamic double null configuration. Consisting of an in-out and up-down sharing scheme. The model considers the separatrix power flow going to each divertor rather than that which reaches the target plate after dissipation. The power split assumes that the fraction of power to the side walls is negligibly small and there is no radiation in the SOL. We thus have:

$$
P_{\text{sep}} \approx P_{\text{div,tot}} = P_{\text{lower,in}} + P_{\text{lower,out}} + P_{\text{upper,in}} + P_{\text{upper,out}}
$$

The main determining variable for the power split is $\delta R_{\text{sep}}$ representing the radial distance between the first and second separatrixes at the plasma outboard mid-plane.

A negative value of $\delta R_{\text{sep}}$ represents a dominant lower null configuration, A $\delta R_{\text{sep}}$ of 0 represents a perfect double null configuration and a positive $\delta R_{\text{sep}}$ is a upper null configuration.

It is observed that the majority of the power crosses the plasma boundary at the outer mid-plane. Assuming that	the	resultant heat flux profile	is exponential with	a single e-folding width equal to the scrape off layer power decay length, $\lambda_{\text{q,o}}$.

For the outboard, side the fraction of the outboard power going to the lower target is:

$$
\frac{P_{\text{lower,out}}}{P_{\text{lower,out}}+P_{\text{upper,out}}} = \frac{1}{1+e^{\frac{\delta R_{\text{sep}}}{\lambda_{\text{q,o}}}}}
$$

For the upper target:

$$
\frac{P_{\text{upper,out}}}{P_{\text{lower,out}}+P_{\text{upper,out}}} = \frac{1}{1+e^{-\frac{\delta R_{\text{sep}}}{\lambda_{\text{q,o}}}}}
$$

For the inboard divertors it was seen that the power sharing follows a similar logistic function trend but with a distinctly different power decay length, $\lambda_{\text{q,i}}$.

For the inboard side the fraction of the inboard power going to the lower target is:

$$
\frac{P_{\text{lower,in}}}{P_{\text{lower,in}}+P_{\text{upper,in}}} = \frac{1}{1+e^{\frac{\delta R_{\text{sep}}}{\lambda_{\text{q,i}}}}}
$$

The total inboard split compared to the outboard was found to give a Gaussian like dependence as follows:

$$
\frac{P_{\text{lower,in}}+P_{\text{lower,out}}}{P_{\text{lower,in}}+P_{\text{upper,in}}+P_{\text{lower,out}}+P_{\text{upper,out}}} = \\ P_{\text{in,0}}+(P_{\text{in,0}}-P_{\text{in,}\infin})\times \left(1-\frac{2}{1+e^{-\left(\frac{\delta R_{\text{sep}}}{\lambda_{\text{q,io}}}\right)^2}}\right)
$$

where $P_{\text{in,0}}$ is the fraction of power to the inner divertors at $\delta R_{\text{sep}} = 0$, $P_{\text{in,}\infin}$ is the fraction of the power to the inner divertors at $R_{\text{sep}} = \infin$

Fitted values to data from DIII-D[^d3_d_values] have shown the values of $P_{\text{in,0}}$ and $P_{\text{in,}\infin}$ to be 0.16 and 0.41 respectively which is what we use here. We also assume that out Gaussian width in the total power sharing is equal to the outboard power decay length, $\lambda_{\text{q,io}} = \lambda_{\text{q,o}}$.

In `PROCESS` we normally want to control how much of the separatrix power is directed to either the upper or lower divertors, this is done through changing `f_p_div_lower_separatrix`. Fitting to the DIII-D[^d3_d_values] we calculate the value of $\delta R_{\text{sep}}$ via:

$$
\delta R_{\text{sep}} = -3\times 10^{-3} \tanh^{-1}\left(2(\texttt{f_p_div_lower_separatrix}-0.5)\right)
$$

Truncated values of the $\delta R_{\text{sep}}$ are set to $1.5\times 10^{-2}\ \text{m}$ and $-1.5\times 10^{-2} \ \text{m}$ respectively based on when the power split gets close to a perfect upper or lower null.

--------------------


[^1] M. Siccinio, G. Federici, R. Kembleton, H. Lux, F. Maviglia, and J. Morris,
"Figure of merit for divertor protection in the preliminary design of the
EU-DEMO reactor," Nuclear Fusion, vol. 59, no. 10, pp. 106026-106026,
Jul. 2019, doi: https://doi.org/10.1088/1741-4326/ab3153.

[^2] H. Zohm et al., "A stepladder approach to a tokamak fusion power plant,"
Nuclear Fusion, vol. 57, no. 8, pp. 086002-086002, May 2017, doi: https://doi.org/10.1088/1741-4326/aa739e.

[^brunner_model] D. Brunner, A. Q. Kuang, B. LaBombard, and J. L. Terry, “The dependence of
divertor power sharing on magnetic flux balance in near double-null configurations on 
Alcator C-Mod,” Nuclear Fusion, vol. 58, no. 7, p. 076010, May 2018, doi: https://doi.org/10.1088/1741-4326/aac006.

[^d3_d_values] T. W. Petrie et al., “The effect of divertor magnetic balance on H-mode performance in DIII-D,” Journal of Nuclear Materials, vol. 290-293, pp. 935-939, Mar. 2001, doi: https://doi.org/10.1016/S0022-3115(00)00492-X