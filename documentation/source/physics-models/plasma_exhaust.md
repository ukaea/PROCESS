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

## Divertor Heat Splits

### Brunner Model | `calculate_brunner_divertor_power_splits()`


[^1] M. Siccinio, G. Federici, R. Kembleton, H. Lux, F. Maviglia, and J. Morris,
"Figure of merit for divertor protection in the preliminary design of the
EU-DEMO reactor," Nuclear Fusion, vol. 59, no. 10, pp. 106026-106026,
Jul. 2019, doi: https://doi.org/10.1088/1741-4326/ab3153.

[^2] H. Zohm et al., "A stepladder approach to a tokamak fusion power plant,"
Nuclear Fusion, vol. 57, no. 8, pp. 086002-086002, May 2017, doi: https://doi.org/10.1088/1741-4326/aa739e.