# Fast Alpha Pressure Contribution | `fast_alpha_beta()`

The pressure contribution from the fast alpha particles can be controlled using switch `ifalphap`.
This sets the value of the physics variable, `beta_fast_alpha`.

A contribution from fast alphas to the total plasma pressure can be relatively high (~ 10-30%) for fusion temperatures of interest Because the maximum volume-averaged beta achievable in a tokamak is limited by  magnetohydrrdynamic (MHD) instabilities (e.g.,ballooning and kink modes), the presence of fast alphas, if $\langle \beta_{\text{tot}} \rangle$ is constant, reduces the background thermal plasma pressure. Furthermore, the energetic alpha population can influence (favorably or unfavorably) the bulk plasma ballooning mode stability boundaries[^2].

----------------------

## IPDG89 model

This model can be used by setting: `ifalphap` = 0[^1]

### Derivation

Below is the derivation given by Uckan, N. A. et.al. [^2].

$$
\beta_{\alpha} = \frac{2\mu_0 \left(\frac{2n_{\alpha}\langle E_{\alpha}\rangle }{3}\right)}{B^2}
$$

Normalizing to the plasma thermal beta $\beta_{\text{th}} = \beta_{\text{i}}+ \beta_{\text{e}}$

$$
\frac{\beta_{\alpha}}{\beta_{\text{th}}} = \frac{\left(\frac{2E_{\alpha,0}}{3T_{\text{e}}}\right) \left(\frac{n_{\alpha}}{n_{\text{e}}}\right) \left(\frac{\langle E_{\alpha} \rangle}{E_{\alpha,0}}\right)}{1+f_{\text{nT}}}
$$

where $f_{\text{nT}} = f_{\text{n}}f_{\text{T}} = \left(\frac{n_{\text{i}}}{n_{\text{e}}}\right)\left(\frac{T_{\text{i}}}{T_{\text{e}}}\right)$. For $T_{\text{i}} \approx T_{\text{e}}$ and $Z_{\text{eff}} \approx 1.5$ (with $Z = 6$, carbon), $f_{\text{nT}} \approx 0.9$ ad typical local values of fractional fast alpha density , beta and energy are given in the table below.

| $T \ [\mathrm{keV}]$ | $\langle\sigma v\rangle_{{D T}} \ \left[\mathrm{m}^3 / \mathrm{s}\right]$ | $n_{\mathrm{\alpha}} / n_{\text{e}} \ [\%]$ | $\beta_{\mathrm{\alpha}} / \beta_{\mathrm{th}} \ [\%]$ | $\langle E_{\alpha}\rangle/ E_{\alpha,0}$ |
| :---: | :---: | :---: | :---: | :---: |
| $5$ | $1.35 \times 10^{-23}$ | $0.01$ | $0.73$ | $0.3$ |
| $10$ | $1.13 \times 10^{-22}$ | $0.9$ | $4.2$ | $0.34$ |
| $20$ | $4.31 \times 10^{-22}$ | $0.8$ | $19$ | $0.39$ |
| $30$ | $6.65 \times 10^{-22}$ | $1.8$ | $31$ | $0.41$ |
| $40$ | $7.93 \times 10^{-22}$ | $2.7$ | $34$ | $0.41$ |
| $50$ | $8.54 \times 10^{-22}$ | $3.45$ | $34$ | $0.4$ |

Assuming a parabolic profile for temperature and density where $\alpha_{\text{n}} \approx 0.0 - 0.5$ (relatively flat density profile) and $\alpha_{\text{T}} \approx 1$ the above beta ratio becomes:

$$
\gamma_{\alpha} = \frac{\beta_{\alpha}}{\beta_{\text{th}}} \approx\frac{ 0.32 f_{\text{DT}}^2\left(\frac{T_{\text{i}}}{T_{\text{e}}}\right) \langle T_{\text{e,10}} \rangle^{5/2} \langle U_{\alpha \text{e}} \rangle}{1+ f_{\text{nT}}} \\
= 0.32 f_{\text{DT}}^2\left(\frac{T_{\text{i}}}{T_{\text{e}}}\right) \langle T_{\text{e,10}} \rangle^{5/2} \left[\frac{2^{5/2}}{(1+f_{\text{nT}})^{7/2}}\right]
$$

where $\langle U_{\alpha \text{e}} \rangle$ is the fraction of alpha energy given to the electrons, $\langle T \rangle = \langle n_{\text{e}}T_{\text{e}}+n_{\text{i}}T_{\text{i}} / \langle 2n_{\text{e}} \rangle = \langle T_{\text{e}} \rangle \left(1+f_{\text{nT}}\right)/2$ is the density-weighted average  temperature. For analytical simplicity, $\langle \sigma v \rangle_{\text{DT}}$ (the fusion reaction-rate parameter) is approximated as $\langle \sigma v \rangle_{\text{DT}} \approx 1.1 x 10^{-22} \left(T_{\text{i,10}}\right)^2$, which is accurate enough for $T \approx 7-20 \  \mathrm{keV}$. 

For the chosen profiles and $Z_{\text{eff}} \approx 1.5$, the average pressure contribution from fast alphas is $\gamma_{\alpha} \approx 5-20 \%$ for $\langle T \rangle \approx 6-15 \  \mathrm{keV}.$ Direct comparison between the predictions and a large number of 1-1/2-D $\mathtt{WHIST}$ transport code calculations (having similar profile shapes and $Z_{\text{eff}}$ values) shows good agreement (within ±15%) over the temperature range $\left(\langle T \rangle \approx 5-20 \  \mathrm{keV}\right)$ considered. A benchmark between above and $\mathtt{WHIST}$ has resulted in a simple functional fit that is more convenient to use in global analyses:

$$
\gamma_{\alpha} \approx 0.2\left(T_{\text{10}}-0.37\right), \ \text{for} \  Z_{\text{eff}} \approx 1.5, T_{\text{i}}/T_{\text{e}} \approx 1, \langle T \rangle \approx 5-20 \text{keV}
$$

!!! quote "Fit validity"

    "*To zeroth order, the assumption of different profiles $\left(\alpha_{\text{n}} \approx 0-1.0, \alpha_{\text{T}} \approx 0.5-2.0\right)$ did 
    not appear to have any significant effect on this simple fit. As expected, significant 
    deviations from the $\gamma_{\alpha}$ were seen in simulations for anomalous fast alpha diffusion 
    and energy relaxation. In such cases, however, global analysis is not adequate to describe 
    the fast alpha behavior*[^2]"

-------------------

The above derivation form is not that explicitly given in IPDG89[^1] and `PROCESS`. Terms for the electron and DT fuel ion species have been added back in. The $\left(\langle T_{\text{10}}\rangle-0.37\right)$ term still remains.

$$
\frac{\beta_{\alpha}}{\beta_{th}} = 0.29 \, \left( \langle T_{10} \rangle -
0.37 \right) \, \left( \frac{n_{\text{DT}}}{n_{\text{e}}} \right)^2
$$

For $Z_{\text{eff}} \approx 1.5, T_{\text{i}}/T_{\text{e}} \approx 1, \langle T \rangle \approx 5-20 \text{keV}$


-----------------------

## H.Lux model

This model can be used by setting: `ifalphap` = 1 (default)[^3]

$$
\frac{\beta_{\alpha}}{\beta_{th}} = 0.26 \, \left( \langle T_{10} \rangle -
0.65 \right)^{0.5} \, \left( \frac{n_{DT}}{n_e} \right)^2
$$


The latter model is a better estimate at higher temperatures.

[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989', https://inis.iaea.org/collection/NCLCollectionStore/_Public/21/068/21068960.pdf


[^2]: Uckan, N. A., Tolliver, J. S., Houlberg, W. A., and Attenberger, S. E. Influence of fast alpha diffusion and thermal alpha buildup on tokamak reactor performance. United States: N. p., 1987. Web.https://www.osti.gov/servlets/purl/5611706


[^3]: H. Lux, R. Kemp, D.J. Ward, M. Sertoli, 'Impurity radiation in DEMO 
systems modelling', Fus. Eng.  | Des. **101**, 42-51 (2015)