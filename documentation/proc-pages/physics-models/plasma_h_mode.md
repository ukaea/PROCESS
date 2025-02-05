## L-H Power Threshold Scalings

Transitions from a standard confinement mode (L-mode) to an improved
confinement regime (H-mode), called L-H transitions, are observed in most
tokamaks. A range of scaling laws are available that provide estimates of the
heating power required to initiate these transitions, via extrapolations
from present-day devices. PROCESS calculates these power threshold values
for the scaling laws listed in the table below, in routine `pthresh`.

For an H-mode plasma, use input parameter `ilhthresh` to
select the scaling to use, and turn on constraint equation no. 15 with
iteration variable no. 103 (`flhthresh`). By default, this will ensure
that the power reaching the divertor is at least equal to the threshold power
calculated for the chosen scaling, which is a necessary condition for
H-mode. 

For an L-mode plasma, use input parameter `ilhthresh` to
select the scaling to use, and turn on constraint equation no. 15 with 
iteration variable no. 103 (`flhthresh`). Set lower and upper bounds for 
the f-value `boundl(103) = 0.001` and `boundu(103) = 1.0` 
to ensure that the power does not exceed the calculated threshold, 
and therefore the machine remains in L-mode.

--------------------

### ITER-1996 Scalings

The general form is:

$$
P_{\text{L-H}} = 0.45 (0.6n_{e,20}R^2)^{\alpha} n^{0.75}_{\text{e},20}B_{\text{T}}R^2
$$

where $\alpha$ lies in the range of $-0.25 \le \alpha \le 0.25$.

------------------

#### ITER-1996 Nominal Scaling

Is selected with `ilhthresh = 1` [^1] [^2]

$$
P_{\text{L-H}} = 0.45 \times n^{0.75}_{\text{e},20}B_{\text{T}}R^2
$$

---------------

#### ITER-1996 Upper Scaling

Is selected with `ilhthresh = 2` [^1] [^2]

$$
P_{\text{L-H}} = 0.3960502816 \times n_{\text{e},20}B_{\text{T}}R^{2.5}
$$

---------------

#### ITER-1996 Lower Scaling

Is selected with `ilhthresh = 3` [^1] [^2]

$$
P_{\text{L-H}} = 0.5112987149 \times n_{\text{e},20}^{0.5}B_{\text{T}}R^{1.5}
$$

---------------

### Martin 2008 Scalings

The general form is:

$$
P_{\text{L-H}} = 0.0488 e^{\pm0.057} \bar{n}_{\text{e},20}^{0.717 \pm 0.035} B_{\text{T}}^{0.803 \pm 0.032} S_{\text{p}}^{0.941 \pm 0.019}
$$

where $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20} \text{m}^{-3}$, $B_{\text{T}}$ is the toroidal magnetic field in Tesla and $S_{\text{p}}$ is the plasma surface area in $\text{m}^2$.

------------------

#### Martin 2008 Nominal Scaling

Is selected with `ilhthresh = 6` 

$$
P_{\text{L-H}} = 0.0488 \bar{n}_{\text{e},20}^{0.717} B_{\text{T}}^{0.803} S_{\text{p}}^{0.941}
$$

---------------

#### Martin 2008 Upper Scaling

Is selected with `ilhthresh = 7` 

$$
P_{\text{L-H}} = 0.05166240355 \times \bar{n}_{\text{e},20}^{0.752} B_{\text{T}}^{0.835} S_{\text{p}}^{0.96}
$$

---------------

#### Martin 2008 Lower Scaling

Is selected with `ilhthresh = 8` 

$$
P_{\text{L-H}} = 0.04609619059 \times \bar{n}_{\text{e},20}^{0.682} B_{\text{T}}^{0.771} S_{\text{p}}^{0.922}
$$

---------------

### Snipes 2000 Scalings

!!! quote "Excerpt fom Snipes et.al"
    *"A number of regression fits were performed on the full data set of all 10 tokamaks for all
    L–H transition points that fit the standard SELDB2 criteria in deuterium plasmas ($N$ = 702).
    A small improvement in the RMSE was obtained by correcting the total ICRF power from
    Alcator C-Mod with density- and toroidal-field-dependent corrections for the absorbed ICRF
    power based on experimental measurements. For hydrogen minority heating between 5 and
    6 $T$, $P_{\text{abs}} = 0.9\bar{n}_{\text{e}}^{-0.6} P_{\text{ICRH}}$ while for He3 minority heating above 6 T, the assumption is made
    that $P_{\text{abs}} = 0.75 P_{\text{ICRH}}$."*

- This scaling has a RMSE of 26.8%

- $P_{\text{L-H}}$ is defined as $\left(P_{\text{in}} - \frac{dW}{dt}\right)$


The general form is:

$$
P_{\text{L-H}} = 1.42\pm 0.127 \times \bar{n}_{\text{e},20}^{0.58 \pm 0.035} B_{\text{T}}^{0.82 \pm 0.031} R^{1.00 \pm 0.089} a^{0.81 \pm 0.066}
$$

where $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20} \text{m}^{-3}$, $B_{\text{T}}$ is the toroidal magnetic field in Tesla, $R$ is the plasma major radius in metres and $a$ is the plasma major radius in metres.

------------------

#### Snipes 2000 Nominal Scaling

Is selected with `ilhthresh = 9` 

$$
P_{\text{L-H}} = 1.42 \times \bar{n}_{\text{e},20}^{0.58} B_{\text{T}}^{0.82} R^{1.00} a^{0.81}
$$

---------------

#### Snipes 2000 Upper Scaling

Is selected with `ilhthresh = 10` 

$$
P_{\text{L-H}} = 1.547 \times \bar{n}_{\text{e},20}^{0.615} B_{\text{T}}^{0.851} R^{1.089} a^{0.876}
$$

---------------

#### Snipes 2000 Lower Scaling

Is selected with `ilhthresh = 11` 

$$
P_{\text{L-H}} = 1.293 \times \bar{n}_{\text{e},20}^{0.545} B_{\text{T}}^{0.789} R^{0.911} a^{0.744}
$$

---------------

### Snipes 2000 Closed Divertor Scalings

!!! quote "Excerpt fom Snipes et.al"
    *"Several machines have reported changes in the H-mode threshold due to divertor geometry modifications. ASDEX-Upgrade saw a 15% increase with a more-closed divertor, attributed to higher edge densities. Conversely, JET and JT-60U experienced a 20% decrease after installing more-closed divertors. Alcator C-Mod, with its inherently closed divertor, showed no change in threshold despite further closure of bypass gaps. Preliminary experiments with C-Mod's divertor bypass indicate no threshold change, even with a significant drop in divertor neutral pressure. These variations contribute to the scatter in H-mode threshold regression fits."*

    *"Only four tokamaks have closed-divertor data in the database (Alcator C-Mod, ASDEX-Upgrade, JET, and JT-60U). Although the resulting data set is limited ($N$ = 169)"*

- This scaling has a RMSE of 22.0%

- $P_{\text{L-H}}$ is defined as $\left(P_{\text{in}} - \frac{dW}{dt}\right)$


The general form is:

$$
P_{\text{L-H}} = 0.8\pm 0.067 \times \bar{n}_{\text{e},20}^{0.50 \pm 0.061} B_{\text{T}}^{0.53 \pm 0.058} R^{1.51 \pm 0.077}
$$

where $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20} \text{m}^{-3}$, $B_{\text{T}}$ is the toroidal magnetic field in Tesla and $R$ is the plasma major radius in metres

------------------

#### Snipes 2000 Closed Divertor Nominal Scaling

Is selected with `ilhthresh = 12` 

$$
P_{\text{L-H}} = 0.8 \times \bar{n}_{\text{e},20}^{0.50} B_{\text{T}}^{0.53} R^{1.51}
$$

---------------

#### Snipes 2000 Closed Divertor Upper Scaling

Is selected with `ilhthresh = 13` 

$$
P_{\text{L-H}} = 0.867 \times \bar{n}_{\text{e},20}^{0.561} B_{\text{T}}^{0.588} R^{1.587}
$$

---------------

#### Snipes 2000 Closed Divertor Lower Scaling

Is selected with `ilhthresh = 14` 

$$
P_{\text{L-H}} = 0.733 \times \bar{n}_{\text{e},20}^{0.439} B_{\text{T}}^{0.472} R^{1.433}
$$

---------------

###  Hubbard 2012 L-I Scalings

The general form is:

$$
P_{\text{L-H}} = 2.11 \times I_{\text{p}}^{0.94\pm 0.24}\bar{n}_{\text{e},20}^{0.65\pm 0.18}
$$

where $I_{\text{p}}$ is the plasma current in $\text{MA}$ and $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20}$.

------------------

#### Hubbard 2012 L-I Nominal Scaling

Is selected with `ilhthresh = 15`

$$
P_{\text{L-H}} = 2.11 \times I_{\text{p}}^{0.94}\bar{n}_{\text{e},20}^{0.65}
$$

---------------

#### Hubbard 2012 L-I Lower Scaling

Is selected with `ilhthresh = 16` 

$$
P_{\text{L-H}} = 2.11 \times I_{\text{p}}^{0.7}\bar{n}_{\text{e},20}^{0.47}
$$

---------------

#### Hubbard 2012 L-I Upper Scaling

Is selected with `ilhthresh = 17` 

$$
P_{\text{L-H}} = 2.11 \times I_{\text{p}}^{1.18}\bar{n}_{\text{e},20}^{0.83}
$$

---------------

###  Hubbard 2017 L-I Scaling

Is selected with `ilhthresh = 18` 

$$
P_{\text{L-H}} = 0.162 \times B_{\text{T}}^{0.26}\bar{n}_{\text{e},20} S_{\text{p}}
$$

where $B_{\text{T}}$ is the toroidal magnetic filed in $\text{T}$, $\bar{n}_{\text{e},20}$ is the line-averaged electron density in units of $10^{20}$ and $S_{\text{p}}$ is the plasma surface area in $\text{m}^2$.

------------------

| `ilhthresh` | Name | Reference |
| :-: | - | - |
| 1 | ITER 1996 nominal | ITER Physics Design Description Document |
| 2 | ITER 1996 upper bound | D. Boucher, p.2-2 |
| 3 | ITER 1996 lower bound | 
| 4 | ITER 1997 excluding elongation | J. A. Snipes, ITER H-mode Threshold Database |
| 5 | ITER 1997 including elongation |  Working Group, Controlled Fusion and Plasma Physics, 24th EPS conference, Berchtesgaden, June 1997, vol.21A, part III, p.961 |
| 6 | Martin 2008 nominal | Martin et al, 11th IAEA Tech. Meeting |
| 7 | Martin 2008 95% upper bound |  H-mode Physics and Transport Barriers, Journal |
| 8 | Martin 2008 95% lower bound |  of Physics: Conference Series **123**, 2008 |
| 9 | Snipes 2000 nominal | J. A. Snipes and the International H-mode |
| 10| Snipes 2000 upper bound | Threshold Database Working Group |
| 11| Snipes 2000 lower bound |  2000, Plasma Phys. Control. Fusion, 42, A299 |
| 12| Snipes 2000 (closed divertor): nominal | 
| 13| Snipes 2000 (closed divertor): upper bound | 
| 14| Snipes 2000 (closed divertor): lower bound | 
| 15| Hubbard 2012 L-I threshold scaling: nominal | [Hubbard et al. (2012; Nucl. Fusion 52 114009)](https://iopscience.iop.org/article/10.1088/0029-5515/52/11/114009) |
| 16| Hubbard 2012 L-I threshold scaling: lower bound | [Hubbard et al. (2012; Nucl. Fusion 52 114009)](https://iopscience.iop.org/article/10.1088/0029-5515/52/11/114009 |
| 17| Hubbard 2012 L-I threshold scaling: upper bound | [Hubbard et al. (2012; Nucl. Fusion 52 114009)](https://iopscience.iop.org/article/10.1088/0029-5515/52/11/114009 |
| 18| Hubbard 2017 L-I threshold scaling | [Hubbard et al. (2017; Nucl. Fusion 57 126039)](https://iopscience.iop.org/article/10.1088/1741-4326/aa8570) |
| 19 | Martin 2008 aspect ratio corrected nominal | Martin et al (2008; J Phys Conf, 123, 012033) |
| 20 | Martin 2008 aspect ratio corrected 95% upper bound | [Takizuka et al. (2004; Plasma Phys. Contol. Fusion, 46, A227)](https://iopscience.iop.org/article/10.1088/0741-3335/46/5A/024)  |
| 21 | Martin 2008 aspect ratio corrected 95% lower bound |  


[^1]: T. Takizuka and International Atomic Energy Agency, Vienna (Austria),"Threshold power and energy confinement for ITER". 1996.
[^2]: J. C. Wesley, “International Thermonuclear Experimental Reactor: Physics issues, capabilities and physics program plans,” Physics of Plasmas, vol. 4, no. 7, pp. 2642-2652, Jul. 1997, doi: https://doi.org/10.1063/1.872406.