# Plasma composition and impurities.

Within `PROCESS` we always assume the plasma is quasi-neutral eg:

$$
n_{\text{e}} = \underbrace{Z_{\text{fuel}}n_{\text{i}}}_{\text{Fuel Ions}} + \underbrace{2n_{\text{e}}f_{\alpha}}_{\text{Alpha particles}} + \underbrace{n_{\text{e}}f_{\text{beam}}}_{\text{Neutral beams}} + \underbrace{\sum_j Z_j n_{\text{e}} f_j}_{\text{Impurities}}
$$

* Since only deuterium and tritium can be placed into the beams the charge coefficient on the $n_{\text{e}}f_{\text{beam}}$ is just 1.

--------------------

## Setting the impurity composition



The impurity number density fractions relative to the electron density are constant and are set 
using input array `f_nd_impurity_electrons(1,...,14)`. The available species along with their `f_nd_impurity_electrons()` index and iteration variable number are as follows:


- `f_nd_impurity_electrons(1)`: Hydrogen isotopes (fraction calculated by code)
- `f_nd_impurity_electrons(2)`: Helium (fraction calculated by code)
- `f_nd_impurity_electrons(3)`: Beryllium, (`ixc = 125`)
- `f_nd_impurity_electrons(4)`: Carbon, (`ixc = 126`)
- `f_nd_impurity_electrons(5)`: Nitrogen, (`ixc = 127`)
- `f_nd_impurity_electrons(6)`: Oxygen, (`ixc = 128`)
- `f_nd_impurity_electrons(7)`: Neon, (`ixc = 129`)
- `f_nd_impurity_electrons(8)`: Silicon, (`ixc = 130`)
- `f_nd_impurity_electrons(9)`: Argon, (`ixc = 131`)
- `f_nd_impurity_electrons(10)`: Iron, (`ixc = 132`)
- `f_nd_impurity_electrons(11)`: Nickel, (`ixc = 133`)
- `f_nd_impurity_electrons(12)`: Krypton, (`ixc = 134`)
- `f_nd_impurity_electrons(13)`: Xenon, (`ixc = 135`)
- `f_nd_impurity_electrons(14)`: Tungsten, (`ixc = 136`)

As stated above, the number density fractions for hydrogen (all isotopes) and
helium should not be set, as they are calculated by the code. This is to ensure 
plasma quasi-neutrality taking into account the fuel ratios
`f_plasma_fuel_deuterium`, `f_plasma_fuel_tritium` and `f_plasma_fuel_helium3`, and the alpha particle fraction `f_nd_alpha_electron` which may 
be input by the user or selected as an iteration variable.

!!! note "Location of impurities"

    All species/impurities are currently assumed to be distributed homogeneously throughout the plasma.
    This is an important assumptions that affects the [plasma radiation](./plasma_radiation.md).

The impurity fraction of any one of the elements listed in array `f_nd_impurity_electrons` (other than hydrogen 
isotopes and helium) may be used as an iteration variable.

**The impurity fraction to be varied can be set simply with `f_nd_impurity_electrons(i) = <value>`, where `i` is the corresponding number value for the desired impurity in the table above.**

----------------

## Plasma Composition Calculation | `plasma_composition()`

This function sets the fractional makeups and determines the relative density and charges of different plasma species.

All of the plasma composites are normally given as a fraction of the volume averaged plasma electron density.

1. **Alpha Ash Portion Calculation**

    - Calculate the number density of alpha particles (`nd_plasma_alphas_vol_avg`) using the electron density (`nd_plasma_electrons_vol_avg`) and the alpha particle to electron ratio (`f_nd_alpha_electron`).
        - `f_nd_alpha_electron` can be set as an iteration variable (`ixc = 109`) or set directly.

    $$
    n_{\alpha} = \mathtt{f\_nd\_alpha\_electron}\times n_{\text{e}}
    $$


2. **Protons Calculation**
    - The calculation of proton density (`nd_plasma_protons_vol_avg`) depends on whether the alpha rate density has been calculated. This should only happen in the first function call as the rates are calculated later on in the code.

    - If the alpha rate density is not yet calculated, use a rough estimate.

    $$
    \texttt{nd_plasma_protons_vol_avg} | n_{\text{p}} = \\
    \text{max}\left[\texttt{f_nd_protium_electrons} \times n_{\text{e}}, n_{\alpha}\times \left(f_{\text{3He}} + 0.001\right)\right]
    $$

    - Otherwise, use the calculated proton rate density.

    $$
    \texttt{nd_plasma_protons_vol_avg} | n_{\text{p}} = \\
    \text{max}\left[\texttt{f_nd_protium_electrons} \times n_{\text{e}}, n_{\alpha}\times \frac{r_{\text{p}}}{r_{\alpha,\text{total}}}\right]
    $$

    where $r_{\text{p}}$ is the rate of proton production and $r_{\alpha,\text{total}}$ is the rate of total alpha particle production, which includes beam fusion (if present).

3. **Beam Hot Ion Component**

    - Calculate the number density of beam ions (`nd_beam_ions`), using the electron density (`nd_plasma_electrons_vol_avg`) and the beam ion to electron ratio (`f_nd_beam_electron`). If the plasma is ignited, set it to zero.
        - `f_nd_beam_electron` can be set as an iteration variable (`ixc = 7`) or set directly.

    $$
    \mathtt{nd\_beam\_ions} | n_{\text{beam}} = \mathtt{f\_nd\_beam\_electron} \times n_{\text{e}}
    $$

4. **Sum of charge number density for all impurity ions**

    - Sum the product of charge number (`Zi`) and number density for all impurity ions with charge greater than helium. 

    $$
    \mathtt{znimp} = \sum_j Z_j n_{\text{e}} f_j
    $$

5. **Fuel Portion - Conserve Charge Neutrality**
    - Calculate the fuel portion (`znfuel`) by conserving charge neutrality.

    $$
   \mathtt{znfuel} | \underbrace{Z_{\text{fuel}}n_{\text{i}}}_{\text{Fuel Ions}} =  n_{\text{e}} - \underbrace{2n_{\text{e}}f_{\alpha}}_{\text{Alpha particles}} - \underbrace{n_{\text{e}}f_{\text{beam}}}_{\text{Neutral beams}} - \underbrace{\sum_j Z_j n_{\text{e}} f_j}_{\text{Impurities}}
    $$

6. **Fuel Ion Density Calculation**

    $$
    \mathtt{nd\_fuel\_ions} | n_{\text{i}} = \frac{\mathtt{znfuel}}{1+f_{\text{3He}}}
    $$

    - Calculate the fuel ion density (`nd_plasma_fuel_ions_vol_avg`).

7. **Set Hydrogen and Helium Impurity Fractions**

    - Set the impurity fractions for hydrogen and helium species.

    $$
    \frac{n_{\text{H}}}{n_{\text{e}}} = n_{\text{protons}} + \left(f_{\text{deuterium}}+f_{\text{tritium}}\right)n_{\text{i}}+ n_{\text{beam}}
    $$

    $$
    \frac{n_{\text{He}}}{n_{\text{e}}} = f_{\text{3He}}n_{\text{i}}+\mathtt{f\_nd\_alpha\_electron}
    $$

8. **Total Impurity Density Calculation**

    - Calculate the total impurity density (`nd_plasma_impurities_vol_avg`).

    $$
    \mathtt{nd\_impurities} | n_{\text{impurities}} = \sum_j n_{\text{e}} f_j
    $$

9. **Total Ion Density Calculation**

    - Calculate the total ion density (`nd_plasma_ions_total_vol_avg`).

    $$
    \mathtt{nd\_ions\_total} | n_{\text{i,total}} = n_{\text{i}} + n_{\alpha}+n_{\text{protons}}+ n_{\text{beam}}+n_{\text{impurities}}
    $$

10. **Set Impurity Fraction Variables**

    - Set global impurity fraction variables for other routines.

    $$
    \mathtt{f_nd_plasma_carbon_electron} = \frac{n_{\text{C}}}{n_{\text{e}}}, \quad \mathtt{f_nd_plasma_oxygen_electron} = \frac{n_{\text{O}}}{n_{\text{e}}} \quad \mathtt{f_nd_plasma_iron_argon_electron} = \frac{n_{\text{Fe}}+ n_{\text{Ar}}}{n_{\text{e}}}
    $$

    The variable above are set as global physics variables to be used in the [`sigbeam()`](../eng-models/heating_and_current_drive/NBI/nbi_overview.md#beam-stopping-cross-section--sigbeam) routine, which calculates the beam stopping cross-section.

11. **Effective Charge Calculation**

    - Calculate the effective charge (`zeff`), which is given by:

    $$
    Z_{\text{eff}} = \frac{\sum_j Z^2_j n_{\text{i}}}{\sum_j Z_j n_{\text{i}}}
    $$

    As we assume the plasma is quai-neutral then:

    $$
    \sum_j Z_j n_{\text{i}} = n_{\text{e}}
    $$

    Thus we can write:

    $$
    \mathtt{zeff} | Z_{\text{eff}} = \frac{\sum_j Z^2_j n_{\text{e}} f_j}{n_{\text{e}}}
    $$

    More info can be found [here](https://wiki.fusion.ciemat.es/wiki/Effective_charge_state).


12. **Fraction of Alpha Energy to Ions and Electrons**

    - Calculate the fraction of alpha energy going to electrons and ions.

    $$
    \mathtt{f\_alpha\_electron} = 0.88155\exp{\left[-\langle T_{\text{e}} \rangle \frac{(1+\alpha_n)(1+\alpha_T)}{67.4036(1+\alpha_T+\alpha_n)}\right]}
    $$

    $$
    \mathtt{f\_alpha\_ion} = (1.0 - \mathtt{f\_alpha\_electron})
    $$

13. **Average Atomic Masses of Injected Fuel Species**

    - Calculate the average atomic masses of injected fuel species.

    $$
    \mathtt{m\_fuel\_amu} | m_{\text{fuel,amu}} = \left(m_{\text{D}}f_{\text{D}}\right) + \left(m_{\text{T}}f_{\text{T}}\right) + \left(m_{\text{3He}}f_{\text{3He}}\right)
    $$

14. **Average Atomic Masses of Injected Fuel Species in Neutral Beams**

    - Calculate the average atomic masses of injected fuel species in the neutral beams.

    $$
    \mathtt{m\_beam\_amu} | m_{\text{beam,amu}} = \left(m_{\text{D}}(1-f_{\text{T,beam}}\right)) + \left(m_{\text{T}}f_{\text{T,beam}}\right)
    $$

15. **Density Weighted Mass Calculation**

    - Calculate the density-weighted mass of ions.

    $$
    \mathtt{m\_ions\_total\_amu} = \\
    \left[\frac{\left(m_{\text{fuel}}n_{\text{i}}\right) + \left(m_{\alpha}n_{\alpha}\right) + \left(m_{\text{p}}n_{p}\right) + \left(m_{\text{beam}}n_{\text{beam}}\right) + \sum_j m_{j} n_{\text{e}} f_j}{n_{\text{i,total}}}\right]
    $$

16. **Mass Weighted Plasma Effective Charge**

    - Calculate the mass-weighted plasma effective charge (`zeffai`). Similar to the calculation of the effective charge except each element is divided by its mass

    $$
    \mathtt{zeffai} | Z_{\text{eff,m}} = \frac{\sum_j \frac{Z^2_j n_{\text{e}} f_j}{m_{\text{j}}}}{n_{\text{e}}}
    $$

---------------


## Key Constraints


### Reinke criterion, divertor impurity fraction lower limit

This constraint can be activated by stating `icc = 78` in the input file.

The minimum impurity fraction required from the Reinke module can be set with, `fzmin`

The scaling value `freinke` can be varied also.


[^1]: H. Lux, R. Kemp, D.J. Ward, M. Sertoli, Impurity radiation in DEMO systems modelling, Fusion Engineering and Design, Volume 101, 2015, Pages 42-51, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2015.10.002.
