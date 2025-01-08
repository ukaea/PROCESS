# Plasma composition and impurities.

Within `PROCESS` we always assume the plasma is quasi-neutral eg:

$$
n_{\text{e}} = \underbrace{n_{\text{i}}}_{\text{Fuel Ions}} + \underbrace{2n_{\text{e}}f_{\alpha}}_{\text{Alpha particles}} + \underbrace{n_{\text{e}}f_{\text{beam}}}_{\text{Neutral beams}} + \underbrace{\sum_j Z_j n_{\text{e}} f_j}_{\text{Impurities}}
$$

## Plasma Composition Calculation | `plasma_composition()`

This function sets the fractional makeups and determines the relative density and charges of different plasma species.

All of the plasma composites are normally given as a fraction of the volume averaged plasma electron density.

1. **Alpha Ash Portion Calculation**

    - Calculate the number density of alpha particles (`nd_alphas`) using the electron density (`dene`) and the alpha particle to electron ratio (`ralpne`).

    `ralpne` can be set as an iteration variable (ixc = 109) or set directly.

2. **Protons Calculation**
    - The calculation of proton density (`nd_protons`) depends on whether the alpha rate density has been calculated. This should only happen in the first function call as the rates are calculated later on in the code.
    ```python
    if physics_variables.alpha_rate_density_total < 1.0e-6:
        physics_variables.nd_protons = max(
            physics_variables.f_nd_protium_electrons * physics_variables.dene,
            physics_variables.nd_alphas * (physics_variables.f_helium3 + 1.0e-3),
        )
    else:
        physics_variables.nd_protons = max(
            physics_variables.f_nd_protium_electrons * physics_variables.dene,
            physics_variables.nd_alphas
            * physics_variables.proton_rate_density
            / physics_variables.alpha_rate_density_total,
        )
    ```
    - If the alpha rate density is not yet calculated, use a rough estimate.
    - Otherwise, use the calculated proton rate density.

3. **Beam Hot Ion Component**
    ```python
    if physics_variables.ignite == 0:
        physics_variables.nd_beam_ions = (
            physics_variables.dene * physics_variables.rnbeam
        )
    else:
        physics_variables.nd_beam_ions = 0.0
    ```
    - Calculate the number density of beam ions (`nd_beam_ions`). If the plasma is ignited, set it to zero.

4. **Sum of Zi.ni for All Impurity Ions**
    ```python
    znimp = 0.0
    for imp in range(impurity_radiation_module.n_impurities):
        if impurity_radiation_module.impurity_arr_z[imp] > 2:
            znimp += impurity_radiation.zav_of_te(
                imp, np.array([physics_variables.te])
            ).squeeze() * (
                impurity_radiation_module.impurity_arr_frac[imp]
                * physics_variables.dene
            )
    ```
    - Sum the product of charge number (`Zi`) and number density (`ni`) for all impurity ions with charge greater than helium.

5. **Fuel Portion - Conserve Charge Neutrality**
    ```python
    znfuel = (
        physics_variables.dene
        - 2.0 * physics_variables.nd_alphas
        - physics_variables.nd_protons
        - physics_variables.nd_beam_ions
        - znimp
    )
    ```
    - Calculate the fuel portion (`znfuel`) by conserving charge neutrality.

6. **Fuel Ion Density Calculation**
    ```python
    physics_variables.deni = znfuel / (1.0 + physics_variables.f_helium3)
    ```
    - Calculate the fuel ion density (`deni`).

7. **Set Hydrogen and Helium Impurity Fractions**
    ```python
    impurity_radiation_module.impurity_arr_frac[
        impurity_radiation.element2index("H_")
    ] = (
        physics_variables.nd_protons
        + (physics_variables.f_deuterium + physics_variables.f_tritium)
        * physics_variables.deni
        + physics_variables.nd_beam_ions
    ) / physics_variables.dene

    impurity_radiation_module.impurity_arr_frac[
        impurity_radiation.element2index("He")
    ] = (
        physics_variables.f_helium3
        * physics_variables.deni
        / physics_variables.dene
        + physics_variables.ralpne
    )
    ```
    - Set the impurity fractions for hydrogen and helium.

8. **Total Impurity Density Calculation**
    ```python
    physics_variables.nd_impurities = 0.0
    for imp in range(impurity_radiation_module.n_impurities):
        if impurity_radiation_module.impurity_arr_z[imp] > 2:
            physics_variables.nd_impurities += (
                impurity_radiation_module.impurity_arr_frac[imp]
                * physics_variables.dene
            )
    ```
    - Calculate the total impurity density (`nd_impurities`).

9. **Total Ion Density Calculation**
    ```python
    physics_variables.nd_ions_total = (
        physics_variables.deni
        + physics_variables.nd_alphas
        + physics_variables.nd_protons
        + physics_variables.nd_beam_ions
        + physics_variables.nd_impurities
    )
    ```
    - Calculate the total ion density (`nd_ions_total`).

10. **Set Obsolescent Impurity Fraction Variables**
    ```python
    physics_variables.rncne = impurity_radiation_module.impurity_arr_frac[
        impurity_radiation.element2index("C_")
    ]
    physics_variables.rnone = impurity_radiation_module.impurity_arr_frac[
        impurity_radiation.element2index("O_")
    ]
    physics_variables.rnfene = (
        impurity_radiation_module.impurity_arr_frac[
            impurity_radiation.element2index("Fe")
        ]
        + impurity_radiation_module.impurity_arr_frac[
            impurity_radiation.element2index("Ar")
        ]
    )
    ```
    - Set some obsolescent impurity fraction variables for other routines.

11. **Effective Charge Calculation**
    ```python
    physics_variables.zeff = 0.0
    for imp in range(impurity_radiation_module.n_impurities):
        physics_variables.zeff += (
            impurity_radiation_module.impurity_arr_frac[imp]
            * impurity_radiation.zav_of_te(
                imp, np.array([physics_variables.te])
            ).squeeze()
            ** 2
        )
    ```
    - Calculate the effective charge (`zeff`).

12. **Define Coulomb Logarithm**
    ```python
    physics_variables.dlamee = (
        31.0
        - (np.log(physics_variables.dene) / 2.0)
        + np.log(physics_variables.te * 1000.0)
    )
    physics_variables.dlamie = (
        31.3
        - (np.log(physics_variables.dene) / 2.0)
        + np.log(physics_variables.te * 1000.0)
    )
    ```
    - Define the Coulomb logarithm for ion-electron and electron-electron collisions.

13. **Fraction of Alpha Energy to Ions and Electrons**
    ```python
    if physics_module.first_call == 1:
        pc = (
            (1.0 + physics_variables.alphan)
            * (1.0 + physics_variables.alphat)
            / (1.0 + physics_variables.alphan + physics_variables.alphat)
        )
        physics_module.first_call = 0
    else:
        pc = physics_variables.pcoef

    physics_variables.f_alpha_electron = 0.88155 * np.exp(
        -physics_variables.te * pc / 67.4036
    )
    physics_variables.f_alpha_ion = 1.0 - physics_variables.f_alpha_electron
    ```
    - Calculate the fraction of alpha energy going to electrons and ions.

14. **Average Atomic Masses of Injected Fuel Species**
    ```python
    physics_variables.m_fuel_amu = (
        (constants.m_deuteron_amu * physics_variables.f_deuterium)
        + (constants.m_triton_amu * physics_variables.f_tritium)
        + (constants.m_helion_amu * physics_variables.f_helium3)
    )
    ```
    - Calculate the average atomic masses of injected fuel species.

15. **Average Atomic Masses of Injected Fuel Species in Neutral Beams**
    ```python
    physics_variables.m_beam_amu = (
        constants.m_deuteron_amu * (1.0 - current_drive_variables.f_tritium_beam)
    ) + (constants.m_triton_amu * current_drive_variables.f_tritium_beam)
    ```
    - Calculate the average atomic masses of injected fuel species in the neutral beams.

16. **Density Weighted Mass Calculation**
    ```python
    physics_variables.m_ions_total_amu = (
        physics_variables.m_fuel_amu * physics_variables.deni
        + (constants.m_alpha_amu * physics_variables.nd_alphas)
        + physics_variables.nd_protons
        + physics_variables.m_beam_amu * physics_variables.nd_beam_ions
    )
    for imp in range(impurity_radiation_module.n_impurities):
        if impurity_radiation_module.impurity_arr_z[imp] > 2:
            physics_variables.m_ions_total_amu += (
                physics_variables.dene
                * impurity_radiation_module.impurity_arr_frac[imp]
                * impurity_radiation_module.impurity_arr_amass[imp]
            )

    physics_variables.m_ions_total_amu = (
        physics_variables.m_ions_total_amu / physics_variables.nd_ions_total
    )
    ```
    - Calculate the density-weighted mass of ions.

17. **Mass Weighted Plasma Effective Charge**
    ```python
    physics_variables.zeffai = (
        physics_variables.f_deuterium * physics_variables.deni / 2.0
        + physics_variables.f_tritium * physics_variables.deni / 3.0
        + 4.0 * physics_variables.f_helium3 * physics_variables.deni / 3.0
        + physics_variables.nd_alphas
        + physics_variables.nd_protons
        + (1.0 - current_drive_variables.f_tritium_beam)
        * physics_variables.nd_beam_ions
        / 2.0
        + current_drive_variables.f_tritium_beam
        * physics_variables.nd_beam_ions
        / 3.0
    ) / physics_variables.dene
    for imp in range(impurity_radiation_module.n_impurities):
        if impurity_radiation_module.impurity_arr_z[imp] > 2:
            physics_variables.zeffai += (
                impurity_radiation_module.impurity_arr_frac[imp]
                * impurity_radiation.zav_of_te(
                    imp, np.array([physics_variables.te])
                ).squeeze()
                ** 2
                / impurity_radiation_module.impurity_arr_amass[imp]
            )
    ```
    - Calculate the mass-weighted plasma effective charge (`zeffai`).








---------------

The impurity radiation model in PROCESS uses a multi-impurity model which
integrates the radiation contributions over an arbitrary choice of density and
temperature profiles[^1]

The impurity number density fractions relative to the electron density are constant and are set 
using input array `fimp(1,...,14)`. The available species are as follows:

| `fimp` | Species |
| :-: | - |
| 1 | Hydrogen isotopes (fraction calculated by code) |
| 2 | Helium (fraction calculated by code) |
| 3 | Beryllium |
| 4 | Carbon |
| 5 | Nitrogen |
| 6 | Oxygen |
| 7 | Neon |
| 8 | Silicon |
| 9 | Argon |
| 10 | Iron |
| 11 | Nickel |
| 12 | Krypton |
| 13 | Xenon |
| 14 | Tungsten |

As stated above, the number density fractions for hydrogen (all isotopes) and
helium need not be set, as they are calculated by the code to ensure 
plasma quasi-neutrality taking into account the fuel ratios
`f_deuterium`, `f_tritium` and `f_helium3`, and the alpha particle fraction `ralpne` which may 
be input by the user or selected as an iteration variable.

The impurity fraction of any one of the elements listed in array `fimp` (other than hydrogen 
isotopes and helium) may be used as an iteration variable.
The impurity fraction to be varied can be set simply with `fimp(i) = <value>`, where `i` is the corresponding number value for the desired impurity in the table above.

[^1]: H. Lux, R. Kemp, D.J. Ward, M. Sertoli, Impurity radiation in DEMO systems modelling, Fusion Engineering and Design, Volume 101, 2015, Pages 42-51, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2015.10.002.
