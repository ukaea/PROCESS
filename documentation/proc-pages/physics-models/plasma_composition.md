# Impurities and Radiation

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
