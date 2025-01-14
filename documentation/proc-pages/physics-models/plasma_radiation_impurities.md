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

The synchrotron radiation power[^2] [^3] is assumed to originate from the 
plasma core. The wall reflection factor `ssync` may be set by the user.

By changing the input parameter `coreradius`, the user may set the normalised 
radius defining the 'core' region. Only the impurity and synchrotron radiation 
from this affects the confinement scaling. Figure 1 below shows the
radiation power contributions.

![Schematic diagram of radiation power contributions](../images/radiation.png "Schematic diagram of radiation power contributions")
*Figure 1: Schematic diagram of the radiation power contributions and how they are split between core and edge radiation*

Constraint equation no. 17 with iteration variable no. 28 (`fradpwr`)
ensures that the calculated total radiation power does not exceed the total
power available that can be converted to radiation (i.e. the sum of the fusion
alpha power, other charged particle fusion power, auxiliary injected power and
the ohmic power). This constraint should always be turned on.

[^1]: H. Lux, R. Kemp, D.J. Ward, M. Sertoli, 'Impurity radiation in DEMO 
systems modelling', Fus. Eng.  | Des. **101**, 42-51 (2015)
[^2]:  Albajar, Nuclear Fusion **41** (2001) 665
[^3]: Fidone, Giruzzi and Granata, Nuclear Fusion **41** (2001) 1755