# Plasma Radiation

The synchrotron radiation power[^1] [^2] is assumed to originate from the 
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

[^1]:  Albajar, Nuclear Fusion **41** (2001) 665
[^2]: Fidone, Giruzzi and Granata, Nuclear Fusion **41** (2001) 1755