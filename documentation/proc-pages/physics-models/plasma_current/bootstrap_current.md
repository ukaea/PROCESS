## Overview

Bootstrap current in tokamaks originates from the pressure gradients within the plasma and the resulting collisions between particles. As the plasma pressure varies radially, it creates a differential in particle velocities, leading to a net drift of charged particles. This drift generates a toroidal current, known as the bootstrap current, which flows parallel to the magnetic field lines. The phenomenon is a consequence of the neoclassical transport theory, where the collisional processes in a magnetically confined plasma lead to a self-sustaining current. This current is particularly advantageous as it reduces the need for external current drive systems, thereby enhancing the efficiency and stability of the tokamak operation. The bootstrap current is proportional to the pressure gradient and the collisionality of the plasma, making it a critical factor in the design and operation of advanced tokamak reactors aiming for steady-state fusion.

Some more info can be found [here](https://wiki.fusion.ciemat.es/wiki/Bootstrap_current)
## Selection

The fraction of the plasma current provided by the bootstrap effect
can be either input into the code directly, or calculated using one of four
methods, as summarised here. Note that methods `i_bootstrap_current = 1-3` do not take into account the 
existence of pedestals, whereas the Sauter et al. scaling 
(`i_bootstrap_current = 4`) allows general profiles to be used. 

| `i_bootstrap_current` | Description |
| :-: | - |
| 1 | ITER scaling -- To use the ITER scaling method for the bootstrap current fraction.  Set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$). This method is valid at high aspect ratio only.
| 2 | General scaling -- To use a more general scaling method, set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$).
| 3 | Numerically fitted scaling [^1] -- To use a numerically fitted scaling method, valid for all aspect ratios, set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$).
| 4 | Sauter, Angioni and Lin-Liu scaling [^2] [^3] -- Set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$).
| 5 | Sakai, Fujita and Okamoto scaling [^4] -- Set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$). The model includes the toroidal diamagnetic current in the calculation due to the dataset, so `idia = 0` can only be used with it

!!! Note "Fixed Bootstrap Current"
    Direct input -- To input the bootstrap current fraction directly, set `bscfmax` 
    to $(-1)$ times the required value (e.g. -0.73 sets the bootstrap faction to 0.73).

[^1]: H.R. Wilson, Nuclear Fusion **32** (1992) 257
[^2]: O. Sauter, C. Angioni and Y.R. Lin-Liu, Physics of Plasmas **6** (1999) 2834 
[^3]: O. Sauter, C. Angioni and Y.R. Lin-Liu, Physics of Plasmas **9** (2002) 5140  
[^4]: Ryosuke Sakai, Takaaki Fujita, Atsushi Okamoto, Derivation of bootstrap current fraction scaling formula for 0-D system code analysis, Fusion Engineering and Design, Volume 149, 2019, 111322, ISSN 0920-3796, https://doi.org/10.1016/j.fusengdes.2019.111322.  