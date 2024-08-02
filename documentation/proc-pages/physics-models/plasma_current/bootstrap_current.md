The fraction of the plasma current provided by the bootstrap effect
can be either input into the code directly, or calculated using one of four
methods, as summarised here. Note that methods `ibss = 1-3` do not take into account the 
existence of pedestals, whereas the Sauter et al. scaling 
(`ibss = 4`) allows general profiles to be used. 

| `ibss` | Description |
| :-: | - |
| 1 | ITER scaling -- To use the ITER scaling method for the bootstrap current fraction.  Set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$). This method is valid at high aspect ratio only.
| 2 | General scaling -- To use a more general scaling method, set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$).
| 3 | Numerically fitted scaling [^1] -- To use a numerically fitted scaling method, valid for all aspect ratios, set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$).
| 4 | Sauter, Angioni and Lin-Liu scaling [^2] [^3] -- Set `bscfmax` to the maximum required bootstrap current fraction ($\leq 1$).

!!! Note "Fixed Bootstrap Current"
    Direct input -- To input the bootstrap current fraction directly, set `bscfmax` 
    to $(-1)$ times the required value (e.g. -0.73 sets the bootstrap faction to 0.73).

[^1]: H.R. Wilson, Nuclear Fusion **32** (1992) 257
[^2]: O. Sauter, C. Angioni and Y.R. Lin-Liu, Physics of Plasmas **6** (1999) 2834 
[^3]: O. Sauter, C. Angioni and Y.R. Lin-Liu, Physics of Plasmas **9** (2002) 5140    