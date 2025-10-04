At the moment `PROCESS` uses raw text for its input files. This guide should be used for users who wish to create their own machine configurations from scratch but dont have the knowedlege of the different models or how to build up the optimisation problem.

It is always recommended to use an pre-existing `*.IN.DAT` file as the template for your work as this removes the majority of the boilerplate in have to re-define key parameters. Example input files can be found in `examples/data/*_IN.DAT`

This guide will not be verbose in describing the model selection but will point to the reference sections of the docs for the user to decide on what models they may wish to use. 

This guide will only give guidance on how to set up a very un-restricted input file where most of the key variables are iteration variables and have wide bounds. This should make the file very easy to converge initially when ran. Further constraining the input id down to the user depending on their aea of interest.

## Plasma

### 1: Geometry

For the plasma we will depict a fixed size and aspect ratio scenario.

1. Select which plasma shape formaultion you wish to have with the [`i_plasma_shape` switch](../physics-models/plasma_geometry.md).

2. Decide on the [aspect ratio ($A$ | `aspect`)](https://euro-fusion.org/glossary/aspect-ratio/) and [major radius ($R_0$ | `rmajor`)](https://euro-fusion.org/glossary/aspect-ratio/) of the plasma you wish to have:

```
>>> IN.DAT

* PLASMA GEOMETRY

i_plasma_shape = 1
aspect = 3.0
rmajor = 8.0

```

3. Select how the elongation ($\kappa$ | `kappa`)($\kappa_{95}$ | `kappa95`) and triangularity ($\delta$, | `triang`)($\delta_{95}$, | `triang95`) is defined by setting the [`i_plasma_geometry` switch](../physics-models/plasma_geometry.md#plasma-geometry-parameters-plasma_geometry)

```
>>> IN.DAT

* PLASMA GEOMETRY

i_plasma_shape = 1
aspect = 3.0
rmajor = 8.0

* ===============================

* Sepatrix elongation and triangularity are input to find the 95% values
i_plasma_geometry = 0
kappa = 2.0
triang = 0.5

```

--------------------

### 2: Profiles

1. Select the type of plasma profile you want with the [`ipedestal` switch](../physics-models/profiles/plasma_profiles.md).

```
>>> IN.DAT

* PLASMA PROFILES

* Plasma profiles have a pedestal (H-mode like plasma)
ipedestal = 1

```

2. Set the volume averaged electron temperature ($T_{\text{e}}$ | `te`) and density ($n_{\text{e}}$ | `dene`) to be wide iteration variables.

```
>>> IN.DAT

* PLASMA PROFILES

* Plasma profiles have a pedestal (H-mode like plasma)
ipedestal = 1

* Sets the volume averaged electron temperature to be an iteration variable that 
* can be found between 7 and 20 keV
ixc = 4 
te = 10.69
boundl(4) = 7.0
boundu(4) = 20.0

* Sets the volume averaged electron density to be an iteration variable that 
* can be found between 0.7E+20 and 2.0E+20
ixc = 6
dene = 1.0E+20
boundl(6) = 0.7E+20
boundu(6) = 2.0E+20

```

3. Decide on peaked the profile should be by setting their indexes. [($\alpha_{\text{T}}$ | `alphat`) and ($\alpha_{\text{n}}$ | `alphan`)](../physics-models/profiles/plasma_profiles.md#pedestal-profile--h-mode)

>>> IN.DAT

* PLASMA PROFILES

* Plasma profiles have a pedestal (H-mode like plasma)
ipedestal = 1

* Sets the volume averaged electron temperature to be an iteration variable that 
* can be found between 7 and 20 keV
ixc = 4 
te = 10.69
boundl(4) = 7.0
boundu(4) = 20.0

* Sets the volume averaged electron density to be an iteration variable that 
* can be found between 0.7E+20 and 2.0E+20
ixc = 6
dene = 1.0E+20
boundl(6) = 0.7E+20
boundu(6) = 2.0E+20


* Peaked temperature profile
alphat = 2.0

* Less peaked and flatter density profile
alphan = 0.5


```

4. Decide on the position of the pedestals 