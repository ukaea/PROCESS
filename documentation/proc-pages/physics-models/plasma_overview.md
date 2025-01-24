# Plasma physics

!!! warning " Documentation re-write"

    The plasma documenation is slowly going through an update to ensure each section is properly decribed and up to date. Each section of the plasma now has its own dedicated page

## Introduction

By default, the plasma is assumed to have an up-down asymmetric, single null
configuration (although this can be changed with user inputs). A great number 
of physics models are coded within PROCESS to describe the behaviour of the 
plasma parameters such as its current, temperature, density, pressure, 
confinement etc., and also the various limits that define the stable operating 
domain. 

More detail is given in [^1], but this webpage is more up to date.


## Other Plasma Physics Options

### Neo-Classical Correction Effects

Neo-classical trapped particle effects are 
included in the calculation of the plasma resistance and ohmic heating power in 
subroutine `plasma_ohmic_heating()`, which is called by routine `physics`.  The scaling used is only valid for aspect 
ratios between 2.5 and 4, and it is possible for the plasma resistance to be 
incorrect or even negative if the aspect ratio is outside this range.  An error is reported if the 
calculated plasma resistance is negative.

### Inverse Quadrature in $\tau_E$ Scaling Laws

Switch `iinvqd` determines whether the energy confinement time scaling
laws due to Kaye-Goldston (`i_confinement_time = 5`) and Goldston (`i_confinement_time = 9`) should include 
an inverse quadrature scaling with the Neo-Alcator result (`i_confinement_time = 1`). A value 
`iinvqd = 1`includes this scaling.

[^1]: M. Kovari, R. Kemp, H. Lux, P. Knight, J. Morris, D.J. Ward, '“PROCESS”: A systems code for fusion power plants—Part 1: Physics' Fusion Engineering and Design 89 (2014) 3054–3069

