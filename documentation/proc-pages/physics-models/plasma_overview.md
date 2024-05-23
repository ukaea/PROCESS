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

More detail is given in [^19], but this webpage is more up to date.


## Other Plasma Physics Options

### Neo-Classical Correction Effects

Neo-classical trapped particle effects are 
included in the calculation of the plasma resistance and ohmic heating power in 
subroutine `pohm`, which is called by routine `physics`.  The scaling used is only valid for aspect 
ratios between 2.5 and 4, and it is possible for the plasma resistance to be 
incorrect or even negative if the aspect ratio is outside this range.  An error is reported if the 
calculated plasma resistance is negative.

### Inverse Quadrature in $\tau_E$ Scaling Laws

Switch `iinvqd` determines whether the energy confinement time scaling
laws due to Kaye-Goldston (`isc = 5`) and Goldston (`isc = 9`) should include 
an inverse quadrature scaling with the Neo-Alcator result (`isc = 1`). A value 
`iinvqd = 1`includes this scaling.

### Plasma-Wall Gap

The region directly outside the last closed flux surface of the core plasma is
known as the scrape-off layer, and contains no structural material. Plasma
entering this region is not confined and is removed by the divertor. PROCESS
treats the scrape-off layer merely as a gap. Switch `iscrp` determines
whether the inboard and outboard gaps should be calculated as 10% of the
plasma minor radius (`iscrp = 0`), or set equal to the input values `scrapli` 
and `scraplo` (`iscrp = 1`).

