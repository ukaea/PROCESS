## Overview

A similar SCENE scaling like that for the diamagnetic current is available for the Pfirsch-Schlüter current fraction $f_{\text{PS}}$.  This is
typically smaller than the diamagnetic current, but is negative.

There is no ability to input the diamagnetic and Pfirsch-Schlüter current
directly.  In this case, it is recommended to turn off these two scalings 
and to use the method of fixing the bootstrap current fraction.

--------------

### No Pfirsch-Schlüter current 

To have it so that the Pfirsch-Schlüter current is not calculated you can set `i_pfirsch_schluter_current = 0`

------------------

### SCENE Fit:

This model can be used by setting: `i_pfirsch_schluter_current = 2`

This model is based off of 108 equilibria from SCENE.
Overall the equilibria cover: 

- $A$ = 1.6 to 3.2
- $\beta$ = 0.5% to 26%
- $\frac{P(0)}{\langle P \rangle}$ = 1.8 to 7.2
- $l_i$(2) = 0.21 to 1.0
- $\frac{q_{95}}{q(0)}$ = 0.9 to 16 \ \ (i.e. deeply hollow current to very peaked)

This fit is derived from the same dataset used to derive the diamagnetic current fraction for  [`i_diamagnetic_current = 2`](diamagnetic_current.md#scene-fit)

$$ f_{\text{PS}} = -0.09 \beta $$

