

## Overview

The diamagnetic current fraction $f_{\text{dia}}$ is strongly related to $\beta$ and is typically small,
hence it is usually neglected.  For high $\beta$ plasmas, such as those at tight
aspect ratio, it should be included and two scalings are offered.  If the diamagnetic
current is expected to be above one per cent of the plasma current, a warning
is issued to calculate it.

$$
J_{\text{dia}} = \frac{RB_{\theta}}{B^2}\frac{dP}{d\psi}
$$

$$
\frac{I_{\text{dia}}}{I_{\text{p}}} = \frac{\beta_0}{2}\int_0^1 \frac{\rho^2q(a)}{q(\rho)}\frac{d\left(\frac{P}{P(0)}\right)}{dx} dx
$$

Where $\rho$ is the normalised radius = $\frac{r}{a}$

-----------------------------------

### No diamagnetic current 

To have it so that the diamagnetic current is not calculated you can set `i_diamagnetic_current = 0`

------------------------------------

### T.Hender fit for ST's:

This model can be used by setting: `i_diamagnetic_current = 1`


$$f_{\text{dia}} = \frac{\beta}{2.8}$$



------------------------------------

###  SCENE fit:

This model can be used by setting: `i_diamagnetic_current = 2`

This model is based off of 108 equilibria from SCENE.
Overall the equilibria cover: 

- $A$ = 1.6 to 3.2
- $\beta$ = 0.5% to 26%
- $\frac{P(0)}{\langle P \rangle}$ = 1.8 to 7.2
- $l_i$(2) = 0.21 to 1.0
- $\frac{q_{95}}{q(0)}$ = 0.9 to 16 \ \ (i.e. deeply hollow current to very peaked)

$$f_{\text{dia}} = 0.414 \space \beta \space \left(\frac{0.1 q_{95}}{q_0} + 0.44\right)$$