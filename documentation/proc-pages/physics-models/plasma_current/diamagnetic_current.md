The diamagnetic current fraction $f_{dia}$ is strongly related to $\beta$ and is typically small,
hence it is usually neglected.  For high $\beta$ plasmas, such as those at tight
aspect ratio, it should be included and two scalings are offered.  If the diamagnetic
current is expected to be above one per cent of the plasma current, a warning
is issued to calculate it.

`i_diamagnetic_current = 0` Diamagnetic current fraction is zero.

`i_diamagnetic_current = 1` Diamagnetic current fraction is calculated using a fit to spherical tokamak calculations by Tim Hender:

$$f_{dia} = \frac{\beta}{2.8}$$

`i_diamagnetic_current = 2` Diamagnetic current fraction is calculated using a SCENE fit for all aspect ratios:

$$f_{dia} = 0.414 \space \beta \space (\frac{0.1 q_{95}}{q_0} + 0.44)$$