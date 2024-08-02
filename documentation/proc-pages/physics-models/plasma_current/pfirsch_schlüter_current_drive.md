A similar scaling is available for the Pfirsch-Schlüter current fraction $f_{PS}$.  This is
typically smaller than the diamagnetic current, but is negative.

`ips = 0` Pfirsch-Schlüter current fraction is set to zero.

`ips = 1` Pfirsch-Schlüter current fraction is calculated using a SCENE fit for all aspect ratios:

$$ f_{PS} = -0.09 \beta $$

There is no ability to input the diamagnetic and Pfirsch-Schlüter current
directly.  In this case, it is recommended to turn off these two scalings 
and to use the method of fixing the bootstrap current fraction.