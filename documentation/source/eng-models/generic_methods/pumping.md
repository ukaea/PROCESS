# Pumping Methods

## Pumping coolant friction | `darcy_friction_haaland()`

 The pressure drop is based on the Darcy fraction factor, using the [Haaland equation](https://en.wikipedia.org/wiki/Darcy_friction_factor_formulae#Haaland_equation), an approximation to the implicit Colebrook–White equation. 

$$
\frac{1}{\sqrt{f}} = -1.8 \log{\left[ \left(\frac{\epsilon / D}{3.7}\right)^{1.11} \frac{6.9}{\text{Re}} \right]}
$$