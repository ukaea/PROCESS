

### Bottura scaling | `bottura_scaling()`

This is a generic scaling proposed for the characterization and production of 
ITER Nb<sub>3</sub>Sn strands. This is also known as the "ITER-2008 parametrization."[^1]

The datasets span 10 years of R&D and production, with 100 to 500 data points each, covering strain (1.5% to 0.4%), temperature (2.35 to 16 K), and field (0.5 to 19 T). The ITER-2008 parameterization achieves an average accuracy error of 3.8 Amps, with the best at 1.5 Amps and the worst at 7.5 Amps. Compared to Durham University's full parameterization, ITER-2008 yields root mean squared errors 1.5 times larger, which is significant but not dramatic.

The strain function is suitable only in the moderate strain region, down to 0.8%. Beyond this, the measured behavior shows an inflection and curvature change not captured in the equations. Fitting strain data beyond an intrinsic compressive strain of 0.5% impacts the moderate strain regime and reduces fit accuracy.


$$
J_{\text{C}} = \frac{C}{B}s(\epsilon)(1-t^{1.52})(1-t^2)b^p (1-b)^q
$$

$$
B_{C2}^*(T,\epsilon) = B_{C20 \text{max}}^* s(\epsilon)(1-t^{1.52})
$$

$$
T_{C}^*(T,\epsilon) = T_{C0 \text{max}}^* \left[s(\epsilon)\right]^{\frac{1}{3}}\left(1-\frac{B}{B_{C2}^*(0,\epsilon)}\right)^{\frac{1}{1.52}}
$$

$$
s(\epsilon) = 1 + \frac{C_{a1}\left(\sqrt{\epsilon_{sh}^2+\epsilon_{0,a}^2}- \sqrt{(\epsilon -\epsilon_{sh})2+\epsilon_{0,a}^2}\right)- C_{a2}\epsilon}{1- C_{a1} \epsilon_{0,a}}
$$

$$
\epsilon_{sh} = \frac{C_{a2}\epsilon_{0,a}}{\sqrt{C_{a1}^2-C_{a2}^2}}
$$

- **\( C \)**: Scaling constant
- **\( J_{\text{C}} \)**: Critical current density, representing the maximum current density a superconductor can carry without resistance.
- **\( B \)**: Magnetic field strength at the conductor
- **\( s(\epsilon) \)**: Strain scaling function, accounting for the effect of strain on superconducting properties.
- **\( t \)**: Reduced temperature, defined as \( T / T_{C0 \text{max}}^* \), where \( T \) is the temperature at the conductor
- **\( b \)**: Reduced magnetic field, defined as \( B / B_{C2}^* \), where \( B_{C2}^* \) is the upper critical field.
- **\( p \)**: Low field exponent of the pinning force ($\approx$ 0.5)
- **\( q \)**: High field exponent of the pinning force ($\approx$ 2)
- **\( B_{C2}^*(T,\epsilon) \)**: Upper critical field as a function of temperature and strain.
- **\( B_{C20 \text{max}}^* \)**: Maximum upper critical field at zero temperature and strain.
- **\( T_{C2}^*(T,\epsilon) \)**: Critical temperature as a function of temperature and strain.
- **\( T_{C0 \text{max}}^* \)**: Maximum critical temperature at zero field and strain.
- **\( \epsilon \)**: Strain applied to the superconductor.
- **\( \epsilon_{sh} \)**: Strain shift parameter, related to the intrinsic strain behavior.
- **\( \epsilon_{0,a} \)**: Residual strain component
- **\( C_{a1}, C_{a2} \)**: Material-specific strain fitting constants used in the strain scaling function.


[^1]: L. Bottura and B. Bordini, “$J_{C}(B,T,\varepsilon)$ Parameterization for the ITER ${\rm Nb}_{3}{\rm Sn}$ Production,” IEEE Transactions on Applied Superconductivity, vol. 19, no. 3, pp. 1521-1524, Jun. 2009, doi: https://doi.org/10.1109/tasc.2009.2018278.