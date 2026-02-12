# Detailed Plasma Physics

It can sometimes be useful to calculate rough values for key plasma paramters that are normally used in higher fidelity codes. The `DetailedPhysics()` class stores functions that are called and the end of the run to show rough values for key plasma behavior parameters. The calculation is done at the end as no other methods currently depend on these values.

## Detailed Plasma Physics | `DetailedPhysics()`


------------------

### Debye length | `calculate_debye_length()`

Calculates the Debye lenght given by:

$$
\lambda_{D} = \sqrt{\frac{\epsilon_0 k_B T_e}{n e^2}}
$$

-------------------

### Relativistic particle speed | `calculate_relativistic_particle_speed()`

$$
v = c \times \sqrt{\left(1- \frac{1}{\left(1+\frac{E}{mc^2}\right)^2}\right)}
$$

------------------

### Coulomb Logarithm | `calculate_coulomb_log_from_impact()`

Calculates the Coulomb logarithm assuming a straight line Landau-Spitzer method

$$
\ln \Lambda = \ln{\left(\frac{b_{\text{max}}}{b_{\text{min}}}\right)}
$$

The maximum impact parameter is given by the Debye length calculated by [`calculate_debye_length()`](#debye-length--calculate_debye_length)
$$
b_{\text{max}} = \lambda_{\text{Debye}}
$$

The minimum impact paramter is the largest of either the classical distance of closest approach or the Debye length. 

$$
\begin{split}b_{\text{min}} ≡
\left\{
    \begin{array}{ll}
               λ_{\text{de Broglie}} & \mbox{if } λ_{\text{de Broglie}} ≥ ρ_⟂ \\
               ρ_⟂         & \mbox{if } ρ_⟂ ≥ λ_{\text{de Broglie}}
    \end{array}
\right.\end{split}
$$

$ρ_⟂$ is the classical distance of closest approach calculated by [`calculate_classical_distance_of_closest_approach()`](#classical-distance-of-closest-approach----calculate_classical_distance_of_closest_approach)

------------------

### Classical distance of closest approach |   `calculate_classical_distance_of_closest_approach()`

$$
\frac{Z_1Z_2e^2}{4\pi \epsilon_0 E_{\text{kinetic}}}
$$

---------------------

### DeBroglie Wavelength | `calculate_debroglie_wavelength()`

$$
\lambda_{\text{DeBroglie}} = \frac{h}{2\pi m v}
$$

----------------------

### Plasma Frequency | `calculate_plasma_frequency()`

$$
\omega_p = \sqrt{\frac{n_ie^2}{\epsilon_0 m_i}}
$$

---------------------

### Larmor Frequency | `calculate_larmor_frequency()`

$$
f_{\text{Larmor}} = \frac{Z_ieB}{2\pi m_i}
$$