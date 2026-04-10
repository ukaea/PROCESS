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

--------------------

### Upper Hybrid Frequency | `calculate_upper_hybrid_frequency()`

$$
f_{\text{UH}} = \sqrt{\omega_{p}^2+ \omega_{\text{Larmor}}^2}
$$

--------------------

### Reduced mass of two particles | `calculate_reduced_mass()`

$$
m_{\text{reduced}} = \frac{m_1 m_2}{m_1+m_2}
$$

-----------------

### Relative average velocity | `calculate_average_relative_velocity()`

$$
v_{\text{rel}} = \sqrt{v_1^2+v_2^2}
$$

-----------------

### Electron-electron collision time | `calculate_electron_electron_collision_time()`

For $T_\text{e}$ in eV

$$
\tau_{\text{ee}} = \frac{12\sqrt{2}\pi^{\frac{3}{2}}\epsilon_0^2 \sqrt{m_{\text{e}}}T_{\text{e}}^{\frac{3}{2}}}{\ln \Lambda_{\text{ee}}\text{e}^4 n_{\text{e}}}
$$

-------------------

### Electron-ion collision time | `calculate_electron_ion_collision_time()`

For $T_\text{e}$ in eV

$$
\tau_{\text{ei}} = \frac{12 \pi^{\frac{3}{2}}\epsilon_0^2 \sqrt{m_{\text{e}}}T_{\text{e}}^{\frac{3}{2}}}{\sqrt{2} Z_i^2 \ln \Lambda_{\text{ei}} \text{e}^4 n_{\text{e}}}
$$

--------------------

### Spitzer ion slowing down time | `calculate_spitzer_ion_slowing_down_time()`

For $T_\text{e}$ in eV

$$
\tau_{\text{spitzer}} = \frac{3 (2\pi)^{\frac{3}{2}}\epsilon_0^2 m_{\text{i}} T_{\text{e}}^{\frac{3}{2}}}{n_{\text{e}} Z_i^2 \ln \Lambda_{\text{ei}} \text{e}^4 n_{\text{e}}}
$$

--------------------

### Spitzer resistivity | `calculate_spitzer_resistivity()`

For $T_\text{e}$ in eV

$$
\eta_{\text{spitzer}} = \frac{4 \sqrt{2 \pi}}{3 }\frac{Z_i e^2 \sqrt{m_{\text{e}}} \ln \Lambda_{\text{ei}}}{\left(4 \pi \epsilon_0 \right)^2 T_{\text{e}}^{\frac{3}{2}}}
$$
