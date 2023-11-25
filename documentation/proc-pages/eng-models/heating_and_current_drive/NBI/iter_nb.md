# ITER Neutral Beam Model

This model calculates the current drive parameters for a neutral beam system, based on the 1990 ITER model.[^1]


Firstly the beam access is checked for such that
$$
\bigg(1+ \frac{1}{A}\bigg) > (R_{\text{tangential}}/R_0)
$$

The beam path length to centre is calculated:

$$
\text{`dpath`} = R_0 \sqrt(\bigg(1+\frac{1}{A}\bigg)^2-`frbeam`^2\bigg)
$$


Beam stopping cross-section ($\sigma_{\text{beam}}$) is calculated using the `sigbeam` method described [here](nbi_overview.md) :


Calculate number of electron decay lengths to centre

$$
\tau_{\text{beam}} = \text{`dpath`} n_e \sigma_{\text{beam}}
$$

Shine-through fraction of beam:
$$
fshine = e^{(-2.0 * dpath * n_e * sigstop)} \\
fshine = max(fshine, 1e-20)
$$

Deuterium and tritium beam densities:
$$
n_D = n_i * (1.0 - current_drive_variables.ftritbm) \\
n_T = n_i * current_drive_variables.ftritbm
$$

Power split to ions / electrons
fpion = self.cfnbi(
    physics_variables.abeam,
    current_drive_variables.enbeam,
    physics_variables.ten,
    physics_variables.dene,
    dend,
    dent,
    physics_variables.zeffai,
    physics_variables.dlamie,
)

# Current drive efficiency
effnbss = current_drive_variables.frbeam * self.etanb(
    physics_variables.abeam,
    physics_variables.alphan,
    physics_variables.alphat,
    physics_variables.aspect,
    physics_variables.dene,
    current_drive_variables.enbeam,
    physics_variables.rmajor,
    physics_variables.ten,
    physics_variables.zeff,
)

return effnbss, fpion, fshine

[^1]: N. A. Uckan and ITER Physics Group, *"ITER Physics Design Guidelines: 1989"*, ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)