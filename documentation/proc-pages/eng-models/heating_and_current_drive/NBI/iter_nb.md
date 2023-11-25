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

        # Calculate beam path length to centre
        dpath = physics_variables.rmajor * np.sqrt(
            (1.0 + physics_variables.eps) ** 2 - current_drive_variables.frbeam**2
        )

Beam stopping cross-section is claculated using the `sigbeam` method described in :


        # Calculate beam stopping cross-section
        sigstop = self.sigbeam(
            current_drive_variables.enbeam / physics_variables.abeam,
            physics_variables.te,
            physics_variables.dene,
            physics_variables.ralpne,
            physics_variables.rncne,
            physics_variables.rnone,
            physics_variables.rnfene,
        )

        # Calculate number of decay lengths to centre
        current_drive_variables.taubeam = dpath * physics_variables.dene * sigstop

        # Shine-through fraction of beam
        fshine = np.exp(-2.0 * dpath * physics_variables.dene * sigstop)
        fshine = max(fshine, 1e-20)

        # Deuterium and tritium beam densities
        dend = physics_variables.deni * (1.0 - current_drive_variables.ftritbm)
        dent = physics_variables.deni * current_drive_variables.ftritbm

        # Power split to ions / electrons
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