Currently in `PROCESS` the inductive current fraction from the CS is not calculated directly but is just equal to ($1 - \mathtt{fvsbrnni}$). Where $\mathtt{fvsbrnni}$ is the sum of the fractions of current driven by non inductive means.

This calculated fraction (`inductive_current_fraction`) is then used in the `vscalc()` and `burn()` functions to calculate the volt-second requirements and the burn time for a pulsed machine.

!!! info "Inductive plasma current fraction refactor"

    It is hoped for the near future to have a more engineering based calculation of the fraction of the plasma current driven by the central solenoid. This would hopefully allow the setting of required ramp and flat-top times for which a given inductive current fraction can be given based on the operational performance and margin in the central solenoid.