# Plant Availability

Switch `iavail` is used to control how the overall plant availability factor `cfactr` is calculated, as follows:

If `iavail = 0`, the input value of `cfactr` is used.

If `iavail = 1`, a model by N. Taylor and D. Ward[^1] is used instead, in which `cfactr` is calculated taking into account the time taken to replace certain components of the fusion power core, and various unplanned unavailability fractions which may be set by the user, as summerised in Table 1.

| Input parameter | description |
| --- | --- | --- |
| `tbktrepl` | time needed to replace blanket (years) |
| `tdivrepl` | time needed to replace divertor (years) |
| `tcomrepl` | time  needed to replace both blanket and divertor (years) |
| `uubop` | unplanned unavailability of balance of plant |
| `uucd` | unplanned unavailability of current drive system |
| `uudiv` | unplanned unavailability of divertor |
| `uufuel` | unplanned unavailability of fuelling system |
| `uufw` | unplanned unavailability of first wall |
| `uumag` | unplanned unavailability of magnets |
| `uuves` | unplanned unavailability of vessel |

Table 3.4: *Summary of the variables in `PROCESS` that relate to the Taylor-Ward availability model (`iavail=1`).*

If `iavail = 2`, the new Morris model is implemented[^2]. It estimates both planned and unplanned unavailability, and the time during which no power is being generated if the reactor is pulsed.

The panned unavailability is linked to the lifetime of the blanket and the time taken to replace them. The lifetime of the blanket is based on the allowable fast neutron fluence. In contrast, the lifetime of the divertor is estimated using the particle and photon load. The time to replace the blanket and divertor have been estimated by Crofts et al, who studied the influence of the number of remote handling systems working in parallel. `PROCESS` uses a simple fit to their results, and adds a month to allow the dose rate to reduce to an acceptable level before remote handling operations start, and a month to allow for pump-down and preparation for operation for operation at the enf of the shutdown.

To estimate the unplanned downtime, each subsystem is represented by a simple model that tries to capture the degradation of reliability when approaching operational and technological limits. This increases the risk of unplanned downtime as design margins are reduced.

For the toroidal field coils, the chance of a quench is likely to be the largest driver of the risk of unplanned unavailability, and this may depend on the temperature margin - the difference between the actual temperature and the critical temperature of the superconductor. The user inputs a minimum temperature margin, below which quench is immediate, and a second, higher value, above which quenches never occur. In between, there is a finite quench rate.

The unplanned downtime for the blanket is based on the number of cycles it experiences before planned replacement. (This model is restricted of pulsed reactors.) The cycle life of the blanket is expressed using a reference lifetime. Before this lifetime, there is a constant but small probability of failure in each pulse. After the reference lifetime the reliability of the blanket starts to decline, reaching zero at the upper lifetime limit.

It is assumed that the vacuum system can be maintained in parallel with blanket replacement, so it does not contribute to the planned downtime. The unplanned downtime is baed on an assumed failure rate for a cryo-pump, and a specified total number pumps, with some of them being redundant. The resulting downtime can be reduced to a negligible level if there are several redundant pumps, but in addition, there is a fixed unavailability to allow for common mode failures affecting several pumps.

[^1]: P. J. Knight, *"PROCESS 3020: Plant Availability Model"*, Work File Note
F/PL/PJK/PROCESS/CODE/<br>
[^2]: M. Kovari, F. Fox, C. Harrington, R. Kembleton, P. Knight, H. Lux, J. Morris *"PROCESS: a systems code for fusion power plants - Part 2: Engineering"*, Fus. Eng. & Des. 104, 9-20 (2016)