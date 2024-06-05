# Plant Availability

Switch `iavail` is used to control how the overall plant availability factor `cfactr` is calculated, as follows:

If `iavail = 0`, the input value of `cfactr` is used.

If `iavail = 1`, a model by N. Taylor and D. Ward[^1] is used instead, in which `cfactr` is calculated taking into account the time taken to replace certain components of the fusion power core, and various unplanned unavailability fractions which may be set by the user, as summerised in Table 1.

| Input parameter | Description |
| :-: | - |
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

If `iavail = 3`, the availability model for Spherical Tokamaks (ST) is implemented. 

!!! Warning "Warning"
		Currently, this model only uses the centrepost to calculate the availability of an ST plant. Other systems/components will be added in the future.

This model takes the user-specified time to replace a centrepost `tmain` and the centrepost lifetime `cplife` (calculated, see below) and calculates the number of maintenance cycles

$$ t_{\text{main}} + t_{\text{CP,life}} = t_{\text{maint cycle}}. $$

The number of maintenance cycles over the lifetime of the plant is calculated and then the ceiling of this value is taken as the number of centreposts required over the lifetime of the plant

$$ n_{\text{cycles}} = t_{\text{life}} / t_{\text{maint cycle}}, $$

$$ n_{\text{CP}} = \lceil n_{\text{cycles}} \rceil. $$

The planned unavailability is then what percent of a maintenance cycle is taken up by the user-specified maintenance time

$$ U_{\text{planned}} = t_{\text{main}} / t_{\text{maint cycle}} $$

and the total operational time is given by

$$ t_{\text{op}} = t_{\text{life}} (1 - U_{\text{planned}}). $$

The total availability of the plant is then given by

$$ A_{\text{tot}} = 1 - (U_{\text{planned}} + U_{\text{unplanned}} + U_{\text{planned}}U_{\text{unplanned}}) $$

where $U_{unplanned}$ is unplanned unavailability which is provided by the user i.e. how often do you expect the centrepost to break over its lifetime. The cross term takes account of overlap between planned and unplanned unavailability.

Finally, the capcity factor is given by

$$ C = A_{\text{tot}} (t_{\text{burn}} / t_{\text{cycle}}) $$

where $t_{\text{burn}}$ is the burn time and $t_{\text{cycle}}$ is the full cycle time.

## Centrepost Lifetime

All availability models in PROCESS require the calculation of the centerpost lifetime, which is detailed here.

!!! Note "Note" 
		The centrepost lifetime is calculated in full-power years (FPY).

For superconducting magnets (`i_tf_sup = 1`), the centrepost lifetime is calculated as

$$ t_{\text{CP,life}} = min(f_{\text{TF,max}}/(\phi_{\text{CP,max}}t_{\text{year}}),t_{\text{life}}) $$

where $f_{\text{TF,max}}$ is the max fast neutron fluence on the TF coil ($\mathrm{m}^{-2} \mathrm{s}$), $\phi_{\text{CP,max}}$ is the centrepost TF fast neutron flux ($\mathrm{m}^{-2}$ $\mathrm{s}^{-1}$) and $t_{\text{year}}$ is the number of seconds in a year. 

For copper or cryogenic aluminium magnets (`i_tf_sup = 0 or 2`), the centrepost lifetime is

$$ t_{\text{CP,life}} = min(f_{\text{CP, allowable}}/P_{\text{wall}}, t_{\text{life}}) $$

where $f_{\text{CP, allowable}}$ is the allowable centrepost neutron fluence and $P_{\text{wall}}$ is the average neutron wall load ($\mathrm{MW} \mathrm{m}^{-2}$).

[^1]: P. J. Knight, *"PROCESS 3020: Plant Availability Model"*, Work File Note
F/PL/PJK/PROCESS/CODE/<br>
[^2]: M. Kovari, F. Fox, C. Harrington, R. Kembleton, P. Knight, H. Lux, J. Morris *"PROCESS: a systems code for fusion power plants - Part 2: Engineering"*, Fus. Eng. & Des. 104, 9-20 (2016)