# Pulsed Plant Operation

If the plasma current is partially or entirely driven by electromagnetic induction, it is 
necessary to operate the plant in a pulsed manner as the current swing in the central solenoid 
cannot be continued indefinitely. `PROCESS` can perform a number of calculations relevant to a 
pulsed power plant, as detailed below.

Switch `lpulse` determines whether the power plant is assumed to be based on steady-state 
(`lpulse = 0`) or pulsed (`lpulse = 1`) operation.  The current ramp calculations apply in both 
cases, as even a steady-state reactor has to be started up.

## Start-up power requirements

The auxiliary power reaching the plasma can be constrained to be more than the minimum allowable 
value `auxmin` by turning on constraint equation no. 40 with iteration variable no. 64 (`fauxmn`). 
The value of `auxmin` is set in the input file.

The auxiliary power required during the start-up and ramp-up phase is not calculated.  (The code 
contains a routine based on a POPCON analysis of access via the Cordey Pass - the path in plasma 
density-temperature space which minimises the power requirement - but this is very CPU-intensive, 
so it is not called at present.)

## Plasma current ramp-up time

When the plasma current is ramped up too quickly an unstable current distribution can develop, 
leading to kink and tearing instabilities.  PROCESS allows an upper limit to be set for the rate 
of change of plasma current. In the default case, the ramp-up time is given by setting the rate of 
change equal to the maximum proposed in [^1], or it can be set by the user.  The physics of this 
constraint is likely to depend on whether the ramp-up is purely inductive or includes current drive, 
but this is not taken ito account.

In the steady-state scenario (`lpulse` = 0), the plasma current ramp-up time `tohs` is determined as follows. 

- If `tohsin` = 0, the rate of change of plasma current is 0.5 MA/s. The PF coil ramp time `tramp` 
  and shutdown time `tqnch` are (arbitrarily) set equal to `tohs`. 
- If `tohsin` $\neq$ 0, the plasma current ramp-up time `tohs` = `tohsin`, and the PF coil ramp 
  and shutdown times are input parameters.

In the pulsed scenario, (`lpulse` = 1), the plasma current ramp-up time `tohs` is an input, and it 
can be set as an iteration variable (65). The ramp-up and shutdown time in the pulsed case are set 
equal to `tohs`. To ensure that the plasma current ramp rate during start-up is prevented from being 
too high, as governed by the requirement to maintain plasma stability by ensuring that the induced 
current has time to diffuse into the body of the plasma, constraint equation no. 41 should be 
turned on with iteration variable no. 66 `ftohs` and input `tohsmn`, the minimum plasma current 
ramp-up time.

## Burn time

The length of the burn time is calculated from the surplus volt-seconds available from the Central 
Solenoid and the other PF coils during the plasma burn phase, after the flux required during the 
plasma start-up is taken into account. A minimum burn time (`tbrnmn`) can be enforced via 
constraint equation no. 13 and iteration variable no 21 (`ftburn`).

## Thermal storage and back-up generation

During every cycle there is a period when no fusion power is produced. If the net electricity 
output from the plant must be maintained this is achieved using thermal storage or back-up 
generation. Three models are available, determined by the value of switch `istore`. If `istore = 1` 
(the default), option 1 of Ref[^1] is assumed, which utilises the thermal storage inherent in the 
machine's steam cycle equipment. This should be used is the machine down time is less than 100 
seconds. If `istore = 2` option 2 of Ref[^2] is assumed, which uses the same method as before, 
but augments it with an additional boiler. This may be used for machine down times of up to 300 
seconds. Finally, if `istore = 3`, a large stainless steel block acts as the thermal storage 
medium. These options are obsolete.

[^1]: *Pulsed Fusion Reactor Study*, AEA Fusion Report AEA FUS 205 (1992)
[^2]: *ITER poloidal field system*, IAEA, Vienna, 1990, 203–30, 56–60.
