# Auxiliary Power Systems

## Current Drive

The use of inductive current drive leads to pulsed plant operation because of the limited flux swing that can be achieved using the central solenoid. This poses problems due to the fact that fatigue failures may result, and there would be a need for thermal storage to maintain output of electricity between pulses, and supply power for starting a new pulse.However, the plasma current can also be produced and maintained (partially or wholly) using non-inductive means which, in principle, removes this restriction. `PROCESS` contains a number of auxiliary current drive schemes, including various RF methods (Lower Hybrid, Electron Cyclotron,Electron Bernstein Wave, and Ion Cyclotron (Fast Wave) current drives) and also Neutral Beam current drive systems. The code calculates the efficiency and the resulting power requirements of the chosen system.

The fraction of the required plasma current to be produced by non-inductive means, `fvsbrnni`, should be set, and flag `irfcd` should be set to 0 for purely inductive scenarios, or 1 otherwise. The current drive efficiency model to be used in this latter case is defined by the value of switch `iefrf`:

- `iefrf` = 1: Fenstermacher Lower Hybrid model
- `iefrf` = 2: Ion cyclotron model[^1],
- `iefrf` = 3: Fernstermacher electron cyclotron resonance model
- `iefrf` = 4: Ehst Lower Hybrid model
- `iefrf` = 5: ITER neutral beam model[^1],[^2],
- `iefrf` = 6: Culham Lower Hybrid model[^2],
- `iefrf` = 7: Culham electron cyclotron model[^2],
- `iefrf` = 8: Culham neutral beam model[^2],
- `iefrf` = 9: Oscillating Field current drive (RFPs only - OBSOLETE-REMOVED),
- `iefrf` = 10: ECRH user input gamma,
- `iefrf` = 11: ECRH "HARE" model [^3],
- `iefrf` = 12: EBW user scaling input[^4],
- `iefrf` = 13: ECRH O-mode cutoff with Zeff and Te [^4], 


(Note that, at present, the neutral beam models do not include the effect of an edge transport barrier (pedestal) in the plasma profile.)

It is sometimes useful to adjust artificially the current drive efficiency values produced by these routines. This can be achieved by setting teh scaling coefficients `feffcd`. The wall plug to plasma efficiencies can also be adjusted, by changing the relevant variable (`etaech`, `etalh`, `etanbi` or `etaof`).

## Plasma heating

In addition to current drive, some auxiliary power can be used to heat the plasma. The value of input parameters `pheat` determines the amount of auxiliary *heating* power (in Watts) to be applied to the plasma. This variable may be used as an iteration variable (no. 11).

## Neutral beam access

If present, a neutral beam injection system needs sufficient space between the TF coils to be able to intercept the plasma tangentially. The major radius `rtanbeam` at which the centre-line of the beam is tangential to the toroidal direction is user-defined using input parameter `frbeam`, which is the ratio of `rtanbeam` to the plasma major radius `rmajor`. The maximum possible tangency radius `rtanmax` is determined by the geometry of the TF coils - see Figure 1, and this can be enforced using constraint equation no. 20 with iteration variable no. 33 (`fportsz`). The thickness of the beam duct walls may be set using input parameter `nbshield`.

<figure>
    <center>
    <img src="../../images/portsize.png" alt="NBI port" 
    title="Neutral beam access geometry" 
    width="550" height="100" />
    <br><br>
    <figcaption><i>Figure 1: Top-down schematic of the neutral beam access geometry. The beam with the maximum possible tangency radius is shown here.
    </i></figcaption>
    <br>
    </center>
</figure>

## Neutral beam losses

Input parameter `forbitloss` can be used to specify the fraction of the net injected neutral beam power that is lost between the beam particles' ionisation and thermalisation (known as the first orbit loss). This quantity cannot easily be calculated as it depends on the field ripple and other three-dimensional effects. The power lost is assumed to be absorbed by the first wall.

The power in the beam atoms that are not ionised as they pass through the plasma (shine-through) is calculated by the code. There are two constraint equations that can be used to control the beam penetration and deposition, as follows:

- It is necessary to use a beam energy that simultaneously gives adequate penetration of the beam to the centre of the plasma and tolerable shine-through of the beam on the wall after the beam has traversed the plasma. The number of exponential decay lengths, $\tau$, for the beam power to fall before it reaches the plasma centre should be in the region of ~ 4-6[^2],. Constraint equation no. 14 may be used to force $\tau$ to be equal to the value given by input parameter `tbeamin`, and is therefore in effect a beam energy consistency equation.
- Alternatively, constraint equation no. 59 with iteration variable no. 105 (`fnbshineef`) may be used to ensure that the beam power fraction emerging from the plasma is no more than the value given by input parameter `nbshinefmax`.

It is recommended that <b>only one</b> of these two constraint equations is used during a run.

## Ignited plasma

Switch `ignite` can be used to denote whether the plasma is ignited, i.e. fully self-sustaining without the need for any injected auxiliary power during the burn. If `ignite` = 1, the calculated injected power does not contribute to the plasma power balance, although the cost of the auxiliary power system is taken into account (the system is then assumed to be required to provide heating etc during the plasma start-up phase only - use `pheat` to indicate the power requirement). If `ignite` = 0, the plasma is not ignited, and the auxiliary power is taken into account in the plasma power balance during the burn phase. Also, constraint equation no. 28 can be turned on to enforce the fusion gain *Q* to be at least `bigqmin`.

[^1]: N. A. Uckan and ITER Physics Group, *"ITER Physics Design Guidelines: 1989"*, ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)

[^2]: T. C. Hender, M. K. Bevir, M. Cox, R. J. Hastie, P. J. Knight, C. N. Lashmore-Davies, B. Lloyd, G. P. Maddison, A. W. Morris, M. R. O'Brien, M.F. Turner abd H. R. Wilson, *"Physics Assessment for the European Reactor Study"*, AEA Fusion Report AEA FUS 172 (1992)

[^3]: E. Poli, M. Müller, H. Zohm, M. Kovari, *"Fast evaluation of the current driven by electron cyclotron waves for reactor studies"*, Physics of Plasmas 1 December 2018; 25 (12): 122501

[^4]: Laqua, H & Maassberg, H & Marushchenko, Nikolai & Volpe, Francesco & Weller, A & Kasparek, W. (2003). *"Electron-Bernstein-Wave Current Drive in an Overdense Plasma at the Wendelstein 7-AS Stellarator"*, Physical review letters. 90. 075003. 10.1103/PhysRevLett.90.075003. 