# Auxiliary Heating & Current Drive Systems | `CurrentDrive`

## Current Drive

The use of inductive current drive leads to pulsed plant operation because of the limited flux swing that can be achieved using the central solenoid. This poses problems due to the fact that fatigue failures may result, and there would be a need for thermal storage to maintain output of electricity between pulses, and supply power for starting a new pulse.However, the plasma current can also be produced and maintained (partially or wholly) using non-inductive means which, in principle, removes this restriction. `PROCESS` contains a number of auxiliary current drive schemes, including various RF methods (Lower Hybrid, Electron Cyclotron, Electron Bernstein Wave, and Ion Cyclotron (Fast Wave) current drives) and also Neutral Beam current drive systems. The code calculates the efficiency and the resulting power requirements of the chosen system.

The fraction of the required plasma current to be produced by non-inductive means, `f_c_plasma_non_inductive`, should be set, and flag `i_hcd_calculations` should be set to 0 for purely inductive scenarios, or 1 otherwise. The current drive efficiency model to be used in this latter case is defined by the value of switch `i_hcd_primary`:

- `i_hcd_primary` = 1: [Fenstermacher Lower Hybrid model](RF/fenstermacher_lower_hybrid.md)
- `i_hcd_primary` = 2: [Ion cyclotron model](RF/ic_model.md)[^1],
- `i_hcd_primary` = 3: [Fenstermacher electron cyclotron resonance model](RF/fenstermacher_electron_cyclotron_resonance.md)
- `i_hcd_primary` = 4: [Ehst Lower Hybrid model](RF/ehst_lower_hybrid.md)
- `i_hcd_primary` = 5: [ITER neutral beam model](NBI/iter_nb.md)[^1],[^2],
- `i_hcd_primary` = 6: [Culham Lower Hybrid model](RF/culham_lower_hybrid.md)[^2],
- `i_hcd_primary` = 7: [Culham electron cyclotron model](RF/culham_electron_cyclotron.md)[^2],
- `i_hcd_primary` = 8: [Culham neutral beam model](NBI/culham_nb.md)[^2],
- `i_hcd_primary` = 9: Oscillating Field current drive :warning: (OBSOLETE-REMOVED),
- `i_hcd_primary` = 10: [ECRH user input gamma](RF/ecrh_gamma.md),
- `i_hcd_primary` = 11: ECRH "HARE" model [^3] :warning: (OBSOLETE-REMOVED),
- `i_hcd_primary` = 12: [EBW user scaling input. Scaling](RF/ebw_freethy.md) (S. Freethy)
- `i_hcd_primary` = 13: [ECRH O-mode cutoff with $Z_{\text{eff}}$ and $T_{\text{e}}$](RF/cutoff_ecrh.md) (S. Freethy) [^4],

!!! Warning "Warning" 
    At present, the neutral beam models do not include the effect of an edge transport barrier (pedestal) in the plasma profile.

It is sometimes useful to adjust artificially the current drive efficiency values produced by these routines. This can be achieved by setting the scaling coefficients `feffcd`. The wall plug to plasma efficiencies can also be adjusted, by changing the relevant variable (`eta_ecrh_injector_wall_plug`, `eta_lowhyb_injector_wall_plug`, `eta_beam_injector_wall_plug` or `etaof`).

### Power limits
The maximum amount of desired heating and current drive power can be set with `p_hcd_injected_max`. This limit can be enforced by activating constraint equation 30 (`icc=30`).
Similarly the lower bound on required heating and current drive power can be set with `p_hcd_injected_min_mw`. This limit can be enforced by activating constraint equation 40 (`icc=40`).

### Secondary current drive

It is possible to have more than one type of heating and current drive system in `PROCESS`. This can be enabled by setting the `i_hcd_secondary` switch to the desired current drive scheme, following the same numbered selection for `i_hcd_primary`.
The power injected by the secondary current drive scheme has to be set to a fixed value. This value can be set with the `p_hcd_secondary_injected_mw` variable.

## Plasma heating only

In addition to current drive, some auxiliary power can be used to only heat the plasma. The value of input parameters `p_hcd_primary_extra_heat_mw` determines the amount of auxiliary heating power (in MW) to be applied to the plasma. This variable may be used as an iteration variable (`ixc = 11`).

### Secondary heating

Like for a current drive and heating system a fixed amount of heating power that does not drive current can be set with the `p_hcd_secondary_extra_heat_mw` variable.

## Ignited plasma

Switch `i_plasma_ignited` can be used to denote whether the plasma is ignited, i.e. fully self-sustaining without the need for any injected auxiliary power during the burn. If `i_plasma_ignited` = 1, the calculated injected power does not contribute to the plasma power balance, although the cost of the auxiliary power system is taken into account (the system is then assumed to be required to provide heating etc during the plasma start-up phase only - use `p_hcd_primary_extra_heat_mw` to indicate the power requirement). If `i_plasma_ignited` = 0, the plasma is not ignited, and the auxiliary power is taken into account in the plasma power balance during the burn phase. Also, constraint equation 28 (`icc = 28`) can be turned on to enforce the fusion gain *Q* to be at least `big_q_plasma_min`.

[^1]: N. A. Uckan and ITER Physics Group, *"ITER Physics Design Guidelines: 1989"*, ITER Documentation Series, No. 10, IAEA/ITER/DS/10 (1990)

[^2]: T. C. Hender, M. K. Bevir, M. Cox, R. J. Hastie, P. J. Knight, C. N. Lashmore-Davies, B. Lloyd, G. P. Maddison, A. W. Morris, M. R. O'Brien, M.F. Turner abd H. R. Wilson, *"Physics Assessment for the European Reactor Study"*, AEA Fusion Report AEA FUS 172 (1992)

[^3]: E. Poli, M. Müller, H. Zohm, M. Kovari, *"Fast evaluation of the current driven by electron cyclotron waves for reactor studies"*, Physics of Plasmas 1 December 2018; 25 (12): 122501

[^4]: Laqua, H & Maassberg, H & Marushchenko, Nikolai & Volpe, Francesco & Weller, A & Kasparek, W. (2003). *"Electron-Bernstein-Wave Current Drive in an Overdense Plasma at the Wendelstein 7-AS Stellarator"*, Physical review letters. 90. 075003. 10.1103/PhysRevLett.90.075003.