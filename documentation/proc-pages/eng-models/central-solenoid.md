# Central Solenoid

The central solenoid (CS) is a PF coil used during start-up and during the burn phase to create and 
maintain the plasma current by electromagnetic induction. Swinging (changing) the current through 
the central solenoid causes a change in the flux linked to the plasma region, inducing a current in 
it. `PROCESS` calculates the amount of flux required to produce the plasma curren, and also the 
amount actually available. The code measures the magnetic flux in units of Volt.seconds (= Webers).

Switch `iohcl` controls whether a central solenoid is present. A value of 1 denotes that this coil 
is present, and should be assigned a non-zero thickness `ohcth`. A value of `iohcl` = 0 denotes 
that no central solenoid is present, in which case the thickness `ohcth` should be zero. No PF 
coils should be located at positions defined by `ipfloc(j)` = 1 if no central solenoid is present.

The central solenoid can be either resistive or superconducting (controlled via switch `ipfres` as 
for the other PF coils), and if superconducting, switch `isumatpf` determines the superconducting 
material to use -  its value is used like `isumattf` and `isumatpf`. The copper fraction (by volume) 
of the superconducting strands is `fcuohsu`.

The hoop stress is calculated using equations 4.10 and 4.11 from "Superconducting magnets", Martin N. 
Wilson (1983).  This is divided by the fraction of the area occupied by steel to obtain the hoop 
stress in the steel, $\sigma_{hoop}$.

The axial stress can be calculated using "Case studies in superconducting magnets", Y. Iwasa, p. 
86, 3.5.2, Special Case 4: Midplane force.  This applies exactly only to a thin-walled solenoid. 
The axial stress in the steel is given by:

$$
\sigma_z = \frac{F_z}{f_z A_z}
$$

where $F_z$ is the axial force, $f_z$ is the fraction of the horizontal cross-section occupied by 
steel, and $A_z$ is the area of the horizontal cross-section.

The fraction of the horizontal cross-section occupied by steel is calculated assuming that the 
conductor is square and has a steel jacket with the same thickness on all four sides, giving:

$$
f_z = \frac{f_V}{2}.
$$

The radial stress is neglected. The hoop and axial stresses are combined to give the maximum shear 
stress, as required by the Tresca stress criterion:

$$
\sigma_{max shear} = max(|\sigma_{hoop} - \sigma_{z}| , |\sigma_z|, |\sigma_{hoop}|)
$$

However, the axial stress is only included if the switch `i_cs_stress` = 1.  The axial stress is 
set to zero if `i_cs_stress` = 0.  This option has no physical justification but can be used if 
there are reasons to be believe that the calculation above gives unrealistically large stresses.

## Fatigue

If the the reactor is assumed to be pulsed, the CS must be assessed against fatigue. 

A simple crack growth model based on Linear Elastic Fracture Mechanics is used to estimate the 
allowable hoop stress in the conduits. The model follows the method described in the ITER Magnet 
Structural Design Criteria, using the Paris law to model the growth of a planar elliptical crack 
across the thickness of a plate with the width and thickness of the conduit wall. The Paris law 
states that the crack growth rate follows a power law:

$$
\frac{da}{dN}=\rm{C}\Delta K^m
$$

where a is the size of the crack, N is the number of cycles, C and m are material constants, and 
$\Delta K$ is the stress intensity factor. The stress intensity factor is, in turn, a function of the crack 
geometry, the residual stress in the conduit, and the alternating tensile stress (i.e. hoop stress 
in the case of the CS coils).

!!! Info "Assumptions"
    1.  The initial defect is a planar half-elliptical surface crack, normal to the long axis of the conductor.  
    2.  Initial aspect ratio of ellipse (semi-major radius \(c) /semi-minor radius (a) = 3).  
    3.  Initial crack dimensions are input.  Defaults: a<sub>0</sub>=2, c<sub>0</sub>=6 mm.  
    4.  The coupled Paris equations for the crack dimensions are integrated using the "Life Cycle" 
        method, in which the crack dimension (either *a* or *c*) is the variable of integration.  
    6.  Only the hoop stress is taken into account.  
    7.  The stress is monotonic (since the hoop stress is always positive), and its minimum value 
        is the residual stress (input).  Default: 240 MPa.  
    8.  The mean stress is taken into account using the [Walker](https://en.wikipedia.org/wiki/Crack_growth_equation#Walker_equation) 
        modification of the Paris equation, with coefficient $\gamma$=0.436  
    9.  Failure occurs when the crack dimension a equals the conduit thickness, or dimension c reaches 
        the conductor width.  
    10. No safety factor is used for the number of cycles.  
    11. A safety factor of 2 is used for the crack size.   

<figure markdown>
![CS defect](../images/conductor_coordinates.png){ width="100%"}
<figcaption>Figure 2: Sketch of CS conductor with a planar defect.</figcaption>
</figure>

An example output follows.  Note that in this example the cycle life is *not* sufficient.

```text
 Residual hoop stress in CS Steel (Pa)                                    (residual_sig_hoop)       2.400E+08     
 Minimum burn time (s)                                                    (tbrnmn)                  7.200E+03     
 Initial vertical crack size (m)                                          (t_crack_vertical)        8.900E-04     
 Initial radial crack size (m)                                            (t_crack_radial)          2.670E-03     
 CS turn area (m)                                                         (a_oh_turn)               1.904E-03     
 CS turn length (m)                                                       (l_cond_cst)              7.557E-02     
 CS turn internal cable space radius (m)                                  (r_in_cst)                6.732E-03     
 CS turn width (m)                                                        (d_cond_cst)              2.519E-02     
 CS structural vertical thickness (m)                                     (t_structural_vertical)   5.863E-03     
 CS structural radial thickness (m)                                       (t_structural_radial)     5.863E-03     
 Allowable number of cycles till CS fracture                              (n_cycle)                 7.529E+02  OP 
 Minimum number of cycles required till CS fracture                       (n_cycle_min)             2.000E+04  OP 
```

The parameters for the Paris law are hard-coded as follows, based on the properties of stainless steel 316LN from
Sarasola et al, IEEE Transactions on Applied Superconductivity, vol. 30, no. 4, pp. 1-5, (2020):

> C = 65e-14 m
> m = 3.5  

The model has some limitations:

1. Cycle life is just an output.  There is no constraint to ensure the cycle life is sufficient.
2. The model only includes hoop stress.

The required cycle life is set in different ways depending on the following switch.

If `bkt_life_csf` = 1 then `n_cycle_min` = `bktcycles`, which is calculated using the blanket life model.
If `bkt_life_csf` = 1 then `n_cycle_min` is an input.

## Pre-compression structure

The central solenoid model in PROCESS consists of a single coil.  In practice, however, a central 
solenoid usually consists of several coils, which can have opposite currents.  This leads to vertical 
forces that tend to separate the coils.  To prevent this, ITER has "tie-plates" which hold the coil 
segments together.  PROCESS has a corresponding structure, known as the pre-compression structure, 
made up of two cylinders, one on the inside and one on the outside, of the same thickness. The 
radii of the two cylinders are `bore` and `bore` + `ohcth`.  The thickness is derived using the 
separation force and the combined cross-sectional area:

$$
p = \frac{F}{2 \pi f \sigma (2r+t) }
$$

where:
$p$ = `precomp`   CS coil precompression structure thickness (m)
$F$ = `fseppc`    Separation force
$f$ = `fcspc`     Fraction of space occupied by pre-compression structure
$\sigma$ = `sigallpc`   allowable stress in pre-compression structure (Pa)  

The central solenoid pre-compression structure is included in the model if and only if `iprecomp` = 1.

## Current density inputs and limits

The absolute value of the central solenoid current density at the end-of-flat-top ('EOF'), `coheof`, 
is specified by the user, and can be used as an iteration variable (no. 37). The current density at 
the beginning-of-pulse ('BOP' - See Figure 1) is specified as a (positive) fraction of `coheof` 
using `fcohbop` (iteration variable no. 41). The current density in the CS at all other times is 
calculated by taking into account the flux swing necessary to initiate and maintain plasma current. 

<figure markdown>
![current-vs-time-plot](../images/current_vs_time.png){ width="100%"}
<figcaption>Figure 2: Plot showing schematically the current waveforms for the plasma, a typical PF 
coil, and the central solenoid. Note that the currents in some of the PF coils may be the opposite 
sign to that shown, and the central solenoid current may remain positive during the I<sub>p</sub> 
ramp-up period, although it will pass through zero during the burn phase.</figcaption>
</figure>

The current density in the central solenoid can be limited at BOP and at EOF. To limit the current 
density at BOP, constraint equation no. 27 is used with iteration variable no. 39 (`fjohc0`). To 
limit the current density at the EOF, constraint equation no. 26 should be turned on with iteration 
variable no. 38 (`fjohc`).

The critical current density *J*<sub>crit</sub> is a function of the temperature of the superconductor. 
The temperature margin $\Delta$*T* is the difference between the current sharing temperature and the 
operating temperature.  The current sharing temperature is the temperature at which *J*<sub>crit</sub> 
is equal to the operating current density *J*<sub>op</sub>. The minimum allowed $\Delta$*T* can be 
set using input parameter `tmargmin` together with constraint equation no. 60 and iteration variable 
no. 106 (`ftmargoh`).

It is recommended that EITHER the temperature margin constraint (60), OR the current density 
constraints (26 and 27) are activated.

!!! tip "Recommended maximum current ratio"
    For engineering feasibility, the centrepost currents at end of flat-top and beginning of pulse (`fjohc` and `fjohc0` respectively) shouldn't be set above 0.7.
