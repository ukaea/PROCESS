# Vacuum vessel

## Stress in vacuum vessel due to fast discharge of superconducting TF coil

If a superconducting TF coil undergoes a quench (rapid loss of superconductivity), a protection circuit will initiate a fast discharge of the entire TF coil set, in which much of the stored magnetic energy is discharged into an external dump resistor. The discharge is characterised by an exponential time constant, determined by the resistance of the dump resistor, but limited by the maximum permissible induced voltage. In a fast discharge a poloidal current is induced in the vessel. In the initial phase of the fast discharge a large fraction of the initial toroidal field is still present. The consequent Lorentz forces act outward on the vessel, analagous to an internal pressure.  Making the vessel stronger by increasing its overall thickness is of limited value, as this also reduces the electrical resistance.  However, the vessel can be strengthened using ribs and ports.

A model has been implemented (not yet merged as of 23/6/23), based on  
[Empirical Formulas for Estimating Self and Mutual Inductances of Toroidal Field Coils and Structures, Itoh et al](https://www.jstage.jst.go.jp/article/pfr/15/0/15_1405078/_article) (2020)

This model takes account of the currents induced in both the vacuum vessel and the steel TF coil structures.

Constraint 65 implements this model, by applying a maximum permitted stress in the vacuum vessel.  
`fmaxvvstress` f-value for constraint 65. Iteration variable 113.   
`theta1_coil` An angle shown in the figure below, relating to the shape of the TF coil conductor centre-line (degrees).   
`theta1_vv` An angle shown in the figure below, relating to the shape of the vacuum vessel centre-line (degrees).   
`max_vv_stress` The maximum permissible maximum shear stress in the vacuum vessel (Pa) (as used in the Tresca criterion)    

Example output:  
```
Minimum allowed quench time due to stress in VV (s)                      (taucq)                   2.831E+01  OP 
Actual quench time (or time constant) (s)                                (tdmptf)                  2.840E+01  ITV
Vacuum Vassel stress on quench (Pa)                                      (vv_stress_quench)        4.589E+07  OP 
```
In this example the stress is much smaller than the permissible value, so the constraint is satisfied but has no effect on the design.


The reference above proceeds as follows.

1.  An analytical formula is derived for the self-inductance of a TF coil composed of 6 circular arcs:

![image](/uploads/ea95aa6adb28e779400c873da63749b9/image.png)

2.  The parameters of the 6 arcs are expressed in terms of the overall parameters parameters κ, δ, A, H and θ1.  (Note these are parameters of the coil or vessel, not the plasma.  Also, θ1 is rather arbitrary, not a meaningful quantity.)

3.  As the analytical formula for the self-inductance is complicated, three surrogate formulae are derived by a regression analysis for the range:  
1.5 ≤ κ ≤ 2.0,  
1.5 ≤ A ≤ 2.0,   
0.22 ≤ δ ≤ 0.5,   
0 < θ1/90◦ ≤ 0.7
H sets the scale and can have any value.

All the inductance formulae are based determining the dimensionless factor ξ in equation 1.  

The aspect ratio A of the vacuum vessel will always be smaller than the aspect ratio of the plasma, but it could easily exceed the maximum figure of 2.0 used for the surrogate formula.

Surrogate 1: Equation 2:  
![image](/uploads/ce2fa0dad1e1358fd6b949c0b3dec8df/image.png)   
![image](/uploads/8b3e4f684b89d77b72905b5fd3a0f313/image.png)  
![image](/uploads/bbd3ce64e55634d7012ba8cbf934d6ce/image.png)  
(Note that ε is not actually part of the formula - it just represents the error.)  
Maximum error 1.5%.

Surrogate 2 (most accurate): Equation 2 and a correction term Δξ given by   
![image](/uploads/5dc97896ebb6afb5f114d8320e454a9f/image.png)![image](/uploads/ddf8fb25438418168f420a59762901d4/image.png)  
![image](/uploads/838e0ddcc0abd57edf80a6e98bc79426/image.png)

Surrogate 3: ξ is given by
![image](/uploads/85d704fded3509233791f3cb6c0bd81d/image.png)  
![image](/uploads/21cb9e7f0d331c9742297e564c698edf/image.png)
![image](/uploads/d210b3d0b7e2ac76b089fb052059afb5/image.png)  
maximum error 1.4%.  
This is the simplest formula and doesn't use θ1, which is useful since θ1 is rather arbitrary.

The surrogate formulae assume up/down symmetry, but the exact formulae in the appendix do not, I think.

4.  The exact and surrogate self-inductance formulae have been checked against the value for the ITER TF coil set.

5.  In a fast discharge, eddy currents are induced in the TF coil structures, consisting of the coil cases and the steel in the winding pack (conduit plus radial plates if any).  Eddy currents are also induced in the vacuum vessel.  The induced currents and voltages are related by two simultaneous circuit equations involving the inductances and resistances.

6.  Based on a conceptual design, we estimate the parameters κ, δ, A, and H for the midlines of the TF coil, the coil structure, and the vacuum vessel.

7.  The self and mutual inductances of these three structures are calculated using the exact or the surrogate formulae.  Note that the mutual inductance of the vacuum vessel and the coils is identical to the self-inductance of the vacuum vessel (multiplied by the number of turns).  The vessel is inside the coils, and with axisymmetry, the field generated by the vessel is zero outside the vessel.  All the field inside the vessel also links the coil set.  Therefore the flux linkage used for self-inductance is identical to that used for mutual inductance.

```math
M_{02} = L_2 N_0
```
where  
throughout, subscript 0 denotes the TF-coil conductor, 1 denotes the coil structure, and 2 denotes the vacuum vessel (VV).  
$` M_{02} = `$ mutual inductance of TF coils and the vacuum vessel,  
$`L_2 = `$ self-inductance of vacuum vessel,    
$`N_0 = N_{TFC} N_C`$  
$`N_{TFC}`$ is the number of TF coils,  
$`N_C`$ is the number of turns per coil.

8.  The poloidal loop resistance of the TF structure is approximately  
```math
R_1 = \frac{\eta_{SSL} l_{CCL}}{N_{TFC} S_{structure}}
```
where, giving some numerical values from the paper,    
$`\eta_{SSL} = `$ resistivity of SS316 (∼0.5 μΩm) at low temperature (4.2K),   
$` l_{CCL} = `$ the poloidal length of the midline of the TF structure,  
$`S_{structure} = `$ the cross-sectional area of the TF structure (including the case and the steel content of the winding pack).

The poloidal loop resistance of the vessel is approximately  
```math
R_2 = \frac{\eta_{SSH}}{\Delta_{VV}} \phi
```
where, giving some numerical values from the paper,    
$`\eta_{SSH} = `$ the resistivity of SS316 (0.84 μΩm) at high temperature (100◦C),  
$`\Delta_{VV} = `$ the vacuum vessel thickness, including an allowance for any poloidally continuous ribs.  
$`\phi = 0.94`$ an approximate numerical factor.

9.  The L/R time constant of the TF coil structure is assumed to be much greater than that for the vacuum vessel.  (This approximation is probably fine, but is not strictly needed.  The equations can be solved without it.)

10.  The current can now be derived as a function of time for each structure.  The TF coil structure has a significant influence on the induced current in the vacuum vessel.  The current in the TF coils is assumed to decay exponentially.

```math
I_0(t) = I_{OP} \space e^{-\lambda_0 t}
```

```math
I_1(t) = \lambda_0 N_0 I_{OP} \frac{e^{-\lambda_1 t} - e^{-\lambda_0 t}}{\lambda_0 - \lambda_1}
```

```math
I_2(t) = \frac{\lambda_1}{\lambda_2} I_1   \space \space \space \space  (assuming \lambda_2>>\lambda_0, \lambda_1)
```

where  
$` I_0 = `$ current in the TF coil conductor  
$` I_1 = `$ current in TF coil structure   
$`I_2 = `$ current in vacuum vessel  

```math
\lambda_0 = 1/\tau_d
```

```math
\lambda_1 = R_1 / L_1
```
```math
\lambda_1 = R_2 / L_2
```
$`\tau_d = `$ decay time constant for the fast discharge of the TF coils,   
$`R_1 = `$ poloidal loop resistance of TF coil structure,  
$`R_2 = `$ poloidal loop resistance of vacuum vessel,  
$`L_1 = `$ self-inductance of the TF coil structure,  
$`L_2 = `$ self-inductance of vacuum vessel,  

The force on the vessel is proportional to the product of the field and the current in the vessel.  

![image](/uploads/bfe3bb4084f50df1eef0e29277e656ed/image.png) 


The field in the inboard VV wall is due to the coils, plus the current in the structure, plus the current in the vessel itself:  
```math
B_{VVI} = μ_0 \frac{I_0 + I_1 + I_2/2} {2πR_{VVI}}
```
$`R_{VVI} = `$ major radius of the inboard VV wall.

This should be evaluated at time $`t_{maxforce}`$.

The current density in the inboard wall of the vessel is
```math
J_{VVI} = \frac{I_2}{2π Δ_{VV} R_{VVI}}
```
This should also be evaluated at time $`t_{maxforce}`$.

The stress used for the Tresco criterion (the maximum shear stress)  is then given by:

![image](/uploads/58847abe2bc0b420ef74da22bda3e752/image.png)  
$`R_{VVI} = `$ major radius of the inboard VV wall.  
$`R_{VVO} = `$ major radius of the outboard VV wall.


Obviously there are some key uncertainties.  The VV will include numerous ports that may tend to interrupt the poloidal flow of current, although the current can still flow around the port tube itself.  The thermal shield and in-vessel components have been ignored.  Only the stress in the inboard VV wall is included, and the possibility of buckling is ignored.

(New 23/6/23)
