# Cryostat
The _vacuum vessel_ provides a toroidal evacuated chamber containing the plasma, first wall, blanket and shield.  The _cryostat_ is a cylindrical chamber enclosing the entire reactor, including the vacuum vessel and all the coils and the intercoil structure.  It provides a vacuum for thermal insulation. 

The top flange of the cryostat will be a large structure taking a considerable load from atmospheric pressure.  PROCESS does not calculate the required thickness, but the vertical distance *h* between the uppermost PF coil and the top flange of the cryostat is set using `clhsf` (default value 4.268 m); a scaling based on ITER is used:
$$
h = \mathtt{clhsf} \left( \frac{2 \times \mathtt{rdewex}}{28.440}\right)
$$

# Cryogenics
The model for the cryogenic cooling power, and the electric power to provide this, is based on D.S. Slack, J.A. Kern, J.R., Miller, Cryogenic system design for a compact tokamak reactor, UCRL-98733, DE89 003176 (1989).  See related issues for comments.

Heat conduction through the gravity support is based on these assumptions:  

| Average thermal conductivity of stainless steel between 300 K and 4.5 K | 10 W/(mK) |
|-------------------------------------------------------------------------|-----------|
| Stress in gravity support                                               | 67 MPa    |
| Length of gravity support                                               | 1 m       |

The power balance for cryogenics is detailed as in the example below.  The calculation of nuclear heating in the coils is selected using switch `inuclear`.  Only the magnet coils are included - no allowance is made for cryopumps.  Resistive current leads are assumed.

``` 
 ************************************************* Cryogenics *************************************************
 
 Conduction and radiation heat loads on cryogenic components (MW)         (qss/1.0D6)               3.246E-02  OP 
 Nuclear heating of cryogenic components (MW)                             (qnuc/1.0D6)              1.292E-02  OP 
 Nuclear heating of cryogenic components is a user input.
 AC losses in cryogenic components (MW)                                   (qac/1.0D6)               3.225E-03  OP 
 Resistive losses in current leads (MW)                                   (qcl/1.0D6)               2.065E-02  OP 
 45% allowance for heat loads in transfer lines, storage tanks etc (MW)   (qmisc/1.0D6)             3.116E-02  OP 
 Sum = Total heat removal at cryogenic temperatures (W)                   (helpow/1.0D6)            1.004E-01  OP 
 Temperature of cryogenic components (K)                                  (temp_tf_coil_cryo)                  4.500E+00     
 Efficiency (figure of merit) of cryogenic plant is 13% of ideal Carnot v                           2.028E-03  OP 
 Electric power for cryogenic plant (MW)                                  (p_cryo_plant)                  4.952E+01  OP 
```

# Vacuum pumping
The vacuum system is used for four different processes. Firstly, before plasma operations the chamber must be evacuated. Secondly, the chamber must be re-evacuated between pulses. Thirdly, helium ash must be removed throughout the burn to prevent it from diluting the fuel. Finally, deuterium and tritium escaping from the confined plasma are removed continuously. `PROCESS` calculates the parameters of a vacuum system that satisfy all four requirements, with the option of either turbo pumps or cryo pumps being used.

Switch `ntype` controls whether a turbo pump (`ntype` = 0) or a cryo pump (`ntype` = 1) is used in the vacuum system.
