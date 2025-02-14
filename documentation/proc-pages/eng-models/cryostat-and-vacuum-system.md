# Cryostat

The _vacuum vessel_ provides a toroidal evacuated chamber containing the plasma, first wall, blanket and shield.  The _cryostat_ is a cylindrical chamber enclosing the entire reactor, including the vacuum vessel and all the coils and the intercoil structure.  It provides a vacuum for thermal insulation.

-------------

## Cryostat geometry | `external_cryo_geometry()`

### Calculate inboard radius

The radius of the inboard side of the cryostat is found by taking the radius of the furthest out PF coil and adding a clearance gap:

$$
\mathtt{r\_cryostat\_inboard}, r_{\text{cryostat}} = \text{max}(r_{\text{PF}}) + \mathtt{dr\_pf\_cryostat}
$$

where $\mathtt{dr\_pf\_cryostat}$ is the radial PF coil to cryostat gap specified by the user at input.

----------------

### Vertical clearance

The top flange of the cryostat will be a large structure taking a considerable load from atmospheric pressure.  The vertical distance $\mathrm{d}z_{\text{PF,cryostat}}$ between the uppermost PF coil and the top inside flange of the cryostat is set using a scaling based on ITER is used:

$$
\mathtt{dz\_pf\_cryostat}, \mathrm{d}z_{\text{PF,cryostat}} = \texttt{f_z_cryostat} \left( \frac{2 \times \texttt{r_cryostat_inboard}}{28.440}\right)
$$

-------------------

### Half-height

The internal half height of the cryostat is then calculated by taking the maximum vertical height of the PF coils and adding the calculated clearance, $\mathtt{dz\_pf\_cryostat}$.

$$
\mathtt{z\_cryostat\_half\_inside} = \text{max}(z_{\text{PF}}) + \mathtt{dz\_pf\_cryostat}
$$

-------------------

### Vertical clearance of TF coil

The vertical clearance between the top of the TF coil and the inside of the cryostat is then calculated:

$$
\mathtt{dz\_tf\_cryostat} = \mathtt{z\_cryostat\_half\_inside} - (z_{\text{TF}} + \mathrm{d}z_{\text{TF}})
$$

where $z_{\text{TF}}$ is the height of the inside of the TF leg and $\mathrm{d}z_{\text{TF}}$ is its thickness.

----------------------


### Calculate cryostat volume

We calculate the cryostat volume by taking the outer dimensions of the cryostat structure and then remove that of the inside structure. This is just subtracting the volumes of two cylinders.

$$
\mathtt{vol\_cryostat}, V_{\text{cryostat}} = \\
\underbrace{\left[\pi \left(r_{\text{cryostat}}+dr_{\text{cryostat}}\right)^2 \times 2\left(\mathtt{z\_cryostat\_half\_inside}+ dr_{\text{cryostat}}\right)\right]}_{\text{Outer shell}} \\
- \underbrace{\left[\pi r_{\text{cryostat}}+^2 \times 2\left(\mathtt{z\_cryostat\_half\_inside}\right) \right]}_{\text{Inner shell}}
$$

where $dr_{\text{cryostat}}$ is the uniform thickness of the cryostat that is set at input by the user with `dr_cryostat =`

-------------------


# Vacuum pumping
The vacuum system is used for four different processes. Firstly, before plasma operations the chamber must be evacuated. Secondly, the chamber must be re-evacuated between pulses. Thirdly, helium ash must be removed throughout the burn to prevent it from diluting the fuel. Finally, deuterium and tritium escaping from the confined plasma are removed continuously. `PROCESS` calculates the parameters of a vacuum system that satisfy all four requirements, with the option of either turbo pumps or cryo pumps being used.

Switch `ntype` controls whether a turbo pump (`ntype` = 0) or a cryo pump (`ntype` = 1) is used in the vacuum system.
