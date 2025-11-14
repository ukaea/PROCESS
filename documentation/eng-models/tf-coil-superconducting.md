# Superconducting TF coil

## Superconducting coil geometry
The TF coils are assumed to be supporting each other against the net centering force.  This can be described as a vaulted or wedged design. Each coil, illustrated in <em>Figure 1</em>, can be separated in two main sections:

- **The winding pack (WP)** : section containing the superconducting cables (Blue area  in <em>Figure 1</em>). The ground insulation and the insertion gap (the clearance required to allow the winding pack to be inserted, shown as the dark grey area in <em>Figure 1</em>) is considered part of the WP by convention.
- **The steel casing**: Section holding the WP providing the necessary structural support (light grey area in <em>Figure 1</em>).

<p style='text-align: justify;'> The next sub-section describes the different parametrization proposed in <em>PROCESS</em>:</p>

<figure>
    <center>
    <img src="../../images/tfcoil_SC_inboard_structure.png" alt="SC_mid_plane_cross_section" 
    title="Inboard mid-plane superconducting TF coil cross-section" 
    width="800" height="200" />
    <br>
    <figcaption><i><p style='text-align: justify;'> 
      Figure 1: Inboard mid-plane cross-section of the superconducting TF coils considered in PROCESS with associated parametrization. The green variables are calculated by the machine_build module subroutine while the blue ones are specific to the sctfcoil module. The light grey area corresponds to the steel casing containing the Winding pack and providing structural support. The dark grey area corresponds to the winding pack insertion gap, and its ground insulation. Finally the light blue area corresponds to the winding pack containing the conductor cables. More details on the parametrization is discussed in this section.</p>
    </i></figcaption>
    <br>
    </center>
</figure>

### TF coil inboard radial size

<p style='text-align: justify;'> 
  Following the geometry and its parametrization presented in <em>Figure 1</em>, the TF total thickness <em>dr_tf_inboard</em> \( \left( \Delta R_\mathrm{TF} \right) \) is related with the inner and outer case radial thicknesses (<em>dr_tf_nose_case</em>, \(  \Delta R_\mathrm{case}^\mathrm{in} \) and <em>dr_tf_plasma_case</em>, \( \Delta R_\mathrm{case}^\mathrm{out} \) respectively) and the WP radial thickness <em>dr_tf_wp_with_insulation</em> \(\Delta R_\mathrm{WP}\) by the following equation :
</p>

$$ 
\Delta R_\mathrm{TF} = R_\mathrm{TF}^\mathrm{in} + \Delta R_\mathrm{WP} + \Delta R_\mathrm{case}^\mathrm{out} + \Delta R_\mathrm{case}^\mathrm{in} - R_\mathrm{TF}^\mathrm{in}
$$

<p style='text-align: justify;'> 
  with \( R_\mathrm{TF}^\mathrm{in} \) (<em>r_tf_inboard_in</em>) the radius of the innermost TF edge, set by the central  solenoid coil size and \( N_\mathrm{TF} \) the number of TF coils. Reverted, to provide the WP thickness, the same equation simply becomes:
</p>

$$  
\Delta R_\mathrm{WP} = \left( R_\mathrm{TF}^\mathrm{in} + \Delta R_\mathrm{TF} \right) - R_\mathrm{TF}^\mathrm{in} - \Delta R_\mathrm{case}^\mathrm{out} - \Delta R_\mathrm{case}^\mathrm{in}
$$

<p style='text-align: justify;'> 
  The TF coil radial thickness (<em>dr_tf_inboard</em>) can parametrized in two ways in <em>PROCESS</em>:
</p>
- <p style='text-align: justify;'> 
    **Direct parametrization**: the TF radial inboard thickness width is set as an input variable : `dr_tf_inboard` (iteration variable 13). The WP radial thickness (`dr_tf_wp_with_insulation`) is calculated from `dr_tf_inboard` and the two case radial thicknesses. This parametrization is used by default.
  </p>
- <p style='text-align: justify;'> 
    **WP thickness parametrization**: the TF inboard radial thickness is calculated from the the case and the WP radial thickness. This option is selected by using the WP thickness (`dr_tf_wp_with_insulation`, iteration variable 140) as an iteration variable. Doing so, any `dr_tf_inboard` values will be overwritten and for this reason `dr_tf_wp_with_insulation` and `dr_tf_inboard` cannot be used as iteration variables simultaneously. Although not set by default for backward compatibility, this parametrization provides a more stable optimization procedure (negative WP area layer cannot be obtained by construction) and is hence encouraged.
  </p>

### Case geometry

<p style='text-align: justify;'> 
  Although not physically divided into pieces, three sections of the case can be considered: 
</p>
- <p style='text-align: justify;'> 
    **The nose casing:** this section corresponds to the case separating the WP with the machine center. Due to the presence of net electromechanical centering forces, this case has a major structural purpose and is often much larger than the other sides. The nose case dimension is set by its radial thickness that the user can specify using the `dr_tf_nose_case`  input variable (iteration variable 57).
  </p>
- <p style='text-align: justify;'> 
    **Sidewall casing:** this section corresponds to the lateral side of the case, separating the WP with the other vaulted coils. As in the WP geometry is generally squared, the sidewall case thickness may vary with the machine radius. For this reason, the user sets its dimensions though its minimal thickness `dx_tf_side_case_min`. The user can either directly specify `dx_tf_side_case_min` or define it as a fraction of the total coil thickness at the inner radius of the WP (`r_tf_wp_inboard_inner`) with the `casths_fraction` input. If `casths_fraction` is set in the input file, the `dx_tf_side_case_min` value will be overwritten.
  </p>
- <p style='text-align: justify;'> 
    **Plasma side casing:** this section corresponds to the case section separating the WP with the plasma. As the geometry of this section is rounded, its thickness is set by its minimal value `dr_tf_plasma_case` (user input). This parameter can also be defined as a fraction of the total TF coil thickness `dr_tf_inboard` using `f_dr_tf_plasma_case`. If the `f_dr_tf_plasma_case` parametrization is used, the `dr_tf_plasma_case` value will be overwritten.
  </p>



### Winding pack geometry

Several Winding pack geometries can chosen with the `i_tf_wp_geom` integer switch as shown in Figure 3: 

- <p style='text-align: justify;'>
    `i_tf_wp_geom = 0` : Rectangular winding pack. It is the only geometry compatible with the integer turn parametrization (`i_tf_turns_integer = 1`).
</p>
- <p style='text-align: justify;'>
    `i_tf_wp_geom = 1` : Double rectangle winding pack. The two rectangles are have the same radial thickness and their width in the toroidal direction is defined with the minimal sidewall thickness at their innermost radius.
  </p>
- <p style='text-align: justify;'> 
    `i_tf_wp_geom = 2` : Trapezoidal WP. The WP area is defined with a trapezoid, keeping the sidewall case thickness constant. This is however probably not a realistic shape as the turns are generally rectangular. This option has been added mostly to allow comparison with simplified FEA analysis configurations.
  </p>

<figure>
    <center>
    <img src="../../images/tfcoil_SC_WP_option.png" alt="WP geometries options"
    title="Superconducting TF coil WP geometry"
    width="650"
    height="100"/>
    <br><br>
    <figcaption><i><p style='text-align: justify;'> 
      Figure 3: Visual illustration of the WP shapes the user can select with the i_tf_wp_geom integer switch. The dx_tf_wp_primary_toroidal and dx_tf_wp_secondary_toroidal parameters, added in option i_tf_wp_geom = 1 are calculated using the minimal sidewall case thickness.
    </p></i></figcaption>
    <br>
    </center>
</figure>

------------------

### Turns geometry

Based on the choice of superconductor used in the coil (`i_tf_sc_mat`) the required turn structure also has to change. NbTi and Nb3Sn superconductors are normally made in cable form and are used in the Cable in Conduit turn configuration. On the other hand high temperature superconductors are normally made in tape format.

The selection of the turn geometry is controlled via the `i_tf_turn_type` switch:

- 1: Cable in conduit conductor
- 2: CroCo
- 3: STEP Vertical stacked tapes

------------------------

#### Cable in conduit conductor

<p style='text-align: justify;'>
  <em>Figure 4</em> illustrates the winding pack internal structure and the
  individual turns structure.
</p>

<figure>
    <center>
    <img src="../../images/tfcoil_SC_turn.png" alt="SC turn geometry" 
    title="Superconducting turn geometry"
    width="650"
    height="100"/>
    <br><br>
    <figcaption><i><p style='text-align: justify;'> 
      Figure 4: Illustration of the winding pack internal structure. The top 
      right diagram shows the inboard mid-plane cross section of a TF coil with
      the steel case in light grey and the winding pack in light blue. The dotted
      lines illustrate the rectangular conductor. The bottom left zoom-in section shows the conductor structure used in the
      <em>PROCESS</em> module. The red section represent the individual turn
      electrical insulation, the turquoise the steel jacket/conduit providing
      structural support, the grey is the area allocated to the superconductor
      material, copper stabiliser and inter-strand coolant, and the white circle the helium cooling channel.  (The superconductor,
      copper stabiliser and coolant together make up the 'cable'.)
    </p></i></figcaption>
    <br>
    </center>
</figure>

The winding pack is assumed to be made of \(N_\mathrm{turn} \) (`n_tf_coil_turns`) 
turns. The number of turns can be parametrized in three different ways :

- <p style='text-align: justify;'>
    **Current per turn parametrization (default):** `i_tf_turns_integer = 0` the
    user sets the value of the current flowing in each turns `c_tf_turn`. The number
    of turns necessary to carry the total TF coil current is then deduced from
    `c_tf_turn`. There is no guarantee that a realistic turn configuration (with all
    the turn geometrically fitting in the allocated space) or even have an
    integer number of turn is used with this parametrization. If the turn
    thickness `dx_tf_turn_general` or the cable space thickness `dx_tf_turn_cable_space_general` is defined by
    the user, this parametrization is not selected.
  </p>   
- <p style='text-align: justify;'>
    **Turn size parametrization:** the dimension of the turn `dx_tf_turn_general` can be
    set by the user. To do so, the user just have to select the following option:
    `i_tf_turns_integer = 0` and to set a value to the variable `dx_tf_turn_general`. The
    area of the corresponding squared turn and the number of turns necessary to
    fill the WP area is deduced. There is no guarantee that a realistic turn
    configuration (with all the turn geometrically fitting in the allocated
    space) or even have an integer number of turns is used with this parametrization.
    The current per turn `c_tf_turn` will be overwitten.
  </p>
- <p style='text-align: justify;'>
    **Cable size parametrization:** the dimension of the SC cable `dx_tf_turn_cable_space_general`
    can be set by the user. To do so, the user just have to select the following
    option: `i_tf_turns_integer = 0` and to set a value to the variable
    `dx_tf_turn_cable_space_general`. The area of the corresponding squared turn is deduced adding
    the steel conduit structure and the turn insulation. The number of turns
    necessary to fill the WP area is then deduced. There is no guarantee that a
    realistic turn configuration (with all the turn geometrically fitting in the
    allocated space) or even have an integer number of turns is used with this
    parametrization. The current per turn `c_tf_turn` will be overwitten.
  </p> 
- <p style='text-align: justify;'> 
    **Integer turn parametrization:** `i_tf_turns_integer = 1` the user sets the
    number of layers in the radial direction (`n_tf_wp_layers`) and the number of turns in the toroidal direction
    (`n_tf_wp_pancakes`). The number of turns is integer. The turn cross-section is not necessarily square, giving different
    averaged structural properties in the radial and toroidal directions. Only a rectangular WP can be used for this parametrization.
  </p>

<p style='text-align: justify;'>  
  The turn internal structure, illustrated in <em>Figure 4</em>, is inspired
  from the cable-in-conduit-conductor (CICC) design, with the main different
  being that a rounded squared cable space is used (grey area in <em>Figure 4
  </em>). The rounding curve radius is take as 0.75 of the steel conduit 
  thickness. The turn geometry is set with with the following thicknesses:
</p>
- <p style='text-align: justify;'>
    **Turn insulation thickness `dx_tf_turn_insulation`:** user input setting the thickness
    of the inter-turn insulation.
  </p>
- <p style='text-align: justify;'>
    **Steel jacket/conduit thickness `dx_tf_turn_steel` (iteration variable 58):** user
    input thickness of the turn steel structures. As it is a crucial variable
    for the TF coil structural properties it is also an iteration variable.
  </p>
- <p style='text-align: justify;'>
    **Helium cooling channel diameter `dia_tf_turn_coolant_channel`:** user input defining the 
    size of the cooling channel.
  </p>


##### Cable composition
<p style='text-align: justify;'>
  As the conductor cable composition is only used to correct the area used to
  compute current density flowing in the superconductor material, to be compared
  with its critical current density, an average material description is enough
  for the <em>PROCESS</em>models. The composition is set with the following
  material fractions:
</p>
- <p style='text-align: justify;'>
    **Cable void fraction (`f_a_tf_turn_cable_space_extra_void`):** user input setting the void fraction
    between the strands. This fraction does not include the helium cooling
    pipe at the cable center.
  </p>
- <p style='text-align: justify;'> 
    **Copper fraction (`f_a_tf_turn_cable_copper`):** user input setting the copper fraction.
    This fraction is applied after the void and helium cooling channels areas
    has been removed from the conductor area. Does not include any copper from
    REBCO tape if used.
  </p>

-----------------

#### STEP Vertical Stacked Tapes



-----------------

## Critical current density for the superconductor 
The minimum conductor cross-section is derived from the critical current density for the superconductor in the operating magnetic field and temperature, and is enforced using constraint 33.

Switch `i_tf_sc_mat` specifies which superconducting material is to be used:

- `i_tf_sc_mat == 1` -- Nb$_3$Sn superconductor, ITER critical surface parameterization[^5], standard critical values
- `i_tf_sc_mat == 2` -- Bi-2212 high temperature superconductor
- `i_tf_sc_mat == 3` -- NbTi superconductor
- `i_tf_sc_mat == 4` -- Nb$_3$Sn superconductor, ITER critical surface parameterization[^5], user-defined critical parameters
- `i_tf_sc_mat == 5` -- WST Nb$_3$Sn parameterization
- `i_tf_sc_mat == 6` -- REBCO HTS tape in CroCo strand
- `i_tf_sc_mat == 7` -- Durham Ginzburg-Landau critical surface model for Nb-Ti
- `i_tf_sc_mat == 8` -- Durham Ginzburg-Landau critical surface model for REBCO
- `i_tf_sc_mat == 9` -- Hazelton experimental data combined with Zhai conceptual model for REBCO

The fraction of copper present in the superconducting filaments is given by `f_a_tf_turn_cable_copper` (iteration variable number 59). For cases where REBCO tape is used this copper fraction does not include the copper within the tape.

For `i_tf_sc_mat = 2`, a technology adjustment factor `fhts` may be used to modify 
the critical current density fit for the Bi-2212 superconductor, to describe the 
level of technology assumed (i.e. to account for stress, fatigue, radiation, 
AC losses, joints or manufacturing variations). The default value for `fhts` is 
0.5 (a value of 1.0 would be very optimistic).

For `i_tf_sc_mat = 4`, important superconductor properties may be input as follows:
- Upper critical field at zero temperature and strain: `bcritsc`,
- Critical temperature at zero field and strain: `tcritsc`.

The toroidal field falls off at a rate $1/R$, with the peak value occurring at the outer edge of the inboard portion of the TF coil winding pack (radius `r_b_tf_inboard_peak`). 

Three constraints are relevant to the operating current density $J_{\mbox{op}}$ in the TF coils.

- Criticial current (`constraint 33`): $f_{\text{iooic}}J_{\mbox{op}}$ must not exceed the critical value $J_{\mbox{crit}}$ where `fiooic` is a margin on the constraint that defaults to `0.7`.

- Temperature margin (`constraint 36`) -- The critical current density $J_{\mbox{crit}}$ falls with 
  the temperature of the superconductor. The temperature margin $\Delta T$ is the difference between the current sharing temperature (at which $J_{\mbox{crit}}$ would be equal to $J_{\mbox{op}}$) and the operating temperature. The minimum allowed $\Delta T$
can be set using `tmargmin` together with constraint equation 36. Note that if the temperature margin is positive, $J_{\mbox{op}}$ is guaranteed to be lower than \jcrit, and so constraints 33 and 36 need not both be turned on. It is recommended that only one of these two constraints is activated.  

---------

## Quench protection

Superconducting quenches in the TF coil are modelled using a simple 0-D adiabatic heat balance in which the heat rise in the copper during a quench is equal to the heat rise required to increase the material temperature by $dT$. This is known as a ``hotspot criterion'' model, as the maximum temperature during a quench is the termination criterion (the point at which the quench is considered irreparably damaging).

The copper is assumed to carry the full current during a quench, and the materials are assumed to be in thermal equilibrium. We have:

$$
P(t)dt = \sum_i c_{P,i}\rho_{i} V_{i} dT
$$

Where $P(t)$ is the power deposited in the copper through Joule heating at time t, and the sum on the RHS is over the conductor consitutents (copper, helium, superconductor), where $V_{i}$ is the volume of said constituent. This can be rewritten as:

$$
\int J(t)^2 dt = \int_{T_{0}}^{T_{max}} \sum_i f_{i}\dfrac{c_{P,i}\rho_{i}}{\nu_{Cu}} dT
$$

The LHS is solved analytically, assuming that the current is at the operational value $J_{op}$ for some quench detection time $t_{detection}$ and then decays exponentially with a time constant $\tau_{discharge}$. The RHS is solved numerically by integration of the temperature-dependent material properties over $dT$.

The resistivity of copper, $\nu_{Cu}$, is an extremely important parameter in this model. The residual-resistance-ratio (RRR) of the copper can be specified, and magneto-resistive and irradiation-induced increases in resistivity are accounted for. The magnetic field at which the copper resisitivity is calculated is kept constant at $B_{TF,peak}$. This is conservative, and to simplify the implementation (keeping the separation of the $t$ and $T$ integrations).

Formally this gives:

$$
\begin{align}
	J_{TF,\mathrm{quench}}^{2} = \frac{1}{\tau_{TF,\mathrm{discharge}}/2 + t_\mathrm{detection}} f_{Cu} \bigg(
	& f_{He} \int_{T_0}^{T_\mathrm{max}} \frac{\rho_{He} C_{p,He}}{\nu_{Cu}}\, dT \nonumber \\
	& + f_{Cu} \int_{T_0}^{T_\mathrm{max}} \frac{\rho_{Cu} C_{p,Cu}}{\nu_{Cu}}\, dT \nonumber \\
	& + f_{sc} \int_{T_0}^{T_\mathrm{max}} \frac{\rho_{sc} C_{p,sc}}{\nu_{Cu}}\, dT \bigg)
\end{align}
$$


- `Constraint 35` -- To ensure that $J_{\mbox{op}}$ does not exceed the quench protection current density limit, $J_{TF,\mathrm{quench}}$, turn on constraint equation no.\ 35.

-----------------------


## Supercoducting TF coil class | `SuperconductingTFCoil(TFCoil)`

----------------

### Winding Pack Geometry | `superconducting_tf_wp_geometry()`

Depending on the value of `i_tf_wp_geom` different WP geometries will be configured.

Initial general dimensions are calculated first as follows:

$$
\overbrace{R_{\text{TF-inboard,WP-inner}}}^{\texttt{r_tf_wp_inboard_inner}} = \overbrace{R_{\text{TF-inboard,in}}}^{\texttt{r_tf_inboard_in}} + \overbrace{\mathrm{d}R_{\text{TF,nose-case}}}^{\texttt{dr_tf_nose_case}}
$$

$$
\overbrace{R_{\text{TF-inboard,WP-outer}}}^{\texttt{r_tf_wp_inboard_outer}} = R_{\text{TF-inboard,WP-inner}} + \overbrace{\mathrm{d}R_{\text{TF,WP}}}^{\texttt{dr_tf_wp_with_insulation}}
$$

$$
\overbrace{R_{\text{TF-inboard,WP-centre}}}^{\texttt{r_tf_wp_inboard_centre}} = \frac{R_{\text{TF-inboard,WP-inner}} + R_{\text{TF-inboard,WP-outer}}}{2}
$$

Find the straight toroidal width of the TF coil at the inside edge of the winding pack:

$$
\mathrm{d}x_{\text{TF,toroidal,WP}} = 2 \times \overbrace{R_{\text{TF,WP-inner}}}^{\texttt{r_tf_wp_inboard_inner}} \times \overbrace{\tan{\left(\phi_{\text{TF,half}}\right)}}^{\texttt{tan_theta_coil}}
$$

To find the straight toroidal length of the winding pack we now take off the side case thicknesses:

$$
\overbrace{\mathrm{d}x_{\text{TF,toroidal,WP-min}}}^{\texttt{dx_tf_wp_toroidal_min}} = \mathrm{d}x_{\text{TF,toroidal,WP}} - \left(2 \times \overbrace{\mathrm{d}x_{\text{TF,side case}}}^{\texttt{dx_tf_side_case_min}}\right)
$$

We also set a commonly used parameter for the radial thickness of the winding pack without the insulation and insertion gap.

$$
\overbrace{\mathrm{d}R_{\text{TF-WP,no-insulation}}}^{\texttt{dr_tf_wp_no_insulation}} = \overbrace{\mathrm{d}R_{\text{TF,WP}}}^{\texttt{dr_tf_wp_with_insulation}} - \\
2\times\left(\overbrace{\mathrm{d}R_{\text{TF,WP-insulation}}}^{\texttt{dx_tf_wp_insulation}} + \overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)
$$

----------

#### Rectangular WP


For a [rectangular winding pack](#winding-pack-geometry) (`i_tf_wp_geom == 0`) of constant shape:

$$
\overbrace{\mathrm{d}x_{\text{TF,toroidal,WP-primary}}}^{\texttt{dx_tf_wp_primary_toroidal}}, \overbrace{\mathrm{d}x_{\text{TF,toroidal,WP-average}}}^{\texttt{dx_tf_wp_toroidal_average}} = \mathrm{d}x_{\text{TF,toroidal,WP-min}}
$$

The full winding pack area with insulation is:

$$
\overbrace{A_{\text{TF,WP}}}^{\texttt{a_tf_wp_with_insulation}} = \mathrm{d}x_{\text{TF,toroidal,WP-primary}} \times \overbrace{\mathrm{d}R_{\text{TF,WP}}}^{\texttt{dr_tf_wp_with_insulation}}
$$

The area of the winding pack with no insulation or gap is:

$$
\overbrace{A_{\text{TF,WP-no-insulation}}}^{\texttt{a_tf_wp_no_insulation}} = 
\\ \left(\mathrm{d}R_{\text{TF,WP}} - 2\times\left(\overbrace{\mathrm{d}R_{\text{TF,WP-insulation}}}^{\texttt{dx_tf_wp_insulation}} + \overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)
\\ \times \left(\mathrm{d}x_{\text{TF,toroidal,WP-primary}} - 2\times\left(\overbrace{\mathrm{d}R_{\text{TF,WP-insulation}}}^{\texttt{dx_tf_wp_insulation}} + \overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)
$$

The area of the surrounding winding pack insulation is:

$$
\overbrace{A_{\text{TF,WP-insulation}}}^{\texttt{a_tf_wp_ground_insulation}} =
\\ \left(\left(\mathrm{d}R_{\text{TF,WP}} - 2\times\overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)
\\ \times \left(\mathrm{d}x_{\text{TF,toroidal,WP-primary}} - 2\times\overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)
\\ - A_{\text{TF,WP-no-insulation}}
$$


--------

#### Double rectangular WP

For a [double rectangular winding pack](#winding-pack-geometry) (`i_tf_wp_geom == 1`):

The straight toroidal width of the primary winding pack is:

$$
\overbrace{\mathrm{d}x_{\text{TF,toroidal,WP-primary}}}^{\texttt{dx_tf_wp_primary_toroidal}} = 2 \times 
\\ \left(\overbrace{R_{\text{TF,WP-centre}}}^{\texttt{r_tf_wp_inboard_centre}} \times \overbrace{\tan{\left(\phi_{\text{TF,half}}\right)}}^{\texttt{tan_theta_coil}} - \overbrace{\mathrm{d}x_{\text{TF,side case}}}^{\texttt{dx_tf_side_case_min}}\right)
$$

The straight toroidal width of the secondary winding pack is:

$$
\overbrace{\mathrm{d}x_{\text{TF,toroidal,WP-secondary}}}^{\texttt{dx_tf_wp_secondary_toroidal}} = 2 \times 
\\ \left(\overbrace{R_{\text{TF,WP-inner}}}^{\texttt{r_tf_wp_inboard_inner}} \times \tan{\left(\phi_{\text{TF,half}}\right)} - \mathrm{d}x_{\text{TF,side case}}\right)
$$

The average toroidal straight width is calculated:

$$
\overbrace{\mathrm{d}x_{\text{TF,toroidal,WP-average}}}^{\texttt{dx_tf_wp_toroidal_average}} = \frac{\left(\mathrm{d}x_{\text{TF,toroidal,WP-secondary}} + \mathrm{d}x_{\text{TF,toroidal,WP-primary}}\right)}{2}
$$

The total winding pack area is calculated from the average:

$$
\overbrace{A_{\text{TF,WP}}}^{\texttt{a_tf_wp_with_insulation}} = \mathrm{d}x_{\text{TF,toroidal,WP-average}} \times \overbrace{\mathrm{d}R_{\text{TF,WP}}}^{\texttt{dr_tf_wp_with_insulation}}
$$

The area of the winding pack with no insulation or gap is:

$$
\overbrace{A_{\text{TF,WP-no-insulation}}}^{\texttt{a_tf_wp_no_insulation}} = 0.5 \times
\\ \left(\mathrm{d}R_{\text{TF,WP}} - 2\times\left(\overbrace{\mathrm{d}R_{\text{TF,WP-insulation}}}^{\texttt{dx_tf_wp_insulation}} + \overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)
\\ \times \left(\mathrm{d}x_{\text{TF,toroidal,WP-primary}} + \mathrm{d}x_{\text{TF,toroidal,WP-secondary}} - 
\\ 4\times\left(\overbrace{\mathrm{d}R_{\text{TF,WP-insulation}}}^{\texttt{dx_tf_wp_insulation}} + \overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)
$$

The area of the surrounding winding pack insulation is:

$$
\overbrace{A_{\text{TF,WP-insulation}}}^{\texttt{a_tf_wp_ground_insulation}} = 0.5 \times
\\ \left(\left(\mathrm{d}R_{\text{TF,WP}} - 2\times\overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)
\\ \times \left(\mathrm{d}x_{\text{TF,toroidal,WP-primary}} + \mathrm{d}x_{\text{TF,toroidal,WP-secondary}} - 4\times\overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)
\\ - A_{\text{TF,WP-no-insulation}}
$$

--------

#### Trapezoidal WP 

For a [trapezoidal winding pack](#winding-pack-geometry) (`i_tf_wp_geom == 2`):

The straight toroidal width of the primary winding pack is (longest side of trapezoid):

$$
\overbrace{\mathrm{d}x_{\text{TF,toroidal,WP-primary}}}^{\texttt{dx_tf_wp_primary_toroidal}} = 2 \times 
\\ \left(\overbrace{R_{\text{TF,WP-outer}}}^{\texttt{r_tf_wp_inboard_outer}} \times \overbrace{\tan{\left(\phi_{\text{TF,half}}\right)}}^{\texttt{tan_theta_coil}} - \overbrace{\mathrm{d}x_{\text{TF,side case}}}^{\texttt{dx_tf_side_case_min}}\right)
$$

The straight toroidal width of the secondary winding pack is (shortest side of trapezoid):

$$
\overbrace{\mathrm{d}x_{\text{TF,toroidal,WP-secondary}}}^{\texttt{dx_tf_wp_secondary_toroidal}} = 2 \times 
\\ \left(\overbrace{R_{\text{TF,WP-inner}}}^{\texttt{r_tf_wp_inboard_inner}} \times \tan{\left(\phi_{\text{TF,half}}\right)} - \mathrm{d}x_{\text{TF,side case}}\right)
$$

The average toroidal straight width is calculated:

$$
\overbrace{\mathrm{d}x_{\text{TF,toroidal,WP-average}}}^{\texttt{dx_tf_wp_toroidal_average}} = \frac{\left(\mathrm{d}x_{\text{TF,toroidal,WP-secondary}} + \mathrm{d}x_{\text{TF,toroidal,WP-primary}}\right)}{2}
$$

The total winding pack area is calculated from the average:

$$
\overbrace{A_{\text{TF,WP}}}^{\texttt{a_tf_wp_with_insulation}} = \frac{\left(\mathrm{d}x_{\text{TF,toroidal,WP-primary}}+\mathrm{d}x_{\text{TF,toroidal,WP-secondary}}\right)}{2} 
\\ \times \overbrace{\mathrm{d}R_{\text{TF,WP}}}^{\texttt{dr_tf_wp_with_insulation}}
$$

The area of the winding pack with no insulation or gap is:

$$
\overbrace{A_{\text{TF,WP-no-insulation}}}^{\texttt{a_tf_wp_no_insulation}} =
\\ \left(\mathrm{d}R_{\text{TF,WP}} - 2\times\left(\overbrace{\mathrm{d}R_{\text{TF,WP-insulation}}}^{\texttt{dx_tf_wp_insulation}} + \overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)
\\ \times \left(\left(\mathrm{d}x_{\text{TF,toroidal,WP-primary}} - 2\times\left(\overbrace{\mathrm{d}R_{\text{TF,WP-insulation}}}^{\texttt{dx_tf_wp_insulation}} + \overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)
\\ + \left(\mathrm{d}x_{\text{TF,toroidal,WP-secondary}} - 2\times\left(\overbrace{\mathrm{d}R_{\text{TF,WP-insulation}}}^{\texttt{dx_tf_wp_insulation}} + \overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)\right)
\\ \times 0.5
$$

The area of the surrounding winding pack insulation is:

$$
\overbrace{A_{\text{TF,WP-insulation}}}^{\texttt{a_tf_wp_ground_insulation}} =
\\ \left(\left(\mathrm{d}R_{\text{TF,WP}} - 2\times\overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)
\\ \times \left(\left(\mathrm{d}x_{\text{TF,toroidal,WP-primary}} - 2\times\overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right) 
\\ + \left(\mathrm{d}x_{\text{TF,toroidal,WP-secondary}} - 2\times\overbrace{\mathrm{d}R_{\text{TF,WP-gap}}}^{\texttt{dx_tf_wp_insertion_gap}}\right)\right)
\\ \times 0.5
\\ - A_{\text{TF,WP-no-insulation}}
$$

-----------------

### Case geometry | `superconducting_tf_case_geometry()`

The areas of the total casing surrounding the winding packs on the inboard and outboard leg are calculated:

$$
\overbrace{A_{\text{TF,case-inboard}}}^{\texttt{a_tf_coil_inboard_case}} = \frac{\overbrace{A_{\text{TF,inboard-total}}}^{\texttt{a_tf_inboard_total}}}{\underbrace{N_{\text{TF,coils}}}_{\texttt{n_tf_coils}}} - \overbrace{A_{\text{TF,WP-no-insulation}}}^{\texttt{a_tf_wp_with_insulation}}
$$

$$
\overbrace{A_{\text{TF,case-outboard}}}^{\texttt{a_tf_coil_outboard_case}} = \overbrace{A_{\text{TF,outboard}}}^{\texttt{a_tf_leg_outboard}} - A_{\text{TF,WP-no-insulation}}
$$

The plasma facing front case area is calculated:

If `i_tf_case_geom == 0` then the front case is circular so:

$$
\overbrace{A_{\text{TF,plasma-case}}}^{\texttt{a_tf_plasma_case}} = \left(\overbrace{\phi_{\text{TF,half}}}^{\texttt{rad_tf_coil_inboard_toroidal_half}} \times \overbrace{R_{\text{TF,inboard-out}}^2}^{\texttt{r_tf_inboard_out}}\right)
\\ - \left(\overbrace{\tan{(\phi_{\text{TF,half}})}}^{\texttt{tan_theta_coil}} \times \overbrace{R_{\text{TF,WP-outer}}^2}^{\texttt{r_tf_wp_inboard_outer}}\right)
$$

The first term is equal to the area of an arc segment of radius $R_{\text{TF,inboard-out}}$. Since the value of `rad_tf_coil_inboard_toroidal_half` is a fraction of $\pi$ for each TF coil it can be substituted as the fraction of a full circle. 

If `i_tf_case_geom == 1` then the front case is straight so:

$$
\overbrace{A_{\text{TF,plasma-case}}}^{\texttt{a_tf_plasma_case}} = \left(\left(\overbrace{R_{\text{TF,WP-outer}}}^{\texttt{r_tf_wp_inboard_outer}} \times \overbrace{\Delta R_{\text{TF,plasma-case}}}^{\texttt{dr_tf_plasma_case}}\right)^2- R_{\text{TF,WP-outer}}^2\right)
\\ \times \tan{(\phi_{\text{TF,half}})}
$$

Next the nose case area is calculated:

$$
\overbrace{A_{\text{TF,nose-case}}}^{\texttt{a_tf_coil_nose_case}} = \left(\tan{(\phi_{\text{TF,half}})} \times \overbrace{R_{\text{TF,WP-inner}}^2}^{\texttt{r_tf_wp_inboard_inner}}\right)
\\ - \left(\phi_{\text{TF,half}} \times \overbrace{R_{\text{TF,inboard-in}}^2}^{\texttt{r_tf_inboard_in}}\right) 
$$

Finally the average side case thickness is calculated:

If `i_tf_wp_geom == 0` then a rectangular casing is:

$$
\overbrace{\Delta x_{\text{TF,side-case-average}}}^{\texttt{dx_tf_side_case_average}} = \overbrace{\Delta x_{\text{TF,side-case-min}}}^{\texttt{dx_tf_side_case_min}} + \frac{1}{2}\left(\tan{(\phi_{\text{TF,half}})} \times \overbrace{\Delta R_{\text{TF,WP}}}^{\texttt{dr_tf_wp_with_insulation}}\right)
$$

This is equal to the sidewall casing thickness at the very centre of the winding pack.

If `i_tf_wp_geom == 1` then a  double rectangular casing is:

$$
\overbrace{\Delta x_{\text{TF,side-case-average}}}^{\texttt{dx_tf_side_case_average}} = \overbrace{\Delta x_{\text{TF,side-case-min}}}^{\texttt{dx_tf_side_case_min}} + \frac{1}{4}\left(\tan{(\phi_{\text{TF,half}})} \times \overbrace{\Delta R_{\text{TF,WP}}}^{\texttt{dr_tf_wp_with_insulation}}\right)
$$

Finally, if `i_tf_wp_geom == 2` then a  trapezoidal casing is:

$$
\overbrace{\Delta x_{\text{TF,side-case-average}}}^{\texttt{dx_tf_side_case_average}} = \overbrace{\Delta x_{\text{TF,side-case-min}}}^{\texttt{dx_tf_side_case_min}}
$$

-------------------

### On coil ripple | `peak_b_tf_inboard_with_ripple()`


The ratio of TF coil magnetic field increase with respect to the axisymmetric
formula has been defined as :

$$
  f_\mathrm{rip}^\mathrm{coil} = \frac{B_\mathrm{rip}}{B_\mathrm{nom}}
$$

with $B_\mathrm{rip}$ being the maximum field measured at the middle of the
plasma facing sides of the winding pack and $B_\mathrm{nom}$ the nominal maximum
field obtained with the axisymmetric formula (see [section TF coils current](../eng-models/tf-coil.md#tf-coil-currents-tf_current)).
The same `FIESTA` runs have been used to estimate the on-coil ripple.
This peaking factor has been fitted separately for 16, 18 and 20 coils using
the following formula:


$$
  f_\mathrm{rip}^\mathrm{coil} = A_0 + A_1 e^{-t} + A_2 z + A_3 zt
$$


with the $A_n$ the fitted coefficients and $t$ the relative winding pack
lateral thickness defined as:


$$
  t = \frac{\Delta R_\mathrm{tWP}^\mathrm{out}}{\Delta R_\mathrm{tWP}^\mathrm{out\ max}}
$$


with $\Delta R_\mathrm{tWP}^\mathrm{out}$ defined in Figure 1 and 3 and
$\Delta R_\mathrm{tWP}^\mathrm{out\ max}$ the same value calculated without
sidewall case as illustrated in Figure 12.

  
<figure>
    <center>
    <img src="../../images/tfcoil_ripple_on_coil_filaments.png" alt="On_coil_ripple_shapes" 
    title="Inboard mid-plane superconducting TF coil stress layers"
    width="800" height="100" />
    <br>
    <figcaption><i>
      <p style='text-align: justify;'> 
        Figure 5: Inboard midplane cross section showing the filaments used for
        the on-coil ripple calculations with `FIESTA` with 16 coils. A
        radial and toroidal winding pack thickness of 0.5 m and 0.8 m for the 
        left graph, while a radial and toroidal thickness of 1.2 m and 1.9 m, 
        respectively, is used in the right figure. The right graph provides
        a visual illustration of the case where the parameter $t$ is close to
        unity.
      </p>
    </i></figcaption>
    <br>
    </center>
</figure>
  
And the relative winding pack radial thickness $z$ given by


$$
z = \frac{\Delta R_\mathrm{WP}}{\Delta R_\mathrm{tWP}^\mathrm{out\ max}}
$$


The three fits (for 16, 18 and 20 coils) are valid for:

- **relative toroidal thickness:**  $t\in[0.35-0.99]$

- **relative radial thickness:** $z\in[0.2-0.7]$

- **Number of TF coils:** Individual fits has been made for 16, 18 and 20.
    <b>For any other number of coils, the *FIESTA* calculations are not used</b>
    and a default ripple increase of 9% is taken (\( f_\mathrm{rip}^\mathrm{coil}
    = 1.09\)). This default value is also used for 17 and 19 coils.


Figure 6 shows a contour plot of the on-coil ripple peaking factor as a
function of the winding pack sizing parameters for 16 coil.


<figure>
    <center>
    <img src="../../images/tfcoil_ripple_on_coil_contours.png" alt="On_coil_ripple_shapes" 
    title="Inboard mid-plane superconducting TF coil stress layers"
    width="520" height="100" />
    <br>
    <figcaption><i>
      <p style='text-align: justify;'> 
        Figure 6: Contour plots of the on-coil peaking factor \(f_\mathrm{rip}^\mathrm{coil}\) obtained with FIESTA and used in the <code>PROCESS</code> scaling for 16 coils. The horizontal and vertical axis corresponds to the relative transverse \(\mathit{t}\) and radial thickness \(\mathit{z}\), respectively.
      </p>
    </i></figcaption>
    <br>
    </center>
</figure>


These ripple calculations are out of the spherical tokamak design range, which generally have
fewer coils (between 10 and 14) and more radially thick winding packs.
It is also worth mentioning that the the ripple must be evaluated layer-by-layer for graded coil designs, to get the genuine B field of each layer
used to quantify the SC cross-section area per layer. Finally, resistive
coils do not suffer from on-coil ripple as there is no radial case present.

