# PF Coil Model

The poloidal field (PF) coils are used initially to cancel the vertical field produced at the 
centre of the plasma by the central solenoid during start-up, and then to maintain the plasma 
position and shape during the flat-top period.

## Positioning

The positions and sizes of te PF coils are partly input, and partly calculated after consideration 
of the required currents and allowable current density.

The PF coil locations are controlled using a set of switched stored in array `i_pf_location` (see 
Figure 1), and are calculated in routine `PFCOIL`. The coils are (usually) organised into groups 
containing two PF coils placed symmetrically above and below the midplane, and each group `j` has 
an element `i_pf_location(j)` assigned to it. Input parameter `n_pf_coil_groups` should be set to the number of groups, 
and `n_pf_coils_in_group(j)` should be assigned the number of coils in each group - which should be 2 in each case.

<figure markdown>
![Machine build](../images/vertical-build.png){ width="100%"}
<figcaption>Figure 1: Machine build for D-shaped major components</figcaption>
</figure>

In the following, all variables are defined in the variable descriptor file `vardes.html`. The 
values for `rpf1`, `rpf2`, `zref(j)` and `routr` should be adjusted by the user to locate the PF 
coils accurately.

The three possible values of `i_pf_location(j)` correspond to the following PF coil positions: (Redo taking 
into account `i_single_null` and other recent changes e.g. rclsnorm)

`i_pf_location(j)` = 1: PF coils are placed above the central solenoid (one group only);
*R* = `r_cs_middle` + `rpf1`<br>
*Z* = $\pm$(`z_tf_inside_half` * `f_z_cs_tf_internal` + 0.1 + 0.5 * (`z_tf_inside_half` * (1 - `f_z_cs_tf_internal`) + `dr_tf_inboard` + 0.1))

`i_pf_location(j)` = 2: PF coils are placed above the TF coils (one group only);<br>
*R* = `rmajor` + `rpf2`<br>
*Z* = $\pm$(`z_tf_inside_half` * `dr_tf_inboard` + 0.86)

`i_pf_location(j)` = 3: PF coils are placed radially outside the TF coils (any number of groups);<br>
*R* = `rtot` + `dr_tf_outboard`/2 + `routr`<br>
*Z* = $\pm$(`rminor` * `zref(j)`

The void fraction (for coolant) in each coil `i`'s winding pack is given by `f_a_pf_coil_void(i)`.

## Coil currents

The peak current per turn, `c_pf_coil_turn_peak_input(i)`, and the winding pack peak current density `j_pf_coil_wp_peak(i)` in 
each PF coil `i` are inputs. The PF coil currents vary as a function of time during the tokamak 
operation as indicated in Figure 2. They contribute part of the flux swing necessary to maintain the plasma current.

<figure markdown>
![Current waveform for Plasma, PF coil and Central Solenoid](../images/current_vs_time.png){ width="100%"}
<figcaption>Figure 2: Plot showing schematically the current waveforms for the plasma, a typical PF 
coil, and the central solenoid. Note that the currents in some of the PF coils may be the opposite 
sign to that shown, and the central solenoid current may remain positive during the I<sub>p</sub> 
ramp-up period, although it will pass through zero during the burn phase.</figcaption>
</figure>

## Materials

The PF coils can be either resistive or superconducting. This is determined from the value of 
`i_pf_conductor`. If `i_pf_conductor` = 0, the PF coils and the central solenoid are assumed to be superconducting. 
If `i_pf_conductor` = 1, they are assumed to be resistive, with their resistivity given by the value of variable `rho_pf_coil`.

If `i_pf_conductor` = 0, switch `i_pf_superconductor` specifies which superconducting material is to be used for the 
PF coils. The values of `i_pf_superconductor` are used in the same way as switch `isumattf` is for the TF coils.

The fraction of copper present in the superconducting filaments if given by the value of 
variable `fcupfsu`.

If the PF coils are superconducting, a steel case is assumed to surround the current-carrying 
winding pack to take the hoop stress. Its cross-sectional area is determined by the *J* $\times$ 
*B* hoop force on the coil divided by the allowable hoop stress, given by input parameter `sigpfcalw`. 
The input parameters `sigpfcf` provides a scale factor (default is 0.666) to adjust the hoop force 
if required, to indicate what proportion of the force is supported by the case.
