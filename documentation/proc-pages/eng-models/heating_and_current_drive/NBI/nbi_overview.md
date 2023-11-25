# Neutral Beam Injection Heating

!!! Warning "Warning" 
    At present, the neutral beam models do not include the effect of an edge transport barrier (pedestal) in the plasma profile.

## Neutral beam access

If present, a neutral beam injection system needs sufficient space between the TF coils to be able to intercept the plasma tangentially. The major radius `rtanbeam` at which the centre-line of the beam is tangential to the toroidal direction is user-defined using input parameter `frbeam`, which is the ratio of `rtanbeam` to the plasma major radius `rmajor`. The maximum possible tangency radius `rtanmax` is determined by the geometry of the TF coils - see Figure 1, and this can be enforced using constraint equation no. 20 with iteration variable no. 33 (`fportsz`). The thickness of the beam duct walls may be set using input parameter `nbshield`.

<figure>
    <center>
    <img src="../../../images/portsize.png" alt="NBI port" 
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