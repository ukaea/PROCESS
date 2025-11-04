# PROCESS Variable Descriptions
---
## Introduction
Variables marked with an **\*** are private variables and cannot be accessed
outside of their Fortran module scope.

The PROCESS convention on inputs dictates that module variables which can be set from the
input file should be initialised in a routine called `init_mod` where `mod` is replaced with the
Fortran module name. However, some variables can also be initialised here too and will register as inputs
when they are not.

Output types signify that these module variables are set by models. This does necessarily mean an "output"
will appear in the MFile/Outfile.
---

## cs_fatigue_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>residual_sig_hoop</td>
        <td>Input</td>
        <td>real</td>
        <td>240000000.0</td>
        <td><p>residual hoop stress in strucutal material (Pa)</p></td>
    </tr>
    
    <tr>
        <td>n_cycle</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Allowable number of cycles for CS stress model</p></td>
    </tr>
    
    <tr>
        <td>n_cycle_min</td>
        <td>Input</td>
        <td>real</td>
        <td>20000.0</td>
        <td><p>Minimum llowable number of cycles for CS stress model</p></td>
    </tr>
    
    <tr>
        <td>t_crack_radial</td>
        <td>Input</td>
        <td>real</td>
        <td>0.006</td>
        <td><p>Initial depth of crack in thickness of conduit (m)</p></td>
    </tr>
    
    <tr>
        <td>t_crack_vertical</td>
        <td>Input</td>
        <td>real</td>
        <td>0.00089</td>
        <td><p>Inital vertical crack size (m)</p></td>
    </tr>
    
    <tr>
        <td>t_structural_radial</td>
        <td>Input</td>
        <td>real</td>
        <td>0.07</td>
        <td><p>Thickness of CS conductor conduit (m)</p></td>
    </tr>
    
    <tr>
        <td>t_structural_vertical</td>
        <td>Input</td>
        <td>real</td>
        <td>0.022</td>
        <td><p>Vertical thickness of CS conductor conduit (m)</p></td>
    </tr>
    
    <tr>
        <td>bkt_life_csf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Switch to pass bkt_life cycles to n_cycle_min</p></td>
    </tr>
    
    <tr>
        <td>sf_vertical_crack</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>Safety factor for vertical crack size (-)</p></td>
    </tr>
    
    <tr>
        <td>sf_radial_crack</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>Safety factor for radial crack size (-)</p></td>
    </tr>
    
    <tr>
        <td>sf_fast_fracture</td>
        <td>Input</td>
        <td>real</td>
        <td>1.5</td>
        <td><p>safety factor for stress intensity factor (-)</p></td>
    </tr>
    
    <tr>
        <td>paris_coefficient</td>
        <td>Input</td>
        <td>real</td>
        <td>6.5e-13</td>
        <td><p>Paris equation material coefficient (-)</p></td>
    </tr>
    
    <tr>
        <td>paris_power_law</td>
        <td>Input</td>
        <td>real</td>
        <td>3.5</td>
        <td><p>Paris equation material power law (-)</p></td>
    </tr>
    
    <tr>
        <td>walker_coefficient</td>
        <td>Input</td>
        <td>real</td>
        <td>0.436</td>
        <td><p>walker coefficent (-)</p></td>
    </tr>
    
    <tr>
        <td>fracture_toughness</td>
        <td>Input</td>
        <td>real</td>
        <td>200.0</td>
        <td><p>fracture toughness (MPa m^1/2)</p></td>
    </tr>
    
</table>

## blanket_library
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>volshldi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Volume of inboard and outboard shield (m3)</p></td>
    </tr>
    
    <tr>
        <td>volshldo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Volume of inboard and outboard shield (m3)</p></td>
    </tr>
    
    <tr>
        <td>volvvi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Volume of inboard and outboard Vacuum Vessel (m3)</p></td>
    </tr>
    
    <tr>
        <td>volvvo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Volume of inboard and outboard Vacuum Vessel (m3)</p></td>
    </tr>
    
    <tr>
        <td>dz_pf_cryostat</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Clearance between uppermost PF coil and cryostat lid (m)</p></td>
    </tr>
    
    <tr>
        <td>vfblkti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard void fraction of blanket</p></td>
    </tr>
    
    <tr>
        <td>vfblkto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard void fraction of blanket</p></td>
    </tr>
    
    <tr>
        <td>bldepti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket coolant channel length (radial direction) (m)</p></td>
    </tr>
    
    <tr>
        <td>bldepto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket coolant channel length (radial direction) (m)</p></td>
    </tr>
    
    <tr>
        <td>blwidti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket mid-plan toroidal circumference for segment (m)</p></td>
    </tr>
    
    <tr>
        <td>blwidto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket mid-plan toroidal circumference for segment (m)</p></td>
    </tr>
    
    <tr>
        <td>bllengi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket segment poloidal length (m)</p></td>
    </tr>
    
    <tr>
        <td>bllengo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket segment poloidal length (m)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard primary blanket flow lengths (m)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard primary blanket flow lengths (m)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengi_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard secondary blanket flow lengths (m)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengo_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard secondary blanket flow lengths (m)</p></td>
    </tr>
    
    <tr>
        <td>pnucfwi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard first wall nuclear heating (MW)</p></td>
    </tr>
    
    <tr>
        <td>pnucfwo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard first wall nuclear heating (MW)</p></td>
    </tr>
    
    <tr>
        <td>tpeakfwi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard first wall peak temperature (K)</p></td>
    </tr>
    
    <tr>
        <td>tpeakfwo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard first wall peak temperature (K)</p></td>
    </tr>
    
    <tr>
        <td>mffwi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard total mass flow rate to remove inboard FW power (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mffwo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard total mass flow rate to remove inboard FW power (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mffw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard total mass flow rate to remove inboard FW power (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>npfwi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/utboard total number of pipes</p></td>
    </tr>
    
    <tr>
        <td>npfwo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/utboard total number of pipes</p></td>
    </tr>
    
    <tr>
        <td>mffwpi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard mass flow rate per coolant pipe (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mffwpo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard mass flow rate per coolant pipe (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>pnucblkti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Neutron power deposited inboard/outboard blanket blanket (MW)</p></td>
    </tr>
    
    <tr>
        <td>pnucblkto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Neutron power deposited inboard/outboard blanket blanket (MW)</p></td>
    </tr>
    
    <tr>
        <td>mfblkti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket mass flow rate for coolant (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mfblkto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket mass flow rate for coolant (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mfblkt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket mass flow rate for coolant (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mfblkti_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket mass flow rate for liquid breeder (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mfblkto_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket mass flow rate for liquid breeder (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mfblkt_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket mass flow rate for liquid breeder (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mftotal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total mass flow rate for coolant (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>npblkti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard total num of pipes</p></td>
    </tr>
    
    <tr>
        <td>npblkto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard total num of pipes</p></td>
    </tr>
    
    <tr>
        <td>mfblktpi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard mass flow rate per coolant pipe (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>mfblktpo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard mass flow rate per coolant pipe (kg/s)</p></td>
    </tr>
    
    <tr>
        <td>velblkti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard coolant velocity in blanket (m/s)</p></td>
    </tr>
    
    <tr>
        <td>velblkto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard coolant velocity in blanket (m/s)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_fwi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard first wall pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_fwo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard first wall pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_blkti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_blkto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard blanket pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_fw_blkti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard fw and blanket pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_fw_blkto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard fw and blanket pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>hblnkt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Blanket internal half-height (m)</p></td>
    </tr>
    
    <tr>
        <td>hshld</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Shield internal half-height (m)</p></td>
    </tr>
    
    <tr>
        <td>hvv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Vacuum vessel internal half-height (m)</p></td>
    </tr>
    
    <tr>
        <td>icomponent</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch used to specify selected component: blanket=0, shield=1, vacuum vessel=2</p></td>
    </tr>
    
</table>

## build_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>aplasmin</td>
        <td>Input</td>
        <td>real</td>
        <td>0.25</td>
        <td><p>minimum minor radius (m)</p></td>
    </tr>
    
    <tr>
        <td>available_radial_space</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Minimal radial space between plasma and coils (m)</p></td>
    </tr>
    
    <tr>
        <td>blarea</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>blanket total surface area (m2)</p></td>
    </tr>
    
    <tr>
        <td>blareaib</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>inboard blanket surface area (m2)</p></td>
    </tr>
    
    <tr>
        <td>blareaob</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>outboard blanket surface area (m2)</p></td>
    </tr>
    
    <tr>
        <td>blbmith</td>
        <td>Input</td>
        <td>real</td>
        <td>0.17</td>
        <td><p>inboard blanket box manifold thickness (m) (<code>blktmodel&gt;0</code>)</p></td>
    </tr>
    
    <tr>
        <td>blbmoth</td>
        <td>Input</td>
        <td>real</td>
        <td>0.27</td>
        <td><p>outboard blanket box manifold thickness (m) (<code>blktmodel&gt;0</code>)</p></td>
    </tr>
    
    <tr>
        <td>blbpith</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>inboard blanket base plate thickness (m) (<code>blktmodel&gt;0</code>)</p></td>
    </tr>
    
    <tr>
        <td>blbpoth</td>
        <td>Input</td>
        <td>real</td>
        <td>0.35</td>
        <td><p>outboard blanket base plate thickness (m) (<code>blktmodel&gt;0</code>)</p></td>
    </tr>
    
    <tr>
        <td>blbuith</td>
        <td>Input</td>
        <td>real</td>
        <td>0.365</td>
        <td><p>inboard blanket breeding zone thickness (m) (<code>blktmodel&gt;0</code>) (<code>iteration variable 90</code>)</p></td>
    </tr>
    
    <tr>
        <td>blbuoth</td>
        <td>Input</td>
        <td>real</td>
        <td>0.465</td>
        <td><p>outboard blanket breeding zone thickness (m) (<code>blktmodel&gt;0</code>) (<code>iteration variable 91</code>)</p></td>
    </tr>
    
    <tr>
        <td>dr_blkt_inboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.115</td>
        <td><p>inboard blanket thickness (m); (calculated if <code>blktmodel&gt;0</code>) (=0.0 if <code>i_blkt_inboard=0</code>)</p></td>
    </tr>
    
    <tr>
        <td>dr_blkt_outboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.235</td>
        <td><p>outboard blanket thickness (m); calculated if <code>blktmodel&gt;0</code></p></td>
    </tr>
    
    <tr>
        <td>blnktth</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>top blanket thickness (m), = mean of inboard and outboard blanket thicknesses</p></td>
    </tr>
    
    <tr>
        <td>dr_bore</td>
        <td>Input</td>
        <td>real</td>
        <td>1.42</td>
        <td><p>central solenoid inboard radius (m) (<code>iteration variable 29</code>)</p></td>
    </tr>
    
    <tr>
        <td>f_z_cryostat</td>
        <td>Input</td>
        <td>real</td>
        <td>4.268</td>
        <td><p>cryostat lid height scaling factor (tokamaks)</p></td>
    </tr>
    
    <tr>
        <td>dr_cryostat</td>
        <td>Input</td>
        <td>real</td>
        <td>0.07</td>
        <td><p>cryostat thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>dr_vv_inboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.07</td>
        <td><p>vacuum vessel inboard thickness (TF coil / shield) (m)</p></td>
    </tr>
    
    <tr>
        <td>dr_vv_outboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.07</td>
        <td><p>vacuum vessel outboard thickness (TF coil / shield) (m)</p></td>
    </tr>
    
    <tr>
        <td>d_vv_top</td>
        <td>Input</td>
        <td>real</td>
        <td>0.07</td>
        <td><p>vacuum vessel topside thickness (TF coil / shield) (m) (= d_vv_bot if double-null)</p></td>
    </tr>
    
    <tr>
        <td>d_vv_bot</td>
        <td>Input</td>
        <td>real</td>
        <td>0.07</td>
        <td><p>vacuum vessel underside thickness (TF coil / shield) (m)</p></td>
    </tr>
    
    <tr>
        <td>f_avspace</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>F-value for stellarator radial space check (<code>constraint equation 83</code>)</p></td>
    </tr>
    
    <tr>
        <td>fcspc</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6</td>
        <td><p>Fraction of space occupied by CS pre-compression structure</p></td>
    </tr>
    
    <tr>
        <td>fseppc</td>
        <td>Input</td>
        <td>real</td>
        <td>350000000.0</td>
        <td><p>Separation force in CS coil pre-compression structure</p></td>
    </tr>
    
    <tr>
        <td>a_fw_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>First wall total surface area [m^2]</p></td>
    </tr>
    
    <tr>
        <td>a_fw_inboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard first wall surface area [m^2]</p></td>
    </tr>
    
    <tr>
        <td>a_fw_outboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Outboard first wall surface area [m^2]</p></td>
    </tr>
    
    <tr>
        <td>dr_fw_inboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>inboard first wall thickness, initial estimate as calculated (m)</p></td>
    </tr>
    
    <tr>
        <td>dr_fw_outboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>outboard first wall thickness, initial estimate as calculated (m)</p></td>
    </tr>
    
    <tr>
        <td>dr_shld_vv_gap_inboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.155</td>
        <td><p>gap between inboard vacuum vessel and thermal shield (m) (<code>iteration variable 61</code>)</p></td>
    </tr>
    
    <tr>
        <td>dr_cs_tf_gap</td>
        <td>Input</td>
        <td>real</td>
        <td>0.08</td>
        <td><p>gap between central solenoid and TF coil (m) (<code>iteration variable 42</code>)</p></td>
    </tr>
    
    <tr>
        <td>gapomin</td>
        <td>Input</td>
        <td>real</td>
        <td>0.234</td>
        <td><p>minimum gap between outboard vacuum vessel and TF coil (m) (<code>iteration variable 31</code>)</p></td>
    </tr>
    
    <tr>
        <td>dr_shld_vv_gap_outboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>gap between outboard vacuum vessel and TF coil (m)</p></td>
    </tr>
    
    <tr>
        <td>hmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum (half-)height of TF coil (inside edge) (m)</p></td>
    </tr>
    
    <tr>
        <td>hpfdif</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>difference in distance from midplane of upper and lower portions of TF
 legs (non-zero for single-null devices) (m)</p></td>
    </tr>
    
    <tr>
        <td>hpfu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>height to top of (upper) TF coil leg (m)</p></td>
    </tr>
    
    <tr>
        <td>hr1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>half-height of TF coil inboard leg straight section (m)</p></td>
    </tr>
    
    <tr>
        <td>iohcl</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for existence of central solenoid:</p>
<ul>
<li>=0 central solenoid not present</li>
<li>=1 central solenoid exists</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_cs_precomp</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for existence of central solenoid pre-compression structure:</p>
<ul>
<li>=0 no pre-compression structure</li>
<li>=1 calculated pre-compression structure</li>
</ul></td>
    </tr>
    
    <tr>
        <td>tf_in_cs</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for placing the TF coil inside the CS</p>
<ul>
<li>= 0 TF coil is outside the CS (default)</li>
<li>= 1 TF coil is inside the CS</li>
</ul></td>
    </tr>
    
    <tr>
        <td>dr_cs</td>
        <td>Input</td>
        <td>real</td>
        <td>0.811</td>
        <td><p>Central solenoid thickness (m) (<code>iteration variable 16</code>)</p></td>
    </tr>
    
    <tr>
        <td>dr_cs_precomp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>CS coil precompression structure thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>rbld</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>sum of thicknesses to the major radius (m)</p></td>
    </tr>
    
    <tr>
        <td>required_radial_space</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Required space between coil and plasma for blanket shield wall etc (m)</p></td>
    </tr>
    
    <tr>
        <td>rinboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.651</td>
        <td><p>plasma inboard radius (m) (<code>consistency equation 29</code>)</p></td>
    </tr>
    
    <tr>
        <td>rsldi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radius to inboard shield (inside point) (m)</p></td>
    </tr>
    
    <tr>
        <td>rsldo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radius to outboard shield (outside point) (m)</p></td>
    </tr>
    
    <tr>
        <td>r_vv_inboard_out</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial plasma facing side position of inboard vacuum vessel [m]</p></td>
    </tr>
    
    <tr>
        <td>r_sh_inboard_in</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial inner side position of inboard neutronic shield [m]</p></td>
    </tr>
    
    <tr>
        <td>r_sh_inboard_out</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial plasma facing side position of inboard neutronic shield [m]</p></td>
    </tr>
    
    <tr>
        <td>r_tf_inboard_in</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Mid-plane inboard TF coil leg radius at the centre-machine side [m]</p></td>
    </tr>
    
    <tr>
        <td>r_tf_inboard_mid</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Mid-plane inboard TF coil leg radius at middle of the coil [m]</p></td>
    </tr>
    
    <tr>
        <td>r_tf_inboard_out</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Mid-plane inboard TF coil leg radius at the plasma side [m]</p></td>
    </tr>
    
    <tr>
        <td>r_tf_outboard_mid</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Mid-plane outboard TF coil leg radius at the middle of the coil [m]</p></td>
    </tr>
    
    <tr>
        <td>i_r_cp_top</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch selecting the he parametrization of the outer radius of the top of the CP part of the TF coil
  0 : <code>r_cp_top</code> is set by the plasma shape
  1 : <code>r_cp_top</code> is a user input
  2 : <code>r_cp_top</code> is set using the CP top and midplane CP radius ratio</p></td>
    </tr>
    
    <tr>
        <td>r_cp_top</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Top outer radius of the centropost (ST only) (m)</p></td>
    </tr>
    
    <tr>
        <td>f_r_cp</td>
        <td>Input</td>
        <td>real</td>
        <td>1.4</td>
        <td><p>Ratio between the top and the midplane TF CP outer radius [-]
 Not used by default (-1) must be larger than 1 otherwise</p></td>
    </tr>
    
    <tr>
        <td>dr_tf_inner_bore</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil horizontal inner dr_bore (m)</p></td>
    </tr>
    
    <tr>
        <td>dh_tf_inner_bore</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil vertical inner dr_bore (m)</p></td>
    </tr>
    
    <tr>
        <td>dr_fw_plasma_gap_inboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.14</td>
        <td><p>Gap between plasma and first wall, inboard side (m) (if <code>i_plasma_wall_gap=1</code>)
 Iteration variable: ixc = 73
 Scan variable: nsweep = 58</p></td>
    </tr>
    
    <tr>
        <td>dr_fw_plasma_gap_outboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.15</td>
        <td><p>Gap between plasma and first wall, outboard side (m) (if <code>i_plasma_wall_gap=1</code>)
 Iteration variable: ixc = 74
 Scan variable: nsweep = 59</p></td>
    </tr>
    
    <tr>
        <td>sharea</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>shield total surface area (m2)</p></td>
    </tr>
    
    <tr>
        <td>shareaib</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>inboard shield surface area (m2)</p></td>
    </tr>
    
    <tr>
        <td>shareaob</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>outboard shield surface area (m2)</p></td>
    </tr>
    
    <tr>
        <td>dr_shld_inboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.69</td>
        <td><p>inboard shield thickness (m) (<code>iteration variable 93</code>)</p></td>
    </tr>
    
    <tr>
        <td>shldlth</td>
        <td>Input</td>
        <td>real</td>
        <td>0.7</td>
        <td><p>lower (under divertor) shield thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>dr_shld_outboard</td>
        <td>Input</td>
        <td>real</td>
        <td>1.05</td>
        <td><p>outboard shield thickness (m) (<code>iteration variable 94</code>)</p></td>
    </tr>
    
    <tr>
        <td>shldtth</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6</td>
        <td><p>upper/lower shield thickness (m); calculated if <code>blktmodel &gt; 0</code> (= shldlth if double-null)</p></td>
    </tr>
    
    <tr>
        <td>sigallpc</td>
        <td>Input</td>
        <td>real</td>
        <td>300000000.0</td>
        <td><p>allowable stress in CSpre-compression structure (Pa)</p></td>
    </tr>
    
    <tr>
        <td>dr_tf_inboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>inboard TF coil thickness, (centrepost for ST) (m)
 (input, calculated or <code>iteration variable 13</code>)</p></td>
    </tr>
    
    <tr>
        <td>tfoffset</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vertical distance between centre of TF coils and centre of plasma (m)</p></td>
    </tr>
    
    <tr>
        <td>tfootfi</td>
        <td>Input</td>
        <td>real</td>
        <td>1.19</td>
        <td><p>TF coil outboard leg / inboard leg radial thickness
 ratio (<code>i_tf_sup=0</code> only) (<code>iteration variable 75</code>)</p></td>
    </tr>
    
    <tr>
        <td>dr_tf_outboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Outboard TF coil thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>dr_tf_shld_gap</td>
        <td>Input</td>
        <td>real</td>
        <td>0.05</td>
        <td><p>Minimum metal-to-metal gap between TF coil and thermal shield (m)</p></td>
    </tr>
    
    <tr>
        <td>dr_shld_thermal_inboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.05</td>
        <td><p>TF-VV thermal shield thickness, inboard (m)</p></td>
    </tr>
    
    <tr>
        <td>dr_shld_thermal_outboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.05</td>
        <td><p>TF-VV thermal shield thickness, outboard (m)</p></td>
    </tr>
    
    <tr>
        <td>thshield_vb</td>
        <td>Input</td>
        <td>real</td>
        <td>0.05</td>
        <td><p>TF-VV thermal shield thickness, vertical build (m)</p></td>
    </tr>
    
    <tr>
        <td>vgap_vv_thermalshield</td>
        <td>Input</td>
        <td>real</td>
        <td>0.163</td>
        <td><p>vertical gap between vacuum vessel and thermal shields (m)</p></td>
    </tr>
    
    <tr>
        <td>vgap_xpoint_divertor</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vertical gap between x-point and divertor (m) (if = 0, it is calculated)</p></td>
    </tr>
    
    <tr>
        <td>vgaptop</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6</td>
        <td><p>vertical gap between top of plasma and first wall (m) (= vgap_xpoint_divertor if double-null)</p></td>
    </tr>
    
    <tr>
        <td>dr_shld_blkt_gap</td>
        <td>Input</td>
        <td>real</td>
        <td>0.05</td>
        <td><p>gap between vacuum vessel and blanket (m)</p></td>
    </tr>
    
    <tr>
        <td>plleni</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>length of inboard divertor plate (m)</p></td>
    </tr>
    
    <tr>
        <td>plleno</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>length of outboard divertor plate (m)</p></td>
    </tr>
    
    <tr>
        <td>plsepi</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>poloidal length, x-point to inboard strike point (m)</p></td>
    </tr>
    
    <tr>
        <td>plsepo</td>
        <td>Input</td>
        <td>real</td>
        <td>1.5</td>
        <td><p>poloidal length, x-point to outboard strike point (m)</p></td>
    </tr>
    
    <tr>
        <td>rspo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>outboard strike point radius (m)</p></td>
    </tr>
    
</table>

## buildings_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>admv</td>
        <td>Input</td>
        <td>real</td>
        <td>100000.0</td>
        <td><p>administration building volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>admvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of administration buildings (m3)</p></td>
    </tr>
    
    <tr>
        <td>aux_build_l</td>
        <td>Input</td>
        <td>real</td>
        <td>60.0</td>
        <td><p>aux building supporting tokamak processes length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>aux_build_w</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>aux building supporting tokamak processes length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>aux_build_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>aux building supporting tokamak processes length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>auxcool_l</td>
        <td>Input</td>
        <td>real</td>
        <td>20.0</td>
        <td><p>Site-Wide Auxiliary Cooling Water facility length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>auxcool_w</td>
        <td>Input</td>
        <td>real</td>
        <td>20.0</td>
        <td><p>Site-Wide Auxiliary Cooling Water facility length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>auxcool_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>Site-Wide Auxiliary Cooling Water facility length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>bioshld_thk</td>
        <td>Input</td>
        <td>real</td>
        <td>2.5</td>
        <td><p>Radial thickness of bio-shield around reactor (m)</p></td>
    </tr>
    
    <tr>
        <td>chemlab_l</td>
        <td>Input</td>
        <td>real</td>
        <td>50.0</td>
        <td><p>Chemistry labs and treatment buldings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>chemlab_w</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>Chemistry labs and treatment buldings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>chemlab_h</td>
        <td>Input</td>
        <td>real</td>
        <td>6.0</td>
        <td><p>Chemistry labs and treatment buldings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>dz_tf_cryostat</td>
        <td>Input</td>
        <td>real</td>
        <td>2.5</td>
        <td><p>vertical clearance from TF coil to cryostat (m) (calculated for tokamaks)</p></td>
    </tr>
    
    <tr>
        <td>clh2</td>
        <td>Input</td>
        <td>real</td>
        <td>15.0</td>
        <td><p>clearance beneath TF coil to foundation (including basement) (m)</p></td>
    </tr>
    
    <tr>
        <td>control_buildings_l</td>
        <td>Input</td>
        <td>real</td>
        <td>80.0</td>
        <td><p>control building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>control_buildings_w</td>
        <td>Input</td>
        <td>real</td>
        <td>60.0</td>
        <td><p>control building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>control_buildings_h</td>
        <td>Input</td>
        <td>real</td>
        <td>6.0</td>
        <td><p>control building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>conv</td>
        <td>Input</td>
        <td>real</td>
        <td>60000.0</td>
        <td><p>control building volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>convol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of control, protection and i&amp;c building (m3)</p></td>
    </tr>
    
    <tr>
        <td>crane_arm_h</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>vertical dimension of crane arm, operating over reactor (m)</p></td>
    </tr>
    
    <tr>
        <td>crane_clrnc_h</td>
        <td>Input</td>
        <td>real</td>
        <td>4.0</td>
        <td><p>horizontal clearance to building wall for crane operation (m)</p></td>
    </tr>
    
    <tr>
        <td>crane_clrnc_v</td>
        <td>Input</td>
        <td>real</td>
        <td>3.0</td>
        <td><p>vertical clearance for crane operation (m)</p></td>
    </tr>
    
    <tr>
        <td>cryomag_l</td>
        <td>Input</td>
        <td>real</td>
        <td>120.0</td>
        <td><p>Cryogenic Buildings for Magnet and Fuel Cycle length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>cryomag_w</td>
        <td>Input</td>
        <td>real</td>
        <td>90.0</td>
        <td><p>Cryogenic Buildings for Magnet and Fuel Cycle length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>cryomag_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>Cryogenic Buildings for Magnet and Fuel Cycle length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>cryostore_l</td>
        <td>Input</td>
        <td>real</td>
        <td>160.0</td>
        <td><p>Magnet Cryo Storage Tanks length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>cryostore_w</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>Magnet Cryo Storage Tanks length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>cryostore_h</td>
        <td>Input</td>
        <td>real</td>
        <td>20.0</td>
        <td><p>Magnet Cryo Storage Tanks length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>cryostat_clrnc</td>
        <td>Input</td>
        <td>real</td>
        <td>2.5</td>
        <td><p>vertical clearance from TF coil to cryostat (m)</p></td>
    </tr>
    
    <tr>
        <td>cryvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of cryoplant building (m3)</p></td>
    </tr>
    
    <tr>
        <td>efloor</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>effective total floor space (m2)</p></td>
    </tr>
    
    <tr>
        <td>elecdist_l</td>
        <td>Input</td>
        <td>real</td>
        <td>380.0</td>
        <td><p>Transformers and electrical distribution facilities length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>elecdist_w</td>
        <td>Input</td>
        <td>real</td>
        <td>350.0</td>
        <td><p>Transformers and electrical distribution facilities length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>elecdist_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>Transformers and electrical distribution facilities length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>elecload_l</td>
        <td>Input</td>
        <td>real</td>
        <td>100.0</td>
        <td><p>Electric (eesential and non-essential) load centres length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>elecload_w</td>
        <td>Input</td>
        <td>real</td>
        <td>90.0</td>
        <td><p>Electric (eesential and non-essential) load centres length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>elecload_h</td>
        <td>Input</td>
        <td>real</td>
        <td>3.0</td>
        <td><p>Electric (eesential and non-essential) load centres length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>elecstore_l</td>
        <td>Input</td>
        <td>real</td>
        <td>100.0</td>
        <td><p>Energy Storage facilities length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>elecstore_w</td>
        <td>Input</td>
        <td>real</td>
        <td>60.0</td>
        <td><p>Energy Storage facilities length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>elecstore_h</td>
        <td>Input</td>
        <td>real</td>
        <td>12.0</td>
        <td><p>Energy Storage facilities length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>elevol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of electrical equipment building (m3)</p></td>
    </tr>
    
    <tr>
        <td>esbldgm3</td>
        <td>Input</td>
        <td>real</td>
        <td>1000.0</td>
        <td><p>volume of energy storage equipment building (m3) (not used if <code>i_pulsed_plant=0</code>)</p></td>
    </tr>
    
    <tr>
        <td>fc_building_l</td>
        <td>Input</td>
        <td>real</td>
        <td>60.0</td>
        <td><p>Fuel Cycle facilities length, width (m)</p></td>
    </tr>
    
    <tr>
        <td>fc_building_w</td>
        <td>Input</td>
        <td>real</td>
        <td>60.0</td>
        <td><p>Fuel Cycle facilities length, width (m)</p></td>
    </tr>
    
    <tr>
        <td>fndt</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>foundation thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>gas_buildings_l</td>
        <td>Input</td>
        <td>real</td>
        <td>25.0</td>
        <td><p>air &amp; gas supply (amalgamated) buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>gas_buildings_w</td>
        <td>Input</td>
        <td>real</td>
        <td>15.0</td>
        <td><p>air &amp; gas supply (amalgamated) buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>gas_buildings_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>air &amp; gas supply (amalgamated) buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>ground_clrnc</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>clearance beneath TF coil (m)</p></td>
    </tr>
    
    <tr>
        <td>hcd_building_l</td>
        <td>Input</td>
        <td>real</td>
        <td>70.0</td>
        <td><p>HCD building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>hcd_building_w</td>
        <td>Input</td>
        <td>real</td>
        <td>40.0</td>
        <td><p>HCD building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>hcd_building_h</td>
        <td>Input</td>
        <td>real</td>
        <td>25.0</td>
        <td><p>HCD building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>hccl</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>clearance around components in hot cell (m)</p></td>
    </tr>
    
    <tr>
        <td>hcwt</td>
        <td>Input</td>
        <td>real</td>
        <td>1.5</td>
        <td><p>hot cell wall thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>heat_sink_l</td>
        <td>Input</td>
        <td>real</td>
        <td>160.0</td>
        <td><p>heat sinks length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>heat_sink_w</td>
        <td>Input</td>
        <td>real</td>
        <td>80.0</td>
        <td><p>heat sinks length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>heat_sink_h</td>
        <td>Input</td>
        <td>real</td>
        <td>12.0</td>
        <td><p>heat sinks length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>hot_sepdist</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>hot cell storage component separation distance (m)</p></td>
    </tr>
    
    <tr>
        <td>hotcell_h</td>
        <td>Input</td>
        <td>real</td>
        <td>12.0</td>
        <td><p>hot cell storage and maintenance facility height (m)</p></td>
    </tr>
    
    <tr>
        <td>hw_storage_l</td>
        <td>Input</td>
        <td>real</td>
        <td>20.0</td>
        <td><p>hazardous waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>hw_storage_w</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>hazardous waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>hw_storage_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>hazardous waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>i_bldgs_size</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch between routines estimating building sizes (0 = default; 1 = updated)</p></td>
    </tr>
    
    <tr>
        <td>i_bldgs_v</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch to select verbose output for buildings (1 = verbose)</p></td>
    </tr>
    
    <tr>
        <td>ilw_smelter_l</td>
        <td>Input</td>
        <td>real</td>
        <td>50.0</td>
        <td><p>radioactive waste smelting facility length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>ilw_smelter_w</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>radioactive waste smelting facility length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>ilw_smelter_h</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>radioactive waste smelting facility length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>ilw_storage_l</td>
        <td>Input</td>
        <td>real</td>
        <td>120.0</td>
        <td><p>ILW waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>ilw_storage_w</td>
        <td>Input</td>
        <td>real</td>
        <td>100.0</td>
        <td><p>ILW waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>ilw_storage_h</td>
        <td>Input</td>
        <td>real</td>
        <td>8.0</td>
        <td><p>ILW waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>llw_storage_l</td>
        <td>Input</td>
        <td>real</td>
        <td>45.0</td>
        <td><p>LLW waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>llw_storage_w</td>
        <td>Input</td>
        <td>real</td>
        <td>20.0</td>
        <td><p>LLW waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>llw_storage_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>LLW waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>magnet_pulse_l</td>
        <td>Input</td>
        <td>real</td>
        <td>105.0</td>
        <td><p>pulsed magnet power building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>magnet_pulse_w</td>
        <td>Input</td>
        <td>real</td>
        <td>40.0</td>
        <td><p>pulsed magnet power building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>magnet_pulse_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>pulsed magnet power building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>magnet_trains_l</td>
        <td>Input</td>
        <td>real</td>
        <td>120.0</td>
        <td><p>steady state magnet power trains building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>magnet_trains_w</td>
        <td>Input</td>
        <td>real</td>
        <td>90.0</td>
        <td><p>steady state magnet power trains building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>magnet_trains_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>steady state magnet power trains building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>maint_cont_l</td>
        <td>Input</td>
        <td>real</td>
        <td>125.0</td>
        <td><p>maintenance control building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>maint_cont_w</td>
        <td>Input</td>
        <td>real</td>
        <td>100.0</td>
        <td><p>maintenance control building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>maint_cont_h</td>
        <td>Input</td>
        <td>real</td>
        <td>6.0</td>
        <td><p>maintenance control building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>mbvfac</td>
        <td>Input</td>
        <td>real</td>
        <td>2.8</td>
        <td><p>maintenance building volume multiplication factor</p></td>
    </tr>
    
    <tr>
        <td>nbi_sys_l</td>
        <td>Input</td>
        <td>real</td>
        <td>225.0</td>
        <td><p>NBI system length, width (m)</p></td>
    </tr>
    
    <tr>
        <td>nbi_sys_w</td>
        <td>Input</td>
        <td>real</td>
        <td>185.0</td>
        <td><p>NBI system length, width (m)</p></td>
    </tr>
    
    <tr>
        <td>pfbldgm3</td>
        <td>Input</td>
        <td>real</td>
        <td>20000.0</td>
        <td><p>volume of PF coil power supply building (m3)</p></td>
    </tr>
    
    <tr>
        <td>pibv</td>
        <td>Input</td>
        <td>real</td>
        <td>20000.0</td>
        <td><p>power injection building volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>qnty_sfty_fac</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>quantity safety factor for component use during plant lifetime</p></td>
    </tr>
    
    <tr>
        <td>rbvfac</td>
        <td>Input</td>
        <td>real</td>
        <td>1.6</td>
        <td><p>reactor building volume multiplication factor</p></td>
    </tr>
    
    <tr>
        <td>rbrt</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>reactor building roof thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>rbvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>reactor building volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>rbwt</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>reactor building wall thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>reactor_clrnc</td>
        <td>Input</td>
        <td>real</td>
        <td>4.0</td>
        <td><p>clearance around reactor (m)</p></td>
    </tr>
    
    <tr>
        <td>reactor_fndtn_thk</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>reactor building foundation thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>reactor_hall_l</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>reactor building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>reactor_hall_w</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>reactor building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>reactor_hall_h</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>reactor building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>reactor_roof_thk</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>reactor building roof thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>reactor_wall_thk</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>reactor building wall thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>rmbvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of maintenance and assembly building (m3)</p></td>
    </tr>
    
    <tr>
        <td>robotics_l</td>
        <td>Input</td>
        <td>real</td>
        <td>50.0</td>
        <td><p>robotics buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>robotics_w</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>robotics buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>robotics_h</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>robotics buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>row</td>
        <td>Input</td>
        <td>real</td>
        <td>4.0</td>
        <td><p>clearance to building wall for crane operation (m)</p></td>
    </tr>
    
    <tr>
        <td>rxcl</td>
        <td>Input</td>
        <td>real</td>
        <td>4.0</td>
        <td><p>clearance around reactor (m)</p></td>
    </tr>
    
    <tr>
        <td>sec_buildings_l</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>security &amp; safety buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>sec_buildings_w</td>
        <td>Input</td>
        <td>real</td>
        <td>25.0</td>
        <td><p>security &amp; safety buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>sec_buildings_h</td>
        <td>Input</td>
        <td>real</td>
        <td>6.0</td>
        <td><p>security &amp; safety buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>shmf</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>fraction of shield mass per TF coil to be moved in the maximum shield lift</p></td>
    </tr>
    
    <tr>
        <td>shov</td>
        <td>Input</td>
        <td>real</td>
        <td>100000.0</td>
        <td><p>shops and warehouse volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>shovol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of shops and buildings for plant auxiliaries (m3)</p></td>
    </tr>
    
    <tr>
        <td>staff_buildings_area</td>
        <td>Input</td>
        <td>real</td>
        <td>480000.0</td>
        <td><p>footprint of staff buildings (m2)</p></td>
    </tr>
    
    <tr>
        <td>staff_buildings_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>staff buildings height (m)</p></td>
    </tr>
    
    <tr>
        <td>stcl</td>
        <td>Input</td>
        <td>real</td>
        <td>3.0</td>
        <td><p>clearance above crane to roof (m)</p></td>
    </tr>
    
    <tr>
        <td>tfcbv</td>
        <td>Input</td>
        <td>real</td>
        <td>20000.0</td>
        <td><p>volume of TF coil power supply building (m3) (calculated if TF coils are superconducting)</p></td>
    </tr>
    
    <tr>
        <td>transp_clrnc</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>transportation clearance between components (m)</p></td>
    </tr>
    
    <tr>
        <td>trcl</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>transportation clearance between components (m)</p></td>
    </tr>
    
    <tr>
        <td>triv</td>
        <td>Input</td>
        <td>real</td>
        <td>40000.0</td>
        <td><p>volume of tritium, fuel handling and health physics buildings (m3)</p></td>
    </tr>
    
    <tr>
        <td>turbine_hall_l</td>
        <td>Input</td>
        <td>real</td>
        <td>109.0</td>
        <td><p>turbine hall length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>turbine_hall_w</td>
        <td>Input</td>
        <td>real</td>
        <td>62.0</td>
        <td><p>turbine hall length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>turbine_hall_h</td>
        <td>Input</td>
        <td>real</td>
        <td>15.0</td>
        <td><p>turbine hall length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>tw_storage_l</td>
        <td>Input</td>
        <td>real</td>
        <td>90.0</td>
        <td><p>tritiated waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>tw_storage_w</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>tritiated waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>tw_storage_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>tritiated waste storage building length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>volrci</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>internal volume of reactor building (m3)</p></td>
    </tr>
    
    <tr>
        <td>volnucb</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>sum of nuclear buildings volumes (m3)</p></td>
    </tr>
    
    <tr>
        <td>warm_shop_l</td>
        <td>Input</td>
        <td>real</td>
        <td>100.0</td>
        <td><p>warm shop length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>warm_shop_w</td>
        <td>Input</td>
        <td>real</td>
        <td>50.0</td>
        <td><p>warm shop length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>warm_shop_h</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>warm shop length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>water_buildings_l</td>
        <td>Input</td>
        <td>real</td>
        <td>110.0</td>
        <td><p>water, laundry &amp; drainage buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>water_buildings_w</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>water, laundry &amp; drainage buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>water_buildings_h</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>water, laundry &amp; drainage buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>wgt</td>
        <td>Input</td>
        <td>real</td>
        <td>500000.0</td>
        <td><p>reactor building crane capacity (kg) (calculated if 0 is input)</p></td>
    </tr>
    
    <tr>
        <td>wgt2</td>
        <td>Input</td>
        <td>real</td>
        <td>100000.0</td>
        <td><p>hot cell crane capacity (kg) (calculated if 0 is input)</p></td>
    </tr>
    
    <tr>
        <td>workshop_l</td>
        <td>Input</td>
        <td>real</td>
        <td>150.0</td>
        <td><p>[cold] workshop buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>workshop_w</td>
        <td>Input</td>
        <td>real</td>
        <td>125.0</td>
        <td><p>[cold] workshop buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>workshop_h</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>[cold] workshop buildings length, width, height (m)</p></td>
    </tr>
    
    <tr>
        <td>wrbi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>distance from centre of machine to building wall (m)</p></td>
    </tr>
    
    <tr>
        <td>wsvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of warm shop building (m3)</p></td>
    </tr>
    
    <tr>
        <td>wsvfac</td>
        <td>Input</td>
        <td>real</td>
        <td>1.9</td>
        <td><p>warm shop building volume multiplication factor</p></td>
    </tr>
    
    <tr>
        <td>a_reactor_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>8320.0</td>
        <td><p>Floor area of reactor building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_ee_ps_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>21330.0</td>
        <td><p>Floor area of electrical equipment and power supply building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_aux_services_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>1000.0</td>
        <td><p>Floor area of auxiliary services building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_hot_cell_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>8430.0</td>
        <td><p>Floor area of hot cell building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_reactor_service_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>2440.0</td>
        <td><p>Floor area of reactor service building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_service_water_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>1567.0</td>
        <td><p>Floor area of service water building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_fuel_handling_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>1670.0</td>
        <td><p>Floor area of fuel handling and storage building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_control_room_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>2880.0</td>
        <td><p>Floor area of controlroom building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_ac_ps_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>6423.0</td>
        <td><p>Floor area of AC power supply building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_admin_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>25674.0</td>
        <td><p>Floor area of admin building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_site_service_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>8300.0</td>
        <td><p>Floor area of site service building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_cryo_inert_gas_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>18380.0</td>
        <td><p>Floor area of cryogenics and inert gas storage building in m^2</p></td>
    </tr>
    
    <tr>
        <td>a_security_bldg</td>
        <td>Input</td>
        <td>real</td>
        <td>4552.0</td>
        <td><p>Floor area of security building in m^2</p></td>
    </tr>
    
</table>

## constants
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>iotty</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>6</td>
        <td></td>
    </tr>
    
    <tr>
        <td>nout</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>11</td>
        <td></td>
    </tr>
    
    <tr>
        <td>nplot</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>12</td>
        <td></td>
    </tr>
    
    <tr>
        <td>mfile</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>13</td>
        <td></td>
    </tr>
    
    <tr>
        <td>vfile</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>14</td>
        <td></td>
    </tr>
    
    <tr>
        <td>opt_file</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>15</td>
        <td></td>
    </tr>
    
    <tr>
        <td>sig_file</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>16</td>
        <td></td>
    </tr>
    
    <tr>
        <td>degrad</td>
        <td>Parameter</td>
        <td>real</td>
        <td>0.01745329251D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>electron_charge</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.602176634D-19</td>
        <td></td>
    </tr>
    
    <tr>
        <td>electron_volt</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.602176634D-19</td>
        <td></td>
    </tr>
    
    <tr>
        <td>kiloelectron_volt</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.602176634D-16</td>
        <td></td>
    </tr>
    
    <tr>
        <td>atomic_mass_unit</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.66053906892D-27</td>
        <td></td>
    </tr>
    
    <tr>
        <td>electron_mass</td>
        <td>Parameter</td>
        <td>real</td>
        <td>9.1093837139D-31</td>
        <td></td>
    </tr>
    
    <tr>
        <td>proton_mass</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.67262192595D-27</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_proton_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.0072764665789</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_protium_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.00782503223</td>
        <td></td>
    </tr>
    
    <tr>
        <td>deuteron_mass</td>
        <td>Parameter</td>
        <td>real</td>
        <td>3.3435837768D-27</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_deuteron_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>2.013553212544</td>
        <td></td>
    </tr>
    
    <tr>
        <td>triton_mass</td>
        <td>Parameter</td>
        <td>real</td>
        <td>5.0073567512D-27</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_triton_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>3.01550071597</td>
        <td></td>
    </tr>
    
    <tr>
        <td>neutron_mass</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.67492750056D-27</td>
        <td></td>
    </tr>
    
    <tr>
        <td>alpha_mass</td>
        <td>Parameter</td>
        <td>real</td>
        <td>6.6446573450D-27</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_alpha_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>4.001506179129</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_helium_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>4.002602</td>
        <td></td>
    </tr>
    
    <tr>
        <td>helion_mass</td>
        <td>Parameter</td>
        <td>real</td>
        <td>5.0064127862D-27</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_helion_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>3.014932246932</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_beryllium_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>9.0121831</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_carbon_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>12.0096</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_nitrogen_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>14.00643</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_oxygen_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>15.99903</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_neon_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>20.1797</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_silicon_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>28.084</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_argon_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>39.948</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_iron_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>55.845</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_nickel_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>58.6934</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_krypton_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>83.798</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_xenon_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>131.293</td>
        <td></td>
    </tr>
    
    <tr>
        <td>m_tungsten_amu</td>
        <td>Parameter</td>
        <td>real</td>
        <td>183.84</td>
        <td></td>
    </tr>
    
    <tr>
        <td>speed_light</td>
        <td>Parameter</td>
        <td>real</td>
        <td>299792458D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>d_t_energy</td>
        <td>Parameter</td>
        <td>real</td>
        <td>(((deuteron_mass+triton_mass)-(alpha_mass+neutron_mass))*speed_light**2)</td>
        <td></td>
    </tr>
    
    <tr>
        <td>d_helium_energy</td>
        <td>Parameter</td>
        <td>real</td>
        <td>(((deuteron_mass+helion_mass)-(alpha_mass+proton_mass))*speed_light**2)</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dd_helium_energy</td>
        <td>Parameter</td>
        <td>real</td>
        <td>(((deuteron_mass+deuteron_mass)-(helion_mass+neutron_mass))*speed_light**2)</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dd_triton_energy</td>
        <td>Parameter</td>
        <td>real</td>
        <td>(((deuteron_mass+deuteron_mass)-(triton_mass+proton_mass))*speed_light**2)</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dt_neutron_energy_fraction</td>
        <td>Parameter</td>
        <td>real</td>
        <td>(alpha_mass/(neutron_mass+alpha_mass))</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dt_alpha_energy</td>
        <td>Parameter</td>
        <td>real</td>
        <td>(1.0D0-dt_neutron_energy_fraction)*d_t_energy</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dd_neutron_energy_fraction</td>
        <td>Parameter</td>
        <td>real</td>
        <td>(helion_mass/(neutron_mass+helion_mass))</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dd_proton_energy_fraction</td>
        <td>Parameter</td>
        <td>real</td>
        <td>(triton_mass/(proton_mass+triton_mass))</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dhelium_proton_energy_fraction</td>
        <td>Parameter</td>
        <td>real</td>
        <td>(alpha_mass/(proton_mass+alpha_mass))</td>
        <td></td>
    </tr>
    
    <tr>
        <td>den_tungsten</td>
        <td>Parameter</td>
        <td>real</td>
        <td>19250.0D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>pi</td>
        <td>Parameter</td>
        <td>real</td>
        <td>3.1415926535897932D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>rmu0</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.256637062D-6</td>
        <td></td>
    </tr>
    
    <tr>
        <td>twopi</td>
        <td>Parameter</td>
        <td>real</td>
        <td>6.2831853071795862D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>umass</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.660538921D-27</td>
        <td></td>
    </tr>
    
    <tr>
        <td>epsilon0</td>
        <td>Parameter</td>
        <td>real</td>
        <td>8.85418781D-12</td>
        <td></td>
    </tr>
    
    <tr>
        <td>cph2o</td>
        <td>Parameter</td>
        <td>real</td>
        <td>4180.0D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dcopper</td>
        <td>Input</td>
        <td>real</td>
        <td>8900.0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dalu</td>
        <td>Input</td>
        <td>real</td>
        <td>2700.0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>denh2o</td>
        <td>Parameter</td>
        <td>real</td>
        <td>985.0D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>k_copper</td>
        <td>Parameter</td>
        <td>real</td>
        <td>330.0D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>kh2o</td>
        <td>Parameter</td>
        <td>real</td>
        <td>0.651D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>muh2o</td>
        <td>Parameter</td>
        <td>real</td>
        <td>4.71D-4</td>
        <td></td>
    </tr>
    
    <tr>
        <td>n_day_year</td>
        <td>Parameter</td>
        <td>real</td>
        <td>365.2425D0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>acceleration_gravity</td>
        <td>Parameter</td>
        <td>real</td>
        <td>9.81D0</td>
        <td><p>Acceleration due to gravity [m/s2]</p></td>
    </tr>
    
</table>

## constraint_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>auxmin</td>
        <td>Input</td>
        <td>real</td>
        <td>0.1</td>
        <td><p>minimum auxiliary power (MW) (<code>constraint equation 40</code>)</p></td>
    </tr>
    
    <tr>
        <td>beta_poloidal_max</td>
        <td>Input</td>
        <td>real</td>
        <td>0.19</td>
        <td><p>maximum poloidal beta (<code>constraint equation 48</code>)</p></td>
    </tr>
    
    <tr>
        <td>bigqmin</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>minimum fusion gain Q (<code>constraint equation 28</code>)</p></td>
    </tr>
    
    <tr>
        <td>bmxlim</td>
        <td>Input</td>
        <td>real</td>
        <td>12.0</td>
        <td><p>maximum peak toroidal field (T) (<code>constraint equation 25</code>)</p></td>
    </tr>
    
    <tr>
        <td>fauxmn</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for minimum auxiliary power (<code>constraint equation 40</code>, <code>iteration variable 64</code>)</p></td>
    </tr>
    
    <tr>
        <td>fbeta_poloidal_eps</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for epsilon beta-poloidal (<code>constraint equation 6</code>, <code>iteration variable 8</code>)</p></td>
    </tr>
    
    <tr>
        <td>fbeta_poloidal</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for poloidal beta (<code>constraint equation 48</code>, <code>iteration variable 79</code>)</p></td>
    </tr>
    
    <tr>
        <td>fbeta_max</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for beta limit (<code>constraint equation 24</code>, <code>iteration variable 36</code>)</p></td>
    </tr>
    
    <tr>
        <td>fbeta_min</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for (lower) beta limit (<code>constraint equation 84</code>, <code>iteration variable 173</code>)</p></td>
    </tr>
    
    <tr>
        <td>fcpttf</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for TF coil current per turn upper limit
 (<code>constraint equation 77</code>, <code>iteration variable 146</code>)</p></td>
    </tr>
    
    <tr>
        <td>fr_conducting_wall</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for conducting wall radius / rminor limit
 (<code>constraint equation 23</code>, <code>iteration variable 104</code>)</p></td>
    </tr>
    
    <tr>
        <td>fdene</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for density limit (<code>constraint equation 5</code>, <code>iteration variable 9</code>)
 (invalid if <code>ipedestal=3</code>)</p></td>
    </tr>
    
    <tr>
        <td>fdivcol</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for divertor collisionality (<code>constraint equation 22</code>, <code>iteration variable 34</code>)</p></td>
    </tr>
    
    <tr>
        <td>fdtmp</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for first wall coolant temperature rise
 (<code>constraint equation 38</code>, <code>iteration variable 62</code>)</p></td>
    </tr>
    
    <tr>
        <td>fecrh_ignition</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for ecrh ignition constraint
 (<code>constraint equation 91</code>, <code>iteration variable 168</code>)</p></td>
    </tr>
    
    <tr>
        <td>fflutf</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for neutron fluence on TF coil (<code>constraint equation 53</code>, <code>iteration variable 92</code>)</p></td>
    </tr>
    
    <tr>
        <td>ffuspow</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum fusion power (<code>constraint equation 9</code>, <code>iteration variable 26</code>)</p></td>
    </tr>
    
    <tr>
        <td>fgamcd</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for current drive gamma (<code>constraint equation 37</code>, <code>iteration variable 40</code>)</p></td>
    </tr>
    
    <tr>
        <td>fhldiv</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for divertor heat load (<code>constraint equation 18</code>, <code>iteration variable 27</code>)</p></td>
    </tr>
    
    <tr>
        <td>fiooic</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>f-value for TF coil operating current / critical current ratio
 (<code>constraint equation 33</code>, <code>iteration variable 50</code>)</p></td>
    </tr>
    
    <tr>
        <td>fipir</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for Ip/Irod upper limit
 constraint equation icc = 46
 iteration variable ixc = 72</p></td>
    </tr>
    
    <tr>
        <td>fjohc</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for central solenoid current at end-of-flattop
 (<code>constraint equation 26</code>, <code>iteration variable 38</code>)</p></td>
    </tr>
    
    <tr>
        <td>fjohc0</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for central solenoid current at beginning of pulse
 (<code>constraint equation 27</code>, <code>iteration variable 39</code>)</p></td>
    </tr>
    
    <tr>
        <td>fjprot</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for TF coil winding pack current density
 (<code>constraint equation 35</code>, <code>iteration variable 53</code>)</p></td>
    </tr>
    
    <tr>
        <td>fl_h_threshold</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for L-H power threshold (<code>constraint equation 15</code>, <code>iteration variable 103</code>)</p></td>
    </tr>
    
    <tr>
        <td>fmva</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum MVA (<code>constraint equation 19</code>, <code>iteration variable 30</code>)</p></td>
    </tr>
    
    <tr>
        <td>fnbshinef</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum neutral beam shine-through fraction
 (<code>constraint equation 59</code>, <code>iteration variable 105</code>)</p></td>
    </tr>
    
    <tr>
        <td>fncycle</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for minimum CS coil stress load cycles
 (<code>constraint equation 90</code>, <code>iteration variable 167</code>)</p></td>
    </tr>
    
    <tr>
        <td>fnesep</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for Eich critical separatrix density
 (<code>constraint equation 76</code>, <code>iteration variable 144</code>)</p></td>
    </tr>
    
    <tr>
        <td>foh_stress</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for Tresca yield criterion in Central Solenoid
 (<code>constraint equation 72</code>, <code>iteration variable 123</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpeakb</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum toroidal field (<code>constraint equation 25</code>, <code>iteration variable 35</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpinj</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for injection power (<code>constraint equation 30</code>, <code>iteration variable 46</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpnetel</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for net electric power (<code>constraint equation 16</code>, <code>iteration variable 25</code>)</p></td>
    </tr>
    
    <tr>
        <td>fportsz</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for neutral beam tangency radius limit
 (<code>constraint equation 20</code>, <code>iteration variable 33</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpsepbqar</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum Psep*Bt/qAR limit (<code>constraint equation 68</code>, <code>iteration variable 117</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpsepr</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum Psep/R limit (<code>constraint equation 56</code>, <code>iteration variable 97</code>)</p></td>
    </tr>
    
    <tr>
        <td>fptemp</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for peak centrepost temperature (<code>constraint equation 44</code>, <code>iteration variable 68</code>)</p></td>
    </tr>
    
    <tr>
        <td>fptfnuc</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum TF coil nuclear heating (<code>constraint equation 54</code>, <code>iteration variable 95</code>)</p></td>
    </tr>
    
    <tr>
        <td>fq</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for edge safety factor (<code>constraint equation 45</code>, <code>iteration variable 71</code>)</p></td>
    </tr>
    
    <tr>
        <td>fqval</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for Q (<code>constraint equation 28</code>, <code>iteration variable 45</code>)</p></td>
    </tr>
    
    <tr>
        <td>fradpwr</td>
        <td>Input</td>
        <td>real</td>
        <td>0.99</td>
        <td><p>f-value for core radiation power limit (<code>constraint equation 17</code>, <code>iteration variable 28</code>)</p></td>
    </tr>
    
    <tr>
        <td>fradwall</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for upper limit on radiation wall load (<code>constr. equ. 67</code>, <code>iteration variable 116</code>)</p></td>
    </tr>
    
    <tr>
        <td>freinke</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for Reinke detachment criterion (<code>constr. equ. 78</code>, <code>iteration variable 147</code>)</p></td>
    </tr>
    
    <tr>
        <td>frminor</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for minor radius limit (<code>constraint equation 21</code>, <code>iteration variable 32</code>)</p></td>
    </tr>
    
    <tr>
        <td>fstrcase</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum TF coil case Tresca yield criterion
 (<code>constraint equation 31</code>, <code>iteration variable 48</code>)</p></td>
    </tr>
    
    <tr>
        <td>fstrcond</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maxiumum TF coil conduit Tresca yield criterion
 (<code>constraint equation 32</code>, <code>iteration variable 49</code>)</p></td>
    </tr>
    
    <tr>
        <td>fstr_wp</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maxiumum TF coil strain absolute value
 (<code>constraint equation 88</code>, <code>iteration variable 165</code>)</p></td>
    </tr>
    
    <tr>
        <td>fmaxvvstress</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum permitted stress of the VV
 (<code>constraint equation 65</code>, <code>iteration variable 113</code>)</p></td>
    </tr>
    
    <tr>
        <td>ftbr</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for minimum tritium breeding ratio (<code>constraint equation 52</code>, <code>iteration variable 89</code>)</p></td>
    </tr>
    
    <tr>
        <td>ft_burn</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for minimum burn time (<code>constraint equation 13</code>, <code>iteration variable 21</code>)</p></td>
    </tr>
    
    <tr>
        <td>ftcycl</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for cycle time (<code>constraint equation 42</code>, <code>iteration variable 67</code>)</p></td>
    </tr>
    
    <tr>
        <td>ftmargoh</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for central solenoid temperature margin
 (<code>constraint equation 60</code>, <code>iteration variable 106</code>)</p></td>
    </tr>
    
    <tr>
        <td>ftmargtf</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for TF coil temperature margin (<code>constraint equation 36</code>, <code>iteration variable 54</code>)</p></td>
    </tr>
    
    <tr>
        <td>ft_current_ramp_up</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for plasma current ramp-up time (<code>constraint equation 41</code>, <code>iteration variable 66</code>)</p></td>
    </tr>
    
    <tr>
        <td>ftpeak</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for first wall peak temperature (<code>constraint equation 39</code>, <code>iteration variable 63</code>)</p></td>
    </tr>
    
    <tr>
        <td>fvdump</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for dump voltage (<code>constraint equation 34</code>, <code>iteration variable 51</code>)</p></td>
    </tr>
    
    <tr>
        <td>fvs</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for flux-swing (V-s) requirement (STEADY STATE)
 (<code>constraint equation 12</code>, <code>iteration variable 15</code>)</p></td>
    </tr>
    
    <tr>
        <td>fvvhe</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for vacuum vessel He concentration limit (<code>i_blanket_type = 2</code>)
 (<code>constraint equation 55</code>, <code>iteration variable 96</code>)</p></td>
    </tr>
    
    <tr>
        <td>fwalld</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum wall load (<code>constraint equation 8</code>, <code>iteration variable 14</code>)</p></td>
    </tr>
    
    <tr>
        <td>fzeffmax</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum zeff (<code>constraint equation 64</code>, <code>iteration variable 112</code>)</p></td>
    </tr>
    
    <tr>
        <td>gammax</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>maximum current drive gamma (<code>constraint equation 37</code>)</p></td>
    </tr>
    
    <tr>
        <td>maxradwallload</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>Maximum permitted radiation wall load (MW/m^2) (<code>constraint equation 67</code>)</p></td>
    </tr>
    
    <tr>
        <td>mvalim</td>
        <td>Input</td>
        <td>real</td>
        <td>40.0</td>
        <td><p>maximum MVA limit (<code>constraint equation 19</code>)</p></td>
    </tr>
    
    <tr>
        <td>nbshinefmax</td>
        <td>Input</td>
        <td>real</td>
        <td>0.001</td>
        <td><p>maximum neutral beam shine-through fraction (<code>constraint equation 59</code>)</p></td>
    </tr>
    
    <tr>
        <td>nflutfmax</td>
        <td>Input</td>
        <td>real</td>
        <td>1e+23</td>
        <td><p>max fast neutron fluence on TF coil (n/m2) (<code>blktmodel&gt;0</code>) (<code>constraint equation 53</code>)
 Also used for demontable magnets (itart = 1) and superconducting coils (i_tf_sup = 1)
 To set the CP lifetime (<code>constraint equation 85</code>)</p></td>
    </tr>
    
    <tr>
        <td>pdivtlim</td>
        <td>Input</td>
        <td>real</td>
        <td>150.0</td>
        <td><p>Minimum pdivt [MW] (<code>constraint equation 80</code>)</p></td>
    </tr>
    
    <tr>
        <td>f_fw_rad_max</td>
        <td>Input</td>
        <td>real</td>
        <td>3.33</td>
        <td><p>peaking factor for radiation wall load (<code>constraint equation 67</code>)</p></td>
    </tr>
    
    <tr>
        <td>pflux_fw_rad_max_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Peak radiation wall load (MW/m^2) (<code>constraint equation 67</code>)</p></td>
    </tr>
    
    <tr>
        <td>pnetelin</td>
        <td>Input</td>
        <td>real</td>
        <td>1000.0</td>
        <td><p>required net electric power (MW) (<code>constraint equation 16</code>)</p></td>
    </tr>
    
    <tr>
        <td>powfmax</td>
        <td>Input</td>
        <td>real</td>
        <td>1500.0</td>
        <td><p>maximum fusion power (MW) (<code>constraint equation 9</code>)</p></td>
    </tr>
    
    <tr>
        <td>psepbqarmax</td>
        <td>Input</td>
        <td>real</td>
        <td>9.5</td>
        <td><p>maximum ratio of Psep*Bt/qAR (MWT/m) (<code>constraint equation 68</code>)</p></td>
    </tr>
    
    <tr>
        <td>pseprmax</td>
        <td>Input</td>
        <td>real</td>
        <td>25.0</td>
        <td><p>maximum ratio of power crossing the separatrix to plasma major radius (Psep/R) (MW/m)
 (<code>constraint equation 56</code>)</p></td>
    </tr>
    
    <tr>
        <td>ptfnucmax</td>
        <td>Input</td>
        <td>real</td>
        <td>0.001</td>
        <td><p>maximum nuclear heating in TF coil (MW/m3) (<code>constraint equation 54</code>)</p></td>
    </tr>
    
    <tr>
        <td>tbrmin</td>
        <td>Input</td>
        <td>real</td>
        <td>1.1</td>
        <td><p>minimum tritium breeding ratio (<code>constraint equation 52</code>)</p></td>
    </tr>
    
    <tr>
        <td>t_burn_min</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>minimum burn time (s) (KE - no longer itv., see issue #706)</p></td>
    </tr>
    
    <tr>
        <td>tcycmn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>minimum cycle time (s) (<code>constraint equation 42</code>)</p></td>
    </tr>
    
    <tr>
        <td>t_current_ramp_up_min</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>minimum plasma current ramp-up time (s) (<code>constraint equation 41</code>)</p></td>
    </tr>
    
    <tr>
        <td>vvhealw</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>allowed maximum helium concentration in vacuum vessel at end of plant life (appm)
 (<code>i_blanket_type =2</code>) (<code>constraint equation 55</code>)</p></td>
    </tr>
    
    <tr>
        <td>walalw</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>allowable neutron wall-load (MW/m2) (<code>constraint equation 8</code>)</p></td>
    </tr>
    
    <tr>
        <td>f_alpha_energy_confinement_min</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>Lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy confinement
 times (<code>constraint equation 62</code>)</p></td>
    </tr>
    
    <tr>
        <td>falpha_energy_confinement</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for lower limit on f_alpha_energy_confinement the ratio of alpha particle to energy
 confinement times (<code>constraint equation 62</code>, <code>iteration variable 110</code>)</p></td>
    </tr>
    
    <tr>
        <td>fniterpump</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for constraint that number of pumps &lt; tfno
 (<code>constraint equation 63</code>, <code>iteration variable 111</code>)</p></td>
    </tr>
    
    <tr>
        <td>zeffmax</td>
        <td>Input</td>
        <td>real</td>
        <td>3.6</td>
        <td><p>maximum value for Zeff (<code>constraint equation 64</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpoloidalpower</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for constraint on rate of change of energy in poloidal field
 (<code>constraint equation 66</code>, <code>iteration variable 115</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpsep</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value to ensure separatrix power is less than value from Kallenbach divertor
 (Not required as constraint 69 is an equality)</p></td>
    </tr>
    
    <tr>
        <td>fcqt</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>TF coil quench temparature remains below tmax_croco
 (<code>constraint equation 74</code>, <code>iteration variable 141</code>)</p></td>
    </tr>
    
</table>

## cost_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>abktflnc</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>allowable first wall/blanket neutron fluence (MW-yr/m2) (<code>blktmodel=0</code>)</p></td>
    </tr>
    
    <tr>
        <td>adivflnc</td>
        <td>Input</td>
        <td>real</td>
        <td>7.0</td>
        <td><p>allowable divertor heat fluence (MW-yr/m2)</p></td>
    </tr>
    
    <tr>
        <td>blkcst</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>blanket direct cost (M$)</p></td>
    </tr>
    
    <tr>
        <td>c221</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total account 221 cost (M$) - first wall, blanket, shield, support structure and div plates</p></td>
    </tr>
    
    <tr>
        <td>c222</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total account 222 cost (M$) - TF coils + PF coils</p></td>
    </tr>
    
    <tr>
        <td>capcost</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total capital cost including interest (M$)</p></td>
    </tr>
    
    <tr>
        <td>cconfix</td>
        <td>Input</td>
        <td>real</td>
        <td>80.0</td>
        <td><p>fixed cost of superconducting cable ($/m)</p></td>
    </tr>
    
    <tr>
        <td>cconshpf</td>
        <td>Input</td>
        <td>real</td>
        <td>70.0</td>
        <td><p>cost of PF coil steel conduit/sheath ($/m)</p></td>
    </tr>
    
    <tr>
        <td>cconshtf</td>
        <td>Input</td>
        <td>real</td>
        <td>75.0</td>
        <td><p>cost of TF coil steel conduit/sheath ($/m)</p></td>
    </tr>
    
    <tr>
        <td>cdcost</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>current drive direct costs (M$)</p></td>
    </tr>
    
    <tr>
        <td>cdirt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total plant direct cost (M$)</p></td>
    </tr>
    
    <tr>
        <td>cdrlife</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Full power year lifetime of heating/current drive system (y)</p></td>
    </tr>
    
    <tr>
        <td>cdrlife_cal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Calendar year lifetime of heating/current drive system (y)</p></td>
    </tr>
    
    <tr>
        <td>cfactr</td>
        <td>Input</td>
        <td>real</td>
        <td>0.75</td>
        <td><p>Total plant availability fraction; input if <code>iavail=0</code></p></td>
    </tr>
    
    <tr>
        <td>cpfact</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total plant capacity factor</p></td>
    </tr>
    
    <tr>
        <td>cfind</td>
        <td>Input</td>
        <td>real</td>
        <td>[0.244 0.244 0.244 0.29 ]</td>
        <td><p>indirect cost factor (func of lsa) (cost model = 0)</p></td>
    </tr>
    
    <tr>
        <td>cland</td>
        <td>Input</td>
        <td>real</td>
        <td>19.2</td>
        <td><p>cost of land (M$)</p></td>
    </tr>
    
    <tr>
        <td>coe</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>cost of electricity ($/MW-hr)</p></td>
    </tr>
    
    <tr>
        <td>coecap</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>capital cost of electricity (m$/kW-hr)</p></td>
    </tr>
    
    <tr>
        <td>coefuelt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>'fuel' (including replaceable components) contribution to cost of electricity (m$/kW-hr)</p></td>
    </tr>
    
    <tr>
        <td>coeoam</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>operation and maintenance contribution to cost of electricity (m$/kW-hr)</p></td>
    </tr>
    
    <tr>
        <td>concost</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plant construction cost (M$)</p></td>
    </tr>
    
    <tr>
        <td>costexp</td>
        <td>Input</td>
        <td>real</td>
        <td>0.8</td>
        <td><p>cost exponent for scaling in 2015 costs model</p></td>
    </tr>
    
    <tr>
        <td>costexp_pebbles</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6</td>
        <td><p>cost exponent for pebbles in 2015 costs model</p></td>
    </tr>
    
    <tr>
        <td>cost_factor_buildings</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>cost scaling factor for buildings</p></td>
    </tr>
    
    <tr>
        <td>cost_factor_land</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>cost scaling factor for land</p></td>
    </tr>
    
    <tr>
        <td>cost_factor_tf_coils</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>cost scaling factor for TF coils</p></td>
    </tr>
    
    <tr>
        <td>cost_factor_fwbs</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>cost scaling factor for fwbs</p></td>
    </tr>
    
    <tr>
        <td>cost_factor_rh</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>cost scaling factor for remote handling</p></td>
    </tr>
    
    <tr>
        <td>cost_factor_vv</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>cost scaling factor for vacuum vessel</p></td>
    </tr>
    
    <tr>
        <td>cost_factor_bop</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>cost scaling factor for energy conversion system</p></td>
    </tr>
    
    <tr>
        <td>cost_factor_misc</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>cost scaling factor for remaining subsystems</p></td>
    </tr>
    
    <tr>
        <td>maintenance_fwbs</td>
        <td>Input</td>
        <td>real</td>
        <td>0.2</td>
        <td><p>Maintenance cost factor: first wall, blanket, shield, divertor</p></td>
    </tr>
    
    <tr>
        <td>maintenance_gen</td>
        <td>Input</td>
        <td>real</td>
        <td>0.05</td>
        <td><p>Maintenance cost factor: All other components except coils, vacuum vessel,
 thermal shield, cryostat, land</p></td>
    </tr>
    
    <tr>
        <td>amortization</td>
        <td>Input</td>
        <td>real</td>
        <td>13.6</td>
        <td><p>amortization factor (fixed charge factor) "A" (years)</p></td>
    </tr>
    
    <tr>
        <td>cost_model</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for cost model:</p>
<ul>
<li>=0 use $ 1990 PROCESS model</li>
<li>=1 use $ 2014 Kovari model</li>
<li>=2 use user-provided model</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_cp_lifetime</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for the centrepost lifetime constraint
  0 : The CP full power year lifetime is set by the user via cplife_input
  1 : The CP lifetime is equal to the divertor lifetime
  2 : The CP lifetime is equal to the breeding blankets lifetime
  3 : The CP lifetime is equal to the plant lifetime</p></td>
    </tr>
    
    <tr>
        <td>cowner</td>
        <td>Input</td>
        <td>real</td>
        <td>0.15</td>
        <td><p>owner cost factor</p></td>
    </tr>
    
    <tr>
        <td>cplife_input</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>User input full power year lifetime of the centrepost (years) (i_cp_lifetime = 0)</p></td>
    </tr>
    
    <tr>
        <td>cplife</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Calculated full power year lifetime of centrepost (years)</p></td>
    </tr>
    
    <tr>
        <td>cplife_cal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Calculated calendar year lifetime of centrepost (years)</p></td>
    </tr>
    
    <tr>
        <td>cpstcst</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ST centrepost direct cost (M$)</p></td>
    </tr>
    
    <tr>
        <td>cpstflnc</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>allowable ST centrepost neutron fluence (MW-yr/m2)</p></td>
    </tr>
    
    <tr>
        <td>crctcore</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>reactor core costs (categories 221, 222 and 223)</p></td>
    </tr>
    
    <tr>
        <td>csi</td>
        <td>Input</td>
        <td>real</td>
        <td>16.0</td>
        <td><p>allowance for site costs (M$)</p></td>
    </tr>
    
    <tr>
        <td>cturbb</td>
        <td>Input</td>
        <td>real</td>
        <td>38.0</td>
        <td><p>cost of turbine building (M$)</p></td>
    </tr>
    
    <tr>
        <td>decomf</td>
        <td>Input</td>
        <td>real</td>
        <td>0.1</td>
        <td><p>proportion of constructed cost required for decommissioning fund</p></td>
    </tr>
    
    <tr>
        <td>dintrt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>diff between borrowing and saving interest rates</p></td>
    </tr>
    
    <tr>
        <td>divcst</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>divertor direct cost (M$)</p></td>
    </tr>
    
    <tr>
        <td>divlife</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Full power lifetime of divertor (y)</p></td>
    </tr>
    
    <tr>
        <td>divlife_cal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Calendar year lifetime of divertor (y)</p></td>
    </tr>
    
    <tr>
        <td>dtlife</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>period prior to the end of the plant life that the decommissioning fund is used (years)</p></td>
    </tr>
    
    <tr>
        <td>fcap0</td>
        <td>Input</td>
        <td>real</td>
        <td>1.165</td>
        <td><p>average cost of money for construction of plant assuming design/construction time of six years</p></td>
    </tr>
    
    <tr>
        <td>fcap0cp</td>
        <td>Input</td>
        <td>real</td>
        <td>1.08</td>
        <td><p>average cost of money for replaceable components assuming lead time for these of two years</p></td>
    </tr>
    
    <tr>
        <td>fcdfuel</td>
        <td>Input</td>
        <td>real</td>
        <td>0.1</td>
        <td><p>fraction of current drive cost treated as fuel (if <code>ifueltyp = 1</code>)</p></td>
    </tr>
    
    <tr>
        <td>fcontng</td>
        <td>Input</td>
        <td>real</td>
        <td>0.195</td>
        <td><p>project contingency factor</p></td>
    </tr>
    
    <tr>
        <td>fcr0</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0966</td>
        <td><p>fixed charge rate during construction</p></td>
    </tr>
    
    <tr>
        <td>fkind</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>multiplier for Nth of a kind costs</p></td>
    </tr>
    
    <tr>
        <td>fwallcst</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>first wall cost (M$)</p></td>
    </tr>
    
    <tr>
        <td>iavail</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Switch for plant availability model:</p>
<ul>
<li>=0 use input value for cfactr</li>
<li>=1 calculate cfactr using Taylor and Ward 1999 model</li>
<li>=2 calculate cfactr using new (2015) model</li>
<li>=3 calculate cfactr using ST model</li>
</ul></td>
    </tr>
    
    <tr>
        <td>ibkt_life</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for fw/blanket lifetime calculation in availability module:</p>
<ul>
<li>=0 use neutron fluence model</li>
<li>=1 use fusion power model (DEMO only)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>life_dpa</td>
        <td>Input</td>
        <td>real</td>
        <td>50.0</td>
        <td><p>Allowable DPA from DEMO fw/blanket lifetime calculation in availability module</p></td>
    </tr>
    
    <tr>
        <td>bktcycles</td>
        <td>Input</td>
        <td>real</td>
        <td>1000.0</td>
        <td><p>Number of fusion cycles to reach allowable DPA from DEMO fw/blanket lifetime calculation</p></td>
    </tr>
    
    <tr>
        <td>avail_min</td>
        <td>Input</td>
        <td>real</td>
        <td>0.75</td>
        <td><p>Minimum availability (<code>constraint equation 61</code>)</p></td>
    </tr>
    
    <tr>
        <td>tok_build_cost_per_vol</td>
        <td>Input</td>
        <td>real</td>
        <td>1283.0</td>
        <td><p>Unit cost for tokamak complex buildings, including building and site services ($/m3)</p></td>
    </tr>
    
    <tr>
        <td>light_build_cost_per_vol</td>
        <td>Input</td>
        <td>real</td>
        <td>270.0</td>
        <td><p>Unit cost for unshielded non-active buildings ($/m3)</p></td>
    </tr>
    
    <tr>
        <td>favail</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>F-value for minimum availability (<code>constraint equation 61</code>)</p></td>
    </tr>
    
    <tr>
        <td>num_rh_systems</td>
        <td>Input</td>
        <td>integer</td>
        <td>4</td>
        <td><p>Number of remote handling systems (1-10)</p></td>
    </tr>
    
    <tr>
        <td>conf_mag</td>
        <td>Input</td>
        <td>real</td>
        <td>0.99</td>
        <td><p>c parameter, which determines the temperature margin at which magnet lifetime starts to decline</p></td>
    </tr>
    
    <tr>
        <td>div_prob_fail</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0002</td>
        <td><p>Divertor probability of failure (per op day)</p></td>
    </tr>
    
    <tr>
        <td>div_umain_time</td>
        <td>Input</td>
        <td>real</td>
        <td>0.25</td>
        <td><p>Divertor unplanned maintenance time (years)</p></td>
    </tr>
    
    <tr>
        <td>div_nref</td>
        <td>Input</td>
        <td>real</td>
        <td>7000.0</td>
        <td><p>Reference value for cycle cycle life of divertor</p></td>
    </tr>
    
    <tr>
        <td>div_nu</td>
        <td>Input</td>
        <td>real</td>
        <td>14000.0</td>
        <td><p>The cycle when the divertor fails with 100% probability</p></td>
    </tr>
    
    <tr>
        <td>fwbs_nref</td>
        <td>Input</td>
        <td>real</td>
        <td>20000.0</td>
        <td><p>Reference value for cycle life of blanket</p></td>
    </tr>
    
    <tr>
        <td>fwbs_nu</td>
        <td>Input</td>
        <td>real</td>
        <td>40000.0</td>
        <td><p>The cycle when the blanket fails with 100% probability</p></td>
    </tr>
    
    <tr>
        <td>fwbs_prob_fail</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0002</td>
        <td><p>Fwbs probability of failure (per op day)</p></td>
    </tr>
    
    <tr>
        <td>fwbs_umain_time</td>
        <td>Input</td>
        <td>real</td>
        <td>0.25</td>
        <td><p>Fwbs unplanned maintenance time (years)</p></td>
    </tr>
    
    <tr>
        <td>redun_vacp</td>
        <td>Input</td>
        <td>real</td>
        <td>25.0</td>
        <td><p>Vacuum system pump redundancy level (%)</p></td>
    </tr>
    
    <tr>
        <td>redun_vac</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Number of redundant vacuum pumps</p></td>
    </tr>
    
    <tr>
        <td>t_operation</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Operational time (yrs)</p></td>
    </tr>
    
    <tr>
        <td>tbktrepl</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>time taken to replace blanket (y) (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>tcomrepl</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>time taken to replace both blanket and divertor (y) (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>tdivrepl</td>
        <td>Input</td>
        <td>real</td>
        <td>0.25</td>
        <td><p>time taken to replace divertor (y) (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>uubop</td>
        <td>Input</td>
        <td>real</td>
        <td>0.02</td>
        <td><p>unplanned unavailability factor for balance of plant (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>uucd</td>
        <td>Input</td>
        <td>real</td>
        <td>0.02</td>
        <td><p>unplanned unavailability factor for current drive (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>uudiv</td>
        <td>Input</td>
        <td>real</td>
        <td>0.04</td>
        <td><p>unplanned unavailability factor for divertor (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>uufuel</td>
        <td>Input</td>
        <td>real</td>
        <td>0.02</td>
        <td><p>unplanned unavailability factor for fuel system (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>uufw</td>
        <td>Input</td>
        <td>real</td>
        <td>0.04</td>
        <td><p>unplanned unavailability factor for first wall (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>uumag</td>
        <td>Input</td>
        <td>real</td>
        <td>0.02</td>
        <td><p>unplanned unavailability factor for magnets (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>uuves</td>
        <td>Input</td>
        <td>real</td>
        <td>0.04</td>
        <td><p>unplanned unavailability factor for vessel (<code>iavail=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>ifueltyp</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for fuel type:</p>
<ul>
<li>=2 treat initial blanket, divertor, first wall
   as capital costs. Treat all later items and
   fraction fcdfuel of CD equipment as fuel costs</li>
<li>=1 treat blanket divertor, first wall and
   fraction fcdfuel of CD equipment as fuel cost</li>
<li>=0 treat these as capital cost</li>
</ul></td>
    </tr>
    
    <tr>
        <td>ipnet</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for net electric power calculation:</p>
<ul>
<li>=0 scale so that always &gt; 0</li>
<li>=1 let go &lt; 0 (no c-o-e)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>ireactor</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for net electric power and cost of electricity calculations:</p>
<ul>
<li>=0 do not calculate MW(electric) or c-o-e</li>
<li>=1 calculate MW(electric) and c-o-e</li>
</ul></td>
    </tr>
    
    <tr>
        <td>lsa</td>
        <td>Input</td>
        <td>integer</td>
        <td>4</td>
        <td><p>Level of safety assurance switch (generally, use 3 or 4):</p>
<ul>
<li>=1 truly passively safe plant</li>
<li>=2,3 in-between</li>
<li>=4 like current fission plant</li>
</ul></td>
    </tr>
    
    <tr>
        <td>moneyint</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>interest portion of capital cost (M$)</p></td>
    </tr>
    
    <tr>
        <td>output_costs</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for costs output:</p>
<ul>
<li>=0 do not write cost-related outputs to file</li>
<li>=1 write cost-related outputs to file</li>
</ul></td>
    </tr>
    
    <tr>
        <td>discount_rate</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0435</td>
        <td><p>effective cost of money in constant dollars</p></td>
    </tr>
    
    <tr>
        <td>startupratio</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>ratio of additional HCD power for start-up to flat-top operational requirements</p></td>
    </tr>
    
    <tr>
        <td>startuppwr</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>cost associated with additional HCD system power required on start-up ($)</p></td>
    </tr>
    
    <tr>
        <td>supercond_cost_model</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for superconductor cost model:</p>
<ul>
<li>=0 use $/kg</li>
<li>=1 use $/kAm</li>
</ul></td>
    </tr>
    
    <tr>
        <td>tlife</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>Full power year plant lifetime (years)</p></td>
    </tr>
    
    <tr>
        <td>tmain</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Maintenance time for replacing CP (years) (iavail = 3)</p></td>
    </tr>
    
    <tr>
        <td>u_unplanned_cp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>User-input CP unplanned unavailability (iavail = 3)</p></td>
    </tr>
    
    <tr>
        <td>ucad</td>
        <td>Parameter</td>
        <td>real</td>
        <td>180.0D0</td>
        <td><p>unit cost for administration buildings (M$/m3)</p></td>
    </tr>
    
    <tr>
        <td>ucaf</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.5D6</td>
        <td><p>unit cost for aux facility power equipment ($)</p></td>
    </tr>
    
    <tr>
        <td>ucahts</td>
        <td>Parameter</td>
        <td>real</td>
        <td>31.0D0</td>
        <td><p>unit cost for aux heat transport equipment ($/W**exphts)</p></td>
    </tr>
    
    <tr>
        <td>ucap</td>
        <td>Parameter</td>
        <td>real</td>
        <td>17.0D0</td>
        <td><p>unit cost of auxiliary transformer ($/kVA)</p></td>
    </tr>
    
    <tr>
        <td>ucblbe</td>
        <td>Input</td>
        <td>real</td>
        <td>260.0</td>
        <td><p>unit cost for blanket beryllium ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucblbreed</td>
        <td>Input</td>
        <td>real</td>
        <td>875.0</td>
        <td><p>unit cost for breeder material ($/kg) (<code>blktmodel&gt;0</code>)</p></td>
    </tr>
    
    <tr>
        <td>ucblli</td>
        <td>Input</td>
        <td>real</td>
        <td>875.0</td>
        <td><p>unit cost for blanket lithium ($/kg) (30% Li6)</p></td>
    </tr>
    
    <tr>
        <td>ucblli2o</td>
        <td>Input</td>
        <td>real</td>
        <td>600.0</td>
        <td><p>unit cost for blanket Li_2O ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucbllipb</td>
        <td>Input</td>
        <td>real</td>
        <td>10.3</td>
        <td><p>unit cost for blanket Li-Pb ($/kg) (30% Li6)</p></td>
    </tr>
    
    <tr>
        <td>ucblss</td>
        <td>Input</td>
        <td>real</td>
        <td>90.0</td>
        <td><p>unit cost for blanket stainless steel ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucblvd</td>
        <td>Input</td>
        <td>real</td>
        <td>200.0</td>
        <td><p>unit cost for blanket vanadium ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucbpmp</td>
        <td>Parameter</td>
        <td>real</td>
        <td>2.925D5</td>
        <td><p>vacuum system backing pump cost ($)</p></td>
    </tr>
    
    <tr>
        <td>ucbus</td>
        <td>Input</td>
        <td>real</td>
        <td>0.123</td>
        <td><p>cost of aluminium bus for TF coil ($/A-m)</p></td>
    </tr>
    
    <tr>
        <td>uccase</td>
        <td>Input</td>
        <td>real</td>
        <td>50.0</td>
        <td><p>cost of superconductor case ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucco</td>
        <td>Parameter</td>
        <td>real</td>
        <td>350.0D0</td>
        <td><p>unit cost for control buildings (M$/m3)</p></td>
    </tr>
    
    <tr>
        <td>uccpcl1</td>
        <td>Input</td>
        <td>real</td>
        <td>250.0</td>
        <td><p>cost of high strength tapered copper ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>uccpclb</td>
        <td>Input</td>
        <td>real</td>
        <td>150.0</td>
        <td><p>cost of TF outboard leg plate coils ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>uccpmp</td>
        <td>Parameter</td>
        <td>real</td>
        <td>3.9D5</td>
        <td><p>vacuum system cryopump cost ($)</p></td>
    </tr>
    
    <tr>
        <td>uccr</td>
        <td>Parameter</td>
        <td>real</td>
        <td>460.0D0</td>
        <td><p>unit cost for cryogenic building (M$/vol)</p></td>
    </tr>
    
    <tr>
        <td>uccry</td>
        <td>Input</td>
        <td>real</td>
        <td>93000.0</td>
        <td><p>heat transport system cryoplant costs ($/W**expcry)</p></td>
    </tr>
    
    <tr>
        <td>uccryo</td>
        <td>Input</td>
        <td>real</td>
        <td>32.0</td>
        <td><p>unit cost for vacuum vessel ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>uccu</td>
        <td>Input</td>
        <td>real</td>
        <td>75.0</td>
        <td><p>unit cost for copper in superconducting cable ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucdgen</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.7D6</td>
        <td><p>cost per 8 MW diesel generator ($)</p></td>
    </tr>
    
    <tr>
        <td>ucdiv</td>
        <td>Input</td>
        <td>real</td>
        <td>280000.0</td>
        <td><p>cost of divertor blade ($)</p></td>
    </tr>
    
    <tr>
        <td>ucdtc</td>
        <td>Parameter</td>
        <td>real</td>
        <td>13.0D0</td>
        <td><p>detritiation, air cleanup cost ($/10000m3/hr)</p></td>
    </tr>
    
    <tr>
        <td>ucduct</td>
        <td>Parameter</td>
        <td>real</td>
        <td>4.225D4</td>
        <td><p>vacuum system duct cost ($/m)</p></td>
    </tr>
    
    <tr>
        <td>ucech</td>
        <td>Input</td>
        <td>real</td>
        <td>3.0</td>
        <td><p>ECH system cost ($/W)</p></td>
    </tr>
    
    <tr>
        <td>ucel</td>
        <td>Parameter</td>
        <td>real</td>
        <td>380.0D0</td>
        <td><p>unit cost for electrical equipment building (M$/m3)</p></td>
    </tr>
    
    <tr>
        <td>uces1</td>
        <td>Parameter</td>
        <td>real</td>
        <td>3.2D4</td>
        <td><p>MGF (motor-generator flywheel) cost factor ($/MVA**0.8)</p></td>
    </tr>
    
    <tr>
        <td>uces2</td>
        <td>Parameter</td>
        <td>real</td>
        <td>8.8D3</td>
        <td><p>MGF (motor-generator flywheel) cost factor ($/MJ**0.8)</p></td>
    </tr>
    
    <tr>
        <td>ucf1</td>
        <td>Input</td>
        <td>real</td>
        <td>22300000.0</td>
        <td><p>cost of fuelling system ($)</p></td>
    </tr>
    
    <tr>
        <td>ucfnc</td>
        <td>Input</td>
        <td>real</td>
        <td>35.0</td>
        <td><p>outer PF coil fence support cost ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucfpr</td>
        <td>Parameter</td>
        <td>real</td>
        <td>4.4D7</td>
        <td><p>cost of 60g/day tritium processing unit ($)</p></td>
    </tr>
    
    <tr>
        <td>ucfuel</td>
        <td>Input</td>
        <td>real</td>
        <td>3.45</td>
        <td><p>unit cost of D-T fuel (M$/year/1200MW)</p></td>
    </tr>
    
    <tr>
        <td>ucfwa</td>
        <td>Parameter</td>
        <td>real</td>
        <td>6.0D4</td>
        <td><p>first wall armour cost ($/m2)</p></td>
    </tr>
    
    <tr>
        <td>ucfwps</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.0D7</td>
        <td><p>first wall passive stabiliser cost ($)</p></td>
    </tr>
    
    <tr>
        <td>ucfws</td>
        <td>Parameter</td>
        <td>real</td>
        <td>5.3D4</td>
        <td><p>first wall structure cost ($/m2)</p></td>
    </tr>
    
    <tr>
        <td>ucgss</td>
        <td>Parameter</td>
        <td>real</td>
        <td>35.0D0</td>
        <td><p>cost of reactor structure ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>uche3</td>
        <td>Input</td>
        <td>real</td>
        <td>1000000.0</td>
        <td><p>cost of helium-3 ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>uchrs</td>
        <td>Input</td>
        <td>real</td>
        <td>87900000.0</td>
        <td><p>cost of heat rejection system ($)</p></td>
    </tr>
    
    <tr>
        <td>uchts</td>
        <td>Input</td>
        <td>real</td>
        <td>[15.3 19.1]</td>
        <td><p>cost of heat transport system equipment per loop ($/W); dependent on coolant type (coolwh)</p></td>
    </tr>
    
    <tr>
        <td>uciac</td>
        <td>Input</td>
        <td>real</td>
        <td>150000000.0</td>
        <td><p>cost of instrumentation, control &amp; diagnostics ($)</p></td>
    </tr>
    
    <tr>
        <td>ucich</td>
        <td>Input</td>
        <td>real</td>
        <td>3.0</td>
        <td><p>ICH system cost ($/W)</p></td>
    </tr>
    
    <tr>
        <td>ucint</td>
        <td>Parameter</td>
        <td>real</td>
        <td>35.0D0</td>
        <td><p>superconductor intercoil structure cost ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>uclh</td>
        <td>Input</td>
        <td>real</td>
        <td>3.3</td>
        <td><p>lower hybrid system cost ($/W)</p></td>
    </tr>
    
    <tr>
        <td>uclv</td>
        <td>Parameter</td>
        <td>real</td>
        <td>16.0D0</td>
        <td><p>low voltage system cost ($/kVA)</p></td>
    </tr>
    
    <tr>
        <td>ucmb</td>
        <td>Parameter</td>
        <td>real</td>
        <td>260.0D0</td>
        <td><p>unit cost for reactor maintenance building (M$/m3)</p></td>
    </tr>
    
    <tr>
        <td>ucme</td>
        <td>Input</td>
        <td>real</td>
        <td>125000000.0</td>
        <td><p>cost of maintenance equipment ($)</p></td>
    </tr>
    
    <tr>
        <td>ucmisc</td>
        <td>Input</td>
        <td>real</td>
        <td>25000000.0</td>
        <td><p>miscellaneous plant allowance ($)</p></td>
    </tr>
    
    <tr>
        <td>ucnbi</td>
        <td>Input</td>
        <td>real</td>
        <td>3.3</td>
        <td><p>NBI system cost ($/W)</p></td>
    </tr>
    
    <tr>
        <td>ucnbv</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1000.0D0</td>
        <td><p>cost of nuclear building ventilation ($/m3)</p></td>
    </tr>
    
    <tr>
        <td>ucoam</td>
        <td>Input</td>
        <td>real</td>
        <td>[68.8 68.8 68.8 74.4]</td>
        <td><p>annual cost of operation and maintenance (M$/year/1200MW**0.5)</p></td>
    </tr>
    
    <tr>
        <td>ucpens</td>
        <td>Input</td>
        <td>real</td>
        <td>32.0</td>
        <td><p>penetration shield cost ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucpfb</td>
        <td>Input</td>
        <td>real</td>
        <td>210.0</td>
        <td><p>cost of PF coil buses ($/kA-m)</p></td>
    </tr>
    
    <tr>
        <td>ucpfbk</td>
        <td>Input</td>
        <td>real</td>
        <td>16600.0</td>
        <td><p>cost of PF coil DC breakers ($/MVA**0.7)</p></td>
    </tr>
    
    <tr>
        <td>ucpfbs</td>
        <td>Input</td>
        <td>real</td>
        <td>4900.0</td>
        <td><p>cost of PF burn power supplies ($/kW**0.7)</p></td>
    </tr>
    
    <tr>
        <td>ucpfcb</td>
        <td>Input</td>
        <td>real</td>
        <td>75000.0</td>
        <td><p>cost of PF coil AC breakers ($/circuit)</p></td>
    </tr>
    
    <tr>
        <td>ucpfdr1</td>
        <td>Input</td>
        <td>real</td>
        <td>150.0</td>
        <td><p>cost factor for dump resistors ($/MJ)</p></td>
    </tr>
    
    <tr>
        <td>ucpfic</td>
        <td>Input</td>
        <td>real</td>
        <td>10000.0</td>
        <td><p>cost of PF instrumentation and control ($/channel)</p></td>
    </tr>
    
    <tr>
        <td>ucpfps</td>
        <td>Input</td>
        <td>real</td>
        <td>35000.0</td>
        <td><p>cost of PF coil pulsed power supplies ($/MVA)</p></td>
    </tr>
    
    <tr>
        <td>ucphx</td>
        <td>Parameter</td>
        <td>real</td>
        <td>15.0D0</td>
        <td><p>primary heat transport cost ($/W**exphts)</p></td>
    </tr>
    
    <tr>
        <td>ucpp</td>
        <td>Parameter</td>
        <td>real</td>
        <td>48.0D0</td>
        <td><p>cost of primary power transformers ($/kVA**0.9)</p></td>
    </tr>
    
    <tr>
        <td>ucrb</td>
        <td>Input</td>
        <td>real</td>
        <td>400.0</td>
        <td><p>cost of reactor building (M$/m3)</p></td>
    </tr>
    
    <tr>
        <td>ucsc</td>
        <td>Input</td>
        <td>real</td>
        <td>[ 600.  600.  300.  600.  600.  600.  300. 1200. 1200.]</td>
        <td><p>cost of superconductor ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>sc_mat_cost_0</td>
        <td>Input</td>
        <td>real</td>
        <td>[ 4.8  2.   1.   4.8  4.8 47.4  1.  47.4 47.4]</td>
        <td><p>cost of superconductor ($/kA m) at 6.4 T, 4.2 K</p></td>
    </tr>
    
    <tr>
        <td>ucsh</td>
        <td>Parameter</td>
        <td>real</td>
        <td>115.0D0</td>
        <td><p>cost of shops and warehouses (M$/m3)</p></td>
    </tr>
    
    <tr>
        <td>ucshld</td>
        <td>Input</td>
        <td>real</td>
        <td>32.0</td>
        <td><p>cost of shield structural steel ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucswyd</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.84D7</td>
        <td><p>switchyard equipment costs ($)</p></td>
    </tr>
    
    <tr>
        <td>uctfbr</td>
        <td>Input</td>
        <td>real</td>
        <td>1.22</td>
        <td><p>cost of TF coil breakers ($/W**0.7)</p></td>
    </tr>
    
    <tr>
        <td>uctfbus</td>
        <td>Input</td>
        <td>real</td>
        <td>100.0</td>
        <td><p>cost of TF coil bus ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>uctfdr</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.75D-4</td>
        <td><p>cost of TF coil dump resistors ($/J)</p></td>
    </tr>
    
    <tr>
        <td>uctfgr</td>
        <td>Parameter</td>
        <td>real</td>
        <td>5000.0D0</td>
        <td><p>additional cost of TF coil dump resistors ($/coil)</p></td>
    </tr>
    
    <tr>
        <td>uctfic</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.0D4</td>
        <td><p>cost of TF coil instrumentation and control ($/coil/30)</p></td>
    </tr>
    
    <tr>
        <td>uctfps</td>
        <td>Input</td>
        <td>real</td>
        <td>24.0</td>
        <td><p>cost of TF coil power supplies ($/W**0.7)</p></td>
    </tr>
    
    <tr>
        <td>uctfsw</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>cost of TF coil slow dump switches ($/A)</p></td>
    </tr>
    
    <tr>
        <td>uctpmp</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.105D5</td>
        <td><p>cost of turbomolecular pump ($)</p></td>
    </tr>
    
    <tr>
        <td>uctr</td>
        <td>Parameter</td>
        <td>real</td>
        <td>370.0D0</td>
        <td><p>cost of tritium building ($/m3)</p></td>
    </tr>
    
    <tr>
        <td>ucturb</td>
        <td>Input</td>
        <td>real</td>
        <td>[2.30e+08 2.45e+08]</td>
        <td><p>cost of turbine plant equipment ($) (dependent on coolant type coolwh)</p></td>
    </tr>
    
    <tr>
        <td>ucvalv</td>
        <td>Parameter</td>
        <td>real</td>
        <td>3.9D5</td>
        <td><p>vacuum system valve cost ($)</p></td>
    </tr>
    
    <tr>
        <td>ucvdsh</td>
        <td>Parameter</td>
        <td>real</td>
        <td>26.0D0</td>
        <td><p>vacuum duct shield cost ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucviac</td>
        <td>Parameter</td>
        <td>real</td>
        <td>1.3D6</td>
        <td><p>vacuum system instrumentation and control cost ($)</p></td>
    </tr>
    
    <tr>
        <td>ucwindpf</td>
        <td>Input</td>
        <td>real</td>
        <td>465.0</td>
        <td><p>cost of PF coil superconductor windings ($/m)</p></td>
    </tr>
    
    <tr>
        <td>ucwindtf</td>
        <td>Input</td>
        <td>real</td>
        <td>480.0</td>
        <td><p>cost of TF coil superconductor windings ($/m)</p></td>
    </tr>
    
    <tr>
        <td>ucws</td>
        <td>Parameter</td>
        <td>real</td>
        <td>460.0D0</td>
        <td><p>cost of active assembly shop ($/m3)</p></td>
    </tr>
    
    <tr>
        <td>ucwst</td>
        <td>Input</td>
        <td>real</td>
        <td>[0.   3.94 5.91 7.88]</td>
        <td><p>cost of waste disposal (M$/y/1200MW)</p></td>
    </tr>
    
</table>

## current_drive_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>beamwd</td>
        <td>Input</td>
        <td>real</td>
        <td>0.58</td>
        <td><p>width of neutral beam duct where it passes between the TF coils (m)
 T Inoue et al, Design of neutral beam system for ITER-FEAT,
 <A HREF=http://dx.doi.org/10.1016/S0920-3796(01)00339-8>
 Fusion Engineering and Design, Volumes 56-57, October 2001, Pages 517-521</A>)</p></td>
    </tr>
    
    <tr>
        <td>bigq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Fusion gain; P_fusion / (P_injection + P_ohmic)</p></td>
    </tr>
    
    <tr>
        <td>bootstrap_current_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>bootstrap current fraction (enforced; see i_bootstrap_current)</p></td>
    </tr>
    
    <tr>
        <td>bootstrap_current_fraction_max</td>
        <td>Input</td>
        <td>real</td>
        <td>0.9</td>
        <td><p>maximum fraction of plasma current from bootstrap; if <code>bootstrap_current_fraction_max &lt; 0</code>,
 bootstrap fraction = abs(bootstrap_current_fraction_max)</p></td>
    </tr>
    
    <tr>
        <td>bscf_iter89</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>bootstrap current fraction, ITER 1989 model</p></td>
    </tr>
    
    <tr>
        <td>bscf_nevins</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>bootstrap current fraction, Nevins et al model</p></td>
    </tr>
    
    <tr>
        <td>bscf_sauter</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>bootstrap current fraction, Sauter et al model</p></td>
    </tr>
    
    <tr>
        <td>bscf_wilson</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>bootstrap current fraction, Wilson et al model</p></td>
    </tr>
    
    <tr>
        <td>bscf_sakai</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Bootstrap current fraction, Sakai et al model</p></td>
    </tr>
    
    <tr>
        <td>bscf_aries</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Bootstrap current fraction, ARIES model</p></td>
    </tr>
    
    <tr>
        <td>bscf_andrade</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Bootstrap current fraction, Andrade et al model</p></td>
    </tr>
    
    <tr>
        <td>bscf_hoang</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Bootstrap current fraction, Hoang et al model</p></td>
    </tr>
    
    <tr>
        <td>bscf_wong</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Bootstrap current fraction, Wong et al model</p></td>
    </tr>
    
    <tr>
        <td>bscf_gi_i</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Bootstrap current fraction, first Gi et al model</p></td>
    </tr>
    
    <tr>
        <td>bscf_gi_ii</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Bootstrap current fraction, second Gi et al model</p></td>
    </tr>
    
    <tr>
        <td>cboot</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>bootstrap current fraction multiplier</p></td>
    </tr>
    
    <tr>
        <td>beam_current</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam current (A)</p></td>
    </tr>
    
    <tr>
        <td>diacf_hender</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>diamagnetic current fraction, Hender fit</p></td>
    </tr>
    
    <tr>
        <td>diacf_scene</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>diamagnetic current fraction, SCENE fit</p></td>
    </tr>
    
    <tr>
        <td>diamagnetic_current_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>diamagnetic current fraction</p></td>
    </tr>
    
    <tr>
        <td>echpwr</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ECH power (MW)</p></td>
    </tr>
    
    <tr>
        <td>echwpow</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ECH wall plug power (MW)</p></td>
    </tr>
    
    <tr>
        <td>effcd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>current drive efficiency (A/W)</p></td>
    </tr>
    
    <tr>
        <td>harnum</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>cyclotron harmonic frequency number, used in cut-off function</p></td>
    </tr>
    
    <tr>
        <td>wave_mode</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for ECRH wave mode :</p>
<ul>
<li>=0 O-mode</li>
<li>=1 X-mode</li>
</ul></td>
    </tr>
    
    <tr>
        <td>beam_energy</td>
        <td>Input</td>
        <td>real</td>
        <td>1000.0</td>
        <td><p>neutral beam energy (keV) (<code>iteration variable 19</code>)</p></td>
    </tr>
    
    <tr>
        <td>etacd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>auxiliary power wall plug to injector efficiency</p></td>
    </tr>
    
    <tr>
        <td>etacdfix</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>secondary auxiliary power wall plug to injector efficiency</p></td>
    </tr>
    
    <tr>
        <td>etaech</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>ECH wall plug to injector efficiency</p></td>
    </tr>
    
    <tr>
        <td>etalh</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>lower hybrid wall plug to injector efficiency</p></td>
    </tr>
    
    <tr>
        <td>etanbi</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>neutral beam wall plug to injector efficiency</p></td>
    </tr>
    
    <tr>
        <td>fpion</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>fraction of beam energy to ions</p></td>
    </tr>
    
    <tr>
        <td>pnbitot</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam power entering vacuum vessel</p></td>
    </tr>
    
    <tr>
        <td>pscf_scene</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Pfirsch-Schlüter current fraction, SCENE fit</p></td>
    </tr>
    
    <tr>
        <td>nbshinemw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam shine-through power</p></td>
    </tr>
    
    <tr>
        <td>feffcd</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>current drive efficiency fudge factor (<code>iteration variable 47</code>)</p></td>
    </tr>
    
    <tr>
        <td>forbitloss</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fraction of neutral beam power lost after ionisation but before
 thermalisation (orbit loss fraction)</p></td>
    </tr>
    
    <tr>
        <td>frbeam</td>
        <td>Input</td>
        <td>real</td>
        <td>1.05</td>
        <td><p>R_tangential / R_major for neutral beam injection</p></td>
    </tr>
    
    <tr>
        <td>f_tritium_beam</td>
        <td>Input</td>
        <td>real</td>
        <td>1e-06</td>
        <td><p>fraction of beam that is tritium</p></td>
    </tr>
    
    <tr>
        <td>gamcd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>normalised current drive efficiency (1.0e20 A/(W m^2))</p></td>
    </tr>
    
    <tr>
        <td>gamma_ecrh</td>
        <td>Input</td>
        <td>real</td>
        <td>0.35</td>
        <td><p>User input ECRH gamma (1.0e20 A/(W m^2))</p></td>
    </tr>
    
    <tr>
        <td>xi_ebw</td>
        <td>Input</td>
        <td>real</td>
        <td>0.8</td>
        <td><p>User scaling input for EBW plasma heating. Default 0.43</p></td>
    </tr>
    
    <tr>
        <td>iefrf</td>
        <td>Input</td>
        <td>integer</td>
        <td>5</td>
        <td><p>Switch for current drive efficiency model:</p>
<ul>
<li>=1 Fenstermacher Lower Hybrid</li>
<li>=2 Ion Cyclotron current drive</li>
<li>=3 Fenstermacher ECH</li>
<li>=4 Ehst Lower Hybrid</li>
<li>=5 ITER Neutral Beam</li>
<li>=6 new Culham Lower Hybrid model</li>
<li>=7 new Culham ECCD model</li>
<li>=8 new Culham Neutral Beam model</li>
<li>=9 RFP option removed in PROCESS (issue #508)</li>
<li>=10 ECRH user input gamma</li>
<li>=11 ECRH "HARE" model (E. Poli, Physics of Plasmas 2019). Removed in #1811.</li>
<li>=12 EBW user scaling input. Scaling (S. Freethy)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>iefrffix</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for 2nd current drive efficiency model:</p>
<ul>
<li>=0 No fixed current drive</li>
<li>=1 Fenstermacher Lower Hybrid</li>
<li>=2 Ion Cyclotron current drive</li>
<li>=3 Fenstermacher ECH</li>
<li>=4 Ehst Lower Hybrid</li>
<li>=5 ITER Neutral Beam</li>
<li>=6 new Culham Lower Hybrid model</li>
<li>=7 new Culham ECCD model</li>
<li>=8 new Culham Neutral Beam model</li>
<li>=9 RFP option removed in PROCESS (issue #508)</li>
<li>=10 ECRH user input gamma</li>
<li>=11 ECRH "HARE" model (E. Poli, Physics of Plasmas 2019). Removed in #1811.</li>
<li>=12 EBW user scaling input. Scaling (S. Freethy)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>irfcd</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for current drive calculation:</p>
<ul>
<li>=0 turned off</li>
<li>=1 turned on</li>
</ul></td>
    </tr>
    
    <tr>
        <td>nbshinef</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam shine-through fraction</p></td>
    </tr>
    
    <tr>
        <td>nbshield</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>neutral beam duct shielding thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>pheat</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>heating power not used for current drive (MW) (<code>iteration variable 11</code>)</p></td>
    </tr>
    
    <tr>
        <td>pheatfix</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>secondary fixed heating power not used for current drive (MW)</p></td>
    </tr>
    
    <tr>
        <td>pinjalw</td>
        <td>Input</td>
        <td>real</td>
        <td>150.0</td>
        <td><p>maximum allowable value for injected power (MW) (<code>constraint equation 30</code>)</p></td>
    </tr>
    
    <tr>
        <td>pinjemw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>auxiliary injected power to electrons (MW)</p></td>
    </tr>
    
    <tr>
        <td>pinjimw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>auxiliary injected power to ions (MW)</p></td>
    </tr>
    
    <tr>
        <td>pinjmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total auxiliary injected power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pinjfixmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>secondary total fixed auxiliary injected power (MW)</p></td>
    </tr>
    
    <tr>
        <td>plasma_current_internal_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma current fraction driven internally (Bootstrap + Diamagnetic + PS)</p></td>
    </tr>
    
    <tr>
        <td>plhybd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>lower hybrid injection power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pnbeam</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam injection power (MW)</p></td>
    </tr>
    
    <tr>
        <td>porbitlossmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam power lost after ionisation but before thermalisation (orbit loss power) (MW)</p></td>
    </tr>
    
    <tr>
        <td>ps_current_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Pfirsch-Schlüter current fraction</p></td>
    </tr>
    
    <tr>
        <td>pwplh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>lower hybrid wall plug power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pwpnb</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam wall plug power (MW)</p></td>
    </tr>
    
    <tr>
        <td>rtanbeam</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam centreline tangency radius (m)</p></td>
    </tr>
    
    <tr>
        <td>rtanmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum tangency radius for centreline of beam (m)</p></td>
    </tr>
    
    <tr>
        <td>taubeam</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam e-decay lengths to plasma centre</p></td>
    </tr>
    
    <tr>
        <td>tbeamin</td>
        <td>Input</td>
        <td>real</td>
        <td>3.0</td>
        <td><p>permitted neutral beam e-decay lengths to plasma centre</p></td>
    </tr>
    
</table>

## dcll_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>r_fci</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial BZ thickness [m]</p></td>
    </tr>
    
    <tr>
        <td>r_backwall</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial BZ thickness [m]</p></td>
    </tr>
    
    <tr>
        <td>bz_r_ib</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Structure/coolant compositional fractions</p></td>
    </tr>
    
    <tr>
        <td>bz_r_ob</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Structure/coolant compositional fractions</p></td>
    </tr>
    
    <tr>
        <td>f_vol_stff_plates</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>MF/BSS compositional fractions</p></td>
    </tr>
    
    <tr>
        <td>f_vol_stl_bz_struct</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>MF/BSS compositional fractions</p></td>
    </tr>
    
    <tr>
        <td>f_vol_stl_back_wall</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>MF/BSS compositional fractions</p></td>
    </tr>
    
    <tr>
        <td>f_vol_stl_fw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>MF/BSS compositional fractions</p></td>
    </tr>
    
    <tr>
        <td>f_vol_mfbss_stl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Volume of FCIs, other BZ structure, liquid channels, backwall and MF/BSS [m^3]</p></td>
    </tr>
    
    <tr>
        <td>f_vol_mfbss_he</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Volume of FCIs, other BZ structure, liquid channels, backwall and MF/BSS [m^3]</p></td>
    </tr>
    
    <tr>
        <td>f_vol_mfbss_pbli</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Volume of FCIs, other BZ structure, liquid channels, backwall and MF/BSS [m^3]</p></td>
    </tr>
    
    <tr>
        <td>vol_fci</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>BZ masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>vol_bz_struct</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>BZ masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>vol_bz_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>BZ masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>vol_bz_liq_ib</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>BZ masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>vol_bz_liq_ob</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>BZ masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>vol_bw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>BZ masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>vol_bss</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>BZ masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>wht_cer</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Backwall masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>wht_stl_struct</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Backwall masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>wht_cool_struct</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Backwall masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>wht_bw_stl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>MF/BSS masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>wht_bw_cool</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>MF/BSS masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>wht_mfbss_stl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>FW masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>wht_mfbss_cool</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>FW masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>wht_mfbss_pbli</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>FW masses by composition [kg]</p></td>
    </tr>
    
    <tr>
        <td>fwmass_stl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total masses of material in blanket [kg]</p></td>
    </tr>
    
    <tr>
        <td>fwmass_cool</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total masses of material in blanket [kg]</p></td>
    </tr>
    
    <tr>
        <td>mass_cool_blanket</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total mass for an inboard/outboard reactor segment [kg]</p></td>
    </tr>
    
    <tr>
        <td>mass_liq_blanket</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total mass for an inboard/outboard reactor segment [kg]</p></td>
    </tr>
    
    <tr>
        <td>mass_stl_blanket</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total mass for an inboard/outboard reactor segment [kg]</p></td>
    </tr>
    
    <tr>
        <td>mass_segm_ib</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>mass_segm_ob</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## divertor_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>adas</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>area divertor / area main plasma (along separatrix)</p></td>
    </tr>
    
    <tr>
        <td>anginc</td>
        <td>Input</td>
        <td>real</td>
        <td>0.262</td>
        <td><p>angle of incidence of field line on plate (rad)</p></td>
    </tr>
    
    <tr>
        <td>beta_div</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>field line angle wrt divertor target plate (degrees)</p></td>
    </tr>
    
    <tr>
        <td>betai</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>poloidal plane angle between divertor plate and leg, inboard (rad)</p></td>
    </tr>
    
    <tr>
        <td>betao</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>poloidal plane angle between divertor plate and leg, outboard (rad)</p></td>
    </tr>
    
    <tr>
        <td>bpsout</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6</td>
        <td><p>reference B_p at outboard divertor strike point (T)</p></td>
    </tr>
    
    <tr>
        <td>c1div</td>
        <td>Input</td>
        <td>real</td>
        <td>0.45</td>
        <td><p>fitting coefficient to adjust ptpdiv, ppdiv</p></td>
    </tr>
    
    <tr>
        <td>c2div</td>
        <td>Input</td>
        <td>real</td>
        <td>-7.0</td>
        <td><p>fitting coefficient to adjust ptpdiv, ppdiv</p></td>
    </tr>
    
    <tr>
        <td>c3div</td>
        <td>Input</td>
        <td>real</td>
        <td>0.54</td>
        <td><p>fitting coefficient to adjust ptpdiv, ppdiv</p></td>
    </tr>
    
    <tr>
        <td>c4div</td>
        <td>Input</td>
        <td>real</td>
        <td>-3.6</td>
        <td><p>fitting coefficient to adjust ptpdiv, ppdiv</p></td>
    </tr>
    
    <tr>
        <td>c5div</td>
        <td>Input</td>
        <td>real</td>
        <td>0.7</td>
        <td><p>fitting coefficient to adjust ptpdiv, ppdiv</p></td>
    </tr>
    
    <tr>
        <td>c6div</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fitting coefficient to adjust ptpdiv, ppdiv</p></td>
    </tr>
    
    <tr>
        <td>delld</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>coeff for power distribution along main plasma</p></td>
    </tr>
    
    <tr>
        <td>dendiv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma density at divertor (10**20 /m3)</p></td>
    </tr>
    
    <tr>
        <td>densin</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>density at plate (on separatrix) (10**20 /m3)</p></td>
    </tr>
    
    <tr>
        <td>divclfr</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>divertor coolant fraction</p></td>
    </tr>
    
    <tr>
        <td>divdens</td>
        <td>Input</td>
        <td>real</td>
        <td>10000.0</td>
        <td><p>divertor structure density (kg/m3)</p></td>
    </tr>
    
    <tr>
        <td>divdum</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for divertor Zeff model:</p>
<ul>
<li>=0 calc</li>
<li>=1 input</li>
</ul></td>
    </tr>
    
    <tr>
        <td>divfix</td>
        <td>Input</td>
        <td>real</td>
        <td>0.2</td>
        <td><p>divertor structure vertical thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>divmas</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>divertor plate mass (kg)</p></td>
    </tr>
    
    <tr>
        <td>divplt</td>
        <td>Input</td>
        <td>real</td>
        <td>0.035</td>
        <td><p>divertor plate thickness (m) (from Spears, Sept 1990)</p></td>
    </tr>
    
    <tr>
        <td>divsur</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>divertor surface area (m2)</p></td>
    </tr>
    
    <tr>
        <td>fdfs</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>radial gradient ratio</p></td>
    </tr>
    
    <tr>
        <td>fdiva</td>
        <td>Input</td>
        <td>real</td>
        <td>1.11</td>
        <td><p>divertor area fudge factor (for ITER, Sept 1990)</p></td>
    </tr>
    
    <tr>
        <td>fhout</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fraction of power to outboard divertor (for single null)</p></td>
    </tr>
    
    <tr>
        <td>fififi</td>
        <td>Input</td>
        <td>real</td>
        <td>0.004</td>
        <td><p>coefficient for gamdiv</p></td>
    </tr>
    
    <tr>
        <td>flux_exp</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>The plasma flux expansion in the divertor (default 2; Wade 2020)</p></td>
    </tr>
    
    <tr>
        <td>frrp</td>
        <td>Input</td>
        <td>real</td>
        <td>0.4</td>
        <td><p>fraction of radiated power to plate</p></td>
    </tr>
    
    <tr>
        <td>hldiv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>divertor heat load (MW/m2)</p></td>
    </tr>
    
    <tr>
        <td>i_hldiv</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for user input hldiv:</p>
<ul>
<li>= 0: divtart model turned off and user inputs hldiv</li>
<li>= 1: divtart model calculates hldiv</li>
<li>= 2: divwade model calculates hldiv</li>
</ul></td>
    </tr>
    
    <tr>
        <td>hldivlim</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>heat load limit (MW/m2)</p></td>
    </tr>
    
    <tr>
        <td>ksic</td>
        <td>Input</td>
        <td>real</td>
        <td>0.8</td>
        <td><p>power fraction for outboard double-null scrape-off plasma</p></td>
    </tr>
    
    <tr>
        <td>lamp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>power flow width (m)</p></td>
    </tr>
    
    <tr>
        <td>minstang</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>minimum strike angle for heat flux calculation</p></td>
    </tr>
    
    <tr>
        <td>omegan</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>pressure ratio (nT)_plasma / (nT)_scrape-off</p></td>
    </tr>
    
    <tr>
        <td>omlarg</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>power spillage to private flux factor</p></td>
    </tr>
    
    <tr>
        <td>ppdivr</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>peak heat load at plate (with radiation) (MW/m2)</p></td>
    </tr>
    
    <tr>
        <td>prn1</td>
        <td>Input</td>
        <td>real</td>
        <td>0.285</td>
        <td><p>n-scrape-off / n-average plasma; (input for <code>ipedestal=0</code>, = nesep/dene if <code>ipedestal&gt;=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>ptpdiv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>peak temperature at the plate (eV)</p></td>
    </tr>
    
    <tr>
        <td>rconl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>connection length ratio, outboard side</p></td>
    </tr>
    
    <tr>
        <td>rlclolcn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ratio of collision length / connection length</p></td>
    </tr>
    
    <tr>
        <td>rlenmax</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>maximum value for length ratio (rlclolcn) (<code>constraintg eqn 22</code>)</p></td>
    </tr>
    
    <tr>
        <td>rsrd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>effective separatrix/divertor radius ratio</p></td>
    </tr>
    
    <tr>
        <td>tconl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>main plasma connection length (m)</p></td>
    </tr>
    
    <tr>
        <td>tdiv</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>temperature at divertor (eV) (input for stellarator only, calculated for tokamaks)</p></td>
    </tr>
    
    <tr>
        <td>tsep</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>temperature at the separatrix (eV)</p></td>
    </tr>
    
    <tr>
        <td>xparain</td>
        <td>Input</td>
        <td>real</td>
        <td>2100.0</td>
        <td><p>parallel heat transport coefficient (m2/s)</p></td>
    </tr>
    
    <tr>
        <td>xpertin</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>perpendicular heat transport coefficient (m2/s)</p></td>
    </tr>
    
    <tr>
        <td>zeffdiv</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>Zeff in the divertor region (if <code>divdum/=0</code>)</p></td>
    </tr>
    
</table>

## error_handling
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>errors_on</td>
        <td>Output</td>
        <td>logical</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error_okay*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error_info*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>1</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error_warn*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>2</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error_severe*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>3</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error_id*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td><p>error_id : identifier for final message encountered</p></td>
    </tr>
    
    <tr>
        <td>error_status</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>error_status : overall status flag for a run; on exit:<UL>
                 <LI> 0  all okay
                 <LI> 1  informational messages have been encountered
                 <LI> 2  warning (non-fatal) messages have been encountered
                 <LI> 3  severe (fatal) errors have occurred</UL></p></td>
    </tr>
    
    <tr>
        <td>int_default*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>-999999</td>
        <td></td>
    </tr>
    
    <tr>
        <td>flt_default*</td>
        <td>Parameter</td>
        <td>real</td>
        <td>real(INT_DEFAULT, kind(1.0D0))</td>
        <td></td>
    </tr>
    
    <tr>
        <td>idiags</td>
        <td>Input</td>
        <td>integer</td>
        <td>[-999999 -999999 -999999 -999999 -999999 -999999 -999999 -999999]</td>
        <td></td>
    </tr>
    
    <tr>
        <td>fdiags</td>
        <td>Input</td>
        <td>real</td>
        <td>[-999999. -999999. -999999. -999999. -999999. -999999. -999999. -999999.]</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error_head*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error_tail*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error_type*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## fson_string_m
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>block_size*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>32</td>
        <td></td>
    </tr>
    
</table>

## fson_value_m
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>type_unknown</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>-1</td>
        <td></td>
    </tr>
    
    <tr>
        <td>type_null</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>type_object</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>1</td>
        <td></td>
    </tr>
    
    <tr>
        <td>type_array</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>2</td>
        <td></td>
    </tr>
    
    <tr>
        <td>type_string</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>3</td>
        <td></td>
    </tr>
    
    <tr>
        <td>type_integer</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>4</td>
        <td></td>
    </tr>
    
    <tr>
        <td>type_real</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>5</td>
        <td></td>
    </tr>
    
    <tr>
        <td>type_logical</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>6</td>
        <td></td>
    </tr>
    
</table>

## fson_library
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>end_of_file*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>-1</td>
        <td></td>
    </tr>
    
    <tr>
        <td>end_of_record*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>-2</td>
        <td></td>
    </tr>
    
    <tr>
        <td>state_looking_for_value*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>1</td>
        <td></td>
    </tr>
    
    <tr>
        <td>state_in_object*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>2</td>
        <td></td>
    </tr>
    
    <tr>
        <td>state_in_pair_name*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>3</td>
        <td></td>
    </tr>
    
    <tr>
        <td>state_in_pair_value*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>4</td>
        <td></td>
    </tr>
    
    <tr>
        <td>pushed_index*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>pushed_char*</td>
        <td>Variable</td>
        <td>character</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## fwbs_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>bktlife</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Full power blanket lifetime (years)</p></td>
    </tr>
    
    <tr>
        <td>bktlife_cal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Calendar year blanket lifetime (years)</p></td>
    </tr>
    
    <tr>
        <td>coolmass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of water coolant (in shield, blanket, first wall, divertor) [kg]</p></td>
    </tr>
    
    <tr>
        <td>vvmass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vacuum vessel mass [kg]</p></td>
    </tr>
    
    <tr>
        <td>denstl</td>
        <td>Input</td>
        <td>real</td>
        <td>7800.0</td>
        <td><p>density of steel [kg m^-3]</p></td>
    </tr>
    
    <tr>
        <td>denwc</td>
        <td>Input</td>
        <td>real</td>
        <td>15630.0</td>
        <td><p>density of tungsten carbide [kg m^-3]</p></td>
    </tr>
    
    <tr>
        <td>dewmkg</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total mass of vacuum vessel + cryostat [kg] (calculated if blktmodel&gt;0)</p></td>
    </tr>
    
    <tr>
        <td>emult</td>
        <td>Input</td>
        <td>real</td>
        <td>1.269</td>
        <td><p>energy multiplication in blanket and shield</p></td>
    </tr>
    
    <tr>
        <td>emultmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>power due to energy multiplication in blanket and shield [MW]</p></td>
    </tr>
    
    <tr>
        <td>fblss</td>
        <td>Input</td>
        <td>real</td>
        <td>0.09705</td>
        <td><p>KIT blanket model: steel fraction of breeding zone</p></td>
    </tr>
    
    <tr>
        <td>fdiv</td>
        <td>Input</td>
        <td>real</td>
        <td>0.115</td>
        <td><p>Solid angle fraction taken by one divertor</p></td>
    </tr>
    
    <tr>
        <td>fhcd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>area fraction covered by heating/current drive apparatus plus diagnostics</p></td>
    </tr>
    
    <tr>
        <td>fhole</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>area fraction taken up by other holes (IFE)</p></td>
    </tr>
    
    <tr>
        <td>fwbsshape</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>switch for first wall, blanket, shield and vacuum vessel shape:</p>
<ul>
<li>=1 D-shaped (cylinder inboard + ellipse outboard)</li>
<li>=2 defined by two ellipses</li>
</ul></td>
    </tr>
    
    <tr>
        <td>life_fw_fpy</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>first wall full-power year lifetime (y)</p></td>
    </tr>
    
    <tr>
        <td>m_fw_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>first wall mass [kg]</p></td>
    </tr>
    
    <tr>
        <td>fw_armour_mass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>first wall armour mass [kg]</p></td>
    </tr>
    
    <tr>
        <td>fw_armour_thickness</td>
        <td>Input</td>
        <td>real</td>
        <td>0.005</td>
        <td><p>first wall armour thickness [m]</p></td>
    </tr>
    
    <tr>
        <td>fw_armour_vol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>first wall armour volume [m^3]</p></td>
    </tr>
    
    <tr>
        <td>i_blanket_type</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for blanket model:</p>
<ul>
<li>=1 CCFE HCPB model</li>
<li>=2 KIT HCPB model  # REMOVED, no longer usable</li>
<li>=3 CCFE HCPB model with Tritium Breeding Ratio calculation</li>
<li>=4 KIT HCLL model  # REMOVED, no longer usable</li>
<li>=5 DCLL model -  no nutronics model included (in development) please check/choose values for
                      'dual-coolant blanket' fractions (provided in this file).
                 -  please use primary_pumping = 0 or 1.</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_blkt_inboard</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for inboard blanket:</p>
<ul>
<li>=0 No inboard blanket (dr_blkt_inboard=0.0)</li>
<li>=1 Inboard blanket present</li>
</ul></td>
    </tr>
    
    <tr>
        <td>inuclear</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for nuclear heating in the coils:</p>
<ul>
<li>=0 Frances Fox model (default)</li>
<li>=1 Fixed by user (qnuc)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>qnuc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating in the coils (W) (<code>inuclear=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>li6enrich</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>lithium-6 enrichment of breeding material (%)</p></td>
    </tr>
    
    <tr>
        <td>pnucblkt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating in the blanket [MW]</p></td>
    </tr>
    
    <tr>
        <td>pnuc_cp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total nuclear heating in the ST centrepost [MW]</p></td>
    </tr>
    
    <tr>
        <td>pnuc_cp_sh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Neutronic shield nuclear heating in the ST centrepost [MW]</p></td>
    </tr>
    
    <tr>
        <td>pnuc_cp_tf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF neutronic nuclear heating in the ST centrepost [MW]</p></td>
    </tr>
    
    <tr>
        <td>pnucdiv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating in the divertor [MW]</p></td>
    </tr>
    
    <tr>
        <td>p_fw_nuclear_heat_total_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating in the first wall [MW]</p></td>
    </tr>
    
    <tr>
        <td>pnuchcd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating in the HCD apparatus and diagnostics [MW]</p></td>
    </tr>
    
    <tr>
        <td>pnucloss</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating lost via holes [MW]</p></td>
    </tr>
    
    <tr>
        <td>pnucvvplus</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating to vacuum vessel and beyond [MW]</p></td>
    </tr>
    
    <tr>
        <td>pnucshld</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating in the shield [MW]</p></td>
    </tr>
    
    <tr>
        <td>whtblkt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of blanket [kg]</p></td>
    </tr>
    
    <tr>
        <td>whtblss</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of blanket - steel part [kg]</p></td>
    </tr>
    
    <tr>
        <td>armour_fw_bl_mass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total mass of armour, first wall and blanket [kg]</p></td>
    </tr>
    
    <tr>
        <td>breeder_f</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>Volume ratio: Li4SiO4/(Be12Ti+Li4SiO4) (<code>iteration variable 108</code>)</p></td>
    </tr>
    
    <tr>
        <td>breeder_multiplier</td>
        <td>Input</td>
        <td>real</td>
        <td>0.75</td>
        <td><p>combined breeder/multipler fraction of blanket by volume</p></td>
    </tr>
    
    <tr>
        <td>vfcblkt</td>
        <td>Input</td>
        <td>real</td>
        <td>0.05295</td>
        <td><p>He coolant fraction of blanket by volume (<code>i_blanket_type= 1,3</code> (CCFE HCPB))</p></td>
    </tr>
    
    <tr>
        <td>vfpblkt</td>
        <td>Input</td>
        <td>real</td>
        <td>0.1</td>
        <td><p>He purge gas fraction of blanket by volume (<code>i_blanket_type= 1,3</code> (CCFE HCPB))</p></td>
    </tr>
    
    <tr>
        <td>whtblli4sio4</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of lithium orthosilicate in blanket [kg] (<code>i_blanket_type=1,3</code> (CCFE HCPB))</p></td>
    </tr>
    
    <tr>
        <td>whtbltibe12</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of titanium beryllide in blanket [kg] (<code>i_blanket_type=1,3</code> (CCFE HCPB))</p></td>
    </tr>
    
    <tr>
        <td>neut_flux_cp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Centrepost TF fast neutron flux (E &gt; 0.1 MeV) [m^(-2).^(-1)]
 This variable is only calculated for superconducting (i_tf_sup = 1 )
 spherical tokamal magnet designs (itart = 0)</p></td>
    </tr>
    
    <tr>
        <td>f_neut_shield</td>
        <td>Input</td>
        <td>real</td>
        <td>-1.0</td>
        <td><p>Fraction of nuclear power shielded before the CP magnet (ST)
 ( neut_absorb = -1 --&gt; a fit on simplified MCNP neutronic
 calculation is used assuming water cooled (13%) tungesten carbyde )</p></td>
    </tr>
    
    <tr>
        <td>f_a_fw_coolant_inboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard FW coolant cross-sectional area void fraction</p></td>
    </tr>
    
    <tr>
        <td>f_a_fw_coolant_outboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard/outboard FW coolant cross-sectional area void fraction</p></td>
    </tr>
    
    <tr>
        <td>psurffwi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Surface heat flux on first wall [MW] (sum = pradfw)</p></td>
    </tr>
    
    <tr>
        <td>psurffwo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Surface heat flux on first wall [MW] (sum = pradfw)</p></td>
    </tr>
    
    <tr>
        <td>vol_fw_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>First wall volume [m3]</p></td>
    </tr>
    
    <tr>
        <td>fblss_ccfe</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Fractions of blanket by volume: steel, lithium orthosilicate, titanium beryllide</p></td>
    </tr>
    
    <tr>
        <td>fblli2sio4</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Fractions of blanket by volume: steel, lithium orthosilicate, titanium beryllide</p></td>
    </tr>
    
    <tr>
        <td>fbltibe12</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Fractions of blanket by volume: steel, lithium orthosilicate, titanium beryllide</p></td>
    </tr>
    
    <tr>
        <td>breedmat</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>breeder material switch (i_blanket_type=2 (KIT HCPB)):</p>
<ul>
<li>=1 Lithium orthosilicate</li>
<li>=2 Lithium methatitanate</li>
<li>=3 Lithium zirconate</li>
</ul></td>
    </tr>
    
    <tr>
        <td>densbreed</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>density of breeder material [kg m^-3] (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>fblbe</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6</td>
        <td><p>beryllium fraction of blanket by volume (if <code>i_blanket_type=2</code>, is Be fraction of breeding zone)</p></td>
    </tr>
    
    <tr>
        <td>fblbreed</td>
        <td>Input</td>
        <td>real</td>
        <td>0.154</td>
        <td><p>breeder fraction of blanket breeding zone by volume (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>fblhebmi</td>
        <td>Input</td>
        <td>real</td>
        <td>0.4</td>
        <td><p>helium fraction of inboard blanket box manifold by volume (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>fblhebmo</td>
        <td>Input</td>
        <td>real</td>
        <td>0.4</td>
        <td><p>helium fraction of outboard blanket box manifold by volume (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>fblhebpi</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6595</td>
        <td><p>helium fraction of inboard blanket back plate by volume (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>fblhebpo</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6713</td>
        <td><p>helium fraction of outboard blanket back plate by volume (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>hcdportsize</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for size of heating/current drive ports (<code>i_blanket_type=2</code> (KIT HCPB)):</p>
<ul>
<li>=1 'small'</li>
<li>=2 'large'</li>
</ul></td>
    </tr>
    
    <tr>
        <td>nflutf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>peak fast neutron fluence on TF coil superconductor [n m^-2] (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>npdiv</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>number of divertor ports (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>nphcdin</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>number of inboard ports for heating/current drive (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>nphcdout</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>number of outboard ports for heating/current drive (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>tbr</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>tritium breeding ratio (<code>i_blanket_type=2,3</code> (KIT HCPB/HCLL))</p></td>
    </tr>
    
    <tr>
        <td>tritprate</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>tritium production rate [g day^-1] (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>wallpf</td>
        <td>Input</td>
        <td>real</td>
        <td>1.21</td>
        <td><p>neutron wall load peaking factor (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>whtblbreed</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of blanket - breeder part [kg] (<code>i_blanket_type=2</code> (KIT HCPB))</p></td>
    </tr>
    
    <tr>
        <td>whtblbe</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of blanket - beryllium part [kg]</p></td>
    </tr>
    
    <tr>
        <td>iblanket_thickness</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Blanket thickness switch (Do not set dr_blkt_inboard, dr_blkt_outboard, dr_fw_inboard or dr_fw_outboard when <code>i_blanket_type=3</code>):</p>
<ul>
<li>=1 thin    0.53 m inboard, 0.91 m outboard</li>
<li>=2 medium  0.64 m inboard, 1.11 m outboard</li>
<li>=3 thick   0.75 m inboard, 1.30 m outboard</li>
</ul></td>
    </tr>
    
    <tr>
        <td>primary_pumping</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Switch for pumping power for primary coolant (mechanical power only and peak first wall
 temperature is only calculated if <code>primary_pumping=2</code>):</p>
<ul>
<li>=0 User sets pump power directly (htpmw_blkt, htpmw_fw, htpmw_div, htpmw_shld)</li>
<li>=1 User sets pump power as a fraction of thermal power (fpumpblkt, fpumpfw, fpumpdiv, fpumpshld)</li>
<li>=2 Mechanical pumping power is calculated</li>
<li>=3 Mechanical pumping power is calculated using specified pressure drop</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_shield_mat</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for shield material - <em>currently only applied in costing routines</em> <code>cost_model = 2</code></p>
<ul>
<li>=0 Tungsten (default)</li>
<li>=1 Tungsten carbide</li>
</ul></td>
    </tr>
    
    <tr>
        <td>secondary_cycle</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for power conversion cycle:</p>
<ul>
<li>=0 Set efficiency for chosen blanket, from detailed models (divertor heat not used)</li>
<li>=1 Set efficiency for chosen blanket, from detailed models (divertor heat used)</li>
<li>=2 user input thermal-electric efficiency (etath)</li>
<li>=3 steam Rankine cycle</li>
<li>=4 supercritical CO2 cycle</li>
</ul></td>
    </tr>
    
    <tr>
        <td>secondary_cycle_liq</td>
        <td>Input</td>
        <td>integer</td>
        <td>4</td>
        <td><p>Switch for power conversion cycle for the liquid breeder component of the blanket:</p>
<ul>
<li>=2 user input thermal-electric efficiency (etath)</li>
<li>=4 supercritical CO2 cycle</li>
</ul></td>
    </tr>
    
    <tr>
        <td>coolwh</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for blanket coolant (set via blkttype):</p>
<ul>
<li>=1 helium</li>
<li>=2 pressurized water</li>
</ul></td>
    </tr>
    
    <tr>
        <td>afwi</td>
        <td>Input</td>
        <td>real</td>
        <td>0.008</td>
        <td><p>inner radius of inboard first wall/blanket coolant channels (stellarator only) [m]</p></td>
    </tr>
    
    <tr>
        <td>afwo</td>
        <td>Input</td>
        <td>real</td>
        <td>0.008</td>
        <td><p>inner radius of outboard first wall/blanket coolant channels (stellarator only) [m]</p></td>
    </tr>
    
    <tr>
        <td>i_fw_coolant_type</td>
        <td>Input</td>
        <td>character</td>
        <td>b'helium'</td>
        <td><p>switch for first wall coolant (can be different from blanket coolant):</p>
<ul>
<li>'helium'</li>
<li>'water'</li>
</ul></td>
    </tr>
    
    <tr>
        <td>dr_fw_wall</td>
        <td>Input</td>
        <td>real</td>
        <td>0.003</td>
        <td><p>wall thickness of first wall coolant channels [m]</p></td>
    </tr>
    
    <tr>
        <td>radius_fw_channel</td>
        <td>Input</td>
        <td>real</td>
        <td>0.006</td>
        <td><p>radius of first wall cooling channels [m]</p></td>
    </tr>
    
    <tr>
        <td>dx_fw_module</td>
        <td>Input</td>
        <td>real</td>
        <td>0.02</td>
        <td><p>Width of a FW module containing a cooling channel [m]</p></td>
    </tr>
    
    <tr>
        <td>temp_fw_coolant_in</td>
        <td>Input</td>
        <td>real</td>
        <td>573.0</td>
        <td><p>inlet temperature of first wall coolant [K]</p></td>
    </tr>
    
    <tr>
        <td>temp_fw_coolant_out</td>
        <td>Input</td>
        <td>real</td>
        <td>823.0</td>
        <td><p>outlet temperature of first wall coolant [K]</p></td>
    </tr>
    
    <tr>
        <td>pres_fw_coolant</td>
        <td>Input</td>
        <td>real</td>
        <td>15500000.0</td>
        <td><p>first wall coolant pressure [Pa] (<code>secondary_cycle&gt;1</code>)</p></td>
    </tr>
    
    <tr>
        <td>temp_fw_peak</td>
        <td>Input</td>
        <td>real</td>
        <td>873.0</td>
        <td><p>peak first wall temperature [K]</p></td>
    </tr>
    
    <tr>
        <td>roughness</td>
        <td>Input</td>
        <td>real</td>
        <td>1e-06</td>
        <td><p>first wall channel roughness epsilon [m]</p></td>
    </tr>
    
    <tr>
        <td>len_fw_channel</td>
        <td>Input</td>
        <td>real</td>
        <td>4.0</td>
        <td><p>Length of a single first wall channel (all in parallel) [m]
 (<code>iteration variable 114</code>, useful for <code>constraint equation 39</code>)</p></td>
    </tr>
    
    <tr>
        <td>f_fw_peak</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>peaking factor for first wall heat loads. (Applied separately to inboard and outboard loads.
 Applies to both neutron and surface loads. Only used to calculate peak temperature - not
 the coolant flow rate.)</p></td>
    </tr>
    
    <tr>
        <td>blpressure</td>
        <td>Input</td>
        <td>real</td>
        <td>15500000.0</td>
        <td><p>blanket coolant pressure [Pa] (<code>secondary_cycle&gt;1</code>)</p></td>
    </tr>
    
    <tr>
        <td>inlet_temp</td>
        <td>Input</td>
        <td>real</td>
        <td>573.0</td>
        <td><p>inlet temperature of blanket coolant  [K] (<code>secondary_cycle&gt;1</code>)</p></td>
    </tr>
    
    <tr>
        <td>outlet_temp</td>
        <td>Input</td>
        <td>real</td>
        <td>823.0</td>
        <td><p>Outlet temperature of blanket coolant [K] (<code>secondary_cycle&gt;1</code>)</p>
<ul>
<li>input if <code>coolwh=1</code> (helium)</li>
<li>calculated if <code>coolwh=2</code> (water)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>coolp</td>
        <td>Input</td>
        <td>real</td>
        <td>15500000.0</td>
        <td><p>blanket coolant pressure [Pa] (stellarator only)</p></td>
    </tr>
    
    <tr>
        <td>nblktmodpo</td>
        <td>Input</td>
        <td>integer</td>
        <td>8</td>
        <td><p>number of outboard blanket modules in poloidal direction (<code>secondary_cycle&gt;1</code>)</p></td>
    </tr>
    
    <tr>
        <td>nblktmodpi</td>
        <td>Input</td>
        <td>integer</td>
        <td>7</td>
        <td><p>number of inboard blanket modules in poloidal direction (<code>secondary_cycle&gt;1</code>)</p></td>
    </tr>
    
    <tr>
        <td>nblktmodto</td>
        <td>Input</td>
        <td>integer</td>
        <td>48</td>
        <td><p>number of outboard blanket modules in toroidal direction (<code>secondary_cycle&gt;1</code>)</p></td>
    </tr>
    
    <tr>
        <td>nblktmodti</td>
        <td>Input</td>
        <td>integer</td>
        <td>32</td>
        <td><p>number of inboard blanket modules in toroidal direction (<code>secondary_cycle&gt;1</code>)</p></td>
    </tr>
    
    <tr>
        <td>temp_fw_max</td>
        <td>Input</td>
        <td>real</td>
        <td>823.0</td>
        <td><p>maximum temperature of first wall material [K] (<code>secondary_cycle&gt;1</code>)</p></td>
    </tr>
    
    <tr>
        <td>fw_th_conductivity</td>
        <td>Input</td>
        <td>real</td>
        <td>28.34</td>
        <td><p>thermal conductivity of first wall material at 293 K (W/m/K) (Temperature dependence
 is as for unirradiated Eurofer)</p></td>
    </tr>
    
    <tr>
        <td>fvoldw</td>
        <td>Input</td>
        <td>real</td>
        <td>1.74</td>
        <td><p>area coverage factor for vacuum vessel volume</p></td>
    </tr>
    
    <tr>
        <td>fvolsi</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>area coverage factor for inboard shield volume</p></td>
    </tr>
    
    <tr>
        <td>fvolso</td>
        <td>Input</td>
        <td>real</td>
        <td>0.64</td>
        <td><p>area coverage factor for outboard shield volume</p></td>
    </tr>
    
    <tr>
        <td>fwclfr</td>
        <td>Input</td>
        <td>real</td>
        <td>0.15</td>
        <td><p>first wall coolant fraction (calculated if <code>i_pulsed_plant=1</code> or <code>ipowerflow=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>praddiv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radiation power incident on the divertor (MW)</p></td>
    </tr>
    
    <tr>
        <td>pradfw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radiation power incident on the first wall (MW)</p></td>
    </tr>
    
    <tr>
        <td>pradhcd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radiation power incident on the heating and current drive system (MW)</p></td>
    </tr>
    
    <tr>
        <td>pradloss</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radiation power lost through holes (eventually hits shield) (MW)
 Only used for stellarator</p></td>
    </tr>
    
    <tr>
        <td>ptfnuc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating in the TF coil (MW)</p></td>
    </tr>
    
    <tr>
        <td>ptfnucpm3</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>nuclear heating in the TF coil (MW/m3) (<code>blktmodel&gt;0</code>)</p></td>
    </tr>
    
    <tr>
        <td>r_cryostat_inboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>cryostat radius [m]</p></td>
    </tr>
    
    <tr>
        <td>z_cryostat_half_inside</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>cryostat height [m]</p></td>
    </tr>
    
    <tr>
        <td>dr_pf_cryostat</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>Radial distance between outer edge of furthest away PF coil (or stellarator
 modular coil) and cryostat [m]</p></td>
    </tr>
    
    <tr>
        <td>vol_cryostat</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Cryostat structure volume [m^3]</p></td>
    </tr>
    
    <tr>
        <td>vol_cryostat_internal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Internal volume of the cryostat [m^3]</p></td>
    </tr>
    
    <tr>
        <td>vdewin</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vacuum vessel volume [m^3]</p></td>
    </tr>
    
    <tr>
        <td>vfshld</td>
        <td>Input</td>
        <td>real</td>
        <td>0.25</td>
        <td><p>coolant void fraction in shield</p></td>
    </tr>
    
    <tr>
        <td>volblkt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of blanket [m^3]</p></td>
    </tr>
    
    <tr>
        <td>volblkti</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of inboard blanket [m^3]</p></td>
    </tr>
    
    <tr>
        <td>volblkto</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of outboard blanket [m^3]</p></td>
    </tr>
    
    <tr>
        <td>volshld</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of shield [m^3]</p></td>
    </tr>
    
    <tr>
        <td>whtshld</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of shield [kg]</p></td>
    </tr>
    
    <tr>
        <td>wpenshld</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of the penetration shield [kg]</p></td>
    </tr>
    
    <tr>
        <td>wtshldi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of inboard shield [kg]</p></td>
    </tr>
    
    <tr>
        <td>wtshldo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of outboard shield [kg]</p></td>
    </tr>
    
    <tr>
        <td>irefprop</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch to use REFPROP routines (stellarator only)</p></td>
    </tr>
    
    <tr>
        <td>fblli</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>lithium fraction of blanket by volume (stellarator only)</p></td>
    </tr>
    
    <tr>
        <td>fblli2o</td>
        <td>Input</td>
        <td>real</td>
        <td>0.08</td>
        <td><p>lithium oxide fraction of blanket by volume (stellarator only)</p></td>
    </tr>
    
    <tr>
        <td>fbllipb</td>
        <td>Input</td>
        <td>real</td>
        <td>0.68</td>
        <td><p>lithium lead fraction of blanket by volume (stellarator only)</p></td>
    </tr>
    
    <tr>
        <td>fblvd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vanadium fraction of blanket by volume (stellarator only)</p></td>
    </tr>
    
    <tr>
        <td>wtblli2o</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of blanket - Li_2O part [kg]</p></td>
    </tr>
    
    <tr>
        <td>wtbllipb</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of blanket - Li-Pb part [kg]</p></td>
    </tr>
    
    <tr>
        <td>whtblvd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of blanket - vanadium part [kg]</p></td>
    </tr>
    
    <tr>
        <td>whtblli</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of blanket - lithium part [kg]</p></td>
    </tr>
    
    <tr>
        <td>vfblkt</td>
        <td>Input</td>
        <td>real</td>
        <td>0.25</td>
        <td><p>coolant void fraction in blanket.</p></td>
    </tr>
    
    <tr>
        <td>blktmodel</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for blanket/tritium breeding model (see i_blanket_type):</p>
<ul>
<li>=0 original simple model</li>
<li>=1 KIT model based on a helium-cooled pebble-bed blanket (HCPB) reference design</li>
</ul></td>
    </tr>
    
    <tr>
        <td>declblkt</td>
        <td>Input</td>
        <td>real</td>
        <td>0.075</td>
        <td><p>neutron power deposition decay length of blanket structural material [m] (stellarators only)</p></td>
    </tr>
    
    <tr>
        <td>declfw</td>
        <td>Input</td>
        <td>real</td>
        <td>0.075</td>
        <td><p>neutron power deposition decay length of first wall structural material [m] (stellarators only)</p></td>
    </tr>
    
    <tr>
        <td>declshld</td>
        <td>Input</td>
        <td>real</td>
        <td>0.075</td>
        <td><p>neutron power deposition decay length of shield structural material [m] (stellarators only)</p></td>
    </tr>
    
    <tr>
        <td>blkttype</td>
        <td>Input</td>
        <td>integer</td>
        <td>3</td>
        <td><p>Switch for blanket type:</p>
<ul>
<li>=1 WCLL;</li>
<li>=2 HCLL; efficiency taken from M. Kovari 2016
 "PROCESS": A systems code for fusion power plants - Part 2: Engineering
 https://www.sciencedirect.com/science/article/pii/S0920379616300072
 Feedheat &amp; reheat cycle assumed</li>
<li>=3 HCPB; efficiency taken from M. Kovari 2016
 "PROCESS": A systems code for fusion power plants - Part 2: Engineering
 https://www.sciencedirect.com/science/article/pii/S0920379616300072
 Feedheat &amp; reheat cycle assumed</li>
</ul></td>
    </tr>
    
    <tr>
        <td>etaiso</td>
        <td>Input</td>
        <td>real</td>
        <td>0.85</td>
        <td><p>isentropic efficiency of FW and blanket coolant pumps</p></td>
    </tr>
    
    <tr>
        <td>etahtp</td>
        <td>Input</td>
        <td>real</td>
        <td>0.95</td>
        <td><p>electrical efficiency of primary coolant pumps</p>
<hr>
<p>BLANKET REFACTOR
 For DCLL, but to be used by all mods that share blanket library after testing.
 Thermodynamic Model for primary_pumping == 2</p>
<hr></td>
    </tr>
    
    <tr>
        <td>ipump</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for whether the FW and BB are on the same pump system
 i.e. do they have the same primary coolant or not
  - =0    FW and BB have the same primary coolant, flow = FWin-&gt;FWout-&gt;BBin-&gt;BBout
  - =1    FW and BB have the different primary coolant and are on different pump systems</p></td>
    </tr>
    
    <tr>
        <td>i_bb_liq</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for Liquid Metal Breeder Material
  - =0   PbLi
  - =1   Li</p></td>
    </tr>
    
    <tr>
        <td>icooldual</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch to specify whether breeding blanket is single-cooled or dual-coolant.
  - =0    Single coolant used for FW and Blanket (H2O or He). Solid Breeder.
  - =1    Single coolant used for FW and Blanket (H2O or He). Liquid metal breeder
          circulted for tritium extraction.
  - =2    Dual coolant: primary coolant (H2O or He) for FW and blanket structure;
          secondary coolant is self-cooled liquid metal breeder.</p></td>
    </tr>
    
    <tr>
        <td>ifci</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for Flow Channel Insert (FCI) type if liquid metal breeder blanket.
  - =0    Thin conducting walls, default electrical conductivity (bz_channel_conduct_liq) is Eurofer
  - =1    Insulating Material, assumed perfect electrical insulator, default density (den_ceramic) is for SiC
  - =2    Insulating Material, electrical conductivity (bz_channel_conduct_liq) is input (default Eurofer), default density (den_ceramic) is for SiC</p></td>
    </tr>
    
    <tr>
        <td>ims</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for Multi Module Segment (MMS) or Single Modle Segment (SMS)
  - =0    MMS
  - =1    SMS</p></td>
    </tr>
    
    <tr>
        <td>n_liq_recirc</td>
        <td>Input</td>
        <td>integer</td>
        <td>10</td>
        <td><p>Number of liquid metal breeder recirculations per day, for use with icooldual=1</p></td>
    </tr>
    
    <tr>
        <td>r_f_liq_ib</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>Radial fraction of BZ liquid channels</p></td>
    </tr>
    
    <tr>
        <td>r_f_liq_ob</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>Radial fraction of BZ liquid channels</p></td>
    </tr>
    
    <tr>
        <td>w_f_liq_ib</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>Toroidal fraction of BZ liquid channels</p></td>
    </tr>
    
    <tr>
        <td>w_f_liq_ob</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>Toroidal fraction of BZ liquid channels</p></td>
    </tr>
    
    <tr>
        <td>den_ceramic</td>
        <td>Input</td>
        <td>real</td>
        <td>3210.0</td>
        <td><p>FCI material density</p></td>
    </tr>
    
    <tr>
        <td>th_wall_secondary</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0125</td>
        <td><p>Liquid metal coolant/breeder wall thickness thin conductor or FCI [m]</p></td>
    </tr>
    
    <tr>
        <td>bz_channel_conduct_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>833000.0</td>
        <td><p>Liquid metal coolant/breeder thin conductor or FCI wall conductance [A V^-1 m^-1]</p></td>
    </tr>
    
    <tr>
        <td>a_bz_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>0.2</td>
        <td><p>Toroidal width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone</p></td>
    </tr>
    
    <tr>
        <td>b_bz_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>0.2</td>
        <td><p>Radial width of the rectangular cooling channel [m] for long poloidal sections of blanket breeding zone</p></td>
    </tr>
    
    <tr>
        <td>nopol</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Number of poloidal sections in a liquid metal breeder/coolant channel for module/segment</p></td>
    </tr>
    
    <tr>
        <td>nopipes</td>
        <td>Input</td>
        <td>integer</td>
        <td>4</td>
        <td><p>Number of Liquid metal breeder/coolant channels per module/segment</p></td>
    </tr>
    
    <tr>
        <td>den_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>9500.0</td>
        <td><p>Liquid metal breeder/coolant density [kg m^-3]</p></td>
    </tr>
    
    <tr>
        <td>wht_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Liquid metal</p></td>
    </tr>
    
    <tr>
        <td>wht_liq_ib</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Liquid metal</p></td>
    </tr>
    
    <tr>
        <td>wht_liq_ob</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Liquid metal</p></td>
    </tr>
    
    <tr>
        <td>specific_heat_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>190.0</td>
        <td><p>Liquid metal breeder/coolant specific heat [J kg^-1 K^-1]</p></td>
    </tr>
    
    <tr>
        <td>thermal_conductivity_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>Liquid metal breeder/coolant thermal conductivity [W m^-1 K^-1]</p></td>
    </tr>
    
    <tr>
        <td>dynamic_viscosity_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Liquid metal breeder/coolant dynamic viscosity [Pa s]</p></td>
    </tr>
    
    <tr>
        <td>electrical_conductivity_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Liquid metal breeder/coolant electrical conductivity [Ohm m]</p></td>
    </tr>
    
    <tr>
        <td>hartmann_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Hartmann number</p></td>
    </tr>
    
    <tr>
        <td>b_mag_blkt</td>
        <td>Input</td>
        <td>real</td>
        <td>[5. 5.]</td>
        <td><p>Toroidal Magnetic feild strength for IB/OB blanket [T]</p></td>
    </tr>
    
    <tr>
        <td>etaiso_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>0.85</td>
        <td><p>Isentropic efficiency of blanket liquid breeder/coolant pumps</p></td>
    </tr>
    
    <tr>
        <td>blpressure_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>1700000.0</td>
        <td><p>blanket liquid metal breeder/coolant pressure [Pa]</p></td>
    </tr>
    
    <tr>
        <td>inlet_temp_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>570.0</td>
        <td><p>Inlet (scan var 68) and Outlet (scan var 69) temperature of the liquid breeder/coolant [K]</p></td>
    </tr>
    
    <tr>
        <td>outlet_temp_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>720.0</td>
        <td><p>Inlet (scan var 68) and Outlet (scan var 69) temperature of the liquid breeder/coolant [K]</p></td>
    </tr>
    
    <tr>
        <td>den_fw_coolant</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Density of the FW primary coolant</p></td>
    </tr>
    
    <tr>
        <td>visc_fw_coolant</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Viscosity of the FW primary coolant</p></td>
    </tr>
    
    <tr>
        <td>rhof_bl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Density of the blanket primary coolant</p></td>
    </tr>
    
    <tr>
        <td>visc_bl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Viscosity of the blanket primary coolant</p></td>
    </tr>
    
    <tr>
        <td>cp_fw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Spesific heat for FW and blanket primary coolant(s)</p></td>
    </tr>
    
    <tr>
        <td>cv_fw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Spesific heat for FW and blanket primary coolant(s)</p></td>
    </tr>
    
    <tr>
        <td>cp_bl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Spesific heat for FW and blanket primary coolant(s)</p></td>
    </tr>
    
    <tr>
        <td>cv_bl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Spesific heat for FW and blanket primary coolant(s)</p></td>
    </tr>
    
    <tr>
        <td>f_nuc_pow_bz_struct</td>
        <td>Input</td>
        <td>real</td>
        <td>0.34</td>
        <td><p>For a dual-coolant blanket, fraction of BZ power cooled by primary coolant</p></td>
    </tr>
    
    <tr>
        <td>f_nuc_pow_bz_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>0.66</td>
        <td><p>For a dual-coolant blanket, fraction of BZ self-cooled power (secondary coolant)</p></td>
    </tr>
    
    <tr>
        <td>pnuc_fw_ratio_dcll</td>
        <td>Input</td>
        <td>real</td>
        <td>0.14</td>
        <td><p>For a dual-coolant blanket, ratio of FW/Blanket nuclear power as fraction of total</p></td>
    </tr>
    
    <tr>
        <td>pnuc_blkt_ratio_dcll</td>
        <td>Input</td>
        <td>real</td>
        <td>0.86</td>
        <td><p>For a dual-coolant blanket, ratio of FW/Blanket nuclear power as fraction of total</p></td>
    </tr>
    
    <tr>
        <td>bzfllengi_n_rad</td>
        <td>Input</td>
        <td>integer</td>
        <td>4</td>
        <td><p>Number of radial and poloidal sections that make up the total primary coolant flow
 length in a blanket module (IB and OB)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengi_n_pol</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Number of radial and poloidal sections that make up the total primary coolant flow
 length in a blanket module (IB and OB)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengo_n_rad</td>
        <td>Input</td>
        <td>integer</td>
        <td>4</td>
        <td><p>Number of radial and poloidal sections that make up the total primary coolant flow
 length in a blanket module (IB and OB)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengo_n_pol</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Number of radial and poloidal sections that make up the total primary coolant flow
 length in a blanket module (IB and OB)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengi_n_rad_liq</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Number of radial and poloidal sections that make up the total secondary coolant/breeder
 flow length in a blanket module (IB and OB)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengi_n_pol_liq</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Number of radial and poloidal sections that make up the total secondary coolant/breeder
 flow length in a blanket module (IB and OB)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengo_n_rad_liq</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Number of radial and poloidal sections that make up the total secondary coolant/breeder
 flow length in a blanket module (IB and OB)</p></td>
    </tr>
    
    <tr>
        <td>bzfllengo_n_pol_liq</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Number of radial and poloidal sections that make up the total secondary coolant/breeder
 flow length in a blanket module (IB and OB)</p></td>
    </tr>
    
</table>

## global_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>icase</td>
        <td>Input</td>
        <td>character</td>
        <td>b'Steady-state tokamak model                      '</td>
        <td><p>power plant type</p></td>
    </tr>
    
    <tr>
        <td>runtitle</td>
        <td>Input</td>
        <td>character</td>
        <td>b"Run Title (change this line using input variable 'runtitle')                                                                                                                        "</td>
        <td><p>short descriptive title for the run</p></td>
    </tr>
    
    <tr>
        <td>verbose</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for turning on/off diagnostic messages</p>
<ul>
<li>=0 turn off diagnostics</li>
<li>=1 turn on diagnostics</li>
</ul></td>
    </tr>
    
    <tr>
        <td>run_tests</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>turns on built-in tests if set to 1</p></td>
    </tr>
    
    <tr>
        <td>maxcal</td>
        <td>Input</td>
        <td>integer</td>
        <td>200</td>
        <td><p>maximum number of VMCON iterations</p></td>
    </tr>
    
    <tr>
        <td>fileprefix</td>
        <td>Input</td>
        <td>character</td>
        <td>b'                                                                                                                                                                                                                                                                                                                                                                                                                '</td>
        <td><p>input file prefix</p></td>
    </tr>
    
    <tr>
        <td>output_prefix</td>
        <td>Input</td>
        <td>character</td>
        <td>b'                                                                                                                                                                                                                                                                                                                                                                                                                '</td>
        <td><p>output file prefix</p></td>
    </tr>
    
    <tr>
        <td>xlabel</td>
        <td>Input</td>
        <td>character</td>
        <td>b'                         '</td>
        <td><p>scan parameter description label</p></td>
    </tr>
    
    <tr>
        <td>vlabel</td>
        <td>Input</td>
        <td>character</td>
        <td>b'                         '</td>
        <td><p>scan value name label</p></td>
    </tr>
    
    <tr>
        <td>xlabel_2</td>
        <td>Input</td>
        <td>character</td>
        <td>b'                         '</td>
        <td><p>scan parameter description label (2nd dimension)</p></td>
    </tr>
    
    <tr>
        <td>vlabel_2</td>
        <td>Input</td>
        <td>character</td>
        <td>b'                         '</td>
        <td><p>scan value name label (2nd dimension)</p></td>
    </tr>
    
    <tr>
        <td>iscan_global</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Makes iscan available globally.</p></td>
    </tr>
    
    <tr>
        <td>convergence_parameter</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>VMCON convergence parameter "sum"</p></td>
    </tr>
    
</table>

## ccfe_hcpb_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>ip</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>ofile</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>armour_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>FW armour density [kg/m3]</p></td>
    </tr>
    
    <tr>
        <td>fw_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>FW density [kg/m3]</p></td>
    </tr>
    
    <tr>
        <td>blanket_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Blanket density [kg/m3]</p></td>
    </tr>
    
    <tr>
        <td>shield_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Shield density [kg/m3]</p></td>
    </tr>
    
    <tr>
        <td>vv_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Vacuum vessel density [kg/m3]</p></td>
    </tr>
    
    <tr>
        <td>x_blanket</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Blanket exponent (tonne/m2)</p></td>
    </tr>
    
    <tr>
        <td>x_shield</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Shield exponent (tonne/m2)</p></td>
    </tr>
    
    <tr>
        <td>tfc_nuc_heating</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Unit nuclear heating in TF coil (W per W of fusion power)</p></td>
    </tr>
    
    <tr>
        <td>fw_armour_u_nuc_heating</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Unit heating of FW and armour in FW armour (W/kg per W of fusion power)</p></td>
    </tr>
    
    <tr>
        <td>shld_u_nuc_heating</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Unit nuclear heating in shield (W per W of fusion power)</p></td>
    </tr>
    
    <tr>
        <td>pnuc_tot_blk_sector</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total nuclear power deposited in blanket covered sector (FW, BLKT, SHLD, TF) (MW)</p></td>
    </tr>
    
    <tr>
        <td>exp_blanket</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Exponential factors in nuclear heating calcs</p></td>
    </tr>
    
    <tr>
        <td>exp_shield1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Exponential factors in nuclear heating calcs</p></td>
    </tr>
    
    <tr>
        <td>exp_shield2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Exponential factors in nuclear heating calcs</p></td>
    </tr>
    
</table>

## heat_transport_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>baseel</td>
        <td>Input</td>
        <td>real</td>
        <td>5000000.0</td>
        <td><p>base plant electric load (W)</p></td>
    </tr>
    
    <tr>
        <td>crypmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>cryogenic plant power (MW)</p></td>
    </tr>
    
    <tr>
        <td>crypmw_max</td>
        <td>Input</td>
        <td>real</td>
        <td>50.0</td>
        <td><p>Maximum cryogenic plant power (MW)
 Constraint equation icc = 87
 Scan variable nwseep = 56</p></td>
    </tr>
    
    <tr>
        <td>f_crypmw</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum cryogenic plant power
 Iteration variable ixc = 164
 Constraint equation icc = 87</p></td>
    </tr>
    
    <tr>
        <td>etatf</td>
        <td>Input</td>
        <td>real</td>
        <td>0.9</td>
        <td><p>AC to resistive power conversion for TF coils</p></td>
    </tr>
    
    <tr>
        <td>etath</td>
        <td>Input</td>
        <td>real</td>
        <td>0.35</td>
        <td><p>thermal to electric conversion efficiency if <code>secondary_cycle=2</code>; otherwise calculated.</p></td>
    </tr>
    
    <tr>
        <td>etath_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>0.35</td>
        <td></td>
    </tr>
    
    <tr>
        <td>fachtmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>facility heat removal (MW)</p></td>
    </tr>
    
    <tr>
        <td>fcsht</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total baseline power required at all times (MW)</p></td>
    </tr>
    
    <tr>
        <td>fgrosbop</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>scaled fraction of gross power to balance-of-plant</p></td>
    </tr>
    
    <tr>
        <td>fmgdmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>power to mgf (motor-generator flywheel) units (MW) (ignored if <code>iscenr=2</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpumpblkt</td>
        <td>Input</td>
        <td>real</td>
        <td>0.005</td>
        <td><p>fraction of total blanket thermal power required to drive the blanket
 coolant pumps (default assumes water coolant) (<code>secondary_cycle=0</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpumpdiv</td>
        <td>Input</td>
        <td>real</td>
        <td>0.005</td>
        <td><p>fraction of total divertor thermal power required to drive the divertor
 coolant pumps (default assumes water coolant)</p></td>
    </tr>
    
    <tr>
        <td>fpumpfw</td>
        <td>Input</td>
        <td>real</td>
        <td>0.005</td>
        <td><p>fraction of total first wall thermal power required to drive the FW coolant
 pumps (default assumes water coolant) (<code>secondary_cycle=0</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpumpshld</td>
        <td>Input</td>
        <td>real</td>
        <td>0.005</td>
        <td><p>fraction of total shield thermal power required to drive the shield coolant
 pumps (default assumes water coolant)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_min</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Minimum total electrical power for primary coolant pumps (MW) (NOT RECOMMENDED)</p></td>
    </tr>
    
    <tr>
        <td>helpow</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Heat removal at cryogenic temperature tmpcry (W)</p></td>
    </tr>
    
    <tr>
        <td>helpow_cryal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Heat removal at cryogenic temperature tcoolin (W)</p></td>
    </tr>
    
    <tr>
        <td>htpmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>heat transport system electrical pump power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_blkt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>blanket primary coolant mechanical pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_blkt_liq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>blanket secondary coolant mechanical pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_blkt_tot</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>blanket primary + secondary coolant mechanical pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_div</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>divertor coolant mechanical pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_fw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>first wall coolant mechanical pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_shld</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>shield and vacuum vessel coolant mechanical pumping power (MW)</p></td>
    </tr>
    
    <tr>
        <td>htpsecmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Waste power lost from primary coolant pumps (MW)</p></td>
    </tr>
    
    <tr>
        <td>ipowerflow</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for power flow model:</p>
<ul>
<li>=0 pre-2014 version</li>
<li>=1 comprehensive 2014 model</li>
</ul></td>
    </tr>
    
    <tr>
        <td>iprimshld</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for shield thermal power destiny:</p>
<ul>
<li>=0 does not contribute to energy generation cycle</li>
<li>=1 contributes to energy generation cycle</li>
</ul></td>
    </tr>
    
    <tr>
        <td>nphx</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>number of primary heat exchangers</p></td>
    </tr>
    
    <tr>
        <td>pacpmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total pulsed power system load (MW)</p></td>
    </tr>
    
    <tr>
        <td>peakmva</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>peak MVA requirement</p></td>
    </tr>
    
    <tr>
        <td>pfwdiv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>heat removal from first wall/divertor (MW)</p></td>
    </tr>
    
    <tr>
        <td>pgrossmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>gross electric power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pinjht</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>power dissipated in heating and current drive system (MW)</p></td>
    </tr>
    
    <tr>
        <td>pinjmax</td>
        <td>Input</td>
        <td>real</td>
        <td>120.0</td>
        <td><p>maximum injector power during pulse (heating and ramp-up/down phase) (MW)</p></td>
    </tr>
    
    <tr>
        <td>pinjwp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>injector wall plug power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pinjwpfix</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>secondary injector wall plug power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pnetelmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>net electric power (MW)</p></td>
    </tr>
    
    <tr>
        <td>precircmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>recirculating electric power (MW)</p></td>
    </tr>
    
    <tr>
        <td>priheat</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total thermal power removed from fusion core (MW)</p></td>
    </tr>
    
    <tr>
        <td>psecdiv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Low-grade heat lost in divertor (MW)</p></td>
    </tr>
    
    <tr>
        <td>psechcd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Low-grade heat lost into HCD apparatus (MW)</p></td>
    </tr>
    
    <tr>
        <td>psechtmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Low-grade heat (MW)</p></td>
    </tr>
    
    <tr>
        <td>pseclossmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Low-grade heat (VV + lost)(MW)</p></td>
    </tr>
    
    <tr>
        <td>psecshld</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Low-grade heat deposited in shield (MW)</p></td>
    </tr>
    
    <tr>
        <td>pthermmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>High-grade heat useful for electric production (MW)</p></td>
    </tr>
    
    <tr>
        <td>pwpm2</td>
        <td>Input</td>
        <td>real</td>
        <td>150.0</td>
        <td><p>base AC power requirement per unit floor area (W/m2)</p></td>
    </tr>
    
    <tr>
        <td>tfacpd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total steady state TF coil AC power demand (MW)</p></td>
    </tr>
    
    <tr>
        <td>tlvpmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>estimate of total low voltage power (MW)</p></td>
    </tr>
    
    <tr>
        <td>trithtmw</td>
        <td>Input</td>
        <td>real</td>
        <td>15.0</td>
        <td><p>power required for tritium processing (MW)</p></td>
    </tr>
    
    <tr>
        <td>tturb</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>coolant temperature at turbine inlet (K) (<code>secondary_cycle = 3,4</code>)</p></td>
    </tr>
    
    <tr>
        <td>vachtmw</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>vacuum pump power (MW)</p></td>
    </tr>
    
</table>

## ife_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>maxmat</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>8</td>
        <td><p>Total number of materials in IFE device. Material numbers are as follows:</p>
<ul>
<li>=0 void</li>
<li>=1 steel</li>
<li>=2 carbon cloth</li>
<li>=3 FLiBe</li>
<li>=4 lithium oxide Li2O</li>
<li>=5 concrete</li>
<li>=6 helium</li>
<li>=7 xenon</li>
<li>=8 lithium</li>
</ul></td>
    </tr>
    
    <tr>
        <td>bldr</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>radial thickness of IFE blanket (m; calculated <code>if ifetyp=4</code>)</p></td>
    </tr>
    
    <tr>
        <td>bldrc</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>radial thickness of IFE curtain (m; <code>ifetyp=4</code>)</p></td>
    </tr>
    
    <tr>
        <td>bldzl</td>
        <td>Input</td>
        <td>real</td>
        <td>4.0</td>
        <td><p>vertical thickness of IFE blanket below chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>bldzu</td>
        <td>Input</td>
        <td>real</td>
        <td>4.0</td>
        <td><p>vertical thickness of IFE blanket above chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>blmatf</td>
        <td>Input</td>
        <td>real</td>
        <td>[[0.05 0.   0.45 0.   0.2  0.   0.3  0.   0.  ]
 [0.05 0.   0.45 0.   0.2  0.   0.3  0.   0.  ]
 [0.05 0.   0.45 0.   0.2  0.   0.3  0.   0.  ]]</td>
        <td><p>IFE blanket material fractions</p></td>
    </tr>
    
    <tr>
        <td>blmatm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE blanket material masses (kg)</p></td>
    </tr>
    
    <tr>
        <td>blmatv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE blanket material volumes (m3)</p></td>
    </tr>
    
    <tr>
        <td>blvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE blanket volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>cdriv0</td>
        <td>Input</td>
        <td>real</td>
        <td>154.3</td>
        <td><p>IFE generic/laser driver cost at edrive=0 (M$)</p></td>
    </tr>
    
    <tr>
        <td>cdriv1</td>
        <td>Input</td>
        <td>real</td>
        <td>163.2</td>
        <td><p>IFE low energy heavy ion beam driver cost extrapolated to <code>edrive=0</code> (M$)</p></td>
    </tr>
    
    <tr>
        <td>cdriv2</td>
        <td>Input</td>
        <td>real</td>
        <td>244.9</td>
        <td><p>IFE high energy heavy ion beam driver cost extrapolated to <code>edrive=0</code> (M$)</p></td>
    </tr>
    
    <tr>
        <td>cdriv3</td>
        <td>Input</td>
        <td>real</td>
        <td>1.463</td>
        <td><p>IFE driver cost ($/J wall plug) (<code>ifedrv==3</code>)</p></td>
    </tr>
    
    <tr>
        <td>chdzl</td>
        <td>Input</td>
        <td>real</td>
        <td>9.0</td>
        <td><p>vertical thickness of IFE chamber below centre (m)</p></td>
    </tr>
    
    <tr>
        <td>chdzu</td>
        <td>Input</td>
        <td>real</td>
        <td>9.0</td>
        <td><p>vertical thickness of IFE chamber above centre (m)</p></td>
    </tr>
    
    <tr>
        <td>chmatf</td>
        <td>Input</td>
        <td>real</td>
        <td>[1. 0. 0. 0. 0. 0. 0. 0. 0.]</td>
        <td><p>IFE chamber material fractions</p></td>
    </tr>
    
    <tr>
        <td>chmatm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE chamber material masses (kg)</p></td>
    </tr>
    
    <tr>
        <td>chmatv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE chamber material volumes (m3)</p></td>
    </tr>
    
    <tr>
        <td>chrad</td>
        <td>Input</td>
        <td>real</td>
        <td>6.5</td>
        <td><p>radius of IFE chamber (m) (<code>iteration variable 84</code>)</p></td>
    </tr>
    
    <tr>
        <td>chvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE chamber volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>dcdrv0</td>
        <td>Input</td>
        <td>real</td>
        <td>111.4</td>
        <td><p>IFE generic/laser driver cost gradient (M$/MJ)</p></td>
    </tr>
    
    <tr>
        <td>dcdrv1</td>
        <td>Input</td>
        <td>real</td>
        <td>78.0</td>
        <td><p>HIB driver cost gradient at low energy (M$/MJ)</p></td>
    </tr>
    
    <tr>
        <td>dcdrv2</td>
        <td>Input</td>
        <td>real</td>
        <td>59.9</td>
        <td><p>HIB driver cost gradient at high energy (M$/MJ)</p></td>
    </tr>
    
    <tr>
        <td>drveff</td>
        <td>Input</td>
        <td>real</td>
        <td>0.28</td>
        <td><p>IFE driver wall plug to target efficiency (<code>ifedrv=0,3</code>) (<code>iteration variable 82</code>)</p></td>
    </tr>
    
    <tr>
        <td>edrive</td>
        <td>Input</td>
        <td>real</td>
        <td>5000000.0</td>
        <td><p>IFE driver energy (J) (<code>iteration variable 81</code>)</p></td>
    </tr>
    
    <tr>
        <td>etadrv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE driver wall plug to target efficiency</p></td>
    </tr>
    
    <tr>
        <td>etali</td>
        <td>Input</td>
        <td>real</td>
        <td>0.4</td>
        <td><p>IFE lithium pump wall plug efficiency (<code>ifetyp=4</code>)</p></td>
    </tr>
    
    <tr>
        <td>etave</td>
        <td>Input</td>
        <td>real</td>
        <td>[0.082 0.079 0.076 0.073 0.069 0.066 0.062 0.059 0.055 0.051]</td>
        <td><p>IFE driver efficiency vs driver energy (<code>ifedrv=-1</code>)</p></td>
    </tr>
    
    <tr>
        <td>fauxbop</td>
        <td>Input</td>
        <td>real</td>
        <td>0.06</td>
        <td><p>fraction of gross electric power to balance-of-plant (IFE)</p></td>
    </tr>
    
    <tr>
        <td>fbreed</td>
        <td>Input</td>
        <td>real</td>
        <td>0.51</td>
        <td><p>fraction of breeder external to device core</p></td>
    </tr>
    
    <tr>
        <td>fburn</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3333</td>
        <td><p>IFE burn fraction (fraction of tritium fused/target)</p></td>
    </tr>
    
    <tr>
        <td>flirad</td>
        <td>Input</td>
        <td>real</td>
        <td>0.78</td>
        <td><p>radius of FLiBe/lithium inlet (m) (<code>ifetyp=3,4</code>)</p></td>
    </tr>
    
    <tr>
        <td>frrmax</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for maximum IFE repetition rate (<code>constraint equation 50</code>, <code>iteration variable 86</code>)</p></td>
    </tr>
    
    <tr>
        <td>fwdr</td>
        <td>Input</td>
        <td>real</td>
        <td>0.01</td>
        <td><p>radial thickness of IFE first wall (m)</p></td>
    </tr>
    
    <tr>
        <td>fwdzl</td>
        <td>Input</td>
        <td>real</td>
        <td>0.01</td>
        <td><p>vertical thickness of IFE first wall below chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>fwdzu</td>
        <td>Input</td>
        <td>real</td>
        <td>0.01</td>
        <td><p>vertical thickness of IFE first wall above chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>fwmatf</td>
        <td>Input</td>
        <td>real</td>
        <td>[[0.05 0.   0.95 0.   0.   0.   0.   0.   0.  ]
 [0.05 0.   0.95 0.   0.   0.   0.   0.   0.  ]
 [0.05 0.   0.95 0.   0.   0.   0.   0.   0.  ]]</td>
        <td><p>IFE first wall material fractions</p></td>
    </tr>
    
    <tr>
        <td>fwmatm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE first wall material masses (kg)</p></td>
    </tr>
    
    <tr>
        <td>fwmatv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE first wall material volumes (kg)</p></td>
    </tr>
    
    <tr>
        <td>fwvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE first wall volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>gain</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE target gain</p></td>
    </tr>
    
    <tr>
        <td>gainve</td>
        <td>Input</td>
        <td>real</td>
        <td>[ 60.  95. 115. 125. 133. 141. 152. 160. 165. 170.]</td>
        <td><p>IFE target gain vs driver energy (<code>ifedrv=-1</code>)</p></td>
    </tr>
    
    <tr>
        <td>htpmw_ife</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE heat transport system electrical pump power (MW)</p></td>
    </tr>
    
    <tr>
        <td>ife</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for IFE option:</p>
<ul>
<li>=0 use tokamak, RFP or stellarator model</li>
<li>=1 use IFE model</li>
</ul></td>
    </tr>
    
    <tr>
        <td>ifedrv</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Switch for type of IFE driver:</p>
<ul>
<li>=-1 use gainve, etave for gain and driver efficiency</li>
<li>=0 use tgain, drveff for gain and driver efficiency</li>
<li>=1 use laser driver based on SOMBRERO design</li>
<li>=2 use heavy ion beam driver based on OSIRIS</li>
<li>=3 Input pfusife, rrin and drveff</li>
</ul></td>
    </tr>
    
    <tr>
        <td>ifetyp</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for type of IFE device build:</p>
<ul>
<li>=0 generic (cylindrical) build</li>
<li>=1 OSIRIS-like build</li>
<li>=2 SOMBRERO-like build</li>
<li>=3 HYLIFE-II-like build</li>
<li>=4 2019 build</li>
</ul></td>
    </tr>
    
    <tr>
        <td>lipmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE lithium pump power (MW; <code>ifetyp=4</code>)</p></td>
    </tr>
    
    <tr>
        <td>mcdriv</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>IFE driver cost multiplier</p></td>
    </tr>
    
    <tr>
        <td>mflibe</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total mass of FLiBe (kg)</p></td>
    </tr>
    
    <tr>
        <td>pdrive</td>
        <td>Input</td>
        <td>real</td>
        <td>23000000.0</td>
        <td><p>IFE driver power reaching target (W) (<code>iteration variable 85</code>)</p></td>
    </tr>
    
    <tr>
        <td>pfusife</td>
        <td>Input</td>
        <td>real</td>
        <td>1000.0</td>
        <td><p>IFE input fusion power (MW) (<code>ifedrv=3 only</code>; <code>itv 155</code>)</p></td>
    </tr>
    
    <tr>
        <td>pifecr</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>IFE cryogenic power requirements (MW)</p></td>
    </tr>
    
    <tr>
        <td>ptargf</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>IFE target factory power at 6 Hz repetition rate (MW)</p></td>
    </tr>
    
    <tr>
        <td>r1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE device radial build (m)</p></td>
    </tr>
    
    <tr>
        <td>r2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE device radial build (m)</p></td>
    </tr>
    
    <tr>
        <td>r3</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE device radial build (m)</p></td>
    </tr>
    
    <tr>
        <td>r4</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE device radial build (m)</p></td>
    </tr>
    
    <tr>
        <td>r5</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE device radial build (m)</p></td>
    </tr>
    
    <tr>
        <td>r6</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE device radial build (m)</p></td>
    </tr>
    
    <tr>
        <td>r7</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE device radial build (m)</p></td>
    </tr>
    
    <tr>
        <td>reprat</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE driver repetition rate (Hz)</p></td>
    </tr>
    
    <tr>
        <td>rrin</td>
        <td>Input</td>
        <td>real</td>
        <td>6.0</td>
        <td><p>Input IFE repetition rate (Hz) (<code>ifedrv=3 only</code>; <code>itv 156</code>)</p></td>
    </tr>
    
    <tr>
        <td>rrmax</td>
        <td>Input</td>
        <td>real</td>
        <td>20.0</td>
        <td><p>maximum IFE repetition rate (Hz)</p></td>
    </tr>
    
    <tr>
        <td>shdr</td>
        <td>Input</td>
        <td>real</td>
        <td>1.7</td>
        <td><p>radial thickness of IFE shield (m)</p></td>
    </tr>
    
    <tr>
        <td>shdzl</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>vertical thickness of IFE shield below chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>shdzu</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>vertical thickness of IFE shield above chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>shmatf</td>
        <td>Input</td>
        <td>real</td>
        <td>[[0.05  0.19  0.    0.    0.    0.665 0.095 0.    0.   ]
 [0.05  0.19  0.    0.    0.    0.665 0.095 0.    0.   ]
 [0.05  0.19  0.    0.    0.    0.665 0.095 0.    0.   ]]</td>
        <td><p>IFE shield material fractions</p></td>
    </tr>
    
    <tr>
        <td>shmatm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE shield material masses (kg)</p></td>
    </tr>
    
    <tr>
        <td>shmatv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE shield material volumes (kg)</p></td>
    </tr>
    
    <tr>
        <td>shvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE shield volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>sombdr</td>
        <td>Input</td>
        <td>real</td>
        <td>2.7</td>
        <td><p>radius of cylindrical blanket section below chamber (<code>ifetyp=2</code>)</p></td>
    </tr>
    
    <tr>
        <td>somtdr</td>
        <td>Input</td>
        <td>real</td>
        <td>2.7</td>
        <td><p>radius of cylindrical blanket section above chamber (<code>ifetyp=2</code>)</p></td>
    </tr>
    
    <tr>
        <td>taufall</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Lithium Fall Time (s)</p></td>
    </tr>
    
    <tr>
        <td>tdspmw</td>
        <td>Input</td>
        <td>real</td>
        <td>0.01</td>
        <td><p>IFE target delivery system power (MW)</p></td>
    </tr>
    
    <tr>
        <td>tfacmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE target factory power (MW)</p></td>
    </tr>
    
    <tr>
        <td>tgain</td>
        <td>Input</td>
        <td>real</td>
        <td>85.0</td>
        <td><p>IFE target gain (if <code>ifedrv = 0</code>) (<code>iteration variable 83</code>)</p></td>
    </tr>
    
    <tr>
        <td>uccarb</td>
        <td>Input</td>
        <td>real</td>
        <td>50.0</td>
        <td><p>cost of carbon cloth ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucconc</td>
        <td>Input</td>
        <td>real</td>
        <td>0.1</td>
        <td><p>cost of concrete ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>ucflib</td>
        <td>Input</td>
        <td>real</td>
        <td>84.0</td>
        <td><p>cost of FLiBe ($/kg)</p></td>
    </tr>
    
    <tr>
        <td>uctarg</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>cost of IFE target ($/target)</p></td>
    </tr>
    
    <tr>
        <td>v1dr</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radial thickness of IFE void between first wall and blanket (m)</p></td>
    </tr>
    
    <tr>
        <td>v1dzl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vertical thickness of IFE void 1 below chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>v1dzu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vertical thickness of IFE void 1 above chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>v1matf</td>
        <td>Input</td>
        <td>real</td>
        <td>[[1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0.]]</td>
        <td><p>IFE void 1 material fractions</p></td>
    </tr>
    
    <tr>
        <td>v1matm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE void 1 material masses (kg)</p></td>
    </tr>
    
    <tr>
        <td>v1matv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE void 1 material volumes (kg)</p></td>
    </tr>
    
    <tr>
        <td>v1vol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE void 1 volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>v2dr</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>radial thickness of IFE void between blanket and shield (m)</p></td>
    </tr>
    
    <tr>
        <td>v2dzl</td>
        <td>Input</td>
        <td>real</td>
        <td>7.0</td>
        <td><p>vertical thickness of IFE void 2 below chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>v2dzu</td>
        <td>Input</td>
        <td>real</td>
        <td>7.0</td>
        <td><p>vertical thickness of IFE void 2 above chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>v2matf</td>
        <td>Input</td>
        <td>real</td>
        <td>[[1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0.]]</td>
        <td><p>IFE void 2 material fractions</p></td>
    </tr>
    
    <tr>
        <td>v2matm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE void 2 material masses (kg)</p></td>
    </tr>
    
    <tr>
        <td>v2matv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE void 2 material volumes (kg)</p></td>
    </tr>
    
    <tr>
        <td>v2vol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE void 2 volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>v3dr</td>
        <td>Input</td>
        <td>real</td>
        <td>43.3</td>
        <td><p>radial thickness of IFE void outside shield (m)</p></td>
    </tr>
    
    <tr>
        <td>v3dzl</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>vertical thickness of IFE void 3 below chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>v3dzu</td>
        <td>Input</td>
        <td>real</td>
        <td>20.0</td>
        <td><p>vertical thickness of IFE void 3 above chamber (m)</p></td>
    </tr>
    
    <tr>
        <td>v3matf</td>
        <td>Input</td>
        <td>real</td>
        <td>[[1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0.]
 [1. 0. 0. 0. 0. 0. 0. 0. 0.]]</td>
        <td><p>IFE void 3 material fractions</p></td>
    </tr>
    
    <tr>
        <td>v3matm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE void 3 material masses (kg)</p></td>
    </tr>
    
    <tr>
        <td>v3matv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE void 3 material volumes (kg)</p></td>
    </tr>
    
    <tr>
        <td>v3vol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE void 3 volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>zl1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build below centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zl2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build below centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zl3</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build below centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zl4</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build below centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zl5</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build below centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zl6</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build below centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zl7</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build below centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zu1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build above centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zu2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build above centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zu3</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build above centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zu4</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build above centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zu5</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build above centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zu6</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build above centre (m)</p></td>
    </tr>
    
    <tr>
        <td>zu7</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>IFE vertical build above centre (m)</p></td>
    </tr>
    
</table>

## impurity_radiation_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>n_impurities</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>14</td>
        <td><p>n_impurities /14/ FIX : number of ion species in impurity radiation model</p></td>
    </tr>
    
    <tr>
        <td>coreradius</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6</td>
        <td><p>coreradius /0.6/ : normalised radius defining the 'core' region</p></td>
    </tr>
    
    <tr>
        <td>coreradiationfraction</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>coreradiationfraction /1.0/ : fraction of radiation from 'core' region that is subtracted from the loss power</p>
<p>fimp(n_impurities) /1.0,0.1,0.02,0.0,0.0,0.0,0.0,0.0,0.0016,0.0,0.0,0.0,0.0,0.0/ :
        impurity number density fractions relative to electron density</p></td>
    </tr>
    
    <tr>
        <td>fimp</td>
        <td>Input</td>
        <td>real</td>
        <td>[1.  0.1 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0. ]</td>
        <td></td>
    </tr>
    
    <tr>
        <td>imp_label</td>
        <td>Input</td>
        <td>character</td>
        <td>[b'H_' b'He' b'Be' b'C_' b'N_' b'O_' b'Ne' b'Si' b'Ar' b'Fe' b'Ni' b'Kr'
 b'Xe' b'W_']</td>
        <td><p>imp_label(n_impurities) : impurity ion species names:<UL>
 
<LI> ( 1)  Hydrogen  (fraction calculated by code)
 <LI> ( 2)  Helium
 <LI> ( 3)  Beryllium
 <LI> ( 4)  Carbon
 <LI> ( 5)  Nitrogen
 <LI> ( 6)  Oxygen
 <LI> ( 7)  Neon
 <LI> ( 8)  Silicon
 <LI> ( 9)  Argon
 <LI> (10)  Iron
 <LI> (11)  Nickel
 <LI> (12)  Krypton
 <LI> (13)  Xenon
 <LI> (14)  Tungsten</UL>

























</p></td>
    </tr>
    
    <tr>
        <td>all_array_hotfix_len*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>200</td>
        <td></td>
    </tr>
    
    <tr>
        <td>impurity_arr_label</td>
        <td>Input</td>
        <td>character</td>
        <td>[b'  ' b'  ' b'  ' b'  ' b'  ' b'  ' b'  ' b'  ' b'  ' b'  ' b'  ' b'  '
 b'  ' b'  ']</td>
        <td></td>
    </tr>
    
    <tr>
        <td>impurity_arr_z</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>impurity_arr_amass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>impurity_arr_frac</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>impurity_arr_len_tab</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>impurity_arr_temp_kev</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>impurity_arr_lz_wm3</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>impurity_arr_zav</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>toolow</td>
        <td>Output</td>
        <td>logical</td>
        <td>-</td>
        <td><p>Used for reporting error in function pimpden</p></td>
    </tr>
    
</table>

## process_input
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>nin</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>10</td>
        <td></td>
    </tr>
    
    <tr>
        <td>maxlen*</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>2000</td>
        <td></td>
    </tr>
    
    <tr>
        <td>line*</td>
        <td>Variable</td>
        <td>character</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>linelen*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>lineno*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>iptr*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>infile*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>outfile*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>report_changes*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>icode*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>subscript_present*</td>
        <td>Variable</td>
        <td>logical</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error*</td>
        <td>Variable</td>
        <td>logical</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>error_message*</td>
        <td>Variable</td>
        <td>character</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>show_changes*</td>
        <td>Variable</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>constraints_exist*</td>
        <td>Variable</td>
        <td>logical</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## define_iteration_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>dummy</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## neoclassics_constants
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>no_roots</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>30</td>
        <td></td>
    </tr>
    
</table>

## neoclassics_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>species</td>
        <td>Variable</td>
        <td>character</td>
        <td>(/"e", "D", "T", "a"/)</td>
        <td></td>
    </tr>
    
    <tr>
        <td>densities</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>temperatures</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dr_densities</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>dr_temperatures</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>roots</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>weights</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>nu</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>nu_star</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>nu_star_averaged</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>vd</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>kt</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>er</td>
        <td>Variable</td>
        <td>real</td>
        <td>0.0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>iota</td>
        <td>Variable</td>
        <td>real</td>
        <td>1.0d0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>d11_mono</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>d11_plateau</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>d111</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>d112</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>d113</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>q_flux</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>gamma_flux</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>d31_mono</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
    <tr>
        <td>eps_eff</td>
        <td>Variable</td>
        <td>real</td>
        <td>1d-5</td>
        <td></td>
    </tr>
    
    <tr>
        <td>r_eff</td>
        <td>Variable</td>
        <td>real</td>
        <td>0</td>
        <td></td>
    </tr>
    
</table>

## numerics
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>ipnvars</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>175</td>
        <td><p>ipnvars FIX : total number of variables available for iteration</p></td>
    </tr>
    
    <tr>
        <td>ipeqns</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>92</td>
        <td><p>ipeqns  FIX : number of constraint equations available</p></td>
    </tr>
    
    <tr>
        <td>ipnfoms</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>19</td>
        <td><p>ipnfoms FIX : number of available figures of merit</p></td>
    </tr>
    
    <tr>
        <td>ipvlam</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>ipeqns+2*ipnvars+1</td>
        <td></td>
    </tr>
    
    <tr>
        <td>iptnt</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>(ipeqns*(3*ipeqns+13))/2</td>
        <td></td>
    </tr>
    
    <tr>
        <td>ipvp1</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>ipnvars+1</td>
        <td></td>
    </tr>
    
    <tr>
        <td>ioptimz</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>ioptimz /1/ : code operation switch:<UL>
           <LI> = -2 for no optimisation, no VMCOM or HYBRD;
           <LI> = -1 for no optimisation, HYBRD only;
           <LI> = 0  for HYBRD and VMCON (not recommended);
           <LI> = 1  for optimisation, VMCON only</UL></p>
<p>minmax /7/ : switch for figure-of-merit (see lablmm for descriptions)
               negative =&gt; maximise, positive =&gt; minimise</p></td>
    </tr>
    
    <tr>
        <td>minmax</td>
        <td>Input</td>
        <td>integer</td>
        <td>7</td>
        <td></td>
    </tr>
    
    <tr>
        <td>lablmm</td>
        <td>Input</td>
        <td>character</td>
        <td>[b'major radius          ' b'not used              '
 b'neutron wall load     ' b'P_tf + P_pf           '
 b'fusion gain           ' b'cost of electricity   '
 b'capital cost          ' b'aspect ratio          '
 b'divertor heat load    ' b'toroidal field        '
 b'total injected power  ' b'H plant capital cost  '
 b'H production rate     ' b'pulse length          '
 b'plant availability    ' b'min R0, max tau_burn  '
 b'net electrical output ' b'Null figure of merit  '
 b'max Q, max t_burn     ']</td>
        <td><p>lablmm(ipnfoms) : labels describing figures of merit:<UL>
<br>
<LI> ( 1) major radius
  <LI> ( 2) not used
  <LI> ( 3) neutron wall load
  <LI> ( 4) P_tf + P_pf
  <LI> ( 5) fusion gain Q
  <LI> ( 6) cost of electricity
  <LI> ( 7) capital cost (direct cost if ireactor=0,
                          constructed cost otherwise)
  <LI> ( 8) aspect ratio
  <LI> ( 9) divertor heat load
  <LI> (10) toroidal field
  <LI> (11) total injected power
  <LI> (12) hydrogen plant capital cost OBSOLETE
  <LI> (13) hydrogen production rate OBSOLETE
  <LI> (14) pulse length
  <LI> (15) plant availability factor (N.B. requires
            iavail=1 to be set)
  <LI> (16) linear combination of major radius (minimised) and pulse length (maximised)
              note: FoM should be minimised only!
  <LI> (17) net electrical output
  <LI> (18) Null Figure of Merit
  <LI> (19) linear combination of big Q and pulse length (maximised)
              note: FoM should be minimised only!</UL>


</p></td>
    </tr>
    
    <tr>
        <td>ncalls</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>ncalls : number of function calls during solution</p></td>
    </tr>
    
    <tr>
        <td>neqns</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>neqns /0/ : number of equality constraints to be satisfied</p></td>
    </tr>
    
    <tr>
        <td>nfev1</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>nfev1 : number of calls to FCNHYB (HYBRD function caller) made</p></td>
    </tr>
    
    <tr>
        <td>nfev2</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>nfev2 : number of calls to FCNVMC1 (VMCON function caller) made</p></td>
    </tr>
    
    <tr>
        <td>nineqns</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>nineqns /0/ : number of inequality constraints VMCON must satisfy
                (leave at zero for now)</p></td>
    </tr>
    
    <tr>
        <td>nvar</td>
        <td>Input</td>
        <td>integer</td>
        <td>16</td>
        <td><p>nvar /16/ : number of iteration variables to use</p></td>
    </tr>
    
    <tr>
        <td>nviter</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>nviter : number of VMCON iterations performed</p>
<p>icc(ipeqns) /0/ :
           array defining which constraint equations to activate
           (see lablcc for descriptions)</p></td>
    </tr>
    
    <tr>
        <td>icc</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>active_constraints</td>
        <td>Output</td>
        <td>logical</td>
        <td>-</td>
        <td><p>active_constraints(ipeqns) : Logical array showing which constraints are active</p></td>
    </tr>
    
    <tr>
        <td>lablcc</td>
        <td>Input</td>
        <td>character</td>
        <td>[b'Beta consistency                 ' b'Global power balance consistency '
 b'Ion power balance                ' b'Electron power balance           '
 b'Density upper limit              ' b'(Epsilon x beta-pol) upper limit '
 b'Beam ion density consistency     ' b'Neutron wall load upper limit    '
 b'Fusion power upper limit         ' b'Toroidal field 1/R consistency   '
 b'Radial build consistency         ' b'Volt second lower limit          '
 b'Burn time lower limit            ' b'NBI decay lengths consistency    '
 b'L-H power threshold limit        ' b'Net electric power lower limit   '
 b'Radiation fraction upper limit   ' b'Divertor heat load upper limit   '
 b'MVA upper limit                  ' b'Beam tangency radius upper limit '
 b'Plasma minor radius lower limit  ' b'Divertor collisionality upper lim'
 b'Conducting shell radius upper lim' b'Beta upper limit                 '
 b'Peak toroidal field upper limit  ' b'CS coil EOF current density limit'
 b'CS coil BOP current density limit' b'Fusion gain Q lower limit        '
 b'Inboard radial build consistency ' b'Injection power upper limit      '
 b'TF coil case stress upper limit  ' b'TF coil conduit stress upper lim '
 b'I_op / I_critical (TF coil)      ' b'Dump voltage upper limit         '
 b'J_winding pack/J_protection limit' b'TF coil temp. margin lower limit '
 b'Current drive gamma limit        ' b'1st wall coolant temp rise limit '
 b'First wall peak temperature limit' b'Start-up inj. power lower limit  '
 b'Plasma curr. ramp time lower lim ' b'Cycle time lower limit           '
 b'Average centrepost temperature   ' b'Peak centrepost temp. upper limit'
 b'Edge safety factor lower limit   ' b'Ip/Irod upper limit              '
 b'TF coil tor. thickness upper lim ' b'Poloidal beta upper limit        '
 b'RFP reversal parameter < 0       ' b'IFE repetition rate upper limit  '
 b'Startup volt-seconds consistency ' b'Tritium breeding ratio lower lim '
 b'Neutron fluence on TF coil limit ' b'Peak TF coil nucl. heating limit '
 b'Vessel helium concentration limit' b'Psep / R upper limit             '
 b'TF coil leg rad width lower limit' b'TF coil leg rad width lower limit'
 b'NB shine-through frac upper limit' b'CS temperature margin lower limit'
 b'Minimum availability value       ' b'f_alpha_energy_confinement       '
 b'number of ITER-like vacuum pumps ' b'Zeff limit                       '
 b'Dump time set by VV stress       ' b'Rate of change of energy in field'
 b'Upper Lim. on Radiation Wall load' b'Upper Lim. on Psep * Bt / q A R  '
 b'pdivt < psep_kallenbach divertor ' b'Separatrix temp consistency      '
 b'Separatrix density consistency   ' b'CS Tresca yield criterion        '
 b'Psep >= Plh + Paux               ' b'TFC quench < tmax_croco          '
 b'TFC current/copper area < Max    ' b'Eich critical separatrix density '
 b'TFC current per turn upper limit ' b'Reinke criterion fZ lower limit  '
 b'Peak CS field upper limit        ' b'pdivt lower limit                '
 b'ne0 > neped                      ' b'toroidalgap >  tftort            '
 b'available_space > required_space ' b'beta > beta_min                  '
 b'CP lifetime                      ' b'TFC turn dimension               '
 b'Cryogenic plant power            ' b'TF coil strain absolute value    '
 b'CS current/copper area < Max     ' b'CS stress load cycles            '
 b'ECRH ignitability                ' b'Fuel composition consistency     ']</td>
        <td><p>lablcc(ipeqns) : labels describing constraint equations (corresponding itvs)<UL>
<br>
<LI> ( 1) Beta (consistency equation) (itv 5)
  <LI> ( 2) Global power balance (consistency equation) (itv 10,1,2,3,4,6,11)
  <LI> ( 3) Ion power balance DEPRECATED (itv 10,1,2,3,4,6,11)
  <LI> ( 4) Electron power balance DEPRECATED (itv 10,1,2,3,4,6,11)
  <LI> ( 5) Density upper limit (itv 9,1,2,3,4,5,6)
  <LI> ( 6) (Epsilon x beta poloidal) upper limit (itv 8,1,2,3,4,6)
  <LI> ( 7) Beam ion density (NBI) (consistency equation) (itv 7)
  <LI> ( 8) Neutron wall load upper limit (itv 14,1,2,3,4,6)
  <LI> ( 9) Fusion power upper limit (itv 26,1,2,3,4,6)
  <LI> (10) Toroidal field 1/R (consistency equation) (itv 12,1,2,3,13 )
  <LI> (11) Radial build (consistency equation) (itv 3,1,13,16,29,42,61)
  <LI> (12) Volt second lower limit (STEADY STATE) (itv 15,1,2,3)
  <LI> (13) Burn time lower limit (PULSE) (itv 21,1,16,17,29,42,44,61)
            (itv 19,1,2,3,6)
  <LI> (14) Neutral beam decay lengths to plasma centre (NBI) (consistency equation)
  <LI> (15) LH power threshold limit (itv 103)
  <LI> (16) Net electric power lower limit (itv 25,1,2,3)
  <LI> (17) Radiation fraction upper limit (itv 28)
  <LI> (18) Divertor heat load upper limit (itv 27)
  <LI> (19) MVA upper limit (itv 30)
  <LI> (20) Neutral beam tangency radius upper limit (NBI) (itv 33,31,3,13)
  <LI> (21) Plasma minor radius lower limit (itv 32)
  <LI> (22) Divertor collisionality upper limit (itv 34,43)
  <LI> (23) Conducting shell to plasma minor radius ratio upper limit
            (itv 104,1,74)
  <LI> (24) Beta upper limit (itv 36,1,2,3,4,6,18)
  <LI> (25) Peak toroidal field upper limit (itv 35,3,13,29)
  <LI> (26) Central solenoid EOF current density upper limit (ipfres=0)
            (itv 38,37,41,12)
  <LI> (27) Central solenoid BOP current density upper limit (ipfres=0)
            (itv 39,37,41,12)
  <LI> (28) Fusion gain Q lower limit (itv 45,47,40)
  <LI> (29) Inboard radial build consistency (itv 3,1,13,16,29,42,61)
  <LI> (30) Injection power upper limit (itv 46,47,11)
  <LI> (31) TF coil case stress upper limit (SCTF) (itv 48,56,57,58,59,60,24)
  <LI> (32) TF coil conduit stress upper limit (SCTF) (itv 49,56,57,58,59,60,24)
  <LI> (33) I_op / I_critical (TF coil) (SCTF) (itv 50,56,57,58,59,60,24)
  <LI> (34) Dump voltage upper limit (SCTF) (itv 51,52,56,57,58,59,60,24)
  <LI> (35) J_winding pack/J_protection upper limit (SCTF) (itv 53,56,57,58,59,60,24)
  <LI> (36) TF coil temperature margin lower limit (SCTF) (itv 54,55,56,57,58,59,60,24)
  <LI> (37) Current drive gamma upper limit (itv 40,47)
  <LI> (38) First wall coolant temperature rise upper limit (itv 62)
  <LI> (39) First wall peak temperature upper limit (itv 63)
  <LI> (40) Start-up injection power lower limit (PULSE) (itv 64)
  <LI> (41) Plasma current ramp-up time lower limit (PULSE) (itv  66,65)
  <LI> (42) Cycle time lower limit (PULSE) (itv 17,67,65)
  <LI> (43) Average centrepost temperature
            (TART) (consistency equation) (itv 13,20,69,70)
  <LI> (44) Peak centrepost temperature upper limit (TART) (itv 68,69,70)
  <LI> (45) Edge safety factor lower limit (TART) (itv 71,1,2,3)
  <LI> (46) Equation for Ip/Irod upper limit (TART) (itv 72,2,60)
  <LI> (47) NOT USED
  <LI> (48) Poloidal beta upper limit (itv 79,2,3,18)
  <LI> (49) NOT USED
  <LI> (50) IFE repetition rate upper limit (IFE)
  <LI> (51) Startup volt-seconds consistency (PULSE) (itv 16,29,3,1)
  <LI> (52) Tritium breeding ratio lower limit (itv 89,90,91)
  <LI> (53) Neutron fluence on TF coil upper limit (itv 92,93,94)
  <LI> (54) Peak TF coil nuclear heating upper limit (itv 95,93,94)
  <LI> (55) Vacuum vessel helium concentration upper limit i_blanket_type =2 (itv 96,93,94)
  <LI> (56) Pseparatrix/Rmajor upper limit (itv 97,1,3)
  <LI> (57) NOT USED
  <LI> (58) NOT USED
  <LI> (59) Neutral beam shine-through fraction upper limit (NBI) (itv 105,6,19,4 )
  <LI> (60) Central solenoid temperature margin lower limit (SCTF) (itv 106)
  <LI> (61) Minimum availability value (itv 107)
  <LI> (62) f_alpha_energy_confinement the ratio of particle to energy confinement times (itv 110)
  <LI> (63) The number of ITER-like vacuum pumps niterpump < tfno (itv 111)
  <LI> (64) Zeff less than or equal to zeffmax (itv 112)
  <LI> (65) Dump time set by VV loads (itv 56, 113)
  <LI> (66) Limit on rate of change of energy in poloidal field
            (Use iteration variable 65(t_current_ramp_up), 115)
  <LI> (67) Simple Radiation Wall load limit (itv 116, 4,6)
  <LI> (68) Psep * Bt / qAR upper limit (itv 117)
  <LI> (69) ensure separatrix power = the value from Kallenbach divertor (itv 118)
  <LI> (70) ensure that teomp = separatrix temperature in the pedestal profile,
            (itv 119 (tesep))
  <LI> (71) ensure that neomp = separatrix density (nesep) x neratio
  <LI> (72) central solenoid shear stress limit (Tresca yield criterion) (itv 123 foh_stress)
  <LI> (73) Psep >= Plh + Paux (itv 137 (fplhsep))
  <LI> (74) TFC quench < tmax_croco (itv 141 (fcqt))
  <LI> (75) TFC current/copper area < Maximum (itv 143 f_coppera_m2)
  <LI> (76) Eich critical separatrix density
  <LI> (77) TF coil current per turn upper limit
  <LI> (78) Reinke criterion impurity fraction lower limit (itv  147 freinke)
  <LI> (79) Peak CS field upper limit (itv  149 fbmaxcs)
  <LI> (80) Divertor power lower limit pdivt (itv  153 fpdivlim)
  <LI> (81) Ne(0) > ne(ped) constraint (itv  154 fne0)
  <LI> (82) toroidalgap >  tftort constraint (itv  171 ftoroidalgap)
  <LI> (83) Radial build consistency for stellarators (itv 172 f_avspace)
  <LI> (84) Lower limit for beta (itv 173 fbeta_min)
  <LI> (85) Constraint for CP lifetime
  <LI> (86) Constraint for TF coil turn dimension
  <LI> (87) Constraint for cryogenic power
  <LI> (88) Constraint for TF coil strain absolute value
  <LI> (89) Constraint for CS coil quench protection
  <LI> (90) Lower Limit on number of stress load cycles for CS (itr 167 fncycle)
  <LI> (91) Checking if the design point is ECRH ignitable (itv 168 fecrh_ignition)
  <LI> (92) D/T/He3 ratio in fuel sums to 1 </UL>



</p></td>
    </tr>
    
    <tr>
        <td>ixc</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>ixc(ipnvars) /0/ :
               array defining which iteration variables to activate
               (see lablxc for descriptions)</p></td>
    </tr>
    
    <tr>
        <td>lablxc</td>
        <td>Input</td>
        <td>character</td>
        <td>[b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ' b'                              '
 b'                              ']</td>
        <td><p>lablxc(ipnvars) : labels describing iteration variables<UL>
<br>
<LI> ( 1) aspect
  <LI> ( 2) bt
  <LI> ( 3) rmajor
 <LI> ( 4) te
 <LI> ( 5) beta
 <LI> ( 6) dene
 <LI> ( 7) f_nd_beam_electron
 <LI> ( 8) fbeta_poloidal_eps (f-value for equation 6)
 <LI> ( 9) fdene (f-value for equation 5)
 <LI> (10) hfact
 <LI> (11) pheat
 <LI> (12) oacdcp
 <LI> (13) dr_tf_inboard (NOT RECOMMENDED)
 <LI> (14) fwalld (f-value for equation 8)
 <LI> (15) fvs (f-value for equation 12)
 <LI> (16) dr_cs
 <LI> (17) t_between_pulse
 <LI> (18) q
 <LI> (19) beam_energy
 <LI> (20) temp_cp_average
 <LI> (21) ft_burn (f-value for equation 13)
 <LI> (22) NOT USED
 <LI> (23) fcoolcp
 <LI> (24) NOT USED
 <LI> (25) fpnetel (f-value for equation 16)
 <LI> (26) ffuspow (f-value for equation 9)
 <LI> (27) fhldiv (f-value for equation 18)
 <LI> (28) fradpwr (f-value for equation 17), total radiation fraction
 <LI> (29) dr_bore
 <LI> (30) fmva (f-value for equation 19)
 <LI> (31) gapomin
 <LI> (32) frminor (f-value for equation 21)
 <LI> (33) fportsz (f-value for equation 20)
 <LI> (34) fdivcol (f-value for equation 22)
 <LI> (35) fpeakb (f-value for equation 25)
 <LI> (36) fbeta_max (f-value for equation 24)
 <LI> (37) coheof
 <LI> (38) fjohc (f-value for equation 26)
 <LI> (39) fjohc0 (f-value for equation 27)
 <LI> (40) fgamcd (f-value for equation 37)
 <LI> (41) fcohbop
 <LI> (42) dr_cs_tf_gap
 <LI> (43) NOT USED
 <LI> (44) fvsbrnni
 <LI> (45) fqval (f-value for equation 28)
 <LI> (46) fpinj (f-value for equation 30)
 <LI> (47) feffcd
 <LI> (48) fstrcase (f-value for equation 31)
 <LI> (49) fstrcond (f-value for equation 32)
 <LI> (50) fiooic (f-value for equation 33)
 <LI> (51) fvdump (f-value for equation 34)
 <LI> (52) NOT USED
 <LI> (53) fjprot (f-value for equation 35)
 <LI> (54) ftmargtf (f-value for equation 36)
 <LI> (55) NOT USED
 <LI> (56) tdmptf
 <LI> (57) thkcas
 <LI> (58) thwcndut
 <LI> (59) fcutfsu
 <LI> (60) cpttf
 <LI> (61) dr_shld_vv_gap_inboard
 <LI> (62) fdtmp (f-value for equation 38)
 <LI> (63) ftpeak (f-value for equation 39)
 <LI> (64) fauxmn (f-value for equation 40)
 <LI> (65) t_current_ramp_up
 <LI> (66) ft_current_ramp_up (f-value for equation 41)
 <LI> (67) ftcycl (f-value for equation 42)
 <LI> (68) fptemp (f-value for equation 44)
 <LI> (69) rcool
 <LI> (70) vcool
 <LI> (71) fq (f-value for equation 45)
 <LI> (72) fipir (f-value for equation 46)
 <LI> (73) dr_fw_plasma_gap_inboard
 <LI> (74) dr_fw_plasma_gap_outboard
 <LI> (75) tfootfi
 <LI> (76) NOT USED
 <LI> (77) NOT USED
 <LI> (78) NOT USED
 <LI> (79) fbeta_poloidal (f-value for equation 48)
 <LI> (80) NOT USED
 <LI> (81) edrive
 <LI> (82) drveff
 <LI> (83) tgain
 <LI> (84) chrad
 <LI> (85) pdrive
 <LI> (86) frrmax (f-value for equation 50)
 <LI> (87) NOT USED
 <LI> (88) NOT USED
 <LI> (89) ftbr (f-value for equation 52)
 <LI> (90) blbuith
 <LI> (91) blbuoth
 <LI> (92) fflutf (f-value for equation 53)
 <LI> (93) dr_shld_inboard
 <LI> (94) dr_shld_outboard
 <LI> (95) fptfnuc (f-value for equation 54)
 <LI> (96) fvvhe (f-value for equation 55)
 <LI> (97) fpsepr (f-value for equation 56)
 <LI> (98) li6enrich
 <LI> (99) NOT USED
 <LI> (100) NOT USED
 <LI> (101) NOT USED
 <LI> (102) fimpvar # OBSOLETE
 <LI> (103) fl_h_threshold (f-value for equation 15)
 <LI> (104)fr_conducting_wall (f-value for equation 23)
 <LI> (105) fnbshinef (f-value for equation 59)
 <LI> (106) ftmargoh (f-value for equation 60)
 <LI> (107) favail (f-value for equation 61)
 <LI> (108) breeder_f: Volume of Li4SiO4 / (Volume of Be12Ti + Li4SiO4)
 <LI> (109) f_nd_alpha_electron: thermal alpha density / electron density
 <LI> (110) falpha_energy_confinement: Lower limit on f_alpha_energy_confinement the ratio of alpha
 <LI> (111) fniterpump: f-value for constraint that number
 <LI> (112) fzeffmax: f-value for max Zeff (f-value for equation 64)
 <LI> (113) ftaucq: f-value for minimum quench time (f-value for equation 65)
 <LI> (114) len_fw_channel: Length of a single first wall channel
 <LI> (115) fpoloidalpower: f-value for max rate of change of
 <LI> (116) fradwall: f-value for radiation wall load limit (eq. 67)
 <LI> (117) fpsepbqar: f-value for  Psep*Bt/qar upper limit (eq. 68)
 <LI> (118) fpsep: f-value to ensure separatrix power is less than
 <LI> (119) tesep:  separatrix temperature calculated by the Kallenbach divertor model
 <LI> (120) ttarget: Plasma temperature adjacent to divertor sheath [eV]
 <LI> (121) neratio: ratio of mean SOL density at OMP to separatrix density at OMP
 <LI> (122) oh_steel_frac : streel fraction of Central Solenoid
 <LI> (123) foh_stress : f-value for CS coil Tresca yield criterion (f-value for eq. 72)
 <LI> (124) qtargettotal : Power density on target including surface recombination [W/m2]
 <LI> (125) fimp(3) :  Beryllium density fraction relative to electron density
 <LI> (126) fimp(4) :  Carbon density fraction relative to electron density
 <LI> (127) fimp(5) :  Nitrogen fraction relative to electron density
 <LI> (128) fimp(6) :  Oxygen density fraction relative to electron density
 <LI> (129) fimp(7) :  Neon density fraction relative to electron density
 <LI> (130) fimp(8) :  Silicon density fraction relative to electron density
 <LI> (131) fimp(9) :  Argon density fraction relative to electron density
 <LI> (132) fimp(10) :  Iron density fraction relative to electron density
 <LI> (133) fimp(11) :  Nickel density fraction relative to electron density
 <LI> (134) fimp(12) :  Krypton density fraction relative to electron density
 <LI> (135) fimp(13) :  Xenon density fraction relative to electron density
 <LI> (136) fimp(14) :  Tungsten density fraction relative to electron density
 <LI> (137) fplhsep (f-value for equation 73)
 <LI> (138) rebco_thickness : thickness of REBCO layer in tape (m)
 <LI> (139) copper_thick : thickness of copper layer in tape (m)
 <LI> (140) dr_tf_wp : radial thickness of TFC winding pack (m)
 <LI> (141) fcqt : TF coil quench temperature < tmax_croco (f-value for equation 74)
 <LI> (142) nesep : electron density at separatrix [m-3]
 <LI> (143) f_copperA_m2 : TF coil current / copper area < Maximum value
 <LI> (144) fnesep : Eich critical electron density at separatrix
 <LI> (145) fgwped :  fraction of Greenwald density to set as pedestal-top density
 <LI> (146) fcpttf : F-value for TF coil current per turn limit (constraint equation 77)
 <LI> (147) freinke : F-value for Reinke detachment criterion (constraint equation 78)
 <LI> (148) fzactual : fraction of impurity at SOL with Reinke detachment criterion
 <LI> (149) fbmaxcs : F-value for max peak CS field (con. 79, itvar 149)
 <LI> (150) REMOVED
 <LI> (151) REMOVED
 <LI> (152) fgwsep : Ratio of separatrix density to Greenwald density
 <LI> (153) fpdivlim : F-value for minimum pdivt (con. 80)
 <LI> (154) fne0 : F-value for ne(0) > ne(ped) (con. 81)
 <LI> (155) pfusife : IFE input fusion power (MW) (ifedrv=3 only)
 <LI> (156) rrin : Input IFE repetition rate (Hz) (ifedrv=3 only)
 <LI> (157) fvssu : F-value for available to required start up flux (con. 51)
 <LI> (158) croco_thick : Thickness of CroCo copper tube (m)
 <LI> (159) ftoroidalgap : F-value for toroidalgap >  tftort constraint (con. 82)
 <LI> (160) f_avspace (f-value for equation 83)
 <LI> (161) fbeta_min (f-value for equation 84)
 <LI> (162) r_cp_top : Top outer radius of the centropost (ST only) (m)
 <LI> (163) f_t_turn_tf : f-value for TF coils WP trurn squared dimension constraint
 <LI> (164) f_crypmw : f-value for cryogenic plant power
 <LI> (165) fstr_wp : f-value for TF coil strain absolute value
 <LI> (166) f_copperaoh_m2 : CS coil current /copper area < Maximum value
 <LI> (167) fncycle : f-value for minimum CS coil stress load cycles
 <LI> (168) fecrh_ignition: f-value for equation 91
 <LI> (169) te0_ecrh_achievable: Max. achievable electron temperature at ignition point
 <LI> (170) beta_div : field line angle wrt divertor target plate (degrees)
 <LI> (171) casths_fraction : TF side case thickness as fraction of toridal case thickness
 <LI> (172) casths : TF side case thickness [m]
 <LI> (173) f_deuterium : Deuterium fraction in fuel
 <LI> (174) EMPTY : Description
 <LI> (175) EMPTY : Description



</p></td>
    </tr>
    
    <tr>
        <td>name_xc</td>
        <td>Variable</td>
        <td>character</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>sqsumsq</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>sqsumsq : sqrt of the sum of the square of the constraint residuals</p></td>
    </tr>
    
    <tr>
        <td>objf_name</td>
        <td>Input</td>
        <td>character</td>
        <td>b'                                        '</td>
        <td><p>Description of the objective function</p></td>
    </tr>
    
    <tr>
        <td>norm_objf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Normalised objective function (figure of merit)</p></td>
    </tr>
    
    <tr>
        <td>epsfcn</td>
        <td>Input</td>
        <td>real</td>
        <td>0.001</td>
        <td><p>epsfcn /1.0e-3/ : finite difference step length for HYBRD/VMCON derivatives</p></td>
    </tr>
    
    <tr>
        <td>epsvmc</td>
        <td>Input</td>
        <td>real</td>
        <td>1e-06</td>
        <td><p>epsvmc /1.0e-6/ : error tolerance for VMCON</p></td>
    </tr>
    
    <tr>
        <td>factor</td>
        <td>Input</td>
        <td>real</td>
        <td>0.1</td>
        <td><p>factor /0.1/ : used in HYBRD for first step size</p></td>
    </tr>
    
    <tr>
        <td>ftol</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0001</td>
        <td><p>ftol /1.0e-4/ : error tolerance for HYBRD</p></td>
    </tr>
    
    <tr>
        <td>boundl</td>
        <td>Input</td>
        <td>real</td>
        <td>[9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99
 9.e-99 9.e-99 9.e-99 9.e-99 9.e-99]</td>
        <td><p>boundl(ipnvars) /../ : lower bounds used on ixc variables during
                         VMCON optimisation runs</p></td>
    </tr>
    
    <tr>
        <td>boundu</td>
        <td>Input</td>
        <td>real</td>
        <td>[9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99
 9.e+99 9.e+99 9.e+99 9.e+99 9.e+99]</td>
        <td></td>
    </tr>
    
    <tr>
        <td>itv_scaled_lower_bounds</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Lower bound of the ixc variables scaled to (divided by)
 the initial value of the corresponding ixc</p></td>
    </tr>
    
    <tr>
        <td>itv_scaled_upper_bounds</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Upper bound of the ixc variables scaled to (divided by)
 the initial value of the corresponding ixc</p></td>
    </tr>
    
    <tr>
        <td>rcm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>resdl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>scafc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>The initial value of each ixc variable</p></td>
    </tr>
    
    <tr>
        <td>scale</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>The reciprocal of the initial value of each ixc variable</p></td>
    </tr>
    
    <tr>
        <td>xcm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>xcs</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>vlam</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## pf_power_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>acptmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>average of currents in PF circuits (kA)</p></td>
    </tr>
    
    <tr>
        <td>ensxpfm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum stored energy in the PF circuits (MJ)</p></td>
    </tr>
    
    <tr>
        <td>iscenr</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>Switch for PF coil energy storage option:</p>
<ul>
<li>=1 all power from MGF (motor-generator flywheel) units</li>
<li>=2 all pulsed power from line</li>
<li>=3 PF power from MGF, heating from line</li>
</ul></td>
    </tr>
    
    <tr>
        <td>pfckts</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>number of PF coil circuits</p></td>
    </tr>
    
    <tr>
        <td>spfbusl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total PF coil circuit bus length (m)</p></td>
    </tr>
    
    <tr>
        <td>spsmva</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>sum of PF power supply ratings (MVA)</p></td>
    </tr>
    
    <tr>
        <td>srcktpm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>sum of resistive PF coil power (kW)</p></td>
    </tr>
    
    <tr>
        <td>vpfskv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>PF coil voltage (kV)</p></td>
    </tr>
    
    <tr>
        <td>peakpoloidalpower</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Peak absolute rate of change of stored energy in poloidal field (MW)</p></td>
    </tr>
    
    <tr>
        <td>maxpoloidalpower</td>
        <td>Input</td>
        <td>real</td>
        <td>1000.0</td>
        <td><p>Maximum permitted absolute rate of change of stored energy in poloidal field (MW)</p></td>
    </tr>
    
    <tr>
        <td>poloidalpower</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Poloidal power usage at time t (MW)</p></td>
    </tr>
    
</table>

## pfcoil_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>nef</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>nfxf</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>ricpf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>ssq0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>sig_axial</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>sig_hoop</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>axial_force</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>rfxf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>zfxf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>cfxf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>xind</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>rcls</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>zcls</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>ccls</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>ccl0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>bpf2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>vsdum</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>first_call</td>
        <td>Input</td>
        <td>logical</td>
        <td>1</td>
        <td></td>
    </tr>
    
    <tr>
        <td>cslimit</td>
        <td>Output</td>
        <td>logical</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductorpf*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>croco_strand*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## pfcoil_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>ngrpmx</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>10</td>
        <td><p>maximum number of groups of PF coils</p></td>
    </tr>
    
    <tr>
        <td>nclsmx</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>2</td>
        <td><p>maximum number of PF coils in a given group</p></td>
    </tr>
    
    <tr>
        <td>nptsmx</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>32</td>
        <td><p>maximum number of points across the midplane of the plasma at which the field from
 the PF coils is fixed</p></td>
    </tr>
    
    <tr>
        <td>nfixmx</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>64</td>
        <td><p>maximum number of fixed current PF coils</p></td>
    </tr>
    
    <tr>
        <td>ngc</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>ngrpmx*nclsmx</td>
        <td><p>maximum total number of coils across all groups</p></td>
    </tr>
    
    <tr>
        <td>ngc2</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>ngc+2</td>
        <td><p>new variable to include 2 additional circuits: plasma and central solenoid</p></td>
    </tr>
    
    <tr>
        <td>alfapf</td>
        <td>Input</td>
        <td>real</td>
        <td>5e-10</td>
        <td><p>smoothing parameter used in PF coil current calculation at the beginning of pulse (BoP)</p></td>
    </tr>
    
    <tr>
        <td>alstroh</td>
        <td>Input</td>
        <td>real</td>
        <td>400000000.0</td>
        <td><p>allowable hoop stress in Central Solenoid structural material (Pa)</p></td>
    </tr>
    
    <tr>
        <td>i_cs_stress</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for CS stress calculation:</p>
<ul>
<li>=0 Hoop stress only</li>
<li>=1 Hoop + Axial stress</li>
</ul></td>
    </tr>
    
    <tr>
        <td>areaoh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Central solenoid vertical cross-sectional area (m2)</p></td>
    </tr>
    
    <tr>
        <td>a_oh_turn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Central solenoid (OH) trun cross-sectional area (m2)</p></td>
    </tr>
    
    <tr>
        <td>awpoh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central solenoid conductor+void area with area of steel subtracted (m2)</p></td>
    </tr>
    
    <tr>
        <td>bmaxoh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum field in central solenoid at end of flat-top (EoF) (T)</p></td>
    </tr>
    
    <tr>
        <td>bmaxoh0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum field in central solenoid at beginning of pulse (T)</p></td>
    </tr>
    
    <tr>
        <td>bpf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>peak field at coil i (T)</p></td>
    </tr>
    
    <tr>
        <td>ccl0_ma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>PF group current array, flux-swing cancellation current (MA)
 Input if i_pf_current=0, computed otherwise</p></td>
    </tr>
    
    <tr>
        <td>ccls_ma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>PF group current array, equilibrium current (MA)
 Input if i_pf_current=0, computed otherwise</p></td>
    </tr>
    
    <tr>
        <td>cohbop</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Central solenoid overall current density at beginning of pulse (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>coheof</td>
        <td>Input</td>
        <td>real</td>
        <td>18500000.0</td>
        <td><p>Central solenoid overall current density at end of flat-top (A/m2) (<code>iteration variable 37</code>) (<code>sweep variable 62</code>)</p></td>
    </tr>
    
    <tr>
        <td>cpt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>current per turn in coil i at time j (A)</p></td>
    </tr>
    
    <tr>
        <td>cptdin</td>
        <td>Input</td>
        <td>real</td>
        <td>[40000. 40000. 40000. 40000. 40000. 40000. 40000. 40000. 40000. 40000.
 40000. 40000. 40000. 40000. 40000. 40000. 40000. 40000. 40000. 40000.
 40000. 40000.]</td>
        <td><p>peak current per turn input for PF coil i (A)</p></td>
    </tr>
    
    <tr>
        <td>curpfb</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>PF coil current array, at beginning of pulse (MA)
 Indexed by coil number, not group number</p></td>
    </tr>
    
    <tr>
        <td>curpff</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>PF coil current array, at flat top (MA)
 Indexed by coil number, not group number</p></td>
    </tr>
    
    <tr>
        <td>curpfs</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>PF coil current array, at end of pulse (MA)
 Indexed by coil number, not group number</p></td>
    </tr>
    
    <tr>
        <td>etapsu</td>
        <td>Input</td>
        <td>real</td>
        <td>0.9</td>
        <td><p>Efficiency of transfer of PF stored energy into or out of storage.</p></td>
    </tr>
    
    <tr>
        <td>fcohbof</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ratio of central solenoid overall current density at beginning of flat-top / end of flat-top</p></td>
    </tr>
    
    <tr>
        <td>fcohbop</td>
        <td>Input</td>
        <td>real</td>
        <td>0.9</td>
        <td><p>ratio of central solenoid overall current density at beginning of pulse / end of flat-top
 (<code>iteration variable 41</code>)</p></td>
    </tr>
    
    <tr>
        <td>fcuohsu</td>
        <td>Input</td>
        <td>real</td>
        <td>0.7</td>
        <td><p>copper fraction of strand in central solenoid</p></td>
    </tr>
    
    <tr>
        <td>fcupfsu</td>
        <td>Input</td>
        <td>real</td>
        <td>0.69</td>
        <td><p>copper fraction of cable conductor (PF coils)</p></td>
    </tr>
    
    <tr>
        <td>fvssu</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>F-value for <code>constraint equation 51</code></p></td>
    </tr>
    
    <tr>
        <td>ipfloc</td>
        <td>Input</td>
        <td>integer</td>
        <td>[2 2 3 0 0 0 0 0 0 0]</td>
        <td><p>Switch for location of PF coil group i:</p>
<ul>
<li>=1 PF coil on top of central solenoid (flux ramp only)</li>
<li>=2 PF coil on top of TF coil (flux ramp only)</li>
<li>=3 PF coil outside of TF coil (equilibrium coil)</li>
<li>=4 PF coil, general location (equilibrium coil)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>ipfres</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for PF &amp; CS coil conductor type:</p>
<ul>
<li>=0 superconducting PF coils</li>
<li>=1 resistive PF coils</li>
</ul></td>
    </tr>
    
    <tr>
        <td>itr_sum</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total sum of I x turns x radius for all PF coils and CS (Am)</p></td>
    </tr>
    
    <tr>
        <td>isumatoh</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for superconductor material in central solenoid:</p>
<ul>
<li>=1 ITER Nb3Sn critical surface model with standard
   ITER parameters</li>
<li>=2 Bi-2212 high temperature superconductor (range of
   validity T &lt; 20K, adjusted field b &lt; 104 T, B &gt; 6 T)</li>
<li>=3 NbTi</li>
<li>=4 ITER Nb3Sn model with user-specified parameters</li>
<li>=5 WST Nb3Sn parameterisation</li>
<li>=6 REBCO HTS tape in CroCo strand</li>
<li>=7 Durham Ginzburg-Landau critical surface model for Nb-Ti</li>
<li>=8 Durham Ginzburg-Landau critical surface model for REBCO</li>
<li>=9 Hazelton experimental data + Zhai conceptual model for REBCO</li>
</ul></td>
    </tr>
    
    <tr>
        <td>isumatpf</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for superconductor material in PF coils:</p>
<ul>
<li>=1 ITER Nb3Sn critical surface model with standard
   ITER parameters</li>
<li>=2 Bi-2212 high temperature superconductor (range of
   validity T &lt; 20K, adjusted field b &lt; 104 T, B &gt; 6 T)</li>
<li>=3 NbTi</li>
<li>=4 ITER Nb3Sn model with user-specified parameters</li>
<li>=5 WST Nb3Sn parameterisation</li>
<li>=6 REBCO HTS tape in CroCo strand</li>
<li>=7 Durham Ginzburg-Landau critical surface model for Nb-Ti</li>
<li>=8 Durham Ginzburg-Landau critical surface model for REBCO</li>
<li>=9 Hazelton experimental data + Zhai conceptual model for REBCO</li>
</ul></td>
    </tr>
    
    <tr>
        <td>j_crit_str_cs</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>superconductor strand critical current density under operating
 conditions in central solenoid (A/m2). Necessary for the cost calculation in $/kA m</p></td>
    </tr>
    
    <tr>
        <td>j_crit_str_pf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>superconductor strand critical current density under operating
 conditions in PF coils (A/m2). Necessary for the cost calculation in $/kA m</p></td>
    </tr>
    
    <tr>
        <td>i_pf_current</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for controlling the current of the PF coils:</p>
<ul>
<li>=0 Input via the variables curpfb, curpff, curpfs</li>
<li>=1 SVD targets zero field across midplane (flux swing
   coils) and the correct vertical field at the plasma
   center (equilibrium coils)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_sup_pf_shape</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for the placement of Location 3 (outboard) PF coils
 when the TF coils are superconducting (i_tf_sup = 1)</p>
<ul>
<li>=0 (Default) Outboard PF coils follow TF shape
   in an ellipsoidal winding surface</li>
<li>=1 Outboard PF coils all have same radius, cylindrical
   winding surface</li>
</ul></td>
    </tr>
    
    <tr>
        <td>jscoh_bop</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central solenoid superconductor critical current density (A/m2) at beginning-of-pulse</p></td>
    </tr>
    
    <tr>
        <td>jscoh_eof</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central solenoid superconductor critical current density (A/m2) at end-of-flattop</p></td>
    </tr>
    
    <tr>
        <td>jcableoh_bop</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central solenoid cable critical current density (A/m2) at beginning-of-pulse</p></td>
    </tr>
    
    <tr>
        <td>jcableoh_eof</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central solenoid cable critical current density (A/m2) at end-of-flattop</p></td>
    </tr>
    
    <tr>
        <td>ncirt</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>number of PF circuits (including central solenoid and plasma)</p></td>
    </tr>
    
    <tr>
        <td>ncls</td>
        <td>Input</td>
        <td>integer</td>
        <td>[1 1 2 0 0 0 0 0 0 0 0 0]</td>
        <td><p>number of PF coils in group j</p></td>
    </tr>
    
    <tr>
        <td>nfxfh</td>
        <td>Input</td>
        <td>integer</td>
        <td>7</td>
        <td><p>number of filaments the top and bottom of the central solenoid should be broken
 into during scaling (5 - 10 is good)</p></td>
    </tr>
    
    <tr>
        <td>ngrp</td>
        <td>Input</td>
        <td>integer</td>
        <td>3</td>
        <td><p>number of groups of PF coils. Symmetric coil pairs should all be in the same group</p></td>
    </tr>
    
    <tr>
        <td>nohc</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>number of PF coils (excluding the central solenoid) + 1</p></td>
    </tr>
    
    <tr>
        <td>ohhghf</td>
        <td>Input</td>
        <td>real</td>
        <td>0.71</td>
        <td><p>Central solenoid height / TF coil internal height</p></td>
    </tr>
    
    <tr>
        <td>oh_steel_frac</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>central solenoid steel fraction (<code>iteration variable 122</code>)</p></td>
    </tr>
    
    <tr>
        <td>pf_current_safety_factor</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>Ratio of permissible PF coil conductor current density to critical conductor
 current density based on short-sample DC measurements</p></td>
    </tr>
    
    <tr>
        <td>pfcaseth</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>steel case thickness for PF coil i (m)</p></td>
    </tr>
    
    <tr>
        <td>pfclres</td>
        <td>Input</td>
        <td>real</td>
        <td>2.5e-08</td>
        <td><p>PF coil resistivity (if ipfres=1) (Ohm-m)</p></td>
    </tr>
    
    <tr>
        <td>rhopfbus</td>
        <td>Input</td>
        <td>real</td>
        <td>3.93e-08</td>
        <td><p>Resistivity of CS and PF coil bus bars (irrespective of
 whether the coils themselves are superconducting or resistive) (Ohm-m)</p></td>
    </tr>
    
    <tr>
        <td>pfmmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of heaviest PF coil (tonnes)</p></td>
    </tr>
    
    <tr>
        <td>pfrmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radius of largest PF coil (m)</p></td>
    </tr>
    
    <tr>
        <td>pfwpmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total mean wall plug power dissipated in PFC and CS power supplies (MW) (issue #713)</p></td>
    </tr>
    
    <tr>
        <td>powohres</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central solenoid resistive power during flattop (W)</p></td>
    </tr>
    
    <tr>
        <td>powpfres</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total PF coil resistive losses during flattop (W)</p></td>
    </tr>
    
    <tr>
        <td>ra</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>inner radius of coil i (m)</p></td>
    </tr>
    
    <tr>
        <td>rb</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>outer radius of coil i (m)</p></td>
    </tr>
    
    <tr>
        <td>ric</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>peak current in coil i (MA-turns)</p></td>
    </tr>
    
    <tr>
        <td>rjconpf</td>
        <td>Input</td>
        <td>real</td>
        <td>[30000000. 30000000. 30000000. 30000000. 30000000. 30000000. 30000000.
 30000000. 30000000. 30000000. 30000000. 30000000. 30000000. 30000000.
 30000000. 30000000. 30000000. 30000000. 30000000. 30000000. 30000000.
 30000000.]</td>
        <td><p>average winding pack current density of PF coil i (A/m2) at time of peak
 current in that coil (calculated for <code>ipfloc=1</code> coils)</p></td>
    </tr>
    
    <tr>
        <td>rjohc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>allowable central solenoid current density at end of flat-top (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>rjohc0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>allowable central solenoid current density at beginning of pulse (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>rjpfalw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>allowable winding pack current density of PF coil i (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>rohc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radius to the centre of the central solenoid (m)</p></td>
    </tr>
    
    <tr>
        <td>routr</td>
        <td>Input</td>
        <td>real</td>
        <td>1.5</td>
        <td><p>radial distance (m) from outboard TF coil leg to centre of <code>ipfloc=3</code> PF coils</p></td>
    </tr>
    
    <tr>
        <td>rpf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radius of PF coil i (m)</p></td>
    </tr>
    
    <tr>
        <td>rpf1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>offset (m) of radial position of <code>ipfloc=1</code> PF coils from being directly above
 the central solenoid</p></td>
    </tr>
    
    <tr>
        <td>rpf2</td>
        <td>Input</td>
        <td>real</td>
        <td>-1.63</td>
        <td><p>offset (m) of radial position of <code>ipfloc=2</code> PF coils from being at
 rmajor (offset = rpf2<em>triang</em>rminor)</p></td>
    </tr>
    
    <tr>
        <td>rref</td>
        <td>Input</td>
        <td>real</td>
        <td>[7. 7. 7. 7. 7. 7. 7. 7. 7. 7.]</td>
        <td><p>PF coil radial positioning adjuster:</p>
<ul>
<li>for groups j with ipfloc(j) = 1; rref(j) is ignored</li>
<li>for groups j with ipfloc(j) = 2; rref(j) is ignored</li>
<li>for groups j with ipfloc(j) = 3; rref(j) is ignored</li>
<li>for groups j with ipfloc(j) = 4; rref(j) is radius of
   the coil in units of minor radii from the major radius
   (r = rmajor + rref*rminor)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>s_tresca_oh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Maximum shear stress (Tresca criterion) coils/central solenoid [MPa]</p></td>
    </tr>
    
    <tr>
        <td>sigpfcalw</td>
        <td>Input</td>
        <td>real</td>
        <td>500.0</td>
        <td><p>maximum permissible tensile stress (MPa) in steel coil cases for superconducting
 PF coils (<code>ipfres=0</code>)</p></td>
    </tr>
    
    <tr>
        <td>sigpfcf</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>fraction of JxB hoop force supported by steel case for superconducting PF coils (<code>ipfres=0</code>)</p></td>
    </tr>
    
    <tr>
        <td>sxlg</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mutual inductance matrix (H)</p></td>
    </tr>
    
    <tr>
        <td>tmargoh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Central solenoid temperature margin (K)</p></td>
    </tr>
    
    <tr>
        <td>turns</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>number of turns in PF coil i</p></td>
    </tr>
    
    <tr>
        <td>vf</td>
        <td>Input</td>
        <td>real</td>
        <td>[0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3
 0.3 0.3 0.3 0.3]</td>
        <td><p>winding pack void fraction of PF coil i for coolant</p></td>
    </tr>
    
    <tr>
        <td>vfohc</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>void fraction of central solenoid conductor for coolant</p></td>
    </tr>
    
    <tr>
        <td>vsbn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total flux swing available for burn (Wb)</p></td>
    </tr>
    
    <tr>
        <td>vsefbn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>flux swing from PF coils for burn (Wb)</p></td>
    </tr>
    
    <tr>
        <td>vsefsu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>flux swing from PF coils for startup (Wb)</p></td>
    </tr>
    
    <tr>
        <td>vseft</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total flux swing from PF coils (Wb)</p></td>
    </tr>
    
    <tr>
        <td>vsoh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total flux swing from the central solenoid (Wb)</p></td>
    </tr>
    
    <tr>
        <td>vsohbn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central solenoid flux swing for burn (Wb)</p></td>
    </tr>
    
    <tr>
        <td>vsohsu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central solenoid flux swing for startup (Wb)</p></td>
    </tr>
    
    <tr>
        <td>vssu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total flux swing for startup (<code>constraint eqn 51</code> to enforce vssu=vs_plasma_res_ramp+vs_plasma_ind_ramp) (Wb)</p></td>
    </tr>
    
    <tr>
        <td>vstot</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total flux swing for pulse (Wb)</p></td>
    </tr>
    
    <tr>
        <td>waves</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>used in current waveform of PF coils/central solenoid</p></td>
    </tr>
    
    <tr>
        <td>whtpf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total mass of the PF coil conductor (kg)</p></td>
    </tr>
    
    <tr>
        <td>whtpfs</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total mass of the PF coil structure (kg)</p></td>
    </tr>
    
    <tr>
        <td>wtc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>conductor mass for PF coil i (kg)</p></td>
    </tr>
    
    <tr>
        <td>wts</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>structure mass for PF coil i (kg)</p></td>
    </tr>
    
    <tr>
        <td>zh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>upper point of PF coil i (m)</p></td>
    </tr>
    
    <tr>
        <td>zl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>lower point of PF coil i (m)</p></td>
    </tr>
    
    <tr>
        <td>zpf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>z (height) location of PF coil i (m)</p></td>
    </tr>
    
    <tr>
        <td>zref</td>
        <td>Input</td>
        <td>real</td>
        <td>[3.6 1.2 2.5 1.  1.  1.  1.  1.  1.  1. ]</td>
        <td><p>PF coil vertical positioning adjuster:</p>
<ul>
<li>for groups j with ipfloc(j) = 1; zref(j) is ignored</li>
<li>for groups j with ipfloc(j) = 2 AND itart=1 (only);
   zref(j) is distance of centre of PF coil from inside
   edge of TF coil (remember that PF coils for STs lie
   within the TF coil)</li>
<li>for groups j with ipfloc(j) = 3; zref(j) = ratio of
   height of coil group j to plasma minor radius</UL></li>
<li>for groups j with ipfloc(j) = 4; zref(j) = ratio of
   height of coil group j to plasma minor radius</UL></li>
</ul></td>
    </tr>
    
    <tr>
        <td>bmaxcs_lim</td>
        <td>Input</td>
        <td>real</td>
        <td>13.0</td>
        <td><p>Central solenoid max field limit [T]</p></td>
    </tr>
    
    <tr>
        <td>fbmaxcs</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>F-value for CS mmax field (<code>cons. 79</code>, <code>itvar 149</code>)</p></td>
    </tr>
    
    <tr>
        <td>ld_ratio_cst</td>
        <td>Input</td>
        <td>real</td>
        <td>3.0</td>
        <td><p>Ratio of CS coil turn conduit length to depth</p></td>
    </tr>
    
    <tr>
        <td>l_cond_cst</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Length of CS of CS coil turn conduit</p></td>
    </tr>
    
    <tr>
        <td>d_cond_cst</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Depth/width of CS of CS coil turn conduit</p></td>
    </tr>
    
    <tr>
        <td>r_out_cst</td>
        <td>Input</td>
        <td>real</td>
        <td>0.003</td>
        <td><p>Length of CS of CS coil turn conduit length</p></td>
    </tr>
    
    <tr>
        <td>r_in_cst</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Length of CS of CS coil turn conduit length</p></td>
    </tr>
    
</table>

## physics_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>iscz</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>err242</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>err243</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>rad_fraction_lcfs</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>e_plasma_beta</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>total_loss_power</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>t_energy_confinement_beta</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>ptarmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>lambdaio</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>drsep</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>fio</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>fli</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>flo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>fui</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>fuo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>plimw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>plomw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>puimw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>puomw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>rho_star</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>nu_star</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>beta_mcdonald</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>itart_r</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>first_call</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td></td>
    </tr>
    
</table>

## physics_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>n_confinement_scalings</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>51</td>
        <td><p>number of energy confinement time scaling laws</p></td>
    </tr>
    
    <tr>
        <td>m_beam_amu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>beam ion mass (amu)</p></td>
    </tr>
    
    <tr>
        <td>m_fuel_amu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>average mass of fuel portion of ions (amu)</p></td>
    </tr>
    
    <tr>
        <td>m_ions_total_amu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>average mass of all ions (amu)</p></td>
    </tr>
    
    <tr>
        <td>alphaj</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>current profile index (calculated from q_0 and q if <code>iprofile=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>alphan</td>
        <td>Input</td>
        <td>real</td>
        <td>0.25</td>
        <td><p>density profile index</p></td>
    </tr>
    
    <tr>
        <td>alphap</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>pressure profile index</p></td>
    </tr>
    
    <tr>
        <td>alpha_rate_density_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Alpha particle production rate per unit volume, from plasma and beams [particles/m3/sec]</p></td>
    </tr>
    
    <tr>
        <td>alpha_rate_density_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Alpha particle production rate per unit volume, just from plasma [particles/m3/sec]</p></td>
    </tr>
    
    <tr>
        <td>alphat</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>temperature profile index</p></td>
    </tr>
    
    <tr>
        <td>aspect</td>
        <td>Input</td>
        <td>real</td>
        <td>2.907</td>
        <td><p>aspect ratio (<code>iteration variable 1</code>)</p></td>
    </tr>
    
    <tr>
        <td>beamfus0</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>multiplier for beam-background fusion calculation</p></td>
    </tr>
    
    <tr>
        <td>beta</td>
        <td>Input</td>
        <td>real</td>
        <td>0.042</td>
        <td><p>total plasma beta (<code>iteration variable 5</code>) (calculated if stellarator)</p></td>
    </tr>
    
    <tr>
        <td>beta_fast_alpha</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fast alpha beta component</p></td>
    </tr>
    
    <tr>
        <td>beta_max</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Max allowable beta</p></td>
    </tr>
    
    <tr>
        <td>beta_min</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>allowable lower beta</p></td>
    </tr>
    
    <tr>
        <td>beta_beam</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutral beam beta component</p></td>
    </tr>
    
    <tr>
        <td>beta_poloidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>poloidal beta</p></td>
    </tr>
    
    <tr>
        <td>beta_poloidal_eps</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Poloidal beta and inverse aspcet ratio product</p></td>
    </tr>
    
    <tr>
        <td>beta_toroidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>toroidal beta</p></td>
    </tr>
    
    <tr>
        <td>beta_thermal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>thermal beta</p></td>
    </tr>
    
    <tr>
        <td>beta_thermal_poloidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>poloidal thermal beta</p></td>
    </tr>
    
    <tr>
        <td>beta_thermal_toroidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>poloidal thermal beta</p></td>
    </tr>
    
    <tr>
        <td>beta_norm_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>normaised total beta</p></td>
    </tr>
    
    <tr>
        <td>beta_norm_thermal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>normaised thermal beta</p></td>
    </tr>
    
    <tr>
        <td>beta_norm_toroidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>normaised toroidal beta</p></td>
    </tr>
    
    <tr>
        <td>beta_norm_poloidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>normaised poloidal beta</p></td>
    </tr>
    
    <tr>
        <td>e_plasma_beta_thermal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Plasma thermal energy derived from thermal beta</p></td>
    </tr>
    
    <tr>
        <td>betbm0</td>
        <td>Input</td>
        <td>real</td>
        <td>1.5</td>
        <td><p>leading coefficient for NB beta fraction</p></td>
    </tr>
    
    <tr>
        <td>bp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>poloidal field (T)</p></td>
    </tr>
    
    <tr>
        <td>bt</td>
        <td>Input</td>
        <td>real</td>
        <td>5.68</td>
        <td><p>toroidal field on axis (T) (<code>iteration variable 2</code>)</p></td>
    </tr>
    
    <tr>
        <td>btot</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total toroidal + poloidal field (T)</p></td>
    </tr>
    
    <tr>
        <td>burnup</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fractional plasma burnup</p></td>
    </tr>
    
    <tr>
        <td>burnup_in</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fractional plasma burnup user input</p></td>
    </tr>
    
    <tr>
        <td>bvert</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vertical field at plasma (T)</p></td>
    </tr>
    
    <tr>
        <td>c_beta</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>Destabalisation parameter for iprofile=6 beta limit</p></td>
    </tr>
    
    <tr>
        <td>csawth</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>coeff. for sawteeth effects on burn V-s requirement</p></td>
    </tr>
    
    <tr>
        <td>f_vol_plasma</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>multiplying factor for the plasma volume (normally=1)</p></td>
    </tr>
    
    <tr>
        <td>f_r_conducting_wall</td>
        <td>Input</td>
        <td>real</td>
        <td>1.35</td>
        <td><p>maximum ratio of conducting wall distance to plasma minor radius for
 vertical stability (<code>constraint equation 23</code>)</p></td>
    </tr>
    
    <tr>
        <td>dene</td>
        <td>Input</td>
        <td>real</td>
        <td>9.8e+19</td>
        <td><p>electron density (/m3) (<code>iteration variable 6</code>)</p></td>
    </tr>
    
    <tr>
        <td>nd_fuel_ions</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fuel ion density (/m3)</p></td>
    </tr>
    
    <tr>
        <td>dlamee</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>electron-electron coulomb logarithm</p></td>
    </tr>
    
    <tr>
        <td>dlamie</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ion-electron coulomb logarithm</p></td>
    </tr>
    
    <tr>
        <td>dlimit</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>density limit (/m3) as calculated using various models</p></td>
    </tr>
    
    <tr>
        <td>nd_alphas</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>thermal alpha density (/m3)</p></td>
    </tr>
    
    <tr>
        <td>nd_beam_ions</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>hot beam ion density, variable (/m3)</p></td>
    </tr>
    
    <tr>
        <td>beam_density_out</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>hot beam ion density from calculation (/m3)</p></td>
    </tr>
    
    <tr>
        <td>beta_norm_max</td>
        <td>Input</td>
        <td>real</td>
        <td>3.5</td>
        <td><p>Troyon-like coefficient for beta scaling</p></td>
    </tr>
    
    <tr>
        <td>dnelimt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>density limit (/m3)</p></td>
    </tr>
    
    <tr>
        <td>nd_ions_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total ion density (/m3)</p></td>
    </tr>
    
    <tr>
        <td>dnla</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>line averaged electron density (/m3)</p></td>
    </tr>
    
    <tr>
        <td>nd_protons</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>proton ash density (/m3)</p></td>
    </tr>
    
    <tr>
        <td>ntau</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Fusion double product (s/m3)</p></td>
    </tr>
    
    <tr>
        <td>nttau</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Lawson triple product [keV s / m3]</p></td>
    </tr>
    
    <tr>
        <td>nd_impurities</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>high Z ion density (/m3)</p></td>
    </tr>
    
    <tr>
        <td>gradient_length_ne</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Max. normalized gradient length in el. density (ipedestal==0 only)</p></td>
    </tr>
    
    <tr>
        <td>gradient_length_te</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Max. normalized gradient length in el. temperature (ipedestal==0 only)</p></td>
    </tr>
    
    <tr>
        <td>beta_poloidal_eps_max</td>
        <td>Input</td>
        <td>real</td>
        <td>1.38</td>
        <td><p>maximum (eps*beta_poloidal) (<code>constraint equation 6</code>). Note: revised issue #346
 "Operation at the tokamak equilibrium poloidal beta-limit in TFTR", 1992 Nucl. Fusion 32 1468</p></td>
    </tr>
    
    <tr>
        <td>eps</td>
        <td>Input</td>
        <td>real</td>
        <td>0.34399724802</td>
        <td><p>inverse aspect ratio</p></td>
    </tr>
    
    <tr>
        <td>aux_current_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fraction of plasma current produced by auxiliary current drive</p></td>
    </tr>
    
    <tr>
        <td>inductive_current_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fraction of plasma current produced inductively</p></td>
    </tr>
    
    <tr>
        <td>f_alpha_electron</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fraction of alpha energy to electrons</p></td>
    </tr>
    
    <tr>
        <td>f_alpha_plasma</td>
        <td>Input</td>
        <td>real</td>
        <td>0.95</td>
        <td><p>Fraction of alpha power deposited in plasma. Default of 0.95 taken from https://doi.org/10.1088/0029-5515/39/12/305.</p></td>
    </tr>
    
    <tr>
        <td>f_alpha_ion</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fraction of alpha power to ions</p></td>
    </tr>
    
    <tr>
        <td>f_deuterium</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>deuterium fuel fraction</p></td>
    </tr>
    
    <tr>
        <td>ftar</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>fraction of power to the lower divertor in double null configuration
 (<code>i_single_null = 0</code> only) (default assumes SN)</p></td>
    </tr>
    
    <tr>
        <td>ffwal</td>
        <td>Input</td>
        <td>real</td>
        <td>0.92</td>
        <td><p>factor to convert plasma surface area to first wall area in neutron wall
 load calculation (<code>iwalld=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>fgwped</td>
        <td>Input</td>
        <td>real</td>
        <td>0.85</td>
        <td><p>fraction of Greenwald density to set as pedestal-top density. If <code>&lt;0</code>, pedestal-top
 density set manually using neped (<code>ipedestal==1</code>).
 (<code>iteration variable 145</code>)</p></td>
    </tr>
    
    <tr>
        <td>fgwsep</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>fraction of Greenwald density to set as separatrix density. If <code>&lt;0</code>, separatrix
 density set manually using nesep (<code>ipedestal==1</code>).
 (<code>iteration variable 152</code>)</p></td>
    </tr>
    
    <tr>
        <td>f_helium3</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>helium-3 fuel fraction</p></td>
    </tr>
    
    <tr>
        <td>figmer</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>physics figure of merit (= plasma_current<em>aspect</em>*sbar, where <code>sbar=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>fkzohm</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>Zohm elongation scaling adjustment factor (<code>i_plasma_geometry=2, 3</code>)</p></td>
    </tr>
    
    <tr>
        <td>fplhsep</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>F-value for Psep &gt;= Plh + Paux (<code>constraint equation 73</code>)</p></td>
    </tr>
    
    <tr>
        <td>fpdivlim</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>F-value for minimum pdivt (<code>constraint equation 80</code>)</p></td>
    </tr>
    
    <tr>
        <td>fne0</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for the constraint ne(0) &gt; ne(ped) (<code>constraint equation 81</code>)
 (<code>Iteration variable 154</code>)</p></td>
    </tr>
    
    <tr>
        <td>f_tritium</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>tritium fuel fraction</p></td>
    </tr>
    
    <tr>
        <td>fusion_rate_density_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fusion reaction rate, from beams and plasma (reactions/m3/sec)</p></td>
    </tr>
    
    <tr>
        <td>fusion_rate_density_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fusion reaction rate, just from plasma (reactions/m3/sec)</p></td>
    </tr>
    
    <tr>
        <td>fvsbrnni</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>fraction of the plasma current produced by non-inductive means (<code>iteration variable 44</code>)</p></td>
    </tr>
    
    <tr>
        <td>ejima_coeff</td>
        <td>Input</td>
        <td>real</td>
        <td>0.4</td>
        <td><p>Ejima coefficient for resistive startup V-s formula</p></td>
    </tr>
    
    <tr>
        <td>f_beta_alpha_beam_thermal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ratio of (fast alpha + neutral beam beta) to thermal beta</p></td>
    </tr>
    
    <tr>
        <td>hfac</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>H factors for an ignited plasma for each energy confinement time scaling law</p></td>
    </tr>
    
    <tr>
        <td>hfact</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>H factor on energy confinement times, radiation corrected (<code>iteration variable 10</code>).</p></td>
    </tr>
    
    <tr>
        <td>taumax</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>Maximum allowed energy confinement time (s)</p></td>
    </tr>
    
    <tr>
        <td>i_bootstrap_current</td>
        <td>Input</td>
        <td>integer</td>
        <td>3</td>
        <td><p>switch for bootstrap current scaling</p>
<ul>
<li>=1 ITER 1989 bootstrap scaling (high R/a only)</li>
<li>=2 for Nevins et al general scaling</li>
<li>=3 for Wilson et al numerical scaling</li>
<li>=4 for Sauter et al scaling</li>
<li>=5 for Sakai et al scaling</li>
<li>=6 for ARIES scaling</li>
<li>=7 for Andrade et al scaling</li>
<li>=8 for Hoang et al scaling</li>
<li>=9 for Wong et al scaling</li>
<li>=10 for Gi-I et al scaling</li>
<li>=11 for Gi-II et al scaling</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_beta_component</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for beta limit scaling (<code>constraint equation 24</code>)</p>
<ul>
<li>=0 apply limit to total beta</li>
<li>=1 apply limit to thermal beta</li>
<li>=2 apply limit to thermal + neutral beam beta</li>
<li>=3 apply limit to toroidal beta</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_plasma_current</td>
        <td>Input</td>
        <td>integer</td>
        <td>4</td>
        <td><p>switch for plasma current scaling to use</p>
<ul>
<li>=1 Peng analytic fit</li>
<li>=2 Peng double null divertor scaling (ST)</li>
<li>=3 simple ITER scaling (k = 2.2, d = 0.6)</li>
<li>=4 later ITER scaling, a la Uckan</li>
<li>=5 Todd empirical scaling I</li>
<li>=6 Todd empirical scaling II</li>
<li>=7 Connor-Hastie model</li>
<li>=8 Sauter scaling allowing negative triangularity</li>
<li>=9 FIESTA ST fit</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_diamagnetic_current</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for diamagnetic current scaling</p>
<ul>
<li>=0 Do not calculate</li>
<li>=1 Use original TART scaling</li>
<li>=2 Use SCENE scaling</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_density_limit</td>
        <td>Input</td>
        <td>integer</td>
        <td>8</td>
        <td><p>switch for density limit to enforce (<code>constraint equation 5</code>)</p>
<ul>
<li>=1 old ASDEX</li>
<li>=2 Borrass model for ITER (I)</li>
<li>=3 Borrass model for ITER (II)</li>
<li>=4 JET edge radiation</li>
<li>=5 JET simplified</li>
<li>=6 Hugill-Murakami Mq limit</li>
<li>=7 Greenwald limit</li>
<li>=8 ASDEX New</li>
</ul></td>
    </tr>
    
    <tr>
        <td>idivrt</td>
        <td>Input</td>
        <td>integer</td>
        <td>2</td>
        <td><p>number of divertors (calculated from <code>i_single_null</code>)</p></td>
    </tr>
    
    <tr>
        <td>i_beta_fast_alpha</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for fast alpha pressure calculation</p>
<ul>
<li>=0 ITER physics rules (Uckan) fit</li>
<li>=1 Modified fit (D. Ward) - better at high temperature</li>
</ul></td>
    </tr>
    
    <tr>
        <td>ignite</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for ignition assumption. Obviously, ignite must be zero if current drive
 is required. If ignite is 1, any auxiliary power is assumed to be used only during
 plasma start-up, and is excluded from all steady-state power balance calculations.</p>
<ul>
<li>=0 do not assume plasma ignition</li>
<li>=1 assume ignited (but include auxiliary power in costs)&lt;/UL</li>
</ul></td>
    </tr>
    
    <tr>
        <td>ipedestal</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for pedestal profiles:</p>
<ul>
<li>=0 use original parabolic profiles</li>
<li>=1 use pedestal profile</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_pfirsch_schluter_current</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for Pfirsch-Schlüter current scaling (issue #413):</p>
<ul>
<li>=0 Do not calculate</li>
<li>=1 Use SCENE scaling</li>
</ul></td>
    </tr>
    
    <tr>
        <td>neped</td>
        <td>Input</td>
        <td>real</td>
        <td>4e+19</td>
        <td><p>electron density of pedestal [m-3] (`ipedestal==1)</p></td>
    </tr>
    
    <tr>
        <td>nesep</td>
        <td>Input</td>
        <td>real</td>
        <td>3e+19</td>
        <td><p>electron density at separatrix [m-3] (`ipedestal==1)</p></td>
    </tr>
    
    <tr>
        <td>alpha_crit</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>critical ballooning parameter value</p></td>
    </tr>
    
    <tr>
        <td>nesep_crit</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>critical electron density at separatrix [m-3]</p></td>
    </tr>
    
    <tr>
        <td>plasma_res_factor</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>plasma resistivity pre-factor</p></td>
    </tr>
    
    <tr>
        <td>rhopedn</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>r/a of density pedestal (<code>ipedestal==1</code>)</p></td>
    </tr>
    
    <tr>
        <td>rhopedt</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>r/a of temperature pedestal (<code>ipedestal==1</code>)</p></td>
    </tr>
    
    <tr>
        <td>rho_te_max</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>r/a where the temperature gradient is largest (<code>ipedestal==0</code>)</p></td>
    </tr>
    
    <tr>
        <td>rho_ne_max</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>r/a where the density gradient is largest (<code>ipedestal==0</code>)</p></td>
    </tr>
    
    <tr>
        <td>tbeta</td>
        <td>Input</td>
        <td>real</td>
        <td>2.0</td>
        <td><p>temperature profile index beta  (`ipedestal==1)</p></td>
    </tr>
    
    <tr>
        <td>teped</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>electron temperature of pedestal (keV) (<code>ipedestal==1</code>)</p></td>
    </tr>
    
    <tr>
        <td>tesep</td>
        <td>Input</td>
        <td>real</td>
        <td>0.1</td>
        <td><p>electron temperature at separatrix (keV) (<code>ipedestal==1</code>) calculated if reinke
 criterion is used (<code>icc=78</code>)</p></td>
    </tr>
    
    <tr>
        <td>iprofile</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for current profile consistency:</p>
<ul>
<li>=0 use input values for alphaj, ind_plasma_internal_norm, beta_norm_max</li>
<li>=1 make these consistent with input q95, q_0 values (recommend <code>i_plasma_current=4</code> with this option)</li>
<li>=2 use input values for alphaj, ind_plasma_internal_norm. Scale beta_norm_max with aspect ratio (original scaling)</li>
<li>=3 use input values for alphaj, ind_plasma_internal_norm. Scale beta_norm_max with aspect ratio (Menard scaling)</li>
<li>=4 use input values for alphaj, beta_norm_max. Set ind_plasma_internal_norm from elongation (Menard scaling)</li>
<li>=5 use input value for alphaj.  Set ind_plasma_internal_norm and beta_norm_max from Menard scaling</li>
<li>=6 use input values for alphaj, c_beta.  Set ind_plasma_internal_norm from Menard and beta_norm_max from Tholerus</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_rad_loss</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for radiation loss term usage in power balance (see User Guide):</p>
<ul>
<li>=0 total power lost is scaling power plus radiation</li>
<li>=1 total power lost is scaling power plus core radiation only</li>
<li>=2 total power lost is scaling power only, with no additional
   allowance for radiation. This is not recommended for power plant models.</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_confinement_time</td>
        <td>Input</td>
        <td>integer</td>
        <td>34</td>
        <td><p>switch for energy confinement time scaling law (see description in <code>labels_confinement_scalings</code>)</p>
<p>labels_confinement_scalings(n_confinement_scalings) : labels describing energy confinement scaling laws</p></td>
    </tr>
    
    <tr>
        <td>labels_confinement_scalings</td>
        <td>Parameter</td>
        <td>character</td>
        <td>(/'User input electron confinement&nbsp;&nbsp;&nbsp;', 'Neo-Alcator&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Ohmic)', 'Mirnov&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'Merezkhin-Muhkovatov&nbsp;&nbsp;&nbsp;&nbsp;(Ohmic)(L)', 'Shimomura&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'Kaye-Goldston&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'ITER 89-P&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'ITER 89-O&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'Rebut-Lallia&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'Goldston&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'T10&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'JAERI / Odajima-Shimomura&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'Kaye-Big Complex&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'ITER H90-P&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'ITER 89-P & 89-O min&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'Riedel&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'Christiansen&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'Lackner-Gottardi&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'Neo-Kaye&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'Riedel&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'ITER H90-P amended&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'LHD&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Stell)', 'Gyro-reduced Bohm&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Stell)', 'Lackner-Gottardi&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Stell)', 'ITER-93H&nbsp;&nbsp;ELM-free&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'TITAN RFP OBSOLETE&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;', 'ITER H-97P ELM-free&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'ITER H-97P ELMy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'ITER-96P (ITER-97L)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'Valovic modified ELMy&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'Kaye 98 modified&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(L)', 'ITERH-PB98P(y)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'IPB98(y)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'IPB98(y,1)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'IPB98(y,2)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'IPB98(y,3)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'IPB98(y,4)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'ISS95&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Stell)', 'ISS04&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(Stell)', 'DS03 beta-independent&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'Murari "Non-power law"&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'Petty 2008&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(ST)(H)', 'Lang high density&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'Hubbard 2017 - nominal&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(I)', 'Hubbard 2017 - lower&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(I)', 'Hubbard 2017 - upper&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(I)', 'Menard NSTX&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(ST)(H)', 'Menard NSTX-Petty08 hybrid (ST)(H)', 'Buxton NSTX gyro-Bohm&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(ST)(H)', 'ITPA20&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)', 'ITPA20-IL&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(H)'/)</td>
        <td></td>
    </tr>
    
    <tr>
        <td>i_plasma_wall_gap</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for plasma-first wall clearances at the mid-plane:</p>
<ul>
<li>=0 use 10% of plasma minor radius</li>
<li>=1 use input (<code>dr_fw_plasma_gap_inboard</code> and <code>dr_fw_plasma_gap_outboard</code>)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_plasma_geometry</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for plasma elongation and triangularity calculations:</p>
<ul>
<li>=0 use input kappa, triang to calculate 95% values</li>
<li>=1 scale q95_min, kappa, triang with aspect ratio (ST)</li>
<li>=2 set kappa to the natural elongation value (Zohm ITER scaling), triang input</li>
<li>=3 set kappa to the natural elongation value (Zohm ITER scaling), triang95 input</li>
<li>=4 use input kappa95, triang95 to calculate separatrix values</li>
<li>=5 use input kappa95, triang95 to calculate separatrix values based on MAST scaling (ST)</li>
<li>=6 use input kappa, triang to calculate 95% values based on MAST scaling (ST)</li>
<li>=7 use input kappa95, triang95 to calculate separatrix values based on fit to FIESTA (ST)</li>
<li>=8 use input kappa, triang to calculate 95% values based on fit to FIESTA (ST)</li>
<li>=9 set kappa to the natural elongation value, triang input</li>
<li>=10 set kappa to maximum stable value at a given aspect ratio (2.6&lt;A&lt;3.6)), triang input (#1399)</li>
<li>=11 set kappa Menard 2016 aspect-ratio-dependent scaling, triang input (#1439)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_plasma_shape</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for plasma boundary shape:</p>
<ul>
<li>=0 use original PROCESS 2-arcs model</li>
<li>=1 use the Sauter model</li>
</ul></td>
    </tr>
    
    <tr>
        <td>itart</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for spherical tokamak (ST) models:</p>
<ul>
<li>=0 use conventional aspect ratio models</li>
<li>=1 use spherical tokamak models</li>
</ul></td>
    </tr>
    
    <tr>
        <td>itartpf</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for Spherical Tokamak PF models:</p>
<ul>
<li>=0 use Peng and Strickler (1986) model</li>
<li>=1 use conventional aspect ratio model</li>
</ul></td>
    </tr>
    
    <tr>
        <td>iwalld</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for neutron wall load calculation:</p>
<ul>
<li>=1 use scaled plasma surface area</li>
<li>=2 use first wall area directly</li>
</ul></td>
    </tr>
    
    <tr>
        <td>plasma_square</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma squareness used by Sauter plasma shape</p></td>
    </tr>
    
    <tr>
        <td>kappa</td>
        <td>Input</td>
        <td>real</td>
        <td>1.792</td>
        <td><p>plasma separatrix elongation (calculated if <code>i_plasma_geometry = 1-5, 7 or 9-10</code>)</p></td>
    </tr>
    
    <tr>
        <td>kappa95</td>
        <td>Input</td>
        <td>real</td>
        <td>1.6</td>
        <td><p>plasma elongation at 95% surface (calculated if <code>i_plasma_geometry = 0-3, 6, or 8-10</code>)</p></td>
    </tr>
    
    <tr>
        <td>kappa_ipb</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Separatrix elongation calculated for IPB scalings</p></td>
    </tr>
    
    <tr>
        <td>ne0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central electron density (/m3)</p></td>
    </tr>
    
    <tr>
        <td>ni0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central ion density (/m3)</p></td>
    </tr>
    
    <tr>
        <td>m_s_limit</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>margin to vertical stability</p></td>
    </tr>
    
    <tr>
        <td>p0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central total plasma pressure (Pa)</p></td>
    </tr>
    
    <tr>
        <td>j_plasma_0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Central plasma current density (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>vol_avg_pressure</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Volume averaged plasma pressure (Pa)</p></td>
    </tr>
    
    <tr>
        <td>f_dd_branching_trit</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>branching ratio for DD -&gt; T</p></td>
    </tr>
    
    <tr>
        <td>alpha_power_density_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Alpha power per volume just from plasma [MW/m3]</p></td>
    </tr>
    
    <tr>
        <td>alpha_power_density_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Alpha power per volume from plasma and beams [MW/m3]</p></td>
    </tr>
    
    <tr>
        <td>alpha_power_electron_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Alpha power per volume to electrons [MW/m3]</p></td>
    </tr>
    
    <tr>
        <td>palpfwmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>alpha power escaping plasma and reaching first wall (MW)</p></td>
    </tr>
    
    <tr>
        <td>alpha_power_ions_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>alpha power per volume to ions (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>alpha_power_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Alpha power from only the plasma (MW)</p></td>
    </tr>
    
    <tr>
        <td>alpha_power_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total alpha power from plasma and beams (MW)</p></td>
    </tr>
    
    <tr>
        <td>alpha_power_beams</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>alpha power from hot neutral beam ions (MW)</p></td>
    </tr>
    
    <tr>
        <td>non_alpha_charged_power</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>non-alpha charged particle fusion power (MW)</p></td>
    </tr>
    
    <tr>
        <td>charged_particle_power</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total charged particle fusion power [MW]</p></td>
    </tr>
    
    <tr>
        <td>charged_power_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Non-alpha charged particle fusion power per volume [MW/m3]</p></td>
    </tr>
    
    <tr>
        <td>pcoef</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>profile factor (= n-weighted T / average T)</p></td>
    </tr>
    
    <tr>
        <td>p_plasma_inner_rad_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radiation power from inner zone (MW)</p></td>
    </tr>
    
    <tr>
        <td>pcoreradpv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total core radiation power per volume (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>dd_power</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>deuterium-deuterium fusion power (MW)</p></td>
    </tr>
    
    <tr>
        <td>dhe3_power</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>deuterium-helium3 fusion power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pdivt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>power to conducted to the divertor region (MW)</p></td>
    </tr>
    
    <tr>
        <td>pdivl</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>power conducted to the lower divertor region (calculated if <code>i_single_null = 0</code>) (MW)</p></td>
    </tr>
    
    <tr>
        <td>pdivu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>power conducted to the upper divertor region (calculated if <code>i_single_null = 0</code>) (MW)</p></td>
    </tr>
    
    <tr>
        <td>pdivmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>power conducted to the divertor with most load (calculated if <code>i_single_null = 0</code>) (MW)</p></td>
    </tr>
    
    <tr>
        <td>dt_power_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total deuterium-tritium fusion power, from plasma and beams [MW]</p></td>
    </tr>
    
    <tr>
        <td>dt_power_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Deuterium-tritium fusion power, just from plasma [MW]</p></td>
    </tr>
    
    <tr>
        <td>p_plasma_outer_rad_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radiation power from outer zone (MW)</p></td>
    </tr>
    
    <tr>
        <td>pedgeradpv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>edge radiation power per volume (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>vs_plasma_internal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>internal plasma V-s</p></td>
    </tr>
    
    <tr>
        <td>pflux_fw_rad_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Nominal mean radiation load on inside surface of reactor (MW/m2)</p></td>
    </tr>
    
    <tr>
        <td>piepv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ion/electron equilibration power per volume (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>plasma_current</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma current (A)</p></td>
    </tr>
    
    <tr>
        <td>neutron_power_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Neutron fusion power from just the plasma [MW]</p></td>
    </tr>
    
    <tr>
        <td>neutron_power_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total neutron fusion power from plasma and beams [MW]</p></td>
    </tr>
    
    <tr>
        <td>neutron_power_density_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutron fusion power per volume from beams and plasma (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>neutron_power_density_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neutron fusion power per volume just from plasma (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>p_plasma_ohmic_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ohmic heating power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pden_plasma_ohmic_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ohmic heating power per volume (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>p_plasma_loss_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>heating power (= transport loss power) (MW) used in confinement time calculation</p></td>
    </tr>
    
    <tr>
        <td>fusion_power</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fusion power (MW)</p></td>
    </tr>
    
    <tr>
        <td>len_plasma_poloidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma poloidal perimeter (m)</p></td>
    </tr>
    
    <tr>
        <td>p_plasma_rad_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total radiation power from inside LCFS (MW)</p></td>
    </tr>
    
    <tr>
        <td>pden_plasma_rad_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total radiation power per volume (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>pradsolmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radiation power from SoL (MW)</p></td>
    </tr>
    
    <tr>
        <td>proton_rate_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Proton production rate [particles/m3/sec]</p></td>
    </tr>
    
    <tr>
        <td>psolradmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>SOL radiation power (MW) (<code>stellarator only</code>)</p></td>
    </tr>
    
    <tr>
        <td>pden_plasma_sync_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>synchrotron radiation power per volume (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>i_l_h_threshold</td>
        <td>Input</td>
        <td>integer</td>
        <td>19</td>
        <td><p>switch for L-H mode power threshold scaling to use (see l_h_threshold_powers for list)</p></td>
    </tr>
    
    <tr>
        <td>p_l_h_threshold_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>L-H mode power threshold (MW) (chosen via i_l_h_threshold, and enforced if
 constraint equation 15 is on)</p></td>
    </tr>
    
    <tr>
        <td>l_h_threshold_powers</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>L-H power threshold for various scalings (MW)</p>
<ul>
<li>=1 ITER 1996 scaling: nominal</li>
<li>=2 ITER 1996 scaling: upper bound</li>
<li>=3 ITER 1996 scaling: lower bound</li>
<li>=4 ITER 1997 scaling: excluding elongation</li>
<li>=5 ITER 1997 scaling: including elongation</li>
<li>=6 Martin 2008 scaling: nominal</li>
<li>=7 Martin 2008 scaling: 95% upper bound</li>
<li>=8 Martin 2008 scaling: 95% lower bound</li>
<li>=9 Snipes 2000 scaling: nominal</li>
<li>=10 Snipes 2000 scaling: upper bound</li>
<li>=11 Snipes 2000 scaling: lower bound</li>
<li>=12 Snipes 2000 scaling (closed divertor): nominal</li>
<li>=13 Snipes 2000 scaling (closed divertor): upper bound</li>
<li>=14 Snipes 2000 scaling (closed divertor): lower bound</li>
<li>=15 Hubbard et al. 2012 L-I threshold scaling: nominal</li>
<li>=16 Hubbard et al. 2012 L-I threshold scaling: lower bound</li>
<li>=17 Hubbard et al. 2012 L-I threshold scaling: upper bound</li>
<li>=18 Hubbard et al. 2017 L-I threshold scaling</li>
<li>=19 Martin 2008 aspect ratio corrected scaling: nominal</li>
<li>=20 Martin 2008 aspect ratio corrected scaling: 95% upper bound</li>
<li>=21 Martin 2008 aspect ratio corrected scaling: 95% lower bound</li>
</ul></td>
    </tr>
    
    <tr>
        <td>p_electron_transport_loss_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>electron transport power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pden_electron_transport_loss_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>electron transport power per volume (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>p_ion_transport_loss_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ion transport power (MW)</p></td>
    </tr>
    
    <tr>
        <td>pscalingmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total transport power from scaling law (MW)</p></td>
    </tr>
    
    <tr>
        <td>pden_ion_transport_loss_mw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ion transport power per volume (MW/m3)</p></td>
    </tr>
    
    <tr>
        <td>q0</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>Safety factor on axis</p></td>
    </tr>
    
    <tr>
        <td>q95</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Safety factor at 95% flux surface (iteration variable 18) (unless icurr=2 (ST current scaling),
 in which case q95 = mean edge safety factor qbar)</p></td>
    </tr>
    
    <tr>
        <td>qfuel</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma fuelling rate (nucleus-pairs/s)</p></td>
    </tr>
    
    <tr>
        <td>tauratio</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>tauratio /1.0/ : ratio of He and pellet particle confinement times</p></td>
    </tr>
    
    <tr>
        <td>q95_min</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>lower limit for edge safety factor</p></td>
    </tr>
    
    <tr>
        <td>qstar</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>cylindrical safety factor</p></td>
    </tr>
    
    <tr>
        <td>rad_fraction_sol</td>
        <td>Input</td>
        <td>real</td>
        <td>0.8</td>
        <td><p>SoL radiation fraction</p></td>
    </tr>
    
    <tr>
        <td>rad_fraction_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radiation fraction total = SoL + LCFS radiation / total power deposited in plasma</p></td>
    </tr>
    
    <tr>
        <td>f_nd_alpha_electron</td>
        <td>Input</td>
        <td>real</td>
        <td>0.1</td>
        <td><p>thermal alpha density/electron density (<code>iteration variable 109</code>)</p></td>
    </tr>
    
    <tr>
        <td>f_nd_protium_electrons</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Seeded f_nd_protium_electrons density / electron density.</p></td>
    </tr>
    
    <tr>
        <td>ind_plasma_internal_norm</td>
        <td>Input</td>
        <td>real</td>
        <td>0.9</td>
        <td><p>Plasma normalised internal inductance (calculated from alphaj if <code>iprofile=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>ind_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma inductance (H)</p></td>
    </tr>
    
    <tr>
        <td>rmajor</td>
        <td>Input</td>
        <td>real</td>
        <td>8.14</td>
        <td><p>plasma major radius (m) (<code>iteration variable 3</code>)</p></td>
    </tr>
    
    <tr>
        <td>rminor</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma minor radius (m)</p></td>
    </tr>
    
    <tr>
        <td>f_nd_beam_electron</td>
        <td>Input</td>
        <td>real</td>
        <td>0.005</td>
        <td><p>hot beam density / n_e (<code>iteration variable 7</code>)</p></td>
    </tr>
    
    <tr>
        <td>rncne</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>n_carbon / n_e</p></td>
    </tr>
    
    <tr>
        <td>rndfuel</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>fuel burnup rate (reactions/second)</p></td>
    </tr>
    
    <tr>
        <td>rnfene</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>n_highZ / n_e</p></td>
    </tr>
    
    <tr>
        <td>rnone</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>n_oxygen / n_e</p></td>
    </tr>
    
    <tr>
        <td>f_res_plasma_neo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>neo-classical correction factor to res_plasma</p></td>
    </tr>
    
    <tr>
        <td>res_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma resistance (ohm)</p></td>
    </tr>
    
    <tr>
        <td>t_plasma_res_diffusion</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma current resistive diffusion time (s)</p></td>
    </tr>
    
    <tr>
        <td>a_plasma_surface</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma surface area</p></td>
    </tr>
    
    <tr>
        <td>a_plasma_surface_outboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>outboard plasma surface area</p></td>
    </tr>
    
    <tr>
        <td>i_single_null</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for single null / double null plasma:</p>
<ul>
<li>=0 for double null</li>
<li>=1 for single null (diverted side down)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>f_sync_reflect</td>
        <td>Input</td>
        <td>real</td>
        <td>0.6</td>
        <td><p>synchrotron wall reflectivity factor</p></td>
    </tr>
    
    <tr>
        <td>t_electron_energy_confinement</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>electron energy confinement time (sec)</p></td>
    </tr>
    
    <tr>
        <td>tauee_in</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Input electron energy confinement time (sec) (<code>i_confinement_time=48 only</code>)</p></td>
    </tr>
    
    <tr>
        <td>t_energy_confinement</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>global thermal energy confinement time (sec)</p></td>
    </tr>
    
    <tr>
        <td>t_ion_energy_confinement</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>ion energy confinement time (sec)</p></td>
    </tr>
    
    <tr>
        <td>t_alpha_confinement</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>alpha particle confinement time (sec)</p></td>
    </tr>
    
    <tr>
        <td>f_alpha_energy_confinement</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>alpha particle to energy confinement time ratio</p></td>
    </tr>
    
    <tr>
        <td>te</td>
        <td>Input</td>
        <td>real</td>
        <td>12.9</td>
        <td><p>volume averaged electron temperature (keV) (<code>iteration variable 4</code>)</p></td>
    </tr>
    
    <tr>
        <td>te0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central electron temperature (keV)</p></td>
    </tr>
    
    <tr>
        <td>ten</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>density weighted average electron temperature (keV)</p></td>
    </tr>
    
    <tr>
        <td>ti</td>
        <td>Input</td>
        <td>real</td>
        <td>12.9</td>
        <td><p>volume averaged ion temperature (keV). N.B. calculated from te if <code>tratio &gt; 0.0</code></p></td>
    </tr>
    
    <tr>
        <td>ti0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>central ion temperature (keV)</p></td>
    </tr>
    
    <tr>
        <td>tin</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>density weighted average ion temperature (keV)</p></td>
    </tr>
    
    <tr>
        <td>tratio</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>ion temperature / electron temperature(used to calculate ti if <code>tratio &gt; 0.0</code></p></td>
    </tr>
    
    <tr>
        <td>triang</td>
        <td>Input</td>
        <td>real</td>
        <td>0.36</td>
        <td><p>plasma separatrix triangularity (calculated if <code>i_plasma_geometry = 1, 3-5 or 7</code>)</p></td>
    </tr>
    
    <tr>
        <td>triang95</td>
        <td>Input</td>
        <td>real</td>
        <td>0.24</td>
        <td><p>plasma triangularity at 95% surface (calculated if <code>i_plasma_geometry = 0-2, 6, 8 or 9</code>)</p></td>
    </tr>
    
    <tr>
        <td>vol_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma volume (m3)</p></td>
    </tr>
    
    <tr>
        <td>vs_plasma_burn_required</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>V-s needed during flat-top (heat + burn times) (Wb)</p></td>
    </tr>
    
    <tr>
        <td>v_plasma_loop_burn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Plasma loop voltage during flat-top (V)</p></td>
    </tr>
    
    <tr>
        <td>vshift</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma/device midplane vertical shift - single null</p></td>
    </tr>
    
    <tr>
        <td>vs_plasma_ind_ramp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total plasma inductive flux consumption for plasma current ramp-up (Vs)(Wb)</p></td>
    </tr>
    
    <tr>
        <td>vs_plasma_res_ramp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Plasma resistive flux consumption for plasma current ramp-up (Vs)(Wb)</p></td>
    </tr>
    
    <tr>
        <td>vs_plasma_total_required</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total V-s needed (Wb)</p></td>
    </tr>
    
    <tr>
        <td>wallmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>average neutron wall load (MW/m2)</p></td>
    </tr>
    
    <tr>
        <td>wtgpd</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of fuel used per day (g)</p></td>
    </tr>
    
    <tr>
        <td>a_plasma_poloidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma poloidal cross-sectional area [m^2]</p></td>
    </tr>
    
    <tr>
        <td>zeff</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>plasma effective charge</p></td>
    </tr>
    
    <tr>
        <td>zeffai</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass weighted plasma effective charge</p></td>
    </tr>
    
</table>

## primary_pumping_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>gamma_he</td>
        <td>Input</td>
        <td>real</td>
        <td>1.667</td>
        <td><p>ratio of specific heats for helium (<code>primary_pumping=3</code>)</p></td>
    </tr>
    
    <tr>
        <td>t_in_bb</td>
        <td>Input</td>
        <td>real</td>
        <td>573.13</td>
        <td><p>temperature in FW and blanket coolant at blanket entrance (<code>primary_pumping=3</code>) [K]</p></td>
    </tr>
    
    <tr>
        <td>t_out_bb</td>
        <td>Input</td>
        <td>real</td>
        <td>773.13</td>
        <td><p>temperature in FW and blanket coolant at blanket exit (<code>primary_pumping=3</code>) [K]</p></td>
    </tr>
    
    <tr>
        <td>p_he</td>
        <td>Input</td>
        <td>real</td>
        <td>8000000.0</td>
        <td><p>pressure in FW and blanket coolant at pump exit (<code>primary_pumping=3</code>) [Pa]</p></td>
    </tr>
    
    <tr>
        <td>dp_he</td>
        <td>Input</td>
        <td>real</td>
        <td>550000.0</td>
        <td><p>pressure drop in FW and blanket coolant including heat exchanger and pipes (<code>primary_pumping=3</code>) [Pa]</p></td>
    </tr>
    
    <tr>
        <td>dp_fw_blkt</td>
        <td>Input</td>
        <td>real</td>
        <td>150000.0</td>
        <td><p>pressure drop in FW and blanket coolant including heat exchanger and pipes (<code>primary_pumping=3</code>) [Pa]</p></td>
    </tr>
    
    <tr>
        <td>dp_fw</td>
        <td>Input</td>
        <td>real</td>
        <td>150000.0</td>
        <td><p>pressure drop in FW coolant including heat exchanger and pipes (<code>primary_pumping=3</code>) [Pa]</p></td>
    </tr>
    
    <tr>
        <td>dp_blkt</td>
        <td>Input</td>
        <td>real</td>
        <td>3500.0</td>
        <td><p>pressure drop in blanket coolant including heat exchanger and pipes (<code>primary_pumping=3</code>) [Pa]</p></td>
    </tr>
    
    <tr>
        <td>dp_liq</td>
        <td>Input</td>
        <td>real</td>
        <td>10000000.0</td>
        <td><p>pressure drop in liquid metal blanket coolant including heat exchanger and pipes (<code>primary_pumping=3</code>) [Pa]</p></td>
    </tr>
    
    <tr>
        <td>htpmw_fw_blkt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mechanical pumping power for FW and blanket including heat exchanger and
 pipes (<code>primary_pumping=3</code>) [MW]</p></td>
    </tr>
    
</table>

## pulse_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>bctmp</td>
        <td>Input</td>
        <td>real</td>
        <td>320.0</td>
        <td><p>first wall bulk coolant temperature (C)</p></td>
    </tr>
    
    <tr>
        <td>dtstor</td>
        <td>Input</td>
        <td>real</td>
        <td>300.0</td>
        <td><p>maximum allowable temperature change in stainless steel thermal storage block (K) (<code>istore=3</code>)</p></td>
    </tr>
    
    <tr>
        <td>istore</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for thermal storage method:</p>
<ul>
<li>=1 option 1 of Electrowatt report, AEA FUS 205</li>
<li>=2 option 2 of Electrowatt report, AEA FUS 205</li>
<li>=3 stainless steel block</li>
</ul></td>
    </tr>
    
    <tr>
        <td>itcycl</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for first wall axial stress model:</p>
<ul>
<li>=1 total axial constraint, no bending</li>
<li>=2 no axial constraint, no bending</li>
<li>=3 no axial constraint, bending</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_pulsed_plant</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for reactor model:</p>
<ul>
<li>=0 continuous operation</li>
<li>=1 pulsed operation</li>
</ul></td>
    </tr>
    
</table>

## rebco_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>rebco_thickness</td>
        <td>Input</td>
        <td>real</td>
        <td>1e-06</td>
        <td><p>thickness of REBCO layer in tape (m) (<code>iteration variable 138</code>)</p></td>
    </tr>
    
    <tr>
        <td>copper_thick</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0001</td>
        <td><p>thickness of copper layer in tape (m) (<code>iteration variable 139</code>)</p></td>
    </tr>
    
    <tr>
        <td>hastelloy_thickness</td>
        <td>Input</td>
        <td>real</td>
        <td>5e-05</td>
        <td><p>thickness of Hastelloy layer in tape (m)</p></td>
    </tr>
    
    <tr>
        <td>tape_width</td>
        <td>Input</td>
        <td>real</td>
        <td>0.004</td>
        <td><p>Mean width of tape (m)</p></td>
    </tr>
    
    <tr>
        <td>tape_thickness</td>
        <td>Input</td>
        <td>real</td>
        <td>6.5e-05</td>
        <td><p>thickness of tape, inc. all layers (hts, copper, substrate, etc.) (m)</p></td>
    </tr>
    
    <tr>
        <td>croco_od</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Outer diameter of CroCo strand (m)</p></td>
    </tr>
    
    <tr>
        <td>croco_id</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inner diameter of CroCo copper tube (m)</p></td>
    </tr>
    
    <tr>
        <td>croco_thick</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0025</td>
        <td><p>Thickness of CroCo copper tube (m) (<code>iteration variable 158</code>)</p></td>
    </tr>
    
    <tr>
        <td>copper_rrr</td>
        <td>Input</td>
        <td>real</td>
        <td>100.0</td>
        <td><p>residual resistivity ratio copper in TF superconducting cable</p></td>
    </tr>
    
    <tr>
        <td>coppera_m2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil current / copper area (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>coppera_m2_max</td>
        <td>Input</td>
        <td>real</td>
        <td>100000000.0</td>
        <td><p>Maximum TF coil current / copper area (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>f_coppera_m2</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for constraint 75: TF coil current / copper area &lt; copperA_m2_max</p></td>
    </tr>
    
    <tr>
        <td>copperaoh_m2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>CS coil current / copper area (A/m2) (<code>sweep variable 61</code>)</p></td>
    </tr>
    
    <tr>
        <td>copperaoh_m2_max</td>
        <td>Input</td>
        <td>real</td>
        <td>100000000.0</td>
        <td><p>Maximum CS coil current / copper area (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>f_copperaoh_m2</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for constraint 88: CS coil current / copper area &lt; copperA_m2_max</p></td>
    </tr>
    
    <tr>
        <td>stack_thickness</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>tapes</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>rebco_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>copper_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>hastelloy_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>solder_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>croco_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## reinke_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>vcritx</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## reinke_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>impvardiv</td>
        <td>Input</td>
        <td>integer</td>
        <td>9</td>
        <td><p>Index of impurity to be iterated for Reinke divertor detachment criterion</p></td>
    </tr>
    
    <tr>
        <td>lhat</td>
        <td>Input</td>
        <td>real</td>
        <td>4.33</td>
        <td><p>Connection length factor L|| = lhat qstar R for Reinke criterion, default value from
 Post et al. 1995 J. Nucl. Mat.  220-2 1014</p></td>
    </tr>
    
    <tr>
        <td>fzmin</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Minimum impurity fraction necessary for detachment. This is the impurity at the SOL/Div.</p></td>
    </tr>
    
    <tr>
        <td>fzactual</td>
        <td>Input</td>
        <td>real</td>
        <td>0.001</td>
        <td><p>Actual impurity fraction of divertor impurity (impvardiv) in the SoL (taking
 impurity_enrichment into account) (<code>iteration variable 148</code>)</p></td>
    </tr>
    
    <tr>
        <td>reinke_mode</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for Reinke criterion H/I mode:</p>
<ul>
<li>=0 H-mode</li>
<li>=1 I-mode</li>
</ul></td>
    </tr>
    
</table>

## scan_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>ipnscns</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>1000</td>
        <td><p>Maximum number of scan points</p></td>
    </tr>
    
    <tr>
        <td>ipnscnv</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>81</td>
        <td><p>Number of available scan variables</p></td>
    </tr>
    
    <tr>
        <td>noutvars</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>84</td>
        <td></td>
    </tr>
    
    <tr>
        <td>width</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>110</td>
        <td></td>
    </tr>
    
    <tr>
        <td>scan_dim</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>1-D or 2-D scan switch (1=1D, 2=2D)</p></td>
    </tr>
    
    <tr>
        <td>isweep</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Number of scan points to calculate</p></td>
    </tr>
    
    <tr>
        <td>isweep_2</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Number of 2D scan points to calculate</p></td>
    </tr>
    
    <tr>
        <td>nsweep</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch denoting quantity to scan:<UL>
         <LI> 1  aspect
         <LI> 2  hldivlim
         <LI> 3  pnetelin
         <LI> 4  hfact
         <LI> 5  oacdcp
         <LI> 6  walalw
         <LI> 7  beamfus0
         <LI> 8  fqval
         <LI> 9  te
         <LI> 10 boundu(15: fvs)
         <LI> 11 beta_norm_max
         <LI> 12 bootstrap_current_fraction_max
         <LI> 13 boundu(10: hfact)
         <LI> 14 fiooic
         <LI> 15 fjprot
         <LI> 16 rmajor
         <LI> 17 bmxlim
         <LI> 18 gammax
         <LI> 19 boundl(16: dr_cs)
         <LI> 20 t_burn_min
         <LI> 21 not used
         <LI> 22 cfactr (N.B. requires iavail=0)
         <LI> 23 boundu(72: fipir)
         <LI> 24 powfmax
         <LI> 25 kappa
         <LI> 26 triang
         <LI> 27 tbrmin (for blktmodel &gt; 0 only)
         <LI> 28 bt
         <LI> 29 coreradius
         <LI> 30 fimpvar # OBSOLETE
         <LI> 31 f_alpha_energy_confinement_min
         <LI> 32 epsvmc
         <LI> 33 ttarget
         <LI> 34 qtargettotal
         <LI> 35 lambda_q_omp
         <LI> 36 lambda_target
         <LI> 37 lcon_factor
         <LI> 38 Neon upper limit
         <LI> 39 Argon upper limit
         <LI> 40 Xenon upper limit
         <LI> 41 dr_blkt_outboard
         <LI> 42 Argon fraction fimp(9)
         <LI> 43 normalised minor radius at which electron cyclotron current drive is maximum
         <LI> 44 Allowable maximum shear stress (Tresca) in tf coil structural material
         <LI> 45 Minimum allowable temperature margin ; tf coils
         <LI> 46 boundu(150) fgwsep
         <LI> 47 impurity_enrichment(9) Argon impurity enrichment
         <LI> 48 TF coil - n_pancake (integer turn winding pack)
         <LI> 49 TF coil - n_layer (integer turn winding pack)
         <LI> 50 Xenon fraction fimp(13)
         <LI> 51 Power fraction to lower DN Divertor ftar
         <LI> 52 SoL radiation fraction
         <LI> 54 GL_nbti upper critical field at 0 Kelvin
         <LI> 55 <code>dr_shld_inboard</code> : Inboard neutron shield thickness
         <LI> 56 crypmw_max: Maximum cryogenic power (ixx=164, ixc=87)
         <LI> 57 <code>bt</code> lower boundary
         <LI> 58 <code>dr_fw_plasma_gap_inboard</code> : Inboard plasma-first wall gap
         <LI> 59 <code>dr_fw_plasma_gap_outboard</code> : Outboard plasma-first wall gap
         <LI> 60 sig_tf_wp_max: Allowable stress in TF Coil conduit (Tresca)
         <LI> 61 copperaoh_m2_max : CS coil current / copper area
         <LI> 62 coheof : CS coil current density at EOF
         <LI> 63 dr_cs : CS thickness (m)
         <LI> 64 ohhghf : CS height (m)
         <LI> 65 n_cycle_min : Minimum cycles for CS stress model constraint 90
         <LI> 66 oh_steel_frac: Steel fraction in CS coil
         <LI> 67 t_crack_vertical: Initial crack vertical dimension (m) </UL>
         <LI> 68 <code>inlet_temp_liq' : Inlet temperature of blanket liquid metal coolant/breeder (K)
         &lt;LI&gt; 69</code>outlet_temp_liq' : Outlet temperature of blanket liquid metal coolant/breeder (K)
         <LI> 70 <code>blpressure_liq' : Blanket liquid metal breeder/coolant pressure (Pa)
         &lt;LI&gt; 71</code>n_liq_recirc' : Selected number of liquid metal breeder recirculations per day
         <LI> 72 <code>bz_channel_conduct_liq' : Conductance of liquid metal breeder duct walls (A V-1 m-1)
         &lt;LI&gt; 73</code>pnuc_fw_ratio_dcll' : Ratio of FW nuclear power as fraction of total (FW+BB)
         <LI> 74 `f_nuc_pow_bz_struct' : Fraction of BZ power cooled by primary coolant for dual-coolant balnket
         <LI> 75 dx_fw_module : pitch of first wall cooling channels (m)
         <LI> 76 etath : Thermal conversion eff.
         <LI> 77 startupratio : Gyrotron redundancy
         <LI> 78 fkind : Multiplier for Nth of a kind costs
         <LI> 79 etaech : ECH wall plug to injector efficiency</p></td>
    </tr>
    
    <tr>
        <td>nsweep_2</td>
        <td>Input</td>
        <td>integer</td>
        <td>3</td>
        <td><p>nsweep_2 /3/ : switch denoting quantity to scan for 2D scan:</p></td>
    </tr>
    
    <tr>
        <td>sweep</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>sweep(ipnscns) /../: actual values to use in scan</p></td>
    </tr>
    
    <tr>
        <td>sweep_2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>sweep_2(ipnscns) /../: actual values to use in 2D scan</p></td>
    </tr>
    
    <tr>
        <td>first_call_1d</td>
        <td>Input</td>
        <td>logical</td>
        <td>1</td>
        <td></td>
    </tr>
    
    <tr>
        <td>first_call_2d</td>
        <td>Input</td>
        <td>logical</td>
        <td>1</td>
        <td></td>
    </tr>
    
</table>

## sctfcoil_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>tf_fit_t</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Dimensionless winding pack width</p></td>
    </tr>
    
    <tr>
        <td>tf_fit_z</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Dimensionless winding pack radial thickness</p></td>
    </tr>
    
    <tr>
        <td>tf_fit_y</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Ratio of peak field with ripple to nominal axisymmetric peak field</p></td>
    </tr>
    
    <tr>
        <td>tfc_current</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Current in each TF coil</p></td>
    </tr>
    
    <tr>
        <td>awpc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total cross-sectional area of winding pack including
 GW insulation and insertion gap [m2]</p></td>
    </tr>
    
    <tr>
        <td>awptf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total cross-sectional area of winding pack without
 ground insulation and insertion gap [m2]</p></td>
    </tr>
    
    <tr>
        <td>a_tf_steel</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard coil steel coil cross-sectional area [m2]</p></td>
    </tr>
    
    <tr>
        <td>a_tf_ins</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard coil insulation cross-section per coil [m2]</p></td>
    </tr>
    
    <tr>
        <td>f_tf_steel</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard coil steel fraction [-]</p></td>
    </tr>
    
    <tr>
        <td>f_tf_ins</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard coil insulation fraction [-]</p></td>
    </tr>
    
    <tr>
        <td>h_cp_top</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Vertical distance from the midplane to the top of the tapered section [m]</p></td>
    </tr>
    
    <tr>
        <td>r_tf_outboard_in</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial position of plasma-facing edge of TF coil outboard leg [m]</p></td>
    </tr>
    
    <tr>
        <td>r_tf_outboard_out</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial position of outer edge of TF coil inboard leg [m]</p></td>
    </tr>
    
    <tr>
        <td>r_wp_inner</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial position of inner edge and centre of winding pack [m]</p></td>
    </tr>
    
    <tr>
        <td>r_wp_outer</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial position of outer edge and centre of winding pack [m]</p></td>
    </tr>
    
    <tr>
        <td>r_wp_centre</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial position of centre and centre of winding pack [m]</p></td>
    </tr>
    
    <tr>
        <td>dr_tf_wp_top</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Conductor layer radial thickness at centercollumn top [m]
 Ground insulation layer included, only defined for itart = 1</p></td>
    </tr>
    
    <tr>
        <td>vol_ins_cp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>CP turn insulation volume [m3]</p></td>
    </tr>
    
    <tr>
        <td>vol_gr_ins_cp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>CP ground insulation volume [m3]</p></td>
    </tr>
    
    <tr>
        <td>vol_case_cp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Volume of the CP outer casing cylinder</p></td>
    </tr>
    
    <tr>
        <td>t_wp_toroidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Minimal toroidal thickness of of winding pack [m]</p></td>
    </tr>
    
    <tr>
        <td>t_wp_toroidal_av</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Averaged toroidal thickness of of winding pack [m]</p></td>
    </tr>
    
    <tr>
        <td>t_lat_case_av</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Average lateral casing thickness [m]</p></td>
    </tr>
    
    <tr>
        <td>a_case_front</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Front casing area [m2]</p></td>
    </tr>
    
    <tr>
        <td>a_case_nose</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Nose casing area [m2]</p></td>
    </tr>
    
    <tr>
        <td>a_ground_ins</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Inboard mid-plane cross-section area of the WP ground insulation [m2]</p></td>
    </tr>
    
    <tr>
        <td>a_leg_ins</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF ouboard leg turn insulation area per coil [m2]</p></td>
    </tr>
    
    <tr>
        <td>a_leg_gr_ins</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF outboard leg ground insulation area per coil [m2]</p></td>
    </tr>
    
    <tr>
        <td>a_leg_cond</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Exact TF ouboard leg conductor area [m2]</p></td>
    </tr>
    
    <tr>
        <td>theta_coil</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Half toroidal angular extent of a single TF coil inboard leg</p></td>
    </tr>
    
    <tr>
        <td>tan_theta_coil</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Tan half toroidal angular extent of a single TF coil inboard leg</p></td>
    </tr>
    
    <tr>
        <td>t_conductor_radial</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Conductor area radial and toroidal dimension (integer turn only) [m]</p></td>
    </tr>
    
    <tr>
        <td>t_conductor_toroidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Conductor area radial and toroidal dimension (integer turn only) [m]</p></td>
    </tr>
    
    <tr>
        <td>t_cable_radial</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Cable area radial and toroidal dimension (integer turn only) [m]</p></td>
    </tr>
    
    <tr>
        <td>t_cable_toroidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Cable area radial and toroidal dimension (integer turn only) [m]</p></td>
    </tr>
    
    <tr>
        <td>t_turn_radial</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Turn radial and toroidal dimension (integer turn only) [m]</p></td>
    </tr>
    
    <tr>
        <td>t_turn_toroidal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Turn radial and toroidal dimension (integer turn only) [m]</p></td>
    </tr>
    
    <tr>
        <td>t_cable</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Cable area averaged dimension (square shape) [m]</p></td>
    </tr>
    
    <tr>
        <td>vforce_inboard_tot</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Total inboard vertical tension (all coils) [N]</p></td>
    </tr>
    
    <tr>
        <td>vv_stress_quench</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>The Tresca stress experienced by the Vacuum Vessel when the SCTF coil quenches [Pa]</p></td>
    </tr>
    
    <tr>
        <td>copper*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>hastelloy*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>solder*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>jacket*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>helium*</td>
        <td>Variable</td>
        <td>type</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>croco_strand_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>croco_strand_critical_current</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_copper_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_copper_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_copper_bar_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_hastelloy_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_hastelloy_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_helium_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_helium_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_solder_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_solder_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_jacket_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_jacket_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_rebco_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_rebco_fraction</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_critical_current</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_acs</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>conductor_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Area of cable space inside jacket</p></td>
    </tr>
    
    <tr>
        <td>t1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>time2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>tau2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>estotft</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>is_leg_cp_temp_same</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## stellarator_module
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>f_n</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>f_r</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>f_aspect</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>f_b</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>f_i</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>f_a</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>first_call</td>
        <td>Variable</td>
        <td>logical</td>
        <td>.true.</td>
        <td></td>
    </tr>
    
    <tr>
        <td>first_call_stfwbs</td>
        <td>Variable</td>
        <td>logical</td>
        <td>.true.</td>
        <td></td>
    </tr>
    
</table>

## stellarator_configuration
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>stella_config_name</td>
        <td>Output</td>
        <td>character</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_symmetry</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_coilspermodule</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_rmajor_ref</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_rminor_ref</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_coil_rmajor</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_coil_rminor</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_aspect_ref</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_bt_ref</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_wp_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_wp_bmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_i0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_a1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_a2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_dmin</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_inductance</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_coilsurface</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_coillength</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_max_portsize_width</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_maximal_coil_height</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_min_plasma_coil_distance</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_derivative_min_lcfs_coils_dist</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_vol_plasma</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_plasma_surface</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_wp_ratio</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_max_force_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_max_force_density_mnm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_min_bend_radius</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_epseff</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_max_lateral_force_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_max_radial_force_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_centering_force_max_mn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_centering_force_min_mn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_centering_force_avg_mn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>stella_config_neutron_peakfactor</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## stellarator_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>istell</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for stellarator option (set via <code>device.dat</code>):</p>
<ul>
<li>=0 use tokamak model</li>
<li>=1 use stellarator model: Helias5</li>
<li>=2 use stellarator model: Helias4</li>
<li>=3 use stellarator model: Helias3</li>
<li>=4 use stellarator model: Wendelstein 7-X with 50 Coils</li>
<li>=5 use stellarator model: Wendelstein 7-X with 30 Coils</li>
<li>=6 use stellarator model: Use stella_conf.json file (any modulear stellarator, see documentation)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>bmn</td>
        <td>Input</td>
        <td>real</td>
        <td>0.001</td>
        <td><p>relative radial field perturbation</p></td>
    </tr>
    
    <tr>
        <td>f_asym</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>divertor heat load peaking factor</p></td>
    </tr>
    
    <tr>
        <td>f_rad</td>
        <td>Input</td>
        <td>real</td>
        <td>0.85</td>
        <td><p>radiated power fraction in SOL</p></td>
    </tr>
    
    <tr>
        <td>f_w</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>island size fraction factor</p></td>
    </tr>
    
    <tr>
        <td>fdivwet</td>
        <td>Input</td>
        <td>real</td>
        <td>0.333333333333333</td>
        <td><p>wetted fraction of the divertor area</p></td>
    </tr>
    
    <tr>
        <td>flpitch</td>
        <td>Input</td>
        <td>real</td>
        <td>0.001</td>
        <td><p>field line pitch (rad)</p></td>
    </tr>
    
    <tr>
        <td>hportamax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum available area for horizontal ports (m2)</p></td>
    </tr>
    
    <tr>
        <td>hportpmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum available poloidal extent for horizontal ports (m)</p></td>
    </tr>
    
    <tr>
        <td>hporttmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum available toroidal extent for horizontal ports (m)</p></td>
    </tr>
    
    <tr>
        <td>iotabar</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>rotational transform (reciprocal of tokamak q) for stellarator confinement time scaling laws</p></td>
    </tr>
    
    <tr>
        <td>isthtr</td>
        <td>Input</td>
        <td>integer</td>
        <td>3</td>
        <td><p>Switch for stellarator auxiliary heating method:</p>
<ul>
<li>= 1electron cyclotron resonance heating</li>
<li>= 2lower hybrid heating</li>
<li>= 3neutral beam injection</li>
</ul></td>
    </tr>
    
    <tr>
        <td>m_res</td>
        <td>Input</td>
        <td>integer</td>
        <td>5</td>
        <td><p>poloidal resonance number (1)</p></td>
    </tr>
    
    <tr>
        <td>max_gyrotron_frequency</td>
        <td>Input</td>
        <td>real</td>
        <td>1000000000.0</td>
        <td><p>Maximal available gyrotron frequency (input parameter) (Hz)</p></td>
    </tr>
    
    <tr>
        <td>n_res</td>
        <td>Input</td>
        <td>integer</td>
        <td>5</td>
        <td><p>toroidal resonance number (1)</p></td>
    </tr>
    
    <tr>
        <td>shear</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>magnetic shear, derivative of iotabar (1)</p></td>
    </tr>
    
    <tr>
        <td>te0_ecrh_achievable</td>
        <td>Input</td>
        <td>real</td>
        <td>100.0</td>
        <td><p>maximal central electron temperature as achievable by the ECRH, input. (keV)</p></td>
    </tr>
    
    <tr>
        <td>vportamax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum available area for vertical ports (m2)</p></td>
    </tr>
    
    <tr>
        <td>vportpmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum available poloidal extent for vertical ports (m)</p></td>
    </tr>
    
    <tr>
        <td>vporttmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>maximum available toroidal extent for vertical ports (m)</p></td>
    </tr>
    
    <tr>
        <td>powerht_constraint</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>powerscaling_constraint</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
</table>

## structure_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>aintmass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>intercoil structure mass (kg)</p></td>
    </tr>
    
    <tr>
        <td>clgsmass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>gravity support structure for TF coil, PF coil and intercoil support systems (kg)</p></td>
    </tr>
    
    <tr>
        <td>coldmass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total mass of components at cryogenic temperatures (kg)</p></td>
    </tr>
    
    <tr>
        <td>fncmass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>PF coil outer support fence mass (kg)</p></td>
    </tr>
    
    <tr>
        <td>gsmass</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>reactor core gravity support mass (kg)</p></td>
    </tr>
    
</table>

## tfcoil_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>acasetf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>external case area per coil (inboard leg) (m2)</p></td>
    </tr>
    
    <tr>
        <td>acasetfo</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>external case area per coil (outboard leg) (m2)</p></td>
    </tr>
    
    <tr>
        <td>acndttf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>area of the cable conduit (m2)</p></td>
    </tr>
    
    <tr>
        <td>acond</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Winding pack conductor area [m2]
 Does not include the area of voids and central helium channel</p></td>
    </tr>
    
    <tr>
        <td>acstf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Cable space area (per turn)  [m2]
 Includes the area of voids and central helium channel</p></td>
    </tr>
    
    <tr>
        <td>insulation_area</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>single turn insulation area (m2)</p></td>
    </tr>
    
    <tr>
        <td>aiwp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>winding pack turn insulation area per coil (m2)</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_case_max</td>
        <td>Input</td>
        <td>real</td>
        <td>600000000.0</td>
        <td><p>Allowable maximum shear stress (Tresca criterion) in TF coil case (Pa)</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_wp_max</td>
        <td>Input</td>
        <td>real</td>
        <td>600000000.0</td>
        <td><p>Allowable maximum shear stress (Tresca criterion) in TF coil conduit (Pa)</p>
<p>Allowable Tresca stress in TF coil structural material (Pa)</p></td>
    </tr>
    
    <tr>
        <td>a_tf_leg_outboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>outboard TF leg area (m2)</p></td>
    </tr>
    
    <tr>
        <td>aswp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>winding pack structure area (m2)</p></td>
    </tr>
    
    <tr>
        <td>avwp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>winding pack void (He coolant) area (m2)</p></td>
    </tr>
    
    <tr>
        <td>awphec</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>winding pack He coil area (m2)</p></td>
    </tr>
    
    <tr>
        <td>bcritsc</td>
        <td>Input</td>
        <td>real</td>
        <td>24.0</td>
        <td><p>upper critical field (T) for Nb3Sn superconductor at zero temperature and
 strain (<code>i_tf_sc_mat=4, =bc20m</code>)</p></td>
    </tr>
    
    <tr>
        <td>bmaxtf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mean peak field at TF coil (T)</p></td>
    </tr>
    
    <tr>
        <td>bmaxtfrp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>peak field at TF conductor with ripple (T)</p></td>
    </tr>
    
    <tr>
        <td>casestr</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>case strain</p></td>
    </tr>
    
    <tr>
        <td>casthi</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>inboard TF coil case plasma side thickness (m) (calculated for stellarators)</p></td>
    </tr>
    
    <tr>
        <td>casthi_fraction</td>
        <td>Input</td>
        <td>real</td>
        <td>0.05</td>
        <td><p>inboard TF coil case plasma side thickness as a fraction of dr_tf_inboard</p></td>
    </tr>
    
    <tr>
        <td>casthi_is_fraction</td>
        <td>Output</td>
        <td>logical</td>
        <td>-</td>
        <td><p>logical switch to make casthi a fraction of TF coil thickness (<code>casthi_fraction</code>)</p></td>
    </tr>
    
    <tr>
        <td>casths</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>inboard TF coil sidewall case thickness (m) (calculated for stellarators)</p></td>
    </tr>
    
    <tr>
        <td>casths_fraction</td>
        <td>Input</td>
        <td>real</td>
        <td>0.06</td>
        <td><p>inboard TF coil sidewall case thickness as a fraction of tftort</p></td>
    </tr>
    
    <tr>
        <td>tfc_sidewall_is_fraction</td>
        <td>Output</td>
        <td>logical</td>
        <td>-</td>
        <td><p>logical switch to make casths a fraction of TF coil thickness (<code>casths_fraction</code>)</p></td>
    </tr>
    
    <tr>
        <td>t_conductor</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Conductor (cable + steel conduit) area averaged dimension [m]</p></td>
    </tr>
    
    <tr>
        <td>t_turn_tf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil turn edge length including turn insulation [m]
   If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
   equivelent size is use to calculated this quantity
   If the t_turn_tf is non zero, cpttf is calculated</p></td>
    </tr>
    
    <tr>
        <td>t_turn_tf_is_input</td>
        <td>Output</td>
        <td>logical</td>
        <td>-</td>
        <td><p>Boolean switch to activated when the user set the TF coil turn dimensions
 Not an input</p></td>
    </tr>
    
    <tr>
        <td>f_t_turn_tf</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>f-value for TF turn edge length constraint
  If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
  equivelent size is use for this constraint
  iteration variable ixc = 175
  constraint equation icc = 86</p></td>
    </tr>
    
    <tr>
        <td>t_turn_tf_max</td>
        <td>Input</td>
        <td>real</td>
        <td>0.05000000074505806</td>
        <td><p>TF turn edge length including turn insulation upper limit [m]
 If the turn is not a square (i_tf_turns_integer = 1) a squared turn of
 equivelent size is use for this constraint
 constraint equation icc = 86</p></td>
    </tr>
    
    <tr>
        <td>t_cable_tf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil superconducting cable squared/rounded dimensions [m]
   If the turn is not a square (i_tf_turns_integer = 1) a squared cable of
   equivelent size is use to calculated this quantity
   If the t_cable_tf is non zero, cpttf is calculated</p></td>
    </tr>
    
    <tr>
        <td>t_cable_tf_is_input</td>
        <td>Output</td>
        <td>logical</td>
        <td>-</td>
        <td><p>Boolean switch to activated when the user set the TF coil cable dimensions
 Not an input</p></td>
    </tr>
    
    <tr>
        <td>acs</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Area of space inside conductor (m2)</p></td>
    </tr>
    
    <tr>
        <td>cdtfleg</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF outboard leg current density (A/m2) (resistive coils only)</p></td>
    </tr>
    
    <tr>
        <td>cforce</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>centering force on inboard leg (per coil) (N/m)</p></td>
    </tr>
    
    <tr>
        <td>cplen</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>length of TF coil inboard leg ('centrepost') (<code>i_tf_sup = 1</code>)</p></td>
    </tr>
    
    <tr>
        <td>cpttf</td>
        <td>Input</td>
        <td>real</td>
        <td>70000.0</td>
        <td><p>TF coil current per turn (A). (calculated for stellarators) (calculated for
 integer-turn TF coils <code>i_tf_turns_integer=1</code>) (<code>iteration variable 60</code>)</p></td>
    </tr>
    
    <tr>
        <td>cpttf_max</td>
        <td>Input</td>
        <td>real</td>
        <td>90000.0</td>
        <td><p>Max TF coil current per turn [A]. (for stellarators and <code>i_tf_turns_integer=1</code>)
 (<code>constraint equation 77</code>)</p></td>
    </tr>
    
    <tr>
        <td>dcase</td>
        <td>Input</td>
        <td>real</td>
        <td>8000.0</td>
        <td><p>density of coil case (kg/m3)</p></td>
    </tr>
    
    <tr>
        <td>dcond</td>
        <td>Input</td>
        <td>real</td>
        <td>[6080. 6080. 6070. 6080. 6080. 8500. 6070. 8500. 8500.]</td>
        <td><p>density of superconductor type given by i_tf_sc_mat/isumatoh/isumatpf (kg/m3)</p></td>
    </tr>
    
    <tr>
        <td>dcondins</td>
        <td>Input</td>
        <td>real</td>
        <td>1800.0</td>
        <td><p>density of conduit + ground-wall insulation (kg/m3)</p></td>
    </tr>
    
    <tr>
        <td>dhecoil</td>
        <td>Input</td>
        <td>real</td>
        <td>0.005</td>
        <td><p>diameter of central helium channel in TF winding (m)</p></td>
    </tr>
    
    <tr>
        <td>estotftgj</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total stored energy in the toroidal field (GJ)</p></td>
    </tr>
    
    <tr>
        <td>b_crit_upper_nbti</td>
        <td>Input</td>
        <td>real</td>
        <td>14.86</td>
        <td><p>upper critical field of GL_nbti</p></td>
    </tr>
    
    <tr>
        <td>t_crit_nbti</td>
        <td>Input</td>
        <td>real</td>
        <td>9.04</td>
        <td><p>critical temperature of GL_nbti</p></td>
    </tr>
    
    <tr>
        <td>max_force_density</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Maximal (WP averaged) force density in TF coils at 1 point. (MN/m3)</p></td>
    </tr>
    
    <tr>
        <td>fcutfsu</td>
        <td>Input</td>
        <td>real</td>
        <td>0.69</td>
        <td><p>copper fraction of cable conductor (TF coils)
 (iteration variable 59)</p></td>
    </tr>
    
    <tr>
        <td>fhts</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>technology adjustment factor for critical current density fit for isumat..=2
 Bi-2212 superconductor, to describe the level of technology assumed (i.e. to
 account for stress, fatigue, radiation, AC losses, joints or manufacturing
 variations; 1.0 would be very optimistic)</p></td>
    </tr>
    
    <tr>
        <td>insstrain</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radial strain in insulator</p></td>
    </tr>
    
    <tr>
        <td>i_tf_stress_model</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for the TF coil stress model
   0 : Generalized plane strain formulation, Issues #977 and #991, O(n^3)
   1 : Old plane stress model (only for SC)
   2 : Axisymmetric extended plane strain, Issues #1414 and #998, O(n)</p></td>
    </tr>
    
    <tr>
        <td>i_tf_tresca</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for TF coil conduit Tresca stress criterion:
   0 : Tresca (no adjustment);
   1 : Tresca with CEA adjustment factors (radial+2%, vertical+60%) </UL></p></td>
    </tr>
    
    <tr>
        <td>i_tf_wp_geom</td>
        <td>Input</td>
        <td>integer</td>
        <td>-1</td>
        <td><p>Switch for TF WP geometry selection
   0 : Rectangular geometry
   1 : Double rectangular geometry
   2 : Trapezoidal geometry (constant lateral casing thickness)
 Default setting for backward compatibility
   if i_tf_turns_integer = 0 : Double rectangular
   if i_tf_turns_integer = 1 : Rectangular</p></td>
    </tr>
    
    <tr>
        <td>i_tf_case_geom</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for TF case geometry selection
   0 : Circular front case (ITER design)
   1 : Straight front case</p></td>
    </tr>
    
    <tr>
        <td>i_tf_turns_integer</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for TF coil integer/non-integer turns:
   0 : non-integer turns
   1 : integer turns</p></td>
    </tr>
    
    <tr>
        <td>i_tf_sc_mat</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for superconductor material in TF coils:</p>
<ul>
<li>=1 ITER Nb3Sn critical surface model with standard
   ITER parameters</li>
<li>=2 Bi-2212 high temperature superconductor (range of
   validity T &lt; 20K, adjusted field b &lt; 104 T, B &gt; 6 T)</li>
<li>=3 NbTi</li>
<li>=4 ITER Nb3Sn model with user-specified parameters</li>
<li>=5 WST Nb3Sn parameterisation</li>
<li>=6 REBCO HTS tape in CroCo strand</li>
<li>=7 Durham Ginzburg-Landau critical surface model for Nb-Ti</li>
<li>=8 Durham Ginzburg-Landau critical surface model for REBCO</li>
<li>=9 Hazelton experimental data + Zhai conceptual model for REBCO</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_tf_sup</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for TF coil conductor model:</p>
<ul>
<li>=0 copper</li>
<li>=1 superconductor</li>
<li>=2 Cryogenic aluminium</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_tf_shape</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for TF coil toroidal shape:</p>
<ul>
<li>=0  Default value : Picture frame coil for TART / PROCESS D-shape for non itart</li>
<li>=1  PROCESS D-shape : parametrise with 2 arcs</li>
<li>=2  Picture frame coils</li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_tf_cond_eyoung_axial</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for the behavior of the TF coil conductor elastic axial properties</p>
<ul>
<li>=0  Young's modulus is set to zero, and the conductor is not considered
       in the stress calculation. This corresponds to the case that the
       conductor is much less stiff than the conduit, or the case that the
       conductor is prevented (isolated) from taking axial loads.</li>
<li>=1  Elastic properties are set by user input, using the variable
       <code>eyoung_cond_axial</code></li>
<li>=2  Elastic properties are set to reasonable defaults taking into
       account the superconducting material <code>i_tf_sc_mat</code></li>
</ul></td>
    </tr>
    
    <tr>
        <td>i_tf_cond_eyoung_trans</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for the behavior of the elastic properties of the TF coil
 conductorin the transverse direction. Only active if
 <code>i_tf_cond_eyoung_axial == 2</code></p>
<ul>
<li>=0  Cable not potted in solder. Transverse Young's modulus set to zero.</li>
<li>=1  Cable potted in solder. If <code>i_tf_cond_eyoung_axial == 2</code>, the
       transverse Young's modulus of the conductor is equal to the axial,
       which is set to a sensible material-dependent default.</li>
</ul></td>
    </tr>
    
    <tr>
        <td>n_pancake</td>
        <td>Input</td>
        <td>integer</td>
        <td>10</td>
        <td><p>Number of pancakes in TF coil. Only used if <code>i_tf_turns_integer=1</code></p></td>
    </tr>
    
    <tr>
        <td>n_layer</td>
        <td>Input</td>
        <td>integer</td>
        <td>20</td>
        <td><p>Number of layers in TF coil. Only used if <code>i_tf_turns_integer=1</code></p></td>
    </tr>
    
    <tr>
        <td>n_rad_per_layer</td>
        <td>Input</td>
        <td>integer</td>
        <td>100</td>
        <td><p>Size of the arrays per layers storing the radial dependent stress
 quantities (stresses, strain displacement etc..)</p></td>
    </tr>
    
    <tr>
        <td>i_tf_bucking</td>
        <td>Input</td>
        <td>integer</td>
        <td>-1</td>
        <td><p>Switch for TF inboard suport structure design:</p>
<p>Default setting for backward compatibility
     - if copper resistive TF (i_tf_sup = 0) : Free standing TF without bucking structure
     - if Superconducting TF  (i_tf_sup = 1) : Free standing TF with a steel casing
     - if aluminium  TF       (i_tf_sup = 2) : Free standing TF with a bucking structure
     Rem : the case is a bucking structure
 - =0 : Free standing TF without case/bucking cyliner (only a conductor layer)
 - =1 : Free standing TF with a case/bucking cylinder made of
     - if copper resistive     TF (i_tf_sup = 0) : used defined bucking cylinder
     - if Superconducting      TF (i_tf_sup = 1) : Steel casing
     - if aluminium resisitive TF (i_tf_sup = 2) : used defined bucking cylinder
 - =2 : The TF is in contact with the CS : "bucked and wedged design"
       Fast version : thin TF-CS interface neglected in the stress calculations (3 layers)
                      The CS is frictionally decoupled from the TF, does not carry axial tension
 - =3 : The TF is in contact with the CS : "bucked and wedged design"
       Full version : thin TF-CS Kapton interface introduced in the stress calculations (4 layers)
                      The CS and kaptop are frictionally decoupled from the TF, do not carry
                      axial tension</p></td>
    </tr>
    
    <tr>
        <td>n_tf_graded_layers</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Number of layers of different stress properties in the WP. If <code>n_tf_graded_layers &gt; 1</code>,
 a graded coil is condidered</p></td>
    </tr>
    
    <tr>
        <td>n_tf_stress_layers</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Number of layers considered for the inboard TF stress calculations
 set in initial.f90 from i_tf_bucking and n_tf_graded_layers</p></td>
    </tr>
    
    <tr>
        <td>n_tf_wp_layers</td>
        <td>Input</td>
        <td>integer</td>
        <td>5</td>
        <td><p>Maximum number of layers that can be considered in the TF coil composited/smeared
 stress analysis. This is the layers of one turn, not the entire WP.
 Default: 5. void, conductor, copper, conduit, insulation.</p></td>
    </tr>
    
    <tr>
        <td>j_tf_bus</td>
        <td>Input</td>
        <td>real</td>
        <td>1250000.0</td>
        <td><p>bussing current density (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>j_crit_str_tf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>j_crit_str : superconductor strand critical current density under operating
 conditions (A/m2). Necessary for the cost calculation in $/kAm</p></td>
    </tr>
    
    <tr>
        <td>j_crit_str_0</td>
        <td>Input</td>
        <td>real</td>
        <td>[5.96905476e+08 1.92550153e+09 7.24544683e+08 5.49858624e+08
 6.69284510e+08 0.00000000e+00 8.98964415e+08 1.15875300e+09
 8.65652123e+08]</td>
        <td><p>j_crit_str_pf_0 : superconductor strand critical current density at 6 T and 4.2 K (A/m2)
 Necessary for the cost calculation in $/kAm</p></td>
    </tr>
    
    <tr>
        <td>jwdgcrt</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>critical current density for winding pack (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>jwdgpro</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>allowable TF coil winding pack current density, for dump temperature rise protection (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>jwptf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>winding pack engineering current density (A/m2)</p></td>
    </tr>
    
    <tr>
        <td>oacdcp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Overall current density in TF coil inboard legs midplane (A/m2)
 Rem SK : Not used in tfcoil to set the current any more. Should not be used as
 iteration variable 12 any more. It is now calculated.</p></td>
    </tr>
    
    <tr>
        <td>eyoung_ins</td>
        <td>Input</td>
        <td>real</td>
        <td>100000000.0</td>
        <td><p>Insulator Young's modulus [Pa]. Default value (1.0D8) setup the following values
  - SC TF, eyoung_ins = 20 Gpa (default value from DDD11-2 v2 2 (2009))
  - Al TF, eyoung_ins = 2.5 GPa (Kapton polymer)</p></td>
    </tr>
    
    <tr>
        <td>eyoung_steel</td>
        <td>Input</td>
        <td>real</td>
        <td>205000000000.0</td>
        <td><p>Steel case Young's modulus (Pa) (default value from DDD11-2 v2 2 (2009))</p></td>
    </tr>
    
    <tr>
        <td>eyoung_cond_axial</td>
        <td>Input</td>
        <td>real</td>
        <td>660000000.0</td>
        <td><p>SC TF coil conductor Young's modulus in the parallel (along the wire/tape)
 direction [Pa]
 Set by user input only if <code>i_tf_cond_eyoung_axial == 1</code>; otherwise
 set by the behavior of that switch.</p></td>
    </tr>
    
    <tr>
        <td>eyoung_cond_trans</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>SC TF coil conductor Young's modulus in the transverse direction [Pa]
 Set by user input only if <code>i_tf_cond_eyoung_axial == 1</code>; otherwise
 set by the behavior of that switch.</p></td>
    </tr>
    
    <tr>
        <td>eyoung_res_tf_buck</td>
        <td>Input</td>
        <td>real</td>
        <td>150000000000.0</td>
        <td><p>Resistive TF magnets bucking cylinder young modulus (Pa)</p></td>
    </tr>
    
    <tr>
        <td>eyoung_copper</td>
        <td>Input</td>
        <td>real</td>
        <td>117000000000.0</td>
        <td><p>Copper young modulus. Default value taken from wikipedia</p></td>
    </tr>
    
    <tr>
        <td>eyoung_al</td>
        <td>Input</td>
        <td>real</td>
        <td>69000000000.0</td>
        <td><p>Aluminium young modulus.  Default value taken from wikipedia</p></td>
    </tr>
    
    <tr>
        <td>poisson_steel</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>Steel Poisson's ratio, Source : https://www.engineeringtoolbox.com/metals-poissons-ratio-d_1268.html</p></td>
    </tr>
    
    <tr>
        <td>poisson_copper</td>
        <td>Input</td>
        <td>real</td>
        <td>0.35</td>
        <td><p>Copper Poisson's ratio. Source : https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html</p></td>
    </tr>
    
    <tr>
        <td>poisson_al</td>
        <td>Input</td>
        <td>real</td>
        <td>0.35</td>
        <td><p>Aluminium Poisson's ratio.
 Source : https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html</p></td>
    </tr>
    
    <tr>
        <td>poisson_ins</td>
        <td>Input</td>
        <td>real</td>
        <td>0.34</td>
        <td><p>Insulation Poisson's ratio. Default: Kapton.
 Source : DuPont™ Kapton® HN datasheet.</p></td>
    </tr>
    
    <tr>
        <td>poisson_cond_axial</td>
        <td>Input</td>
        <td>real</td>
        <td>0.30000001192092896</td>
        <td><p>SC TF coil conductor Poisson's ratio in the parallel-transverse direction</p></td>
    </tr>
    
    <tr>
        <td>poisson_cond_trans</td>
        <td>Input</td>
        <td>real</td>
        <td>0.30000001192092896</td>
        <td><p>SC TF coil conductor Poisson's ratio in the transverse-transverse direction</p></td>
    </tr>
    
    <tr>
        <td>rbmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Radius of maximum TF B-field (m)</p></td>
    </tr>
    
    <tr>
        <td>res_tf_leg</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil leg resistance (ohm)</p></td>
    </tr>
    
    <tr>
        <td>toroidalgap</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>Minimal distance between two toroidal coils. (m)</p></td>
    </tr>
    
    <tr>
        <td>ftoroidalgap</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>F-value for minimum tftort (<code>constraint equation 82</code>)</p></td>
    </tr>
    
    <tr>
        <td>ripmax</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>aximum allowable toroidal field ripple amplitude at plasma edge (%)</p></td>
    </tr>
    
    <tr>
        <td>ripple</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>peak/average toroidal field ripple at plasma edge (%)</p></td>
    </tr>
    
    <tr>
        <td>c_tf_total</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total (summed) current in TF coils (A)</p></td>
    </tr>
    
    <tr>
        <td>n_radial_array</td>
        <td>Parameter</td>
        <td>integer</td>
        <td>50</td>
        <td><p>Size of the radial distribution arrays per layers
 used for stress, strain and displacement distibution</p></td>
    </tr>
    
    <tr>
        <td>radial_array</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Array refining the radii of the stress calculations arrays</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_r</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF Inboard leg radial stress in steel r distribution at mid-plane [Pa]</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_t</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF Inboard leg tangential stress in steel r distribution at mid-plane [Pa]</p></td>
    </tr>
    
    <tr>
        <td>deflect</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil radial deflection (displacement) radial distribution [m]</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_z</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF Inboard leg vertical tensile stress in steel at mid-plane [Pa]</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_vmises</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF Inboard leg Von-Mises stress in steel r distribution at mid-plane [Pa]</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_tresca</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF Inboard leg maximum shear stress (Tresca criterion) in steel r distribution at mid-plane [Pa]</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_cs_bucked</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Maximum shear stress (Tresca criterion) in CS structures at CS flux swing [Pa]:</p>
<ul>
<li>If superconducting CS (ipfres = 0): turn steel conduits stress</li>
<li>If resistive       CS (ipfres = 1): copper conductor stress</li>
</ul>
<p>Quantity only computed for bucked and wedged design (<code>i_tf_bucking &gt;= 2</code>)
 Def : CS Flux swing, instant when the current changes sign in CS (null current)</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_case</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Maximum shear stress (Tresca criterion) in TF casing steel structures (Pa)</p></td>
    </tr>
    
    <tr>
        <td>sig_tf_wp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td></td>
    </tr>
    
    <tr>
        <td>str_cs_con_res</td>
        <td>Input</td>
        <td>real</td>
        <td>-0.005</td>
        <td><p>Residual manufacturing strain in CS superconductor material</p></td>
    </tr>
    
    <tr>
        <td>str_pf_con_res</td>
        <td>Input</td>
        <td>real</td>
        <td>-0.005</td>
        <td><p>Residual manufacturing strain in PF superconductor material</p></td>
    </tr>
    
    <tr>
        <td>str_tf_con_res</td>
        <td>Input</td>
        <td>real</td>
        <td>-0.005</td>
        <td><p>Residual manufacturing strain in TF superconductor material
 If <code>i_str_wp == 0</code>, used to compute the critical surface.
 Otherwise, the self-consistent winding pack <code>str_wp</code> is used.</p></td>
    </tr>
    
    <tr>
        <td>str_wp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Axial (vertical) strain in the TF coil winding pack found by
 self-consistent stress/strain calculation.
 if <code>i_str_wp == 1</code>, used to compute the critical surface.
 Otherwise, the input value <code>str_tf_con_res</code> is used.
 Constrain the absolute value using <code>constraint equation 88</code>
 You can't have constraint 88 and i_str_wp = 0 at the same time</p></td>
    </tr>
    
    <tr>
        <td>str_wp_max</td>
        <td>Input</td>
        <td>real</td>
        <td>0.007</td>
        <td><p>Maximum allowed absolute value of the strain in the TF coil
 (<code>Constraint equation 88</code>)</p></td>
    </tr>
    
    <tr>
        <td>i_str_wp</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>Switch for the behavior of the TF strain used to compute
 the strain-dependent critical surface:</p>
<ul>
<li>=0  str_tf_con_res is used</li>
<li>=1  str_wp is used</li>
</ul></td>
    </tr>
    
    <tr>
        <td>quench_model</td>
        <td>Input</td>
        <td>character</td>
        <td>b'exponential '</td>
        <td><p>switch for TF coil quench model (Only applies to REBCO magnet at present, issue #522):</p>
<ul>
<li>='exponential' exponential quench with constant discharge resistor</li>
<li>='linear' quench with constant voltage</li>
</ul></td>
    </tr>
    
    <tr>
        <td>time1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Time at which TF quench is detected (s)</p></td>
    </tr>
    
    <tr>
        <td>tcritsc</td>
        <td>Input</td>
        <td>real</td>
        <td>16.0</td>
        <td><p>critical temperature (K) for superconductor at zero field and strain (<code>i_tf_sc_mat=4, =tc0m</code>)</p></td>
    </tr>
    
    <tr>
        <td>tdmptf</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>fast discharge time for TF coil in event of quench (s) (<code>iteration variable 56</code>)</p>
<p>For REBCO model, meaning depends on quench_model:</p>
<ul>
<li>exponential quench : e-folding time (s)`</li>
<li>linear quench : discharge time (s)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>tfareain</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Area of inboard midplane TF legs (m2)</p></td>
    </tr>
    
    <tr>
        <td>len_tf_bus</td>
        <td>Input</td>
        <td>real</td>
        <td>300.0</td>
        <td><p>TF coil bus length (m)</p></td>
    </tr>
    
    <tr>
        <td>m_tf_bus</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil bus mass (kg)</p></td>
    </tr>
    
    <tr>
        <td>tfckw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>available DC power for charging the TF coils (kW)</p></td>
    </tr>
    
    <tr>
        <td>tfcmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Peak power per TF power supply (MW)</p></td>
    </tr>
    
    <tr>
        <td>tfcpmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Peak resistive TF coil inboard leg power (MW)</p></td>
    </tr>
    
    <tr>
        <td>tfjtsmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF joints resistive power losses (MW)</p></td>
    </tr>
    
    <tr>
        <td>tfcryoarea</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>surface area of toroidal shells covering TF coils (m2)</p></td>
    </tr>
    
    <tr>
        <td>tficrn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil half-width - inner dr_bore (m)</p></td>
    </tr>
    
    <tr>
        <td>tfind</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil inductance (H)</p></td>
    </tr>
    
    <tr>
        <td>tfinsgap</td>
        <td>Input</td>
        <td>real</td>
        <td>0.01</td>
        <td><p>TF coil WP insertion gap (m)</p></td>
    </tr>
    
    <tr>
        <td>tflegmw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil outboard leg resistive power (MW)</p></td>
    </tr>
    
    <tr>
        <td>rho_cp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil inboard leg resistivity [Ohm-m]. If <code>itart=0</code>, this variable is the
 average resistivity over the whole magnet</p></td>
    </tr>
    
    <tr>
        <td>rho_tf_leg</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Resistivity of a TF coil leg (Ohm-m)</p></td>
    </tr>
    
    <tr>
        <td>rho_tf_bus</td>
        <td>Input</td>
        <td>real</td>
        <td>1.86e-08</td>
        <td><p>Resistivity of a TF coil bus (Ohm-m). Default values is for that of GLIDCOP AL-15 (C15715) at 293K</p></td>
    </tr>
    
    <tr>
        <td>frhocp</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>Centrepost resistivity enhancement factor. For <code>itart=0</code>, this factor
 is used for the whole magnet</p></td>
    </tr>
    
    <tr>
        <td>frholeg</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>Ouboard legs resistivity enhancement factor. Only used for <code>itart=1</code>.</p></td>
    </tr>
    
    <tr>
        <td>i_cp_joints</td>
        <td>Input</td>
        <td>integer</td>
        <td>-1</td>
        <td><p>Switch for CP demoutable joints type
  -= 0 : Clampled joints
  -= 1 : Sliding joints
 Default value (-1) choses :
   Sliding joints for resistive magnets (i_tf_sup = 0, 2)
   Clampled joints for superconducting magents (i_tf_sup = 1)</p></td>
    </tr>
    
    <tr>
        <td>rho_tf_joints</td>
        <td>Input</td>
        <td>real</td>
        <td>2.5e-10</td>
        <td><p>TF joints surfacic resistivity [ohm.m]. Feldmetal joints assumed.</p></td>
    </tr>
    
    <tr>
        <td>n_tf_joints_contact</td>
        <td>Input</td>
        <td>integer</td>
        <td>6</td>
        <td><p>Number of contact per turn</p></td>
    </tr>
    
    <tr>
        <td>n_tf_joints</td>
        <td>Input</td>
        <td>integer</td>
        <td>4</td>
        <td><p>Number of joints
 Ex: n_tf_joints = 2 for top and bottom CP joints</p></td>
    </tr>
    
    <tr>
        <td>th_joint_contact</td>
        <td>Input</td>
        <td>real</td>
        <td>0.03</td>
        <td><p>TF sliding joints contact pad width [m]</p></td>
    </tr>
    
    <tr>
        <td>pres_joints</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Calculated TF joints resistive power losses [W]</p></td>
    </tr>
    
    <tr>
        <td>len_tf_coil</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil circumference (m)</p></td>
    </tr>
    
    <tr>
        <td>eff_tf_cryo</td>
        <td>Input</td>
        <td>real</td>
        <td>-1.0</td>
        <td><p>TF cryoplant efficiency (compared to pefect Carnot cycle).
 Using -1 set the default value depending on magnet technology:</p>
<ul>
<li>i_tf_sup = 1 : SC magnet, eff_tf_cryo = 0.13 (ITER design)</li>
<li>i_tf_sup = 2 : Cryo-aluminium, eff_tf_cryo = 0.4</li>
</ul></td>
    </tr>
    
    <tr>
        <td>n_tf_coils</td>
        <td>Input</td>
        <td>real</td>
        <td>16.0</td>
        <td><p>Number of TF coils (default = 50 for stellarators). Number of TF coils outer legs for ST</p></td>
    </tr>
    
    <tr>
        <td>tfocrn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil half-width - outer dr_bore (m)</p></td>
    </tr>
    
    <tr>
        <td>tfsai</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>area of the inboard TF coil legs (m2)</p></td>
    </tr>
    
    <tr>
        <td>tfsao</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>area of the outboard TF coil legs (m2)</p></td>
    </tr>
    
    <tr>
        <td>tftmp</td>
        <td>Input</td>
        <td>real</td>
        <td>4.5</td>
        <td><p>peak helium coolant temperature in TF coils and PF coils (K)</p></td>
    </tr>
    
    <tr>
        <td>tftort</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>TF coil toroidal thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>thicndut</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0008</td>
        <td><p>conduit insulation thickness (m)</p></td>
    </tr>
    
    <tr>
        <td>layer_ins</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Additional insulation thickness between layers (m)</p></td>
    </tr>
    
    <tr>
        <td>thkcas</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>inboard TF coil case outer (non-plasma side) thickness (m) (<code>iteration variable 57</code>)
 (calculated for stellarators)</p></td>
    </tr>
    
    <tr>
        <td>dr_tf_wp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radial thickness of winding pack (m) (<code>iteration variable 140</code>) (issue #514)</p></td>
    </tr>
    
    <tr>
        <td>thwcndut</td>
        <td>Input</td>
        <td>real</td>
        <td>0.008</td>
        <td><p>TF coil conduit case thickness (m) (<code>iteration variable 58</code>)</p></td>
    </tr>
    
    <tr>
        <td>tinstf</td>
        <td>Input</td>
        <td>real</td>
        <td>0.018</td>
        <td><p>Thickness of the ground insulation layer surrounding (m)
   - Superconductor TF (<code>i_tf_sup == 1</code>) : The TF coil Winding packs
   - Resistive magnets (<code>i_tf_sup /= 1</code>) : The TF coil wedges
 Rem : Thickness calculated for stellarators.</p></td>
    </tr>
    
    <tr>
        <td>tmargmin_tf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>minimum allowable temperature margin : TF coils (K)</p></td>
    </tr>
    
    <tr>
        <td>tmargmin_cs</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>minimum allowable temperature margin : CS (K)</p></td>
    </tr>
    
    <tr>
        <td>tmargmin</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>minimum allowable temperature margin : TFC AND CS (K)</p></td>
    </tr>
    
    <tr>
        <td>temp_margin</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>temperature margin (K)</p></td>
    </tr>
    
    <tr>
        <td>tmargtf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil temperature margin (K)</p></td>
    </tr>
    
    <tr>
        <td>tmaxpro</td>
        <td>Input</td>
        <td>real</td>
        <td>150.0</td>
        <td><p>maximum temp rise during a quench for protection (K)</p></td>
    </tr>
    
    <tr>
        <td>tmax_croco</td>
        <td>Input</td>
        <td>real</td>
        <td>200.0</td>
        <td><p>CroCo strand: maximum permitted temp during a quench (K)</p></td>
    </tr>
    
    <tr>
        <td>croco_quench_temperature</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>CroCo strand: Actual temp reached during a quench (K)</p></td>
    </tr>
    
    <tr>
        <td>tmpcry</td>
        <td>Input</td>
        <td>real</td>
        <td>4.5</td>
        <td><p>coil temperature for cryogenic plant power calculation (K)</p></td>
    </tr>
    
    <tr>
        <td>n_tf_turn</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>number of turns per TF coil</p></td>
    </tr>
    
    <tr>
        <td>vdalw</td>
        <td>Input</td>
        <td>real</td>
        <td>20.0</td>
        <td><p>max voltage across TF coil during quench (kV) (<code>iteration variable 52</code>)</p></td>
    </tr>
    
    <tr>
        <td>vforce</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vertical tension on inboard leg/coil (N)</p></td>
    </tr>
    
    <tr>
        <td>f_vforce_inboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.5</td>
        <td><p>Fraction of the total vertical force taken by the TF inboard leg tension
 Not used for resistive <code>itart=1</code> (sliding joints)</p></td>
    </tr>
    
    <tr>
        <td>vforce_outboard</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Vertical tension on outboard leg/coil (N)</p></td>
    </tr>
    
    <tr>
        <td>vftf</td>
        <td>Input</td>
        <td>real</td>
        <td>0.4</td>
        <td><p>coolant fraction of TFC 'cable' (<code>i_tf_sup=1</code>), or of TFC leg (<code>i_tf_ssup=0</code>)</p></td>
    </tr>
    
    <tr>
        <td>voltfleg</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume of each TF coil outboard leg (m3)</p></td>
    </tr>
    
    <tr>
        <td>vtfkv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil voltage for resistive coil including bus (kV)</p></td>
    </tr>
    
    <tr>
        <td>vtfskv</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>voltage across a TF coil during quench (kV)</p></td>
    </tr>
    
    <tr>
        <td>whtcas</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass per coil of external case (kg)</p></td>
    </tr>
    
    <tr>
        <td>whtcon</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>TF coil conductor mass per coil (kg/coil).
 For <code>itart=1</code>, coil is return limb plus centrepost/n_tf_coils</p></td>
    </tr>
    
    <tr>
        <td>whtconcu</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>copper mass in TF coil conductor (kg/coil).
 For <code>itart=1</code>, coil is return limb plus centrepost/n_tf_coils</p></td>
    </tr>
    
    <tr>
        <td>whtconal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Aluminium mass in TF coil conductor (kg/coil).
 For <code>itart=1</code>, coil is return limb plus centrepost/n_tf_coils</p></td>
    </tr>
    
    <tr>
        <td>whtconin</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>conduit insulation mass in TF coil conductor (kg/coil)</p></td>
    </tr>
    
    <tr>
        <td>whtconsc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>superconductor mass in TF coil cable (kg/coil)</p></td>
    </tr>
    
    <tr>
        <td>whtconsh</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>steel conduit mass in TF coil conductor (kg/coil)</p></td>
    </tr>
    
    <tr>
        <td>whtgw</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of ground-wall insulation layer per coil (kg/coil)</p></td>
    </tr>
    
    <tr>
        <td>whttf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total mass of the TF coils (kg)</p></td>
    </tr>
    
    <tr>
        <td>wwp1</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>width of first step of winding pack (m)</p></td>
    </tr>
    
    <tr>
        <td>wwp2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>width of second step of winding pack (m)</p></td>
    </tr>
    
    <tr>
        <td>dthet</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>angle of arc i (rad)</p></td>
    </tr>
    
    <tr>
        <td>radctf</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>radius of arc i (m)</p></td>
    </tr>
    
    <tr>
        <td>xarc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>x location of arc point i on surface (m)</p></td>
    </tr>
    
    <tr>
        <td>xctfc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>x location of arc centre i (m)</p></td>
    </tr>
    
    <tr>
        <td>yarc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>y location of arc point i on surface (m)</p></td>
    </tr>
    
    <tr>
        <td>yctfc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>y location of arc centre i (m)</p></td>
    </tr>
    
    <tr>
        <td>tfa</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Horizontal radius of inside edge of TF coil (m)</p></td>
    </tr>
    
    <tr>
        <td>tfb</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Vertical radius of inside edge of TF coil (m)</p></td>
    </tr>
    
    <tr>
        <td>drtop</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>centrepost taper maximum radius adjustment (m)</p></td>
    </tr>
    
    <tr>
        <td>dztop</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>centrepost taper height adjustment (m)</p></td>
    </tr>
    
    <tr>
        <td>etapump</td>
        <td>Input</td>
        <td>real</td>
        <td>0.8</td>
        <td><p>centrepost coolant pump efficiency</p></td>
    </tr>
    
    <tr>
        <td>fcoolcp</td>
        <td>Input</td>
        <td>real</td>
        <td>0.3</td>
        <td><p>coolant fraction of TF coil inboard legs (<code>iteration variable 23</code>)</p></td>
    </tr>
    
    <tr>
        <td>f_a_tf_cool_outboard</td>
        <td>Input</td>
        <td>real</td>
        <td>0.2</td>
        <td><p>coolant fraction of TF coil outboard legs</p></td>
    </tr>
    
    <tr>
        <td>a_cp_cool</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Centrepost cooling area toroidal cross-section (constant over the whole CP)</p></td>
    </tr>
    
    <tr>
        <td>ncool</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>number of centrepost coolant tubes</p></td>
    </tr>
    
    <tr>
        <td>ppump</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>centrepost coolant pump power (W)</p></td>
    </tr>
    
    <tr>
        <td>p_cp_resistive</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>resistive power in the centrepost (itart=1) [W].
 If <code>itart=0</code>, this variable is the ressitive power on the whole magnet</p></td>
    </tr>
    
    <tr>
        <td>p_tf_leg_resistive</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Summed resistive power in the TF coil legs [W]. Remain 0 if <code>itart=0</code>.</p></td>
    </tr>
    
    <tr>
        <td>ptempalw</td>
        <td>Input</td>
        <td>real</td>
        <td>473.15</td>
        <td><p>maximum peak centrepost temperature (K) (<code>constraint equation 44</code>)</p></td>
    </tr>
    
    <tr>
        <td>rcool</td>
        <td>Input</td>
        <td>real</td>
        <td>0.005</td>
        <td><p>average radius of coolant channel (m) (<code>iteration variable 69</code>)</p></td>
    </tr>
    
    <tr>
        <td>tcoolin</td>
        <td>Input</td>
        <td>real</td>
        <td>313.15</td>
        <td><p>centrepost coolant inlet temperature (K)</p></td>
    </tr>
    
    <tr>
        <td>dtiocool</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>inlet / outlet TF coil coolant temperature rise (K)</p></td>
    </tr>
    
    <tr>
        <td>temp_cp_average</td>
        <td>Input</td>
        <td>real</td>
        <td>373.15</td>
        <td><p>Average temperature of centrepost called CP (K). Only used for resistive coils
 to compute the resisitive heating. Must be an iteration variable for
 ST (<code>itart=1</code>) (<code>iteration variable 20</code>)</p></td>
    </tr>
    
    <tr>
        <td>tcpav2</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Computed centrepost average temperature (K) (for consistency)</p></td>
    </tr>
    
    <tr>
        <td>temp_tf_legs_outboard</td>
        <td>Input</td>
        <td>real</td>
        <td>-1.0</td>
        <td><p>Average temperature of the TF outboard legs [K]. If <code>temp_tf_legs_outboard=-1.0</code>, the ouboard
 legs and CP temperatures are the same. Fixed for now, should use a contraints eq like temp_cp_average</p></td>
    </tr>
    
    <tr>
        <td>tcpmax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>peak centrepost temperature (K)</p></td>
    </tr>
    
    <tr>
        <td>vcool</td>
        <td>Input</td>
        <td>real</td>
        <td>20.0</td>
        <td><p>inlet centrepost coolant flow speed at midplane (m/s) (<code>iteration variable 70</code>)</p></td>
    </tr>
    
    <tr>
        <td>vol_cond_cp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Exact conductor volume in the centrepost (m3)</p></td>
    </tr>
    
    <tr>
        <td>whtcp</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of TF coil inboard legs (kg)</p></td>
    </tr>
    
    <tr>
        <td>whttflgs</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of the TF coil legs (kg)</p></td>
    </tr>
    
    <tr>
        <td>cryo_cool_req</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>Cryo cooling requirement at helium temp 4.5K (kW)</p></td>
    </tr>
    
    <tr>
        <td>theta1_coil</td>
        <td>Input</td>
        <td>real</td>
        <td>45.0</td>
        <td><p>The angle of the outboard arc forming the TF coil current center line [deg]</p></td>
    </tr>
    
    <tr>
        <td>theta1_vv</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>The angle of the outboard arc forming the Vacuum Vessel current center line [deg]</p></td>
    </tr>
    
    <tr>
        <td>max_vv_stress</td>
        <td>Input</td>
        <td>real</td>
        <td>143000000.0</td>
        <td><p>The allowable peak maximum shear stress in the vacuum vessel due to quench and fast discharge of the TF coils [Pa]</p></td>
    </tr>
    
</table>

## times_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>pulsetimings</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>Switch for pulse timings (if i_pulsed_plant=1):</p>
<ul>
<li>=0, t_current_ramp_up = Ip(MA)/0.1 t_precharge, t_ramp_down = input</li>
<li>=1, t_current_ramp_up = iteration var or input. t_precharge/t_ramp_down max of input or t_current_ramp_up</li>
</ul></td>
    </tr>
    
    <tr>
        <td>t_burn</td>
        <td>Input</td>
        <td>real</td>
        <td>1000.0</td>
        <td><p>flat-top duration (s) (calculated if <code>i_pulsed_plant=1</code>)</p></td>
    </tr>
    
    <tr>
        <td>t_burn_0</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>burn time (s) - used for internal consistency</p></td>
    </tr>
    
    <tr>
        <td>t_cycle</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>full cycle time (s)</p></td>
    </tr>
    
    <tr>
        <td>tdown</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>down time (s)</p></td>
    </tr>
    
    <tr>
        <td>t_between_pulse</td>
        <td>Input</td>
        <td>real</td>
        <td>1800.0</td>
        <td><p>time between pulses in a pulsed reactor (s) (<code>iteration variable 17</code>)</p></td>
    </tr>
    
    <tr>
        <td>t_fusion_ramp</td>
        <td>Input</td>
        <td>real</td>
        <td>10.0</td>
        <td><p>time for plasma temperature and density rise to full values (s)</p></td>
    </tr>
    
    <tr>
        <td>tim</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>array of time points during plasma pulse (s)</p></td>
    </tr>
    
    <tr>
        <td>timelabel</td>
        <td>Input</td>
        <td>character</td>
        <td>[b'Start      ' b'BOP        ' b'EOR        ' b'BOF        '
 b'EOF        ' b'EOP        ']</td>
        <td><p>array of time labels during plasma pulse (s)</p></td>
    </tr>
    
    <tr>
        <td>intervallabel</td>
        <td>Input</td>
        <td>character</td>
        <td>[b't_precharge        ' b't_current_ramp_up  ' b't_fusion_ramp      '
 b't_burn             ' b't_ramp_down        ']</td>
        <td><p>time intervals - as strings (s)</p></td>
    </tr>
    
    <tr>
        <td>t_current_ramp_up</td>
        <td>Input</td>
        <td>real</td>
        <td>30.0</td>
        <td><p>time for plasma current to ramp up to approx. full value (s) (calculated if <code>i_pulsed_plant=0</code>)
 (<code>iteration variable 65</code>)</p></td>
    </tr>
    
    <tr>
        <td>i_t_current_ramp_up</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>Switch for plasma current ramp-up time (if i_pulsed_plant=0):</p>
<ul>
<li>= 0, t_current_ramp_up = t_precharge = t_ramp_down = Ip(MA)/0.5</li>
<li>= 1, t_current_ramp_up, t_precharge, t_ramp_down are input</li>
</ul></td>
    </tr>
    
    <tr>
        <td>t_pulse_repetition</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>pulse length = t_current_ramp_up + t_fusion_ramp + t_burn + t_ramp_down</p></td>
    </tr>
    
    <tr>
        <td>t_ramp_down</td>
        <td>Input</td>
        <td>real</td>
        <td>15.0</td>
        <td><p>time for plasma current, density, and temperature to ramp down to zero, simultaneously (s); if pulsed, = t_current_ramp_up
 the CS and PF coil currents also ramp to zero at the same time</p></td>
    </tr>
    
    <tr>
        <td>t_precharge</td>
        <td>Input</td>
        <td>real</td>
        <td>15.0</td>
        <td><p>the time for the central solenoid and PF coils to ramp from zero to max current (s); if pulsed, = t_current_ramp_up</p></td>
    </tr>
    
</table>

## vacuum_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>vacuum_model</td>
        <td>Input</td>
        <td>character</td>
        <td>b'old   '</td>
        <td><p>switch for vacuum pumping model:</p>
<ul>
<li>='old' for old detailed ETR model</li>
<li>='simple' for simple steady-state model with comparison to ITER cryopumps</li>
</ul></td>
    </tr>
    
    <tr>
        <td>niterpump</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>number of high vacuum pumps (real number), each with the throughput of one
 ITER cryopump (50 Pa m3 s-1), all operating at the same time (<code>vacuum_model='simple'</code>)</p></td>
    </tr>
    
    <tr>
        <td>ntype</td>
        <td>Input</td>
        <td>integer</td>
        <td>1</td>
        <td><p>switch for vacuum pump type:</p>
<ul>
<li>=0 - for turbomolecular pump (magnetic bearing) with speed of 2.0 m3/s
   (1.95 for N2, 1.8 for He, 1.8 for DT)</li>
<li>=1 - for compound cryopump with nominal speed of 10.0 m3/s
   (9.0 for N2, 5.0 for He and 25.0 for DT)</li>
</ul></td>
    </tr>
    
    <tr>
        <td>nvduct</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>number of ducts (torus to pumps)</p></td>
    </tr>
    
    <tr>
        <td>dlscal</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>vacuum system duct length scaling</p></td>
    </tr>
    
    <tr>
        <td>pbase</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0005</td>
        <td><p>base pressure during dwell before gas pre-fill(Pa)</p></td>
    </tr>
    
    <tr>
        <td>prdiv</td>
        <td>Input</td>
        <td>real</td>
        <td>0.36</td>
        <td><p>divertor chamber pressure during burn (Pa)</p></td>
    </tr>
    
    <tr>
        <td>pumptp</td>
        <td>Input</td>
        <td>real</td>
        <td>1.2155e+22</td>
        <td><p>Pump throughput (molecules/s) (default is ITER value)</p></td>
    </tr>
    
    <tr>
        <td>rat</td>
        <td>Input</td>
        <td>real</td>
        <td>1.3e-08</td>
        <td><p>plasma chamber wall outgassing rate (Pa-m/s)</p></td>
    </tr>
    
    <tr>
        <td>tn</td>
        <td>Input</td>
        <td>real</td>
        <td>300.0</td>
        <td><p>neutral gas temperature in chamber (K)</p></td>
    </tr>
    
    <tr>
        <td>vacdshm</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>mass of vacuum duct shield (kg)</p></td>
    </tr>
    
    <tr>
        <td>vcdimax</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>diameter of duct passage (m)</p></td>
    </tr>
    
    <tr>
        <td>vpumpn</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>number of high vacuum pumps</p></td>
    </tr>
    
    <tr>
        <td>dwell_pump</td>
        <td>Output</td>
        <td>integer</td>
        <td>-</td>
        <td><p>switch for dwell pumping options:</p>
<ul>
<li>=0 pumping only during t_between_pulse</li>
<li>=1 pumping only during t_precharge</li>
<li>=2 pumping during t_between_pulse + t_precharge</li>
</ul></td>
    </tr>
    
    <tr>
        <td>pumpareafraction</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0203</td>
        <td><p>area of one pumping port as a fraction of plasma surface area</p></td>
    </tr>
    
    <tr>
        <td>pumpspeedmax</td>
        <td>Input</td>
        <td>real</td>
        <td>27.3</td>
        <td><p>maximum pumping speed per unit area for deuterium &amp; tritium, molecular flow</p></td>
    </tr>
    
    <tr>
        <td>pumpspeedfactor</td>
        <td>Input</td>
        <td>real</td>
        <td>0.167</td>
        <td><p>effective pumping speed reduction factor due to duct impedance</p></td>
    </tr>
    
    <tr>
        <td>initialpressure</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>initial neutral pressure at the beginning of the dwell phase (Pa)</p></td>
    </tr>
    
    <tr>
        <td>outgasindex</td>
        <td>Input</td>
        <td>real</td>
        <td>1.0</td>
        <td><p>outgassing decay index</p></td>
    </tr>
    
    <tr>
        <td>outgasfactor</td>
        <td>Input</td>
        <td>real</td>
        <td>0.0235</td>
        <td><p>outgassing prefactor kw: outgassing rate at 1 s per unit area (Pa m s-1)</p></td>
    </tr>
    
</table>

## water_usage_variables
<table class="vardes">
    <tr>
        <th>Name</th>
        <th>Type</th>
        <th>Datatype</th>
        <th>Default Value</th>
        <th>Description</th>
    </tr>
    
    <tr>
        <td>airtemp</td>
        <td>Input</td>
        <td>real</td>
        <td>15.0</td>
        <td><p>ambient air temperature (degrees Celsius)</p></td>
    </tr>
    
    <tr>
        <td>watertemp</td>
        <td>Input</td>
        <td>real</td>
        <td>5.0</td>
        <td><p>water temperature (degrees Celsius)</p></td>
    </tr>
    
    <tr>
        <td>windspeed</td>
        <td>Input</td>
        <td>real</td>
        <td>4.0</td>
        <td><p>wind speed (m/s)</p></td>
    </tr>
    
    <tr>
        <td>waterdens</td>
        <td>Input</td>
        <td>real</td>
        <td>998.02</td>
        <td><p>density of water (kg/m3)
   for simplicity, set to static value applicable to water at 21 degC</p></td>
    </tr>
    
    <tr>
        <td>latentheat</td>
        <td>Input</td>
        <td>real</td>
        <td>2257000.0</td>
        <td><p>latent heat of vaporization (J/kg)
   for simplicity, set to static value applicable at 1 atm (100 kPa) air pressure</p></td>
    </tr>
    
    <tr>
        <td>volheat</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volumetric heat of vaporization (J/m3)</p></td>
    </tr>
    
    <tr>
        <td>evapratio</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>evaporation ratio: ratio of the heat used to evaporate water
   to the total heat discharged through the tower</p></td>
    </tr>
    
    <tr>
        <td>evapvol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>evaporated volume of water (m3)</p></td>
    </tr>
    
    <tr>
        <td>energypervol</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>input waste (heat) energy cooled per evaporated volume (J/m3)</p></td>
    </tr>
    
    <tr>
        <td>volperenergy</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>volume evaporated by units of heat energy (m3/MJ)</p></td>
    </tr>
    
    <tr>
        <td>waterusetower</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total volume of water used in cooling tower (m3)</p></td>
    </tr>
    
    <tr>
        <td>wateruserecirc</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total volume of water used in recirculating system (m3)</p></td>
    </tr>
    
    <tr>
        <td>wateruseonethru</td>
        <td>Output</td>
        <td>real</td>
        <td>-</td>
        <td><p>total volume of water used in once-through system (m3)</p></td>
    </tr>
    
</table>
