"""Module containing routines to perform a parameter scan

None
This module contains routines to perform a parameter scan
over a range of values of a particular scanning variable.
"""

import numpy as np

IPNSCNS: int = 1000
"""Maximum number of scan points"""


IPNSCNV: int = 81
"""Number of available scan variables"""


NOUTVARS: int = 84


WIDTH: int = 110


scan_dim: int = None
"""1-D or 2-D scan switch (1=1D, 2=2D)"""


isweep: int = None
"""Number of scan points to calculate"""


isweep_2: int = None
"""Number of 2D scan points to calculate"""


nsweep: int = None
"""Switch denoting quantity to scan:<UL>
<LI> 1  aspect
<LI> 2  pflux_div_heat_load_max_mw
<LI> 3  p_plant_electric_net_required_mw
<LI> 4  hfact
<LI> 5  oacdcp
<LI> 6  pflux_fw_neutron_max_mw
<LI> 7  beamfus0
<LI> 8  NOT USED
<LI> 9  temp_plasma_electron_vol_avg_kev
<LI> 10 NOT USED
<LI> 11 beta_norm_max
<LI> 12 f_c_plasma_bootstrap_max
<LI> 13 boundu(10: hfact)
<LI> 14 fiooic
<LI> 16 rmajor
<LI> 15 NOT USED
<LI> 17 b_tf_inboard_max
<LI> 18 eta_cd_norm_hcd_primary_max
<LI> 19 boundl(16: dr_cs)
<LI> 20 t_burn_min
<LI> 21 NOT USED
<LI> 22 f_t_plant_available (N.B. requires i_plant_availability=0)
<LI> 23 NOT USED
<LI> 24 p_fusion_total_max_mw
<LI> 25 kappa
<LI> 26 triang
<LI> 27 tbrmin (for blktmodel > 0 only)
<LI> 28 b_plasma_toroidal_on_axis
<LI> 29 radius_plasma_core_norm
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
<LI> 42 Argon fraction f_nd_impurity_electrons(9)
<LI> 43 normalised minor radius at which electron cyclotron current drive is maximum
<LI> 44 Allowable maximum shear stress (Tresca) in tf coil structural material
<LI> 45 Minimum allowable temperature margin ; tf coils
<LI> 46 boundu(150) f_nd_plasma_separatrix_greenwald
<LI> 47 impurity_enrichment(9) Argon impurity enrichment
<LI> 48 TF coil - n_tf_wp_pancakes (integer turn winding pack)
<LI> 49 TF coil - n_tf_wp_layers (integer turn winding pack)
<LI> 50 Xenon fraction f_nd_impurity_electrons(13)
<LI> 51 Power fraction to lower DN Divertor f_p_div_lower
<LI> 52 SoL radiation fraction
<LI> 54 GL_nbti upper critical field at 0 Kelvin
<LI> 55 `dr_shld_inboard` : Inboard neutron shield thickness
<LI> 56 p_cryo_plant_electric_max_mw: Maximum cryogenic power (ixx=164, ixc=87)
<LI> 57 `b_plasma_toroidal_on_axis` lower boundary
<LI> 58 `dr_fw_plasma_gap_inboard` : Inboard plasma-first wall gap
<LI> 59 `dr_fw_plasma_gap_outboard` : Outboard plasma-first wall gap
<LI> 60 sig_tf_wp_max: Allowable stress in TF Coil conduit (Tresca)
<LI> 61 copperaoh_m2_max : CS coil current / copper area
<LI> 62 j_cs_flat_top_end : CS coil current density at EOF
<LI> 63 dr_cs : CS thickness (m)
<LI> 64 f_z_cs_tf_internal : CS height (m)
<LI> 65 n_cycle_min : Minimum cycles for CS stress model constraint 90
<LI> 66 f_a_cs_turn_steel: Steel fraction in CS coil
<LI> 67 t_crack_vertical: Initial crack vertical dimension (m) </UL>
<LI> 68 `inlet_temp_liq' : Inlet temperature of blanket liquid metal coolant/breeder (K)
<LI> 69 `outlet_temp_liq' : Outlet temperature of blanket liquid metal coolant/breeder (K)
<LI> 70 `blpressure_liq' : Blanket liquid metal breeder/coolant pressure (Pa)
<LI> 71 `n_liq_recirc' : Selected number of liquid metal breeder recirculations per day
<LI> 72 `bz_channel_conduct_liq' : Conductance of liquid metal breeder duct walls (A V-1 m-1)
<LI> 73 `pnuc_fw_ratio_dcll' : Ratio of FW nuclear power as fraction of total (FW+BB)
<LI> 74 `f_nuc_pow_bz_struct' : Fraction of BZ power cooled by primary coolant for dual-coolant balnket
<LI> 75 dx_fw_module : pitch of first wall cooling channels (m)
<LI> 76 eta_turbine : Thermal conversion eff.
<LI> 77 startupratio : Gyrotron redundancy
<LI> 78 fkind : Multiplier for Nth of a kind costs
<LI> 79 eta_ecrh_injector_wall_plug : ECH wall plug to injector efficiency
"""


nsweep_2: int = None
"""nsweep_2 /3/ : switch denoting quantity to scan for 2D scan:"""


sweep: list[float] = None
"""sweep(IPNSCNS) /../: actual values to use in scan"""


sweep_2: list[float] = None
"""sweep_2(IPNSCNS) /../: actual values to use in 2D scan"""


# Vars in subroutines scan_1d and scan_2d requiring re-initialising before
# each new run

first_call_1d: bool = None

first_call_2d: bool = None


def init_scan_variables():
    """Initialise the scan module"""
    global \
        scan_dim, \
        isweep, \
        isweep_2, \
        nsweep, \
        nsweep_2, \
        sweep, \
        sweep_2, \
        first_call_1d, \
        first_call_2d

    scan_dim = 1
    isweep = 0
    isweep_2 = 0
    nsweep = 1
    nsweep_2 = 3
    sweep = np.zeros(1000, dtype=np.float64)
    sweep_2 = np.zeros(1000, dtype=np.float64)
    first_call_1d = True
    first_call_2d = True
