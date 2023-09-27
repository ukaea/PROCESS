from dataclasses import dataclass


@dataclass
class VariableMetadata:
    latex: str
    description: str


var_dicts = {
    "shldith": VariableMetadata(
        latex=r"$\Delta R_\mathrm{sh}$ [$m$]", description="Inboard shield"
    ),
    "rmajor": VariableMetadata(
        latex=r"$R_\mathrm{maj}$ [$m$]", description="Major radius"
    ),
    "crypmw": VariableMetadata(
        latex=r"$P_\mathrm{cryo}$ [$MW$]", description="Cryogenic power"
    ),
    "bt": VariableMetadata(
        latex=r"$B_\mathrm{T}$ [$T$]", description="Toroidal magnetic field"
    ),
    "tfcth": VariableMetadata(
        latex=r"$\Delta R_\mathrm{TF}$ [$m$]", description="TF coil thickness"
    ),
    "powfmw": VariableMetadata(
        latex=r"$P_\mathrm{fus}$ [$MW$]", description="Fusion power"
    ),
    "pinjemw": VariableMetadata(
        latex=r"$P_\mathrm{inj}$ [$MW$]", description="Injected power"
    ),
    "pnetelmw": VariableMetadata(
        latex=r"$P_\mathrm{Net\ elec}$ [$MW$]", description="Net electrical power"
    ),
    "taueff": VariableMetadata(
        latex=r"$\tau_\mathrm{E}$ [s]", description="Effective energy confinement time"
    ),
    "ralpne": VariableMetadata(
        latex=r"$f_\mathrm{\alpha}$", description="Alpha particle heating fraction"
    ),
    "te": VariableMetadata(
        latex=r"$\left< T_\mathrm{e} \right>$",
        description="Average electron temperature",
    ),
    "taulimit": VariableMetadata(
        latex=r"$max : \frac{\tau_\mathrm{\alpha}}{\tau_\mathrm{E}}$",
        description="Ratio of alpha heating time to energy confinement time",
    ),
    "scrapli": VariableMetadata(
        latex=r"$\Delta R_\mathrm{FW-sep}$ [$m$]", description="Scrape-off layer width"
    ),
    "scraplo": VariableMetadata(
        latex=r"$\Delta R_\mathrm{FW-sep}^\mathrm{out}$ [$m$]",
        description="Outer scrape-off layer width",
    ),
    "vforce": VariableMetadata(
        latex=r"$F_\mathrm{z}^\mathrm{in}$ [$N$]", description="Vertical force"
    ),
    "thkcas": VariableMetadata(
        latex=r"$\Delta R_\mathrm{TF}^\mathrm{buck}$ [$m$]",
        description="Bucking mode TF coil thickness",
    ),
    "bmaxtf": VariableMetadata(
        latex=r"$B_\mathrm{TF}^\mathrm{max}$ [$T$]",
        description="Maximum TF magnetic field",
    ),
    "ritfc": VariableMetadata(
        latex=r"$I_\mathrm{TF}^\mathrm{tot}$ [$A$]",
        description="Total TF coil current",
    ),
    "dr_tf_wp": VariableMetadata(
        latex=r"$\Delta R_\mathrm{TF}^\mathrm{WP}$ [$m$]",
        description="TF coil winding pack width",
    ),
    "aspect": VariableMetadata(latex=r"$A$", description="Aspect ratio"),
    "rminor": VariableMetadata(
        latex=r"$a_\mathrm{min}$ [$m$]", description="Minor radius"
    ),
    "capcost": VariableMetadata(
        latex=r"$C_\mathrm{cap}$ [$M\$ $]", description="Capital cost"
    ),
    "r_tf_outboard_mid": VariableMetadata(
        latex=r"$\Delta R_\mathrm{TF}^\mathrm{out\ mid}$ [$m$]",
        description="Outboard midplane TF coil thickness",
    ),
    "pgrossmw": VariableMetadata(
        latex=r"$P_\mathrm{gross}^\mathrm{elec}$ [$MW$]",
        description="Gross electrical power",
    ),
    "htpmw": VariableMetadata(
        latex=r"$P_\mathrm{Primary\ coolant}^\mathrm{elec}$ [$MW$]",
        description="Primary coolant power",
    ),
    "ppfmw": VariableMetadata(
        latex=r"$P_\mathrm{PF}^\mathrm{elec}$ [$MW$]",
        description="Poloidal field power",
    ),
    "hmax": VariableMetadata(
        latex=r"$z_\mathrm{TF}^\mathrm{pl\ side}$ [$m$]",
        description="Height of TF poloidal field side",
    ),
    "thicndut": VariableMetadata(
        latex=r"\Delta l_\mathrm{steel\ jacket}^\mathrm{turn}",
        description="Thickness of steel jacket turn",
    ),
    "cpttf": VariableMetadata(
        latex=r"$I_\mathrm{TF}^\mathrm{turn}$ [$A$]",
        description="TF turn current",
    ),
    "boundl(2)": VariableMetadata(
        latex=r"$B_\mathrm{T}^\mathrm{min}$ [$A$]",
        description="Minimum toroidal field",
    ),
    "pinjmw": VariableMetadata(
        latex=r"$P_\mathrm{inj}$ [$MW$]", description="Injected power"
    ),
    "hldivlim": VariableMetadata(
        latex=r"$q_\mathrm{div}^\mathrm{max}$ [$MW.m^{-2}$]",
        description="Maximum divertor heat load",
    ),
    "hfact": VariableMetadata(latex=r"$f_\mathrm{H}$", description="H-factor"),
    "kappa": VariableMetadata(
        latex=r"$\kappa_\mathrm{sep}$", description="Separation factor"
    ),
    "triang": VariableMetadata(
        latex=r"$\delta_\mathrm{sep}$", description="Triangularity"
    ),
    "f_tf_steel": VariableMetadata(
        latex=r"f_\mathrm{steel}^\mathrm{TF}", description="TF steel fraction"
    ),
    "plascur/1d6": VariableMetadata(
        latex=r"$I_{\mathrm{p}}$[$MA$]", description="Plasma current"
    ),
    "alstroh": VariableMetadata(
        latex=r"$\sigma_{\mathrm{oh}}^{\mathrm{max}}$[$Pa$]",
        description="Maximum allowable stress",
    ),
    "ohcth": VariableMetadata(
        latex=r"$\Delta R_{\mathrm{CS}}$[$m$]", description="CS coil thickness"
    ),
    "bore": VariableMetadata(
        latex=r"$\Delta R_{\mathrm{bore}}$[$m$]", description="Bore radius"
    ),
    "dnla": VariableMetadata(
        latex=r"$\bar{n}_{\mathrm{e}}$[$m^{-3}$]",
        description="Average electron density",
    ),
    "dnla_gw": VariableMetadata(
        latex=r"$f_{\mathrm{GW}}$", description="Greenwald fraction"
    ),
    "normalised_toroidal_beta": VariableMetadata(
        latex=r"$\beta_{N,\mathrm{tor}}$", description="Normalized toroidal beta"
    ),
    "copperaoh_m2": VariableMetadata(
        latex=r"$\frac{I_{\mathrm{CS}}}{CuA} [$A m$^{-2}$$]$",
        description="Copper-to-area ratio",
    ),
    "copperaoh_m2_max": VariableMetadata(
        latex=r"$max\frac{I_{\mathrm{CS}}}{CuA} [$A m$^{-2}$$]$",
        description="Maximum copper-to-area ratio",
    ),
    "coreradius": VariableMetadata(latex=r"$r_{core} [M]$", description="Core radius"),
    "fcuohsu": VariableMetadata(
        latex=r"$f_{\mathrm{Cu}}^{\mathrm{CS}}$",
        description="Copper fraction of strand in central solenoid",
    ),
    "coheof": VariableMetadata(latex=r"$J [A M^{-2}]$", description="Current density"),
    "ohhghf": VariableMetadata(
        latex=r"$Thickness_{\mathrm{CS}}[m]$",
        description="High field CS coil thickness",
    ),
    "pheat": VariableMetadata(
        latex=r"$ P_{\mathrm{heat}}$ [$MW$]", description="Heat power"
    ),
    "effcd": VariableMetadata(
        latex=r"$\eta_{\mathrm{CD}}$[$A/W$]", description="CD efficiency"
    ),
    "bigq": VariableMetadata(latex=r"$Q$", description="Quality factor"),
    "faccd": VariableMetadata(latex=r"$f_{\mathrm{CD}}$", description="CD factor"),
    "facoh": VariableMetadata(
        latex=r"$f_{\mathrm{CD,ind}}$", description="Inductive CD factor"
    ),
    "bootipf": VariableMetadata(
        latex=r"$f_{\mathrm{BS}}$", description="Bootstrap current fraction"
    ),
    "pdivt": VariableMetadata(
        latex=r"$P_{\mathrm{sep}}$ [$MW$]", description="Separation power"
    ),
    "pradmw": VariableMetadata(
        latex=r"$P_{\mathrm{rad}}$ [$MW$]", description="Radiation power"
    ),
    "pdivtbt/qar": VariableMetadata(
        latex=r"$\frac{P_{\mathrm{sep}}B_T}{q_{95}AR_{\mathrm{maj}}}$ [$MWTm^{-1}$]",
        description="Normalized separation power",
    ),
    "iooic": VariableMetadata(
        latex=r"$I_{\mathrm{TF}}/I_{\mathrm{TF},\mathrm{crit}}$",
        description="Normalized TF current",
    ),
    "bktlife": VariableMetadata(
        latex=r"$T_{\mathrm{blk}}$", description="Blanket lifetime"
    ),
    "bktcycles": VariableMetadata(
        latex=r"$N_{\mathrm{blk},\mathrm{cycle}}$", description="Blanket cycles"
    ),
    "zeff": VariableMetadata(
        latex=r"$Z_{\mathrm{eff}}$", description="Effective charge"
    ),
    "tburn": VariableMetadata(
        latex=r"$t_{\mathrm{burn}}$[$s$]", description="Burn time"
    ),
    "vburn": VariableMetadata(
        latex=r"$V_{\mathrm{loop}}$ [$V$]", description="Loop voltage"
    ),
    "sig_tf_wp_max": VariableMetadata(
        latex=r"$\sigma_{TP,wp}^{max}$", description="Maximum TF winding pack stress"
    ),
    "rli": VariableMetadata(
        latex=r"$l_i$", description="Normalized internal inductance"
    ),
    "n_cycle_min": VariableMetadata(
        latex=r"$MinCycles_{\mathrm{Stress.min}}^{\mathrm{CS}}$",
        description="Minimum cycles for stress",
    ),
    "n_cycle": VariableMetadata(
        latex=r"$Cycles_{\mathrm{Stress}}^{\mathrm{CS}}$",
        description="Cycles for stress",
    ),
    "a_oh_turn": VariableMetadata(
        latex=r"$Turn_{\mathrm{area}}^{\mathrm{CS}}[$m$^{2}]$",
        description="Cross-sectional area of CS coil turns",
    ),
    "tbrnmn": VariableMetadata(
        latex=r"$t_{\mathrm{burn.min}}$[$s$]", description="Minimum burn time"
    ),
    "pfv.oh_steel_frac": VariableMetadata(
        latex=r"$f_{\mathrm{Steel}}^{\mathrm{CS}}$",
        description="Steel fraction in CS coil",
    ),
    "csfv.t_structural_radial": VariableMetadata(
        latex=r"$Turn_{\mathrm{radial}}^{\mathrm{CS}}[$m$]$",
        description="Radial turn length",
    ),
    "csfv.t_crack_vertical": VariableMetadata(
        latex=r"$Crack_{\mathrm{vertical}}^{\mathrm{CS}}[$m$]$",
        description="Vertical crack length",
    ),
    "inlet_temp_liq": VariableMetadata(
        latex=r"Breeder/coolant inlet T [K]",
        description="Breeder/coolant inlet temperature",
    ),
    "outlet_temp_liq": VariableMetadata(
        latex=r"Breeder/coolant outlet T [K]",
        description="Breeder/coolant outlet temperature",
    ),
    "blpressure_liq": VariableMetadata(
        latex=r"Breeder/coolant pressure [Pa]",
        description="Breeder/coolant pressure",
    ),
    "n_liq_recirc": VariableMetadata(
        latex=r"Breeder/coolant recirculations",
        description="Number of breeder/coolant recirculations",
    ),
    "bz_channel_conduct_liq": VariableMetadata(
        latex=r"channel wall conductivity [A V-1 m-1]",
        description="Liquid breeder/coolant channel wall conductivity",
    ),
    "pnuc_fw_ratio_dcll": VariableMetadata(
        latex=r"FW nuclear power fraction",
        description="Nuclear power fraction in the first wall",
    ),
    "f_nuc_pow_bz_struct": VariableMetadata(
        latex=r"BZ structure fraction",
        description="Fraction of nuclear power deposited in the blanket structure",
    ),
    "pitch": VariableMetadata(
        latex=r"FW pitch [m]",
        description="Pitch of the first wall",
    ),
    "coe": VariableMetadata(
        latex=r"$\mathrm{LCOE}$ [$m\$/kWh$]",
        description="Levelized Cost of Electricity",
    ),
    "beta": VariableMetadata(
        latex=r"$\beta$",
        description="Beta (plasma confinement parameter)",
    ),
    "fimp(13)": VariableMetadata(
        latex=r"$Xe_{\mathrm{f}}$",
        description="Impurity fraction (Xenon)",
    ),
    "pdivmax/rmajor": VariableMetadata(
        latex=r"$P_{\mathrm{div}}/R_\mathrm{maj}$ [MW/m]",
        description="Divertor power per major radius",
    ),
    "etath": VariableMetadata(
        latex=r"Thermal to Electric efficiency",
        description="",
    ),
    "fkind": VariableMetadata(
        latex=r"N$^\mathrm{th}$ of a kind factor",
        description="Factor representing the number of distinct scenarios",
    ),
    "startupratio": VariableMetadata(
        latex=r"Gyrotron Redundancy",
        description="Redundancy factor for gyrotrons",
    ),
    "etaech": VariableMetadata(
        latex=r"ECH wall plug to injector efficiency",
        description="Efficiency of electron cyclotron heating",
    ),
}
