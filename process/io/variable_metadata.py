from dataclasses import dataclass


@dataclass
class VariableMetadata:
    latex: str
    description: str
    units: str


var_dicts = {
    "shldith": VariableMetadata(
        latex=r"$\Delta R_\mathrm{sh}$ [$m$]",
        description="Inboard shield thickness",
        units="m",
    ),
    "rmajor": VariableMetadata(
        latex=r"$R_\mathrm{major}$ [$m$]", description="Major radius", units="m"
    ),
    "crypmw": VariableMetadata(
        latex=r"$P_\mathrm{cryo}$ [$MW$]", description="Cryogenic power", units="MW"
    ),
    "bt": VariableMetadata(
        latex=r"$B_\mathrm{T}$ [$T$]", description="Toroidal magnetic field", units="T"
    ),
    "tfcth": VariableMetadata(
        latex=r"$\Delta R_\mathrm{TF}$ [$m$]",
        description="TF coil thickness",
        units="m",
    ),
    "powfmw": VariableMetadata(
        latex=r"$P_\mathrm{fus}$ [$MW$]", description="Fusion power", units="MW"
    ),
    "pinjemw": VariableMetadata(
        latex=r"$P_\mathrm{inj}$ [$MW$]", description="Injected power", units="MW"
    ),
    "pnetelmw": VariableMetadata(
        latex=r"$P_\mathrm{Net\ elec}$ [$MW$]",
        description="Net electrical power",
        units="MW",
    ),
    "taueff": VariableMetadata(
        latex=r"$\tau_\mathrm{E}$ [s]",
        description="Effective energy confinement time",
        units="s",
    ),
    "ralpne": VariableMetadata(
        latex=r"$f_\mathrm{\alpha}$",
        description="Alpha particle heating fraction",
        units="",
    ),
    "te": VariableMetadata(
        latex=r"$\left< T_\mathrm{e} \right>$",
        description="Average electron temperature",
        units="keV",
    ),
    "taulimit": VariableMetadata(
        latex=r"$max : \frac{\tau_\mathrm{\alpha}}{\tau_\mathrm{E}}$",
        description="Ratio of alpha heating time to energy confinement time",
        units="",
    ),
    "scrapli": VariableMetadata(
        latex=r"$\Delta R_\mathrm{FW-sep}$ [$m$]",
        description="Inner FW to plasma gap",
        units="m",
    ),
    "scraplo": VariableMetadata(
        latex=r"$\Delta R_\mathrm{FW-sep}^\mathrm{out}$ [$m$]",
        description="Outer FW to plasma gap",
        units="m",
    ),
    "vforce": VariableMetadata(
        latex=r"$F_\mathrm{z}^\mathrm{in}$ [$N$]",
        description="TF coil vertical force",
        units="N",
    ),
    "thkcas": VariableMetadata(
        latex=r"$\Delta R_\mathrm{TF}^\mathrm{buck}$ [$m$]",
        description="Inboard TF coil case outer (non-plasma side) thickness",
        units="m",
    ),
    "bmaxtf": VariableMetadata(
        latex=r"$B_\mathrm{TF}^\mathrm{max}$ [$T$]",
        description="Mean peak field at TF coil",
        units="T",
    ),
    "ritfc": VariableMetadata(
        latex=r"$I_\mathrm{TF}^\mathrm{tot}$ [$A$]",
        description="Total TF coil current",
        units="A",
    ),
    "dr_tf_wp": VariableMetadata(
        latex=r"$\Delta R_\mathrm{TF}^\mathrm{WP}$ [$m$]",
        description="TF coil winding pack width",
        units="m",
    ),
    "aspect": VariableMetadata(latex=r"$A$", description="Aspect ratio", units=""),
    "rminor": VariableMetadata(
        latex=r"$a_\mathrm{min}$ [$m$]", description="Minor radius", units="m"
    ),
    "capcost": VariableMetadata(
        latex=r"$C_\mathrm{cap}$ [$M\$ $]", description="Capital cost", units="M$"
    ),
    "r_tf_outboard_mid": VariableMetadata(
        latex=r"$\Delta R_\mathrm{TF}^\mathrm{out\ mid}$ [$m$]",
        description="Mid-plane outboard TF coil leg radius at the middle of the coil",
        units="m",
    ),
    "pgrossmw": VariableMetadata(
        latex=r"$P_\mathrm{gross}^\mathrm{elec}$ [$MW$]",
        description="Gross electrical power",
        units="MW",
    ),
    "htpmw": VariableMetadata(
        latex=r"$P_\mathrm{Primary\ coolant}^\mathrm{elec}$ [$MW$]",
        description="Primary coolant power",
        units="MW",
    ),
    "ppfmw": VariableMetadata(
        latex=r"$P_\mathrm{PF}^\mathrm{elec}$ [$MW$]",
        description="Poloidal field power",
        units="MW",
    ),
    "hmax": VariableMetadata(
        latex=r"$z_\mathrm{TF}^\mathrm{pl\ side}$ [$m$]",
        description="maximum (half-)height of TF coil (inside edge)",
        units="m",
    ),
    "thicndut": VariableMetadata(
        latex=r"\Delta l_\mathrm{steel\ jacket}^\mathrm{turn}",
        description="Thickness of steel jacket turn",
        units="",
    ),
    "cpttf": VariableMetadata(
        latex=r"$I_\mathrm{TF}^\mathrm{turn}$ [$A$]",
        description="TF turn current",
        units="A",
    ),
    "boundl(2)": VariableMetadata(
        latex=r"$B_\mathrm{T}^\mathrm{min}$ [$A$]",
        description="Toroidal field lower bound",
        units="A",
    ),
    "pinjmw": VariableMetadata(
        latex=r"$P_\mathrm{inj}$ [$MW$]", description="Injected power", units="MW"
    ),
    "hldivlim": VariableMetadata(
        latex=r"$q_\mathrm{div}^\mathrm{max}$ [$MW.m^{-2}$]",
        description="Maximum divertor heat load",
        units="MW.m^{-2}",
    ),
    "hfact": VariableMetadata(
        latex=r"$f_\mathrm{H}$", description="H-factor", units=""
    ),
    "kappa": VariableMetadata(
        latex=r"$\kappa_\mathrm{sep}$", description="Elongation", units=""
    ),
    "triang": VariableMetadata(
        latex=r"$\delta_\mathrm{sep}$", description="Triangularity", units=""
    ),
    "f_tf_steel": VariableMetadata(
        latex=r"f_\mathrm{steel}^\mathrm{TF}", description="TF steel fraction", units=""
    ),
    "plascur/1d6": VariableMetadata(
        latex=r"$I_{\mathrm{p}}$[$MA$]", description="Plasma current", units="MA"
    ),
    "n_cycle": VariableMetadata(
        latex=r"$N_{\mathrm{CS},\mathrm{cycle}}$",
        description="Number of cycles",
        units="",
    ),
    "alstroh": VariableMetadata(
        latex=r"$\sigma_{\mathrm{oh}}^{\mathrm{max}}$[$Pa$]",
        description="Maximum allowable stress",
        units="Pa",
    ),
    "ohcth": VariableMetadata(
        latex=r"$\Delta R_{\mathrm{CS}}$[$m$]",
        description="CS coil thickness",
        units="m",
    ),
    "bore": VariableMetadata(
        latex=r"$\Delta R_{\mathrm{bore}}$[$m$]", description="Bore radius", units="m"
    ),
    "dnla": VariableMetadata(
        latex=r"$\bar{n}_{\mathrm{e}}$[$m^{-3}$]",
        description="Average electron density",
        units="m^{-3}",
    ),
    "dnla_gw": VariableMetadata(
        latex=r"$f_{\mathrm{GW}}$", description="Greenwald fraction", units=""
    ),
    "normalised_toroidal_beta": VariableMetadata(
        latex=r"$\beta_{N,\mathrm{tor}}$",
        description="Normalized toroidal beta",
        units="",
    ),
    "copperaoh_m2": VariableMetadata(
        latex=r"$\frac{I_{\mathrm{CS}}}{CuA} [$A m$^{-2}$$]$",
        description="Copper-to-area ratio",
        units="A m^{-2}",
    ),
    "copperaoh_m2_max": VariableMetadata(
        latex=r"$max\frac{I_{\mathrm{CS}}}{CuA} [$A m$^{-2}$$]$",
        description="Maximum copper-to-area ratio",
        units="A m^{-2}",
    ),
    "coreradius": VariableMetadata(
        latex=r"$r_{core} [M]$", description="Core radius", units="M"
    ),
    "fcuohsu": VariableMetadata(
        latex=r"$f_{\mathrm{Cu}}^{\mathrm{CS}}$",
        description="Copper fraction of strand in central solenoid",
        units="",
    ),
    "coheof": VariableMetadata(
        latex=r"$J [A m^{-2}]$",
        description="central solenoid overall current density at end of flat-top",
        units="A m^{-2}",
    ),
    "ohhghf": VariableMetadata(
        latex=r"$Thickness_{\mathrm{CS}}[m]$",
        description="Central solenoid height / TF coil internal height",
        units="m",
    ),
    "pheat": VariableMetadata(
        latex=r"$ P_{\mathrm{heat}}$ [$MW$]", description="Heat power", units="MW"
    ),
    "effcd": VariableMetadata(
        latex=r"$\eta_{\mathrm{CD}}$[$A/W$]", description="CD efficiency", units="A/W"
    ),
    "bigq": VariableMetadata(latex=r"$Q$", description="Plasma Q value", units=""),
    "faccd": VariableMetadata(
        latex=r"$f_{\mathrm{CD}}$", description="CD factor", units=""
    ),
    "facoh": VariableMetadata(
        latex=r"$f_{\mathrm{CD,ind}}$", description="Inductive CD factor", units=""
    ),
    "bootipf": VariableMetadata(
        latex=r"$f_{\mathrm{BS}}$", description="Bootstrap current fraction", units=""
    ),
    "pdivt": VariableMetadata(
        latex=r"$P_{\mathrm{sep}}$ [$MW$]", description="Power to divertor", units="MW"
    ),
    "pradmw": VariableMetadata(
        latex=r"$P_{\mathrm{rad}}$ [$MW$]", description="Radiation power", units="MW"
    ),
    "pdivtbt/qar": VariableMetadata(
        latex=r"$\frac{P_{\mathrm{sep}}B_T}{q_{95}AR_{\mathrm{maj}}}$ [$MWTm^{-1}$]",
        description="",
        units="MWTm^{-1}",
    ),
    "iooic": VariableMetadata(
        latex=r"$I_{\mathrm{TF}}/I_{\mathrm{TF},\mathrm{crit}}$",
        description="Normalized TF current",
        units="",
    ),
    "bktlife": VariableMetadata(
        latex=r"$T_{\mathrm{blk}}$", description="Blanket lifetime", units=""
    ),
    "bktcycles": VariableMetadata(
        latex=r"$N_{\mathrm{blk},\mathrm{cycle}}$",
        description="Blanket cycles",
        units="",
    ),
    "zeff": VariableMetadata(
        latex=r"$Z_{\mathrm{eff}}$", description="Effective charge", units=""
    ),
    "tburn": VariableMetadata(
        latex=r"$t_{\mathrm{burn}}$[$s$]", description="Burn time", units="s"
    ),
    "vburn": VariableMetadata(
        latex=r"$V_{\mathrm{loop}}$ [$V$]",
        description="Loop voltage during burn",
        units="V",
    ),
    "sig_tf_wp_max": VariableMetadata(
        latex=r"$\sigma_{TP,wp}^{max}$",
        description="Maximum TF winding pack stress",
        units="",
    ),
    "rli": VariableMetadata(
        latex=r"$l_i$", description="Normalized internal inductance", units=""
    ),
    "n_cycle_min": VariableMetadata(
        latex=r"$MinCycles_{\mathrm{Stress.min}}^{\mathrm{CS}}$",
        description="Minimum cycles for stress",
        units="",
    ),
    "a_oh_turn": VariableMetadata(
        latex=r"$Turn_{\mathrm{area}}^{\mathrm{CS}}[$m$^{2}]$",
        description="Cross-sectional area of CS coil turns",
        units="m^2",
    ),
    "tbrnmn": VariableMetadata(
        latex=r"$t_{\mathrm{burn.min}}$[$s$]",
        description="Minimum burn time",
        units="s",
    ),
    "pfv.oh_steel_frac": VariableMetadata(
        latex=r"$f_{\mathrm{Steel}}^{\mathrm{CS}}$",
        description="Steel fraction in CS coil",
        units="",
    ),
    "csfv.t_structural_radial": VariableMetadata(
        latex=r"$Turn_{\mathrm{radial}}^{\mathrm{CS}}[$m$]$",
        description="Radial turn length",
        units="m",
    ),
    "csfv.t_crack_vertical": VariableMetadata(
        latex=r"$Crack_{\mathrm{vertical}}^{\mathrm{CS}}[$m$]$",
        description="Vertical crack length",
        units="m",
    ),
    "inlet_temp_liq": VariableMetadata(
        latex=r"Breeder/coolant inlet T [K]",
        description="Breeder/coolant inlet temperature",
        units="K",
    ),
    "outlet_temp_liq": VariableMetadata(
        latex=r"Breeder/coolant outlet T [K]",
        description="Breeder/coolant outlet temperature",
        units="K",
    ),
    "blpressure_liq": VariableMetadata(
        latex=r"Breeder/coolant pressure [Pa]",
        description="Breeder/coolant pressure",
        units="Pa",
    ),
    "n_liq_recirc": VariableMetadata(
        latex=r"Breeder/coolant recirculations",
        description="Number of breeder/coolant recirculations",
        units="",
    ),
    "bz_channel_conduct_liq": VariableMetadata(
        latex=r"channel wall conductivity [A V-1 m-1]",
        description="Liquid breeder/coolant channel wall conductivity",
        units="A V-1 m-1",
    ),
    "pnuc_fw_ratio_dcll": VariableMetadata(
        latex=r"FW nuclear power fraction",
        description="Nuclear power fraction in the first wall",
        units="",
    ),
    "f_nuc_pow_bz_struct": VariableMetadata(
        latex=r"BZ structure fraction",
        description="Fraction of nuclear power deposited in the blanket structure",
        units="",
    ),
    "pitch": VariableMetadata(
        latex=r"FW pitch [m]", description="Pitch of the first wall", units="m"
    ),
    "coe": VariableMetadata(
        latex=r"$\mathrm{LCOE}$ [$m\$/kWh$]",
        description="Levelized Cost of Electricity",
        units="m$/kWh",
    ),
    "beta": VariableMetadata(
        latex=r"$\beta$", description="Total plasma beta", units=""
    ),
    "fimp(13)": VariableMetadata(
        latex=r"$Xe_{\mathrm{f}}$", description="Impurity fraction (Xenon)", units=""
    ),
    "pdivmax/rmajor": VariableMetadata(
        latex=r"$P_{\mathrm{div}}/R_\mathrm{maj}$ [MW/m]",
        description="Divertor power per major radius",
        units="MW/m",
    ),
    "etath": VariableMetadata(
        latex=r"Thermal to Electric efficiency",
        description="Thermal to electric efficiency",
        units="",
    ),
    "fkind": VariableMetadata(
        latex=r"N$^\mathrm{th}$ of a kind factor",
        description="Multiplier for Nth of a kind costs",
        units="",
    ),
    "startupratio": VariableMetadata(
        latex=r"Gyrotron Redundancy",
        description="Redundancy factor for gyrotrons",
        units="",
    ),
    "etaech": VariableMetadata(
        latex=r"ECH wall plug to injector efficiency",
        description="Efficiency of electron cyclotron heating",
        units="",
    ),
    "tauee": VariableMetadata(
        latex=r"$\tau_E$",
        description="Electron energy confinement time (sec)",
        units="s",
    ),
    "dene": VariableMetadata(
        latex=r"$n_e$",
        description="Volume-averaged electron density (/m3)",
        units="m-3",
    ),
}
