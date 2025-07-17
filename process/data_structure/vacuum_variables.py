vacuum_model: str = None
"""switch for vacuum pumping model:

 - ='old' for old detailed ETR model
 - ='simple' for simple steady-state model with comparison to ITER cryopumps
  !#TODO: old and simple not suitable names."""

niterpump: float = None
"""number of high vacuum pumps (real number), each with the throughput of one
ITER cryopump (50 Pa m3 s-1), all operating at the same time (`vacuum_model='simple'`)"""

ntype: int = None
"""switch for vacuum pump type:

   - =0 - for turbomolecular pump (magnetic bearing) with speed of 2.0 m3/s
     (1.95 for N2, 1.8 for He, 1.8 for DT)
   - =1 - for compound cryopump with nominal speed of 10.0 m3/s
     (9.0 for N2, 5.0 for He and 25.0 for DT)"""

nvduct: int = None
"""number of ducts (torus to pumps)"""

dlscal: float = None
"""vacuum system duct length scaling"""

pbase: float = None
"""base pressure during dwell before gas pre-fill(Pa)"""

prdiv: float = None
"""divertor chamber pressure during burn (Pa)"""

pumptp: float = None
"""Pump throughput (molecules/s) (default is ITER value)"""

rat: float = None
"""plasma chamber wall outgassing rate (Pa-m/s)"""

tn: float = None
"""neutral gas temperature in chamber (K)"""

vacdshm: float = None
"""mass of vacuum duct shield (kg)"""

vcdimax: float = None
"""diameter of duct passage (m)"""

vpumpn: int = None
"""number of high vacuum pumps"""

dwell_pump: int = None
"""switch for dwell pumping options:

   - =0 pumping only during t_between_pulse
   - =1 pumping only during t_precharge
   - =2 pumping during t_between_pulse + t_precharge

   The following are used in the Battes, Day and Rohde pump-down model
   See "Basic considerations on the pump-down time in the dwell phase of a pulsed fusion DEMO"
   http://dx.doi.org/10.1016/j.fusengdes.2015.07.011)(vacuum_model=simple')"""

pumpareafraction: float = None
"""area of one pumping port as a fraction of plasma surface area"""

pumpspeedmax: float = None
"""maximum pumping speed per unit area for deuterium & tritium, molecular flow"""

pumpspeedfactor: float = None
"""effective pumping speed reduction factor due to duct impedance"""

initialpressure: float = None
"""initial neutral pressure at the beginning of the dwell phase (Pa)"""

outgasindex: float = None
"""outgassing decay index"""

outgasfactor: float = None
"""outgassing prefactor kw: outgassing rate at 1 s per unit area (Pa m s-1)"""


def init_vacuum_variables():
    """Initialise Vacuum variables"""
    global vacuum_model
    global niterpump
    global ntype
    global nvduct
    global dlscal
    global pbase
    global prdiv
    global pumptp
    global rat
    global tn
    global vacdshm
    global vcdimax
    global vpumpn
    global dwell_pump
    global pumpareafraction
    global pumpspeedmax
    global pumpspeedfactor
    global initialpressure
    global outgasindex
    global outgasfactor

    vacuum_model = "old"
    niterpump = 0.0
    ntype = 1
    nvduct = 0
    dlscal = 0.0
    pbase = 5.0e-4
    prdiv = 0.36
    pumptp = 1.2155e22
    rat = 1.3e-8
    tn = 300.0
    vacdshm = 0.0
    vcdimax = 0.0
    vpumpn = 0
    dwell_pump = 0
    pumpareafraction = 0.0203
    pumpspeedmax = 27.3
    pumpspeedfactor = 0.167
    initialpressure = 1.0
    outgasindex = 1.0
    outgasfactor = 0.0235
