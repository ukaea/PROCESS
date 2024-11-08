module constraint_variables
  !! author: J. Morris (UKAEA)
  !!
  !! This module contains global variables relating to the constraint
  !! equations (f-values, limits, etc.).
  !!
  !!### References
  !!
  !! -
#ifndef dp
  use, intrinsic :: iso_fortran_env, only: dp=>real64
#endif

  implicit none

  public

  real(dp) :: auxmin
  !! minimum auxiliary power (MW) (`constraint equation 40`)

  real(dp) :: beta_poloidal_max
  !! maximum poloidal beta (`constraint equation 48`)

  real(dp) :: bigqmin
  !! minimum fusion gain Q (`constraint equation 28`)

  real(dp) :: bmxlim
  !! maximum peak toroidal field (T) (`constraint equation 25`)

  real(dp) :: fauxmn
  !! f-value for minimum auxiliary power (`constraint equation 40`, `iteration variable 64`)

  real(dp) :: fbeta
  !! f-value for epsilon beta-poloidal (`constraint equation 6`, `iteration variable 8`)

  real(dp) :: fbetap
  !! f-value for poloidal beta (`constraint equation 48`, `iteration variable 79`)

  real(dp) :: fbetatry
  !! f-value for beta limit (`constraint equation 24`, `iteration variable 36`)

  real(dp) :: fbetatry_lower
  !! f-value for (lower) beta limit (`constraint equation 84`, `iteration variable 173`)

  real(dp) :: fcpttf
  !! f-value for TF coil current per turn upper limit
  !! (`constraint equation 77`, `iteration variable 146`)

  real(dp) :: fcwr
  !! f-value for conducting wall radius / rminor limit
  !! (`constraint equation 23`, `iteration variable 104`)

  real(dp) :: fdene
  !! f-value for density limit (`constraint equation 5`, `iteration variable 9`)
  !! (invalid if `ipedestal=3`)

  real(dp) :: fdivcol
  !! f-value for divertor collisionality (`constraint equation 22`, `iteration variable 34`)

  real(dp) :: fdtmp
  !! f-value for first wall coolant temperature rise
  !! (`constraint equation 38`, `iteration variable 62`)

  real(dp) :: fecrh_ignition
  !! f-value for ecrh ignition constraint
  !! (`constraint equation 91`, `iteration variable 168`)

  real(dp) :: fflutf
  !! f-value for neutron fluence on TF coil (`constraint equation 53`, `iteration variable 92`)

  real(dp) :: ffuspow
  !! f-value for maximum fusion power (`constraint equation 9`, `iteration variable 26`)

  real(dp) :: fgamcd
  !! f-value for current drive gamma (`constraint equation 37`, `iteration variable 40`)

  real(dp) :: fhldiv
  !! f-value for divertor heat load (`constraint equation 18`, `iteration variable 27`)

  real(dp) :: fiooic
  !! f-value for TF coil operating current / critical current ratio
  !! (`constraint equation 33`, `iteration variable 50`)

  real(dp) :: fipir
  !! f-value for Ip/Irod upper limit
  !! constraint equation icc = 46
  !! iteration variable ixc = 72

  real(dp) :: fjohc
  !! f-value for central solenoid current at end-of-flattop
  !! (`constraint equation 26`, `iteration variable 38`)

  real(dp) :: fjohc0
  !! f-value for central solenoid current at beginning of pulse
  !! (`constraint equation 27`, `iteration variable 39`)

  real(dp) :: fjprot
  !! f-value for TF coil winding pack current density
  !! (`constraint equation 35`, `iteration variable 53`)

  real(dp) :: flhthresh
  !! f-value for L-H power threshold (`constraint equation 15`, `iteration variable 103`)

  real(dp) :: fmva
  !! f-value for maximum MVA (`constraint equation 19`, `iteration variable 30`)

  real(dp) :: fnbshinef
  !! f-value for maximum neutral beam shine-through fraction
  !! (`constraint equation 59`, `iteration variable 105`)

  real(dp) :: fncycle
  !! f-value for minimum CS coil stress load cycles
  !! (`constraint equation 90`, `iteration variable 167`)

  real(dp) :: fnesep
  !! f-value for Eich critical separatrix density
  !! (`constraint equation 76`, `iteration variable 144`)

  real(dp) :: foh_stress
  !! f-value for Tresca yield criterion in Central Solenoid
  !! (`constraint equation 72`, `iteration variable 123`)

  real(dp) :: fpeakb
  !! f-value for maximum toroidal field (`constraint equation 25`, `iteration variable 35`)

  real(dp) :: fpinj
  !! f-value for injection power (`constraint equation 30`, `iteration variable 46`)

  real(dp) :: fpnetel
  !! f-value for net electric power (`constraint equation 16`, `iteration variable 25`)

  real(dp) :: fportsz
  !! f-value for neutral beam tangency radius limit
  !! (`constraint equation 20`, `iteration variable 33`)

  real(dp) :: fpsepbqar
  !! f-value for maximum Psep*Bt/qAR limit (`constraint equation 68`, `iteration variable 117`)

  real(dp) :: fpsepr
  !! f-value for maximum Psep/R limit (`constraint equation 56`, `iteration variable 97`)

  real(dp) :: fptemp
  !! f-value for peak centrepost temperature (`constraint equation 44`, `iteration variable 68`)

  real(dp) :: fptfnuc
  !! f-value for maximum TF coil nuclear heating (`constraint equation 54`, `iteration variable 95`)

  real(dp) :: fq
  !! f-value for edge safety factor (`constraint equation 45`, `iteration variable 71`)

  real(dp) :: fqval
  !! f-value for Q (`constraint equation 28`, `iteration variable 45`)

  real(dp) :: fradpwr
  !! f-value for core radiation power limit (`constraint equation 17`, `iteration variable 28`)

  real(dp) :: fradwall
  !! f-value for upper limit on radiation wall load (`constr. equ. 67`, `iteration variable 116`)

  real(dp) :: freinke
  !! f-value for Reinke detachment criterion (`constr. equ. 78`, `iteration variable 147`)

  real(dp) :: frminor
  !! f-value for minor radius limit (`constraint equation 21`, `iteration variable 32`)

  real(dp) :: fstrcase
  !! f-value for maximum TF coil case Tresca yield criterion
  !! (`constraint equation 31`, `iteration variable 48`)

  real(dp) :: fstrcond
  !! f-value for maxiumum TF coil conduit Tresca yield criterion
  !! (`constraint equation 32`, `iteration variable 49`)

  real(dp) :: fstr_wp
  !! f-value for maxiumum TF coil strain absolute value
  !! (`constraint equation 88`, `iteration variable 165`)

  real(dp) :: fmaxvvstress
  !! f-value for maximum permitted stress of the VV
  !! (`constraint equation 65`, `iteration variable 113`)

  real(dp) :: ftbr
  !! f-value for minimum tritium breeding ratio (`constraint equation 52`, `iteration variable 89`)

  real(dp) :: ftburn
  !! f-value for minimum burn time (`constraint equation 13`, `iteration variable 21`)

  real(dp) :: ftcycl
  !! f-value for cycle time (`constraint equation 42`, `iteration variable 67`)

  real(dp) :: ftmargoh
  !! f-value for central solenoid temperature margin
  !! (`constraint equation 60`, `iteration variable 106`)

  real(dp) :: ftmargtf
  !! f-value for TF coil temperature margin (`constraint equation 36`, `iteration variable 54`)

  real(dp) :: ftohs
  !! f-value for plasma current ramp-up time (`constraint equation 41`, `iteration variable 66`)

  real(dp) :: ftpeak
  !! f-value for first wall peak temperature (`constraint equation 39`, `iteration variable 63`)

  real(dp) :: fvdump
  !! f-value for dump voltage (`constraint equation 34`, `iteration variable 51`)

  real(dp) :: fvs
  !! f-value for flux-swing (V-s) requirement (STEADY STATE)
  !! (`constraint equation 12`, `iteration variable 15`)

  real(dp) :: fvvhe
  !! f-value for vacuum vessel He concentration limit (`iblanket = 2`)
  !! (`constraint equation 55`, `iteration variable 96`)

  real(dp) :: fwalld
  !! f-value for maximum wall load (`constraint equation 8`, `iteration variable 14`)

  real(dp) :: fzeffmax
  !! f-value for maximum zeff (`constraint equation 64`, `iteration variable 112`)

  real(dp) :: gammax
  !! maximum current drive gamma (`constraint equation 37`)

  real(dp) :: maxradwallload
  !!  Maximum permitted radiation wall load (MW/m^2) (`constraint equation 67`)

  real(dp) :: mvalim
  !! maximum MVA limit (`constraint equation 19`)

  real(dp) :: nbshinefmax
  !! maximum neutral beam shine-through fraction (`constraint equation 59`)

  real(dp) :: nflutfmax
  !! max fast neutron fluence on TF coil (n/m2) (`blktmodel>0`) (`constraint equation 53`)
  !! Also used for demontable magnets (itart = 1) and superconducting coils (i_tf_sup = 1)
  !! To set the CP lifetime (`constraint equation 85`)

  real(dp) :: pdivtlim
  !! Minimum pdivt [MW] (`constraint equation 80`)

  real(dp) :: peakfactrad
  !! peaking factor for radiation wall load (`constraint equation 67`)

  real(dp) :: peakradwallload
  !! Peak radiation wall load (MW/m^2) (`constraint equation 67`)

  real(dp) :: pnetelin
  !! required net electric power (MW) (`constraint equation 16`)

  real(dp) :: powfmax
  !! maximum fusion power (MW) (`constraint equation 9`)

  real(dp) :: psepbqarmax
  !! maximum ratio of Psep*Bt/qAR (MWT/m) (`constraint equation 68`)

  real(dp) :: pseprmax
  !! maximum ratio of power crossing the separatrix to plasma major radius (Psep/R) (MW/m)
  !! (`constraint equation 56`)

  real(dp) :: ptfnucmax
  !! maximum nuclear heating in TF coil (MW/m3) (`constraint equation 54`)

  real(dp) :: tbrmin
  !! minimum tritium breeding ratio (`constraint equation 52`)

  real(dp) :: tbrnmn
  !! minimum burn time (s) (KE - no longer itv., see issue #706)

  real(dp) :: tcycmn
  !! minimum cycle time (s) (`constraint equation 42`)

  real(dp) :: tohsmn
  !! minimum plasma current ramp-up time (s) (`constraint equation 41`)

  real(dp) :: vvhealw
  !! allowed maximum helium concentration in vacuum vessel at end of plant life (appm)
  !! (`iblanket =2`) (`constraint equation 55`)

  real(dp) :: walalw
  !! allowable neutron wall-load (MW/m2) (`constraint equation 8`)

  real(dp) :: taulimit
  !! Lower limit on taup/taueff the ratio of alpha particle to energy confinement
  !! times (`constraint equation 62`)

  real(dp) :: ftaulimit
  !! f-value for lower limit on taup/taueff the ratio of alpha particle to energy
  !! confinement times (`constraint equation 62`, `iteration variable 110`)

  real(dp) :: fniterpump
  !! f-value for constraint that number of pumps < tfno
  !! (`constraint equation 63`, `iteration variable 111`)

  real(dp) :: zeffmax
  !! maximum value for Zeff (`constraint equation 64`)

  real(dp) :: fpoloidalpower
  !! f-value for constraint on rate of change of energy in poloidal field
  !! (`constraint equation 66`, `iteration variable 115`)

  real(dp) :: fpsep
  !! f-value to ensure separatrix power is less than value from Kallenbach divertor
  !! (Not required as constraint 69 is an equality)

  real(dp) :: fcqt
  !! TF coil quench temparature remains below tmax_croco
  !! (`constraint equation 74`, `iteration variable 141`)

  contains

  subroutine init_constraint_variables
    !! Initialise module variables
    implicit none

    auxmin = 0.1D0
    beta_poloidal_max = 0.19D0
    bigqmin = 10.0D0
    bmxlim = 12.0D0
    fauxmn = 1.0D0
    fbeta = 1.0D0
    fbetap = 1.0D0
    fbetatry = 1.0D0
    fbetatry_lower = 1.0D0
    fcpttf = 1.0D0
    fcwr = 1.0D0
    fdene = 1.0D0
    fdivcol = 1.0D0
    fdtmp = 1.0D0
    fflutf = 1.0D0
    ffuspow = 1.0D0
    fgamcd = 1.0D0
    fhldiv = 1.0D0
    fiooic = 0.5D0
    fipir = 1.0D0
    fjohc = 1.0D0
    fjohc0 = 1.0D0
    fjprot = 1.0D0
    flhthresh = 1.0D0
    fmva = 1.0D0
    fnbshinef = 1.0D0
    fncycle = 1.0D0
    fnesep = 1.0D0
    foh_stress = 1.0D0
    fpeakb = 1.0D0
    fpinj = 1.0D0
    fpnetel = 1.0D0
    fportsz = 1.0D0
    fpsepbqar = 1.0D0
    fpsepr = 1.0D0
    fptemp = 1.0D0
    fptfnuc = 1.0D0
    fq = 1.0D0
    fqval = 1.0D0
    fradpwr = 0.99D0
    fradwall = 1.0D0
    freinke = 1.0D0
    frminor = 1.0D0
    fstrcase = 1.0D0
    fstrcond = 1.0D0
    fstr_wp = 1.0D0
    fmaxvvstress = 1.0D0
    ftbr = 1.0D0
    ftburn = 1.0D0
    ftcycl = 1.0D0
    ftmargoh = 1.0D0
    ftmargtf = 1.0D0
    ftohs = 1.0D0
    ftpeak = 1.0D0
    fvdump = 1.0D0
    fvs = 1.0D0
    fvvhe = 1.0D0
    fwalld = 1.0D0
    fzeffmax = 1.0D0
    gammax = 2.0D0
    maxradwallload = 1.0D0
    mvalim = 40.0D0
    nbshinefmax = 1.0D-3
    nflutfmax = 1.0D23
    pdivtlim = 150.0D0
    peakfactrad = 3.33D0
    peakradwallload = 0.0D0
    pnetelin = 1.0D3
    powfmax = 1.5D3
    psepbqarmax = 9.5D0
    pseprmax = 25.0D0
    ptfnucmax = 1.0D-3
    tbrmin = 1.1D0
    tbrnmn = 1.0D0
    tcycmn = 0.0D0
    tohsmn = 1.0D0
    vvhealw = 1.0D0
    walalw = 1.0D0
    taulimit = 5.0D0
    ftaulimit = 1.0D0
    fniterpump = 1.0D0
    zeffmax = 3.6D0
    fpoloidalpower = 1.0D0
    fpsep = 1.0D0
    fcqt = 1.0D0
    fecrh_ignition = 1.0D0
  end subroutine init_constraint_variables
end module constraint_variables
