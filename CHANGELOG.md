# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.4.0] - 2022-05-18
### Added
- Support for debug builds using `-DCMAKE_BUILD_TYPE=Debug` flag (#1277)

### Changed
- Converted `pfcoil.f90` to Python (#1456)
- ChangeLog update procedure (#1574)
- Converted `pulse.f90` from Fortran to Python
- Converted `buildings_module.f90` to pure Python (#1552)
- Converted `water_usage.f90` to Python (#1570)
- Converted `machine_build.f90` to Python (#1576)
- Converted `kallenbach_module.f90` to Python, with the non-scanning part becoming utilities (#1579)
- Converted `structure.f90` to Python (#1538)
- Converted `costs.f90` to Python (#1637)

### Fixed
- Bug where fpdivlim was being called for iteration variable 154 when it should be fne0 (#1659)
- thwcndut may now be input as zero. This is a conduit-free winding pack.
- Divertor's `run` method was not being called with `output=True` from `output.py` (#1553)
- Tracker's plotting inner loop should run over `set(titles)` to avoid running in `O(n)` opposed to `O(1)` (#1520)
- Fixed a typo in the implementation of `extended_plane_strain` from `sctfcoil.f90` (#1565)

### Removed
- Removed Python 3.6 support in CI (#1490)
- Majority of `fispact` code (`lib/fispact/` and `fispact.f90`) as it is no longer used in PROCESS (#1650)


## [2.3.0] - 2022-01-20
### Added
- Create scan notebook template (#1497)
- Create example Jupyter notebook (#1487)
- CS fatigue model (#1400)
- Tracking for cost model (#1344)

### Changed
- Beta-norm and elongation aspect ratio dependence from Menard et al. 2016 (#1439)
- Convert tfcoil into Python/Hybrid Python-Fortran (#1452)
- Update LTS/HTS cost model (#1343)
- Convert costs step to pure python (#1448)
- Update HCD costs in cost model 2 (#1304)
- Elongation scaling aspect ratio stability margin (#1399)
- General PF coil placement (#1418)
- Axisymmetric extended plane strain (#998, #1414)

### Fixed
- 5% regression test job missing values (#1503)
- N_cycle_fix (#1506)
- Plasma-side case of inboard TF coil (`casthi`) is neglected in computing radial stress distribution (#1509)
- Generalized plane strain TF coil stress model is not regression tested with a bucking cylinder (#1442)
- Resolve "REBCO current density incorrectly calculated (#1494)
- Resolve "Incorrect calculation of building sizes (#1395)

## [2.2.0] - 2021-10-26
### Added
- Add pumping variables as inputs for varying the efficiency of thermal to electrical conversion in develop-stable (#1374)
- Added xi_ebw input variable and set default value (#1371)
- Remove f90wrap from build process (#1397)
- Unit tests for costs_step added (#1372)
- Tracking for cost model (#1344)
- Convert caller.f90 to Python (#1354)
- EBW heating system (#1262)
- Call fcnvmc1 and 2 from Python (#1079)

## Changed
- Fix the blnkith variable so that an Inboard blanket can be used whilst itart = 1 (Spherical Tokamak) (#1388)
- Elongation should be scaleable with aspect ratio (#1345)
- PROCESS TF geometry output issue (#1328)
- TF coil stress constraints using different limits for case and WP (#1327)
- Add centrepost nuclear heating to TF nuclear heating (#1331)
- Updated proc_plot.py for EBW plasma heating fix (#1387)
- change the formulae for the dittus-boelter eqn in tfcoil.f90 as it is wrong in develop-stable (#1375)
- Higher current density required in TF coil REBCO tape" - new sc material (#1350)
- Update of Shield Cost Calculations (#1407)
- PROCESS TF geometry output issue (#1328)
- Remove Cost Model 2 link to LSA (#1361)
- Issue 935 power plant water use (#935)
- Secondary EBW Heating Missing (#1324)

## [2.1.2] - 2021-07-01
### Added
- Add costs and temperature margin to Durham REBCO model (#1333)
- Parameterised regression test names added in pytest reporting (#1315)
- BoP cost for Cost model 2 (#1295)
- Add cryo-aluminium costing (#1272)
- Option to use toroidal beta in beta limit (#1290)
- Add array bounds check to Fortran compilation flags (#1294)

### Changed
- cmake configuration step made more robust (#1326)
- Pf coil and BB cost model 2 update
- Enforced up/down symmetry in the vertical build for double-null - (#1309)
- Update indirect costs in cost model 2 (#1296)
- Update plot_scan tool labelling
- dcond dimension fix (#1297)
- Updated EBW gamma
- CMake build step fixed when using numpy v1.20.0 (#1258)

## [2.1.1] - 2021-03-17
### Added
- Adding the EBW scaling (#1262)
- Cost model 2 adds remote handing costs
- Cryopower constraint (#1120)
- Find a way to only review regression test diffs >5% (#1242)
- Adding the insulation effect to the Young modulus smearing in the - vertical direction
- Adding the scan plot utility
- Create tool for updating test references (#1239)

### Changed
- Move to semantic versioning (#1241)
- Move obsolete_vars.py into Python package for Blueprint use (#1279)
- Joints heating error (#1264)
- Corrected some power accounting issues in creating the sankey diagram - (#1265)
- Plot ploc plotting (#1217)
- Break vmcon into smaller subroutines (#1078)
- Python-Fortran dictionaries not updating (#1235)
- Updated core references in documentation (#1247)

## [2.1] - 2021-01-25
### Added
- VMCON unit tests now run in test suite (#1077)
- Added ground insulation layer to resistive magnets (#1122)
- Implemented code coverage reporting for the regression tests (#?)
- Python 3.6 CI test jobs added (#1223)
- Added Flinter FORTRAN Code Quality Scoring (#?)

### Changed
- Inboard and outboard vacuum vessel thickness variables changed (#433)
- New fixed turn cable size formulation (#1182)

### Removed
- Removed Freia support: now untested in CI (#1211)

## [2.0] - 2020-12-17
### Bug Fixes
- Issue 1029 compiler warning on n contact tot (#1029)
- Fixed bug in TorGA interface which breaks under GFortran-10 (#?)
- Handle zero value in IN.DAT writer (#1101)
- Latest markdown version breaks Ford on develop (#1145)
- CMakeLists.txt not working with ; in commit message (#1089)
- HTML introduced into IN.DAT comments (#1124)
- Fixes Failing Tests on Freia (#?)

### Features
- Connecting process output with other codes using json files (#1017)
- Allow Process to be installed as a Python package (#1013)
- Cryo-aluminium magnet model updates (#1016)
- Issue 1021 install script (#1021)
- Issue 999 New coil module for stellarator.f90 (#999)
- Issue 1015 costing (#1015)
- Issue 1010 LH (#1010)
- Issue 1035 JSON DN (#1035)
- Issue 1031 bucked and wedged (#1031)
- Issue 1044 wp geometry (#1044)
- Durham Nb-Ti model based on Ginsburg-Landau theory (#?)
- Issue 1053 resistive magnet updates (#1053)
- Issue 1063 ip warn (#1063)
- Issue 1054 sc cp neutronic shielding (#1054)
- Issue 1086 Elongation and Triangularity (#1086)
- Add Steady-State DEMO documentation and test suite (#?)
- Issue-1114 (#1114)
- Issue 1132 resistive model stability (#1132)
- Issue 1085 durham rebco v1 (#1085)
- Issue 1168 plot proc with rebco (#1168)
- issue-1171 : Adding the shield (#1171)
- Issue 1167 cost model 2 scalings (#1167)
- Adding a line setting the TF coil (#1165)
- Issue 866 ripple for resistive TF coils & 1143 sidewall thickness parametrization (#866)
- Convert project to Python-wrapped Fortran (#1092)

### Minor Changes
- Real declarations in PROCESS (#1027)
- Make subroutines/functions explicitly use module variables (#980)
- Updated base docker image Dockerfile for PROCESS (#1108)
- Use find_package macro to locate GTest, removing need for GTEST_DIR Env variable (#?)
- Issue 1175 (#1175)
- Issue 899 error handling (and other issues detailed on commit messages) (#899)

### Documentation updates
- Tidy docstrings to Ford format (#1026)
- Issue 944 updating gitpages (#944)
- Issue 1002 pl scenario docs (#1002)
- Get variable descriptions from inside the Ford site to display in gitlab pages (#948)
- Remove unwanted and empty modules from vardes (#1040)
- Issue 944 updating gitpages (#944)
- Large commented section in `initial.f90` (#1103)
- Issue 1125 st documentation update (#1125)
- Correcting relative path for included image (#1163)
- Issue 1173 update developper documentation (#1173)
- Issue 1179 update tf ripple doc (#1179)

## 1.0.17

### Bug Fixes

### Features
- Restored IFE (Issue #901)
- GitLab pages now hosts the autodoc html output (only develop) (issue #418)
- Git branch now in output file (#912)
- CI jobs now run for cmake custom targets
- Added diamagnetic and Pfirsch-Schlüter current scalings #992

### Minor Changes
- Renamed 'test_files' to 'unit_tests' #972

### Documentation update


## 1.0.16

### Bug Fixes

- Fixed units issue with Lang et al. (2012) confinement scaling (#821)
- Fixed issue with error numbering (#826)
- Fixed issue with costing of TF coil dump resistors (#847)

### Features

- new command line argument `./process.exe help` provides help info
- new CMake option `-Ddll=ON/OFF`. Default is `ON`. Useful for profiling with gprof
- Updated version of Kallenbach testing (now can run test case that matches Kallenbach
  paper or user defined inputs). See Userguide for more info.
- Updated version of Kallenbach scanning (now can specify the variable to scan
  from a selection, number of scan points etc.). See Userguide for more info.
- Added NSTX and NSTX-Petty08 confinement time scalings (#820)
- Added option to input the confinement time
- Added a new spherical tokamak plasma current relation based on FIESTA fitting
- CI now runs on all branches named "issue-*"
- Unit tests incorporated into main branch

### Minor Changes

- Increased the number of scan points to 1000 (issue #809)
- For issue #379 constraint 52 now gives warning for iblanket=1
- Updated FNSF test case (#822)
- Removed obsolete variable estotf (#199 #847)

### Documentation update

- Developper documentation update (code description/compilation/git instructions)

## 1.0.15

### Bug Fixes
- Added emultmw calculation to stellarator and fixed power balance errors (Issue #783)
- Amended fpump* output to match with primary_pumping options.
- Corrected power crossing the separatrix for stellarators (Issue #787)
- Changed Connor-Hastie plasma current model to kappa95 and triang95% (Issue #791)

### Features
- HTS REBCO model final version implemented
- Can now limit the CS peak field to be below set maximum
- Added Hubbard 2012 and 2017 I-mode threshold scaling
- Added Hubbard I mode confinement time scaling
- Added I-mode version of Reinke criterion (fzmin)
- New figure of linear combination figure of merit. Linear combination (50/50
  weighted) of $`Q`$ and $`t_{burn}`$.
- I mode scalings for confinement time and L-I power threshold from Hubbard 2017.
- New utility called `plot_profiles.py`. Plots T and n profiles for a list of given MFILES.
- Can now setup the repo in `debug` mode for compilation. See `README.md` for instructions.
- New scan variables - `impurity_ratio(9)` and `fgwsep`.
- New constraint on CS peak field.

### Minor Changes
- Explicitly state 1990 $ for old cost model
- Made photon_wall and rad_fraction global variables, and added calculations to stellarator.
- TF coil documentation now in repository and makefile target `tfdoc`.

## 1.0.14

### Bug Fixes

- Wrong pedestal position used in plot_proc temperature plot (Issue #653) ([957f94a7](https://git.ccfe.ac.uk/process/process/commit/957f94a723b026f67544fa46548bc8a1be062d35))
- Removed hardwired Martin scaling for L-H threshold in plot_proc.py (Issue #679 and #680)
- Fixed error in spherical tokamak divertor geometry calculation (Issue #697)
- Fixed error in spherical tokamak radial build calculation (Issue #704)
- Fixed error in current drive fractions adding to > 1 (Issue #705).
- Fixed issues with uncertainty python utility (Issue #716 #746)
- Fixed issues when there is no inboard blanket (Issue #722 #732)
- Fixed incorrect cross sectional area calculation in resistive TF coils (Issue #727)
- Fixed constraint equations plot in diagnose_process.py (Issue #738)
- Corrected units in resistive TF coil stress output
- Corrected units on ucme and uciac in global variables.
- Fixed issue with plot_proc.py scan counting (Issue #748)
- Fixed issue with run_process.py not working (Issue #766)
- Switched obsolete estotf for estotftgj in stellarator
- Corrected ztot calculation in tfpwr subroutine for resistive TF coils (#773)
- Corrected deltf in sctfcoil.f90 (#779)

### Features

PLASMOD
 - PLASMOD is a 1D transport model which replaces many of the plasma physics calculations in PROCESS. The previous set up remains available.
 - See reference: E. Fable et al., Fusion Engineering and Design, Volume 130, May 2018, Pages 131-136
 - PLASMOD can be run during every PROCESS iteration by setting ipedestal to 3. It can be run just once, at the end of a PROCESS run by setting ipedestal to 2.

 - Created a new file 'physics_functions.f90' to store code moved from physics.f90 which may be used by PLASMOD and other semi-independent models.
 - This is to prevent circular dependencies.
 - Subroutines include: beamcalc, beamfus, imprad, palph, palph2, prad_ipdg89, psync_albajar_fidone, pthresh, radpwr
 - Functions include: bosch_hale, fsv, p_eped_scaling, t_eped_scaling,

 - New user-defined inputs have been added, which all have the prefix 'plasmod_'. These are specific controls and inputs to PLASMOD.
 - For a complete list, see the vardes file. Where appropriate, previously-existing PROCESS input parameters still apply.
 - Certain constraints and iterations variables cannot be used with PLASMOD - see the User Guide for more information.

2D Scan
 - implemented a basic 2-D scan feature in PROCESS.
 - new inputs `scan_dim`, `isweep_2`, `nsweep_2` and `sweep_2`.
 - Does the scan in a basic grid like manner (i.e. jumps to start of next row from end of previous). Would be nice to upgrade to 'zig-zag'-like approach.

Utilities
 - New script compare_radials.py to plot two radial profiles on the same chart for comparison. Takes input columns of data representing the profiles, with the first column being the x-axis, e.g. radial position.
 - evaluate_uncertainties.py now outputs and additional file to allow analysis of failed PROCESS Runs.
 - New script plot_sankey.py to plot a Sankey diagram of the PROCESS power flow
 - New scripts cost_pie.py and cost_bar.py to analyse cost data.
 - New script popcon.py to plot POPCON plot from MFILE.

Miscellaneous
- TF stress in conduit Tresca criterion can now have regular and CEA adjusted options
  (adjustment from [Torre et al. 2016](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7390035)
  paper). (Issue #678)
- [Snipes et al.](2000; http://iopscience.iop.org/article/10.1088/0741-3335/42/5A/336) H-mode threshold scaling options added (Issue #680)
- Initial version of `.gitlab-ci.yml` created for GitLab CI.
- Added Spherical Tokamak and Stellarator examples to the test suite (Issues #715 and #718)
- Output to MFILE variable names for cost models
- Added [Reinke detachment criterion](http://iopscience.iop.org/article/10.1088/1741-4326/aa5145/meta) as constraint equation and formula for tesep (Issue #707)

### Minor changes

- Changed upper bound on `coheof` from 1.0e8 to 5.0e8 (Issue #668).
- A number of changes to `plot_proc.py` and outputs in the fortran associated
  with vertical build. (Merge request !18)
- Update utilities guide for a number of Python utilities
    - `cad_output.py` (Issue #671)
    - `convert_in_dat.py` (Issue #672)
    - `mcnp_output.py` (Issue #674)
    - `output_summary.py` and `output_detailed.py` (Issue #675)
    - `plot_comparison.py` (Merge request !21)
- New Python utility
    - `plot_comparison.py` (Merge request !21)
- File prefixes for input files now works as intended. For example input file called `my_input_IN.DAT`
  will be outputted as `my_input_OUT.DAT` etc.
- `tbrnmn` no longer iteration variable as there is constraint equation 13 and f-value `ftburn` already. `tbrnmn` will act as the constraint limit input value.
- `cdtfleg` no longer an iteration variable.  The outboard leg current density is now calculated for resistive TF coils. (Issue #727)
- `tfacpd` is now calculted for resistive TF coils so is no longer an input.
- Reset test_suite files (Issue #719)
- Added error reporting to function ncore (Issue #735)
- Added input `plasma_res_factor` for adjustment factor for plasma resistivity. Default is 1.0   to preserve old behaviour.
- Added additional scaling factor 'eped_sf' for the EPED pedestal model (pressure and temperature versions).
- Slight change to functionality of utilities/write_new_in_dat.py: This script will no longer create a new IN.DAT from a non-feasible solution. If a scan is run, it will take by default the last feasible solution. If required there is also an option to use the first feasible solution from a scan (Issue #752).
- Made more robust the reading of input files - comments are now denoted only via an asterisk (*), and if a comment is present without an asterisk the reading of the input file will stop (previously it simply ignored constraint equations and iteration variables that could not be read). It is no longer permissible to write an input over multiple lines. Users can now use punctuation in comments as they wish, including full stops and commas.
- Stellarator radial build is output to MFILE (Issue #770)
