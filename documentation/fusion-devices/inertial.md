# Inertial Fusion Energy Model

As well as magnetic confinement devices, `PROCESS` has the ability to model inertial fusion plants, in which a laser or ion beam is used to ignite a target pellet containing the fusion fuel.

To activate the inertial fusion energy (IFE) coding, it is necessary to create a file `device.dat`, containing the single character 3 in the first row, in the working directory. This has the effect of setting the internally-used switch `ife = 1`. If the file is absent, or its first character is set to something other than 3, the IFE model is not used, and `ife` is set to 0.

The IFE model[^1] is controlled using two additional switches.

`ifetyp = 0` : Generic device build

`ifetyp = 1` : OSIRIS-type device build[^2] [^3] [^4]

`ifetyp = 2` : SOMBRERO-type device build[^5] [^6]

`ifetyp = 3` : HYLIFE-II-type device build[^7] [^8] [^9]

`iftype = 4` : 2019 device build

Switch `ifetyp` defines the type of device that is assumed; this varies widely between different conceptual designs. The generic type assumes a cylindrical symmetric device, while the other types are approximations to the builds of the given conceptual machines[^10]. In general, the build from the centre of the device (at the target ignition location) is in the order: chamber, first wall, gap, blanket, gap, shield, gap, building wall. The user specifies the thicknesses of these regions, and also the materials that are present and in what proportions.

`ifedrv = -1` : Driver efficiency and target gain are input as functions of driver energy

`ifedrv = 0` : Driver efficiency and target gain are input

`ifedrv = 1` : SOMBRERO laser drive efficiency and target gin assumed

`ifedrv = 2` : OSIRIS heavy ion beam driver efficiency and target gain are assumed[^11]

Switch `ifedrv` defines how the code calculates the drivers efficiency and target gain - these are the primary outputs required from the physics part of the model. For the SOMBRERO and OSIRIS cases (`ifedrv = 1` and `ifedrv = 2`, respectively) the driver efficiency and gain are calculated from curves of these parameters as functions of the driver energy, via the two arrays`etaxe(1:10)` and `gainve(1:10)` respectively; the element number corresponds to the driver energy in MJ, and outside the range 1-10 MJ the curves are extrapolated linearly. Finally, for the `ifedrv = 0` case, the user inputs single values for the driver efficiency (`drveff`) and target gain (`tgain`).

Constraint equation no. 50 can be turned on to enable the ignition repetition rate to remain below a user-specified upper limit (`rrmax`); iteration variable no. 86 (`frrmax`) is the associated f-value. The other iteration variables relevant for the IFE model are nos. 81-85 (`edrive`, `drveff`, `tgain`, `chrad` and `pdrive`).

[^1]: P. J. Knight, *"PROCESS 3009: Incorporation of Inertial Fusion Energy Model"*, Work File Note F/MI/PJK/PROCESS/CODE/032
[^2]: Bourque et al., *"Overview of the OSIRIS IFE Reactor Conceptual Design"*, Fusion Technology **21** (1992) 1465
[^3]: Meier and Bieri, *"Economic Modeling and Parametric Studies for OSIRIS - a HIB-Driven IFE Power Plant"*, Fusion Technology **21** (1992) 1547
[^4]: Ghose et al., *"BOP Designs for OSIRIS and SOMBRERO IFE Reactor Plants"*, Fusion Technology **21** (1992) 1501
[^5]: Sviatoslavsky et al., *"A KrF Laser Driven Inertial Fusion Reactor SOMBRERO"*, Fusion Technology **21** (1992) 1470
[^6]: Meier and Bieri, *"Economic Modeling and Parametric Studies for SOMBRERO - a Laser-Driven IFE Power Plant"*, Fusion Technology **21** (1992) 1552
[^7]: Moir et al., *"HYLIFE-II: A Molten-Salt Inertial Fusion Energy Power Plant Design | Final Report"*, Fusion Technology **25** (1994) 5
[^8]: Moir, *"HYLIFE-II Inertial Fusion Energy Power Plant Design"*, Fusion Technology **21** (1992) 1475
[^9]: Homan and Lee, *"Performance and Cost of the HYLIFE-II Balance of Plant"*, Fusion Technology **21** (1992) 1475
[^10]: P. J. Knight, *"PROCESS IFE Build Details"*, F/MI/PJK/LOGBOOK12, pp. 52, 53, 56, 57
[^11]: Bieri and Meier, *"Heavy-Ion Driver Parametric Studies and Choice of a Base 5 MJ Driver Design"*, Fusion Technology **21** (1992) 1557