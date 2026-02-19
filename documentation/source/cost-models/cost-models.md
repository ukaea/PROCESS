# Cost Models

Two cost models are available, determined by the switch `cost_model`.

## 1990 cost model (`cost_model = 0`)

This combines methods[^1] used in the TETRA code [^2] and the Generomak[^3] scheme. The costs are split into accounting categories[^4]. The best references for the algorithms used are[^5], and source file `costs.f90` in the code itself. The majority of the costed items have a unit cost associated with them. These values scale with (for example) power output, volume, component mass etc., and many are available to be changed via the input file. All costs and their algorithms correspond to 1990 dollars.

The unit costs of the components of the fusion power core are relevant to "first-of-a-kind" items. That is to say, the items are assumed to be relatively expensive to build as they are effectively prototypes and specialised tools and machines have perhaps been made specially to create them. However, if a "production line" has been set up, and R & D progress has allowed more experience to be gained in constructing the power core components, the cost will be reduced as a result. Variable `fkind` may be used to multiply the raw unit costs of the fusion power core items (by a factor less than one) to simulate this cost reduction for an *N<sup>th</sup>*-of-a-kind device. In other systems studies of fusion power plants[^6], values for this multiplier have ranged from 0.5 to 0.8.

Many of the unit costs have four possible choices, relating to the level of safety assurance[^7] flag `lsa`. A value `lsa = 1` corresponds to a plant with a full safety credit (i.e. is truly passively safe). Levels 2 and 3 lie between the two extremes, and level 4 corresponds to a present day fission reactor, with no safety credit.

The first wall, blanket, divertor, centrepost (if present) and current drive system have relatively short lifetimes because of their hostile environment, after which they must be replaced. Because of this frequent renewal they can be regarded as though they are "fuel" items, and be costed accordingly. Switch `ifueltyp` is used to control whether this option is used in the code. If `ifueltyp = 1`, the costs of the first wall, blanket, divertor and a fraction `fcdfuel` of the cost of the current drive system are treated as fuel costs. If `ifueltyp = 0`, these are treated as capital costs.

If the switch `ireactor = 0`, no cost of electricity calculation is performed. If `ireactor = 1`, then the cost of electricity is evaluated, with the value quoted in units of $/MWh.

The net electric power is calculated in routine `POWER` It is possible that the net electric power can become negative due to a high recirculating power. Switch `ipnet` determines whether the net electric power is scaled to always remain positive (`ipnet = 0`, or whether it is allowed to become negative (`ipnet = 1`), in which case no cost of electricity calculation is performed.

## 2015 Kovari model (`cost_model = 1`)

This model[^8] provides only capital cost, and it is not currently suitable for estimating the cost of electricity. *N<sup>th</sup>*-of-a-kind factors, level of safety assurance factors, and blanket replacement costs are not included. The mean electric output is calculated using the capacity factor, which takes account of the availability and the dwell time for a pulsed reactor. The capital cost divided by the mean electric output is a useful comparison parameter.

[^1]: R. L. Reid and Y-K. M. Peng, *"Potential Minimum Cost of Electricity of Superconducting Coil Tokamak Power Reactors"*, Proceedings of 13th IEEE Symposium on Fusion Engineering, Knoxville, Tennessee, October 1989, p. 258
[^2]: R. L. Reid et al., *"ETR/ITER Systems Code"*, Oak Ridge Report ORNL/FEDC-87/7 (1988)
[^3]: J. Sheeld et al., *"Cost Assessment of a Generic Magnetic Fusion Reactor"*, Fusion Technology **9** (1986) 199
[^4]: S. Thompson, *"Systems Code Cost Accounting"*, memo FEDC-M-88-SE,-004 (1988)
[^5]: J. D. Galambos, *"STAR Code : Spherical Tokamak Analysis and Reactor Code"*, Unpublished internal Oak Ridge document. A copy exists in the `PROCESS` Project Work File [^9].
[^6]: J. D. Galambos, L. J. Perkins, S. W. Haney and J. Mandrekas, *Nuclear Fusion*, **35** (1995) 551
[^7]: J. P. Holdren et al., *"Report of the Senior Committee on Environmental Safety and Economic Aspects of Magnetic Fusion Energy"*, Fusion Technology, **13** (1988) 7
[^8]: M. Kovari et al., *"The cost of a fusion power plant: extrapolation from ITER"*, in preparation (2015)
[^9]: P. J. Knight, *"PROCESS Reactor Systems Code"*, AEA Fusion Project Work File, F/RS/CIRE5523/PWF (1992)