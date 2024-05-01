# Culham Lower Hybrid 

- `iefrf` = 6: Culham Lower Hybrid model[^1]

This routine calculates the current drive parameters for a
lower hybrid system, based on the AEA FUS 172 model.
AEA FUS 251: A User's Guide to the PROCESS Systems Code
AEA FUS 172: Physics Assessment for the European Reactor Study

        rratio = self.lhrad()
        rpenet = rratio * physics_variables.rminor

        # Local density, temperature, toroidal field at this minor radius

        dlocal = 1.0e-19 * profiles_module.nprofile(
            rratio,
            physics_variables.rhopedn,
            physics_variables.ne0,
            physics_variables.neped,
            physics_variables.nesep,
            physics_variables.alphan,
        )
        tlocal = profiles_module.tprofile(
            rratio,
            physics_variables.rhopedt,
            physics_variables.te0,
            physics_variables.teped,
            physics_variables.tesep,
            physics_variables.alphat,
            physics_variables.tbeta,
        )
        blocal = (
            physics_variables.bt
            * physics_variables.rmajor
            / (physics_variables.rmajor - rpenet)
        )  # Calculated on inboard side

        # Parallel refractive index needed for plasma access

        frac = np.sqrt(dlocal) / blocal
        nplacc = frac + np.sqrt(1.0e0 + frac * frac)

        # Local inverse aspect ratio

        epslh = rpenet / physics_variables.rmajor

        # LH normalised efficiency (A/W m**-2)

        x = 24.0e0 / (nplacc * np.sqrt(tlocal))

        term01 = 6.1e0 / (nplacc * nplacc * (physics_variables.zeff + 5.0e0))
        term02 = 1.0e0 + (tlocal / 25.0e0) ** 1.16e0
        term03 = epslh**0.77e0 * np.sqrt(12.25e0 + x * x)
        term04 = 3.5e0 * epslh**0.77e0 + x

        if term03 > term04:
            eh.fdiags[0] = term03
            eh.fdiags[1] = term04
            eh.report_error(129)

        gamlh = term01 * term02 * (1.0e0 - term03 / term04)

        # Current drive efficiency (A/W)

        return gamlh / ((0.1e0 * dlocal) * physics_variables.rmajor)


    def lhrad(self):

    """Routine to calculate Lower Hybrid wave absorption radius
    author: P J Knight, CCFE, Culham Science Centre
    rratio  : output real : minor radius of penetration / rminor
    This routine determines numerically the minor radius at which the
    damping of Lower Hybrid waves occurs, using a Newton-Raphson method.
    AEA FUS 251: A User's Guide to the PROCESS Systems Code
    AEA FUS 172: Physics Assessment for the European Reactor Study
    """
    #  Correction to refractive index (kept within valid bounds)
    drfind = min(0.7e0, max(0.1e0, 12.5e0 / physics_variables.te0))

    #  Use Newton-Raphson method to establish the correct minor radius
    #  ratio. g is calculated as a function of r / r_minor, where g is
    #  the difference between the results of the two formulae for the
    #  energy E given in AEA FUS 172, p.58. The required minor radius
    #  ratio has been found when g is sufficiently close to zero.

    #  Initial guess for the minor radius ratio

    rat0 = 0.8e0

    for _ in range(100):
        #  Minor radius ratios either side of the latest guess

        r1 = rat0 - 1.0e-3 * rat0
        r2 = rat0 + 1.0e-3 * rat0

        #  Evaluate g at rat0, r1, r2

        g0 = self.lheval(drfind, rat0)
        g1 = self.lheval(drfind, r1)
        g2 = self.lheval(drfind, r2)

        #  Calculate gradient of g with respect to minor radius ratio

        dgdr = (g2 - g1) / (r2 - r1)

        #  New approximation

        rat1 = rat0 - g0 / dgdr

        #  Force this approximation to lie within bounds

        rat1 = max(0.0001e0, rat1)
        rat1 = min(0.9999e0, rat1)

        if abs(g0) <= 0.01e0:
            break
        rat0 = rat1

    else:
        eh.report_error(16)
        rat0 = 0.8e0

    return rat0

[^1]: T. C. Hender, M. K. Bevir, M. Cox, R. J. Hastie, P. J. Knight, C. N. Lashmore-Davies, B. Lloyd, G. P. Maddison, A. W. Morris, M. R. O'Brien, M.F. Turner abd H. R. Wilson, *"Physics Assessment for the European Reactor Study"*, AEA Fusion Report AEA FUS 172 (1992)