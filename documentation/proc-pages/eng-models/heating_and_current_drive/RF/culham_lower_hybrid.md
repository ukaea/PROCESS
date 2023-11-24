# Culham Lower Hybrid 

- `iefrf` = 6: Culham Lower Hybrid model[^2]




def cullhy(self):
        """Routine to calculate Lower Hybrid current drive efficiency
        author: P J Knight, CCFE, Culham Science Centre
        effrfss : output real : lower hybrid current drive efficiency (A/W)
        This routine calculates the current drive parameters for a
        lower hybrid system, based on the AEA FUS 172 model.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
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



[^2]: T. C. Hender, M. K. Bevir, M. Cox, R. J. Hastie, P. J. Knight, C. N. Lashmore-Davies, B. Lloyd, G. P. Maddison, A. W. Morris, M. R. O'Brien, M.F. Turner abd H. R. Wilson, *"Physics Assessment for the European Reactor Study"*, AEA Fusion Report AEA FUS 172 (1992)