# Culham Electron Cyclotron Model

def culecd(self):
        """Routine to calculate Electron Cyclotron current drive efficiency
        author: M R O'Brien, CCFE, Culham Science Centre
        author: P J Knight, CCFE, Culham Science Centre
        effrfss : output real : electron cyclotron current drive efficiency (A/W)
        This routine calculates the current drive parameters for a
        electron cyclotron system, based on the AEA FUS 172 model.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        AEA FUS 172: Physics Assessment for the European Reactor Study
        """
        rrr = 1.0e0 / 3.0e0

        #  Temperature
        tlocal = profiles_module.tprofile(
            rrr,
            physics_variables.rhopedt,
            physics_variables.te0,
            physics_variables.teped,
            physics_variables.tesep,
            physics_variables.alphat,
            physics_variables.tbeta,
        )

        #  Density (10**20 m**-3)
        dlocal = 1.0e-20 * profiles_module.nprofile(
            rrr,
            physics_variables.rhopedn,
            physics_variables.ne0,
            physics_variables.neped,
            physics_variables.nesep,
            physics_variables.alphan,
        )

        #  Inverse aspect ratio
        epsloc = rrr * physics_variables.rminor / physics_variables.rmajor

        #  Effective charge (use average value)
        zlocal = physics_variables.zeff

        #  Coulomb logarithm for ion-electron collisions
        #  (From J. A. Wesson, 'Tokamaks', Clarendon Press, Oxford, p.293)
        coulog = 15.2e0 - 0.5e0 * np.log(dlocal) + np.log(tlocal)

        #  Calculate normalised current drive efficiency at four different
        #  poloidal angles, and average.
        #  cosang = cosine of the poloidal angle at which ECCD takes place
        #         = +1 outside, -1 inside.
        cosang = 1.0e0
        ecgam1 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)
        cosang = 0.5e0
        ecgam2 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)
        cosang = -0.5e0
        ecgam3 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)
        cosang = -1.0e0
        ecgam4 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)

        #  Normalised current drive efficiency (A/W m**-2)
        ecgam = 0.25e0 * (ecgam1 + ecgam2 + ecgam3 + ecgam4)

        #  Current drive efficiency (A/W)
        return ecgam / (dlocal * physics_variables.rmajor)