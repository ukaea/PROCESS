# Culham Electron Cyclotron Model | `culecd()`

- `i_hcd_primary/i_hcd_secondary` = 7

This routine calculates the current drive parameters for a electron cyclotron system, based on the AEA FUS 172 model[^1]

1. Local electron temperature $(\mathtt{tlocal})$ is calculated using the `teprofile` method  
2. Local electron density $(\mathtt{dlocal})$ is calculated using the `neprofile` method

3. Calculate the inverse aspect ratio `epsloc`.

$$
\mathtt{epsloc} = \frac{1}{3A}
$$

4. Calculate the Coulomb logarithm for ion-electron collisions `coulog`.[^2]

$$
\mathtt{coulog} = 15.2 - 0.5\log({\mathtt{dlocal}}) + \log({\mathtt{tlocal}})
$$

Calculate normalised current drive efficiency at four different poloidal angles, and average.
cosang = cosine of the poloidal angle at which ECCD takes place = +1 outside, -1 inside.

## Normalised current drive efficiency

Uses the `eccdef` model found [here](ec_overview.md)

5. Calculate the normalised current drive efficiency at four different poloidal angles, and average.
         - Set `cosang` to 1.0 and calculate `ecgam1`.
         - Set `cosang` to 0.5 and calculate `ecgam2`.
         - Set `cosang` to -0.5 and calculate `ecgam3`.
         - Set `cosang` to -1.0 and calculate `ecgam4`.

        cosang = 1.0e0
        ecgam1 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)
        cosang = 0.5e0
        ecgam2 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)
        cosang = -0.5e0
        ecgam3 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)
        cosang = -1.0e0
        ecgam4 = self.eccdef(tlocal, epsloc, zlocal, cosang, coulog)    

6. Calculate the normalised current drive efficiency `ecgam` as the average of `ecgam1`, `ecgam2`, `ecgam3`, and `ecgam4`.

$$
\mathtt{ecgam} = 0.25(\mathtt{ecgam1} + \mathtt{ecgam2} +\mathtt{ecgam3} + \mathtt{ecgam4})
 $$

7. Calculate the current drive efficiency by dividing `ecgam` by `(dlocal * physics_variables.rmajor)`.

$$
\text{Current drive efficiency [A/W]} = \frac{\mathtt{ecgam}}{\mathtt{dlocal} \times R_0}
$$

Note: The `eccdef` method is called to calculate the current drive efficiency at each poloidal angle.

[^1]: Hender, T.C., Bevir, M.K., Cox, M., Hastie, R.J., Knight, P.J., Lashmore-Davies, C.N., Lloyd, B., Maddison, G.P., Morris, A.W., O’Brien, M.R. and Turner, M.F., 1992. *"Physics assessment for the European reactor study."* AEA FUS, 172.

[^2]: Wesson, J. and Campbell, D.J., 2011. *"Tokamaks"* (Vol. 149). Oxford university press.
