# Electron Cyclotron Heating

Electron cyclotron resonance heating is the simplest of the radio frequency heating methods. In contrast to ion cyclotron and lower hybrid heating, there is no evanescent region between the antenna and the plasma,although cut-offs can exist within the plasma. 
As a result, the antenna can be retracted into a less hostile environment than for the other schemes. Electron cyclotron heating has been made possible by the invention of the gyrotron millimetre wave source. This and related devices have only emerged since the mid-1970s. Gyrotron tubes capable of 0.5-1 MW out- put and with frequencies in the range 100-200 GHz are now under intense development. The frequency required for a reactor would be in the range 100-200 GHz which corresponds to vacuum wavelengths of 1-2 mm.
Since $\omega_{ce}\le \omega_{pe}$ only the electrons respond to electron cyclotron waves and only the electrons are heated directly. However, under reactor- like conditions the ions will be heated collisionally by the electrons. As the density is increased, a limit is encountered above which electron cyclotron waves cannot penetrate to the central regions of a tokamak. Propagation is again described below

$$
n_{\perp}^2 = 1 - \frac{\omega_{pe}^2}{\omega^2} \ \  (\text{O-mode})
$$

$$
n_{\perp}^2=\frac{(1-\frac{\omega_{pe}^2}{\omega^2}-\frac{\omega_{ce}}{\omega})(1-\frac{\omega_{pe}^2}{\omega^2}+\frac{\omega_{ce}}{\omega})}{(1-\frac{\omega_{pe}^2}{\omega^2}-\frac{\omega_{ce}^2}{\omega^2})} \ \ \text{(X-mode)}
$$

## Normalised current drive efficiency `eccdef`

One of the methods for calculating the normalised current drive efficiency is the `eccdef` method found below.

| Input       | Description                          |
| :---------- | :----------------------------------- |
| $\mathtt{tlocal}$       |      Local electron temperature (keV)  |
| $\mathtt{epsloc}$       |  Local inverse aspect ratio |
| $\mathtt{zlocal}$    |     Local plasma effective charge |
| $\mathtt{cosang}$       |  Cosine of the poloidal angle at which ECCD takes place (+1 outside, -1 inside) 
| $\mathtt{coulog}$       | Local coulomb logarithm for ion-electron collisions |

This routine calculates the current drive parameters for a electron cyclotron system, based on the AEA FUS 172 model.
It works out the ECCD efficiency using the formula due to Cohen quoted in the ITER Physics Design Guidelines : 1989 (but including division by the Coulomb Logarithm omitted from IPDG89). 

We have assumed $\gamma^2$-1 << 1, where gamma is the relativistic factor. 
The notation follows that in IPDG89.
The answer ECGAM is the normalised efficiency $n_{\text{e}}IR/P$ with $n_{\text{e}}$ the local density in [$10^{20} / \text{m}^3$], I the driven current in [$\text{MA}$], $R$ the major radius in [$\text{m}$], and $P$ the absorbed power in  [$\text{MW}$].
        

$$
\mathtt{mcsq} = 9.1095\times10^{-31} \frac{c^2}{1  \text{keV}}
$$

$$
\mathtt{f} = 16\left(\frac{\mathtt{tlocal}}{\mathtt{mcsq}}\right)^2
$$

$\mathtt{fp}$ is the derivative of $\mathtt{f}$ with respect to gamma, the relativistic factor, taken equal to $1 + \frac{2T_{\text{e}}}{(m_{\text{e}}c^2)}$

$$
\mathtt{fp} = 16 \left(\frac{\mathtt{tlocal}}{\mathtt{mcsq}}\right)
$$
        
$\mathtt{lam}$ is IPDG89's lambda. `legend` calculates the Legendre function of order $\alpha$ and argument `lam`, `palpha`, and its derivative, `palphap`.
Here `alpha` satisfies $\alpha(\alpha+1) = \frac{-8}{(1+z_{\text{local}})}$. $\alpha$ is of the form  $(-1/2 + ix)$, with x a real number and $i = \sqrt{-1}$.

$$
\mathtt{lam} = 1.0
$$

$$
\mathtt{palpha, palphap} = \text{legend}(\mathtt{zlocal, lam})
$$

$$
\mathtt{lams} = \sqrt{\frac{2 \times \mathtt{epsloc}}{1 + \mathtt{epsloc}}}
$$

$$
\mathtt{palphas} = \text{legend}(\mathtt{zlocal, lams})
$$

$\mathtt{hp}$ is the derivative of IPDG89's \mathtt{h} function with respect to $\mathtt{lam}$

$$
\mathtt{h} = \frac{-4 \times \mathtt{lam}}{\mathtt{zlocal + 5}} \times \frac{1- (\mathtt{lams \times \mathtt{palpha}})}{\mathtt{lam \times \mathtt{palphas}}}
$$

$$
\mathtt{hp} = \frac{-4}{\mathtt{zlocal}+5} \times \left(1- \mathtt{lams} \times \frac{\mathtt{palphap}}{\mathtt{palphas}}
\right)
$$

$\mathtt{facm}$ is IPDG89's momentum conserving factor

$$
\mathtt{facm} = 1.5
$$

$$
\mathtt{y} = \frac{\mathtt{mcsq}}{2 \times \mathtt{tlocal}} \times (1+ \mathtt{epsloc} \times \mathtt{cosang})
$$

We take the negative of the IPDG89 expression to get a positive number

$$
\mathtt{ecgam} = \left(\frac{-7.8 \times \mathtt{facm \times \sqrt{\frac{1 + \mathtt{epsloc}}{1- \mathtt{epsloc}}}}}{\mathtt{coulog}}\right) \times (\mathtt{h} \times \mathtt{fp} -0.5 \times \mathtt{y} \times \mathtt{f} \times \mathtt{hp})
$$

----------------------------------------------------------------------------------
### Legendre function and its derivative `legend`


| Input       | Description                          |
| :---------- | :----------------------------------- |
| $\mathtt{zlocal}$       |  Local plasma effective charge  |
| $\mathtt{arg}$       |  Argument of Legendre function |

``` py
def legend(self, zlocal, arg):
        """Routine to calculate Legendre function and its derivative
        author: M R O'Brien, CCFE, Culham Science Centre
        author: P J Knight, CCFE, Culham Science Centre
        zlocal  : input real : local plasma effective charge
        arg     : input real : argument of Legendre function
        palpha  : output real : value of Legendre function
        palphap : output real : derivative of Legendre function
        This routine calculates the Legendre function `palpha`
        of argument `arg` and order
        `alpha = -0.5 + i sqrt(xisq)``,
        and its derivative `palphap`.
        This Legendre function is a conical function and we use the series
        in `xisq`` given in Abramowitz and Stegun. The
        derivative is calculated from the derivative of this series.
        The derivatives were checked by calculating `palpha` for
        neighboring arguments. The calculation of `palpha` for zero
        argument was checked by comparison with the expression
        `palpha(0) = 1/sqrt(pi) * cos(pi*alpha/2) * gam1 / gam2`
        (Abramowitz and Stegun, eqn 8.6.1). Here `gam1`` and
        `gam2`` are the Gamma functions of arguments
        `0.5*(1+alpha)`` and `0.5*(2+alpha)` respectively.
        Abramowitz and Stegun, equation 8.12.1
        """
        if abs(arg) > (1.0e0 + 1.0e-10):
            eh.fdiags[0] = arg
            eh.report_error(18)

        arg2 = min(arg, (1.0e0 - 1.0e-10))
        sinsq = 0.5e0 * (1.0e0 - arg2)
        xisq = 0.25e0 * (32.0e0 * zlocal / (zlocal + 1.0e0) - 1.0e0)
        palpha = 1.0e0
        pold = 1.0e0
        pterm = 1.0e0
        palphap = 0.0e0
        poldp = 0.0e0

        for n in range(10000):
            #  Check for convergence every 20 iterations

            if (n > 1) and ((n % 20) == 1):
                term1 = 1.0e-10 * max(abs(pold), abs(palpha))
                term2 = 1.0e-10 * max(abs(poldp), abs(palphap))

                if (abs(pold - palpha) < term1) and (abs(poldp - palphap) < term2):
                    return palpha, palphap

                pold = palpha
                poldp = palphap

            pterm = (
                pterm
                * (4.0e0 * xisq + (2.0e0 * n - 1.0e0) ** 2)
                / (2.0e0 * n) ** 2
                * sinsq
            )
            palpha = palpha + pterm
            palphap = palphap - n * pterm / (1.0e0 - arg2)
        else:
            eh.report_error(19)
        
```

[^1]: Abramowitz, Milton. *"Abramowitz and stegun: Handbook of mathematical functions."* US Department of Commerce 10 (1972).        