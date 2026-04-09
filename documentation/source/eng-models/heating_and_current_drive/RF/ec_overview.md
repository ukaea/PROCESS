# Electron Cyclotron Heating | `ElectronCyclotron`

Electron cyclotron resonance heating is the simplest of the radio frequency heating methods. In contrast to ion cyclotron and lower hybrid heating, there is no evanescent region between the antenna and the plasma,although cut-offs can exist within the plasma. 
As a result, the antenna can be retracted into a less hostile environment than for the other schemes. Electron cyclotron heating has been made possible by the invention of the gyrotron millimetre wave source. This and related devices have only emerged since the mid-1970s. Gyrotron tubes capable of 0.5-1 MW out- put and with frequencies in the range 100-200 GHz are now under intense development. The frequency required for a reactor would be in the range 100-200 GHz which corresponds to vacuum wavelengths of 1-2 mm.
Since $\omega_{ce}\le \omega_{pe}$ only the electrons respond to electron cyclotron waves and only the electrons are heated directly. However, under reactor- like conditions the ions will be heated collisionally by the electrons. As the density is increased, a limit is encountered above which electron cyclotron waves cannot penetrate to the central regions of a tokamak. Propagation is again described below

$$
n_{\perp}^2 = 1 - \frac{\omega_{pe}^2}{\omega^2} \ \  (\text{O-mode})
$$

$$
n_{\perp}^2=\frac{(1-\frac{\omega_{pe}^2}{\omega^2}-\frac{\omega_{ce}}{\omega})(1-\frac{\omega_{pe}^2}{\omega^2}+\frac{\omega_{ce}}{\omega})}{(1-\frac{\omega_{pe}^2}{\omega^2}-\frac{\omega_{ce}^2}{\omega^2})} \ \ \text{(X-mode)}
$$

## Normalised current drive efficiency | `eccdef()`

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
\mathtt{mcsq} = m_{\text{e}} \frac{c^2}{1  \text{keV}}
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
### Legendre function and its derivative | `legend()`


| Input       | Description                          |
| :---------- | :----------------------------------- |
| $\mathtt{zlocal}$       |  Local plasma effective charge  |
| $\mathtt{arg}$       |  Argument of Legendre function |

The `legend()` function is a routine that calculates the Legendre function and its derivative. It takes two input parameters: `zlocal` (local plasma effective charge) and `arg` (argument of the Legendre function). The function returns two output values: `palpha` (value of the Legendre function) and `palphap` (derivative of the Legendre function).

Here is the explanation of the `legend()` function:

1. Check if the absolute value of `arg` is greater than `1.0 + 1.0e-10`. If it is, set `eh.fdiags[0]` to `arg` and report an error (error code 18).

2. Set `arg2` to the minimum value between `arg` and `1.0 - 1.0e-10`.

3. Calculate `sinsq` as `0.5 * (1.0 - arg2)`.

4. Calculate `xisq` as `0.25 * (32.0 * zlocal / (zlocal + 1.0) - 1.0)`.

5. Initialize `palpha` to `1.0`, `pold` to `1.0`, `pterm` to `1.0`, `palphap` to `0.0`, and `poldp` to `0.0`.

6. Start a loop that iterates up to 10000 times.

7. Check for convergence every 20 iterations:
    - If `n > 1` and `(n % 20) == 1`, calculate `term1` as `1.0e-10 * max(abs(pold), abs(palpha))` and `term2` as `1.0e-10 * max(abs(poldp), abs(palphap))`.
    - If the absolute difference between `pold` and `palpha` is less than `term1` and the absolute difference between `poldp` and `palphap` is less than `term2`, return `palpha` and `palphap`.

8. Update `pold` to `palpha` and `poldp` to `palphap`.

9. Calculate `pterm` as `pterm * (4.0 * xisq + (2.0 * n - 1.0) ** 2) / (2.0 * n) ** 2 * sinsq`.

10. Update `palpha` as `palpha + pterm`.

11. Update `palphap` as `palphap - n * pterm / (1.0 - arg2)`.

12. If the loop completes without returning, report an error (error code 19).

[^1]: Abramowitz, Milton. *"Abramowitz and stegun: Handbook of mathematical functions."* US Department of Commerce 10 (1972).        