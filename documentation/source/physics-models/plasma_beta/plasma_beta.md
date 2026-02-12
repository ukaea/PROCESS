# Plasma Beta

## Overview

The efficiency of confinement of plasma pressure by the magnetic field is represented by the ratio:

$$
\beta(\rho) = \frac{2\mu_0p(\rho)}{\left(B(\rho)\right)^2}
$$

Where $\beta$ is a function of normalised minor radius across the plasma ($\rho$), due to the change in pressure and magnetic field strength. 

The standard $\beta$ term used for comparison and to represent the plasma as a whole in many calculations is the volume averaged value given by:

$$
\langle \beta \rangle = \frac{2\mu_0 \langle p \rangle}{\langle B \rangle^2}
$$

Where $\langle p \rangle$ is the volume averaged plasma pressure and $\langle B \rangle$ is the average field.

There are several different measures of this type, arising from different choices of definition and from the need to quantify different equilibrium properties.

In its expanded form of magnetic field components:

$$
\beta = \frac{2\mu_0p}{B^2_{\text{T}}+B^2_{\text{P}}}
$$

This is often broken down into separate related quantities knows as the toroidal $\beta_{\text{t}}$ and the poloidal beta $\beta_{\text{p}}$ whose definitions are:

$$
\beta_{\text{t}} = \frac{2\mu_0p}{B^2_{\text{T}}} \\
\beta_{\text{p}} = \frac{2\mu_0p}{B^2_{\text{P}}}
$$

The relationship between these quantities can be written as:

$$
\frac{1}{\beta} = \frac{1}{\beta_{\text{t}}} + \frac{1}{\beta_{\text{p}}}
$$

implying that the total $\beta$ is dominated by the smaller of the two contributions.Observe that by definition $\beta \le 1$. However,either $\beta_{\text{t}}$ or $\beta_{\text{p}}$, but not both, can be greater than unity.

The above in its simplest form is known as the total thermal beta $\beta_{\text{th}}$ as it only takes into account the pressure produced by the ions and electrons such that:

$$
\beta_{\text{th}} = \frac{2\mu_0 \langle n_{\text{e}}T_{\text{e}}+n_{\text{i}}T_{\text{i}}\rangle}{B^2}
$$

In future reactors there will be a notable pressure contribution due to the productions of fast alpha particles and that from plasma-beam interactions if used. So the true total beta can be defined as:

$$
\beta_{\text{tot}} = \frac{2\mu_0 \langle n_{\text{e}}T_{\text{e}}+n_{\text{i}}T_{\text{i}} + p_{\alpha}+ p_{\text{beam}}\rangle}{B^2} \\
\beta_{\text{tot}} = \beta_{\text{th}}+\beta_{\alpha}+ \beta_{\text{beam}}
$$

Models for the fast alpha particle pressure contribution can be found [here](plasma_alpha_beta_contribution.md).

The calculation of the beta component given by neutral beams is calculated in the neutral beam fusion calculations in [`beam_fusion()`](../fusion_reactions/beam_reactions.md#beam-slowing-down-properties--beam_fusion). A description can be found [here](../fusion_reactions/beam_reactions.md#derivation-of-beam-slowing-down-rate-and-critical-energy).

----------------

## Derivation of plasma beta parameter

From the lowest order momentum balance equation:

$$
\nabla p = \mathbf{j} \times \mathbf{B}
$$

as it pertains to plasma equilibrium. The current and the field must also satisfy Maxwell’s equations

$$
\mu_0\mathbf{j} = \nabla \times \mathbf{B},  \quad \nabla \cdot \mathbf{B} = 0
$$

A consequence of the above is that the current and  magnetic field lie on isobaric surfaces. This follows directly from the observation that $\mathbf{j} \cdot \nabla p = \mathbf{B} \cdot \nabla p = 0$ and $\nabla p$ is everywhere normal to the surface $p = \text{const}$. The force is everywhere normal to the isobaric surface and just balances the pressure gradient force, $-\nabla p$. Although the current and the magnetic field lie in a common flux surface, they can only be parallel in regions where the pressure gradient vanishes. Since currents that are parallel to the field do not contribute to the $\mathbf{j} \times \mathbf{B}$ force, they are known as “force-free currents.”

Taking the divergence of $\mu_0\mathbf{j} = \nabla \times \mathbf{B}$  yields $\nabla \cdot \mathbf{j} = 0$ which is consistent with the quasineutrality assumption. The vector product of $\mathbf{B}$ with the momentum balance equation leads to an expression for the current perpendicular to $\mathbf{B}$,

$$
\mathbf{j}_{\perp} = \frac{\mathbf{B} \times \nabla p}{B^2}
$$

We see that there can be no perpendicular current in the absence of a pressure gradient. This leads to a momentum balance expressed in terms of the pressure and magnetic field.

$$
\nabla \left(p+\frac{B^2}{2\mu_0}\right) = \frac{1}{\mu_0}\left(\mathbf{B} \cdot \nabla \right)\mathbf{B}
$$

The right hand side vanishes when the field lines are straight and parallel, in which case above reduces to a simple statement that the total (kinetic plus magnetic) pressure is constant everywhere within a confined plasma, 

$$
p + \frac{B^2}{2\mu_0} \equiv \frac{B^2}{2\mu_0}\left(1+\beta \right) = \text{const}
$$

where we have introduced the quantity

$$
\beta \equiv \frac{2\mu_0p}{B^2}
$$

------------------------

## Definitions

### Volume averaged thermal toroidal beta

We define $B_{\text{T,on-axis}}$ as the toroidal field at the plasma major radius, $R_0$

$$
\overbrace{\langle \beta_t \rangle_{\text{V}}}^{\texttt{beta_toroidal_thermal_vol_avg}} = \frac{2\mu_0 \overbrace{\langle p_{\text{thermal}} \rangle}^{\texttt{pres_plasma_thermal_vol_avg}}}{\underbrace{B_{\text{T,on-axis}}^2}_{\texttt{b_plasma_toroidal_on_axis}}}
$$

### Volume averaged thermal poloidal beta



$$
\overbrace{\langle \beta_p \rangle_{\text{V}}}^{\texttt{beta_poloidal_thermal_vol_avg}} = \frac{2\mu_0 \overbrace{\langle p_{\text{thermal}} \rangle}^{\texttt{pres_plasma_thermal_vol_avg}}}{\underbrace{\langle B_{\text{P,average}} \rangle^2}_{\texttt{b_plasma_poloidal_average}}}
$$

### Volume averaged thermal beta

$$
\overbrace{\langle \beta \rangle_{\text{V}}}^{\texttt{beta_thermal_vol_avg}} = \frac{2\mu_0 \overbrace{\langle p_{\text{thermal}} \rangle}^{\texttt{pres_plasma_thermal_vol_avg}}}{\sqrt{\langle B_{\text{P,average}} \rangle^2+B_{\text{T,on-axis}}^2}}
$$

------------------

## Troyon Beta Limit

The Troyon plasma beta limit is given by[^0][^1]:

$$\begin{aligned}
\beta < 0.01\, g \, \frac{I \ [\mbox{MA}]}{a \ [\mbox{m}] \, B_0 \ [\mbox{T}]}
\end{aligned}$$

where $B_0$ is the axial vacuum toroidal field. The beta
coefficient $g$ is set using input parameter `beta_norm_max`. To apply the Troyon limit please see [Beta upper limit](#beta-upper-limit).

By default, $\beta$ is defined with respect to the total equilibrium B-field. This can be changed depending on the setting of `i_beta_component` seen below:

| `i_beta_component` | Description |
| :-: | - |
| 0 (default) | Apply the $\beta$ limit to the total plasma beta (including the contribution from fast alphas and  neutral beams) |
| 1 | Apply the $\beta$ limit to only the thermal component of beta |
| 2 | Apply the $\beta$ limit to only the thermal plus neutral beam contributions to beta |
| 3 | Apply the $\beta$ limit to the total toroidal beta (including the contribution from fast alphas and neutral beams) |

------------

### Setting the Beta g Coefficient

Switch `i_beta_norm_max` determines how the beta $g$ coefficient `beta_norm_max` should
be calculated. The following switch options are available below:

#### User Input

The user can specify the maximum allowed value of $\beta_{\text{N}}$ directly by stating `i_beta_norm_max = 0` in the input file.

```python
IN.DAT

i_beta_norm_max = 0
beta_norm_max = 3.0


```

---------

#### Wesson Relation | `calculate_beta_norm_max_wesson()`

This can be activated by stating `i_beta_norm_max = 1` in the input file.

`beta_norm_max` is set to `beta_norm_max_wesson` using:  

$$
\texttt{beta_norm_max_wesson} = g = 4l_i
$$

This relation is based off of data taken from DIII-D shots for $\beta_{\text{N}} \ge 2.5$[^7].

This is only recommended for high aspect ratio tokamaks[^3].

**It is recommended to use this switch with [`i_alphaj = 1`](../plasma_current/plasma_current.md#wesson-relation) and [`i_ind_plasma_internal_norm = 1`](../plasma_current/plasma_inductance.md#wesson-relation) as they are self-consistent with each other.**

---------

#### Original Scaling Law | `calculate_beta_norm_max_original()`

This can be activated by stating `i_beta_norm_max = 2` in the input file.

`beta_norm_max` is set to `beta_norm_max_original_scaling` using:  

$$
\texttt{beta_norm_max_original_scaling} = g =2.7(1+5\epsilon^{3.5})
$$ 

<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>Original Normalized Beta Limit</title>
    <style>
      html, body {
        box-sizing: border-box;
        display: flow-root;
        height: 100%;
        margin: 0;
        padding: 0;
      }
    </style>
    <script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-3.6.0.min.js"></script>
    <script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-mathjax-3.6.0.min.js"></script>
    <script type="text/javascript">
        Bokeh.set_log_level("info");
    </script>
  </head>
  <body>
    <div id="cf29dd69-dc18-43d8-9476-146578601b91" data-root-id="p1004" style="display: contents;"></div>
  
    <script type="application/json" id="ed3c7b85-c675-406a-a3ff-f1b0785e8266">
      {"62884a47-c32c-4b55-a144-365023ac52d1":{"version":"3.6.0","title":"Bokeh Application","roots":[{"type":"object","name":"Figure","id":"p1004","attributes":{"width":400,"height":400,"x_range":{"type":"object","name":"Range1d","id":"p1014","attributes":{"start":1,"end":5}},"y_range":{"type":"object","name":"Range1d","id":"p1015","attributes":{"start":2,"end":15}},"x_scale":{"type":"object","name":"LinearScale","id":"p1016"},"y_scale":{"type":"object","name":"LinearScale","id":"p1017"},"title":{"type":"object","name":"Title","id":"p1007","attributes":{"text":"Original Normalized Beta Limit"}},"renderers":[{"type":"object","name":"GlyphRenderer","id":"p1047","attributes":{"data_source":{"type":"object","name":"ColumnDataSource","id":"p1001","attributes":{"selected":{"type":"object","name":"Selection","id":"p1002","attributes":{"indices":[],"line_indices":[]}},"selection_policy":{"type":"object","name":"UnionRenderers","id":"p1003"},"data":{"type":"map","entries":[["x",{"type":"ndarray","array":{"type":"bytes","data":"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"},"shape":[500],"dtype":"float64","order":"little"}],["y",{"type":"ndarray","array":{"type":"bytes","data":"NDMzMzMzMEC5AilM7KcvQJVzml0k8C5ARd2ThsU+LkDCBkBaipMtQPJQkN0w7ixAqzkkVXpOLEDTYmQWK7QrQJIFl1sKHytAoD64GuKOKkA1w+TefgMqQPhXKaSvfClAn9yMtUX6KEAd7i2NFHwoQDkOULbxAShAYu42sbSLJ0AX9LDYNhknQGtRNElTqiZApCN0yeY+JkCv9VS0z9YlQKjJKeTtcSVA7m4knyIQJUCJa+WEULEkQJEjGX1bVSRAmjERpyj8I0DtEkpKnqUjQItizsejUSNAoeFpjCEAI0AOfJ8DAbEiQO5VVossZCJA6rwyaI8ZIkDplpC6FdEhQKeYFXSsiiFA1TLSTUFGIUCuvOi+wgMhQGTwsvMfwyBAU1FexUiEIECtjfixLUcgQHpf5dS/CyBAbrNvv+GjH0DO4tImZjMfQPva0nPyxR5AHHz/P21bHkBaFf0evvMdQNgjfJPNjh1AO3K9BIUsHUAQ2Jq0zswcQMNQDbaVbxxATpIp5MUUHEActYzZS7wbQIPgMugUZhtAXkuyEQ8SG0D6NdX/KMAaQCbTjf1RcBpAUGE/8HkiGkA3+1dRkdYZQIvpNiiJjBlA+HxbBFNEGUD9sdj34P0YQCcWCpIluRhAwpmG2hN2GECGKE5MnzQYQAYRL9G79BdAzGxgvV22F0DD4k7LeXkXQMFCmRcFPhdAT5w6HfUDF0Bgk9+xP8sWQCXVZALbkxZA+Kt8j71dFkAGzXkq3igWQDOUPfIz9RVA8fpHULbCFUAbruf1XJEVQO67iNkfYRVANWggNPcxFUChx7R+2wMVQJHV/m/F1hRAnsgl+q2qFEB2epJIjn8UQMbG2b1fVRRAfdO88RssFEAbQj6vvAMUQNFWy/I73BNASy546JO1E0BxJk7qvo8TQHOpq363ahNAk5S0VnhGE0DvfdJM/CITQN4lRGM+ABNA72i7wjneEkAlEAm56bwSQK3k1bdJnBJAvXNoU1V8EkBt9nZBCF0SQObXBFhePhJAmVpLjFMgEkBR46zx4wISQGZ2srgL5hFA+/gSLsfJEUAVzcO5Eq4RQGtkEt7qkhFARGrGNkx4EUAvKUx4M14RQIvV5m6dRBFArmnq/YYrEUCAxPwe7RIRQI6+XeHM+hBAZe01aSPjEED6z+vu7csQQB0if74ptRBAyRfqNtSeEEAFRIjJ6ogQQJvxgvlqcxBAfbZCW1JeEEATDeaTnkkQQPrAvFhNNRBA4v7HblwhEEBZ2T6qyQ0QQIssLtwl9Q9AWzUkV2zPD0CdmJvDYqoPQE2EsjoFhg9AWj2a709iD0CT4M4uPz8PQEnrVV3PHA9AkkoD+Pz6DkBfssSSxNkOQA8A89ciuQ5APm+phxSZDkAGaSJ3lnkOQDS5GZClWg5A+/Uz0D48DkD46WpIXx4OQAvRfhwEAQ5AqztsgirkDUADbebBz8cNQOUK1jPxqw1AV/fbQYyQDUD5LdhlnnUNQA+AdCklWw1AiQyzJR5BDUB5UYAChycNQBe3SHZdDg1AeXSRRZ/1DEBQsJRCSt0MQIDA4ExcxQxA/Gz6UNOtDEDvGgJIrZYMQK/GWzfofwxAm7NZMIJpDEAuuulPeVMMQE0dRb7LPQxAz9CirncoDEAUHOxeexMMQCeEcxfV/gtA0+mtKoPqC0DByO30g9YLQE+FINzVwgtArLeNT3evC0A2Y5jHZpwLQOEJgsWiiQtA6osv0yl3C0C0xO+C+mQLQD/WQ28TUwtAKhWpOnNBC0Cuh2SPGDALQIHqTx8CHwtAFS+ooy4OC0D0Zt3cnP0KQJ0QZJJL7QpAeLqHkjndCkAU8T6yZc0KQA5uAM3OvQpAln2ZxHOuCkCgkAWBU58KQH7yRvBskApAgphABr+BCkAUA5G8SHMKQIonbhIJZQpArVqCDP9WCkDSM8q0KUkKQA1hcxqIOwpA6GS8URkuCkCTNtVz3CAKQIi9wJ7QEwpAIyI39fQGCkBv7IieSPoJQEDrgsbK7QlAO91SnXrhCUBC1mxXV9UJQGpbcS1gyQlAITAUXJS9CUAuzwMk87EJQI6K0cl7pglADE3alS2bCUAL+S/UB5AJQKJfg9QJhQlA2MsO6jJ6CUBlHYFrgm8JQANv6bL3ZAlAIUSjHZJaCUAYOkMMUVAJQAs5hOIzRglA3SA1Bzo8CUCI7ibkYjIJQJFWG+atKAlAENGzfBofCUBDFGEaqBUJQFn6UjRWDAlAmM9oQiQDCUDSBSK/EfoIQFJKjyce8QhAkvtD+0joCEDv+0e8kd8IQOPeCe/31ghALm5RGnvOCEByhDLHGsYIQBk7AIHWvQhA9GhA1a21CECYb59ToK0IQCxV5I2tpQhAoCjlF9WdCEBRrnuHFpYIQPhSenRxjghAQ2OheOWGCED1hZQvcn8IQPl20DYXeAhAhAGhLdRwCEC3Nxe1qGkIQP3l/2+UYghAq0DaApdbCEA4ys4TsFQIQLdwpkrfTQhA5ODBUCRHCEChDRHRfkAIQDvqCnjuOQhAYVal83IzCEBTOk3zCy0IQCnS3ie5JghA1iaeQ3ogCEDqsy/6ThoIQK04kQA3FAhAsLMSDTIOCECTh0/XPwgIQAfIJxhgAghABa65iZL8B0AqMlvn1vYHQFDMk+0s8QdAc1cWWpTrB0DoGLvrDOYHQBXqeWKW4AdAuoNkfzDbB0AF6qAE29UHQI74Y7WV0AdAgQ3sVWDLB0Ah03urOsYHQOcmVXwkwQdAiB20jx28B0AmI8qtJbcHQPU2uZ88sgdAt0GPL2KtB0BYhkEolqgHQAssqFXYowdAT+F5hCifB0A1l0eChpoHQGBUeB3ylQdAGx9FJWuRB0AM/rRp8YwHQNsOmbuEiAdAZ7KI7CSEB0Dozd3O0X8HQJ0gsTWLewdAYa3W9FB3B0DTN9rgInMHQJHU+84AbwdA8osslepqB0APDwsK4GYHQGd+4AThYgdA7UGdXe1eB0D28dXsBFsHQKxQwIsnVwdAoVMwFFVTB0A1PZVgjU8HQE3F9kvQSwdAHVHysR1IB0CpObhudUQHQJggCV/XQAdAFFMzYEM9B0BYOhBQuTkHQKjZAQ05NgdAZFnwdcIyB0DRnkdqVS8HQHvw9MnxKwdAwaZkdZcoB0BP6H9NRiUHQE9yqjP+IQdA+GvACb8eB0BRRRSyiBsHQMigbA9bGAdAj0cCBTYVB0BFKH52GRIHQPNf90cFDwdA20zxXfkLB0Arq1md9QgHQCG7huv5BQdAg3A1LgYDB0BIq4dLGgAHQBF5Aio2/QZAcF+MsFn6BkCyrmvGhPcGQP/cRFO39AZApukYP/HxBkBpyENyMu8GQJjUetV67AZA2UvLUcrpBkBd0JjQIOcGQHfymzt+5AZAW8HgfOLhBkDeYsV+Td8GQAuy+Cu/3AZAi+R4bzfaBkB9NpI0ttcGQNyc3WY71QZAJ34/8sbSBkApcebCWNAGQOcBSsXwzQZAWHwp5o7LBkD+vIoSM8kGQCAHuTfdxgZAluBDQ43EBkD98v0iQ8IGQEjy+8T+vwZAdoiTF8C9BkCARloJh7sGQCOaJIlTuQZArMgEhiW3BkCD7knv/LQGQGcDf7TZsgZAUeNpxbuwBkDCWwoSo64GQJo9mYqPrAZAI3OHH4GqBkBzGn3Bd6gGQOWjWGFzpgZArfQt8HOkBkBijUVfeaIGQHO0G6CDoAZAeqRfpJKeBkA2vvJdppwGQFe+576+mgZAv/aBuduYBkBpizRA/ZYGQLSyoUUjlQZAIvmZvE2TBkBZiBuYfJEGQHNxUcuvjwZAhvqSSeeNBkA572IGI4wGQIn0bvViigZAct+OCqeIBkCrDsQ574YGQCvHOHc7hQZAn5M/t4uDBkCVplLu34EGQG0/ExE4gAZA9xFJFJR+BkCvsOHs83wGQIv5749XewZAV4Wr8r55BkCCGXAKKngGQGYcvcyYdgZA4Qs1Lwt1BkBc9pwngXMGQAT226v6cQZAXa76sXdwBkDxyyIw+G4GQEKGnhx8bQZAzCPYbQNsBkAjgFkajmoGQBmUyxgcaQZA9P/1X61nBkCKl77mQWYGQFXwKKTZZAZAcPFVj3RjBkByZYOfEmIGQAuOC8yzYAZAjblkDFhfBkAb2iBY/10GQJoe7aapXAZAY42R8FZbBkCOoPAsB1oGQOPjBlS6WAZAcZTqXXBXBkCuQctCKVYGQCFw8frkVAZAnj2+fqNTBkDyBqvGZFIGQBIPScsoUQZAqydBhe9PBkAvW1PtuE4GQDSYVvyETQZAPF44q1NMBkDDa/zyJEsGQLhtvMz4SQZAJrCnMc9IBkAq0AIbqEcGQDJvJ4KDRgZAZueDYGFFBkBKAZuvQUQGQJSqA2kkQwZAHq5ohglCBkAKbYgB8UAGQPWYNNTaPwZAUe9R+MY+BkDF9ddntT0GQLG30BymPAZAooRYEZk7BkDir50/jjoGQABR4KGFOQZAVQVyMn84BkCCsrXrejcGQOhJH8h4NgZAAo0zwng1BkC+0ofUejQGQKrNwfl+MwZAFlOXLIUyBkAJI85njTEGQBuxO6aXMAZAIu7E4qMvBkC1El4Ysi4GQH5qCkLCLQZAYiDcWtQsBkBkC/Rd6CsGQF18gUb+KgZAcgzCDxYqBkA/bAG1LykGQNgzmTFLKAZAYrPwgGgnBkCFxHyehyYGQHacv4WoJQZAvZ5IMsskBkCqMLSf7yMGQG6Nq8kVIwZA45rkqz0iBkDwviFCZyEGQJO1MYiSIAZAkmfveb8fBkC8wUET7h4GQMiMG1AeHgZA1EV7LFAdBkBr92qkgxwGQCsTALS4GwZA80tbV+8aBkCfcKiKJxoGQFVHHkphGQZAVmn+kZwYBkBaH5Ve2RcGQGs+OawXFwZASgVMd1cWBkBM+ji8mBUGQMfJdXfbFAZA2ySCpR8UBkDcoOdCZRMGQBiXOUysEgZAHQUVvvQRBkB5bSCVPhEGQOK4C86JEAZA2xeQZdYPBkC55G9YJA8GQCmGdqNzDgZAC1J4Q8QNBkDKcFI1Fg0GQBTB6nVpDAZA97svAr4LBkBiWRjXEwsGQBb1o/FqCgZA6DPaTsMJBkBq6crrHAkGQPr9jcV3CAZAHlVD2dMHBkBNtBIkMQcGQAaqK6OPBgZASHXFU+8FBkBf7R4zUAUGQPxpfj6yBAZAvKsxcxUEBkDoxI3OeQMGQJIC703fAgZACda47kUCBkCSvlWurQEGQHUzN4oWAQZAWI7Vf4AABkDi9a+M6/8FQK9ITK5X/wVAhAg34sT+BUDcRQMmM/4FQKiLSnei/QVAaMus0xL9BUB8SdA4hPwFQA=="},"shape":[500],"dtype":"float64","order":"little"}]]}}},"view":{"type":"object","name":"CDSView","id":"p1048","attributes":{"filter":{"type":"object","name":"AllIndices","id":"p1049"}}},"glyph":{"type":"object","name":"Line","id":"p1044","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.6,"line_width":3}},"nonselection_glyph":{"type":"object","name":"Line","id":"p1045","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.1,"line_width":3}},"muted_glyph":{"type":"object","name":"Line","id":"p1046","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.2,"line_width":3}}}}],"toolbar":{"type":"object","name":"Toolbar","id":"p1013","attributes":{"tools":[{"type":"object","name":"PanTool","id":"p1028"},{"type":"object","name":"WheelZoomTool","id":"p1029","attributes":{"renderers":"auto"}},{"type":"object","name":"BoxZoomTool","id":"p1030","attributes":{"overlay":{"type":"object","name":"BoxAnnotation","id":"p1031","attributes":{"syncable":false,"line_color":"black","line_alpha":1.0,"line_width":2,"line_dash":[4,4],"fill_color":"lightgrey","fill_alpha":0.5,"level":"overlay","visible":false,"left":{"type":"number","value":"nan"},"right":{"type":"number","value":"nan"},"top":{"type":"number","value":"nan"},"bottom":{"type":"number","value":"nan"},"left_units":"canvas","right_units":"canvas","top_units":"canvas","bottom_units":"canvas","handles":{"type":"object","name":"BoxInteractionHandles","id":"p1037","attributes":{"all":{"type":"object","name":"AreaVisuals","id":"p1036","attributes":{"fill_color":"white","hover_fill_color":"lightgray"}}}}}}}},{"type":"object","name":"SaveTool","id":"p1038"},{"type":"object","name":"ResetTool","id":"p1039"},{"type":"object","name":"HelpTool","id":"p1040"}]}},"left":[{"type":"object","name":"LinearAxis","id":"p1023","attributes":{"ticker":{"type":"object","name":"BasicTicker","id":"p1024","attributes":{"mantissas":[1,2,5]}},"formatter":{"type":"object","name":"BasicTickFormatter","id":"p1025"},"axis_label":"Normalized beta limit, \\  $$[\\beta_N]$$","major_label_policy":{"type":"object","name":"AllLabels","id":"p1026"}}}],"below":[{"type":"object","name":"LinearAxis","id":"p1018","attributes":{"ticker":{"type":"object","name":"BasicTicker","id":"p1019","attributes":{"mantissas":[1,2,5]}},"formatter":{"type":"object","name":"BasicTickFormatter","id":"p1020"},"axis_label":"Aspect ratio, \\  $$[A]$$","major_label_policy":{"type":"object","name":"AllLabels","id":"p1021"}}}],"center":[{"type":"object","name":"Grid","id":"p1022","attributes":{"axis":{"id":"p1018"}}},{"type":"object","name":"Grid","id":"p1027","attributes":{"dimension":1,"axis":{"id":"p1023"}}}]}}]}}
    </script>
    <script type="text/javascript">
      (function() {
        const fn = function() {
          Bokeh.safely(function() {
            (function(root) {
              function embed_document(root) {
              const docs_json = document.getElementById('ed3c7b85-c675-406a-a3ff-f1b0785e8266').textContent;
              const render_items = [{"docid":"62884a47-c32c-4b55-a144-365023ac52d1","roots":{"p1004":"cf29dd69-dc18-43d8-9476-146578601b91"},"root_ids":["p1004"]}];
              root.Bokeh.embed.embed_items(docs_json, render_items);
              }
              if (root.Bokeh !== undefined) {
                embed_document(root);
              } else {
                let attempts = 0;
                const timer = setInterval(function(root) {
                  if (root.Bokeh !== undefined) {
                    clearInterval(timer);
                    embed_document(root);
                  } else {
                    attempts++;
                    if (attempts > 100) {
                      clearInterval(timer);
                      console.log("Bokeh: ERROR: Unable to run BokehJS code because BokehJS library is missing");
                    }
                  }
                }, 10, root)
              }
            })(window);
          });
        };
        if (document.readyState != "loading") fn();
        else document.addEventListener("DOMContentLoaded", fn);
      })();
    </script>
  </body>
</html>

---------

#### Menard Beta Relation | `calculate_beta_norm_max_menard()`

This can be activated by stating `i_beta_norm_max = 3` in the input file.

`beta_norm_max` is set to `beta_norm_max_menard` using[^4]:


$$
\texttt{beta_norm_max_menard} = g =3.12+3.5\epsilon^{1.7}
$$

Found as a reasonable fit to the computed no wall limit at $f_{\text{BS}} \approx 50%$. Uses maximum $\kappa$ data from NSTX at $A = 1.45, A = 1.75.$ Along with record $\beta_{\text{T}}$ data from DIII-D at $A = 2.9$ and high $\kappa$.

**This is only recommended for spherical tokamaks**

**It is recommended to use this switch with [`i_ind_plasma_internal_norm = 2`](../plasma_current/plasma_inductance.md#menard-inductance-relation) as they are self-consistent with each other.**

<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <title>Menard Normalized Beta Limit</title>
    <style>
      html, body {
        box-sizing: border-box;
        display: flow-root;
        height: 100%;
        margin: 0;
        padding: 0;
      }
    </style>
    <script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-3.6.0.min.js"></script>
    <script type="text/javascript" src="https://cdn.bokeh.org/bokeh/release/bokeh-mathjax-3.6.0.min.js"></script>
    <script type="text/javascript">
        Bokeh.set_log_level("info");
    </script>
  </head>
  <body>
    <div id="eaf741e9-fe70-489d-96a0-bfc0abcdfa36" data-root-id="p1004" style="display: contents;"></div>
  
    <script type="application/json" id="c6bbe5cd-254d-4dc8-8c3c-700334e74f6d">
      {"84526480-8179-4d70-b025-efb598715583":{"version":"3.6.0","title":"Bokeh Application","roots":[{"type":"object","name":"Figure","id":"p1004","attributes":{"width":400,"height":400,"x_range":{"type":"object","name":"Range1d","id":"p1014","attributes":{"start":1,"end":5}},"y_range":{"type":"object","name":"Range1d","id":"p1015","attributes":{"start":2,"end":8}},"x_scale":{"type":"object","name":"LinearScale","id":"p1016"},"y_scale":{"type":"object","name":"LinearScale","id":"p1017"},"title":{"type":"object","name":"Title","id":"p1007","attributes":{"text":"Menard Normalized Beta Limit"}},"renderers":[{"type":"object","name":"GlyphRenderer","id":"p1047","attributes":{"data_source":{"type":"object","name":"ColumnDataSource","id":"p1001","attributes":{"selected":{"type":"object","name":"Selection","id":"p1002","attributes":{"indices":[],"line_indices":[]}},"selection_policy":{"type":"object","name":"UnionRenderers","id":"p1003"},"data":{"type":"map","entries":[["x",{"type":"ndarray","array":{"type":"bytes","data":"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"},"shape":[500],"dtype":"float64","order":"little"}],["y",{"type":"ndarray","array":{"type":"bytes","data":"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"},"shape":[500],"dtype":"float64","order":"little"}]]}}},"view":{"type":"object","name":"CDSView","id":"p1048","attributes":{"filter":{"type":"object","name":"AllIndices","id":"p1049"}}},"glyph":{"type":"object","name":"Line","id":"p1044","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.6,"line_width":3}},"nonselection_glyph":{"type":"object","name":"Line","id":"p1045","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.1,"line_width":3}},"muted_glyph":{"type":"object","name":"Line","id":"p1046","attributes":{"x":{"type":"field","field":"x"},"y":{"type":"field","field":"y"},"line_color":"#1f77b4","line_alpha":0.2,"line_width":3}}}}],"toolbar":{"type":"object","name":"Toolbar","id":"p1013","attributes":{"tools":[{"type":"object","name":"PanTool","id":"p1028"},{"type":"object","name":"WheelZoomTool","id":"p1029","attributes":{"renderers":"auto"}},{"type":"object","name":"BoxZoomTool","id":"p1030","attributes":{"overlay":{"type":"object","name":"BoxAnnotation","id":"p1031","attributes":{"syncable":false,"line_color":"black","line_alpha":1.0,"line_width":2,"line_dash":[4,4],"fill_color":"lightgrey","fill_alpha":0.5,"level":"overlay","visible":false,"left":{"type":"number","value":"nan"},"right":{"type":"number","value":"nan"},"top":{"type":"number","value":"nan"},"bottom":{"type":"number","value":"nan"},"left_units":"canvas","right_units":"canvas","top_units":"canvas","bottom_units":"canvas","handles":{"type":"object","name":"BoxInteractionHandles","id":"p1037","attributes":{"all":{"type":"object","name":"AreaVisuals","id":"p1036","attributes":{"fill_color":"white","hover_fill_color":"lightgray"}}}}}}}},{"type":"object","name":"SaveTool","id":"p1038"},{"type":"object","name":"ResetTool","id":"p1039"},{"type":"object","name":"HelpTool","id":"p1040"}]}},"left":[{"type":"object","name":"LinearAxis","id":"p1023","attributes":{"ticker":{"type":"object","name":"BasicTicker","id":"p1024","attributes":{"mantissas":[1,2,5]}},"formatter":{"type":"object","name":"BasicTickFormatter","id":"p1025"},"axis_label":"Normalized beta limit, \\  $$[\\beta_N]$$","major_label_policy":{"type":"object","name":"AllLabels","id":"p1026"}}}],"below":[{"type":"object","name":"LinearAxis","id":"p1018","attributes":{"ticker":{"type":"object","name":"BasicTicker","id":"p1019","attributes":{"mantissas":[1,2,5]}},"formatter":{"type":"object","name":"BasicTickFormatter","id":"p1020"},"axis_label":"Aspect ratio, \\  $$[A]$$","major_label_policy":{"type":"object","name":"AllLabels","id":"p1021"}}}],"center":[{"type":"object","name":"Grid","id":"p1022","attributes":{"axis":{"id":"p1018"}}},{"type":"object","name":"Grid","id":"p1027","attributes":{"dimension":1,"axis":{"id":"p1023"}}}]}}]}}
    </script>
    <script type="text/javascript">
      (function() {
        const fn = function() {
          Bokeh.safely(function() {
            (function(root) {
              function embed_document(root) {
              const docs_json = document.getElementById('c6bbe5cd-254d-4dc8-8c3c-700334e74f6d').textContent;
              const render_items = [{"docid":"84526480-8179-4d70-b025-efb598715583","roots":{"p1004":"eaf741e9-fe70-489d-96a0-bfc0abcdfa36"},"root_ids":["p1004"]}];
              root.Bokeh.embed.embed_items(docs_json, render_items);
              }
              if (root.Bokeh !== undefined) {
                embed_document(root);
              } else {
                let attempts = 0;
                const timer = setInterval(function(root) {
                  if (root.Bokeh !== undefined) {
                    clearInterval(timer);
                    embed_document(root);
                  } else {
                    attempts++;
                    if (attempts > 100) {
                      clearInterval(timer);
                      console.log("Bokeh: ERROR: Unable to run BokehJS code because BokehJS library is missing");
                    }
                  }
                }, 10, root)
              }
            })(window);
          });
        };
        if (document.readyState != "loading") fn();
        else document.addEventListener("DOMContentLoaded", fn);
      })();
    </script>
  </body>
</html>

---------

#### Tholerus Relation | `calculate_beta_norm_max_thloreus()`

This can be activated by stating `i_beta_norm_max = 4` in the input file.

`beta_norm_max` is set to `beta_norm_max_tholerus` using[^5]:


$$
C_{\beta}\approx\frac{(g-3.7)F_p}{12.5-3.5 F_p}
$$

where $F_p$ is the pressure peaking, $F_p = p_{\text{ax}} / \langle p \rangle$ and $C_{\beta}$ is the destabilization parameter (default 0.5)[^5].  

**This is only recommended for spherical tokamaks**

-------------

#### Stambaugh Relation | `calculate_beta_norm_max_stambaugh()`

This can be activated by stating `i_beta_norm_max = 5` in the input file.

`beta_norm_max` is set to `beta_norm_max_stambaugh` using[^8] [^9]:


$$
g=\frac{f_\beta \times 10 \times\left(-0.7748+1.2869 \kappa-0.2921 \kappa^2+0.0197 \kappa^3\right)}{A^{0.5523} \times \tanh \left[(1.8524+0.2319 \kappa) / A^{0.6163}\right]}
$$

This fit was done for $A = 1.2 -7.0, \kappa = 1.5-6.0$ with $\delta = 0.5$ for nearly 100% bootstrap current

---------

## Stored energy | `calculate_plasma_energy_from_beta()`

As the $\beta$ metric is simply the ratio between the kinetic pressure of the plasma and the magnetic pressure, $\beta$ can be used to get the total kinetic energy of the plasma (assuming the total $\beta$ is used).


$$
E_{\text{plasma}} \approx \frac{3}{2}\frac{\beta B^2}{2\mu_0}V_{\text{plasma}}
$$

---------------

## Key Constraints

### Beta consistency

This constraint can be activated by stating `icc = 1` in the input file.

Ensures the relationship between $\beta$, density, temperature and total magnetic field is withheld by checking the fixed input or iteration variable $\texttt{beta}$ is consistent in value with the rest of the physics parameters

$$
\texttt{beta_total_vol_avg} \equiv \frac{2\mu_0 \langle n_{\text{e}}T_{\text{e}}+n_{\text{i}}T_{\text{i}}\rangle}{B^2} + \beta_{\alpha} + \beta_{\text{beam}}
$$

Here the calculation of the volume averaged pressure of the ions and electrons has to use the density weighted temperature for each. This is because $\langle nT \rangle_{\text{V}} \neq \langle n \rangle_{\text{V}} \langle T \rangle_{\text{V}}$, where $\text{V}$ denotes the volume averaged value. The true value is, $\langle nT \rangle_{\text{V}} = \langle n \rangle_{\text{V}} \langle T \rangle_{\text{n}}$, where $\text{n}$ is the density weighted averaged. For example:

$$
\langle n_{\text{e}}T_{\text{e}} \rangle_{\text{V}} = \overbrace{\langle n_{\text{e}} \rangle_{\text{V}}}^{\texttt{nd_plasma_electrons_vol_avg}} \times  \overbrace{\langle T_{\text{e}} \rangle_{\text{n}}}^{\texttt{temp_plasma_electron_density_weighted_kev}}
$$

**It is highly recommended to always have this constraint on as it is a global consistency checker**

----------------

### Poloidal beta and inverse aspect upper limit

This constraint can be activated by stating `icc = 6` in the input file [^6].

The limiting value of $\epsilon\beta_p$ is be set using input parameter `beta_poloidal_eps_max`.

!!! note "Origin of the $\epsilon\beta_p$ limit"

    High poloidal beta shots in TFTR were performed[^6] and it was found that as $\beta_p$,
    exceeds approximately 1.2 times the aspect ratio, a separatrix with
    an inside poloidal field null is observed to limit the outer boundary
    of the plasma. Since the curvature of TFTR’s applied vertical field
    is constant, the appearance of the poloidal field null corresponds to
    the equilibrium poloidal beta limit.

--------------------

### Beta upper limit

This constraint can be activated by stating `icc = 24` in the input file.

It is the general setting of the $\beta$ limit depending on the $\beta_{\text{N}}$ value calculated in the [beta limit](#beta-limit) calculations.

The upper limit value of beta is calculated by `calculate_beta_limit()`. The beta
coefficient $g$ can be set using `beta_norm_max`, depending on the setting of [`i_beta_norm_max`](#setting-the-beta--coefficient). It can be set directly or follow some relation.

**It is recommended to have this constraint on as it is a plasma stability model**

--------------------

### Poloidal upper limit

This constraint can be activated by stating `icc = 48` in the input file.

The value of `beta_poloidal_max` can be set to the desired maximum poloidal beta.

-------------------

### Beta lower limit

This constraint can be activated by stating `icc = 84` in the input file.

The value of `beta_vol_avg_min` can be set to the desired minimum total beta.

[^0]: F. Troyon et.al,  “Beta limit in tokamaks. Experimental and computational status,” Plasma Physics and Controlled Fusion, vol. 30, no. 11, pp. 1597–1609, Oct. 1988, doi: https://doi.org/10.1088/0741-3335/30/11/019.

[^1]: N.A. Uckan and ITER Physics Group, 'ITER Physics Design Guidelines: 1989',

[^2]: D.J. Ward, 'PROCESS Fast Alpha Pressure', Work File Note F/PL/PJK/PROCESS/CODE/050

[^3]: Tokamaks 4th Edition, Wesson, page 116

[^4]: J. E. Menard et al., “Fusion nuclear science facilities and pilot plants based on the spherical tokamak,” Nuclear Fusion, vol. 56, no. 10, p. 106023, Aug. 2016, doi: https://doi.org/10.1088/0029-5515/56/10/106023.

[^5]: E. Tholerus et al., “Flat-top plasma operational space of the STEP power plant,” Nuclear Fusion, Aug. 2024, doi: https://doi.org/10.1088/1741-4326/ad6ea2.
‌
[^6]: M. E. Mauel et al., “Operation at the tokamak equilibrium poloidal beta-limit in TFTR,” Nuclear Fusion, vol. 32, no. 8, pp. 1468–1473, Aug. 1992. doi:https://dx.doi.org/10.1088/0029-5515/32/8/I14

[^7]: T. T. S et al., “Profile Optimization and High Beta Discharges and Stability of High Elongation Plasmas in the DIII-D Tokamak,” Osti.gov, Oct. 1990. https://www.osti.gov/biblio/6194284 (accessed Dec. 19, 2024).

[^8]: R. D. Stambaugh et al., “Fusion Nuclear Science Facility Candidates,” Fusion Science and Technology, vol. 59, no. 2, pp. 279–307, Feb. 2011, doi: https://doi.org/10.13182/fst59-279.

[^9]: Y. R. Lin-Liu and R. D. Stambaugh, “Optimum equilibria for high performance, steady state tokamaks,” Nuclear Fusion, vol. 44, no. 4, pp. 548–554, Mar. 2004, doi: https://doi.org/10.1088/0029-5515/44/4/009.
‌
‌
‌
