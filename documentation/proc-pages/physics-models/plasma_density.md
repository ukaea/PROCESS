Several density limit models[^1] are available in PROCESS. These are
calculated in routine `culdlm`, which is called by `physics`. To enforce any of 
these limits, turn on constraint equation no. 5 with iteration variable no. 9 
(`fdene`). In addition, switch `idensl` must be set to the relevant value, as 
follows:

| `idensl` | Description |
| :-: | - |
| 1 | ASDEX model |
| 2 | Borrass model for ITER, I |
| 3 | Borrass model for ITER, II |
| 4 | JET edge radiation model |
| 5 | JET simplified model |
| 6 | Hugill-Murakami $M.q$ model |
| 7 | Greenwald model: $n_G=10^{14} \frac{I_p}{\pi a^2}$ where the units are m and ampere. For the Greenwald model the limit applies to the line-averaged electron density, not the volume-averaged density. |



[^1]: T. C. Hender et al., 'Physics Assessment for the European Reactor Study',
AEA Fusion Report AEA FUS 172 (1992)