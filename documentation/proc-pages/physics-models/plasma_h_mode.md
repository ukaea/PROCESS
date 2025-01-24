## L-H Power Threshold Scalings

Transitions from a standard confinement mode (L-mode) to an improved
confinement regime (H-mode), called L-H transitions, are observed in most
tokamaks. A range of scaling laws are available that provide estimates of the
heating power required to initiate these transitions, via extrapolations
from present-day devices. PROCESS calculates these power threshold values
for the scaling laws listed in the table below, in routine `pthresh`.

For an H-mode plasma, use input parameter `ilhthresh` to
select the scaling to use, and turn on constraint equation no. 15 with
iteration variable no. 103 (`flhthresh`). By default, this will ensure
that the power reaching the divertor is at least equal to the threshold power
calculated for the chosen scaling, which is a necessary condition for
H-mode. 

For an L-mode plasma, use input parameter `ilhthresh` to
select the scaling to use, and turn on constraint equation no. 15 with 
iteration variable no. 103 (`flhthresh`). Set lower and upper bounds for 
the f-value `boundl(103) = 0.001` and `boundu(103) = 1.0` 
to ensure that the power does not exceed the calculated threshold, 
and therefore the machine remains in L-mode.


| `ilhthresh` | Name | Reference |
| :-: | - | - |
| 1 | ITER 1996 nominal | ITER Physics Design Description Document |
| 2 | ITER 1996 upper bound | D. Boucher, p.2-2 |
| 3 | ITER 1996 lower bound | 
| 4 | ITER 1997 excluding elongation | J. A. Snipes, ITER H-mode Threshold Database |
| 5 | ITER 1997 including elongation |  Working Group, Controlled Fusion and Plasma Physics, 24th EPS conference, Berchtesgaden, June 1997, vol.21A, part III, p.961 |
| 6 | Martin 2008 nominal | Martin et al, 11th IAEA Tech. Meeting |
| 7 | Martin 2008 95% upper bound |  H-mode Physics and Transport Barriers, Journal |
| 8 | Martin 2008 95% lower bound |  of Physics: Conference Series **123**, 2008 |
| 9 | Snipes 2000 nominal | J. A. Snipes and the International H-mode |
| 10| Snipes 2000 upper bound | Threshold Database Working Group |
| 11| Snipes 2000 lower bound |  2000, Plasma Phys. Control. Fusion, 42, A299 |
| 12| Snipes 2000 (closed divertor): nominal | 
| 13| Snipes 2000 (closed divertor): upper bound | 
| 14| Snipes 2000 (closed divertor): lower bound | 
| 15| Hubbard 2012 L-I threshold scaling: nominal | [Hubbard et al. (2012; Nucl. Fusion 52 114009)](https://iopscience.iop.org/article/10.1088/0029-5515/52/11/114009) |
| 16| Hubbard 2012 L-I threshold scaling: lower bound | [Hubbard et al. (2012; Nucl. Fusion 52 114009)](https://iopscience.iop.org/article/10.1088/0029-5515/52/11/114009 |
| 17| Hubbard 2012 L-I threshold scaling: upper bound | [Hubbard et al. (2012; Nucl. Fusion 52 114009)](https://iopscience.iop.org/article/10.1088/0029-5515/52/11/114009 |
| 18| Hubbard 2017 L-I threshold scaling | [Hubbard et al. (2017; Nucl. Fusion 57 126039)](https://iopscience.iop.org/article/10.1088/1741-4326/aa8570) |
| 19 | Martin 2008 aspect ratio corrected nominal | Martin et al (2008; J Phys Conf, 123, 012033) |
| 20 | Martin 2008 aspect ratio corrected 95% upper bound | [Takizuka et al. (2004; Plasma Phys. Contol. Fusion, 46, A227)](https://iopscience.iop.org/article/10.1088/0741-3335/46/5A/024)  |
| 21 | Martin 2008 aspect ratio corrected 95% lower bound |  