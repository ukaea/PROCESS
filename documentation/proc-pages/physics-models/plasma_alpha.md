# Fast Alpha Pressure Contribution

The pressure contribution from the fast alpha particles can be controlled using switch `ifalphap`. 
There are two options 1[^1] and 2[^2]:

$$\begin{aligned}
\frac{\beta_{\alpha}}{\beta_{th}} & = 0.29 \, \left( \langle T_{10} \rangle -
  0.37 \right) \, \left( \frac{n_{DT}}{n_e} \right)^2
\hspace{20mm} \mbox{if alphap = 0} \\
\frac{\beta_{\alpha}}{\beta_{th}} & = 0.26 \, \left( \langle T_{10} \rangle -
  0.65 \right)^{0.5} \, \left( \frac{n_{DT}}{n_e} \right)^2
\hspace{16mm} \mbox{if alphap = 1 (default)}
\end{aligned}$$

The latter model is a better estimate at higher temperatures.

[^1]: T. C. Hender et al., 'Physics Assessment for the European Reactor Study',

[^2]: H. Lux, R. Kemp, D.J. Ward, M. Sertoli, 'Impurity radiation in DEMO 
systems modelling', Fus. Eng.  | Des. **101**, 42-51 (2015)