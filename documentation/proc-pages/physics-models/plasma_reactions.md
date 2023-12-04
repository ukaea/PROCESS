# Fusion Reactions

The most likely fusion reaction to be utilised in a power plant is the
deuterium-tritium reaction:

$$
\mathrm{D + T} \Longrightarrow \mathrm{^{4}He + n + 17.6 \,MeV}
$$

20% of the energy produced is given to the alpha particles (\(^4\)He). The remaining 80% is carried 
away by the neutrons, which deposit their energy within the blanket and shield and other reactor components. 
The fraction of the alpha energy deposited in the plasma is `falpha`.

PROCESS can also model D-\(^3\)He power plants, which utilise the following 
primary fusion reaction:

$$
\mathrm{D + \text{$^3$He}} \Longrightarrow \mathrm{^{4}He + p + 18.3 \,MeV}
$$

The fusion reaction rate is significantly different to that for D-T fusion,
and the power flow from the plasma is modified since charged particles are
produced rather than neutrons. Because only charged particles (which remain in
the plasma) are produced by this reaction, the whole of the fusion power is
used to heat the plasma. Useful energy is extracted from the plasma since the
radiation power produced is very high, and this, in theory, can be converted to
electricity without using a thermal cycle.

Since the temperature required to ignite the D-\(^3\)He reaction is considerably
higher than that for D-T, it is necessary to take into account the following
D-D reactions, which have significant reaction rates at such temperatures:

$$\begin{aligned}
\mathrm{D + D}  & \Longrightarrow \mathrm{^{3}He + n + 3.27 \,MeV} \\
\mathrm{D + D}  & \Longrightarrow \mathrm{T + p + 4.03 \,MeV}
\end{aligned}$$

Also, as tritium is produced by the latter reaction, D-T fusion also occurs. 
As a result, there is still a small amount of neutron power
extracted from the plasma.

Pure D-\(^3\)He tokamak power plants do not include breeding blankets, because 
no tritium needs to be produced for fuel.

The contributions from all four of the above fusion reactions are included in
the total fusion power production calculation. The fusion reaction rates are
calculated using the parameterizations in [^4], integrated over the plasma 
profiles (correctly, with or without pedestals).

The fractional composition of the 'fuel' ions (D, T and \(^3\)He) is
controlled using the three variables `fdeut`, `ftrit` and `fhe3`, respectively:

$$\begin{aligned}
n_{\mbox{fuel}}  & = n_D + n_T + n_{\mathrm{^{3}He}}  \;\;\; \mbox{particles/m$^3$} \\
n_D  & = \mathtt{fdeut} \, n_{\mbox{fuel}} \\
n_T  & = \mathtt{ftrit} \, n_{\mbox{fuel}} \\
n_{\mathrm{^{3}He}} & = \mathtt{fhe3} \, n_{\mbox{fuel}}
\end{aligned}$$

PROCESS checks that $fdeut + ftrit + fhe3 = 1.0$, and stops with an error message otherwise.

Constraint equation no. 28 can be turned on to enforce the fusion gain *Q* to be at 
least equal to `bigqmin`.