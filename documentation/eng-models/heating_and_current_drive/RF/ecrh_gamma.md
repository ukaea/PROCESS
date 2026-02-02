# ECRH User Input Gamma Model

- `i_hcd_primary/i_hcd_secondary` = 10 

This model allows the user to input a scaling factor to the current drive efficiency with the variable `eta_cd_norm_ecrh`. The value of this variable should follow the value and form of the expression below:

$$
\gamma_{CD} = \frac{\langle n_{e,20} \rangle I_{CD}R_0}{P_{CD}}
$$

The current drive efficiency is then calculated with $\gamma_{CD}$ in the form below:

$$
\text{Current drive efficiency [A/W]} =  \frac{\gamma_{CD}}{R_{0} n_{e,20}}
$$                