# Shield

The shield reduces the neutron flux reaching the vacuum vessel, TF coils and 
beyond. This minimises the radiological impact of the neutrons, and their 
heating of the TF coils which, if superconducting, need to remain at liquid 
helium temperatures. As with the blanket the energy deposited in the coolant 
may or may not be used to produce electricity, depending on the value of the 
switch `i_shld_primary_heat`. The shield coolant fraction by volume is set via the input 
parameter `vfshld`.

The inboard and outboard shield thicknesses (`dr_shld_inboard` and `dr_shld_outboard`, 
respectively) may be used as iteration variables.
