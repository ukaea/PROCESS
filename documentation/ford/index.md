project: PROCESS
src_dir: ../../source
output_dir: ../../ford_site
summary: PROCESS is a systems code at CCFE that calculates in a 
         self-consistent manner the parameters of a fusion power
		 plant with a specified performance, ensuring that 
		 its operating limits are not violated, and with the
		 option to optimise a given function of these parameters.
display: public
         protected
         private
fpp_extensions: f90
graph: true
source: true
search: false
author: UKAEA

During the course of studies into a proposed fusion power plant, there
may be times when questions of the following type arise:

- Are the machine’s physics and engineering parameters consistent with
  one another?
- Which machine of a given size and shape produces the cheapest
  electricity?
- What is the effect of a more optimistic limit on the maximum plasma
  density on the amount of auxiliary power required?

Questions such as these are extremely difficult to answer, since the
large number of parameters involved are highly dependent on one
another. Fortunately, computer programs have been written to address
these issues, and PROCESS is one of them.

Suppose that an outline power plant design calls for a machine with a
given size and shape, which will produce a certain net electric
power. There may be a vast number of different conceptual machines
that satisfy the problem as stated so far, and PROCESS can be used in
“non-optimisation” mode to find one of these whose physics and
engineering parameters are self-consistent. However, the machine found
by PROCESS in this manner may not be possible to build in practice —
the coils may be overstressed, for instance, or the plasma pressure
may exceed the maximum possible value. PROCESS contains a large number
of constraints to prevent the code from finding a machine with such
problems, and running the code in so-called “optimisation” mode forces
these constraints to be met. The number of possible conceptual
machines is thus considerably reduced, and optimisation of the
parameters with respect to (say) the cost of electricity will reduce
this number to a minimum (possibly one).

Formally then, PROCESS is a systems code that calculates in a
self-consistent manner the parameters of a fusion power plant with a
specified performance, ensuring that its operating limits are not
violated, and with the option to optimise a given function of these
parameters.
