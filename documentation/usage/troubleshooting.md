# Troubleshooting

Experience has shown that the first few attempts at running PROCESS with a new file tends to 
produce infeasible results; that is, the code will not find a consistent set of machine parameters. 
The highly non-linear nature of the numerics of PROCESS is the reason for this difficulty, and 
it often requires a great deal of painstaking adjustment of the input file to overcome.

## Error handling

In general, errors detected during a run are handled in a consistent manner, with the code 
producing useful diagnostic messages to help the user understand what has happened.

In the case of an unrecoverable error, PROCESS will fail with a Python exception that will detail 
the nature of the error and its location (line in the code).

In addition, PROCESS will log recoverable errors and warnings (and other information) to the process.log 
file and, for serious logs, the terminal. Logs (errors and warnings) that occur during the output phase are captured
and written to the OUTFile and terminal in the run summary. The ouput phase is when the models are run to write the 
results into the MFile, and occurs at the final iteration of an optimising run or when running an evaluation. These logs
can be classified into two types:
1. **Warning**: provides information that the user should be aware about but that is not an error. For example, if the users
input is modified by the code or if an extrapolation is used. 
2. **Errors**: provides information about an error that PROCESS is able to recover from, but which could impact the validity 
of the results. For example, an input to a model is unphysical causing an unreliable output. 

### General problems

A code of the size and complexity of PROCESS contains myriads of equations and variables. Virtually 
everything depends indirectly on everything else because of the nature of the code structure, so perhaps 
it is not surprising that it is often difficult to achieve a successful outcome.

Naturally, problems will occur if some of the parameters become un-physical. For example, if the aspect 
ratio becomes less than or equal to one, then we must expect problems to appear. For this reason, 
the bounds on the iteration variables should be selected with care.

Occasionally arithmetic ("NaN") errors are reported. They usually occur when the code is exploring 
un-physical values of the parameters, and often suggest that no feasible solution exists for the 
input file used.

The error messages produced by the code attempt to provide diagnostic information, telling the user 
where the problems occurs, and also suggest a possible solution. These messages are out of 
necessity brief, and so cannot promise to lead to a more successful outcome.

The is the option to turn on extra debugging output; to do this, set `verbose = 1` in the input file.

### Optimisation problems

On reflection it is perhaps surprising that PROCESS ever does manage to find the global minimum 
figure of merit value, if there are `nvar` iteration variables active the search is 
over `nvar`-dimensional parameter space, in which there may be many shallow minima of approximately 
equal depth. Remember that `nvar` is usually of the order of twenty.

The machine found by PROCESS may not, therefore, be the absolute optimal device. It is quite easy 
to have two or more solutions, with results only a few percent different, but a long way apart in 
parameter space. The technique of "stationary" scans is sometimes used in this situation: a scan is 
requested, but the same value of the scan variable is listed repeatedly.

Scans should be started in the middle of a range of values, to try to keep the scan within the same 
family of machines. The optimum machine found may otherwise suddenly jump to a new region of 
parameter space, causing the output variables to seem to vary unpredictably with the scanning variable.

It should be noted that in general the machine produced by PROCESS will always sit against one or 
more operation limits. If, during a scan, the limit being leant upon changes (i.e. if the machine 
jumps from leaning on the beta limit to leaning on the density limit) the output parameters may well 
become discontinuous in gradient, and trends may suddenly change direction.

### Unfeasible results

In the numerics section of the output file, the code indicates whether the run produced a feasible 
or unfeasible result.

The former implies a successful outcome, although it is always worth checking that the sum of the 
squares of the constraint residuals (`sqsumsq`) is small ($~10^{-3}$ or less); the code will issue 
a warning if the solver reports convergence but the value of `sqsumsq` exceeds $10^{-2}$. If this 
occurs, reducing the value of the HYBRD tolerance `ftol` or `VMCON` tolerance `epsvmc` as appropriate 
should indicate whether the result is valid ot not; the output can usually be trusted of (1) the 
constraint residuals[^1] fall as the tolerance is reduced to about $10^{-8}$, and (2) the code 
indicates that a feasible solution is still found.

An unfeasible result occurs if PROCESS cannot find a set of values for the iteration variables 
which satisfies all the given constraints. In this case, the values of the constraint residues 
shown in the output give some indication of which constraint equations are not being 
satisfied - those with the highest residues should be examined further. In optimisation mode, 
the code also indicates which iteration variables lie at the edge of their allowed range.

Unfeasible runs can be caused by specifying physical incompatible input parameters, using 
insufficient iteration variables, or by starting the problem with unsuitable values of the iteration variables.

The utility `run_process` carries out many runs, changing the starting values of the iteration 
variables randomly. It stops once a feasible solution is found.

It is important to choose the right number of *useful* iteration variables for the problem to be 
solved - it is possible to activate too many iteration variables as well as too few, some of which may be redundant.

Both optimisation and non-optimisation runs can fail with an error message suggesting that the 
iteration process is not making good progress. This is likely to be due to the code itself unable 
to escape a region of the parameters space where the minimum in the residuals is significantly 
above zero. In this situation, there is either no solution possible (the residuals can therefore 
never approach zero), or the topology of the local minimum makes it difficult for the code to 
escape to the global minimum. Again, a helpful technique os to wither change the list of iteration 
variables in use, or to simply modify their initial values to try to help the code avoid such regions.

A technique that occasionally removes problems due to unfeasible results, particularly if an error 
code `ifail = 3` is encountered during an optimisation run, is to adjust slightly one of the limits 
imposed on the iteration variables, even if the limit in question has not been reached. This subtly 
alters the gradients computed by the code during the iteration process and may tip the balance so 
that the code decides that the device produced is feasible after all. For instance, a certain 
component's temperature might be 400 K, and its maximum allowable temperature is 1000 K. Adjusting 
this limit to 900 K (which will make no different to the *actual* temperature) may be enough to 
persuade the code that it has found a feasible solution.

Similarly, the order in which the constraint equations and iteration variables are stored in the 
`icc` and `ixc` arrays can make the difference between a feasible and unfeasible result. This 
seemingly illogical behaviour is  typical of the way in which the code works.

Another technique in such situations may be to change the finite difference step length `epsfcn`, 
as this might subtly change the path taken in the approach towards a solution.

It may be the case that the act of satisfying all the required constraints is impossible. No 
machine can exist if the allowed operating regime is too restrictive, or if two constraints 
require conflicting (non-overlapping) parameters spaces. In this case some relaxation of the 
requirements is needed for the code to produce a successful machine design.

### Hints

The above sections should indicate that it is the complex inter-play between the constraint 
equations and the iteration variables that determines whether the code will eb successful at 
producing a useful result. It can be somewhat laborious process to arrive at a working vase, 
and experience is often of great value in this situation.

The lower and upper bounds of the iteration variables are all available to be changed in the input 
file. Constraints can be relaxed in a controlled manner by moving these bounds, although in some 
cases care should be taken to ensure that un-physical values cannot occur. The code indicates which 
iteration variables lie at the edge of their range.

It is suggested that constraint equations should be added one at a time, with required new iteration 
variables activated at each step. If the situation becomes unfeasible it can be helpful to reset 
the initial iteration variable values to those shown in the output from a previous feasible case 
and rerun the code.

[^1]: The constraint residuals are the final values of $c_i$ in the constraint equations. The value `sqsumsq` is the square root of the sum of the squares of these residuals.
