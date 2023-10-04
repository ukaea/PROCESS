from process.caller import Caller
from process.fortran import global_variables as gv
from process.fortran import constraints
from process.fortran import cost_variables as cv
from process.fortran import numerics
from process.fortran import physics_variables as pv
from process.fortran import stellarator_variables as sv
from process.fortran import times_variables as tv
from process.fortran import function_evaluator
import numpy as np
import math
import logging

logger = logging.getLogger(__name__)


class Evaluators:
    """Calls models to evaluate function and gradient functions."""

    def __init__(self, models, x):
        """Instantiate Caller with model objects.

        :param models: physics and engineering model objects
        :type models: process.main.Models
        :param x: optimisation parameters
        :type x: np.ndarray
        """
        self.caller = Caller(models, x)

    def fcnvmc1(self, n, m, xv, ifail):
        """Function evaluator for VMCON.

        This routine is the function evaluator for the VMCON
        maximisation/minimisation routine.

        It calculates the objective and constraint functions at the
        n-dimensional point of interest xv.
        Note that the equality constraints must precede the inequality
        constraints in conf.
        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        :param n: number of variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param xv: scaled variable values, length n
        :type xv: numpy.array
        :param ifail: ifail error flag
        :type ifail: int
        :return: tuple containing: objfn objective function, conf(m) constraint
        functions
        :rtype: tuple
        """
        # Output array for constraint functions
        conf = np.zeros(m, dtype=np.float64, order="F")

        # Evaluate machine parameters at xv
        self.caller.call_models(xv)

        # Convergence loop to ensure burn time consistency
        if sv.istell == 0:
            loop = 0
            while (loop < 10) and (
                abs((tv.tburn - tv.tburn0) / max(tv.tburn, 0.01)) > 0.001
            ):
                loop += 1
                self.caller.call_models(xv)
                if gv.verbose == 1:
                    print("Internal tburn consistency check: ", tv.tburn, tv.tburn0)

            if loop >= 10:
                print(
                    "Burn time values are not consistent in iteration: ",
                    numerics.nviter,
                )
                print("tburn, tburn0: ", tv.tburn, tv.tburn0)

        # Evaluate figure of merit (objective function)
        objf = function_evaluator.funfom()

        # Evaluate constraint equations
        conf, _, _, _, _ = constraints.constraint_eqns(m, -1)

        # Verbose diagnostics
        if gv.verbose == 1:
            summ = 0.0
            for i in range(m):
                summ = summ + conf[i] ** 2

            sqsumconfsq = math.sqrt(summ)
            logger.debug("Key evaluator values:")
            logger.debug(f"{numerics.nviter = }")
            logger.debug(f"{(1 - (ifail % 7)) - 1 = }")
            logger.debug(f"{(numerics.nviter % 2) - 1 = }")
            logger.debug(f"{pv.te = }")
            logger.debug(f"{cv.coe = }")
            logger.debug(f"{pv.rmajor = }")
            logger.debug(f"{pv.powfmw = }")
            logger.debug(f"{pv.bt = }")
            logger.debug(f"{tv.tburn = }")
            logger.debug(f"{sqsumconfsq = }")
            logger.debug(f"{xv = }")

        return objf, conf

    def fcnvmc2(self, n, m, xv, lcnorm):
        """Gradient function evaluator for VMCON.

        This routine is the gradient function evaluator for the VMCON
        maximisation/minimisation routine. It calculates the gradients of the
        objective and constraint functions at the n-dimensional point of interest
        xv. Note that the equality constraints must precede the inequality
        constraints in conf. The constraint gradients or normals are returned as the
        columns of cnorm.

        AEA FUS 251: A User's Guide to the PROCESS Systems Code
        :param n: number of variables
        :type n: int
        :param m: number of constraints
        :type m: int
        :param xv: scaled variable names, size n
        :type xv: numpy.array
        :param lcnorm: number of columns in cnorm
        :type lcnorm: int
        :return: fgrdm (numpy.array (n)) gradient of the objective function
        cnorm (numpy.array (lcnorm, m)) constraint gradients, i.e. cnorm[i, j] is
        the derivative of constraint j w.r.t. variable i
        :rtype: tuple
        """
        xfor = np.zeros(n, dtype=np.float64, order="F")
        xbac = np.zeros(n, dtype=np.float64, order="F")
        cfor = np.zeros(m, dtype=np.float64, order="F")
        cbac = np.zeros(m, dtype=np.float64, order="F")
        fgrd = np.zeros(n, dtype=np.float64, order="F")
        cnorm = np.zeros((lcnorm, m), dtype=np.float64, order="F")

        ffor = 0.0
        fbac = 0.0

        for i in range(n):
            for j in range(n):
                xfor[j] = xv[j]
                xbac[j] = xv[j]
                if i == j:
                    xfor[i] = xv[j] * (1.0 + numerics.epsfcn)
                    xbac[i] = xv[j] * (1.0 - numerics.epsfcn)

            # Evaluate at (x+dx)
            self.caller.call_models(xfor)
            ffor = function_evaluator.funfom()
            cfor, _, _, _, _ = constraints.constraint_eqns(m, -1)

            # Evaluate at (x-dx)
            self.caller.call_models(xbac)
            fbac = function_evaluator.funfom()
            cbac, _, _, _, _ = constraints.constraint_eqns(m, -1)

            # Calculate finite difference gradients
            fgrd[i] = (ffor - fbac) / (xfor[i] - xbac[i])

            for j in range(m):
                cnorm[i, j] = (cfor[j] - cbac[j]) / (xfor[i] - xbac[i])

        # Additional evaluation call to ensure that final result is consistent
        # with the correct iteration variable values.
        # If this is not done, the value of the nth (i.e. final) iteration
        # variable in the solution vector is inconsistent with its value
        # shown elsewhere in the output file, which is a factor (1-epsfcn)
        # smaller (i.e. its xbac value above).
        self.caller.call_models(xv)

        return fgrd, cnorm
