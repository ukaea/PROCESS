import logging
import math

import numpy as np

from process.caller import Caller
from process.data_structure import cost_variables as cv
from process.data_structure import global_variables as gv
from process.data_structure import numerics
from process.data_structure import physics_variables as pv
from process.data_structure import times_variables as tv

logger = logging.getLogger(__name__)


class Evaluators:
    """Calls models to evaluate function and gradient functions."""

    def __init__(self, models, _x):
        """Instantiate Caller with model objects.

        :param models: physics and engineering model objects
        :type models: process.main.Models
        :param x: optimisation parameters
        :type x: np.ndarray
        """
        self.caller = Caller(models)

    def fcnvmc1(self, _n, m, xv, ifail):
        """Function evaluator for VMCON.

        This routine is the function evaluator for the VMCON
        maximisation/minimisation routine.

        It calculates the objective and constraint functions at the
        n-dimensional point of interest xv.
        Note that the equality constraints must precede the inequality
        constraints in conf.
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
        objf, conf = self.caller.call_models(xv, m)

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
            logger.debug(f"{pv.temp_plasma_electron_vol_avg_kev = }")
            logger.debug(f"{cv.coe = }")
            logger.debug(f"{pv.rmajor = }")
            logger.debug(f"{pv.p_fusion_total_mw = }")
            logger.debug(f"{pv.b_plasma_toroidal_on_axis = }")
            logger.debug(f"{tv.t_plant_pulse_burn = }")
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
            ffor, cfor = self.caller.call_models(xfor, m)

            # Evaluate at (x-dx)
            fbac, cbac = self.caller.call_models(xbac, m)

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
        self.caller.call_models(xv, m)

        return fgrd, cnorm
