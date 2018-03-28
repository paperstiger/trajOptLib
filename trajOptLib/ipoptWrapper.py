#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
ipoptWrapper.py

A wrapper for ipopt solver
"""


import numpy as np
try:
    import pyipopt
except:
    print("pyipopt not detected, do not use ipopt")
from . import trajOptProblem
from .libsnopt import result

class ipOption(object):
    """A class containing settings for ipopt solver."""
    def __init__(self):
        self.intOpt = []
        self.strOpt = []
        self.floatOpt = []
        self.tol = 1e-6
        self.max_iter = 1000
        self.dual_inf_tol = 10
        self.constr_vio_tol = 1e-6

    def addIntOption(self, key, value):
        self.intOpt.append((key, value))

    def addFloatOption(self, key, value):
        self.floatOpt.append((key, value))

    def addStrOption(self, key, value):
        self.strOpt.append((key, value))

    def addOption(self, key, value):
        if isinstance(value, int):
            self.addIntOption(key, value)
        elif isinstance(value, float):
            self.addFloatOption(key, value)
        elif isinstance(value, str):
            self.addStrOption(key, value)

class ipSolver(object):
    """A solver class for ipopt. It accepts my conventional problem type."""
    def __init__(self, prob, option=None):
        assert isinstance(prob, trajOptProblem)
        nvar = prob.numSol
        x_L = prob.xlb
        x_U = prob.xub
        ncon = prob.numF
        g_L = prob.lb
        g_U = prob.ub
        nnzj = prob.nG
        nnzh = 1  # this is not reasonable
        self.prob = prob
        self.nlp = pyipopt.create(nvar, x_L, x_U, ncon, g_L, g_U, nnzj, nnzh, prob.ipEvalF,
                prob.ipEvalGradF, prob.ipEvalG, prob.ipEvalJacG)
        if option is None:
            option = ipOption()
        self._parseOption(option)

    def _parseOption(self, option):
        """Parse an option object."""
        self.nlp.num_option('tol', option.tol)
        self.nlp.int_option('max_iter', option.max_iter)
        self.nlp.num_option('dual_inf_tol', option.dual_inf_tol)
        self.nlp.num_option('constr_viol_tol', option.constr_vio_tol)
        for key, value in option.intOpt:
            self.nlp.int_option(key, value)
        for key, value in option.floatOpt:
            self.nlp.num_option(key, value)
        for key, value in option.strOpt:
            self.nlp.str_option(key, value)

    def solveRand(self):
        """Use problem to generate a random guess and solve it from there.

        It generates a random initial guess (which might be far off for unknown problem) and use it to solve.

        """
        x0 = self.prob.randomGenX()
        return self.solveGuess(x0)

    def solveGuess(self, x0):
        """Solve the problem with some user-supplied guess.

        :param x0: ndarray, solution/guess to the problem.

        """
        x, zl, zu, lmd, obj, status = self.nlp.solve(x0)
        self.nlp.close()
        rst = result()
        rst.flag = status + 1  # such that 1 means successful
        rst.sol = x
        rst.obj = obj
        return rst
