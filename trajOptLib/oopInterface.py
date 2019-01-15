#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
oopInterface.py

Construct a object oriented interface for problem solving.
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
from .trajOptCollocationProblem import trajOptCollocProblem
from .trajOptBase import lqrObj
from . import solver, snoptConfig
from . import ipSolver, ipOption
from .trajOptBase import nonLinearPointConstr
from .utility import showSol, getInf


class AbstractSolver(object):
    """A class that can provides abstract methods for a general solver."""
    def __init__(self, system, N, tf, ip=False, config=None):
        """Constructor for the problem.

        Parameters
        ----------
        system : a dae system
        N : knot point size
        tf : float/array-like, final time
        ip : bool, if use interior point solver
        config : a configuration object for either snopt or ipopt

        """
        self.N = N
        self.system = system
        self.tf = tf
        self.ip = ip
        self.prob = trajOptCollocProblem(system, N, 0, tf)
        self.solver_is_up_to_data = False
        self.config = config

    def constructSolver(self, ip=False):
        if ip:
            if self.config is None:
                option = ipOption()
                option.max_iter = 1000
                option.tol = 1e-4
                option.dual_inf_tol = 1e-4
                option.print_level = 0
            else:
                option = self.config
            self.option = option
            self.slv = ipSolver(self.prob, option)
        else:
            if self.config is None:
                cfg = snoptConfig()
                cfg.printLevel = 0
                cfg.verifyLevel = 0
                cfg.majorIterLimit = 500
            else:
                cfg = self.config
            self.option = cfg
            self.slv = solver(self.prob, cfg)

    def addLQRObj(self, *args, **kwargs):
        """Add a lqr objective."""
        obj = lqrObj(*args, **kwargs)
        self.prob.addLQRObj(obj)

    def preProcess(self, *args, **kwargs):
        self.prob.preProcess(*args, **kwargs)

    def guessGen(self, *args, **kwargs):
        raise NotImplementedError

    def parseSol(self, x):
        """Parse a solution."""
        return self.prob.parseSol(x)

    def solve(self, x0=None, *args, **kwargs):
        if not self.solver_is_up_to_data:
            self.constructSolver(self.ip)
        if x0 is None:
            x0 = self.guessGen(*args, **kwargs)
        return self.slv.solveGuess(x0)

    def setx0(self, x0):
        """Set bounds for initial states.

        :param x0: ndarray, (ndyn*order,), initial state

        """
        x0 = np.atleast_1d(x0)
        self.x0 = x0
        nf = self.system.nf
        nx0 = self.system.order * self.system.nf
        nx = self.system.nx
        if len(x0) == nx0:
            x0lb = np.concatenate((x0, -getInf(nf)))
            x0ub = np.concatenate((x0, getInf(nf)))
            self.prob.x0bd = [x0lb, x0ub]  # it has to be fixed
        else:
            if len(x0) == nx:
                self.prob.x0bd = [x0, x0]
            else:
                ninf = nx - len(x0)
                self.prob.x0bd = [np.concatenate((x0, -getInf(ninf))), np.concatenate((x0, getInf(ninf)))]

    def setxf(self, xf):
        """Set bounds for final states.

        :param xf: ndarray, (4,) final state

        """
        xf = np.atleast_1d(xf)
        self.xf = xf
        nf = self.system.nf
        nx0 = self.system.order * self.system.nf
        nx = self.system.nx
        if len(xf) == nx0:
            xflb = np.concatenate((xf, -getInf(nf)))
            xfub = np.concatenate((xf, getInf(nf)))
            self.prob.xfbd = [xflb, xfub]  # it has to be fixed
        else:
            if len(xf) == nx:
                self.prob.xfbd = [xf, xf]
            else:
                ninf = nx - len(xf)
                self.prob.xfbd = [np.concatenate((xf, -getInf(ninf))), np.concatenate((xf, getInf(ninf)))]

    def setXbound(self, xlb, xub):
        """Set bounds for state variables."""
        nf = self.system.nf
        nx0 = self.system.order * self.system.nf
        nx = self.system.nx
        assert len(xlb) == len(xub)
        assert len(xlb) == nx0 or len(xlb) == nx
        if len(xlb) == nx0:
            xlb = np.concatenate((xlb, -getInf(nf)))
            xub = np.concatenate((xub, getInf(nf)))
            self.prob.xbd = [xlb, xub]  # it has to be fixed
        else:
            self.prob.xbd = [xlb, xub]

    def setUbound(self, ulb, uub=None):
        """Set bounds on control.

        If uub is None, we are using -ulb and ulb as bounds

        """
        if uub is None:
            ubd = np.array(ulb)
            self.prob.ubd = [-ubd, ubd]
        else:
            self.prob.ubd = [np.array(ulb), np.array(uub)]

    def update(self):
        """Update bounds on x0 and xf"""
        self.prob.__setXbound__()
        if self.ip:
            self.solver_is_up_to_data = False
