#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
simpleDemo.py

Use a 1D problem to show how it works.
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
sys.path.append('../')
from trajOptLib.io import getOnOffArgs
from trajOptLib.trajOptBase import system, nonPointObj, lqrObj
from trajOptLib.trajOptProblem import trajOptProblem
from trajOptLib.libsnopt import snoptConfig, probFun, solver
from trajOptLib.utility import showSol
from trajOptLib import ipSolver
from scipy.sparse import coo_matrix


if False:
    import pydevd
    pydevd.settrace('10.197.84.153', port=10000, stdoutToServer=True, stderrToServer=True)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class oneDcase(system):
    """Test second order system."""
    def __init__(self):
        system.__init__(self, 2, 1, 0, 'Euler')

    def dyn(self, t, x, u, p=None):
        y1 = x[1]
        y2 = u[0]
        return np.array([y1, y2])

    def Jdyn(self, t, x, u, p=None):
        y1 = x[1]
        y2 = u[0]
        J = np.zeros((2, 4))
        J[0, 2] = 1.0
        J[1, 3] = 1.0
        return np.array([y1, y2]), coo_matrix(J)


class pendulum(system):
    """Test pendulum nonlinearity."""
    def __init__(self):
        system.__init__(self, 2, 1, 0, 'Euler')

    def Jdyn(self, t, x, u, p=None):
        y1 = x[1]
        y2 = u[0] / 5 - 5 * np.sin(x[0])
        y = np.array([y1, y2])
        J = coo_matrix(([1, -5 * np.cos(x[0]), 0.2], ([0, 1, 1], [2, 1, 3])))
        return y, J


class quadCost(nonPointObj):
    """A quadratic cost on control."""
    def __init__(self):
        nonPointObj.__init__(self, -1, 2, 1, 0, 'user', nG=1)
        self.R = 1.0

    def __callf__(self, x, y):
        u = x[3]
        y[0] = u * self.R * u

    def __callg__(self, x, y, G, row, col, rec, needg):
        u = x[3]
        y[0] = u * self.R * u
        if needg:
            G[0] = 2 * self.R * u
            if rec:
                row[0] = 0
                col[0] = 3


def main():
    args = getOnOffArgs('fd', 'grad', 'pen', 'lqr', 'ip')
    if args.fd:
        fdmode(args.lqr)
    if args.grad:
        gradmode(args.lqr)
    if args.pen:
        penMode(args.lqr)
    if args.ip:
        ipMode(args.lqr)


def penMode(lqr):
    """Run pendulum swing-up problem"""
    sys = pendulum()
    cost = quadCost()
    N = 40
    t0 = 0.0
    tf = 9.0
    prob = trajOptProblem(sys, N, t0, tf, gradmode=True)
    prob.xbd = [np.array([-1e20, -1e20]), np.array([1e20, 1e20])]
    prob.ubd = [np.array([-1.0]), np.array([1.0])]
    prob.x0bd = [np.array([0, 0]), np.array([0, 0])]
    prob.xfbd = [np.array([np.pi, 0]), np.array([np.pi, 0])]
    if not lqr:
        prob.addNonPointObj(cost, True)  # add a path cost
    else:
        lqr = lqrObj(R=np.ones(1))
        prob.addLQRObj(lqr)
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = snoptConfig()
    cfg.printLevel = 1
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    print(rst.flag)
    if rst.flag == 1:
        print(rst.sol)
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        showSol(sol)


def gradmode(lqr):
    """Solve the simple problem with gradient."""
    sys = oneDcase()
    cost = quadCost()
    N = 20
    t0 = 0.0
    tf = 2.0
    prob = trajOptProblem(sys, N, t0, tf, gradmode=True)
    prob.xbd = [np.array([-1e20, -1e20]), np.array([1e20, 1e20])]
    prob.ubd = [np.array([-1e20]), np.array([1e20])]
    prob.x0bd = [np.array([0, 0]), np.array([0, 0])]
    prob.xfbd = [np.array([1, 0]), np.array([1, 0])]
    if not lqr:
        prob.addNonPointObj(cost, True)  # add a path cost
    else:
        lqr = lqrObj(R=np.ones(1))
        prob.addLQRObj(lqr)
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = snoptConfig()
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    print(rst.flag, rst.sol)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        showSol(sol)


def ipMode(lqr):
    """Use ipopt to solve this simple problem and see what is happening."""
    sys = oneDcase()
    cost = quadCost()
    N = 20
    t0 = 0.0
    tf = 2.0
    prob = trajOptProblem(sys, N, t0, tf, gradmode=True)
    prob.xbd = [np.array([-1e20, -1e20]), np.array([1e20, 1e20])]
    prob.ubd = [np.array([-1e20]), np.array([1e20])]
    prob.x0bd = [np.array([0, 0]), np.array([0, 0])]
    prob.xfbd = [np.array([1, 0]), np.array([1, 0])]
    if not lqr:
        prob.addNonPointObj(cost, True)  # add a path cost
    else:
        lqr = lqrObj(R=np.ones(1))
        prob.addLQRObj(lqr)
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    slv = ipSolver(prob)
    rst = slv.solveRand()
    print(rst.flag, rst.sol)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        showSol(sol)


def fdmode(lqr):
    """Solve the simple problem with finite difference."""
    sys = oneDcase()
    cost = quadCost()
    x0 = np.random.random(2)
    u = np.random.random(1)
    print(sys.dyn(0, x0, u))
    N = 20
    t0 = 0.0
    tf = 2.0
    prob = trajOptProblem(sys, N, t0, tf, gradmode=False)
    prob.xbd = [np.array([-1e20, -1e20]), np.array([1e20, 1e20])]
    prob.ubd = [np.array([-1e20]), np.array([1e20])]
    prob.x0bd = [np.array([0, 0]), np.array([0, 0])]
    prob.xfbd = [np.array([1, 0]), np.array([1, 0])]
    if not lqr:
        prob.addNonPointObj(cost, True)  # add a path cost
    else:
        lqr = lqrObj(R=np.ones(1))
        prob.addLQRObj(lqr)
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = snoptConfig()
    cfg.printLevel = 1
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    print(rst.flag, rst.sol)
    # parse the sulution
    sol = prob.parseSol(rst.sol.copy())
    showSol(sol)


if __name__ == '__main__':
    main()
