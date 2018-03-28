#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
collocDemo.py

Use the collocation version of problem.
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
sys.path.append('../')
from trajOptLib.io import getOnOffArgs
from trajOptLib import daeSystem, trajOptCollocProblem, lqrObj
from trajOptLib import snoptConfig, solver
from trajOptLib.utility import showSol
from scipy.sparse import coo_matrix


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


class oneDcase(daeSystem):
    def __init__(self):
        daeSystem.__init__(self, 3, 1, 0, 1, 2)  # ddx = u

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[2] - u[0]
        if needg:
            G[0] = 1
            G[1] = -1
            if rec:
                row[0] = 0
                row[1] = 0
                col[0] = 3
                col[1] = 4


class pendulum(daeSystem):
    """Test pendulum nonlinearity."""
    def __init__(self):
        daeSystem.__init__(self, 3, 1, 0, 1, 3)  # ddq = u/5 - .5*sin(q)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[2] - u[0] / 5. + 0.5 * np.sin(x[0])
        if needg:
            G[0] = 0.5 * np.cos(x[0])
            G[1] = 1
            G[2] = -0.2
            if rec:
                row[:3] = 0
                col[:3] = [1, 3, 4]


def main():
    args = getOnOffArgs('oned', 'pen', 'lqr', 'ip')
    if args.oned:
        testOneD()
    if args.pen:
        testPen()


def testOneD():
    """Test solving one-dim problem using collocation approach"""
    sys = oneDcase()
    N = 10
    t0 = 0.0
    tf = 2.0
    prob = trajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1e20]), np.array([1e20])]
    prob.x0bd = [np.array([0, 0, -1e20]), np.array([0, 0, 1e20])]
    prob.xfbd = [np.array([1, 0, -1e20]), np.array([1, 0, 1e20])]
    lqr = lqrObj(R=np.ones(1))
    prob.addLQRObj(lqr)
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = snoptConfig()
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 3
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    print(rst.flag, rst.sol)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        showSol(sol)


def testPen():
    """Test solving pendulum swing up problem using collocation approach"""
    sys = pendulum()
    N = 40
    t0 = 0.0
    tf = 20.0
    prob = trajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1.5]), np.array([1.5])]
    prob.x0bd = [np.array([0, 0, -1e20]), np.array([0, 0, 1e20])]
    prob.xfbd = [np.array([np.pi, 0, -1e20]), np.array([np.pi, 0, 1e20])]
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


if __name__ == '__main__':
    main()
