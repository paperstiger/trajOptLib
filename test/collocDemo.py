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
from trajOptLib import daeSystem, trajOptCollocProblem
from trajOptLib import nonLinearPointObj, linearPointObj, linearPointConstr
from trajOptLib import lqrObj
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


class orderOneOneD(daeSystem):
    def __init__(self):
        daeSystem.__init__(self, 4, 1, 0, 2, 4)  # ddx = u

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[2] - x[1]
        y[1] = x[3] - u[0]
        if needg:
            G[0] = 1
            G[1] = -1
            G[2] = 1
            G[3] = -1
            if rec:
                row[:2] = 0
                row[2:] = 1
                col[:] = [3, 2, 4, 5]


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


class orderOnePendulum(daeSystem):
    """Pendulum with order 1"""
    def __init__(self):
        daeSystem.__init__(self, 4, 1, 0, 2, 5)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        theta, omega, dtheta, domega = x
        y[0] = dtheta - omega
        y[1] = domega - u[0] / 5.0 + 0.5 * np.sin(x[0])
        if needg:
            G[0] = 1
            G[1] = -1
            G[2] = 1
            G[3] = 0.5 * np.cos(x[0])
            G[4] = -0.2
            if rec:
                row[:2] = 0
                row[2:] = 1
                col[:] = [3, 2, 4, 1, 5]


def main():
    args = getOnOffArgs('oned', 'pen', 'lqr', 'ip', 'linear', 'orderone')
    if args.oned:
        testOneD()
    if args.pen:
        testPen()
    if args.linear:
        testLinear()
    if args.orderone:
        testOrderOne(args)


def testOrderOne(args):
    """Test order one pendulum case, this is seen everywhere."""
    if args.pen:
        sys = orderOnePendulum()
    else:
        sys = orderOneOneD()
    N = 20
    t0 = 0.0
    tf = 20.0
    prob = trajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1.5]), np.array([1.5])]
    prob.x0bd = [np.array([0, 0, -1e20, -1e20]), np.array([0, 0, 1e20, 1e20])]
    prob.xfbd = [np.array([np.pi, 0, -1e20, -1e20]), np.array([np.pi, 0, 1e20, 1e20])]
    lqr = lqrObj(R=np.ones(1))
    prob.addLQRObj(lqr)
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 3
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    print(rst.flag)
    if rst.flag == 1:
        print(rst.sol)
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        showSol(sol)


def testLinear():
    """Test 1d problem with linear constraints and linear objective"""
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
    A = np.zeros(5)
    A[1] = 1
    A[2] = 1  # so it basically does nothing
    linPntObj = linearPointObj(0, A, 3, 1, 0)
    prob.addLinearPointObj(linPntObj)
    # add linear constraint that x is increasing
    A = np.zeros(5)
    A[1] = 1
    lb = np.zeros(1)
    ub = np.ones(1)
    linPntCon = linearPointConstr(-1, A, lb, ub)
    prob.addLinearPointConstr(linPntCon, True)
    # we want mid point to be close to 0.8
    wantState = np.array([0.8, 0])
    pntObj = pointObj(N, wantState)
    prob.addObj(pntObj)
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


class pointObj(nonLinearPointObj):
    """A objective function to make mid point close to a selected point"""
    def __init__(self, N, state):
        nonLinearPointObj.__init__(self, 15, 3, 1, 0, 'user', 2)
        self.state = state
        self.weight = 100

    def __callg__(self, x, F, G, row, col, rec, needg):
        dx = x[1:3] - self.state
        F[0] = self.weight * np.sum(dx ** 2)
        if needg:
            G[:2] = self.weight * 2 * dx
            if rec:
                row[:2] = 0
                col[:2] = np.arange(1, 3)


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
    N = 20
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