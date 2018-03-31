#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
twoPhaseDemo.py

Test the multiple-phase solver using this example.
We consider a 2D linear car problem.
Basically, the car has to reach a certain point with free velocity.
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
sys.path.append('../')
from trajOptLib.io import getOnOffArgs
from trajOptLib import daeSystem, trajOptCollocProblem
from trajOptLib import nonLinearPointObj, linearPointObj, linearPointConstr
from trajOptLib import nonLinearPointConstr
from trajOptLib import lqrObj
from trajOptLib import snoptConfig, solver
from trajOptLib import TrajOptMultiPhaseCollocProblem, LinearConnectConstr, NonLinearConnectConstr
from trajOptLib.utility import showSol
from scipy.sparse import coo_matrix


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


def main():
    args = getOnOffArgs('one', 'two', 'oneleg')
    if args.one:
        testOne()
    if args.oneleg:
        testOneLeg()
    if args.two:
        testTwo()


class orderTwoModel(daeSystem):
    """A class with dynamics :math:`\ddot{x}=u; \ddot{y}=v`"""
    def __init__(self):
        daeSystem.__init__(self, 6, 2, 0, 2, 4)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[5] - u[0]
        y[1] = x[6] - u[1]
        if needg:
            G[0:2] = 1
            G[2:4] = -1
            if rec:
                row[:] = [0, 1, 0, 1]
                col[:] = [5, 6, 7, 8]


class orderOneModel(daeSystem):
    """The 2D model in first order :math:`\dot{x}=v_x, \dot{y}=v_y, \dot{v_x}=u, \dot{v_y}=v`"""
    def __init__(self):
        daeSystem.__init__(self, 8, 2, 0, 4, 8)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[:2] = x[4:6] - x[2:4]
        y[2:4] = x[6:] - u
        if needg:
            G[:4] = 1
            G[4:8] = -1
            if rec:
                row[:4] = np.arange(4)
                row[4:] = np.arange(4)
                col[:] = [5, 6, 7, 8, 3, 4, 9, 10]


class connectStateConstr(LinearConnectConstr):
    """Linear constraint such that time is continuous; x, y continuous."""
    def __init__(self, N, dimq):
        a1 = np.eye(1 + dimq)
        a2 = -np.eye(1 + dimq)
        LinearConnectConstr.__init__(self, 0, 1, a1, a2)  # this makes time continuous


class connectVelocityConstr(NonLinearConnectConstr):
    """Nonlinear constraint that u must change direction. x1 and x2 are both of length 11"""
    def __init__(self):
        NonLinearConnectConstr.__init__(self, 0, 1, 2, nG=4)

    def __callg__(self, x1, x2, F, G, row, col, rec, needg, addx=None):
        """Override this function."""
        u1 = x1[9:]
        u2 = x2[9:]
        F[:] = u1 + u2
        if needg:
            G[:] = [1, 1, 1, 1]
            if rec:
                row[:] = [0, 1, 0, 1]
                col[:] = [9, 10, 20, 21]

    def __callhumang__(self, x1, x2, needg, addx=None):
        u1 = x1[9:]
        u2 = x2[9:]
        F = u1 + u2
        if needg:
            G1 = coo_matrix(([1, 1], ([0, 1], [9, 10])))
            G2 = coo_matrix(([1, 1], ([0, 1], [9, 10])))
            return F, G1, G2
        else:
            return F, None, None


class objAvoidConstr(nonLinearPointConstr):
    """A constraint that at end of first phase, it is at least 0.2 distance w.r.t center"""
    def __init__(self, N):
        self.midpoint = np.array([0.5, 0.5])
        self.distance = 0.2
        nonLinearPointConstr.__init__(self, N - 1, 1, 8, 2, np=0, lb=[0], ub=[1e20], nG=2)

    def __callg__(self, x, y, G, row, col, rec, needg):
        dis = x[1:3] - self.midpoint
        y[0] = np.dot(dis, dis) - self.distance ** 2
        if needg:
            G[:2] = 2 * dis
            if rec:
                row[:2] = [0, 0]
                col[:2] = [1, 2]


def testOneLeg():
    """Test order one, leg one"""
    sys = orderOneModel()
    N = 20
    t0 = 0.0
    tf = 10.0
    prob1 = trajOptCollocProblem(sys, N, t0, tf)
    xlb = -1e20 * np.ones(8)
    xub = 1e20 * np.ones(8)
    ulb = -1.5 * np.ones(2)
    uub = 1.5 * np.ones(2)
    x0lb = np.concatenate((np.zeros(4), -1e20*np.ones(4)))
    x0ub = np.concatenate((np.zeros(4), 1e20*np.ones(4)))
    xflb = np.concatenate((np.ones(2), np.zeros(2), -1e20*np.ones(4)))
    xfub = np.concatenate((np.ones(2), np.zeros(2), 1e20*np.ones(4)))
    prob1.xbd = [xlb, xub]
    prob1.ubd = [ulb, uub]
    prob1.x0bd = [x0lb, x0ub]
    prob1.xfbd = [xflb, xfub]
    # define objective function
    lqr = lqrObj(R=np.ones(2))
    prob1.addLQRObj(lqr)
    # add several constraints
    obj_avoid = objAvoidConstr(N)
    prob1.addConstr(obj_avoid)
    # ready to construct this problem
    prob1.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 3
    slv = solver(prob1, cfg)
    rst = slv.solveRand()
    print(rst.flag)
    if rst.flag == 1:
        print(rst.sol)
        # parse the solution
        sol = prob1.parseSol(rst.sol.copy())
        showSol(sol)


def testOne():
    """Test order one pendulum case, this is seen everywhere."""
    sys = orderOneModel()
    N = 20
    t0 = 0.0
    tf = 10.0
    tmid = 5.0
    prob1 = trajOptCollocProblem(sys, N, t0, [t0+1, tmid])
    prob2 = trajOptCollocProblem(sys, N, tmid, tf) # [t0, tf], tf)
    xlb = -1e20 * np.ones(8)
    xub = 1e20 * np.ones(8)
    ulb = -1.5 * np.ones(2)
    uub = 1.5 * np.ones(2)
    x0lb = np.concatenate((np.zeros(4), -1e20*np.ones(4)))
    x0ub = np.concatenate((np.zeros(4), 1e20*np.ones(4)))
    xflb = np.concatenate((np.ones(2), np.zeros(2), -1e20*np.ones(4)))
    xfub = np.concatenate((np.ones(2), np.zeros(2), 1e20*np.ones(4)))
    prob1.xbd = [xlb, xub]
    prob1.ubd = [ulb, uub]
    prob1.x0bd = [x0lb, x0ub]
    prob1.xfbd = [xlb, xub]
    prob2.xbd = [xlb, xub]
    prob2.ubd = [ulb, uub]
    prob2.x0bd = [xlb, xub]
    prob2.xfbd = [xflb, xfub]
    # define objective function
    lqr = lqrObj(R=np.ones(2))
    prob1.addLQRObj(lqr)
    prob2.addLQRObj(lqr)
    # add several constraints
    obj_avoid = objAvoidConstr(2*N - 2)
    prob1.addConstr(obj_avoid)
    # ready to construct this problem
    prob = TrajOptMultiPhaseCollocProblem([prob1, prob2], addx=None)
    # add connect constraints
    constr1 = connectStateConstr(N, 4)  # since dimq = 4
    constr2 = connectVelocityConstr()
    prob.add_constr(constr1)
    prob.add_constr(constr2)
    # ready to solve this problem. It is quite complicated to construct this problem, how come?
    # TODO: add support for easy-to-construct constraints. 
    # For example, linear constraints can be constructed by giving matrix
    # nonlinear constraints can be constructed by functions. These functions can be auto-diffed
    prob.pre_process()  # construct the problem
    # construct a solver for the problem
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 3
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    print(rst.flag)
    if rst.flag == 1:
        # print(rst.sol)
        # parse the solution
        sol = prob.parse_sol(rst)
        showSol(sol)
        fig, ax = plt.subplots()
        phase1 = sol['phases'][0]
        phase2 = sol['phases'][1]
        ax.plot(phase1['x'][:, 0], phase1['x'][:, 1])
        ax.plot(phase1['x'][-1, 0], phase1['x'][-1, 1], marker='*', markersize=5)
        ax.plot(phase2['x'][:, 0], phase2['x'][:, 1])
        ax.plot(phase2['x'][-1, 0], phase2['x'][-1, 1], marker='*', markersize=5)
        plt.show()


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


if __name__ == '__main__':
    main()
