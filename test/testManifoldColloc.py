#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testManifoldColloc.py

Test the collocation approach with manifold support.
Basically I will define a second-order omni-directional car. 
It is constrained to move on a circle. We move from (0, 0) to (1, 1) on a circle (x-0.5)^2+(y-0.5)^2=0.5
"""
import sys, os, time
import numpy as np
import matplotlib.pyplot as plt
import logging
sys.path.append('../')
from trajOptLib.fakeio import getOnOffArgs
from trajOptLib import daeSystem, trajOptCollocProblem
from trajOptLib import manifoldConstr, trajOptManifoldCollocProblem
from trajOptLib import nonLinearPointConstr
from trajOptLib import nonLinearPointObj, linearPointObj, linearPointConstr
from trajOptLib import lqrObj
from trajOptLib import snoptConfig, solver
from trajOptLib.utility import showSol
from scipy.sparse import coo_matrix


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SimpleOmniCar(daeSystem):
    """A car without constraint. you use control to guarantee it"""
    def __init__(self):
        daeSystem.__init__(self, 6, 2, 0, 2, 4)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[4] - u[0]
        y[1] = x[5] - u[1]
        if needg:
            G[0:2] = 1.0
            G[2:4] = -1.0
            if rec:
                row[:] = [0, 1, 0, 1]
                col[:] = [5, 6, 7, 8]


class OmniCar(daeSystem):
    """A omni-car system in 2D space.

    It has simple dynamics \ddot{x}=u_x; \ddot{y}=u_y

    """
    def __init__(self):
        daeSystem.__init__(self, 6, 2, 1, 2, 8)

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[4] - u[0] - 2 * p[0] * (x[0] - 0.5)
        y[1] = x[5] - u[1] - 2 * p[0] * (x[1] - 0.5)
        if needg:
            G[0:2] = 1.0
            G[2:4] = -1.0
            G[4:6] = -p[0]  # w.r.t. x[0], x[1]
            G[6:8] = [-2*(x[0] - 0.5), -2*(x[1] - 0.5)]  # w.r.t p
            if rec:
                row[:] = [0, 1, 0, 1, 0, 1, 0, 1]
                col[:] = [5, 6, 7, 8, 1, 2, 9, 9]


class collocCircleConstr(nonLinearPointConstr):
    def __init__(self, with_vel=False, with_acce=False):
        nG = 2
        nf = 1
        self.with_vel = with_vel
        self.with_acc = with_acce
        if with_vel:
            nG += 4
            nf += 1
            if with_acce:
                nG += 6
                nf += 1
        nonLinearPointConstr.__init__(self, -1, nf, 6, 2, nG=nG)

    def __callg__(self, X, F, G, row, col, rec, needg):
        x, y, dx, dy, ddx, ddy = X[1:7]
        F[0] = (x-0.5)**2 + (y-0.5)**2 - 0.5
        if self.with_vel:
            F[1] = (x-0.5)*dx + (y-0.5)*dy
            if self.with_acc:
                F[2] = dx*dx + (x-0.5)*ddx + dy*dy + (y-0.5)*ddy
        if needg:
            G[:2] = [2*(x-0.5), 2*(y-0.5)]
            if self.with_vel:
                G[2:6] = [dx, dy, x - 0.5, y - 0.5]
                if self.with_acc:
                    G[6:12] = [ddx, ddy, 2*dx, 2*dy, x - 0.5, y - 0.5]
            if rec:
                row[:2] = [0, 0]
                col[:2] = [1, 2]
                if self.with_vel:
                    row[2:6] = [1, 1, 1, 1]
                    col[2:6] = [1, 2, 3, 4]
                    if self.with_acc:
                        row[6:12] = 2
                        col[6:12] = np.arange(6) + 1


class circleConstr(manifoldConstr):
    def __init__(self, with_vel=False):
        manifoldConstr.__init__(self, 6, 1, 2, constr_order=0, nG=12, nnzJx=2, nnzJg=2)
        self.with_vel = with_vel

    def __callg__(self, X, F, G, row, col, rec, needg):
        """Calculate the constraints at knot points."""
        x, y, dx, dy, ddx, ddy = X
        F[0] = (x-0.5)**2 + (y-0.5)**2 - 0.5
        F[1] = (x-0.5)*dx + (y-0.5)*dy
        F[2] = dx*dx + (x-0.5)*ddx + dy*dy + (y-0.5)*ddy
        if needg:
            G[:] = [2*(x-0.5), 2*(y-0.5), dx, dy, x - 0.5, y - 0.5, ddx, ddy, 2*dx, 2*dy, x - 0.5, y-0.5]
            if rec:
                row[:] = [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
                col[:] = [0, 1, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5]

    def __calc_correction__(self, X, Gamma, Y, Gq, rowq, colq, Gg, rowg, colg, rec, needg):
        """Calculate the constraint J(q)^T \gamma and corresponding Jacobians"""
        x, y, dx, dy, ddx, ddy = X
        gamma = Gamma[0]
        Y[0] = 2 * (x - 0.5) * gamma
        Y[1] = 2 * (y - 0.5) * gamma
        if needg:
            Gg[:] = [2*(x-0.5), 2*(y-0.5)]
            Gq[:] = [2*gamma, 2*gamma]
            if rec:
                rowg[:] = [0, 1]
                colg[:] = [0, 0]
                rowq[:] = [0, 1]
                colq[:] = [0, 1]


def main():
    args = getOnOffArgs('car', 'debug', 'naive')
    if args.debug:
        sys.path.append('/home/motion/GAO/ThirdPartyLibs/pycharm-debug.egg')
        import pydevd
        pydevd.settrace('10.197.47.90', port=10000, stdoutToServer=True, stderrToServer=True)
    if args.naive:
        testNaiveCar()
    if args.car:
        testCar()


def testCar():
    """Test the omni-directional car problem."""
    # sys = OmniCar()
    sys = SimpleOmniCar()
    N = 10
    t0 = 0
    tf = 1
    man_constr = circleConstr()
    prob = trajOptManifoldCollocProblem(sys, N, t0, tf, man_constr)
    prob.xbd = [-1e20*np.ones(6), 1e20*np.ones(6)]
    prob.ubd = [-1e20*np.ones(2), 1e20*np.ones(2)]
    # prob.pbd = [np.zeros(1), 1e20*np.ones(1)]
    x0bd = np.array([0, 0, 0, 0, -1e20, -1e20])
    prob.x0bd = [x0bd, -x0bd]
    xflb = np.array([1, 1, 0, 0, -1e20, -1e20])
    xfub = np.array([1, 1, 0, 0, 1e20, 1e20])
    prob.xfbd = [xflb, xfub]
    lqr = lqrObj(R=np.ones(2))
    prob.addLQRObj(lqr)
    prob.preProcess()
    # construct solver
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 3
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    if True:  #rst.flag == 1 or rst.flag == 2:
        psf = prob.parseF(rst.sol.copy())
        print(psf['Gamma'])
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        showSol(sol)


def testNaiveCar():
    """Use ordinary approach and see"""
    sys = OmniCar()
    N = 10
    t0 = 0
    tf = 1
    man_constr = circleConstr()
    prob = trajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [-1e20*np.ones(6), 1e20*np.ones(6)]
    prob.ubd = [-1e20*np.ones(2), 1e20*np.ones(2)]
    x0bd = np.array([0, 0, 0, 0, -1e20, -1e20])
    prob.x0bd = [x0bd, -x0bd]
    xflb = np.array([1, 1, 0, 0, -1e20, -1e20])
    xfub = np.array([1, 1, 0, 0, 1e20, 1e20])
    prob.xfbd = [xflb, xfub]
    lqr = lqrObj(R=np.ones(2))
    prob.addLQRObj(lqr)
    constr = collocCircleConstr(with_vel=True, with_acce=True)
    prob.addConstr(constr, path=True)
    prob.preProcess()
    # construct solver
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 3
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    sol = prob.parseF(rst.sol.copy())
    if rst.flag == 1:
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        x = sol['x']
        fig, ax = plt.subplots()
        ax.plot(x[:, 0], x[:, 1])
        # check velocity
        velvio = (x[:, 0] - 0.5) * x[:, 2] + (x[:, 1] - 0.5) * x[:, 3]
        print(velvio)
        accevio = x[:, 2] ** 2 + (x[:, 0] - 0.5) * x[:, 4] + x[:, 3] ** 2 + (x[:, 1] - 0.5) * x[:, 5]
        print(accevio)
        plt.show()
        # showSol(sol)


if __name__ == '__main__':
    main()
