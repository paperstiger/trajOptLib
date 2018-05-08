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
from trajOptLib.io import getOnOffArgs
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
RX = 0
RY = 0
RADIUS = 1


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
        self.rx = RX
        self.ry = RY
        self.radius = RADIUS

    def dyn(self, t, x, u, p, y, G, row, col, rec, needg):
        y[0] = x[4] - u[0] - p[0] * (x[0] - self.rx)
        y[1] = x[5] - u[1] - p[0] * (x[1] - self.ry)
        if needg:
            G[0:2] = 1.0
            G[2:4] = -1.0
            G[4:6] = -p[0]  # w.r.t. x[0], x[1]
            G[6:8] = [-(x[0] - self.rx), -(x[1] - self.ry)]  # w.r.t p
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
        self.x = 0
        self.y = 0
        self.r = 1

    def __callg__(self, X, F, G, row, col, rec, needg):
        x, y, dx, dy, ddx, ddy = X[1:7]
        F[0] = (x-self.x)**2 + (y-self.y)**2 - self.r**2
        if self.with_vel:
            F[1] = (x-self.x)*dx + (y-self.y)*dy
            if self.with_acc:
                F[2] = dx*dx + (x-self.x)*ddx + dy*dy + (y-self.y)*ddy
        if needg:
            G[:2] = [2*(x-self.x), 2*(y-self.y)]
            if self.with_vel:
                G[2:6] = [dx, dy, x - self.x, y - self.y]
                if self.with_acc:
                    G[6:12] = [ddx, ddy, 2*dx, 2*dy, x - self.x, y - self.y]
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
        self.x = RX
        self.y = RY
        self.r = RADIUS

    def __callg__(self, X, F, G, row, col, rec, needg):
        """Calculate the constraints at knot points."""
        x, y, dx, dy, ddx, ddy = X
        F[0] = (x-self.x)**2 + (y-self.y)**2 - self.r**2
        F[1] = (x-self.x)*dx + (y-self.y)*dy
        F[2] = dx*dx + (x-self.x)*ddx + dy*dy + (y-self.y)*ddy
        if needg:
            G[:] = [2*(x-self.x), 2*(y-self.y), dx, dy, x - self.x, y - self.y, ddx, ddy, 2*dx, 2*dy, x - self.x, y-self.y]
            if rec:
                row[:] = [0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
                col[:] = [0, 1, 0, 1, 2, 3, 0, 1, 2, 3, 4, 5]

    def __calc_correction__(self, X, Gamma, Y, Gq, rowq, colq, Gg, rowg, colg, rec, needg):
        """Calculate the constraint J(q)^T \gamma and corresponding Jacobians"""
        x, y, dx, dy, ddx, ddy = X
        gamma = Gamma[0]
        Y[0] = 2 * (x - self.x) * gamma
        Y[1] = 2 * (y - self.x) * gamma
        if needg:
            Gg[:] = [2*(x-self.x), 2*(y-self.y)]
            Gq[:] = [2*gamma, 2*gamma]
            if rec:
                rowg[:] = [0, 1]
                colg[:] = [0, 0]
                rowq[:] = [0, 1]
                colq[:] = [0, 1]


def main():
    args = getOnOffArgs('car', 'debug', 'naive', 'analy', 'warm')
    if args.debug:
        sys.path.append('/home/motion/GAO/ThirdPartyLibs/pycharm-debug.egg')
        import pydevd
        pydevd.settrace('10.197.52.193', port=10000, stdoutToServer=True, stderrToServer=True)
    if args.naive:
        testNaiveCar()
    if args.car:
        testCar(warm=args.warm)
    if args.analy:
        testAnalytic()


def testAnalytic():
    """Test car problem with analytic solution. See what's missing."""
    sys = SimpleOmniCar()  # simple car has no p
    N = 10
    t0 = 0
    tf = 1
    man_constr = circleConstr()
    prob = trajOptManifoldCollocProblem(sys, N, t0, tf, man_constr)
    setBound(prob)

    lqr = lqrObj(R=np.ones(2))
    prob.addLQRObj(lqr)
    prob.preProcess(defect_u=False, defect_p=False)
    # construct solver
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 1
    slv = solver(prob, cfg)
    # generate an guess
    x0 = np.zeros(prob.nx)
    # generate the random initial guess
    t = np.linspace(0, 1, 2*N - 1)
    theta = -np.pi * t ** 3 + 1.5 * np.pi * t ** 2
    dtheta = -3*np.pi*t**2 + 3*np.pi*t
    ddtheta = -6*np.pi*t + 3*np.pi
    # assign them to it
    useX, useU, useP = prob.__parseX__(x0)
    radius = 1
    stheta = np.sin(theta)
    ctheta = np.cos(theta)
    useX[:, 0] = radius * ctheta
    useX[:, 1] = radius * stheta
    useX[:, 2] = -radius * stheta * dtheta
    useX[:, 3] = radius * ctheta * dtheta
    useX[:, 4] = -radius * ctheta * dtheta ** 2 - radius * stheta * ddtheta
    useX[:, 5] = -radius * stheta * dtheta ** 2 + radius * ctheta * ddtheta
    useU[:, 0] = useX[:, 4]
    useU[:, 1] = useX[:, 5]
    # gamma is assumed zero
    # set auxv
    useAuxV = prob.__parseAuxVc__(x0)
    useAuxV[:] = useX[1::2, 2:4]
    psf0 = prob.parseF(x0)
    rst = slv.solveGuess(x0)
    psf = prob.parseF(rst.sol)
    # parse the solution
    sol = prob.parseSol(rst.sol)
    showSol(sol)


def testCar(warm=False):
    """Test the omni-directional car problem."""
    sys = OmniCar()
    # sys = SimpleOmniCar()
    N = 10
    t0 = 0
    tf = 0.9
    man_constr = circleConstr()
    prob = trajOptManifoldCollocProblem(sys, N, t0, tf, man_constr)
    setBound(prob)

    lqr = lqrObj(R=np.ones(2))
    prob.addLQRObj(lqr)
    prob.preProcess(defect_u=True, defect_p=False)
    # construct solver
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 1
    slv = solver(prob, cfg)
    if warm:
        x0 = np.zeros(prob.numSol)
        data = np.load('first_order_traj.npz')
        useX, useU, useP = prob.__parseX__(x0)
        useGamma = prob.__parseGamma__(x0)
        useAuxVc = prob.__parseAuxVc__(x0)
        useX[:] = data['X']
        useU[:] = data['U']
        useP[:] = data['P']
        useGamma[:] = data['gamma']
        useAuxVc[:] = useX[1::2, 2:4]
        psf = prob.parseF(x0)
        rst = slv.solveGuess(x0)
    else:
        rst = slv.solveRand()
    if rst.flag == 1 or rst.flag == 2:
        psf = prob.parseF(rst.sol)
        np.savez('tmp.npz', **psf)
        print(psf.keys())
        print(psf['gamma'])
        # parse the solution
        sol = prob.parseSol(rst.sol)
        showSol(sol)


def testNaiveCar():
    """Use ordinary approach and see"""
    sys = SimpleOmniCar()
    N = 10
    t0 = 0
    tf = 1
    prob = trajOptCollocProblem(sys, N, t0, tf)
    setBound(prob)

    lqr = lqrObj(R=np.ones(2))
    prob.addLQRObj(lqr)

    constr = collocCircleConstr(with_vel=True, with_acce=True)
    prob.addConstr(constr, path=True)
    prob.preProcess(colloc_constr_is_on=False, defect_p=False, defect_u=False)
    # construct solver
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 1
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    psf = prob.parseF(rst.sol.copy())
    if rst.flag == 1:
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        x = sol['x']
        u = sol['u']
        fig, ax = plt.subplots(2, 2)
        ax[0][0].plot(x[:, 0], x[:, 1])
        ax[0][1].plot(u[:, 0])
        ax[0][1].plot(u[:, 1])
        # check velocity
        velvio = (x[:, 0] - constr.x) * x[:, 2] + (x[:, 1] - constr.y) * x[:, 3]
        print(velvio)
        accevio = x[:, 2] ** 2 + (x[:, 0] - constr.x) * x[:, 4] + x[:, 3] ** 2 + (x[:, 1] - constr.y) * x[:, 5]
        print(accevio)
        plt.show()
        # showSol(sol)


def setBound(prob):
    x0lb = np.array([1, 0, 0, 0, -1e20, -1e20])
    x0ub = np.array([1, 0, 0, 0, 1e20, 1e20])
    prob.x0bd = [x0lb, x0ub]
    xflb = np.array([0, 1, 0, 0, -1e20, -1e20])
    xfub = np.array([0, 1, 0, 0, 1e20, 1e20])
    prob.xfbd = [xflb, xfub]
    prob.xbd = [-1e20*np.ones(6), 1e20*np.ones(6)]
    prob.ubd = [-1e20*np.ones(2), 1e20*np.ones(2)]
    if prob.numP > 0:
        prob.pbd = [-1e20*np.ones(1), 1e20*np.ones(1)]


if __name__ == '__main__':
    main()
