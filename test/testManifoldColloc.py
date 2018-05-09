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

from carCommon import CircleConstr, SecondOrderOmniCar


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
RX = 0
RY = 0
RADIUS = 1


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
    sys = SecondOrderOmniCar()  # simple car has no p
    N = 10
    t0 = 0
    tf = 1
    man_constr = CircleConstr()
    prob = trajOptManifoldCollocProblem(sys, N, t0, tf, man_constr)
    setBound(prob)

    lqr = lqrObj(R=np.ones(2), P=np.ones(1))
    prob.addLQRObj(lqr)
    prob.preProcess(defect_u=True, defect_p=False)
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
    psf0 = prob.parseF(x0)
    rst = slv.solveGuess(x0)
    psf = prob.parseF(rst.sol)
    # solve the problem
    rst = slv.solveGuess(x0)
    print(rst.flag)
    # parse the solution
    sol = prob.parseSol(rst.sol)
    showSol(sol)


def testCar(warm=False):
    """Test the omni-directional car problem."""
    sys = SecondOrderOmniCar()
    # sys = SimpleOmniCar()
    N = 10
    t0 = 0
    tf = 0.9
    man_constr = CircleConstr()
    prob = trajOptManifoldCollocProblem(sys, N, t0, tf, man_constr)
    setBound(prob)

    lqr = lqrObj(R=np.ones(2), P=1*np.ones(1))
    prob.addLQRObj(lqr)
    prob.preProcess(defect_u=True, defect_p=False, gamma_bound=1)
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
        useX[:] = data['X']
        useU[:] = data['U']
        useP[:] = data['P']
        useGamma[:] = data['gamma']
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
