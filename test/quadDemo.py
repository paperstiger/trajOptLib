#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
quadDemo.py

A demo of quadrotor which is quite challenging.
Use the classical model used everywhere else.
"""
import sys, os, time
import numpy as np
import logging
from pyLib.io import getOnOffArgs
sys.path.append('../')
from trajOptLib.trajOptBase import system, nonLinearPointObj, nonLinearObj, lqrObj
from trajOptLib.trajOptProblem import trajOptProblem
from trajOptLib.libsnopt import snoptConfig, probFun, solver
from trajOptLib.utility import showSol
from scipy.sparse import coo_matrix
from libRotor import Rotor


if False:
    import pydevd
    pydevd.settrace('10.197.84.153', port=10000, stdoutToServer=True, stderrToServer=True)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class quadrotor(system, Rotor):
    """A class derived from system and Rotor"""
    def __init__(self):
        system.__init__(self, 12, 4, 0, 'Euler')
        Rotor.__init__(self)

    def Jdyn(self, t, x, u, p=None):
        f = np.zeros(self.nx)
        J = np.zeros((self.nx, self.nx + self.nu + 1 + self.np), order='F')
        Rotor.dyn(self, t, x, u, f, J)
        J = np.ascontiguousarray(J)
        return f, J


class quadCost(nonLinearObj):
    """A quadratic cost on control."""
    def __init__(self, N, dimx, dimu):
        lenSol = N * (dimx + dimu)
        nonLinearObj.__init__(self, lenSol, 'user', nG=N*dimu)
        self.R = 1.0
        self.N = N
        self.dimx = dimx
        self.dimu = dimu

    def __callf__(self, x, y):
        u = x[3]
        y[0] = u * self.R * u

    def __callg__(self, x, y, G, row, col, rec, needg):
        u = np.reshape(x[self.N*self.dimx:], (self.N, self.dimu))
        y[0] = np.sum(u**2)
        if needg:
            G[:self.N*self.dimu] = 2.0 * u.flatten()
            if rec:
                row[:self.N*self.dimu] = 0
                col[:self.N*self.dimu] = np.arange(self.N*self.dimx, self.N*(self.dimx+self.dimu))


def main():
    sys = quadrotor()
    N = 40
    dimx, dimu = sys.nx, sys.nu
    cost = quadCost(N, sys.nx, sys.nu)
    t0 = 0.0
    tf = 5.0
    prob = trajOptProblem(sys, N, t0, tf, gradmode=True)
    prob.xbd = [-1e20*np.ones(sys.nx), 1e20*np.ones(sys.nx)]
    prob.ubd = [0*np.ones(sys.nu), 4*np.ones(sys.nu)]
    prob.x0bd = [np.zeros(sys.nx), np.zeros(sys.nx)]
    prob.xfbd = [np.zeros(sys.nx), np.zeros(sys.nx)]
    prob.xfbd[0][:3] = 5
    prob.xfbd[1][:3] = 5
    if False:
        prob.addNonLinearObj(cost)
    else:
        lqr = lqrObj(R=np.ones(4))
        prob.addLQRObj(lqr)
    prob.preProcess()
    # construct a solver for the problem
    cfg = snoptConfig()
    # cfg.printFile = 'Toy0.out'
    cfg.printLevel = 1
    cfg.addIntOption('Verify Level', 0)
    slv = solver(prob, cfg)
    guessx = np.zeros(prob.nx)
    straightx = np.reshape(guessx[:N*dimx], (N, dimx))
    for i in range(3):
        straightx[:, i] = np.linspace(0, prob.xfbd[0][i], N)
    guessx[N*dimx:-1] = np.random.random(N * dimu)
    rst = slv.solveGuess(guessx)
    print(rst.flag)
    if False: #rst.flag == 1:
        print(rst.sol)
        # parse the solution
        sol = prob.parseSol(rst.sol.copy())
        showSol(sol)


if __name__ == '__main__':
    main()
