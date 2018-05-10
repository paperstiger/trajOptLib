#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2018 Gao Tang <gt70@duke.edu>
#
# Distributed under terms of the MIT license.

"""
testMagicWrapper.py

We test if autograd can be used and use the toy problem to demonstrate so.

"""
import sys, os, time
import numpy as np
import autograd.numpy as numpy
import matplotlib.pyplot as plt
import logging
import autograd.numpy as numpy
sys.path.append('../')
from trajOptLib import daeSystemWrapper
from trajOptLib import daeSystem, trajOptCollocProblem
from trajOptLib import lqrObj
from trajOptLib.utility import showSol
from trajOptLib import snoptConfig, solver


def sysFunOrder2(t, X, u, p):
    x, dx, ddx = X
    return numpy.array([ddx - u[0]])


def sysFunOrder1(t, X, u, p):
    x, v, dx, dv = X
    y0 = dx - v
    y1 = dv - u[0]
    return numpy.array([y0, y1])


def main():
    # prob = constructOrderOne()
    prob = constructOrderTwo()
    # construct a solver for the problem
    cfg = snoptConfig()
    cfg.printLevel = 1
    cfg.printFile = 'test.out'
    cfg.verifyLevel = 3
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    print(rst.flag)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parseSol(rst.sol)
        showSol(sol)


def constructOrderOne():
    """Test the wrapper class for this naive problem"""
    sys = daeSystemWrapper(sysFunOrder1, 4, 1, 0, 2)
    N = 20
    t0 = 0.0
    tf = 10.0
    prob = trajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1.5]), np.array([1.5])]
    prob.x0bd = [np.array([0, 0, -1e20, -1e20]), np.array([0, 0, 1e20, 1e20])]
    prob.xfbd = [np.array([np.pi, 0, -1e20, -1e20]), np.array([np.pi, 0, 1e20, 1e20])]
    lqr = lqrObj(R=np.ones(1))
    prob.addLQRObj(lqr)
    prob.preProcess()  # construct the problem
    return prob


def constructOrderTwo():
    """Test the wrapper class for yet another naive problem."""
    sys = daeSystemWrapper(sysFunOrder2, 3, 1, 0, 1)
    N = 20
    t0 = 0.0
    tf = 10.0
    prob = trajOptCollocProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20, -1e20]), np.array([1e20, 1e20, 1e20])]
    prob.ubd = [np.array([-1.5]), np.array([1.5])]
    prob.x0bd = [np.array([0, 0, -1e20]), np.array([0, 0, 1e20])]
    prob.xfbd = [np.array([np.pi, 0, -1e20]), np.array([np.pi, 0, 1e20])]
    lqr = lqrObj(R=np.ones(1))
    prob.addLQRObj(lqr)
    prob.preProcess()  # construct the problem
    return prob


if __name__ == '__main__':
    main()
