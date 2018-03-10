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
from trajOptBase import system, pointObj
from trajOptProblem import trajOptProblem
from libsnopt import snoptConfig, probFun, solver
from utility import showSol


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


class quadCost(pointObj):
    """A quadratic cost on control."""
    def __init__(self):
        self.R = 1.0

    def __callf__(self, x, y):
        u = x[3]
        y[0] = u * self.R * u


def main():
    sys = oneDcase()
    cost = quadCost()
    N = 20
    t0 = 0.0
    tf = 2.0
    prob = trajOptProblem(sys, N, t0, tf)
    prob.xbd = [np.array([-1e20, -1e20]), np.array([1e20, 1e20])]
    prob.ubd = [np.array([-1e20, -1e20]), np.array([1e20, 1e20])]
    prob.x0bd = [np.array([0, 0]), np.array([0, 0])]
    prob.xfbd = [np.array([1, 0]), np.array([1, 0])]
    prob.addNonPointObj(cost, True)  # add a path cost
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = snoptConfig()
    slv = solver(prob, cfg)
    rst = slv.solveRand()
    print(rst.flag, rst.sol)
    # parse the sulution
    sol = prob.parseSol(rst.sol.copy())
    showSol(sol)


if __name__ == '__main__':
    main()
