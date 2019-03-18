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
from trajoptlib.io import get_onoff_args
from trajoptlib import System, NonLinearPointObj, LqrObj
from trajoptlib import TrajOptProblem
from trajoptlib import OptConfig, OptSolver
from trajoptlib.utility import show_sol
from scipy.sparse import coo_matrix


class OneDcase(System):
    """Test second order system."""
    def __init__(self):
        System.__init__(self, 2, 1, 0, 'Euler')

    def dyn(self, t, x, u, p=None):
        y1 = x[1]
        y2 = u[0]
        return np.array([y1, y2])

    def jac_dyn(self, t, x, u, p=None):
        y1 = x[1]
        y2 = u[0]
        J = np.zeros((2, 4))
        J[0, 2] = 1.0
        J[1, 3] = 1.0
        return np.array([y1, y2]), coo_matrix(J)


class Pendulum(System):
    """Test pendulum nonlinearity."""
    def __init__(self):
        System.__init__(self, 2, 1, 0, 'Euler')

    def jac_dyn(self, t, x, u, p=None):
        y1 = x[1]
        y2 = u[0] / 5 - 5 * np.sin(x[0])
        y = np.array([y1, y2])
        J = coo_matrix(([1, -5 * np.cos(x[0]), 0.2], ([0, 1, 1], [2, 1, 3])))
        return y, J


class QuadCost(NonLinearPointObj):
    """A quadratic cost on control."""
    def __init__(self):
        NonLinearPointObj.__init__(self, -1, 2, 1, 0, 'user', nG=1)
        self.R = 1.0

    def __callf__(self, x, y):
        u = x[3]
        y[0] = u * self.R * u
        return 0

    def __callg__(self, x, y, G, row, col, rec, needg):
        u = x[3]
        y[0] = u * self.R * u
        if needg:
            G[0] = 2 * self.R * u
            if rec:
                row[0] = 0
                col[0] = 3
        return (0, 0)


def main():
    args = get_onoff_args('grad', 'pen', 'lqr', 'backend ipopt')
    if args.grad:
        gradmode(args)
    if args.pen:
        penMode(args)


def penMode(args):
    """Run pendulum swing-up problem"""
    sys = Pendulum()
    cost = QuadCost()
    N = 40
    t0 = 0.0
    tf = 9.0
    prob = TrajOptProblem(sys, N, t0, tf, gradmode=True)
    prob.xbd = [np.array([-1e20, -1e20]), np.array([1e20, 1e20])]
    prob.ubd = [np.array([-1.0]), np.array([1.0])]
    prob.x0bd = [np.array([0, 0]), np.array([0, 0])]
    prob.xfbd = [np.array([np.pi, 0]), np.array([np.pi, 0])]
    if not args.lqr:
        prob.add_obj(cost, True)  # add a path cost
    else:
        lqr = LqrObj(R=np.ones(1))
        prob.add_lqr_obj(lqr)
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = OptConfig(args.backend)
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print(rst.flag)
    if rst.flag == 1:
        print(rst.sol)
        # parse the solution
        sol = prob.parse_sol(rst.sol.copy())
        show_sol(sol)


def gradmode(args):
    """Solve the simple problem with gradient."""
    sys = OneDcase()
    cost = QuadCost()
    N = 20
    t0 = 0.0
    tf = 2.0
    prob = TrajOptProblem(sys, N, t0, tf, gradmode=True)
    prob.xbd = [np.array([-1e20, -1e20]), np.array([1e20, 1e20])]
    prob.ubd = [np.array([-1e20]), np.array([1e20])]
    prob.x0bd = [np.array([0, 0]), np.array([0, 0])]
    prob.xfbd = [np.array([1, 0]), np.array([1, 0])]
    if not args.lqr:
        prob.add_obj(cost, True)  # add a path cost
    else:
        lqr = LqrObj(R=np.ones(1))
        prob.add_lqr_obj(lqr)
    prob.preProcess()  # construct the problem
    # construct a solver for the problem
    cfg = OptConfig(args.backend)
    slv = OptSolver(prob, cfg)
    rst = slv.solve_rand()
    print(rst.flag, rst.sol)
    if rst.flag == 1:
        # parse the solution
        sol = prob.parse_sol(rst.sol.copy())
        show_sol(sol)


if __name__ == '__main__':
    main()
